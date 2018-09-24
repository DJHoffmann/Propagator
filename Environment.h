/***************************************************************************************************************/
// Header file to contain enviromental definitions/functions for Propagator.cpp 

/***************************************************************************************************************/
// Notes:

/***************************************************************************************************************/
#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include<Constants_Conversions.h>
#include<Input.h>

/***************************************************************************************************************/
// Class to create spatial axes in cylindrically symmetric space:
class spatial_axes {
 private:
  /***************************************************************************************************************/
  // Private Class Function Prototypes:
  void Set_and_Convert();
  void Initialise();
  void Make_Axes();
 public:
  /***************************************************************************************************************/
  // Public Class Members:
  int Nr, Nz;
  int Nr_mpi;
  double dr, dz, z0;
  double *r, *z;
  double *AB; // radial absorbing boundary
  /***************************************************************************************************************/
  // Public Class Function Prototypes:
  spatial_axes();
  spatial_axes(double, double, int, int, int);
  spatial_axes(double, double, double, int, int, int);
  spatial_axes(double, double, double, int, int, double, int);
  spatial_axes(std::ifstream &, int);
  ~spatial_axes();
  void Radial_Absorbing_Boundary(double);
  void Output(double *, int, std::ofstream &);
  void Output(std::ofstream &, std::ofstream &);
};

/***************************************************************************************************************/
// Class Functions:

// Default constructor
spatial_axes::spatial_axes() {
  Nr = 10;
  dr = 0.1;
  Nz = 20;
  dz = 0.1;
  z0 = -0.5*Nz*dz;
  Nr_mpi = Nr;
  Set_and_Convert();
  Initialise();
  Make_Axes();
}

// Constructor for axes with z-axis centred at zero
spatial_axes::spatial_axes(double dr_, double dz_, int Nr_, int Nz_, int mpi_size=1) {
  Nr = Nr_;
  dr = dr_;
  Nz = Nz_;
  dz = dz_;
  z0 = -0.5*Nz*dz;
  Nr_mpi = Nr/mpi_size;
  Set_and_Convert();
  Initialise();
  Make_Axes();
}

// Constructor for axes with z-axis centred at specified z0
spatial_axes::spatial_axes(double z0_, double dr_, double dz_, int Nr_, int Nz_, int mpi_size=1) {
  Nr = Nr_;
  dr = dr_;
  Nz = Nz_;
  dz = dz_;
  z0 = z0_;
  Nr_mpi = Nr/mpi_size;
  Set_and_Convert();
  Initialise();
  Make_Axes();
}

// Constructor for axes with z-axis centred at specified z0 and including radial absorbing boundary
spatial_axes::spatial_axes(double z0_, double dr_, double dz_, int Nr_, int Nz_, double frac=1.0, int mpi_size=1) {
  Nr = Nr_;
  dr = dr_;
  Nz = Nz_;
  dz = dz_;
  z0 = z0_;
  Nr_mpi = Nr/mpi_size;
  Set_and_Convert();
  Initialise();
  Make_Axes();
  Radial_Absorbing_Boundary(frac);
}

// Constructor as above with parameters drawn from specified filestream
spatial_axes::spatial_axes(std::ifstream &INSTREAM, int mpi_size=1) {
  Nr = Read_Variable<int>(INSTREAM, "Nr");
  dr = Read_Variable<double>(INSTREAM, "rMax")/(Nr-1);
  char zRel = Read_Variable<char>(INSTREAM, "zRel");
  double zj = Read_Variable<double>(INSTREAM, "zj");
  Nz = Read_Variable<int>(INSTREAM, "Nz");
  if(zRel=='y') {
    //std::cout << "Adjusting for relative zMin, zMax\n";
    z0 = zj+Read_Variable<double>(INSTREAM, "zMin");
    dz = (zj+Read_Variable<double>(INSTREAM, "zMax")-z0)/(Nz-1); 
  }
  else {
    //std::cout << "Taking absolute zMin, zMax\n";
    z0 = Read_Variable<double>(INSTREAM, "zMin");
    dz = (Read_Variable<double>(INSTREAM, "zMax")-z0)/(Nz-1); 
  }
  Nr_mpi = Nr/mpi_size;
  Set_and_Convert();
  Initialise();
  Make_Axes();
  Radial_Absorbing_Boundary(Read_Variable<double>(INSTREAM, "rAB"));
}

// Destructor
spatial_axes::~spatial_axes() {
  delete[] r;
  delete[] z;
  delete[] AB;
}

// Converts parameters to a.u.
void spatial_axes::Set_and_Convert() {
  dr *= m_to_au;
  dz *= m_to_au;
  z0 *= m_to_au; 
}

// Initialise memory for spatial axes variables
void spatial_axes::Initialise() {
  r = new double[Nr];
  z = new double[Nz];
  AB = new double[Nr];
}

// Make spatial axes
void spatial_axes::Make_Axes() {
  for(int nr=0; nr<Nr; nr++)
    r[nr] = nr*dr;
  for(int nz=0; nz<Nz; nz++)
    z[nz] = z0+nz*dz;
}

// Function to define radially absorbing boundary
void spatial_axes::Radial_Absorbing_Boundary(double frac) {
  double barrier = frac*(Nr-1);
  for(int nr=0; nr<Nr; nr++)
    if(nr<barrier)
      AB[nr] = 1.0;
    else
      AB[nr] = exp(-1.0*pow((nr-barrier)/(Nr-1-barrier), 2));
  //AB[nr] = pow(cos((nr-barrier)/(Nr-1-barrier)*pi/2), 2);
}

// Output specified spatial axis
void spatial_axes::Output(double *Ax, int N, std::ofstream &OUTSTREAM) {
  for(int n=0; n<N; n++)
    OUTSTREAM << Ax[n] << '\n';
}

// Output both spatial axes 
void spatial_axes::Output(std::ofstream &OUTSTREAM_R, std::ofstream &OUTSTREAM_Z) {
  for(int nr=0; nr<Nr; nr++) 
    OUTSTREAM_R << r[nr] << '\n';
  for(int nz=0; nz<Nz; nz++) 
    OUTSTREAM_Z << z[nz] << '\n';
}

/***************************************************************************************************************/
// Inherited spatial class to include also far-field axis information:
class spatial_nearfar_axes : public spatial_axes {
 public:
  /***************************************************************************************************************/
  // Public Class Members:
  double drp;
  double L;
  /***************************************************************************************************************/
  // Public Class Function Prototypes:
  spatial_nearfar_axes();
  spatial_nearfar_axes(std::ifstream &, int);
  void Output(std::ofstream &, std::ofstream &, std::ofstream &);
};

/***************************************************************************************************************/
// Class Functions:

// Default constructor
spatial_nearfar_axes::spatial_nearfar_axes() : spatial_axes() {
  L = 1.0*m_to_au;
  /* drp = std::max(dr,L*tan(20e-3)/Nr); // set for max divergence of 20mrad  */
  drp = std::max(dr,L*tan(30e-3)/Nr); 
}

// Constructor with parameters drawn from specified filestream
spatial_nearfar_axes::spatial_nearfar_axes(std::ifstream &INSTREAM, int mpi_size=1) : spatial_axes(INSTREAM, mpi_size) {
  L = Read_Variable<double>(INSTREAM, "zFar")*m_to_au;
  /* drp = std::max(dr,L*tan(20e-3)/Nr); */
  drp = std::max(dr,L*tan(30e-3)/Nr);
}

// Output all spatial axes 
void spatial_nearfar_axes::Output(std::ofstream &OUTSTREAM_RP, std::ofstream &OUTSTREAM_R, std::ofstream &OUTSTREAM_Z) {
  double scaling = drp/dr;
  for(int nr=0; nr<Nr; nr++) 
    OUTSTREAM_RP << r[nr]*scaling << '\t' << atan2(r[nr]*scaling,L) << '\n';
  spatial_axes::Output(OUTSTREAM_R, OUTSTREAM_Z);
}

/***************************************************************************************************************/
// Class to create the temporal (t & omega) axes:
class temporal_axes {
 private:
  /***************************************************************************************************************/
  // Private Class Function Prototypes:
  void Set_and_Convert();
  void Initialise();
  void Make_Axes();
 public:
  /***************************************************************************************************************/
  // Public Class Members:
  int Nt;
  int Nt_mpi;
  double dt, dw;
  double *t, *w;
  double *AB; // can use for time or frequency
  /***************************************************************************************************************/
  // Public Class Function Prototypes:
  temporal_axes();
  temporal_axes(double, int, int);
  temporal_axes(double, int, double, int);
  temporal_axes(std::ifstream &, int);
  ~temporal_axes();
  void Temporal_Absorbing_Boundary(double);
  void Output(double *, std::ofstream &);
  void Output(std::ofstream &, std::ofstream &);
};

/***************************************************************************************************************/
// Class Functions:

// Default constructor
temporal_axes::temporal_axes() {
  Nt = 10;
  dt = 0.1;
  dw = (2.0*pi)/(Nt*dt);
  Nt_mpi = Nt;
  Set_and_Convert();
  Initialise();
  Make_Axes();
}

// Constructor for temporal axes without specifying absorbing boundary
temporal_axes::temporal_axes(double dt_, int Nt_, int mpi_size=1) {
  Nt = Nt_;
  dt = dt_;
  dw = (2.0*pi)/(Nt*dt);
  Nt_mpi = Nt/mpi_size;
  Set_and_Convert();
  Initialise();
  Make_Axes();
}

// Constructor for temporal axes including absorbing boundary
temporal_axes::temporal_axes(double dt_, int Nt_, double frac=1.0, int mpi_size=1) {
  Nt = Nt_;
  dt = dt_;
  dw = (2.0*pi)/(Nt*dt);
  Nt_mpi = Nt/mpi_size;
  Set_and_Convert();
  Initialise();
  Make_Axes();
  Temporal_Absorbing_Boundary(frac);
}

// Constructor as above with parameters drawn from specified filestream
temporal_axes::temporal_axes(std::ifstream &INSTREAM, int mpi_size=1) {
  Nt = Read_Variable<int>(INSTREAM, "Nt");  // currently must be divisible by mpi_size
  dt = Read_Variable<double>(INSTREAM, "tRange")/(Nt-1); 
  dw = (2.0*pi)/(Nt*dt);
  Nt_mpi = Nt/mpi_size;
  Set_and_Convert();
  Initialise();
  Make_Axes();
  Temporal_Absorbing_Boundary(Read_Variable<double>(INSTREAM, "tAB"));
}

// Destructor
temporal_axes::~temporal_axes() {
  delete[] t;
  delete[] w;
  delete[] AB;
}

// Converts parameters to a.u.
void temporal_axes::Set_and_Convert() {
  dt *= s_to_au;
  dw /= s_to_au;
}

// Initialise memory time and frequency axes variables
void temporal_axes::Initialise() {
  t = new double[Nt];
  w = new double[Nt];
  AB = new double[Nt];
}

// Makes time and frequency axes
void temporal_axes::Make_Axes() {
  for(int nt=0; nt<Nt; nt++)
    t[nt] = (nt-Nt/2)*dt;
  for(int nt=0; nt<Nt; nt++)
    if(nt>=0 && nt<=Nt/2)
      w[nt] = nt*dw;
    else if(nt>Nt/2 && nt<Nt)
      w[nt] = (nt-Nt)*dw;
    else
      w[nt] = 0.0;
}

// Temporal dampening
void temporal_axes::Temporal_Absorbing_Boundary(double frac) {
  double barrier = std::max(0.5,frac)*(Nt-1);
  if(frac>=1.0)
    for(int nt=0; nt<Nt; nt++)
      AB[nt] = 1.0;
  else
    for(int nt=0; nt<Nt; nt++)
      if(nt<=Nt-1-barrier)
	AB[nt] = pow(sin(nt/(Nt-1-barrier)*pi/2.0), 2);
      else if(nt>=barrier)
	AB[nt] = pow(cos((nt-barrier)/(Nt-1-barrier)*pi/2.0), 2);
      else
	AB[nt] = 1.0; 

/*   double barrier = frac*(Nt-1); */
/*   for(int nt=Nt/2; nt<Nt; nt++) */
/*     if(nt<barrier) */
/*       AB[nt] = 1.0; */
/*     else */
/*       AB[nt] = exp(-10.0*pow(3.0*(nt-barrier)/(Nt/2-1-barrier), 2)); */
/*   for(int nt=1; nt<Nt/2; nt++) */
/*     AB[nt] = AB[Nt-nt]; */
/*   AB[0] = AB[1]; */
}

// Output specified temporal axis
void temporal_axes::Output(double *Ax, std::ofstream &OUTSTREAM) {
  for(int nt=0; nt<Nt; nt++)
    OUTSTREAM << Ax[nt] << '\n';
}

// Output both temporal axes 
void temporal_axes::Output(std::ofstream &OUTSTREAM_T, std::ofstream &OUTSTREAM_W) {
  for(int nt=0; nt<Nt; nt++) {
    OUTSTREAM_T << t[nt] << '\n';
    OUTSTREAM_W << w[nt] << '\n';
  }
}

/***************************************************************************************************************/
#endif
