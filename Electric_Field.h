/***************************************************************************************************************/
// Header file to contain electric field class and methods for Propagator.cpp 

/***************************************************************************************************************/
// Notes:
// in Make_spacetime for two-dimensional field the IFT performed for both components (new should equal old for x-component)

/***************************************************************************************************************/
#ifndef ELECTRIC_FIELD_H
#define ELECTRIC_FIELD_H

#include<cmath>
#include<fftw3.h>
#ifdef useMPI
#ifdef oldMPI
#undef SEEK_SET  // undef needed due to a bug with MPI-2, bug fixed in openmpi >= 1.3.2
#undef SEEK_END
#undef SEEK_CUR
#endif
#include<mpi.h>
#endif
#include<Constants_Conversions.h>
#include<Misc_Functions.h>
#include<Environment.h>
#include<Laser_Parameters.h>
#include<Gas_Jet.h>

/***************************************************************************************************************/
// Governing electric_field abstract class to allow for selection of the following one-/two-dimensional field classes at runtime
class electric_field {
  protected:
  /***************************************************************************************************************/
  // Protected Class Members:
  temporal_axes *pT;
  spatial_axes *pX;
  //laser_parameters **pLP;
  bool spatial, FFTW; 
  bool mpi; // whether we need E_mpi (will often import directly as mpi_active from main.cpp)
  dcmplx *Er; // not currently needed outside of class
  fftw_complex *fftwArrayT;
  fftw_plan time2freq, freq2time;
 public:
  /***************************************************************************************************************/
  // Public Class Members:
  int Dim;
  laser_parameters **pLP;
  int colours;
  double maxEt;
  dcmplx *Et, *Ew, *Ert, *Erw, *Et2, *Ew2, *Ert2, *Erw2;
  dcmplx *E_mpi, *E_mpi2; // use for segments of both Erw,Erw2 and Ertp,Ertp2
  /***************************************************************************************************************/
  // Public Virtual Class Functions:
  virtual double Pulse_Energy(char) = 0; 
  virtual void Transform_Field_FT() = 0;
  virtual void Transform_Field_FT(fftw_complex *, fftw_plan) = 0;
  virtual void Transform_Field_IFT() = 0;
  virtual void Transform_Field_IFT(fftw_complex *, fftw_plan) = 0;
  virtual void Update_Et_at_RZ(int) = 0;
  virtual void Update_Ew_at_RZ(int) = 0;
  virtual void Transpose_Overwrite(char) = 0;
  virtual void Combine_Fields(char) = 0;
  virtual void Check_TemporalAxes_Pointer(temporal_axes &) = 0;
  virtual void Check_SpatialAxes_Pointer(spatial_axes &) = 0;
  virtual void Check_LaserParameters_Pointer(laser_parameters **) = 0;
  virtual void Analyse_Field(laser_parameters **&, bool) = 0;
  virtual void Analyse_Field(laser_parameters **&, int, bool, bool) = 0;
  virtual void Analyse_Field(laser_parameters **&, dcmplx *, bool) = 0;
  virtual void Analyse_Field(laser_parameters **&, dcmplx *, dcmplx *, bool) = 0;
  virtual void Output_Temporal(dcmplx *, std::ofstream &, char='r') = 0;
  virtual void Output_Radial(dcmplx *, int, char, std::ofstream &, char='r') = 0;
  virtual void Output_RadialTemporal(dcmplx *, char, std::ofstream &, char='r', bool=false, int=1, int=1) = 0;
};

/***************************************************************************************************************/
#include<Diffraction.h>

/***************************************************************************************************************/
// electric_field_1d class for electric fields linearly polarised along one axis:
class electric_field_1d : public electric_field {
 protected:
  /***************************************************************************************************************/
  // Protected Class Function Prototypes:
  void Initialise_Time();
  void Make_Time(fftw_complex *, fftw_plan);
  void Temporal_Field(dcmplx *, const int, char, bool, bool);
  void Temporal_Field(bool, bool, bool);
  void Calculate_MaxValue();
  void Reset_Field();
  void Initialise_SpaceTime();
  void Make_SpaceTime(fftw_complex *, fftw_plan, fftw_plan);
  void Make_SpaceTime(fftw_complex *, fftw_plan, fftw_plan, double);
  void Gaussian_Beam(double, double, int, bool);
  void Vector_Potential(dcmplx *, dcmplx *);
  void Analyse_SingleColour_Field(laser_parameters **&, dcmplx *, bool);
  void Analyse_MultiColour_Field(laser_parameters **&, dcmplx *, bool);
  void Destroy_Time();
  void Destroy_SpaceTime();
 public:
  /***************************************************************************************************************/
  // Public Class Function Prototypes:
  electric_field_1d(temporal_axes &);
  electric_field_1d(temporal_axes &, laser_parameters **, int);
  electric_field_1d(temporal_axes &, laser_parameters **, fftw_complex *, fftw_plan, int);
  electric_field_1d(temporal_axes &, spatial_axes &, bool);
  electric_field_1d(temporal_axes &, spatial_axes &, laser_parameters **, int, bool);
  electric_field_1d(temporal_axes &, spatial_axes &, laser_parameters **, double, int, bool);
  electric_field_1d(temporal_axes &, spatial_axes &, laser_parameters **, fftw_complex *, fftw_plan, fftw_plan, int, bool);
  ~electric_field_1d();
  double Pulse_Energy(char); 
  void Transform_Field_FT();
  void Transform_Field_FT(fftw_complex *, fftw_plan);
  void Transform_Field_IFT();
  void Transform_Field_IFT(fftw_complex *, fftw_plan);
  void Update_Et_at_RZ(int);
  void Update_Ew_at_RZ(int);
  void Transpose_Overwrite(char);
  void Combine_Fields(char);
  void Check_TemporalAxes_Pointer(temporal_axes &);
  void Check_SpatialAxes_Pointer(spatial_axes &);
  void Check_LaserParameters_Pointer(laser_parameters **);
  void Analyse_Field(laser_parameters **&, bool);
  void Analyse_Field(laser_parameters **&, int, bool, bool);
  void Analyse_Field(laser_parameters **&, dcmplx *, bool);
  void Analyse_Field(laser_parameters **&, dcmplx *, dcmplx *, bool);
  void Output_Temporal(dcmplx *, std::ofstream &, char);
  void Output_Radial(dcmplx *, int, char, std::ofstream &, char);
  void Output_RadialTemporal(dcmplx *, char, std::ofstream &, char, bool, int, int);
};

/***************************************************************************************************************/
// Class Functions:

// Constructor for to-be-calculated temporal fields Et, Ew
electric_field_1d::electric_field_1d(temporal_axes &T_) {
  Dim = 1;
  pT = &T_;
  pX = NULL;
  pLP = NULL;
  colours = 0;
  spatial = false;
  FFTW = true;
  mpi = false;
  Initialise_Time();
}

// Constructor for temporal fields Et, Ew calculated at creation with parameters LP
electric_field_1d::electric_field_1d(temporal_axes &T_, laser_parameters **pLP_, int colours_=1) {
  Dim = 1;
  pT = &T_;
  pX = NULL;
  pLP = pLP_;
  colours = colours_;
  spatial = false;
  FFTW = true;
  mpi = false;
  Initialise_Time();
  Make_Time(fftwArrayT, time2freq);
}

// Constructor as previous but takes in FFTW arguments
electric_field_1d::electric_field_1d(temporal_axes &T_, laser_parameters **pLP_, fftw_complex *fftwArrayT_, fftw_plan time2freq_, int colours_=1) {
  Dim = 1;
  pT = &T_;
  pX = NULL;
  pLP = pLP_;
  colours = colours_;
  spatial = false;
  FFTW = false;
  mpi = false;
  Initialise_Time();
  Make_Time(fftwArrayT_, time2freq_);
}

// Constructor for to-be-calculated temporal fields Et, Ew; spatial field Er; temporo-spatial fields Ert, Erw
electric_field_1d::electric_field_1d(temporal_axes &T_, spatial_axes &X_, bool mpi_=false) {
  Dim = 1;
  pT = &T_;
  pX = &X_;
  pLP = NULL;
  spatial = true;
  FFTW = true;
  mpi = mpi_;
  Initialise_SpaceTime();
}

// Constructor for temporal fields Et, Ew; spatial field Er; temporo-spatial fields Ert, Erw defined via parameters L
electric_field_1d::electric_field_1d(temporal_axes &T_, spatial_axes &X_, laser_parameters **pLP_, int colours_=1, bool mpi_=false) {
  Dim = 1;
  pT = &T_;
  pX = &X_;
  pLP = pLP_;
  colours = colours_;
  spatial = true;
  FFTW = true;
  mpi = mpi_;
  Initialise_SpaceTime();
  Make_SpaceTime(fftwArrayT, time2freq, freq2time);
}

// As above but for laser parameters defined at a given z=zP (ie. different from the focus, usually centre of gas jet)
electric_field_1d::electric_field_1d(temporal_axes &T_, spatial_axes &X_, laser_parameters **pLP_, double zP_, int colours_=1, bool mpi_=false) {
  Dim = 1;
  pT = &T_;
  pX = &X_;
  pLP = pLP_;
  colours = colours_;
  spatial = true;
  FFTW = true;
  mpi = mpi_;
  Initialise_SpaceTime();
  Make_SpaceTime(fftwArrayT, time2freq, freq2time, zP_);
}

// Constructor as previous but takes in FFTW arguments
electric_field_1d::electric_field_1d(temporal_axes &T_, spatial_axes &X_, laser_parameters **pLP_, fftw_complex *fftwArrayT_, fftw_plan time2freq_, fftw_plan freq2time_, int colours_=1, bool mpi_=false) {
  Dim = 1;
  pT = &T_;
  pX = &X_;
  pLP = pLP_;
  colours = colours_;
  spatial = true;
  FFTW = false;
  mpi = mpi_;
  Initialise_SpaceTime();
  Make_SpaceTime(fftwArrayT_, time2freq_, freq2time_);
}

// Destructor
electric_field_1d::~electric_field_1d() {
  if(spatial==false)  Destroy_Time();
  else  Destroy_SpaceTime();
}

// Initialisation of variable memory for temporal-only
void electric_field_1d::Initialise_Time() {
  Et = new dcmplx[pT->Nt];
  Ew = new dcmplx[pT->Nt];
  if(FFTW) {
#ifdef useMPI 
    int rank = 0, size = 1;  // force sequential for thread-safety 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);  MPI_Comm_size(MPI_COMM_WORLD,&size);
    for(int n=0; n<size; n++) {
      if(rank==n)  {
	fftwArrayT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*pT->Nt);
	time2freq = fftw_plan_dft_1d(pT->Nt, fftwArrayT, fftwArrayT, FFTW_FORWARD, FFTW_MEASURE);
	freq2time = fftw_plan_dft_1d(pT->Nt, fftwArrayT, fftwArrayT, FFTW_BACKWARD, FFTW_MEASURE);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
#else
    fftwArrayT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*pT->Nt);
    time2freq = fftw_plan_dft_1d(pT->Nt, fftwArrayT, fftwArrayT, FFTW_FORWARD, FFTW_MEASURE);
    freq2time = fftw_plan_dft_1d(pT->Nt, fftwArrayT, fftwArrayT, FFTW_BACKWARD, FFTW_MEASURE);
#endif
  }
  else {
    fftwArrayT = NULL;
    time2freq = NULL;  freq2time = NULL;
  }
  Et2 = Ew2 = NULL;
}

// Calculation of temporal fields
void electric_field_1d::Make_Time(fftw_complex *fftwArrayT_, fftw_plan time2freq_) {
  Temporal_Field(true, false, true);
  if(norm(Et[0])>1.0e-4)  std::cout << "WARNING! |Et[0]|^2>1e-4, killing tails through absorbing boundary\n";
  for(int nt=0; nt<pT->Nt; nt++) {
    Et[nt] *= pT->AB[nt]; // apply temporal absorbing boundary
    fftwArrayT_[nt][0] = real(Et[nt]);
    fftwArrayT_[nt][1] = imag(Et[nt]);
    //fftwArrayT_[nt][1] = 0.0;
  }
  fftw_execute(time2freq_);
  double prefac = pT->dt/sqrt(2.0*pi);
  for(int nt=0; nt<pT->Nt; nt++)
    Ew[nt] = prefac * (fftwArrayT_[nt][0]+i*fftwArrayT_[nt][1]);
  Calculate_MaxValue();
}

// Creates single colour temporal field using specified envelope type and method (direct or via vector potential)
void electric_field_1d::Temporal_Field(dcmplx *E, const int nc, char PolSwitch='n', bool Adefn=false, bool realise=false) {
  double env_range;   // Pol_Switch = 'x/c' to project onto x-axis (using phi), 'y/s' for y-axis
  double pol_factor = (PolSwitch=='x'||PolSwitch=='c') ? cos(pLP[nc]->phi) : ((PolSwitch=='y'||PolSwitch=='s') ? sin(pLP[nc]->phi) : 1.0);
  switch(pLP[nc]->envelope) {
    case 'g':
      if(Adefn==true)
        for(int nt=0; nt<pT->Nt; nt++)
	  E[nt] += pLP[nc]->E0*(1.0-2.0*i*(pT->t[nt]-pLP[nc]->t0)*pLP[nc]->B/pLP[nc]->w0)*exp(-pow(pT->t[nt]-pLP[nc]->t0, 2)*pLP[nc]->B)*exp(i*(pLP[nc]->w0*(pT->t[nt]-pLP[nc]->t0)+pLP[nc]->cep)) * pol_factor;
      else
        for(int nt=0; nt<pT->Nt; nt++)
 	  E[nt] += pLP[nc]->E0 * exp(-pow(pT->t[nt]-pLP[nc]->t0, 2)*pLP[nc]->B) * exp(i*(pLP[nc]->w0*(pT->t[nt]-pLP[nc]->t0)+pLP[nc]->cep)) * pol_factor;
      break;
    case 'c':
      env_range = 0.5*pi/pLP[nc]->B;
      if(Adefn==true)
        for(int nt=0; nt<pT->Nt; nt++)
          if((pT->t[nt]-pLP[nc]->t0)>=-env_range && (pT->t[nt]-pLP[nc]->t0)<=env_range) 
	    E[nt] += pLP[nc]->E0/pLP[nc]->w0 * cos((pT->t[nt]-pLP[nc]->t0)*pLP[nc]->B) * (pLP[nc]->w0*cos((pT->t[nt]-pLP[nc]->t0)*pLP[nc]->B)+i*2.0*pLP[nc]->B*sin((pT->t[nt]-pLP[nc]->t0)*pLP[nc]->B)) * exp(i*(pLP[nc]->w0*(pT->t[nt]-pLP[nc]->t0)+pLP[nc]->cep)) * pol_factor;
          else
            E[nt] = 0.0; 
      else
        for(int nt=0; nt<pT->Nt; nt++)
          if((pT->t[nt]-pLP[nc]->t0)>=-env_range && (pT->t[nt]-pLP[nc]->t0)<=env_range) 
	    E[nt] += pLP[nc]->E0 * pow(cos((pT->t[nt]-pLP[nc]->t0)*pLP[nc]->B), 2) * exp(i*(pLP[nc]->w0*(pT->t[nt]-pLP[nc]->t0)+pLP[nc]->cep)) * pol_factor;
          else
            E[nt] = 0.0; 
      break;
    case 's':
      if(Adefn==true)
        for(int nt=0; nt<pT->Nt; nt++)
	  E[nt] += pLP[nc]->E0/pLP[nc]->w0/cosh((pT->t[nt]-pLP[nc]->t0)*pLP[nc]->B) * (pLP[nc]->w0+i*pLP[nc]->B*tanh((pT->t[nt]-pLP[nc]->t0)*pLP[nc]->B)) * exp(i*(pLP[nc]->w0*(pT->t[nt]-pLP[nc]->t0)+pLP[nc]->cep)) * pol_factor;
      else
        for(int nt=0; nt<pT->Nt; nt++)
	  E[nt] += pLP[nc]->E0 / cosh((pT->t[nt]-pLP[nc]->t0)*pLP[nc]->B) * exp(i*(pLP[nc]->w0*(pT->t[nt]-pLP[nc]->t0)+pLP[nc]->cep)) * pol_factor;
      break;
    default:
      std::cerr << "Incorrect envelope type for colour " << nc+1 << '\n';
#ifdef useMPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(1);
  }
  if(realise==true)  
    for(int nt=0; nt<pT->Nt; nt++)  E[nt] = 0.5*(E[nt]+conj(E[nt]));
}

// Function to create multicolour temporal field
void electric_field_1d::Temporal_Field(bool Adefn=false, bool realise=false, bool Ereset=false) {
  if(Ereset==true)  Reset_Field();
  for(int nc=0;nc<colours; nc++)
    Temporal_Field(Et, nc, 'x', true, false);
}

// Calculates temporal field maximum
void electric_field_1d::Calculate_MaxValue() {
  maxEt = 0;
  for(int nt=0; nt<pT->Nt; nt++)
    if(abs(Et[nt])>maxEt)
      maxEt = abs(Et[nt]);
}

// Function to reset Et to zero
void electric_field_1d::Reset_Field() {
  for(int nt=0; nt<pT->Nt; nt++) {
    Et[nt] = 0.0;
    Ew[nt] = 0.0;
  }
}

// Initialisation of variable memory for temporal-spatial (which ones depends on the class member "spatial")
void electric_field_1d::Initialise_SpaceTime() {
  Initialise_Time();
  Er = new dcmplx[pX->Nr];
  Erw = new dcmplx[pX->Nr*pT->Nt];
  Ert = new dcmplx[pX->Nr*pT->Nt];
  E_mpi = (mpi==true) ? new dcmplx[pX->Nr_mpi*pT->Nt] : NULL; // could also be pX->Nr*pT->Nt_mpi
}

// Calculation of temporal-spatial fields (which ones depends on the class member "spatial")
void electric_field_1d::Make_SpaceTime(fftw_complex *fftwArrayT_, fftw_plan time2freq_, fftw_plan freq2time_) {
  double k;
  double pre_factor = pT->dt/sqrt(2.0*pi); 
  for(int nt=0; nt<pT->Nt; nt++)
    for(int nr=0; nr<pX->Nr; nr++) 
      Erw[Index2V(nt,nr,pX->Nr)] = 0.0;
  for(int nc=0;nc<colours; nc++) {
    Reset_Field();
    Temporal_Field(Et, nc, 'x', true, false);
    for(int nt=0; nt<pT->Nt; nt++) {
      fftwArrayT_[nt][0] = real(Et[nt]);
      fftwArrayT_[nt][1] = imag(Et[nt]);
      //fftwArrayT_[nt][1] = 0.0;
    }
    fftw_execute(time2freq_);
    for(int nt=0; nt<pT->Nt; nt++) {
      Ew[nt] = pre_factor * (fftwArrayT_[nt][0]+i*fftwArrayT_[nt][1]);
      k = pT->w[nt]/c_au; 
      Gaussian_Beam(k, pX->z0, nc, true);
      for(int nr=0; nr<pX->Nr; nr++)
	Erw[Index2V(nt,nr,pX->Nr)] += Er[nr]*Ew[nt];
    }
  }
  Transform_Field_IFT(fftwArrayT_, freq2time_);
  Make_Time(fftwArrayT_, time2freq_);  // calculates combined multicolour Et, Ew independent of space (on-axis temporal field)
  Ert2 = Erw2 = E_mpi2 = NULL;
}

// As above but for laser parameters defined at a given z=zP (ie. different from the focus, usually centre of gas jet)
void electric_field_1d::Make_SpaceTime(fftw_complex *fftwArrayT_, fftw_plan time2freq_, fftw_plan freq2time_, double zP) {
  bool diffract = true;
  if(zP<pX->z0 || zP>pX->z0+(pX->Nz-1)*pX->dz) {  // check validity of zP 
    std::cerr << "Error in choice of z position where laser parameters are defined, resetting to laser focus\n";
    zP = pX->z0;
    diffract = false;
  }
  double k;
  double pre_factor = pT->dt/sqrt(2.0*pi); 
  for(int nt=0; nt<pT->Nt; nt++)
    for(int nr=0; nr<pX->Nr; nr++) 
      Erw[Index2V(nt,nr,pX->Nr)] = 0.0;
  for(int nc=0;nc<colours; nc++) {
    Reset_Field();
    Temporal_Field(Et, nc, 'x', true, false);
    for(int nt=0; nt<pT->Nt; nt++) {
      fftwArrayT_[nt][0] = real(Et[nt]);
      fftwArrayT_[nt][1] = imag(Et[nt]);
      //fftwArrayT_[nt][1] = 0.0;
    }
    fftw_execute(time2freq_);
    for(int nt=0; nt<pT->Nt; nt++) {
      Ew[nt] = pre_factor * (fftwArrayT_[nt][0]+i*fftwArrayT_[nt][1]);
      k = pT->w[nt]/c_au; 
      Gaussian_Beam(k, zP, nc, !diffract);
      for(int nr=0; nr<pX->Nr; nr++)
	Erw[Index2V(nt,nr,pX->Nr)] += Er[nr]*Ew[nt];
    }
  }
  if(diffract) {
    int mpi_rank = 0; 
    dcmplx *A = new dcmplx[pX->Nr], *B = new dcmplx[pX->Nr], *C = new dcmplx[pX->Nr], *D = new dcmplx[pX->Nr];
    // Diffract laser field backwards to start of propagation axis (frequency domain):
    for(int nz=(int)((zP-pX->z0)/pX->dz); nz>=0; nz--) {
#ifdef useMPI
      MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
      MPI_Scatter(this->Erw, pX->Nr*pT->Nt_mpi, MPI_DOUBLE_COMPLEX, this->E_mpi, pX->Nr*pT->Nt_mpi, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
      //Diffract_Field(*this, A, B, C, D, *pX, *pT, -pX->dz, mpi_rank);
      Diffract_Field(this->E_mpi, A, B, C, D, *pX, *pT, -pX->dz, mpi_rank);
      MPI_Allgather(this->E_mpi, pX->Nr*pT->Nt_mpi, MPI_DOUBLE_COMPLEX, this->Erw, pX->Nr*pT->Nt_mpi, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
#else
      //Diffract_Field(*this, A, B, C, D, *pX, *pT, -pX->dz, mpi_rank);
      Diffract_Field(this->Erw, A, B, C, D, *pX, *pT, -pX->dz, mpi_rank);
#endif
    }
    delete[] A, B, C, D;
  }
  Transform_Field_IFT(fftwArrayT_, freq2time_);
  Make_Time(fftwArrayT_, time2freq_);  // calculates combined multicolour Et, Ew independent of space (on-axis temporal field)
  Ert2 = Erw2 = E_mpi2 = NULL;
}

// Gaussian beam function for spatial field
void electric_field_1d::Gaussian_Beam(double k, double z0, int nc, bool wFac=true) {
  double r2;
  double zr = pow(pLP[nc]->wz0,2)*fabs(k)/2.0;
  double z = z0-pLP[nc]->zf;
  double wz = (z==0) ? pLP[nc]->wz0 : pLP[nc]->wz0*sqrt(1.0+pow(z/zr,2));
  double oRz = (z==0) ? 0.0 : z/(z*z+zr*zr); // something very strange here when z = 0
  dcmplx prefac = wFac ? pLP[nc]->wz0/wz*exp(i*atan2(z,zr)) : 1.0;  // ignore wz0/wz & Gouy phase when parameters not at focus
  for(int nr=0; nr<pX->Nr; nr++) {
    r2 = pow(pX->r[nr],2);
    Er[nr] = prefac * exp(-r2/(wz*wz)) * exp(-i*0.5*k*r2*oRz);
  }
}

// Function to output the (temporo-spatial) pulse energy in Joules (removes intensity ambiguity in a.u.)
double electric_field_1d::Pulse_Energy(char method='n') {
  double Energy = 0.0;
  switch(method) { 
  case 'n':  // numerical integration
    for(int nt=0; nt<pT->Nt; nt++) {
      Energy += 0.5*(pow(real(Ert[Index2V(0,nt,pT->Nt)]),2)+pow(real(Ert[Index2V(1,nt,pT->Nt)]),2))*pX->dr/8.0;
      for(int nr=1; nr<pX->Nr; nr++)
	Energy += pow(real(Ert[Index2V(nr,nt,pT->Nt)]),2)*pX->r[nr];
    }
    Energy *= 2.0*pi*(au_to_Wcm2*1.0e4)*au_to_m*(pX->dr*au_to_m)*(pT->dt*au_to_s) * 2.0; // last 2.0 to convert to cycle averaged power (v accurate to even 0.5 cycle FWHMI)
    break;
  case 'e':  // experimental estimate 
    double factor;
    for(int nc=0; nc<colours; nc++) {
      if(pLP[nc]->envelope=='s') // sech
	factor = 0.88;
      else
	factor = 0.94;  // need to do for cos2 too
      Energy += (pLP[nc]->I0*au_to_Wcm2*1.0e4) * 0.5*pi*(pLP[nc]->wz0*au_to_m)*(pLP[nc]->wz0*au_to_m) * (pLP[nc]->fwhm*au_to_s)/factor;
    }
    break;
  default:
    std::cerr << "Error in method for calculation of pulse energy\n";
  }
  return Energy;
}

//  FFT(W): Ert -> Erw
void electric_field_1d::Transform_Field_FT() {
  double prefac = pT->dt/sqrt(2.0*pi);
  for(int nr=0; nr<pX->Nr; nr++) {
    for(int nt=0; nt<pT->Nt; nt++) {
      fftwArrayT[nt][0] = real(Ert[Index2V(nr,nt,pT->Nt)]);
      fftwArrayT[nt][1] = imag(Ert[Index2V(nr,nt,pT->Nt)]);
      //fftwArrayT[nt][1] = 0.0; // imaginary component of Ert should be zero for Et=real(...) definition
    }
    fftw_execute(time2freq);
    for(int nt=0; nt<pT->Nt; nt++)
      Erw[Index2V(nt,nr,pX->Nr)] = prefac * (fftwArrayT[nt][0]+i*fftwArrayT[nt][1]);
  }
}

//  As above taking fftw array and plan as arguments
void electric_field_1d::Transform_Field_FT(fftw_complex *fftwArrayT_, fftw_plan time2freq_) {
  double prefac = pT->dt/sqrt(2.0*pi);
  for(int nr=0; nr<pX->Nr; nr++) {
    for(int nt=0; nt<pT->Nt; nt++) {
      fftwArrayT_[nt][0] = real(Ert[Index2V(nr,nt,pT->Nt)]);
      fftwArrayT_[nt][1] = imag(Ert[Index2V(nr,nt,pT->Nt)]);
      //fftwArrayT_[nt][1] = 0.0; // imaginary component of Ert should be zero for Et=real(...) definition
    }
    fftw_execute(time2freq_);
    for(int nt=0; nt<pT->Nt; nt++)
      Erw[Index2V(nt,nr,pX->Nr)] = prefac * (fftwArrayT_[nt][0]+i*fftwArrayT_[nt][1]);
  }
}

// Inverse FFT(W): Erw -> Ert
void electric_field_1d::Transform_Field_IFT() {
  double prefac = sqrt(2.0*pi)/(pT->Nt*pT->dt);
  for(int nr=0; nr<pX->Nr; nr++) {
    for(int nt=0; nt<pT->Nt; nt++) {
      fftwArrayT[nt][0] = real(Erw[Index2V(nt,nr,pX->Nr)]);
      fftwArrayT[nt][1] = imag(Erw[Index2V(nt,nr,pX->Nr)]);
    }
    fftw_execute(freq2time);
    for(int nt=0; nt<pT->Nt; nt++)
      Ert[Index2V(nr,nt,pT->Nt)] = prefac * (fftwArrayT[nt][0]+i*fftwArrayT[nt][1]);
    //Ert[Index2V(nr,nt,pT->Nt)] = prefac * fftwArrayT[nt][0]; // imaginary component of Ert should be zero for Et=real(...) definition
  }
}

//  As above taking fftw array and plan as arguments
void electric_field_1d::Transform_Field_IFT(fftw_complex *fftwArrayT_, fftw_plan freq2time_) {
  double prefac = sqrt(2.0*pi)/(pT->Nt*pT->dt);
  for(int nr=0; nr<pX->Nr; nr++) {
    for(int nt=0; nt<pT->Nt; nt++) {
      fftwArrayT_[nt][0] = real(Erw[Index2V(nt,nr,pX->Nr)]);
      fftwArrayT_[nt][1] = imag(Erw[Index2V(nt,nr,pX->Nr)]);
    }
    fftw_execute(freq2time_);
    for(int nt=0; nt<pT->Nt; nt++)
      Ert[Index2V(nr,nt,pT->Nt)] = prefac * (fftwArrayT_[nt][0]+i*fftwArrayT_[nt][1]);
    //Ert[Index2V(nr,nt,pT->Nt)] = prefac * fftwArrayT_[nt][0]; // imaginary component of Ert should be zero for Et=real(...) definition
  }
}

// Updates to Et spatial field Ert at current z and specified r
void electric_field_1d::Update_Et_at_RZ(int nr=0) {
  for(int nt=0; nt<pT->Nt; nt++)
    Et[nt] = Ert[Index2V(nr,nt,pT->Nt)];
}

// Updates to Ew spatial field Erw at current z and specified r
void electric_field_1d::Update_Ew_at_RZ(int nr) {
  for(int nt=0; nt<pT->Nt; nt++)
    Ew[nt] = Erw[Index2V(nr,nt,pT->Nt)];
}

// Function to copy-and-transpose one temporo-spatial field onto the other (use for harmonic field where Et not needed during propagation)
void electric_field_1d::Transpose_Overwrite(char direction) {
  if(spatial==true && direction=='f')
    for(int nr=0; nr<pX->Nr; nr++)
      for(int nt=0; nt<pT->Nt; nt++)
	Erw[Index2V(nt,nr,pX->Nr)] = Ert[Index2V(nr,nt,pT->Nt)];
  else if(spatial==true && direction=='b')
    for(int nt=0; nt<pT->Nt; nt++)
      for(int nr=0; nr<pX->Nr; nr++)
	Ert[Index2V(nr,nt,pT->Nt)] = Erw[Index2V(nt,nr,pX->Nr)];
}

// Function to add one temporo-spatial field onto the other (specified by destination) - used in harmonic propagation
void electric_field_1d::Combine_Fields(char destination) {
  if(spatial==true && destination=='t')
    for(int nr=0; nr<pX->Nr; nr++)
      for(int nt=0; nt<pT->Nt; nt++)
	Ert[Index2V(nr,nt,pT->Nt)] += Erw[Index2V(nt,nr,pX->Nr)];
  else if(spatial==true && destination=='w')
    for(int nt=0; nt<pT->Nt; nt++)
      for(int nr=0; nr<pX->Nr; nr++)
	Erw[Index2V(nt,nr,pX->Nr)] += Ert[Index2V(nr,nt,pT->Nt)];
}

// Function to check temporal pointer
void electric_field_1d::Check_TemporalAxes_Pointer(temporal_axes &T_) {
  if(pT==&T_)
    std::cout << "Matching temporal pointers\n";
  else 
    std::cout << "Temporal pointers do not match!\n";
}

// Function to check spatial pointer
void electric_field_1d::Check_SpatialAxes_Pointer(spatial_axes &X_) {
  if(pX==&X_)
    std::cout << "Matching spatial pointers\n";
  else 
    std::cout << "Spatial pointers do not match!\n";
}

// Function to check laser parameters pointer
void electric_field_1d::Check_LaserParameters_Pointer(laser_parameters **pLP_) {
  if(pLP==pLP_)
    std::cout << "Matching laser parameters pointers\n";
  else 
    std::cout << "Laser parameter pointers do not match!\n";
}

// Function to calculate vector potential from temporal electric field
void electric_field_1d::Vector_Potential(dcmplx *At, dcmplx *Efield) {
  At[0] = 0.0;
  At[1] = -pT->dt*0.5*(Efield[0]+Efield[1]);
  for(int nt=2; nt<pT->Nt; nt++)  At[nt] = At[nt-2]-pT->dt*(Efield[nt-2]+4.0*Efield[nt-1]+Efield[nt])/3.0;
}

// Governing interface to analyse electric field and write parameters found to specified laser_parameters object for temporal-only fields
void electric_field_1d::Analyse_Field(laser_parameters **&pLP_, bool findAll) {
  if(colours==1)  Analyse_SingleColour_Field(pLP_, Et, findAll);
  else if(colours>1)  Analyse_MultiColour_Field(pLP_, Et, findAll);
}

// As above for temporo-spatial fields
void electric_field_1d::Analyse_Field(laser_parameters **&pLP_, int nr, bool findAll, bool mpi) {
  if(colours==1) 
    (mpi==true) ? Analyse_SingleColour_Field(pLP_, &E_mpi[Index2V(nr,0,pT->Nt)], findAll) : Analyse_SingleColour_Field(pLP_, &Ert[Index2V(nr,0,pT->Nt)], findAll);
  else if(colours>1) 
    (mpi==true) ? Analyse_MultiColour_Field(pLP_, &E_mpi[Index2V(nr,0,pT->Nt)], findAll) : Analyse_MultiColour_Field(pLP_, &Ert[Index2V(nr,0,pT->Nt)], findAll);
}

// As above but taking fields in explicitly
void electric_field_1d::Analyse_Field(laser_parameters **&pLP_, dcmplx *Efield, bool findAll) {
  if(colours==1)  Analyse_SingleColour_Field(pLP_, Efield, findAll);
  else if(colours>1)  Analyse_MultiColour_Field(pLP_, Efield, findAll);
}

// Exists purely for compatibility with electric_field_2d
void electric_field_1d::Analyse_Field(laser_parameters **&pLP_, dcmplx *Efield1, dcmplx *Efield2, bool findAll) {
  if(colours==1)  Analyse_SingleColour_Field(pLP_, Efield1, findAll);
  else if(colours>1)  Analyse_MultiColour_Field(pLP_, Efield1, findAll);
}

// Function to analyse single-colour electric field and write parameters found to specified laser_parameters object
void electric_field_1d::Analyse_SingleColour_Field(laser_parameters **&pLP_, dcmplx *Efield, bool findAll) {
  double I0, lambda, fwhm, t0, cep;
  char envelope = 'c';  // fix this as cos^2 envelope for now (for QO)
  int nt_max;
  double Enorm_max = 0.0;
  for(int nt=0; nt<pT->Nt; nt++)
  //for(int nt=pT->Nt/8; nt<7*pT->Nt/8; nt++)  // assume within central 3/4 of the range
    if(norm(Efield[nt])>Enorm_max)  {
      nt_max = nt;
      Enorm_max = norm(Efield[nt]);
    }
  double Enorm_maxM1 = norm(Efield[nt_max-1]);
  double Enorm_maxP1 = norm(Efield[nt_max+1]);
  I0 = Enorm_max - 0.125*pow(Enorm_maxP1-Enorm_maxM1,2)/(Enorm_maxP1-2.0*Enorm_max+Enorm_maxM1); 
  double nt_shift = -0.5*(Enorm_maxP1-Enorm_maxM1)/(Enorm_maxP1-2.0*Enorm_max+Enorm_maxM1); 
  t0 = pT->t[nt_max]+nt_shift*pT->dt;
  double arg_max = atan2(imag(Efield[nt_max]),real(Efield[nt_max]));
  double arg_near = (nt_shift>0.0) ? atan2(imag(Efield[nt_max+1]),real(Efield[nt_max+1])) : atan2(imag(Efield[nt_max-1]),real(Efield[nt_max-1]));
  if(arg_near-arg_max>pi)  arg_near-=2.0*pi;
  else if(arg_near-arg_max<-pi)  arg_near+=2.0*pi;
  cep = (1.0-fabs(nt_shift))*arg_max+fabs(nt_shift)*arg_near;
  if(cep>pi)  cep-=2.0*pi;
  else if(cep<-pi)  cep+=2.0*pi;
  if(findAll==true) {
    dcmplx *At = new dcmplx[pT->Nt];
    Vector_Potential(At, Efield);
    double t_minM=0.0, t_minP=0.0;
    for(int nt=nt_max+2; nt<pT->Nt; nt++)  
      if(real(At[nt-1]*exp(-i*cep))<0.0&&real(At[nt]*exp(-i*cep))>0.0 || real(At[nt-1]*exp(-i*cep))>0.0&&real(At[nt]*exp(-i*cep))<0.0) {
	t_minP = pT->t[nt-1]+pT->dt*fabs(real(At[nt-1]*exp(-i*cep))/real((At[nt]-At[nt-1])*exp(-i*cep)));  
	break;
      }
    lambda = 2.0*(t_minP-t0)*c_au;
    double halfAmodsq = 0.5 * (norm(At[nt_max])-0.125*pow(norm(At[nt_max+1])-norm(At[nt_max-1]),2)/(norm(At[nt_max+1])-2.0*norm(At[nt_max])+norm(At[nt_max-1]))); 
    for(int nt=nt_max+2; nt<pT->Nt; nt++)
      if(norm(At[nt-1])>halfAmodsq && norm(At[nt])<halfAmodsq) {
	fwhm = 2.0*(pT->t[nt-1]+pT->dt*fabs((norm(At[nt-1])-halfAmodsq)/(norm(At[nt])-norm(At[nt-1])))-t0);
	break;
      }
    delete[] At;
  }
  else {
    lambda = pLP[0]->lambda;
    fwhm = pLP[0]->fwhm;
  }    
  //if(pLP_[0])  delete[] pLP_[0];
  pLP_[0] = new laser_parameters(I0*au_to_Wcm2, lambda*au_to_m, fwhm*au_to_s, t0*au_to_s, cep/pi, envelope);
}

// Function to analyse multi-colour electric field and write parameters found to specified laser_parameters object
void electric_field_1d::Analyse_MultiColour_Field(laser_parameters **&pLP_, dcmplx *Efield, bool findAll) {
  // both Ert and Erw should be up-to-date & consistant here
  double I0, lambda, fwhm, t0, cep;
  char envelope = 'c';  // fix this as cos^2 envelope for now (for QO)
  int nt_max;
  double nt_shift;
  double w_lower, w_upper, w_shift = pT->dt*(double)pT->Nt;
  double prefac = sqrt(2.0*pi)/(pT->Nt*pT->dt);
  double arg_max, arg_near;
  for(int nt=0; nt<pT->Nt; nt++) {
    fftwArrayT[nt][0] = real(Efield[nt]);
    fftwArrayT[nt][1] = imag(Efield[nt]); 
  }
  fftw_execute(time2freq); 
  if(findAll==false) {
    dcmplx *E_filter = new dcmplx[pT->Nt];
    double Enorm_max, Enorm_maxM1, Enorm_maxP1;
    for(int nc=0; nc<colours; nc++) {
      w_lower = pLP[nc]->w0-0.035*0.5;  w_upper = pLP[nc]->w0+0.035*0.5;
      for(int nt=0; nt<pT->Nt; nt++) {  
	if(fabs(pT->w[nt])<w_lower || fabs(pT->w[nt])>w_upper) { //  killed beyond one half-harmonic of 1300nm from w0
	  fftwArrayT[nt][0] = 0.0;
	  fftwArrayT[nt][1] = 0.0;
	}
      }
      fftw_execute(freq2time); 
      for(int nt=0; nt<pT->Nt; nt++)  E_filter[nt] = prefac*(fftwArrayT[nt][0]+i*fftwArrayT[nt][1]);
      Enorm_max = 0.0;
      for(int nt=0; nt<pT->Nt; nt++)
	//for(int nt=pT->Nt/4; nt<3*pT->Nt/4; nt++)  // assume within central 1/2 of the range
	if(norm(E_filter[nt])>Enorm_max)  {
	  nt_max = nt;
	  Enorm_max = norm(E_filter[nt]);
	}
      Enorm_maxM1 = norm(E_filter[nt_max-1]);
      Enorm_maxP1 = norm(E_filter[nt_max+1]);
      I0 = Enorm_max - 0.125*pow(Enorm_maxP1-Enorm_maxM1,2)/(Enorm_maxP1-2.0*Enorm_max+Enorm_maxM1); 
      nt_shift = -0.5*(Enorm_maxP1-Enorm_maxM1)/(Enorm_maxP1-2.0*Enorm_max+Enorm_maxM1); 
      t0 = pT->t[nt_max]+nt_shift*pT->dt;
      arg_max = atan2(imag(E_filter[nt_max]),real(E_filter[nt_max]));
      arg_near = (nt_shift>0.0) ? atan2(imag(E_filter[nt_max+1]),real(E_filter[nt_max+1])) : atan2(imag(E_filter[nt_max-1]),real(E_filter[nt_max-1]));
      if(arg_near-arg_max>pi)  arg_near-=2.0*pi;
      else if(arg_near-arg_max<-pi)  arg_near+=2.0*pi;
      cep = (1.0-fabs(nt_shift))*arg_max+fabs(nt_shift)*arg_near;
      if(cep>pi)  cep-=2.0*pi;
      else if(cep<-pi)  cep+=2.0*pi;
      lambda = pLP[nc]->lambda;
      fwhm = pLP[nc]->fwhm;
      //if(pLP_[nc])  delete[] pLP_[nc];
      pLP_[nc] = new laser_parameters(I0*au_to_Wcm2, lambda*au_to_m, fwhm*au_to_s, t0*au_to_s, cep/pi, envelope);
    }
    delete[] E_filter;
  }
  else {
    dcmplx *A_full = new dcmplx[pT->Nt];
    Vector_Potential(A_full, Efield);
    for(int nt=0; nt<pT->Nt; nt++) {
      fftwArrayT[nt][0] = real(A_full[nt]);
      fftwArrayT[nt][1] = imag(A_full[nt]); 
    }
    fftw_execute(time2freq); 
    for(int nt=0; nt<pT->Nt; nt++)  A_full[nt] = pT->dt/sqrt(2.0*pi)*(fftwArrayT[nt][0]+i*fftwArrayT[nt][1]) * exp(i*pT->w[nt]*w_shift);    
    double nw_shift, w_shift;  // not directly connected
    int nw_lower, nw_upper;
    double w0_temp;
    double Anorm_max, Anorm_maxM1, Anorm_maxP1, halfAmodsq;
    dcmplx *A_filter = new dcmplx[pT->Nt];
    for(int nc=0; nc<colours; nc++) {
      nw_lower = (int)((pLP[nc]->w0-0.035)/pT->dw);  nw_upper = (int)((pLP[nc]->w0+0.035)/pT->dw);
      Anorm_max = 0.0;
      for(int nw=nw_lower; nw<nw_upper; nw++) 
	if(norm(A_full[nw])>Anorm_max && norm(A_full[nw])>norm(A_full[nw-1]) && norm(A_full[nw])>norm(A_full[nw+1]))  {
	  nt_max = nw;
	  Anorm_max = norm(A_full[nw]);
	}
      Anorm_maxM1 = norm(A_full[nt_max-1]);
      Anorm_maxP1 = norm(A_full[nt_max+1]);
      Anorm_max = Anorm_max-0.125*pow(Anorm_maxP1-Anorm_maxM1,2)/(Anorm_maxP1-2.0*Anorm_max+Anorm_maxM1);
      nw_shift = -0.5*(Anorm_maxP1-Anorm_maxM1)/(Anorm_maxP1-2.0*Anorm_max+Anorm_maxM1); 
      w0_temp = pT->w[nt_max]+nw_shift*pT->dw;      
      lambda = 2.0*pi/w0_temp*c_au;
      w_lower = w0_temp-0.035*0.5;  w_upper = w0_temp+0.035*0.5;
      for(int nt=0; nt<pT->Nt; nt++) {  
	if(fabs(pT->w[nt])>=w_lower && fabs(pT->w[nt])<=w_upper)  {
	  fftwArrayT[nt][0] = real(A_full[nt]*exp(-i*pT->w[nt]*w_shift));  
	  fftwArrayT[nt][1] = imag(A_full[nt]*exp(-i*pT->w[nt]*w_shift));
	}
	else {  //  killed beyond one half-harmonic of 1300nm from w0_1
	  fftwArrayT[nt][0] = 0.0;
	  fftwArrayT[nt][1] = 0.0;
	}
      }
      fftw_execute(freq2time); 
      for(int nt=0; nt<pT->Nt; nt++)  A_filter[nt] = prefac*(fftwArrayT[nt][0]+i*fftwArrayT[nt][1]);
      Anorm_max = 0.0;
      //for(int nt=Nt/4; nt<Nt/4*3; nt++)  // assume within central 1/2 of the range
      for(int nt=0; nt<pT->Nt; nt++)
	if(norm(A_filter[nt])>Anorm_max)  {
	  nt_max = nt;
	  Anorm_max = norm(A_filter[nt]);
	}
      Anorm_maxM1 = norm(A_filter[nt_max-1]);
      Anorm_maxP1 = norm(A_filter[nt_max+1]);
      Anorm_max = Anorm_max-0.125*pow(Anorm_maxP1-Anorm_maxM1,2)/(Anorm_maxP1-2.0*Anorm_max+Anorm_maxM1);
      I0 = pow(w0_temp,2)*Anorm_max;
      nt_shift = -0.5*(Anorm_maxP1-Anorm_maxM1)/(Anorm_maxP1-2.0*Anorm_max+Anorm_maxM1); 
      t0 = pT->t[nt_max]+nt_shift*pT->dt;     
      arg_max = atan2(imag(A_filter[nt_max]),real(A_filter[nt_max])); 
      arg_near = (nt_shift>0.0) ? atan2(imag(A_filter[nt_max+1]),real(A_filter[nt_max+1])) : atan2(imag(A_filter[nt_max-1]),real(A_filter[nt_max-1]));
      if(arg_near-arg_max>pi)  arg_near-=2.0*pi;
      else if(arg_near-arg_max<-pi)  arg_near+=2.0*pi;
      cep = ((1.0-fabs(nt_shift))*arg_max+fabs(nt_shift)*arg_near)-pi/2.0;
      if(cep>pi)  cep-=2.0*pi;
      else if(cep<-pi)  cep+=2.0*pi;
      halfAmodsq = 0.5*Anorm_max; 
      for(int nt=nt_max+2; nt<pT->Nt; nt++)
	if(norm(A_filter[nt-1])>halfAmodsq && norm(A_filter[nt])<halfAmodsq) {
	  fwhm = 2.0*(pT->t[nt-1]+pT->dt*fabs((norm(A_filter[nt-1])-halfAmodsq)/(norm(A_filter[nt])-norm(A_filter[nt-1])))-t0);  
	  break;
	} 
      //if(pLP_[nc])  delete[] pLP_[nc];
      pLP_[nc] = new laser_parameters(I0*au_to_Wcm2, lambda*au_to_m, fwhm*au_to_s, t0*au_to_s, cep/pi, envelope);
    }
    delete[] A_full, A_filter;
  }
}

// Output temporal pulse in specified domain
void electric_field_1d::Output_Temporal(dcmplx *Etw, std::ofstream &OUTSTREAM, char RoIoS='r') {
  switch(RoIoS) {
  case 'r':
    for(int nt=0; nt<pT->Nt; nt++)
      OUTSTREAM << real(Etw[nt]) << '\t';
    OUTSTREAM << '\n';
    break;
  case 'i':
    for(int nt=0; nt<pT->Nt; nt++)
      OUTSTREAM << imag(Etw[nt]) << '\t';
    OUTSTREAM << '\n';
    break;
  case 'n':
    for(int nt=0; nt<pT->Nt; nt++)
      OUTSTREAM << norm(Etw[nt]) << '\t';
    OUTSTREAM << '\n';
    break;
  case 'p':
    for(int nt=0; nt<pT->Nt; nt++)
      OUTSTREAM << arg(Etw[nt]) << '\t';
    OUTSTREAM << '\n';
    break;
  default:
    std::cerr << "Error in real/imaginary/spectrum choice for output function";
  }
}

// Output radial beam at specified nt (t0=t[Nt/2]) for radially-and temporally-dependent array
void electric_field_1d::Output_Radial(dcmplx *Ertw, int nt, char ToF, std::ofstream &OUTSTREAM, char RoIoS='r') {
  switch(ToF) {
  case 't':
    switch(RoIoS) {
    case 'r':
      for(int nr=0; nr<pX->Nr; nr++) 
	OUTSTREAM << real(Ertw[Index2V(nr,nt,pT->Nt)]) << '\t';
      OUTSTREAM << '\n';
      break;
    case 'i':
      for(int nr=0; nr<pX->Nr; nr++) 
	OUTSTREAM << imag(Ertw[Index2V(nr,nt,pT->Nt)]) << '\t';
      OUTSTREAM << '\n';
      break;
    case 'n':
      for(int nr=0; nr<pX->Nr; nr++) 
	OUTSTREAM << norm(Ertw[Index2V(nr,nt,pT->Nt)]) << '\t';
      OUTSTREAM << '\n';
      break;
    case 'p':
      for(int nr=0; nr<pX->Nr; nr++) 
	OUTSTREAM << arg(Ertw[Index2V(nr,nt,pT->Nt)]) << '\t';
      OUTSTREAM << '\n';
      break;
    default:
      std::cerr << "Error in real/imaginary/spectrum choice for output function\n";
    }
    break;
  case 'w':
    switch(RoIoS) {
    case 'r':
      for(int nr=0; nr<pX->Nr; nr++) 
	OUTSTREAM << real(Ertw[Index2V(nt,nr,pX->Nr)]) << '\t';
      OUTSTREAM << '\n';
      break;
    case 'i':
      for(int nr=0; nr<pX->Nr; nr++) 
	OUTSTREAM << imag(Ertw[Index2V(nt,nr,pX->Nr)]) << '\t';
      OUTSTREAM << '\n';
      break;
    case 'n':
      for(int nr=0; nr<pX->Nr; nr++) 
	OUTSTREAM << norm(Ertw[Index2V(nt,nr,pX->Nr)]) << '\t';
      OUTSTREAM << '\n';
    break;
    case 'p':
      for(int nr=0; nr<pX->Nr; nr++) 
	OUTSTREAM << arg(Ertw[Index2V(nt,nr,pX->Nr)]) << '\t';
      OUTSTREAM << '\n';
      break;
    default:
      std::cerr << "Error in real/imaginary/spectrum choice for output function\n";
    }
    break;
  default:
    std::cerr << "Error in time/frequency choice for output function\n";   
  }
}

// Output radial-temporal pulse
void electric_field_1d::Output_RadialTemporal(dcmplx *Ertw, char ToF, std::ofstream &OUTSTREAM, char RoIoS='r', bool mpi=false, int n1step=1, int n2step=1) { // currently n1step, n2step must be powers of 2
  int N1, N2;
  switch(ToF) {
  case 't':
    if(mpi==true)  N1 = pX->Nr_mpi;
    else  N1 = pX->Nr;
    N2 = pT->Nt;
    break;
  case 'w':
    if(mpi==true)  N1 = pT->Nt_mpi;
    else  N1 = pT->Nt;
    N2 = pX->Nr;
    break;
  default:
    std::cerr << "Error in time/freq choice for output function\n";
  }
  switch(RoIoS) {
  case 'r':
    for(int n1=0; n1<N1; n1+=n1step) {
      for(int n2=0; n2<N2; n2+=n2step)
	OUTSTREAM << real(Ertw[Index2V(n1,n2,N2)]) << '\t';
      OUTSTREAM << '\n';
    }
    break;
  case 'i':
    for(int n1=0; n1<N1; n1+=n1step) {
      for(int n2=0; n2<N2; n2+=n2step)
	OUTSTREAM << imag(Ertw[Index2V(n1,n2,N2)]) << '\t';
      OUTSTREAM << '\n';
    }
    break;
  case 'n':
    for(int n1=0; n1<N1; n1+=n1step) {
      for(int n2=0; n2<N2; n2+=n2step)
	OUTSTREAM << norm(Ertw[Index2V(n1,n2,N2)]) << '\t';
      OUTSTREAM << '\n';
    }
    break;
  case 'p':
    for(int n1=0; n1<N1; n1+=n1step) {
      for(int n2=0; n2<N2; n2+=n2step)
	OUTSTREAM << arg(Ertw[Index2V(n1,n2,N2)]) << '\t';
      OUTSTREAM << '\n';
    }
    break;
  default:
    std::cerr << "Error in real/imaginary/spectrum choice for output function\n";
  }
}

// Function to de-assign memory for temporal-only functions
void electric_field_1d::Destroy_Time() {
  delete[] Et;
  delete[] Ew;
  if(FFTW) {
    fftw_free(fftwArrayT);  // don't force sequential as can't use MPI functions (destructors run after MPI_finalize)
    fftw_destroy_plan(time2freq);  fftw_destroy_plan(freq2time);
  }
}

// Function to de-assign memory for temporal-spatial functions (which ones depends on the class member "spatial")
void electric_field_1d::Destroy_SpaceTime() {
  Destroy_Time();
  delete[] Er, Ert, Erw;
  if(mpi==true)  delete[] E_mpi;
}

/***************************************************************************************************************/
// electric_field_2d class for two-dimensional electric fields, inherited from electric_field_1d:
class electric_field_2d : public electric_field_1d {
 private:
  /***************************************************************************************************************/
  // Private Class Function Prototypes:
  void Initialise_Time();
  void Make_Time(fftw_complex *, fftw_plan);
  void Temporal_Field(bool, bool);
  void Calculate_MaxValue();
  void Reset_Field();
  void Initialise_SpaceTime();
  void Make_SpaceTime(fftw_complex *, fftw_plan, fftw_plan);
  void Make_SpaceTime(fftw_complex *, fftw_plan, fftw_plan, double);
  void Gaussian_Beam(double, double, int);
  void Vector_Potential(dcmplx *, dcmplx *, dcmplx *, dcmplx *);
  void Analyse_TwoColour_Field(laser_parameters **&, dcmplx *, dcmplx *, bool);
  void Destroy_Time();
  void Destroy_SpaceTime();
 public:
  /***************************************************************************************************************/
  // Public Class Function Prototypes:
  electric_field_2d(temporal_axes &);
  electric_field_2d(temporal_axes &, laser_parameters **, int);
  electric_field_2d(temporal_axes &, laser_parameters **, fftw_complex *, fftw_plan, int);
  electric_field_2d(temporal_axes &, spatial_axes &, bool);
  electric_field_2d(temporal_axes &, spatial_axes &, laser_parameters **, int, bool);
  electric_field_2d(temporal_axes &, spatial_axes &, laser_parameters **, double, int, bool);
  electric_field_2d(temporal_axes &, spatial_axes &, laser_parameters **, fftw_complex *, fftw_plan, fftw_plan, int, bool);
  ~electric_field_2d();
  double Pulse_Energy(char); 
  void Transform_Field_FT();
  void Transform_Field_FT(fftw_complex *, fftw_plan);
  void Transform_Field_IFT();
  void Transform_Field_IFT(fftw_complex *, fftw_plan);
  void Update_Et_at_RZ(int);
  void Update_Ew_at_RZ(int);
  void Transpose_Overwrite(char);
  void Combine_Fields(char);
  void Analyse_Field(laser_parameters **&, bool);
  void Analyse_Field(laser_parameters **&, int, bool, bool);
  void Analyse_Field(laser_parameters **&, dcmplx *, dcmplx *, bool);
  

  void Analyse_Field(laser_parameters **&, dcmplx *, bool);
};

/***************************************************************************************************************/
// Class Functions:

// Constructor for to-be-calculated temporal fields Et,Ew,Et2,Ew2
electric_field_2d::electric_field_2d(temporal_axes &T_) : electric_field_1d(T_) {
  Dim = 2;
  Initialise_Time();
}

// Constructor for temporal fields Et,Ew,Et2,Ew2 calculated at creation with parameters LP
electric_field_2d::electric_field_2d(temporal_axes &T_, laser_parameters **pLP_, int colours_=1) : electric_field_1d(T_, pLP_, colours_) {
  Dim = 2;
  Initialise_Time();
  Make_Time(fftwArrayT, time2freq);
}

// Constructor as previous but takes in FFTW arguments
electric_field_2d::electric_field_2d(temporal_axes &T_, laser_parameters **pLP_, fftw_complex *fftwArrayT_, fftw_plan time2freq_, int colours_=1) : electric_field_1d(T_, pLP_, fftwArrayT_, time2freq_, colours_) {
  Dim = 2;
  Initialise_Time();
  Make_Time(fftwArrayT_, time2freq_);
}

// Constructor for to-be-calculated temporal fields Et,Ew,Et2,Ew2; spatial field Er; temporo-spatial fields Ert,Erw,Ert2,Erw2
electric_field_2d::electric_field_2d(temporal_axes &T_, spatial_axes &X_, bool mpi_=false) : electric_field_1d(T_, X_, mpi_) {
  Dim = 2;
  Initialise_SpaceTime();
}

// Constructor as above with fields defined via parameters L
electric_field_2d::electric_field_2d(temporal_axes &T_, spatial_axes &X_, laser_parameters **pLP_, int colours_=1, bool mpi_=false) : electric_field_1d(T_, X_, pLP_, colours_, mpi_) {
  Dim = 2;
  Initialise_SpaceTime();
  Make_SpaceTime(fftwArrayT, time2freq, freq2time);
}

// As above but for laser parameters defined at a given z=zP (ie. different from the focus, usually centre of gas jet)
electric_field_2d::electric_field_2d(temporal_axes &T_, spatial_axes &X_, laser_parameters **pLP_, double zP_, int colours_=1, bool mpi_=false) : electric_field_1d(T_, X_, pLP_, zP_, colours_, mpi_) {
  Dim = 2;
  Initialise_SpaceTime();
  Make_SpaceTime(fftwArrayT, time2freq, freq2time, zP_);
}

// Constructor as previous but takes in FFTW arguments
electric_field_2d::electric_field_2d(temporal_axes &T_, spatial_axes &X_, laser_parameters **pLP_, fftw_complex *fftwArrayT_, fftw_plan time2freq_, fftw_plan freq2time_, int colours_=1, bool mpi_=false) : electric_field_1d(T_, X_, pLP_, fftwArrayT_, time2freq_, freq2time_, colours_, mpi_) {
  Dim = 2;
  Initialise_SpaceTime();
  Make_SpaceTime(fftwArrayT_, time2freq_, freq2time_);
}

// Destructor
electric_field_2d::~electric_field_2d() {
  if(spatial==false)  Destroy_Time();
  else  Destroy_SpaceTime();
}

// Initialisation of variable memory for temporal-only
void electric_field_2d::Initialise_Time() {
  Et2 = new dcmplx[pT->Nt];
  Ew2 = new dcmplx[pT->Nt];  
}

// Calculation of temporal fields
void electric_field_2d::Make_Time(fftw_complex *fftwArrayT_, fftw_plan time2freq_) {
  Temporal_Field(true, true);
  for(int nt=0; nt<pT->Nt; nt++) {
    Et2[nt] *= pT->AB[nt]; // apply temporal absorbing boundary
    fftwArrayT_[nt][0] = real(Et2[nt]);
    fftwArrayT_[nt][1] = imag(Et2[nt]);
    //fftwArrayT_[nt][1] = 0.0;
  }
  fftw_execute(time2freq_);
  double prefac = pT->dt/sqrt(2.0*pi);
  for(int nt=0; nt<pT->Nt; nt++)
    Ew2[nt] = prefac * (fftwArrayT_[nt][0]+i*fftwArrayT_[nt][1]);
  Calculate_MaxValue();
}

// Function to create multicolour temporal field
void electric_field_2d::Temporal_Field(bool Adefn=false, bool Ereset=false) {
  if(Ereset==true)  Reset_Field();
  for(int nc=0;nc<colours; nc++)
    electric_field_1d::Temporal_Field(Et2, nc, 'y', true);
}

// Calculates temporal field maximum
void electric_field_2d::Calculate_MaxValue() {
  maxEt = 0;
  for(int nt=0; nt<pT->Nt; nt++)
    if(sqrt(norm(Et[nt])+norm(Et2[nt]))>maxEt)
      maxEt = sqrt(norm(Et[nt])+norm(Et2[nt]));
}

// Function to reset just Et2 to zero
void electric_field_2d::Reset_Field() {
  for(int nt=0; nt<pT->Nt; nt++)  {
    Et2[nt] = 0.0;  Ew2[nt] = 0.0;
  }
}

// Initialisation of variable memory for temporal-spatial (which ones depends on the class member "spatial")
void electric_field_2d::Initialise_SpaceTime() {
  Initialise_Time();
  Erw2 = new dcmplx[pX->Nr*pT->Nt];
  Ert2 = new dcmplx[pX->Nr*pT->Nt];
  E_mpi2 = (mpi==true) ? new dcmplx[pX->Nr_mpi*pT->Nt] : NULL; // could also be pX->Nr*pT->Nt_mpi
}

// Calculation of temporal-spatial fields (which ones depends on the class member "spatial")
void electric_field_2d::Make_SpaceTime(fftw_complex *fftwArrayT_, fftw_plan time2freq_, fftw_plan freq2time_) {
  double k;
  double pre_factor = pT->dt/sqrt(2.0*pi);
  for(int nt=0; nt<pT->Nt; nt++)
    for(int nr=0; nr<pX->Nr; nr++) 
      Erw2[Index2V(nt,nr,pX->Nr)] = 0.0;
  for(int nc=0;nc<colours; nc++) {
    Reset_Field();
    electric_field_1d::Temporal_Field(Et2, nc, 'y', true);
    for(int nt=0; nt<pT->Nt; nt++) {
      fftwArrayT_[nt][0] = real(Et2[nt]);
      fftwArrayT_[nt][1] = imag(Et2[nt]);
      //fftwArrayT_[nt][1] = 0.0;
    }
    fftw_execute(time2freq_);
    for(int nt=0; nt<pT->Nt; nt++) {
      Ew2[nt] = pre_factor * (fftwArrayT_[nt][0]+i*fftwArrayT_[nt][1]);
      k = pT->w[nt]/c_au; 
      electric_field_1d::Gaussian_Beam(k, pX->z0, nc);
      for(int nr=0; nr<pX->Nr; nr++)
	Erw2[Index2V(nt,nr,pX->Nr)] += Er[nr]*Ew2[nt];
    }
  }
  Transform_Field_IFT(fftwArrayT_, freq2time_); // does x-component aswell (should not change values(?))
  Make_Time(fftwArrayT_, time2freq_);  // calculates combined multicolour Et, Ew independent of space (on-axis temporal field)
}

// As above but for laser parameters defined at a given z=zP (ie. different from the focus, usually centre of gas jet)
void electric_field_2d::Make_SpaceTime(fftw_complex *fftwArrayT_, fftw_plan time2freq_, fftw_plan freq2time_, double zP) {
  bool diffract = true;
  if(zP<pX->z0 || zP>pX->z0+(pX->Nz-1)*pX->dz) {  // check validity of zP 
    std::cerr << "Error in choice of z position where laser parameters are defined, resetting to laser focus\n";
    zP = pX->z0;
    diffract = false;
  }
  double k;
  double pre_factor = pT->dt/sqrt(2.0*pi);
  for(int nt=0; nt<pT->Nt; nt++)
    for(int nr=0; nr<pX->Nr; nr++) 
      Erw2[Index2V(nt,nr,pX->Nr)] = 0.0;
  for(int nc=0;nc<colours; nc++) {
    Reset_Field();
    electric_field_1d::Temporal_Field(Et2, nc, 'y', true);
    for(int nt=0; nt<pT->Nt; nt++) {
      fftwArrayT_[nt][0] = real(Et2[nt]);
      fftwArrayT_[nt][1] = imag(Et2[nt]);
      //fftwArrayT_[nt][1] = 0.0;
    }
    fftw_execute(time2freq_);
    for(int nt=0; nt<pT->Nt; nt++) {
      Ew2[nt] = pre_factor * (fftwArrayT_[nt][0]+i*fftwArrayT_[nt][1]);
      k = pT->w[nt]/c_au; 
      electric_field_1d::Gaussian_Beam(k, pX->z0, nc, !diffract);
      for(int nr=0; nr<pX->Nr; nr++)
	Erw2[Index2V(nt,nr,pX->Nr)] += Er[nr]*Ew2[nt];
    }
  }
  if(diffract) {
    int mpi_rank = 0; 
    dcmplx *A = new dcmplx[pX->Nr], *B = new dcmplx[pX->Nr], *C = new dcmplx[pX->Nr], *D = new dcmplx[pX->Nr];
    // Diffract laser field backwards to start of propagation axis (frequency domain):
    for(int nz=(int)((zP-pX->z0)/pX->dz); nz>=0; nz--) {
#ifdef useMPI
      MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
      MPI_Scatter(this->Erw2, pX->Nr*pT->Nt_mpi, MPI_DOUBLE_COMPLEX, this->E_mpi2, pX->Nr*pT->Nt_mpi, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
      //Diffract_Field(*this, A, B, C, D, *pX, *pT, -pX->dz, mpi_rank);
      Diffract_Field(this->E_mpi2, A, B, C, D, *pX, *pT, -pX->dz, mpi_rank);
      MPI_Allgather(this->E_mpi2, pX->Nr*pT->Nt_mpi, MPI_DOUBLE_COMPLEX, this->Erw2, pX->Nr*pT->Nt_mpi, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
#else
      //Diffract_Field(*this, A, B, C, D, *pX, *pT, -pX->dz, mpi_rank);
      Diffract_Field(this->Erw2, A, B, C, D, *pX, *pT, -pX->dz, mpi_rank);
#endif
    }
    delete[] A, B, C, D;
  }
  Transform_Field_IFT(fftwArrayT_, freq2time_); // does x-component aswell (should not change values(?))
  Make_Time(fftwArrayT_, time2freq_);  // calculates combined multicolour Et, Ew independent of space (on-axis temporal field)
}

// Function to output the (temporo-spatial) pulse energy in Joules (removes intensity ambiguity in a.u.)
double electric_field_2d::Pulse_Energy(char method='n') {
  double Energy = 0.0;
  switch(method) { 
  case 'n':  // numerical integration
    for(int nt=0; nt<pT->Nt; nt++) {
      Energy += 0.5*(pow(real(Ert[Index2V(0,nt,pT->Nt)]),2)+pow(real(Ert2[Index2V(0,nt,pT->Nt)]),2))*pX->dr/8.0;  
      Energy += 0.5*(pow(real(Ert[Index2V(1,nt,pT->Nt)]),2)+pow(real(Ert2[Index2V(1,nt,pT->Nt)]),2))*pX->dr/8.0; 
      for(int nr=1; nr<pX->Nr; nr++)
	Energy += (pow(real(Ert[Index2V(nr,nt,pT->Nt)]),2)+pow(real(Ert2[Index2V(nr,nt,pT->Nt)]),2))*pX->r[nr];
    }
    Energy *= 2.0*pi*(au_to_Wcm2*1.0e4)*au_to_m*(pX->dr*au_to_m)*(pT->dt*au_to_s) * 2.0; // last 2.0 to convert to cycle averaged power (v accurate to even 0.5 cycle FWHMI)
    break;
  // experimental estimate removed (not needed since numerical calculation verified as accurate)
  default:
    std::cerr << "Error in method for calculation of pulse energy\n";
  }
  return Energy;
}

// FFT(W): Ert2 -> Erw2
void electric_field_2d::Transform_Field_FT() {
  electric_field_1d::Transform_Field_FT();
  double prefac = pT->dt/sqrt(2.0*pi);
  for(int nr=0; nr<pX->Nr; nr++) {
    for(int nt=0; nt<pT->Nt; nt++) {
      fftwArrayT[nt][0] = real(Ert2[Index2V(nr,nt,pT->Nt)]);
      fftwArrayT[nt][1] = imag(Ert2[Index2V(nr,nt,pT->Nt)]);
      //fftwArrayT[nt][1] = 0.0; // imaginary component of Ert2 should be zero for Et2=real(...) definition
    }
    fftw_execute(time2freq);
    for(int nt=0; nt<pT->Nt; nt++)
      Erw2[Index2V(nt,nr,pX->Nr)] = prefac * (fftwArrayT[nt][0]+i*fftwArrayT[nt][1]);
  }
}

// As above taking fftw array and plan as arguments
void electric_field_2d::Transform_Field_FT(fftw_complex *fftwArrayT_, fftw_plan time2freq_) {
  electric_field_1d::Transform_Field_FT(fftwArrayT_, time2freq_);
  double prefac = pT->dt/sqrt(2.0*pi);
  for(int nr=0; nr<pX->Nr; nr++) {
    for(int nt=0; nt<pT->Nt; nt++) {
      fftwArrayT_[nt][0] = real(Ert2[Index2V(nr,nt,pT->Nt)]);
      fftwArrayT_[nt][1] = imag(Ert2[Index2V(nr,nt,pT->Nt)]);
      //fftwArrayT_[nt][1] = 0.0; // imaginary component of Ert2 should be zero for Et2=real(...) definition
    }
    fftw_execute(time2freq_);
    for(int nt=0; nt<pT->Nt; nt++)
      Erw2[Index2V(nt,nr,pX->Nr)] = prefac * (fftwArrayT_[nt][0]+i*fftwArrayT_[nt][1]);
  }
}

// Inverse FFT(W): Erw2 -> Ert2
void electric_field_2d::Transform_Field_IFT() {
  electric_field_1d::Transform_Field_IFT();
  double prefac = sqrt(2.0*pi)/(pT->Nt*pT->dt);
  for(int nr=0; nr<pX->Nr; nr++) {
    for(int nt=0; nt<pT->Nt; nt++) {
      fftwArrayT[nt][0] = real(Erw2[Index2V(nt,nr,pX->Nr)]);
      fftwArrayT[nt][1] = imag(Erw2[Index2V(nt,nr,pX->Nr)]);
    }
    fftw_execute(freq2time);
    for(int nt=0; nt<pT->Nt; nt++)
      Ert2[Index2V(nr,nt,pT->Nt)] = prefac * (fftwArrayT[nt][0]+i*fftwArrayT[nt][1]);
    //Ert2[Index2V(nr,nt,pT->Nt)] = prefac * fftwArrayT[nt][0]; // imaginary component of Ert2 should be zero for Et2=real(...) definition
  }
}

// As above taking fftw array and plan as arguments
void electric_field_2d::Transform_Field_IFT(fftw_complex *fftwArrayT_, fftw_plan freq2time_) {
  electric_field_1d::Transform_Field_IFT(fftwArrayT_, freq2time_);
  double prefac = sqrt(2.0*pi)/(pT->Nt*pT->dt);
  for(int nr=0; nr<pX->Nr; nr++) {
    for(int nt=0; nt<pT->Nt; nt++) {
      fftwArrayT_[nt][0] = real(Erw2[Index2V(nt,nr,pX->Nr)]);
      fftwArrayT_[nt][1] = imag(Erw2[Index2V(nt,nr,pX->Nr)]);
    }
    fftw_execute(freq2time_);
    for(int nt=0; nt<pT->Nt; nt++)
      Ert2[Index2V(nr,nt,pT->Nt)] = prefac * (fftwArrayT_[nt][0]+i*fftwArrayT_[nt][1]);
    //Ert2[Index2V(nr,nt,pT->Nt)] = prefac * fftwArrayT_[nt][0]; // imaginary component of Ert2 should be zero for Et2=real(...) definition
  }
}

// Updates to Et2 spatial field Ert2 at current z and specified r
void electric_field_2d::Update_Et_at_RZ(int nr) {
  electric_field_1d::Update_Et_at_RZ(nr);
  for(int nt=0; nt<pT->Nt; nt++)
    Et2[nt] = Ert2[Index2V(nr,nt,pT->Nt)];
}

// Updates to Ew2 spatial field Erw2 at current z and specified r
void electric_field_2d::Update_Ew_at_RZ(int nr) {
  electric_field_1d::Update_Ew_at_RZ(nr);
  for(int nt=0; nt<pT->Nt; nt++)
    Ew2[nt] = Erw2[Index2V(nr,nt,pT->Nt)];
}

// Function to copy-and-transpose one temporo-spatial field onto the other (use for harmonic field where Et2 not needed during propagation)
void electric_field_2d::Transpose_Overwrite(char direction) {
  electric_field_1d::Transpose_Overwrite(direction);
  if(spatial==true && direction=='f')
    for(int nr=0; nr<pX->Nr; nr++)
      for(int nt=0; nt<pT->Nt; nt++)
	Erw2[Index2V(nt,nr,pX->Nr)] = Ert2[Index2V(nr,nt,pT->Nt)];
  else if(spatial==true && direction=='b')
    for(int nt=0; nt<pT->Nt; nt++)
      for(int nr=0; nr<pX->Nr; nr++)
	Ert2[Index2V(nr,nt,pT->Nt)] = Erw2[Index2V(nt,nr,pX->Nr)];
}

// Function to add one temporo-spatial field onto the other (specified by destination) - used in harmonic propagation
void electric_field_2d::Combine_Fields(char destination) {
  electric_field_1d::Combine_Fields(destination);
  if(spatial==true && destination=='t')
    for(int nr=0; nr<pX->Nr; nr++)
      for(int nt=0; nt<pT->Nt; nt++)
	Ert2[Index2V(nr,nt,pT->Nt)] += Erw2[Index2V(nt,nr,pX->Nr)];
  else if(spatial==true && destination=='w')
    for(int nt=0; nt<pT->Nt; nt++)
      for(int nr=0; nr<pX->Nr; nr++)
	Erw2[Index2V(nt,nr,pX->Nr)] += Ert2[Index2V(nr,nt,pT->Nt)];
}

// Function to calculate 2D vector potential from temporal electric fields
void electric_field_2d::Vector_Potential(dcmplx *At, dcmplx *At2, dcmplx *Efield, dcmplx *Efield2) {
  electric_field_1d::Vector_Potential(At, Efield);
  At2[0] = 0.0;
  At2[1] = -pT->dt*0.5*(Efield2[0]+Efield2[1]);
  for(int nt=2; nt<pT->Nt; nt++)  At2[nt] = At2[nt-2]-pT->dt*(Efield2[nt-2]+4.0*Efield2[nt-1]+Efield2[nt])/3.0;
}

// Governing interface to analyse electric fields and write parameters found to specified laser_parameters object for temporal-only fields
void electric_field_2d::Analyse_Field(laser_parameters **&pLP_, bool findAll) {
  // currently only consider 2 colours in the non-parallel case
  Analyse_TwoColour_Field(pLP_, Et, Et2, findAll); 
}

// As above for temporo-spatial fields
void electric_field_2d::Analyse_Field(laser_parameters **&pLP_, int nr, bool findAll, bool mpi) {
  // currently only consider 2 colours in the non-parallel case
  (mpi==true) ? Analyse_TwoColour_Field(pLP_, &E_mpi[Index2V(nr,0,pT->Nt)], &E_mpi2[Index2V(nr,0,pT->Nt)], findAll) : Analyse_TwoColour_Field(pLP_, &Ert[Index2V(nr,0,pT->Nt)], &Ert2[Index2V(nr,0,pT->Nt)], findAll); 
}

// As above but taking in fields explicitly
void electric_field_2d::Analyse_Field(laser_parameters **&pLP_, dcmplx *EfieldX, dcmplx *EfieldY, bool findAll) {
  // currently only consider 2 colours in the non-parallel case
  //(mpi==true) ? Analyse_TwoColour_Field(pLP_, EfieldX, EfieldY, findAll) : Analyse_TwoColour_Field(pLP_, EfieldX, EfieldY, findAll); 
  Analyse_TwoColour_Field(pLP_, EfieldX, EfieldY, findAll);
}

// As above but taking fields in explicitly
void electric_field_2d::Analyse_Field(laser_parameters **&pLP_, dcmplx *Efield, bool findAll) {
  electric_field_1d::Analyse_SingleColour_Field(pLP_, Efield, findAll);
  //delete[] pLP_[1];
  pLP_[1] = new laser_parameters(0.0, 0.5*pLP_[0]->lambda*au_to_m, pLP_[0]->fwhm*au_to_s, pLP_[0]->t0*au_to_s, pLP_[0]->cep/pi, pLP_[0]->envelope);
}

// Function to analyse TWO-colour, two-dimensional electric field and write parameters found to specified laser_parameters object
void electric_field_2d::Analyse_TwoColour_Field(laser_parameters **&pLP_, dcmplx *Efield, dcmplx *Efield2, bool findAll) {
  // use class member pLP to get polAxis for both colours
  // need to store polAxis in new pLP_ then!!!!
  double sin_phi0 = sin(pLP[0]->phi), cos_phi0 = cos(pLP[0]->phi);
  double sin_phi1 = sin(pLP[1]->phi), cos_phi1 = cos(pLP[1]->phi);
  double inv_sin_diff = 1.0/sin(pLP[0]->phi-pLP[1]->phi);
  dcmplx *E_split = new dcmplx[2*pT->Nt];
  for(int nt=0; nt<pT->Nt; nt++) {
    E_split[nt] = (-Efield[nt]*sin_phi1+Efield2[nt]*cos_phi1)*inv_sin_diff;
    E_split[nt+pT->Nt] = (Efield[nt]*sin_phi0-Efield2[nt]*cos_phi0)*inv_sin_diff;
  } 
  double I0, lambda, fwhm, t0, cep, phi;
  char envelope = 'c';  // fix this as cos^2 envelope for now (for QO)
  int nt_max, nt_shift;
  double Enorm_max, Enorm_maxM1, Enorm_maxP1;
  double arg_max, arg_near;
  double t_minM, t_minP;
  dcmplx *A_split = NULL;
  if(findAll==true) {
    A_split = new dcmplx[2*pT->Nt];
    Vector_Potential(&A_split[0], &A_split[pT->Nt], Efield, Efield2);
    dcmplx A_temp;
    for(int nt=0; nt<pT->Nt; nt++) {
      A_temp = A_split[nt];
      A_split[nt] = (-A_temp*sin_phi1+A_split[nt+pT->Nt]*cos_phi1)*inv_sin_diff;
      A_split[nt+pT->Nt] = (A_temp*sin_phi0-A_split[nt+pT->Nt]*cos_phi0)*inv_sin_diff;
    }
  }
  for(int nax=0; nax<2; nax++) {
    Enorm_max = 0.0;
    for(int nt=0; nt<pT->Nt; nt++)
      //for(int nt=pT->Nt/4; nt<3*pT->Nt/4; nt++)  // assume within central 1/2 of the range
      if(norm(E_split[nt+nax*pT->Nt])>Enorm_max)  {
	nt_max = nt;
	Enorm_max = norm(E_split[nt+nax*pT->Nt]);
      }
    Enorm_maxM1 = norm(E_split[nt_max-1+nax*pT->Nt]);
    Enorm_maxP1 = norm(E_split[nt_max+1+nax*pT->Nt]);
    I0 = Enorm_max - 0.125*pow(Enorm_maxP1-Enorm_maxM1,2)/(Enorm_maxP1-2.0*Enorm_max+Enorm_maxM1); 
    nt_shift = -0.5*(Enorm_maxP1-Enorm_maxM1)/(Enorm_maxP1-2.0*Enorm_max+Enorm_maxM1); 
    t0 = pT->t[nt_max]+nt_shift*pT->dt;
    arg_max = atan2(imag(E_split[nt_max+nax*pT->Nt]),real(E_split[nt_max+nax*pT->Nt]));
    arg_near = (nt_shift>0.0) ? atan2(imag(E_split[nt_max+1+nax*pT->Nt]),real(E_split[nt_max+1+nax*pT->Nt])) : atan2(imag(E_split[nt_max-1+nax*pT->Nt]),real(E_split[nt_max-1+nax*pT->Nt]));
    if(arg_near-arg_max>pi)  arg_near-=2.0*pi;
    else if(arg_near-arg_max<-pi)  arg_near+=2.0*pi;
    cep = (1.0-fabs(nt_shift))*arg_max+fabs(nt_shift)*arg_near;
    if(cep>pi)  cep-=2.0*pi;
    else if(cep<-pi)  cep+=2.0*pi;
    if(findAll==true) {
      t_minM=0.0, t_minP=0.0;
      for(int nt=nt_max+2; nt<pT->Nt; nt++)  
	if(real(A_split[nt-1+nax*pT->Nt]*exp(-i*cep))<0.0&&real(A_split[nt+nax*pT->Nt]*exp(-i*cep))>0.0 || real(A_split[nt-1+nax*pT->Nt]*exp(-i*cep))>0.0&&real(A_split[nt+nax*pT->Nt]*exp(-i*cep))<0.0) {
	  t_minP = pT->t[nt-1]+pT->dt*fabs(real(A_split[nt-1+nax*pT->Nt]*exp(-i*cep))/real((A_split[nt+nax*pT->Nt]-A_split[nt-1+nax*pT->Nt])*exp(-i*cep)));  
	  break;
	}
      lambda = 2.0*(t_minP-t0)*c_au;
      double halfAmodsq = 0.5 * (norm(A_split[nt_max+nax*pT->Nt])-0.125*pow(norm(A_split[nt_max+1+nax*pT->Nt])-norm(A_split[nt_max-1+nax*pT->Nt]),2)/(norm(A_split[nt_max+1+nax*pT->Nt])-2.0*norm(A_split[nt_max+nax*pT->Nt])+norm(A_split[nt_max-1+nax*pT->Nt]))); 
      for(int nt=nt_max+2; nt<pT->Nt; nt++)
	if(norm(A_split[nt-1+nax*pT->Nt])>halfAmodsq && norm(A_split[nt+nax*pT->Nt])<halfAmodsq) {
	  fwhm = 2.0*(pT->t[nt-1]+pT->dt*fabs((norm(A_split[nt-1+nax*pT->Nt])-halfAmodsq)/(norm(A_split[nt+nax*pT->Nt])-norm(A_split[nt-1+nax*pT->Nt])))-t0);
	  break;
	}
    }
    else {
      lambda = pLP[nax]->lambda;
      fwhm = pLP[nax]->fwhm;
    } 
    //if(pLP_[nax])  delete[] pLP_[nax];
    pLP_[nax] = new laser_parameters(I0*au_to_Wcm2, lambda*au_to_m, fwhm*au_to_s, t0*au_to_s, cep/pi, envelope, pLP[nax]->phi/pi);
  }
  delete[] E_split;
  if(findAll==true)  delete[] A_split;
}

// Function to de-assign memory for temporal-only functions
void electric_field_2d::Destroy_Time() {
  delete[] Et2;
  delete[] Ew2;
}

// Function to de-assign memory for temporal-spatial functions (which ones depends on the class member "spatial")
void electric_field_2d::Destroy_SpaceTime() {
  Destroy_Time();
  Ert2, Erw2;
  if(mpi==true)  delete[] E_mpi2;
}

/***************************************************************************************************************/
// Universal functions on electric_field objects

// Function to calculate (re)absorption of electric field by the neutral gas in freq domain (so doesn't include depletion of ground state with time)
void Absorption_By_Gas(dcmplx *Erw, gas_jet &Jet, int nz, spatial_axes &X, temporal_axes &T, int mpi_rank) {
  double abs_factor;
  for(int nt=0; nt<T.Nt_mpi; nt++) {
    abs_factor = exp(-Jet.nd[nz]*Jet.absorb[mpi_rank*T.Nt_mpi+nt]);
    for(int nr=0; nr<X.Nr; nr++)
      Erw[Index2V(nt,nr,X.Nr)] *= abs_factor;
  }
}

// As above taking electric_field class as argument
void Absorption_By_Gas(electric_field &Efield, gas_jet &Jet, int nz, spatial_axes &X, temporal_axes &T, int mpi_rank) {
#ifdef useMPI  
  Absorption_By_Gas(Efield.E_mpi, Jet, nz, X, T, mpi_rank);
  if(Efield.Dim==2)  Absorption_By_Gas(Efield.E_mpi2, Jet, nz, X, T, mpi_rank);
#else
  Absorption_By_Gas(Efield.Erw, Jet, nz, X, T, mpi_rank);
  if(Efield.Dim==2)  Absorption_By_Gas(Efield.Erw2, Jet, nz, X, T, mpi_rank);
#endif
}

/***************************************************************************************************************/
#endif
