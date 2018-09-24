/***************************************************************************************************************/
// Header file to contain dipole_response class for use in Propagator.cpp 

/***************************************************************************************************************/
// Notes:
// Removed constructor which didn't take an electric field object as an argument - what was this even for?

/***************************************************************************************************************/
#ifndef DIPOLE_RESPONSE_H
#define DIPOLE_RESPONSE_H

#include<iostream>
#include<fstream>
#include<Environment.h>
#include<Misc_Functions.h>
#include<Gas_Jet.h>
#include<Electric_Field.h>
//#include<SFA_Response.h> - must be included below (after dipole_response class)
//#include<QO_Response.h> - must be included below (after dipole_response class)

/***************************************************************************************************************/
// Abstract base class for general dipole response (doesn't include atomic density - only population fraction):
class dipole_response {
 private:
  /***************************************************************************************************************/
  // Private Class Members:
  char method;
  /***************************************************************************************************************/
  // Private Class Function Prototypes:
  void Initialise();
 protected:
  /***************************************************************************************************************/
  // Protected Class Members:
  temporal_axes *pT;
  electric_field *pLaser;
  gas_jet *pJet;
 public:
  /***************************************************************************************************************/
  // Public Class Members:
  dcmplx *DA, *DA2;  // (freq domain)
  /***************************************************************************************************************/
  // Public Class Function Prototypes:
  dipole_response(temporal_axes &, electric_field &, gas_jet &, char);
  ~dipole_response();
  virtual void Dipole_Acceleration() = 0; // pure virtual function (no implementation)
  virtual void Dipole_Acceleration(int, bool) = 0; 
  virtual void Dipole_Acceleration(dcmplx *) = 0; 
  virtual void Dipole_Acceleration(dcmplx *, int) = 0; 
  virtual void Dipole_Acceleration(dcmplx *, dcmplx *) = 0; 
  virtual void Dipole_Acceleration(dcmplx *, dcmplx *, int) = 0; 
};

/***************************************************************************************************************/
// Class Functions:
dipole_response::dipole_response(temporal_axes &T_, electric_field &Laser_, gas_jet &Jet_, char method_='S') {
  pT = &T_;
  pLaser = &Laser_;
  pJet = &Jet_;
  method = method_;
  Initialise();
}

dipole_response::~dipole_response() {
  delete[] DA;
  if(pLaser->Dim==2)  delete[] DA2;
}

void dipole_response::Initialise() {
  DA = new dcmplx[pT->Nt];
  DA2 = (pLaser->Dim==2) ? new dcmplx[pT->Nt] : NULL;
}

/***************************************************************************************************************/
#include<SFA_Response.h>
#include<QO_Response.h>

/***************************************************************************************************************/
// Universal functions on dipole_response objects

// Function to assign method to calculate dipole response taking method as argument
void Assign_Method(dipole_response *&Dipole, electric_field &Laser, gas_jet &Jet, temporal_axes &T, char method) {
  if(method=='S')  
    Dipole = (Laser.Dim==2) ? new sfa_response_2d(T, Laser, Jet) : new sfa_response_1d(T, Laser, Jet);
  else if(method=='Q')  
    Dipole = (Laser.Dim==2) ? new qo_response_2d(T, Laser, Jet) : new qo_response_1d(T, Laser, Jet); 

/*     if(Laser.Dim==2)  */
/*       Dipole = new qo_response_2d(T, Laser, Jet); */
/*     else if(Laser.Dim==1)  */
/*       Dipole = new qo_response_1d(T, Laser, Jet); */

  else {
    std::cerr << "Error in choice of method for calculating dipole response\n";
#ifdef useMPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(1);
  }
}

// Function to assign method to calculate dipole response taking filestream containing method as argument
void Assign_Method(dipole_response *&Dipole, electric_field &Laser, gas_jet &Jet, temporal_axes &T, std::ifstream &INSTREAM) {
  char method = Read_Variable<char>(INSTREAM, "DRMethod");
  if(method=='S')  
    Dipole = (Laser.Dim==2) ? new sfa_response_2d(T, Laser, Jet, INSTREAM) : new sfa_response_1d(T, Laser, Jet, INSTREAM);
  else if(method=='Q')  
    Dipole = (Laser.Dim==2) ? new qo_response_2d(T, Laser, Jet, INSTREAM) : new qo_response_1d(T, Laser, Jet, INSTREAM);

/*     if(Laser.Dim==2) */
/*       Dipole = new qo_response_2d(T, Laser, Jet, INSTREAM); */
/*     else if(Laser.Dim==1) */
/*       Dipole = new qo_response_1d(T, Laser, Jet, INSTREAM);  */

  else {
    std::cerr << "Error in choice of method for calculating dipole response\n";
#ifdef useMPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(1);
  }
}

// Function to generate new harmonics at this z-position
void Generate_Harmonics(electric_field &Hfield, electric_field &Efield, double minI, dipole_response *&Dipole, double nd, spatial_axes &X, temporal_axes &T) {
  int nmax;
  double Emax, shift, Elocalmax, Emag;
  double minE = sqrt(minI*Wcm2_to_au);
  dcmplx *Ex, *Ey, *Hx, *Hy;
  bool mpi_active;
#ifdef useMPI  
  Ex = Efield.E_mpi;  Ey = Efield.E_mpi2;
  Hx = Hfield.E_mpi;  Hy = Hfield.E_mpi2;
  mpi_active = true;
#else
  Ex = Efield.Ert;  Ey = Efield.Ert2;
  Hx = Hfield.Ert;  Hy = Hfield.Ert2;
  mpi_active = false;
#endif
  for(int nr=0; nr<X.Nr_mpi; nr++) {
    Elocalmax=0.0;
    for(int nt=T.Nt/4; nt<3*T.Nt/4; nt++) {
      Emag = (Efield.Dim==2) ? sqrt(pow(abs(Ex[Index2V(nr,nt,T.Nt)]),2)+pow(abs(Ey[Index2V(nr,nt,T.Nt)]),2)) : abs(Ex[Index2V(nr,nt,T.Nt)]);
      if(Emag>Elocalmax)  Elocalmax = Emag;
    }
    if(Elocalmax>minE) {
      /* (Efield.Dim==2) ? Dipole->Dipole_Acceleration(&Ex[Index2V(nr,0,T.Nt)], &Ey[Index2V(nr,0,T.Nt)]) : Dipole->Dipole_Acceleration(&Ex[Index2V(nr,0,T.Nt)]); */
      Dipole->Dipole_Acceleration(nr, mpi_active);  // changed to sync up SFA and QO methods
      for(int nt=0; nt<T.Nt; nt++) {
	Hx[Index2V(nr,nt,T.Nt)] = nd * Dipole->DA[nt];
	/* if(Efield.Dim==2)  Hy[Index2V(nr,nt,T.Nt)] = nd * Dipole->DA2[nt]; */
	if(Hfield.Dim==2)  Hy[Index2V(nr,nt,T.Nt)] = nd * Dipole->DA2[nt];
      }
    }
    else
      for(int nt=0; nt<T.Nt; nt++) {
	Hx[Index2V(nr,nt,T.Nt)] = 0.0;
	/* if(Efield.Dim==2)  Hy[Index2V(nr,nt,T.Nt)] = 0.0; */
	if(Hfield.Dim==2)  Hy[Index2V(nr,nt,T.Nt)] = 0.0;
      }
  }
}

/***************************************************************************************************************/
#endif


