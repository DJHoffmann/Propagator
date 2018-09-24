/***************************************************************************************************************/
// Header file to contain sfa_response_1/2d classes for use in Propagator.cpp 

/***************************************************************************************************************/
// Notes:
// Removed constructor which didn't take an electric field object as an argument - what was this even for?

/***************************************************************************************************************/
#ifndef SFA_RESPONSE_H
#define SFA_RESPONSE_H

#include<fftw3.h>
#ifdef useMPI
#ifdef oldMPI
#undef SEEK_SET  // undef needed due to a bug with MPI-2, bug fixed in openmpi >= 1.3.2
#undef SEEK_END
#undef SEEK_CUR
#endif
#include <mpi.h>
#endif
#include<Constants_Conversions.h>
#include<Misc_Functions.h>
#include<Input.h>
#include<Environment.h>
#include<Electric_Field.h>
#include<Dipole_Response.h>

/***************************************************************************************************************/
// Class to create one-dimensional dipole response via SFA method:
class sfa_response_1d : public dipole_response {  
 protected:    
  /***************************************************************************************************************/
  // Protected Class Members:
  int min_ntau, max_ntau;                
  double p, Sp;
  dcmplx dtme[2], pre_dtme;
  double *A, *intA, *intAsq, *population;
  dcmplx *spread;
  dcmplx *D;
  fftw_complex *fftwArrayT;
  fftw_plan time2freq;
  /***************************************************************************************************************/
  // Protected Class Function Prototypes:
  void Initialise();
  void Set_MinMaxTau();
  void Set_MinMaxTau(std::ifstream &);
  void Calculate_A(dcmplx *);
  void Calculate_intA();
  void Calculate_intAsq();
  void Wavepacket_Spreading();
  void DTME_Prefactor();
  dcmplx Dipole_Transition_Matrix_Element(double);
  dcmplx Dipole_Transition_Matrix_Element(double, double);
  double Tau_Smooth(int, int, int);
  void Dipole_Moment(dcmplx *, char);
  void Dipole_Moment_to_Acceleration();
 public:
  /***************************************************************************************************************/
  // Public Class Function Prototypes:
  sfa_response_1d(temporal_axes &T_, electric_field &Laser_, gas_jet &Jet_);
  sfa_response_1d(temporal_axes &T_, electric_field &Laser_, gas_jet &Jet_, std::ifstream &);
  ~sfa_response_1d();
  void Dipole_Acceleration(); 
  void Dipole_Acceleration(int, bool); 
  void Dipole_Acceleration(dcmplx *); 
  void Dipole_Acceleration(dcmplx *, int); 
  void Dipole_Acceleration(dcmplx *, dcmplx *); 
  void Dipole_Acceleration(dcmplx *, dcmplx *, int); 
};

/***************************************************************************************************************/
// Class Functions:
  
// Constructor
sfa_response_1d::sfa_response_1d(temporal_axes &T_, electric_field &Laser_, gas_jet &Jet_) : dipole_response(T_, Laser_, Jet_, 'S') {
  Initialise();
  Set_MinMaxTau();
  DTME_Prefactor();
  Wavepacket_Spreading();
}  

// Constructor
sfa_response_1d::sfa_response_1d(temporal_axes &T_, electric_field &Laser_, gas_jet &Jet_, std::ifstream &INSTREAM) : dipole_response(T_, Laser_, Jet_, 'S') {
  Initialise();
  Set_MinMaxTau(INSTREAM);
  DTME_Prefactor();
  Wavepacket_Spreading();
}  

// Destructor
sfa_response_1d::~sfa_response_1d() {
  delete[] A;
  delete[] intA;
  delete[] intAsq;
  delete[] population;
  delete[] spread;
  delete[] D;
  fftw_free(fftwArrayT);  // don't force sequential as can't use MPI functions (destructors run after MPI_finalize)
  fftw_destroy_plan(time2freq);
}

// Initialise memory
void sfa_response_1d::Initialise() {
  A = new double[pT->Nt];
  intA = new double[pT->Nt];
  intAsq = new double[pT->Nt];
  population = new double[pT->Nt];
  spread = new dcmplx[pT->Nt];
  D = new dcmplx[pT->Nt];
#ifdef useMPI  
  int rank = 0, size = 1;  // force sequential for thread-safety 
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);  MPI_Comm_size(MPI_COMM_WORLD,&size);
  for(int n=0; n<size; n++) {
    if(rank==n) {
      fftwArrayT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*pT->Nt);
      time2freq = fftw_plan_dft_1d(pT->Nt, fftwArrayT, fftwArrayT, FFTW_FORWARD, FFTW_MEASURE);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#else
  fftwArrayT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*pT->Nt);
  time2freq = fftw_plan_dft_1d(pT->Nt, fftwArrayT, fftwArrayT, FFTW_FORWARD, FFTW_MEASURE);
#endif
}

// Function to set max_ntau for SFA calculations
void sfa_response_1d::Set_MinMaxTau() {
  double max_tau_cycles = 2.0;
  min_ntau = 1; // must be >=1
  max_ntau = 0;
  for(int nc=0; nc<pLaser->colours; nc++)
    if((int)(max_tau_cycles*(2.0*pi/pLaser->pLP[nc]->w0)/pT->dt)>max_ntau)
      max_ntau = (int)(max_tau_cycles*(2.0*pi/pLaser->pLP[nc]->w0)/pT->dt);
}

// Function to set max_ntau for SFA calculations using values taken in from file
void sfa_response_1d::Set_MinMaxTau(std::ifstream &INSTREAM) {
  double max_tau_cycles = Read_Variable<double>(INSTREAM, "SFAcycles");
  min_ntau = 1; // must be >=1
  max_ntau = 0;
  for(int nc=0; nc<pLaser->colours; nc++)
    if((int)(max_tau_cycles*(2.0*pi/pLaser->pLP[nc]->w0)/pT->dt)>max_ntau)
      max_ntau = (int)(max_tau_cycles*(2.0*pi/pLaser->pLP[nc]->w0)/pT->dt);
}

// Calculate A
void sfa_response_1d::Calculate_A(dcmplx *Et) {
  A[0] = 0.0;
  for(int nt=1; nt<pT->Nt-1; nt++)
    A[nt] = A[nt-1] - real(Et[nt])*pT->dt;
}

// Calculate intA
void sfa_response_1d::Calculate_intA() {
  intA[0] = 0.0;
  for(int nt=1; nt<pT->Nt-1; nt++)
    intA[nt] = intA[nt-1] + A[nt]*pT->dt;
}

// Calculate intAsq
void sfa_response_1d::Calculate_intAsq() {
  intAsq[0] = 0.0;
  for(int nt=1; nt<pT->Nt-1; nt++)
    intAsq[nt] = intAsq[nt-1] + A[nt]*A[nt]*pT->dt;
}

// Wavepacket spreading
void sfa_response_1d::Wavepacket_Spreading() {
  double eps = 0.0001;
  for(int nt=0; nt<pT->Nt; nt++) 
    spread[nt] = sqrt(pow(2.0*pi/(eps+i*(double)nt*pT->dt), 3));
  //spread[nt] = pow(2.0*pi/(eps+i*(double)nt*pT->dt), 1.5); // don't use pT->t[nt] as representing tau
}

// Dipole transition matrix element
void sfa_response_1d::DTME_Prefactor() {
  if(pJet->n==1)
    pre_dtme = i*8.0*sqrt(2.0)*pow(2.0*pJet->Ip, 1.25)/pi;
  else if(pJet->n==2)
    if(pJet->l==0)
      pre_dtme = i*pow(2.0, 4.5)*pow(2.0*pJet->Ip, 1.25)/pi;
    else if(pJet->l==1)
      pre_dtme = pow(2.0, 3.5)*pow(2.0*pJet->Ip, 1.75)/pi; 
  else if(pJet->n==3)
    if(pJet->l==0)
      pre_dtme = i*pow(2.0,3.5)*pow(2.0*pJet->Ip, 1.25)/pi;
    else
      pre_dtme = 16.0*sqrt(3.0)*pow(2.0*pJet->Ip, 1.75)/pi;
  else
    pre_dtme = 0.0;    
}

// Dipole transition matrix element
dcmplx sfa_response_1d::Dipole_Transition_Matrix_Element(double k) {
  dcmplx element;
  double kSign = (k<0) ? -1.0 : 1.0;
  if(pJet->n==1)
    element = pre_dtme*k/pow(k*k+2.0*pJet->Ip, 3);
  else if(pJet->n==2)
    if(pJet->l==0)
      element = pre_dtme*k*(k*k-4.0*pJet->Ip)/pow(k*k+2.0*pJet->Ip, 4);
    else if(pJet->l==1)
      element = pre_dtme*k*(5.0*k*k-2.0*pJet->Ip)/pow(k*k+2.0*pJet->Ip, 4);  // assuming p-state aligned along (1d) laser polarisation axis
  else if(pJet->n==3)
    if(pJet->l==0)
      element = pre_dtme*(3.0*pow(k, 4)-18.0*k*k*2.0*pJet->Ip+11.0*pow(2.0*pJet->Ip, 2))/pow(k*k+2.0*pJet->Ip, 5);
    else if(pJet->l==1)
      element = pre_dtme*(2.0*k*k*(3.0*k*k-10.0*pJet->Ip)-(pow(k, 4)-pow(2.0*pJet->Ip, 2)))/pow(k*k+2.0*pJet->Ip, 5);
    else
      element = pre_dtme*0.0;
  else 
    element = pre_dtme*0.0;     
  return element;
}

// Smoothing function for application of mintau and maxtau
double sfa_response_1d::Tau_Smooth(int ntau, int min_ntau, int max_ntau) {
  double smooth;
  /*
  double factor_low = 0.95;
  double factor_high = 0.95;
  if(ntau<min_ntau)
    smooth = 0.0;
  else if(ntau<factor_low*min_ntau+(1.0-factor_low)*max_ntau) 
    smooth = pow(sin(1.0/(1.0-factor_low)*(tau-min_ntau)/(max_ntau-min_ntau)*pi/2), 2.0);
  else if(ntau<(1.0-factor_high)*min_ntau+factor_high*max_ntau)  
    smooth = 1.0;
  else if(ntau<max_ntau) 
    smooth = pow(cos(1.0/(1.0-factor_high)*(tau-(1.0-factor_high)*min_ntau-factor_high*max_ntau)/(max_ntau-min_ntau)*pi/2),2.0);
  else 
    smooth = 0.0;
  return smooth;
  */
  if(ntau>0.8*max_ntau)
    smooth = exp(-pow(5.0*(ntau-0.8*max_ntau)/(double)max_ntau, 2));
  return smooth;
}

// Dipole moment in time-domain
void sfa_response_1d::Dipole_Moment(dcmplx *Et, char method='b') {
  double smooth;
  Calculate_A(Et);
  Calculate_intA();
  Calculate_intAsq(); 
  pJet->ADK_Ionisation(population, Et);
  switch(method) {
    case 'b':  // backward temporal integration
      for(int ntr=0; ntr<min_ntau; ntr++)
        D[ntr] = 0.0;
      for(int ntr=min_ntau; ntr<pT->Nt; ntr++) { 
        D[ntr] = 0.0;
        for(int nti=ntr-min_ntau; nti>=std::max(ntr-max_ntau,0); nti--) { 
	  // what about a "if E>Emin" statement around these calculations?
          p = (intA[nti]-intA[ntr])/(pT->t[ntr]-pT->t[nti]);
          Sp = (pJet->Ip+0.5*p*p)*(pT->t[ntr]-pT->t[nti]) + 0.5*(intAsq[ntr]-intAsq[nti]) + p*(intA[ntr]-intA[nti]);
          dtme[0] = Dipole_Transition_Matrix_Element(p+A[nti]);
          dtme[1] = Dipole_Transition_Matrix_Element(p+A[ntr]);  // p is ti & tr dependent so this line must be inside nested loop
          smooth = 1.0;
          if(ntr-nti>0.8*max_ntau)
	    smooth = Tau_Smooth(ntr-nti,min_ntau,max_ntau);
          if(nti<max_ntau)
	    smooth *= Tau_Smooth(max_ntau-nti-1,min_ntau,max_ntau);
          D[ntr] += smooth*spread[ntr-nti]*conj(dtme[1])*population[nti]*dtme[0]*Et[nti]*exp(-i*Sp);
        }
        D[ntr] *= i*pT->dt*population[ntr];
        D[ntr] *= pT->AB[ntr];
      }
      break;
    case 'f':  // forward temporal integration
      for(int nt=0; nt<pT->Nt; nt++)
        D[nt] = 0.0;
      for(int nti=0; nti<pT->Nt; nti++) {
        for(int ntr=nti+min_ntau; ntr<nti+max_ntau && ntr<pT->Nt; ntr++) {
          p = (intA[nti]-intA[ntr])/(pT->t[ntr]-pT->t[nti]);
          Sp = (pJet->Ip+0.5*p*p)*(pT->t[ntr]-pT->t[nti]) + 0.5*(intAsq[ntr]-intAsq[nti]) + p*(intA[ntr]-intA[nti]);
          smooth = 1.0;
          if(ntr-nti>0.8*max_ntau)
            smooth = Tau_Smooth(ntr-nti,min_ntau,max_ntau);
          if(nti<max_ntau)
	    smooth *= Tau_Smooth(max_ntau-nti-1,min_ntau,max_ntau);
          dtme[0] = Dipole_Transition_Matrix_Element(p+A[nti]);  // p is ti & tr dependent so this line must be inside nested loop
          dtme[1] = Dipole_Transition_Matrix_Element(p+A[ntr]);  
          D[ntr] += smooth*spread[ntr-nti]*conj(dtme[1])*population[ntr]*Et[nti]*dtme[0]*exp(-i*Sp);
        } 
        D[nti] *= i*pT->dt*population[nti];
        D[nti] *= pT->AB[nti];
      }
      break;
    default:
      std::cerr << "Error in SFA method (forwards / backwards integration)\n";
#ifdef useMPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(1);
  }  
}

// Function to calculate dipole acceleration from dipole moment
void sfa_response_1d::Dipole_Moment_to_Acceleration() {
  for(int nt=0; nt<pT->Nt; nt++) {
    fftwArrayT[nt][0] = real(D[nt]);
    fftwArrayT[nt][1] = 0.0;
  }
  fftw_execute(time2freq);
  double pre_factor = pT->dt/sqrt(2.0*pi);
  for(int nt=0; nt<pT->Nt; nt++)
    DA[nt] = pre_factor*(pT->w[nt]*pT->w[nt])*(fftwArrayT[nt][0]+i*fftwArrayT[nt][1]);
}

// Dipole acceleration in freq-domain via pointer to electric_field class
void sfa_response_1d::Dipole_Acceleration() {
  Dipole_Moment(pLaser->Et);
  Dipole_Moment_to_Acceleration();
}

// Dipole acceleration in freq-domain importing row number of Ert
void sfa_response_1d::Dipole_Acceleration(int nr, bool mpi=false) {
  (mpi==true) ? Dipole_Moment(&pLaser->E_mpi[Index2V(nr,0,pT->Nt)]) : Dipole_Moment(&pLaser->Ert[Index2V(nr,0,pT->Nt)]);
  Dipole_Moment_to_Acceleration();
}

// Dipole acceleration in freq-domain
void sfa_response_1d::Dipole_Acceleration(dcmplx *Et) {
  Dipole_Moment(Et);
  Dipole_Moment_to_Acceleration();
}

// Dipole acceleration in freq-domain
void sfa_response_1d::Dipole_Acceleration(dcmplx *Ert, int nr) {
  Dipole_Moment(&Ert[Index2V(nr,0,pT->Nt)]);
  Dipole_Moment_to_Acceleration();
}

// Dipole acceleration in freq-domain from specified temporal fields
void sfa_response_1d::Dipole_Acceleration(dcmplx *EtX, dcmplx *EtY) {
  Dipole_Acceleration(EtX);  // just to give form to this virtual function in dipole_response class
}

// Dipole acceleration in freq-domain importing row number of specified temporospatial fields
void sfa_response_1d::Dipole_Acceleration(dcmplx *ErtX, dcmplx *ErtY, int nr) {
  Dipole_Acceleration(ErtX, nr);  // just to give form to this virtual function in dipole_response class
}  

/***************************************************************************************************************/
// Class to create two-dimensional dipole response via SFA method:
class sfa_response_2d : public sfa_response_1d {  
 private:    
  /***************************************************************************************************************/
  // Private Class Members:           
  double p2;
  dcmplx dtme2[2];
  double *A2, *intA2, *intAsq2;
  dcmplx *D2;
  /***************************************************************************************************************/
  // Private Class Function Prototypes:
  void Initialise();
  void Calculate_A2(dcmplx *);
  void Calculate_intA2();
  void Calculate_intAsq2();
  dcmplx Dipole_Transition_Matrix_Element(double, double);
  dcmplx Dipole_Transition_Matrix_Element(double, double, double, char);
  void Dipole_Moment(dcmplx *, dcmplx *, char);
  void Dipole_Moment_to_Acceleration();
 public:
  /***************************************************************************************************************/
  // Public Class Function Prototypes:
  sfa_response_2d(temporal_axes &T_, electric_field &Laser_, gas_jet &Jet_);
  sfa_response_2d(temporal_axes &T_, electric_field &Laser_, gas_jet &Jet_, std::ifstream &);
  ~sfa_response_2d();
  void Dipole_Acceleration(); 
  void Dipole_Acceleration(int, bool); 
  void Dipole_Acceleration(dcmplx *, dcmplx *); 
  void Dipole_Acceleration(dcmplx *, dcmplx *, int); 
};

/***************************************************************************************************************/
// Class Functions:
  
// Constructor
sfa_response_2d::sfa_response_2d(temporal_axes &T_, electric_field &Laser_, gas_jet &Jet_) : sfa_response_1d(T_, Laser_, Jet_) { 
  Initialise();
}  

// Constructor
sfa_response_2d::sfa_response_2d(temporal_axes &T_, electric_field &Laser_, gas_jet &Jet_, std::ifstream &INSTREAM) : sfa_response_1d(T_, Laser_, Jet_, INSTREAM) { 
  Initialise();
}  

// Destructor
sfa_response_2d::~sfa_response_2d() {
  delete[] A2;
  delete[] intA2;
  delete[] intAsq2;
  delete[] D2;
}

// Initialise memory
void sfa_response_2d::Initialise() {
  A2 = new double[pT->Nt];
  intA2 = new double[pT->Nt];
  intAsq2 = new double[pT->Nt];
  D2 = new dcmplx[pT->Nt];
}

// Calculate A2
void sfa_response_2d::Calculate_A2(dcmplx *Et) {
  A2[0] = 0.0;
  for(int nt=1; nt<pT->Nt-1; nt++)
    A2[nt] = A2[nt-1] - real(Et[nt])*pT->dt;
}

// Calculate intA2
void sfa_response_2d::Calculate_intA2() {
  intA2[0] = 0.0;
  for(int nt=1; nt<pT->Nt-1; nt++)
    intA2[nt] = intA2[nt-1] + A2[nt]*pT->dt;
}

// Calculate intAsq2
void sfa_response_2d::Calculate_intAsq2() {
  intAsq2[0] = 0.0;
  for(int nt=1; nt<pT->Nt-1; nt++)
    intAsq2[nt] = intAsq2[nt-1] + A2[nt]*A2[nt]*pT->dt;
}

// Dipole transition matrix element (order of arguments specifies polarisation along which axes - usually x,y)
dcmplx sfa_response_2d::Dipole_Transition_Matrix_Element(double k_par, double k_perp) {
  dcmplx element;
  double k_modsq = k_par*k_par+k_perp*k_perp;
  if(pJet->n==1)
    element = pre_dtme*k_par/pow(k_modsq+2.0*pJet->Ip, 3);
  else if(pJet->n==2)
    if(pJet->l==0)
      element = pre_dtme*k_par*(k_modsq-4.0*pJet->Ip)/pow(k_modsq+2.0*pJet->Ip, 4);
    else if(pJet->l==1)  
      element = pre_dtme*(5.0*k_modsq-2.0*pJet->Ip)/pow(k_modsq+2.0*pJet->Ip, 4);  
  else if(pJet->n==3)
    if(pJet->l==0)
      element = pre_dtme*(3.0*pow(k_modsq, 2)-18.0*k_modsq*2.0*pJet->Ip+11.0*pow(2.0*pJet->Ip, 2))/pow(k_modsq+2.0*pJet->Ip, 5);
    else if(pJet->l==1)  // this is not correct for arbitrary phi (only phi = 0) - recalculate if & when necessary 
      element = pre_dtme*(2.0*k_modsq*(3.0*k_modsq-10.0*pJet->Ip)-(pow(k_modsq, 2)-pow(2.0*pJet->Ip, 2)))/pow(k_modsq+2.0*pJet->Ip, 5);
    else
      element = pre_dtme*0.0;
  else 
    element = pre_dtme*0.0;     
  return element;
}

// Dipole transition matrix element 
dcmplx sfa_response_2d::Dipole_Transition_Matrix_Element(double kx, double ky, double phi, char axis) {
  dcmplx element;
  double k_modsq = kx*kx+ky*ky;
  bool yAxis;
  if(axis=='y')  yAxis = true;
  else  yAxis = false;  // defaults to x-axis if neither x or y entered
  if(pJet->n==1)
    element = (yAxis==true) ? pre_dtme*ky/pow(k_modsq+2.0*pJet->Ip, 3) : pre_dtme*kx/pow(k_modsq+2.0*pJet->Ip, 3);
  else if(pJet->n==2)
    if(pJet->l==0)
      element = (yAxis==true) ? pre_dtme*ky*(k_modsq-4.0*pJet->Ip)/pow(k_modsq+2.0*pJet->Ip, 4) : pre_dtme*kx*(k_modsq-4.0*pJet->Ip)/pow(k_modsq+2.0*pJet->Ip, 4);
    else if(pJet->l==1) {
      double phiP = atan2(ky,kx);
      dcmplx factorX = pre_dtme*((2.0*pJet->Ip-5.0*k_modsq)*pow(cos(phiP),2)+k_modsq*(k_modsq+2.0*pJet->Ip)*pow(sin(phiP),2))/pow(k_modsq+2.0*pJet->Ip, 4);
      dcmplx factorY = pre_dtme*((2.0*pJet->Ip-k_modsq*(5.0+k_modsq+2.0*pJet->Ip))*0.5*sin(2*(phiP)))/pow(k_modsq+2.0*pJet->Ip, 4);
      if(yAxis==true)  element = (sin(phi)*factorX+cos(phi)*factorY)*ky;
      else  element = (cos(phi)*factorX-sin(phi)*factorY)*kx;
    }
  else if(pJet->n==3)  
    if(pJet->l==0)
      element = pre_dtme*(3.0*pow(k_modsq, 2)-18.0*k_modsq*2.0*pJet->Ip+11.0*pow(2.0*pJet->Ip, 2))/pow(k_modsq+2.0*pJet->Ip, 5);
      else if(pJet->l==1)  // this is not correct for arbitrary phi (only phi = 0) - recalculate if & when necessary
      element = (yAxis==true) ? pre_dtme*(2.0*k_modsq*(3.0*k_modsq-10.0*pJet->Ip)-(pow(k_modsq, 2)-pow(2.0*pJet->Ip, 2)))*sin(phi)/pow(k_modsq+2.0*pJet->Ip, 5) : pre_dtme*(2.0*k_modsq*(3.0*k_modsq-10.0*pJet->Ip)-(pow(k_modsq, 2)-pow(2.0*pJet->Ip, 2)))*cos(phi)/pow(k_modsq+2.0*pJet->Ip, 5);
    else
      element = pre_dtme*0.0;
  else 
    element = pre_dtme*0.0;     
  // fill out for additional states if & when necessary
  return element;
}

// Dipole moment in time-domain
void sfa_response_2d::Dipole_Moment(dcmplx *EtX, dcmplx *EtY, char method='b') {
  double smooth, tau;
  dcmplx int_factor;
  double phi;
  Calculate_A(EtX);  Calculate_A2(EtY);
  Calculate_intA();  Calculate_intA2();
  Calculate_intAsq();  Calculate_intAsq2();
  pJet->ADK_Ionisation(population, EtX, EtY);
  switch(method) {
    case 'b':  // backward temporal integration
/*       for(int nt=0; nt<min_ntau; nt++) { */
/*         D[nt] = 0.0;  D2[nt] = 0.0; */
/*       } */
      for(int nt=0; nt<pT->Nt; nt++) {
        D[nt] = 0.0;  D2[nt] = 0.0;
      }
      for(int ntr=min_ntau; ntr<pT->Nt; ntr++) { 
        D[ntr] = 0.0;  D2[ntr] = 0.0;
        for(int nti=ntr-min_ntau; nti>=std::max(ntr-max_ntau,0); nti--) { 
	  // what about a "if E>Emin" statement around these calculations?
	  tau = pT->t[ntr]-pT->t[nti];
          p = (intA[nti]-intA[ntr])/tau;
          p2 = (intA2[nti]-intA2[ntr])/tau;
	  Sp = pJet->Ip*tau;
	  Sp += 0.5*p*p*tau + 0.5*(intAsq[ntr]-intAsq[nti]) + p*(intA[ntr]-intA[nti]);
	  Sp += 0.5*p2*p2*tau + 0.5*(intAsq2[ntr]-intAsq2[nti]) + p2*(intA2[ntr]-intA2[nti]);

	  phi = atan2(real(EtY[nti]), real(EtX[nti]));
          dtme[0] = Dipole_Transition_Matrix_Element(p+A[nti], p2+A2[nti], phi, 'x');
          dtme2[0] = Dipole_Transition_Matrix_Element(p+A[nti], p2+A2[nti], phi, 'y');
          dtme[1] = Dipole_Transition_Matrix_Element(p+A[ntr], p2+A2[ntr], phi, 'x');
          dtme2[1] = Dipole_Transition_Matrix_Element(p+A[ntr], p2+A2[ntr], phi, 'y');

/*           dtme[0] = Dipole_Transition_Matrix_Element(p+A[nti], p2+A2[nti]); */
/*           dtme2[0] = Dipole_Transition_Matrix_Element(p2+A2[nti], p+A[nti]); */
/*           dtme[1] = Dipole_Transition_Matrix_Element(p+A[ntr], p2+A2[ntr]);   */
/*           dtme2[1] = Dipole_Transition_Matrix_Element(p2+A2[ntr], p+A[ntr]);  */

          smooth = 1.0;
          if(ntr-nti>0.8*max_ntau)
	    smooth = Tau_Smooth(ntr-nti,min_ntau,max_ntau);
          if(nti<max_ntau)
	    smooth *= Tau_Smooth(max_ntau-nti-1,min_ntau,max_ntau);
	  int_factor = smooth*spread[ntr-nti]*population[nti]*(dtme[0]*EtX[nti]+dtme2[0]*EtY[nti])*exp(-i*Sp);
          D[ntr] += conj(dtme[1])*int_factor;
	  D2[ntr] += conj(dtme2[1])*int_factor;
	}
        D[ntr] *= i*pT->dt*population[ntr];  D2[ntr] *= i*pT->dt*population[ntr];
	D[ntr] *= pT->AB[ntr];  D2[ntr] *= pT->AB[ntr];
      }
      break;
    case 'f':  // forward temporal integration
      for(int nt=0; nt<pT->Nt; nt++) {
        D[nt] = 0.0;  D2[nt] = 0.0;
      }
      for(int nti=0; nti<pT->Nt; nti++) {
        for(int ntr=nti+min_ntau; ntr<nti+max_ntau && ntr<pT->Nt; ntr++) {
	  tau = pT->t[ntr]-pT->t[nti];
          p = (intA[nti]-intA[ntr])/tau;
	  p2 = (intA2[nti]-intA2[ntr])/tau;
	  Sp = pJet->Ip*tau;
	  Sp += 0.5*p*p*tau + 0.5*(intAsq[ntr]-intAsq[nti]) + p*(intA[ntr]-intA[nti]);
	  Sp += 0.5*p2*p2*tau + 0.5*(intAsq2[ntr]-intAsq2[nti]) + p2*(intA2[ntr]-intA2[nti]);

	  phi = atan2(real(EtY[nti]), real(EtX[nti]));
          dtme[0] = Dipole_Transition_Matrix_Element(p+A[nti], p2+A2[nti], phi, 'x');
          dtme2[0] = Dipole_Transition_Matrix_Element(p+A[nti], p2+A2[nti], phi, 'y');
          dtme[1] = Dipole_Transition_Matrix_Element(p+A[ntr], p2+A2[ntr], phi, 'x');  
          dtme2[1] = Dipole_Transition_Matrix_Element(p+A[ntr], p2+A2[ntr], phi, 'y'); 

          smooth = 1.0;
          if(ntr-nti>0.8*max_ntau)
            smooth = Tau_Smooth(ntr-nti,min_ntau,max_ntau);
          if(nti<max_ntau)
	    smooth *= Tau_Smooth(max_ntau-nti-1,min_ntau,max_ntau);
   	  int_factor = smooth*spread[ntr-nti]*population[ntr]*(dtme[0]*EtX[nti]+dtme2[0]*EtY[nti])*exp(-i*Sp);
          D[ntr] += conj(dtme[1])*int_factor;
	  D2[ntr] += conj(dtme2[1])*int_factor;
        } 
        D[nti] *= i*pT->dt*population[nti];  D2[nti] *= i*pT->dt*population[nti];
        D[nti] *= pT->AB[nti];  D2[nti] *= pT->AB[nti];
      }
      break;
    default:
      std::cerr << "Error in SFA method (forwards / backwards integration)\n";
#ifdef useMPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(1);
  }  
}

// Function to calculate dipole acceleration from dipole moment
void sfa_response_2d::Dipole_Moment_to_Acceleration() {
  sfa_response_1d::Dipole_Moment_to_Acceleration();
  for(int nt=0; nt<pT->Nt; nt++) {
    fftwArrayT[nt][0] = real(D2[nt]);
    fftwArrayT[nt][1] = 0.0;
  }
  fftw_execute(time2freq);
  double pre_factor = pT->dt/sqrt(2.0*pi);
  for(int nt=0; nt<pT->Nt; nt++)
    DA2[nt] = pre_factor*(pT->w[nt]*pT->w[nt])*(fftwArrayT[nt][0]+i*fftwArrayT[nt][1]);
}

// Dipole acceleration in freq-domain via pointer to electric_field class
void sfa_response_2d::Dipole_Acceleration() {
  Dipole_Moment(pLaser->Et, pLaser->Et2);
  Dipole_Moment_to_Acceleration();
}

// Dipole acceleration in freq-domain importing row number of Ert
void sfa_response_2d::Dipole_Acceleration(int nr, bool mpi=false) {
  (mpi==true) ? Dipole_Moment(&pLaser->E_mpi[Index2V(nr,0,pT->Nt)], &pLaser->E_mpi2[Index2V(nr,0,pT->Nt)]) : Dipole_Moment(&pLaser->Ert[Index2V(nr,0,pT->Nt)], &pLaser->Ert2[Index2V(nr,0,pT->Nt)]);
  Dipole_Moment_to_Acceleration();
}

// Dipole acceleration in freq-domain from specified temporal fields
void sfa_response_2d::Dipole_Acceleration(dcmplx *EtX, dcmplx *EtY) {
  Dipole_Moment(EtX, EtY);
  Dipole_Moment_to_Acceleration();
}

// Dipole acceleration in freq-domain importing row number of specified temporospatial fields
void sfa_response_2d::Dipole_Acceleration(dcmplx *ErtX, dcmplx *ErtY, int nr) {
  Dipole_Moment(&ErtX[Index2V(nr,0,pT->Nt)], &ErtY[Index2V(nr,0,pT->Nt)]);
  Dipole_Moment_to_Acceleration();
}  

/***************************************************************************************************************/
#endif


