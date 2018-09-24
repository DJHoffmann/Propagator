/***************************************************************************************************************/
// Header file to contain qo_response_1/2d classes for use in Propagator.cpp 

/***************************************************************************************************************/
// Notes:

// qo_response_2d assumes exactly 2 colours

// In principle maybe no reason to ignore high free-electron density: 
// (i) drop off in amplitude could be added through t-dependent decay and some chain-rule like solution
// (ii) red/blue-shifting through t-dependence in frequency???

//  need to change Characterise_Orbits() for the 2D case

/***************************************************************************************************************/
#ifndef QO_RESPONSE_H
#define QO_RESPONSE_H

#include<iostream>
#include<vector>
#include<Constants_Conversions.h>
#include<Laser_Parameters.h>
#include<Electric_Field.h>

#define __airy_functions_complex_double_airy_aic __airy_functions_complex_double_MOD_airy_aic
extern "C" void __airy_functions_complex_double_airy_aic(dcmplx*, dcmplx*, dcmplx*, int*, int*, int*);

/***************************************************************************************************************/
// New Definitions:
typedef std::vector<int> intvec;
typedef std::vector<double> dblvec;
typedef std::vector<dcmplx> dcmplxvec;

/***************************************************************************************************************/
// Class to create one-dimensional dipole response via QO method:
class qo_response_1d : public dipole_response {  
 protected:    
  /***************************************************************************************************************/
  // Protected Class Members:
  bool analyse, findAll;
  char method;  // S=saddle-point approx, U=uniform approx
  char traj;  // S=short, L=long, B=both(short&long), C=combined(short+long), A=all
  int colourRef;
  int min_ntp, max_ntp;
  int min_ntau, max_ntau;
  double min_tp, max_tp;
  double min_tau, max_tau;
  double t_shift;
  double inc_frac;
  laser_parameters **pLP;
  int Nguess;  // max number of guesses
  int Nit;     // number of iterations
  int NEn;     // number of energy points
  double dEn;
  dblvec En;
  /* double intA2s[pLaser->colours], intA2f[pLaser->colours], env_range[pLaser->colours]; */
  double *intA2s, *intA2f, *env_range;
  dcmplxvec t_saddle, tp_saddle;
  dcmplxvec t_onepath, tp_onepath;
  std::vector<dcmplxvec> t_orbits, tp_orbits;
  std::vector<intvec> orbit_labels;
  dblvec population;
  dcmplxvec spectrum_short, spectrum_long, spectrum_pair;
  dcmplxvec DA_local; 
  /***************************************************************************************************************/
  // Protected Class Function Prototypes:
  void Initialise();
  void Specify_Orbits();
  void Specify_Orbits(std::ifstream &);
  void Set_Bounds();
  void Build_Energy_Axis();
  void Calculate_FieldFactors();
  dcmplx E(dcmplx);
  dcmplx A(dcmplx);
  dcmplx intA(dcmplx);
  dcmplx intA2(dcmplx);  
  dcmplx p(dcmplx, dcmplx);
  dcmplx S(dcmplx, dcmplx, double);
  dcmplx dtdS(dcmplx, dcmplx, double);
  dcmplx dtpdS(dcmplx, dcmplx);
  dcmplx dttd2S(dcmplx, dcmplx);
  dcmplx dtptpd2S(dcmplx, dcmplx);
  dcmplx dtptd2S(dcmplx, dcmplx);
  dcmplx Determinant_Hessian(dcmplx, dcmplx);
  bool Find_Initial_Saddle(dcmplx &, dcmplx &);
  bool Find_Saddle(dcmplx &, dcmplx &, double, bool);
  double Calculate_Increments(dcmplx &, dcmplx &, dcmplx, dcmplx, double);
  bool Check_Convergence(dcmplx, dcmplx, double, bool &, bool);
  bool Add_Saddle(dcmplx, dcmplx);
  void Sort_Saddle_Time();
  void Saddlepoint_Search();
  bool Track_Orbits();
  void Characterise_Orbits();
  void Calculate_Population();
  void Fit_Population(double &, double &, dcmplx, dcmplx, double, double);
  dcmplx DTME_Prefactor();
  dcmplx DTME(dcmplx, dcmplx);
  dcmplx DTME_singular(dcmplx, dcmplx); 
  void Dipole_Acceleration_SaddlePoint();
  void Dipole_Acceleration_Uniform();
  void Dipole_Acceleration_Combined();
  void Fit_to_wAxis();
 public:
  /***************************************************************************************************************/
  // Public Class Function Prototypes:
  //qo_response_1d(laser_parameters &);
  qo_response_1d(temporal_axes &, electric_field &, gas_jet &, bool);
  qo_response_1d(temporal_axes &, electric_field &, gas_jet &, std::ifstream &);
  ~qo_response_1d();
  void Dipole_Acceleration(dcmplx *);
  void Dipole_Acceleration(); 
  void Dipole_Acceleration(int, bool);
  void Dipole_Acceleration(dcmplx *, int); 
  void Dipole_Acceleration(dcmplx *, dcmplx *); 
  void Dipole_Acceleration(dcmplx *, dcmplx *, int); 
};

/***************************************************************************************************************/
// Class Functions:
 
/* // Constructor using imported laser parameter specification - DO LATER */
/* qo_response_1d::qo_response_1d(laser_parameters &LP_) : dipole_response(.................) {  */
/*   analyse = false; */
/*   findAll = false; */
/*   Initialise(); */
/*   //DTME_Prefactor(); */
/* }  */ 

// Constructor diagnosing laser parameters from imported electric field (creates new laser_parameters object then copies this to pLP)
qo_response_1d::qo_response_1d(temporal_axes &T_, electric_field &Laser_, gas_jet &Jet_, bool findAll_=false) : dipole_response(T_, Laser_, Jet_, 'Q') {
  analyse = true;  
  findAll = findAll_;
  method = 'U';
  traj = 'B';
  colourRef = 0; 
  Nguess = 10000; //1000; 
  Nit = 20; //100;
  NEn = 1024;
  inc_frac = 0.1;//0.25;
  Initialise();
  Specify_Orbits();
  // DTME_Prefactor();  // put this within Dipole_Acceleration() - could feasibly change if >1 jet
}  

// As above but takes in control parameters from specified ifstream
qo_response_1d::qo_response_1d(temporal_axes &T_, electric_field &Laser_, gas_jet &Jet_, std::ifstream &INSTREAM) : dipole_response(T_, Laser_, Jet_, 'Q') { 
  analyse = true;  
  findAll = (Read_Variable<char>(INSTREAM, "QOtrackAll")=='y') ? true : false;
  method = Read_Variable<char>(INSTREAM, "QOmethod");
  traj = Read_Variable<char>(INSTREAM, "QOtraj");
  colourRef = Read_Variable<int>(INSTREAM, "QOcolour");  
  Nguess = 1000;  //10000; //1000; 
  Nit = 100;  //20; //100;
  NEn = 1024;
  inc_frac = 0.3;//0.1;//0.25;
  Initialise();
  Specify_Orbits(INSTREAM);
  // DTME_Prefactor();  // put this within Dipole_Acceleration() - could feasibly change if >1 jet
}  

// Destructor
qo_response_1d::~qo_response_1d() {
  for(int nc=0; nc<pLaser->colours; nc++)  delete[] pLP[nc];
  delete[] pLP;
  delete[] intA2s, intA2f, env_range;
  DA_local.clear();
  if(method=='U'||method=='C')  spectrum_pair.clear();
  if(method=='C') {spectrum_short.clear();  spectrum_long.clear();}  
}

// Initialise memory
void qo_response_1d::Initialise() {
  pLP = new laser_parameters*[pLaser->colours];
  t_onepath.resize(NEn);  tp_onepath.resize(NEn);  
  En.resize(NEn);  DA_local.resize(NEn);
  intA2s = new double[pLaser->colours]; 
  intA2f = new double [pLaser->colours]; 
  env_range = new double[pLaser->colours];
/*   if(method=='C') {spectrum_short.resize(NEn);  spectrum_long.resize(NEn);  spectrum_pair.resize(NEn);} */
  if(method=='U'||method=='C')  spectrum_pair.resize(NEn);
  if(method=='C') {spectrum_short.resize(NEn);  spectrum_long.resize(NEn);}
}

// Sets which orbits will be calculated subsequently - default choices
void qo_response_1d::Specify_Orbits() {
  min_ntp = -2;  max_ntp = 2;
  min_ntau = 0;  max_ntau = 1;
}

// Sets which orbits will be calculated subsequently based on filestream taken in by argument
void qo_response_1d::Specify_Orbits(std::ifstream &INSTREAM) {
  min_ntp = Read_Variable<int>(INSTREAM, "QOminIon");  max_ntp = Read_Variable<int>(INSTREAM, "QOmaxIon");
  min_ntau = Read_Variable<int>(INSTREAM, "QOminTau");  max_ntau = Read_Variable<int>(INSTREAM, "QOmaxTau");
  int temp;
  if(min_ntp>max_ntp) {
    std::cout << "Ordering error in QO ionisation time selection ... switching\n";
    temp = min_ntp;  min_ntp = max_ntp;  max_ntp = temp;
  }
  if(min_ntau>max_ntau) {
    std::cout << "Time-ordering error in QO excursion time selection ... switching\n";
    temp = min_ntau;  min_ntau = max_ntau;  max_ntau = temp;
  }
}

// Sets bounds on ionisation time and excursion duration based on half-cycle integer labels
void qo_response_1d::Set_Bounds() {
  min_tau = (min_ntau==0 && pLaser->colours>1) ? 0.0 : 0.45*pLP[colourRef]->T+min_ntau*(0.5*pLP[colourRef]->T);
  max_tau = 0.45*pLP[colourRef]->T+max_ntau*(0.5*pLP[colourRef]->T); 
  t_shift = -fmod(pLP[colourRef]->cep,pi)*pLP[colourRef]->T/(2.0*pi);
  if(t_shift<-0.25)  t_shift+=0.5;
  else if(t_shift>=0.25)  t_shift-=0.5;
  min_tp = pLP[colourRef]->t0+t_shift+pLP[colourRef]->T*(-0.25+0.5*min_ntp); 
  max_tp = pLP[colourRef]->t0+t_shift+pLP[colourRef]->T*(0.25+0.5*max_ntp);
}

// Function to build the energy axis appropriate to these field values (& Ip)
void qo_response_1d::Build_Energy_Axis() {
  double En_lower = pJet->Ip+0.1*eV_to_au;//pJet->Ip; 
  double En_upper;
  //if(pLaser->colours==1)  En_upper = 1.5*(1.32*pJet->Ip+3.17*pLP[0]->Ponderomotive_Potential());  // 1.5 * classical cut-off
  //if(pLaser->colours==1)  En_upper = std::max(3.0*pJet->Ip,1.5*(1.32*pJet->Ip+3.17*pLP[0]->Ponderomotive_Potential()));
  if(pLaser->colours==1)  En_upper = std::max(3.0*pJet->Ip,3.0*(1.32*pJet->Ip+3.17*pLP[0]->Ponderomotive_Potential()));
  else {
    double E0_mix = 0.0, w0_mix = 0.0;
    for(int nc=0; nc<pLaser->colours; nc++) {
      E0_mix += pLP[nc]->E0;
      w0_mix += pLP[nc]->w0*pLP[nc]->E0;
    }
    w0_mix /= E0_mix;
    //En_upper = 1.5*(1.32*pJet->Ip+3.17*pow(0.5*E0_mix/w0_mix,2));
    En_upper = 3.0*(1.32*pJet->Ip+3.17*pow(0.5*E0_mix/w0_mix,2));
  }
  //En[0] = En_upper;
  En[0] = std::max(3.0*pJet->Ip,En_upper);
  dEn = (En_upper-En_lower)/(double)(NEn-1);
  for(int nEn=1; nEn<NEn; nEn++)  En[nEn] = En[nEn-1]-dEn;
}

// Function to precalculate initial/final values of intA2 and envelope boundary
void qo_response_1d::Calculate_FieldFactors() {
  for(int nc=0; nc<pLaser->colours; nc++) {
    env_range[nc] = 0.5*pi/pLP[nc]->B;
    intA2s[nc] = 0.5*pow(0.5*pLP[nc]->E0/pLP[nc]->w0,2) * (1.5*(-env_range[nc])+sin(2.0*pLP[nc]->B*(-env_range[nc]))/pLP[nc]->B+0.125*sin(4.0*pLP[nc]->B*(-env_range[nc]))/pLP[nc]->B-0.25*(3.0/pLP[nc]->w0+4.0*pLP[nc]->w0/(pLP[nc]->w0*pLP[nc]->w0-pLP[nc]->B*pLP[nc]->B)*cos(2.0*pLP[nc]->B*(-env_range[nc]))+pLP[nc]->w0/(pLP[nc]->w0*pLP[nc]->w0-4.0*pLP[nc]->B*pLP[nc]->B)*cos(4.0*pLP[nc]->B*(-env_range[nc])))*sin(2.0*(pLP[nc]->w0*(-env_range[nc])+pLP[nc]->cep))+0.5*pLP[nc]->B*(2.0/(pLP[nc]->w0*pLP[nc]->w0-pLP[nc]->B*pLP[nc]->B)*sin(2.0*pLP[nc]->B*(-env_range[nc]))+1.0/(pLP[nc]->w0*pLP[nc]->w0-4.0*pLP[nc]->B*pLP[nc]->B)*sin(4.0*pLP[nc]->B*(-env_range[nc])))*cos(2.0*(pLP[nc]->w0*(-env_range[nc])+pLP[nc]->cep)));
    intA2f[nc] = 0.5*pow(0.5*pLP[nc]->E0/pLP[nc]->w0, 2) * (1.5*env_range[nc]+sin(2.0*pLP[nc]->B*env_range[nc])/pLP[nc]->B+0.125*sin(4.0*pLP[nc]->B*env_range[nc])/pLP[nc]->B-0.25*(3.0/pLP[nc]->w0+4.0*pLP[nc]->w0/(pLP[nc]->w0*pLP[nc]->w0-pLP[nc]->B*pLP[nc]->B)*cos(2.0*pLP[nc]->B*env_range[nc])+pLP[nc]->w0/(pLP[nc]->w0*pLP[nc]->w0-4.0*pLP[nc]->B*pLP[nc]->B)*cos(4.0*pLP[nc]->B*env_range[nc]))*sin(2.0*(pLP[nc]->w0*env_range[nc]+pLP[nc]->cep))+0.5*pLP[nc]->B*(2.0/(pLP[nc]->w0*pLP[nc]->w0-pLP[nc]->B*pLP[nc]->B)*sin(2.0*pLP[nc]->B*env_range[nc])+1.0/(pLP[nc]->w0*pLP[nc]->w0-4.0*pLP[nc]->B*pLP[nc]->B)*sin(4.0*pLP[nc]->B*env_range[nc]))*cos(2.0*(pLP[nc]->w0*env_range[nc]+pLP[nc]->cep))) - intA2s[nc];
  }
}

// Electric field
dcmplx qo_response_1d::E(dcmplx t) {
  dcmplx Et = 0.0;
  for(int nc=0; nc<pLaser->colours; nc++)
    if(real(t-pLP[nc]->t0)>=-env_range[nc] && real(t-pLP[nc]->t0)<=env_range[nc]) 
      Et += pLP[nc]->E0/pLP[nc]->w0 * cos((t-pLP[nc]->t0)*pLP[nc]->B) * (pLP[nc]->w0*cos((t-pLP[nc]->t0)*pLP[nc]->B)*cos(pLP[nc]->w0*(t-pLP[nc]->t0)+pLP[nc]->cep)-2.0*pLP[nc]->B*sin((t-pLP[nc]->t0)*pLP[nc]->B)*sin(pLP[nc]->w0*(t-pLP[nc]->t0)+pLP[nc]->cep));
  return Et;
}

// Vector Potential
dcmplx qo_response_1d::A(dcmplx t) {
  dcmplx At = 0.0;
  for(int nc=0; nc<pLaser->colours; nc++)
    if(real(t-pLP[nc]->t0)>=-env_range[nc] && real(t-pLP[nc]->t0)<=env_range[nc]) 
      At += - pLP[nc]->E0/pLP[nc]->w0 * pow(cos((t-pLP[nc]->t0)*pLP[nc]->B), 2) * sin(pLP[nc]->w0*(t-pLP[nc]->t0)+pLP[nc]->cep);
  return At;
}

// Integral of Vector Potential
dcmplx qo_response_1d::intA(dcmplx t) {
  dcmplx intAt = 0.0;
  for(int nc=0; nc<pLaser->colours; nc++)
    if(real(t-pLP[nc]->t0)>=-env_range[nc] && real(t-pLP[nc]->t0)<=env_range[nc]) 
      intAt += 0.5*pLP[nc]->E0/pLP[nc]->w0 * (cos(pLP[nc]->w0*(t-pLP[nc]->t0)+pLP[nc]->cep)/pLP[nc]->w0+1.0/(pLP[nc]->w0*pLP[nc]->w0-4.0*pLP[nc]->B*pLP[nc]->B)*(pLP[nc]->w0*cos(2.0*pLP[nc]->B*(t-pLP[nc]->t0))*cos(pLP[nc]->w0*(t-pLP[nc]->t0)+pLP[nc]->cep)+2.0*pLP[nc]->B*sin(2.0*pLP[nc]->B*(t-pLP[nc]->t0))*sin(pLP[nc]->w0*(t-pLP[nc]->t0)+pLP[nc]->cep)));
  return intAt;
}

// Integral of (Vector Potential Squared)
dcmplx qo_response_1d::intA2(dcmplx t) {
  dcmplx intAt2 = 0.0;
  for(int nc=0; nc<pLaser->colours; nc++)
    if(real(t-pLP[nc]->t0)>=-env_range[nc] && real(t-pLP[nc]->t0)<=env_range[nc]) 
      intAt2 += 0.5*pow(0.5*pLP[nc]->E0/pLP[nc]->w0, 2) * (1.5*(t-pLP[nc]->t0)+sin(2.0*pLP[nc]->B*(t-pLP[nc]->t0))/pLP[nc]->B+0.125*sin(4.0*pLP[nc]->B*(t-pLP[nc]->t0))/pLP[nc]->B-0.25*(3.0/pLP[nc]->w0+4.0*pLP[nc]->w0/(pLP[nc]->w0*pLP[nc]->w0-pLP[nc]->B*pLP[nc]->B)*cos(2.0*pLP[nc]->B*(t-pLP[nc]->t0))+pLP[nc]->w0/(pLP[nc]->w0*pLP[nc]->w0-4.0*pLP[nc]->B*pLP[nc]->B)*cos(4.0*pLP[nc]->B*(t-pLP[nc]->t0)))*sin(2.0*(pLP[nc]->w0*(t-pLP[nc]->t0)+pLP[nc]->cep))+0.5*pLP[nc]->B*(2.0/(pLP[nc]->w0*pLP[nc]->w0-pLP[nc]->B*pLP[nc]->B)*sin(2.0*pLP[nc]->B*(t-pLP[nc]->t0))+1.0/(pLP[nc]->w0*pLP[nc]->w0-4.0*pLP[nc]->B*pLP[nc]->B)*sin(4.0*pLP[nc]->B*(t-pLP[nc]->t0)))*cos(2.0*(pLP[nc]->w0*(t-pLP[nc]->t0)+pLP[nc]->cep))) - intA2s[nc];
    else if(real(t-pLP[nc]->t0)>env_range[nc])  intAt2 += intA2f[nc];
  return intAt2;
}

// Complex Saddle-point Momentum
dcmplx qo_response_1d::p(dcmplx t, dcmplx tp) {
  return (intA(tp)-intA(t))/(t-tp);
}

// Action
dcmplx qo_response_1d::S(dcmplx t, dcmplx tp, double Ef) {
  dcmplx ps = p(t,tp); 
  return -(pJet->Ip*(t-tp)-0.5*ps*ps*(t-tp)+0.5*(intA2(t)-intA2(tp))-Ef*t);
}

// dS/dt
dcmplx qo_response_1d::dtdS(dcmplx t, dcmplx tp, double Ef) {
  dcmplx ps = p(t,tp); 
  return -(pJet->Ip-Ef+0.5*pow(ps+A(t), 2.0));
}

// dS/dtp
dcmplx qo_response_1d::dtpdS(dcmplx t, dcmplx tp) {
  dcmplx ps = p(t,tp); 
  return -(-pJet->Ip-0.5*pow(ps+A(tp), 2.0));
}

// d2S/dtdt
dcmplx qo_response_1d::dttd2S(dcmplx t, dcmplx tp) {
  dcmplx ps = p(t,tp); 
  return -(-(ps+A(t))*(E(t)+(ps+A(t))/(t-tp)));
}

// d2S/dtpdtp
dcmplx qo_response_1d::dtptpd2S(dcmplx t, dcmplx tp) {
  dcmplx ps = p(t,tp); 
  return -((ps+A(tp))*(E(tp)-(ps+A(tp))/(t-tp)));
}

// d2S/dtpdt
dcmplx qo_response_1d::dtptd2S(dcmplx t, dcmplx tp) {
  dcmplx ps = p(t,tp); 
  return -((ps+A(t))*(ps+A(tp))/(t-tp));
}

// Determinant of Hessian Matrix
dcmplx qo_response_1d::Determinant_Hessian(dcmplx t, dcmplx tp) {
  dcmplx DtptD2S = dtptd2S(t,tp);
  return dttd2S(t,tp)*dtptpd2S(t,tp)-DtptD2S*DtptD2S;
}

// Find saddle points at En[0]
bool qo_response_1d::Find_Initial_Saddle(dcmplx &t, dcmplx &tp) {
  dcmplx t_guess(0.0), tp_guess(0.0), tau_guess(0.0);
  bool found = false;
  real(tp_guess) = fmod(rand(), max_tp-min_tp) + min_tp;
  real(tau_guess) = fmod(rand(), max_tau-min_tau) + min_tau;
  real(t_guess) = real(tp_guess) + real(tau_guess);
  imag(tp_guess) = fmod(rand(), 0.7*pLP[colourRef]->T);
  imag(t_guess) = fmod(rand(), 1.4*pLP[colourRef]->T) - 0.7*pLP[colourRef]->T;  
  found = Find_Saddle(t_guess, tp_guess, En[0], true);
  if(found==true) {t = t_guess;  tp = tp_guess;}
  return found;
}

// Find single saddle point at a specifed energy
bool qo_response_1d::Find_Saddle(dcmplx &t, dcmplx &tp, double Energy, bool boundcheck=true) {
  dcmplx dt(0.0), dtp(0.0);
  double error = 1.0e10;//0.0;
  bool found = false, stop = false;
  for(int nit=0; nit<Nit; nit++) {
    error = Calculate_Increments(dt, dtp, t, tp, Energy);
    found = Check_Convergence(t, tp, error, stop, boundcheck);
    if(found==true)  return true;
    else if(stop==true)  return false;
    // t += dt;  tp += dtp;
    if(nit!=Nit-1) {t += dt;  tp += dtp;}
  }
  //cout << stop << '\t' << found << '\t' << error << '\n';
  return found;
}

// Function to calculate increments to t,t' and return a measure of convergene in 1st derivatives of action
double qo_response_1d::Calculate_Increments(dcmplx &dt, dcmplx &dtp, dcmplx t, dcmplx tp, double Energy) {
  dcmplx DtDS = dtdS(t,tp,Energy);
  dcmplx DtpDS = dtpdS(t,tp);
  dcmplx DttD2S = dttd2S(t,tp);
  dcmplx DtptpD2S = dtptpd2S(t,tp);
  dcmplx DtptD2S = dtptd2S(t,tp);
  dcmplx DetHess = Determinant_Hessian(t,tp);
  dt = (DtpDS*DtptD2S-DtDS*DtptpD2S)/DetHess;
  dtp = (DtDS*DtptD2S-DtpDS*DttD2S)/DetHess;
  double error = abs(DtDS) + abs(DtpDS);
  //double error = abs(dt) + abs(dtp);
  return error;
}

// Function to check convergence to saddle-point solution
bool qo_response_1d::Check_Convergence(dcmplx t, dcmplx tp, double error, bool &stop, bool boundcheck) {
  dcmplx tau = t-tp;
  if(error<1.0e-8) {
    if(boundcheck==true) {
      if((real(tau)>min_tau) && (real(tau)<max_tau)) { 
	if((real(tp)>=min_tp) && (real(tp)<=max_tp)) {stop = false;  return true;} 
	else {stop = true;  return false;}
      } 
      else {stop = true;  return false;}
    }
    else {stop = false;  return true;}
  }
  else {stop = false;  return false;}
}

// Add found saddle to vector if not already present and satifying  Im(t')>=0
bool qo_response_1d::Add_Saddle(dcmplx t, dcmplx tp) {
  bool unique = false;
  if(imag(tp)>=0.0) {
    unique = true;
    for(int ns=0; ns<t_saddle.size(); ns++) 
      //if((abs(t_saddle[ns]-t)<1.0e-5) && (abs(tp_saddle[ns]-tp)<1.0e-5))  unique = false;
      if((abs(t_saddle[ns]-t)<1.0e-5*pLP[colourRef]->T) && (abs(tp_saddle[ns]-tp)<1.0e-5*pLP[colourRef]->T))  unique = false;
    if(unique==true) {
      t_saddle.push_back(t);
      tp_saddle.push_back(tp);
    }
  }
  return unique;
}

// Sort final element into the vector according to ascending t
void qo_response_1d::Sort_Saddle_Time() {
  int ns = 0;
  dcmplx temp;
  while(ns<t_saddle.size() && (real(t_saddle.back())>real(t_saddle[ns])))  ns++;
  if(ns<t_saddle.size()-1) {
    temp = t_saddle.back();
    t_saddle.insert(t_saddle.begin()+ns, temp);
    t_saddle.pop_back();
    temp = tp_saddle.back();
    tp_saddle.insert(tp_saddle.begin()+ns, temp);
    tp_saddle.pop_back();
  }
}

// Find all saddle points within specified region (and up to Nguess)
void qo_response_1d::Saddlepoint_Search() {
  t_saddle.clear();  tp_saddle.clear();
  dcmplx t(0,0), tp(0,0);
  bool found = false;
  for(int nguess=0; nguess<Nguess; nguess++) {
    found = Find_Initial_Saddle(t, tp);
    if(found==true) {
      found = Add_Saddle(t, tp);
      if(found==true) {
	Sort_Saddle_Time();
	//std::cout << "Sorted saddle tp=(" << real(tp)/pLP[colourRef]->T << ',' << imag(tp)/pLP[colourRef]->T << "); t=(" << real(t)/pLP[colourRef]->T << ',' << imag(t)/pLP[colourRef]->T << "), spread = " << arg(1.0e-10+i*(t-tp))/pi << " pi rad\n";
      }
      found = false; 
    }
  }
  //std::cout << "Found " << t_saddle.size() << " unique saddle-points and sorted according to ascending re(t)\n";
}

// Function to follow saddle-point orbits with decreasing energy
bool qo_response_1d::Track_Orbits() {
  t_orbits.clear();  tp_orbits.clear();
  dcmplx t, tp;
  dcmplx t_prev, tp_prev;
  dcmplx dt, dtp;
  bool found;
  int nEn = 1;
  double Energy = En[0];
  dcmplx S_track[3];
  double arg_S = 0.0;
  const int Nstep = 10, Nlevel = 10;//5;
  const double step_factor = 1.0/(double)(Nstep);
  int step_level, step_count[Nlevel], old_level;
  t_orbits.clear();  tp_orbits.clear();
  t_orbits.resize(t_saddle.size(), t_onepath);  tp_orbits.resize(tp_saddle.size(), tp_onepath);
  for(int ns=0; ns<t_saddle.size(); ns++) {
    nEn = 1;  Energy = En[0];
    t_prev = (t_onepath[0] = t_saddle[ns]);
    tp_prev = (tp_onepath[0] = tp_saddle[ns]);
    dt = 0.0;  dtp = 0.0;
    S_track[0] = 0.0;  S_track[1] = 0.0;  S_track[2] = 0.0;
    arg_S = 0.0;
    step_level = 0;
    for(int nlevel=0; nlevel<Nlevel; nlevel++)  step_count[nlevel] = 0;
    while(nEn<NEn) {
/*       Energy -= dEn; // this seems to work better....? */
/*       t = t_prev+dt; */
/*       tp = tp_prev+dtp; */
      
      Energy -= dEn*pow(step_factor, step_level);
      t = t_prev+dt*pow(step_factor, step_level);
      tp = tp_prev+dtp*pow(step_factor, step_level);

      found = Find_Saddle(t, tp, Energy, false);
      if(found==true) {
	S_track[2] = S(t, tp, Energy);
	if(nEn==1)  S_track[0] = S_track[2];
	else if(nEn==2)  S_track[1] = S_track[2];
	else {
/* 	  arg_S = abs(arg(S_track[2]-S_track[1])-arg(S_track[1]-S_track[0])); */
/* 	  S_track[0] = S_track[1];  S_track[1] = S_track[2];  */
/* 	  if(arg_S>pi)  arg_S = 2.0*pi-arg_S; */

	  if(fabs(real(S_track[2])-real(S_track[1]))>3.0*fabs(real(S_track[1])-real(S_track[0])))  found = false;
	  else {
	    arg_S = abs(arg(S_track[2]-S_track[1])-arg(S_track[1]-S_track[0]));
	    S_track[0] = S_track[1];  S_track[1] = S_track[2]; 
	    if(arg_S>pi)  arg_S = 2.0*pi-arg_S;
	  }
	}
      }

      // if(found==true && arg_S<0.5*pi) { 
      if(found==true) {
	if(step_level>0)  step_count[step_level]++;
	old_level = step_level;
	for(int nlevel=step_level; nlevel>0; nlevel--) {
	  if(step_count[nlevel]==Nstep) {
	    step_level--;  // if reached a more major point, raise back 1 level
	    //std::cout << "Rejoined level " << step_level << " after " << step_count[step_level+1] << " steps.\n";
	    step_count[step_level]++;  // and advance that level one, recheck
	  }
	  else  break;
	}	
	if(step_level==0) {
	  t_onepath[nEn] = t;
	  tp_onepath[nEn] = tp;
	  nEn++;
	}
	dt = (t-t_prev)*pow(step_factor, step_level-old_level);
	dtp = (tp-tp_prev)*pow(step_factor, step_level-old_level);
	t_prev = t;  tp_prev = tp;
      }
      else {
	//std::cout << "Reverting back from " << En*27.212;
	Energy += dEn*pow(step_factor, step_level);  // revert back
	//std::cout << "eV to " << En*27.212 << "eV\n";
	step_level++;  step_count[step_level] = 0;  // drop down level & reset step_count at new level

	if(step_level==Nlevel) {std::cout << "Nlevel reached (ns = " << ns << ", nEn = " << nEn << ") - total failure!\n";  return false;}
	//if(step_level==Nlevel) {std::cout << "Nlevel reached (ns = " << ns << ", nEn = " << nEn << ") - total failure!\n";  continue;} 

	dt *= step_factor;  dtp *= step_factor;  // reduce dt & dtp according to dropped level
	//std::cout << "Track_Orbit failed for ns = " << ns << " on nEn = " << nEn << ',' << step_level-1 << ',' << step_count[step_level-1] << " (En = " << En*27.212 << "eV vs En_cutoff = " << En_cutoff*27.212 << "eV) ... descending to level " << step_level << " of " << Nlevel << ".\n";
      }
    }
    t_orbits[ns] = t_onepath;
    tp_orbits[ns] = tp_onepath;
  }
  return true;
}

// Function to ascribe 3 labels for each orbit (t position wrt t0[=0] vs T, tau=#cycles vs T, short[=-1]/long[=1])
void qo_response_1d::Characterise_Orbits() {
  const int Nl = 7;
  dcmplx tau;
  intvec one_label;
  double T_1 = pLP[colourRef]->T;
  double t0_1 = pLP[colourRef]->t0;
  int colourRefAlt = (colourRef==0) ? 1 : 0;
  double T_2 = (pLaser->colours==1 || pLP[colourRefAlt]->E0==0.0) ? T_1 : pLP[colourRefAlt]->T;
  double t0_2 = (pLaser->colours==1 || pLP[colourRefAlt]->E0==0.0) ? t0_1 : pLP[colourRefAlt]->t0;
  double t_shift_2 = -fmod(pLP[colourRefAlt]->cep,pi)*T_2/(2.0*pi);
  if(t_shift_2<-0.25)  t_shift_2+=0.5;
  else if(t_shift_2>=0.25)  t_shift_2-=0.5;
  int min_ntp_2 = floor(pLP[colourRefAlt]->w0/pLP[colourRef]->w0*min_ntp)-1;
  int max_ntp_2 = floor(pLP[colourRefAlt]->w0/pLP[colourRef]->w0*max_ntp)+1;
  int min_ntau_2 = floor(pLP[colourRefAlt]->w0/pLP[colourRef]->w0*min_ntau)-1;
  int max_ntau_2 = floor(pLP[colourRefAlt]->w0/pLP[colourRef]->w0*max_ntau)+1;
  for(int nl=0; nl<Nl; nl++)  one_label.push_back(0);
  orbit_labels.resize(t_saddle.size(), one_label);
  //  const int nEnRef = (int)(t_orbits[0].size()*0.9);  // 90% of way from E_upper to E_lower
  const int nEnRef = (int)(t_orbits[0].size()-10);  // just shy of full way from E_upper to E_lower
  
  //  const int nEcut = t_orbits[0].size()/2;
  double E_cutoff; 
  if(pLaser->colours>1)  E_cutoff = 1.32*pJet->Ip + 3.17*(pLP[0]->I0+pLP[1]->I0)/(4.0*pow((pLP[0]->w0*pLP[0]->I0+pLP[1]->w0*pLP[1]->I0)/(pLP[0]->I0+pLP[1]->I0),2));
  else  E_cutoff = 1.32*pJet->Ip + 3.17*pLP[0]->I0/(4.0*pow(pLP[0]->w0,2));

  const int nEcut = (int)(0.5*(En[0]-E_cutoff)/dEn);  

  for(int ns=0; ns<t_orbits.size(); ns++) {
    for(int nl=0; nl<Nl; nl++)  one_label[nl] = -100;
    for(int nion=min_ntp; nion<=max_ntp+max_ntau; nion++)
      if(real(t_orbits[ns][0])>=t0_1+t_shift-0.4*T_1+nion*(0.5*T_1) && real(t_orbits[ns][0])<t0_1+t_shift+0.4*T_1+nion*(0.5*T_1)) {
	one_label[0] = nion;
	break;
      }
    for(int nion2=min_ntp_2; nion2<=max_ntp_2+max_ntau_2; nion2++)
      //if(real(t_orbits[ns][0])>=t0_2+t_shift_2-0.4*T_2+nion2*(0.5*T_2) && real(t_orbits[ns][0])<t0_2+t_shift_2+0.4*T_2+nion2*(0.5*T_2)) {
      if(real(t_orbits[ns][0])>=t0_2+t_shift_2-0.5*T_2+nion2*(0.5*T_2) && real(t_orbits[ns][0])<t0_2+t_shift_2+0.5*T_2+nion2*(0.5*T_2)) {
	one_label[2] = nion2;
	break;
      }
    tau = t_orbits[ns][0]-tp_orbits[ns][0];
    if(real(tau)<0.45*T_1)  one_label[1] = 0;
    for(int ntau=1; ntau<=max_ntau; ntau++)
      if(real(tau)>=0.45*T_1+(ntau-1)*(0.5*T_1) && real(tau)<0.45*T_1+ntau*(0.5*T_1)) {
	one_label[1] = ntau;
	break;
      }
    if(real(t_orbits[ns][nEcut]-tp_orbits[ns][nEcut])<0.5*T_1) one_label[1] = 0;  // possible correction
    if(real(tau)<0.45*T_2)  one_label[3] = 0;
    for(int ntau2=1; ntau2<=max_ntau_2; ntau2++)
      if(real(tau)>=0.45*T_2+(ntau2-1)*(0.5*T_2) && real(tau)<0.45*T_2+ntau2*(0.5*T_2)) {
	one_label[3] = ntau2;
	break;
      }
    if(real(t_orbits[ns][nEcut]-tp_orbits[ns][nEcut])<0.5*T_2) one_label[3] = 0;  // possible correction
    orbit_labels[ns] = one_label;
  }
  for(int ns=0; ns<t_orbits.size(); ns++)
    for(int ns2=ns+1; ns2<t_orbits.size(); ns2++)
      if(orbit_labels[ns2][0]==orbit_labels[ns][0] && orbit_labels[ns2][1]==orbit_labels[ns][1] && orbit_labels[ns2][2]==orbit_labels[ns][2] && orbit_labels[ns2][3]==orbit_labels[ns][3]) {
	if(real(t_orbits[ns][nEnRef]-tp_orbits[ns][nEnRef])<real(t_orbits[ns2][nEnRef]-tp_orbits[ns2][nEnRef])) {orbit_labels[ns][4] = -1;  orbit_labels[ns2][4] = 1;}
	else {orbit_labels[ns][4] = 1;  orbit_labels[ns2][4] = -1;}
      }
  for(int ns=0; ns<orbit_labels.size(); ns++) {
    if(imag(S(t_orbits[ns][10], tp_orbits[ns][10], En[10]))<imag(S(t_orbits[ns][0], tp_orbits[ns][0], En[0])))  orbit_labels[ns][5] = -1;
    else  orbit_labels[ns][5] = 1;
  }
  double T_min = 0.01*std::max(T_1,T_2); //0.05*std::max(T_1,T_2);
  for(int ns=0; ns<t_orbits.size(); ns++)
    for(int ns2=ns+1; ns2<t_orbits.size(); ns2++)
/*       if(orbit_labels[ns2][0]==orbit_labels[ns][0] && orbit_labels[ns2][1]==orbit_labels[ns][1] && orbit_labels[ns2][2]==orbit_labels[ns][2] && orbit_labels[ns2][3]==orbit_labels[ns][3] && orbit_labels[ns][6]==-100 && orbit_labels[ns2][6]==-100) { */

      if(orbit_labels[ns2][0]==orbit_labels[ns][0] && orbit_labels[ns2][1]==orbit_labels[ns][1] && orbit_labels[ns2][2]==orbit_labels[ns][2] && orbit_labels[ns2][3]==orbit_labels[ns][3] && orbit_labels[ns2][4]==-orbit_labels[ns][4] && orbit_labels[ns2][5]==-orbit_labels[ns][5] && orbit_labels[ns][6]==-100 && orbit_labels[ns2][6]==-100) {


	if(fabs(real(t_orbits[ns][nEnRef]-t_orbits[ns][0]))<T_min) {orbit_labels[ns][6] = 100000;  orbit_labels[ns2][6] = 100000;}
	else if(fabs(real(t_orbits[ns2][nEnRef]-tp_orbits[ns2][nEnRef])-real(t_orbits[ns][nEnRef]-tp_orbits[ns][nEnRef]))<=inc_frac*T_1) {orbit_labels[ns][6] = 200000;  orbit_labels[ns2][6] = 200000;}
	else {orbit_labels[ns][6] = 1;  orbit_labels[ns2][6] = 1;}
	break;
      }
  for(int ns=0; ns<t_orbits.size(); ns++)
    if(real(t_orbits[ns][0]-tp_orbits[ns][0])<0.1*std::max(T_1,T_2))  orbit_labels[ns][6] = 300000; 
}

// Calculates population vector along real time axis
void qo_response_1d::Calculate_Population() {
  const double min_t = std::min(real(tp_orbits[0][NEn-1]),real(tp_orbits[1][NEn-1]));
  const double max_t = std::max(real(t_orbits[t_orbits.size()-2][NEn-1]),real(t_orbits[t_orbits.size()-1][NEn-1]));
  const int Nt = 2*NEn; //1024;
  const double dt = (max_t-min_t)/Nt;
  double t;
  population.clear();  population.resize(Nt);
  double Eabs, ADK; 
  double exp_factor = 0.0;
  for(int nt=0; nt<Nt; nt++) {
    t = min_t+nt*dt;
    exp_factor -= dt*pJet->ADK_Rate(E(t));
    population[nt] = exp(exp_factor);
  }
}

// Function to fit t,tp onto calculated population axis
void qo_response_1d::Fit_Population(double &pop_t, double &pop_tp, dcmplx t, dcmplx tp, double min_t, double dt) {
  double t_frac, t_int;
  t_frac = modf((real(t)-min_t)/dt, &t_int);
  if(t_int<0.0)  pop_t = population[0];
  else if(t_int>=population.size()-1)  pop_t = population[population.size()-1];
  else  pop_t = (1.0-t_frac)*population[t_int]+t_frac*population[t_int+1];
  t_frac = modf((real(tp)-min_t)/dt, &t_int);
  if(t_int<0.0)  pop_tp = population[0];
  else if(t_int>=population.size()-1)  pop_tp = population[population.size()-1];
  else  pop_tp = (1.0-t_frac)*population[t_int]+t_frac*population[t_int+1];
} 

// Dipole transition matrix element
dcmplx qo_response_1d::DTME_Prefactor() {
  dcmplx pre_dtme; 
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
  return pre_dtme;
}

// Dipole transition matrix element
dcmplx qo_response_1d::DTME(dcmplx k, dcmplx pre_dtme) {
  dcmplx element;
  if(pJet->n==1)
    element = pre_dtme*k/pow(k*k+2.0*pJet->Ip, 3);
  else if(pJet->n==2)
    if(pJet->l==0)
      element = pre_dtme*k*(k*k-4.0*pJet->Ip)/pow(k*k+2.0*pJet->Ip, 4);
    else if(pJet->l==1)
      element = pre_dtme*(5.0*k*k-2.0*pJet->Ip)/pow(k*k+2.0*pJet->Ip, 4);
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

// Dipole transition matrix element without k^2 in denominator (crude way to remove singlarity)
dcmplx qo_response_1d::DTME_singular(dcmplx k, dcmplx pre_dtme) {
  dcmplx element;
  if(pJet->n==1)
    element = pre_dtme*k/pow(2.0*pJet->Ip, 3);  // must have k factor for polarity purposes (to give odd- rather than even-harmonics)
else if(pJet->n==2)
    if(pJet->l==0)
      element = pre_dtme*k*(k*k-4.0*pJet->Ip)/pow(2.0*pJet->Ip, 4);
    else if(pJet->l==1)
      element = pre_dtme*(5.0*k*k-2.0*pJet->Ip)/pow(2.0*pJet->Ip, 4);
  else if(pJet->n==3)
    if(pJet->l==0)
      element = pre_dtme*(3.0*pow(k, 4)-18.0*k*k*2.0*pJet->Ip+11.0*pow(2.0*pJet->Ip, 2))/pow(2.0*pJet->Ip, 5);
    else if(pJet->l==1)
      element = pre_dtme*(2.0*k*k*(3.0*k*k-10.0*pJet->Ip)-(pow(k, 4)-pow(2.0*pJet->Ip, 2)))/pow(2.0*pJet->Ip, 5);
    else
      element = pre_dtme*0.0;
  else 
    element = pre_dtme*0.0;     
  return element;
}

// Calculates harmonic spectrum resulting from previously specified orbits using the Saddle-Point Approx
void qo_response_1d::Dipole_Acceleration_SaddlePoint() { 
  bool found_pair;
  int ns_pair;
  dcmplx t, tp;
  double epsilon = 1.0e-10;
  double Energy; 
  intvec ns_store;
  dcmplx DetHess, DetHess_prev[2];
  dcmplx spread, spread_prev[2];
  double sign_spreadHess[2];
  dcmplx spreadHess; 
  double min_t = std::min(real(tp_orbits[0][NEn-1]),real(tp_orbits[1][NEn-1]));
  double max_t = std::max(real(t_orbits[t_orbits.size()-2][NEn-1]),real(t_orbits[t_orbits.size()-1][NEn-1]));
  double dt = (max_t-min_t)/population.size();
  double pop_t, pop_tp;
  dcmplx cmplxAmp;
  dcmplx action;
  for(int nEn=0; nEn<NEn; nEn++)  DA_local[nEn] = 0.0;
  dcmplx pre_dtme = DTME_Prefactor();
  dcmplx dtme[2];
  for(int ns=0; ns<t_orbits.size(); ns++) {
    if(traj=='A' || ((traj=='S'||traj=='B') && orbit_labels[ns][4]==-1) || ((traj=='L'||traj=='B') && orbit_labels[ns][4]==1)) {
      ns_store.push_back(ns);  // records position in orbit_labels
      sign_spreadHess[0] = 1.0;  sign_spreadHess[1] = 1.0;
      for(int nEn=0; nEn<t_orbits[0].size(); nEn++) {
	Energy = En[nEn];
	t = t_orbits[ns][nEn];
	tp = tp_orbits[ns][nEn];
	DetHess = Determinant_Hessian(t,tp);
	spread = epsilon+i*(t-tp);
	if(nEn>0) {
	  if(abs(arg(DetHess)-arg(DetHess_prev[0]))>pi)  sign_spreadHess[0] *= -1.0;
	  if(abs(arg(spread)-arg(spread_prev[0]))>pi)  sign_spreadHess[0] *= -1.0;
	}
	DetHess_prev[0] = DetHess;  spread_prev[0] = spread;
	spreadHess = sign_spreadHess[0]/std::polar(sqrt(pow(abs(spread),3)*abs(DetHess)), 0.5*(3.0*arg(spread)+arg(DetHess)));
	//spreadHess = sign_spreadHess[0]/pow(spread,1.5)/sqrt(DetHess);
	Fit_Population(pop_t, pop_tp, t, tp, min_t, dt);
	cmplxAmp = i*2.0*pow(pi*Energy,2)*spreadHess*(pop_t*conj(DTME(p(t,tp)+A(t),pre_dtme)))*(pop_tp*E(tp)*DTME_singular(p(t,tp)+A(tp),pre_dtme));	  
	action = S(t, tp, Energy);
 	DA_local[nEn] += conj(cmplxAmp*exp(i*action));  // conj() as was previously giving time-reversed emission 
      }
    }
  }     
}

// Calculates harmonic spectrum resulting from previously specified orbits using the Uniform Approx
void qo_response_1d::Dipole_Acceleration_Uniform() { 
  bool found_pair;
  int ns_pair;
  dcmplx t, tp;
  double epsilon = 1.0e-10;
  dcmplx DetHess, DetHess_prev[2];
  dcmplx spread, spread_prev[2];
  double sign_spreadHess[2];
  dcmplx spreadHess[2];
  dcmplx pre_dtme = DTME_Prefactor();
  dcmplx dtmes[2];
  double Energy;
  double min_t = std::min(real(tp_orbits[0][NEn-1]),real(tp_orbits[1][NEn-1]));
  double max_t = std::max(real(t_orbits[t_orbits.size()-2][NEn-1]),real(t_orbits[t_orbits.size()-1][NEn-1]));
  double dt = (max_t-min_t)/population.size();
  double pop_t, pop_tp;
  dcmplx cmplxAmp[2], cmplxAmp_plus, cmplxAmp_minus;
  double short_shift=0.0, long_shift=0.0;
  dcmplx action[2], action_plus, action_minus;
  double sign_spectrum, eta;
  dcmplx sqrt_z, minus_z;
  dcmplx Airy, dzAiry;
  int ierr = 0, argflag = 0, modflag = 0;
  dcmplx pre_factor = sqrt(6.0*pi)*exp(i*0.25*pi);
  double diff_im_S, prev_diff_im_S;
  bool antistokes_passed;

  int nEn_1, nEn_2;
  double x1, x2, y1, y2, k1, k2, a2, b2;
  double x_frac;

  spectrum_pair.clear();
  spectrum_pair.resize(NEn);

  for(int nEn=0; nEn<NEn; nEn++)  DA_local[nEn] = 0.0;

  for(int ns=0; ns<t_orbits.size(); ns++) {
    if(orbit_labels[ns][4]==-1 && orbit_labels[ns][6]==1) {
    //if(orbit_labels[ns][4]==-1) {
      found_pair = false;
      for(ns_pair=0; ns_pair<t_orbits.size(); ns_pair++)
	if(orbit_labels[ns_pair][0]==orbit_labels[ns][0] && orbit_labels[ns_pair][1]==orbit_labels[ns][1] && orbit_labels[ns_pair][2]==orbit_labels[ns][2] && orbit_labels[ns_pair][3]==orbit_labels[ns][3] && orbit_labels[ns_pair][4]==1 && orbit_labels[ns_pair][6]==1) {
	found_pair = true;
	  break;
	}
      if(found_pair==false) {
	std::cout << "Failed to find pair for orbit " << orbit_labels[ns][0] << ',' << orbit_labels[ns][1] << ',' << orbit_labels[ns][2] << ',' << orbit_labels[ns][3] << ',' << orbit_labels[ns][4] << '\n';
	return;
      }
    }
    else  continue;  // skips back to start of for-loop and advances ns
    eta = (double)(abs(orbit_labels[ns][4]+orbit_labels[ns][5]))-1.0;
    sign_spreadHess[0] = 1.0;  sign_spreadHess[1] = 1.0;
    prev_diff_im_S = 0.0;  antistokes_passed = false;
    for(int nEn=0; nEn<NEn; nEn++) {
      Energy = En[nEn];
      // short orbit
      t = t_orbits[ns][nEn];
      tp = tp_orbits[ns][nEn];
      dtmes[0] = conj(DTME(p(t,tp)+A(t),pre_dtme))*DTME_singular(p(t,tp)+A(tp),pre_dtme);
      DetHess = Determinant_Hessian(t,tp);
      spread = epsilon+i*(t-tp);
      if(nEn>0) {
	if(abs(arg(DetHess)-arg(DetHess_prev[0]))>pi)  sign_spreadHess[0] *= -1.0;
	if(abs(arg(spread)-arg(spread_prev[0]))>pi)  sign_spreadHess[0] *= -1.0;
      }
      DetHess_prev[0] = DetHess;  spread_prev[0] = spread;
      spreadHess[0] = sign_spreadHess[0]/std::polar(sqrt(pow(abs(spread),3)*abs(DetHess)), 0.5*(3.0*arg(spread)+arg(DetHess)));
      //spreadHess[0] = sign_spreadHess[0]/pow(spread,1.5)/sqrt(DetHess);
      Fit_Population(pop_t, pop_tp, t, tp, min_t, dt);
      cmplxAmp[0] = i*2.0*pow(pi*Energy,2)*spreadHess[0]*pop_t*dtmes[0]*pop_tp*E(tp);
      action[0] = S(t, tp, Energy);
      // long orbit
      t = t_orbits[ns_pair][nEn];
      tp = tp_orbits[ns_pair][nEn];
      dtmes[1] = conj(DTME(p(t,tp)+A(t),pre_dtme))*DTME_singular(p(t,tp)+A(tp),pre_dtme);
      DetHess = Determinant_Hessian(t,tp);
      spread = epsilon+i*(t-tp);
      if(nEn>0) {
	if(abs(arg(DetHess)-arg(DetHess_prev[1]))>pi)  sign_spreadHess[1] *= -1.0;
	if(abs(arg(spread)-arg(spread_prev[1]))>pi)  sign_spreadHess[1] *= -1.0;
      }
      DetHess_prev[1] = DetHess;  spread_prev[1] = spread;
      spreadHess[1] = sign_spreadHess[1]/std::polar(sqrt(pow(abs(spread),3)*abs(DetHess)), 0.5*(3.0*arg(spread)+arg(DetHess)));
      //spreadHess[1] = sign_spreadHess[1]/pow(spread,1.5)/sqrt(DetHess);
      Fit_Population(pop_t, pop_tp, t, tp, min_t, dt);
      cmplxAmp[1] = i*2.0*pow(pi*Energy,2)*spreadHess[1]*pop_t*dtmes[1]*pop_tp*E(tp);
      action[1] = S(t, tp, Energy);
      if(nEn==0) {
	// potentially a problem here when both short and long orbit 'groups' form >1 band each
	short_shift = 0.0;  long_shift = 0.0;  // reset values
/* 	if(orbit_labels[ns][1]!=1 && orbit_labels[ns][3]!=1) */
/* 	  if(arg(spreadHess[1])-arg(spreadHess[0])>pi)  long_shift -= pi;  // if phase(long)>phase(short)+pi subract pi (but keep phase(long)>phase(short)) */
/* 	  else if(arg(spreadHess[1])-arg(spreadHess[0])<-pi) {short_shift -=pi;  long_shift += pi;}  // if phase(long)+pi<phase(short) subract pi from phase(short) and add pi to phase(long) */
/* 	  else if(arg(spreadHess[1])-arg(spreadHess[0])<0.0)  long_shift += pi;  // if instead phase(long)<phase(short) just add pi to phase(long) */
/* 	//if(orbit_labels[ns][1]==1 || orbit_labels[ns][1]==3)  {short_shift += pi;  long_shift -= pi;} */
/* 	if(arg(spreadHess[1])-arg(spreadHess[0])>pi)  long_shift -= pi;  // if phase(long)>phase(short)+pi subract pi (but keep phase(long)>phase(short)) */
/* 	else if(arg(spreadHess[1])-arg(spreadHess[0])<-pi) {short_shift -=pi;  long_shift += pi;}  // if phase(long)+pi<phase(short) subract pi from phase(short) and add pi to phase(long) */
/* 	else if(arg(spreadHess[1])-arg(spreadHess[0])<0.0)  long_shift += pi;  // if instead phase(long)<phase(short) just add pi to phase(long) */

	if(pLP[0]->I0>1e10*Wcm2_to_au && orbit_labels[ns][1]>1)
	  if(arg(spreadHess[1])-arg(spreadHess[0])>pi)  long_shift -= pi;  // if phase(long)>phase(short)+pi subract pi (but keep phase(long)>phase(short))
	  else if(arg(spreadHess[1])-arg(spreadHess[0])<-pi) {short_shift -=pi;  long_shift += pi;}  // if phase(long)+pi<phase(short) subract pi from phase(short) and add pi to phase(long)
	  else if(arg(spreadHess[1])-arg(spreadHess[0])<0.0)  long_shift += pi;  // if instead phase(long)<phase(short) just add pi to phase(long)

      }
      cmplxAmp[0] *= exp(i*short_shift);
      cmplxAmp[1] *= exp(i*long_shift);
      // uniform approx
      cmplxAmp_plus = 0.5 * (cmplxAmp[0]+i*cmplxAmp[1]);
      cmplxAmp_minus = 0.5 * (cmplxAmp[0]-i*cmplxAmp[1]);
      diff_im_S = imag(action[0])-imag(action[1]);
      action_plus = 0.5 * (action[0]+action[1]);
      action_minus = 0.5 * (action[0]-action[1]);
      if(antistokes_passed==false && ((diff_im_S>0.0&&prev_diff_im_S<0.0) || (diff_im_S<0.0&&prev_diff_im_S>0.0)))  antistokes_passed = true;
      if(antistokes_passed==false) {
	sign_spectrum = -1.0;
	sqrt_z = -pow(1.5*action_minus, 1.0/3.0)*exp(i*eta*pi/3.0);
	minus_z = -pow(sqrt_z, 2);
      }
      else {
	sign_spectrum = 1.0;
	sqrt_z = pow(1.5*action_minus, 1.0/3.0);
	minus_z = -pow(sqrt_z, 2);
      }
      __airy_functions_complex_double_airy_aic(&minus_z, &Airy, &dzAiry, &ierr, &argflag, &modflag);
      spectrum_pair[nEn] = sign_spectrum*pre_factor*sqrt(action_minus)*exp(i*action_plus)*(cmplxAmp_minus/sqrt_z*Airy-i*cmplxAmp_plus/minus_z*dzAiry);
      prev_diff_im_S = diff_im_S;
    }

    nEn_1 = (int)std::min(0.1*eV_to_au/dEn, t_orbits[0].size()-3.0);  // -2eV from E_C
    nEn_2 = (int)std::max(-0.1*eV_to_au/dEn, 3.0);  // +2eV from E_C
/*     nEn_1 = (int)(0.1*eV_to_au/dEn);  */
/*     nEn_2 = 3.0;  // +2eV from E_C */

    //std::cout << nEn_1 << '\t' << nEn_2 << '\t' << (int)(2.0*eV_to_au/dEn) << '\n';

    for(int nEn=nEn_1-1; nEn>nEn_2; nEn--) {
      nEn; //DA_local[nEn]; //spectrum_long[nEn];
    }

    //for(int nEn=nEn_1-1; nEn>=nEn_2; nEn--);

    for(int nEn=0; nEn<NEn; nEn++) {
      if(traj=='B')  DA_local[nEn] += conj(spectrum_pair[nEn]);
    }

/*     if(norm(spectrum_pair[NEn-1])<10.0*norm(spectrum_pair[NEn-31])) */
/*       for(int nEn=0; nEn<NEn; nEn++)  */
/* 	if(traj=='B')  DA_local[nEn] += conj(spectrum_pair[nEn]); */
    
  }
}

// Calculates harmonic spectrum resulting from previously specified orbits using a combination of the Saddle-Point and Uniform Approx
void qo_response_1d::Dipole_Acceleration_Combined() {
  bool found_pair;
  int ns_pair;
  dcmplx t, tp;
  double epsilon = 1.0e-10;
  dcmplx DetHess, DetHess_prev[2];
  dcmplx spread, spread_prev[2];
  double sign_spreadHess[2];
  dcmplx spreadHess[2];
  dcmplx pre_dtme = DTME_Prefactor();
  dcmplx dtmes[2];
  double Energy;
  double min_t = std::min(real(tp_orbits[0][NEn-1]),real(tp_orbits[1][NEn-1]));
  double max_t = std::max(real(t_orbits[t_orbits.size()-2][NEn-1]),real(t_orbits[t_orbits.size()-1][NEn-1]));
  double dt = (max_t-min_t)/population.size();
  double pop_t, pop_tp;
  dcmplx cmplxAmp[2], cmplxAmp_plus, cmplxAmp_minus;
  double short_shift=0.0, long_shift=0.0;
  dcmplx action[2], action_plus, action_minus;
  double sign_spectrum, eta;
  dcmplx sqrt_z, minus_z;
  dcmplx Airy, dzAiry;
  int ierr = 0, argflag = 0, modflag = 0;
  dcmplx pre_factor = sqrt(6.0*pi)*exp(i*0.25*pi);
  double diff_im_S, prev_diff_im_S;
  bool antistokes_passed;
  int nEn_aS;
  double phase_wrap[2], phase_wrap_prev[2], phase_unwrap[2], phase_unwrap_prev[2];
  bool crossing_passed;
  int nEn_C;
  dcmplx shift_factor;
  for(int nEn=0; nEn<NEn; nEn++)  DA_local[nEn] = 0.0;
  for(int ns=0; ns<t_orbits.size(); ns++) {
    if(orbit_labels[ns][4]==-1 && orbit_labels[ns][6]==1) {
    //if(orbit_labels[ns][4]==-1) {
      found_pair = false;
      for(ns_pair=0; ns_pair<t_orbits.size(); ns_pair++)
	if(orbit_labels[ns_pair][0]==orbit_labels[ns][0] && orbit_labels[ns_pair][1]==orbit_labels[ns][1] && orbit_labels[ns_pair][2]==orbit_labels[ns][2] && orbit_labels[ns_pair][3]==orbit_labels[ns][3] && orbit_labels[ns_pair][4]==1 && orbit_labels[ns_pair][6]==1) {
	/* if(orbit_labels[ns_pair][0]==orbit_labels[ns][0] && orbit_labels[ns_pair][1]==orbit_labels[ns][1] && orbit_labels[ns_pair][4]==1) { */
	found_pair = true;
	  break;
	}
      if(found_pair==false) {
	std::cout << "Failed to find pair for orbit " << orbit_labels[ns][0] << ',' << orbit_labels[ns][1] << ',' << orbit_labels[ns][2] << ',' << orbit_labels[ns][3] << ',' << orbit_labels[ns][4] << '\n';
	return;
      }
    }
    else  continue;  // skips back to start of for-loop and advances ns
    eta = (double)(abs(orbit_labels[ns][4]+orbit_labels[ns][5]))-1.0;
    sign_spreadHess[0] = 1.0;  sign_spreadHess[1] = 1.0;
    prev_diff_im_S = 0.0;  antistokes_passed = false;
    phase_wrap_prev[0] = 0.0;  phase_wrap_prev[1] = 0.0;  phase_wrap[0] = 0.0;  phase_wrap[1] = 0.0;
    phase_unwrap_prev[0] = 0.0;  phase_unwrap_prev[1] = 0.0;  phase_unwrap[0] = 0.0;  phase_unwrap[1] = 0.0;
    crossing_passed = false;
    nEn_aS = 0;  nEn_C = t_orbits[0].size()-3; //nE_C = 0;
    for(int nEn=0; nEn<NEn; nEn++) {
      Energy = En[nEn];
      // short orbit
      t = t_orbits[ns][nEn];
      tp = tp_orbits[ns][nEn];
      dtmes[0] = conj(DTME(p(t,tp)+A(t),pre_dtme))*DTME_singular(p(t,tp)+A(tp),pre_dtme);
      DetHess = Determinant_Hessian(t,tp);
      spread = epsilon+i*(t-tp);
      if(nEn>0) {
	if(abs(arg(DetHess)-arg(DetHess_prev[0]))>pi)  sign_spreadHess[0] *= -1.0;
	if(abs(arg(spread)-arg(spread_prev[0]))>pi)  sign_spreadHess[0] *= -1.0;
      }
      DetHess_prev[0] = DetHess;  spread_prev[0] = spread;
      spreadHess[0] = sign_spreadHess[0]/std::polar(sqrt(pow(abs(spread),3)*abs(DetHess)), 0.5*(3.0*arg(spread)+arg(DetHess)));
      //spreadHess[0] = sign_spreadHess[0]/pow(spread,1.5)/sqrt(DetHess);
      Fit_Population(pop_t, pop_tp, t, tp, min_t, dt);
      cmplxAmp[0] = i*2.0*pow(pi*Energy,2)*spreadHess[0]*pop_t*dtmes[0]*pop_tp*E(tp);
      action[0] = S(t, tp, Energy);
      // long orbit
      t = t_orbits[ns_pair][nEn];
      tp = tp_orbits[ns_pair][nEn];
      dtmes[1] = conj(DTME(p(t,tp)+A(t),pre_dtme))*DTME_singular(p(t,tp)+A(tp),pre_dtme);
      DetHess = Determinant_Hessian(t,tp);
      spread = epsilon+i*(t-tp);
      if(nEn>0) {
	if(abs(arg(DetHess)-arg(DetHess_prev[1]))>pi)  sign_spreadHess[1] *= -1.0;
	if(abs(arg(spread)-arg(spread_prev[1]))>pi)  sign_spreadHess[1] *= -1.0;
      }
      DetHess_prev[1] = DetHess;  spread_prev[1] = spread;
      spreadHess[1] = sign_spreadHess[1]/std::polar(sqrt(pow(abs(spread),3)*abs(DetHess)), 0.5*(3.0*arg(spread)+arg(DetHess)));
      //spreadHess[1] = sign_spreadHess[1]/pow(spread,1.5)/sqrt(DetHess);
      Fit_Population(pop_t, pop_tp, t, tp, min_t, dt);
      cmplxAmp[1] = i*2.0*pow(pi*Energy,2)*spreadHess[1]*pop_t*dtmes[1]*pop_tp*E(tp);
      action[1] = S(t, tp, Energy);
      if(nEn==0) {
	// potentially a problem here when both short and long orbit 'groups' form >1 band each
	short_shift = 0.0;  long_shift = 0.0;  // reset values
/* 	if(orbit_labels[ns][1]!=1 && orbit_labels[ns][3]!=1) */
/* 	  if(arg(spreadHess[1])-arg(spreadHess[0])>pi)  long_shift -= pi;  // if phase(long)>phase(short)+pi subract pi (but keep phase(long)>phase(short)) */
/* 	  else if(arg(spreadHess[1])-arg(spreadHess[0])<-pi) {short_shift -=pi;  long_shift += pi;}  // if phase(long)+pi<phase(short) subract pi from phase(short) and add pi to phase(long) */
/* 	  else if(arg(spreadHess[1])-arg(spreadHess[0])<0.0)  long_shift += pi;  // if instead phase(long)<phase(short) just add pi to phase(long) */
/* 	//if(orbit_labels[ns][1]==1 || orbit_labels[ns][1]==3)  {short_shift += pi;  long_shift -= pi;} */
/* 	if(arg(spreadHess[1])-arg(spreadHess[0])>pi)  long_shift -= pi;  // if phase(long)>phase(short)+pi subract pi (but keep phase(long)>phase(short)) */
/* 	else if(arg(spreadHess[1])-arg(spreadHess[0])<-pi) {short_shift -=pi;  long_shift += pi;}  // if phase(long)+pi<phase(short) subract pi from phase(short) and add pi to phase(long) */
/* 	else if(arg(spreadHess[1])-arg(spreadHess[0])<0.0)  long_shift += pi;  // if instead phase(long)<phase(short) just add pi to phase(long) */

	if(pLP[0]->I0>1e10*Wcm2_to_au && orbit_labels[ns][1]>1)
	  if(arg(spreadHess[1])-arg(spreadHess[0])>pi)  long_shift -= pi;  // if phase(long)>phase(short)+pi subract pi (but keep phase(long)>phase(short))
	  else if(arg(spreadHess[1])-arg(spreadHess[0])<-pi) {short_shift -=pi;  long_shift += pi;}  // if phase(long)+pi<phase(short) subract pi from phase(short) and add pi to phase(long)
	  else if(arg(spreadHess[1])-arg(spreadHess[0])<0.0)  long_shift += pi;  // if instead phase(long)<phase(short) just add pi to phase(long)
      }
      cmplxAmp[0] *= exp(i*short_shift);
      cmplxAmp[1] *= exp(i*long_shift);
      spectrum_short[nEn] = cmplxAmp[0]*exp(i*action[0]);
      spectrum_long[nEn] = cmplxAmp[1]*exp(i*action[1]);
      // uniform approx
      cmplxAmp_plus = 0.5 * (cmplxAmp[0]+i*cmplxAmp[1]);
      cmplxAmp_minus = 0.5 * (cmplxAmp[0]-i*cmplxAmp[1]);
      diff_im_S = imag(action[0])-imag(action[1]);
      action_plus = 0.5 * (action[0]+action[1]);
      action_minus = 0.5 * (action[0]-action[1]);
      if(antistokes_passed==false && ((diff_im_S>0.0&&prev_diff_im_S<0.0) || (diff_im_S<0.0&&prev_diff_im_S>0.0))) {
	antistokes_passed = true;
	nEn_aS = nEn;
      }
      if(antistokes_passed==false && ((diff_im_S>0.0&&prev_diff_im_S<0.0) || (diff_im_S<0.0&&prev_diff_im_S>0.0)))  antistokes_passed = true;
      if(antistokes_passed==false) {
	sign_spectrum = -1.0;
	sqrt_z = -pow(1.5*action_minus, 1.0/3.0)*exp(i*eta*pi/3.0);
	minus_z = -pow(sqrt_z, 2);
      }
      else {
	sign_spectrum = 1.0;
	sqrt_z = pow(1.5*action_minus, 1.0/3.0);
	minus_z = -pow(sqrt_z, 2);
      }
      __airy_functions_complex_double_airy_aic(&minus_z, &Airy, &dzAiry, &ierr, &argflag, &modflag);
      spectrum_pair[nEn] = sign_spectrum*pre_factor*sqrt(action_minus)*exp(i*action_plus)*(cmplxAmp_minus/sqrt_z*Airy-i*cmplxAmp_plus/minus_z*dzAiry);
      prev_diff_im_S = diff_im_S;

      // find crossing - energy at which final hump of half-cycle spectra is at a maximum
      if(crossing_passed==false) {
	phase_wrap_prev[0] = phase_wrap[0];  phase_unwrap_prev[0] = phase_unwrap[0];
  	phase_wrap[0] = arg(spectrum_short[nEn]);
	if(nEn==0) {phase_unwrap[0] = phase_wrap[0];  phase_unwrap_prev[0] = 0.0;}
	else if(phase_wrap[0]-phase_wrap_prev[0]<-pi)  phase_unwrap[0] = phase_unwrap_prev[0]+(phase_wrap[0]+2.0*pi)-phase_wrap_prev[0];
	else if(phase_wrap[0]-phase_wrap_prev[0]>pi)  phase_unwrap[0] = phase_unwrap_prev[0]+(phase_wrap[0]-2.0*pi)-phase_wrap_prev[0];
	else  phase_unwrap[0] = phase_unwrap_prev[0]+phase_wrap[0]-phase_wrap_prev[0];
	phase_wrap_prev[1] = phase_wrap[1];  phase_unwrap_prev[1] = phase_unwrap[1];
	phase_wrap[1] = arg(spectrum_long[nEn]);
	if(nEn==0) {
	  phase_unwrap[1] = phase_wrap[1];  phase_unwrap_prev[1] = 0.0;
/* 	  if(phase_wrap[1]-phase_wrap[0]>pi)  phase_unwrap[1] -= 2.0*pi; */
/* 	  else if (phase_wrap[1]-phase_wrap[0]<-pi)  phase_unwrap[1] += 2.0*pi */;
	  if(phase_wrap[1]<phase_wrap[0])  phase_unwrap[1] += 2.0*pi;
	}
	else if(phase_wrap[1]-phase_wrap_prev[1]<-pi)  phase_unwrap[1] = phase_unwrap_prev[1]+(phase_wrap[1]+2.0*pi)-phase_wrap_prev[1];
	else if(phase_wrap[1]-phase_wrap_prev[1]>pi)  phase_unwrap[1] = phase_unwrap_prev[1]+(phase_wrap[1]-2.0*pi)-phase_wrap_prev[1];
	else  phase_unwrap[1] = phase_unwrap_prev[1]+phase_wrap[1]-phase_wrap_prev[1];
	//if((phase_unwrap[1]-phase_unwrap[0]>0.0&&phase_unwrap_prev[1]-phase_unwrap_prev[0]<0.0) || (phase_unwrap[1]-phase_unwrap[0]<0.0&&phase_unwrap_prev[1]-phase_unwrap_prev[0]>0.0)) {

	if(antistokes_passed==true && ((phase_unwrap[1]-phase_unwrap[0]>0.0&&phase_unwrap_prev[1]-phase_unwrap_prev[0]<0.0) || (phase_unwrap[1]-phase_unwrap[0]<0.0&&phase_unwrap_prev[1]-phase_unwrap_prev[0]>0.0))) {
	  crossing_passed = true;
// 	    if(fabs(phase_unwrap[1]-phase_unwrap[0])<fabs(phase_unwrap_prev[1]-phase_unwrap_prev[0]))  nEn_C = nEn;
// 	    else  nEn_C = nEn-1;
	  nEn_C = nEn;
	}
      }
      // Above the crossing, phase of both average of SPA short&long gradually shifts away from UA. This removes the shift by comparing SPA average to UA.
      if(crossing_passed==false) {
	spectrum_short[nEn] *= exp(i*(arg(spectrum_pair[nEn])-0.5*(phase_unwrap[0]+phase_unwrap[1])));
	spectrum_long[nEn] *= exp(i*(arg(spectrum_pair[nEn])-0.5*(phase_unwrap[0]+phase_unwrap[1])));
      }
    }
    //std::cout << "nEn_C = " << nEn_C << ", energy = " << En[nEn_C]*27.212 << "eV\n";
    double pair_shift;
    for(int nEn=0; nEn<nEn_C; nEn++) {
      // At the crossing, phase of average of SPA short&long not quite equal to UA. This performs slight adjustment to match them up.
      spectrum_short[nEn] *= exp(i*(arg(spectrum_short[nEn_C+1]+spectrum_long[nEn_C+1])-arg(spectrum_pair[nEn_C+1])));
      spectrum_long[nEn] *= exp(i*(arg(spectrum_short[nEn_C+1]+spectrum_long[nEn_C+1])-arg(spectrum_pair[nEn_C+1])));
      // Adjust the magnitude of SPA short&long above crossing to match UA (with universal scaling factor taken at crossing) -> derivative wrt E discontinuous at crossing
      spectrum_long[nEn] = abs(spectrum_long[nEn_C])/abs(spectrum_pair[nEn_C]) * abs(spectrum_pair[nEn]) * exp(i*arg(spectrum_long[nEn]));
      spectrum_short[nEn] = abs(spectrum_short[nEn_C])/abs(spectrum_pair[nEn_C]) * abs(spectrum_pair[nEn]) * exp(i*arg(spectrum_short[nEn]));
      // Further adjust both together such that magnitude of their sum matches with UA at the crossing -> s/l spectra discontinuous at crossing (along with their derivatives)
      pair_shift = abs(spectrum_pair[nEn])/abs(spectrum_short[nEn]+spectrum_long[nEn]);
      spectrum_long[nEn] *= pair_shift;
      spectrum_short[nEn] *= pair_shift;
    }
    // Smooth magnitude (and also its derivative wrt energy) using spline interpolation

    double dEn = En[1]-En[2];
    int nEn_1 = (int)std::min(nEn_C+2.0*eV_to_au/dEn, t_orbits[0].size()-3.0);  // -2eV from E_C
    int nEn_2 = (int)std::max(nEn_C-2.0*eV_to_au/dEn, 3.0);  // +2eV from E_C
    //if(nEn_2<1)  std::cout << "\nnEn_1 = " << nEn_1 << ", nEn_2 = " << nEn_2 << " on saddle " << ns << " of " << t_saddle.size() << '\n';

//     int nEn_1 = nEn_C+1;
//     int nEn_2 = nEn_aS-1;
    double x1 = En[nEn_1];
    double x2 = En[nEn_2];
    // short orbit
    double y1 = abs(spectrum_short[nEn_1]);
    double y2 = abs(spectrum_short[nEn_2]);
    double k1 = (abs(spectrum_short[nEn_1-1])-y1)/(En[nEn_1-1]-x1);
    double k2 = (abs(spectrum_short[nEn_2-1])-y2)/(En[nEn_2-1]-x2);
    double a2 = k1*(x2-x1)-(y2-y1);
    double b2 = -k2*(x2-x1)+(y2-y1);
    double x_frac;
    for(int nEn=nEn_1-1; nEn>nEn_2; nEn--) {
      x_frac = (En[nEn]-x1)/(x2-x1);
      spectrum_short[nEn] *= ((1-x_frac)*y1 + x_frac*y2 + x_frac*(1.0-x_frac)*(a2*(1.0-x_frac)+b2*x_frac))/abs(spectrum_short[nEn]);
    }
    // long orbit
    y1 = abs(spectrum_long[nEn_1]);
    y2 = abs(spectrum_long[nEn_2]);
    k1 = (abs(spectrum_long[nEn_1-1])-y1)/(En[nEn_1-1]-x1);
    k2 = (abs(spectrum_long[nEn_2-1])-y2)/(En[nEn_2-1]-x2);
    a2 = k1*(x2-x1)-(y2-y1);
    b2 = -k2*(x2-x1)+(y2-y1);
    for(int nEn=nEn_1-1; nEn>nEn_2; nEn--) {
      x_frac = (En[nEn]-x1)/(x2-x1);
      spectrum_long[nEn] *= ((1-x_frac)*y1 + x_frac*y2 + x_frac*(1.0-x_frac)*(a2*(1.0-x_frac)+b2*x_frac))/abs(spectrum_long[nEn]);
    }
    for(int nEn=0; nEn<NEn; nEn++) {
      if(traj=='S'||traj=='B'||traj=='A')  DA_local[nEn] += conj(spectrum_short[nEn]);  // conj() as was previously giving time-reversed emission
      if(traj=='L'||traj=='B'||traj=='A')  DA_local[nEn] += conj(spectrum_long[nEn]);
    }

/*     for(int nEn=0; nEn<NEn; nEn++) { */
/*       if(traj=='B')  DA_local[nEn] += conj(spectrum_short[nEn]); */
/*     } */

  }
}

// Function to translate calculated spectrum onto w axis specified in pT object
void qo_response_1d::Fit_to_wAxis() {
  double wEn2, wEn1, wEn0, wEnm1; 
  for(int nw=0; nw<=pT->Nt/2; nw++) {
    if(pT->w[nw]>En[1])  DA[nw] = 0.0;  // remember En goes high->low
    else if(pT->w[nw]<En[NEn-1])  DA[nw] = 0.0;
    else if(pT->w[nw]<En[NEn-2])  DA[nw] = DA_local[NEn-1]*(En[NEn-2]-pT->w[nw])/(En[NEn-2]-En[NEn-1])+DA_local[NEn-2]*(pT->w[nw]-En[NEn-1])/(En[NEn-2]-En[NEn-1]);
    else if(pT->w[nw]<En[NEn-3])  DA[nw] = DA_local[NEn-2]*(En[NEn-3]-pT->w[nw])/(En[NEn-3]-En[NEn-2])+DA_local[NEn-3]*(pT->w[nw]-En[NEn-2])/(En[NEn-3]-En[NEn-2]);
    else     
      for(int nEn=NEn-3; nEn>=1; nEn--)
	if(pT->w[nw]>En[nEn+1] && pT->w[nw]<=En[nEn]) {
	  wEn2 = pT->w[nw]-En[nEn+2];
	  wEn1 = pT->w[nw]-En[nEn+1];
	  wEn0 = pT->w[nw]-En[nEn];
	  wEnm1 = pT->w[nw]-En[nEn-1];
	  // 4 point (Lagrange-)polynomial interpolation
	  DA[nw] = DA_local[nEn+2]*wEn1/(En[nEn+2]-En[nEn+1])*wEn0/(En[nEn+2]-En[nEn])*wEnm1/(En[nEn+2]-En[nEn-1]);
	  DA[nw] += DA_local[nEn+1]*wEn2/(En[nEn+1]-En[nEn+2])*wEn0/(En[nEn+1]-En[nEn])*wEnm1/(En[nEn+1]-En[nEn-1]);
	  DA[nw] += DA_local[nEn]*wEn2/(En[nEn]-En[nEn+2])*wEn1/(En[nEn]-En[nEn+1])*wEnm1/(En[nEn]-En[nEn-1]);
	  DA[nw] += DA_local[nEn-1]*wEn2/(En[nEn-1]-En[nEn+2])*wEn1/(En[nEn-1]-En[nEn+1])*wEn0/(En[nEn-1]-En[nEn]);
	  break;
	}
  }
  for(int nw=pT->Nt/2+1; nw<pT->Nt; nw++)  DA[nw] = conj(DA[pT->Nt-nw]);
}

// Governing interface to calculate dipole response according to specified method
void qo_response_1d::Dipole_Acceleration(dcmplx *Et) {
  if(analyse==true)  pLaser->Analyse_Field(pLP, Et, findAll);
  Calculate_FieldFactors();
  Set_Bounds();
  Build_Energy_Axis();
  Saddlepoint_Search();
  bool tracked = Track_Orbits();
  if(tracked==false) {
    std::cout << "Orbit tracking failed. Exiting...\n";
    return;//exit(1);
  }
/*   if(t_saddle.size()==0) std::cout << "Catastrophic error here 1 !!!!\n";
/*   if(t_orbits.size()==0) std::cout << "Catastrophic error here 2 !!!!\n"; */
/*   if(t_orbits.size()==0) std::cout << pLP[0]->I0*au_to_Wcm2 << '\t' << pLP[0]->lambda*au_to_m*1e9 << '\t' << pLP[0]->t0*au_to_s*1e15 << '\t' << pLP[0]->fwhm*au_to_s*1e15 << '\t' << pLP[0]->cep << '\n'; */
/*   // changed findAll to n */
  Characterise_Orbits();
  Calculate_Population();
  if(method=='S')  Dipole_Acceleration_SaddlePoint();
  else if(method=='U')  Dipole_Acceleration_Uniform();
  else if(method=='C')  Dipole_Acceleration_Combined();
  else {
    std::cerr << "Error in method for calculating quantum orbit response, setting to uniform approximation.\n";
    Dipole_Acceleration_Uniform();
  }
  Fit_to_wAxis();
  //for(int nw=0; nw<pT->Nt; nw++)  DA[nw]*=exp(-i*pT->w[nw]*pLP[colourRef]->t0);
}

// As above but for default Et
void qo_response_1d::Dipole_Acceleration() {
  Dipole_Acceleration(pLaser->Et);
}

// As above but for default E_mpi / Ert (depending on flag) with specified row
void qo_response_1d::Dipole_Acceleration(int nr, bool mpi=false) {
  (mpi==true) ? Dipole_Acceleration(&pLaser->E_mpi[Index2V(nr,0,pT->Nt)]) : Dipole_Acceleration(&pLaser->Ert[Index2V(nr,0,pT->Nt)]);
}

// As above but for a specifed row of temporo-spatial field
void qo_response_1d::Dipole_Acceleration(dcmplx *Ert, int nr) {
  Dipole_Acceleration(&Ert[Index2V(nr,0,pT->Nt)]);
}  

// Dipole acceleration in freq-domain from specified temporal fields
void qo_response_1d::Dipole_Acceleration(dcmplx *EtX, dcmplx *EtY) {
  Dipole_Acceleration(EtX);  // just to give form to this virtual function in dipole_response class
}

// Dipole acceleration in freq-domain importing row number of specified temporospatial fields
void qo_response_1d::Dipole_Acceleration(dcmplx *ErtX, dcmplx *ErtY, int nr) {
  Dipole_Acceleration(ErtX, nr);  // just to give form to this virtual function in dipole_response class
}  

/***************************************************************************************************************/
// Class to create two-dimensional dipole response via QO method: (assumes exactly 2 colours) 
class qo_response_2d : public qo_response_1d {  
 private:    
  /***************************************************************************************************************/
  // Private Class Members:
  dcmplxvec spectrum_shortY, spectrum_longY, spectrum_pairY;
  dcmplxvec DA_localY; 

  /***************************************************************************************************************/
  // Private Class Function Prototypes:
  void Initialise();
  void Build_Energy_Axis();
  dcmplx Ex(dcmplx);
  dcmplx Ax(dcmplx);
  dcmplx intAx(dcmplx);
  dcmplx intA2x(dcmplx);  
  dcmplx Ey(dcmplx);
  dcmplx Ay(dcmplx);
  dcmplx intAy(dcmplx);
  dcmplx intA2y(dcmplx);
  dcmplx px(dcmplx, dcmplx);
  dcmplx py(dcmplx, dcmplx);
  dcmplx S(dcmplx, dcmplx, double);
  dcmplx dtdS(dcmplx, dcmplx, double);
  dcmplx dtpdS(dcmplx, dcmplx);
  dcmplx dttd2S(dcmplx, dcmplx);
  dcmplx dtptpd2S(dcmplx, dcmplx);
  dcmplx dtptd2S(dcmplx, dcmplx);
  dcmplx Determinant_Hessian(dcmplx, dcmplx);
  bool Find_Initial_Saddle(dcmplx &, dcmplx &);
  bool Find_Saddle(dcmplx &, dcmplx &, double, bool);
  double Calculate_Increments(dcmplx &, dcmplx &, dcmplx, dcmplx, double);
  bool Check_Convergence(dcmplx, dcmplx, double, bool &, bool);
  bool Add_Saddle(dcmplx, dcmplx);
  void Sort_Saddle_Time();
  void Saddlepoint_Search();
  bool Track_Orbits();
  void Calculate_Population();
  void DTME(dcmplx &, dcmplx &, dcmplx, dcmplx, dcmplx, double);
  void DTME_singular(dcmplx &, dcmplx &, dcmplx, dcmplx, dcmplx, double);
  void Dipole_Acceleration_SaddlePoint();
  void Dipole_Acceleration_Uniform();
  void Dipole_Acceleration_Combined();
  void Fit_to_wAxis();

 public:
  /***************************************************************************************************************/
  // Public Class Function Prototypes:
  //qo_response_2d(laser_parameters);
  qo_response_2d(temporal_axes &, electric_field &, gas_jet &, bool);
  qo_response_2d(temporal_axes &, electric_field &, gas_jet &, std::ifstream &);
  void Dipole_Acceleration(dcmplx *, dcmplx *);
  void Dipole_Acceleration(); 
  void Dipole_Acceleration(int, bool);
  void Dipole_Acceleration(dcmplx *, dcmplx *, int); 
};

/***************************************************************************************************************/
// Class Functions:
 
/* // Constructor using imported laser parameter specification - DO LATER */
/* qo_response_2d::qo_response_2d(laser_parameters &LP_) : qo_response_1d(laser_parameters &LP_) {  */
/*   Initialise(); */
/* }  */

// Constructor diagnosing laser parameters from imported electric field (creates new laser_parameters object then copies this to pLP)
qo_response_2d::qo_response_2d(temporal_axes &T_, electric_field &Laser_, gas_jet &Jet_, bool findAll_=false) : qo_response_1d(T_, Laser_, Jet_, findAll_) {
  Initialise();
  // DTME_Prefactor();  // put this within Dipole_Acceleration() - could feasibly change if >1 jet
}

// As above but takes in control parameters from specified ifstream
qo_response_2d::qo_response_2d(temporal_axes &T_, electric_field &Laser_, gas_jet &Jet_, std::ifstream &INSTREAM) : qo_response_1d(T_, Laser_, Jet_, INSTREAM) {
  Initialise();
  if(colourRef>1)  colourRef = 0;  // resets erroneous colourRef specification (2 colours max for qo_response_2d; not needed in above constructor)
  // DTME_Prefactor();  // put this within Dipole_Acceleration() - could feasibly change if >1 jet
}

// Initialise memory
void qo_response_2d::Initialise() {
  DA_localY.resize(NEn);
  //if(method=='C') {spectrum_shortY.resize(NEn);  spectrum_longY.resize(NEn);  spectrum_pairY.resize(NEn);}

  if(method=='C') {spectrum_shortY.resize(NEn);  spectrum_longY.resize(NEn);  spectrum_pairY.resize(NEn);}
  else if(method=='U')  spectrum_pairY.resize(NEn);

}

// Function to build the energy axis appropriate to these field values (& Ip)
void qo_response_2d::Build_Energy_Axis() {
  double En_lower = pJet->Ip+0.1*eV_to_au;//pJet->Ip;
  double E0_max = sqrt(pow(pLP[0]->E0,2)+pow(pLP[1]->E0,2));
  double w0_mix = (pLP[0]->w0*pLP[0]->E0+pLP[1]->w0*pLP[1]->E0)/E0_max;
  double En_upper = 3.0*(1.32*pJet->Ip+3.17*pow(0.5*E0_max/w0_mix,2));
  En[0] = En_upper;
  dEn = (En_upper-En_lower)/(double)(NEn-1);
  for(int nEn=1; nEn<NEn; nEn++)  En[nEn] = En[nEn-1]-dEn;
}

// Electric field along x-axis
dcmplx qo_response_2d::Ex(dcmplx t) {
  dcmplx Et = 0.0;
  if(real(t-pLP[0]->t0)>=-env_range[0] && real(t-pLP[0]->t0)<=env_range[0])
    Et = pLP[0]->E0/pLP[0]->w0 * cos((t-pLP[0]->t0)*pLP[0]->B) * (pLP[0]->w0*cos((t-pLP[0]->t0)*pLP[0]->B)*cos(pLP[0]->w0*(t-pLP[0]->t0)+pLP[0]->cep)-2.0*pLP[0]->B*sin((t-pLP[0]->t0)*pLP[0]->B)*sin(pLP[0]->w0*(t-pLP[0]->t0)+pLP[0]->cep));
  return Et;
}

// Vector Potential along x-axis
dcmplx qo_response_2d::Ax(dcmplx t) {
  dcmplx At = 0.0;
  if(real(t-pLP[0]->t0)>=-env_range[0] && real(t-pLP[0]->t0)<=env_range[0])
    At = - pLP[0]->E0/pLP[0]->w0 * pow(cos((t-pLP[0]->t0)*pLP[0]->B), 2) * sin(pLP[0]->w0*(t-pLP[0]->t0)+pLP[0]->cep);
  return At;
}

// Integral of Vector Potential along x-axis
dcmplx qo_response_2d::intAx(dcmplx t) {
  dcmplx intAt = 0.0;
  if(real(t-pLP[0]->t0)>=-env_range[0] && real(t-pLP[0]->t0)<=env_range[0])
      intAt = 0.5*pLP[0]->E0/pLP[0]->w0 * (cos(pLP[0]->w0*(t-pLP[0]->t0)+pLP[0]->cep)/pLP[0]->w0+1.0/(pLP[0]->w0*pLP[0]->w0-4.0*pLP[0]->B*pLP[0]->B)*(pLP[0]->w0*cos(2.0*pLP[0]->B*(t-pLP[0]->t0))*cos(pLP[0]->w0*(t-pLP[0]->t0)+pLP[0]->cep)+2.0*pLP[0]->B*sin(2.0*pLP[0]->B*(t-pLP[0]->t0))*sin(pLP[0]->w0*(t-pLP[0]->t0)+pLP[0]->cep)));
  return intAt;
}

// Integral of (Vector Potential Squared) along x-axis
dcmplx qo_response_2d::intA2x(dcmplx t) {
  dcmplx intAt2 = 0.0;
  if(real(t-pLP[0]->t0)>=-env_range[0] && real(t-pLP[0]->t0)<=env_range[0])
      intAt2 = 0.5*pow(0.5*pLP[0]->E0/pLP[0]->w0, 2) * (1.5*(t-pLP[0]->t0)+sin(2.0*pLP[0]->B*(t-pLP[0]->t0))/pLP[0]->B+0.125*sin(4.0*pLP[0]->B*(t-pLP[0]->t0))/pLP[0]->B-0.25*(3.0/pLP[0]->w0+4.0*pLP[0]->w0/(pLP[0]->w0*pLP[0]->w0-pLP[0]->B*pLP[0]->B)*cos(2.0*pLP[0]->B*(t-pLP[0]->t0))+pLP[0]->w0/(pLP[0]->w0*pLP[0]->w0-4.0*pLP[0]->B*pLP[0]->B)*cos(4.0*pLP[0]->B*(t-pLP[0]->t0)))*sin(2.0*(pLP[0]->w0*(t-pLP[0]->t0)+pLP[0]->cep))+0.5*pLP[0]->B*(2.0/(pLP[0]->w0*pLP[0]->w0-pLP[0]->B*pLP[0]->B)*sin(2.0*pLP[0]->B*(t-pLP[0]->t0))+1.0/(pLP[0]->w0*pLP[0]->w0-4.0*pLP[0]->B*pLP[0]->B)*sin(4.0*pLP[0]->B*(t-pLP[0]->t0)))*cos(2.0*(pLP[0]->w0*(t-pLP[0]->t0)+pLP[0]->cep))) - intA2s[0];
  else if(real(t-pLP[0]->t0)>env_range[0])  intAt2 = intA2f[0];
  return intAt2;
}

// Electric field along y-axis
dcmplx qo_response_2d::Ey(dcmplx t) {
  dcmplx Et = 0.0;
  if(real(t-pLP[1]->t0)>=-env_range[1] && real(t-pLP[1]->t0)<=env_range[1])
    Et = pLP[1]->E0/pLP[1]->w0 * cos((t-pLP[1]->t0)*pLP[1]->B) * (pLP[1]->w0*cos((t-pLP[1]->t0)*pLP[1]->B)*cos(pLP[1]->w0*(t-pLP[1]->t0)+pLP[1]->cep)-2.0*pLP[1]->B*sin((t-pLP[1]->t0)*pLP[1]->B)*sin(pLP[1]->w0*(t-pLP[1]->t0)+pLP[1]->cep));
  return Et;
}

// Vector Potential along x-axis
dcmplx qo_response_2d::Ay(dcmplx t) {
  dcmplx At = 0.0;
  if(real(t-pLP[1]->t0)>=-env_range[1] && real(t-pLP[1]->t0)<=env_range[1])
    At = - pLP[1]->E0/pLP[1]->w0 * pow(cos((t-pLP[1]->t0)*pLP[1]->B), 2) * sin(pLP[1]->w0*(t-pLP[1]->t0)+pLP[1]->cep);
  return At;
}

// Integral of Vector Potential along x-axis
dcmplx qo_response_2d::intAy(dcmplx t) {
  dcmplx intAt = 0.0;
  if(real(t-pLP[1]->t0)>=-env_range[1] && real(t-pLP[1]->t0)<=env_range[1])
      intAt = 0.5*pLP[1]->E0/pLP[1]->w0 * (cos(pLP[1]->w0*(t-pLP[1]->t0)+pLP[1]->cep)/pLP[1]->w0+1.0/(pLP[1]->w0*pLP[1]->w0-4.0*pLP[1]->B*pLP[1]->B)*(pLP[1]->w0*cos(2.0*pLP[1]->B*(t-pLP[1]->t0))*cos(pLP[1]->w0*(t-pLP[1]->t0)+pLP[1]->cep)+2.0*pLP[1]->B*sin(2.0*pLP[1]->B*(t-pLP[1]->t0))*sin(pLP[1]->w0*(t-pLP[1]->t0)+pLP[1]->cep)));
  return intAt;
}

// Integral of (Vector Potential Squared) along x-axis
dcmplx qo_response_2d::intA2y(dcmplx t) {
  dcmplx intAt2 = 0.0;
  if(real(t-pLP[1]->t0)>=-env_range[1] && real(t-pLP[1]->t0)<=env_range[1])
      intAt2 = 0.5*pow(0.5*pLP[1]->E0/pLP[1]->w0, 2) * (1.5*(t-pLP[1]->t0)+sin(2.0*pLP[1]->B*(t-pLP[1]->t0))/pLP[1]->B+0.125*sin(4.0*pLP[1]->B*(t-pLP[1]->t0))/pLP[1]->B-0.25*(3.0/pLP[1]->w0+4.0*pLP[1]->w0/(pLP[1]->w0*pLP[1]->w0-pLP[1]->B*pLP[1]->B)*cos(2.0*pLP[1]->B*(t-pLP[1]->t0))+pLP[1]->w0/(pLP[1]->w0*pLP[1]->w0-4.0*pLP[1]->B*pLP[1]->B)*cos(4.0*pLP[1]->B*(t-pLP[1]->t0)))*sin(2.0*(pLP[1]->w0*(t-pLP[1]->t0)+pLP[1]->cep))+0.5*pLP[1]->B*(2.0/(pLP[1]->w0*pLP[1]->w0-pLP[1]->B*pLP[1]->B)*sin(2.0*pLP[1]->B*(t-pLP[1]->t0))+1.0/(pLP[1]->w0*pLP[1]->w0-4.0*pLP[1]->B*pLP[1]->B)*sin(4.0*pLP[1]->B*(t-pLP[1]->t0)))*cos(2.0*(pLP[1]->w0*(t-pLP[1]->t0)+pLP[1]->cep))) - intA2s[1];
  else if(real(t-pLP[1]->t0)>env_range[1])  intAt2 = intA2f[1];
  return intAt2;
}

// Complex Saddle-point Momentum along x-axis
dcmplx qo_response_2d::px(dcmplx t, dcmplx tp) {
  return (intAx(tp)-intAx(t))/(t-tp);
}

// Complex Saddle-point Momentum along y-axis
dcmplx qo_response_2d::py(dcmplx t, dcmplx tp) {
  return (intAy(tp)-intAy(t))/(t-tp);
}

// Action
dcmplx qo_response_2d::S(dcmplx t, dcmplx tp, double Ef) {
  dcmplx psx = px(t,tp);
  dcmplx psy = py(t,tp);
  return -(pJet->Ip*(t-tp)-0.5*(psx*psx+psy*psy)*(t-tp)+0.5*(intA2x(t)-intA2x(tp)+intA2y(t)-intA2y(tp))-Ef*t);
}

// dS/dt
dcmplx qo_response_2d::dtdS(dcmplx t, dcmplx tp, double Ef) {
  dcmplx psx = px(t,tp);
  dcmplx psy = py(t,tp);
  return -(pJet->Ip-Ef+0.5*(pow(psx+Ax(t),2.0)+pow(psy+Ay(t),2.0)));
}

// dS/dtp
dcmplx qo_response_2d::dtpdS(dcmplx t, dcmplx tp) {
  dcmplx psx = px(t,tp);
  dcmplx psy = py(t,tp);
  return -(-pJet->Ip-0.5*(pow(psx+Ax(tp),2.0)+pow(psy+Ay(tp),2.0)));
}

// d2S/dtdt
dcmplx qo_response_2d::dttd2S(dcmplx t, dcmplx tp) {
  dcmplx psx = px(t,tp);
  dcmplx psy = py(t,tp);
  return -(-(psx+Ax(t))*(Ex(t)+(psx+Ax(t))/(t-tp))-(psy+Ay(t))*(Ey(t)+(psy+Ay(t))/(t-tp)));
}

// d2S/dtpdtp
dcmplx qo_response_2d::dtptpd2S(dcmplx t, dcmplx tp) {
  dcmplx psx = px(t,tp);
  dcmplx psy = py(t,tp);
  return -((psx+Ax(tp))*(Ex(tp)-(psx+Ax(tp))/(t-tp))+(psy+Ay(tp))*(Ey(tp)-(psy+Ay(tp))/(t-tp)));
}

// d2S/dtpdt
dcmplx qo_response_2d::dtptd2S(dcmplx t, dcmplx tp) {
  dcmplx psx = px(t,tp);
  dcmplx psy = py(t,tp);
  return -((psx+Ax(t))*(psx+Ax(tp))/(t-tp)+(psy+Ay(t))*(psy+Ay(tp))/(t-tp));
}

// Determinant of Hessian Matrix
dcmplx qo_response_2d::Determinant_Hessian(dcmplx t, dcmplx tp) {
  dcmplx DtptD2S = dtptd2S(t,tp);
  return dttd2S(t,tp)*dtptpd2S(t,tp)-DtptD2S*DtptD2S;
}

// Find saddle points at En[0]
bool qo_response_2d::Find_Initial_Saddle(dcmplx &t, dcmplx &tp) {
  dcmplx t_guess(0.0), tp_guess(0.0), tau_guess(0.0);
  bool found = false;
  real(tp_guess) = fmod(rand(), max_tp-min_tp) + min_tp;
  real(tau_guess) = fmod(rand(), max_tau-min_tau) + min_tau;
  real(t_guess) = real(tp_guess) + real(tau_guess);
  imag(tp_guess) = fmod(rand(), 0.7*pLP[colourRef]->T);
  imag(t_guess) = fmod(rand(), 1.4*pLP[colourRef]->T) - 0.7*pLP[colourRef]->T;
  found = Find_Saddle(t_guess, tp_guess, En[0], true);
  if(found==true) {t = t_guess;  tp = tp_guess;}
  return found;
}

// Find single saddle point at a specifed energy
bool qo_response_2d::Find_Saddle(dcmplx &t, dcmplx &tp, double Energy, bool boundcheck=true) {
  dcmplx dt(0.0), dtp(0.0);
  double error = 1.0e10;//0.0;
  bool found = false, stop = false;
  for(int nit=0; nit<Nit; nit++) {
    error = Calculate_Increments(dt, dtp, t, tp, Energy);
    found = Check_Convergence(t, tp, error, stop, boundcheck);
    if(found==true)  return true;
    else if(stop==true)  return false;
    // t += dt;  tp += dtp;
    if(nit!=Nit-1) {t += dt;  tp += dtp;}
  }
  //cout << stop << '\t' << found << '\t' << error << '\n';
  return found;
}

// Function to calculate increments to t,t' and return a measure of convergene in 1st derivatives of action
double qo_response_2d::Calculate_Increments(dcmplx &dt, dcmplx &dtp, dcmplx t, dcmplx tp, double Energy) {
  dcmplx DtDS = dtdS(t,tp,Energy);
  dcmplx DtpDS = dtpdS(t,tp);
  dcmplx DttD2S = dttd2S(t,tp);
  dcmplx DtptpD2S = dtptpd2S(t,tp);
  dcmplx DtptD2S = dtptd2S(t,tp);
  dcmplx DetHess = Determinant_Hessian(t,tp);
  dt = (DtpDS*DtptD2S-DtDS*DtptpD2S)/DetHess;
  dtp = (DtDS*DtptD2S-DtpDS*DttD2S)/DetHess;
  double error = abs(DtDS) + abs(DtpDS);
  //double error = abs(dt) + abs(dtp);
  return error;
}

// Function to check convergence to saddle-point solution
bool qo_response_2d::Check_Convergence(dcmplx t, dcmplx tp, double error, bool &stop, bool boundcheck) {
  dcmplx tau = t-tp;
  if(error<1.0e-8) {
    if(boundcheck==true) {
      if((real(tau)>min_tau) && (real(tau)<max_tau)) {
	if((real(tp)>=min_tp) && (real(tp)<=max_tp)) {stop = false;  return true;}
	else {stop = true;  return false;}
      }
      else {stop = true;  return false;}
    }
    else {stop = false;  return true;}
  }
  else {stop = false;  return false;}
}

// Add found saddle to vector if not already present and satifying  Im(t')>=0
bool qo_response_2d::Add_Saddle(dcmplx t, dcmplx tp) {
  bool unique = false;
  if(imag(tp)>=0.0) {
    unique = true;
    for(int ns=0; ns<t_saddle.size(); ns++)
      //if((abs(t_saddle[ns]-t)<1.0e-5) && (abs(tp_saddle[ns]-tp)<1.0e-5))  unique = false;
      if((abs(t_saddle[ns]-t)<1.0e-5*pLP[colourRef]->T) && (abs(tp_saddle[ns]-tp)<1.0e-5*pLP[colourRef]->T))  unique = false;
    if(unique==true) {
      t_saddle.push_back(t);
      tp_saddle.push_back(tp);
    }
  }
  return unique;
}

// Sort final element into the vector according to ascending t
void qo_response_2d::Sort_Saddle_Time() {
  int ns = 0;
  dcmplx temp;
  while(ns<t_saddle.size() && (real(t_saddle.back())>real(t_saddle[ns])))  ns++;
  if(ns<t_saddle.size()-1) {
    temp = t_saddle.back();
    t_saddle.insert(t_saddle.begin()+ns, temp);
    t_saddle.pop_back();
    temp = tp_saddle.back();
    tp_saddle.insert(tp_saddle.begin()+ns, temp);
    tp_saddle.pop_back();
  }
}

// Find all saddle points within specified region (and up to Nguess)
void qo_response_2d::Saddlepoint_Search() {
  t_saddle.clear();  tp_saddle.clear();
  dcmplx t(0,0), tp(0,0);
  bool found = false;
  for(int nguess=0; nguess<Nguess; nguess++) {
    found = Find_Initial_Saddle(t, tp);
    if(found==true) {
      found = Add_Saddle(t, tp);
      if(found==true) {
	Sort_Saddle_Time();
	//std::cout << "Sorted saddle tp=(" << real(tp)/pLP[colourRef]->T << ',' << imag(tp)/pLP[colourRef]->T << "); t=(" << real(t)/pLP[colourRef]->T << ',' << imag(t)/pLP[colourRef]->T << "), spread = " << arg(1.0e-10+i*(t-tp))/pi << " pi rad\n";
      }
      found = false;
    }
  }
  //std::cout << "Found " << t_saddle.size() << " unique saddle-points and sorted according to ascending re(t)\n";
}

// Function to follow saddle-point orbits with decreasing energy
bool qo_response_2d::Track_Orbits() {
  t_orbits.clear();  tp_orbits.clear();
  dcmplx t, tp;
  dcmplx t_prev, tp_prev;
  dcmplx dt, dtp;
  bool found;
  int nEn = 1;
  double Energy = En[0];
  dcmplx S_track[3];
  double arg_S = 0.0;
  const int Nstep = 10, Nlevel = 10;//5;
  const double step_factor = 1.0/(double)(Nstep);
  int step_level, step_count[Nlevel], old_level;
  t_orbits.clear();  tp_orbits.clear();
  t_orbits.resize(t_saddle.size(), t_onepath);  tp_orbits.resize(tp_saddle.size(), tp_onepath);
  for(int ns=0; ns<t_saddle.size(); ns++) {
    nEn = 1;  Energy = En[0];
    t_prev = (t_onepath[0] = t_saddle[ns]);
    tp_prev = (tp_onepath[0] = tp_saddle[ns]);
    dt = 0.0;  dtp = 0.0;
    S_track[0] = 0.0;  S_track[1] = 0.0;  S_track[2] = 0.0;
    arg_S = 0.0;
    step_level = 0;
    for(int nlevel=0; nlevel<Nlevel; nlevel++)  step_count[nlevel] = 0;
    while(nEn<NEn) {

/*       Energy -= dEn; // this seems to work better....? */
/*       t = t_prev+dt; */
/*       tp = tp_prev+dtp; */

      Energy -= dEn*pow(step_factor, step_level);
      t = t_prev+dt*pow(step_factor, step_level);
      tp = tp_prev+dtp*pow(step_factor, step_level);

      found = Find_Saddle(t, tp, Energy, false);

/*       if(found==true) { */
/* 	S_track[2] = S(t, tp, Energy); */
/* 	if(nEn==1)  S_track[0] = S_track[2]; */
/* 	else if(nEn==2)  S_track[1] = S_track[2]; */
/* 	else { */
/* 	  arg_S = abs(arg(S_track[2]-S_track[1])-arg(S_track[1]-S_track[0])); */
/* 	  S_track[0] = S_track[1];  S_track[1] = S_track[2]; */
/* 	  if(arg_S>pi)  arg_S = 2.0*pi-arg_S; */
/* 	} */
/*       } */

      if(found==true) {
	S_track[2] = S(t, tp, Energy);
	if(nEn==1)  S_track[0] = S_track[2];
	else if(nEn==2)  S_track[1] = S_track[2];
	else {
	  if(fabs(real(S_track[2])-real(S_track[1]))>3.0*fabs(real(S_track[1])-real(S_track[0])))  found = false;
	  else {
	    arg_S = abs(arg(S_track[2]-S_track[1])-arg(S_track[1]-S_track[0]));
	    S_track[0] = S_track[1];  S_track[1] = S_track[2];
	    if(arg_S>pi)  arg_S = 2.0*pi-arg_S;
	  }
	}
      }

      //if(found==true && arg_S<0.5*pi) {
      if(found==true) {
	if(step_level>0)  step_count[step_level]++;
	old_level = step_level;
	for(int nlevel=step_level; nlevel>0; nlevel--) {
	  if(step_count[nlevel]==Nstep) {
	    step_level--;  // if reached a more major point, raise back 1 level
	    //std::cout << "Rejoined level " << step_level << " after " << step_count[step_level+1] << " steps.\n";
	    step_count[step_level]++;  // and advance that level one, recheck
	  }
	  else  break;
	}
	if(step_level==0) {
	  t_onepath[nEn] = t;
	  tp_onepath[nEn] = tp;
	  nEn++;
	}
	dt = (t-t_prev)*pow(step_factor, step_level-old_level);
	dtp = (tp-tp_prev)*pow(step_factor, step_level-old_level);
	t_prev = t;  tp_prev = tp;
      }
      else {
	//std::cout << "Reverting back from " << En*27.212;
	Energy += dEn*pow(step_factor, step_level);  // revert back
	//std::cout << "eV to " << En*27.212 << "eV\n";
	step_level++;  step_count[step_level] = 0;  // drop down level & reset step_count at new level

	if(step_level==Nlevel) {std::cerr << "Nlevel reached (ns = " << ns << ", nEn = " << nEn << ") - total failure!\n";  return false;}
	//if(step_level==Nlevel) {std::cout << "Nlevel reached (ns = " << ns << ", nEn = " << nEn << ") - total failure!\n";  continue;}

	dt *= step_factor;  dtp *= step_factor;  // reduce dt & dtp according to dropped level
	//std::cout << "Track_Orbit failed for ns = " << ns << " on nEn = " << nEn << ',' << step_level-1 << ',' << step_count[step_level-1] << " (En = " << En*27.212 << "eV vs En_cutoff = " << En_cutoff*27.212 << "eV) ... descending to level " << step_level << " of " << Nlevel << ".\n";
      }
    }
    t_orbits[ns] = t_onepath;
    tp_orbits[ns] = tp_onepath;
  }
  return true;
}

// Calculates population vector along real time axis
void qo_response_2d::Calculate_Population() {
  const double min_t = std::min(real(tp_orbits[0][NEn-1]),real(tp_orbits[1][NEn-1]));
  const double max_t = std::max(real(t_orbits[t_orbits.size()-2][NEn-1]),real(t_orbits[t_orbits.size()-1][NEn-1]));
  const int Nt = 2*NEn; //1024;
  const double dt = (max_t-min_t)/Nt;
  double t;
  population.clear();  population.resize(Nt);
  double Eabs, ADK;
  double exp_factor = 0.0;
  for(int nt=0; nt<Nt; nt++) {
    t = min_t+nt*dt;
    exp_factor -= dt*pJet->ADK_Rate(Ex(t), Ey(t));
    population[nt] = exp(exp_factor);
  }
}

// Dipole transition matrix element
void qo_response_2d::DTME(dcmplx &dtmeX, dcmplx &dtmeY, dcmplx kX, dcmplx kY, dcmplx pre_dtme, double phiE) {
  double kNorm = norm(kX)+norm(kY);
  if(pJet->n==1) {
    dtmeX = pre_dtme*kX/pow(kNorm+2.0*pJet->Ip, 3);
    dtmeY = pre_dtme*kY/pow(kNorm+2.0*pJet->Ip, 3);
  }
  else if(pJet->n==2)
    if(pJet->l==0) {
      dtmeX = pre_dtme*kX*(kNorm-4.0*pJet->Ip)/pow(kNorm+2.0*pJet->Ip, 4);
      dtmeY = pre_dtme*kY*(kNorm-4.0*pJet->Ip)/pow(kNorm+2.0*pJet->Ip, 4);
    }
    else if(pJet->l==1) {
      double phiP = atan2(real(kY),real(kX));
      dcmplx factorX = pre_dtme*((2.0*pJet->Ip-5.0*kNorm)*pow(cos(phiP),2)+kNorm*(kNorm+2.0*pJet->Ip)*pow(sin(phiP),2))/pow(kNorm+2.0*pJet->Ip, 4);
      dcmplx factorY = pre_dtme*((2.0*pJet->Ip-kNorm*(5.0+kNorm+2.0*pJet->Ip))*0.5*sin(2*(phiP)))/pow(kNorm+2.0*pJet->Ip, 4);
      dtmeX = (cos(phiE)*factorX-sin(phiE)*factorY)*kX;
      dtmeY = (sin(phiE)*factorX+cos(phiE)*factorY)*kY;
    }
  else {
    dtmeX = pre_dtme*0.0;
    dtmeY = pre_dtme*0.0;
  }
}

// Dipole transition matrix element
void qo_response_2d::DTME_singular(dcmplx &dtmeX, dcmplx &dtmeY, dcmplx kX, dcmplx kY, dcmplx pre_dtme, double phiE) {
  double kNorm = norm(kX)+norm(kY);
  if(pJet->n==1) {
    dtmeX = pre_dtme*kX/pow(2.0*pJet->Ip, 3);
    dtmeY = pre_dtme*kY/pow(2.0*pJet->Ip, 3);
  }
  else if(pJet->n==2)
    if(pJet->l==0) {
      dtmeX = pre_dtme*kX*(kNorm-4.0*pJet->Ip)/pow(2.0*pJet->Ip, 4);
      dtmeY = pre_dtme*kY*(kNorm-4.0*pJet->Ip)/pow(2.0*pJet->Ip, 4);
    }
    else if(pJet->l==1) {
      double phiP = atan2(real(kY),real(kX));
      dcmplx factorX = pre_dtme*((2.0*pJet->Ip-5.0*kNorm)*pow(cos(phiP),2)+kNorm*(kNorm+2.0*pJet->Ip)*pow(sin(phiP),2));
      dcmplx factorY = pre_dtme*((2.0*pJet->Ip-kNorm*(5.0+kNorm+2.0*pJet->Ip))*0.5*sin(2*(phiP)));
      dtmeX = (cos(phiE)*factorX-sin(phiE)*factorY)*kX;
      dtmeY = (sin(phiE)*factorX+cos(phiE)*factorY)*kY;
    }
  else {
    dtmeX = pre_dtme*0.0;
    dtmeY = pre_dtme*0.0;
  }
}

// Calculates harmonic spectrum resulting from previously specified orbits using the Saddle-Point Approx
void qo_response_2d::Dipole_Acceleration_SaddlePoint() {
  bool found_pair;
  int ns_pair;
  dcmplx t, tp;
  double epsilon = 1.0e-10;
  double Energy;
  intvec ns_store;
  dcmplx DetHess, DetHess_prev[2];
  dcmplx spread, spread_prev[2];
  double sign_spreadHess[2];
  dcmplx spreadHess;
  double min_t = std::min(real(tp_orbits[0][NEn-1]),real(tp_orbits[1][NEn-1]));
  double max_t = std::max(real(t_orbits[t_orbits.size()-2][NEn-1]),real(t_orbits[t_orbits.size()-1][NEn-1]));
  double dt = (max_t-min_t)/population.size();
  double pop_t, pop_tp;
  dcmplx cmplxAmp;
  dcmplx action;
  for(int nEn=0; nEn<t_orbits[0].size(); nEn++)  {DA_local[nEn] = 0.0;  DA_localY[nEn] = 0.0;}
  dcmplx pre_dtme = qo_response_1d::DTME_Prefactor();
  dcmplx dtmex, dtmey, dtmeSx, dtmeSy;
  for(int ns=0; ns<t_orbits.size(); ns++) {
    if(traj=='A' || ((traj=='S'||traj=='B'||traj=='C') && orbit_labels[ns][4]==-1) || ((traj=='L'||traj=='B') && orbit_labels[ns][4]==1)) {
      ns_store.push_back(ns);  // records position in orbit_labels
      sign_spreadHess[0] = 1.0;  sign_spreadHess[1] = 1.0;
      for(int nEn=0; nEn<t_orbits[0].size(); nEn++) {
	Energy = En[nEn];
	t = t_orbits[ns][nEn];
	tp = tp_orbits[ns][nEn];
	DetHess = Determinant_Hessian(t,tp);
	spread = epsilon+i*(t-tp);
	if(nEn>0) {
	  if(abs(arg(DetHess)-arg(DetHess_prev[0]))>pi)  sign_spreadHess[0] *= -1.0;
	  if(abs(arg(spread)-arg(spread_prev[0]))>pi)  sign_spreadHess[0] *= -1.0;
	}
	DetHess_prev[0] = DetHess;  spread_prev[0] = spread;
	spreadHess = sign_spreadHess[0]/std::polar(sqrt(pow(abs(spread),3)*abs(DetHess)), 0.5*(3.0*arg(spread)+arg(DetHess)));
	//spreadHess = sign_spreadHess[0]/pow(spread,1.5)/sqrt(DetHess);
	qo_response_1d::Fit_Population(pop_t, pop_tp, t, tp, min_t, dt);
	DTME(dtmex, dtmey, px(t,tp)+Ax(tp), py(t,tp)+Ay(tp), pre_dtme, atan2(real(Ey(tp)),real(Ex(tp))));
	DTME_singular(dtmeSx, dtmeSy, px(t,tp)+Ax(tp), py(t,tp)+Ay(tp), pre_dtme, atan2(real(Ey(tp)),real(Ex(tp))));
	cmplxAmp = i*2.0*pow(pi*Energy,2)*spreadHess*pop_t*pop_tp*(Ex(tp)*dtmeSx+Ey(tp)*dtmeSy);
 	action = S(t, tp, Energy);
 	DA_local[nEn] += conj(cmplxAmp*exp(i*action)*conj(dtmex));  // conj() as was previously giving time-reversed emission
 	DA_localY[nEn] += conj(cmplxAmp*exp(i*action)*conj(dtmey));
      }
    }
  }
}

// Calculates harmonic spectrum resulting from previously specified orbits using the Uniform Approx
void qo_response_2d::Dipole_Acceleration_Uniform() {
  bool found_pair;
  int ns, ns_pair;
  dcmplx t, tp;
  double epsilon = 1.0e-10;
  dcmplx DetHess, DetHess_prev[2];
  dcmplx spread, spread_prev[2];
  double sign_spreadHess[2];
  dcmplx spreadHess[2];
  double min_t = std::min(real(tp_orbits[0][NEn-1]),real(tp_orbits[1][NEn-1]));
  double max_t = std::max(real(t_orbits[t_orbits.size()-2][NEn-1]),real(t_orbits[t_orbits.size()-1][NEn-1]));
  double dt = (max_t-min_t)/population.size();
  double pop_t, pop_tp;
  dcmplx cmplxAmp[2];
  bool antistokes_passed;
  double Energy;
  dcmplx action[2];
  dcmplx cmplxAmp_plus, cmplxAmp_minus;
  double short_shift=0.0, long_shift=0.0;
  double diff_im_S, prev_diff_im_S;
  dcmplx action_plus, action_minus;
  double sign_spectrum, eta;
  dcmplx sqrt_z, minus_z;
  dcmplx Airy, dzAiry;
  int ierr = 0, argflag = 0, modflag = 0;
  dcmplx pre_factor = sqrt(6.0*pi)*exp(i*0.25*pi);
  for(int nEn=0; nEn<t_orbits[0].size(); nEn++)  {DA_local[nEn] = 0.0;  DA_localY[nEn] = 0.0;}
  dcmplx pre_dtme = qo_response_1d::DTME_Prefactor();
  dcmplx dtmex[2], dtmey[2], dtmeSx, dtmeSy;
  double E0_max = sqrt(pow(pLP[0]->E0,2)+pow(pLP[1]->E0,2));
  for(ns=0; ns<t_orbits.size(); ns++) {
    if(orbit_labels[ns][4]==-1 && orbit_labels[ns][6]==1) {
      found_pair = false;
      for(ns_pair=0; ns_pair<t_orbits.size(); ns_pair++)
	if(orbit_labels[ns_pair][0]==orbit_labels[ns][0] && orbit_labels[ns_pair][1]==orbit_labels[ns][1] && orbit_labels[ns_pair][2]==orbit_labels[ns][2] && orbit_labels[ns_pair][3]==orbit_labels[ns][3] && orbit_labels[ns_pair][4]==1 && orbit_labels[ns_pair][6]==1) {
	  found_pair = true;
	  break;
	}
      if(found_pair==false) {
	//std::cout << "Failed to find pair for orbit " << orbit_labels[ns][0] << ',' << orbit_labels[ns][1] << ',' << orbit_labels[ns][2] << ',' << orbit_labels[ns][3] << ',' << orbit_labels[ns][4] << '\n';

	std::cout << "\nFailed to find pair for orbit " << orbit_labels[ns][0] << ',' << orbit_labels[ns][1] << ',' << orbit_labels[ns][2] << ',' << orbit_labels[ns][3] << ',' << orbit_labels[ns][4] << '\n';
	for(int nss=0; nss<t_orbits.size(); nss++)
	  std::cout << nss << ":\t" << orbit_labels[nss][0] << ',' << orbit_labels[nss][1] << ',' << orbit_labels[nss][2] << ',' << orbit_labels[nss][3] << ',' << orbit_labels[nss][4] << ',' << orbit_labels[nss][5] << ',' << orbit_labels[nss][6] << '\n';
	std::cout << '\n';
	return;
      }
      eta = (double)(abs(orbit_labels[ns][4]+orbit_labels[ns][5]))-1.0;
      //std::cout << "Orbit pair (" << orbit_labels[ns][0] << ',' << orbit_labels[ns][1] << ',' << orbit_labels[ns][2] << ")+(" << orbit_labels[ns_pair][0] << ',' << orbit_labels[ns_pair][1] << ',' << orbit_labels[ns_pair][2] << ") has eta = " << eta;// << '\n';
      sign_spreadHess[0] = 1.0;  sign_spreadHess[1] = 1.0;
      prev_diff_im_S = 0.0;
      antistokes_passed = false;
      for(int nEn=0; nEn<t_orbits[0].size(); nEn++) {
	Energy = En[nEn];
	// short orbit
	t = t_orbits[ns][nEn];
	tp = tp_orbits[ns][nEn];
 	DTME(dtmex[0], dtmey[0], px(t,tp)+Ax(tp), py(t,tp)+Ay(tp), pre_dtme, atan2(real(Ey(tp)),real(Ex(tp)))); 
	DTME_singular(dtmeSx, dtmeSy, px(t,tp)+Ax(tp), py(t,tp)+Ay(tp), pre_dtme, atan2(real(Ey(tp)),real(Ex(tp))));
	DetHess = Determinant_Hessian(t,tp);
	spread = epsilon+i*(t-tp);
	if(nEn>0) {
	  if(abs(arg(DetHess)-arg(DetHess_prev[0]))>pi)  sign_spreadHess[0] *= -1.0;
	  if(abs(arg(spread)-arg(spread_prev[0]))>pi)  sign_spreadHess[0] *= -1.0;
	}
	DetHess_prev[0] = DetHess;  spread_prev[0] = spread;
	spreadHess[0] = sign_spreadHess[0]/std::polar(sqrt(pow(abs(spread),3)*abs(DetHess)), 0.5*(3.0*arg(spread)+arg(DetHess)));
	//spreadHess[0] = sign_spreadHess[0]/pow(spread,1.5)/sqrt(DetHess);
	qo_response_1d::Fit_Population(pop_t, pop_tp, t, tp, min_t, dt);
	cmplxAmp[0] = i*2.0*pow(pi*Energy,2) * spreadHess[0]*pop_t*pop_tp*(Ex(tp)*dtmeSx+Ey(tp)*dtmeSy);
	action[0] = S(t, tp, Energy);
	// long orbit
	t = t_orbits[ns_pair][nEn];
	tp = tp_orbits[ns_pair][nEn];
	DTME(dtmex[1], dtmey[1], px(t,tp)+Ax(tp), py(t,tp)+Ay(tp), pre_dtme, atan2(real(Ey(tp)),real(Ex(tp))));
	DTME_singular(dtmeSx, dtmeSy, px(t,tp)+Ax(tp), py(t,tp)+Ay(tp), pre_dtme, atan2(real(Ey(tp)),real(Ex(tp))));
	DetHess = Determinant_Hessian(t,tp);
	spread = epsilon+i*(t-tp);
	if(nEn>0) {
	  if(abs(arg(DetHess)-arg(DetHess_prev[1]))>pi)  sign_spreadHess[1] *= -1.0;
	  if(abs(arg(spread)-arg(spread_prev[1]))>pi)  sign_spreadHess[1] *= -1.0;
	}
	DetHess_prev[1] = DetHess;  spread_prev[1] = spread;
	spreadHess[1] = sign_spreadHess[1]/std::polar(sqrt(pow(abs(spread),3)*abs(DetHess)), 0.5*(3.0*arg(spread)+arg(DetHess)));
	//spreadHess[1] = sign_spreadHess[1]/pow(spread,1.5)/sqrt(DetHess);
	qo_response_1d::Fit_Population(pop_t, pop_tp, t, tp, min_t, dt);
	cmplxAmp[1] = i*2.0*pow(pi*Energy,2) * spreadHess[1]*pop_t*pop_tp*(Ex(tp)*dtmeSx+Ey(tp)*dtmeSy);
	action[1] = S(t, tp, Energy);
	// adjusting orbits if necessary
	if(nEn==0) {
	  // potentially a problem here when both short and long orbit 'groups' form >1 band each
	  short_shift = 0.0;  long_shift = 0.0;  // reset values
/* 	  if(orbit_labels[ns][1]!=0 && orbit_labels[ns][3]!=0) */
/* 	    if(arg(spreadHess[1])-arg(spreadHess[0])>pi)  long_shift -= pi;  // if phase(long)>phase(short)+pi subract pi (but keep phase(long)>phase(short)) */
/* 	    else if(arg(spreadHess[1])-arg(spreadHess[0])<-pi) {short_shift -=pi;  long_shift += pi;}  // if phase(long)+pi<phase(short) subract pi from phase(short) and add pi to phase(long) */
/* 	    else if(arg(spreadHess[1])-arg(spreadHess[0])<0.0)  long_shift += pi;  // if instead phase(long)<phase(short) just add pi to phase(long) */
/* 	  if(arg(spreadHess[1])-arg(spreadHess[0])>pi)  long_shift -= pi;  // if phase(long)>phase(short)+pi subract pi (but keep phase(long)>phase(short)) */
/* 	  else if(arg(spreadHess[1])-arg(spreadHess[0])<-pi) {short_shift -=pi;  long_shift += pi;}  // if phase(long)+pi<phase(short) subract pi from phase(short) and add pi to phase(long) */
/* 	  else if(arg(spreadHess[1])-arg(spreadHess[0])<0.0)  long_shift += pi;  // if instead phase(long)<phase(short) just add pi to phase(long) */

	  if((pLP[0]->I0>1e10*Wcm2_to_au && orbit_labels[ns][1]>1) || (pLP[1]->I0>1e10*Wcm2_to_au && orbit_labels[ns][3]>1))
	    if(arg(spreadHess[1])-arg(spreadHess[0])>pi)  long_shift -= pi;  // if phase(long)>phase(short)+pi subract pi (but keep phase(long)>phase(short))
	    else if(arg(spreadHess[1])-arg(spreadHess[0])<-pi) {short_shift -=pi;  long_shift += pi;}  // if phase(long)+pi<phase(short) subract pi from phase(short) and add pi to phase(long)
	    else if(arg(spreadHess[1])-arg(spreadHess[0])<0.0)  long_shift += pi;  // if instead phase(long)<phase(short) just add pi to phase(long)
	}
	// uniform approx - note different E_tp (short) and E(tp) (long)
	diff_im_S = imag(action[0])-imag(action[1]);
	action_plus = 0.5 * (action[0]+action[1]);
	action_minus = 0.5 * (action[0]-action[1]);
	if(antistokes_passed==false && ((diff_im_S>0.0&&prev_diff_im_S<0.0) || (diff_im_S<0.0&&prev_diff_im_S>0.0))) {
	  antistokes_passed = true;
	  //std::cout << ", E_aS = " << En[nEn]*27.212 << "eV\n";
	}
 	if(antistokes_passed==false) {
	  sign_spectrum = -1.0;
	  sqrt_z = -pow(1.5*action_minus, 1.0/3.0)*exp(i*(double)eta*pi/3.0);
	  minus_z = -pow(sqrt_z, 2);
	}
	else {
	  sign_spectrum = 1.0;
	  sqrt_z = pow(1.5*action_minus, 1.0/3.0);
	  minus_z = -pow(sqrt_z, 2);
	}
 	__airy_functions_complex_double_airy_aic(&minus_z, &Airy, &dzAiry, &ierr, &argflag, &modflag);
	cmplxAmp_plus = 0.5 * (cmplxAmp[0]*conj(dtmex[0])*exp(i*short_shift)+i*cmplxAmp[1]*conj(dtmex[1])*exp(i*long_shift));
	cmplxAmp_minus = 0.5 * (cmplxAmp[0]*conj(dtmex[0])*exp(i*short_shift)-i*cmplxAmp[1]*conj(dtmex[1])*exp(i*long_shift));
	//DA_local[nEn] += conj(sign_spectrum*pre_factor*sqrt(action_minus)*exp(i*action_plus)*(cmplxAmp_minus/sqrt_z*Airy-i*cmplxAmp_plus/minus_z*dzAiry));  // conj() as was previously giving time-reversed emission
	spectrum_pair[nEn] = conj(sign_spectrum*pre_factor*sqrt(action_minus)*exp(i*action_plus)*(cmplxAmp_minus/sqrt_z*Airy-i*cmplxAmp_plus/minus_z*dzAiry));
/* 	if(real(cmplxAmp_plus)<=real(cmplxAmp_minus)) */
/* 	  spectrum_pair[nEn] = conj(sign_spectrum*pre_factor*sqrt(action_minus)*exp(i*action_plus)*(cmplxAmp_minus/sqrt_z*Airy-i*cmplxAmp_plus/minus_z*dzAiry)); */
/* 	else */
/* 	  spectrum_pair[nEn] = conj(sign_spectrum*pre_factor*sqrt(action_minus)*exp(i*action_plus)*(cmplxAmp_plus/sqrt_z*Airy-i*cmplxAmp_minus/minus_z*dzAiry)); */


	cmplxAmp_plus = 0.5 * (cmplxAmp[0]*conj(dtmey[0])*exp(i*short_shift)+i*cmplxAmp[1]*conj(dtmey[1])*exp(i*long_shift));
	cmplxAmp_minus = 0.5 * (cmplxAmp[0]*conj(dtmey[0])*exp(i*short_shift)-i*cmplxAmp[1]*conj(dtmey[1])*exp(i*long_shift));
	//DA_localY[nEn] += conj(sign_spectrum*pre_factor*sqrt(action_minus)*exp(i*action_plus)*(cmplxAmp_minus/sqrt_z*Airy-i*cmplxAmp_plus/minus_z*dzAiry));  // conj() as was previously giving time-reversed emission
	spectrum_pairY[nEn] = conj(sign_spectrum*pre_factor*sqrt(action_minus)*exp(i*action_plus)*(cmplxAmp_minus/sqrt_z*Airy-i*cmplxAmp_plus/minus_z*dzAiry));
/* 	if(real(cmplxAmp_plus)<=real(cmplxAmp_minus)) */
/* 	  spectrum_pairY[nEn] = conj(sign_spectrum*pre_factor*sqrt(action_minus)*exp(i*action_plus)*(cmplxAmp_minus/sqrt_z*Airy-i*cmplxAmp_plus/minus_z*dzAiry)); */
/* 	else  */
/* 	  spectrum_pairY[nEn] = conj(sign_spectrum*pre_factor*sqrt(action_minus)*exp(i*action_plus)*(cmplxAmp_plus/sqrt_z*Airy-i*cmplxAmp_minus/minus_z*dzAiry)); */

	prev_diff_im_S = diff_im_S;
      }

      //if(norm(spectrum_pair[NEn-1])<10.0*norm(spectrum_pair[NEn-31]))
      //if((norm(spectrum_pair[NEn-1])<10.0*norm(spectrum_pair[NEn-11])||norm(spectrum_pair[NEn-1])<10.0*norm(spectrum_pair[NEn-21])||norm(spectrum_pair[NEn-1])<10.0*norm(spectrum_pair[NEn-31])) && norm(spectrum_pair[0])<norm(spectrum_pair[NEn-10]))
      //      if(norm(spectrum_pair[NEn-1])<10.0*norm(spectrum_pair[NEn-11])||norm(spectrum_pair[NEn-1])<10.0*norm(spectrum_pair[NEn-21])||norm(spectrum_pair[NEn-1])<10.0*norm(spectrum_pair[NEn-31]))

      //      if((norm(spectrum_pair[NEn-1])<10.0*norm(spectrum_pair[NEn-11])||norm(spectrum_pair[NEn-1])<10.0*norm(spectrum_pair[NEn-21])||norm(spectrum_pair[NEn-1])<10.0*norm(spectrum_pair[NEn-31])) && norm(spectrum_pair[NEn-11])<1.0e-6)

/*       if((norm(spectrum_pair[NEn-1])<10.0*norm(spectrum_pair[NEn-11])||norm(spectrum_pair[NEn-1])<10.0*norm(spectrum_pair[NEn-21])||norm(spectrum_pair[NEn-1])<10.0*norm(spectrum_pair[NEn-31])) && norm(spectrum_pair[NEn-11])<1.0e-6*E0_max) */
/* 	//if(norm(spectrum_pair[2])<1.0e-6*E0_max && norm(spectrum_pair[NEn-2])<1.0e-6*E0_max) */
/* 	for(int nEn=0; nEn<t_orbits[0].size(); nEn++)	DA_local[nEn] += spectrum_pair[nEn]; */

      //if(norm(spectrum_pairY[NEn-1])<10.0*norm(spectrum_pairY[NEn-31]))
      //if((norm(spectrum_pairY[NEn-1])<10.0*norm(spectrum_pairY[NEn-11])||norm(spectrum_pairY[NEn-1])<10.0*norm(spectrum_pairY[NEn-21])||norm(spectrum_pairY[NEn-1])<10.0*norm(spectrum_pairY[NEn-31])) && norm(spectrum_pairY[0])<norm(spectrum_pairY[NEn-10]))
     //      if(norm(spectrum_pairY[NEn-1])<10.0*norm(spectrum_pairY[NEn-11])||norm(spectrum_pairY[NEn-1])<10.0*norm(spectrum_pairY[NEn-21])||norm(spectrum_pairY[NEn-1])<10.0*norm(spectrum_pairY[NEn-31]))

      //      if((norm(spectrum_pairY[NEn-1])<10.0*norm(spectrum_pairY[NEn-11])||norm(spectrum_pairY[NEn-1])<10.0*norm(spectrum_pairY[NEn-21])||norm(spectrum_pairY[NEn-1])<10.0*norm(spectrum_pairY[NEn-31])) && norm(spectrum_pairY[NEn-11])<1.0e-6)
     
/*       if((norm(spectrum_pairY[NEn-1])<10.0*norm(spectrum_pairY[NEn-11])||norm(spectrum_pairY[NEn-1])<10.0*norm(spectrum_pairY[NEn-21])||norm(spectrum_pairY[NEn-1])<10.0*norm(spectrum_pairY[NEn-31])) && norm(spectrum_pairY[NEn-11])<1.0e-6*E0_max) */
/*       //if(norm(spectrum_pairY[2])<1.0e-6*E0_max && norm(spectrum_pairY[NEn-2])<1.0e-6*E0_max) */
/* 	for(int nEn=0; nEn<t_orbits[0].size(); nEn++)	DA_localY[nEn] += spectrum_pairY[nEn]; */

/*       if(antistokes_passed==false) */
/* 	std::cout << "Error in locating anti-Stokes transition in cycle pair (" << orbit_labels[ns][0] << ',' << orbit_labels[ns][1] << ',' << orbit_labels[ns][2] << ") & (" << orbit_labels[ns_pair][0] << ',' << orbit_labels[ns_pair][1] << ',' << orbit_labels[ns_pair][2] << ")\n";   */



      if(norm(spectrum_pair[2])<1.0e-6*E0_max && norm(spectrum_pair[NEn-2])<1.0e-6*E0_max && norm(spectrum_pairY[2])<1.0e-6*E0_max && norm(spectrum_pairY[NEn-2])<1.0e-6*E0_max)
	for(int nEn=0; nEn<t_orbits[0].size(); nEn++) {
	  DA_local[nEn] += spectrum_pair[nEn];
	  DA_localY[nEn] += spectrum_pairY[nEn];
	}

    }
  }
}

// Calculates harmonic spectrum resulting from previously specified orbits using a combination of the Saddle-Point and Uniform Approx
void qo_response_2d::Dipole_Acceleration_Combined() {
  bool found_pair;
  int ns_pair;
  dcmplx t, tp;
  double epsilon = 1.0e-10;
  dcmplx DetHess, DetHess_prev[2];
  dcmplx spread, spread_prev[2];
  double sign_spreadHess[2];
  dcmplx spreadHess[2];
  dcmplx pre_dtme = qo_response_1d::DTME_Prefactor();
  dcmplx dtmex[2], dtmey[2], dtmeSx, dtmeSy;
  double Energy;
  double min_t = std::min(real(tp_orbits[0][NEn-1]),real(tp_orbits[1][NEn-1]));
  double max_t = std::max(real(t_orbits[t_orbits.size()-2][NEn-1]),real(t_orbits[t_orbits.size()-1][NEn-1]));
  double dt = (max_t-min_t)/population.size();
  double pop_t, pop_tp;
  dcmplx cmplxAmp[2], cmplxAmp_plus, cmplxAmp_minus;
  double short_shift=0.0, long_shift=0.0;
  dcmplx action[2], action_plus, action_minus;
  double sign_spectrum, eta;
  dcmplx sqrt_z, minus_z;
  dcmplx Airy, dzAiry;
  int ierr = 0, argflag = 0, modflag = 0;
  dcmplx pre_factor = sqrt(6.0*pi)*exp(i*0.25*pi);
  double diff_im_S, prev_diff_im_S;
  bool antistokes_passed;
  int nEn_aS;
  double phase_wrap[2], phase_wrap_prev[2], phase_unwrap[2], phase_unwrap_prev[2];
  double phase_wrapY[2], phase_wrap_prevY[2], phase_unwrapY[2], phase_unwrap_prevY[2];
  double pair_shift, pair_shiftY;
  bool crossing_passed;
  int nEn_C;
  dcmplx shift_factor;
  double dEn;
  int nEn_1, nEn_2;
  double x1, x2, y1, y2, k1, k2, a2, b2;
  double x_frac;
  double E0_max = sqrt(pow(pLP[0]->E0,2)+pow(pLP[1]->E0,2));
  for(int nEn=0; nEn<NEn; nEn++)  {DA_local[nEn] = 0.0;  DA_localY[nEn] = 0.0;}
  for(int ns=0; ns<t_orbits.size(); ns++) {
    if(orbit_labels[ns][4]==-1 && orbit_labels[ns][6]==1) {
    //if(orbit_labels[ns][4]==-1) {
      found_pair = false;
      for(ns_pair=0; ns_pair<t_orbits.size(); ns_pair++)
	if(orbit_labels[ns_pair][0]==orbit_labels[ns][0] && orbit_labels[ns_pair][1]==orbit_labels[ns][1] && orbit_labels[ns_pair][2]==orbit_labels[ns][2] && orbit_labels[ns_pair][3]==orbit_labels[ns][3] && orbit_labels[ns_pair][4]==1 && orbit_labels[ns_pair][6]==1) {
	/* if(orbit_labels[ns_pair][0]==orbit_labels[ns][0] && orbit_labels[ns_pair][1]==orbit_labels[ns][1] && orbit_labels[ns_pair][4]==1) { */
	found_pair = true;
	  break;
	}
      if(found_pair==false) {
	std::cout << "Failed to find pair for orbit " << orbit_labels[ns][0] << ',' << orbit_labels[ns][1] << ',' << orbit_labels[ns][2] << ',' << orbit_labels[ns][3] << ',' << orbit_labels[ns][4] << '\n';
	return;
      }
    }
    else  continue;  // skips back to start of for-loop and advances ns
    eta = (double)(abs(orbit_labels[ns][4]+orbit_labels[ns][5]))-1.0;
    sign_spreadHess[0] = 1.0;  sign_spreadHess[1] = 1.0;
    prev_diff_im_S = 0.0;  antistokes_passed = false;
    phase_wrap_prev[0] = 0.0;  phase_wrap_prev[1] = 0.0;  phase_wrap[0] = 0.0;  phase_wrap[1] = 0.0;
    phase_unwrap_prev[0] = 0.0;  phase_unwrap_prev[1] = 0.0;  phase_unwrap[0] = 0.0;  phase_unwrap[1] = 0.0;
    phase_wrap_prevY[0] = 0.0;  phase_wrap_prevY[1] = 0.0;  phase_wrapY[0] = 0.0;  phase_wrapY[1] = 0.0;
    phase_unwrap_prevY[0] = 0.0;  phase_unwrap_prevY[1] = 0.0;  phase_unwrapY[0] = 0.0;  phase_unwrapY[1] = 0.0;
    crossing_passed = false;
    nEn_aS = 0;  nEn_C = t_orbits[0].size()-3; //nE_C = 0;
    for(int nEn=0; nEn<NEn; nEn++) {
      Energy = En[nEn];
      // short orbit
      t = t_orbits[ns][nEn];
      tp = tp_orbits[ns][nEn];
      DTME(dtmex[0], dtmey[0], px(t,tp)+Ax(tp), py(t,tp)+Ay(tp), pre_dtme, atan2(real(Ey(tp)),real(Ex(tp))));
      DTME_singular(dtmeSx, dtmeSy, px(t,tp)+Ax(tp), py(t,tp)+Ay(tp), pre_dtme, atan2(real(Ey(tp)),real(Ex(tp))));
      DetHess = Determinant_Hessian(t,tp);
      spread = epsilon+i*(t-tp);
      if(nEn>0) {
	if(abs(arg(DetHess)-arg(DetHess_prev[0]))>pi)  sign_spreadHess[0] *= -1.0;
	if(abs(arg(spread)-arg(spread_prev[0]))>pi)  sign_spreadHess[0] *= -1.0;
      }
      DetHess_prev[0] = DetHess;  spread_prev[0] = spread;
      spreadHess[0] = sign_spreadHess[0]/std::polar(sqrt(pow(abs(spread),3)*abs(DetHess)), 0.5*(3.0*arg(spread)+arg(DetHess)));
      //spreadHess[0] = sign_spreadHess[0]/pow(spread,1.5)/sqrt(DetHess);
      qo_response_1d::Fit_Population(pop_t, pop_tp, t, tp, min_t, dt);
      cmplxAmp[0] = i*2.0*pow(pi*Energy,2)*spreadHess[0]*pop_t*pop_tp*(Ex(tp)*dtmeSx+Ey(tp)*dtmeSy);
      action[0] = S(t, tp, Energy);
      // long orbit
      t = t_orbits[ns_pair][nEn];
      tp = tp_orbits[ns_pair][nEn];
      DTME(dtmex[1], dtmey[1], px(t,tp)+Ax(tp), py(t,tp)+Ay(tp), pre_dtme, atan2(real(Ey(tp)),real(Ex(tp))));
      DTME_singular(dtmeSx, dtmeSy, px(t,tp)+Ax(tp), py(t,tp)+Ay(tp), pre_dtme, atan2(real(Ey(tp)),real(Ex(tp))));
      DetHess = Determinant_Hessian(t,tp);
      spread = epsilon+i*(t-tp);
      if(nEn>0) {
	if(abs(arg(DetHess)-arg(DetHess_prev[1]))>pi)  sign_spreadHess[1] *= -1.0;
	if(abs(arg(spread)-arg(spread_prev[1]))>pi)  sign_spreadHess[1] *= -1.0;
      }
      DetHess_prev[1] = DetHess;  spread_prev[1] = spread;
      spreadHess[1] = sign_spreadHess[1]/std::polar(sqrt(pow(abs(spread),3)*abs(DetHess)), 0.5*(3.0*arg(spread)+arg(DetHess)));
      //spreadHess[1] = sign_spreadHess[1]/pow(spread,1.5)/sqrt(DetHess);
      qo_response_1d::Fit_Population(pop_t, pop_tp, t, tp, min_t, dt);
      cmplxAmp[1] = i*2.0*pow(pi*Energy,2)*spreadHess[1]*pop_t*pop_tp*(Ex(tp)*dtmeSx+Ey(tp)*dtmeSy);
      action[1] = S(t, tp, Energy);
      if(nEn==0) {
	// potentially a problem here when both short and long orbit 'groups' form >1 band each
	short_shift = 0.0;  long_shift = 0.0;  // reset values
/* 	if(orbit_labels[ns][1]!=1 && orbit_labels[ns][3]!=1) */
/* 	  if(arg(spreadHess[1])-arg(spreadHess[0])>pi)  long_shift -= pi;  // if phase(long)>phase(short)+pi subract pi (but keep phase(long)>phase(short)) */
/* 	  else if(arg(spreadHess[1])-arg(spreadHess[0])<-pi) {short_shift -=pi;  long_shift += pi;}  // if phase(long)+pi<phase(short) subract pi from phase(short) and add pi to phase(long) */
/* 	  else if(arg(spreadHess[1])-arg(spreadHess[0])<0.0)  long_shift += pi;  // if instead phase(long)<phase(short) just add pi to phase(long) */
/* 	//if(orbit_labels[ns][1]==1 || orbit_labels[ns][1]==3)  {short_shift += pi;  long_shift -= pi;} */
/* 	if(arg(spreadHess[1])-arg(spreadHess[0])>pi)  long_shift -= pi;  // if phase(long)>phase(short)+pi subract pi (but keep phase(long)>phase(short)) */
/* 	else if(arg(spreadHess[1])-arg(spreadHess[0])<-pi) {short_shift -=pi;  long_shift += pi;}  // if phase(long)+pi<phase(short) subract pi from phase(short) and add pi to phase(long) */
/* 	else if(arg(spreadHess[1])-arg(spreadHess[0])<0.0)  long_shift += pi;  // if instead phase(long)<phase(short) just add pi to phase(long) */

	if((pLP[0]->I0>1e10*Wcm2_to_au && orbit_labels[ns][1]>1) || (pLP[1]->I0>1e10*Wcm2_to_au && orbit_labels[ns][3]>1))
	  if(arg(spreadHess[1])-arg(spreadHess[0])>pi)  long_shift -= pi;  // if phase(long)>phase(short)+pi subract pi (but keep phase(long)>phase(short))
	  else if(arg(spreadHess[1])-arg(spreadHess[0])<-pi) {short_shift -=pi;  long_shift += pi;}  // if phase(long)+pi<phase(short) subract pi from phase(short) and add pi to phase(long)
	  else if(arg(spreadHess[1])-arg(spreadHess[0])<0.0)  long_shift += pi;  // if instead phase(long)<phase(short) just add pi to phase(long)
      }
      cmplxAmp[0] *= exp(i*short_shift);
      cmplxAmp[1] *= exp(i*long_shift);
      spectrum_short[nEn] = cmplxAmp[0]*conj(dtmex[0])*exp(i*action[0]);
      spectrum_long[nEn] = cmplxAmp[1]*conj(dtmex[1])*exp(i*action[1]);
      spectrum_shortY[nEn] = cmplxAmp[0]*conj(dtmey[0])*exp(i*action[0]);
      spectrum_longY[nEn] = cmplxAmp[1]*conj(dtmey[1])*exp(i*action[1]);
      // uniform approx
      diff_im_S = imag(action[0])-imag(action[1]);
      action_plus = 0.5 * (action[0]+action[1]);
      action_minus = 0.5 * (action[0]-action[1]);
      if(antistokes_passed==false && ((diff_im_S>0.0&&prev_diff_im_S<0.0) || (diff_im_S<0.0&&prev_diff_im_S>0.0))) {
	antistokes_passed = true;
	nEn_aS = nEn;
      }
      if(antistokes_passed==false && ((diff_im_S>0.0&&prev_diff_im_S<0.0) || (diff_im_S<0.0&&prev_diff_im_S>0.0)))  antistokes_passed = true;
      if(antistokes_passed==false) {
	sign_spectrum = -1.0;
	sqrt_z = -pow(1.5*action_minus, 1.0/3.0)*exp(i*eta*pi/3.0);
	minus_z = -pow(sqrt_z, 2);
      }
      else {
	sign_spectrum = 1.0;
	sqrt_z = pow(1.5*action_minus, 1.0/3.0);
	minus_z = -pow(sqrt_z, 2);
      }
      __airy_functions_complex_double_airy_aic(&minus_z, &Airy, &dzAiry, &ierr, &argflag, &modflag);
      cmplxAmp_plus = 0.5 * (cmplxAmp[0]*conj(dtmex[0])+i*cmplxAmp[1]*conj(dtmex[1]));
      cmplxAmp_minus = 0.5 * (cmplxAmp[0]*conj(dtmex[0])-i*cmplxAmp[1]*conj(dtmex[1]));
      spectrum_pair[nEn] = sign_spectrum*pre_factor*sqrt(action_minus)*exp(i*action_plus)*(cmplxAmp_minus/sqrt_z*Airy-i*cmplxAmp_plus/minus_z*dzAiry);
      cmplxAmp_plus = 0.5 * (cmplxAmp[0]*conj(dtmey[0])+i*cmplxAmp[1]*conj(dtmey[1]));
      cmplxAmp_minus = 0.5 * (cmplxAmp[0]*conj(dtmey[0])-i*cmplxAmp[1]*conj(dtmey[1]));
      spectrum_pairY[nEn] = sign_spectrum*pre_factor*sqrt(action_minus)*exp(i*action_plus)*(cmplxAmp_minus/sqrt_z*Airy-i*cmplxAmp_plus/minus_z*dzAiry);
      prev_diff_im_S = diff_im_S;

      // find crossing - energy at which final hump of half-cycle spectra is at a maximum
      if(crossing_passed==false) {
	phase_wrap_prev[0] = phase_wrap[0];  phase_unwrap_prev[0] = phase_unwrap[0];
  	phase_wrap[0] = arg(spectrum_short[nEn]);
	if(nEn==0) {phase_unwrap[0] = phase_wrap[0];  phase_unwrap_prev[0] = 0.0;}
	else if(phase_wrap[0]-phase_wrap_prev[0]<-pi)  phase_unwrap[0] = phase_unwrap_prev[0]+(phase_wrap[0]+2.0*pi)-phase_wrap_prev[0];
	else if(phase_wrap[0]-phase_wrap_prev[0]>pi)  phase_unwrap[0] = phase_unwrap_prev[0]+(phase_wrap[0]-2.0*pi)-phase_wrap_prev[0];
	else  phase_unwrap[0] = phase_unwrap_prev[0]+phase_wrap[0]-phase_wrap_prev[0];
	phase_wrap_prev[1] = phase_wrap[1];  phase_unwrap_prev[1] = phase_unwrap[1];
	phase_wrap[1] = arg(spectrum_long[nEn]);
	if(nEn==0) {
	  phase_unwrap[1] = phase_wrap[1];  phase_unwrap_prev[1] = 0.0;
	  if(phase_wrap[1]<phase_wrap[0])  phase_unwrap[1] += 2.0*pi;
	}
	else if(phase_wrap[1]-phase_wrap_prev[1]<-pi)  phase_unwrap[1] = phase_unwrap_prev[1]+(phase_wrap[1]+2.0*pi)-phase_wrap_prev[1];
	else if(phase_wrap[1]-phase_wrap_prev[1]>pi)  phase_unwrap[1] = phase_unwrap_prev[1]+(phase_wrap[1]-2.0*pi)-phase_wrap_prev[1];
	else  phase_unwrap[1] = phase_unwrap_prev[1]+phase_wrap[1]-phase_wrap_prev[1];
	phase_wrap_prevY[0] = phase_wrapY[0];  phase_unwrap_prevY[0] = phase_unwrapY[0];
  	phase_wrapY[0] = arg(spectrum_shortY[nEn]);
	if(nEn==0) {phase_unwrapY[0] = phase_wrapY[0];  phase_unwrap_prevY[0] = 0.0;}
	else if(phase_wrapY[0]-phase_wrap_prevY[0]<-pi)  phase_unwrapY[0] = phase_unwrap_prevY[0]+(phase_wrapY[0]+2.0*pi)-phase_wrap_prevY[0];
	else if(phase_wrapY[0]-phase_wrap_prevY[0]>pi)  phase_unwrapY[0] = phase_unwrap_prevY[0]+(phase_wrapY[0]-2.0*pi)-phase_wrap_prevY[0];
	else  phase_unwrapY[0] = phase_unwrap_prevY[0]+phase_wrapY[0]-phase_wrap_prevY[0];
	phase_wrap_prevY[1] = phase_wrapY[1];  phase_unwrap_prevY[1] = phase_unwrapY[1];
	phase_wrapY[1] = arg(spectrum_longY[nEn]);
	if(nEn==0) {
	  phase_unwrapY[1] = phase_wrapY[1];  phase_unwrap_prevY[1] = 0.0;
	  if(phase_wrapY[1]<phase_wrapY[0])  phase_unwrapY[1] += 2.0*pi;
	}
	else if(phase_wrapY[1]-phase_wrap_prevY[1]<-pi)  phase_unwrapY[1] = phase_unwrap_prevY[1]+(phase_wrapY[1]+2.0*pi)-phase_wrap_prevY[1];
	else if(phase_wrapY[1]-phase_wrap_prevY[1]>pi)  phase_unwrapY[1] = phase_unwrap_prevY[1]+(phase_wrapY[1]-2.0*pi)-phase_wrap_prevY[1];
	else  phase_unwrapY[1] = phase_unwrap_prevY[1]+phase_wrapY[1]-phase_wrap_prevY[1];

	if(antistokes_passed==true && ((phase_unwrap[1]-phase_unwrap[0]>0.0&&phase_unwrap_prev[1]-phase_unwrap_prev[0]<0.0) || (phase_unwrap[1]-phase_unwrap[0]<0.0&&phase_unwrap_prev[1]-phase_unwrap_prev[0]>0.0))) {
	  crossing_passed = true;
	  nEn_C = nEn;
	}
      }
      // can use nE_C for both polarisations?

      // Above the crossing, phase of both average of SPA short&long gradually shifts away from UA. This removes the shift by comparing SPA average to UA.
      if(crossing_passed==false) {
	spectrum_short[nEn] *= exp(i*(arg(spectrum_pair[nEn])-0.5*(phase_unwrap[0]+phase_unwrap[1])));
	spectrum_long[nEn] *= exp(i*(arg(spectrum_pair[nEn])-0.5*(phase_unwrap[0]+phase_unwrap[1])));
	spectrum_shortY[nEn] *= exp(i*(arg(spectrum_pairY[nEn])-0.5*(phase_unwrapY[0]+phase_unwrapY[1])));
	spectrum_longY[nEn] *= exp(i*(arg(spectrum_pairY[nEn])-0.5*(phase_unwrapY[0]+phase_unwrapY[1])));
      }
    }

    if(norm(spectrum_pair[2])<1.0e-6*E0_max && norm(spectrum_pair[NEn-2])<1.0e-6*E0_max && norm(spectrum_pairY[2])<1.0e-6*E0_max && norm(spectrum_pairY[NEn-2])<1.0e-6*E0_max) {

      for(int nEn=0; nEn<NEn; nEn++) {
	  DA_local[nEn] += conj(spectrum_pair[nEn]);
	  DA_localY[nEn] += conj(spectrum_pairY[nEn]);
      }

      /* for(int nEn=0; nEn<nEn_C; nEn++) { */
/* 	// At the crossing, phase of average of SPA short&long not quite equal to UA. This performs slight adjustment to match them up. */
/* 	spectrum_short[nEn] *= exp(i*(arg(spectrum_short[nEn_C+1]+spectrum_long[nEn_C+1])-arg(spectrum_pair[nEn_C+1]))); */
/* 	spectrum_long[nEn] *= exp(i*(arg(spectrum_short[nEn_C+1]+spectrum_long[nEn_C+1])-arg(spectrum_pair[nEn_C+1]))); */
/* 	// Adjust the magnitude of SPA short&long above crossing to match UA (with universal scaling factor taken at crossing) -> derivative wrt E discontinuous at crossing */
/* 	spectrum_short[nEn] = abs(spectrum_short[nEn_C])/abs(spectrum_pair[nEn_C]) * abs(spectrum_pair[nEn]) * exp(i*arg(spectrum_short[nEn])); */
/* 	spectrum_long[nEn] = abs(spectrum_long[nEn_C])/abs(spectrum_pair[nEn_C]) * abs(spectrum_pair[nEn]) * exp(i*arg(spectrum_long[nEn])); */
/* 	// Further adjust both together such that magnitude of their sum matches with UA at the crossing -> s/l spectra discontinuous at crossing (along with their derivatives) */
/* 	pair_shift = abs(spectrum_pair[nEn])/abs(spectrum_short[nEn]+spectrum_long[nEn]); */
/* 	spectrum_short[nEn] *= pair_shift; */
/* 	spectrum_long[nEn] *= pair_shift; */
/* 	// At the crossing, phase of average of SPA short&long not quite equal to UA. This performs slight adjustment to match them up. */
/* 	spectrum_shortY[nEn] *= exp(i*(arg(spectrum_shortY[nEn_C+1]+spectrum_longY[nEn_C+1])-arg(spectrum_pairY[nEn_C+1]))); */
/* 	spectrum_longY[nEn] *= exp(i*(arg(spectrum_shortY[nEn_C+1]+spectrum_longY[nEn_C+1])-arg(spectrum_pairY[nEn_C+1]))); */
/* 	// Adjust the magnitude of SPA short&long above crossing to match UA (with universal scaling factor taken at crossing) -> derivative wrt E discontinuous at crossing */
/* 	spectrum_shortY[nEn] = abs(spectrum_shortY[nEn_C])/abs(spectrum_pairY[nEn_C]) * abs(spectrum_pairY[nEn]) * exp(i*arg(spectrum_shortY[nEn])); */
/* 	spectrum_longY[nEn] = abs(spectrum_longY[nEn_C])/abs(spectrum_pairY[nEn_C]) * abs(spectrum_pairY[nEn]) * exp(i*arg(spectrum_longY[nEn])); */
/* 	// Further adjust both together such that magnitude of their sum matches with UA at the crossing -> s/l spectra discontinuous at crossing (along with their derivatives) */
/* 	pair_shiftY = abs(spectrum_pairY[nEn])/abs(spectrum_shortY[nEn]+spectrum_longY[nEn]); */
/* 	spectrum_shortY[nEn] *= pair_shiftY; */
/* 	spectrum_longY[nEn] *= pair_shiftY; */
/*       } */
/*       // Smooth magnitude (and also its derivative wrt energy) using spline interpolation */
/*       dEn = En[1]-En[2]; */
/*       nEn_1 = (int)std::min(nEn_C+2.0*eV_to_au/dEn, t_orbits[0].size()-3.0);  // -2eV from E_C */
/*       nEn_2 = (int)std::max(nEn_C-2.0*eV_to_au/dEn, 3.0);  // +2eV from E_C */
/*       x1 = En[nEn_1]; */
/*       x2 = En[nEn_2]; */
/*       // short orbit along x */
/*       y1 = abs(spectrum_short[nEn_1]); */
/*       y2 = abs(spectrum_short[nEn_2]); */
/*       k1 = (abs(spectrum_short[nEn_1-1])-y1)/(En[nEn_1-1]-x1); */
/*       k2 = (abs(spectrum_short[nEn_2-1])-y2)/(En[nEn_2-1]-x2); */
/*       a2 = k1*(x2-x1)-(y2-y1); */
/*       b2 = -k2*(x2-x1)+(y2-y1); */
/*       for(int nEn=nEn_1-1; nEn>nEn_2; nEn--) { */
/* 	x_frac = (En[nEn]-x1)/(x2-x1); */
/* 	spectrum_short[nEn] *= ((1-x_frac)*y1 + x_frac*y2 + x_frac*(1.0-x_frac)*(a2*(1.0-x_frac)+b2*x_frac))/abs(spectrum_short[nEn]); */
/*       } */
/*       // long orbit along x */
/*       y1 = abs(spectrum_long[nEn_1]); */
/*       y2 = abs(spectrum_long[nEn_2]); */
/*       k1 = (abs(spectrum_long[nEn_1-1])-y1)/(En[nEn_1-1]-x1); */
/*       k2 = (abs(spectrum_long[nEn_2-1])-y2)/(En[nEn_2-1]-x2); */
/*       a2 = k1*(x2-x1)-(y2-y1); */
/*       b2 = -k2*(x2-x1)+(y2-y1); */
/*       for(int nEn=nEn_1-1; nEn>nEn_2; nEn--) { */
/* 	x_frac = (En[nEn]-x1)/(x2-x1); */
/* 	spectrum_long[nEn] *= ((1-x_frac)*y1 + x_frac*y2 + x_frac*(1.0-x_frac)*(a2*(1.0-x_frac)+b2*x_frac))/abs(spectrum_long[nEn]); */
/*       } */
/*       // short orbit along y */
/*       y1 = abs(spectrum_shortY[nEn_1]); */
/*       y2 = abs(spectrum_shortY[nEn_2]); */
/*       k1 = (abs(spectrum_shortY[nEn_1-1])-y1)/(En[nEn_1-1]-x1); */
/*       k2 = (abs(spectrum_shortY[nEn_2-1])-y2)/(En[nEn_2-1]-x2); */
/*       a2 = k1*(x2-x1)-(y2-y1); */
/*       b2 = -k2*(x2-x1)+(y2-y1); */
/*       for(int nEn=nEn_1-1; nEn>nEn_2; nEn--) { */
/* 	x_frac = (En[nEn]-x1)/(x2-x1); */
/* 	spectrum_shortY[nEn] *= ((1-x_frac)*y1 + x_frac*y2 + x_frac*(1.0-x_frac)*(a2*(1.0-x_frac)+b2*x_frac))/abs(spectrum_shortY[nEn]); */
/*       } */
/*       // long orbit along y */
/*       y1 = abs(spectrum_longY[nEn_1]); */
/*       y2 = abs(spectrum_longY[nEn_2]); */
/*       k1 = (abs(spectrum_longY[nEn_1-1])-y1)/(En[nEn_1-1]-x1); */
/*       k2 = (abs(spectrum_longY[nEn_2-1])-y2)/(En[nEn_2-1]-x2); */
/*       a2 = k1*(x2-x1)-(y2-y1); */
/*       b2 = -k2*(x2-x1)+(y2-y1); */
/*       for(int nEn=nEn_1-1; nEn>nEn_2; nEn--) { */
/* 	x_frac = (En[nEn]-x1)/(x2-x1); */
/* 	spectrum_longY[nEn] *= ((1-x_frac)*y1 + x_frac*y2 + x_frac*(1.0-x_frac)*(a2*(1.0-x_frac)+b2*x_frac))/abs(spectrum_longY[nEn]); */
/*       } */
/*       for(int nEn=0; nEn<NEn; nEn++) { */
/* 	if(traj=='S'||traj=='B'||traj=='A') { */
/* 	  DA_local[nEn] += conj(spectrum_short[nEn]);  // conj() as was previously giving time-reversed emission */
/* 	  DA_localY[nEn] += conj(spectrum_shortY[nEn]); */
/* 	} */
/* 	if(traj=='L'||traj=='B'||traj=='A') { */
/* 	  DA_local[nEn] += conj(spectrum_long[nEn]); */
/* 	  DA_localY[nEn] += conj(spectrum_longY[nEn]); */
/* 	} */
/*       } */
    }
  }
}

// Function to translate calculated spectrum onto w axis specified in pT object
void qo_response_2d::Fit_to_wAxis() {
  double wEn2, wEn1, wEn0, wEnm1;
  for(int nw=0; nw<=pT->Nt/2; nw++) {
    if(pT->w[nw]>En[1]) {   //if(pT->w[nw]>En[0]) {
      DA[nw] = 0.0;  // remember En goes high->low
      DA2[nw] = 0.0;
    }
    else if(pT->w[nw]<En[NEn-1]) {
      DA[nw] = 0.0;
      DA2[nw] = 0.0;
    }
    else if(pT->w[nw]<En[NEn-2]) {
      DA[nw] = DA_local[NEn-1]*(En[NEn-2]-pT->w[nw])/(En[NEn-2]-En[NEn-1])+DA_local[NEn-2]*(pT->w[nw]-En[NEn-1])/(En[NEn-2]-En[NEn-1]);
      DA2[nw] = DA_localY[NEn-1]*(En[NEn-2]-pT->w[nw])/(En[NEn-2]-En[NEn-1])+DA_localY[NEn-2]*(pT->w[nw]-En[NEn-1])/(En[NEn-2]-En[NEn-1]);
    }
    else if(pT->w[nw]<En[NEn-3]) {
      DA[nw] = DA_local[NEn-2]*(En[NEn-3]-pT->w[nw])/(En[NEn-3]-En[NEn-2])+DA_local[NEn-3]*(pT->w[nw]-En[NEn-2])/(En[NEn-3]-En[NEn-2]);
      DA2[nw] = DA_localY[NEn-2]*(En[NEn-3]-pT->w[nw])/(En[NEn-3]-En[NEn-2])+DA_localY[NEn-3]*(pT->w[nw]-En[NEn-2])/(En[NEn-3]-En[NEn-2]);
    }
    else {
      for(int nEn=NEn-3; nEn>=1; nEn--)   //for(int nEn=nEn_track; nEn>=0; nEn--)
	if(pT->w[nw]>En[nEn+1] && pT->w[nw]<=En[nEn]) {
	  wEn2 = pT->w[nw]-En[nEn+2];
	  wEn1 = pT->w[nw]-En[nEn+1];
	  wEn0 = pT->w[nw]-En[nEn];
	  wEnm1 = pT->w[nw]-En[nEn-1];
	  // 4 point (Lagrange-)polynomial interpolation
	  DA[nw] = DA_local[nEn+2]*wEn1/(En[nEn+2]-En[nEn+1])*wEn0/(En[nEn+2]-En[nEn])*wEnm1/(En[nEn+2]-En[nEn-1]);
	  DA[nw] += DA_local[nEn+1]*wEn2/(En[nEn+1]-En[nEn+2])*wEn0/(En[nEn+1]-En[nEn])*wEnm1/(En[nEn+1]-En[nEn-1]);
	  DA[nw] += DA_local[nEn]*wEn2/(En[nEn]-En[nEn+2])*wEn1/(En[nEn]-En[nEn+1])*wEnm1/(En[nEn]-En[nEn-1]);
	  DA[nw] += DA_local[nEn-1]*wEn2/(En[nEn-1]-En[nEn+2])*wEn1/(En[nEn-1]-En[nEn+1])*wEn0/(En[nEn-1]-En[nEn]);
	  DA2[nw] = DA_localY[nEn+2]*wEn1/(En[nEn+2]-En[nEn+1])*wEn0/(En[nEn+2]-En[nEn])*wEnm1/(En[nEn+2]-En[nEn-1]);
	  DA2[nw] += DA_localY[nEn+1]*wEn2/(En[nEn+1]-En[nEn+2])*wEn0/(En[nEn+1]-En[nEn])*wEnm1/(En[nEn+1]-En[nEn-1]);
	  DA2[nw] += DA_localY[nEn]*wEn2/(En[nEn]-En[nEn+2])*wEn1/(En[nEn]-En[nEn+1])*wEnm1/(En[nEn]-En[nEn-1]);
	  DA2[nw] += DA_localY[nEn-1]*wEn2/(En[nEn-1]-En[nEn+2])*wEn1/(En[nEn-1]-En[nEn+1])*wEn0/(En[nEn-1]-En[nEn]);
	  break;
	}
    }
  }
  for(int nw=pT->Nt/2+1; nw<pT->Nt/2; nw++) {
    DA[nw] = conj(DA[pT->Nt/2-nw]);
    DA2[nw] = conj(DA2[pT->Nt/2-nw]);
  }
}

/* // Function to translate calculated spectrum onto w axis specified in pT object */
/* void qo_response_2d::Fit_to_wAxis() { */
/*   int nEn_track = NEn-2; */
/*   double wEn2, wEn1, wEn0, wEnm1; */
/*   for(int nw=0; nw<=pT->Nt/2; nw++) { */
/*     if(pT->w[nw]>En[0]) { */
/*       DA[nw] = 0.0;  // remember En goes high->low */
/*       DA2[nw] = 0.0; */
/*     } */
/*     else if(pT->w[nw]<En[NEn-1]) { */
/*       DA[nw] = 0.0; */
/*       DA2[nw] = 0.0; */
/*     } */
/*     else { */
/*       for(int nEn=nEn_track; nEn>=0; nEn--) */
/* 	if(pT->w[nw]>En[nEn+1] && pT->w[nw]<=En[nEn]) { */
/* 	  nEn_track = nEn+1; */
/* 	  wEn2 = pT->w[nw]-En[nEn+2]; */
/* 	  wEn1 = pT->w[nw]-En[nEn+1]; */
/* 	  wEn0 = pT->w[nw]-En[nEn]; */
/* 	  wEnm1 = pT->w[nw]-En[nEn-1]; */
/* 	  // 4 point (Lagrange-)polynomial interpolation */
/* 	  DA[nw] = DA_local[nEn+2]*wEn1/(En[nEn+2]-En[nEn+1])*wEn0/(En[nEn+2]-En[nEn])*wEnm1/(En[nEn+2]-En[nEn-1]); */
/* 	  DA[nw] += DA_local[nEn+1]*wEn2/(En[nEn+1]-En[nEn+2])*wEn0/(En[nEn+1]-En[nEn])*wEnm1/(En[nEn+1]-En[nEn-1]); */
/* 	  DA[nw] += DA_local[nEn]*wEn2/(En[nEn]-En[nEn+2])*wEn1/(En[nEn]-En[nEn+1])*wEnm1/(En[nEn]-En[nEn-1]); */
/* 	  DA[nw] += DA_local[nEn-1]*wEn2/(En[nEn-1]-En[nEn+2])*wEn1/(En[nEn-1]-En[nEn+1])*wEn0/(En[nEn-1]-En[nEn]); */
/* 	  DA2[nw] = DA_localY[nEn+2]*wEn1/(En[nEn+2]-En[nEn+1])*wEn0/(En[nEn+2]-En[nEn])*wEnm1/(En[nEn+2]-En[nEn-1]); */
/* 	  DA2[nw] += DA_localY[nEn+1]*wEn2/(En[nEn+1]-En[nEn+2])*wEn0/(En[nEn+1]-En[nEn])*wEnm1/(En[nEn+1]-En[nEn-1]); */
/* 	  DA2[nw] += DA_localY[nEn]*wEn2/(En[nEn]-En[nEn+2])*wEn1/(En[nEn]-En[nEn+1])*wEnm1/(En[nEn]-En[nEn-1]); */
/* 	  DA2[nw] += DA_localY[nEn-1]*wEn2/(En[nEn-1]-En[nEn+2])*wEn1/(En[nEn-1]-En[nEn+1])*wEn0/(En[nEn-1]-En[nEn]); */
/* 	  break; */
/* 	} */
/*     } */
/*   } */
/*   for(int nw=pT->Nt/2+1; nw<pT->Nt/2; nw++) { */
/*     DA[nw] = conj(DA[pT->Nt/2-nw]); */
/*     DA2[nw] = conj(DA2[pT->Nt/2-nw]); */
/*   } */
/* } */

// Governing interface to calculate dipole response according to specified method
void qo_response_2d::Dipole_Acceleration(dcmplx *EtX, dcmplx *EtY) {
  if(analyse==true)  pLaser->Analyse_Field(pLP, EtX, EtY, findAll);

  //if(analyse==true)  pLaser->Analyse_Field(pLP, EtX, findAll);

  qo_response_1d::Calculate_FieldFactors();
  qo_response_1d::Set_Bounds();
  Build_Energy_Axis();
  Saddlepoint_Search();
  bool tracked = Track_Orbits();
  if(tracked==false) {
    std::cout << "Orbit tracking failed. Exiting...\n";
    return;//exit(1);
  }
  qo_response_1d::Characterise_Orbits();
  Calculate_Population(); 
  if(method=='S')  Dipole_Acceleration_SaddlePoint();
  else if(method=='U')  Dipole_Acceleration_Uniform();
  else if(method=='C')  Dipole_Acceleration_Combined();
  else {
    std::cerr << "Error in method for calculating quantum orbit response, setting to uniform approximation.\n";
    Dipole_Acceleration_Uniform();
  }
  Fit_to_wAxis();
  //for(int nw=0; nw<pT->Nt; nw++)  DA[nw]*=exp(-i*pT->w[nw]*pLP[colourRef]->t0);
}

// As above but for default Et, Et2
void qo_response_2d::Dipole_Acceleration() {
  Dipole_Acceleration(pLaser->Et, pLaser->Et2);
}

// As above but for default E_mpi / Ert (depending on flag) with specified row
void qo_response_2d::Dipole_Acceleration(int nr, bool mpi=false) {
  (mpi==true) ? Dipole_Acceleration(&pLaser->E_mpi[Index2V(nr,0,pT->Nt)], &pLaser->E_mpi2[Index2V(nr,0,pT->Nt)]) : Dipole_Acceleration(&pLaser->Ert[Index2V(nr,0,pT->Nt)], &pLaser->Ert2[Index2V(nr,0,pT->Nt)]);
  //(mpi==true) ? qo_response_1d::Dipole_Acceleration(&pLaser->E_mpi[Index2V(nr,0,pT->Nt)]) : qo_response_1d::Dipole_Acceleration(&pLaser->Ert[Index2V(nr,0,pT->Nt)]);
}

// Dipole acceleration in freq-domain importing row number of specified temporospatial fields
void qo_response_2d::Dipole_Acceleration(dcmplx *ErtX, dcmplx *ErtY, int nr) {
  Dipole_Acceleration(&ErtX[Index2V(nr,0,pT->Nt)], &ErtY[Index2V(nr,0,pT->Nt)]);
}

/***************************************************************************************************************/
#endif
