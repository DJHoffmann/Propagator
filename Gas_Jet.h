/***************************************************************************************************************/
// Header file to contain gas_jet class for use in Propagator.cpp 

/***************************************************************************************************************/
/* Notes:
   Need to add in molar mass (private) and absorption array (public; will also need T.w)
   IMPORTANT: X-Ray absorption data taken from /Inputs/Absorb_Z.dat where Z=particle name/symbol as entered in inputs file
   IMPORTANT: Henke file does not record selected temperature, so always use 295K (which is assumed below)
*/

/***************************************************************************************************************/
#ifndef GAS_JET_H
#define GAS_JET_H

#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<Constants_Conversions.h>
#include<Misc_Functions.h>
#include<Environment.h>
#include<Input.h>

using std::abs;

/***************************************************************************************************************/
// Class to create gas jet number density distribution:
class gas_jet {
 private:
  /***************************************************************************************************************/
  // Private Class Members:
  double nd0;
  double fwhm;
  double gamma;
  double pre_ADK, power_ADK, exp_ADK;
  std::string particle;
  spatial_axes *pX;
  temporal_axes *pT;
  /***************************************************************************************************************/
  // Private Class Function Prototypes:
  void Set_and_Convert();
  void Initialise();
  void Make_Jet();
  void ADK_Factors();
  void Particle_Data(std::string);
  void Check_Particle(std::ifstream &);
  void Import_Particle(std::ifstream &);
  void Import_XRay_Transmission(double *&, double *&, double &, double &, int &, std::ifstream &);
  void Fit_XRay_Absorption(double *, double *, double, double, int);
 public:
  /***************************************************************************************************************/
  // Public Class Members:
  double Ip;
  double zj;
  int n, l, m;
  double *nd;
  double *ionised; // will hold on-axis ionised fraction of sample (integrated over time)
  double *absorb;
  /***************************************************************************************************************/
  // Public Class Function Prototypes:
  gas_jet(double, double, double, spatial_axes &);
  gas_jet(double, int, int, int, spatial_axes &);
  gas_jet(double, double, double, double, int, int, int, spatial_axes &);
  gas_jet(double, double, double, int, int, int, spatial_axes &, temporal_axes &, char *);
  gas_jet(std::ifstream &, spatial_axes &, temporal_axes &, char *);
  ~gas_jet();
  double ADK_Rate(dcmplx);
  double ADK_Rate(dcmplx, dcmplx);
  double ADK_Ionisation(dcmplx *);
  double ADK_Ionisation(dcmplx *, dcmplx *);
  double ADK_Ionisation(dcmplx *, temporal_axes &);
  void ADK_Ionisation(double *, dcmplx *);
  void ADK_Ionisation(double *, dcmplx *, dcmplx *);
  void ADK_Ionisation(double *, dcmplx *, temporal_axes &);
  void Import_Data(const char *);
  void Output(double *, std::ofstream &);
};

/***************************************************************************************************************/
// Class Functions:

// Contructor for gas jet taking only jet data as arguments and assuming 2p-state of argon
gas_jet::gas_jet(double nd0_, double zj_, double fwhm_, spatial_axes &X_) {
  pX = &X_;
  pT = NULL;
  nd0 = nd0_;
  zj = zj_;
  fwhm = fwhm_;
  Ip = 15.76;
  n = 2;  l = 1;  m = 0;
  Initialise();
  Set_and_Convert();
  Make_Jet();
}

// General gas jet constructor taking only atomic data as arguments
gas_jet::gas_jet(double Ip_, int n_, int l_, int m_, spatial_axes &X_) {
  pX = &X_;
  pT = NULL;
  nd0 = 1.0e18;
  zj = 2.0e-3;
  fwhm = 1.0e-3;
  Ip = Ip_;
  n = n_;  l = l_;  m = m_;
  Initialise();
  Set_and_Convert();
  Make_Jet();
}

// General gas jet constructor taking jet and atomic data as arguments
gas_jet::gas_jet(double nd0_, double zj_, double fwhm_, double Ip_, int n_, int l_, int m_, spatial_axes &X_) {
  pX = &X_;
  pT = NULL;
  nd0 = nd0_;
  zj = zj_;
  fwhm = fwhm_;
  Ip = Ip_;
  n = n_;  l = l_;  m = m_;
  Initialise();
  Set_and_Convert();
  Make_Jet();
}

// General gas jet constructor taking jet and atomic data as arguments, also with temporal axes for absorption data
gas_jet::gas_jet(double nd0_, double zj_, double fwhm_, int n_, int l_, int m_, spatial_axes &X_, temporal_axes &T_, char *ConfigPath) {
  pX = &X_;
  pT = &T_;
  nd0 = nd0_;
  zj = zj_;
  fwhm = fwhm_;
  n = n_;  l = l_;  m = m_;
  Initialise();
  Import_Data(ConfigPath);
  Set_and_Convert();
  Make_Jet();
}

// Constructor as above and ConfigPath imported
gas_jet::gas_jet(std::ifstream &INSTREAM, spatial_axes &X_, temporal_axes &T_, char *ConfigPath) {
  pX = &X_;
  pT = &T_;
  
  //double rj = Read_Variable<double>(SOURCE, "rj");

  nd0 = Read_Variable<double>(INSTREAM, "nd0"); // replace with ? a:b statement if backing pressure and distance supplied
  zj = Read_Variable<double>(INSTREAM, "zj");
  fwhm = Read_Variable<double>(INSTREAM, "zw");
  particle = Read_Variable<std::string>(INSTREAM, "GasName");
  Particle_Data(particle);
  n = Read_Variable<int>(INSTREAM, "nGS");  
  l = Read_Variable<int>(INSTREAM, "lGS");
  m = Read_Variable<int>(INSTREAM, "mGS");
  Initialise();
  Import_Data(ConfigPath);
  Set_and_Convert();
  Make_Jet();
}

// Destructor
gas_jet::~gas_jet() {
  delete[] nd;
  delete[] ionised;
  if(pT!=NULL)
    delete[] absorb;
}

// Converts input values to a.u.
void gas_jet::Set_and_Convert() {
  Ip *= eV_to_au;
  nd0 *= pow(100.0/m_to_au, 3);
  zj *= m_to_au;
  fwhm *= m_to_au;
  ADK_Factors();
}

// Initialise memory of gas jet variables
void gas_jet::Initialise() {
  nd = new double[pX->Nz];
  ionised = new double[pX->Nz];
  if(pT!=NULL)
    absorb = new double[pT->Nt];
  else
    absorb = NULL;
}

// Makes the gas jet density profile
void gas_jet::Make_Jet() {
  double env_factor = 4.0*log(2.0)/(fwhm*fwhm);
  for(int nz=0; nz<pX->Nz; nz++)
    nd[nz] = nd0 * exp(-(pX->z[nz]-zj)*(pX->z[nz]-zj)*env_factor);
}

// Provides three atomic numerical factors for ADK calculation
void gas_jet::ADK_Factors() {
  m = 0; // m should be a projection of l onto E field polarisation vector - set to 0 for now
  double n_eff = 1.0/sqrt(2.0*fabs(Ip));
  double l_eff = n_eff-1.0; //0.0;
  double Anl = pow(2.0, 2.0*n_eff) / (n_eff*Gamma(n_eff+l_eff+1) * Gamma(n_eff-l_eff));
  double Blm = (2*l+1) * Factorial(l+abs(m)) / (pow(2.0, abs(m)) * Factorial(abs(m)) * Factorial(l-abs(m)));
  pre_ADK = Anl * Blm * fabs(Ip) * pow(2.0*pow(2.0*fabs(Ip), 1.5), 2.0*n_eff-abs(m)-1); 
  power_ADK = -(2.0*n_eff-abs(m)-1);
  exp_ADK = -2.0/3.0*pow(2.0*fabs(Ip), 1.5); 
}

// Calculates ADK ionisation rate
double gas_jet::ADK_Rate(dcmplx Et) {
  double ADK; 
  double Eabs = fabs(real(Et)); //abs(Et);
  if(Eabs>0.0)  ADK = pre_ADK * pow(Eabs, power_ADK) * exp(exp_ADK/Eabs);
  else  ADK = 0.0;
  return ADK;
}

// Calculates ADK ionisation rate for 2D field
double gas_jet::ADK_Rate(dcmplx EtX, dcmplx EtY) {
  double ADK; 
  double Eabs = sqrt(pow(real(EtX),2)+pow(real(EtY),2));
  if(Eabs>0.0)  ADK = pre_ADK * pow(Eabs, power_ADK) * exp(exp_ADK/Eabs);
  else  ADK = 0.0;
  return ADK;
}

// Calculates total ionisation arising from ADK ionisation rate
double gas_jet::ADK_Ionisation(dcmplx *Et) {
  double exp_factor = 0.0;
  for(int nt=0; nt<pT->Nt; nt++)
    exp_factor -= pT->dt * ADK_Rate(Et[nt]);
  return 1.0 - exp(exp_factor);
}

// Calculates total ionisation arising from ADK ionisation rate
double gas_jet::ADK_Ionisation(dcmplx *EtX, dcmplx *EtY) {
  double exp_factor = 0.0;
  for(int nt=0; nt<pT->Nt; nt++)
    exp_factor -= pT->dt * ADK_Rate(sqrt(pow(abs(EtX[nt]),2)+pow(abs(EtY[nt]),2)));
  return 1.0 - exp(exp_factor);
}

// Calculates total ionisation arising from ADK ionisation rate, taking in temporal axes as argument
double gas_jet::ADK_Ionisation(dcmplx *Et, temporal_axes &T) {
  double exp_factor = 0.0;
  for(int nt=0; nt<T.Nt; nt++)
    exp_factor -= T.dt * ADK_Rate(Et[nt]);
  return 1.0 - exp(exp_factor);
}

// Calculates total ionisation arising from ADK ionisation rate
void gas_jet::ADK_Ionisation(double *population, dcmplx *Et) {
  double exp_factor = 0.0;
  for(int nt=0; nt<pT->Nt; nt++) {
    exp_factor -= pT->dt * ADK_Rate(Et[nt]);
    population[nt] = exp(exp_factor);
  }
}

// Calculates total ionisation arising from ADK ionisation rate
void gas_jet::ADK_Ionisation(double *population, dcmplx *EtX, dcmplx *EtY) {
  double exp_factor = 0.0;
  for(int nt=0; nt<pT->Nt; nt++) {
    exp_factor -= pT->dt * ADK_Rate(sqrt(pow(abs(EtX[nt]),2)+pow(abs(EtY[nt]),2)));
    population[nt] = exp(exp_factor);
  }
}

// Calculates total ionisation arising from ADK ionisation rate, taking in temporal axes as argument
void gas_jet::ADK_Ionisation(double *population, dcmplx *Et, temporal_axes &T) {
  double exp_factor = 0.0;
  for(int nt=0; nt<T.Nt; nt++) {
    exp_factor -= T.dt * ADK_Rate(Et[nt]);
    population[nt] = exp(exp_factor);
  }
}

// Function to take in gas type from file and thus provide ionisation potential & specific heat ratio
void gas_jet::Particle_Data(std::string particle) {
  if(particle.compare("H")==0 || particle.compare("Hydrogen")==0) {
    Ip = 13.5984;
    gamma = 1.41;
  }
  else if(particle.compare("He")==0 || particle.compare("Helium")==0) {
    Ip = 24.5874;
    gamma = 1.63;
  }
  else if(particle.compare("Ne")==0 || particle.compare("Neon")==0) {
    Ip = 21.5645;
    gamma = 1.642;
  }
  else if(particle.compare("Ar")==0 || particle.compare("Argon")==0) {
    Ip = 15.7596;
    gamma = 1.667;
  }
  else if(particle.compare("Kr")==0 || particle.compare("Krypton")==0) {
    Ip = 13.9996;
    gamma = 1.689;
  }
  else if(particle.compare("Xe")==0 || particle.compare("Xenon")==0) {
    Ip = 12.12984;
    gamma = 1.666;
  }
  else {
    std::cerr << "Error in particle selection\n";
#ifdef useMPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(1);
  }
}

// Function to take in gas type from file and thus provide ionization potential & molar mass
void gas_jet::Import_Particle(std::ifstream &INSTREAM) {
  INSTREAM.seekg(0, std::ios::beg);
  INSTREAM.ignore(100,' ');
  INSTREAM >> particle;
  Particle_Data(particle);
}

// Function to check correct absorption file in use
void gas_jet::Check_Particle(std::ifstream &INSTREAM) {
  std::string check;
  INSTREAM.seekg(0, std::ios::beg);
  INSTREAM.ignore(100,' ');
  INSTREAM >> check;
  if(check.compare(particle)!=0) {
    std::cerr << "Error in choice of absorption file\n";
#ifdef useMPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(1);
  }
}

// Function to import x-ray transmission & corresponding energy axis from file (from http://henke.lbl.gov/optical_constants/)
void gas_jet::Import_XRay_Transmission(double *&EnF, double *&transF, double &pressure, double &path, int &Ntrans, std::ifstream &INSTREAM) {
  INSTREAM.seekg(0, std::ios::beg);
  INSTREAM.ignore(100,'=');
  INSTREAM >> pressure;
  INSTREAM.ignore(100,'=');
  INSTREAM >> path;
  INSTREAM.ignore(1000,'\n');  INSTREAM.ignore(1000,'\n');
  while(INSTREAM.good()) {
    INSTREAM.ignore(1000,'\n');
    Ntrans++;
  }
  if(Ntrans==0) {
    std::cerr << "Error in reading in x-ray transmission data\n";
#ifdef useMPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(1);
  }
  INSTREAM.clear();
  INSTREAM.seekg(0, std::ios::beg);
  INSTREAM.ignore(1000,'\n');  INSTREAM.ignore(1000,'\n');
  EnF = new double[Ntrans];
  transF = new double[Ntrans];
  for(int ntrans=0; ntrans<Ntrans; ntrans++) {
    INSTREAM >> EnF[ntrans];
    INSTREAM >> transF[ntrans];
  }
}

// Function to calculate x-ray absorption in a.u. fitted to T.w frequency axis
void gas_jet::Fit_XRay_Absorption(double *EnF, double *transF, double pressure, double path, int Ntrans) {
  int ntrans;
  double En, frac;
  for(int nt=0; nt<1+pT->Nt/2; nt++) {  // first fit transmission profile onto our frequency axes
    En = pT->w[nt]*au_to_eV;
    ntrans = 0;
    if(En<=EnF[0])
      absorb[nt] = transF[0];
    else {
      while(ntrans<Ntrans && EnF[ntrans]<En)
	ntrans++;
      if(ntrans==Ntrans-1)
	absorb[nt] = transF[ntrans];
      else {
	frac = (En-EnF[ntrans-1])/(EnF[ntrans]-EnF[ntrans-1]);
	absorb[nt] = (1.0-frac)*transF[ntrans-1] + frac*transF[ntrans];
      }
    }
  }
  for(int nt=1+pT->Nt/2; nt<pT->Nt; nt++)  // negative frequencies
    absorb[nt] = absorb[pT->Nt-nt];
  double T = 295.0*K_to_au; // Henke file does not record selected temperature - always pick 295K which we use here
  double ndF = pressure/T; // Boltzmann constant = 1 in atomic units
  for(int nt=0; nt<pT->Nt; nt++)  // now convert to absorption per dz per unit number density (test 100 factor for zero transmission)
    absorb[nt] = (absorb[nt]==0) ? 100.0 : -log(absorb[nt])*(pX->dz/path)/ndF;    
}

// Overriding function to import and convert all data 
void gas_jet::Import_Data(const char *ConfigPath) {
  std::string FilePath = ConfigPath;
  FilePath.append("Absorb_"); FilePath.append(particle); FilePath.append(".conf");
  //std::cout << FilePath << '\n';
  std::ifstream ABSORB(FilePath.c_str());
  Import_Particle(ABSORB);
  int Ntrans = 0;
  double pressure, path;
  double *EnF, *transF;
  Import_XRay_Transmission(EnF, transF, pressure, path, Ntrans, ABSORB);
  pressure *= Torr_to_Pa*Pa_to_au;
  path *= 1.0e-2*m_to_au;
  Fit_XRay_Absorption(EnF, transF, pressure, path, Ntrans);
  delete[] EnF;
  delete[] transF;
  ABSORB.close();
}

// Outputs specified array along z-axis
void gas_jet::Output(double *A, std::ofstream &OUTSTREAM) {
  for(int nz=0; nz<pX->Nz; nz++)
    OUTSTREAM << A[nz] << '\t';
  OUTSTREAM << '\n';
}

/***************************************************************************************************************/
#endif


