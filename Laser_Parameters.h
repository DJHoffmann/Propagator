/***************************************************************************************************************/
// Header file to contain laser_parameters class for Propagator.cpp 

/***************************************************************************************************************/
// Notes:

/***************************************************************************************************************/
#ifndef LASER_PARAMETERS_H
#define LASER_PARAMETERS_H

#include<fstream>
#include<Constants_Conversions.h>
#include<Input.h>

/***************************************************************************************************************/
// laser_parameters class:
class laser_parameters {
 private:
  /***************************************************************************************************************/
  // Private Class Function Prototypes:
  void Set_and_Convert();
 public:
  /***************************************************************************************************************/
  // Public Class Members:
  double I0, lambda, fwhm, t0, cep, zf, wz0, phi;
  char envelope;
  double E0, T, w0, k0, cycles, B;
  /***************************************************************************************************************/
  // Public Class Function Prototypes:
  laser_parameters();
  laser_parameters(double, double);
  laser_parameters(double, double, double, double, double, char, double);
  laser_parameters(double, double, double, double);
  laser_parameters(double, double, double, double, double, char, double, double);
  laser_parameters(double, double, double, double, double, char, double, double, double);
  void Set_Polarisation_Angle(double);
  double Ponderomotive_Potential();
};

/***************************************************************************************************************/
// Class Functions:

// Default constructor
laser_parameters::laser_parameters() {
  I0 = 0.0;
  lambda = 0.0;
  fwhm = 0.0;
  t0 = 0.0;
  cep = 0.0;
  envelope = 'g';
  zf = 0.0;
  wz0 = 0.0;
  phi = 0.0;
  E0 = 0.0;
  T = 0.0;
  w0 = 0.0;
  k0 = 0.0;
  cycles = 0.0;
  B = 0.0;
}
  
// Constructor for CW temporal field
laser_parameters::laser_parameters(double I0_, double lambda_) {
  I0 = I0_;
  lambda = lambda_;
  fwhm = 1.0e30;
  t0 = 0.0;
  cep = 0.0;
  envelope = 'g';
  zf = 0.0;
  wz0 = 0.0;
  phi = 0.0;
  Set_and_Convert();
}

// Constructor for pulsed temporal field
laser_parameters::laser_parameters(double I0_, double lambda_, double fwhm_, double t0_, double cep_, char envelope_, double phi_=1.0) {
  I0 = I0_;
  lambda = lambda_;
  fwhm = fwhm_;
  t0 = t0_;
  cep = cep_;
  envelope = envelope_;
  zf = 0.0;
  wz0 = 0.0;
  phi = phi_;
  Set_and_Convert();
}

// Constructor for CW temporo-spatial field
laser_parameters::laser_parameters(double I0_, double lambda_, double zf_, double wz0_) {
  I0 = I0_;
  lambda = lambda_;
  fwhm = 1.0e30;
  t0 = 0.0;
  cep = 0.0;
  envelope = 'g';
  zf = zf_;
  wz0 = wz0_;
  phi = 0.0;
  Set_and_Convert();
}

// Constructor for pulsed temporo-spatial field
laser_parameters::laser_parameters(double I0_, double lambda_, double fwhm_, double t0_, double cep_, char envelope_, double zf_, double wz0_) {
  I0 = I0_;
  lambda = lambda_;
  fwhm = fwhm_;
  t0 = t0_;
  cep = cep_;
  envelope = envelope_;
  zf = zf_;
  wz0 = wz0_;
  phi = 0.0;  
  Set_and_Convert();
}

// Constructor for pulsed temporo-spatial field
laser_parameters::laser_parameters(double I0_, double lambda_, double fwhm_, double t0_, double cep_, char envelope_, double zf_, double wz0_, double phi_) {
  I0 = I0_;
  lambda = lambda_;
  fwhm = fwhm_;
  t0 = t0_;
  cep = cep_;
  envelope = envelope_;
  zf = zf_;
  wz0 = wz0_;
  phi = phi_;  
  Set_and_Convert();
}

// Function to convert inputted laser parameters to a.u. and calculate derived parameters
void laser_parameters::Set_and_Convert() {
  I0 *= Wcm2_to_au;
  lambda *= m_to_au;
  fwhm *= s_to_au;
  t0 *= s_to_au;
  cep *= pi;
  zf *= m_to_au;
  wz0 *= m_to_au;
  phi *= pi;
  E0 = sqrt(I0);
  T = lambda/c_au;
  w0 = 2.0*pi/T;
  k0 = w0*c_au;
  cycles = fwhm/lambda*c_au;
  switch(envelope) {
  case 'g':
    B = 2.0*log(2.0)/(fwhm*fwhm);
    break;
  case 'c':
    B = 2.0*acos(pow(0.5,0.25))/fwhm;
    break;
  case 's':
    B = 2.0*acosh(sqrt(2.0))/fwhm;
    break;
  default:
    std::cerr << "Incorrect envelope specifed.\n";
    exit(1);
  }
}

// Function to set (or reset) polarisation angle
void laser_parameters::Set_Polarisation_Angle(double phi_) {
  phi = phi_;
}

// Function to return the ponderomotive potential
double laser_parameters::Ponderomotive_Potential() {
  return pow(E0/(2.0*w0), 2);
}

/***************************************************************************************************************/
// Universal functions on laser_parameters objects

// Function to read in laser parameters from file located at FilePath (not class function as acts on array of such classes)
void Import_Parameters(laser_parameters **&pLP, int &colours, std::ifstream &INSTREAM) { // need & on **LP to maintain memory allocation upon function exit
  double I0, lambda, fwhm, t0, cep, zf, wz0, phi;
  char envelope;
  colours = Read_Variable<int>(INSTREAM, "colours");
  pLP = new laser_parameters*[colours]; // allocate memory based on number of colours found
  for(int nc=1; nc<=colours; nc++) {
    I0 = Read_Variable<double>(INSTREAM, "I0", nc);
    lambda = Read_Variable<double>(INSTREAM, "lambda", nc);
    fwhm = Read_Variable<double>(INSTREAM, "Ifwhm", nc);
    t0 = Read_Variable<double>(INSTREAM, "t0", nc);
    cep = Read_Variable<double>(INSTREAM, "cep", nc);
    envelope = Read_Variable<char>(INSTREAM, "envshape", nc);
    zf = Read_Variable<double>(INSTREAM, "zFocus", nc);
    wz0 = Read_Variable<double>(INSTREAM, "beamwaist", nc);
    phi = Read_Variable<double>(INSTREAM, "polAxis", nc);
    pLP[nc-1] = new laser_parameters(I0, lambda, fwhm, t0, cep, envelope, zf, wz0, phi);
    //std::cout << I0 << '\t' << lambda <<'\t' << fwhm << '\t' << cep << '\t' << envelope << '\t' << zf << '\t' << wz0 << '\n';
  }
}

/* // Function to check polarisation angles and return true/false to set one-/two-dimensional fields */
/* bool Set_Dimension(laser_parameters **pLP, int colours) { */
/*   bool TwoDim = false; */
/*   int nc0; */
/*   double angle; */
/*   for(int nc=1; nc<=colours; nc++)   */
/*     if(pLP[nc-1]->I0!=0.0) { */
/*       angle = pLP[nc-1]->phi; */
/*       nc0 = nc; */
/*       break; */
/*     } */
/*   for(int nc=nc0+1; nc<=colours; nc++)   */
/*     if(pLP[nc-1]->phi!=angle && pLP[nc-1]->I0!=0.0)  TwoDim = true;   */
/*   if(TwoDim==false) { */
/* /\*     if(angle!=0.0) *\/ */
/* /\*       std::cout << "Single laser polarisation axis found, remapping each colour onto x-axis (phi=0)\n"; *\/ */
/*     for(int nc=1; nc<=colours; nc++)   */
/*       pLP[nc-1]->Set_Polarisation_Angle(0.0); */
/*   } */
/*   return TwoDim; */
/* } */

// Function to check polarisation angles and return true/false to set one-/two-dimensional fields
bool Set_Dimension(laser_parameters **pLP, int &colours) {
  int nc0;
  bool contract = true;
  bool TwoDim = false;
  double angle;
  for(int nc=1; nc<=colours; nc++)  
    if(pLP[nc-1]->I0!=0.0) {
      angle = pLP[nc-1]->phi;
      nc0 = nc;
      break;
    }
  for(int nc=nc0+1; nc<=colours; nc++)  
    if(pLP[nc-1]->I0!=0.0) {
      contract = false;  
      break;
    }
  if(contract==true) {
    if(nc0!=1) {
      pLP[0]->I0 = pLP[nc0-1]->I0;
      pLP[0]->lambda = pLP[nc0-1]->lambda;
      pLP[0]->fwhm = pLP[nc0-1]->fwhm;
      pLP[0]->t0 = pLP[nc0-1]->t0;
      pLP[0]->cep = pLP[nc0-1]->cep;
      pLP[0]->zf = pLP[nc0-1]->zf;
      pLP[0]->wz0 = pLP[nc0-1]->wz0;
      pLP[0]->phi = pLP[nc0-1]->phi;
      pLP[0]->envelope = pLP[nc0-1]->envelope;
      pLP[0]->E0 = pLP[nc0-1]->E0;
      pLP[0]->T = pLP[nc0-1]->T;
      pLP[0]->w0 = pLP[nc0-1]->w0;
      pLP[0]->k0 = pLP[nc0-1]->k0;
      pLP[0]->cycles = pLP[nc0-1]->cycles;
      pLP[0]->B = pLP[nc0-1]->B;
    }
    colours = 1;
  }
  else
    for(int nc=nc0+1; nc<=colours; nc++)  
      if(pLP[nc-1]->phi!=angle && pLP[nc-1]->I0!=0.0) {
	TwoDim = true;
	break;
      }
  if(TwoDim==false) {
/*     if(angle!=0.0) */
/*       std::cout << "Single laser polarisation axis found, remapping each colour onto x-axis (phi=0)\n"; */
    for(int nc=1; nc<=colours; nc++)  
      pLP[nc-1]->Set_Polarisation_Angle(0.0);
  }
  return TwoDim;
}

/***************************************************************************************************************/
#endif

