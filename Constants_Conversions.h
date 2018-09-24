/***************************************************************************************************************/
// Header file to contain universal constants and conversion factors for Propagator_Vacuum_CrankNicolson.cpp 

/***************************************************************************************************************/
// Notes:

/***************************************************************************************************************/
#ifndef CONSTANTS_CONVERSIONS_H
#define CONSTANTS_CONVERSIONS_H

#include<complex>

/***************************************************************************************************************/
// New Definitions:
typedef std::complex<double> dcmplx;

/***************************************************************************************************************/
// Functions Prototypes:
double FS_to_Cycles(double, double);

/***************************************************************************************************************/
// Universal Constants & Conversions Factors:
const dcmplx i(0.0,1.0);
const double pi = acos(-1.0);

const double mole = 6.023e23;
const double kB = 1.38065e-23;

const double Torr_to_atm = 1.3158e-3; 
const double Torr_to_Pa = 133.322;

const double c_au = 137.036;
const double au_to_m = 5.29177e-11;
const double m_to_au = 1.0/au_to_m;
const double au_to_s = 2.4189e-17;
const double s_to_au = 1.0/au_to_s;
const double au_to_eV = 27.212;
const double eV_to_au = 1.0/au_to_eV;
const double au_to_Wcm2 = 3.509e16;
const double Wcm2_to_au = 1.0/au_to_Wcm2;
const double au_to_kg = 9.109e-31;
const double kg_to_au = 1.0/au_to_kg;
const double au_to_K = 3.158e5;
const double K_to_au = 1.0/au_to_K;
const double au_to_Pa = 2.942e13;
const double Pa_to_au = 1.0/au_to_Pa;

/***************************************************************************************************************/
// Functions:

// Function to calculate number of cycles from femtosecond FWHMI
double FS_to_Cycles(double fs, double lambda) {
  return fs/0.024189*5.29177e-2*c_au/lambda;
}

/***************************************************************************************************************/
#endif
