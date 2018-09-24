/***************************************************************************************************************/
// Header file to contain plasma defocussing / blue-shifting functions due to free electrons for Propagator.cpp 

/***************************************************************************************************************/
// Notes:

/***************************************************************************************************************/
#ifndef PLASMA_DEFOCUSSING_H
#define PLASMA_DEFOCUSSING_H

#include<Misc_Functions.h>
#include<Gas_Jet.h>

/***************************************************************************************************************/
// Function Prototypes:
double Plasma_Frequency_Squared(double, double);
void RungeKutta_FirstOrder(dcmplx *, double &, dcmplx *, gas_jet &, double, const double, const int, const int);
void RungeKutta_FirstOrder(dcmplx *, double &, dcmplx *, gas_jet &, double, const double, const int, const int, const int);
void RungeKutta_SecondOrder(dcmplx **, double *, dcmplx *, gas_jet &, double, const double, const double, const int, const int);
void RungeKutta_FourthOrder(dcmplx **, double *, dcmplx *, gas_jet &, double, const double, const double, const int, const int);
void Gills_FourthOrder(dcmplx **, double *, dcmplx *, gas_jet &, double, const double, const double, const int, const int);
void Plasma_Defocussing(electric_field &, gas_jet &, int, int, spatial_axes &, temporal_axes &, const int);

/***************************************************************************************************************/
// Functions:

// Calculates plasma frequency for a given free electron density
double Plasma_Frequency_Squared(double ionised, double nd) {
  double wp2 = 4.0*pi*nd*ionised; // plasma frequency squared
  return wp2;
}

// Function to obtain coefficients for 1st order Runge-Kutta evolution
void RungeKutta_FirstOrder(dcmplx *integral, double &exp_factor, dcmplx *E, gas_jet &Jet, double pre_factor, const double dt, const int nz, const int Dim) {
  dcmplx Emag = (Dim==2) ? sqrt(pow(abs(E[0]),2)+pow(abs(E[1]),2)) : E[0];
  exp_factor -= dt * Jet.ADK_Rate(Emag);
  double ionised = 1.0 - exp(exp_factor);
  for(int nd=0; nd<Dim; nd++)
    integral[nd] += pre_factor * Plasma_Frequency_Squared(ionised, Jet.nd[nz]) * E[nd];
}

// Function to obtain coefficients for 1st order Runge-Kutta evolution after half-step advancement
void RungeKutta_FirstOrder(dcmplx *integral, double &exp_factor, dcmplx *E, gas_jet &Jet, double pre_factor, const double dt, const int nz1, const int nz2, const int Dim) {
  dcmplx Emag = (Dim==2) ? sqrt(pow(abs(E[0]),2)+pow(abs(E[1]),2)) : E[0];
  exp_factor -= dt * Jet.ADK_Rate(Emag);
  double ionised = 1.0 - exp(exp_factor);
  for(int nd=0; nd<Dim; nd++)
    integral[nd] += pre_factor * Plasma_Frequency_Squared(ionised, 0.5*(Jet.nd[nz1]+Jet.nd[nz2])) * E[nd];
}

// Function to obtain coefficients for 2nd order Runge-Kutta evolution
void RungeKutta_SecondOrder(dcmplx **integral, double *exp_factor, dcmplx *E, gas_jet &Jet, double pre_factor, const double dz, const double dt, const int nz, const int Dim) {
  dcmplx Estep[2];
  RungeKutta_FirstOrder(integral[0], exp_factor[0], E, Jet, pre_factor, dt, nz, Dim);
  for(int nd=0; nd<Dim; nd++)  Estep[nd] = E[nd]+(0.5*dz)*integral[0][nd];  
  RungeKutta_FirstOrder(integral[1], exp_factor[1], Estep, Jet, pre_factor, dt, nz, nz+1, Dim);
}

// Function to obtain coefficients for 4th order Runge-Kutta evolution
void RungeKutta_FourthOrder(dcmplx **integral, double *exp_factor, dcmplx *E, gas_jet &Jet, double pre_factor, const double dz, const double dt, const int nz, const int Dim) {
  dcmplx Estep[2];
  RungeKutta_FirstOrder(integral[0], exp_factor[0], E, Jet, pre_factor, dt, nz, Dim);
  for(int nd=0; nd<Dim; nd++)  Estep[nd] = E[nd]+(0.5*dz)*integral[0][nd]; 
  RungeKutta_FirstOrder(integral[1], exp_factor[1], Estep, Jet, pre_factor, dt, nz, nz+1, Dim);
  for(int nd=0; nd<Dim; nd++)  Estep[nd] = E[nd]+(0.5*dz)*integral[1][nd]; 
  RungeKutta_FirstOrder(integral[2], exp_factor[2], Estep, Jet, pre_factor, dt, nz, nz+1, Dim);
  for(int nd=0; nd<Dim; nd++)  Estep[nd] = E[nd]+dz*integral[2][nd]; 
  RungeKutta_FirstOrder(integral[3], exp_factor[3], Estep, Jet, pre_factor, dt, nz+1, Dim);
}

// Function to obtain coefficients for 4th order Runge-Kutta-Gill evolution
void Gills_FourthOrder(dcmplx **integral, double *exp_factor, dcmplx *E, gas_jet &Jet, double pre_factor, const double dz, const double dt, const int nz, const int Dim) {
  dcmplx Estep[2];
  RungeKutta_FirstOrder(integral[0], exp_factor[0], E, Jet, pre_factor, dt, nz, Dim);
  for(int nd=0; nd<Dim; nd++)  Estep[nd] = E[nd]+(0.5*dz)*integral[0][nd]; 
  RungeKutta_FirstOrder(integral[1], exp_factor[1], Estep, Jet, pre_factor, dt, nz, nz+1, Dim);
  for(int nd=0; nd<Dim; nd++)  Estep[nd] = E[nd]+0.5*(sqrt(2.0)-1.0)*dz*integral[0][nd]+(1.0-0.5*sqrt(2.0))*dz*integral[1][nd];
  RungeKutta_FirstOrder(integral[2], exp_factor[2], Estep, Jet, pre_factor, dt, nz, nz+1, Dim);
  for(int nd=0; nd<Dim; nd++)  Estep[nd] = E[nd]-0.5*sqrt(2.0)*dz*integral[1][nd] + (1.0+0.5*sqrt(2.0))*dz*integral[2][nd];
  RungeKutta_FirstOrder(integral[3], exp_factor[3], Estep, Jet, pre_factor, dt, nz+1, Dim);
} 

// Function to apply plasma dispersion step due to free-electrons 
void Plasma_Defocussing(electric_field &Efield, gas_jet &Jet, int nz, int order, spatial_axes &X, temporal_axes &T, const int mpi_rank) {
  double pre_factor = -T.dt*0.5/c_au; //-dt*2.0*pi/c_au;
  double exp_factor[order];
  dcmplx E[Efield.Dim];
  dcmplx **integral = new dcmplx*[order];
  for(int no=0; no<order; no++)  *(integral+no) = new dcmplx[Efield.Dim];
  dcmplx *Ex, *Ey;
#ifdef useMPI  
  Ex = Efield.E_mpi;  Ey = Efield.E_mpi2;
#else
  Ex = Efield.Ert;  Ey = Efield.Ert2;
#endif
  if(mpi_rank==0) // this is just for later output to ionised fraction vs propagation distance (not used here; expecting Et as argument so on-axis only)
    Jet.ionised[nz+1] = (Efield.Dim==2) ? Jet.ADK_Ionisation(Ex, Ey) : Jet.ADK_Ionisation(Ex);
  for(int nr=0; nr<X.Nr_mpi; nr++) {
    for(int no=0; no<order; no++) {
      exp_factor[no] = 0.0;
      for(int nd=0; nd<Efield.Dim; nd++)  
	integral[no][nd] = 0.0;
    }
    for(int nt=0; nt<T.Nt; nt++) {
      E[0] = Ex[Index2V(nr,nt,T.Nt)];
      if(Efield.Dim==2)  E[1] = Ey[Index2V(nr,nt,T.Nt)];
      switch(order) {
      case 1: // first-order Runge-Kutta
	RungeKutta_FirstOrder(integral[0], exp_factor[0], E, Jet, pre_factor, T.dt, nz, Efield.Dim);
	Ex[Index2V(nr,nt,T.Nt)] += X.dz*integral[0][0];
	if(Efield.Dim==2)  Ey[Index2V(nr,nt,T.Nt)] += X.dz*integral[0][1];
	break;
      case 2: // second-order Runge-Kutta
	RungeKutta_SecondOrder(integral, exp_factor, E, Jet, pre_factor, X.dz, T.dt, nz, Efield.Dim);
	Ex[Index2V(nr,nt,T.Nt)] += X.dz*integral[1][0];
	if(Efield.Dim==2)  Ey[Index2V(nr,nt,T.Nt)] += X.dz*integral[1][1];
	break;
      case 4: // fourth-order Runge-Kutta
	RungeKutta_FourthOrder(integral, exp_factor, E, Jet, pre_factor, X.dz, T.dt, nz, Efield.Dim);
	Ex[Index2V(nr,nt,T.Nt)] += X.dz * (integral[0][0]/6.0 + integral[1][0]/3.0 + integral[2][0]/3.0 + integral[3][0]/6.0);
	if(Efield.Dim==2)  Ey[Index2V(nr,nt,T.Nt)] += X.dz * (integral[0][1]/6.0 + integral[1][1]/3.0 + integral[2][1]/3.0 + integral[3][1]/6.0);
	break;
      case 5: // fourth-order Runge-Kutta-Gill (call case 5 for sake of ease/laziness)
	Gills_FourthOrder(integral, exp_factor, E, Jet, pre_factor, X.dz, T.dt, nz, Efield.Dim);
	Ex[Index2V(nr,nt,T.Nt)] += X.dz * (integral[0][0]/6.0 + (2.0-sqrt(2.0))*integral[1][0]/6.0 + (2.0+sqrt(2.0))*integral[2][0]/6.0 + integral[3][0]/6.0);
	if(Efield.Dim==2)  Ey[Index2V(nr,nt,T.Nt)] += X.dz * (integral[0][1]/6.0 + (2.0-sqrt(2.0))*integral[1][1]/6.0 + (2.0+sqrt(2.0))*integral[2][1]/6.0 + integral[3][1]/6.0);
	break;
      default:
	std::cerr << "Error in order specified for Runge-Kutta";
#ifdef useMPI
	MPI_Abort(MPI_COMM_WORLD, 1);
	MPI_Finalize();
#endif
	exit(1);
      }
    }
  }
  for(int no=0; no<order; no++)  delete[] *(integral+no);
  delete[] integral;
}

/***************************************************************************************************************/
#endif
