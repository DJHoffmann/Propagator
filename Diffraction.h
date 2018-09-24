/***************************************************************************************************************/
// Header file to contain diffraction functions based on Crank-Nicolson finite difference for Propagator.cpp 

/***************************************************************************************************************/
// Notes:

/***************************************************************************************************************/
#ifndef DIFFRACTION_H
#define DIFFRACTION_H

#ifdef useMPI
#ifdef oldMPI
#undef SEEK_SET  // undef needed due to a bug with MPI-2, bug fixed in openmpi >= 1.3.2
#undef SEEK_END
#undef SEEK_CUR
#endif
#include<mpi.h>
#endif

/***************************************************************************************************************/
// Function Prototypes:
void Set_DC_to_Zero(dcmplx *, const int);
void CrankNicolson_RHS(dcmplx *, dcmplx *, dcmplx, spatial_axes &);
void CrankNicolson_LHS(dcmplx *, dcmplx *, dcmplx *, dcmplx, spatial_axes &);
void Invert_Tridiagonal(dcmplx *, dcmplx *, dcmplx *, dcmplx *, dcmplx *, const int);
void Diffract_Field(dcmplx *, dcmplx *, dcmplx *, dcmplx *, dcmplx *, spatial_axes &, temporal_axes &, const int);
void Diffract_Field(electric_field &, dcmplx *, dcmplx *, dcmplx *, dcmplx *, spatial_axes &, temporal_axes &, const int);
void Diffract_Field(dcmplx *, dcmplx *, dcmplx *, dcmplx *, dcmplx *, spatial_axes &, temporal_axes &, const double, const int);
void Diffract_Field(electric_field &, dcmplx *, dcmplx *, dcmplx *, dcmplx *, spatial_axes &, temporal_axes &, const double, const int);

/***************************************************************************************************************/
// Functions:

// Sets DC to component to zero (must be physical and would otherwise diverge to infinity by the following calculation)
void Set_DC_to_Zero(dcmplx *Erw, const int Nr) {  // only ever call on mpi_rank holding the zero frequency!
  for(int nr=0; nr<Nr; nr++)
    Erw[nr] = 0.0; // since only the first row don't need Index2V(0,nr,Nr)
}

// Calculates D solution vector of C-N (RHS; length Ntp)
void CrankNicolson_RHS(dcmplx *D, dcmplx *Er, dcmplx alpha, spatial_axes &X) {
  double beta;
  D[0] = (1.0+2.0*alpha)*Er[0] - 2.0*alpha*Er[1]; // counts +1 element for missing -1 element (symmetric about 0)
  for(int nr=1; nr<X.Nr-1; nr++) {
    beta = X.dr/(2.0*X.r[nr]);
    D[nr] = -alpha*(1.0-beta)*Er[nr-1]*X.AB[nr-1] + (1.0+2.0*alpha)*Er[nr]*X.AB[nr] - alpha*(1.0+beta)*Er[nr+1]*X.AB[nr+1]; //check Vr part
  }
  D[X.Nr-1] = 0.0; //D[X.Nr-2];  // is this correct???
}

// Calculates A,B,C coefficient vectors of C-N (LHS; length Ntp)
void CrankNicolson_LHS(dcmplx *A, dcmplx *B, dcmplx *C, dcmplx alpha, spatial_axes &X) {
  double beta;
  A[0] = 0.0;
  B[0] = 1.0 - 2.0*alpha;
  C[0] = 2.0*alpha; // counts +1 element for missing -1 element (symmetric about 0)
  for(int nr=1; nr<X.Nr; nr++) {
    beta = X.dr/(2.0*X.r[nr]);
    A[nr] = alpha*(1.0-beta);
    B[nr] = 1.0-2.0*alpha;
    C[nr] = alpha*(1.0+beta);
  }
  C[X.Nr-1] = 0.0;
}

// Inverts the tridiagonal matrix to find z+dz field (see Numerical Recipes in C++)
void Invert_Tridiagonal(dcmplx *Er, dcmplx *A, dcmplx *B, dcmplx *C, dcmplx *D, const int Nr) {
  C[0] /= B[0];
  D[0] /= B[0];
  for(int nr=1; nr<Nr; nr++) {
    C[nr] /= (B[nr]-C[nr-1]*A[nr]);
    D[nr] = (D[nr]-D[nr-1]*A[nr])/(B[nr]-C[nr-1]*A[nr]);
  }
  Er[Nr-1] = D[Nr-1];
  for(int nr=Nr-2; nr>=0; nr--)
    Er[nr] = D[nr]-C[nr]*Er[nr+1];
}

// Diffract field using Crank-Nicholson scheme
void Diffract_Field(dcmplx *Erw, dcmplx *A, dcmplx *B, dcmplx *C, dcmplx *D, spatial_axes &X, temporal_axes &T, const int mpi_rank) {
  dcmplx prefac = i*c_au/4.0*X.dz/(X.dr*X.dr);
  for(int nw=0; nw<T.Nt_mpi; nw++) {
    CrankNicolson_RHS(D, &Erw[Index2V(nw,0,X.Nr)], prefac/T.w[mpi_rank*T.Nt_mpi+nw], X);
    CrankNicolson_LHS(A, B, C, prefac/T.w[mpi_rank*T.Nt_mpi+nw], X);
    Invert_Tridiagonal(&Erw[Index2V(nw,0,X.Nr)], A, B, C, D, X.Nr);
  }
  if(mpi_rank==0)
    Set_DC_to_Zero(Erw, X.Nr);
}

// As above taking electric_field class as argument
void Diffract_Field(electric_field &Efield, dcmplx *A, dcmplx *B, dcmplx *C, dcmplx *D, spatial_axes &X, temporal_axes &T, const int mpi_rank) {
#ifdef useMPI  
  Diffract_Field(Efield.E_mpi, A, B, C, D, X, T, mpi_rank);
  if(Efield.Dim==2)  Diffract_Field(Efield.E_mpi2, A, B, C, D, X, T, mpi_rank);
#else
  Diffract_Field(Efield.Erw, A, B, C, D, X, T, mpi_rank);
  if(Efield.Dim==2)  Diffract_Field(Efield.Erw2, A, B, C, D, X, T, mpi_rank);
#endif
}

// Diffract field using Crank-Nicholson scheme but taking in dz as an argument
void Diffract_Field(dcmplx *Erw, dcmplx *A, dcmplx *B, dcmplx *C, dcmplx *D, spatial_axes &X, temporal_axes &T, const double dz, const int mpi_rank) {
  dcmplx prefac = i*c_au/4.0*dz/(X.dr*X.dr);
  for(int nw=0; nw<T.Nt_mpi; nw++) {
    CrankNicolson_RHS(D, &Erw[Index2V(nw,0,X.Nr)], prefac/T.w[mpi_rank*T.Nt_mpi+nw], X);
    CrankNicolson_LHS(A, B, C, prefac/T.w[mpi_rank*T.Nt_mpi+nw], X);
    Invert_Tridiagonal(&Erw[Index2V(nw,0,X.Nr)], A, B, C, D, X.Nr);
  }
  if(mpi_rank==0)
    Set_DC_to_Zero(Erw, X.Nr);
}

// As above taking electric_field class as argument
void Diffract_Field(electric_field &Efield, dcmplx *A, dcmplx *B, dcmplx *C, dcmplx *D, spatial_axes &X, temporal_axes &T, const double dz, const int mpi_rank) {
#ifdef useMPI  
  Diffract_Field(Efield.E_mpi, A, B, C, D, X, T, dz, mpi_rank);
  if(Efield.Dim==2)  Diffract_Field(Efield.E_mpi2, A, B, C, D, X, T, dz, mpi_rank);
#else
  Diffract_Field(Efield.Erw, A, B, C, D, X, T, dz, mpi_rank);
  if(Efield.Dim==2)  Diffract_Field(Efield.Erw2, A, B, C, D, X, T, dz, mpi_rank);
#endif
}

/***************************************************************************************************************/
#endif
