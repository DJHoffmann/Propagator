/***************************************************************************************************************/
// Header file to contain Hankel transform to farfield for use in Propagator.cpp 

/***************************************************************************************************************/
// Notes:

/***************************************************************************************************************/
#ifndef FAR_FIELD_H
#define FAR_FIELD_H

#include<iostream>
#include<complex>
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

/***************************************************************************************************************/
// Function Prototypes:
void Hankel_Propagation(electric_field &, spatial_nearfar_axes &, temporal_axes &, const char, const int);
void Hankel_Abel(dcmplx *, dcmplx *, spatial_nearfar_axes &, temporal_axes &, fftw_complex *, fftw_plan, const int);
void Hankel_Abel(dcmplx *, dcmplx *, spatial_nearfar_axes &, temporal_axes &, const int);
dcmplx Abel_Integral_by_Parts(dcmplx *, dcmplx *, double, spatial_axes &);
void Hankel_Direct(dcmplx *, dcmplx *, spatial_nearfar_axes &, temporal_axes &, const int);
double BESSJ0 (double);

/***************************************************************************************************************/
// Functions:

// Governing interface for Hankel transform
void Hankel_Propagation(electric_field &Hfield, spatial_nearfar_axes &X, temporal_axes &T, const char FPmethod, const int mpi_rank) {
  switch(FPmethod) {
  case 'A':
#ifdef useMPI 
    Hankel_Abel(Hfield.E_mpi, Hfield.E_mpi, X, T, mpi_rank);
    if(Hfield.Dim==2)  Hankel_Abel(Hfield.E_mpi2, Hfield.E_mpi2, X, T, mpi_rank);
#else
    Hankel_Abel(Hfield.Erw, Hfield.Erw, X, T, mpi_rank);
    if(Hfield.Dim==2)  Hankel_Abel(Hfield.Erw2, Hfield.Erw2, X, T, mpi_rank);
#endif
    break;
  case 'D':
#ifdef useMPI 
    Hankel_Direct(Hfield.E_mpi, Hfield.E_mpi, X, T, mpi_rank);
    if(Hfield.Dim==2)  Hankel_Direct(Hfield.E_mpi2, Hfield.E_mpi2, X, T, mpi_rank);
#else
    Hankel_Direct(Hfield.Erw, Hfield.Erw, X, T, mpi_rank);
    if(Hfield.Dim==2)  Hankel_Direct(Hfield.Erw2, Hfield.Erw2, X, T, mpi_rank);
#endif
    break;
  default:
    std::cerr << "Error in choice of method for far-field propagation\n";
#ifdef useMPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(1);
  }
}

// Function to perform Hankel integration looped over multiple frequency components using the Abel integral method
void Hankel_Abel(dcmplx *EF, dcmplx *EN, spatial_nearfar_axes &X, temporal_axes &T, fftw_complex *fftwArrayY, fftw_plan y2rp, const int mpi_rank) {
  const int Ny = 2*X.Nr; // forms axis by reflecting rp through zero (Nrp=Nr) - synchronises with length of fftwArrayY and y2rp (don't change!)
  const double dy = pi/(X.drp*X.Nr);
  double *y = new double[Ny];
  for(int ny=0; ny<Ny; ny++)
    y[ny] = (ny-Ny/2)*dy;
  double Bessfac, rp, fftw_factor = dy/sqrt(2.0*pi), PoM;
  dcmplx Expfac, tempC;
  dcmplx *F = new dcmplx[X.Nr], *drF = new dcmplx[X.Nr];
  const int nw0 = (mpi_rank==0) ? 1 : 0;  // to prevent division by zero below
  for(int nw=nw0; nw<T.Nt_mpi; nw++) {
    Bessfac = T.w[mpi_rank*T.Nt_mpi+nw]/(c_au*X.L);
    Expfac = -0.5*i*Bessfac;
    for(int nr=0; nr<X.Nr; nr++)
      F[nr] = EN[Index2V(nw,nr,X.Nr)]*exp(Expfac*X.r[nr]*X.r[nr]);
    Differentiate_Central(drF, F, X.dr, X.Nr);
    for(int ny=0; ny<Ny; ny++) {
      tempC = Abel_Integral_by_Parts(F, drF, y[ny]/fabs(Bessfac), X);
      fftwArrayY[ny][0] = real(tempC);
      fftwArrayY[ny][1] = imag(tempC);
    }
    fftw_execute(y2rp);
    for(int nrp=0; nrp<X.Nr; nrp++) {  // selects positive values of y
      rp = X.r[nrp]*X.drp/X.dr;
      PoM = (T.w[mpi_rank*T.Nt_mpi+nw]<0.0) ? -1.0 : 1.0;
      // EF[Index2V(nw,nrp,X.Nr)] = 1.0/sqrt(2.0*pi) * (fftw_factor/Bessfac)*(fftwArrayY[nrp][0]+i*fftwArrayY[nrp][1]) * (PoM*i)*Bessfac * exp(Expfac*rp*rp);
      // --------> simplifies as below --------->
      EF[Index2V(nw,nrp,X.Nr)] = 1.0/sqrt(2.0*pi) * fftw_factor*(fftwArrayY[nrp][0]+i*fftwArrayY[nrp][1]) * (PoM*i) * exp(Expfac*rp*rp);
    }
  }
  if(mpi_rank==0)  // DC component
    for(int nrp=0; nrp<X.Nr; nrp++)
      EF[Index2V(0,nrp,X.Nr)] = 0.0;
  delete[] y, F, drF;
}

// As above but creating FFTW arra & plan within function
void Hankel_Abel(dcmplx *EF, dcmplx *EN, spatial_nearfar_axes &X, temporal_axes &T, const int mpi_rank) {
  const int Ny = 2*X.Nr; // forms axis by reflecting rp through zero (Nrp=Nr) [repeated within function above]
  fftw_complex *fftwArrayY;
  fftw_plan plan_y2rp;
#ifdef useMPI
  int size = 1;  // force sequential for thread-safety
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  for(int n=0; n<size; n++) {
    if(mpi_rank==n) {
      fftwArrayY = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*Ny);
      plan_y2rp = fftw_plan_dft_1d(Ny, fftwArrayY, fftwArrayY, FFTW_FORWARD, FFTW_MEASURE);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#else
  fftwArrayY = (fftw_complex *)fftw_malloc(sizeof(fftw_complex)*Ny);
  plan_y2rp = fftw_plan_dft_1d(Ny, fftwArrayY, fftwArrayY, FFTW_FORWARD, FFTW_MEASURE);
#endif
  Hankel_Abel(EF, EN, X, T, fftwArrayY, plan_y2rp, mpi_rank);
#ifdef useMPI
  for(int n=0; n<size; n++) {
    if(mpi_rank==n) {
      fftw_free(fftwArrayY);
      fftw_destroy_plan(plan_y2rp);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#else
  fftw_free(fftwArrayY);
  fftw_destroy_plan(plan_y2rp);
#endif
}

// Performs Abel integral for single y value by parts (to eliminate singularity)
dcmplx Abel_Integral_by_Parts(dcmplx *EN, dcmplx *drEN, double y, spatial_axes &X) {
  double y2 = y*y;
  dcmplx EF = 0.0;
  int nr0 = (int)ceil(fabs(y)/X.dr); //std::min((int)(fabs(y)/X.dr),X.Nr-1);
  for(int nr=nr0; nr<X.Nr; nr++) 
    EF -= sqrt(X.r[nr]*X.r[nr]-y2)*drEN[nr];
  EF *= X.dr;
  EF += sqrt(X.r[X.Nr-1]*X.r[X.Nr-1]-y2)*EN[X.Nr-1];
  EF *= 2.0; 
  if(fabs(y)>=(X.Nr-1)*X.dr) EF = 0.0;
  return EF;
}

/* // Function to perform Hankel integration looped over multiple frequency components */
/* void Hankel_Direct(dcmplx *EF, dcmplx *EN, spatial_nearfar_axes &X, temporal_axes &T, const int mpi_rank) { */
/*   double Bessfac, r, rp, PoM; */
/*   dcmplx Expfac; */
/*   dcmplx *F = new dcmplx[X.Nr]; */
/*   for(int nw=0; nw<T.Nt_mpi; nw++) { */
/*     Bessfac = T.w[mpi_rank*T.Nt_mpi+nw]/(c_au*X.L); */
/*     Expfac = -0.5*i*Bessfac; */
/*     for(int nrp=0; nrp<X.Nr; nrp++) { */
/*       rp = X.r[nrp]*X.drp/X.dr; */
/*       F[nrp] = 0.0; */
/*       for(int nr=0; nr<X.Nr; nr++) { */
/* 	r = X.r[nr]; */
/* 	F[nrp] += EN[Index2V(nw,nr,X.Nr)]*exp(Expfac*r*r)*BESSJ0(fabs(Bessfac)*r*rp)*r; */
/*       } */
/*       PoM = (T.w[mpi_rank*T.Nt_mpi+nw]<0.0) ? -1.0 : 1.0; */
/*       F[nrp] *= X.dr * (PoM*i)*fabs(Bessfac)*exp(Expfac*rp*rp); */

/*     } */
/*     for(int nrp=0; nrp<X.Nr; nrp++)   */
/*       EF[Index2V(nw,nrp,X.Nr)] = F[nrp]; */
/*   } */
/* } */

// Function to perform Hankel integration looped over multiple frequency components
void Hankel_Direct(dcmplx *EF, dcmplx *EN, spatial_nearfar_axes &X, temporal_axes &T, const int mpi_rank) {
  double Bessfac, r, rp, PoM;
  dcmplx Expfac;
  dcmplx *F = new dcmplx[X.Nr];
  for(int nw=0; nw<T.Nt_mpi; nw++) {
    Bessfac = T.w[mpi_rank*T.Nt_mpi+nw]/(c_au*X.L);
    Expfac = -0.5*i*Bessfac;
    for(int nrp=0; nrp<X.Nr; nrp++) {
      rp = X.r[nrp]*X.drp/X.dr;
      F[nrp] = 0.0;
      for(int nr=0; nr<X.Nr; nr++) {
	r = X.r[nr];
	F[nrp] += EN[Index2V(nw,nr,X.Nr)]*exp(Expfac*r*r)*BESSJ0(Bessfac*r*rp)*r;
      }
      F[nrp] *= X.dr * i*Bessfac*exp(Expfac*rp*rp);

    }
    for(int nrp=0; nrp<X.Nr; nrp++)  
      EF[Index2V(nw,nrp,X.Nr)] = F[nrp];
  }
}

// see http://jean-pierre.moreau.pagesperso-orange.fr/c_bessel.html
double BESSJ0 (double X) {
/***********************************************************************
      This subroutine calculates the First Kind Bessel Function of
      order 0, for any real number X. The polynomial approximation by
      series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
      REFERENCES:
      M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
      C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
      VOL.5, 1962.
************************************************************************/
  const double
    P1=1.0, P2=-0.1098628627E-2, P3=0.2734510407E-4,
    P4=-0.2073370639E-5, P5= 0.2093887211E-6,
    Q1=-0.1562499995E-1, Q2= 0.1430488765E-3, Q3=-0.6911147651E-5,
    Q4= 0.7621095161E-6, Q5=-0.9349451520E-7,
    R1= 57568490574.0, R2=-13362590354.0, R3=651619640.7,
    R4=-11214424.18, R5= 77392.33017, R6=-184.9052456,
    S1= 57568490411.0, S2=1029532985.0, S3=9494680.718,
    S4= 59272.64853, S5=267.8532712, S6=1.0;
  double AX,FR,FS,Z,FP,FQ,XX,Y, TMP;
  if (X==0.0) return 1.0;
  AX = fabs(X);
  if (AX < 8.0) {
    Y = X*X;
    FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))));
    FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))));
    TMP = FR/FS;
  }
  else {
    Z = 8./AX;
    Y = Z*Z;
    XX = AX-0.785398164;
    FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)));
    FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)));
    TMP = sqrt(0.636619772/AX)*(FP*cos(XX)-Z*FQ*sin(XX));
  }
  return TMP;
}

/***************************************************************************************************************/
#endif
