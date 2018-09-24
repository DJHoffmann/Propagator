/***************************************************************************************************************/
// Header file to contain miscellaneous functions for Propagator.cpp 

/***************************************************************************************************************/
// Notes:

/***************************************************************************************************************/
#ifndef MISC_FUNCTIONS_H
#define MISC_FUNCTIONS_H

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
int Index2D(int, int, const int);
int Factorial(int);
double Gamma(double);
template<class T> void Differentiate_Central(T *, T *, double, int);
template<class T> char * Combine_String(char *, const char *, const T);
template<class A, class B> char * Combine_String(char *, const char *, const A, const B);
char Dimension_Label(int);

/***************************************************************************************************************/
// Functions:

// Returns the correct 1D array index for 2D data stored in a 1D array
int Index2V(int nr, int nc, const int Nc) {
  return nr*Nc+nc;
}

// Factorial function
int Factorial(int knd) {
  static double F[180];
  static int sk=1;
  double answer;
  int k;
  F[0]=F[1]=1.0;
  if(knd>=169)
    std::cout << "x too large for x!" << '\n';
  else if(knd<0)
    std::cout << "negative input to x!" << '\n';
  else if(knd<=sk)
    answer=F[knd];
  else {
    answer=F[sk];
    for(k=sk;k<=knd;k++) {
      answer *= k;
      F[k]=answer;
    }
  }
  return answer;
}

// Gamma function
double Gamma(double z) {
  double G = 1.0/z, k;
  for(k = 1; k <= 10000; k++)
    G *= pow(1.0+1.0/(double)k, z) / (1.0+(double)z/k);
  return G;
}

// Simple calulcation of differential using central difference method
template<class T> void Differentiate_Central(T *drE, T *E, double dr, int Nr) {
  drE[0] = 0.0;
  for(int nr=1; nr<Nr-1; nr++)
    drE[nr] = (E[nr+1]-E[nr-1])*0.5/dr;
  drE[Nr-1] = drE[Nr-2];
}

// Function to allow single line filename specification
template<class A> char * Combine_String(char *buffer, const char *name, const A a) {
  sprintf(buffer, name, a);
  return buffer;
}

// Function to allow single line filename specification
template<class A, class B> char * Combine_String(char *buffer, const char *name, const A a, const B b) {
  sprintf(buffer, name, a, b);
  return buffer;
}

// Function to return X,Y,Z label based on input integer
char Dimension_Label(int N) {
  char label;
  switch(N) {
  case 0:
    label = 'X';
    break;
  case 1:
    label = 'Y';
    break;
  case 2:
    label = 'Z';
    break;
  default:
    label = 'E';
    std::cerr << "Error in input for Dimension_Label function \n";
  }
  return label;
}

/***************************************************************************************************************/
#ifdef useMPI

#ifdef oldMPI
#undef SEEK_SET  // undef needed due to a bug with MPI-2, bug fixed in openmpi >= 1.3.2
#undef SEEK_END
#undef SEEK_CUR
#endif
#include<mpi.h>
#include<Electric_Field.h>

/***************************************************************************************************************/
// MPI Function Prototypes:
void MPI_Scatter_Interface(electric_field &, const int, const char, const int, const int);
void MPI_Allgather_Interface(electric_field &, const int, const char, const int, const int);

/***************************************************************************************************************/
// MPI Functions:

// Function to govern MPI_Scatter commands for electric fields
void MPI_Scatter_Interface(electric_field &Efield, const int N, const char domain, const int n1=0, const int n2=0) {
  switch(domain) {
  case 't':
    MPI_Scatter(&Efield.Ert[n1], N, MPI_DOUBLE_COMPLEX, &Efield.E_mpi[n2], N, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    if(Efield.Dim==2)  MPI_Scatter(&Efield.Ert2[n1], N, MPI_DOUBLE_COMPLEX, &Efield.E_mpi2[n2], N, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    break;
  case 'w':
    MPI_Scatter(&Efield.Erw[n1], N, MPI_DOUBLE_COMPLEX, &Efield.E_mpi[n2], N, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    if(Efield.Dim==2)  MPI_Scatter(&Efield.Erw2[n1], N, MPI_DOUBLE_COMPLEX, &Efield.E_mpi2[n2], N, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    break;
  default:
    std::cerr << "Error in domain choice for MPI_Scatter_2D interface\n";
  }
}

// Function to govern MPI_Allgather commands for electric fields
void MPI_Allgather_Interface(electric_field &Efield, const int N, const char domain, const int n1=0, const int n2=0) {
  switch(domain) {
  case 't':
    MPI_Allgather(&Efield.E_mpi[n1], N, MPI_DOUBLE_COMPLEX, &Efield.Ert[n2], N, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
    if(Efield.Dim==2)  MPI_Allgather(&Efield.E_mpi2[n1], N, MPI_DOUBLE_COMPLEX, &Efield.Ert2[n2], N, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
    break;
  case 'w':
    MPI_Allgather(&Efield.E_mpi[n1], N, MPI_DOUBLE_COMPLEX, &Efield.Erw[n2], N, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
    if(Efield.Dim==2)  MPI_Allgather(&Efield.E_mpi2[n1], N, MPI_DOUBLE_COMPLEX, &Efield.Erw2[n2], N, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
    break;
  default:
    std::cerr << "Error in domain choice for MPI_Allgather_2D interface\n";
  }
}

#endif

/***************************************************************************************************************/
#endif
