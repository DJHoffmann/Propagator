/***************************************************************************************************************/
// Header file to contain universal output functions for Propagator.cpp 

/***************************************************************************************************************/
// Notes:

/***************************************************************************************************************/
#ifndef OUTPUT_H
#define OUTPUT_H

#include<fstream>
#include<Environment.h>
#include<Laser_Parameters.h>
#include<Electric_Field.h>
#include<Gas_Jet.h>

/***************************************************************************************************************/
// Function Prototypes:
void Output_Rank_Selector(int *, const int, const int);
void Output_Spec(electric_field *, laser_parameters **, spatial_axes &, temporal_axes &, const int, const int, std::ofstream&);
void Output_Spec(electric_field *, laser_parameters **, spatial_nearfar_axes &, temporal_axes &, const int, const int, const int, std::ofstream&);
void Get_Frequency_Components(double *&, int *&, int &, gas_jet &, laser_parameters **, temporal_axes &T, std::ifstream &);

/***************************************************************************************************************/
// Functions:

// Function to correctly distribute output functions across mpi_ranks
void Output_Rank_Selector(int *rank_choice, const int streams, const int mpi_size) {
  for(int ns=0; ns<streams; ns++)
    //rank_choice[ns] = ns%mpi_size;
    rank_choice[ns] = ns%std::min(mpi_size,8);
}

// Output spec
void Output_Spec(electric_field *Efield, laser_parameters **pLP, spatial_axes &X, temporal_axes &T, const int mpi_rank, const int mpi_size, std::ofstream &OUTSTREAM) {
  OUTSTREAM << pLP[0]->E0 << '\t' << pLP[0]->w0 << '\t' << pLP[0]->cycles << '\t' << pLP[0]->cep << '\t' << pLP[0]->wz0 << '\t' << Efield->Dim << '\n';
  OUTSTREAM << X.dr << '\t' << X.dz << '\t' << T.dt << '\t' << X.Nr << '\t' << X.Nz << '\t' << T.Nt << '\n'; 
  OUTSTREAM << mpi_rank << '\t' << mpi_size << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\n'; 
}

// Output spec
void Output_Spec(electric_field *Efield, laser_parameters **pLP, spatial_nearfar_axes &X, temporal_axes &T, const int mpi_rank, const int mpi_size, const int filter, std::ofstream &OUTSTREAM) {
  OUTSTREAM << pLP[0]->E0 << '\t' << pLP[0]->w0 << '\t' << pLP[0]->cycles << '\t' << pLP[0]->cep << '\t' << pLP[0]->wz0 << '\t' << Efield->Dim << '\n';
  OUTSTREAM << X.dr << '\t' << X.dz << '\t' << T.dt << '\t' << X.Nr << '\t' << X.Nz << '\t' << T.Nt << '\n'; 
  OUTSTREAM << mpi_rank << '\t' << mpi_size << '\t' << X.L << '\t' << X.drp << '\t' << filter << '\t' << 0 << '\n'; 
}

// Calculates the frequency components (as indices of frequncy vector) for outputting spatially through interaction region
void Get_Frequency_Components(double *&FreqCmpts, int *&FreqInds, int &Nfreqs, gas_jet &Jet, laser_parameters **pLP, temporal_axes &T, std::ifstream &INSTREAM) {
  char DoThis = Read_Variable<char>(INSTREAM, "OutFreqCmpts");
  if(DoThis=='y') {
    Nfreqs = Read_Variable<int>(INSTREAM, "NFreqs");
    FreqCmpts = new double[Nfreqs];  FreqInds = new int[Nfreqs];
    int FreqRef = Read_Variable<int>(INSTREAM, "FreqRef");
    double wRef = (FreqRef>0) ? 1.32*Jet.Ip + 3.17*pLP[FreqRef-1]->Ponderomotive_Potential() : T.w[T.Nt/2];
    double wMin = Read_Variable<double>(INSTREAM, "minFreq") * wRef;
    double wMax = Read_Variable<double>(INSTREAM, "maxFreq") * wRef;
    if(wMin>=T.w[T.Nt/2]) {
      wMin = T.w[T.Nt/2]/Nfreqs;
      wMax = T.w[T.Nt/2-1];
      std::cout << "Both minFreq & maxFreq (wrt FreqRef) set too high for defined frequency axis, resetting to " << 1.0/Nfreqs << "* and 1*max defined frequency, respectively\n";
    }
    else if(wMax>=T.w[T.Nt/2]) {
      wMax = T.w[T.Nt/2-1];
      std::cout << "maxFreq (wrt FreqRef) set too high for defined frequency axis, resetting to max defined frequency\n";
    }
    for(int nfreqs=0; nfreqs<Nfreqs; nfreqs++) 
      FreqCmpts[nfreqs] = wMin+(double)nfreqs/(Nfreqs-1)*(wMax-wMin);
    char FreqShift = Read_Variable<char>(INSTREAM, "FreqShift");
    double wHO;
    switch(FreqShift) {  // shifts energies according to following options wrt colour specified by FreqRef (or 1st if FreqRef=0) 
    case 'h':  // nearest harmonic
      for(int nfreqs=0; nfreqs<Nfreqs; nfreqs++) {
	wHO = FreqCmpts[nfreqs]/pLP[FreqRef-1]->w0;
	FreqCmpts[nfreqs] = (wHO-(int)wHO>0.5) ? ((int)wHO+1)*pLP[std::max(0,FreqRef-1)]->w0 : ((int)wHO)*pLP[std::max(0,FreqRef-1)]->w0;
      } 
      break;
    case 'o':  // nearest odd harmonic
      for(int nfreqs=0; nfreqs<Nfreqs; nfreqs++) {
	wHO = FreqCmpts[nfreqs]/(2.0*pLP[std::max(0,FreqRef-1)]->w0);
	FreqCmpts[nfreqs] = (2*(int)wHO+1)*pLP[std::max(0,FreqRef-1)]->w0; 
      } 
      break;
    case 'e':  // nearest even harmonic
      for(int nfreqs=0; nfreqs<Nfreqs; nfreqs++) {
	wHO = FreqCmpts[nfreqs]/(2.0*pLP[std::max(0,FreqRef-1)]->w0);
	FreqCmpts[nfreqs] = (wHO-(int)wHO>0.5) ? ((int)wHO+1)*(2.0*pLP[std::max(0,FreqRef-1)]->w0) : ((int)wHO)*(2.0*pLP[std::max(0,FreqRef-1)]->w0);
      } 
      break;
    default:
      break;
    }
    for(int nfreqs=0; nfreqs<Nfreqs; nfreqs++) {
      wHO = FreqCmpts[nfreqs]/T.dw;
      FreqInds[nfreqs] = (wHO-(int)wHO>0.5) ? 1+(int)wHO : (int)wHO; 
    }
  }
  else
    Nfreqs = 0;
}

/***************************************************************************************************************/
#endif
