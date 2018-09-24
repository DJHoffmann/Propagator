/***************************************************************************************************************/
// Propagator.cpp

// Code to proagated an intense, multicolour laser pulse through a gas volume in cylindrically symmetric coordinates.
// The generated harmonic field is provided via the nonlinear dipole response calculated by SFA or Quantum Orbit models.

// David Hoffmann
// Imperial College 
// 2013

// Use of this code and the header files listed below may be granted upon request.
// Any contribution towards published work should be properly acknowledged.

/***************************************************************************************************************/
// Libraries:

#include<iostream>
#include<fstream>
#include<complex>
#include<time.h>
#include<fftw3.h>

#ifdef useMPI
#ifdef oldMPI
#undef SEEK_SET  // undef needed due to a bug with MPI-2, bug fixed in openmpi >= 1.3.2
#undef SEEK_END
#undef SEEK_CUR
#endif
#include<mpi.h>
#endif

using namespace std;

/***************************************************************************************************************/
// Header Files:
#include<Constants_Conversions.h>
#include<Misc_Functions.h>
#include<Input.h>
#include<Output.h>
#include<Environment.h>
#include<Gas_Jet.h>
#include<Laser_Parameters.h>
#include<Electric_Field.h>
#include<Diffraction.h>
#include<Plasma_Defocussing.h>
#include<Dipole_Response.h>
#include<Far_Field.h>


/***************************************************************************************************************/
// Main function:
int main(int argc, char **argv) {

  /***************************************************************************************************************/
  // MPI parameters & initialisation:
  int mpi_size = 1;  // set defaults
  int mpi_rank = 0;  // have these 2 outside useMPI as can use to direct program as single/parallel
  int mpi_sum = 1;  // equal to mpi_rank plus mpi_size
#ifdef useMPI  
  int mpi_err;
  mpi_err = MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  mpi_sum = mpi_rank + mpi_size;
  MPI_Status mpistatus;
  const bool mpi_active(true);
#else
  const bool mpi_active(false);
#endif

  /***************************************************************************************************************/
  //Clock variables & start:
  time_t global_start, global_finish, local_start, local_finish;
  double elapsed;
  time(&global_start);

  /***************************************************************************************************************/
  // Input path and stream
  char ConfigDefault[100] = "./Config/", *ConfigPath = ConfigDefault;
  Read_Command_Options(ConfigPath, ConfigDefault, argc, argv);
  string ControlPath = ConfigPath;  ControlPath.append("ControlBox.conf");
  ifstream CONTROL(ControlPath.c_str());
  
  /***************************************************************************************************************/
  // Program control parameters:
  const char Propagate = Read_Variable<char>(CONTROL, "Propagate");       // perform propagation?
  const char TransLaser = Read_Variable<char>(CONTROL, "WhereLaserDef");  // specifies where along propagation axis the laser parameters apply
  const int order = Read_Variable<int>(CONTROL, "RKOrder");               // order for Runge-Kutta advancement in free-electron dispersion step
  const double minI = Read_Variable<double>(CONTROL, "minI0");            // specifies min (temporal-peak) intensity for which dipole response is calculated
  const char FPmethod = Read_Variable<char>(CONTROL, "FPMethod");         // method for calculating far-field propagation
  const int filter = Read_Variable<int>(CONTROL, "OutFilter");            // specifies density of time/freq-axis for output through z (n outputs every nth point)
  const char PulseTime = Read_Variable<char>(CONTROL, "PulseTime");       // dictates whether driving pulse is outputted in time domain running through z
  const char PulseFreq = Read_Variable<char>(CONTROL, "PulseFreq");       // dictates whether driving pulse is outputted in freq domain running through z
  const char HarmFreq = Read_Variable<char>(CONTROL, "HarmFreq");         // dictates whether harmonics are outputted in freq domain running through z

  /***************************************************************************************************************/
  // Define space & time:
  spatial_nearfar_axes X(CONTROL, mpi_size);
  temporal_axes T(CONTROL, mpi_size);

  /***************************************************************************************************************/
  // Construct gas jet:
  gas_jet Jet(CONTROL, X, T, ConfigPath);

  /***************************************************************************************************************/
  // Construct laser field:
  int colours;
  laser_parameters **LP = NULL;
  Import_Parameters(LP, colours, CONTROL);
  electric_field *Laser = NULL;
  if(TransLaser=='j'||TransLaser=='i') {
    double zDef = (TransLaser=='j') ? Jet.zj : X.z0;
    Laser = (Set_Dimension(LP, colours)==true) ? new electric_field_2d(T, X, LP, zDef, colours, mpi_active) : new electric_field_1d(T, X, LP, zDef, colours, mpi_active);
    //Laser = new electric_field_1d(T, X, LP, zDef, colours, mpi_active);
  }
  else  Laser = (Set_Dimension(LP, colours)==true) ? new electric_field_2d(T, X, LP, colours, mpi_active) : new electric_field_1d(T, X, LP, colours, mpi_active);
  if(mpi_rank==0) {
    cout << "\nFound " << colours << " colour(s) ... creating " << Laser->Dim << "-dimensional laser field\n";
    cout << "Pulse energy at z = " << X.z0*au_to_m*1e3 << "mm is calculated as E = " << Laser->Pulse_Energy('n')*1.0e3 << " mJ\n";
  }

  /***************************************************************************************************************/
  if(mpi_rank==0) {
    cout << "colours = " << colours << '\n';
    for(int nc=0; nc<colours; nc++)
      cout << "I0 = " << pow(LP[nc]->E0,2)*au_to_Wcm2*1e-14 << "e14Wcm^-2, lambda = " << LP[nc]->lambda*au_to_m*1e9 << "nm, t0 = " << LP[nc]->t0*au_to_s*1e15 << "fs, fwhm = " << LP[nc]->fwhm*au_to_s*1e15 << "fs, cep = " << LP[nc]->cep/pi << "pi rad\n";
  }

  laser_parameters **LP_analysed = new laser_parameters*[Laser->colours];
  Laser->Analyse_Field(LP_analysed, 0, false, false);
  if(mpi_rank==0) {
    for(int nc=0; nc<colours; nc++)
      cout << "I0 = " << pow(LP_analysed[nc]->E0,2)*au_to_Wcm2*1e-14 << "e14Wcm^-2, lambda = " << LP_analysed[nc]->lambda*au_to_m*1e9 << "nm, t0 = " << LP_analysed[nc]->t0*au_to_s*1e15 << "fs, fwhm = " << LP_analysed[nc]->fwhm*au_to_s*1e15 << "fs, cep = " << LP_analysed[nc]->cep/pi << "pi rad\n";
  }

  Laser->Analyse_Field(LP_analysed, 0, true, false);
  if(mpi_rank==0) {
    for(int nc=0; nc<colours; nc++)
      cout << "I0 = " << pow(LP_analysed[nc]->E0,2)*au_to_Wcm2*1e-14 << "e14Wcm^-2, lambda = " << LP_analysed[nc]->lambda*au_to_m*1e9 << "nm, t0 = " << LP_analysed[nc]->t0*au_to_s*1e15 << "fs, fwhm = " << LP_analysed[nc]->fwhm*au_to_s*1e15 << "fs, cep = " << LP_analysed[nc]->cep/pi << "pi rad\n";
  }

  /***************************************************************************************************************/
  // Construct harmonic field and single-atom dipole response objects:
  electric_field *Harmonics = NULL;
  Harmonics = (Laser->Dim==2) ? new electric_field_2d(T, X, mpi_active) : new electric_field_1d(T, X, mpi_active);
  //Harmonics = new electric_field_1d(T, X, mpi_active);
  dipole_response *Dipole = NULL;
  Assign_Method(Dipole, *Laser, Jet, T, CONTROL);  

  /***************************************************************************************************************/
  // Define remaining variables:
  dcmplx *A = new dcmplx[X.Nr];
  dcmplx *B = new dcmplx[X.Nr];
  dcmplx *C = new dcmplx[X.Nr];
  dcmplx *D = new dcmplx[X.Nr];

  /***************************************************************************************************************/
  // Determine which frequency components to output through interaction region:
  int Nfreq = 0;
  int *FreqInds = NULL;
  double *FreqCmpts = NULL;
  Get_Frequency_Components(FreqCmpts, FreqInds, Nfreq, Jet, LP, T, CONTROL);

  /***************************************************************************************************************/
  // Output filestreams:
  ofstream SPEC("Spec.dat");
  ofstream AXISR("AxisR.dat");
  ofstream AXISZ("AxisZ.dat");
  ofstream AXISTP("AxisTp.dat");
  ofstream AXISW("AxisW.dat");
  ofstream AXISRP("AxisRp.dat");
  ofstream GASJET("GasJet.dat");
  ofstream IONISED("IonisedOA.dat");
  ofstream *PULSE_W_NORM = new ofstream[Laser->Dim], *PULSE_W_PHASE = new ofstream[Laser->Dim];
  ofstream *PULSE_TP = new ofstream[Laser->Dim];
  ofstream *PULSE_RZ = new ofstream[Laser->Dim];
  ofstream *PULSE_RW_NORM = new ofstream[Laser->Dim], *PULSE_RW_PHASE = new ofstream[Laser->Dim];
  ofstream *PULSE_RTP = new ofstream[Laser->Dim];
  ofstream *PULSE_RZW_NORM = new ofstream[mpi_size*Laser->Dim], *PULSE_RZW_PHASE = new ofstream[mpi_size*Laser->Dim]; 
  ofstream *PULSE_RZTP = new ofstream[mpi_size*Laser->Dim];
  ofstream *PULSE_ZTP = new ofstream[Laser->Dim];
  ofstream *HARM_NEAR_RW_NORM = new ofstream[Harmonics->Dim], *HARM_NEAR_RW_PHASE = new ofstream[Harmonics->Dim];
  ofstream *HARM_NEAR_RTP = new ofstream[Harmonics->Dim];
  ofstream *HARM_NEAR_RZW_NORM = new ofstream[mpi_size*Harmonics->Dim], *HARM_NEAR_RZW_PHASE = new ofstream[mpi_size*Harmonics->Dim]; 
  ofstream *HARM_FAR_RW_NORM = new ofstream[Harmonics->Dim], *HARM_FAR_RW_PHASE = new ofstream[Harmonics->Dim];
  ofstream *HARM_FAR_RTP = new ofstream[Harmonics->Dim];
  ofstream *HARM_NEAR_RZ_NORM = new ofstream[Nfreq*Harmonics->Dim], *HARM_NEAR_RZ_PHASE = new ofstream[Nfreq*Harmonics->Dim];
  char str[100];
  for(int nd=0; nd<Laser->Dim; nd++) {
    PULSE_W_NORM[nd].open(Combine_String(str, "PulseWnorm%c.dat", Dimension_Label(nd)));
    PULSE_W_PHASE[nd].open(Combine_String(str, "PulseWphase%c.dat", Dimension_Label(nd)));
    PULSE_TP[nd].open(Combine_String(str, "PulseTp%c.dat", Dimension_Label(nd)));    
    PULSE_RZ[nd].open(Combine_String(str, "PulseRZ%c.dat", Dimension_Label(nd)));
    PULSE_RW_NORM[nd].open(Combine_String(str, "PulseRWnorm%c.dat", Dimension_Label(nd)));
    PULSE_RW_PHASE[nd].open(Combine_String(str, "PulseRWphase%c.dat", Dimension_Label(nd)));
    PULSE_RTP[nd].open(Combine_String(str, "PulseRTp%c.dat", Dimension_Label(nd)));
    if(PulseFreq=='y') {
      PULSE_RZW_NORM[mpi_rank+nd*mpi_size].open(Combine_String(str, "PulseRZWnorm%c_C%d.dat", Dimension_Label(nd), mpi_rank));
      PULSE_RZW_PHASE[mpi_rank+nd*mpi_size].open(Combine_String(str, "PulseRZWphase%c_C%d.dat", Dimension_Label(nd), mpi_rank));
    }
    if(PulseTime=='y')  PULSE_RZTP[mpi_rank+nd*mpi_size].open(Combine_String(str, "PulseRZTp%c_C%d.dat", Dimension_Label(nd), mpi_rank)); 
    PULSE_ZTP[nd].open(Combine_String(str, "PulseZTp%c.dat", Dimension_Label(nd)));  
  }
  for(int nd=0; nd<Harmonics->Dim; nd++) {
    HARM_NEAR_RW_NORM[nd].open(Combine_String(str, "HarmNearRWnorm%c.dat", Dimension_Label(nd)));
    HARM_NEAR_RW_PHASE[nd].open(Combine_String(str, "HarmNearRWphase%c.dat", Dimension_Label(nd)));
    HARM_NEAR_RTP[nd].open(Combine_String(str, "HarmNearRTp%c.dat", Dimension_Label(nd)));
    if(HarmFreq=='y') {
      HARM_NEAR_RZW_NORM[mpi_rank+nd*mpi_size].open(Combine_String(str, "HarmNearRZWnorm%c_C%d.dat", Dimension_Label(nd), mpi_rank));
      HARM_NEAR_RZW_PHASE[mpi_rank+nd*mpi_size].open(Combine_String(str, "HarmNearRZWphase%c_C%d.dat", Dimension_Label(nd), mpi_rank));
    }
    HARM_FAR_RW_NORM[nd].open(Combine_String(str, "HarmFarRWnorm%c.dat", Dimension_Label(nd)));
    HARM_FAR_RW_PHASE[nd].open(Combine_String(str, "HarmFarRWphase%c.dat", Dimension_Label(nd)));
    HARM_FAR_RTP[nd].open(Combine_String(str, "HarmFarRTp%c.dat", Dimension_Label(nd)));
    for(int nfreq=0; nfreq<Nfreq; nfreq++) {
      HARM_NEAR_RZ_NORM[nfreq+nd*Nfreq].open(Combine_String(str, "HarmNearRZ%2$deVnorm%1$c.dat", Dimension_Label(nd), (int)(FreqCmpts[nfreq]*27.212)));
      HARM_NEAR_RZ_PHASE[nfreq+nd*Nfreq].open(Combine_String(str, "HarmNearRZ%2$deVphase%1$c.dat", Dimension_Label(nd), (int)(FreqCmpts[nfreq]*27.212)));
    }
  }
  ofstream HARM_NEAR_ZW("HarmNearZW.dat");
  ofstream HARM_NEAR_ZWsum("HarmNearZWsum.dat");

  /*
  ofstream ABST("AbsT.dat");
  for(int nt=0; nt<T.Nt; nt++)
    ABST << T.AB[nt] << '\n';
  */

  /***************************************************************************************************************/
  // Initial output of spatial fields:
  for(int nfreq=0; nfreq<Nfreq; nfreq++) {
    if(mpi_rank==0) {
      Harmonics->Output_Radial(Harmonics->Erw, FreqInds[nfreq], 'w', HARM_NEAR_RZ_NORM[nfreq], 'n');
      Harmonics->Output_Radial(Harmonics->Erw, FreqInds[nfreq], 'w', HARM_NEAR_RZ_PHASE[nfreq], 'p');  
    } 
    if(mpi_rank==min(1,mpi_size-1) && Harmonics->Dim==2) {
      Harmonics->Output_Radial(Harmonics->Erw2, FreqInds[nfreq], 'w', HARM_NEAR_RZ_NORM[nfreq+Nfreq], 'n');
      Harmonics->Output_Radial(Harmonics->Erw2, FreqInds[nfreq], 'w', HARM_NEAR_RZ_PHASE[nfreq+Nfreq], 'p');  
    }
  }

  /***************************************************************************************************************/
  // Initial output of temporospatial fields:
  if(mpi_rank==0) {
    Laser->Update_Et_at_RZ(0);
    Laser->Output_Temporal(Laser->Et, PULSE_ZTP[0]);  
    if(Laser->Dim==2)  Laser->Output_Temporal(Laser->Et2, PULSE_ZTP[1]); 
  }
  if(PulseTime=='y') {
    Laser->Output_RadialTemporal(&Laser->Ert[Index2V(mpi_rank*X.Nr_mpi,0,T.Nt)], 't', PULSE_RZTP[mpi_rank], 'r', mpi_active, 1, filter);
    if(Laser->Dim==2)  Laser->Output_RadialTemporal(&Laser->Ert2[Index2V(mpi_rank*X.Nr_mpi,0,T.Nt)], 't', PULSE_RZTP[mpi_sum], 'r', mpi_active, 1, filter);
  }
  if(PulseFreq=='y') { 
    Laser->Output_RadialTemporal(&Laser->Erw[Index2V(mpi_rank*T.Nt_mpi,0,X.Nr)], 'w', PULSE_RZW_NORM[mpi_rank], 'n', mpi_active, filter, 1);
    Laser->Output_RadialTemporal(&Laser->Erw[Index2V(mpi_rank*T.Nt_mpi,0,X.Nr)], 'w', PULSE_RZW_PHASE[mpi_rank], 'p', mpi_active, filter, 1);
    if(Laser->Dim==2) {
      Laser->Output_RadialTemporal(&Laser->Erw2[Index2V(mpi_rank*T.Nt_mpi,0,X.Nr)], 'w', PULSE_RZW_NORM[mpi_sum], 'n', mpi_active, filter, 1);
      Laser->Output_RadialTemporal(&Laser->Erw2[Index2V(mpi_rank*T.Nt_mpi,0,X.Nr)], 'w', PULSE_RZW_PHASE[mpi_sum], 'p', mpi_active, filter, 1);
    }
  }
  if(HarmFreq=='y') { 
    Harmonics->Output_RadialTemporal(&Harmonics->Erw[Index2V(mpi_rank*T.Nt_mpi,0,X.Nr)], 'w', HARM_NEAR_RZW_NORM[mpi_rank], 'n', mpi_active, filter, 1);
    Harmonics->Output_RadialTemporal(&Harmonics->Erw[Index2V(mpi_rank*T.Nt_mpi,0,X.Nr)], 'w', HARM_NEAR_RZW_PHASE[mpi_rank], 'p', mpi_active, filter, 1);
    if(Harmonics->Dim==2) {
      Harmonics->Output_RadialTemporal(&Harmonics->Erw2[Index2V(mpi_rank*T.Nt_mpi,0,X.Nr)], 'w', HARM_NEAR_RZW_NORM[mpi_sum], 'n', mpi_active, filter, 1);
      Harmonics->Output_RadialTemporal(&Harmonics->Erw2[Index2V(mpi_rank*T.Nt_mpi,0,X.Nr)], 'w', HARM_NEAR_RZW_PHASE[mpi_sum], 'p', mpi_active, filter, 1);
    }
  }

//   if(mpi_rank==0)
//     for(int nc=0; nc<colours; nc++)
//       cout << LP[nc]->I0 << '\t' << LP[nc]->lambda << '\t' << LP[nc]->fwhm << '\t' << LP[nc]->t0 << '\t' << LP[nc]->cep << '\t' << LP[nc]->envelope << '\t' << LP[nc]->zf << '\t' << LP[nc]->wz0 << '\t' << LP[nc]->phi << '\n';


// #ifdef useMPI
//   MPI_Barrier(MPI_COMM_WORLD);
// #endif


  if(Propagate=='y') {
    /***************************************************************************************************************/
    // Begin propagation:
    time(&local_start);
    if(mpi_rank==0)  cout << "\nCommencing propagation through target medium:\n";
    for(int nz=0; nz<X.Nz-1; nz++) {   // nz advances E(nz)->E(nz+1) and outputs E(nz+1)
      if(mpi_rank==0 && nz==3) {
	time(&local_finish);  elapsed = difftime(local_finish,local_start);
	cout << "Propagation runtime initially estimated as " << (int)((double)(X.Nz-nz)/nz*elapsed/60) << 'm' << (int)((double)(X.Nz-nz)/nz*elapsed)%60 << "s\n";
      }
      else if(mpi_rank==0 && nz>1 && (nz-1)%((int)(X.Nz/20.0))==0) {
	time(&local_finish);  elapsed = difftime(local_finish,local_start);
	cout << "Propagation " << 5*(int)((nz-1)/(X.Nz/20.0)) << "% complete in " << (int)(elapsed/60) << 'm' << (int)(elapsed)%60 << "s\t...\ttime remaining estimated as " << (int)((double)(X.Nz-nz)/nz*elapsed/60) << 'm' << (int)((double)(X.Nz-nz)/nz*elapsed)%60 << "s\n";
      }
      
      /***************************************************************************************************************/
      // Diffract laser field (frequency domain):
#ifdef useMPI
      MPI_Scatter_Interface(*Laser, X.Nr*T.Nt_mpi, 'w'); 
      Diffract_Field(*Laser, A, B, C, D, X, T, mpi_rank);
      MPI_Allgather_Interface(*Laser, X.Nr*T.Nt_mpi, 'w'); 
#else
      Diffract_Field(*Laser, A, B, C, D, X, T, mpi_rank);
#endif
    
    /***************************************************************************************************************/
    // Transform laser field into time-domain:
      Laser->Transform_Field_IFT();
      
      /***************************************************************************************************************/
    // Plasma defocussing & blue-shifting due to free-electrons (time domain):
#ifdef useMPI
      MPI_Scatter_Interface(*Laser, X.Nr_mpi*T.Nt, 't');  
      Plasma_Defocussing(*Laser, Jet, nz, order, X, T, mpi_rank);
      MPI_Allgather_Interface(*Laser, X.Nr_mpi*T.Nt, 't');  
#else
      Plasma_Defocussing(*Laser, Jet, nz, order, X, T, mpi_rank);
#endif 
    
      /***************************************************************************************************************/
      // Transform laser field back into frequency-domain:
      Laser->Transform_Field_FT();
      
      /***************************************************************************************************************/
      // Laser field outputs:
      if(mpi_rank==0) {
	Laser->Update_Et_at_RZ(0);
	Laser->Output_Temporal(Laser->Et, PULSE_ZTP[0]);  
	if(Laser->Dim==2)  Laser->Output_Temporal(Laser->Et2, PULSE_ZTP[1]); 
      }
      if(PulseTime=='y') {
	Laser->Output_RadialTemporal(&Laser->Ert[Index2V(mpi_rank*X.Nr_mpi,0,T.Nt)], 't', PULSE_RZTP[mpi_rank], 'r', mpi_active, 1, filter);
	if(Laser->Dim==2)  Laser->Output_RadialTemporal(&Laser->Ert2[Index2V(mpi_rank*X.Nr_mpi,0,T.Nt)], 't', PULSE_RZTP[mpi_sum], 'r', mpi_active, 1, filter);
      }
      if(PulseFreq=='y') {
	Laser->Output_RadialTemporal(&Laser->Erw[Index2V(mpi_rank*T.Nt_mpi,0,X.Nr)], 'w', PULSE_RZW_NORM[mpi_rank], 'n', mpi_active, filter, 1);
	Laser->Output_RadialTemporal(&Laser->Erw[Index2V(mpi_rank*T.Nt_mpi,0,X.Nr)], 'w', PULSE_RZW_PHASE[mpi_rank], 'p', mpi_active, filter, 1);
	if(Laser->Dim==2) {
	  Laser->Output_RadialTemporal(&Laser->Erw2[Index2V(mpi_rank*T.Nt_mpi,0,X.Nr)], 'w', PULSE_RZW_NORM[mpi_sum], 'n', mpi_active, filter, 1);
	  Laser->Output_RadialTemporal(&Laser->Erw2[Index2V(mpi_rank*T.Nt_mpi,0,X.Nr)], 'w', PULSE_RZW_PHASE[mpi_sum], 'p', mpi_active, filter, 1);
	}
      }   
      if((nz==0||fabs(X.z[nz-1]-Jet.zj-X.dz/10.0)<=0.5*X.dz) && mpi_rank==0) {  // -X.dz/10 to select the z>zj if two z's equidistant from zj
      	Laser->Update_Et_at_RZ(0);
	Laser->Output_Temporal(Laser->Et, PULSE_TP[0]);  
	if(Laser->Dim==2)  Laser->Output_Temporal(Laser->Et2, PULSE_TP[1]); 
      }  
      if(mpi_rank==0)   Laser->Output_Radial(Laser->Ert, T.Nt/2, 't', PULSE_RZ[0], 'r');
      if(mpi_rank==min(1,mpi_size-1) && Laser->Dim==2) 	Laser->Output_Radial(Laser->Ert2, T.Nt/2, 't', PULSE_RZ[1], 'r');


      if((nz==0||fabs(X.z[nz-1]-X.dz/10.0)<=0.5*X.dz||fabs(X.z[nz-1]-Jet.zj-X.dz/10.0)<=0.5*X.dz||nz==X.Nz-2) && mpi_rank==0) {
	cout << "Pulse energy at z = " << X.z[nz-1]*au_to_m*1e3 << "mm is calculated as E = " << Laser->Pulse_Energy('n')*1.0e3 << " mJ (numerically) & " << Laser->Pulse_Energy('e')*1.0e3 << " mJ (estimate)\n";
	Laser->Analyse_Field(LP_analysed, 0, false, mpi_active);
 	  for(int nc=0; nc<colours; nc++)
	    cout << "I0 = " << pow(LP_analysed[nc]->E0,2)*au_to_Wcm2*1e-14 << "e14Wcm^-2, lambda = " << LP_analysed[nc]->lambda*au_to_m*1e9 << "nm, t0 = " << LP_analysed[nc]->t0*au_to_s*1e15 << "fs, fwhm = " << LP_analysed[nc]->fwhm*au_to_s*1e15 << "fs, cep = " << LP_analysed[nc]->cep/pi << "pi rad\n";
	Laser->Analyse_Field(LP_analysed, 0, true, mpi_active);
 	  for(int nc=0; nc<colours; nc++)
	    cout << "I0 = " << pow(LP_analysed[nc]->E0,2)*au_to_Wcm2*1e-14 << "e14Wcm^-2, lambda = " << LP_analysed[nc]->lambda*au_to_m*1e9 << "nm, t0 = " << LP_analysed[nc]->t0*au_to_s*1e15 << "fs, fwhm = " << LP_analysed[nc]->fwhm*au_to_s*1e15 << "fs, cep = " << LP_analysed[nc]->cep/pi << "pi rad\n";
      }
      
      /***************************************************************************************************************/
      // Reabsorb and diffract propagating harmonic field (frequency domain):
#ifdef useMPI
      MPI_Scatter_Interface(*Harmonics, X.Nr*T.Nt_mpi, 'w'); 
      Absorption_By_Gas(*Harmonics, Jet, nz, X, T, mpi_rank);
      Diffract_Field(*Harmonics, A, B, C, D, X, T, mpi_rank);
      MPI_Allgather_Interface(*Harmonics, X.Nr*T.Nt_mpi, 'w'); 
#else
      Absorption_By_Gas(*Harmonics, Jet, nz, X, T, mpi_rank);
      Diffract_Field(*Harmonics, A, B, C, D, X, T, mpi_rank);
#endif
     
      /***************************************************************************************************************/
      // Generation of new harmonics (use Ert class member here for convenience even though in frequency domain):
#ifdef useMPI // scatter in strips (to enable efficient use of max(Et)>minE to determine whether to calculate harmonics)
      for(int nr=0; nr<X.Nr_mpi; nr++)
	MPI_Scatter_Interface(*Laser, T.Nt, 't', Index2V(nr*mpi_size,0,T.Nt), Index2V(nr,0,T.Nt));
      Generate_Harmonics(*Harmonics, *Laser, minI, Dipole, Jet.nd[nz], X, T);
      for(int nr=0; nr<X.Nr_mpi; nr++)
	MPI_Allgather_Interface(*Harmonics, T.Nt, 't', Index2V(nr,0,T.Nt), Index2V(nr*mpi_size,0,T.Nt)); 
#else
      Generate_Harmonics(*Harmonics, *Laser, minI, Dipole, Jet.nd[nz], X, T);
#endif
      
      /***************************************************************************************************************/
      // Combine newly generated harmonics with propagating harmonic field (in frequency domain)
      Harmonics->Combine_Fields('w');  
      
//       for(int nt=0; nt<T.Nt; nt++) {
// 	HARM_NEAR_ZW << norm(Harmonics->Ert[Index2V(0,nt,T.Nt)]) << '\t';
// 	HARM_NEAR_ZWsum << norm(Harmonics->Erw[Index2V(nt,0,X.Nr)]) << '\t';
//       }
//       HARM_NEAR_ZW << '\n';
//       HARM_NEAR_ZWsum << '\n';
      
      /***************************************************************************************************************/
      // Harmonic field outputs:
      if(HarmFreq=='y') {
	Harmonics->Output_RadialTemporal(&Harmonics->Erw[Index2V(mpi_rank*T.Nt_mpi,0,X.Nr)], 'w', HARM_NEAR_RZW_NORM[mpi_rank], 'n', mpi_active, filter, 1);
	Harmonics->Output_RadialTemporal(&Harmonics->Erw[Index2V(mpi_rank*T.Nt_mpi,0,X.Nr)], 'w', HARM_NEAR_RZW_PHASE[mpi_rank], 'p', mpi_active, filter, 1);
	if(Harmonics->Dim==2) {
 	  Harmonics->Output_RadialTemporal(&Harmonics->Erw2[Index2V(mpi_rank*T.Nt_mpi,0,X.Nr)], 'w', HARM_NEAR_RZW_NORM[mpi_sum], 'n', mpi_active, filter, 1);
 	  Harmonics->Output_RadialTemporal(&Harmonics->Erw2[Index2V(mpi_rank*T.Nt_mpi,0,X.Nr)], 'w', HARM_NEAR_RZW_PHASE[mpi_sum], 'p', mpi_active, filter, 1);
	}
      }

      /***************************************************************************************************************/
      // Spatial (r vs z) outputs at specific frequency components (given by FreqCmpts array):
      for(int nfreq=0; nfreq<Nfreq; nfreq++) {
	if(mpi_rank==0) {
	  Harmonics->Output_Radial(Harmonics->Erw, FreqInds[nfreq], 'w', HARM_NEAR_RZ_NORM[nfreq], 'n');
	  Harmonics->Output_Radial(Harmonics->Erw, FreqInds[nfreq], 'w', HARM_NEAR_RZ_PHASE[nfreq], 'p');  
	} 
	if(mpi_rank==min(1,mpi_size-1) && Harmonics->Dim==2) {
	  Harmonics->Output_Radial(Harmonics->Erw2, FreqInds[nfreq], 'w', HARM_NEAR_RZ_NORM[nfreq+Nfreq], 'n');
	  Harmonics->Output_Radial(Harmonics->Erw2, FreqInds[nfreq], 'w', HARM_NEAR_RZ_PHASE[nfreq+Nfreq], 'p');  
	}
      }
    }
    // end propagation
    if(mpi_rank==0) {
      time(&local_finish);  elapsed = difftime(local_finish,local_start);
      cout << "Propagation completed in " << (int)(elapsed/60) << 'm' << (int)(elapsed)%60 << "s\n";
    }
  }
  
  /***************************************************************************************************************/
  // Transform harmonic field into time-domain and final output of near-field:
  Harmonics->Transform_Field_IFT();

  if(mpi_rank==0)  Harmonics->Output_RadialTemporal(Harmonics->Erw, 'w', HARM_NEAR_RW_NORM[0], 'n');
  if(mpi_rank==1%mpi_size)  Harmonics->Output_RadialTemporal(Harmonics->Erw, 'w', HARM_NEAR_RW_PHASE[0], 'p');
  if(mpi_rank==2%mpi_size)  Harmonics->Output_RadialTemporal(Harmonics->Ert, 't', HARM_NEAR_RTP[0], 'r');
  if(Harmonics->Dim==2) {
    if(mpi_rank==3%mpi_size)  Harmonics->Output_RadialTemporal(Harmonics->Erw2, 'w', HARM_NEAR_RW_NORM[1], 'n');
    if(mpi_rank==4%mpi_size)  Harmonics->Output_RadialTemporal(Harmonics->Erw2, 'w', HARM_NEAR_RW_PHASE[1], 'p');
    if(mpi_rank==5%mpi_size)  Harmonics->Output_RadialTemporal(Harmonics->Ert2, 't', HARM_NEAR_RTP[1], 'r');
  }

  /***************************************************************************************************************/
  // Propagate harmonic near-field through vacuum to far-field:
  time(&local_start);
  if(mpi_rank==0)  cout << "\nCommencing propagation to far-field ... ";
#ifdef useMPI
  MPI_Scatter_Interface(*Harmonics, X.Nr*T.Nt_mpi, 'w'); 
  Hankel_Propagation(*Harmonics, X, T, FPmethod, mpi_rank);
  MPI_Allgather_Interface(*Harmonics, X.Nr*T.Nt_mpi, 'w');
#else
  Hankel_Propagation(*Harmonics, X, T, FPmethod, mpi_rank);
#endif 

  time(&local_finish);  elapsed = difftime(local_finish,local_start);
  if(mpi_rank==0)  cout << "completed in " << (int)(elapsed/60) << 'm' << (int)(elapsed)%60 << "s\n";

  /***************************************************************************************************************/
  // Transform harmonic far-field into time-domain for output:
  Harmonics->Transform_Field_IFT();

  /***************************************************************************************************************/
  // Post-Propagation Bodges:
  Jet.ionised[0] = Jet.ionised[1];

  /***************************************************************************************************************/
  // Final Outputs:
  if(mpi_rank==0)  cout << "\nOutputting remaining data files\n";
  int streams = (Laser->Dim==2) ? 16 : 10; // simply number of output functions below (counting sum of rank0 outputs as 1)
  int *rank_choice = new int[streams];
  Output_Rank_Selector(rank_choice, streams, mpi_size);
  if(mpi_rank==rank_choice[0]) {
    Output_Spec(Laser, LP, X, T, mpi_rank, mpi_size, filter, SPEC);
    Laser->Output_Radial(Laser->Ert, T.Nt/2, 't', PULSE_RZ[0]);
    Laser->Update_Et_at_RZ(0);
    Laser->Output_Temporal(Laser->Et, PULSE_TP[0]);  
    if(Laser->Dim==2)  Laser->Output_Temporal(Laser->Et2, PULSE_TP[1]);
  }
  if(mpi_rank==rank_choice[1]) {
    X.Output(AXISRP, AXISR, AXISZ);  T.Output(AXISTP, AXISW);
    if(Laser->Dim==2)  Laser->Output_Radial(Laser->Ert2, T.Nt/2, 't', PULSE_RZ[1]);
  }
  if(mpi_rank==rank_choice[2]) {
    Jet.Output(Jet.nd, GASJET);
    Jet.Output(Jet.ionised, IONISED);
  }
  if(mpi_rank==rank_choice[3]) {
    Laser->Output_Temporal(Laser->Ew, PULSE_W_NORM[0], 'n');  
    Laser->Output_Temporal(Laser->Ew, PULSE_W_PHASE[0], 'p'); 
    if(Laser->Dim==2) {
      Laser->Output_Temporal(Laser->Ew2, PULSE_W_NORM[1], 'n');  
      Laser->Output_Temporal(Laser->Ew2, PULSE_W_PHASE[1], 'p'); 
    }
  }
  if(mpi_rank==rank_choice[4])  Laser->Output_RadialTemporal(Laser->Erw, 'w', PULSE_RW_NORM[0], 'n');
  if(mpi_rank==rank_choice[5])  Laser->Output_RadialTemporal(Laser->Erw, 'w', PULSE_RW_PHASE[0], 'p');
  if(mpi_rank==rank_choice[6])  Laser->Output_RadialTemporal(Laser->Ert, 't', PULSE_RTP[0], 'r');
  if(mpi_rank==rank_choice[7])  Harmonics->Output_RadialTemporal(Harmonics->Erw, 'w', HARM_FAR_RW_NORM[0], 'n');
  if(mpi_rank==rank_choice[8])  Harmonics->Output_RadialTemporal(Harmonics->Erw, 'w', HARM_FAR_RW_PHASE[0], 'p');
  if(mpi_rank==rank_choice[9])  Harmonics->Output_RadialTemporal(Harmonics->Ert, 't', HARM_FAR_RTP[0], 'r');
  if(Laser->Dim==2) {
    if(mpi_rank==rank_choice[10])  Laser->Output_RadialTemporal(Laser->Erw2, 'w', PULSE_RW_NORM[1], 'n');
    if(mpi_rank==rank_choice[11])  Laser->Output_RadialTemporal(Laser->Erw2, 'w', PULSE_RW_PHASE[1], 'p');
    if(mpi_rank==rank_choice[12])  Laser->Output_RadialTemporal(Laser->Ert2, 't', PULSE_RTP[1], 'r');
  }
  if(Harmonics->Dim==2) {
    if(mpi_rank==rank_choice[13])  Harmonics->Output_RadialTemporal(Harmonics->Erw2, 'w', HARM_FAR_RW_NORM[1], 'n');
    if(mpi_rank==rank_choice[14])  Harmonics->Output_RadialTemporal(Harmonics->Erw2, 'w', HARM_FAR_RW_PHASE[1], 'p');
    if(mpi_rank==rank_choice[15])  Harmonics->Output_RadialTemporal(Harmonics->Ert2, 't', HARM_FAR_RTP[1], 'r');
  }

  /***************************************************************************************************************/
  // Tidy Up & Exit:
  if(mpi_rank==0)  cout << "Freeing up memory & exiting\n";
  CONTROL.close();
  SPEC.close();
  AXISR.close();  AXISZ.close();  AXISTP.close();  AXISW.close();  AXISRP.close();
  GASJET.close();
  IONISED.close();
  for(int nd=0; nd<Laser->Dim; nd++) {
    PULSE_W_NORM[nd].close();  PULSE_W_PHASE[nd].close();
    PULSE_TP[nd].close();  
    PULSE_RZ[nd].close();
    PULSE_RW_NORM[nd].close();  PULSE_RW_PHASE[nd].close();
    PULSE_RTP[nd].close(); 
    if(PulseFreq=='y') {
      PULSE_RZW_NORM[mpi_rank+nd*mpi_size].close();  PULSE_RZW_PHASE[mpi_rank+nd*mpi_size].close();
    }
    if(PulseTime=='y')  PULSE_RZTP[mpi_rank+nd*mpi_size].close();
    PULSE_ZTP[nd].close();
  }
  for(int nd=0; nd<Harmonics->Dim; nd++) {
    HARM_NEAR_RW_NORM[nd].close();  HARM_NEAR_RW_PHASE[nd].close();
    HARM_NEAR_RTP[nd].close();
    if(HarmFreq=='y') { 
      HARM_NEAR_RZW_NORM[mpi_rank+nd*mpi_size].close();  HARM_NEAR_RZW_PHASE[mpi_rank+nd*mpi_size].close();
    }
    HARM_FAR_RW_NORM[nd].close();  HARM_FAR_RW_PHASE[nd].close();
    HARM_FAR_RTP[nd].close();
  }
  HARM_NEAR_ZW.close();
  HARM_NEAR_ZWsum.close();
  for(int nfreq=0; nfreq<Nfreq; nfreq++) {
    HARM_NEAR_RZ_NORM[nfreq].close();  HARM_NEAR_RZ_PHASE[nfreq].close();
  }

  delete[] A, B, C, D;
  delete[] rank_choice;
  for(int nc=0; nc<colours; nc++)  delete[] LP[nc];
  delete[] LP;
  delete[] Laser, Harmonics;
  if(Nfreq!=0)  delete[] FreqCmpts, FreqInds;
  delete[] PULSE_TP, PULSE_W_NORM, PULSE_W_PHASE;
  delete[] PULSE_RTP, PULSE_RW_NORM, PULSE_RW_PHASE;
  delete[] PULSE_RZ, PULSE_RZTP, PULSE_RZW_NORM, PULSE_RZW_PHASE;
  delete[] PULSE_ZTP;
  //delete[] HARM_NEAR_RW_NORM, HARM_NEAR_RW_PHASE, HARM_NEAR_RTP;
  delete[] HARM_NEAR_RZW_NORM, HARM_NEAR_RZW_PHASE;
  delete[] HARM_FAR_RW_NORM, HARM_FAR_RW_PHASE, HARM_FAR_RTP;
  delete[] HARM_NEAR_RZ_NORM, HARM_NEAR_RZ_PHASE;

  if(mpi_rank==0) {
    time(&global_finish);  elapsed = difftime(global_finish,global_start);
    cout << "Total time elapsed = " << (int)(elapsed/60) << 'm' << (int)elapsed%60 << "s\n\n";
  }
#ifdef useMPI
  MPI_Finalize();  // runs before destructors so can't use mpi functions to force sequential
#endif

  return 0;
}

/***************************************************************************************************************/
