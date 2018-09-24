/***************************************************************************************************************/
// Header file to contain universal input functions for Propagator.cpp 

/***************************************************************************************************************/
// Notes:

/***************************************************************************************************************/
#ifndef INPUT_H
#define INPUT_H

#include<fstream>
#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include <ctype.h>
#include<unistd.h>
#include<Misc_Functions.h>

/***************************************************************************************************************/
// Function Prototypes:
void Read_Command_Options(char *&, char *, int, char **);
template<class T> T Read_Variable(std::ifstream &, const char *, const int=1);

/***************************************************************************************************************/
// Functions:

// Function to handle Getopt arguments
void Read_Command_Options(char *&ConfigPath, char *ConfigDefault, int argc, char **argv) {
  ConfigPath = ConfigDefault;
  int opt;
  opterr = 0;
  while((opt=getopt(argc,argv,"C:"))!=-1)
    switch(opt) {
    case 'C':
      ConfigPath = optarg;
      break;
    case '?':
      if(optopt=='C')
	fprintf(stderr,"Option -%c requires an argument.\n",optopt);
      else if(isprint(optopt))  // ctype.h only used for isprint
	fprintf(stderr,"Unknown option `-%c'.\n",optopt);
      else
	fprintf(stderr,"Unknown option character `\\x%x'.\n",optopt);
      abort();
    default:
      abort();
    }
  for(int index=optind; index<argc; index++)
    printf("Non-option argument %s\n",argv[index]);

}

// Function to read in a specified variable from a filestream 
template<class T> T Read_Variable(std::ifstream &INSTREAM, const char *varname, const int Ninst) {
  T var;
  INSTREAM.clear();  
  INSTREAM.seekg(std::ios::beg);
  if(INSTREAM.good()) {
    const int Nch = strlen(varname);
    char teststr[Nch];
    char testch;
    bool match(false);
    int ninst = 0;
    while(!INSTREAM.eof() && !match) {
      for(int nch=0; nch<Nch; nch++) {
	INSTREAM.get(testch);
	if(INSTREAM.eof())  break;
	else if(testch==('%'||'='||';')) {
	  nch = -1;
	  INSTREAM.ignore(1000,'\n');
	}
	else if(isspace(testch) || (ispunct(testch)&&testch!='_'))  nch = -1;
	else  teststr[nch] = testch;
      }
      if(!INSTREAM.eof() && strncmp(varname, teststr, Nch)==0) {
	INSTREAM >> testch;
	if(isspace(testch) || testch=='=') {  // checking this variable name not contained within another variable name
	  ninst++;
	  if(ninst==Ninst) {
	    INSTREAM >> testch;
	    while(isspace(testch) || testch=='=')  INSTREAM >> testch;
	    INSTREAM.putback(testch);
	    INSTREAM >> var;
	    match = true;
	    // break; ---> there is an issue with variable names repeated in the comments: maybe this would fix it?
	  }
	}
	//else {}
      }
    }
    if(!match) {
      char str[50];
      std::cerr << Combine_String(str, "Cannot find variable %s in file!\n", varname);
#ifdef useMPI
      MPI_Abort(MPI_COMM_WORLD, 1);
      MPI_Finalize();
#endif
      exit(1);
    }
  }
  else {
    std::cerr << "Cannot find source file!\n";
#ifdef useMPI
    MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_Finalize();
#endif
    exit(1);
  }
  return var;
}

/***************************************************************************************************************/
#endif
