#include <stdlib.h>
#include <stdio.h>
#include <cpgplot.h>
#include <math.h>
#include <fitsio.h>
#include <string.h>
#include <time.h>
#include <tclap/CmdLine.h>
#include "cpeds-consts.h"
#include "cpeds-math.h"
#include "Mscs-map.h"
#include "Mscs-global-defs.h"


#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif


//declaration of the global variables
void parseOptions(int argc, char** argv);
string _input_file,_mask_file;
long _lmax;
bool masked=false;
//-----------------------------------


int main(int argc, char **argv) {

  //----------------------------------------------------------------------------------------------------
  //PARSE THE COMMAND LINE ARGUMENTS
  //----------------------------------------------------------------------------------------------------
  Mscs_initiate_global_variables();
  mscsMap map("map"),mask("mask");
  parseOptions(argc,argv);
  mscsAlms alms("alms");
  //----------------------------------------------------------------------------------------------------

  alms.lmax(_lmax);
  map.loadtxtF(_input_file.c_str(),1);
  map.calculate_C_l(0,map.lmax,1);
  map.savetxtC_l(_input_file.c_str(),1);
  map.killalmsF();

}


void parseOptions(int argc, char** argv) {
  try {
	  
	CmdLine cmd("calculate_alms\n calculates alms from a bin map file", ' ', Mscs_version );

	UnlabeledValueArg<string> infile("infile", "input file name (prefix)","", "string");     	cmd.add( infile );
/* 	ValueArg<string> mask("m","mask","wheather to mask the file before doing statistics (prefix)",false,"","string"); cmd.add(mask); */
	ValueArg<long> lmax("l","lmax","maximal multipole number for fourier decomposition",false,512,"long"); cmd.add(lmax);
	//
	// Parse the command line.
	//
	cmd.parse(argc,argv);
	//
	// Set variables
	//
	_input_file = infile.getValue();
/* 	_mask_file = mask.getValue(); if (_mask_file == "") { masked = false; } else {masked = true; } */
	_lmax = lmax.getValue();
 
  } catch ( ArgException& e )
      { cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
}



