/*!
* \file gadget-snapshot-dump-data.cpp
*
*  Created on: Feb 18, 2011
*      Author: blew
*/


#include <stdlib.h>
#include <stdio.h>
#include <cpgplot.h>
#include <math.h>
#include <string.h>
#include <tclap/CmdLine.h>
#include "cpeds-math.h"
#include "Nbody_io.h"
#include "cpeds-msgs.h"
#include "cpeds-rng.h"


#ifndef _NO_NAMESPACE
using namespace std;
//using namespace math;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif

bool _sodTube;
long _Npart;
double _SodNumberDensityRatio, _BoxSizeX, _BoxSizeY, _BoxSizeZ;
string _outputFile;
void parseOptions(int argc, char** argv);
//-----------------------------------

int main(int argc, char **argv) {
	parseOptions(argc,argv);
	matrix<double> ic(_Npart,7);
	Nbody_io snapshot;
	Nbody_io::gadget_IC_blocks block;
	cpedsRNG rnsX;
	cpedsRNG rnsY;
	cpedsRNG rnsZ;

	if (_sodTube) {		
		long NpartLow=_Npart/_SodNumberDensityRatio;
		long NpartHigh=_Npart-NpartLow;
		
		rnsX.seedOffset(0);
		rnsY.seedOffset(1);
		rnsZ.seedOffset(2);
		rnsX.setMinMax(0,_BoxSizeX);
		rnsY.setMinMax(0,_BoxSizeY);
		rnsZ.setMinMax(0,_BoxSizeZ/2);
		for (long i = 0; i < NpartHigh; i++) {
			ic(i,0)=rnsX.getRN();
			ic(i,1)=rnsY.getRN();
			ic(i,2)=rnsZ.getRN();
			ic(i,3)=0;
			ic(i,4)=0;
			ic(i,5)=0;
			ic(i,6)=1;
		}
		
		rnsX.setMinMax(0, _BoxSizeX);
		rnsY.setMinMax(0, _BoxSizeY);
		rnsZ.setMinMax(_BoxSizeZ/2, _BoxSizeZ);
		for (long i = NpartHigh; i < _Npart; i++) {
			ic(i,0)=rnsX.getRN();
			ic(i,1)=rnsY.getRN();
			ic(i,2)=rnsZ.getRN();
			ic(i,3)=0;
			ic(i,4)=0;
			ic(i,5)=0;
			ic(i,6)=1;
		}

		snapshot.writeGadget2data(ic,_outputFile,0,0);
		cpeds_matrix_save(ic,_outputFile+".mat");
	}
	return 0;
}




void parseOptions(int argc, char** argv) {
	long i;
	string::size_type j;
	
	try {
		
		CmdLine cmd("cosmological calculator\n calculates vatious cosmological distances and other useful quantities from well known formulae", ' ', "" );
		
		// 
		// Define arguments
		//
		
		
//		UnlabeledValueArg<string> input_file("input_file","snapshot files prefix","nofile","string");	cmd.add( input_file );
		ValueArg<long> Npart("n", "Npart", "number of particles ( default: 1000 )", false,1000,"long");	cmd.add( Npart);
		ValueArg<double> SodNumberDensityRatio("", "sodNDR", "number densities ratio for the sod's tube test, given by number density in high pressure region to number density in low pressure region( default: 10 )", false,10,"double");	cmd.add(SodNumberDensityRatio );
		ValueArg<double> BoxSizeX("X", "bX", "box size in X direction ( default: 1 )", false,1,"double");	cmd.add(BoxSizeX);
		ValueArg<double> BoxSizeY("Y", "bY", "box size in Y direction ( default: 1 )", false,1,"double");	cmd.add(BoxSizeY);
		ValueArg<double> BoxSizeZ("Z", "bZ", "box size in Z direction  ( default: 100 )", false,100,"double");	cmd.add(BoxSizeZ);
		
		SwitchArg sodTube("","sod", "generate ICs for Sod's tube test", false);	cmd.add( sodTube);
		
//		vector<string> allowedStr;
//		allowedStr.push_back("Coordinates");
//		allowedStr.push_back("Velocities");
//		allowedStr.push_back("IDs");
//		ValueArg<string> block("b", "block", "block to dump ( default: \"Coordinates\" )", false,"Coordinates",allowedStr);	cmd.add( block );
		
		
		ValueArg<string> outputFile("o", "", "name of the output file ( default: ic.out )", false,"ic.out","string");	cmd.add( outputFile );
		
		//
		// Parse the command line.
		//
		cmd.parse(argc,argv);
		
		//
		// Set variables
		//
		
		
//		_inputFile=input_file.getValue();
//		_files=files.getValue();
//		_type=type.getValue();
		_Npart=Npart.getValue();
		_sodTube=sodTube.getValue();
		_SodNumberDensityRatio=SodNumberDensityRatio.getValue();
		_BoxSizeX=BoxSizeX.getValue();
		_BoxSizeY=BoxSizeY.getValue();
		_BoxSizeZ=BoxSizeZ.getValue();
		_outputFile=outputFile.getValue();
		
	} catch ( ArgException& e )
	{ cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
}



