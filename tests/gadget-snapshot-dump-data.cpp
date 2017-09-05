/*!
* \file gadget-snapshot-dump-data.cpp
*
*  Created on: Feb 16, 2011
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
#include  "cpeds-msgs.h"


#ifndef _NO_NAMESPACE
using namespace std;
//using namespace math;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif

string _inputFile;
long _files, _type;
string _block;
bool _saveAsType2;

void parseOptions(int argc, char** argv);
//-----------------------------------

int main(int argc, char **argv) {
	parseOptions(argc,argv);
	
	Nbody_io snapshot;
	Nbody_io::gadget_IC_blocks block;

	if (_saveAsType2) {		
		snapshot.convertGadgetToGadget2(_inputFile,_inputFile+".type2");
		
	}
	else {
		if (_block=="Coordinates") { block=Nbody_io::Coordinates; }
		if (_block=="Velocities") { block=Nbody_io::Velocities; }
		if (_block=="IDs") { block=Nbody_io::ParticleIDs; }		
		snapshot.dumpGadgetSnapshotData(_inputFile,_files,_type,block);
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
		
		
		UnlabeledValueArg<string> input_file("input_file","snapshot files prefix","nofile","string");	cmd.add( input_file );
		ValueArg<long> files("n", "files", "number of files in the snapshot( default: 1 )", false,1,"long");	cmd.add( files );
		ValueArg<long> type("t", "type", "gadget snapshot file format ( default: 1 )", false,1,"long");	cmd.add( type );
		
		SwitchArg saveAsType2("","saveAsType2", "saves the input gadget type 1 snapshot file as type gadget type 2 file", false);	cmd.add( saveAsType2);
		
		vector<string> allowedStr;
		allowedStr.push_back("Coordinates");
		allowedStr.push_back("Velocities");
		allowedStr.push_back("IDs");
		ValueArg<string> block("b", "block", "block to dump ( default: \"Coordinates\" )", false,"Coordinates",allowedStr);	cmd.add( block );
		
		
		
		//
		// Parse the command line.
		//
		cmd.parse(argc,argv);
		
		//
		// Set variables
		//
		
		
		_inputFile=input_file.getValue();
		_files=files.getValue();
		_type=type.getValue();
		_block=block.getValue();
		_saveAsType2=saveAsType2.getValue();
		
	} catch ( ArgException& e )
	{ cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
}



