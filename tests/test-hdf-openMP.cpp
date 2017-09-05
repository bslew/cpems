/*!
 * \file test-hdf-openMP.cpp
 *
 *  Project: Mscs
 *  Created on: Jan 31, 2016 2:38:34 PM
 *  Author: blew
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <tclap/CmdLine.h>
#include <QtCore/QStringList>
#include <QtCore/QString>
#include "cpeds-msgs.h"
#ifdef ENABLE_OMP
#include <omp.h>
#endif

#include "Mscs-function3dregc.h"


#ifndef _NO_NAMESPACE
using namespace std;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif

string _programVersionString;

void parseOptions(int argc, char** argv);
string getProgramVersionString();

int main(int argc, char **argv) {
	
	cpedsMsgs msgs(".test-hdf-openMP",false,"",High);
	msgs.setSaveRunWriteMode('a');
	string s;
	//----------------------------------------------------------------------------------------------------
	//PARSE THE COMMAND LINE ARGUMENTS
	//----------------------------------------------------------------------------------------------------
	parseOptions(argc,argv);
	long i,k=0;
	mscsFunction3dregc f;
	f.setSize(10,10);
	f.allocFunctionSpace();
	mscsFunction3dregc g=f;
	
	
#ifdef ENABLE_OMP
	omp_lock_t  _lock;
	omp_init_lock(&_lock);
#endif

	double _h0=44;
	
#ifdef ENABLE_OMP
#pragma omp parallel for schedule(guided) default(none) shared(k,_lock) private(i) firstprivate(f,g,_h0,msgs) num_threads(100)
#endif
	for (i = 0; i < 1000; ++i) {
		mscsFunction h;
		// make some time consuming calculation
		for (unsigned long j = 0; j < 1000000; j++) {
#pragma omp atomic
			k++;
		}
		
		// save to the same hdf file from different threads but in a serial manner
#ifdef ENABLE_OMP
		omp_set_lock(&_lock);
#endif
		f.saveHDF5("test.hdf5","f_"+msgs.toStr(i));
		f.setHDF5_scalarStringAttribute("test.hdf5","f_"+msgs.toStr(i),"hubble","uadasdadsnitless Hubble parameter unitless Hubble parameterunitless Hubble parameterunitless Hubble parameterunitless Hubble parameter");

//		g.saveHDF5("test.hdf5","g_"+msgs.toStr(i));
//		g.setHDF5_scalarStringAttribute("test.hdf5","f_"+msgs.toStr(i),"hubble","uadasdadsnitless Hubble parameter unitless Hubble parameterunitless Hubble parameterunitless Hubble parameterunitless Hubble parameter");

//		h.HDF5_scalarStringAttribute("test.hdf5","f_"+msgs.toStr(i),"hubble2","uadasdadsnitless Hubble parameter unitless Hubble parameterunitless Hubble parameterunitless Hubble parameterunitless Hubble parameter");
		
#ifdef ENABLE_OMP
		omp_unset_lock(&_lock);
#endif
	}
#ifdef ENABLE_OMP
	omp_destroy_lock(&_lock);
#endif

	return 0;
}


void parseOptions(int argc, char** argv) {
//	long i;
	_programVersionString=getProgramVersionString();
	try {
		
		//     CmdLine cmd("RT4\n\n. ",' ', "0.1" );
		CmdLine cmd("test-hdf-openMP\nTorun Centre for Astronomy, UMK, POLAND.\n "
				"This program ...\n\n. "
				"\n"
				"example usage: test-hdf-openMP"
				"",' ', _programVersionString.c_str() );
		//     CmdLine cmd("RT4 azimutal scan simulator for OCRA-p\n", 'This program simulates the RT4 scan motion, sky rotation during \n\n. ', 0.1 );
		
		// 
		// Define arguments
		//
		
//		ValueArg<double> Ndec("", "Ndec", "number of periods in Lissajous trajectory [21]", false,21,"double");	cmd.add( Ndec );
//		ValueArg<long> nLissajous("", "nLj", "number of timesteps to take for Lissajous trajectory (5000)", false,5000,"long");	cmd.add( nLissajous );
		
//		SwitchArg horizTrack("","trackHoriz", "follow horizontal sky motion (default: true)", false);	cmd.add( horizTrack );
		
//		ValueArg<string> outfile("o","","output file name",false,"scan-trajectory.txt","string"); cmd.add(outfile);
//		ValueArg<string> outdir("","odir","output directory name",false,"./","string"); cmd.add(outdir);
		
//		std::vector<string> allowedStr;
//		allowedStr.push_back("none");
//		allowedStr.push_back("incAIncH");
//		allowedStr.push_back("incADecH");
		
//		ValuesConstraint<string>* allowedStrNew;
//		allowedStrNew = new ValuesConstraint<string>(allowedStr);
//		ValueArg<string> traj("", "traj", "trajectory type:", false,"none", allowedStrNew);	cmd.add( traj );
//		allowedStr.clear();

		//     UnlabeledValueArg<string> input_file("input_file","file to plot",true,"nofile","string");	cmd.add( input_file );
		//     ValueArg<string> mask("m","mask","wheather to mask the file before doing statistics (prefix)",false,"","string"); cmd.add(mask);
		//     SwitchArg mask_from_here("","MM", "take mask from current directory rather than from the default directory", false);	cmd.add( mask_from_here );
		//     SwitchArg dont_save_mask("","dont_save_mask", "do not print to file masked pixels; useful for -Tn-txt savings", false);	cmd.add( dont_save_mask );
		//     ValueArg<string> proj("p","proj","projection to use [mall - Mollweide/ aitoff etc] (default: mall)",false,"mall","string"); cmd.add(proj);
		//     ValueArg<string> ft("f","ft","input file type [bin, txt] (default: bin)",false,"bin","string"); cmd.add(ft);
		//     ValueArg<string> ord("","ord","ordering of the map [ring, nest] default: nest",false,"nest","string"); cmd.add(ord);
		//     SwitchArg mink("M","Mink", "plot also the minkowski functionals for the map", false);	cmd.add( mink );
		//     ValueArg<long> mlevels("", "mlev", "number of levels to calculate the minkowski functionals (default: 100)", false,100,"long");	cmd.add( mlevels );
		
		cmd.parse(argc,argv);
		
		//
		// Set variables
		//
//		_Ndec = Ndec.getValue();
//		_nLissajous=nLissajous.getValue();
//		_trackHoriz = horizTrack.getValue();
		
//		_outfile = outfile.getValue();
//		_outdir = outdir.getValue(); 	if (_outdir=="") _outdir=".";
		
		
	} catch ( ArgException& e )
	{ cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
	
	
}

string getProgramVersionString() {
	string rev;
#ifdef GENERATE_CVS_VERSION_STRING
	rev="$Revision: 1.1 $";
	QString qrev=rev.c_str();
	QStringList qrevl=qrev.split(" ");
	if (qrevl.size()==3) rev=qrevl[1].toStdString();
	else rev="could not retrieve version number";
#else
	rev="0.1.1";
#endif
	
#ifdef GOMPI
	rev+=" (MPI)";
#endif

	return rev;
}
