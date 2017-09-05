/*!
 * \file mkFractal.cpp
 *
 *  Project: Mscs
 *  Created on: Oct 18, 2012 10:54:50 PM
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
#include "mandelbrotSet.h"



#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif

string _programVersionString;

void parseOptions(int argc, char** argv);
string getProgramVersionString();

double _xmin,_xmax, _ymin, _ymax;
long _Nx, _Ny, _Nz, _maxIterations;
bool _mk3d;
string _outfile;

int main(int argc, char **argv) {
	
	cpedsMsgs msgs(".mkFractal",false,"",High);
	msgs.setSaveRunWriteMode('a');
	string s;
	//----------------------------------------------------------------------------------------------------
	//PARSE THE COMMAND LINE ARGUMENTS
	//----------------------------------------------------------------------------------------------------
	parseOptions(argc,argv);

	mandelbrotSet ms(_Nx,_Ny,_xmin,_xmax,_ymin,_ymax,_maxIterations);
	ms.mkSet();
	if (_mk3d) {
		long ik=0;
		double maxV=ms.getMaxValue(NULL);
		double minV=ms.getMinValue(NULL);

		mscsFunction3dregc out("fractal3D", CPEDS_defaultVerbosityLevel);
		out.setSizeRange(_Nx,_Ny,_Nz,_xmin,_ymin,0,_xmax,_ymax,1);
		out.allocFunctionSpace();
		out=double(0); // initialize function

		for (long i = 0; i < _Nx; i++) {
			for (long j = 0; j < _Ny; j++) {
				for (long k = 0; k < _Nz; k++) {
					ik=long(double(_Nz-1)/(maxV - minV)*(ms.fRe(i,j,0) - minV));
					if (ik<0) ik=0; 
					if (ik>_Nz-1) ik=_Nz-1;
					out.fRe(i,j,ik)=ms.fRe(i,j,0);
				}
			}
		}
		out.saveDF3(_outfile+".df3",0);
		out.savetxtlin(_outfile+".txt",false);
	}
	else {
		ms.saveSlice(2,0,_outfile,0);
		ms.saveSlice(2,0,_outfile+".im",1);		
	}
	
	return 0;
}


void parseOptions(int argc, char** argv) {
//	long i;
	_programVersionString=getProgramVersionString();
	try {
		
		//     CmdLine cmd("RT4\n\n. ",' ', "0.1" );
		CmdLine cmd("mkFractal - a fractal explorer.\n"
				"\n"
				"example usage: mkFractal"
				"",' ', _programVersionString.c_str() );
		
		// 
		// Define arguments
		//
		
		ValueArg<double> xmin("", "xmin", "xmin (-2)", false,-2,"double");	cmd.add( xmin );
		ValueArg<double> xmax("", "xmax", "xmax (1)", false,1,"double");	cmd.add( xmax );
		ValueArg<double> ymin("", "ymin", "ymin (-1)", false,-1,"double");	cmd.add( ymin );
		ValueArg<double> ymax("", "ymax", "ymax (1)", false,1,"double");	cmd.add( ymax );
//		ValueArg<double> zmin("", "zmin", "specifies Z direction domain range (only mk3d) (0)", false,0,"double");	cmd.add( zmin );
//		ValueArg<double> zmax("", "zmax", "specifies Z direction domain range (only mk3d) (1)", false,1,"double");	cmd.add( zmax );
		
		ValueArg<long> Nx("", "Nx", "Nx (1000)", false,1000,"long");	cmd.add( Nx );
		ValueArg<long> Ny("", "Ny", "Ny (1000)", false,1000,"long");	cmd.add( Ny );
		ValueArg<long> Nz("", "Nz", "Nz (1000)", false,1000,"long");	cmd.add( Nz );
		ValueArg<long> maxIter("i", "iter", "maximal number of iterations (2000)", false,2000,"long");	cmd.add( maxIter );
		SwitchArg mk3d("","mk3d", "saves the output fractal in format readable by povray as a 3d array with the number of cells along"
				"Z direction by Nz. The Z coordinate represents the convergence of sequences and so are the values in the cells (default: false)", false);	cmd.add( mk3d );
		
		ValueArg<string> outfile("o","","output file name",false,"mandelbrot.txt","string"); cmd.add(outfile);
		
//		std::vector<string> allowedStr;
//		allowedStr.push_back("none");
//		allowedStr.push_back("incAIncH");
//		allowedStr.push_back("incADecH");
//		
//		ValuesConstraint<string>* allowedStrNew;
//		allowedStrNew = new ValuesConstraint<string>(allowedStr);
//		ValueArg<string> traj("", "traj", "trajectory type:", false,"none", allowedStrNew);	cmd.add( traj );
//		allowedStr.clear();

		//     UnlabeledValueArg<string> input_file("input_file","file to plot","nofile","string");	cmd.add( input_file );
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
		_xmin = xmin.getValue();
		_xmax = xmax.getValue();
		_ymin = ymin.getValue();
		_ymax = ymax.getValue();
//		_zmin = zmin.getValue();
//		_zmax = zmax.getValue();
		_Nx = Nx.getValue();
		_Ny = Ny.getValue();
		_Nz = Nz.getValue();
		_maxIterations=maxIter.getValue();
		_mk3d=mk3d.getValue();
		
		_outfile = outfile.getValue();
//		_outdir = outdir.getValue(); 	if (_outdir=="") _outdir=".";
		
		
	} catch ( ArgException& e )
	{ cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
	
	
}

string getProgramVersionString() {
	string rev;
#ifdef GENERATE_CVS_VERSION_STRING
	rev="$Revision: 1.1 $";
	QString qrev=rev.c_str();
	QStringList qrevl=qrev.split(' ');
	if (qrevl.size()==3) rev=qrevl[1].toStdString();
#else
	rev="0.1.1";
#endif
	
#ifdef GOMPI
	rev+=" (MPI)";
#endif

	return rev;
}
