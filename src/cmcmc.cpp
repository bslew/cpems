/*!
 * \file cmcmc.cpp
 *
 *  Created on: Nov 19, 2010, 6:59:22 PM
 *      Author: blew
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <tclap/CmdLine.h>
#include "cpedsMCMC.h"

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
string _input_file,_outfile;
cpedsList<long> _bintab;
cpedsList<double> _w;
/* long _bin_num; */
bool _plot, _eqw,_geo;
long _fsize, _from;
double _dl,_gm;
//-----------------------------------


int main(int argc, char **argv) {
  long i;
  string outfile;
/*   long bintabs; */

  //----------------------------------------------------------------------------------------------------
  //PARSE THE COMMAND LINE ARGUMENTS
  //----------------------------------------------------------------------------------------------------
  mscsFunction f;

  parseOptions(argc,argv);
  f.load(_input_file);
//  f.mkSin(0,1,0.001);
//  f.save("sin");
  _fsize=f.pointsCount();


  if (_bintab.size()==0){
	  if (_geo) {
		  double b=_dl;
		  i=_from;
		  while (i<_fsize) {
			  i+=long(round(b));
			  if (i<=_fsize) { _bintab.append(long(round(b))); } else { _bintab.append(_fsize-i+long(round(b))); }
			  b*=_gm;
		  }
	  }
	  else {
		  i=_from;
		  while (i<_fsize) {
			  i+=_dl;
			  if (i<=_fsize) _bintab.append(_dl); 
		  }		  
	  }
  }

  
  f=f.binFunctionLin2(_from,_bintab,_w);

  f.save(_outfile);

}


void parseOptions(int argc, char** argv) {
  long i;
  try {

	CmdLine cmd("bin_function: program for binning the input function read from the 2-column text file", ' ', "");

	UnlabeledValueArg<string> file("files", "input binary file name (prefix) in Mscs format ","", "string");     	cmd.add( file );
	ValueArg<string> outfile("o","outfile","name of the output file (postfix)",false,"","string"); cmd.add(outfile);
	ValueArg<long> from("f","from","bin from l (default: 10)",false,0,"long"); cmd.add(from);
//	ValueArg<long> to("t","to","bin to l (default: 1024)",false,1024,"long"); cmd.add(to);
	ValueArg<long> dl("d","dl","bin size in l (default: 10)",false,10,"long"); cmd.add(dl);
	MultiArg<long> b("b","bin","list of arbitrary bins to use for binning (if defined, the f,t,d args will be ignored)",false,"long"); cmd.add(b);
/* 	ValueArg<long> nsigma("o","output","",false,4,"long"); cmd.add(nsigma); */
	SwitchArg eqw("","eqw","equal weights for binning (default: false). NOT IMPLEMENTED YET",false); cmd.add(eqw);
	SwitchArg geo("","geo","make bins in geometric sequence using multiplier gm and dl (d) parameters (default: false).",false); cmd.add(geo);
	ValueArg<double> gm("","gm","geometric sequence multiplied (default: 2)",false,2,"double"); cmd.add(gm);
//	SwitchArg interpolate("","interpolate","after the binning is done do the cspline interpolation in the range of initial multipoles  (default: false)",false); cmd.add(interpolate);

/* 	SwitchArg hisamplingathires("","hisamathires","do finer sampling at smaller angular scales (default: false)",false); cmd.add(hisamplingathires); */

	//
	// Parse the command line.
	//
	cmd.parse(argc,argv);
	//
	// Set variables
	//
	_input_file = file.getValue(); 
	_outfile = outfile.getValue(); 

	_dl = dl.getValue();
	_from = from.getValue();
//	_to = to.getValue();
	_eqw = eqw.getValue();  
	_eqw=true;

	_geo = geo.getValue();  
	_gm=gm.getValue();
//	_w.append(0);

//	_interpolate = interpolate.getValue();


	vector<long> bins = b.getValue();
	if (bins.size()>0) {
		for (i=0;i<bins.size();i++) { _bintab.append(bins[i]); }
	}

 
  } catch ( ArgException& e )
      { cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
}



