#include <stdlib.h>
#include <stdio.h>
#include <cpgplot.h>
#include <math.h>
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
filenamestr mask_file;
string _input_file,_mask_file,_outfile, _CthType;
double _minang,_maxang,_res;
/* long _bin_num; */
bool _masked=false,_plot, _hisamplingathires;
//-----------------------------------


int main(int argc, char **argv) {
  long i;
  string outfile;

  //----------------------------------------------------------------------------------------------------
  //PARSE THE COMMAND LINE ARGUMENTS
  //----------------------------------------------------------------------------------------------------
  Mscs_initiate_global_variables();
  mscsMap map("map");
  
  parseOptions(argc,argv);
  //----------------------------------------------------------------------------------------------------
  map.setVerbosityLevel(High);
  map.loadbinT(_input_file);
//  map.clean_mask();
  if (_mask_file!="") { // load mask if required
//    sprintf(mask_file,"%s%s",program_dir,_mask_file.c_str());
//    mask.loadfitsT(mask_file,2); 
//	  mscsMap mask("mask");
//	  mask.loadbinm(mask_file);
//	  map.import_map_data(mask,"m",2);    
//	  map.check_mask();    
//	  map.mask_map_merge();
	  map.loadbinm(_mask_file);
  }
  map.calculate_map_stats(1);

  mscsCorrelationFunction Cth;
  
  if (_CthType=="hist")  Cth=map.calculate_C_th(_minang,_maxang,_res);
  if (_CthType=="histNormVar")  {
	  Cth=map.calculate_C_th(_minang,_maxang,_res);
	  Cth/=map.get_varianceT();
  }
  if (_CthType=="Sstat") Cth=map.calculate_C_th(_minang,_maxang,_res);
//  outfile = _input_file+_outfile;
  Cth.save(_outfile);


}


void parseOptions(int argc, char** argv) {
  try {

	CmdLine cmd("calculate_Cth: calculating 2-pt correlation function on a healpix map", ' ', Mscs_version );

	UnlabeledValueArg<string> file("files", "input binary file name (prefix) in Mscs format ",true,"", "string");     	cmd.add( file );
	ValueArg<string> CthType("","type","type of correlation function. Possible values are:"
			"hist - w(th=ang(n1,n2)) = f(n1)*f(n2)/count_in_bin,"
			"histNormVar - like hist but normalized by variance"
			"DDRR - w(th=ang(n1,n2)) = DD/RR-1 - NOT IMPLEMENTED,"
			"Sstat - S(th_ij) = 2*<fi*fj*mi*mj>/(<fi^2*mi*mj> + <fj^2*mi*mj>) "
			"where i,j are map pixels, m is the mask and f are the map values."
			"",false,"out","string"); cmd.add(CthType);
	ValueArg<string> outfile("o","outfile","name of the output file",false,"out","string"); cmd.add(outfile);
	ValueArg<string> mask("m","mask","mask file",false,"","string"); cmd.add(mask);
	ValueArg<double> res("r","res","resolution [deg] (5)",false,5,"double"); cmd.add(res);
	ValueArg<double> minang("f","from","minimal angular scale [deg] (0)",false,0,"double"); cmd.add(minang);
	ValueArg<double> maxang("t","to","maximal angular scale [deg] (180)",false,180,"double"); cmd.add(maxang);
/* 	ValueArg<long> nsigma("o","output","",false,4,"long"); cmd.add(nsigma); */
//	SwitchArg plot("p","plot","plot the result (default: false)",false); cmd.add(plot);
//	SwitchArg hisamplingathires("","hisamathires","do finer sampling at smaller angular scales (default: false)",false); cmd.add(hisamplingathires);

	//
	// Parse the command line.
	//
	cmd.parse(argc,argv);
	//
	// Set variables
	//
	_input_file = file.getValue(); 
	_outfile = outfile.getValue(); 
	_mask_file = mask.getValue(); 
	_res = res.getValue();
//	_plot = plot.getValue();
//	_hisamplingathires = hisamplingathires.getValue();
	_minang = minang.getValue();
	_maxang = maxang.getValue();
/* 	_nsigma = nsigma.getValue(); */
/* 	_show_plots = show_plots.getValue(); */
	_CthType=CthType.getValue();
 
  } catch ( ArgException& e )
      { cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
}



