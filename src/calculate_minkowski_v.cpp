#include <stdlib.h>
#include <stdio.h>
#include <cpgplot.h>
#include <math.h>
#include <fitsio.h>
#include <string.h>
#include <tclap/CmdLine.h>
#include "cpeds-consts.h"
#include "cpeds-math.h"
//#include "chealpix.h"
#include "Mscs-map.h"
#include "Mscs-global-defs.h"
#include "Mscs-minkowski-f.h"


#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
//using namespace MATPACK; // use all classes and functions in the MATPACK namespace
using namespace TCLAP;
#define STD std
#else
#define STD
#endif


//declaration of the global variables
void parseOptions(int argc, char** argv);
filenamestr _input_file_my;
string _input_file,_output_file,_mask_file,_ft,_ord,_functype;
double _nsigmaMin,_nsigmaMax,_maxFskyMask;
long _mlevels;
bool _col_ord, _no_normdev;
bool masked=false,_mink=false;
//-----------------------------------


int main(int argc, char **argv) {
  filenamestr file2;
  mscsMap::mapOrderings mapordering;
  string fname;
  minkowski_fs *mkfs;

  //----------------------------------------------------------------------------------------------------
  //PARSE THE COMMAND LINE ARGUMENTS
  //----------------------------------------------------------------------------------------------------
  Mscs_initiate_global_variables();
  mscsMap map("map"),map2("mask");
  parseOptions(argc,argv);
  //----------------------------------------------------------------------------------------------------

/*   strcpy(mappref, ARGV[1]); */
/*   mapordering = (int)(strtol(ARGV[2],NULL,10)); */
/*   projection = (int)(strtol(ARGV[3],NULL,10)); */

  //if ((ARGC != 10) && (ARGC != 5)) { print_program_usage(); exit(1);}
/*  if (ARGC == 5) { */
/*     strcpy(txtbin,ARGV[4]); */
/*     //if (strcmp(txtbin,"") == 0) { strcpy(txtbin,"bin"); } */
/*   } */

/*   if (ARGC == 11) { */
/*     start = (int)(strtol(ARGV[4],NULL,10)); */
/*     stop = (int)(strtol(ARGV[5],NULL,10)); */
/*     step = (int)(strtol(ARGV[6],NULL,10)); */
/*     minka = (int)(strtol(ARGV[7],NULL,10)); */
/*     minkc = (int)(strtol(ARGV[8],NULL,10)); */
/*     minks = (int)(strtol(ARGV[9],NULL,10)); */
/*     strcpy(txtbin,ARGV[10]); */
/*   } */
/*   else { start = 1; stop = 1; step = 1; minka = minkc = minks = 1; } */

/*   sim1.minkowski_level_num_circ = 100; */
/*   for (i=start;i<=stop;i=i+step) { */
  //sprintf(out,"simulated_map_sum_lmax%i",i);
/*     if (ARGC == 11) { sprintf(out,"%s%i",mappref,i); } else  { sprintf(out,"%s",mappref); } */
/*     if (strcmp(txtbin,"bin") == 0) { sim1.loadbinT(out,mapordering); } */
/*     if (strcmp(txtbin,"txt") == 0) { sim1.loadtxtT(out,mapordering); } */
  if (_ord == "nest") { mapordering = mscsMap::nested; } 
  if (_ord == "ring") { mapordering = mscsMap::ring; } 

  if (_ft == "bin") { map.loadbinT(_input_file,mapordering); }
  if (_ft == "txt") { map.loadtxtT(_input_file,mapordering); }

  if (masked) {   // load mask if required
//    fname=program_dir+_mask_file;
    fname=_mask_file;
    if (fname.find("-bin",fname.size()-4) != string::npos) { 
    	map2.loadbinm(fname); 
    	map2.change_map_resolution(map.nside());
    } 
/*
    else  //map2.change_map_resolution(map1.nside); 
      //    sprintf(file2,"%s%s",program_dir,_mask_file.c_str());
      map2.loadfitsT(fname.c_str(),2); //map2.change_map_resolution(map.nside); 
*/
    map.import_map_data(map2,"m",2);    
//    map.m()=map2.m();
    map2.killAll();
    map.check_mask();
//    map.multimask_map_merge();
  }

  map.calculate_map_stats(1);

  if ((_functype == "v0") ) {
    map.calculate_minkowski_v0(_mlevels,_nsigmaMin,_nsigmaMax).save(_output_file);
//    map.savetxtM(_output_file.c_str(),10);
  }

  if ((_functype == "v1") ) {
    map.calculate_minkowski_v1(_mlevels,_nsigmaMin,_nsigmaMax).save(_output_file);
//    if (_plotMF) map.plot_minkowski_v(1,1);
/*     fname = _input_file; */
//    map.savetxtM(_output_file.c_str(),11);
  }
  if ((_functype == "v2") ) {
    map.calculate_minkowski_v2(_mlevels,_nsigmaMin,_nsigmaMax).save(_output_file);
//    if (_plotMF) map.plot_minkowski_v(2,1);
/*     fname = _input_file; */
//    map.savetxtM(_output_file.c_str(),12);
  }

  if (_functype == "all") {
//	    mscsFunction v0=map.calculate_minkowski_v0(_mlevels,_nsigmaMin,_nsigmaMax);
//	    mscsFunction v1=map.calculate_minkowski_v1(_mlevels,_nsigmaMin,_nsigmaMax);
//	    mscsFunction v2=map.calculate_minkowski_v2(_mlevels,_nsigmaMin,_nsigmaMax);
	    mscsFunction3dregc mink=map.calculate_minkowski_v0v1v2(_mlevels,_nsigmaMin,_nsigmaMax,!_no_normdev);
	    mink.saveSlice(2,0,_output_file);
//    mkfs=map.calculate_minkowski_v0v1v2(_mlevels,-_nsigma,_nsigma,true,true,NULL,NULL);
/*     mkfs=map.calculate_minkowski_v0v1v2(_mlevels,map.minT,map.maxT,true,true,NULL,NULL); */

/*     map.mink_v0 = new minkowski_f("v0",_mlevels); */
/*     map.mink_v1 = new minkowski_f("v1",_mlevels); */
/*     map.mink_v2 = new minkowski_f("v2",_mlevels); */

/*     map.plot_minkowski_v(2,1); */
/*     map.savetxtM(fname.c_str(),10); */
/*     map.savetxtM(fname.c_str(),11); */
/*     map.savetxtM(fname.c_str(),12); */

/*     mkfs[0].mfs->v0->save(_output_file+"v0"); */
/*     mkfs[0].mfs->v1->save(_output_file+"v1"); */
/*     mkfs[0].mfs->v2->save(_output_file+"v2"); */

//    mkfs->save(_output_file,1);
  }

//  map.kill_map();
}


void parseOptions(int argc, char** argv) {
  //  long i;
	try {

	CmdLine cmd("Calculate the minkowski functionals on map", ' ', Mscs_version );

	// 
	// Define arguments
	//

/* 	SwitchArg btest("B","existTestB", "tests for the existence of B", false);	cmd.add( btest ); */
/* 	SwitchArg ctest("C","existTestC", "tests for the existence of C", false);	cmd.add( ctest ); */
/* 	SwitchArg atest("A","existTestA", "tests for the existence of A", false);	cmd.add( atest ); */
/* 	ValueArg<string> stest("s","stringTest","string test",true,"homer","string");	cmd.add( stest ); */
/* 	ValueArg<int> itest("i", "intTest", "integer test", true, 5, "int");	        cmd.add( itest );  */
/* 	ValueArg<double> ftest("lf", "LookFrom", "position of the camera", false, 3.7, "float");	cmd.add( ftest ); */
/* 	UnlabeledValueArg<string> utest("unTest","unlabeld test","default","string");	cmd.add( utest ); */
/* 	UnlabeledMultiArg<string> mtest("fileName", "file names", "string");     	cmd.add( mtest ); */
/* 	MultiArg<float> ftest("f", "floatTest", "multi float test", false,"float" );	cmd.add( ftest ); */
	//MultiArg<double> camera("l", "LookFrom", "position of the camera", false,"double");	cmd.add( camera );
	//SwitchArg col_ord("C","ColOrdering", "Column major ordering", false);	cmd.add( col_ord );
	//ValueArg<double> size("s","size","Y size of the plots",false,600,"double"); cmd.add(size);

	UnlabeledValueArg<string> input_file("input_file","map file",true, "","string");	cmd.add( input_file );
	ValueArg<string> mask("m","mask","wheather to mask the file before doing statistics (prefix)",false,"","string"); cmd.add(mask);
	ValueArg<string> ft("f","ft","input file type [bin, txt] (default: bin)",false,"bin","string"); cmd.add(ft);
	ValueArg<string> ord("","ord","ordering of the map [ring, nest] default: nest",false,"nest","string"); cmd.add(ord);
 	SwitchArg no_normdev("","nonormdev", "normalize the map by standard deviation before calculating the functionals (false)", false);	cmd.add( no_normdev); 
	ValueArg<long> mlevels("l", "mlev", "number of levels to calculate the minkowski functionals (default: 100)", false,100,"long");	cmd.add( mlevels );
	ValueArg<double> nsigmaMin("", "nsigMin", "number of std.dev. to calculate the minkowski functionals from (default: -4)", false,-4.0,"double");	cmd.add( nsigmaMin );
	ValueArg<double> nsigmaMax("", "nsigMax", "number of std.dev. to calculate the minkowski functionals to (default: 4)", false,4.0,"double");	cmd.add( nsigmaMax );
	ValueArg<string> functype("t","type","which functional to calculate [v0,v1,v2,all]",false,"all","string"); cmd.add(functype);
	ValueArg<string> output("o","out","outfile name (prefix), by default input file name",false,"","string"); cmd.add(output);
//	SwitchArg plotMF("p","plot", "plot the functional", false);	cmd.add( plotMF );
	ValueArg<double> maxFskyMask("", "maxFskyMask", "Maximal allowed masked sky fraction."
			"If the provided sky mask contains fraction of pixels masked equal or above this threshold"
			"then the program will return zeros. By default (1) the program will return"
			"zeros only when whole sky is masked.", false,1.0,"double");	cmd.add( maxFskyMask );

	//
	// Parse the command line.
	//
	cmd.parse(argc,argv);

	//
	// Set variables
	//
	_input_file = input_file.getValue();
	_mask_file = mask.getValue(); if (_mask_file.length() == 0) { masked = false; } else {masked = true; }
	_output_file = output.getValue(); if (_output_file=="") { _output_file=_input_file+"mink"; }
	_ft = ft.getValue();
	_ord = ord.getValue();
 	_no_normdev = no_normdev.getValue(); 
	_mlevels = mlevels.getValue();
	_nsigmaMin = nsigmaMin.getValue();
	_nsigmaMax = nsigmaMax.getValue();
	_functype = functype.getValue();
	_maxFskyMask = maxFskyMask.getValue();
//	_plotMF = plotMF.getValue();


	} catch ( ArgException& e )
	{ cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
}



