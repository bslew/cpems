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
string _input_file,_mask_file,_outfile,_smB;
long _lmax,_Fmet, _renside,_loaduptol;
double _smG;
filenamestr smooth;
bool masked=false,_pixtf,_save_alms,_almslong,_Folmax,_dont_save_Cl;
//-----------------------------------


int main(int argc, char **argv) {
  filenamestr mask_file;
  double pixtf;
  long l,m;
  //  a_lm a;

  //----------------------------------------------------------------------------------------------------
  //PARSE THE COMMAND LINE ARGUMENTS
  //----------------------------------------------------------------------------------------------------
  Mscs_initiate_global_variables(256,_lmax);
  mscsMap map("map"),mask("mask");
  parseOptions(argc,argv);
  mscsAlms alms("alms");
  //----------------------------------------------------------------------------------------------------

  map.setVerbosityLevel(High);
  map.loadbinT(_input_file,map.nested);
  printf("map loaded\n");
  if (masked) { // load mask if required
    
//    if (_mask_file.find("-Tn-bin",_mask_file.size()-7) != string::npos) { _mask_file.erase(_mask_file.size()-7,7); sprintf(mask_file,"%s%s",program_dir,_mask_file.c_str()); mask.loadbinT(mask_file,4); } else
//      if (_mask_file.find(".fits",_mask_file.size()-5) != string::npos) { _mask_file.erase(_mask_file.size()-5,5);  sprintf(mask_file,"%s%s",program_dir,_mask_file.c_str()); mask.loadfitsT(mask_file,2); } else { sprintf(mask_file,"%s%s",program_dir,_mask_file.c_str()); 
//	mask.loadfits(_mask_file,"T","m"); //}
	  mask.loadbinm(_mask_file);
	  printf("mask loaded\n");
//	  mask.calculate_map_stats(1);

/*     mask.loadfitsT(mask_file,2);  */
    map.import_map_data(mask,"m",2);    mask.killAll();
    map.check_mask();    
    map.mask_map_merge();
//	  exit(0);
//    mask.killAll();
  }

  if (_renside != 0 ) {
	  printf("changing resolution\n");
	  map.change_map_resolution(_renside);
  }


  //smoothing = -1; // no (de)smoothing
  //make_detailed_file_name(&outfile,MAPcalc_dir,256,512,"map",DA,i,DA,i,"nomask","nosmooth",8);
  //alms.lmax(_lmax);
//  Mscs_initiate_fourier_variables(map.nside,2*map.nside);
  if(_pixtf) pixtf = (double)map.nside(); else pixtf = 1;
//  map.SH_analysis(_lmax,)calculate_transformF(GLOBAL_fourier_method_num,pixtf,_smB.c_str(),_smG,0,"",GLOBAL_Plmsfile_forward); 
  alms=map.SH_analysis(_lmax);

/*   // EXPERIMENTAL START */
/*   printf("experimental start\n"); */
/*   for (l=0;l<=map.lmax;l++) { map.set_F(l,0,0,0); } */
/*   printf("experimental end\n"); */
/*   // EXPERIMENTAL END */

  if (_save_alms) { 
    if (_almslong) { alms.savetxtAlms(_outfile,"lmRI"); }
//    if (_Folmax) { alms.savetxtAlms(_outfile,"uptoLmaxOnly",_Folmax);  } 
    if (_Folmax) { alms.savetxtAlms(_outfile,"uptoLmaxOnly");  } 
    else { alms.savebinAlms(_outfile,"RI");   alms.savetxtAlms(_outfile,"RI"); }
  }
  if (!_dont_save_Cl) {
    mscsAngularPowerSpectrum cl=alms.get_Cl(0,_lmax,0);
    cl.save(_outfile+".cl");
    cl.multiply_llpotwoPI();
    cl.save(_outfile+".Cl");
  }
}


void parseOptions(int argc, char** argv) {
  string::size_type j;
  try {
	  
	CmdLine cmd("calculate_alms\n calculates alms from a bin map file", ' ', Mscs_version );

	UnlabeledValueArg<string> infile("infile", "input file name (prefix)",true,"", "string");     	cmd.add( infile );
	ValueArg<string> mask("m","mask","wheather to mask the file before doing statistics (prefix)",false,"","string"); cmd.add(mask);
	ValueArg<string> outfile("o","","name for the output file (prefix)",false,"","string"); cmd.add(outfile);
	ValueArg<long> lmax("l","lmax","maximal multipole number for fourier decomposition",false,512,"long"); cmd.add(lmax);
	ValueArg<long> loaduptol("", "loaduptol", "l up to which load alms from the indicated alms file. NOTE: This works only with the -Folxxx- files now (not with Falxxx files) (default: lmax)", false,-1,"long");	cmd.add( loaduptol );
	ValueArg<long> Fmet("","Fmet","fourier method number: default 1",false,1,"long"); cmd.add(Fmet);
	ValueArg<double> smG("s","smG","fwhm of the Gaussian (de)smoothing beam. -1 for no desmoothing",false,-1,"double"); cmd.add(smG);
	ValueArg<string> smB("","smB","file name with the beam transfer function. default none",false,"","string"); cmd.add(smB);
	ValueArg<long> renside("","renside","change the resolution of the map before doing SHT",false,0,"long"); cmd.add(renside);
	SwitchArg pixtf("","pixtf", "use pix.tf. for given resolution", false);	cmd.add(pixtf);
	SwitchArg save_alms("","save_alms", "saves also the alms (default: true)", false);	cmd.add(save_alms);
	SwitchArg almslong("","almslong", "saves the alms in long format: with the l,m information (default: false)", false);	cmd.add(almslong);
	SwitchArg dont_save_Cl("","dont_save_Cl", "do not save it :) (default: false)", false);	cmd.add(dont_save_Cl);

	SwitchArg Folmax("","Folmax", "indicates the preferred file format to be saved. (default false) The default format is Falmax.", false);	cmd.add(Folmax);

	//
	// Parse the command line.
	//
	cmd.parse(argc,argv);
	//
	// Set variables
	//
	_input_file = infile.getValue(); //j=_input_file.find("-Tn-bin",0); if (j != string::npos) { _input_file.erase(j,(string::size_type)7); }
	_mask_file = mask.getValue(); if (_mask_file == "") { masked = false; } else {masked = true; }
	_lmax = lmax.getValue(); 
	_loaduptol = loaduptol.getValue(); if (_lmax != _loaduptol) {  _loaduptol = _lmax; }
	_smG = smG.getValue(); 
	_smB = smB.getValue(); 
	_outfile = outfile.getValue(); if (_outfile == "") { _outfile = _input_file; }
	_Fmet = Fmet.getValue(); GLOBAL_fourier_method_num = _Fmet;
	_renside = renside.getValue(); 
	_pixtf = pixtf.getValue();
	_almslong = almslong.getValue();
	_save_alms = save_alms.getValue();
	_dont_save_Cl = dont_save_Cl.getValue();

	_Folmax = Folmax.getValue();
 
  } catch ( ArgException& e )
      { cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
}



