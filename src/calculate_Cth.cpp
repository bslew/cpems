#include <stdlib.h>
#include <stdio.h>
#include <cpgplot.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <tclap/CmdLine.h>
#include <omp.h>
#include "cpeds-consts.h"
#include "cpeds-math.h"
#include "Mscs-map.h"
#include "Mscs-global-defs.h"
#include <chrono>


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
double _minang,_maxang,_res, _wisdom_frac;
int _Nthreads,_mapres;
/* long _bin_num; */
bool _masked=false,_plot, _hisamplingathires;
bool _useWisdom,_test_wisdom;
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
	
	if (_test_wisdom) {
		map.set_nside(_mapres);
		map.makekill_space_manager("make","T");
		map.set_map_coord();
		map.make_gaussian_map(0,1);
		
		mscsCorrelationFunction Cth;
		
		typedef std::chrono::high_resolution_clock Clock;
		auto t1 = Clock::now();
		Cth=map.calculate_C_th(_minang,_maxang,_res,1);
		auto t2 = Clock::now();
		double duration_wisdom=double(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count())/1000;
		cout << "Total calculation time (wisdom) [s]: " 
				<< duration_wisdom
				<< std::endl;
		Cth.save("Cth_wisdom.txt");
		
		
		t1=Clock::now();
		double wisdom_frac=_wisdom_frac;
		Cth=map.calculate_C_th(_minang,_maxang,_res,wisdom_frac);
		t2=Clock::now();
		double duration_wisdom_frac=double(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count())/1000;
		cout << "Total calculation time (wisdom fract " << wisdom_frac << ") [s]: " 
				<< duration_wisdom_frac
				<< std::endl;
		Cth.save("Cth_wisdom_frac.txt");
		
		t1=Clock::now();
		Cth=map.calculate_C_th(_minang,_maxang,_res,false);
		t2=Clock::now();
		double duration_normal=double(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count())/1000;
		cout << "Total calculation time (normal) [s]: " 
				<< duration_normal
				<< std::endl;
		//	  cout << "Time increase"
		Cth.save("Cth.txt");
		return 0;
	}
	
	
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
	
	if (_Nthreads>0) omp_set_num_threads(_Nthreads);
	
	if (_CthType=="Cth")  {
		Cth=map.calculate_C_th(_minang,_maxang,_res,_useWisdom);
	}
	if (_CthType=="CthNormVar")  {
		Cth=map.calculate_C_th(_minang,_maxang,_res);
		Cth/=map.get_varianceT();
	}
	if (_CthType=="Sth") Cth=map.calculate_Sth(_minang,_maxang,_res);
	//  outfile = _input_file+_outfile;
	Cth.save(_outfile);
	
	return 0;
}


void parseOptions(int argc, char** argv) {
	try {
		
		CmdLine cmd("calculate_Cth: calculating 2-pt correlation function on a healpix map", ' ', Mscs_version );
		
		UnlabeledValueArg<string> file("files", "input binary file name (prefix) in Mscs format ",true,"", "string");     	cmd.add( file );
		ValueArg<string> CthType("","type","type of correlation function. Possible values are:"
				"Cth - w(th=ang(n1,n2)) = f(n1)*f(n2)/count_in_bin,"
				"CthNormVar - like hist but normalized by variance"
				"DDRR - w(th=ang(n1,n2)) = DD/RR-1 - NOT IMPLEMENTED,"
				"Sth - S(th_ij) = 2*<fi*fj*mi*mj>/(<fi^2*mi*mj> + <fj^2*mi*mj>) "
				"where i,j are map pixels, m is the mask and f are the map values."
				"",false,"out","string"); cmd.add(CthType);
		ValueArg<string> outfile("o","outfile","name of the output file",false,"out","string"); cmd.add(outfile);
		ValueArg<string> mask("m","mask","mask file",false,"","string"); cmd.add(mask);
		ValueArg<int> Nthreads("","Nthreads","Number of threads to use (1). "
				"-1 to use all CPUs.",false,1,"int"); cmd.add(Nthreads);
		ValueArg<double> res("r","res","resolution [deg] (5)",false,5,"double"); cmd.add(res);
		ValueArg<double> minang("f","from","minimal angular scale [deg] (0)",false,0,"double"); cmd.add(minang);
		ValueArg<double> maxang("t","to","maximal angular scale [deg] (180)",false,180,"double"); cmd.add(maxang);
		/* 	ValueArg<long> nsigma("o","output","",false,4,"long"); cmd.add(nsigma); */
		SwitchArg test_wisdom("","test_wisdom","test wisdom (default: false)",false); cmd.add(test_wisdom);
		ValueArg<double> wisdom_frac("","wisdom_frac","use wisdom fraction",false,1,"double"); cmd.add(wisdom_frac);
		ValueArg<int> mapres("","mapres","healpix map resolution (test_wisdom option)",false,8,"int"); cmd.add(mapres);
		SwitchArg record("","record","record and store the map pixel idexes for "
				"faster calculation of multiple maps of the same format (default: false)",false); cmd.add(record);
		
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
		_Nthreads=Nthreads.getValue();
		/* 	_nsigma = nsigma.getValue(); */
		/* 	_show_plots = show_plots.getValue(); */
		_CthType=CthType.getValue();
		_useWisdom=record.getValue();
		_test_wisdom=test_wisdom.getValue();
		_mapres=mapres.getValue();
		_wisdom_frac=wisdom_frac.getValue();
		
	} catch ( ArgException& e )
	{ cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
}



