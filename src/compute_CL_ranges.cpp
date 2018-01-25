#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <tclap/CmdLine.h>
#include "cpeds-math.h"
#include "cpeds-pdf.h"
#include "Mscs-function.h"

#ifndef _NO_NAMESPACE
using namespace std;
//using namespace math;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif


void parseOptions(int argc, char** argv);
cpeds_queue<double> _CL;
double _CLr,_CLacc;
bool _testGauss;
bool _doCLr,_nogaps,_save_CLr,_show_info,_save_PDF_info;
string _input_file,_CLr_file,_save_PDF_info_file;
long _modesnum;

void makeGaussTest(cpeds_queue<double>& X, cpeds_queue<double>& Y);

int main(int argc, char** argv) {
  double x,y;
  cpeds_queue<double> X,Y;
  FILE* f;
  cpeds_PDF *pdf=NULL;
  long i;
  bool gaps;
  cpeds_queue_double * modes;
  cpeds_queue_long * modesi;


  parseOptions(argc,argv);
  _CL.sort(12);
  
  if (_testGauss) {
	  makeGaussTest(X,Y);
  }
  else {
	  printf("* reading the PDF function\n");
	  f=fopen(_input_file.c_str(),"r");
	  if (f==NULL) { printf("ERROR: the input file not found. exiting... \n"); exit(0); }
	  while (fscanf(f,"%lE %lE",&x,&y)!=EOF) { X.addq(x); Y.addq(y); }
	  fclose(f);
  }

  //
  // initiate the PDF object
  //
  if (_CL.get_size() > 0) {
    pdf=new cpeds_PDF(X.export_array(),Y.export_array(),X.get_size(),&_CL);  
  } 
  else {
    pdf=new cpeds_PDF(X.export_array(),Y.export_array(),X.get_size()); 
  }

  //
  // derive the CL ranges
  //
  printf("* analysing the PDF\n");
  if (_CL.get_size() > 0) {
    pdf->get_info();


    //
    // print out the information
    //
    double pdfMode=pdf->get_MODE();
    double pdfEX=pdf->get_EX();
    printf("MOD: %lE\n",pdfMode);
    printf("EX: %lE\n",pdfEX);

    for (i=0;i<_CL.get_size();i++) {
      printf("CLr( %lE ): %lE %lE ",_CL(i),pdf->get_CL_range(true,i,&gaps),pdf->get_CL_range(false,i,NULL));
      if (gaps) printf(" possible gaps\n"); else printf("\n");
    }

    for (i=0;i<_CL.get_size();i++) {
      printf("errors w.r.t MOD( %lE ): %lE %lE ",_CL(i),pdfMode-pdf->get_CL_range(true,i,&gaps),pdf->get_CL_range(false,i,NULL)-pdfMode);
      if (gaps) printf(" possible gaps\n"); else printf("\n");
    }

    if (_save_PDF_info) {
      f=fopen(_save_PDF_info_file.c_str(),"a"); 
      fprintf(f,"%lE %lE ",pdf->get_EX(),pdf->get_MODE()); 
      for (i=0;i<_CL.get_size();i++) { fprintf(f,"%lE %lE %lE ",pdf->get_CL(i),pdf->get_CL_range(true,i,NULL),pdf->get_CL_range(false,i,NULL)); } 
      fprintf(f,"\n"); 
      fclose(f);
    }
  }


  if (_show_info) {
    pdf->get_info();
    printf("|compute_CL_ranges> * MOD: %lE\n",pdf->get_MODE());
    printf("|compute_CL_ranges> * EX: %lE\n",pdf->get_EX());

    pdf->get_PDF_modes(&modes,&modesi);
    printf("|compute_CL_ranges> * the mode values are at:\n");
    modes->printq_double();
    printf("|compute_CL_ranges> * the mode values indexes are:\n");
    modesi->printq_long();
    printf("|compute_CL_ranges> * number of modes: %li\n",modes->get_size());
  }

  //
  // derive the CL corresponding to the requested lower or upper limit CLr (the other limit will be returned and the corresponding CL)
  //

  if (_doCLr) {
/*     _CL.delq(); */
/*     _CL.mk_list(_CLacc,1-_CLacc,_CLacc); */
/*     _CL.printq_double(); */

/*     if (pdf!=NULL) { delete pdf; pdf=NULL; } */
/*     pdf=new cpeds_PDF(X.export_array(),Y.export_array(),X.get_size()); */
/*     pdf->get_info(); */

    pdf->get_PDF_modes(&modes,&modesi);
    pdf->integrate_to_CLr(_CLr,!(_nogaps),_modesnum);
    printf("CLr( %lE ): %lE %lE ",pdf->get_CL(0),pdf->get_CL_range(true,0,&gaps),pdf->get_CL_range(false,0,NULL));
    if (gaps) printf(" possible gaps\n"); else printf("\n");
    if (_save_CLr) { pdf->save_CLr(_CLr_file,0); }

    if (_save_PDF_info) {
      f=fopen(_save_PDF_info_file.c_str(),"a"); 
      fprintf(f,"%lE %lE ",pdf->get_EX(),pdf->get_MODE()); 
      fprintf(f,"%lE %lE %lE ",pdf->get_CL(0),pdf->get_CL_range(true,0,NULL),pdf->get_CL_range(false,0,NULL)); 
      fprintf(f,"\n"); 
      fclose(f);
    }

  }
  

  
  if (pdf!=NULL) { delete pdf; pdf=NULL; }

}





void parseOptions(int argc, char** argv) {
  long i;
  string::size_type j;


	try {

	  CmdLine cmd("Given a PDF function in the input file it computes the expectation value, mod value and CL ranges corresponding to given CL", ' ', "" );

	// 
	// Define arguments
	//

//	  UnlabeledValueArg<string> input_file("input_file","file with input PDF","","string");	cmd.add( input_file );
	  UnlabeledValueArg<string> input_file("input_file","file with input PDF",false,"","string");	cmd.add( input_file );
	  
	  MultiArg<double>CL("","CL","CL up to which integrate the PDF sideways from the mode value (must be beween 0 and 1)",false,"double"); cmd.add(CL);
	  SwitchArg doCLr("","doCLr", "whether or not to proceed with integrations as defined in options CLr", false);	cmd.add( doCLr );
	  ValueArg<double> CLr("","CLr","value of the lower CL range down to which to the sideways integration; will return the corresponding upper CL range and CL; if the CLr option is not enabled then this option will not be used",false,0,"double"); cmd.add(CLr);
	  ValueArg<long> modesnum("m","modes","numer of modes to enforce to run over in the integration (default: 0 - don't enforce anything)",false,0,"long"); cmd.add(modesnum);

	  ValueArg<double> CLacc("","CLacc","accuracy in integrating the CL in case of CLr stuff (default: 0.0001)",false,0.0001,"double"); cmd.add(CLacc);

	  SwitchArg nogaps("","nogaps", "whether or not to care about possible gaps in the derived CLRs (for eg. multi-modal PDFs) when eg. integrating the PDF according to some given values of CLrlow/high", false);	cmd.add( nogaps );
	  SwitchArg save_CLr("","save_CLr", "whether or not to save a file the points belonging to the current confidence range (the arguments of the PDF used for the integration)", false);	cmd.add( save_CLr );
	  SwitchArg show_info("","show_info", "whether or not to print the info on the modes of the PDF", false);	cmd.add( show_info);
	  SwitchArg save_PDF_info("","save_PDF_info", "whether or not to save the PDF info to file: EX MODE 68 95 99 percentile (integrated from the MODE)", false);	cmd.add( save_PDF_info);
	  ValueArg<string> save_PDF_info_file("","save_PDF_info_file","filename for save_PDF_info option  (default: inputfile.pdfinfo) ",false,"","string"); cmd.add(save_PDF_info_file);
	  ValueArg<string> CLr_file("","CLr_file","filename for the save_CLr option (default: inputfile.Clr) ",false,"CLr","string"); cmd.add(CLr_file);
	  SwitchArg testGauss("","testGauss", "generate gaussian sample from N(0,1) and calculate 1,2,3 sigma confidence ranges.", false);	cmd.add( testGauss );

/* 	ValueArg<string> overplot_region_type("","plot_reg_type","region type for plotting on map (dot: default, emptydot) (to be connected with ml and mb parameters) ",false,"dot","string"); cmd.add(overplot_region_type); */
/* 	ValueArg<long> overplot_region_type_dot_points("","plot_reg_type_circle_points","numer of points to plot in the width of circle (default: 1) (to be connected with ml and mb parameters) ",false,1,"long"); cmd.add(overplot_region_type_dot_points); */

/* /\* 	ValueArg<string> save_overplot_as("","save_over_plot_as","saves the overplots into a txt file with pixel numbers instead of l,b ",false,"dot","string"); cmd.add(save_overplot_as); *\/ */

/* 	SwitchArg show_ecliptic("","show_ecliptic", "overplots ecliptic and poles from data/512-poles_ecliptic file", false);	cmd.add( show_ecliptic ); */
/* 	SwitchArg show_equator("","show_equator", "overplots equator and poles from data/512-poles_equator file", false);	cmd.add( show_equator ); */

/* 	ValueArg<long> Yl("","Yl","l'th harmonic to plot (default: 4)",false,4,"long"); cmd.add(Yl); */
/* 	ValueArg<long> Ym("","Ym","m'th harmonoc to plot (default: 0)",false,0,"long"); cmd.add(Ym); */
/* 	ValueArg<double> reim_ratio("","reim_ratio","ratio of re to im part in requested Ylm to plot (default: 0)",false,0,"double"); cmd.add(reim_ratio); */

	//SwitchArg pyth("","pyth", "whether or not plot the data with python script", false);	cmd.add( pyth );
	//
	// Parse the command line.
	//
	  cmd.parse(argc,argv);

	//
	// Set variables
	//
	  _input_file = input_file.getValue(); 
	  _CLr = CLr.getValue(); 
	  _doCLr = doCLr.getValue(); 

	  _save_PDF_info = save_PDF_info.getValue(); 
	  _save_PDF_info_file = save_PDF_info_file.getValue(); if (_save_PDF_info_file == "") _save_PDF_info_file = _input_file+".pdfinfo";
	  _nogaps = nogaps.getValue(); 
	  _show_info = show_info.getValue(); 
	  _save_CLr = save_CLr.getValue(); 
	  _CLr_file = CLr_file.getValue(); 
	  _CLacc = CLacc.getValue(); 
	  _modesnum = modesnum.getValue(); 
	  _testGauss=testGauss.getValue();



	// read in the CLs and check
	  vector<double>	v3;
	  v3 = CL.getValue(); if (v3.size() > 0) { for (i=0;i<(long)v3.size();i++) { _CL.addq(v3[i]);  if (_CL(i) > 1) { printf("ERROR: one of the CL values given is larger than unity.\n"); exit(1); } } }
	  


	} catch ( ArgException& e ) { cout << "ERROR: " << e.error() << " " << e.argId() << endl;  }

}



/***************************************************************************************/
void makeGaussTest(cpeds_queue<double>& X, cpeds_queue<double>& Y) {
	mscsFunction sample;
	
	double res=0.001;
	double xmin=-5,xmax=5;
	printf("Generating gassian PDF: resolution %lf, range: %lf, %lf\n",res,xmin,xmax);
	sample.mkGauss(xmin,xmax,res,1,1,1,0);

	printf("Setting CLs\n");
	_CL.delete_all_queue();
	_CL.addq(0.682689);
	_CL.addq(0.9545);
	_CL.addq(0.9973);
	_CL.printq_double();
	
	printf("The confidence ranges should be:\n");
	printf("(-1,1)\n");
	printf("(-2,2)\n");
	printf("(-3,3)\n");
	
	double* tmp=sample.extractArguments();
	X.import_array(tmp,sample.pointsCount());
	delete [] tmp;
	tmp=sample.extractValues();
	Y.import_array(tmp,sample.pointsCount());
}
