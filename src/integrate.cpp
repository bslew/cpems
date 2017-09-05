#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <tclap/CmdLine.h>
#include "cpeds-math.h"
//#include "cpeds-templates.h"
#include "Mscs-function.h"

using namespace std;
using namespace TCLAP;

void parseOptions(int argc, char** argv);
string _input_file;
double _from,_to;
bool _normalize;
bool _xminAuto,_xmaxAuto;

int main(int argc, char** argv) {
	parseOptions(argc,argv);
	
	mscsFunction f;
	
	
	//
	// load data
	//
	f.load(_input_file);
	
	//
	// check ranges
	//
	
	f.checkRanges();
	if (_xminAuto) _from=f.getMinArg();
	if (_xmaxAuto) _to=f.getMaxArg();
	
	//
	// operations before integration
	//
	if (_normalize) f.normalize();
	
	
	//
	// integrate
	//
	double I=f.integrate(_from,_to);
	
	printf("%lE\n",I);
	
	//  double x,y,I;
	//  cpeds_queue<double> X,Y;
	//  double *px,*py;
	//  FILE* f;
	//
	//
	//  parseOptions(argc,argv);
	//
	//  if (_input_file == "stdin") f=stdin; else  f=fopen(_input_file.c_str(),"r");
	//  if (f==NULL) { printf("ERROR: the input file not found. exiting... \n"); exit(0); }
	//  while (fscanf(f,"%lE %lE",&x,&y)!=EOF) { X.addq(x); Y.addq(y); }
	//  if (_input_file != "stdin") fclose(f);
	//
	//  // 
	//  // derive the integral
	//  //
	//
	//  px=X.export_array();
	//  py=Y.export_array();
	//  I=cpeds_integrate_1d(X.get_size(),px,py);
	//  printf("%lE\n",I);
	//  delete [] px;
	//  delete [] py;
	
}





void parseOptions(int argc, char** argv) {
	/*   long i; */
	/*   string::size_type j; */
	
	
	try {
		
		CmdLine cmd("integrate: Given a function in format X Y  in the input file it computes the integral", ' ', "" );
		
		// 
		// Define arguments
		//
		
		UnlabeledValueArg<string> input_file("input_file","file with input PDF",true,"","string");	cmd.add( input_file );
		SwitchArg normalize("","normalize", "normalize function before integration (useful for integrating PDF)", false);	cmd.add( normalize );
		ValueArg<double> from("f","from","integrate from argument value (0)",false,0,"double"); cmd.add(from);
		ValueArg<double> to("t","to","integrate to argument value (0)",false,0,"double"); cmd.add(to);
		
		
		//
		// Parse the command line.
		//
		cmd.parse(argc,argv);
		
		//
		// Set variables
		//
		_input_file = input_file.getValue(); 
		_xminAuto=true;
		_xmaxAuto=true;
		if (from.isSet()) { _from=from.getValue(); _xminAuto=false; }
		if (to.isSet()) { _to=to.getValue(); _xmaxAuto=false; }
		_normalize=normalize.getValue();
		
		
		
		
	} catch ( ArgException& e ) { cout << "ERROR: " << e.error() << " " << e.argId() << endl;  }
	
}



