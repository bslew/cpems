#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <tclap/CmdLine.h>
#include "cpeds-math.h"
#include "cpeds-templates.h"

#ifndef _NO_NAMESPACE
using namespace std;
//using namespace math;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif

void print_info();
void print_warning_and_exit();
cpeds_queue<double>* interpolate_stuff(cpeds_queue<double>* qx, cpeds_queue<double>* qy, cpeds_queue<double>* arg, cpeds_queue<long>* from, cpeds_queue<long>* to, cpeds_queue<long>* type);

void parseOptions(int argc, char** argv);
string _input_file,_output_file,_xargs,_conf;
double _from,_to,_step;
bool _pdf,_from_conf,_from_args;
long fromprev,toprev;
long _pointsNum, _interpolationType;

int main(int argc, char** argv) {
  double x,y;
  cpeds_queue<double> qx("x"),qy("y"),arg("arg"),*Qy; // 
  cpeds_queue<long> from("from"),to("to"),type("type"); // interpolation configuration
  cpeds_queue<double> tmpx,tmpy,tmparg,*tmpQy; 

  FILE* f;

  long i,j,k,i1,i2,j1,j2;
  bool pdfOK,rangeOK;

  parseOptions(argc,argv);

  //
  // read in the function to interpolate
  //
  printf("* reading the function to interpolate...\n");
  if (_input_file == "stdin") f=stdin; else  f=fopen(_input_file.c_str(),"r");
  if (f==NULL) { printf("ERROR: the input file not found. exiting... \n"); exit(0); }
  while (fscanf(f,"%lE %lE",&x,&y)!=EOF) { qx.addq(x); qy.addq(y); }
  if (_input_file != "stdin") fclose(f);


  //
  // set up the interpolation arguments
  // 
  if (_from == _to) {
	  _from=qx(0);
	  _to=qx(qx.get_size()-1);
  }
  if (_step==0 and _pointsNum==0) {
	  printf("ERROR in configuration. Set either the step or the number of points\n");
	  exit(0);
  }
  else {
	  if (_step==0) {
		  _step=(_to-_from)/_pointsNum;
	  }
  }
  
  printf("* Configuring interpolation arguments...\n");
  if (_from_args) arg.loadq(_xargs,"double");
  else { for (x=_from;x<=_to;x+=_step) { arg.addq(x); } }
  if (arg(arg.get_size()-1) < _to) arg.addq(_to);

  //
  // setup the configuration for the interpolation
  //
  printf("* Configuring the interpolation type...\n");
  if (_from_conf) {
    f=fopen(_conf.c_str(),"r");
    while (fscanf(f,"%li %li %li",&i,&j,&k)!= EOF) { from.addq(i); to.addq(j); type.addq(k); }
    fclose(f);
  }
  else {
    from.addq(0L); to.addq(qx.get_size()-1); type.addq(_interpolationType);
  }

/*   qx.printq_double(); */
/*   qy.printq_double(); */
/*   arg.printq_double(); */

/*   from.printq_long(); */
/*   to.printq_long(); */
/*   type.printq_long(); */
    


  //
  // PIECE-WISE INTERPOLATE
  //

  Qy=interpolate_stuff(&qx,&qy,&arg,&from,&to,&type);



  //
  // PDF part
  //
  // this checks if there are negative values in the interpolated function, and if so replaces that part of the interpolated function with a new interpolation
  // done only for the concerned range. this is for special use only. eg. it doesn't make sense if the input function has any negative values
  if (_pdf) {
    pdfOK = Qy->all_nonnegative();
/*     rangeOK=true; */

/*     while ( pdfOK==false || rangeOK==true ) { */
    while ( pdfOK==false ) {
      //
      // copy the first bad part onto a new queue
      //
      i=Qy->get_first_negative_idx();
      i1=qx.find_closest(arg(i));

      if (qx(i1) < arg(i)) i2=i1+1; else { i1--; i2=i1+1; }


      // set interpolation function
      tmpx.delete_all_queue(); tmpy.delete_all_queue();
      tmpx.addq(qx(i1));      tmpx.addq(qx(i2));
      tmpy.addq(qy(i1));      tmpy.addq(qy(i2));

      // set interpolation configuration
/*       fromprev=from(from.get_size()-1); */
/*       toprev=to(to.get_size()-1); */
/*       printf("fromprev %li toprev %li\n",fromprev,toprev); */

      from.delete_all_queue();      to.delete_all_queue();      type.delete_all_queue();
      from.addq(0L);       to.addq(1L);  type.addq(1L);      
/*       printf("from %li to %li\n",from(0),to(0)); */
/*       if (fromprev  != from(0) && toprev!=to(0)) rangeOK=true; else rangeOK=false; */
      
      // set interpolation arguments
      j1=arg.find_closest(qx(i1)); if (arg(j1)<qx(i1)) j1++;
      j2=arg.find_closest(qx(i2)); if (arg(j2)>qx(i2)) j2--;
      tmparg.delete_all_queue();
      for (j=j1;j<=j2;j++) { tmparg.addq(arg(j)); }

      // interpolate this bit
      tmpQy=interpolate_stuff(&tmpx,&tmpy,&tmparg,&from,&to,&type);

      // copy the results
/*       if (rangeOK) */
	for (j=j1;j<=j2;j++) { Qy->setq(j,(*tmpQy)(j-j1)); }
/*       else */
/* 	for (j=j1;j<=j2;j++) { Qy->setq(j,0); } */

      pdfOK = Qy->all_nonnegative();
/*       if (pdfOK) printf("pdfok\n"); else printf("printf pdfnot ok\n"); */
    }

    Qy->zero_all_below(0);
    
  }

/* SAVE THE INTERPOLATION  */  
  printf("* Saving results to: %s\n",_output_file.c_str());
  if (_output_file != "") {
    f=fopen(_output_file.c_str(),"w");
  
    for (i=0;i<Qy->get_size();i++) { fprintf(f,"%lE %lE\n",arg(i),(*Qy)(i));    }
    fclose(f);
  }
  else {
    for (i=0;i<Qy->get_size();i++) { printf("%lE %lE\n",arg(i),(*Qy)(i));    }
  }


  delete Qy;  
  return 0;
}

void print_info() {
  printf("USAGE: cpeds-interpolate data_file arg_file out_file [ negative: 0 - don't want negative values we're doing PDF, 1 - I don't care, just do it] [confing file] \n\n");
  printf("This is a test program for testing GSL interpolation implemented in CPEDS library\n");
  printf("The program requires two inputd date files that can be given in the command line.\n");

  printf("The first one is the file with numbers defining pairs of X and Y data values used to interpolate on (default: interpolate-data).\n");
  printf("The second file contains X arguments used to derive the interpolated value (default: interpolate-onargs).\n");
  printf("If no files are given at the command line, the default filesnames are used.\n");
  printf("The result is by default stored on the interpolate-results file.\n\n");
  printf("The arguments in the data_file are sorted in the rising order (The y values aren't).\n");
  printf("The arguments in the arg_file need to be within the rage covered by the arguments given in the data_file.\n\n\n");
  printf("The optional config file contains information on how the interpolation should be perfrormed. Ea ch row of the\n");
  printf("file should be as follows\n");
  printf("0 j type\n");
  printf("j k type\n");
  printf("where i and j are from and to indexes of the points to be interpolated and type is the type of the interpolation to be used: 1 - linear, 3 - cubic");
  printf("If no file is given then default interpolation - cubic - is applied to all data points\n");

  printf("Using the default values.\n");

}



void parseOptions(int argc, char** argv) {
/*   long i; */
/*   string::size_type j; */


	try {

	  CmdLine cmd("integrate: Given a function in format X Y  in the input file it computes the integral using different methods  from the cpeds package", ' ', "" );

	// 
	// Define arguments
	//

//	  UnlabeledValueArg<string> input_file("input_file","file with input PDF","","string");	cmd.add( input_file );
	  UnlabeledValueArg<string> input_file("input_file","file with input PDF",true,"","string");	cmd.add( input_file );
	  
	  ValueArg<string> output("o","out","output file",false,"","string"); cmd.add(output);
	  ValueArg<double> from("f", "from", "interpolate from", false,0,"double");	cmd.add( from );
	  ValueArg<double> to("t", "to", "interpolate to", false,0,"double");	cmd.add( to );
	  ValueArg<double> step("s", "step", "interpolate with step ", false,0,"double");	cmd.add( step );
	  ValueArg<long> pointsNum("n", "", "number of interpolated points (defined alternatively to step) ", false,0,"long");	cmd.add( pointsNum );
	  ValueArg<long> interpolationType("", "tid", "interpolation type id: 1- linear, 3 - cspline, 33 - akima", false,3,"long");	cmd.add( interpolationType );

	  ValueArg<string> xargs("","xargs","file name with the x argumets data to interpolate on if the f,t,s parameters were not used",false,"","string"); cmd.add(xargs);
	  ValueArg<string> conf("","conf","interpolation config file with defiend scheme for this interpolation",false,"","string"); cmd.add(conf);
	  SwitchArg pdf("","pdf", "triggers the interpolation mode assuming that the function is the PDF and we don't want any oscilations about zero etc. (default: false)", false);	cmd.add( pdf );
	//
	// Parse the command line.
	//
	  cmd.parse(argc,argv);

	//
	// Set variables
	//
	  _input_file = input_file.getValue(); 
	  _output_file = output.getValue(); 

	  _from = from.getValue(); 
	  _to = to.getValue(); 
	  _step = step.getValue(); 
	  _xargs = xargs.getValue(); if (_conf != "") _from_args=true;  
	  _conf = conf.getValue();  if (_conf != "") _from_conf=true;
	  _pdf = pdf.getValue(); 
	  _pointsNum= pointsNum.getValue(); 
	  _interpolationType=interpolationType.getValue();



	  


	} catch ( ArgException& e ) { cout << "ERROR: " << e.error() << " " << e.argId() << endl;  }

}


// returns an address to the list object with interpolated function values defined by qx,qy on arguments arg configured as def. by from, to, type
cpeds_queue<double>* interpolate_stuff(cpeds_queue<double>* qx, cpeds_queue<double>* qy, cpeds_queue<double>* arg, cpeds_queue<long>* from, cpeds_queue<long>* to, cpeds_queue<long>* type) {

  double *X,*Y,*Xint,*Yint,*xx,*yy,*xxint,*yyint;
  long N,Nint,i,j,l,n,nint;
  string interp;
  cpeds_queue<double>* Yq = new cpeds_queue<double>;

  X=qx->export_array();
  Y=qy->export_array();
  N=qx->get_size();

  Xint=arg->export_array();
  Yint = new double[(*arg).get_size()];
  Nint=arg->get_size();

  for (i=0;i<from->get_size();i++) { // loop over ranges of interpolation
    printf("* Interpolating from %lE to %lE\n", X[(*from)(i)], X[(*to)(i)]);
    // prepare arrays for this piece
    n=(*to)(i)-(*from)(i)+1; //printf("from %li to %li, n=%li\n",from(i),to(i),n);
    xx=new double[n]; for (j=(*from)(i);j<=(*to)(i);j++) xx[j-(*from)(i)]=X[j];
    yy=new double[n]; for (j=(*from)(i);j<=(*to)(i);j++) { yy[j-(*from)(i)]=Y[j]; } //printf("X:%lE Y:%lE\n",xx[j-from(i)],yy[j-from(i)]); }
    nint=0;
    for (j=0;j<Nint;j++) if (Xint[j]>=X[(*from)(i)] && Xint[j]<=X[(*to)(i)]) nint++;
    xxint = new double[nint]; l=0; for (j=0;j<Nint;j++) if (Xint[j]>=X[(*from)(i)] && Xint[j]<=X[(*to)(i)]) { xxint[l]=Xint[j]; /* printf("xxint[l]=%lE\n",xxint[l]); */ l++; }
  
    if ((*type)(i)==1) interp="linear";
    if ((*type)(i)==3) interp="cspline";
    if ((*type)(i)==33) interp="akima";
    
    // interpolate
    yyint=cpeds_interpolate(xx,yy,n,xxint,nint,interp.c_str());


    // copy the results onto the output array
    l=0; for (j=0;j<Nint;j++) if (Xint[j]>=X[(*from)(i)] && Xint[j]<=X[(*to)(i)]) { Yint[j]=yyint[l]; /* printf("%lE, ",yyint[l]); */ l++; }
    
    // release memory
    delete [] xx;
    delete [] yy;
    delete [] xxint;
    delete [] yyint;
  }

  for (j=0;j<Nint;j++) { Yq->addq(Yint[j]); }
/*   arg->printq_double(); */
/*   Yq->printq_double(); */

/* CLEAR THE MEMORY */
  printf("* Clearing up\n");

  delete [] X;
  delete [] Y;
  delete [] Xint;
  delete [] Yint;

  return Yq;

}
