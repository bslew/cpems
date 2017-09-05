/*!
 * \file bin_function.cpp
 *
 *  Created on: Aug 3, 2010
 *      Author: blew
 */

#include <stdlib.h>
#include <stdio.h>
#include <cpgplot.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <tclap/CmdLine.h>
#include "cpeds-consts.h"
#include "cpeds-math.h"
#include "Mscs-function.h"
#include "cpeds-msgs.h"


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
bool _plot, _eqw,_geo, _geoinv,_fromAsArg, _sort,_sortBinsDesc;
long _fsize, _colx,_coly;
double _dl,_gm, _from;
bool _testName, _calcBinAvg;
//-----------------------------------

void calcBinAvg(mscsFunction& f, double from, double d);

int main(int argc, char **argv) {
	long i;
	string outfile;
	/*   long bintabs; */
	cpedsMsgs msgs("bin_function");
	msgs.saveThisRun(argc,argv);

	
	//----------------------------------------------------------------------------------------------------
	//PARSE THE COMMAND LINE ARGUMENTS
	//----------------------------------------------------------------------------------------------------
	mscsFunction f;
	mscsFunction g;
	long from;
	double x;
	
	parseOptions(argc,argv);
	msgs.setLogFileName(_outfile+".info");
	msgs.setVerbosity(Zero);
	msgs.setLogVerbosity(High);
	
	if (_testName) {
		f.mkSin(0,1,0.001);
		f.save("sin");
		g=f;
	}
	else {
		f.load(_input_file,true,_colx,_coly);	  
	}
	if (_sort) {
//		if (_sortDesc) f.sortFunctionArgDescending();
//		else 
			f.sortFunctionArgAscending();
	}
	
	
	if (_calcBinAvg) {
		calcBinAvg(f,_from,_dl);
	}
	_fsize=f.pointsCount();
	
	if (_fromAsArg) {
		f.f(_from,&from);
		_from=from;
		x=f.getX(_from);
		printf("binning will start from: %lE, idx: %li\n",x,_from);
	}
	
	if (_bintab.size()==0){
		if (_geo) {
			double b=_dl;
			i=_from;
			while (i<_fsize) {
				i+=long(round(b));
				if (i<=_fsize) { _bintab.append(long(round(b))); } else { _bintab.append(_fsize-i+long(round(b))); }
				if (_geoinv) {
					b/=_gm;
					if (b<1) b=1;
				}
				else
					b*=_gm;
			}
		}
		else {
			i=_from;
			while (i<_fsize) {
				if (_fromAsArg) {
					x+=_dl;
					f.f(x,&from);
					if (from<_fsize-1) { _bintab.append(from-i); i=from; }
					else { _bintab.append(_fsize-from); i=_fsize; }
					
					//				  printf("binning will start from: %lE, idx: %li\n",x,i);
				}
				else {
					i+=_dl;				  
					if (i<=_fsize) _bintab.append(_dl); 
					else _bintab.append(_fsize-i+_dl);
				}
			}		  
		}
	}
	
	
	if (_sortBinsDesc) {
		mscsFunction tmpf;
		for (unsigned long i = 0; i < _bintab.size(); i++) {
			tmpf.newPoint(_bintab[i],_bintab[i]);
		}
		tmpf.sortFunctionArgDescending();
		for (unsigned long i = 0; i < _bintab.size(); i++) {
			_bintab[i]=tmpf.getX(i);
		}
	}
	
	cpedsList<double> varInBin;
	cpedsList<double> alpha;
	// select quantile values to calculate
	alpha << 0.01/2 << 1.0-0.01/2; // 99% CL
	alpha << 0.0455003/2 << 1.0-0.0455003/2 ; // 95% CL
	alpha << 0.317311/2 <<  1.0-0.317311/2; // 68% CL
	alpha << 0.5 ; // 2nd quartile
	
	vector< cpedsList<double> > stats;
	double xmin,xmax;
	f.checkRanges();
	xmin=f.getMinArg();
	xmax=f.getMaxArg();
	printf("xmin: %lE, xmax: %lE\n",xmin,xmax);
	f.setVerbosityLevel(High);
	f=f.binFunctionLin2(_from,_bintab,_w,&varInBin,&alpha,&stats);
	
	
	if (_outfile!="") {
		f.save(_outfile);
		_bintab.save(_outfile+".bins",false,"double");	  
		//	  varInBin.save(_outfile+".var",false,"double",10);
		matrix<double> m(_bintab.size(),7+alpha.size()+8);
		for (long i = 0; i < m.RowNo(); i++) {
			m(i,0)=f.getx(i); // binned function arg
			m(i,1)=f.f(i); // binned function value
			if (i==0) {
				m(i,2)=(f.getx(i)-xmin); // lower bin half width
				m(i,3)=(f.getx(i+1)-f.getx(i))/2; // upper bin half width
			}
			else {
				if (i==m.RowNo()-1) {
					m(i,2)=(f.getx(i)-f.getx(i-1))/2; // lower bin half width
					m(i,3)=(xmax-f.getx(i)); // upper bin half width
				}
				else {
					m(i,2)=(f.getx(i)-f.getx(i-1))/2; // lower bin half width
					m(i,3)=(f.getx(i+1)-f.getx(i))/2; // upper bin half width
				}
			}
			m(i,4)=sqrt(varInBin[i]/2); // st.dev. in bin
			m(i,5)=sqrt(varInBin[i]/2); // st.dev. in bin
			m(i,6)=_bintab[i];  // number of binned data points

			//
			// append extra columns from stats structure
			//
			m(i,7)=stats[0][i]; // mean x 
			m(i,8)=stats[1][i]; // mean y 
			m(i,9)=stats[2][i]; // sigma x 
			m(i,10)=stats[3][i]; // sigma y
			m(i,11)=stats[4][i]; // Distance from the binned argument to lower bin boundary (should be the same as col 2)
			m(i,12)=stats[5][i]; // Distance from the binned argument to upper bin boundary (should be the same as col 3)
			m(i,13)=stats[6][i]; // Distance from the binned value to the minimal data value 
			m(i,14)=stats[7][i]; // Distance from the binned value to the maximal data value 
			
			for (long j = 0; j < alpha.size(); j++) {
				m(i,15+j)=stats[8+j][i]; // quantile function value for a given alpha threshold
			}
			
		}
		cpeds_matrix_save(m,_outfile+".all");
		
		msgs.setLogFileName(_outfile+".all.info");
		msgs.loggingOn();
		msgs.say(".all matrix description:",High);
		msgs.say("col0: binned function argument (according to weights)",Medium);
		msgs.say("col1: binned function value (according to weights)",Medium);
		msgs.say("col2: disntace from the binned argument to lower bin boundary",Medium);
		msgs.say("col3: disntace from the binned argument to upper bin boundary",Medium);
		msgs.say("col4: variance in the bin",Medium);
		msgs.say("col5: variance in the bin",Medium);
		msgs.say("col6: number of binned data points",Medium);
		msgs.say("col7: mean x",Medium);
		msgs.say("col8: mean y",Medium);
		msgs.say("col9: sigma x",Medium);
		msgs.say("col10: sigma y",Medium);
		msgs.say("col11: Distance from the binned argument to lower bin boundary (should be the same as col 2)",Medium);
		msgs.say("col12: Distance from the binned argument to upper bin boundary (should be the same as col 3)",Medium);
		msgs.say("col13: Distance from the binned value to the minimal data value ",Medium);
		msgs.say("col14: Distance from the binned value to the maximal data value ",Medium);
		for (long j = 0; j < alpha.size(); j++) {
			msgs.say("col%.0lf: a bound of two-sided confidence interval corresponding to a given alpha threshold (%lf)",double(j+15),alpha[j],Medium);
		}
		msgs.say("The quoted thresholds correspond to two-sided confidence intervals, so e.g. a 2sigma CR"
				"corresponds to the threshold value erfc(2/sqrt(2))/2 = 0.0455/2 = 0.02275",Medium);
//		msgs.say("The bounds are specified such that the ")
		msgs.loggingOff();
//		msgs.setLogFileName(_outfile+".stats.info");
	}
	else {
		f.print();
	}
	
	//  if (_test) {
	//	  double* xin=g.extractArguments();
	//	  double* yin=g.extractValues();
	//	  long* binSize=_bintab.toCarray();
	//	  double* xout;
	//	  double* yout;
	//	  long Nout;
	//	  yout=cpeds_bin_function(xin,yin,NULL,g.pointsCount(),binSize,_bintab.count(),&xout,&Nout,_from);
	//	  mscsFunction test("test",xout,yout,Nout);
	//	  test.save("test.txt");
	//	  delete [] xin;
	//	  delete [] yin;
	//	  delete [] xout;
	//	  delete [] yout;
	//	  delete [] binSize;
	//  }
}


void parseOptions(int argc, char** argv) {
	long i;
	try {
		
		CmdLine cmd("bin_function: program for binning the input function read from the 2-column text file", ' ', "");
		
		UnlabeledValueArg<string> file("files", "input binary file name (prefix) in Mscs format ",false,"", "string");     	cmd.add( file );
		ValueArg<string> outfile("o","outfile","name of the output file (postfix)",false,"","string"); cmd.add(outfile);
		//	ValueArg<long> from("f","from","bin from l (default: 10)",false,0,"long"); cmd.add(from);
		ValueArg<double> from("f","from","bin from l (default: 10)",false,0,"double"); cmd.add(from); // BLmodification (Jul 14, 2011, 12:28:08 PM): this is to allow using this parameter also when option fromAsArg is true
		ValueArg<long> colx("x","colX","x column number, starting from 0 (default: 0)",false,0,"long"); cmd.add(colx);
		ValueArg<long> coly("y","colY","y column number, starting from 0 (default: 1)",false,1,"long"); cmd.add(coly);
		ValueArg<long> dl("d","dl","bin size in l (default: 10)",false,10,"long"); cmd.add(dl);
		MultiArg<long> b("b","bin","list of arbitrary bins to use for binning (if defined, the f,t,d args will be ignored)",false,"long"); cmd.add(b);
		/* 	ValueArg<long> nsigma("o","output","",false,4,"long"); cmd.add(nsigma); */
		SwitchArg eqw("","eqw","equal weights for binning (default: false). NOT IMPLEMENTED YET",false); cmd.add(eqw);
		SwitchArg geo("","geo","make bins in geometric sequence using multiplier gm and dl (d) parameters (default: false).",false); cmd.add(geo);
		SwitchArg geoinv("","geoinv","make bins in decreading geometric sequence using dividor gm and dl (d) parameters (default: false).",false); cmd.add(geoinv);
		SwitchArg sort("","sort","sort the input data first in ascending order. "
				"This option should be used if he data are not sorted. (default: false).",false); cmd.add(sort);
		SwitchArg sortBinsDesc("","sortBinsDesc","sort generated bins in descending order before using "
				"them for binning. (default: false).",false); cmd.add(sortBinsDesc);
		ValueArg<double> gm("","gm","geometric sequence multiplied (default: 2)",false,2,"double"); cmd.add(gm);
		SwitchArg test("","test","test the program (default: false)",false); cmd.add(test);
		SwitchArg fromAsArg("","fromAsArg","treat from argument as argument of the function, and not as index of the point at which the binning starts (default: false). NOT IMPLEMENTED YET",false); cmd.add(fromAsArg);
		
		SwitchArg calcBinAvg("","calcBinAvg","calculate and print out average value within a bin defined with -f and -d options (default: false). ",false); cmd.add(calcBinAvg);
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
		_geoinv = geoinv.getValue();   if (_geoinv) _geo=true;
		_sort=sort.getValue();
		_sortBinsDesc=sortBinsDesc.getValue();
		_gm=gm.getValue();
		_testName=test.getValue();
		_calcBinAvg=calcBinAvg.getValue();
		//	_w.append(0);
		
		//	_interpolate = interpolate.getValue();
		
		_fromAsArg=fromAsArg.getValue();
		_colx=colx.getValue();
		_coly=coly.getValue();
		
		vector<long> bins = b.getValue();
		if (bins.size()>0) {
			for (i=0;i<bins.size();i++) { _bintab.append(bins[i]); }
		}
		
		
	} catch ( ArgException& e )
	{ cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
}



/***************************************************************************************/
void calcBinAvg(mscsFunction& f, double from, double d) {
	double st,en;
	st=from;
	en=from+d;
	printf("%lE\n",f.cut(st,en).meanf());
	exit(0);
}
