/*!
 * \file traffic_lights_mcmc.cpp
 *
 *  Project: cpeds
 *  Created on: Feb 12, 2017 9:20:52 PM
 *  Author: blew
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <tclap/CmdLine.h>
#include "qtbase.h"
#include "cpeds-msgs.h"
#include "cpeds-list.h"



#ifndef _NO_NAMESPACE
using namespace std;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif

string _programVersionString;

void parseOptions(int argc, char** argv);
string getProgramVersionString();
long _ns,_nm,_tMax,_N,_M;
double _Pred;
string _f;
bool _debug;

cpedsList<double> getRiskFactors(string fs) {
	cpedsList<double> f;
	QString qf=fs.c_str();
	QStringList l=qf.split(",");
	for (long i = 0; i < l.size(); i++) {
		f.append(l[i].toDouble());
	}
	
	if (f.size()<_ns) {
		for (long i = f.size(); i < _ns; i++) {
			f.append(f.last());
		}
	}
	//	printf("assuming factors:\n");
	//	f.print();
	
	return f;
}


int main(int argc, char **argv) {
	
	cpedsMsgs msgs(".traffic_lights_mcmc",false,"",High);
	msgs.setSaveRunWriteMode('a');
	string s;
	//----------------------------------------------------------------------------------------------------
	//PARSE THE COMMAND LINE ARGUMENTS
	//----------------------------------------------------------------------------------------------------
	parseOptions(argc,argv);
	
	double* lights;
	double* counter;
	cpedsList<double> mean_wait_timeTOT;
	
	for (long m = 0; m < _M; m++) {
		
		counter=cpeds_random_uniform_numbersGSL(0.5,2*_tMax+0.5,_N*_ns,0,false); // values >80 mean green light
		//	cpeds_save_matrix(counters,_N*_ns,1,"counters");
		
		cpedsList<double> mean_wait_time;
		cpedsList<double> mean_Nwands;
		cpedsList<double> mean_Nred;
		cpedsList<double> mean_Ngreen;
		
		cpedsList<double> f=getRiskFactors(_f);
		
		double EXt=20.25;
		long Nwands;
		long idx,idx_c=0;
		double c;
		double r;
		double wait_time=0;
		double Nred,Ngreen;
		double n1=0,n2=0,n3=0,n4=0;
		for (long i = 0; i < _N; i++) {
			
			Nwands=_nm;
			wait_time=0;
			Nred=0;
			Ngreen=0;
			
			for (long is = 0; is < _ns; is++) {
				idx=i*_ns+is;
				// signal lights >0.5 means red
				//			r=lights[idx]; 
				
				c=round(counter[idx]);
				if (c<=_tMax) { // if red
					Nred++;
					if (Nwands>=_ns-is) { // use magic on every red light if we have enough wands and don't wait 
						Nwands--;
						n4++;
					}
					else {
						
						// look at the counter value
						if (Nwands>0) { // it we have magic wands
							if (c>f[is]*EXt) { // use magic wand now, and don't wait at all
								Nwands--;
								n3++;
							}
							else { // save magic wand for latter use, and wait time c.
								wait_time+=c;					
								n2++;
							}
						}
						else { // if not, then simply wait time c
							wait_time+=c;
							n1++;
						}
					}
				}
				else { // if green, we don't wait
					Ngreen++;
				}
			}
			//		cout << wait_time << "\n";
			mean_wait_time.append(wait_time);
			mean_Nwands.append(Nwands);
			mean_Nred.append(Nred);
			mean_Ngreen.append(Ngreen);
		}
		
		cout << m << " EX(t) = " << mean_wait_time.mean() << endl;
		if (_debug) {
			cout << "mean N wands = " << mean_Nwands.mean() << endl;
			cout << "mean N red = " << mean_Nred.mean() << endl;
			cout << "mean N green = " << mean_Ngreen.mean() << endl;
			cout << "n1 = " << n1 << endl;
			cout << "n2 = " << n2 << endl;
			cout << "n3 = " << n3 << endl;
			cout << "n4 = " << n4 << endl;
			cout << "n3+n4 = " << n3+n4 << endl;
			
			cout << "n1/n2 = " << n1/n2 << endl;
		}
		delete [] counter;
		mean_wait_timeTOT.append(mean_wait_time.mean());
		cout << "EX(t) all = " << mean_wait_timeTOT.mean() << endl;
	}
	return 0;
}

/***************************************************************************************/
void parseOptions(int argc, char** argv) {
	long i;
	_programVersionString=getProgramVersionString();
	try {
		
		//     CmdLine cmd("RT4\n\n. ",' ', "0.1" );
		CmdLine cmd("traffic_lights_mcmc\nTorun Centre for Astronomy, UMK, POLAND.\n "
				"This program ...\n\n. "
				"\n"
				"example usage: traffic_lights_mcmc"
				"",' ', _programVersionString.c_str() );
		//     CmdLine cmd("RT4 azimutal scan simulator for OCRA-p\n", 'This program simulates the RT4 scan motion, sky rotation during \n\n. ', 0.1 );
		
		// 
		// Define arguments
		//
		
		ValueArg<long> M("", "M", "repeat (1)", false,1,"long");	cmd.add( M );
		ValueArg<long> N("", "N", "number of trials (1000)", false,1000,"long");	cmd.add( N );
		ValueArg<long> ns("", "ns", "number of traffic lights (2)", false,2,"long");	cmd.add( ns );
		ValueArg<long> nm("", "nm", "number of magic wands (1)", false,1,"long");	cmd.add( nm );
		ValueArg<double> Pred("", "Pred", "red lights probability (0.5)", false,0.5,"double");	cmd.add( Pred );
		ValueArg<string> f("f", "", "coma separated risk factors for each lights (1)", false,"1","string");	cmd.add( f );
		ValueArg<long> tMax("", "tMax", "Maximal waiting time [s] (80)", false,80,"long");	cmd.add( tMax );
		//		ValueArg<long> nLissajous("", "nLj", "number of timesteps to take for Lissajous trajectory (5000)", false,5000,"long");	cmd.add( nLissajous );
		
		//		SwitchArg horizTrack("","trackHoriz", "follow horizontal sky motion (default: true)", false);	cmd.add( horizTrack );
		
		//		ValueArg<string> outfile("o","","output file name",false,"scan-trajectory.txt","string"); cmd.add(outfile);
		//		ValueArg<string> outdir("","odir","output directory name",false,"./","string"); cmd.add(outdir);
		
		//		std::vector<string> allowedStr;
		//		allowedStr.push_back("none");
		//		allowedStr.push_back("incAIncH");
		//		allowedStr.push_back("incADecH");
		
		//		ValuesConstraint<string>* allowedStrNew;
		//		allowedStrNew = new ValuesConstraint<string>(allowedStr);
		//		ValueArg<string> traj("", "traj", "trajectory type:", false,"none", allowedStrNew);	cmd.add( traj );
		//		allowedStr.clear();
		
		//     UnlabeledValueArg<string> input_file("input_file","file to plot",true,"nofile","string");	cmd.add( input_file );
		//     ValueArg<string> mask("m","mask","wheather to mask the file before doing statistics (prefix)",false,"","string"); cmd.add(mask);
		//     SwitchArg mask_from_here("","MM", "take mask from current directory rather than from the default directory", false);	cmd.add( mask_from_here );
		//     SwitchArg dont_save_mask("","dont_save_mask", "do not print to file masked pixels; useful for -Tn-txt savings", false);	cmd.add( dont_save_mask );
		//     ValueArg<string> proj("p","proj","projection to use [mall - Mollweide/ aitoff etc] (default: mall)",false,"mall","string"); cmd.add(proj);
		//     ValueArg<string> ft("f","ft","input file type [bin, txt] (default: bin)",false,"bin","string"); cmd.add(ft);
		//     ValueArg<string> ord("","ord","ordering of the map [ring, nest] default: nest",false,"nest","string"); cmd.add(ord);
		SwitchArg debug("","debug", "pring extra info", false);	cmd.add( debug );
		//     ValueArg<long> mlevels("", "mlev", "number of levels to calculate the minkowski functionals (default: 100)", false,100,"long");	cmd.add( mlevels );
		
		cmd.parse(argc,argv);
		
		//
		// Set variables
		//
		_M=M.getValue();
		_N=N.getValue();
		_ns=ns.getValue();
		_nm=nm.getValue();
		_Pred=Pred.getValue();
		_tMax=tMax.getValue();
		_f=f.getValue();
		_debug=debug.getValue();
		
		//		_Ndec = Ndec.getValue();
		//		_nLissajous=nLissajous.getValue();
		//		_trackHoriz = horizTrack.getValue();
		
		//		_outfile = outfile.getValue();
		//		_outdir = outdir.getValue(); 	if (_outdir=="") _outdir=".";
		
		
	} catch ( ArgException& e )
	{ cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
	
	
}

string getProgramVersionString() {
	string rev;
#ifdef GENERATE_CVS_VERSION_STRING
	rev="$Revision: 1.1 $";
	QString qrev=rev.c_str();
	QStringList qrevl=qrev.split(" ");
	if (qrevl.size()==3) rev=qrevl[1].toStdString();
	else rev="could not retrieve version number";
#else
	rev="0.1.1";
#endif
	
#ifdef GOMPI
	rev+=" (MPI)";
#endif
	
	return rev;
}
