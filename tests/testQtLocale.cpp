/*!
 * \file testQtLocale.cpp
 *
 *  Project: cpeds
 *  Created on: Jan 14, 2013 10:39:50 AM
 *  Author: blew
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <tclap/CmdLine.h>
#include <QtCore/QStringList>
#include <QtCore/QString>
#include "cpeds-msgs.h"
#include "cpeds-consts.h"
#include "QtCore/QLocale"


#ifndef _NO_NAMESPACE
using namespace std;
//using namespace math;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif

string _programVersionString;

void parseOptions(int argc, char** argv);
string getProgramVersionString();

int main(int argc, char **argv) {
	
	cpedsMsgs msgs(".testQtLocale",false,"",High);
	msgs.setSaveRunWriteMode('a');
	string s;
	//----------------------------------------------------------------------------------------------------
	//PARSE THE COMMAND LINE ARGUMENTS
	//----------------------------------------------------------------------------------------------------
	parseOptions(argc,argv);


	 QLocale polish(QLocale::Polish, QLocale::Poland);
	 QLocale us(QLocale::English, QLocale::UnitedStates);
	 QLocale::setDefault(polish);
	 QString s1 = polish.toString(1.571429E+07, 'e');
	 QString s2 = us.toString(1.571429E+07, 'e');
	 printf("s1:%s\n",s1.toStdString().c_str());
	 printf("s2:%s\n",s2.toStdString().c_str());

//	 double d = egyptian.toDouble(s1);
//	 int i = egyptian.toInt(s2);
	 
	QString tmpstr(".");
//	printf("tmpstr: %s \n",tmpstr.toStdString().c_str());
	tmpstr.sprintf("to jest liczba: %.10lE\n",PI);
	printf("tmpstr: %s \n",tmpstr.toStdString().c_str());
	
	
	
	return 0;
}


void parseOptions(int argc, char** argv) {
	long i;
	_programVersionString=getProgramVersionString();
	try {
		
		//     CmdLine cmd("RT4\n\n. ",' ', "0.1" );
		CmdLine cmd("testQtLocale\nTorun Centre for Astronomy, UMK, POLAND.\n "
				"This program ...\n\n. "
				"\n"
				"example usage: testQtLocale"
				"",' ', _programVersionString.c_str() );
		//     CmdLine cmd("RT4 azimutal scan simulator for OCRA-p\n", 'This program simulates the RT4 scan motion, sky rotation during \n\n. ', 0.1 );
		
		// 
		// Define arguments
		//
		
//		ValueArg<double> Ndec("", "Ndec", "number of periods in Lissajous trajectory [21]", false,21,"double");	cmd.add( Ndec );
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

		//     UnlabeledValueArg<string> input_file("input_file","file to plot","nofile","string");	cmd.add( input_file );
		//     ValueArg<string> mask("m","mask","wheather to mask the file before doing statistics (prefix)",false,"","string"); cmd.add(mask);
		//     SwitchArg mask_from_here("","MM", "take mask from current directory rather than from the default directory", false);	cmd.add( mask_from_here );
		//     SwitchArg dont_save_mask("","dont_save_mask", "do not print to file masked pixels; useful for -Tn-txt savings", false);	cmd.add( dont_save_mask );
		//     ValueArg<string> proj("p","proj","projection to use [mall - Mollweide/ aitoff etc] (default: mall)",false,"mall","string"); cmd.add(proj);
		//     ValueArg<string> ft("f","ft","input file type [bin, txt] (default: bin)",false,"bin","string"); cmd.add(ft);
		//     ValueArg<string> ord("","ord","ordering of the map [ring, nest] default: nest",false,"nest","string"); cmd.add(ord);
		//     SwitchArg mink("M","Mink", "plot also the minkowski functionals for the map", false);	cmd.add( mink );
		//     ValueArg<long> mlevels("", "mlev", "number of levels to calculate the minkowski functionals (default: 100)", false,100,"long");	cmd.add( mlevels );
		
		cmd.parse(argc,argv);
		
		//
		// Set variables
		//
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
	QStringList qrevl=qrev.split(' ');
	if (qrevl.size()==3) rev=qrevl[1].toStdString();
#else
	rev="0.1.1";
#endif
	
#ifdef GOMPI
	rev+=" (MPI)";
#endif

	return rev;
}
