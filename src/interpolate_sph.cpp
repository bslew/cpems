/*!
 * \file interpolate_sph.cpp
 *
 * This program interpolates an input Healpix map to higher resolutions
 *
 *  Project: Mscs
 *  Created on: Nov 7, 2012 8:29:27 PM
 *  Author: blew
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <tclap/CmdLine.h>
#include <QtCore/QStringList>
#include <QtCore/QString>
#include "cpeds-msgs.h"
#include "Mscs-map.h"



#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif

string _programVersionString, _outfile, _outdir, _inputFile, _rbf;
long _ns, _tons, _seed;
double _r0,_r1;
bool _testName, _maskSH;

void parseOptions(int argc, char** argv);
string getProgramVersionString();
void makeSelfTest();
void interpolateMap();
void interpolateDirections(cpedsDirectionSet& ds);

int main(int argc, char **argv) {
	
	cpedsMsgs msgs(".interpolate_sph",false,"",High);
	msgs.setTimeSinceStartOn(true);
	msgs.setSaveRunWriteMode('a');
	msgs.say("PARSE THE COMMAND LINE ARGUMENTS",High);
	string s;
	QStringList qsl;
	QString qs;
	//----------------------------------------------------------------------------------------------------
	//PARSE THE COMMAND LINE ARGUMENTS
	//----------------------------------------------------------------------------------------------------
	parseOptions(argc,argv);


	if (_testName) {
		makeSelfTest();
		msgs.say("FINISHING",High);
		return 0;		
	}
	qs=_inputFile.c_str();
	if (qs.contains("-Tn-bin",Qt::CaseSensitive)) {
		// interpolate a map file
		interpolateMap();
	}
	if (_inputFile=="stdin") {
		string input_line;
		cpedsDirectionSet ds;
		double l,b,T;
		bool ok;
		while(cin) {
			getline(cin, input_line);
			qs=input_line.c_str();
			printf("read: %s\n",qs.toStdString().c_str());
			qsl=qs.split(' ');
			if (qsl.count()==3) { //msgs.criticalError("wrong number of data",Top);
				l=qsl[0].toDouble(&ok); if (not ok) msgs.criticalError("incorrect conversion: l",Top);
				b=qsl[1].toDouble(&ok); if (not ok) msgs.criticalError("incorrect conversion: b",Top);
				T=qsl[2].toDouble(&ok); if (not ok) msgs.criticalError("incorrect conversion: T",Top);
				ds.append(cpedsDirection(l*PI180, b*PI180, T));
			}
		}
		interpolateDirections(ds);
	}
	
	
	msgs.say("FINISHING",High);
	return 0;
}


void parseOptions(int argc, char** argv) {
	long i;
	_programVersionString=getProgramVersionString();
	try {
		
		//     CmdLine cmd("RT4\n\n. ",' ', "0.1" );
		CmdLine cmd("interpolate_sph\nTorun Centre for Astronomy, UMK, POLAND.\n "
				"This program interpolates healpix maps to at higher resolution pixels using radial function interpolation method on a tangent plane.\n\n. "
				"\n"
				"example usage: interpolate_sph"
				"",' ', _programVersionString.c_str() );
		
		// 
		// Define arguments
		//
		
		ValueArg<long> seed("", "seed", "seed to generate a rangom map (testing stage only)", false,1,"long");	cmd.add( seed );
		ValueArg<long> ns("", "ns", "healpix resolution parameter for the source map", false,2,"long");	cmd.add( ns );
		ValueArg<long> tons("", "tons", "target healpix resolution parameter", false,4,"long");	cmd.add( tons );
		ValueArg<double> r0("", "r0", "interpolation scale parameter [deg]", false,10,"double");	cmd.add( r0 );
		ValueArg<double> r1("", "r1", "interpolation range parameter. Neighboring points will be selected from within this search radii [deg]", false,20,"double");	cmd.add( r1 );
		
		SwitchArg test("","test", "perform a self-test (default: false)", false);	cmd.add( test );
		SwitchArg maskSH("","maskSH", "mask southern hemisphere before the interpolation (default: false)", false);	cmd.add( maskSH );
//		
		ValueArg<string> outfile("o","","output file name",false,"dstMap-Tn-bin","string"); cmd.add(outfile);
		ValueArg<string> outdir("","odir","output directory name",false,"","string"); cmd.add(outdir);
		
		std::vector<string> allowedStr;
		allowedStr.push_back("gaussian");
		allowedStr.push_back("multiquad");
		
		ValuesConstraint<string>* allowedStrNew;
		allowedStrNew = new ValuesConstraint<string>(allowedStr);
		ValueArg<string> rbf("", "rbf", "radial basis function name (default: gaussian)", false,"gaussian", allowedStrNew);	cmd.add( rbf );
//		allowedStr.clear();

		UnlabeledValueArg<string> input_file("input_file","data file to interpolate",false, "stdin","string");	cmd.add( input_file );
		
		cmd.parse(argc,argv);
		
		//
		// Set variables
		//
		_ns = ns.getValue();
		_tons = tons.getValue();
		_r0=r0.getValue();
		_r1=r1.getValue();
		_seed=seed.getValue();
		_testName=test.getValue();
		_inputFile=input_file.getValue();
		_maskSH=maskSH.getValue();

		_rbf=rbf.getValue();
		
		_outfile = outfile.getValue();
		_outdir = outdir.getValue(); 	if (_outdir=="") _outdir=".";
		
		
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



void makeSelfTest() {
	mscsMap srcMap("srcMap", _ns);
	mscsMap dstMap("dstMap", _tons);
	dstMap.make_equatorial_mask(20);
	
	srcMap.makekill_space_manager("make","T");
	cpedsRNG rns("gaussian_circle");
	rns.seed(_seed);
	srcMap.make_gaussian_map(0,1,-1,&rns);
	srcMap.savebinT("srcMap-Tn-bin");
	srcMap.interpolate(dstMap, _r1,_r0, "gaussian");
	dstMap.savebinT("dstMap-Tn-bin");	
}
void interpolateMap() {
	mscsMap srcMap("srcMap", _ns);
	mscsMap dstMap("dstMap", _tons);
	
	srcMap.loadbinT(_inputFile);
	srcMap.interpolate(dstMap, _r1, _r0, "gaussian");
	dstMap.savebinT(_outdir+"/"+_outfile);		
}

// ds in [rad]
void interpolateDirections(cpedsDirectionSet& ds) {
	mscsMap srcMap("srcMap", _ns);
	mscsMap dstMap("dstMap", _tons);
	srcMap.makekill_space_manager("make","T");
	srcMap.makekill_space_manager("make","m");
	dstMap.makekill_space_manager("make","m");
	srcMap.clear_map();
	dstMap.clear_map();
	srcMap.clean_mask();
	dstMap.clean_mask();
	srcMap.m()=double(0);
	srcMap.set_map_coord();
	dstMap.set_map_coord();
	long j;
	cpedsDirection n;
	for (long i = 0; i < ds.size(); i++) {
		j=srcMap.get_Ci(ds[i]);
		srcMap.set_T(j,ds[i].val());
		srcMap.set_m(j,1);
		n=ds[i]; n.lat()=-n.lat();

		j=srcMap.get_Ci(n);
		srcMap.set_T(j,ds[i].val());
		srcMap.set_m(j,1);
		
//		srcMap.make_circle_dot(PI180inv*ds[i].l(),PI180inv*ds[i].b(),_r1,1,"m","dot",-1,"");
//		dstMap.make_circle_dot(PI180inv*ds[i].l(),PI180inv*ds[i].b(),_r1,1,"m","dot",-1,"");
	}
	srcMap.savebinT(_outdir+"/"+_outfile+".src");
	srcMap.savebinm(_outdir+"/"+_outfile+".srcm");
	if (_maskSH) {
//		srcMap.make_circle_dot(0,-90,90,0,"m","dot",-1,"");
		for (long i = 0; i < dstMap.pixNum(); i++) {
			if (dstMap.get_C(i).b()<0) dstMap.set_m(i,0);
		}
//		dstMap.make_circle_dot(0,-90,90,0,"m","dot",-1,"");
	}
	dstMap.savebinm(_outdir+"/"+_outfile+".mask");
	srcMap.interpolate(dstMap, _r1, _r0, _rbf);
	dstMap.savebinT(_outdir+"/"+_outfile);		
}
