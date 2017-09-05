/*!
 * \file test-cpeds-smooth.cpp
 *
 *  Project: cpeds
 *  Created on: Jan 4, 2013 9:22:55 AM
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
#include "cpeds-direction_set.h"
#include "cpeds-project.h"
#include "cpeds-smooth.h"
#include "cpeds-math.h"
#include "Mscs-function3dregc.h"



#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif

string _programVersionString;
bool _populateTest;

void populateTest();
void parseOptions(int argc, char** argv);
string getProgramVersionString();

int main(int argc, char **argv) {
	
	cpedsMsgs msgs(".test-cpeds-smooth",false,"",High);
	msgs.setSaveRunWriteMode('a');
	string s;
	//----------------------------------------------------------------------------------------------------
	//PARSE THE COMMAND LINE ARGUMENTS
	//----------------------------------------------------------------------------------------------------
	parseOptions(argc,argv);

	if (_populateTest) { populateTest(); exit(0); }
	
	cpedsDirectionSet ds;
	cpedsList<double> v;
	ds.load("corrections2.ds");
	for (long i = 0; i < ds.size(); i++) {		v.append(ds[i].value());	}
	double m=v.mean();
	printf("the mean is: %lf\n", m );
	for (long i = 0; i < ds.size(); i++) {	ds[i].setVal(ds[i].value()-m);	}
	
	double res=0.002*PI180;
	double fwhm=0.5*PI180;
	matrix<double> smo;
	/*
	  //Jan 4, 2013 3:29:37 PM - DEBUG BEGIN there is something wrong with this smoothing :(
	  
	  cpedsSmooth smoF(ds,ds[0],cpedsDirection(0,0), cpedsDirection(res,res),fwhm,fwhm);
	  smo=smoF.smoothGauss();
	  cpeds_matrix_save(smo,"test-cpeds-smooth-corrections.mat");
	  
	  //DEBUG - END
	*/

//	ds.print();

	
	mscsFunction3dregc f;
	f.setSize(1000,1000,1,1,1,1,0,0,0);
	f.allocFunctionSpace();
	cpedsProject pro(ds);
	pro.setInvertLon(false);
	cpedsPointSet3D ps=pro.projectOnPlane(ds[0]*PI180inv,"+proj=%s +lat_0=%lE +lon_0=%lEE +R=%lE");
	ps.save("corrections.ps");
//		ps.clear();
//		ps.append(cpedsPoint3D(0,0,2));
//		ps.append(cpedsPoint3D(-1,-1,3));
//		ps.append(cpedsPoint3D(1,1,4));
//		ps.print();
//		pro.points()=ps;
	
//	ds=pro.projectOnSphere(ds[0]*PI180inv);
	ds.print();
	f.populateField(ps,true,0,0,false);
	f.printInfo();
	long iMax,i,j,k;
		smo=f.getSlice(2,0,0);
		cpeds_matrix_save(smo,"test-cpeds-smooth-corrections.matfnin");
		f.getMaxValue(&iMax,true);
		f.idx2ijk(iMax,i,j,k);
		printf("maxX: %lE, maxY: %lE maxVal: %lE, maxi: %li\n",f.getX(i),f.getY(j), f.fRe(iMax), iMax);
		printf("maxXi: %li, maxYi: %li maxZi: %li\n",i,j,k);
//		exit(0);
	fwhm=pro.getLengthScaleFromAng(0.12*PI180);
//	f.loadtxtLin("cpeds-smooth.mat.in",0);
	printf("fwhm: %lE\n",fwhm);
	f.smooth3DGauss(fwhm,fwhm,1);
		smo=f.getSlice(2,0,0);
		cpeds_matrix_save(smo,"test-cpeds-smooth-corrections.matfn");
		f.getMaxValue(&iMax,true);
		f.idx2ijk(iMax,i,j,k);
		printf("maxX: %lE, maxY: %lE maxVal: %lE, maxi: %li\n",f.getX(i),f.getY(j), f.fRe(iMax), iMax);
		printf("maxXi: %li, maxYi: %li maxZi: %li\n",i,j,k);
	ps.clear();
	ps.append(cpedsPoint3D(f.getX(i),f.getY(j),0.0));
	pro.points()=ps;
	ds=pro.projectOnSphere(ds[0]*PI180inv,"+proj=%s +lat_0=%lE +lon_0=%lEE +R=%lE");
	ds[0].print_direction("maximal signal direction",true,2);
	f.printInfo();
	return 0;
}


void parseOptions(int argc, char** argv) {
	long i;
	_programVersionString=getProgramVersionString();
	try {
		
		//     CmdLine cmd("RT4\n\n. ",' ', "0.1" );
		CmdLine cmd("test-cpeds-smooth\nTorun Centre for Astronomy, UMK, POLAND.\n "
				"This program ...\n\n. "
				"\n"
				"example usage: test-cpeds-smooth"
				"",' ', _programVersionString.c_str() );
		//     CmdLine cmd("RT4 azimutal scan simulator for OCRA-p\n", 'This program simulates the RT4 scan motion, sky rotation during \n\n. ', 0.1 );
		
		// 
		// Define arguments
		//
		
//		ValueArg<double> Ndec("", "Ndec", "number of periods in Lissajous trajectory [21]", false,21,"double");	cmd.add( Ndec );
//		ValueArg<long> nLissajous("", "nLj", "number of timesteps to take for Lissajous trajectory (5000)", false,5000,"long");	cmd.add( nLissajous );
//		
		SwitchArg populateTest("","populateTest", "test method of populating field with a points set (default: true)", false);	cmd.add( populateTest );
//		
//		ValueArg<string> outfile("o","","output file name",false,"scan-trajectory.txt","string"); cmd.add(outfile);
//		ValueArg<string> outdir("","odir","output directory name",false,"./","string"); cmd.add(outdir);
//		
//		std::vector<string> allowedStr;
//		allowedStr.push_back("none");
//		allowedStr.push_back("incAIncH");
//		allowedStr.push_back("incADecH");
//		
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
//		
		_populateTest=populateTest.getValue();
		
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
void populateTest() {
	mscsFunction3dregc f(10,10,1,1,1,1,0,0,0);
	f.allocFunctionSpace();
	
	cpedsPointSet3D ps;
	
	ps.append(cpedsPoint3D(-1,-2,1));
	f.populateField(ps,true,0,0,false,false);
	f.saveSlice(2,0,"popTest.mat",0);
	
	
	
}
