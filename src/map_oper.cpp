/*!
 * \file map_oper.cpp
 *
 *  Project: Mscs
 *  Created on: May 18, 2016 5:23:02 PM
 *  Author: blew
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <tclap/CmdLine.h>
//#include <QtCore/QStringList>
//#include <QtCore/QString>
#include "cpeds-msgs.h"
#include "Mscs-map.h"



#ifndef _NO_NAMESPACE
using namespace std;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif

string _programVersionString;
bool _add;
vector<string> _infiles;
string _outfile;
QString _cmd;
string _hduName,_mask_file;
double _const_value;

void parseOptions(int argc, char** argv);
string getProgramVersionString();
void add_files(vector<string> _infiles);
void subtract_files(vector<string> _infiles);
void multiply_files(vector<string> _infiles);
void multiply_value(vector<string> _infiles, double val);

/*!
	\brief 
	\details 
	@param fwhm [deg]
	@return

	\date May 18, 2016, 6:12:53 PM
	\author Bartosz Lew
*/
void smooth_gaussian(mscsMap& map, double fwhm);
void smoothG(vector<string> _infiles);
void thresold_map(vector<string> _infiles,QString cmd);
void print_map_stat(vector<string> _infiles, string mask_file);
void count_values(vector<string> _infiles, double val);
void convPLANCKfits2Tnbin(vector<string> _infiles);
void convPLANCKfits_1sthdu1col_2Tnbin(vector<string> _infiles,string hdu);
void change_nside(vector<string> _infiles);
void smoothB(vector<string> _infiles);
void convolveSH(mscsMap& map, mscsFunction b);

int main(int argc, char **argv) {
	
	cpedsMsgs msgs(".map_oper",false,"",High);
	msgs.setSaveRunWriteMode('a');
	msgs.saveThisRun(argc,argv);
	string s;
	//----------------------------------------------------------------------------------------------------
	//PARSE THE COMMAND LINE ARGUMENTS
	//----------------------------------------------------------------------------------------------------
	parseOptions(argc,argv);

	if (_add or _cmd=="add")  { add_files(_infiles);	}
	if (_cmd.contains("smoothG"))  { smoothG(_infiles);	}
	if (_cmd.contains("smoothB"))  { smoothB(_infiles);	}
	if (_cmd=="sub")  { subtract_files(_infiles);	}
	if (_cmd=="mul")  { multiply_files(_infiles);	}
	if (_cmd=="mulV")  { multiply_value(_infiles,_const_value);	}
	if (_cmd=="stat")  { print_map_stat(_infiles,_mask_file);	}
	if (_cmd=="convPLANCKfits2Tnbin")  { convPLANCKfits2Tnbin(_infiles);	}
	if (_cmd=="convPLANCKfits_1sthdu1col_2Tnbin")  { convPLANCKfits_1sthdu1col_2Tnbin(_infiles,_hduName);	}
	if (_cmd.contains("renside"))  { change_nside(_infiles);	}
	if (_cmd.contains("thres"))  { thresold_map(_infiles,_cmd);	}
	if (_cmd=="count")  { count_values(_infiles,_const_value);	}

	return 0;
}


void parseOptions(int argc, char** argv) {
	long i;
	_programVersionString=getProgramVersionString();
	try {
		
		//     CmdLine cmd("RT4\n\n. ",' ', "0.1" );
		CmdLine cmd("map_oper\nTorun Centre for Astronomy, UMK, POLAND.\n "
				"This program ...\n\n. "
				"\n"
				"example usage: map_oper"
				"",' ', _programVersionString.c_str() );
		//     CmdLine cmd("RT4 azimutal scan simulator for OCRA-p\n", 'This program simulates the RT4 scan motion, sky rotation during \n\n. ', 0.1 );
		
		// 
		// Define arguments
		//
		
//		ValueArg<double> Ndec("", "Ndec", "number of periods in Lissajous trajectory [21]", false,21,"double");	cmd.add( Ndec );
//		ValueArg<long> nLissajous("", "nLj", "number of timesteps to take for Lissajous trajectory (5000)", false,5000,"long");	cmd.add( nLissajous );
		
		SwitchArg add("","add", "sum input files (default: false)", false);	cmd.add( add );
		
		ValueArg<double> const_value("","val","const value - sub-option for operations on maps",false,0.0,"double"); cmd.add(const_value);
		ValueArg<string> hdu("","hdu","hdu name - sub-option for reading fits files",false,"","string"); cmd.add(hdu);
		ValueArg<string> outfile("o","","output file name",false,"","string"); cmd.add(outfile);
		ValueArg<string> oper("c","cmd","command name (default: '')."
				"Possible values are:\n"
				"add,sub\n"
				"mulV - with value specified with --val option\n"
				"smoothG,fwhp [']\n"
				"smoothB,beam_transfer_function_file\n"
				"thres,ctrVal,low,high - reset maps values to low/high if T<thres/T>=thres respectively \n"
				"stat\n"
				"convPLANCKfits2Tnbin - to extract maps form PLANCK fits files from FULL SKY MAP hdus\n"
				"convPLANCKfits_1sthdu1col_2Tnbin - to extract a single column data from hdu specified with --hdu option\n"
				"",false,"","string"); cmd.add(oper);
		
		ValueArg<string> mask("m","mask","mask file",false,"","string"); cmd.add(mask);
//		std::vector<string> allowedStr;
//		allowedStr.push_back("none");
//		allowedStr.push_back("incAIncH");
//		allowedStr.push_back("incADecH");
		
//		ValuesConstraint<string>* allowedStrNew;
//		allowedStrNew = new ValuesConstraint<string>(allowedStr);
//		ValueArg<string> traj("", "traj", "trajectory type:", false,"none", allowedStrNew);	cmd.add( traj );
//		allowedStr.clear();

		UnlabeledMultiArg<string> infiles("infiles", "input file names ","", "string");     	cmd.add( infiles );
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
		_infiles = infiles.getValue();  if (_infiles.size() == 0) { printf("ERROR: not enough input files specified\n"); exit(0);}
//		_Ndec = Ndec.getValue();
//		_nLissajous=nLissajous.getValue();
		_add = add.getValue();
		string operLoc=oper.getValue(); _cmd=operLoc.c_str();
		_outfile = outfile.getValue();
//		_outdir = outdir.getValue(); 	if (_outdir=="") _outdir=".";
		_hduName=hdu.getValue();
		_const_value=const_value.getValue();
		_mask_file = mask.getValue(); 
		
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
/***************************************************************************************/
void add_files(vector<string> _infiles) {
	mscsMap m1,m2;
	printf("loading %s\n",_infiles[0].c_str());
	m1.loadbinT(_infiles[0]);
	for (unsigned long i = 1; i < _infiles.size(); i++) {
		printf("loading %s\n",_infiles[i].c_str());
		m2.loadbinT(_infiles[i]);
		printf("adding\n");
		m1.T()+=m2.T();
	}
	printf("saving to file: %s\n",_outfile.c_str());
	m1.savebinT(_outfile);
}
/***************************************************************************************/
void subtract_files(vector<string> _infiles) {
	mscsMap m1,m2;
	printf("loading %s\n",_infiles[0].c_str());
	m1.loadbinT(_infiles[0]);
	for (unsigned long i = 1; i < _infiles.size(); i++) {
		printf("loading %s\n",_infiles[i].c_str());
		m2.loadbinT(_infiles[i]);
		printf("subtracting\n");
		m1.T()-=m2.T();
	}
	printf("saving to file: %s\n",_outfile.c_str());
	m1.savebinT(_outfile);
}
/***************************************************************************************/
void multiply_files(vector<string> _infiles) {
	mscsMap m1,m2;
	printf("loading %s\n",_infiles[0].c_str());
	m1.loadbinT(_infiles[0]);
	for (unsigned long i = 1; i < _infiles.size(); i++) {
		printf("loading %s\n",_infiles[i].c_str());
		m2.loadbinT(_infiles[i]);
		printf("multiplying\n");
		m1.T()*=m2.T();
	}
	printf("saving to file: %s\n",_outfile.c_str());
	m1.savebinT(_outfile);
}
/***************************************************************************************/
void multiply_value(vector<string> _infiles, double val) {
	mscsMap m1;
	for (unsigned long i = 0; i < _infiles.size(); i++) {
		printf("loading %s\n",_infiles[i].c_str());
		m1.loadbinT(_infiles[i]);
		printf("multiplying\n");
		m1.T()*=val;
	}
	printf("saving to file: %s\n",_outfile.c_str());
	m1.savebinT(_outfile);
}
/***************************************************************************************/
void smoothG(vector<string> _infiles) {
	mscsMap m1;
	QStringList qsl=_cmd.split(",");
	double fwhm=qsl[1].toDouble()/60; // convert from ' to deg
	
	for (unsigned long i = 0; i < _infiles.size(); i++) {
		printf("loading %s\n",_infiles[i].c_str());
		m1.loadbinT(_infiles[i]);
		printf("gaussian smoothing (fwhm=%lf [deg])\n",fwhm);
		
		smooth_gaussian(m1,fwhm);
		printf("saving to file: %s\n",_outfile.c_str());
		m1.savebinT(_outfile);
		return;
	}
}
/***************************************************************************************/
void smoothB(vector<string> _infiles) {
	mscsMap m1;
	QStringList qsl=_cmd.split(",");
	QString beamFile=qsl[1]; 
	mscsFunction b;
	b.load(beamFile.toStdString());
	
	for (unsigned long i = 0; i < _infiles.size(); i++) {
		printf("loading %s\n",_infiles[i].c_str());
		m1.loadbinT(_infiles[i]);
		printf("convolving with transfer function  %s\n",beamFile.toStdString().c_str());
		
		convolveSH(m1,b);
		printf("saving to file: %s\n",_outfile.c_str());
		m1.savebinT(_outfile);
		return;
	}
}
/***************************************************************************************/
void convolveSH(mscsMap& map, mscsFunction b) {
	mscsWindowFunction bl;
	long lmax=map.nside()*2;
	for (long l = 0; l < lmax; l++) {
		bl.newPoint(l,b.getY(l));
	}

	mscsWindowFunction wf;
	wf.make_unit_kernel(lmax);
//	map.setVerbosityLevel(High);
	printf("analysis\n");
	mscsAlms alm=map.SH_analysis(lmax,wf,wf);
	alm.savetxtAlms("alms","RI");
	printf("synthesis\n");
	map.SH_synthesis(alm,lmax,bl,wf);
	
}
/***************************************************************************************/
void smooth_gaussian(mscsMap& map, double fwhm) {
	mscsWindowFunction bl;
	long lmax=map.nside()*2;
	bl.make_unit_kernel(lmax);

	mscsWindowFunction wf=bl;
	bl.make_gaussian_kernel_kspace(fwhm*PI180);
//	map.setVerbosityLevel(High);
	printf("analysis\n");
	mscsAlms alm=map.SH_analysis(lmax,wf,wf);
	alm.savetxtAlms("alms","RI");
	printf("synthesis\n");
	map.SH_synthesis(alm,lmax,bl,wf);
}
/***************************************************************************************/
void print_map_stat(vector<string> _infiles, string mask_file) {
	mscsMap m1;
	m1.setVerbosityLevel(High);
	for (unsigned long i = 0; i < _infiles.size(); i++) {
		printf("loading %s\n",_infiles[i].c_str());
		m1.loadbinT(_infiles[i]);
		if (mask_file!="") m1.loadbinm(mask_file);
		m1.calculate_map_stats(1);
	}
}
/* ******************************************************************************************** */
void count_values(vector<string> _infiles, double val) {
	mscsMap m1;
	long n=0;
	m1.setVerbosityLevel(High);
	for (unsigned long i = 0; i < _infiles.size(); i++) {
		printf("loading %s\n",_infiles[i].c_str());
		m1.loadbinT(_infiles[i]);
		m1.calculate_map_stats(1);
		for (long i = 0; i < m1.pixNum(); i++) {
			if (m1.T(i)==val) n++;
		}
		printf("number of %lE values in the map: %li\n",val,n);
	}
}
/***************************************************************************************/
void convPLANCKfits2Tnbin(vector<string> _infiles) {
	mscsMap m1;
	m1.setVerbosityLevel(High);
	for (unsigned long i = 0; i < _infiles.size(); i++) {
		printf("loading %s\n",_infiles[i].c_str());
		m1.loadPLANCK_temp_fits(_infiles[i],"\'FULL SKY MAP\'");
		m1.savebinT(_outfile);
	}
}
/***************************************************************************************/
void convPLANCKfits_1sthdu1col_2Tnbin(vector<string> _infiles,string hduName) {
	mscsMap m1;
	m1.setVerbosityLevel(High);
	for (unsigned long i = 0; i < _infiles.size(); i++) {
		printf("loading %s\n",_infiles[i].c_str());
		m1.loadPLANCK_temp_fits(_infiles[i],hduName);
		m1.savebinT(_outfile);
	}	
}
/***************************************************************************************/
void change_nside(vector<string> _infiles) {
	mscsMap m1;
	QStringList qsl=_cmd.split(",");
	long nside=qsl[1].toLong(); 
	
	for (unsigned long i = 0; i < _infiles.size(); i++) {
		printf("loading %s\n",_infiles[i].c_str());
		m1.loadbinT(_infiles[i]);
		printf("changing resolution to nside=%li \n",nside);
		
		m1.change_map_resolution(nside);
		printf("saving to file: %s\n",_outfile.c_str());
		m1.savebinT(_outfile);
		return;
	}
}
/* ******************************************************************************************** */
void thresold_map(vector<string> _infiles,QString cmd) {
	mscsMap m1;
	QStringList qsl=_cmd.split(",");
	double ctr=qsl[1].toDouble(); 
	double l=qsl[2].toDouble(); 
	double h=qsl[3].toDouble(); 

	for (unsigned long i = 0; i < _infiles.size(); i++) {
		printf("loading %s\n",_infiles[i].c_str());
		m1.loadbinT(_infiles[i]);
		printf("thresholding map\n");
		
		for (long j = 0; j < m1.pixNum(); j++) {
			if (m1.T()[j]>=ctr) m1.T()[j]=h;
			else m1.T()[j]=l;
		}
		
		printf("saving to file: %s\n",_outfile.c_str());
		m1.savebinT(_outfile);
		return;
	}

}
