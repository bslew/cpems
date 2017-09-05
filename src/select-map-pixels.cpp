/*!
 * \file select-map-pixels.cpp - selects pixels from loaded healpix map and saves to txt file
 *
 *  Project: Mscs
 *  Created on: Aug 2, 2016 10:07:43 PM
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
#include "Mscs-map.h"
#include "Mscs-function.h"


#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif

string _programVersionString, _mapFile, _outfile, _sel;
double _angPhi,_angTheta;

void parseOptions(int argc, char** argv);
string getProgramVersionString();

/***************************************************************************************/
/*!
	\brief selects directions from map from within circular disks with locations
	defined in file diskFileName all of which have radius diskRadius
	
	\details 
	@param map - nested map to select data from
	@param diskFileName - file specifying directions l,b [deg] of the centres of disks
	@param diskRadius - angular radius of each disk [deg]
	@return set of selected directions with values containing map temperatures

	WARNING: this implementation should work very fast but currently has certain flaw:
	if the disk location happens to be at the boundary of one of the 12 
	healpix base regions then incomplete disk may be selected.
	
	
	In the current implementation if the defined disks are partially overlapping
	then the corresponding pixels will be selected only once,
	but there is also a parallel implementation that selects the pixels twice 
	if this should be a feature then the commented code can be used instead.

	\date Aug 3, 2016, 10:47:17 AM
	\author Bartosz Lew
*/
cpedsDirectionSet selectSmallDisks_file(mscsMap& map,string diskFileName, double diskRadius);

/***************************************************************************************/
/*!
	\brief calculate nside such that the pixel radius is minimally larger than requested radius 
	but the nside cannot be larger than nsideMax
	\details 
	@param radius - requested radius [deg]
	@param nsideMax - maximal nside that defines the minimal pixel size
	@return nside 

	\date Aug 2, 2016, 11:07:35 PM
	\author Bartosz Lew
*/
long get_nside_pixRadius(double radius, long nsideMax);
/***************************************************************************************/
int main(int argc, char **argv) {
	
	cpedsMsgs msgs(".select-map-pixels",false,"",High);
	msgs.setSaveRunWriteMode('a');
	msgs.saveThisRun(argc,argv);
	string s;
	//----------------------------------------------------------------------------------------------------
	//PARSE THE COMMAND LINE ARGUMENTS
	//----------------------------------------------------------------------------------------------------
	parseOptions(argc,argv);

	QString selection=_sel.c_str();
	QStringList selectionL=selection.split(",");
	
	cpedsDirectionSet selectedDirs;
	mscsMap map;
	map.setVerbosityLevel(High);
	
	msgs.say("loading map: "+_mapFile,High);
	map.loadbinT(_mapFile);
	
	msgs.say("selecting directions",High);
	if (selectionL[0]=="disks") {
//		map.set_map_coord();
//		exit(0);
		selectedDirs=selectSmallDisks_file(map,selectionL[1].toStdString(),selectionL[2].toDouble());
	}
	msgs.say("saving selection to file: "+_outfile,High);
	selectedDirs.save(_outfile,"",true);
	return 0;
}


void parseOptions(int argc, char** argv) {
	long i;
	_programVersionString=getProgramVersionString();
	try {
		
		//     CmdLine cmd("RT4\n\n. ",' ', "0.1" );
		CmdLine cmd("select-map-pixels\nTorun Centre for Astronomy, UMK, POLAND.\n "
				"This program ...\n\n. "
				"\n"
				"example usage: select-map-pixels"
				"",' ', _programVersionString.c_str() );
		//     CmdLine cmd("RT4 azimutal scan simulator for OCRA-p\n", 'This program simulates the RT4 scan motion, sky rotation during \n\n. ', 0.1 );
		
		// 
		// Define arguments
		//
		
//		ValueArg<double> Ndec("", "Ndec", "number of periods in Lissajous trajectory [21]", false,21,"double");	cmd.add( Ndec );
//		ValueArg<long> nLissajous("", "nLj", "number of timesteps to take for Lissajous trajectory (5000)", false,5000,"long");	cmd.add( nLissajous );
		
//		SwitchArg horizTrack("","trackHoriz", "follow horizontal sky motion (default: true)", false);	cmd.add( horizTrack );
		
		ValueArg<string> inputfile("m","map","input file name (should be a -Tn-bin file) - healpix nexted 8bit float binary file.",false,"","string"); cmd.add(inputfile);
		ValueArg<string> outfile("o","outfile","output file name (default: selected-pixels.txt)",false,"selected-pixels.txt","string"); cmd.add(outfile);
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
		ValueArg<string> sel("","sel","pixel selection function:\n"
				"- rand_disks,N,size - selects N randomly selected disks of angular radius given by size [deg]\n"
				"- disks,file,size - disks of angular radius given by size using directions (l,b [deg]) specified "
				"in file file",false,"","string"); cmd.add(sel);
		ValueArg<double> angPhi("","angPhi","rotation angle that is added to the list of directions"
				"read from file before selecting map pixels (default: 0 [deg])",false,0.0,"double"); cmd.add(angPhi);
		ValueArg<double> angTheta("","angTheta","rotation angle. The rotation is performed after --angPhi, "
				"before selecting map pixels (default: 0 [deg])",false,0.0,"double"); cmd.add(angTheta);
		//     ValueArg<string> ord("","ord","ordering of the map [ring, nest] default: nest",false,"nest","string"); cmd.add(ord);
		//     SwitchArg mink("M","Mink", "plot also the minkowski functionals for the map", false);	cmd.add( mink );
		//     ValueArg<long> mlevels("", "mlev", "number of levels to calculate the minkowski functionals (default: 100)", false,100,"long");	cmd.add( mlevels );
		
		cmd.parse(argc,argv);
		
		//
		// Set variables
		//
		_outfile = outfile.getValue();
		_mapFile = inputfile.getValue();
		_sel= sel.getValue();
		_angPhi=angPhi.getValue();
		_angTheta=angTheta.getValue();
//		_trackHoriz = horizTrack.getValue();
		
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

/***************************************************************************************/
long get_nside_pixRadius(double radius, long nsideMax) {
	double res0;
	long nsi=1;
	while (nsi<=nsideMax) {
		res0=cpeds_pix_size_healpix(nsi)*PI180inv/2; // coarse pixel radius
		printf("nsi: %li, pix radius [deg]: %lE, radius: %lE\n",nsi, res0,radius);
		if (res0 < radius) return nsi/2;
		nsi*=2;
	}
	return nsi;
}
/***************************************************************************************/
cpedsDirectionSet selectSmallDisks_file(mscsMap& map,string diskFileName, double diskRadius) {
	cpedsDirectionSet selected;
	cpedsDirectionSet dirs0,dummyd;
	cpedsDirection d;
	dirs0.load(diskFileName);
	cpedsList<long> pix0,dummy;
	QList< cpedsList<long> > pix;
	QList< cpedsDirectionSet > pixd;
	mscsFunction selected_pix;
	double th,phi;

	//
	// rotate directions if requested
	//
	if (_angPhi!=0) {
		double angPhi=_angPhi; // deg
		for (long i = 0; i < dirs0.size(); i++) {
			dirs0[i].setLon(cpeds_check_phi( (dirs0[i].lon()+angPhi)*PI180 )*PI180inv);
		}
	}
	if (_angTheta!=0) {
		double ang=_angTheta*PI180; // deg
		for (long i = 0; i < dirs0.size(); i++) {
			dirs0[i]*=PI180;
			dirs0[i].Ry(ang);
			dirs0[i]*=PI180inv;
		}
	}
	
	
	
	//
	// create low resolution map containing the pixel size larger than the requested diskRadius
	//
	
	// calculate map resolution
	long nsideCoarse = get_nside_pixRadius(diskRadius,map.nside());
	printf("nside coarse: %li\n",nsideCoarse);
	
	// decreasing nsideCoarse to improve positioning
	nsideCoarse/=4;

	// making coarse map
	mscsMap cmap;
	cmap.set_nside(nsideCoarse);
	cmap.makekill_space_manager("make", "T");
	cmap.T()=0;
	cmap.set_map_coord();
	
	long idx;
	for (long i = 0; i < dirs0.size(); i++) {
		dirs0[i]*=PI180; // convert directions to radians
//		cmap.set_T(dirs0[i]*PI180,1.0);
		cmap.get_T(dirs0[i],&idx);
		pix0.append(idx);
		cpeds_pix2ang_healpix(cmap.nside(),idx,&th,&phi,1);
#ifdef DEBUG_SELECT_MAP_PIXELS
		printf("th, phi: %lf %lf\n",th*PI180inv,phi*PI180inv);
#endif
		
		pix.append(dummy); // prepare space for storing pixel indexes in high-res map individually for each disk 
		pixd.append(dummyd);
	}
//	// prograde map
//	cmap.change_map_resolution(map.nside());
//	
//	// select directions from the map from within the pixel regions
//	for (long i = 0; i < map.pixNum(); i++) {
//		if (cmap.T(i)==1) pix.append(i);
//	}

	// calculate how many progrades we need
	int prograde=int(log2(map.nside())) - int(log2(nsideCoarse));
	int shift=prograde*2;
	long NpixProg=pow(2.0,shift);
	
#ifdef DEBUG_SELECT_MAP_PIXELS
	printf("shift bits: %i\n",shift);
	printf("NpixProg: %li\n",NpixProg);
#endif
	
	// calculate the 1st pixel ids in the prograded map
	for (long i = 0; i < pix0.size(); i++) {
#ifdef DEBUG_SELECT_MAP_PIXELS
		printf("coarse pixel id: %li\n",pix0[i]);
#endif
		idx=pix0[i];
		pix0[i]=idx << shift;
#ifdef DEBUG_SELECT_MAP_PIXELS
		printf("base high-res pixel id: %li\n",pix0[i]);
#endif
	}
	
/*
	// discard pixels that are outside of the requested diskRadius from the center directions
	double diskRadiusRad=diskRadius*PI180;
	for (long i = 0; i < dirs0.size(); i++) {
		for (long j = pix0[i]; j < pix0[i]+NpixProg; j++) {
			cpeds_pix2ang_healpix(map.nside(),j,&th,&phi,1);
#ifdef DEBUG_SELECT_MAP_PIXELS
			printf("th, phi: %lf %lf\n",th*PI180inv,phi*PI180inv);
#endif
			d.set(phi,PIsnd-th);
			if (d.angle(dirs0[i])<diskRadiusRad) {
				d.setVal(map.get_T(j));
				selected.append(d*PI180inv);
			}
		}
	}
*/

/*
	// discard pixels that are outside of the requested diskRadius from the center directions
	// this parallel implementation generates pixels of overlapping disks selected twice
	double diskRadiusRad=diskRadius*PI180;
	long i;
#pragma omp parallel for default(shared) private (i,th,phi,d) 
	for (i = 0; i < dirs0.size(); i++) {
		for (long j = pix0[i]; j < pix0[i]+NpixProg; j++) {
			cpeds_pix2ang_healpix(map.nside(),j,&th,&phi,1);
			d.set(phi,PIsnd-th);
			if (d.angle(dirs0[i])<diskRadiusRad) {
				d.setVal(map.get_T(j));
				pix[i].append(j);
				pixd[i].append(d);
			}
		}
	}

	// combine results into output list
	for (long i = 0; i < dirs0.size(); i++) {
		for (long j = 0; j < pix[i].size(); j++) {
			selected.append(pixd[i][j]*PI180inv);
		}
	}
*/

	// discard pixels that are outside of the requested diskRadius from the center directions
	double diskRadiusRad=diskRadius*PI180;
	long i;
#pragma omp parallel for default(shared) private (i,th,phi,d) 
	for (i = 0; i < dirs0.size(); i++) {
		for (long j = pix0[i]; j < pix0[i]+NpixProg; j++) {
			cpeds_pix2ang_healpix(map.nside(),j,&th,&phi,1);
			d.set(phi,PIsnd-th);
			if (d.angle(dirs0[i])<diskRadiusRad) {
//				d.setVal(map.get_T(j));
				pix[i].append(j);
//				pixd[i].append(d);
			}
		}
	}

	// combine results into output list
	for (long i = 0; i < dirs0.size(); i++) {
		for (long j = 0; j < pix[i].size(); j++) {
			selected_pix.newPoint(pix[i][j],0.0);
		}
	}
	selected_pix.sortFunctionArgAscending();
	selected_pix.average_sameArgs();

	for (i = 0; i < selected_pix.pointsCount(); i++) {
		long j=selected_pix.getx(i);
		cpeds_pix2ang_healpix(map.nside(),j,&th,&phi,1);
		d.set(phi,PIsnd-th);
		d.setVal(map.get_T(j));
		selected.append(d*PI180inv);
	}

	return selected;
}
