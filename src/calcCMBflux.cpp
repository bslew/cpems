/*!
 * \file calcCMBflux.cpp
 *
 *  Project: Mscs
 *  Created on: Apr 18, 2013 8:54:20 PM
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
#include "Mscs-function3dregc.h"
#include "Mscs-global-defs.h"


#ifndef _NO_NAMESPACE
using namespace std;
//using namespace math;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif

string _programVersionString;
string _mapFile, _hdf5dset,_fpblist, _outFile;
string _includeCols; // _fluxPerBeamRelList sub-option
//string _paralacticAngleRange;
double _paralacticAngleRange;
vector<double> _phi; // paralactic angle range
double _pixSize, _beamSize, _freq, _aperture, _thSt,_thEn,_dth,_pointingAcc, _refBeamSep;
long _i,_j,_NrefBeam, _numTh;
bool _fluxPerBeam,_fluxPerAperture,_calculateCMBflux_vs_angSize,_exactl, _fluxPerBeamList,_fluxPerBeamRelList, _mJy;
bool _paralacticAngleRange_oneSided;
double _geo;
double _haloMassThres;

void parseOptions(int argc, char** argv);
string getProgramVersionString();


/*!
	\brief average relative flux densities
	\details 
	@param S - array with relative flux densities for all halos
	@param colSt - start averaging S matrix from column colSt
	@return 3xNhalo column array containing halo id, mean flux density and flux density standard deviation 
	calculated from S
	
	S array contains Ny rows, each for every halos
	Each row contains Nx columns where the first column has stored haloID and 
	the remaining columns relative flux density measurements from different paralactic angles.

	\date Apr 20, 2016, 7:16:47 PM
	\author Bartosz Lew
*/
mscsFunction3dregc averageRelativeFluxDensities(mscsFunction3dregc& S, long colSt);

void calculateFluxPerBeam(mscsFunction3dregc& f);
void calculateFluxPerAperture(mscsFunction3dregc& f);
/*!
	\brief calculate flux densities from map f at a number of locations
	\details 
	@param f - field with specific intensity
	@param xyList - an object containing a list of pixel coordinates (x,y) at which a gaussian beam should be placed
	@param coli - a column in xyList object containing the pixel indexes along x axis in the field f
	@param colj - a column in xyList object containing the pixel indexes along y axis in the field f	
	@return

	\date Apr 19, 2016, 10:57:44 AM
	\author Bartosz Lew
*/
mscsFunction3dregc calculateFluxPerBeamList(mscsFunction3dregc& f,mscsFunction3dregc& xyList, long coli, long colj);
/*!
	\brief calculate flux densities from map f at a number of locations relative to the flux density at some other location
	\details 
	@param f - field with specific intensity
	@param xyList - an object containing a list of pixel coordinates (x,y) at which a gaussian beam should be placed
	The columns in this list should be compatible with the output txt file from the lss_pyramidLOSintegrate program.
	In particular column 32 is used for screen halos agains virial mass.
	@param coli - a column in xyList object containing the pixel indexes along x axis in the field f
	@param colj - a column in xyList object containing the pixel indexes along y axis in the field f	
	@param refBeamDist - distance of the reference beam. It's parameters are defined in the same way as the primary beam
	(i.e. via the global variables _aperture and _freq)
	@param includeCols - list of columns to be copied from xyList to the output array as the first columns 
	@return

	\date Apr 19, 2016, 10:57:44 AM
	\author Bartosz Lew
*/
mscsFunction3dregc calculateFluxPerBeamRelList(mscsFunction3dregc& f,mscsFunction3dregc xyList, long coli, long colj, cpedsList<long> includeCols);
mscsFunction calculateCMBflux_vs_angSize(mscsFunction3dregc& f);

int main(int argc, char **argv) {
	
	cpedsMsgs msgs(".calcCMBflux",false,"",High);
	msgs.setSaveRunWriteMode('a');
	msgs.saveThisRun(argc,argv);
	string s;
	//----------------------------------------------------------------------------------------------------
	//PARSE THE COMMAND LINE ARGUMENTS
	//----------------------------------------------------------------------------------------------------
	parseOptions(argc,argv);
	
	mscsFunction3dregc f;

	if (Mscs_getExtension(_mapFile)=="hdf" or Mscs_getExtension(_mapFile)=="hdf5") {
		cpedsStatusCodes res=f.loadHDF5(_mapFile,_hdf5dset);
		if (res!=cpedsSuccess) {
			QString hdf5dset=_hdf5dset.c_str();
			if (hdf5dset.contains("DLcmb-")) {
				msgs.say("Procedding special case for DLcmb- dataset",High);
				// we typically don't store such dataset, so try to construct it from other datasets
				mscsFunction3dregc f1,f2;
				QString hdf5dset1=hdf5dset;
				hdf5dset1.replace("DLcmb-","DLcmbTSZ-");
				QString hdf5dset2=hdf5dset;
				hdf5dset2.replace("DLcmb-","DLTSZ-");

				msgs.say("Loading dataset: "+hdf5dset1.toStdString(),High);
				cpedsStatusCodes res1=f1.loadHDF5(_mapFile,hdf5dset1.toStdString());
				msgs.say("Loading dataset: "+hdf5dset2.toStdString(),High);
				cpedsStatusCodes res2=f2.loadHDF5(_mapFile,hdf5dset2.toStdString());
				
				printf("max value1: %lE\n",f1.getMaxValue());
				printf("max value2: %lE\n",f2.getMaxValue());
				if (res1==cpedsSuccess and res2==cpedsSuccess) {
					msgs.say("Calculating difference",High);
					f=f1-f2;
					printf("max value: %lE\n",f.getMaxValue());
				}
			}
		}
	}
	else 
		f.loadMatrix(_mapFile);

	if (_fluxPerBeamList) {
		mscsFunction3dregc xyList,S;
		QStringList qsl;
		QString s=_fpblist.c_str();
		qsl=s.split(",");
		msgs.say("assumed list file name: "+qsl[0].toStdString(),Low);
		msgs.say("assumed X column: "+qsl[1].toStdString(),Low);
		msgs.say("assumed Y column: "+qsl[2].toStdString(),Low);
		if (qsl.size()!=3) {
			msgs.criticalError("wrong syntax in --fpblist",Top);
		}
		long listX=qsl[1].toDouble();
		long listY=qsl[2].toDouble();
		xyList.loadMatrix(qsl[0].toStdString());

		msgs.say("Calculating flux per beam around a list of locations",Top);
		// reset the coordinates
		f.setSize(f.Nx(),f.Ny(),f.Nz(),1,1,1,0,0,0);
		S=calculateFluxPerBeamList(f,xyList,listX,listY);

		S.saveSlice(2,0,_outFile,0);
		
	}

	if (_fluxPerBeamRelList) {
		mscsFunction3dregc xyList,S;
		QStringList qsl;
		QString s=_fpblist.c_str();
		qsl=s.split(",");
		msgs.say("assumed list file name: "+qsl[0].toStdString(),Low);
		msgs.say("assumed X column: "+qsl[1].toStdString(),Low);
		msgs.say("assumed Y column: "+qsl[2].toStdString(),Low);
		if (qsl.size()!=3) {
			msgs.criticalError("wrong syntax in --fpbrlist",Top);
		}
		long listX=qsl[1].toLong();
		long listY=qsl[2].toLong();
		xyList.loadMatrix(qsl[0].toStdString());

		cpedsList<long> includeCols;
		s=_includeCols.c_str();
		qsl=s.split(",");
		for (long i = 0; i < qsl.size(); i++) {			includeCols.append(qsl[i].toLong());		}

		msgs.say("Calculating flux per beam around a list of locations",Top);
		// reset the coordinates
		f.setSize(f.Nx(),f.Ny(),f.Nz(),1,1,1,-0.5,-0.5,-0.5);
		S=calculateFluxPerBeamRelList(f,xyList,listX,listY,includeCols);

//		S.saveSlice(2,0,_outFile+".full",0); // debug purpose only
		
		// calculate mean difference flux density and its standard deviation
		S=averageRelativeFluxDensities(S,includeCols.size()+1);
		S.saveSlice(2,0,_outFile,0);
		
	}

//	if (_fluxPerBeamRelList_cmdLineParams) {
//		mscsFunction3dregc xyList,S;
//		QStringList qsl;
//		QString s=_fpblist.c_str();
//		qsl=s.split(",");
//		msgs.say("assumed list file name: "+qsl[0].toStdString(),Low);
//		msgs.say("assumed X column: "+qsl[1].toStdString(),Low);
//		msgs.say("assumed Y column: "+qsl[2].toStdString(),Low);
//		if (qsl.size()!=3) {
//			msgs.criticalError("wrong syntax in --fpbrlist",Top);
//		}
//		long listX=qsl[1].toLong();
//		long listY=qsl[2].toLong();
//		xyList.loadMatrix(qsl[0].toStdString());
//
//		cpedsList<long> includeCols;
//		s=_includeCols.c_str();
//		qsl=s.split(",");
//		for (long i = 0; i < qsl.size(); i++) {			includeCols.append(qsl[i].toLong());		}
//
//		msgs.say("Calculating flux per beam around a list of locations",Top);
//		// reset the coordinates
//		f.setSize(f.Nx(),f.Ny(),f.Nz(),1,1,1,-0.5,-0.5,-0.5);
//		S=calculateFluxPerBeamRelList(f,xyList,listX,listY,includeCols);
//
////		S.saveSlice(2,0,_outFile+".full",0); // debug purpose only
//		
//		// calculate mean difference flux density and its standard deviation
//		S=averageRelativeFluxDensities(S,includeCols.size()+1);
//		S.saveSlice(2,0,_outFile,0);
//		
//	}
	
	
	if (_fluxPerBeam) {
		msgs.say("Calculating flux per beam: ",Top);
		calculateFluxPerBeam(f);
	}

	if (_fluxPerAperture) {
		msgs.say("Calculating flux per aperture: ",Top);
		calculateFluxPerAperture(f);
	}

	if (_calculateCMBflux_vs_angSize) {
		msgs.say("Calculating flux per beam from CMB power spectrum: ",Top);
		mscsFunction Sth=calculateCMBflux_vs_angSize(f);
		string fname="Sth-f_"+msgs.toStrf(_freq,0)+"GHz-D_"+msgs.toStrf(_aperture,0)+"m";
		if (_exactl)
			fname+="-exactl";
		Sth.save(fname);
	}
	
	
	return 0;
}


void parseOptions(int argc, char** argv) {
	long i;
	_programVersionString=getProgramVersionString();
	try {
		
		//     CmdLine cmd("RT4\n\n. ",' ', "0.1" );
		CmdLine cmd("calcCMBflux\nTorun Centre for Astronomy, UMK, POLAND.\n "
				"This program ...\n\n. "
				"\n"
				"example usage: calcCMBflux"
				"",' ', _programVersionString.c_str() );
		//     CmdLine cmd("RT4 azimutal scan simulator for OCRA-p\n", 'This program simulates the RT4 scan motion, sky rotation during \n\n. ', 0.1 );
		
		// 
		// Define arguments
		//
		
		ValueArg<double> pixSize("", "pixSize", "angular pixel size in the fpb file [deg]", false,1,"double");	cmd.add( pixSize );
		//		ValueArg<double> beamSize("", "beamSize", "angular gaussian beam fwhm [deg]", true,1,"double");	cmd.add( beamSize );
		ValueArg<double> freq("", "freq", "beam frequency [GHz] (default: 30)", false,30,"double");	cmd.add( freq );
		ValueArg<double> aperture("D", "aperture", "telescope aperture. Used to derive theoretical gaussian beam size [m] (default: 32)", false,32,"double");	cmd.add( aperture );
		
		ValueArg<double> thSt("", "thSt", "minimal angular size to probe the CMB flux per beam vs theta relation [deg] (default: 0.001 )", false,0.001,"double");	cmd.add( thSt );
		ValueArg<double> thEn("", "thEn", "minimal angular size to probe the CMB flux per beam vs theta relation [deg] (default: 1 )", false,1,"double");	cmd.add( thEn );
		ValueArg<double> dth("",  "dth", "minimal angular size to probe the CMB flux per beam vs theta relation [deg] (default: 0.001 )", false,0.001,"double");	cmd.add( dth );
		SwitchArg exactl("","exactl", "use for option --fCl. Calculate flux for temperature rms fluctuation calculated exactly at the corresponding multipole. (default: false)", false);	cmd.add( exactl );
		ValueArg<double> geo("",  "geo", "use for option --fCl. Increase dth in geomethi sequence using this multiplied (default: 1.0 )", false,1.0,"double");	cmd.add( geo );
		
		ValueArg<long> beam_i("", "bi", "field matrix index along x of the center of the beam ", false,0,"long");	cmd.add( beam_i );
		ValueArg<long> beam_j("", "bj", "field matrix index along y of the center of the beam ", false,0,"long");	cmd.add( beam_j );
		
		SwitchArg mJy("","mJy", "output flux-densities in mJy (default: false)", false);	cmd.add( mJy );
		SwitchArg fluxPerBeam("","fpb", "calculate field flux density per beam (default: false)", false);	cmd.add( fluxPerBeam );
		SwitchArg fluxPerAperture("","fpa", "integrate field within aperture --aperture given in deg at --bi --bj pixel in the field."
				"The pixel size needs to be also specified in deg. The output will be in arcmin^2. (default: false)", false);	cmd.add( fluxPerAperture );
		SwitchArg calcCMBflux_vs_angSize("","fCl", "calculate CMB flux per beam vs angular scale relation using provided input CMB power spectra. (default: false). "
				"In this case the -f option must be the input power spectrum in units of [K] not multiplied by l(l+1)/2pi."
				"The -D and -freq options are needed as CMB flux is frequnecy dependent and beam size dependent (i.e. apterture dependent)."
				"Use --thSt --thEn --dth for this option to specify the interesting range to probe the output function.", false);	cmd.add( calcCMBflux_vs_angSize );

		SwitchArg fluxPerBeamList("","fpbl", "calculate field flux density per beam using a list of i,j pixel positions in the provided"
				"map. The list should be provided using --list option (default: false)", false);	cmd.add( fluxPerBeamList );
		SwitchArg fluxPerBeamRelList("","fpbrl", "calculate field flux density per beam relative to a reference beam for "
				"a list of i,j pixel positions in the provided map. "
				"The list should be provided using --list option (default: false). "
				"For this option you should also consider specifying the pointing accuracy with --pointingAcc option and "
				"relative beam angular separation with --relBeamSep option and "
				"the number of reference beam pointings using option --NrelBeam", false);	cmd.add( fluxPerBeamRelList );
		ValueArg<long> NrefBeam("", "NrelBeam", "fpbrl sub-option: Number of reference beam pointings. These are chosen randomly"
				"around the main beam location. (default: 500)", false,500,"long");	cmd.add( NrefBeam );
		ValueArg<double> pointingAcc("", "pointingAcc", "fpbrl sub-option: telescope pointing accuracy [fraction of fwhm defined by --freq and -D]. "
				"If non-zero then the provided list of pointings will be offset in a random direction by "
				"pointingAcc. This value could be associated with the pointing/tracking accuracy "
				"over the time scales relevant for the ON-OFF SZ observations. (default: 0.000)", false,0.0,"double");	cmd.add( pointingAcc );
		ValueArg<double> refBeamSep("", "refBeamSep", "fpbrl sub-option: reference beam angular separation from the main beam [deg]."
				"(default: 0.0526)", false,0.0526,"double");	cmd.add( refBeamSep);
		ValueArg<string> mapFile("f","field","input file name containing matrix with the field of spectral radiance in units of [Jy/sr]",false,"radiance.mat","string"); cmd.add(mapFile);
		ValueArg<string> hdf5dset("","hdf5dset","name of the hdf5 dset to use (default: '')",false,"","string"); cmd.add(hdf5dset);
		ValueArg<string> fpblist("","fpblist","file name containing the x,y locations of the pixels around which the flux should be calculated with --fpbl option."
				"After the file name the two coma-separated integers should follow indicating the columns inside of the file that should represent the"
				"x coordinate and y coordinate. (default: ''). "
				"For example:\n "
				"--fpblist listFile.txt,10,13\n"
				"which indicates that column 10 should be used for x coordinate and 13'th for y coordinate counting from 0."
				"",false,"","string"); cmd.add(fpblist);
		ValueArg<double> haloMassThres("",  "haloMassThres", "sub-option for fpblist. Defines minimal halo mass required for "
				"flux density calculation [10^10 Msol/h] (default: 1e4 )", false,1.0e4,"double");	cmd.add( haloMassThres );
		ValueArg<long> numTh("",  "numTh", "Number of OMP threads to use for the run."
				"(default: 0 -- auto )", false,0,"long");	cmd.add( numTh);
		
		ValueArg<string> outFile("o","","output file name (default: out)",false,"out","string"); cmd.add(outFile);
		ValueArg<string> includeCols("","includeCols","--fpbrl sub-option. Coma-separated list of columns to be"
				"included in the output file and copied from the file specified with --fpblist option. (default: 0,9,32,34)."
				"The default choice demands that the 4 first columns of the output file should be "
				"filled with the data from columns 0,9,32 and 34 respectively from the input file given with --fpblist option."
				"Since typically this program takes the outputs from the lss_pyramidLOSintegrate program as input,"
				"this choice corresponds to haloID,redshift,virial mass and angular size of virial radius.",false,"out","string"); cmd.add(includeCols);
//		ValueArg<string> paralacticAngleRange("","paralacticAngleRange","--fpbrl sub-option. Coma-separated ranges of paralactic angles [deg]"
//				"to be considered when calculating difference flux density. (default: 0,180,180,360)."
//				"Exactly two ranges should be specified. The default choice corresponds to two complementary pieces of the full angle"
//				"",false,"0,180,180,360","string"); cmd.add(paralacticAngleRange);
		ValueArg<double> paralacticAngleRange("",  "paralacticAngleRange", "sub-option for fpblist. "
				"Defines range of paralactic angles at which the reference flux density calculation can be done. "
				"E.g. The value of 90 means that allowed PAs are from within [0,90] sum [180,270]."
				"(default: 180 - i.e. the full ring around the source position is available for reference flux density calculation )", false,180.0,"double");	cmd.add( paralacticAngleRange );
		SwitchArg paralacticAngleRange_oneSided("","paralacticAngleRange_oneSided", "--paralacticAngleRange sub-option."
				"It forces that the paralactic angle range is not symmetrical (default) but defined only within [0,paralacticAngleRange] range."
				"(default: false)", false);	cmd.add( paralacticAngleRange_oneSided );
		
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
		
		
		cmd.parse(argc,argv);
		
		//
		// Set variables
		//
		_fluxPerAperture = fluxPerAperture.getValue();
		_fluxPerBeam = fluxPerBeam.getValue();
		_calculateCMBflux_vs_angSize = calcCMBflux_vs_angSize.getValue();
		_thSt=thSt.getValue();
		_thEn=thEn.getValue();
		_dth=dth.getValue();
		_pixSize=pixSize.getValue();
		//		_beamSize=beamSize.getValue();
		_freq=freq.getValue();
		_aperture=aperture.getValue();
		_mapFile= mapFile.getValue();
		_i=beam_i.getValue();
		_j=beam_j.getValue();
		_exactl=exactl.getValue();
		_geo=geo.getValue();
		_haloMassThres=haloMassThres.getValue();
		_hdf5dset=hdf5dset.getValue();
		_fluxPerBeamList=fluxPerBeamList.getValue();
		_fluxPerBeamRelList=fluxPerBeamRelList.getValue();
		_refBeamSep=refBeamSep.getValue();
		_pointingAcc=pointingAcc.getValue();
		_NrefBeam=NrefBeam.getValue();
		_fpblist=fpblist.getValue();
		_outFile=outFile.getValue();
		_includeCols=includeCols.getValue();
		_mJy=mJy.getValue();
		_numTh=numTh.getValue();
		_paralacticAngleRange=paralacticAngleRange.getValue();
		_paralacticAngleRange_oneSided=paralacticAngleRange_oneSided.getValue();
		
		_phi.push_back(0.0); 
		_phi.push_back(_paralacticAngleRange*PI180 ); 
		if (_paralacticAngleRange_oneSided==false) {
			_phi.push_back( (180.0+0.0)*PI180 ); 
			_phi.push_back( ( 180.0+_paralacticAngleRange)*PI180 ); 
		}
		
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
void calculateFluxPerBeam(mscsFunction3dregc& f) {
	double fwhm=cpeds_calculateAngularResolution(_aperture,_freq) ; //deg
	fwhm/=_pixSize;
	//	double fwhm=_beamSize/(_pixSize);
	double s=fwhm/(2.0*sqrt(2.0*log(2.0)));
	printf("fwhm [pix]: %lf\n",fwhm);
	//	msgs.say("lambda [cm]: %lf\n",lambda*100, Medium);
	
	
	mscsFunction3dregc beam=f;	
	beam.mkGauss3D(s,s,s);
	beam.shift(_i-f.Nx()/2,_j-f.Ny()/2,0);
	beam/=beam.getMaxValue();
	beam.saveSlice(2,0,"beamPosition.mat",0);
	
	f*=beam;
	f.saveSlice(2,0,"CMB_times_beam.mat",0);
	double Wpix=_pixSize*_pixSize*pow(PI180,2);
	double F=f.sumRe()*Wpix;
	printf("pixSize [sr]: %lE\n",Wpix);
	printf("flux per beam [Jy]: %lE\n",F);
	
}
/***************************************************************************************/
void calculateFluxPerAperture(mscsFunction3dregc& f) {
	double fwhm=_aperture; //deg
//	fwhm/=_pixSize;
	printf("fwhm [deg]: %lf\n",fwhm);
	double xmax,ymax;
	xmax=f.Nx()*_pixSize/2;
	ymax=f.Ny()*_pixSize/2;
	f.setSizeRange(f.Nx(),f.Ny(),f.Nz(),-xmax,-ymax,0,xmax,ymax,0);
	
	mscsFunction3dregc beam=f;	
	beam=0.0;
	beam.mkBall3D(fwhm/2);
	beam.shift(_i-f.Nx()/2,_j-f.Ny()/2,0);
	beam/=beam.getMaxValue();
	beam.saveSlice(2,0,"beamPosition.mat",0);
	
	f*=beam;
	f.saveSlice(2,0,"field_times_beam.mat",0);
	double Wpix=_pixSize*_pixSize*3600.0; // convert to arcmin^2 from deg^2
	double F=f.sumRe()*Wpix;
//	printf("pixSize [sr]: %lE\n",Wpix);
	printf("integral inside the aperture [fieldUnit * arcmin^2]: %lE\n",F);
	
}
/***************************************************************************************/
mscsFunction3dregc calculateFluxPerBeamList(mscsFunction3dregc& f,mscsFunction3dregc& xyList, long coli, long colj) {
	double fwhm=cpeds_calculateAngularResolution(_aperture,_freq) ; //deg
	fwhm/=_pixSize; // convert to pixels
	//	double fwhm=_beamSize/(_pixSize);
	double s=fwhm/(2.0*sqrt(2.0*log(2.0)));
	printf("fwhm [pix]: %lf\n",fwhm);
	//	msgs.say("lambda [cm]: %lf\n",lambda*100, Medium);
	long cutRadius=5*fwhm; // in pixels
	long cutSize=cutRadius*2+1;
	long di,dj;
	mscsFunction3dregc beam,source;
	mscsFunction3dregc S;
	S.setSize(1,xyList.Ny());
	S.allocFunctionSpace();
	
	double FluxDensityUnitConv=1; // unit conversion factor: 1 for Jy
	if (_mJy) FluxDensityUnitConv=1000;
	
	double Wpix=_pixSize*_pixSize*pow(PI180,2)*FluxDensityUnitConv;
	printf("pixSize [sr]: %lE\n",Wpix);
	long i,j;
	for (unsigned long row = 0; row < xyList.Ny(); row++) {
		i=xyList(coli,row);
		j=xyList(colj,row);
#ifdef DEBUG_CALCCMBFLUX
		printf("i: %li, j:%li,  cutx: %li, %li, cuty: %li, %li\n",i,j,i-cutRadius,i+cutRadius,j-cutRadius,j+cutRadius);
#endif
		source=f.cutAwayBlock(i-cutRadius,j-cutRadius,0,i+cutRadius,j+cutRadius,0);	
#ifdef DEBUG_CALCCMBFLUX
		source.saveSlice(2,0,"source.mat",0);
#endif
		beam=source; // setting size
		beam.mkGauss2D(i,j,s,s,1,2,0);
#ifdef DEBUG_CALCCMBFLUX
		beam.setVerbosityLevel(Top);
		beam.printInfo();
		beam.saveSlice(2,0,"beamPosition.mat",0);
#endif
		source*=beam;

#ifdef DEBUG_CALCCMBFLUX
		source.saveSlice(2,0,"source_times_beam.mat",0);
#endif
		S(0,row)=source.sumRe()*Wpix;
		printf("%li) flux per beam [Jy]: %lE\n",row,S(0,row));
#ifdef DEBUG_CALCCMBFLUX
		exit(0);
#endif
	}

#ifdef DEBUG_CALCCMBFLUX
	S.saveSlice(2,0,"fpbl.flux",0);
#endif
	
	return S;
}
/***************************************************************************************/
mscsFunction3dregc calculateFluxPerBeamRelList(mscsFunction3dregc& f,mscsFunction3dregc xyList, long coli, long colj, cpedsList<long> includeCols) {
	double fwhm=cpeds_calculateAngularResolution(_aperture,_freq) ; //deg
	fwhm/=_pixSize; // convert to pixels
	//	double fwhm=_beamSize/(_pixSize);
	double s=fwhm/(2.0*sqrt(2.0*log(2.0)));
	//	msgs.say("lambda [cm]: %lf\n",lambda*100, Medium);
	double Wpix=_pixSize*_pixSize*pow(PI180,2);
	printf("pixSize [sr]: %lE\n",Wpix);
	printf("fwhm [pix]: %lf\n",fwhm);
	
	
	double FluxDensityUnitConv=1; // unit conversion factor: 1 for Jy
	if (_mJy) FluxDensityUnitConv=1000;
	Wpix*=FluxDensityUnitConv; // put flux density unit conversion here to avoid extra multiplication
	
	long Nref=_NrefBeam;
	double cutRadius=2.5*fwhm; // in pixels
	if (cutRadius<=1) { printf("curRadius [pix]: %lf\n. This probably doesn't make sense. exiting",cutRadius); exit(-1); }
	double MaxRefBeamSep=0.2; // deg
	long refBeamRegionRadius=cutRadius+MaxRefBeamSep/_pixSize + 2*fwhm; /* [pix]
	 * Comment: we add + 2*fwhm to make sure we have the same halos processed for different beam pointing accuracy
	 * and through all other possible variations.
	 * 
	 * author: blew
	 * date: May 31, 2016 6:18:32 PM
	 *
	 */
	
//	long refBeamSepPix=_refBeamSep/_pixSize; // in pixels
	cpedsRNG rns("uniform");
	rns.setMinMax(0.0,twoPI); // we will randomly choose offset direction
	cpedsList<double> phiAll,phi,dX,dY;

	mscsFunction3dregc S;
	long colSt=includeCols.size()+1;  // add 1 extra column to store source flux density (not relative flux density)
	long colSsource=includeCols.size(); 
//	S.setSize(1+Nref,xyList.Ny()); // source flux density and Nref reference flux densities from annulus around
	S.setSize(colSt+Nref,xyList.Ny()); // source flux density and Nref reference flux densities from annulus around
	S.allocFunctionSpace();
	long N=xyList.Ny();
	long row;
	long haloDone=0;
	
	// generate/load random referece beam positions	
	if (cpeds_fileExists("angles")==false) {	
		phiAll=rns.getRNs(Nref*100); // we use 100x more Nref to accommodate for cases when not all paralactic angles will be considered
		// save random numbers for another runs
		phiAll.save("angles");
		phiAll.clear();
	}
	phiAll.load("angles");
	
	for (long i = 0; i < phiAll.size(); i++) {
		if ((phiAll[i]>_phi[0] and phiAll[i]<=_phi[1]) or (phiAll[i]>_phi[2] and phiAll[i]<=_phi[3]) ) {
			phi.push_back(phiAll[i]);
			if (phi.size()==Nref) i=phiAll.size(); // load at most Nref suitable angles
		}
	}
	
	// generate random pointing imperfections
	if (_pointingAcc>0) {
		if (cpeds_fileExists("pointingErrX")==false) {	
//			rns.setRNsType("gaussian_circle"); // we will randomly choose offset direction
//			rns.setMeanVariance(0,1);
			rns.setRNsType("uniform"); // we will randomly choose offset direction
			rns.setMinMax(-1,1); // we will randomly choose offset direction
			dX=rns.getRNs(Nref);
			dY=rns.getRNs(Nref);
			// save random numbers for another runs
			dX.save("pointingErrX");
			dY.save("pointingErrY");
		}
		dX.load("pointingErrX");
		dY.load("pointingErrY");
	}
	else {
		dX.makeLength(Nref);
		dY.makeLength(Nref);
		dX=0.0;
		dY=0.0;
	}
	
	// convert pointing offsets to pixel space
	for (long h= 0; h< Nref; h++) {
		// gaussian case
//		dX[h]=dX[h]*s*_pointingAcc; // s is already in pixel space
//		dY[h]=dY[h]*s*_pointingAcc;
//		printf("s=%lE, ptAcc: %lE, _pixSize: %lE\n",s,_pointingAcc,_pixSize);
//		exit(0);
		// uniform case
		dX[h]=dX[h]*fwhm*_pointingAcc; // fwhm is already in pixel space
		dY[h]=dY[h]*fwhm*_pointingAcc;
	}

	//
	// calculate which halos will be processed
	//
	
	// maybe in future
	
	if (_numTh!=0) {
		omp_set_num_threads(_numTh);		
	}
//	firstprivate(includeCols,colSsource, cutRadius,_refBeamSep, refBeamRegionRadius, Wpix,s,colSt,_pixSize,coli,colj,Nref,N,_haloMassThres )
//#pragma omp parallel for \
//	default (none) \
//	private(row) \
//	shared(haloDone,S,f,phi,dX,dY,xyList, _pointingAcc, \
//			includeCols,colSsource, cutRadius,_refBeamSep, refBeamRegionRadius, Wpix,s,colSt,_pixSize,coli,colj,Nref,N,_haloMassThres \			
//	) 
	for (row = 0; row < N; row++) {
		bool processHalo=true;
		double Ssrc,Sbg;
		long i,j,iOff,jOff;
		double dx,dy;
		mscsFunction3dregc beam,source;
		mscsFunction3dregc bg; // reference beam pointing at offset location _relBeamSep away in a random direction
		i=xyList(coli,row); //+dX[row % Nref];
		j=xyList(colj,row); //+dY[row % Nref];
		
		//				i-cutRadius<0 or i+cutRadius >= f.Nx() or
		//				j-cutRadius<0 or j+cutRadius >= f.Ny() or
		if (
				i-refBeamRegionRadius<0 or i+refBeamRegionRadius >= f.Nx() or
				j-refBeamRegionRadius<0 or j+refBeamRegionRadius >= f.Ny()
		) {
			processHalo=false;
		}
		
		if (xyList(32,row)< _haloMassThres) processHalo=false; // neglect the least massive halos
		
		long cutRadiusL=long(cutRadius);
		if (processHalo) {
			source=f.cutAwayBlock(i-cutRadiusL,j-cutRadiusL,0,i+cutRadiusL,j+cutRadiusL,0);	
			beam=source; // setting size
			beam.mkGauss2D(i,j,s,s,1,2,0);
//			printf("id: %li i=%li, j=%li, ist: %li jst: %li ien: %li jen: %li, cutRdius: %li\n",
//					long(xyList(0,row)),i,j,
//					i-cutRadiusL,j-cutRadiusL,
//					i+cutRadiusL,j+cutRadiusL,
//					cutRadiusL);

//			source.saveSlice(2,0,"source.tmp"); 
//			beam.saveSlice(2,0,"beam.tmp"); 
//			f.saveSlice(2,0,"f.tmp"); 
//			exit(0);
			
			source*=beam;
			// source flux density
			Ssrc=source.sumRe()*Wpix;
			
			for (long col = 0; col < includeCols.size(); col++) {
				S(col,row)=xyList(includeCols[col],row); // set the requested included columns data				
			}
			
			S(colSsource,row)=Ssrc;
			// calculate reference flux densities
			long iref;
#pragma omp parallel for default (shared) private (iref, iOff,jOff,bg,Sbg,source) firstprivate (Ssrc)
			for (iref = 0; iref < Nref; iref++) {
				iOff=i+dX[(row+iref) % Nref]+_refBeamSep*cos(phi[iref])/_pixSize;
				jOff=j+dY[(row+iref) % Nref]+_refBeamSep*sin(phi[iref])/_pixSize;
//				printf("halo: %li iref: %li, i:%li, j:%li, dX: %li dY:%li, %li %li %li %li\n",haloDone,iref,i,j,(long)dX[row],(long)dY[row],iOff-cutRadiusL,iOff+cutRadiusL,jOff-cutRadiusL,jOff+cutRadiusL);
				bg=f.cutAwayBlock(iOff-cutRadiusL,jOff-cutRadiusL,0,iOff+cutRadiusL,jOff+cutRadiusL,0);
				bg*=beam;
				// background flux density
				Sbg=bg.sumRe()*Wpix;
						
				// if pointing accuracy is > 0 then we need to recalculate source flux density as well
				if (_pointingAcc>0) {
					iOff=i+dX[(row+iref) % Nref];
					jOff=j+dY[(row+iref) % Nref];
					source=f.cutAwayBlock(iOff-cutRadiusL,jOff-cutRadiusL,0,iOff+cutRadiusL,jOff+cutRadiusL,0);
					source*=beam;
					// source flux density
					Ssrc=source.sumRe()*Wpix;
				}
				// difference flux density
				S(colSt+iref,row)=Ssrc-Sbg;
				
			}
		}
//#pragma omp critical
		{
		haloDone++;
		printf("halo %li/%li\n",haloDone,N);
		}
	}

	return S;	
}

/***************************************************************************************/
mscsFunction calculateCMBflux_vs_angSize(mscsFunction3dregc& f) {
	mscsFunction Cl,Sth,beam,sinTh,sig;
	
	double fwhm=cpeds_calculateAngularResolution(_aperture,_freq)*PI180 ; //rad
	double s=fwhm/(2.0*sqrt(2.0*log(2.0)));
	double dx=fwhm/1000;
	beam.mkGauss(0,3*fwhm,dx,1,0,s);
	sinTh.mkSin(0,3*fwhm,dx,twoPI,0,1);
	sig=beam*sinTh;
	double ang=sig.integrate()*twoPI; // 2*PI int b(theta) sin(theta) d theta
	printf("fwhm [deg]: %lf\n",fwhm);
	printf("freq [GHz]: %lf\n",_freq);
	printf("D [m]: %lf\n",_aperture);
	printf("twoPI int beam(th,D,freq)*sin(th) dth [sr]: %lE\n",ang);
	printf("thSt [deg]: %lf\n",_thSt);
	printf("thEn [deg]: %lf\n",_thEn);
	printf("dth [deg]: %lf\n",_dth);
	
	Cl.load(_mapFile);
	Cl.checkRanges();
	
	double sigmaCMB=0;
	double dl=100;
	double l;
	long lst,len,li;
	double th=_thSt;
	double Bnu0,S0,S;
	double freq=_freq*1e9;
	Bnu0=cpeds_black_body_radiance_Bnu(freq,T_CBR0);
	S0=Bnu0*ang;
	printf("Bnu0 [Jy/sr]: %lE\n",Bnu0);
	printf("S0 [Jy]: %lE\n",S0);
	
	len=Cl.getMaxArg();
	while (th<_thEn) {
		l=round(0.5*(360.0/th-1.0));
		//		lst=Cl.getx
		sigmaCMB=0;
		lst=l;
		if (_exactl) {
			li=lst;
			if (li>=Cl.getMinArg() and li<=Cl.getMaxArg())
				sigmaCMB=Cl.f(li)*(2*li+1);
			else
				sigmaCMB=0;				
		}
		else {
			for (li = lst; li <= len; li++) {
				if (li>=Cl.getMinArg() and li<=Cl.getMaxArg())
					sigmaCMB+=Cl.f(li)*(2*li+1);
				else
					sigmaCMB=0;
			}
			
		}
		sigmaCMB=sqrt(sigmaCMB/fourPI); // rms fluctuation value below angular scale th
		
		// integrate CMB flux per beam
		//		2 \[Nu]^2/c^2  h \[Nu]/(Exp[h \[Nu]/(kB T)] - 1)/Jy * 2\[Pi] Integrate[	      beam[\[Theta], \[Theta]rtRTH] Sin[\[Theta]], {\[Theta], 0, 	        3\[Theta]rtRTH}];
		//		printf("f [GHz]: %f,  Bnu [Jy/sr]: %lE\n",)
		S=fabs(S0-cpeds_black_body_radiance_Bnu(freq,T_CBR0+sigmaCMB)*ang)*1000; // store in mJy
		
		Sth.newPoint(th,S ); 
		
		
		th+=_dth;
		_dth*=_geo;
	}
	
	return Sth;
	
}

/***************************************************************************************/
mscsFunction3dregc averageRelativeFluxDensities(mscsFunction3dregc& S, long colSt) {
	mscsFunction3dregc Savg;
	mscsFunction Srel;

	// calculate the processed halo count
	long count=0;
	for (long j = 0; j < S.Ny(); j++) {		if (S(0,j)!=0) count++;	}

	Savg.setSize(colSt+2,count);
	Savg.allocFunctionSpace();
	count=0;
	for (long j = 0; j < S.Ny(); j++) {
		if (S(0,j)!=0) {
			Srel=S.getSlice1Dfn(0,j,0,0); // get row
			for (long col = 0; col < colSt; col++) {
				Savg(col,count)=S(col,j); // copy 'included' columns
				Srel.deletePoint(0); // remove 'included' columns from averaging
			}
			Savg(colSt,count)=Srel.meanf();
			Savg(colSt+1,count)=Srel.stdev();
	//		Savg(colSt+2,j)=Srel.stdev();	
			count++;
		}
	}
	return Savg;
}
