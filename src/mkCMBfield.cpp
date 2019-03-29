/*!
 * \file mkCMBfield.cpp
 *
 *  Project: ocraToolkit
 *  Created on: Mar 18, 2013 3:50:26 PM
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
//#include "qtbase.h"

#include "cpeds-msgs.h"
#include "cpeds-math.h"
#include "Mscs-function3dregc.h"
#include "Mscs-map.h"
#include "cpeds-project.h"


#ifndef _NO_NAMESPACE
using namespace std;
//using namespace math;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif

string _programVersionString;
string _outfile, _input_file, _outDir, _freq; //, _noise_file;
double _field_l,_field_b, _field_size, _dipoleAmp, _dipolel, _dipoleb, _monopoleAmp, _NneighTol;
double _fieldMarginFactor;
double _aperture;
long _Npix,_Nneigh,_NneighMax;
bool _dipole, _monopole, _convolveBeam, _save_all;// _addNoise;


void parseOptions(int argc, char** argv);
string getProgramVersionString();


mscsFunction3dregc getRawField(mscsFunction3dregc& f,mscsMap& map);
/*!
	\brief selects a direction set from the full sphere CMB simulation 
	\details 
	@param map - reference to the map from which the directions are taken
	@param msgs - messanger object
	@return a direction set selected from the full sphere CMB simulation located within angular distance smaller than the field size defined by the _field_size parameter.
	The actular direction set is created from larger spherical area in order to make SPH interpolations work correctly on the edges of the field.
	

	\date Mar 19, 2013, 1:58:08 PM
	\author Bartosz Lew
 */
cpedsDirectionSet selectDirections(mscsMap& map, bool inRectField, cpedsMsgs& msgs);
/*!
	\brief projects the selected directions onto a tangent plane using stareographic projection
	\details 
	@param ds - direction set to be projected [rad]
	@param pixSize - reference at which a linear size on the projection plane of the pixel will be assigned (based on _Npix and _field_size parameters).
	The pixel size is defined on sphere as _field_size/_Npix - the pixSize quantity is to be understood in the flat sky approximation as an projected linear size along a big circle.
	@param msgs - messanger object
	@return projected points set.

	\date Mar 19, 2013, 2:01:16 PM
	\author Bartosz Lew
 */
cpedsPointSet3D projectOnPlane(cpedsDirectionSet& ds, double& pixSize, string fname, cpedsMsgs& msgs);
/*!
	\brief generates interpolated CMB field on a tangent plane
	\details 
	@param ps - porjected point set
	@param pixSize - single pixel size on the tangent plane
	@param ds - direction set (keeps info on function values)
	@param msgs - messanger object
	@return SPH interpolated field

	\date Mar 19, 2013, 2:09:56 PM
	\author Bartosz Lew
 */
mscsFunction3dregc makeInterpolatedField(cpedsPointSet3D& ps, double pixSize, cpedsMsgs& msgs);
void saveDirectionSet(cpedsDirectionSet ds, string fname, cpedsMsgs& msgs);
cpedsDirectionSet projectOnSphere(mscsFunction3dregc& f,cpedsMsgs& msgs);
void addDipoleComponent(mscsMap& map);
cpedsDirectionSet mkDipoleDirections(cpedsDirectionSet& ds,cpedsMsgs& msgs);
void addMonopoleComponent(mscsMap& map);

/*!
	\brief convolves input function with gaussian beam on the tangent plane
	\details 
	@param freq - frequency for which the theoretical gaussian beam will be derived [GHz]
	@param aperture - Cassegrain telescope apterture [m]
	@param pixSize - pixel size on the tangent plane [arbitrary units]
	@param fieldSize - field size on tangent plane [deg]
	@return convoled field
	
	this approach has a drawback because CMB in the field is not periodic and fft introduces 
	some edge-effects (but very minor). In any case it is better the convolve the full sky map with intrumental beam 
	in SH space. THerefore this function is not used.
	
	BTW. smoothing CMB with ocra beam has a negligible effects on the fluctuations because there are almost no 
	fluctuations at these angular scales so the smoothed map looks almost the same as the input map

	\date Apr 17, 2013, 11:03:26 PM
	\author Bartosz Lew
 */
mscsFunction3dregc convolveWithGaussianBeam(mscsFunction3dregc& f, double pixSize, double fieldSize, double freq, double aperture, cpedsMsgs& msgs);

mscsMap convolveWithGaussianBeam(mscsMap& f, double freq, double aperture, cpedsMsgs& msgs);

int main(int argc, char **argv) {
	
	cpedsMsgs msgs("mkCMBfield",false,"",High);
	msgs.setSaveRunWriteMode('a');
	msgs.saveThisRun(argc,argv);
	//----------------------------------------------------------------------------------------------------
	//PARSE THE COMMAND LINE ARGUMENTS
	//----------------------------------------------------------------------------------------------------
	parseOptions(argc,argv);
	msgs.setLogFileName(_outDir+"/"+"mkCMBfield-readme.txt");
	msgs.loggingOn();
	string frequencySuffix;
	QString freqstr=_freq.c_str();
	QStringList qfl=freqstr.split(",");
	cpedsList<double> freq;
	bool ok;
	double num;
	double frequency;
	for (unsigned long i = 0; i < qfl.size(); i++) {
		num=qfl[i].toDouble(&ok);
		if (ok) freq.append(num);
		else msgs.criticalError("could not convert frequency do double. Better stop now.",Top);
	}
	

	
	//
	// make output directory
	//
	if (_outDir=="")
		_outDir=_input_file+".data";
	mkdir(_outDir.c_str(),0755);
	

	//
	// load CMB map
	//
	
	mscsMap map, mapb;
	map.loadbinT(_input_file);
	map.set_map_coord();

	//
	// convolve with gaussian beam
	//
	if (_convolveBeam) {
		mapb=convolveWithGaussianBeam(map,frequency,_aperture,msgs);
		mapb.savebinT("cmb-beam-Tn-bin");
		map=mapb;
	}

//	if (_addNoise) {
//		mscsMap noise;
//		noise.loadbinT(_noise_file);
//		noise.set_map_coord();
//	}

	//
	// extract list of directions from the map within the requested field
	//
	cpedsDirectionSet dsf = selectDirections(map,true,msgs);
	if (_save_all) saveDirectionSet(dsf,"selected_dirs_rect_field.ds",msgs);
	
	cpedsDirectionSet ds = selectDirections(map,false,msgs);
	msgs.say("Found %li directions (including the margin space)",ds.size(),High);
	if (_save_all) saveDirectionSet(ds,"selected_dirs.ds",msgs);
	
	//
	// project direction set on a tangent plane
	//
	double pixSize;
	cpedsPointSet3D psf=projectOnPlane(dsf,pixSize,_outDir+"/"+"projected_rect_field.ps",msgs);
	cpedsPointSet3D ps =projectOnPlane(ds,pixSize,_outDir+"/"+"projected.ps",msgs);
	//	cpedsPointSet3D psll =projectOnPlane(dsll,pixSize,_outDir+"/"+"projected_lowell.ps",msgs);
	
	//
	// make field
	//
	mscsFunction3dregc f=makeInterpolatedField(ps,pixSize,msgs);
	mscsFunction3dregc fSum=f;
	mscsFunction3dregc DTcmb=f;
	f.saveHDF5(_outDir+"/"+_outfile+".hdf5","DTcmb");
	//	f.saveSlice(2,0,_outDir+"/DTcmb.mat",0);
	f.setHDF5_scalarDoubleAttribute(_outDir+"/"+_outfile+".hdf5","DTcmb","field_l",_field_l,"field center galactic longitude [deg]");
	f.setHDF5_scalarDoubleAttribute(_outDir+"/"+_outfile+".hdf5","DTcmb","field_b",_field_b,"field center galactic latitude [deg]");
	f.setHDF5_scalarDoubleAttribute(_outDir+"/"+_outfile+".hdf5","DTcmb","field_sizeX",_field_size,"full with of the field [deg]");
	f.setHDF5_scalarDoubleAttribute(_outDir+"/"+_outfile+".hdf5","DTcmb","field_sizeY",_field_size,"full with of the field [deg]");

//	if (_save_all) {
		mscsFunction3dregc g=getRawField(f,map);
		g.saveHDF5(_outDir+"/"+_outfile+".hdf5","DTcmbCoarse");
		
//	}
	
	//
	// iterate over frequencies
	//
	for (unsigned long i = 0; i < freq.size(); i++) {
		frequency=freq[i];
		msgs.say("PROCESSING FREQUENCY: %lf GHz",frequency,High);
		f=DTcmb;
		fSum=f;
		frequencySuffix="-"+msgs.toStrf(frequency,1)+"GHz";
		//
		// precalculate conversion factor from K to Jy/sr
		//
		double Tcmb=_monopoleAmp;
		double nu=frequency*1.0e9; // conversion from GHz to Hz
		double A=2.0*(CPEDS_kB*Tcmb)*(CPEDS_kB*Tcmb)*(CPEDS_kB*Tcmb)/( (CPEDS_h*CPEDS_c)*(CPEDS_h*CPEDS_c) );
		double x=CPEDS_h*nu/(CPEDS_kB*Tcmb);
		double B=(x*x*x*x*exp(x))/( (exp(x)-1.0)* (exp(x)-1.0) );
		double dT2JyoSr=A*B/Tcmb/CPEDS_Jy; // conversion factor from K to Jy
		double Bnu=2.0*CPEDS_h*nu*nu*nu/(CPEDS_c*CPEDS_c)/(exp(x)-1)/CPEDS_Jy; // CMB monopole spectral radiance [Jy/sr]
		
		
		
		f*=dT2JyoSr; 
		//	f.saveSlice(2,0,_outDir+"/DLcmb.mat",0);
		f.saveHDF5(_outDir+"/"+_outfile+".hdf5","DLcmb"+frequencySuffix);
		f.setHDF5_scalarDoubleAttribute(_outDir+"/"+_outfile+".hdf5","DLcmb"+frequencySuffix,"field_l",_field_l,"field center galactic longitude [deg]");
		f.setHDF5_scalarDoubleAttribute(_outDir+"/"+_outfile+".hdf5","DLcmb"+frequencySuffix,"field_b",_field_b,"field center galactic latitude [deg]");
		f.setHDF5_scalarDoubleAttribute(_outDir+"/"+_outfile+".hdf5","DLcmb"+frequencySuffix,"dT2JyoSr",dT2JyoSr,"Kelvin to Jansky/sr conversion factor for this frequency");
		f.setHDF5_scalarDoubleAttribute(_outDir+"/"+_outfile+".hdf5","DLcmb"+frequencySuffix,"freq",frequency,"map frequency [GHz]");
		
		//	exit(0);
		//
		// add monopole component
		//
		
		if (_monopole) {
			f=double(_monopoleAmp);
			//		f.saveSlice(2,0,_outDir+"/Tmonopole.mat",0);
			if (i==0) {
				f.saveHDF5(_outDir+"/"+_outfile+".hdf5","Tmonopole");
				f.setHDF5_scalarDoubleAttribute(_outDir+"/"+_outfile+".hdf5","Tmonopole","monopole_amp",_monopoleAmp,"dipole amplitude [K]");
			}
			f=double(Bnu);
			//		f.saveSlice(2,0,_outDir+"/Lmonopole.mat",0);
			f.saveHDF5(_outDir+"/"+_outfile+".hdf5","Lmonopole"+frequencySuffix);
			f.setHDF5_scalarDoubleAttribute(_outDir+"/"+_outfile+".hdf5","Lmonopole"+frequencySuffix,"Bnu",Bnu,"CMB monopole radiance [Jy/sr]");
			//		fSum+=f;
		}
		//
		// project interpolated function values back on sphere 
		//
		ds=projectOnSphere(f,msgs);
		if (i==0 and _save_all) saveDirectionSet(ds,"DTcmb.ds",msgs);
		 
		//
		// add dipole component
		//
		if (_dipole) {
			cpedsDirectionSet ds_dipole;
			ds_dipole = mkDipoleDirections(ds,msgs); // the values ordering in ds is the same is in function object f
			cpedsPointSet3D ps_dipole =projectOnPlane(ds_dipole,pixSize,_outDir+"/"+"projected_dipole.ps",msgs);
			for (long i = 0; i < ds.size(); i++) {
				ds[i].value()+=ds_dipole[i].val();
			}
			if (i==0) {
				if (_save_all) saveDirectionSet(ds,_outfile+"_interpolated_all.ds",msgs);
			}
			for (long i = 0; i < ds_dipole.size(); i++) {
				f.fRe(i)=ds_dipole[i].val();
			}
			//		f.saveSlice(2,0,_outDir+"/"+_outfile+"-Tdipole.mat",0);
			f.saveHDF5(_outDir+"/"+_outfile+".hdf5","Tdipole");
			f.setHDF5_scalarDoubleAttribute(_outDir+"/"+_outfile+".hdf5","Tdipole","dipole_l",_dipolel,"dipole axis galactic longitude [deg]");
			f.setHDF5_scalarDoubleAttribute(_outDir+"/"+_outfile+".hdf5","Tdipole","dipole_b",_dipoleb,"dipole axis galactic latitude [deg]");
			f.setHDF5_scalarDoubleAttribute(_outDir+"/"+_outfile+".hdf5","Tdipole","dipole_amp",_dipoleAmp,"dipole amplitude [K]");
			fSum+=f;
			//		fSum.saveSlice(2,0,_outDir+"/"+_outfile+"-DTcmb_all.mat",0);
			fSum.saveHDF5(_outDir+"/"+_outfile+".hdf5","DTcmb_all");
			f*=dT2JyoSr; 
			//		f.saveSlice(2,0,_outDir+"/"+_outfile+"-Ldipole.mat",0);
			f.saveHDF5(_outDir+"/"+_outfile+".hdf5","Ldipole"+frequencySuffix);
		}
		
		//
		// convert to spectral radiance units [Jy/sr]
		//
		fSum*=dT2JyoSr; // this is DeltaS = S-S0 = dBnu/dT * (DeltaT)
		printf("temperature to Jansky/sr conversion factor: %lE\n",dT2JyoSr);
		//	printf("dT2Jy: %lE\n",dT2Jy);
		//	fSum.saveSlice(2,0,_outDir+"/"+_outfile+"-DLcmb_all.mat",0);
		fSum.saveHDF5(_outDir+"/"+_outfile+".hdf5","DLcmb_all"+frequencySuffix);
		fSum.setHDF5_scalarDoubleAttribute(_outDir+"/"+_outfile+".hdf5","DLcmb_all"+frequencySuffix,"dT2JyoSr",dT2JyoSr,"Kelvin to Jansky/sr conversion factor for this frequency");
		
		
		//
		// convert to flux density units [Jy] 
		//
		double deg2tosr=pow(PI180,2);
		double Wpix=pow(_field_size/_Npix,2)*deg2tosr; /*
		 * Comment:
		 *  The pixels in the spherical rectangle have different areas and so these factors should be pixel dependent
		 * 
		 * author: blew
		 * date: Nov 8, 2013 7:29:54 PM
		 *
		 */
		
		printf("Wpix: %lE\n",Wpix);
		fSum*=Wpix;
		//	fSum.saveSlice(2,0,_outDir+"/"+_outfile+"-DFcmb_all.mat",0);
		fSum.saveHDF5(_outDir+"/"+_outfile+".hdf5","DFcmb_all"+frequencySuffix);

		
		msgs.say("Files maming conventions:",High);
		msgs.say("DTcmb/DFcmb/DLcmb - means that the file contains cmb temperature fluctuations (DT=T-Tcmb0)/ spectral flux density (DF=F-Fcmb0)/radiance respectively (DL=L-Lcmb0)",Medium);
		msgs.say("L/T/F - means that the file contains radiance ,temperature or spectral flux density respectively ",Medium);
		msgs.say("_beam.mat - contains the field convolved with gaussian beam",Medium);
		msgs.say("Units used",High);
		msgs.say("temperature units T: [K]",Medium);
		msgs.say("spectral flux density F: [Jy]",Medium);
		msgs.say("radiance L: [Jy/sr]",Medium);
		msgs.say("Parameters used for this run: ",High);
		msgs.say("Tcmb0 [K]: %lf",_monopoleAmp,Medium);
		msgs.say("Dipole [K]: %lE",_dipoleAmp,Medium);
		msgs.say("Dipole l [deg] %lf: ",_dipolel,Medium);
		msgs.say("Dipole b [deg] %lf: ",_dipoleb,Medium);
		msgs.say("Frequency [GHz]: %lf",_freq,Medium);
		msgs.say("Aperture [m]: %lf",_aperture,Medium);
		msgs.say("'*all' files contain monopole component: ",_monopole,Medium);
		msgs.say("'*all' files contain dipole component: ",_dipole,Medium);
		msgs.say("'*all' beam convolution was done: ",_convolveBeam,Medium);
		msgs.say("Temperature to Jansky/sr conversion factor [Jy/sr/K]: %lE",dT2JyoSr,Medium);
		msgs.say("Bnu [Jy/sr]: %lE",Bnu,Medium);
		msgs.say("Bnu - spectral CMB monopole radiance",Medium);
		msgs.say(string("The main products are saved in: ")+_outfile+"-D*cmb_all.mat",High);

	}
	
	
	
	
	return 0;
}

/***************************************************************************************/
void parseOptions(int argc, char** argv) {
	long i;
	_programVersionString=getProgramVersionString();
	try {
		
		//     CmdLine cmd("RT4\n\n. ",' ', "0.1" );
		CmdLine cmd("mkCMBfield\nTorun Centre for Astronomy, UMK, POLAND.\n "
				"This program takes input CMB map in healpix nested binary format and produces an SPH interpolated temperature fluctuations field towards requested direction in the sky."
				"The output is given in the ascii matrix form of the requested size\n\n. "
				"\n"
				"example usage: "
				"mkCMBfield SZbg-0-sig-WMAP-Tn-bin -o cmb-l_90-b_85 --size 5 --outdir field-5x5 --lon 90 --lat 85 --Npix 256 --Nneigh 30 --NneighTol 0 --monopole --dipole"
				"",' ', _programVersionString.c_str() );
		//     CmdLine cmd("RT4 azimutal scan simulator for OCRA-p\n", 'This program simulates the RT4 scan motion, sky rotation during \n\n. ', 0.1 );
		
		// 
		// Define arguments
		//
		
		//		ValueArg<double> field_l("l", "", "galactic longitude of the field center  [deg] (default: 0)", false,0,"double");	cmd.add( field_l );
		//		ValueArg<double> field_b("b", "", "galactic latitude of the field center [deg] (default: 0)", false,0,"double"); cmd.add( field_b );
		ValueArg<double> field_l("", "lon", "galactic longitude of the field center  [deg] (default: 90)", false,90,"double");	cmd.add( field_l );
		ValueArg<double> field_b("", "lat", "galactic latitude of the field center [deg] (default: 0)", false,0,"double"); cmd.add( field_b );
		ValueArg<double> field_size("s", "size", "size of the field [deg]. The field is defined as an spherical rectangle with center at (l,b) and "
				"with the full field span of 'size'. (default: 5)", false,5,"double"); cmd.add( field_size );
		ValueArg<long> Npix("", "Npix", "number of pixelx in the field along side (1024)", false,512,"long");	cmd.add( Npix );
		ValueArg<long> Nneigh("", "Nneigh", "minimal number of neighbors for SPH smoothing (300)", false,300,"long");	cmd.add( Nneigh );
		ValueArg<double> NneighTol("", "NneighTol", "neighbors number tolerance for SPH smoothing (5 \%)", false,5,"double");	cmd.add( NneighTol );
		
		SwitchArg monopole("","monopole", "enable monopole component (default: false)", false);	cmd.add( monopole);
		ValueArg<double> monopoleAmp("", "mAmp", "monopole compoment amplitude. [thermodynamic K]", false,2.726,"double"); cmd.add( monopoleAmp );
		SwitchArg dipole("","dipole", "enable dipole component (default: false)", false);	cmd.add( dipole );
		ValueArg<double> dipoleAmp("", "dAmp", "dipole compoment amplitude. [thermodynamic K]", false,0.003346,"double"); cmd.add( dipoleAmp );
		ValueArg<double> dipolel("", "dl", "dipole compoment direction l [deg]", false,263.85,"double"); cmd.add( dipolel );
		ValueArg<double> dipoleb("", "db", "dipole compoment direction b [deg]", false,48.25,"double"); cmd.add( dipoleb );
		ValueArg<string> freq("f", "freq", "comma separated list of frequencies for which to calculate cmb maps [GHz] (default: 15,22,30)", false,"15,22,30","string"); cmd.add( freq );
		ValueArg<double> aperture("D", "aperture", "telescope aperture. Used to derive theoretical gaussian beam size [m] (default: 32)", false,32,"double");	cmd.add( aperture );
		SwitchArg convolveBeam("","convolveBeam", "enable theoretical gaussian beam convolution (default: false)", false);	cmd.add( convolveBeam );
		SwitchArg save_all("","save_all", "if used then much more outputs is saved. Useful to not to use for large simulations. (default: false)", false);	cmd.add( save_all );
//		SwitchArg addNoise("","addNoise", "add noise to the generated m(default: false)", false);	cmd.add( convolveBeam );

		
		ValueArg<double> fieldMarginFactor("", "fieldMarginFactor", "the projected field will be larger by this factor for make SPH interpolations "
				"behave better within the field. (default: 0.1 - increases the field size by 10 per-cent)", false,0.1,"double");	cmd.add( fieldMarginFactor );
		
		
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
		
		UnlabeledValueArg<string> input_file("input_file","file to plot",true,"nofile","string");	cmd.add( input_file );
		//		ValueArg<string> mask("m","mask","wheather to mask the file before doing statistics (prefix)",false,"","string"); cmd.add(mask);
		//     SwitchArg mask_from_here("","MM", "take mask from current directory rather than from the default directory", false);	cmd.add( mask_from_here );
		//     SwitchArg dont_save_mask("","dont_save_mask", "do not print to file masked pixels; useful for -Tn-txt savings", false);	cmd.add( dont_save_mask );
		ValueArg<string> outfile("o","outfile","name of the output file (default: cmb)",false,"cmb","string"); cmd.add(outfile);
		ValueArg<string> outdir("","outdir","name of the output directory (default: )",false,"","string"); cmd.add(outdir);
		//     ValueArg<string> ft("f","ft","input file type [bin, txt] (default: bin)",false,"bin","string"); cmd.add(ft);
		//     ValueArg<string> ord("","ord","ordering of the map [ring, nest] default: nest",false,"nest","string"); cmd.add(ord);
		//     SwitchArg mink("M","Mink", "plot also the minkowski functionals for the map", false);	cmd.add( mink );
		//     ValueArg<long> mlevels("", "mlev", "number of levels to calculate the minkowski functionals (default: 100)", false,100,"long");	cmd.add( mlevels );
		
		cmd.parse(argc,argv);
		
		//
		// Set variables
		//
		_Npix = Npix.getValue();
		_field_size=field_size.getValue();
		_field_l = field_l.getValue();
		_field_b = field_b.getValue();
		_Nneigh=Nneigh.getValue();
		_NneighTol=NneighTol.getValue();
		_NneighMax=_Nneigh+long(double(_Nneigh)*_NneighTol/100);
		
		_outfile = outfile.getValue();
		_input_file=input_file.getValue();
		_outDir = outdir.getValue(); 	if (_outDir=="") _outDir=".";
		
		
		_monopole=monopole.getValue();
		_monopoleAmp=monopoleAmp.getValue();
		_dipole=dipole.getValue();
		_dipoleAmp=dipoleAmp.getValue();
		_dipolel=dipolel.getValue();
		_dipoleb=dipoleb.getValue();
		_freq=freq.getValue();
		_aperture=aperture.getValue();
		_convolveBeam=convolveBeam.getValue();
		_save_all=save_all.getValue();
		_fieldMarginFactor=fieldMarginFactor.getValue();
		
		
	} catch ( ArgException& e )
	{ cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
	
	
}
/***************************************************************************************/
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
cpedsDirectionSet selectDirections(mscsMap& map, bool inRectField, cpedsMsgs& msgs) {
	//	
	// calculate filed
	//
	double r=_field_size;
	
	double fieldSizeFactor=_fieldMarginFactor; // this gives the total filed size 30x the scan field size
	//	msgs.say("POINT SET RANGES:",High);
	//	printf("minx: %lE maxx: %lE, miny: %lE, maxy: %lE, zmin: %lE, zmax: %lE\n",minx,maxx,miny,maxy,minz,maxz);
	r+=r*fieldSizeFactor;
	
	//	sizeX=maxx-minx;
	//	sizeY=maxy-miny;
	
	double l1,l2,b1,b2;
	l1=_field_l-_field_size/2;
	l2=_field_l+_field_size/2;
	b1=_field_b-_field_size/2;
	b2=_field_b+_field_size/2;
	
	l1*=PI180;
	l2*=PI180;
	b1*=PI180;
	b2*=PI180;
	
	cpedsDirection lb(_field_l,_field_b);
	cpedsDirection n;
	lb*=PI180;
	r*=PI180;
	cpedsDirectionSet ds;
	for (long i = 0; i < map.pixNum(); i++) {
		n=map.n(i);
		if (inRectField) {
			if (n.lon()>=l1 and n.lon()<l2 and n.lat()>=b1 and n.lat()<b2) {
				n.setVal(map.T(i));
				ds.append(n);
			}
		}
		else {
			if (lb.angle(n)<r) {
				n.setVal(map.T(i));
				ds.append(n);
			}
		}
	}
	return ds;
}
/***************************************************************************************/
cpedsPointSet3D projectOnPlane(cpedsDirectionSet& ds, double& pixSize, string fname,cpedsMsgs& msgs) {
	msgs.say("PROJECTING ON PLANE",High);
	//
	// project on the central direction
	//
	cpedsDirection n(_field_l,_field_b);
	cpedsProject pro(ds,"stere");
	pro.setInvertLon(false);
	cpedsPointSet3D ps=pro.projectOnPlane(n,"+proj=%s +lat_0=%lE +lon_0=%lEE +R=%lE");
	//	ps.printRanges("projected ps0");
	if (_save_all) ps.save(fname);
	double dang=_field_size/_Npix; //deg -- defines the resolution of the field
	pixSize=pro.getLengthScaleFromAng(dang*PI180);
	
	msgs.say("LINEAR PIXEL SIZE IS: %lE\n",pixSize,High);
	
	return ps;
	
}
/***************************************************************************************/
mscsFunction3dregc makeInterpolatedField(cpedsPointSet3D& ps, double pixSize, cpedsMsgs& msgs) {
	msgs.say("INTERPOLATING FIELD ON GRID CELLS",High);
	
	mscsFunction3dregc f;
	double minx,maxx,miny,maxy,minz,maxz;
	ps.getRanges(minx,maxx,miny,maxy,minz,maxz);
	
	subDomain_region_t sd;
	subDomain_region_t td;
	//
	// define tree domain
	//
	td.xmin=minx;
	td.xmax=maxx;
	td.ymin=miny;
	td.ymax=maxy;
	td.zmin=0;
	td.zmax=0;
	td.subx=2;
	td.suby=2;
	td.subz=1;
	
	//
	// define function domain
	//
	ps.printRanges("projected ps");
	double fSize=pixSize*_Npix/2.0;
	
	/*
	 * Comment: this below definition of the function domain - defined in the projective plane - effectively translates onto
	 * a spherical region that has the same angular dimensions all over the sphere, regardless of the requested spherical rectangle
	 * field definition according to which field directions are selected in the main program.
	 * 
	 * author: blew
	 * date: Mar 20, 2013 10:57:49 PM
	 *
	 */
	sd.xmin=-fSize;
	sd.xmax= fSize;
	sd.ymin=-fSize;
	sd.ymax= fSize;
	
	sd.zmin=0;
	sd.zmax=0;
	sd.subx=_Npix;
	sd.suby=_Npix;
	sd.subz=1;
	
	mscsVector<double> vals;
	vals.setSize(ps.size());
	for (long i = 0; i < ps.size(); i++) {
		vals[i]=ps[i].z();
		ps[i].setZ(0);
	}
	//
	// interpolate filed
	//
	ps.printRanges("projected ps before interpolation");
	f.mkInterpolatedFieldScatter(sd,td,ps,vals,"gadget2",_Nneigh,_NneighMax);
	f.printInfo();
	if (_save_all) {
		f.saveSlice(2,0,_outDir+"/interpolated.mat",0);
		f.saveSlice(2,0,_outDir+"/interpolated-Neff.mat",1);
	}
	return f;
}
/***************************************************************************************/
void saveDirectionSet(cpedsDirectionSet ds, string fname, cpedsMsgs& msgs) {
	msgs.say("saving direction set to file: ",_outDir+"/"+fname,High);
	for (long i = 0; i < ds.size(); i++) {
		ds[i]*=double(PI180inv); // convert to deg for saving
	}
	ds.save(_outDir+"/"+fname,"",true);
}
/***************************************************************************************/
void addDipoleComponent(mscsMap& map) {
	mscsMap dipole;
	cpedsDirection m(_dipolel,_dipoleb);
	m.setVal(_dipoleAmp);
	dipole.set_nside(map.nside());
	dipole.makeDipole(m);
	
	map.T()+=dipole.T();
	dipole.savebinT(_outDir+"/"+"dipole-Tn-bin");
	map.savebinT(_outDir+"/"+"CMB_and_dipole-Tn-bin");
	
}
/***************************************************************************************/
void addMonopoleComponent(mscsMap& map) {
	map.T()+=double(_monopoleAmp);
}

/***************************************************************************************/
cpedsDirectionSet projectOnSphere(mscsFunction3dregc& f, cpedsMsgs& msgs) {
	msgs.say("PROJECTING ON SPHERE",High);
	cpedsPointSet3D ps;
	//
	// project on the central direction
	//
	cpedsDirection n(_field_l,_field_b);
	
	/*
	 * Comment: this should be changed for a single loop over all function values in case 
	 * when something is changed (unlikely) internally in the ordering of the data inside of 
	 * the mscsFunction3dregc class.
	 * 
	 * author: blew
	 * date: Apr 17, 2013 4:59:36 PM
	 *
	 */
	
	n.print_direction("central projection direction",true,2);
	ps.clear();
	cpedsDirectionSet ds;
	double sigMax=f.fRe(0,0,0);
	for (long i = 0; i < f.Nx(); i++) {
		for (long j = 0; j < f.Ny(); j++) {
			ps.append(cpedsPoint3D(f.getX(i),f.getY(j),f.fRe(i,j,0)));
		}
	}
	cpedsProject pro(ps,"stere");
	pro.setInvertLon(false);
	//	pro.points()=ps;
	ds=pro.projectOnSphere(n,"+proj=%s +lat_0=%lE +lon_0=%lEE +R=%lE");				
	return ds;
}
/***************************************************************************************/
cpedsDirectionSet mkDipoleDirections(cpedsDirectionSet& ds,cpedsMsgs& msgs) {
	cpedsDirectionSet dipole_ds;
	cpedsDirection m(_dipolel,_dipoleb);
	m*=PI180;
	cpedsDirection p;
	for (long i = 0; i < ds.size(); i++) {
		p=ds[i];
		p.setVal(_dipoleAmp*cos(p.angle(m)));
		dipole_ds.append(p);
	}
	if (_save_all) saveDirectionSet(dipole_ds,"dipole.ds",msgs);
	
	return dipole_ds;
}
/***************************************************************************************/
//mscsFunction3dregc convolveWithGaussianBeam(mscsFunction3dregc& f, double pixSize, double fieldSize, double freq, double aperture, cpedsMsgs& msgs) {
//	
//	double D=aperture;
//	double lambda=CPEDS_c/(freq*1.0e9);
//	double fwhm=1.22*lambda/D * PI180inv; // deg
//	msgs.say("fwhm [deg]: %lf\n",fwhm,Medium);	
//	fwhm*=double(pixSize*f.Nx()/fieldSize); // convert to fwhm in units on the tangent plane
//	msgs.say("fwhm [tangent plane length units]: %lf\n",fwhm,Medium);	
//	msgs.say("lambda [cm]: %lf\n",lambda*100, Medium);
//
//	mscsFunction3dregc fsm;
//	
//	fsm=f;
//	double s=fwhm/(2.0*sqrt(2.0*log(2.0)));
//	fsm.smooth3DGauss(s,s,s);
//	
//	return fsm;
//}
/***************************************************************************************/
mscsMap convolveWithGaussianBeam(mscsMap& map, double freq, double aperture, cpedsMsgs& msgs) {
	
	double D=aperture;
	double lambda=CPEDS_c/(freq*1.0e9);
	double fwhm=1.22*lambda/D; // rad
	msgs.say("fwhm [deg]: %lf\n",fwhm * PI180inv,Medium);	
	msgs.say("lambda [cm]: %lf\n",lambda*100, Medium);
	
	mscsAlms a;
	mscsMap mapb;
	mscsWindowFunction bl,pixtf;
	bl.mkConst(0,10000,1,0);
	pixtf.mkConst(0,10000,1,1);
	bl.make_gaussian_kernel(fwhm);
	long lmax=map.nside()*2;
	a=map.SH_analysis(lmax);
	mscsAngularPowerSpectrum Cl=a.get_Cl(0,lmax,0);
	//	Cl.save("Cl");
	mapb.set_nside(map.nside());
	//	mapb.SH_synthesis(a,lmax,pixtf,pixtf,0);
	mapb.SH_synthesis(a,lmax,bl,pixtf,0);
	
	string fname=string("beam-lmax_")+msgs.toStr(lmax)+"-fwhm_"+msgs.toStrf(fwhm*PI180inv,3);
	bl.save(fname);
	
	return mapb;
}
/***************************************************************************************/
mscsFunction3dregc getRawField(mscsFunction3dregc& f,mscsMap& map) {
	printf("calculating raw field\n");
	mscsFunction3dregc g=f;
	cpedsPointSet3D ps=f.exportAsPoints(false,1.0e120);
//	string fname=_outDir+"/"+"exported.ps";
//	ps.save(fname);
	//
	// project on the central direction
	//
	cpedsDirection n(_field_l,_field_b);
//	cpedsDirection n(0.0,0.0);
	cpedsProject pro(ps,"stere");
	pro.setInvertLon(false);
	printf("projecting on spehre\n");
	cpedsDirectionSet ds=pro.projectOnSphere(n,"+proj=%s +lat_0=%lE +lon_0=%lEE +R=%lE");
	//	ps.printRanges("projected ps0");
//	fname=_outDir+"/"+"backprojected.ds";
//	ds.save(fname);

	long ii=0;
	for (long i = 0; i < f.Nx(); i++) {
		for (long j = 0; j < f.Ny(); j++) {
			for (long k = 0; k < f.Nz(); k++) {
				g.fRe(i,j,k)=map.get_T(ds[ii].lon(),ds[ii].lat());
				ii++;
			}
		}
	}
	
	return g;
	
}
