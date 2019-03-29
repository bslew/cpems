/*!
 * \file matrix_oper.cpp
 *
 *  Project: Mscs
 *  Created on: Nov 6, 2013 11:34:12 PM
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

vector<string> _infiles;
vector<string> _dsets;
vector<string> _dsetsPtSrc;
string _outFile, _operation, _hdf5dset, _freq, _hdf5_SZcomponent_output_dset_prefix;
bool _hdf5,_resizeCMB;
double _val,_val2;

void parseOptions(int argc, char** argv);
string getProgramVersionString();

mscsFunction3dregc make_TcmbtszMap(cpedsMsgs& msgs);
mscsFunction3dregc average_matrices(cpedsMsgs& msgs);
mscsFunction3dregc dumpNonZero(cpedsMsgs& msgs);
mscsFunction3dregc replaceBelow(cpedsMsgs& msgs);

int main(int argc, char **argv) {
	
	cpedsMsgs msgs(".matrix_oper",false,"",High);
	msgs.setSaveRunWriteMode('a');
	msgs.saveThisRun(argc,argv);
	string s;
	//----------------------------------------------------------------------------------------------------
	//PARSE THE COMMAND LINE ARGUMENTS
	//----------------------------------------------------------------------------------------------------
	parseOptions(argc,argv);

	
	vector<mscsFunction3dregc> inf;
	mscsFunction3dregc outf;

//	inf.resize(_infiles.size());
//	if (_hdf5) {
//		for (unsigned long i = 0; i < _infiles.size(); i++) {
//			inf[i].loadHDF5(_infiles[i],_hdf5dset,0);
//		}		
//	}

	if (_operation=="add" or _operation=="sub" or _operation=="mul" or _operation=="divide") {
		inf.resize(1);
		for (unsigned long i = 0; i < _infiles.size(); i++) {
			msgs.say("loading file: "+_infiles[i],High);
			if (cpeds_fileExists(_infiles[i])) {
				if (_hdf5) {					inf[0].loadHDF5(_infiles[i],_dsets[i % _dsets.size()],0);				}
				else {
					inf[0].loadMatrix(_infiles[i]);
				}
				if (i==0) outf=inf[0];
				else {
					if (_operation=="add") {
						msgs.say("adding",Medium);
						outf+=inf[0];
					}
					if (_operation=="sub") {
						msgs.say("subtracting",Medium);
						outf-=inf[0];
					}
					if (_operation=="mul") {
						msgs.say("multiplying",Medium);
						outf*=inf[0];
					}
					if (_operation=="divide") {
						msgs.say("dividing",Medium);
						outf.divide(inf[0]);
					}
				}
			}
			else {
				msgs.say("file does not exist, skipping ("+_infiles[i]+")",High);
			}
		}
	}

	if (_operation=="mean") {		outf=average_matrices(msgs); }

	if (_operation=="dumpNonZero") {		outf=dumpNonZero(msgs); }
	if (_operation=="replaceBelow") {		outf=replaceBelow(msgs); }


	if (_operation=="mkTcmbtsz") { 		outf=make_TcmbtszMap(msgs); 	exit(0); }

	
	if (_hdf5)
		outf.saveHDF5(_outFile,_dsets[_dsets.size()-1]);
	else {
		msgs.say("In non hdf5 mode we assume a 2D matrix and only save the 0'th slice",Top);
		outf.saveSlice(2,0,_outFile,0);		
	}
	
	

	return 0;
}


void parseOptions(int argc, char** argv) {
	long i;
	_programVersionString=getProgramVersionString();
	try {
		
		//     CmdLine cmd("RT4\n\n. ",' ', "0.1" );
		CmdLine cmd("matrix_oper\nTorun Centre for Astronomy, UMK, POLAND.\n "
				"This program ...\n\n. "
				"\n"
				"example usage: matrix_oper"
				"",' ', _programVersionString.c_str() );
		//     CmdLine cmd("RT4 azimutal scan simulator for OCRA-p\n", 'This program simulates the RT4 scan motion, sky rotation during \n\n. ', 0.1 );
		
		// 
		// Define arguments
		//
		UnlabeledMultiArg<string> infiles("infiles", "input file names (<100)","", "string");     	cmd.add( infiles );
		ValueArg<string> outfile("o","outfile","name of output file",false,"","string"); cmd.add(outfile);

//		SwitchArg abs("","abs", "load absolute values of the maps (default: false)", false);	cmd.add(abs);
		ValueArg<string> comm("c","comm","operation to perform [sum, mkTcmbtsz ]\n"
				"\n"
				"add/sub/mul/divide - add/sub/mul/divide all input files and save to output file\n"
				"mean - calculate mean of all input files and save to output file. They must have the same dimensions. \n"
				"mkTcmbtsz - combine cmb map from the first input file [K] with the TSZE compton-y map taken from the second file."
				"The second file should be the hdf5 file containing comptonY dataset.",true,"","string"); cmd.add(comm);
//		ValueArg<double> val("","val","numerical value for various commands",false,1e300,"double"); cmd.add(val);
		SwitchArg hdf5("","hdf5", "intput files are hdf5 files (default: false)", false);	cmd.add(hdf5);
		ValueArg<string> freq("f", "freq", "comma separated list of frequencies for which to calculate cmb maps [GHz] (default: 15,22,30)", false,"15,22,30","string"); cmd.add( freq );
//		SwitchArg resizeCMB("","resizeCMB", "This options works only with operation -c mkTcmbtsz. This interpolates CMB map if it was provided in lower resolution than the comtonY map. (default: false)", false);	cmd.add(resizeCMB);

		ValueArg<double> val("", "val", "value for numerical manipulations (default: 0)", false,0,"double");	cmd.add( val );
		ValueArg<double> val2("", "val2", "value2 for numerical manipulations (default: 0)", false,0,"double");	cmd.add( val2 );
		
		ValueArg<string> hdf5dset("","hdf5dset","name of the hdf5 dset to be used. or comma separated list of datasets for each input file and one additional for the output file (default: '') ",false,"","string"); cmd.add(hdf5dset);
		ValueArg<string> hdf5_SZcomponent_output_dset_prefix("","hdf5_SZ_output_dset_prefix","name of the hdf5 dset to be used "
				"in all maps containing the converted compton Y-parameter component (e.g. TSZ or SZ maps)."
				"When processing the compton Y-parameter map containing converted kSZ component, you might want to "
				"change this name from TSZ to SZ to indicate that the output datasets contain both TSZ and KSZ contributions "
				"(default: 'TSZ') ",false,"","string"); cmd.add(hdf5_SZcomponent_output_dset_prefix);

		
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
		//     SwitchArg mink("M","Mink", "plot also the minkowski functionals for the map", false);	cmd.add( mink );
		
		cmd.parse(argc,argv);
		
		//
		// Set variables
		//
		_infiles = infiles.getValue();  if (_infiles.size() == 0) { printf("ERROR: not enough input files specified\n"); exit(0);}

//		_Ndec = Ndec.getValue();
//		_nLissajous=nLissajous.getValue();
//		_trackHoriz = horizTrack.getValue();
		
		_outFile = outfile.getValue();
		_operation=comm.getValue();
		_hdf5=hdf5.getValue();
		_hdf5dset=hdf5dset.getValue();
		_dsets=Mscs_stringToList(_hdf5dset);
		
		_freq=freq.getValue();
//		_resizeCMB=resizeCMB.getValue();
		_val=val.getValue();
		_val2=val2.getValue();
		_hdf5_SZcomponent_output_dset_prefix=hdf5_SZcomponent_output_dset_prefix.getValue();
		
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
mscsFunction3dregc make_TcmbtszMap(cpedsMsgs& msgs) {

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
	string frequencySuffix;

	
	if (_dsets[0]=="comptonY") { msgs.criticalError("the first file should be the CMB file, not compton-y. This will fail.",Top); }
	
	vector<mscsFunction3dregc> inf;
	vector<mscsFunction3dregc> infPtSrc;
	mscsFunction3dregc outf,DLTSZ;
	double FOVx,FOVy,x,y, cmb_interpolated,gnu;
	if (_infiles.size()<2 or _infiles.size()>3) 
		msgs.criticalError("There should be exactly two input files for this operation: \n"
				"1) cmb map \n"
				"2) compton-Y map.\n\n"
				"or\n\n"
				"exactly three if also point sources should be included. \n"
				"3) point source map file with flux densities expressed in Jy and contain datasets for each considered frequency"
				"named fluxDensity-%%.2fGHz",Top);

	inf.resize(2);
	if (_infiles.size()==3) {
		infPtSrc.resize(freq.size());
	}

	if (_hdf5) {
		for (unsigned long i = 0; i < _infiles.size(); i++) {
			if (i==0) {
				msgs.say("loading CMB map from file (assuming units [K]): "+_infiles[i],High);
				if (cpeds_fileExists(_infiles[i])) { inf[i].loadHDF5(_infiles[i],_dsets[i],0);	}
				else {	msgs.criticalError("file does not exist, skipping ("+_infiles[i]+")",Top);	}				
			}
			if (i==1) {
				msgs.say("loading compton-Y map from file: "+_infiles[i],High);
				if (cpeds_fileExists(_infiles[i])) { inf[i].loadHDF5(_infiles[i],_dsets[i],0);	}
				else {	msgs.criticalError("file does not exist, skipping ("+_infiles[i]+")",Top);	}				
			}
			if (i==2) {
				msgs.say("loading point source maps from file: "+_infiles[i],High);
				if (cpeds_fileExists(_infiles[i])) {
					string dset;
					for (unsigned long j = 0; j < freq.size(); j++) {
						dset="fluxDensity-"+msgs.toStrf(freq[j],2)+"GHz";
						msgs.say("loading map for frequency %lf GHz ",freq[j],Medium);
						msgs.say("using dataset: "+dset,Medium);
						infPtSrc[j].loadHDF5(_infiles[i],dset,0);
						_dsetsPtSrc.push_back(dset);
					}
				}
				else {	msgs.criticalError("file does not exist, skipping ("+_infiles[i]+")",Top);	}				
			}
		}
		
		// resize the CMB map region to match with FOV size using the internal key field_size
		int errCode;
		FOVx=inf[0].getHDF5_scalarDoubleAttribute(_infiles[0],"DTcmb","field_sizeX",&errCode);
		FOVy=inf[0].getHDF5_scalarDoubleAttribute(_infiles[0],"DTcmb","field_sizeY",&errCode);
		inf[0].setSizeRange(inf[0].Nx(),inf[0].Ny(),inf[0].Nz(),inf[1].getMinX(),inf[1].getMinY(),inf[1].getMinZ(),inf[1].getMaxX(),inf[1].getMaxY(),inf[1].getMaxZ());
		msgs.say("Resetting CMB map FOV range to: %lf x %lf",inf[0].lengthX(),inf[0].lengthY(),Medium);
//		if (inf[0].Nx()!=inf[1].Nx() or inf[0].Ny()!=inf[1].Ny()) _resizeCMB=true; else _resizeCMB=false;
//		inf[0].setVerbosityLevel(High);
//		inf[0].printInfo();
//		inf[1].setVerbosityLevel(High);
//		inf[1].printInfo();
		
		//
		// interpolate CMB map at the locations of the TSZ map grid cells and add to TSZ map
		//
		
		for (unsigned long i = 0; i < freq.size(); i++) {
			frequency=freq[i];
			msgs.say("PROCESSING FREQUENCY: %lf GHz",frequency,High);

			gnu=cpeds_TSZEgnu_factor(frequency*1e9,T_CBR0);
			msgs.say("gnu factor: %lf",gnu,Medium);
			outf=inf[1];
			outf*=double(gnu*T_CBR0);
			DLTSZ=outf; // store DTTSZ map

			long xi;
			if (inf[0].Nx()!=inf[1].Nx() or inf[0].Ny()!=inf[1].Ny()) {
				msgs.say("Interpolating input cmb map",Medium);
				mscsFunction3dregc cmb=outf;
				cmb.setVerbosityLevel(High);
				cmb.printInfo();
#pragma omp parallel for default(shared) private (xi,x,y)
				for (xi = 0; xi < cmb.Nx(); xi++) {
					x=cmb.getX(xi);
//					if (xi % 100==0) msgs.say("done: %lf \%",double(xi)/outf.Nx()*100,Low);
					for (long yj = 0; yj < cmb.Ny(); yj++) {
						y=cmb.getY(yj);
						cmb.fRe(xi,yj,0)=inf[0].fxy(x,y,0,0);
					}
				}
				inf[0]=cmb;
			}

			// add CMB
			outf+=inf[0]; 
			
			// save CMB+TSZ
			frequencySuffix="-"+msgs.toStrf(frequency,1)+"GHz";
			_hdf5dset="DTcmb"+_hdf5_SZcomponent_output_dset_prefix+frequencySuffix;
			outf.saveHDF5(_outFile,_hdf5dset);
			outf.setHDF5_scalarDoubleAttribute(_outFile,_hdf5dset,"gnu",gnu);
			outf.setHDF5_scalarDoubleAttribute(_outFile,_hdf5dset,"freq",frequency,"[GHz]");
			outf.setHDF5_scalarDoubleAttribute(_outFile,_hdf5dset,"Tcmb0",T_CBR0,"[K]");
			
			
			//
			// calculate spectral radiance map
			//
			double deg2toSr=PI180*PI180;
			double pixSize=FOVx/outf.Nx()*FOVy/outf.Ny()*deg2toSr; // sr
			double Jyosr2K=cpeds_radiance_to_temp_conversionFactor(frequency*1e9,T_CBR0);
			outf/=Jyosr2K;
			_hdf5dset="DLcmb"+_hdf5_SZcomponent_output_dset_prefix+frequencySuffix;
			outf.saveHDF5(_outFile,_hdf5dset);
			outf.setHDF5_scalarDoubleAttribute(_outFile,_hdf5dset,"freq",frequency,"[GHz]");
			outf.setHDF5_scalarDoubleAttribute(_outFile,_hdf5dset,"Tcmb0",T_CBR0,"[K]");
			msgs.say("pixel size [sr]: %lE",pixSize,Medium);
			msgs.say("T2Jyosr [Jy/sr/K]: %lE",1.0/Jyosr2K,Medium);
			
			// save TSZ radiance alone
			_hdf5dset="DL"+_hdf5_SZcomponent_output_dset_prefix+frequencySuffix;
			DLTSZ/=Jyosr2K;
			DLTSZ.saveHDF5(_outFile,_hdf5dset);
			DLTSZ.setHDF5_scalarDoubleAttribute(_outFile,_hdf5dset,"gnu",gnu);
			DLTSZ.setHDF5_scalarDoubleAttribute(_outFile,_hdf5dset,"freq",frequency,"[GHz]");
			DLTSZ.setHDF5_scalarDoubleAttribute(_outFile,_hdf5dset,"Tcmb0",T_CBR0,"[K]");
			DLTSZ.setHDF5_scalarDoubleAttribute(_outFile,_hdf5dset,"Jyosr2K",Jyosr2K,"[K/Jy/sr]");
			DLTSZ.setHDF5_scalarDoubleAttribute(_outFile,_hdf5dset,"pixelSize",pixSize,"[sr]");


			//
			// include point sources if requested
			//
			if (_infiles.size()==3) {
				if (infPtSrc[i].Nx()!=outf.Nx()) { msgs.say("The Nx():%li of the point source map doesn't match the Nx(): %li of the compton-y map",infPtSrc[i].Nx(),outf.Nx(),Top); exit(-1); }
				if (infPtSrc[i].Ny()!=outf.Ny()) { msgs.say("The Ny():%li of the point source map doesn't match the Ny(): %li of the compton-y map",infPtSrc[i].Nx(),outf.Nx(),Top); exit(-1); }
				infPtSrc[i]/=pixSize; // convert to spectral radiance from Jy
				outf+=infPtSrc[i]; // combine CMB+TSZ and ptSrc maps
				DLTSZ+=infPtSrc[i]; // combine TSZ and ptSrc maps

				_hdf5dset="DLcmb"+_hdf5_SZcomponent_output_dset_prefix+"ptSrc"+frequencySuffix;
				outf.saveHDF5(_outFile,_hdf5dset);
				outf.setHDF5_scalarDoubleAttribute(_outFile,_hdf5dset,"freq",frequency,"[GHz]");
				outf.setHDF5_scalarDoubleAttribute(_outFile,_hdf5dset,"Tcmb0",T_CBR0,"[K]");
				outf.setHDF5_scalarDoubleAttribute(_outFile,_hdf5dset,"pixelSize",pixSize,"[sr]");
				outf.setHDF5_scalarDoubleAttribute(_outFile,_hdf5dset,"Jyosr2K",Jyosr2K,"[K/Jy/sr]");

				_hdf5dset="DL"+_hdf5_SZcomponent_output_dset_prefix+"ptSrc"+frequencySuffix;
				DLTSZ.saveHDF5(_outFile,_hdf5dset);
				DLTSZ.setHDF5_scalarDoubleAttribute(_outFile,_hdf5dset,"freq",frequency,"[GHz]");
				DLTSZ.setHDF5_scalarDoubleAttribute(_outFile,_hdf5dset,"Tcmb0",T_CBR0,"[K]");
				DLTSZ.setHDF5_scalarDoubleAttribute(_outFile,_hdf5dset,"pixelSize",pixSize,"[sr]");
				DLTSZ.setHDF5_scalarDoubleAttribute(_outFile,_hdf5dset,"Jyosr2K",Jyosr2K,"[K/Jy/sr]");
			}
		}
	}
	else {			msgs.criticalError("this operation only works with hdf5 files",Top);		}
	
	// set the name for the output dataset
	return outf;
}
/***************************************************************************************/
mscsFunction3dregc average_matrices(cpedsMsgs& msgs) {
	vector<mscsFunction3dregc> inf;
	mscsFunction3dregc res;
	long count=0;
	inf.resize(1);
	for (unsigned long i = 0; i < _infiles.size(); i++) {
		msgs.say("loading file: "+_infiles[i],High);
		if (cpeds_fileExists(_infiles[i])) {
			if (_hdf5) {					inf[0].loadHDF5(_infiles[i],_dsets[i % _dsets.size()],0);				}
			else {
				inf[0].loadMatrix(_infiles[i]);
			}
			msgs.say("averaging",Medium);
			if (i==0) res=inf[0];
			else
				res+=inf[0];					
			count++;
		}
		else {
			msgs.say("file does not exist, skipping ("+_infiles[i]+")",High);
		}
	}
	
	res/=double(count);

	return res;
}
/***************************************************************************************/
mscsFunction3dregc replaceBelow(cpedsMsgs& msgs) {
	vector<mscsFunction3dregc> inf;
	mscsFunction3dregc res;
	inf.resize(1);
	for (unsigned long i = 0; i < _infiles.size(); i++) {
		msgs.say("loading file: "+_infiles[i],High);
		if (cpeds_fileExists(_infiles[i])) {
			if (_hdf5) {					inf[0].loadHDF5(_infiles[i],_dsets[i % _dsets.size()],0);				}
			else {
				inf[0].loadMatrix(_infiles[i]);
				_dsets.push_back("");
			}
			msgs.say("replacing values below val to val2",Medium);
			for (unsigned long i = 0; i < inf[0].size(); i++) {
				if (inf[0].fRe(i)<_val) inf[0].fRe(i)=_val2;
			}
			
		}
		else {
			msgs.say("file does not exist, skipping ("+_infiles[i]+")",High);
		}
	}


	return inf[0];	
}
/***************************************************************************************/
mscsFunction3dregc dumpNonZero(cpedsMsgs& msgs) {
	vector<mscsFunction3dregc> inf;
	mscsFunction3dregc res;
	cpedsList<double> l;
	inf.resize(1);
	for (unsigned long i = 0; i < _infiles.size(); i++) {
		msgs.say("loading file: "+_infiles[i],High);
		if (cpeds_fileExists(_infiles[i])) {
			if (_hdf5) {					inf[0].loadHDF5(_infiles[i],_dsets[i % _dsets.size()],0);				}
			else {
				inf[0].loadMatrix(_infiles[i]);
				_dsets.push_back("");
			}
			msgs.say("dumping non-zero values",Medium);
			l.clear();
			for (unsigned long i = 0; i < inf[0].size(); i++) {
				if (inf[0].fRe(i)!=0) l.append(inf[0].fRe(i));
			}
			l.save(_infiles[i]+"_"+_dsets[i % _dsets.size()]+".nonZeroDump",false,"double",15,false);
			
		}
		else {
			msgs.say("file does not exist, skipping ("+_infiles[i]+")",High);
		}
	}
	

	return res;
}
