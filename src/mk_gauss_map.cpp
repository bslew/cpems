/*!
 * \file mk_gauss_map.cpp
 *
 *  Project: Mscs
 *  Created on: Mar 26, 2013 11:19:40 AM
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
#include "cpeds-math.h"
#include "Mscs-map.h"
//#include "Mscs-global-defs.h"
#include "Mscs-angular_power_spectrum.h"
#include "Mscs-alms.h"
#include "Mscs-map-window_function.h"



#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif


string _outFile,_outDir,_mask_file, _Cl, _alms, _smf, _lowls, _noiseFile;
bool _from_Cl, _from_alms,_save_alms,_save_Cl,_masked, _nol01, _unit_alms,_check;
bool _white_noise,_extract_l,_tmp_block,_pixtf,_fastsims,_Folmax,_rrange;
long _lmax,_nside,_lmaxP,_wn_method,_extract_l_size,*_extract_ls,_loaduptol,_rr1,_rr2,_seed_offset;
double _smG, _ATfactor,_wns,_wnm,_Clcalib, _gaussianBeamFWHP;
string _uniform_noise, _uniform_noise_int;
bool _Clconv;

string _programVersionString;

void parseOptions(int argc, char** argv);
string getProgramVersionString();
void mkUniformNoise(string uniform_noise, mscsMap& sim, cpedsMsgs& msgs);
void mkUniformNoiseMap(string uniform_noise, mscsMap& sim, cpedsMsgs& msgs);
void mkUniformIntegerNoiseMap(string uniform_noise, mscsMap& sim, cpedsMsgs& msgs);
/*!
	\brief adapts resolution of noise map to the specified nside
	\details 
	@param map - noise map
	@param target_nside - requested target nside parameter

	If nside of the map > nside then dowongrade of the noise map is done
	and the noise is averaged between the pixels.
	If map.nside < nside, the map is prograded and new pixes are populated with
	gaussian random noise with the mean value of the initial pixel value
	and with standard deviation equal the standard deviation calculated from the pixels 
	that tile the pixel in the map resolution nside/2.

	\date May 16, 2016, 2:28:28 PM
	\author Bartosz Lew
*/
void fitNoiseToResolution(mscsMap& map, long target_nside);

int main(int argc, char **argv) {
	
	cpedsMsgs msgs(".mk_gauss_map",false,"",High);
	msgs.setSaveRunWriteMode('a');
	msgs.saveThisRun(argc,argv);
	string s;
	//----------------------------------------------------------------------------------------------------
	//PARSE THE COMMAND LINE ARGUMENTS
	//----------------------------------------------------------------------------------------------------
	parseOptions(argc,argv);
	
	//
	// prepare output directory
	//
	if (_outDir!="" and _outDir!="." and _outDir!="./") {
		mkdir(_outDir.c_str(),0755);
		chdir(_outDir.c_str());		
		_Cl=string("../"+_Cl);
	}
	
	
	mscsMap sim("simulation",_nside);
	//	mscsMap mask("mask",_nside);
	mscsAngularPowerSpectrum Cl,Clgen;
	mscsAlms a;
	mscsWindowFunction bl;
	
	//	  if(_pixtf) pixtf = sim.nside; else pixtf = 1;
	
	//	  if (_masked) {   // load mask if required
	//	    sprintf(file2,"%s%s",program_dir,_mask_file.c_str()); fname=file2; // fname=_mask_file;
	//	    if (fname.find("-bin",fname.size()-4) != string::npos) { fname.erase(fname.size()-7,7); mask.loadbinT(fname.c_str(),4); } else  //mask.change_map_resolution(sim.nside); 
	//	      //if (fname.find(".fits",fname.size()-5) != string::npos) mask.loadfitsT(file2,2);  //mask.change_map_resolution(sim.nside); 
	//	      mask.loadfitsT(file2,2); // assume that it's fits if there is no -bin extension
	//	    sim.import_map_data(mask,"m",2);
	//	    sim.check_mask();    //sim.mask_map_merge();
	//	  }
	
	//	  if (_white_noise) {
	//	    wnmap = cpeds_random_gauss_numbers(_wnm,_wns,sim.pix_num,_wn_method,cpeds_seed_offset,_fastsims); sim.makekill_space_manager("make","T",1);
	//	    for (i=0;i<sim.pix_num;i++) sim.set_T(i,wnmap[i]); sim.savebinT(_outFile.c_str(),1);
	//	    return 0;
	//	  }
	
	if (_uniform_noise!="") {
		mkUniformNoiseMap(_uniform_noise, sim,msgs);
		exit(0);
	}
	if (_uniform_noise_int!="") {
		mkUniformIntegerNoiseMap(_uniform_noise_int, sim,msgs);
		exit(0);
	}
	
	// check out about the power spectrum of the map
	
	// load C_l file or alms file
	if (_from_Cl) {     
		Cl.load(_Cl,true); 
		msgs.say("loading Cl done",Low);
		if (_Clconv) Cl.divide_llpotwoPI();
		//		if (_Clcalib != "nocalib") {
		//			if (_Clcalib=="Tcmb") { ClunitConversionFactor=double(T_CBR0*T_CBR0); }
		//			if (_Clcalib=="TcmbWMAP") {
		//				if (_yr=="WMAP1yr" or _yr=="WMAP3yrs" or _yr=="WMAP5yrs" or _yr=="WMAP7yrs") {
		//					ClCalibrationFactor=wmapSpecs.get_Cl_calibration(wmapVersion);
		//					ClunitConversionFactor=double(T_CBR0*T_CBR0);
		//				}
		//			}
		//		}
		Cl*=_Clcalib;
	}
	

	if (_from_alms) {
		msgs.say("loading alms",High);
		a.setVerbosityLevel(High);
		a.loadbinAlms(_outFile+".alms","RI");
		if (a.lmax()!=_lmax) { msgs.say("loaded alms have lmax different from requested. Cannot use it.",High); }
	//	_save_alms=false;
//		printf("done lmax=%li %li\n",a.lmax(),_lmax);
	}
	//		if (_Folmax) { sim.loadbinF(_alms.c_str(),5,_loaduptol); _lmax=_loaduptol; } else sim.loadbinF(_alms.c_str(),1); 
	//		if (sim.lmax > _lmax) sim.zero_alms_multipole_range(sim.lmax+1,_lmax); 
	//		if (_nol01) { printf("mk_gauss_map> REMOVING MULTIPOLES l=0 and l=1 \n"); sim.zero_alms_multipole_range(0,1); }
	//		if (_rrange) { printf("mk_gauss_map> REMOVING MULTIPOLE RANGE\n"); sim.zero_alms_multipole_range(_rr1,_rr2); }
	//		if (_nolq.get_size() > 0) for (i=0; i<_nolq.get_size();i++) sim.zero_alms_multipole_range(_nolq(i),_nolq(i));
	//		sim.calculate_C_l(0,_lmax,1); 
	//	}
	//	else {  
	//		printf("mk_gauss_map> DOING FLAT POWER SPECTRUM\n"); 
	//		//if (_rr1!=-1 && _rr2!=0) sim.zero_alms_multipole_range(_rr1,_rr2);  
	//		sim.set_Cl_multipole_range(0, _lmax,1); 
	//	} // doing flat power spectrum
	
	// generate alms if needed
	//	a.makekill_space_manager("make","F",1);
	if (a.lmax()!=_lmax) { // TODO this is stupid: we should not regenerate alms unless the power spectrum changed, not lmax. This is a very easy change.
		msgs.say("Generating gaussian alms",High);
		a.clear();
		a.lmax(_lmax);
		a.generate_gaussian_alms(Cl,_lmax,0,false,false,-1,NULL);
		a.antysymmetrize();
		msgs.say("done",Low);
	}
	if (_save_alms) { 
		msgs.say("Saving alms to file",High);
		a.savebinAlms(_outFile+".alms","RI"); 
		msgs.say("done",Low);
	}
	if (_save_Cl) { 
		msgs.say("Saving generated Cl to file",High);
		Clgen=a.get_Cl(0,_lmax,0);  Clgen.save(_outFile+".cl"); 
		Cl.save("input.cl");
		msgs.say("done",Low);
	}
	msgs.say("Making unit kernel",High);
	bl.make_unit_kernel(_lmax);
	msgs.say("done",Low);
	msgs.say("Making SH synthesis",High);
	sim.SH_synthesis(a,_lmax,bl,bl);
	msgs.say("done",Low);

	msgs.say("Saving temperature map",High);
	sim.savebinT(_outFile+"-Tn-bin");
	msgs.say("done",Low);

	
	if (_gaussianBeamFWHP>0) {	
		mscsWindowFunction wf=bl;
		bl.make_gaussian_kernel_kspace(_gaussianBeamFWHP/60.0*PI180);	
		bl.save("beam_tf.txt");
		msgs.say("Making SH synthesis for beam-smoothed map",High);
		sim.SH_synthesis(a,_lmax,bl,wf);
		msgs.say("done",Low);

		msgs.say("Saving temperature map",High);
		_outFile+="-beam_"+msgs.toStrf(_gaussianBeamFWHP,2);
		sim.savebinT(_outFile+"-Tn-bin");
		msgs.say("done",Low);
	
	}
	if (_noiseFile!="") {
		msgs.say("Making beam-smoothed map with noise",High);
		mscsMap noise("noise");
		if (_noiseFile.substr(_noiseFile.size()-4,4)=="fits") {
			noise.loadPLANCK_temp_fits(_noiseFile,"\'FULL SKY MAP\'");
		}
		else noise.loadbinT(_noiseFile);
		msgs.say("calculating noise",High);
		fitNoiseToResolution(noise, sim.nside());
		
//		noise.savebinT("noise-Tn-bin");
//		exit(0);
//		sim.setVerbosityLevel(High);
//		sim.calculate_map_stats(1);
//		noise.calculate_map_stats(1);
//		exit(0);
		// add noise
		sim.T()+=noise.T();

		msgs.say("Saving temperature map",High);
		sim.savebinT(_outFile+"-noise-Tn-bin");
		msgs.say("done",Low);
		
	}
	
	//	  if ( _from_alms == false) {
	//	    if (_unit_alms) sim.generate_gaussian_alms("RI","unnormCl",_lowls.c_str(),0,1,_lmax,4,-1,_fastsims);
	//	    if (_from_Cl) sim.generate_gaussian_alms("RI","unnormCl",_lowls.c_str(),0,1,_lmax,5,-1,_fastsims);
	//	  }
	
	//	  if (_extract_l) {
	//	    lmap = new map_class("lmap",data_dir); 
	//	    for (i=0;i<_extract_l_size;i++) {
	//	      (*lmap).set_map_nside(sim.nside); (*lmap).set_alms_lmax(sim.lmax);  (*lmap).import_map_data(sim,"alms",1);
	//	      (*lmap).zero_alms_multipole_range(0,_extract_ls[i]-1);
	//	      (*lmap).zero_alms_multipole_range(_extract_ls[i]+1,sim.lmax);
	//	      (*lmap).set_alms_lmax(_extract_ls[i],"copy"); 
	//	      (*lmap).calculate_inverse_transformF(1,pixtf,_smf.c_str(),_smG,0,"","","");
	//	      sprintf(tmpch,"%s-l%li",_outFile.c_str(),_extract_ls[i]); (*lmap).savebinT(tmpch,1);
	//	      (*lmap).makekill_space_manager("kill","F",1);
	//	    }
	//	    delete lmap;
	//	    exit(0);
	//	  }
	
	//	  if (_masked) sim.mask_map_merge();
	//	  sim.savebinT(_outFile.c_str(),1);
	
	//	// do check if required
	//	if (_check) {
	//		/*     _lmax=2048;  _lmaxP=256; _nside=512; */
	//		/*     sim.change_map_resolution(_nside);    sim.savebinT("downgraded",1); */
	//		/*     sim.makekill_space_manager("kill", "F",1); sim.set_alms_lmax(_lmax); sim.makekill_space_manager("make", "F",1);   Mscs_initiate_global_variables(_nside,_lmaxP); */
	//		
	//		sim.set_alms_multipole_range(0,_lmax,"RI",0);
	//		//sim.calculate_transformF(8,1,_smf.c_str(),_smG,0,"",GLOBAL_Plmsfile_forward);
	//		sim.calculate_transformF(1,pixtf,_smf.c_str(),_smG,0,"","");
	//		fname="check_"+_outFile;    sim.savebinF(fname.c_str(),1); sim.savetxtF(fname.c_str(),1);   sim.calculate_C_l(0,_lmax,1);  sim.savetxtC_l(fname.c_str(),1); 
	//	}
	//	
	//	if (_extract_l) delete [] _extract_ls;
	
	
	
	
	
	return 0;
}


void parseOptions(int argc, char** argv) {
	long i;
	_programVersionString=getProgramVersionString();
	
	try {
		
		//     CmdLine cmd("RT4\n\n. ",' ', "0.1" );
		CmdLine cmd("mk_gauss_map\nTorun Centre for Astronomy, UMK, POLAND.\n "
				"This program ...\n\n. "
				"\n"
				"example usage: mk_gauss_map"
				"",' ', _programVersionString.c_str() );
		
		ValueArg<string> outFile("o","outFile","outfile prefix (default: cmb)",false,"cmb","string"); cmd.add(outFile);
		ValueArg<string> outDir("","outDir","output main directory (default: current dir)",false,"","string"); cmd.add(outDir);
		ValueArg<long> nside("n", "nside", "nside of the map (default: 512)", false,512,"long");	cmd.add( nside );
		ValueArg<long> lmax("l", "lmax", "lmax of the alms (default: 1024)", false,1024,"long");	cmd.add( lmax );
		
		
		ValueArg<string> Cl("","Cl", "file with Cls in CMBfast format (no header)", false,"","string");	cmd.add( Cl );
		SwitchArg from_alms("","from_alms", "use saved alms in file outFile.alms for this run. (default: false)", false);	cmd.add(from_alms);
//		ValueArg<string> alms("","alms", "file with almss in RI format (no header)", false,"","string");	cmd.add( alms );
		
		//		allowedStr.clear();
		//		allowedStr.push_back("nocalib");
		//		allowedStr.push_back("Tcmb");
		//		allowedStr.push_back("TcmbWMAP");
		//		allowedStr.push_back("other");
		//		allowedClCal = new ValuesConstraint<string>(allowedStr);
		//		ValueArg<string> Clcalib("","clcal","multiply theoretical cl by a calibration factor",false,"nocalib",allowedClCal); cmd.add(Clcalib);
		ValueArg<double> Clcalib("","cal","calibration factor. Cl will be multiplied by this value. (default: 1)",false,1,"double"); cmd.add(Clcalib);
		//		ValueArg<double> ClcalibVal("","cal","calibration factor. Cl will be multiplied by this value of clcal \"other\" was used",false,0,"double"); cmd.add(ClcalibVal);
		ValueArg<double> gaussianBeamFWHP("","gaussianBeamFWHP","Use this value to create beam transfer function other "
				"than unity for the map generation process. "
				"The beam transfer function is generated in SH space and convolved in SH space with map."
				"The beamwidth is given in arcmin. (default: -1 = not used)",false,-1,"double"); cmd.add(gaussianBeamFWHP);
		ValueArg<string> noiseFile("","noiseFile","Heapix nested format noise file."
				"Should be of the same ns as the generated cmb map. "
				"If the resolution of the noise map is lower than the generated cmb map then"
				"while noise component will be generated in the prograded sub-pixels,"
				"otherwise downgrade is applied and the noise is averaged. (default: '' )",false,"","string"); cmd.add(noiseFile);

//		ValueArg<string> tf("","tf","File name containing a transfer function."
//				"If provided an extra file will contain the generated cmb map convolved with this function in SH space. "
//				"If the resolution of the noise map is lower than the generated cmb map then"
//				"while noise component will be generated in the prograded sub-pixels,"
//				"otherwise downgrade is applied and the noise is averaged. (default: '' )",false,"","string"); cmd.add(noiseFile);

		SwitchArg Clconv("","clconv", "whether or not to divide the input Cl by l(l+1)/twoPI factor before generating alms, default: false)", false);	cmd.add(Clconv);
		
		
		
		//		ValueArg<string> alms("","alms", "file with alms in Mscs bin format ", false,"","string");	cmd.add( alms );
		SwitchArg save_alms("","save_alms", "whether or not to save the alms file (default: false)", false);	cmd.add(save_alms);
		SwitchArg save_Cl("","save_Cl", "whether or not to save the Cl file (default: false)", false);	cmd.add(save_Cl);
		SwitchArg nol01("","nol01", "whether or not erase 0,1 multipole (default: false)", false);	cmd.add(nol01);
		//		ValueArg<long> rr1("","rr1","remove range of multipoles: from multipole in SH space",false,-1,"long"); cmd.add(rr1);
		//		ValueArg<long> rr2("","rr2","remove range of multipoles: to multipole in SH space",false,-1,"long"); cmd.add(rr2);
		//		MultiArg<long> nol("r","r","remove the individual multipole number in SH space",false,"long"); cmd.add(nol);
		
		//		SwitchArg unit_alms("","unit_alms", "whether or not make unitary all alms (default: false)", false);	cmd.add(unit_alms);
		//		SwitchArg check("","check", "whether or not make fwd STH to check resluts (default: false)", false);	cmd.add(check);
//		SwitchArg white_noise("","white_noise", "makes gaussian real space white noise map with variace --wns (default: false)", false);	cmd.add(white_noise);
		ValueArg<string> uniform_noise("","uniform_noise", "makes uniform real space white noise map. "
				"You should specify"
				"min and max values in the map (default: '')", false, "", "string");	cmd.add(uniform_noise);
		ValueArg<string> uniform_noise_int("","uniform_noise_int", "makes uniform real space white noise map of integer "
				"values. You should specify"
				"min and max values in the map (default: '')", false, "","string");	cmd.add(uniform_noise_int);
		//		SwitchArg tmp_block("","tmp", "perform the special temporary block in the program and quit", false);	cmd.add(tmp_block);
		//		SwitchArg pixtf("","pixtf", "use pix.tf. for given resolution", false);	cmd.add(pixtf);
		//		SwitchArg fastsims("","fastsims", "use when the simsulation will be repeatedly done faster than a second to help avoind dwawing the same set of numbers from GRF", false);	cmd.add(fastsims);
		//		ValueArg<double> wns ("","wns","white noise variance (default: 1)",false,1,"double"); cmd.add(wns);
		//		ValueArg<double> wnm ("","wnm","white noise mean (default: 0)",false,0,"double"); cmd.add(wnm);
		//		ValueArg<long> wn_method ("","wn_method","white noise mean (default: 3 - cpeds numbering)",false,3,"long"); cmd.add(wn_method);
		//		ValueArg<long> seed_offset ("","seed_offset","offset of the seed for the GRNG for parallel runs",false,3,"long"); cmd.add(seed_offset);
		
		
		//
		// Parse the command line.
		//
		cmd.parse(argc,argv);
		
				
		//
		// Set variables
		//
		//		_mask_file = mask.getValue(); if (_mask_file.length() == 0) { _masked = false; } else {_masked = true; }
		_outFile = outFile.getValue(); 
		_outDir=outDir.getValue();
		_nside = nside.getValue(); 
		_lmax = lmax.getValue(); 
		//		_loaduptol = loaduptol.getValue(); 
		//		_lmaxP = lmaxP.getValue(); 
		//		_smG = smG.getValue(); if (_smG == 0) _smG=-1;
		//		_smf = smf.getValue();
		//		_ATfactor = ATfactor.getValue();
		
		_Cl = Cl.getValue(); if (_Cl.size() == 0) { _from_Cl = false; } else { _from_Cl = true; }
//		_alms = alms.getValue(); 
		_from_alms=from_alms.isSet();
		_Clconv = Clconv.getValue();
		_Clcalib = Clcalib.getValue();
		_gaussianBeamFWHP=gaussianBeamFWHP.getValue();
		_noiseFile=noiseFile.getValue();
		//		_alms = alms.getValue(); if (_alms.size() == 0) { _from_alms = false; } else { _from_alms = true; }
		_save_alms = save_alms.getValue();
		_save_Cl = save_Cl.getValue();
		_nol01 = nol01.getValue(); if (_nol01 == true) { _lowls = "nol0-1"; }  else { _lowls = ""; }
		//		_unit_alms = unit_alms.getValue();
		//		_check = check.getValue();
		
		//		_rr1 = rr1.getValue();
		//		_rr2 = rr2.getValue();
		//		if ((_rr1 == -1) || (_rr2 == -1)) { _rrange = false; } else { _rrange = true; }
		//		vector<long> noltab = nol.getValue();	
		//		if (noltab.size() > 0) 
		//		  for ( i = 0; (unsigned int)i < noltab.size(); i++ ) { _nolq.addq(noltab[i]); }
		
		
		//		_white_noise = white_noise.getValue();
		_uniform_noise = uniform_noise.getValue();
		_uniform_noise_int = uniform_noise_int.getValue();
		//		_wns = wns.getValue();
		//		_wnm = wnm.getValue();
		//		_wn_method = wn_method.getValue();
		
		//		_tmp_block = tmp_block.getValue();
		//		_pixtf = pixtf.getValue();
		//		_fastsims = fastsims.getValue();
		//		_seed_offset = seed_offset.getValue(); cpeds_seed_offset=_seed_offset;
		
		//		_Folmax = Folmax.getValue();
		//
		//		vector<long>	vec;
		//		vec = extract_l.getValue(); 
		//		if (vec.size() > 0) { 
		//		  _extract_l = true; _extract_l_size = (long)vec.size(); _extract_ls = new long[_extract_l_size]; 
		//		  for (i=0;i<(long)vec.size();i++) { _extract_ls[i] = vec[i]; }	}
		
		
		
		
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
void fitNoiseToResolution(mscsMap& map, long target_nside) {	
	if (map.nside()==target_nside) return;
	if (map.nside()>target_nside) {
//		map.setVerbosityLevel(High);
		map.change_map_resolution(target_nside);
	}
	else {
		if (map.nside()<=4) { printf("findNoiseToResolution:: nside too low\n"); exit(1); }
		mscsMap mmmap=map;
		mmmap.make_multi_maskHP(target_nside/4); // will calculate variance statistics from 16 pixels
		stat_info *stat=mmmap.calculate_map_variance_with_multi_mask();
		for (long i = 0; i < mmmap.pixNum(); i++) { mmmap.T(i)=stat[i].s;	}
		mmmap.change_map_resolution(target_nside); // distribute the same variance among pixels in the target resolution
		
		// generate white noise in the prograded map
		cpedsRNG rns("gaussian_circle");
		int output;
		map.change_map_resolution(target_nside);
		for (long i = 0; i < map.pixNum(); i++) {
			rns.setMeanVariance(0.0,mmmap.T(i));
			map.T(i)=rns.getRN();
		}
		
	}
	
	
}
/* ******************************************************************************************** */
void mkUniformNoise(string uniform_noise, mscsMap& sim, cpedsMsgs& msgs) {
	cpedsRNG rns("uniform");
	QString qs=uniform_noise.c_str();
	QStringList qsl=qs.split(',');
	if (qsl.size()!=2) {
		msgs.error("You should provide range",High);
		exit(-1);
	}
	rns.setMinMax(qsl[0].toDouble(),qsl[1].toDouble());
	
	sim.T()=rns.getRNs(sim.pixNum());
	
}

/* ******************************************************************************************** */
void mkUniformNoiseMap(string uniform_noise, mscsMap& sim, cpedsMsgs& msgs) {
	mkUniformNoise(uniform_noise,sim,msgs);
	msgs.say("Saving temperature map",High);
	sim.savebinT(_outFile+"-Tn-bin");
	msgs.say("done",Low);
}

/* ******************************************************************************************** */
void mkUniformIntegerNoiseMap(string uniform_noise, mscsMap& sim, cpedsMsgs& msgs) {
	mkUniformNoise(uniform_noise,sim,msgs);
	for (long i = 0; i < sim.pixNum(); i++) {
		sim.T()[i]=round(sim.T()[i]);
	}
	msgs.say("Saving temperature map",High);
	sim.savebinT(_outFile+"-Tn-bin");
	msgs.say("done",Low);
}

