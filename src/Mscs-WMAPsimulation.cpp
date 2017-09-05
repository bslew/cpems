#include "Mscs-WMAPsimulation.h"

//************************************************************************
mscsWMAPsimulation::~mscsWMAPsimulation() {
	for (int i = 0; i < sim.count(); ++i) {
		delete sim.value(i);
	}
	sim.clear();
//	printf("jestem w mscsWMAPsimulation destruktorze\n"); exit(0);
}

//************************************************************************
mscsWMAPsimulation::mscsWMAPsimulation(DAnames DA, WMAPversion year, const long ns, const long& lmax, const mscsAngularPowerSpectrum& C, string pixtfType, cpedsRNG *rns, bool precalib, double precalibrate, bool calib, double calibrate, mscsMap* mask) {
//mscsWMAPsimulation::mscsWMAPsimulation(DAnames DA, WMAPversion year, cpedsRNG *rns=NULL) {

	_DA=DA;
	_dr=year;
	_ns=ns;
	_lmax=lmax;
	_cl=C;
	_pixtfType=pixtfType;
	rng=rns;
	_precalib=precalib;
	_precalibrate=precalibrate;
	_calib=calib;
	_calibrate=calibrate;
	_mask=mask;

	_simNo=-1;


//	makeWMAPsimulation();

}

void mscsWMAPsimulation::makeWMAPsimulation() {

	bool doK1,doK2,doQ1, doQ2, doQ, doV1, doV2, doV, doW1, doW2, doW3, doW4, doW, doQVW, doVW, doSignalMap;
//	bool doDA[15];
	doK1=doK2=doQ1=doQ2=doQ=doV1=doV2=doV=doW1=doW2=doW3=doW4=doW=doQVW=doVW=doSignalMap=false;
//	for (int i = 0; i < 15; ++i) {		doDA=false;	}

	switch (_DA) {
		case K1:   doK1=true;		break;
		case K2:   doK2=true;		break;
		case Q1:   doQ1=true;		break;
		case Q2:   doQ2=true;		break;
		case V1:   doV1=true;		break;
		case V2:   doV2=true;		break;
		case W1:   doW1=true;		break;
		case W2:   doW2=true;		break;
		case W3:   doW3=true;		break;
		case W4:   doW4=true;		break;

		case Qinc: doQ1=true; doQ2=true;		break;
		case Vinc: doV1=true; doV2=true;		break;
		case Winc: doW1=true; doW2=true; doW3=true; doW4=true;		break;

		case QVWinc: doQ1=true; doQ2=true; doV1=true; doV2=true; doW1=true; doW2=true; doW3=true; doW4=true;		break;
		case VWinc:                        doV1=true; doV2=true; doW1=true; doW2=true; doW3=true; doW4=true;		break;

		case Q_V_W_inc: doQ1=true; doQ2=true; doV1=true; doV2=true; doW1=true; doW2=true; doW3=true; doW4=true;		break;
		case V_W_inc:                         doV1=true; doV2=true; doW1=true; doW2=true; doW3=true; doW4=true;		break;
		case V_W_VW_inc:                      doV1=true; doV2=true; doW1=true; doW2=true; doW3=true; doW4=true;		break;
//		case unitBeamTf:   doSignalMap=true;		break;
		default:
			break;
	}

//	msgs->say("doV1 doV2 is:\n"+msgs->toStr(doV1)+msgs->toStr(doV2),Top );

	mscsWMAPspecifications WMAPspec;

//	double sigma0;
//	mscsWindowFunction b;
	DAnames da;
	long N=0; // number of maps to be generated //TODO remove this variable it's not needed
	// define how many maps are to be generated,
	// the proper beam and noise amplitude for the requested DA
	
//	if (doSignalMap) {	N=1; da=unitBeamTf;
//		_beamtf.make_unit_kernel();
//		generateGaussianSimulation(_ns,_lmax,_cl,WMAPspec.get_beamtf(da,_dr),_pixtfType,WMAPspec.get_Nobs(da,_dr),WMAPspec.get_sigma0(da,_dr),rng,getPartialSimulationFileName(da)); setName("K1");
//		inc.useMap(this, WMAPspec.get_sigma0(da,_dr));		clear_map();	}

	
	if (doK1) {	N=1; da=K1;
		generateGaussianSimulation(_ns,_lmax,_cl,WMAPspec.get_beamtf(da,_dr),_pixtfType,WMAPspec.get_Nobs(da,_dr),WMAPspec.get_sigma0(da,_dr),rng,getPartialSimulationFileName(da)); setName("K1");
		inc.useMap(this, WMAPspec.get_sigma0(da,_dr));		clear_map();	}
	if (doK2) {	N=1; da=K2;
		generateGaussianSimulation(_ns,_lmax,_cl,WMAPspec.get_beamtf(da,_dr),_pixtfType,WMAPspec.get_Nobs(da,_dr),WMAPspec.get_sigma0(da,_dr),rng,getPartialSimulationFileName(da)); setName("K2");
		inc.useMap(this, WMAPspec.get_sigma0(da,_dr));		clear_map();	}
	if (doQ1) {	N=1; da=Q1;
		generateGaussianSimulation(_ns,_lmax,_cl,WMAPspec.get_beamtf(da,_dr),_pixtfType,WMAPspec.get_Nobs(da,_dr),WMAPspec.get_sigma0(da,_dr),rng,getPartialSimulationFileName(da)); setName("Q1");
		inc.useMap(this, WMAPspec.get_sigma0(da,_dr));		clear_map();	}
	if (doQ2) {	N=1; da=Q2;
		generateGaussianSimulation(_ns,_lmax,_cl,WMAPspec.get_beamtf(da,_dr),_pixtfType,WMAPspec.get_Nobs(da,_dr),WMAPspec.get_sigma0(da,_dr),rng,getPartialSimulationFileName(da)); setName("Q2");
		inc.useMap(this, WMAPspec.get_sigma0(da,_dr));		clear_map();	}
	if (doV1) {	N=1; da=V1;
		generateGaussianSimulation(_ns,_lmax,_cl,WMAPspec.get_beamtf(da,_dr),_pixtfType,WMAPspec.get_Nobs(da,_dr),WMAPspec.get_sigma0(da,_dr),rng,getPartialSimulationFileName(da)); setName("V1");
		inc.useMap(this, WMAPspec.get_sigma0(da,_dr));		clear_map();	}
	if (doV2) {	N=1; da=V2;
		generateGaussianSimulation(_ns,_lmax,_cl,WMAPspec.get_beamtf(da,_dr),_pixtfType,WMAPspec.get_Nobs(da,_dr),WMAPspec.get_sigma0(da,_dr),rng,getPartialSimulationFileName(da)); setName("V2");
		inc.useMap(this, WMAPspec.get_sigma0(da,_dr));		clear_map();	}
	if (doW1) {	N=1; da=W1;
		generateGaussianSimulation(_ns,_lmax,_cl,WMAPspec.get_beamtf(da,_dr),_pixtfType,WMAPspec.get_Nobs(da,_dr),WMAPspec.get_sigma0(da,_dr),rng,getPartialSimulationFileName(da)); setName("W1");
		inc.useMap(this, WMAPspec.get_sigma0(da,_dr));		clear_map();	}
	if (doW2) {	N=1; da=W2;
		generateGaussianSimulation(_ns,_lmax,_cl,WMAPspec.get_beamtf(da,_dr),_pixtfType,WMAPspec.get_Nobs(da,_dr),WMAPspec.get_sigma0(da,_dr),rng,getPartialSimulationFileName(da)); setName("W2");
		inc.useMap(this, WMAPspec.get_sigma0(da,_dr));		clear_map();	}
	if (doW3) {	N=1; da=W3;
		generateGaussianSimulation(_ns,_lmax,_cl,WMAPspec.get_beamtf(da,_dr),_pixtfType,WMAPspec.get_Nobs(da,_dr),WMAPspec.get_sigma0(da,_dr),rng,getPartialSimulationFileName(da)); setName("W3");
		inc.useMap(this, WMAPspec.get_sigma0(da,_dr));		clear_map();	}
	if (doW4) {	N=1; da=W4;
		generateGaussianSimulation(_ns,_lmax,_cl,WMAPspec.get_beamtf(da,_dr),_pixtfType,WMAPspec.get_Nobs(da,_dr),WMAPspec.get_sigma0(da,_dr),rng,getPartialSimulationFileName(da)); setName("W4");
		inc.useMap(this, WMAPspec.get_sigma0(da,_dr));		clear_map();	}

	if (_mask!=NULL and (_precalib or _calib) ) inc.import_map_data(*_mask,"m",1);

	switch (_DA) {
		case K1:   { inc.getMap(0).makekill_space_manager("kill","N"); sim.append(new mscsMap("K1",inc.getMap(0))); inc.removeMap(0); break; }
		case K2:   { inc.getMap(0).makekill_space_manager("kill","N"); sim.append(new mscsMap("K2",inc.getMap(0))); inc.removeMap(0); break; }
		case Q1:   { inc.getMap(0).makekill_space_manager("kill","N"); sim.append(new mscsMap("Q1",inc.getMap(0))); inc.removeMap(0); break; }
		case Q2:   { inc.getMap(0).makekill_space_manager("kill","N"); sim.append(new mscsMap("Q2",inc.getMap(0))); inc.removeMap(0); break; }
		case V1:   { inc.getMap(0).makekill_space_manager("kill","N"); sim.append(new mscsMap("V1",inc.getMap(0))); inc.removeMap(0); break; }
		case V2:   { inc.getMap(0).makekill_space_manager("kill","N"); sim.append(new mscsMap("V2",inc.getMap(0))); inc.removeMap(0); break; }
		case W1:   { inc.getMap(0).makekill_space_manager("kill","N"); sim.append(new mscsMap("W1",inc.getMap(0))); inc.removeMap(0); break; }
		case W2:   { inc.getMap(0).makekill_space_manager("kill","N"); sim.append(new mscsMap("W2",inc.getMap(0))); inc.removeMap(0); break; }
		case W3:   { inc.getMap(0).makekill_space_manager("kill","N"); sim.append(new mscsMap("W3",inc.getMap(0))); inc.removeMap(0); break; }
		case W4:   { inc.getMap(0).makekill_space_manager("kill","N"); sim.append(new mscsMap("W4",inc.getMap(0))); inc.removeMap(0); break; }

		case Kinc:   { sim.append(new mscsMap("Kinc",*inc.makeINC(-1,-1,"inv_noise",_precalib,_precalibrate,_calib, _calibrate))); sim.value(0)->makekill_space_manager("kill","N"); break; }
		case Qinc:   { sim.append(new mscsMap("Qinc",*inc.makeINC(-1,-1,"inv_noise",_precalib,_precalibrate,_calib, _calibrate))); sim.value(0)->makekill_space_manager("kill","N"); break; }
		case Vinc:   { sim.append(new mscsMap("Vinc",*inc.makeINC(-1,-1,"inv_noise",_precalib,_precalibrate,_calib, _calibrate))); sim.value(0)->makekill_space_manager("kill","N"); break; }
		case Winc:   { sim.append(new mscsMap("Winc",*inc.makeINC(-1,-1,"inv_noise",_precalib,_precalibrate,_calib, _calibrate))); sim.value(0)->makekill_space_manager("kill","N"); break; }
		case QVWinc: { sim.append(new mscsMap("QVWinc",*inc.makeINC(-1,-1,"inv_noise",_precalib,_precalibrate,_calib, _calibrate))); sim.value(0)->makekill_space_manager("kill","N"); break; }
		case VWinc:  { sim.append(new mscsMap("VWinc",*inc.makeINC(-1,-1,"inv_noise",_precalib,_precalibrate,_calib, _calibrate))); sim.value(0)->makekill_space_manager("kill","N"); break; }

		case Q_V_W_inc:     { sim.append(new mscsMap("Qinc",*inc.makeINC(0,1,"inv_noise",_precalib,_precalibrate,_calib, _calibrate))); sim.append(new mscsMap("Vinc",*inc.makeINC(2,3,"inv_noise",_precalib,_precalibrate,_calib, _calibrate))); sim.append(new mscsMap("Winc",*inc.makeINC(4,7,"inv_noise",_precalib,_precalibrate,_calib, _calibrate))); sim.value(0)->makekill_space_manager("kill","N"); sim.value(1)->makekill_space_manager("kill","N"); sim.value(2)->makekill_space_manager("kill","N");break; }
		case Q_V_W_QVW_inc: { sim.append(new mscsMap("Qinc",*inc.makeINC(0,1,"inv_noise",_precalib,_precalibrate,_calib, _calibrate))); sim.append(new mscsMap("Vinc",*inc.makeINC(2,3,"inv_noise",_precalib,_precalibrate,_calib, _calibrate))); sim.append(new mscsMap("Winc",*inc.makeINC(4,7,"inv_noise",_precalib,_precalibrate,_calib, _calibrate))); sim.append(new mscsMap("QVWinc",*inc.makeINC())); sim.value(0)->makekill_space_manager("kill","N"); sim.value(1)->makekill_space_manager("kill","N"); sim.value(2)->makekill_space_manager("kill","N"); sim.value(3)->makekill_space_manager("kill","N"); break; }
		case V_W_inc:       { sim.append(new mscsMap("Vinc",*inc.makeINC(0,1,"inv_noise",_precalib,_precalibrate,_calib, _calibrate))); sim.append(new mscsMap("Winc",*inc.makeINC(2,5,"inv_noise",_precalib,_precalibrate,_calib, _calibrate))); sim.value(0)->makekill_space_manager("kill","N"); sim.value(1)->makekill_space_manager("kill","N"); break; }
		case V_W_VW_inc:    { sim.append(new mscsMap("Vinc",*inc.makeINC(0,1,"inv_noise",_precalib,_precalibrate,_calib, _calibrate))); sim.append(new mscsMap("Winc",*inc.makeINC(2,5,"inv_noise",_precalib,_precalibrate,_calib, _calibrate))); sim.append(new mscsMap("VWinc",*inc.makeINC(0,5,"inv_noise",_precalib,_precalibrate,_calib, _calibrate))); sim.value(0)->makekill_space_manager("kill","N"); sim.value(1)->makekill_space_manager("kill","N"); sim.value(2)->makekill_space_manager("kill","N"); break; }
		default:

			break;
	}

	savePartialSimulations();
	inc.clearMaps();

}

//************************************************************************
// CAUTION!!! THESE PROCEDURES ARE FOR SPECIAL DEDICATED USE ONLY. BE SURE YOU KNOW WHAT YOU'RE DOING.
// wmap_data - defines what kind of INC simulation is to be generated, KINC, QINC,VINC,WINC,allQVWINC,
// save - 2 - save the files
//        1 - save the alms and C_l signal files only  (all stuff will be cleaned) ,
//        0 - don't save (the map will not be killed and alms and coords and T map will be present ); this enforces the quiet! option in the called routines
// in case of save=0 the return value will contain an array of INC sims of size  as requested in the wmap_data parameter and the ordering
// of the sims will be according to the rising frequency with the INC as the last one if requested

//mscsMap* mscsWMAPsimulation::generate_gaussian_maps_WMAP_QVW(long sim_idx,string how,package_dirs dirs,string wmap_data, string mask, double precalibrate, double calibrate,string what, string work_scheme ) {
//  mscsMap *s = new mscsMap[10];
//  mscsMap *sims=NULL;
//  long i,j,save,sim_num=0;
//  filenamestr tmp;
////  string danames[10] = {"k1","k2","q1","q2","v1","v2","w1","w2","w3","w4"};
//  bool Q=false,V=false,W=false,INC=false, K=false, KQVW=false;
//
//  if (Mscs_contains_str(work_scheme,"quiet!",6)) save=0 ;
//  else {
//    if (Mscs_contains_str(work_scheme,"quiet",5)) save=1;
//    else save=2;
//  }
//
//  if (what == "allQVWINC") { Q=true; V=true; W=true; INC=true;}
//  if (what == "KINC") { K=true; sim_num++;}
//  if (what == "QINC") { Q=true; sim_num++;}
//  if (what == "VINC") { V=true; sim_num++;}
//  if (what == "WINC") { W=true; sim_num++;}
//  if (what == "INC") { INC=true; sim_num++;}
//  if (what == "KQVWINC") { KQVW=true; INC=true;}
//  if (what == "KQVW") { KQVW=true; }
///*   if (what == "QVINC") { Q=true; V=true;} */
///*   if (what == "QWINC") { Q=true; W=true;} */
///*   if (what == "VWINC") { V=true; W=true;} */
//
//  if (save==0 && sim_num>0) sims = new mscsMap[sim_num];
//
//  // generating D/A maps
//  if (K || KQVW ) { // to be extended for ILC maps in future
///*     generate_gaussian_map("k1",sim_idx,0,how.c_str(),dirs,wmap_data,work_scheme); */
//    printf("!!!!!!!!! WARNING !!!!!!!!!!! THIS IS a special case where the work_scheme is replaced with \"\" so that k1 and k2 were saved but also other parameters (if any) will not be passed !!!!!!!!!!!!!!!\n");
//    printf("!!!!!!!!! WARNING !!!!!!!!!!! THIS IS a special case where the work_scheme is replaced with \"\" so that k1 and k2 were saved but also other parameters (if any) will not be passed !!!!!!!!!!!!!!!\n");
//    printf("!!!!!!!!! WARNING !!!!!!!!!!! THIS IS a special case where the work_scheme is replaced with \"\" so that k1 and k2 were saved but also other parameters (if any) will not be passed !!!!!!!!!!!!!!!\n");
//    printf("!!!!!!!!! WARNING !!!!!!!!!!! THIS IS a special case where the work_scheme is replaced with \"\" so that k1 and k2 were saved but also other parameters (if any) will not be passed !!!!!!!!!!!!!!!\n");
//
//    generate_gaussian_map("k1",sim_idx,0,how.c_str(),dirs,wmap_data,"");
//    s[0].set_map_nside(512);  s[0].import_map_data(*this,"T",2); s[0].import_map_data(*this,"N",2); makekill_space_manager("kill","N",1);
///*     generate_gaussian_map("k2",sim_idx,0,how.c_str(),dirs,wmap_data,work_scheme); */
//    printf("!!!!!!!!! WARNING !!!!!!!!!!! THIS IS a special case where the work_scheme is replaced with \"\" so that k1 and k2 were saved but also other parameters (if any) will not be passed !!!!!!!!!!!!!!!\n");
//    printf("!!!!!!!!! WARNING !!!!!!!!!!! THIS IS a special case where the work_scheme is replaced with \"\" so that k1 and k2 were saved but also other parameters (if any) will not be passed !!!!!!!!!!!!!!!\n");
//    printf("!!!!!!!!! WARNING !!!!!!!!!!! THIS IS a special case where the work_scheme is replaced with \"\" so that k1 and k2 were saved but also other parameters (if any) will not be passed !!!!!!!!!!!!!!!\n");
//    printf("!!!!!!!!! WARNING !!!!!!!!!!! THIS IS a special case where the work_scheme is replaced with \"\" so that k1 and k2 were saved but also other parameters (if any) will not be passed !!!!!!!!!!!!!!!\n");
//
//    generate_gaussian_map("k2",sim_idx,0,how.c_str(),dirs,wmap_data,"");
//    s[1].set_map_nside(512);  s[1].import_map_data(*this,"T",2); s[1].import_map_data(*this,"N",2); makekill_space_manager("kill","N",1);
//    if (save!=0) { sprintf(tmp,"if [ ! -d %s ]; then mkdir -p %s; fi\n",PROGRAM_DIRS[25],PROGRAM_DIRS[25]); system(tmp); } // create the INCK directory
//  }
//
//  if (Q  || KQVW || INC) {
//    generate_gaussian_map("q1",sim_idx,0,how.c_str(),dirs,wmap_data,work_scheme);
//    s[2].set_map_nside(512);  s[2].import_map_data(*this,"T",2); s[2].import_map_data(*this,"N",2); makekill_space_manager("kill","N",1);
//    generate_gaussian_map("q2",sim_idx,0,how.c_str(),dirs,wmap_data,work_scheme);
//    s[3].set_map_nside(512);  s[3].import_map_data(*this,"T",2); s[3].import_map_data(*this,"N",2); makekill_space_manager("kill","N",1);
//    if (save!=0) { sprintf(tmp,"if [ ! -d %s ]; then mkdir -p %s; fi\n",PROGRAM_DIRS[22],PROGRAM_DIRS[22]); system(tmp); }// create the INCQ directory
//  }
//
//  if (V  || KQVW || INC) {
//    generate_gaussian_map("v1",sim_idx,0,how.c_str(),dirs,wmap_data,work_scheme); // TEMPORARY CHANGE HERE: normally should be quiet
//    s[4].set_map_nside(512);  s[4].import_map_data(*this,"T",2); s[4].import_map_data(*this,"N",2); makekill_space_manager("kill","N",1);
//    generate_gaussian_map("v2",sim_idx,0,how.c_str(),dirs,wmap_data,work_scheme); // TEMPORARY CHANGE HERE: normally should be quiet
//    s[5].set_map_nside(512);  s[5].import_map_data(*this,"T",2); s[5].import_map_data(*this,"N",2); makekill_space_manager("kill","N",1);
//    if (save!=0) { sprintf(tmp,"if [ ! -d %s ]; then mkdir -p %s; fi\n",PROGRAM_DIRS[23],PROGRAM_DIRS[23]); system(tmp); } // create the INCV directory
//  }
//
//  if (W  || KQVW || INC) {
//    generate_gaussian_map("w1",sim_idx,0,how.c_str(),dirs,wmap_data,work_scheme);
//    s[6].set_map_nside(512);  s[6].import_map_data(*this,"T",2); s[6].import_map_data(*this,"N",2); makekill_space_manager("kill","N",1);
//    generate_gaussian_map("w2",sim_idx,0,how.c_str(),dirs,wmap_data,work_scheme);
//    s[7].set_map_nside(512);  s[7].import_map_data(*this,"T",2); s[7].import_map_data(*this,"N",2); makekill_space_manager("kill","N",1);
//    generate_gaussian_map("w3",sim_idx,0,how.c_str(),dirs,wmap_data,work_scheme);
//    s[8].set_map_nside(512);  s[8].import_map_data(*this,"T",2); s[8].import_map_data(*this,"N",2); makekill_space_manager("kill","N",1);
//    generate_gaussian_map("w4",sim_idx,0,how.c_str(),dirs,wmap_data,work_scheme);
//    s[9].set_map_nside(512);  s[9].import_map_data(*this,"T",2); s[9].import_map_data(*this,"N",2); makekill_space_manager("kill","N",1);
//    if (save!=0) { sprintf(tmp,"if [ ! -d %s ]; then mkdir -p %s; fi\n",PROGRAM_DIRS[24],PROGRAM_DIRS[24]); system(tmp); } // create the INCW directory
//  }
//
//
//  // load Nobs data
///*   for (i=0;i<8;i++) { */
///*     if (wmap_data == "wmap1" ) { sprintf(tmp,"%swmap-%s",dirs[3],danames[i].c_str()); loadfitsT(tmp,12); } */
///*     if (wmap_data == "wmap3" ) { sprintf(tmp,"%swmap3sum-%s",dirs[3],danames[i].c_str()); loadfitsT(tmp,14); } */
///*     if (wmap_data == "wmap5" ) { sprintf(tmp,"%swmap5sum-%s",dirs[3],danames[i].c_str()); loadfitsT(tmp,12); } */
///*     s[i].import_map_data(*this,"N",2); */
///*   } */
//
//  // load mask
//  if (mask == "kp03" || mask == "kp2" || mask == "kq75" || mask == "kq85") {    sprintf(tmp,"%s%s",dirs[2],mask.c_str()); loadfitsT(tmp,2);   }
//
// // make shift <T> = precalibrate (= 0)  // performing in this object to avoid loading the same mask into 10 objects or so
//  // precalibration accounts here for the statistics done outside the loaded mask
//  if (K || KQVW)        for (j=0;j<=1;j++) { import_map_data(s[j],"T",2);  shift_mean_to(precalibrate,true);   s[j].import_map_data(*this,"T",2); }
//  if (Q || INC || KQVW) for (j=2;j<=3;j++) { import_map_data(s[j],"T",2);  shift_mean_to(precalibrate,true);   s[j].import_map_data(*this,"T",2); }
//  if (V || INC || KQVW) for (j=4;j<=5;j++) { import_map_data(s[j],"T",2);  shift_mean_to(precalibrate,true);   s[j].import_map_data(*this,"T",2); }
//  if (W || INC || KQVW) for (j=6;j<=9;j++) { import_map_data(s[j],"T",2);  shift_mean_to(precalibrate,true);   s[j].import_map_data(*this,"T",2); }
//
//  //making INC maps
//  sim_num=0;
//  if (K || KQVW)        { mk_INC_DA(0,1,s,"inv_noise"); shift_mean_to(calibrate,true); sprintf(tmp,"%s512-1024-sim%li_%s_INC-k0-k0-%s-sm_no-Fmet_1",dirs[25],sim_idx,wmap_data.c_str(),"nomask"); if (save>=1) savebinT(tmp,1); if (save==0) { sims[sim_num].import_map_data(*this,"T",3); sim_num++; }}
//  if (Q || KQVW || INC) { mk_INC_DA(2,3,s,"inv_noise"); shift_mean_to(calibrate,true); sprintf(tmp,"%s512-1024-sim%li_%s_INC-q0-q0-%s-sm_no-Fmet_1",dirs[22],sim_idx,wmap_data.c_str(),"nomask"); if (save>=1) savebinT(tmp,1); if (save==0) { sims[sim_num].import_map_data(*this,"T",3); sim_num++; }}
//  if (V || KQVW || INC) { mk_INC_DA(4,5,s,"inv_noise"); shift_mean_to(calibrate,true); sprintf(tmp,"%s512-1024-sim%li_%s_INC-v0-v0-%s-sm_no-Fmet_1",dirs[23],sim_idx,wmap_data.c_str(),"nomask"); if (save>=1) savebinT(tmp,1); if (save==0) { sims[sim_num].import_map_data(*this,"T",3); sim_num++; }}
//  if (W || KQVW || INC) { mk_INC_DA(6,9,s,"inv_noise"); shift_mean_to(calibrate,true); sprintf(tmp,"%s512-1024-sim%li_%s_INC-w0-w0-%s-sm_no-Fmet_1",dirs[24],sim_idx,wmap_data.c_str(),"nomask"); if (save>=1) savebinT(tmp,1); if (save==0) { sims[sim_num].import_map_data(*this,"T",3); sim_num++; }}
//  if (INC)              { mk_INC_DA(2,9,s,"inv_noise"); shift_mean_to(calibrate,true); sprintf(tmp,"%s512-1024-sim%li_%s_INC-0-0-%s-sm_no-Fmet_1",dirs[9],sim_idx,wmap_data.c_str(),"nomask");    if (save>=1) savebinT(tmp,1); if (save==0) { sims[sim_num].import_map_data(*this,"T",3); sim_num++; }}
//
//  // free space
//  delete [] s;
//
//  if (save>=1) {
//    kill_map();  makekill_space_manager("kill","F",1);
//  }
//
//  return sims;
//}

// ********************************************************************
const mscsMap* mscsWMAPsimulation::getSimulation(string name) {

	for (int i = 0; i < sim.count(); ++i) {
		if (sim.value(i)->getName()==name) return sim.value(i);
	}
	return NULL;
}
/***************************************************************************************/
void mscsWMAPsimulation::saveSimulation(string outDir, string pref, long num, string suff, bool saveSignalAlms) {
	mscsWMAPspecifications ws;
	string fname=pref+msgs->toStr(num)+suff;
	string subDir=outDir+"/"+ws.getVersion(_dr)+"/";
	string dir;
	string cmd;

	for (int i = 0; i < sim.count(); ++i) {
		dir=subDir+sim.value(i)->getName()+"/";
		cmd="mkdir -p "+dir; system(cmd.c_str());
		sim.value(i)->savebinT(dir+fname);
	}
	
	if (saveSignalAlms) {
		mscsWMAPspecifications ws;
		string fname=pref+msgs->toStr(num)+suff;
		string subDir=_outDir+"/"+ws.getVersion(_dr)+"/";
		string dir=subDir+"alms/";
		cmd="mkdir -p "+dir; system(cmd.c_str());
		
		mscsAlms alms("alms");
		alms=getSignalAlms(); alms.savebinAlms(dir+fname+"-signal-alms-bin","RI");
		mscsAngularPowerSpectrum cl=alms.get_Cl(0,alms.lmax(),0);
		cl.save(dir+fname+"-signal-cl");		
	}

}
/***************************************************************************************/
void mscsWMAPsimulation::saveSimulation() {
	saveSimulation(_outDir,_outFilePref,_simNo,_outFileSuff);
}
/***************************************************************************************/
void mscsWMAPsimulation::setNames(string outDir, string pref, long num, string suff) {
	_outDir=outDir;
	_outFilePref=pref;
	_outFileSuff=suff;
	_simNo=num;
}
/***************************************************************************************/
void mscsWMAPsimulation::savePartialSimulations() {
	savePartialSimulations(_outDir,_outFilePref,_simNo,_outFileSuff);
}
/***************************************************************************************/
void mscsWMAPsimulation::savePartialSimulations(string outDir, string pref, long num, string suff) {
	if (savePartialSimulations_) {
		mscsWMAPspecifications ws;
		string fname=pref+msgs->toStr(num)+suff;
		string subDir=outDir+"/"+ws.getVersion(_dr)+"/";
		string dir;
		string cmd;

		for (int i = 0; i < inc.getMaps().count(); ++i) {
			dir=subDir+inc.getMaps().value(i)->getName()+"/";
			cmd="mkdir -p "+dir; system(cmd.c_str());
			inc.getMaps().value(i)->savebinT(dir+fname);
		}
	}
}

/***************************************************************************************/
string mscsWMAPsimulation::getPartialSimulationFileName(DAnames da) {
	
	mscsWMAPspecifications ws;
	if (ws.getDA(da)!="" && savePartialSimulations_) {
		string fname=_outFilePref+msgs->toStr(_simNo)+"-sig"+_outFileSuff;
		string subDir=_outDir+"/"+ws.getVersion(_dr)+"/";
		string dir;
		string cmd;
		
		
		dir=subDir+ws.getDA(da)+"/";
		cmd="mkdir -p "+dir; system(cmd.c_str());
		return dir+fname;
	}
	return "";
}
/***************************************************************************************/
//const vector<string> mscsWMAPsimulation::getSubDirName(DAnames DA, WMAPversion year) {
//	mscsWMAPspecifications ws;
//	string s=ws.getVersion(year)+"/";
//	string tmp=ws.getDA(DA);
//	vector<string> sd; // list of subdirectories
//
//	if (tmp.find(""))
//}
