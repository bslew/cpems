#include "Mscs-gaussian_simulation.h"

/***************************************************************************************/
mscsGaussianSimulation::mscsGaussianSimulation() : a("CMBsigmalAlms",0) { 	
	saveSignalMap_=false;
}

//************************************************************************
mscsGaussianSimulation::mscsGaussianSimulation(const long& ns, const long& lmax, const mscsAngularPowerSpectrum& C, const mscsWindowFunction& b, string pixtfType, mscsMap& Nobs, const double& sigma0, cpedsRNG *rns) : a("CMBsigmalAlms",0) {
	generateGaussianSimulation(ns,lmax,C,b,pixtfType,Nobs,sigma0,rns);
	saveSignalMap_=false;
}


//************************************************************************
void mscsGaussianSimulation::generateGaussianSimulation(const long& ns, const long& lmax, const mscsAngularPowerSpectrum& C, const mscsWindowFunction& b, string pixtfType, mscsMap Nobs, const double& sigma0, cpedsRNG *rns, string signalMap) {
	set_nside(ns);
	makeCMBsignalMap(lmax,C,b,pixtfType, rns);
	if (signalMap!="") { savebinT(signalMap); }
	if (Nobs.nside()!=ns) {
		Nobs.change_map_resolution(ns);
	}
	N()=Nobs.get_N();
	addWhiteGaussianNoise(0.0,sigma0,true,rns);
}

//************************************************************************
void mscsGaussianSimulation::makeCMBsignalMap(const long& lmax, const mscsAngularPowerSpectrum & C, const mscsWindowFunction & b, string pixtfType, cpedsRNG *rns) {
	if (a.lmax()!=lmax) { // TODO this is stupid: we should not regenerate alms unless the power spectrum changed, not lmax. This is a very easy change.
		a.clear();
		a.lmax(lmax);
		a.generate_gaussian_alms(C,lmax,0,false,false,-1,rns);
		a.antysymmetrize();
	}
	SH_synthesis(a,lmax,b,pixtfType);
//	savebinT("test-signal");
}

//************************************************************************
void mscsGaussianSimulation::makeCMBsignalMap(const long& lmax, const mscsWindowFunction& b, string pixtfType) {
	SH_synthesis(a,lmax,b,pixtfType);
}

//************************************************************************
void mscsGaussianSimulation::addWhiteGaussianNoise(double m, double s, bool useNobs, cpedsRNG *rns) {
	mscsMap map("tmpMap",nside());
	map.makekill_space_manager("make","T");
	bool rnsWasNull;
	cpedsRNG::distrType RNGtype;

	if (useNobs) {
		long n=map.pixNum();
		double sigma0=s;
		if (rns==NULL) { rns = new cpedsRNG("gaussian","double"); rnsWasNull=true; } 
		else { 
			rnsWasNull=false; 
			// store the current generator type;
			RNGtype=rns->getRNsType();
			rns->setRNsType("gaussian");
		}
		for (int i = 0; i < n; ++i) {
			s=sigma0/sqrt(N(i));
			rns->setMeanStd(m,s);
			map.T(i)=rns->getRN();
		}
		if (rnsWasNull) delete rns;
	}
	else {
		map.make_gaussian_map(m,s,1,rns);
	}
//	map.T().save("noise-Tn-bin",true);
	T()+=map.T();
	if (rnsWasNull) delete rns;
	else {
		rns->setRNsType(RNGtype);		
	}
}


//************************************************************************

////************************************************************************
//// generates a full gaussian map according to stored model power spectrum (needs to be loaded first - then it's replaced with the one generated)
//// for a requested DA (differential assembly) of WMAP radiometer, with according gaussian noice patterns, smoothed to the requested beam (by default
//// the DA beam size gaussian smoothing)
//
///* this is the abbreviated version of the procedure, if you want the have a full control on the parameters required for the generation just call the right  */
///* function straight away with a buchn of parameters. This is actually only to setup the parameters for map simulation */
//void mscsGaussianSimulation::generate_gaussian_map(strarg DA, long sim_num, int save_partial, strarg how,package_dirs dirs, string wmap_data,string work_scheme) {
//  long lmax_loc,nside_loc;
//  filenamestr Clfile;
//  double sigma0=-1;
//  double smooth_user=-1;
//  filenamestr WMAPNobsfile,WMAPNobsfile3,bl_file;
//  int gaussgenerator;
//  filenamestr sim_dir,alms_file;
//  filenamestr files_prefix, files_prefix1,files_prefix3;
//  filenamestr plms_file;
//  filenamestr show_process,tmp;
//  filenamestr DA1;
//  long DAi1=0;
//  int fourier_method;
//  long sigma0offset;
//
//
///*----------------------  definition of the initial resolution parameter   -----------------*/
//  nside_loc = 512;
//
///*----------------------  definition of the range of the alms to be generated   ------------*/
//  lmax_loc = 1024;
//  //lmax_loc = 80;
//  fourier_method = 1;
//  //lmax_loc = 2;
//
///*----------------------  definition of the underlying power spectra model -----------------*/
//   // theoretical models
///*   sprintf(Clfile,"%sbestWMAP",program_dir); */
//  //sprintf(Clfile,"%sbestWMAP",dirs[3]);
//
//  if (wmap_data == "wmap13" || wmap_data == "wmap3") {
//    sprintf(Clfile,"%sbestLCDMflat6par.cl",dirs[3]); // WMAP3 one of the  LCDM best fits power spectrum as in Spergel (http://lambda.gsfc.nasa.gov/product/map/dr2/pub_papers/threeyear/parameters/64897.web.pdf) p.388 tab 5 col 2 for wmap3 only ; use this for normal GRF sims
//  }
///*   if (wmap_data == "wmap5") */
///*     sprintf(Clfile,"%sbestLCDMflat6par-cmb-wmap5.cl",dirs[3]); // WMAP5 power spectrum LCDM best fit to CMB only max likelihood Cl; use this for normal GRF sims */
///*   if (wmap_data == "wmap5")  */
///*     sprintf(Clfile,"%swmap_lcdm_sz_lens_wmap5_cl_v3.dat",dirs[3]); // WMAP5 power spectrum LCDM best fit ; use this for normal GRF sims */
//  if (wmap_data == "wmap5") {
//    sprintf(Clfile,"%sbestLCDMflat6par-cmb-wmap5-meanL.cl",dirs[3]); // WMAP5 power spectrum LCDM best fit to CMB only (Dunkley et al. 2008 p.16 tab 2 "5 Year Mean"), marginalized mean  likelihood Cl; use this for normal GRF sims ; see /home/blew/BIG/gaussian_simulations/wmap5-tests/README for more info on this
//  }
//
//  if (wmap_data == "wmap5ILCfit") {
//    sprintf(Clfile,"%sILC/wmap_ilc_5yr_v3-binup-inter-join",dirs[3]); // the cubic spline fit to the binned full sky ILC5 map for l>30 and bestLCDM as in Dunkley (see above) for l<=30 is used. see README ( /home/blew/programy/Mscs/WMAPstuff/ILC/README ) file for more info on this.
//  }
//
//
///*   sprintf(Clfile,"%scomb-wmap3.cl",dirs[3]); // wmap reconstructed Cls from wmap3 */
///*   sprintf(Clfile,"%swmap3-sim-mix100.cl",dirs[3]); // wmap reconstructed Cls from wmap3 upto l=100 and from l=101... simulated from lcdm best fit */
///*   sprintf(Clfile,"%swmap3-lowltune2",dirs[3]); // wmap reconstructed Cls from wmap3 upto l=1024; the negative values exist for high l end but are put 0 when doing alms and l=3,5,7 are changed so that better match the pseudo Cls and variance in the sims. to the data */
//
//// EXPERIMANTALL -  stuff for Naoshi END
///*   sprintf(Clfile,"%snaoshi-open.cl",dirs[3]); */
///*   sprintf(Clfile,"%snaoshi-closed.cl",dirs[3]); */
///*   sprintf(Clfile,"%snaoshi-highB.cl",dirs[3]); */
///*   sprintf(Clfile,"%snaoshi-lowB.cl",dirs[3]); */
///*   sprintf(Clfile,"%snaoshi-noCDM-flat-ObeqOm.cl",dirs[3]); */
///*   sprintf(Clfile,"%snaoshi-noCDM-flat-ObSTDOlBIG.cl",dirs[3]); */
//// EXPERIMANTALL -  stuff for Naoshi END
//
//  //sprintf(Clfile,"%sbestWMAP-hitau",dirs[2]);
//
//  //  dirty power spectrum from WMAP data according to the channel
//  //sprintf(Clfile,"%s256-512-wmap-%s-%s-kp0-sm_nosmooth-Fmet_8",dirs[3],DA,DA);
//
//  // NASA power spectrum all the same for all channels
//  //sprintf(Clfile,"%swmap_INCNASA-power.cl",dirs[3]);
///*------------------------------------------------------------------------------------------*/
//  //
//  // definition of sigma0 for noice generation for given DA
//  //
//  if (wmap_data == "wmap1") sigma0offset=0;
//  if (wmap_data == "wmap13" || wmap_data == "wmap3") sigma0offset=10;
//  if (wmap_data == "wmap5") sigma0offset=20;
//
//  if (wmap_data == "wmap5ILCfit") sigma0=0.0;
//  if (strcmp(DA,"k1") == 0) { sigma0 = sigma0_WMAP[0+sigma0offset]; }
//  if (strcmp(DA,"k2") == 0) { sigma0 = sigma0_WMAP[1+sigma0offset]; }
//  if (strcmp(DA,"q1") == 0) { sigma0 = sigma0_WMAP[2+sigma0offset]; }
//  if (strcmp(DA,"q2") == 0) { sigma0 = sigma0_WMAP[3+sigma0offset]; }
//  if (strcmp(DA,"v1") == 0) { sigma0 = sigma0_WMAP[4+sigma0offset]; }
//  if (strcmp(DA,"v2") == 0) { sigma0 = sigma0_WMAP[5+sigma0offset]; }
//  if (strcmp(DA,"w1") == 0) { sigma0 = sigma0_WMAP[6+sigma0offset]; }
//  if (strcmp(DA,"w2") == 0) { sigma0 = sigma0_WMAP[7+sigma0offset]; }
//  if (strcmp(DA,"w3") == 0) { sigma0 = sigma0_WMAP[8+sigma0offset]; }
//  if (strcmp(DA,"w4") == 0) { sigma0 = sigma0_WMAP[9+sigma0offset]; }
//
//  //
//  // definition of the beam tf^2 of the experiment i.e. the window function due to beam
//  //
//  if (wmap_data == "wmap1" || wmap_data == "wmap3" || wmap_data == "wmap13") {  sprintf(bl_file,"%sbl-%s.txt",dirs[3],DA); }
//  if (wmap_data == "wmap5") {  sprintf(bl_file,"%swmap_%s_ampl_bl_5yr_v3.txt",dirs[3],DA); }
//  if (wmap_data == "wmap5ILCfit") {  sprintf(bl_file,""); }
//
//  /*   if ((strcmp(DA,"q1") == 0) || (strcmp(DA,"q2") == 0)) { smooth_user = thFWHM_Q; } */
//  /*   if ((strcmp(DA,"v1") == 0) || (strcmp(DA,"v2") == 0)) { smooth_user = thFWHM_V; } */
//  /*   if ((strcmp(DA,"w1") == 0)  || (strcmp(DA,"w2") == 0)  || (strcmp(DA,"w3") == 0)  || (strcmp(DA,"w4") == 0) )  { smooth_user = thFWHM_W; } */
//
//  //
//  // definition of the file with number of observations of a given pixel
//  //
//  if (wmap_data == "wmap1" || wmap_data == "wmap3" || wmap_data == "wmap13") {
//    sprintf(WMAPNobsfile,"%swmap",dirs[3]);  sprintf(WMAPNobsfile3,"%swmap3sum",dirs[3]);
//    if (strcmp(DA,"k1") == 0) { strcat(WMAPNobsfile,"-k1"); strcat(WMAPNobsfile3,"-k1"); strcpy(DA1,"k"); DAi1=1; }
//    if (strcmp(DA,"k2") == 0) { strcat(WMAPNobsfile,"-k2"); strcat(WMAPNobsfile3,"-k2"); strcpy(DA1,"k"); DAi1=2; }
//    if (strcmp(DA,"q1") == 0) { strcat(WMAPNobsfile,"-q1"); strcat(WMAPNobsfile3,"-q1"); strcpy(DA1,"q"); DAi1=1; }
//    if (strcmp(DA,"q2") == 0) { strcat(WMAPNobsfile,"-q2"); strcat(WMAPNobsfile3,"-q2"); strcpy(DA1,"q"); DAi1=2; }
//    if (strcmp(DA,"v1") == 0) { strcat(WMAPNobsfile,"-v1"); strcat(WMAPNobsfile3,"-v1"); strcpy(DA1,"v"); DAi1=1; }
//    if (strcmp(DA,"v2") == 0) { strcat(WMAPNobsfile,"-v2"); strcat(WMAPNobsfile3,"-v2"); strcpy(DA1,"v"); DAi1=2; }
//    if (strcmp(DA,"w1") == 0) { strcat(WMAPNobsfile,"-w1"); strcat(WMAPNobsfile3,"-w1"); strcpy(DA1,"w"); DAi1=1; }
//    if (strcmp(DA,"w2") == 0) { strcat(WMAPNobsfile,"-w2"); strcat(WMAPNobsfile3,"-w2"); strcpy(DA1,"w"); DAi1=2; }
//    if (strcmp(DA,"w3") == 0) { strcat(WMAPNobsfile,"-w3"); strcat(WMAPNobsfile3,"-w3"); strcpy(DA1,"w"); DAi1=3; }
//    if (strcmp(DA,"w4") == 0) { strcat(WMAPNobsfile,"-w4"); strcat(WMAPNobsfile3,"-w4"); strcpy(DA1,"w"); DAi1=4; }
//  }
//  if (wmap_data == "wmap5") {
//    sprintf(WMAPNobsfile,"%swmap5sum",dirs[3]);
//    if (strcmp(DA,"k1") == 0) { strcat(WMAPNobsfile,"-k1"); strcat(WMAPNobsfile3,"-k1"); strcpy(DA1,"k"); DAi1=1; }
//    if (strcmp(DA,"k2") == 0) { strcat(WMAPNobsfile,"-k2"); strcat(WMAPNobsfile3,"-k2"); strcpy(DA1,"k"); DAi1=2; }
//    if (strcmp(DA,"q1") == 0) { strcat(WMAPNobsfile,"-q1"); strcat(WMAPNobsfile3,"-q1"); strcpy(DA1,"q"); DAi1=1; }
//    if (strcmp(DA,"q2") == 0) { strcat(WMAPNobsfile,"-q2"); strcat(WMAPNobsfile3,"-q2"); strcpy(DA1,"q"); DAi1=2; }
//    if (strcmp(DA,"v1") == 0) { strcat(WMAPNobsfile,"-v1"); strcat(WMAPNobsfile3,"-v1"); strcpy(DA1,"v"); DAi1=1; }
//    if (strcmp(DA,"v2") == 0) { strcat(WMAPNobsfile,"-v2"); strcat(WMAPNobsfile3,"-v2"); strcpy(DA1,"v"); DAi1=2; }
//    if (strcmp(DA,"w1") == 0) { strcat(WMAPNobsfile,"-w1"); strcat(WMAPNobsfile3,"-w1"); strcpy(DA1,"w"); DAi1=1; }
//    if (strcmp(DA,"w2") == 0) { strcat(WMAPNobsfile,"-w2"); strcat(WMAPNobsfile3,"-w2"); strcpy(DA1,"w"); DAi1=2; }
//    if (strcmp(DA,"w3") == 0) { strcat(WMAPNobsfile,"-w3"); strcat(WMAPNobsfile3,"-w3"); strcpy(DA1,"w"); DAi1=3; }
//    if (strcmp(DA,"w4") == 0) { strcat(WMAPNobsfile,"-w4"); strcat(WMAPNobsfile3,"-w4"); strcpy(DA1,"w"); DAi1=4; }
//  }
//
//  if (wmap_data == "wmap5ILCfit") {      sprintf(WMAPNobsfile,""); }
//
//  //
//  // definition of the gaussian RNG to be used for the run
//  //
//  gaussgenerator = 0; // this parameter is not used
//  //save_partial = 2; // 0 - do not save the partial files 1- save sum l, 2 - save both sum and multipoles
///*   sprintf(sim_dir,"%sDA",WMAPsim_dir); */
//
//  //
//  // definition and check of the output directory
//  //
//  sprintf(sim_dir,"%sDA",dirs[8]);
//  if (strcmp(DA,"k1") == 0) { strcat(sim_dir,"-k1"); }     if (strcmp(DA,"k2") == 0) { strcat(sim_dir,"-k2"); }
//  if (strcmp(DA,"q1") == 0) { strcat(sim_dir,"-q1"); }     if (strcmp(DA,"q2") == 0) { strcat(sim_dir,"-q2"); }
//  if (strcmp(DA,"v1") == 0) { strcat(sim_dir,"-v1"); }     if (strcmp(DA,"v2") == 0) { strcat(sim_dir,"-v2"); }
//  if (strcmp(DA,"w1") == 0) { strcat(sim_dir,"-w1"); }     if (strcmp(DA,"w2") == 0) { strcat(sim_dir,"-w2"); }
//  if (strcmp(DA,"w3") == 0) { strcat(sim_dir,"-w3"); }     if (strcmp(DA,"w4") == 0) { strcat(sim_dir,"-w4"); }
//  strcat(sim_dir,"/");
//  if (wmap_data == "wmap5ILCfit") {  sprintf(sim_dir,"%s",dirs[28]); }
//  if (!Mscs_contains_str(work_scheme,"quiet!",6)) { sprintf(tmp,"if [ ! -d %s ]; then mkdir -p %s; fi\n",sim_dir,sim_dir); system(tmp); } // create the simulation directory for DAs
//
//  //
//  // definition of the output file name for a new simulated map ( this should be connected with the make_detailed_file_name routine from Mscs-common )
//  //
//  if (wmap_data == "wmap1" || wmap_data == "wmap3" || wmap_data == "wmap13") {
//    sprintf(tmp,"sim%li_wmap1",sim_num);
//    sprintf(files_prefix1,"%li-%li-%s-%s%li-%s%li-%s-sm_%s-Fmet_%i",nside_loc,lmax_loc,tmp,DA1,DAi1,DA1,DAi1,"nomask","no",fourier_method);
//    sprintf(tmp,"sim%li_wmap3",sim_num);
//    sprintf(files_prefix3,"%li-%li-%s-%s%li-%s%li-%s-sm_%s-Fmet_%i",nside_loc,lmax_loc,tmp,DA1,DAi1,DA1,DAi1,"nomask","no",fourier_method);
//  }
//  if (wmap_data == "wmap5") {
//    sprintf(tmp,"sim%li_wmap5",sim_num);
//    sprintf(files_prefix1,"%li-%li-%s-%s%li-%s%li-%s-sm_%s-Fmet_%i",nside_loc,lmax_loc,tmp,DA1,DAi1,DA1,DAi1,"nomask","no",fourier_method);
//  }
//  if (wmap_data == "wmap5ILCfit") {
//    sprintf(tmp,"sim%li_wmap5_ILCfit",sim_num);
//    sprintf(files_prefix1,"%li-%li-%s-%s%li-%s%li-%s-sm_%s-Fmet_%i",nside_loc,lmax_loc,tmp,"",0,"",0,"nomask","no",fourier_method);
//  }
//
//  sprintf(tmp,"sim%li",sim_num);
//  sprintf(files_prefix,"%li-%li-%s-%s%li-%s%li-%s-sm_%s-Fmet_%i",nside_loc,lmax_loc,tmp,DA1,DAi1,DA1,DAi1,"nomask","no",fourier_method); // this thing is actually not needed --- to be removed some day (below it's used only to name the file with the signal only map if requested but  some of the files_prefix[13] filenames could be used instead
//
//  //
//  // definition of the file with alms for the new map (existing or not)
//  //
//  if (wmap_data == "wmap5ILCfit") {    sprintf(tmp,"sim%li_ILC5fit",sim_num); } else   { sprintf(tmp,"sim%li",sim_num); }
//  if (strcmp(how,"same_alms") == 0)  sprintf(alms_file,"%s%li-%li-%s-%s%i-%s%i-%s-sm_%s-Fmet_%i",dirs[21],nside_loc,lmax_loc,tmp,"",0,"",0,"nomask","no",0);
//  else { printf("Mscs-map.c file> !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n THIS OPTION (SOMETHING ELSE THAN same_alms) IS NOT USUALLY USED; BE CAREFUL OF WHAT YOU'RE DOING AND JUST IN CASE CHECK THE CODE\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
//    sprintf(alms_file,"%s%li-%li-%s-%s%li-%s%li-%s-sm_%s-Fmet_%i",dirs[21],nside_loc,lmax_loc,tmp,DA1,DAi1,DA1,DAi1,"nomask","no",0); }
//  if (!Mscs_contains_str(work_scheme,"quiet!",6)) { sprintf(tmp,"if [ ! -d %s ]; then mkdir -p %s; fi\n",dirs[21],dirs[21]); system(tmp); }// create the alms directory
//
//  //sprintf(tmp,"if [ ! -d %s%s_partial_files ]; then mkdir -p %s%s_partial_files; fi\n",sim_dir,files_prefix,sim_dir,files_prefix); system(tmp); // create the destination directory for partial files
//
//  //
//  // definition of the file with precomputed Ylm functions
//  //
///*   sprintf(plms_file,"%sPmthls-%li-%li",data_dir,nside_loc,2*nside_loc); */
//  if (fourier_method == 7) sprintf(plms_file,"%sPlms-%li-%li",dirs[7],nside_loc,2*nside_loc);
//  if ((fourier_method == 8) ||  (fourier_method == 9)) sprintf(plms_file,"%sPmthls-%li-%li",dirs[7],nside_loc,2*nside_loc);
//  if ((fourier_method == 7) || (fourier_method == 8) ||  (fourier_method == 9)) printf("-- Plms file: %s\n",plms_file);
//  //sprintf(plms_file,"");
//  strcpy(show_process,"");
//  //strcpy(show_process,"show");
//  printf("-- simulation directory: %s\n",sim_dir);
//  printf("-- simulation files prefix: %s\n",files_prefix);
//
//  //
//  // run the thing
//  //
//  generate_full_signal_noice_gaussian_map(nside_loc,lmax_loc,Clfile,sigma0,bl_file,smooth_user,WMAPNobsfile,WMAPNobsfile3,gaussgenerator,save_partial,sim_dir,alms_file,files_prefix,files_prefix1,files_prefix3,plms_file,show_process,how,fourier_method,wmap_data,work_scheme);
//
//}
//
//// nside_loc - the nside of the map
//// lmax_loc - the maximal multipole included in the map
//// Clfile - the name of the file containing the power spectrum for the map (in the cmbfast format)
//// sigma0 - noice per DA information
//// smooth_beam - the size of the [gaussian or other] smoothing beam of the instrument ( smooth_user = -1  means don't do the smoothing)
//// smooth_user - gaussian smoothing required by user
///* // pixtf - pixel transfer function to be used [ normally for given nside ] -- no such thing for now */
//// WMAPNobsfile - the file containing the nuber of observations in each pixel in the map corresponding to the WMAP observations (depends on the DA) - fits file format
//// gaussgenerator -- accepted but not used !! the method 5 of generate_gaussian_alms is used
//// save_partial [1/0] whether to save the partial files of individial multipoles or not
//// partial_files_prefix - the prefix name of the partial files and output files
//// plms_file - file containing the precomputerd Legendre associated polynominals (in format of the program precalculate_legendre_polynominals
//// show_process ["show"/"sth else"] - wheter to show or not the process of building the map
//// how - the info on how the simulation process should be realized. if how is empty then by default the simulations are done independently from DAs
//// but if it's not empty then it is treated as the path to the filename with alms generated for some particular simulation for all DAs.
//// thi is usefull since the shape of the CMB genus is independent from the DA they are watched by. The difrerences arise only from the instrumental reasons.
//// work_scheme - this defines various ways of what is supposed to be done during this run using keywords passed as strings.
////               use eg. _ or " " space, to separate keywords
//// work_scheme - is dedicated for cluster use, where the I/O traffic should be maximally minimalized. no partial files are stored -
////               quiet - makes that nothing is stored on the hard disk (except the alms and signal power spectrum). but careful - if you choose quiet then the defauld data to generate is wmap3. !!
////               quiet! - makes things really quiet and the alms and powr spectrum are not stored either
////               savesignal -- use to save signal map to the disk as well (regardless of the quientness parameters now)
////               savenoise1,savenoise3,savenoise5 - to save noise maps to disk
////               pixtf - keyword to use the pixel tf. for SHT
////               almsonly - generate and save the alms only ; dont generate the maps
///* for wmap5 it is not possible to generate more than just wmap5 simulation at one run because the beams are different and sigma0 are also  */
///* different so WMAPNobsfile is used to keep the informaiton on the file with Nobs. Also there will be a different best fit Cl  model */
//
///* in fact for the simultaneous generation of sims for wmap1 and wmap3 the same sigma0 was used - i.e. sigma0 for wmap3  for noise generation*/
//
///* for wmap5 the output files prefix files_prefix1 is used */
//void mscsGaussianSimulation::generate_full_signal_noice_gaussian_map(long nside_loc, long lmax_loc, strarg Clfile, double sigma0, strarg smooth_beam, double smooth_user, strarg WMAPNobsfile, strarg WMAPNobsfile3, int gaussgenerator, int save_partial, strarg sim_dir, strarg alms_file, strarg files_prefix, strarg files_prefix1, strarg files_prefix3, strarg plms_file, strarg show_process, strarg how, int fourier_method, string wmap_data,string work_scheme) {
//  long i,l,m;
//  a_lm a,b;
//  double s,x,s3,x3,s5,x5;
//  double * Nobs=NULL, *Smap=NULL, *Nobs3=NULL, *Smap3=NULL, *noise3=NULL, *Nobs5=NULL, *Smap5=NULL, *noise5=NULL;
//  filenamestr full_path, full_path_partial,  Nobs_map_file, noise_map_file,signal_noise_map_file;
//  filenamestr tmp;
//  FILE* f=NULL;
//  string fname;
//  bool wmap1,wmap3, wmap5;
//  double tmpClratio,tmpCl;
//  double pixtf;
//
//  if (Mscs_contains_str(work_scheme,"pixtf",5)) pixtf=(double)nside_loc; else pixtf=1.0;
//
///*   if (work_scheme == "quiet" && wmap_data == "wmap13") wmap_data = "wmap3"; // can't work quietly with more than one set of data */
//  if (Mscs_contains_str(work_scheme,"quiet",5) && wmap_data == "wmap13") wmap_data = "wmap3"; // can't work quietly with more than one set of data
//
//  if (wmap_data == "wmap1" || wmap_data == "wmap13")  wmap1 = true;  else wmap1 = false;
//  if (wmap_data == "wmap3" || wmap_data == "wmap13")  wmap3 = true;  else wmap3 = false;
//  if (wmap_data == "wmap5")  { wmap5 = true; wmap3 = false; wmap1 = false; }  else wmap5 = false;
//
//  //if (smooth_user > 0) { strcpy(sm,"smooth"); } else { strcpy(sm,"nosmooth"); }
//
//
///*  OLD STUFF -- BEGIN */
//  sprintf(full_path_partial,"%s%s_partial_files/%s",sim_dir,files_prefix,files_prefix);
//  if (save_partial != 0) { sprintf(tmp,"if [ ! -d %s%s_partial_files ]; then mkdir -p %s%s_partial_files; fi\n",sim_dir,files_prefix,sim_dir,files_prefix); system(tmp); } // create the destination directory for partial files
///*  OLD STUFF -- END */
//
//
//  //
//  //map signal generation
//  //
//  printf("|%s> * generation of full signal and noise gaussian map\n\n",object_name.c_str());
//  printf("|%s>  -- generating signal gaussian map\n",object_name.c_str());
//
///*   kill_map(); */
//  makekill_space_manager("kill","T",1);  makekill_space_manager("kill","m",1);
//
//  if (!Mscs_contains_str(work_scheme,"quiet",5) || alms_loaded == 0) {
//    set_map_nside(nside_loc);
//    if (alm!=NULL) delete alm;
//    alm=new mscsAlms;
//    alm->set_lmax(lmax_loc);
//
//    if (strcmp(how,"same_alms") == 0) {
//      sprintf(tmp,"%s-Fal%li-bin",alms_file,alm->get_lmax());
//      printf("|%s>  -- looking if there is a alms file for this simulation (%s)...",object_name.c_str(),tmp);
//      f = fopen(tmp,"r"); }
//    //else { sprintf(tmp,"%s",full_path); }
//
//    // CURRENT VERSION
//    if (f != NULL) {
//      if (strcmp(how,"same_alms") == 0) {  printf(" yes, there is\n"); fclose(f); }
//      if (Mscs_contains_str(work_scheme,"almsonly",8)) { printf(" WARNING!! the alms file already exist, and almsonly requested: skipping to the next simulation\n"); return; }
//      alm->loadbinAlms(alms_file,1);
//    }
//    // CURRENT VERSION
//
//// EXPERIMANTALL -  stuff for Naoshi BRGIN -- for keeping the phases while changing the power spectrum
///*     if (f != NULL) { if (strcmp(how,"same_alms") == 0) {  printf(" yes, there is\n"); fclose(f);} loadbinAlms(alms_file,1); */
///*       printf("*** PERFORMING EXPERIMENTAL BLOCK: recalibrating the power spectrum while keeping phases unchanged ***\n"); */
///*       // new stuff here */
///*       // load the new power spectrum to which we should fit with the same phases */
///*       //loadtxtC_l("/home/blew/programy/cmbfast-new/cmbfast4.5.1/naoshi-open.cl",2); */
///*       sprintf(tmp,"%sbestLCDMflat6par.cl",PROGRAM_DIRS[3]); // current power spectrum LCDM best fit ; use this for normal GRF sims */
///*       loadtxtC_l(tmp,2); */
///*       // use proper calibration if needed */
///*       for (i=0;i<=lmax;i++) { C_l[i][1] = C_l[i][1] * 23.7e-10* 2.726*2.726; } // we simulate WMAP in [K] temperature maps so we need to multiply the Cl by T_CMB^2 and there is a correction factor between the simulated from best fit LCDM 6params model of the Hinshlaw 2006 et al. and the obtained value from the cmbfast of the sigma8 parameter. this is the raio of the two by which the Cl is corrected to better fit the data. -- this is current one */
//
///*       // recalibrate the alms accordingly */
//
///*       for (l=0; l<=lmax; l++) { */
///* 	tmpCl=calculate_C_l(l); */
///* 	if ( tmpCl != 0) tmpClratio=sqrt(C_l[l][1]/tmpCl); else tmpClratio = sqrt(C_l[l][1]); // if the original is zero it will stay zero cause there's no information about phases */
///* 	printf("tmpCLratio: sqrt(new/original)=%lE\n",tmpClratio); */
///* 	for ( m=-l; m<=l; m++) { */
///* 	  a = get_F(l,m); a.R*=tmpClratio; a.I*=tmpClratio; */
///* 	  if (l==0 || l==1) { a.R=0; a.I=0; } // remove the l=0,1 if you like */
///* 	  set(l,m,a);	   */
///* 	} */
///* 	printf("l=%li\n",l); */
///*       } */
//
///*       // check the result; */
///*       calculate_C_l(0,lmax,1); */
///*       savetxtC_l("check.tmp",2); */
///*       printf("*** PERFORMING EXPERIMENTAL BLOCK DONE:  ***\n"); */
///*     } */
//// EXPERIMANTALL -  stuff for Naoshi END
//
//    else {   if (strcmp(how,"same_alms") == 0)   printf("no, there isn't -- I will generate it\n");
//
//      if (wmap_data == "wmap5ILCfit" ) {
//	C_l->load(Clfile); // load the power spectrum in Mscs format format - unnormalized
//	//loadtxtC_l(Clfile,1); // load the power spectrum in Mscs format format - unnormalized   //commented out during transition into version-1.0
//      }
//      else {
//	C_l->load(Clfile); // load the power spectrum in CMBfast format
//	C_l->divide_llpotwoPI();
//	//loadtxtC_l(Clfile,2); // load the power spectrum in CMBfast format   //commented out during transition into version-1.0
//      }
//      //*****************EXPERIMENTAL STUFF ****************
//      // the power spectrum in different chanels ...
//
//      /*     if ((sigma0 == sigma0_Q1) || (sigma0 == sigma0_Q2))        */
//      /*       { for (i=0;i<=1000;i++) { C_l[i][1] *= 14.251; } // value estimated from fit to WMAP in q1 band */
//      /* 	printf("----------------> EXPERIMENTAL : Cl x for Q channel: sigma0 = %lE\n",sigma0); } */
//      /*     if ((sigma0 == sigma0_V1) || (sigma0 == sigma0_V2))        */
//      /*       { for (i=0;i<=1000;i++) { C_l[i][1] *= 15.622; } // value estimated from fit to WMAP in v1 band */
//      /* 	printf("----------------> EXPERIMENTAL : Cl x for V channel: sigma0 = %lE\n",sigma0); } */
//      /*     if ((sigma0 == sigma0_W1) || (sigma0 == sigma0_W2) || (sigma0 == sigma0_W3) || (sigma0 == sigma0_W4))  */
//      /*       { for (i=0;i<=1000;i++) { C_l[i][1] *= 31.8279; } // value estimated from fit to WMAP in w2 band */
//      /* 	printf("----------------> EXPERIMENTAL : Cl x for W channel:  sigma0 = %lE\n",sigma0); } */
//
//
//
//      // CHOOSE CALIBRATION OF THE POWER SPECTRUM
//
//      //for (i=0;i<=1000;i++) { C_l[i][1] *= 1e-12; } // changing to K^2 from muK^2 // this is for the rough WMAP Cl
//      //for (i=0;i<=lmax;i++) { C_l[i][1] *= 8.0; } // this is for the model Cl best WMAP but calibrated by eye to fit the data. I use this one to hava a statistical approach to low quadrupol estimation; x8 factor is due to the T_cmb^2 * A uncertainty
///*       for (i=0;i<=lmax;i++) { C_l[i][1] = C_l[i][1] * 0.84/0.744 * 2.726*2.726; } // we simulate WMAP in [K] temperature maps so we need to multiply the Cl by T_CMB^2 and there is a correction factor between the simulated from best fit LCDM 6params model of the Hinshlaw 2006 et al. and the obtained value from the cmbfast of the sigma8 parameter. this is the raio of the two by which the Cl is corrected to better fit the data. -- this is old one */
//
//// CURRENT VERSIONS
//      if (wmap_data == "wmap13" || wmap_data == "wmap3") {
//	//commented out during transition into version-1.0
//	// 	for (i=0;i<=lmax;i++) { C_l[i][1] = C_l[i][1] * 23.7e-10 * 2.726*2.726; } // we simulate WMAP in [K] temperature maps so we need to multiply the Cl by T_CMB^2 and there is a correction factor between the simulated from best fit LCDM 6params model of the Hinshlaw 2006 et al. and the obtained value from the cmbfast of the sigma8 parameter. this is the raio of the two by which the Cl is corrected to better fit the data. -- this is current one
//	(*C_l)*= (23.7e-10 * 2.726*2.726);
//      }
//
//      if (wmap_data == "wmap5" ) {
//	//commented out during transition into version-1.0
//	// 	for (i=0;i<=lmax;i++) { C_l[i][1] = C_l[i][1] * 24.1e-10 * 2.726*2.726; } // we simulate WMAP in [K] temperature maps so we need to multiply the Cl by T_CMB^2 and there is a correction factor between the simulated from best fit LCDM 6params model WMAP3 best fit parameters table from lambda including lensing and SZ power. this particular deriviation of the cl from cmbfast does not include the lensing however. INCLUDE IT LATTER   -- this is current one
//	(*C_l)*= (24.1e-10 * 2.726*2.726);
//      }
//
///*       if (wmap_data == "wmap5ILCfit" ) { dont do anything } */
//// CURRENT VERSIONS
//
///*       if (wmap_data == "wmap5" ) { */
///* 	for (i=0;i<=lmax;i++) { C_l[i][1] = C_l[i][1] * 1e-12; } // we simulate WMAP in [K] temperature maps so we need to multiply the Cl by T_CMB^2 and there is a correction factor between the simulated from best fit LCDM 6params model WMAP3 best fit parameters table from lambda including lensing and SZ power. this particular deriviation of the cl from cmbfast does not include the lensing however. INCLUDE IT LATTER   -- this is current one */
///*       } */
//
///*       for (i=0;i<=lmax;i++) { C_l[i][1] = C_l[i][1] * 1e-12; } // this is for the WAMP COMP POW. SP. -- ONLY CONVERSION FROM muK^2 --> K^2 */
//
//
//
//// POWER SPECTRUM CORRECTIONS BEGIN - NORMALLY COMMENTED OUT
//
//// fit to the power reconstructed from the ILC 5yr map
///*       C_l[2][1]=2.540796E-10; */
///*       C_l[3][1]=5.511079E-10; */
//// fit to the power reconstructed from the TOH 5yr map
//
//  //commented out during transition into version-1.0
////       C_l[2][1]=2.119171E-10;
////       C_l[3][1]=3.920900E-10;
//      C_l->set_Cl(2, 2.119171E-10);
//      C_l->set_Cl(3, 3.920900E-10);
//
///*       C_l[2][1]=2.23053E-10; // these values are taken from the reconstructed power spectrum using ML  */
///*       C_l[3][1]=5.44019E-10; */
//// POWER SPECTRUM CORRECTIONS END - NORMALLY COMMENTED OUT
//
//
//
//
//// EXPERIMANTALL -  stuff for Naoshi BEGIN
//// for open model
////      for (i=0;i<=lmax;i++) { C_l[i][1] = C_l[i][1] * 23.7e-10 * 2.726*2.726 * 0.4737/0.83113; } // we simulate WMAP in [K] temperature maps so we need to multiply the Cl by T_CMB^2 and there is a correction factor between the simulated from best fit LCDM 6params model of the Hinshlaw 2006 et al. and the obtained value from the cmbfast of the sigma8 parameter. this is the raio of the two by which the Cl is corrected to better fit the data. -- this is for Naoshi
//
//// for closed model
////      for (i=0;i<=lmax;i++) { C_l[i][1] = C_l[i][1] * 23.7e-10 * 2.726*2.726 * 0.4737/0.57021; } // we simulate WMAP in [K] temperature maps so we need to multiply the Cl by T_CMB^2 and there is a correction factor between the simulated from best fit LCDM 6params model of the Hinshlaw 2006 et al. and the obtained value from the cmbfast of the sigma8 parameter. this is the raio of the two by which the Cl is corrected to better fit the data. -- this is for Naoshi
//
//// for high baryon model
////      for (i=0;i<=lmax;i++) { C_l[i][1] = C_l[i][1] * 23.7e-10 * 2.726*2.726 * 0.4737/0.48351; } // we simulate WMAP in [K] temperature maps so we need to multiply the Cl by T_CMB^2 and there is a correction factor between the simulated from best fit LCDM 6params model of the Hinshlaw 2006 et al. and the obtained value from the cmbfast of the sigma8 parameter. this is the raio of the two by which the Cl is corrected to better fit the data. -- this is for Naoshi
//
//// for low baryon model
////      for (i=0;i<=lmax;i++) { C_l[i][1] = C_l[i][1] * 23.7e-10 * 2.726*2.726 * 0.4737/0.4766; } // we simulate WMAP in [K] temperature maps so we need to multiply the Cl by T_CMB^2 and there is a correction factor between the simulated from best fit LCDM 6params model of the Hinshlaw 2006 et al. and the obtained value from the cmbfast of the sigma8 parameter. this is the raio of the two by which the Cl is corrected to better fit the data. -- this is for Naoshi
//
//// for noCDM ObeqOm
////      for (i=0;i<=lmax;i++) { C_l[i][1] = C_l[i][1] * 23.7e-10 * 2.726*2.726 * 0.4737/0.4939; } // we simulate WMAP in [K] temperature maps so we need to multiply the Cl by T_CMB^2 and there is a correction factor between the simulated from best fit LCDM 6params model of the Hinshlaw 2006 et al. and the obtained value from the cmbfast of the sigma8 parameter. this is the raio of the two by which the Cl is corrected to better fit the data. -- this is for Naoshi
//
//// for noCDM ObSTDOlBIG
///*       for (i=0;i<=lmax;i++) { C_l[i][1] = C_l[i][1] * 23.7e-10 * 2.726*2.726 * 0.4737/0.62839; } // we simulate WMAP in [K] temperature maps so we need to multiply the Cl by T_CMB^2 and there is a correction factor between the simulated from best fit LCDM 6params model of the Hinshlaw 2006 et al. and the obtained value from the cmbfast of the sigma8 parameter. this is the raio of the two by which the Cl is corrected to better fit the data. -- this is for Naoshi */
//
//
//// EXPERIMANTALL -  stuff for Naoshi END
//
//      //*****************EXPERIMENTAL STUFF ****************
//
//
//
//
//      // CHOOSE THE WAY OF GENERATION OF ALMS
//      alm->generate_gaussian_alms("RI","unnormCl","nol0-1",0,1,0,alm->get_lmax(),5,-1); // current one - random gaussian realization of the power spectrum
//
///*       C_l[0][0] = 6.8e-11; C_l[0][1] = 8.6e-11; // this is an "by hand" fix to the WMAP l=0,1 power estimation for dodec project only - calibration as in the combtt C_l */
///*       generate_gaussian_alms("RI","unnormCl","",0,1,lmax,4,-1); // stick to the exact power spectrum */
///*       generate_gaussian_alms("RI","unnormCl","nol0-1",0,1,lmax,4,-1); // stick to the exact power spectrum but remove l=0,1 */
//
//
//      // EXPERIMENTAL - BEGIN -- making aligned and planar alms by hand
///*       printf("********************* PERFORMING THE EXPERIMENTAL BLOCK: START -- making aligned multipoles *********************\n"); */
///*       // make planar and aligned l=2,3 with the same power as those generated/loaded or whatever */
///*       for (l=2;l<=25;l++) { */
///* 	tmpCl = calculate_C_l(l); */
///* 	// clean all alms for a given l. */
///* 	for (m=-l; m<l; m++) {	a.R=0;a.I=0; set(l,m,a); } b=get_F(l,l); */
//
///* 	// set one alm with the correct power and antisymmetrize it */
///* 	a.R=sqrt(tmpCl*(2.0*l+1.0)/2.0)*cos(b.R/b.I); a.I=sqrt(tmpCl*(2.0*l+1.0)/2.0)*sin(b.R/b.I); set(l,l,a); // put all the power into the a_ll yet keep random phases in it */
///* 	a.R=pow(-1.0,(double)l)*a.R; a.I=-pow(-1.0,(double)l)*a.I; set(l,-l,a); // antisymmetrize */
///*       } */
///*       printf("********************* PERFORMING THE EXPERIMENTAL BLOCK: END  *********************\n"); */
//      // REMEMBER TO UNCOMMENT ALSO THE BELOW ROTATE MAP IF YOU LIKE THE MAP TO HAVE SOME PARTICULAR ORIENTATION
//      // EXPERIMENTAL - END
//
//
//
///*       zero_alms_multipole_range(0,5); */
///*       zero_alms_multipole_range(7,1024); */
//      //generate_gaussian_alms("RI","noCl","nol0-1",0,1,lmax,5,-1);
//      if (C_l!=NULL) delete C_l;
//      C_l=alm->calculate_C_l(0,alm->get_lmax(),1);
//      if (!Mscs_contains_str(work_scheme,"quiet!",6)) { // save alms unless you want it really quiet: "quiet!"
//	if (!Mscs_contains_str(work_scheme,"dontsavealms",12)) alm->savebinAlms(alms_file,1);
//	C_l->save(alms_file);   //changed during transition into version-1.0
//      }
//
//      if (Mscs_contains_str(work_scheme,"almsonly",8)) { return; }
//      //exit(0);
//    } // close condition for generation from alms from power spectrum
//  }
//
//
//
//
//
//
//
//
//
//    //*****************EXPERIMENTAL STUFF ****************
//    //calculate_inverse_transformF(8,1,"nosmooth",-1,0,"",plms_file,"");
//  //calculate_inverse_transformF(fourier_method,1,"",smooth_user,save_partial,full_path_partial,plms_file,show_process);
///*   printf("********** REMOVING L=2 **************\n"); */
///*   l=2; */
///*   for (m=-l; m<=l; m++) {	a.R=0;a.I=0; set(l,m,a); } */
///*   printf("********** DONE **************\n"); */
//    //*****************EXPERIMENTAL STUFF ****************
//  SHT_alm2map(fourier_method,alm,alm->get_lmax(),pixtf,smooth_beam,smooth_user,save_partial,full_path_partial,plms_file,show_process);
///*   calculate_inverse_transformF(fourier_method,512,smooth_beam,smooth_user,save_partial,full_path_partial,plms_file,show_process); // includes pixlel tf. */
//
//  // EXPERIMENTAL - BEGIN -- rotating map anisotropical map to a given orienation: useful with "making aligned and planar alms by hand"
///*   rotate_map(30,260,0,true,"T"); */
//  // EXPERIMENTAL - END
//
//  sprintf(full_path,"%s%s-signal",sim_dir,files_prefix);
//  if (Mscs_contains_str(work_scheme,"savesignal",10)) savebinT(full_path,1); // saving the signal map
//
//  if (wmap_data == "wmap5ILCfit" ) {
//    sprintf(full_path,"%s%s",sim_dir,files_prefix1);
//    savebinT(full_path,1);
//  } // saving the signal only ILC5fit map
//
//  //*****************ADDITIONAL DEBUG STUFF ****************
//  // SIGNAL
//  // # var signal, var noise, var S+N, var S theoretical ,var S+N theoretical
///*   f = fopen("debug-test","a"); */
///*   fprintf(f,"%.10lE ",cpeds_variance((*map).T,pix_num)); */
///*   fclose(f); */
//  //*****************ADDITIONAL DEBUG STUFF ****************
//
//  if (wmap_data != "wmap5ILCfit" ) {
//
//  if (wmap1) Smap = new double[pix_num];
//  if (wmap3) Smap3 = new double[pix_num];
//  if (wmap5) Smap5 = new double[pix_num];
//  if (wmap1) for (i=0;i<pix_num;i++) {   Smap[i] = map->T[i]; } // hold the signal map
//  if (wmap3) for (i=0;i<pix_num;i++) {   Smap3[i] = map->T[i]; } // hold the signal map
//  if (wmap5) for (i=0;i<pix_num;i++) {   Smap5[i] = map->T[i]; } // hold the signal map
//
///*   kill_map(); */
//  makekill_space_manager("kill","T",1);  makekill_space_manager("kill","m",1);
//
//  //exit(0);
//
//
//
//
//
///**********************************************************************************************************/
///* noise generation for both WMAP 1 and WMAP 3, the difference comes only from the number of observations */
///**********************************************************************************************************/
//
//  printf("|%s>  -- generating noise gaussian map for requested WMAP DA\n",object_name.c_str());
//  if (wmap1) {
//    loadfitsT(WMAPNobsfile,12); // for I maps: WMAP1
//    change_map_resolution(nside_loc);
//    Nobs = new double[pix_num];  for (i=0;i<pix_num;i++) { Nobs[i] = map->N[i]; map->T[i]=map->N[i];} // hold the Nobs info
///*   sprintf(Nobs_map_file,"%s%s-Nobs",sim_dir,files_prefix1); */
///*   savebinT(Nobs_map_file,1); // saving file with pixel number observations - we don't need this */
//  }
//
//  if (wmap3) {
//    loadfitsT(WMAPNobsfile3,14); // for IQU maps: WMAP3
//    change_map_resolution(nside_loc);
//    Nobs3 = new double[pix_num];  for (i=0;i<pix_num;i++) { Nobs3[i] = map->N[i]; map->T[i]=map->N[i];} // hold the Nobs3 info
///*   sprintf(Nobs_map_file,"%s%s-Nobs",sim_dir,files_prefix3); */
///*   savebinT(Nobs_map_file,1); // saving file with pixel number observations - we don't need this */
//  }
//
//  if (wmap5) {
//    loadfitsT(WMAPNobsfile,12); // for I maps: WMAP5
//    change_map_resolution(nside_loc);
//    Nobs5 = new double[pix_num];  for (i=0;i<pix_num;i++) { Nobs5[i] = map->N[i]; map->T[i]=map->N[i];} // hold the Nobs info
///*   sprintf(Nobs_map_file,"%s%s-Nobs",sim_dir,files_prefix5); */
///*   savebinT(Nobs_map_file,1); // saving file with pixel number observations - we don't need this */
//  }
//
///*   kill_map(); */
//  makekill_space_manager("kill","T",1);  makekill_space_manager("kill","m",1);
//
//  set_map_nside(nside_loc);
//
///*   if (Mscs_contains_str(work_scheme,"savenoise",9)) { */
//    if (wmap3) noise3 = new double[pix_num]; // uncomment if you want to save the noise wmap3
//    if (wmap5) noise5 = new double[pix_num]; // uncomment if you want to save the noise wmap5
///*   } */
//
//  makekill_space_manager("make","T",1); clear_map(map);      clean_mask();
//  for (i=0;i<pix_num;i++) { // keeping the generated noise on the Nobs structure after it's generated (Nobs is no longer needed then)
//    if (wmap1) { s = sigma0/sqrt(Nobs[i]);    x = cpeds_random_gauss_number(0.0,s,100,2); Smap[i] += x; map->T[i]=x; } // making the SN 1 map here
//    if (wmap3) { s3 = sigma0/sqrt(Nobs3[i]);   if (wmap1)  x3 = x*s3/s;  else     x3 = cpeds_random_gauss_number(0.0,s3,100,2);  Smap3[i] += x3; noise3[i] = x3; }// making the SN 3 map here }
//    if (wmap5) { s5 = sigma0/sqrt(Nobs5[i]);    x5 = cpeds_random_gauss_number(0.0,s5,100,2); Smap5[i] += x5;  noise5[i] = x5; } // making the SN 5 map here
//    //*****************EXPERIMENTAL STUFF ****************
//    //s = 200*sigma0*sigma0/sqrt(Nobs[i]);
//    //*****************EXPERIMENTAL STUFF ****************
//
//    //xp = cpeds_random_gauss_numbers(0.0,s,1,2); x = xp[0];
//
//    // recording the noise only in case we want so save it to a file
///*     if (Mscs_contains_str(work_scheme,"savenoise",9)) { */
///*       if (wmap1) map->T[i] = x; //  uncomment if you want to save the noise wmap1 */
///*       if (wmap3) noise3[i] = x3; // uncomment if you want to save the noise wmap3 */
///*       if (wmap5) noise5[i] = x5; // uncomment if you want to save the noise wmap5 */
///*     } */
//  }
//
//
//  // saving noise maps
//  sprintf(noise_map_file,"%s%s-noise",sim_dir,files_prefix1);
//  if (Mscs_contains_str(work_scheme,"savenoise1",10)) savebinT(noise_map_file,1); // saving the noise 1
//
//  if (Mscs_contains_str(work_scheme,"savenoise3",10)) {
//    for (i=0;i<pix_num;i++) { map->T[i] = noise3[i]; }
//    sprintf(noise_map_file,"%s%s-noise",sim_dir,files_prefix3);
//    savebinT(noise_map_file,1); // saving the noise 3
//  }
//  if (wmap3) delete [] noise3;
//
//
//  if (Mscs_contains_str(work_scheme,"savenoise5",10)) {
//    for (i=0;i<pix_num;i++) { map->T[i] = noise5[i]; }
//    sprintf(noise_map_file,"%s%s-noise",sim_dir,files_prefix1);
//    savebinT(noise_map_file,1); // saving the noise 5
//  }
//  if (wmap5) delete [] noise5;
//
//  //*****************ADDITIONAL DEBUG STUFF ****************
//  // NOISE
//  // # var signal, var noise, var S+N, var S theoretical ,var S+N theoretical
///*   f = fopen("debug-test","a"); */
///*   s3=cpeds_variance((*map).T,pix_num); */
///*   fprintf(f,"%.10lE ",s3); */
///*   fclose(f); */
//  //*****************ADDITIONAL DEBUG STUFF ****************
//
//
//  //adding noise to the map (actually it's only a copying to a proper map structure since it's already added on Smap
//  printf("|%s>  -- adding generated noise to the map\n",object_name.c_str());
//  //sprintf(smoothed_map_file,"%s%s"files_prefix); // this line is here since the smoothing is done in the fourier transform at the moment
//
//  if (wmap1) {
//    for (i=0;i<pix_num;i++) {     map->T[i] = Smap[i];  }
//    sprintf(signal_noise_map_file,"%s%s",sim_dir,files_prefix1);
//    if (!Mscs_contains_str(work_scheme,"quiet",5)) savebinT(signal_noise_map_file,1); // saves the signal and noise map
//  }
//  if (wmap3) {
//    for (i=0;i<pix_num;i++) {     map->T[i] = Smap3[i];  }
//    sprintf(signal_noise_map_file,"%s%s",sim_dir,files_prefix3);
//    if (!Mscs_contains_str(work_scheme,"quiet",5)) savebinT(signal_noise_map_file,1); // saves the signal and noise map
//
//  //*****************ADDITIONAL DEBUG STUFF ****************
//  // S+N
//  // # var noise, var signal, var S+N, var S theoretical ,var S+N theoretical
///*   f = fopen("debug-test","a"); */
///*   fprintf(f,"%.10lE ",cpeds_variance((*map).T,pix_num)); */
///*   fclose(f); */
//  //*****************ADDITIONAL DEBUG STUFF ****************
//
//  }
//
//  if (wmap5) {
//    for (i=0;i<pix_num;i++) {     map->T[i] = Smap5[i];  }
//    sprintf(signal_noise_map_file,"%s%s",sim_dir,files_prefix1);
//    if (!Mscs_contains_str(work_scheme,"quiet",5)) savebinT(signal_noise_map_file,1); // saves the signal and noise map
//  }
//
//  if (wmap1) {  delete Nobs;  delete Smap; }
//  if (wmap3) {  delete Nobs3;  delete Smap3; }
//  if (wmap5) {  delete Nobs5;  delete Smap5; }
//
//  //*****************ADDITIONAL DEBUG STUFF ****************
//  // THEORETICAL S, SN
//  // # var noise, var signal, var S+N, var S theoretical ,var S+N theoretical, var Sth_nol012, var S+Nth_nol012
///*   f = fopen("debug-test","a");  */
///*   s=0; */
///*   calculate_C_l(0,lmax,1); */
///*   for (i=0; i<=lmax; i++) { */
///*     s+=(2*C_l[i][0]+1)*C_l[i][1]*bl_tab[8]->get_wl(i) * bl_tab[9]->get_wl(i) * bl_tab[7]->get_wl(i)*bl_tab[8]->get_wl(i) * bl_tab[9]->get_wl(i) * bl_tab[7]->get_wl(i); */
///*   } */
///*   s=s/(4*PI); */
///*   fprintf(f,"%.10lE %.10lE ",s,s+s3); */
//
///*   s=0; */
///*   for (i=3; i<=lmax; i++) { */
///*     s+=(2*C_l[i][0]+1)*C_l[i][1]*bl_tab[8]->get_wl(i) * bl_tab[9]->get_wl(i) * bl_tab[7]->get_wl(i)*bl_tab[8]->get_wl(i) * bl_tab[9]->get_wl(i) * bl_tab[7]->get_wl(i); */
///*   } */
///*   s=s/(4*PI); */
///*   fprintf(f,"%.10lE %.10lE\n",s,s+s3); */
///*   fprintf(f,"#var_N var_S var_SN var_Sth var_SNth var_Sth_nol012 var_SNth_nol012\n"); */
//
///*   fclose(f); */
///*   kill_window_functions(); */
//  //*****************ADDITIONAL DEBUG STUFF ****************
//
//
//  }
//
//}
//
//
//
////************************************************************************
