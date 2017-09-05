#include "Mscs-minkowski-f.h"

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// IMPLEMENTATION OF THE CLASS class minkowski_f
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

// CONSTRUCTOR
minkowski_f::minkowski_f(string ftype, long th_num) {
  long i;
  // initiate the functional space
  f = new v[th_num];
  // initiate functional parameters
  size = th_num;
  min_i = 0;
  max_i = th_num-1;
  func_type = ftype;

  // zero other info
  min_nu = max_nu = min_vnu = max_vnu = 0;
  for (i=min_i;i<=max_i;i++) {     f[i].nu=0.0; f[i].vnu=0.0;  }

  printf("|%s> * initiating minkowski functional object. size: %li\n",func_type.c_str(),th_num);
}

// DESTRUCTOR
minkowski_f::~minkowski_f() {
  delete [] f;
}


// HANDLERS
double minkowski_f::get_nu(long i) {
  double ret;
  if (check_irange(i))  ret = f[i].nu; else ret = 0;
  return ret;
}

double minkowski_f::get_vnu(long i) {
  double ret;
  if (check_irange(i)) ret = f[i].vnu; else ret = 0;
  return ret;
}

void minkowski_f::set_nu(long i, double nu) {
  f[i].nu = nu;
}

void minkowski_f::set_vnu(long i, double vnu) {
  f[i].vnu = vnu;
}

string minkowski_f::get_type() {
  return func_type;
}

void minkowski_f::set_type(string ftype) {
  func_type  = ftype;
}

// this returns the range of table indexes that are allowed
void minkowski_f::get_irange(long *_min_i, long *_max_i) {
  *_min_i = min_i; *_max_i = max_i;
}
// returns the number of thresholds
long minkowski_f::get_size() {
  return size;
}

void minkowski_f::clear() {
/*  to be implemented */
}

void minkowski_f::get_frange(double * _min_nu, double * _max_nu, double * _min_vnu, double * _max_vnu) {
  calculate_functional_minmax_values();
  *_min_nu = min_nu;   *_max_nu = max_nu;
  *_min_vnu = min_vnu;   *_max_vnu = max_vnu;
}

// normalize the minkowski functional
void minkowski_f::normalize_V0() {
  long i;
  long sizele1=size-1;
  for (i=0;i<size;i++) f[i].vnu/=f[sizele1].vnu;

}

void minkowski_f::save(string name) {
  FILE *f;
  long i;
  f=fopen(name.c_str(),"w");
  for (i=0;i<size;i++) {    fprintf(f,"%lE %lE\n",get_nu(i),get_vnu(i)); }
  fclose(f);
}

//--------------------------------------------------------------------------------
// prrivate methods
//--------------------------------------------------------------------------------

// checking index range passed to handlers
bool minkowski_f::check_irange(long i) {
  bool ret;
  if ((i <= max_i) && (i >= min_i)) ret = true;   else ret = false;
  return ret;
}

void minkowski_f::calculate_functional_minmax_values() {
  long i;

  min_nu = f[0].nu;  max_nu = f[0].nu;
  min_vnu = f[0].vnu;  max_vnu = f[0].vnu;

  for (i=min_i;i<=max_i;i++) { 
    if (f[i].nu < min_nu) { min_nu = f[i].nu; min_nui = i; }
    if (f[i].nu > max_nu) { max_nu = f[i].nu; max_nui = i; }
    if (f[i].vnu < min_vnu) { min_vnu = f[i].vnu; min_vnui = i; }
    if (f[i].vnu > max_vnu) { max_vnu = f[i].vnu; max_vnui = i; }
  }

}






















////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// IMPLEMENTATION OF THE CLASS class minkowski_fs
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

// general constructor
minkowski_fs::minkowski_fs() {
  fnum=0;
  mfs=NULL;
}

// some other general constructor
minkowski_fs::minkowski_fs(long num) {
  long i;
  fnum = num;
  mfs=new mfs_type[fnum];
  
  for (i=0;i<fnum;i++) {
    mfs[i].v0=NULL;
    mfs[i].v1=NULL;
    mfs[i].v2=NULL;
    mfs[i].regid=0;
    mfs[i].size=0;
  }
  
}

// constructor for a given common number of thresholds and number of set of functionals
minkowski_fs::minkowski_fs(long num, long thnum) {
  long i;
  string name;
  filenamestr tmpch;

  fnum = num;
  mfs=new mfs_type[fnum];
  
  for (i=0;i<fnum;i++) {
    sprintf(tmpch,"F%li_%s",i,"v0"); name=tmpch; mfs[i].v0 = new minkowski_f(name,thnum);
    sprintf(tmpch,"F%li_%s",i,"v1"); name=tmpch; mfs[i].v1 = new minkowski_f(name,thnum);
    sprintf(tmpch,"F%li_%s",i,"v2"); name=tmpch; mfs[i].v2 = new minkowski_f(name,thnum);
    mfs[i].regid=i;
    mfs[i].size=0;
  }

  object_name="minkowski_functionals";
}

/****************************************************************************************************/


minkowski_fs::~minkowski_fs() {
  fnum=0;
  delete [] mfs;
  mfs=NULL;
}

// returns the number of sets of minkowski functionals stored in the object
long minkowski_fs::get_size() { return fnum; }
void minkowski_fs::set_name(string name) { object_name=name; }
string minkowski_fs::get_name() { return object_name; }

// saves the functionals to a file
void minkowski_fs::save(string name, long how) {
  FILE * f;
  long thres_num,i,j;

  if (how == 1) {
    name+=".mink";
    printf("|%s> * saving the minkowski functionals to file %s\n",object_name.c_str(),name.c_str());
    thres_num=mfs[0].v0->get_size(); // assume the same number of threshods for all functional types and all regions
    f = fopen(name.c_str(),"w");
    fprintf(f,"# %li %li\n",fnum, thres_num);

    for (j=0;j<fnum;j++) { 
      for (i=0;i<thres_num;i++) {
	fprintf(f,"%li %li %lE %lE %lE %lE\n",mfs[j].regid,mfs[j].size, mfs[j].v0->get_nu(i),mfs[j].v0->get_vnu(i),mfs[j].v1->get_vnu(i),mfs[j].v2->get_vnu(i));
      }
    } 
    
    fclose(f);
  }

}

void minkowski_fs::normalize_V0s() {
  long i;
  for (i=0;i<fnum;i++) mfs[i].v0->normalize_V0();

}

double minkowski_fs::get_mf(long i,long imf,long inu) {
  double ret=0;
  if (imf==0) { ret=mfs[i].v0->get_vnu(inu); }
  if (imf==1) { ret=mfs[i].v1->get_vnu(inu); }
  if (imf==2) { ret=mfs[i].v2->get_vnu(inu); }
  return ret;
}
























////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// IMPLEMENTATION OF THE CLASS class minkowski_fs_MFstat
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


/* // constructor for a given common number of thresholds and number of set of functionals */
/* minkowski_fs_MFstat::minkowski_fs_MFstat(long num, long thnum) { */
/*   long i; */
/*   string name; */
/*   filenamestr tmpch; */

/*   fnum = num; */
/*   mfs=new mfs_type[fnum]; */
  
/*   for (i=0;i<fnum;i++) { */
/*     sprintf(tmpch,"F%li_%s",i,"v0"); name=tmpch; mfs[i].v0 = new minkowski_f(name,thnum); */
/*     sprintf(tmpch,"F%li_%s",i,"v1"); name=tmpch; mfs[i].v1 = new minkowski_f(name,thnum); */
/*     sprintf(tmpch,"F%li_%s",i,"v2"); name=tmpch; mfs[i].v2 = new minkowski_f(name,thnum); */
/*     mfs[i].regid=i; */
/*     mfs[i].size=0; */
/*   } */

/*   object_name="minkowski_functionals"; */
/* } */

/****************************************************************************************************/
// This constructor is dedicated to initiate the structures to derive the covariance matrices of the regionally calculated minkowski 
// functionals for a GRF wmap5 simulations. The final product will be therefore: a set of 
// averages <f_nu,i,j,k,l>_j and inverse covariance matrices C_nu1nu2,i,k,l
// where 
// nu1/2 - temperature threshold
// i - MF type (0,1,2)
// j - simulation number (theese can be split between different nodes)
// k - region number (in the multi-mask)
// l - smoothing ID of the simulation (dont use too many smoothing scales to avoid memory allovation error)
// m - multi mask index
// d - dataset index
//
minkowski_fs_MFstat::minkowski_fs_MFstat(long node_no, string sim_type, long sim_stnum, long sim_ennum, long sm_num, double* sm, long mm_num, string* mmf, long th_num, string mask_name) {
  long nu,k,l,m,d,sim_num,reg_num;
  long sim_type_num,cov_size;
  filenamestr tmpch;
  //cpeds_queue<string> qstr;
  string tstr[10];
  cpeds_queue<long> qlong;
  
  //
  // initiate the object
  //
  MF = new mfsd_type;
  node = node_no;
  sprintf(tmpch,"mfsMFstat_node:%li",node);
  object_name=tmpch;
  mask = mask_name;
  simtype=sim_type;
  mmsf=mmf;
  thnum=th_num;

  // resolve the simulation type
  sim_type_num=0;
  if (Mscs_contains_str(sim_type,"Q",1))  { tstr[sim_type_num]="Q";  sim_type_num++;  }
  if (Mscs_contains_str(sim_type,"V",1))  { tstr[sim_type_num]="V";  sim_type_num++;  }
  if (Mscs_contains_str(sim_type,"W",1))  { tstr[sim_type_num]="W";  sim_type_num++;  }
  if (Mscs_contains_str(sim_type,"VW",2)) { tstr[sim_type_num]="VW"; sim_type_num++;  }

  if (sim_type_num <=0) { printf("|%s> * nothing to do; sim_type_num: %li\n",object_name.c_str(),sim_type_num); exit(0); } 
  else { printf("|%s> * initating the data level; sim_type_num: %li\n",object_name.c_str(),sim_type_num); } 
    
  //
  // initiate the data level
  //
  MF->data_num = sim_type_num;
  MF->dstr = new string[sim_type_num]; 
  MF->mfd = new mfssm_type[sim_type_num]; 
  for (d=0;d<sim_type_num;d++) { 
    MF->dstr[d] = tstr[d]; 
    //MF->mfd[d] = NULL; 
  }


  // check the simulations number
  simst=sim_stnum; 
  simen=sim_ennum;
  sim_num=sim_ennum-sim_stnum+1;
  if (sim_num <=0) { printf("|%s> * nothing to do; sim_num: %li\n",object_name.c_str(),sim_num); exit(0); }
  else { printf("|%s> * will process %li simulations\n",object_name.c_str(),sim_num); } 



  // check the smoothing 
  if (sm_num <= 0) { printf("|%s> * nothing to do; sm_num: %li\n",object_name.c_str(),sm_num); exit(0); }
  else { printf("|%s> * initating the smoothing level; sm_num: %li\n",object_name.c_str(),sm_num); }
 
  //
  // initiate the smoothing level
  //
  for (d=0;d<sim_type_num;d++) { 
    MF->mfd[d].sm_num = sm_num;
    MF->mfd[d].sm = new double[sm_num];
    MF->mfd[d].mfs = new mfsmm_type[sm_num];
    for (l=0;l<sm_num;l++) { 
      MF->mfd[d].sm[l] = sm[l];
      MF->mfd[d].mfs[l].mfm = NULL;
    }
  }
  
  // checking the number and types of the multimasks
  if (mm_num <= 0) { printf("|%s> * nothing to do; mm_num: %li\n",object_name.c_str(),mm_num); exit(0); }
  else { printf("|%s> * initating the multi-mask level; mm_num: %li\n",object_name.c_str(),mm_num); }

  for (m=0;m<mm_num;m++) { 
    if (Mscs_contains_str(mmf[m],"HP_2",4)) qlong.addq(48); else
      if (Mscs_contains_str(mmf[m],"HP_4",4)) qlong.addq(192); else
	if (Mscs_contains_str(mmf[m],"HP_8",4)) qlong.addq(768); 
	else { printf("ERROR: some of the multimasks defined is not either HP_2, nor HP_4, nor HP_8, check this out\n\n"); exit(0); }
  }
	  

  //
  // initiate the multi-mask  level
  //
  for (d=0;d<sim_type_num;d++) { 
    for (l=0;l<sm_num;l++) { 
      MF->mfd[d].mfs[l].mm_num = mm_num;
      MF->mfd[d].mfs[l].mfm = new mfsreg_type[mm_num];
      for (m=0;m<mm_num;m++) { 
	MF->mfd[d].mfs[l].mfm[m].mfr = NULL;
 	MF->mfd[d].mfs[l].mfm[m].reg_num = qlong(m);
	MF->mfd[d].mfs[l].mfm[m].Cm = NULL;
	MF->mfd[d].mfs[l].mfm[m].Cs = NULL;
	MF->mfd[d].mfs[l].mfm[m].CS = NULL;
	MF->mfd[d].mfs[l].mfm[m].CK = NULL;
      }
    }
  }
  qlong.delete_all_queue();


  //
  // initiate the regions level
  //
  for (d=0;d<sim_type_num;d++) { 
    for (l=0;l<sm_num;l++) { 
      for (m=0;m<mm_num;m++) { 
	MF->mfd[d].mfs[l].mfm[m].mfr = new mfs_type;
	reg_num = MF->mfd[d].mfs[l].mfm[m].reg_num;
	for (k=0;k<reg_num;k++) { 
	  MF->mfd[d].mfs[l].mfm[m].mfr[k].th = NULL;
	  MF->mfd[d].mfs[l].mfm[m].mfr[k].v0 = NULL;
	  MF->mfd[d].mfs[l].mfm[m].mfr[k].v1 = NULL;
	  MF->mfd[d].mfs[l].mfm[m].mfr[k].v2 = NULL;
	  MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv0 = NULL;
	  MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv1 = NULL;
	  MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv2 = NULL;
	  MF->mfd[d].mfs[l].mfm[m].mfr[k].regid = k;
	  MF->mfd[d].mfs[l].mfm[m].mfr[k].size = 0;
	}
      }
    }
  }

  // define the covariance matrix size
  cov_size = (th_num+1)*th_num/2; 

  //
  // initiate the thresholds and covariance matrices level
  //
  for (d=0;d<sim_type_num;d++) { 
    for (l=0;l<sm_num;l++) { 
      for (m=0;m<mm_num;m++) { 
	reg_num = MF->mfd[d].mfs[l].mfm[m].reg_num;
	for (k=0;k<reg_num;k++) { 
	  MF->mfd[d].mfs[l].mfm[m].mfr[k].th = new double[th_num];
	  MF->mfd[d].mfs[l].mfm[m].mfr[k].v0 = new double[th_num];
	  MF->mfd[d].mfs[l].mfm[m].mfr[k].v1 = new double[th_num];
	  MF->mfd[d].mfs[l].mfm[m].mfr[k].v2 = new double[th_num];
	  MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv0 = new float[cov_size];
	  MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv1 = new float[cov_size];
	  MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv2 = new float[cov_size];
	  for (nu=0;nu<th_num;nu++) { 
	    MF->mfd[d].mfs[l].mfm[m].mfr[k].th[nu] = 0;
	    MF->mfd[d].mfs[l].mfm[m].mfr[k].v0[nu] = 0;
	    MF->mfd[d].mfs[l].mfm[m].mfr[k].v1[nu] = 0;
	    MF->mfd[d].mfs[l].mfm[m].mfr[k].v2[nu] = 0;
	  }
	  for (nu=0;nu<cov_size;nu++) { 
	    MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv0[nu] = 0;
	    MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv1[nu] = 0;
	    MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv2[nu] = 0;
	  }
	}
      }
    }
  }

}

/****************************************************************************************************/
minkowski_fs_MFstat::~minkowski_fs_MFstat() {
  long d,l,m,k;
  delete [] MF->dstr;

  for (d=0;d<MF->data_num;d++) { 
    delete [] MF->mfd[d].sm;
    for (l=0;l<MF->mfd[d].sm_num;l++) { 
      for (m=0;m<MF->mfd[d].mfs[l].mm_num;m++) { 	
	for (k=0;k<MF->mfd[d].mfs[l].mfm[m].reg_num;k++) { 
	  delete [] MF->mfd[d].mfs[l].mfm[m].mfr[k].th;
	  delete [] MF->mfd[d].mfs[l].mfm[m].mfr[k].v0;
	  delete [] MF->mfd[d].mfs[l].mfm[m].mfr[k].v1;
	  delete [] MF->mfd[d].mfs[l].mfm[m].mfr[k].v2;
	  delete [] MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv0;
	  delete [] MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv1;
	  delete [] MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv2;
	}
      }
    }
  } 

  delete MF;
  MF=NULL;

}
/****************************************************************************************************/
// Routine that for a given number of strting simulation (sim_stnum) and ending simulation (sim_ennum) to be generated on fly,
// types of simulations to be generated (sim_type) = [KINC, QINC,VINC,WINC,allQVWINC, etc. see the comments of the 
// map_class::generate_gaussian_maps_WMAP_QVW routine for more info]. 
// (the size of the corresponding structure will be set automatically)
// for a given common number of thresholds at which functionals will be calculated (thnum)
// and for a given list of files to multimasks of size (mm_num) stored on the array (mmf)
// and for sm_num requested smoothing scales of the generated datasets, stored in sm double array 
// will calculate the quantity proportional to the
// covariance matrix of the MFs in each region independently for each simulation and dataset
// C = Sum_j (delta_nu1ijkl)(delta_nu2ijkl)
// delta_nu1/2ijkl = f_nu1/2ijkl - <f_nu1/2ijkl>

// this routine will start generating requested simulations one by one and calculating minkowski functionals 
// on the simulated maps. It will not store the simulations.

void minkowski_fs_MFstat::make_cov_and_mean_part(long node) {
  long j,l,m,d;
  map_class *sims,*sim,**derivs;
  minkowski_fs* fs;
  double **grads;
  bool calc_derivs;
  long seed_offset_size=8000000; // [s] - roughly 3 months of time
  
  //
  // initiate the simulation space
  //
  sim = new map_class("a sim");
  
  //
  // loop over the simulations
  //
  for (j=simst;j<=simen;j++) {
    // set seed offset different for different nodes to make sure there will be no same seed in the simulations
    // sim->set_seed_offset(-seed_offset_size*node); // commented out on 2009-01-22 -- moved to cpeds and individual applications
    // generate simulations
    sims=sim->generate_gaussian_maps_WMAP_QVW(j,"same_alms",PROGRAM_DIRS,"wmap5",mask.c_str(),0,0,simtype,"quiet!"); // the mask "mask" will be loaded in the "sim" object
    
    //
    // loop over the datasets
    //
    for (d=0;d<MF->data_num;d++) { 
      sim->import_map_data(sims[d],"T",1);
      sim->calculate_transformF(1,0,"",-1,0,"",""); // this will keep the alms of the map for further smoothing
      
      //
      // loop over the smoothings lengths
      //
      for (l=0;l<MF->mfd[d].sm_num;l++) {

	// smooth the simulation
	sim->calculate_inverse_transformF(1,(double)0,"",MF->mfd[d].sm[l],0,"","",0);
	sim->calculate_map_stats(0);
	sims[d].import_map_data((*sim),"T",1);
	calc_derivs=true;
	derivs = new map_class*[6];
	grads = new double*[3];

	//
	// loop over the multimasks ( for the moment these will be loaded each time a simulation is produced)
	//
	for (m=0;m<MF->mfd[d].mfs[l].mm_num;m++) { 
	  // load the multi-mask
	  sims[d].loadbinT(mmsf[m].c_str(),4);
	  if (calc_derivs) {
	    fs = sims[d].calculate_minkowski_v0v1v2(thnum,sim->minT,sim->maxT,true,true,derivs,grads);
	    calc_derivs=false; }
	  else 
	    fs = sims[d].calculate_minkowski_v0v1v2(thnum,sim->minT,sim->maxT,false,false,derivs,grads);

	  correlate_and_add_mfs(d,l,m,fs);
	  delete fs;
	}

	// you can save the thresholds used (grads[0]) here

	// remove the grads and derivs
	for (m=0;m<3;m++) { delete grads[m]; } delete [] grads;
	for (m=0;m<6;m++) { delete derivs[m]; } delete [] derivs;
      }
    }

    delete [] sims;
  }
  
  
}
/****************************************************************************************************/

void minkowski_fs_MFstat::make_cov_and_mean(long _simst, long _simen) {
  long data_num,sm_num,mm_num,reg_num,th_num;
  long d,l,k,m,idx,nu1,nu2;
  long Nsim = _simen-_simst+1;
  long Nsimleo=Nsim-1;

  data_num = get_data_num();
  th_num = get_thres_num();
  for (d=0;d<data_num;d++) {
    sm_num = get_sm_num(d);
    for (l=0;l<sm_num;l++) {
      mm_num = get_mm_num(d,l);
      for (m=0;m<mm_num;m++) {
	reg_num = get_reg_num(d,l,m);
	for (k=0;k<reg_num;k++) {
	  idx=0;
	  for (nu1=0;nu1<thnum;nu1++) { // for each threshold
	    for (nu2=0;nu2<=nu1;nu2++) { // subtract the means from the correlations for each two levels combinations
	      MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv0[idx] = ( MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv0[idx] - MF->mfd[d].mfs[l].mfm[m].mfr[k].v0[nu1] * MF->mfd[d].mfs[l].mfm[m].mfr[k].v0[nu2] / (double)Nsim)  / (double)Nsimleo;
	      MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv1[idx] = ( MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv1[idx] - MF->mfd[d].mfs[l].mfm[m].mfr[k].v1[nu1] * MF->mfd[d].mfs[l].mfm[m].mfr[k].v1[nu2] / (double)Nsim)  / (double)Nsimleo;
	      MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv2[idx] = ( MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv2[idx] - MF->mfd[d].mfs[l].mfm[m].mfr[k].v2[nu1] * MF->mfd[d].mfs[l].mfm[m].mfr[k].v2[nu2] / (double)Nsim)  / (double)Nsimleo;
	      idx++;
	    }
	    // derive the average MFs
	    MF->mfd[d].mfs[l].mfm[m].mfr[k].v0[nu1] /= (double)Nsim;
	    MF->mfd[d].mfs[l].mfm[m].mfr[k].v1[nu1] /= (double)Nsim;
	    MF->mfd[d].mfs[l].mfm[m].mfr[k].v2[nu1] /= (double)Nsim;
	  }
	}
      }
    }
  }

  simst = _simst;
  simen = _simen;
}
/****************************************************************************************************/

void minkowski_fs_MFstat::correlate_and_add_mfs(long d,long l,long m, minkowski_fs* fs) {
  long nu1,nu2;
  long k,idx;
  double mf0nu1,mf1nu1,mf2nu1;

  for (k=0;k<MF->mfd[d].mfs[l].mfm[m].reg_num;k++) { // for each region in mm
    idx=0;
    for (nu1=0;nu1<thnum;nu1++) { // for each threshold
      // sum the MFs (for average)
      mf0nu1 = fs->get_mf(k,0,nu1);
      mf1nu1 = fs->get_mf(k,1,nu1);
      mf2nu1 = fs->get_mf(k,2,nu1);
      MF->mfd[d].mfs[l].mfm[m].mfr[k].v0[nu1] += mf0nu1;
      MF->mfd[d].mfs[l].mfm[m].mfr[k].v1[nu1] += mf1nu1;
      MF->mfd[d].mfs[l].mfm[m].mfr[k].v2[nu1] += mf2nu1;
      for (nu2=0;nu2<=nu1;nu2++) { // correlate the MFs between thresholds fnu1*fnu2
	MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv0[idx] += mf0nu1 * fs->get_mf(k,0,nu2);
	MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv1[idx] += mf1nu1 * fs->get_mf(k,1,nu2);
	MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv2[idx] += mf2nu1 * fs->get_mf(k,2,nu2);
	idx++;
      }
    }
  }

}
/****************************************************************************************************/
void minkowski_fs_MFstat::invert_covariance_matrices() {
  long data_num,sm_num,mm_num,reg_num,th_num;
  long d,l,k,m,idx,nu1,nu2;
  matrix <double> C0(thnum,thnum);
  matrix <double> C1(thnum,thnum);
  matrix <double> C2(thnum,thnum);

  data_num = get_data_num();
  th_num = get_thres_num();
  for (d=0;d<data_num;d++) {
    sm_num = get_sm_num(d);
    for (l=0;l<sm_num;l++) {
      mm_num = get_mm_num(d,l);
      for (m=0;m<mm_num;m++) {
	reg_num = get_reg_num(d,l,m);
	for (k=0;k<reg_num;k++) {

	  // copy covs onto the matrix object
	  idx=0;
	  for (nu1=0;nu1<thnum;nu1++) { 
	    for (nu2=0;nu2<=nu1;nu2++) {
	      C0(nu1,nu2)=(double)(MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv0[idx]);
	      C0(nu2,nu1)=(double)(MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv0[idx]);
	      C1(nu1,nu2)=(double)(MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv1[idx]);
	      C1(nu2,nu1)=(double)(MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv1[idx]);
	      C2(nu1,nu2)=(double)(MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv2[idx]);
	      C2(nu2,nu1)=(double)(MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv2[idx]);
	      idx++;
	    }
	  }

	  // inverse covs
	  !C0; // here in general will have to be inserted the condition number check of some sort
	  !C1; // here in general will have to be inserted the condition number check of some sort
	  !C2; // here in general will have to be inserted the condition number check of some sort

	  // copy covs back from the matrix object
	  idx=0;
	  for (nu1=0;nu1<thnum;nu1++) { 
	    for (nu2=0;nu2<=nu1;nu2++) {
	      MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv0[idx]=(float)C0(nu1,nu2);
	      MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv1[idx]=(float)C1(nu1,nu2);
	      MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv2[idx]=(float)C2(nu1,nu2);
	      idx++;
	    }
	  }

	}
      }
    }
  }


}
/****************************************************************************************************/
// returns the number of sets of minkowski functionals stored in the object
/* long minkowski_fs_MFstat::get_size() { return fnum; } */

void minkowski_fs_MFstat::set_name(string name) { object_name=name; }
string minkowski_fs_MFstat::get_name() { return object_name; }

double* minkowski_fs_MFstat::get_mf(long d, long l, long m, long k, long mftype ) { 
  double * ret=NULL;
  if (mftype == 0) { ret = MF->mfd[d].mfs[l].mfm[m].mfr[k].v0; }
  if (mftype == 1) { ret = MF->mfd[d].mfs[l].mfm[m].mfr[k].v1; }
  if (mftype == 2) { ret = MF->mfd[d].mfs[l].mfm[m].mfr[k].v2; }
  if (mftype == 3) { ret = MF->mfd[d].mfs[l].mfm[m].mfr[k].th; }
  return ret;
}

float* minkowski_fs_MFstat::get_mfcov(long d, long l, long m, long k, long covtype ) { 
  float * ret=NULL;
  if (covtype == 0) { ret = MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv0; }
  if (covtype == 1) { ret = MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv1; }
  if (covtype == 2) { ret = MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv2; }
  return ret;
}

void minkowski_fs_MFstat::set_mf(long d, long l, long m, long k, long mftype, double* data ) {
  long nu;
  if (mftype == 0)  for (nu=0;nu<thnum;nu++) { MF->mfd[d].mfs[l].mfm[m].mfr[k].v0[nu]=data[nu]; }
  if (mftype == 1)  for (nu=0;nu<thnum;nu++) { MF->mfd[d].mfs[l].mfm[m].mfr[k].v1[nu]=data[nu]; }
  if (mftype == 2)  for (nu=0;nu<thnum;nu++) { MF->mfd[d].mfs[l].mfm[m].mfr[k].v2[nu]=data[nu]; }
  if (mftype == 3)  for (nu=0;nu<thnum;nu++) { MF->mfd[d].mfs[l].mfm[m].mfr[k].th[nu]=data[nu]; }
}

void minkowski_fs_MFstat::set_mfcov(long d, long l, long m, long k, long covtype, float* data ) {
  long covsize=get_cov_size();
  long nu;
  if (covtype == 0)  for (nu=0;nu<covsize;nu++) { MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv0[nu]=data[nu]; }
  if (covtype == 1)  for (nu=0;nu<covsize;nu++) { MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv1[nu]=data[nu]; }
  if (covtype == 2)  for (nu=0;nu<covsize;nu++) { MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv2[nu]=data[nu]; }
}

void minkowski_fs_MFstat::add_mf   (long d, long l, long m, long k, long mftype, double* data) {
  long nu;
  if (mftype == 0)  for (nu=0;nu<thnum;nu++) { MF->mfd[d].mfs[l].mfm[m].mfr[k].v0[nu]+=data[nu]; }
  if (mftype == 1)  for (nu=0;nu<thnum;nu++) { MF->mfd[d].mfs[l].mfm[m].mfr[k].v1[nu]+=data[nu]; }
  if (mftype == 2)  for (nu=0;nu<thnum;nu++) { MF->mfd[d].mfs[l].mfm[m].mfr[k].v2[nu]+=data[nu]; }
  if (mftype == 3)  for (nu=0;nu<thnum;nu++) { MF->mfd[d].mfs[l].mfm[m].mfr[k].th[nu]+=data[nu]; }
}

void minkowski_fs_MFstat::add_mfcov(long d, long l, long m, long k, long covtype, float*  data) {
  long covsize=get_cov_size();
  long nu;
  if (covtype == 0)  for (nu=0;nu<covsize;nu++) { MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv0[nu]+=data[nu]; }
  if (covtype == 1)  for (nu=0;nu<covsize;nu++) { MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv1[nu]+=data[nu]; }
  if (covtype == 2)  for (nu=0;nu<covsize;nu++) { MF->mfd[d].mfs[l].mfm[m].mfr[k].Cv2[nu]+=data[nu]; }
}


long minkowski_fs_MFstat::get_data_num() { return MF->data_num; }
long minkowski_fs_MFstat::get_sm_num(long d) { return MF->mfd[d].sm_num; }
long minkowski_fs_MFstat::get_mm_num(long d,long l) { return MF->mfd[d].mfs[l].mm_num; }
long minkowski_fs_MFstat::get_reg_num(long d,long l,long m) { return MF->mfd[d].mfs[l].mfm[m].reg_num; }
long minkowski_fs_MFstat::get_thres_num() { return thnum; }
long minkowski_fs_MFstat::get_reg_id(long d,long l,long m, long k) { return MF->mfd[d].mfs[l].mfm[m].mfr[k].regid; }
long minkowski_fs_MFstat::get_reg_size(long d,long l,long m, long k) { return MF->mfd[d].mfs[l].mfm[m].mfr[k].size; }
long minkowski_fs_MFstat::get_cov_size() { return (thnum+1)*thnum/2; }



// saves the functionals to a file
void minkowski_fs_MFstat::save(string name, long how) {
  long data_num,sm_num,mm_num,reg_num;
  long d,l,k,m;
  FILE * f;
  size_t sd=sizeof(double);
  size_t sf=sizeof(float);
  long covsize=get_cov_size();

  if (how == 1) {
    name+=".minkbin";
    printf("|%s> * saving the minkowski functionals to a binary file %s\n",object_name.c_str(),name.c_str());
    f = fopen(name.c_str(),"wb");

    // object name
    fprintf(f,"%s\n",object_name.c_str());
    // save simtype
    fprintf(f,"%s\n",simtype.c_str());
    // save mask name
    fprintf(f,"%s\n",mask.c_str());
    // sim start sim end
    fprintf(f,"%li %li\n",simst,simen);
    // thresholds number
    fprintf(f,"%li\n",thnum);
    // save data_num
    data_num = get_data_num();
    fprintf(f,"%li\n",data_num);
    // save dataset name
    for (d=0;d<data_num;d++) {
      fprintf(f,"%s\n",MF->dstr[d].c_str());

      sm_num = get_sm_num(d);
      fprintf(f,"%li\n",sm_num);
      for (l=0;l<sm_num;l++) {
	fprintf(f,"%lE\n",MF->mfd[d].sm[l]);

	mm_num = get_mm_num(d,l);
	fprintf(f,"%li\n",mm_num);
	for (m=0;m<mm_num;m++) {

	  reg_num = get_reg_num(d,l,m);
	  fprintf(f,"%li\n",reg_num);
	  for (k=0;k<reg_num;k++) {
	    fprintf(f,"%li %li\n",get_reg_id(d,l,m,k),get_reg_size(d,l,m,k));

	    fwrite(get_mf(d,l,m,k,3),sd,thnum,f);
	    fwrite(get_mf(d,l,m,k,0),sd,thnum,f);
    	    fwrite(get_mf(d,l,m,k,1),sd,thnum,f);
	    fwrite(get_mf(d,l,m,k,2),sd,thnum,f);
	    fwrite(get_mfcov(d,l,m,k,0),sf,covsize,f);
    	    fwrite(get_mfcov(d,l,m,k,1),sf,covsize,f);
	    fwrite(get_mfcov(d,l,m,k,2),sf,covsize,f);
	  }
	}
      }
    }
    fclose(f);
  }

}

/* void minkowski_fs_MFstat::normalize_V0s() { */
/*   long i; */
/*   for (i=0;i<fnum;i++) mfs[i].v0->normalize_V0(); */

/* } */

