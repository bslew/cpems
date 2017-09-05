// This class implements handling of radiative transfer functions for eg. NG map generation etc.
//#include <gsl/gsl_sf_bessel.h>
#include "Mscs-transfer_function.h"
#include "cpeds-templates.h"


//*********************************************************************************
transfer_function::transfer_function() { // minimal constructor 
  t.f=NULL;
  t.k=NULL;
  //lmax_all=0;
  t.kmin=t.kmax=0;
  t.k_num=0;
  tf_loaded=false;
  real_space=false;
  tf_load_file="";
  object_name="tf";
  printf("|%s> -- initiating empty transfer function\n", object_name.c_str());
}

//*********************************************************************************
transfer_function::transfer_function(long ksize, long lsize) { // constructor with size in l and k directions
  t.f=NULL;
  t.k=NULL;

  t.kmin=t.kmax=0;
  t.k_num=ksize;
  t.lmin=0; t.lmax=0;
  t.l_num=lsize;
  object_name="tf";

  printf("|%s> -- initiating transfer function of ksize: %li and lsize: %li\n", object_name.c_str(), ksize,lsize);
  initiate_tf(0,0,ksize,0,0,lsize);
}

//*********************************************************************************
transfer_function::transfer_function(transfer_function* src,string clonename) { // constructor for cloning tfs
  string srcname=src->get_name();

  t.f=NULL;
  t.k=NULL;
  t.kmin=t.kmax=0;
  t.k_num=0;
  tf_loaded=false;
  real_space=false;
  tf_load_file="";

  if (clonename=="") object_name="clone of "+src->get_name(); else object_name=clonename;
  printf("|%s> -- cloning transfer function from %s\n", object_name.c_str(), srcname.c_str());

  if (!src->is_loaded()) { 
    t.f=NULL;  t.k=NULL; 
    t.kmin=t.kmax=0;
    t.k_num=0;
    t.lmin=0; t.lmax=0;
    t.l_num=0;
  } 
  else {
    // initiate space
    //initiate_tf(1.0/src->get_kmax(),1.0/src->get_kmin(),src->get_size(),0,src->get_lmax(),src->get_lnum()); 
    initiate_tf(src->get_kmin(),src->get_kmax(),src->get_knum(),0,src->get_lmax(),src->get_l_num()); 

    // copy the content and set all the variables
    *this=*src;    
  }
    
}

//*********************************************************************************
void transfer_function::kill_k() {
  printf("|%s> -- deleting table for k\n", object_name.c_str());
  if (t.k !=NULL) { delete [] t.k; t.k=NULL; }
  tf_loaded=false;
  tf_load_file="";
}
//*********************************************************************************
void transfer_function::kill_f() {
  printf("|%s> -- deleting table for f\n", object_name.c_str());
  if (t.f !=NULL) { delete [] t.f; t.f=NULL; }
  tf_loaded=false;
  tf_load_file="";
}
//*********************************************************************************
transfer_function::~transfer_function() { 
  kill_k();
  kill_f();
}

//*********************************************************************************
void transfer_function::set_name(string name) { 
  printf("|%s> -- changing name of object to: %s\n", object_name.c_str(), name.c_str());
  object_name=name; 
}
//*********************************************************************************
string transfer_function::get_name() { return object_name; }

//*********************************************************************************
void transfer_function::initiate_tf(double kmin_loc, double kmax_loc, long k_num_loc, long lmin_loc, long lmax_loc, long l_num_loc) {
  long i;
  t.kmin=kmin_loc;  
  t.kmax=kmax_loc;
  t.k_num=k_num_loc;
  t.lmin=lmin_loc;
  t.lmax=lmax_loc;
  t.l_num=l_num_loc; 
  kill_f(); t.f=new power_spectrum[t.k_num]; 
  for (i=0;i<k_num_loc;i++) { t.f[i].set_quietness(0); t.f[i].initiate_power_spectrum(l_num_loc); }
  kill_k(); t.k=new double[t.k_num]; 
  clear_tf(); 
  tf_loaded=false;
  tf_load_file="";

}


//*********************************************************************************
// returns the tf for requested l and k; assumes that k vales are stored in increasing order
double transfer_function::get_tf(double k, long l, int *result) {
  long ik;
  ik=cpeds_find_value(k,t.k,t.k_num,0,t.k_num);
  if (result != NULL) *result=0; // exit status code; not used yet.
  return t.f[ik].get_Cl(l);
}

//*********************************************************************************
// returns the tf for requested l and k; assumes that k vales are stored in increasing order
double transfer_function::get_tf(long ik, long l, int *result) {
  if (result != NULL) *result=0; // exit status code; not used yet.
  return t.f[ik].get_Cl(l);
}
//*********************************************************************************
int transfer_function::set_tf(long ik, double k, long l, double tf) {
  t.k[ik]=k;
  t.f[ik].set_Cl(l,tf);
  return 0;
}

//*********************************************************************************
long transfer_function::get_lmax(double k) { 
  long ik;
  ik=cpeds_find_value(k,t.k,t.k_num,0,t.k_num);
  return t.f[ik].get_lmax();
}
//*********************************************************************************
long transfer_function::get_lmax() { 
  return t.lmax;
}

//*********************************************************************************
long transfer_function::get_lsize(long ik) { 
  return t.f[ik].get_size();
}

//*********************************************************************************
// to be implemented when needed
long transfer_function::get_lsize(double k) { 
  return 0;
}
//*********************************************************************************
// returns the maximal size of the tf. in l direction l_num=lmax+1 to include l=0
// this might be different from the stored value of t.l_num since t.l_num keeps the number of ls that were eg. loaded from file with useful info. other l < lmin were set to zero by default. this is theoretically stupid since the C_l object is capable of storing also l not only Cl and hence the arrays might start filling from index 0 yet keeping the information on higher ls.
long transfer_function::get_l_num() { 
  return t.lmax+1;
}

//*********************************************************************************
// NOT IMPLEMENTED YET; DEFINE HOW IT IS TO BE IMPLEMENTED...
long transfer_function::set_lmax(double k) {
  return 0;
}

//*********************************************************************************
// NOT IMPLEMENTED YET; DEFINE HOW IT IS TO BE IMPLEMENTED...
long transfer_function::set_all_lmax() {
  return 0;
}

//*********************************************************************************
double transfer_function::get_kmax() { return t.kmax; }

//*********************************************************************************
/* double transfer_function::set_kmax() { t.kmax= */
/* } */

//*********************************************************************************
double transfer_function::get_kmin() { return t.kmin; }
//*********************************************************************************
long transfer_function::get_knum() { return t.k_num; }

//*********************************************************************************
// returns the address of an array with ks values of the transfer function for various purposes
double* transfer_function::get_kstab() { 
  double * tmp = new double[t.k_num];
  long ik;
  for (ik=0;ik<t.k_num;ik++) { tmp[ik]=t.k[ik]; /* printf("%lE\n",tmp[ik]); */ }
  return tmp; 
}
//*********************************************************************************
// returns the address of an array with 2/PI*k^2 values of the transfer function for various purposes
double* transfer_function::get_ks2tab() { 
  double * tmp = new double[t.k_num];
  long ik;
  for (ik=0;ik<t.k_num;ik++) { tmp[ik]=twooverPI*t.k[ik]*t.k[ik]; /* printf("%lE\n",tmp[ik]); */ }
/*   for (ik=0;ik<t.k_num;ik++) { tmp[ik]=4*PI*t.k[ik]*t.k[ik]; /\* printf("%lE\n",tmp[ik]); *\/ } */
  return tmp; 
}

//*********************************************************************************
// -1 val means an error
double transfer_function::get_k(long ik) { 
  if ( ik >=0 && ik < t.k_num) return t.k[ik]; 
  return -1;
}
//*********************************************************************************
void transfer_function::set_k(long ik, double k) { 
  if ( ik >=0 && ik < t.k_num) t.k[ik]=k; 
}

//*********************************************************************************
long transfer_function::get_ikmin_non_zero(long l) { 
  long ik,r=0;
  
  for (ik=0;ik<t.k_num;ik++) { r=ik;
    if (get_tf(ik,l,NULL) != 0) ik=t.k_num; }
  return r;
}
//*********************************************************************************
long transfer_function::get_ikmax_non_zero(long l) { 
  long ik,ikst=t.k_num-1,r=ikst;
  
  for (ik=ikst;ik>=0;ik--) { r=ik;
    if (get_tf(ik,l,NULL) != 0) ik=-1; }
  return r;
}

//*********************************************************************************
// returns the addres to the power_spectrum object containing a transfer function 
//for a mode ik
power_spectrum* transfer_function::get_tf_ik(long ik) {
  return &t.f[ik];
}

//*********************************************************************************
/* double transfer_function::set_kmin() { */
/* } */
//*********************************************************************************
void transfer_function::clear_tf() {
  long i;

  printf("|%s> -- zeroing transfer function\n", object_name.c_str());
  for (i=0;i<t.k_num;i++) {
    t.f[i].reset_Cl();
    t.k[i]=0;
  }

  tf_load_file="";
}

//*********************************************************************************
int transfer_function::load_tf(string fname) {
  FILE *f,*finfo;

  long il,ik,tff_cols_num;
  long l,lmin,lmax,k_num,l_num;
  double k,kprev,kmin,kmax,tf;  
  string info_file; info_file=fname+".info";
  cpeds_queue <double> ksl; // cpeds_queue list of ks in the transfer function file
  cpeds_queue <long> lmaxl; // cpeds_queue of lmax values

  tff_cols_num = cpeds_get_cols_num_first_ln(fname.c_str(),NULL);

  printf("|%s> * loading the transfer function from file: %s with %li columns\n",object_name.c_str(),fname.c_str(),tff_cols_num);
  if (tf_loaded) { if (tf_load_file==fname) {
      printf("|%s> NOTE: -- transfer function already loaded from this file: %s, not loading; to force loading kill the tf. function first\n",object_name.c_str(),fname.c_str());
      return 0;
    } 
    else {
      printf("|%s> NOTE: -- some other transfer function is already loaded from file: %s; killing the old tf. first\n",object_name.c_str(),fname.c_str());
      kill_k(); kill_f();
    }
  }

  f = fopen(fname.c_str(),"r");
  if (f==NULL) { return -1; } else 

  finfo = fopen(info_file.c_str(),"r");
  if (finfo==NULL) {  
  

    // scanning the file structure to initiate the tf structures
    printf("|%s>  -- scanning the file structure\n",object_name.c_str());
    if  (!feof(f)) {
      if (tff_cols_num == 3) fscanf(f,"%lE %li %lE",&k,&l,&tf); // for the T only tf
      if (tff_cols_num == 4) fscanf(f,"%lE %li %lE %*E",&k,&l,&tf); // for the T E tf
      /*     printf("%lE %li %lE\n",k,l,tf); */
      lmin=lmax=l;
      kmin=kmax=kprev=k; ksl.addq(k);
      k_num=1;
    }
    else { printf("|%s>  -- the tf file: %s is empty. tf not loaded\n",object_name.c_str(),fname.c_str()); return -1; }
    
    while (!feof(f)) {
      if (tff_cols_num == 3) fscanf(f,"%lE %li %lE",&k,&l,&tf); // for the T only tf
      if (tff_cols_num == 4) fscanf(f,"%lE %li %lE %*E",&k,&l,&tf); // for the T E tf
      /*     printf("%lE %li %lE\n",k,l,tf); */
      if (k!=kprev) { ksl.addq(k); k_num++; kprev=k; printf("|%s>  -- number of ks: %li\r",object_name.c_str(),k_num); }
      if (l < lmin) lmin=l;
      if (l > lmax) lmax=l;
      if (k < kmin) kmin=k;
      if (k > kmax) kmax=k;
    };
    fclose(f); f=NULL;
    printf("\n");
    l_num=lmax-lmin+1;
    printf("|%s>  -- scanning the file structure done.\n",object_name.c_str());
    
    
    printf("|%s>  -- saving info file and loading\n",object_name.c_str());
    finfo = fopen(info_file.c_str(),"w");
    fprintf(finfo,"kmin: %lE\n",kmin);
    fprintf(finfo,"kmax: %lE\n",kmax);
    fprintf(finfo,"k_num: %li\n",k_num);
    fprintf(finfo,"lmin: %li\n",lmin);
    fprintf(finfo,"lmax: %li\n",lmax);
    fprintf(finfo,"l_num: %li\n\nks tab: k lmax\n",l_num);
    for (ik=0;ik<k_num;ik++) { 
/*       fprintf(finfo,"%lE %li\n",ksl(ik),(long)0);  */
      fprintf(finfo,"%lE %li\n",0.0,(long)0); 
    }
    fclose(finfo);
  } 
  else {
    printf("|%s>  -- reading info file and loading\n",object_name.c_str());
    finfo = fopen(info_file.c_str(),"r");
    fscanf(finfo,"%*s %lE",&kmin);
    fscanf(finfo,"%*s %lE",&kmax);
    fscanf(finfo,"%*s %li",&k_num);
    fscanf(finfo,"%*s %li",&lmin);
    fscanf(finfo,"%*s %li",&lmax);
    fscanf(finfo,"%*s %li",&l_num);
    fclose(finfo);
  }



  // print file info

  printf("|%s>  -- tf file %s: lmin: %li lmax:%li l_num: %li   kmin: %lE kmax:%lE k_num: %li \n",object_name.c_str(),fname.c_str(),lmin,lmax,l_num,kmin,kmax,k_num);

  initiate_tf(kmin,kmax,k_num,0,lmax,l_num+lmin);


  printf("|%s>  -- tf file %s: lmin: %li lmax:%li l_num: %li   kmin: %lE kmax:%lE k_num: %li \n",object_name.c_str(),fname.c_str(),t.lmin,t.lmax,t.l_num,t.kmin,t.kmax,t.k_num);

  // load the tf 

  if (f==NULL) f = fopen(fname.c_str(),"r");


  for (ik=0;ik<t.k_num;ik++) {
    //printf("|%s>  -- loading k NO %li of k_num: %li\r",object_name.c_str(),ik,k_num);
    for (il=lmin;il<=lmax;il++) { 
      if (feof(f)!=0) {printf("ERROR: check the loading tf: eof occured"); exit(0); }
      if (tff_cols_num == 3) fscanf(f,"%lE %li %lE",&k,&l,&tf); // for the T only tf
      if (tff_cols_num == 4) fscanf(f,"%lE %li %lE %*E",&k,&l,&tf); // for the T E tf

      set_tf(ik,k,l,tf);
      //printf("%li %lE %li %lE\n",ik,k,l,get_tf(ik,l,NULL));
    }
  }
  printf("\n");
  fclose(f);

  t.lmin=lmin; t.lmax=lmax; t.l_num=l_num; t.kmin=kmin; t.kmax=kmax; t.k_num=k_num;
  tf_loaded=true;
  tf_load_file=fname;

  shrink_in_lmax();

  printf("|%s> * loading the transfer function done\n",object_name.c_str());
  return 0;
}
//*********************************************************************************
bool transfer_function::is_loaded() {
  return tf_loaded;
}
//*********************************************************************************
int transfer_function::save_tf(string fname, string how) {
  long l,ik;
  FILE * f;
  filenamestr tmpch;

  if (how == "matrix") {
    printf("|%s> * saving the transfer function from file: %s in a %s format\n",object_name.c_str(),fname.c_str(),how.c_str());

    sprintf(tmpch,"%s.mat",fname.c_str());
    f=fopen(tmpch,"w");
    for (l=0;l<=t.lmax;l++) {
      for (ik=0;ik<t.k_num;ik++) { fprintf(f,"%lE ",t.f[ik].get_Cl(l)); }
      fprintf(f,"\n");
    }
    fclose(f);
  }

  if (how == "matrixl") {
    printf("|%s> * saving the transfer function from file: %s in a %s format\n",object_name.c_str(),fname.c_str(),how.c_str());

    sprintf(tmpch,"%s.mat",fname.c_str());
    f=fopen(tmpch,"w");
    for (l=0;l<=t.lmax;l++) {
      fprintf(f,"%li ",(long)t.f[0].get_l(l)); 
      for (ik=0;ik<t.k_num;ik++) { fprintf(f,"%lE ",t.f[ik].get_Cl(l)); }
      fprintf(f,"\n");
    }
    fclose(f);
  }

  if (how == "tf") {
    printf("|%s> * saving the transfer function from file: %s in a %s format\n",object_name.c_str(),fname.c_str(),how.c_str());

    sprintf(tmpch,"%s",fname.c_str());
    f=fopen(tmpch,"w");
    for (ik=0;ik<t.k_num;ik++) { 
      for (l=0;l<=t.lmax;l++) { fprintf(f,"%lE %li %lE\n",t.k[ik],(long)t.f[0].get_l(l),t.f[ik].get_Cl(l)); }
    }
    fclose(f);
  }

  sprintf(tmpch,"%s.mat.info",fname.c_str());
  f=fopen(tmpch,"w");
  fprintf(f,"kmin: %lE\n",t.kmin);
  fprintf(f,"kmax: %lE\n",t.kmax);
  fprintf(f,"k_num: %li\n",t.k_num);
  fprintf(f,"lmin: %li\n",t.lmin);
  fprintf(f,"lmax: %li\n",t.lmax);
  fprintf(f,"l_num: %li\n\nks tab: k lmax\n",t.l_num);
  for (ik=0;ik<t.k_num;ik++) { fprintf(f,"%lE %li\n",get_k(ik),get_lsize(ik)-1); }
  fclose(f);

  
  return 0;
}
//*********************************************************************************
// imports al power spectrum to i'th position on it's own structure
void transfer_function::import_ithCl(power_spectrum *al, long i) {
  string tmp =al->get_name();
  printf("|%s> * importing the part of transfer function in l direction at ik position %li from: %s\n",object_name.c_str(),i,tmp.c_str());
  t.f[i] = *al;
}
//*********************************************************************************
void transfer_function::update_self() {
  long ik;

  if (t.f==NULL) { t.kmin=t.kmax=t.lmin=t.lmax=0; return; }

  // update kmin and kmax from tab k and lmin and lmax from C_l structures
  t.kmin=t.k[0]; t.kmax=t.k[0];
  t.lmin=t.lmax=t.f[0].get_lmin();

  for (ik=0;ik<t.k_num;ik++) { 
    if (t.kmin > t.k[ik]) t.kmin=t.k[ik];
    if (t.kmax < t.k[ik]) t.kmax=t.k[ik];
    if (t.lmin > t.f[ik].get_lmin()) t.lmin = t.f[ik].get_lmin(); // CAVEAT !! make sure you know how special cases are treated: eg. when there are gaps in the initiation of the C_l in the f table then eg. lmin can be returned as -1 which is a control values innitiated in the power_spectrum object
    if (t.lmax < t.f[ik].get_lmax()) t.lmax = t.f[ik].get_lmax();
  }

  t.l_num=t.lmax-t.lmin+1;
  printf("|%s> -- updating: transfer function: k_num %li kmin: %lE max: %lE l_num: %li lmin: %li lmax: %li \n", object_name.c_str(), t.k_num, t.kmin, t.kmax, t.l_num, t.lmin,t.lmax );

}
//*********************************************************************************
// sorts the ks table and t.f table accordingly
// NOT IMPLEMENTED YET
void transfer_function::sort_in_ks(long ord) {
  long ik;

  for (ik=0;ik<t.k_num;ik++) {
    

  }
  
}
//*********************************************************************************
// compute trhe real transfer function for a given comoving distance r
// Delta_l(r) = 2/pi int k^2 dk Delta_l(k) j_l(kr)
// distance should be given in Mpc (to match the transfer function converntions of the cmbfast
// can be given in H^-1

power_spectrum* transfer_function::get_real_tf(double r, double * k, double *kk) {
  long ik,l;
  long lnum;
  power_spectrum * Deltar;
  double *fn = new double[t.k_num];
/*   transfer_function  Dlr(this,"Xlr"); // clone this transfer function */

  long ikminnz,ikmaxnz;

  //x=get_kstab();
  lnum=get_lmax(); //printf("lnum %li\n",lnum); 
  Deltar = new power_spectrum(lnum,0);

  // precompute the 2/PI * k^2 values
  //for (ik=0;ik<knum;ik++) {     kk[ik]= twooverPI * x[ik]*x[ik];  } // remove this outside this function or implement calculations for range of r values here without many calls to this routine
/* printf("ik %li k %lE k^2 %lE\n",ik,x[ik],kk[ik]); */
/*   exit(0); */

  for (l=0;l<lnum;l++) {
/*   for (l=150;l<151;l++) { */
    
/*     EXPERIMENTALL BEG */
    //printf("l %li\n",l);

/*     EXPERIMENTALL END */

    // define the k integration limits
    ikminnz=get_ikmin_non_zero(l);
    ikmaxnz=get_ikmax_non_zero(l);

/*     EXPERIMENTALL BEG */
    //printf("ikminmax %li %li\n",ikminnz,ikmaxnz);
/*     EXPERIMENTALL END */
    

    // integrating over k
    for (ik=ikminnz;ik<=ikmaxnz;ik++) {
/*       fn[ik] = kk[ik] * t.f[ik].get_Cl(l) * cpeds_spherical_bessel_fn(l,k[ik] * r); */
      fn[ik] = t.f[ik].get_Cl(l) * cpeds_spherical_bessel_fn(l,k[ik] * r) / k[ik] *4*PI;
/*     EXPERIMENTALL BEG */
      //printf("tf %lE jl: %lE\n", t.f[ik].get_Cl(l), cpeds_spherical_bessel_fn(l,x[ik] * r));
/*     EXPERIMENTALL END */

    }
    
    Deltar->set_Cl(l, cpeds_integrate_1d(t.k_num,k,fn,ikminnz,ikmaxnz));
  }
  //  exit(0);
/*   *(get_tf_ik(knum-ir-1))=*Deltar; */
/*     set_k(knum-ir-1, 1.0/Dlr.get_k(ir)); */

  //  sort_in_ks(12);
/*   update_self(); */

  delete [] fn;
  //  delete [] kk;
  //delete [] x;
/*   Dlr.kill_k(); Dlr.kill_f(); */
  return Deltar;
}

//*********************************************************************************
// builds a real transfer function object containinig the real transfer functions in selected 
// range of the distances sampled with  Nsteps transfer functions

// the alghorithm for testing how the resulting real space transfer functions really transfer should be applied
// so that given sufficiently large range and fine sampling the resulting transfer function object won't grow too much
// since in most of the regions the tfs are zero and there's no point in keeping them.

// So the alghorithm is: THIS IS NOT IMPLEMENTED YET -- perhaps it will not be needed
// * get range of distances and suggested sampling density in number of steps
// * obtain each transfer functions for each distance value
// * decide whether the transfer function is usefull or not (by testing it's valies vs some critiria)
// * drop it if it's useless or store it in the object if it is useful.
// 
// the sampling of the real space should be no worse than the thickness of the SLS which is ~ 71 Mpc comoving (for redshifts 1100 to 800)
// 
transfer_function* transfer_function::get_real_space_transfer_functions(double rmin,double rmax, long Nsteps) {
  double dr,r;
  long ir,lmax=get_lmax();
  transfer_function* tfr = new transfer_function(Nsteps,lmax);
  double *k = get_kstab();
  double *kk = get_ks2tab();

  // step with which sample the tf
  dr=(rmax-rmin)/(double)Nsteps;

  r=rmin;
  

  for (ir=0;ir<Nsteps;ir++) {
    tfr->import_ithCl(get_real_tf(r,k,kk),ir);
    tfr->set_k(ir,r);
    r+=dr;
  }
  
  delete [] k;
  delete [] kk;
  return tfr;
}
//*********************************************************************************
bool transfer_function::real_space_tf() {
  return real_space;
}
//*********************************************************************************
// checks for larges nonzero l in each k in tf and rescales the tf accordingly
void transfer_function::shrink_in_lmax() {
  long ik;

  printf("|%s> * shrinking transfer function in lmax\n",object_name.c_str());

  for (ik=0;ik<t.k_num;ik++) {     t.f[ik].shrink_in_lmax(1);  }  
  update_self();
}
//*********************************************************************************
// these routines calculate the C_l = 4pi \int Delta^2_l(k) * P_R(k) dk/k  
// which is the power spectrum of the initial curvature perturbation
// P_R(k) = A(k/k0)^(ns-1)

/* power_spectrum* transfer_function::derive_Cl(double A, double ns, double k0) { // derives the C_l given power-law primordial power spectrum with given A and ns */
/*   power_spectrum* Cl = new power_spectrum(t.lmax+1); */
/*   long l,ik,knum=t.k_num-1; */
/*   double tf,dk; */

/*   printf("|%s> * deriving the CMB angular power spectrum from transfer function and given model P(k)\n",object_name.c_str()); */

/*   for (l=t.lmin;l<=t.lmax;l++) { */
/* /\*     printf("integrating l=%li, lmax=%li\r",l,t.lmax); *\/ */
/*     for (ik=0;ik<knum;ik++) { */
/*       dk=t.k[ik+1]-t.k[ik];  tf=t.f[ik].get_Cl(l);  */
/*       Cl->set_Cl(l,Cl->get_Cl(l) + tf*tf * dk/t.k[ik] * pow(t.k[ik]/k0,ns-1.0)  ); */
/* /\*       Cl->set_Cl(l,Cl->get_Cl(l) + tf*tf * t.k[ik]*t.k[ik]*dk * pow(t.k[ik]/k0,ns-1.0)  ); *\/ */
/*     } */
/* /\*     Cl->set_Cl(l,Cl->get_Cl(l) * 4*PI * A ); *\/ */
/* /\*     Cl->set_Cl(l,Cl->get_Cl(l) * 2.0/PI * A ); *\/ */
/*   } */
/*   (*Cl)*= 4*PI*A; // normalization */

/*   return Cl; */
/* } */

//*********************************************************************************
// these routines calculate the C_l = 4pi \int Delta^2_l(k) * Delta_zeta(k) dk/k  
// where Delta_zeta(k) = k^3/2pi^2 P_zeta(k) -- power per log interval of the comoving curvatrure pert.
// where P_zeta(k) = 2pi^2/k^3 (k/k0)^(ns-1) - primordial power spectrum of the comoving curvatrure pert.
// the unnorm option of the new cmbfast calibrates Delta_zeta(k) = 1 for scale invariant power spectrum ns=1
// and normalized to unity at pivot point k0 for ns!=1
// which is the power spectrum of the initial curvature perturbation
// P_R(k) = A(k/k0)^(ns-1)

power_spectrum* transfer_function::derive_Cl(matter_power_spectrum* Pk) { // derives the C_l given power-law primordial power spectrum with given P(k) in object matter_power_spectrum 

  power_spectrum* Cl = new power_spectrum(t.lmax+1);
  long l,ik;
  double tf,*fn=new double[t.k_num],*x;
  long ikminnz,ikmaxnz;

  printf("|%s> * deriving the CMB angular power spectrum from transfer function and given P(k)\n",object_name.c_str());
  x=get_kstab();

  for (l=t.lmin;l<=t.lmax;l++) {    

    // define the k integration limits
    ikminnz=get_ikmin_non_zero(l);
    ikmaxnz=get_ikmax_non_zero(l);

    // integrating over k
    for (ik=ikminnz;ik<=ikmaxnz;ik++) {
/*     for (ik=0;ik<t.k_num;ik++) { */
      tf=t.f[ik].get_Cl(l); 
      fn[ik]= tf*tf * Pk->get_P(x[ik]) / x[ik];
/*       Cl->set_Cl(l,Cl->get_Cl(l) + tf*tf * dk * Pk->get_P(t.k[ik])/t.k[ik] ); */
    }
    Cl->set_Cl(l, cpeds_integrate_1d(t.k_num,x,fn,ikminnz,ikmaxnz));
  }

  (*Cl)*=4*PI;

  delete [] fn;
  delete [] x;
  return Cl;
}

//*********************************************************************************
// derives the alpha_l(r) or beta_l(r) as defined in astro-ph/0305189 depending on what goes into P(k)
// alpha_l(r) = 2/PI int k^2 dk      T_l(k) jl(kr)
// beta_l(r)  = 2/PI int k^2 dk P(k) T_l(k) jl(kr)
// 
power_spectrum* transfer_function::integrate_alphabeta_l(double r, matter_power_spectrum* Pk) { 
  power_spectrum *Xlr = new power_spectrum(get_l_num()); Xlr->set_name("Xl");
  long ik,l,lmaxloc=get_lmax();

  double *fn = new double[t.k_num],*k,*kk;
  string tmp=Pk->get_name();
  //  FILE *f1,*f;
  long ikminnz,ikmaxnz;

  printf("|%s> * deriving alpha_l(r) or beta_l(r) for r=%lE given transfer function and P(k) given in %s\n",object_name.c_str(), r, tmp.c_str());


  k=get_kstab();
  kk=get_ks2tab();

  // precompute the 2/PI* k^2 values
/*   for (ik=0;ik<t.k_num;ik++) {  */
/*     fn[ik]=t.k[ik]*t.k[ik]*twooverPI; */
/*     //printf("k=%lE kr=%lE\n",t.k[ik],t.k[ik]*r); */
/*   } */

/*   f=fopen("/home/blew/programy/cmbfast/jlkrs.conv","r"); */


/* EXPERIMENTALL -- DUMPING BESSELS FUNCTIONS TO FILE FOR TESTS  BEGIN*/
/*   jlkr = new double[1501]; */
/*   f1 = fopen("bessel-sph-CDM.dump","w");  */
/*   for (ik=0;ik<t.k_num;ik++) { */
/*     gsl_sf_bessel_jl_steed_array((int)1500,t.k[ik]*r,jlkr); // these are the spherical Bessel functions: this is the good one */
/*     for (l=0;l<=1500;l++) { fprintf(f1,"%lE\n",jlkr[l]); printf("%lE %li %lE\n",t.k[ik]*r,l,jlkr[l]); }  */
/*   } */
/*   delete [] jlkr; */
/*   fclose(f1); // dump bessels */
/*   exit(0); */

/* ------------------------------------------------------------------ */
/*  NOTE ON THE ACCURACY OF THE SPHERICAL BESSEL FUNCTION DERIVIATION */
/* ------------------------------------------------------------------ */
/* Comparizon of the bessel functions derived from the jlgen of the cmbfast and the gsl */
/* gsl_sf_bessel_jl gives ~ 14% of results inconsistent @ level > 1% */
/* for CDM model run where kr ranges from 1.617783E-01 to 2.934378E+03 from lmin=0 to  lmax=1500 */
/* Of those the top 1000 most inconsistent were checked against results from Mathematica */
/* and definitelly the jlgen alghoritm returns far better results than the gsl routines. */
/* The relative errors reach vals as large as 82000 ! */
/* Most of the larges rel errros is concentrated around l=100 and k>1000 */
/* for k<1000 the consistency is @ level ~1 */
/* the strongest rel errors are concentrated at |jl| values of ~ (1e-4,1e-3) */
/* see plots in the cmbfast directory */

/* HENCE THOSE OF THE JLGEN OF THE CMBFAST WILL BE USED IN THIS PROGRAM */

/* EXPERIMENTALL -- DUMPING BESSELS FUNCTIONS TO FILE  END*/


  // integrating over k
/*   for (ik=0;ik<k_numleo;ik++) { */
/*   for (ik=0;ik<t.k_num;ik++) { */

/*     lmaxloc=t.f[ik].get_lmax(); */

/*     // precompute the spherical bessel functions */
/*     jlkr = new double[lmaxloc+1]; */
/*     for (l=0;l<=lmaxloc;l++) {  */
/*       if (go) jl= gsl_sf_bessel_jl(l,t.k[ik]*r); else jl=0; // this breaks down for large ls as well l~111 underflow; so I workaround it. */
/*       if (fabs(jl)<1e-50) go=false; */
/*       jlkr[l]=jl; */
/*     }  */

/*     printf("ik %li kr: %lE, lmaxloc: %li \n",ik, t.k[ik]*r, lmaxloc); */
/*     gsl_sf_bessel_jl_steed_array((int)lmaxloc,t.k[ik]*r,jlkr); // these are the spherical Bessel functions: this is the good one */

/*     gsl_sf_bessel_Jn_array(0,(int)lmaxloc,t.k[ik]*r,jlkr); // these are the cylindrical Bessel functions */
/*     for (l=0;l<=lmaxloc;l++) { jlkr[l]=bessel_jn(l,t.k[ik]*r); } */
/*     for (l=0;l<=lmaxloc;l++) { fscanf(f,"%lE",&jlkr[l]);  } */


/*     // dump bessels */
/*     for (l=0;l<=lmaxloc;l++) { fprintf(f1,"%lE\n",jlkr[l]); }  */
/*     for (l=0;l<=lmaxloc;l++) {jlkr[l]=(double)l; } */

/*     kkdkPk = kkdk[ik] * Pk->get_Pi(ik) * twooverPI ; */
/*     // we integrate for all ls at the same time */
/*     for (l=0;l<=lmaxloc;l++) { */
      //printf("l=%li jl=%lE\n",l,jlkr[l]);
/*       Xlr->set_Cl(l, Xlr->get_Cl(l) + get_tf(ik,l,NULL) * cpeds_spherical_bessel_fn(l,t.k[ik]*r) * kkdkPk); */
/*     } */
/*     printf("Pk: %lE\n",Pk->get_Pi(ik)); */

/*     delete [] jlkr; */
/*   } */










  for (l=0;l<=lmaxloc;l++) {

    // define the k integration limits
    ikminnz=get_ikmin_non_zero(l);
    ikmaxnz=get_ikmax_non_zero(l);

    // integrating over k
    for (ik=ikminnz;ik<=ikmaxnz;ik++) {
      fn[ik] = kk[ik] * Pk->get_Pi(ik) * get_tf(ik,l,NULL) * cpeds_spherical_bessel_fn(l,t.k[ik]*r);
    }
    Xlr->set_Cl(l, cpeds_integrate_1d(t.k_num,k,fn,ikminnz,ikmaxnz) );
  }

  delete [] fn;
  delete [] k;
  delete [] kk;
  
  return Xlr;
}






//*********************************************************************************
// operator definitions
//*********************************************************************************
// designed for copying the content of the identical tfs. 
// the transfer function must be initiated
// NOT TESTED IN OTHER CASES...

transfer_function & transfer_function::operator = (transfer_function &src) {
  
  long ik;

  for (ik=0;ik<t.k_num;ik++) {
    t.k[ik]=src.get_k(ik);
    t.f[ik]=*src.get_tf_ik(ik);
  }
  t.kmin=src.get_kmin();
  t.kmax=src.get_kmax();
  t.k_num=src.get_knum();

  t.lmin=0;
  t.lmax=src.get_lmax();
  t.l_num=src.get_l_num();

  shrink_in_lmax();
  tf_loaded=true;
  real_space = src.real_space_tf();
  return *this;

}
