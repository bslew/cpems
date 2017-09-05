#include "Mscs-matter_power_spectrum.h"

  //  power_spectrum(long lmax_loc) {
matter_power_spectrum::matter_power_spectrum(long Pksize) {
  //long i;
  P=NULL;
  k=NULL;
  initiate_power_spectrum(Pksize);
  object_name="Pk";
}

matter_power_spectrum:: matter_power_spectrum() {
  k_num=0;
  P=NULL;
  k=NULL;
  object_name="Pk";
}

void matter_power_spectrum::initiate_power_spectrum(long Pksize) {
  if (Pksize >= 0) { k_num = Pksize-1; } else { kmax = 0; Pksize=1; }
  make_Pk(Pksize);
  reset_Pk();
    
}

matter_power_spectrum::~matter_power_spectrum() {
  kill_Pk();
}

double matter_power_spectrum::get_kmin() { return k[0]; }
double matter_power_spectrum::get_kmax() { return k[k_num-1]; }
double matter_power_spectrum::get_P(double kk) { return P[cpeds_find_value(kk,k,k_num,0,k_num)]; } // assumes that k is ascending ordered 
double matter_power_spectrum::get_Pi(long ik) { return P[ik]; }
double matter_power_spectrum::get_k(long ik) { return k[ik]; }
//void set_Pk(double k,double P) { if (l >= l_Clmin && l <= l_Clmax) { C_l[l].C = C; if (C_l[l].set == false) C_l[l].l = l; C_l[l].set = true; }} // here development could be to resize the object if a call outsize the size is given
void matter_power_spectrum::set_Pi(long ik,double PP) { P[ik]=PP; }
void matter_power_spectrum::set_k(long ik,double kk) { k[ik]=kk; }
void matter_power_spectrum::set_P(double kk,double PP) { P[cpeds_find_value(kk,k,k_num,0,k_num)]=PP; }
void matter_power_spectrum::set_k(double kk,double kk_new) { k[cpeds_find_value(kk,k,k_num,0,k_num)]=kk_new; }
void matter_power_spectrum::set_P_to(double val) { long i; for (i=0;i<k_num;i++)  P[i] = val;  }
void matter_power_spectrum::reset_Pk() { long i;   for (i=0;i<k_num;i++) { P[i]=0; k[i] = 0; } }
void matter_power_spectrum::set_name(string name) {
  printf("|%s> -- changing name of object to: %s\n", object_name.c_str(), name.c_str());
  object_name=name; 
}
string matter_power_spectrum::get_name() {  return object_name; }


//void set_err(long l,double err) { if (l >= l_Clmin && l <= l_Clmax) { C_l[l].err = err; } }
//void set_cv(long l,double cv) { if (l >= l_Clmin && l <= l_Clmax) { C_l[l].cv = cv; } }
/*   void compute_cosmic_variance() { */
/*     long i; */
/*     for (i=0;i<=lmax;i++) { if (C_l[i].set) { C_l[i].cv = cosmic_variance(i); }} */
/*   } */

void matter_power_spectrum::save(string filename) {
  long i;
  FILE* f;
  f = fopen(filename.c_str(),"w");
  printf("|%s>  -- saving the matter power spectrum to file %s\n",object_name.c_str(),filename.c_str());
  for (i=0;i<k_num;i++) {      fprintf(f,"%lE %lE\n",k[i],P[i]);     }
  fclose(f);
}

/*   void save_all(string filename) { */
/*     long l; */
/*     FILE* f; */
/*     f = fopen(filename.c_str(),"w"); */
/*     printf("  -- saving the power spectrum with errors to file %s\n",filename.c_str()); */
/*     for (l=l_Clmin;l<=l_Clmax;l++) {      fprintf(f,"%lE %lE %lE %lE\n",C_l[l].l,C_l[l].C,C_l[l].err,C_l[l].cv);     } */
/*     fclose(f); */
/*   } */

void matter_power_spectrum::kill_Pk() { 
  printf("|%s>  -- killing the matter power spectrum structures\n",object_name.c_str());
  if (P!=NULL) { delete [] P; P=NULL; }  
  if (k!=NULL) { delete [] k; k=NULL; } 
}

void matter_power_spectrum::make_Pk(long size) { 
  printf("|%s>  -- making the matter power spectrum structure of size %li\n",object_name.c_str(),size);
  kill_Pk(); 
  P = new double[size];    
  k = new double[size];    
  k_num=size;
    
  reset_Pk();
}

    

void matter_power_spectrum::load(string filename) {
  long i,n=0;
  FILE* f;
  f = fopen(filename.c_str(),"r"); if (f != NULL) {  n=0; while (fscanf(f,"%*E %*E") != EOF) { n++;}         fclose(f); } else { printf("|%s>ERROR: no such file %s\n",object_name.c_str(),filename.c_str()); }
  printf("|%s>  -- the file to load has %li lines\n",object_name.c_str(),n);
  if (n != k_num) { printf("|%s>  -- resizing the existing matter power spectra to size: %li\n",object_name.c_str(),n);  make_Pk(n);  }
  printf("|%s>  -- reading the power spectrum from file %s\n",object_name.c_str(),filename.c_str());
  kill_Pk();
  f = fopen(filename.c_str(),"r"); for (i=0;i<k_num;i++) {      fscanf(f,"%lE %lE",&k[i],&P[i]);     }
  fclose(f);
}

void matter_power_spectrum::flush() {
  long i;
  printf("\n");
  for (i=0;i<k_num;i++) {      printf("|%s>%lE %lE\n",object_name.c_str(),k[i],P[i]);     }
  printf("\n");
}


/******************************************************************************************/
//seeds the power spectrum with the power law spectrum P(k) = A(k/k0)^(ns-1)
// ks should be in increasing order
void matter_power_spectrum::make_power_law(double A, double ns, double k0, double* ks) {
  long ik;
  double nsleo=ns-1.0;

  printf("|%s>  -- generating the matter power law power spectrum A=%lE ns=%lE k0=%lE for given tab of ks\n",object_name.c_str(),A,ns,k0);
  for (ik=0;ik<k_num;ik++) {
    set_k(ik,ks[ik]);
    set_Pi(ik,A*pow(ks[ik]/k0,nsleo));
  }

}
/******************************************************************************************/
//seeds the power spectrum with the power law spectrum P(k) = A(k/k0)^(ns-1)
// ks should be in increasing order
void matter_power_spectrum::make_power_law_k_3(double A, double ns, double k0, double* ks) {
  long ik;
  double nsleo=ns-1,k3;

  printf("|%s>  -- generating the matter power law power spectrum A=%lE ns=%lE k0=%lE for given tab of ks\n",object_name.c_str(),A,ns,k0);
  for (ik=0;ik<k_num;ik++) {
    set_k(ik,ks[ik]);
    k3=ks[ik]/k0; k3=k3*k3*k3;
    set_Pi(ik,A*pow(ks[ik]/k0,nsleo)/k3);
  }

}

/******************************************************************************************/
//seeds the power spectrum with ones
// ks should be in increasing order -- define the k wavenumber for which to set the power spectrum
void matter_power_spectrum::make_ones(double* ks) {
  long ik;

  printf("|%s>  -- setting matter power spectrum to ones\n",object_name.c_str());
  for (ik=0;ik<k_num;ik++) {
    set_k(ik,ks[ik]);
    set_Pi(ik,1.0);
  }

}
