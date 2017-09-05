#include "cpeds-consts.h"
#include "cpeds-math.h"
#include "Mscs-common.h"
#include "Mscs-global-defs.h"
#include "Mscs-power_spectrum.h"

#ifndef MSCS_MATTER_POWER_SPECTRUM
#define MSCS_MATTER_POWER_SPECTRUM

class matter_power_spectrum {

 public:
/*   typedef struct { */
/*     double k; // effective multipol number */
/*     double P; // power  */
/*     double err; // error */
/*     double cv; // cosmic vatiance */
/*     bool set; // if set flag */
/*   } Pktype; */

/*   Pktype *P; */

  // power spectrum specifications variables
  double kmin,kmax;
  double Pkmax, Pklmin;
  long k_num;
  double *P;
  double *k;
  string object_name;

  // and so on


  //  power_spectrum(long lmax_loc) {
  matter_power_spectrum(long Pksize);
  matter_power_spectrum();
  void initiate_power_spectrum(long Pksize);
  ~matter_power_spectrum();
  double get_kmin();
  double get_kmax();
  double get_P(double kk);
  double get_Pi(long ik);
  double get_k(long ik);
  //void set_Pk(double k,double P) { if (l >= l_Clmin && l <= l_Clmax) { C_l[l].C = C; if (C_l[l].set == false) C_l[l].l = l; C_l[l].set = true; }} // here development could be to resize the object if a call outsize the size is given
  void set_Pi(long ik,double PP);
  void set_k(long ik,double kk);
  void set_P(double kk,double PP);
  void set_k(double kk,double kk_new);
  void set_P_to(double val);
  void reset_Pk();
  string get_name();
  void set_name(string name);


  //void set_err(long l,double err) { if (l >= l_Clmin && l <= l_Clmax) { C_l[l].err = err; } }
  //void set_cv(long l,double cv) { if (l >= l_Clmin && l <= l_Clmax) { C_l[l].cv = cv; } }
/*   void compute_cosmic_variance() { */
/*     long i; */
/*     for (i=0;i<=lmax;i++) { if (C_l[i].set) { C_l[i].cv = cosmic_variance(i); }} */
/*   } */

  void save(string filename);

/*   void save_all(string filename) { */
/*     long l; */
/*     FILE* f; */
/*     f = fopen(filename.c_str(),"w"); */
/*     printf("  -- saving the power spectrum with errors to file %s\n",filename.c_str()); */
/*     for (l=l_Clmin;l<=l_Clmax;l++) {      fprintf(f,"%lE %lE %lE %lE\n",C_l[l].l,C_l[l].C,C_l[l].err,C_l[l].cv);     } */
/*     fclose(f); */
/*   } */

  void kill_Pk();
  void make_Pk(long size);
    

  void load(string filename);
  void flush();


  void make_power_law(double A, double ns, double k0, double* ks);
  void make_power_law_k_3(double A, double ns, double k0, double* ks);
  void make_ones(double* ks);

 private:


};

#endif
