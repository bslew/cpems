// This class implements handling of radiative transfer functions for eg. NG map generation etc.
#ifndef _ISOC99_SOURCE
#define _ISOC99_SOURCE
#endif
#include <features.h>
#include <stdlib.h>
#include <stdio.h>
#include "Mscs-common.h"
#include "Mscs-global-defs.h"
#include "Mscs-power_spectrum.h"
#include "Mscs-matter_power_spectrum.h"
#include "cpeds-math.h"
//#include "boost/math/special_functions/bessel.hpp"
//#include "bessel_yv.hpp"


class transfer_function {
 public:

  // data types definitions

  typedef struct {
    power_spectrum *f; // array keeping the tf. vals.
    double *k; // array keeping the k-vals. @ which tf is probed
    double kmin,kmax;
    long k_num;    
    long lmin,lmax,l_num; // l_num=lmax-lmin+1; // lmax is the globally largest used l in tf
  } tf_type;

  // data structures definitions
  tf_type t;
  string object_name;
  bool tf_loaded,real_space;
  string tf_load_file;

  //long lmax_all;
  

  // methods definitions

  transfer_function(); // minimal constructor 
  transfer_function(long ksize, long lsize); // constructor with size in l and k directions
  transfer_function(double kmin_loc, double kmax_loc, long k_num, long lmax ); // constructor for allocating memory of requested size with equal spacing in k
  transfer_function(transfer_function* src,string clonename); // constructor for cloning tfs

  ~transfer_function();

  // IO
  int load_tf(string fname);
  bool is_loaded();
  int save_tf(string fname, string how);
  void import_ithCl(power_spectrum *al, long i);
  
  // handlers
  void set_name(string name);
  string get_name();
  double get_tf(double k, long l, int *result);
  double get_tf(long ik, long l, int *result);
  int set_tf(long ik, double k, long l, double tf);
  long get_lmax(double k);
  long get_lmax();
  long get_lsize(long ik);
  long get_l_num();
  long get_lsize(double k);
  long set_lmax(double k);
  long set_all_lmax();
  double get_kmax();
  //  void set_kmax();
  double get_kmin();
  long get_knum();
  double* get_kstab();
  double* get_ks2tab();
  double get_k(long ik);
  void set_k(long ik, double k);
  long get_ikmin_non_zero(long l);
  long get_ikmax_non_zero(long l);
  power_spectrum* get_tf_ik(long ik);
  //  double set_kmin();
  void clear_tf();
  void kill_k();
  void kill_f();


  // other methods
  void update_self(); // updates the proper ranges of kmin and kmax and lmin and lmax
  void sort_in_ks(long ord);
  //power_spectrum* get_real_tf(double r); // compute trhe real transfer function for a given comoving distance r
  power_spectrum* get_real_tf(double r, double * k, double *kk); // compute trhe real transfer function for a given comoving distance r
  transfer_function* get_real_space_transfer_functions(double rmin,double rmax, long Nsteps);
  bool real_space_tf(); // tells whether the tf was converted into the real space or not
  void shrink_in_lmax(); // checks for larges nonzero l in each k in tf and rescales the tf accordingly

/*   power_spectrum* derive_Cl(double A, double ns, double k0); // derives the C_l given power-law primordial power spectrum with given A and ns */
  power_spectrum* derive_Cl(matter_power_spectrum* Pk); // derives the C_l given power-law primordial power spectrum with given P(k) in object matter_power_spectrum 
  power_spectrum* integrate_alphabeta_l(double r, matter_power_spectrum* Pk); // derives the alpha_l(r) or beta_l(r) as defined in astro-ph/0305189 depending on what goes into P(k)



//*********************************************************************************
// operator definitions
//*********************************************************************************

  transfer_function & operator = (transfer_function &src);
 private:
  void initiate_tf(double kmin_loc, double kmax_loc, long k_num, long lmin_loc, long lmax_loc, long l_num_loc ); 

};
