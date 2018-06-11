#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "Mscs-common.h"
#include "matrix.h"
#include "Mscs-function.h"

#ifndef MSCS_MINKOWSKI_F
#define MSCS_MINKOWSKI_F


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// this is a basic class for minkowski functionals
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

using namespace std;

class minkowski_f : public mscsFunction {

 public:
  minkowski_f(string ftype, long th_num=0);
  ~minkowski_f();

  double get_nu(long i) { return getX(i); }
  double get_vnu(long i) { return getY(i); }
  void set_nu(long i, double nu) { setarg(i,nu); }
  void set_vnu(long i, double vnu) { setf(i,vnu); }
  string get_type();
  void set_type(string ftype);

//  void get_irange(long *_min_i, long *_max_i);
  long get_size();
  void clear();
//  void get_frange(double * _min_nu, double * _max_nu, double * _min_vnu, double * _max_vnu);
  void normalize_V0();
//  void save(string name);


 private:

//  typedef struct {
//    double nu;
//    double vnu;
//  } v;

//  v * f; // pointer to array holding the functional
  string func_type; // the type of the functional eg. "genus", "circ", "area" etc.
  double min_nu, max_nu, min_vnu, max_vnu; // minimal and maximal values of the functional values in threshold and the value respectively
  long min_i, max_i, size; // minimal and maximal indexes of the functional array and its size
  long min_nui,max_nui,min_vnui,max_vnui; // the indexes of the minimal and maximal nu and vnu values respectively in the object


//  bool check_irange(long i);
//  void calculate_functional_minmax_values();

};

#endif
















////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// this class stores a set of three minkowski functionals for each region 
// of multi-mask. uses the basic minkowski_f class
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

#ifndef MSCS_MINKOWSKI_FS
#define MSCS_MINKOWSKI_FS

class minkowski_fs {

  
 public:

  typedef struct {
    minkowski_f* v0;
    minkowski_f* v1;
    minkowski_f* v2;
    long regid; // stores information on the region number 
    long size; // stores information on the size of the region (number of pixels)
    long dataid; // stores information on the data number 
    //    char[3] datatype; // stores information on the data type (like the type of the simulation Q,W etc)
    long mmid; // stores information on the multi-mask number 
    //    string datastr; // string describing the map data on which the functionals were calculated
    //    string mmstr; // string describing the mm on which the functionals were calculated
  } mfs_type;

  long fnum;
  mfs_type* mfs; // pointer to and array of mink. func.
  //double** mfcovs; // pointer to array of the MFs covariance matrix (cross-correlations actually: f_nu1ijkl*f_nu2ijkl for all combinations of nu1/2
  string object_name;

  minkowski_fs();
  minkowski_fs(long num);
  minkowski_fs(long num, long thnum);
  ~minkowski_fs();

  long get_size();
  void set_name(string name);
  string get_name();
  void save(string name, long how);
  void normalize_V0s();
  double get_mf(long i,long imf,long inu);

 private:

};

#endif









