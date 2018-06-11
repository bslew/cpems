/*!
  \file Mscs-minkowski-stat.h - 
*/


#ifndef SRC_MSCS_MINKOWSKI_STAT_H_
#define SRC_MSCS_MINKOWSKI_STAT_H_




////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// this class stores a set of three minkowski functionals for each region 
// of multi-mask.
// dedicated class space optimized for MFstat project
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

  // ASSUMED INDEX CONVENTIONS

  // i - MF type {0,1,2}
  // j - sim number {1...10^4}
  // k - reg number {1..48}, {1..192}, {1..768}
  // l - smoothing length index {0..2} for eg. {0.3, 1, 3} deg
  // m - mm number {1..150} for HP2,4,8; 50 mms in each regsk. type
  // d - data type {0,1,2,3,4,5} eg for {Q,V,W,VW} for wmap5 etc.


class minkowski_fs_MFstat {

  
 public:

  // 
  // structures types
  //
  typedef struct {
    double* th; // array of MF thresholds for all thress MFs 
    double* v0; // array of MF values in consecutive thresholds for v0
    double* v1; // array of MF values in consecutive thresholds for v1
    double* v2; // array of MF values in consecutive thresholds for v2
    float* Cv0; // covaraince matrix for v0; (only the triangular part of the matrix is stored)
    float* Cv1; // covaraince matrix for v1; (only the triangular part of the matrix is stored)
    float* Cv2; // covaraince matrix for v2; (only the triangular part of the matrix is stored)
    double m,s,S,K; // mean, variance, skewness, kurtosis in the region
    long regid; // stores information on the region number 
    long size; // stores information on the size of the region (number of pixels)
    //long dataid; // stores information on the data number 
    //    char[3] datatype; // stores information on the data type (like the type of the simulation Q,W etc)
    //long mmid; // stores information on the multi-mask number 
    //    string datastr; // string describing the map data on which the functionals were calculated
    //    string mmstr; // string describing the mm on which the functionals were calculated
  } mfs_type;

  typedef struct {
    mfs_type *mfr; // MFs in a reigons of mm array
    long reg_num; // number of regions in a single mm
    double* Cm;
    double* Cs;
    double* CS;
    double* CK;
  } mfsreg_type;

  typedef struct {
    mfsreg_type *mfm; // MFs in multi-masks in all regions array
    long mm_num; // number of multimasks
  } mfsmm_type;

  typedef struct {
    mfsmm_type *mfs; // MFs for all multi-masks array
    long sm_num; // size of sm
    double *sm; // smoothing lenghts
  } mfssm_type;

  typedef struct {
    mfssm_type *mfd; // MFs for all datasets  array
    long data_num; // size of d
    string *dstr; // data set name
  } mfsd_type;

  //
  // data structures of the object
  //
  mfsd_type* MF; // pointer to a structure containing the info on the MFs for each dataset, each, smoothing length, enach multimask, each region and each MF and it's covariance matrix at each threshold.
  string object_name, mask,simtype;
  long simst,simen,node,thnum;
  string* mmsf;
  

  //
  // methods
  //
  minkowski_fs_MFstat(long node_no, string sim_type, long sim_stnum, long sim_ennum, long sm_num, double* sm, long mm_num, string* mmf, long th_num, string mask);
  ~minkowski_fs_MFstat();

  void make_cov_and_mean_part(long node);
  void make_cov_and_mean(long _simst, long _simen);
  void correlate_and_add_mfs(long d,long l,long m, minkowski_fs* f);
  void invert_covariance_matrices();
/*   long get_size(); */

  void set_name(string name);
  string get_name();

  // data handlers
  double* get_mf(long d, long l, long m, long k, long mftype );
  float* get_mfcov(long d, long l, long m, long k, long covtype );
  void set_mf(long d, long l, long m, long k, long mftype, double* data );
  void set_mfcov(long d, long l, long m, long k, long covtype, float* data );
  void add_mf   (long d, long l, long m, long k, long mftype, double* data);
  void add_mfcov(long d, long l, long m, long k, long covtype, float*  data);

  long get_data_num();
  long get_sm_num(long d);
  long get_mm_num(long d,long l);
  long get_reg_num(long d,long l,long m);
  long get_thres_num();
  long get_reg_id(long d,long l,long m, long k);
  long get_reg_size(long d,long l,long m, long k);
  long get_cov_size();
/*   double get_mf(long d, long l, long m, long k, long nu); */
/*   double get_mfcov(long d, long l, long m, long k, long nu1, long nu2); */


  void save(string name, long how);
/*   void normalize_V0s(); */

 private:

};


#endif /* SRC_MSCS_MINKOWSKI_STAT_H_ */ 



