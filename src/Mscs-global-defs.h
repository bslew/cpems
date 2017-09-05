#include "Mscs-common.h"
#include <vector>
//#include <QtCore/QString>
//#include <QtCore/QList>

#ifndef MSCS_GLOBAL_DEFS
#define MSCS_GLOBAL_DEFS
//! release version information
extern string Mscs_version;

//commented out during transition into version-1.0
/* //  */
/* //  */
/* /\*! */
/*   \struct f_stat */
/*   \brief container to keep information on the minkowski functionals */
/*   \details  */

/*   \date 2009/05/27 11:28:31  */
/*   \author Bartosz Lew */
/* *\/ */
/* typedef struct { */
/*   long k_num; //!< number of temperature thresholds */
/*   double *k; //!< the pointer to array of thresholds at which mink. f. is calculated */
/*   double *fkm; //!< the pointer to array of median values of mink. functional  */
/*   double *fkv; //!< the pointer to array of variance values of mink. functional  */
/*   double *fks; //!< the pointer to array of skewness values of mink. functional  */
/*   double *fkk; //!< the pointer to array of kurtosis values of mink. functional  */

/*   // the histogram of values in the temperature thresholds (for visual gaussianity check) */
/*   long bin_num; //!< number of bins */
/*   double *bin; //!< pointer to array defining bins */
/*   double *Nbin; //!< pointer to array of pointers to arrays with counts in the bins (in row major order) */
/* } f_stat; */


/*******************************************************************************/
/*  declaration of shared memory variables that hold running tasks information */
/*******************************************************************************/

//! shared variable that holds information on where to find the current information on the running processes
extern Mscs_run_control_t* MSCS_GLOBAL__RUN_SIMULATION_INFO;
//! key to access the shared memory area for simulations run
extern char MSCS_GLOBAL__RUN_SIMULATION_SHM_KEY[255];




/******************************************************************/
/* //! declaration of variables that hold naming Mscs conventions */
/******************************************************************/
// heaplix pix system transfer function file prefix
extern string MSCS_GLOBAL__HEALPIX_PIXTF_PREF;
// heaplix pix system transfer function file suffix
extern string MSCS_GLOBAL__HEALPIX_PIXTF_SUFF;

// nested coordinates file prefix
extern string MSCS_GLOBAL__NESTED_COORDINATES_PREF;
// nested coordinates file suffix
extern string MSCS_GLOBAL__NESTED_COORDINATES_SUFF;
// ring coordinates file prefix
extern string MSCS_GLOBAL__RING_COORDINATES_PREF;
// ring coordinates file suffix
extern string MSCS_GLOBAL__RING_COORDINATES_SUFF;

// nested 2 ring conversion array file prefix
extern string MSCS_GLOBAL__N2R_CONV_TAB_PREF;
// nested 2 ring conversion array file suffix
extern string MSCS_GLOBAL__N2R_CONV_TAB_SUFF_BIN;
extern string MSCS_GLOBAL__N2R_CONV_TAB_SUFF_TXT;
// ring 2 nested conversion array file prefix
extern string MSCS_GLOBAL__R2N_CONV_TAB_PREF;
// ring 2 nested conversion array file prefix
extern string MSCS_GLOBAL__R2N_CONV_TAB_SUFF_BIN;
extern string MSCS_GLOBAL__R2N_CONV_TAB_SUFF_TXT;



/*************************************************************/
/* // definitions of directories where usefull stuff is kept */
/*************************************************************/
extern string HOME_DIR; //!< users home directory
extern string MSCS_PROGRAM_DIR;//!< mscs package directory - full path
extern string MSCS_WMAP_DATA_DIR;//!< wmap cosmo data directory
extern string MSCS_WMAP_LOCAL_DIR;//!< derived wmap stuff directory
extern string MSCS_DATA_DIR;//!< data directory (with pixel transfer functions, precalculated coordinates, conversion arrays, etc.)

//commented out during transition into version-1.0
/* //! Declarations of various differential assembly receiver  names;  */
/* /\*! \note FIXME: This is really stupid - for each of DA string size memory is allocated ! This should be done by enum type  *\/ */
/* extern string DAnames_nums[];// =  {"k1" "k2" "q1","q2","v1","v2","w1","w2","w3","w4"}; */

/* //! Declarations of various differential assemblies (DAs) names;  */
/* /\*! \note FIXME: This is really stupid - for each of DA string size memory is allocated ! This should be done by enum type  *\/ */
/* extern string DAnames[];// =  {"q","q","v","v","w","w","w","w"}; */

/* extern long DAnums[];// =  {1,2,1,2,1,2,3,4}; */

//! Some global definitions of extern files to be used in programs
//! \note I guess this should be gradually removed from the package
extern string GLOBAL_path;
extern long GLOBAL_nside;
extern long GLOBAL_pix_num;
extern long GLOBAL_lmax;
extern string GLOBAL_file_prefix;
extern string GLOBAL_DA1;
extern int GLOBAL_DA1i;
extern string GLOBAL_DA2;
extern int GLOBAL_DA2i;
extern string GLOBAL_mask_file;
extern string GLOBAL_smoothing_method;
extern int GLOBAL_fourier_method_num;
extern int GLOBAL_fourier_forward_method_num;
extern int GLOBAL_fourier_backward_method_num;

//commented out during transition into version-1.0
/* extern string GLOBAL_Plmsfile_forward; */
/* extern string GLOBAL_Plmsfile_inverse; */

//! initialization of some global variables
void Mscs_initiate_global_variables();
//! prints the current state of global variables to the screen
void Mscs_print_global_variables();
//! definitions of directories where usefull stuff is kept; additionally initiates the global lmax and nside
void Mscs_initiate_global_variables(long _nside, long _lmax);
//! initiates the global directories variables to their default values
void Mscs_initiate_directories();

//commented out during transition into version-1.0
/* //! initiates the global WMAP DAs  variables to their default values */
/* extern void Mscs_initiate_WMAP_DAs(); */


//commented out during transition into version-1.0
/* /\*! initiates the global WMAP DAs  variables to their default values */
/*    \note this name is misleading (fourier) this is about the spherical  */
/*    harmonic transformation, not fourier transformation. */
/*    It should be changed */
/* *\/ */
/* extern void Mscs_initiate_fourier_variables(long _nside, long _lmax); */



/* //! generates an elaborated file name */
/* /\*! \note This is rather to be removed. I don't see it to be useful *\/ */
/* extern string*  make_detailed_file_name(string *name, */
/* 				    strarg pathP, */
/* 				    long nsideP, */
/* 				    long lmaxP, */
/* 				    strarg file_prefixP, */
/* 				    strarg DA1P, */
/*     				    int DA1iP, */
/* 				    strarg DA2P, */
/* 				    int DA2iP, */
/* 				    strarg mask_fileP, */
/* 				    strarg smoothing_methodP, */
/* 			      int fourier_method_numP); */

/* //! generates an elaborated file name with additional comment part */
/* /\*! \note This is rather to be removed. I don't see it to be useful *\/ */
/* extern void  make_detailed_file_name_comment(string *name, strarg comment); */

/* extern void make_covariance_matrix_name(string* name, strarg dir, long start_num, long end_num, long sim_num, long th_num, strarg comment); */

//---------------------
//commented out during transition into version-1.0
/* //! declaration of common routines for programs of Mscs package */
/* //! matrix load/save routines */
void Mscs_matrix_print(matrix<double> *M);
/* extern void Mscs_matrix_save(matrix <double> *M, strarg where); */
/* extern void Mscs_matrix_save(matrix <double> *M, strarg where, string how); */
/* extern matrix <double> * Mscs_matrix_load(strarg where); */
/* extern matrix <double> * Mscs_matrix_load(strarg where, string how); */
/* extern void initiate_f_stat(f_stat *fstat, long th_num,long mink_num,long nkbin); */
/* extern void delete_f_stat(f_stat *fstat,long mink_num); */
//commented out during transition into version-1.0

matrix<double>* Mscs_matrix_bin(matrix<double> *M, long dx, long dy, long work_mode);
/* extern matrix<double>* Mscs_matrix_bin(matrix <double> *M, long dx, long dy); */

/*!
  \brief checks if the string s contains c characters
  \details
  @param s - string to be checed for containing nc caracters from c array
  @return true if contains, false if doesn't
*/
bool Mscs_contains_str(string s, const char* c, long nc);

/*!
  \brief removes the extension from the input file string
  \details
  @param
  @return - input string without the extension
  \note Removes all repetitions of the ext from the file stirng.
*/
string Mscs_basename_noext(string file,string ext);
string Mscs_basename(string name);
string Mscs_getExtension(string name);
vector<string> Mscs_stringToList(string str, string sep=",");
vector<double> Mscs_stringToDoubleList(string str, string sep);
vector<long> Mscs_stringToLongList(string str, string sep);
//QList<QString> Mscs_stringToList(string str, string sep=",");

#endif
