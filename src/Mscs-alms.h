/*!
  \file Declares a class for storing and managing alms (coefficients of the spherical harmonic transformation)
*/

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */

#ifndef MSCS_ALMS
#define MSCS_ALMS

#include "Mscs-object.h"
#include "Mscs-alm.h"
#include "cpeds-list.h"
#include "cpeds-rng.h"
/* #include "Mscs-common.h" */
/* #include "Mscs-global-defs.h" */
#include "Mscs-angular_power_spectrum.h"

// interdependent headers
//#include "Mscs-map.h"

// forward declarations
//class mscsMap;

using namespace std;
/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class mscsAlms
  \brief Encapsulates the spherical harmonic (SH) transformation coefficients
  \details
  This class is needed whenever a SH transformation is done, or any other
  operation are performed in the SH space.

  \date 2009/05/28 15:34:19
  \author Bartosz Lew
*/
class mscsAlms : public mscsObject {

 public:
/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PUBLIC MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
  typedef struct {
    long lmax;
    long int num;
    bool halfTriangle;
  } almsInfoStructure;

/* ------------- */
/* CLASS FRIENDS */
/* ------------- */

/* ---------------------------- */
/* CONSTRUCTORS AND DESTRUCTORS */
/* ---------------------------- */
  mscsAlms();
  mscsAlms(const mscsAlms& parent);
  mscsAlms(const cpedsList<mscsAlm>& al);
  mscsAlms(long lmaxl);
  mscsAlms(string name, long lmaxl=-1);
  //! A cloning constructor
  ~mscsAlms() {}


/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */

  // HANDLERS
  long lmax() const { return almsInfo.lmax; } //!< returns the maximal multipole number in the fourier space
  void lmax(long l) { set_lmax(l);  _alms.makeLength(almsNum()); } //!< sets the alms number lmax and allocates the  alms space accordingly
/*   long get_alms_lmax(); //!< returns the maximal multipole number in the fourier space */   //commented out during transition into version-1.0
  long almsNum() const { return almsInfo.num; } //!< returns the total number of alm complex coefficients
  void info() { msgs->say("lmax: "+msgs->toStr(lmax())+", alms count: "+msgs->toStr(almsNum()),High); }
  const cpedsList<mscsAlm>& alms() const { return _alms; }
  mscsAlm& a(long l,long m);
  mscsAlm& a(long i) { return _alms[i]; }
  cpedsList<mscsAlm>& a() { return _alms; }
  const mscsAlm& get(long l,long m) const;
  const mscsAlm& get(long i) const; //!< returns the i'th alm coefficient as stored in the array  
  const mscsAlms& set(const cpedsList<mscsAlm>& cl) { _alms=cl; return *this; }
  void set(long l, long m, const mscsAlm& aa)  { a(l,m)=aa; } // sets alm value
  void set(long l, long m, double aR, double aI) { a(l,m)=mscsAlm(aR,aI); } // sets alm value
  /* void setR(long i, double v) { a(i).r()=v; } // sets alm value  */
  /* void setI(long i, double v) { a(i).i()=v; } // sets alm value  */


  // DATA OPERATIONS
  void clear();
  void num2alm(long int num, long *l, long *m) const;
  long alm2num(long l, long m) const; // move to private
  long alm2num_norm(long l, long m) const;
  void num2alm_norm(long int num, long *l, long *m) const;
  //! returns the alms in the lm ordering. The current alms are not changed
  mscsAlms tolmOrdering() const;
  void zero_alms_multipole_range(long l1, long l2);
  /* void set_alms_multipole_range(long l1, long l2,string what, double val); */


  /* mscsAngularPowerSpectrum calculate_C_l(long lminloc, long lmaxloc, int how);  // calculates the power spectrum from alms   */
  /*!
    \brief derives the power spectrum from the alms for the requested range of multipoles
    \details
    @param lmin - start calculations at this multipole (default 0)
    @param lmax - end calculations at this multipole (default -1) \n
    if the default value is used the true lmax value of the allocated alms is used for calculations
    @param what - defines how to calculate to power spectrum\n
    what=0 - (default) the whole multipole power is calculated
    what=1 - only the real part power is calculated
    what=2 - only the imaginary part power is calculated
    @return power spectrum object (note that the index in the new object need not the coincide with the corresponding multipole number

    \date 2010/02/22 13:39:03
    \author Bartosz Lew
  */
  mscsAngularPowerSpectrum get_Cl(long lminl=0,long lmaxl=-1,long what=0) const;
  double calculate_single_C_l(long l) const;
  /*!
    \brief calculates power in a single multipole from alms
    \details
    @param l - multipole number
    @param what - defines how to calculate the power:\n
    what = 0 - whole multipole power is calculated
    what = 1 - only real part power
    what = 2 - only imaginary part power
    what = 3 - only m>=0 part power
    what = 4 - only m<0 part power
    @return power in the multipole
  */
  double calculate_single_C_l(long l, long what) const;
  double calculate_single_C_l(long l, long what, long ordering) const;

  /*!
    \brief  calculates the cross power spectrum
    \details
    @param a2 - second set of alms for cross power spectra computation

    The details are as in the get_Cl member function.
    \note what=1,2 are not implemented yet
    \date 2010/02/22 13:44:08
    \author Bartosz Lew
  */
  mscsAngularPowerSpectrum get_cross_Cl(const mscsAlms& a2, long lmin=0,long lmax=-1,long what=0);  // calculates the cross power spectrum on F - i.e. from the alms

  //! calibrates the alms to value cl at some fiducial multipole number l
  const mscsAlms& calibrateOn(long ll, double cl);
  //! removes monopole and dipole
  const mscsAlms& removel01() { zero_alms_multipole_range(0,1); return *this; }

  //SPHERICAL HARMONIC SPACE METHODS

  /*!
    \brief generates random gaussian alms all with the same variance and mean (white noise)
    \details
    @param llmax - defines maximal multipole up to which the random numbers are generated
    @param fm - mean of the generated numbers
    @param fs - standard deviation of the generated numbers
    @param method - method for generating random gaussian numbers (see the code) -- CURRENTLY THIS IS NOT USED - USED TO BE ONE OF THE BELOW

      method: 1 -  the invers CDF from cpeds is used for generagion of random nubers (iterative) ( very old and slow )\n

      method: 2 -  the invers CDF from cpeds is used for generagion of random nubers (by half division) ( here you can play with variace as a function of C_l or set it to a cosmic variance)\n

      method: 3 -  gsl library is used   for generagion of random nubers ( these RNG I haven't harmed yet )\n

      method: 4 -  the invers CDF from cpeds is used for generagion of random nubers (search by half division) but make just one call to the cpeds function for all alms_num at once. here the variance is equal to fs.\n
      in this method the variance throughout all multipoles is the same. nevertheless the underlying power spectrum is respected (exactly as in previous methods )

      method: 5 - just like in method 4 but the power spectrum is trully a random realization - the exact calibration to the underlying Cl is not performed.
      also it is dedicated to fourier transform not Re(a_lm*Y_lm) but simply a_lm * Y_lm - this requires that a_l(-m) = (-1)^m * a_lm* that must be satisfied here -
      the reality condition.
      this makes more than half of the alms redundant. for the DI kind this means that the all phases of a_l0 = 0 and for RI kind it means that all I parts of a_l0 = 0.
      The previous methods were self consistent within this program but probably would not be consistent with external programls calculating fourier transforms, hence
      method 5 is the best for practical applications. this method performs the C_l calibration at l=10 to fit the uncerlying spectrum; for l<0 no calibration is performed
      alm = x*C_l^orig where x is a complex gaussian number (C_l^gen allways is 1 ---> fs = 1) . More precisely:
      In method 5:
      for m > 0, a_lm = R_lm + i * I_lm = zeta_1lm * sqrt(C_l/2) + i * zeta_2lm(C_l/2)
      where zeta_1/2 is driven from N(0,1) so that the expectancy value of the simulated C_l - EX(C_l^sim) = C_l^M - the model power spectrum
      this is granted by the properties of the chisq distribution with n degrees of freedom
      for m = 0, a_lm = R_lm = zeta * sqrt(C_l)
      where zeta is driven from N(0,1)
      for m < 0, a_l-m = (-1)^m * a_lm^* - which is the reality condition for the C -> R fourier transformation
      for this method to work the operative distribution MUST be the normal distribution ---> N(0,1)

      method 6 is like in 5 but fs is a function of C_l or cosmic variance - fs = sqrt(C_l)

    @param half - indicates whether or not only the m>=0 part of the alms should be generated for further antysymmetrization (for real temperature maps)
    By default this is false;
    If false - the full alms will be generated -- as needed for polarization maps.

    @param rns - random numbers object pointer; if null given (default) than a new one will be allocated and destroyed when it is done.\n
    if non-null is given then it will be used and will not be destroyed and can be re-used
    for further calls. In this case after all is done you need to delete it from outside of this method.
    This is useful for very frequent calls to this method (faster than 1Hz) to avoid generation of random numbers from
    the same seed.

    All alms are cleared prior random number generation.
    @return *this is returned as const reference

    \date 2010/02/22 20:17:10
    \author Bartosz Lew
  */
  const mscsAlms& generate_gaussian_alms(long llmax=-1, double fm=0, double fs=1, long method=0, bool half=false, cpedsRNG* rns=NULL); //generates the gaussian alms either with or without the underlying Cl

  /*!
    \brief generates random gaussian alms according to some power spectrum
    \details
    @param cl - the model power spectrun that will be statistically realized by the alms.
    @param exact - indicates whether the power spectrum should be realized with the exact power as defined in the power spectrum (true)
    or within the cosmic variance uncertainty (false - default).
    @param calibrateon - indicates multipole on which the alms are calibrated (if -1 is given (default) - no calibration is done)
    @param method - this is fed to another generate_gaussian_alms (seek there)
    @return *this is returned as const reference

    All alms are cleared prior random number generation.
    All alms are generated (inlcuding those with m<0)
  */
  const mscsAlms& generate_gaussian_alms(const mscsAngularPowerSpectrum& cl, long llmax, long method=0, bool half=false, bool exact=false, long calibrateon=-1, cpedsRNG* rns=NULL); //generates the gaussian alms either with or without the underlying Cl
  //! make alms antysymmetrical; only the information  for m>=0 is used; alm with m<0 are overwritten
  void antysymmetrize();
  /* double calculate_cosmic_variance(double l) const; */








  // IO METHODS

  /*!
    \brief save alms to binary file
    \details
    @param alms_file - file name
    @param how - specifies how the alms should be saved in the file. Possible strings are:\n
    "RI" - saves only real and imaginary part\n
    "lmRI" - saves l m R I --- NOT IMPLEMENTED YET\n
    "lmRIMP" - saves l m R I modulus phase --- NOT IMPLEMENTED YET\n
    "lmMP" - saves l m modulus phase --- NOT IMPLEMENTED YET\n
    @return return code: 0 if O.K., -1 if error occured

    The file format is 2 double numbers (8byte each) for each alm number
    The defult structure for saving alms in file is as follows:\n

    Binary file of double numbers pairs.
    Each complete alm is by default written as RE(alm) IM(alm) for each l and -l<m<l
    The alms ordering is from left to right and from top to bottom of a
    triangular set of alm numbers eg.

                                     (4)\n
                                     a00\n
                              (2)    (5)     (7)\n
                             a1-1    a10     a11\n
                      (1)     (3)    (6)     (8)    (9)\n
                     a2-2    a2-1    a20     a21    a22\n
    and so on

    \date 2010/02/18 23:43:44
    \author Bartosz Lew
  */
  long savebinAlms(string  alms_file, string how) const; // saves alms to a file in a requested manner of which depends the output file name ending
  //! saves alms to txt file: description as for savebinAlms
  long savetxtAlms(string  alms_file, string how) const; // how - 1 -- reads/saves all alms to file ending with "-Falxxx" .bin or .txt

  /*!
    \brief loads alms from binary file
    \details
    @param alms_file - file name
    @param how - specifies how the alms should be saved in the file. Possible strings are:\n
    "RI" - saves only real and imaginary part\n
    "lmRI" - saves l m R I --- NOT IMPLEMENTED YET\n
    "lmRIMP" - saves l m R I modulus phase --- NOT IMPLEMENTED YET\n
    "lmMP" - saves l m modulus phase --- NOT IMPLEMENTED YET\n
    @return return code: 0 if O.K., -1 if error occured

    Details as in savebinAlms
    \date 2010/02/18 23:43:44
    \author Bartosz Lew
  */
  long loadbinAlms(string  alms_file, string how);
  //! loades alms from txt file: description as for savebinAlms
  long loadtxtAlms(string  alms_file, string how);

  //! saves alms to txt file: description as for savebinAlms -- NOT IMPLEMENTED YET
  long savetxtAlms(string  alms_file, string how, long l) const; // how - 1 -- reads/saves all alms to file ending with "-Falxxx" .bin or .txt
  //! loades alms from binary file: description as for savebinAlms -- NOT IMPLEMENTED YET
  long loadbinAlms(string  alms_file, string how, long l); // read in alms from a file

  void printtxtAlms(string what="") const;
  //! clears all alms and deallocates space
  /* void killAlms(); // kills the alms from the memory */

  //commented out during transition into version-1.0
  void loadsave_manager(string whattodo, string fileformat, string what, int how, string where, long whatmultipole);


  mscsAlm& operator[](long i) { return a(i); }
  const mscsAlms& operator=(const cpedsList<mscsAlm>& cl) { set(cl); return *this; }
  const mscsAlms& operator=(const mscsAlms& rhs) {   if (this != &rhs) { almsInfo=rhs.almsInfo;  a()=rhs.alms();  this->mscsObject::operator=(rhs); } return *this; }


 protected:
/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PROTECTED MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */

/* ---------------------------- */
/* PROTECTED METHODS */
/* ---------------------------- */

  //! sets the maximal multipole number and adjusts the number of alms accordingly
  void set_lmax(long l) { almsInfo.lmax=l; almsInfo.num=alm2num(l,l)+1;  }
  //! sets the number of alms - does not alter the memory allocation
  void set_alms_num(long alms_num);

  const cpedsList<double> toDoubleList() const;
  void read_binalm_parameters(string alm_file); //!< reads the alms information from a binary alms file
  void read_txtalm_parameters(string alm_file); //!< reads the alms information from a text alms file


/* ---------------------------- */
/* PROTECTED STRUCTURES */
/* ---------------------------- */
  //SH SPACE RELATED VARIABLES
  cpedsList<mscsAlm> _alms;
  almsInfoStructure almsInfo;

/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PRIVATE MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 private:
/* ---------------------------- */
/* PRIVATE METHODS */
/* ---------------------------- */


/* ---------------------------- */
/* PRIVATE STRUCTURES */
/* ---------------------------- */

  };
#endif



/*   void precalculate_fourier_stuff(double* COS, double * SIN, pixel * mapRING,long int ring_num, long int * ringtab, long int *ringpixtab); */
/*   void precalculate_fourier_stuff(double *COS, pixel *mapRING,long int ring_num,long int *ringtab); */
/*   void precalculate_fourier_stuff(long ring_num, pixel *mapRING, double *ThVals,long *ThBreaks, double* phi0, double * dphi); */


/*   void calculate_Plm_harmonics(double* COS, long int ring_num, long int rmax,double* Ptab,long l, long m); */
/*   void calculate_Plm_harmonics(long int ring_num,double* Ptab, Plms_tab * Plms); */
/*   void calculate_Pmth_harmonics(long m, long skip, double * Ptab, Plms_tab * Plms); */
/*   void calculate_Pmth_harmonics(double * Ptab, long m, double costh); */
/*   void precalculate_Legendre_polynominals(strarg ordering); */


  /* double get_Plm_harmonic(double* Ptab, long int pix); */


  /* double Klm(long l, long m); */
  /* double Nlm_harmonic(long l, long m); */
  /* double Plm_harmonic(long l, long m, double x); */
  /* double Plmp_harmonic(long l, long m, double x); */


  //commented out during transition into version-1.0
/*   void set_alms_lmax(long l); //sets the maximal multipole number in the fourier space */
/*   void set_alms_lmax(long l,string how); //sets the maximal multipole number in the fourier space */

  /* void set_lmax(long l,string how); //sets the maximal multipole number in the fourier space */

  /* void calculate_C_th(double theta_min, double theta_max, double resolution); // calculates the correlation function from the power spectrum  */

  //  a_lm get_F(long i); // returns alm  value   //commented out during transition into version-1.0
  //a_lm get_F(long l, long m); // returns alm  value   //commented out during transition into version-1.0
  /* complex get(long l, long m); //!< returns the alm  value */

  //void set_F(long l, long m, a_lm a); // sets alm value    //commented out during transition into version-1.0

  //void set_F(long l, long m, double aR, double aI); // sets alm value   //commented out during transition into version-1.0
