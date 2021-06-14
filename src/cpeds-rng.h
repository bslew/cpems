/*!
  \file cpeds-rng.h declares the interface for the random number generator class
*/

#ifndef CPEDSRNGS
#define CPEDSRNGS

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
//#include "cpeds-consts.h"
//#include "cpeds-common.h"
#include "cpeds-math.h"
#include "cpeds-msgs.h"
#include "cpeds-list.h"

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */


/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class cpedsRNG
  \brief Encapsulates the random number generation
  \details

  \date 2009/08/18 13:22:23
  \author Bartosz Lew
*/
class cpedsRNG {


/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PUBLIC MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 public:

	 typedef enum { uniform, gaussian, gaussian_invcdf, gaussian_circle, gaussian_power_law, gaussian_power_law_t, gaussian_power_law_t2, gaussian_power_law_fft, from_array_invCDF, from_array_invCDF2d,
		 from_array_invCDF_sample} distrType;

/* ------------- */
/* CLASS FRIENDS */
/* ------------- */


/* ---------------------------- */
/* CONSTRUCTORS AND DESTRUCTORS */
/* ---------------------------- */

  /*!
    \brief initiate the rng
    \details
    @param distr - desired distribution: \n
    "uniform" - uniform distribution,\n
    "gaussian" - gaussian distribution from central limit theorem,\n
    "gaussian_invcdf" - gaussian distribution using method of inverting gaussian cdf. Note that it will be very slow to use this method with signle RN draws vis getRN()\n
    because at every call a new inverse CDF will be generated - which is not necessary\n
    "gaussian_circle" - draws from unitary circle are used to convert two uniform distribution numbers into a gaussian random number\n
    "gaussian_polar" - polar method: NOT IMPLEMENTED YET\n
    "gaussianPowerLaw" - gaussian with defined power law spectrum. The 1/f numbers are generated in real space\n
    "gaussianPowerLaw_t" - gaussian with defined power law spectrum. Solution rewritten on arrays for speed up\n
    "gaussianPowerLaw_t2" - gaussian with defined power law spectrum. Second solution on arrays but slower than the first one.\n
    "gaussianPowerLaw_fft" - gaussian with defined power law spectrum but generated in Fourier space. Once calculated, all numbers are stored in memory.\n
    "invCDF" - random numbers from provided tabulated CDF distribution.
    "invCDF2d" - random numbers from provided tabulated CDF distribution.
	"invCDF_sample" - as invCDF but the random numbers for sampling CDF are taken from the provided sample set with setRandomSample

    @param rn - indicates data type of random numbers: now only "double" supported
    @param generator - specifies the random number generator to be used (default: gsl_rng_mt19937)
    @param seed_ini - initial value for the seed - default: 0
    @param seed_offset - offset of the seed - will be added to the seed before initializing RNG (default: 0)
  */
  cpedsRNG(string distr="uniform", string rn="double", const gsl_rng_type* generator=gsl_rng_mt19937, long seed_ini=0, long seed_offset=0) {
    _seed=seed_ini;
    _seed=random();
    _seed_offset=seed_offset;
    initiateRNG(distr, rn, generator);
  }

  cpedsRNG(const cpedsRNG& parent) { 
	  initiateRNG("uniform", string("double"), parent._generator_type);
	  clone(parent);
  }

/*   ~cpedsRNG() { gsl_rng_free(_state); killGCDF();  saveSeed(string("/home/blew/nao4/tmp/seed."+cpedsMsgs("").toStr(seed()))); } //NOTE: this is temporary debug stuff */
  ~cpedsRNG() { 
	  gsl_rng_free(_state); 
	  killGCDF();  /* saveSeed(string("seed."+cpedsMsgs("").toStr(seed()))); */ 
	  killGCDF2d(); 
	  if (_bt!=NULL) delete [] _bt;
	  if (_oofnst!=NULL) delete [] _oofnst;
  }
/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */
  /*!
    \brief initiates the RN generator with some parameters
    \details
    @param distr - tells what kind of distribution is needed: can be either  "uniform" or "gaussian"
    @param rn - tells what kind of random numbers should be returned - can be "integer" or "double"
    @param generator - tells which algorithm for the RN generation should be used
    It initiates the state of the RNG with the seed and seed offset. if the seed is 0 then it will be initialized with the current system time + seed_offset.
    The convention in the entire object is that the seed+seedOffset should give the right value with which the RNG was initialized
  */
  void initiateRNG(string distr, string rn, const gsl_rng_type* generator);

  //! clone the parent object
  void clone(const cpedsRNG& parent);

  //! define the distribution type
  void setRNsType(string distr);
  
  //! define the distribution type
  void setRNsType(distrType distr);

  void setRandomSample(cpedsList<double> sample);
  
  /*!
    \brief defines the PDF for RNs generation
  \details 
  @param size - size of the x and p arrays
  @param x - defines the range of tabulated PDF
  @param p - defines the pdf tabulated at values x

  The RNs are drawn with inverse CDF method with half division search
  
  The provided arrays will be used inside the object so do not delete them.
  They will be deleted after the object is destroyed
  \date Nov 22, 2010, 12:59:34 PM
  \author Bartosz Lew
  */
  void setPDF(long size, double* x, double *p);
  /*!
	\brief alternative to setPDF
	\details 
	@param as in setPDF

	\date Jun 12, 2021, 5:22:37 PM
*/
  void setCDF(long size, double* x, double *p);

  /*!
	\brief defines the 2-d PDF for RNs generation
	\details 
	@param sizeX - size of the x arguments array
	@param sizeY - size of the y arguments array
	@param x - defines the range of tabulated PDF in x dimension
	@param y - defines the range of tabulated PDF in y dimension
	@param p - defines the pdf tabulated at values x,y. The size is sizeX*sizeY and ordering is X-major
	@return

	  The RNs are drawn with inverse CDF method with half division search
	  
	  The provided arrays will be used inside the object so do not delete them.
	  They will be deleted after the object is destroyed
	\date Apr 10, 2013, 10:03:09 AM
	\author Bartosz Lew
   */
  void setPDF2d(long sizeX, double* x, long sizeY, double *y, double *p);

//   returns an array of size n of the generated random points from 2D PDF
//  cpedsPoint2D* getRN2d(long n);
  //! returns an array of size n of the generated random numbers
  double* getRN(long n);
  //! returns a list of size n of the generated random numbers
  const cpedsList<double> getRNs(long n);
  /*!
	\brief returns a random number with previously given specifications.
	\details 
	@return a RN
	
	Do not use this method when you chose gaussianPowerLaw_fft as a distribution it will not work.

	\date Nov 18, 2010, 6:39:12 PM
	\author Bartosz Lew
  */
  double getRN();

  //! sets the minimal and maximal value for the uniform RNs
  void setMinMax(double min, double max) { _min=min; _max=max; }
  //! returns the theoretical requested minimal value for the uniform RNs
  double Min() const {return _min; }
  //! returns the theoretical requested maximal value for the uniform RNs
  double Max() const {return _max; }
  //! sets the mean and standard deviation (not variance!!) value for the gaussian RNs
  void setMeanStd(double m, double s) { _mean=m; _stdev=s; }
  //! returns the theoretical requested mean for the gaussian RNs
  double mean() const { return _mean; }
  //! returns the theoretical variance mean for the gaussian RNs
  double gauss_std() const { return _stdev; }
  //! set the number of numbers uniformly distributed that will be used to create the gaussian distributed variate
  void setCentralLimitNumbers(long n) { _centralLimitNumbers=n; }
  //! get the number of numbers uniformly distributed that will be used to create the gaussian distributed variate
  long centralLimitNumbers() const { return _centralLimitNumbers; }
  
  //! set the spectral index for the generation of the power law gaussian noise
  void setSpectralIndex(double alpha=-1.0) { _spectralIndex=alpha; }
  double getSpectralIndex() const { return _spectralIndex; }
  //! set the pivot point (k0) for the generation of the power law gaussian noise
  void setPivotPoint(double k0=1) { _pivotPoint=k0; }
  double getPivotPoint() const { return _pivotPoint; }
  //! set the amplitude of the power law spectrum used for the generation of the power law gaussian noise
  void setOofAmplitude(double A=1) { _oofAmplitude=A; }
  double getOofAmplitude() const { return _oofAmplitude; }
  void setOofKminKmax(double kmin,double kmax) { _oofKmin=kmin; _oofKmax=kmax; }
  
  double getPowerLawNoiseMinimalCoefficient() { return _KasdinPowerLawNoiseMinimalCoefficient; }
  void setPowerLawNoiseMinimalCoefficient(double v) { _KasdinPowerLawNoiseMinimalCoefficient=v; }
  long getPowerLawNoiseCoefficiensCount() const;
  //! returns the total number of generated 1/f numbers
  long long getOOfTot() const {return _oofTot; }
  
  //! returns the initial seed of the generator
  long seed() const { return _seed; }
  //! sets the initial seed of the generator;
  //! \details if this is set to 0 then the current time is used as a number of seconds from 1980 combined with th seedOffset
  void seed(long s);
  //! returns the seed offset
  long seedOffset() const { return _seed_offset; }
  //! sets the seed offset
  void seedOffset(long seed_offset);
  long long drawsCount() const { return _drawsCount; }

  //! return the current RNG type
  distrType getRNsType() { return _distr; }
  
  double* getCDF() { return gCDF; }
  double* getCDFargs() { return CDFargs; }
  
  void saveSeed(string fname) {
    FILE* f;
    f=fopen(fname.c_str(),"a");
    fprintf(f,"%li\n",seed()+seedOffset());
    fclose(f);
  }


  //! a copy operator
  const cpedsRNG& operator=(const cpedsRNG& rhs) { clone(rhs); return *this; }

/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PROTECTED MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 protected:


/* ---------------------------- */
/* PROTECTED METHODS */
/* ---------------------------- */
  gsl_rng* getState() const { return _state; }
  void setState(gsl_rng* s) {   gsl_rng_free(_state);    _state=gsl_rng_clone(s);  }


/* ---------------------------- */
/* PROTECTED STRUCTURES */
/* ---------------------------- */

  //
  // power noise generation variables
  // 
  double _spectralIndex; //!< spectral index n of the f^n power spectrum random numbers. -1 for 1/f noise
  double _pivotPoint; //!< pivot point k0 for definition of the power law for generation of 1/f gaussian power law noise
  double _oofAmplitude; //!< power law spectrum amplitude used for generation of 1/f noise in fourier space
  double _oofKmin; //!< kmin - defines the range of the power spectrum used to generate the 1/f noise
  double _oofKmax; //!< kmax - defines the range of the power spectrum used to generate the 1/f noise
  double _KasdinPowerLawNoiseMinimalCoefficient; //!< the minimal value of the coefficient
  bool _KasdinCoefficientsCalculated;
  cpedsList<double> _b,_oofns; //!< coefficients of the real-space filter and 1/f numbers
  long long _oofTot;
  
  // speed-up tables for generation of 1/f noise in real space	
  double* _bt;
  double* _oofnst;
  long _bt_size;
  long _oofnst_curidx;
  long _oofnst_stfrom;

  
  typedef enum { Double, Integer } rnDataType;
  distrType _distr;
  rnDataType _rnDataType;
  long _seed;
  long _seed_offset;
  double _min,_max;
  double _mean, _stdev;
  long _centralLimitNumbers;
  double *gCDF;
  double *gCDF2d;
  double *CDFargs;
  double *CDFargsX,*CDFargsY; // for the case of 2d PDF
/*   double *gCDFarg; */
  long gCDFsize;
  long gCDFsize2d;
  long gCDFsizeX, gCDFsizeY; // for the case of 2d PDF
  long long _drawsCount;
  cpedsList<double> random_sample;
  
  gsl_rng* _state;
  const gsl_rng_type* _generator_type;
  
  
/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PRIVATE MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */

 private:
  void copyGCDF(const cpedsRNG& parent);
  void killGCDF();
  void killGCDF2d();
/* ---------------------------- */
/* PRIVATE METHODS */
/* ---------------------------- */

  /*!
	\brief the method performs a sum over generated 1/f numbers to derive the next one
	\details 
	@param n - n'th sum; must be greater than 0 since sum(n=0)=0
	@return
	  The sum is:
	  s_n = \sum_{i!=0} b_i*X_{n-i} \n
	  where\n
	  b_i - is the coefficient defined as: b_i=(i-1-ns/2) b_{i-1}/i\n
	  n - defines the number of terms in the summation
	  X_{n-i} - n-i'th generated 1/f number stored in the _oofns object
	  
	\date Nov 4, 2010, 3:42:11 PM
	\author Bartosz Lew
  */
  double sumOOfNumbers(long n);
  /*!
	\brief same as above but operating on arrays for speed-up. 
	\details 
	Currently this is the fastest method of generation in real space

	\date Nov 5, 2010, 12:56:20 PM
	\author Bartosz Lew
  */
  double sumOOfNumbers_t(long n);
  //! same as above but supposed to be faster but it isn't
  double sumOOfNumbers_t2(long n);

/* ---------------------------- */
/* PRIVATE STRUCTURES */
/* ---------------------------- */

  
  
  
  
 };
#endif /* CPEDSRNGS */

