/*!
  \file a general gaussian CMB simulations functionalities
*/

#ifndef MSCS_GAUSSIAN_SIMULATION
#define MSCS_GAUSSIAN_SIMULATION

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
#include "Mscs-map.h"
#include "Mscs-alms.h"

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */


/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class mscsGaussianSimulation
  \brief A build brick for gaussian CMB simulations generation.
  \details
  Encapsulates basic functionalities for generating gaussian CMB simulations.
  This includes eg. generation of white noise maps,
  gaussian maps that realize a particular underlying and provided power spectrum,
  generate noise according to some number-of-observations pattern, and experimental beam transfer function.
  Adding experimental noise to the existing map according to Nobs pattern and sigma_0 - noise amplitude
  for a given receiver.

  \date 2009/06/05 10:36:57
  \author Bartosz Lew
*/
class mscsGaussianSimulation : public mscsMap {


/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PUBLIC MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 public:


/* ------------- */
/* CLASS FRIENDS */
/* ------------- */


/* ---------------------------- */
/* CONSTRUCTORS AND DESTRUCTORS */
/* ---------------------------- */
  mscsGaussianSimulation();
  
  /*!
  \brief generates a gaussian simulation according to given power spectrum instrumental beam transfer function, and instrumental and systematical noise properties.
  \details
  @param ns - healpix resolution parameter of the simulated map
  @param C - angular power spectrum
  @param b - beam transfer function
  @param pixtfType - name of the pixel transfer function to be used. See documentation of SH_synthesis method of the mscsMap class for the possible string values
  @param Nobs - a map holding the number of observations of individual pixels
  @param sigma0 - noise per DA in K*sqrt(N_obs) (as defined in http://arxiv.org/abs/astro-ph/0603452 )
  @param rns - pointer to the random number generator. The generator should be initiated, but if NULL is given it will be initiated. This argument is given to allow
  for fast simulations/generations, and reduce the number of seed initializations from the system time to minimum. The state of the generator is kept in that object.

  \date Apr 22, 2010, 11:59:02 PM
  \author Bartosz Lew
  */
  mscsGaussianSimulation(const long& ns, const long& lmax, const mscsAngularPowerSpectrum& C, const mscsWindowFunction& b, string pixtfType, mscsMap& Nobs, const double& sigma0, cpedsRNG *rns=NULL);


  virtual ~mscsGaussianSimulation() {}

/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */
  /*!
	\brief generate CMB signal map according to the given Cl power spectrum and beam window function
	\details
	@param C - theoretical CMB power spectrum
	@param b - experimental beam transfer function

	The map is stored on the map structure (this is currently only for temperature map)
	The generated CMB signal only alms are stored inside of the object for possible further reuse for map generation with different beam transfer function, with eg.
	void makeCMBsignalMap(const long& lmax, const mscsWindowFunction& b, string pixtfType, cpedsRNG *rns=NULL);\n
	Under current implementation the new alms are not generated in the new call if the lmax number matches with the lmax of the stored alms.
	This is possibly a subject to change in future.

	\date Apr 23, 2010, 4:06:13 PM
  */
  void makeCMBsignalMap(const long& lmax, const mscsAngularPowerSpectrum& C, const mscsWindowFunction& b, string pixtfType, cpedsRNG *rns=NULL);

  /*!
	\brief generates CMB signal map from the alms previously generated and stored inside of the object.
	\details
	  The parameters are the same as in the method above. Note that in the current implementation calling this method is nothing else as calling the inherited SH_synthesis

	\date Apr 24, 2010, 9:12:36 AM
  */
  void makeCMBsignalMap(const long& lmax, const mscsWindowFunction& b, string pixtfType);

  /*!
    \brief Adds a white Gaussian noise to the existing map
    \details
    @param m - mean value of the gaussian noise
    @param s - standard deviation of the gaussian noise
    @param useNobs - indicates whether the number of observations field should
    be used for the noise generation. If true, then the s parameter will be divided by sqrt(Nobs) for a given pixel; otherwise
    it will not be divided.

    The map must be allocated for this routine.
  */
  void addWhiteGaussianNoise(double m, double s, bool useNobs=true, cpedsRNG *rns=NULL);


  void generateGaussianSimulation(const long& ns, const long& lmax, const mscsAngularPowerSpectrum& C, const mscsWindowFunction& b, string pixtfType, mscsMap Nobs, const double& sigma0, cpedsRNG *rns=NULL, string signalMap="");

  mscsAlms& getSignalAlms() { return a; }
  
  //! set flag indicating whether to save the signal only map; not used at the moment
  void setSaveSignalMap(bool v) { saveSignalMap_=v; }


  //  void generate_gaussian_map(strarg DA, long sim_num, int save_partial, strarg how,package_dirs dirs, string wmap_data,string work_scheme);
  //void generate_full_signal_noice_gaussian_map(long nside_loc, long lmax_loc, strarg Clfile, double sigma0, double smooth_beam, strarg WMAPNobsfile, int gaussgenerator, int save_partial, strarg sim_dir, strarg alms_dir, strarg files_prefix, strarg plms_file, strarg show_process, strarg how, int fourier_method);
  //  void generate_full_signal_noice_gaussian_map(long nside_loc, long lmax_loc, strarg Clfile, double sigma0, strarg smooth_beam, double smooth_user, strarg WMAPNobsfile, strarg WMAPNobsfile3, int gaussgenerator, int save_partial, strarg sim_dir, strarg alms_file, strarg files_prefix, strarg files_prefix1, strarg files_prefix3, strarg plms_file, strarg show_process, strarg how, int fourier_method, string wmap_data,string work_scheme);

/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PROTECTED MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 protected:


/* ---------------------------- */
/* PROTECTED METHODS */
/* ---------------------------- */


/* ---------------------------- */
/* PROTECTED STRUCTURES */
/* ---------------------------- */


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
	 mscsAlms a;
	 bool saveSignalMap_;

 };
#endif /* MSCS_GAUSSIAN_SIMULATION */

