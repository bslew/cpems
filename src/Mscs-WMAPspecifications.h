/*!
  \file stores WMAP technical specifications
*/

#ifndef MSCS_WMAP_SPECIFICATIONS
#define MSCS_WMAP_SPECIFICATIONS

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
#include <string>
#include <complex>
#include "MscsWMAPdata.h"
#include "Mscs-map.h"

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */


/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class mscsWMAPspecifications
  \brief Encapsulates various technical details of the WMAP satellite.
  \details
  Contains details on WMAP beam transfer functions, noise specifications,
  frequency coverage, receivers, number of receivers, beam resolutions etc.


  \date 2009/06/05 13:16:57
  \author Bartosz Lew
*/
class mscsWMAPspecifications {


/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PUBLIC MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 public:


	 /*!
	  \enum DAnames
	  \brief names for DAs used in the code
	  \details

	 inc - refers to the inverse noise coadded maps.
	 all - refers to a property common for all DAs of the same band. Don't use these for gaussian map generation. These can be used
	 eg for choosing FWHM
	 unitBeamTf - refers to a simulation with unit beam transfer function

	  \date Apr 27, 2010, 1:55:35 PM
	 */
	 typedef enum { K1, K2, Kall, Q1, Q2, Qall, V1, V2, Vall, W1, W2, W3, W4, Wall, Kinc, Qinc, Vinc, Winc, VWinc, QVWinc, Q_V_W_inc, Q_V_W_QVW_inc, V_W_inc, V_W_VW_inc, unknownDA } DAnames;
	 typedef enum { WMAP1yr, WMAP3yrs, WMAP5yrs, WMAP7yrs } WMAPversion;

/* ------------- */
/* CLASS FRIENDS */
/* ------------- */


/* ---------------------------- */
/* CONSTRUCTORS AND DESTRUCTORS */
/* ---------------------------- */
  mscsWMAPspecifications();
  virtual ~mscsWMAPspecifications() {}

/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */

  /*!
    \brief returns the FWHM for a given DA in degrees
    \details
    @param da - the requested DA
    @return the FWHM in [deg]; -1 when wrong argument is given

    \date 2009/06/05 13:42:36
    \author Bartosz Lew
  */
  double get_FWHM(DAnames da);

  /*!
    \brief returns the sigma_0 values for a given DA and WMAP release
    \details
    @param
    @return sigma_0 noise parameter; -1 when wrong arguments are provided

    \date 2009/06/05 13:43:40
    \author Bartosz Lew
  */
  double get_sigma0(DAnames da, WMAPversion v);
  /*!
	\brief returns a map with initiated Nobs for requested DA and for requested year of obserwation
	\details

	\date Apr 27, 2010, 12:36:05 PM
  */
  mscsMap get_Nobs(DAnames da, WMAPversion v);
  /*!
	\brief returns a name of the file where the Nobs information is stored for a given DA and year of observation
	\details
	@param
	@return

	\date Apr 27, 2010, 12:36:57 PM
  */
  const string get_Nobs_fileName(DAnames da, WMAPversion v);

  const mscsWindowFunction get_beamtf(DAnames da, WMAPversion v);

  string get_beamtf_fileName(DAnames da, WMAPversion v);


  /*!
	\brief returns the calibration of the power spectrum for requested data release
	\details
	@param v - WMAP data release
	@return power spectrum calibration in Kelvins^2.\n Remember to multiply the power spectrum by T_cmb^2 when converting from CMBfast format into
	Kelvins

	\date Apr 28, 2010, 12:15:08 AM
  */
  double get_Cl_calibration(WMAPversion v);

  //! converts the WMAP version information into corresponding string
  string getVersion(WMAPversion v);
  //! converts the WMAP version given in string into corresponding WMAPversion type value
  WMAPversion getVersion(string v);

  //! returns the name of the DA corresponding to the given string: NOT IMPLEMENTED YET
  DAnames getDA(string da);
  //! returns a string corresponding to a given  DA
  string getDA(DAnames da);

/*   double operator() (DAnames */

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

  //! Normalization of the power spectrum at quadrupole
  double C_2WMAP;


  //! Definition of sigma0 per each frequency channel for the most actual data set (wmap5 for now)
  double sigma0_K1;
  double sigma0_K2;
  double sigma0_Q1;
  double sigma0_Q2;
  double sigma0_V1;
  double sigma0_V2;
  double sigma0_W1;
  double sigma0_W2;
  double sigma0_W3;
  double sigma0_W4;

  /*!
    the array sigma0_WMAP[] stores all the sigma0 values for all
    years and all channels such that first ten is the first year
    data from k1 to w4 second 10 is the three year data and third
    ten i.e. from 20..29 is for WMAP5yrs and from 30..39 for WMAP7yrs
  */
  double sigma0_WMAP[40];

  //! WMAP Beam sizes - the Full Width at Half Maximum in degrees for different DA
  double thFWHM_K,thFWHM_Q, thFWHM_V, thFWHM_W; // K is not initialized yet


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
#endif /* MSCS_WMAP_SPECIFICATIONS */




