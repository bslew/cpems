/*!
  \file handles WMAP-like simulations
*/

#ifndef MSCS_WMAP_SIMULATION
#define MSCS_WMAP_SIMULATION

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
#include "Mscs-map.h"
#include "Mscs-alms.h"
#include "Mscs-WMAPspecifications.h"
#include "Mscs-gaussian_simulation.h"
#include "mscsMapINC.h"

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */


/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class mscsWMAPsimulation
  \brief Encapsulates a WMAP-like maps and simulations
  \details

  \date 2009/06/05 16:23:40
  \author Bartosz Lew
*/
class mscsWMAPsimulation : public mscsGaussianSimulation, public mscsWMAPspecifications {


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
  mscsWMAPsimulation() { }
  mscsWMAPsimulation(DAnames DA, WMAPversion year, const long ns, const long& lmax, const mscsAngularPowerSpectrum& C, string pixtfType, cpedsRNG *rns=NULL, bool precalib=true, double precalibrate=0.0, bool _calib=true, double _calibrate=0.0, mscsMap* mask=NULL);

  virtual ~mscsWMAPsimulation();

/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */
  void makeWMAPsimulation();
  const mscsMap* getSimulation(string name);
  /*!
	\brief defines the directory and naming conventions for the current simulation
	\details
	@param outDir - name of the output sub-directory
	@param pref - file name prefix
    @param suff - file name suffix
    @param num - simulation id number

	\date May 10, 2010, 1:28:42 PM
  */
  void setNames(string outDir, string pref, long num, string suff);
  //! saves the simulation to disk using the previously set naming conventions with setNames method
  void saveSimulation();
  //! saves the simulation to disk using naming conventions given in arguments
  void saveSimulation(string outDir, string pref, long num, string suff, bool saveSignalAlms=false);

  void setSavePartialSimulations(bool val) { savePartialSimulations_=val; }

//  void setBeamTf(mscsWindowFunction& wf) { _beamtf=wf; }
  
  /*!
	\brief saves to disk partial simulations if the option is enabled with setSavePartialSimulations(true)
	\details
	@param outDir - name of the output sub-directory
	@param pref - file name prefix
    @param suff - file name suffix
    @param num - simulation id number

	\date May 10, 2010, 1:28:42 PM
  */
  void savePartialSimulations(string outDir, string pref, long num, string suff);
  /*!
	\brief saves to disk partial simulations if the option is enabled with setSavePartialSimulations(true)
	\details
	  The values for the prefix suffix and simulaion id number are those previously set with setNames command

	\date May 10, 2010, 1:28:42 PM
  */
  void savePartialSimulations();
  
  string getPartialSimulationFileName(DAnames da);

//  /*!
//    \brief generates gaussian simulation of the WMAP requested channels
//    \details
//    @param
//    @return
//
//    \date 2009/06/05 17:14:02
//    \author Bartosz Lew
//  */
//  mscsMap* generate_gaussian_maps_WMAP_QVW(long sim_idx,string how,package_dirs dirs,string wmap_data, string mask, double precalibrate, double calibrate,string what,string work_scheme);

//  const string getSubDirName(DAnames DA, WMAPversion year);

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
	 cpedsRNG *rng; //!< RNG for the simulations
	 QList<mscsMap*> sim;
	 mscsMapINC inc;
	 DAnames _DA;
	 WMAPversion _dr;
	 long _lmax, _ns, _simNo;
	 string _pixtfType, _outDir, _outFilePref, _outFileSuff;
	 bool savePartialSimulations_, _precalib, _calib;
	 mscsAngularPowerSpectrum _cl;
//	 mscsWindowFunction _beamtf;
	 double _precalibrate, _calibrate;
	 mscsMap* _mask;


 };
#endif /* MSCS_WMAP_SIMULATION */

