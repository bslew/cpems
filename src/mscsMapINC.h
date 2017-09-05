/*!
  \file mscsMapINC.h - implements generation of inverse noise co-added maps.
*/


#ifndef MSCS_MAP_INC_H_
#define MSCS_MAP_INC_H_

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
#include <QtCore/QList>
#include <string>
#include "Mscs-map.h"
/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */


/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
  \class mscsMapINC
  \brief Implements generation of inverse noise co-added maps
  \details

  \date Apr 23, 2010, 2:29:26 PM
  \author Bartosz Lew
*/
class mscsMapINC : public mscsMap {


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
	mscsMapINC();
	~mscsMapINC();

/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */

	/*!
		\brief adds a map pointed by map to local storage for the INC map making
		\details
		@param map - pointer to the mscsMap object that should be used for making INC
		@param sigma0 - corresponding sigma0 value for that map's DA; Doesn't have to be given for "simple" INC
		@return

		The map is copied so the original map need not to be held in memory. After INC is done
		these copies are removed from the memory.
		\date Apr 24, 2010, 10:05:59 AM
		\author Bartosz Lew
	*/
	void useMap(const mscsMap* map, double sigma0=0);

	mscsMap& getMap(long DA) { return *maps.value(DA); }
	void removeMap(long DA) { delete maps.value(DA); maps.removeAt(DA); }
	QList<mscsMap*> getMaps() { return maps; }

	/*!
		\brief makes INC co-added map from the indicated maps (via useMap method)
		\details
		@param DAst - start DA for making INC; if default value is left the INC will process over all stored maps
		@param DAen - end DA for making INC; if default value is left the INC will process over all stored maps
		@param weight - defines how the INC is to be weighted.\n
		if weight=="simple" then Nobs information is not used
		if weight=="inv_noise" (default) then the INC map is constructed as:
		@param precalib - if true - maps that will be used for INC will be shifted so that their mean was equal the precalibrate value
		@param precalibrate - the value to which the mean of the maps will be set to before doing INC\n
		The mask can be used for precalibration and in that case it should be stored in this object. It will be copied onto invididual maps for precalibration
		and removed just after it.
		@param calib - if true the INC map's mean will be shifted to the value of "calibrate"
		@param calibrate - the value to which the mean if the INC map will be set after INC process. If the mask is loaded then this will be the mean outside of the
		mask.
		@return a const pointer to the INC map
		The INC map is stored in this object.

		\date Apr 24, 2010, 10:43:29 AM
		\author Bartosz Lew
	*/
	const mscsMap* makeINC(long DAst=-1, long DAen=-1, string weight="inv_noise", bool precalib=true, double precalibrate=0, bool calib=true, double calibrate=0);


	void clearMaps();
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
	 QList<mscsMap*> maps; //!< stores maps from individual DAs
	 QList<double> sigmas; //!< keeps information on sigma0 of a given DA

 };
#endif /* MSCS_MAP_INC_H_ */



