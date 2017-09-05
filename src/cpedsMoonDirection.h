/*!
 * \file cpedsMoonDirection.h
 *
 *  Created on: Dec 5, 2012
 *      Author: blew
 */

#ifndef CPEDSMOONDIRECTION_H_
#define CPEDSMOONDIRECTION_H_

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */

/* STANDALONE HEADERS */

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */

/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
 \class 
 \brief Encapsulates the direction to the moon
 \details 
 
 \date Dec 5, 2012, 1:01:36 PM 
 \author Bartosz Lew
 */

#include "cpeds-direction.h"

class cpedsMoonDirection: public cpedsDirection {
	public:
		
		typedef struct {
			double JD; //!< the moment for which the topocentric position was calculated
			double epochJD; //!< epoch in which the coordinates should be stored
			double phase; //!< current phase angle sun-moon-earth [0..180] [rad]
			double distance; //!< center of moon to center of Earth distance [m]
			double angularSize; //!< angular size of the moon in rad
			double equatorialRadius; //!< [m]
			double polarRadius;  //!< [m]
			double bodyTemperature; //!< [K]
			cpedsDirection observer; //!< observer geodetic coordinates on the surface of ellipsoid [deg]
			double observerAltitude; //!< observer altitude above the ellipsoid [m]
			double temperature; //!< atmospheric temperature at observation stop [Celsius]
			double pressure; //!< observer altitude above the ellipsoid [millibars]
			string ephemeridesFile;
			string COname; //!< name of the celestial object
			int novas_obj_type; //! type of the object for make_object function as defined by novas library
			int novas_body_number; //! number of the body within that type
			
		} properties_t;
		
		
		cpedsMoonDirection(cpedsDirection obsLoc=cpedsDirection(18.56406,53.09546), double altitude=133, double pressure=1012, double temperature=20);
		cpedsMoonDirection(const cpedsMoonDirection& parent);
		virtual ~cpedsMoonDirection();
		
		/*!
			\brief returns the topocentric direction at which the moon currently is
			\details 
			@param JD - JD UT time to calculate direction for. for 0 - the current UT time is taken. 
			@param epochJD - the epoch for which the returned ra-dec coordinates are specified: A value of 2000 corresponds to JD2000
			@param geocentric - coordinates modifier. if false (default) then topocentric coordinates are returned. If true, the geocentric ra,dec are returned instead.
			@return
			
			
		
			\date Dec 5, 2012, 1:04:37 PM
			\author Bartosz Lew
		*/
		virtual DirectionRaDec nowAtRaDec(double JD=0, double epochJD=2000, double ut1_utc=0, bool geocentric=false);
		/*!
			\brief returns optical phase angle in radians for requested JD
			\details 
			@param
			@return
			
			angle between the Sun-Moon-Earth
		
			\date Dec 6, 2012, 9:51:40 AM
			\author Bartosz Lew
		*/
		virtual double opticalPhase(double JD=0);
		virtual double brightLimbPA(double JD=0);
		virtual double illuminatedDiskFraction(double JD=0);
		/*!
			\brief returns the angular size of the moon at JD as seen from the center of the Earth
			\details 
			@param JD - JD for which the distance is calculated. if JD=0 then the current JD for UT time is calculated from system time.
			@param which - flag to indicate which angular size you want: \n
			0 - for the average polar - equatorial radii
			1 - for maximal radii
			-1 - for minimal radii
			@return angular size in radians
		
			The angular distance is based upon the averaged (polar + equatorial)/2 radius of the moon.
			
			\date Dec 5, 2012, 5:18:28 PM
			\author Bartosz Lew
		*/
		virtual double angularSize(double JD=0, int which=0);
		/*!
			\brief center of the Earth - center of the Moon distance in meters for requested JD
			\details 
			@param
			@return distance in meters
			This function relies on the libnova implementation. Should be reimplemented to use JPLEPH ephemerices from novas library
		
			\date Dec 5, 2012, 5:24:29 PM
			\author Bartosz Lew
		*/
		virtual double distance(double JD=0);
		/*!
			\brief get precalculated distance to the object
			\details 
			@return distance to object in [m]
			Returns the actual value of the distance value, set previously by some other command (eg. nowAtRaDec)
		
			\date Dec 13, 2012, 11:12:23 AM
			\author Bartosz Lew
		*/
		virtual double getDistance() const { return _COproperties.distance; }
		virtual double getJD() const { return _COproperties.JD; }
		
		virtual properties_t getProperties() { return _COproperties; }
		
		/*!
			\brief get the direction to the center of the center of the dark part of the moon disk
			\details 
			@param JD - JD in UT time for which the direction is to be calculated
			@param epochJD - epoch in which the coordinates are to be provided (default= J2000 - epochJD=2000, other values indicate particular JD)
			@param dsCtrSize - pointer to a variable to hold the maximal angular size of the dark side part of the moon disk [deg] (not set if null given)
			@return get the direction to the center of the center of the dark part of the moon disk in radians
		
			\date Dec 6, 2012, 11:18:46 AM
			\author Bartosz Lew
		*/
		virtual DirectionRaDec getDarkSideCenter(double JD, double* dsCtrSize=NULL);

		/*!
			\brief calculates an estimated flux density per beam of the planet
			\details 
			@param freq - requested frequency [Hz]
			@param fwhm - gaussian beam size assumed for calculations [rad]
			@param T - planet temperature [K] to be used for calculation (if not given a default value is assumed)
			@return flux density per beam [Jy]
			
			The calculation assumes black body spectral radiance distribution and thermodynamic temperature for a given
		
			\date Jan 7, 2015, 1:52:07 PM
			\author Bartosz Lew
		*/
		virtual double fluxPerBeam(double freq, double fwhm, double T=0);

		
		cpedsMoonDirection& operator=(const cpedsMoonDirection& rhs);
	protected:
		properties_t _COproperties;
};

#endif /* CPEDSMOONDIRECTION_H_ */
