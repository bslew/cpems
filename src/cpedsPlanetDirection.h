/*!
 * \file cpedsPlanetDirection.h
 *
 *  Created on: Dec 13, 2012
 *      Author: blew
 */

#ifndef CPEDSPLANETDIRECTION_H_
#define CPEDSPLANETDIRECTION_H_

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
 \brief Encapsulates directions to planets of the Solar System
 \details 
 
 \date Dec 13, 2012, 9:07:23 AM 
 \author Bartosz Lew
 */

#include "cpedsMoonDirection.h"

class cpedsPlanetDirection: public cpedsMoonDirection {
	public:

		/*!
			\brief generates an object which can be queried for position
			\details 
			@param planetName - name of the planet. The possible values are:\n
			Merkury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, Sun, Moon
		
			\date Dec 13, 2012, 10:58:27 AM
			\author Bartosz Lew
		*/
		cpedsPlanetDirection(string planetName="Jupiter", cpedsDirection obsLoc=cpedsDirection(18.56406,53.09546), double altitude=133, double pressure=1012, double temperature=20);
		cpedsPlanetDirection(const cpedsPlanetDirection& parent);
		cpedsPlanetDirection(const cpedsMoonDirection& parent);
		virtual ~cpedsPlanetDirection();
		
//		/*!
//			\brief returns the topocentric direction at which the planet currently is
//			\details 
//			@param JD - JD UT time to calculate direction for. for 0 - the current UT time is taken. 
//			@param epochJD - the epoch for which the returned ra-dec coordinates are specified: A value of 2000 corresponds to JD2000
//			@param geocentric - coordinates modifier. if false (default) then topocentric coordinates are returned. If true, the geocentric ra,dec are returned instead.
//			@return
//			
//			
//		
//			\date Dec 5, 2012, 1:04:37 PM
//			\author Bartosz Lew
//		*/
//		virtual DirectionRaDec nowAtRaDec(double JD=0, double epochJD=2000, double ut1_utc=0, bool geocentric=false);

//		/*!
//			\brief calculates an estimated flux density per beam of the planet
//			\details 
//			@param freq - requested frequency [Hz]
//			@param fwhm - gaussian beam size assumed for calculations [rad]
//			@param T - planet temperature [K] to be used for calculation (if not given a default value is assumed)
//			@return flux density per beam [Jy]
//			
//			The calculation assumes black body spectral radiance distribution and thermodynamic temperature for a given
//		
//			\date Jan 7, 2015, 1:52:07 PM
//			\author Bartosz Lew
//		*/
//		double fluxPerBeam(double freq, double fwhm, double T=0);
		
		const cpedsPlanetDirection& operator=(const cpedsPlanetDirection& rhs);
		const cpedsPlanetDirection& operator=(const cpedsMoonDirection& rhs);

		
		
		virtual double angularSize(double JD=0, int which=0);
	protected:
		void setPlanetProperties(string planetName);
		
};

#endif /* CPEDSPLANETDIRECTION_H_ */
