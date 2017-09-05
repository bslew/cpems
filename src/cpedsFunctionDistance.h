/*!
 * \file cpedsFunctionDistance.h
 *
 *  Created on: Oct 19, 2011
 *      Author: blew
 */

#ifndef CPEDSFUNCTIONDISTANCE_H_
#define CPEDSFUNCTIONDISTANCE_H_

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */

/* STANDALONE HEADERS */
#include "Mscs-function.h"

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */

/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
 \class 
 \brief Encapsulates calculations of distances between two sets of points defined in 2D space
 \details 
 Each set of points is given as an mscsFunction object where y values define the coordinate 
 while x values parameterize it. Therefore two functions are needed to define one point set,
 and 4 function to define two points sets.
 
 The distance between the two point sets is calculated as follows.
	 ...(finish this description)..."Do I look like a guy who's not lazy ?" :)
 
 \date Oct 19, 2011, 12:59:12 PM 
 \author Bartosz Lew
 */

class cpedsFunctionDistance {
	public:
		cpedsFunctionDistance();
		virtual ~cpedsFunctionDistance();
		
		/*!
			\brief calculates distance between two point sets
			\details 
			@param s1x - function holding the parameterized x coordinate from the first set
			@param s1y - function holding the parameterized y coordinate from the first set 
			@param s2x - function holding the parameterized x coordinate from the second set
			@param s2y - function holding the parameterized y coordinate from the second set 
			@return
		
			\date Oct 19, 2011, 2:07:40 PM
			\author Bartosz Lew
		*/
		mscsFunction distance(mscsFunction& s1x,mscsFunction& s1y, mscsFunction& s2x, mscsFunction& s2y);
		
	private:
		mscsFunction calculateDistance(mscsFunction &f1, mscsFunction &f2, mscsFunction &s2interp);
		mscsFunction _dist;
		mscsFunction _s2xinterp,_s2yinterp;
		
//		options
		bool _smallerSet;
};

#endif /* CPEDSFUNCTIONDISTANCE_H_ */
