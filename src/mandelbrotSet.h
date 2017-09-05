/*!
 * \file mandelbrotSet.h
 *
 *  Created on: Oct 18, 2012
 *      Author: blew
 */

#ifndef MANDELBROTSET_H_
#define MANDELBROTSET_H_

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */

/* STANDALONE HEADERS */
#include "Mscs-function3dregc.h"
/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */

/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
 \class 
 \brief Encapsulates generation of Mandelbrot set with arbitrary precision
 \details 
 
 \date Oct 18, 2012, 9:59:31 PM 
 \author Bartosz Lew
 */

class mandelbrotSet : public mscsFunction3dregc {
	public:
		typedef struct {
			double x,y;
			int isIN;
			long iter;
			//double C;
		} MSpoint;
		
		
		typedef struct {
				double R,I;
		} MScomplex;
		
		mandelbrotSet(long Nx=1000, long Ny=500, double xmin=-1, double xmax=1,double ymin=-0.5 , double ymax=0.5, long maxIterations=1000, cpeds_VerbosityLevel verbosityLoc=CPEDS_defaultVerbosityLevel);
		//		mandelbrotSet(string name, cpeds_VerbosityLevel verbosityLoc=CPEDS_defaultVerbosityLevel);
		
		
		virtual ~mandelbrotSet();
		
		void mkSet();
		
	protected:
				
		//this procedure tests a given rectangular area on the complex plane and finds out which points belong to the set and which doesn't
		void find_set();
		
		int check_point(double x, double y, long iter_num, long *color);		
		
		long _maxIterations; //!< maximal amount of iterations after which calculation is stopped
		
	private:
		
		MScomplex complex_sum(MScomplex C1, MScomplex C2);
		MScomplex complex_multiply(MScomplex C1, MScomplex C2);
		double complex_abs2(MScomplex C);		
};

#endif /* MANDELBROTSET_H_ */
