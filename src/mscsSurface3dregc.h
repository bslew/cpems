/*!
 * \file mscsSurface3dregc.h
 *
 *  Created on: Jan 7, 2012
 *      Author: blew
 */

#ifndef MSCSSURFACE3DREGC_H_
#define MSCSSURFACE3DREGC_H_

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */

/* STANDALONE HEADERS */
#include "Mscs-function3dregc.h"
#include "mscsLine.h"
/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */

/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
 \class 
 \brief Encapsulates generation of variety of surfaces defined in 3d
 \details 
 
 \date Jan 7, 2012, 11:01:25 PM 
 \author Bartosz Lew
 */

class mscsSurface3dregc : public mscsFunction3dregc {
	public:
		mscsSurface3dregc();
		virtual ~mscsSurface3dregc();
		
		/*!
			\brief generates a paraboloid defined within given range and consisting of requested number of points
			\details 
			@param fromX - range specification of the paraboloid for X coordinate
			@param fromY - range specification of the paraboloid for Y coordinate
			@param toX - range specification of the paraboloid for X coordinate
			@param toY - range specification of the paraboloid for Y coordinate
			@param Nx - number of points along X direction
			@param Ny - number of points along Y direction
			@param A - paraboloid shape parameter in x direction
			@param B - paraboloid shape parameter in y direction
			@param vx - x coordinate of the tip of the paraboloid
			@param vy - y coordinate of the tip of the paraboloid
			@param vz - z coordinate of the tip of the paraboloid
			@param hyper - parameter allowing to create a hyperbolic-paraboloid

			The paraboloid is defined as: z = ( (x-vx)^2/A^2 + (-1)^hyper (y-vy)^2/B^2 ) + vz
			
			IMPORTANT: this function allocates only the space corresponding to a 2-d complex grid of points in correspondence to the 
			area over which the paraboloid is to be defined. The z-coodrinates of the paraboloid are stored in the real part of the
			function values. 
		
			\date Jan 7, 2012, 11:08:21 PM
			\author Bartosz Lew
		*/
		void mkParaboloid(double fromX, double toX, double fromY, double toY, long Nx, long Ny, double A=1, double B=1, double vx=0, double vy=0, double vz=0, bool hyper=false);
		
		/*!
			\brief generates a paraboloid defined within given range and consisting of requested number of points
			\details 
			@param fromX - range specification of the hyperboloid for X coordinate
			@param fromY - range specification of the hyperboloid for Y coordinate
			@param toX - range specification of the hyperboloid for X coordinate
			@param toY - range specification of the hyperboloid for Y coordinate
			@param Nx - number of points along X direction
			@param Ny - number of points along Y direction
			@param A - hyperboloid shape parameter in x direction
			@param B - hyperboloid shape parameter in y direction
			@param C - hyperboloid shape parameter in z direction
			@param vx - x coordinate of the tip of the hyperboloid
			@param vy - y coordinate of the tip of the hyperboloid
			@param vz - z coordinate of the tip of the hyperboloid
			@param hyper - parameter allowing to create a hyperbolic-paraboloid

			The hyperboloid is defined as: 
			z = C * Sqrt(1+ (x-vx)^2/A^2 + (y-vy)^2/B^2 ) - C + vz
			
			This gives the two-sheated "upper" part of the hyperboloid with the tip located at (vx,vy,vz).
			
			IMPORTANT: this function allocates only the space corresponding to a 2-d complex grid of points in correspondence to the 
			area over which the paraboloid is to be defined. The z-coodrinates of the paraboloid are stored in the real part of the
			function values. 
		
			\date Jan 8, 2012, 9:14:21 PM
			\author Bartosz Lew
		*/
		void mkHyperboloid(double fromX, double toX, double fromY, double toY, long Nx, long Ny, double A=1, double B=1, double C=1, double vx=0, double vy=0, double vz=0);


		void mkPlane(double fromX, double toX, double fromY, double toY, long Nx, long Ny, double A=1, double B=1, double C=1, double vx=0, double vy=0, double vz=0);

		
		mscsSurface3dregc& operator=(const mscsSurface3dregc& rhs);
		mscsSurface3dregc& operator=(const mscsFunction3dregc& rhs);

		mscsFunction3dregc& mkSin3D(double from, double to, double dr, double T, double phi, int dir);
		void mkGauss2D(double fromX, double toX, double fromY, double toY, long Nx, long Ny, double sx=1, double sy=1, double sz=1, double vx=0, double vy=0, double vz=0);

		/*!
			\brief generates exponential function
			\details 
			@param

			The function is of the type:
			
			f(x,y,r,base)  = base^( ((x-x0)^2+(y-y0)^2)/scale^2 ) 
			
			THe function allocates 2d sheet i.e. the size of z-dimension is 1;
		
			\date Jan 16, 2012, 10:18:07 AM
			\author Bartosz Lew
		*/
		void mkExponentialSymmetrical2D(double fromX, double toX, double fromY, double toY, long Nx, long Ny, double base=10, double scale=1, double vx=0, double vy=0);
		
		/*!
			\brief returns the normal line to the surface at the requested point
			\details 
			@param 
			@return a normal line directed downwards and hooked at P1 at the surface
		
			\date Jan 19, 2012, 1:45:20 PM
			\author Bartosz Lew
		*/
		mscsLine getNormalLine(double x, double y);


	protected:
		
	private:
		
		
};

#endif /* MSCSSURFACE3DREGC_H_ */
