/*
 * mscsFunction3d.h
 *
 *  Created on: Sep 27, 2010
 *      Author: blew
 */

/*!
\file mscsFunction3dreg.h - 
*/


#ifndef MSCSFUNCTION3DREG
#define MSCSFUNCTION3DREG

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
#include <fftw3.h>
#include "Mscs-function.h"
/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */


/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
\class mscsFunction3dreg
\brief Encapsulates a 3d function defined on a regular, equally-spaced grid.
\details This is a simple extension of the one dimensional mscsFunction for the regular 3d grid of points.
The regularity of the grid makes points storing simpler. Only the spacing in each dimension of the grid is stored and the
initial value for each dimension. The actual coordinate values for the stored points are calculated once at the initialization and stored in
linear array, and are calculated as:

(x,y,z) = (i*dx+dx/2, j*dy+dy/2, k*dz+dz/2),

where:
x,y,z - are the coordinate values
i,j,k - are the indexes of the grid cell
dx,dy,dz - are the grid separations in each dimension.

The 3d function is stored on an linear array of QPointF objects, hence the function can store complex values.

The ordering of points in the array is X-major: i.e. the index of the last dimension (Z) changes 
most frequently. - consistently with fftw conventions, hence the N'th element in the array corresponds to
i,j,k indexes in the 3d array as:\n

N=i*Ny*Nz+j*Nz+k\n

where i,j,k iterate dimensions x,y,z respectively and can vary from 0...Mx-1, 0...My-1, 0...Mz-1 respectively
where Mx,My,Mz define the size of respective dimension.

\date created: Sep 27, 2010, 8:08:43 PM 
\author Bartosz Lew
*/
class mscsFunction3dreg : public mscsFunction {
	
	
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
		mscsFunction3dreg();
		mscsFunction3dreg(string name);
		
		mscsFunction3dreg(long Nx, long Ny, long Nz, double dx,double dy, double dz);
		
		
		
		~mscsFunction3dreg();
		
		
		
		/* ---------------------------- */
		/* PUBLIC METHODS */
		/* ---------------------------- */
		
		
		
		
		/*****************************/
		/* VARIOUS FUNCTION HANDERLS */
		/*****************************/
		
		/*!
		\brief sets i,j,k'th point of the function
		\details
		If you want do substitution faster then use f(i)=... function
		\author Bartosz Lew
		*/
		void setf(long i, long j, long k, double x, double y);

		
		//!
		/*!
		\brief returns the value of an argument which is  closest to the argument x;
		\details
		
		if the function isn't allocated zero is returned
		*/
//		double f(double x, double y, double z) const;
		//! same as above but the index of the closest argument is returned under i pointer
//		double f(double x, double y, double z, long* i, long* j, long* k) const;

		//! returns the real part reference at the ijk'th position; use this if you are sure that the ijk'th value exists in the function
		double& f(long i, long j, long k) { return _f[ijk2idx(i,j,k)].rx(); }
		double& fRe(long i, long j, long k) { return _f[ijk2idx(i,j,k)].rx(); }
		//! returns the complex part reference at the ijk'th position; use this if you are sure that the ijk'th value exists in the function
		double& fIm(long i, long j, long k) { return _f[ijk2idx(i,j,k)].ry(); }

		//! returns the coordinate value at i'th position
//		double getX(long i) const { return _X[i]; }
		double getX(long i) const { return i*_dx+_dxo2; }
		//! returns the coordinate value at j'th position
//		double getY(long j) const { return _Y[j]; };
		double getY(long j) const { return j*_dy+_dyo2; }
		//! returns the coordinate value at k'th position
//		double getZ(long k) const { return _Z[k]; };
		double getZ(long k) const { return k*_dz+_dzo2; }
//		double& X(long i) { return _f[i].rx(); }
		
		//! returns a modifiable reference to the list of X coordinates list
		cpedsList<double>& getX() {return _X; }
		//! returns a modifiable reference to the list of Y coordinates list
		cpedsList<double>& getY() {return _Y; }
		//! returns a modifiable reference to the list of Z coordinates list
		cpedsList<double>& getZ() {return _Z; }
		
		//! returns the cell separation in X direction
		double getDx() const { return _dx; }
		//! returns the cell separation in Y direction
		double getDy() const { return _dy; }
		//! returns the cell separation in Z direction
		double getDz() const { return _dz; }

		long Nx() const { return _Nx; }
		long Ny() const { return _Ny; }
		long Nz() const { return _Nz; }		
		//! returns the theoretical number of points that the function should be made of, as derived from the grid size;
		long getN() { return Nx()*Ny()*Nz(); }

		void setSize(long Nx, long Ny, long Nz, double dx, double dy, double dz, double x0, double y0, double z0);
		void setSize(long Nx, long Ny, long Nz);
		
		/*!
			\brief extracts a slice of the 3D scalar function as matrix 
			\details 
			@param plane - defines the plane of the matrix. If x is associated with coordinate 0, y with 1 and z with 2, then
				plane=0 defines y-z plane, plane=1 defines x-z plane and plane=2 defines x-y plane.
			@param coord - defines the index across the plane at which to make the slice
			@param part - defines which part of the function to store (0 - real or X, and 1 - imaginary or Y)
			@return
		
			\date Oct 6, 2010, 6:27:20 PM
			\author Bartosz Lew
		*/
		matrix<double> getSlice(long plane, long coord, long part=0) const;
		
		void printInfo() const;
		
		/*****************************/
		/* FUNCTION GENERATORS		 */
		/*****************************/

		mscsFunction3dreg& mkSin3Drad(double from, double to, double dr, double T, double phi);
		mscsFunction3dreg& mkSin3D(double from, double to, double dr, double T, double phi);

		
		
		/*!
			\brief generate gaussian random field with mean m and sigma s, on the defined function domain
			\details 
			@param m - mean
			@param s - standard deviation
			@return this
		
			\date Sep 27, 2010, 9:07:45 PM
			\author Bartosz Lew
		*/
		mscsFunction3dreg& generateRandomGaussianField(double m, double s, long seed);
		
		/*!
			\brief generate gaussian random field in fourier space according to given power spectrum P(k)
			\details 
			@param Pk - power spectrum to be realized in the field
			@return this
			
			The k = sqrt(kx^2+ky^2+kz^2)
		
			\date Sep 27, 2010, 9:09:05 PM
			\author Bartosz Lew
		*/
		mscsFunction3dreg& generateRandomGaussianFField(const mscsFunction& Pk);
		
		/*!
			\brief generate gaussian random field in real space according to given Pk
			\details 
			@param
			@return
			This performs the generateRandomGaussianFField(const mscsFunction& Pk) followed by antysymmetrization and fft to the real space
		
			\date Sep 27, 2010, 9:11:15 PM
			\author Bartosz Lew
		*/
		mscsFunction3dreg& generateRandomGaussianField(const mscsFunction& Pk);

		/*****************************/
		/* FUNCTION OPERATIONS 		 */
		/*****************************/

		/*!
			\brief  derive the power spectrum of the 3D scalar real function
			\details 
			@param
			@return
		
			\date Oct 6, 2010, 3:05:09 PM
			\author Bartosz Lew
		*/
		mscsFunction powerSpectrum();
		
		
		/*!
			\brief make the function antysymmetrical in fourier space
			\details 
			@param
			@return
			
			complexConjugate(f(-(x,y,z))) = f((x,y,z))
		
			\date Sep 27, 2010, 9:13:22 PM
			\author Bartosz Lew
		*/
		mscsFunction& antysymmetrize();
		
		
		/*****************************/
		/* IO OPERATIONS 		 */
		/*****************************/
		/*!
			\brief load function data from file
			\details 
			@param
			@return
			The implementation is very similar to the load function from the mscsFunction class but allows for
			single column read. It reads the points one by one from the data file in pairs or as single values
			depending on how many colums there is in the file:\n
			one - real valued function or\n
			two - complex valued function\n
			
			Larger number of columns will result in unpredictable results. Can be extended to
			make it possible to read regardless of the number of columns in the data file and
			simply add a flag indicating whether the function is to be real or complex which
			will result in appropriate way of reading the numbers from the file: singlets for real 
			and doublets for complex.
		
			\date Oct 4, 2010, 2:58:24 PM
			\author Bartosz Lew
		*/
		virtual cpedsStatusCodes load(string filename, bool commentedFile=false);

		
		
		
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PROTECTED MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	protected:
		
		
		/* ---------------------------- */
		/* PROTECTED METHODS */
		/* ---------------------------- */
		/*!
			\brief fast Fourier transform the field
			\details 
			@param dk - address to the address of the beginning of the Fourier components of the function
			@return
		
			\date Sep 27, 2010, 9:12:19 PM
			\author Bartosz Lew
		*/
		mscsFunction3dreg& fft_r2c(fftw_complex** dk);
		
		mscsFunction3dreg& fft_c2c(fftw_complex** dk);
		
		mscsFunction3dreg slowft_r2c(const cpedsList<double>& kx, const cpedsList<double>& ky, const cpedsList<double>& kz);

		
		
		/* ---------------------------- */
		/* PROTECTED STRUCTURES */
		/* ---------------------------- */
		double _sizeX,_sizeY,_sizeZ; //!< these are the sizes of the function domain in each of the direction. They are used to calculate the cell in number in the qlist
		long _Nx,_Ny,_Nz; //!< defines the size of the regular grid in each dimension
		long _NXY,_NYZ; //!< = sizeX*sizeY, and sizeZ*sizeY respectively
		cpedsList<double> _X,_Y,_Z; // this is not needed in the "reg" class as the points spacing is assumed to be regular, however using this allows for non-regular grids to be stored as well.
		double _dx,_dy,_dz; //!< the separation of points
		double _dxo2,_dyo2,_dzo2; //!< the half-separation of points
		double _x0,_y0,_z0; //!< coordinate of the 3D box corner with smallest coordinates
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PRIVATE MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	private:
		
		
		/* ---------------------------- */
		/* PRIVATE METHODS */
		/* ---------------------------- */
		long ijk2idx(long i, long j, long k) const;
		
		/* ---------------------------- */
		/* PRIVATE STRUCTURES */
		/* ---------------------------- */
		
		
};
#endif /* MSCSFUNCTION3DREG */ 










