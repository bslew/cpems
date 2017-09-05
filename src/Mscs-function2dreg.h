/*
 * mscsFunction3d.h
 *
 *  Created on: Dec 21, 2010, 3:23:59 PM
 *      Author: blew
 */

/*!
\file mscsFunction2dreg.h - 
*/


#ifndef MSCSFUNCTION2DREG
#define MSCSFUNCTION2DREG

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
\class mscsFunction2dreg
\brief Encapsulates a 2d function defined on a regular, equally-spaced grid.
\details This is a simple extension of the one dimensional mscsFunction for the regular 3d grid of points.
The regularity of the grid makes points storing simpler. Only the spacing in each dimension of the grid is stored and the
initial value for each dimension. The actual coordinate values for the stored points are calculated once at the initialization and stored in
linear array, and are calculated as:

(x,y) = (i*dx+dx/2, j*dy+dy/2),

where:
x,y - are the coordinate values
i,j - are the indexes of the grid cell
dx,dy - are the grid separations in each dimension.

The 2d function is stored on an linear array of QPointF objects, hence the function can store complex values.

The ordering of points in the array is X-major: i.e. the index of the last dimension (y) changes 
most frequently. - consistently with fftw conventions, hence the N'th element in the array corresponds to
i,j indexes in the 2d array as:\n

N=i*Ny+j\n

where i,j iterate dimensions x,y respectively and can vary from 0...Mx-1, 0...My-1, 0...Mz-1 respectively
where Mx,My,Mz define the size of respective dimension.

\date created: Sep 27, 2010, 8:08:43 PM 
\author Bartosz Lew
*/
class mscsFunction2dreg : public mscsFunction {
	
	
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
		mscsFunction2dreg();
		mscsFunction2dreg(string name);
		
		mscsFunction2dreg(long Nx, long Ny, double dx,double dy);
		
		
		
		~mscsFunction2dreg();
		
		
		
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
		void setf(long i, long j, double x, double y);

		
		//!
		/*!
		\brief returns the value of an argument which is  closest to the argument x;
		\details
		
		if the function isn't allocated zero is returned
		*/
//		double f(double x, double y, double z) const;
		//! same as above but the index of the closest argument is returned under i pointer
//		double f(double x, double y, double z, long* i, long* j, long* k) const;

		//! returns the real part reference at the ij'th position; use this if you are sure that the ijk'th value exists in the function
		double& f(long i, long j) { return _f[ij2idx(i,j)].rx(); }
		double& fRe(long i, long j) { return _f[ij2idx(i,j)].rx(); }
		//! returns the complex part reference at the ij'th position; use this if you are sure that the ijk'th value exists in the function
		double& fIm(long i, long j) { return _f[ij2idx(i,j)].ry(); }

		//! returns the coordinate value at i'th position
//		double getX(long i) const { return _X[i]; }
		double getX(long i) const { return i*_dx+_dxo2; }
		//! returns the coordinate value at j'th position
//		double getY(long j) const { return _Y[j]; };
		double getY(long j) const { return j*_dy+_dyo2; }
		//! returns the coordinate value at k'th position
//		double& X(long i) { return _f[i].rx(); }
		
		//! returns a modifiable reference to the list of X coordinates list
		cpedsList<double>& getX() {return _X; }
		//! returns a modifiable reference to the list of Y coordinates list
		cpedsList<double>& getY() {return _Y; }
		
		//! returns the cell separation in X direction
		double getDx() const { return _dx; }
		//! returns the cell separation in Y direction
		double getDy() const { return _dy; }

		long Nx() const { return _Nx; }
		long Ny() const { return _Ny; }
		//! returns the theoretical number of points that the function should be made of, as derived from the grid size;
		long getN() { return Nx()*Ny(); }

		void setSize(long Nx, long Ny, double dx, double dy, double x0, double y0);
		void setSize(long Nx, long Ny);
		
		
		void printInfo() const;
		
		/*****************************/
		/* FUNCTION GENERATORS		 */
		/*****************************/

		mscsFunction2dreg& mkSin2Drad(double from, double to, double dr, double T, double phi);
		mscsFunction2dreg& mkSin2D(double from, double to, double dr, double T, double phi);
		mscsFunction2dreg& mkSin2Dx(double from, double to, double dr, double T, double phi);
		mscsFunction2dreg& mkSin2Dy(double from, double to, double dr, double T, double phi);
		
		
		/*!
			\brief generate gaussian random field with mean m and sigma s, on the defined function domain
			\details 
			@param m - mean
			@param s - standard deviation
			@return this
		
			\date Sep 27, 2010, 9:07:45 PM
			\author Bartosz Lew
		*/
		mscsFunction2dreg& generateRandomGaussianField(double m, double s, long seed);
		
		/*!
			\brief generate gaussian random field in fourier space according to given power spectrum P(k)
			\details 
			@param Pk - power spectrum to be realized in the field
			@return this
			
			The k = sqrt(kx^2+ky^2)
		
			\date Sep 27, 2010, 9:09:05 PM
			\author Bartosz Lew
		*/
		mscsFunction2dreg& generateRandomGaussianFField(const mscsFunction& Pk);
		
		/*!
			\brief generate gaussian random field in real space according to given Pk
			\details 
			@param
			@return
			This performs the generateRandomGaussianFField(const mscsFunction& Pk) followed by antysymmetrization and fft to the real space
		
			\date Sep 27, 2010, 9:11:15 PM
			\author Bartosz Lew
		*/
		mscsFunction2dreg& generateRandomGaussianField(const mscsFunction& Pk);

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
			
			complexConjugate(f(-(x,y))) = f((x,y))
		
			\date Sep 27, 2010, 9:13:22 PM
			\author Bartosz Lew
		*/
		void antysymmetrize();
		
		
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

		
		virtual cpedsStatusCodes saveAsMatrix(string filename);
		
		mscsFunction2dreg& fft(bool fwd=true);
		
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
		mscsFunction2dreg& fft_r2c(fftw_complex** dk);
		
		mscsFunction2dreg& fft_c2c(fftw_complex** dk,bool fwd=true);
		
		mscsFunction2dreg slowft_r2c(const cpedsList<double>& kx, const cpedsList<double>& ky);

		
		
		/* ---------------------------- */
		/* PROTECTED STRUCTURES */
		/* ---------------------------- */
		double _sizeX,_sizeY; //!< these are the sizes of the function domain in each of the direction. They are used to calculate the cell in number in the qlist
		long _Nx,_Ny; //!< defines the size of the regular grid in each dimension
		cpedsList<double> _X,_Y; // this is not needed in the "reg" class as the points spacing is assumed to be regular, however using this allows for non-regular grids to be stored as well.
		double _dx,_dy; //!< the separation of points
		double _dxo2,_dyo2; //!< the half-separation of points
		double _x0,_y0; //!< coordinate of the 2D box corner with smallest coordinates
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PRIVATE MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	private:
		
		
		/* ---------------------------- */
		/* PRIVATE METHODS */
		/* ---------------------------- */
		long ij2idx(long i, long j) const;
		
		/* ---------------------------- */
		/* PRIVATE STRUCTURES */
		/* ---------------------------- */
		
		
};
#endif /* MSCSFUNCTION2DREG */ 










