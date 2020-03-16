/*
 * mscsFunction3d.h
 *
 *  Created on: Dec 22, 2010, 10:43:06 PM
 *      Author: blew
 */

/*!
\file mscsFunction3dregc.h - reimplementation of the mscsFunction3dreg class but based on C arrays instead of
QLists
*/


#ifndef MSCSFUNCTION3DREGC
#define MSCSFUNCTION3DREGC

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <vector>
#include <tuple>
#ifndef NO_HDF5
#include <hdf5.h>
#endif
//#include "Mscs-object.h"
//#include "matrix.h"
#include "Mscs-function.h"
#include "cpeds-point3d.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "Mscs-map-window_function.h"
#include "mscsVector.h"
#include "cpeds-point3d.h"
#include "cpeds-point_set.h"
#include "subdomain.h"


/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */


/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
\class mscsFunction3dregcc
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
class mscsFunction3dregc : public mscsObject {
	
	
	/* ---------------------------------------------------------------------------------------------------- */
	/* CLASS PUBLIC MEMBERS */
	/* ---------------------------------------------------------------------------------------------------- */
	public:
		
		typedef struct {
			double sizeX,sizeY,sizeZ; //!< these are the sizes of the function domain in each of the direction. They are used to calculate the cell in number in the qlist
			long Nx,Ny,Nz; //!< defines the size of the regular grid in each dimension
			long NXY,NYZ,NXYZ; //!< = sizeX*sizeY, and sizeZ*sizeY respectively
			double dx,dy,dz; //!< the separation of points
			double dxo2,dyo2,dzo2; //!< the half-separation of points
			double x0,y0,z0; //!< coordinate of the 3D box corner with smallest coordinates
			double xMax,yMax,zMax; //!< coordinate of the 3D box corner with largest coordinates; this is provided explicitly to help avoid numerical round-off errors in some cases;
			
			bool derivativeXperiodic,derivativeYperiodic;
			int savePrecision;
		} function_parameters_t;
		
		/* ------------- */
		/* CLASS FRIENDS */
		/* ------------- */
		
		/* ---------------------------- */
		/* CONSTRUCTORS AND DESTRUCTORS */
		/* ---------------------------- */
		mscsFunction3dregc();
		mscsFunction3dregc(string name, cpeds_VerbosityLevel verbosityLoc=CPEDS_defaultVerbosityLevel);
		
		mscsFunction3dregc(long Nx, long Ny, long Nz, double dx=1,double dy=1, double dz=1,double x0=0, double y0=0, double z0=0, cpeds_VerbosityLevel verbosityLoc=CPEDS_defaultVerbosityLevel);
		mscsFunction3dregc(const mscsFunction3dregc& parent);
		
		
		
		virtual ~mscsFunction3dregc();
		
		void initiateParams();
		int destroyFFTPlans();
		
		/* ---------------------------- */
		/* PUBLIC METHODS */
		/* ---------------------------- */
		
		
		
		
		/*****************************/
		/* VARIOUS FUNCTION HANDERLS */
		/*****************************/
		
		/*!
		\brief sets i,j,k'th point of the function
		@param x - real part
		@param y - imaginary part
		\details
		\author Bartosz Lew
		*/
		void setf(long i, long j, long k, double x, double y);
		mscsFunction3dregc& setf(const double re,const  double im);
		void setf(long i, double x, double y) { _data[i][0]=x; _data[i][1]=y; }
		void setf(double x, double y,double z, double re, double im);
		void setRe(double value);
		void setIm(double value);
		double& fReCoord(cpedsPoint3D &p) { return fReCoord(p.x(),p.y(),p.z()); }
		double& fImCoord(cpedsPoint3D &p) { return fImCoord(p.x(),p.y(),p.z()); }
		double& fReCoord(double x, double y,double z);
		double& fImCoord(double x, double y,double z);
		double& fReCoordPeriodic(double x, double y,double z);
		double& fImCoordPeriodic(double x, double y,double z);

		/*!
			\brief convert coordinate to array index number for the current function setup and for requested axis
			\details 
			@param coord - coordinate value to be converted to array index
			@param ax - requested dimension: 0 - X, 1 - Y, 2 - Z
			@param overflow - if non NULL given then returns if crossing the function space happened
			@param idxd - pointer to an allocated variable which if different from NULL will be set to the floating point index of the cell.
			It can be used  to get information on how the overflow happened
			@return array index. May range from 0 to Nx()-1 for ax=0 and correspondingly for other axes
		
			\date Jan 4, 2013, 1:47:58 PM
			\author Bartosz Lew
		*/
		long coord2idx(double coord,long ax=0, bool* overflow=NULL, double *idxd=NULL);
		
		/*!
			\brief convert coordinate to array index number for the current function setup and for requested axis periodically
			\details 
			@param coord - coordinate value to be converted to array index
			@param ax - requested dimension: 0 - X, 1 - Y, 2 - Z
			@return array index. May range from 0 to Nx()-1 for ax=0 and correspondingly for other axes
		
			\date Jun 20, 2013, 6:33:18 PM
			\author Bartosz Lew
		*/
		long coord2idxPeriodic(double coord,long ax);

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
		fftw_complex& f(long i, long j, long k) { return _data[ijk2idx(i,j,k)]; }
		const fftw_complex& f(long i, long j, long k) const { return _data[ijk2idx(i,j,k)]; }
		double& fRe(long i, long j, long k) { return _data[ijk2idx(i,j,k)][0]; }		
		double fRe(long i, long j, long k) const { return _data[ijk2idx(i,j,k)][0]; }
		//! returns the complex part reference at the ijk'th position; use this if you are sure that the ijk'th value exists in the function
		double& fIm(long i, long j, long k) { return _data[ijk2idx(i,j,k)][1]; }
		double fIm(long i, long j, long k) const { return _data[ijk2idx(i,j,k)][1]; }

		double& fRe(long i) { return _data[i][0]; }
		double& fIm(long i) { return _data[i][1]; }
		double fRe(long i) const { return _data[i][0]; }
		double fIm(long i) const { return _data[i][1]; }
		/*!
			\brief memcopies data of the assumed size taken from this object to from src to this object.
			\details 
			@param src - pointer to the source data
		
			\date Aug 28, 2013, 10:44:48 AM
			\author Bartosz Lew
		*/
		void dataHardCopy(fftw_complex* src);
		void dataHardCopy(const fftw_complex* src);
		
		/*!
			\brief returns an interpolated function value for the requested x,y coordinates in the slice through the function at grid cell k
			\details 
			@param x - x coordinate
			@param y - y coordinate 
			@param k - z coordinate index
			@param extrapolate -- NOT USED 

			IMPORTANT:
			This function uses the internal variable controlling how the derivatives are calculated in the function field.
			Use
			non-periodic - derivatives are calculated in a non-periodic way
			
			@return interpolated function value is calculated using the bicubic interpolation 
		
			\date Jan 12, 2012, 5:19:09 PM
			\author Bartosz Lew
		*/
		double fxy(double x, double y, int k=0, int part=0, double* dfdx=NULL, double* dfdy=NULL, bool extrapolate=false);

		/*!
			\brief returns an interpolated function value for the requested x,y coordinates in the slice through the function at grid cell k
			\details 
			@param x - x coordinate
			@param y - y coordinate 
			@param k - z coordinate index
			@param extrapolate -- NOT USED 

			IMPORTANT:
			This function works periodically. 
			Consider re-implementation to use
			the internal variables to control how to go about the boundaries.
			
			@return interpolated function value is calculated using the bilinear interpolation 
		
			\date Feb 5, 2018, 4:22:07 PM
		*/
		double fxyLin(double x, double y, int k=0, int part=0, bool extrapolate=false);

		/*!
			\brief fill the holes in the real part in the slice through the function at grid cell k using imaginary part as a mask
			\details 
			@param k - z coordinate index
			@return returns *this
		
			NOT IMPLEMENTED YET
		
			\date Feb 5, 2018, 4:39:07 PM
		*/
		mscsFunction3dregc& interpolateXYholes(int k=0);
		
		mscsFunction3dregc& multiply(const mscsFunction3dregc& g);
		mscsFunction3dregc& multiply(const double v);
		mscsFunction3dregc& divide(const double v);
		/*!
			\brief divide by function f but avoiding division by zero
			\details 
			@param f dividor
			@param part - 0 = Re, 1 = Im, 2 = both
			@return self
		
			\date Aug 20, 2014, 2:26:43 PM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& divide(const mscsFunction3dregc& f, int part=0);
		mscsFunction3dregc& subtract(const double v);
		mscsFunction3dregc& add(const double v);
		mscsFunction3dregc rotateSlice(int plane, int coord, double angle, double x, double y, double z);
		/*!
			\brief center the maximal value of the function by periodic wrap
			\details 
			@param part - 0 = re, 1 = im
			@return *this
		
			\date Oct 29, 2014, 9:54:08 AM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& centerMax(bool re=true); 
		
		
		double& operator()(long i, long j);
		double operator()(long i, long j) const;
		fftw_complex& operator()(long i, long j, long k);
		fftw_complex& operator[](long i);
		const fftw_complex& operator[](long i) const;
		
		mscsFunction3dregc& operator=(const double& v);
		mscsFunction3dregc& operator=(const fftw_complex& v);
		mscsFunction3dregc& operator=(const mscsFunction3dregc& rhs);
//		const mscsFunction3dregc& operator=(const mscsFunction3dregc& rhs);

		mscsFunction3dregc& operator+=(const mscsFunction3dregc& f);
		mscsFunction3dregc& operator-=(const mscsFunction3dregc& f);
		mscsFunction3dregc& operator*=(const mscsFunction3dregc& f);
		mscsFunction3dregc& operator/=(const mscsFunction3dregc& f);
		mscsFunction3dregc& operator+=(const double v);
		mscsFunction3dregc& operator-=(const double v);
		mscsFunction3dregc& operator*=(const double v);
		mscsFunction3dregc& operator/=(const double v);

		const mscsFunction3dregc operator + (const mscsFunction3dregc& f) const { mscsFunction3dregc tmp(*this); return tmp+=f; }
		const mscsFunction3dregc operator - (const mscsFunction3dregc& f) const { mscsFunction3dregc tmp(*this); return tmp-=f; }
		const mscsFunction3dregc operator * (const mscsFunction3dregc& f) const { mscsFunction3dregc tmp(*this); return tmp*=f; }
		const mscsFunction3dregc operator / (const mscsFunction3dregc& f) const { mscsFunction3dregc tmp(*this); return tmp/=f; }


		void re2im();
		void im2re();
	
		
		//! returns the coordinate value at i'th position
//		double getX(long i) const { return _X[i]; }
		double getX(long i) const { return _param.x0+i*_param.dx+_param.dxo2; }
		long idxX(double v) const { return (v-_param.x0)/_param.dx; }
		long idxY(double v) const { return (v-_param.y0)/_param.dy; }
		long idxZ(double v) const { return (v-_param.z0)/_param.dz; }
		double cellX(double v) const { return (v-_param.x0)/_param.dx; }
		double cellY(double v) const { return (v-_param.y0)/_param.dy; }
		double cellZ(double v) const { return (v-_param.z0)/_param.dz; }
		//! returns the coordinate value at j'th position
//		double getY(long j) const { return _Y[j]; };
		double getY(long j) const { return _param.y0+j*_param.dy+_param.dyo2; }
		//! returns the coordinate value at k'th position
//		double getZ(long k) const { return _Z[k]; };
		double getZ(long k) const { return _param.z0+k*_param.dz+_param.dzo2; }
//		double& X(long i) { return _f[i].rx(); }
		
		double getMinX() const { return _param.x0; }
		double getMinY() const { return _param.y0; }
		double getMinZ() const { return _param.z0; }
		double getMaxX() const { return getX(Nx()-1)+_param.dxo2; }
		double getMaxY() const { return getY(Ny()-1)+_param.dyo2; }
		double getMaxZ() const { return getZ(Nz()-1)+_param.dzo2; }

		/*!
			\brief returns the maximal value of the function
			\details 
			@param iMax - the pointer is set to the index (in the linear array) of the maximal value element
			@param Re - if true the real part is searched, otherwise imaginary part is searched
			@return
		
			\date Mar 15, 2013, 12:53:19 PM
			\author Bartosz Lew
		*/
		double getMaxValue(long* iMax=NULL, bool Re=true) const;
		double getMinValue(long* iMin=NULL, bool Re=true) const;
		cpedsPoint3D getMinValueCell(bool Re=true) const;
		/*!
			\brief return function space coordinates of the maximal function value
			\details 
			@param Re - if true the maximal value corresponds to real part, otherwise
				it concerns the imaginary part.
			@return x,y,z point in function space of the maximal value.
		
			\date Mar 13, 2020, 12:52:29 PM
		*/
		cpedsPoint3D getMaxValueCell(bool Re=true) const;

//		//! returns a modifiable reference to the list of X coordinates list
//		cpedsList<double>& getX() {return _X; }
//		//! returns a modifiable reference to the list of Y coordinates list
//		cpedsList<double>& getY() {return _Y; }
//		//! returns a modifiable reference to the list of Z coordinates list
//		cpedsList<double>& getZ() {return _Z; }
		
		//! returns the cell separation in X direction
		double getDx() const { return _param.dx; }
		//! returns the cell separation in Y direction
		double getDy() const { return _param.dy; }
		//! returns the cell separation in Z direction
		double getDz() const { return _param.dz; }
		//! returns the cell center offset in X direction from the cell's lowest coordinate value
		double getDxo2() const { return _param.dxo2; }
		//! returns the cell center offset in Y direction from the cell's lowest coordinate value
		double getDyo2() const { return _param.dyo2; }
		//! returns the cell center offset in Z direction from the cell's lowest coordinate value
		double getDzo2() const { return _param.dzo2; }

		long Nx() const { return _param.Nx; }
		long Ny() const { return _param.Ny; }
		long Nz() const { return _param.Nz; }		
		//! returns the theoretical number of points that the function should be made of, as derived from the grid size;
		long getN() const { return _param.NXYZ; }
		long size() const { return Nx()*Ny()*Nz(); }
		double lengthX() const { return _param.sizeX; }
		double lengthY() const { return _param.sizeY; }
		double lengthZ() const { return _param.sizeZ; }
		/*!
			\brief i,j,k to array index converter in i- then j- major ordering
			\details 
			@param
			@return
		
			\date Nov 3, 2015, 5:27:52 PM
			\author Bartosz Lew
		*/
		long ijk2idx(long i, long j, long k) const;
		void idx2ijk(long idx, long& i, long& j, long& k) const;
		/*!
			\brief i,j,k to array index converter in k- then j- major ordering
			\details 
			@param i,j,k indexes defining a grid cell
			@return linear array cell index mapped using an i-minor order i.e.
			idx=i+j*Nx+k*Nx*Ny
			The actual ordering of the cells in the linear array is consistent with
			the fftw conventions: i.e. the last dimension index changes fastest.
		
			\date Nov 3, 2015, 5:27:52 PM
			\author Bartosz Lew
		*/
		long kji2idx(long i, long j, long k) const;
		
		
		/*!
			\brief calculates a volume of a given pixel
			\details 
			@param i - index along x
			@param j - index along y
			@param k - index along z
			@param spherical - if false then volume is simply dx*dy if Nz=1 and dx*dy*dz if Nz>1 and in spherical case
			this is a surface element on the sphere in case of Nz=1 and volume element of the sphere in case of Nz>1
			@return area of volume of a given pixel [srad]
		
			\date Oct 18, 2012, 9:05:54 AM
			\author Bartosz Lew
		*/
		double getPixelVolume(long i, long j, long k, bool spherical=false);
		double getVolume(bool spherical=false);

		/*!
			\brief privides the full specification on the function and grid sizes 
			\details This function does not allocate nor reallocate the memory space 
			Use realloc for this purpose.
			
			@param x0 - x coordinate of the 3D box corner with smallest coordinates
			@param y0 - y coordinate of the 3D box corner with smallest coordinates
			@param z0 - z coordinate of the 3D box corner with smallest coordinates
			
			\date Dec 23, 2010, 8:48:31 PM
			\author Bartosz Lew
		*/
		void setSize(long Nx, long Ny, long Nz, double dx, double dy, double dz, double x0, double y0, double z0);

		void setSizeRange(long Nx, long Ny, long Nz, double x0, double y0, double z0, double xMax, double yMax, double zMax);

		/*!
			\brief sets the location of the cell center in each direction
			\details 
		
			\date Jan 10, 2012, 12:27:09 PM
			\author Bartosz Lew
		*/
		void setDo2(double dxo2, double dyo2, double dzo2);
		/*!
			\brief sets the size parameters of the function space
			\details 
			This assumes that dx,dy, and dz are already defined. If not, use the another setSize function
		
			\date Dec 23, 2010, 8:47:18 PM
			\author Bartosz Lew
		*/
		void setSize(long Nx, long Ny, long Nz);
		void setSize(long Nx, long Ny);
		
		/*!
			\brief sets internal parameters that alter some calculations details
			\details 
			@param paramName - name of the parameter to alter. Could be one of:\n
				derivativeXperiodic = [0/1]
				
				derivativeYperiodic = [0/1]
				
		
			\date Jul 11, 2014, 5:28:13 PM
			\author Bartosz Lew
		*/
		void setInternalParam(string paramName, double paramVal);
		
		/*!
			\brief move the box domain 
			\details 
			@param
			@return
		
			\date Jun 21, 2013, 8:23:53 AM
			\author Bartosz Lew
		*/
		void moveBox(double dx, double dy, double dz);
		/*!
			\brief reallocate the space according to the current size settings
			\details 
			This deletes the memory and allocates the space again
		
			\date Dec 22, 2010, 10:55:50 PM
			\author Bartosz Lew
		*/
		void realloc() { freeFunctionSpace(); allocFunctionSpace(); }
		
		//! free the dynamically allocated space of the function
		void freeFunctionSpace();
		/*!
			\brief allocate the space according to the current size settings
			\details 
			@param n - number of cells to allocate; if 0 then it will allocate according to the current settings
			This does free the function memory space first and then reallocates again according to the new size.
			It does not alter the function size settings if n parameter is inconsistent with the current size settings.
		
			\date Dec 22, 2010, 10:55:50 PM
			\author Bartosz Lew
		*/
		void allocFunctionSpace(long n=0);
		
		/*!
			\brief extracts a 2D slice of the 3D scalar complex function and returns as matrix 
			\details 
			@param plane - defines the plane of the matrix. If x is associated with coordinate 0, y with 1 and z with 2, then
				plane=0 defines y-z plane, plane=1 defines x-z plane and plane=2 defines x-y plane.
			@param coord - defines the index across the plane at which to make the slice
			@param part - defines which part of the function to store (0 - real or X, and 1 - imaginary or Y, 2 - real weighted by imaginary)
			@return
		
			The output matrix contains the slice. The ordering of the data in the slice stored in the matrix is worth mentioning.
			Typically we like to think of second dimension (y) to be something vertical and first dimension (x) something horizontal on the plots.
			But in the matrix representation is it customary to have exactly the opposite: i.e. the first dimension goes across rows while the second
			across columns. 
			
			Here we stick to the visual interpretation of rows and columns so, 
			for the case of plane=2 (the x-y slice) the matrix rows hold the pixels along X axis while columns hold pixels along Y axis.

			For the case of plane=1 (the x-z slice) the matrix rows hold the pixels along x axis and cols hold the pixels along z axis.

			Finally for the plane=0 case (the y-z slice) the matrix rows hold the pixels along y axis and columns hold the pixels along z axis.
		
			\date Oct 6, 2010, 6:27:20 PM
			\author Bartosz Lew
		*/
		matrix<double> getSlice(long plane, long coord, long part=0) const;

		/*!
			\brief extracts 1D slice of the 3D scalar complex function and returns as vector
			\details 
			@param dir - defines the direction along which the slice will be taken. If x is associated with coordinate 0, y with 1 and z with 2, then
				dir indicates these directions respectively
			@param coord1 - defines the index of the first coordinate in the plane perpendicular to the direction dir
			@param coord2 - defines the index of the second coordinate in the plane perpendicular to the direction dir
			@param part - defines which part of the function to store 
				(0 - real or X, and 1 - imaginary or Y, and 2 - function arguments)
			@return pointer to the newly allocated array of doubles holding the requested part of the data
			
			For example:\n 
			if dir=0 (slice along X), then coord1 and coord2 will define the y and z coordinates on the perpendicular plane respectively.
			
			if dir=1 (slice along Y), then coord1 and coord2 will define the z and x coordinates on the perpendicular plane respectively.
			
			if dir=2 (slice along Z), then coord1 and coord2 will define the x and y coordinates on the perpendicular plane respectively.
		
			\date Nov 11, 2011, 3:04:29 PM
			\author Bartosz Lew
		*/
		double* getSlice1D(long dir, long coord1, long coord2, long part=0) const;
		mscsFunction getSlice1Dfn(long dir, long coord1, long coord2, long part=0) const;
		
		/*!
			\brief inserts the XY function X and Y values into two neighboring columns of this function starting at location i0,j0,k0
			\details 
			@param
			@return returns this object

			The function should have allocated appropriate size in each dimension
		
			\date Jan 28, 2015, 11:01:07 AM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& insertFunctionXY(mscsFunction& XY, long i0, long j0, long k0=0);
		
		/*!
			\brief inserts the list values into column of this function at location i0,j0,k0
			\details 
			@param
			@return returns this object

			The function should have allocated appropriate size in each dimension
		
			\date Jan 28, 2015, 11:30:35 AM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& insertList(cpedsList<double>& l, long i0, long j0, long k0=0);
		
		/*!
			\brief take LOS 1D slice through the box
			\details 
			@param lx,ly,lz - define the LOS vector in the function coordinates
			@param x0,y0,z0 - define the starting point
			@param dt - step when moving along LOS
			@param tmax - maximal value of LOS parameter
			@param zmax - optional pointer to zmax value; If not NULL, then LOS vector will be cut upon reaching tmax or zmax whichever comes first
			@param options - behaviour modifiers. Valid keys are:\n
				nearest - nearest cell values are taken - this is currencly the only supported method\n
				\n
			@return
		
			THIS IS NOT IMPLEMENTED YET
			\date Jun 20, 2013, 6:08:56 PM
			\author Bartosz Lew
		*/
		mscsVector<double> getLOSvec(double lx,double ly,double lz,double x0, double y0, double z0, double dt, double tmax, double* zmax=NULL, mscsVector<double>* X=NULL, mscsVector<double>* Y=NULL, mscsVector<double>* Z=NULL, string options="nearest");
		mscsVector<double> getLOSvecZmax(double lx,double ly,double lz,double x0, double y0, double z0, double dt, double Zmax, mscsVector<double>* X=NULL, mscsVector<double>* Y=NULL, mscsVector<double>* Z=NULL, string options="nearest");
		mscsVector<double> getLOSvec(double lx,double ly,double lz,double x0, double y0, double z0, double dt, double Lmax, mscsVector<double>* X=NULL, mscsVector<double>* Y=NULL, mscsVector<double>* Z=NULL, string options="nearest");

		/*!
			\brief cuts away a block of cells from [st,en] (inclusively)
			\details 
			@param iSt - cut away start cell id (X axis)
			@param jSt - cut away start cell id (Y axis)
			@param kSt - cut away start cell id (Z axis)
			@param iEn - cut away end cell id (X axis)
			@param jEn - cut away end cell id (Y axis)
			@param kEn - cut away end cell id (Z axis)
			@return a new cut-away function with domain defined based on the cut-away domain range from the original function
		
			\date Jan 15, 2016, 9:55:29 AM
			\author Bartosz Lew
		*/
		mscsFunction3dregc cutAwayBlock(long iSt, long jSt, long kSt,long iEn, long jEn, long kEn);
		mscsFunction3dregc pasteAdd(mscsFunction3dregc& block, long iSt, long jSt, long kSt);
//		const mscsFunction3dregc get_cutAwayBlock(long iSt, long jSt, long kSt,long iEn, long jEn, long kEn);


		/*!
			\brief averages the function values on planes perpendicular to requested coordinate direction
			\details 
			@param coord - 0: averages Y-Z plane, 1: averages Z-X plane, 2: averages X-Y plane
			@return mscsFunction with averaged values. The arguments are set according to the function domain in the coord direction
		
			\date Apr 9, 2013, 9:57:21 PM
			\author Bartosz Lew
		*/
		mscsFunction averagePlane(long coord);
		
		void printInfo() const;
		void printFunction() const;
		
		mscsFunction3dregc& fft(bool fwd=true);
		
		bool intersects(mscsFunction3dregc& f);
		
		/*!
			\brief set value at grid cells that intersect with function f
			\details 
			@param
			@return
		
			\date Oct 14, 2013, 4:31:29 PM
			\author Bartosz Lew
		*/
		void intersection_set(mscsFunction3dregc& f, double value, int part=0);

		/*****************************/
		/* FUNCTION GENERATORS		 */
		/*****************************/
		/*!
			\brief populates the initialized function space with points as defined by points set
			\details 
			@param ps - point set object with the position and values data
			@param ZasVals - if true then it is assumed that the function should be 2D array and only x,y coordinates from the point set are used and z coordinate is used 
			as a function value. If false - then all x,y,z coordinates are used to populate 3D array and the values are taken from the points set values structure.
			@param initialValueRe - initial real part function value for the array cells where no points are available
			@param initialValueIm - initial imaginary part function value for the array cells where no points are available
			@param redefineFunctionDomain - if false then the current function domain will be used for populating the filed with ps data and the additional option overflow is used to further adjust 
			the way the populating is done. If true, then the function space will be redefined anew to fit the volume defined by the ps. In this case overflow parameter is not used.
			@param overflow - if false (default), then no points will be allowed to fall outside of the defined function domain and will be snapped to the nearest
			function cell. If true - then the points that do not fit into the defined function domain will be ignored.
			@param how - defines how the field is populated. Accepted values are:\n 
				mean - will calculate the mean value of the points that fall into the same grid cell (default)\n
				min - will take the minimal value of the points that fall into the same grid cell\n
				max - will take the maximal value of the points that fall into the same grid cell\n
				last - will reset the cell value each time the cell contains the input point. As a result the cell value is defined by the last point that fell into any given cell.
			@return returns this object
			The imaginary part will hold the map of counts of values in different grid cells
		
			\date Jan 4, 2013, 2:33:17 PM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& populateField(cpedsPointSet3D &ps, bool ZasVals=true, double initialValueRe=0, double initialValueIm=0, bool redefineFunctionDomain=false, bool overflow=false, string how="mean");
		
		mscsFunction3dregc& mk3DRegularGrid(double fromX, double toX, double fromY, double toY, double fromZ, double toZ, double dx, double dy, double dz, double vre, double vim);
		/*!
			\brief generates a density field out of provided set of particle positions and their masses (weights)
			\details 
			@param r - definition of the domain of the function after the operation
			@param positions - a vector of positions of the partiles (points)
			@param weights - alternatively masses of the particles (if NULL, equal weights will be used)
			@param int part - 0 - store the field on the real part; 1 - store the field on imaginary part of the function
			@param smoothingKernel - string defining the SPH smoothing kernel to be used for calculating density field.  currently only "gadget2" is valid.
			@param NeighboursMin - minimal number of neighbors for the smoothing
			Allowed values are: "gadget2"
			@return returns this
		
			This function implies that the function space will be completely reallocated.
			If not allocated, it will be allocated anew based on the provided region specifications.
			In the region specification the filds with number of subdivisions should defined the number of cells in the 
			domain along each of the dimension.
			
			This function uses octTree algorithm in domain decompositioning for N log(N) complexity.
			
			The smoothingKernel should be defined in terms of r/h
			This method works both for 3D and for 2D cases. In order to recognize 2D case dz should be 0.
			
			This density calculation proceeds according to the "gather" scheme: i.e. the smoothing length is calculated for every grid point
			so that there is at least Nmin particles within hsml. Then the density is calculated according to the weights defined by the kernel
			defined with this smoothing length.
			
			The imaginary part stores the calculated hsml for each grid point.
		
		
		
			\date Jan 27, 2012, 3:38:00 PM
			\author Bartosz Lew
		*/
//		mscsFunction3dregc& mkDensityField(subDomain_region_t r, const mscsVector<cpedsPoint3D>& positions, const mscsVector<cpedsPoint3D>* weights=NULL, int part=0, mscsFunction& smoothingKernel=mscsWindowFunction().mkSPHkernelGadget2());
		mscsFunction3dregc& mkDensityFieldGather(subDomain_region_t r, subDomain_region_t treeScheme, const mscsVector<cpedsPoint3D>& positions, const mscsVector<double>* weights=NULL, int part=0, string smKernel="gadget2", long NeighborsMin=30);
		/*!
			\brief generates a density field out of provided set of particle positions and their masses (weights)
			\details 
			@param r - region that defines the function space and resolution. If you want 2D case, then set zmin==zmax and subz=1. 
			In 2D case all positions must have the same z coordinate.
			@param treeScheme - defines the tree type and domain. The domain must be defined. It is used for making tree. For 2D case use zmin==zmax
			@param positions - vector of data points positions
			@param weights - particle masses (or weights to be applied for density calculation)
			@param smKernel - defines which smoothing kernel should be used - default is "gadget2"
			@param NeighborsMin - minimal number of neighbors
			@param NeighborsMax - maximal number of neighbors
			@return returns this object
			
			Here the density calculation is done according to the "scatter" scheme: i.e. smoothing lengths are calculated only for the given particles
			so that there is between [Nmin,Nmax] neighboring particles within hsml and then the density at the grid point is obtained by summing
			contributions from individual particles according to the kernel defined by smoothing lengths of those particles and their distances to the 
			grid point.

			The density field is generated on the real part of the function. Im part contains the effective number of particles that contributed
			to the density at particular grid cell.
		
			\date Feb 2, 2012, 2:02:13 AM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& mkDensityFieldScatter2(subDomain_region_t r, subDomain_region_t treeScheme, const mscsVector<cpedsPoint3D>& positions, const mscsVector<double>* weights=NULL, string smKernel="gadget2", long NeighborsMin=30, long NeighborsMax=33);
		//! depreciated routine for doing the same as mkDensityFieldScatter2. This will be reimplemented by calling mkDensityFieldScatter2. NOT fully tested
		mscsFunction3dregc& mkDensityFieldScatter(subDomain_region_t r, subDomain_region_t treeScheme, const mscsVector<cpedsPoint3D>& positions, const mscsVector<double>* weights=NULL, int part=0, string smKernel="gadget2", long NeighborsMin=30);

		/*!
			\brief create a density field out of set of points with requested number of points in the output grid
			\details 
			@param positions - point set with the positions of the points from which the density is to be calculated
			@param maxDepth - maximal depth of the domain decomposition (-1 - decomposes the domain without using maxDepth criteria)
			@param spherical - indicates to treat x coordinate as phi, y latitude in element surface (volume) calculations (angles are assumed to be given in radians and so the resulting number density will be in sr^-1 and density in mass_unit/sr)
			@param calcMass - this is a calculation behaviour modifier. If false then the weights of the particles located in given grid cell will be summed and multiplied by the
			density of points distribution. If true, then the weights of all points located in the grid cell will be averaged i.e.
			if the weights correspond to particle masses then option true would return average mass density per particle in given grid cell rather than mass-weighted number density.
			
			@return this function returns *this
			
			This first performs a domain decomposition with domain splits as defined in the treeScheme and provides full flexibility on 
			settings of the domain sizes and grid sizes. All other routines for calculating density in the "Cube" approach use this method.
		
			\date Aug 8, 2012, 12:11:46 PM
			
			revision: Feb 22, 2013, 9:25:43 PM
			This routine has been hugely speed up by correcting tree walking in the subdomain class. The full power of the tree algorithm
			significantly reduced the numerical complexity providing exactly the same results.
			The OLD version of this routine is still in the code with the suffix _old but is now depreciated.
			It is kept for algorithm comparizon purposes that can be done using gadget-snapshot-dump-data program
			
			TODO: rethink the way 2d case is distinguished from 3d case. It should be decided upon whether dz=0 or not. Reimplement this and check
			results. Will have to also change couple of other methods as well (and possibly programs too).
			
			\author Bartosz Lew
		*/
		mscsFunction3dregc& mkDensityFieldCube(subDomain_region_t r, subDomain_region_t treeScheme, const mscsVector<cpedsPoint3D>& positions, const mscsVector<double>* weights=NULL, long maxDepth=-1, bool spherical=false, bool calcMass=false, bool unitVolume=false);
		mscsFunction3dregc& mkDensityFieldCube_old(subDomain_region_t r, subDomain_region_t treeScheme, const mscsVector<cpedsPoint3D>& positions, const mscsVector<double>* weights=NULL, int part=0, long maxDepth=-1, bool spherical=false, bool calcMass=false, bool unitVolume=false);
		/*!
			\brief create a density field out of set of points with requested number of points in the output grid
			\details 
			@param Nx - requested number of points in the output grid along X direction
			@param Ny - requested number of points in the output grid along Y direction
			@param Nz - requested number of points in the output grid along Z direction
			@param positions - point set with the positions of the points from which the density is to be calculated
			@param maxDepth - maximal depth of the domain decomposition
			@param spherical - indicates to treat x coordinate as phi, y as polar angle (theta) in element surface (volume) calculations
			@return this function returns *this
			
			This first performs a domain decomposition with domain splits 2,2,2 (the octTree scheme) down to the requested maximal depth. The root domain has depth 0.
			Decompositioning stops at the maxDepth. Then a density field is generated out of the numers of particles located in the leafs and the number density is 
			calculated as the number of particles in a leaf divided by the volume of the leaf domain. The corresponding cell in the grid is assigned with the calculated
			density value.
			As a result the domain subdivision scheme the effective resolution of the density calculations is 2^maxDepth. which is casted on to the grid of requested
			resolution (the size function grid covers exactly the root domain size).
			
			If the values are set in the positions point set, then they are used as weights in density calculation. Otherwise equal weights are used.
		
			\date Aug 8, 2012, 12:11:46 PM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& mkDensityFieldCube(int Nx, int Ny, int Nz, const cpedsPointSet3D& positions, long maxDepth, bool is2D, bool spherical, bool calcMass, bool unitVolume=false);
		/*!
			\brief create a density field out of set of points with requested number of points in the output grid
			\details 
			@param Nx - requested number of points in the output grid along X direction
			@param Ny - requested number of points in the output grid along Y direction
			@param Nz - requested number of points in the output grid along Z direction
			@param positions - point set with the positions of the points from which the density is to be calculated
			@param spherical - indicates to treat x coordinate as phi, y as polar angle (theta) in element surface (volume) calculations
			@return this function returns *this

			The difference here from the density calculations with maxDepth parameter is that the domain decomposition is always done with subdomain scheme of 
			Nx,Ny,Nz and always down to the depth of 1. This allows for matching the size of the domain resolution to the requested grid resolution.
					
			If the values are set in the positions point set, then they are used as weights in density calculation. Otherwise equal weights are used.

			\date Aug 8, 2012, 12:11:46 PM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& mkDensityFieldCube(int Nx, int Ny, int Nz, const cpedsPointSet3D& positions, bool is2D, bool spherical, bool calcMass, bool unitVolume=false);

		
		/*!
			\brief 
			\details 
			@param r - region that defines the function space and resolution. If you want 2D case, then set zmin==zmax and subz=1. 
			In 2D case all positions must have the same z coordinate.
			@param treeScheme - defines the tree type and domain. The domain must be defined. It is used for making tree. For 2D case use zmin==zmax
			@param positions - vector of data points positions
			@param vals - field values to be interpolated on the grid cells
			@param smKernel - defines which smoothing kernel should be used - default is "gadget2"
			@param NeighborsMin - minimal number of neighbors
			@param NeighborsMax - maximal number of neighbors
			@param providedHSML - a list of smoothing lengths that can be provided (if the size is
				equal the size of positions) or requested if the allocated vector size is 0,
				in which case the vector will be modified with the calculated smoothing lengths.
				
			@return returns this object
			
			This is SPH interpolation. The interpolated field is hold on real part of the function. The imaginary part holds the effective 
			number of particles that contributed to the calculation at given grid cell as a result of the neighbors finding algorithm.
			
			
		
			\date Feb 27, 2013, 9:53:42 AM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& mkInterpolatedFieldScatter(subDomain_region_t r, subDomain_region_t treeScheme, const mscsVector<cpedsPoint3D>& positions, const mscsVector<double>& vals, string smKernel="gadget2", long NeighborsMin=30, long NeighborsMax=33, mscsVector<double>* providedHSML=NULL, double MassMin=0, double MassMax=0);

//		mscsFunction3dregc& mkInterpolatedFieldScatterMass(subDomain_region_t r, subDomain_region_t treeScheme, const mscsVector<cpedsPoint3D>& positions, const mscsVector<double>& vals, string smKernel="gadget2", double MassMin, double MassMax, mscsVector<double>* providedHSML=NULL);

//		mscsFunction3dregc& mkInterpolatedFieldScatter(subDomain_region_t r, subDomain_region_t treeScheme, const mscsVector<cpedsPoint3D>& positions, const mscsVector<double>& vals, mscsVector<double>& providedHSML, mscsVector<double>& providedDens, string smKernel="gadget2");

		/*!
			\brief 
			\details 
			@param maxHSML - defines the maximal hsml that is used to define region within which the neighboring particles are sought for. if -1 then not used and maximal hsml from the 
			@param D3 - an optional tree, that will be used for optimization reasons. It should have a different domain split scheme than D tree.
			provided hsmls is used.
			@param providedDens - points number density at positions in the same length units as positions.
			@param providedHSML - smoothing lengths of the particles in the same length units as positions.
			@return
		
			\date Oct 10, 2013, 9:25:44 AM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& mkInterpolatedFieldScatter(subDomain_region_t r, subDomain* D, const mscsVector<cpedsPoint3D>& positions, const mscsVector<double>& vals, mscsVector<double>& providedHSML, mscsVector<double>& providedDens, double maxHSML=-1, string smKernel="gadget2", subDomain* D3=NULL);

		/*!
			\brief generate 2D field on an allocated function space using 2D triangulation of input points
			\details 
			@param positions - vector of x,y positions with values stored at z coordinate
			@return *this - interpolated real function values at function grid points
			
			The function space and range need to be allocated and initialized with default value.
			If the interpolation fails the fuction value is not set.
			The imaginary part of the function stores 1 on successful interpolation and 0 otherwise.

			THe returned function has Nz()=1
		
			\date May 29, 2017, 5:57:44 PM
		*/
		mscsFunction3dregc& mkInterpolatedFieldTriangSibsonGrad2D(const mscsVector<cpedsPoint3D>& positions);
		
		mscsFunction3dregc& mkInterpolatedFieldTriangLinear2D(const mscsVector<cpedsPoint3D>& positions);
		/*!
			\brief convenience function for mkInterpolatedFieldTriangLinear2D(const mscsVector<cpedsPoint3D>& positions);
			\details 
			@param positions - values are not used. Z-coordinate is used as function values.
			@return *this
			
		
			\date May 29, 2017, 8:07:00 PM
		*/
		mscsFunction3dregc& mkInterpolatedFieldTriangLinear2D(const cpedsPointSet3D& positions);

		mscsFunction3dregc& mkSin3Drad(double from, double to, double dr, double T, double phi);
		mscsFunction3dregc& mkSin3D(double from, double to, double dr, double T, double phi);
		/*!
			\brief creates a 3-d gauss function on the defined grid. 
			\details 
			@param sx - defines the variance of the gauss in units of Lx in x direction
			@param sy - defines the variance of the gauss in units of Ly in y direction
			@param sz - defines the variance of the gauss in units of Lz in z direction
			@return this
			
			The gaussian function peak is centered i.e it is located at: Nx()/2,Ny()/2, Nz()/2
			It can be periodically shifted if need be with shiftN method.
			
			The Gauss function is
			
			exp(- (x^2/(2sx^2)) - (y^2/(2sy^2)) - (z^2/(2sz^2)) )
			
			but if the sizeZ of the function domain is 0 then the function is
			
			exp(- (x^2/(2sx^2)) - (y^2/(2sy^2)) )
		
			\date Nov 11, 2011, 2:38:29 PM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& mkGauss3D(double sx, double sy, double sz);
		
		/*!
			\brief generates top hat function of radii R in the center of the field.
			\details 
			@param R - radii of the top-hat in the same units as the box size.
			@param acc - accuracy to decide between 0 and 1 when creating top-hat function. Units [the same units as the box size].
			@return *this
		
			The imaginary part is set to 0.
		
			\date Aug 28, 2013, 12:00:43 PM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& mkTopHat3D(double R, double acc=0);
		/*!
			\brief create a 3D ball - a filled sphere
			\details 
			@param R - ball radius in the same units as the box size.
			@param acc - accuracy to decide between 0 and 1 when creating top-hat function. Units [the same units as the box size].
			@return
			
			Only the real side is populated with 1s within the ball and zeros outside. 
			The ball is placed in the center of the function domain.
		
			\date Sep 3, 2014, 6:54:54 PM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& mkBall3D(double R, double acc=0);
		mscsFunction3dregc& mkBall3Dat(double x, double y, double z, double R, double acc=0);
		/*!
			\brief generates 2D gauss function across requested slice
			\details 
			@param c1 - position of the center of gauss along first coordinate [function units]. Values from within [Xmin,Xmax]
			@param c2 - position of the center of gauss along second coordinate [function units]. Values from within [Ymin,Ymax]
			@param s1 - sigma parameter along X axis [function units]. 
			@param s1 - sigma parameter along Y axis [function units]. 
			@param amplitude - gauss amplitude
			@param plane - 2 will generate gauss on X-Y plane, 1 -- will generate gauss on Z-X plane, 0 -- will generate gauss on Y-Z plane
			@param coord - coordinate of the plane along third axis
			@return *this
		
			\date Apr 9, 2013, 2:44:47 PM
		*/
		mscsFunction3dregc& mkGauss2D(double c1, double c2, double s1, double s2, double amplitude, long plane, long coord);
		/*!
			\brief generates 2D gauss function across requested slice
			\details 
			@param s11 - 1st diagonal term of the covariance matrix (this is variance, not stdev)
			@param s22 - 2nd diagonal term of the covariance matrix (this is variance, not stdev)
			@param s12 - off-diagonal term of the covariance matrix (this is covariance)
			@return *this
		
			\date Apr 9, 2013, 2:44:47 PM
		*/
		mscsFunction3dregc& mkGauss2D(double c1, double c2, double s11, double s22, double s12, double amplitude, long plane, long coord);
		mscsFunction3dregc mkKernel(mscsFunction& Pk, string interpolationType="cspline");
		
		
		/*!
			\brief generate gaussian random field with mean m and sigma s, on the defined function domain
			\details 
			@param m - mean
			@param s - standard deviation
			@param seed - seed to use for the RNs generation (used if other than 0)
			@param seedOffset - seed offset to use for the RNs generation (used if other than 0)
			@param re - if true populate the real part of the function
			@param im - if true populate the imaginary part of the function
			@return this
		
			\date Sep 27, 2010, 9:07:45 PM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& generateRandomGaussianField(double m, double s, bool re=true, bool im=true, long seed=0, long seedOffset=0);
		
		/*!
			\brief generate gaussian random field in fourier space according to given power spectrum P(k)
			\details 
			@param Pk - power spectrum to be realized in the field
			@return this
			
			The k = sqrt(kx^2+ky^2+kz^2)
		
			\date Sep 27, 2010, 9:09:05 PM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& generateRandomGaussianField(mscsFunction& Pk, long seed=0);
		mscsFunction3dregc& generateRandomGaussianField2(mscsFunction& Pk);
		mscsFunction3dregc& generateRandomGaussianField3(mscsFunction& Pk);		
		
		
		/*!
			\brief calculate gravitational potential on a grid for a given set of points
			\details 
			@param ps - point set with distribution of masses stored as values [ some length and mass units]
			@param resolution - spatial resolution with which the potential is to be resolved [length units as in points set]
			@param softening - gravitational softening [length units as in points set]
			@return returns this
			
			The function is allocated and centered on the points set.
			If softening is given -1 then it is calculated as a diagonal of the function domain divided by the number of particles in the points set.
			
		
			\date Oct 25, 2013, 9:41:05 AM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& calculateGravitationalPotential(cpedsPointSet3D& ps, double resolution, double softening=-1);		
		
//		/*!
//			\brief generate gaussian random field in real space according to given Pk
//			\details 
//			@param
//			@return
//			This performs the generateRandomGaussianFField(const mscsFunction& Pk) followed by antysymmetrization and fft to the real space
//		
//			\date Sep 27, 2010, 9:11:15 PM
//			\author Bartosz Lew
//		*/
//		mscsFunction3dregc& generateRandomGaussianField(const mscsFunction& Pk);

		/*****************************/
		/* FUNCTION OPERATIONS 		 */
		/*****************************/

		/*!
			\brief  derive the power spectrum of the 3D scalar real function
			\details 
			@param kRes - resolution factor. If 1 then the output spectrum will be of the same resolution
			as the X resolution of this function. The output spectrum resolution is defined by
			this function resolution along X box side divided by kRes, so if kRes=0.1 the output spectrum
			is tabulated with 10x larger k resolution.
			@return returns angular averaged 3d power spectrum of this function.
		
			This routine was tested with test 34 and 35 in the test-mscsfn program.
			It returns the power spectrum amplitude for a white noise N(0,sigma) map equal to sigma^2,
			as expected. Consistent results are also obtained when white noise power spectrum is used 
			to generate a real space map.
			
			Nov 3, 2015, 11:00:57 PM
			However the definition of fft_c2c has changed: I added 1/sqrt(N) factors before forward and 
			backward transformations.
			
			
			\date Oct 6, 2010, 3:05:09 PM
			\author Bartosz Lew
		*/
		mscsFunction powerSpectrum(double kRes=1);
		
		/*!
			\brief calculates correlation function in real space
			\details 
			@param Lmin - minimal value of the correlation length
			@param Lmax - maximal value of the correlation length
			@param dr - correlation function resolution
			@return correlation function.
			The range of calculation is R=(Lmin,Lmax) with resolution dr.
			The default values for Lmin, Lmax, dr are: cell size along X, box size along X, cell size along X
		
			\date Nov 2, 2015, 2:54:05 PM
			\author Bartosz Lew
		*/
		mscsFunction correlationFunctionR(double Lmin=-1, double Lmax=-1, double dr=-1);
		
		/*!
			\brief make the function antysymmetrical in fourier space
			\details 
			@param
			@return
			
			complexConjugate(f(-(x,y,z))) = f((x,y,z))
		
			\date Sep 27, 2010, 9:13:22 PM
			\author Bartosz Lew
		*/
		void antisymmetrize();
		
		void flipXYZ();
		
		/*!
			\brief shift periodically the function
			\details 
			@param nx - shift along x
			@param ny - shift along y
			@param nz - shift along y
			@return
			
			Negative shifts are also possible
		
			\date Nov 11, 2011, 7:08:36 PM
			\author Bartosz Lew
		*/
		void shift(long nx, long ny, long nz);
		
		
		/*!
			\brief gaussian smoothing in fourier space by convolution with gaussian kernel
			\details 
			@param sxG - sigma of the gaussian kernel given in the same units as Lx
			@param syG - sigma of the gaussian kernel given in the same units as Ly
			@param szG - sigma of the gaussian kernel given in the same units as Lz
			@return returns this
		
			BLcomment (Aug 27, 2013, 12:22:22 AM): this routine is compatible with the 
		 	test3DGaussSmooth from function_fft.
		 	Both of these were tested for giving the correct amplitude and width of fwhm of the input gaussian signal that was smoothed
		 	with gaussian kernel via fourier space convolution. (see also test #10 of test-mscsfn program)
		 	
		 	This smoothing conserves the function integral as the convolution is done using a normalized Gaussian kernel.

			\date Nov 11, 2011, 8:57:36 PM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& smooth3DGauss(double sxG,double syG,double szG);
		/*!
			\brief convolve real-space function with kernel ker function (given in fourier space)
			\details 
			@param ker - convolution kernel in Fourier space. Same size as this function
			@return returns convolved copy of this function
		
			\date Oct 21, 2014, 5:43:19 PM
			\author Bartosz Lew
		*/
		mscsFunction3dregc convolve_fft(mscsFunction3dregc& ker);
		
		/*!
			\brief calculate real-space convolution of this function with g convolution kernel and return the result.
			\details 
			@param
			@return
			
			 This is slow real-space implementation for testing purposes.
			 The input kernel should be normalized and it should be centered in the field. 
			 In case of eg. gaussian kernel the gauss peak should be in the middle of the field.
		
			\date Aug 28, 2013, 11:52:45 AM
			\author Bartosz Lew
		*/
		mscsFunction3dregc convolve(mscsFunction3dregc g);

		/*!
			\brief performs an iterative Richadrson-Lucy deconvolution on the current function using provided kernel function
			\details 
			@param psf - real-space point spread function (centred in the field)
			@param iter - number of iterations to perform
			@return this runction returns *this
		
			\date Oct 21, 2014, 5:13:23 PM
			\author Bartosz Lew
		*/
		mscsFunction3dregc RLdeconvolution(mscsFunction3dregc data, mscsFunction3dregc psf, long iter);
//		/*!
//			\brief performs an iterative Richadrson-Lucy deconvolution on the current function using provided kernel function
//			\details 
//			@param psf - real-space point spread function (centred in the field)
//			@param MaxIter - number of iterations to perform
//			@param loopGain - typically (0.03, 0.3)
//			@return this function returns *this
//		
//			\date Oct 21, 2014, 5:13:23 PM
//			\author Bartosz Lew
//		*/
//		mscsFunction3dregc CLEANdevonvolution(mscsFunction3dregc data, mscsFunction3dregc psf, long MaxIter, double loopGain);

		/*!
			\brief calculates derivative of the function in requested direction
			\details 
			@param dir - specifies the type of defivative: 0 - df/dx, 1 - df/dy, 2- df/dz
			@param part - 0 - re, 1 - im
			@return the differentiated function.

			Differentiation uses the GSL gsl_deriv_central algorithm.
			\date Jan 11, 2012, 2:11:33 PM
			\author Bartosz Lew
		*/
		mscsFunction3dregc derivative(int dir, int part=0);
		
		/*!
			\brief returns partial derivatives at point i,j,k calculated in the plane xy at k'th grid point
			\details 
			@param fi - pointer to the variable where df/dx will be stored
			@param fj - pointer to the variable where df/dy  will be stored
			@param fij - pointer to the variable where d^2f/dx/dy will be stored
			@part - defines which part 0 - real, 1- imaginary should be taken for calculating derivatives

			The derivatives are done using 5-point stencil and periodic boundary conditions.
			
			You can change the default behavior using the setInternalParameter function.
		
			\date Jan 13, 2012, 12:50:10 PM
			\author Bartosz Lew
		*/
		void get_derivativesXY(int i, int j, int k, double* fi, double* fj, double* fij, int part=0);
		
		
		
		/*!
			\brief calculates integral over the whole real function volume
			\details 
			@return
		
			In case of 2D function (i.e. when Dz==0) the integral is Dx*Dy*sumRe()
			otherwise it is Dx*Dy*Dz*sumRe()
		
			\date Apr 18, 2013, 9:25:06 AM
			\author Bartosz Lew
		*/
		double integrateRe();
		/*!
			\brief integrate along selected axis
			\details 
			@param axis - 0 - integrate along X, 1 - integrate along Y, 2 - integrate along Z
			@return an integrated 2D function
		
			\date Oct 11, 2013, 12:13:56 AM
			\author Bartosz Lew
		*/
		mscsFunction3dregc integrateRe(long axis);
		mscsFunction3dregc& normalizeRe();
		
		/*!
			\brief calculates sum over the whole real function space
			\details 
			@return sum
		
			\date Apr 18, 2013, 9:25:06 AM
			\author Bartosz Lew
		*/
		double sumRe();
		double averageRe();
		double varianceRe();
		double RMSre();
		/*!
			\brief calculate absolute value and substitute it into re and im parts
			\details 
			@return returns *this
		
			\date Aug 27, 2013, 9:39:39 AM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& absoluteValue();
		
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
		
		virtual cpedsStatusCodes loadMatrix(string filename);
		
		
		/*!
			\brief load ordered data from txt file
			\details 
			@param filename - file name
			@param how - 0 - only real part is read and imaginary part is set to 0; 1- both re and im parts are read from file
			@return standard status code
			
			The structure in the txt file doens't matter. Tha data are read in the default function order which is x-major then y-major
		
			\date Sep 28, 2012, 8:45:15 PM
			\author Bartosz Lew
		*/
		virtual cpedsStatusCodes loadtxtLin(string filename, int how=0);
		/*!
			\brief save the values to a file
			\details 
			@param filename - prefix for the filenames to be generated. Suffix .re and .im will be added for real and imaginary parts respectively
			@param all - false only values are stored in the file in a single column; true - x,y,z coods are in first 3 cols. and re im values in 4th and 5th col.
			@return
		
			\date Jan 19, 2012, 3:56:52 PM
			\author Bartosz Lew
		*/
		virtual cpedsStatusCodes savetxtlin(string filename, bool all=false);
		cpedsStatusCodes save(string filename);
		/*!
			\brief saves the function data in df3 format readable by pov-ray
			\details 
			@param filename - name of the output file
			@param part - 0 - real part, 1 - imaginary part
			
			The output binary data file saves the truncated 32-bit integer values obtained from scaled float function values.
			The ordering of the data in the output file is Z-major.
		
			\date Aug 23, 2012, 8:53:23 PM
			\author Bartosz Lew
		*/
		cpedsStatusCodes saveDF3(string filename, int part=0);
		cpedsStatusCodes saveAllSlices(string dirPref, string filenamePref, int plane=0, int part=0);
		
		void setSavePrecision(int precision);
		int getSavePrecision() const;
		
		/*!
			\brief  save chosen slice as matrix
			\details 
			@param part - 0 - real part, 1 - imaginary part, 2 - real part weighted by imaginary (useful for masking)

			The slice in case of plane=2 (the x-y plane) is saved so that X coordinate is saved horizontally (same y-coordinate pixels in the same row)
			and Y coordinate is saved vertically (same x-coordinate pixels in the same column)
			TODO: change this implementation for using c-arrays rather than matrix object. It will be faster
			\date Jan 13, 2012, 2:17:06 PM
			\author Bartosz Lew
		*/
		cpedsStatusCodes saveSlice(int plane, int coord, string filename, int part=0);
#ifndef NO_HDF5
		/*!
			\brief save function into hdf5 file
			\details 
			@param hdf5 file name
			@param datasetName - name for the new dataset - the file should not contain this dataset
			@param ix - element index along X to start saving with
			@param iy - element index along Y to start saving with
			@param iz - element index along Z to start saving with
			@param nx - number of elements to save along X dimension
			@param ny - number of elements to save along Y dimension
			@param nz - number of elements to save along Z dimension
			@param part - which part to save: 0 - real , 1 -imaginary
			@return exit status code
		
		
			Important: This method as many other HDF5 related are not thread-safe. You must use omp_set_lock and omp_unset_lock pairs to 
			enforce serialized access to HDF5 internals.
		
			\date Mar 23, 2012, 4:44:30 PM
			\author Bartosz Lew
		*/
		cpedsStatusCodes saveHDF5(string filename, string datasetName, int ix=0, int iy=0, int iz=0, int nx=-1, int ny=-1, int nz=-1, int part=0);

		/*!
			\brief checks if a dataset exists in the hdf5 file
			\details 
			@param
			@return
		
			\date Dec 22, 2015, 6:20:02 PM
			\author Bartosz Lew
		*/
		bool hdf5DatasetExists(string filename, string dsetName);

		vector<string> getHDF5dataSets(string filename);
		cpedsStatusCodes loadHDF5(string filename, string datasetName, int part=0);
//#pragma omp threadprivate (HDF5_stringAttributeMaxLength)
		/*
		 * Comment: Is this correct to make threadprivate a static const ?
		 * 
		 * author: blew
		 * date: Mar 25, 2016 11:19:57 AM
		 *
		 */
		
//		omp_lock_t *omp_hdf_lock;
//		void setOMPlock(omp_lock_t **lock);
#endif
		
		/*!
			\brief imports the N values from the array a onto the real or imaginary part of the function as indicated
			by the re parameter
			\details 
			
			The values are imported according to the ordering of the data in the linear array that keeps the 3d function data.
			The ordering is C-ordering - row major ordering or (x-major ordering)
			
			IMPORTANT: The memory for this function must be correctly allocated before importing the data.
		
			\date Dec 23, 2010, 7:07:41 PM
			\author Bartosz Lew
		*/
		void importFunction(double* a, long N, bool re=true);
		void importFunction(fftw_complex* a, long N);
		void importFunction(double* re, double* im, long N, bool del=true);
		void importFunction(mscsFunction& inf);
		mscsFunction3dregc concatenate(mscsFunction3dregc& inf, int axis);
		
		
		/*!
			\brief import an YZ slice
			\details 
			@param a - Y-major array of data in the YZ slice
			@param N - size of the a array - N should be  _param.Ny x _param.Nz
			@param slice - number of slice into which import the data - ranges from 0 to _param.Nx-1

			The function data space should be allocated.
			The array a is deleted inside of the mathod.
			
			\date Jan 6, 2011, 6:37:46 PM
			\author Bartosz Lew
		*/
		void importSlice(double* a, long N, int slice, bool re=true);
//		void importSlice(cpedsList<double> a,  int slice, bool re=true);
		/*!
			\brief exports a copy of the real part of the function; 
			\details 
			@param ordering - ordering of the data in the returned array. Allowed values are: "XYZmajor" (default), "YXZmojor"
			@return array with real part of the data
			
			The default ordering is the same as the this function ordering (X-major then Y-major Z-minor)
		
			\date Aug 8, 2012, 11:47:04 AM
			\author Bartosz Lew
		*/
		double* exportRe(string ordering="XYZmajor") const;
		//! exports a copy of the imaginary part of the function; The ordering is the same as the this function ordering
		double* exportIm() const;
		//! returns a i-mojor ordered list of all grid coordinates as a vector of 3-d points
		vector<cpedsPoint3D> positions2Vec();
		/*!
			\brief export the field to a point set
			\details 
			@param ZasVals - if true, then the points z coordinate and values have the same value of the function value.
				if false, then x,y,z coordinates are transferred onto the point's x,y,z coordinates and the value is set into point's value
			@param mask - this parameter can be used to indicate which data is to be exported using the information stored on the imaginary part of the function.
			Only the grid cells having the value different than the mask value as an imaginary part will be exported. If you want all cells to be exported then the function
			imaginary part need to be preset to a value different than the mask value.
			@return
		
			\date May 7, 2014, 9:13:47 AM
			REVISION:
			 May 7, 2014, 9:38:47 AM - added mask stuff so the default behavior may change wrt the previous implementation
			\author Bartosz Lew
		*/
		cpedsPointSet3D exportAsPoints(bool ZasVals, double mask=0);

		/*!
			\brief fast Fourier transform the field
			\details 
			@return fourier components of the function
			This transform is done out-of-place. New structures are created to store the result. The input is not changed.
			\date Sep 27, 2010, 9:12:19 PM
			\author Bartosz Lew
		*/
		mscsFunction3dregc fft_r2c() const;
		
		/*!
			\brief FFT of the complex field
			\details 
			@param fwd - direction parameter: true - r to c, false - c to r
			@return this function returns this.
			
			The transform is done in place. No new structures are created.
			The forward transform is divided by the number of cells so that the fundamental mode (F[0][0][0]) corresponds to average
		
			\date Dec 23, 2010, 8:01:14 PM
			\author Bartosz Lew
		*/
		mscsFunction3dregc& fft_c2c(bool fwd=true);
		
//		mscsFunction3dregc slowft_r2c(const cpedsList<double>& kx, const cpedsList<double>& ky, const cpedsList<double>& kz);
		
		void getXrange(double* xmin, double* xmax) const;
		void getYrange(double* ymin, double* ymax) const;
		void getZrange(double* zmin, double* zmax) const;
		
		function_parameters_t getFunctionParameters() const { return _param; }
		
#ifndef NO_HDF5
		void initiatie_hdf5_params();
		void setHDF5_scalarDoubleAttribute(string fname, string dsetName, string attributeName, double value, string attributeComment="");
		void setHDF5_scalarStringAttribute(string fname, string dsetName, string attributeName, string attributeValue);
		double getHDF5_scalarDoubleAttribute(string fname, string dsetName, string attributeName, int* errCode);
		string getHDF5_stringAttribute(string fname, string dsetName, string attributeName, int* errCode);

#endif
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PROTECTED MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	protected:
		


		
		/* ---------------------------- */
		/* PROTECTED METHODS */
		/* ---------------------------- */
		
		/*!
			\brief initializator for internal variables. This function should be called by the constructor ONLY !
			\details 
			@param alloc - parameter that specifies if the function space should be allocated 
			@return
		
			\date Oct 29, 2013, 5:19:54 PM
			\author Bartosz Lew
		*/
		void initiate(bool alloc=true);
		
		const fftw_complex* data() const { return _data; }
		fftw_complex* data() { return _data; }

		void calcLengths();

		/*!
			\brief get the subDomain_region definition for the requested function grid cell
			\details 
			@param i,j,k - grid cell coordinates
			@return
		
			\date Feb 22, 2013, 2:42:39 PM
			\author Bartosz Lew
		*/
		subDomain_region_t getCellSDregion(long i, long j, long k);

		
#ifndef NO_HDF5
		double getHDF5attribute(hid_t dset);
		double getHDF5_scalarDoubleAttribute(hid_t& file, string dsetName, string attributeName, int* errCode=NULL);
		string getHDF5_stringAttribute(hid_t& file, string dsetName, string attributeName, int* errCode);
		void setHDF5_scalarDoubleAttribute(hid_t& dset, string attributeName, double value, string attributeComment="");
		void setHDF5_scalarStringAttribute(hid_t& dset, string attributeName, string attributeValue);
		bool hdf5DatasetExists(hid_t& file, string dsetName);
		long getHDF5dsetCount(hid_t& file);
		bool hdf5createGroup(hid_t& file,string linkName);
#endif		

		
		
		/* ---------------------------- */
		/* PROTECTED STRUCTURES */
		/* ---------------------------- */
		
		function_parameters_t _param;
		cpedsList<double> _X,_Y,_Z; // this is not needed in the "reg" class as the points spacing is assumed to be regular, however using this allows for non-regular grids to be stored as well.

		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PRIVATE MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	private:
		
		
		/* ---------------------------- */
		/* PRIVATE METHODS */
		/* ---------------------------- */
		mscsFunction3dregc cutAwayBlock_priv(long iSt, long jSt, long kSt,long iEn, long jEn, long kEn);
		
		
		
		/*
		  //Feb 4, 2012 9:35:44 AM - DEBUG BEGIN
		  
		  typedef struct {
				  cpedsPoint3D p0;
				  long Nmin;
				  double* hsml;
				  double acc;
				  subDomain* sD;
		  } getNeighbors_fnargs_t;
		  
		  void sDgetNeighborsThreads(void* args) {
			  args->sD->getNeighbors(args->p0,args->Nmin, args->hsml,args->acc);
			  pthread_exit(0);
		  }
		  
		  //DEBUG - END
		*/
		
//		int _cur_i, _cur_j, _cur_k;
//		double gsl_function_handlerX(double x, void * params)   { return f(,_cur_j,_cur_k); }
		
		/* ---------------------------- */
		/* PRIVATE STRUCTURES */
		/* ---------------------------- */
		fftw_complex* _data;
		fftw_plan *plan_r2c_forward;
		fftw_plan *plan_r2c_backward;
		fftw_plan *plan_c2c_forward;
		fftw_plan *plan_c2c_backward;
		cpedsRNG *_rns;
};
#endif /* MSCSFUNCTION3DREGC */ 










