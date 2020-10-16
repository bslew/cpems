//#include "cpeds-templates.h"
//#include "invCDFGauss.h"

//#include <lidia/bigfloat.h>
//#include <lidia/bigint.h> // commented out for a while from compatibility reasons with MATPACK ; possible fix would be to move the lidia using procedures to a separate cpeds-extra.h/c files and not include them if not needed, when compiling something with MATPACK routines

#ifndef CPEDS_MATH_DEFS
#define CPEDS_MATH_DEFS

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
//#include <iostream>
//#include <iomanip>
#include <math.h>
//#include <sstream>
#include <complex>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <tuple>

#include "cpeds-common.h"
#include "cpeds-consts.h"
//#include "cpeds-templates.h"
//#include "cpeds-pdf.h"
#include "matrix.h"
//#include <fftw3.h>

#ifndef _NO_NAMESPACE
/* using namespace LiDIA; */
using namespace std;
using namespace math;
#define STD std
#else
#define STD
#endif




typedef double fftwComplex[2];
typedef char cpeds_filenamestr[500];
typedef const char cpeds_strarg[500];

/* gsl_rng_type* CPEDS_gsl_rng_type; */
/* gsl_rng * CPEDS_gsl_rng; */

// *********************************************************************************************
// VARIOUS HELP DEFINITIONS FOR VARIOUS FUNCTIONS
// *********************************************************************************************
/* Coefficients in rational approximations. */
static const double _invCDFGauss_a[] =
{
	-3.969683028665376e+01,
	 2.209460984245205e+02,
	-2.759285104469687e+02,
	 1.383577518672690e+02,
	-3.066479806614716e+01,
	 2.506628277459239e+00
};

static const double _invCDFGauss_b[] =
{
	-5.447609879822406e+01,
	 1.615858368580409e+02,
	-1.556989798598866e+02,
	 6.680131188771972e+01,
	-1.328068155288572e+01
};

static const double _invCDFGauss_c[] =
{
	-7.784894002430293e-03,
	-3.223964580411365e-01,
	-2.400758277161838e+00,
	-2.549732539343734e+00,
	 4.374664141464968e+00,
	 2.938163982698783e+00
};

static const double _invCDFGauss_d[] =
{
	7.784695709041462e-03,
	3.224671290700398e-01,
	2.445134137142996e+00,
	3.754408661907416e+00
};

#define _invCDFGauss_LOW 0.02425
#define _invCDFGauss_HIGH 0.97575
// *********************************************************************************************

// *********************************************************************************************
// *******************   COORDINATES CHANGE STUFF ************************************************
// *********************************************************************************************
//! Conversion from cartesian coordinates to spherical coordinates
/*! @param coord: defines which spherical coordine is to be returned (in radians): 0 - theta, 1 - phi */
/*! @param x,y,z - Cartesian coordinates */
double cpeds_cart2sph(int coord, double x, double y,double z); // a conversion from cartesian to spherical system, only for directions   // coord - 0 - theta, 1 - phi

//! Conversion from spherical coordinates to cartesian coordinates
/*! @param coord: defines which cartesian coordine is to be returned: 0 - x, 1 - y, 2 - z */
/*! @param th, phi - spherical coordinates (in radians), theta is the polar angle  */
double cpeds_sph2cart(int coord, double th, double phi); // a conversion from cartesian to spherical system, only for directions   // coord - 0 - x, 1 - y 2 - z

//! Returns the cosine of the angle between two directions given in radians in spherical coordinate system with the polar angle taken from the north pole
/*! @param th1,phi1: define the first direction */
/*! @param th2,phi2: define the 2nd direction */
double cpeds_cosang_n1n2(double th1, double phi1,double th2,double phi2);

//! Returns the angle between two directions given in radians in spherical coordinate system with the polar angle taken from the north pole
/*! @param th1,phi1: define the first direction */
/*! @param th2,phi2: define the 2nd direction */
double cpeds_ang_n1n2(double th1, double phi1,double th2,double phi2);
//! Returns the angle between two directions given in radians in spherical coordinate system with the polar angle taken from the equator, just as is defined in the cpeds_direction structure
/*! @param n1: define the first direction */
/*! @param n2: define the 2nd direction */
double cpeds_ang_n1n2(cpeds_direction n1, cpeds_direction n2);

//! Verifies whether the input direction values are legal: i.e.  \f$ \theta \in <0,\pi>\f$ and \f$ \phi \in <0,2\pi>\f$.
/*! If these directions are not legal then they are corrected to be legal. */
/*! @param th1,phi1: are the pointers to exising allocated variables */
void cpeds_check_thphi(double *th, double *phi);

//! Verifies whether the input direction values are legal: i.e.  \f$ b \in <-\pi,\pi>\f$ and \f$ l \in <0,2\pi>\f$.
/*! If these directions are not legal then they are corrected to be legal. */
/*! @param b,l: are the pointers to exising allocated variables */
void cpeds_check_bl(double *l, double *b);

//! checks the b coordinate for the correct value and returns the value from within range <-PIsnd...PIsnd>
double cpeds_check_b(double b);
/*!
	\brief checks the phi coordinate for the correct value and returns the value from within range <0...2\pi>
	\details 
	@param pointer to allocated variable; the variable is modified inside the function
	@return copy of the modified variable.

	\date Jun 18, 2020, 2:35:26 PM
*/
double cpeds_check_phi(double *phi);
double cpeds_check_phi(double phi);

void cpeds_cross_product_3d(double x1, double y1, double z1,double x2, double y2, double z2, double* x3, double* y3, double* z3);

double cpeds_dot_product_3d(cpeds_direction n1, cpeds_direction n2);
double cpeds_dot_product_3d(double x1, double y1, double z1,double x2, double y2, double z2);

void cpeds_swap(double& d1, double& d2);
void cpeds_swap(int& d1, int& d2);
void cpeds_swap(long& d1, long& d2);







// *********************************************************************************************
// *******************   RANDOM NUMBERS STUFF ************************************************
// *********************************************************************************************

unsigned int cpeds_get_devrandom_seed();

//! returns a random number taken from uniform distribution from range <min,max>
/*! For this to return different results for each run the srand(0) function should be previously invoked */
/*! The RNG is the generic linux system C libraries generator */
double cpeds_random_uniform_number(double min, double max);

//! returns a pointer to a newly allocated array with pix_num random numbers taken from a uniform distribution from range <min,max>
/*! For this to return different results for each run the srand(0) function should be previously invoked */
/*! The RNG is from the implementation in as in the GSL library */
/*! @param - fast: option is now NOT IMPLEMENTED  -- will be used to indicate a very freuqent calls to this routine (faster than once per second) to avoid taking the same seeds from the system time. */
double* cpeds_random_uniform_numbersGSL(double min, double max, long pix_num,long seed_offset=0, bool fast=false);

/*!
  \brief generates a random uniformly distributed numbers
  \details
  @param min,max - range of numbers to generate
  @param num - number of the RNs
  @param t - pointer to a pointer of doubles. the pointer pointed by this pointer will be initiated with newly allocated array of doubles of size num and filled with the random numbers generated
  @param rng - pointer to the structure holding the current state of the raussian random number generato.
  @return rng - returns the pointer to the updated structure with the state of the RNG after the numbers were drawn.

  This routine must be initiated before use with code that might look like this:
  \code
  gsl_rng_env_setup();
  const gsl_rng_type * T = gsl_rng_default;
  gsl_rng* r = gsl_rng_alloc (T);
  gsl_rng_set(r,initial_seed);
  \endcode
  After all is done the structure indicated by the r variable should be released by
  issuing
  \code
  gsl_rng_free(rng);
  \endcode

  \date 2009/08/18 21:58:22
  \author Bartosz Lew
*/
gsl_rng* cpeds_random_uniform_numbersGSL(double min, double max, long num,  double **t, gsl_rng* rng);


/*!
  \brief initiates the GSL uniform RNG and returns its state
  \details
  @param seed - seed value; if 0 then taked from time
  @param seed_offset - offset to be appended to the seed
  @return state of the RNG
  The default generator is gsl_rng_mt19937

  \date 2009/08/19 00:42:36
  \author Bartosz Lew
*/
gsl_rng* cpeds_random_uniform_numbersGSL_init(long* seed, long seed_offset=0, const gsl_rng_type* generator=gsl_rng_default);

/*!
  \brief generates a random uniformly distributed numbers
  \details
  @param min,max - range of numbers to generate
  @param num - number of the RNs
  @param rng - pointer to the initalized GLS rng
  @return t - pointer to a pointer of doubles. the pointer pointed by this pointer will be initiated with newly allocated array of doubles of size num and filled with the random numbers generated

	\details see above for more details

  \date 2009/08/18 21:58:22
  \author Bartosz Lew
*/
double* cpeds_random_uniform_numbersGSL(double min, double max, long num, gsl_rng* rng);


//! same as above by returns only one number; this helps avoid frequend new and delete operation for small number of drwas
double cpeds_random_uniform_numberGSL(double min, double max, gsl_rng* rng);

//! Returns a pointer to a newly allocated array with pix_num random numbers taken from a uniform distribution from range <min,max>
/*! For this to return different results for each run the srand(0) function should be previously invoked */
double* cpeds_random_uniform_numbers(double min, double max, long pix_num);



//! Returns a single value from gaussian distribution with mean m and variance s.
/*! @param method - defines which method should be used for the RN generation.  */
/*! Method 1 - this is totaly old slow and crappy implementation that uses inverse CDF method with linear search and no interpolation. god damn  !!!! - this should be removed or developed */
/*! Method 2 - this method gives a single gaussian variate value based on the central limit theorem. pix_num is the number of points summed for one value. */
/*! Method 3 - this method is just as in the method 2 but the uniform generator is the gsl default generator -- NOT IMPLEMENTED YET see double cpeds_random_gauss_number(double m,double s, long int pix_num, int method, gsl_rng** rng) */
/*! Method 5 - this is for the internal use only; for testing purposes only; (allowes to steer the a value to find the proper arange through m value. the real m is 0 here by default */
double cpeds_random_gauss_number(double m,double s, long int pix_num, int method);

/*!
  \brief generate gaussian random number with mean m and variance s
  \details
  @param m,s - mean and standard deviation of the GRN
  @param num - number of numbers
  @param rng - pointer pointing to the structure with the state of the gaussian random number generator
  The value of pointed by the rng will change when the numbers are drawn
  @return returns the  number from the gaussian distribution

  This routine uses the central limit theorem to generate the gaussian numbers and uses
  sum of the num numbers from the uniform distribution drawn from the GSL uniform generator

  \date 2009/08/18 23:07:12
  \author Bartosz Lew
*/
double cpeds_random_gauss_number(double m,double s, long int num, gsl_rng* rng);

/*!
  \brief generate n gaussian random numbers with mean m and variance s
  \details
  A newly allocated and populated array is returned. The rest is as in the routine cpeds_random_gauss_number(double m,double s, long int num, gsl_rng* rng)

  \date 2009/08/18 23:28:43
  \author Bartosz Lew
*/
double* cpeds_random_gauss_numbers(double m,double s, long n, long int num, gsl_rng* rng);


//! Returns a pointer to a newly allocated array with pix_num random gaussian numbers from distribution with mean m and variance s
/*! @param method - defines which method should be used for the RN generation.  */
/*! Method 1 - this is totaly old slow and crappy and slow implementation that uses inverse CDF method with linear search and no interpolation. god damn  !!!! - this should be removed or developed. Remember to initialize the srand(time(0)) */
/*! Method 2 - ivnersion of the CDF but search algorithm is improved (search is done by range division by 2 method) */
/*! remember to initialize the srand(time(0)) for this method. After finding the closest mathing point, the linear interpolation is performed for the requested point between the adjacent grid points in iversion process */
/*! Method 3 - uses the uniform gsl generator but performs inverse gaussCDF transformaiton to get gaussian distr. */
/*! Method 4 - gsl gaussian distribution */
/*! Method 5 - uses the default GSL uniform generator and central limit theorem to make gaussian numbers from 20 uniform numbers  */
/*! @param sims_fast - this is a workaround for the very frequent calls to the routine (faster then 1Hz) to avoid getting the same numbers. This in principle could also be done by saving somehow the state of the RNG - which probably would be a better idea: by default: false */
/*! @param seed_offset - the offset to add to the seed taken from the current system time based on the number of seconts that elapsed from Jan 1st 1980. (by default 0*/
double* cpeds_random_gauss_numbers(double m, double s,long int pix_num,int method, long seed_offset=0, bool sims_fast=false);


//! Assigns value to the passed pointer that is a pointer to the newly generated array th and phi of size pix_num
//! that define directions (in radians) uniformly distributed on the sphere. The directions are in the CPEDS coordinate system where the polar angle is maesured from 0 (on the NP) to PI (on the SP)
void cpeds_random_uniform_directions_on_sphere(long pix_num, double **th, double **phi );
//! Same as above but only for a part of the sphere; the angles are given (in radians) defining the region for directions generation
void cpeds_random_uniform_directions_on_sphere(long pix_num, double **th, double **phi, double fromPhi, double toPhi, double fromTh, double toTh);


//! Returns the array of size size with the gaussian CDF of variance s and mean x0
/*! The CDF is stored in a table twice as big as the "size" and the ordering goes as: 1 double number - x1, 2 double number - CDF(x1) 3 double number is x2 and so on */
double* cpeds_generate_gaussian_distribuant_function(long int size, double s, double x0);





/* ******************************************************************************** */
/* ******************   PIXELIZATION SYSTEM STUFF (HEALPIX COORDINATES ARE BY DEFALUT) *********/
/* ******************************************************************************** */
//! Routine for converting pixel number to its direction in the Healpix pixelization system
/*! ordering -- 1 - when nested ordering */
/*! ordering -- 0 - when ring ordering */
/*! ordering -- 3 - iterative ordering in rows and cols - this allowes for any Nside parameter to work */
long cpeds_pix2ang_healpix(long int nside, long int pix, double *th, double *phi, int ordering);
/**********************************************************************************/
//! Routine for converting pixel number to its direction in the Healpix pixelization system
/*! ordering -- 1 - when nested ordering */
/*! ordering -- 0 - when ring ordering */
/*! ordering -- 3 - iterative ordering in rows and cols - this allowes for any Nside parameter to work */
long cpeds_ang2pix_healpix(long int nside, long int *pix, double th, double phi, int ordering); // this works only for nested at the moment

//! Routine for converting from pixel number to x,y coordinates. For internal use only.
long cpeds_pix2xy(long int nside, long int pix, double *x, double *y, int ordering);

//! Routine for converting from pixel number to column and row number in ring ordering. For internal use only.
/*! Fills the tab array with the binary representation of the pixel number for the ring ordering and returns the col and row of the pixel in the region */
void cpeds_pix2colrow_healpix_ring(long nside,long pix, long *row, long *col);

//! Routine for converting from pixel number to column and row number. For internal use only.
/*! Fills the tab array with the binary representation of the pixel number for the nested ring ordering and returns the col and row of the pixel in the region */
void cpeds_pix2colrow_healpix_nest(long nside,long pix, long *row, long *col);

//! Returns an estimated size of the pixel in the Healpix pixelization scheme for a given nside parameter in radians
double cpeds_pix_size_healpix(long int nside);

//! For a given pixel number and nside parameter returns the region of the Healpix pixelization scheme (this is only for nested scheme I guess)
long int cpeds_get_healpix_region(long int nside, long int pix);

long int cpeds_get_healpix_pixnum_in_region(long int nside, long int pix);
long int cpeds_get_healpix_colnum_in_region(long int nside, long int pix);
long int cpeds_get_healpix_rownum_in_region(long int nside, long int pix);

//! For a given nside parameter, returns the total number of rings
long int cpeds_get_ring_num_healpix(long int nside);

//! For a given nside and pixel numer i parameters and for nested ordering, returns the ring in whch the i'th pixel lies (The rings numbering stars from 0).
long int cpeds_get_ring_healpix_nest(long int nside,long int i);

//! For a given nside and pixel numer i parameters and for ring ordering, returns the ring in whch the i'th pixel lies (The rings numbering stars from 0).
long int cpeds_get_ring_healpix_ring(long int nside,long int i);

//! Converts the nested ordered map into ring ordered map for a given nside.
void cpeds_conv_nest2ring_healpix(long int nside, double * mapNEST, double * mapRING);

//! Converts the nested ordering map of directions into ring ordering map for a given nside.
void cpeds_conv_nest2ring_healpix(long int nside, cpeds_direction * mapNEST, cpeds_direction * mapRING);

//! Converts the ring ordered map into nested ordered map for a given nside.
void cpeds_conv_ring2nest_healpix(long int nside, double * mapRING, double * mapNEST);

//! Converts the ring ordered map of directions into nested ordered map for a given nside.
void cpeds_conv_ring2nest_healpix(long int nside, cpeds_direction * mapRING, cpeds_direction * mapNEST);


long cpeds_i2j_healpix(long nside, long i, long* ringpixtab); // internal procedure

//! Returns the number of pixels in a given ring for a given nside
long int cpeds_get_pixnum_in_ring_healpix(long int nside,long int ring);

//! Returns the number of pixels in all rings for a given nside; the space needs to be deallocated manually
long* cpeds_get_pixnum_in_ring_healpix(long int nside, long* N);

/*!
	\brief Returns the colatitudes of the rings for a given nside of the Heaplix pixelization scheme;
	\details
	@param
	@return pointer to the newly allocated array

	Colatitudes in [rad] are returned in increasing order;\n
	The allocated space needs to be deallocated manually
	\date May 21, 2010, 3:34:35 PM
*/
double* cpeds_get_ring_colatitudes_healpix(long int nside, long* N);

//! Returns the number of pixels above a given ring; ring numbering starts with 0
long cpeds_get_pix_num_above_ring_healpix(long nside, long ring);

//! Returns the number of pixels in the healpix pixelization scheme for a given nside parameter
/*! \f$ N_{pix} = 12*n_s^2 \f$ */
long cpeds_get_healpix_pix_num(long nside);


/*!
	\brief get list of nested pixel IDs in ns resolution that are contained in the coarser 
	parent pixel with id=pixID in the coarse_ns resolution
	\details 
	@param ns - HP ns parameter
	@param coarse_ns - parent pixel ns parameter (must be < ns)
	@param pixID - pixel ID in nested map of coarse_ns resolution
	@return returns a pair of long values: 
	
	offset and Npix 
	
	where:
	
	 offset is the pixel id in the ns resolution of the first pixel 
	 	 that belong to the coarse pixel and 
	 	 
	 Npix is the number of pixels that belong 
	 	 to the coarse pixel. The numbering of the pixels that belong to the coarse
	 	 pixel is continuous.

	\date Jun 4, 2020, 5:03:03 PM
*/
std::pair<long,long> cpeds_get_healpix_nested_pixels(long ns, long coarse_ns, long pixID);












/**********************************************************************************/
/*******************   MATH STUFF *************************************************/
/**********************************************************************************/
//! Rotates the point p about the Ox axix by Ax angle given in radians
/*! The coordinate convention is that: */
/*! Ox: (1,0,0) */
/*! Oy: (0,1,0) */
/*! Oz: (0,0,1) */
/*! The polar angle of the sphericl coordinate system is calculated from */
/*! the xy plane and not from from the pole, Convension consistent with the Mscs package, but not with the Healpix pakage */
cpeds_direction cpeds_Rx(cpeds_direction p,double Ax);
void cpeds_Rx(double x, double y, double z,double Ax, double* xr, double* yr, double* zr);
//! Rotates the point p about the Oy axix by Ay angle given in radians
cpeds_direction cpeds_Ry(cpeds_direction p,double Ay);
void cpeds_Ry(double x, double y, double z,double Ay, double* xr, double* yr, double* zr);
//! Rotates the point p about the Oz axix by Az angle given in radians
cpeds_direction cpeds_Rz(cpeds_direction p,double Az);
void cpeds_Rz(double x, double y, double z,double Az, double* xr, double* yr, double* zr);

























// ****************************************************************************************************************
// ****************************************************************************************************************
// ****************************************************************************************************************
// ****************************************************************************************************************
// Astro stuff
// ****************************************************************************************************************
// ****************************************************************************************************************
// ****************************************************************************************************************
// ****************************************************************************************************************
//! Converts the a direction from horizontal coordinate system into a direction in the first equatorial coordinate system
/*! All angles are given and returned in radians in the same coordinate system.  */
/*! @param A,h the azimuth and elevation above the horizon in radians. The Azimuth is calculated from the north eastwards. */
/*! @param th - the latitude  of the observation point in radians */
/*! The return value is the hour angle and declination of the input direction in the equatorial CS  */
/*! The routine does not account for the precession effects. The math is: */
/*! \f$ sin d = sin(th) sin(h) + cos(th) cos(h) cos(A) \f$ */
/*! \f$ sin ha = - cos(h) sin(A)/cos(d) \f$ */
/*! \f$cos(ha)=(sin(h) cos(th)-cos(h) sin(th) cos(A))/cos(d) \f$ */
/*! where h,A are the elevation and azimuth, and */
/*! th - is the geodetic latitude and */
/*! ha,d are hour angle and declination */
cpeds_direction cpeds_horizontalToEquatorialFirst(double A, double h, double th);
cpeds_direction cpeds_horizontalToEquatorialFirst(cpeds_direction n, double th);


// ****************************************************************************************************************
//! Converts the a direction from horizontal coordinate system into a direction in the first equatorial coordinate system
/*! All angles are given and returned in radians in the same coordinate system.  */
/*! @param ha,d the hour angle and declination in radians.  */
/*! @param th - the lattitude  of the observation point in radians */
/*! The return value is the Azimuth (The Azimuth is calculated from the north eastwards) and elevation of the input direction in the horizontal CS  */
/*! The routine does not account for the precession effects. The math is: */
/*! \f$ sin h = sin(th) sin(d) + cos(th) cos(d) cos(ha) \f$ */
/*! \f$ sin A = - cos(d) sin(ha)/cos(h) \f$ */
/*! \f$cos(A)=(sin(d) cos(th)-cos(d) cos(ha) sin(th))/cos(h) \f$ */
/*! where h,A are the elevation and azimuth, and */
/*! th - is the geographical lattitude and */
/*! ha,d are hour angle and declination */
cpeds_direction cpeds_equatorialFirstToHorizontal(double ha, double d, double th);
cpeds_direction cpeds_equatorialFirstToHorizontal(cpeds_direction n, double th);

// ****************************************************************************************************************
//! Converts from first equatorial to second equatorial system.
//! @param lst - the local star time is given in radians.
//! @param n - direction given in radians
cpeds_direction cpeds_equatorialFirstToEquatorialSecond(cpeds_direction n, double lst);


// ****************************************************************************************************************
//! Converts from first equatorial to second equatorial system.
//! @param lst - the local star time is given in radians.
//! @param n - direction given in radians
cpeds_direction cpeds_equatorialSecondToEquatorialFirst(cpeds_direction n, double lst);


// ****************************************************************************************************************
/*! Returns the local sidereal time in radians
 *! Nutation is not included.
 *
 *! @param localJulianTime - JD  Julian date/time given in julian days
 *! @return  local sidereal time, in radians
 *!
 *! @version 16.02.2006
 *! @author E. Pazderski, modified by Bartosz Lew
 */
// this doesn't work correctly
double cpeds_local_sidereal_time2(double localJulianTime);



// ****************************************************************************************************************
/*! Returns the local mean sidereal time in radians
 *! Nutation is not included.
 *
 *! @param localJulianTime - the local JD Julian date/time given in julian days on the place of observation at longitude lon
 *! @param lon - observers longitude in degrees (from within \pm 180)
 *! @return  local sidereal time, in radians
 */
/*! The algorithm behind this routine is based on book */
/*!  "Practical astronomy with your calculator" */
/*!  By Peter Duffett-Smith 
 * 
 * This procedure is highly (<0.01 s) compatible with novas mean sidereal time with equonix method.
 * if localJulianTime is given with dUT1 correction.
 * The two are strongly (~2s) inconsistent with the ln_get_mean_sidereal_time of libnova.
 * 
 * */
double cpeds_local_sidereal_time(double localJulianTime, double lon);

/***************************************************************************************/
/*!
	\brief wrapper for novas version of the local sidereal time with nutation implemented
	\details 
	@param jd_ut1 - current ut1 time [JD]
	@param ut1_utc - UT1-UTC = dUT [s]
	@param DeltaAT - TAI-UTC [integer s] - number of leap seconds since the beginning of TAI (this is now provided by the astro_time server)
	@param longitude - RT32 longitude [s] = RTLON[deg]/15*3600
	@param LSTtype - 0 - mean, 1 - true
	@return mean (or apparent) local sidereal time [seconds]

	TT-TAI = 32.184 s
	
	
	\date Jun 14, 2013, 7:08:13 PM
	\author Bartosz Lew
*/
double cpeds_local_sidereal_time_novas(double jd_ut1, double ut1_utc, double DeltaAT,double longitude, int LSTtype=0 );


// ****************************************************************************************************************
/*! Converts julian date to date format
 *!
 *! @param JD - julian date in days
 */
void cpeds_JDToYMDhms(double JD,long* year, long* month, long* day, long* hour, long* minute, double* sec);

/***************************************************************************************/
/*!
	\brief obtain the string corresponding to given JD UT time
	\details 
	@param JD - julian date/time (default current UT time)
	@param fmt - format to convert to string
	@return

	WARNING: this routine is dangerous. If you truncate sedonds to an integer then
	you may get number 60 eg. which is wrong. This should be corrected.

	\date Nov 20, 2011, 3:46:21 PM
	\author Bartosz Lew
*/
string cpeds_JDToYMDhms(double JD=-1e9, string fmt="%li-%02li-%02li, %02li:%02li:%06.3lf");

/*!
	\brief convert unix time stamp to date
	\details 
	@param u - unix time stamp (assumed to be localtime in seconds)
	@return date string 

	\date Jun 4, 2018, 9:44:42 AM
*/
string cpeds_unixToYMDhms(double u=-1, string fmt="%Y-%m-%d %H:%M:%S");
/***************************************************************************************/
/*!
	\brief converts the given JD date for the UT time to the corresponding year.
	\details 
	@param JD - julian date/time (default current UT time)
	@return

	\date Nov 21, 2011, 2:45:40 PM
	\author Bartosz Lew
*/
long cpeds_JDToYear(double JD=-1e9);
///*!
//	\brief converts the given JD date for the UT time to the corresponding year and fraction of year.
//	\details 
//	@param JD - julian date/time (default current UT time)
//	@return
//
//	\date Nov 21, 2011, 2:45:40 PM
//	\author Bartosz Lew
//*/
//double cpeds_JDToYearFrac(double JD=-1e9);

// ****************************************************************************************************************
/*!   PURPOSE: */
/*!      This function will compute the Julian time for a given calendar */
/*!      date (year, month, day, hour). */

/*   REFERENCES:  */
/*      Fliegel & Van Flandern, Comm. of the ACM, Vol. 11, No. 10, October */
/*      1968, p. 657.. */

/*   VER./DATE/ */
/*   PROGRAMMER: */
/*      V1.0/06-98/JAB (USNO/AA) */

/*   NOTES: */
/*      1. This function is the "C" version of Fortran NOVAS routine */
/*      'juldat'. */
/*      2. This function makes no checks for a valid input calendar */
/*      date. */

/* @param year	Year. */
/* @param month	Month number. */
/* @param day	Day-of-month.    */
/* @param hour	Hour-of-day. BL. Can have minutes and seconds included in it. The hour is for the winter time (does not include the daylight saving time */
/* @return        Julian date. The 0 hour allways means the midday i.e. the noon*/
/* @author NOVAS adapted by Bartosz Lew */
/* @version 14.02.2006  */
/*! Definition from Wikipedia: The Julian date (JD) is the interval of time in days and fractions of a day, since January 1, 4713 BC Greenwich noon, 
 * Julian proleptic calendar.[1] In precise work, the timescale, e.g., Terrestrial Time (TT) or Universal Time (UT), should be specified.[2] 
 * 
 * This is the time in Julian callendar time as measured at longitude 0.
 * hour is in hours:from [0, 24) for lon=0
 * */
double cpeds_julian_time(long year, long month, long day, double hour);

//! same as cpeds_julian_time (long year, long month, long day, double hour) but for the current system time expressed in UTC time: i.e. for the current UTC time
double cpeds_julian_time();

/*!
	\brief converts time given in seconds since 1970 to Julian date for UTC time.
	\details 
	@param s - seconds of utc time elapsed since Jan 1st 1970
	@return

	\date Dec 17, 2012, 11:29:49 PM
	\author Bartosz Lew
*/
double cpeds_timeSec_to_julian_time(double s);

//! Returns the local julian time for a given gregorian date and hour (at lon=0) and given longitude of the observer
//! THe longitude is given in degrees and can range from -180 (W) to 180 (E)
//! hour is in hours:from [0, 24)
double cpeds_julian_local_time(long year, long month, long day, double hour, double longitude);



//! Returns the modified julian date (time) given juian date (time).
//! Modified JD has time arrangenent to start a new day at midnight not midday.
double cpeds_modified_julian_time (double JD);

/*!
  \brief converts from JD local mean solar time, to zonal time.
  \details
  @param JD - local JD time (i.e. UTC corrected for the longitude - DST state should not affect given JD to this routine. DST exists only in h,m,s space. In other words we do not consider such a construction as a "julian local DST")
  @param timeZone - time zone (ranges from -12...+12)
  @param longitude - geographic longitude in deg
  @param dst - indicates whether the local time should be DST or not.
  @return The date/time values are returned via corresponding pointers.

  This accounts for possible change due to DST if needed.
  \date 2010/04/14 20:02:57
*/
void cpeds_JDToLocalZonalTime(double JD,long timeZone, double longitude, bool dst, long* YY, long* MM, long* DD, long* h, long* m, double* s);

/*!
  \brief converts from zonal civil time to local JD time (in julian days)
  \details
  @param JD - pointer to returned local JD time (i.e. UTC corrected for the longitude - DST state should not affect given JD to this routine. DST exists only in h,m,s space. In other words we do not consider such a construction as a "julian local DST")
  @param timeZone - time zone (ranges from -12...+12)
  @param longitude - geographic longitude in deg
  @param dst - indicates whether the local time should be DST or not.
  @param YY - year
  @param MM - month
  @param DD - day
  @param h - year
  @param m - minute
  @param s - second

  This accounts for possible change due to DST if needed.
  \date 2010/04/29 18:10:40
*/
double cpeds_ZonalTimeToLocalJD(long YY, long MM, long DD, long h, long m, double s, long timeZone, double longitude, bool dst);

// ****************************************************************************************************************
//! Converts an angle in degrees to degrees,minutes,seconds convention
void cpeds_angToDMS(double ang, double *d, double *m, double *s, double acc=1e7);

//! Converts an angle in degrees to hours,minutes,seconds convention
void cpeds_angToHMS(double ang, double *h, double *m, double *s);


/*!
	\brief Converts an angle in degrees minutes and seconds to a double valued angle in degrees
	\details 
	@param d - degrees [deg]
	@param m - minutes [']
	@param s - seconds ['']
	@return angle corresponding to d,m,s
	
	The convention for the negative coordinates is that if d!=0 then the sign for the hemisphere should be carried by m, then if m=0
	then the sign for the hemisphere is indicated by s.
	
	if d<0 then m and s are assumed to be positive. If they are negative then they are converted to positive numbers so that -1 30 0 means
	-1.5 deg below equator and -1 -30 0 means the same thing.

	\date Jun 14, 2012, 1:12:27 AM
	\author Bartosz Lew
*/
double cpeds_DMSToAng(double d, double m, double s);

//! Converts an hour angle in hours, minutes and seconds to a double valued angle in degrees
double cpeds_HMSToAng(double h, double m, double s);

//double cpeds_HMSstrToAng(string HHMMSS,string separator=":");

/*!
  \brief derives azimuth of the star at horizon from which the star will reach the zenith.
  \details
  @param lat - local latitude [rad]
  @return - azimuth [rad]

  The formula for the aziumuth of the rising star that will cross local zenith only makes sense for latitudes from +-45 deg because
  for all higher latitudes the stars that cross the local zenith are circum-polar so they don't rise.
  The formula is:
  A = arccos(tan(lat))

  \date 2009/10/28 14:19:44
  \author Bartosz Lew
*/
double cpeds_zenithStarAzimuthAtRising(double lat);

/*!
  \brief derives the angle at which a star will rise w.r.t horizon at azimuth A for observer at latitude lat.
  \details
  @param A - azimuth [rad]
  @param lat - latitude [rad]
  @return angle [rad]
  The formula is:
  alpha = ArcCos[ Sin[lat] / Sin[ arccos[cos(A) * cos(lat) ] ] ]

  The azimuth is calculated from the north eastwards.

  \date 2009/10/28 14:24:35
  \author Bartosz Lew
*/
double cpeds_angleToHorizonAtRising(double A, double lat);

/*!
	\brief calculates length of a section through a spherical shell
	\details 
	@param zd - zenith angle [deg]
	@param R0 - radius of the innermost shell - radius of the Earth [km]
	@param Ri - radius of the i'th shell base [km]
	@param hi - height of the i'th shell above its base [km]
	@return distance from the shell base to the i+1 shell base as looking at the angle zd. [km]

	todo: it'd be nice to implement the refraction corrections to this formula depending on 
	the given temperature and RH profiles.
	
	\date Apr 17, 2015, 12:17:49 AM
	\author Bartosz Lew
*/
double cpeds_nonplanar_atmospheric_layer_path(double zd, double Ri, double hi, double R0=6371);

/*!
	\brief calculates air mass for requested elevation
	\details 
	@param elev - elevation [deg]
	@param fit - UNUSED at the moment, leave the default value
	@return air mass - a ratio of how much more atmosphere there is at elevation elev than there is towards zenith.
	
	The formula is a fit to the analytical solution of the air mass for Earth atmosphere taken from Tools of radio astronomy. page 186 (book page 175). 
	Comes from paper Schoenberg 1929)
	
	The fitting formula gives reasonable results within range of elevations from ~2.7 to 90 deg. Below ~3 deg the fit doesn't make sense anymore as it has maximum there.
	Error is reported to be less than: 6.4 × 10−4 for fitting formula.
	
	\date Sep 13, 2012, 12:26:10 PM
	\author Bartosz Lew
*/
double cpeds_air_mass(double elev, bool fit=true);

/*!
	\brief converts temperature given in Kelvins to temperature in eV
	\details 
	@param temp - temperature in [K]
	@return temperature in eV

	this is simply conversion as:
	temp [eV] = temp [K] * CPEDS_kB/CPEDS_e

	\date Apr 25, 2013, 3:47:59 PM
	\author Bartosz Lew
*/
double cpeds_K2eV(double temp);
double cpeds_K2eV();
double cpeds_K2keV();
double cpeds_K2keV(double temp);


/*!
	\brief calculates dew point from given relative humidity and temperature
	\details 
	@param T - temperature [Celsius]
	@param RH - relative humidity [%]
	@param alpha - fitting parameter [hPa]
	@param beta - fitting parameter
	@param lambda - fitting parameter [degC]
	@return dew temperature in Celsius
	
	Based on document:
	http://irtfweb.ifa.hawaii.edu/~tcs3/tcs3/Misc/Dewpoint_Calculation_Humidity_Sensor_E.pdf
	see also the references therein.
	
	The accuracy of this formula is to within 0.3 degC within temperature range -50 - 50 degC and RH 0 -- 100 %

	\date Oct 1, 2012, 10:53:53 AM
	\author Bartosz Lew
*/
double cpeds_dew_point(double T, double RH, double alpha=6.112, double beta=17.62, double lambda=243.12 );

/*!
	\brief calculate relative humidity from T and Tdew by inverting Magnus approximation
	\details 
	@param T - temperature of the mixture of gases [degC]
	@param Tdew - dew point temperature [degC]
	@param alpha - fitting parameter [hPa]
	@param beta - fitting parameter
	@param lambda - fitting parameter [degC]
	@return

	\date Nov 5, 2014, 9:16:06 AM
	\author Bartosz Lew
*/
double cpeds_h2o_RH(double T, double Tdew, double alpha=6.112, double beta=17.62, double lambda=243.12 );
/*!
	\brief calculates saturation pressure for water vapor at temperature T
	\details 
	@param T - temp [K]
	@return saturation pressure [mbar]

	H20 saturation pressure for abundance 1.0 of HITRAN isotopologues of H2O.

	From AM program of S.Paine
   Computes H2O saturation vapor pressure over a flat liquid
   surface at temperature T.  This is the customary meteorological
   definition for reporting RH, even at sub-freezing temperature,
   and for vapor in equilibrium with small droplets.  The
   formulation used in this function is  Eq. 10 of

     D. M. Murphy and T. Koop 2005, "Review of the vapour
     pressures of ice and supercooled water for atmospheric
     applications."  Q. J. R. Meteorol. Soc. 131:1539.

   which has a stated range of applicability of 123 K < T < 332 K.

   Standard isotopic composition is assumed, and changes in
   isotopic fraction between phases are ignored.  This is an
   approximation which can be off by several percent in nature.
   See for example

     M. Kakiuchi and S. Matsuo 1979, "Direct measurements of
     D/H and 18O/16O fractionation factors between vapor and
     liquid water in the temperature range from 10 to 40 C."
     Geochemical Journal 13:307.
	
	\date Jun 26, 2013, 3:11:52 PM
	\author Bartosz Lew
*/
double cpeds_h2o_Psat(double T);
/*!
	\brief converts relative humidity to vmr for water vapor
	\details 
	@param RH - relative humidity [%]
	@param T - temperature [C]
	@param P - atmospheric pressure [mbar]
	@return vmr - volume mixing ratio

	\date Jun 26, 2013, 3:11:48 PM
	\author Bartosz Lew
*/
double cpeds_h2o_RH2vmr(double RH, double P,double T);


double cpeds_h2o_RH(double T, double P, double vmr);

/*!
	\brief generates atmospheric pressure profile
	\details 
	@param z - altitude [m]
	@param P0 - pressure at the reference level at altutude z0 [hPa]
	@param T0 - temperature at the reference level z0 [K]
	@param z0 - reference altutude [m]
	@param Lr - lapse rate [K/m]
	@param g0 - gravitational constant [m/s^2]
	@param mu - molar mass of the atmosphere [kg/mol]
	@return returns this

	\date Oct 4, 2012, 11:32:36 PM
	\author Bartosz Lew
*/
double cpeds_pressure_vs_altitude(double z, double P0, double T0, double z0, double Lr=-6.5, double g0=9.80665, double mu=0.0289644);

/*!
	\brief calculate factor to convert from flux density to thermodynamic temperature
	\details 
	@param freq - frequency [Hz]
	@param T0 - reference temperature for which approximation works
	@return conversion factor [K/Jy]

	\date Apr 16, 2013, 3:04:39 PM
	\author Bartosz Lew
*/
double cpeds_radiance_to_temp_conversionFactor(double freq, double T0=2.726);

/*!
	\brief calculate back body radiance at given temperature and frequnecy
	\details 
	@param freq - frequency [Hz]
	@param T0 - black body temperature [K]
	@return 2 h nu^2/c^2  *  h nu /(exp(x)-1) / Jy  [Jy/sr]
	
	where 
	
	x=h*nu/(kB*T0)

	\date Jan 29, 2014, 3:11:50 PM
	\author Bartosz Lew
*/
double cpeds_black_body_radiance_Bnu(double freq, double T0=2.726);

/*!
	\brief calculates TSZE frequency factor
	\details 
	@param freq - frequency [Hz]
	@param T0 - cmb temperature today [K]
	@return

	\date Nov 8, 2013, 10:03:00 PM
	\author Bartosz Lew
*/
double cpeds_TSZEgnu_factor(double freq, double T0);

/*!
	\brief 
	\details 
	@param ZDobs - observed (refracted) zenith distance [deg]
	@param alt - altitude above sea level [m]
	@param T - temperature [C]
	@param P - pressure [mbar]
	@param hum - relative humidity [%]
	@param lambda - wavelength [cm] 
	@param lat - latitude of the observer [deg]
	@param Tlapse temperature lapse in the troposphere [K/m] (suggested value: 0.0065)
	@param acc - accuracy to terminate the the iteration [rad] (suggested value 1e-8)
	@return returns the unrefracted (in space) zenith distance [deg]

	\date Jun 3, 2013, 3:23:05 PM
	\author Bartosz Lew
*/
double cpeds_refraction(double ZDobs, double alt, double T, double P, double H, double lambda, double lat, double Tlapse, double acc);

/*!
	\brief calculate refraction for the ZDspace direction
	\details 
	@param ZDspace - true (in space) zenith distance [deg]
	@param alt - altitude above sea level [m]
	@param T - temperature [C]
	@param P - pressure [mbar]
	@param hum - relative humidity [%]
	@param lambda - wavelength [cm] 
	@param lat - latitude of the observer [deg]
	@param Tlapse temperature lapse in the troposphere [K/m] (suggested value: 0.0065)
	@param acc - accuracy to terminate the the iteration [rad] (suggested value 1e-8)
			Value passed to cpeds_refraction.
	@return returns the refracted (observed) zenith distance [deg]
	
	The accuracy of this routine is as follows:
	for ZDobs=80 deg error is ~3.3e-05 deg
	for ZDobs=85 deg error is ~0.00026 deg
	for ZDobs=89 deg error is ~0.0012 deg
	

	\date Apr 29, 2020, 5:54:53 PM
*/
double cpeds_refraction_space(double ZDspace, double alt, double T, double P, double H, double lambda, double lat, double Tlapse, double acc);

extern "C" {
	extern void* slarefro_(double *zobs, double* alt,double* T, double* P, double* hum, double* lambda, double* lat, double* Tlapse, double* acc, double* ref);

/*!
	\brief radial velocity of RT32 at given direction and time
	\details 
	Input quantities:
		dj, rar,decr
	Outout quantities:
		aLMST, Vsun, Vobs, Vdop

	Jednostki:            d rad   rad  rad    km/s  km/s km/s

	Oblicza rzut predkosci, Vtot [km/s], teleskopu RT32 (Piwnice
	k. Torunia) na kierunek obserwowanego zrodla o wspolrzednych
	    ra, dec - rektascensja i deklinacja [rad] przeprecesowane 
	na zadany moment czasu:
	    dj - data julianska (UTC)
	    aLMST - sredni czas gwiazdowy dla polozenia RT32 [rad]
	    Vsun - predkosc Slonca (ku apeksowi) jest uwzgledniona w Vtot
	dlatego, jesli nie chcemy jej uwzgledniac, trzeba uzyc Vtot-Vsun.
	Maksymalna roznica Vtot w stosunku do obliczen scislych z
	wykorzystaniem numerycznych efemeryd JPL w latach 1990 - 2019
	jest mniejsza niz 0.57 m/s

	\date Jan 2, 2015, 3:24:46 PM
	\author KB, included by Bartosz Lew
*/
	extern void* cpeds_vdrt32_(double* dj,double* rar,double* decr,double* aLMST,double* Vsun,double* Vobs,double* Vdop);
}


/*!
	\brief this refraction is taken from the RT32 control system for tests which way it works
	\details 
	@param h - elevation [rad]
	@param iObs - parameter defining if the refraction is to be added or removed
	@return elevation [rad]

	\date Jun 5, 2013, 9:14:40 AM
	\author Bartosz Lew
*/
double cpeds_refractionKB(double h, int iObs);


double cpeds_RT4_arcsin(double x);

/*!
	\brief calculate theoretical beam fwhm for cassegrain telescope of given aperture at requested frequency freq.
	\details 
	@param aperture - [m]
	@param freq - frequency [GHz]
	@return beam FWHM [deg]

	The returned angular resolution quoted as the full beam width at half maximum of the beam power pattern
	is suitable for the telescope with primary 10 times bigger than the secondary mirror and 
	with gaussian secondary illumination with 12 dB taper at the edges.

	\date Nov 5, 2013, 1:20:10 PM
	\author Bartosz Lew
	
	Apr 27, 2014, 12:17:20 AM - REVISION
		updated the fitting formula value to 1.80102 due to corrections in the theoretical beam pattern derivations
		provided by KB.
*/
double cpeds_calculateAngularResolution(double aperture, double freq);


typedef struct {
	int day; 			/**< current day*/
	int month;			/**< current month*/
	int year;			/**< current year*/
	int hour;			/**< current hour*/
	int min;				/**< current min*/
	int sec;				/**< current sec*/
	double usec;		/**< current micro sec*/
	int dayofyear;		/**< day of year*/
	double JD;			/**< Julian Date*/
	double JDhour;		/**< Julian Data with hours*/
	double ST;			/**<local sidereal time [sec] */
	int Day_of_Year;	/**< number of days since 1994.01.01 (GP )*/
	double UT;			/*<< UT1 time as a number of seconds since midnight*/
	double date;		/*<< date as double: 1e6*day+1e4*month+year*/

} rt4_time_date;

rt4_time_date cpeds_RT4_cal_date(double tjd);

double cpeds_rt4_JulianDay(int tye,int tmo,int tda);
void cpeds_rt4_nutation_and_aberration(double DJ,double *dRA,double *dDEC);
/*!
	\brief nutation as implemented in RT4 control system by KB
	\details 
	@param DJ - jd_tt
	@param dRa - input mean RA [deg] for the epoch DJ
	@param dDec - input mean DEC [deg] for the epoch DJ
	@param dPsi - nutation in longitude ["]
	@param dEps - nutation in obliquity ["]
	@param eObliquity - ecliptic obliquity ["]
	@return
	
	The output true ra,dec coordinates: dRa and dDec are also in degrees

	\date Jun 13, 2016, 11:36:15 AM
	\author Bartosz Lew
*/
void cpeds_rt4_nutation(double DJ,double *dRA,double *dDEC,double* dPsi=NULL, double* dEps=NULL,double* eObliquity=NULL);
/*!
	\brief calculate RA precession of mean equatorial coordinates from juda0 to juda
	\details 
	@param alpha0 [deg]
	@param delta0 [deg]
	@return ra processed to epoch juda [deg]

	\date Jun 13, 2016, 1:45:52 PM
	\author Bartosz Lew
*/
double RT4_PrecessionRA(double juda0,double juda,double alpha0,double delta0);
double RT4_PrecessionDEC(double juda0,double juda,double alpha0,double delta0);

// ****************************************************************************************************************
// ****************************************************************************************************************
// ****************************************************************************************************************
// ****************************************************************************************************************
// ********************* MATH STUFF ******************
// ****************************************************************************************************************
// ****************************************************************************************************************
// ****************************************************************************************************************
// ****************************************************************************************************************

//! returns true if d>=0  and false if d<0
bool cpeds_positive(double d);
//! Returns the index of the element in 2D array of elements x,y
/*! y can only by 0 or 1 -- this is conversion from linear table to rectangular table of size nx2 */
long int  cpeds_xy2num(long int x, long int y);

//! returns 1 if the number is the correct real number, and 0 if the number if NAN
int cpeds_isnan (double x);

//! Returns the base 2 logarithm  of the argument
double cpeds_log2(double x);

//! Returns the factorial of the integer argument n. The result is a double number
double cpeds_factorial(long int n);

//! Same as above but the return type is long double
long double cpeds_factorial_ld(long int n);

/* bigint cpeds_factorial_lidia(long int n); */

//! Returns the ratio n!/k! of the factorial values of the arguments: i.e. k*(k+1)*...*n (k<=n)
double cpeds_factorial_frac(long int n, long k);

//! Same as aobove but for the long double type
long double cpeds_factorial_frac_ld(long int n, long k);

/* bigint cpeds_factorial_frac_lidia(long int n, long k); */
//int cpeds_find_minimal_value(double *a, long int size);
//int cpeds_find_minimal_value(float *a, long int size);

//! Finds the minimal and maximal values in the input array a of size size and assigns values of max, max, mini, maxi
/*! @param - mini and maxi are the indexes of the min and max values in the array */
void cpeds_find_minmax_value(long *a, long int size, long *min, long *max, long int *mini, long int *maxi);

//! Finds the minimal and maximal values in the input array a of size size and assigns values of max, max, mini, maxi
/*! @param - mini and maxi are the indexes of the min and max values in the array */
void cpeds_find_minmax_value(double *a, long int size, double *min, double *max, long int *mini, long int *maxi);

//! Same as above but for the float type.
void cpeds_find_minmax_value(float *a, long int size, float *min, float *max, long int *mini, long int *maxi);

//! Same as above but for the complex numbers; minimal and maximal numbers are found by their absolute values
void cpeds_find_minmax_value(complex<double> *a, long int size, complex<double> *min, complex<double> *max, long int *mini, long int *maxi);

//! Returns the last maximal value in the array a of size size starting from index st; The maxi variable is assigned to the index of that point in the array.
double cpeds_find_max_value(double *a, long int size, long st, long int *maxi=NULL);

//! Same as above but for the float type.
long cpeds_find_max_value(long *a, long int size, long st, long int *maxi);

//! Returns the last maximal value in the cpeds_pixel structure of size size, found in the first field starting from st index. The maxi variable is assigned to the index of that point in the array.
cpeds_point cpeds_find_max_value(cpeds_point *a, long int size, long st, long int *maxi);

/*!
	\brief Find the minimal value in the array a of doubles of size size starting from index st; 
	\details 
	@param a array pointer
	@param size - total size of the array
	@param st - index the start with
	@return minimal value within indexes range [st,size-1]
	
	The maxi variable is assigned to the index of that point in the array.
	
	IMPORTANT: note that the st parameter has a different meaning from the one used
	in cpeds_mean_value and similar routines.

	\date Jan 26, 2016, 1:46:10 PM
	\author Bartosz Lew
*/
double cpeds_find_min_value(double *a, long int size, long st, long int *mini=NULL);

//! Returns the last minimal value in the array a of longs of size size starting from index st; The mini variable is assigned to the index of that point in the array.
long cpeds_find_min_value(long *a, long int size, long st, long int *mini);

//! Returns the last minimal value in the cpeds_pixel structure of size size, found in the first field starting from st index. The mini variable is assigned to the index of that point in the array.
cpeds_point cpeds_find_min_value(cpeds_point *a, long int size, long st, long int *mini);


//! Returns the larger of the two values
long cpeds_get_max(long v1, long v2);

//! Returns the larger of the two values and indicates which one is bigger
/*! @param i - is set to 1 when v1 is bigger and to 2 when v2 is bigger */
long cpeds_get_max(long v1, long v2, long* i);

//! Returns the smaller of the two values
long cpeds_get_min(long v1, long v2);

//! Returns the smaller of the two values and indicates which one is smaller
/*! @param i - is set to 1 when v1 is smaller and to 2 when v2 is smaller */
long cpeds_get_min(long v1, long v2, long* i);

//! Returns the smaller of the two double values
double cpeds_get_min(double v1, double v2);

//! Returns the bigger of the two double values
double cpeds_get_max(double v1, double v2);

//! Returns the bigger of the two double values
/*! @param i - is set to 1 when v1 is smaller and to 2 when v2 is smaller */
double cpeds_get_max(double v1, double v2, long* i);

//! Returns the bigger of the two double values
/*! @param i - is set to 1 when v1 is smaller and to 2 when v2 is smaller */
double cpeds_get_min(double v1, double v2, long* i);

//! Returns the sum of the values in the array tab of long values of size size
long cpeds_sum(long * tab,long int size, bool deleteTab=false);

//! Returns the sum of the values in the array tab of double values of size size
double cpeds_sum(double * tab,long int size, bool deleteTab=false);
long double cpeds_sum(long double * tab,long int size, bool deleteTab=false);
complex<double> cpeds_sum(complex<double> * tab,long int size, bool deleteTab=false);



/*!
	\brief Returns the mean value of the values in an array
	\details 
	@param tab - array pointer
	@param size - number of cells to be averaged
	@param startFrom - starting index
	@return mean value calculated on tab range [startFrom,startFrom+size)

	It is user's responsibility to avoid running outside of the tab.
	In case of startFrom=0, size can be interpreted as the total size of tab array
	(if the whole array is the be used).

	\date Jan 26, 2016, 4:10:41 PM
	\author Bartosz Lew
*/
double cpeds_mean_value(const double * tab,long size, long startFrom=0);
double cpeds_median_value(const double * tab,long size);
complex<double> cpeds_mean_value(complex<double> * tab,long int size);


//! Returns the central moment value of the values in the tab array of size size
/*! @param moment - defines which central moment to calculate; The value should be integer 
 *  @param size - the size of the data on which to calculate the moment, this doesn't have to be the size of the array
 *  @param startFrom - index of the element to start the calculation with
 * */
double cpeds_central_moment(double * tab,long int size, double moment, long startFrom=0);
complex<double> cpeds_central_moment(complex<double> * tab,long int size, double moment);

//! Returns the unbiased estimator of the variance value of the values in the tab array of size size;
double cpeds_variance(double * tab, long int size, long startFrom=0);

//! returns the covariance betewwn two vectors of length n.
double cpeds_covariance(const double* d1,const double* d2,long n);

//! Returns the root mean square;
double cpeds_rms(double * tab, long int size);
complex<double> cpeds_rms(complex<double> * tab, long int size);

//! Returns the skewness value of the values in the tab array of size size
double cpeds_skewness(double * tab, long int size);

//! Returns the kurtosis value of the values in the tab array of size size
double cpeds_kurtosis(double * tab, long int size);

//! Returns the value of the N[m,s](x-x0) function i.e. \f$  \frac{A}{\sqrt{2\pi}*s}*exp(-\frac{(x-m)^2}{2*s^2}) \f$
double cpeds_gaussN(double x, double A, double m, double s);

//! I'm not sure about the significance or validity of this routine. Check the sources. Basically there's no such thing as Gauss on sphere...
double cpeds_gauss_on_sphere(double thc, double phic, double s, double th, double phi);

//! I'm not sure about the significance or validity of this routine. Check the sources. Basically there's no such thing as Gauss on sphere...
double cpeds_normalize_gauss_on_sphere(double s, long int pix_num);

//! Returns the \f$ \chi^2 \f$ value for the data array of k measurements each of which done with statistical error given in the array err.
/*!  The model used for the \f$ \chi^2 \f$ calculation is the gaussian curve as defined by overall aplitude A, standard deviation sigma0 and the mean mean0 */
double cpeds_test_chisq(long int k, double *data, double *err, double A, double mean0, double sigma0);

//! Something about the pearson test ... see the code for more details
double cpeds_Pearson_test(long int k, double *t, double A, double mean0, double sigma0, double bin, double * bin_num_real, double * Cstat);

//! True if the given chisq value of a chi-squared test with degr_of_freedom degrees of freedom was successful at the confidence level 1-alph
bool cpeds_chisq_accepted(double alph, double degr_of_freedom, double chisq);

//! True if the given chisq value of a chi-squared test with degr_of_freedom degrees of freedom was successful at the confidence level 1-alph; Qx is set to the actual p-value in the test.
/*! So the test is successful if the Qx (the p-value i.e.  the probability to accept) is greater then alpha. */
bool cpeds_chisq_accepted(double alph, double degr_of_freedom, double chisq, double *Qx);

int cpeds_print_chisq_confidence(double alph, double degr_of_freedom, double chisq);
double cpeds_print_chisq_accepted_confidence(double degr_of_freedom, double chisq);
//! this doesn't exist
double cpeds_chisq_prob(double x);
double cpeds_calculate_chisq(double * tab1, double * tab2, long int num);


//! tabulates the normalized but scaled gauss function with given normalization A, sigma s, mean m, in num point within range min < x < max; the return value is the poiner to the table
double* tabulate_gauss_function(int f, double A, double s, double m, long int num, double min, double max);


//! Returns the number of counts of data points in the array t of size k that fall into a range <x,x+bin)
/*! method : 1 - do this without sorting the array */
/*! method : 2 - do this upon sorting the data first, using different method */
long int cpeds_count_numbers_in_bin(long int k, double * t, double bin, double x, int method);


/*!
  \brief bins the data in he t array of size k into num equal bins which values are returned in the xbin array
  \details
  @param k -size of the t array
  @param t - array of numbers to bin
  @param num - number of bins
  @param xbin - pointer to the allocated array of size num; It will be filled by this routine with the bin values
  @param align - alignment of the bins.<br>
  -1 - the xbin values indicate the left edge of the bin - this is default<Br>
  0 -  the xbin values indicate the center of the bin<br>
  1 -  the xbin values indicate the right edge of the bin<br>
  @param from - if < 'to' then num bins will be defined within [from,to] range
  @return Returns pointer to a newly allocated array of binned input array of data t of size num.


  The space for xbin must be preallocated.<br>
  In case of -1, the first bin starts with the smallest value in t array and is of size bin=(max-min)/num: i.e <min,bin)
  The next bin is  <bin,2bin) and so on. The returned xbin values are min,bin,2bin...<br>
  The last bin is right inclusive: <(num-1)*bin,num*bin>

  This routine is heavily used in the histogram program

  \date 2009/08/28 15:43:04
  \author Bartosz Lew
*/
double * cpeds_bin_data(long int k, double * t, long num, double * xbin, long align=-1, double from=0, double to=0, double geometricMultiplier=1.0);

double* cpeds_bin_function(const double* xin, double* yin, const double* w, long Nin, long* binSize, long Nbins, double** xout, long* Nout, long imin);


///*!
//	\brief bin data in bins of size increasing geometrically
//	\details 
//	@param k -size of the t array
//	@param t - array of numbers to bin
//	@param from - start binning from this index
//	@param to - end binning on this index
//	@param binSize - initial bin size
//	@param gm - geometrical sequence multiplier
//	@param xbin - array containing derived bins. It is always allocated inside of this function
//	@param n - size of the xbin array
//	@return
//
//	\date Jan 31, 2011, 8:48:06 AM
//	\author Bartosz Lew
//*/
//double* cpeds_bin_data_geo(long k, double* t, long from, long to, double binSize=3, double gm=1.01, cpedsList<double>& xbin );

//! Same as above but for floats
float * cpeds_bin_data_float(long int k, double * t, long num, float * xbin);

//! copy the input array and return a pointer to it
double* cpeds_copy_array(const double* t, long N);

/*!
  \brief - adds a double value to an array
  \details
  @param t - double array
  @param N - size of the array
  @param val - value to be added
*/
void cpeds_add_value(double val, double* t, long N);
void cpeds_add_value(double val, fftwComplex* t, long N);

void cpeds_sub_value(double val, double* t, long N);
void cpeds_sub_value(double val, fftwComplex* t, long N);

void cpeds_mul_value(double val, double* t, long N);
void cpeds_mul_value(double val, fftwComplex* t, long N);

void cpeds_divide_value(double val, double* t, long N);
void cpeds_divide_value(double val, fftwComplex* t, long N);

//template <class Type> Type *cpeds_bin_data(double *t,long int k, Type *xbin, long num );
//#include "cpeds-templates.c"
//double *cpeds_bin_data(double *t,long int k, double *xbin, long num );
//float *cpeds_bin_data(double *t,long int k, float *xbin, long num );

//! Sorts the t array of data of size k
/*! @param direction - 12 - increasing order; 21 - decreasing order */
void cpeds_sort_data(long int k, double * t, int direction);

//! Sorts the t array of data of size k
/*! @param direction - 12 - increasing order; 21 - decreasing order */
void cpeds_sort_data(long k, long * t, int direction);

//! Sorts the t array of data of size k according to the first field in the point structure
/*! @param direction - 12 - increasing order; 21 - decreasing order */
void cpeds_sort_data(long int k, cpeds_point * t, int direction);


/*!
	\brief convert fwhm to sigma for a gaussian distribution
	\details 
	@param fwhm - arbitrary units
	@return sigma - same units as fwhm
	
	This function returns fwhm/(2 sqrt(2 ln(2)))

	\date Feb 4, 2014, 9:29:01 AM
	\author Bartosz Lew
*/
double cpeds_fwhm2sigma(double fwhm);


/*!
	\brief 2D data smoothing in Fourier space
	\details 
	@param t - data array in row-major ordering
	@param N - total size of the t array
	@param n1 - size of the last dimension (consistent with fftw naming conventions)
	@param sl0 - smoothing length in first direction defined as fwhm in units of cells count
	@param sl1 - smoothing length in second direction defined as fwhm in units of cells count

	The smoothed field is returned on the same array t.\n
	The smoothing is done via two successive fourier trnasforms in each direction followed by
	convolution with a gaussian kernel.
	
	\date Jan 12, 2011, 12:19:51 PM
	\author Bartosz Lew
*/
void cpeds_smoothGauss_data2D(double* data, long N, long n1, double sl0, double sl1);

void cpeds_smoothGauss_data3D(double* data, long n0, long n1, long n2, double sl0, double sl1);

/*!
	\brief fft real 3D data along selected axis.
	\details 
	@param in - linear row-major real array representing 3D array of size n0*n1*n2
	@param out - linear row-major comples array representing 3D array of size n0*n1*(n2/2+1)
	@param n0 - size of the first dimention
	@param n1 - size of the second dimention
	@param n2 - size of the third dimention	
	@param dim - fft direction indicator: 0 - along n0, 1 - along n1, 2 - along n2
	
	Depending on the dim value the corresponding dimension will be reduced down to n/2+1
	after fft.\n
	The fft complex values are stored on out array.

	\date Jan 18, 2011, 10:58:19 AM
	\author Bartosz Lew
*/
void cpeds_fft_data3D(double* in, fftwComplex* out, long n0, long n1, long n2, int dim);

/*!
	\brief fft of complex 3D data along selected axis.
	\details 
	@param in - linear row-major complex array representing 3D array of size n0*n1*(n2/2+1)
	@param out - linear row-major real array representing 3D array of size n0*n1*n2
	@param n0 - size of the first dimention; 
	@param n1 - size of the second dimention
	@param n2 - size of the third dimention	
	@param dim - fft direction indicator: 0 - along n0, 1 - along n1, 2 - along n2
	
	The n0,n1,n2 - always refer to the size of the array in the real space. Eg. for
	real array of size 2 x 3 x 4 the input n parameters for this routine should be
	2,3 4.
	Depending on the dim value the corresponding dimension will be reduced down to n/2+1
	after fft.\n
	The fft complex values are stored on out array.

	\date Jan 18, 2011, 10:58:19 AM
	\author Bartosz Lew
*/
void cpeds_fft_data3D(fftwComplex* in, double* out, long n0, long n1, long n2, int dim);

/*!
	\brief convolve the 3d data with gaussian smoothing kernel in fourier space
	\details 
	@param data - the array of data consistent with the format obtained from 
	cpeds_fft_data3D(double* in, fftwComplex* out, long n0, long n1, long n2, int dim) routine.
	@param fwhm - fwhm of the gaussian kernel.
	@param scale - if true than the data will be scaled by n where n is the number of cells along the dimention dim.
	This is to make sure that after inverse fft the data will be correctly scalled as fftw doesn't do that.
	If scale is false then the scaling is not performed.

	The convolution is done in place.
	
	\date Jan 18, 2011, 11:06:08 AM
	\author Bartosz Lew
*/
void cpeds_convolve_gauss(fftwComplex* data, long n0, long n1, long n2, int dim, double fwhm, bool scale);

/*!
	\brief convolve the 3d data with gaussian smoothing kernel in fourier space
	\details 
	@param data - the array of data consistent with the format of the input data passed to this
	cpeds_fft_data3D(double* in, fftwComplex* out, long n0, long n1, long n2, int dim) routine.
	@param fwhm - fwhm of the gaussian kernel.
	@param scale - if true than the data will be scaled by n where n is the number of cells along the dimention dim.
	This is to make sure that after inverse fft the data will be correctly scalled as fftw doesn't do that.
	If scale is false then the scaling is not performed.

	The convolution is done in place.
	The inpud data are first transformed into the fourier space, multiplied with the gaussian kernel and then transformed back
	to the real space. The convolution is done only along dimension indicated by dim. See info on this parameter from
	the cpeds_fft_data3D routine.
	
	\date Jan 18, 2011, 11:06:08 AM
	\author Bartosz Lew
*/
void cpeds_convolve_gauss(double* data, long n0, long n1, long n2, int dim, double fwhm, bool scale);

//void cpeds_convolve_data3D(double* data, long n0, long n1, long n2, double* Pk);


//! Finds the value closest to the value val in the array t of double values of size ts searching in num values from array index starting at start.
/*! Assumes that the array is an array of values SORTED  in rising order. */
/*! The seartch is performed with half division method */
/*! Returns the array index of the closest value */
long cpeds_find_value(double val,double * t,long ts, long start, long num);
long cpeds_find_value(long val,long * t,long ts, long start, long num);

//! Checks if the value val exists in the array t of size size
bool cpeds_val_in_array(double val, double* t, long size);

//! Returns an array of coordinates of the list of the points given in the Pxy array of size N
/*! @param coord: 0 - returns the x coordinates; 1 - returns the y coordinates */
double* cpeds_point2array(long N, cpeds_point *Pxy, long coord);

//! Saves the matrix of elements given in Dvec array. The matrix is of size vec_size x vec_num and it is stored in the vector_num major format.
/*! This means that the vecors v1, v2, ... are stored as: v11, v12, ..., v1vec_size, v21,...v2vec_size,...,vvec_num1,...,vvec_numvec_size */
/*! The output file (named file) format is in a single column format  */
void cpeds_save_matrix_lin(double *Dvec, long vec_size, long vec_num, string file, bool deleteDvec=false, bool append=false);

//! Prints out the Dvec array as a matrix matrix of size vec_size x vec_num. The ordering of the array Dvec is vec_num major
void cpeds_print_matrix(double *Dvec, long vec_size, long vec_num);

//! Prints out the matrix object to the screen
void cpeds_print_matrix(const matrix<double>& m);

//! Saves the M array as a matrix in a rows x cols format to a file with name name.
/*! The ordering of the M array is rows-major, meaning that the adjacent elements in the array will be placed
 *! in the neighboring columns in the file within the same row.  */
void cpeds_save_matrix(double * M, long rows, long cols, string name, bool deleteM=false,bool append=false);
void cpeds_save_matrix(float * M, long rows, long cols, string name, bool deleteM=false,bool append=false);

/*!
	\brief  save the complex data to file into two column (re, im) format
	\details 
	@param M array of complex data
	@param rows - size of the array
	@param name - file name
	@param deleteM - flag to destroy M after saving

	\date Nov 11, 2011, 10:05:27 PM
	\author Bartosz Lew
*/
void cpeds_save(fftwComplex* M, long rows, string name, bool deleteM=false);


/*!
  \brief shift matrix elements by Nrows and Ncols
  \details
  @param
  @return returns the reference to the input matrix

  \date 2009/11/20 12:56:36
  \author Bartosz Lew
*/
matrix<double>& cpeds_shift_matrix(matrix<double>& m, long Nrows, long Ncols);

/*!
  \brief extends the matrix adding padding in rows and cols as defined in the parameters
  \details
  @param NrowPaddingAbove - number of rows to be added above the first row
  @param NrowPaddingBelow - number of rows to be added below the last row
  @param NcolPaddingBefore - number of cols to be added before the first col
  @param NcolPaddingAfter - number of rows to be added after the last col
  @return returns the reference to the input matrix

  All padding parameters must be greater than zero. If they are not - they are set to zero.

  \date 2009/11/20 12:56:36
  \author Bartosz Lew
*/
matrix<double>& cpeds_pad_matrix(matrix<double>& m, long NrowPaddingAbove, long NrowPaddingBelow, long NcolPaddingBefore, long NcolPaddingAfter);
/*!
  \brief converts matrix to C linear array
  \details
  @param size - gives the size of the resulting array
  @param rowMajor - defines the ordering of the elements in the array; if false the ordering will be column major; if true it will be row major
  @return newly allocated array of size m.RowNum()*m.ColNum() in case of no padded cols/rows

  The convention for m matrix object is that the first coordinate refer to row number and second to the column number

  \date 2009/11/12 17:18:49
  \author Bartosz Lew
*/
  /* @param rowPadding - number of cells at the end of each column to be added to the output matrix (the padded values will be zero) */
  /* @param colPadding - number of cells at the end of each row to be added to the output matrix (the padded values will be zero) */
double* cpeds_matrix2array(const matrix<double>& m, long& size, bool rowMajor=true);

//! convert double array to float array; new array is allocated, old array is not deleted.
float* float_double_array(double* t, long size);


//! shifts array by N cells according to the fwd direction flag (default: forward) and returns 
void cpeds_shift_array(double*t, long size, long N, bool fwd=true);
void cpeds_shift_array(fftwComplex*t, long size, long N, bool fwd=true);

/*!
  \brief converts array to matrix
  \details
  The parameters are generally as for matrix2array function
  @param vecSize - is the size of the first vector stored in the array - for the rowMajor ordering this will be the number of columns in the matrix and for the colMajor ordering it will be the rows number in the matrix
  @return the rectangular matrix corresponding to the input array

  \date 2009/11/23 22:41:58
  \author Bartosz Lew
*/
const matrix<double> cpeds_array2matrix(const double* t, long size, long vecSize, bool rowMajor=true);

//! The name is self explanatory.
void cpeds_change_matrix_ordering_from_rows_to_cols_major(double *M, long rows, long cols);

/*! 
	\brief Calculates the covariance matrix of the measured data in Dvec table
	\details 
	@param Dvec - pointer to linear array of size vec_size*vec_num
	@param vec_num - number of vectors in the Dvec array
	@param vec_size - size of the vector
	@param diagonal - If diagonal is true than offdiagonal terms are not calculated and set zero (this is faster). In this case the return array is of size vec_size
	otherwise it is of size vec_size x vec_size and the corresponding covatiance matrix is symmetrical
	@return pointer to the covariance matrix
	
	Calculates the covariance matrix of the measured data in Dvec table, organized in vec_num vectors 1-row vectors, each of size vec_size and form (1___x___vec_size) 
	If i=0..vec_num-1 iterates vector index and j=0..vec_size-1 iterates variate in i'th vector then the ordering of the Dvec array is i-major: i.e. vector index major.
	Hence the Dvec is a set of vectors where each row-vector is a single measurement of all variates
	The resulting cov matrix is a square var_num___x___var_num size symmetric matrix 
	given by pointer cov.
		
	The array pointed by the returned pointer need not be allocated. It is allocated in this function.
		
	TODO: add option of calculate un-biased estimator of the covariance matrix

	\date Jan 29, 2011, 1:31:33 PM
	\author Bartosz Lew
*/
double * cpeds_calculate_covariance_matrix(double *Dvec, long vec_size, long var_num, bool diagonal=false);

//! Returns the quantile probability of occurance of x value in the t array of data of size ts
/*! Calculates the quantile probability of getting value x from distr given in t of size ts */
/*! The probability is of getting a measurment deviating in abs value from second quartile (Q24) by more than the measured value. */
/*! It's not cumulative so it calculates relaive to the 2nd  quartile, so very small and very big deviations will be unprobable. */
/*! Unlike in the Q probability case (defined as in the GSL library). */
/*! This routine assumes that the t array is sorted in rising order */
double cpeds_quantile_P(double x,double * t,long ts);

//! Same as above but this doesn't require data to be sorted. The sorting is done first.
double cpeds_quantile_P_sort(double x,double * t,long ts);

//! As above, but the sign of the returned probability indicates the measurment position with regard to 2nd quartile (Q24): - value below the Q24, + value above Q24
/*!  CAUTION !! assumes that the data in t are sorted in rising order */
double cpeds_quantile_P_signed(double x,double * t,long ts);

//! As above but does the sorting of the t array first.
double cpeds_quantile_P_sort_signed(double x,double * t,long ts);

//! As above but the the data are not sorted (you need to do this first).
/*! The improvement is that linear interpolation is done between the sorted values of the probed PDF  */
double cpeds_quantile_P_signed_interpolated(double x,double * t,long ts);

//! As above but the sorting is performed before the probability calculation.
/*! P=2 is a control value to detect in case something unpredicted happened. */
/*! Routine can also extrapolate probability using gaussian exponential lower/upper tail extrapolation */
/*! CAUTION !! assumes that the data in t are sorted in rising order */
double cpeds_quantile_P_sort_signed_interpolated(double x,double * t,long ts);

/*!
	\brief calculate a quantile value for a given array t of size ts
	\details 
	@param t - pointer to an array
	@param ts - array size
	@param quantile - requested quantile value [0..1]
	@param deleteAfter - if true, the array will be deleted when done
	@return quantile value

	\date Sep 23, 2014, 8:28:39 AM
	\author Bartosz Lew
*/
double cpeds_quantile(double* t, long ts, double quantile, bool deleteAfter=true);

//! The signed routine for calculation gaussian extrapolation (interpolation is also possible) given set of data.
/*! the tail probabilities are caclulated separatelly depening on which tail the data sits. and  */
/*! to that tail the gauss curve is fitted. the sewing assumes that the CDF should reach 2/N value at the x[0] */
/*! (and also at x[ts-1]). the sign of the probabilitiy is only indication on which tail the data sits. */
/*! SO THIS ROUTUNE RETURNS THE TAIL PROBABILITY NOT UPPER OR LOWER TAIL PROB. */
/*! IF YOU WANT A TAIL PROBABILITY -- IT'S JUST 0.5 OF THE RETURNED VALUE SINCE THE GAUSSIAN PDF IS SYMMETRICAL */
/*! WARNING ! if this routine is used as fitting gauss curve within the MCPDF probed region, some |P|>1 */
/*! are possible around Q24. to be fixed ! */
/*! CAUTION !! assumes that the data in t are sorted in rising order */
/*! the input data MUST BE SORTED IN ASCENDING ORDER */
double cpeds_extrapolate_gauss(double x, double* t, long ts);

//! Same as above but the sorting is first done
double cpeds_sort_extrapolate_gauss(double x, double* t, long ts);

double cpeds_extrapolate_linear(double x, double x1, double y1, double x2, double y2);
//! Calculates the lower tail inverse gaussian CDF
/*! Invoke with /f$ P \in \{0,p/2\} /f$ for finding the numer of sigmas deviation */
/*! E.g. invCDFGauss ( 0.0027/2) =~ -3 (sigma) (and remember about the signs) */
double invCDFGauss(double p);

//! do the power on integers
/*! returns /f$ a^b /f$; a and b must be natural numbers and a > 0  */
long cpeds_pow_int(long a, long b);

//! Calculate the Wigner3J coefficient with m=0
double cpeds_W3jm0(long l1, long l2, long l3);
/* double cpeds_W3jm0_lidia(long l1, long l2, long l3); */

//! Returns true when the Wj3 (all m1=m2=m3=0) is zero for given l1,l2,l3
bool cpeds_W3jm0_non_zero(long l1, long l2, long l3);

bool cpeds_W3jm0_lcond(long l1, long l2, long l3, long i);
/* double cpeds_CDFGauss(double u); */


//! calculates the spherical bessel functions of x -- adopted from cmbfast code
double cpeds_spherical_bessel_fn(long l,double x);

//! Calculates integral of the function given in array fx sampled on grid points given in array x of size size
/*! The x table and fx table must correspond to each other and the arguments of the functions must be given in increasing order */
double cpeds_integrate_1d(long size, double *x, double *fx);

//! As above but the integration is done from index imin to imax.
double cpeds_integrate_1d(long size, double *x, double *fx, long imin, long imax);


//! Performs GSL cubic spline interpolation with or without periodic boundary conditions
/*! using the input data passed in arrays X and Y of size N over the arguments values */
/*! passed in the Xint array of size Nint. The result is returned in the Yint array of size Nint. */

/*! IMPORTANT: the supplied points must be sorted according to the increasing X value */

/*! IMPORTANT: if the periodic boundary conditions were requested then the last value in the Y array */
/*! must be the same one as the first value. This is needed to know the length of the period. */
/*! Then if the interpolated points are requested formally from outside of the range defined in the */
/*! X array then they will be transformed accordingly to make sense out of the results. */

/*! WARNING: The thing I just said only works within +- one period from the given range defined in the */
/*! X array. It doesn't make sense to make this thing work ab infinity because if the signal is periodic */
/*! then one can just copy the results as many times as needed. */
double * cpeds_interpolate(double *X,double *Y, long N, double *Xint, long Nint, string type, bool check_points);
double * cpeds_interpolate(double *X,double *Y, long N, double *Xint, long Nint, string type);

double cpeds_bilinear_interpolation(double x1, double x2, double y1, double y2, double f11, double f12, double f21, double f22, double x, double y);
/*!
	\brief routine to calculate c matrix coefficients for bicubic interpolation it is called by bicubic_interpolation function
	\details 
	@param f - function values at four points around the square within which the interpolated value is calculated
	@param fx - df/dx at four points around the square within which the interpolated value is calculated
	@param fy - df/dy at four points around the square within which the interpolated value is calculated
	@param fxy - d2f/dxdy at four points around the square within which the interpolated value is calculated
	@param dx - length of the square along x direction
	@param dy - length of the square along y direction
	@param c - 2d array in row-major format of size 4x4. This is calculated.
		
	All arrays must be initialized. The size of y,y1,y2,y12 arrays is 4.

	\date Jan 12, 2012, 1:43:05 PM
	\author Bartosz Lew
*/
void cpeds_bicubic_interpolation_ccoef(double* y, double* y1, double* y2, double* y12, double d1, double d2, double c [][4]);
/*!
	\brief routine to calculate c matrix coefficients for bicubic interpolation it is called by bicubic_interpolation function
	\details 
	@param f - function values at four points around the square within which the interpolated value is calculated
	@param fx - df/dx at four points around the square within which the interpolated value is calculated
	@param fy - df/dy at four points around the square within which the interpolated value is calculated
	@param fxy - d2f/dxdy at four points around the square within which the interpolated value is calculated
	@param dx - length of the square along x direction
	@param dy - length of the square along y direction
	@param c - 2d array in row-major format of size 4x4. This is calculated.
	@param x1l - lower bound of the grid cell in direction 1
	@param x1u - upper bound of the grid cell in direction 1
	@param x2l - lower bound of the grid cell in direction 2
	@param x2u - upper bound of the grid cell in direction 2
	@param x1 - x coordinate of the interpolation point
	@param x2 - y coordinate of the interpolation point
	@param ansy - reference to the variable that will hold the interpolated value
	@param ansy1 - reference to the variable that will hold the derivative along direction 1 at the interpolated point
	@param ansy2 - reference to the variable that will hold the derivative along direction 2 at the interpolated point		
	
	All arrays must be initialized. The size of y,y1,y2,y12 arrays is 4.

	The code is based on the NUMERICAL RECIPES for c++
	
	This function is used in mscsfn classes for calculating interpolated function values between grid points.

	\date Jan 12, 2012, 1:43:05 PM
	\author Bartosz Lew
*/
void cpeds_bicubic_interpolation(double* y, double* y1, double* y2, double* y12, 
		const double x1l, const double x1u, const double x2l, const double x2u,
		const double x1, const double x2, double &ansy, double &ansy1, double &ansy2);


/*!
  \brief Calculates derivative of the 2D function stored in the matrix m.
  \details

  @param  yvals - The abscissa values for each row of the matrix are stored in yvals linear array
  @param  m - matrix of values of the function to be differentiated

  The derivative is calculated  in vertical direction - i.e. along columns of the matrix;
  The data stored in each column of the matrix represents values of the function to to be differentiated.

  The yvals don't have to be equally spaced, but they should be in increasing order.
  The derivative is done in place so the original function values in the matrix are replaced with its derivative
*/
void cpeds_matrix_derivative(matrix<double>* m, double* yvals);

/*!
	\brief calculate derivative 
	\details 
	@param
	@return
	
	In order to minimize numerical noise the function should 
	be sampled on equal grid. In case of periodic calculations, this means that
	the distance between the first and the last point should not be  equal to the period,
	bacause they are formally the same point.

	\date Nov 12, 2017, 5:48:24 PM
*/
void cpeds_derivative(double* x, double* y, long size, double* periodX=NULL, double *periodY=NULL);
double cpeds_getMinAbs3(double v1, double v2, double v3);
//cpeds_PDF_info* cpeds_get_PDF_info(double *X,double *Y, long N, cpeds_queue_double *CL);

/**********************************************************************************/
/*******************   IO STUFF  *************************************************/
/**********************************************************************************/


//! this functions measures the number of cols in the first line of the file fn
/*! it was tested and it works for the normal linux txt files */
/*! the field delimiter is single space */
/*! it works good whether or not the last line in file ends with \n (this doesn't matter here since it's only for the first line) */
/*! if you run it on some strange stuff then you better check the results. */
long cpeds_get_cols_num_first_ln(strarg fn,char * lastc);

/*!
  \brief returns the number of columns in the first line of the file
  \details
  @param fn - file name to examine
  @param commentedFile - indicates whetehr the file is commented or not\n
  if true: then the routing will skip to the first non-commented line.\n
  The default comment symbol is #. This is not implemented yet.
  @param maxLineSize - maximal size of line in bytes
  @return number of columns. -1 if file doesn't  exist
  This is a smarter implementation of the routine above: long cpeds_get_cols_num_first_ln(strarg fn,char * lastc)
  the file doesn't have to have the same field separators, it can have both the
  numerical and non-numerical entries and there can be many field separators between two neighboring columns.


  \date 2010/04/30 13:17:15
  \author Bartosz Lew
*/
long cpeds_get_cols_num_first_ln(strarg fn, bool commentedFile=false, long maxLineSize=1000 );

//! This function measures the number of columns starting from a position in the file indicated by the f pointer
/*! It works in the same way as the cpeds_get_cols_num_first_ln but does not close the file. */
/*! The position of the f pointer is increased by the line width since this routine reads the whole line untill the */
/*! \n character */
long cpeds_get_cols_num(FILE *f,char * lastc);


/*!
	\brief return the size of the file in bytes
	\details 
	@param fName 
	@return

	\date Dec 13, 2013, 11:00:55 AM
	\author Bartosz Lew
*/
long long cpeds_get_file_size(string fName);

long long cpeds_get_txt_file_lines_count(string fName);

////! This routine checks the txt file given by fn, and returns the cpeds_queue object that contains the
////! number of columns (space separated words) in each row of the file and much other useful info.
//cpeds_queue<long>* cpeds_get_txt_file_cols_rows(strarg fn);

//! checks if the file exists
bool cpeds_fileExists(string fname);


/*!
  \brief load matrix from file
  \details
  @param fileName - file name
  @param how - describes various options for loading\n
  In particular:\n
  if "header" keyword is included in the how string then it is assumed that the first row in the file contains two integer numbers  - the number of rows and columns of the matrix\n
  if "binary" keyword is included then the file is read as a double valued binary file\n
  if no keyword is given then it is assumed that the file is a text file and the file format will be detected automatically in case of the text files\n
  In case of the binary files the format will be assumed to be single column vector
  if "float" is also given to "binary" then 4byte floats are read from the file

  @param result - pointer to the operation status variable. If not needed, send null of skip parameter. If given 0 result indicate success

  @return returns a matrix containing the loaded data

  The format for saving is the raw-major format.

  \date 2009/11/06 16:50:00
  \author Bartosz Lew
*/
const matrix<double> cpeds_matrix_load(string fileName, string how="", long * result=NULL);
/*!
  \brief  saves matrix data to file
  \details

  See above for more details.
  
  The matrix M is saved to a txt file in a way that the first dimension (i) increases with row number
  and second dimension (j) increases with column number

  \date 2009/11/06 16:54:06
  \author Bartosz Lew
*/
long cpeds_matrix_save(const matrix<double>& M, string fileName, string how="",int precision=20);

double cpeds_get_memory_usage();

std::tuple<string, string, string> cpeds_get_dirname_filebase_ext(string s);

#endif




#ifndef __CHEALPIX_H__
#define __CHEALPIX_H__

/* -------------------- */
/* Constant Definitions */
/* -------------------- */

#ifndef HEALPIX_NULLVAL
#define HEALPIX_NULLVAL (-1.6375e30)
#endif /* HEALPIX_NULLVAL */

/* --------------------- */
/* Function Declarations */
/* --------------------- */

/* pixel operations */
/* ---------------- */
void ang2pix_nest(const long nside, double theta, double phi, long *ipix);
void ang2pix_ring(const long nside, double theta, double phi, long *ipix);

void pix2ang_nest(long nside, long ipix, double *theta, double *phi);
void pix2ang_ring(long nside, long ipix, double *theta, double *phi);

void nest2ring(long nside, long ipnest, long *ipring);
void ring2nest(long nside, long ipring, long *ipnest);

void mk_pix2xy(int *pix2x, int *pix2y);
void mk_xy2pix(int *x2pix, int *y2pix);


#endif /* __CHEALPIX_H__ */
