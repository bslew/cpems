/*!
  \mainpage Mscs - Map Statistics package

  \section Introduction
  This package was initially started probably around mid 2005 when I started
  working on my PhD in Japan and was developed to ease numer of numerical
  problems that needed to be solved.

  The package basic usage is to generate and analyse CMB high resolution
  spherical maps from experiments like WMAP.
  The package allows for map generation,
  loading, saving, many useful operations on it, as well as numer of
  more ellaborated statistics (like, tests of gaussianity, power spectrum
  reconstructions,  various statistics etc.

  The functionality of the package is continously extended (eg. work on
  non-Gaussian simulations and bispectrum statistics are under way.
  Polarization maps management and
  operations are also in plans as well as needlets decomplsition).

  The core of the package is the Mscs library on which several programs operate.
  These programs will be documented elsewhere.

  The general structure of the Mscs package can be viewed <a href="../../Mscs-structure-by-UML/Mscs-structure-by-UML.html">here</a>

  \section advanced Compilation and Requirements
  To compile this code several other libraries are required, see
  INSTALLATION.txt file in the distribution directory for details on that.

  \section Conventions
  There are some conventions that are generally followed in the code
  about few things. Knowing them might sometimes make life a bit easier.

  \li \b angles - generally kept in radians unless otherwise noted
  \li <b> coordinate systems</b> - The generic convention is that Mscs
  operates in the usual galactic coordinate system (l,b) with
  \f$ b\in <\pi/2,\pi/2> \f$ and \f$ l\in <0,2\pi) \f$.
  \li

  \section Glossary
  \li DA - differential assembly
  \li nside - healpix pixelization system number of pixels on the side of the
  fundamental pixelization region (see some healpix papers for more info)

  \section Licence
  This code is not yet enough smoothed out to exist as a release version
  (no time :-( ) and so as for the moment exists only as a pretty much
  private code. But in future an official release is planned of course.

  \date 2009/05/27 10:42:27
  \author Bartosz Lew
*/

//#include <complex>
#include <string>
#include "matrix.h"

//#include "Mscs-map.h"

#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
#define STD std
#else
#define STD
#endif

#ifndef MSCS_COMMON_DEFS
#define MSCS_COMMON_DEFS
#define _Mscs_version_ "Mscs-1.0"

#define MSCS_NUMERICAL_ACCURACY 1e-8


/* #ifndef __MscsTypeDefs__ */
/* #define __MscsTypeDefs__ */
//typedef char filenamestr[1000];
//typedef const char strarg[500];
//typedef filenamestr package_dirs[30];
//typedef strarg strarg_dirs[30];
typedef struct { long r,b,g; } MscsColor; //!< color structure defined on tripet of longs; The intensities should be within range of [0,255]
typedef struct { double r,b,g; } MscsColorD; //!< color structure defined on tripet of doubles; The intensities should be within range of [0,1]

typedef struct {
	int procNum; //!< number of threads that is running
	int st_thread;//!< starting number of the thread
	int en_thread; //!< ending number of the thread
	char sim_run_file_pref[500]; //!< prefix of the run file that holds the current status information (includes absolute path)
	char sim_run_file_suff[100]; //!< suffix of the run file that holds the current status information
} Mscs_run_control_t;


/* #endif */




#endif

