#include "Mscs-function.h"

#ifndef MSCS_MAP_WINDOW_FUNCTION
#define MSCS_MAP_WINDOW_FUNCTION

/*!
  \class mscsWindowFunction
  \brief Implements the window function functionality
  \details 
  This may be either a pixel window function of beam window (transfer) function.

  \date 2009/06/05 17:30:21 
  \author Bartosz Lew
*/
class mscsWindowFunction : public mscsFunction {

 public:
  mscsWindowFunction();
  mscsWindowFunction(const mscsFunction& rhs) : mscsFunction(rhs) { }

  //! another constructor but probably not much useful at the moment
  mscsWindowFunction(string _wf_name, double _kmax);
  mscsWindowFunction(string _wf_name, cpeds_VerbosityLevel verbosity=CPEDS_defaultVerbosityLevel);
  /*!
    \brief Initiates with zeros the function of size size.
    \details If _lmax is given -1 then the function will not be initiated
    
    \date 2009/10/20 21:34:22 
    \author Bartosz Lew
  */
  mscsWindowFunction(string _wf_name, long _lmax=-1, cpeds_VerbosityLevel verbosity=CPEDS_defaultVerbosityLevel);
  mscsWindowFunction(const mscsWindowFunction& rhs) : mscsFunction(rhs) { }
  virtual ~mscsWindowFunction();


  /*!
    \brief creates a gaussian kernel in multipole space
    \details 
    @param FWHM - the full width at half maximum of the beam in [rad]
    
    \date 2009/06/05 17:29:56 
    \author Bartosz Lew
  */
  virtual const mscsWindowFunction& make_gaussian_kernel(double FWHM);

  /*!
    \brief creates a gaussian kernel in k-space
    \details 
    @param FWHM - the full width at half maximum of the beam in [rad]
    
    \date 2009/06/05 17:29:56 
    \author Bartosz Lew
  */
  virtual const mscsWindowFunction& make_gaussian_kernel_kspace(double FWHM);

  /*!
	\brief populates function with the gaussian filter using the current arguments of the function as k modes
	\details 
	@param R - scale parameter in units [1/k] - i.e. inverse argument unit so that product is unitless
	@return
	
	The gaussian - window function will have form:
		f[i] = exp ( - k^2 R^2 / 2 * 4.0*PI^2 )
		where k = k(i) and R is the provided scale

	Note that this is fourier transform definition dependent. For the fourier transform definition that is used in
	fftw we need the extra 4PI^2 factor (applied in this implementation) to the usual definition of this trnasfer function.

	\date Oct 8, 2010, 11:36:25 AM
	\author Bartosz Lew
  */
  mscsFunction& mkGaussianKernel(double R, double k0=0);


  /*!
	\brief populates function with the top-hat filter using the current arguments of the function as k modes
	\details 
	@param R - scale parameter in units [1/k] - i.e. inverse argument unit so that product is unitless
	@param fftConvention - fft definition convention. Possible values are:\n
		0 - fftw convenction
		1 - physical fft convenction
	
	@return
	
	The top hat - this window function ussually has the form (convention 1):
		f[i] = 3 * (sin(kR)/(kR)^3 - cos(kR)/(kR)^2 )
		where k = k(i) and R is the provided scale

	But,
	Note that this is fourier transform definition dependent. For the fourier transform definition that is used in
	fftw we need a different definition (convention 0):
	
		f[i] = sin(2pi kR)/(2pi kR)
	

	The type of fft convenction is defined by fftConvention parameter

	\date Oct 8, 2010, 11:36:25 AM
	\author Bartosz Lew
  */
  mscsFunction& mkTopHatKernel(double R, int fftConvention=0);
  
//  virtual const mscsWindowFunction& make_gaussian_kernel_kspace(double kmin, double kmax, double dk, double R);

  /*!
    \brief creates a unitary kernel
  */
  const mscsWindowFunction& make_unit_kernel(long lmax=-1);


  /*!
	\brief generate Hann window function for the data of size N
	\details 
	@param N - number of points in the window
	@return returns this

	\date Feb 28, 2011, 5:14:40 PM
	\author Bartosz Lew
  */
  mscsFunction& mkHannWindow(long N);


  /*!
	\brief make exponential window in frequency space for time series convolutions in Fourier space
	\details 
	@param k0 - frequency below which the window function is 1.
	@param tC - time constant that defines the exponential slope
	@return

	  The window is defined as follows:\n
	  w_i = 1 for k < k0
	  w_i = exp(-(k-k0)/tC) for k >=k0
	  
	\date May 17, 2011, 3:14:44 PM
	\author Bartosz Lew
  */
  mscsFunction& mkExpWindow(double kmin, double kmax, long Nk, double k0, double tC);
  /*!
	\brief make exponential high pass window in frequency space for time series convolutions in Fourier space
	\details 
	@param k0 - frequency below which the window function is 1.
	@return

	  The window is defined as follows:\n
	  w_i = 1 for k >= k0
	  w_i = exp(-(k0-k)/tC) for k >=k0

	\date Mar 10, 2013, 3:35:49 AM
	\author Bartosz Lew
   */
  mscsFunction& mkHighPassExpWindow(double kmin, double kmax, long Nk, double k0, double tC);

  /*!
	\brief a low-pass filter transfer function
	\details 
	@param fmin - omega_min
	@param fmax - omega_max
	@param Nf - number of points in min max range
	@param RC - time constant
	@return returns the transfer function in frequency space
	
	The transfer function is defined as:
	tf(w) = 1/Sqrt(1+ (w RC)^2 )

	The time constant RC = 1/(2pi f_c), where f_c is the cut-off frequency where the amplitide of the input signal drops by a factor of 2.

	\date May 24, 2011, 12:14:44 PM
	\author Bartosz Lew
  */
  mscsFunction& mkLowPassWindow(double fmin, double fmax, long Nf, double RC);

  /*!
	\brief a low-pass filter transfer function
	\details 
	@param fmin - omega_min
	@param fmax - omega_max
	@param Nf - number of points in min max range
	@param RC - time constant
	@return returns the transfer function in frequency space
	
	The transfer function is defined as:
	tf(w) = 1 for f <= f_c 
	tf(w) = 0 for f >  f_c 

	where f_c is the cut-off frequency defined as f_c = 1/(2pi RC) 

	\date Mar 28, 2014, 9:35:25 AM
	\author Bartosz Lew
  */
  mscsFunction& mkLowPassStepWindow(double fmin, double fmax, long Nf, double RC);

  
  /*!
	\brief make smooting kernel function used in gadget-2 for SPH
	\details 
	@param Rmin - minimal radial value
	@param Rmax - maximal radial value
	@param hsml - smoothing length; if 0 then not used and the kernel is unitless
	@param N - number of points probing this kernel
	@return returns this

	The Rmin and Rmax are the normalized radial values of the kernel expressed in units of hsml
	

	\date Jan 27, 2012, 3:54:35 PM
	\author Bartosz Lew
   */
  mscsFunction& mkSPHkernelGadget2(long N=100, double Rmin=0, double Rmax=2, double hsml=0);

  
  /*!
	\brief generate air mass- elevation relation for realistic profile, spherical atmosphere model
	\details 
	@param from - minimal elevation [deg]
	@param to - maximal elevation [deg]
	@param dh - step in elevation [deg]
	@return returns this
	
	The small elevations ~<3 deg are extrapolated using linear relation down to h=0 

	\date Sep 13, 2012, 2:26:42 PM
	\author Bartosz Lew
   */
  mscsFunction& mkAirMassElevationRelation(double from=0, double to=90, double dh=0.1);

  
  
  static double kernelGadget(double R);

  
  
  
  
  const mscsWindowFunction& operator=(const mscsWindowFunction& rhs) { if (this != &rhs) { this->mscsFunction::operator=(rhs); } return *this; }
  const mscsWindowFunction& operator=(const mscsFunction& rhs) { mscsFunction::operator=(rhs);   return *this; }

 protected:


};

#endif
