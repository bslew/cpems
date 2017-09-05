#include "Mscs-function.h"

#ifndef MSCS_MAP_BEAM_FUNCTION
#define MSCS_MAP_BEAM_FUNCTION

/*!
  \class mscsBeamFunction
  \brief Implements the circular beam profile function functionality
  \details 

  The beam is assumed the be a function b(theta) where theta is in radians
  and b is calibrated to unity at the main lobe.

  \date 2009/10/22 22:22:04 
  \author Bartosz Lew
*/
class mscsBeamFunction : public mscsFunction {

 public:
  mscsBeamFunction();
  mscsBeamFunction(string _wf_name);

  /*!
    \brief Initiates with zeros the function of size size.
    \details If _lmax is given -1 then the function will not be initiated
    @param _wf_name - name of the beam function
    @param FWHM - full width at half maximum [rad]
    @param nHWHM - the number of HWHMs to cover in the half-beam digitization from one side of the beam
    @param _points_per_HWHM - self explanatory
    The actual angle values of the function points are located in the middle of the bin.
    
    \date 2009/10/20 21:34:22 
    \author Bartosz Lew
  */
  mscsBeamFunction(string _wf_name, double fwhm, long nFWHM, long _points_per_FWHM);
  mscsBeamFunction(const mscsBeamFunction& rhs) : mscsFunction(rhs) {  } 
  ~mscsBeamFunction();


  /************/
  /* HANDLERS */
  /************/
  double FWHM() const { return _FWHM; }
  


  /*!
    \brief creates a gaussian kernel
    \details 
    @param FWHM - the full width at half maximum of the beam in [rad]
    
    \date 2009/06/05 17:29:56 
    \author Bartosz Lew
  */
  void make_gaussian_beam(double fwhm);
  void make_gaussian_beam();

  const mscsBeamFunction& operator=(const mscsBeamFunction& rhs) { _f=rhs._f; range=rhs.range;  return *this; }
  const mscsBeamFunction& operator=(const mscsFunction& rhs) { mscsFunction::operator=(rhs);  return *this; }

  /*************/
  /* OPERATORS */
  /*************/

  
 protected:

  void setFWHM(double fwhm) { _FWHM=fwhm; }

 private:
  double _FWHM;


};

#endif
