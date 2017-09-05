/*!
  \file implements the smoothing of the given directions field in the sky
*/

#ifndef CPEDSSMOOTH
#define CPEDSSMOOTH

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
#include "cpeds-math.h"
#include "cpeds-point_set.h"
#include "cpeds-direction_set.h"
#include "cpeds-project.h"

#include <fftw3.h>

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */


/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class cpedsSmooth
  \brief Encapsulates smoothing of a set of directions in the sky
  \details 
  This class makes gaussian smoothing of requested strength of the set of directions in the sky by

  first: projecting them onto a tangent plane via orthogonal projection\n
  second: FFTing it \n
  third: combining convolving it with the requested gaussian kernel and ffting back\n
  

  The original units that appear in the documentation mean refer to the initial units in which the points coordinates were
  provided (for the x,y,z coordinates case). Unless otherwise specified these are the basic units to use in this object.


  Convention clarifications:\n
  Whenever there's a field concerned (stored in a matrix object) it has rows and columns and the first dimension is rows (corresponds to X although 
  usually we write them vertically one under another),
  and second dimension (and the last one) are columns (which correspond to Y - usually we write them horizontally one after another). 
  This is the situation within the object and FFTW3.

  The situation from the user point of view is different where we have fields or coordinate systems and got used to the situation where X is something
  horizontal and Y is vertical axis. So whenever the smoothing is concerned and fwhm of the smoothing the fwhmX and fwhmY correspond to the user space,
  and the same applies to resolution parameters and anything that can be accesses by the user via handlers.


  IMPORTANT NOTE: If you use this class in order to smooth the directions in the sky via first projecting them on a plane and then gaussian smoothing 
  and back-projecting, then keep this in mind. The field generated on plane from a set of projected directions is associated with the true on-plane coordinates
  using a single point (from the projection direction) and the grid resolution information and a linear transformation. Since the stereographic projection
  on plane is a non-linear transformation, the back projection of cells of the field will not generally coorespond to the correct celestial coordinates.
  The further the cell from the projection direction is located, the larger distortion of the field will occur. To mitigate this effect,  one can restrict
  to smoothing only small portions of the sky (where the flat sky approximation works fine) or by developing a better method for associating the 
  field matrix with the true on-plane coordinates for every cell of the matrix. This is not implemented yet. For more info  see docs for associateMatrix2Field()
  method.
  
  

  \date 2009/11/05 23:16:51 
  \author Bartosz Lew
*/
class cpedsSmooth : public cpedsPointSet3D {


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
  cpedsSmooth();
  cpedsSmooth(cpedsSmooth& p);
  cpedsSmooth(const cpedsPointSet2D& p,double alpha=1.0);
  /*!
    \brief initiate and perform smoothing of the input point set
    \details 
    @param p - point set to be smoothed
    @param borderSize - size of the border of the field to be used (in original units of the point set)
    @param resolution - resolution with which to probe point set in order to make a rectangular field (in original units of the point set)
    @param fwhmX - size of the smoothing in X direction (in original units of the point set)
    @param fwhmY - size of the smoothing in Y direction (in original units of the point set)
    @param alpha - currently rather not used ( so don't bother for the moment)
    @return
  */
  cpedsSmooth(const cpedsPointSet3D& p, const cpedsPoint3D& borderSize=cpedsPoint3D(), const cpedsPoint3D& resolution=cpedsPoint3D(1,1,0), double fwhmX=0, double fwhmY=0, double alpha=1.0);
  cpedsSmooth(const cpedsDirectionSet& d, const cpedsDirection& n, const cpedsDirection& borderSize=cpedsDirection(), const cpedsDirection& resolution=cpedsDirection(1,1,0), double fwhmX=0, double fwhmY=0, double alpha=1.0);


/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */

  void initiate_smoothing_parameters();

  /*!
    \brief defines various smoothing options; this method defaults to the smoothing of the field that is given as a set of directions in the sky
    \details 
    @param resolution - defines the angular resolution of the smoothed field in radians for longitude and latitude.\n
    Eg. in order to have angular resolution of 0.001 deg x 0.001 deg in each direction use cpedsDirection(0.001,0.001)*PI180 for this argument

    @param smoothingLength - defines the angular smoothing lenght in radians of the field in each direction.\n
    Eg. in order to have smoothing of size of 1 deg per 1 deg  use cpedsDirection(1,1)*PI180 for this argument

    @param fieldBorder - defines the angular size in radians of the field border in each direction.\n
    The resulting field will be larger than the one resulting from the input directions set by this given amount in every direction along longitude and latitude

    There are two smoothing types to be supported:\n
    "fft" - smoothing performed by convolving with the gaussian kernel in fourier space - the set of points is first cast onto a dense grid with resolution defined by the resolution parameters\n
    "convolve" - real space based convolotion with gaussian beam profile. Currently NOT IMPLEMENTED

    \date 2009/11/12 11:55:42 
    \author Bartosz Lew
  */
  void set_smoothing_parameters(const cpedsDirection& resolution, const  cpedsDirection& smoothingLength, const cpedsDirection& fieldBorder);




















  /*!
    \brief defines the global resolution parameger of the smoothing
    \details
    @param a - resolution parameter value; the larger value results in smaller number of pixels in the smoothed map; the smaller value (than unity) results in larger number of pixels in the smoothed map

    The smoothing is performed via convolution of the rectangular field with a smoothing kernel (which can be eg. a gaussian beam).
    The field is defined in a form of grid with nods with spacing defined as alpha*rms_X,Y,Z ( for X,Y and Z dierection respectively) where rms is the rms value
    of the point set in X,Y and Z direction respectively. Hence alpha=0.1 (the default value) will result in field that resolves the points 10 times better than
    rms value of the points in the set.

    \date 2009/11/10 23:02:42
    \author Bartosz Lew
  */
  void setAlpha(double a) { _alpha=a; }
  double alpha() const { return _alpha; }

  //! set or get the resolution to resolve the points of the field in original units; careful - this returns the reference to the resolution object
  cpedsPoint3D& res() { return _resXY; }
  //! set the resolution to resolve the points of the field in original units
  void setRes(const cpedsPoint3D r)  { _resXY=r; }
  cpedsPoint3D getRes() const { return _resXY; }

  //! set or get the border size of the field in original units; careful - this returns the reference to the resolution object
  cpedsPoint3D& border() { return _fieldBorder; }
  //! set the border size of the field in original units
  void setBorder(const cpedsPoint3D b) { _fieldBorder=b; }
  cpedsPoint3D getBorder() const { return _fieldBorder; }

  //! returns the smoothing length object for the field in original units
  const cpedsPoint3D smoothingLength() { return cpedsPoint3D(_fwhmX,_fwhmY,0); }

  /*!
    \brief defines the smoothing length in each direction
    \details 
    @param fwhmX - smoothing lendth in X direction [original units - in the plane - not on the sphere]
    @param fwhmY - smoothing lendth in Y direction [original units - in the plane - not on the sphere]

    If you want to initiate the smoothing length when smoothing a set of directions on the sphere then use set_smoothing_parameters method instead
    or apropriate constructor. Don't use this to is you want to smooth directions and don't know anything about the on-plane scale onto which these directions
    will be translated upon projection.
    
    \date 2009/12/09 12:19:19 
    \author Bartosz Lew
  */
  void setSmoothingLength(double fwhmX, double fwhmY) { _fwhmX=long(fwhmX/getRes().x());  _fwhmY=long(fwhmY/getRes().y()); }
  //! returns the smoothing parameters in the original units.
  cpedsPoint3D getSmoothingLength() const { return cpedsPoint3D(_fwhmX,_fwhmY,0)*getRes(); }
  //! returns the smoothing parameters as a number of pixels in the matrix, that prepresent fwhm in each direction 
  cpedsPoint3D getSmoothingLengthPix() const { return cpedsPoint3D(_fwhmX,_fwhmY,0); }

  //! clears all the structures of the object
  void clearAll() { _m.SetSize(0,0); _pssm.clear(); clear(); }

  /* double& rangeX() { return _rX; } */
  /* double& rangeY() { return _rY; } */

  /* /\*! */
  /*   \brief  */
  /*   \details  */
  /*   @param b - resolution parameter of the convolution */

  /*   This parameters defines how many points will be used to generate the convolution kernel around given grid nod of the field. */
  /* *\/ */
  /* void setBeta(double b) { _beta=b; } */
  /* double beta() const { return beta; } */


  /* /\*! */
  /*   \brief defines the smoothing angular scale for each direction of the field in radians */
  /*   \details  */
  /*   @param fwhmX - smoothing scale in radians for X direction */
  /*   @param fwhmY - smoothing scale in radians for Y direction */
    
  /*   This scale is translated onto the linear scale via the proper projection of the corresponding arc onto the plane tangent to the center point of the arc length */
  /*   situated in the center of the field. */
    
  /*   \date 2009/11/10 23:35:35  */
  /*   \author Bartosz Lew */
  /* *\/ */
  /* void setSmoothingAngScale(double fwhmLon, double fwhmLat) { _fwhmLon=fwhmLon; _fwhmLat=fwhmLat; deriveSmoothingScale(); } */

  /* //! sets the linear on-plane smoothing scale for each direction */
  /* void setSmoothingScale(double fwhmX, double fwhmY) { _fwhmX=fwhmX; _fwhmY=fwhmY;  } */
  

  /* //! derives the _fwhmX and _fwhmY parameters based on set _fwhmLon and _fwhmLat parameters */
  /* void deriveSmoothingScale(); */

  

  void printInfo() const;


  /*!
    \brief smooth field with gaussian beam
    \details 
    @param nfwhmX - smoothing length in number of points per fwhm
    @param nfwhmY - smoothing length in number of points per fwhm
    @param m - an alternative matrix that can be smoothed instead of the default one -- THIS IS NOT FULLY IMPLEMENTED YET.\n
    @return returns the smoothed field in the matrix form
    
    \date 2009/11/06 17:32:25 
    \author Bartosz Lew
  */
  const matrix<double> smoothGauss(long nfwhmX, long nfwhmY, matrix<double>* m=NULL);
  /*!
    \brief smooth field with gaussian beam
    \details
    @param fwhmX - number of pixels in the matrix that represent the fwhm in X direction
    @param fwhmY - number of pixels in the matrix that represent the fwhm in Y direction
    @param m - an alternative matrix that can be smoothed instead of the default one -- THIS IS NOT FULLY IMPLEMENTED YET.\n
    The way it can be used is to eg. pass the resulting smoothed matrix from one smoothing as an imput to another smoothing since each smoothing 
    is done from the beginning - ie. starting from the initial point set, so smoothing first in eg. X direction and then in Y direction gives you the same result
    as smoothing in Y directions only. If you want to smooth in both directions you can use smoothGauss or apply smoothing in X direction and then in Y direction
    supplying the resulting matrix from first smoothing to the second smoothing.
    @return returns the smoothed field in the matrix form

    The fwhmX and fwhmY are supposed to be integers. If you had set all the smoothing parameters before, you don't need to supply any parameters to this method.\n
    If you want to supply the smoothing values in original units rather than in number of pixels the use another smoothGauss method.
  */
  const matrix<double> smoothGauss(double fwhmX=0, double fwhmY=0, matrix<double>* m=NULL);


  const matrix<double> smoothGaussX(double fwhmX=0, matrix<double>* m=NULL);
  const matrix<double> smoothGaussY(double fwhmY=0, matrix<double>* m=NULL);
  /*!
    \brief smooth field with gaussian beam along X direction only
    \details 
    @param nfwhmX - smoothing length in number of points per fwhm
    @param m - pointer to a matrix containing an alternative field to smooth
    @return returns the smoothed field in the matrix form
    
    \date 2009/11/06 17:32:25 
    \author Bartosz Lew
  */
  const matrix<double> smoothGaussX(long nfwhmX=0,matrix<double>* m=NULL);
  /*!
    \brief smooth field with gaussian beam along Y direction only
    \details 
    @param nfwhmY - smoothing length in number of points per fwhm
    @param m - pointer to a matrix containing an alternative field to smooth
    @return returns the smoothed field in the matrix form
    
    \date 2009/11/06 17:32:25 
    \author Bartosz Lew
  */
  const matrix<double> smoothGaussY(long nfwhmY=0, matrix<double>* m=NULL);


  /*!
    \brief generates a rectangular matrix populated with data points
    \details 
    @param resx, resy - resolution parameters; if the smoothing parametres are not pre-set with the set_smoothing_parameters function
    then you should provide the resolution information via these parameters. If the default values are given then the function will use the 
    parameters set by the set_smoothing_parameters method.
    
    The data points do not have to be found in the matrix grid nods at their exact original location.
    The input data points will be snapped to the grid 

    \date 2009/11/09 10:41:47 
  \author Bartosz Lew
  */
  const matrix<double> make_field(double resx=0, double resy=0);
  const matrix<double> make_field(long sizeX, long sizeY);

  const matrix<double> getSmoothedField() const {  return _m; }
  const cpedsPointSet3D& getSmoothedPointSet() const { return _pssm; }
  cpedsDirectionSet getSmoothedDirectionSet() const { return cpedsProject(_pssm).projectOnSphere(getProjectionPlane()); }


  /* virtual const cpedsSmooth& operator=(const cpedsPointSet3D &rhs); */
  virtual cpedsSmooth& operator=(const cpedsSmooth& rhs);

  cpedsDirection& projectionPlane() { return _projectionPlane; }
  cpedsDirection getProjectionPlane() const { return _projectionPlane; }
/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PROTECTED MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 protected:


/* ---------------------------- */
/* PROTECTED METHODS */
/* ---------------------------- */
  /*!
    \brief associates (currently) first point from the input direction set with its location in the matrix accounting for the size of the border.
    \details 
    This is needed in case of the directions smoothing in order to generate the correct coordinates in the 
    output point set that will be back converted into a direction set.
    For regular point sets smoothing this is also needed in order to assign the correct coordinates to the interpolated (smoothed) new points
    in the matrix.
  */
  void associateMatrix2Field();
  

  /* //! sets the on-plane field ranges excluding the current definition of the field margins */
  /* void setFieldRanges(); */



  /*!
    \brief imposes the Hermitian symmetry on the 
    \details 

    \warning
    something is wrong with this routine -- there's some missing symmetry here or something else screwed up !!!
    FIX THIS BEFORE USING - OR USE make_field_symmetricalCol and make_field_symmetricalRow together instead
  */
  void make_field_symmetrical(fftw_complex* t, long vecNum, long vecSize);
  void make_field_symmetricalCol(fftw_complex* t, long vecNum, long vecSize);
  void make_field_symmetricalRow(fftw_complex* t, long vecNum, long vecSize);


  /*!
    \brief FFT along X
    \details 
    This transforms the field along the first direction (X) along which the rows change. They change slower than columns in the table t since the array is in row-major ordering - consistently with the fftw3 standards.
    Each column is FFTed independently.
  */
  void field_fftX(fftw_complex* t, long vecNum, long vecSize, bool fwd=true);
  /*!
    \brief FFT along Y
    \details 

    Parameters is in 

    This transforms the field along the last direction (Y) along which the column index changes faster than row index since the array is in row-major ordering
    Each row is FFTed independently.

  */
  void field_fftY(fftw_complex* t, long vecNum, long vecSize, bool fwd=true);

  /*!
    \brief FFT along X and Y
    \details 
    @param t - pointer to 2D ffw_conmplex array that is in X-major ordering, (Y-minor ordering) - X is vertical, Y - horizontal
    @param vecNum - number of rows (number of lines along X direction)
    @param vecSize - number of cells in each row (along Y direction)
    @param fwd - fft direction: true - forward, false - backward

    Performs a 2D FFT on the matrix data - stored as row-major ordered array with vecNum number of rows and vecSize number of columns
  */
  void field_fft(fftw_complex* t, long vecNum, long vecSize, bool fwd=true);

  /* void field_fft(fftw_complex* t, long vecNum, long vecSize); */
  /* void field_inv_fft(fftw_complex* t, long vecNum, long vecSize); */
  void field_convolve_Gaussian_kernel(fftw_complex* t, long vecNum, long vecSize, double fwhmX, double fwhmY);

/* ---------------------------- */
/* PROTECTED STRUCTURES */
/* ---------------------------- */


/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PRIVATE MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 private:


/* ---------------------------- */
/* PRIVATE METHODS */
/* ---------------------------- */


/* ---------------------------- */
/* PRIVATE STRUCTURES */
/* ---------------------------- */

  /* double _Xmin, _Xmax, _Ymin, _Ymax, _Zmin, _Zmax; //!< ranges of the projected field */

  /* THESE PARAMETERS ARE NOT USED IN FFT SMOOTHING MODE - BEGIN*/
  double _dx, _dy; //!< resolution parameters to resolve the smoothing kernel (not used in fft smoothing mode)
  double _rX, _rY; //!< convolution range parameter for each direction
  /* double _beta; //!< convolution integral resolution parameter  */
  /* THESE PARAMETERS ARE NOT USED IN FFT SMOOTHING MODE - END */

  //! resolution parameter [dimensionless] 
  //! \details The parameter defines the factor by which the output resolution of the fied will be multiplied in order to better resolve the input point set before the smoothing is done. Then, after smoothing appropriate recasting (downgrading) to the requested resolution is done. Eg. for alpha= 0.1 the input point set will be probed with 0.1*res() resolution in order to generte field more accurately. Then smoothing is done and then the field is downgraded to the requested resolution res.
  double _alpha; 
  double _fwhmLon, _fwhmLat; //!< angular smoothing scales for each direction
  double _fwhmX,_fwhmY; //!< number of points per gaussian fwhm in each direction of the smoothing \details This is in the matrix pixel space - not in the point set space


  double _Lx; //!< size of the box in x direction in the original units
  double _Ly; //!< size of the box in y direction in the original units
  

  cpedsPoint3D _resXY; // (on-plane) resolution parameters of the field to be smoothed
  cpedsPoint3D _smoothingLength; // (on-plane) resolution parameters of the field to be smoothed
  cpedsPoint3D _fieldBorder; // (on-plane) size of the border around the field to be smoothed


  //cpedsPointSet3D _ps; //!< the input point set to be smoothed
  matrix<double> _m; //!< contains the field to be smoothed
  matrix<long> _mask; //!< contains the id of the point(direction); -1 if doesn't have point(direction) counterpart
  cpedsPointSet3D _pssm; //!< smoothed point set - to be generated after the smoothing is done
  cpedsDirection _projectionPlane; //!< the direction of the tangent plane for projections
 protected:
  long _referencePointRow; //!< a row in the matrix _m in which the first point of the point set is stored
  long _referencePointCol; //!< a column in the matrix _m in which the first point of the point set is stored
  long _referencePointNumber; //! an index of the point in the list that corresponds to the _referencePointCol and _referencePointRow
};
#endif /* CPEDSSMOOTH */ 

