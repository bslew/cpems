/*!
  \file encapsulates a set of cpeds point objects in 2D
*/

#ifndef CPEDSPOINTSET2D
#define CPEDSPOINTSET2D

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
#include <QtCore/qlist.h>
#include "cpeds-math.h"
#include "cpeds-point2d.h"
#include <vector>
#include "mscsVector.h"
#include "subdomain.h"

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */


/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class cpedsPointSet2D
  \brief Encapsulates a set of cpeds point objects in 2D
  \details

  \date 2009/11/06 12:51:40
  \author Bartosz Lew
*/
//class cpedsPointSet2D : public QList<cpedsPoint2D> {
class cpedsPointSet2D : public mscsVector<cpedsPoint2D> {


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
  cpedsPointSet2D();
  cpedsPointSet2D(const QList<cpedsPoint2D> &q);
  cpedsPointSet2D(const cpedsPointSet2D &q);
  cpedsPointSet2D(long N, double *lon, double *lat);
  cpedsPointSet2D(const vector<cpedsPoint2D> &v);
  virtual ~cpedsPointSet2D();

/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */
  cpedsPointSet2D& append(const cpedsPoint2D& p) { push_back(p); return *this; }
  cpedsPointSet2D& append(const cpedsPointSet2D& ps) { for (unsigned long i = 0; i < ps.size(); i++) { append(ps[i]); }  return *this;  }
  long count() const { return size(); }
  cpedsPoint2D value(long i) { return cpedsPoint2D(at(i)); }

  double* getXvals(long* size) const;
  double* getYvals(long* size) const;

  void getRanges(double &minx, double &maxx, double &miny, double &maxy) const;
  double minX() const;
  double maxX() const;
  double minY() const;
  double maxY() const;

  virtual cpedsPointSet2D& operator=(const cpedsPointSet2D &rhs);
  virtual cpedsPointSet2D& operator=(const QList<cpedsPoint2D> &rhs);

  void print() const;
  virtual void save(string fname) const;

  virtual const matrix<double> exportAsMatrix() const;

/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PROTECTED MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 protected:


/* ---------------------------- */
/* PROTECTED METHODS */
/* ---------------------------- */
  virtual double extremal(int coord, bool max) const;


/* ---------------------------- */
/* PROTECTED STRUCTURES */
/* ---------------------------- */
//  QList<cpedsPoint2D> _ps; // points set

  mscsVector<double> _vals;
  
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


 };
#endif /* CPEDSPOINTSET2D */




















#ifndef CPEDSPOINTSET3D
#define CPEDSPOINTSET3D

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
#include <QtCore/qlist.h>
#include <vector>
#include "cpeds-math.h"
#include "cpeds-point3d.h"
#include "mscsVector.h"

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */


/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class cpedsPointSet3D
  \brief Encapsulates a set of cpeds point objects in 3D
  \details

  \date 2009/11/06 12:51:40
  \author Bartosz Lew
*/
//class cpedsPointSet3D  : public QList<cpedsPoint3D> {
class cpedsPointSet3D  : public mscsVector<cpedsPoint3D> {


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
  cpedsPointSet3D();
  cpedsPointSet3D(const QList<cpedsPoint3D> &q);
  cpedsPointSet3D(const cpedsPointSet3D &q);
  cpedsPointSet3D(const cpedsPointSet2D &q);
  //cpedsPointSet3D(long N, double *lon, double *lat);
  cpedsPointSet3D(long N, double *lon, double *lat, double* z);
  cpedsPointSet3D(long N, double *lon, double *lat, double z=1, bool deleteInside=true);
  /*!
	\brief constructor
	\details 
	@param N - number of points to generate
	@param len - longitudes
	@param lat - latitudes
	@param z - z values (if NULL - then 0 is set)
	@param v - values (if NULL - then values are not being set)

	\date Aug 8, 2012, 2:16:29 PM
	\author Bartosz Lew
*/
  cpedsPointSet3D(long N, double *lon, double *lat, double* z, double* v, bool deleteInside=true);
  cpedsPointSet3D(const vector<cpedsPoint3D> &v);
  cpedsPointSet3D(const mscsVector<cpedsPoint3D> &v);
  cpedsPointSet3D(const mscsVector<double>& x,const mscsVector<double>& y, const mscsVector<double>& z);
  cpedsPointSet3D(const mscsVector<double>& x,const mscsVector<double>& y, const mscsVector<double>& z, const mscsVector<double>& v);

  virtual ~cpedsPointSet3D();

/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */
  cpedsPointSet3D& append(const cpedsPoint3D& p) { push_back(p); return *this; }
  cpedsPointSet3D& append(const cpedsPoint3D& p, const double val) { push_back(p); values().push_back(val); return *this; }
  cpedsPointSet3D& append(const cpedsPointSet3D& ps) { for (unsigned long i = 0; i < ps.size(); i++) { append(ps[i]); }  return *this;  }
  cpedsPointSet3D& append(long N, double *x, double *y, double* z, double* v, bool deleteInside=true);
  
  
  long count() const { return size(); }
  cpedsPoint3D value(long i) { return cpedsPoint3D(at(i)); }
  cpedsPoint3D& point(long i) { return at(i); }
  double val(long i) { return _vals[i]; }
  double val(long i) const { return _vals[i]; }
  
  void clear() { vector<cpedsPoint3D>::clear(); _vals.clear(); }
  
  mscsVector<double>& values() { return _vals; }
  const mscsVector<double>& values() const { return _vals; }
  mscsVector<double>* valuesPtr() { if (_vals.size()==0) return NULL; return &_vals; }
  const mscsVector<double>* getValuesPtr() const { if (_vals.size()==0) return NULL; return &_vals; }
  void setSize(long n) { mscsVector<cpedsPoint3D>::setSize(n); _vals.setSize(n); }
  
  /*!
	\brief returns array with coordinates
	\details 
	@param coord - 0 - x, 1-y, 2-z
	@param size - pointer to an allocated long where the size of the returned array is stored.
	@return pointer to the newly allocated array with coordinates

	\date Feb 2, 2012, 9:01:50 PM
	\author Bartosz Lew
   */
  double* getVals(int coord, long* size) const;
  /*!
    \brief get the C array of X coordinats
    \details
    @param size - pointer which will be filed with the number of elements of the set
    @return
  */
  double* getXvals(long* size) const;
  //! see above
  double* getYvals(long* size) const;
  //! see above
  double* getZvals(long* size) const;

  void getRanges(double &minx, double &maxx, double &miny, double &maxy, double &minz, double &maxz) const;
  double lengthX() const;
  double lengthY() const;
  double lengthZ() const;
  void shiftX(double delta);
  void shiftY(double delta);
  void shiftZ(double delta);
  void makePeriodic(double periodX,double periodY,double periodZ);
  //! get the minimal value of in X direction from the point set
  double minX() const;
  //! get the maximal value of in X direction from the point set
  double maxX() const;
  //! get the minimal value of in Y direction from the point set
  double minY() const;
  //! get the maximal value of in Y direction from the point set
  double maxY() const;
  //! get the minimal value of in Z direction from the point set
  double minZ() const;
  //! get the maximal value of in Z direction from the point set
  double maxZ() const;

  //! calculate a dispersion of points for the requested coordinate
  double disp(int coord) const;
  //! get X dispersion of the points (rms value)
  double dispersionX() const;
  //! get Y dispersion of the points (rms value)
  double dispersionY() const;
  //! get Z dispersion of the points (rms value)
  double dispZ() const;

  void multiply(double val);
  void add(double val);

  virtual const cpedsPointSet3D& operator=(const cpedsPointSet3D &rhs);
  virtual const cpedsPointSet3D& operator=(const QList<cpedsPoint3D> &rhs);
  virtual const cpedsPointSet3D& operator=(const cpedsPointSet2D &rhs);
  virtual const cpedsPointSet3D& operator+=(const cpedsPointSet3D &rhs);
  virtual cpedsPointSet3D operator+(const cpedsPointSet3D &rhs);

  /*!
    \brief assigns values to the points in the set
    \details
    @param coord - defines the which coordinate should be used to assign the values to the newly generated points.
    The binary system is used for numbering the dimentions. \n
    1 - X axis \n
    2 - Y axis \n
    4 - Z axis \n
    3 - X and Y axis \n
    5 - X and Z axis \n
    6 - Y and Z axis \n
    7 - X and Y and Z axis \n
    The values to be assigned are not used for assignment for axes other than chosen for the value assignment

    @return The newly generated set is returned by reference

    \date 2009/11/10 23:56:32
    \author Bartosz Lew
  */
  const cpedsPointSet3D& setVals(int coord, double vX, double vY, double vZ );
  const cpedsPointSet3D& set(int i, double x, double y, double z, double v);
  const cpedsPointSet3D& set(int i, const cpedsPoint3D p, double v);
  //! sets the values according to the given vector
  void setVals(mscsVector<double> v);
  void setVals(long N, double* v, bool deleteInside=false);
  

  /*!
    \brief generates n points in the set and assigns values to them
    \details
    @param n - number of points to generate

    @return The newly generated set is returned by reference

    \date 2009/11/10 23:58:45
    \author Bartosz Lew
  */
  const cpedsPointSet3D& generatePoints(long n, double vX, double vY, double vZ);

//  /*!
//	\brief generate a random set of points structured hierarchically
//	\details 
//	@param
//	@return
//
//	\date Nov 4, 2015, 8:31:04 PM
//	\author Bartosz Lew
//   */
//  const cpedsPointSet3D& generateRandomPointsHierarchical(long nLevels, long nPtsLev, subDomain_region_t sd);

  const cpedsPointSet3D& generateRandomPoints(long nPts, subDomain_region_t sd);
  
  void print() const;
  /*!
	\brief calculates an prints the current ranges
	\details 
	This method calculates the ranges in X,Y,Z directions over which this points sets spans and prints the result to stdout.

	\date Mar 20, 2013, 9:38:49 AM
	\author Bartosz Lew
   */
  void printRanges(string comment="") const;
  virtual void save(string fname) const;
  /*!
	\brief loads points from file
	\details 
	@param fname - input file name with point set data

	\date Feb 2, 2012, 10:22:18 PM
	\author Bartosz Lew
   */
  virtual void load(string fname);

  cpedsPointSet3D getPoints() const { return cpedsPointSet3D(*this); }

  //! exports the points as a column of points three coordinates in each of the rows
  virtual const matrix<double> exportAsMatrix() const;

  /*!
    \brief exports the directions as 2-D field matrix with z values stored in them.
    \details
    @param sizeX - number of columns in the matrix
    @param sizeY - number of rows in the matrix
    @return field in a form of matrix

    The x and y coordinates information is not stored but the distribution of points in the matrix makes a map of points corresponding
    to those stored in the point set.

    There is no free space around the projected field in the resulting matrix - i.e there will be some data
    at first/last row/column in the matrix. Increasing resolution will increase the number of cells in between the extremal values
    and the sparsity of the matrix

    \date 2010/03/11 12:20:33
  */
  virtual const matrix<double> exportAsField(long sizeX, long sizeY) const;
  //!
  /*!
    \brief same as above but also populates the mask matrix with id numbers of the points from the set
    \details
    If some cell in the matrix doesn't have its counterpart in the point set, then the mask has -1 value in that cell
    \note The mask will only hold the information on the last id number of a point that fell in particular cell.

    \date 2010/03/09 20:52:25
  */
  virtual const matrix<double> exportAsField(long sizeX, long sizeY, matrix<long>& mask) const;

//  virtual const mscsFunction3dregc exportAsField3d(long sizeX, long sizeY, matrix<long>& mask) const;

  /*!
    \brief exports the set as a matrix and fills the empty cells with interpolated values
    \details
    @param mask - holds information on which cells were interpolated (-1 value)
    and which have a counterpart in the point set (valid number of the point in the set)
    @param psint - reference to the interpolated point set. This object will hold the point set with the interpolated
    points kept in the interpolated matrix generated from the original point set, stored in the row-major ordering.

    The points are snapped to the matrix grid.

    The interpolation of the points is done row by row (column by column) for the first coordinate (for the second
    coordinate), using data from the same row (column). For linear interpolation at least two points in the same
    row (column) are needed. If more than two points are found then cubic interpolation is used. If there's only
    one point in given row (column) then first (second) coordinate is sought for in second row (column) and so forth,
    until the sufficient number of points for interpolation are found. If not sufficient number of points is found then
    mask flag is set to -2 to indicate that this particular point was not interpolated and the corresponding psint
    value is set to the default zero value.
    Flag value of -4 indicates that both coordinates were successfully interpolated;
    Flag value of -3 indicates that only first coordinates was successfully interpolated;

    \warning This code is experimental and should be checked before use. Currently the interpolation is done only along the rows
    of the field
    @return populated interpolated matrix
    \date 2010/03/09 21:36:57
  */
  virtual const matrix<double> exportAsInterpolatedField(long sizeX, long sizeY, matrix<long>& mask, cpedsPointSet3D& psint) const;

  virtual vector<cpedsPoint3D> exportAsVector() const;
  
  /*!
    \brief generates a set from the m matrix.
    \details
    @param m - matrix reference used as a source to generate the set
    @param refRow - reference row
    @param refCol - reference column. Look below for more info.
    @param p - real coordinates of the element in the refRow,refCol cell.
    @param resolution - resolution parameter of the matrix grid

    The coordinates of the points in the newly generated set are calculated based on the reference cell in the matrix refRow and refCol which is known
    to hava a true coordinates as defined by point p. The remaining points coordinates are calculated by linear transformation using the matrix grid resulution
    information.

    \date 2009/12/09 15:53:46
  */
  virtual const cpedsPointSet3D field2set(const matrix<double>& m, long refRow, long refCol, cpedsPoint3D p, cpedsPoint3D resolution) const;

  virtual const cpedsPoint3D massCenter() const;
  virtual cpedsPointSet3D calculateGravitationalPotential(double gravSoft);
  
  cpedsPoint3D getMinValuePoint();
  
  
/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PROTECTED MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 protected:


/* ---------------------------- */
/* PROTECTED METHODS */
/* ---------------------------- */
  virtual double extremal(int coord, bool max) const;


/* ---------------------------- */
/* PROTECTED STRUCTURES */
/* ---------------------------- */
//  QList<cpedsPoint2D> _ps; // points set
  mscsVector<double> _vals;

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
  /* double _Xmin, _Xmax, _Ymin, _Ymax, _Zmin, _Zmax; //!< extremal values in the set */


 };
#endif /* CPEDSPOINTSET3D */

