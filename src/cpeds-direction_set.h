/*!
  \file encapsulates a set of cpeds point objects in 2D
*/

#ifndef CPEDSDIRECTIONSET
#define CPEDSDIRECTIONSET

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
//#include <QtCore/qlist.h>
#include "cpeds-direction.h"

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */


/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class cpedsDirectionSet
  \brief Encapsulates a set of cpeds direction objects
  \details

  \date 2009/11/06 12:51:40
  \author Bartosz Lew
*/
class cpedsDirectionSet : public QList<cpedsDirection> {


/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PUBLIC MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 public:

	using  QList<cpedsDirection>::append;
	using  QList<cpedsDirection>::size;
	using  QList<cpedsDirection>::count;
	using  QList<cpedsDirection>::clear;
	using  QList<cpedsDirection>::value;
	using  QList<cpedsDirection>::at;
	using  QList<cpedsDirection>::first;
	using  QList<cpedsDirection>::last;
	using  QList<cpedsDirection>::operator[];
	using  QList<cpedsDirection>::operator<<;
	using  QList<cpedsDirection>::takeFirst;
	using  QList<cpedsDirection>::takeLast;

/* ------------- */
/* CLASS FRIENDS */
/* ------------- */


/* ---------------------------- */
/* CONSTRUCTORS AND DESTRUCTORS */
/* ---------------------------- */
  cpedsDirectionSet();
  cpedsDirectionSet(const QList<cpedsDirection> &q);
  cpedsDirectionSet(const cpedsDirectionSet &q);
  cpedsDirectionSet(long N, double *lon, double *lat);
  cpedsDirectionSet(long N, double *lon, double *lat, double *val);
  cpedsDirectionSet(const QList<double>& lon, const QList<double>& lat, const QList<double>& val);
  ~cpedsDirectionSet();

/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */

  double* getLonVals(long* size) const;
  double* getLatVals(long* size) const;
  double* getVals(long* size) const;

  void getRanges(double &minLon, double &maxLon, double &minLat, double &maxLat) const;
  double minLon() const;
  double maxLon() const;
  double minLat() const;
  double maxLat() const;
  /*!
	\brief calculate the mean direction out of the set
	\details 
	@return mean direction in radians
	
	It is assumed that the set is given in radians
	
	The mean is calculated as the mean vector after converting coordinates to rectangular coordinates

	\date Jun 9, 2017, 10:28:38 AM
*/
  cpedsDirection mean() const;

  virtual cpedsDirectionSet& operator=(const cpedsDirectionSet &rhs);
  virtual cpedsDirectionSet& operator=(const QList<cpedsDirection> &rhs);

  /*!
    \brief sets the length of the set to n using default values;
    \details
    If n is larger than the current length of the set then new directions will be appended.
    If n is smaller then the current length of the set then  the elements will be deleted from the beginning of the set.
  */
  void setLength(long n);

  void sort();


  void print() const;
  long save(string fname, string how="", bool withVals=false) const;
  /*!
    \brief loads set from file
    \details
    @param how - "binary" indicates that the file is binary, "float" - 4byte binary file
    @para withVals - true: three columns are assumed, false: two columns are assumed
    @return load status: 0 - ok

    \date 2010/03/12 00:35:21
    \author Bartosz Lew
  */
  long load(string fname, string how="", bool withVals=false);

  /*!
	\brief exports the list as 2 dimentional matrix with 2 or 3 columns.
	\details
	@param withVals: true - 3 columns (lon, lat, val), false - two columns export (lon, lat) only
	@return

	\date May 17, 2010, 3:01:03 PM
	\author Bartosz Lew
  */
  const matrix<double> exportAsMatrix(bool withVals=false) const;
  /*!
    \brief imports data from matrix m
    \details
    @param withVals - only used in case of column vector matrices\n
    if true - every triplet of rows will be treated as lon,lat,val for generating the direction set.\n
    if false - every duplet of rows will be treated as lon,lat for generating the direction set. the value will be assigned zero value\n
  */
  void importMatrix(const matrix<double>& m, bool withVals=false) ;


/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PROTECTED MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 protected:


/* ---------------------------- */
/* PROTECTED METHODS */
/* ---------------------------- */


/* ---------------------------- */
/* PROTECTED STRUCTURES */
/* ---------------------------- */
//  QList<cpedsDirection> _ps; // points set


/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PRIVATE MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 private:

/* ---------------------------- */
/* PRIVATE METHODS */
/* ---------------------------- */
  double extremal(bool lon, bool max) const;


/* ---------------------------- */
/* PRIVATE STRUCTURES */
/* ---------------------------- */


 };
#endif /* CPEDDIRECTIONSET */

