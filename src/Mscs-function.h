/*!
  \file Defines a 1-dimentional scalar function class f(x) and the like
*/
#ifndef MSCS_FUNCTION
#define MSCS_FUNCTION

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */
#include <QtCore/QtAlgorithms>
#include <QtCore/QPointF>
#include <QtCore/QList>
#include <string>
#include <iostream>
#include <sstream>
#ifndef NO_HDF5
#include <hdf5.h>
#endif
#include "Mscs-object.h"
#include "cpeds-math.h"
#include "cpeds-list.h"
#include "cpeds-rng.h"





extern string _HDF5currentGroup_tmp;
extern vector<string> _HDF5_datasets;

#ifndef NO_HDF5
herr_t op_func (hid_t loc_id, const char *name, const H5L_info_t *info,void *operator_data);

/*
 * Define operator data structure type for H5Literate callback.
 * During recursive iteration, these structures will form a
 * linked list that can be searched for duplicate groups,
 * preventing infinite recursion.
 */
struct opdata {
		unsigned        recurs;         /* Recursion level.  0=root */
		struct opdata   *prev;          /* Pointer to previous opdata */
		haddr_t         addr;           /* Group address */
};
int group_check (struct opdata *od, haddr_t target_addr);
//static const int HDF5_stringAttributeMaxLength=500;
#endif





using namespace std;
/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class mscsFunction
  \brief Encapsulates a 1-D scalar function f(x) functionality
  \details
  The function is realized as a queue of reala valued points.
	Some complex or vector functionality is also supported

  \date 2009/05/28 15:50:18
  \author Bartosz Lew
*/
class mscsFunction : public mscsObject {

 private:
  /*!
    \struct
    \brief a ranges within which a function is embedded
    \details

    \date 2009/05/28 20:15:55
    \author Bartosz Lew
  */
  typedef struct {
    double xmin; //!< minimal value of the function argument
    double xmax; //!< maximal value of the function argument
    double ymin; //!< minimal value of the function value
    double ymax; //!< maximal value of the function value
    long ixmin; //!< index of the minimal value of the function argument
    long ixmax; //!< index of the maximal value of the function argument
    long iymin; //!< index of the minimal value of the function value
    long iymax; //!< index of the maximal value of the function value
  } frange;


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

  //! an empty constructor
  mscsFunction();

  //! an empty constructor giving a name to this function
  mscsFunction(const string name, cpeds_VerbosityLevel verbosity=CPEDS_defaultVerbosityLevel);

  /*!
    \brief initiates the function with the array of cpeds_points of length num
    \details
    @param name - name of the function
    @param t - pointer to an array of cpeds_points
    @param num - size of the array

	  If t is NULL and num is given along, then an array of cpeds_points is generated and initiated with zeros
	and used for function initialization.

    \date 2009/05/28 19:20:09
    \author Bartosz Lew
  */
  mscsFunction(const string name, cpeds_point* t, long num);

  /*!
    \brief initiates the function with the array of points given in x and y arrays of length num
    \details
    @param name - name of the function
    @param x - pointer to an array of arguments
    @param y - pointer to an array of values
    @param num - size of the array

    \date 2009/05/28 19:20:09
    \author Bartosz Lew
  */
  mscsFunction(const string name, double *x, double *y, long num, cpeds_VerbosityLevel verbosity=CPEDS_defaultVerbosityLevel);


  /*!
    \brief initiates the function with the array of points given in x and y objects
    \details
    @param name - name of the function
    @param x - QList object reference with function arguments
    @param y - QList object reference with function values

    \date 2009/05/28 19:20:09
    \author Bartosz Lew
  */
  mscsFunction(const string name, const QList<double>& x, const QList<double>& y);

  //! a cloning constructor
  mscsFunction(const mscsFunction& fn);

  //! destructor
  virtual ~mscsFunction();

/****************************************/
/* /\* ---------------------------- *\/ */
/* /\* PUBLIC METHODS *\/	        */
/* /\* ---------------------------- *\/ */
/****************************************/

  /**************************/
  /* RANGE AND DATA TESTERS */
  /**************************/

  //! returns the maximal value of the function
  double getMaxValue() const;
  //! returns the minimal value of the function
  double getMinValue() const;
  //! returns the maximal argument of the function
  double getMaxArg() const;
  //! returns the maximal argument of the function
  double getMinArg() const;
  //! returns index of the maximal value of the function
  long getMaxValueIdx() const { return range.iymax; }
  //! returns the index of the  minimal value of the function
  long getMinValueIdx() const  { return range.iymin; }
  //! returns index of the  maximal argument of the function
  long getMaxArgIdx() const { return range.ixmax; }
  //! returns index of the  maximal argument of the function
  long getMinArgIdx() const { return range.ixmin; }


  /*****************************/
  /* VARIOUS FUNCTION HANDERLS */
  /*****************************/

  /*!
    \brief sets i'th point of the function
    \details
    and updates the information about function ranges.
    If you want do substitution faster then use f(i)=... function
    \author Bartosz Lew
  */
  virtual void setf(long i, double x, double y);
  /*!
    \brief sets i'th value of the function to y
    \details
    and updates the information about function ranges.
    If you want do substitution faster then use f(i)=... function
    \author Bartosz Lew
  */
  void setf(long i,  double y);
  //! sets the entire function to value val
  const mscsFunction& setf(double val);
  //! sets i'th argument of the function to x
  void setarg(long i,  double x);

  //!
  /*!
    \brief returns the value of an argument which is  closest to the argument x;
    \details

    if the function isn't allocated zero is returned
  */
//  double f(double x) const;
  //! same as above but the index of the closest argument is returned under i pointer
  double f(double x, long* i=NULL) const;
  //! same as above but the value at x is interpolated 
  double finter(double x, string type="linear") const;
//  double finter(double x, string type="linear");
  //! returns the value reference at the i'th position; use this if you are sure that the i'th value exists in the function
  double& f(long i) { return _f[i].ry(); }
  //! returns the value reference at the index i treated in the python array index style. This allows to get function values from the end by using negative indexes.
  double& fpsty(long i);
  //! returns the value reference at the i'th position using periodic boundary conditions
  double& fp(long i) { return _f[i % _f.size()].ry(); }
  //! returns function value without a check for validity of range i. functionality changed since Jun 28, 2012, 5:17:57 PM. Now Y has the previous functionality
  double y(long i) const { return _f[i].y(); }
  //! same as getY
  double Y(long i) const;
  //! returns function value after a check for validity of range i; if non valid i is given then defult value of 0 is returned
  double getY(long i) const;
  //! returns the argument value at i'th position
  double getX(long i) const;
  double getx(long i) const { return _f[i].x(); }
  double& X(long i) { return _f[i].rx(); }
  //! Returns the number of points that build up this function
  long pointsCount() const;
  
  /*!
	\brief get i'th point from the functio as const reference
	\details 
	@param i - index of the point
	@return QPointF object is returned
	
	No checking for wrong indexes is done/

	\date Mar 3, 2011, 9:23:06 PM
	\author Bartosz Lew
  */
  const QPointF& get(long i) const { return _f[i]; }
  /*!
	\brief return a range of the data from the function
	\details 
	@param from - starting index
	@param to - ending index
	@return a new function with copied range is returned
	
	The copied range is taken as: [from,to)
	If to is given larger than the number of points in the function it is set to the pointsCount() value.

	\date Mar 3, 2011, 9:21:32 PM
	\author Bartosz Lew
  */
  mscsFunction get(long from, long to);
  //! checks whether the i is in a viable range for usage with this function
  bool argOK(long i) const;


  //! sets function values to v for function values smaller than th
  const mscsFunction& setBelow(double th, double val);
  //! sets function values to v for function values larger than th
  const mscsFunction& setAbove(double th, double val);
  //! sets the function value to value val within range <from,to>
  mscsFunction& setRange(double from, double to, double val);
  

  /*!
    \brief returns a newly allocated array with function argument values
    \details
    @param n - number of arguments to extract; default=0 means all will be extracted
    @return - an array of double values; if there are no points in the function
    the return value is NULL
  */
  double* extractArguments(long n=0) const;
  double* extractArguments(long imin,long imax) const;

  cpedsList<double> toXList() const;

  /*!
    \brief exports the function arguments as a list
    \details
    @param
    @return QList object list  of doubles with arguments
    Use this method for compound commands instead of exractArguments to avoid memory leaks when the poiner is not stored

    \date 2010/01/29 23:14:26
    \author Bartosz Lew
  */
  const QList<double> getArguments() const;

  /*!
    \brief returns a newly allocated array with function values
    \details
    @param n - number of values to extract; default 0 - means all will be extracted
    @return - an array of double values; if there are no points in the function
    the return value is NULL
  */
  double* extractValues(long n=0) const;
  double* extractValues(long imin,long imax) const;
  cpedsList<double> toYList() const;

  
  /*!
    \brief exports the function values as a list
    \details
    @param
    @return QList object list  of doubles with values
    Use this method for compound commands instead of exractValues to avoid memory leaks when the poiner is not stored

    \date 2010/01/29 23:14:26
    \author Bartosz Lew
  */
  const QList<double> getValues() const;

  QPointF* extractPoints() const;
  void importPoints(QPointF* t, long n, bool reverseOrder=false);

  /*!
    \brief adds points to this function from the X,Y array of size size
    \details
    The existing points in the function are not deleted.
	This function can also be used when only X or Y parts need to be imported. 
	In such case if X=NULL arguments will be populated with numbers starting from 0,1,2 ...
	If Y=NULL then values are set to zero.

    \date 2009/05/28 21:09:05
    \author Bartosz Lew
  */
  void importFunction(double* X, double* Y, long size, bool deleteAfter=false);
  void importFunction(mscsVector<double> X, mscsVector<double> Y);
  /*!
    \brief adds points to this function from the array of cpeds_points of size size
    \details
    The existing points in the function are not deleted.

    \date 2009/05/28 21:09:05
    \author Bartosz Lew
  */
  void importFunction(cpeds_point* p, long size);


  /*!
	\brief sets a function using list
	\details
	@param x list of arguments
	@return *this
	
		This operation clears the function first.
		
	changed on: Jun 25, 2012, 7:23:44 PM
		From now on this operation does not clear the existing function and does not change its size.
		Use setPointsNum first if you are not sure about the sizes of the function and x list.
	  
	\date May 21, 2010, 2:17:09 PM
  */
  const mscsFunction& setX(const QList<double> x);

  const mscsFunction& setX(double v);

  /*!
	\brief sets a function using lists
	\details
	@param y list of function values
	@return *this

	  This operation does not clear the function.
	  
	changed on: Jun 25, 2012, 7:23:44 PM
		From now on this operation does not clear the existing function and does not change its size.
		Use setPointsNum first if you are not sure about the sizes of the function and x list.
	  
	\date May 21, 2010, 2:17:09 PM
  */
  const mscsFunction& setY(const QList<double> y);

  /*!
    \brief adds new points (or removes points) so that the total number of points was i
    \details
    The new points will be zeros.\n
    In case of deleting points, the points will be deleted from the end of the queue.
    @param num - total number of required points
  */
  void setPointsNum(long num);

  cpeds_point getPoint(long i) const { cpeds_point p; p.x=getX(i); p.y=f(i); return p; }
  QPointF getQPoint(long i) const { return _f[i]; }
  void setQPoint(long i,QPointF qp) { _f[i]=qp; }
  QPointF last() { return _f.last(); }

  /*!
    \brief Adds a new point to the function
    \details
    It doesn't update the the ranges information.

    \date 2009/05/28 21:38:45
  */
//  void addPoint(double x, double y) {  _f.append(QPointF::QPointF(x,y)); } // commented out when switching to the 64bit system; todo need to find out if it still compiles without problems on 32-bit installations in standard configuration
  void addPoint(double x, double y) {  _f.append(QPointF(x,y)); }

  /*!
    \brief This and following function adds a new point to the function
    \details
    This does not looks whether the ordering of the arguments is in any particular
    order. Keep that in mind; Write another method for doing so.
    It keeps the ranges of the function up to date.

    \date 2009/05/28 21:38:45
    \author Bartosz Lew
  */
  void newPoint(double x=0, double y=0);
  void newPoint(cpeds_point p);
  void newPoint(const QPointF& p);

  /*!
    \brief This and following function inserts a new point to the function
    \details
    @param x - argument 
    @param y - value
    @param i - the index at which a new point will exist.
    The points starting with i'th + one will by pushed forward
    This does not looks whether the ordering of the arguments is in any particular
    order. Keep that in mind; Write another method for doing so.
    @param N - number of points to insert

    \date 2009/05/28 21:38:45
    \author Bartosz Lew
  */
  void insertPoint(double x, double y,long i,long N=1);
  void insertPoint(cpeds_point p,long i);

  //! removes the i'th point from the function
  void deletePoint(long i);
  //! removes the last point from the function
  void deletePoint();

  /*!
	\brief removes points from this function that are in function p
	\details 
	@param p - function with points to be removed from this function
	@return this function returns this

	\date Feb 7, 2011, 12:43:17 PM
	\author Bartosz Lew
  */
  mscsFunction& remove(const mscsFunction& p);
  
  /*!
	\brief removes points that have value v
	\details 
	@param v - value by which all points will be removed
	@return returns this

	\date Jun 25, 2012, 6:22:33 PM
	\author Bartosz Lew
   */
  mscsFunction& removeValue(double v);
  mscsFunction& removeSmaller(double v);
  mscsFunction& removeLarger(double v);

  mscsFunction& removePoints(cpedsList<long> l);
  mscsFunction& removeNans();

  void clearFunction();

  /*!
	\brief export values as matrix
	\details 
	@return

	\date Dec 13, 2010, 8:28:41 PM
	\author Bartosz Lew
  */
  matrix<double> valuesToVector(bool vertical);
  
  /**************/
  /* IO METHODS */
  /**************/
  
  /*!
	\brief load data from file as a function
	\details 
	@param filename - name of the file
	@param commentedFile - if true then the read will disregard the lines starting from #, if false then it is assumed that all rows start with numbers
	@param colx - column number in a multicolumn file that should be regarded as function arguments: if -1 then it is assumed that only y values are in the single column file.
	@param coly - column number in a multicolumn file that should be regarded as function values: if -1 then it is assumed that only x values are in the single column file.
	@param Nrows - number of rows to read; Default -1 - means all rows are read
	@param skipRows - number of rows the be skipped before loading the data
	@return cpeds file io in status code

	\date Nov 17, 2010, 7:02:34 PM
	\author Bartosz Lew
  */
  virtual cpedsStatusCodes load(string filename, bool commentedFile=false, int colx=0, int coly=1, long long Nrows=-1, long long skipRows=0);
  virtual cpedsStatusCodes save(string filename, bool append=false) const;
  void print(bool printPointNo=true) const;
  void printRanges() const;

#ifndef NO_HDF5
		/*!
			\brief save function into hdf5 file
			\details 
			@param hdf5 file name
			@param datasetName - name for the new dataset - the file should not contain this dataset
			@return exit status code
		
			\date Mar 23, 2012, 4:44:30 PM
			\author Bartosz Lew
		*/
		cpedsStatusCodes saveHDF5(string filename, string datasetName, string xUnit="",string yUnit="",string xComment="",string yComment="");
		cpedsStatusCodes loadHDF5(string filename, string datasetName);
		vector<string> getHDF5dataSets(string filename);
		void setHDF5_scalarStringAttribute(string fname, string dsetName, string attributeName, string attributeValue);
		void setHDF5_scalarDoubleAttribute(string fname, string dsetName, string attributeName, double value, string attributeComment="");
		string getHDF5_stringAttribute(string fname, string dsetName, string attributeName, int* errCode);

		
		void initiatie_hdf5_params();
		static const int HDF5_stringAttributeMaxLength=500;
#pragma omp threadprivate(HDF5_stringAttributeMaxLength)

#endif

  
  
#ifdef PGPLOT
  int plotPG(int plottype=0);   //! plots the function using pgplot
#endif


  /***********************/
  /* FUNCTION OPERATIONS */
  /***********************/

  /*!
    \brief sorts the function by arguments in ascending ordering
    \details The function structure is saved of course :)
    This routine is useful when defining the function as a set of
    arbitrary points without any particular care about the ordering of
    arguments.

    \date 2009/06/04 21:29:34
    \author Bartosz Lew
  */
  const mscsFunction& sortFunctionArgAscending();
  const mscsFunction& sortFunctionArgDescending();

/*   power_spectrum* bin_function(long lmin_loc, long lmax_loc, long *bintabs, long* bintab, double w); */
/*   void set_Cl_multipole_range(long l1, long l2, double val); */

  bool isPositive() const;
  
  void Xadd(long i, double v); //<! adds a value to i'th argument
  void xadd(long i, double v); //<! adds a value to i'th argument
  void Xadd(double* v); //<! adds an array of values to arguments
  void Xsubtract(long i, double v); //<! subtracts a value from i'th argument
  void Xmultiply(long i, double v); //<! multiplies i'th argument by value
  void Xdivide(long i, double v); //<! divides i'th argument by a value
  void Yadd(long i, double v); //<! adds a value to i'th function value if the argument exists
  void yadd(long i, double v); //<! adds a value to i'th function value
  void Ysubtract(long i, double v); //<! subtracts a value from i'th function value  if the argument exists
  void ysubtract(long i, double v); //<! subtracts a value from i'th function value
  void Ymultiply(long i, double v); //<! multiplies i'th function value by value  if the argument exists
  void ymultiply(long i, double v); //<! multiplies i'th function value by value
  void Ydivide(long i, double v); //<! divides i'th function value by a value  if the argument exists
  void ydivide(long i, double v); //<! divides i'th function value by a value
  void asVectorMultiply(long i,double v); //!< multiplies i'th vector by v
  double modulus(long i); //!< retunrs modulus of i'th complex number
  double sumX() const; //<! returns the sum of all x values
  double sumY() const; //<! returns the sum of all y values

  mscsFunction& add(double v); //<! adds a value to all function values
  mscsFunction& subtract(double v); //<! subtracts a value from all function values
  mscsFunction& multiply(double v); //<! multiplies by value all function values
  mscsFunction& divide(double v); //<! divides by a value all function values
  mscsFunction& power(double a); //<! derives power a of the function
  mscsFunction& absoluteValue(); //<! calculates absolute value of the function
  mscsFunction& sqroot(); //<! derives square root a of the function (returns 0 on negative function values)
  mscsFunction& logX(double base=10); //<! converts X values into a log X values
  mscsFunction& logY(double base=10); //<! converts Y values into a log Y values
  mscsFunction& lnY(); //<! convenience function  that converts Y values into a ln Y values
  mscsFunction& exponent(); //<! calculates exponent of the current function
  mscsFunction& exponentX(double base=10); //<! calculates exponent of the arguments
  mscsFunction& exponentY(double base=10); //<! calculates exponent of the function
  mscsFunction& add(const mscsFunction& F, long startFrom=0); //!< adds a function according to the index of points stored starting from startFrom index
  mscsFunction& add(const double* t); //!< adds to function values the values from the array t whic must be at least of the size of this function
  double vectorSum(); //!< adds function X and Y values as vectors and returns vector length
  mscsFunction& vectorAdd(const mscsFunction& F); //!< adds function F as vector treating X and Y values as vector coordinates
  mscsFunction& vectorAdd(const double* v); //!< adds function F as vector treating X and Y values as vector coordinates
  mscsFunction& vectorSubtract(const mscsFunction& F); //!< subtracts function F as vector treating X and Y values as vector coordinates
  mscsFunction& scalarMultiply(const mscsFunction& F); //!< multiplies function F as vector treating X and Y values as vector coordinates
  mscsFunction& subtract(const mscsFunction& F);  //!< subtracts a function according to the index of points stored
  mscsFunction& multiply(const mscsFunction& F); //!< multiplies by a function according to the index of points stored
  mscsFunction& divide(const mscsFunction& F); //!< divides by a function according to the index of points stored
  mscsFunction& divide(const cpedsList<double>& l); //!< divides function values by a list element-wise 
  mscsFunction& divide(const cpedsList<long>& l);
  
  mscsFunction& epoweri(double phi, bool phiR=false); //!< complex multiplication by e^(i*phi) - treats arguments as real part and values as imaginary part of the complex number, unless Rphi=true, in which case x corresponds to phase of the complex number of y corresponds to modulus
  mscsFunction& epowerii(long i, double phi, bool phiR=false); //!< complex multiplication of i'th point in the function by e^(i*phi) - treats arguments as real part and values as imaginary part of the complex number, unless Rphi=true, in which case x corresponds to phase of the complex number of y corresponds to modulus
  
  
  mscsFunction& X2Y(bool x2y=true); //!< copies arguments onto the function values

  mscsFunction& shiftX(double x); //!< shifts the function by x according to transformation: f(y) -> f(y-x) (negative x will shift towards smaller x values; positive x will shift towards larger values)
  /*!
	\brief transforms function points arguments to fit into range [from, to)
	\details
	@param acc - tolerance with which the points inside of the range but very close to the upper limit will be folded as if they were exactly on "to" value
	acc should be positive.

	The arguments x of the points from outside of the region [from, to] are transformed as:\n
	  z = z-from\n
	  x-> (z - [z/(to-from)] * (to-from)) + from\n
	  where [y] - is an integer part of y trucated towards zero\n

	@return this method returns "*this"

	\date May 19, 2010, 11:49:28 AM
  */
  mscsFunction& foldXinto(double from, double to, double acc=0);
  mscsFunction& scaleX(double a); //!< scales in the X direction by a; *this is returned

  //! deletes points outside of the <from,to> range
  mscsFunction& deleteOutsideOf(double from, double to);
  
  //! deletes points inside of the <from,to> range
  mscsFunction& deleteRange(double from, double to);
  //! deletes ranges specified by [x_i,y_i] points in fromto function
  mscsFunction& deleteRanges(mscsFunction fromto);
  /*!
	\brief shifts the function arguments and values by n cells periodically
	\details 
	@param n - amount of the shift, can be positive or negative
	@return returns this

	  This is fast because it does not rewrite tables but appends and deletes elemets to/from the list.
	\date Oct 29, 2010, 2:30:31 PM
	\author Bartosz Lew
  */
  mscsFunction& shift(long n);
  mscsFunction& shiftYwrtX(long n);
  
  /***********************/
  /* FUNCTION GENERATION */
  /***********************/
  /*!
	\brief generate sinus function
	\details 
	@param from - from X
	@param to - to X
	@param dx - dx
	@param T - period in units of X
	@param phi - offset in units of T
	@return

	\date Jul 22, 2010, 8:41:35 PM
	\author Bartosz Lew
  */
  mscsFunction& mkSin(double from, double to, double dx, double T=1.0, double phi=0.0, double A=1.0);
  mscsFunction& mkSquareWave(double from, double to, double dx, double T=1.0, double phi=0.0);
  mscsFunction& mkConst(double from, double to, double dx, double v=0);
  mscsFunction& mkTopHat(double from, double to, double dx, double v1, double from2, double to2, double v2=1);
  mscsFunction& mkGumbelDistr(double from, double to, double dx, double beta=1, double m=0, double tol=0);
  mscsFunction& mkWeibullDistr(double from, double to, double dx, double lambda, double k, double tol=0);
  mscsFunction& mkGauss(double from, double to, double dx, double A=1, double m=0, double s=1, double B=0, double tol=0);
  /*!
	\brief generates skew-gaussian function
	\details 
	@param alpha - skewness parameter: >0 gives positive skewness (right tail longer), <0 gives negative skewness (left tail longer)
	@param tol - tolerance
	@return

	The function is:
	gauss * (1.0-gsl_sf_erf(nsig/sq2))
	where
		gauss=A * exp(-(x-m)*(x-m)/(2*s*s));
	and
		nsig=(m-x)*alpha/s;
	

	\date Apr 10, 2013, 11:46:32 PM
	\author Bartosz Lew
   */
  mscsFunction& mkSkewGauss(double from, double to, double dx, double A=1, double m=0, double s=1,double alpha=0, double tol=0);
  mscsFunction& mkLog(double from, double to, double dx, double base=0, double tol=0);
  mscsFunction& mkExponent(double from, double to, double dx, double base, double tol=0);
  mscsFunction& mkLogSpace(double from, double to, long N, double base=0);
  mscsFunction& mkLine(double A, double B);
  mscsFunction& mkLine(double from, double to, double dx, double A=1, double B=0, double tol=0);
  mscsFunction& mkLine(double from, double to, double dx, double A, double B, double x0, double tol);
  mscsFunction& mkGaussCDF(double from, double to, double dx, double A=1, double m=0, double s=1);
  /*!
	\brief create a power law function
	\details 
	@param
	@return
	If the function is already initialized then the power law will be generated on the existing function arguments.
	and dk parameter will be ignored

	\date Nov 14, 2013, 1:14:58 PM
	\author Bartosz Lew
*/
  mscsFunction& mkPowerLaw(double from, double to, double dk, double A, double k0, double ns);
  mscsFunction& mkPowerLaw(double from, double to, double dk, double A, double k0, double k1, double ns, double B);
  mscsFunction& mkPowerLaw(double A, double k0, double ns);

  /*!
	\brief generates beta model of the NFW density profile
	\details 
	@param
	@return this

	n(r) = n0 ( 1+ (r/rc)^2 )^(-3/2*beta)
	
	param1 - n0
	param2 - rc
	param3 - beta

	\date Nov 28, 2013, 3:15:43 PM
	\author Bartosz Lew
   */
  mscsFunction& mkBetaModel(double from, double to, double dr, double n0, double rc, double beta);

  
  /*!
	\brief generates atmospheric pressure profile
	\details 
	@param from - from [km]
	@param to - to [km]
	@param dz - step [km]
	@param P0 - pressure at the reference level at altutude z0 [hPa]
	@param T0 - temperature at the reference level z0 [K]
	@param z0 - reference altutude [km]
	@param Lr - lapse rate [K/km]
	@param g0 - gravitational constant [m/s^2]
	@param mu - molar mass of the atmosphere [kg/mol]
	@return returns this

	\date Oct 4, 2012, 11:32:36 PM
	\author Bartosz Lew
*/
  mscsFunction& mkPressureProfile(double from, double to, double dz, double P0, double T0, double z0, double Lr=-6.5, double g0=9.80665, double mu=0.0289644);
  /*!
	\brief black body spectral energy density
	\details 
	@param from - initial frequency in Hz
	@param to - final frequency in Hz
	@param dnu - frequency increment
	@param black body temperature
	@return 8pi nu^2/c^3 * h nu/(exp(-h nu/(kT))-1) [J/m^3/Hz]

	\date Jan 2, 2016, 1:54:34 PM
	\author Bartosz Lew
   */
  mscsFunction& mkPlank(double from, double to, double dnu, double T);
  /*!
	\brief generate polynomial given polynomil coefficients
	\details 
	@param from
	@param to
	@param dx
	@param order - polynomial order: eg. 2 for polynomial of second order (parabola)
	@param a0 - 0'th order polynomial coefficient
	@param ... arbitrary number of polynomial coefficients.
	@return returns this

	\date Feb 9, 2012, 12:23:16 PM
	\author Bartosz Lew
   */
  mscsFunction& mkPolynomial(double from, double to, double dx, int order, double a0, ...);
  mscsFunction& mkPolynomial(double from, double to, double dx, cpedsList<double> a);
  mscsFunction& mkPolynomial(cpedsList<double> a);
  
  /*!
	\brief tabulates the Boltzmann distribution
	\details 
	@param mu - chemical potential
	@param T - temperature
	@return

	\date Nov 22, 2010, 3:39:13 PM
	\author Bartosz Lew
  */
  mscsFunction& mkBoltzmann(double from, double to, double dE, double T, double mu=0);


  

  /*!
	\brief generate gaussian white noise
	\details 
	@param N - number of points generated
	@param m - mean of the gaussian PDF
	@param s - stdev of the gaussian PDF
	@param seed - seed to be used for generation
	@param rng - random number generator to be used
	@param noisex - make noise on x axix
	@param noisey - make noise on y axis
 	@return this function is returned
 	
 	if noisex and noisey are chosen as:\n
 	true, false then the y values will be assigned 0 value\n
 	false, true - the x values will be enumerated from 0 to N-1
 	true true - both axes are randomized

	\date Dec 22, 2010, 3:09:32 PM
	\author Bartosz Lew
  */
  mscsFunction& mkGaussianNoise(long N, double m, double s, long seed=0, cpedsRNG* rng=NULL, bool noisex=false, bool noisey=true);

  /*!
	\brief generate chisq noise with NDOF=1
	\details 
	@param N - number of points generated
	@param m - mean of the gaussian PDF
	@param s - stdev of the gaussian PDF
	@param seed - seed to be used for generation
	@param rng - random number generator to be used
	@param noisex - make noise on x axix
	@param noisey - make noise on y axis
 	@return this function is returned
 	
 	if noisex and noisey are chosen as:\n
 	true, false then the y values will be assigned 0 value\n
 	false, true - the x values will be enumerated from 0 to N-1
 	true true - both axes are randomized

	This basically generates gaussian noise and takes the square of the output.

	\date Dec 22, 2010, 3:09:32 PM
	\author Bartosz Lew
  */
  mscsFunction& mkChisqNoise(long N, double m, double s, long seed=0, cpedsRNG* rng=NULL, bool noisex=false, bool noisey=true);

  
  /*!
	\brief generate random gaussian, real-space noise with requested power spectrum
	\details 
	@param A - amplitude of noise
	@param alpha - spectral indes of the power law noise: 0 - white, -1 - 1/f red noise, -2 - 1/f^2 noise, 1 - f noise (blue)etc.
	@param k0 - defines pivot scale
	@param kmin - minimal frequency probing the power spectrum
	@param kmax - maximal frequency probing the power spectrum
	@param N - number of points in the noise
	@param seed - seed for the RN sequence generation (only used if rng is not provided)
	@param powerSpectrum - reference to the object at which to append the power spectrum of the generated real signal. If the realSpace=false, the real space signal is not calculated at all, and if powerSpectrum==NULL then the 1/f noise power spectrum is appended to *this object
	@param realSpace - defines if the noise should be defined in real space of Fourier space. If true then Fourier transform is done on the generated data and returned.
	If false, then Fourier transform is not done and the generated power spectrum of the noise is returned via powerSpectrum parameter.
	@return this function returns *this.

	Important !
	
	In fact, the sufficient information to define the spectra resolution both in real space is N and kmax since:
	N dt = T and
	dt = 1/2kmax and
	k_i = i * 2kmax/N where i=0,1,2,...N-1
	so
	for the real space noise generations the kmin argument is not needed. It is only used for the fourier space noise power  spectrum realizations
	in which case the  spectra resolution is given by:
	dk = (kmax-kmin)/N.
	
	For the real space noise realizations the resolution is defined by: dk = 2kmax/N in fourier space and the time resolution is defined by:
	dt = 1/2kmax;
	
	This change was introduced on: Jul 4, 2011, 2:44:10 PM. It may be needed to correct some of the test programs (mostly related to MCMC tests)
	to account for this change.


	\date Aug 6, 2010, 12:50:08 PM
	\author Bartosz Lew
  */
  mscsFunction& mkPowerLawNoise(double A, double alpha, double k0, double kmin, double kmax, long N, long seed, mscsFunction* powerSpectrum=NULL, cpedsRNG* rng=NULL, bool realSpace=true);
//  mscsFunction& mkPowerLawNoise(double A, double alpha, double k0, double *klist, long N, long seed, mscsFunction* powerSpectrum=NULL, cpedsRNG* rng=NULL, bool realSpace=true);


  /*!
	\brief 
	\details 
	@param
	@return

	\date Nov 13, 2011, 10:48:12 PM
	\author Bartosz Lew
  */
  mscsFunction& mkPhaseNoise(long N, long seed, cpedsRNG* rng=NULL);
  
  
  /*!
	\brief generate power law gaussian noise in Fourier space on the allocated arguments
	\details 
	@param A - amplitude of noise
	@param alpha - spectral index of the power law noise: 0 - white, -1 - 1/f red noise, -2 - 1/f^2 noise, 1 - f noise (blue)etc.
	@param k0 - defines pivot scale
	@param seed - seed for the RN sequence generation (only used if rng is not provided)
	@param powerSpectrum - reference to the object at which to append the power spectrum of the generated real signal. If the realSpace=false, the real space not calculated at all, and if powerSpectrum==NULL then the 1/f noise power spectrum is appended to *this object
	@param realSpace - defines if the noise should be defined in real space of Fourier space. If true then Fourier transform is done on the generated data and returned. - THIS IS NOT IMPLEMENTED YET
	@return

	Currently this function does not perform the Fourier transform to the real space.
	It can only generate a given power law noise spectrum random realization for a range of arguments.
	
	
	
	\date Mar 9, 2011, 9:17:38 AM
	\author Bartosz Lew
  */
  mscsFunction& mkPowerLawNoise(double A, double alpha, double k0, long seed, mscsFunction* powerSpectrum=NULL, cpedsRNG* rng=NULL, bool realSpace=true);

  /*!
	\brief generates a histogram from the input data
	\details 
	@param data - list of data
	@param Nbins - number of bins in the histogram
	@param binAlign - bin locations (-1 - x vals taken as beginning of the bin, 0 - center, 1 - end of the bin)
	@return returns *this

	\date Apr 5, 2013, 12:53:13 PM
	\author Bartosz Lew
   */
  mscsFunction& mkHistogram(cpedsList<double>& data, long Nbins=100, long binAlign=0);

  mscsFunction& mkHistogram(cpedsList<double>& data, double fromX, double toX, long Nbins, long binAlign=0);
  
  /*!
	\brief convert arguments of this function from unitx time to JD
	\details 
	@param
	@return returns this

	\date Jan 9, 2013, 11:37:21 PM
	\author Bartosz Lew
   */
  mscsFunction& convertUnixTimeToJD_x();


  
  /********************/
  /* FUNCTIONS BASKET */
  /********************/
  static double FN_gaussian1D(double x, double A, double m, double s) { return A*exp(-(x-m)*(x-m)/(2.0*s*s)); }
  static double FN_gaussian1D_A1m0(double x, double s) { return exp(-(x*x)/(2.0*s*s)); }
  static double FN_multiquadratic1D(double r, double r0) { return sqrt(r*r+r0*r0); }
  

  
  
  /*!
    \brief concatenates two functions together.
    \details
    @param f - function to be added to this function
    @return returns a const reference to this function object

    The new points are added at the end of the function.
    NOTE: No sorting is done.

    \date 2010/01/27 21:51:16
    \author Bartosz Lew
  */
  const mscsFunction& concatenate(const mscsFunction& f);

  /*!
    \brief inserts a function f on k'th position.
    \details
    @param f - function to be inserted in this function
    @param k - insert function at position k pushing forward existing points starting from k'th point. 
    -1 - default indicates the last point and the function f will be pasted before the the last element.
    If you want to paste after the last element then use concatenate or specify the correct positive k
    @param overwrite - whether or not to overwrite the existing data. If the function reaches the end of range a new point is 
    added at the end.
    @return returns this 

    NOTE: No sorting is done.

    \date 2010/11/08 
    \author Bartosz Lew
  */
  mscsFunction& paste(const mscsFunction& f, long k=-1, bool overwrite=false);

  /*!
	\brief cuts out N points from the function starting from k'th element 
	\details 
	@param N - number of points to cut
	@param k - starting index
	@return the cut points are returned
  
	if k < 0 then last N points will be cut out

	\date Nov 8, 2010, 4:16:23 PM
	\author Bartosz Lew
  */
  mscsFunction cut(long N,long k=-1);
  mscsFunction cut(double from, double to);
  mscsFunction copy(long N,long k=-1);
  mscsFunction copy(double from, double to);
  
  /*!
	\brief removes the points that are identical
	\details
	@return returns *this

	  This method has O(N^2/2) complexity so it is very slow. If you need faster unique then see uniqueFast
	\date May 20, 2010, 3:03:41 PM
  */
  const mscsFunction& unique();

  /*!
	\brief sorts the function in argument ascending order and removes neighboring points that are separated in their arguments by less then acc
	\details
	@param acc - maximal separation of points below which points are considered the same.
	@return returns *this

	The point with larger index will be removed.

	\date May 20, 2010, 3:05:47 PM
  */
  const mscsFunction& uniqueXFast(double acc=0);
  const mscsFunction& uniqueYFast(double acc=0);

  /*!
	\brief  averages the function values within the same argument values
	\details 
	@param acc - accuracy below which arguments separations are considered as zero
	@return returns this

	This function assumes that the function is sorted arg ascending
	
	\date Oct 1, 2012, 3:20:42 PM
	\author Bartosz Lew
   */
  mscsFunction& average_sameArgs(double acc=0, mscsFunction* stdev=NULL);
  
  //! inverts the function and returns the modifiable reference to itself
  mscsFunction& invert();

  
  //! convert from phiR convension of X and Y to Re Im convention of X and Y
  mscsFunction& convertPhiRToReIm();
  mscsFunction& convertReImToPhiR();
  
  /*!
	\brief extracts spikes from the function
	\details 
	@param n - size of the window in numbers
	@param m - step by which to move the window
	@param nsigma - number of standard deviations above which threshold the point should be removed
	@param startFrom - value of the argument to start with
	@param direction - the direction of shifting the window: 12 - the windown is shifted towards increasing point indexes; 21 - otherwise
	
	The function extracts the points from the function that outstand by more than nsigma of standard deviations calculated within the
	window of n numbers that is swept across the range where the function is defined by m numbers. This makes sense when the function
	is defined on equi-spaced grid.
	For the last shift the window is set at full width starting at the end of the point indexes and reaching as far as to cover all 
	n points unless the remaining number of points is less the n, in which case n will be set for the remaining number of points.
	
	The spikes found will be removed from this function.
	@return this method returns the function containing the extracted spike points.

	\date Sep 8, 2010, 2:04:53 PM
	\author Bartosz Lew
  */
  mscsFunction extractSpikes(long n, long m, long nsigma, double startFrom=0,long direction=12, bool startFromAsIdx=true);

  /*!
	\brief as above but it can search for spikes in a selected range of data
	\details 

	\date Feb 7, 2011, 3:17:16 PM
	\author Bartosz Lew
  */
  mscsFunction extractSpikesRange(long n, long m, long nsigma, double startFrom=0,double endAt=0, long direction=12, bool startFromAsIdx=true);
  
  /*!
	\brief get mean function value
	\details 
	@return mean value of the function

	\date Jul 23, 2010, 12:04:20 AM
	\author Bartosz Lew
  */
  double meanf() const;
  double meanX(bool weightByY=false) const;

  //! calculate standard deviation of the values of the function
  double stdev();
  //! calculate variance deviation of the values of the function
  double variance();
  //! calculate the covariance with f2
  double covariance(mscsFunction& f2);
  //! calculate rms of the values of the function
  double rms();


  /*!
    \brief calculates the derivative of the function
    \details
    @param inPlace - if true then no extra memory is allocated for derivation and the derivative is stored in the object - THIS IS NOT IMPLEMENTED YET\n
                     if flase then the original function is not changed and derivative is returned. Will allocate memory for another function
    @param periodX - if NULL then non-periodic derivative is calculated and if anything else, then periodic derivative is calculated. No particular value 
    is returned.
    @param periodY - value that specifies the pediod and implicitly the range of allowec y-values of the function [0,pediodY). 
    If function detects that the function values had been wrapped around period
    then the derivative will choose the smaller of the two possible values of the derivative. Special use only.
    @return calculated derivative of the initial function

    \date 2010/01/19 10:09:00
    \author Bartosz Lew
  */
  const mscsFunction derivative(bool inPlace=false, double* periodX=NULL, double* periodY=NULL);

  /*!
    \brief performs interpolation of function
    \details
    @param Xint - pointer to array of x values on which to interpolate
    @param Nint - size of the Xint array
    @param type - type of interpolation (cspline, linear, auto, akima_periodic, cspline_periodic)\n
    if periodic type is chosen then the value of the first point must be the same as the value of the last point
    @param inPlace - if true the interpolation will be inserted in the original function and returned 
    otherwise only the interpolated function is returned. If inplace=false then the original function is not modified.
    @return returns the interpolated function

    \date 2010/01/27 15:12:56
    \author Bartosz Lew
  */
  const mscsFunction interpolate(double *Xint, long Nint, string type, bool inPlace=false);

  /*!
    \brief performs interpolation of function only within indicated range of arguments
    \details
    @param from - abscissa value to start interpolation from
    @param to - abscissa value to end interpolation at
    @param dx - increment used to derive arguments for interpolation
    @param inPlace - if true the interpolation will be interted in the original function, otherwise it will only be retutne in a separate function
    @param exact - defines where the grid of interpolated values should start. By default it starts with the first valid point inside of the range in which
    the function is defined (see below)
    @param acc - range defining maximal separation from the interpolation limits for the points to be snapped to those limits;
    @return - sorted, interpolated frunction

    This method assumes that the function is sorted in ascending order in arguments.\n
    The interpolation - by default (exact=false) - will be done for argument values starting from the argument closest to the "from" value and with step dx
    until the argument closest to the "to" but smaller than "to".\n
    Another option is when "exact=true". Then,
    the interpolation will be done for argument values starting from "from+N *dx" with intenger N yielding condition such that the argument is the first value
    that is inside of the range upon which the function is defined ( [from,to) ),
    and will proceed with step dx until the argument from+M*dx with integer M for which the argument is still inside of the range inside which the function is
    defined.\n
    If acc argument is provided and different from 0, then if the value at which the interpolation will start is closer to "from" than acc then it will be treated
		as equal to from. Similarly, if acc argument is provided and different from 0 then if the value at which the interpolation will end lies closer
		to "to" than acc then it will be set exactly to "to" and the interpolation range will be [from,to].\n
		This is different when acc=0 - then the last interpolated value will necessarily be inside of
		[from, to) range.


    \date 2010/01/27 16:34:36
    \author Bartosz Lew
  */
  const mscsFunction interpolate(double from, double to, double dx, bool inPlace=false, string type="linear", bool exact=false, double acc=0);
  const mscsFunction interpolate(double from, double to, long N, bool inPlace=false, string type="linear");

  //  const mscsFunction interpolate(double *onX, long Nx, bool inPlace=false, string type="linear", bool exact=false, double acc=0);
  /*!
    \brief performs interpolation of function within the whole range of arguments
    \details
    @param dx - separation of the interpolated function arguments
    @param inPlace - if true the interpolation will be inserted in the original function, otherwise it will only be returned in a separate function
    @return returns the interpolated function

	  The inetrpolation will be done for argument values starting from the argument of the first point of the function + dx upto the argument of the last point
	  of the function - dx

    \date 2010/01/27 15:12:56
    \author Bartosz Lew
  */
  const mscsFunction interpolate(double dx, bool inPlace=false, string type="linear");
  const mscsFunction interpolate(double dx, string type, double DX, string type2);


  /*!
    \brief performs extrapolation of function values outside of the allocated region
    \details
    @param inPlace - the extrapolated points will be added to this object, otherwise this object will not be changed
    @return
    This function extrapolates function values from from to minimal
    argument of the function and from to to maximal argument of the function with step dx.
    So from/to should be smaller/larger than minimal/maximal argument of the function in order
    for this function to do anything.

    The extrapolation is linear at the moment using the two first and last values.

    \date 2010/03/10 00:33:00
  */
  const mscsFunction extrapolate(double from, double to, double dx, bool inPlace=false, string type="linear");

  /*!
	\brief find estimates of roots of the tabulated function
	\details 
	@param periodic - if >0 then the root is calculated on periodic boundary conditions assuming the specified period
	@return vector of roots

	\date Nov 12, 2017, 8:03:04 PM
   */
  vector<double> findRoot(double period=0);

  /*!
    \brief integrates the function the returns the result
    \details
    @return returns the integral

    The function needs to be sorted ascending in arguments
  */
  double integrate();

  double integrate(double xmin,double xmax);

  /*!
	\brief convolve the current function with window function w
	\details 
	@param w - window function (defined in real space)
	@param wreal - specifies if the w is defined in real or fourier space; (default- real sapce; fft is done on w before convolving in fourier space)
	@return widowed real space function 
	
	The convolution is done in fourier space by multiplying the function fft with the 
	window function interpolated on appropriate frequencies (k-modes).

	\note: This method was tested to return the correct values under the following tests:
	convolution of a gaussian function with itself. The result of cyclic convolution gave similar results as compared with
	mathematica. 
	The FWHM of the convolved function grew by a factor of sqrt(2) as expected.
  
	\date May 17, 2011, 3:38:59 PM
	\author Bartosz Lew
  */
  mscsFunction convolve_window(mscsFunction w, bool wreal=true);
  
  /*!
    \brief piece-wise average the function then return the result
    \details
    @param imin - bin from this index
    @param binSize - list containing the number of points to average into one bin starting from imin;
    Each element indicates the number of points that should be averaged into a separate bin
    @param w - weights to be used for each point; Each point indicates the weight which is to be applied when doing the binning
    The size of the w list should be the total number of points that will be binned. If the size is zero, then equal weights
    will be applied.
    @return returns the binned function

    Binning is done only from imin to imax, outside of this range the function is copied without any changes;
    Binning is done according to the bintab table of size bintabs that holds a set of long type numbers defining how many ls to bin using weights w=1/bintab[i]
    or weights from the err field (previously the 3rd col) of the C_l structure for each l in which case w should be 0: w=0; otherwise it should be w=1
    the sum of values stored in bintab array should be lmax-lmin+1
    value in array bintabs - 1 corresponds to no binning; 2 - means that two multipoles are binned and so on.

    \date 2009/06/04 23:31:24 \n
    reimplemented on:\n
    \date 2010/02/19 15:59:49
    \author Bartosz Lew
  */
  mscsFunction& binFunction(long imin, cpedsList<long>& binSize, cpedsList<double>& w);

  /*!
	\brief bin function at fixed argument intervals
	\details 
	@param dx - interval defining a fixed bin width
	@param binSize - output parameter. Contains number of elements used for each bin. Data is being appended to this object.
	Object is not cleaned.
	@param methodX - binning method for arguments
	@param methodY - binning method for values

	Can be any of:
	
	bin_center,bin_max,bin_min - only for methodX
	mean - \n
	median- \n
	min - \n
	max - \n

	@return *this
	
	The function should be sorted before calling this method.
	The i'th bin center is defined as:
	
	Xmin+i*dx+dx/2, except for the last bin which may be smaller.

	\date May 25, 2017, 10:11:10 AM
   */
  mscsFunction binFunction(double dx, cpedsList<long>& binSize, string methodX="mean", string methodY="mean");

  /*!
	\brief bin function using space few * O(N)
	\details 
	  Parameters: same as above
	@return same as above

	\date Aug 3, 2010, 1:50:26 PM
	\author Bartosz Lew
  */
  mscsFunction& binFunctionLin(long imin, cpedsList<long>& binSize, cpedsList<double>& w);

  /*!
	\brief bin function using space few * O(N)
	\details 
	@return

	  Input parameters are as above but the implementation doesn't use matrices but linear arrays so it is much faster.
	  
	  NOTE: the binning matrix (or what is left of it) can be reduced in size down to the largest bin size in the binSize list.
	  Currently the space is allocated for the full length of the input data. This not necessary in case of the top-hat binning
	  problem.

	\date Oct 28, 2010, 9:15:37 AM
	\author Bartosz Lew
  */
  mscsFunction& binFunctionLin2(long imin, cpedsList<long>& binSize, cpedsList<double>& w, cpedsList<double>* varInbin=NULL, cpedsList<double>* alpha=NULL, vector< cpedsList<double> >* stats=NULL);
  
  /*!
	\brief same as binFunctionLin2 but the binSizes list is generated using the bs parameter
	\details 
	@param imin - the starting bin
	@param bs - width of the bins. This can be a real number in which case the bin widths are derived as
	b_1 = round(bs)
	b_i = round(i*bs) - b_{i-1}
	@param binSize - a reference to the list of the generated bins. This need not be generated as this list is cleared inside of this method and filled with the
	generated bins that are used in binning.
	@param w - a reference to the list of the weights to be used. This may by empty in which case same unit weights will be used.
	@return returns *this

	\date Jul 4, 2011, 12:45:09 PM
	\author Bartosz Lew
   */
  mscsFunction& binFunctionLin2(long imin, double bs, cpedsList<long>& binSize, cpedsList<double>& w);

  
  /*!
	\brief bin data in bins of increasing width in a geometrical sequence
	\details 
	@param imin - start binning ar this index
	@param firstBin - width of the first bin in number of points
	@param gm - geometrical sequence multiplier
	@param binSize - reference to an empty list that will contain the derived bin sizes (if not empty - then it will be cleaned)
	@param w - weights to apply for individual data points (must be of the length of the data)
	@return binned function

	\date Jan 31, 2011, 11:47:43 AM
	\author Bartosz Lew
  */
  mscsFunction& binFunctionLin2geo(long imin, double firstBin, double gm, cpedsList<long>& binSize, cpedsList<double>& w, cpedsList<double>* varInbin=NULL);

  //  mscsFunction& fft_real(bool fwd);
  
  /*!
	\brief calculate forward fast Fourier transform of the function and return the power spectrum
	\details 
	@param re - pointer to a function to store the real part of the Fourier transform
	@param im - pointer to a function to store the imaginary part of the Fourier transform
	@param regridX - defines how the function should be interpolated before proceeding to the fftw.
	For the default value 0 the function will not be interpolated at all. It will be assumed that the 
	function is defined on equally spaced grid and the separation of the points will be calculated as the mean separation of the points i.e.:\n
	dt=(x_last-x_first)/Npoints. 
	If eg. regridX=0.1 then the function will be probed and interpolated
	at resolution 0.1*medianStep, where medianStep is the median separation of the two neighboring points.
	
	@return power spectrum of the function
	
	If re or im can be set to NULL in which case it will not be returned.
	The resulting power spectrum wave numbers correspond to 1/(X unit), so eg for the time ordered data in seconds
	the result is in Hz and it's inverse correspond to the period for a give wave number.
	
	At the moment the re and im functions store all the non-redundant information up to the nyquist frequency =
	M=N/2+1 where division is rounded down. The i'th point of the function corresponds to the i/T frequency
	where T is the sampling period.
	
	If the data were provided on an non-uniform grid, option regridX=true makes the function to regrid the data
	onto the uniform mesh using linear interpolation with steps defined by the median step in the input function.
	
	The FFT convention is as follows: if x_i is the real signal and f_j fourier expansion coefficients then
	f_j = 1/N sum_{i=0}^{N-1} x_i exp(-2pi sqrt(-1) i j/N)

	and the power spectrum is:

	P_j = f_j^\star * f_j

	where only non-negative frequencies are considered
	
	The inverse fft would be
	
	x_i =  sum_{j=0}^{N-1} x_i exp(2pi sqrt(-1) i j/N)
	
	so, the FFT definition is with and (a,b) = (-1,-1) where the general forward DFT is defined as 
	y_j = 1/(N^[(1-a)/2]) sum_i (x_i exp(2pi b sqrt(-1) ij/N)
	
	This is a bit unusual convention but it is useful because the zero'th mode carries the value of the mean real signal and 
	the sum of the spectra is the RMS^2 of the real signal.
	
	For this definition of the FT the parseval's theorem is:
	
	1/N sum_i (x_i)^2 = sum_j |y_j|^2 
	
	
	With this definition the unit of the power spectra is V^2 if the input signal is a time series of voltage measurements 
	
	If the input signal is a time series in time units then N -> T=N*dt where T is the total sampling period. -- not sure this is valid after the change of FT definition


	\date Jul 22, 2010, 7:16:54 PM
	updated on: Apr 25, 2012, 3:18:25 PM
	\author Bartosz Lew
  */
  mscsFunction powerSpectrum(mscsFunction* re=NULL, mscsFunction* im=NULL, double regridX=0, bool calibrateByPointsCount=true);
  
//  mscsFunction 
  
  /*!
	\brief performs inverse fourier transform using re and im information
	\details 
	@param re - real part of the signal stored on y values
	@param im - imaginary part of the signal stored on y values
	@return
	The re and im are assumed to be non-negative frequencies only. They must be of the same length.

	\date Aug 22, 2011, 9:59:06 PM
	\author Bartosz Lew
   */
  mscsFunction inverseFFT(mscsFunction& re, mscsFunction& im, double* x=NULL, bool deleteX=true);

  /*!
	\brief calculates the slow forward Fourier transform on this function inside defined range of k values.
	\details 
	@param
	@return returns the power spectrum
	THe fourier transform is defined as:
	
	F(k) = int_x exp(-i 2pi k x) f(x) dx
	integrated over the x range where the function is defined.
	

	\date Apr 24, 2012, 5:09:48 PM
	\author Bartosz Lew
   */
  mscsFunction powerSpectrum(double kmin, double kmax, double dk,mscsFunction* Re=NULL, mscsFunction* Im=NULL);
  
  /*!
	\brief calculates Fourier series on a function and optionally returns coefficients
	\details 
	@param N - series order
	@param A - pointer to amplitude coefficients, arguments store period
	@param phi - pointer to phase coefficients, arguments store period
	
	The function is reconstructed by the sum:
	
	f(x) = mean + Sum_i A_i sin(2pi/T_i + phi_i)
	
	where mean is mean value of the function.
	
	@return function that is assembled from the calculated coefficients up to order N

	\date Oct 28, 2016, 11:00:27 PM
	\author Bartosz Lew
*/
  mscsFunction FourierSeries(long N, double regridX=0, mscsFunction* A=NULL, mscsFunction* phi=NULL);  
  
  /*!
	\brief derive the expectation value of the function as \int x f(x) dx
	\details 
	@return expectation value

	\date Nov 23, 2010, 12:00:30 PM
	\author Bartosz Lew
  */
  double expectationValue();
  
  
  /*!
	\brief convert this function into a cumulative distribution function
	\details 
	@param zoomIn - optional parameter. If false this just calculates CDF from the function.
	If true, then extra checks are done for nonsense values in CDF like multiple-ones in case
	of integrating probability values that are too small to hold on double precision.
	If the PDF drops very rapidly then CDF quickly becomes one and it is numerically just one for the rest of the domain
	it is defined upon.
	So if the zoomIn option if true, then checks are done if CDF becomes one and no further integration is done.
    Instead the CDF is interpolated on the values that are <=1 to match the initial number of points. (This is effectively increasing resolution
    in the range where interesting things happen and neglecting other regions-- hence zoomIn)
	If true then the range of arguments will change.
	
	@return returns *this

	\date Mar 14, 2011, 4:55:17 PM
	\author Bartosz Lew
  */
  mscsFunction& mkCDF(bool zoomIn=false, double acc=1e-10);
//  mscsFunction& mkcCDF(bool zoomIn=false, double acc=1e-10);
  
  mscsFunction integrateX_Xmax();

  /*!
	\brief make the function normalized
	\details 
	This calculates the integral of the function and divides its values by the integral.

	\date Dec 13, 2010, 8:57:23 PM
	\author Bartosz Lew
  */
  void normalize();
  
  /*!
	\brief Calibrate the signal by the piece-wise calculated standard deviation in the real space. 
	\details 
	@param windowSize - number of points to calculate standard deviation with in one window.
	@param windowStep - the number of points by which the window is to be moved. By default this is -1 which means that the 
	window will be shifted by the windowSize.
	@return
	The width of the window is defined by this parameter (It must be larger than 3 for this option to be activated). 
	The calibration can be done on the moving window by setting windowStep as required.
	The signal is calibrated in the neighbouring, non-overlapping windows, by default when windowStep=-1.
	
	\date Aug 25, 2011, 10:17:30 AM
	\author Bartosz Lew
   */
  mscsFunction& calibrateByStDev(long windowSize, long windowStep=-1);
  
  
  /*!
	\brief digitize the function with a given number of levels within given range of values
	\details 
	@param n - number of digitization levels
	@param vmin - minimal value after quantization (the lower values in the function will be saturated with this value)
	@param vmax - maximal value after quantization (the higher values in the function will be saturated with this value)
	@return this method returns *this

	
	\date Jul 4, 2011, 10:32:40 AM
	\author Bartosz Lew
   */
  mscsFunction& quantize(long n, double vmin=0, double vmax=0);
  
  
  /*!
	\brief calculate correlation coefficient between this function and function f2
	\details 
	@param f2 - function to calculate the correlation function with
	@return correlation coefficient defined as
	
			c= cov(f, f2) / ( sqrt(var(f)) +sqrt(var(f2)) )

	f2 function must be of the same length. Only the values are taken into account.
	This function changes this function and the parameter function in a way that they
	are shifted to the zero mean.

	\date Oct 1, 2011, 9:51:55 AM
	\author Bartosz Lew
  */
  double correlationCoefficient(mscsFunction& f2);

  /*!
	\brief calculates the correlation coefficient function 
	\details 
	@param f2 - function to calculate the correlation coefficients with
	@param step - number of points by which to shift one function wrt other.
	@param N - number of steps to perform; default 0 - means that all steps until periodic wrap
	If step > 1 and N==0 then N becomes n/step where n is the number of points in the interpolated
	function.	
	
	By default the f2 function is shifted forward. The shift is periodic.
	IMPORTANT: This method assumes that the two functions are defined on a regular grid.
	@return correlation coefficient function 

	\date Oct 1, 2011, 10:39:58 AM
	\author Bartosz Lew
  */
  mscsFunction correlationCoefficientFunction(mscsFunction& f2, long step=1, long N=0);
  
  //! checks the ranges and updates the function ranges structure
  /*! The updated ranges are returned.   */
  void checkRanges();
 

  /*!
	\brief kill the current random number generator if it is allocated
	\details 
	\date Jul 4, 2011, 5:37:56 PM
	\author Bartosz Lew
   */	
  void killRNGs() { if (_rns!=NULL) delete _rns; _rns=NULL; }
  
  
  /*************/
  /* OPERATORS */
  /*************/
/*   double& operator () (long x, long y); */
/*   double operator () (long x, long y); */

  const mscsFunction& operator=(const double v) { setf(v); return *this; }
  const mscsFunction& operator=(const mscsFunction& rhs);
  bool operator==(const mscsFunction& rhs);
  double operator () (long i) { return f(i); }
  double operator () (double x) const { return f(x); }
  const mscsFunction& operator += (const mscsFunction& f) { add(f); return *this; }
  const mscsFunction& operator += (const double* t) { add(t); return *this; }
  const mscsFunction& operator -= (const mscsFunction& f) { subtract(f); return *this; }
  const mscsFunction& operator *= (const mscsFunction& f) { multiply(f); return *this; }
  const mscsFunction& operator /= (const mscsFunction& f) { divide(f); return *this; }
  const mscsFunction& operator /= (const cpedsList<double>& l) { divide(l); return *this; }
  const mscsFunction& operator /= (const cpedsList<long>& l) { divide(l); return *this; }

  const mscsFunction& operator += (const double v) { add(v); return *this; }
  const mscsFunction& operator -= (const double v) { subtract(v); return *this; }
  const mscsFunction& operator *= (const double v) { multiply(v); return *this; }
  const mscsFunction& operator /= (const double v) { divide(v); return *this; }

  const mscsFunction operator + (const mscsFunction& f) const { mscsFunction tmp(*this); return tmp+=f; }
  const mscsFunction operator - (const mscsFunction& f) const { mscsFunction tmp(*this); return tmp-=f; }
  const mscsFunction operator * (const mscsFunction& f) const { mscsFunction tmp(*this); return tmp*=f; }
  const mscsFunction operator / (const mscsFunction& f) const { mscsFunction tmp(*this); return tmp/=f; }

  const mscsFunction operator + (const double v) const { mscsFunction tmp(*this); return tmp+=v; }
  const mscsFunction operator - (const double v) const { mscsFunction tmp(*this); return tmp-=v; }
  const mscsFunction operator * (const double v) const { mscsFunction tmp(*this); return tmp*=v; }
  const mscsFunction operator / (const double v) const { mscsFunction tmp(*this); return tmp/=v; }


  
  
  QPointF& operator[](long i) { return _f[i]; }


/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PROTECTED MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 protected:

/* ---------------------------- */
/* PROTECTED METHODS */
/* ---------------------------- */

#ifndef NO_HDF5
		bool hdf5DatasetExists(hid_t& file, string dsetName);
		long getHDF5dsetCount(hid_t& file);
		void setHDF5_scalarStringAttribute(hid_t& dset, string attributeName, string attributeValue);
		void setHDF5_scalarDoubleAttribute(hid_t& dset, string attributeName, double value, string attributeComment="");
		string getHDF5_stringAttribute(hid_t& file, string dsetName, string attributeName, int* errCode);
		bool hdf5createGroup(hid_t& file,string linkName);
#endif		

  
  
  //! This doens't do anyting at the moment;
  void initiate_variables();
  //! returns modifiable reference to the function structure
  QList<QPointF>& f() { return _f; }



  //! returns the function ranges structure
  frange getRanges() const;

  void clearRanges();

  
  /*!
	\brief derive complex to real fft on the input data. 
	\details 
	@param inRe - array of real part of the complex input array of size size 
	@param inIm - array of imaginary part of the complex input array of size size 
	@param M - size of the input array
	@param N - size of the output array: Must be M = N/2+1 for even N and M=(N-1)/2+1 for odd N (so the division is rounded down)
	@param deleteIn - true - the inRe and inIm will be deleted inside of this function before doing fft and after copying onto fftw_complex array
	
	@return real valued output array

	dL defines the spacings of points in the real space. This is used in combination with size, to derive the k values for the output
	real function. The in array is assumed to hold only the non-negative frequencies of the signal - ie. its size is
	size=N/2+1 where the division is rounded down for odd N. The real valued function will be of size N and will span over a domain:
	dL, 2dL, ..., N*dL.
	
	
	\date Oct 6, 2010, 1:47:22 AM
	\author Bartosz Lew
  */
  double* fft_1D_c2r(double *inRe, double* inIm, long M, long N, bool deleteIn);

  /*!
	\brief DTF of an in array data of size size
	\details 
	@param in - input array
	@param size - size of in
	@param re - pointer to output array containing the real part of the output
	@param im - pointer to output array containing the imaginary part of the output
	@param Nout - size of re and im = size/2+1 where division is rounded down
	
	The output re and im arrays need not (should not) be allocated; They will be allocated within this routine.

	\date May 17, 2011, 8:59:53 PM
	\author Bartosz Lew
  */
  void fft_1D_r2c(double *in, long size, double *re, double *im);


/* ---------------------------- */
/* PROTECTED STRUCTURES */
/* ---------------------------- */
  QList<QPointF> _f;
//  cpedsList<double> _Xerr;
//  cpedsList<double> _Yerr;
  cpedsRNG *_rns;

  frange range;

/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PRIVATE MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 private:


/* ---------------------------- */
/* PRIVATE METHODS */
/* ---------------------------- */

  static int cmp_points_X(const QPointF &p1, const QPointF &p2);
  static int cmp_points_Y(const QPointF &p1, const QPointF &p2);


/* ---------------------------- */
/* PRIVATE STRUCTURES */
/* ---------------------------- */



  typedef int (*cmpPt)(const void*, const void*);
  typedef int (*cmpPt2)(const QPointF&, const QPointF&);
  string object_name;


 };
#endif
