/*!
  \file implements a point class in 2d space
*/

#ifndef CPEDSPOINT2D
#define CPEDSPOINT2D

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
#include "cpeds-common.h"

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */


/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class cpedsPoint2D
  \brief Encapsulates points in 2d space
  \details 
  
  \date 2009/11/05 23:44:33 
  \author Bartosz Lew
*/
class cpedsPoint2D {


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
  cpedsPoint2D() { setX(0); setY(0); }
  cpedsPoint2D(const cpedsPoint2D &p) { setX(p.x());  setY(p.y()); }
  cpedsPoint2D(double x, double y=0) {  setX(x); setY(y); }
  virtual ~cpedsPoint2D() {}


/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */
  void setX(double X) { _p.x=X; }
  void setY(double Y) { _p.y=Y; }  
  void setXY(double X, double Y) { setX(X); setY(Y); }
  double x() const { return _p.x; }
  double y() const { return _p.y; }
  double getX() { return _p.x; }
  double getY() { return _p.y; }
  double& x() { return _p.x; }
  double& y() { return _p.y; }
  double* xP() { return &_p.x; }
  double* yP() { return &_p.y; }

  void add(const cpedsPoint2D &p) { setX(x()+p.x()); setY(y()+p.y()); }
  void subtract(const cpedsPoint2D &p) { setX(x()-p.x()); setY(y()-p.y()); }

  //! return the distance between two points
  double dist(const cpedsPoint2D& p) { 
    double dx=x()-p.x();  
    double dy=y()-p.y();
    return sqrt(dx*dx+dy*dy);
  }
  
  //! true is the point p lies within radius r from *this point.
  bool withinRadius(double r, const cpedsPoint2D &p) {
    if (dist(p)<=r) return true;
    return false;
  }

  virtual void print_point() const { printf("x: %lf h y: %lf\n", x(),y());  }

  virtual bool operator==(const cpedsPoint2D &p) const { const bool r=(x()==p.x() and y()==p.y()); return r; }
  virtual bool operator!=(const cpedsPoint2D &p) const { return !(*this==p); }
  virtual cpedsPoint2D operator=(const cpedsPoint2D &p) { if (&p!=this) { setX(p.x());  setY(p.y()); } return *this; }
  virtual cpedsPoint2D operator+=(const cpedsPoint2D &p) { add(p); return *this; }
  virtual cpedsPoint2D operator-=(const cpedsPoint2D &p) { subtract(p); return *this; }
  virtual const cpedsPoint2D operator+(const cpedsPoint2D &p) const  { cpedsPoint2D o(*this); o+=p; return o; }
  virtual const cpedsPoint2D operator-(const cpedsPoint2D &p) const  { cpedsPoint2D o(*this); o-=p; return o; }
//  virtual cpedsPoint2D& operator<<(const cpedsPoint2D& vec) const { cout << "x: " << vec.x() << " y: " << vec.y() ;  return *this; }
//  virtual void operator<<(const cpedsPoint2D& vec) const { cout << "x: " << vec.x() << " y: " << vec.y() ;   }
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
  typedef struct { 
    double x; 
    double y; 
  } cpeds_point2D;

  cpeds_point2D _p;

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

//void cout::operator<<(const cpedsPoint2D& vec) const { cout << "x: " << vec.x() << " y: " << vec.y() ;   }

#endif /* CPEDSPOINT2D */ 

