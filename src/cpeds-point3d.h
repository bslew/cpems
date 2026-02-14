/*!
  \file implements a point class in 3d space
*/

#ifndef CPEDSPOINT3D
#define CPEDSPOINT3D

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
#include <math.h>
#include "cpeds-common.h"
#include "cpeds-math.h"
#include "cpeds-point2d.h"

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */

using namespace std;

/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class cpedsPoint3D
  \brief Encapsulates points in 3d space
  \details 
  
  \date 2009/11/05 23:44:33 
  \author Bartosz Lew
*/
class cpedsPoint3D : public cpedsPoint2D {


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
 cpedsPoint3D() : cpedsPoint2D() { _z=0; }
 cpedsPoint3D(const cpedsPoint3D &p) : cpedsPoint2D(p) { setZ(p.z()); }
 cpedsPoint3D(double x, double y, double z) : cpedsPoint2D(x,y) { setZ(z); }
 cpedsPoint3D(double x, double y) : cpedsPoint2D(x,y) {  setZ(0); }
 cpedsPoint3D(const cpedsPoint2D &p) : cpedsPoint2D(p) { setZ(0); }
 virtual  ~cpedsPoint3D() {}


/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */
  void setZ(double Z) { _z=Z; }
  double z() const { return _z; }
  double& z() { return _z; }
  double* zP() { return &_z; }
  void set(double X, double Y, double Z) { setX(X); setY(Y); _z=Z; }

  void add(const cpedsPoint3D &p) { setX(x()+p.x()); setY(y()+p.y()); setZ(z()+p.z()); }
  void subtract(const cpedsPoint3D &p) { setX(x()-p.x()); setY(y()-p.y()); setZ(z()-p.z()); }
  void multiply(const cpedsPoint3D &p) { setX(x()*p.x()); setY(y()*p.y()); setZ(z()*p.z()); }
  void divide(const cpedsPoint3D &p) { if (p.x()!=0) setX(x()/p.x()); if (p.y()!=0) setY(y()/p.y()); if (p.z()!=0) setZ(z()/p.z()); }
  cpedsPoint3D& multiply(const double &k) { setX(x()*k); setY(y()*k); setZ(z()*k); return *this; }
  cpedsPoint3D& divide(const double &k) { setX(x()/k); setY(y()/k); setZ(z()/k); return *this; }
  cpedsPoint3D& add(const double &k) { setX(x()+k); setY(y()+k); setZ(z()+k); return *this; }
  cpedsPoint3D& subtract(const double &k) { setX(x()-k); setY(y()-k); setZ(z()-k); return *this; }
  cpedsPoint3D& shiftX(double delta) { setX(x()+delta); return *this; }
  cpedsPoint3D& shiftY(double delta) { setY(x()+delta); return *this; }
  cpedsPoint3D& shiftZ(double delta) { setZ(x()+delta); return *this; }

  /*!
	\brief a dot product between two vectors defined by this and the provided point
	\details 
	@param p - point for the provided vector
	@return dot product

	\date Jan 17, 2012, 11:46:55 AM
	\author Bartosz Lew
   */
  double dot(const cpedsPoint3D &p) const { return x()*p.x()+y()*p.y()+z()*p.z(); }
  /*!
	\brief returns the angle between vectors defined by this and provided point dia their dot product.
	\details 
	@param p - point for the provided vector
	@return angle between the corresponding vectors in radians

	\date Jan 17, 2012, 11:45:51 AM
	\author Bartosz Lew
*/
  double angle(const cpedsPoint3D &p) const { return acos(dot(p)/(dist()*p.dist())); }
  

  //! return the distance between two points
  double dist(const cpedsPoint3D& p) const { 
    double dx=x()-p.x();  
    double dy=y()-p.y();
    double dz=z()-p.z();
    return sqrt(dx*dx+dy*dy+dz*dz);
  }
  //! distance from the origin of the CS
  double dist() const { 
    double dx=x();  
    double dy=y();
    double dz=z();
    return sqrt(dx*dx+dy*dy+dz*dz);
  }
  
  cpedsPoint3D& Rx(double Ax) { cpeds_Rx(x(),y(),z(),Ax,xP(),yP(),zP()); return *this; }
  cpedsPoint3D& Ry(double Ay) { cpeds_Ry(x(),y(),z(),Ay,xP(),yP(),zP()); return *this; }
  cpedsPoint3D& Rz(double Az) { cpeds_Rz(x(),y(),z(),Az,xP(),yP(),zP()); return *this; }
  
  //! true is the point p lies within radius r from *this point.
  bool withinRadius(double r, const cpedsPoint3D &p) {
    if (dist(p)<=r) return true;
    return false;
  }

  virtual cpedsPoint3D roundToInteger() { setX(round(x())); setY(round(y())); setZ(round(z())); return *this; }
  virtual void print_point(string name="") const { printf("%s> x: %lf y: %lf z: %lf\n", name.c_str(), x(),y(),z());  }

  virtual bool operator==(const cpedsPoint3D &p) const { return (x()==p.x() && y()==p.y() && z()==p.z()); }
  virtual const cpedsPoint3D operator=(const cpedsPoint3D &p) { if (&p!=this) { setX(p.x());  setY(p.y()); setZ(p.z()); } return *this; }
  virtual bool operator!=(const cpedsPoint3D &p) const { return !(*this==p); }
  virtual cpedsPoint3D operator+=(const cpedsPoint3D &p) { add(p); return *this; }
  virtual cpedsPoint3D operator-=(const cpedsPoint3D &p) { subtract(p); return *this; }
  virtual cpedsPoint3D operator*=(const cpedsPoint3D &p) { multiply(p); return *this; }
  virtual cpedsPoint3D operator*=(const double k) { multiply(k); return *this; }
  virtual cpedsPoint3D operator/=(const double k) { divide(k); return *this; }
  virtual cpedsPoint3D operator+=(const double k) { add(k); return *this; }
  virtual cpedsPoint3D operator-=(const double k) { subtract(k); return *this; }
  virtual cpedsPoint3D operator/=(const cpedsPoint3D &p) { divide(p); return *this; }
  virtual cpedsPoint3D operator+(const cpedsPoint3D &p) const  { return cpedsPoint3D(*this)+=p; }
  virtual cpedsPoint3D operator-(const cpedsPoint3D &p) const  { return cpedsPoint3D(*this)-=p; }
  virtual cpedsPoint3D operator*(const cpedsPoint3D &p) const  { return cpedsPoint3D(*this)*=p; }
//  virtual cpedsPoint3D operator*(const double k) const { return cpedsPoint3D(*this)*=p; }
  virtual cpedsPoint3D operator/(const cpedsPoint3D &p) const  { return cpedsPoint3D(*this)/=p; }
  
  friend cpedsPoint3D operator * ( double k, cpedsPoint3D& v ) { return cpedsPoint3D(v)*=k; }
  friend cpedsPoint3D operator / ( double k, cpedsPoint3D& v ) { return cpedsPoint3D(v)/=k; }
  friend cpedsPoint3D operator + ( double k, cpedsPoint3D& v ) { return cpedsPoint3D(v)+=k; }
  friend cpedsPoint3D operator - ( double k, cpedsPoint3D& v ) { return cpedsPoint3D(v)-=k; }
  
  cpedsPoint3D operator - (cpedsPoint3D& v) const { cpedsPoint3D r(*this); r-=v; return r; }
//  cpedsPoint3D operator - (cpedsPoint3D& v) const { return v.multiply(double(-1.0)); }
  cpedsPoint3D operator - () { return multiply(double(-1.0)); }
//  friend cpedsPoint3D operator - (cpedsPoint3D& v) const { return v.multiply(double(-1.0)); }

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
  double _z;

 };

ostream& operator<<(ostream& s, cpedsPoint3D& p);

#endif /* CPEDSPOINT3D */ 

