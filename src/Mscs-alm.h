/*!
  \file implements alm - a class for storing and managing single SHT coefficient
*/

#ifndef MSCS_ALM
#define MSCS_ALM

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
#include <complex>
#include <ostream>
#include <fstream>

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */
using namespace std;

/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
  \class mscsAlm
  \brief Encapsulates narrowes alm of double type
  \details 
  
  \date 2010/02/18 12:10:38 
  \author Bartosz Lew
*/

class mscsAlm : public complex<double> {


/* ---------------------------------------------------------------------------------------------------- */
/* CLASS PUBLIC MEMBERS */
/* ---------------------------------------------------------------------------------------------------- */
 public:

  typedef complex<double> complexDouble;

  typedef struct { 
      double R; // real part
      double I; // unreal part
  } mscs_a_lm;

  
/* ------------- */
/* CLASS FRIENDS */
/* ------------- */


/* ---------------------------- */
/* CONSTRUCTORS AND DESTRUCTORS */
/* ---------------------------- */
  mscsAlm() {}
  mscsAlm(double v) : complexDouble(v) {}
  mscsAlm(complexDouble v) : complexDouble(v) {}
  mscsAlm(double r, double i) : complexDouble(r,i) {}
  ~mscsAlm() {}

/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */

  void R(double v) { operator=(complexDouble(v,imag())); }
  void I(double v) { operator=(complexDouble(real(),v)); }
  complexDouble Z() { return complexDouble(real(),imag()); }
  double R() const { return real(); }
  double I() const { return imag(); }
  mscs_a_lm get_a_lm() const; //!< returns the alm in a C form

  /* double& r() { return real(); } */

  mscsAlm& conj() { complex<double>::operator=(std::conj(complex<double>(real(),imag()))); return *this; }
  double abs() const { return std::abs(complex<double>(real(),imag()));  }

  const mscsAlm& operator=(double v) { complex<double>::operator=(complex<double>(v,0)); return *this; }
  const mscsAlm& operator=(const complex<double>& rhs) { complex<double>::operator=(rhs); return *this; }
  const mscsAlm& operator=(const mscsAlm& rhs) { complex<double>::operator=(complex<double>(rhs.real(),rhs.imag())); return *this; }
  operator double() { return real(); }

  const mscsAlm& operator*=(const double v) { complex<double>::operator*=(v); return *this; }
  mscsAlm operator*(const mscsAlm& rhs) const { return mscsAlm(complex<double>(real(),imag())*complex<double>(rhs.real(),rhs.imag())); }

//  ostream& operator << (ostream& out) { return out<<R()<<I();}
  
  void print();
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


 };


ostream& operator<<(ostream& out, mscsAlm& a); 
//istream& operator<<(istream& in, mscsAlm& a); 


#endif /* MSCS_ALM */ 

