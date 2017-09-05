/*!
  \file Defines the correlation function class
*/

#ifndef MSCS_CORRELATION_FUNCTION
#define MSCS_CORRELATION_FUNCTION
/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */
#include "Mscs-function.h"


/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class mscsCorrelationFunction
  \brief Encapsulates a correlation function
  \details 
  
  \date 2009/05/28 15:45:22 
  \author Bartosz Lew
*/
class mscsCorrelationFunction : public mscsFunction {


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
 mscsCorrelationFunction() : mscsFunction("Cth") {}
 mscsCorrelationFunction(string name) : mscsFunction(name) {}
 ~mscsCorrelationFunction() {}

/* ---------------------------- */
/* PUBLIC METHODS */
/* ---------------------------- */
  /* void savetxtC_th(strarg C_th_file); // save the correlation function */
  /* void loadtxtC_th(strarg C_th_file); */
  /* void savebinC_th(strarg C_th_file); // save the correlation function */
  /* void loadbinC_th(strarg C_th_file); */


/*   int plot_C_l(int plottype);   // plots the power spectrum of the map */
/*   int plot_C_th(int plottype);   // plots the power spectrum of the map */

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

#endif
