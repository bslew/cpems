/*!
 \file MscsPDF1D.h - 
 */

#ifndef INCLUDE_MSCSPDF1D_H_
#define INCLUDE_MSCSPDF1D_H_

#include "Mscs-function.h"


/*!
 \class MscsPDF1D
 \brief Encapsulates 
 \details 
 
 \date created: Jun 6, 2017, 5:39:07 PM 
 \author Bartosz Lew
 */
class MscsPDF1D :public mscsFunction {
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PUBLIC MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	public:
		
		/* ------------- */
		/* CLASS TYPES */
		/* ------------- */
		typedef struct {
			mscsFunction CDF_sel;
			double CDF_last;
			double twoSidedIntegral_lastyMin;
		} MscsPDF1D_t;

		/* ---------------------------- */
		/* CONSTRUCTORS AND DESTRUCTORS */
		/* ---------------------------- */
		MscsPDF1D();
		~MscsPDF1D();

		MscsPDF1D(const mscsFunction& parent);
		/* ---------------------------- */
		/* PUBLIC METHODS */
		/* ---------------------------- */

		/*!
			\brief calculate and return confidence region for given confidence level
			\details 
			@param CL - confidence level [0,1). E.g 0.68 for 1-sigma confidence region.
			@param LVL - PDF value that corresponds to the requested CL - or more specifically
			the value of the inverse CDF (calculated sideways from the mode of the PDF) which 
			corresponds (the value) to the requested CL.
			@return returns a list of abscissa values that yield PDF==LVL equation.
			There should be two values for a well behaved unimodal distributions.
		
			\date Jun 6, 2017, 5:44:11 PM
		*/
		cpedsList<double> getCR(double CL, double* LVL=NULL);
		
		
		MscsPDF1D& operator=(const mscsFunction& rhs);

		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PROTECTED MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	protected:
		
		/* ---------------------------- */
		/* PROTECTED METHODS */
		/* ---------------------------- */
		
		/*!
			\brief integrate sideways from the PDF maximal within range defined by yMin value
			\details 
			@param
			@return
		
			\date Jun 6, 2017, 6:10:36 PM
		*/
		double getTwosidedIntegral(double yMin);
		mscsFunction getCDF();

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
		MscsPDF1D_t _MscsPDF1D_data;
		
};
#endif /* INCLUDE_MSCSPDF1D_H_ */ 

