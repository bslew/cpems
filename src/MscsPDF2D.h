/*!
 \file MscsPDF2D.h - 
 */

#ifndef INCLUDE_MSCSPDF2D_H_
#define INCLUDE_MSCSPDF2D_H_

#include "Mscs-function3dregc.h"
#include "cpeds-point2d.h"

/*!
 \class MscsPDF2D
 \brief Encapsulates 2D PDF.
 \details 
 
 \date created: May 31, 2017, 10:19:56 AM 
 \author Bartosz Lew
 */
class MscsPDF2D : public mscsFunction3dregc {
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PUBLIC MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	public:
		
		/* ------------- */
		/* CLASS TYPES */
		/* ------------- */
		typedef struct {
			double contour_density;
			double normalization;
		} MscsPDF2D_t;

		/* ---------------------------- */
		/* CONSTRUCTORS AND DESTRUCTORS */
		/* ---------------------------- */
		MscsPDF2D();
		~MscsPDF2D();
		MscsPDF2D(const mscsFunction3dregc& parent);
		
		/* ---------------------------- */
		/* PUBLIC METHODS */
		/* ---------------------------- */
		
		/*!
			\brief calculate CR contour for a given CL
			\details 
			@param CL - conficence level (e.g. 0.68)
			@return returns set of points approximating vertices of a polygon corresponding to requested CR 
		
			\date May 31, 2017, 11:09:06 AM
		*/
		mscsFunction getContour(double CL, double* LVL=NULL);

		mscsFunction getCDF();
		double getIntegral(double zMin);
		
		double& f(long i, long j) { return fRe(i,j,0); }
		
		/*!
			\brief make the PDF normalized.
			\details 
		
			\date May 31, 2017, 11:44:08 AM
		*/
		void normalize();

		
		
		MscsPDF2D& operator=(const mscsFunction3dregc& rhs);

		
		
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
		MscsPDF2D_t _MscsPDF2D_data;
		
};
#endif /* INCLUDE_MSCSPDF2D_H_ */ 

