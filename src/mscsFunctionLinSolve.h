/*!
 \file mscsFunctionLinSolve.h - 
 */

#ifndef MSCSFUNCTIONLINSOLVE_H_
#define MSCSFUNCTIONLINSOLVE_H_

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */

/* STANDALONE HEADERS */
#include "Mscs-function3dregc.h"
/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */

/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
 \class mscsFunctionLinSolve
 \brief Encapsulates linear algebra
 \details 
 
 \date created: Jan 20, 2015, 11:47:34 PM 
 \author Bartosz Lew
 */
class mscsFunctionLinSolve : public mscsFunction3dregc {
		
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PUBLIC MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	public:
		
		
		/* ------------- */
		/* CLASS TYPES */
		/* ------------- */
		typedef struct {
				
		} mscsFunctionLinSolve_t;

		/* ---------------------------- */
		/* CONSTRUCTORS AND DESTRUCTORS */
		/* ---------------------------- */
		mscsFunctionLinSolve() : mscsFunction3dregc() {}
		virtual ~mscsFunctionLinSolve();
		
		mscsFunctionLinSolve(long Nx, long Ny, cpeds_VerbosityLevel verbosityLoc=CPEDS_defaultVerbosityLevel) : mscsFunction3dregc(Nx,Ny, 1, 1,1,1,0,0,0,CPEDS_defaultVerbosityLevel) {}
		mscsFunctionLinSolve(const mscsFunctionLinSolve& parent) : mscsFunction3dregc(parent) {}


		/* ---------------------------- */
		/* PUBLIC METHODS */
		/* ---------------------------- */

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
		mscsFunctionLinSolve_t _mscsFunctionLinSolve_data;
		
};
#endif /* MSCSFUNCTIONLINSOLVE_H_ */ 

