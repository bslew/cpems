/*!
 \file cpedsMCMCset.h - 
 */

#ifndef CPEDSMCMCSET_H_
#define CPEDSMCMCSET_H_

#include "cpeds-msgs.h"
#include "cpedsMCMC.h"

/*!
 \class cpedsMCMCset
 \brief Encapsulates a collection of cpedsMCMC objects and performs marginalization
 \details 
 
 \date created: May 24, 2017, 10:59:43 AM 
 \author Bartosz Lew
 */
class cpedsMCMCset : public cpedsMsgs {
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PUBLIC MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	public:
		
		/* ------------- */
		/* CLASS TYPES */
		/* ------------- */
		typedef struct {
			cpedsMCMC bflinks; //!< chain containing only the best links from different chains
			cpedsMCMC alllinks; //!< chain containing all links from different chains 
		} cpedsMCMCset_t;

		/* ---------------------------- */
		/* CONSTRUCTORS AND DESTRUCTORS */
		/* ---------------------------- */
		cpedsMCMCset(cpeds_VerbosityLevel verb=Zero);
		~cpedsMCMCset();

		/* ---------------------------- */
		/* PUBLIC METHODS */
		/* ---------------------------- */
		void addChain(cpedsMCMC& chain);

		
		
		
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
		cpedsMCMCset_t _cpedsMCMCset_data;
		
};
#endif /* CPEDSMCMCSET_H_ */ 

