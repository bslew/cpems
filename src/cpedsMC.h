/*!
 \file cpedsMC.h - 
 */

#ifndef SRC_CPEDSMC_H_
#define SRC_CPEDSMC_H_

//#include <QtCore/QList>
#include "MClink.h"


/*!
 \class cpedsMC
 \brief Encapsulates 
 \details 
 
 \date created: Jan 23, 2018, 4:09:22 PM 
 \author Bartosz Lew
 */
class cpedsMC : public QList<MClink> {
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PUBLIC MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	public:
		
		/* ------------- */
		/* CLASS TYPES */
		/* ------------- */
		typedef struct {
				
		} cpedsMC_t;

		/* ---------------------------- */
		/* CONSTRUCTORS AND DESTRUCTORS */
		/* ---------------------------- */
		cpedsMC();
		~cpedsMC();

		/* ---------------------------- */
		/* PUBLIC METHODS */
		/* ---------------------------- */
		void append(MClink l, long idx);
		void append(MClink l);
		void append(cpedsMC& mc);
		void load(string fname);
		const cpedsMC& operator=(const cpedsMC& rhs);

		void save(string fname, long precision=5);
		
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
		cpedsMC_t _cpedsMC_data;
		
};
#endif /* SRC_CPEDSMC_H_ */ 

