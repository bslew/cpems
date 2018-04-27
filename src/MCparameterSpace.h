/*!
 \file MCparameterSpace.h - 
 */

#ifndef SRC_MCPARAMETERSPACE_H_
#define SRC_MCPARAMETERSPACE_H_

#include <QtCore/QList>
#include "Mscs-function.h"

/*!
 \class MCparameterSpace
 \brief Encapsulates parameter space for MCMC run
 \details 
 
 \date created: Jan 24, 2018, 1:15:58 PM 
 \author Bartosz Lew
 */
class MCparameterSpace : public QList<mscsFunction> {
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PUBLIC MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	public:
		
		/* ------------- */
		/* CLASS TYPES */
		/* ------------- */
		typedef struct {
			QList<string> names; //!< parameter name
			QList<string> names_full; //!< parameter name with units in latex typesetting used for automatic plot generation
				
		} MCparameterSpace_t;

		/* ---------------------------- */
		/* CONSTRUCTORS AND DESTRUCTORS */
		/* ---------------------------- */
		MCparameterSpace();
		~MCparameterSpace();

		/* ---------------------------- */
		/* PUBLIC METHODS */
		/* ---------------------------- */

		void addParameter(const mscsFunction& p, string parameter_name="", string parameter_full_name="");
		QList<string>& names_latex() { return _MCparameterSpace_data.names_full; }
		QList<string>& names_full() { return _MCparameterSpace_data.names_full; }
		QList<string>& names() { return _MCparameterSpace_data.names; }
		
		
		const MCparameterSpace& operator=(const MCparameterSpace& rhs);

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
		MCparameterSpace_t _MCparameterSpace_data;
		
};
#endif /* SRC_MCPARAMETERSPACE_H_ */ 

