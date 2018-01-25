/*!
 * \file MCparameterSpace.cpp
 *
 *  Created on: Jan 24, 2018
 *      Author: blew
 */
#include <MCparameterSpace.h>

MCparameterSpace::MCparameterSpace() {
	// TODO Auto-generated constructor stub
	
}

MCparameterSpace::~MCparameterSpace() {
	// TODO Auto-generated destructor stub
}

void MCparameterSpace::addParameter(const mscsFunction& p, string parameter_name, string parameter_full_name) {
	append(p);
	last().checkRanges();

	names_latex().append(parameter_full_name);
	names().append(parameter_name);
	
}

const MCparameterSpace& MCparameterSpace::operator =(const MCparameterSpace& rhs) {
	if (this!=&rhs) {
		QList<mscsFunction>::operator=(rhs);
		_MCparameterSpace_data=rhs._MCparameterSpace_data;
	}
	return *this;
}
