/*!
 * \file cpedsMCMCset.cpp
 *
 *  Created on: May 24, 2017
 *      Author: blew
 */
#include "cpedsMCMCset.h"

/***************************************************************************************/
cpedsMCMCset::cpedsMCMCset(cpeds_VerbosityLevel verb) : cpedsMsgs("cpedsMCMCset", false,"",verb) {
	
}
/***************************************************************************************/
cpedsMCMCset::~cpedsMCMCset() {}
/***************************************************************************************/
void cpedsMCMCset::addChain(cpedsMCMC& chain) {
	_cpedsMCMCset_data.alllinks.append(chain);
	_cpedsMCMCset_data.bflinks.append(chain.bestFitLink());
}
