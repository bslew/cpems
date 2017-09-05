/*!
 * \file readGraficIC.cpp
 *
 *  Created on: Jan 4, 2011
 *      Author: blew
 */

#include "Nbody_io.h"


int main(int argc, char** argv) {
	
	Nbody_io ic;
		

	ic.loadGraficICsDensity(argv[1]);
	matrix<double> slice;
	slice=ic.getICDensityData().getSlice(0,0,0);
	cpeds_matrix_save(slice,"densitySlice");
	
	ic.convertGraficToGadgetIC(argv[1],"gadgetIC");
	
	
}
