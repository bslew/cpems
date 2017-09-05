/*!
 * \file readGraficIC.cpp
 *
 *  Created on: Jan 4, 2011
 *      Author: blew
 */

#include "Nbody_io.h"


int main(int argc, char** argv) {
	
	Nbody_io ic;
		
	if (argc==1) { 
		printf("USAGE: lss_convert_IC_grafic2Gadget dir_with_grafic_IC gadgetICfileName\n");
		printf("\n");
		exit(0);
	}

	ic.loadGraficICsDensity(argv[1]);
	matrix<double> slice;
	slice=ic.getICDensityData().getSlice(0,0,0);
	cpeds_matrix_save(slice,"densityTestSlice-0-0.tmp");
	
	ic.convertGraficToGadgetIC(argv[1],argv[2]);
	return 0;
}
