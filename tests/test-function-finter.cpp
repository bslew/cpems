/*!
 * \file test-function-finter.cpp
 *
 *  Created on: May 22, 2012
 *      Author: blew
 */
#include "Mscs-function.h"
#include <string.h>
#include "cpeds-msgs.h"

int main() {
	mscsFunction _detectorNL;
	cpedsMsgs* msgs = new cpedsMsgs("test-function-finter");
	string NLfile="/home/blew/projects/ocra-f-simulations/nonlinearity-fit/gieniu-data/compilation.txt";	
	_detectorNL.load(NLfile);
	if (_detectorNL.pointsCount()==0) { msgs->criticalError(string("I didn't find the non-linearity function for the detecors (")+NLfile+"). Please correct the class and recompile.",Top); }

	_detectorNL.print();
	_detectorNL.save("test-function-finter.in");
	mscsFunction nlint;
	double p=0;
	double dp=1e-4;
	
	while (p<1) {
//		printf("interpolating: p=%lf\n",p);
		if (p<0.316228) {
			nlint.newPoint(p,_detectorNL.finter(p,"cspline"));
		}
		else
			nlint.newPoint(p,_detectorNL.finter(p,"linear"));
		p+=dp;
//		printf("\n");
	}
	nlint.save("test-function-finter.out");
	delete msgs;
	return 0;
	
}
