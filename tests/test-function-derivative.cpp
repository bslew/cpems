/*!
 * \file test-function-derivative.cpp
 *
 *  Created on: Jul 22, 2010
 *      Author: blew
 */

#include <stdlib.h>
#include <stdio.h>
#include "Mscs-function.h"

void printUsage() {
  printf("USAGE: test-function-derivative ");
}

int main(int argc, char** argv) {
	char tmpch[1000];
	mscsFunction f("f"),g("g");

//	if (argc==1) { printUsage(); exit(0); }
	if (argc==1) {
		
		printUsage();
		for (long i = 0; i < 100; i++) {
			f.newPoint(i,10-(i%10));
		}
//		f.mkSin(0,1,0.01,0.1,0,1);
//		f.absoluteValue();
		f.save("/home/blew/tmp/input_sin_signal");
		double Py=10;
		g=f.derivative(false,NULL,NULL);
		g.print();
		g.save("/home/blew/tmp/input_sin_signal-deriative");
		g=f.derivative(false,NULL,&Py);
		g.print();
		g.save("/home/blew/tmp/input_sin_signal-deriative2");

	}
	
	return 0;
}
