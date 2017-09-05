/*!
 * \file testoofNumbers.cpp
 *
 *  Created on: Nov 4, 2010
 *      Author: blew
 */

#include "cpeds-list.h"
#include "cpeds-rng.h"

int main(int argc, char** argv) {
	cpedsRNG* rns=new cpedsRNG("gaussianPowerLaw_t","double");
	rns->seed(1);
	long N;
	N=strtol(argv[1],NULL,10);
//	long N=10;
//	cpedsList<double> t=rns->getRNs(N);
//	cpedsList<double> t2=rns->getRNs(N);
	long i=0;
	double t;
	while (i<N) {
		t=rns->getRN();
		printf("i: %li, t: %lE\n",i,t);
		i++;
	}

//	printf("saving\n");
//	t.save("1ofnumbers");
//	t2.save("1ofnumbers2");
//	t.print();
//	t2.print();
	
	return 0;
}
