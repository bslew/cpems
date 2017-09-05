/*!
 * \file test-function-fft.cpp
 *
 *  Created on: Jul 22, 2010
 *      Author: blew
 */

#include <stdlib.h>
#include <stdio.h>
#include "Mscs-function.h"

void printUsage() {
	printf("USAGE: test-function-fft XYfile regridFactor outfile\n\n where:\n\n regridFactor is 0 for unifirm spacing assumption, and <1 for fraction of mean X points separations at which to linearly interpolate the function\n");	
}

int main(int argc, char** argv) {
	char tmpch[1000];
			//	cpedsMsgs *msgs = new cpedsMsgs("test");
			//	msgs->say("test",High);
			//	delete msgs;
			//	exit(0);
	mscsFunction f("f"),c("copy");
	mscsFunction re,im,P;

	{ // test polynomials
		f.mkPolynomial(-10,10,0.01,3, 1.0, 2.0, 3.0, -4.0);
		f.save("parabola.dat");
		exit(0);
	}
	
	
//	if (argc==1) { printUsage(); exit(0); }
	if (argc==1) {
		
		printUsage();
		
		c=f.mkSin(0,1,0.0001,0.01);
		
		for (long i=0;i<f.pointsCount();i++) { f.f(i)+=10; 	i+=500; }
		f.save("/home/blew/tmp/input_sin_signal");
		
		f.extractSpikes(500,250,3);
		f.save("/home/blew/tmp/input_sin_signal-no-spikes");
//		exit(0);

		P=f.powerSpectrum(&re,&im,0);		
		P.save("/home/blew/tmp/input_sin_signal.power");
		f=c;
		cpedsList<double> w;
		cpedsList<long> b;
		int i=0,N=f.pointsCount();
		int bs=100;
		while (i<N) {
			b.append(bs); 
//			w.append(1);
			i+=bs;
		}
		f.binFunction(long(0),b, w);
		f.save("/home/blew/tmp/sampled_signal");
		P=f.powerSpectrum(&re,&im,0);		
		P.save("/home/blew/tmp/sampled_signal.power");
		
//		P=f.powerSpectrum(&re,&im,0);		
//		P.save("/home/blew/tmp/test-function-fft-regrid_0.power");
//		f=c;
//		P=f.powerSpectrum(&re,&im,1);		
//		P.save("/home/blew/tmp/test-function-fft-regrid_1.power");
//		re.save("/home/blew/tmp/test-function-fft.re");
//		im.save("/home/blew/tmp/test-function-fft.im");
	}
	else {
		if (argc!=3 and argc!=4) { printUsage(); printf("arg num: %li\n",argc); exit(0); }
		f.load(argv[1]);
		f-=f.meanf();
		double regridX=strtod(argv[2],NULL);
		P=f.powerSpectrum(&re,&im,regridX);
		if (argc==4) P.save(argv[3]);
		else {
			sprintf(tmpch,"%s-regrid_%.2lf.power",argv[1],regridX);
			P.save(tmpch);
		}
		
	}
	
	return 0;
}
