/*!
 * \file test-fft-alias.cpp
 *
 *  Created on: Dec 14, 2010
 *      Author: blew
 */

#include "Mscs-function.h"

int main(int argc, char** argv){
		
	mscsFunction f,p;
	
//	f.load(argv[1]);

	
//	f.mkSin(0,10,0.1,1,0);
//	f.save("sin-0-10-0.1-1.txt");	
//	p=f.powerSpectrum(0.01,20,0.01);
//	p.save("sin-0-10-0.1-1--power--kmin_0.01-kmax_20-dk_0.01.txt");
//
//	f.clearFunction();
//	
//	f.mkSin(0,10,0.2,1,0);
//	f.save("sin-0-10-0.2-1.txt");	
//	p=f.powerSpectrum(0.01,20,0.01);
//	p.save("sin-0-10-0.2-1--power--kmin_0.01-kmax_20-dk_0.01.txt");
//	
//	f.clearFunction();
//	
//	f.mkSin(0,10,0.3,1,0);
//	f.save("sin-0-10-0.3-1.txt");	
//	p=f.powerSpectrum(0.01,20,0.01);
//	p.save("sin-0-10-0.3-1--power--kmin_0.01-kmax_20-dk_0.01.txt");
//	
//	f.clearFunction();
//	
//	f.mkSin(0,10,0.4,1,0);
//	f.save("sin-0-10-0.4-1.txt");	
//	p=f.powerSpectrum(0.01,20,0.01);
//	p.save("sin-0-10-0.4-1--power--kmin_0.01-kmax_20-dk_0.01.txt");
//
//	f.clearFunction();
//	
//	f.mkSin(0,10,0.8,1,0);
//	f.save("sin-0-10-0.8-1.txt");	
//	p=f.powerSpectrum(0.01,20,0.01);
//	p.save("sin-0-10-0.8-1--power--kmin_0.01-kmax_20-dk_0.01.txt");
	

//	f.clearFunction();
//	
//	f.mkSin(0,1,0.01234,0.001,0);
//	f.save("sin-0-1-0.01234-0.001.txt");	
//	p=f.powerSpectrum(1,2000,.1);
//	p.save("sin-0-1-0.01234-0.001--power--kmin_1-kmax_2000-dk_0.1.txt");
	
	
	
	
	
	
	
	
	
	
	
//	mscsFunction g;
//	
	double fs=277.0;
	f.mkSin(0,0.2,1.0/fs,0.02,0);
	f.save("sin-0-0.2-0.00361-0.02.txt");	
	p=f.powerSpectrum(0.0,3*fs,1);
	p.save("sin-0-0.2-0.00361-0.02--power--kmin_0.0-kmax_3fs-dk_1.txt");

//	f.clearFunction();
//
//	f.mkSin(0,.1,0.00361,0.02,0);
//	f.save("sin-0-0.1-0.00361-0.02.txt");	
//	p=f.powerSpectrum(0.00361,277,1);
//	p.save("sin-0-0.1-0.00361-0.02--power--kmin_0.00361-kmax_277-dk_1.txt");
//
//	f.clearFunction();

	
	//	f.mkSin(0,1,0.00361,0.02,0);
//	f.save("sin-0-1-0.00361-0.02.txt");	
//	p=f.powerSpectrum(0.00361,137,0.00361);
//	p.save("sin-0-1-0.00361-0.02--power--kmin_0.00361-kmax_137-dk_0.00361.txt");
//
////	f.clearFunction();
//
//	g.mkSin(0,1,0.00361,0.01,0);
//	g.save("sin-0-1-0.00361-0.01.txt");	
//	p=g.powerSpectrum(0.00361,137,0.00361);
//	p.save("sin-0-1-0.00361-0.01--power--kmin_0.00361-kmax_137-dk_0.00361.txt");
//
////	f.clearFunction();
//	f+=g;
//
////	f.mkSin(0,1,0.00361,0.01,0);
//	f.save("sin-0-1-0.00361-0.01_sum_0.02.txt");	
//	p=f.powerSpectrum(0.00361,137,0.00361);
//	p.save("sin-0-1-0.00361-0.01_sum_0.02--power--kmin_0.00361-kmax_137-dk_0.00361.txt");
//	
}
