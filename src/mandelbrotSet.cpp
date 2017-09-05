/*!
 * \file mandelbrotSet.cpp
 *
 *  Created on: Oct 18, 2012
 *      Author: blew
 */

#include "mandelbrotSet.h"

mandelbrotSet::mandelbrotSet(long Nx, long Ny, double xmin, double xmax, double ymin, double ymax, long maxIterations, cpeds_VerbosityLevel verbosityLoc) : mscsFunction3dregc("mandelbrot set", verbosityLoc) {
	setSizeRange(Nx, Ny, 1, xmin, ymin, -1, xmax, ymax, 1);
	allocFunctionSpace();
	_maxIterations=maxIterations;
}

mandelbrotSet::~mandelbrotSet() {
	freeFunctionSpace();
}
/***************************************************************************************/
void mandelbrotSet::mkSet() {
	find_set();
}
/***************************************************************************************/
void mandelbrotSet::find_set() {
	long i,j,color;
	double x,y;
	
	for (i=0;i<Nx();i++) {
		x = getX(i);
		for (j=0;j<Ny();j++) { // for each point in the plane
			y = getY(j);
			fIm(i,j,0)=double(check_point(x,y,_maxIterations, &color));
			fRe(i,j,0)=double(color);
		}
		msgs->say("%i/%i",i,Nx(),Low);
	}
}
/***************************************************************************************/
int mandelbrotSet::check_point(double x, double y, long iter_num, long *color) {
	MScomplex C,Zn;
	//complex *Z = new complex[iter_num];
	double d;
	long i;
	int ok;
	
	C.R = x; C.I = y;
	//Z[0].R = 0; Z[0].I = 0; i=0;       d=complex_abs2(Z[i]); 
	Zn.R = 0; Zn.I = 0;    d=complex_abs2(Zn);  i=0;
	while ((i < iter_num) && (d < 4)) { 
		i++;
		Zn = complex_sum(complex_multiply(Zn,Zn),C); // for Mandelbrot
		d = complex_abs2(Zn);
	}
	if (d < 4) { ok = 1; } else { ok = 0; }
	*color = i; 
	return ok;
}
/***************************************************************************************/
mandelbrotSet::MScomplex mandelbrotSet::complex_sum(mandelbrotSet::MScomplex C1, mandelbrotSet::MScomplex C2) { 
	MScomplex C; 
	C.R = C1.R+C2.R;    C.I = C1.I+C2.I; 
	return C;
}
mandelbrotSet::MScomplex mandelbrotSet::complex_multiply(mandelbrotSet::MScomplex C1, mandelbrotSet::MScomplex C2) {
	MScomplex C;
	C.R = C1.R*C2.R-C1.I*C2.I;    C.I = C1.R*C2.I+C1.I*C2.R;
	return C;
}

double mandelbrotSet::complex_abs2(mandelbrotSet::MScomplex C) {
	double d;
	d = C.R*C.R+C.I*C.I;    
	return d;
}
