/*!
 * \file cpedsFunctionDistance.cpp
 *
 *  Created on: Oct 19, 2011
 *      Author: blew
 */

#include "cpedsFunctionDistance.h"

cpedsFunctionDistance::cpedsFunctionDistance() : _dist("distance") {

	_smallerSet=false;
}

cpedsFunctionDistance::~cpedsFunctionDistance() {
}

mscsFunction cpedsFunctionDistance::distance(mscsFunction& s1x, mscsFunction& s1y, mscsFunction& s2x, mscsFunction& s2y) {
	mscsFunction dx,dy;
	dx=calculateDistance(s1x,s2x,_s2xinterp);
	_s2xinterp.save("set2x-extra_interp.tmp");
	dy=calculateDistance(s1y,s2y,_s2yinterp);
	_s2yinterp.save("set2y-extra_interp.tmp");
	
	long N=dx.pointsCount();
	double dz,x,y;
	for (long i = 0; i < N; i++) {
		x=dx.f(i);
		y=dy.f(i);
		_dist.newPoint(dx.X(i),sqrt(x*x+y*y));
	}
	
	return _dist;
}

mscsFunction cpedsFunctionDistance::calculateDistance(mscsFunction &f1, mscsFunction &f2, mscsFunction &s2interp) {
	mscsFunction d;
	double *x = f1.extractArguments();
	s2interp=f2.extrapolate(f1.X(0),f1.last().rx(),f1.X(1)-f1.X(0),true,"linear");
	s2interp.interpolate(x,f1.pointsCount(),"linear",true);
	long N=f1.pointsCount();
	for (long i = 0; i < N; i++) {
		d.newPoint(f1.X(i),fabs(s2interp.f(i)-f1.f(i)));
	}
	delete [] x;
	return d;
}
