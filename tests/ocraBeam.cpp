/*!
 * \file ocraBeam.cpp
 *
 *	A very simple beam calculator
 *
 *
 *  All length units are in meters.
 *
 *  Created on: Mar 21, 2011
 *      Author: blew
 */
#include <math.h>
#include "Mscs-function.h"

double dist(const QPointF& p1, const QPointF& p2) {
	double d=sqrt((p1.x()-p2.x())*(p1.x()-p2.x()) + (p1.y()-p2.y())*(p1.y()-p2.y()));
	return d;
}

mscsFunction calcScatterBeam(double h, double lambda, mscsFunction& primary, mscsFunction& plane) {
	long Nprim=primary.pointsCount();
	long Npl=plane.pointsCount();
	mscsFunction sum("sum");
	sum=plane;
	sum=double(0);
	
	for (long i = 0; i < Nprim; i++) {
		for (long j = 0; j < Npl; j++) {
//			Nlambda=long(dist(primary[i],plane[j])/lambda);
//			phase=dist(primary[i],plane[j])/lambda-Nlambda*lambda;
			
			sum.f(j)+=sin(twoPI/lambda*dist(primary[i],plane[j]));
		}
		printf("%li / %li\n",i,Nprim);
	}
	return sum;
}

mscsFunction calcReflectedBeam(double h, double lambda, mscsFunction& primary, mscsFunction& plane) {
	long Nprim=primary.pointsCount();
	long Npl=plane.pointsCount();
	mscsFunction sum("sum");
	sum=plane;
	sum=double(0);
	
	for (long i = 0; i < Nprim; i++) {
		sum.f(i)=sin(twoPI/lambda*dist(primary[i],plane[i]));
	}
	return sum;
}

int main(int argc, char** argv) {
	double dx=0.01;
	double r=16.0;
	double lambda=0.01; // 30 GHz
	double h;//=strtod(argv[1],NULL);
	cpedsMsgs msgs;
	
	mscsFunction primary("primary mirror");
	double f=11.2;
	primary.mkPowerLaw(-r,r,dx,1/(4*f),1,2);
	
	mscsFunction plane("Clouds plane");
	
//	mscsFunction phase("phase");
	mscsFunction sum("sum");
	
	double phase;

//	h=2000;
//	plane.mkConst(-r,r,dx,h);
//	sum=calcReflectedBeam(h,lambda, primary,plane);
//	sum.save("ocra-beam-test-reflect.sum--h_"+msgs.toStr(h));
//	sum.power(2);
//	sum.save("ocra-beam-test-reflect.power--h_"+msgs.toStr(h));
//	for (long i = 0; i < sum.pointsCount(); i++) {		printf("%lE\n",sum.getx(i)/h);		sum.setarg(i,atan(sum.getX(i)/h)*PI180inv);	}
//	sum.save("ocra-beam-test-reflect.power.ang--h_"+msgs.toStr(h));

	
	
//	h=2000;
//	plane.mkConst(-r,r,dx,h);
//	sum=calcScatterBeam(h,lambda, primary,plane);
//	sum.save("ocra-beam-test.sum--h_"+msgs.toStr(h));
//	sum.power(2);
//	sum.save("ocra-beam-test.power--h_"+msgs.toStr(h));
//	for (long i = 0; i < sum.pointsCount(); i++) {		printf("%lE\n",sum.getx(i)/h);		sum.setarg(i,atan(sum.getX(i)/h)*PI180inv);	}
//	sum.save("ocra-beam-test.power.ang--h_"+msgs.toStr(h));
//
//	plane.clearFunction();
//	h=10000;
//	plane.mkConst(-r,r,dx,h);
//	sum=calcScatterBeam(h,lambda, primary,plane);
//	sum.save("ocra-beam-test.sum--h_"+msgs.toStr(h));
//	sum.power(2);
//	sum.save("ocra-beam-test.power--h_"+msgs.toStr(h));
//	for (long i = 0; i < sum.pointsCount(); i++) {		printf("%lE\n",sum.getx(i)/h);		sum.setarg(i,atan(sum.getX(i)/h)*PI180inv);	}
//	sum.save("ocra-beam-test.power.ang--h_"+msgs.toStr(h));

	plane.clearFunction();
	h=1000000;
	plane.mkConst(-10*r,10*r,dx,h);
	sum=calcScatterBeam(h,lambda, primary,plane);
	sum.save("ocra-beam-test.sum--10r-h_"+msgs.toStr(h));
	sum.power(2);
	sum.save("ocra-beam-test.power--10r-h_"+msgs.toStr(h));
	for (long i = 0; i < sum.pointsCount(); i++) {		sum.setarg(i,atan(sum.getX(i)/h)*PI180inv);	}
	sum.save("ocra-beam-test.power.ang--10r-h_"+msgs.toStr(h));

	return 0;
}
