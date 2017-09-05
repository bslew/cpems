/*!
 * \file test-bicubic-interpolation.cpp
 *
 *  Created on: Jan 12, 2012
 *      Author: blew
 */

#include "cpeds-consts.h"
#include "cpeds-math.h"
#include "Mscs-function3dregc.h"

void testPtr(double* p) {
	printf("f+...: %lf\n",*p);
}

int main() {
	/*	We test function z= sin(x+y) on the range (x,y)=(0,1)x(0,1) */
	// case for L=1
	double f[4] = {0,sin(1),sin(2),sin(1) };
	double f1[4] = { 1,cos(1),cos(2),cos(1) };
	double f2[4] = { 1,cos(1),cos(2),cos(1) };
	double f12[4] = { 0,-sin(1),-sin(2),-sin(1) };

	// casse for L=PI
//	double f[4] = {0,0,0,0 };
//	double f1[4] = { 1,-1,1,-1 };
//	double f2[4] = { 1,-1,1,-1 };
//	double f12[4] = {0,0,0,0 };
//	testPtr(f1+0);
//	testPtr(f1+1);
//	testPtr(f1+2);
//	testPtr(f1+3);
//	exit(0);
	
	long k=0;
	long n=100;
//	double Lx=PI;
//	double Ly=PI;
	double Lx=1;
	double Ly=1;
	double dx=Lx/n;
	double dy=Ly/n;
	printf("dx: %lf, dy: %lf\n",dx,dy);
	mscsFunction3dregc fg(n,n,1,dx,dy,1,0,0,0);
	double fint, f1int, f2int;
	double x,y;
//	printf("(-2+10) % 10 = %i\n",(-2L+10)%10L);
//	printf("(2+10) % 10 = %i\n",(2L+10)%10L);
//	printf("(10+10) % 10 = %i\n",(10L+10)%10L);
	
	x=0;
	y=0;
	printf("\ninterpolating at (x,y)=(%lf, %lf)\n",x,y);
	cpeds_bicubic_interpolation(f,f1,f2,f12,0,Lx,0,Ly,x,y,fint,f1int,f2int);
	printf("interpolated value is: %lf\n\n",fint);

	x=0.5;
	y=0.5;
	printf("\ninterpolating at (x,y)=(%lf, %lf)\n",x,y);
	cpeds_bicubic_interpolation(f,f1,f2,f12,0,Lx,0,Ly,x,y,fint,f1int,f2int);
	printf("interpolated value is: %lf\n\n",fint);

	x=0;
	y=1;
	printf("\ninterpolating at (x,y)=(%lf, %lf)\n",x,y);
	cpeds_bicubic_interpolation(f,f1,f2,f12,0,Lx,0,Ly,x,y,fint,f1int,f2int);
	printf("interpolated value is: %lf\n\n",fint);

	x=1;
	y=0;
	printf("\ninterpolating at (x,y)=(%lf, %lf)\n",x,y);
	cpeds_bicubic_interpolation(f,f1,f2,f12,0,Lx,0,Ly,x,y,fint,f1int,f2int);
	printf("interpolated value is: %lf\n\n",fint);

	x=1;
	y=1;
	printf("\ninterpolating at (x,y)=(%lf, %lf)\n",x,y);
	cpeds_bicubic_interpolation(f,f1,f2,f12,0,Lx,0,Ly,x,y,fint,f1int,f2int);
	printf("interpolated value is: %lf\n\n",fint);

	x=0.7;
	y=0.8;
	printf("\ninterpolating at (x,y)=(%lf, %lf)\n",x,y);
	cpeds_bicubic_interpolation(f,f1,f2,f12,0,Lx,0,Ly,x,y,fint,f1int,f2int);
	printf("interpolated value is: %lf\n\n",fint);

	x=0.3;
	y=0.9;
	printf("\ninterpolating at (x,y)=(%lf, %lf)\n",x,y);
	cpeds_bicubic_interpolation(f,f1,f2,f12,0,Lx,0,Ly,x,y,fint,f1int,f2int);
	printf("interpolated value is: %lf\n\n",fint);
	
	//	exit(0);
	for (long i = 0; i < n; i++) {
		x=i*dx;
		for (long j = 0; j < n; j++) {
			y=j*dy;
			printf("interpolating at: %lf %lf\n",x,y);
			cpeds_bicubic_interpolation(f,f1,f2,f12,0,Lx,0,Ly,x,y,fint,f1int,f2int);
			
			fg.setf(i,j,0,fint,double(k++));
		}
	}
	
	cpeds_matrix_save(fg.getSlice(2,0,0),"test-bicubic-interpolation.dat");
	
	
}
