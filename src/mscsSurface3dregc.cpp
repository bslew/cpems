/*!
 * \file mscsSurface3dregc.cpp
 *
 *  Created on: Jan 7, 2012
 *      Author: blew
 */

#include "mscsSurface3dregc.h"

mscsSurface3dregc::mscsSurface3dregc() : mscsFunction3dregc() {


}

mscsSurface3dregc::~mscsSurface3dregc() {
}

/***************************************************************************************/
void mscsSurface3dregc::mkParaboloid(double fromX, double toX, double fromY, double toY, long Nx, long Ny, double A, double B, double vx, double vy, double vz, bool hyper) {
	double dx=fabs(toX-fromX)/Nx;
	double dy=fabs(toY-fromY)/Ny;
	
	setSize(Nx,Ny,1,dx,dy,1,fromX,fromY,0);
	allocFunctionSpace();
	
	double z,x,y;
	double tx,ty;
	
	if (hyper) {
		for (long i = 0; i < Nx; i++) {
//			x=i*dx-fromX;
			x=getX(i);
			tx=(x-vx)/A;
			for (long j = 0; j < Ny; j++) {
//				y=j*dy-fromY;
				y=getY(j);

				ty=(y-vy)/B;
				z= tx*tx - ty*ty + vz;
				setf(i,j,0,z,0);
			}
		}
	}
	else {
		for (long i = 0; i < Nx; i++) {
//			x=i*dx+fromX;
			x=getX(i);
			tx=(x-vx)/A;
			for (long j = 0; j < Ny; j++) {
//				y=j*dy+fromY;
				y=getY(j);

				ty=(y-vy)/B;
				z= tx*tx + ty*ty + vz;
				setf(i,j,0,z,0);
			}
		}		
	}
	
}


/***************************************************************************************/
void mscsSurface3dregc::mkHyperboloid(double fromX, double toX, double fromY, double toY, long Nx, long Ny, double A, double B, double C, double vx, double vy, double vz) {
	double dx=fabs(toX-fromX)/Nx;
	double dy=fabs(toY-fromY)/Ny;
	
	setSize(Nx,Ny,1,dx,dy,1,fromX,fromY,0);
	allocFunctionSpace();
	
	double z,x,y;
	double tx,ty;
	
	for (long i = 0; i < Nx; i++) {
//		x=i*dx+fromX;
		x=getX(i);
		tx=(x-vx)/A;
		for (long j = 0; j < Ny; j++) {
//			y=j*dy+fromY;
			y=getY(j);
			ty=(y-vy)/B;
			z= C*sqrt(1.0 + tx*tx + ty*ty) -C + vz;
			setf(i,j,0,z,0);
		}
	}
}
/***************************************************************************************/


/***************************************************************************************/
mscsFunction3dregc& mscsSurface3dregc::mkSin3D(double from, double to, double dr, double T, double phi, int dir) {
	double x,y,z;
	long i,j,k;

	// check the sizes of the output array
	x=0; i=0;	while (x<=to) {		x+=dr;		i++;	}
//	y=0; j=0;	while (y<=to) {		y+=dr;		j++;	}
//	z=0; k=0;	while (z<=to) {		z+=dr;		k++;	}

	setSize(i,i,i,dr,dr,dr,from,from,from);
	realloc();
	double wx=0,wy=0,wz=0;
	
	if (dir==0) wx=1;
	if (dir==1) wy=1;
	if (dir==2) wz=1;
	// make function
	x=from; i=0;
	while (x<=to) {
		y=from; j=0;
		while (y<=to) {
			z=from; k=0;
			while (z<=to) {
				fRe(i,j,k)=sin(twoPI/T*(wx*x + wy*y + wz*z));
				fIm(i,j,k)=0;
				z+=dr;
				k++;
			}
			y+=dr;
			j++;
		}
		x+=dr;
		i++;
	}
	
	return *this;	
}
/***************************************************************************************/
mscsSurface3dregc& mscsSurface3dregc::operator=(const mscsSurface3dregc& rhs) {
	mscsFunction3dregc::operator =(rhs);
	return *this;
}
/***************************************************************************************/
mscsSurface3dregc& mscsSurface3dregc::operator=(const mscsFunction3dregc& rhs) {
	mscsFunction3dregc::operator =(rhs);
	return *this;	
}
/***************************************************************************************/
void mscsSurface3dregc::mkGauss2D(double fromX, double toX, double fromY, double toY, long  Nx, long  Ny, double sx, double sy, double sz, double vx, double vy, double vz) {
	double dx=fabs(toX-fromX)/Nx;
	double dy=fabs(toY-fromY)/Ny;
	
	setSize(Nx,Ny,1,dx,dy,1,fromX,fromY,0);
	allocFunctionSpace();
	
	double z,x,y;
	double tx,ty;
	double sqrt2=sqrt(2);
	
	for (long i = 0; i < Nx; i++) {
		x=getX(i);
		tx=(x-vx)/sx;
		for (long j = 0; j < Ny; j++) {
			y=getY(j);
			ty=(y-vy)/B;
			z= exp(-(tx*tx/2 + ty*ty/2)) + vz;
			setf(i,j,0,z,0);
		}
	}
	
}

/***************************************************************************************/
void mscsSurface3dregc::mkExponentialSymmetrical2D(double fromX, double toX, double fromY, double toY, long Nx, long Ny, double base, double scale, double vx, double vy) {
	double dx=fabs(toX-fromX)/Nx;
	double dy=fabs(toY-fromY)/Ny;
	
	setSize(Nx,Ny,1,dx,dy,1,fromX,fromY,0);
	allocFunctionSpace();
	
	double z,x,y;
	double tx,ty;
	
	for (long i = 0; i < Nx; i++) {
		x=getX(i);
		tx=(x-vx)/scale;
		for (long j = 0; j < Ny; j++) {
			y=getY(j);
			ty=(y-vy)/scale;
			z= pow(base,(tx*tx + ty*ty));
			setf(i,j,0,z,0);
		}
	}	
}

/***************************************************************************************/
mscsLine mscsSurface3dregc::getNormalLine(double x, double y) {
	double dfdx,dfdy;
	double z=fxy(x,y,0,0,&dfdx,&dfdy);
//	dfdx=0.5968811774295497; // debug: this is a value for ray on the edge of secondary
//	printf("dfdx (x=%.15lf): %.15lf\n",x,dfdx);
//	printf("dfdy (y=%.15lf): %.15lf\n",y,dfdy);
	mscsLine n(x,y,z,x+dfdx,y+dfdy,z-1);
	return n;
}
