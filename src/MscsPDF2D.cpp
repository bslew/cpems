/*!
 * \file MscsPDF2D.cpp
 *
 *  Created on: May 31, 2017
 *      Author: blew
 */
#include "MscsPDF2D.h"

MscsPDF2D::MscsPDF2D() {
	_MscsPDF2D_data.contour_density=0.5;
}

MscsPDF2D::~MscsPDF2D() {
}

MscsPDF2D::MscsPDF2D(const mscsFunction3dregc& parent) : mscsFunction3dregc(parent) {
//	subtract(getMinValue());
//	normalize();
//	divide(getMaxValue());
	_MscsPDF2D_data.normalization=integrateRe();
}

mscsFunction MscsPDF2D::getContour(double CL, double* LVL) {
	mscsFunction c;
	mscsFunction cdf=getCDF();
	mscsFunction cdfi=cdf;
	cdfi.invert();
	cdfi.sortFunctionArgAscending();
	cdfi.average_sameArgs(0);
//	cdfi.save("cdfi");
	double lvl=cdfi.finter(CL);
//	printf("min: %lE, max: %lE\n",getMinValue(),getMaxValue());
//	printf("I(lvl=%lf) = %lE\n",0.0,getIntegral(0.0));
//	printf("I(lvl=%lf) = %lE\n",0.5,getIntegral(0.5));
//	printf("I(lvl=%lf) = %lE\n",lvl,getIntegral(lvl)/getIntegral(0));
	
	double xmin,xmax,ymin,ymax,x,y,dx,dy,dfx,dfy,tmp,df;
	xmin=getMinX();
	xmax=getMaxX();
	ymin=getMinY();
	ymax=getMaxY();
	dx=getDx()/2;
	dy=getDy()/2;
	x=xmin;
	while (x<xmax) {
		y=ymin;
		while (y<ymax) { 
//			printf("x %lf, y %lf, f %lf, dfx %lE, dfy %lE\n",x,y,tmp,dfx,dfy);
			tmp=fxy(x,y,0,0,&dfx,&dfy,false);
			df=cpeds_get_max(fabs(dfx*getDx()),fabs(dfy*getDy()))/2;
			if (fabs(tmp-lvl)<=df) c.newPoint(x,y);
			y+=dy;
		}
		x+=dx;
	}
	
	
	
	
	if (LVL!=NULL) *LVL=lvl;
	
	return c;
}
/***************************************************************************************/
mscsFunction MscsPDF2D::getCDF() {
	mscsFunction cdf;
	
	double zMin=0.0; // getMinValue();
	double zMax=getMaxValue();
	double dz=(zMax-zMin)/1000;
	double z=zMax;
	do {
		z-=dz;
		if (z<zMin) z=zMin;
		cdf.newPoint(z,getIntegral(z));
	} while (z>zMin);
	
//	cdf.invert();
//	cdf.average_sameArgs(0);
//	cdf.invert();
//	cdf.sortFunctionArgAscending();
	cdf.checkRanges();
	cdf/=cdf.getMaxValue();
//	cdf.save("cdf");
	return cdf;
}
/***************************************************************************************/
double MscsPDF2D::getIntegral(double zMin) {
	double I=0;
	long Ngtr=0;
	
	for (long i = 0; i < Nx(); i++) {
		for (long j = 0; j < Ny(); j++) {
			if (fRe(i,j,0)>=zMin) {
				I+=fRe(i,j,0);
//				I++;
				Ngtr++;
			}
		}
	}
	I*=(getDx()*getDy());
//	cout << "z: " << zMin << " Ngtr: " << Ngtr << " int: "<< I << "\n";
	
	return I;
}

/***************************************************************************************/
MscsPDF2D& MscsPDF2D::operator=(const mscsFunction3dregc& rhs) {
	mscsFunction3dregc::operator=(rhs);
	return *this;
}
/***************************************************************************************/
void MscsPDF2D::normalize() {
	normalizeRe();
}
