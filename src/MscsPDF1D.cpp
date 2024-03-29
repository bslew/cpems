/*!
 * \file MscsPDF1D.cpp
 *
 *  Created on: Jun 6, 2017
 *      Author: blew
 */
#include "MscsPDF1D.h"
#include <assert.h>

/***************************************************************************************/
MscsPDF1D::MscsPDF1D()  {
	
}

/***************************************************************************************/
MscsPDF1D::~MscsPDF1D() {
}

/***************************************************************************************/
MscsPDF1D::MscsPDF1D(const mscsFunction& parent) : mscsFunction(parent) {
	
}

/***************************************************************************************/
MscsPDF1D& MscsPDF1D::operator=(const mscsFunction& rhs) {
	mscsFunction::operator=(rhs);
	return *this;
}
/***************************************************************************************/
cpedsList<double> MscsPDF1D::getCR(double CL, double* LVL) {
	mscsFunction c;
	mscsFunction cdf=getCDF();
	mscsFunction cdfi=cdf;
#ifdef DEBUG_MCMC_PDF
	cdf.save("rusty.cdf");
#endif
	cdfi.invert();
	cdfi.sortFunctionArgAscending();
#ifdef DEBUG_MCMC_PDF
	cdfi.save("rusty.icdf");
#endif
	cdfi.average_sameArgs();
	double lvl=cdfi.finter(CL);
	
	if (LVL!=NULL) *LVL=lvl;

	mscsFunction tmp(*this);
	tmp-=lvl;
	vector<double> CR=tmp.findRoot();
	cpedsList<double> CRl(CR);
	checkRanges();
	double modalX=getX(getMaxValueIdx());
	
	//
	// sanity checks
	//
	if (CR.size()==1) {
		/*
		 * we probably have ill defined PDF that is not probed sufficiently for this CL.
		 * We need to find out which tail of the PDF is missing and return the extremal argument
		 * as CR boundary.
		 *
		 */
		
		if (CR[0]-modalX > 0) CRl.prepend(getMinArg());
		else CRl.append(getMaxArg());
	}
//	else {
//		if (f(getMinArg()>lvl) or f(getMaxArg()>lvl)) {
//			if (CR[0]-modalX > 0) CRl.prepend(getMinArg());
//			else CRl.append(getMaxArg());
//		}		
//	}
	
	return CRl;
}
/***************************************************************************************/
mscsFunction MscsPDF1D::getCDF() {
	mscsFunction cdf;
//	_MscsPDF1D_data.CDF_sel=*this;
//	_MscsPDF1D_data.CDF_last=0;
	
	checkRanges();
	double yMin=0.0; // getMinValue();
	double yMax=getMaxValue();
	double dy=(yMax-yMin)/1000;
	double y=yMax;
//	_MscsPDF1D_data.twoSidedIntegral_lastyMin=yMax;
	do {
		y-=dy;
		if (y<yMin) y=yMin;
//		_MscsPDF1D_data.CDF_last=;
		cdf.newPoint(y,getTwosidedIntegral(y));
	} while (y>yMin);
	
	cdf.checkRanges();
	cdf/=cdf.getMaxValue();
	return cdf;
}
/***************************************************************************************/
double MscsPDF1D::getTwosidedIntegral(double yMin) {
	double I=0;
	long N=pointsCount();
	assert(N>2);
//	mscsFunction* tmp;
//	
//	if (yMin<_MscsPDF1D_data.twoSidedIntegral_lastyMin) {
//		_MscsPDF1D_data.CDF_sel.removeLarger(_MscsPDF1D_data.twoSidedIntegral_lastyMin);
//		N=_MscsPDF1D_data.CDF_sel.pointsCount();
//		assert(N>2);
//		for (long i = 0; i < N; i++) {
//			double dx;
//			if (i==0) dx=(getx(1)-getx(0))/2;
//			else {
//				if (i==N-1) {
//					dx=(getx(N-1)-getx(N-2))/2;
//				}
//				else { dx=(getx(i)-getx(i-1))/2 + (getx(i+1)-getx(i))/2; }
//			}
//			I+=f(i)*dx;
//		}
//	}
//	else {
		for (long i = 0; i < N; i++) {
			if (f(i)>=yMin) {
				double dx;
				if (i==0) dx=(getx(1)-getx(0))/2;
				else {
					if (i==N-1) {
						dx=(getx(N-1)-getx(N-2))/2;
					}
					else { dx=(getx(i)-getx(i-1))/2 + (getx(i+1)-getx(i))/2; }
				}
				I+=f(i)*dx;
	//			Ngtr++;
			}
		}		
//	}
	
//	_MscsPDF1D_data.CDF_last=I;
	return I;
}

