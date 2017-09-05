/*!
 * \file Mscs-map-interpolate.cpp
 *
 *  Created on: Nov 7, 2012
 *      Author: blew
 */

#include "Mscs-map.h"
#include "cpeds-project.h"
#include "armadillo"


void mscsMap::interpolate(mscsMap& outMap, double R, double r0,  string rbfName) {
	long toNs=outMap.nside();
	cpedsList<long> srcIdx, destIdx;
	cpedsDirectionSet dsDest,dsSrc;
	cpedsDirection nt; // tangent plane direction
	cpedsPointSet3D psSrc;
	double yInt; // the interpolated value
	
	r0*=PI180;
	R*=PI180;
	
	// initiate source map directions
	set_map_coord();
	// initiate target map directions and values
	outMap.set_map_coord();
	if (outMap.mapLoaded()==false) {
		outMap.makekill_space_manager("make","T");
		outMap.clear_map();
	}
	
	// select target directions 
	if (outMap.maskLoaded()) {
		for (long i = 0; i < outMap.pixNum(); i++) {
			if (outMap.isMasked(i)==false) {
				dsDest.append(outMap.get_C(i));
				destIdx.append(i);
			}
		}
	}
	else {
		for (long i = 0; i < outMap.pixNum(); i++) {
			dsDest.append(outMap.get_C(i));
			destIdx.append(i);
		}
	}
	if (dsDest.size()==0) { msgs->warning("no target directions found. Perhaps all sky is masked ?", High); }
	
	if (maskLoaded()==false) {
		clean_mask();
	}
	//
	// for every target direction perform the interpolation using the unmasked neighboring directions from this object
	//
	bool res=true;
	double normv=0;
	
	for (long i = 0; i < dsDest.count(); i++) {
		nt=dsDest[i];
//		nt.print_direction("nt");
		dsSrc.clear();
		srcIdx.clear();
		psSrc.clear();
		// select the neighboring points from source data within the distance r1 from tangent direction
		for (long j = 0; j < pixNum(); j++) {
			if (isMasked(j)==false) {
				if (cpeds_ang_n1n2(nt.get_direction(),n(j).get_direction()) <= R) {
					dsSrc.append(n(j));
					srcIdx.append(j);
				}
			}
		}
		
		// project nt and dsSrc onto a tangent plane
		nt*=PI180inv;
//		dsSrc.save("dsSrc.txt");
		psSrc.append( cpedsProject(dsSrc, "stere").projectOnPlane(nt) );
//		psDest.save("projectedPoints.txt");
//		exit(0);
		
		//
		//
		// interpolate
		//
		//
		double (*rbfP)(double, double);
		if (rbfName=="gaussian") rbfP=&mscsFunction::FN_gaussian1D_A1m0;
		if (rbfName=="multiquad") rbfP=&mscsFunction::FN_multiquadratic1D;
		
		//
		// derive w_i coefficients
		//
		// srcY * Vnorm = PHI * W
		// where: 
		// srcY = temperature values at selected neighboring points
		// Vnorm - W coefficients normalization, Vnorm = Sum_i=0^N-1 phi_ji, where phi_ji = phi(|vec(x_j) - vec(x_i)|) and phi is a radial basis function
		// PHI - matrix of phi_ji values
		// W - vector of weights to be found
		
		// generate PHI matrix of distances between pairs of points
		arma::mat PHI(psSrc.size(),psSrc.size());
		arma::mat Vnorm(psSrc.size(),1);
		arma::mat srcY(psSrc.size(),1);
		arma::mat W(psSrc.size(),1);
//		cpedsPoint3D pnt(0,0,0);
		double r=0;
		for (long k = 0; k < psSrc.size(); k++) {
			Vnorm(k,0)=0;
			for (long l = 0; l < psSrc.size(); l++) {
				r=psSrc[k].dist(psSrc[l]);
				PHI(k,l)=rbfP(r,r0);
				Vnorm(k,0)+=PHI(k,l);
			}
			srcY(k,0)=get_T(srcIdx[k]);
		}

		// solve the set of linear equations to get the weights
		srcY%=Vnorm;
		res=arma::solve(W,PHI,srcY);
		if (res==false) { printf("solution not found\n"); exit(0); }
		// derive the interpolated value at point (0,0,0)
		yInt=0;
		normv=0;
		for (long j = 0; j < psSrc.size(); j++) {
			yInt+=W(j,0)*rbfP(psSrc[j].dist(),r0);
			normv+=rbfP(psSrc[j].dist(),r0);
		}
		yInt/=normv;
		
		// substitute the interpolated value y to the correct map pixel
		outMap.set_T(destIdx[i],yInt);
	}
	
	
}
