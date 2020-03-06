/*!
 * \file test-SPinterpolation-2d.cpp
 *
 *  Created on: Feb 28, 2020
 *      Author: blew
 */


#include <stdio.h>
#include <string.h>
#include "Mscs-function3dregc.h"
#include "pointsDensity.h"

int main(int argc, char** argv) {
	
	mscsFunction3dregc f("field");

	
	mscsVector<double> vals;
	pointsDensity dens;
	dens.append(cpedsPoint3D(0.5,0.5,0)); 		vals.push_back(1);
	dens.append(cpedsPoint3D(5.5,5.5,0)); 	vals.push_back(1);
	dens.append(cpedsPoint3D(10.5,10.5,0)); 		vals.push_back(1);
	dens.append(cpedsPoint3D(35.5,10.5,0)); 	vals.push_back(1);
	subDomain_region_t r,t;
	r.xmin=0; 		r.xmax=100;
	r.ymin=0; 		r.ymax=100;
	r.zmin=0; 		r.zmax=0;
	r.subx=100;
	r.suby=100;
	r.subz=1;
	t=r;
	t.ymin=-1;
	t.subx=2;
	t.suby=2;
	t.subz=1;

	for (auto& p : dens) {
		for (auto& q: dens) {
			std::cout << p << q << "dist: " << p.dist(q) << std::endl;
		}
	}
	
	std::vector<long> neigh={4};
	for (auto& n : neigh) {
		std::cout << "\n\nTesting interpolation accuracy for Nneigh=" << n << std::endl;
		f.mkInterpolatedFieldScatter(r,t,dens,vals,"gadget2",n,n);
		f.printInfo();

		int i=0;
		for (auto& p : dens) {
			std::cout << "testing interpolation at: " << p << ", " <<
					"should be: " << vals[i] << ", got: " << f.fxy(p.x(),p.y()) << " " <<
					"relerr: " << (f.fxy(p.x(),p.y())-vals[i])/vals[i]*100 << "%" <<
					std::endl;
			i++;
		}
		f.saveSlice(2,0,"xyinter.mat");
	}
	
	return 0;
}
