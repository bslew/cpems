/*!
 * \file test-domain-decomposition.cpp
 *
 *  Created on: Jan 27, 2012
 *      Author: blew
 */

#include <vector>
#include "cpeds-point3d.h"
#include "cpeds-point_set.h"
#include "subdomain.h"
#include "cpeds-math.h"
#include "mscsVector.h"
#include "Mscs-function3dregc.h"

#define domainMIN 0.0
#define domainMAX 10.0

int main(int argc, char** argv) {
//	vector<cpedsPoint3D> v;
	mscsVector<cpedsPoint3D> v;
	cpedsMsgs msgs("test");
	msgs.setTimeSinceStartOn(true);
//	mscsVector<double> vd;
//	vd.setSize(10);
//	vd=10;
//	vd.printVector();

//	cpedsPoint2D p(2,2);
//	cout << p;
//	exit(0);

	msgs.say("generting points",High);
	long Npoints=0;
	long NminPart=25;
	if (argc==1) { 
		printf("\n *** NOTE: if you want to generate N randon points then the usage is: test-domain-decomposition N\n"
				"if you want to specify the maximal number of particles in cell then use: test-domain-decomposition N Nmax\n\n");
		v.push_back(cpedsPoint3D(1,1,5));
		v.push_back(cpedsPoint3D(2,2,5));
		v.push_back(cpedsPoint3D(7,1,5));
		v.push_back(cpedsPoint3D(7,8,5));
		v.push_back(cpedsPoint3D(1.1,1.1,5));
		v.push_back(cpedsPoint3D(1.1,1.2,5));
		v.push_back(cpedsPoint3D(1.1,1.3,5));
		v.push_back(cpedsPoint3D(1.1,1.4,5));
		v.push_back(cpedsPoint3D(1.1,1.5,5));

//		v.push_back(cpedsPoint3D(3,5,5));
//		v.push_back(cpedsPoint3D(7,5,5));
//		v.push_back(cpedsPoint3D(5,5,5));
		Npoints=v.size();
	}
	else {
		Npoints=strtol(argv[1],NULL,10);
		if (argc>2) {
			NminPart=strtol(argv[2],NULL,10);
		}
		srand(time(0));
		double *x=cpeds_random_uniform_numbers(domainMIN,domainMAX,Npoints);
		double *y=cpeds_random_uniform_numbers(domainMIN,domainMAX,Npoints);
		double *z=cpeds_random_uniform_numbers(domainMIN,domainMAX,Npoints);
		for (long i = 0; i < Npoints; i++) {
			v.push_back(cpedsPoint3D(x[i],y[i],z[i]));
		}
		delete [] x; delete [] y; delete [] z;
	}
	cpedsPointSet3D ps(v);	
	if (Npoints<10000) {
		msgs.say("saving points",High);
		ps.save("pointSet");
	}
	
	subDomain_region_t r,subr;
	r.subx=2;
	r.suby=2;
	r.subz=1;
	r.xmin=domainMIN;
	r.xmax=domainMAX;	
	r.ymin=domainMIN;
	r.ymax=domainMAX;
	r.zmin=domainMIN;
	r.zmax=domainMAX;

	subDomain D(&v,&r,NminPart);
	
	msgs.say("doing tree",High);
	D.tree();
	
//	exit(0);
	if (Npoints<10000) {
		msgs.say("saving domains",High);
		D.saveDomains("domains");
	}

	vector<long>* ptr;
	if (Npoints<10000) {
		printf("\n\nsearching for particles in region:\n\n");
	//	subr.xmin=2;
	//	subr.xmax=4;
	//	subr.ymin=0.1;
	//	subr.ymax=1.1;
	//	subr.zmin=0.0;
	//	subr.zmax=10;

		subr.xmin=1.5;
		subr.xmax=2.1;
		subr.ymin=1.5;
		subr.ymax=2.1;
		subr.zmin=0.0;
		subr.zmax=10;
		D.print_domain_range(subr);
		subDomain* sD;
		ptr=D.getContainingDomainPointIdx(subr,&sD);
	
		if (ptr==NULL) { printf("\n\nnie ma cząstek w tym regionie\n\n"); }
		else {
			printf("cząstki w domenie zawierającej:\n");
			for (long i = 0; i < ptr->size(); i++) {
				printf("%lf, %lf, %lf\n",v[ptr->at(i)].x(),v[ptr->at(i)].y(),v[ptr->at(i)].z());
			}
			printf("containing subdomain info:\n");
			sD->print_domain_range();
		}
		printf("\n");
	}
	
	msgs.say("FINDING NEIGHBORS TEST 1",Top);
	cpedsPoint3D p(1,2,3);
	long t0,t1;
	t0=msgs.timeElapsed();
	long Nmin=cpeds_get_min(25,Npoints-1);
	long Nmax=cpeds_get_min(25,Npoints-1)+1;
	if (Npoints==v.size()) {
		Nmin=25;
		Nmax=25;
		p.set(1.5,1.5,5.0);
	}
	printf("Nmin: %li, Max: %li\n",Nmin,Nmax);

	mscsVector<long> neigh=D.getNeighbors(p,Nmin,Nmax);
	t1=msgs.timeElapsed();
	neigh.printVector();
	msgs.say("neigh search time: %lf [s]",double(t1-t0),High);
	exit(0);

	{
		printf("\nSEARCHING PARTICLES IN RECTANGULAR REGION:\n\n");
		D.print_domain_range(subr);
		cpedsList<long> l=D.getPointsIdx(subr);
		if (ptr==NULL) { printf("\n\nnie ma cząstek w tym regionie\n\n"); }
		else {
			printf("cząstki w tym regionie:\n");
			for (long i = 0; i < l.size(); i++) {
				printf("%lf, %lf, %lf\n",v[l[i]].x(),v[l[i]].y(),v[l[i]].z());
			}
		}
	}

	{
		double R=3;
		printf("\nSEARCHING PARTICLES IN SPHERICAL REGION OF RADIUS: %lf\n\n",R);
		double tmp=(domainMAX+domainMIN)/2;
		mscsVector<long> l=D.getPointsIdx(cpedsPoint3D(tmp,tmp,tmp),R);
		if (l.size()==0) { printf("\n\nnie ma cząstek w tym regionie\n\n"); }
		else {
			printf("cząstki w tym regionie:\n");
			for (long i = 0; i < l.size(); i++) {
				printf("%lf, %lf, %lf\n",v[l[i]].x(),v[l[i]].y(),v[l[i]].z());
			}
		}
	}
	
//	{
//		long Nmin=3;
////		long Nmax=1;
//		printf("\nSEARCHING PARTICLES IN SPHERICAL REGION in the first sub-domain with Nmin=%li\n\n",Nmin);	
//		double tmp=(domainMAX+domainMIN)/2;
////		cpedsList<long> l=D.getNeighbors(cpedsPoint3D(tmp/2,tmp/2,tmp/2),Nmin,Nmax);
//		cpedsList<long> l=D.getNeighbors(cpedsPoint3D(0,0,5),Nmin);
//		if (l.size()==0) { printf("\n\nnie ma cząstek w tym regionie\n\n"); }
//		else {
//			printf("cząstki w tym regionie:\n");
//			for (long i = 0; i < l.size(); i++) {
//				printf("%lf, %lf, %lf\n",v[l[i]].x(),v[l[i]].y(),v[l[i]].z());
//			}
//		}
//	}
	
	{
		long Nmin=4;
		printf("\nTEST THE DENSITY CALCULATIONS - gather\n\n");	
		mscsFunction3dregc f;
		subr.subx=2;
		subr.suby=2;
		subr.subz=1;
		subr.xmin=domainMIN;
		subr.xmax=domainMAX;	
		subr.ymin=domainMIN;
		subr.ymax=domainMAX;
		subr.zmin=domainMIN;
		subr.zmax=domainMAX;

		r.xmin=domainMIN; r.xmax=domainMAX;
		r.ymin=domainMIN; r.ymax=domainMAX;
		r.zmin=5; r.zmax=5;
		r.subx=100;
		r.suby=100;
		r.subz=1;
		
		f.mkDensityFieldGather(r,subr,v,NULL,0,"gadget2",Nmin);
		f.saveSlice(2,0,"densityGather4",0);
		f.saveSlice(2,0,"hsmlGather4",1);
//		f.saveSlice(2,0,"neighbors",1);
		
	}

//	{
//		long Nmin=4;
//		printf("\nTEST THE DENSITY CALCULATIONS 2\n\n");	
//		mscsFunction3dregc f;
//		r.xmin=domainMIN; r.xmax=domainMAX;
//		r.ymin=domainMIN; r.ymax=domainMAX;
//		r.zmin=5; r.zmax=5;
//		r.subx=100;
//		r.suby=100;
//		r.subz=1;
//		
//		f.mkDensityFieldScatter(r,v,NULL,0,"gadget2",Nmin);
//		f.saveSlice(2,0,"density4",0);
//		f.saveSlice(2,0,"hsml4",1);
////		f.saveSlice(2,0,"neighbors",1);
//		
//	}

	return 0;
}
