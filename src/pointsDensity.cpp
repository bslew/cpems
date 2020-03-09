/*!
 * \file pointsDensity.cpp
 *
 *  Created on: Feb 26, 2013
 *      Author: blew
 */

#include "pointsDensity.h"
#include "Mscs-map-window_function.h"
#ifdef HAVE_OPENMP
#include <omp.h>
#endif


pointsDensity::pointsDensity() {
	initiateVars();
}

pointsDensity::pointsDensity(const mscsVector<cpedsPoint3D>& points) : cpedsPointSet3D(points) {
	initiateVars();	
}

pointsDensity::~pointsDensity() {
	if (_pointsDensity_data.tree!=NULL) delete _pointsDensity_data.tree;
}

void pointsDensity::initiateVars() {
	_pointsDensity_data.tree=NULL;
	_pointsDensity_data.hsmlMax=0;
	
}
void pointsDensity::calculateDensity(long NeighborsMin, long NeighborsMax, bool is2dcase, string smKernel, subDomain_region_t *treeScheme, double MassMin, double MassMax) {	
	subDomain_region_t treeSchemeLoc;
	if (treeScheme==NULL) {
		treeSchemeLoc.subx=2;
		treeSchemeLoc.suby=2;
		treeSchemeLoc.subz=1;
		getRanges(treeSchemeLoc.xmin,treeSchemeLoc.xmax,treeSchemeLoc.ymin,treeSchemeLoc.ymax,treeSchemeLoc.zmin,treeSchemeLoc.zmax);
	}
	else
		treeSchemeLoc=*treeScheme;	
	
	double norm;
	double (*kernel)(double );
	if (smKernel == "gadget2") { 
		kernel=&mscsWindowFunction::kernelGadget; 	
		if (is2dcase) norm=double(40.0)/(7.0*PI); //2d case
//		if (r.subz==1) norm=double(40.0)/(7.0*PI); //2d case
		else norm=8.0/PI; //3d case
	}
	else {
		if (smKernel == "gadget2b") { 
			kernel=&mscsWindowFunction::kernelGadget2b; 	
			if (is2dcase) norm=double(40.0)/(7.0*PI); //2d case
			else norm=8.0/PI; //3d case
		}
		else {
			msgs->criticalError("mscsFunction3dregc::mkDensityField >> don't know this smoothing kernel function: "+smKernel,High);
		}
	}
	
	// prepare tree
	msgs->say("building tree with minimal number of particles in subdomain: %li",NeighborsMin,Medium);
//	msgs->say("building tree",Medium);
//	vector<cpedsPoint3D> points=exportAsVector();
	_pointsDensity_data.tree = new subDomain(this,&treeSchemeLoc,NeighborsMin); 
	subDomain *D=_pointsDensity_data.tree;
	D->tree();
#ifdef DEBUG_DENSITY
	D->saveDomains("domains");	
#endif
	msgs->say("done",Low);


	
	// prepare initial smoothing length guess
//	subDomain_region_t particle_reg;
	double hsml=0, rho=0,dist=0;
	
//	_pointsDensity_data.hsml.reserve(size());
	
	mscsVector<long> neighbors;
	long acc=13;
	double dx,dy,dz;
	dx=treeSchemeLoc.xmax-treeSchemeLoc.xmin;
	dy=treeSchemeLoc.ymax-treeSchemeLoc.ymin;
	dz=treeSchemeLoc.zmax-treeSchemeLoc.zmin;
	
//	if (treeSchemeLoc.zmax!=treeSchemeLoc.zmin) { // 3D case
//		hsml=sqrt(dx*dx+dy*dy+dz*dz);
//	}
//	else { // 2D case
//		hsml=sqrt(dx*dx+dy*dy);
//	}
	if (is2dcase) { // 2D case
		hsml=sqrt(dx*dx+dy*dy);
	}
	else { // 3D case
		hsml=sqrt(dx*dx+dy*dy+dz*dz);
	}
	double dr2d=sqrt(dx*dx+dy*dy);
	double dr3d=sqrt(dx*dx+dy*dy+dz*dz);
	double dr;
	if (is2dcase) dr=dr2d; else dr=dr3d;
#ifdef DEBUG_DENSITY
	printf("hsml guess: %lf\n",hsml);
#endif
	//
	// calculate smoothing lengths for all particles (if they were not pre-set)
	//
	long i;
	long N=size();
	if (sml().size()==0) {
		msgs->say("calculating smoothing lengths",Medium);
		_pointsDensity_data.hsml.setSize(size());
#ifdef DEBUG_DENSITY
		printf("looking for neighbors around points \n");
#endif
		
		vector< vector<double> > v;
		v=D->getBranchParticlesCount(true);
#ifdef DEBUG_DENSITY
		printf("number of leaf domains: %li\n", v.size());
#endif		
		long k;
		
#ifdef DEBUG_POINTS_DENSITY2
		long pdone=0;
#endif
		//	long foundNmin, foundNmax;
		mscsVector<long> foundNeigh;
//#ifdef HAVE_OPENMP
//		omp_set_num_threads(omp_get_num_procs());
//#endif
//#pragma omp parallel for schedule (guided) firstprivate(D,hsml,acc,NeighborsMin,NeighborsMax,v,foundNeigh)  private(k,i) reduction(+:pdone)
#ifdef DEBUG_POINTS_DENSITY2
#pragma omp parallel for firstprivate(D,hsml,acc,NeighborsMin,NeighborsMax,v,foundNeigh)  private(k,i) \
	reduction(+:pdone)
#else
#pragma omp parallel for schedule(guided) shared(v,D) firstprivate(hsml,acc,NeighborsMin,NeighborsMax)  private(k,i,foundNeigh) default(shared)
#endif
		for (k = 0; k < v.size(); k++) {
#ifdef DEBUG_POINTS_DENSITY2
			pdone=pdone+1;
#ifdef HAVE_OPENMP
			if (omp_get_thread_num()==0 and pdone%100==0) printf("done: %lf\%\n",double(pdone)/v.size()*omp_get_num_threads()*100);
#else
			if (pdone%100==0) printf("done: %lf\%\n",double(pdone)/v.size()*100);
#endif
#endif
			//		printf("%lf %li\n",v[k][6],v[k].size());
			for (long j = 8; j < v[k].size(); j++) {
				i=v[k][j];
				
				foundNeigh=D->getNeighbors(at(i),NeighborsMin, NeighborsMax, &hsml,acc);
#ifdef DEBUG_POINTS_DENSITY
				if (foundNeigh.size()<NeighborsMin or foundNeigh.size()>NeighborsMax) printf("problem in finding neighbors for particle: %li shold be [%li,%li], found: %li, hsml: %lE\n",i,NeighborsMin, NeighborsMax,foundNeigh.size(),hsml);
#endif
				//			printf("hsml: %lE\n",hsml);
				if (hsml==0) hsml=10*dr; // this is only for NeighborsMin=1 case. TODO: This can be improved 
#pragma omp critical
				sml(i)=hsml;
				//			printf("hsml(i=%li): %lE\n",i,hsml);
				//			if (i==0) {  
				//				at(i).print_point("point 0");
				//				foundNeigh.printVector();
				//				for (long m = 0; m < foundNeigh.size(); m++) {
				//					at(foundNeigh[m]).print_point();
				//				}
				//				exit(0); 
				//			}
			}
		}
	}
	else {
		msgs->say("using pre-set smoothing lengths (hsml vector size: %li)",long(sml().size()),Medium);
	}
//#define DEBUG_POINTS_DENSITY
#ifdef DEBUG_POINTS_DENSITY
	printf("check\n");
	for (long i = 0; i < size(); i++) {
		printf("x: %lf y: %lf z: %lf,  hsml: %lf\n",at(i).x(),at(i).y(),at(i).z(),sml(i));
	}
#endif	
	
	
	// find the biggest hsml for density computations
	double hsmlMax=0;
	double hsmlMin=0;
	hsmlMax=sml(0);
	hsmlMin=sml(0);
	for (long i = 0; i < N; i++) {	if (sml(i)>hsmlMax) hsmlMax=sml(i);	if (sml(i)<hsmlMin) hsmlMin=sml(i);}
	msgs->say("hsmlMin: %lE, hsmlMax: %lE, hsmlMax/boxSize: %lE\n",hsmlMin,hsmlMax,hsmlMax/(dx),Low);
	msgs->say("done",Low);
	_pointsDensity_data.hsmlMax=hsmlMax;
	
	//
	// do the density calculation
	//
	msgs->say("calculating density",Medium);
	double x=0,y=0,z=0;
//	double hsmlo2=hsml/2;

	
	density().setSize(size());
//#ifdef HAVE_OPENMP
//	omp_set_num_threads(omp_get_num_procs());
//#endif
#pragma omp parallel for schedule (guided) firstprivate(D,hsml,x,y,z,hsmlMax, neighbors,rho,dist,norm)
//#pragma omp parallel for firstprivate(D,hsml,x,y,z,hsmlMax, neighbors,rho,dist,norm)
	for (long i = 0; i < size(); i++) {
#ifdef DEBUG_POINTS_DENSITY2
#ifdef HAVE_OPENMP
		if (omp_get_thread_num()==0 and (i%1000==0) ) 
			msgs->say("done: "+msgs->toStr(double(i)*omp_get_num_threads()/size()*100)+" %",Low);
#else
		if (i%1000==0) printf("done: %lf\%\n",double(i)/size()*100);
#endif
#endif
		x=at(i).x();
		y=at(i).y();
		z=at(i).z();

		// get particles within the distance of R defining region containing constant number of particles
		neighbors=D->getPointsIdx(cpedsPoint3D(x,y,z),hsmlMax,NULL);
		if (hsml==0) { msgs->criticalError("mkDensityField>>> hsml is 0. Cannot continue",Top); }
		rho=0;
		// calculate density
		if (is2dcase) {
			for (long l = 0; l < neighbors.size(); l++) {
				dist=cpedsPoint3D(x,y,z).dist(at(neighbors[l]));
				rho+=kernel(dist/(sml(neighbors[l])))/(sml(neighbors[l])*sml(neighbors[l]));
			}
		}
		else {
			for (long l = 0; l < neighbors.size(); l++) {
				dist=cpedsPoint3D(x,y,z).dist(at(neighbors[l]));
				rho+=kernel(dist/(sml(neighbors[l])))/(sml(neighbors[l])*sml(neighbors[l])*sml(neighbors[l]));
			}
		}

		density().at(i)=norm*rho;
	}
	msgs->say("calculating density done",Low);
#ifdef DEBUG_POINTS_DENSITY
	printf("check\n");
	for (long i = 0; i < size(); i++) {
		printf("x: %lf y: %lf z: %lf,  hsml: %lf, rho: %lf\n",at(i).x(),at(i).y(),at(i).z(),sml(i), density()[i]);
	}
#endif

	
}
