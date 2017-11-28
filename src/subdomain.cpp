/*!
 * \file subDomain.cpp
 *
 *  Created on: Jan 26, 2012
 *      Author: blew
 */

#include "subdomain.h"
#include "Mscs-common.h"
#include "QtCore/QList"


/***************************************************************************************/
/* top domain initialization */ 
/***************************************************************************************/

subDomain::subDomain(const vector<cpedsPoint3D>* points, const subDomain_region_t *domain, long minPoints) {
//	subDomain::subDomainDef_t _rootDomain;
	_points=points;
	_def.region=*domain;
	_def.subDomainsCount=domain->subx*domain->suby*domain->subz;
	_def.minObjectsCount=minPoints;
	_def.root=true;
	_def.curDepth=0;
	_def.maxDepth=-1;
	_def.rootDomain=this;
//	printf("%li %li %li\n",_def.region.subx,_def.region.suby,_def.region.subz);
	
	initiate(points);
}
subDomain::subDomain(const mscsVector<cpedsPoint3D>* points, const subDomain_region_t *domain, long minPoints) {
	//	subDomain::subDomainDef_t _rootDomain;
		_points=points;
		_def.region=*domain;
		_def.subDomainsCount=domain->subx*domain->suby*domain->subz;
		_def.minObjectsCount=minPoints;
		_def.root=true;
		_def.curDepth=0;
		_def.maxDepth=-1;
		_def.rootDomain=this;
	//	printf("%li %li %li\n",_def.region.subx,_def.region.suby,_def.region.subz);
		
		initiate(points);	
}

void subDomain::initiate(const vector<cpedsPoint3D>* points) {
	if (_points!=NULL) {
//		printf("allocating space for sub-domain %li pointers\n",_def.subDomainsCount);
//		_def.sub = new subDomain*[_def.subDomainsCount];
//		for (long i = 0; i < _def.subDomainsCount; i++) { _def.sub[i]=NULL; }

//		_def.idx=new vector<long>(points->size());
//		for (long i = 0; i < points->size(); i++) { _def.idx->at(i)=i; } // all points are initially in the domain

		// find all particles in this domain and put their indexes on the idx vector
		_def.idx=new vector<long>;
		for (unsigned long i = 0; i < points->size(); i++) { if (inDomain(points->at(i))) _def.idx->push_back(i); } // find points in the domain, not all of them might be in this domain
		if (_def.idx->size()==0) { printf("created subDomain contains no points in it. This is probably not what you wanted. please check the ranges.\n"); exit(0); }
	}
}
/***************************************************************************************/
subDomain::subDomain(const subDomain& rhs) {
	_points=rhs._points;
	_def=rhs.getSubdomainInfo();
//	printf("_def.subDomainsCount: %li\n",_def.subDomainsCount);
	_def.idx=NULL;
	if (rhs.getSubdomainInfo().idx!=NULL) _def.idx = new vector<long>(*(rhs.getSubdomainInfo().idx));
	_def.rootDomain=this;
//	printf("rhs sd count: %li:\n", rhs._def.subDomainsCount);
	if (rhs.getSubdomainInfo().subDomainsCount!=0) { _def.sub = new subDomain*[rhs.getSubdomainInfo().subDomainsCount]; }
	
	for (long i = 0; i < _def.subDomainsCount; i++) {
//		printf("rhs.getSubdomainInfo().sub[i]: %li\n", rhs.getSubdomainInfo().sub[i]);
		_def.sub[i]=NULL;
//		printf("rhs.getSubdomainInfo().sub[i]: %li\n", rhs.getSubdomainInfo().sub[i]);
//		printf("subdomain pass %li\n",i); 
		if (rhs.getSubdomainInfo().sub[i]!=NULL) { 
//			printf("initiating subdomain %li\n",i); 
		_def.sub[i] = new subDomain(*(rhs.getSubdomainInfo().sub[i]), this); }
	}
}
subDomain::subDomain(const subDomain& rhs, subDomain* rootP) {
	_points=rhs._points;
	_def=rhs.getSubdomainInfo();
	_def.idx=NULL;
	if (rhs.getSubdomainInfo().idx!=NULL) _def.idx = new vector<long>(*(rhs.getSubdomainInfo().idx));
	_def.rootDomain=rootP;	
	if (rhs.getSubdomainInfo().subDomainsCount!=0) { _def.sub = new subDomain*[rhs.getSubdomainInfo().subDomainsCount]; }

	for (long i = 0; i < _def.subDomainsCount; i++) {
		_def.sub[i]=NULL;
		if (rhs.getSubdomainInfo().sub[i]!=NULL) { _def.sub[i] = new subDomain(*(rhs.getSubdomainInfo().sub[i]), this); }
	}
}
/***************************************************************************************/
subDomain::~subDomain() {
	if (_def.sub!=NULL) {
		// delete sub-domains
		for (long i = 0; i < _def.subDomainsCount; i++) {
			if (_def.sub[i]!=NULL) { delete _def.sub[i]; }
		}
		delete [] _def.sub; _def.sub=NULL;
	}
	// delete list of particles in this domain
	if (_def.idx!=NULL) { delete _def.idx; _def.idx=NULL; }

	// clean up tree walk stuff
	subDomain_region_t reg;
	getOverlappingDomainPointIdx(reg,NULL,false,true);
}

/***************************************************************************************/
/* sub domain initialization for tree building */
/***************************************************************************************/
subDomain::subDomain(subDomain* parent, const subDomain_region_t& newSubDomain) {
	_def=parent->_def;
//	printf("%li %li %li\n",_def.region.subx,_def.region.suby,_def.region.subz);
	_def.region=newSubDomain; 
	_def.root=false;
	_def.leaf=true;
	parent->_def.leaf=false;
	_def.rootDomain=parent->_def.rootDomain;
	_points=parent->_points;
	_def.maxDepth=parent->_def.maxDepth;
	_def.curDepth=parent->_def.curDepth+1;
	
//	printf("allocating space for sub-domain %li pointers\n",_def.subDomainsCount);
	_def.sub=NULL;
	//	_def.sub = new subDomain*[_def.subDomainsCount];
//	for (long i = 0; i < _def.subDomainsCount; i++) { _def.sub[i]=NULL; }

	// find all particles in this domain and put their indexes on the idx vector
	_def.idx=new vector<long>;	
//	if (_def.idx->size()>_def.minObjectsCount) {
	findPointsInDomain(parent);
	if (_def.maxDepth>0) {
		if (_def.curDepth<_def.maxDepth) tree();
		else {
			/*
			 * Comment: this is a fix for tree initialization
			 * It always should have subdomains initialized even if all are empty.
			 * This can be re-witten a different way someday but for the moment
			 * in many places the code iterates over the subdomains checking if 
			 * they are NULL. So if this structure should be initialized
			 * 
			 * author: blew
			 * date: Nov 27, 2017 4:02:20 PM
			 *
			 */
			_def.sub = new subDomain*[_def.subDomainsCount];
			for (long i = 0; i < _def.subDomainsCount; i++) { _def.sub[i]=NULL; }			
		}
	}
	else {
		tree();
	}
		
//	}
}

/***************************************************************************************/
//void subDomain::tree(const subDomain *domain) {
void subDomain::tree() {
#ifdef DEBUG_SUBDOMAIN
	printf("entering tree()\n");
#endif
	_def.sub = new subDomain*[_def.subDomainsCount];
	for (long i = 0; i < _def.subDomainsCount; i++) { _def.sub[i]=NULL; }

	//	if (domain==NULL) domain=this;
	subDomain_region_t r;
//	double Dx=(domain->_def.region->xmax-domain->_def.region->xmin)/domain->_def.region->subx;
//	double Dy=(domain->_def.region->ymax-domain->_def.region->ymin)/domain->_def.region->suby;
//	double Dz=(domain->_def.region->zmax-domain->_def.region->zmin)/domain->_def.region->subz;

	if (_def.idx->size()>_def.minObjectsCount) { // if the number of points in this domain is larger than the minimal number of points in the domain then do another decomposition
#ifdef DEBUG_SUBDOMAIN
		printf("- number of particles in this domain is > min, decomposing\n");
		print_domain_range(_def.region);
#endif
		double Dx=(_def.region.xmax-_def.region.xmin)/_def.region.subx;
		double Dy=(_def.region.ymax-_def.region.ymin)/_def.region.suby;
		double Dz=(_def.region.zmax-_def.region.zmin)/_def.region.subz;
#ifdef DEBUG_SUBDOMAIN
		printf("- subdomain size: dx=%lf, dy= %lf, dz=%lf\n",Dx,Dy,Dz);
#endif
	
		
		int l=0;
		r.subx=_def.region.subx;
		r.suby=_def.region.suby;
		r.subz=_def.region.subz;
		for (long i = 0; i < _def.region.subx; i++) {
			//		r.xmin=domain->_def.region->xmin + i*Dx;  r.xmax=r.xmin + (i+1)*Dx;
			r.xmin=_def.region.xmin + i*Dx;  
			r.xmax=_def.region.xmin + (i+1)*Dx;   // this better prevents creating subdomains with tiny gaps in space in between subdomain boundaries
			if (i==_def.region.subx-1) r.xmax=_def.region.xmax;
//			r.xmax=r.xmin + Dx;
			for (long j = 0; j < _def.region.suby; j++) {
				//			r.ymin=domain->_def.region->ymin + j*Dy;  r.ymax=r.ymin + (j+1)*Dy;
				r.ymin=_def.region.ymin + j*Dy;
				r.ymax=_def.region.ymin + (j+1)*Dy;
				if (j==_def.region.suby-1) r.ymax=_def.region.ymax;
//				r.ymax=r.ymin + Dy;
				for (long k = 0; k < _def.region.subz; k++) {
								//				r.zmin=domain->_def.region->zmin + k*Dz;  r.zmax=r.zmin + (k+1)*Dz;
								//				_def.sub[l] = new subDomain(domain,r);

					// define a new sub-domain region
					r.zmin=_def.region.zmin + k*Dz;
					r.zmax=_def.region.zmin + (k+1)*Dz;
					if (k==_def.region.subz-1) r.zmax=_def.region.zmax;
//					r.zmax=r.zmin + Dz;
					
//					print_domain_range(r);
					// create a new sub-domain
					_def.sub[l] = new subDomain(this,r);
					if (_def.sub[l]->_def.idx->size()==0) { delete _def.sub[l]; _def.sub[l]=NULL; 
//						delete _def.idx; _def.idx=NULL; 
#ifdef DEBUG_SUBDOMAIN
					printf("deleting subdomain which does not contain any points in it\n");
#endif
					}
					l++;
				}
			}
		}
		
		
	}
	else { 
#ifdef DEBUG_SUBDOMAIN
		printf("- number of particles in this domain is <= min, NOT decomposing\n");	
#endif
		_def.leaf=true;
	}
	
#ifdef DEBUG_SUBDOMAIN
	printf("leaving tree()\n");
#endif

}

/***************************************************************************************/
void subDomain::findPointsInDomain(const subDomain* parent) {
#ifdef DEBUG_SUBDOMAIN
	printf("finding points in domain\n");
	print_domain_range();	
#endif
//	long n=parent->_points->size();
	long n=parent->_def.idx->size();
	for (long i = 0; i < n; i++) {
		if (inDomain(parent->_points->at(parent->_def.idx->at(i)))) { 
			_def.idx->push_back(parent->_def.idx->at(i)); 		
			
#ifdef DEBUG_SUBDOMAIN
			printf("point %li (%lf, %lf,%lf) is in domain\n",i, parent->_points->at(parent->_def.idx->at(i)).x(),parent->_points->at(parent->_def.idx->at(i)).y(),parent->_points->at(parent->_def.idx->at(i)).z());
#endif

		}
	}
#ifdef DEBUG_SUBDOMAIN
	printf("found: %li points \n",_def.idx->size());
#endif
}
/***************************************************************************************/
bool subDomain::inDomain(const cpedsPoint3D& p) const {
//	if (p.x()>=_def.region.xmin and p.x()<=_def.region.xmax and 
//		p.y()>=_def.region.ymin and p.y()<=_def.region.ymax and 
//		p.z()>=_def.region.zmin and p.z()<=_def.region.zmax ) return true;

	if (p.x()>=_def.region.xmin and p.x()<_def.region.xmax and 
		p.y()>=_def.region.ymin and p.y()<_def.region.ymax and 
		p.z()>=_def.region.zmin and p.z()<_def.region.zmax ) return true;

	//
	// check along the subdomain faces
	//
	
	// x coordinate is on the outermost x coordinate of the root domain
	if (p.x()==_def.region.xmax and p.x()==_def.rootDomain->getSubdomainInfo().region.xmax and 
		p.y()>=_def.region.ymin and p.y()<_def.region.ymax and 
		p.z()>=_def.region.zmin and p.z()<_def.region.zmax ) return true;
	
	// y coordinate is on the outermost y coordinate of the root domain
	if (p.x()>=_def.region.xmin and p.x()<_def.region.xmax and 
		p.y()==_def.region.ymax and p.y()==_def.rootDomain->getSubdomainInfo().region.ymax and 
		p.z()>=_def.region.zmin and p.z()<_def.region.zmax ) return true;


	// z coordinate is on the outermost z coordinate of the root domain
	if (p.x()>=_def.region.xmin and p.x()<_def.region.xmax and 
		p.y()>=_def.region.ymin and p.y()<_def.region.ymax and 
		p.z()==_def.region.zmax and p.z()==_def.rootDomain->getSubdomainInfo().region.zmax ) return true;

	//
	// check along the subdomain lines
	//

	// x,y coordinate is on the outermost x and y coordinate of the root domain
	if (p.x()==_def.region.xmax and p.x()==_def.rootDomain->getSubdomainInfo().region.xmax and 
		p.y()==_def.region.ymax and p.y()==_def.rootDomain->getSubdomainInfo().region.ymax and 
		p.z()>=_def.region.zmin and p.z()<_def.region.zmax ) return true;
	// x,z coordinate is on the outermost x and z coordinate of the root domain
	if (p.x()==_def.region.xmax and p.x()==_def.rootDomain->getSubdomainInfo().region.xmax and 
		p.y()>=_def.region.ymin and p.y()<_def.region.ymax and 
		p.z()==_def.region.zmax and p.z()==_def.rootDomain->getSubdomainInfo().region.zmax ) return true;
	// y,z coordinate is on the outermost y and z coordinate of the root domain
	if (p.x()>=_def.region.xmin and p.x()<_def.region.xmax and 
		p.y()==_def.region.ymax and p.y()==_def.rootDomain->getSubdomainInfo().region.ymax and 
		p.z()==_def.region.zmax and p.z()==_def.rootDomain->getSubdomainInfo().region.zmax ) return true;

	//
	// check the far corner of the subdomain
	//
	if (p.x()==_def.region.xmax and p.x()==_def.rootDomain->getSubdomainInfo().region.xmax and 
		p.y()==_def.region.ymax and p.y()==_def.rootDomain->getSubdomainInfo().region.ymax and 
		p.z()==_def.region.zmax and p.z()==_def.rootDomain->getSubdomainInfo().region.zmax ) return true;

#ifdef DEBUG
	if (_def.root)	{
		printf("\nWARNING: ");
		p.print_point("lost point in root domain");
	}
#endif
	return false;
}
/***************************************************************************************/
bool subDomain::inDomain(const cpedsPoint3D& p, double rx,double ry,double rz) const {
	cpedsPoint3D pxlow(p), pxhigh(p), pylow(p), pyhigh(p), pzlow(p), pzhigh(p);
	pxlow.x()-=rx;	pxhigh.x()+=rx;
	pylow.y()-=ry;	pyhigh.y()+=ry;
	pzlow.z()-=rz;	pzhigh.z()+=rz;
//	printf("pxlow: %i pxhigh: %i\n",inDomain(pxlow) ,inDomain(pxhigh));
//	printf("pylow: %i pyhigh: %i\n",inDomain(pylow) ,inDomain(pyhigh));
//	printf("pzlow: %i pzhigh: %i\n",inDomain(pzlow) ,inDomain(pzhigh));
	if (inDomain(pxlow) and inDomain(pxhigh) and
		inDomain(pylow) and inDomain(pyhigh) and
		inDomain(pzlow) and inDomain(pzhigh) ) return true;
	
	return false;
}
/***************************************************************************************/
bool subDomain::partiallyInDomain(const subDomain_region_t& r, const subDomain* D) const {
	
	if (D!=NULL) {
		return D->partiallyInDomain(r,NULL);
	}
	if (r.xmin>_def.region.xmax ) return false;
	if (r.ymin>_def.region.ymax ) return false;
	if (r.zmin>_def.region.zmax ) return false;

	if (r.xmax<_def.region.xmin) return false;
	if (r.ymax<_def.region.ymin) return false;
	if (r.zmax<_def.region.zmin) return false;

	return true;	
}

/***************************************************************************************/
bool subDomain::inDomain(const subDomain_region_t& r, const subDomain* D) const {

	// NEW recursive VERSION: Mar 29, 2012, 1:17:14 PM
	
	if (D!=NULL) {
		return D->inDomain(r,NULL);
	}
	
	if (r.xmin>=_def.region.xmin and r.xmax<=_def.region.xmax and
		r.ymin>=_def.region.ymin and r.ymax<=_def.region.ymax and
		r.zmin>=_def.region.zmin and r.zmax<=_def.region.zmax) return true;
	return false;
	
	
	/***************************************************************************************/
//	cpedsPoint3D p((r.xmax+r.xmin)/2, (r.ymax+r.ymin)/2, (r.zmax+r.zmin)/2);
//
//	if (D!=NULL) {
//		return D->inDomain(p,(r.xmax-r.xmin)/2, (r.ymax-r.ymin)/2, (r.zmax-r.zmin)/2);
//	}
//
////	printf("looking in\n");
////	p.print_point();
////	printf("D: %lf %lf %lf\n",(r.xmax-r.xmin)/2, (r.ymax-r.ymin)/2, (r.zmax-r.zmin)/2);
//	return inDomain(p,(r.xmax-r.xmin)/2, (r.ymax-r.ymin)/2, (r.zmax-r.zmin)/2);
}
/***************************************************************************************/
bool subDomain::inRegion(const subDomain_region_t& r, const cpedsPoint3D& p) {
//	p.print_point("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$");
//	print_domain_range(r);
	if (p.x()>=r.xmin and p.x()<=r.xmax and 
		p.y()>=r.ymin and p.y()<=r.ymax and 
		p.z()>=r.zmin and p.z()<=r.zmax ) return true;
	return false;	
}
/***************************************************************************************/
bool subDomain::inRange(double r,const cpedsPoint3D& p0, const cpedsPoint3D& p) {
	if (p.dist(p0) <=r) return true;
	return false;		
}
/***************************************************************************************/
vector< vector<long>* >* subDomain::getOverlappingDomainPointIdx(const subDomain_region_t& r, subDomain** sd, bool newSearch, bool claenUp) {
	static vector< vector<long>* >* Vptr=NULL; // the value of this variable is shared between all objects that build the tree
//	static vector< vector<long>* >*  ptr=NULL; // the value of this variable is shared between all objects that build the tree
//	static subDomain* containingDomainPtr=NULL; // the value of this variable is shared between all objects that build the tree
	static int treeLevel=0; // 
//	static int containingDomainLevel=0; // 
//#pragma omp threadprivate(Vptr,containingDomainPtr, treeLevel, containingDomainLevel)
#pragma omp threadprivate(Vptr, treeLevel)
	//	vector<long>* ptr=NULL;
	//	subDomain* containingDomainPtr=NULL;
	
	if (claenUp) {
		if (Vptr!=NULL) {
			delete Vptr;
			Vptr=NULL;
		}
		return NULL;
	}
	if (newSearch) {
//		ptr=&Vptr;
		if (Vptr==NULL) {
			Vptr=new vector< vector<long>* >;
//			printf("allocating new vector\n");
		
		}
		else Vptr->clear();
//		containingDomainPtr=NULL;		
		treeLevel=0;
//		containingDomainLevel=0;
	}
	//	printf("ptr BEGIN: %li\n",ptr);
	treeLevel++;
	
		
		//	print_domain_range(_def.region); 
#ifdef DEBUG_SUBDOMAIN2
		printf("\n\n\ngetOverlappingDomainPointIdx, region is\n");
		print_domain_range(r);
		printf("tree level: %i\n",treeLevel);
		printf("checking new domain\n");
		print_domain_range();
#endif
		
		if (partiallyInDomain(r)) { 
			if (_def.leaf) {
				Vptr->push_back(_def.idx);
#ifdef DEBUG_SUBDOMAIN2
				printf("-------------------found leaf domain:\n"); 
				printf("appending indexes\n");
				printf("number of domains containing the region: %li\n",Vptr->size());
				print_domain_range(_def.region); 
#endif
			}
//			ptr=_def.idx;  containingDomainPtr=this;
#ifdef DEBUG_SUBDOMAIN2
			printf("!!!!!!!!!!!!!!!!!!!!!found domain, but it's not a leaf:\n"); 
#endif
//			containingDomainLevel=treeLevel;
		}
		else {
#ifdef DEBUG_SUBDOMAIN2
			printf("now lev: %i\n",treeLevel);
			printf("going level up\n");
#endif
			treeLevel--;
//			if (sd!=NULL) *sd=containingDomainPtr; // TODO this can be improved somehow since eg for 1000 levels when going up this assignment will be done 1000 times with the same value.

			return Vptr;
		}
#ifdef DEBUG_SUBDOMAIN2
		printf("now lev: %i\n",treeLevel);
		printf("checking subdomains\n");
#endif

		for (long i = 0; i < _def.subDomainsCount; i++) { // TODO: this can be improved if after one pass a domain is found then there is no point in further search at the same level because the domains are non-overlapping
#ifdef DEBUG_SUBDOMAIN2
			printf("checking subdomain: %i\n",i);
#endif
//			if (treeLevel>=containingDomainLevel) {
				if (_def.sub[i]!=NULL) {
					_def.sub[i]->getOverlappingDomainPointIdx(r);
			//			printf("ptr: %li\n",ptr);
				}
#ifdef DEBUG_SUBDOMAIN2
				else {
					printf("    no more subdomains\n");
				}
#endif
//			}
//#ifdef DEBUG_SUBDOMAIN2
//			else {
//				printf("    wrong branch\n");				
//			}
//#endif
		}
//		if (sd!=NULL) *sd=containingDomainPtr; // TODO this can be improved somehow since eg for 1000 levels when going up this assignment will be done 1000 times with the same value.
		treeLevel--;
#ifdef DEBUG_SUBDOMAIN2
		printf("now lev: %i\n",treeLevel);
		printf("going level up\n");
#endif

	return Vptr;
	
}
/***************************************************************************************/
//vector<long>* subDomain::getContainingDomainPointIdx(const subDomain_region_t& r, subDomain** sd, bool ignoreRegionSize) {
vector<long>* subDomain::getContainingDomainPointIdx(const subDomain_region_t& r, subDomain** sd, bool newSearch) {
	static vector<long>* ptr=NULL; // the value of this variable is shared between all objects that build the tree
	static subDomain* containingDomainPtr=NULL; // the value of this variable is shared between all objects that build the tree
	static int treeLevel=0; // 
	static int containingDomainLevel=0; // 
#pragma omp threadprivate(ptr,containingDomainPtr, treeLevel, containingDomainLevel)
	//	vector<long>* ptr=NULL;
	//	subDomain* containingDomainPtr=NULL;
	if (newSearch) {
		ptr=NULL;
		containingDomainPtr=NULL;		
		treeLevel=0;
		containingDomainLevel=0;
	}
	//	printf("ptr BEGIN: %li\n",ptr);
	treeLevel++;
	

	long NnonNULLsub=0, NNULLsub=0;
	for (int i=0;i< _def.subDomainsCount;i++) {
		if (_def.sub!=NULL) {
			if (_def.sub[i]!=NULL) NnonNULLsub++;
			else NNULLsub++;
		}
	}
	if (NnonNULLsub==0 and NNULLsub==0) { printf("NnonNULLsub=0 and NNULLsub=0. This should not happen. "
			"This is most likely because _def.sub==NULL i.e. is not initialized at the bottom of the tree. FIX IT !\n"); exit(1); }

		//	print_domain_range(_def.region); 
#ifdef DEBUG_SUBDOMAIN
		printf("\n\n\ngetContainingDomainPointIdx, region is\n");
		print_domain_range(r);
		printf("\n\n going down: containing domain PTR: %li\n",containingDomainPtr);
		printf("lev: %i, CDL: %i\n",treeLevel,containingDomainLevel);
		printf("checking new domain\n");
		print_domain_range();
		printf("this domain has %li non-NULL subdomains and %li NULL subdomains\n-----\n\n",NnonNULLsub, NNULLsub);
//		if (treeLevel>_def.maxDepth) return ptr;
#endif
		
		if (inDomain(r)) { 
			ptr=_def.idx;  containingDomainPtr=this;
#ifdef DEBUG_SUBDOMAIN
			printf("!!!!!!!!!!!!!!!!!!!!!found domain:\n"); 
			print_domain_range(_def.region); 
#endif
			containingDomainLevel=treeLevel;
		}
		else {
#ifdef DEBUG_SUBDOMAIN
		printf("going up (containing domain PTR: %li\n\n\n",containingDomainPtr);		
#endif
			treeLevel--;
			if (sd!=NULL) *sd=containingDomainPtr; // TODO this can be improved somehow since eg for 1000 levels when going up this assignment will be done 1000 times with the same value.

			return ptr;
		}
#ifdef DEBUG_SUBDOMAIN
		printf("now containing domain PTR: %li\n",containingDomainPtr);
		printf("now lev: %i, CDL: %i\n",treeLevel,containingDomainLevel);
		printf("checking subdomains\n");
#endif

		for (long i = 0; i < _def.subDomainsCount; i++) { // TODO: this can be improved if after one pass a domain is found then there is no point in further search at the same level because the domains are non-overlapping
#ifdef DEBUG_SUBDOMAIN
			printf("checking subdomain: %i/%i\n",i+1,_def.subDomainsCount);
#endif
			if (treeLevel>=containingDomainLevel) {
				if (_def.sub[i]!=NULL) {
					//			_def.sub[i]->getContainingDomainPointIdx(r,l);
					ptr=_def.sub[i]->getContainingDomainPointIdx(r,sd);
			//			printf("ptr: %li\n",ptr);
				}
#ifdef DEBUG_SUBDOMAIN
				else {
					printf("    no more subdomains\n");
				}
#endif
			}
#ifdef DEBUG_SUBDOMAIN
			else {
				printf("    wrong branch\n");				
			}
#endif
		}
		//	if (l==1) { // check if the particles are indeed in side of the requested region
		//		if (ptr!=NULL) {
		//			for (long i = 0; i < ptr->size(); i++) {
		//			}
		//		}
		//	}
		//	l--;
		if (sd!=NULL) *sd=containingDomainPtr; // TODO this can be improved somehow since eg for 1000 levels when going up this assignment will be done 1000 times with the same value.
#ifdef DEBUG_SUBDOMAIN
		printf("going up (containing domain PTR: %li\n\n\n",containingDomainPtr);		
#endif
		treeLevel--;
#ifdef DEBUG_SUBDOMAIN
		printf("lev: %i, CDL: %i\n",treeLevel,containingDomainLevel);
#endif

	return ptr;
}
/***************************************************************************************/
cpedsList<long> subDomain::getPointsIdx(const subDomain_region_t& r) {
	cpedsList<long> id;

//#ifdef DEBUG_SUBDOMAIN2
	/*
	 * Comment: this is a possible fix when searching tree that was decomposed down to a fixed depth and not
	 * down to a fixed maximal number of particles per cell.
	 * 
	 * author: blew
	 * date: Nov 27, 2017 4:00:19 PM
	 *
	 */
	
	subDomain* sD;
	vector<long>* rdomainPoints=getContainingDomainPointIdx(r,&sD,true);
//#endif	

	if (rdomainPoints!=NULL) {
		for (unsigned long i = 0; i < rdomainPoints->size(); i++) {
			if (inRegion(r,_points->at( rdomainPoints->at(i) ))) id.append(rdomainPoints->at(i));
		}
	}
	return id;
}
/***************************************************************************************/
mscsVector<long> subDomain::getPointsIdx(const cpedsPoint3D& p0, double r, subDomain* D, bool ignoreRegionSize, bool cubicalRegion) {
//	static mscsVector<long> id;
	mscsVector<long> id;
	id.clear();
	subDomain* tmpsD=D;
	if (tmpsD==NULL) tmpsD=this;
	// create parallel-piped region containing the sphere. We will look for the domain that completely contains this region
	// and then select points from that region that are in S(p0,r)
	subDomain_region_t reg;
	reg.xmin=p0.x()-r;	reg.xmax=p0.x()+r; 
	reg.ymin=p0.y()-r;	reg.ymax=p0.y()+r; 
	reg.zmin=p0.z()-r;	reg.zmax=p0.z()+r; 
	// chop off the parts that are outside of the domain
	if (reg.xmin<_def.rootDomain->getSubdomainInfo().region.xmin) reg.xmin=_def.rootDomain->getSubdomainInfo().region.xmin; if (reg.xmax>_def.rootDomain->getSubdomainInfo().region.xmax) reg.xmax=_def.rootDomain->getSubdomainInfo().region.xmax;
	if (reg.ymin<_def.rootDomain->getSubdomainInfo().region.ymin) reg.ymin=_def.rootDomain->getSubdomainInfo().region.ymin; if (reg.ymax>_def.rootDomain->getSubdomainInfo().region.ymax) reg.ymax=_def.rootDomain->getSubdomainInfo().region.ymax;
	if (reg.zmin<_def.rootDomain->getSubdomainInfo().region.zmin) reg.zmin=_def.rootDomain->getSubdomainInfo().region.zmin; if (reg.zmax>_def.rootDomain->getSubdomainInfo().region.zmax) reg.zmax=_def.rootDomain->getSubdomainInfo().region.zmax;
	vector< vector<long>* >* rdomains=NULL;
	bool deleteRdomains=false;
	vector<long>* rdomainPoints;
	subDomain* sD;
//	if (D!=NULL) sD=D; else sD=this;

#ifdef DEBUG_SUBDOMAIN2
	printf("search guess radius: %lE\n",r);
//	exit(0);
#endif

	if (true) { 
		rdomains=getOverlappingDomainPointIdx(reg,&sD,true);
	}
	else {
		rdomains = new vector< vector<long>* >;
		deleteRdomains=true;
		if (D==NULL) {
			rdomainPoints=getContainingDomainPointIdx(reg,&sD,true);
			if (ignoreRegionSize and sD==NULL) rdomainPoints=getSubdomainInfo().idx;
		}
		else {
			rdomainPoints=D->getContainingDomainPointIdx(reg,&sD,true);
			if (ignoreRegionSize and sD==NULL) rdomainPoints=D->getSubdomainInfo().idx;
		}
		rdomains->push_back(rdomainPoints);			
	}

	
//	rdomains = new vector< vector<long>* >;
//	deleteRdomains=true;
//	if (D==NULL) {
//		rdomainPoints=getContainingDomainPointIdx(reg,&sD,true);
//		if (ignoreRegionSize and sD==NULL) rdomainPoints=getSubdomainInfo().idx;
//	}
//	else {
//		rdomainPoints=D->getContainingDomainPointIdx(reg,&sD,true);
//		if (ignoreRegionSize and sD==NULL) rdomainPoints=D->getSubdomainInfo().idx;
//	}
//	rdomains->push_back(rdomainPoints);
//
//	if (rdomainPoints->size()>= _def.rootDomain->getSubdomainInfo().idx->size()/512) { 
//		/*
//		 * Comment: This causes that if halo is at the tree cells boundary above and inclusively third level (2^3^3)
//		 * then the search is done again to select smaller number of particles by searching not the region containing 
//		 * domain but all domains that partially contain the region.
//		 * 
//		 * If the halo resides on the tree subdomain walls below fourth level, then all particels from the level up
//		 * will be checked in the next loop.
//		 * 
//		 * author: blew
//		 * date: Jan 22, 2014 10:28:00 AM
//		 *
//		 */
//	
//		deleteRdomains=false;
//		delete rdomains;
//		rdomains=NULL;
//		rdomains=getOverlappingDomainPointIdx(reg,&sD,true);
//	}
	
	
//	printf("points in containing domain count: %li\n", rdomainPoints->size());
	long selectedPart=0;
	for (long i = 0; i < rdomains->size(); i++) { selectedPart+=rdomains->at(i)->size(); }
#ifdef DEBUG_SUBDOMAIN2
	printf("preselected particles: %li\n",selectedPart);
//	exit(0);
#endif
	
	// now search within subdomain for points matching criteria (must lie within spherical region or parallel-piped region)
	if (rdomains!=NULL) {
		
		for (unsigned long sDi = 0; sDi < rdomains->size(); sDi++) {
			rdomainPoints=rdomains->at(sDi);
			
			if (rdomainPoints!=NULL) {
				id.reserve(rdomainPoints->size());
				//#pragma omp critical
				if (cubicalRegion) {
					//			printf("%li\n",long(rdomainPoints->size()));
					//			if (cubicalRegion) exit(0);
					
					for (unsigned long i = 0; i < rdomainPoints->size(); i++) {
						if (inRegion(reg,_points->at( rdomainPoints->at(i) ))) { id.push_back(rdomainPoints->at(i)); }
					}
				}
				else {
					
					for (unsigned long i = 0; i < rdomainPoints->size(); i++) {
						if (inRange(r,p0,_points->at( rdomainPoints->at(i) ))) { id.push_back(rdomainPoints->at(i));
						//					printf("adding\n");
						}
					}			
				}
			}
		}
	}
	
	if (deleteRdomains) {
		delete rdomains;
		rdomains=NULL;
	}
	return id;	
}
/***************************************************************************************/
mscsVector<long> subDomain::getShellPointsIdx(const cpedsPoint3D& p0, double r, double dr,subDomain* D) {
		double R=r+dr;
		mscsVector<long> id;
		id.clear();
		subDomain* tmpsD=D;
		if (tmpsD==NULL) tmpsD=this;
		// create parallel-piped region containing the sphere. We will look for the domain that completely contains this region
		// and then select points from that region that are in S(p0,r)
		subDomain_region_t reg;
		reg.xmin=p0.x()-R;	reg.xmax=p0.x()+R; 
		reg.ymin=p0.y()-R;	reg.ymax=p0.y()+R; 
		reg.zmin=p0.z()-R;	reg.zmax=p0.z()+R; 
		// chop off the parts that are outside of the domain
		if (reg.xmin<_def.rootDomain->getSubdomainInfo().region.xmin) reg.xmin=_def.rootDomain->getSubdomainInfo().region.xmin; if (reg.xmax>_def.rootDomain->getSubdomainInfo().region.xmax) reg.xmax=_def.rootDomain->getSubdomainInfo().region.xmax;
		if (reg.ymin<_def.rootDomain->getSubdomainInfo().region.ymin) reg.ymin=_def.rootDomain->getSubdomainInfo().region.ymin; if (reg.ymax>_def.rootDomain->getSubdomainInfo().region.ymax) reg.ymax=_def.rootDomain->getSubdomainInfo().region.ymax;
		if (reg.zmin<_def.rootDomain->getSubdomainInfo().region.zmin) reg.zmin=_def.rootDomain->getSubdomainInfo().region.zmin; if (reg.zmax>_def.rootDomain->getSubdomainInfo().region.zmax) reg.zmax=_def.rootDomain->getSubdomainInfo().region.zmax;
		vector< vector<long>* >* rdomains=NULL;
		bool deleteRdomains=false;
		vector<long>* rdomainPoints;
		subDomain* sD;
	//	if (D!=NULL) sD=D; else sD=this;


		rdomains=getOverlappingDomainPointIdx(reg,&sD,true);


		
		long selectedPart=0;
		for (long i = 0; i < rdomains->size(); i++) { selectedPart+=rdomains->at(i)->size(); }
		
		// now search within subdomain for points matching criteria (must lie within spherical region or parallel-piped region)
		if (rdomains!=NULL) {
			
			for (unsigned long sDi = 0; sDi < rdomains->size(); sDi++) {
				rdomainPoints=rdomains->at(sDi);
				
				if (rdomainPoints!=NULL) {
					id.reserve(rdomainPoints->size());
						
					for (unsigned long i = 0; i < rdomainPoints->size(); i++) {
						if (inRange(R,p0,_points->at( rdomainPoints->at(i) )) and !inRange(r,p0,_points->at( rdomainPoints->at(i) ))) { id.push_back(rdomainPoints->at(i));
						}
					}			
				}
			}
		}
		
		if (deleteRdomains) {
			delete rdomains;
			rdomains=NULL;
		}
		return id;	
	
}

/***************************************************************************************/
void subDomain::print_domain_range(subDomain_region_t r) {
	printf("  --sub-domain region:\n");
	printf("    xmin: %lf, xmax: %lf\n",r.xmin,r.xmax);
	printf("    ymin: %lf, ymax: %lf\n",r.ymin,r.ymax);
	printf("    zmin: %lf, zmax: %lf\n",r.zmin,r.zmax);
	printf("    nx: %li, ny: %li, nz: %li\n",r.subx, r.suby, r.subz);
	printf("\n");
}
/***************************************************************************************/
void subDomain::saveDomains(string fileName) {
	FILE* f=fopen(fileName.c_str(), "w");
	if (f!=NULL) {
//		fprintf(f,"#new domain\n ");
		saveDomain(&f);
		fclose(f);
	}
	else { printf("cound not open the file\n"); exit(0); }
	
}
/***************************************************************************************/
void subDomain::saveDomain(FILE** f) {
//	for (long i = 0; i < _def.region.subx; i++) {
//		for (long j = 0; j < _def.region.suby; j++) {
//			for (long k = 0; k < _def.region.subz; k++) {
				if (_def.idx!=NULL) {
//					if (_def.idx->size() >= _def.minObjectsCount) {
#ifdef DEBUG
//						printf("%lf %lf %lf %lf %lf %lf %li %i\n",_def.region.xmin,_def.region.ymin,_def.region.zmin,_def.region.xmax,_def.region.ymax,_def.region.zmax,_def.idx->size(),int(_def.root));
#endif
						fprintf(*f,"%lf %lf %lf %lf %lf %lf %li %i %i\n",_def.region.xmin,_def.region.ymin,_def.region.zmin,_def.region.xmax,_def.region.ymax,_def.region.zmax,_def.idx->size(),int(_def.root),int(_def.leaf));

						
						if (_def.sub!=NULL){
							for (long i = 0; i < _def.subDomainsCount; i++) {
								if (_def.sub[i]!=NULL) _def.sub[i]->saveDomain(f);	
							}
						}
//					}
				}				
//			}
//		}
//	}	
}
/***************************************************************************************/
vector< vector<double> > subDomain::getBranchParticlesCount(bool leafsOnly) {
	vector< vector<double> > pcount;
	branchParticlesCount(pcount, leafsOnly);
	return pcount;
}
/***************************************************************************************/
void subDomain::branchParticlesCount(vector< vector<double> >& pcount, bool leafsOnly) {
	if (leafsOnly) {
		if (_def.leaf) {
			tmpvec.clear();
			tmpvec.push_back(_def.region.xmin);
			tmpvec.push_back(_def.region.xmax);
			tmpvec.push_back(_def.region.ymin);
			tmpvec.push_back(_def.region.ymax);
			tmpvec.push_back(_def.region.zmin);
			tmpvec.push_back(_def.region.zmax);
			tmpvec.push_back(_def.idx->size());
			tmpvec.push_back(1); // leaf flag
			
			// add a list of indexes of particles located in this subdomain
			for (unsigned long i = 0; i < _def.idx->size(); i++) {
				tmpvec.push_back(_def.idx->at(i)); 
			}
			
			pcount.push_back(tmpvec);
		}
	}
	else {
		//	if (_def.idx->size() >= _def.minObjectsCount) {
		tmpvec.clear();
		tmpvec.push_back(_def.region.xmin);
		tmpvec.push_back(_def.region.xmax);
		tmpvec.push_back(_def.region.ymin);
		tmpvec.push_back(_def.region.ymax);
		tmpvec.push_back(_def.region.zmin);
		tmpvec.push_back(_def.region.zmax);
		tmpvec.push_back(_def.idx->size());
		if (_def.leaf)	tmpvec.push_back(1); // leaf flag
		else tmpvec.push_back(0);
		
		// add a list of indexes of particles located in this subdomain
		for (unsigned long i = 0; i < _def.idx->size(); i++) {
			tmpvec.push_back(_def.idx->at(i)); 
		}
		
		
		pcount.push_back(tmpvec);
		//	}
	}
	if (_def.sub!=NULL)
		for (long i = 0; i < _def.subDomainsCount; i++) {
			if (_def.sub[i]!=NULL) _def.sub[i]->branchParticlesCount(pcount,leafsOnly);
		}
}
/***************************************************************************************/
mscsVector<long> subDomain::getNeighbors(cpedsPoint3D p0, long Nmin, long Nmax, double* hsml, long subDivisionsCount) {
	
//	if (Nmin>Nmax) { printf("mscsFunction3dregc::getNeighbors>> Nmin > Nmax. I better stop now."); exit(0); }
	mscsVector<long> l; 
//	static mscsVector<long> l; 
	l.clear();
	long n=0;
	double sidex,sidey,sidez;
	subDomain_region_t r=getSubdomainInfo().region;
	sidex=(r.xmax-r.xmin);
	sidey=(r.ymax-r.ymin);
	sidez=(r.zmax-r.zmin);
//	double rMax=cpeds_get_min(cpeds_get_min(sidex, sidey),sidez); // get the biggest possible radii in the domain
	double rMax=sqrt(sidex*sidex+sidey*sidey+sidez*sidez); // get the biggest possible radii in the domain
//	printf("rMax=%lf\n",rMax);
	double h, hhigh,hlow;
	hhigh=rMax;
	hlow=0;

//#pragma omp critical
	if (hsml==NULL) {		
//		h=rMax/2;
		h=pow(sidex*sidey*sidez/getSubdomainInfo().idx->size(),0.333)/2; // this is a guess value to begin with
//		h=pow(sidex*sidey/getSubdomainInfo().idx->size(),0.333); // this is a guess value to begin with
//		*hsml=(r.xmax-r.xmin + r.ymax-r.ymin + r.zmax-r.zmin)/12; // get the average side quarter length of the containing domain
//		*hsml=cpeds_get_min(cpeds_get_min(cpeds_get_min(*hsml,sidex), sidey),sidez); // make sure it's not too big
	}
	else h=*hsml;
//	printf("h=%lf\n",h);
//	exit(0);

//#pragma omp critical
	l=getPointsIdx(p0,h,NULL,true); // get indexes of points in S(p0,h)
	n=l.size();
	int iterations=0;
	long iterationsMax=subDivisionsCount;
	bool cont=true;
	double step=h/2;
	bool dirChanged=false;
	bool increasing=true;
//	if (n>Nmax) increasing=false; else increasing=true;
//#pragma omp critical
	while ( (n<Nmin or n>Nmax) and cont==true) {
#ifdef DEBUG_SUBDOMAIN_NEIGHBORS
		printf("-- new loop iteration: %li, h=%lf,  step: %lf, iterationsMax: %li, hhigh: %lf\n",iterations, h, step, iterationsMax, hhigh);
#endif
		if (iterations<iterationsMax) {
			if (n>Nmax) {  
				h-=step; if (increasing) dirChanged=true; 
#ifdef DEBUG_SUBDOMAIN_NEIGHBORS
				printf("decreasing radii\n");
#endif
			}
			if (n<Nmin) { 
				h+=step; if (!increasing) dirChanged=true; 
#ifdef DEBUG_SUBDOMAIN_NEIGHBORS
				printf("increasing radii\n");
#endif
			}
#ifdef DEBUG_SUBDOMAIN_NEIGHBORS
			printf("-- now iteration: %li, h=%lf, step: %lf, iterationsMax: %li, hhigh: %lf\n",iterations, h, step, iterationsMax, hhigh);
#endif
			if (dirChanged) { 
				step/=2; increasing=!increasing; dirChanged=false; iterations++; 
#ifdef DEBUG_SUBDOMAIN_NEIGHBORS
				printf("dirchanged, dividing step\n"); 
#endif
			}
		}
//		if (n>=Nmin) {
//			hhigh=h; 
//			hlow=h;
//		}
//					*hsml=cpeds_get_min(cpeds_get_min(cpeds_get_min(*hsml,sidex), sidey),sidez); // make sure it's not too big
//		if (iterations==iterationsMax or h>rMax) {
		if (iterations==iterationsMax) {
			cont=false;
//#pragma omp critical
			l=getPointsIdx(p0,h,this,true);
			n=l.size();			
		}
		else {
//#pragma omp critical
			l=getPointsIdx(p0,h,this);
			n=l.size();
		}
		
		// in case when it is possible to find a radii enclosing the requested Nmin particles use this radii as the highest radii in the worst case scenario i.e. when bail-out condition will be met (cont=false)
		if (n>=Nmin and h<hhigh) { hhigh=h; }
//		if (n<=Nmax and h>hlow) { hhigh=h; }
#ifdef DEBUG_SUBDOMAIN_NEIGHBORS
		printf("   -- iteration: %li, h=%lf, n: %li, step: %lf, iterationsMax: %li, hhigh: %lf\n",iterations, h,n, step, iterationsMax,hhigh);
#endif
	}
//#pragma omp critical
	if (iterations==iterationsMax) { h=hhigh;	l=getPointsIdx(p0,h,this); }
//	printf("-- iteration end: %li, h=%lf, n: %li, iterationsMax: %li, rMax: %lf, acc: %lf, hhigh: %lf\n",iterations, h,l.size(), iterationsMax,rMax,acc, hhigh);
	
	// if the hsml was given set it to the largest distance between p0 and the set of found Nmin points
	if (hsml!=NULL) { 
		h=0;
		double dist;
//#pragma omp critical
		for (unsigned long i = 0; i < l.size(); i++) {
			dist=p0.dist(_points->at(l[i]));
			if (dist>h) h=dist;
		}
//#pragma omp critical
		*hsml=h; 
//		printf("h: %lf\n",h);
	}
	return l;
}
