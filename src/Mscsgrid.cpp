/*!
 * \file Mscsgrid.cpp
 *
 *  Created on: May 18, 2010
 *      Author: blew
 */

#include "Mscsgrid.h"
#include "Mscs-function.h"

/***************************************************************************************/
mscsGrid::mscsGrid(const string name, cpeds_VerbosityLevel verbosity) : mscsObject(name, verbosity) {
	srcMap=NULL;
}

/***************************************************************************************/
mscsGrid::~mscsGrid() {

}


/***************************************************************************************/
mscsGrid::mscsGrid(mscsMap* map, const string name, cpeds_VerbosityLevel verbosity) : mscsObject(name, verbosity) {
	srcMap=map;
}
/***************************************************************************************/
const matrix<double>& mscsGrid::makeGrid(gridType_t g, long Nphi, long Nth) {

	dth_=PI/(Nth+1);
	dphi_=twoPI/Nphi;
	Nphi_=Nphi;
	Nth_=Nth;
	grid_.SetSize(Nth_,Nphi_);
	gridType=g;

	//
	// define grid
	//
	double phi=0;
	for (int i=0;i<Nphi_;i++) {
		phiVals_.append(phi);
		phi+=dphi_;
	}
	double th=0;
	for (int i=0;i<Nth_;i++) {
		th+=dth_; // do not include the poles
		thVals_.append(th);
	}

	return grid();
}

/***************************************************************************************/
const matrix<double>& mscsGrid::makeGrid(gridType_t g, double dphi, double dth) {
	dth_=dth;
	dphi_=dphi;
	Nth_=long(trunc(PI/dth_));
	dth_=twoPI/Nth_;
	Nphi_=long(trunc(twoPI/dphi_));
	dphi_=twoPI/Nphi_;
	grid_.SetSize(Nth_,Nphi_);
	gridType=g;

	//
	// define grid
	//
	double phi=0;
	for (int i=0;i<Nphi_;i++) {
		phiVals_.append(phi);
		phi+=dphi_;
	}
	double th=0;
	for (int i=0;i<Nth_;i++) {
		thVals_.append(th);
		th+=dth_;
	}

	return grid();
}

/***************************************************************************************/
const matrix<double>& mscsGrid::grid() {
	double acc=1e-10;
	//
	// prepare the source map ordering
	//

	if (!srcMap->coordLoaded()) {
		if ( srcMap->ordering()==srcMap->ring) { srcMap->conv_ring2nest(); }
		srcMap->set_map_coord();
	}

	// convert to ring
	if ( srcMap->ordering()==srcMap->nested) { srcMap->conv_nest2ring(); }

	long ring_num = cpeds_get_ring_num_healpix(srcMap->nside());

	srcMap->mask_map_merge();

	//
	// perform ring by ring interpolation at given longitudes
	//
	matrix<double> m(ring_num,Nphi_);

	mscsFunction ring("ring");

	for (int i=0;i<ring_num;i++) {
		//		printf("RING NUM:%i\n", i);
		ring=srcMap->getRing(i);
		//		ring.scaleX(PI180inv).print();
		//		ring.scaleX(PI180);

		// sort in ring
		ring.sortFunctionArgAscending();
		// make periodic
		if (ring.getX(0)>acc) ring.newPoint(ring.getX(0)+twoPI,ring.getY(0));
		//			ring.save("test-map.txt");
		// interpolate on the whole range and the range that will be folded to complete the period
		ring=ring.interpolate(0,ring.getX(0)+twoPI,dphi_,false,"linear",true);

		// cast function arguments into [0..twoPI) range. Upper limit within accuracy 1e-10
		ring.foldXinto(0,twoPI,1e-10);
//			ring.save("test-inter.txt");
		//		if (i==3) exit(0);
		ring.uniqueXFast(1e-10);
		// check interpolation
		if (ring.pointsCount()!=Nphi_) { msgs->criticalError("wrong number of interpolated points in phi direction: "+msgs->toStr(ring.pointsCount())+" should be "+msgs->toStr(Nphi_)+ " - CHECK THE CODE",High); }

		// store in matrix
		for (int j=0;j<Nphi_;j++) {			m(i,j)=ring.Y(j);		}

	}


	//
	// perform column by column interpolation at given latitudes
	//
	for (int i=0;i<grid_.ColNo();i++) {
		// copy column onto a function for interpolation
		ring.clearFunction();
		for (int j=0;j<ring_num;j++) {
			ring.newPoint(PIsnd-srcMap->get_C(  cpeds_get_pix_num_above_ring_healpix(srcMap->nside(),j)  ).b(),m(j,i));
		}
		// sort in ring
		ring.sortFunctionArgAscending();
		//		ring.print();
		// interpolate
		ring=ring.interpolate(thVals_.first(),thVals_.last(),dth_,false,"linear",true,1e-10);
		// check if function is defined from thVals_.first() to thVals_.last()
//		ring.setVerbosityLevel(High);
//		thVals_.print();
//		ring.print();
//		printf("minArg: %lE maxArg %lE, %lE %lE\n",ring.getMinArg(), ring.getMaxArg(), PIsnd-srcMap->get_C(0).b(),PIsnd-srcMap->get_C(srcMap->pixNum()-1).b());
//		printf("iminArg: %li \n",ring.getMinArgIdx());

		//
		// extrapolate if needed outside of the range where map data is defined (will only copy the extremal values around the poles)
		//
		long N1=(Nth_-ring.pointsCount());
		if (thVals_.first()<PIsnd-srcMap->get_C(0).b()) {
			if (N1 % 2 != 0) { msgs->criticalError("upper co-latitude extrapolation: number of missing points is odd. This shouldn't happen. CHECK THE CODE",High); }
			N1/=2;
//			printf("number of extra points from the top: %li\n", N1);
			for (int j=0;j<N1;j++) { ring.insertPoint(thVals_[i],ring.f(0),0);  }
		}
		long N2=(Nth_-ring.pointsCount());
		if (thVals_.last()>PIsnd-srcMap->get_C( srcMap->pixNum()-1 ).b()) {
			if (N1!=N2) msgs->criticalError("lower co-latitude extrapolation: number of missing points is different than in upper co-latitude extrapolation. This shouldn't happen with healpix map. CHECK THE CODE",High);
//			printf("number of extra points from the bottom: %li\n", N2);
			for (int j=0;j<N2;j++) { ring.newPoint(thVals_.last(),ring.f(ring.pointsCount()-1));  }
		}
		if (ring.pointsCount()!=Nth_) { msgs->criticalError("wrong number of interpolated points in th direction: "+msgs->toStr(ring.pointsCount())+" should be "+msgs->toStr(Nth_)+ " - CHECK THE CODE",High); }
		// store back onto the matrix
		for (int j=0;j<Nth_;j++) {			grid_(j,i)=ring.Y(j);		}
	}

	m.SetSize(0,0);

	return grid_;

}
/***************************************************************************************/
double mscsGrid::getVal(double th, double phi, bool check) {
	long i,j;
	if (check) cpeds_check_thphi(&th,&phi);

	i=long(trunc(th/dth_));
	j=long(trunc(phi/dphi_));

	return grid_(i,j);
}
/***************************************************************************************/
mscsMap mscsGrid::toMap(long nside) {
	mscsMap map("map");

	map.set_nside(nside);
	long N=map.pixNum();
	map.set_map_coord();
	map.makekill_space_manager("make","T");
	cpedsDirection n;
	for (long i=0;i<N;i++) { n=map.get_C(i); map.set_T(i, getVal( PIsnd-n.b(),n.l() )); }
	map.makekill_space_manager("kill","C");
	return map;
}
/***************************************************************************************/
// this function is based on extrapolation in order to restore the power which would
// otherwise be lost via averaging

// THIS IS NOT FINISHED
mscsMap mscsGrid::makeMap() {
	mscsMap map("map");
	long ring_num = cpeds_get_ring_num_healpix(srcMap->nside());
	matrix<double> m(ring_num,Nphi_);

	//
	// extrapolate values on the original rings
	//
	mscsFunction merid("merid");
	merid.setX(thVals_);
	long Nth;
	double * th=cpeds_get_ring_colatitudes_healpix(srcMap->nside(),&Nth);
	mscsFunction extra("extra");
	extra.setX(cpedsList<double>(th,Nth));
	long N,k1,k2;
	double x1,x2,y1,y2;
	for (int i=0;i<grid_.ColNo();i++) {
		// copy column onto a function for extrapolation
		for (int j=0;j<Nth_;j++) { merid.setf(j,grid_(j,i)); }
		//
		// extrapolate
		//
		N=extra.pointsCount();
		for (int j = 0; j < N; j++) {
			merid.f(extra.X(j),&k1);
			if (k1==0) k2=1; else k2=k1-1;
			x1=merid.X(k1);
			x2=merid.X(k2);
			y1=merid.Y(k1);
			y2=merid.Y(k2);

		}
		// find suitable points
//		x1=merid()
	}

//	map.set_nside(nside);
//	long N=map.pixNum();
//	map.set_map_coord();
//	map.makekill_space_manager("make","T");
//	cpedsDirection n;
//	for (long i=0;i<N;i++) { n=map.get_C(i); map.set_T(i, getVal( PIsnd-n.b(),n.l() )); }
//	map.makekill_space_manager("kill","C");
	return map;

}
