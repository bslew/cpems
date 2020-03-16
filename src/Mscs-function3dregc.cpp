/*
 * mscsFunction3d.cpp
 *
 *  Created on: Dec 22, 2010, 10:43:15 PM
 *      Author: blew
 */
//#define _REENTRANT
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <sys/dir.h>
//#include <pthread.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <limits>
#include <assert.h>
//#include <H5Cpp.h>
//#include <fftw3.h>
#include "Mscs-function3dregc.h"
#include "Mscs-map-window_function.h"
#include "pointsDensity.h"
#include "cpeds-math.h"
#include "cpeds-point3d.h"

// dependences for 2d sibson triangulation based interpolation
#include <CGAL/basic.h>
#include <utility>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/Interpolation_gradient_fitting_traits_2.h>
//#include <CGAL/sibson_gradient_fitting.h>
#include <CGAL/interpolation_functions.h>

// dependences for 2d linear triangulation based interpolation
#include <CGAL/Interpolation_traits_2.h>

struct K : CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef std::vector< std::pair< K::Point_2, K::FT  > > Point_coordinate_vector;
/***************************************************************************************/
mscsFunction3dregc::mscsFunction3dregc() : mscsObject("function3Dregc",Zero) {
	//	plan_r2c=new fftw_plan;
	initiateParams();
	initiate(false);
}
/***************************************************************************************/
mscsFunction3dregc::mscsFunction3dregc(string name, cpeds_VerbosityLevel verbosityLoc) : mscsObject(name,verbosityLoc) { 
	//	plan_r2c=new fftw_plan;
	initiateParams();
	initiate(false);
}
/***************************************************************************************/

mscsFunction3dregc::mscsFunction3dregc(long Nx, long Ny, long Nz, double dx,double dy, double dz, double x0, double y0, double z0, cpeds_VerbosityLevel verbosityLoc) : mscsObject("function3Dregc",verbosityLoc) {
	initiateParams();
	
	setSize(Nx,Ny,Nz,dx,dy,dz,x0,y0,z0);
	// initiate the function domain
	//	double tmp=_param.dxo2;
	//	for (int i = 0; i < Nx; i++) {
	//		_X.append(tmp);
	//		tmp+=dx;
	//	}
	//	
	//	tmp=_param.dyo2;
	//	for (int i = 0; i < Ny; i++) {
	//		_Y.append(tmp);
	//		tmp+=dy;
	//	}
	//	
	//	tmp=_param.dzo2;
	//	for (int i = 0; i < Nz; i++) {
	//		_Z.append(tmp);
	//		tmp+=dz;
	//	}
	initiate(true);
	//	allocFunctionSpace();
}
/***************************************************************************************/
mscsFunction3dregc::mscsFunction3dregc(const mscsFunction3dregc& parent) : mscsObject(parent) {
	initiateParams();
	initiate(false);
	
	long n=parent.size();
	_param=parent.getFunctionParameters();
	if (parent.data()!=NULL) {
		allocFunctionSpace();
		//		for (long i = 0; i < n; i++) { _data[i][0]=parent[i][0]; _data[i][1]=parent[i][1];	}
		dataHardCopy(parent.data());
	}
	else {
		freeFunctionSpace();
	}
	
}

/***************************************************************************************/
mscsFunction3dregc::~mscsFunction3dregc() {
	freeFunctionSpace();
	delete _rns;
	destroyFFTPlans();
	//#ifdef HAVE_OPENMP
	//	if (omp_hdf_lock!=NULL)
	//		omp_destroy_lock(omp_hdf_lock);
	//
	//#endif
}
/***************************************************************************************/
void mscsFunction3dregc::initiateParams() {
#ifndef NO_HDF5
	initiatie_hdf5_params();
#endif
	
	_param.NXY=0;
	_param.NXYZ=0;
	_param.NYZ=0;
	_param.Nx=0;
	_param.Ny=0;
	_param.Nz=0;
	_param.dx=0;
	_param.dy=0;
	_param.dz=0;
	_param.dxo2=0;
	_param.dyo2=0;
	_param.dzo2=0;
	_param.sizeX=0;
	_param.sizeY=0;
	_param.sizeZ=0;
	_param.x0=0;
	_param.y0=0;
	_param.z0=0;
	_param.xMax=0;
	_param.yMax=0;
	_param.zMax=0;	
	
	_param.derivativeXperiodic=true;
	_param.derivativeYperiodic=true;
}
/***************************************************************************************/
void mscsFunction3dregc::initiate(bool alloc) {
	_rns = new cpedsRNG("gaussian_circle");
	_data=NULL;
	plan_r2c_forward=NULL;
	plan_r2c_backward=NULL;
	plan_c2c_forward=NULL;
	plan_c2c_backward=NULL;
	
	setSavePrecision(5);
	
	if (alloc) allocFunctionSpace();
}
/***************************************************************************************/
int mscsFunction3dregc::destroyFFTPlans() {
	int ret=0;
	if (plan_c2c_forward!=NULL) {  ret=1; 	}
	if (plan_c2c_backward!=NULL) {  ret=1; }
	if (plan_r2c_forward!=NULL) {  ret=1; }
	if (plan_r2c_backward!=NULL) {  ret=1; }
#ifdef DEBUG_FUNC3DREGC_FFTWPLANS
	if (ret) msgs->say("FFTW PLANS WILL BE DELETED",Top);
#endif
	
	if (plan_c2c_forward!=NULL) { 
		fftw_destroy_plan(*plan_c2c_forward); 
		delete plan_c2c_forward; plan_c2c_forward=NULL; ret=1; 
	}
	if (plan_c2c_backward!=NULL) { fftw_destroy_plan(*plan_c2c_backward); 
	delete plan_c2c_backward; plan_c2c_backward=NULL; ret=1; }
	if (plan_r2c_forward!=NULL) { fftw_destroy_plan(*plan_r2c_forward); 
	delete plan_r2c_forward; plan_r2c_forward=NULL; ret=1; }
	if (plan_r2c_backward!=NULL) { fftw_destroy_plan(*plan_r2c_backward); 
	delete plan_r2c_backward; plan_r2c_backward=NULL; ret=1; }

//	fftw_cleanup();
	return ret;
}
/***************************************************************************************/
void mscsFunction3dregc::freeFunctionSpace() {
#ifdef DEBUG_FUNC3DREGC_FFTWPLANS
	msgs->say("FREEING FUNCTION SPACE",Top);
#endif
	int r=0;
	if (_data!=NULL) { 
		r=destroyFFTPlans(); // we need to destroy plans too in case they were already used, because the plan stores the location of the _data
		fftw_free(_data); 
		_data=NULL; 
	}
#ifdef DEBUG_FUNC3DREGC_FFTWPLANS
	if (r ) { printf("fftw plans were destroyed\n");  }
#endif

}
/***************************************************************************************/
void mscsFunction3dregc::allocFunctionSpace(long n) {
	freeFunctionSpace();
	
	if (n==0) {
		//		printf("getN: %li\n",getN());
		if (getN()>0) {
			msgs->say("allocating "+msgs->toStr(getN())+" cells to store the function.",Low);
			_data=(fftw_complex*)fftw_malloc((getN())*sizeof(fftw_complex)); 
		}
		else {
			msgs->criticalError("cannot allocate the function space. No or wrong size specification was given (n="+msgs->toStr(n)+").",High);
		}
	}
	else {
		_data=(fftw_complex*)fftw_malloc(n*sizeof(fftw_complex));
	}
}

/***************************************************************************************/
void mscsFunction3dregc::setf(long i, long j, long k, double x, double y) {
	long idx=ijk2idx(i,j,k);
	setf(idx,x,y);
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::setf(const double re, const double im) {
	long n=getN();
	
	for (long i = 0; i < n; i++) {
		fRe(i)=re;
		fIm(i)=im;
	}
	
	return *this;
	
}
/***************************************************************************************/
void mscsFunction3dregc::getYrange(double *ymin, double *ymax) const {
	*ymin=_param.y0+_param.dyo2;
	*ymax=*ymin+Ny()*_param.dy;
}

/***************************************************************************************/
void mscsFunction3dregc::getZrange(double *zmin, double *zmax) const {
	*zmin=_param.z0+_param.dzo2;
	*zmax=*zmin+Nz()*_param.dz;
}
/***************************************************************************************/
void mscsFunction3dregc::getXrange(double *xmin, double *xmax) const {
	*xmin=_param.x0+_param.dxo2;
	*xmax=*xmin+Nx()*_param.dx;
}
/***************************************************************************************/
mscsFunction3dregc & mscsFunction3dregc::mk3DRegularGrid(double fromX, double toX, double fromY, double toY, double fromZ, double toZ, double dx, double dy, double dz, double vre, double vim) {
	setSize((toX-fromX)/dx,(toY-fromY)/dy,(toZ-fromZ)/dz,dx,dy,dz,fromX,fromY,fromZ);
	allocFunctionSpace();
	long nx=Nx();
	long ny=Ny();
	long nz=Nz();
	
	for (long i = 0; i < nx; i++) {
		for (long j = 0; j < ny; j++) {
			for (long k = 0; k < nz; k++) {
				fRe(i,j,k)=vre;
				fIm(i,j,k)=vim;
			}
		}
	}
	
	return *this;
}

/***************************************************************************************/
long mscsFunction3dregc::ijk2idx(long i, long j, long k) const {
	//	return k*_param.sizeXY+j*_param.sizeY+i;
	return k+j*_param.Nz+i*_param.NYZ;
}
/***************************************************************************************/
long mscsFunction3dregc::kji2idx(long i, long j, long k) const {
	return i+j*_param.Nx+k*_param.NXY;	
}
/***************************************************************************************/
void mscsFunction3dregc::idx2ijk(long idx, long& i, long& j, long& k) const {	
	i=idx / _param.NYZ;
	j=(idx-i*_param.NYZ) / _param.Nz;
	k=idx-i*_param.NYZ-j*_param.Nz;
}

/***************************************************************************************/
cpedsStatusCodes mscsFunction3dregc::load(string filename, bool commentedFile) {
	mscsFunction tmp;
	tmp.load(filename,commentedFile,0,1);
	if (tmp.pointsCount()!=size()) { msgs->error("mscsFunction3dregc::load> Something is WRONG ! The file size and the allocated function size do NOT match.",Top); }
	if (tmp.pointsCount()>size()) { 
		msgs->say("found %li values in file",tmp.pointsCount(),High);
		msgs->criticalError("mscsFunction3dregc::load> The file size is too big. Cannot load it onto the function.",Top); }
	importFunction(tmp);
	//	cpedsList<
	
	//	FILE* f;
	//	string format;
	//	double x,y;
	//	int res;
	//	msgs->say("Loading function from file "+filename,High);
	//	long cols=cpeds_get_cols_num_first_ln(filename.c_str(),commentedFile);
	//	msgs->say("Number of columns in the input file: "+msgs->toStr(cols),Medium);
	//	if (cols==1) { msgs->say("Assuming real-valued function.", Medium);	}
	//	if (cols==2) { msgs->say("Assuming complex-valued function.", Medium);	}
	//
	//	switch (cols) {
	//		case 1:
	//			format="%lE";
	//			break;
	//		case 2:
	//			format="%lE %lE";
	//
	//			break;
	//		default:
	//			break;
	//	}
	////	for (int i = 2; i < cols; ++i) {	format+=" %*E";  }
	//	if (commentedFile) {
	//		ifstream ifs(filename.c_str());
	//		if (ifs.is_open()) {
	//			string s;
	//			switch (cols) {
	//				case 1:
	//					while (getline(ifs,s)) {
	//						if (s[0]!= '#') { sscanf(s.c_str(),format.c_str(),&x); newPoint(x,0); }
	//					}
	//					break;
	//				case 2:
	//					while (getline(ifs,s)) {
	//						if (s[0]!= '#') { sscanf(s.c_str(),format.c_str(),&x,&y); newPoint(x,y); }
	//					}
	//
	//					break;
	//			}
	//		}
	//		else { msgs->error("ERROR: cannot load from the file: "+filename+". No such file", High);	  return cpedsNoSuchFile;	  }
	//	}
	//	else {
	//		f= fopen(filename.c_str(),"r");
	//		if (f!=NULL) {
	//			switch (cols) {
	//				case 1:
	//					while (!feof(f)) { res=fscanf(f,format.c_str(),&x); if (res!=EOF) newPoint(x,0); }
	//					break;
	//				case 2:
	//					while (!feof(f)) { res=fscanf(f,format.c_str(),&x,&y); if (res!=EOF) newPoint(x,y); }
	//
	//					break;
	//			}
	//			fclose(f);
	//			msgs->say("Read "+msgs->toStr((long)_f.count())+" lines",Medium);
	//		}
	//		else { msgs->error("ERROR: cannot load from the file: "+filename+". No such file", High);	  return cpedsNoSuchFile;	  }
	//	}
	//	
	//		
	//	long n=pointsCount();
	//	msgs->say("Loaded "+msgs->toStr(pointsCount())+" points.",Medium);
	//	if (n!=getN()) {
	//		msgs->warning("The expected number of loaded points is different from the declared grid size: "+msgs->toStr(getN()),High);
	//	}
	return cpedsSuccess;
	
}
/***************************************************************************************/
cpedsStatusCodes mscsFunction3dregc::loadMatrix(string filename) {
	matrix<double> m;
	m=cpeds_matrix_load(filename);
	//	exit(0);
	if (m.ColNo()==0 or m.RowNo()==0) {
		msgs->say("Incorrect matrix size detected, will not load. The size of the function remains unchanged.",High);
		return cpedsSuccess;
	}
	
	setSize(m.ColNo(),m.RowNo(),1,1,1,0,0,0,0);
	allocFunctionSpace();
	
	for (long i = 0; i < Nx(); i++) {
		for (long j = 0; j < Ny(); j++) {
			fRe(i,j,0)=m(j,i);
		}
	}
	
	return cpedsSuccess;
}
/***************************************************************************************/
cpedsStatusCodes mscsFunction3dregc::loadtxtLin(string filename,int how) {
	FILE* f=fopen(filename.c_str(),"r");
	if (f==NULL) {msgs->error("mscsFunction3dregc::loadtxtlin> cannot open file for reading",High); return cpedsError; }
	long nx,ny,nz;
	double re,im=0;
	nx=Nx();
	ny=Ny();
	nz=Nz();
	for (long i = 0; i < nx; i++) {
		for (long j = 0; j < ny; j++) {
			for (long k = 0; k < nz; k++) {
				fscanf(f,"%lE",&re);
				if (how==0) im=0;
				if (how==1) fscanf(f,"%lE",&im);
				setf(i,j,k,re,im);
			}
		}
	}
	fclose(f);
	//}
	//else {
	//	cpeds_save_matrix(exportRe(),getN(),1,filename+".re",true);
	//	cpeds_save_matrix(exportIm(),getN(),1,filename+".im",true);		
	//}
	//	
	return cpedsSuccess;
}
/***************************************************************************************/
cpedsStatusCodes mscsFunction3dregc::savetxtlin(string filename, bool all) {
	//	msgs->say("savetxtlin is not implemented yet");
	if (all) {
		FILE* f=fopen(filename.c_str(),"w");
		if (f==NULL) {msgs->error("mscsFunction3dregc::savetxtlin> cannot open file for saving",High); return cpedsError; }
		long nx,ny,nz;
		nx=Nx();
		ny=Ny();
		nz=Nz();
		for (long i = 0; i < nx; i++) {
			for (long j = 0; j < ny; j++) {
				for (long k = 0; k < nz; k++) {
					fprintf(f,"%lE %lE %lE %lE %lE\n",getX(i),getY(j),getZ(k),fRe(i,j,k),fIm(i,j,k));
				}
			}
		}
		fclose(f);
	}
	else {
		cpeds_save_matrix(exportRe(),getN(),1,filename+".re",true);
		cpeds_save_matrix(exportIm(),getN(),1,filename+".im",true);		
	}
	
	return cpedsSuccess;
	
}
/***************************************************************************************/
cpedsStatusCodes mscsFunction3dregc::save(string filename) {
	//	msgs->say("savetxtlin is not implemented yet");
	//	cpeds_save(reinterpret_cast<double**>(_data),size(),filename,false);
	cpeds_save((fftwComplex*)_data,size(),filename,false);
	filename+=".fparams";
	FILE* f=fopen(filename.c_str(),"w");
	fprintf(f,"%li %li %li\n",_param.Nx,_param.Ny,_param.Nz);
	fprintf(f,"%lf %lf %lf\n",_param.sizeX,_param.sizeY,_param.sizeZ);
	fprintf(f,"%lf %lf %lf\n",_param.x0,_param.y0,_param.z0);
	fclose(f);
	return cpedsSuccess;
}
/***************************************************************************************/
cpedsStatusCodes mscsFunction3dregc::saveDF3(string filename, int part) {
	//	filename+=".df3";
	FILE* f=fopen(filename.c_str(),"wb");
	unsigned short int n;
	n=(unsigned short int)(Nx()); 
	printf("field size X: %i\n",n);
	n=(n << 8) | (n >> 8);
	fwrite(&n,sizeof(n),1,f);
	n=(unsigned short int)(Ny()); 
	printf("field size Y: %i\n",n);
	n=(n << 8) | (n >> 8);
	fwrite(&n,sizeof(n),1,f);
	n=(unsigned short int)(Nz()); 
	printf("field size Z: %i\n",n);
	n=(n << 8) | (n >> 8);
	fwrite(&n,sizeof(n),1,f);
	int v;
	long l=0;
	double* d;
	double minv,maxv;
	double minv2,maxv2;
	
	
	if (part==0) {
		d=exportRe();
	}
	if (part==1) {
		d=exportIm();
	}
	cpeds_find_minmax_value(d,getN(),&minv,&maxv,NULL,NULL);
#ifdef DEBUG_DF3
	printf("field minv: %lE maxV: %lE\n",minv,maxv);
#endif
	cpeds_sub_value(minv,d,getN());
#ifdef DEBUG_DF3
	cpeds_find_minmax_value(d,getN(),&minv2,&maxv2,NULL,NULL);
	printf("minv: %lE maxV: %lE\n",minv2,maxv2);
#endif
	cpeds_mul_value(double(2.0*numeric_limits<int>::max())/(maxv-minv),d,getN());
	//	cpeds_mul_value(double(1.0)/(maxv-minv),d,getN());
#ifdef DEBUG_DF3
	cpeds_find_minmax_value(d,getN(),&minv2,&maxv2,NULL,NULL);
	printf("minv: %lE maxV: %lE\n",minv2,maxv2);
#endif
	cpeds_sub_value(numeric_limits<int>::max(),d,getN());
#ifdef DEBUG_DF3
	cpeds_find_minmax_value(d,getN(),&minv2,&maxv2,NULL,NULL);
	printf("corrected minv: %lE maxV: %lE\n",minv2,maxv2);
#endif
	
	for (long k = 0; k < Nz(); k++) {
		for (long j = 0; j < Ny(); j++) {
			for (long i = 0; i < Nx(); i++) {
				v=int(d[l]); // convert to int
#ifdef DEBUG_DF3
				printf("f: %lE i: %i, factor:%li\n",d[l],v, numeric_limits<int>::max());
#endif
				v= (v << 24) | ((v << 8) & 0x00ff0000) | ((v >> 8) & 0x0000ff00) | (v >> 24); // change to big-endian format
				fwrite(&v,sizeof(v),1,f);
				l++;
#ifdef DEBUG_DF3
				//				setf(i,j,k,d[l],0);				
#endif
			}
		}
	}
	delete [] d;
	fclose(f);
	return cpedsSuccess;
}

/***************************************************************************************/
cpedsStatusCodes mscsFunction3dregc::saveAllSlices(string dirPref, string filenamePref, int plane, int part) {
	if (getN()==0) return cpedsError;
	long N;
	if (plane==0) N=_param.Nx;
	if (plane==1) N=_param.Ny;
	if (plane==2) N=_param.Nz;
	
	int ret=mkdir(dirPref.c_str(),0755);
	if (ret==0) {
		for (long i = 0; i < N; i++) {
			cpeds_matrix_save(getSlice(plane,i,part),dirPref+"/"+filenamePref+"-plane_"+msgs->toStr(plane)+"."+msgs->toStr(i,"%04li"));
		}
	}
	else {
		msgs->warning("could not create directory "+dirPref+". It may already exist. Will not save.", Top);
	}
	return cpedsSuccess;
}
/***************************************************************************************/
cpedsStatusCodes mscsFunction3dregc::saveSlice(int plane, int coord, string filename, int part) {
	if (getN()==0) { return cpedsError; }
	cpeds_matrix_save(getSlice(plane,coord,part),filename,"",_param.savePrecision);
	return cpedsSuccess;
}
/* ******************************************************************************************** */
void mscsFunction3dregc::setSavePrecision(int precision) {
	_param.savePrecision=precision;
}
/* ******************************************************************************************** */
int mscsFunction3dregc::getSavePrecision() const {
	return _param.savePrecision;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::generateRandomGaussianField(double m, double s, bool re, bool im, long seed, long seedOffset) {
	//	cpedsRNG *rns = new cpedsRNG("gaussian_circle");
	if (_rns->getRNsType()!=_rns->gaussian_circle) {
		_rns->setRNsType(_rns->gaussian_circle);
	}
	_rns->setMeanStd(m,s);
	
	if (seed!=0) {
		_rns->seed(seed);
	}
	if (seedOffset!=0) {
		_rns->seedOffset(seedOffset);
	}
	
	double *a;
	if (re) {
		a = _rns->getRN(getN());
		importFunction(a,getN(),true);
		delete [] a;
	}
	if (im) {
		a = _rns->getRN(getN());
		importFunction(a,getN(),false);
		delete [] a;
	}
	
	return *this;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::generateRandomGaussianField(mscsFunction& Pk, long seed) {
	msgs->warning("mscsFunction3dregc::generateRandomGaussianField: THIS IMPLEMENTATION HAS NOT BEEN TESTED. PLEASE PERFORM TESTS. In particular re2im is suspicious (see function_fft on how to do this correctly)",Top);
	allocFunctionSpace();
	
	mscsFunction3dregc ker(*this);
	ker.mkKernel(Pk,"linear");
	
	msgs->say("generating gaussian field",High);
	generateRandomGaussianField(0,1,true,true,seed,0);
	//	long N=size();
	//	exit(0);
	
	msgs->say("Forming the power spectrum",High);
	//	l=0;
	
	ker/=double(sqrt(2.0)); 
	(*this)*=ker;
	
	//	for (long i = 0; i < nx; i++) {
	//		for (long j = 0; j < ny; j++) {
	//			for (long k = 0; k <= nz; k++) {
	//				
	//				modk=kk[l];
	//				if (Pkint[l]<0) {  // this is to get rid of possible residual oscillations from non-linear interpolation
	//					msgs->warning("Pkint is "+msgs->toStr(Pkint[l])+" < 0. Assuming power 0 for k-mode: (i,j,k)= ("+msgs->toStr(i)+", "+msgs->toStr(j)+", "+msgs->toStr(k)+")",High);
	//					Pkint[l]=0;
	//				}
	//
	//				//				Kmax=sqrt(pow(nx-1,2)+pow(ny-1,2)+pow(nz,2));
	////				delta=sqrt(Pkint[l])/sqrt2;
	//				if (modk<Kmax) delta=sqrt(Pkint[l])/sqrt2; //sqrt(Pk.f(modk))/sqrt2;
	//				else {
	////					printf("zeroing power: kx, ky, kz = %lE %lE %lE\n",kx,ky,kz);
	//					delta=0;
	//				}
	//				l++;
	////				printf("delta: %lE\n",delta);
	//				fRe(i,j,k)*=delta;
	//				fIm(i,j,k)*=delta;
	////				if ((nx-i) % nx==i and (ny-j) % ny==j and (nz-k) % nz==k) { fRe(i,j,k)*=sqrt2; fIm(i,j,k)=0.0; }
	////				if ((nx-i) % nx==i and (ny-j) % ny==j and (Nz()-k) % Nz()==k) { fRe(i,j,k)*=sqrt2; fIm(i,j,k)=0.0; }
	////				if ((Nx()-i) % Nx()==i and (Ny()-j) % Ny()==j and (Nz()-k) % Nz()==k) { fRe(i,j,k)*=sqrt2; fIm(i,j,k)=0.0; }
	//				
	////				Pkintf.newPoint(modk,fRe(i,j,k)*fRe(i,j,k)+fIm(i,j,k)*fIm(i,j,k));
	//				
	//			}
	//		}
	//	}
	//	delete [] Pkint;
	//	Pkintf.save("Pkint.rand");
	
	msgs->say("Antisymmetrizing",High);
	antisymmetrize();
	
	msgs->say("backward FFT",High);
	fft(false);
	msgs->say("generateRandomGaussianField...done",High);
	
	return *this;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::generateRandomGaussianField2(mscsFunction& Pk) {
	msgs->warning("mscsFunction3dregc::generateRandomGaussianField2: THIS IMPLEMENTATION HAS NOT BEEN TESTED. PLEASE PERFORM TESTS.  In particular re2im is suspicious (see function_fft on how to do this correctly)",Top);
	freeFunctionSpace();
	
	double delta,kx,ky,kz,modk;
	double sqrt2=sqrt(double(2.0));
	long nx=Nx();
	long ny=Ny();
	long nz=Nz()/2;
	nz=Nz();
	long n=nx*ny*(nz+1);
	n=getN();
	n=nx;
	double* kk = new double[n];
	long l=0;
	double kmin,kmax;
	Pk.sortFunctionArgAscending();
	kmin=Pk.getMinArg();
	kmax=Pk.getMaxArg();
	
	msgs->say("Precalculating the k-modes and interpolating the power spectrum",High);
	for (long i = 0; i < nx; i++) {
		kx=double(i)/double(_param.sizeX);
		modk=kx;
		if (modk<kmin) { msgs->say("P(k) is not sufficiently tabulated for small k -- using the closest, minimal k value: "+msgs->toStr(kmin)+"  (value requested: "+msgs->toStr(modk)+")",High); modk=kmin; }
		if (modk>kmax) { msgs->say("P(k) is not sufficiently tabulated for large k -- using the closest, maximal k value"+msgs->toStr(kmax)+"  (value requested: "+msgs->toStr(modk)+")",High); modk=kmax; }
		
		kk[l]=modk;
		l++;
	}
	double *kin=Pk.extractArguments();
	double *Pkin=Pk.extractValues();
	double *Pkint=cpeds_interpolate(kin,Pkin,Pk.pointsCount(),kk,n,"cspline");
	//	mscsFunction Pkintf("Pkint",kk,Pkint,n);
	//	Pkintf.save("Pkint");
	//	Pk.save("Pk");
	delete [] kin;
	delete [] Pkin;
	delete [] kk;
	//	Pkintf.clearFunction();
	
	msgs->say("generating gaussian field",High);
	allocFunctionSpace();
	generateRandomGaussianField(0,1,true,true,0,0);
	
	msgs->say("Forming the power spectrum",High);
	for (long i = 0; i < nx; i++) {
		for (long j = 0; j < ny; j++) {
			l=0;
			for (long k = 0; k < nz; k++) {
				
				if (Pkint[l]<0) { 
					msgs->warning("Pkint is "+msgs->toStr(Pkint[l])+" < 0. Assuming power 0 for k-mode: (i,j,k)= ("+msgs->toStr(i)+", "+msgs->toStr(j)+", "+msgs->toStr(k)+")",High);
					Pkint[l]=0;
				}
				delta=1;
				delta*=powf(sqrt(Pkint[i])/sqrt2,1.0/3.0); 
				delta*=powf(sqrt(Pkint[j])/sqrt2,1.0/3.0); 
				delta*=powf(sqrt(Pkint[k])/sqrt2,1.0/3.0); 
				
				fRe(i,j,k)*=delta;
				fIm(i,j,k)*=delta;
				if ((nx-i) % nx==i and (ny-j) % ny==j and (Nz()-k) % Nz()==k) { fRe(i,j,k)*=sqrt2; fIm(i,j,k)=0.0; }				
			}
		}
	}
	delete [] Pkint;
	
	msgs->say("Antisymmetrizing",High);
	antisymmetrize();
	
	msgs->say("backward FFT",High);
	fft_c2c(false);
	msgs->say("generateRandomGaussianField...done",High);
	
	return *this;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::generateRandomGaussianField3(mscsFunction& Pk) {
	msgs->warning("mscsFunction3dregc::generateRandomGaussianField3: THIS IMPLEMENTATION HAS NOT BEEN TESTED. PLEASE PERFORM TESTS. In particular re2im is suspicious (see function_fft on how to do this correctly)",Top);
	freeFunctionSpace();
	
	double delta,kx,ky,kz,modk;
	double sqrt2=sqrt(double(2.0));
	long nx=Nx();
	long ny=Ny();
	long nz=Nz()/2;
	nz=Nz();
	long n=nx*ny*(nz+1);
	n=getN();
	double* kk = new double[n];
	long l=0;
	double kmin,kmax;
	Pk.sortFunctionArgAscending();
	kmin=Pk.getMinArg();
	kmax=Pk.getMaxArg();
	double Kmax=double(nx-1)/double(_param.sizeX);
	
	double p;
	
	msgs->say("generating gaussian field",High);
	allocFunctionSpace();
	for (long i = 0; i < getN(); i++) {		fRe(i)=0; fIm(i)=0;	}
	fRe(nx/2,ny/2,nz/2)=1;
	fIm(nx/2,ny/2,nz/2)=1;
	cpeds_matrix_save(getSlice(0,ny/2,0),"slice1.in");
	fft_c2c();
	
	msgs->say("Forming the power spectrum",High);
	l=0;
	for (long i = 0; i < nx; i++) {
		//		printf("%li/%li\n \%",i,nx);
		for (long j = 0; j < ny; j++) {
			//			for (long k = 0; k <= nz; k++) {
			for (long k = 0; k < nz; k++) {
				//			for (long k = 0; k < nz/2; k++) {
				
				kx=double(i)/double(_param.sizeX);
				ky=double(j)/double(_param.sizeY);
				kz=double(k)/double(_param.sizeZ);
				modk=sqrt(kx*kx + ky*ky + kz*kz);
				p=exp(-modk/0.025);
				//				if (Pkint[l]<0) { 
				//					msgs->warning("Pkint is "+msgs->toStr(Pkint[l])+" < 0. Assuming power 0 for k-mode: (i,j,k)= ("+msgs->toStr(i)+", "+msgs->toStr(j)+", "+msgs->toStr(k)+")",High);
				//					Pkint[l]=0;
				//				}
				
				//				Kmax=sqrt(pow(nx-1,2)+pow(ny-1,2)+pow(nz,2));
				delta=p;
				//				if (modk<Kmax) delta=sqrt(Pkint[l])/sqrt2; //sqrt(Pk.f(modk))/sqrt2;
				//				else {
				////					printf("zeroing power: kx, ky, kz = %lE %lE %lE\n",kx,ky,kz);
				//					delta=0;
				//				}
				//				l++;
				//				printf("delta: %lE\n",delta);
				fRe(i,j,k)*=delta;
				fIm(i,j,k)*=delta;
				//				if ((nx-i) % nx==i and (ny-j) % ny==j and (nz-k) % nz==k) { fRe(i,j,k)*=sqrt2; fIm(i,j,k)=0.0; }
				if ((nx-i) % nx==i and (ny-j) % ny==j and (Nz()-k) % Nz()==k) { fRe(i,j,k)*=sqrt2; fIm(i,j,k)=0.0; }
				//				if ((Nx()-i) % Nx()==i and (Ny()-j) % Ny()==j and (Nz()-k) % Nz()==k) { fRe(i,j,k)*=sqrt2; fIm(i,j,k)=0.0; }
				
				//				Pkintf.newPoint(modk,fRe(i,j,k)*fRe(i,j,k)+fIm(i,j,k)*fIm(i,j,k));
				
			}
		}
	}
	msgs->say("Antisymmetrizing",High);
	antisymmetrize();
	
	msgs->say("backward FFT",High);
	fft_c2c(false);
	msgs->say("generateRandomGaussianField...done",High);
	
	return *this;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::calculateGravitationalPotential(cpedsPointSet3D& ps, double resolution, double softening) {
	double xmin,xmax,ymin,ymax,zmin,zmax;
	ps.getRanges(xmin,xmax,ymin,ymax,zmin,zmax);
	long nx,ny,nz;
	double V,d;
	nx=(xmax-xmin)/resolution;
	ny=(ymax-ymin)/resolution;
	nz=(zmax-zmin)/resolution;
	cpedsPoint3D p,p2;
	
	if (nx==0) nx++;
	if (ny==0) ny++;
	if (nz==0) nz++;
	
	if (softening==-1) {
		softening=sqrt(pow((xmax-xmin),2)+pow((ymax-ymin),2)+pow((zmax-zmin),2))/ps.size();
		//		printf("gravitational softening: %lE\n",gravSoft);
	}
	//	ps.printRanges();
	//	printf("nx: %li ny: %li, nz: %li soft: %lE, size: %li\n",nx,ny,nz, softening,ps.size());
	setSizeRange(nx,ny,nz,xmin,ymin,zmin,xmax,ymax,zmax);
	allocFunctionSpace();
	setf(0,0);
	
	for (long i = 0; i < nx; i++) {
		p.setX(getX(i));
		for (long j = 0; j < ny; j++) {
			p.setY(getY(j));
			for (long k = 0; k < nz; k++) {
				p.setZ(getZ(k));				
				
				V=0;
				for (unsigned long l = 0; l < ps.size(); l++) {
					p2=ps[l];
					d=p.dist(p2);
					V-=ps.val(l)/(d+softening);
				}
				V*=CPEDS_G;
				fRe(i,j,k)=V;
			}
		}
	}
	
	
	return *this;
}
/***************************************************************************************/
//mscsFunction3dregc& mscsFunction3dregc::generateRandomGaussianField(const mscsFunction& Pk) {
//		
//	return *this;
//}

/***************************************************************************************/
mscsFunction3dregc mscsFunction3dregc::fft_r2c() const {
	double* in=exportRe();
	
	long M=Nx()*Ny()*(long(Nz()/2)+1);
	fftw_complex* out=(fftw_complex*)fftw_malloc((M)*sizeof(fftw_complex)); 
	
	fftw_plan p;
	p=fftw_plan_dft_r2c_3d(Nx(),Ny(),Nz(), in,out, FFTW_ESTIMATE); 
	fftw_execute(p);
	fftw_destroy_plan(p);
	//	fftw_free(in);
	delete [] in;
	mscsFunction3dregc tmpout(Nx(),Ny(),Nz()/2+1);
	tmpout.importFunction(out,M);
	return tmpout;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::fft_c2c(bool fwd) {
	//	fftw_complex* out=(fftw_complex*)fftw_malloc((getN())*sizeof(fftw_complex)); 
	/*
	 * Comment: REMEMBER: WHILE FFTW_ESTIMATE does not touch the _data array FFTW_MEASURE DOES !
	 * so if you ever thing of using this flag then the _data array would need to be reinitiated
	 * before fftw_execute.
	 * 
	 * author: blew
	 * date: Aug 29, 2014 1:16:53 PM
	 *
	 */
	
	fftw_plan plan=NULL;

	/*
	 * Comment: the block below does lead to some segfault errors when destroying plans in destructor
	 * of the class. I still don't know what the problem is. The temporary solution is that
	 * the plan is created and destroyed locally and it is not stored throughout the lifetime 
	 * of the object.
	 * 
	 * author: blew
	 * date: Aug 31, 2014 4:13:27 PM
	 *
	 */
	
	/*
	  //Aug 31, 2014 4:13:02 PM - DEBUG BEGIN
	  
//	  if (fwd) {
//		  if (plan_c2c_forward==NULL) {
//			  plan_c2c_forward=new fftw_plan;
//			  *plan_c2c_forward=fftw_plan_dft_3d(Nx(),Ny(),Nz(), _data,_data, FFTW_FORWARD,FFTW_ESTIMATE); 
//			  fftw_execute(*plan_c2c_forward);	
//		  }
//		  else {
//			  fftw_execute(*plan_c2c_forward);	
//		  }
//		  // divide by number of cells
//		  cpeds_divide_value(double(getN()),_data,getN());
//		  //		divide(getVolume(false));
//	  }
//	  else {
//		  if (plan_c2c_backward==NULL) {
//			  plan_c2c_backward=new fftw_plan;
//			  *plan_c2c_backward=fftw_plan_dft_3d(Nx(),Ny(),Nz(), _data,_data, FFTW_BACKWARD,FFTW_ESTIMATE); 
//			  fftw_execute(*plan_c2c_backward);	
//		  }
//		  else {
//			  fftw_execute(*plan_c2c_backward);	
//		  }
//	  }
	  
	  //DEBUG - END
	*/

	if (fwd) {
		if (lengthZ()!=0) 
			plan=fftw_plan_dft_3d(Nx(),Ny(),Nz(), _data,_data, FFTW_FORWARD,FFTW_ESTIMATE); 
		else
			plan=fftw_plan_dft_2d(Nx(),Ny(), _data,_data, FFTW_FORWARD,FFTW_ESTIMATE); 
		fftw_execute(plan);	
		// divide by number of cells
		cpeds_divide_value(sqrt(double(getN())),_data,getN());
		//		divide(getVolume(false));
	}
	else {
		if (lengthZ()!=0) 
			plan=fftw_plan_dft_3d(Nx(),Ny(),Nz(), _data,_data, FFTW_BACKWARD,FFTW_ESTIMATE); 
		else
			plan=fftw_plan_dft_2d(Nx(),Ny(), _data,_data, FFTW_BACKWARD,FFTW_ESTIMATE); 
		fftw_execute(plan);	
		cpeds_divide_value(sqrt(double(getN())),_data,getN());
	}

	fftw_destroy_plan(plan);
	fftw_cleanup();
	
	
	return *this;
}
/***************************************************************************************/
//mscsFunction3dregc mscsFunction3dregc::slowft_r2c(const cpedsList<double>& kx, const cpedsList<double>& ky, const cpedsList<double>& kz) {
////	fftw_complex* in=(fftw_complex*)fftw_malloc((getN())*sizeof(fftw_complex)); 
//	long nx=Nx();
//	long ny=Ny();
//	long nz=Nz();
//	long nkx=kx.size();
//	long nky=ky.size();
//	long nkz=kz.size();
//	mscsFunction3dregc g("FT");
//	double re,im;
//	
//	for (long ki = 0; ki < nkx; ki++) {
//		for (long kj = 0; kj < nky; kj++) {
//			for (long kk = 0; kk < nkz; kk++) {
//				re=im=0;
//				
//				for (long i = 0; i < nx; i++) {
//					for (long j = 0; j < ny; j++) {
//						for (long k = 0; k < nz; k++) {
//							re+= fRe(i,j,k)*cos(   twoPI * ( kx[ki]*getX(i) + ky[kj]*getY(j) + kz[kk]*getZ(k) )   );
//							im+=-fIm(i,j,k)*sin(   twoPI * ( kx[ki]*getX(i) + ky[kj]*getY(j) + kz[kk]*getZ(k) )   );
//						}
//					}
//				}
//				g.newPoint(re,im);
//			}
//		}
//	}
//	
//	g._X=kx;
//	g._Y=ky;
//	g._Z=kz;
//	g.setSize(nx,ny,nz);
//	
//	return g;
//}

/***************************************************************************************/
mscsFunction mscsFunction3dregc::powerSpectrum(double kRes) {
	mscsFunction Pk("power spectrum");
	mscsVector<double> kvec;
	mscsVector<double> Nk;
	mscsVector<double> P;
	mscsVector<double> kk;
	
	// derive the fourier modes
	//	fft_r2c(&dk);
	fft_c2c(true);
//	multiply(double(getN()));
	// derive the power spectrum
	
	long nx=Nx();
	long ny=Ny();
	long nz=Nz();

	if (lengthZ()!=0) { // 3D case
		nz=Nz()/2+1;
	}
	else { // 2D case
		ny=Ny()/2+1;
		assert(nz==1);
	}
	
	printf("nz %li ny %li nz %li\n",nx,ny,nz);
	mscsVector<double> kx,ky,kz;
	double modk,modi;
	double kxMin=1.0/_param.sizeX;
	double kxMax=double((nx/2+1))/_param.sizeX;

	//
	// calculate power P(vec k) 
	//

	// pre-calculate kx, ky, kz
	kx.setSize(nx);
	ky.setSize(ny);
	kz.setSize(nz);
	
	// define the k-space of the output P(k)
	for (long i = 0; i < nx; i++) { kx[i]=double(i)/double(_param.sizeX); }
	for (long j = 0; j < ny; j++) { ky[j]=double(j)/double(_param.sizeY); }
	if (lengthZ()==0) kz[0]=0;
	else
		for (long k = 0; k < nz; k++) { kz[k]=double(k)/double(_param.sizeZ); }

	// calculate the size of the output P(k)
	long l1=0;
//	long l2=0;
	long m=0;
	double tmp;
//	long sizek=kAverageAcc/kxMin;
	long sizek=Nx()/kRes;

	printf("kmin: %lE, kmax: %lE, dk: %lE\n",kxMin,kxMax,kRes*kxMin);
	printf("sizek: %li\n",sizek);
	// allocate and initialize with 0
	Nk.setSize(sizek); 
	P.setSize(sizek);
	kk.setSize(sizek);
//	for (unsigned long i = 0; i < sizek; i++) {
//		Nk[i]=0;
//		P[i]=0;
//		kk[i]=0;
//	}
	
	//
	// calculate power P(k) 
	//
	for (long i = 0; i < nx; i++) {
		for (long j = 0; j < ny; j++) {
			for (long k = 0; k < nz; k++) {
				
				l1=ijk2idx(i,j,k);
//				l2=kji2idx((nx-i) % nx,(ny-j) % ny,(nz-k) % nz);
//				if (l1<=l2) {
					// derive the power
					tmp=(_data[l1][0]*_data[l1][0] + _data[l1][1]*_data[l1][1]); //tmp*=tmp;
					if (_param.sizeZ==0) {
						modk=sqrt(kx[i]*kx[i] + ky[j]*ky[j]);
						modi=sqrt(i*i + j*j);
					}
					else {
						modi=sqrt(double(i*i + j*j + k*k));
						modk=sqrt(kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k]);
					}
					m=modi/_param.Nx*sizek;
					if (m<sizek) { 
						kk[m]+=modk;
						P[m]+=tmp;
						Nk[m]+=1.0;
					}
//				}
			}
		}
	}
	
	
////	long i;
////	printf("starting parallel for\n");
////#pragma omp parallel for shared(kk,P,Nk) firstprivate(nx,ny,nz,kx,ky,kz,sizek) private(l,m,tmp,modk)
//	for (long k = 0; k < nz; k++) {
//		if (k==0 or k==Nz()/2) { nx=Nx()/2; ny=Ny()/2; } else { nx=Nx(); ny=Ny(); }
//		for (long i = 0; i < nx; i++) {
//			for (long j = 0; j < ny; j++) {
//				l=ijk2idx(i,j,k);
//				// derive the power
//				tmp=(_data[l][0]*_data[l][0] + _data[l][1]*_data[l][1]); //tmp*=tmp;
//				
////				if (_param.sizeZ==0) modk=sqrt(kx[i]*kx[i] + ky[j]*ky[j]);
////				else modk=sqrt(kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k]);
////				m=modk/kxMax*sizek;
////				if (m<sizek) { // shpere mode
////					kk[m]+=modk;
////					P[m]+=tmp;
////					Nk[m]+=1.0;
////				}
//				if (_param.sizeZ==0) {
//					modk=sqrt(kx[i]*kx[i] + ky[j]*ky[j]);
//					modi=sqrt(i*i + j*j);
//				}
//				else {
//					modi=sqrt(double(i*i + j*j + k*k));
//					modk=sqrt(kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k]);
//				}
//				m=modi/_param.Nx*sizek;
//				if (m<sizek) { 
//					kk[m]+=modk;
//					P[m]+=tmp;
//					Nk[m]+=1.0;
//				}
//			}
//		}
//	}
//	printf("end of parallel block\n");
	
	//
	// transcript the results onto output function
	//
	for (unsigned long i = 0; i < sizek; i++) {
		if (Nk[i]>0) {
			Pk.newPoint(kk[i]/Nk[i],P[i]/Nk[i]);
		}
	}
	Pk.sortFunctionArgAscending();
//	cpeds_save_matrix(Nk.toArray(),Nk.size(),1,"Nk",true,false);
//	cpeds_save_matrix(kk.toArray(),Nk.size(),1,"k",true,false);
//	cpeds_save_matrix(P.toArray(),Nk.size(),1,"P",true,false);
//	Pk*=double(getN());
	
//	//
//	// calculate power P(vec k) 
//	//
//	long l=0;
//	double tmp;
//	//	Pveck.setSize(nx*ny*nz);
//	Pveck.setSize(getN());
//	for (long i = 0; i < nx; i++) {
//		for (long j = 0; j < ny; j++) {
//			for (long k = 0; k < nz; k++) {
//				l=ijk2idx(i,j,k);
//				// derive the power
//				tmp=(_data[l][0]*_data[l][0] + _data[l][1]*_data[l][1]); //tmp*=tmp;
//				//				Pveck.append(tmp);
//				Pveck[l]=tmp;
//				//				l++;
//			}
//		}
//	}
//	//
//	// calculate power P(k) = integrate sqrt(kx^2+ky^2+kz^2) f(vec k) dkx dky dkz
//	//
//	l=0;
//	for (long i = 0; i < nx; i++) {
//		for (long j = 0; j < ny; j++) {
//			for (long k = 0; k < nz; k++) {
//				// derive k=sqrt(kx^2+ky^2+kz^2)
//				kx=double(i)/double(_param.sizeX);
//				ky=double(j)/double(_param.sizeY);
//				kz=double(k)/double(_param.sizeZ);
//				
//				if (_param.sizeZ==0) modk=sqrt(kx*kx + ky*ky);
//				else modk=sqrt(kx*kx + ky*ky + kz*kz);
//
//				if (modk<=kxMax) {
//					l=ijk2idx(i,j,k);
//					Pk.newPoint(modk,Pveck[l]);
//				}
//				//				l++;
//			}
//		}
//	}
//	//	Pk/=double(_param.sizeX*_param.sizeY*_param.sizeZ); // BLcomment (Jun 25, 2013, 6:28:21 PM): commented out since this division is done in fft_c2c(true)
//	msgs->say("sorting Pk",Low);
//	msgs->say("Pk has %li points",Pk.pointsCount(),Low);
//	msgs->say("kx_max = %lE",kxMax,Low);
//	Pk.sortFunctionArgAscending();
//	msgs->say("averaging within same k values",Low);
//	Pk.average_sameArgs(kAverageAcc);
	return Pk;	
}
/***************************************************************************************/
mscsFunction mscsFunction3dregc::correlationFunctionR(double Lmin, double Lmax, double dr) {
	mscsFunction corr;
	double sep=0;
	double mu,sigma2;
	double x1,y1,z1,x2,y2,z2;
	long di,dj,dk;
	long dx,dy,dz;
	
	long nC=0;
	long l;
	if (Lmin==-1) Lmin=getDx();
	if (Lmax==-1) Lmax=lengthX();
	if (dr==-1) dr=getDx();
	
	nC=(Lmax-Lmin)/dr;
	mscsVector<double> C;
	mscsVector<double> dist;
	mscsVector<double> num;
	
	C.setSize(nC);
	dist.setSize(nC);
	num.setSize(nC);

	mu=averageRe();
	mu=0;
	sigma2=varianceRe();
	
	for (long i1 = 0; i1 < Nx(); i1++) {
		printf("%li/%li\n",i1,Nx());
		for (long j1 = 0; j1 < Ny(); j1++) {
			for (long k1 = 0; k1 < Nz(); k1++) {

//	for (unsigned long i = 0; i < getN(); i++) {
		
				
				for (long i2 = i1; i2 < Nx(); i2++) {
					if (i1!=i2) {
						di=i2-i1;
						dx=di*_param.dx;
						
						for (long j2 = j1; j2 < Ny(); j2++) {
							if (j1!=j2) {
								dj=j2-j1;
								dy=dj*_param.dy;
								for (long k2 = k1; k2 < Nz(); k2++) {
									if (k1!=k2) {
										dk=k2-k1;
										dz=dk*_param.dz;
										sep=sqrt(dx*dx+dy*dy+dz*dz);

										if (sep>Lmin and sep<=Lmax) {
											l=(sep-Lmin)/(Lmax-Lmin)*nC;
											if (l>=0 and l<nC) {
												num[l]++;
												C[l]+=(fRe(i1,j1,k1)-mu)*(fRe(i2,j2,k2)-mu);
												dist[l]+=sep;
											}
										}
									}
									
								}
							}
						}
					}
				}
				
				
			}
		}
	}
	
	for (unsigned long i = 0; i < nC; i++) {
		if (num[i]>0) {
			corr.newPoint(dist[i]/num[i], C[i]/num[i]/sigma2);
		}
	}
	corr.sortFunctionArgAscending();
	
	return corr;
}
/***************************************************************************************/
//mscsFunction3dregc& mscsFunction3dregc::mkDensityField(subDomain_region_t r, const mscsVector<cpedsPoint3D>& positions, const mscsVector<cpedsPoint3D>* weights, int part, mscsFunction& smoothingKernel) {
mscsFunction3dregc& mscsFunction3dregc::mkDensityFieldGather(subDomain_region_t r, subDomain_region_t treeScheme, const mscsVector<cpedsPoint3D>& positions, const mscsVector<double>* weights, int part, string smKernel, long NeighborsMin) {
	double dx,dy,dz;
	
	// define the function space
	dx=(r.xmax-r.xmin)/r.subx;
	dy=(r.ymax-r.ymin)/r.suby;
	dz=(r.zmax-r.zmin)/r.subz;
	setSize(r.subx,r.suby,r.subz,dx,dy,dz,r.xmin,r.ymin,r.zmin);
	allocFunctionSpace();
	
	// define the smoothing kernel - this kernel is good for any smoothing length as it is not normalized (it is not in inverse volume units)
	// distance is in R=r/hsml
	double (*kernel)(double );
	double norm;
	//	mscsWindowFunction wfn;
	if (smKernel == "gadget2") { 
		kernel=&mscsWindowFunction::kernelGadget; 	
		if (dz==0) norm=double(40.0)/(7.0*PI); //2d case
		else norm=8.0/PI; //3d case
	}
	else {
		if (smKernel == "gadget2b") { 
			kernel=&mscsWindowFunction::kernelGadget2b; 	
			if (dz==0) norm=double(40.0)/(7.0*PI); //2d case
			else norm=8.0/PI; //3d case
		}
		else {
			msgs->criticalError("mscsFunction3dregc::mkDensityField >> don't know this smoothing kernel function: "+smKernel,High);
		}
	}
	
	//	mscsWindowFunction w("SPH kernel");
	//	w.mkSPHkernelGadget2(100);
	//	double *wd=w.extractValues();
	//	double *rd=w.extractArguments();
	
	// prepare tree
	subDomain D(&positions,&treeScheme,1L); D.tree();

	double hsml, rho,W,dist;
//	static mscsVector<long> *neighbors;
//#pragma omp threadprivate(neighbors)
	mscsVector<long> *neighbors;
	neighbors = new mscsVector<long>;
	
	long acc=13;
	//	if (r.zmax!=r.zmin) acc=pow((r.xmax-r.xmin)*(r.ymax-r.ymin)*(r.zmax-r.zmin)/positions.size(),0.3333333)/10;
	//	else acc=sqrt((r.xmax-r.xmin)*(r.ymax-r.ymin)/positions.size())/10;
	printf("acc: %li\n",acc);
	hsml=sqrt(dx*dx+dy*dy+dz*dz);
	double dV;
	// do the density calculation
	double x,y,z;
	for (long i = 0; i < r.subx; i++) {
		x=getX(i);
		for (long j = 0; j < r.suby; j++) {
			y=getY(j);
			for (long k = 0; k < r.subz; k++) {
				z=getZ(k);
				
				// get particles within the distance of R defining region containing constant number of particles
				printf("--------------------\ngetting neighbors\n");
				(*neighbors)=D.getNeighbors(cpedsPoint3D(x,y,z),NeighborsMin,NeighborsMin, &hsml,acc);
				printf("ijk: %li %li %li,   x: %lf, y: %lf z: %lf, Nneigh: %li, hsml: %lf\n",i,j,k, x,y,z,neighbors->size(),hsml);
				//				neighbors.print();
				rho=0;
				// define hsml as the distance to the maximally distant point -- this in principle excludes this most distance point from the density calculations
				// because the kernel becomes 0 at r/h=1
				hsml=0;
				for (long l = 0; l < neighbors->size(); l++) {
					dist=cpedsPoint3D(x,y,z).dist(positions[(*neighbors)[l]]);
					if (dist>hsml) hsml=dist;
				}
				if (hsml==0) { msgs->criticalError("mkDensityField>>> hsml is 0. Cannot continue",Top); }
				printf("calculating density: hsml=%lf n=%li\n\n",hsml,neighbors->size());
				// calculate density
				W=0;
				if (weights!=NULL) {
					for (long l = 0; l < neighbors->size(); l++) {
						dist=cpedsPoint3D(x,y,z).dist(positions[(*neighbors)[l]]);
						rho+=weights->at((*neighbors)[l])*kernel(dist/hsml);
					}
				}
				else {
					for (long l = 0; l < neighbors->size(); l++) {
						dist=cpedsPoint3D(x,y,z).dist(positions[(*neighbors)[l]]);
						rho+=(*kernel)(dist/(hsml));
						//						if (i==20) { printf("rho: %lf\n",rho); }
					}			
					//					if (hsml==0) { printf("**************hsml=0\n"); }
				}
				dV=hsml*hsml;
				if (dz!=0) dV*=hsml;
				f(i,j,k)[part]=norm*rho/dV;
				//				f(i,j,k)[part]=rho;
				f(i,j,k)[1]=hsml;
				//				f(i,j,k)[1]=neighbors.size();
				//				if (cpeds_isnan(f(i,j,k)[part])) printf("********************we have NAN\n");
			}
		}
	}
	
	delete neighbors;
	
	return *this;
}

/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::mkDensityFieldScatter(subDomain_region_t r, subDomain_region_t treeScheme, const mscsVector<cpedsPoint3D>& positions, const mscsVector<double>* weights, int part, string smKernel, long NeighborsMin) {
	double dx,dy,dz;
	
	// define the function space
	dx=(r.xmax-r.xmin)/r.subx;
	dy=(r.ymax-r.ymin)/r.suby;
	dz=(r.zmax-r.zmin)/r.subz;
	setSize(r.subx,r.suby,r.subz,dx,dy,dz,r.xmin,r.ymin,r.zmin);
	allocFunctionSpace();
	
	// define the smoothing kernel - this kernel is good for any smoothing length as it is not normalized (it is not in inverse volume units)
	// distance is in R=r/hsml
	double (*kernel)(double );
	double norm;
	//	mscsWindowFunction wfn;
	if (smKernel == "gadget2") { 
		kernel=&mscsWindowFunction::kernelGadget; 	
	//		if (dz==0) norm=double(40.0)/(7.0*PI); //2d case
		if (r.subz==1) norm=double(40.0)/(7.0*PI); //2d case
		else norm=8.0/PI; //3d case
	}
	else {
		if (smKernel == "gadget2b") { 
			kernel=&mscsWindowFunction::kernelGadget2b; 	
			if (dz==0) norm=double(40.0)/(7.0*PI); //2d case
			else norm=8.0/PI; //3d case
		}
		else {
			msgs->criticalError("mscsFunction3dregc::mkDensityField >> don't know this smoothing kernel function: "+smKernel,High);
		}
	}
	
	// prepare tree
	msgs->say("building tree",Medium);
	subDomain D(&positions,&treeScheme,1L); D.tree();
	//	D.saveDomains("domains");
	msgs->say("done",Low);
	
	//	{
	//	subDomain* E= new subDomain(&positions,&treeScheme,1L); E->tree();
	//	E->saveDomains("orig");
	//	msgs->say("copying building tree",Medium);
	//	subDomain DD(*E); 
	//	delete E;
	//	DD.saveDomains("orig-copy");
	//	msgs->say("done",Low);
	//	}
	//exit(0);
	
	treeScheme.zmin=-1; 	treeScheme.zmax=1;
	D.print_domain_range(treeScheme);
	
	double hsml, rho,dist;
	mscsVector<long> neighbors;
	mscsVector<double> sml;
	long acc=13;
	//	if (r.zmax!=r.zmin) acc=pow(dx*dy*dz,0.3333333)/10;
	if (r.subz>1) { 
		//		acc=pow(dx*dy*dz,0.3333333)/10;
		hsml=sqrt(dx*dx+dy*dy+dz*dz);
	}
	else { 
		//		acc=sqrt(dx*dy)/10;
		hsml=sqrt(dx*dx+dy*dy);
	}
	//	printf("acc: %lE\n",acc);
	//	printf("hsml: %lE\n",hsml);
	double dr2d=sqrt(dx*dx+dy*dy);
	double dr3d=sqrt(dx*dx+dy*dy+dz*dz);
	double dr;
	if (r.subz==1) dr=dr2d; else dr=dr3d;
	//
	// calculate smoothing lengths for all particles
	//
	msgs->say("calculating smoothing lengths",Medium);
	
	long N=positions.size();
	sml.resize(N,0.0);
	double hsmlMax=0;
	long i;
#pragma omp parallel for firstprivate(D,hsml,acc,NeighborsMin,N) shared(positions, sml) private(i)
	for (i = 0; i < N; i++) {
		printf("looking for neighbors around point %li\n",i);
		//		positions[i].print_point();
		//		neighbors=D.getNeighbors(positions[i],NeighborsMin, &hsml,acc);
		D.getNeighbors(positions[i],NeighborsMin,NeighborsMin, &hsml,acc).size();
		
		//		int n=D.getNeighbors(positions[i],NeighborsMin, &(sml[i]),acc).size();
		
		//		printf("aaa %li %lf %i\n",i,hsml,n);
		//		printf("aaa %li %lf\n",i,hsml);
		if (hsml==0) hsml=10*dr; // this is only for NeighborsMin=1 case. TODO: This can be improved 
		//		if (sml[i]==0) sml[i]=20*sqrt(dx*dx+dy*dy+dz*dz); // this is only for NeighborsMin=1 case. TODO: This can be improved 
		//		if (i==40) printf("hsml: %lf\n",hsml);
		sml[i]=hsml;
		//		sml[i]=positions[i].x();
		//		d[i]=hsml;
		//		if (hsml>hsmlMax) hsmlMax=hsml;
		//		printf("particle: %li hsml: %lf, n=%li\n",i, hsml,neighbors.size());
		//		positions[i].print_point();
	}
	for (long i = 0; i < N; i++) {	if (sml[i]>hsmlMax) hsmlMax=sml[i];	}
	//	cpedsList<double> tmp(sml); tmp.save("testSML");
	cpedsList<double> tmp(sml); tmp.save("testSML-para");
	msgs->say("done",Low);
	//	exit(0);
	
	//
	// do the density calculation
	//
	msgs->say("calculating density",Medium);
	double x,y,z;
	for (long i = 0; i < r.subx; i++) {
		msgs->say("done: "+msgs->toStr(double(i)/r.subx*100)+" %",Low);
		x=getX(i);
		//#pragma omp parallel for firstprivate(D,hsml,neighbors,x,y,z,i,dist) reduction(+: rho)
		//#pragma omp parallel for firstprivate(D,hsml,x,i,acc,NeighborsMin) private(y,z,dist,neighbors) reduction(+: rho)
		for (long j = 0; j < r.suby; j++) {
			y=getY(j);
			for (long k = 0; k < r.subz; k++) {
				z=getZ(k);
				
				// get particles within the distance of R defining region containing constant number of particles
				//				printf("--------------------\ngetting neighbors\n");
				neighbors=D.getNeighbors(cpedsPoint3D(x,y,z),NeighborsMin,NeighborsMin, &hsml,acc);
				//				printf("ijk: %li %li %li,   x: %lf, y: %lf z: %lf, Nneigh: %li, hsml: %lf\n",i,j,k, x,y,z,neighbors.size(),hsml);
				if (hsml==0) { msgs->criticalError("mkDensityField>>> hsml is 0. Cannot continue",Top); }
				//				neighbors.print();
				rho=0;
				// define hsml as the distance to the maximally distant point -- this in principle excludes this most distance point from the density calculations
				// because the kernel becomes 0 at r/h=1
				//				printf("calculating density: hsml=%lf n=%li\n\n",hsml,neighbors.size());
				// calculate density
				if (dz==0) {
					if (weights!=NULL) {
						for (long l = 0; l < neighbors.size(); l++) {
							dist=cpedsPoint3D(x,y,z).dist(positions[neighbors[l]]);
							rho+=weights->at(neighbors[l])*kernel(dist/sml[neighbors[l]])/(sml[neighbors[l]]*sml[neighbors[l]]);
						}
					}
					else {
						//#pragma omp for 
						for (long l = 0; l < neighbors.size(); l++) {
							dist=cpedsPoint3D(x,y,z).dist(positions[neighbors[l]]);
							rho+=kernel(dist/(sml[neighbors[l]]))/(sml[neighbors[l]]*sml[neighbors[l]]);
						}			
					}					
				}
				else {
					if (weights!=NULL) {
						for (long l = 0; l < neighbors.size(); l++) {
							dist=cpedsPoint3D(x,y,z).dist(positions[neighbors[l]]);
							rho+=weights->at(neighbors[l])*kernel(dist/sml[neighbors[l]])/(sml[neighbors[l]]*sml[neighbors[l]]*sml[neighbors[l]]);
						}
					}
					else {
						for (long l = 0; l < neighbors.size(); l++) {
							dist=cpedsPoint3D(x,y,z).dist(positions[neighbors[l]]);
							rho+=kernel(dist/(sml[neighbors[l]]))/(sml[neighbors[l]]*sml[neighbors[l]]*sml[neighbors[l]]);
						}
					}					
				}
				f(i,j,k)[part]=norm*rho;
				//				f(i,j,k)[part]=rho;
				f(i,j,k)[1]=hsml;
				//				f(i,j,k)[1]=neighbors.size();
				//				if (cpeds_isnan(f(i,j,k)[part])) printf("********************we have NAN\n");
			}
		}
	}
	msgs->say("calculating density done",Low);
	
	
	
#ifdef DEBUG
	// mass test
	{
		double M=0;
		dx=getDx();
		dy=getDy();
		dz=getDz();
		D.print_domain_range(r);
		double dV;
		if (dz==0) dV=dx*dy; else dV=dx*dy*dz;
		for (long i = 0; i < r.subx; i++) {
			x=getX(i);
			for (long j = 0; j < r.suby; j++) {
				y=getY(j);
				for (long k = 0; k < r.subz; k++) {
					z=getZ(k);
					M+=f(i,j,k)[part]*dV;
				}
			}
		}
		printf("test: M=%lf, dV=%lE\n",M,dV);
	}
#endif
	return *this;
}

/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::mkDensityFieldScatter2(subDomain_region_t r, subDomain_region_t treeScheme, const mscsVector<cpedsPoint3D>& positions, const mscsVector<double>* weights, string smKernel, long NeighborsMin, long NeighborsMax) {
	double dx,dy,dz;
	
	if (positions.size()<NeighborsMin) {
		NeighborsMin=positions.size();
		msgs->warning("The provided number of particles is smaller than the requested minimal number of neighbors. Decreasing NeighborsMin to the actual number of particles.", Top);
	}
	if (positions.size()<NeighborsMax) {
		NeighborsMax=positions.size();
		msgs->warning("The provided number of particles is smaller than the requested maximal number of neighbors. Decreasing NeighborsMax to the actual number of particles.", Top);
	}
	if (NeighborsMax<NeighborsMin) {
		NeighborsMax=NeighborsMin;
		msgs->warning("The maximal number of neighbors is smaller than the minimal number of neighbors. Will make them equal.", High);
	}
	
	//
	// define the function space
	//
	dx=(r.xmax-r.xmin)/r.subx;
	dy=(r.ymax-r.ymin)/r.suby;
	dz=(r.zmax-r.zmin)/r.subz;
	setSize(r.subx,r.suby,r.subz,dx,dy,dz,r.xmin,r.ymin,r.zmin);
	allocFunctionSpace();
	//	printInfo();
	
	//
	// define the smoothing kernel - this kernel is good for any smoothing length as it is not normalized (it is not in inverse volume units)
	// distance is in R=r/hsml
	double (*kernel)(double );
	double norm;
	if (smKernel == "gadget2") { 
		kernel=&mscsWindowFunction::kernelGadget; 	
		if (dz==0) norm=double(40.0)/(7.0*PI); //2d case
		else norm=8.0/PI; //3d case
	}
	else {
		if (smKernel == "gadget2b") { 
			kernel=&mscsWindowFunction::kernelGadget2b; 	
			if (dz==0) norm=double(40.0)/(7.0*PI); //2d case
			else norm=8.0/PI; //3d case
		}
		else {
			msgs->criticalError("mscsFunction3dregc::mkDensityField >> don't know this smoothing kernel function: "+smKernel,High);
		}
	}
	
	//
	// prepare tree
	//
	msgs->say("building tree",Medium);
	subDomain D(&positions,&treeScheme,NeighborsMin); D.tree();
	//	D.saveDomains("domains");
	msgs->say("done",Low);
	//	treeScheme.zmin=-1; 	treeScheme.zmax=1;
	//	D.print_domain_range(treeScheme);
	
	//
	// prepare some important variables
	//
	
	double hsml, rho,dist;
	mscsVector<long> neighbors;
	mscsVector<double> sml;
	long acc=13;
	if (dz!=0) { // 3D case
		hsml=sqrt(dx*dx+dy*dy+dz*dz);
	}
	else { // 2D case
		hsml=sqrt(dx*dx+dy*dy);
	}
	//	printf("acc: %li\n",acc);
	//	printf("hsml: %lE\n",hsml);
	double dr2d=sqrt(dx*dx+dy*dy);
	double dr3d=sqrt(dx*dx+dy*dy+dz*dz);
	double dr;
	if (dz==0) dr=dr2d; else dr=dr3d;
	
	//
	// calculate smoothing lengths for all particles
	//
	msgs->say("calculating smoothing lengths",Medium);
	
	long N=positions.size();
	sml.resize(N,0.0);
	double hsmlMax=0;
	long i;
	//	printf("looking for neighbors around points \n");
#ifdef DEBUG_DENSITY
	
#pragma omp parallel for firstprivate(D,hsml,acc,NeighborsMin,NeighborsMax,N) shared(positions, sml) private(i)
	for (i = 0; i < N; i++) {
		if (i%1000==0) printf("done: %lf\%\n",float(i)/N*100);
		//		positions[i].print_point();
		//		neighbors=D.getNeighbors(positions[i],NeighborsMin, &hsml,acc);
		D.getNeighbors(positions[i],NeighborsMin, NeighborsMax, &hsml,acc).size();
		
		//		int n=D.getNeighbors(positions[i],NeighborsMin, &(sml[i]),acc).size();
		
		//		printf("aaa %li %lf %i\n",i,hsml,n);
		//		printf("aaa %li %lf\n",i,hsml);
		if (hsml==0) hsml=10*dr; // this is only for NeighborsMin=1 case. TODO: This can be improved 
		//		if (sml[i]==0) sml[i]=20*sqrt(dx*dx+dy*dy+dz*dz); // this is only for NeighborsMin=1 case. TODO: This can be improved 
		//		if (i==40) printf("hsml: %lf\n",hsml);
		sml[i]=hsml;
		//		sml[i]=positions[i].x();
		//		d[i]=hsml;
		//		if (hsml>hsmlMax) hsmlMax=hsml;
		//		printf("particle: %li hsml: %lf, n=%li\n",i, hsml,neighbors.size());
		//		positions[i].print_point();
	}
	
#else
	
	/*
	 * Comment: this block intends to optimize smoothing length calculations by (possibly) more efficient walk from one particle to another by exploiting 
	 * the information from the tree. The idea is that it is possible that the initial guess value for the final hsml given to the getNeighbors
	 * method will be pretty much the same for the particles spatially located close to each other. This is why intuitively it is better to walk
	 * over the idenxes of particles from the same leaf domain, and then to move to particles from another leaf domain.
	 * 
	 * The comparizon test with the regular iterative walk shows consistenty at the level 1e-8 in the final density computations for NeighborsMin==NeighborsMax, but
	 * it also shows that the discrepancy can be as large as O(10) % when NeighborsMin=30 NeighborsMax=33. This is because it is possible that different particles 
	 * will get slightly different smoothing lengths.
	 * 
	 * In the case of NeighborsMin==NeighborsMax the results are not exactly the same most likely due in cases when it is not possible to find exactly the requested number
	 * of neighbors. In such cases the getNeighbors routine returns the biggest hsml found that yields condition Nneighbors>NeighborsMin
	 * 
	 * The previous block is now available with DEBUG_DENSITY compiler flag - but since this one gives compatible results the other one will 
	 * be depreciated
	 * 
	 * COMMENT Aug 23, 2013, 2:05:57 PM - passing D as first private is somewhat ineffective because it leads to copying omp_get_thread_num() times this structure
	 * and since the structure can be very big - it is very memory consuming and also time consuming.
	 * 
	 * author: blew
	 * date: Feb 25, 2013 8:28:49 PM
	 *
	 */
	
	vector< vector<double> > v;
	v=D.getBranchParticlesCount(true);
	msgs->say("number of leaf domains: %li\n", long(v.size()),Low);
	
	long idx,k,pdone=0;
	
#pragma omp parallel for firstprivate(D,hsml,acc,NeighborsMin,NeighborsMax,v) shared(positions, sml) private(i,k) reduction(+:pdone)
	for (k = 0; k < v.size(); k++) {
		pdone=pdone+1;
#ifdef HAVE_OPENMP
		if (omp_get_thread_num()==0)
			msgs->say("done: %lf",double(double(pdone)/v.size()*omp_get_num_threads()*100),Low);
#else
		if (pdone%100==0) msgs->say("done: %lf\%\n",double(double(pdone)/v.size()*100),Low);
#endif
		for (long j = 8; j < v[k].size(); j++) {
			i=v[k][j];
			
			D.getNeighbors(positions[i],NeighborsMin, NeighborsMax, &hsml,acc);
			
			if (hsml==0) hsml=10*dr; // this is only for NeighborsMin=1 case. TODO: This can be improved 
			sml[i]=hsml;
		}
		
	}
	
#endif
	
	//
	// find the biggest hsml for density computations
	//
	hsmlMax=0;
	for (long i = 0; i < N; i++) {	if (sml[i]>hsmlMax) hsmlMax=sml[i];	}
	msgs->say("hsmlMax: %lE, hsmlMax/boxSize: %lE\n",hsmlMax,hsmlMax/(r.xmax-r.xmin),Low);
	//	cpedsList<double> tmp(sml); tmp.save("testSML");
	//	cpedsList<double> tmp(sml); tmp.save("testSML-para");
	msgs->say("done",Low);
	//	exit(0);
	
	//
	// do the density calculation
	//
	msgs->say("calculating density",Medium);
	double x,y,z,Neff;
	
	/*
	 * Comment: the below density calculation block was tested for giving the same results as the serial implementation, so now
	 * it is a standard implementation.
	 * 
	 * author: blew
	 * date: Feb 26, 2013 1:12:15 AM
	 *
	 */
	
#pragma omp parallel for firstprivate(D,hsml,x,y,z,hsmlMax, neighbors,r,rho,Neff,dist,norm) shared(positions, sml) 
	for (long i = 0; i < r.subx; i++) {
#ifdef HAVE_OPENMP
		if (omp_get_thread_num()==0) 
			msgs->say("done: "+msgs->toStr(double(i)*omp_get_num_threads()/r.subx*100)+" %",Low);
#else
		msgs->say("done: %lf\%",double(i)/r.subx*100,Low);
#endif
		x=getX(i);
		for (long j = 0; j < r.suby; j++) {
			y=getY(j);
			for (long k = 0; k < r.subz; k++) {
				z=getZ(k);
				
				// get particles within the distance of R defining region containing constant number of particles
				neighbors=D.getPointsIdx(cpedsPoint3D(x,y,z),hsmlMax,NULL);
				
				//				printf("ijk: %li %li %li,   x: %lf, y: %lf z: %lf, Nneigh: %li, hsml: %lf\n",i,j,k, x,y,z,neighbors.size(),hsml);
				//				neighbors.print();
				rho=0;
				// define hsml as the distance to the maximally distant point -- this in principle excludes this most distance point from the density calculations
				// because the kernel becomes 0 at r/h=1
				//				printf("calculating density: hsml=%lf n=%li\n\n",hsml,neighbors.size());
				// calculate density
				Neff=0;
				if (dz==0) {
					if (weights!=NULL) {
						for (long l = 0; l < neighbors.size(); l++) {
							dist=cpedsPoint3D(x,y,z).dist(positions[neighbors[l]]);
							if (kernel(dist/sml[neighbors[l]])>0) Neff++; // debug
							rho+=weights->at(neighbors[l])*kernel(dist/sml[neighbors[l]])/(sml[neighbors[l]]*sml[neighbors[l]]);
						}
					}
					else {
						for (long l = 0; l < neighbors.size(); l++) {
							dist=cpedsPoint3D(x,y,z).dist(positions[neighbors[l]]);
							rho+=kernel(dist/(sml[neighbors[l]]))/(sml[neighbors[l]]*sml[neighbors[l]]);
						}			
					}					
				}
				else {
					if (weights!=NULL) {
						for (long l = 0; l < neighbors.size(); l++) {
							dist=cpedsPoint3D(x,y,z).dist(positions[neighbors[l]]);
							rho+=weights->at(neighbors[l])*kernel(dist/sml[neighbors[l]])/(sml[neighbors[l]]*sml[neighbors[l]]*sml[neighbors[l]]);
						}
					}
					else {
						for (long l = 0; l < neighbors.size(); l++) {
							dist=cpedsPoint3D(x,y,z).dist(positions[neighbors[l]]);
							rho+=kernel(dist/(sml[neighbors[l]]))/(sml[neighbors[l]]*sml[neighbors[l]]*sml[neighbors[l]]);
						}
					}					
				}
				f(i,j,k)[0]=norm*rho;
				f(i,j,k)[1]=Neff;
			}
		}
	}
	msgs->say("calculating density done",Low);
	
	
	
#ifdef DEBUG_TEST_DENSITY_CALCULATIONS
	// mass test
	{
		double M=0;
		dx=getDx();
		dy=getDy();
		dz=getDz();
		D.print_domain_range(r);
		double dV;
		if (dz==0) dV=dx*dy; else dV=dx*dy*dz;
		for (long i = 0; i < r.subx; i++) {
			x=getX(i);
			for (long j = 0; j < r.suby; j++) {
				y=getY(j);
				for (long k = 0; k < r.subz; k++) {
					z=getZ(k);
					M+=f(i,j,k)[part]*dV;
				}
			}
		}
		printf("test: M=%lE (integrated mass in function domain, not in tree), dV=%lE\n",M,dV);
		cpedsPoint3D p(_param.x0+_param.sizeX/2,_param.y0+_param.sizeY/2,_param.z0+_param.sizeZ/2);
		p.print_point();
		printf("radius: %lf\n",cpeds_get_max(_param.sizeZ,_param.sizeX)/2);
		subDomain_region_t rr;
		rr.xmin=_param.x0; 	rr.xmax=_param.x0+_param.sizeX;
		rr.ymin=_param.y0; 	rr.ymax=_param.y0+_param.sizeY;
		rr.zmin=_param.z0; 	rr.zmax=_param.z0+_param.sizeZ;
		long NpartInSlice=D.getPointsIdx(rr).size();
		
		printf("mass in function region in the tree: M=%lE\n",	(double(NpartInSlice)));
		printf("mass in the tree: M=%lE\n",	(double(positions.size())));
		//		printf("mass in function region in the tree: M=%lE\n",	(double(NpartInSlice)));
		//		printf("mass loss: Mlost=%lE \%\n",	(double(positions.size())-M)/double(positions.size()));
		//		exit(0);
	}
#endif
	return *this;
}

/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::mkDensityFieldCube_old(subDomain_region_t r, subDomain_region_t treeScheme, const mscsVector<cpedsPoint3D>& positions, const mscsVector<double>* weights, int part, long maxDepth, bool spherical, bool calcMass, bool unitVolume) {
	
	double dx,dy,dz;
	
	// define the function space
	dx=(r.xmax-r.xmin)/r.subx;
	dy=(r.ymax-r.ymin)/r.suby;
	dz=(r.zmax-r.zmin)/r.subz;
	//	setSize(r.subx,r.suby,r.subz,dx,dy,dz,r.xmin,r.ymin,r.zmin);
	setSizeRange(r.subx,r.suby,r.subz, r.xmin, r.ymin, r.zmin, r.xmax, r.ymax, r.zmax);
	allocFunctionSpace();
	
	// prepare tree
	msgs->say("building tree",Medium);
	subDomain D(&positions,&treeScheme,1L); 
	D.setMaxDepth(maxDepth);
	D.tree();
	
	//	D.saveDomains("domains");
	msgs->say("done",Low);
	
	long N=getN();
	for (long i = 0; i < N; i++) {		setf(i,0,0);	}
	
	vector< vector<double> > v=D.getBranchParticlesCount(true);
	//	printf("number of domains: %li\n",v.size());
	long leafsCount,totalParticleCount=0;
	long particlesOnTree=0;
	double x,y,z,dv,w, dvTest=0;
	double massCount, mass;
	double drho;
	//	printInfo();
	
	dv=_param.dx*_param.dy;
	if (_param.Nz>=1) dv*=_param.dz;
	
	for (long i = 0; i < _param.Nx; i++) {
		x=getX(i);
		msgs->say("%li / %li\n",i,_param.Nx,Medium);
		for (long j = 0; j < _param.Ny; j++) {
			y=getY(j);
			//			printf("Ny=%li\n",_param.Ny);
			//			printInfo();
			//			msgs->say("j: %li / %li\n",j,_param.Ny,Medium);
			for (long k = 0; k < _param.Nz; k++) { // for every grid cell
				//				printf("Nz=%li\n",_param.Nz);
				z=getZ(k);
				leafsCount=0;
				massCount=0;
				mass=0;
				for (long l = 0; l < v.size(); l++) { // go through every subdomain
					if (v[l][7]==1) { // check if it is a leaf subdomain (it is redundant now as we now use leaf domains only)
						leafsCount++;
						w=0;
						
						// check if this subdomain has anything in common with current grid cell
						if (v[l][0] <= x+_param.dxo2 and v[l][1] >= x-_param.dxo2 and 
								v[l][2] <= y+_param.dyo2 and v[l][3] >= y-_param.dyo2 and
								v[l][4] <= z+_param.dzo2 and v[l][5] >= z-_param.dzo2 ) {
							
							// if so then for every particle in this subdomain, 
							particlesOnTree+=v[l].size()-8;
							for (long ii = 8; ii < v[l].size(); ii++) {	
								//recheck if it is inside of the current grid cell 
								// since the cell includes when it is inside  |- - -|
								//                                            |
								//                                            |     |     
								//                                            |
								//                                            |-----|
								//
								
								if (positions[v[l][long(ii)]].x() < x+_param.dxo2 and positions[v[l][long(ii)]].x() >= x-_param.dxo2 and 
										positions[v[l][long(ii)]].y() < y+_param.dyo2 and positions[v[l][long(ii)]].y() >= y-_param.dyo2 and
										positions[v[l][long(ii)]].z() < z+_param.dzo2 and positions[v[l][long(ii)]].z() >= z-_param.dzo2 ) {
									
									// and count the density contribution
									//									dv=(v[l][1]-v[l][0])*(v[l][3]-v[l][2]);
									//									if (_param.Nz>1) dv*=(v[l][5]-v[l][4]);
									dv=_param.dx;
									if (spherical) dv*=cos(PIsnd-(y+_param.dyo2))-cos(PIsnd-(y-_param.dyo2)); else dv*=_param.dy;
									if (_param.Nz>1) dv*=_param.dz;
									if (weights==NULL) w=1;
									else {
										w=weights->at(v[l][long(ii)]);
									}
									
									massCount+=1;
									mass+=w;
								}
								
								// special check for the last grid cells along given direction
								// but if it is the last cell along x, y or z then it must be like this
								//                                            |----|
								//                                            |    |
								//                                            |----|
								// so
								if (i==_param.Nx-1) {
									if (positions[v[l][long(ii)]].x() == _param.xMax  and 
											positions[v[l][long(ii)]].y() < y+_param.dyo2 and positions[v[l][long(ii)]].y() >= y-_param.dyo2 and
											positions[v[l][long(ii)]].z() < z+_param.dzo2 and positions[v[l][long(ii)]].z() >= z-_param.dzo2 ) {
										// and count the density contribution
										dv=_param.dx;
										if (spherical) dv*=cos(PIsnd-(y+_param.dyo2))-cos(PIsnd-(y-_param.dyo2)); else dv*=_param.dy;
										if (_param.Nz>1) dv*=_param.dz;
										if (weights==NULL) w=1;
										else {
											w=weights->at(v[l][long(ii)]);
										}
										massCount+=1;
										mass+=w;
									}
								}
								if (j==_param.Ny-1) {
									if (positions[v[l][long(ii)]].x() < x+_param.dxo2 and positions[v[l][long(ii)]].x() >= x-_param.dxo2 and 
											positions[v[l][long(ii)]].y() == _param.yMax and 
											positions[v[l][long(ii)]].z() < z+_param.dzo2 and positions[v[l][long(ii)]].z() >= z-_param.dzo2 ) {
										// and count the density contribution
										dv=_param.dx;
										if (spherical) dv*=cos(PIsnd-(y+_param.dyo2))-cos(PIsnd-(y-_param.dyo2)); else dv*=_param.dy;
										if (_param.Nz>1) dv*=_param.dz;
										if (weights==NULL) w=1;
										else {
											w=weights->at(v[l][long(ii)]);
										}
										massCount+=1;
										mass+=w;
									}
								}
								if (k==_param.Nz-1) {
									if (positions[v[l][long(ii)]].x() < x+_param.dxo2 and positions[v[l][long(ii)]].x() >= x-_param.dxo2 and 
											positions[v[l][long(ii)]].y() < y+_param.dyo2 and positions[v[l][long(ii)]].y() >= y-_param.dyo2 and
											positions[v[l][long(ii)]].z() == _param.zMax ) {
										// and count the density contribution
										dv=_param.dx;
										if (spherical) dv*=cos(PIsnd-(y+_param.dyo2))-cos(PIsnd-(y-_param.dyo2)); else dv*=_param.dy;
										if (_param.Nz>1) dv*=_param.dz;
										if (weights==NULL) w=1;
										else {
											w=weights->at(v[l][long(ii)]);
										}
										massCount+=1;
										mass+=w;
									}
								}
								// special check for the last grid cell along given directions at the same time
								// but if it is the last cell along x, y or z then it must be like this
								//                                            |----*
								//                                            |    |
								//                                            |----|
								// so
								if (i==_param.Nx-1 and j==_param.Ny-1) {
									if (positions[v[l][long(ii)]].x() == _param.xMax  and 
											positions[v[l][long(ii)]].y() == _param.yMax  and
											positions[v[l][long(ii)]].z() < z+_param.dzo2 and positions[v[l][long(ii)]].z() >= z-_param.dzo2 ) {
										// and count the density contribution
										dv=_param.dx;
										if (spherical) dv*=cos(PIsnd-(y+_param.dyo2))-cos(PIsnd-(y-_param.dyo2)); else dv*=_param.dy;
										if (_param.Nz>1) dv*=_param.dz;
										if (weights==NULL) w=1;
										else {
											w=weights->at(v[l][long(ii)]);
										}
										massCount+=1;
										mass+=w;
									}
								}
								
								
							}
							
						}
						
					}
					//					printf("%li %li %li %lE\n",i,j,k,v[l][6]);
				}
				totalParticleCount+=long(massCount);
				//				printf("total  number of particles accounted for till now: %li\n",totalParticleCount);
				if (calcMass) mass/=double(massCount);
				if (unitVolume) setf(i,j,k,mass,0);
				else
					setf(i,j,k,mass/dv,0);
#ifdef DEBUG
				dvTest+=dv;
#endif
				//				setf(i,j,k,mass,0);
				//				printf("dv: %lE\n",dv);
			}
		}
	}
	
#ifdef DEBUG
	
	
	particlesOnTree=0;
	for (long i = 0; i < v.size(); i++) {
		if (v[i][7]==1)
			particlesOnTree+=v[i].size()-8;
	}
	
	
	printf("total  number of particles checked out from tree: %li\n",particlesOnTree);
	printf("total  number of particles accounted for: %li\n",totalParticleCount);
	printf("total  number of particles : %li\n",positions.size());
	printf("leafs count: %li\n",v.size());
	
	
	printf("dvTest: %lE\n",dvTest);
	printf("\nnumber of leafs: %li\n",leafsCount);
#endif
	
	return *this;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::mkDensityFieldCube(subDomain_region_t r, subDomain_region_t treeScheme, const mscsVector<cpedsPoint3D>& positions, const mscsVector<double>* weights, long maxDepth, bool spherical, bool calcMass, bool unitVolume) {
	//	if (positions.size()==0) return *this;
	double dx,dy,dz;
	//
	// define the function space
	//
	dx=(r.xmax-r.xmin)/r.subx;
	dy=(r.ymax-r.ymin)/r.suby;
	dz=(r.zmax-r.zmin)/r.subz;
	//	setSize(r.subx,r.suby,r.subz,dx,dy,dz,r.xmin,r.ymin,r.zmin);
	setSizeRange(r.subx,r.suby,r.subz, r.xmin, r.ymin, r.zmin, r.xmax, r.ymax, r.zmax);
	allocFunctionSpace();
	setf(0,0);
	subDomain_region_t cell_reg;
	//
	// prepare tree
	//
	msgs->say("building tree",Medium);
	subDomain D(&positions,&treeScheme,2L); 
	D.setMaxDepth(maxDepth);
	D.tree();
	
	//	D.saveDomains("domains");
	msgs->say("done",Low);
	vector<long>* v; // pointer to a vector to containing particles in subdomain that contains the requested function grid cell
	
	long leafsCount,totalParticleCount=0;
	long particlesOnTree=0;
	double x,y,z,dv,w, dvTest=0;
	double massCount, mass;
	double drho;
	
	dv=_param.dx*_param.dy;
	if (dz!=0) dv*=_param.dz;
	
	for (long i = 0; i < _param.Nx; i++) {
		x=getX(i);
		msgs->say("i: %li / %li\n",i,_param.Nx,Medium);
		for (long j = 0; j < _param.Ny; j++) {
			y=getY(j);
			//			msgs->say("j: %li / %li\n",j,_param.Ny,Medium);
			for (long k = 0; k < _param.Nz; k++) { // for every grid cell
				//				msgs->say("k: %li / %li\n",k,_param.Nz,Medium);
				z=getZ(k);
				leafsCount=0;
				massCount=0;
				mass=0;
				cell_reg=getCellSDregion(i,j,k);
				//				msgs->say("obtaining getContainingDomainPointIdx",Medium);
				v=D.getContainingDomainPointIdx(cell_reg,NULL,true); // obtain the vector containing particles within spatial region enclosing i,j,k function grid cell
#ifdef DEBUG_SUBDOMAIN
				printf("containing SD has %li particles\n",long(v->size()));
#endif
				if (v!=NULL) { // 
					w=0;
					
					// for every particle in this subdomain, 
					particlesOnTree+=v->size();
					for (long ii = 0; ii < v->size(); ii++) {
						//recheck if it is inside of the current grid cell 
						// since the cell includes when it is inside  |- - -|
						//                                            |
						//                                            |     |     
						//                                            |
						//                                            |-----|
						//
						
						if (positions[v->at(ii)].x() < x+_param.dxo2 and positions[v->at(ii)].x() >= x-_param.dxo2 and 
								positions[v->at(ii)].y() < y+_param.dyo2 and positions[v->at(ii)].y() >= y-_param.dyo2 and
								positions[v->at(ii)].z() < z+_param.dzo2 and positions[v->at(ii)].z() >= z-_param.dzo2 ) {
							
							// and count the density contribution
							dv=_param.dx;
							if (spherical) dv*=cos(PIsnd-(y+_param.dyo2))-cos(PIsnd-(y-_param.dyo2)); else dv*=_param.dy;
							if (dz!=0) dv*=_param.dz;
							if (weights==NULL) w=1;
							else {
								w=weights->at(v->at(ii));
							}
							
							massCount+=1;
							mass+=w;
						}
						
						// special check for the last grid cells along given direction
						// but if it is the last cell along x, y or z then it must be like this
						//                                            |----|
						//                                            |    |
						//                                            |----|
						// so
						if (i==_param.Nx-1) {
							if (positions[v->at(ii)].x() == _param.xMax  and 
									positions[v->at(ii)].y() < y+_param.dyo2 and positions[v->at(ii)].y() >= y-_param.dyo2 and
									positions[v->at(ii)].z() < z+_param.dzo2 and positions[v->at(ii)].z() >= z-_param.dzo2 ) {
								// and count the density contribution
								dv=_param.dx;
								if (spherical) dv*=cos(PIsnd-(y+_param.dyo2))-cos(PIsnd-(y-_param.dyo2)); else dv*=_param.dy;
								if (dz!=0) dv*=_param.dz;
								if (weights==NULL) w=1;
								else {
									w=weights->at(v->at(ii));
								}
								massCount+=1;
								mass+=w;
							}
						}
						if (j==_param.Ny-1) {
							if (positions[v->at(ii)].x() < x+_param.dxo2 and positions[v->at(ii)].x() >= x-_param.dxo2 and 
									positions[v->at(ii)].y() ==  _param.yMax and 
									positions[v->at(ii)].z() < z+_param.dzo2 and positions[v->at(ii)].z() >= z-_param.dzo2 ) {
								// and count the density contribution
								dv=_param.dx;
								if (spherical) dv*=cos(PIsnd-(y+_param.dyo2))-cos(PIsnd-(y-_param.dyo2)); else dv*=_param.dy;
								if (dz!=0) dv*=_param.dz;
								if (weights==NULL) w=1;
								else {
									w=weights->at(v->at(ii));
								}
								massCount+=1;
								mass+=w;
							}
						}
						if (k==_param.Nz-1) {
							if (positions[v->at(ii)].x() < x+_param.dxo2 and positions[v->at(ii)].x() >= x-_param.dxo2 and 
									positions[v->at(ii)].y() < y+_param.dyo2 and positions[v->at(ii)].y() >= y-_param.dyo2 and
									positions[v->at(ii)].z() == _param.zMax ) {
								// and count the density contribution
								dv=_param.dx;
								if (spherical) dv*=cos(PIsnd-(y+_param.dyo2))-cos(PIsnd-(y-_param.dyo2)); else dv*=_param.dy;
								if (dz!=0) dv*=_param.dz;
								if (weights==NULL) w=1;
								else {
									w=weights->at(v->at(ii));
								}
								massCount+=1;
								mass+=w;
							}
						}
						// special check for the last grid cell along given directions at the same time
						//                                            |----*
						//                                            |    |
						//                                            |----|
						// so
						if (i==_param.Nx-1 and j==_param.Ny-1) {
							if (positions[v->at(ii)].x() == _param.xMax  and 
									positions[v->at(ii)].y() == _param.yMax  and
									positions[v->at(ii)].z() < z+_param.dzo2 and positions[v->at(ii)].z() >= z-_param.dzo2 ) {
								// and count the density contribution
								dv=_param.dx;
								if (spherical) dv*=cos(PIsnd-(y+_param.dyo2))-cos(PIsnd-(y-_param.dyo2)); else dv*=_param.dy;
								if (dz!=0) dv*=_param.dz;
								if (weights==NULL) w=1;
								else {
									w=weights->at(v->at(ii));
								}
								massCount+=1;
								mass+=w;
							}
						}
						
						
					} // loop over particles in the containing subdomain
					
				}
				totalParticleCount+=long(massCount);
				if (calcMass) 
					if (massCount==0) mass=0; else mass/=double(massCount); 
				if (unitVolume) setf(i,j,k,mass,0);
				else
					setf(i,j,k,mass/dv,0);
#ifdef DEBUG
				dvTest+=dv;
#endif
				//				setf(i,j,k,mass,0);
				//				printf("dv: %lE\n",dv);
			} // k
			//			exit(0);
		} // j
	} // i
	
#ifdef DEBUG
	
	
	particlesOnTree=0;
	for (long i = 0; i < v.size(); i++) {
		if (v[i][7]==1)
			particlesOnTree+=v[i].size()-8;
	}
	
	
	printf("total  number of particles checked out from tree: %li\n",particlesOnTree);
	printf("total  number of particles accounted for: %li\n",totalParticleCount);
	printf("total  number of particles : %li\n",positions.size());
	printf("leafs count: %li\n",v.size());
	
	
	printf("dvTest: %lE\n",dvTest);
	printf("\nnumber of leafs: %li\n",leafsCount);
#endif
	
	return *this;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::mkDensityFieldCube(int Nx, int Ny, int Nz, const cpedsPointSet3D& positions, long maxDepth, bool is2D, bool spherical, bool calcMass, bool unitVolume) {
	//	setSize(Nx,Ny,Nz);
	//	allocFunctionSpace();
	subDomain_region_t treeSubDomain;
	subDomain_region_t reg;
	//
	// set subdomain to match the positions maximal span in all directions
	//
	treeSubDomain.xmin=positions[0].x();	
	treeSubDomain.ymin=positions[0].y();
	treeSubDomain.zmin=positions[0].z();
	treeSubDomain.xmax=positions[0].x();
	treeSubDomain.ymax=positions[0].y();
	treeSubDomain.zmax=positions[0].z();
	for (long i = 0; i < positions.size(); i++) {
		if (positions[i].x()<treeSubDomain.xmin) treeSubDomain.xmin=positions[i].x();
		if (positions[i].y()<treeSubDomain.ymin) treeSubDomain.ymin=positions[i].y();
		if (positions[i].z()<treeSubDomain.zmin) treeSubDomain.zmin=positions[i].z();
		
		if (positions[i].x()>treeSubDomain.xmax) treeSubDomain.xmax=positions[i].x();
		if (positions[i].y()>treeSubDomain.ymax) treeSubDomain.ymax=positions[i].y();
		if (positions[i].z()>treeSubDomain.zmax) treeSubDomain.zmax=positions[i].z();
	}
	// use oct-tree
	treeSubDomain.subx=2;
	treeSubDomain.suby=2;
	if (is2D) {
		treeSubDomain.subz=1; 
		treeSubDomain.zmin-=1;
		treeSubDomain.zmax+=1;
	}
	else treeSubDomain.subz=2;
	
#ifdef DEBUG
	printf("xmin: %lE xmax: %lE\n",treeSubDomain.xmin,treeSubDomain.xmax);
	printf("ymin: %lE ymax: %lE\n",treeSubDomain.ymin,treeSubDomain.ymax);
	printf("zmin: %lE zmax: %lE\n",treeSubDomain.zmin,treeSubDomain.zmax);
#endif
	
	//
	// use the same area in function space but with possibly different resolution
	//
	reg=treeSubDomain;
	reg.subx=Nx; reg.suby=Ny; reg.subz=Nz;	
	
	return mkDensityFieldCube(reg,treeSubDomain,positions,positions.getValuesPtr(),maxDepth,spherical,calcMass,unitVolume);
}

/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::mkDensityFieldCube(int Nx, int Ny, int Nz, const cpedsPointSet3D& positions, bool is2D, bool spherical, bool calcMass, bool unitVolume) {
	subDomain_region_t treeSubDomain;
	subDomain_region_t reg;
	
	treeSubDomain.xmin=positions[0].x();	
	treeSubDomain.ymin=positions[0].y();
	treeSubDomain.zmin=positions[0].z();
	treeSubDomain.xmax=positions[0].x();
	treeSubDomain.ymax=positions[0].y();
	treeSubDomain.zmax=positions[0].z();
	for (long i = 0; i < positions.size(); i++) {
		if (positions[i].x()<treeSubDomain.xmin) treeSubDomain.xmin=positions[i].x();
		if (positions[i].y()<treeSubDomain.ymin) treeSubDomain.ymin=positions[i].y();
		if (positions[i].z()<treeSubDomain.zmin) treeSubDomain.zmin=positions[i].z();
		
		if (positions[i].x()>treeSubDomain.xmax) treeSubDomain.xmax=positions[i].x();
		if (positions[i].y()>treeSubDomain.ymax) treeSubDomain.ymax=positions[i].y();
		if (positions[i].z()>treeSubDomain.zmax) treeSubDomain.zmax=positions[i].z();
	}
	treeSubDomain.subx=Nx;
	treeSubDomain.suby=Ny;
	if (is2D) {
		treeSubDomain.subz=1; 
		treeSubDomain.zmin-=1;
		treeSubDomain.zmax+=1;
	}
	else treeSubDomain.subz=Nz;
	
	reg=treeSubDomain;
	reg.subx=Nx; reg.suby=Ny; reg.subz=Nz;	
	
	return mkDensityFieldCube(reg,treeSubDomain,positions,positions.getValuesPtr(),1,spherical,calcMass, unitVolume);
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::mkInterpolatedFieldScatter(subDomain_region_t r, subDomain_region_t treeScheme, const mscsVector<cpedsPoint3D>& positions, const mscsVector<double>& vals, string smKernel, long NeighborsMin, long NeighborsMax, mscsVector<double>* providedHSML, double MassMin, double MassMax) {
	pointsDensity dens(positions);
	dens.setVerbosityLevel(getVerbosityLevel());
#ifdef DEBUG_DENSITY
	dens.printRanges("positions ranges");
#endif
	bool is2Dcase=false;
	if (r.zmin==r.zmax) is2Dcase=true;
	if (is2Dcase)
		msgs->warning("2D case detected.", High);
	
	//
	// calculate points number density first
	//
	if (providedHSML!=NULL) {
		if (providedHSML->size()==positions.size())	dens.sml()=(*providedHSML);
		else {
			if (providedHSML->size()!=0) { 
				std::cerr << "Wrong providedHSML->size(). "
						"Should be either 0 or " << positions.size()<<std::endl;
				exit(1);
				
			}
//			else {
				// this means we don't know them 
				// but we want to know them when the calculation is done
//			}
		}
	}
	if (positions.size()<NeighborsMin) {
		NeighborsMin=positions.size();
		msgs->warning("The provided number of particles is smaller than the requested minimal number of neighbors. Decreasing NeighborsMin to the actual number of particles.", Top);
	}
	if (positions.size()<NeighborsMax) {
		NeighborsMax=positions.size();
		msgs->warning("The provided number of particles is smaller than the requested maximal number of neighbors. Decreasing NeighborsMax to the actual number of particles.", Top);
	}
	if (NeighborsMax<NeighborsMin) {
		NeighborsMax=NeighborsMin;
		msgs->warning("The maximal number of neighbors is smaller than the minimal number of neighbors. Will make them equal.", High);
	}
	
	dens.calculateDensity(NeighborsMin,NeighborsMax,is2Dcase,smKernel,&treeScheme,MassMin,MassMax);
	
	// we now have hsml calculated
	if (providedHSML!=0 and providedHSML->size()==0) { 
		*providedHSML=dens.sml();
	}
	
#ifdef DEBUG_DENSITY
	double *v=dens.density().toArray();
	cpeds_save_matrix(v,dens.size(),1,"mkInterpolatedFieldScatter_ptsDensity",true,false);
	v=dens.sml().toArray();
	cpeds_save_matrix(v,dens.size(),1,"mkInterpolatedFieldScatter_hsml",true,false);
#endif
	//
	// define the function space
	//
	double dx,dy,dz;
	dx=(r.xmax-r.xmin)/r.subx;
	dy=(r.ymax-r.ymin)/r.suby;
	dz=(r.zmax-r.zmin)/r.subz;
	setSize(r.subx,r.suby,r.subz,dx,dy,dz,r.xmin,r.ymin,r.zmin);
	allocFunctionSpace();
#ifdef DEBUG_DENSITY
	printInfo();
#endif
	
	//
	// define the smoothing kernel - this kernel is good for any smoothing length as it is not normalized (it is not in inverse volume units)
	// distance is in R=r/hsml
	double (*kernel)(double );
	double norm; //select kernel normalization depending on dimensionality
	if (smKernel == "gadget2") { 
		kernel=&mscsWindowFunction::kernelGadget; 	
		if (dz==0) norm=double(40.0)/(7.0*PI); //2d case
		else norm=8.0/PI; //3d case
	}
	else {
		if (smKernel == "gadget2b") { 
			kernel=&mscsWindowFunction::kernelGadget2b; 	
			if (dz==0) norm=double(40.0)/(7.0*PI); //2d case
			else norm=8.0/PI; //3d case
		}
		else {
			msgs->criticalError("mscsFunction3dregc::mkDensityField >> don't know this smoothing kernel function: "+smKernel,High);
		}
	}
	
	//
	// get the tree
	//
	subDomain *D=dens.getTree();
	
	//
	// set variables
	//
	double hsmlMax=dens.getMaxHSML();
	subDomain_region_t particle_reg;
	double hsml, vInt,dist;
	mscsVector<long> neighbors;
	mscsVector<double> sml;
	long acc=13;
	if (is2Dcase) { // 2D case
		hsml=sqrt(dx*dx+dy*dy);
	}
	else { // 3D case
		hsml=sqrt(dx*dx+dy*dy+dz*dz);
	}
	double dr2d=sqrt(dx*dx+dy*dy);
	double dr3d=sqrt(dx*dx+dy*dy+dz*dz);
	double dr;
	if (is2Dcase) dr=dr2d; else dr=dr3d;
	
	//
	// do the interpolation
	//
	msgs->say("calculating interpolated values at grid cells",Medium);
	msgs->say("hsmlMax: %lE",hsmlMax,Medium);
	double x,y,z,Neff,kernVal;
	//#ifdef HAVE_OPENMP
	//	omp_set_num_threads(omp_get_num_procs());
	//#endif
	//#pragma omp parallel num_threads(omp_get_num_procs())
#pragma omp parallel for schedule (guided) firstprivate(D,hsml,x,y,z,hsmlMax, neighbors,r,Neff,kernVal, vInt, dist,norm) shared(positions, sml) 
	//#pragma omp parallel for firstprivate(D,hsml,x,y,z,hsmlMax, neighbors,r,Neff,kernVal, vInt, dist,norm) shared(positions, sml) 
	for (long i = 0; i < r.subx; i++) {
#ifdef HAVE_OPENMP
		if (omp_get_thread_num()==0) 
			msgs->say("done: "+msgs->toStr(double(i)*omp_get_num_threads()/r.subx*100)+" %",Low);
#else
		msgs->say("done: %lf\%",double(i)/r.subx*100,Low);
#endif
		x=getX(i);
		for (long j = 0; j < r.suby; j++) {
			y=getY(j);
			for (long k = 0; k < r.subz; k++) {
				z=getZ(k);
				
				cpedsPoint3D p0=cpedsPoint3D(x,y,z);
#ifdef DEBUG_DENSITY
				//				p0=cpedsPoint3D(1.09855697694289025591E-15,-1.51987663191803165641E-02,0);
				for (long o=0;o<positions.size();o++) {
					p0=positions[o];
#endif
					
					//				cpedsPoint3D p(x,y,z); p.print_point("current point");
					// get particles within the distance of R defining region containing constant number of particles
					neighbors=D->getPointsIdx(p0,hsmlMax,NULL);
					vInt=0;
					Neff=0;
					if (is2Dcase) { //2d case
						for (long l = 0; l < neighbors.size(); l++) {
							dist=p0.dist(positions[neighbors[l]]);
							hsml=dens.sml(neighbors[l]);
							kernVal=kernel(dist/hsml);
							if (kernVal>0) Neff++; // debug
							vInt+=vals.at(neighbors[l])*kernVal/(hsml*hsml)  / dens.rho(neighbors[l]);
#ifdef DEBUG_DENSITY
							//						printf("%li) %li %li %li, %lf, %lf, %lf, vInt: %lE, kern: %lE, hsml: %lE\n",l,i,j,k,x,y,z,vInt,kernVal,hsml);
#endif
						}
					}
					else { // 3d case
						for (long l = 0; l < neighbors.size(); l++) {
							dist=cpedsPoint3D(x,y,z).dist(positions[neighbors[l]]);
							hsml=dens.sml(neighbors[l]);
							kernVal=kernel(dist/hsml);
							if (kernVal>0) Neff++; // debug
							vInt+=vals.at(neighbors[l])*kernVal/(hsml*hsml*hsml) / dens.rho(neighbors[l]);
						}
					}
					f(i,j,k)[0]=norm*vInt;
					f(i,j,k)[1]=Neff;
#ifdef DEBUG_DENSITY
					printf("--- %li %li %li, %lf, %lf, %lf, fInt: %lE, checkVal: %lE, relErr: %lf\n",i,j,k,x,y,z,fRe(i,j,k),vals.at(o),(vals.at(o)-fRe(i,j,k))/vals.at(o)*100);
				}
#endif
			}
		}
	}
	msgs->say("interpolation done",Low);
	
	return *this;
}
/***************************************************************************************/
//mscsFunction3dregc& mkInterpolatedFieldScatterMass(subDomain_region_t r, subDomain_region_t treeScheme, const mscsVector<cpedsPoint3D>& positions, const mscsVector<double>& vals, string smKernel="gadget2", double MassMin, double MassMax, mscsVector<double>* providedHSML) {
//	pointsDensity dens(positions);
//	dens.setVerbosityLevel(getVerbosityLevel());
//#ifdef DEBUG_DENSITY
//	dens.setVerbosityLevel(getVerbosityLevel());
//	dens.printRanges("positions ranges");
//#endif
//	bool is2Dcase=false;
//	if (r.zmin==r.zmax) is2Dcase=true;
//	if (is2Dcase)
//		msgs->warning("2D case detected.", High);
//	
//	//
//	// calculate points number density first
//	//
//	if (providedHSML!=NULL) {
//		dens.sml()=(*providedHSML);
//	}
//	if (positions.size()<NeighborsMin) {
//		NeighborsMin=positions.size();
//		msgs->warning("The provided number of particles is smaller than the requested minimal number of neighbors. Decreasing NeighborsMin to the actual number of particles.", Top);
//	}
//	if (positions.size()<NeighborsMax) {
//		NeighborsMax=positions.size();
//		msgs->warning("The provided number of particles is smaller than the requested maximal number of neighbors. Decreasing NeighborsMax to the actual number of particles.", Top);
//	}
//	if (NeighborsMax<NeighborsMin) {
//		NeighborsMax=NeighborsMin;
//		msgs->warning("The maximal number of neighbors is smaller than the minimal number of neighbors. Will make them equal.", High);
//	}
//
//	dens.calculateDensity(NeighborsMin,NeighborsMax,is2Dcase,smKernel,&treeScheme);
//#ifdef DEBUG_DENSITY
//	double *v=dens.density().toArray();
//	cpeds_save_matrix(v,dens.size(),1,"mkInterpolatedFieldScatter_ptsDensity",true,false);
//	v=dens.sml().toArray();
//	cpeds_save_matrix(v,dens.size(),1,"mkInterpolatedFieldScatter_hsml",true,false);
//#endif
//	//
//	// define the function space
//	//
//	double dx,dy,dz;
//	dx=(r.xmax-r.xmin)/r.subx;
//	dy=(r.ymax-r.ymin)/r.suby;
//	dz=(r.zmax-r.zmin)/r.subz;
//	setSize(r.subx,r.suby,r.subz,dx,dy,dz,r.xmin,r.ymin,r.zmin);
//	allocFunctionSpace();
//#ifdef DEBUG_DENSITY
//	printInfo();
//#endif
//	
//	//
//	// define the smoothing kernel - this kernel is good for any smoothing length as it is not normalized (it is not in inverse volume units)
//	// distance is in R=r/hsml
//	double (*kernel)(double );
//	double norm; //select kernel normalization depending on dimensionality
//	if (smKernel == "gadget2") { 
//		kernel=&mscsWindowFunction::kernelGadget; 	
//		if (dz==0) norm=double(40.0)/(7.0*PI); //2d case
//		else norm=8.0/PI; //3d case
//	}
//	else msgs->criticalError("mscsFunction3dregc::mkDensityField >> don't know this smoothing kernel function: "+smKernel,High);
//
//	//
//	// get the tree
//	//
//	subDomain *D=dens.getTree();
//
//	//
//	// set variables
//	//
//	double hsmlMax=dens.getMaxHSML();
//	subDomain_region_t particle_reg;
//	double hsml, vInt,dist;
//	mscsVector<long> neighbors;
//	mscsVector<double> sml;
//	long acc=13;
//	if (is2Dcase) { // 2D case
//		hsml=sqrt(dx*dx+dy*dy);
//	}
//	else { // 3D case
//		hsml=sqrt(dx*dx+dy*dy+dz*dz);
//	}
//	double dr2d=sqrt(dx*dx+dy*dy);
//	double dr3d=sqrt(dx*dx+dy*dy+dz*dz);
//	double dr;
//	if (is2Dcase) dr=dr2d; else dr=dr3d;
//
//	//
//	// do the interpolation
//	//
//	msgs->say("calculating interpolated values at grid cells",Medium);
//	msgs->say("hsmlMax: %lE",hsmlMax,Medium);
//	double x,y,z,Neff,kernVal;
////#ifdef HAVE_OPENMP
////	omp_set_num_threads(omp_get_num_procs());
////#endif
////#pragma omp parallel num_threads(omp_get_num_procs())
//#pragma omp parallel for schedule (guided) firstprivate(D,hsml,x,y,z,hsmlMax, neighbors,r,Neff,kernVal, vInt, dist,norm) shared(positions, sml) 
////#pragma omp parallel for firstprivate(D,hsml,x,y,z,hsmlMax, neighbors,r,Neff,kernVal, vInt, dist,norm) shared(positions, sml) 
//	for (long i = 0; i < r.subx; i++) {
//#ifdef HAVE_OPENMP
//		if (omp_get_thread_num()==0) 
//			msgs->say("done: "+msgs->toStr(double(i)*omp_get_num_threads()/r.subx*100)+" %",Low);
//#else
//		msgs->say("done: %lf\%",double(i)/r.subx*100,Low);
//#endif
//		x=getX(i);
//		for (long j = 0; j < r.suby; j++) {
//			y=getY(j);
//			for (long k = 0; k < r.subz; k++) {
//				z=getZ(k);
//
//				cpedsPoint3D p0=cpedsPoint3D(x,y,z);
//#ifdef DEBUG_DENSITY
////				p0=cpedsPoint3D(1.09855697694289025591E-15,-1.51987663191803165641E-02,0);
//				for (long o=0;o<positions.size();o++) {
//					p0=positions[o];
//#endif
//				
//				//				cpedsPoint3D p(x,y,z); p.print_point("current point");
//				// get particles within the distance of R defining region containing constant number of particles
//				neighbors=D->getPointsIdx(p0,hsmlMax,NULL);
//				vInt=0;
//				Neff=0;
//				if (is2Dcase) { //2d case
//					for (long l = 0; l < neighbors.size(); l++) {
//						dist=p0.dist(positions[neighbors[l]]);
//						hsml=dens.sml(neighbors[l]);
//						kernVal=kernel(dist/hsml);
//						if (kernVal>0) Neff++; // debug
//						vInt+=vals.at(neighbors[l])*kernVal/(hsml*hsml)  / dens.rho(neighbors[l]);
//#ifdef DEBUG_DENSITY
////						printf("%li) %li %li %li, %lf, %lf, %lf, vInt: %lE, kern: %lE, hsml: %lE\n",l,i,j,k,x,y,z,vInt,kernVal,hsml);
//#endif
//					}
//				}
//				else { // 3d case
//					for (long l = 0; l < neighbors.size(); l++) {
//						dist=cpedsPoint3D(x,y,z).dist(positions[neighbors[l]]);
//						hsml=dens.sml(neighbors[l]);
//						kernVal=kernel(dist/hsml);
//						if (kernVal>0) Neff++; // debug
//						vInt+=vals.at(neighbors[l])*kernVal/(hsml*hsml*hsml) / dens.rho(neighbors[l]);
//					}
//				}
//				f(i,j,k)[0]=norm*vInt;
//				f(i,j,k)[1]=Neff;
//#ifdef DEBUG_DENSITY
//				printf("--- %li %li %li, %lf, %lf, %lf, fInt: %lE, checkVal: %lE, relErr: %lf\n",i,j,k,x,y,z,fRe(i,j,k),vals.at(o),(vals.at(o)-fRe(i,j,k))/vals.at(o)*100);
//				}
//#endif
//			}
//		}
//	}
//	msgs->say("interpolation done",Low);
//
//	return *this;
//	
//}
/***************************************************************************************/
//mscsFunction3dregc& mscsFunction3dregc::mkInterpolatedFieldScatter(subDomain_region_t r, subDomain_region_t treeScheme, const mscsVector<cpedsPoint3D>& positions, const mscsVector<double>& vals, mscsVector<double>& providedHSML, mscsVector<double>& providedDens, string smKernel) {
mscsFunction3dregc& mscsFunction3dregc::mkInterpolatedFieldScatter(subDomain_region_t r, subDomain* D, const mscsVector<cpedsPoint3D>& positions, const mscsVector<double>& vals, mscsVector<double>& providedHSML, mscsVector<double>& providedDens, double maxHSML, string smKernel, subDomain* D3) {
	bool is2Dcase=false;
	if (r.zmin==r.zmax) is2Dcase=true;
	
	//	long D3activationPartNumThres=positions.size();
	long D3activationPartNumThres=positions.size()/512;
	
	//
	// define the function space
	//
	double dx,dy,dz;
	dx=(r.xmax-r.xmin)/r.subx;
	dy=(r.ymax-r.ymin)/r.suby;
	dz=(r.zmax-r.zmin)/r.subz;
	setSize(r.subx,r.suby,r.subz,dx,dy,dz,r.xmin,r.ymin,r.zmin);
	allocFunctionSpace();
	
	
	cpedsPoint3D p0ctr((r.xmax+r.xmin)/2,(r.ymax+r.ymin)/2,(r.zmax+r.zmin)/2);
	double regRadii=sqrt((r.xmax-r.xmin)*(r.xmax-r.xmin)+(r.ymax-r.ymin)*(r.ymax-r.ymin)+(r.zmax-r.zmin)*(r.zmax-r.zmin))/2;
	//
	// define the smoothing kernel - this kernel is good for any smoothing length as it is not normalized (it is not in inverse volume units)
	// distance is in R=r/hsml
	double (*kernel)(double );
	double norm; //select kernel normalization depending on dimensionality
	if (smKernel == "gadget2") { 
		kernel=&mscsWindowFunction::kernelGadget; 	
		if (dz==0) norm=double(40.0)/(7.0*PI); //2d case
		else norm=8.0/PI; //3d case
	}
	else {
		if (smKernel == "gadget2b") { 
			kernel=&mscsWindowFunction::kernelGadget2b; 	
			if (dz==0) norm=double(40.0)/(7.0*PI); //2d case
			else norm=8.0/PI; //3d case
		}
		else {
			msgs->criticalError("mscsFunction3dregc::mkDensityField >> don't know this smoothing kernel function: "+smKernel,High);
		}
	}
	
	
	//
	// set variables
	//
	double hsmlMax=cpeds_find_max_value(providedHSML.data(),providedHSML.size(),0);
	double hsmlMaxProvided=hsmlMax;
	//	if (hsmlMax>maxHSML and maxHSML!=-1) hsmlMax=maxHSML;
	if (maxHSML!=-1) hsmlMax=maxHSML;
	subDomain_region_t particle_reg;
	double hsml, vInt,dist;
	mscsVector<long> neighbors;
	mscsVector<long> neighbors3;
	mscsVector<double> sml;
	long acc=13;
	if (is2Dcase) { // 2D case
		hsml=sqrt(dx*dx+dy*dy);
	}
	else { // 3D case
		hsml=sqrt(dx*dx+dy*dy+dz*dz);
	}
	double dr2d=sqrt(dx*dx+dy*dy);
	double dr3d=sqrt(dx*dx+dy*dy+dz*dz);
	double dr;
	if (is2Dcase) dr=dr2d; else dr=dr3d;
	
	providedHSML.printVector();
	providedDens.printVector();
	printf("region: xmin: %lf, xmax: %lf\n",r.xmin,r.xmax);
	printf("region: ymin: %lf, ymax: %lf\n",r.ymin,r.ymax);
	printf("region: zmin: %lf, zmax: %lf\n",r.zmin,r.zmax);
	D->print_domain_range();
	//
	// do the interpolation
	//
	msgs->say("calculating interpolated values at grid cells",Medium);
	msgs->say("hsmlMax: %lE",hsmlMax,Medium);
	double x=0,y=0,z=0,Neff=0,kernVal=0;
	
	cpedsPoint3D p0prev=p0ctr;
	
	if (regRadii>hsmlMax/2) {
		msgs->say("will look for neighbors at every grid cell (regRadii: %lE, hsmlMax: %lE",regRadii,hsmlMax,Low);
		//		hsmlMax=hsmlMaxProvided;
		msgs->say("falling back to maximal provided hsml (hsmlMaxProvided: %lE)",hsmlMaxProvided,Low);
	}
	else {
		neighbors=D->getPointsIdx(p0ctr,hsmlMax,NULL);  
		msgs->say("neighbors wrt region center: %li",long(neighbors.size()),Low);		
#ifdef DEBUG_DENSITY
		FILE* tmpf=fopen("neigh.tmp","w");
		for (unsigned long i = 0; i < neighbors.size(); i++) {
			fprintf(tmpf,"%lE %lE %lE %lE\n",positions[neighbors[i]].x(),positions[neighbors[i]].y(),positions[neighbors[i]].z(),vals.at(neighbors[i]));
		}
		fclose(tmpf);
#endif
	}
	
	//#ifdef HAVE_OPENMP
	//	omp_set_num_threads(omp_get_num_procs());
	//#endif
	//#pragma omp parallel for 
#pragma omp parallel for schedule (guided) \
		firstprivate(hsmlMaxProvided,regRadii,dr,D,D3,hsml,x,y,z,hsmlMax, neighbors,neighbors3,r,Neff,kernVal, vInt, dist,norm,D3activationPartNumThres) \
private(p0prev) \
shared(positions, sml) 
	for (long i = 0; i < r.subx; i++) {
#ifdef HAVE_OPENMP
		if (omp_get_thread_num()==0) 
			msgs->say("done: "+msgs->toStr(double(i)*omp_get_num_threads()/r.subx*100)+" %",Low);
#else
		msgs->say("done: %lf\%",double(i)/r.subx*100,Low);
#endif
		x=getX(i);
		for (long j = 0; j < r.suby; j++) {
			//#ifdef HAVE_OPENMP
			//		if (omp_get_thread_num()==0) 
			//			msgs->say("j:"+msgs->toStr(double(j)),Low);
			//#endif
			y=getY(j);
			for (long k = 0; k < r.subz; k++) {
				z=getZ(k);
				
				cpedsPoint3D p0=cpedsPoint3D(x,y,z);
				
				// get particles within the distance of R defining region containing constant number of particles
				//				if (p0.dist(p0prev) > hsmlMax or neighbors.size()==0)
				if (regRadii>hsmlMax/2)
					neighbors=D->getPointsIdx(p0,hsmlMaxProvided,NULL);  
				/*
				 * Comment: getting neighbors for every grid cell is sub-optimal.
				 * If hsml is large enough (and we used hsmlMax) then it is possible that
				 * the same set of particles will suffice to calculate interpolated value for the
				 * neighboring (next) grid cell. The optimization here could include
				 * a) searching for particles on the tree only if the current point p0 is located outside of the
				 * sphere (p0_prev, hsmlMax) where p0_prev is the previous grid point
				 * b) walk over the grid in a nested healpix-type manner and not along rows,cols etc.
				 * This is possible only for the grid sizes of 2^n on side.
				 * c) if the region of interpolation is defined across one of the top-tree domains
				 * then at every pass the full set of particles is taken and the algorithm is highly inefective
				 * in selecting particles. This can be improved by using two-trees shifted wrt each other
				 * by eg hsmlMax. or by using a combination of oct-tree and 3^3 type tree -i.e. with 3 generic subdivisions.
				 * 
				 * author: blew
				 * date: Oct 10, 2013 9:14:39 AM
				 *
				 */
				if (D3!=NULL) { // this fixes the problem indicated in point c) - comment above, but the code is not tested and something works wrong, DEBUG needed
					//					printf("oct-tree neighbors: %li\n",neighbors.size());
					if (neighbors.size()>D3activationPartNumThres) {
						msgs->say("searching tree3 neighbors", Low);
						neighbors3=D3->getPointsIdx(p0,hsmlMax,NULL);
						if (neighbors3.size()<neighbors.size()) {
							msgs->say("using tree3 neighbors: new neighbors %li, old oct-tree neighbors %li",long(neighbors3.size()),long(neighbors.size()),Low);
							neighbors=neighbors3;
						}
					}
				}
				
				
				vInt=0;
				Neff=0;
				if (is2Dcase) { //2d case
					for (long l = 0; l < neighbors.size(); l++) {
						dist=p0.dist(positions[neighbors[l]]);
						hsml=providedHSML[neighbors[l]];
						kernVal=kernel(dist/hsml);
						if (kernVal>0) Neff++; // debug
						vInt+=vals.at(neighbors[l])*kernVal/(hsml*hsml)  / providedDens[neighbors[l]];
					}
				}
				else { // 3d case
					for (long l = 0; l < neighbors.size(); l++) {
						dist=cpedsPoint3D(x,y,z).dist(positions[neighbors[l]]);
						hsml=providedHSML[neighbors[l]];
						kernVal=kernel(dist/hsml);
						if (kernVal>0) Neff++; // debug
						vInt+=vals.at(neighbors[l])*kernVal/(hsml*hsml*hsml) / providedDens[neighbors[l]];
					}
				}
				f(i,j,k)[0]=norm*vInt;
				f(i,j,k)[1]=Neff;
			}
		}
	}
	msgs->say("interpolation done",Low);
	
	return *this;	
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::mkInterpolatedFieldTriangSibsonGrad2D(const mscsVector<cpedsPoint3D>& positions) {
//	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//	typedef CGAL::Delaunay_triangulation_2<K>               Delaunay_triangulation;
//	typedef CGAL::Interpolation_gradient_fitting_traits_2<K> Traits;
//	typedef K::FT                                            Coord_type;
//	typedef K::Point_2                                       Point;
//	typedef std::map<Point, Coord_type, K::Less_xy_2>        Point_value_map ;
//	typedef std::map<Point, K::Vector_2 , K::Less_xy_2 >     Point_vector_map;
//
//	Delaunay_triangulation T;
//	std::map<K::Point_2, Coord_type, K::Less_xy_2> function_values;
//	typedef CGAL::Data_access< std::map<K::Point_2, Coord_type, K::Less_xy_2 > > Value_access;
//
//	for (long i = 0; i < positions.size(); i++) {
//		K::Point_2 p(positions[i].x(),positions[i].y());
//		T.insert(p);
//		function_values.insert(std::make_pair(p,positions[i].z()));
//	}
//
	return *this;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::mkInterpolatedFieldTriangLinear2D(const mscsVector<cpedsPoint3D>& positions) {
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//	typedef CGAL::Exact_predicates_exact_constructions_kernel K;
	typedef CGAL::Delaunay_triangulation_2<K>             Delaunay_triangulation;
	typedef CGAL::Interpolation_traits_2<K>               Traits;
	typedef K::FT                                         Coord_type;
	typedef K::Point_2                                    Point;
	typedef std::vector< std::pair< Point, Coord_type  > > Point_coordinate_vector;
	
	Delaunay_triangulation T;
	std::map<Point, Coord_type, K::Less_xy_2> function_values;
	typedef CGAL::Data_access< std::map<Point, Coord_type, K::Less_xy_2 > > Value_access;
	
	//import data
	for (long i = 0; i < positions.size(); i++) {
		K::Point_2 p(positions[i].x(),positions[i].y());
		T.insert(p);
		function_values.insert(std::make_pair(p,positions[i].z()));
	}

	//coordiante computation
	Point_coordinate_vector coords;
	
//	for (long k = 0; k< Nz(); k++) {
	long k=0,i,j;
	for (i = 0; i < Nx(); i++) { 
#pragma omp parallel for private(j) 
		for (j = 0; j < Ny(); j++) {
			//coordinate computation
			K::Point_2 p(getX(i),getY(j));
			Point_coordinate_vector coords;

/*
			CGAL::Triple< std::back_insert_iterator<Point_coordinate_vector>, Coord_type, bool> result;
			result=CGAL::natural_neighbor_coordinates_2(T, p,std::back_inserter(coords));
			Coord_type norm =result.second;
*/

			Coord_type norm =
					CGAL::natural_neighbor_coordinates_2
					(T, p,std::back_inserter(coords)).second;
//			if (result.third) {
				Coord_type res = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(function_values));
				fRe(i,j,k)=res;
//			}
		}
	}
//	}


	return *this;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::mkInterpolatedFieldTriangLinear2D(const cpedsPointSet3D& positions) {
	mscsVector<cpedsPoint3D> pos=positions.exportAsVector();
	mkInterpolatedFieldTriangLinear2D(pos);
	return *this;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::mkSin3Drad(double from, double to, double dr, double T, double phi) {
	double x,y,z;
	long i,j,k;
	
	// check the sizes of the output array
	x=0; i=0;	while (x<=to) {		x+=dr;		i++;	}
	//	y=0; j=0;	while (y<=to) {		y+=dr;		j++;	}
	//	z=0; k=0;	while (z<=to) {		z+=dr;		k++;	}
	
	setSize(i,i,i,dr,dr,dr,from,from,from);
	realloc();
	
	// make function
	x=from; i=0;
	while (x<=to) {
		y=from; j=0;
		while (y<=to) {
			z=from; k=0;
			while (z<=to) {
				fRe(i,j,k)=sin(twoPI/T*sqrt(x*x + y*y + z*z));
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
mscsFunction3dregc& mscsFunction3dregc::mkSin3D(double from, double to, double dr, double T, double phi) {
	double x,y,z;
	long i,j,k;
	
	// check the sizes of the output array
	x=0; i=0;	while (x<=to) {		x+=dr;		i++;	}
	//	y=0; j=0;	while (y<=to) {		y+=dr;		j++;	}
	//	z=0; k=0;	while (z<=to) {		z+=dr;		k++;	}
	
	setSize(i,i,i,dr,dr,dr,from,from,from);
	realloc();
	
	// make function
	x=from; i=0;
	while (x<=to) {
		y=from; j=0;
		while (y<=to) {
			z=from; k=0;
			while (z<=to) {
				fRe(i,j,k)=sin(twoPI/T*(x + y + z));
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
mscsFunction3dregc& mscsFunction3dregc::mkGauss3D(double sx, double sy, double sz) {
	long nx=Nx();
	long ny=Ny();
	long nz=Nz();
	double sx2,sy2,sz2, sx2inv, sy2inv, sz2inv;
	// convert variances to units of the grid
	sx*=double(nx)/_param.sizeX;
	sy*=double(ny)/_param.sizeY;
	sz*=double(nz)/_param.sizeZ;
	
	sx2=2.0*sx*sx;
	sy2=2.0*sy*sy;
	sz2=2.0*sz*sz;

	sx2inv=1.0/sx2;
	sy2inv=1.0/sy2;
	sz2inv=1.0/sz2;
	
	if (_param.sizeZ==0)  sz2inv=0;

	
	double x,y,z;
	
	for (long i = 0; i < nx; i++) {
		x=double(i)-double(nx)/2;
		for (long j = 0; j < ny; j++) {
			y=double(j)-double(ny)/2;
			for (long k = 0; k < nz; k++) {
				z=double(k)-double(nz)/2;
				f(i,j,k)[0]=exp(- ( 
						x*x*sx2inv +  
						y*y*sy2inv +
						z*z*sz2inv
				)
				);
				f(i,j,k)[1]=0.0;
			}
		}
	}
	return *this;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::mkTopHat3D(double R, double acc) {
	long nx=Nx();
	long ny=Ny();
	long nz=Nz();
	double x0=lengthX()/2;
	double y0=lengthY()/2;
	double z0=lengthZ()/2;
	double x,y,z;
	
	for (long i = 0; i < nx; i++) {
		x=getX(i);
		for (long j = 0; j < ny; j++) {
			y=getY(j);
			for (long k = 0; k < nz; k++) {
				z=getZ(k);
				
				if (fabs(x-x0)-R<=acc and 
						fabs(y-y0)-R<=acc and 
						fabs(z-z0)-R<=acc) 
					f(i,j,k)[0]=1.0; 
				else f(i,j,k)[0]=0.0;
				
				f(i,j,k)[1]=0.0;
			}
		}
	}
	return *this;	
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::mkBall3D(double R, double acc) {
	long nx=Nx();
	long ny=Ny();
	long nz=Nz();
	double x0=getFunctionParameters().x0+lengthX()/2;
	double y0=getFunctionParameters().y0+lengthY()/2;
	double z0=getFunctionParameters().z0+lengthZ()/2;
	double x,y,z;
	
	for (long i = 0; i < nx; i++) {
		x=getX(i);
		for (long j = 0; j < ny; j++) {
			y=getY(j);
			for (long k = 0; k < nz; k++) {
				z=getZ(k);
				
				if (pow(x-x0,2)+pow(y-y0,2)+pow(z-z0,2)-pow(R,2)<=acc)
					f(i,j,k)[0]=1.0; 
				else f(i,j,k)[0]=0.0;
				
				f(i,j,k)[1]=0.0;
			}
		}
	}
	return *this;		
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::mkBall3Dat(double x0, double y0, double z0, double R, double acc) {
	long nx=Nx();
	long ny=Ny();
	long nz=Nz();
	double x,y,z;
	
	for (long i = 0; i < nx; i++) {
		x=getX(i);
		for (long j = 0; j < ny; j++) {
			y=getY(j);
			for (long k = 0; k < nz; k++) {
				z=getZ(k);
				
				if (pow(x-x0,2)+pow(y-y0,2)+pow(z-z0,2)-pow(R,2)<=acc)
					f(i,j,k)[0]=1.0; 
				else f(i,j,k)[0]=0.0;
				
				f(i,j,k)[1]=0.0;
			}
		}
	}
	return *this;		
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::mkGauss2D(double c1, double c2, double s1, double s2, double amplitude, long plane, long coord) {
	if (plane==0) {
		msgs->criticalError("mkGauss2D: not implemented yet",Top);
	}
	if (plane==1) {		
		msgs->criticalError("mkGauss2D: not implemented yet",Top);
	}
	if (plane==2) {
		double x,y;
		for (long i = 0; i < Nx(); i++) {
			x=getX(i);
			for (long j = 0; j < Ny(); j++) {
				y=getY(j);
				fRe(i,j,coord)=amplitude*exp(-(x-c1)*(x-c1)/(s1*s1*2.0)-(y-c2)*(y-c2)/(s2*s2*2.0));
			}
		}
	}
	return *this;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::mkGauss2D(double c1, double c2, double s11, double s22, double s12, double amplitude, long plane, long coord) {
	matrix<double> C(2,2);
	C(0,0)=s11;
	C(0,1)=s12;
	C(1,0)=s12;
	C(1,1)=s22;
	C=C.Inv();
	
	if (plane==0) {
		msgs->criticalError("mkGauss2D: not implemented yet",Top);
	}
	if (plane==1) {		
		msgs->criticalError("mkGauss2D: not implemented yet",Top);
	}
	if (plane==2) {
		double x,y;
		for (long i = 0; i < Nx(); i++) {
			x=getX(i);
			for (long j = 0; j < Ny(); j++) {
				y=getY(j);
				fRe(i,j,coord)=amplitude*exp(-0.5*( (x-c1)*(x-c1)*C(0,0) + (y-c2)*(y-c2)*C(1,1) + 2.0*(x-c1)*(y-c2)*C(0,1) ));
			}
		}
	}
	return *this;	
}
/***************************************************************************************/
mscsFunction3dregc mscsFunction3dregc::mkKernel(mscsFunction& Pk, string interpolationType) {

	allocFunctionSpace();
	double kx,ky,kz,modk;
	long nx=Nx();
	long ny=Ny();
	long nz=Nz();
	
	if (lengthZ()!=0) { // 3D case
		nz=Nz()/2+1;
	}
	else { // 2D case
		ny=Ny()/2+1;
		assert(nz==1);
	}

	
	double* kk = new double[size()];
	long l=0;
//	long l1=0,l2=0;
	double kmin,kmax;
	Pk.sortFunctionArgAscending();
	kmin=Pk.getMinArg();
	kmax=Pk.getMaxArg();
//	long ik0=0;
	long Nbelow=0, Nabove=0;
	//	double Kmax=double(nx-1)/double(_param.sizeX);
	
	msgs->say("Pre-calculating the k-modes and interpolating the power spectrum",High);
	for (long i = 0; i < nx; i++) {
		for (long j = 0; j < ny; j++) {
			for (long k = 0; k < nz; k++) {	

				kx=double(i)/double(_param.sizeX);
				ky=double(j)/double(_param.sizeY);
				if (_param.sizeZ==0) kz=0; else kz=double(k)/double(_param.sizeZ);
				
				modk=sqrt(kx*kx + ky*ky + kz*kz);
				if (modk<kmin) { msgs->warning("P(k) is not sufficiently tabulated for small k -- using the closest, minimal k value: "+msgs->toStr(kmin)+"  (value requested: "+msgs->toStr(modk)+")",High); modk=kmin; Nbelow++; }
				if (modk>kmax) { msgs->warning("P(k) is not sufficiently tabulated for large k -- using the closest, maximal k value"+msgs->toStr(kmax)+"  (value requested: "+msgs->toStr(modk)+")",High); modk=kmax; Nabove++; }
				kk[l]=modk;
				l++;
			}
		}
	}
	printf("l: %li\n",l);
	printf("Nbelow: %li\n",Nbelow);
	printf("Nabove: %li\n",Nabove);
	//
	// calculate the power in the order needed by the 3d-kernel structure
	//
	double *kin=Pk.extractArguments();
	double *Pkin=Pk.extractValues();
	double *Pkint=cpeds_interpolate(kin,Pkin,Pk.pointsCount(),kk,l,interpolationType);
	delete [] kin;
	delete [] Pkin;
	delete [] kk;
	
	
	
	//
	// generate kernel values
	//	
	l=0;
	for (long i = 0; i < nx; i++) {
		for (long j = 0; j < ny; j++) {
			for (long k = 0; k < nz; k++) {
				if (Pkint[l]>0)	f(i,j,k)[0]=sqrt(Pkint[l]); else f(i,j,k)[0]=0; // this is to get rid of possible residual oscillations in case of non-linear interpolation
				f(i,j,k)[1]=f(i,j,k)[0];
				l++;
			}
		}
	}
	delete [] Pkint;

	// remove the mean
	f(0,0,0)[0]=0; 
	f(0,0,0)[1]=0;
	
	antisymmetrize(); // populate the other half of the kernel that we didn't really walk over
	re2im();
	
	
	
//	allocFunctionSpace();
//	double kx,ky,kz,modk;
//	long nx=Nx();
//	long ny=Ny();
//	long nz=Nz();
//	//	long n=nx*ny*(nz+1);
//	
//	double* kk = new double[size()];
//	long l=0;
//	long l1=0,l2=0;
//	long Nbelow=0, Nabove=0;
//	double kmin,kmax;
//	Pk.sortFunctionArgAscending();
//	kmin=Pk.getMinArg();
//	kmax=Pk.getMaxArg();
//	//	double Kmax=double(nx-1)/double(_param.sizeX);
//	
//	msgs->say("Pre-calculating the k-modes and interpolating the power spectrum",High);
//	for (long i = 0; i < nx; i++) {
//		for (long j = 0; j < ny; j++) {
//			for (long k = 0; k < nz; k++) {	
//
//				kx=double(i)/double(_param.sizeX)-_param.Nx/(2.0*_param.sizeX);
//				ky=double(j)/double(_param.sizeY)-_param.Ny/(2.0*_param.sizeY);
//				kz=double(k)/double(_param.sizeZ)-_param.Nz/(2.0*_param.sizeZ);
//				modk=sqrt(kx*kx + ky*ky + kz*kz);
//				
//				if (modk<kmin) { msgs->warning("P(k) is not sufficiently tabulated for small k -- using the closest, minimal k value: "+msgs->toStr(kmin)+"  (value requested: "+msgs->toStr(modk)+")",High); modk=kmin; Nbelow++; }
//				if (modk>kmax) { msgs->warning("P(k) is not sufficiently tabulated for large k -- using the closest, maximal k value"+msgs->toStr(kmax)+"  (value requested: "+msgs->toStr(modk)+")",High); modk=kmax; Nabove++; }
//				kk[l]=modk;
//				l++;
//			}
//		}
//	}
//	printf("l: %li\n",l);
//	printf("Nbelow: %li\n",Nbelow);
//	printf("Nabove: %li\n",Nabove);
//	//
//	// calculate the power in the order needed by the 3d-kernel structure
//	//
//	double *kin=Pk.extractArguments();
//	double *Pkin=Pk.extractValues();
//	double *Pkint=cpeds_interpolate(kin,Pkin,Pk.pointsCount(),kk,l,interpolationType);
//	delete [] kin;
//	delete [] Pkin;
//	delete [] kk;
//	
//	
//	
//	//
//	// generate kernel values
//	//	
//	l=0;
//	for (long i = 0; i < nx; i++) {
//		for (long j = 0; j < ny; j++) {
//			for (long k = 0; k < nz; k++) {
//				//				if (Pkint[l]>0)	f(i,j,k)[0]=sqrt(Pkint[l]/2); else f(i,j,k)[0]=0; // this is to get rid of possible residual oscillations in case of non-linear interpolation
//				if (Pkint[l]>0)	f(i,j,k)[0]=sqrt(Pkint[l]); else f(i,j,k)[0]=0; // this is to get rid of possible residual oscillations in case of non-linear interpolation
//				f(i,j,k)[1]=f(i,j,k)[0];
//				l++;
//			}
//		}
//	}
//	delete [] Pkint;
//	
//	//#ifndef NO_HDF5
//	//	saveHDF5("kernel.tmp.hdf5","ker");
//	//#endif
//	shift(nx/2, ny/2, nz/2);
//	// remove the mean
//	f(0,0,0)[0]=0; 
//	f(0,0,0)[1]=0;
//	//#ifndef NO_HDF5
//	//	saveHDF5("kernel.tmp.hdf5","kersh");
//	//#endif
//	//	divide(getMaxValue());

	
	
	// OLD BELOW
	//	antisymmetrize();
	//	saveHDF5("kernel.tmp.hdf5","kershan");
	//	absoluteValue();
	//	cpeds_matrix_save(getSlice(0,nx/2,0),"kerFFT.mid");
	//	cpeds_matrix_save(getSlice(0,0,0),"kerFFT.0");
	//	exit(0);
	
	//	for (long i = 0; i < nx; i++) {
	//		for (long j = 0; j < ny; j++) {
	//			for (long k = 0; k <= nz; k++) {	
	//				
	//				modk=kk[l];
	//				if (Pkint[l]<0) {  // this is to get rid of possible residual oscillations from non-linear interpolation
	//					msgs->warning("Pkint is "+msgs->toStr(Pkint[l])+" < 0. Assuming power 0 for k-mode: (i,j,k)= ("+msgs->toStr(i)+", "+msgs->toStr(j)+", "+msgs->toStr(k)+")",High);
	//					Pkint[l]=0;
	//				}
	//
	//				if (modk<Kmax) delta=sqrt(Pkint[l])/sqrt2; //sqrt(Pk.f(modk))/sqrt2;
	//				else {
	//					delta=0;
	//				}
	//				l++;
	//			}
	//		}
	//	}
	
	return *this;
}
/***************************************************************************************/
void mscsFunction3dregc::setSize(long Nx, long Ny, long Nz, double dx, double dy, double dz, double x0, double y0, double z0) {
	_param.dx=dx; _param.dxo2=dx/2.0;
	_param.dy=dy; _param.dyo2=dy/2.0;
	_param.dz=dz; _param.dzo2=dz/2.0;
	_param.x0=x0;
	_param.y0=y0;
	_param.z0=z0;
	
	_param.xMax=_param.x0+Nx*_param.dx;
	_param.yMax=_param.y0+Ny*_param.dy;
	_param.zMax=_param.z0+Nz*_param.dz;
	calcLengths();
	
#ifdef DEBUG
	printf("mscsFunction3dregc::setSize>> setting new size Nx, Ny, Nz: %li %li %li\n",Nx,Ny,Nz);
#endif
	setSize(Nx,Ny,Nz);
	
}

/***************************************************************************************/
void mscsFunction3dregc::setSizeRange(long Nx, long Ny, long Nz, double x0, double y0, double z0, double xMax, double yMax, double zMax) {
	_param.x0=x0;
	_param.y0=y0;
	_param.z0=z0;
	
	_param.xMax=xMax;
	_param.yMax=yMax;
	_param.zMax=zMax;
	
	_param.dx=(_param.xMax-_param.x0)/double(Nx); _param.dxo2=_param.dx/2.0;
	_param.dy=(_param.yMax-_param.y0)/double(Ny); _param.dyo2=_param.dy/2.0;
	_param.dz=(_param.zMax-_param.z0)/double(Nz); _param.dzo2=_param.dz/2.0;
	
#ifdef DEBUG
	printf("mscsFunction3dregc::setSize>> setting new size Nx, Ny, Nz: %li %li %li\n",Nx,Ny,Nz);
#endif
	setSize(Nx,Ny,Nz);	
#ifdef DEBUG
	printf("mscsFunction3dregc::setSize>> setting new size Nx, Ny, Nz: %li %li %li\n",_param.Nx,_param.Ny,_param.Nz);
#endif
	
	calcLengths();
}
/***************************************************************************************/
void mscsFunction3dregc::calcLengths() {
	_param.sizeX=_param.xMax-_param.x0;
	_param.sizeY=_param.yMax-_param.y0;
	_param.sizeZ=_param.zMax-_param.z0;
}
/***************************************************************************************/
void mscsFunction3dregc::setSize(long Nx, long Ny, long Nz) {
	_param.Nx=Nx;
	_param.Ny=Ny;
	_param.Nz=Nz;
	
	//	_param.sizeX=_param.Nx*_param.dx;
	//	_param.sizeY=_param.Ny*_param.dy;
	//	_param.sizeZ=_param.Nz*_param.dz;
	//	_param.sizeX=Nx;
	//	_param.sizeY=Ny;
	//	_param.sizeZ=Nz;
	
	
	_param.NYZ=_param.Ny*_param.Nz;
	_param.NXY=_param.Nx*_param.Ny;
	_param.NXYZ=_param.NXY*_param.Nz;
}
/***************************************************************************************/
void mscsFunction3dregc::setSize(long Nx, long Ny) {
	setSize(Nx,Ny,1);
}
/***************************************************************************************/
void mscsFunction3dregc::moveBox(double dx, double dy, double dz) {
	_param.x0+=dx;
	_param.y0+=dy;
	_param.z0+=dz;
	
	_param.xMax+=dx;
	_param.yMax+=dy;
	_param.zMax+=dz;
	
}
/***************************************************************************************/
void mscsFunction3dregc::setDo2(double dxo2, double dyo2, double dzo2) {
	_param.dxo2=dxo2;
	_param.dyo2=dyo2;
	_param.dzo2=dzo2;
}

/***************************************************************************************/
matrix<double> mscsFunction3dregc::getSlice(long plane, long coord, long part) const {
	
	matrix<double> m;
	long n1,n2;
	
	switch (plane) {
		case 0:
			n1=Ny();
			n2=Nz();
			break;
		case 1:
			n1=Nx();
			n2=Nz();
			break;
		case 2:
			n1=Nx();
			n2=Ny();
			break;
	}
	if (plane>2 or plane < 0) { msgs->criticalError("mscsFunction3dregc::getSlice: the plane parameter value out of range",Top); }
	
	m.SetSize(n2,n1);
	for (long i1 = 0; i1 < n1; i1++) {
		for (long i2 = 0; i2 < n2; i2++) {
			//			printf("idx: %li, i1 %li, i2 %li, n1: %li, n2: %li, _param.NyZ: %li\n",ijk2idx(coord,i1,i2),i1,i2,n1,n2,_param.NyZ);
			
			if (plane==0) {
				if (part==0) m(i2,i1)=fRe(coord,i1,i2); 
				else  { 
					if (part==1) m(i2,i1)=fIm(coord,i1,i2);  
					else { m(i2,i1)=fRe(coord,i1,i2)*fIm(coord,i1,i2); }
				}
			}
			if (plane==1) {
				if (part==0) m(i2,i1)=fRe(i1,coord,i2); 
				else { 
					if (part==1) m(i2,i1)=fIm(i1,coord,i2); 
					else { m(i2,i1)=fRe(i1,coord,i2)*fIm(i1,coord,i2); }
				}
			}
			if (plane==2) {
				if (part==0) m(i2,i1)=fRe(i1,i2,coord); 
				else  {
					if (part==1) m(i2,i1)=fIm(i1,i2,coord); 
					else { m(i2,i1)=fRe(i1,i2,coord)*fIm(i1,i2,coord); }
				}
			}
		}
	}
	return m;
}
/***************************************************************************************/
double* mscsFunction3dregc::getSlice1D(long dir, long coord1, long coord2, long part) const {
	if (dir>2 or dir < 0) { msgs->criticalError("mscsFunction3dregc::getSlice: the dir parameter value out of range",Top); }
	long N=0;
	if (dir==0) { N=_param.Nx;}
	if (dir==1) { N=_param.Ny;}
	if (dir==2) { N=_param.Nz;}
	//	printf("N=%li\n",N);
	double* d = new double[N];
	
	if (part<2 and part>=0) {
		if (dir==0) { for (long i = 0; i < N; i++) { d[i]=f(i,coord1,coord2)[part]; }	}
		if (dir==1) { for (long i = 0; i < N; i++) { d[i]=f(coord2,i,coord1)[part]; }	}
		if (dir==2) { for (long i = 0; i < N; i++) { d[i]=f(coord1,coord2,i)[part]; }	}
		return d;
	}
	if (part==2) {
		if (dir==0) { for (long i = 0; i < N; i++) { d[i]=getX(i); }	}
		if (dir==1) { for (long i = 0; i < N; i++) { d[i]=getY(i); }	}
		if (dir==2) { for (long i = 0; i < N; i++) { d[i]=getZ(i); }	}
		return d;
	}
	return NULL;
}
/***************************************************************************************/
mscsFunction mscsFunction3dregc::getSlice1Dfn(long dir, long coord1, long coord2, long part) const {
	double *d=getSlice1D(dir,coord1,coord2,part);
	double *arg=getSlice1D(dir,coord1,coord2,2);
	long N=0;
	if (dir==0) { N=_param.Nx;}
	if (dir==1) { N=_param.Ny;}
	if (dir==2) { N=_param.Nz;}
	mscsFunction f("f",arg,d,N,Zero);
	delete [] d;
	delete [] arg;
	return f;
}
/***************************************************************************************/
mscsVector<double> mscsFunction3dregc::getLOSvec(double lx,double ly,double lz,double x0, double y0, double z0, double dt, double tmax, double* zmax, mscsVector<double>* X, mscsVector<double>* Y, mscsVector<double>* Z, string options) {
	msgs->warning("mscsFunction3dregc::getLOSvec>> has not been tested", Top);
	mscsVector<double> los;
	double t;
	
	t=0;
	double x,y,z;
	double Zcond=z0;
	
	if (zmax == NULL) Zcond=_param.zMax;
	else Zcond=*zmax;
	
	if (X!=NULL and Y!=NULL and Z!=NULL) {
		while (t<=tmax and z<=Zcond) {
			x=lx*t+x0;
			y=ly*t+y0;
			z=lz*t+z0;
			los.push_back(fReCoordPeriodic(x,y,z));
			X->push_back(x);
			Y->push_back(y);
			Z->push_back(z);
			t+=dt;
		}
	}	
	else {
		while (t<=tmax and z<=Zcond) {
			x=lx*t+x0;
			y=ly*t+y0;
			z=lz*t+z0;		
			los.push_back(fReCoordPeriodic(x,y,z));
			t+=dt;
		}
	}
	
	//	msgs->criticalError("mscsFunction3dregc::getLOSvec>> THIS IS NOT IMPLEMENTED YET",Top);
	return los;
}
/***************************************************************************************/
mscsVector<double> mscsFunction3dregc::getLOSvecZmax(double lx,double ly,double lz,double x0, double y0, double z0, double dt, double Zmax, mscsVector<double>* X, mscsVector<double>* Y, mscsVector<double>* Z, string options) {
	//	msgs->warning("mscsFunction3dregc::getLOSvecZmax>> has not been tested", Top);
	mscsVector<double> los;
	double t;
	
	t=0;
	double x,y,z;
	
	if (X!=NULL and Y!=NULL and Z!=NULL) {
		while (z<=Zmax) {
			x=lx*t+x0;
			y=ly*t+y0;
			z=lz*t+z0;
			los.push_back(fReCoordPeriodic(x,y,z));
			X->push_back(x);
			Y->push_back(y);
			Z->push_back(z);
			t+=dt;
		}
	}	
	else {
		while (z<=Zmax) {
			x=lx*t+x0;
			y=ly*t+y0;
			z=lz*t+z0;		
			los.push_back(fReCoordPeriodic(x,y,z));
			t+=dt;
		}
	}
	
	return los;	
}
/***************************************************************************************/
mscsVector<double> mscsFunction3dregc::getLOSvec(double lx,double ly,double lz,double x0, double y0, double z0, double dt, double Lmax, mscsVector<double>* X, mscsVector<double>* Y, mscsVector<double>* Z, string options) {
	msgs->warning("mscsFunction3dregc::getLOSvec>> has not been tested", Top);
	mscsVector<double> los;
	double t=0;
	double L=0;
	double x,y,z;
	cpedsPoint3D p1(x0,y0,z0),p2;
	
	if (X!=NULL and Y!=NULL and Z!=NULL) {
		while (L<=Lmax) {
			x=lx*t+x0;
			y=ly*t+y0;
			z=lz*t+z0;
			los.push_back(fReCoordPeriodic(x,y,z));
			X->push_back(x);
			Y->push_back(y);
			Z->push_back(z);
			t+=dt;
			p2.set(x,y,z);
			L=p1.dist(p2);
		}
	}	
	else {
		while (L<=Lmax) {
			x=lx*t+x0;
			y=ly*t+y0;
			z=lz*t+z0;		
			los.push_back(fReCoordPeriodic(x,y,z));
			t+=dt;
			L=p1.dist(p2);
		}
	}
	
	return los;
	
}
/***************************************************************************************/
void mscsFunction3dregc::printInfo() const {
	msgs->say("sizeX sizeY sizeZ: %lE %lE %lE",_param.sizeX,_param.sizeY,_param.sizeZ,Low);
	msgs->say("Nx Ny Nz: %li %li %li",_param.Nx,_param.Ny,_param.Nz,Low);	
	msgs->say("dx dy dz: %lE %lE %lE\n",_param.dx,_param.dy,_param.dz,Low);	
	msgs->say("x0 y0 z0: %lE %lE %lE\n",_param.x0,_param.y0,_param.z0,Low);	
	msgs->say("xMax yMax zMax: %lE %lE %lE\n",_param.xMax,_param.yMax,_param.zMax,Low);	
	msgs->say("DOMAIN RANGES:",Medium);
	msgs->say("X: %lE %lE\n",getMinX(),getMaxX(),Low);
	msgs->say("Y: %lE %lE\n",getMinY(),getMaxY(),Low);
	msgs->say("Z: %lE %lE\n",getMinZ(),getMaxZ(),Low);
}

/***************************************************************************************/
void mscsFunction3dregc::printFunction() const {
	for (long i = 0; i < Nx(); i++) {
		for (long j = 0; j < Ny(); j++) {
			for (long k = 0; k < Nz(); k++) {
				printf("i: %li, j: %li k: %li   Re: %.10lE Im: %.10lE\n", i,j,k,fRe(i,j,k),fIm(i,j,k));
			}
		}
	}
}

/***************************************************************************************/
void mscsFunction3dregc::antisymmetrize() {
	msgs->warning("mscsFunction3dregc::antisymmetrize() check this method for the case of even number of grid cells.", Top);
	
	if (lengthZ()==0) { // 2D
		for (long i = 0; i < Nx(); i++) {
			for (long j = 0; j <= Ny()/2; j++) {
				setf((Nx()-i) % Nx(),(Ny()-j) % Ny(), 0, fRe(i,j,0),-fIm(i,j,0));
				if ((Nx()-i) % Nx()==i and (Ny()-j) % Ny()==j) { fIm(i,j,0)=0.0; }
			}
		}	
	}
	else {
		for (long i = 0; i < Nx(); i++) {
			for (long j = 0; j < Ny(); j++) {
				for (long k = 0; k <= Nz()/2; k++) {
					setf((Nx()-i) % Nx(),(Ny()-j) % Ny(), (Nz()-k) % Nz(), fRe(i,j,k),-fIm(i,j,k));
					if ((Nx()-i) % Nx()==i and (Ny()-j) % Ny()==j and (Nz()-k) % Nz()==k) { fIm(i,j,k)=0.0; }
				}
			}
		}

	}
	// set the imaginary part at the Nyquist frequency to zero if the size of the transform is even.
	if (Nx()%2==0) { fIm(Nx()/2,0,0)=0; }
	if (Ny()%2==0) { fIm(0,Ny()/2,0)=0; }
	if (Nz()%2==0) { fIm(0,0,Nz()/2)=0; }
}
/***************************************************************************************/
void mscsFunction3dregc::flipXYZ() {
	mscsFunction3dregc tmp(*this);
	for (long i = 0; i < Nx(); i++) {
		for (long j = 0; j < Ny(); j++) {
			for (long k = 0; k < Nz(); k++) {
				setf(Nx()-i-1 ,Ny()-j-1, Nz()-k-1, tmp.fRe(i,j,k),-tmp.fIm(i,j,k));
			}
		}
	}
}
/***************************************************************************************/
void mscsFunction3dregc::shift(long nx, long ny, long nz) {
	fftw_complex* data=new fftw_complex[size()];
	nx = nx % _param.Nx;
	ny = ny % _param.Ny;
	nz = nz % _param.Nz;
	if (nx<0) nx=_param.Nx+nx;
	if (ny<0) ny=_param.Ny+ny;
	if (nz<0) nz=_param.Nz+nz;
	
	for (long i = 0; i < _param.Nx; i++) {
		for (long j = 0; j < _param.Ny; j++) {
			for (long k = 0; k < _param.Nz; k++) {
				data[ijk2idx((i+nx) % _param.Nx,(j+ny) % _param.Ny,(k+nz) % _param.Nz)][0]=f(i,j,k)[0];
				data[ijk2idx((i+nx) % _param.Nx,(j+ny) % _param.Ny,(k+nz) % _param.Nz)][1]=f(i,j,k)[1];
				//				printf("shift: %li %li %li,   i,j,k: %li %li %li --> %li %li %li\n",nx,ny,nz,i,j,k,(i+nx) % _param.Nx,(j+ny) % _param.Ny,(k+nz) % _param.Nz);
			}
		}
	}
	long n=size();
	dataHardCopy(data);
	//	for (long i = 0; i < n; i++) { _data[i][0]=data[i][0]; _data[i][1]=data[i][1]; }
	delete [] data;
	
}

/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::fft(bool fwd) {
	fft_c2c(fwd);
	
	//	if (fwd==false) {		cpeds_divide_value(double(getN()),_data,getN());	} // BLcomment (Nov 11, 2011, 10:07:24 AM): this is done inside of fft_c2c now
	//	if (fwd) multiply(getPixelVolume(0,0,0,false)); 
	/*
	 * Comment: the fft_c2c function is left unchanged so that FT can be done on raw function data without
	 * messing with the function ranges.
	 * But, since the function domain can in general be specified in various units (other than 1)
	 * when making FFT we need to account for the pixel volume so that forward and backward transform
	 * give correct results in case of eg. smoothing.
	 * 
	 * author: blew
	 * date: Apr 18, 2013 12:18:41 PM
	 *
	 */
	
	return *this;
}


/***************************************************************************************/
void mscsFunction3dregc::importFunction(double* a, long N, bool re) {
	if (re) {
		for (long i = 0; i < N; i++) {
			_data[i][0]=a[i];
		}
	}
	else {
		for (long i = 0; i < N; i++) {
			_data[i][1]=a[i];
		}
	}
}
/***************************************************************************************/
void mscsFunction3dregc::importFunction(fftw_complex* a, long N) {
	freeFunctionSpace();
	allocFunctionSpace(N);
	
	//	for (long i = 0; i < N; i++) {
	//		_data[i][0]=a[i][0];
	//		_data[i][1]=a[i][1];
	//	}
	dataHardCopy(a);
	delete [] a;
}
/***************************************************************************************/
void mscsFunction3dregc::importFunction(double* re, double* im, long N, bool del) {
	for (long i = 0; i < N; i++) {
		_data[i][0]=re[i];
		_data[i][1]=im[i];
	}
	if (del) { delete [] re; delete [] im; }
}
/***************************************************************************************/
void mscsFunction3dregc::importFunction(mscsFunction& inf) {
	double* re=inf.extractArguments();
	double* im=inf.extractValues();
	importFunction(re,im,inf.pointsCount(),true);
}
/***************************************************************************************/
mscsFunction3dregc mscsFunction3dregc::concatenate(mscsFunction3dregc& inf, int axis) {
	msgs->warning("concatenate is not full implemented, the domain ranges are not updated in the concatenated function",High);
	mscsFunction3dregc f;
	long nx,ny,nz,offset_i=0,offset_j=0,offset_k=0;
	
	if (size()==0) {f=inf; return f; }
	
	if (axis==0) {
		nx=Nx()+inf.Nx();
		ny=Ny();
		nz=Nz();
		offset_i=Nx();
		if (ny!=inf.Ny()) msgs->criticalError("mscsFunction3dregc mscsFunction3dregc::concatenate: ny!=inf.Ny()",Top);
		if (nz!=inf.Nz()) msgs->criticalError("mscsFunction3dregc mscsFunction3dregc::concatenate: nz!=inf.Nz()",Top);
	}
	if (axis==1) {
		nx=Nx();
		ny=Ny()+inf.Ny();
		nz=Nz();
		offset_j=Ny();
		if (nx!=inf.Nx()) msgs->criticalError("mscsFunction3dregc mscsFunction3dregc::concatenate: nx!=inf.Nx()",Top);
		if (nz!=inf.Nz()) msgs->criticalError("mscsFunction3dregc mscsFunction3dregc::concatenate: nz!=inf.Nz()",Top);
	}
	if (axis==2) {
		nx=Nx();
		ny=Ny();
		nz=Nz()+inf.Nz();
		offset_k=Nz();
		if (nx!=inf.Nx()) msgs->criticalError("mscsFunction3dregc mscsFunction3dregc::concatenate: nx!=inf.Nx()",Top);
		if (ny!=inf.Ny()) msgs->criticalError("mscsFunction3dregc mscsFunction3dregc::concatenate: ny!=inf.Ny()",Top);
	}
	
	f.setSize(nx,ny,nz);
	f.allocFunctionSpace();
	for (long i = 0; i < nx; i++) {
		for (long j = 0; j < ny; j++) {
			for (long k = 0; k < nz; k++) {
				if (i<Nx() and j<Ny() and k<Nz())
					f.setf(i,j,k,
							fRe(i,j,k),
							fIm(i,j,k));
				else {
					f.setf(i,j,k,
							inf.fRe((i-offset_i) % inf.Nx(),(j-offset_j) % inf.Ny(),(k-offset_k) % inf.Nz()),
							inf.fIm((i-offset_i) % inf.Nx(),(j-offset_j) % inf.Ny(),(k-offset_k) % inf.Nz())
					);
					//					printf("i: %li, j: %li, k: %li    i2: %li, j2: %li, k2: %li  \n",i,j,k,   (i-offset_i) % inf.Nx(), (j-offset_j) % inf.Ny(),  (k-offset_k) % inf.Nz());
				}
			}
		}
	}
	
	//	f.saveHDF5("concatenate-test.hdf5","joint0");
	
	return f;
}
/***************************************************************************************/
void mscsFunction3dregc::importSlice(double* a, long N, int slice, bool re) {
	long ist,ien;
	
	ist=slice*_param.NYZ;
	ien=ist+N;
	
	if (re) {
		for (long i = ist; i < ien; i++) {
			_data[i][0]=a[i-ist];
		}
	}
	else {
		for (long i = ist; i < ien; i++) {
			_data[i][1]=a[i-ist];
		}
	}
}
/***************************************************************************************/
double* mscsFunction3dregc::exportRe(string ordering) const {
	long N=getN();
	if (N==0) return NULL;
	
	double* t = new double[N];
	
	if (ordering=="XYZmajor") {
		for (long i = 0; i < N; i++) {
			t[i]=_data[i][0];
		}		
	}
	if (ordering=="YXZmajor") {
		long l=0;
		for (long j = 0; j < _param.Ny; j++) {
			for (long i = 0; i < _param.Nx; i++) {
				for (long k = 0; k < _param.Nz; k++) {
					t[l]=_data[i*_param.NYZ+j*_param.Nz+k][0];
					l++;
				}
			}
		}
	}
	if (ordering=="ZYXmajor") {
		long l=0;
		for (long k = 0; k < _param.Nz; k++) {
			for (long j = 0; j < _param.Ny; j++) {
					for (long i = 0; i < _param.Nx; i++) {
					t[l]=_data[i*_param.NYZ+j*_param.Nz+k][0];
					l++;
				}
			}
		}
	}
	
	return t;	
}
/***************************************************************************************/
double* mscsFunction3dregc::exportIm() const {
	long N=getN();
	if (N==0) return NULL;
	
	double* t = new double[N];
	
	for (long i = 0; i < N; i++) {
		t[i]=_data[i][1];
	}
	return t;	
}
/***************************************************************************************/
vector<cpedsPoint3D> mscsFunction3dregc::positions2Vec() {
	long nx=Nx();
	long ny=Ny();
	long nz=Nz();
	long N=getN();
	long l=0;
	vector<cpedsPoint3D> v(N);
	for (long i = 0; i < nx; i++) {
		for (long j = 0; j < ny; j++) {
			for (long k = 0; k < nz; k++) {
				v[l]=cpedsPoint3D(getX(i),getY(j),getZ(k));
				l++;
			}
		}
	}
	return v;
}
/***************************************************************************************/
cpedsPointSet3D mscsFunction3dregc::exportAsPoints(bool ZasVals, double mask) {
	long nx=Nx();
	long ny=Ny();
	long nz=Nz();
	long N=getN();
	cpedsPointSet3D ps;
	if (ZasVals) {
		for (long i = 0; i < nx; i++) {
			for (long j = 0; j < ny; j++) {
				for (long k = 0; k < nz; k++) {
					if (fIm(i,j,k)!=mask)
						ps.append(cpedsPoint3D(getX(i),getY(j),fRe(i,j,k)),fRe(i,j,k));
				}
			}
		}		
	}
	else {
		for (long i = 0; i < nx; i++) {
			for (long j = 0; j < ny; j++) {
				for (long k = 0; k < nz; k++) {
					if (fIm(i,j,k)!=mask)
						ps.append(cpedsPoint3D(getX(i),getY(j),getZ(k)),fRe(i,j,k));
				}
			}
		}
	}
	return ps;	
}
/***************************************************************************************/
void mscsFunction3dregc::dataHardCopy(fftw_complex* src) {
	memcpy(_data,src,size()*sizeof(fftw_complex));
}
/***************************************************************************************/
void mscsFunction3dregc::dataHardCopy(const fftw_complex* src) {
	memcpy(_data,src,size()*sizeof(fftw_complex));
}

/***************************************************************************************/
double mscsFunction3dregc::fxy(double x, double y, int k, int part, double* dfdx, double* dfdy, bool extrapolate) {
	
	//
	// derive the grid cells
	//
	long i=(x-_param.x0-_param.dxo2)/_param.dx;
	long j=(y-_param.y0-_param.dyo2)/_param.dy;
	//	long i=(x-_param.x0)/_param.dx;
	//	long j=(y-_param.y0)/_param.dy;
#ifdef DEBUG
	printf("derived i,j are: %i, %i, (x,y)=%lE, %lE\n",i,j,x,y);
#endif
	//
	// check the range
	//
	//	if (i<0) msgs->criticalError("fRe(double x, double y, k):: coordinate out of range (x)",High);
	//	if (i>=_param.Nx) msgs->criticalError("fRe(double x, double y, k):: coordinate out of range (x)",High);
	//	if (j<0) msgs->criticalError("fRe(double x, double y, k):: coordinate out of range (y)",High);
	//	if (j>=_param.Ny) msgs->criticalError("fRe(double x, double y, k):: coordinate out of range (y)",High);
	if (k<0) msgs->criticalError("fRe(double x, double y, k):: coordinate out of range (z)",High);
	if (k>=_param.Nz) msgs->criticalError("fRe(double x, double y, k):: coordinate out of range (z)",High);
	
	if (i<0) i=0;
	if (i>=_param.Nx-1) i=_param.Nx-2;
	if (j<0) j=0;
	if (j>=_param.Ny-1) j=_param.Ny-2;
	//	if (k<0) k=0;
	//	if (k>=_param.Nz-1) k=_param.Nz-2;
	//	printf("(x,y)=(%lf, %lf), corrected i,j, k are: %li, %li, %li\n",x,y,i,j,k);
	
	//
	// derive function values and periodic derivatives for the 4 surrounding grid points
	//
//	static double ff[4],ff1[4],ff2[4],ff12[4];
//#pragma omp threadprivate(ff,ff1,ff2,ff12)
	/*
	 * Comment: Removed static declaration because in that form it was not thread save anyway
	 * (only pointer was threadprivate) and implementation with static declarations of
	 * vector objects did not compile with gcc for some reason.
	 * 
	 * So, for the moment we will use local arrays at the penalty of 
	 * possibly slower execution time at frequent calls to this routine.
	 * 
	 * author: blew
	 * date: Jan 31, 2016 11:14:15 AM
	 *
	 */
	
	double ff[4],ff1[4],ff2[4],ff12[4];

//	extern static cpedsList<double> ff;
//	extern static vector<double> ff1;
//	extern static vector<double> ff2;
//	extern static vector<double> ff12;
//	vector<double> ff;
//	vector<double> ff1;
//	vector<double> ff2;
//	vector<double> ff12;
//	
//	ff.resize(4);
//	ff1.resize(4);
//	ff2.resize(4);
//	ff12.resize(4);
	if (part==0) {
		ff[0]=fRe(i, j, k);
		ff[1]=fRe(i+1, j, k);
		ff[2]=fRe(i+1, j+1, k);
		ff[3]=fRe(i, j+1, k);
	}
	else {
		ff[0]=fIm(i, j, k);
		ff[1]=fIm(i+1, j, k);
		ff[2]=fIm(i+1, j+1, k);
		ff[3]=fIm(i, j+1, k);		
	}
	get_derivativesXY(i,j,k,     ff1,   ff2,   ff12  , part);
	get_derivativesXY(i+1,j,k,   ff1+1, ff2+1, ff12+1, part);
	get_derivativesXY(i+1,j+1,k, ff1+2, ff2+2, ff12+2, part);
	get_derivativesXY(i,j+1,k,   ff1+3, ff2+3, ff12+3, part);

//	get_derivativesXY(i,j,k,     ff1.data(),   ff2.data(),   ff12.data()  , part);
//	get_derivativesXY(i+1,j,k,   ff1.data()+1, ff2.data()+1, ff12.data()+1, part);
//	get_derivativesXY(i+1,j+1,k, ff1.data()+2, ff2.data()+2, ff12.data()+2, part);
//	get_derivativesXY(i,j+1,k,   ff1.data()+3, ff2.data()+3, ff12.data()+3, part);

	//	for (long i = 0; i < 4; i++) {
	//		printf("point i=%li: f: %lf, f1: %lf, f2: %lf, f12: %lf\n",i,ff[i],ff1[i],ff2[i], ff12[i]);
	//	}
	
	static double fint, fint1,fint2;
#pragma omp threadprivate(fint, fint1,fint2)
//	cpeds_bicubic_interpolation(ff.data(),ff1.data(),ff2.data(),ff12.data(),getX(i),getX(i+1),getY(j),getY(j+1),x,y,fint,fint1,fint2);
	cpeds_bicubic_interpolation(ff,ff1,ff2,ff12,getX(i),getX(i+1),getY(j),getY(j+1),x,y,fint,fint1,fint2);
	
	if (dfdx!=NULL) *dfdx=fint1;
	if (dfdy!=NULL) *dfdy=fint2;
	return fint;
}
/* ******************************************************************************************** */
double mscsFunction3dregc::fxyLin(double x, double y, int k, int part, bool extrapolate) {
	//
	// derive the grid cells
	//
	long i=(x-_param.x0-_param.dxo2)/_param.dx;
	long j=(y-_param.y0-_param.dyo2)/_param.dy;
	//
	// check the range
	//
	if (k<0) msgs->criticalError("fRe(double x, double y, k):: coordinate out of range (z)",High);
	if (k>=_param.Nz) msgs->criticalError("fRe(double x, double y, k):: coordinate out of range (z)",High);
	
	if (i<0) i=0;
	if (i>_param.Nx-1) i=_param.Nx-1;
	if (j<0) j=0;
	if (j>_param.Ny-1) j=_param.Ny-1;

	
	double fint;
	double x1,x2,y1,y2;
	int i1,i2,j1,j2;

	if (x<getX(i))
		i1=(i-1+_param.Nx) % _param.Nx;
	else
		i1=i;
	i2=(i1+1) % _param.Nx;
	
	if (y<getY(j))
		j1=(j-1+_param.Ny) % _param.Ny;
	else
		j1=j;
	j2=(j1+1) % _param.Ny;
	
	x1=getX(i1);//-_param.dx;
	x2=getX(i2);//+_param.dx;
	y1=getY(j1);//-_param.dy;
	y2=getY(j2);//+_param.dy;
//	printf("x: %lf, y: %lf,    x1: %lf, x2: %lf, y1: %lf, y2: %lf\n",x,y,x1,x2,y1,y2);
	if (x1>=x2) x1=getX(i1)-_param.dx;
	if (y1>=y2) y1=getY(j1)-_param.dy;
//	printf("i: %li, j: %li,    i1: %li, i2: %li, j1: %li, j2: %li\n",i,j,i1,i2,j1,j2);
//	printf("x: %lf, y: %lf,    x1: %lf, x2: %lf, y1: %lf, y2: %lf\n",x,y,x1,x2,y1,y2);

	
	fint=cpeds_bilinear_interpolation(x1,x2,y1,y2,
			fRe(i1,j1,k),
			fRe(i1,j2,k),
			fRe(i2,j1,k),
			fRe(i2,j2,k),
			x,y);
	return fint;
	
}
/* ******************************************************************************************** */
double& mscsFunction3dregc::operator()(long i, long j) { return _data[ijk2idx(i,j,0)][0]; }
/* ******************************************************************************************** */
double mscsFunction3dregc::operator()(long i, long j) const { return _data[ijk2idx(i,j,0)][0]; }
/***************************************************************************************/
fftw_complex& mscsFunction3dregc::operator()(long i, long j, long k) { return _data[ijk2idx(i,j,k)]; }
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::operator=(const fftw_complex& v) { 
	long n=size();
	for (long i = 0; i < n; i++) { _data[i][0]=v[0]; _data[i][1]=v[1];	}
	return *this;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::operator=(const double& v) {
	setf(v,v);
	return *this;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::operator=(const mscsFunction3dregc& rhs) {
	this->mscsObject::operator =(rhs);
	long n=rhs.size();
	if (rhs.data()!=NULL) {
		if (size()!=rhs.size()) {
			_param=rhs.getFunctionParameters();
			allocFunctionSpace();
		}
		dataHardCopy(rhs.data());
		//		memcpy(_data,rhs.data(),size()*sizeof(fftw_complex));
		//		for (long i = 0; i < n; i++) { _data[i][0]=rhs[i][0]; _data[i][1]=rhs[i][1];	}
	}
	else {
		freeFunctionSpace();
	}
	_param=rhs.getFunctionParameters();
	
	return *this;
}
/***************************************************************************************/
//const mscsFunction3dregc& mscsFunction3dregc::operator=(const mscsFunction3dregc& rhs) {
//	this->mscsObject::operator =(rhs);
//	long n=rhs.size();
//	_param=rhs.getFunctionParameters();
//	allocFunctionSpace();
//	for (long i = 0; i < n; i++) { _data[i][0]=rhs[i][0]; _data[i][1]=rhs[i][1];	}
//
//	return *this;
//}
/***************************************************************************************/
fftw_complex& mscsFunction3dregc::operator[](long i) { return _data[i]; };	
const fftw_complex& mscsFunction3dregc::operator[](long i) const { return _data[i]; };	
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::operator+=(const mscsFunction3dregc& f) {
	long n=size();
	long i;
#pragma omp parallel for default(shared) private (i)
	for (i = 0; i < n; i++) { _data[i][0]+=f[i][0]; _data[i][1]+=f[i][1];	}
	return *this;	
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::operator-=(const mscsFunction3dregc& f) {
	long n=size();
	long i;
#pragma omp parallel for default(shared) private (i)
	for (i = 0; i < n; i++) { _data[i][0]-=f[i][0]; _data[i][1]-=f[i][1];	}
	return *this;	
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::operator*=(const mscsFunction3dregc& f) {
	long n=size();
	long i;
#pragma omp parallel for default(shared) private (i)
	for (long i = 0; i < n; i++) { _data[i][0]*=f[i][0]; _data[i][1]*=f[i][1];	}
	return *this;	
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::operator/=(const mscsFunction3dregc& f) {
	long n=size();
	long i;
#pragma omp parallel for default(shared) private (i)
	for (i = 0; i < n; i++) { _data[i][0]/=f[i][0]; _data[i][1]/=f[i][1];	}
	return *this;	
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::operator+=(const double v) {
	return add(v);
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::operator-=(const double v) {
	return subtract(v);
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::operator*=(const double v) {
	return multiply(v);
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::operator/=(const double v) {
	return divide(v);
}

/***************************************************************************************/
void mscsFunction3dregc::re2im() {
	long n=size();
	for (long i = 0; i < n; i++) { _data[i][1]=_data[i][0];	}	
}
/***************************************************************************************/
void mscsFunction3dregc::im2re() {
	long n=size();
	for (long i = 0; i < n; i++) { _data[i][0]=_data[i][1];	}	
}

/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::smooth3DGauss(double sxG,double syG,double szG) {
	//
	// make convolution kernel in real space
	//
	mscsFunction3dregc ker(_param.Nx,_param.Ny,_param.Nz, _param.sizeX/_param.Nx,_param.sizeY/_param.Ny,_param.sizeZ/_param.Nz);
	ker.setName("kernel");
	ker.mkGauss3D(sxG,syG,szG);
	// normalize kernel
	//	printf("kernel normalization factor: %lE\n",ker.integrateRe());
	ker/=ker.integrateRe(); 
	ker.shift(_param.Nx/2,_param.Ny/2, _param.Nz/2);
	//	printf("kernel normalization check: %lE\n",ker.integrateRe());
	
	// make kernel
	ker.fft(true);
	//	ker.re2im();
	ker.absoluteValue();
//	ker*=double(size());
	/*
	 * Comment: we have 2 forward transforms and 1 backward and since after forward transform 
	 * function is divided by the number of cells, so that after backward transform obtain the same values 
	 * as at the beginning, we need to multiply once here by the number of cells
	 * 
	 * author: blew
	 * date: Apr 18, 2013 11:44:09 AM
	 *
	 * REVISION Dec 29, 2015, 9:17:00 PM:
	 * since we have recently changed the fft normalization, we need to reflect that change here as well.
	 * Hence, we now multiply by sqrt(size())
	 */
	ker*=double(sqrt(double(getN())));
	
	//
	// convolve
	//
	fft(true);
	multiply(getPixelVolume(0,0,0,false));
	(*this)*=ker;
	//	f.antisymmetrize();
	fft(false);
	
	return *this;
}
/***************************************************************************************/
mscsFunction3dregc mscsFunction3dregc::convolve_fft(mscsFunction3dregc& ker) {
	mscsFunction3dregc conv=(*this);
	//
	// convolve
	//
	conv.fft(true);
	conv.multiply(conv.getPixelVolume(0,0,0,false));
	conv*=ker;
	//	f.antisymmetrize();
	conv.fft(false);
	return conv;
}
/***************************************************************************************/
mscsFunction3dregc mscsFunction3dregc::convolve(mscsFunction3dregc g) {
	mscsFunction3dregc fg=(*this);
	fg.setName("fg");
	mscsFunction3dregc tmp=(*this);
	//	tmp.setVerbosityLevel(Zero);
	tmp.setName("tmp");
	g.setVerbosityLevel(Zero);
	
	g.shift(g.Nx()/2,g.Ny()/2,g.Nz()/2);
	
	for (long i = 0; i < Nx(); i++) {
		msgs->say("done: %li/%li",i,Nx(),Low);
		for (long j = 0; j < Ny(); j++) {
			for (long k = 0; k < Nz(); k++) {
				
				tmp=g;
				tmp.shift(i,j,k);
				tmp.multiply(*this);
				fg(i,j,k)[0]=tmp.integrateRe();
				fg(i,j,k)[1]=0;
				
				//				for (long t1 = 0; t1 < Nx(); t1++) {
				//					for (long t2 = 0; t2 < Ny(); t2++) {
				//						for (long t3 = 0; t3 < Nz(); t3++) {
				//							fg(i,j,k)[0]+=f(t1,t2,t3)*g(i-t1,j-t2,k-t3);
				//						}
				//					}
				//				}
				
				
			}
		}
	}
	
	return fg;
}
/***************************************************************************************/
mscsFunction3dregc mscsFunction3dregc::RLdeconvolution(mscsFunction3dregc data, mscsFunction3dregc psf, long iter) {
	msgs->say("RL deconvolution",Medium);
	mscsFunction3dregc Ci,estimate,tmp,psfhat,R;
	estimate=data;
	R=data;
	estimate=0.5;
	double sum=0;

	
	//
	// real space implementation
	//
	// make kernels
//#define testwithHDF
#ifdef testwithHDF
	data.saveHDF5("RLestimateRS.hdf5","data");
	psfhat=psf;
//	psf.shift(psf.Nx()/2,psf.Ny()/2, psf.Nz()/2);
	psf.saveHDF5("RLestimateRS.hdf5","psf");
//	printf("psf integral: %lE\n",psf.integrateRe());

	psfhat.flipXYZ();
//	psfhat.shift(psfhat.Nx()/2,psfhat.Ny()/2, psfhat.Nz()/2);
	psfhat.saveHDF5("RLestimateRS.hdf5","psfhat");

	// make iterations
	mscsFunction conv;
	for (long i = 0; i < iter; i++) {
		msgs->say("iteration %li / %li",i+1,iter,Low);
		Ci=psf.convolve(estimate);
		Ci.saveHDF5("Ci.hdf5",msgs->toStr(i+1));
		tmp=data/Ci;
		estimate=estimate*psfhat.convolve(tmp);
		estimate.saveHDF5("RLestimateRS.hdf5",msgs->toStr(i+1));
		sum=tmp.sumRe();
		printf("sum: %lE\n",sum);
		conv.newPoint(double(i),sum);
	}

	conv.save("convergence");
#endif	


	
	// 
	// fourier space implementation
	//
#ifndef NO_HDF5
	data.saveHDF5("RLestimate.hdf5","data");
	psf.saveHDF5("RLestimate.hdf5","psf");
#endif	
	// make kernels
	psfhat=psf;
//	psf/=psf.integrateRe(); 
	psf.shift(psf.Nx()/2,psf.Ny()/2, psf.Nz()/2);
	psf.fft(true);
	psf.absoluteValue();
	psf*=double(psf.size()); 

	psfhat.flipXYZ();
#ifndef NO_HDF5
	psfhat.saveHDF5("RLestimate.hdf5","psfhat");
#endif	
//	psfhat/=psfhat.integrateRe(); 
	psfhat.shift(psfhat.Nx()/2,psfhat.Ny()/2, psfhat.Nz()/2);
	psfhat.fft(true);
	psfhat.absoluteValue();
	psfhat*=double(psfhat.size()); 

	// make iterations
	for (long i = 0; i < iter; i++) {
		msgs->say("iteration %li / %li",i+1,iter,Low);
		Ci=estimate.convolve_fft(psf);
		tmp=data/Ci;
		estimate=estimate*tmp.convolve_fft(psfhat);
#ifndef NO_HDF5
		estimate.saveHDF5("RLestimate.hdf5",msgs->toStr(i+1));
#endif	
	}

	return estimate;
}
/***************************************************************************************/
//mscsFunction3dregc mscsFunction3dregc::CLEANdevonvolution(mscsFunction3dregc data, mscsFunction3dregc psf, long MaxIter, double loopGain) {
//	
//}

/***************************************************************************************/
mscsFunction3dregc mscsFunction3dregc::derivative(int dir, int part) {
	long nx=Nx();
	long ny=Ny();
	long nz=Nz();
	long n;
	double dx=getDx();
	double dy=getDy();
	double dz=getDz();
	mscsFunction3dregc fp(nx,ny,nz,dx,dy,dz);
	double* t;
	double h12;
	
	if (dir==0) {
		t = new double[nx];
		h12=dx*12;
		n=nx-2;
		
		for (long j = 0; j < ny; j++) {
			for (long k = 0; k < nz; k++) {
				
				// copy the column of data for speed
				if (part==0) { 
					for (long i = 0; i < nx; i++) {	t[i]=fRe(i,j,k); }
					
					
					// calculate derivative for 0'th cell
					fp.fRe(0,j,k)= (-t[2]+8*t[1]-8*t[nx-1]+t[nx-2])/h12;
					// calculate derivative for 1'st cell
					fp.fRe(1,j,k)= (-t[3]+8*t[2]-8*t[0]+t[nx-1])/h12;
					
					// calculate derivative for 2'nd through N-3'th cell
					for (long i = 2; i < n; i++) {
						fp.fRe(i,j,k)= (-t[i+2]+8*t[i+1]-8*t[i-1]+t[i-2])/h12;
					}
					
					// calculate derivative for N-2'th cell
					fp.fRe(nx-2,j,k)= (-t[0]+8*t[nx-1]-8*t[nx-3]+t[nx-4])/h12;
					// calculate derivative for N-1'th cell
					fp.fRe(nx-1,j,k)= (-t[1]+8*t[0]-8*t[nx-2]+t[nx-3])/h12;
				}
				
				if (part==1) { 
					for (long i = 0; i < nx; i++) {	t[i]=fIm(i,j,k); }
					
					// calculate derivative for 0'th cell
					fp.fIm(0,j,k)= (-t[2]+8*t[1]-8*t[nx-1]+t[nx-2])/h12;
					// calculate derivative for 1'st cell
					fp.fIm(1,j,k)= (-t[3]+8*t[2]-8*t[0]+t[nx-1])/h12;
					
					// calculate derivative for 2'nd through N-3'th cell
					for (long i = 2; i < n; i++) {
						fp.fIm(i,j,k)= (-t[i+2]+8*t[i+1]-8*t[i-1]+t[i-2])/h12;
					}
					
					// calculate derivative for N-2'th cell
					fp.fIm(nx-2,j,k)= (-t[0]+8*t[nx-1]-8*t[nx-3]+t[nx-4])/h12;
					// calculate derivative for N-1'th cell
					fp.fIm(nx-1,j,k)= (-t[1]+8*t[0]-8*t[nx-2]+t[nx-3])/h12;
				}
				
			}
		}
		delete [] t;
	}
	
	if (dir==1) {
		t = new double[ny];
		h12=dy*12;
		n=ny-2;
		
		for (long i = 0; i < nx; i++) {
			for (long k = 0; k < nz; k++) {
				
				// copy the column of data for speed
				if (part==0) { 
					for (long j = 0; j < ny; j++) {	t[j]=fRe(i,j,k); }
					
					
					// calculate derivative for 0'th cell
					fp.fRe(i,0,k)= (-t[2]+8*t[1]-8*t[ny-1]+t[ny-2])/h12;
					// calculate derivative for 1'st cell
					fp.fRe(i,1,k)= (-t[3]+8*t[2]-8*t[0]+t[ny-1])/h12;
					
					// calculate derivative for 2'nd through N-3'th cell
					for (long j = 2; j < n; j++) {
						fp.fRe(i,j,k)= (-t[j+2]+8*t[j+1]-8*t[j-1]+t[j-2])/h12;
					}
					
					// calculate derivative for N-2'th cell
					fp.fRe(i,ny-2,k)= (-t[0]+8*t[ny-1]-8*t[ny-3]+t[ny-4])/h12;
					// calculate derivative for N-1'th cell
					fp.fRe(i,ny-1,k)= (-t[1]+8*t[0]-8*t[ny-2]+t[ny-3])/h12;
				}
				
				if (part==1) { 
					for (long j = 0; j < ny; j++) {	t[j]=fIm(i,j,k); }
					
					
					// calculate derivative for 0'th cell
					fp.fIm(i,0,k)= (-t[2]+8*t[1]-8*t[ny-1]+t[ny-2])/h12;
					// calculate derivative for 1'st cell
					fp.fIm(i,1,k)= (-t[3]+8*t[2]-8*t[0]+t[ny-1])/h12;
					
					// calculate derivative for 2'nd through N-3'th cell
					for (long j = 2; j < n; j++) {
						fp.fIm(i,j,k)= (-t[j+2]+8*t[j+1]-8*t[j-1]+t[j-2])/h12;
					}
					
					// calculate derivative for N-2'th cell
					fp.fIm(i,ny-2,k)= (-t[0]+8*t[ny-1]-8*t[ny-3]+t[ny-4])/h12;
					// calculate derivative for N-1'th cell
					fp.fIm(i,ny-1,k)= (-t[1]+8*t[0]-8*t[ny-2]+t[ny-3])/h12;
				}
				
			}
		}
		delete [] t;
	}
	
	if (dir==2) {
		t = new double[nz];
		h12=dz*12;
		n=nz-2;
		
		for (long i = 0; i < nx; i++) {
			for (long j = 0; j < ny; j++) {
				
				// copy the column of data for speed
				if (part==0) { 
					for (long k = 0; k < nz; k++) {	t[k]=fRe(i,j,k); }
					
					
					// calculate derivative for 0'th cell
					fp.fRe(i,j,0)= (-t[2]+8*t[1]-8*t[nz-1]+t[nz-2])/h12;
					// calculate derivative for 1'st cell
					fp.fRe(i,j,1)= (-t[3]+8*t[2]-8*t[0]+t[nz-1])/h12;
					
					// calculate derivative for 2'nd through N-3'th cell
					for (long k = 2; k < n; k++) {
						fp.fRe(i,j,k)= (-t[k+2]+8*t[k+1]-8*t[k-1]+t[k-2])/h12;
					}
					
					// calculate derivative for N-2'th cell
					fp.fRe(i,j,nz-2)= (-t[0]+8*t[nz-1]-8*t[nz-3]+t[nz-4])/h12;
					// calculate derivative for N-1'th cell
					fp.fRe(i,j,nz-1)= (-t[1]+8*t[0]-8*t[nz-2]+t[nz-3])/h12;
				}
				
				if (part==1) { 
					for (long k = 0; k < nz; k++) {	t[k]=fIm(i,j,k); }
					
					
					// calculate derivative for 0'th cell
					fp.fIm(i,j,0)= (-t[2]+8*t[1]-8*t[nz-1]+t[nz-2])/h12;
					// calculate derivative for 1'st cell
					fp.fIm(i,j,1)= (-t[3]+8*t[2]-8*t[0]+t[nz-1])/h12;
					
					// calculate derivative for 2'nd through N-3'th cell
					for (long k = 2; k < n; k++) {
						fp.fIm(i,j,k)= (-t[k+2]+8*t[k+1]-8*t[k-1]+t[k-2])/h12;
					}
					
					// calculate derivative for N-2'th cell
					fp.fIm(i,j,nz-2)= (-t[0]+8*t[nz-1]-8*t[nz-3]+t[nz-4])/h12;
					// calculate derivative for N-1'th cell
					fp.fIm(i,j,nz-1)= (-t[1]+8*t[0]-8*t[nz-2]+t[nz-3])/h12;
				}
				
			}
		}
		delete [] t;
	}
	
	return fp;
}
/***************************************************************************************/
void mscsFunction3dregc::get_derivativesXY(int i, int j, int k, double* fi, double* fj, double* fij, int part) {
//	static long n[2];
//	static double t[5][5];
//	static double hx12,hy12,hxy12;
//	static long ii,jj,jjj;
//#pragma omp threadprivate(n, t,hx12,hy12,hxy12,ii,jj,jjj)
	/*
	 * Comment: Removed static declaration because in that form it was not thread save anyway
	 * (only pointer was threadprivate) and implementation with static declarations of
	 * vector objects did not compile with gcc for some reason.
	 * 
	 * So, for the moment we will use local arrays at the penalty of 
	 * possibly slower execution time at frequent calls to this routine.
	 * 
	 * author: blew
	 * date: Jan 31, 2016 11:14:15 AM
	 *
	 */
	long n[2];
	double t[5][5];
	static double hx12,hy12,hxy12;
	static long ii,jj,jjj;
#pragma omp threadprivate(hx12,hy12,hxy12,ii,jj,jjj)

	hx12=_param.dx*12;
	hy12=_param.dy*12;
	hxy12=hx12*hy12;
	n[0]=_param.Nx-2;
	n[1]=_param.Ny-2;
	
	
	// copy a patch of of data 5x5 in size with periodic boundary conditions
	if (part==0) { 
		for (jj = 0; jj < 5; jj++) {
			jjj=(jj+j-2+_param.Ny) % _param.Ny;
			for (ii = 0; ii < 5; ii++) {	t[ii][jj]=fRe((ii+i-2+_param.Nx) % _param.Nx,jjj,k); }
		}
	}
	if (part==1) { 
		for (jj = 0; jj < 5; jj++) {
			jjj=(jj+j-2+_param.Ny) % _param.Ny;
			for (ii = 0; ii < 5; ii++) {	t[ii][jj]=fIm((ii+i-2+_param.Nx) % _param.Nx,jjj,k); }			
		}
	}
		
	if (_param.derivativeXperiodic==0) { // remove data that stand out
		if (i<2 or i>=_param.Nx-2) {
			if (i<2) {
				for (unsigned long ii = 0; ii < 2-i; ii++) {
					for (unsigned long jj = 0; jj < 5; jj++) {
						t[ii][jj]=0;
					}
				}
			}
			else {
				for (unsigned long ii = _param.Nx-i+2; ii < 5; ii++) {
					for (unsigned long jj = 0; jj < 5; jj++) {
						t[ii][jj]=0;
					}
				}			
			}	
		}	
	}
	if (_param.derivativeYperiodic==0) { // remove data that stand out
		if (j<2 or j>=_param.Ny-2) {
			if (j<2) {
				for (unsigned long jj = 0; jj < 2-j; jj++) {
					for (unsigned long ii = 0; ii < 5; ii++) {
						t[ii][jj]=0;
					}
				}
			}
			else {
				for (unsigned long jj = _param.Ny-j+2; jj < 5; jj++) {
					for (unsigned long ii = 0; ii < 5; ii++) {
						t[ii][jj]=0;
					}
				}			
			}	
		}	
	}
	
	// calculate derivative at cell i,j
	jj=2;
	*fi  = (-t[4][jj]+8*t[3][jj]-8*t[1][jj]+t[0][jj])/hx12; // (df/dx)_2
	*fj  = (-t[jj][4]+8*t[jj][3]-8*t[jj][1]+t[jj][0])/hy12; // (df/dy)_2
	
	for (jj = 0; jj < 5; jj++) { // (df/dx)_i
		t[2][jj]  = (-t[4][jj]+8*t[3][jj]-8*t[1][jj]+t[0][jj])/hx12; // df/dx
	}
	*fij  = (-t[2][4]+8*t[2][3]-8*t[2][1]+t[2][0])/hy12; // d2f/dx/dy
	
}
/***************************************************************************************/
double mscsFunction3dregc::integrateRe() {
	return sumRe()*getPixelVolume(0,0,0,false);
}
/***************************************************************************************/
mscsFunction3dregc mscsFunction3dregc::integrateRe(long axis) {
	mscsFunction3dregc fint;
	if (size()==0) return fint;
	
	double tmpint=0;
	if (axis==0) {
		fint.setSizeRange(1,Ny(),Nz(),0,getDy(),getDz(),0,getMinY(),getMinZ());
		fint.allocFunctionSpace();
		fint.setf(0,0);
		for (long j = 0; j < Ny(); j++) {
			for (long k = 0; k < Nz(); k++) {
				tmpint=0;
				for (long i = 0; i < Nx(); i++) {
					tmpint+=fRe(i,j,k);	
				}
				fint.fRe(0,j,k)+=tmpint*getDx();
			}
		}
	}
	if (axis==1) {
		fint.setSize(Nx(),1,Nz(),getDx(),0,getDz(),getMinX(),0,getMinZ());
		fint.allocFunctionSpace();
		fint.setf(0,0);
		for (long i = 0; i < Nx(); i++) {
			for (long k = 0; k < Nz(); k++) {
				tmpint=0;
				for (long j = 0; j < Ny(); j++) {
					tmpint+=fRe(i,j,k);					
				}
				fint.fRe(i,0,k)+=tmpint*getDy();
			}
		}
	}
	if (axis==2) {
		fint.setSize(Nx(),Ny(),1,getDx(),getDy(),0,getMinX(),getMinY(),0);
		fint.allocFunctionSpace();
		fint.setf(0,0);
		if (Nz()>1) { // trapezoid rule
			for (long i = 0; i < Nx(); i++) {
				for (long j = 0; j < Ny(); j++) {
					tmpint=0;
					for (long k = 0; k < Nz()-1; k++) {
						tmpint+=(fRe(i,j,k)+fRe(i,j,k+1))/2.0;
					}
					fint.fRe(i,j,0)+=tmpint*getDz();
				}
			}			
		}
		else {
			for (long i = 0; i < Nx(); i++) {
				for (long j = 0; j < Ny(); j++) {
					tmpint=0;
					for (long k = 0; k < Nz(); k++) {
						tmpint+=fRe(i,j,k);
					}
					fint.fRe(i,j,0)+=tmpint*getDz();
				}
			}			
		}
	
	}
	return fint;
}
mscsFunction3dregc& mscsFunction3dregc::normalizeRe() {
	divide(integrateRe());
	return *this;
}
/***************************************************************************************/
double mscsFunction3dregc::sumRe() {
	double res=0;
	for (long i = 0; i < size(); i++) {
		res+=fRe(i);
	}
	return res;
}
/***************************************************************************************/
double mscsFunction3dregc::averageRe() {
	return sumRe()/getFunctionParameters().NXYZ;
}
/***************************************************************************************/
double mscsFunction3dregc::varianceRe() {
	long n=size();
	double v,tmp,av=averageRe();
	v=0;
	for (long i = 0; i < n; i++) { 
		tmp=_data[i][0]-av;
		v+=tmp*tmp;
	}
	return v/(n-1);
}
/***************************************************************************************/
double mscsFunction3dregc::RMSre() {
	long n=size();
	double v,tmp;
	v=0;
	for (long i = 0; i < n; i++) { 
		tmp=_data[i][0];
		v+=tmp*tmp;
	}
	return sqrt(v/n);
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::absoluteValue() {
	long n=size();
	double absVal;
	for (long i = 0; i < n; i++) { 
		absVal=sqrt(_data[i][0]*_data[i][0]+_data[i][1]*_data[i][1]);	
		_data[i][0]=absVal; _data[i][1]=absVal;
	}
	return *this;			
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::multiply(const mscsFunction3dregc& g) {
	long n=size();
	for (long i = 0; i < n; i++) { _data[i][0]*=g[i][0]; _data[i][1]*=g[i][1];	}
	return *this;		
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::multiply(const double v) {
	long n=size();
	for (long i = 0; i < n; i++) { _data[i][0]*=v; _data[i][1]*=v;	}
	return *this;		
}
mscsFunction3dregc& mscsFunction3dregc::divide(const double v) {
	long n=size();
	for (long i = 0; i < n; i++) { _data[i][0]/=v; _data[i][1]/=v;	}
	return *this;		
}
mscsFunction3dregc& mscsFunction3dregc::divide(const mscsFunction3dregc& f, int part) {
	long n=size();
	if (part==0 or part==2)
		for (long i = 0; i < n; i++) { if (f.fRe(i)!=0) _data[i][0]/=f.fRe(i); }
	if (part==1 or part==2)
		for (long i = 0; i < n; i++) { if (f.fIm(i)!=0) _data[i][1]/=f.fIm(i); }

	return *this;			
}
mscsFunction3dregc& mscsFunction3dregc::subtract(const double v) {
	long n=size();
	for (long i = 0; i < n; i++) { _data[i][0]-=v; _data[i][1]-=v;	}
	return *this;		
}
mscsFunction3dregc& mscsFunction3dregc::add(const double v) {
	long n=size();
	for (long i = 0; i < n; i++) { _data[i][0]+=v; _data[i][1]+=v;	}
	return *this;		
}
/***************************************************************************************/
mscsFunction3dregc mscsFunction3dregc::cutAwayBlock(long iSt, long jSt, long kSt,long iEn, long jEn, long kEn) {
	return cutAwayBlock_priv(iSt,jSt,kSt,iEn,jEn,kEn);
}
/***************************************************************************************/
mscsFunction3dregc mscsFunction3dregc::pasteAdd(mscsFunction3dregc& block, long iSt, long jSt, long kSt) {
	for (long i = iSt; i < iSt+block.Nx(); i++) {
		for (long j = jSt; j < jSt+block.Ny(); j++) {
			for (long k = kSt; k < kSt+block.Nz(); k++) {
				fRe(i,j,k)+=block.fRe(i-iSt,j-jSt,k-kSt);
				fIm(i,j,k)+=block.fIm(i-iSt,j-jSt,k-kSt);
			}
		}
	}
	return *this;
}
/***************************************************************************************/
//const mscsFunction3dregc mscsFunction3dregc::get_cutAwayBlock(long iSt, long jSt, long kSt,long iEn, long jEn, long kEn) {
//	return cutAwayBlock_priv(iSt,jSt,kSt,iEn,jEn,kEn);
//}
/***************************************************************************************/
mscsFunction3dregc mscsFunction3dregc::cutAwayBlock_priv(long iSt, long jSt, long kSt,long iEn, long jEn, long kEn) {
	if (iSt<0) iSt=0;
	if (jSt<0) jSt=0;
	if (kSt<0) kSt=0;
	if (iEn>=Nx()) iEn=Nx()-1;
	if (jEn>=Ny()) jEn=Ny()-1;
	if (kEn>=Nz()) kEn=Nz()-1;
	
	long nx=iEn-iSt+1;
	long ny=jEn-jSt+1;
	long nz=kEn-kSt+1;
	long newSize=nx*ny*nz;
	
	mscsFunction3dregc cut("", getVerbosityLevel());
	if (newSize==0) return cut;
	
	fftw_complex* data = new fftw_complex[newSize];
	long idx=0;
	for (long i = iSt; i <= iEn; i++) {
		for (long j = jSt; j <= jEn; j++) {
			for (long k = kSt; k <= kEn; k++) {
				data[idx][0]=fRe(i,j,k);
				data[idx][1]=fIm(i,j,k);
				idx++;
			}
		}
	}
	//	mscsFunction3dregc cut(nx,ny,nz,getDx(),getDy(),getDz(), getX(iSt)-_param.dxo2, getY(jSt)-_param.dyo2, getZ(kSt)-_param.dzo2,getVerbosityLevel());
	cut.setSize(nx,ny,nz,getDx(),getDy(),getDz(),getX(iSt)-_param.dxo2, getY(jSt)-_param.dyo2, getZ(kSt)-_param.dzo2);
	cut.importFunction(data,newSize);
	return cut;
}
/***************************************************************************************/
double mscsFunction3dregc::getPixelVolume(long i, long j, long k, bool spherical) {
	double dv=0;
	double y=getY(j);
	dv=getDx();	
	if (spherical) dv*=cos(PIsnd-(y+_param.dyo2))-cos(PIsnd-(y-_param.dyo2)); else dv*=getDy();
	if (getDz()>0) dv*=getDz();
	
	return dv;
}
/***************************************************************************************/
double mscsFunction3dregc::getVolume(bool spherical) {
	//	printInfo();
	double dv=lengthX();
	//	printf("V: %lE\n",dv);
	if (spherical) dv*=cos(PIsnd-(getMaxY()))-cos(PIsnd-(getMinY())); else dv*=lengthY();
	//	printf("V: %lE\n",dv);
	if (Nz()>1) dv*=lengthZ();
	//	printf("V: %lE\n",dv);
	
	return dv;	
}
/***************************************************************************************/
double mscsFunction3dregc::getMaxValue(long* iMax, bool Re) const {
	double vMax;
	long imax;
	if (Re) {
		vMax=fRe(0);
		for (long i = 0; i < _param.NXYZ; i++) {
			if (_data[i][0] > vMax) { vMax=_data[i][0]; imax=i; }
		}
	}
	else {
		vMax=fIm(0);
		for (long i = 0; i < _param.NXYZ; i++) {
			if (_data[i][1] > vMax) { vMax=_data[i][1]; imax=i; }
		}		
	}
	
	if (iMax!=NULL) *iMax=imax;
	return vMax;
}
/***************************************************************************************/
double mscsFunction3dregc::getMinValue(long* iMin, bool Re) const {
	double vMin;
	long imin;
	if (Re) {
		vMin=fRe(0);
		for (long i = 0; i < _param.NXYZ; i++) {
			if (_data[i][0] < vMin) { vMin=_data[i][0]; imin=i; }
		}
	}
	else {
		vMin=fIm(0);
		for (long i = 0; i < _param.NXYZ; i++) {
			if (_data[i][1] < vMin) { vMin=_data[i][1]; imin=i; }
		}		
	}
	if (iMin!=NULL) *iMin=imin;
	return vMin;
}
/***************************************************************************************/
cpedsPoint3D mscsFunction3dregc::getMinValueCell(bool Re) const {
	long imin,i,j,k;
	getMinValue(&imin,Re);
	cpedsPoint3D p;
	idx2ijk(imin,i,j,k);
	p.set(getX(i),getY(j),getZ(k));
	return p;
}
/***************************************************************************************/
cpedsPoint3D mscsFunction3dregc::getMaxValueCell(bool Re) const {
	long imin,i,j,k;
	getMaxValue(&imin,Re);
	cpedsPoint3D p;
	idx2ijk(imin,i,j,k);
	p.set(getX(i),getY(j),getZ(k));
	return p;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::populateField(cpedsPointSet3D &ps, bool ZasVals, double initialValueRe, double initialValueIm, bool redefineFunctionDomain, bool overflow, string how) {
	double x,y,z;	
	long N=ps.size();
	double minx,maxx,miny,maxy,minz,maxz;
	ps.getRanges(minx,maxx,miny,maxy,minz,maxz);
	//	printf("minx %lE maxx %lE miny %lE maxy %lE\n",minx,maxx,miny,maxy);
	
	int howLoc=-1;
	if (how=="mean") howLoc=0;
	if (how=="min") howLoc=1;
	if (how=="last") howLoc=2;
	if (how=="max") howLoc=3;
	if (howLoc==-1) return *this;
	
	mscsFunction3dregc counts;
	if (redefineFunctionDomain) {
		if (ZasVals) {
			setSizeRange(Nx(),Ny(),1,minx,miny,-1,maxx,maxy,1);
			counts.setSize(Nx(),Ny(),1);
		}
		else {
			setSizeRange(Nx(),Ny(),Nz(),minx,miny,minz,maxx,maxy,maxz);
			counts.setSize(Nx(),Ny(),Nz());
		}
		allocFunctionSpace();		
		setf(initialValueRe,initialValueIm);
		counts.allocFunctionSpace();
		counts=double(0);
	}
	else {
		allocFunctionSpace();
		setf(initialValueRe,initialValueIm);
		counts=(*this);
		counts=double(0);
	}
	
	
	// generate the field
	long ix,iy,iz;
	bool of, ovflw;
	
	if (ZasVals) {
		
		for (long k=0;k<N;k++) { // loop over points in point set
			//		ps[k].print_point("point");
			ix=coord2idx(ps[k].x(),0,&of,&x); 
			ovflw=of;
			//		msgs->say("ix: %li",ix,High);
			//		msgs->say("x: %lf",x,High);
			//		msgs->say("of: "+msgs->toStr(of),High);
			//		msgs->say("ovflw: "+msgs->toStr(ovflw),High);
			
			iy=coord2idx(ps[k].y(),1,&of,&y);  
			ovflw=ovflw | of;
			//		msgs->say("iy: %li",iy,High);
			//		msgs->say("y: %lf",y,High);
			//		msgs->say("of: "+msgs->toStr(of),High);
			//		msgs->say("ovflw: "+msgs->toStr(ovflw),High);
			
			
			
			iz=0;
			if (overflow) {
				if (ovflw==false) {
					//					msgs->say("ovflw==false",High);
					
					switch (howLoc) {
						case 0:
							fRe(ix,iy,iz)+=ps[k].z();
							break;
						case 1:
							if (ps[k].z() < fRe(ix,iy,iz) )	fRe(ix,iy,iz)=ps[k].z();
							break;
						case 2:
							fRe(ix,iy,iz)=ps[k].z();
							break;
						case 3:
							if (ps[k].z() > fRe(ix,iy,iz) )	fRe(ix,iy,iz)=ps[k].z();
							break;
						default:
							break;
					}
					counts.fRe(ix,iy,iz)++;
					fIm(ix,iy,iz)++;
				}
			}
			else {
				switch (howLoc) {
					case 0:
						fRe(ix,iy,iz)+=ps[k].z();
						break;
					case 1:
						if (ps[k].z() < fRe(ix,iy,iz) )	fRe(ix,iy,iz)=ps[k].z();
						break;
					case 2:
						fRe(ix,iy,iz)=ps[k].z();
						break;
					case 3:
						if (ps[k].z() > fRe(ix,iy,iz) )	fRe(ix,iy,iz)=ps[k].z();
						break;
					default:
						break;
				}
				
				counts.fRe(ix,iy,iz)++;
				fIm(ix,iy,iz)++;
			}
		}
	}
	else {
		
		
		
		for (long k=0;k<N;k++) { // loop over points in point set
			ix=coord2idx(ps[k].x(),0,&of,&x); 
			ovflw=of;
			iy=coord2idx(ps[k].y(),1,&of,&y);  
			ovflw=ovflw | of;
			
			iz=coord2idx(ps[k].z(),2,&of,&z); 
			ovflw=ovflw | of;
			if (overflow) {
				if (ovflw==false) {
					switch (howLoc) {
						case 0:
							fRe(ix,iy,iz)+=ps.val(k);
							break;
						case 1:
							if (ps.val(k) < fRe(ix,iy,iz) )	fRe(ix,iy,iz)=ps.val(k);
							break;
						case 2:
							fRe(ix,iy,iz)=ps.val(k);
							break;
						case 3:
							if (ps.val(k) > fRe(ix,iy,iz) )	fRe(ix,iy,iz)=ps.val(k);
							break;
						default:
							break;
					}
					counts.fRe(ix,iy,iz)++;
					fIm(ix,iy,iz)++;
				}
			}
			else {
				switch (howLoc) {
					case 0:
						fRe(ix,iy,iz)+=ps.val(k);
						break;
					case 1:
						if (ps.val(k) < fRe(ix,iy,iz) )	fRe(ix,iy,iz)=ps.val(k);
						break;
					case 2:
						fRe(ix,iy,iz)=ps.val(k);
						break;
					case 3:
						if (ps.val(k) > fRe(ix,iy,iz) )	fRe(ix,iy,iz)=ps.val(k);
						break;
					default:
						break;
				}
				counts.fRe(ix,iy,iz)++;
				fIm(ix,iy,iz)++;
			}
		}
	}
	//		printf("ix: %li iy: %li iz: %li\n\n\n",ix,iy,iz);
	
	
	//derive the mean of the field for the points that feel on the same grid cell
	if (ZasVals) {
		if (howLoc==0) {
			for (long i=0;i<Nx();i++) {
				for (long j=0;j<Ny();j++) {
					if (counts(i,j,0)[0]!=0) { fRe(i,j,0)/=counts.fRe(i,j,0);  }
				}
			}
		}
	}
	else {
		if (howLoc==0) {
			for (long i=0;i<Nx();i++) {
				for (long j=0;j<Ny();j++) {
					for (long k=0;k<Nz();k++) {
						if (counts.fRe(i,j,k)!=0) { fRe(i,j,k)/=counts.fRe(i,j,k);  }
					}
				}
			}
		}
	}
	
	return *this;
	
}
void mscsFunction3dregc::setf(double x, double y,double z, double re, double im) {
	long ix,iy,iz;
	//	bool overflw=false;
	
	ix=coord2idx(x,0);
	iy=coord2idx(y,1);
	iz=coord2idx(z,2);
	printf("ix: %li, iy: %li, iz: %li\n",ix,iy,iz);
	setf(ix,iy,iz,re,im);
}
/***************************************************************************************/
void mscsFunction3dregc::setRe(double value) {
	for (long i=0;i<Nx();i++) {
		for (long j=0;j<Ny();j++) {
			for (long k=0;k<Nz();k++) {
				fRe(i,j,k)=value;
			}
		}
	}	
}
void mscsFunction3dregc::setIm(double value) {
	for (long i=0;i<Nx();i++) {
		for (long j=0;j<Ny();j++) {
			for (long k=0;k<Nz();k++) {
				fIm(i,j,k)=value;
			}
		}
	}	
}
/***************************************************************************************/
double& mscsFunction3dregc::fReCoord(double x, double y,double z) {
	long ix,iy,iz;
	ix=coord2idx(x,0);
	iy=coord2idx(y,1);
	iz=coord2idx(z,2);
	return fRe(ix,iy,iz);
}
double& mscsFunction3dregc::fImCoord(double x, double y,double z) {
	long ix,iy,iz;
	ix=coord2idx(x,0);
	iy=coord2idx(y,1);
	iz=coord2idx(z,2);
	return fIm(ix,iy,iz);	
}
/***************************************************************************************/
double& mscsFunction3dregc::fReCoordPeriodic(double x, double y,double z) {
	long ix,iy,iz;
	ix=coord2idxPeriodic(x,0);
	iy=coord2idxPeriodic(y,1);
	iz=coord2idxPeriodic(z,2);
	return fRe(ix,iy,iz);
}
double& mscsFunction3dregc::fImCoordPeriodic(double x, double y,double z) {
	long ix,iy,iz;
	ix=coord2idxPeriodic(x,0);
	iy=coord2idxPeriodic(y,1);
	iz=coord2idxPeriodic(z,2);
	return fIm(ix,iy,iz);	
}

/***************************************************************************************/
long mscsFunction3dregc::coord2idx(double coord,long ax, bool *overflow, double *idxd) {
	double n, coordMin, totalLength;
	if (ax==0) { n=Nx(); coordMin=_param.x0; totalLength=_param.sizeX; }
	if (ax==1) { n=Ny(); coordMin=_param.y0; totalLength=_param.sizeY; }
	if (ax==2) { n=Nz(); coordMin=_param.z0; totalLength=_param.sizeZ; }
	double idx=(coord-coordMin) / (totalLength) * n;
	long idxI=long(idx);
	//	printf("idxI: %li\n",idxI);
	bool of=false;
	if (idxI<0) { idxI=0; of=true; }
	if (idxI>=n) { idxI=n-1; of=true; }
	if (overflow!=NULL) (*overflow)=of; 
	if (idxd!=NULL) (*idxd)=idx;
	return idxI;
}
/***************************************************************************************/
long mscsFunction3dregc::coord2idxPeriodic(double coord,long ax) {
	double n, coordMin, totalLength;
	if (ax==0) { n=Nx(); coordMin=_param.x0; totalLength=_param.sizeX; }
	if (ax==1) { n=Ny(); coordMin=_param.y0; totalLength=_param.sizeY; }
	if (ax==2) { n=Nz(); coordMin=_param.z0; totalLength=_param.sizeZ; }
	double idx=(coord-coordMin) / (totalLength) * n;
	long idxI=long(idx) % long(n);
	if (idxI<0) idxI+=n;
	return idxI;
}
/***************************************************************************************/
subDomain_region_t mscsFunction3dregc::getCellSDregion(long i, long j, long k) {
	subDomain_region_t sdr;
	sdr.xmin=getX(i)-getDxo2();
	sdr.xmax=getX(i+1)-getDxo2();
	sdr.ymin=getY(j)-getDyo2();
	sdr.ymax=getY(j+1)-getDyo2();
	sdr.zmin=getZ(k)-getDzo2();
	sdr.zmax=getZ(k+1)-getDzo2();
	return sdr;
}
/***************************************************************************************/
mscsFunction mscsFunction3dregc::averagePlane(long coord) {
	mscsFunction avf("av",Zero);
	avf=getSlice1Dfn(coord,0,0,0);
	double av;
	if (coord==0) {
		for (long i = 0; i < Nx(); i++) {
			av=0;
			for (long k = 0; k < Nz(); k++) {
				for (long j = 0; j < Ny(); j++) {
					av+=fRe(i,j,k);
				}
			}
			av/=_param.NYZ;
			avf.f(i)=av;
		}
	}
	if (coord==1) {
		for (long j = 0; j < Ny(); j++) {
			av=0;
			for (long i = 0; i < Nx(); i++) {
				for (long k = 0; k < Nz(); k++) {
					av+=fRe(i,j,k);
				}
			}
			av/=(Nx()*Nz());
			avf.f(j)=av;
		}
	}
	if (coord==2) {
		for (long k = 0; k < Nz(); k++) {
			av=0;
			for (long i = 0; i < Nx(); i++) {
				for (long j = 0; j < Ny(); j++) {
					av+=fRe(i,j,k);
				}
			}
			av/=_param.NXY;
			avf.f(k)=av;
		}
	}
	return avf;
}
/***************************************************************************************/
bool mscsFunction3dregc::intersects(mscsFunction3dregc& f) {
	if (
			getMaxX()>f.getMinX() and getMinX()<f.getMaxX() and
			getMaxY()>f.getMinY() and getMinY()<f.getMaxY() and
			getMaxZ()>f.getMinZ() and getMinZ()<f.getMaxZ()
	) return true;
	return false;
}
/***************************************************************************************/
void mscsFunction3dregc::intersection_set(mscsFunction3dregc& f, double value, int part) {
	if (intersects(f)) {
		long iSt,iEn,jSt,jEn,kSt,kEn;
		
		iSt=idxX(cpeds_get_max(getMinX(),f.getMinX())); //if (iSt==Nx() iSt)
		iEn=idxX(cpeds_get_min(getMaxX(),f.getMaxX())); if (iEn==Nx()) iEn--;
		
		jSt=idxY(cpeds_get_max(getMinY(),f.getMinY()));
		jEn=idxY(cpeds_get_min(getMaxY(),f.getMaxY())); if (jEn==Ny()) jEn--;
		
		kSt=idxZ(cpeds_get_max(getMinZ(),f.getMinZ()));
		kEn=idxZ(cpeds_get_min(getMaxZ(),f.getMaxZ())); if (kEn==Nz()) kEn--;
		
		//		printf("iSt: %li, iEn: %li\n",iSt,iEn);
		//		printf("jSt: %li, jEn: %li\n",jSt,jEn);
		//		printf("kSt: %li, kEn: %li\n",kSt,kEn);
		
		//		if (value==141) exit(0);
		//		if (iSt<0) exit(0);
		//		if (jSt<0) exit(0);
		//		if (kSt<0) exit(0);
		//		if (iEn>Nx()) exit(0);
		//		if (jEn>Ny()) exit(0);
		//		if (kEn>Nz()) exit(0);
		
		if (part==0) {
			for (long i = iSt; i <= iEn; i++) {
				for (long j = jSt; j <= jEn; j++) {
					for (long k = kSt; k <= kEn; k++) {
						fRe(i,j,k)=value;
					}
				}
			}			
		}
		
		if (part==1) {
			for (long i = iSt; i <= iEn; i++) {
				for (long j = jSt; j <= jEn; j++) {
					for (long k = kSt; k <= kEn; k++) {
						fIm(i,j,k)=value;						
					}
				}
			}			
		}
		
		
	}
}
/***************************************************************************************/
void mscsFunction3dregc::setInternalParam(string paramName, double paramVal) {
	if (paramName=="derivativeXperiodic") {		_param.derivativeXperiodic=paramVal;	}
	if (paramName=="derivativeYperiodic") {		_param.derivativeYperiodic=paramVal;	}
}
/***************************************************************************************/
mscsFunction3dregc mscsFunction3dregc::rotateSlice(int plane, int coord, double angle, double x0, double y0, double z0) {
	mscsFunction3dregc r(*this);
	if (angle==0 or angle==twoPI) return r;
	
	double v,x,y;
	long iSrc,jSrc;
	double xSrc,ySrc;
	bool ovflwX, ovflwY;
	if (plane==2) {
		for (long i = 0; i < Nx(); i++) {
			x=getX(i);
			for (long j = 0; j < Ny(); j++) {
				y=getY(j);
//				xSrc=(x+x0)*cos(angle) + (y+y0)*sin(angle) - x0;
//				ySrc=-(x+x0)*sin(angle) + (y+y0)*cos(angle) - y0;
				xSrc=(x-x0)*cos(angle) + (y-y0)*sin(angle) + x0;
				ySrc=-(x-x0)*sin(angle) + (y-y0)*cos(angle) + y0;
				iSrc=coord2idx(xSrc,0,&ovflwX);
				if (ovflwX) v=0;
				else {
					jSrc=coord2idx(ySrc,1,&ovflwY);
					if (ovflwY) v=0; 
					else v=fRe(iSrc,jSrc,coord);					
				}
				r.setf(i,j,0,v,0);
			}
		}
	}
	else {
		msgs->criticalError("mscsFunction3dregc::rotateSlice:: plane!=2 is not implemented yet",Top);
	}
	return r;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::centerMax(bool re) {
	long id;
	getMaxValue(&id,re);
	long i,j,k;
	long di,dj,dk;
	idx2ijk(id,i,j,k);
	di=Nx()/2-i;
	dj=Ny()/2-j;
	dk=Nz()/2-k;
	shift(di,dj,dk);
	return *this;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::insertFunctionXY(mscsFunction& XY, long i0, long j0, long k0) {
	assert(Ny()>=j0+XY.pointsCount());
	for (unsigned long j = 0; j < XY.pointsCount(); j++) {
		fRe(i0  ,j+j0,k0)=XY.getX(j);
		fRe(i0+1,j+j0,k0)=XY.f(j);
	}
	return *this;
}
/***************************************************************************************/
mscsFunction3dregc& mscsFunction3dregc::insertList(cpedsList<double>& l, long i0, long j0, long k0) {
	assert(Ny()>=j0+l.size());
	for (unsigned long j = 0; j < l.size(); j++) {
		fRe(i0,j+j0,k0)=l[j];
	}
	return *this;	

}
/* ******************************************************************************************** */
mscsFunction3dregc& mscsFunction3dregc::interpolateXYholes(int k) {
//	mscsFunction3dregc
	
/*
	for (long i = 0; i < Nx(); i++) {
		for (long j = 0; j < Ny(); j++) {
			if (fIm(i,j,k)!=1) {
				
			}
		}
	}
*/
	return *this;
}
