#include "cpeds-point_set.h"
#include "cpeds-math.h"
#include "cpeds-rng.h"
#include "Mscs-function.h"

/****************************************************************************************************/
cpedsPointSet2D::cpedsPointSet2D() { }
/****************************************************************************************************/
cpedsPointSet2D::cpedsPointSet2D(const QList<cpedsPoint2D> &q) {
  clear();
  long N=q.count();
  for (long i=0;i<N;i++) { append(q.at(i));  }
}
/****************************************************************************************************/
cpedsPointSet2D::cpedsPointSet2D(const cpedsPointSet2D &q) {
  clear();
  *this=q;
}
/****************************************************************************************************/
cpedsPointSet2D::cpedsPointSet2D(long N, double *lon, double *lat) {
  for (long i=0;i<N;i++) { append(cpedsPoint2D(lon[i],lat[i]));  }
}
/***************************************************************************************/
cpedsPointSet2D::cpedsPointSet2D(const vector<cpedsPoint2D> &v) {  
	long N=v.size(); 
	for (long i=0;i<N;i++) { append(v[i]);  } 
}

/****************************************************************************************************/
cpedsPointSet2D::~cpedsPointSet2D() {

}

/****************************************************************************************************/
double* cpedsPointSet2D::getXvals(long* size) const {
  *size=count();
  double *t = new double[*size];

  for (long i=0;i<*size;i++) { t[i]=at(i).x();  }
  return t;
}
/****************************************************************************************************/
double* cpedsPointSet2D::getYvals(long* size) const {
  *size=count();
  double *t = new double[*size];

  for (long i=0;i<*size;i++) { t[i]=at(i).y();  }
  return t;
}


/****************************************************************************************************/
void cpedsPointSet2D::getRanges(double &minx, double &maxx, double &miny, double &maxy) const {
  minx=minX();
  maxx=maxX();
  miny=minY();
  maxy=maxY();
}
/****************************************************************************************************/
double cpedsPointSet2D::extremal(int coord, bool max) const {
  long N;
  double* t=NULL;
  double v;
  if (coord==0) { t=getXvals(&N); }
  else {
    if (coord==1) { t=getYvals(&N);
      //      else { t=getZvals(&N); }
    }
  }
  if (max) v=cpeds_find_max_value(t,N,0);
  else v=cpeds_find_min_value(t,N,0);
  delete t;
  return v;
}
/****************************************************************************************************/
double cpedsPointSet2D::minX() const { return extremal(0, false); }
/****************************************************************************************************/
double cpedsPointSet2D::maxX() const{ return extremal(0, true); }
/****************************************************************************************************/
double cpedsPointSet2D::maxY() const { return extremal(1, true); }
/****************************************************************************************************/
double cpedsPointSet2D::minY() const { return extremal(1, false); }

/****************************************************************************************************/
cpedsPointSet2D& cpedsPointSet2D::operator=(const cpedsPointSet2D &rhs) {
  if (this!=&rhs) {
//    this->QList<cpedsPoint2D>::operator=(rhs);
    this->mscsVector<cpedsPoint2D>::operator=(rhs);
  }
  return *this;
}
/****************************************************************************************************/
cpedsPointSet2D& cpedsPointSet2D::operator=(const QList<cpedsPoint2D> &rhs) {
  *this=cpedsPointSet2D(rhs);
  return *this;
}

/****************************************************************************************************/
void cpedsPointSet2D::print() const {
  long N=count();
  for (long i=0;i<N;i++) { at(i).print_point();  }
}
/****************************************************************************************************/
void cpedsPointSet2D::save(string fname) const {
  cpeds_matrix_save(exportAsMatrix(), fname);
}

/****************************************************************************************************/
const matrix<double> cpedsPointSet2D::exportAsMatrix() const {
  long N=count();
  matrix<double> M(N,2);
  for (long i=0;i<N;i++) { M(i,0)=at(i).x(); M(i,1)=at(i).y(); }
  return M;
}

/****************************************************************************************************/
/****************************************************************************************************/
/****************************************************************************************************/




















/****************************************************************************************************/
/****************************************************************************************************/
/****************************************************************************************************/
/****************************************************************************************************/
cpedsPointSet3D::cpedsPointSet3D() { }
/****************************************************************************************************/
cpedsPointSet3D::cpedsPointSet3D(const QList<cpedsPoint3D> &q) {
  clear();
  long N=q.count();
  for (long i=0;i<N;i++) { append(q.at(i));  }
}
/****************************************************************************************************/
cpedsPointSet3D::cpedsPointSet3D(const cpedsPointSet3D &q) {
  clear();
  *this=q;
}
/****************************************************************************************************/
cpedsPointSet3D::cpedsPointSet3D(const cpedsPointSet2D &q) {
  clear();
  *this=q;
}

/****************************************************************************************************/
cpedsPointSet3D::cpedsPointSet3D(long N, double *lon, double *lat, double* z) {
  for (long i=0;i<N;i++) { append(cpedsPoint3D(lon[i],lat[i],z[i]));  }
}
/***************************************************************************************/
cpedsPointSet3D::cpedsPointSet3D(long N, double *lon, double *lat, double z, bool deleteInside) {
	for (long i=0;i<N;i++) { append(cpedsPoint3D(lon[i],lat[i],z));  }	
	if (deleteInside) {
		delete [] lon;
		delete [] lat;
	}
}
/***************************************************************************************/
cpedsPointSet3D::cpedsPointSet3D(long N, double *lon, double *lat, double* z, double* v, bool deleteInside) {
	append(N,lon,lat,z,v,deleteInside);
}

/****************************************************************************************************/
cpedsPointSet3D::~cpedsPointSet3D() {
}
/***************************************************************************************/
cpedsPointSet3D& cpedsPointSet3D::append(long N, double *x, double *y, double* z, double* v, bool deleteInside) {
	if (z!=NULL) {
		for (long i=0;i<N;i++) { append(cpedsPoint3D(x[i],y[i],z[i]));  }	
		if (deleteInside) delete [] z;
	}
	else {
		for (long i=0;i<N;i++) { append(cpedsPoint3D(x[i],y[i],0));  }			
	}
	if (v!=NULL) {
		for (long i=0;i<N;i++) { _vals.push_back(v[i]); }		
		if (deleteInside) delete [] v;
	}
	
	if (deleteInside) {
		delete [] x;
		delete [] y;
	}
	return *this;
}
/****************************************************************************************************/
double* cpedsPointSet3D::getXvals(long* size) const {
  return getVals(0,size);
}

/****************************************************************************************************/
double* cpedsPointSet3D::getYvals(long* size) const {
  return getVals(1,size);
}
/****************************************************************************************************/
double* cpedsPointSet3D::getZvals(long* size) const {
  return getVals(2,size);
}
/****************************************************************************************************/
double* cpedsPointSet3D::getVals(int coord, long* size) const {
  *size=count();
  double *t = new double[*size];

  if (coord==0) {
    for (long i=0;i<*size;i++) { t[i]=at(i).x();  }
  }
  else {
    if (coord==1) {
      for (long i=0;i<*size;i++) { t[i]=at(i).y();  } }
    else {
      for (long i=0;i<*size;i++) { t[i]=at(i).z();  }
    }
  }

  return t;
}

/****************************************************************************************************/
void cpedsPointSet3D::getRanges(double &minx, double &maxx, double &miny, double &maxy, double &minz, double &maxz) const {
  minx=minX();
  maxx=maxX();
  miny=minY();
  maxy=maxY();
  minz=minZ();
  maxz=maxZ();
  // _Xmin=minx;
  // _Xmax=maxx;
  // _Ymin=miny;
  // _Ymax=maxy;
  // _Zmin=minz;
  // _Zmax=maxz;
}
/***************************************************************************************/
double cpedsPointSet3D::lengthX() const {
	return maxX()-minX();
}
/****************************************************************************************************/
double cpedsPointSet3D::lengthY() const {
	return maxY()-minY();
	
}
/****************************************************************************************************/
double cpedsPointSet3D::lengthZ() const {
	return maxZ()-minZ();
}
/***************************************************************************************/
void cpedsPointSet3D::shiftX(double delta) {
	for (unsigned long i = 0; i < size(); i++) {
		point(i).shiftX(delta);
	}
}
/****************************************************************************************************/
void cpedsPointSet3D::shiftY(double delta) {
	for (unsigned long i = 0; i < size(); i++) {
		point(i).shiftY(delta);
	}	
}
/****************************************************************************************************/
void cpedsPointSet3D::shiftZ(double delta) {
	for (unsigned long i = 0; i < size(); i++) {
		point(i).shiftZ(delta);
	}

}
/***************************************************************************************/
void cpedsPointSet3D::makePeriodic(double periodX,double periodY,double periodZ) {
//	cpedsPoint3D delta(0,0,0);
	if (periodX!=0) {
		for (unsigned long i = 0; i < size(); i++) {
			if (point(i).x()>periodX) point(i).x()-=long(point(i).x()/periodX)*periodX;
			if (point(i).x()<0) point(i).x()-=long(point(i).x()/periodX)*periodX;
		}		
	}
	if (periodY!=0) {
		for (unsigned long i = 0; i < size(); i++) {
			if (point(i).y()>periodY) point(i).y()-=long(point(i).y()/periodY)*periodY;
			if (point(i).y()<0) point(i).y()-=long(point(i).y()/periodY)*periodY;
		}		
	}
	if (periodZ!=0) {
		for (unsigned long i = 0; i < size(); i++) {
			if (point(i).z()>periodZ) point(i).z()-=long(point(i).z()/periodZ)*periodZ;
			if (point(i).z()<0) point(i).z()-=long(point(i).z()/periodZ)*periodZ;
		}		
	}
}
/****************************************************************************************************/
double cpedsPointSet3D::extremal(int coord, bool max) const {
  long N;
  double* t=NULL;
  double v;
  if (coord==0) { t=getXvals(&N); }
  else {
    if (coord==1) { t=getYvals(&N); }
    else { t=getZvals(&N); }
  }

  if (max) v=cpeds_find_max_value(t,N,0);
  else v=cpeds_find_min_value(t,N,0);
  delete [] t;
  return v;
}
/****************************************************************************************************/
double cpedsPointSet3D::maxX() const { return extremal(0, true); }
/****************************************************************************************************/
double cpedsPointSet3D::minX() const { return extremal(0, false); }
/****************************************************************************************************/
double cpedsPointSet3D::maxY() const { return extremal(1, true); }
/****************************************************************************************************/
double cpedsPointSet3D::minY() const { return extremal(1, false); }
/****************************************************************************************************/
double cpedsPointSet3D::maxZ() const { return extremal(2, true); }
/****************************************************************************************************/
double cpedsPointSet3D::minZ() const { return extremal(2, false); }
/****************************************************************************************************/
double cpedsPointSet3D::disp(int coord) const {
  long n;
  double *tmp;

  if (coord==0) { tmp=getXvals(&n); }
  else {
    if (coord==1) { tmp=getYvals(&n); }
    else {
      tmp=getZvals(&n); }
  }

  double rms = sqrt(cpeds_variance(tmp,n));
  delete tmp;
  return rms;
}
/****************************************************************************************************/
double cpedsPointSet3D::dispersionX() const {  return disp(0); }
/****************************************************************************************************/
double cpedsPointSet3D::dispersionY() const {  return disp(1); };
/****************************************************************************************************/
double cpedsPointSet3D::dispZ() const {  return disp(2); };
/***************************************************************************************/
void cpedsPointSet3D::multiply(double val) {
	long N=count();
	for (long i=0;i<N;i++) { at(i).multiply(val);  }
}
/***************************************************************************************/
void cpedsPointSet3D::add(double val) {
	long N=count();
	for (long i=0;i<N;i++) { at(i).add(val);  }
}

/****************************************************************************************************/
void cpedsPointSet3D::print() const {
  long N=count();
  for (long i=0;i<N;i++) { at(i).print_point();  }
}
/****************************************************************************************************/
void cpedsPointSet3D::save(string fname) const {
  cpeds_matrix_save(exportAsMatrix(), fname);
}
/***************************************************************************************/
void cpedsPointSet3D::load(string fname) {
	long cols=cpeds_get_cols_num_first_ln(fname.c_str(),true);
	msgs->say("detected "+msgs->toStr(cols)+" columns in the input file.",Low);
	FILE* f=fopen(fname.c_str(),"r");
	if (f==NULL) msgs->criticalError("cannot open file: "+fname,High);
	double x,y,z,v;
	string fmt;
	int res;
	char fmtch[20];
	if (cols==3) fmt="%lE %lE %lE";
	if (cols==4) fmt="%lE %lE %lE %lE";
	strcpy(fmtch,fmt.c_str());
	if (cols==3) {
		while (!feof(f)) {	res=fscanf(f,fmtch,&x,&y,&z); if (res!=EOF)	append(cpedsPoint3D(x,y,z));		/*printf("vec size: %li\n",size());*/ }
	}
	if (cols==4) {
		while (!feof(f)) {	res=fscanf(f,fmtch,&x,&y,&z,&v); if (res!=EOF) {		append(cpedsPoint3D(x,y,z));	_vals.push_back(v);	}	}
	}
	msgs->say("Loaded "+msgs->toStr(long(size()))+" points from the input file.",Medium);
	
	fclose(f);			
}
/****************************************************************************************************/
const cpedsPointSet3D& cpedsPointSet3D::operator=(const cpedsPointSet3D &rhs) {
  if (this!=&rhs) {
//    this->QList<cpedsPoint3D>::operator=(rhs);
    this->mscsVector<cpedsPoint3D>::operator=(rhs);
    values()=rhs.values();
  }
  return *this;
}
/****************************************************************************************************/
const cpedsPointSet3D& cpedsPointSet3D::operator=(const QList<cpedsPoint3D> &rhs) {
  *this=cpedsPointSet3D(rhs);
  return *this;
}
/****************************************************************************************************/
const cpedsPointSet3D& cpedsPointSet3D::operator=(const cpedsPointSet2D &rhs) {
  clear();
  long N=rhs.count();
  for (long i=0;i<N;i++) { append(cpedsPoint3D(rhs.at(i)));  }
  return *this;
}

/****************************************************************************************************/
const cpedsPointSet3D& cpedsPointSet3D::operator+=(const cpedsPointSet3D &rhs) {
  long N=cpeds_get_min(long(rhs.count()),long(count()));
  for (long i=0;i<N;i++) { value(i)+=rhs.at(i);  }
  return *this;
}
/****************************************************************************************************/
cpedsPointSet3D cpedsPointSet3D::operator+(const cpedsPointSet3D &rhs) {
  cpedsPointSet3D set(*this);
  set+=rhs;
  return set;
}

/****************************************************************************************************/
const matrix<double> cpedsPointSet3D::exportAsMatrix() const {
  long N=count();
  matrix<double> M;
  
  if (values().size()==size()) {
	  M.SetSize(N,4);
	  for (long i=0;i<N;i++) { M(i,0)=at(i).x(); M(i,1)=at(i).y(); M(i,2)=at(i).z(); M(i,3)=val(i); }
  }
  else {
	  M.SetSize(N,3);
	  for (long i=0;i<N;i++) { M(i,0)=at(i).x(); M(i,1)=at(i).y(); M(i,2)=at(i).z(); }
  }
  
  return M;
}
/****************************************************************************************************/
const matrix<double> cpedsPointSet3D::exportAsField(long sizeX, long sizeY) const {
  matrix<long> mask;
  return exportAsField(sizeX,sizeY,mask);
}

/****************************************************************************************************/
// this is depreciated
const matrix<double> cpedsPointSet3D::exportAsField(long sizeX, long sizeY, matrix<long>& mask) const {
  long N=count();
  matrix<double> M(sizeY,sizeX);
  matrix<double> counts(sizeY,sizeX);
  mask.SetSize(sizeY,sizeX);
  double x,y;
  mask-=1;

  double minx,maxx,miny,maxy,minz,maxz;
  getRanges(minx,maxx,miny,maxy,minz,maxz);
  printf("minx %lE maxx %lE miny %lE maxy %lE\n",minx,maxx,miny,maxy);
  // double resxo2=(maxx-minx)/sizeX/2;
  // double resyo2=(maxy-miny)/sizeY/2;
  // generate the field
  for (long k=0;k<N;k++) {
    if (maxx==minx) { x=0; }
    else {
      // x=(at(k).x()-minx-resxo2) / (maxx-minx) * sizeX;
      x=(at(k).x()-minx) / (maxx-minx) * sizeX - 0.5;
    }
    if (maxy==miny) { y=0; }
    else {
      // y=(at(k).y()-miny-resyo2) / (maxy-miny) * sizeY;
      y=(at(k).y()-miny) / (maxy-miny) * sizeY - 0.5;
    }

    // printf("sizex %li sizeY %li, x: %li  y: %li\n",sizeX, sizeY, long(x),long(y));
    M(long(y),long(x))+=at(k).z();
    mask(long(y),long(x))=k; // here the assumption is that the mask will store only the last value that falls into this cell
    counts(long(y),long(x))++;
  }

  //derive the mean of the field
  for (long i=0;i<sizeX;i++) {
    for (long j=0;j<sizeY;j++) {
      if (counts(j,i)!=0) { M(j,i)/=counts(j,i); }
    }
  }

  return M;
}
/***************************************************************************************/
//const mscsFunction3dregc cpedsPointSet3D::exportAsField3d(long sizeX, long sizeY, matrix<long>& mask) const {
//	  long N=count();
//	  mscsFunction3dregc M(sizeX,sizeY,1);
//	  mscsFunction3dregc counts(sizeX,sizeY,1);
//	  M.allocFunctionSpace();
//	  counts.allocFunctionSpace();
//	  M=double(0);
//	  counts=double(0);
//	  mask.SetSize(sizeY,sizeX);
//	  double x,y;
//	  mask-=1;
//
//	  double minx,maxx,miny,maxy,minz,maxz;
//	  getRanges(minx,maxx,miny,maxy,minz,maxz);
//	  printf("minx %lE maxx %lE miny %lE maxy %lE\n",minx,maxx,miny,maxy);
//	  // double resxo2=(maxx-minx)/sizeX/2;
//	  // double resyo2=(maxy-miny)/sizeY/2;
//	  // generate the field
//	  for (long k=0;k<N;k++) {
//	    if (maxx==minx) { x=0; }
//	    else {
//	      // x=(at(k).x()-minx-resxo2) / (maxx-minx) * sizeX;
//	      x=(at(k).x()-minx) / (maxx-minx) * sizeX - 0.5;
//	    }
//	    if (maxy==miny) { y=0; }
//	    else {
//	      // y=(at(k).y()-miny-resyo2) / (maxy-miny) * sizeY;
//	      y=(at(k).y()-miny) / (maxy-miny) * sizeY - 0.5;
//	    }
//
//	    // printf("sizex %li sizeY %li, x: %li  y: %li\n",sizeX, sizeY, long(x),long(y));
//	    M(long(x),long(y),0)[0]+=at(k).z();
//	    mask(long(y),long(x))=k; // here the assumption is that the mask will store only the last value that falls into this cell
//	    counts(long(x)long(y),0)[0]++;
//	  }
//
//	  //derive the mean of the field
//	  for (long i=0;i<sizeX;i++) {
//	    for (long j=0;j<sizeY;j++) {
//	      if (counts(i,j,0)[0]!=0) { M(i,j,0)[0]/=counts(i,j,0); }
//	    }
//	  }
//
//	  return M;
//	
//}

/****************************************************************************************************/
const cpedsPointSet3D cpedsPointSet3D::field2set(const matrix<double>& m, long refRow, long refCol, cpedsPoint3D p, cpedsPoint3D resolution) const {
  cpedsPointSet3D tmpps;
  long rows=m.RowNo();
  long cols=m.ColNo();
  double x,y,z;

  for (long i=0;i<rows;i++) {
    for (long j=0;j<cols;j++) {
      x=(j-refRow)*resolution.x()+p.x();
      y=(i-refCol)*resolution.y()+p.y();
      z=m(i,j);
      tmpps.append(cpedsPoint3D(x,y,z));
    }
  }
  return tmpps;
}
/****************************************************************************************************/
const cpedsPointSet3D& cpedsPointSet3D::setVals(int coord, double vX, double vY, double vZ ) {
  long N=count();

  for (long i=0;i<N;i++) {
    if (coord & 1) value(i).setX(vX);
    if (coord & 2) value(i).setY(vY);
    if (coord & 3) value(i).setZ(vZ);
  }
  return *this;
}
/***************************************************************************************/
const cpedsPointSet3D& cpedsPointSet3D::set(int i, double x, double y, double z, double v) {
	mscsVector<cpedsPoint3D>::at(i)=cpedsPoint3D(x,y,z);
	_vals[i]=v;
	return *this;
}
/***************************************************************************************/
const cpedsPointSet3D& cpedsPointSet3D::set(int i,const  cpedsPoint3D p, double v) {
	mscsVector<cpedsPoint3D>::at(i)=p;
	_vals[i]=v;
	return *this;	
}
/***************************************************************************************/
void cpedsPointSet3D::setVals(mscsVector<double> v) {
	_vals=v;	
}
/***************************************************************************************/
void cpedsPointSet3D::setVals(long N, double* v, bool deleteInside) {
	_vals.setSize(N);
	for (long i = 0; i < N; i++) {		_vals[i]=v[i];	}
	if (deleteInside) delete [] v;
}
/****************************************************************************************************/
const cpedsPointSet3D& cpedsPointSet3D::generatePoints(long n, double vX, double vY, double vZ) {

  for (long i=0;i<n;i++) { append(cpedsPoint3D(vX,vY,vZ));  }
  return *this;
}
/***************************************************************************************/
const cpedsPointSet3D& cpedsPointSet3D::generateRandomPoints(long nPts, subDomain_region_t sd) {
	cpedsRNG rns("uniform");
	double* x;
	double* y;
	double* z;
	setSize(nPts);
	rns.setMinMax(sd.xmin,sd.xmax);
	x=rns.getRN(nPts);
	rns.setMinMax(sd.ymin,sd.ymax);
	y=rns.getRN(nPts);
	rns.setMinMax(sd.zmin,sd.zmax);
	z=rns.getRN(nPts);
	for (long i = 0; i < nPts; i++) {
		set(i,x[i],y[i],z[i],1);
	}
	delete [] x;
	delete [] y;
	delete [] z;
	
	return *this;
}
/****************************************************************************************************/
const matrix<double> cpedsPointSet3D::exportAsInterpolatedField(long sizeX, long sizeY, matrix<long>& mask, cpedsPointSet3D& psint) const {
  matrix<double> m;
  printf("exporting as field\n");
  m=exportAsField(sizeX, sizeY, mask);
  // printf("saving mask\n");
  // cpeds_matrix_save(mask,"ilc5.mask.mat");
  mscsFunction f1("coord1",Zero);
  mscsFunction f2("coord2",Zero);

  //  mscsFunction f("coord");
  mscsFunction v("value",Zero);
  long pNum;
  long k;

  //printf("---------------------------- deriving interpolated coordinates ------------------------------\n");
  for (long j=0;j<sizeY;j++) {
    k=j;
    // printf("-------------- row number: %li\n",j);
    f1.clearFunction();
    f2.clearFunction();
    v.clearFunction();
    do {
    	for (long i=0;i<sizeX;i++) {
    		if (mask(k,i)>=0) {	// extract only original point coordinates from matrix
    			// printf("mask: row %li col %li %li\n",k,i,mask(k,i) );
    			f1.newPoint(double(i),at(mask(k,i)).x()); //! NOTE: here I only store the last id of the coordinate that falls in this matrix cell. Information about all previous points that also happen to fall in this cell is  not stored
    			f2.newPoint(double(i),at(mask(k,i)).y()); //! NOTE: here I only store the last id of the coordinate that falls in this matrix cell. Information about all previous points that also happen to fall in this cell is  not stored
    			v.newPoint(double(i),m(k,i));    																																						      //
    		}
    	}
    	// printf("k = %li\n",k);
    	if (f1.pointsCount() == 0) { if (k<sizeY) k++; else { printf("error on row %li in the map. This shouldn't happen, please check the code in functions exportAsInterpolatedField and exportAsField",k); return -1; } } // the requested resolution supersedes the map resolution -- will borrow some data from the next row
    } while (f1.pointsCount() == 0 and k<sizeY);

    // f1.save("coor1.rowlast");
    // f2.save("coor2.rowlast");
    //pNum=;
    // printf("number of original points before interpolation: %li %li\n", f1.pointsCount(), f2.pointsCount());
    // f.print();
    if (f1.pointsCount()>1) {       // interpolate
      f1.interpolate(1.0,true,"auto"); f2.interpolate(1.0,true,"auto");  v.interpolate(1.0,true,"linear");           }
    // f.print();
    // f.save("interp.txt");
    // printf("===extrapolation:\n");
    if (f1.pointsCount() < sizeX) {      // extrapolate
      f1.extrapolate(0,sizeX-1,1.0,true); f2.extrapolate(0,sizeX-1,1.0,true);       v.extrapolate(0,sizeX-1,1.0,true);     }
    // f.save("extrap.txt");
    // exit(0);
    // f1.save("coor1.rowlast.extra");
    // f2.save("coor2.rowlast.extra");
    // printf("====points count after extrapolation f1= %li f2=%li\n",f1.pointsCount(),f2.pointsCount());
    //f.print();
    if (f1.pointsCount()==sizeX and f2.pointsCount()==sizeX) {
      // printf("%li ******** points count f == sizeX; flagging  mask OK\n",j);
      // store inter/extra/polated x coordinates in row-major ordering to psint point set
      // store interpolated values onto the matrix m
      for (long i=0;i<sizeX;i++) {
	// if (mask(j,i)<0) mask(j,i)=-3;
	psint.append(cpedsPoint3D(f1.Y(i),f2.Y(i),v.Y(i)));
	m(j,i)=v.Y(i); } // -3 indicates that first coordinate was interpolated fine
    }
    else {
      printf("%li ******** points count f != sizeX; flagging mask NOK -- SOMETHING IS WRONG WITH THE CODE -- THIS MESSAGE SHOULD NOT APPEAR\n",j);
      f1.print();
      f2.print();
      exit(0);
      for (long i=0;i<sizeX;i++) { mask(j,i)=-2;  psint.append(cpedsPoint3D(0,0,0)); }
    } // -2 indicates that interpolated values were not calculated
  }



  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // //printf("---------------------------- deriving first coordinates ------------------------------\n");																														      //
  // // interpolate first coordinate																																							      //
  // for (long j=0;j<sizeY;j++) {																																							      //
  //   k=j;																																										      //
  //   // printf("-------------- row number: %li\n",j);																																					      //
  //   f.clearFunction();																																								      //
  //   v.clearFunction();																																								      //
  //   do {																																										      //
  //     for (long i=0;i<sizeX;i++) {																																							      //
  // 	if (mask(k,i)>=0) {	// extract only original point coordinates from matrix																																	      //
  // 	  // printf("mask: row %li col %li %li\n",k,i,mask(k,i) ); 																																			      //
  // 	  f.newPoint(double(i),at(mask(k,i)).x()); //! NOTE: here I only store the last id of the coordinate that falls in this matrix cell. Information about all previous points that also happen to fall in this cell is  not stored															      //
  // 	  v.newPoint(double(i),m(k,i));    																																						      //
  // 	}																																										      //
  //     }																																										      //
  //     // printf("k = %li\n",k);																																							      //
  //     if (f.pointsCount() == 0) { if (k<sizeY) k++; else { printf("error on row %li in the map. Cannot plot the map. This shouldn't happen, please check the code in functions exportAsInterpolatedField and exportAsField",k); return -1; } } // the requested resolution supersedes the map resolution -- will borrow some data from the next row	      //
  //   } while (f.pointsCount() == 0 and k<sizeY);																																					      //
  //   																																											      //
  //   // f.save("coor1.rowlast");																																							      //
  //   pNum=f.pointsCount();																																								      //
  //   // printf("number of original points before interpolation: %li\n", f.pointsCount());																																      //
  //   // f.print();																																									      //
  //   if (pNum>1) {       // interpolate																																						      //
  //     f.interpolate(1.0,true,"auto");  v.interpolate(1.0,true,"linear");           }																																	      //
  //   // f.print();																																									      //
  //   // f.save("interp.txt");																																								      //
  //   // printf("===extrapolation:\n");																																						      //
  //   if (pNum>=1 and pNum < sizeX) {      // extrapolate																																				      //
  //     f.extrapolate(0,sizeX-1,1.0,true);       v.extrapolate(0,sizeX-1,1.0,true);     }																																      //
  //   // f.save("extrap.txt");																																								      //
  //   // exit(0);																																									      //
  //   // f.save("coor1.rowlast.extra");																																						      //
  //   // printf("====points count after extrapolation f= %li\n",f.pointsCount());																																	      //
  //   //f.print();																																									      //
  //   if (f.pointsCount()==sizeX) {																																							      //
  //     // printf("%li ******** points count f == sizeX; flagging  mask OK\n",j);																																	      //
  //     // store inter/extra/polated x coordinates in row-major ordering to psint point set																																      //
  //     // store interpolated values onto the matrix m																																					      //
  //     for (long i=0;i<sizeX;i++) { if (mask(j,i)<0) mask(j,i)=-3;  psint.append(cpedsPoint3D(f.Y(i),0,v.Y(i)));       m(j,i)=v.Y(i); } // -3 indicates that first coordinate was interpolated fine																			      //
  //   }																																										      //
  //   else {																																										      //
  //     printf("%li ******** points count f != sizeX; flagging mask NOK -- SOMETHING IS WRONG WITH THE CODE -- THIS MESSAGE SHOULD NOT APPEAR\n",j);																									      //
  //     exit(0);																																									      //
  //     for (long i=0;i<sizeX;i++) { mask(j,i)=-2;  psint.append(cpedsPoint3D(0,0,0)); } 																																      //
  //   } // -2 indicates that interpolated values were not calculated																																			      //
  // }																																											      //
  // 																																											      //
  // // matrix<double> m1=m;																																								      //
  // // cpeds_matrix_save(m1,"interX.mat");																																						      //
  // 																																											      //
  // //printf("---------------------------- interpolating second coordinate ------------------------------\n");																														      //
  // // interpolate second coordinate																																							      //
  // for (long i=0;i<sizeX;i++) {																																							      //
  //   k=i;																																										      //
  //   f.clearFunction();																																								      //
  //   v.clearFunction();																																								      //
  //   do {																																										      //
  //     for (long j=0;j<sizeY;j++) {																																							      //
  // 	if (mask(j,i)>=0) {	 // extract only original point coordinates from matrix																																	      //
  // 	  f.newPoint(double(j),at(mask(j,i)).y());    																																					      //
  // 	  v.newPoint(double(j),m(j,i));    																																						      //
  // 	}																																										      //
  //     }																																										      //
  //     if (f.pointsCount() == 0) { if (k<sizeX) k++; else { printf("error on col %li in the map. Cannot plot the map. This shouldn't happen, please check the code in functions exportAsInterpolatedField and exportAsField",k); return -1; } } // the requested resolution supersedes the map resolution -- will borrow some data from the next row	      //
  //   } while (f.pointsCount() == 0 and k<sizeX);																																					      //
  // 																																											      //
  //   pNum=f.pointsCount();																																								      //
  //   // printf("x=%li number of points in the first column for interpolation: %li\n", i, f.pointsCount());																														      //
  //   if (pNum>1) {    																																								      //
  //     // interpolate																																									      //
  //     f.interpolate(1.0,true,"auto");  v.interpolate(1.0,true,"linear");           }																																	      //
  // 																																											      //
  //   if (pNum>=1 and pNum < sizeY) {      // extrapolate																																				      //
  //     // extrapolate																																									      //
  //     f.extrapolate(0,sizeY-1,1.0,true);       v.extrapolate(0,sizeY-1,1.0,true); }																																	      //
  // 																																											      //
  //   if (f.pointsCount()==sizeY) {																																							      //
  //     // printf("%li ******** points count f == sizeY; flagging  mask OK\n",i);																																	      //
  //     // store inter/extra/polated y coordinates in row-major ordering to psint point set																																      //
  //     // store interpolated values onto the matrix m																																					      //
  //     for (long j=0;j<sizeY;j++) { 																																							      //
  // 	if (mask(j,i)==-3) mask(j,i)=-4; // interpolation ok for second and ok for first 																																      //
  // 	if (mask(j,i)==-2) mask(j,i)=-5; // interpolation ok for second coord, nok for first																																      //
  // 	psint[j*sizeX+i].setY(f.Y(j));  																																						      //
  // 	if (mask(j,i)==-4) m(j,i)=(m(j,i)+v.Y(j))/2.0;																																					      //
  // 	if (mask(j,i)==-5) m(j,i)=v.Y(j);																																						      //
  // 	psint[j*sizeX+i].setZ(m(j,i));																																							      //
  //     }																																										      //
  //   }																																										      //
  //   else {																																										      //
  //     printf("%li ******** points count f != sizeY; flagging mask NOK -- SOMETHING IS WRONG WITH THE CODE -- THIS MESSAGE SHOULD NOT APPEAR\n",i);																									      //
  //     for (long j=0;j<sizeY;j++) { mask(j,i)=-2;  psint[j*sizeX+i].setY(0.0); } } // -2 indicates that interpolated values were not calculated																									      //
  // }																																											      //
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // // clean up the flags
  // for (long i=0;i<sizeX;i++) {
  //   for (long j=0;j<sizeY;j++) {
  //     if (mask(j,i)==-4;  psint[j*sizeX+i].setY(0.0); } } // -2 indicates that interpolated values were not calculated

  return m;
}
/***************************************************************************************/
cpedsPointSet3D::cpedsPointSet3D(const vector<cpedsPoint3D> &v) {
	for (long i = 0; i < v.size(); i++) {
		append(v.at(i));		
	}
}
/***************************************************************************************/
cpedsPointSet3D::cpedsPointSet3D(const mscsVector<cpedsPoint3D> &v) : mscsVector<cpedsPoint3D>(v) {
	
}
/***************************************************************************************/
cpedsPointSet3D::cpedsPointSet3D(const mscsVector<double>& x,const mscsVector<double>& y, const mscsVector<double>& z, const mscsVector<double>& v) {
	if ((x.size()==y.size()) and (x.size()==z.size()) and (x.size()==v.size())) {
		long N=x.size();
		for (long i = 0; i < N; i++) {
			append(cpedsPoint3D(x.at(i),y.at(i),z.at(i)));
			_vals.push_back(v.at(i));
		}
	}
}
/***************************************************************************************/
cpedsPointSet3D::cpedsPointSet3D(const mscsVector<double>& x,const mscsVector<double>& y, const mscsVector<double>& z) {
	if (x.size()==y.size() and x.size()==z.size()) {
		long N=x.size();
		for (long i = 0; i < N; i++) {		append(cpedsPoint3D(x.at(i),y.at(i),z.at(i)));		}
	}
}

/***************************************************************************************/
vector<cpedsPoint3D> cpedsPointSet3D::exportAsVector() const {
	vector<cpedsPoint3D> v;
	for (long i = 0; i < count(); i++) {		v.push_back(at(i));	}
	return v;
}

/***************************************************************************************/
void cpedsPointSet3D::printRanges(string comment) const {
	double xmin,xmax,ymin,ymax,zmin,zmax;
	getRanges(xmin,xmax,ymin,ymax,zmin,zmax);
	printf("%s> xmin: %lE, xmax: %lE, ymin: %lE, ymax: %lE, zmin: %lE, zmax: %lE\n",comment.c_str(),xmin,xmax,ymin,ymax,zmin,zmax);
}
/***************************************************************************************/
const cpedsPoint3D cpedsPointSet3D::massCenter() const {
	cpedsPoint3D c;
	double x=0,y=0,z=0;
	double M=0;
	for (unsigned long i = 0; i < size(); i++) {
		x+=val(i)*at(i).x();
		y+=val(i)*at(i).y();
		z+=val(i)*at(i).z();
		M+=val(i);
	}
	x/=M;
	y/=M;
	z/=M;
	c.set(x,y,z);
	return c;
}
/***************************************************************************************/
cpedsPointSet3D cpedsPointSet3D::calculateGravitationalPotential(double gravSoft) {
	cpedsPointSet3D gravPot=(*this);
	double V,d;
	cpedsPoint3D p,p2;
	
	for (unsigned long i = 0; i < size(); i++) {
		V=0;
		p=at(i);
		for (unsigned long j = 0; j < size(); j++) {
			if (i!=j) {
				p2=at(j);
				d=p.dist(p2);
				V-=val(j)/(d+gravSoft);
			}
		}
		V*=CPEDS_G;
		gravPot.values()[i]=V;
	}
	
	return gravPot;
}
/***************************************************************************************/
cpedsPoint3D cpedsPointSet3D::getMinValuePoint() {
	long iMin=0;
	double minv=val(iMin);
	
	for (unsigned long i = 0; i < size(); i++) {
		if (val(i)< minv) {			iMin=i; minv=val(i);		}
	}
	
	return point(iMin);
}
