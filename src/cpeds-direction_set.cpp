#include "cpeds-direction_set.h"
#include "cpeds-math.h"
#include "cpeds-list.h"


// ****************************************************************************************************
cpedsDirectionSet::cpedsDirectionSet() {}
// ****************************************************************************************************
cpedsDirectionSet::cpedsDirectionSet(const QList<cpedsDirection> &q) {
  clear();
  long N=q.count();
  for (long i=0;i<N;i++) { append(q.at(i));  }
}
// ****************************************************************************************************
cpedsDirectionSet::cpedsDirectionSet(const cpedsDirectionSet &q) {
  clear();
  *this=q;
}
// ****************************************************************************************************
cpedsDirectionSet::cpedsDirectionSet(long N, double *lon, double *lat) {
  for (long i=0;i<N;i++) { append(cpedsDirection(lon[i],lat[i]));  }
}
// ****************************************************************************************************
cpedsDirectionSet::cpedsDirectionSet(long N, double *lon, double *lat, double *val) {
  for (long i=0;i<N;i++) { append(cpedsDirection(lon[i],lat[i],val[i]));  }
}
// ****************************************************************************************************
cpedsDirectionSet::cpedsDirectionSet(const QList<double>& lon, const QList<double>& lat, const QList<double>& val) {
  long N=lon.count();
  for (long i=0;i<N;i++) { append(cpedsDirection(lon[i],lat[i],val[i]));  }
}
// ****************************************************************************************************
cpedsDirectionSet::~cpedsDirectionSet() {

}

// ****************************************************************************************************
double* cpedsDirectionSet::getLonVals(long* size) const {
  *size=count();
  double *t = new double[*size];

  for (long i=0;i<*size;i++) { t[i]=at(i).lon();  }
  return t;
}
// ****************************************************************************************************
double* cpedsDirectionSet::getLatVals(long* size) const {
  *size=count();
  double *t = new double[*size];

  for (long i=0;i<*size;i++) { t[i]=at(i).lat();  }
  return t;
}
// ****************************************************************************************************
double* cpedsDirectionSet::getVals(long* size) const {
  *size=count();
  double *t = new double[*size];

  for (long i=0;i<*size;i++) { t[i]=at(i).val();  }
  return t;
}


// ****************************************************************************************************
void cpedsDirectionSet::getRanges(double &minlon, double &maxlon, double &minlat, double &maxlat) const {
  minlon=minLon();
  maxlon=maxLon();
  minlat=minLat();
  maxlat=maxLat();
}
// ****************************************************************************************************
double cpedsDirectionSet::extremal(bool lon, bool max)  const{
  long N;
  double* t;
  double v;
  if (lon) t=getLonVals(&N); 
  else t=getLatVals(&N);
  if (max) v=cpeds_find_max_value(t,N,0); 
  else v=cpeds_find_min_value(t,N,0); 
  delete t;
  return v;
}
// ****************************************************************************************************
double cpedsDirectionSet::minLon() const { return extremal(true, false); }
// ****************************************************************************************************
double cpedsDirectionSet::maxLon() const{ return extremal(true, true); }
// ****************************************************************************************************
double cpedsDirectionSet::maxLat() const { return extremal(false, true); }
// ****************************************************************************************************
double cpedsDirectionSet::minLat() const { return extremal(false, false); }

// ****************************************************************************************************
cpedsDirectionSet& cpedsDirectionSet::operator=(const cpedsDirectionSet &rhs) {
  if (this!=&rhs) { QList<cpedsDirection>::operator=(rhs); }
  // clear();
  // long N=rhs.count();
  // for (long i=0;i<N;i++) { append(rhs.at(i));  }
  return *this;
}
// ****************************************************************************************************
cpedsDirectionSet& cpedsDirectionSet::operator=(const QList<cpedsDirection> &rhs) {
  *this=cpedsDirectionSet(rhs);
  return *this;
}

// ****************************************************************************************************
void cpedsDirectionSet::print() const {
  long N=count();
  for (long i=0;i<N;i++) { at(i).print_direction();  } 
}
// ****************************************************************************************************
long cpedsDirectionSet::save(string fname, string how, bool withVals) const {
  return cpeds_matrix_save(exportAsMatrix(withVals), fname,how);
}

// ****************************************************************************************************
long cpedsDirectionSet::load(string fname, string how, bool withVals) {
  long res;
  importMatrix(cpeds_matrix_load(fname,how,&res),withVals);
  return res;
}

// ****************************************************************************************************
const matrix<double> cpedsDirectionSet::exportAsMatrix(bool withVals) const {
  long N=count();
  matrix<double> M;
  if (withVals) {
    M.SetSize(N,3);
    for (long i=0;i<N;i++) { M(i,0)=at(i).lon(); M(i,1)=at(i).lat(); M(i,2)=at(i).val(); }
  }
  else {
    M.SetSize(N,2);
    for (long i=0;i<N;i++) { M(i,0)=at(i).lon(); M(i,1)=at(i).lat(); }
  }
  return M;  
}

// ****************************************************************************************************
void cpedsDirectionSet::importMatrix(const matrix<double>& m, bool withVals) {
  long N=m.RowNo();
  if (m.ColNo()==1) {
    if (withVals) for (long i=0;i<N;i+=3) { append(cpedsDirection(m(i,0),m(i+1,0),m(i+2,0))); }
    else for (long i=0;i<N;i+=2) { append(cpedsDirection(m(i,0),m(i+1,0)));  }
  }
  if (m.ColNo()==2) {
    for (long i=0;i<N;i++) { append(cpedsDirection(m(i,0),m(i,1))); }
  }
  if (m.ColNo()==3) {
    for (long i=0;i<N;i++) { append(cpedsDirection(m(i,0),m(i,1),m(i,2))); }
  }
}

// ****************************************************************************************************
void cpedsDirectionSet::setLength(long n) {
   if (n<0) return;
   if (n==0) { clear(); return; }
   if (count() < n) { while (count() < n) { append(cpedsDirection()); }   }
   if (count() > n) { while (count() > n) { erase(0);  }   }
}
/***************************************************************************************/
cpedsDirection cpedsDirectionSet::mean() const {
	cpedsList<double> X,Y,Z;
	cpedsDirection av;
	double x,y;
	for (unsigned long i = 0; i < size(); i++) {
		X.append(cpeds_sph2cart(0,PIsnd-at(i).lat(),at(i).lon()));
		Y.append(cpeds_sph2cart(1,PIsnd-at(i).lat(),at(i).lon()));
		Z.append(cpeds_sph2cart(2,PIsnd-at(i).lat(),at(i).lon()));
	}
	av.setLat(PIsnd-cpeds_cart2sph(0,X.mean(),Y.mean(),Z.mean()));
	av.setLon(      cpeds_cart2sph(1,X.mean(),Y.mean(),Z.mean()));
	return av;
}
