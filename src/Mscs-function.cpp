#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <fstream>
#include <sstream>
#include <fftw3.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_errno.h>
#ifdef PGPLOT
#include <cpgplot.h>
#endif

#include "Mscs-function.h"


/* **************************************************************************************************** */
mscsFunction::mscsFunction() : mscsObject("function",Zero) {
	initiate_variables();
}

/* **************************************************************************************************** */
mscsFunction::mscsFunction(const string name, cpeds_VerbosityLevel verbosity) : mscsObject(name, verbosity) {
	initiate_variables();
}
/* **************************************************************************************************** */
mscsFunction::mscsFunction(const string name, cpeds_point* t, long num) : mscsObject(name) {
	initiate_variables();
	if (t!=NULL)  importFunction(t,num);
	else {
		t=new cpeds_point[num];
		for (int i = 0; i < num; ++i) {		t[i].x=0; t[i].y=0;	}
		importFunction(t,num);
		delete [] t; t=NULL;
	}
}
/* **************************************************************************************************** */
mscsFunction::mscsFunction(const string name, double *x, double *y, long num, cpeds_VerbosityLevel verbosity) : mscsObject(name, verbosity) {
	initiate_variables();
	importFunction(x,y,num);
}
/* **************************************************************************************************** */
mscsFunction::mscsFunction(const string name, const QList<double>& x, const QList<double>& y) : mscsObject(name) {
	initiate_variables();
	long num=x.count();
	for (long i=0;i<num;i++) { newPoint(x.value(i),y.value(i));  }
}
/* **************************************************************************************************** */
mscsFunction::mscsFunction(const mscsFunction& fn ) : mscsObject(fn) {
	// clearFunction();
	// long N=fn.pointsCount();
	// for (long i=0;i<N;i++) {
	//   newPoint(fn.getPoint(i));
	// }
	initiate_variables();
	*this=fn;
}
/* **************************************************************************************************** */
mscsFunction::~mscsFunction() {
	if (_rns!=NULL) delete _rns; _rns=NULL;
}

/* **************************************************************************************************** */
void mscsFunction::initiate_variables() {
	_rns=NULL;
	clearRanges();
	
}



/* **************************************************************************************************** */
mscsFunction::frange mscsFunction::getRanges() const { return range; }
/* **************************************************************************************************** */
double mscsFunction::getMinValue() const { return range.ymin; }
/* **************************************************************************************************** */
double mscsFunction::getMaxValue() const { return range.ymax; }
/* **************************************************************************************************** */
double mscsFunction::getMinArg() const { return range.xmin; }
/* **************************************************************************************************** */
double mscsFunction::getMaxArg() const { return range.xmax; }
/* ******************************************************************************************** */
std::pair<double,double> mscsFunction::getMinMaxValues() const {
	return std::make_pair(getMinValue(),getMaxValue());
}

/* **************************************************************************************************** */
mscsFunction& mscsFunction::checkRanges() {
	double* tmp;
	clearRanges();
	
	tmp=extractArguments();
	if (tmp!=NULL) {
		cpeds_find_minmax_value(tmp,_f.count(),&range.xmin,&range.xmax,&range.ixmin,&range.ixmax); 
		delete [] tmp;
	}
	else { range.xmin=range.xmax=0.0; range.ixmin=range.ixmax=0; }
	tmp=extractValues();
	if (tmp!=NULL) {
		cpeds_find_minmax_value(tmp,_f.count(),&range.ymin,&range.ymax,&range.iymin,&range.iymax);
		delete [] tmp;
	}
	else { range.ymin=range.ymax=0.0; range.iymin=range.iymax=0; }
	
	return *this;
}

void mscsFunction::clearRanges() {
	range.xmin=range.ymin=range.xmax=range.ymax=0;
	range.ixmax=range.ixmin=range.iymax=range.iymin=0;
}
/* **************************************************************************************************** */
double* mscsFunction::extractArguments(long n) const {
	if (n==0) n=pointsCount();
	if (n>0) {
		long i;
		double * tmp = new double[n];
		
		for (i=0;i<n;i++) {
			tmp[i]=_f[i].x();
		}
		return tmp;
	}
	return NULL;
}
/***************************************************************************************/
double* mscsFunction::extractArguments(long imin,long imax) const {
	long n=imax-imin+1;
	if (n>0) {
		long i;
		double * tmp = new double[n];
		
		for (i=imin;i<=imax;i++) {
			tmp[i-imin]=_f[i].x();
		}
		return tmp;
	}
	return NULL;	
}
/***************************************************************************************/
double* mscsFunction::extractValues(long imin,long imax) const {
	long n=imax-imin+1;
	if (n>0) {
		long i;
		double * tmp = new double[n];
		
		for (i=imin;i<=imax;i++) {
			tmp[i-imin]=_f[i].y();
		}
		return tmp;
	}
	return NULL;	
}
/***************************************************************************************/
cpedsList<double> mscsFunction::toXList() const {
	cpedsList<double> tmp;
	if (_f.count()>0) {
		long i;
		
		for (i=0;i<_f.count();i++) {
			tmp.append(_f[i].x());
		}
	}
	return tmp;
}
/* **************************************************************************************************** */
const QList<double> mscsFunction::getArguments() const {
	QList<double> l;
	long num=_f.count();
	for (long i=0;i<num;i++) { l.append(_f[i].x());  }
	return l;
}
/* **************************************************************************************************** */
const QList<double> mscsFunction::getValues() const {
	QList<double> l;
	long num=_f.count();
	for (long i=0;i<num;i++) { l.append(_f[i].y());  }
	return l;
}

/* **************************************************************************************************** */
double* mscsFunction::extractValues(long n) const {
	if (n==0) n=pointsCount();
	if (n>0) {
		long i;
		double * tmp = new double[n];
		for (i=0;i<n;i++) {
			tmp[i]=_f[i].y();
		}
		return tmp;
	}
	return NULL;
}
/***************************************************************************************/
cpedsList<double> mscsFunction::toYList() const {
	cpedsList<double> tmp;
	if (_f.count()>0) {
		long i;
		
		for (i=0;i<_f.count();i++) {
			tmp.append(_f[i].y());
		}
	}
	return tmp;	
}
/* **************************************************************************************************** */
void mscsFunction::setf(long i, double x, double y) {
	if (argOK(i)) {
		_f[i].setX(x);  _f[i].setY(y);
		if (x<range.xmin) { range.xmin=x; range.ixmin=i; }
		if (x>range.xmax) { range.xmax=x; range.ixmax=i; }
		if (y<range.ymin) { range.ymin=y; range.iymin=i; }
		if (y>range.ymax) { range.ymax=y; range.iymax=i; }
	}
}
/* **************************************************************************************************** */
void mscsFunction::setf(long i,  double y)  {
	if (argOK(i)) {  _f[i].setY(y);
	if (y<range.ymin) { range.ymin=y; range.iymin=i; }
	if (y>range.ymax) { range.ymax=y; range.iymax=i; }
	}
}
/* **************************************************************************************************** */
void mscsFunction::setarg(long i,  double x) {
	if (argOK(i)) {  _f[i].setX(x);
	if (x<range.xmin) { range.xmin=x; range.ixmin=i; }
	if (x>range.xmax) { range.xmax=x; range.ixmax=i; }
	}
}
/* **************************************************************************************************** */
void mscsFunction::setPointsNum(long num) {
	long num_cur=pointsCount();
	if (num>num_cur) {
		for (long i=num_cur;i<num;i++) { newPoint(0,0); }
	}
	else {
		if (num<num_cur) {
			for (long i=num;i<num_cur;i++) { deletePoint(num); }
		}
	}
}
/***************************************************************************************/
double mscsFunction::Y(long i) const { 
	return getY(i);
}
/* **************************************************************************************************** */
double mscsFunction::getY(long i) const {
	if (argOK(i)) return _f[i].y();
	else {
//		msgs->error("index out of range: "+msgs->toStr(i),Medium);
		msgs->say("index out of range: "+msgs->toStr(i),Medium);
	}
	return 0;
}
/* **************************************************************************************************** */
double mscsFunction::f(double x, long* i) const {
	if (_f.count()>0) {
#ifdef DEBUG
		if (getx(0)>getx(1) and pointsCount()>=2) {
			msgs->warning("mscsFunction::f(double x, long* i) const >> The function seems not to be sorted. THIS OPERATION WILL FAIL!!! ",Top);
		}
#endif
		double* tmp = extractArguments();
		long ii=cpeds_find_value(x,tmp,_f.count(),0,_f.count());
		delete [] tmp;
		if (i!=NULL) *i=ii;
		return _f[ii].y();
	}
	return 0;
}
/***************************************************************************************/
double mscsFunction::finter(double x, string type) const {
	//double mscsFunction::finter(double x, string type) {
	long idx1;
	double a[5],v[5];
	long n=pointsCount();
	
	v[2]=f(x,&idx1);
#ifdef DEBUG
	printf("idx1: %li\n",idx1);
	//	checkRanges();
	setVerbosityLevel(Top);
	printRanges();
#endif
	if (n<2) {
		msgs->say("cannot interpolate for < 2 points. closest point will be returned",High);
		return f(x,NULL);
	}
	
	if (n<5) {
		if (type!="linear") { msgs->say(type+" interpolation is not possible with < 5 points. Changing type to 'linear'",High); type="linear"; }
		//		printf("x: %lf, idx1: %li\n", x,idx1);
		if (x<getx(idx1)) {
			if (idx1>0) {
				a[0]=getx(idx1-1); a[1]=getx(idx1);
				v[0]=y(idx1-1); v[1]=y(idx1); 			
			}
			else { // extrapolate low
				a[0]=getx(idx1); a[1]=getx(idx1+1);
				v[0]=y(idx1); v[1]=y(idx1+1); 			
			}
		}
		else {
			if (idx1<n-1) {
				a[0]=getx(idx1); a[1]=getx(idx1+1);
				v[0]=y(idx1); v[1]=y(idx1+1); 
			}
			else { // extrapolate high
				a[0]=getx(idx1-1); a[1]=getx(idx1);
				v[0]=y(idx1-1); v[1]=y(idx1);
			}
		}
		//		printf("a0: %lf, a1: %lf, v0: %lf, v1: %lf\n",a[0],a[1],v[0],v[1]);
		return cpeds_extrapolate_linear(x,a[0],v[0],a[1],v[1]);
		
	}
	else {
		if (idx1>=2 and idx1<=n-3) {
			a[0]=getx(idx1-2); a[1]=getx(idx1-1); a[2]=getx(idx1); a[3]=getx(idx1+1); a[4]=getx(idx1+2); 
			v[0]=Y(idx1-2); v[1]=Y(idx1-1); v[2]=Y(idx1); v[3]=Y(idx1+1); v[4]=Y(idx1+2); 		
			//		printf("betw: %lf %lf %lf %lf %lf\n",a[0],a[1],a[2],a[3],a[4]);
		}
		else {
			if (idx1==1) {
				a[0]=getx(idx1-1); a[1]=getx(idx1); a[2]=getx(idx1+1); a[3]=getx(idx1+2); a[4]=getx(idx1+3); 
				v[0]=Y(idx1-1); v[1]=Y(idx1); v[2]=Y(idx1+1); v[3]=Y(idx1+2); v[4]=Y(idx1+3); 		
			}
			if (idx1==0) {
				a[0]=getx(idx1); a[1]=getx(idx1+1); a[2]=getx(idx1+2); a[3]=getx(idx1+3); a[4]=getx(idx1+4); 
				v[0]=Y(idx1); v[1]=Y(idx1+1); v[2]=Y(idx1+2); v[3]=Y(idx1+3); v[4]=Y(idx1+4);
			}
			if (idx1==n-1) {
				a[0]=getx(idx1-4); a[1]=getx(idx1-3); a[2]=getx(idx1-2); a[3]=getx(idx1-1); a[4]=getx(idx1); 
				v[0]=Y(idx1-4); v[1]=Y(idx1-3); v[2]=Y(idx1-2); v[3]=Y(idx1-1); v[4]=Y(idx1);
				//			printf("n-1: %lf %lf %lf %lf %lf\n",a[0],a[1],a[2],a[3],a[4]);
			}
			if (idx1==n-2) {
				a[0]=getx(idx1-3); a[1]=getx(idx1-2); a[2]=getx(idx1-1); a[3]=getx(idx1); a[4]=getx(idx1+1); 
				v[0]=Y(idx1-3); v[1]=Y(idx1-2); v[2]=Y(idx1-1); v[3]=Y(idx1); v[4]=Y(idx1+1);
				//			printf("n-2: %lf %lf %lf %lf %lf\n",a[0],a[1],a[2],a[3],a[4]);
			}
		}		
		
		
		if (x>a[4] or x<a[0]) {
#ifdef DEBUG
			printf("linearly extrapolating at: %lE\n",x);
#endif
			if (x>a[4]) return cpeds_extrapolate_linear(x,a[3],v[3],a[4],v[4]);
			else return cpeds_extrapolate_linear(x,a[0],v[0],a[1],v[1]);
		}
		else {
			double* vint;
			vint=cpeds_interpolate(a,v,5,&x,1,type);
#ifdef DEBUG
			cpeds_print_matrix(a,5,1);
			cpeds_print_matrix(v,5,1);
			printf("vint: %lE\n",vint[0]);
#endif
			v[0]=vint[0]; delete [] vint;
		}
	}
	
	return v[0];
}
/* **************************************************************************************************** */
//double mscsFunction::f(double x) const {
//	long i;
//	return f(x,&i);
//}
/* **************************************************************************************************** */
void mscsFunction::newPoint(double x, double y) {
	_f.append(QPointF(x,y));
	if (_f.count()==1) { range.xmin=x; range.ixmin=_f.count()-1; range.xmax=x; range.ixmax=_f.count()-1; range.ymin=y; range.iymin=_f.count()-1; range.ymax=y; range.iymax=_f.count()-1; }
	if (x<range.xmin) { range.xmin=x; range.ixmin=_f.count()-1; }
	if (x>range.xmax) { range.xmax=x; range.ixmax=_f.count()-1; }
	if (y<range.ymin) { range.ymin=y; range.iymin=_f.count()-1; }
	if (y>range.ymax) { range.ymax=y; range.iymax=_f.count()-1; }
}
/* **************************************************************************************************** */
void mscsFunction::newPoint(const QPointF& p) {
	newPoint(p.x(),p.y());
}
/* **************************************************************************************************** */
void mscsFunction::newPoint(cpeds_point p) { 
	printf("WARNING: I am in mscsFunction::newPoint(cpeds_point p) make sure that this routine actually works");
	newPoint(p.x,p.y); }
//   _f.append(QPointF(p.x,p.y));
//  if (p.x<range.xmin) { range.xmin=p.x; range.ixmin=_f.count()-1; }
//  if (p.x>range.xmax) { range.xmax=p.x; range.ixmax=_f.count()-1; }
//  if (p.y<range.ymin) { range.ymin=p.y; range.iymin=_f.count()-1; }
//  if (p.y>range.ymax) { range.ymax=p.y; range.iymax=_f.count()-1; }
// }
/* **************************************************************************************************** */
void mscsFunction::insertPoint(double x, double y,long i,long N) {
	for (long j = 0; j < N; j++) {
		_f.insert(i,QPointF(x,y));
	}
	if (x<range.xmin) { range.xmin=x; range.ixmin=i; }
	if (x>range.xmax) { range.xmax=x; range.ixmax=i; }
	if (y<range.ymin) { range.ymin=y; range.iymin=i; }
	if (y>range.ymax) { range.ymax=y; range.iymax=i; }
}
/* **************************************************************************************************** */
void mscsFunction::insertPoint(cpeds_point p,long i) {
	_f.insert(i,QPointF(p.x,p.y));
	if (p.x<range.xmin) { range.xmin=p.x; range.ixmin=i; }
	if (p.x>range.xmax) { range.xmax=p.x; range.ixmax=i; }
	if (p.y<range.ymin) { range.ymin=p.y; range.iymin=i; }
	if (p.y>range.ymax) { range.ymax=p.y; range.iymax=i; }
}

/* **************************************************************************************************** */
void mscsFunction::deletePoint(long i) { _f.removeAt(i); }
/* **************************************************************************************************** */
void mscsFunction::deletePoint() { if ( pointsCount()>0) _f.removeLast(); }
/* **************************************************************************************************** */
void mscsFunction::clearFunction() { 
	//	qDeleteAll(_f); 
	_f.clear(); 
	clearRanges(); 
}

/***************************************************************************************/
matrix<double> mscsFunction::valuesToVector(bool vertical) {
	matrix<double> m;
	long n=pointsCount();
	if (vertical) {
		m.SetSize(n,1);
		for (long i = 0; i < n; i++) {		m(i,0)=f(i);	}
	}
	else {
		m.SetSize(1,n);
		for (long i = 0; i < n; i++) {		m(0,i)=f(i);	}		
	}
	return m;
}

/* **************************************************************************************************** */
bool mscsFunction::argOK(long i) const {
	if (i>=0 && i<_f.count()) {
		return true;
	}
	return false;
}
/***************************************************************************************/
double& mscsFunction::fpsty(long i) { 
	i=i % pointsCount();
	if (i<0) i=(pointsCount()+i);
	return _f[i].ry(); 
}

/* **************************************************************************************************** */
void mscsFunction::importFunction(double* X, double* Y, long size, bool deleteAfter) {
	msgs->say("Importing new data points: "+msgs->toStr(size),Low);
	if (X!=NULL and Y!=NULL) {
		for (long i=0;i<size;i++) {
			newPoint(X[i],Y[i]);
			// _f.append(QPointF(X[i],Y[i]));
		}
	}
	else {
		if (X==NULL) {			for (long i=0;i<size;i++) {	newPoint(i,Y[i]);	}					}
		if (Y==NULL) {			for (long i=0;i<size;i++) {	newPoint(X[i],0);	}					}
	}
	
	if (deleteAfter) {
		if (X!=NULL) delete [] X;
		if (Y!=NULL) delete [] Y;
	}
}
/***************************************************************************************/
void mscsFunction::importFunction(mscsVector<double> X, mscsVector<double> Y) {
	for (long i=0;i<X.size();i++) {
		newPoint(X[i],Y[i]);
	}	
}
/* **************************************************************************************************** */
void mscsFunction::importFunction(cpeds_point* p, long size) {
	msgs->say("Importing new data points: "+msgs->toStr(size),Low);
	for (long i=0;i<size;i++) {
		newPoint(p[i]);
	}
}


/* **************************************************************************************************** */
const mscsFunction& mscsFunction::setX(const QList<double> x) {
	//	clearFunction();
	//	long N=x.size();
	//	for (long i = 0; i < N; i++) {
	//		newPoint(x.at(i),0);
	//	}
	long N=x.size(); /*
	 * Comment: changed implementation here. 
	 * This change may possibly influence work of some of the programs,
	 * but I believe this implementation makes more sense.
	 * 
	 * If this is not what you want then use setPointsNum first.
	 * 
	 * author: blew
	 * date: Jun 25, 2012 7:19:47 PM
	 *
	 */
	
	for (long i = 0; i < N; i++) {
		setarg(i,x.at(i));
	}
	return *this;
}
/* **************************************************************************************************** */
const mscsFunction& mscsFunction::setX(double v) {
	long N=pointsCount();
	for (long i = 0; i < N; i++) {		setarg(i,v);	}
	return *this;	
}

/* **************************************************************************************************** */
const mscsFunction& mscsFunction::setY(const QList<double> y) {
	//	long N=pointsCount();
	//	for (long i = 0; i < N; i++) {		setf(i,y.at(i));	}
	//	return *this;
	
	long N=y.size();		/*
	 * Comment: changed implementation here. 
	 * This change may possibly influence work of some of the programs,
	 * but I believe this implementation makes more sense.
	 * 
	 * author: blew
	 * date: Jun 25, 2012 7:19:47 PM
	 *
	 */
	
	for (long i = 0; i < N; i++) {		setf(i,y.at(i));	}
	return *this;
}

/* **************************************************************************************************** */
double mscsFunction::getX(long i) const { if (argOK(i)) return _f[i].x(); return 0; }

/* **************************************************************************************************** */
long mscsFunction::pointsCount() const { return _f.count(); }

/***************************************************************************************/
mscsFunction mscsFunction::get(long from, long to) {
	mscsFunction f;
	if (argOK(to)==false) to=pointsCount();
	for (long i = from; i < to; i++) {
		f.newPoint(get(i));
	}
	return f;
}
/* **************************************************************************************************** */
cpedsStatusCodes mscsFunction::load(string filename, bool commentedFile, int colx, int coly, long long Nrows, long long skipRows) {
	FILE* f;
	double x,y;
	int res;
	bool invertFn=false;
	clearFunction();
	msgs->say("Loading function from file "+filename,High);
	long cols=cpeds_get_cols_num_first_ln(filename.c_str(),commentedFile,10000);
	if (cols==-1) { msgs->error("ERROR: cannot load from the file: "+filename+". No such file", High);	  return cpedsNoSuchFile; }
	
	msgs->say("Number of columns in the input file: "+msgs->toStr(cols),Medium);
	bool xonly, yonly;
	if (colx==-1) yonly=true; else yonly=false;
	if (coly==-1) xonly=true; else xonly=false;
	if (xonly && yonly) { msgs->error("ERROR: mscsFunction::load- wrong parameters", High);	  return wrongFormat;	  }
	if (cols<1) { msgs->error("ERROR: cannot load from the file: "+filename+". Too few columns", High);	  return wrongFormat;	  }
	if (colx>coly) { 
		cpeds_swap(colx,coly); 
		//		printf("%i %i\n",colx,coly); 
		invertFn=true; }
	
	//
	// create the right format to read in the function
	//
	string format;
	string skipFormat="  %*s";
	if (cols==1) {
		format="%lE";	
	}
	else {
		//		format="%lE %lE";	
		//		if (colx==0 and coly==1) {
		//			for (int i = 2; i < cols; ++i) {	format+=" %*E";  }
		//		}
		//		else {
		format="";
		for (int i = 0; i < colx; ++i) {	format+=skipFormat;  }
		if (colx>-1) format+=" %lE";
		for (int i = colx+1; i < coly; ++i) {	format+=skipFormat;  }
		format+=" %lE";
		for (int i = coly+1; i < cols; ++i) {	format+=skipFormat;  }
#ifdef DEBUG
		printf("%s\n",format.c_str());
#endif
		//		exit(0);
		//		}		
	}
	
	if (commentedFile) {
		if (Nrows!=-1 or skipRows!=0) { msgs->warning("Nrows!=-1 and  skipRows!=0 are NOT IMPLEMENTED YET FOR COMMENTED FILES. Will ignore these parameters.",Top); }
		long j=0;
		ifstream ifs(filename.c_str());
		if (ifs.is_open()) {
			string s;
			if (xonly) {
				msgs->say("Reading into X; enumerating y: "+msgs->toStr(cols),Medium);
				while (getline(ifs,s)) {
					if (s[0]!= '#') { sscanf(s.c_str(),format.c_str(),&x); newPoint(x,double(j)); }											
					j++;
				}
			}
			else {
				if (yonly) {
					msgs->say("Reading into y; enumerating x: "+msgs->toStr(cols),Medium);
					while (getline(ifs,s)) {
						if (s[0]!= '#') { sscanf(s.c_str(),format.c_str(),&y); newPoint(double(j),y); }												
						j++;
					}
				}
				else {
					while (getline(ifs,s)) {
						if (s[0]!= '#') { sscanf(s.c_str(),format.c_str(),&x,&y); newPoint(x,y); }						
					}
				}
			}
		}
		else { msgs->error("ERROR: cannot load from the file: "+filename+". No such file", High);	  return cpedsNoSuchFile;	  }
	}
	else {
		long j=0;
		f= fopen(filename.c_str(),"r");
		if (f!=NULL) {
			if (xonly) {
				if (Nrows==-1) {
					while (!feof(f)) { res=fscanf(f,format.c_str(),&x); if (res!=EOF and j>=skipRows) newPoint(x,double(j)); j++; }
				}
				else {
					long long i=Nrows+skipRows;;
					while (i>0 and !feof(f)) { res=fscanf(f,format.c_str(),&x); if (res!=EOF and j>=skipRows) newPoint(x,double(j)); j++; i--; }
				}
			}
			else {
				if (yonly) {
					if (Nrows==-1) {
						while (!feof(f)) { res=fscanf(f,format.c_str(),&y); if (res!=EOF and j>=skipRows) newPoint(double(j),y); j++; }
					}
					else {
						long long i=Nrows+skipRows;;
						while (i>0 and !feof(f)) { res=fscanf(f,format.c_str(),&y); if (res!=EOF and j>=skipRows) newPoint(double(j),y); j++; i--; }
					}
				}
				else {
					long long i=Nrows+skipRows;;
					if (Nrows==-1) {
						while (!feof(f)) { res=fscanf(f,format.c_str(),&x,&y); if (res!=EOF and j>=skipRows) newPoint(x,y); }					
					}
					else {
						while (i>0 and !feof(f)) { res=fscanf(f,format.c_str(),&x,&y); if (res!=EOF and j>=skipRows) newPoint(x,y); j++; i--; }
					}
				}
			}
			fclose(f);
			msgs->say("Read "+msgs->toStr((long)_f.count())+" lines",Medium);
			
			if (invertFn) invert();
			return cpedsSuccess;
		}
		else { msgs->error("ERROR: cannot load from the file: "+filename+". No such file", High);	  return cpedsNoSuchFile;	  }
	}
	
	if (invertFn) invert();
	return cpedsSuccess;
}
/* **************************************************************************************************** */
cpedsStatusCodes mscsFunction::save(string filename, bool append) const {
	FILE* F;
	//double x,y;
	
	msgs->say("Saving function to file "+filename,High);
	if (append)
		F= fopen(filename.c_str(),"a");
	else
		F= fopen(filename.c_str(),"w");
	if (F!=NULL) {
		for (long i=0;i<_f.count();i++) {
			fprintf(F,"%.20lE %.20lE\n",getX(i),getY(i));
		}
		fclose(F);
		stringstream ss; ss<<"Saved "<<pointsCount()<<" points.";
		msgs->say(ss.str(),High);
		
		return cpedsSuccess;
	}
	else { msgs->error("ERROR: cannot save to file", High); }
	
	return cpedsCannotWriteToFile;
}
/* **************************************************************************************************** */
//void mscsFunction::print() const {
void mscsFunction::print(bool printPointNo) const {
	msgs->say("printing the function. Points number: "+msgs->toStr(pointsCount()),High);
	printf("\n");
	if (printPointNo) {
		for (long i=0;i<_f.count();i++) {
			printf("%li| %.15lE %.15lE\n",i,getX(i),getY(i));
			//		printf("%li| %.15lE %.15lE\n",i,X(i),f(i));
		}
	}
	else {
		for (long i=0;i<_f.count();i++) {
			printf("%.15lE %.15lE\n",getX(i),getY(i));
		}
	}
	printf("\n");
}
/***************************************************************************************/
double mscsFunction::sumX() const {
	double s=0;
	for (long i=0;i<_f.count();i++) { s+=getx(i); }
	return s;
}
/***************************************************************************************/
double mscsFunction::sumY() const {
	double s=0;
	for (long i=0;i<_f.count();i++) { s+=y(i); }
	return s;	
}

/* **************************************************************************************************** */
mscsFunction& mscsFunction::add(double v) {
	for (long i=0;i<_f.count();i++) { _f[i].setY(_f[i].y()+v); }
	return *this;
}
/* **************************************************************************************************** */
mscsFunction& mscsFunction::subtract(double v) {
	for (long i=0;i<_f.count();i++) { _f[i].setY(_f[i].y()-v); }
	return *this;
}
/* **************************************************************************************************** */
mscsFunction& mscsFunction::multiply(double v) {
	for (long i=0;i<_f.count();i++) { _f[i].setY(_f[i].y()*v); }
	return *this;
}
/* **************************************************************************************************** */
mscsFunction& mscsFunction::divide(double v) {
	if (v!=0.0) {
		for (long i=0;i<_f.count();i++) { _f[i].setY(_f[i].y()/v); }
	}
	else { msgs->error(" cannot divide by 0\n",High); }
	return *this;
}
/* **************************************************************************************************** */
mscsFunction& mscsFunction::power(double a) {
	for (long i=0;i<_f.count();i++) { _f[i].setY(pow(_f[i].y(),a)); }
	return *this;
}
/***************************************************************************************/
mscsFunction& mscsFunction::absoluteValue() { //<! calculates absolute value of the function
	for (long i=0;i<_f.count();i++) { _f[i].setY(fabs(_f[i].y())); }
	checkRanges();
	return *this;	
}
/* **************************************************************************************************** */
mscsFunction& mscsFunction::sqroot() {
	for (long i=0;i<_f.count();i++) { if (_f[i].y()>0) _f[i].setY(sqrt(_f[i].y())); else _f[i].setY(0); }
	checkRanges();
	return *this;
}
/* **************************************************************************************************** */
mscsFunction& mscsFunction::exponent() {
	for (long i=0;i<_f.count();i++) { _f[i].setY(exp(double(_f[i].y()))); }
	checkRanges();
	return *this;
}
/***************************************************************************************/
mscsFunction& mscsFunction::exponentX(double base) {
	for (long i=0;i<_f.count();i++) { _f[i].setX(pow(base,double(_f[i].x()))); }
	checkRanges();
	return *this;	
}
/***************************************************************************************/
mscsFunction& mscsFunction::exponentY(double base) {
	for (long i=0;i<_f.count();i++) { _f[i].setY(pow(base,double(_f[i].y()))); }	
	checkRanges();
	return *this;	
}

/***************************************************************************************/
mscsFunction& mscsFunction::logX(double base) {
	long n=_f.count();
	double logbase=log(base);
	for (long i=0;i<n;i++) { _f[i].setX(log(double(_f[i].x()))/logbase); }
	checkRanges();
	return *this;	
}
/***************************************************************************************/
mscsFunction& mscsFunction::logY(double base) {
	long n=_f.count();
	double logbase=log(base);
	for (long i=0;i<n;i++) { _f[i].setY(log(double(_f[i].y()))/logbase); }
	checkRanges();
	return *this;		
}
/***************************************************************************************/
mscsFunction& mscsFunction::lnY() { //<! convenience function  that converts Y values into a ln Y values
	long n=_f.count();
	for (long i=0;i<n;i++) { _f[i].setY(log(double(_f[i].y()))); }
	checkRanges();
	return *this;			
}


/* **************************************************************************************************** */
void mscsFunction::Xadd(long i, double v) {  if (argOK(i)) { _f[i].setX(_f[i].x()+v); } }
void mscsFunction::xadd(long i, double v) {   _f[i].setX(_f[i].x()+v);  }
/* **************************************************************************************************** */
void mscsFunction::Xadd(double* v) {  
	long n=pointsCount();
	for (long i=0;i<n;i++) { _f[i].rx()+=v[i]; }
}
/* **************************************************************************************************** */
void mscsFunction::Xsubtract(long i, double v) { if (argOK(i)) { _f[i].setX(_f[i].x()-v); } }
/* **************************************************************************************************** */
void mscsFunction::Xmultiply(long i, double v) { if (argOK(i)) { _f[i].setX(_f[i].x()*v); } };
/* **************************************************************************************************** */
void mscsFunction::Xdivide(long i, double v) { if (argOK(i)) { _f[i].setX(_f[i].x()/v); } }
/* **************************************************************************************************** */
void mscsFunction::Yadd(long i, double v) { if (argOK(i)) { _f[i].setY(_f[i].y()+v); } else printf("argument not valid\n");}
void mscsFunction::yadd(long i, double v) { _f[i].setY(_f[i].y()+v); }
//void mscsFunction::Yadd(long i, double v) {  _f[i].setf(_f[i].y()+v); }
/* **************************************************************************************************** */
void mscsFunction::Ysubtract(long i, double v) { if (argOK(i)) { _f[i].setY(_f[i].y()-v); } }
void mscsFunction::ysubtract(long i, double v) { _f[i].setY(_f[i].y()-v); }
/* **************************************************************************************************** */
void mscsFunction::Ymultiply(long i, double v) { if (argOK(i)) { _f[i].setY(_f[i].y()*v); } }
void mscsFunction::ymultiply(long i, double v) {  _f[i].setY(_f[i].y()*v); }
/* **************************************************************************************************** */
void mscsFunction::Ydivide(long i, double v) { if (argOK(i)) { _f[i].setY(_f[i].y()/v); } }
void mscsFunction::ydivide(long i, double v) {  _f[i].setY(_f[i].y()/v); } 
/* **************************************************************************************************** */
void mscsFunction::asVectorMultiply(long i, double v) { Xmultiply(i,v); Ymultiply(i,v); }
/* **************************************************************************************************** */
double mscsFunction::modulus(long i) { double x,y; x=getX(i); y=getY(i); return sqrt(x*x+y*y);  }

/* **************************************************************************************************** */
//! adds a function
mscsFunction& mscsFunction::add(const mscsFunction& F, long startFrom) {
	long N=F.pointsCount();
	for (long i=0;i<N;i++) {  Yadd(i+startFrom,F.y(i));   }
	return *this;
}
/* **************************************************************************************************** */
//! adds a function
mscsFunction& mscsFunction::add(const double*t) {
	long N=pointsCount();
	for (long i=0;i<N;i++) {  Yadd(i,t[i]);   }
	return *this;
}
/***************************************************************************************/
double mscsFunction::vectorSum() {
	double x=sumX();
	double y=sumY();
	return sqrt(x*x+y*y);
}
/***************************************************************************************/
mscsFunction& mscsFunction::vectorAdd(const mscsFunction& F) {
	long N=F.pointsCount();
	for (long i=0;i<N;i++) {  Xadd(i,F.getx(i)); Yadd(i,F.y(i));   }
	return *this;	
}
/***************************************************************************************/
mscsFunction& mscsFunction::vectorAdd(const double* v) {
	long N=pointsCount();
	for (long i=0;i<N;i++) {  Xadd(i,v[i]); Yadd(i,v[i]);   }
	return *this;	
}
/***************************************************************************************/
mscsFunction& mscsFunction::vectorSubtract(const mscsFunction& F) {
	long N=F.pointsCount();
	for (long i=0;i<N;i++) {  Xsubtract(i,F.getx(i)); Ysubtract(i,F.y(i));   }
	return *this;	
}
/***************************************************************************************/
mscsFunction& mscsFunction::scalarMultiply(const mscsFunction& F) {
	long N=F.pointsCount();
	for (long i=0;i<N;i++) {  Xmultiply(i,F.getx(i)); Ymultiply(i,F.y(i));   }
	return *this;	
}

/* **************************************************************************************************** */
//! subtracts a function
mscsFunction& mscsFunction::subtract(const mscsFunction& F) {
	for (long i=0;i<F.pointsCount();i++) {    Ysubtract(i,F.y(i));  }
	return *this;
}
/* **************************************************************************************************** */
//! multiplies by a function
mscsFunction& mscsFunction::multiply(const mscsFunction& F) {
	long N=F.pointsCount();
	for (long i=0;i<N;i++) {    Ymultiply(i,F.y(i));  }
	return *this;
}
/* **************************************************************************************************** */
//! divides by a function
mscsFunction& mscsFunction::divide(const mscsFunction& F) {
	long N=F.pointsCount();
	for (long i=0;i<N;i++) {    Ydivide(i,F.y(i));  }
	return *this;
}
/***************************************************************************************/
mscsFunction& mscsFunction::divide(const cpedsList<double>& l) {
	long N=cpeds_get_min(pointsCount(),l.size());
	for (long i=0;i<N;i++) {    Ydivide(i,l[i]);  }
	return *this;	
}
/***************************************************************************************/
mscsFunction& mscsFunction::divide(const cpedsList<long>& l) {
	long N=cpeds_get_min(pointsCount(),l.size());
	for (long i=0;i<N;i++) {    Ydivide(i,double(l[i]));  }
	return *this;	
}
/* **************************************************************************************************** */
mscsFunction& mscsFunction::epoweri(double phi, bool phiR) {
	long N=_f.count();
	if (phiR) {
		for (long i=0;i<N;i++) { X(i)+=phi;  }				
	}
	else {
		double cosphi=cos(phi);
		double sinphi=sin(phi);
		for (long i=0;i<N;i++) {    X(i)=X(i)*cosphi-f(i)*sinphi; f(i)=X(i)*sinphi+f(i)*cosphi;  }		
	}
	return *this;	
}
/* **************************************************************************************************** */
mscsFunction& mscsFunction::epowerii(long i, double phi, bool phiR) {
	if (phiR) { X(i)+=phi;	}
	else {
		double cosphi=cos(phi);
		double sinphi=sin(phi);
		X(i)=X(i)*cosphi-f(i)*sinphi; f(i)=X(i)*sinphi+f(i)*cosphi;
	}
	return *this;	
}

/* **************************************************************************************************** */
//! shifts the function
mscsFunction& mscsFunction::shiftX(double x) {
	for (long i=0;i<_f.count();i++) {    Xadd(i,x);  }
	return *this;
}

/* **************************************************************************************************** */
mscsFunction& mscsFunction::foldXinto(double from, double to, double acc) {
	long N=_f.count();
	double x,d=to-from;
	long Nfold;
	//	if (acc>0) to-=acc;
	for (long i=0;i<N;i++) {
		x= X(i);
		//		printf("%li >  x: %.15lE, foldX: x-to: %.10lE\n",i,x*PI180inv, (x-to)*PI180inv);
		if (x<to and to-x< acc and acc>0) x=to; // snap to the upper limit to fold points lying very close to the upper limit
		
		if (x>=to) {
			x-=from;
			Nfold=int(trunc(x/d));
			x-=Nfold*d;
			x+=from;
			setarg(i,x);
			//			printf("folding %li >  Nfold: %li, x: %.10lE\n",i,Nfold, x*PI180inv);
		} else
			if (x<from) {
				x-=from;
				Nfold=int(trunc(x/d))-1;
				x-=Nfold*d;
				x+=from;
				setarg(i,x);
			}
	}
	return *this;
}

/* **************************************************************************************************** */
//! scales the function
mscsFunction& mscsFunction::scaleX(double a) {
	for (long i=0;i<_f.count();i++) {    Xmultiply(i,a);  }
	return *this;
}
/* **************************************************************************************************** */

// double& mscsFunction::operator () (long x, long y) { return &f(x)
// double mscsFunction::operator () (long x, long y);


/* **************************************************************************************************** */
#ifdef PGPLOT
int mscsFunction::plotPG(int plottype) {
	double* x=extractArguments();
	double* y=extractValues();
	checkRanges();
	
	cpgbeg(0,"/XWINDOW",1,1);
	cpgask(false);
	cpgenv(getMinValue(), getMaxValue(), getMinArg(), getMaxArg(), 0, 2); //coordinates
	cpglab("x", "y", "title"); //axes
	cpgsci(9);
	cpgbin((int)pointsCount(),x,y,false);
	do {
		cpgcurs(&x,&y,&ch);
	} while ((ch != ' ') && (ch != 'p'));
	delete [] hist; delete [] xbin;
	cpgclos();
	if (ch == ' ') retval = 0;
	if (ch == 'p') retval = -2;
	
	delete [] x;
	delete [] y;
	return retval;
}
#endif


/* **************************************************************************************************** */
const mscsFunction& mscsFunction::operator=(const mscsFunction& rhs) {
#ifdef DEBUG
	msgs->say("entering assignment operator. The random number generator will not be cloned",Low);
#endif
	
	//	clearRanges();
	//	initiate_variables(); // basically initiate_variables();  is a one time call function and should not be called multiple times; so be careful with calling it.
	if (this != &rhs) {
		//		clearFunction();
		this->mscsObject::operator=(rhs);
		this->_f=rhs._f;
		range=rhs.range;
		//		this->mscsObject::operator=(rhs); // I do not want to copy the messanger nor the random number generator
	}
	
	//     clearFunction();
	//     long N=rhs.pointsCount();
	//     for (long i=0;i<N;i++) {
	// #ifdef DEBUG
	//       // msgs->say("assigning point no: "+msgs->toStr(i),Low);
	//       printf("assigning point no: %li\n",i);
	// #endif
	//       newPoint(rhs.getPoint(i));
	//     }
	//   }
	return *this;
}
/***************************************************************************************/
bool mscsFunction::operator==(const mscsFunction& rhs) {
	if (pointsCount()==rhs.pointsCount()) {
		for (long i = 0; i < pointsCount(); i++) {
			if (getx(i)!=rhs.getx(i) or f(i)!=rhs.f(i)) return false;
		}
		return true;
	}
	return false;
}
/***************************************************************************************/
double mscsFunction::meanf() const {
	double *X=extractValues();
	double m=cpeds_mean_value(X,pointsCount());
	delete [] X;
	return m;
}
/***************************************************************************************/
double mscsFunction::meanX(bool weightByY) const {
	double m=0;
	if (weightByY) {
		for (unsigned long i = 0; i < pointsCount(); i++) {			m+=getx(i)*y(i);		}
		m/=sumY();
	}
	else {
		double *X;
		X=extractArguments();
		m=cpeds_mean_value(X,pointsCount());
		delete [] X;		
	}
	return m;
	
}
/***************************************************************************************/
double mscsFunction::stdev() {
	return sqrt(variance());
}
/***************************************************************************************/
double mscsFunction::variance() {
	double *X=extractValues();
	double v=cpeds_variance(X,pointsCount());
	delete [] X;
	return v;
}
/***************************************************************************************/
double mscsFunction::rms() {
	double *X=extractValues();
	double v=cpeds_rms(X,pointsCount());
	delete [] X;
	return v;
}
/* **************************************************************************************************** */
const mscsFunction mscsFunction::derivative(bool inPlace, double* periodX, double* periodY) {
	long N=pointsCount();
	
	if (N>1) {
		long i;
		//		matrix<double> m(N,1);
		double* xvals = extractArguments();
		double* yvals = extractValues();
		//		for (i=0;i<N;i++) { m(i,0)=f(i); }
		//		cpeds_matrix_derivative(&m,xvals);
#ifdef DEBUG_MSCSFN_DERIVATIVES
		//		cpeds_matrix_save(m,"deriv_mat.ref");
		cpeds_derivative(xvals,yvals,pointsCount(),NULL,NULL);		
		cpeds_save_matrix(yvals,pointsCount(),1,"deriv_old",false);
		delete [] yvals;
		yvals = extractValues();
		cpeds_derivative(xvals,yvals,pointsCount(),periodX,periodY);		
		cpeds_save_matrix(yvals,pointsCount(),1,"deriv",false);
#endif
		cpeds_derivative(xvals,yvals,pointsCount(),periodX,periodY);		
		
		if (inPlace) {
			for (i=0;i<N;i++) { setf(i,yvals[i]); }
			delete [] xvals;			delete [] yvals;
			return mscsFunction(*this);
		}
		else {
			mscsFunction fp("derivative");
			for (i=0;i<N;i++) { fp.newPoint(xvals[i],yvals[i]); }
			delete [] xvals;			delete [] yvals;
			return fp;
		}
	}
	return mscsFunction("undefined derivative");
	
}
/* **************************************************************************************************** */
const mscsFunction& mscsFunction::sortFunctionArgAscending() {
	QPointF* t= extractPoints();
	long num=_f.count();
	cmpPt2 fp=&cmp_points_X;
	qsort(t,num,sizeof(QPointF),(cmpPt)fp);
	importPoints(t,num);
	delete [] t;
	return *this;
}
/* **************************************************************************************************** */
const mscsFunction& mscsFunction::sortFunctionArgDescending() {
	QPointF* t= extractPoints();
	long num=_f.count();
	cmpPt2 fp=&cmp_points_X;
	qsort(t,num,sizeof(QPointF),(cmpPt)fp);
	importPoints(t,num,true);
	delete [] t;
	return *this;
}
/* **************************************************************************************************** */
int mscsFunction::cmp_points_X(const QPointF &p1, const QPointF &p2) {
	if (p1.x() < p2.x()) return -1;
	if (p1.x() > p2.x()) return 1;
	return 0;
}
/* **************************************************************************************************** */
int mscsFunction::cmp_points_Y(const QPointF &p1, const QPointF &p2) {
	if (p1.y() < p2.y()) return -1;
	if (p1.y() > p2.y()) return 1;
	return 0;
}

/* **************************************************************************************************** */
QPointF* mscsFunction::extractPoints() const {
	long num=_f.count();
	QPointF* t= new QPointF[num];
	
	for (long i=0;i<num;i++) {
		t[i]=_f[i];
	}
	return t;
}

/* **************************************************************************************************** */
void mscsFunction::importPoints(QPointF* t, long n, bool reverseOrder) {
	clearFunction();
	if (reverseOrder) {
		for (long i=1;i<=n;i++) {
			newPoint(t[n-i].x(),t[n-i].y());
		}
	}
	else {
		for (long i=0;i<n;i++) {
			newPoint(t[i].x(),t[i].y());
		}
	}
}
/* **************************************************************************************************** */
const mscsFunction mscsFunction::interpolate(double *Xint, long Nint, string type, bool inPlace) {
	double *x=extractArguments();
	double *y=extractValues();
	double *Yint=cpeds_interpolate(x,y,_f.count(),Xint,Nint,type);
	
	// for (long i=0;i<Nint;i++) {    insertPoint(Xint[i],Yint[i]);  }
	mscsFunction f("Interpolating function",Xint,Yint,Nint,Zero);
	//  f.sortFunctionArgAscending();
	
	delete [] x;
	delete [] y;
	delete [] Yint;
	if (inPlace) {
		*this=f;
		//		clearFunction();
		//		concatenate(f);
		//		sortFunctionArgAscending(); // BLcomment (Jul 7, 2011, 12:25:18 PM): is this thing really needed here ? I think not, since inPlace does not return the combined result.
	}
	return f;
}
/***************************************************************************************/
const mscsFunction mscsFunction::interpolate(double dx, string type, double DX, string type2) {
	mscsFunction tmpf("interpolatedFunction",getVerbosityLevel());
	mscsFunction tmpf2("interpolatedFunction2",getVerbosityLevel());
	mscsFunction tmpf3("interpolatedFunction3",getVerbosityLevel());
	
	tmpf=(*this);
	tmpf.interpolate(dx,true,type);
	double f,t;
	f=0;
	t=0;
	for (long i = 0; i < pointsCount()-1; i++) {
		if (getX(i+1)-getX(i) > DX) { 
			msgs->say("found region of increased points spacing",High);
			f=getX(i); t=getX(i+1);
			tmpf.cut(f,t);
			tmpf2=copy(f,t);
			tmpf2.interpolate(dx,true,type2);
			tmpf3.concatenate(tmpf2);
		}
	}
	
	tmpf.concatenate(tmpf3);
	tmpf.sortFunctionArgAscending();
	return tmpf;
}
/***************************************************************************************/
mscsFunction mscsFunction::cut(double from, double to) {
	long fromi, toi;
	f(from,&fromi);
	f(to,&toi);
	if (getX(fromi)< from) fromi++;
	if (getX(toi)> to) toi--;
	mscsFunction tmpf("tmpFunction",getVerbosityLevel());
	if (fromi>toi) return tmpf;
	return cut(toi-fromi+1,fromi);
}
/***************************************************************************************/
mscsFunction mscsFunction::copy(double from, double to) {
	long fromi, toi;
	f(from,&fromi);
	f(to,&toi);
	if (getX(fromi)< from) fromi++;
	if (getX(toi)> to) toi--;
	mscsFunction tmpf("tmpFunction",getVerbosityLevel());
	if (fromi>toi) return tmpf;
	return copy(toi-fromi+1,fromi);
}
/* **************************************************************************************************** */
//const mscsFunction mscsFunction::interpolate(double *onX, long Nx, bool inPlace, string type, bool exact, double acc) {
//	long iSt,iEn;
//	double x,st,en;
//	mscsFunction tmpfn("interpolating arguments",Zero);
//
//	//
//	// define ranges and interpolation grid
//	//
//	f(onX[0],&iSt);
//	f(onX[Nx-1],&iEn);
////			printf("iSt %li iEn %li\n", iSt, iEn);
//	if (exact) { // set from and to to the first and the last respectively valid values inside of the region where function is defined.
//		if (iSt==0) { st=X(iSt); } else { st=from; }
//
//		x=onX[0];
//		while (x<st) { if (abs(st-x)<acc) x=st; else x+=dx; }
//		from=x;
//
//		x=from;
//		en=cpeds_get_min(X(iEn),to);
////				printf("from %lE to %lE en %lE\n", from, to,en);
//		double en2=en+acc;
//		while (x<en2) { if (abs(x-en)<acc) { x=en; tmpfn.newPoint(x,0.0); x=en2; printf("acc\n"); } else { tmpfn.newPoint(x,0.0); x+=dx; } }
//		//if (abs(x-en)<acc) { x=en; tmpfn.newPoint(x,0.0); }
//		if (acc==0) x-=dx;
//		if (x==en2) x=en;
//		to=x;
////				printf("from %lE to %lE\n", from, to);
//	}
//	else {
//		from=getX(iSt);
//		to=getX(iEn);
//		// set wider range to include derivatives information for the interpolation
//		if (iSt>=1) { tmpfn.insertPoint(getX(iSt-1),0.0,0); } // commented out on 2010/03/10 11:56:45 -- this messes up the requested interpolation range
//
//		x=from;
//
//		while (x<to) {
//			tmpfn.newPoint(x,0.0);
//			x+=dx;
//		}
//		if (abs(x-to)<acc) tmpfn.newPoint(to,0.0);
//		if (iEn+1<_f.count()) { tmpfn.newPoint(X(iEn+1),0.0); }
//
//	}
//
//	msgs->say("interpolating from "+msgs->toStr(from,15)+" to "+msgs->toStr(to,15)+" | from argument: "+msgs->toStr(iSt)+" to argument "+msgs->toStr(iEn),Low);
//	msgs->say("sampling step is:"+msgs->toStr(dx,15),Low);
//
//	//
//	// define the X values to interpolate on
//	//
////			tmpfn.print();
//	double *Xint=tmpfn.extractArguments();
//
//	//
//	// interpolate function
//	//
//	if (type=="auto") {
//		if (pointsCount()>5) tmpfn=interpolate(Xint,tmpfn.pointsCount(),"cspline", inPlace);
//		else tmpfn=interpolate(Xint,tmpfn.pointsCount(),"linear", inPlace);
//	}
//	if (type=="linear") { tmpfn=interpolate(Xint,tmpfn.pointsCount(),"linear", inPlace); }
//	if (type=="cspline") { tmpfn=interpolate(Xint,tmpfn.pointsCount(),"cspline", inPlace); }
//
//	// delete points outside of the requested interpolation region
//	if (!exact) tmpfn.deleteOutsideOf(from,to);
//
//	delete Xint;
//	tmpfn.checkRanges();
//	msgs->say("points count after interpolation: "+msgs->toStr(pointsCount()),Low);
//	return tmpfn;
//	
//}
/* **************************************************************************************************** */
const mscsFunction mscsFunction::interpolate(double from, double to, double dx, bool inPlace, string type, bool exact, double acc) {
	long iSt,iEn;
	double x,st,en;
	mscsFunction tmpfn("interpolating arguments",Zero);
	double from_orig=from, to_orig=to;
	
	//
	// define ranges and interpolation grid
	//
	f(from,&iSt);
	f(to,&iEn);
	if (from+acc < X(iSt) and iSt==0) msgs->warning("the chosen lower interpolation bound is smaller than the lower range over which the function is defined. I will not interpolate there.\nFrom is: "+msgs->toStr(from)+", first argument is: "+msgs->toStr(X(iSt)),High);
	if (to-acc > X(iEn) and iEn==pointsCount()-1) msgs->warning("the chosen upper interpolation bound is larger than the upper range over which the function is defined. I will not interpolate there.\nTo is: "+msgs->toStr(to)+", last argument is: "+msgs->toStr(X(iEn)),High);
	
	//			printf("iSt %li iEn %li\n", iSt, iEn);
	if (exact) { // set from and to to the first and the last respectively valid values inside of the region where function is defined.
		if (iSt==0) { st=X(iSt); } else { st=from; }
		
		x=from;
		while (x<st) { if (abs(st-x)<acc) x=st; else x+=dx; }
		from=x;
		
		x=from;
		en=cpeds_get_min(X(iEn),to);
		//				printf("from %lE to %lE en %lE\n", from, to,en);
		double en2=en+acc;
		while (x<en2) { if (abs(x-en)<acc) { x=en; tmpfn.newPoint(x,0.0); x=en2; printf("acc\n"); } else { tmpfn.newPoint(x,0.0); x+=dx; } }
		//if (abs(x-en)<acc) { x=en; tmpfn.newPoint(x,0.0); }
		if (acc==0) x-=dx;
		if (x==en2) x=en;
		to=x;
		//				printf("from %lE to %lE\n", from, to);
	}
	else {
		from=getX(iSt);
		if (to-acc > X(iEn) and iEn<pointsCount()-1) { iEn++; }
		to=getX(iEn);
		// set wider range to include derivatives information for the interpolation
		if (iSt>=1) { tmpfn.insertPoint(getX(iSt-1),0.0,0); } // commented out on 2010/03/10 11:56:45 -- this messes up the requested interpolation range
		
		x=from;
		
		while (x<to) {
			tmpfn.newPoint(x,0.0);
			x+=dx;
		}
		if (abs(x-to)<acc) tmpfn.newPoint(to,0.0);
		if (iEn+1<_f.count()) { tmpfn.newPoint(X(iEn+1),0.0); }
		
	}
	
	msgs->say("interpolating from "+msgs->toStr(from,15)+" to "+msgs->toStr(to,15)+" | from argument: "+msgs->toStr(iSt)+" to argument "+msgs->toStr(iEn),Low);
	msgs->say("sampling step is:"+msgs->toStr(dx,15),Low);
	msgs->say(string("interpolation type: ")+type,Low);
	
	//
	// define the X values to interpolate on
	//
	//			tmpfn.print();
	double *Xint=tmpfn.extractArguments();
	
	//
	// interpolate function
	//
	if (type=="auto") {
		if (pointsCount()>5) tmpfn=interpolate(Xint,tmpfn.pointsCount(),"cspline", inPlace);
		else tmpfn=interpolate(Xint,tmpfn.pointsCount(),"linear", inPlace);
	}
	if (type=="linear") { tmpfn=interpolate(Xint,tmpfn.pointsCount(),"linear", inPlace); }
	if (type=="cspline") { tmpfn=interpolate(Xint,tmpfn.pointsCount(),"cspline", inPlace); }
	if (type=="cspline_periodic") { tmpfn=interpolate(Xint,tmpfn.pointsCount(),"cspline_periodic", inPlace); }
	if (type=="akima") { 
		if (pointsCount()>5) tmpfn=interpolate(Xint,tmpfn.pointsCount(),"akima", inPlace); 
		else tmpfn=interpolate(Xint,tmpfn.pointsCount(),"linear", inPlace); 
	}
	if (type=="akima_periodic") { 
		tmpfn=interpolate(Xint,tmpfn.pointsCount(),"akima_periodic", inPlace); 
	}
	
	// delete points outside of the requested interpolation region
	if (!exact) { tmpfn.deleteOutsideOf(from_orig,to_orig);  if (inPlace) deleteOutsideOf(from_orig,to_orig); }
	
	
	delete [] Xint;
	tmpfn.checkRanges();
	msgs->say("points count after interpolation: "+msgs->toStr(tmpfn.pointsCount()),Low);
	return tmpfn;
}
/***************************************************************************************/
const mscsFunction mscsFunction::interpolate(double from, double to, long N, bool inPlace, string type) {
	mscsFunction tmpfn("interpolating arguments",Zero);
	//
	// define the X values to interpolate on
	//
	//			tmpfn.print();
	double *Xint=new double[N];
	double dx=(to-from)/N;
	for (unsigned long i = 0; i < N; i++) {
		Xint[i]=from+dx*i;
	}
	//	printf("dx: %lE\n",dx);
	
	//
	// interpolate function
	//
	if (type=="auto") {
		if (pointsCount()>5) tmpfn=interpolate(Xint,N,"cspline", inPlace);
		else tmpfn=interpolate(Xint,N,"linear", inPlace);
	}
	else {
		tmpfn=interpolate(Xint,N,type, inPlace); 
	}
	
	
	
	delete [] Xint;
	tmpfn.checkRanges();
	msgs->say("points count after interpolation: "+msgs->toStr(tmpfn.pointsCount()),Low);
	return tmpfn;
	
}
/* **************************************************************************************************** */
const mscsFunction mscsFunction::interpolate(double dx, bool inPlace, string type) {
	msgs->say("Interpolating",Medium);
	checkRanges();
	printRanges();
	return interpolate(getMinArg(),getMaxArg(),dx,inPlace,type,false);
}

/* **************************************************************************************************** */
const mscsFunction mscsFunction::extrapolate(double from, double to, double dx, bool inPlace, string type) {
	mscsFunction tmpfn("",Zero);
	if (pointsCount()==0) return tmpfn;
	double x,y,a,b;
	long i;
	checkRanges();
	msgs->say("extrapolating from "+msgs->toStr(from,15)+" to "+msgs->toStr(to,15)+" with step "+msgs->toStr(dx,15)+" type="+type,Low);
	
	if (pointsCount()==1) { // do const extrapolation
		x=from;
		while (x<getMinArg()) {      tmpfn.addPoint(x,Y(0));      x+=dx;    }
		x=getMaxArg()+dx;
		while (x<=to) { tmpfn.addPoint(x,Y(0)); x+=dx;  }
		if (inPlace) { concatenate(tmpfn); sortFunctionArgAscending(); return *this; }
		else { tmpfn.concatenate(*this); tmpfn.sortFunctionArgAscending(); return tmpfn; }
		
	}
	
	if (type=="linear") {
		if (pointsCount() > 1) { // do linear extrapolation
			x=from;
			while (x<getMinArg()) {
				f(x,&i);      tmpfn.addPoint(x,cpeds_extrapolate_linear(x,X(i),Y(i),X(i+1),Y(i+1)));            x+=dx;        }
			x=getMaxArg()+dx;
			while (x<=to) {
				f(x,&i);      tmpfn.addPoint(x,cpeds_extrapolate_linear(x,X(i-1),Y(i-1),X(i),Y(i)));            x+=dx;        }
			
			if (inPlace) { concatenate(tmpfn); sortFunctionArgAscending(); return *this; }
			else { tmpfn.concatenate(*this); tmpfn.sortFunctionArgAscending(); return tmpfn; }
		}
	}
	return tmpfn;
}
/* **************************************************************************************************** */
const mscsFunction& mscsFunction::concatenate(const mscsFunction& f) {
	long num=f.pointsCount();
	for (long i=0;i<num;i++) {    newPoint(f.getQPoint(i));  }
	return *this;
}
/* **************************************************************************************************** */
mscsFunction& mscsFunction::paste(const mscsFunction& f, long k, bool overwrite) {
	long num=f.pointsCount();
	long j;
	if (k<0) {	k=pointsCount()+k;	}
	if (overwrite) {
		for (long i=0;i<num;i++) {    j=k+i; if (j<pointsCount()) _f[j]=f.getQPoint(i); else { newPoint(f.getQPoint(i));  }	}
	}
	else {
		for (long i=0;i<num;i++) {    _f.insert(k+i,f.getQPoint(i));  }
	}
	return *this;
}

/* **************************************************************************************************** */
mscsFunction mscsFunction::cut(long N, long k) {
	mscsFunction c;
	if (N==0) return c;
	
	if (k<0) {	k=pointsCount()+k;	}
	if (k<0) { k=0; }
	if (N>pointsCount()-k) N=pointsCount()-k;
	for (long i=0;i<N;i++) {    c.newPoint(_f.takeAt(k));  }		
	return c;
}
/***************************************************************************************/
mscsFunction mscsFunction::copy(long N, long k) {
	mscsFunction c;
	if (N==0) return c;
	
	if (k<0) {	k=pointsCount()+k;	}
	if (k<0) { k=0; }
	if (N>pointsCount()-k) N=pointsCount()-k;
	for (long i=0;i<N;i++) {    c.newPoint(_f.at(k+i));  }		
	return c;
}
/* **************************************************************************************************** */
const mscsFunction& mscsFunction::unique() {
	long num=_f.count();
	long k;
	for (long i=0;i<num;i++) {
		k=i+1;
		for (long j=k;j<num;j++) {
			if (getQPoint(i)==getQPoint(j)) { deletePoint(j); num--;  j--; }
		}
	}
	checkRanges();
	return *this;
}

/* **************************************************************************************************** */
const mscsFunction& mscsFunction::uniqueXFast(double acc) {
	sortFunctionArgAscending();
	long num=_f.count()-1;
	long ip;
	for (long i=0;i<num;i++) {
		ip=i+1;
		if (X(ip)-X(i)<=acc) { deletePoint(ip); num--; i--; }
	}
	return *this;
}
/***************************************************************************************/
const mscsFunction& mscsFunction::uniqueYFast(double acc) {
	invert();
	uniqueXFast(acc);
	invert();
//	sortFunctionArgAscending();
//	long num=_f.count()-1;
//	long ip;
//	for (long i=0;i<num;i++) {
//		ip=i+1;
//		//		printf("f: %lE, fp: %lE, i: %li\n",f(ip),f(i),i);
//		if (fabs(f(ip)-f(i))<=acc) { deletePoint(ip); num--; i--; } //printf("  deleting\n"); }
//	}
	return *this;	
}
/* **************************************************************************************************** */
void mscsFunction::printRanges() const {
	msgs->say("Function ranges",Low);
	msgs->say("Xmin: "+msgs->toStr(getMinArg(),15)+" Xmax: "+msgs->toStr(getMaxArg(),15)+" Ymin: "+msgs->toStr(getMinValue(),15)+" Ymax: "+msgs->toStr(getMaxValue(),15),Low);
	msgs->say("iXmin: "+msgs->toStr(getMinArgIdx())+" iXmax: "+msgs->toStr(getMaxArgIdx())+" iYmin: "+msgs->toStr(getMinValueIdx())+" iYmax: "+msgs->toStr(getMaxValueIdx()),Low);
}
/* **************************************************************************************************** */
mscsFunction& mscsFunction::invert()  {
	long num=_f.count();
	double tmp;
	for (long i=0;i<num;i++) {
		tmp=getX(i);
		setarg(i,f(i));
		setf(i,tmp);
	}
	return *this;
}
/***************************************************************************************/
mscsFunction& mscsFunction::convertPhiRToReIm() {
	long num=_f.count();
	double phi,z;
	for (long i=0;i<num;i++) {
		phi=X(i);
		z=f(i);
		X(i)=z*cos(phi);
		f(i)=z*sin(phi);
	}
	return *this;
}
/***************************************************************************************/
mscsFunction& mscsFunction::convertReImToPhiR() {
	long num=_f.count();
	double x,y;
	for (long i=0;i<num;i++) {
		x=X(i); y=f(i);
		X(i)=cpeds_cart2sph(1,x,y,0);
		f(i)=sqrt(x*x+y*y);
	}
	return *this;
}

/***************************************************************************************/
mscsFunction mscsFunction::extractSpikes(long n, long m, long nsigma, double startFrom, long direction, bool startFromAsIdx) {
	mscsFunction spikes;
	if (n>pointsCount()) n=pointsCount();
	long N=pointsCount()/m;
	//	double* x=new double[n];
	double* y=new double[n];
	double stdev,thres,mean;
	long ist,ien,i,j,size;
	
	bool go=true,stop=false;
	
	if (direction==12) {
		if (startFromAsIdx) { ist=long(startFrom); } else { f(startFrom,&ist); }
		ien=ist+n-1;
	}
	else {
		if (startFromAsIdx) { ien=pointsCount()-long(startFrom)-1; } else { f(startFrom,&ien); }
		ist=ien-n+1;
	}
	
	mscsFunction ids; // indexes of points to be removed
	
	while (go) {
		if 	(stop) go=false;
		for (i=ist;i<=ien;i++) { j=i-ist; y[j]=f(i); } // copy the window
		// derive stats on the window
		size=ien-ist+1;
		stdev=sqrt(cpeds_variance(y,size));
		mean=cpeds_mean_value(y,size);
		thres=nsigma*stdev;
		// remove the outstanding points
		for (i=ist;i<=ien;i++) { if (fabs(f(i)-mean)>thres) ids.newPoint(double(i),0.0);  }
		// shift the window 
		//		if (direction==12) {	ist+=ien+1;		ien=ist+m; } else { ien=ist-1; ist=ien-m; }
		if (direction==12) {	ist+=m;		ien+=m; } else { ist-=m; ien-=m; }
		// check the ranges
		if (direction==12) { if (ien>=pointsCount()) { ien=pointsCount()-1; /*ist=ien-n+1;*/ stop=true; } }
		else  { if (ist<0) { ist=0; /*ien=ist+n-1;*/ stop=true; } }
#ifdef DEBUG
		printf("shifting window ist: %li ien: %li pointsCount: %li, m: %li, n: %li\n",ist,ien,pointsCount(),m,n);
#endif
	}
	
	//	printf("ids points: %li\n",ids.pointsCount());
	//	ids.uniqueXFast();
	ids.unique();
	//	printf("ids points: %li\n",ids.pointsCount());
	ids.sortFunctionArgAscending();
	ien=ids.pointsCount();
	for (i=0;i<ien;i++) { spikes.newPoint(X(long(ids.X(i))),f(long(ids.X(i)))); }
	for (i=0;i<ien;i++) { deletePoint(long(ids.X(ien-i-1))); }
	
	delete [] y;
	//	delete [] x;
	return spikes;
}
/***************************************************************************************/
mscsFunction mscsFunction::extractSpikesRange(long n, long m, long nsigma, double startFrom,double endAt, long direction, bool startFromAsIdx) {
	mscsFunction spikes;
	if (n>pointsCount()) n=pointsCount();
	long N=pointsCount()/m;
	//	double* x=new double[n];
	double* y=new double[n];
	double stdev,thres,mean;
	long ist,ien,i,j,size;
	long iEnd;
	
	bool go=true,stop=false;
	
	if (direction==12) {
		if (startFromAsIdx) { ist=long(startFrom); } else { f(startFrom,&ist); }
		ien=ist+n-1;
		if (startFromAsIdx) { iEnd=long(endAt); } else { f(endAt,&iEnd); }
	}
	else {
		if (startFromAsIdx) { ien=pointsCount()-long(startFrom)-1; } else { f(startFrom,&ien); }
		ist=ien-n+1;
		if (startFromAsIdx) { iEnd=pointsCount()-long(endAt)-1; } else { f(endAt,&iEnd); }
	}
	
	mscsFunction ids; // indexes of points to be removed
	
	while (go) {
		if 	(stop) go=false;
		for (i=ist;i<=ien;i++) { j=i-ist; y[j]=f(i); } // copy the window
		// derive stats on the window
		size=ien-ist+1;
		stdev=sqrt(cpeds_variance(y,size));
		mean=cpeds_mean_value(y,size);
		thres=nsigma*stdev;
		// remove the outstanding points
		for (i=ist;i<=ien;i++) { if (fabs(f(i)-mean)>thres) ids.newPoint(double(i),0.0);  }
		// shift the window 
		//		if (direction==12) {	ist+=ien+1;		ien=ist+m; } else { ien=ist-1; ist=ien-m; }
		if (direction==12) {	ist+=m;		ien+=m; } else { ist-=m; ien-=m; }
		// check the ranges
		if (direction==12) { if (ien>=iEnd) { ien=iEnd-1; /*ist=ien-n+1;*/ stop=true; } }
		else  { if (ist<iEnd) { ist=iEnd; /*ien=ist+n-1;*/ stop=true; } }
#ifdef DEBUG
		printf("shifting window ist: %li ien: %li pointsCount: %li, m: %li, n: %li iEnd: %li\n",ist,ien,pointsCount(),m,n,iEnd);
#endif
	}
	
	//	printf("ids points: %li\n",ids.pointsCount());
	//	ids.uniqueXFast();
	ids.unique();
	//	printf("ids points: %li\n",ids.pointsCount());
	ids.sortFunctionArgAscending();
	ien=ids.pointsCount();
	for (i=0;i<ien;i++) { spikes.newPoint(X(long(ids.X(i))),f(long(ids.X(i)))); }
	for (i=0;i<ien;i++) { deletePoint(long(ids.X(ien-i-1))); }
	
	delete [] y;
	//	delete [] x;
	return spikes;
	
}


/* **************************************************************************************************** */
double mscsFunction::integrate() {
	double *X=extractArguments();
	double *Y=extractValues();
	double integral = cpeds_integrate_1d(pointsCount(),X,Y);
	delete [] X;
	delete [] Y;
	return integral;
}
/***************************************************************************************/
double mscsFunction::integrate(double xmin,double xmax) {
	long imin, imax;
	
	f(xmin,&imin);
	f(xmax,&imax);
	
	//	printf("imin: %li, imax: %li\n",imin,imax);
	double *X=extractArguments(imin,imax);
	double *Y=extractValues(imin,imax);
	//	if (imin==2) exit(0);
	//	if (X==NULL) {printf("X null\n"); exit(0); }
	double integral = cpeds_integrate_1d(imax-imin+1,X,Y);
	delete [] X;
	delete [] Y;
	return integral;	
}

/* **************************************************************************************************** */
const mscsFunction& mscsFunction::setBelow(double th, double val) {
	long num=_f.count();
	for (long i=0;i<num;i++) {
		if (f(i) < th) { f(i)=val; }
	}
	return *this;
}
/* **************************************************************************************************** */
const mscsFunction& mscsFunction::setAbove(double th, double val) {
	long num=_f.count();
	for (long i=0;i<num;i++) {
		if (f(i) > th) { f(i)=val; }
	}
	return *this;
}
/***************************************************************************************/
mscsFunction& mscsFunction::setRange(double from, double to, double val) {
	for (long i=0;i<_f.count();i++) { if (X(i)>=from and X(i)<=to) { f(i)=val; }  }
	return *this;
}
/* **************************************************************************************************** */
const mscsFunction& mscsFunction::setf(double val) {
	long num=_f.count();
	for (long i=0;i<num;i++) { f(i)=val; }
	return *this;
}
/* **************************************************************************************************** */
mscsFunction& mscsFunction::binFunction(long imin, cpedsList<long>& binSize, cpedsList<double>& w) {
	//  mscsAngularPowerSpectrum bin_C_l(mscsAngularPowerSpectrum* c_l, long lmin, long lmax, long *bintabs, long* bintab, double w) {
	long i,j,k,binPointsNum=binSize.sum(),imax=imin+binPointsNum-1;
	long Cls, Cbs,cls,cbs;
	bool EQweights;
	double wsum,W=0;
	stringstream ss; ss<<"Binning the function from index: "<<imin<<" to index: "<<imax<<", together: "<<imax-imin+1<<" points.";
	msgs->say(ss.str(),High);
	msgs->say("  -- bins are:",Low);
	binSize.print();
	
	// check safety conditions
	if (imin<0) { msgs->error("negative imin arument: "+msgs->toStr(imin)+". Not binning",High);     return *this; }
	// requirement that binSize should be positive numbers
	if (!(binSize>0)) { msgs->error("negative bin size given. Check the binSize vector. Not binning",High);     return *this; }
	// if there's more bins than data
	if (imax>pointsCount()-1) {
		msgs->warning("The binning array is incompatible with requested binning multipole range. Will bin upto the end of the data dropping some bins from the binning vector",High);
		// printf("number of bins is: %li.\n",*bintabs);
		
		// fit the length of the binning vector to the actual data size if the supplied vector was too large.
		while (imax>pointsCount()-1) {
			if (binSize.last()>1) binSize.last()--; else { binSize.removeLast(); }
			if (w.size()!=0) w.removeLast();
			binPointsNum=binSize.sum();
			imax=imin+binPointsNum-1;
		}
		// printf("number of bins is : %li.\n",*bintabs);
		// (*bintabs)++;
		// bintab[*bintabs-1]=lmax_loc-lmin_loc+2-binl_sum;
		// printf("adding last bin of size: %li and increasing number of bins to: %li.\n",bintab[*bintabs-1], *bintabs);
		// binl_sum=0;
		// for (i=0;i<*bintabs;i++) { binl_sum+=bintab[i]; } printf("shrinking binl_sum to: %li\n",binl_sum);
		// printf("bin sum is now: %li.\n",binl_sum);
		
		msgs->say("Binning the function from index: "+msgs->toStr(imin)+" to index: "+msgs->toStr(imax)+", together: "+msgs->toStr(imax-imin+1),High);
		if (imin==imax) return *this; // no binning is done
	}
	
	
	// // find out the size of the binned C_b
	cls = pointsCount(); // original size of the C_l
	cbs = pointsCount()  - binPointsNum + binSize.size(); // total size of binned C_l: bs
	Cls = imax-imin+1; // size of original C_l for binning
	Cbs = binSize.size(); // size of binned C_l: bs
	if (w.size() == 0) { EQweights = true;  } else { EQweights = false; }
	msgs->say("Using same weights for all points: "+msgs->toStr(EQweights),Medium);
	if (w.size() !=0 && w.size() != Cls) {   msgs->error("Weights table size ("+msgs->toStr(long(w.size()))+") doesn't match the size of the data to be binned ("+msgs->toStr(Cls),Medium); }
	
	
	// printf("  -- total size of cl vector before binning: %li\n",cls);
	// printf("  -- size of Cl vector for binning: %li\n",Cls);
	// printf("  -- total size of cb - binned Cl: %li\n",cbs);
	// printf("  -- size of Cb binned vector: %li\n",Cbs);
	
	
	// define the structures needed
	matrix <double> Cl(Cls,1); // original C_l vector for binning
	matrix <double> x(Cls,1); // original function arguments in the binning range
	matrix <double> Cb(Cbs,1); // binned C_b vector
	matrix <double> M(Cbs,Cls); // binning operator
	cpedsList<double> effl;
	
	// // copy the power spectra outside the reigon for binning
	// for (i=0;i<lmin_loc;i++) { cb->newPoint(getX(i), getY(i)); } // lower end
	// tmp = lmax_loc+1; j=lmin_loc+(*bintabs);
	// for (i=tmp;i<=lmax;i++)  { cb->set_l(j,c_l->get_l(i)); cb->set_Cl(j,c_l->get_Cl(i)); j++;} // upper end shifted as to fit the binned Cb in range lmin_loc, lmax_loc
	
	// prepare the power spectrum vector for binning
	for (j=Cls;j>0;j--) {  Cl(j-1,0) = getY(j-1+imin); x(j-1,0)=getX(j-1+imin); deletePoint(j-1+imin); } // copy the power spectra part for binning and remove the points that will be binned
	
	
	// prepare the binning matrix operator
	
	for (i=0;i<Cbs;i++)
		for (j=0;j<Cls;j++)
			M(i,j) = 0;  // zero to binning matrix
	
	k=0;
	for (i=0;i<Cbs;i++) {
		//if (i==0) jst = i; else jst = bintab[i]-1;
		if (EQweights) W = 1/(double)binSize[i];
		effl.append(0); wsum=0;
		for (j=k;j<k+binSize[i];j++) {
			if (!EQweights) { W = w[j]; wsum+=W; }
			M(i,j) = W;
			effl[i]+=x(j,0)*W;
		}
		if (!EQweights) { for (j=k;j<k+binSize[i];j++) { M(i,j)=M(i,j)/wsum; } effl[i]/=wsum; }
		// printf("effective ls: %lf\n",effl[i]);
		k+=binSize[i];
	}
	
	// do the binning
	Cb = M*Cl;
	
	// rewrite the Cb onto the object
	
	for (i=0;i<Cbs;i++)  {
		newPoint(effl[i], Cb(i,0));
		// printf("l %lf Cl %lE            binned Cb: l %lf Cb %lE\n",effl[i],Cb(i,0),cb->get_l(j),cb->get_Cl(j));
	}
	
	// delete [] effl;
	
	// return cb;
	
	// sort function
	sortFunctionArgAscending();
	return *this;
}
/***************************************************************************************/
mscsFunction mscsFunction::binFunction(double dx, cpedsList<long>& binSize, string methodX, string methodY) {
	mscsFunction fbinned;
	long i,ist, ien, binidx;
	double x,binst,binen, xmin,xmax,binArg, binVal;
	cpedsList<double> binX,binY;
	
	
	checkRanges();
	binidx=0;
	ist=0;
	ien=pointsCount();
	x=getx(ist);
	i=ist;

	do {
		binX.clear();
		binY.clear();

		binst=getMinArg()+binidx*dx;
		binen=binst+dx;
		x=getx(i);
//		printf("binst: %lf binen: %lf, x=%lf i=%li\n",binst,binen,x,i);
		while (x>=binst and x<binen) {
			binX.append(x);
			binY.append(f(i));
			i++;
			if (i==pointsCount()) break; 
			x=getx(i);
		}
		
		if (binX.size()>0) {
			if (methodX=="mean") binArg=binX.mean();
			else {
				if (methodX=="median") binArg=binX.median();
				else {
					if (methodX=="min") binArg=binX.min();
					else {
						if (methodX=="max") binArg=binX.max();
						else {
							if (methodX=="bin_center") binArg=binst+dx/2;
							else {
								if (methodX=="bin_max") binArg=binen;
								else {
									if (methodX=="bin_min") binArg=binst;
									else {
										cout << "Wrong methodX in binFunction\n"; exit(1);
									}
								}
							}
						}
					}
				}
			}
	
			if (methodY=="mean") binVal=binY.mean();
			else {
				if (methodY=="median") binVal=binY.median();
				else {
					if (methodY=="min") binVal=binY.min();
					else {
						if (methodY=="max") binVal=binY.max();
						else {
							cout << "Wrong methodY in binFunction\n"; exit(1);
						}
					}
				}
			}
		
			binSize.append(binX.size());
			fbinned.newPoint(binArg,binVal);
		}
		binidx++;
		
	} while (binen<=getMaxArg());
	
	return fbinned;
}
/* ******************************************************************************************** */
mscsFunction mscsFunction::binFunction(cpedsList<double>& binEdge, cpedsList<long>& binCounts, string methodX, string methodY) {
	mscsFunction fbinned;
	if (binEdge.size()<2) return fbinned; 
	long i,ist, ien, binidx;
	double x,binst,binen, binctr, xmin,xmax,binArg, binVal;
	cpedsList<double> binX,binY;
	
	
	checkRanges();
	binidx=0;
	ist=0;
	ien=pointsCount();
	x=getx(ist);
	i=ist;

	do {
		binX.clear();
		binY.clear();

		binst=binEdge[binidx];
		binen=binEdge[binidx+1];
		binctr=(binst+binen)/2;
		x=getx(i);
		while (x>=binst and x<binen) {
			binX.append(x);
			binY.append(f(i));
			i++;
			if (i==pointsCount()) break; 
			x=getx(i);
//			printf("----------\n");
//			binY.print();
		}
//		printf("----------\n");
//		binY.print();
		
		if (binX.size()>0) {
			if (methodX=="mean") binArg=binX.mean();
			else {
				if (methodX=="median") binArg=binX.median();
				else {
					if (methodX=="min") binArg=binX.min();
					else {
						if (methodX=="max") binArg=binX.max();
						else {
							if (methodX=="bin_center") binArg=binctr;
							else {
								if (methodX=="bin_max") binArg=binen;
								else {
									if (methodX=="bin_min") binArg=binst;
									else {
										cout << "Wrong methodX in binFunction\n"; exit(1);
									}
								}
							}
						}
					}
				}
			}
	
			if (methodY=="mean") binVal=binY.mean();
			else {
				if (methodY=="median") binVal=binY.median();
				else {
					if (methodY=="min") binVal=binY.min();
					else {
						if (methodY=="max") binVal=binY.max();
						else {
							cout << "Wrong methodY in binFunction\n"; exit(1);
						}
					}
				}
			}
//			binY.print();
			binCounts.append(binX.size());
			fbinned.newPoint(binArg,binVal);
		}
		else {
			// make sure binCounts has the right size when there is no data in the bin
			binCounts.append(0);
		}
		
		if (x>=binen) binidx++;
		if (x<binst) i++;
		if (i==pointsCount()) break; 
		
	} while (binidx<binEdge.size()-1);
	
	return fbinned;
	
}
/***************************************************************************************/
mscsFunction& mscsFunction::binFunctionLin(long imin, cpedsList<long>& binSize, cpedsList<double>& w) {
	
	//  mscsAngularPowerSpectrum bin_C_l(mscsAngularPowerSpectrum* c_l, long lmin, long lmax, long *bintabs, long* bintab, double w) {
	long i,j,k,binPointsNum=binSize.sum(),imax=imin+binPointsNum-1;
	long Cls, Cbs,cls,cbs;
	bool EQweights;
	double wsum,W=0;
	stringstream ss; ss<<"Binning the function from index: "<<imin<<" to index: "<<imax<<", together: "<<imax-imin+1<<" points.";
	
	msgs->say(ss.str(),High);
	msgs->say("  -- bins are:",Low);
	binSize.print();
	
	// check safety conditions
	if (imin<0) { msgs->error("negative imin arument: "+msgs->toStr(imin)+". Not binning",High);     return *this; }
	// requirement that binSize should be positive numbers
	if (!(binSize>0)) { msgs->error("negative bin size given. Check the binSize vector. Not binning",High);     return *this; }
	// if there's more bins than data
	if (imax>pointsCount()-1) {
		msgs->warning("The binning array is incompatible with requested binning multipole range. Will bin upto the end of the data dropping some bins from the binning vector",High);
		// printf("number of bins is: %li.\n",*bintabs);
		
		// fit the length of the binning vector to the actual data size if the supplied vector was too large.
		while (imax>pointsCount()-1) {
			if (binSize.last()>1) binSize.last()--; else { binSize.removeLast(); }
			if (w.size()!=0) w.removeLast();
			binPointsNum=binSize.sum();
			imax=imin+binPointsNum-1;
		}
		// printf("number of bins is : %li.\n",*bintabs);
		// (*bintabs)++;
		// bintab[*bintabs-1]=lmax_loc-lmin_loc+2-binl_sum;
		// printf("adding last bin of size: %li and increasing number of bins to: %li.\n",bintab[*bintabs-1], *bintabs);
		// binl_sum=0;
		// for (i=0;i<*bintabs;i++) { binl_sum+=bintab[i]; } printf("shrinking binl_sum to: %li\n",binl_sum);
		// printf("bin sum is now: %li.\n",binl_sum);
		
		msgs->say("Binning the function from index: "+msgs->toStr(imin)+" to index: "+msgs->toStr(imax)+", together: "+msgs->toStr(imax-imin+1),High);
		if (imin==imax) return *this; // no binning is done
	}
	
	
	// // find out the size of the binned C_b
	cls = pointsCount(); // original size of the C_l
	cbs = pointsCount()  - binPointsNum + binSize.size(); // total size of binned C_l: bs
	Cls = imax-imin+1; // size of original C_l for binning
	Cbs = binSize.size(); // size of binned C_l: bs
	if (w.size() == 0) { EQweights = true;  } else { EQweights = false; }
	msgs->say("Using same weights for all points: "+msgs->toStr(EQweights),Medium);
	if (w.size() !=0 && w.size() != Cls) {   msgs->error("Weights table size ("+msgs->toStr(long(w.size()))+") doesn't match the size of the data to be binned ("+msgs->toStr(Cls),Medium); }
	
	
	// printf("  -- total size of cl vector before binning: %li\n",cls);
	// printf("  -- size of Cl vector for binning: %li\n",Cls);
	// printf("  -- total size of cb - binned Cl: %li\n",cbs);
	// printf("  -- size of Cb binned vector: %li\n",Cbs);
	
	
	// define the structures needed
	matrix <double> Cl(Cls,1); // original C_l vector for binning
	matrix <double> x(Cls,1); // original function arguments in the binning range
	matrix <double> Cb(Cbs,1); // binned C_b vector
	matrix <double> Cbtmp(Cbs,1); // binned C_b vector
	matrix <double> M(1,Cls); // binning operator
	cpedsList<double> effl;
	
	// // copy the power spectra outside the reigon for binning
	// for (i=0;i<lmin_loc;i++) { cb->newPoint(getX(i), getY(i)); } // lower end
	// tmp = lmax_loc+1; j=lmin_loc+(*bintabs);
	// for (i=tmp;i<=lmax;i++)  { cb->set_l(j,c_l->get_l(i)); cb->set_Cl(j,c_l->get_Cl(i)); j++;} // upper end shifted as to fit the binned Cb in range lmin_loc, lmax_loc
	
	// prepare the power spectrum vector for binning
	for (j=Cls;j>0;j--) {  Cl(j-1,0) = getY(j-1+imin); x(j-1,0)=getX(j-1+imin); deletePoint(j-1+imin); } // copy the power spectra part for binning and remove the points that will be binned
	
	
	// prepare the binning matrix operator
	
	//	for (i=0;i<Cbs;i++)
	
	k=0;
	for (i=0;i<Cbs;i++) {
		for (j=0;j<Cls;j++)	M(0,j) = 0;  // zero the binning matrix
		//if (i==0) jst = i; else jst = bintab[i]-1;
		if (EQweights) W = 1/(double)binSize[i];
		effl.append(0); wsum=0;
		for (j=k;j<k+binSize[i];j++) {
			if (!EQweights) { W = w[j]; wsum+=W; }
			M(0,j) = W;
			effl[i]+=x(j,0)*W;
		}
		if (!EQweights) { for (j=k;j<k+binSize[i];j++) { M(0,j)=M(0,j)/wsum; } effl[i]/=wsum; }
		// printf("effective ls: %lf\n",effl[i]);
		k+=binSize[i];
		
		// do the binning
		Cbtmp = M*Cl;
		Cb(i,0)=Cbtmp(0,0);
		printf("binning bin: %li of %li\n",i,Cbs);
	}
	
	
	// rewrite the Cb onto the object
	
	for (i=0;i<Cbs;i++)  {
		newPoint(effl[i], Cb(i,0));
		// printf("l %lf Cl %lE            binned Cb: l %lf Cb %lE\n",effl[i],Cb(i,0),cb->get_l(j),cb->get_Cl(j));
	}
	
	// delete [] effl;
	
	// return cb;
	
	// sort function
	sortFunctionArgAscending();
	return *this;
	
	
}
/***************************************************************************************/
mscsFunction& mscsFunction::binFunctionLin2(long imin, cpedsList<long>& binSize, cpedsList<double>& w, cpedsList<double>* varInbin, cpedsList<double>* alpha, vector< cpedsList<double> >* stats) {
	
	//  mscsAngularPowerSpectrum bin_C_l(mscsAngularPowerSpectrum* c_l, long lmin, long lmax, long *bintabs, long* bintab, double w) {
	long i,j,k,binPointsNum=binSize.sum(),imax=imin+binPointsNum-1;
	long Cls, Cbs;//,cls,cbs;
	bool EQweights;
	double wsum,W=0;
	stringstream ss; ss<<"Binning the function from index: "<<imin<<" to index: "<<imax<<", together: "<<imax-imin+1<<" points.";
	
	msgs->say(ss.str(),High);
	//msgs->say("  -- bins are:",Low);
	//binSize.print();
	
	// check safety conditions
	if (imin<0) { msgs->error("negative imin arument: "+msgs->toStr(imin)+". Not binning",High);     return *this; }
	// requirement that binSize should be positive numbers
	if (!(binSize>0)) { msgs->error("negative bin size given. Check the binSize vector. Not binning",High);     return *this; }
	// if there's more bins than data
	if (imax>pointsCount()-1) {
		msgs->warning("The binning array is incompatible with requested binning multipole range. Will bin upto the end of the data dropping some bins from the binning vector",High);
		// printf("number of bins is: %li.\n",*bintabs);
		
		// fit the length of the binning vector to the actual data size if the supplied vector was too large.
		while (imax>pointsCount()-1) {
			if (binSize.last()>1) binSize.last()--; else { binSize.removeLast(); }
			if (w.size()!=0) w.removeLast();
			binPointsNum=binSize.sum();
			imax=imin+binPointsNum-1;
		}
		// printf("number of bins is : %li.\n",*bintabs);
		// (*bintabs)++;
		// bintab[*bintabs-1]=lmax_loc-lmin_loc+2-binl_sum;
		// printf("adding last bin of size: %li and increasing number of bins to: %li.\n",bintab[*bintabs-1], *bintabs);
		// binl_sum=0;
		// for (i=0;i<*bintabs;i++) { binl_sum+=bintab[i]; } printf("shrinking binl_sum to: %li\n",binl_sum);
		// printf("bin sum is now: %li.\n",binl_sum);
		
		msgs->say("Binning the function from index: "+msgs->toStr(imin)+" to index: "+msgs->toStr(imax)+", together: "+msgs->toStr(imax-imin+1),High);
		if (imin==imax) return *this; // no binning is done
	}
	
	
	// // find out the size of the binned C_b
	//	cls = pointsCount(); // original size of the C_l
	//	cbs = pointsCount()  - binPointsNum + binSize.size(); // total size of binned C_l: bs
	Cls = imax-imin+1; // size of original C_l for binning
	Cbs = binSize.size(); // size of binned C_l: bs
	if (w.size() == 0) { EQweights = true;  } else { EQweights = false; }
	msgs->say("Using same weights for all points: "+msgs->toStr(EQweights),Medium);
	if (w.size() !=0 && w.size() != Cls) {   msgs->error("Weights table size ("+msgs->toStr(long(w.size()))+") doesn't match the size of the data to be binned ("+msgs->toStr(Cls)+")",Medium); }
	
	
	// printf("  -- total size of cl vector before binning: %li\n",cls);
	// printf("  -- size of Cl vector for binning: %li\n",Cls);
	// printf("  -- total size of cb - binned Cl: %li\n",cbs);
	// printf("  -- size of Cb binned vector: %li\n",Cbs);
	
	
	// define the structures needed
	double *Cl = new double[Cls]; // original C_l vector for binning
	double *x = new double[Cls]; // original function arguments in the binning range
	double *Cb = new double[Cbs]; // binned C_b vector
	double Cbtmp; // binned range value
	double *M = new double[Cls]; // binning operator
	double *effl = new double[Cbs]; // effective l-value
	
	
	// // copy the power spectra outside the reigon for binning
	// for (i=0;i<lmin_loc;i++) { cb->newPoint(getX(i), getY(i)); } // lower end
	// tmp = lmax_loc+1; j=lmin_loc+(*bintabs);
	// for (i=tmp;i<=lmax;i++)  { cb->set_l(j,c_l->get_l(i)); cb->set_Cl(j,c_l->get_Cl(i)); j++;} // upper end shifted as to fit the binned Cb in range lmin_loc, lmax_loc
	
	// prepare the power spectrum vector for binning
	for (j=Cls;j>0;j--) {  Cl[j-1] = getY(j-1+imin); x[j-1]=getX(j-1+imin); deletePoint(j-1+imin); } // copy the power spectra part for binning and remove the points that will be binned
	

	cpedsList<double> significance;
	if (alpha!=NULL and stats!=NULL) {
		significance=*alpha;
		cpedsList<double> dummy;
		stats->resize(significance.size()+8,dummy);
	}

	
	
	//
	// prepare the binning matrix operator
	//
	
	for (j=0;j<Cls;j++)	M[j] = 0;  // zero the binning matrix
	k=0;
	double kmin,kmax;
	for (i=0;i<Cbs;i++) {
		if (EQweights) W = 1/(double)binSize[i];
		effl[i]=0; wsum=0;
		kmax=k+binSize[i];
		for (j=k;j<kmax;j++) { M[j]=0; }
		for (j=k;j<kmax;j++) {
			if (!EQweights) { W = w[j]; wsum+=W; }
			M[j] = W;
			effl[i]+=x[j]*W;
		}
		if (wsum==0) wsum=1; // BLmodification (Jul 13, 2011, 10:35:43 PM): protection against zero weights in all bins
		if (!EQweights) { for (j=k;j<kmax;j++) { M[j]/=wsum; } effl[i]/=wsum; }
		
		// do the binning
		Cbtmp=0;
		for (j=k;j<kmax;j++) { Cbtmp+= M[j]*Cl[j]; }
		Cb[i]=Cbtmp;
		
		// calculate the variance within the bin from the original vector
		if (varInbin!=NULL) varInbin->append(cpeds_variance(Cl,kmax-k,k));
		//		printf("binning bin: %li of %li\n",i+1,Cbs);
		//		printf(" - effective l: %lf\n",effl[i]);

		if (alpha!=NULL and stats!=NULL) {
			mscsFunction tmpbin;
			cpedsList<double> tmpl;
			long si=0;
			tmpl.fromCarray(&Cl[k],kmax-k);
			qSort(tmpl);

			stats->at(si++).append(cpeds_mean_value(x,kmax-k,k)); // mean x
			stats->at(si++).append(cpeds_mean_value(Cl,kmax-k,k)); // mean y
			stats->at(si++).append(sqrt(cpeds_variance(x,kmax-k,k))); // sigma x
			stats->at(si++).append(sqrt(cpeds_variance(Cl,kmax-k,k))); // sigma y
			stats->at(si++).append(cpeds_mean_value(x,kmax-k,k) - cpeds_find_min_value(x,kmax,k)); // delta x-
			stats->at(si++).append(cpeds_find_max_value(x,kmax,k) - cpeds_mean_value(x,kmax-k,k)); // delta x+
			stats->at(si++).append(cpeds_mean_value(Cl,kmax-k,k) - cpeds_find_min_value(Cl,kmax,k)); // delta y-
			stats->at(si++).append(cpeds_find_max_value(Cl,kmax,k) - cpeds_mean_value(Cl,kmax-k,k)); // delta y+
			
			for (long si2= 0; si2 < significance.size(); si2++) {
				long idx=tmpl.size()*significance[si2];
				stats->at(si++).append(tmpl[idx]);
			}
			
		}
		
		k+=binSize[i];
	}
	
	//	printf("finishing");
	
	//
	// rewrite the Cb onto *this function
	//
	for (i=0;i<Cbs;i++)  {
		insertPoint(effl[i], Cb[i],imin+i);
		// printf("l %lf Cl %lE            binned Cb: l %lf Cb %lE\n",effl[i],Cb(i,0),cb->get_l(j),cb->get_Cl(j));
	}
	
	
	// 
	// clean up
	//
	delete [] effl;
	delete [] x;
	delete [] Cl;
	delete [] Cb;
	delete [] M;
	
	
	// sort function
	//	sortFunctionArgAscending();
	return *this;
	
	
}
/***************************************************************************************/
mscsFunction& mscsFunction::binFunctionLin2(long imin, double bs, cpedsList<long>& binSize, cpedsList<double>& w) {
	
	binSize.clear();
	long i,j,k;
	double b=bs;
	i=imin;
	j=0;
	k=0;
	long fsize=pointsCount();
	while (i<fsize) {
		j++;
		i=long(round(j*b));
		if (i+imin<=fsize) { binSize.append(i-k); } else { binSize.append(fsize-k); }
		k=i;
	}
	
	return binFunctionLin2(imin,binSize,w);
	
}


/***************************************************************************************/
mscsFunction& mscsFunction::binFunctionLin2geo(long imin, double firstBin, double gm, cpedsList<long>& binSize, cpedsList<double>& w, cpedsList<double>* varInbin) {
	binSize.clear();
	long i;
	double b=firstBin;
	i=imin;
	long fsize=pointsCount();
	while (i<fsize) {
		i+=long(round(b));
		if (i<=fsize) { binSize.append(long(round(b))); } else { binSize.append(fsize-i+long(round(b))); }
		b*=gm;
	}
	
	return binFunctionLin2(imin,binSize,w,varInbin);
}

/***************************************************************************************/
mscsFunction mscsFunction::powerSpectrum(mscsFunction* re, mscsFunction* im, double regridX, bool calibrateByPointsCount) {
	mscsFunction P("power spectrum", getVerbosityLevel());
	if (re!=NULL) re->clearFunction();
	if (im!=NULL) im->clearFunction();
	long N;
	double dX;
	N=pointsCount();
	
	//
	// regridding
	//
	if (regridX==0) dX=(last().rx()-X(0))/pointsCount(); // assume that function is equally spaced - take separation as an average from the whole signal span
	if (regridX>0) { // calculate the mean separation of arguments
		double *X=extractArguments();
		long M=N-1;
		for (int i = 0; i < M; i++) { X[i]=X[i+1]-X[i]; }
		dX=regridX*cpeds_mean_value(X,M);
		interpolate(dX,true,"linear");
		delete [] X;
	}
	
	//
	// fft part
	//
	N=pointsCount();
	fftw_plan p;
	long i;
#ifdef DEBUG_USING_FFTW2AND3
	fftw_complex* t=new fftw_complex[N*sizeof(fftw_complex)]; 
#else
	fftw_complex* t=(fftw_complex*)fftw_malloc(N*sizeof(fftw_complex)); 
#endif
	int dir=FFTW_FORWARD;
#ifdef DEBUG
	struct timeval startTime;
	struct timeval endTime;
	struct rusage ru;
	getrusage(RUSAGE_SELF, &ru);
	startTime = ru.ru_utime;
	//	clock_t debug_time1,debug_time2;
	//	debug_time1=clock();
#endif
	
//#pragma omp critical
//	{
	p=fftw_plan_dft_1d(N, t,t, dir, FFTW_ESTIMATE); 
//	}
	for (i=0;i<N;i++) { t[i][0]=f(i); t[i][1]=0.0; }
	fftw_execute(p);
#ifdef DEBUG
	getrusage(RUSAGE_SELF, &ru);
	endTime = ru.ru_utime;
	getrusage(RUSAGE_SELF, &ru);
	double tS = startTime.tv_sec*1000000 + (startTime.tv_usec);
	double tE = endTime.tv_sec*1000000  + (endTime.tv_usec);
	FILE* F=fopen("mscsFunction_powerSpectrumFFT1d.timing","w"); fprintf(F,"%lf\n",(tE-tS)*1e-6); fclose(F);
	printf("FFT exec time [s]: %lf\n",(tE-tS)*1e-6);
	//	debug_time2=clock();
	//	FILE* F=fopen("mscsFunction_powerSpectrumFFT1d.timing","w"); fprintf(F,"%lf\n",double(debug_time2-debug_time1)/CLOCKS_PER_SEC); fclose(F);
	//	printf("FFT exec time [s]: %lf\n",double(debug_time2-debug_time1)/CLOCKS_PER_SEC);
#endif
	
	//
	// power spectrum
	//
	
	if (regridX==0) regridX=1;
	double tmp=double(N)*double(dX); // sampling period
	
	long M;
	double N2=N; N2*=N2;
	if (N % 2==0) M=N/2+1;
	else M=(N-1)/2+1;
	
	if (calibrateByPointsCount)
		for (i=0;i<M;i++) { t[i][0]/=N; t[i][1]/=N; }
	//			for (i=0;i<M;i++) { t[i][0]/=tmp; t[i][1]/=tmp; } 
	
	if (re!=NULL) {		for (i=0;i<M;i++) { re->newPoint(double(i)/tmp,t[i][0]); }	}
	if (im!=NULL) {		for (i=0;i<M;i++) { im->newPoint(double(i)/tmp,t[i][1]); }	}
	
	for (i=0;i<M;i++) { P.newPoint(double(i)/tmp,  (t[i][0]*t[i][0]  +  t[i][1]*t[i][1]) ); } 
	
	// this is an implementation that sums over the positive and negative frequencies
	//	for (i=0;i<M;i++) { P.newPoint(double(i)/tmp,  2.0*(t[i][0]*t[i][0]  +  t[i][1]*t[i][1]) ); } //May 19, 2011, 11:56:58 AM: the  factor of 2 was be added due to summation over positive and negative frequencies
	//	P.f(0)/=double(2.0); // zero'th frequency does not get the factor of 2 bacause its a real number
	//	if (N%2==0) // divide the Nyquist frequency by 2 when it doesn't have its negative counterpart (i.e. when N is even)
	//		P.f(M-1)/=double(2.0);
#pragma omp critical
	{
	fftw_destroy_plan(p);
	}
#ifdef DEBUG_USING_FFTW2AND3
	delete [] t;
#else
	fftw_free(t);
#endif
	fftw_cleanup();
	
	return P;
}

/***************************************************************************************/
mscsFunction mscsFunction::inverseFFT(mscsFunction& re, mscsFunction& im, double* x, bool deleteX) {
	long pointsNum;
	pointsNum=2*(re.pointsCount()-1);
	if (abs(im.getY(re.pointsCount()-1))>1e-13) pointsNum++;
	
	//	long M=re.pointsCount();
	//	double dX=(x[pointsNum-1]-x[0])/pointsNum; // assume that function is equally spaced - take separation as an average from the whole signal span
	//	if (calibrateBySamplingPeriod) {
	//		re*=dX;
	//		im*=dX;
	//	}
	
	double *sig=fft_1D_c2r(re.extractValues(),im.extractValues(),re.pointsCount(),pointsNum,true);
	clearFunction();
	importFunction(x,sig,pointsNum,false);
	if (deleteX and x!=0) delete [] x;
	delete [] sig;
	return *this;
}

/***************************************************************************************/
mscsFunction mscsFunction::powerSpectrum(double kmin, double kmax, double dk,mscsFunction* Re, mscsFunction* Im) {
	mscsFunction re,im;
	double k=kmin;
	long j=0;
	long n=pointsCount();
	double *x=extractArguments();
	double *y=extractValues();
	
	double *yy=new double[n];
#ifdef DEBUG
	struct timeval startTime;
	struct timeval endTime;
	struct rusage ru;
	getrusage(RUSAGE_SELF, &ru);
	startTime = ru.ru_utime;
	//	clock_t debug_time1,debug_time2;
	//	debug_time1=clock();
#endif
	while (k<=kmax) {
		re.newPoint(k,0);
		im.newPoint(k,0);
#pragma omp parallel for shared(y,yy,x,n)
		for (long i = 0; i < n; i++) {	yy[i]=cos(twoPI*k*x[i])*y[i]; }
		re.f(j)=cpeds_integrate_1d(n,x,yy);
#pragma omp parallel for shared(y,yy,x,n)
		for (long i = 0; i < n; i++) {	yy[i]=-sin(twoPI*k*x[i])*y[i]; }
		im.f(j)=cpeds_integrate_1d(n,x,yy);
		j++;
		k+=dk;
#ifndef MSCS_QUIET
		printf("finished in: %lf %%              \r",(k-kmin)/(kmax-kmin)*100);
#endif
	}
#ifdef DEBUG
	getrusage(RUSAGE_SELF, &ru);
	endTime = ru.ru_utime;
	getrusage(RUSAGE_SELF, &ru);
	double tS = startTime.tv_sec*1000000 + (startTime.tv_usec);
	double tE = endTime.tv_sec*1000000  + (endTime.tv_usec);
	FILE* F=fopen("mscsFunction_powerSpectrumSFT1d.timing","w"); fprintf(F,"%lf\n",(tE-tS)*1e-6); fclose(F);
	printf("SFT exec time [s]: %lf\n",(tE-tS)*1e-6);
	//	debug_time2=clock();
	//	FILE* F=fopen("mscsFunction_powerSpectrumSFT1d.timing","w"); fprintf(F,"%lf\n",double(debug_time2-debug_time1)/CLOCKS_PER_SEC); fclose(F);
	//	printf("SFT exec time [s]: %lf\n",double(debug_time2-debug_time1)/CLOCKS_PER_SEC);
#endif
	
	delete [] x;
	delete [] y;
	delete [] yy;
	
	if (Re!=NULL) (*Re)=re;
	if (Im!=NULL) (*Im)=im;
	re.power(2.0);
	im.power(2.0);
	return re+im;	
}
/***************************************************************************************/
mscsFunction mscsFunction::FourierSeries(long N, double regridX, mscsFunction* A, mscsFunction* phi) {
	mscsFunction amp,phase;
	mscsFunction re,im,tmp;
	double T,a,p,T0;
	
	if (N>pointsCount()/2) { 
		msgs->warning("N in Fourier series is too big. Will use maximal possible value",Medium); 
		N=pointsCount()/2;
	}
	powerSpectrum(&re,&im,regridX,1);
//	re.save("re");
//	im.save("im");
	T0=1.0/re.getx(1);
	for (unsigned long i = 1; i <= N; i++) {
		T=1.0/re.getx(i);
		a=2.0*sqrt(pow(re.f(i),2)+pow(im.f(i),2));
		p=atan2(re.f(i),im.f(i));
		amp.newPoint(T,a);
		phase.newPoint(T,p);
//		printf("T=%lf A=%lf phi=%lf\n",T,A,phi);
	}

	if (A!=NULL) *A=amp;
	if (phi!=NULL) *phi=phase;

	tmp=*this;
	tmp=sqrt(pow(re.f(0),2)+pow(im.f(0),2));
	for (unsigned long k = 0; k < amp.pointsCount(); k++) {
		T=-amp.getx(k);
		p=phase.f(k); //-tmp.getx(0);
		a=amp.f(k);
		for (unsigned long i = 0; i < tmp.pointsCount(); i++) {
			tmp.f(i)+=a*sin(twoPI/T*(tmp.getx(i)-tmp.getx(0)) + (p));
		}
	}
	return tmp;
	
}
/***************************************************************************************/
double mscsFunction::expectationValue() {
	double *X=extractArguments();
	double *Y=extractValues();
	long N=pointsCount();
	for (long i = 0; i < N; i++) {		Y[i]*=X[i];	}
	double integral = cpeds_integrate_1d(pointsCount(),X,Y);
	delete [] X;
	delete [] Y;
	return integral;	
}

/***************************************************************************************/
void mscsFunction::normalize() {
	divide(integrate());
}

/***************************************************************************************/
mscsFunction& mscsFunction::calibrateByStDev(long windowSize, long windowStep) {
	long i=0;
	long N=pointsCount()-windowSize;
	double* y=extractArguments();
	double stdev;
	if (windowStep==-1) windowStep=windowSize;
	
	while (i<N) {
		stdev=cpeds_variance(y,windowSize,i);
		for (long j = 0; j < windowSize; j++) {			f(i+j)/=stdev;		}
		i+=windowStep;
	}
	
	return *this;
}

/***************************************************************************************/
mscsFunction& mscsFunction::quantize(long n, double vmin, double vmax) {
	double v;
	double delta=(vmax-vmin)/n;
	
	long N=pointsCount();
	for (long i = 0; i < N; i++) {
		v= round(Y(i)/delta)*delta;
		if (v<vmin) v=vmin; 
		else
			if (v>vmax) v=vmax;
		//		printf("%li: from %lE to %lE\n",i,Y(i),v);
		f(i)=v;
	}
	return *this;
}

// ****************************************************************************************************
mscsFunction& mscsFunction::mkSin(double from, double to, double dx, double T, double phi, double A) {
	double x=from;
	double phiT=phi*T;
	while (x<=to) {
		newPoint(x,A*sin(twoPI/T*(x+phiT)));
		x+=dx;
	}
	return *this;
}
// ****************************************************************************************************
mscsFunction& mscsFunction::mkSquareWave(double from, double to, double dx, double T, double phi) {
	mkSin(from,to,dx,T,phi);
	long N=pointsCount();
	for (long i = 0; i < N; i++) {
		if (f(i)>=0) f(i)=1; else f(i)=-1;
	}
	return *this;
}
// ****************************************************************************************************
mscsFunction& mscsFunction::mkConst(double from, double to, double dx, double v) {
	msgs->say("generating constant function from: "+msgs->toStr(from)+" to: "+msgs->toStr(to)+" with step: "+msgs->toStr(dx)+" and value: "+msgs->toStr(v),Low);
	double x=from;
	while (x<=to) {
		newPoint(x,v);
		x+=dx;
		// printf("%lE\n",x);
	}
	msgs->say("done",Low);
	return *this;
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkTopHat(double from, double to, double dx, double v1, double from2, double to2, double v2) {
	mkConst(from,to,dx,v1);
	msgs->say("generating top-hat function.",Low);
	for (long i = 0; i < pointsCount(); i++) {
		if (getx(i)>from2 and getx(i)<to2) {
			f(i)=v2;
		}
	}
	return *this;
}

// ****************************************************************************************************
mscsFunction& mscsFunction::mkGumbelDistr(double from, double to, double dx, double beta, double m, double tol) {
	double x=from;
	double v;
	double z;
	if (pointsCount()==0) {
		while (x<=to+tol) {
			z=(x-m)/beta;
			v=1.0/beta * exp(-(z+exp(-z)));
			newPoint(x,v);
			x+=dx;
		}
	}
	else {
		for (unsigned long i = 0; i < pointsCount(); i++) {
			x=getx(i);
			z=(x-m)/beta;
			v=1.0/beta * exp(-(z+exp(-z)));
			setf(i,v);
		}
	}
	
	return *this;	
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkWeibullDistr(double from, double to, double dx, double lambda, double k, double tol) {
	double x=from;
	double v;
	if (pointsCount()==0) {
		while (x<=to+tol) {
			v=k/lambda * pow(x/lambda,k-1) * exp(-pow(x/lambda,k));
			newPoint(x,v);
			x+=dx;
		}
	}
	else {
		for (unsigned long i = 0; i < pointsCount(); i++) {
			x=getx(i);
			v=k/lambda * pow(x/lambda,k-1) * exp(-pow(x/lambda,k));
			setf(i,v);
		}
	}
	
	return *this;		
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkGauss(double from, double to, double dx, double A, double m, double s, double B, double tol) {
	double x=from;
	double v;
	if (pointsCount()==0) {
		while (x<=to+tol) {
			v=A * exp(-(x-m)*(x-m)/(2*s*s))+B;
			newPoint(x,v);
			x+=dx;
		}
	}
	else {
		for (unsigned long i = 0; i < pointsCount(); i++) {
			x=getx(i);
			v=A * exp(-(x-m)*(x-m)/(2*s*s))+B;
			setf(i,v);
		}
	}
	
	return *this;
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkSkewGauss(double from, double to, double dx, double A, double m, double s,double alpha, double tol) {
	double x=from;
	double v,nsig;
	double sq2=1.41421356237310;
	
	if (pointsCount()==0) {
		while (x<=to+tol) {
			v=A * exp(-(x-m)*(x-m)/(2*s*s));
			nsig=(m-x)*alpha/s;
			newPoint(x,v*(1.0-gsl_sf_erf(nsig/sq2)));
			x+=dx;
		}
	}
	else {
		for (unsigned long i = 0; i < pointsCount(); i++) {
			x=getx(i);
			v=A * exp(-(x-m)*(x-m)/(2*s*s));
			nsig=(m-x)*alpha/s;
			setf(i,v*(1.0-gsl_sf_erf(nsig/sq2)));
		}
	}
	
	return *this;	
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkLog(double from, double to, double dx, double base, double tol) {
	if (from==0) from=dx;
	if (from<0) msgs->criticalError("mkLog: from cannot be negative",Top);
	
	mkLine(from,to,dx,1,0,tol);
	if (base==0) lnY();
	else logY(base);
	return *this;		
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkExponent(double from, double to, double dx, double base, double tol) {
	if (from==0) from=dx;
	if (from<0) msgs->criticalError("mkLog: from cannot be negative",Top);
	
	mkLine(from,to,dx,1,0,tol);
	if (base==0) exponent();
	else exponentY(base);
	return *this;		
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkLogSpace(double from, double to, long N, double base) {
	double f,t;
	if (base==0) base=exp(1.0);
	f=log10(from)/log10(base);
	t=log10(to)/log10(base);
	//	printf("from: %lE\n",from);
	//	printf("to: %lE\n",to);
	//	printf("f: %lE\n",f);
	//	printf("t: %lE\n",t);
	double dx=(t-f)/N;
	
	
	mkLine(f,t,dx,1,0,dx/10);
	exponentY(base);
	
	return *this;		
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkLine(double A, double B) {
	double x;
	for (unsigned long i = 0; i < pointsCount(); i++) {
		x=getx(i);
		setf(i,A * x + B);
	}
	return *this;		
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkLine(double from, double to, double dx, double A, double B, double tol) {
	if (pointsCount()==0) {
		double x=from;
		while (x<=to+tol) {
			newPoint(x,A * x + B);
			x+=dx;
		}
	}
	else {
		double x;
		for (unsigned long i = 0; i < pointsCount(); i++) {
			x=getx(i);
			setf(i,A * x + B);
		}
	}
	return *this;	
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkLine(double from, double to, double dx, double A, double B, double x0, double tol) {
	if (pointsCount()==0) {
		double x=from;
		while (x<=to+tol) {
			newPoint(x,A * (x-x0) + B);
			x+=dx;
		}
	}
	else {
		double x;
		for (unsigned long i = 0; i < pointsCount(); i++) {
			x=getx(i);
			setf(i,A * (x-x0) + B);
		}
	}
	return *this;		
}
// ****************************************************************************************************
mscsFunction& mscsFunction::mkGaussCDF(double from, double to, double dx, double A, double m, double s) {
	double x=from;
	while (x<=to) {
		if (pointsCount()>0) newPoint(x,A * exp(-(x-m)*(x-m)/(2*s*s)));
		else newPoint(x,_f.last().y()+A * exp(-(x-m)*(x-m)/(2*s*s)));
		x+=dx;
	}
	return *this;
}

/***************************************************************************************/
mscsFunction& mscsFunction::mkPowerLaw(double from, double to, double dk, double A, double k0, double ns) {
	if (pointsCount()==0) {
		double k=from;
		while (k<=to) {
			newPoint(k,A * pow((k/k0),ns));
			k+=dk;
		}
	}
	else {
		double k;
		for (unsigned long i = 0; i < pointsCount(); i++) {
			k=getx(i);
			setf(i,A * pow((k/k0),ns));
		}
	}
	return *this;	
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkPowerLaw(double from, double to, double dk, double A, double k0, double k1, double ns, double B) {
	if (pointsCount()==0) {
		double k=from;
		while (k<=to) {
			newPoint(k,A * pow(((k-k0)/k1),ns)+B);
			k+=dk;
		}
	}
	else {
		double k;
		for (unsigned long i = 0; i < pointsCount(); i++) {
			k=getx(i);
			setf(i,A * pow(((k-k0)/k1),ns)+B);
		}
	}
	return *this;		
}

/***************************************************************************************/
mscsFunction& mscsFunction::mkPowerLaw(double A, double k0, double ns) {
	long n=pointsCount();
	for (long i = 0; i < n; i++) {
		f(i)=A * pow((getx(i)/k0),ns);
	}
	return *this;	
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkBetaModel(double from, double to, double dr, double n0, double rc, double beta) {
	double alpha=-1.5*beta;
	if (pointsCount()==0) {
		double r=from;
		while (r<=to) {
			newPoint(r, n0 * pow(1.0+(r/rc)*(r/rc),alpha));
			r+=dr;
		}
	}
	else {
		double r;
		for (unsigned long i = 0; i < pointsCount(); i++) {
			r=getx(i);
			setf(i,n0 * pow(1.0+(r/rc)*(r/rc),alpha));
		}
	}
	return *this;		
}


/***************************************************************************************/
mscsFunction& mscsFunction::mkPressureProfile(double from, double to, double dz, double P0, double T0, double z0, double Lr, double g0, double mu) {
	double z=from;
	while (z<=to) {
		newPoint(z, cpeds_pressure_vs_altitude(z*1000,P0,T0,z0*1000,Lr/1000,g0,mu));
		z+=dz;
	}
	return *this;		
	
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkPlank(double from, double to, double dnu, double T) {
	double nu=from;
	double c3=CPEDS_c*CPEDS_c*CPEDS_c;
	while (nu<=to) {
		newPoint(nu, 8*PI*nu*nu/c3 * CPEDS_h*nu / (exp(CPEDS_h*nu/(CPEDS_kB*T))-1.0) );
		nu+=dnu;
	}
	return *this;		
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkPolynomial(double from, double to, double dx, int order, double a0, ...) {
	cpedsList<double> a;
	va_list args;
	va_start(args,a0);
	a.append(a0);
	for (long i = 1; i <= order; i++) {
		a.append(va_arg(args,double));
	}
	va_end(args);
	
	mkPolynomial(from,to,dx,a);
	return *this;
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkPolynomial(double from, double to, double dx, cpedsList<double> a) {
//	while (x<=to+tol) {
//		v=A * exp(-(x-m)*(x-m)/(2*s*s));
//		nsig=(m-x)*alpha/s;
//		newPoint(x,v*(1.0-gsl_sf_erf(nsig/sq2)));
//		x+=dx;
//	}

	if (a.size()>0) {
		mkConst(from,to,dx,a[0]);
		long N=pointsCount();
		double x;
		for (long i = 1; i <a.size(); i++) {
			for (long j = 0; j < N; ++j) {
				x=getx(j);
				f(j)+=a[i]*pow(x,double(i));
			}
		}
		
	}
	return *this;	
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkPolynomial(cpedsList<double> a) {
	long N=pointsCount();
	double x;
	if (a.size()>0) {
		setf(a[0]);
		for (long i = 1; i <a.size(); i++) {
			for (long j = 0; j < N; ++j) {
				x=getx(j);
				f(j)+=a[i]*pow(x,double(i));
			}
		}			
	}
	return *this;	
}

/***************************************************************************************/
mscsFunction& mscsFunction::mkBoltzmann(double from, double to, double dE, double T, double mu) {
	double E=from;
	while (E<=to) {
		newPoint(E, exp(-(E-mu)/(CPEDS_kB*T)) );
		E+=dE;
	}
	return *this;		
}

/***************************************************************************************/
mscsFunction& mscsFunction::mkGaussianNoise(long N, double m, double s, long seed, cpedsRNG* rng, bool noisex, bool noisey) {
	cpedsRNG *rns;
	//
	// setup the random number generator
	//
	if (rng!=NULL) {
		rns=rng;
	}
	else {
		if (_rns==NULL) {
			_rns = new cpedsRNG("gaussian_circle");
			rns=_rns;
			rns->seed(seed);
		}
		else {
			rns=_rns;		
		}
	}	
	rns->setMeanStd(m,s);
	
	double *a=NULL;
	
	if (noisex and !noisey) {
		a= rns->getRN(N);
		for (long i=0;i<N;i++) { newPoint(a[i],0); }
	}
	if (noisey and !noisex) {
		a= rns->getRN(N);
		for (long i=0;i<N;i++) { newPoint(i, a[i]); }		
	}
	if (noisex and noisey) {
		long M=2*N;
		a= rns->getRN(M);
		for (long i=0;i<M;i++) { newPoint(a[i], a[i+N]); }		
	}
	
	if (a!=NULL) delete [] a;
	
	return *this;
}

/***************************************************************************************/
mscsFunction& mscsFunction::mkChisqNoise(long N, double m, double s, long seed, cpedsRNG* rng, bool noisex, bool noisey) {
	
	mkGaussianNoise(N,m,s,seed,rng,noisex,noisey);
	power(2.0);
	
	return *this;
}

/***************************************************************************************/
mscsFunction& mscsFunction::mkPowerLawNoise(double A, double alpha, double k0, double kmin, double kmax, long N, long seed, mscsFunction* powerSpectrum, cpedsRNG* rng, bool realSpace) {
	cpedsRNG *rns;
	if (powerSpectrum!=NULL) powerSpectrum->clearFunction();
	
	//
	// setup the random number generator
	//
	
	if (rng!=NULL) {
		rns=rng;
	}
	else {
		if (_rns==NULL) {
			_rns = new cpedsRNG("gaussian_circle");
			rns=_rns;
			rns->seed(seed);
		}
		else {
			rns=_rns;		
		}
	}
	
	//
	// generate gaussian noise
	//
	long M;
	
	if (N%2==1) M=(N-1)/2+1; else M=N/2+1;
	
	double *a= rns->getRN(M);
	double *b= rns->getRN(M);
	
	//
	// form the right power spectrum in Fourier space (positive frequencies only)
	//
	double Pk,sqrPj;
	//	double Lmin=1/kmax; // Jul 4, 2011, 2:05:26 PM commented out
	//	double Lmax=1/kmin; // Jul 4, 2011, 2:05:26 PM commented out
	double T; // this is the total sampling period T = N dt = N/(2 kmax)
	//	double Lj;
	//	L=fabs(Lmax-Lmin); // Jul 4, 2011, 2:05:45 PM - changed
	
	long j; // iterates k_j; j=0..M-1
	double k;
	double dk;
	
	// changes on  Jul 4, 2011, 2:07:40 PM --  begin
	if (realSpace) { 
		dk=double(2)*kmax/double(N);
		T=double(N)/(2.0 * kmax);
		kmin=0;
		//		printf("N: %li\n",N);
		//		printf("kmax: %lE\n",kmax);
		//		printf("dk: %lE\n",dk);
		//		printf("T: %lE\n",T);
	}
	else { 
		dk=(kmax-kmin)/double(M);
		T=double(M)/(kmax-kmin); // Jul 4, 2011, 2:07:40 PM
	}
	// changes on  Jul 4, 2011, 2:07:40 PM --  end
	
	for (j=1;j<M;j++) {
		k=double(j)*dk+kmin;
		Pk=A*powf( k / k0,alpha); 
		//		printf("k: %lE, dk: %lE, pk: %lE\n",k, dk, Pk);
		sqrPj=sqrt(Pk/2.0);
		a[j]*=sqrPj;
		b[j]*=sqrPj;
	}
	
	// set power for k=0;
	a[0]=0;
	b[0]=0;
	
	//
	//	store power spectrum realization if requested
	//
	if (powerSpectrum!=NULL) {		
		//		for (j=1;j<M;j++) {	powerSpectrum->newPoint(double(j)/L,a[j]*a[j]+b[j]*b[j]);	}
		k=kmin;
		for (j=0;j<M;j++) {	powerSpectrum->newPoint(k,a[j]*a[j]+b[j]*b[j]); k+=dk;	}
	}
	else {
		if (realSpace==false) {
			k=kmin;
			for (j=0;j<M;j++) {	newPoint(k,a[j]*a[j]+b[j]*b[j]); k+=dk;	}			
		}
	}
	
	if (realSpace) {
		//
		// do fftw
		//
		long i; // iterates x_i, i=0..N-1 
		//		for (i=0;i<M;i++) { printf("a: %lE, b: %lE\n",a[i],b[i]); }
		double* out=fft_1D_c2r(a,b,M,N,true);
		
		//
		// make a new function
		//
		clearFunction();
		double dt=T/double(N);
		for (i=0;i<N;i++) { newPoint(dt*i, out[i]); }
		delete [] out;
		//		fftw_free(out);
		
	}
	else {
		delete [] a;
		delete [] b;
	}
	return *this;
	
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkPhaseNoise(long N, long seed, cpedsRNG* rng) {
	cpedsRNG *rns;
	//	if (powerSpectrum!=NULL) powerSpectrum->clearFunction();
	
	//
	// setup the random number generator
	//
	
	if (rng!=NULL) {
		rns=rng;
	}
	else {
		if (_rns==NULL) {
			_rns = new cpedsRNG("gaussian_circle");
			rns=_rns;
			rns->seed(seed);
		}
		else {
			rns=_rns;		
		}
	}
	
	//
	// generate uniform phase noise
	//
	long M;
	
	if (N%2==1) M=(N-1)/2+1; else M=N/2+1;
	
	double *a= rns->getRN(M);
	double *b= rns->getRN(M);
	a[0]=0; // set mean
	b[0]=0; // reality condition
	if (N%2==0) b[M-1]=0; // reality condition
	
	
	long j;
	double p;
	for (j=0;j<M;j++) {	
		p=sqrt(a[j]*a[j]+b[j]*b[j]);
		if (p==0) {
			a[j]=p;
			b[j]=p;			
		}
		else {
			a[j]/=p;
			b[j]/=p;
		}
		//		printf("a,b: %lE, %lE\n",a[j],b[j]);
	}
	
	//	if (powerSpectrum!=NULL) { 
	//		for (j=0;j<M;j++) {	
	//			powerSpectrum->newPoint(double(j)/N,a[j]*a[j]+b[j]*b[j]);
	//		}
	
	double* out=fft_1D_c2r(a,b,M,N,true);
	
	clearFunction();
	importFunction(NULL,out,N);
	delete [] out;
	return *this;
}

/***************************************************************************************/
mscsFunction& mscsFunction::mkPowerLawNoise(double A, double alpha, double k0, long seed, mscsFunction* powerSpectrum, cpedsRNG* rng, bool realSpace) {
	if (powerSpectrum!=NULL) powerSpectrum->clearFunction();
	
	cpedsRNG *rns;
	//
	// setup the random number generator
	//
	
	if (rng!=NULL) {
		rns=rng;
	}
	else {
		if (_rns==NULL) {
			_rns = new cpedsRNG("gaussian_circle");
			rns=_rns;
			rns->seed(seed);
		}
		else {
			rns=_rns;		
		}
	}
	
	//
	// generate gaussian noise
	//
	long M;
	M=pointsCount();
	
	
	double *a= rns->getRN(M);
	double *b= rns->getRN(M);
	
	//
	// form the right power spectrum in Fourier space (positive frequencies only)
	//
	double Pk,sqrPj;
	long j; // iterates k_j; j=0..M-1
	double k;
	for (j=0;j<M;j++) {
		k=getX(j);
		Pk=A*powf( k / k0,alpha);
		sqrPj=sqrt(Pk/2.0);
		a[j]*=sqrPj;
		b[j]*=sqrPj;
	}
	
	//
	//	store power spectrum realization if requested
	//
	if (powerSpectrum!=NULL) {		
		for (j=0;j<M;j++) {	powerSpectrum->newPoint(getX(j),a[j]*a[j]+b[j]*b[j]); }
	}
	else {
		if (realSpace==false) {
			for (j=0;j<M;j++) {	f(j)=a[j]*a[j]+b[j]*b[j]; }
		}
	}
	
	delete [] a;
	delete [] b;
	
	
	return *this;
}

/***************************************************************************************/
mscsFunction& mscsFunction::mkHistogram(cpedsList<double>& data, long Nbins, long binAlign) {
	if (data.size()>0) {
		clearFunction();
		double *d=data.toCarray();
		long Nbin=Nbins;
		double * bins=new double[Nbin];
		double *counts;
		counts=cpeds_bin_data(data.size(),d,Nbin,bins,binAlign);
		
		importFunction(bins,counts,Nbin);
		
		delete [] counts;
		delete [] bins;
		delete [] d;		
	}
	return *this;
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkHistogram(cpedsList<double>& data, double fromX, double toX, long Nbins, long binAlign) {
	if (data.size()>0) {
		clearFunction();
		double *d=data.toCarray();
		long j=0;
		for (long i = 0; i < data.size(); i++) {
			if (data[i]>=fromX and data[i]<=toX) { d[j]=data[i]; j++; }
		}
		
		long Nbin=Nbins;
		double * bins=new double[Nbin];
		double *counts;
		counts=cpeds_bin_data(j,d,Nbin,bins,binAlign);
		
		importFunction(bins,counts,Nbin);
		
		delete [] counts;
		delete [] bins;
		delete [] d;		
	}
	return *this;
}


// ****************************************************************************************************
mscsFunction& mscsFunction::deleteOutsideOf(double from, double to) {
	for (long i=0;i<_f.count();i++) { if (X(i)<from or X(i)>to) { deletePoint(i); i--;} }
	return *this;
}
// ****************************************************************************************************
mscsFunction& mscsFunction::deleteRange(double from, double to) {
	for (long i=0;i<_f.count();i++) { if (X(i)>=from and X(i)<=to) { deletePoint(i); i--; }  } //printf("deleteing frequency: %lf\n",X(i)); 
	return *this;
}
/***************************************************************************************/
mscsFunction& mscsFunction::deleteRanges(mscsFunction fromto) {
	for (long i = 0; i < fromto.pointsCount(); i++) {
		deleteRange(fromto.getx(i),fromto.f(i));
	}
	return *this;
}
/***************************************************************************************/
mscsFunction& mscsFunction::shift(long n) {
	long N=_f.count();
	
	if (n!=0 && n!=N) {
		
		if (n>0) {
			for (long j = 1; j <= n; j++) {				_f.prepend(_f[_f.count()-j]);			}
			for (long j = 0; j < n; j++) {				_f.takeLast();			}
		}
		else {
			//			printf("n=%li\n",n);
			n=-n;
			for (long j = 0; j < n; j++) {				_f.append(_f[j]);	}
			for (long j = 0; j < n; j++) {				_f.takeFirst();	 }
		}
		
		//		mscsFunction tmp(*this);
		//		long ip;
		//		for (long i=0;i<N;i++) { 
		//			ip=(i+n) % N;
		//			f(ip)=tmp.f(i);		
		//		}
		
	}
	return *this;	
}
mscsFunction& mscsFunction::shiftYwrtX(long n) {
	mscsFunction tmp=(*this);
	tmp.shift(n);
	this->setY(tmp.getValues());
	return *this;	
}

/***************************************************************************************/
double* mscsFunction::fft_1D_c2r(double *inRe, double* inIm, long M, long N, bool deleteIn) {
#ifdef DEBUG_USING_FFTW2AND3
	fftw_complex* in=new fftw_complex[(M)*sizeof(fftw_complex)]; 
#else
	fftw_complex* in=(fftw_complex*)fftw_malloc((M)*sizeof(fftw_complex)); 
#endif
	for (long j=0;j<M;j++) { in[j][0]=inRe[j]; in[j][1]=inIm[j]; }
	if (deleteIn) {	delete [] inRe;	delete [] inIm; }
	
	//	double* out=(double*)fftw_malloc((N)*sizeof(double)); 
	double* out=new double[N];
	fftw_plan p;
	p=fftw_plan_dft_c2r_1d(N, in,out, FFTW_ESTIMATE); 
	fftw_execute(p);
	fftw_destroy_plan(p);
#ifdef DEBUG_USING_FFTW2AND3
	delete [] in;
#else
	fftw_free(in);
#endif
	fftw_cleanup();
	//	for (long j=0;j<N;j++) { out[j]/=N; } // BLcomment (Nov 14, 2011, 12:59:53 AM): this should be in r2c function
	return out;
}
/***************************************************************************************/
void mscsFunction::fft_1D_r2c(double *in, long size, double *re, double *im) {
	long N=size;
	fftw_plan p;
	long i;
#ifdef DEBUG_USING_FFTW2AND3
	fftw_complex* t=new fftw_complex[N*sizeof(fftw_complex)];
#else
	fftw_complex* t=(fftw_complex*)fftw_malloc(N*sizeof(fftw_complex)); 
#endif
	int dir=FFTW_FORWARD;
	for (long j = 0; j < N; j++) {
		t[j][0]=in[j];
		t[j][1]=0.0;
	}
	p=fftw_plan_dft_1d(N, t,t, dir, FFTW_ESTIMATE); 
	fftw_execute(p);
	
	fftw_destroy_plan(p);
	long M=N/2+1;
	re=new double[M];
	im=new double[M];
	for (long j = 0; j < M; j++) {
		//		re[j]=t[j][0];
		//		im[j]=t[j][1];
		re[j]=t[j][0]/N;
		im[j]=t[j][1]/N;
	}
#ifdef DEBUG_USING_FFTW2AND3
	delete [] t;
#else
	fftw_free(t);
#endif
	fftw_cleanup();
}
/***************************************************************************************/
mscsFunction& mscsFunction::X2Y(bool x2y) {
	long N=pointsCount();
	if (x2y)
		for (long i=0;i<N;i++) { f(i)=getx(i); }
	else
		for (long i=0;i<N;i++) { X(i)=Y(i); }
	
	return *this;
}
/***************************************************************************************/
mscsFunction& mscsFunction::remove(const mscsFunction& p) {
	long j;
	for (long i = 0; i < p.pointsCount(); i++) {
		j=0;
		while (j<pointsCount()) {
			if (p.getQPoint(i)==_f[i]) deletePoint(j);
			else
				j++;
		}
		msgs->say("removing point: "+msgs->toStr(i)+"from: "+msgs->toStr(p.pointsCount()),Low);
	}
	return *this;
}

/***************************************************************************************/
mscsFunction& mscsFunction::removeValue(double v) {
	long N=pointsCount()-1;
	for (long i=N;i>=0;i--) { 
		if (f(i)==v) {
			deletePoint(i); 
		}
	}
	return *this;	
}
/***************************************************************************************/
mscsFunction& mscsFunction::removeSmaller(double v) {
//	long N=pointsCount()-1;
//	for (long i=N;i>=0;i--) { 
//		if (f(i)<v) {
//			deletePoint(i); 
//		}
//	}
	
	mscsFunction tmp;
	long N=pointsCount();
	for (long i=0;i<N;i++) { 
		if (f(i)>=v) {
			tmp.newPoint(getx(i),f(i));
		}
	}
	
	*this=tmp;
	
	return *this;		
}
/***************************************************************************************/
mscsFunction& mscsFunction::removeLarger(double v) {
	
	mscsFunction tmp;
	long N=pointsCount();
	for (long i=0;i<N;i++) { 
		if (f(i)<=v) {
			tmp.newPoint(getx(i),f(i));
		}
	}
	
	*this=tmp;
	
	return *this;			
}
/***************************************************************************************/
mscsFunction& mscsFunction::removePoints(cpedsList<long> l) {
	mscsFunction tmpf("tmpf",getVerbosityLevel());
	long* idx=l.toCarray();
	cpeds_sort_data(l.size(),idx,21);
	for (long i = 0; i < l.size(); i++) {
		deletePoint(idx[i]);
	}
	delete [] idx;
	return *this;
}
/***************************************************************************************/
mscsFunction& mscsFunction::removeNans() {
	for (long i = 0; i < pointsCount(); i++) {
		if (cpeds_isnan(getx(i))) { deletePoint(i); i--; }
		else {
			if (cpeds_isnan(f(i))) { deletePoint(i); i--; }
		}
	}
	return *this;	
}
/***************************************************************************************/
mscsFunction& mscsFunction::mkCDF(bool zoomIn, double acc) {
	long N=pointsCount();
	long Nleo=N-1;
	for (long i=1;i<N;i++) { 
		f(i)+=f(i-1); 
	}
	double norm=f(N-1);
	for (long i=0;i<N;i++) { f(i)/=norm; }		
	
	if (zoomIn) {
		for (long i=1;i<N;i++) { 
			if (abs(f(i)-1.0)<acc and i<Nleo) {
				f(i)=double(1.0);
				cut(N-i-1,i+1);
				msgs->warning("The CDF calculation exceeds the double accuracy. Will zoom-in into interesting regions, neglecting non-computable one.",High);
				mscsFunction tmp("tmp",Zero);
				tmp=interpolate(getX(0),getX(i),double(getX(i)-getX(0))/N,false,"cspline",false);
				*this=tmp;
				i=N;
			}
		}
	}
	return *this;
	
}
/***************************************************************************************/
//mscsFunction& mscsFunction::mkcCDF(bool zoomIn, double acc) {
//	mscsFunction tmp=(*this);
//	long N=pointsCount();
//	long Nleo=N-1;
//	for (long i=1;i<N;i++) { 
//		f(i)+=f(i-1); 
//	}
//	double norm=f(N-1);
//	for (long i=0;i<N;i++) { f(i)/=norm; }		
//
//	if (zoomIn) {
//		for (long i=1;i<N;i++) { 
//			if (abs(f(i)-1.0)<acc and i<Nleo) {
//				f(i)=double(1.0);
//				cut(N-i-1,i+1);
//				msgs->warning("The CDF calculation exceeds the double accuracy. Will zoom-in into interesting regions, neglecting non-computable one.",High);
//				mscsFunction tmp("tmp",Zero);
//				tmp=interpolate(getX(0),getX(i),double(getX(i)-getX(0))/N,false,"cspline",false);
//				*this=tmp;
//				i=N;
//			}
//		}
//	}
//	return *this;
//		
//}
/***************************************************************************************/
mscsFunction mscsFunction::integrateX_Xmax() {
	long N=pointsCount();
	mscsFunction intxxmax;
	checkRanges();
	double xmax=getMaxArg();
	printf("xmax: %lE\n",xmax);
	//	exit(0);
	double x;
	for (long i=0;i<N;i++) { 
		x=getx(i);
		intxxmax.newPoint(x,integrate(x,xmax));
	}
	return intxxmax;
}
/***************************************************************************************/
mscsFunction mscsFunction::convolve_window(mscsFunction w, bool wreal) {
	//	double *in=extractValues();
	//	long Nout=pointsCount()/2+1;
	//	double *re;
	//	double *im;
	//	fft_1D_r2c(in,pointsCount(),re,im);
	mscsFunction re,im,wre,wim;
	mscsFunction wf;	
	if (wreal) { 
		wf=w.powerSpectrum(&wre,&wim,0); 
	} 
	else wf=w;
	
	//	wf.power(0.5);
	wf.sqroot();
	
	
	mscsFunction p=powerSpectrum(&re,&im,0);
	
	
	wf.interpolate(p.getMinArg(),p.getMaxArg(),(p.getMaxArg()-p.getMinArg())/(p.pointsCount()-1),true,"linear",true);
	
	//	wf.save("window_interpolated_linear.p");
	//	wf.interpolate(p.getMinArg(),p.getMaxArg(),(p.getMaxArg()-p.getMinArg())/p.pointsCount(),true,"cspline",true);
	//	wf.save("window_interpolated_cspline.p");
	//	delete [] xint;
	re*=wf;
	im*=wf;
	
	wf.clearFunction();
	double *x=extractArguments();
	double *conv=fft_1D_c2r(re.extractValues(),im.extractValues(),p.pointsCount(),pointsCount(),true);
	mscsFunction convolved("convolved",x,conv,pointsCount(),getVerbosityLevel());
	//		mscsFunction convolved("convolved",x,x,pointsCount(),getVerbosityLevel());
	delete [] x;
	delete [] conv;
	return convolved;
	
	//	return *this;
}
/***************************************************************************************/
double mscsFunction::correlationCoefficient(mscsFunction& f2) {
	double corr=covariance(f2)/(stdev() * f2.stdev());	
	return corr;
}

/***************************************************************************************/
double mscsFunction::covariance(mscsFunction& f2) {
	long n1=pointsCount();
	long n2=f2.pointsCount();
	long n=cpeds_get_min(n1,n2);
	if (n1!=n2) { msgs->warning("covariance>> the functions are of different length. Will calculate only for the length of the shorter vector.",High);  }
	double* d1=extractValues(n);
	double* d2=f2.extractValues(n);
	
	double cov=cpeds_covariance(d1,d2,n);
	delete [] d1;
	delete [] d2;
	return cov;
}
/***************************************************************************************/
mscsFunction mscsFunction::correlationCoefficientFunction(mscsFunction& f2, long step, long N) {
	long n1=pointsCount();
	long n2=f2.pointsCount();
	long n=cpeds_get_min(n1,n2);
	if (n1!=n2) { msgs->warning("correlationCoefficientFunction>> the functions are of different length. Will calculate only for the length of the shorter vector.",High);  }
	if (cpeds_get_min(n1,n2)<2) { msgs->error("correlationCoefficientFunction>> the function length must be larger than 1 point.",High);  }
	double* d1=extractValues(n);
	double* d2=f2.extractValues(n);
	double s1=stdev();
	double s2=f2.stdev();
	double s=s1*s2;
	double dx=(X(1)-X(0))*step;
	if (N==0) N=n;
	mscsFunction corr("correlation function");
	//	long i=0;
	//	while (i<N) {
	for (long i = 0; i < N; i+=step) {
//		printf("i:%li, x:%lf cov: %lf\n ",i,i*dx,cpeds_covariance(d1,d2,n)/s);
		corr.newPoint(i*dx,cpeds_covariance(d1,d2,n)/s);
		// shift array in f2
//		printf("i:%li, x:%lf\n ",i,i*dx);
		cpeds_shift_array(d2,n,step,true);
		//		i+=step;
	}

	delete [] d1;
	delete [] d2;
	return corr;
}
/***************************************************************************************/
bool mscsFunction::isPositive() const {
	for (long i = 0; i < pointsCount(); i++) {
		if (f(i)<0) return false;
	}
	return true;	
}
/***************************************************************************************/
mscsFunction& mscsFunction::average_sameArgs(double acc, mscsFunction* stdev) {
	if (pointsCount()>1) {
		
		bool sameArg=true;
		double arg=getx(0);
		long toIdx=0;
		long fromIdx=0;
		mscsFunction tmpf("tmpf",getVerbosityLevel());
		
		if (stdev!=NULL) stdev->clearFunction();
		while (fromIdx<pointsCount()) {
			arg=getx(fromIdx);
			toIdx=fromIdx;
			sameArg=true;
			while (sameArg) {
				toIdx++;
				if (toIdx<pointsCount()) {
					if (fabs(getx(toIdx)-arg)>acc ) sameArg=false;				
				}
				else sameArg=false;
			}
			toIdx--;
			//		printf("from: %li, to: %li\n",fromIdx,toIdx);
			if (toIdx!=fromIdx) {
				// average
				tmpf=cut(toIdx-fromIdx+1, fromIdx);
				insertPoint(tmpf.getx(0),tmpf.meanf(),fromIdx,1);
				if (stdev!=NULL) stdev->newPoint(tmpf.getx(0),tmpf.stdev());
			}
			fromIdx++;
			//		printf("   from: %li, pointsCount:%li\n",fromIdx,pointsCount());
		}
	}
	
	return *this;
}
/***************************************************************************************/
mscsFunction& mscsFunction::convertUnixTimeToJD_x() {
	for (long i = 0; i < pointsCount(); i++) {
		setarg(i,cpeds_timeSec_to_julian_time(getx(i)));
	}
	return *this;	
}
/***************************************************************************************/
vector<double> mscsFunction::findRoot(double period) {
	vector<double> r;
	double x1,x2,y1,y2;

	if (period>0) {
		newPoint(getx(0)+period,getY(0));
	}
	
	x1=getX(0);
	y1=getY(0);
	if (pointsCount()<2) return r;
	for (long i = 1; i < pointsCount(); i++) {
		x2=getX(i);
		y2=getY(i);
		
		if (y1<0 and y2>=0) {
			// found root
			mscsFunction tmp;
			tmp.newPoint(y1,x1);
			tmp.newPoint(y2,x2);
			r.push_back(tmp.finter(0,"linear"));
		}
		else {
			
			if (y1>0 and y2<=0) {
				// found root
				mscsFunction tmp;
				tmp.newPoint(y1,x1);
				tmp.newPoint(y2,x2);
				r.push_back(tmp.finter(0,"linear"));
			}
		}
		
		x1=getX(i);
		y1=getY(i);
		
	}

	if (period>0) {
		for (unsigned long i = 0; i < r.size(); i++) {
			if (r[i]>x1+period) r[i]-=period;
		}
		deletePoint();
	}

	return r;
}
