/*!
 * \file MClink.cpp
 *
 *  Created on: Nov 19, 2010
 *      Author: blew
 */

#include "MClink.h"
#include "cpeds-templates.h"
#include <assert.h>

/***************************************************************************************/
MClink::MClink() {
	_theta=NULL;
	_Nparam=0;
	_likelihood=-1;
	_chisq=-1;
	_accepted=false;
	_idx=0;
}
/***************************************************************************************/
MClink::MClink(long Nparams) {
	_theta=new double[Nparams];	
	for (long i = 0; i < Nparams; i++) {		_theta[i]=0;	}
	_Nparam=Nparams;
	_likelihood=-1;
	_chisq=-1;
	_accepted=false;
	_idx=0;
}
/***************************************************************************************/
MClink::MClink(const MClink& parent) {
//	*this=parent;
	
	_Nparam=parent.dims();
	_likelihood=parent.L();
	_chisq=parent.chisq();
	_theta=cpeds_copy_array(parent._theta,_Nparam);
	_accepted=parent.isAccepted();
	_idx=parent.getIdx();
}
/***************************************************************************************/
MClink::~MClink() {
	clear();
}
/***************************************************************************************/
void MClink::set(int n, ...) {
	if (n==_Nparam) {
		va_list Numbers;
		va_start(Numbers, n); 
		for(int i = 0; i < n; ++i )
			_theta[i]=va_arg(Numbers, double);
		va_end(Numbers);

	}
}
/***************************************************************************************/
void MClink::set(int n, double *t, bool deleteInside) {
	setNpar(n);
	for (long i = 0; i < n; i++) {		_theta[i]=t[i];	}
	if (deleteInside) delete [] t;
}

/***************************************************************************************/
void MClink::setNpar(long n) {
	if (n!=_Nparam) { 
		clear();
		if (n!=0) {	
			_theta=new double[n];
			_Nparam=n;
		}
		else { printf("MClink ERROR: setNpar n cannot be 0. I better stop now.\n"); exit(0); }
	}
}

/***************************************************************************************/
const MClink& MClink::operator=(const MClink& rhs) {
//	printf("MClink::operator= in>>> _theta: %li, rhs.theta: %li, dims: %li, rhs.dims: %li\n",_theta,rhs._theta,dims(), rhs.dims());
	if (this!=&rhs) {
		if (dims()!=rhs.dims()) {
			clear();
			if (rhs.dims()>0) _theta=rhs.getParams();
			_Nparam=rhs.dims();
		}
		else {
			for (long i = 0; i < _Nparam; i++) {				_theta[i]=rhs[i];			}
		}
		setL(rhs.L());
		setChisq(rhs.chisq());
		setIdx(rhs.getIdx());
		setAccepted(rhs.isAccepted());
	}
//	printf("MClink::operator= out>>> _theta: %li, rhs.theta: %li, dims: %li, rhs.dims: %li\n",_theta,rhs._theta,dims(), rhs.dims());
	return *this;
}
/* ******************************************************************************************** */
const MClink& MClink::operator=(const double v) {
	for (long i = 0; i < _Nparam; i++) {				_theta[i]=v;			}
	return *this;
}
/* ******************************************************************************************** */
bool MClink::operator==(const double v) {
	for (long i = 0; i < _Nparam; i++) {
		if (_theta[i]!=v) return false;	}
	return true;
}
/* ******************************************************************************************** */
MClink MClink::operator-(const MClink& rhs) {
	MClink mc(*this);
	return mc-=rhs;
}
/* ******************************************************************************************** */
MClink MClink::operator*(const double v) {
	MClink mc(*this);
	return mc*=v;	
}
/* ******************************************************************************************** */
const MClink& MClink::operator*=(const double v) {
	for (long i=0;i<_Nparam;i++) { setParam(i,getParam(i)*v); }
	return *this;
}
/* ******************************************************************************************** */
MClink& MClink::operator-=(const MClink& rhs) {
	assert(rhs.dims()==dims());
	for (long i=0;i<_Nparam;i++) { setParam(i,getParam(i)-rhs.getParam(i)); }
	setChisq(chisq()-rhs.chisq());
	setIdx(getIdx()-rhs.getIdx());
	_likelihood-=rhs.L();
	return *this;
}
/***************************************************************************************/
const double& MClink::operator[](const long i) const {
//	printf("operator[]:: theta: %li, i: %li, dims: %li\n",_theta,i,dims());
	return _theta[i];
}
/***************************************************************************************/
double& MClink::operator[](const long i) {
	return _theta[i];
}
/***************************************************************************************/
void MClink::printParams() const {
//	printf("\n------\nMClink\n");
	for (long i = 0; i < _Nparam; i++) {
		printf("param %li:  %lE\n",i,_theta[i]);
	}
	if (_chisq!=-1) {
		printf("chisq:  %lE\n",_chisq);
	}
	if (_likelihood!=-1) {
		printf("likelihood:  %lE\n",_likelihood);
	}
	printf("idx:  %li\n",_idx);
	if (isAccepted()) printf("The link is accepted\n"); else printf("The link is rejected\n");
//	printf("------\n");
	
}
/***************************************************************************************/
void MClink::save(string fname) {
	FILE* f=fopen(fname.c_str(),"w");
	for (long i = 0; i < dims(); i++) {
		fprintf(f,"%.10lE ", getParam(i));
	}
	fprintf(f,"%.10lE ", chisq());
	fprintf(f,"%.10lE ", L());
	fprintf(f,"%lf ", _idx);
	long acc;
	if (isAccepted()) acc=1; else acc=0;
	fprintf(f,"%li ", acc);
	fprintf(f,"\n");
	fclose(f);
}
/* ******************************************************************************************** */
MClink& MClink::load(ifstream& ifs, long Nparam) {
	double data[Nparam];
	double X2,L,idx,accepted;
	
	for (long i = 0; i < Nparam; i++) {
		ifs >> data[i];
//		cout << data[i] << " ";
	}	
	
	ifs >> X2;
	ifs >> L;
	ifs >> idx;
	ifs >> accepted;
	
//	cout << X2 << " " << L << " " << idx << " " << accepted << "\n";
	
	set(dims(),data,false);
	setChisq(X2);
	setL(L);
	setAccepted(accepted);
	return *this;
}
/* ******************************************************************************************** */
MClink& MClink::load(string fname) {
	cpeds_queue<long>* q;
	q=cpeds_get_txt_file_cols_rows(fname.c_str());
	if (q==0) { return *this; }	

	long Nrows=q->get_size();
	long Ncols=q->getq(0);
	long Nparam=Ncols-4;
		
	ifstream ifs;
	ifs.open(fname.c_str());
	
	setNpar(Nparam);
	load(ifs,Nparam);
	ifs.close();
	return *this;
}
/***************************************************************************************/
cpedsList<double> MClink::getParameters() const {
	cpedsList<double> p;
	p.fromCarray(_theta,dims());
	return p;
}
/* ******************************************************************************************** */
void MClink::setParam(int n, double t) {
	assert(n<dims());
	_theta[n]=t;
}
/* ******************************************************************************************** */
void MClink::printLink() const {
	cout << "dims: " << dims() << "\n";
	printParams();
}
/* ******************************************************************************************** */
ostream& operator<< (ostream& output, const MClink& l ) {
	for (long i=0;i<l.dims();i++) {
		output << l.getParam(i) << " ";
	}
	output << l.chisq();
	output << std::endl;
	return output;
}
/* ******************************************************************************************** */
/* ******************************************************************************************** */
MClink& operator*(double v,MClink& l) {
	for (long i=0;i<l.dims();i++) { l.setParam(i,l.getParam(i)*v); }
	return l;	

}

void MClink::printParamsLine(string info) const {
	std::cout << info;
	for (long i=0;i<dims();i++) {
		printf("p%i: %lE ",i,getParam(i));
	}
	printf("\n");

}
//ostream& operator<<(MClink& link) {
//	for (long i=0;i<link.dims();i++) {
//		*sb << link.getParam(i) << " ";
//	}
//	
//}
