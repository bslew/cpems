/*!
 * \file MClink.cpp
 *
 *  Created on: Nov 19, 2010
 *      Author: blew
 */

#include "MClink.h"

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
void MClink::set(int n, double *t) {
	setNpar(n);
	for (long i = 0; i < n; i++) {		_theta[i]=t[i];	}
	delete [] t;
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
	if (this!=&rhs) {
		if (_Nparam!=rhs.dims()) {
			clear();
			_theta=rhs.getParams();
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
	return *this;
}
/***************************************************************************************/
const double& MClink::operator[](const long i) const {
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
	fprintf(f,"%li ", _idx);
	long acc;
	if (isAccepted()) acc=1; else acc=0;
	fprintf(f,"%li ", acc);
	fprintf(f,"\n");
	fclose(f);
}
/***************************************************************************************/
cpedsList<double> MClink::getParameters() const {
	cpedsList<double> p;
	p.fromCarray(_theta,dims());
	return p;
}
