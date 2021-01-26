/*!
 * \file cpedsMC.cpp
 *
 *  Created on: Jan 23, 2018
 *      Author: blew
 */
#include <cpedsMC.h>
#include "cpeds-templates.h"
#include <iostream>
#include <fstream>

using namespace cpems;

cpedsMC::cpedsMC() {
	// TODO Auto-generated constructor stub
	
}

cpedsMC::~cpedsMC() {
	// TODO Auto-generated destructor stub
}

void cpedsMC::load(string fname) {
	cpeds_queue<long>* q;
	q=cpeds_get_txt_file_cols_rows(fname.c_str());
	if (q==0) { return; }	

	long Nrows=q->get_size();
	long Ncols=q->getq(0);
	long Nparam=Ncols-4;
	
//	printf("Ncols: %li, Nrows: %li\n",Ncols,Nrows);
	
	ifstream ifs;
	ifs.open(fname.c_str());
	
	MClink l(Nparam);
	double data[Nparam];
	double X2,L,idx,accepted;
	
	for (long j = 0; j < Nrows; j++) {
		l.load(ifs,Nparam);
/*
		for (long i = 0; i < Nparam; i++) {
			ifs >> data[i];
			cout << data[i] << " ";
		}	

		ifs >> X2;
		ifs >> L;
		ifs >> idx;
		ifs >> accepted;

		cout << X2 << " " << L << " " << idx << " " << accepted << "\n";
		
		l.set(Nparam,data,false);
		l.setChisq(X2);
		l.setL(L);
		l.setAccepted(accepted);
*/
		append(l);
	}
	ifs.close();

}
/***************************************************************************************/
const cpedsMC& cpedsMC::operator=(const cpedsMC& rhs) {
	if (this!=&rhs) {
		QList<MClink>::operator=(rhs);
	}
	return *this;
}
/* ******************************************************************************************** */
void cpedsMC::save(string fname, long precision) {
	ofstream f;
	f.open(fname);
	f << std::setprecision(precision);
	for (auto it=this->begin();it!=this->end();it++) {
		f << it->getIdx() << " " << *it;
	}
	f.close();
}
/* ******************************************************************************************** */
void cpedsMC::append(MClink l, long idx) {
	l.setIdx(idx);
	QList::append(l);
}
/* ******************************************************************************************** */
void cpedsMC::append(MClink l) {
	QList::append(l);
}
/* ******************************************************************************************** */
void cpedsMC::append(cpedsMC &mc) {
	QList::append(mc);
}
