/*!
 * \file MClink.h - a Markov chain link container
 *
 *  Created on: Nov 19, 2010
 *      Author: blew
 */

#ifndef MCLINK_H_
#define MCLINK_H_

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */

/* STANDALONE HEADERS */
#include <stdio.h>
#include <stdarg.h>
#include "cpeds-math.h"
#include "cpeds-list.h"
/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */

/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
 \class 
 \brief Encapsulates a link in a Monte-Carlo Markov Chain 
 \details 
 
This class can store the information about a likelihood value at a given point in parameter space.
It also stores the values of parameters for that likelihood value.
 
 \date Nov 19, 2010, 7:29:03 PM 
 \author Bartosz Lew
 */

class MClink {
	public:
		MClink();
		MClink(long Nparams);
		MClink(const MClink& parent);
		virtual ~MClink();
		void setNpar(long n);
		void clear() { 	if (_theta!=NULL) { delete [] _theta; _theta=NULL; } }
		//! set the parameter values from provided values; 
		void set(int n, ...);
//		void set(long p, double v) { _theta[p]=v; }
		//! set the parameter values from array; The t array will be deleted here.
		void set(int n, double *t, bool deleteInside=true);
		//! set the likelihood value
		void setL(double v) {	_likelihood=v; }
		void setChisq(double v) {	_chisq=v; _likelihood=(exp(-_chisq/2)); }
		
		//! set both the parameter values and the likelihood value for the MC link
		void setLink(int n, double *t,double L) { setL(L); set(n,t); }

		void setAccepted(bool tf) { _accepted=tf; }
		bool isAccepted() const { return _accepted; }
		
		void setIdx(long i) { _idx=i; }
		long getIdx() const { return _idx; }
		
		//! get the likelihood values
		double L() const { return _likelihood; }
		double chisq() const { return _chisq; }
		
		//! get the parameters array pointer. This is the internal structure of this object.
		double* params() { return _theta; }
		//! return a copy of the parameter array
		double* getParams() const { return cpeds_copy_array(_theta,dims()); }
		cpedsList<double> getParameters() const;
		double getParam(long i) const { return _theta[i]; }
		long dims() const { return _Nparam; }
		
		void printParams() const;
		void save(string fname);
		MClink& load(string fname);
		MClink& load(ifstream& ifs, long Nparam);
		
		const MClink& operator=(const MClink& rhs);
		const double& operator[](const long i) const;
		
		
	private:
		double* _theta; //!< point in parameter space
		int _Nparam;
		double _likelihood; //!< likelihood value at that point
		double _chisq;
		bool _accepted;
		long _idx;
		
};
#endif /* MCLINK_H_ */
