/*!
 * \file mscsVector.h
 *
 *  Created on: Jan 27, 2012
 *      Author: blew
 */

#ifndef MSCSVECTOR_H_
#define MSCSVECTOR_H_

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */

/* STANDALONE HEADERS */
#include <vector>
#include <numeric>
#include "Mscs-object.h"
/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */

/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
 \class 
 \brief Encapsulates 
 \details 
 
 \date Jan 27, 2012, 11:57:18 AM 
 \author Bartosz Lew
 */

template <typename T>
class mscsVector : public vector<T>, public mscsObject {
	public:
		mscsVector() : vector<T>(), mscsObject("mscsVector",Zero) {};
		mscsVector(const vector<T>& v) : vector<T>(v), mscsObject("mscsVector",Zero) {};
		virtual ~mscsVector() {};
		
		using vector<T>::reserve;
		using vector<T>::size;
		using vector<T>::resize;
		using vector<T>::at;
		
		
		mscsVector<T>& pushBackVec(mscsVector<T>& rhs);
		virtual T* toArray();
		virtual void setSize(long N);
		void printVector();
		
		mscsVector<T>& operator=(const T& rhs);
		mscsVector<T>& operator=(const mscsVector<T>& rhs);
		mscsVector<T>& operator=(const vector<T>& rhs);
		mscsVector<T>& operator+=(const T& rhs);
		mscsVector<T>& operator-=(const T& rhs);
		mscsVector<T>& operator*=(const T& rhs);
		mscsVector<T>& operator/=(const T& rhs);
//		vector<double> tmp;
};


template <typename T> mscsVector<T>& mscsVector<T>::pushBackVec(mscsVector<T>& rhs) { 
	long st=size();
	setSize(size()+rhs.size());
	for (unsigned long i = 0; i < rhs.size(); i++) {
		at(st+i)=rhs[i];
	}
	return *this;
}

template <typename T> void mscsVector<T>::setSize(long N) { resize(N,T(0.0)); }
template <typename T> mscsVector<T>& mscsVector<T>::operator=(const T& rhs) {
	long N=size();
//	this->resize(N,rhs);
	for (long i = 0; i < N; i++) {		at(i)=rhs;	}
	return *this; 
}
template <typename T> mscsVector<T>& mscsVector<T>::operator=(const mscsVector<T>& rhs) {
	vector<T>::operator=(rhs);
	return *this; 
}
template <typename T> mscsVector<T>& mscsVector<T>::operator=(const vector<T>& rhs) {
	vector<T>::operator=(rhs);
	return *this; 
}
template <typename T> mscsVector<T>& mscsVector<T>::operator+=(const T& rhs) {
	long N=size();
	for (long i = 0; i < N; i++) {		at(i)+=rhs;	}
	return *this; 
}
template <typename T> mscsVector<T>& mscsVector<T>::operator-=(const T& rhs) {
	long N=size();
	for (long i = 0; i < N; i++) {		at(i)-=rhs;	}
	return *this; 
}

template <typename T> mscsVector<T>& mscsVector<T>::operator*=(const T& rhs) {
	long N=size();
	for (long i = 0; i < N; i++) {		at(i)*=rhs;	}
	return *this; 
}
template <typename T> mscsVector<T>& mscsVector<T>::operator/=(const T& rhs) {
	long N=size();
	for (long i = 0; i < N; i++) {		at(i)/=rhs;	}
	return *this; 
}

template <typename T> T* mscsVector<T>::toArray() {
	T* a=new T[size()];
	for (unsigned long i = 0; i < size(); i++) {
		a[i]=at(i);
	}
	return a;
}
template <typename T> void mscsVector<T>::printVector() { 
	long N=size();
	msgs->say("printing vector",High);
	for (long i = 0; i < N; i++) {		cout << at(i) << "\n"; 	}	
}

#endif /* MSCSVECTOR_H_ */
