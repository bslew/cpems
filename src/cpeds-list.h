/*!
  \file extends the QList for operations and IO methods
 */

#ifndef CPEDSLIST
#define CPEDSLIST

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */


/* STANDALONE HEADERS */
#include <stdio.h>
#include<stdlib.h>
#include <iostream>
#include <fstream>
#include<string.h>
#include <iomanip>
#include <complex>
//#include <QtCore/QList>
#include "qtbase.h"
#include "cpeds-math.h"
#include "Mscs-alm.h"
#include "mscsVector.h"
#include <sys/stat.h>
/* #include <sys/types.h> */
/* #include <unistd.h> */

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */
using namespace std;

/* **************************************************************************************************** */
/* CLASS DEFINITION */
/* **************************************************************************************************** */
/*!
  \class
  \brief Encapsulates a list
  \details


  \date 2010/02/16 11:24:31
  \author Bartosz Lew
 */
template <typename T>
class cpedsList : public QList<T> {
		
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PUBLIC MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	public:
		
		using  QList<T>::append;
		using  QList<T>::size;
		using  QList<T>::count;
		using  QList<T>::clear;
		using  QList<T>::value;
		using  QList<T>::at;
		using  QList<T>::mid;
		using  QList<T>::first;
		using  QList<T>::last;
		using  QList<T>::operator[];
		using  QList<T>::operator<<;
		using  QList<T>::takeFirst;
		using  QList<T>::takeLast;
		
		/* ------------- */
		/* CLASS FRIENDS */
		/* ------------- */
		
		
		/* ---------------------------- */
		/* CONSTRUCTORS AND DESTRUCTORS */
		/* ---------------------------- */
		cpedsList() : QList<T>() {}
		cpedsList(T * tab, long N) : QList<T>() { fromCarray(tab,N); }
		cpedsList(const QList<T>& q) : QList<T>(q) {}
		cpedsList(const mscsVector<T>& v);
		~cpedsList() {}
		
		
		
		/* ---------------------------- */
		/* PUBLIC METHODS */
		/* ---------------------------- */
		
		/*!
   \brief makes the list to have N cells.
   \details
   @param N - requested number of cells in the list

   If the current number of cells is larger than N then elements from the beginning of the list
   will be removed until the requested number of elements is reached.

   If the current number of cells is smaller than N then default value (0) elements will be
   added to the end of the list.

   If N is negaitve then nothing will happen
		 */
		void makeLength(long N) {
			if (N<0) return;
			if (N==0) { QList<T>::clear(); return; }
			if (count() < N) { while (count() < N) { QList<T>::append(0); }   }
			if (count() > N) { while (count() > N) { QList<T>::pop_front();  }   }
		}
		
		/*!
   \brief derives the extreme values in the list and their positions in the list
   \details
   @param min - pointer pointing where to store the minimal value
   @param max - pointer pointing where to store the maximal value
   @param mini - pointer pointing where to store the minimal value index
   @param maxi - pointer pointing where to store the maximal value index

   For complex type the extreme values are selected by their abs() value
   and the returned complex numbers have that value stored in the real parts
   
   If the list is of size 0 then mini and maxi are set to -1 and the function
   returns.
   \date 2010/02/18 12:52:51
   \author Bartosz Lew
		 */
		void getMinMaxValues(T* min, T* max, long* mini, long* maxi) const {
			if (size()==0) {
				*mini=-1;
				*maxi=-1;
				return;
			}
			T *t=toCarray();
			cpeds_find_minmax_value(t,QList<T>::size(),min,max,mini,maxi);
			delete t;
		}
		
		//! calculates the mininal value of the list
		T max() const { T* t=toCarray(); T m=(T)cpeds_find_max_value(t,(long)size(),0); delete t; return m; }
		//! calculates the mininal value of the list
		T min() const { T* t=toCarray(); T m=(T)cpeds_find_min_value(t,(long)size(),0); delete t; return m; }
		//! calculates the mean value of the list
		T mean() const { T* t=toCarray(); T m=(T)cpeds_mean_value(t,(long)size()); delete t; return m; }
		//! calculates the unbiased estimator of the variance of the list
		T variance() const { T* t=toCarray(); T v=(T)cpeds_variance(t,size()); delete t; return v; }
		//! calculates the standard deviation of the list
		T std() const { T* t=toCarray(); T v=(T)sqrt(cpeds_variance(t,size())); delete t; return v; }
		//! calculates the rms of the list
		T rms() const { T* t=toCarray(); T v=(T)cpeds_rms(t,size()); delete t; return v; }
		//! calculates the sum of the list
		T sum() const { T* t=toCarray(); T m=(T)cpeds_sum(t,size(),false); delete t; return m; }
		//! calculates the skewness of the list
		T skewness() const { T* t=toCarray(); T v=(T)cpeds_skewness(t,size()); delete t; return v; }
		//! calculates the kurtosis of the list
		T kurtosis() const { T* t=toCarray(); T v=(T)cpeds_kurtosis(t,size()); delete t; return v; }
		
		//! calculates the mean value of the list
		T median() const { 
			T* t=toCarray(); 
			T m=(T)cpeds_median_value(t,(long)size()); 
			delete t; 
			return m; 
		}
		
		
		cpedsList<T>& sort(int dir=12) { 
			T* a= toCarray();
			long n=size();
			cpeds_sort_data(n,a,dir); 
			clear();
			fromCarray(a,n);
			delete [] a;
			return *this; 
		}
		
		
		cpedsList<T>& add(T v) { for (long i=0;i<count();i++) { (*this)[i]+=v; } return *this; }
		cpedsList<T>& subtract(T v) { for (long i=0;i<count();i++) { (*this)[i]-=v; } return *this; }
		cpedsList<T>& multiply(T v) { for (long i=0;i<count();i++) { (*this)[i]*=v; } return *this; }
		cpedsList<T>& divide(T v) { for (long i=0;i<count();i++) { (*this)[i]/=v; } return *this; }
		
		/* cpedsList<T>& add(long v) { for (long i=0;i<count();i++) { at(i)+=v; } return *this; } */
		/* cpedsList<T>& subtract(long v) { for (long i=0;i<count();i++) { at(i)-=v; } return *this; } */
		/* cpedsList<T>& multiply(long v) { for (long i=0;i<count();i++) { at(i)*=v; } return *this; } */
		/* cpedsList<T>& divide(long v) { for (long i=0;i<count();i++) { at(i)/=v; } return *this; } */
		
		cpedsList<T>& add(const cpedsList<T>& cl) { for (long i=0;i<count();i++) { (*this)[i]+=cl[i]; } return *this; }
		cpedsList<T>& subtract(const cpedsList<T>& cl) { for (long i=0;i<count();i++) { (*this)[i]-=cl[i]; } return *this; }
		cpedsList<T>& multiply(const cpedsList<T>& cl) { for (long i=0;i<count();i++) { (*this)[i]*=cl[i]; } return *this; }
		cpedsList<T>& divide(const cpedsList<T>& cl) { for (long i=0;i<count();i++) { (*this)[i]/=cl[i]; } return *this; }
		
		cpedsList<T>& invert() { for (long i=0;i<count();i++) { (*this)[i]=1.0/(*this)[i]; } return *this; }
		
		
		cpedsList<T>& operator+=(const cpedsList<T>& cl) { add(cl); return *this; }
		cpedsList<T>& operator-=(const cpedsList<T>& cl) { subtract(cl); return *this; }
		cpedsList<T>& operator*=(const cpedsList<T>& cl) { multiply(cl); return *this; }
		cpedsList<T>& operator/=(const cpedsList<T>& cl) { divide(cl); return *this; }
		
		cpedsList<T>& operator+=(T v) { add(v); return *this; }
		cpedsList<T>& operator-=(T v) { subtract(v); return *this; }
		cpedsList<T>& operator*=(T v) { multiply(v); return *this; }
		cpedsList<T>& operator/=(T v) { divide(v); return *this; }
		
		cpedsList<T>& operator=(const cpedsList<T>& cl) { if (this!=&cl) { QList<T>::operator=(cl); }  return *this; }
		cpedsList<T>& operator=(T v) { for (long i=0;i<count();i++) { (*this)[i]=v; } return *this; }
		/* cpedsList<T>& operator=(const cpedsList& cl) {  (QList<T>::this); return *this; } */
		cpedsList<T>& operator==(const cpedsList<T>& cl) { if (count()!=cl.count()) return false; for (long i=0;i<cl.count();i++) { if (at(i)!=cl[i]) return false; } return true; }
		
		cpedsList<T> operator+(const cpedsList<T>& cl) { cpedsList<T> tmp(*this); tmp.add(cl); return tmp; }
		cpedsList<T> operator-(const cpedsList<T>& cl) { cpedsList<T> tmp(*this); tmp.subtract(cl); return tmp; }
		cpedsList<T> operator*(const cpedsList<T>& cl) { cpedsList<T> tmp(*this); tmp.multiply(cl); return tmp; }
		cpedsList<T> operator/(const cpedsList<T>& cl) { cpedsList<T> tmp(*this); tmp.divide(cl); return tmp; }
		
		cpedsList<T> operator+(T v) const { cpedsList<T> tmp(*this); tmp.add(v); return tmp; }
		cpedsList<T> operator-(T v) const { cpedsList<T> tmp(*this); tmp.subtract(v); return tmp; }
		cpedsList<T> operator*(T v) const { cpedsList<T> tmp(*this); tmp.multiply(v); return tmp; }
		cpedsList<T> operator/(T v) const { cpedsList<T> tmp(*this); tmp.divide(v); return tmp; }
		
		bool operator>(T v) { for (long i=0;i<count();i++) { if (value(i)<=v) return false; } return true; }
		bool operator<(T v) { for (long i=0;i<count();i++) { if (value(i)>=v) return false; } return true; }
		bool operator>=(T v) { for (long i=0;i<count();i++) { if (value(i)<v) return false; } return true; }
		bool operator<=(T v) { for (long i=0;i<count();i++) { if (value(i)>v) return false; } return true; }
		
		
		
		long save(string file, bool binary=false, string type="double", long precision=10, bool append=false) const;
		long load(string file, bool binary=false, string type="double", bool Append=false);
		
		void print() const { long i,N=size();  for (i=0;i<N;i++) { cout <<value(i)<<endl; } }
		
		cpedsList<T>& fromMscsVector(const mscsVector<T>& t);
		mscsVector<T> toMscsVector();
		
		void fromCarray(T* t, long N) {
			for (long i=0;i<N;i++) { QList<T>::append(t[i]); }
		}
		
		T* toCarray() const {
			T * t = new T[size()];
			for (long i=0;i<size();i++) { t[i]=value(i); }
			return t;
		}
		
		T* toArray() const {
			return toCarray();
		}
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PROTECTED MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	protected:
		
		
		/* ---------------------------- */
		/* PROTECTED METHODS */
		/* ---------------------------- */
		
		
		/* ---------------------------- */
		/* PROTECTED STRUCTURES */
		/* ---------------------------- */
		
		
		/* ---------------------------------------------------------------------------------------------------- */
		/* CLASS PRIVATE MEMBERS */
		/* ---------------------------------------------------------------------------------------------------- */
	private:
		
		
		/* ---------------------------- */
		/* PRIVATE METHODS */
		/* ---------------------------- */
		
		
		/* ---------------------------- */
		/* PRIVATE STRUCTURES */
		/* ---------------------------- */
		
		
};



/* **************************************************************************************************** */
/* **************************************************************************************************** */
/* **************************************************************************************************** */
/* **************************************************************************************************** */
template <typename T> cpedsList<T>::cpedsList(const mscsVector<T>& v) {
	long i,N=v.size();
	//	printf("v size: %li\n",N);
	for (i=0;i<N;i++) { append(v[i]);  }
}

template <typename T> cpedsList<T>& cpedsList<T>::fromMscsVector(const mscsVector<T>& v) {
	long i,N=v.size();
	for (i=0;i<N;i++) { append(v[i]);  }
	return *this;
}
template <typename T> mscsVector<T> cpedsList<T>::toMscsVector() {
	mscsVector<T> v;
	v.resize(size());
	for (unsigned long i = 0; i < size(); i++) {
		v[i]=at(i);
	}
	return v;
}

template <typename T> long cpedsList<T>::save(string file, bool binary, string type, long precision, bool append) const {
	ofstream F;
	long i;
	long N=size();
	F.precision(precision);
	
	if (binary) {
		if (append) F.open(file.c_str(),ios_base::out | ios_base::binary | ios_base::app);	
		else F.open(file.c_str(),ios_base::out | ios_base::binary);	
		
		if (F.is_open()==false) return -1;
		
		T v;
		//		for (i=0;i<N;i++) { v=value(i); F << v;  }
		for (i=0;i<N;i++) { v=value(i); F.write((char*)(&v),sizeof(v));  }
		
	}
	else {
		if (append) F.open(file.c_str(),ios::out  | ios::app);
		else F.open(file.c_str(),ios::out);	
		
		if (F.is_open()==false) return -1;
		
		T v;
		for (i=0;i<N;i++) { v=value(i); F << v << "\n";  }
		
	}
	
	//	print();
	
	F.close();
	return 0;
}

//template <typename T> long cpedsList<T>::save(string file, bool binary, string type, long precision, bool append) const {
//  FILE* F;
//  long i;
//  long N=size();
//  if (binary) {
//	  if (append) F=fopen(file.c_str(),"ab");
//	  else F=fopen(file.c_str(),"wb");
//	  
//	  if (F==NULL) { return -1; }
//	  
//	  if (type == "complex") {
//		  mscsAlm v;
//		  struct {double  RE, IM; } c;
//		  size_t s=sizeof(c);
//		  for (i=0;i<N;i++) { v=value(i); c.RE=v.real(); c.IM=v.imag();	fwrite(&c,s,1,F);  }
//	  }
//	  else {
//		  /* if (type == "double" || type == "long") { */
//		  T *v = toCarray();
//		  fwrite(v,sizeof(T),N,F);
//		  delete v;
//	  }
//  }
//  else {
//	  string s;
//	  stringstream ss;
//	  if (append) F=fopen(file.c_str(),"ab");
//	  else F=fopen(file.c_str(),"wb");
//	  //    F=fopen(file.c_str(),"w");
//	  if (F==NULL) { return -1; }
//	  
//	  if (type == "double") { ss<<"%."<<precision<<"lE\n"; s=ss.str(); }
//	  if (type == "long") { s="%li\n"; }
//	  if (type == "complex") { ss<<"%."<<precision<<"lE %."<<precision<<"lE\n"; s=ss.str(); }
//	  
//	  if (type == "complex" || type == "mscsAlm" ) {
//		  mscsAlm v;
//		  for (i=0;i<N;i++) { v=value(i); fprintf(F,s.c_str(),v.real(),v.imag()); }
//	  }
//	  else {
//		  T v;
//		  for (i=0;i<N;i++) { 
//			  v=at(i); 
////			  fprintf(F,s.c_str(),v.real(),v.imag());
//		  	  fprintf(F,s.c_str(),v); 
//		  }
//	  }
//	  
//  }
//  fclose(F);
//  return 0;
//}

/* **************************************************************************************************** */
template <typename T> long cpedsList<T>::load(string file, bool binary, string type, bool Append) {
	FILE* F;
	long i;
	long N=0;
	
	if (!Append) clear();
	
	if (binary) {
		struct stat info;
		stat(file.c_str(),&info);
		if (type=="double")  N = (long)(info.st_size)/(long)sizeof(double);
		if (type=="long")  N = (long)(info.st_size)/(long)sizeof(long);
		if (type=="complex")  N = (long)(info.st_size)/(long)sizeof(complex<double>);
		if (type=="mscsAlm")  N = (long)(info.st_size)/(2*(long)sizeof(double));
#ifdef DEBUG_CPEDS_LIST
		printf("**** DEBUG **** reading list from binary file: number of items to read: %li\n",N);
#endif
		
		if (N==0) return -1;
		//makeLength(N);
		F=fopen(file.c_str(),"rb");    if (F==NULL) { return -1; }
		
		if (type == "mscsAlm") {
			struct {double  RE, IM; } c;
			size_t s=sizeof(c);
			/* for (i=0;i<N;i++) { fread(&c,s,1,F); QList<T>::append(mscsAlm(c.RE,c.IM));  } */
			long n;
			for (i=0;i<N;i++) { n=fread(&c,s,1,F); QList<T>::append(T( mscsAlm(c.RE,c.IM) ));  if (i<N and n==0) return -1; }
		} //else
		/* if (type == "complex") { */
		/* 	struct {double  RE, IM; } c; */
		/* 	size_t s=sizeof(c); */
		/* 	for (i=0;i<N;i++) { fread(&c,s,1,F); QList<T>::append(complex<double>(c.RE,c.IM));  } */
		/* } */
		else {
			
			T *v = new T[N];
			fread(v,sizeof(T),N,F);
			fromCarray(v,N);
			delete [] v;
			//			printf("size of T is : %i\n",int(sizeof(T)));
			//			T tmp;
			//			while (fread(&tmp,sizeof(T),1,F)>0) {				append(tmp);			printf("read: %lf\n",tmp); }
			
		}
		fclose(F);
	}
	else {
		string s;
		stringstream ss;
		int ret;
		F=fopen(file.c_str(),"r");
		if (F==NULL) { return -1; }
		
		if (type == "double") { s="%lE"; }
		if (type == "long") { s="%li"; }
		if (type == "complex") { s="%lE %lE"; }
		
		if (type == "mscsAlm") {
			double RE,IM;
			while (!feof(F)) {  ret=fscanf(F,s.data(),&RE,&IM);        if (ret!=EOF) append(mscsAlm(RE,IM)); }
		}
		else
			if (type == "complex") {
				double RE,IM;
				while (!feof(F)) {  ret=fscanf(F,s.data(),&RE,&IM);  	if (ret!=EOF) append(mscsAlm(RE,IM)); }
			}
			else {
				T v;
				while (!feof(F)) {  ret=fscanf(F,s.data(),&v);  	if (ret!=EOF) append(v); }
			}
		
		fclose(F);
		
	}
	return 0;
}
/* **************************************************************************************************** */
/* **************************************************************************************************** */
/* **************************************************************************************************** */
/* **************************************************************************************************** */
#endif /* CPEDSLIST */
