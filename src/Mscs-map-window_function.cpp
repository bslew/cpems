// in fact the structure with l and wl is highly redundant, should be enough to keep wl only
// in the right cell in memory; this can be changed at some point


/*   two different definitions of gaussian kernels used in this program: */
/*   G_l = exp[-0.5 * sigma^2 * l * (l+1) ]  as in J. Schmalzing and K.M. Gorski, Mon.Not. R. Astron. Soc. 297,355-365 (1988) */
/*   or */
/*   used at the moment for simulations */
/*   G_l = exp[ - 0.5 * sigma^2 * (l+1/2)^2 ] as in Changbom Park, et al. astro-ph/9711057 */
/* THESE TWO DEFINITIONS ARE CONSISTENT TO WITHIN A FACTOR 1.6*10^-6 !!! */




#include <math.h>
#include "Mscs-map-window_function.h"


/* **************************************************************************************************** */
// CONSTRUCTORS
/* **************************************************************************************************** */

mscsWindowFunction::mscsWindowFunction() : mscsFunction() {
	setVerbosityLevel(Zero);
}

/* **************************************************************************************************** */

mscsWindowFunction::mscsWindowFunction(string _wf_name, long _lmax,cpeds_VerbosityLevel verbosity) : mscsFunction(_wf_name,verbosity) {
	
	for (long i=0;i<=_lmax;i++) {
		newPoint(double(i),double(0.0));
	}
	
	/*   // initiate the functional space */
	/*   size = _size+1; */
	/*   wf = new w_f[size]; */
	
	/*   // initiate functional parameters */
	/*   min_i = 0; */
	/*   max_i = size; */
	/*   wf_name = _wf_name; */
	
	/* /\*   if (wf_name == "unit_kernel") make_unit_kernel(); *\/ */
	/* /\*   if (wf_name == "unit_kernel") make_unit_kernel(); *\/ */
	
	/*   // zero other info */
	/*   min_l = max_l = min_wl = max_wl = 0; */
}


/* **************************************************************************************************** */
mscsWindowFunction::mscsWindowFunction(string _wf_name, double _kmax) : mscsFunction(_wf_name) {
	long kmax=long(_kmax);
	for (long i=0;i<=kmax;i++) {
		newPoint(double(i),double(0.0));
	}
}
mscsWindowFunction::mscsWindowFunction(string _wf_name, cpeds_VerbosityLevel verbosity) : mscsFunction(_wf_name,verbosity) {
	
}
/* **************************************************************************************************** */
// DESTRUCTOR
mscsWindowFunction::~mscsWindowFunction() {
	/*   delete [] wf; */
}

/* **************************************************************************************************** */

// HANDLERS
/* double mscsWindowFunction::get_l(long i) { */
/*   double ret; */
/*   if (check_irange(i))  ret = wf[i].l; else ret = 0; */
/*   return ret; */
/* } */

/* double mscsWindowFunction::get_wl(long i) { */
/*   double ret; */
/*   if (check_irange(i)) ret = wf[i].wl; else ret = 0; */
/*   return ret; */
/* } */

/* void mscsWindowFunction::set_l(long i, double l) { */
/*   wf[i].l = l; */
/* } */

/* void mscsWindowFunction::set_wl(long i, double wl) { */
/*   wf[i].wl = wl; */
/* } */

/* string mscsWindowFunction::get_name() { */
/*   return wf_name; */
/* } */

/* void mscsWindowFunction::set_name(string _wf_name) { */
/*   wf_name = _wf_name; */
/* } */

/* // this returns the range of table indexes that are allowed */
/* void mscsWindowFunction::get_irange(long *_min_i, long *_max_i) { */
/*   *_min_i = min_i; *_max_i = max_i; */
/* } */

/* long mscsWindowFunction::get_size() { */
/*   return size; */
/* } */

/* void mscsWindowFunction::clear() { */
/* /\*  to be implemented *\/ */
/* } */

/* void mscsWindowFunction::get_frange(double * _min_l, double * _max_l, double * _min_wl, double * _max_wl) { */
/*   calculate_functional_minmax_values(); */
/*   *_min_l = min_l;   *_max_l = max_l; */
/*   *_min_wl = min_wl;   *_max_wl = max_wl; */
/* } */


/* //-------------------------------------------------------------------------------- */
/* // prrivate methods */
/* //-------------------------------------------------------------------------------- */

/* // checking index range passed to handlers */
/* bool mscsWindowFunction::check_irange(long i) { */
/*   bool ret; */
/*   if ((i <= max_i) && (i >= min_i)) ret = true;   else ret = false; */
/*   return ret; */
/* } */

/* void mscsWindowFunction::calculate_functional_minmax_values() { */
/*   long i; */

/*   min_l = wf[0].l;  max_l = wf[0].l; */
/*   min_wl = wf[0].wl;  max_wl = wf[0].wl; */

/*   for (i=min_i;i<max_i;i++) {  */
/*     if (wf[i].l < min_l) { min_l = wf[i].l; min_li = i; } */
/*     if (wf[i].l > max_l) { max_l = wf[i].l; max_li = i; } */
/*     if (wf[i].wl < min_wl) { min_wl = wf[i].wl; min_wli = i; } */
/*     if (wf[i].wl > max_wl) { max_wl = wf[i].wl; max_wli = i; } */
/*   } */

/* } */


/* **************************************************************************************************** */

const mscsWindowFunction&  mscsWindowFunction::make_gaussian_kernel(double FWHM) {
	long i;
	double s = FWHM/(2*sqrt(2*log(2.0)));
	double l;
	
	msgs->say("making gaussian window function kernel. Points:"+msgs->toStr(pointsCount()),Medium);
	long num=pointsCount();
	for (i=0;i<num;i++) {
		setf(i,double(i),exp(-pow((double)i+0.5,2)*s*s/2));
	}
	
	checkRanges();
	printRanges();
	
	/*   for (i=min_i;i<max_i;i++) {  */
	/*     wf[i].l = (double)i; */
	/*     wf[i].wl = exp(-pow((double)i+0.5,2)*s*s/2); */
	/*   } */
	/*   calculate_functional_minmax_values(); */
	return *this;
}

/* **************************************************************************************************** */
const mscsWindowFunction&  mscsWindowFunction::make_gaussian_kernel_kspace(double FWHM) {
	long i;
	double s = FWHM/(2*sqrt(2*log(2.0)));
	double s2o2=s*s/2.0;
	double l;
	
	msgs->say("making gaussian window function kernel in k-space\n",Medium);
	long num=pointsCount();
	for (i=0;i<num;i++) {
		setf(i,i,exp(-double(i)*double(i)*s2o2));
	}
	checkRanges();
	
	
	/*   for (i=min_i;i<max_i;i++) {  */
	/*     wf[i].l = (double)i; */
	/*     wf[i].wl = exp(-pow((double)i+0.5,2)*s*s/2); */
	/*   } */
	/*   calculate_functional_minmax_values(); */
	return *this;
}
/* **************************************************************************************************** */
mscsFunction& mscsWindowFunction::mkTopHatKernel(double R, int fftConvention) {
	msgs->say("making top-hat kernel window function in k-space for scale R\n",Medium);
	long num=pointsCount();
	double k;
	double kR;
	if (fftConvention==0) {
		for (long i=0;i<num;i++) {
			k=getx(i);
			kR=k*R;
			if (kR==0) 
				setf(i,1.0 );
			else
				setf(i, sin(twoPI*kR)/(twoPI*kR) );
		}		  
	}
	else {
		if (fftConvention==1) {
			for (long i=0;i<num;i++) {
				k=getx(i);
				kR=k*R;
				if (kR==0) 
					setf(i,1.0 );
				else
					setf(i,3.0* ( sin(kR)/(kR*kR*kR)-cos(kR)/(kR*kR) ) );
			}		  
		}
	}
	checkRanges();	
	return *this;
}
/* **************************************************************************************************** */
mscsFunction& mscsWindowFunction::mkGaussianKernel(double R, double k0) {
	msgs->say("making gaussian kernel window function in k-space for scale R\n",Medium);
	long num=pointsCount();
	double k;
	double kR;
	for (long i=0;i<num;i++) {
		k=getx(i);
		if (k<k0)	setf(i,1.0);
		else {
			kR=(k-k0)*R;
			setf(i,exp(- kR*kR /2.0 * 4.0*PI*PI));			  
		}
	}
	checkRanges();	
	return *this;
}

/* **************************************************************************************************** */
const mscsWindowFunction& mscsWindowFunction::make_unit_kernel(long lmax) {
	long i;
	double l;
	
	msgs->say("making unit window function kernel with pointsCount: "+msgs->toStr(pointsCount()),Medium);
	
	long num;
	if (lmax!=-1) setPointsNum(lmax+1);
	num=pointsCount();
	for (i=0;i<num;i++) {
		setf(i,i,double(1));
	}
	checkRanges();
	
	/*   for (i=min_i;i<max_i;i++) {  */
	/*     wf[i].l = (double)i; */
	/*     wf[i].wl = (double)1.0; */
	/*   } */
	/*   calculate_functional_minmax_values(); */
	return *this;
}
/* **************************************************************************************************** */
mscsFunction& mscsWindowFunction::mkHannWindow(long N) {
	long Nleo=N-1;
	clearFunction();
	for (long i = 0; i < N; i++) {
		newPoint(i,0.5*(1.0-cos(twoPI*i/Nleo)));
	}
	return *this;
}
/***************************************************************************************/
mscsFunction& mscsWindowFunction::mkExpWindow(double kmin, double kmax, long Nk, double k0, double tC) {
	double k=kmin;
	double dk=(kmax-kmin)/double(Nk);
	while (k<kmax) {
		if (k<k0)
			newPoint(k, 1.0);
		else
			newPoint(k, exp(-(k-k0)/tC) );
		k+=dk;
	}
	return *this;		
}
/***************************************************************************************/
mscsFunction& mscsWindowFunction::mkHighPassExpWindow(double kmin, double kmax, long Nk, double k0, double tC) {
	double k=kmin;
	double dk=(kmax-kmin)/double(Nk);
	while (k<kmax) {
		if (k>k0)
			newPoint(k, 1.0);
		else
			newPoint(k, exp(-(k0-k)/tC) );
		k+=dk;
	}
	return *this;			
}
/***************************************************************************************/
mscsFunction& mscsWindowFunction::mkLowPassWindow(double fmin, double fmax, long Nf, double RC) {
	double nu=fmin;
	double dnu=(fmax-fmin)/double(Nf-1);
	double tmp;
	while (nu<fmax) {
		tmp=twoPI*nu*RC;
		//		newPoint(nu, 1.0/sqrt(1.0+tmp*tmp) );
		newPoint(nu, 1.0/(1.0+tmp*tmp) );
		nu+=dnu;
	}
	return *this;		
}
/***************************************************************************************/
mscsFunction& mscsWindowFunction::mkLowPassStepWindow(double fmin, double fmax, long Nf, double RC) {
	double nu=fmin;
	double dnu=(fmax-fmin)/double(Nf-1);
	double fc=1.0/(twoPI*RC);
	
	while (nu<fmax) {
		if (nu<=fc)	newPoint(nu, 1.0 );
		else newPoint(nu, 0.0 );
		nu+=dnu;
	}
	return *this;			
}
/***************************************************************************************/
mscsFunction & mscsWindowFunction::mkSPHkernelGadget2(long  N, double Rmin, double Rmax, double hsml) {
	if (Rmin <0) {
		msgs->criticalError("mscsWindowFunction::mkSPHkernelGadget2>>> Rmin cannot be smaller than 0,", Top);
	}
	double R=Rmin;
	double dR=(Rmax-Rmin)/double(N-1);
	double tmp;
	double norm=1.0;
	long i=0;
	if (hsml!=0) {		norm=8.0/(PI*hsml*hsml*hsml);	}
	
	while (R<Rmax) {
		if (R<=0.5) { tmp=1.0-6*(R*R)+6*(R*R*R); }
		else {
			if (R<=1.0) { tmp=1.0-R; tmp=tmp*tmp*tmp; tmp=2.0*tmp; }
			else tmp=0;
		}
		
		newPoint(R, norm*tmp);
		i++;
		R=dR*i;
	}
	return *this;			
}


/***************************************************************************************/
double mscsWindowFunction::kernelGadget(double R) {
	double tmp;
//#pragma omp threadprivate(tmp)
	if (R<=0.5) { tmp=1.0-6*(R*R)+6*(R*R*R); }
	else {
		if (R<=1.0) { tmp=1.0-R; tmp=tmp*tmp*tmp; tmp=2.0*tmp; }
		else tmp=0;
	}
	return tmp;
}
/***************************************************************************************/
double mscsWindowFunction::kernelGadget2b(double R) {
	return kernelGadget(R/2);
}
/***************************************************************************************/
mscsFunction& mscsWindowFunction::mkAirMassElevationRelation(double from, double to, double dh) {
	double h=from;
	double fromLoc=from;
	double minimalElevation=4.0;
	if (from<minimalElevation) fromLoc=minimalElevation; else fromLoc=from;
	h=fromLoc;
	while (h<=to) {		
		newPoint(h,cpeds_air_mass(h));
		h+=dh;
	}
	if (from<minimalElevation) {
		extrapolate(from,minimalElevation,dh,true,"linear");
	}
	return *this;				
}
