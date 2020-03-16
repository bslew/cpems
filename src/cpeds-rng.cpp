#include "cpeds-rng.h"
#include <fftw3.h>

// ****************************************************************************************************
void cpedsRNG::initiateRNG(string distr, string rn, const gsl_rng_type* generator) {
	_state=cpeds_random_uniform_numbersGSL_init(&_seed, seedOffset(),generator); // this is the default one in gsl: gsl_rng_mt19937
	_generator_type=generator;
	setMinMax(0,1);
	setMeanStd(0,1);
	setCentralLimitNumbers(100);
	gCDF=NULL;
	gCDF2d=NULL;
	CDFargs=NULL;
	CDFargsX=NULL;
	CDFargsY=NULL;
	gCDFsize=0;
	gCDFsize2d=0;
	gCDFsizeX=0;
	gCDFsizeY=0;
	_drawsCount=0;
	if (rn=="double") { _rnDataType=Double; }

	setRNsType(distr);

	_KasdinPowerLawNoiseMinimalCoefficient=2e-7;
//	_KasdinPowerLawNoiseMinimalCoefficient=0.121;
	_KasdinCoefficientsCalculated=false;
	_oofTot=0;
	_spectralIndex=-1.0;
	_oofnst=NULL;
	_bt=NULL;

	_pivotPoint=1.0;
	_oofAmplitude=1.0;
	_oofKmin=-1; // this will be converted to 1/N when N will be defined
	_oofKmax=0.5; // this defines the larges frequency in the signal
	
	
}

// ****************************************************************************************************
void cpedsRNG::clone(const cpedsRNG& parent) {
	_generator_type=parent._generator_type;
	seed(parent.seed());
//	exit(0);
	seedOffset(parent.seedOffset());
//	printf("rng clone: rng: %li, state: %li, type: %li, seed:%li off: %li\n",_state,
//			_state->state,_state->type,_seed,_seed_offset);
	setState(parent.getState());
	_rnDataType=parent._rnDataType;
	_distr=parent._distr;
	setMinMax(parent.Min(),parent.Max());
	setMeanStd(parent.mean(),parent.gauss_std());
	setCentralLimitNumbers(parent.centralLimitNumbers());
	copyGCDF(parent);
	_drawsCount=parent._drawsCount;
	_KasdinCoefficientsCalculated=parent._KasdinCoefficientsCalculated;
	_KasdinPowerLawNoiseMinimalCoefficient=parent._KasdinPowerLawNoiseMinimalCoefficient;
	_oofTot=parent.getOOfTot();
	_spectralIndex=parent.getSpectralIndex();
	_bt=NULL;
	_bt_size=0;
	_oofnst=NULL;
	_oofnst_curidx=0;
	
	_pivotPoint=parent.getPivotPoint();
	_oofAmplitude=parent.getOofAmplitude();
	_oofKmin=parent._oofKmin;
	_oofKmax=parent._oofKmax;
}

/***************************************************************************************/
void cpedsRNG::setRNsType(string distr) {

	if (distr=="uniform") { _distr=uniform; } else
		if (distr=="gaussian") { _distr=gaussian; } else
			if (distr=="gaussian_invcdf") {		_distr=gaussian_invcdf;	} else 
				if (distr=="gaussian_circle") {		_distr=gaussian_circle;	} else
					if (distr=="gaussianPowerLaw") {		_distr=gaussian_power_law;	} else 
						if (distr=="gaussianPowerLaw_t") {		_distr=gaussian_power_law_t;	} else 
							if (distr=="gaussianPowerLaw_t2") {		_distr=gaussian_power_law_t2;	} else 
								if (distr=="gaussianPowerLaw_fft") {		_distr=gaussian_power_law_fft;	} else 
									if (distr=="invCDF") {		_distr=from_array_invCDF;	} else 
										if (distr=="invCDF2d") {		_distr=from_array_invCDF2d;	} else 
										{						
											printf("unknown distribution: %s\n. Please check the code. Will exit now - this is to force debugging.\n",distr.c_str());
											exit(0);
										}
}
/***************************************************************************************/
void cpedsRNG::setRNsType(distrType distr) {
	switch (distr) {
		case from_array_invCDF:
			_distr=distr;
			break;
		case from_array_invCDF2d:
			_distr=distr;
			break;
		case uniform:
			_distr=distr;
			break;
		case gaussian:
			_distr=distr;
			break;
		case gaussian_invcdf:
			_distr=distr;
			break;
		case gaussian_circle:
			_distr=distr;
			break;
		case gaussian_power_law:
			_distr=distr;
			break;
		case gaussian_power_law_t:
			_distr=distr;
			break;
		case gaussian_power_law_t2:
			_distr=distr;
			break;
		case gaussian_power_law_fft:
			_distr=distr;
			break;
		default:
			printf("setRNsType: unknown distribution: n. Please check the code. Will exit now - this is to force debugging.\n");			
			break;
	}
}

/***************************************************************************************/
void cpedsRNG::setPDF(long size, double* x, double *p) {
	long i;
	if (size!=0) { killGCDF(); }
	gCDFsize=size-1;
	gCDF=p;
	CDFargs=x;
	double dxo2=0.5*(x[1]-x[0]);
//	double* cdf=new double[gCDFsize];
//	for (i=0;i<gCDFsize;i++) {	    cdf[i]=p[i];	  }
//	for (i=1;i<gCDFsize;i++) {	    cdf[i]+=cdf[i-1];	  }
	
	// convert to CDF
	for (i=1;i<gCDFsize;i++) {	    p[i]+=p[i-1];	  }
//	for (i=1;i<gCDFsize;i++) {	    p[i]=p[i-1]+0.5*(cdf[i-1]+cdf[i]);	  }
	// shift xargx
	for (i=0;i<gCDFsize;i++) {	    x[i]+=dxo2;	  }
	
//	for (i=0;i<gCDFsize;i++) {	    p[i]+=cdf[i];	  }
	
	// normalization of the CDF
	long last=gCDFsize-1;
	for (i=0;i<gCDFsize;i++) {	    p[i]/= p[last];	  }
//	delete [] cdf;
}
/***************************************************************************************/
void cpedsRNG::setPDF2d(long sizeX, double* x, long sizeY, double *y, double *p) {
	long i,j;
	if (sizeX!=0 and sizeY!=0) { killGCDF2d(); }
	gCDFsize2d=sizeX*sizeY;
	gCDFsizeX=sizeX;
	gCDFsizeY=sizeY;
	gCDF2d=p;
	CDFargsX=x;
	CDFargsY=y;
	// convert to CDF
	for (i=0;i<gCDFsizeX;i++) {
		for (j=1;j<gCDFsizeY;j++) {	    
			p[i*sizeY+j]+=p[i*sizeY+j-1];
		}
	}
	for (j=0;j<gCDFsizeY;j++) {	    
		for (i=1;i<gCDFsizeX;i++) {
			p[i*sizeY+j]+=p[(i-1)*sizeY+j];
		}
	}
	
	// normalization of the CDF
	long last=gCDFsize2d-1;
	for (i=0;i<gCDFsize2d;i++) {	    p[i]/= p[last];	  }
	
}

// ****************************************************************************************************
const cpedsList<double> cpedsRNG::getRNs(long n) {
	cpedsList<double> cl;
	if (n==0) return cl;
	double *t=getRN(n);
	for (long i=0;i<n;i++) { cl.append(t[i]); }
	delete [] t;
	return cl;
}
// ****************************************************************************************************
double* cpedsRNG::getRN(long n) {
	double *t=NULL;
	long i;
	/***************************************************************************************/
	if (_distr==uniform) {
		double delta=Max()-Min();
		double min=Min();
		t=new double[n];
		for (i=0;i<n;i++) { t[i]=  gsl_rng_uniform(_state) * delta + min; }

		_drawsCount+=n;
	}
	/***************************************************************************************/
	if (_distr==gaussian) {
		t=cpeds_random_gauss_numbers(mean(),gauss_std(),n, centralLimitNumbers(),getState());
		//for (i=0;i<n;i++) { t[i]= cpeds_random_gauss_numbers(mean(),variance(),n, centralLimitNumbers(),getState()); }

		_drawsCount+=n*centralLimitNumbers();
	}
	/***************************************************************************************/
	if (_distr==gaussian_invcdf) {
		t=new double[n];
		long int j,np=0; // the size of the tabulating array is to be twice as large as the number of points drown from this distribution
		double x,delta,xp,xf,yf,xi,yi;
		long istart,iend;
		if (gCDF==NULL || 2*n!=gCDFsize) {
			killGCDF();
			gCDFsize=2*n;
			gCDF = cpeds_generate_gaussian_distribuant_function(gCDFsize, gauss_std(), mean());
		}
		for (j=0;j<n;j++) {
			x = gsl_rng_uniform(getState());

			istart = 0; iend = gCDFsize;
			np = (long int)round(((double)iend - (double)istart)/2); xp=gCDF[cpeds_xy2num(np,1)];
			delta=xp-x;
//			printf("x: %lE ist: %li, ien: %li np: %li, delta: %lE, xp: %lE\n",x,istart,iend,np,delta,xp);
			while (np+1 != iend ) {
//				printf("   beg> x: %lE ist: %li, ien: %li np: %li, delta: %lE, xp: %lE\n",x,istart,iend,np,delta,xp);
				if (delta < 0) { istart = np; }
				if (delta > 0) { iend = np; }
				np = (int)round((double)istart+((double)iend - (double)istart)/2); xp=gCDF[cpeds_xy2num(np,1)];
				delta = xp - x;
//				printf("   end> x: %lE ist: %li, ien: %li np: %li, delta: %lE, xp: %lE\n",x,istart,iend,np,delta,xp);
			}

			// linear interpolation between the selected points
			if (x > xp) {xi = gCDF[cpeds_xy2num(np,0)]; xf = gCDF[cpeds_xy2num(np+1,0)]; yi = gCDF[cpeds_xy2num(np,1)]; yf = gCDF[cpeds_xy2num(np+1,1)];}
			if (x < xp) {xi = gCDF[cpeds_xy2num(np-1,0)]; xf = gCDF[cpeds_xy2num(np,0)]; yi = gCDF[cpeds_xy2num(np-1,1)]; yf = gCDF[cpeds_xy2num(np,1)];}
			t[j] = (x-yi)*(xf-xi)/(yf-yi) + xi;
//			printf("tj: %lE\n",t[j]);
//			t[j]=xp;
		}

		_drawsCount+=n;
	}
	/***************************************************************************************/
	if (_distr==gaussian_circle) {
		t=new double[n];
		double v=gauss_std();
		double m=mean();
		double v1,v2;
		double S;
		for (long j=0;j<n;j++) {
			do {
				v1=cpeds_random_uniform_numberGSL(-1,1,getState());
				v2=cpeds_random_uniform_numberGSL(-1,1,getState());
				S=v1*v1+v2*v2;
				_drawsCount+=2;
			} while (S>1 or S==0);
			t[j]=v * v1*sqrt(-2.0*log(S)/S) + m;
		}
	}

	/***************************************************************************************/
	if (_distr==gaussian_power_law) {
		_distr=gaussian_circle;
		setMeanStd(0,1);
		double* gauss=getRN(n); // _drawsCount will be increased here
		_distr=gaussian_power_law;
		
		//
		// calculate the coefficients
		//
		if (!_KasdinCoefficientsCalculated) {
			double b=1.0;
			double bp;
			double i=0;
			double ns=0.5*getSpectralIndex();
			while (fabs(b)>_KasdinPowerLawNoiseMinimalCoefficient) {
				_b.append(b);
				bp=b; 
				i++;
				b=(i-1.0+ns)*bp/i;
			}
			_KasdinCoefficientsCalculated=true;
		}
		
		//
		// calculate 1/f numbers
		//
		long bmax,N,Nb;
		Nb=_b.size();
		bmax=Nb-1;

		for (long i = 0; i < n; i++) {
			N=cpeds_get_min(bmax, _oofTot);
			if (_oofTot>bmax) { _oofns.takeFirst(); }

			_oofns.append(gauss[i]-sumOOfNumbers(N));
			gauss[i]=_oofns.last();
			_oofTot++;
		}

		//
		// finalize
		//
		t=gauss;
	}
	
	/***************************************************************************************/
	if (_distr==gaussian_power_law_t) {
		_distr=gaussian_circle;
		setMeanStd(0,1);
		double* gauss=getRN(n); // _drawsCount will be increased here
		_distr=gaussian_power_law_t;
		
		//
		// calculate the coefficients
		//
		if (!_KasdinCoefficientsCalculated) {
			double b=1.0;
			double bp;
			double i=0;
			double ns=0.5*getSpectralIndex();
			while (fabs(b)>_KasdinPowerLawNoiseMinimalCoefficient) {
				_b.append(b);
				bp=b; 
				i++;
				b=(i-1.0+ns)*bp/i;
			}
			_KasdinCoefficientsCalculated=true;
			
			// rewrite this onto a table
			_bt_size=_b.size();
			_bt=_b.toCarray();
			_b.clear();
			
			// setup an array for previos 1/f numbers
			_oofnst=new double[_bt_size];
			_oofnst_curidx=0;
		}
		
		//
		// calculate 1/f numbers
		//
		long bmax,N;
		bmax=_bt_size-1;
		for (long i = 0; i < n; i++) {
			N=cpeds_get_min(bmax, _oofTot);
			if (_oofTot>bmax) { 
				cpeds_shift_array(_oofnst,_bt_size,1,false); // commenting out this shift reduced the time generation time by factor of ~5 to 6
			}

			_oofnst[_oofnst_curidx]=gauss[i]-sumOOfNumbers_t(N);
			gauss[i]=_oofnst[_oofnst_curidx];
			if (_oofnst_curidx<bmax) _oofnst_curidx++;
			_oofTot++;
		}

		//
		// finalize
		//
		t=gauss;
	}
	/***************************************************************************************/
	if (_distr==gaussian_power_law_t2) {
		_distr=gaussian_circle;
		setMeanStd(0,1);
		double* gauss=getRN(n); // _drawsCount will be increased here
		_distr=gaussian_power_law_t2;
		
		//
		// calculate the coefficients
		//
		if (!_KasdinCoefficientsCalculated) {
			double b=1.0;
			double bp;
			double i=0;
			double ns=0.5*getSpectralIndex();
			while (fabs(b)>_KasdinPowerLawNoiseMinimalCoefficient) {
				_b.append(b);
				bp=b; 
				i++;
				b=(i-1.0+ns)*bp/i;
			}
			_KasdinCoefficientsCalculated=true;
			
			// rewrite this onto a table
			_bt_size=_b.size();
			_bt=_b.toCarray();
			_b.clear();
			
			// setup an array for previos 1/f numbers
			_oofnst=new double[_bt_size];
			_oofnst_curidx=0;
		}
		
		//
		// calculate 1/f numbers
		//
		long bmax,N;
		bmax=_bt_size-1;
		_oofnst_stfrom=0;
		long idx;
		for (long i = 0; i < n; i++) {
			N=cpeds_get_min(bmax, _oofTot);
			if (_oofTot>bmax) { 
//				_oofnst=cpeds_shift_array(_oofnst,_bt_size,1,false); // commenting out this shift reduced the time generation time by factor of ~5 to 6
																	/* but the solution with this _oofnst_stfrom business takes much longer time than the _t 
																	 * approach with the shifting array method - suprisingly even though the results returned
																	 * are identical*/
				_oofnst_stfrom++;
				_oofnst_stfrom=_oofnst_stfrom % _bt_size;
			}
			idx=(_oofnst_curidx + _oofnst_stfrom) % _bt_size;
			_oofnst[idx]=gauss[i]-sumOOfNumbers_t2(N);
			gauss[i]=_oofnst[idx];
			if (_oofnst_curidx<bmax) _oofnst_curidx++;
//			printf("bmax:%li curidx: %li\n",bmax,_oofnst_curidx);
			_oofTot++;
		}

		//
		// finalize
		//
		t=gauss;
	}
	/***************************************************************************************/	
	if (_distr==gaussian_power_law_fft) {
#ifdef DEBUG_CPEDS
		printf("WARNING !! double* cpedsRNG::getRN(long n); if (_distr==gaussian_power_law_fft): THIS ROUTUNE IS NOW CORRECTED BUT NOT TESTED. \n");	
#endif
		
//		mscsFunction& mscsFunction::mkPowerLawNoise(double A, double alpha, double k0, double kmin, double kmax, long N, long seed, mscsFunction* powerSpectrum, cpedsRNG* rng) {
//			cpedsRNG *rns;
		//
		// setup the random number generator
		//
		_distr=gaussian_circle;
		setMeanStd(0,1);
		//
		// generate gaussian noise
		//
		long M;		
		if (n%2==1) M=(n-1)/2+1; else M=n/2+1;
		double *a= getRN(M);
		double *b= getRN(M);
		_distr=gaussian_power_law_fft;
				
		//
		// form the right power spectrum in Fourier space (positive frequencies only)
		//
		double Pk,sqrPj;
		if (_oofKmin==-1) { _oofKmin=double(1)/double(n); } // define _oofKmin if it wans't defined yet
//		double Lmin=1/kmax;
//		double Lmax=1/kmin;
//		double L;
		double T; // this is the total sampling period T = N dt = N/(2 kmax)

//		L=fabs(Lmax-Lmin);
		long j; // iterates k_j; j=0..M-1
		double k;
		double dk;
		/*
		  //Jul 19, 2011 2:47:51 PM - DEBUG BEGIN - changed when fixing a bug
		  
		  double k=_oofKmin;
		  double dk=(_oofKmax-_oofKmin)/double(M);
		  
		  //DEBUG - END
		*/
		dk=double(2)*_oofKmax/double(n);
		T=double(n)/(2.0 * _oofKmax);
		_oofKmin=0;

		double alpha=getSpectralIndex();
		double k0=getPivotPoint();
		double A=getOofAmplitude();
		for (j=1;j<M;j++) {
			k=double(j)*dk+_oofKmin;
			Pk=A*powf( k / k0,alpha); // L=N dL; dL=L/N
			sqrPj=sqrt(Pk/2.0);
			a[j]*=sqrPj;
			b[j]*=sqrPj;
			k+=dk;
		}
		// set power for k=0;
		a[0]=0;
		b[0]=0;

		//
		// do fftw
		//

		fftw_complex* in=(fftw_complex*)fftw_malloc((M)*sizeof(fftw_complex)); 
		for (long j=0;j<M;j++) { in[j][0]=a[j]; in[j][1]=b[j]; }
		delete [] a; delete [] b; 
		
		double* out=new double[n];
		fftw_plan p;
		p=fftw_plan_dft_c2r_1d(n, in,out, FFTW_ESTIMATE); 
		fftw_execute(p);
		fftw_destroy_plan(p);
		fftw_free(in);
		fftw_cleanup();
		t=out;		
	}
	

	/***************************************************************************************/
	if (_distr==from_array_invCDF) {
	  double x,xp,xi,xf,yi,yf; 

	  // generate uniformly distributed numbers
	  _distr=uniform;
	  setMinMax(gCDF[0],gCDF[gCDFsize-1]);
	  t=getRN(n); // draws count is increased here
	  _distr=from_array_invCDF;
	  
	  // convert to requested PDF according to the CDF table
	  long k;
	  
//	  cpeds_save_matrix(gCDF,gCDFsize,1,"gCDF",false,false);
	  for (long j=0;j<n;j++) {
		  
		  // find the point in the cdf array
		  k=cpeds_find_value(t[j],gCDF,gCDFsize,0,gCDFsize);		  
		  x=t[j];
		  xp=gCDF[k];
//		  if (k==0 or k==1) printf("x: %lE k: %li, xp: %lE\n",x,k,xp);
//		  printf("x: %lE k: %li, xp: %lE\n",x,k,xp);
		  
//		  // linear interpolation between the selected points
		  if (x > xp) { xi = CDFargs[k]; xf = CDFargs[k+1]; yi = xp; yf = gCDF[k+1];}
		  if (x < xp) { 
			  xi = CDFargs[k-1]; xf = CDFargs[k]; yi = gCDF[k-1]; yf = gCDF[k];
//			  if (k==0) xi=0;
//		  	  if (k==0) printf("warning x: %lE k: %li, xp: %lE, xi:%lE, xf:%lE, yi:%lE yd:%lE\n\n",t[j],k,xp,  xi,xf,yi,yf);		  
		  }
		  
		  t[j] = (x-yi)*(xf-xi)/(yf-yi) + xi;
//		  if (k==0 or k==1) 		  printf("int x: %lE k: %li, xp: %lE, xi:%lE, xf:%lE, yi:%lE yd:%lE\n\n",t[j],k,xp,  xi,xf,yi,yf);
//		  t[j]=CDFargs[k]; // case for no interpolation
	  }
	}
	/***************************************************************************************/
	
	return t;
}
// ****************************************************************************************************
double cpedsRNG::getRN() {
	double *dp=getRN(1);
	double d=*dp;
	delete dp;
	return d;
}
/***************************************************************************************/
//cpedsPoint2D* cpedsRNG::getRN2d(long n) {
//	cpedsPoint2D *t=NULL;
//	long i;
//	msgs->criticalError("cpedsRNG::getRN2d> not implemented yet",Top);
//	/***************************************************************************************/
//	if (_distr==from_array_invCDF2d) {
//	  double x,xp,xi,xf,yi,yf; 
//
//	  // generate uniformly distributed numbers
//	  _distr=uniform;
//	  setMinMax(gCDF[0],gCDF[gCDFsize-1]);
//	  t=getRN(n); // draws count is increased here
//	  _distr=from_array_invCDF2d;
//	  
//	  // convert to requested PDF according to the CDF table
//	  long k;
//	  for (long j=0;j<n;j++) {
//		  
//		  // find the point in the cdf array
//		  k=cpeds_find_value(t[j],gCDF,gCDFsize,0,gCDFsize);		  
//		  x=t[j];
//		  xp=gCDF[k];
//		  
////		  // linear interpolation between the selected points
//		  if (x > xp) { xi = CDFargs[k]; xf = CDFargs[k+1]; yi = xp; yf = gCDF[k+1];}
//		  if (x < xp) { xi = CDFargs[k-1]; xf = CDFargs[k]; yi = gCDF[k-1]; yf = gCDF[k];}
//		  
//		  t[j] = (x-yi)*(xf-xi)/(yf-yi) + xi;
////		  t[j]=CDFargs[k]; // case for no interpolation
//	  }
//	  
//	  
//	  
//	}
//
//	return t;
//}

// ****************************************************************************************************
// THESE METHODS SHOULD BE IMPROVED CAUSE THEY ARE WRITTEN QUITE BADLY
void cpedsRNG::copyGCDF(const cpedsRNG& parent) {
	if (gCDF!=NULL) killGCDF();
	if (gCDF2d!=NULL) killGCDF2d();
	if (parent.CDFargs==NULL) { CDFargs=parent.CDFargs; }
	else {
		CDFargs=new double[parent.gCDFsize];
		for (long i=0;i<parent.gCDFsize;i++) { CDFargs[i]=parent.CDFargs[i]; }		
		gCDFsize=parent.gCDFsize;
	}
	if (parent.gCDF==NULL) { gCDF=parent.gCDF; gCDFsize=parent.gCDFsize; }
	else {
		gCDF=new double[parent.gCDFsize];
		for (long i=0;i<parent.gCDFsize;i++) { gCDF[i]=parent.gCDF[i]; }
		gCDFsize=parent.gCDFsize;
	}

//2d
	if (parent.CDFargsX==NULL) { CDFargsX=parent.CDFargsX; }
	else {
		CDFargsX=new double[parent.gCDFsizeX];
		for (long i=0;i<parent.gCDFsizeX;i++) { CDFargsX[i]=parent.CDFargsX[i]; }
		gCDFsizeX=parent.gCDFsizeX;
	}
	if (parent.CDFargsY==NULL) { CDFargsY=parent.CDFargsY; }
	else {
		CDFargsY=new double[parent.gCDFsizeY];
		for (long i=0;i<parent.gCDFsizeY;i++) { CDFargsY[i]=parent.CDFargsY[i]; }
		gCDFsizeY=parent.gCDFsizeY;
	}
	if (parent.gCDF2d==NULL) { gCDF2d=parent.gCDF2d; gCDFsize2d=parent.gCDFsize2d; }
	else {
		gCDF2d=new double[parent.gCDFsize2d];
		for (long i=0;i<parent.gCDFsize2d;i++) { gCDF2d[i]=parent.gCDF2d[i]; }
		gCDFsize2d=parent.gCDFsize2d;
	}


}
// ****************************************************************************************************
void cpedsRNG::killGCDF() {
	if (gCDF!=NULL) { delete gCDF; gCDF=NULL; }
	if (CDFargs!=NULL) { delete CDFargs; CDFargs=NULL; }
	gCDFsize=0;
}
/***************************************************************************************/
void cpedsRNG::killGCDF2d() {
	if (gCDF2d!=NULL) { delete gCDF2d; gCDF2d=NULL; }
	if (CDFargsX!=NULL) { delete CDFargsX; CDFargsX=NULL; }
	if (CDFargsY!=NULL) { delete CDFargsY; CDFargsY=NULL; }
	gCDFsize2d=0;
	gCDFsizeX=0;
	gCDFsizeY=0;
}
/***************************************************************************************/
long cpedsRNG::getPowerLawNoiseCoefficiensCount() const { 
	if (_distr==gaussian_power_law) return _b.size(); 
	if (_distr==gaussian_power_law_t) return _bt_size; 
	else return 0;	
}

// ****************************************************************************************************
void cpedsRNG::seedOffset(long seed_offset) {
	_seed_offset=seed_offset;
	gsl_rng_free(_state);
	_state=cpeds_random_uniform_numbersGSL_init(&_seed, seedOffset(),_generator_type); // this is the default one in gsl: gsl_rng_mt19937
}
// ****************************************************************************************************
void cpedsRNG::seed(long s) {
	_seed=s;
	gsl_rng_free(_state);
	_state=cpeds_random_uniform_numbersGSL_init(&_seed, seedOffset(),_generator_type); // this is the default one in gsl: gsl_rng_mt19937
}
// ****************************************************************************************************
double cpedsRNG::sumOOfNumbers(long n) {
//	printf("entering sum: n=%li\n",n);
	if (n==0) return 0;
	double s=0;
	n=cpeds_get_min(n,_b.count());
//	printf("entering sum: n=%li\n",n);
	for (long i = 1; i <= n; i++) {
		s+=_b[i]*_oofns[n-i];
//		printf("i: %li, n: %li, n-i: %li\n",i,n,n-i);
	}
//	printf("leaving sum s: %lE\n",s);
	return s;
}
// ****************************************************************************************************
double cpedsRNG::sumOOfNumbers_t(long n) {
//	printf("entering sum: n=%li\n",n);
	if (n==0) return 0;
	double s=0;
	n=cpeds_get_min(n,_bt_size);
//	printf("entering sum: n=%li\n",n);
	for (long i = 1; i <= n; i++) {
		s+=_bt[i]*_oofnst[n-i];
//		printf("i: %li, n: %li, n-i: %li\n",i,n,n-i);
	}
//	printf("leaving sum s: %lE\n",s);
	return s;
}
// ****************************************************************************************************
double cpedsRNG::sumOOfNumbers_t2(long n) {
	if (n==0) return 0;
	double s=0;
	long idx;
	n=cpeds_get_min(n,_bt_size);
	for (long i = 1; i <= n; i++) {
		idx=(n-i-_oofnst_stfrom) % _bt_size;
		if (idx<0) idx=_bt_size-idx;
		s+=_bt[i]*_oofnst[idx];
	}
	return s;
}
// ****************************************************************************************************
