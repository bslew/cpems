/*!
 * \file test-function-fft.cpp
 *
 *  Created on: Jul 22, 2010
 *      Author: blew
 */

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <tclap/CmdLine.h>
#include <fftw3.h>
#include <iomanip>
#include <iostream>
#include "Mscs-function.h"
#include "Mscs-function2dreg.h"
#include "Mscs-function3dregc.h"
#include "Mscs-map-window_function.h"
#include <QtCore/QList>
#include <QtCore/QString>
#include <QtCore/QStringList>
#include <assert.h>
#include "cpedsTime.h"
#include "cpeds-list.h"
#include "cpedsFunctionDistance.h"
#include <time.h>
#include <sys/resource.h>

#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
#endif
//using namespace MATPACK; // use all classes and functions in the MATPACK namespace
using namespace TCLAP;

typedef struct {
		long long i;
		double ch[12];
} ocrafDAQd_t;



// switches
bool _doPowerSpectrum=true;
bool  _mkNoise, _mkPhaseNoise, _mkGnoise, _mkSignal, _calcDistance1D, _calcDistance2D, _smooth3DGauss, _saveReIm,_nofftw, _shift0,_hann, _saveSlices, _save3Dpower, _calibrateByDuration, _psd;
bool _save3DboxSlices, _equalizePower, _printData, _saveData, _averageFns, _printStats, _fromTo;
bool _commentedFile;
cpedsList<double> _extract_f;
cpedsList<double> _extract_w;

// tests
bool _testName, _filter, _test_sin,  _test3DGaussSmooth,_test3DGaussSmooth2, _test1of,_testKasdin,_testGaussoofNormalized, _testCorrelation, _testParseval, _testPtSrcDistr, _testFFTcoefdistrchisq, _testFFTcoefdistrchisq2, _testPk3d;
double _convolve;
//bool _testModulation;


long _Nx,_Ny,_Nz, _dimension, _sliceAndAverage, _calibrateByStdev, _movingWindowStep, _corrStep, _corrNum;
long _skipFirst, _Nrows, _calculate1Dto3DpowerCalibFact;
int _colx,_coly, _colx2, _coly2, _ref1, _ref2, _save3DboxSlicesPlane;
string _outfile, _infile, _infile2, _format, _filterData, _prewhiten, _correlation, _outputFormat;
double _regridFactor, _Lx,_Ly,_Lz,_A,_k0,_ns,_seed,_kmin,_kmax,_dk,_samplingRate,_samplingRate2;
double _mkSignal_x1,_mkSignal_x2;
double _sxG,_syG,_szG;
mscsFunction _removeRanges;
vector<string> _inputFns;
string _mkSignalType;
mscsFunction _extractFromTo;









// tests
void mkTestParseval(int testNo);
void printStatistics(mscsFunction& f);
void averageFunctions(cpedsMsgs* msgs);
void save_run_config();
void parseOptions(int argc, char** argv);
void test3DGaussSmooth(cpedsMsgs* msgs);
void test3DGaussSmooth2(cpedsMsgs* msgs);
void loadInputData(mscsFunction& f,mscsFunction& f2, cpedsMsgs* msgs);
void calculateFunctionDistance1D(mscsFunction& f,mscsFunction& f2, cpedsMsgs* msgs);
void calculateFunctionDistance2D(mscsFunction& f,mscsFunction& f2, cpedsMsgs* msgs);
void loadDAQdFormatData(mscsFunction& f, cpedsMsgs* msgs);
void extractRanges(mscsFunction& re, mscsFunction& im, mscsFunction& P, cpedsMsgs* msgs);
//void removeRanges(mscsFunction& re, mscsFunction& im, mscsFunction& P, cpedsMsgs* msgs);

mscsFunction calculateCorrelationFunction(mscsFunction &f1, mscsFunction &f2,cpedsMsgs* msgs);
mscsFunction calculateDistance1D(mscsFunction &f1, mscsFunction &f2,cpedsMsgs* msgs);
void smooth3DGauss(mscsFunction3dregc& f, cpedsMsgs* msgs);


//void printUsage() {
//	printf("USAGE: function-fft XYfile regridFactor outfile\n\n where:\n\n regridFactor is 0 for uniform spacing assumption, and <1 for fraction of mean X points separations at which to linearly interpolate the function\n");	
//}

int main(int argc, char** argv) {
#ifdef DEBUG
	cpedsTime timer;
	timer.start(1);
	struct timeval startTime;
	struct timeval endTime;
	struct timeval startTime2;
	struct timeval endTime2;
	struct rusage ru;
	getrusage(RUSAGE_SELF, &ru);
	startTime = ru.ru_utime;
//
//			clock_t t0,t10;
//			t0=clock();
#endif
	cpedsMsgs* msgs = new cpedsMsgs("function_fft");
	msgs->setSaveRunWriteMode('a');
	msgs->saveThisRun(argc,argv);
	save_run_config();
	parseOptions(argc, argv);
	mscsFunction re,im,P;
	mscsWindowFunction wfn;

	if (_test_sin) {
		mscsFunction f;
		f.mkSin(1,10,10.0/_Nx,1,0);
		f.save("sin");
//		f.shift(_Ny);
//		f.save("sinShifted");
		exit(0);
	}
	
	if (_testGaussoofNormalized) {
		mscsFunction f("1/f"),g("gauss"), s("sum");
		mscsFunction fp("1/f"),gp("gauss"), sp("sum");
		
		f.mkPowerLawNoise(_A,_ns,_k0,_kmin,_kmax,_Nx,1,&fp);
		g.mkPowerLawNoise(1,0,_k0,_kmin,_kmax,_Nx,2,&gp);
		
		f/=f.stdev();
		g/=g.stdev();
		f.save(_outfile+".1of");
		//		fp.save(_outfile+".1of.p");
		g.save(_outfile+".gauss");
		//		gp.save(_outfile+".gauss.p");
		
		s=f+g;
		s.save(_outfile+".gauss1of");
		cpedsRNG rns;
		rns.seed(3);
		rns.setRNsType("uniform");
		rns.setMinMax(-1,1);
		double dx=s.getX(1)-s.getX(0);
		for (long i = 0; i < s.pointsCount(); i++) {
			s.X(i)=f.getx(i)+rns.getRN()*dx/3.0;
		}
		s.save(_outfile+".gauss1of0.3.pert");
		
		for (long i = 0; i < s.pointsCount(); i++) {
			s.X(i)=f.getx(i)+rns.getRN()*dx/2.0;
		}
		s.save(_outfile+".gauss1of0.5.pert");
		
		
		exit(0);
	}
	
	if (_testKasdin) {
		//
		// single run test
		//
		//		cpedsRNG* rns=new cpedsRNG("gaussianPowerLaw_t");
		//		rns->setPowerLawNoiseMinimalCoefficient(2e-8);
		//		rns->getRNs(68000000);//.save("1oftest3");
		//		rns->seed(1);
		//		rns->seedOffset(0);
		//		printf("Value of the minimal Kasdin coefficient: %lE\n",rns->getPowerLawNoiseMinimalCoefficient());
		//		printf("Number of Kasdin coefficients used: %li\n",rns->getPowerLawNoiseCoefficiensCount());
		//		delete rns;
		//		exit(0);
		
		// multiple runs with partial files save
		cpedsRNG* rns=new cpedsRNG("gaussianPowerLaw_t");
		rns->setPowerLawNoiseMinimalCoefficient(1e-8);
		rns->seed(1);
		rns->seedOffset(0);
		
		long N=70000000;
		long i=0;
		while (i<N) {
			rns->getRNs(1000000).save("1oftest4-70000000",false,"double",10,true);
			printf("Value of the minimal Kasdin coefficient: %lE\n",rns->getPowerLawNoiseMinimalCoefficient());
			printf("Number of Kasdin coefficients used: %li\n",rns->getPowerLawNoiseCoefficiensCount());
		}
		delete rns;
		exit(0);
		
		
	}
	if (_testName) {
		mscsFunction f("f"),c("copy");
		
		
		c=f.mkSin(0,1,0.0001,0.1);
		f=double(0);
		f.f(5000)=1;
		//		f.clearFunction();
		//		f.mkSin(0,1,0.0001,0.2);
		//		f+=c;
		//		for (long i=0;i<f.pointsCount();i++) { f.f(i)+=10; 	i+=500; }
		//		f.save("/home/blew/tmp/input_sin_signal");
		f.save("/home/blew/tmp/input_impulse_signal");
		
		//		f.extractSpikes(500,250,3);
		//		f.save("/home/blew/tmp/input_sin_signal-no-spikes");
		
		P=f.powerSpectrum(&re,&im,0);		
		//		P.save("/home/blew/tmp/input_sin_signal.power");
		P.save("/home/blew/tmp/input_impulse_signal.power");
		exit(0);
		f=c;
		cpedsList<double> w;
		cpedsList<long> b;
		int i=0,N=f.pointsCount();
		int bs=100;
		while (i<N) {
			b.append(bs); 
			//			w.append(1);
			i+=bs;
		}
		f.binFunction(long(0),b, w);
		f.save("/home/blew/tmp/sampled_signal");
		P=f.powerSpectrum(&re,&im,0);		
		P.save("/home/blew/tmp/sampled_signal.power");
		
		//		P=f.powerSpectrum(&re,&im,0);		
		//		P.save("/home/blew/tmp/test-function-fft-regrid_0.power");
		//		f=c;
		//		P=f.powerSpectrum(&re,&im,1);		
		//		P.save("/home/blew/tmp/test-function-fft-regrid_1.power");
		//		re.save("/home/blew/tmp/test-function-fft.re");
		//		im.save("/home/blew/tmp/test-function-fft.im");
		
		return 0;
	}
	
	if (_test1of) {
		mscsFunction f("1/f"),g("gauss"), s("sum");
		mscsFunction fp("1/f"),gp("gauss"), sp("sum");
		
		f.mkPowerLawNoise(_A,_ns,_k0,_kmin,_kmax,_Nx,0,&fp);
		g.mkPowerLawNoise(1,0,_k0,_kmin,_kmax,_Nx,2,&gp);
		
		double stdf,stdg;
		stdf=f.stdev();
		stdg=g.stdev();
		f/=stdf;
		g/=stdg;
		f.scaleX(1.0/_Nx);
		g.scaleX(1.0/_Nx);
		f.save(_outfile+"1of");
		fp.save(_outfile+"1of.p");
		g.save(_outfile+"gauss");
		gp.save(_outfile+"gauss.p");
		
		//		s=f+(g+double(1200));
		
		//		s.clearFunction();
		//		for (int i = 0; i < f.pointsCount(); i++) {
		//			s.newPoint(f.getX(i),f.getY(i)+g.getY(i));
		//		}
		//		s.save(_outfile+"1of_plus_gauss");
		
		//		s=(f+double(100*f.stdev()));
		//		s.save(_outfile+"1of_shift");
		
		//		s=(f+double(10*f.stdev()))*(g+double(1200));
		s=g*f;
		s.save(_outfile+"1of_times_gauss");
		sp=s.powerSpectrum(&re,&im,0);
		sp.save(_outfile+"1of_times_gauss.p");
		s.power(2);
		s.save(_outfile+"1of_times_gauss2");
		sp=s.powerSpectrum(&re,&im,0);
		sp.save(_outfile+"1of_times_gauss2.p");
		// simulate video amplifiers smoothing effect
		mscsWindowFunction w;
		//		w.mkExpWindow(sp.getMinArg(),sp.getMaxArg(),(sp.getMaxArg()-sp.getMinArg())/(sp.pointsCount()-1),sp.getX(long(sp.pointsCount())/2),1);
		w.mkExpWindow(sp.getMinArg(),sp.getMaxArg(),_Nx,1500,1000);
		w.save("window.txt");
		s=s.convolve_window(w,false);
		s.save(_outfile+"1of_times_gauss2-video");
		sp=s.powerSpectrum(&re,&im,0);
		sp.save(_outfile+"1of_times_gauss2-video.p");
		
		//		g*=g;
		//		g.save(_outfile+"gauss2");
		
		return 0;
	}
	
	
	if (_testPk3d) {
		mscsFunction3dregc f(_Nx,_Ny,_Nz, _Lx/_Nx,_Ly/_Ny,_Lz/_Nz);
		mscsFunction Pk;
		//		Pk.mkPowerLaw(_kmin,_kmax,1.0/(_Lx/_Nx),_A,_k0,_ns);
		//		Pk.mkPowerLaw(_kmin,_kmax,(_kmax-_kmin)/10000,_A,_k0,_ns);
		Pk.mkGauss(_kmin,_kmax,(_kmax-_kmin)/10000,1,0,.10);
		Pk.printRanges();
		Pk.save("Pk");
		f.printInfo();
		//		exit(0);
		f.generateRandomGaussianField3(Pk);
		msgs->say("saving the files...",High);
		cpeds_matrix_save(f.getSlice(0,50,0),"slice1");
		matrix<double> tmpm=f.getSlice(0,0,0);
		cpeds_shift_matrix(tmpm,_Nx/2,_Ny/2);
		cpeds_matrix_save(tmpm,"slice1-shift");
		cpeds_matrix_save(f.getSlice(1,50,0),"slice2");
		cpeds_matrix_save(f.getSlice(1,50,1),"slice2.im");
		cpeds_matrix_save(f.getSlice(2,50,0),"slice3");
		cpeds_matrix_save(f.getSlice(2,50,1),"slice3.im");
		cpeds_matrix_save(f.getSlice(2,20,0),"slice4");
		msgs->say("done",High);
		
		//		char tmpch[100];
		//		for (long i = 0; i < _Nx; i++) {
		//			msgs->say("saving slice:"+msgs->toStr(i),Low);
		//			sprintf(tmpch,"all-slices/slice0i0.%03li",i);
		//			cpeds_matrix_save(f.getSlice(0,i,0),string(tmpch));			
		//		}
		
		//		FILE* tmpf=fopen("raw","w");
		//		for (long i = 0; i < _Nx*_Ny; i++) {
		//			fprintf(tmpf,"%lE ", f.fRe(i));
		//		}
		//		fclose(tmpf);
		
		//		f.fft_c2c(false);
		msgs->say("saving all",High);
		//		f.savetxtlin("all");
		msgs->say("done",High);
		exit(0);
	}

	
	if (_testParseval) {
		mkTestParseval(0);
		mkTestParseval(1);
		mkTestParseval(3);
		mkTestParseval(5);
		exit(0);		
	}
	
	
	if (_testPtSrcDistr) {
		
		mscsFunction f("1/f");
		mscsFunction fp("1/f");
		
		//		f.mkPowerLawNoise(_A,_ns,_k0,_kmin,_kmax,_Nx,0,&fp);
		f.mkChisqNoise(_Nx,0,1);
		
		//		f/=f.stdev();
		f.save(_outfile+".ptsrcDistrText");
		
		// correct for the distance from the beginning of the coordinate system
		long n = f.pointsCount();
		for (long i = 0; i < n; i++) {
			f.f(i)/=(f.getX(i)+1.0)*(f.getX(i)+1.0);
		}
		f.save(_outfile+".ptsrcDistrText.obs");
		
		// randomize the generated fluxes so that they have random position wrt observer
		
		// calculate the power spectrum of this
		mscsFunction re,im;
		fp=f.powerSpectrum(&re,&im,0);
		
		fp.save(_outfile+".ptsrcDistrText.obs.p");
		exit(0);
		
	}
	
	if (_testFFTcoefdistrchisq) {
		mscsFunction f("1/f"),re,im;
		mscsFunction fp("1/f");
		
		f.mkChisqNoise(_Nx,0,1);
		//		f-=f.meanf();
		f.save(_outfile+"_testFFTcoefdistrchisq.signal");
		fp=f.powerSpectrum(&re,&im,0);
		re.deletePoint(0);
		im.deletePoint(0);
		fp.deletePoint(0);
		re.save(_outfile+"_testFFTcoefdistrchisq.re");
		im.save(_outfile+"_testFFTcoefdistrchisq.im");
		fp.save(_outfile+"_testFFTcoefdistrchisq.p");
		
		exit(0);
	}
	
	if (_testFFTcoefdistrchisq2) {
		mscsFunction f("1/f");
		mscsFunction fp("1/f");
		mscsFunction *re=NULL,*im=NULL;
		mscsFunction vars, varav,bins;
		long Nsim=1000;
		cpedsList<long> bs;
		cpedsList<double> w;
		cpedsList<double> dcovL, varInBin;
		double *dcov,*cov, *x,*y;
		long bsini=10;
		double geogm=1.1;
		
		double spidx=_ns;
		
		for (long i = 0; i < Nsim; i++) {
			msgs->say("SIMULATION: "+msgs->toStr(i),Top);
			bs.clear();
			fp.clearFunction();
			fp.mkPowerLawNoise(1,spidx,1,0.1,1000,100000,0,NULL,NULL,false);
			fp.binFunctionLin2geo(0,bsini,geogm,bs,w,&varInBin);
			
			// store and average the variances in the bins
			y=varInBin.toCarray();			vars.importFunction(NULL,y,varInBin.size()); delete [] y;
			if (i==0) { varav=vars; } else { varav+=vars; }
			vars.clearFunction();
			varInBin.clear();
			
			f.clearFunction();
			f.setX(fp.toXList());
			f.mkPowerLawNoise(1,spidx,1,0,&fp);
			
			// form a dcov of the binned spectra to calculate the covariance matrix
			// add the sparsely generated power spectrum values to dcovv
			y=fp.extractValues();			dcovL.fromCarray(y,fp.pointsCount()); delete [] y;
			f.clearFunction();
		}
		// calculate average variance per bin
		varav/=double(Nsim);
		// save dcov, bins and averaged variances in bins
		//		dcovv.save(_outfile+"_testFFTcoefdistrchisq2.dcov.vecSize"+msgs->toStr(fp.pointsCount())+"-Nsim_"+msgs->toStr(Nsim));
		for (long i = 0; i < fp.pointsCount(); i++) {			
			bins.newPoint(fp.X(i),bs[i]);
			varav.setarg(i,fp.X(i));
		}
		bins.save(_outfile+"_binSize_geo-bsini_"+msgs->toStrf(bsini,0)+"-gm_"+msgs->toStrf(geogm,2));
		varav.save(_outfile+"_highres.varAV-alpha_"+msgs->toStrf(spidx,1));
		
		// calculate the variance on the bins from the power spectra generated on the bins
		dcov=dcovL.toCarray();
		cov=cpeds_calculate_covariance_matrix(dcov,fp.pointsCount(),Nsim,true);
		delete [] dcov;
		
		vars.clearFunction();
		vars=varav;
		// multiply the variances calculated from the power spectra generated on the bins by the size of the bins to test against the variances in the bins calculated on the non-binned spectr
		for (long i = 0; i < vars.pointsCount(); i++) {			vars.f(i)=cov[i]; }//*bs[i];		}
		vars.save(_outfile+"_cov-alpha_"+msgs->toStrf(spidx,1));
		//		cpeds_save_matrix(cov,fp.pointsCount(),1,"cov.highres",true);
		delete [] cov;
		
		exit(0);
	}
	
	
	
	/***************************************************************************************/	
	// MAIN PART OF THE PROGRAM	
	/***************************************************************************************/	
	
	if (_averageFns) {
		averageFunctions(msgs);
		return 0;
	}
	
	
	if (_dimension==1) {
		mscsFunction f("f"),c("copy");
		mscsFunction f2("f2");
		mscsFunction re,im,P,Pav, Pcal, Pcalint;
		char tmpch[1000];
		long sliceSize;
		
		if (_mkPhaseNoise) {
			f.mkPhaseNoise(_Nx,_seed);
			f.save(_outfile);
			exit(0);
		}
		
		if (_mkNoise) { 
			msgs->say("Generating noise data...",High);
			f.mkPowerLawNoise(_A,_ns,_k0,_kmin,_kmax,_Nx,_seed);
			msgs->say("Generating noise data done",Medium);
			msgs->say("Saving noise data...",High);
			string tmps; if (_outfile=="") tmps="signal.txt"; else tmps="-signal.txt";
			f.save(_outfile+tmps);
			printf("stdev: %lE\n",f.stdev());
			msgs->say("Generating noise data...done",Medium);
		}
		else {
			if (_mkSignal) {
				if (_mkSignalType=="topHat" or _mkSignalType=="topHatNoise") {
					msgs->say("Generating signal data...",High);
					f.mkTopHat(-_mkSignal_x2,_mkSignal_x2,(2*_mkSignal_x2)/_Nx,0,-_mkSignal_x1,_mkSignal_x1,1);
				}
				if (_mkSignalType=="topHatNoise") {
					cpedsRNG rns("uniform");
					rns.setMinMax(-0.1,0.1);
					for (long i = 0; i < f.pointsCount(); i++) {
						if (f.getx(i)>-_mkSignal_x1 and f.getx(i)<_mkSignal_x1) {
							f.f(i)+=rns.getRN();
						}
					}
				}
				
				if (_mkSignalType=="sineSineModulated") {
					msgs->say("Generating signal data: modulated sine function",High);
					f.mkSin(0,10,0.0001,1.0/50,0);
					mscsFunction fmod;
					fmod.mkSin(0,10,0.0001,1.0,0);
//					fmod*=double(0.4);
					fmod+=double(1);
					f*=fmod;
				}
				if (_mkSignalType=="sineSqModulated") {
					msgs->say("Generating signal data: modulated sine function",High);
					f.mkSin(0,10,0.0001,1.0/50,0);
					mscsFunction fmod;
//					fmod.mkGaussianNoise(f.pointsCount(),0,0.1,0);
//					f+=fmod;
					double Mamp=0.2;
					double famp=2;
					fmod.mkSquareWave(0,10,0.0001,1.0/famp,0);
					fmod*=double(0.2);
					fmod+=double(1.0-Mamp);
					f*=fmod;
				}
				if (_mkSignalType=="twoSines") {
					msgs->say("Generating signal data: multiple sines function added together",High);
					mscsFunction f2;
					f.mkSin(0,_Lx,1.0/_samplingRate,1.0/50,0.12);
					
					f2.clearFunction();
					f2.mkSin(0,_Lx,1.0/_samplingRate,1.0/50.01,0.12);
					f2*=double(0.5);
					f+=f2;					

					f2.clearFunction();
					f2.mkSin(0,_Lx,1.0/_samplingRate,1.0/100,0.6);
					f2*=double(0.5);
					f+=f2;					
					
					f2.clearFunction();
					f2.mkSin(0,_Lx,1.0/_samplingRate,1.0/150,0.87);
					f2*=double(0.3);
					f+=f2;					
				}
				if (QString(_mkSignalType.c_str()).split(':')[0]=="sinesFromDphiFile") {
					string fname=QString(_mkSignalType.c_str()).split(':')[1].toStdString();
					msgs->say("Generating signal data from D,phi information stored in file.",High);
					mscsFunction D,phi,tmp,tmp2;
					D.load(fname,false,0,1);
					phi.load(fname,false,0,2);
					for (long i = 0; i < D.pointsCount(); i++) {
						cout << "simulating sine: " << i << " / " << D.pointsCount() << "\r";
						tmp.clearFunction();
						tmp.mkSin(0,_Lx,1.0/_samplingRate,1.0/D.X(i),PIsnd-phi.f(i)); // this is NOT doing what you think it suppose to.
						tmp*=double(D.f(i));
						if (i==0) f=tmp; else f+=tmp;
					}
				}
				if (QString(_mkSignalType.c_str()).split(':')[0]=="signalFromReImFile") {
					string fname=QString(_mkSignalType.c_str()).split(':')[1].toStdString();
					msgs->say("Generating signal data from re,im information stored in file.",High);
					mscsFunction re,im,tmp,tmp2;
					re.load(fname,false,0,1);
					im.load(fname,false,0,2);
					for (long i = 0; i < re.pointsCount(); i++) {
						cout << "simulating component: " << i << " / " << re.pointsCount() << "\r";
						tmp.clearFunction();
						tmp.mkSin( 0,_Lx,1.0/_samplingRate,1.0/re.X(i),0.25,re.f(i));
						tmp2.clearFunction();
						tmp2.mkSin(0,_Lx,1.0/_samplingRate,1.0/re.X(i),0,im.f(i));
						tmp-=tmp2;
						if (i==0) f=tmp; else f+=tmp;
					}
					if (_hann) { // EXPERIMENTALL 
//						wfn.mkHannWindow(_Ly*_samplingRate);
//						wfn=wfn.cut(f.pointsCount(),0);
//						wfn.print();
//						f/=wfn;
						f*=double(4.0); // this factor comes from tests; I'm not sure exactly where does it come from - probably some DFT vs fft definition differences						
					}
					else {
						f*=double(2.0); // this factor comes from tests; I'm not sure exactly where does it come from - probably some DFT vs fft definition differences						
					}
				}

				msgs->say("Saving signal data...",High);
				string tmps; 
				if (_infile=="") {
					if (_outfile!="signal.txt") tmps=_outfile+".generated_signal"; else tmps="input_signal.txt";
				}
				else {					tmps=_infile;				}
				f.save(tmps);
			}
			else {
				//
				// load input data
				//
				msgs->say("Loading input data",High);
				loadInputData(f,f2,msgs);
				msgs->say("Loading done",Medium);				
			}
		}
		
		if (_printStats) { printStatistics(f); }
		
		//
		// extra stuff Begin
		//
		if (_calcDistance1D) { calculateFunctionDistance1D(f,f2,msgs); }
		if (_calcDistance2D) { calculateFunctionDistance2D(f,f2,msgs); }
		
		
		if (_correlation!="") { 
			mscsFunction corr=calculateCorrelationFunction(f,f2,msgs); 
			corr.setName("correlationFunction");
			corr.save(_outfile);
			exit(0);
		}
		
		if (_convolve>0) {
			f.sortFunctionArgAscending();
			mscsWindowFunction w=f;
//			w.mkGauss(0,w.getMaxArg(),w.getX(1)-w.getX(0),1,0,_convolve);
			w.mkGauss(0,0,0,1,0,_convolve);
			f=f.convolve_window(w,true);

//			w.make_gaussian_kernel_kspace(_convolve);
//			double T=f.getMaxArg();
//			w.mkExpWindow(2.0/T,f.pointsCount()/T,f.pointsCount()/2,2.0/T,_convolve);
			w.save("w");
//			f=f.convolve_window(w,false);
			f.save(_outfile);
			exit(0);
		}
		
		//
		// extra stuff End
		//
				

		//
		// FFT stuff begin
		//
		
		// calibrate by standard deviation
		if (_calibrateByStdev>3) {
			f.calibrateByStDev(_calibrateByStdev,_movingWindowStep);
		}
		
		if (_psd) {
			_sliceAndAverage=long(sqrt(f.pointsCount()));
			_saveSlices=false;
		}
		
		if (_prewhiten!="") {
			Pcal.load(_prewhiten);
		}
		
		sliceSize=f.pointsCount()/_sliceAndAverage;
		c=f.get(0,sliceSize*_sliceAndAverage);
		msgs->say("Calculating FFT",High);
		msgs->say(string("Will average from ")+msgs->toStr(_sliceAndAverage)+" slices.",Low);
//		f.print();
		for (long i = 0; i < _sliceAndAverage; i++) {
			msgs->say("Calculating FFT on slice: "+msgs->toStr(i),Medium);
			// prepare slice
			f=c.get(i*sliceSize,(i+1)*sliceSize);
			
			if (_shift0) f-=f.meanf();
			if (_hann) {
				wfn.mkHannWindow(f.pointsCount());
				f*=wfn;
			}
			
			if (_saveSlices) { f.save(_outfile+".slice_"+msgs->toStr(i)); }
			// calculate power spectrum
#ifdef DEBUG
			getrusage(RUSAGE_SELF, &ru);
			startTime2 = ru.ru_utime;
//			clock_t t1,t2;
//			t1=clock();
#endif
			if (_nofftw) {
				P=f.powerSpectrum(_kmin,_kmax,_dk,&re,&im);
			}
			else {
				P=f.powerSpectrum(&re,&im,_regridFactor);
				//			P=f.powerSpectrum(NULL,NULL,_regridFactor);
			}
			
#ifdef DEBUG
			getrusage(RUSAGE_SELF, &ru);
			endTime2 = ru.ru_utime;
			getrusage(RUSAGE_SELF, &ru);
			double tS = startTime2.tv_sec*1000000 + (startTime2.tv_usec);
			double tE = endTime2.tv_sec*1000000  + (endTime2.tv_usec);
			printf("power spectrum exec time [s]: %lf\n",(tE-tS)*1e-6);
//			t2=clock();
//			printf("power spectrum exec time [s]: %lf\n",double(t2-t1)/CLOCKS_PER_SEC);
#endif

			extractRanges(re,im,P,msgs);
//			removeRanges(re,im,P,msgs);
			
			if (_calibrateByDuration) { P/= (f.pointsCount()/_samplingRate); } // here we divide by sampling period but this is now done (by default) inside of the power spectrum routine
			if (_saveSlices) { P.save(_outfile+".slice_"+msgs->toStr(i)+".p"); }
			
			// average power spectra
			if (i==0) Pav=P;
			else Pav+=P;
			
			if (_prewhiten!="") {
				msgs->say("Prewhitening",High);
				Pcalint=Pcal.interpolate(P.extractArguments(),P.pointsCount(),"linear",false);
				P/=Pcalint;
				Pcalint.sqroot();
				re/=Pcalint;
				im/=Pcalint;
			}
			if (_equalizePower) {
				msgs->say("equalizing power",High);
				c=re.power(2)+im.power(2);
				c.sqroot();
				re/=c;
				im/=c;
			}
			
			if (_filter) {
				msgs->say("filtering",High);
				mscsFunction filter("filter");
				if (_filterData!="") {
					QString qs=_filterData.c_str();
					bool ok;
					QStringList qsl=qs.split(',');
					if (qsl.size()>1) {
						for (long i = 0; i < qsl.size(); i++) {
							filter.newPoint( qsl[i].toDouble(&ok)-_extract_w[i % _extract_w.size()], qsl[i].toDouble(&ok)+_extract_w[i % _extract_w.size()] );
							assert(ok==true);
						}
					}
					else {
						filter.load(_filterData);						
					}
					filter.print();
					for (long i = 0; i < filter.pointsCount(); i++) {
						re.setRange(filter.getX(i),filter.getY(i),0);
						im.setRange(filter.getX(i),filter.getY(i),0);
					}
				}
				msgs->say("inverse fft",High);
				f.inverseFFT(re,im,f.extractArguments(),true);
				msgs->say("saving filtered data",High);
				if (_saveSlices) { f.save(_outfile+".slice_"+msgs->toStr(i)+".sig"); }
//				if (_filterData!="") {	f.save(_outfile+".filtered--"+_filterData); }
//				else 
					f.save(_outfile+".filtSig");
			}
					
		}
		if (_sliceAndAverage>1) {
			P=Pav/double(_sliceAndAverage);
		}
		msgs->say("FFT done.",Medium);
		//			if (_calibrateByDuration) { P/= (f.pointsCount()/_samplingRate); }
		msgs->say("Saving data...",High);
		P.save(_outfile);
		if (_saveReIm) {
			if (_outputFormat=="reim2") {
				re.save(_outfile+".re");
				im.save(_outfile+".im");				
			}
			if (_outputFormat=="reim2nonzero") {
				re.removeValue(0);
				re.save(_outfile+".renz");
				im.removeValue(0);
				im.save(_outfile+".imnz");				
			}
			if (_outputFormat=="reim1nonzero") {
				long n=re.pointsCount()-1;
				for (long i = n; i >= 0; i--) {
					if (re.f(i)==0 and im.f(i)==0) { re.deletePoint(i); im.deletePoint(i); }				}
//				re.removeValue(0);
//				im.removeValue(0);
								
				cpedsList<double> l;
				l.append(re.getArguments());
				l.append(re.getValues());
				l.append(im.getValues());
				double* p=l.toCarray();
				matrix<double> m=cpeds_array2matrix(p,l.size(),re.pointsCount(),false);
				cpeds_matrix_save(m,_outfile+".reimnz");
				delete p;
			}

			if (_outputFormat=="Dphi2nonzero") {
				re.removeValue(0);
				im.removeValue(0);
				
//				mscsFunction D,phi;				
//				D=re;
//				D.invert();
//				D.setY(im.getValues());
//				D.convertReImToPhiR();
//				phi=D;
//				phi.invert();
//				phi.setX(re.getArguments());
//				D.setX(re.getArguments());
//				
//				D.save(_outfile+".Dnz");
//				phi.save(_outfile+".phinz");				
				
				
				
				re.setX(im.getValues()); // re: im, re
				re.invert(); // re: re, im
				re.convertReImToPhiR(); // re: phi, R
				im.setY(re.getArguments()); // im: f, phi
				re.setX(im.getArguments());
				
				re.save(_outfile+".Dnz");
				im.save(_outfile+".phinz");				
			}
			if (_outputFormat=="Dphi1nonzero") {
				long n=re.pointsCount()-1;
				for (long i = n; i >= 0; i--) {
					if (re.f(i)==0 and im.f(i)==0) { re.deletePoint(i); im.deletePoint(i); }				}
//				re.removeValue(0);
//				im.removeValue(0);
								
				re.setX(im.getValues()); // re: im, re
				re.invert(); // re: re, im
				re.convertReImToPhiR(); // re: phi, R
				cpedsList<double> l;
				l.append(im.getArguments());
				l.append(re.getValues());
				l.append(re.getArguments());
				double* p=l.toCarray();
				matrix<double> m=cpeds_array2matrix(p,l.size(),re.pointsCount(),false);
				cpeds_matrix_save(m,_outfile+".Dphinz");
				delete p;
			}
			
		}
		msgs->say("Saving done",Medium);
	}
	
	if (_dimension==2) {
		mscsFunction2dreg f(_Nx,_Ny, _Lx/_Nx,_Ly/_Ny);
		f.clearFunction();
		//		f.mkSin2Drad(0,1,0.1,1,0);
		//		f.mkSin2D(0,1,0.05000001,1,0);
		
		//		f.mkGaussianNoise(16,0,1,1,NULL,true,false);
		//		f.save("sig");
		//		//		f.invert();
		//		//		f.mscsFunction::setf(double(0.0));
		//
		//		f.saveAsMatrix("sig.real");
		//		f.fft(true);
		//		f.saveAsMatrix("sig.kspace");		
		//		f.antisymmetrize();		
		//		f.saveAsMatrix("sig.kspace.check");		
		//		f.fft(false);
		//		f.saveAsMatrix("sig.real.check");		
		//		exit(0);
		
		f.mkGaussianNoise(16,0,1,1,NULL,true,true);
		f.saveAsMatrix("sig.kspace");
		
		f.antysymmetrize();		
		f.saveAsMatrix("sig.kspace.asym");
		//		exit(0);
		//		f.printInfo();
		f.fft(false);
		f.saveAsMatrix("sigreal");
		exit(0);
		//		f.fft(true);
		//		f.saveAsMatrix("check");
		
		//		cpeds_matrix_save(f.getSlice(0,0,0),"slice");
		mscsFunction Pk=f.powerSpectrum();
		Pk.save("power_spectrum.tmp");
	}
	
	if (_dimension==3) {
		if (_test3DGaussSmooth) { test3DGaussSmooth(msgs); exit(0); }
		if (_test3DGaussSmooth2) { test3DGaussSmooth2(msgs); exit(0); }
		
		mscsFunction3dregc f(_Nx,_Ny,_Nz, _Lx/_Nx,_Ly/_Ny,_Lz/_Nz);
		if (_infile!="") {	f.load(_infile); }
		else {
			if (_mkNoise) {
				msgs->say("Generating power law noise field in 3D box",High);
				_doPowerSpectrum=false;		if (_save3Dpower) _doPowerSpectrum=true;
				
				mscsFunction Pk("Pk");
				if (_calculate1Dto3DpowerCalibFact) {
					Pk.mkPowerLaw(_kmin, sqrt(3.0)*_kmax,_dk, _A,_k0,_ns);
					Pk.save("Pk1d");
					_ns-=2;
					Pk.clearFunction();
				}
				Pk.mkPowerLaw(_kmin, sqrt(3.0)*_kmax,_dk, _A,_k0,_ns);
				Pk.save("Pk3d");					

				f.generateRandomGaussianField(Pk,_seed);
				if (_calculate1Dto3DpowerCalibFact>0) {
					mscsFunction slice1d,sliceAvP;
					cpedsRNG rnsY, rnsZ;
					rnsY.setMinMax(0,f.Nx()); rnsZ.setMinMax(0,f.Nz()); rnsY.seedOffset(1); rnsZ.seedOffset(2);
					
					slice1d.clearFunction();
					sliceAvP.clearFunction();
					double P1dAv=0, calibFact;
					for (long i = 0; i < _calculate1Dto3DpowerCalibFact; i++) {
						slice1d=f.getSlice1Dfn(0,long(rnsY.getRN()),long(rnsZ.getRN()),0);
						Pk=slice1d.powerSpectrum();
						if (i==0) sliceAvP=Pk; else sliceAvP+=Pk;
						P1dAv+=Pk.f(_k0,NULL);
					}
					P1dAv/=_calculate1Dto3DpowerCalibFact;
					sliceAvP/=_calculate1Dto3DpowerCalibFact;
					calibFact=sqrt(P1dAv/_A); 
					calibFact=P1dAv/_A; 
												/*
												 * Comment: This is an arbitrary choice based on experiments.
												 * Eventually this calibration should be done in terms of fitting the amplitude
												 * of the reconstructed average power spectra to the desired one.
												 * Or ideally via deconvolution methods to mitigate the effects of excessive
												 * power leaking to large scales and suppressed power in small scales.
												 * 
												 * author: blew
												 * date: Jun 26, 2013 12:59:49 PM
												 *
												 */
					
					sliceAvP.save("1DsliceAv.p");
					sliceAvP/=calibFact;
					sliceAvP.save("1DsliceAv.corrected.p");
					msgs->setLogFileName("readme");
					msgs->loggingOn();
					msgs->say("1D to 3D calibration factor: %lE",calibFact,Medium);
					f/=sqrt(calibFact);
					

					// re-test

									slice1d.clearFunction();
									sliceAvP.clearFunction();
									for (long i = 0; i < _calculate1Dto3DpowerCalibFact; i++) {
										slice1d=f.getSlice1Dfn(0,long(rnsY.getRN()),long(rnsZ.getRN()),0);
										Pk=slice1d.powerSpectrum();
										if (i==0) sliceAvP=Pk; else sliceAvP+=Pk;
										P1dAv+=Pk.f(_k0,NULL);
									}
									P1dAv/=_calculate1Dto3DpowerCalibFact;
									sliceAvP/=_calculate1Dto3DpowerCalibFact;
									sliceAvP.save("1DsliceAv.corrected.reconstructed.p");
									
					msgs->say("This is a WV field simulation run. The principal files are:",Medium);
					msgs->say("1DsliceAv.p - 1D power averaged spectrum from %li draws from the generated 3D field before amplitude corrections",_calculate1Dto3DpowerCalibFact,Low);
					msgs->say("1DsliceAv.corrected.p - same as 1DsliceAv.p but scaled by calibration factor: %lE, which is power ratio at k0 of what was derived with what was wanted",calibFact,Low);
					msgs->say("1DsliceAv.corrected.reconstructed.p - 1D power averaged spectrum from %li draws from the generated 3D field after amplitude corrections to the field in real space by factor sqrt(fact): %lE",sqrt(calibFact),Low);
					msgs->say("Pk1D - 1D input model power spectrum",Low);
					msgs->say("Pk3D - 3D input model power spectrum steeper by 2 than the input 1D power spectrum",Low);
					msgs->say("The current approach is not optimal, due to 3D power distortions which lead to suppression of reconstructed 1D power wrt model at small scales and "
							"boosting power at large scales. Some sort of fourier space filtering could be applied to mitigate this effect."
							"Currently the WV fluctuations at small scales will be underestimated. Use:\n\n"
							"plot_function.py --logX --logY Pk1d Pk3d 1DsliceAv.p 1DsliceAv.corrected.p 1DsliceAv.corrected.reconstructed.p --mkLabels "
							"\n\nto plot.",Medium);
					msgs->loggingOff();
				}
			}
			if (_mkGnoise) {
				msgs->say("Generating gaussian noise field in 3D box",High);
				_doPowerSpectrum=false;
				f.generateRandomGaussianField(0,1,true,0,_seed,0);
				f.saveHDF5(_outfile+".hdf5","gnoise");

			}
			if (_mkSignal) {
				if (_mkSignalType=="gauss3d") {
					f.mkGauss3D(_sxG,_syG,_szG);
				}
			}

		}
		if (_smooth3DGauss) { 
			_doPowerSpectrum=false;
			smooth3DGauss(f,msgs); 
		}
		if (_outputFormat=="df3") {
			msgs->say("Saving 3D box to a df3 files: "+_outfile,High);
			f.saveDF3(_outfile+"-re.df3",0);
			f.saveDF3(_outfile+"-im.df3",1);
		}
		if (_outputFormat=="raw") {
			msgs->say("Saving 3D box to a raw file: "+_outfile+".raw",High);
			f.save(_outfile+".raw");			
		}
#ifndef NO_HDF5
		if (_outputFormat=="hdf5") {
			msgs->say("Saving 3D box to a hdf5 file: "+_outfile+".hdf5",High);
			f.saveHDF5(_outfile+".hdf5","field");
		}
#endif
		if (_save3DboxSlices) { 
			msgs->say("Saving all slices of the 3D box directory with slices: "+_outfile,High);
			f.saveAllSlices(_outfile,_outfile,_save3DboxSlicesPlane,0); 
		}
		
		if (_doPowerSpectrum) {	
			mscsFunction Pk=f.powerSpectrum(); 
			Pk.save(_outfile+".raw");				
		}
	}
	
	
	//		f.mkSin3D(0,1,0.05,0.5,0);
	//			f.generateRandomGaussianField(0,1,true,false,1,0);
	//		f.antisymmetrize();
	
	//		f.save("sin3d");
	//		f.printInfo();
	//		cpeds_matrix_save(f.getSlice(0,0,0),"slice");
	
	//			f.fft(true);
	//		exit(0);
	//		f.print();
	//			f.savetxtlin("field.re");

#ifdef DEBUG
//			t10=clock();
//			printf("function_fft exec time [s]: %lf\n",double(t10-t0)/CLOCKS_PER_SEC);
			getrusage(RUSAGE_SELF, &ru);
			endTime = ru.ru_utime;
			getrusage(RUSAGE_SELF, &ru);
			double tS = startTime.tv_sec*1000000 + (startTime.tv_usec);
			double tE = endTime.tv_sec*1000000  + (endTime.tv_usec);
//			printf("function_fft exec time [s]: %lf\n",(tE-tS)*1e-6);
			timer.stop(1);
			printf("function_fft exec time [s]: %lf\n",timer.timeDiff(1));

#endif
	return 0;

}















void parseOptions(int argc, char** argv) {
	long i;
	string::size_type j;
	char tmpch[1000];
	
	try {
		
		CmdLine cmd("function_fft: calculate the fft on an input data. Can also do other stuff...\n", ' ', "" );
		
		//
		// Define arguments
		//
		
		UnlabeledMultiArg<string> inputFns("inputFunctions","files to be used for various operations unrelated to fft: eg. averaging", false,"string"); cmd.add(inputFns);
		
		ValueArg<string> input_file("","inf","file XY[Z] function file. The XY format is for 1D FFT of a function like structure. "
				"The spacings of the points may not be equal."
				"For the case of dims>1 the format is one of those described in the format parameter.",
				false,"","string");	cmd.add( input_file );
		
		ValueArg<string> input_file2("","inf2","file Y-values function file. This for for 1D type analyses (1d fft, cross-correlation analysis)"
				"The spacings of the points may not be equal. If this file is given then the program will use it to read the y values from it for "
				"the analysis using the column defined by --coly parameter while the file given in --inf option will be used for reading the x values"
				"from the column indicated by the --colx option",
				false,"","string");	cmd.add( input_file2 );
		
		//		ValueArg<string> nofile("","","use if no input file is needed",false,"","string"); cmd.add(nofile);
		ValueArg<string> output("o","out","outfile ",false,"out","string"); cmd.add(output);
		
		vector<string> allowedFormats;
		allowedFormats.push_back("1col");
//		allowedFormats.push_back("2col");
		allowedFormats.push_back("ascii");
		allowedFormats.push_back("binDAQd");
		ValuesConstraint<string>* allowedFormatsNew = new ValuesConstraint<string>(allowedFormats);
		
		ValueArg<string> format("f","input_format","format of the input data for the n>1 case\n"
				"The following formats are available for the n>1 case:\n\n"
				"* 1col - a set of numbers is read from the input file until the end of the file and an n-dimensional array is filled with values"
				"in the X-major ordering: i.e. n=Nz*Ny*i+Nz*j+k where {i,j,k} in 0..{Nx,Ny,Nz} for 3d case (currently implemented only this one)\n"
				"Parameters Nx,Ny,Nz - define the size of the array for each dimension. As a result the fft is performed on an equally spaced points. No interpolation is performed in this case\n"
					"* 2col - a set of doublets is read as in 1col case but in this case the function is assumed to be complex-valued.NOT YET IMPLEMENTED.\n\n"
				"The following file formats are available for the 1D case:\n"
				"* ascii - default format for reading ordinary ascii files formatted in columns of data"
				"* binDAQd - binary file format with one row defined as: long long, 11xdouble. With such defined row, it is possible to"
				"use the colx and coly parameters to choose the columns that are required. The coly parameter cannot be set to 0. The column number identifies"
				"the column of data stored in the binary file. The columns in the binary files are:\n"
				"i, time, ch1,..., ch8, timeCal, switch signal, calibration signal, so that column 1 (coly 1) identifies the time column and coly 2 identifies "
				"the first data channel."
				,false,"ascii",allowedFormatsNew); cmd.add(format);
		
		vector<long> allowedDims;
		allowedDims.push_back(1);
		allowedDims.push_back(2);
		allowedDims.push_back(3);
		
		ValuesConstraint<long>* allowedDimsNew = new ValuesConstraint<long>(allowedDims);
		ValueArg<long> dimension("n","dims","number of dimensions (default:1)",false,1,allowedDimsNew); cmd.add(dimension);
		
		//		MultiArg<long> N("N","boxSize","number of points. Defines Nx when single number is given. Provided for convinience ",false,"long"); cmd.add(N);
		ValueArg<long> Nx("","Nx","size of grid in n>1 case in X direction ",false,1,"long"); cmd.add(Nx);
		ValueArg<long> Ny("","Ny","size of grid in n>1 case in Y direction ",false,1,"long"); cmd.add(Ny);
		ValueArg<long> Nz("","Nz","size of grid in n>1 case in Z direction ",false,1,"long"); cmd.add(Nz);
		
		//		MultiArg<long> L("L","boxLength","length of the box in each direction. Defines Lx when single number is given. Lx, Ly if two values are passed and Lx Ly Lz when three values are passed",false,"long"); cmd.add(L);
		ValueArg<double> Lx("","Lx","physical size of grid in n>1 case in X direction (default: 1)",false,1,"double"); cmd.add(Lx);
		ValueArg<double> Ly("","Ly","physical size of grid in n>1 case in Y direction (default: 1)",false,1,"double"); cmd.add(Ly);
		ValueArg<double> Lz("","Lz","physical size of grid in n>1 case in Z direction (default: 1)",false,1,"double"); cmd.add(Lz);
		
		ValueArg<double> regridFact("r","regridFact","regrid factor: defines how many time denser the sampling of the initial function must be "
				"(using linear interpolation) when casting the input data onto a equally spaced grid. Only used for 1d ffts. The number is given "
				"wrt the mean spacing in the provided data set -- corresponding to regridFact=1. (default 0)",false,0,"double"); cmd.add(regridFact);
		
		SwitchArg saveReIm("","saveReIm", "Save the real and imaginary part of the real signal fft.", false);	cmd.add(saveReIm);

		ValueArg<long> Nrows("","Nrows","number of rows to read from the input file for further analysis (default: all rows are read)",false,-1,"long"); cmd.add(Nrows);
		ValueArg<long> skipFirst("","skipFirst","skip first N forws when reading the data from the input file (default: 0)",false,0,"long"); cmd.add(skipFirst);
		SwitchArg printData("","printData", "print the data that was loaded before further analysis", false);	cmd.add(printData);
		SwitchArg commentedFile("","commentedFile", "indicate that the input files may have commented rows - i.e starting from #", false);	cmd.add(commentedFile);

		SwitchArg saveData("","saveData", "save the input data to a file: _infile-col_x-col_y.txt (default: false)", false);	cmd.add(saveData);
		
		SwitchArg save3DboxSlices("","save3DboxSlices", "Save all slices through 3D box. These will be kept in a subdirectory", false);	cmd.add(save3DboxSlices);
		ValueArg<long> save3DboxSlicesPlane("","3DboxSlicesPlane","Defines in which planes the 3d box will be sliced. 0 - for yz, 1 - for xz and 2 for xy (default: 0)",false,0,"long"); cmd.add(save3DboxSlicesPlane);
		
		SwitchArg mkNoise("","mkNoise", "generates gaussian noise in real space according to some power law function defined by parameters A,ns. (Can also be used for 3D simulations of noise via fft)", false);	cmd.add(mkNoise);
		SwitchArg mkGNoise("","mkGNoise", "generates white gaussian random noise in real space with unit variance and zero mean. It can be smoothed latter on with smooth3DGauss option. 3D use only", false);	cmd.add(mkGNoise);
		SwitchArg mkPhaseNoise("","mkPhaseNoise", "generates white gaussian phase noise. The power spectrum is exactly flat", false);	cmd.add(mkPhaseNoise);
		SwitchArg mkSignal("","mkSignal", "generates signal instead from reading from file. The signal type is defined by the signalType option.", false);	cmd.add(mkSignal);
		ValueArg<string> signalType("","signalType","defines what function is generated instead of being read from file.\n"
				"These are the possible values for the mkSignal option:\n"
				"	topHat - generates the top hat function consisting of Nx points in total, with Ny points sitting on the top of the hat. "
				"The values on the top of the hat are 1 and the values outside of it are 0."
				"The function is defined symmetrically about x=0 with the width of the hat defined by mkSignal_x2 parameter and "
				"the width of the top of the hat defined by the mkSignal_x1 parameter\n"
				"	sineSineModulated - a 1 Hz modulated 50Hz sine function.\n"
				"	sineSqModulated - a 1 Hz square wave modulated 50Hz sine function.\n"
				"	twoSines - generates 50 Hz sine + 100 Hz sine offset in phase by 0.66 of its period. The --Lx parameter controls the lendth of the signal and --sampling controls the sampling."
				"	sinesFromDphiFile - a signal defined by a sum of sines functions with frequencies, amplitudes and phases"
				"defined by triplets of numbers saved in rows of a file. The file name should be specified after a colon: eg.\n"
				"		sinesFromDphiFile:fileName.txt."
				"The length of the signal generated is defined by --Lx parameter,"
				"and resolution is defined by _samplingRate parameter\n"
				"	signalFromReImFile - a signal defined by a sum of sines functions with frequencies, re and im"
				"defined by triplets of numbers saved in rows of a file. The file name should be specified after a colon: eg."
				"signalFromReImFile:fileName.txt. "
				"The length of the signal generated is defined by --Lx parameter,"
				"and resolution is defined by _samplingRate parameter\n"
				"   gauss3d - 3d gauss is generated as signal and placed in the center of the box. amplitude is 1 and sigma is defined by _s[xyz]G parameters. Only to be used with n=3.\n"
				"", false,"","string");	cmd.add( signalType);
		ValueArg<double> mkSignalX1("","mkSignal_x1","fist parameter defining the input signal to be generated (used with mkSignal option)",false,0,"double"); cmd.add(mkSignalX1);
		ValueArg<double> mkSignalX2("","mkSignal_x2","second parameter defining the input signal to be generated (used with mkSignal option)",false,0,"double"); cmd.add(mkSignalX2);
		ValueArg<double> A("A","A","amplitude of the noise (default: 1)",false,1,"double"); cmd.add(A);
		ValueArg<double> k0("","k0","pivot frequency (default: 1)",false,1,"double"); cmd.add(k0);
		ValueArg<double> ns("","ns","spectral index of the noise: -1 for 1/f noise (default: -1)",false,-1,"double"); cmd.add(ns);
		ValueArg<double> kmin("","kmin","minimal k value for the noise generation: (default: 1/N). This parameter is also used to define the range of frequencies for 1d FT calculation without fftw3.",false,0.001,"double"); cmd.add(kmin);
		ValueArg<double> kmax("","kmax","maximal k value for the noise generation: (default: 1). This parameter is also used to define the range of frequencies for 1d FT calculation without fftw3.",false,1,"double"); cmd.add(kmax);
		ValueArg<double> dk("","dk","step at which to probe the kmin,kmax range in 1d FT that doesn't use the fftw3: (default: 1)",false,1,"double"); cmd.add(dk);
		SwitchArg nofftw("","nofftw", "if used the power spectrum will be calculated using the slow O(N^2) method, but for arbitrary range of k values as defined by kmin,kmax,dk parameters.", false);	cmd.add(nofftw);
		ValueArg<double> seed("","seed","seed for RNs generation: (default: from date)",false,-1,"double"); cmd.add(seed);
		
		ValueArg<double> convolveGauss("","convolveGauss","convolve the input function with gaussian kernel of given FWHM (default -1 - not done). "
				"",false,-1,"double"); cmd.add(convolveGauss);
		SwitchArg averageFns("","averageFns", "Average functions given as unlabeled parameters. All files must be of the same length and format", false);	cmd.add(averageFns);
		ValueArg<long> colx("","colx","number of column for function arguments (starting from 0). If -1 given then it is ignored and arguments will be enumerated from 0 as rows in the input file.",false,0,"long"); cmd.add(colx);
		ValueArg<long> coly("","coly","number of column for function values (starting from 0 - must be greater than colx)",false,1,"long"); cmd.add(coly);
		ValueArg<long> colx2("","colx2","number of column for function arguments (starting from 0) for second input file. If -1 given then it is ignored and arguments will be enumerated from 0 as rows in the input file.",false,0,"long"); cmd.add(colx2);
		ValueArg<long> coly2("","coly2","number of column for function values (starting from 0 - must be greater than colx) for second input file",false,1,"long"); cmd.add(coly2);
		ValueArg<long> ref1("","ref1","reference column number for the first input file. Only used with calcDistance2D option.",false,2,"long"); cmd.add(ref1);
		ValueArg<long> ref2("","ref2","reference column number for the second input file. Only used with calcDistance2D option.",false,2,"long"); cmd.add(ref2);
		ValueArg<string> correlation("","correlation","calculate correlation coefficients between the inf and inf2. "
				"Types of possible correlation calculations are:\n"
				"cutBeg - cut out from inf and inf2 a common number of samples from the beginning stopping at whichever data set is shorter. \n"
				"cutEnd - cut out from inf and inf2 a common number of samples from the end stopping at whichever data set is shorter.\n"
				"match - use x and x2 information to interpolate the two functions to a common argument values. This uses the -r option for interpolations.\n"
				"If --Nx option is given then only Nx samples are used.\n"
				"(default: none)",false,"","string"); cmd.add(correlation);
		ValueArg<long> corrStep("","corrStep","offset in number of samples for correlation function calculations.",false,0,"long"); cmd.add(corrStep);
		ValueArg<long> corrNum("","corrNum","Number of points in the correlation coefficient function - corresponds to the "
				"number of times of shifting one function wrt other."
				"(default: 0 - means until period wrap.)",false,0,"long"); cmd.add(corrNum);
		
		SwitchArg calcDistance1D("","calcDistance1D", "Calculate distance between point sets defined in the two files."
				"If the point sets are defined for different x values interpolation is done in the second set to the x values from the"
				"first set.", false);	cmd.add(calcDistance1D);
		SwitchArg calcDistance2D("","calcDistance2D", "Calculate distance between point sets defined in the two files."
				"The x and y data are parameterized. Each file reads in the parameter value column (as x values) and x coordinates column (as y values). The y coordinate column is"
				"indicated by the option ref and ref2 respectively for the first and the second set."
				"If the point sets are defined for different x values interpolation is done in the second set to the x values from the"
				"first set. This option allows to specify one additional column in each of the the two input files that holds"
				"the y-values.", false);	cmd.add(calcDistance2D);
		
		SwitchArg shift0("","shift0", "shift the input signal to zero mean before fftw. This only applies now to 1d transforms", false);	cmd.add(shift0);
		SwitchArg hann("","hann", "Apply the hann window before fftw.", false);	cmd.add(hann);
		ValueArg<double> sampling("","sampling","sampling rate of the input signal in case when colx=-1 was given [Hz] (default: 1)",false,1,"double"); cmd.add(sampling);
		ValueArg<double> sampling2("","sampling2","sampling rate of the input signal in case when colx2=-1 was given [Hz] (default: 1)",false,1,"double"); cmd.add(sampling2);
		
		ValueArg<long> sliceAndAverage("","sliceAv","slice the input signal into sliceAv equal slices calculate the spectra of each slice and produce averaged spectra from all slices. "
				"For values approaching infinity the result approaches to power spectral density (PSD) of the input signal. "
				"A typically used value here is sqrt(N_tot) which is a trade-off between noise and spectral resolution, where N_tot is the total number of samples in the input signal"
				"The input signal is chopped off at the end so as to fit its length to the integer multiplicity of the slice length (default: 1 i.e. no slicing is done)",false,1,"long"); cmd.add(sliceAndAverage);
		SwitchArg psd("","psd", "same as setting sliceAndAverage=sqrt(N_tot). In this case saveSlices option is ignored", false);	cmd.add(psd);
		
		SwitchArg filter("","filter", "Indicates that the input signal should be filtered in fourier space by removing ranges of freqiencies defined in "
				"text file specified by the filterData.", false);	cmd.add(filter);
		ValueArg<string> filterData("","filterData","the name of a text file with two columns defining ranges of frequencies to be removed in fourier space."
				"The first column defines the central frequency to be removed and the second column defines the half-width of the range to be removed."
				"This is so unless the fromTo modifier is used. The string can also be a comma separated list of central frequencies (with widths defined by --w option). (default: )",false,"","string"); cmd.add(filterData);
		
		SwitchArg printStats("","printStats", "print some statistics on the input signal before calculations.", false);	cmd.add(printStats);
		SwitchArg saveSlices("","saveSlices", "Save slices in separate files.", false);	cmd.add(saveSlices);
		SwitchArg save3Dpower("","save3Dpower", "Save 3D power spectra calculated from generated 3-D field.", false);	cmd.add(save3Dpower);
		ValueArg<long> calculate1Dto3DpowerCalibFact("","calculate1Dto3DpowerCalibFact","3D OPTION. If the provided power spectrum parameters (A,k0,ns) are for"
				"1D power spectrum, then 1D to 3D power spectrum conversions may be needed. This is done numerically. "
				"This option specifies the number of retrieved 1D power spectra from the 3D box sampled along X axis and calculated using 1D FFT."
				"If this option is > 0 then it is assumed that these corrections should be applied and so ns parameter is dectreased by 2 for 3D field generation."
				"Then slices are drawn and power spectrum calculated. The power calibration factor is calculated at k0 frequency as an average power from the "
				"reconstructed 1D spectra. It is known that due to projection effects the 1D reconstructed spectra will have more power on large scales than "
				"the input 1D power, so this procedure may help in mitigating this effect. Finally, the simulated field is divided in real space"
				"by sqrt(<P_1drec>/P_1din). If this option is 0, then nothing of this sort is done. (default: 0)",false,0,"long"); cmd.add(calculate1Dto3DpowerCalibFact);
		
		
		
		SwitchArg calibrateByDuration("","calibrateByDuration", "The output spectra will be calibrated by the duration of the recording (or slice): i.e the power spectrum will be divided by ( (points count in real space)/(sampling rate) )", false);	cmd.add(calibrateByDuration);
		ValueArg<long> calibrateByStDev("","calibrateByStDev","Calibrate the signal by the piece-wise calculated standard deviation in the real space. "
				"The width of the window is defined by this parameter (It must be larger than 3 for this option to be activated). By default no calibration is done. The calibration is done not on the moving window."
				"The signal is calibrated in the neighbouring, non-overlapping windows.",false,-1,"long"); cmd.add(calibrateByStDev);
		ValueArg<long> movingWindowStep("","windowShift","number of samples by which to move the window."
				"This is used in combination with option calibrateByStDev"
				"The default value is the size of the window.",false,-1,"long"); cmd.add(movingWindowStep);
		
		
		ValueArg<string> prewhiten("","prewhiten","load the power spectrum from file indicated by this option and use it to calibrate the re and im and P components. "
				"This simply gives re/power, im/power, P/power where power is the power spectrum loaded from file and interpolated onto the needed frequency: (default: )",false,"","string"); cmd.add(prewhiten);		
		SwitchArg equalizePower("","equalizePower", "Indicates that the input signal should be filtered in fourier space by equalizing the power in all frequencies.", false);	cmd.add(equalizePower);
		
		SwitchArg smooth3DGauss("","smooth3DGauss", "Smooth the input field with gaussian kernel defined by smoothing length parameters: sx,sy,sz.", false);	cmd.add(smooth3DGauss);
		ValueArg<double> sxG("","sxG","x smoothing length for gaussian smoothing given in units of Lx (default: 1)",false,1,"double"); cmd.add(sxG);
		ValueArg<double> syG("","syG","y smoothing length for gaussian smoothing given in units of Ly (default: 1)",false,1,"double"); cmd.add(syG);
		ValueArg<double> szG("","szG","z smoothing length for gaussian smoothing given in units of Lz (default: 1)",false,1,"double"); cmd.add(szG);

		vector<string> outputFormats;
		outputFormats.push_back("reim2");
		outputFormats.push_back("reim2nonzero");
		outputFormats.push_back("reim1nonzero");
		outputFormats.push_back("Dphi2nonzero");
		outputFormats.push_back("Dphi1nonzero");
		outputFormats.push_back("df3");
		outputFormats.push_back("hdf5");
		ValuesConstraint<string>* outputFormatsNew = new ValuesConstraint<string>(outputFormats);
		ValueArg<string> outputFormat("","outputFormat","output format of the re and im data. Possible values are:\n\n"
				"reim2 - (default) - output in two files one for re part and one for im part, each storing all frequencies\n"
				"reim2nonzero - output in two files one for re part and one for im part, each storing only non-zero frequencies\n"
				"reim1nonzero - output in one file with 3 columns: frequency, re, im. Storing only non-zero frequencies\n"
				"Dphi2nonzero - output in two files, one for modulus and one for phase. Storing only non-zero frequencies.\n"
				"Dphi1nonzero - output in one file, with three columns: frequency, modulus, phase. Storing only non-zero frequencies.\n"
				"df3 - df3 file format for 3D density data for POV-RAY ray-tracing.\n"
				"raw - raw file format for storing 3D density data (default for -n 3 option).\n"
				"hdf5 - hdf5 file format for storing 3D density data. The data are saved in dataset 'field'\n"
				"",false,"reim2",outputFormatsNew); cmd.add(outputFormat);

		//		MultiArg<double> extract_f("e","extractf","list of arbitrary bins to use for binning (if defined, the f,t,d args will be ignored)",false,"long"); cmd.add(b);
		ValueArg<string> extract_f("e","extract_f","coma-separated list of frequencies to extract from the input signal "
				"along with frequencies laying within range [e-w, e+w), where w is the half width of the range defined by w parameter."
				"If no comas are provided and no w or extract_w option is provided then this string is interpreted as file name with ranges in a form center of the bin in the first"
				"column and half width of the range in the second column. This convention can be changed with option fromTo in which case the first column is"
				"considered to hold lower bound of the range and second column is considered to hold the upper bound of the range. (default: none)",false,"","string"); cmd.add(extract_f);
		ValueArg<string> extract_w("w","extract_w","coma-separated list of ranges half-widths corresponding to the central frequencies defined by extract_f parameter"
				"It need not be of the same length as the extract_f parameter as the program will cycle through the array as i % N where N is the size of the array. (default: none)",false,"","string"); cmd.add(extract_w);

		SwitchArg fromTo("","fromTo", "Switch to alter the convention of the columns in the input file for extracting data ranges from the spectra in the fourier domain.", false);	cmd.add(fromTo);

//		ValueArg<string> removeRanges("","removeRanges","file name with frequency range in each row to be removed in the Fourier space from the input signal. "
//				"The filtered signal can be obtained by using invfft option."
//				"(default: none)",false,"","string"); cmd.add(removeRanges);
//		SwitchArg invfft("","invfft", "Calculate inverse fft on the filtered signal in the fourier space. The output is stored in the output file with suffix .invfft", false);	cmd.add(invfft);
		
		
		SwitchArg testSin("","testSin", "test of sine function generation; useful for other tests", false);	cmd.add(testSin);
		
		SwitchArg testParseval("","testParseval", "Run Parseval theorem test of the FFT for a sine input signal. "
				"sum_i=0^N-1 |x_i|^2 = 1/N sum_k=0^N-1 |X_k|^2 ,"
				"where X_k is DFT if x_i both of length N (i.e. the summation is over positive and negative frequencies.", false);	cmd.add(testParseval);
		
		SwitchArg testPtSrcDistr("","testPtSrcDistr", "Run a simulation of the point source distribution in 1D"
				"and calculate its power spectrum", false);	cmd.add(testPtSrcDistr);
		
		SwitchArg testFFTcoefdistrchisq("","testFFTcoefdistrchisq", "Test what is the distribution of the FFT coefficients of the"
				"input signal wihch is a realization of white chisq noise.", false);	cmd.add(testFFTcoefdistrchisq);
		
		SwitchArg testFFTcoefdistrchisq2("","testFFTcoefdistrchisq2", "Test the averaging of the power spectra "
				"to speed up the calculation of the covariance matrix in the MCMC OCRA noise parameter estimation", false);	cmd.add(testFFTcoefdistrchisq2);
		SwitchArg testCorrelation("","testCorrelation", "Test the calculation of correlation coefficient function", false);	cmd.add(testCorrelation);
		SwitchArg test3Dsmooth("","test3Dsmooth", "Test the 3D smoothing using the impulse signal in the center of the field. The smoothing length is defined by sx,sy,sz parameters.", false);	cmd.add(test3Dsmooth);
		SwitchArg test3Dsmooth2("","test3Dsmooth2", "Test the 3D smoothing for the case of convolution of the gaussian input signal itself. Test for the fwhm increase.", false);	cmd.add(test3Dsmooth2);
		
		SwitchArg test("","test", "starts a demo test of the program.", false);	cmd.add(test);
		SwitchArg test1of("","test1of", "starts a demo test of the program of noise properties.", false);	cmd.add(test1of);
		SwitchArg testKasdin("","testKasdin", "starts a test of the 1/f rng in real space.", false);	cmd.add(testKasdin);
		SwitchArg testGaussoofNormalized("","testGaussoofNormalized", "starts a test of generation gaussian noise then 1/f noise with the same number of samples: 1e6, "
				"saves both and saves their sum. The parameters required for this test is Nx -- the number of samples to generate", false);	cmd.add(testGaussoofNormalized);
		SwitchArg testPk3d("","testPk3d", "starts test of the generation of the 3d random gaussian filed with initial power spectrum defined by A,k0, and ns.", false);	cmd.add(testPk3d);
//		SwitchArg testModulation("","testModulation", "This test provides a modulated sine signal and calculates its power spectrum according to normal parameters set.", false);	cmd.add(testModulation);

		//
		// Parse the command line.
		//
		cmd.parse(argc,argv);
		delete allowedDimsNew;
		delete outputFormatsNew;
		
		//
		// Set variables
		//
		_inputFns=inputFns.getValue();
		_infile = input_file.getValue(); /*if ( _input_file == "nofile" || _input_file == "Ylm" ) { _not_from_file = true; } else { _not_from_file = false; }*/
		_infile2 = input_file2.getValue(); 
		_outfile = output.getValue(); if (!output.isSet()) { sprintf(tmpch,"%s-regrid_%.2lf.power",_infile.c_str(),_regridFactor); _outfile=tmpch; }
		_outputFormat=outputFormat.getValue();
		_dimension = dimension.getValue(); if (_dimension==3 and outputFormat.isSet()==false) { _outputFormat="raw"; }
		_format=format.getValue();
		_averageFns=averageFns.getValue();
		_printStats=printStats.getValue();
		_testName = test.getValue();
		_testKasdin = testKasdin.getValue();
		_test1of = test1of.getValue();
		_testPk3d = testPk3d.getValue();
		_testGaussoofNormalized= testGaussoofNormalized.getValue();
		_testCorrelation = testCorrelation.getValue();
//		_testModulation=testModulation.getValue();
		_calcDistance1D=calcDistance1D.getValue();
		_calcDistance2D=calcDistance2D.getValue();
		
		_dimension=dimension.getValue();
		_shift0=shift0.getValue();
		_hann=hann.getValue();
		
		//		std::vector<long> Nv=N.getValue();
		//		switch (Nv.size()) {
		//			case 0:
		//				_Nx=_Ny=_Nz=0;
		//				break;
		//			case 1:
		//				_Nx=Nv[0];
		//				break;
		//			case 2:
		//				_Nx=Nv[0];
		//				_Ny=Nv[1];
		//				break;
		//			case 3:
		//				_Nx=Nv[0];
		//				_Ny=Nv[1];
		//				_Nx=Nv[2];
		//				break;
		//		}
		//		_Nx=N.getValue();
		_Nx=Nx.getValue();
		_Ny=Ny.getValue();
		_Nz=Nz.getValue();
		
		
		//		std::vector<long> Lv=L.getValue();
		//		switch (Lv.size()) {
		//			case 1:
		//				_Lx=Nv[0];
		//				break;
		//			case 2:
		//				_Lx=Nv[0];
		//				_Ly=Nv[1];
		//				break;
		//			case 3:
		//				_Lx=Lv[0];
		//				_Ly=Lv[1];
		//				_Lx=Lv[2];
		//				break;
		//		}
		_Lx=Lx.getValue();
		_Ly=Ly.getValue();
		_Lz=Lz.getValue();
		
		_regridFactor=regridFact.getValue();
		
		_save3DboxSlices=save3DboxSlices.getValue();
		_mkNoise=mkNoise.getValue();
		_mkSignal=mkSignal.getValue();
		_mkSignalType=signalType.getValue();
		_mkSignal_x1=mkSignalX1.getValue();
		_mkSignal_x2=mkSignalX2.getValue();
		
		_mkGnoise=mkGNoise.getValue();
		_mkPhaseNoise=mkPhaseNoise.getValue();
		_A=A.getValue();
		_k0=k0.getValue();
		_ns=ns.getValue();
		_kmin=kmin.getValue(); if (!kmin.isSet()) { if (_Nx!=0) _kmin=1.0/double(_Nx); else _kmin=0; }
		_kmax=kmax.getValue(); if (!kmax.isSet()) { _kmax=1.0; }
		_dk=dk.getValue(); if (!dk.isSet()) { _dk=1.0; }
		_nofftw=nofftw.getValue(); 		
		_seed=seed.getValue(); if (!seed.isSet()) { _seed=0; }
		
		
		_convolve=convolveGauss.getValue();
		_colx=colx.getValue();
		_coly=coly.getValue();
		_colx2=colx2.getValue();
		_coly2=coly2.getValue();
		_ref1=ref1.getValue();
		_ref2=ref2.getValue();
		_correlation=correlation.getValue();
		_corrStep=corrStep.getValue();
		_corrNum=corrNum.getValue();
		_sliceAndAverage=sliceAndAverage.getValue();
		_samplingRate=sampling.getValue();
		_samplingRate2=sampling2.getValue();
		_saveReIm=saveReIm.getValue();
		_Nrows=Nrows.getValue();
		_skipFirst=skipFirst.getValue();
		_printData=printData.getValue();
		_commentedFile=commentedFile.getValue();
		_saveData=saveData.getValue();
		
		_saveSlices=saveSlices.getValue();
		_save3Dpower=save3Dpower.getValue();
		_calculate1Dto3DpowerCalibFact=calculate1Dto3DpowerCalibFact.getValue();
		_save3DboxSlicesPlane=save3DboxSlicesPlane.getValue();
		_calibrateByDuration=calibrateByDuration.getValue();
		_calibrateByStdev=calibrateByStDev.getValue();
		_movingWindowStep=movingWindowStep.getValue(); if (_movingWindowStep==-1) _movingWindowStep=_calibrateByStdev;
		_prewhiten=prewhiten.getValue();
		_equalizePower=equalizePower.getValue();
		_psd=psd.getValue();
		_filter=filter.getValue();
		_filterData=filterData.getValue();
		_test_sin=testSin.getValue();
		_testParseval=testParseval.getValue();
		_testPtSrcDistr=testPtSrcDistr.getValue();
		_testFFTcoefdistrchisq=testFFTcoefdistrchisq.getValue();
		_testFFTcoefdistrchisq2=testFFTcoefdistrchisq2.getValue();
		_test3DGaussSmooth=test3Dsmooth.getValue();
		_test3DGaussSmooth2=test3Dsmooth2.getValue();
		_smooth3DGauss=smooth3DGauss.getValue();
		_sxG=sxG.getValue();
		_syG=syG.getValue();
		_szG=szG.getValue();
		
		QString extract;
		QStringList qsl;
		bool ok;

		_fromTo=fromTo.getValue();
		extract=extract_w.getValue().c_str();
		if (extract.size()>0) {
			qsl=extract.split(',');
			//		printf("qsl size: %li\n",qsl.size());
			//		cout << qsl[0].toStdString();
			//		exit(0);
			for (long i = 0; i < qsl.size(); i++) {
				_extract_w.append(qsl[i].toDouble(&ok));
				assert(ok==true);
			}
		}
		extract=extract_f.getValue().c_str();
		if (extract.size()>0) {
			if (extract_w.getValue().length()==0) {
				mscsFunction r;
				r.load(extract_f.getValue());
				for (long i = 0; i < r.pointsCount(); i++) {
					if (_fromTo) { _extract_f.append(0.5*(r.getX(i)+r.f(i))); _extract_w.append(0.5*fabs(r.getX(i)-r.f(i))); }
					else {
						_extract_f.append(r.getX(i)); _extract_w.append(r.f(i)); 
					}
				}
			}
			else {
				qsl=extract.split(',');
				for (long i = 0; i < qsl.size(); i++) {
					cout << "Will extract frequency range: " << qsl[i].toStdString() << " +_ "<< _extract_w[i % _extract_w.size()] << "\n";
					_extract_f.append(qsl[i].toDouble(&ok));
					assert(ok==true);
				}
				if (_extract_f.size()>0) assert(_extract_w.size()>0);				
			}
		}		
//		_invfft=invfft.getValue();
//		_removeRanges.load(removeRanges.getValue());
		
		
	} catch ( ArgException& e )
	{ cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
	
}



mscsFunction calculateCorrelationFunction(mscsFunction &f1, mscsFunction &f2,cpedsMsgs* msgs) {
	
	if (_testCorrelation) {
		msgs->say("testing correlations",Top);
		f1.clearFunction();
		f2.clearFunction();
		f1.mkSin(0,1,0.01,1,0);
		f2.mkSin(0.05,1,0.012,1,0.5);
		//		f1.mkSin(0,1,0.1,1,0);
		//		f2.mkSin(0,1,0.1,1,0.25);
		f1.save("corrTest.f1");
		f2.save("corrTest.f2");
	}
	
	
	mscsFunction f;
	long n;
	long N1,N2;
	N1=f1.pointsCount();
	N2=f2.pointsCount();
	long iSt,iEn;
	double X1st = f1.X(0);
	double X1en = f1.X(N1-1);
	double X2st = f2.X(0);
	double X2en = f2.X(N2-1);
	double xSt,xEn,dx;
	
	if (_Nx!=1) { n=_Nx; }
	else {
		n=cpeds_get_min(N1,N2);
	}
	
	if (_correlation=="cutBeg") {
	}
	if (_correlation=="cutEnd") {
	}
	if (_correlation=="match") {
		//
		// sanity checks
		//
		if (X1st>X1en) { msgs->error("Arguments are not rising in f1",High); }
		if (X2st>X2en) { msgs->error("Arguments are not rising in f2",High); }
		if (X1en<X2st) { msgs->error("The two functions are not defined over a common range",High); }
		if (X2en<X1st) { msgs->error("The two functions are not defined over a common range",High); }
		
		
		//
		// find common range of arguments
		//
		xSt=cpeds_get_max(f1.X(0),f2.X(0));
		xEn=cpeds_get_min(f1.X(N1-1),f2.X(N2-1));
		//
		printf("xSt: %lf, xEn: %lf\n",xSt,xEn);
		
		//
		// interpolate 
		//
		double dx;
		if (N1>N2) dx=(xEn-xSt)/f2.pointsCount(); else dx=(xEn-xSt)/f1.pointsCount();
		dx/=_regridFactor;
		
		//		f1.interpolate(xSt,xEn,dx,true,"linear",true);
		//		f2.interpolate(xSt,xEn,dx,true,"linear",true);
		f1.interpolate(xSt,xEn,dx,true,"cspline",true);
		f2.interpolate(xSt,xEn,dx,true,"cspline",true);
		
		f1.deleteOutsideOf(xSt,xEn);
		f2.deleteOutsideOf(xSt,xEn);
		
		if (_testCorrelation) {
			f1.save("corrTest.f1.int");
			f2.save("corrTest.f2.int");
		}
		
		
	}
	//
	// calculate correlation function
	//
	
	f=f1.correlationCoefficientFunction(f2,_corrStep,_corrNum);
	
	return f;
}

mscsFunction  calculateDistance1D(mscsFunction &f1, mscsFunction &f2,cpedsMsgs* msgs) {
	mscsFunction d;
	double *x = f1.extractArguments();
	mscsFunction s2=f2.extrapolate(f1.X(0),f1.last().rx(),f1.X(1)-f1.X(0),true,"linear");
	//	s2.print();
	s2.interpolate(x,f1.pointsCount(),"linear",true);
	//	s2.print();
	long N=f1.pointsCount();
	for (long i = 0; i < N; i++) {
		d.newPoint(f1.X(i),fabs(s2.f(i)-f1.f(i)));
	}
	//	d.print();
	delete [] x;
	return d;
}


void test3DGaussSmooth(cpedsMsgs* msgs) {
	msgs->say("Testing 3D gaussian smoothing",High);
	mscsFunction3dregc f(_Nx,_Ny,_Nz, _Lx/_Nx,_Ly/_Ny,_Lz/_Nz);
	f.setName("f");
	
	// make real space signal
	//	fftw_complex z={0,0};
	//	f=z;
	//	f(_Nx/2, _Ny/2, _Nz/2)[0]=1.0;
	f.generateRandomGaussianField(0,1,true,false,_seed,0);
	
	
	// check
	cpeds_matrix_save(f.getSlice(0,_Nx/2,0),"slice");
	
	//	f.smooth3DGauss(_sxG,_syG,_szG);
	//
	// make convolution kernel in real space
	//
	mscsFunction3dregc ker(_Nx,_Ny,_Nz, _Lx/_Nx,_Ly/_Ny,_Lz/_Nz);
	ker.setName("kernel");
	ker.mkGauss3D(_sxG,_syG,_szG);
	ker.shift(_Nx/2,_Ny/2, _Nz/2);
	// check
	ker.getSlice1D(0,_Ny/2,_Nz/2,0);
	cpeds_save_matrix(ker.getSlice1D(0,_Ny/2,_Nz/2,0),_Nx,1,"ker1d",true);
	cpeds_matrix_save(ker.getSlice(0,0,0),"ker-slice");
	
	// make kernel
	ker.fft(true);
	ker.absoluteValue();
//	ker.re2im();
	ker/=ker.getMaxValue();
	cpeds_save_matrix(ker.getSlice1D(0,_Ny/2,_Nz/2,0),_Nx,1,"ker1d",true);
	cpeds_matrix_save(ker.getSlice(0,0,0),"kerFFT-slice");
	
	//
	// convolve
	//
	f.fft(true);
//	ker*=double(f.size()); 
	f*=ker;
	//	f.antisymmetrize();
	f.fft(false);
	
	// check
	cpeds_matrix_save(f.getSlice(0,_Nx/2,0),"slice-sm");
	
	
	exit(0);
}
void test3DGaussSmooth2(cpedsMsgs* msgs) {
	msgs->say("Testing 3D gaussian smoothing",High);
	_Nx=101;
	_Ny=_Nx;
	_Nz=_Nx;
	_Lx=1;
	_Ly=1;
	_Lz=1;
	mscsFunction3dregc f(_Nx,_Nx,_Nx, _Lx/_Nx,_Lx/_Ny,_Lx/_Nz);
	f.setName("f");
	
	// make real space signal
	//	fftw_complex z={0,0};
	//	f=z;
	//	f(_Nx/2, _Ny/2, _Nz/2)[0]=1.0;
	double sigma=0.1;
	double fwhm=2.0*sqrt(2*log(2.0))*sigma;
	f.mkGauss3D(sigma,sigma,sigma);
	
	printf("will be convolving with the same gauss function\n");
	printf("input fwhm: %lE\n",fwhm);
	printf("output fwhm should be: %lE\n",sqrt(2)*fwhm);
	// check
	cpeds_matrix_save(f.getSlice(0,_Nx/2,0),"slice");
	
	//	f.smooth3DGauss(_sxG,_syG,_szG);
	//
	// make convolution kernel in real space
	//
	mscsFunction3dregc ker(_Nx,_Ny,_Nz, _Lx/_Nx,_Ly/_Ny,_Lz/_Nz);
	ker.setName("kernel");
	ker.mkGauss3D(sigma,sigma,sigma);
	ker.shift(_Nx/2,_Ny/2, _Nz/2);
	// check
	ker.getSlice1D(0,_Ny/2,_Nz/2,0);
	cpeds_save_matrix(ker.getSlice1D(0,_Ny/2,_Nz/2,0),_Nx,1,"ker1d",true);
	cpeds_matrix_save(ker.getSlice(0,0,0),"ker-slice");
	
	// make kernel
	ker.fft(true);
	ker.absoluteValue();
	ker/=ker.getMaxValue();
	cpeds_save_matrix(ker.getSlice1D(0,_Ny/2,_Nz/2,0),_Nx,1,"ker1d",true);
	cpeds_matrix_save(ker.getSlice(0,0,0),"kerFFT-slice");
	
	//
	// convolve
	//
	f.fft(true);
//	ker*=double(f.size()); 
	f*=ker;
	//	f.antisymmetrize();
	f.fft(false);
	
	// check
	cpeds_matrix_save(f.getSlice(0,_Nx/2,0),"slice-sm");
	
	
	exit(0);
}

void smooth3DGauss(mscsFunction3dregc& f, cpedsMsgs* msgs) {
	msgs->say("Smoothing 3D box with gaussian kernel",High);
	f.smooth3DGauss(_sxG,_syG,_szG);
}

/***************************************************************************************/
void loadInputData(mscsFunction& f,mscsFunction& f2, cpedsMsgs* msgs) {

	//
	// load f
	//
	
	if (_format=="ascii") {
			f.load(_infile,_commentedFile,_colx,_coly,_Nrows,_skipFirst);
	} 
	else {
		if (_format=="binDAQd") {
			loadDAQdFormatData(f,msgs);
		}
		else {
			msgs->criticalError("this format is not dedicated for 1d data. Will quit now.",Top);
		}
	}
	if (_colx==-1 and _samplingRate!=1) {
		f.scaleX(1.0/_samplingRate);
	}
	else {
		if (_samplingRate!=1) printf("WARNING: samplingRate was given but will be ignored since the column X timings are used. Use --colx=-1 to override it\n");
	}
	if (_printData) f.print();

	//
	// load f2
	//
	
	if (_infile2!="") {
		f2.load(_infile2,_commentedFile,_colx2,_coly2,_Nrows,_skipFirst);
		if (_colx2==-1 and _samplingRate2!=1) {
			f2.scaleX(1.0/_samplingRate2);
		}
		else {
			if (_samplingRate2!=1) printf("WARNING: samplingRate2 was given but will be ignored since the column X timings are used. Use --colx2=-1 to override it\n");
		}
		if (_printData) f2.print();
	}

	
}

/***************************************************************************************/
void calculateFunctionDistance1D(mscsFunction& f,mscsFunction& f2, cpedsMsgs* msgs) {
	mscsFunction distance = calculateDistance1D(f,f2,msgs);
	distance.setName("distance1D");
	msgs->say("Average between two point sets is: "+msgs->toStr(distance.meanf()),High);
	msgs->say("distance st.dev: "+msgs->toStr(distance.stdev()),High);
	
	distance.save(_outfile);
	delete msgs;
	exit(0);
}
/***************************************************************************************/
void calculateFunctionDistance2D(mscsFunction& f,mscsFunction& f2, cpedsMsgs* msgs) {
	cpedsFunctionDistance fd;
	mscsFunction s1y("s1y"), s2y("s2y");
	s1y.load(_infile,false,_colx,_ref1,_Nrows,_skipFirst);
	s2y.load(_infile2,false,_colx2,_ref2,_Nrows,_skipFirst);
	
	mscsFunction distance = fd.distance(f,s1y,f2,s2y);
	msgs->say("Average between two point sets is: "+msgs->toStr(distance.meanf()),High);
	msgs->say("distance st.dev: "+msgs->toStr(distance.stdev()),High);
	
	distance.save(_outfile);
	delete msgs;
	exit(0);
}

/***************************************************************************************/
void loadDAQdFormatData(mscsFunction& f, cpedsMsgs* msgs) {
	long long N=0;
	double Nd;
	if (_coly < 0 or _coly>11) { msgs->criticalError("for this format the coly can assume values only from range 0..10.",Top); }

	FILE* F=fopen(_infile.c_str(),"rb"); if (F==NULL) { msgs->criticalError("cannot open file: "+_infile, Top); }
	fclose(F);
	
	
	struct stat info;
	stat(_infile.c_str(),&info);
	Nd = ((double)(info.st_size))/(long)sizeof(ocrafDAQd_t);
	if (Nd - (long long)Nd != 0) msgs->warning("The size of the input file is incorrect.",High);
	N=(long long)(Nd);
	
	msgs->say("There are "+msgs->toStr(N)+" valid rows in the input file.",High);

	if (N==0) msgs->criticalError("there is N=0 rows in the input file. Will exit now.",Top);
	if (_Nrows!=-1) {
		if (_skipFirst+_Nrows>N) {
			msgs->warning("The _skipFirst+_Nrows>N. WIll skip skipFirst and read until the end of the file",Top);
			_Nrows=N-_skipFirst;
		}
	}
	if (_Nrows==-1)	_Nrows=N;
	F=fopen(_infile.c_str(),"rb"); 

	ocrafDAQd_t row;
	size_t s=sizeof(ocrafDAQd_t);
	long n;
	if (_skipFirst>0 and _skipFirst<N) n=fseek(F,_skipFirst*sizeof(ocrafDAQd_t),SEEK_SET);
	assert(_coly>0);
	if (_colx==-1) {
		for (long long i=0;i<_Nrows;i++) { n=fread(&row,s,1,F);  f.newPoint(double(i),row.ch[_coly-1]); 	}
	}
	else {
		if (_colx==0) {
			for (long long i=0;i<_Nrows;i++) { n=fread(&row,s,1,F);  f.newPoint(double(row.i),row.ch[_coly-1]); 	}		
		}
		else {
			for (long long i=0;i<_Nrows;i++) { n=fread(&row,s,1,F);  f.newPoint(row.ch[_colx-1],row.ch[_coly-1]); 	}
		}
	}
	fclose(F);
	if (_saveData) {
		if (_Nrows!=-1)
			f.save(_infile+"-Nrows_"+msgs->toStr(_Nrows)+"-skipFirst_"+msgs->toStr(_skipFirst)+"-col_"+msgs->toStr(_colx)+"-"+msgs->toStr(_coly)+".txt");
		else
			f.save(_infile+"-col_"+msgs->toStr(_colx)+"-"+msgs->toStr(_coly)+".txt");
	}
}
/***************************************************************************************/
void save_run_config() {
	
}
/***************************************************************************************/
void averageFunctions(cpedsMsgs* msgs) {
	mscsFunction f("f"),av("average");
	for (long i = 0; i < _inputFns.size(); i++) {
		msgs->say("loading function: "+_inputFns[i],High);
		f.load(_inputFns[i],false,_colx,_coly);
		if (i==0) av=f; else av+=f;
	}
	msgs->say("averaging: ",High);
	av/=double(_inputFns.size());
	av.save(_outfile);
}
/***************************************************************************************/
void mkTestParseval(int testNo) {
	mscsFunction f;
	
	printf("--------------------- NEW TEST ---------------------\n");
	if (testNo==0) {
		printf("--------------------- signal: white nosise m:0 s:1---------------------\n");
		f.mkGaussianNoise(1000,0,1);
		f-=f.meanf();
	}
	if (testNo==1) {
		printf("--------------------- signal: white nosise m:1 s:1---------------------\n");
		f.mkGaussianNoise(1000,1,1);
	}
	if (testNo==2) {
		f.mkGaussianNoise(1001,1,1);
	}
	if (testNo==3) {
		printf("--------------------- signal: sine, period: 1 --------------------\n");
		f.mkSin(0,10,0.01,1,0);
	}
	if (testNo==5) {
		printf("--------------------- signal: white nosise m:0 s:1, arbitrary sampling---------------------\n");
		f.mkGaussianNoise(1000,1,1);
		f.scaleX(0.1);
	}
	printStatistics(f);
	mscsFunction re,im;
	mscsFunction p=f.powerSpectrum(&re,&im,0);
	//check of Parseval's theorem sum_i=0^N-1 |x_i|^2 = N dt sum_k=0^N-1 |X_k|^2 ,
	// where X_k is DFT of x_i and both are of length N (i.e. the summation is over positive and negative frequencies)
	// converting this so that LHS stands for st.dev. for <x>=0 signals we have
	double LHS=f.rms();  LHS*=LHS; //LHS*=double(f.pointsCount()); 
	double RHS=0;
	
	//		// for implementation when power spectrum takes sum of positive and negative frequencies into one k mode
	//		for (long i = 0; i < p.pointsCount(); i++) { // sum zero and positive frequencies power
	//			RHS+=p.f(i);
	//		}
	
	// for old (current) definition of the power spectrum where only positive frequencies participated
	for (long i = 0; i < p.pointsCount(); i++) { // sum zero and positive frequencies power
		RHS+=p.f(i);
	}
	for (long i = 1; i < p.pointsCount(); i++) { // sum negative frequencies power
		RHS+=p.f(i);
	}
	double dt=f.X(1)-f.X(0);
	printf("dt=%lE\n",dt);
	RHS*=dt*dt;
	/* BLcomment (Apr 24, 2012, 6:39:16 PM): this is commented out since the refefinition of the power spectrum where we divide the spectrum by the 1/N^2 factor
	 * to make the spectrum return the same power for the same signals that are sampled with different frequency. 
	 * The only effect should be the k_max in the output power spectra. Hence if the 1/N^2 factor goes in the the power spectra then we cannot divide by it here.
	 * This also fixes the issue of the mean value (DC mode). Now the sqrt(DC mode of the power spectrum) is the mean of the input signal
//	RHS/=f.pointsCount();  
//	RHS/=f.pointsCount();
	 	 	 	 	 	 	 */
	//		printf("%li\n",(2*p.pointsCount()-2));
	printf("power spectrum points count: %li\n",p.pointsCount());
	
	printf("LHS: %.15lE\n",LHS);
	printf("RHS: %.15lE\n",RHS);	
	printf("LHS-RHS: %.15lE\n",LHS-RHS);	
	printf("\n\nmeans test\n");
	printf("DC term: %.15lE\n",p.f(0));
	printf("--------------------- END ---------------------\n");

}
void printStatistics(mscsFunction& f) {
	printf("signal points count: %li\n",f.pointsCount());
	printf("signal mean value: %.15lE\n",f.meanf());
	printf("signal mean value^2: %.15lE\n",pow(f.meanf(),2));
	printf("signal variance: %.15lE\n",f.variance());
	printf("signal stdev: %.15lE\n",f.stdev());
	printf("signal rms: %.15lE\n",f.rms());
	printf("signal rms^2: %.15lE\n",pow(f.rms(),2));
	printf("mind that variance is calculated with N-1 and rms with N\n");
}
/***************************************************************************************/
void extractRanges(mscsFunction& re, mscsFunction& im, mscsFunction& P, cpedsMsgs* msgs) {
	if (_extract_f.size()>0) {
		msgs->say("Extracting requested frequency ranges",High);
		msgs->say("There frequency ranges count: "+msgs->toStr(_extract_f.size()),Low);
		bool zeroRange;
		for (long j = 0; j < P.pointsCount(); j++) {
			cout << j << "/" << P.pointsCount() << "\r";
			zeroRange=true;
			for (long i = 0; i < _extract_f.size(); i++) {
				zeroRange=zeroRange and (P.getX(j)<_extract_f[i]-_extract_w[i % _extract_w.size()] or P.getX(j)>=_extract_f[i]+_extract_w[i % _extract_w.size()]);
			}
			if (zeroRange) {
				re.f(j)=0; 
				im.f(j)=0;
				P.f(j)=0;
			}
			
		}		
		
		for (long i = 0; i < _extract_f.size(); i++) {
			_extractFromTo.newPoint(_extract_f[i]-_extract_w[i % _extract_w.size()],_extract_f[i]+_extract_w[i % _extract_w.size()]);
		}
		_extractFromTo.save(_outfile+".extractedRanges.FromTo");

		msgs->say("Done",Low);
	}
}
/***************************************************************************************/
//void removeRanges(mscsFunction& re, mscsFunction& im, mscsFunction& P, cpedsMsgs* msgs) {
//	msgs->say("Removing requested frequency ranges",High);
//	for (long j = 0; j < P.pointsCount(); j++) {
//		for (long i = 0; i < _removeRanges.pointsCount(); i++) {
//			if (P.getX(j)>=_removeRanges[i].rx() and P.getX(j)<_removeRanges[i].ry()) {
//				re.f(j)=0; 
//				im.f(j)=0;
//				P.f(j)=0;				
//			}
//		}
//	}
//	msgs->say("Done",Low);
//}
