/*!
 * \file test-cpedsMCMC-1ofnoiseV2.cpp
 *
 * see comments from the first version of this program
 * 
 * This program implements binning of the data and simulations to reduce noise and make data distribution more uniform
 * It also implements usage of masks - that remove some data from the analysis (eg. due to spikes in the power spectra)
 *
 *  Created on: Feb 2, 2011
 *      Author: blew
 */

#include "cpeds-rng.h"
#include "Mscs-function.h"
#include "cpedsMCMC.h"

class fitoofPowerSpectrum : public cpedsMCMC {
	public:
		fitoofPowerSpectrum(long Ndim=0, long id=0, long runOffset=0) : cpedsMCMC(Ndim,id, runOffset) { msgs->setSender("OCRApower MCMC"); }
		
		double modelFunction(double k0, double n, double fknee, double f) const { 
			double k=pow(10,f);
			double knee=pow(10,fknee);
			k0=pow(10,k0);
			return pow( k/k0,n ) + pow(knee/k0,n);
		}

		double chisq(const MClink &th) {
			printf("calling reimplemented chisq\n"
					"data size is : %li\n",_data.pointsCount());
			return chisq(_data,th);
		}

		
		virtual double chisq(const mscsFunction& data, const MClink& th) {
			msgs->say("Entering chisq calculations",High);

			
			//
			// DERIVE THE COVARIANCE MATRIX FOR THIS MODEL
			//
			msgs->say("Calculating covariance matrix",Low);

			mscsFunction oof("oof",Zero),gauss("gauss",Zero),sim("sim",Zero);
			double dk=1.7732e-3;
			long m=137.0/dk;
			long n;
			n=data.pointsCount();

			long Nrealizations=5;
			double *Dvec=new double[n*Nrealizations];
			double *Cov;
			long k=0;
			
			for (long i = 0; i < Nrealizations; i++) {
				oof.clearFunction();
				gauss.clearFunction();
				sim.clearFunction();
				sim.mkPowerLawNoise(1,-1,1,0.01,137,m,1,&oof,&_rnCov,false);
				sim.mkPowerLawNoise(1,0,1,0.01,137,m,2,&gauss,&_rnCov,false); 
				sim=oof+gauss;
				sim.deletePoint(0);

				// convert simulation into a log-log space
				sim.logX();
				sim.logY();
				
				for (long j = 0; j < n; j++) {
					Dvec[k]=sim.f(j);
					k++;
				}
			}
			
			
			
			Cov=cpeds_calculate_covariance_matrix(Dvec,n,Nrealizations,true);
			delete [] Dvec;

			sim.clearFunction();
			// bin and re-interpolate to reduce noise in Cov
			sim.importFunction(NULL,Cov,n);
//			sim.save("cov-N_"+msgs->toStr(Nrealizations));

			cpedsList<double> w;//=sim.extractValues();
//			w.invert();
			cpedsList<long> bins;
			sim.binFunctionLin2geo(0,10,1.05,bins, w);
//			sim.save("cov-N_"+msgs->toStr(Nrealizations)+".bin");
			sim.extrapolate(0,n,1,true,"linear");
//			sim.save("cov-N_"+msgs->toStr(Nrealizations)+".bin.extra");
			sim.interpolate(0,n-1,1,true,"linear",true,1e-5);
//			sim.save("cov-N_"+msgs->toStr(Nrealizations)+".bin.extra.inter");
			//
			// return the covatiance matrix if needed on cov object
			//
//			for (long i = 0; i < m; i++) {				cov(i,i)=Cov[i];			}
//			cpeds_save_matrix(Cov,n,1,"Cov20");
//			printf("m: %li\n",n);
//			exit(0);
			for (long i = 0; i < n; i++) {				Cov[i]=sim.f(i);			}
			if (_covDiagonal) {
				if (_covDiag==NULL) _covDiag=new double[n];
				for (long i = 0; i < n; i++) {				_covDiag[i]=Cov[i];			}	
			}
			msgs->say("done",Low);
			
			//
			// calculate chisq for this model
			//
			msgs->say("calculating chisq",Low);
			mscsFunction model;
			_chisqPerDOF.clear();
			double X2=0;
			long N=_data.pointsCount();
			double tmp,x,y;
			for (long i = 0; i < N; i++) {
				x=data.getX(i);
				y=log(modelFunction(th.getParam(0), th.getParam(1), th.getParam(2), x));
				model.newPoint(x,y);
				tmp=(  data.getY(i) - y  ); tmp*=tmp;
				_chisqPerDOF.append(tmp/Cov[i]);
				X2+=tmp/Cov[i];
			}
			
			_model=model;
		
			delete [] Cov;
			msgs->say("done",Low);
			return X2;
		}
		
		
	private:
		// data
//		mscsFunction _data;
};


/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/
int main(int argc, char** argv) {
	
	
//	//
//	// READ THE OBSERVATIONAL DATA
//	//
//	
	mscsFunction oof,gauss,data;
//	data.load(argv[1]);
	
	//
	// GENERATE OBSERVATIONAL DATA
	//
	double dk=1.7732e-3;
	long m=137.0/dk;
	data.mkPowerLawNoise(1,-1,1,0.01,137,m,1,&oof);
	data.mkPowerLawNoise(1,0,1,0.01,137,m,2,&gauss); 
	data=oof+gauss;
	data.deletePoint(0); // remove the zero'th k mode
	data.save("noisePower-data.txt");

	// normalize data to one for large k-modes
	//	data/=double(data.f(data.pointsCount()-1)); // no need to normalize in case of the fiducial data
	
	// convert data into a log-log space
	data.logX();
	data.logY();
	// 
	
	// save the fiducial data
	data.save("noisePower-data.loglog.txt");
	
	//
	// MODEL
	//
	
	
//	// define the model parameter and priors
//	/*The model: 
//	 * The model of the data is a piece-wise combination of two linear functions that represent 
//	 * two limiting cases of noise domains: the 1/f dominated part of the spectrum and 
//	 * gaussian fluctuations dominated part of the spectrum i.e.:
//	 * p(k) = A k + B ; k <= fknee and 
//	 * p(k) = C       ; k > fknee and 
//	 * 
//	 * But A fknee + B == b2 so 
//	 * the model has 3 parameters: 
//	 * A, B, and fknee.
//	 * Typically the values of interests lie in ranges:
//	 * A in (-3,0)
//	 * B in (0,10)
//	 * fknee in (-2.5, 2.5) */
//		
//	mscsFunction A("A"); 
//	mscsFunction B("B"); 
//	mscsFunction fknee("fknee"); 
//
//	// make a flat prior on parameter A
//	A.mkConst(-3.0,0.0,0.01,1.0);
//	A.normalize();
//	// make a flat prior on parameter B
//	B.mkConst(0.0,10.0,0.01,1.0);
//	B.normalize();
//	// make a flat prior on parameter fknee
//	fknee.mkConst(-2.5,2.5,0.01,1.0);
//	fknee.normalize();

	
	// define a new model parameter and priors
	/*The model: 
	 * The model of the data is a combination of two functions that represent 
	 * two limiting cases of noise domains: the 1/f dominated part of the spectrum and 
	 * gaussian fluctuations dominated part of the spectrum i.e.:
	 * p(k) = (k/k0)^n + G ;
	 * 
	 * But (fknee/k0) + G == 0 so 
	 * the model has 3 parameters: 
	 * k0, n, and fknee.
	 * It can be rewriten as:
	 * p(k) = (k/k0)^n + (fknee/k0)^n ;
	 * Typically the values of interests lie in ranges:
	 * k0 in (-2.5,2.5)
	 * n in (-3,0)
	 * fknee in (-2.5, 2.5) */

	mscsFunction k0("k0"); 
	mscsFunction n("n"); 
	mscsFunction fknee("fknee"); 

	// make a flat prior on parameter A
	k0.mkConst(-2.5,2.5,0.01,1.0);
	k0.normalize();
	// make a flat prior on parameter B
	n.mkConst(-3.0,0.0,0.01,1.0);
	n.normalize();
	// make a flat prior on parameter fknee
	fknee.mkConst(-2.5,2.5,0.01,1.0);
	fknee.normalize();

	
	//
	// define the covariance matrix for the input data within the defined model
	//
	
	// number of simulations
	
	
	//
	// FIT 
	//
	
	// set the dimensionality of the problem
	long idx=strtol(argv[1],NULL,10);
	fitoofPowerSpectrum fitSpectra(3,idx);

	
	// add parameters
//	fitSpectra.addParameter(A);
//	fitSpectra.addParameter(B);
//	fitSpectra.addParameter(fknee);
	fitSpectra.addParameter(k0);
	fitSpectra.addParameter(n);
	fitSpectra.addParameter(fknee);

//	fitSpectra.saveStepPDFs("stepPDFini");
	//
	// set the observational data
	//
	printf("data size is :%li\n",data.pointsCount());
	fitSpectra.setData(data);
	printf("data size is :%li\n",data.pointsCount());
	
	fitSpectra.setSaveParameters(20,20,20);
	fitSpectra.setDumpEveryMCMCstate(true);
	fitSpectra.setCoolingRate(1);
	fitSpectra.setBurnInLength(100);
	fitSpectra.setCovarianceMatrixDiagonal(true);
	//
	// RUN THE MCMC CHAIN
	//
	fitSpectra.startChain();
	fitSpectra.saveResults();
	
	//
	// SAVE THE CHAIN TO FILE
	//
	
	// save the chain of chisq values 
//	fitSpectra.saveChisq("chisq.mc");
//	fitSpectra.saveParams("params.mc");
//	fitSpectra.saveTemperature("temperature.mc");
//	fitSpectra.saveStepPDFs("stepPDFend");
	
	return 0;
}





