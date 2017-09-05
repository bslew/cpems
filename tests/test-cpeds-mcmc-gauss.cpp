/*!
 * \file test-cpeds-mcmc-gauss.cpp
 *
 *  Created on: Feb 6, 2012
 *      Author: blew
 */

#include "cpeds-rng.h"
#include "Mscs-function.h"
#include "cpedsMCMC.h"
#include "cpeds-math.h"


class gaussMCMC : public cpedsMCMC {
	public:
		gaussMCMC(long Ndim=0, long id=0, long runOffset=1000) : cpedsMCMC(Ndim,id,runOffset) { 
			msgs->setSender("gauss2dMCMC"); 
			chisqPtr=NULL;
		}
		
		/*!
			\brief defines appropriate chisq function to call to evaluate chisq value for a given model
			\details 
			@param name of the model

			\date Feb 8, 2011, 11:05:48 AM
			\author Bartosz Lew
		*/
		void setChisqModel(string model) {
//			if (model=="chisqCov") chisqPtr=&chisqCov;
//			if (model=="chisqNoCov") chisqPtr=&chisqNoCov;
//			if (model=="chisqNoCovNoLog") chisqPtr=&chisqNoCovNoLog;
//			if (chisqPtr==NULL) { msgs->criticalError("wrong chisq address. Please check the code",High); }
		}
		
		void setCovRNG() {
			_rnCov.setRNsType("gaussian_circle");
		}
				
		double chisq(const MClink &th) {
//			printf("calling reimplemented chisq\n ****************** data size is : %li\n",_data.pointsCount());
//			return (*chisqPtr)(_data,th);
//			return chisqNoCov1Log3Lin4Param(_data,th);
//			return chisqNoCovNoLog(_data,th);
			return chisqMeanSigma(_data,th);
//			return chisqCovDiag(_data,th);
		}
		
		
		/***************************************************************************************/
		/*!
		\brief chisq implementation without covariance matrix simulations
		\details 
		@param
		@return
		
		\date Feb 3, 2011, 8:21:42 PM
		\author Bartosz Lew
		*/
		double chisqMeanSigma(const mscsFunction& data, const MClink& th) {
//			msgs->say("Entering chisq calculations",High);
			_chisqSignature="chisqMeanSigma";
			
			
			//
			// calculate chisq for this model
			//
			mscsFunction model;
			_chisqPerDOF.clear();

			model.newPoint(1,th[0]);
			model.newPoint(2,th[1]);
			double X2=0;

			X2=pow(_data.f(0)-model.f(0),2)/getCov(0) + pow(_data.f(1)-model.f(1),2)/getCov(1);
			
			_model=model;
			
			return X2;
		}

		/***************************************************************************************/
		/*!
		\brief chisq implementation with diagonal covariance matrix simulations.
		\details 
		@param
		@return
		
		\date Feb 3, 2011, 8:21:42 PM
		\author Bartosz Lew
		*/
		double chisqCovDiag(const mscsFunction& data, const MClink& th) {
			msgs->say("Entering chisq calculations",Low);
			_chisqSignature="chisqCovDiag";
			
			
			//
			// derive covariance matrix for this calculation
			//
			long vecSize=_data.pointsCount();
			long vecNum=100; // we want 100x more measurements sets than the length of the data set
			long N=vecSize*vecNum; 
//			double* dcov=new double[N * vecSize];
			double* x=_data.extractArguments();
			_rnCov.setMeanVariance(th[0],th[1]);			
			double* dcov=_rnCov.getRN(N);
			double* sim=_rnCov.getRN(vecSize);

			double* cov=cpeds_calculate_covariance_matrix(dcov,vecSize,vecNum,true);
			delete [] dcov;
//			double covav=cpeds_sum(cov,vecSize,true);
			//
			// calculate chisq for this model
			//
			msgs->say("calculating chisq",Low);
			mscsFunction model;
			_chisqPerDOF.clear();
			double X2=0;
			double tmp,y;
			for (long i = 0; i < vecSize; i++) {
				y=sim[i];
				model.newPoint(x[i],y);
				tmp=(  data.getY(i) - y  ); tmp*=tmp; tmp/=cov[i];
				_chisqPerDOF.append(tmp);
				X2+=tmp;
			}
			
			_model=model;
			
			delete [] x;
			delete [] cov;
			delete [] sim;
			msgs->say("done",Low);
			return X2;
		}
		

//		double _acc; //!< measurements accuracy
//		double _binw; //!< start binning at
//		double _bingm; //!< geometrical sequence multiplier
//		mscsFunction _simData; //!< stores arguments at which the initial data were tabulated before binning
	private:
		double (*chisqPtr)(const mscsFunction&, const MClink& );
};



/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/
int main() {
	
	
	//
	// PREPARE THE OBSERVATIONAL DATA
	//

	cout << "\n";
	cout << "***************************************************\n";
	cout << "Generating random sample from gaussian distribution\n";

	cpedsRNG rns;
	rns.setRNsType(rns.gaussian_circle);
	double model_mean=0;
	double model_variance=1;
	long sampleSize=10000;
	rns.setMeanVariance(model_mean,model_variance);
	rns.seed(1);
	cpedsList<double> sample=rns.getRNs(sampleSize);
	mscsFunction data("data");
	data.newPoint(1,sample.mean());
	data.newPoint(2,sample.variance());

	cout << "theoretical mean: " << model_mean << "\n";
	cout << "theoretical variance: " << model_variance << "\n";
	cout << "sample size: " << sampleSize << "\n";

	//
	// test covariance
	//
	cout << "\n";
	cout << "***************************************************\n";
	cout << "Calculating sample covariance matrix\n";
	long Nsim=1000;
	cout << "Number of simulations: " << Nsim << "\n";

	cpedsList<double> sim_m;
	cpedsList<double> sim_s;
	for (long i = 0; i < Nsim; i++) {
		cpedsList<double> sim=rns.getRNs(sampleSize);
		sim_m.append(sim.mean());
		sim_s.append(sim.variance());
	}
	cpedsList<double> cov;
	cov.append(sim_m.variance());
	cov.append(sim_s.variance());
	cout << "Calculated diagonal terms of the covariance matrix\n";
	cov.print();

	cout << "Theoretical diagonal terms of the covariance matrix\n";
	double variance_of_mean=1.0/sampleSize;
	double variance_of_variance=pow(sampleSize-1,2)/pow(sampleSize,3)*3-(sampleSize-1)*(sampleSize-3)*model_variance/pow(sampleSize,3);
	cout << variance_of_mean << "\n";
	cout << variance_of_variance << "\n";

	cout << "Using theoretical covariance matrix\n";
	cov[0]=variance_of_mean;
	cov[1]=variance_of_variance;
	

	//
	// Starting MCMC
	//
	cout << "\n";
	cout << "***************************************************\n";
	cout << "Starting MCMC\n";
	cout << "Search parameters: mean, variance\n";
	long idx=3;
	cpedsMsgs msgs;
	cpedsList<double> bestFitValues;
	QList<mscsFunction> bestFitPosteriors;
	for (long i = 7; i < 8; i++) {
		{
			// set the dimensionality of the problem
			gaussMCMC fitgauss(2,i,1000);
			fitgauss.setOutputDir("gauss2d/mcmc");
			
//			fitgauss._acc=dev;
			fitgauss.setCovRNG();
			fitgauss.setCovarianceMatrixDiagonal(true);
			fitgauss.setDiagonalCovarianceMatrix(cov);
			
			
			// define the parameter priors
			mscsFunction m("mean");
			mscsFunction s("sigma");
			
			// make a flat prior on parameter a
			m.mkConst(-3.0,3.0,0.01,1.0);
			m.normalize();
			s.mkConst(1e-2,30.0,0.01,1.0);
			s.normalize();
			
			fitgauss.addParameter(m);
			fitgauss.addParameter(s);
			
			
			fitgauss.setData(data);
			//
			// set MCMC parameters
			//
			fitgauss.setSaveParameters(0,0,0);
			fitgauss.setDumpEveryMCMCstate(false);
			fitgauss.setInitialStepSize(0.5);
			fitgauss.setInitialWalkStepSize(0.05);
			fitgauss.setCoolingRate(0.1); // medium cooling rate 
			fitgauss.setTemperatures(1000,1);
			fitgauss.setBurnInLength(20000);
			fitgauss.setMaximalRejectionsCount(100);
			fitgauss.setCovarianceMatrixDiagonal(true);
			fitgauss.setMaximalChainLength(1000000);
			fitgauss.setUpHillClimbing(true);

			fitgauss.addImplementationNote("");
			fitgauss.addImplementationNote(string("DATA OPERATION OPTIONS: "));
			fitgauss.addImplementationNote(string("gaussian distributed measurements with stdev=: ")+msgs.toStr(model_variance));
			fitgauss.addImplementationNote(string("PROGRAM NAME: test-cpeds-mcmc-gauss.cpp: "));
			
			
			//
			// RUN THE MCMC CHAIN
			//
			fitgauss.startChain();
			fitgauss.saveResults();
			bestFitPosteriors=fitgauss.posteriors();
			bestFitValues=fitgauss.bestFitLink().getParameters();
		}
	}

	cout << "\n";
	cout << "***************************************************\n";
	cout << "Provided sample parameters:\n";
	cout << "mean: " << sample.mean() << "\n";
	cout << "variance: " << sample.variance() << "\n";

	cout << "\n";
	cout << "***************************************************\n";
	cout << "Reconstructed parameters:\n";
	cout << "mean: " << bestFitValues[0] << "\n";
	cout << "variance: " << bestFitValues[1] << "\n";

	cout << "\n";
	cout << "***************************************************\n";
	cout << "Theoretical width of the reconstructed PDF:\n";
	double sigma2fwhm=2.0*sqrt(2.0*log(2.0));
	cout << "width(mean): " << sqrt(variance_of_mean)*sigma2fwhm << "\n";
	cout << "width(variance): " << sqrt(variance_of_variance)*sigma2fwhm << "\n";

	cout << "\n";
	cout << "***************************************************\n";
	cout << "Width of reconstructed PDFs:\n";
	bestFitPosteriors[0]-=0.5;
	vector<double> v;
	v=bestFitPosteriors[0].findRoot();
	double width_mean=v[1]-v[0];
	bestFitPosteriors[1]-=0.5;
	v=bestFitPosteriors[1].findRoot();
	double width_variance=v[1]-v[0];
	cout << "width(mean): " << width_mean << "\n";
	cout << "width(variance): " << width_variance << "\n";
	return 0;
}





