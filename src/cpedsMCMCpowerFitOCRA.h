/*!
* \file cpedsMCMCpowerFitOCRA.h
*
*  Created on: Feb 3, 2011
*      Author: blew
*/

#ifndef CPEDSMCMCPOWERFITOCRA_H_
#define CPEDSMCMCPOWERFITOCRA_H_

#include "cpeds-rng.h"
#include "Mscs-function.h"
#include "cpedsMCMC.h"
#include "cpeds-math.h"

class OCRAPowerMCMC : public cpedsMCMC {
	public:
		OCRAPowerMCMC(long Ndim=0, long id=0, long runOffset=1000, double k0=1.0, long runSeed=0) : cpedsMCMC(Ndim,id,runOffset, runSeed),_simData("sim data",Zero) { 
			msgs->setSender("OCRApower MCMC"); 
			chisqPtr=NULL;
			_x=NULL;
			_X=NULL;
			_noise1of=NULL;
			_noiseG=NULL;
			_noiseSim=NULL;
			_binSize=NULL;
			_Nbins=0;
//			_binnedSim=NULL;
			
			_k0=k0;
			_chisqCalcInLogSpace=true;
		}
		
		~OCRAPowerMCMC() {
//			if (_binnedSim!=NULL) { delete [] _binnedSim; _binnedSim=NULL; }
			cleanCovCalc();
		}
		
		void setCovRNG() {
			_rnCov.setRNsType("gaussian_circle");
		}
		
		/*!
			\brief set if the chisq and covariance matrix should be calculated in linear of logaritmic space for a given model (set of models)
			\details 
			@param tf - if true (default) then every simulation will go into log space before storing for cov calculation and every model will go into log space
			before chisq calculation. If false, then it will not go into log space.
			@return
		
			\date Jul 15, 2011, 3:43:48 PM
			\author Bartosz Lew
		*/
		void chisqCalcInLogSpace(bool tf) {			_chisqCalcInLogSpace=tf;		}
		
		void cleanCovCalc() {
			if (_x!=NULL) { delete [] _x; _x=NULL; }
			if (_X!=NULL) { delete [] _X; _X=NULL; }
			if (_noise1of!=NULL) { delete [] _noise1of; _noise1of=NULL; }
			if (_noiseG!=NULL) { delete [] _noiseG; _noiseG=NULL; }
			if (_noiseSim!=NULL) { delete [] _noiseSim; _noiseSim=NULL; }
			if (_Dcov!=NULL) { delete [] _Dcov; _Dcov=NULL; }
			if (_binSize!=NULL) { delete [] _binSize; _binSize=NULL; }			
		}
		void prepareCovCalc(const cpedsList<long>& binSize) {
			_x=_data.extractArguments(); 
			_X=_simData.extractArguments();			
			_noise1of=new double[_simData.pointsCount()];
			_noiseG=new double[_simData.pointsCount()];
			_noiseSim=new double[_simData.pointsCount()];
			_Dcov=new double[_data.pointsCount()*_NcovSim];
//			_binnedSim=new double[_data.pointsCount()];
			
			// prepare bins
//			cpedsList<long> binSize;
//			long i;
//			double b=_binw;
//			i=_binst;
//			long fsize=_simData.pointsCount();
//			while (i<fsize) {
//				i+=long(round(b));
//				if (i<=fsize) { binSize.append(long(round(b))); } else { binSize.append(fsize-i+long(round(b))); }
//				b*=_bingm;
//			}
			_binSize=binSize.toCarray();
			_Nbins=binSize.count();

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
		
		double modelFunction(double k0, double n, double fknee, double f) const { 
			double k=pow(10,f);
			double knee=pow(10,fknee);
			k0=pow(10,k0);
			return pow( k/k0,n ) + pow(knee/k0,n);
		}
		
		double chisq(const MClink &th) {
//			printf("calling reimplemented chisq\n ****************** data size is : %li\n",_data.pointsCount());
//			return (*chisqPtr)(_data,th);
			if (_NcovSim==0) {
				return chisqNoCov1Log3Lin4Param(_data,th);
			}
			else {
#ifdef COVFAST
				return chisqCovFast(_data,th);
#else
				return chisqCov(_data,th);
#endif
			}
			return 0; // we never return this

//			return chisqNoCovNoLog(_data,th);
//			return chisqNoCov(_data,th);
		}
		
		
		/***************************************************************************************/
		/*!
			\brief derive chisq value with the appropriate covariance matrix
			\details 
			@return chisq value
			
			We assume here that the signal is a gaussian noise. This isn't true for most of the power spectra.
			Only in the gaussian dominated part of the spectra (high frequencies) this becomes justified.
			However calculating the covariance matrix for the correctly simulated non-gaussian data would be VERY
			computationally intensive so here we assume that the data are gaussian. This is more correct in the case of the 
			double-difference power spectra.
			Since we generated the noise in the fourier space as gaussian realization, the covariance matrix is diagonal.
			We do not assume though that all diagonal terms of the covariance matrix are the same. 
			
			We calculate the diagonal terms from a limited number of noise realizations, and average them by binning the same way
			as the data are binned for noise reduction. Such derived covariances are then used for the calculation of the chisq value.
		
			The model is the 4-parameter model as in chisqNoCov1Log3Lin4Param
			
			BLcomment (Jul 8, 2011, 10:11:05 PM): the model is now changed to 3-parameter model since there is an exact degeneracy between k0 and A,
			so k0 is from on on fixed at k0=1 Hz
			
			\date Mar 8, 2011, 2:40:32 PM
			updated on
			\date Mar 8, 2011, 2:46:11 PM
			\author Bartosz Lew
		*/ 
		double chisqCov(const mscsFunction& data, const MClink& th) {
//			msgs->say("Entering chisq calculations",High);
			double *Cov=NULL;
			long vecSize=_simData.pointsCount();
			double *dcov;//=new double[_NcovSim*vecSize];
			_chisqSignature="chisqCov";
			
			//
			// DERIVE THE COVARIANCE MATRIX FOR THIS MODEL
			//

			mscsFunction oof("oof",Zero),gauss("gauss",Zero),sim("sim",Zero), sims("sims",Zero);
			cpedsList<long> bins;
			double A;
			A=pow(10,th[3]); // convert to linear space
			
			// generate simulations for covariance matrix calculation
			for (long j = 0; j < _NcovSim; j++) {
				oof=_simData; // set the correct arguments to calculate the model on
				oof.mkPowerLawNoise(A,th[1],th[0],0,NULL,&_rnCov,false);
				gauss=_simData; // set the correct arguments to calculate the model on
				gauss.mkPowerLawNoise(A*pow(th[2]/th[0],th[1]),0,1,0,NULL,&_rnCov,false);
				sim=oof+gauss; // model in linear space
				// bin simulation
				sim.binFunctionLin2geo(_binst,_binw,_bingm,bins,_dataMask); //here _dataMask should reflect the mask applied to the data
				// convert simulation into a log-log space
//				sim.logX();
				sim.logY();
				sims.concatenate(sim);

				oof.clearFunction();
				gauss.clearFunction();
			}
//			sim.save("sim-test-arguments-values.txt");
			// calculate covariance matrix
			dcov=sims.extractValues();
			vecSize=sim.pointsCount();
			
			Cov=cpeds_calculate_covariance_matrix(dcov,vecSize,_NcovSim,true);
			delete [] dcov;
//			cpedsList<double> ld; ld.fromCarray(Cov,vecSize);
//			sim.setY(ld);
//			sim.save("cov");
			//
			// return the covatiance matrix if needed on cov object
			//
			if (_covDiagonal) {
				if (_covDiag!=NULL) delete [] _covDiag;
				_covDiag=cpeds_copy_array(Cov,vecSize);
			}
			
			//
			// calculate chisq for this model
			//
			msgs->say("calculating chisq",Low);
			
			mscsFunction model;
			oof=_simData; // set the correct arguments to calculate the model on
			A=pow(10,th[3]); // convert to linear space
			oof.mkPowerLaw(A,th[0],th[1]);
			gauss=_simData; // set the correct arguments to calculate the model on
			gauss=double(A*pow(th[2]/th[0],th[1]));
			model=oof+gauss; // model in linear space
			
			// bin data
			model.binFunctionLin2geo(_binst,_binw,_bingm,bins,_dataMask);
			model.logY();
			_model=model;
//			model.save("model-test-arguments-values.txt");
//			exit(0);

			_chisqPerDOF.clear();
			double X2=0;
			long N=model.pointsCount();
			double tmp,x,y;
			for (long i = 0; i < N; i++) {
				x=model.getX(i);
				y=model.f(i);
				tmp=(  data.getY(i) - y  ); tmp*=tmp; tmp/=Cov[i];
				_chisqPerDOF.append(tmp);
				X2+=tmp;
			}
			
			delete [] Cov;
			
			msgs->say("done",Low);
			return X2;
		}
		/***************************************************************************************/
		/*!
			\brief same as chisqCov but without mscsFunction class which in theory should be faster
			\details 
		
				In the BINNED mode there is a tremendous speed-up in calculation of the covariance matrix
				by assuming that the independent fourier components of the power spectrum are independent
				(assumption of gaussianity of the signal). This is btw. not true in the real signal.
				The only quasi-gaussian part of the signal is in the high frequency limit. Low frequencies are
				trully non-gaussian and so modes must be coupled.
				By the assumption of gaussianity we derive the covariance matrix only at the frequencies
				that exist after the binning of the spectra. Then we divide it by the number of points
				that were binned in each bin. This would be correct if the individual non-binned fourier
				components had the same variance. But they don't. ! This is one mistake in this calculation.
				It leads to systematically larger estimates of diagonal covariance matrix terms,
				and possibly only reduces sensitivity.
				
				If the BINNED option is switched off, then the covariance matrix is calculated exactly from
				simulations of the full signal, which are next binned. But this is too slow for large practical runs.

			
			\date Mar 9, 2011, 6:29:42 PM
			\author Bartosz Lew
		*/
		double chisqCovFast(const mscsFunction& data, const MClink& th) {
			double *Cov=NULL;
			long vecSize=_data.pointsCount();
			_chisqSignature="chisqCovfast";
			
//			MClink th(4);
//			th.set(4,9.7717762087E+01,-7.0628572069E-01,8.5888987128E+01,9.4590996741E-01);
			//
			// DERIVE THE COVARIANCE MATRIX FOR THIS MODEL
			//
#ifdef BINNED
			mkDcov2(vecSize,_NcovSim,th);
			Cov=cpeds_calculate_covariance_matrix(_Dcov,vecSize,_NcovSim,true);
			//			for (long i = 0; i < vecSize; i++) { Cov[i]/=_binSize[i]; } // BLmodification (Jul 8, 2011, 9:47:30 PM): commented out this will possibly solve the problem of too small covariances. See /home/blew/programy/test-FFTcoefficients-of-chisq-distributed-variate/readme for more info and tests
#else
			mkDcov(vecSize,_NcovSim,th);			
			Cov=cpeds_calculate_covariance_matrix(_Dcov,vecSize,_NcovSim,true);
#endif
#ifdef DEBUG			
			{
				cpedsList<double> ld; ld.fromCarray(Cov,vecSize);
				mscsFunction sim=_data;
				sim.setY(ld);
				sim.save("cov");
				exit(0);
			}
#endif

			//
			// return the covatiance matrix if needed on cov object
			//
			if (_covDiagonal) {
				if (_covDiag!=NULL) delete [] _covDiag;
				_covDiag=cpeds_copy_array(Cov,vecSize);
			}
			
			//
			// calculate chisq for this model
			//
			msgs->say("calculating chisq",Low);
			
			mscsFunction oof("oof",Zero),gauss("gauss",Zero),sim("sim",Zero);
			mscsFunction model;
			cpedsList<long> bins;
			double A;
			oof=_simData; // set the correct arguments to calculate the model on
			A=pow(10,th[2]); // convert to linear space
			oof.mkPowerLaw(A,_k0,th[0]); // BLmodification (Jul 8, 2011, 10:58:41 PM): change to 3-param model
			gauss=_simData; // set the correct arguments to calculate the model on
			gauss=double(A*pow(th[1]/_k0,th[0]));
			model=oof+gauss; // model in linear space
			
		// bin data
			model.binFunctionLin2geo(_binst,_binw,_bingm,bins,_dataMask);
			if (_chisqCalcInLogSpace) model.logY();
			_model=model;
			//			model.save("model-test-arguments-values.txt");
			//			exit(0);
			
			_chisqPerDOF.clear();
			double X2=0;
			long N=model.pointsCount();
			double tmp,x,y;
			for (long i = 0; i < N; i++) {
				x=model.getX(i);
				y=model.f(i);
				tmp=(  data.getY(i) - y  ); tmp*=tmp; tmp/=Cov[i];
				_chisqPerDOF.append(tmp);
				X2+=tmp;
			}
			
			delete [] Cov;
			
			msgs->say("done",Low);
			return X2;
			
			
			
		}

		
		/***************************************************************************************/
		/*!
		\brief chisq implementation without covariance matrix simulations
		\details 
		@param
		@return
		
		CURRENTLY NOT USED
		\date Feb 3, 2011, 8:21:42 PM
		\author Bartosz Lew
		*/
		double chisqNoCov(const mscsFunction& data, const MClink& th) {
			msgs->say("Entering chisq calculations",High);
			_chisqSignature="chisqNoCov";
			
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
				y=log(modelFunction(th[0], th[1], th[2], x));
				model.newPoint(x,y);
				tmp=(  data.getY(i) - y  ); tmp*=tmp;
				_chisqPerDOF.append(tmp);
				X2+=tmp;
			}
			
			_model=model;
			
			msgs->say("done",Low);
			return X2;
		}
		/***************************************************************************************/		
		/***************************************************************************************/
		/*!
		\brief chisq implementation without covariance matrix simulations
		\details 
		@param
		@return
		This model doesn't go into log-log space with the fiducial model.
		The model is still 3 parameter model:
		P(k) = (k/k0)^n + (knee/k0)^n
		
		CURRENTLY NOT USED
		\date Feb 3, 2011, 8:21:42 PM
		\author Bartosz Lew
		*/
		double chisqNoCovNoLog(const mscsFunction& data, const MClink& th) {
			msgs->say("Entering chisq calculations",High);
			_chisqSignature="chisqNoCovNoLog";
			
			//
			// calculate chisq for this model
			//
			msgs->say("calculating chisq",Low);
			
			mscsFunction oof,gauss,model;
			double dk=1.7732e-3;
			long m=137.0/dk;
			oof.mkPowerLaw(0.01,137,dk,1,th[0],th[1]);
//			gauss.mkPowerLaw(0.01,137,dk,pow(th[2]/th[0],th[1]),1,0);
			gauss.mkConst(0.01,137,dk,pow(th[2]/th[0],th[1]));
			model=oof+gauss;
			model.deletePoint(0); // remove the zero'th k mode
			
			// bin data
			cpedsList<long> bins;
			cpedsList<double> dataMask;
			
			//here load the dataMask. It should be the mask removing the spikes from the signal: TODO use approriate data mask to remove spikes from the signal
			model.binFunctionLin2geo(0,10,1.1,bins,dataMask);
			
			_model=model;

			_chisqPerDOF.clear();
			double X2=0;
			long N=model.pointsCount();
			double tmp,x,y;
			for (long i = 0; i < N; i++) {
				x=model.getX(i);
				y=model.f(i);
				tmp=(  data.getY(i) - y  ); tmp*=tmp;
				_chisqPerDOF.append(tmp);
				X2+=tmp;
			}
			
			
			msgs->say("done",Low);
			return X2;
		}
		
		/***************************************************************************************/
		/*!
		\brief chisq implementation without covariance matrix simulations in linear space 4-parameters
		\details 
		@param
		@return
		This model doesn't go into log-log space with the fiducial model.
		The model is a 4 parameter model:
		P(k) = A [ (k/k0)^n + (knee/k0)^n ]
		
		CURRENTLY NOT USED
		\date Feb 3, 2011, 8:21:42 PM
		\author Bartosz Lew
		*/
		double chisqNoCovNoLog4Param(const mscsFunction& data, const MClink& th) {
			msgs->say("Entering chisq calculations",High);
			_chisqSignature="chisqNoCovNoLog4Param";
			
			//
			// calculate chisq for this model
			//
			msgs->say("calculating chisq",Low);
			
			mscsFunction oof,gauss,model;
			oof=_simData;
			oof.mkPowerLaw(th[3],th[0],th[1]);
			gauss=_simData;
			gauss=double(th[3]*pow(th[2]/th[0],th[1]));
			model=oof+gauss;
			
			// bin data
			cpedsList<long> bins;
			cpedsList<double> dataMask;
			
			//here load the dataMask. It should be the mask removing the spikes from the signal: TODO use approriate data mask to remove spikes from the signal
			model.binFunctionLin2geo(0,10,1.05,bins,dataMask);
			
			_model=model;

			_chisqPerDOF.clear();
			double X2=0;
			long N=model.pointsCount();
			double tmp,x,y;
			for (long i = 0; i < N; i++) {
				x=model.getX(i);
				y=model.f(i);
				tmp=(  data.getY(i) - y  ); tmp*=tmp;
				_chisqPerDOF.append(tmp);
				X2+=tmp;
			}
			
			
			msgs->say("done",Low);
			return X2;
		}
		
		/***************************************************************************************/
		/*!
		\brief chisq implementation without covariance matrix simulations in log-linear space 4-parameters
		\details 
		@param
		@return
		This model doesn't go into log space with the fiducial model for parameters 0,1,2
		Parameter 4 is taken as log.\n
		The model is a 4 parameter model:
		P(k) = log( A [ (k/k0)^n + (knee/k0)^n ] )
		
		\date Feb 3, 2011, 8:21:42 PM
		\author Bartosz Lew
		*/
		double chisqNoCov1Log3Lin4Param(const mscsFunction& data, const MClink& th) {
			msgs->say("Entering chisq calculations",High);
			_chisqSignature="chisqNoCov1Log3Lin4Param";
			
			//
			// calculate chisq for this model
			//
			msgs->say("calculating chisq",Low);
			
			mscsFunction oof,gauss,model;
			double A;
			oof=_simData; // set the correct arguments to calculate the model on
			A=pow(10,th[3]); // convert to linear space
			oof.mkPowerLaw(A,th[0],th[1]);
			gauss=_simData; // set the correct arguments to calculate the model on
			gauss=double(A*pow(th[2]/th[0],th[1]));
			model=oof+gauss; // model in linear space
			
			// bin data
			cpedsList<long> bins;
			
			model.binFunctionLin2geo(_binst,_binw,_bingm,bins,_dataMask);
			model.logY();
			_model=model;

			_chisqPerDOF.clear();
			double X2=0;
			long N=model.pointsCount();
			double tmp,x,y;
			for (long i = 0; i < N; i++) {
				x=model.getX(i);
				y=model.f(i);
				tmp=(  data.getY(i) - y  ); tmp*=tmp;
				_chisqPerDOF.append(tmp);
				X2+=tmp;
			}
			
			
			msgs->say("done",Low);
			return X2;
		}
		
		
		virtual void saveBestFitCovSimulations() {
			cpedsList<long> binSize(_binSize,_Nbins);
			cleanCovCalc();
			_NcovSim=_bestFitNcovSimNum;
			prepareCovCalc(binSize);
			long vecSize=_data.pointsCount();
//#ifdef BINNED
//			mkDcov2(_data.pointsCount(),_NcovSim,_bestFit);
//#else
			mkDcov(_data.pointsCount(),_NcovSim,_bestFit);
//#endif
			
			cpeds_save_matrix(_Dcov,_NcovSim,vecSize,_partialFilesDir+"/bestFit-sims.mc",false);	
		}
		
		

		double _binst; //!< bins size
		double _binw; //!< start binning at
		double _bingm; //!< geometrical sequence multiplier
//		double _dk; //!< use this delta k to generate
		long _NcovSim; //!< number of simulations generated to derive covariance matrix
		mscsFunction _simData; //!< stores arguments at which the initial data were tabulated before binning
	private:
		double (*chisqPtr)(const mscsFunction&, const MClink& );
		
		// PRIVATE OPTIMIZATION C STRUCTURES AND ROUTINES USED TO SPEED UP COV CALCULATION
		double* _x; // binned argument values
		double* _X; // original argument values with gaps
		double* _noise1of;
		double* _noiseG;
		double* _noiseSim;
		double* _Dcov;
		long* _binSize;
		long _Nbins;
		bool _chisqCalcInLogSpace;
//		double* _binnedSim;
		
		double _k0; //!< a k-nod fixed for the 3-parameter model (A,ns,fknee). This is only used on 3-parameter model runs. It is not used otherwise.
		
		
		/*!
			\brief prepare the simulations for covariance matrix calculation
			\details 
			@param vecSize - number of variates in the binned spectra (size of the cov matrix)
			@param vecNum - number of simulations to generate
			@param th - model for which the covariance matrix is to be calculated

			In this approach in each simulation the full input unbinned noise signal is
			generated within gaussian noise model. Then the simulation is binned as was the 
			data and placed in dcov array.
			
			\date Mar 10, 2011, 12:00:38 PM
			\author Bartosz Lew
		*/
		void mkDcov(long vecSize, long vecNum, const MClink& th) {
			double* binnedSim;
			double *xout;
			long Nout;
			double logbase=log(10.0);
			long VecSize=_simData.pointsCount();

			for (long i = 0; i < vecNum; i++) {
				printf("%li\n",i);
				mkPowerLawNoise(th,0); // make a new 1/f realization
				mkPowerLawNoise(th,1); // make a new gaussian realization
				for (long j = 0; j < VecSize; j++) { _noiseSim[j]=_noise1of[j]+_noiseG[j]; 	}
				// bin noise simulation
				binnedSim=cpeds_bin_function(_X,_noiseSim,NULL,_simData.pointsCount(),_binSize,_Nbins,&xout,&Nout,_binst);

				for (long j = 0; j < Nout; j++) { binnedSim[j]=log(binnedSim[j])/logbase; }
				for (long j = 0; j < Nout; j++) { _Dcov[i*vecSize+j]=binnedSim[j]; }	
#ifdef DEBUG
				// debug
				cpeds_save_matrix(binnedSim,Nout,1,string("sim")+msgs->toStr(i),false);
#endif
				delete [] binnedSim;
			}
		}

		/*!
			\brief prepare the simulations for covariance matrix calculation
			\details 
			@param vecSize - number of variates in the binned spectra (size of the cov matrix)
			@param vecNum - number of simulations to generate
			@param th - model for which the covariance matrix is to be calculated


			In this approach we make a significant speed up due to the number of random numbers generated.
			This was the main computational burden in the mkDcov() implementation.
			But since we operate within the gaussian noise model for the covariance matrix calculations anyway,
			a significant speed-up is possible because it is known how the variance of a gaussian noise averages.
			
			So in this implementation we generate the noise simulations that are already binned
			(small number of random numbers to be generated- hence speed-up), and then 
			in the covariance matrix calculation we divide the covariance matrix by the 
			factors that are equal binSize as the variance of the averaged signal over N independent (as they are in our model)
			samples goes down by a factor of N. BLcomment (Jul 15, 2011, 3:35:13 PM): -- we don't do that anymore as this is wrong. 
			

			Important thing is that the data for the Cov calculation (the power spectra values) are in log space.
		
			\date Mar 10, 2011, 12:03:13 PM
			\author Bartosz Lew
		*/
		void mkDcov2(long vecSize, long vecNum, const MClink& th) {
			double logbase=log(10.0);

			for (long i = 0; i < vecNum; i++) {
				mkPowerLawNoise(th,0,1); // make a new 1/f realization
				mkPowerLawNoise(th,1,1); // make a new gaussian realization
				for (long j = 0; j < vecSize; j++) { _noiseSim[j]=_noise1of[j]+_noiseG[j]; 	}

				if (_chisqCalcInLogSpace) for (long j = 0; j < vecSize; j++) { _noiseSim[j]=log(_noiseSim[j])/logbase; }
				for (long j = 0; j < vecSize; j++) { _Dcov[i*vecSize+j]=_noiseSim[j]; }	
#ifdef DEBUG
				cpeds_save_matrix(_noiseSim,vecSize,1,string("sim")+msgs->toStr(i),false);
#endif
			}
		}
		
		/*!
			\brief generate the noise realization fast (don't use the standar function class for that)
			\details 
			@param MCMC link defining the noise properties in the assumed model
			@param noiseType - type of the noise to generate:\n
			0 - 1/f
			 
			1 - gaussian
			
			@param binned - 0 - generate noise on _X grid, 1 - generate noise on _x grid
			@return
		
			\date Mar 9, 2011, 5:47:41 PM
			\author Bartosz Lew
		*/
		void mkPowerLawNoise(const MClink& th, long noiseType, long binned=0) {
			//
			// generate gaussian noise
			//
			long M;
			double* x;
			if (binned==0) {
				M=_simData.pointsCount();
				x=&_X[0];
			}
			if (binned==1) {
				M=_data.pointsCount();
				x=&_x[0];
			}
			
			double *a= _rnCov.getRN(M);
			double *b= _rnCov.getRN(M);
			//
			// form the right power spectrum in Fourier space (positive frequencies only)
			//
			double Pk,sqrPj;
			long j; // iterates k_j; j=0..M-1
			double k;
			double A;
			
			if (dims()==3) {
				if (noiseType==0) { // 1/f
					A=pow(10.0,th[2]);
					for (j=0;j<M;j++) {
						Pk=A*pow( x[j] / _k0,th[0]);
						sqrPj=sqrt(Pk/2.0);
						a[j]*=sqrPj;
						b[j]*=sqrPj;
					}
				}
				if (noiseType==1) { // gaussian
					A=pow(10.0,th[2]);
					Pk=A*pow( th[1] / _k0,th[0]);
					for (j=0;j<M;j++) {
						sqrPj=sqrt(Pk/2.0);
						a[j]*=sqrPj;
						b[j]*=sqrPj;
					}
				}				
			}
			
			if (dims()==4) {
				
				if (noiseType==0) { // 1/f
					A=pow(10.0,th[3]);
					for (j=0;j<M;j++) {
						Pk=A*pow( x[j] / th[0],th[1]);
						sqrPj=sqrt(Pk/2.0);
						a[j]*=sqrPj;
						b[j]*=sqrPj;
					}
				}
				if (noiseType==1) { // gaussian
					A=pow(10.0,th[3]);
					Pk=A*pow( th[2] / th[0],th[1]);
					for (j=0;j<M;j++) {
						sqrPj=sqrt(Pk/2.0);
						a[j]*=sqrPj;
						b[j]*=sqrPj;
					}
				}

			}
			
			
			//
			//	store power spectrum realization
			//
			if (noiseType==0) {
				for (j=0;j<M;j++) {	_noise1of[j]=a[j]*a[j]+b[j]*b[j];  }
			}
			if (noiseType==1) {
				for (j=0;j<M;j++) {	_noiseG[j]=a[j]*a[j]+b[j]*b[j]; }					
			}
			delete [] a;
			delete [] b;
		}
};

#endif /* CPEDSMCMCPOWERFITOCRA_H_ */
