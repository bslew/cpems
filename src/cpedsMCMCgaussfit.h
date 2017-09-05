/*!
 * \file cpedsMCMCgaussfit.h
 * 
 * test cpeds MCMC with 3 parameters. fitting gauss data
 *
 *  Created on: Feb 3, 2011
 *      Author: blew
 */

#ifndef CPEDSMCMCGAUSSFIT_H_
#define CPEDSMCMCGAUSSFIT_H_


#include "cpeds-rng.h"
#include "Mscs-function.h"
#include "cpedsMCMC.h"

class gaussMCMC : public cpedsMCMC {
	public:
		gaussMCMC(long Ndim=0, long id=0) : cpedsMCMC(Ndim,id) { msgs->setSender("gauss MCMC"); }
		
		double modelFunction(double A, double m, double s, double x) const { 
			return A*exp(-((x-m)*(x-m))/(2*s*s));
		}

		double chisq(const MClink &th) {
			msgs->say("calling reimplemented chisq, data size is: "+msgs->toStr(_data.pointsCount()),High);
			return chisq(_data,th);
		}

		
		virtual double chisq(const mscsFunction& data, const MClink& th) {
			msgs->say("Entering chisq calculations",High);
			
			double *Cov=NULL;
//			//
//			// DERIVE THE COVARIANCE MATRIX FOR THIS MODEL
//			//
//			msgs->say("Calculating covariance matrix",Low);
//
//			mscsFunction oof("oof",Zero),sim("sim",Zero);
//			double dk=1.7732e-3;
//			long m=137.0/dk;
//			long n;
//			n=data.pointsCount();
//
//			long Nrealizations=2;
//			double *Dvec=new double[n*Nrealizations];
//			long k=0;
//			
//			for (long i = 0; i < Nrealizations; i++) {
//				oof.clearFunction();
//				sim.clearFunction();
//				sim.mkPowerLawNoise(1,-1,1,0.01,137,m,1,&oof,&_rnCov,false);
//				sim=oof;
//				sim.deletePoint(0);
//
//				// convert simulation into a log-log space
//				sim.logX();
//				sim.logY();
//				
//				for (long j = 0; j < n; j++) {
//					Dvec[k]=sim.f(j);
//					k++;
//				}
//			}
//			
//			
//			
//			Cov=cpeds_calculate_covariance_matrix(Dvec,n,Nrealizations,true);
//			delete [] Dvec;
//
//			sim.clearFunction();
//			// bin and re-interpolate to reduce noise in Cov
//			sim.importFunction(NULL,Cov,n);
////			sim.save("cov-N_"+msgs->toStr(Nrealizations));
//
//			cpedsList<double> w;//=sim.extractValues();
////			w.invert();
//			cpedsList<long> bins;
//			sim.binFunctionLin2geo(0,10,1.05,bins, w);
//			sim.extrapolate(0,n,1,true,"linear");
//			sim.interpolate(0,n-1,1,true,"linear",true,1e-5);
//			//
//			// return the covatiance matrix if needed on cov object
//			//
//			for (long i = 0; i < n; i++) {				Cov[i]=sim.f(i);			}
//			if (_covDiagonal) {
//				if (_covDiag==NULL) _covDiag=new double[n];
//				for (long i = 0; i < n; i++) {				_covDiag[i]=Cov[i];			}	
//			}
//			msgs->say("done",Low);
			
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
				y=modelFunction(th.getParam(0), th.getParam(1),th.getParam(2), x);
				model.newPoint(x,y);
				tmp=(  data.getY(i) - y  ); tmp*=tmp;
//				_chisqPerDOF.append(tmp/Cov[i]);
				_chisqPerDOF.append(tmp);
//				X2+=tmp/Cov[i];
				X2+=tmp;
			}
			
			_model=model;
		
			if (Cov!=NULL) delete [] Cov;
			msgs->say("done",Low);
			return X2;
		}
		
		
	private:
};


#endif /* CPEDSMCMCGAUSSFIT_H_ */
