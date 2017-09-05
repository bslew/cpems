#include "cpeds-rng.h"
#include "Mscs-function.h"
#include "cpedsMCMC.h"
#include "cpeds-math.h"


class parabola2dMCMC : public cpedsMCMC {
	public:
		parabola2dMCMC(long Ndim=0, long id=0, long runOffset=1000) : cpedsMCMC(Ndim,id,runOffset) { 
			msgs->setSender("parabola2dMCMC"); 
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
			_rnCov.setRNsType("uniform");
			_rnCov.setMinMax(-_acc,_acc);			
		}
				
		double chisq(const MClink &th) {
			printf("calling reimplemented chisq\n ****************** data size is : %li\n",_data.pointsCount());
//			return (*chisqPtr)(_data,th);
//			return chisqNoCov1Log3Lin4Param(_data,th);
//			return chisqNoCovNoLog(_data,th);
//			return chisqNoCov(_data,th);
			return chisqCovDiag(_data,th);
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
				y=th[0]*x*x+th[1];
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
		/*!
		\brief chisq implementation with diagonal covariance matrix simulations.
		\details 
		@param
		@return
		
		\date Feb 3, 2011, 8:21:42 PM
		\author Bartosz Lew
		*/
		double chisqCovDiag(const mscsFunction& data, const MClink& th) {
			msgs->say("Entering chisq calculations",High);
			_chisqSignature="chisqCovDiag";
			
			
			//
			// derive covariance matrix for this calculation
			//
			long N=100;
			long vecSize=_data.pointsCount();
			double* dcov=new double[N*vecSize];
			double* x=_data.extractArguments();
			long j;
			for (long i = 0; i < N; i++) {	
				j=i*vecSize;
				for (long k = 0; k < vecSize; k++) {		
					dcov[j+k]=th[0]*x[k]*x[k]+th[1]+_rnCov.getRN();
				}
			}
			double* cov=cpeds_calculate_covariance_matrix(dcov,vecSize,N,true);
			delete [] dcov;
			double covav=cpeds_sum(cov,vecSize,true);
//			mscsFunction covint;
//			covint.importFunction(NULL,cov);
//			covint.binFunction
//			cpeds_save_matrix(cov,vecSize,1,"cov",false);
//			cpeds_save_matrix(dcov,N,vecSize,"dcov",false);
//			exit(0);
			//
			// calculate chisq for this model
			//
			msgs->say("calculating chisq",Low);
			mscsFunction model;
			_chisqPerDOF.clear();
			double X2=0;
			double tmp,y;
			for (long i = 0; i < vecSize; i++) {
				y=th[0]*x[i]*x[i]+th[1];
				model.newPoint(x[i],y);
				tmp=(  data.getY(i) - y  ); tmp*=tmp; tmp/=covav;
				_chisqPerDOF.append(tmp);
				X2+=tmp;
			}
			
			_model=model;
			
			delete [] x;
//			delete [] cov;
			msgs->say("done",Low);
			return X2;
		}
		

		double _acc; //!< measurements accuracy
//		double _binw; //!< start binning at
//		double _bingm; //!< geometrical sequence multiplier
//		mscsFunction _simData; //!< stores arguments at which the initial data were tabulated before binning
	private:
		double (*chisqPtr)(const mscsFunction&, const MClink& );
};

//class parabolaMCMC : public cpedsMCMC {
//	public:
//		parabolaMCMC() { msgs->setSender("parabolaMCMC"); }
//		parabolaMCMC(long Ndim) : cpedsMCMC(Ndim) { msgs->setSender("parabolaMCMC"); }
//		parabolaMCMC(mscsFunction& data) { _data=data; }
//		
//		//		virtual double chisq(matrix<double>& data, const MClink &th, matrix<double>& cov) {
//		virtual double chisq(const mscsFunction& data, const MClink& th,  matrix<double>& cov) {
////			msgs->say("calling the right chisq function",High);
////			th.printParams();
////			msgs->say("a parameter is: "+msgs->toStr(th[0]),High);
//			
////			if (_data.pointsCount()==0) setData(data);
//			//			
//			//
//			// define the model data
//			//
//			mscsFunction model;
//			model.mkPowerLaw(-1,1,0.01,th[0],1,2);
//			_model=model;
////			model.print();
////			printf("a is: %lE\n",th[0]);
//			//
//			// calculate chisq for this model
//			//
//			mscsFunction X=data-model;
//			X.power(2.0);
////			X.print();
////			matrix<double> x=X.valuesToVector(false);
////			x.
////			printf("chisq is: %lE\n",cpeds_sum(X.extractValues(),X.pointsCount(),true));
////			exit(0);
//			return cpeds_sum(X.extractValues(),X.pointsCount(),true);
//			
////			double X=0;
////			long N=_data.pointsCount();
////			double *y=_data.extractValues();
////			double var=cpeds_variance(y,N);
////			double mean=cpeds_mean_value(y,N);
////			
////			for (long i = 0; i < N; i++) {
////				y[i]-=link[0]*_data.getX(i)*_data.getX(i);
////			}
////			
////			for (long i = 0; i < N; i++) {				X+=y[i]*y[i];			}
////			
////			return X;
//		}
//		
//		
//	private:
//		// data
//		mscsFunction _data;
//};


/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/
int main() {
	
	
	//
	// PREPARE THE OBSERVATIONAL DATA
	//
	
	mscsFunction data;
	data.mkPowerLaw(-1,1,0.01,1,1,2);

	cpedsRNG rns;
	rns.setRNsType(rns.uniform);
	rns.seed(1);
	double dev;
//	dev=0.5;
	dev=0.01;
	rns.setMinMax(-dev,dev);
	double B=1;
	for (long i = 0; i < data.pointsCount(); i++) {		data.f(i)+=rns.getRN()+B;	}
	
	long idx=3;
	cpedsMsgs msgs;
	
	for (long i = 0; i < 10; i++) {
		{
			// set the dimensionality of the problem
			parabola2dMCMC fitparabola(2,i,1000);
			fitparabola.setOutputDir(string("dev_")+msgs.toStr(dev)+msgs.toStr("-mcmc.")+msgs.toStr(idx));
			
			fitparabola._acc=dev;
			fitparabola.setCovRNG();
			
			
			// define the parameter priors
			mscsFunction a("a");
			mscsFunction b("b");
			
			// make a flat prior on parameter a
			a.mkConst(-2.0,3.0,0.01,1.0);
			a.normalize();
			b.mkConst(-2.0,3.0,0.01,1.0);
			b.normalize();
			
			fitparabola.addParameter(a);
			fitparabola.addParameter(b);
			
			
			fitparabola.setData(data);
			//
			// set MCMC parameters
			//
			fitparabola.setSaveParameters(0,0,0);
			fitparabola.setDumpEveryMCMCstate(false);
			//		fitparabola.setCoolingRate(0.1); // very slow cooling rate
			fitparabola.setCoolingRate(0.5); // medium cooling rate 
//			fitparabola.setCoolingRate(1.0); // normal cooling rate
			//	fitparabola.setCoolingRate(_coolingRate); // normal cooling rate
			fitparabola.setBurnInLength(1000);
			fitparabola.setMaximalRejectionsCount(50);
			fitparabola.setCovarianceMatrixDiagonal(true);
			//	fitparabola.setDataMask(dataMask);
			fitparabola.setMaximalChainLength(50000);
			fitparabola.setUpHillClimbing(true);

			fitparabola.addImplementationNote("");
			fitparabola.addImplementationNote(string("DATA OPERATION OPTIONS: "));
			fitparabola.addImplementationNote(string("uniformly distributed deviation from the model set in the data (dev): ")+msgs.toStr(dev));
			fitparabola.addImplementationNote(string("PROGRAM NAME: test-parabolaMCMC-2d.cpp: "));
			
			
			//
			// RUN THE MCMC CHAIN
			//
			fitparabola.startChain();
			fitparabola.saveResults();
		}
	}

	return 0;
}





