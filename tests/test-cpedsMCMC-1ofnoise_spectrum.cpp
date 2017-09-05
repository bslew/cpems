/*
 * This is a test program to test the cpeds MCMCs and to find the optimal solution for fiting the power spectra
 * comming from the ocra-f raw data.
 * 
 * The explored problems are as follows:
 * - the unequal distribution of points in the power spectra makes the chisq statistic insensitive to the 1/f part of the spectra resulting in a poor fit, especially for low fknee 
 * - the chisq takes most of its value from the high frequency region which is noisier but has thousands times more points than low freq. region.
 * - the fits are generally poor or bad for the assumed 3-parameter model even if the fit is done in log-log space
 * - the covariance matrix isn't noisy (its smooth) but the data is noisy and the resulting fit is poor.
 * - in this implementation the covariance matrix changes from model to model
 * 
 * - due to the weak dependence of the chisq on the 1/f part, some fits tend to find solutions where fknee=k0 which is consistent with the input data
 * but individual values are wrong. This is because 1/f part is insignificant to the chisq estimator so that apropriate k0 value cannot be estabilished
 * 
 * POSSIBLE SOLUTIIONS:
 * - estabilished a model with different parametrization and/or
 * - bin data to reduce noise and reduce the effects of unequal distribution
 * 
 * This will be addressed in test-cpedsMCMC-1ofnoiseV2.cpp program
 *  
 * BUGS:
 * - the generated simulated data for cov calculation were actually not dependent on the model (Feb 3, 2011, 11:23:59 AM) correcting the file
 * 
 * 
 * CHANGELOG:
 * 
 * 
 * CREATED ON:
 * Feb 3, 2011, 8:55:26 AM
 * 
 */

#include "cpedsMCMCpowerFitOCRA.h"
#ifdef GOMPI
#include "mpi.h"
#endif

/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/
int main(int argc, char** argv) {
	int ierr,node_num=1;
	int world_id=0;
	
//#ifdef GOMPI
//	MPI_Status mpi_stat;
//#endif  
	
#ifdef GOMPI
	ierr = MPI_Init ( &argc, &argv ); // initiate MPI
	ierr = MPI_Comm_size ( MPI_COMM_WORLD, &node_num ); // Get the number of processes
	ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &world_id ); // Get the individual process ID
#endif
	
	
	mscsFunction oof,gauss,data;
	//	data.load(argv[1]);
	
	//
	// GENERATE OBSERVATIONAL DATA
	//
	double dk=1.7732e-3;
	long m=137.0/dk;
//	data.mkPowerLawNoise(1,-1,1,0.01,137,m,3,&oof);
//	data.mkPowerLawNoise(1,0,1,0.01,137,m,4,&gauss); 
	oof.mkPowerLaw(0.01,137,dk,1,1,-1);
//	gauss.mkPowerLaw(0.01,137,dk,1,1,0);
	gauss.mkConst(0.01,137,dk,1);
	data=oof+gauss;
	data.deletePoint(0); // remove the zero'th k mode
	if (world_id==0) {
		data.save("noisePower-data.txt");
	}
	
	// bin data
	cpedsList<long> bins;
	cpedsList<double> dataMask;
	
	//here load the dataMask. It should be the mask removing the spikes from the signal: TODO use approriate data mask to remove spikes from the signal
	data.binFunctionLin2geo(0,10,1.1,bins,dataMask);
	
	// normalize data to one for large k-modes
	//	data/=double(data.f(data.pointsCount()-1)); // no need to normalize in case of the fiducial data
	
	// convert data into a log-log space
//	data.logX();
//	data.logY();
	// 
	
	// save the fiducial data
	if (world_id==0) {
		data.save("noisePower-data.loglog.txt");
	}	
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
	 * But fknee is defined as (fknee/k0)^n == G so 
	 * the model has 3 parameters: 
	 * k0, n, and fknee.
	 * It can be rewriten as:
	 * p(k) = (k/k0)^n + (fknee/k0)^n ;
	 * Typically the values of interests lie in ranges (in log space for k0 and fknee):
	 * k0 in (-2.5,2.5)
	 * n in (-3,0)
	 * fknee in (-2.5, 2.5) */
	
	mscsFunction k0("k0"); 
	mscsFunction n("n"); 
	mscsFunction fknee("fknee"); 
	
//	// make a flat prior on parameter A
//	k0.mkConst(-2.5,2.5,0.01,1.0);
//	k0.normalize();
//	// make a flat prior on parameter B
//	n.mkConst(-3.0,0.0,0.01,1.0);
//	n.normalize();
//	// make a flat prior on parameter fknee
//	fknee.mkConst(-2.5,2.5,0.01,1.0);
//	fknee.normalize();
	
//	// make a flat prior on parameter A
//	k0.mkConst(0,137,0.01,1.0);
//	k0.normalize();
//	// make a flat prior on parameter B
//	n.mkConst(-3.0,0.0,0.01,1.0);
//	n.normalize();
//	// make a flat prior on parameter fknee
//	fknee.mkConst(0,137,0.01,1.0);
//	fknee.normalize();

	// make a flat prior on parameter A
	k0.mkConst(0,3,0.01,1.0);
	k0.normalize();
	// make a flat prior on parameter B
	n.mkConst(-3.0,0.0,0.01,1.0);
	n.normalize();
	// make a flat prior on parameter fknee
	fknee.mkConst(0,3,0.01,1.0);
	fknee.normalize();
	
	//
	// define the covariance matrix for the input data within the defined model
	//
	
	// number of simulations
	
	
	//
	// FIT 
	//
	
	// set the dimensionality of the problem
	cpedsMsgs msgs;
	long idx=strtol(argv[1],NULL,10);
#ifdef GOMPI
	idx=world_id;
#endif	
	OCRAPowerMCMC fitSpectra(3,idx,1000);
	fitSpectra.msgs->setSender("MCMC."+msgs.toStr(world_id));
	if (node_num>1) {
		fitSpectra.setOutputDir(msgs.toStr("mcmc.")+msgs.toStr(argv[1]));
	}
	else 
		fitSpectra.setOutputDir("mcmc");
	
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
	
	fitSpectra.setSaveParameters(0,0,0);
	fitSpectra.setDumpEveryMCMCstate(true);
	fitSpectra.setCoolingRate(1);
	fitSpectra.setBurnInLength(1000);
	fitSpectra.setMaximalRejectionsCount(50);
	fitSpectra.setCovarianceMatrixDiagonal(true);
	fitSpectra.setDataMask(dataMask);
	fitSpectra.setMaximalChainLength(20000);
	fitSpectra.setUpHillClimbing(true);
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
	
	
#ifdef GOMPI
	ierr = MPI_Finalize ( ); // Shut down MPI.
#endif
	
	return 0;
}





