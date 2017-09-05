#include "cpedsMCMCgaussfit.h"

int main(int argc, char** argv) {
	
	
//	//
//	// READ THE OBSERVATIONAL DATA
//	//
//	
	mscsFunction oof,data;
//	data.load(argv[1]);
	
	//
	// GENERATE OBSERVATIONAL DATA
	//
	data.mkGauss(-3,3,0.1,1,0,1);
	data.save("gauss-data.txt");

	
	//
	// MODEL
	//
	
	

	mscsFunction A("A"); 
	mscsFunction m("m"); 
	mscsFunction s("s"); 

	// make a flat prior on parameter k0
	A.mkConst(-10,10,0.1,1.0);
	A.normalize();
	// make a flat prior on parameter m
	m.mkConst(-10.0,10.0,0.1,1.0);
	m.normalize();
	// make a flat prior on parameter s
	s.mkConst(0,10.0,0.1,1.0);
	s.normalize();
	
	//
	// define the covariance matrix for the input data within the defined model
	//
	
	// number of simulations
	
	
	//
	// FIT 
	//
	
	// set the dimensionality of the problem
	long idx=strtol(argv[1],NULL,10);
	gaussMCMC fitGauss(3,idx);

	
	// add parameters
	fitGauss.addParameter(A);
	fitGauss.addParameter(m);
	fitGauss.addParameter(s);

//	fitSpectra.saveStepPDFs("stepPDFini");
	//
	// set the observational data
	//
	printf("data size is :%li\n",data.pointsCount());
	fitGauss.setData(data);
	printf("data size is :%li\n",data.pointsCount());
	
	fitGauss.setSaveParameters(0,0,0);
	fitGauss.setDumpEveryMCMCstate(true);
	fitGauss.setCoolingRate(1);
	fitGauss.setBurnInLength(1000);
	fitGauss.setCovarianceMatrixDiagonal(true);
	fitGauss.setMaximalRejectionsCount(50);
	fitGauss.setMaximalChainLength(20000);
	//
	// RUN THE MCMC CHAIN
	//
	fitGauss.startChain();
	fitGauss.saveResults();
	
	return 0;
}





