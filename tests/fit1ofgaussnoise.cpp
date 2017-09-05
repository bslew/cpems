/*
 * A program that should fit 1/f + gauss noise model to the OCRA-f data using MCMC
 * 
 * CREATED ON:
 * Feb 3, 2011, 8:55:26 AM
 * 
 */
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <tclap/CmdLine.h>
#include "cpedsMCMCpowerFitOCRA.h"
#ifdef GOMPI
#include "mpi.h"
#endif

/***************************************************************************************/
#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
#endif
//using namespace MATPACK; // use all classes and functions in the MATPACK namespace
using namespace TCLAP;
/***************************************************************************************/

/******************************************************************************************/
/* GLOBAL VARIABLES                                                                                        */
/******************************************************************************************/
void parseOptions(int argc, char** argv);

bool _testName, _cutBelowFlag, _cutAboveFlag;
double _cutBelow, _cutAbove,_bingm, _coolingRate,_BIfrac,_ABIfrac;
double _log10Amin,_log10Amax,_k0min,_k0max,_nMin,_nMax, _fkneeMin,_fkneeMax;
long _idx,_binw,_binst, _burnInLength, _dirStIdx,_NcovSim, _forcedCoolingsCount;
string _outfile, _input_file, _maskfile, _mcmcDirsPrefix;
vector<string> _msg;
/***************************************************************************************/








/***************************************************************************************/
int main(int argc, char** argv) {
	int ierr,node_num=1;
	int world_id=0;
	cpedsMsgs msgs("cpedsMCMC");
	parseOptions(argc, argv);
	long idx=0;
	
#ifdef GOMPI
	ierr = MPI_Init ( &argc, &argv ); // initiate MPI
	ierr = MPI_Comm_size ( MPI_COMM_WORLD, &node_num ); // Get the number of processes
	ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &world_id ); // Get the individual process ID
	idx=world_id;
#endif	

	
	mscsFunction oof,gauss,data,maskRanges,binnedData;

	
	/******************************************************************************************/
	/* DATA                                                                                       */
	/******************************************************************************************/
	if (_testName) {
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
		_maskfile="";
	}
	else {
		//
		// LOAD THE INPUT POWER SPECTRUM DATA FROM FILE
		//
		data.load(_input_file);
		data.deletePoint(0); // remove the zero'th k mode
		
		
		//
		// REMOVE DATA OUTSIDE OF REQUESTED RANGE
		//
		if (_cutBelowFlag or _cutAboveFlag) {
			if (_cutBelowFlag==false) _cutBelow=data.getMinArg();
			if (_cutAboveFlag==false) _cutAbove=data.getMaxArg();
			data.deleteOutsideOf(_cutBelow,_cutAbove);
		}

		binnedData=data;
	
	}


	/******************************************************************************************/
	/* MASK                                                                                   */
	/******************************************************************************************/
	
	// bin data
	cpedsList<long> bins;
	cpedsList<double> dataMask;
	
	
	//
	// load the masked data ranges here. It should be the mask removing the spikes from the signal.
	//
	if (_maskfile!="") {
		msgs.say("processing mask",High);
		maskRanges.load(_maskfile);

		// generate mask for the data set

		long n=data.pointsCount();
		long m=maskRanges.pointsCount();
		double min,max,x;
		dataMask.makeLength(n);
		for (long j = 0; j < m; j++) {
			if (world_id==0) { msgs.say("processing mask range: "+msgs.toStr(j)+" of "+msgs.toStr(m),Medium); }
			min=maskRanges.getX(j);
			max=maskRanges.getY(j);
			for (long i = 0; i < n; i++) {
				x=data.getx(i);
				if (x>min and x<max) dataMask[i]=0;
				else dataMask[i]=1;
			}
		}
	}
		
	//
	// bin data 
	//
	binnedData.binFunctionLin2geo(_binst,_binw,_bingm,bins,dataMask);
	
	// normalize data to one for large k-modes
	//	data/=double(data.f(data.pointsCount()-1)); // no need to normalize in case of the fiducial data
	
	// convert data into a log-log space
//	data.logX();
	binnedData.logY();
	// 
	
	// save the fiducial data
//	if (world_id==0) {
//		binnedData.save("noisePower-data.txt");
//	}	
	
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
	
	
	/******************************************************************************************/
	/* PARAMETER SPACE                                                                        */
	/******************************************************************************************/
	
	mscsFunction A("log A"); 
	mscsFunction k0("f_0"); 
	mscsFunction n("n"); 
	mscsFunction fknee("f_{knee}"); 
	

	if (_testName) {
		// make a flat prior on parameter k0
		k0.mkConst(0,3,0.01,1.0);
		k0.normalize();
		// make a flat prior on parameter n
		n.mkConst(-3.0,0.0,0.01,1.0);
		n.normalize();
		// make a flat prior on parameter fknee
		fknee.mkConst(0,3,0.01,1.0);
		fknee.normalize();
	}
	else {
		// make a flat prior on parameter A
//		A.mkConst(0,1e5,10,1.0); // constraint from fknee >0.1 Hz // run 0
//		A.normalize();
		// make a flat prior on parameter log10 A
		A.mkConst(_log10Amin,_log10Amax,.01,1.0); 
		A.normalize();
		
		// make a flat prior on parameter k0
//		k0.mkConst(0,data.getMaxArg(),0.01,1.0); // run 0
		k0.mkConst(_k0min,_k0max,0.01,1.0); // run 1
		k0.normalize();

		// make a flat prior on parameter n
		n.mkConst(_nMin,_nMax,0.01,1.0);
		n.normalize();

		// make a flat prior on parameter fknee
//		fknee.mkConst(0,data.getMaxArg(),0.01,1.0); // run 0
		fknee.mkConst(_fkneeMin,_fkneeMax,0.01,1.0); // run 1
		fknee.normalize();		
	}
	
	
	//
	// define the covariance matrix for the input data within the defined model
	//
	
	// number of simulations
	
	
	/******************************************************************************************/
	/* MCMC SETUP                                                                                       */
	/******************************************************************************************/
	
	
	// set the dimensionality of the problem
	if (_testName) {
		OCRAPowerMCMC fitSpectra(3,idx,1000);
		if (_mcmcDirsPrefix!="") _mcmcDirsPrefix+="-";
		fitSpectra.msgs->setSender("MCMC."+msgs.toStr(world_id));
		if (node_num>1) {
			fitSpectra.setOutputDir(_mcmcDirsPrefix+msgs.toStr("mcmc.")+msgs.toStr(_idx));
		}
		else 
			fitSpectra.setOutputDir(_mcmcDirsPrefix+msgs.toStr("mcmc.")+msgs.toStr(_idx));
//			fitSpectra.setOutputDir(msgs.toStr("mcmc.")+msgs.toStr(argv[1]));
//			fitSpectra.setOutputDir("mcmc");
		
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
		
		
		//	double (*chisqPtr)(const mscsFunction&, const MClink& );
		//	chisqPtr=&(fitSpectra.chisqNoCovNoLog);
		fitSpectra.setSaveParameters(0,0,0);
		fitSpectra.setDumpEveryMCMCstate(false);
		fitSpectra.setCoolingRate(_coolingRate);
		fitSpectra.setBurnInLength(_burnInLength);
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
	}
	else {
		OCRAPowerMCMC fitSpectra(4,idx,1000);
		fitSpectra.msgs->setSender("MCMC."+msgs.toStr(_idx));
		fitSpectra.setDirStIdx(_dirStIdx);
		if (node_num>1) {
			fitSpectra.setOutputDir(_mcmcDirsPrefix+msgs.toStr("mcmc.")+msgs.toStr(_idx));
		}
		else 
			fitSpectra.setOutputDir(_mcmcDirsPrefix+msgs.toStr("mcmc.")+msgs.toStr(_idx));
//			fitSpectra.setOutputDir("mcmc");
		//
		// add parameters
		//
		fitSpectra.setInitialStepSize(_BIfrac);
		fitSpectra.setInitialWalkStepSize(_ABIfrac);
		fitSpectra.addParameter(k0);
		fitSpectra.addParameter(n);
		fitSpectra.addParameter(fknee);
		fitSpectra.addParameter(A);
		
		//
		// set the observational data
		//
		fitSpectra.setData(binnedData);
		
		//
		// set chisq related variables
		//
		fitSpectra._simData=data;
		fitSpectra._binst=_binst;
		fitSpectra._binw=_binw;
		fitSpectra._bingm=_bingm;
		fitSpectra._NcovSim=_NcovSim;
		fitSpectra.setCovRNG();
#ifdef COVFAST
		fitSpectra.prepareCovCalc(bins);
#endif
//		if (_NcovSim==0 and fitSpectra.
		
		//
		// WRITE IMPLEMENTATION NOTES FOR THE RUN
		//
		fitSpectra.addImplementationNote("");
		fitSpectra.addImplementationNote(string("DATA OPERATION OPTIONS: "));
		fitSpectra.addImplementationNote(string("mask file: ")+_maskfile);
		fitSpectra.addImplementationNote(string("cut data below flag: ")+msgs.toStr(_cutBelowFlag));
		fitSpectra.addImplementationNote(string("cut data below frequency: ")+msgs.toStr(_cutBelow));
		fitSpectra.addImplementationNote(string("cut data above flag: ")+msgs.toStr(_cutAboveFlag));
		fitSpectra.addImplementationNote(string("cut data above frequency: ")+msgs.toStr(_cutAbove));

		fitSpectra.addImplementationNote("");
		fitSpectra.addImplementationNote(string("PROGRAM OPTIONS"));
		fitSpectra.addImplementationNote(string("program name: fit1ofgaussnoise.cpp"));
		fitSpectra.addImplementationNote(string("argv[0] is: ")+argv[0]);
		fitSpectra.addImplementationNote(string("Number of simulations for the covariance matrix estimation: ")+msgs.toStr(_NcovSim));
#ifdef COVFAST
		fitSpectra.addImplementationNote(string("The program was compiled with COVFAST compiler option."));
#ifdef BINNED 
		fitSpectra.addImplementationNote(string("The program was also compiled with BINNED compiler option. The noise realizations for covariance matrix calculation will be done on binned grid of points and the covariance matrix will be rescalled by binsSize factors"));
#endif
#endif
		fitSpectra.addImplementationNote("");
		fitSpectra.addImplementationNote(string("BINNING OPTIONS: "));
		fitSpectra.addImplementationNote(string("binning starts at: ")+msgs.toStr(_binst));
		fitSpectra.addImplementationNote(string("initial bin width: ")+msgs.toStr(_binw));		
		fitSpectra.addImplementationNote(string("bin width multiplier: ")+msgs.toStr(_bingm));
		if (_msg.size()>0) {
			fitSpectra.addImplementationNote("");
			fitSpectra.addImplementationNote(string("RUN COMMENT: "));
			for (long i = 0; i < _msg.size(); i++) {
				fitSpectra.addImplementationNote(_msg[i]);				
			}
		}
		//
		// set MCMC parameters
		//
		fitSpectra.setSaveParameters(0,0,0);
		fitSpectra.setDumpEveryMCMCstate(false);
//		fitSpectra.setCoolingRate(0.1); // very slow cooling rate
//		fitSpectra.setCoolingRate(0.5); // medium cooling rate 
//		fitSpectra.setCoolingRate(1.0); // normal cooling rate
		fitSpectra.setCoolingRate(_coolingRate); // normal cooling rate
		fitSpectra.setBurnInLength(_burnInLength);
		fitSpectra.setCovarianceMatrixDiagonal(true);
		fitSpectra.setDataMask(dataMask);
		fitSpectra.setMaximalChainLength(25000);
		fitSpectra.setUpHillClimbing(true);
		double p=0.05;
//		fitSpectra.setAcceptWorsePvalues(0.3,0.05);
//		fitSpectra.setMaximalRejectionsCount(1.0/p);
		fitSpectra.setMaximalRejectionsCount(100);
		fitSpectra.setAcceptWorseUniformScheme(false);
		fitSpectra.setAcceptWorseBoltzmannSchemeFactors(1000,10);
		fitSpectra.setBestFitCovSimNum(5000);
		fitSpectra.setMaximalNumberOfForcedCoolings(_forcedCoolingsCount);
		//
		// RUN THE MCMC CHAIN
		//
		fitSpectra.startChain();
		fitSpectra.saveResults();
	}
		
#ifdef GOMPI
		ierr = MPI_Finalize ( ); // Shut down MPI.
#endif
		
		return 0;
}




/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/
/***************************************************************************************/

void parseOptions(int argc, char** argv) {
	long i;
	string::size_type j;
	
	try {
		
		CmdLine cmd("test-cpedsMCMC: a test program for reconstructing power spectra parameters\n", ' ', "" );
		
		//
		// Define arguments
		//
		
		
		UnlabeledValueArg<string> input_file("input_file","file XY function file with raw power spectra","nofile","string");	cmd.add( input_file );
//		
//		ValueArg<string> output("o","out","outfile ",false,"out","string"); cmd.add(output);
		ValueArg<string> maskfile("m","mask","mask file",false,"","string"); cmd.add(maskfile);
		ValueArg<string> mcmcDirsPrefix("","mcmcDirsPrefix","prefix for the name of the output mcmc directories",false,"","string"); cmd.add(mcmcDirsPrefix);

//		vector<string> runNotes;
//		runNotes.push_back(12);
//		runNotes.push_back(21);
//		ValueArg<long> direction("d","direction","shift direction (default:12)",false,12,allowedDirs); cmd.add(direction);
		
		MultiArg<string> msgNote("","msg","add a description note for this run; it will be saved in setup.mc files.",false,"string"); cmd.add(msgNote);
//		
		ValueArg<long> idx("i","idx","unique index of the mcmc run (default: 1)",false,1,"long"); cmd.add(idx);
//		ValueArg<long> shift("s","shift","Number of points by which to shift the window (default: 10)",false,10,"long"); cmd.add(shift);
		ValueArg<long> dirStIdx("","dirStIdx","starting index for the first chain used to name the directory for the run. Other mcmc chains will have indexes according to their world_id+dirStIdx (default: 1)",false,1,"long"); cmd.add(dirStIdx);

		ValueArg<long> binst("","bst","data bin starting from (default: 1)",false,1,"long"); cmd.add(binst);
		ValueArg<long> binw("","bw","data first bin width(default: 10)",false,10,"long"); cmd.add(binw);
		ValueArg<double> bingm("","gm","data bin geometrical sequence multiplier (default: 1.05)",false,1.05,"double"); cmd.add(bingm);
		ValueArg<long> NcovSim("","Ncov","number of simulations used to calculate covariance matrix for a model (default: 0)",false,0,"long"); cmd.add(NcovSim);
		
		ValueArg<double> log10Amin("","log10Amin","lower bound on log_10 A (default: 0)",false,0,"double"); cmd.add(log10Amin);
		ValueArg<double> log10Amax("","log10Amax","upper bound on log_10 A (default: 5)",false,5,"double"); cmd.add(log10Amax);

		ValueArg<double> k0min("","k0min","lower bound on k0 (default: 0)",false,0,"double"); cmd.add(k0min);
		ValueArg<double> k0max("","k0max","upper bound on k0 (default: 200)",false,200,"double"); cmd.add(k0max);
		ValueArg<double> nMin("","nMin","lower bound on n (default: -3)",false,-3,"double"); cmd.add(nMin);
		ValueArg<double> nMax("","nMax","upper bound on n (default: 0)",false,0,"double"); cmd.add(nMax);
		ValueArg<double> fkneeMin("","fkneeMin","lower bound on 1/f knee (default: 0)",false,0,"double"); cmd.add(fkneeMin);
		ValueArg<double> fkneeMax("","fkneeMax","upper bound on 1/f knee (default: 200)",false,200,"double"); cmd.add(fkneeMax);

		ValueArg<double> coolingRate("","coolingRate","cooling rate; 0.1 - very slow, 0.5 - slow, 1 - normal and so on (default: 1)",false,1,"double"); cmd.add(coolingRate);
		ValueArg<long> burnInLength("","burnInLength","minimal length of the MC chain before the system starts to cool down; (default: 1000)",false,1000,"long"); cmd.add(burnInLength);
		
		ValueArg<double> BIfrac("","BIfrac","fraction of the size of the confidence range corresponding to the confidence level equal to expectation value used to calibrate step generating PDF and RNG during burn-in period (default 0.5)",false,0.5,"double"); cmd.add(BIfrac);
		ValueArg<double> ABIfrac("","ABIfrac","fraction of the size of the confidence range corresponding to the confidence level equal to expectation value used to calibrate step generating PDF and RNG after burn-in period  (default 0.1)",false,0.1,"double"); cmd.add(ABIfrac);
		//		
		ValueArg<double> cutbelow("","cutBelow","low frequencies cut-off threshold (default 0.1)",false,0.1,"double"); cmd.add(cutbelow);
		ValueArg<double> cutabove("","cutAbove","high frequencies cut-off threshold (default 137)",false,0.1,"double"); cmd.add(cutabove);

		ValueArg<long> forcedCoolingsCount("","coolingTimes","maximal number of forced coolings; (default: 10)",false,10,"long"); cmd.add(forcedCoolingsCount);
		
		
		//		ValueArg<double> startFrom("f","startFrom","In the asIds=true mode the flag indicates the point index to start from (from the left for direction 12 and from the right for direction 21)."
//				"In the asIdx=false mode the flag indicates the value of the argument closest to the given number. The window will start/end on this number and will be  shifted to the right/left depending on the value of the direction flag 12/21. "
//				"If not set, the full range will be considered"
//				"(default 0)",false,0,"double"); cmd.add(startFrom);
//		ValueArg<double> endAt("t","endAt","In the asIds=true mode the flag indicates the point index to end at (from the left for direction 12 and from the right for direction 21)."
//				"In the asIdx=false mode the flag indicates the value of the argument closest to the given number. The window will start/end on this number and will be  shifted to the right/left depending on the value of the direction flag 12/21. "
//						"(default 0)",false,0,"double"); cmd.add(endAt);
//		SwitchArg startFromAsIdx("","asIdx", "whether or not to treat the startFrom parameter as the point index or the value of the X argument of the function", false);	cmd.add(startFromAsIdx);
				
		SwitchArg test("","test", "starts a demo test of the program.", false);	cmd.add(test);
				
//		SwitchArg loop("","loop", "loops the spikes removal untill no points are removed for the current window and shift parameters."
//				"In case of the looped run, the output file with spikes will contain the removed points only from the last pass of the loop", false);	cmd.add(loop);
//
//		SwitchArg removeInLoop("","removeInLoop", "Remove window of width w also at every pass over the loop (relevant only for loop calculations).", false);	cmd.add(removeInLoop);
//
//				
//		ValueArg<double> width("w","width","The half width of the window centered on the detected spike position within which the data are to be considered as contaminated and treated as spikes. "
//				"For 0 (default 0) value only the detected spike will be in the output file.",false,0,"double"); cmd.add(width);
				
		//
		// Parse the command line.
		//
		cmd.parse(argc,argv);
		
		//
		// Set variables
		//
		_input_file = input_file.getValue(); /*if ( _input_file == "nofile" || _input_file == "Ylm" ) { _not_from_file = true; } else { _not_from_file = false; }*/
		_mcmcDirsPrefix=mcmcDirsPrefix.getValue();
		_msg=msgNote.getValue();


		_log10Amin=log10Amin.getValue();
		_log10Amax=log10Amax.getValue();
		_k0min=k0min.getValue();
		_k0max=k0max.getValue();
		_nMin=nMin.getValue();
		_nMax=nMax.getValue();
		_fkneeMin=fkneeMin.getValue();
		_fkneeMax=fkneeMax.getValue();
		
		_coolingRate=coolingRate.getValue();
		_burnInLength=burnInLength.getValue();
		_BIfrac=BIfrac.getValue();
		_ABIfrac=ABIfrac.getValue();
		
		//		_outfile = output.getValue();
		_idx= idx.getValue();
		_dirStIdx= dirStIdx.getValue();
		_forcedCoolingsCount=forcedCoolingsCount.getValue();


		// PARAMETERS STORED AS NOTES IN SETUP
		_maskfile= maskfile.getValue();
		_cutBelow= cutbelow.getValue(); if (cutbelow.isSet()) _cutBelowFlag=true; else _cutBelowFlag=false;
		_cutAbove= cutabove.getValue(); if (cutabove.isSet()) _cutAboveFlag=true; else _cutAboveFlag=false;
		_binst=binst.getValue();
		_binw=binw.getValue();
		_bingm=bingm.getValue();
		_NcovSim=NcovSim.getValue();		
		
		
		//		_shiftSize = shift.getValue();
//		_direction = direction.getValue();
//		// _loffset = loffset.getValue();
//		_nsigma = nsigma.getValue();
//		_startFrom = startFrom.getValue();
//		_endAt = endAt.getValue(); if (endAt.isSet()) _endAtSet=true; else _endAtSet=false; 
//		_startFromAsIdx = startFromAsIdx.getValue();
		_testName = test.getValue();
//		_width=width.getValue();
		
//		_loop= loop.getValue();
//		_removeInLoop= removeInLoop.getValue();
				
	} catch ( ArgException& e )
	{ cout << "ERROR: " << e.error() << " " << e.argId() << endl; }
}



