/*!
 * \file cpedsMCMC.cpp
 *
 *  Created on: Nov 19, 2010
 *      Author: blew
 */

#include <math.h>
#include <assert.h>
#include "cpedsMCMC.h"
#include "cpeds-list.h"
#include "cpeds-point_set.h"
#include "subdomain.h"
#include "MscsPDF2D.h"
#include <iostream>
#include <assert.h>

/***************************************************************************************/
//cpedsMCMC::cpedsMCMC(int runIdx) {
//	initialize(runIdx);
//}
/***************************************************************************************/
//cpedsMCMC::cpedsMCMC(long Npar,int runIdx, long runOffset, long runSeed, cpeds_VerbosityLevel verb) : current()(Npar) {
cpedsMCMC::cpedsMCMC(long Npar,int runIdx, long runOffset, long runSeed, cpeds_VerbosityLevel verb)  {
	current().setNpar(Npar);
	msgs=NULL;
	initialize(runIdx,runOffset,Npar, runSeed);
	msgs->setVerbosity(verb);
}
/***************************************************************************************/
cpedsMCMC::cpedsMCMC(cpedsMCMC &parent) {
	msgs = new cpedsMsgs("cpedsMCMC");
	this->operator =(parent);			

}
/***************************************************************************************/
cpedsMCMC::~cpedsMCMC() {
	delete msgs;
//	if (_covDiag!=NULL) delete [] _covDiag;
//	delete [] _walk.walk.currentDirection;
}

/***************************************************************************************/
void cpedsMCMC::initialize(int runIdx, long runOffset, long Npar, long runSeed) {
	//
	// run/debug/save control stuff
	//
	_IOcontrol.mcmcRunIdx=runIdx;
	_IOcontrol.runOffset=runOffset*runIdx;
	_IOcontrol.dirStartIndex=0;
	if (_IOcontrol.mcmcRunIdx==-1) {
		_IOcontrol.runFilesDir="mcmc_run";
	}
	else {
		_IOcontrol.runFilesDir="mcmc_run."+msgs->toStr(_IOcontrol.mcmcRunIdx+_IOcontrol.dirStartIndex);
	}
	_IOcontrol.partialFilesDir=getOutputDir()+"/"+getPartialFilesDir();
	_IOcontrol.saveCoolingPDFEveryNthState=0; // do not save anything
	_IOcontrol.saveStepPDFsEveryNthState=0;
	_IOcontrol.saveModelEveryNthState=0;
	_IOcontrol.chisqSignature="abstract";
	
	_IOcontrol.saveCoolingPDF_lastSavedNumber=0;
	_IOcontrol.saveStepPDFs_lastSavedNumber=0;
	_IOcontrol.saveModel_lastSavedNumber=0;
	_IOcontrol.dumpEveryMCMCstate=false;
	_IOcontrol.storeBestFitNow=false;
//	_covDiag=NULL;
	_chisqData.covDiagonal=false;
	
	_IOcontrol.PDF1d_pointsCount=250;
	_IOcontrol.PDF2d_pointsCount=50;
	_IOcontrol.calculate2Dpdf=false;
	//
	// walk parameters for step PDFs and RNGs
	// 
	_cooling.stepGeneratingPDFPoints=1000;
	_walk.runRNGseed=runSeed; // all other rngs will use this seed to initialize themselves (with different offsets of course)
	_walk.rnSgn.seed(_walk.runRNGseed); 
	_walk.rnSgn.setRNsType(_walk.rnSgn.uniform);
	_walk.rnSgn.setMinMax(-1,1);
//	_walk.rnSgn.seedOffset(100+_IOcontrol.runOffset);
	_walk.rnSgn.seedOffset(cpeds_get_devrandom_seed());
	_walk.Nparam=Npar;
	_walk.maximalChainLength=200000;
	_walk.statesTot=0;
	_walk.initialPDFstepCRfraction=0.5;
	_walk.initialPDFstepCRfractionAfterBurnIn=0.05;
//	current()Direction=new double[Npar];
	_walk.currentDirection.makeLength(Npar);
	cpedsList<double> init_step;
	_walk.steps.clear();
	for (unsigned long i = 0; i < Npar; i++) { _walk.steps.push_back(init_step); }
	_walk.keepDirectionNow=false;
	_walk.keepDirectionStepGenerationFailuresNumber=50;
	_walk.rejectedLinkNumber=0;
	_walk.uphillClimbing=true;
	_walk.uphillGradient=false;
	_walk.userOutputFrequency=1000;
	_walk.likelihoodRejectionThreshold=1.0e-10;
	_walk.convergenceTestMinimalLength=300; 
	_walk.burnOutLength=_walk.convergenceTestMinimalLength; 
	
	//
	// convergence variables
	//
	_convergence.relChisqChangeThres=0.0001;

	//
	// set the cooling PDF and corresponding RNG for "down-hill" walk
	//
	_cooling.maximalRejectionsCount=30.0*sqrt(Npar);
	_cooling.maximalNumberOfForcedCoolings=10;
	_cooling.initialTemperature=1000;
	_cooling.temperature=getInitialTemperature();
	_cooling.forcedCoolingTemperatureDecrement=-1;
	//	_temperatureHistory.append(_temperature);
	_cooling.finalTemperature=getInitialTemperature()*_convergence.relChisqChangeThres;
	_cooling.initialMaximalEnergy=10; // i.e. 10 kT - there should be no need to change this
	setCoolingScheme(cpedsMCMC::coolingScheme_chisqPropLinLog);
	
	
	
	//
	// convergence variables again
	//
	
	_convergence.last_chisq=-1;
	_convergence.minimalLengthToTestConvergence=getBurnInLength();
	_convergence.nextConvergenceCheck=-1;
	_convergence.statesToConvergenceCheck=_cooling.maximalRejectionsCount/2;

	//
	// cooling variables
	//
	_cooling.acceptWorseStateFromUniformDistr=true; // this means we don't use boltzmann distribution
	mscsFunction f("coolingPDF",Zero);
	f.mkBoltzmann(0,_cooling.initialMaximalEnergy*CPEDS_kB*_cooling.initialTemperature,_cooling.initialMaximalEnergy*CPEDS_kB*_cooling.initialTemperature/_cooling.stepGeneratingPDFPoints,_cooling.initialTemperature);
	if (_cooling.acceptWorseStateFromUniformDistr==false) { // i.e. actually use the boltzmann distribution
		// scale the distribution so it make sense to compare the arguments with delta chisq
		f.scaleX(1.0/(CPEDS_kB*getTemperature()));
	}
	f.normalize();
	_cooling.coolingDistr=f;
	_cooling.coolingDistr.setVerbosityLevel(Zero);
	_cooling.coolingDistrInvCDF=_cooling.coolingDistr;
	_cooling.coolingDistrInvCDF.setVerbosityLevel(Zero);
	_cooling.coolingDistrInvCDF.mkCDF(true);
	_cooling.coolingDistrInvCDF.invert();
	_cooling.coolingRNG.setRNsType(_cooling.coolingRNG.from_array_invCDF);
	_cooling.coolingRNG.setPDF(_cooling.coolingDistr.pointsCount(),_cooling.coolingDistr.extractArguments(),_cooling.coolingDistr.extractValues());
//	_cooling.coolingRNGseedOffset=101+_IOcontrol.runOffset;
	_cooling.coolingRNGseedOffset=cpeds_get_devrandom_seed();
	_cooling.coolingRNG.seedOffset(_cooling.coolingRNGseedOffset);
	_cooling.coolingRNG.seed(_walk.runRNGseed);
	setAcceptWorsePvalues(0.3,0.05);
	setAcceptWorseBoltzmannSchemeFactors(1,1);
	_cooling.coolingRate=0.5;
	
	_chisqData.data2dCol=-1;
	_chisqData.bestFitNcovSimNum=5000;
	_chisqData.rnCov.setRNsType("gaussian_circle");
//	_chisqData.rnCov.seedOffset(102+_IOcontrol.runOffset);
	_chisqData.rnCov.seedOffset(cpeds_get_devrandom_seed());
	_chisqData.rnCov.seed(_walk.runRNGseed);
	_chisqData.bestFitData.confidenceLevels.append(0.68);
	_chisqData.bestFitData.confidenceLevels.append(0.95);
	_chisqData.bestFitData.confidenceLevels.append(0.99);
	
	
	// messanger
	if (msgs==NULL) msgs = new cpedsMsgs("cpedsMCMC");
	
}
/***************************************************************************************/
double cpedsMCMC::chisq(const mscsFunction& data, const MClink &th, matrix<double>& cov) {
	msgs->criticalError("you are calling the wrong chisq function. This is the abstract method that needs to be reimplemented in the derived class",High);
	return 0;
}
/***************************************************************************************/
double cpedsMCMC::chisq(const MClink &th) {
	//	msgs->say("calling short native chisq function",High);
	matrix<double> m;
//	matrix<double> m(_chisqData.data.pointsCount(),_chisqData.data.pointsCount());
//	for (long i = 0; i < _chisqData.data.pointsCount(); i++) {
//		for (long j = 0; j < _chisqData.data.pointsCount(); j++) {
//			if (i==j) m(i,i)=1; else m(i,j)=0;
//		}
//	}
	msgs->criticalError("you are calling the wrong chisq function. This is the abstract method that needs to be reimplemented in the derived class",High);
	
//	return chisq(_chisqData.data,th,m);
	return -1;
}

/***************************************************************************************/
void cpedsMCMC::addParameter(const mscsFunction& p, string parameter_full_name) { 
//	_parameterSpace.append(p); 
//	_parameterSpace.last().checkRanges();
	_parameterSpace.addParameter(p,p.getName(),parameter_full_name);

	_walk.Nparam=_parameterSpace.size();
	// prepare RNGs for drawing starting point for the chain
	_walk.rns.append(cpedsRNG());
	_walk.rns.last().setRNsType(_walk.rns.last().from_array_invCDF);
	_walk.rns.last().setPDF(p.pointsCount(),p.extractArguments(),p.extractValues());
//	_walk.rns.last().seedOffset(_walk.rns.size()+1+_IOcontrol.runOffset);
	_walk.rns.last().seedOffset(cpeds_get_devrandom_seed());
	_walk.rns.last().seed(_walk.runRNGseed);
	
	// define the grid size
	_walk.delta0.append(p.getX(1)-p.getX(0));
	
	//
	// define step size generating distribution
	//
	mscsFunction f(p.getName()+"PDF",Zero);
	f.mkBoltzmann(0,_cooling.initialMaximalEnergy*CPEDS_kB*_cooling.initialTemperature,_cooling.initialMaximalEnergy*CPEDS_kB*_cooling.initialTemperature/_cooling.stepGeneratingPDFPoints,_cooling.initialTemperature);
	f.scaleX(_walk.initialPDFstepCRfraction*(_parameterSpace.last().getMaxArg()-_parameterSpace.last().getMinArg())/f.getMaxArg()); // make the maximal step size equal to half of parameter space prior - this is because the step can be positive or negative
	f.normalize(); // normalize PDF
	_walk.stepGeneratingDistr.append(f);
	_walk.stepGeneratingDistr.last().setVerbosityLevel(Zero);
	
	
	// define step size generating RNG
	_walk.stepGeneratingRNG.append(cpedsRNG());
	_walk.stepGeneratingRNG.last().setRNsType(_walk.stepGeneratingRNG.last().from_array_invCDF);
	_walk.stepGeneratingRNG.last().setPDF(f.pointsCount(),f.extractArguments(),f.extractValues());
//	_walk.stepGeneratingRNG.last().seedOffset(-_walk.stepGeneratingRNG.size()-1+_IOcontrol.runOffset);
	_walk.stepGeneratingRNG.last().seedOffset(cpeds_get_devrandom_seed());
	_walk.stepGeneratingRNG.last().seed(_walk.runRNGseed);
	
	
}
/***************************************************************************************/
void cpedsMCMC::addParameter(string param_name, double from, double to, long Npts, string param_full_name) {
	mscsFunction param(param_name);
	double dx=(to-from)/Npts;
	param.mkConst(from,to,dx,1.0); 
	param.normalize();
	addParameter(param,param_full_name);
}
/* ******************************************************************************************** */
void cpedsMCMC::addParameter(string param_name, double from, double to, double delta, string param_full_name) {
	mscsFunction param(param_name);
	double dx=delta;
	param.mkConst(from,to,dx,1.0); 
	param.normalize();
	addParameter(param,param_full_name);	
}
/***************************************************************************************/
void cpedsMCMC::walkNstates(long Nstates, MClink startingLink, int how) {
	bool walk=true;
	long burnOutStates=0;
	if (Nstates==0) walk=false;
	current()=startingLink;
	
	while (walk) {
		_walk.statesTot++;

		
		nextCandidate()=getNextPoint(current());
		nextCandidate().setIdx(_walk.statesTot);
		nextCandidate().setChisq(chisq(nextCandidate()));

		nextCandidate().setAccepted(true);
		_MCaccepted.append(nextCandidate());
		
		if (nextCandidate().chisq()<_bestFit.chisq()) {
			addNewBestFitLink(nextCandidate());
		}

		switch (how) {
			case 0:
				current()=nextCandidate(); 				
				break;
			case 1:
				current()=bestFitLink();
				break;
			default:
//				current()=startingLink();
				break;
		}
		//
		// saving section
		//
//		dumpAll();

		// control burn-out length
		burnOutStates++;
		if (burnOutStates>=getBurnOutLength()) walk=false;
		
		//
		// control the maximal chain length
		//		
		if (length()==_walk.maximalChainLength) walk=false;

		
		//
		// output
		//
		if (length() % getWalkInfoOutputFrequency() == 0 or walk==false or msgs->getVerbosity()>=High)
			msgs->say("B-out: MC len.: %.0f, next.X2: %lE, bf.X2: %lE, left: %.0f ",
					double(length()),
					nextCandidate().chisq(),
					_bestFit.chisq(),
					double(getBurnOutLength()-burnOutStates),
					Medium);

	}
}
/* ******************************************************************************************** */
void cpedsMCMC::startChain(bool followPriors) {
	mkOutputDir();

//	MClink next(dims());
	double deltaX,deltaX2perDOF;
	_cooling.forcedCoolingTemperatureDecrement=-1;
	double acceptWorseThreshold;
	double acceptWorseBoltzmannFactor;
	double acceptWorseRN;
	
	_walk.rejectionsCount=0;
	
	
	//
	// define the starting point of the walk in the parameter space 
	//
	current().setNpar(dims());
	_walk.statesTot=0;
	
	// randomly choose the starting location in the parameter space
	msgs->say("Generating random starting point in the prameter space",Medium);
	current().set(dims(),getStartingPoint());
	current().setAccepted(true);
	current().setIdx(0);
	_startingLink=current();
	// calculate the chisq there
	current().setChisq(chisq(current())); // get first chisq and assign chisq signature
	_convergence.last_chisq=current().chisq();
//	printInfo(); // print and store setup here to get the chisq assigned correctly before storing
	
	_MCaccepted.append(current());
	appendLinkForConvergenceTest(current());
	_bestFit=current();
	_cooling.temperatureHistory.append(getTemperature());
	
	bool walk=true;
	bool walkDone=false;
	bool burnOutDone=false;
	long burnOutStates=0;
	long statesTotBeforeBurnOut=0;
	
	_startingLink.save(getRunDependentFileName(getPartialFilesDirFull()+"/startingLink")); 

	msgs->say("STARTING THE MC CHAIN",High);
	while (walk) {
		//
		// generate next step
		//
		_walk.statesTot++;
		if (length()==getBurnInLength()) { 
			deltaX=current().chisq()-nextCandidate().chisq();
			nextCandidate()=_bestFit;  // go back to the best fit link in the chain after the burn-in period and continue from there
			nextCandidate().setIdx(_walk.statesTot);
			msgs->say("Reached end of burin-in period",Medium);
			
			if (_cooling.acceptWorseAfterBurnInBoltzmannFact==-1) {
				_cooling.acceptWorseAfterBurnInBoltzmannFact=_bestFit.chisq()/dataSize();
			}
			
			if (_walk.initialPDFstepCRfractionAfterBurnIn>0) { // regenerate the step size PDFs for all parameters according to the new initial step size to be used after burn-in period
				generateInitialStepSizePDFsAfterBurnIn();
			}
			//			if (deltaX > 0) { // cool down at the same time to reduce chance of going away from the best fit solution
			//				coolDownLogLin(deltaX);
			//			}
		}
		else {
			if (getInitialStepSize()==0 and length()<getBurnInLength()) { // we are in the burn-in phase
				MClink randomLink;
				randomLink.set(dims(),getStartingPoint());
				randomLink.setAccepted(true);
				randomLink.setIdx(_walk.statesTot);

				nextCandidate()=randomLink;
				msgs->say("Current temperature is: "+msgs->toStr(getTemperature())+", final temperature is:"+msgs->toStr(getFinalTemperature()),Zero);
				nextCandidate().setChisq(chisq(nextCandidate()));
			}
			else { // we are no longer in the burn-in phase
				if (getUphillGradient()) {
					nextCandidate()=getNextPointUphillGradient(current());
				}
				else {
					nextCandidate()=getNextPoint(current());
					msgs->say("Current temperature is: "+msgs->toStr(getTemperature())+", final temperature is:"+msgs->toStr(getFinalTemperature()),Zero);
				}
				nextCandidate().setIdx(_walk.statesTot);
				nextCandidate().setChisq(chisq(nextCandidate()));
			}
		}
		
		//
		// decide what to do with the new state
		//
		
		// define the current accept worse threshold
		if (length()<=getBurnInLength()) { acceptWorseThreshold=_cooling.acceptWorseThreshold; acceptWorseBoltzmannFactor=_cooling.acceptWorseBoltzmannFact; }
		else  { acceptWorseThreshold=_cooling.acceptWorseThresholdAfterBurnIn; acceptWorseBoltzmannFactor=_cooling.acceptWorseAfterBurnInBoltzmannFact; }
		
		
//		if (nextCandidate().chisq()<bestFitLink().chisq()) {
		if (nextCandidate().chisq()<current().chisq()) {
//		if ((current().chisq()-next.chisq())/_bestFit.chisq()>_convergence.relChisqChangeThres) {
			msgs->say("the new chisq is < than the current one: ACCEPTING the step",Zero);
			// take the step
			nextCandidate().setAccepted(true);
			_MCaccepted.append(nextCandidate());
			appendLinkForConvergenceTest(nextCandidate());
			if (nextCandidate().chisq()<_bestFit.chisq()) {
				addNewBestFitLink(nextCandidate());
			}

			_walk.poorX2improvement=(current().chisq()-nextCandidate().chisq())/_bestFit.chisq()<_convergence.relChisqChangeThres;
			
			if (_walk.poorX2improvement) { 
//				if (!getUphillGradient()) {
					_walk.rejectionsCount++;
					msgs->say("Poor X2 improvement detected. Counting as rejection.",Low);
//				}
			}
			else {
				/*
				 * Comment: if the improvement is very small cancel the counter; 
				 * This could be because the convergence cannot be reached because the step size is still
				 * too big, system too hot, and the resulting chisq variance does not meeet the assumed 
				 * convergence criterion. If this is the case then counting the small improvements as
				 * rejections will speed up the forced coolings, and thus speed up smaller steps.
				 * 
				 * author: blew
				 * date: Jan 9, 2018 2:54:09 PM
				 *
				 */
				
				_walk.rejectionsCount=0;
			}

			//
			// cool down
			//
			if (!_walk.poorX2improvement) {
				msgs->say("Cooling down",Low);
				if (length() >= getBurnInLength()) { // do not cool at all until the chain reaches minimal length suitable for testing the convergence
					if (getTemperature() > getFinalTemperature()) {
						//						deltaX=deltaChisq(length()-1); // commented out on Mar 17, 2011, 3:50:41 PM since we have accepted and rejected on separate lists.
						deltaX=current().chisq()-nextCandidate().chisq();
						//						msgs->say("Delta chisq is: "+msgs->toStr(deltaX),Medium);
						//
						// here is the cooling scheme
						//
						if (deltaX > 0) {
							switch (getCoolingScheme()) {
								case cpedsMCMC::coolingScheme_geometric:
									coolDownGeometric();
									break;
								case cpedsMCMC::coolingScheme_chisqPropLinLog:
									coolDownLogLin(deltaX);
									break;
								case cpedsMCMC::coolingScheme_chisqPerDoFprop:
									coolDownProp(deltaX);
									break;
								case cpedsMCMC::coolingScheme_none:
									break;
								default:
									break;
							}
							_acceptWorseRNs.clear();
						}
					}
					else { msgs->say("Cannot cool more. Reached the minimal temperature already Tmin= "+msgs->toStr(_cooling.finalTemperature),Low); }
					//				}
				}
				else { msgs->say("MC chain too short to cool down (accepted links: "+msgs->toStr(_MCaccepted.length())+", rejected links: "+msgs->toStr(_MCrej.length())+", together: "+msgs->toStr(length())+"). The minimal length to start cooling is: "+msgs->toStr(getBurnInLength()),Low); }
			}
			
			current()=nextCandidate();
		}
		else {
			msgs->say("the new chisq is > than the current one",Zero);
			// take the step with probability according to the cooling PDF encoded in coolingRNG
#ifdef DEBUGMCMC
			// DEBUG STUFF - BEGIN
			if (length()>getBurnInLength()) {
				for (long i = 0; i < 1000; i++) {
					acceptWorseRN=_coolingRNG.getRN();
					_acceptWorseRNs.append(acceptWorseRN);					
				}
				_acceptWorseRNs.save("testRNs");
				_acceptWorseRNs.clear();
				_acceptWorseRNs.append(acceptWorseThreshold);
				_acceptWorseRNs.save("testRNs-thres");
				exit(0);
			}
			else {
				acceptWorseRN=_coolingRNG.getRN();
				_acceptWorseRNs.append(acceptWorseRN);
				
			}
#else
			acceptWorseRN=_cooling.coolingRNG.getRN();
			_acceptWorseRNs.append(acceptWorseRN);
#endif
			if (_cooling.acceptWorseStateFromUniformDistr) {
				if (acceptWorseRN > acceptWorseThreshold) { // take the step
					acceptWorseState(true,nextCandidate(),_walk.rejectionsCount);
					msgs->say("Accepting worse state (chisq: %lE)",nextCandidate().chisq(),Low);
				}
				else acceptWorseState(false,nextCandidate(),_walk.rejectionsCount);
			}
			else {
				deltaX2perDOF=(nextCandidate().chisq()-current().chisq())/dataSize();
				if (acceptWorseBoltzmannFactor*acceptWorseRN > deltaX2perDOF) { // take the step
					acceptWorseState(true,nextCandidate(),_walk.rejectionsCount);
					msgs->say("Accepting worse state (chisq: %lE)",nextCandidate().chisq(),Low);
				}
				else acceptWorseState(false,nextCandidate(),_walk.rejectionsCount);
			}
		}
		
		_cooling.temperatureHistory.append(getTemperature());
		//		printInfo();
		
		//
		// saving section
		//
//#pragma omp critical
		{
		dumpAll();
		}
		
		//
		// test convergence section 
		//
		if (length()==getBurnInLength()) { 
			msgs->say("Starting testing convergence at level: "+msgs->toStr(_convergence.relChisqChangeThres),Medium);
			_convergence.last_chisq=_bestFit.chisq(); 
			_convergence.nextConvergenceCheck=length()+_convergence.statesToConvergenceCheck;
//			msgs->say("Next convergence check at: %li",_convergence.nextConvergenceCheck,Medium);
		}
		if (length()==_convergence.nextConvergenceCheck) {
//			printf("bfX2: %lE bfX2last: %lE\n",_bestFit.chisq(),_convergence.last_chisq);
#ifdef DEBUGMCMC_CONVERGENCE
			_convergence.chisq.save("convergence.chisq");
			printf("chisq std/bf_chisq: %lE, bf chisq: %lE\n",fabs(sqrt(_convergence.chisq.variance()))/_bestFit.chisq(), _convergence.relChisqChangeThres);
#endif			
			if (getX2Convergence() < _convergence.relChisqChangeThres) {
				if (coolingPossible()) { 
					/* if the system is too hot then we can be here because the generated 
						steps are too big and MCMC departs too far from the bestFit solution and 
						the system won't have a chance to cool down. Cooling down here will force smaller steps 
						and will possibly improve the fit */
					msgs->say("Convergence (deltaX2<%lf) reached at T>Tfinal",_convergence.relChisqChangeThres,Medium);
					coolForcibly();
				}
				else {
//					walk=false;
					walkDone=true;
					msgs->say("Breaking the chain. Convergence (deltaX2<%lf) reached at T=Tfinal",_convergence.relChisqChangeThres,Medium);
				}
			}
			_convergence.last_chisq=_bestFit.chisq();
			_convergence.nextConvergenceCheck+=_convergence.statesToConvergenceCheck;
		}

		// test for maximal rejections count
		if (_walk.rejectionsCount>=_cooling.maximalRejectionsCount and length()>getBurnInLength()) {
			if (forcedCoolingPossible()) { 
				msgs->say("Maximal rejections count (%li) reached at T>Tfinal",_walk.rejectionsCount,Medium);
				/* if the system is too hot then we can be here because the generated 
					steps are too big and MCMC departs too far from the bestFit solution and 
					the system won't have a chance to cool down. Cooling down here will force smaller steps 
					and will possibly further improve the fit */
				coolForcibly();
			}
			else {
//				walk=false;
				walkDone=true;
			}
		}
		

		
		if (walkDone) {
			walk=false;
/*
			//
			// burn-out
			//
			if (statesTotBeforeBurnOut==0) {
				msgs->say("Starting burn-out states",Medium);
				statesTotBeforeBurnOut=length();
				setUpHillClimbing(false);
				setInitialWalkStepSize(1);
			}
			burnOutStates++;


			if (burnOutStates==getBurnOutLength()) walk=false;
*/
		}

		
		//
		// control the maximal chain length
		//		
		if (length()==_walk.maximalChainLength) walk=false;

		
		
		//
		// output
		//
//		printf("step param0: %lE\n",getMCsteps(0).last());
		
		if (length() % getWalkInfoOutputFrequency() == 0 or walkDone==true or walk==false or msgs->getVerbosity()>=High)
			msgs->say("MC len.: %.0f, next.X2: %lE, bf.X2: %lE, cur./min.T: %lf/%f, next conv.chk: %.0f ",double(length()),nextCandidate().chisq(),_bestFit.chisq(),getTemperature(), getFinalTemperature(),double(_convergence.nextConvergenceCheck),Medium);
	
	}
	if (_walk.rejectionsCount>=_cooling.maximalRejectionsCount) { msgs->say("Breaking the chain. _maximalRejectionsCount reached: "+msgs->toStr(_walk.rejectionsCount),Medium); }
	if (length()==_walk.maximalChainLength) { msgs->say("Breaking the chain. The chain length reached: "+msgs->toStr(_walk.maximalChainLength),Medium); }
	
	//
	// print best fit in this chain
	//
	bestFitLink().printLink();
	
	
	//
	// save step PDF after walk finishes
	//
	saveStepPDFs(getPartialFilesDirFull()+"/stepPDF-walkEnd");				
	saveTemperature(getPartialFilesDirFull()+"/temperature-walkEnd");
	
	//
	// burn-out
	//
	statesTotBeforeBurnOut=length();
	_walk.statesTotBeforeBurnOut=length();
	
	burnOut();
	
	//
	// save step sizes
	//
	saveMCsteps(getPartialFilesDirFull()+"/step-sizes");
	saveBestFit();
	
	saveChain(_MCrej,getPartialFilesDirFull()+"/rejected.mc"); 
	saveChain(_MCaccepted,getPartialFilesDirFull()+"/accepted.mc"); 
//	QList<MClink> tmp=_MCrej; tmp.append(_MCaccepted);
//	saveChain(tmp,getOutputDir()+"/all.mc"); 


	
	msgs->say("Best fit link is:", Medium);
	printLink(_bestFit);
	
	// save everything we've got now
}

/***************************************************************************************/
void cpedsMCMC::updateStepGeneratingPDFs() {
	msgs->say("Updating step generating PDFs",Low);
	mscsFunction f("step generating PDF",Zero);
	double frac;
	frac=_walk.initialPDFstepCRfraction;
	if (length()>=getBurnInLength()) {
		if (_walk.initialPDFstepCRfractionAfterBurnIn>0) frac=_walk.initialPDFstepCRfractionAfterBurnIn;
	}
	
	for (long i = 0; i < dims(); i++) {
		f.setName(_parameterSpace[i].getName()+"PDF");
		// define step size generating distribution
//		f.mkBoltzmann(0,_cooling.initialMaximalEnergy*CPEDS_kB*_cooling.initialTemperature,_cooling.initialMaximalEnergy*CPEDS_kB*getInitialTemperature()/_cooling.stepGeneratingPDFPoints,getTemperature());
		f.mkBoltzmann(0,_cooling.initialMaximalEnergy*CPEDS_kB*getTemperature(),_cooling.initialMaximalEnergy*CPEDS_kB*getTemperature()/_cooling.stepGeneratingPDFPoints,getTemperature());
//		f.scaleX(frac*(_parameterSpace[i].getMaxArg()-_parameterSpace[i].getMinArg())/f.getMaxArg()); // make the maximal step size equal to half of parameter space prior - this is because the step can be positive or negative
		f.scaleX(frac*getTemperature()/getInitialTemperature()*(_parameterSpace[i].getMaxArg()-_parameterSpace[i].getMinArg())/f.getMaxArg()); 
		/*
		 * Comment: we generate step PDF basically from the same initial distribution, so there is no need to re-generate this every time 
		 * and for every parameter - but it's not a big deal in terms of efficiency.
		 * The improvement to the commented lines is that new we only probe the part of the PDF that is useful (i.e. non-zero) which
		 * is useful when temperatures become very small. In the previous implementation the PDFs became more like a step functions.
		 * 
		 * So, now we essentially rescale the same PDF according to the requested expected step size fraction (frac) and a proportion
		 * of the current temperature to the initial temperature.
		 * 
		 * author: blew
		 * date: Jan 9, 2018 7:33:11 PM
		 *
		 */
		
		f.normalize(); // normalize
		_walk.stepGeneratingDistr[i]=f;

/*
		f.checkRanges();
		
		mscsFunction g=f;
		g/=g.getMaxValue();
		g.save("step-param_"+msgs->toStr(i)+".mc");
		g.invert(); g.sortFunctionArgAscending();
		long xi;
		double tmp=g.f(0.25,&xi);
		printf("expected step size in direction: %li is %lE, xi: %li\n",i,f.getX(f.pointsCount()-xi),f.pointsCount()-xi);
*/
		
		// define step generating RNG
		_walk.stepGeneratingRNG[i].setPDF(f.pointsCount(),f.extractArguments(),f.extractValues());	
		f.clearFunction();
	}
}

/***************************************************************************************/
void cpedsMCMC::updateCoolingPDF() {
	if (getAcceptWorseScheme()==cpedsMCMC::awSchemeBoltzmann) {
		// update cooling PDF
		msgs->say("Updating the cooling PDF",Low);
		mscsFunction f("coolingPDF",Zero);
		//	f.mkBoltzmann(0,_initialMaximalEnergy*CPEDS_kB*_initialTemperature,_initialMaximalEnergy*CPEDS_kB*_initialTemperature/_stepGeneratingPDFPoints,_temperature);
		//		f.mkBoltzmann(0,_initialMaximalEnergy*CPEDS_kB*_temperature,_initialMaximalEnergy*CPEDS_kB*_temperature/_stepGeneratingPDFPoints,_temperature);
		f.mkBoltzmann(0,_cooling.initialMaximalEnergy*CPEDS_kB*_cooling.initialTemperature,_cooling.initialMaximalEnergy*CPEDS_kB*getInitialTemperature()/_cooling.stepGeneratingPDFPoints,getTemperature());
		//		if (_acceptWorseStateFromUniformDistr==false) { // i.e. actually use the boltzmann distribution
		// scale the distribution so it make sense to compare the arguments with delta chisq
		f.scaleX(1.0/(CPEDS_kB*getInitialTemperature()));
		//		}
		f.normalize();
		_cooling.coolingDistr=f;
		_cooling.coolingDistrInvCDF=_cooling.coolingDistr;
		_cooling.coolingDistrInvCDF.setVerbosityLevel(Zero);
		_cooling.coolingDistrInvCDF.mkCDF(true);
		_cooling.coolingDistrInvCDF.invert();
		
		// update cooling RNG
		_cooling.coolingRNG.setPDF(f.pointsCount(),f.extractArguments(),f.extractValues());
		//		if (getAcceptWorseScheme()==cpedsMCMC::awSchemeUniform) {
		//			updateAcceptWorseThresholds();
		//		}
	}
	
}
/***************************************************************************************/
double* cpedsMCMC::getStartingPoint() {
	double* p=new double[dims()];
//	printf("_walk.rns size is: %li and dims: %li\n",_walk.rns.size(),dims());
	for (long i = 0; i < dims(); i++) {		
		p[i]=_walk.rns[i].getRN();	}
	return p;
}
/***************************************************************************************/
MClink cpedsMCMC::getNextPoint(const MClink& current) {
	msgs->say("generating the next point", Zero);
#ifdef DEBUGMCMC
	msgs->say("current state is:", Low);
	//	if (msgs->getVerbosity()>
	printLink(current);
	msgs->say("Best fit link is:", Low);
	printLink(_bestFit);	
#endif
	MClink newLink(dims());
	//	double f; // step shrink factor
	double delta;
	double sign;
	double tmp;
	long failures;
	for (long i = 0; i < dims(); i++) {
		//		f=_temperature/_initialTemperature;
		failures=-1;
//		do {
			failures++;
			delta=_walk.stepGeneratingRNG[i].getRN();
//			printf("step size along %li: %lE\n",i,delta);
//			getWalk().steps.append(delta);
			
			if (_walk.uphillClimbing and _walk.keepDirectionNow and failures<_walk.keepDirectionStepGenerationFailuresNumber) {
				sign=_walk.currentDirection[i];
			}
			else {
				sign=getSign();
				_walk.currentDirection[i]=sign;
				if (_walk.uphillClimbing and _walk.keepDirectionNow) msgs->say("Reached maximal rejections count in walk along parameter: "+msgs->toStr(i)+". Will walk in any direction from now on.",Low);
			}
			
			//			msgs->say("generating new step for dimension: "+msgs->toStr(i)+",delta is: "+msgs->toStr(delta), Low);
			//		p[i]=current[i] + getSign()*(_delta0[i]*f+delta);	
			tmp = current[i] + sign*(delta);
			
			MCsteps(i).append(sign*(delta)); // store the step
			//			msgs->say("new candidate: "+msgs->toStr(tmp)+", boundL: "+msgs->toStr(_parameterSpace[i].getX(0))+" boundU: "+msgs->toStr(_parameterSpace[i].getX(_parameterSpace[i].pointsCount()-1)), Low);

		
			/*
			 * Comment: This loop is supposed to ensure that the chain won't walk outside of the predefined
			 * parameter space. But the implementation is wrong and sometimes leads to never ending cycles.
			 * So, temporarily commented out.
			 * Besides, it's not clear whether this functionality is really useful. Perhaps it is, so
			 * it could be considered for implementing as an optional switchable parameter
			 * 
			 * author: blew
			 * date: Jan 9, 2018 3:53:50 PM
			 *
			 */
			
//		} while (tmp < _parameterSpace[i].getX(0) || tmp > _parameterSpace[i].getX(_parameterSpace[i].pointsCount()-1));
//		newLink.setParam(i,tmp);
		newLink[i]=tmp;
		//		msgs->say("accepted", Low);
	}
	_walk.keepDirectionNow=false;
	
	return newLink;	
}



/***************************************************************************************/
double cpedsMCMC::getSign() {
	double s=_walk.rnSgn.getRN();
	if (s<0) return -1;
	return 1;	
}
/* ******************************************************************************************** */
double cpedsMCMC::getX2Convergence() {
	return fabs(sqrt(_convergence.chisq.variance()))/_bestFit.chisq();
}
/* ******************************************************************************************** */
double cpedsMCMC::getRelX2Improvement() {
//	return current().chisq()-next.chisq()/_bestFit.chisq();
}
/***************************************************************************************/
void cpedsMCMC::setData(const mscsFunction& data) {	
	msgs->say("Setting data vector of length: "+msgs->toStr(data.pointsCount()),Medium);
	_chisqData.data=data; 	
	_chisqData.model=data; // this is to set the size and arguments
}
/* ******************************************************************************************** */
void cpedsMCMC::setData2D(const mscsFunction3dregc& data, int colNo) {
	_chisqData.data=data.getSlice1Dfn(1,0,colNo,0);	
	_chisqData.data2dCol=colNo;
	_chisqData.data2d.setSize(data.Nx()+2,data.Ny()); // two extra columns for best fit model and residuals
	_chisqData.data2d.allocFunctionSpace();
	for (long i = 0; i < data.Nx(); i++) {
		for (long j = 0; j < data.Ny(); j++) {
			_chisqData.data2d(i,j)=data(i,j);
		}
	}
}

/***************************************************************************************/
const mscsFunction& cpedsMCMC::getData() const {
	return _chisqData.data;
}
/***************************************************************************************/
const mscsFunction3dregc& cpedsMCMC::getData2D() const {
	return _chisqData.data2d;
}
/***************************************************************************************/
const double cpedsMCMC::getData2D(long i, long j) const {
	return _chisqData.data2d(i,j);
}
/***************************************************************************************/
const mscsFunction& cpedsMCMC::getModel() const {
	return _chisqData.model;
}
/***************************************************************************************/
void cpedsMCMC::setDiagonalCovarianceMatrix(const cpedsList<double>& cov) {
//	_covDiag=cov.toCarray();
	_chisqData.covarianceMatrixDiag=cov;
}
/***************************************************************************************/
void cpedsMCMC::setDiagonalCovarianceMatrix(double var) { _chisqData.variance=var; }
/***************************************************************************************/
void cpedsMCMC::setCovarianceMatrix(const mscsFunction3dregc& cov) {
	_chisqData.covarianceMatrix=cov;
}
/***************************************************************************************/

void cpedsMCMC::printInfo() const {
	msgs->say("MC setup:",Top);
	msgs->say("MC walk parameters:",High);
	msgs->say("Total number of visited states: "+msgs->toStr(_walk.statesTot),Medium);
	msgs->say("initial temperature: "+msgs->toStr(getInitialTemperature()),Medium);
	msgs->say("final temperature: "+msgs->toStr(getFinalTemperature()),Medium);
	msgs->say("temperature: "+msgs->toStr(getTemperature()),Medium);
	msgs->say("initial maximal energy: "+msgs->toStr(_cooling.initialMaximalEnergy),Medium);
	msgs->say("step generating PDF points: "+msgs->toStr(_cooling.stepGeneratingPDFPoints),Medium);
	msgs->say("Burn-in stage step PDF domain fraction (_walk.initialPDFstepCRfraction): "+msgs->toStr(_walk.initialPDFstepCRfraction),Medium);
	msgs->say("After Burn-in stage step PDF domain fraction (_walk.initialPDFstepCRfractionAfterBurnIn): "+msgs->toStr(_walk.initialPDFstepCRfractionAfterBurnIn),Medium);
	msgs->say("user output frequency: "+msgs->toStr(getWalkInfoOutputFrequency()),Medium);
	msgs->say("Burn-in period length (number of steps to start cooling down): "+msgs->toStr(getBurnInLength(),"%li"),Medium);
	msgs->say("maximal number of times to force cooling: "+msgs->toStr(_cooling.maximalNumberOfForcedCoolings),Medium);
	msgs->say("maximal number of rejections to force cooling or break chain: "+msgs->toStr(_cooling.maximalRejectionsCount),Medium);
	msgs->say("maximal chain length: "+msgs->toStr(_walk.maximalChainLength),Medium);
	msgs->say("burn-in length: "+msgs->toStr(getBurnInLength()),Medium);
	msgs->say("burn-out length: "+msgs->toStr(getBurnOutLength()),Medium);
	msgs->say("cooling rate: "+msgs->toStr(getCoolingRate()),Medium);
	
	msgs->say("accept worse state uniform scheme: "+msgs->toStr(_cooling.acceptWorseStateFromUniformDistr),Medium);
	if (getAcceptWorseScheme()==awSchemeUniform) {
		msgs->say("accept worse state p-value: "+msgs->toStr(_cooling.acceptWorsePvalue),Medium);
		msgs->say("accept worse state threshold: "+msgs->toStr(_cooling.acceptWorseThreshold),Medium);
		msgs->say("accept worse state p-value after burn-in: "+msgs->toStr(_cooling.acceptWorsePvalueAfterBurnIn),Medium);
		msgs->say("accept worse state threshold after burn-in: "+msgs->toStr(_cooling.acceptWorseThresholdAfterBurnIn),Medium);
	}
	if (getAcceptWorseScheme()==awSchemeBoltzmann) {
		msgs->say("accept worse RN multiplicative factor during burn-in period: "+msgs->toStr(_cooling.acceptWorseBoltzmannFact),Medium);
		msgs->say("accept worse RN multiplicative factor after burn-in period: "+msgs->toStr(_cooling.acceptWorseAfterBurnInBoltzmannFact),Medium);
	}
	msgs->say("up-hill climbing: "+msgs->toStr(_walk.uphillClimbing),Medium);
	msgs->say("number of failures to generate the next step in the same direction before changing direction: "+msgs->toStr(_walk.keepDirectionStepGenerationFailuresNumber),Medium);
	msgs->say("chisq computing function signature: "+_IOcontrol.chisqSignature,Medium);
	msgs->say("best fit cov simulations number: "+msgs->toStr(_chisqData.bestFitNcovSimNum),Medium);
	
	
	msgs->say("RNG options:",Top);
	msgs->say("run seed: "+msgs->toStr(_walk.runRNGseed),Medium);
	msgs->say("run offset: "+msgs->toStr(_IOcontrol.runOffset),Medium);
	
	msgs->say("MC save options:",Top);
	msgs->say("MCMC run index: "+msgs->toStr(_IOcontrol.mcmcRunIdx),Medium);
	msgs->say("run files directory: "+getOutputDir(),Medium);
	msgs->say("partial files directory: "+getPartialFilesDirFull(),Medium);
	msgs->say("save cooling PDF every: "+msgs->toStr(_IOcontrol.saveCoolingPDFEveryNthState),Medium);
	msgs->say("save step PDF every: "+msgs->toStr(_IOcontrol.saveStepPDFsEveryNthState),Medium);
	msgs->say("save model every: "+msgs->toStr(_IOcontrol.saveModelEveryNthState),Medium);
	msgs->say("dump every MCMC state: "+msgs->toStr(_IOcontrol.dumpEveryMCMCstate),Medium);
	
	msgs->say("PARAMETER SPACE:",Top);
	msgs->say("parameter names:",High);
	for (long i = 0; i < dims(); i++) {
		msgs->say("name of parameter "+msgs->toStr(i)+": "+_parameterSpace[i].getName(),Medium);
	}
	msgs->say("parameter ranges:",High);
	for (long i = 0; i < dims(); i++) {
		msgs->say("parameter "+msgs->toStr(i)+" range> from: "+msgs->toStr(_parameterSpace[i].getMinArg()),Medium);
		msgs->say("parameter "+msgs->toStr(i)+" range> to  : "+msgs->toStr(_parameterSpace[i].getMaxArg()),Medium);
	}
	
	msgs->say("Starting link:",High);
	printLink(_startingLink);
	msgs->println("",Low);
	
	msgs->say("IMPLEMENTATION NOTES BLOCK:",Top);
	for (long i = 0; i < _IOcontrol.implementationNote.count(); i++) {
		if (_IOcontrol.implementationNote[i]=="") {
			msgs->say(_IOcontrol.implementationNote[i+1],High);			
			i++;
		}
		else {
			msgs->say("note "+msgs->toStr(i)+": "+_IOcontrol.implementationNote[i],Medium);
		}
	}
	
	
}

/***************************************************************************************/
void cpedsMCMC::saveAcceptedChisq(string fname) {
	cpedsList<double> l;
	
	for (long i = 0; i < _MCaccepted.size(); i++) {
		l.append(_MCaccepted[i].chisq());
	}
	
	l.save(fname);
	
}
/* ******************************************************************************************** */
cpedsMC cpedsMCMC::loadMC(string fileName) {
	cpedsMC chain;
	chain.load(fileName);
	return chain;
}
/* ******************************************************************************************** */
void cpedsMCMC::loadMCrun(string dirName) {
	bestFitLink().load(dirName+"/"+getPartialFilesDir()+"/bestFit-link");
#ifdef DEBUG_MCMC_PDF2
	cout << "loading best fit link:\n";
	bestFitLink().printLink();	
#endif
	_walk.Nparam=bestFitLink().dims();

//	printf("dims: %i\n",dims());
	initialize(-1,0,dims(),0);
//	printf("dims: %i\n",dims());
	_chisqData.data.load(dirName+"/inputData.txt");
	_MCaccepted.load(dirName+"/accepted.mc");
	_MCrej.load(dirName+"/rejected.mc");
	_MCbestFitChain.load(dirName+"/bestFitChain.mc");

	loadBestFitStepPDFs(dirName);
	
	
	// load parameter space
	for (long i = 0; i < dims(); i++) {
		mscsFunction pspace;
//		psapce.load(dirName+"/parameterPrior.param."+msgs->toStr(i));
		pspace.loadHDF5(dirName+"/parameterPriors.hdf5",msgs->toStr(i));
		int errCode=0;
		string pname=pspace.getHDF5_stringAttribute(dirName+"/parameterPriors.hdf5",msgs->toStr(i),"paramName",&errCode);
		if (errCode!=0) msgs->say("could not read parameter name from file: "+dirName+"/parameterPriors.hdf5 (param "+msgs->toStr(i)+")",High);
		string pnameFull=pspace.getHDF5_stringAttribute(dirName+"/parameterPriors.hdf5",msgs->toStr(i),"paramNameFull",&errCode);
		if (errCode!=0) msgs->say("could not read full parameter name from file: "+dirName+"/parameterPriors.hdf5 (param "+msgs->toStr(i)+")",High);
		parameterSpace().addParameter(pspace,pname,pnameFull);
//		printf("adding parameter prior -- LOADING PARAMETER SPACE IS NOT FULLY IMPLEMENTED, THIS WILL NOT LOAD PARAMETER NAMES, FIX IT\n");

	
		mscsFunction stepPDF;
		stepPDF.load(dirName+"/"+getPartialFilesDir()+"/stepPDFend-param_"+msgs->toStr(i)+".mc");
		_walk.stepGeneratingDistr.append(stepPDF);

	}
	
	
/*
	
	// append parts of structures 
	for (long i = 0; i < getNparam(); i++) {
		_walk.steps[i].append(chain.getMCsteps(i));
	}
*/

}
/***************************************************************************************/
void cpedsMCMC::saveParams(string fname) {
	if (_MCaccepted.size()>0) {
		matrix<double> m(_MCaccepted.size(),dims());
		
		for (long i = 0; i < _MCaccepted.size(); i++) {
			for (long j = 0; j < dims(); j++) {
				m(i,j)=_MCaccepted[i].getParam(j);			
			}
		}
		
		cpeds_matrix_save(m,fname);
	}
}
/***************************************************************************************/
void cpedsMCMC::saveData(string fname) {
	_chisqData.data.save(fname);
}
/***************************************************************************************/
void cpedsMCMC::savePriors() {
	for (long i = 0; i < dims(); i++) {
//		_parameterSpace[i].save(getOutputDir()+"/parameterPrior.param."+msgs->toStr(i));
		_parameterSpace[i].saveHDF5(getOutputDir()+"/parameterPriors.hdf5",msgs->toStr(i));
		_parameterSpace[i].setHDF5_scalarStringAttribute(getOutputDir()+"/parameterPriors.hdf5",msgs->toStr(i),"paramName",_parameterSpace.names()[i]);
		_parameterSpace[i].setHDF5_scalarStringAttribute(getOutputDir()+"/parameterPriors.hdf5",msgs->toStr(i),"paramNameFull",_parameterSpace.names_full()[i]);		
	}
}

/* ******************************************************************************************** */
string cpedsMCMC::getRunDependentFileName(string fname) {
	string outfile=fname+".mc";
	
	if (getRunIdx()!=-1) outfile+=msgs->toStr(getRunIdx()); 
	return outfile;
}
/* ******************************************************************************************** */
string cpedsMCMC::getParameterDependentFileName(string fname, long paramID) {
	return fname+"-param_"+msgs->toStr(paramID);
}
/* ******************************************************************************************** */
string cpedsMCMC::get2ParameterDependentFileName(string fname, long paramID1, long paramID2) {
	return fname+"-params"+msgs->toStr(paramID1)+"_"+msgs->toStr(paramID2);
}
/***************************************************************************************/
void cpedsMCMC::saveTemperature(string fname) {
	if (_cooling.temperatureHistory.size()>0) _cooling.temperatureHistory.save(getRunDependentFileName(fname));
}

/***************************************************************************************/
void cpedsMCMC::saveStepPDFs(string fname) {
	string outfile;
	outfile=getRunDependentFileName(fname);
	msgs->say("Saving step PDF to file: "+outfile,Medium);
#pragma omp critical 
	for (long j = 0; j < dims(); j++) {
//		outfile=getParameterDependentFileName(fname,j);
//		outfile=getRunDependentFileName(outfile);
//		_walk.stepGeneratingDistr[j].save(outfile);

		
		_walk.stepGeneratingDistr[j].saveHDF5(outfile+".hdf5",msgs->toStr(j));
	}
	_IOcontrol.saveStepPDFs_lastSavedNumber++;

}
/* ******************************************************************************************** */
cpedsList<double> cpedsMCMC::getMCsteps(long param) {
	assert(param < _walk.steps.size());
	return _walk.steps[param];
}
/* ******************************************************************************************** */
cpedsList<double>& cpedsMCMC::MCsteps(long param) {
	assert(param < _walk.steps.size());
	return _walk.steps[param];	
}
/* ******************************************************************************************** */
void cpedsMCMC::saveMCsteps(string fname) {
	long irej=0,iacc=0;
	string outfile;

	outfile=getRunDependentFileName(fname);
#pragma omp critical 	
	for (long i = 0; i < getNparam(); i++) {
//		getMCsteps(i).save(fname+"param_"+msgs->toStr(i,"%li")+".mc");
		
		mscsFunction3dregc steps;
		steps.setSize(1,getMCsteps(i).size());
		steps.allocFunctionSpace();
		
		steps.saveHDF5(outfile+".hdf5",msgs->toStr(i));
	}
}
/***************************************************************************************/
void cpedsMCMC::saveBestFitStepPDFs(string fname) {
	string outfile;
	outfile=getRunDependentFileName(fname);
	msgs->say("Saving step PDF to file: "+outfile,Medium);
#pragma omp critical 
	for (long j = 0; j < dims(); j++) {
		
//		outfile=getParameterDependentFileName(fname,j);
//		outfile=getRunDependentFileName(outfile);
		
		//		_chisqData.bestFitData.stepGeneratingDistr[j].save(outfile);
		_chisqData.bestFitData.stepGeneratingDistr[j].saveHDF5(outfile+".hdf5",msgs->toStr(j));
	}
}
/***************************************************************************************/
void cpedsMCMC::loadBestFitStepPDFs(string dirName) {
	mscsFunction pdf;
	string outfile;
	for (long j = 0; j < dims(); j++) {
/*
		outfile=getParameterDependentFileName(dirName+"/"+getPartialFilesDir()+"/bestFit-stepPDF",j);
		outfile=getRunDependentFileName(outfile);
		msgs->say("Loading step PDF from file: "+outfile,Medium);
		pdf.load(outfile);
*/
		
		outfile=getRunDependentFileName(dirName+"/"+getPartialFilesDir()+"/bestFit-stepPDF");
		msgs->say("Loading step PDF from file: "+outfile,Medium);
		pdf.loadHDF5(outfile,msgs->toStr(j));
		
		_chisqData.bestFitData.stepGeneratingDistr.append(pdf);
	}
}
/***************************************************************************************/
void cpedsMCMC::setSaveParameters(long saveCoolingPDFEvery, long saveStepPDFsEvery, long saveModelEvery) {
	_IOcontrol.saveCoolingPDFEveryNthState=saveCoolingPDFEvery; 
	_IOcontrol.saveStepPDFsEveryNthState=saveStepPDFsEvery;
	_IOcontrol.saveModelEveryNthState=saveModelEvery;	
	msgs->say("Setting save parameters to store Cooling PDF every "+msgs->toStr(_IOcontrol.saveCoolingPDFEveryNthState)+" states",Low);
	msgs->say("Setting save parameters to store Step PDFs every "+msgs->toStr(_IOcontrol.saveStepPDFsEveryNthState)+" states",Low);
	msgs->say("Setting save parameters to store Model every "+msgs->toStr(_IOcontrol.saveModelEveryNthState)+" states",Low);
}
/***************************************************************************************/
void cpedsMCMC::saveParameterNames(string fname) {
	string s;
	if (fname=="") s="parameter_names.mc";
	else s=fname;
	
	FILE* f=fopen(s.c_str(),"w");
	for (long j = 0; j < dims(); j++) {
		fprintf(f,"%li %s\n",j,_parameterSpace[j].getName().c_str());
	}
	fclose(f);
}

/***************************************************************************************/
const mscsFunction& cpedsMCMC::getParameter(string Pname, int* ctrl) const {
	for (long i = 0; i < dims(); i++) {	
		if (_parameterSpace[i].getName()==Pname) { if (ctrl!=NULL) *ctrl = 0; return _parameterSpace[i]; }
	}
	if (ctrl!=NULL) *ctrl=-1;
	return _parameterSpace[0];
}
/***************************************************************************************/
const int cpedsMCMC::getParameterByName(string Pname) const {
	for (long i = 0; i < dims(); i++) {	
		if (_parameterSpace[i].getName()==Pname) { return i; }
	}	
	return -1;
}
/* ******************************************************************************************** */
void cpedsMCMC::saveCooling() {
	string tmps;

	cpedsList<double> tmpl;
	tmps=msgs->toStr(_IOcontrol.saveCoolingPDF_lastSavedNumber,"%05li");
	_cooling.coolingDistr.save(getPartialFilesDirFull()+"/"+tmps+"-coolingPDF.step_"+msgs->toStr(_MCaccepted.size())); 
	_cooling.coolingDistrInvCDF.save(getPartialFilesDirFull()+"/"+tmps+"-coolingInvCDF.step_"+msgs->toStr(_MCaccepted.size())); 
	_acceptWorseRNs.save(getPartialFilesDirFull()+"/"+tmps+"-acceptWorseRNs.step_"+msgs->toStr(_MCaccepted.size())); 
	tmpl.append(_cooling.acceptWorseThreshold);
	tmpl.append(_cooling.acceptWorseThresholdAfterBurnIn);
	tmpl.save(getPartialFilesDirFull()+"/"+tmps+"-acceptWorseThresholds.step_"+msgs->toStr(_MCaccepted.size())); 
	
	_IOcontrol.saveCoolingPDF_lastSavedNumber++;	
}
/* ******************************************************************************************** */
void cpedsMCMC::saveModel() {
	string tmps;

	tmps=msgs->toStr(_IOcontrol.saveModel_lastSavedNumber,"%05li");
	_chisqData.model.save(getPartialFilesDirFull()+"/"+tmps+"-model.step_"+msgs->toStr(_MCaccepted.size())); 
	_chisqData.chisqPerDOF.save(getPartialFilesDirFull()+"/"+tmps+"-chisqPerDOF.step_"+msgs->toStr(_MCaccepted.size()));
	_IOcontrol.saveModel_lastSavedNumber++;				

	
}

/***************************************************************************************/
/* NOTE: The indexing of partial files will be smaller by 1 wrt states number (statesTotal variable)
 * (for the case of dumping at every state), because the 
 * zero'th state is not saved, and the indexing variables used in this method work 
 * independently */

void cpedsMCMC::dumpAll() {
	char tmpch[10];
	string tmps;
	if (_cooling.acceptWorseStateFromUniformDistr==false) {
		if (_IOcontrol.saveCoolingPDFEveryNthState!=0) {
			if (_MCaccepted.size() % _IOcontrol.saveCoolingPDFEveryNthState == 0) {
				saveCooling();
			}					
		}
	}

	if (_IOcontrol.saveStepPDFsEveryNthState!=0)
		if (_MCaccepted.size() % _IOcontrol.saveStepPDFsEveryNthState== 0) { 
			tmps=msgs->toStr(_IOcontrol.saveStepPDFs_lastSavedNumber,"%05li");
			saveStepPDFs(getPartialFilesDirFull()+"/"+tmps+"-stepPDFs.step_"+msgs->toStr(_MCaccepted.size())); 
		}
	
	if (_IOcontrol.saveModelEveryNthState!=0)
		if (_MCaccepted.size() % _IOcontrol.saveModelEveryNthState== 0) { 
			saveModel();
		}

	if (_IOcontrol.dumpEveryMCMCstate) {
//		saveAcceptedChisq(getOutputDir()+"/chisq.mc"); 
//		saveParams(getOutputDir()+"/params.mc");
//		saveChain(_MCrej,getOutputDir()+"/rejected.mc"); 
//		saveChain(_MCaccepted,getOutputDir()+"/accepted.mc"); 
//		QList<MClink> tmp=_MCrej; tmp.append(_MCaccepted);
//		saveChain(tmp,getOutputDir()+"/all.mc"); 
		saveTemperature(getPartialFilesDirFull()+"/temperature");
		saveStepPDFs(getPartialFilesDirFull()+"/stepPDFend");				
		saveMCsteps(getPartialFilesDirFull()+"/step-sizes");
		if (_cooling.acceptWorseStateFromUniformDistr==false)
			_cooling.coolingDistr.save(getOutputDir()+"/coolingPDFend.mc");
		
		if (_chisqData.covarianceMatrixDiag.size()>0 and _chisqData.covDiagonal) {
			_chisqData.covarianceMatrixDiag.save(getOutputDir()+"/current-covarianceMatrixDiagonal");
//			cpeds_save_matrix(_covDiag,_data.pointsCount(),1,getOutputDir()+"/current-covarianceMatrixDiagonal",false);
		}
	}
	
	if (_IOcontrol.storeBestFitNow) {
		// this is now moved to storeBestFit method
		/*
		_bestFitData.coolingDistr=_coolingDistr;
		_bestFitData.stepGeneratingDistr=_stepGeneratingDistr;
		_bestFitData.model=_model;
		
		if (_covDiag!=NULL and _covDiagonal) {
			_bestFitData.covDiag.clear();
			_bestFitData.covDiag.fromCarray(_covDiag,_data.pointsCount());
		}
		 */
		
		/*
		_coolingDistr.save(_partialFilesDir+"/bestFit-coolingPDF"); 
		saveStepPDFs(_partialFilesDir+"/bestFit-stepPDFs"); 
		_model.save(_partialFilesDir+"/bestFit-model"); 			
//		cout << "saving best fit\n";
//		exit(0);
//		_temperatureHistory.save(_partialFilesDir+"/temperature-history"); 			
		_bestFit.save(_partialFilesDir+"/bestFit-link");
		saveChain(_MCbestFitChain,_runFilesDir+"/bestFitChain");
		if (_covDiag!=NULL and _covDiagonal) {
			cpeds_save_matrix(_covDiag,_data.pointsCount(),1,_partialFilesDir+"/bestFit-covarianceMatrixDiagonal",false);
		}
		 */
		_IOcontrol.storeBestFitNow=false;
	}
	
	//	saveCurrentLength();
	//	saveCurrentTemperature();
	//	saveFirstLinkLocation();
}

/***************************************************************************************/
double cpedsMCMC::coolDownGeometric() {
	msgs->say("COOLING: current temperature: "+msgs->toStr(getTemperature()),Low);
	double oldTemperature=_cooling.temperature;
//	_cooling.temperature/=pow(2.0,getCoolingRate());
	_cooling.temperature/=2.0;
	
	if (getTemperature()<getFinalTemperature()) _cooling.temperature=getFinalTemperature();
	
	msgs->say("COOLING: new temperature: "+msgs->toStr(getTemperature()),Low);
	
	updateStepGeneratingPDFs();		
	updateCoolingPDF();

	return oldTemperature-_cooling.temperature;
}
/* ******************************************************************************************** */
void cpedsMCMC::coolDownLogLin(double deltaX) {
	msgs->say("COOLING: current temperature: "+msgs->toStr(getTemperature()),Low);
	//	double deltaX2perDOF=deltaX/_data.pointsCount();
	double dT;
	if (log(deltaX) > 0) {
		dT= getCoolingRate() * log(deltaX); // the cooling is proportional to log delta chisq
	}
	else {
		dT= getCoolingRate() * deltaX; // the cooling is proportional to delta chisq							
	}
	
//	printf("Cooling loglin: deltaX: %lE, log(deltaX): %lE, cooling by: %lE, MClen: %li\n",deltaX, log(deltaX),dT,length());
	_cooling.temperature-= dT;
	
	if (getTemperature()<getFinalTemperature()) _cooling.temperature=getFinalTemperature();
	
	msgs->say("COOLING: new temperature: "+msgs->toStr(getTemperature()),Low);
	
	if (getTemperature()>0)  {
		updateStepGeneratingPDFs();		
		updateCoolingPDF();		
	}
	
}
/***************************************************************************************/
void cpedsMCMC::coolDownProp(double deltaX) {
	msgs->say("COOLING: current temperature: "+msgs->toStr(getTemperature()),Low);
	double deltaX2perDOF=deltaX/dataSize();
	double chisqIni=_MCaccepted[_walk.rejectedLinkNumber].chisq()/dataSize();
	//	double dT=_finalTemperature/_initialTemperature*deltaX2perDOF*chisqIni;
	double dT=(_cooling.finalTemperature-_cooling.initialTemperature)*deltaX2perDOF/chisqIni; // this is under construction - what goes here ?
	
	_cooling.temperature-=getCoolingRate()*dT;
	
	
	if (getTemperature()<getFinalTemperature()) _cooling.temperature=getFinalTemperature();
	
	msgs->say("COOLING: new temperature: "+msgs->toStr(getTemperature()),Low);
	
	if (getTemperature()>0)  {
		updateStepGeneratingPDFs();		
		updateCoolingPDF();
	}
	
}
/***************************************************************************************/
void cpedsMCMC::saveSetup() {
	mkOutputDir();
	msgs->loggingOn();
	msgs->messagingOff();
	msgs->setLogFileName(getOutputDir()+"/setup.mc");
	msgs->setLogVerbosity(msgs->getVerbosity());
	printInfo();
	msgs->messagingOn();
	msgs->loggingOff();
	
	_startingLink.save(getRunDependentFileName(getPartialFilesDirFull()+"/startingLink"));
	saveData(getOutputDir()+"/inputData.txt"); 
	savePriors();
	if (_cooling.acceptWorseStateFromUniformDistr==false)	
		_cooling.coolingDistr.save(getOutputDir()+"/coolingPDFini.mc");

}
/***************************************************************************************/
void cpedsMCMC::saveChain(QList<MClink>& chain, string fname) {	
	if (chain.size()>0) {
		string output=getRunDependentFileName(fname);
		matrix<double> m(chain.size(),dims()+4);
		
		for (long i = 0; i < chain.size(); i++) {
			for (long j = 0; j < dims(); j++) {
				m(i,j)=chain[i].getParam(j);			
			}
			m(i,dims())=chain[i].chisq();
			m(i,dims()+1)=chain[i].L();
			m(i,dims()+2)=chain[i].getIdx();
			m(i,dims()+3)=double((int)(chain[i].isAccepted()));
		}
		
		cpeds_matrix_save(m,output);
	}
}
/***************************************************************************************/
void cpedsMCMC::printLink(const MClink& link) const {
	msgs->say("Link Information:",Medium);
	for (long i = 0; i < link.dims(); i++) {
		msgs->say("param "+msgs->toStr(i)+": "+msgs->toStr(link[i]),Medium);
	}
	if (link.chisq()!=-1) {
		msgs->say("chisq: "+msgs->toStr(link.chisq()),Medium);
	}
	if (link.L()!=-1) {
		msgs->say("likelihood: "+msgs->toStr(link.L()),Medium);
	}
}
/***************************************************************************************/
void cpedsMCMC::mkOutputDir() {
	
	string cmd="if [ ! -d "+getPartialFilesDirFull()+" ]; then mkdir -p "+getPartialFilesDirFull()+"; fi";
	system(cmd.c_str());
/*
	string dname=getPartialFilesDirFull();
	mkdir(dname.c_str(),0755);
*/
}
/***************************************************************************************/
void cpedsMCMC::setOutputDir(string dirName) { 
	if (_IOcontrol.mcmcRunIdx==-1) {
		_IOcontrol.runFilesDir=dirName;
	}
	else _IOcontrol.runFilesDir=dirName+"."+msgs->toStr(_IOcontrol.mcmcRunIdx+_IOcontrol.dirStartIndex);
	_IOcontrol.partialFilesDir=_IOcontrol.runFilesDir+"/"+getPartialFilesDir();
}
/***************************************************************************************/
void cpedsMCMC::saveFirstLinkLocation() {
	cpedsList<long> l;
	l.append(_walk.rejectedLinkNumber);
	l.save(getPartialFilesDirFull()+"/first-link-index",false,"long");
}

/***************************************************************************************/
void cpedsMCMC::saveCurrentLength() {
	cpedsList<long> l;
	l.append(_MCaccepted.length());
	l.save(getPartialFilesDirFull()+"/current-length",false,"long");
}
/***************************************************************************************/
void cpedsMCMC::saveCurrentTemperature() {
	cpedsList<double> ct;
	ct.append(_cooling.temperatureHistory.last());
	ct.save(getPartialFilesDirFull()+"/current-temperature"); 			
}
/***************************************************************************************/
// currently not used
void cpedsMCMC::generateInitialStepSizePDFsAfterBurnIn() {
	mscsFunction f("PDF",Zero);
	
	for (long i = 0; i < _walk.stepGeneratingDistr.count(); i++) {		
		//
		// define step size generating distribution
		//
		f.mkBoltzmann(0,_cooling.initialMaximalEnergy*CPEDS_kB*getTemperature(),_cooling.initialMaximalEnergy*CPEDS_kB*_cooling.temperature/_cooling.stepGeneratingPDFPoints,getTemperature());
		f.scaleX(_walk.initialPDFstepCRfractionAfterBurnIn*(_parameterSpace.last().getMaxArg()-_parameterSpace.last().getMinArg())/f.getMaxArg()); // make the maximal step size equal to the requested fraction of the domain size 
		f.normalize(); // normalize PDF
		_walk.stepGeneratingDistr[i]=f;
		_walk.stepGeneratingDistr[i].setVerbosityLevel(Zero);
		
		// define step size generating RNG
		_walk.stepGeneratingRNG[i].setRNsType(_walk.stepGeneratingRNG[i].from_array_invCDF);
		_walk.stepGeneratingRNG[i].setPDF(f.pointsCount(),f.extractArguments(),f.extractValues());	
		f.clearFunction();
	}
	
}
/***************************************************************************************/
void cpedsMCMC::setAcceptWorsePvalues(double bipv, double abipv) {
	_cooling.acceptWorsePvalue=bipv;
	_cooling.acceptWorsePvalueAfterBurnIn=abipv;
	//	_acceptWorseThreshold=_coolingDistr.expectationValue();
	updateAcceptWorseThresholds();
	//	_maximalRejectionsCount=long(1.0/_acceptWorsePvalueAfterBurnIn);
	
}
/***************************************************************************************/
void cpedsMCMC::updateAcceptWorseThresholds() {
	_cooling.acceptWorseThreshold=_cooling.coolingDistrInvCDF.f(1.0-_cooling.acceptWorsePvalue,NULL);
	_cooling.acceptWorseThresholdAfterBurnIn=_cooling.coolingDistrInvCDF.f(1.0-_cooling.acceptWorsePvalueAfterBurnIn,NULL);
}
/***************************************************************************************/
void cpedsMCMC::acceptWorseState(bool acceptWorse, MClink& next, long& rejectionsCount) {
	if (acceptWorse) {
		next.setAccepted(true);
		_MCaccepted.append(next);
		appendLinkForConvergenceTest(next);
		current()=next;			
		msgs->say("ACCEPTING the step",Zero);
		rejectionsCount=0;
	}
	else {
		msgs->say("REJECTING the step",Zero);
		next.setAccepted(false);
		//				_MCchain.prepend(next); // store the calculations that were done as they also probe the likelihood surface, but we are prepending to avoid confusing the delta chisq calculations
		_walk.rejectedLinkNumber++;
		msgs->say("rejected "+msgs->toStr(rejectionsCount)+" points in row. Maximal rejections count is: "+msgs->toStr(_cooling.maximalRejectionsCount),Zero);
		_MCrej.append(next); // store the calculations that were done as they also probe the likelihood surface, but we are prepending to avoid confusing the delta chisq calculations
//		appendLinkForConvergenceTest(next);
		rejectionsCount++;		
	}
}
/***************************************************************************************/
void cpedsMCMC::setAcceptWorseUniformScheme(bool uniform) {
	_cooling.acceptWorseStateFromUniformDistr=uniform;
	updateCoolingPDF();
}
/***************************************************************************************/
void cpedsMCMC::setBurnInLength(long l) { 
	_walk.convergenceTestMinimalLength=l; _convergence.minimalLengthToTestConvergence=l; 
	if (l>_walk.maximalChainLength) setMaximalChainLength(l);
}
/* ******************************************************************************************** */
void cpedsMCMC::setBurnOutLength(long l) {
	_walk.burnOutLength=l;
	if (l>_walk.maximalChainLength) setMaximalChainLength(l);	
}
/***************************************************************************************/
void cpedsMCMC::saveBestFitCovSimulations() {
	msgs->warning("you are calling an abstract saveBestFitCovSimulations() method. I hope this is what you want.", High);
}
/***************************************************************************************/
void cpedsMCMC::saveBestFit() {
	if (_cooling.acceptWorseStateFromUniformDistr==false) {
		_chisqData.bestFitData.coolingDistr.save(getRunDependentFileName(getPartialFilesDirFull()+"/bestFit-coolingPDF")); 
	}
	saveBestFitStepPDFs(getPartialFilesDirFull()+"/bestFit-stepPDF"); 

	_chisqData.bestFitData.model.save(getRunDependentFileName(getPartialFilesDirFull()+"/bestFit-model")); 			
	//		_temperatureHistory.save(_partialFilesDir+"/temperature-history"); 			
	_bestFit.save(getRunDependentFileName(getPartialFilesDirFull()+"/bestFit-link"));
	
	
	//
	// check if the best fit is close to the imposed prior boundary and if so, issue a warrning
	//
	for (long i = 0; i < _bestFit.dims(); i++) {
		mscsFunction p=getParameterSpace()[i];
		double min=p.getX(0);
		double max=p.getX(1);
		double delta=max-min;
		double thres=0.05;
		if ((_bestFit[i] - min)/delta < thres or (max-_bestFit[i])/delta > 1.0-thres) {
			msgs->say("WARNING: the best-fit solution for parameter %li is close to the boundary, or beyond the range defined by the prior. "
					"Consider enlarging the parameter space volume to improve likelihood reconstruction.", i,High);
		}
	}
	
	
	saveChain(_MCbestFitChain,getRunDependentFileName(getPartialFilesDirFull()+"/bestFitChain"));
	
	if (_chisqData.covarianceMatrixDiag.size()>0 and _chisqData.covDiagonal) {
		_chisqData.bestFitData.covDiag.save(getRunDependentFileName(getPartialFilesDirFull()+"/bestFit-covarianceMatrixDiagonal.mc"),false,"double");
	}

	
/*
	// save length
	cpedsList<long> l;
	l.append(_chisqData.bestFitData.lengthAtStoring);
	l.save(getPartialFilesDirFull()+"/current-length",false,"long");
	
	//save temperature
	cpedsList<double> ct;
	ct.append(_chisqData.bestFitData.temperatureAtStoring);
	ct.save(getPartialFilesDirFull()+"/current-temperature"); 			
*/
	
	
}
/***************************************************************************************/
void cpedsMCMC::saveResults() {
	saveSetup();
	setDumpEveryMCMCstate(true);
	dumpAll();
	saveBestFit();
	saveResiduals();
	calculate1Dposteriors();
	calculate1DCRs();
	save1Dposteriors();
	if (_IOcontrol.calculate2Dpdf) { //TODO: implement this stuff analogically as 1D posteriors
		save2Dposteriors();
		for (long i = 0; i < chisqData().bestFitData.confidenceLevels.size(); i++) {
			save2DCR(chisqData().bestFitData.confidenceLevels[i]);			
		}
	}

	saveBestFitCovSimulations();
}
/***************************************************************************************/
void cpedsMCMC::calculate1DCRs() {
	_chisqData.bestFitData.confidenceRanges1D.clear();
	// calculate 1d CRs
	
	for (long i = 0; i < getNparam(); i++) {
		msgs->say("Calculating 1-D CRs for parameter: %li",i,Medium);
		
		mscsFunction3dregc CRs;
		CRs.setSize(4,chisqData().bestFitData.confidenceLevels.size());
		CRs.allocFunctionSpace();
		for (long j = 0; j < chisqData().bestFitData.confidenceLevels.size(); j++) {
			double CL=chisqData().bestFitData.confidenceLevels[j];
			double lvl;
			cpedsList<double> CR=get1DCR(i,CL,&lvl);
			
			CRs(0,j)=CL;
			CRs(1,j)=CR[0];
			CRs(2,j)=CR[CR.size()-1];
			CRs(3,j)=lvl;
			
		}
		chisqData().bestFitData.confidenceRanges1D.append(CRs);	
	}
}
/* ******************************************************************************************** */
void cpedsMCMC::calculate1Dposteriors() {
	_chisqData.bestFitData.posteriors1D.clear();
	msgs->say("Calculating 1-D posteriors",Medium);

	for (long i = 0; i < getNparam(); i++) {
		MscsPDF1D p=get1Dposterior(i,_IOcontrol.PDF1d_pointsCount);
		p.checkRanges();
		_chisqData.bestFitData.posteriors1D.append(p);
	}	
}
/***************************************************************************************/
void cpedsMCMC::calculate2Dposteriors() {
	if (_IOcontrol.calculate2Dpdf) {
		_chisqData.bestFitData.posteriors2D.clear();
		for (long i = 0; i < getNparam(); i++) {
			for (long j = i+1; j < getNparam(); j++) {
				msgs->say("Calculating 2-D PDF for parameters: %li %li",i,j,Medium);
				mscsFunction3dregc p=get2Dposterior(i,j,_IOcontrol.PDF2d_pointsCount);
				_chisqData.bestFitData.posteriors2D.append(p);
			}
		}		
	}
}
/***************************************************************************************/
void cpedsMCMC::save1Dposteriors() {
	if (_chisqData.bestFitData.posteriors1D.size()!=dims()) {
		calculate1Dposteriors();
	}
	if (_chisqData.bestFitData.confidenceRanges1D.size()!=dims()) {
		calculate1DCRs();		
	}

	for (long i = 0; i < getNparam(); i++) {
		_chisqData.bestFitData.posteriors1D[i].saveHDF5(getOutputDir()+"/posterior1D.hdf5",msgs->toStr(i));

		// save parameter name
		_chisqData.bestFitData.posteriors1D[i].setHDF5_scalarStringAttribute(
				getOutputDir()+"/posterior1D.hdf5",msgs->toStr(i),
				"full_name",parameterSpace().names_latex()[i]);

		//
		//  save best-fit parameter value
		//
		_chisqData.bestFitData.posteriors1D[i].setHDF5_scalarDoubleAttribute(
				getOutputDir()+"/posterior1D.hdf5",msgs->toStr(i),
				"bestFit",bestFitLink().getParam(i),"best fit parameter value");

		// save CR
		_chisqData.bestFitData.confidenceRanges1D[i].saveHDF5(
				getOutputDir()+"/posterior1D.hdf5",msgs->toStr(i)+"_CR");
	
	}
	
}
/***************************************************************************************/
void cpedsMCMC::save2Dposteriors() {
	if (_chisqData.bestFitData.posteriors2D.size()!=double(dims()-1)/2*dims()) calculate2Dposteriors();
//	long k=0;
	for (long i = 0; i < getNparam(); i++) {
		for (long j = i+1; j < getNparam(); j++) {
			
//			printf("i: %li j:%li, param idx: %li, names_latex_size: %i\n",i,j,paramij2idx(i,j),parameterSpace().names_latex().size());
			_chisqData.bestFitData.posteriors2D[paramij2idx(i,j)].saveHDF5(
					getOutputDir()+"/posterior2D.hdf5",msgs->toStr(i)+"-"+msgs->toStr(j));
			
			// save parameter names
			_chisqData.bestFitData.posteriors2D[paramij2idx(i,j)].setHDF5_scalarStringAttribute(
					getOutputDir()+"/posterior2D.hdf5",msgs->toStr(i)+"-"+msgs->toStr(j),
					"full_name_x",parameterSpace().names_latex()[i]);
			
			_chisqData.bestFitData.posteriors2D[paramij2idx(i,j)].setHDF5_scalarStringAttribute(
					getOutputDir()+"/posterior2D.hdf5",msgs->toStr(i)+"-"+msgs->toStr(j),
					"full_name_y",parameterSpace().names_latex()[j]);
			
			//
			//  save best-fit parameter
			//
			_chisqData.bestFitData.posteriors2D[paramij2idx(i,j)].setHDF5_scalarDoubleAttribute(
					getOutputDir()+"/posterior2D.hdf5",msgs->toStr(i)+"-"+msgs->toStr(j),
					"bestFit_x",bestFitLink().getParam(i),"best fit parameter for X axis");
			
			_chisqData.bestFitData.posteriors2D[paramij2idx(i,j)].setHDF5_scalarDoubleAttribute(
					getOutputDir()+"/posterior2D.hdf5",msgs->toStr(i)+"-"+msgs->toStr(j),
					"bestFit_y",bestFitLink().getParam(j),"best fit parameter for Y axis");

//			exit(0);
		}

	}
	
	
	
}
/* ******************************************************************************************** */
void cpedsMCMC::recalculateLikelihoods() {
	cpedsList<double> p;
	for (long i = 0; i < acceptedStates().size(); i++) {
		p.append(chain(i).chisq());
	}
	
	for (long i = 0; i < rejectedStates().size(); i++) {
		p.append(rejectedStates()[i].chisq());
	}
	
	double chisqOffset,dummy;
	long i,j;
	
	p.getMinMaxValues(&chisqOffset,&dummy,&i,&j);

	for (long i = 0; i < acceptedStates().size(); i++) {
		chain(i).setL(exp(-(chain(i).chisq()-chisqOffset)/2));
	}
	
	for (long i = 0; i < rejectedStates().size(); i++) {
		rejectedStates()[i].setL(exp(-(rejectedStates()[i].chisq()-chisqOffset)/2));
	}

}
/* ******************************************************************************************** */
long cpedsMCMC::statesAboveRejectionThreshold() {
	cpedsList<double> p;
	for (long i = 0; i < acceptedStates().size(); i++) {
		p.append(chain(i).chisq());
	}
	
	for (long i = 0; i < rejectedStates().size(); i++) {
		p.append(rejectedStates()[i].chisq());
	}
	
	double chisqOffset,dummy;
	long i,j;
	
	p.getMinMaxValues(&chisqOffset,&dummy,&i,&j);
	long Nstates=0;
	
	for (long i = 0; i < acceptedStates().size(); i++) {
		if (exp(-(chain(i).chisq()-chisqOffset)/2)> getWalk().likelihoodRejectionThreshold) Nstates++;
	}
	
	for (long i = 0; i < rejectedStates().size(); i++) {
		if (exp(-(rejectedStates()[i].chisq()-chisqOffset)/2) > getWalk().likelihoodRejectionThreshold) Nstates++;
	}
	
	return Nstates;
}
/***************************************************************************************/
/*
void cpedsMCMC::save1DCR(double CL, string fname, string dset) {
	msgs->say("Calculating confidence regions for CL=%lf %%",CL*100,Medium);

	for (long i = 0; i < getNparam(); i++) {
		double lvl;
		stringstream ssCL;
		ssCL << "CL " << setprecision(2) << CL*100;
		

		// save 1D CR contours for given CL
		cpedsList<double> CR=get1DCR(i,CL,&lvl);

		stringstream ssCRlow,ssCRhigh;
		ssCRlow  << "CRlow " << CR[0];
		ssCRhigh << "CRhigh " << CR[0];

		// save CR for given CL
		_chisqData.bestFitData.posteriors1D[i].setHDF5_scalarDoubleAttribute(
				fname,dset,
				ss.str(),lvl,"PDF value corresponding to a confidence region");
		
		// save CR level for given CL
		_chisqData.bestFitData.posteriors1D[i].setHDF5_scalarDoubleAttribute(
				fname,dset,
				ss.str(),lvl,"PDF value corresponding to a confidence region");
		
	}

}
*/
/* ******************************************************************************************** */

void cpedsMCMC::save2DCR(double CL) {
	msgs->say("Calculating confidence regions for CL=%lf %%",CL*100,Medium);
	for (long i = 0; i < getNparam(); i++) {
		for (long j = i+1; j < getNparam(); j++) {
			double lvl;
//			cpedsList<double> lvls;
			stringstream ss;
			ss << "CLcontour" << setprecision(2) << CL*100;

			// save 2D CR contours for given CL
			get2DCR(i,j,CL,&lvl).saveHDF5(getOutputDir()+"/CR_"+msgs->toStrf(CL*100,0)+".hdf5",msgs->toStr(i)+"-"+msgs->toStr(j));
//			get2DCR(i,j,CL,&lvl).save(getOutputDir()+"/CR"+msgs->toStrf(CL,2)+".params"+msgs->toStr(i)+"_"+msgs->toStr(j));

			// save 
//			lvls.append(lvl);
//			lvls.save(getOutputDir()+"/PDF-contour"+msgs->toStrf(CL,2)+".params"+msgs->toStr(i)+"_"+msgs->toStr(j));

			// save 2-D contour level for given CL
			_chisqData.bestFitData.posteriors2D[paramij2idx(i,j)].setHDF5_scalarDoubleAttribute(
					getOutputDir()+"/posterior2D.hdf5",msgs->toStr(i)+"-"+msgs->toStr(j),
					ss.str(),lvl,"contour level corresponding to a confidence region");

			
//			cout << ss.str() << "\n";
		}
	}

}
/***************************************************************************************/
cpedsList<double> cpedsMCMC::get1DCR(int paramID1, double CL, double* LVL) {
	MscsPDF1D pdf=get1Dposterior(paramID1);
	cpedsList<double> CR;
	CR=pdf.getCR(CL,LVL);
	
	return CR;	
}
/* ******************************************************************************************** */
mscsFunction cpedsMCMC::get2DCR(int paramID1, int paramID2, double CL, double* LVL) {
	MscsPDF2D pdf=get2Dposterior(paramID1, paramID2);
	mscsFunction CR;
	CR=pdf.getContour(CL,LVL);
	
	return CR;
}
/***************************************************************************************/
void cpedsMCMC::storeBestFit() {
	_chisqData.bestFitData.lengthAtStoring=_MCaccepted.length();
	_chisqData.bestFitData.coolingDistr=_cooling.coolingDistr;
	_chisqData.bestFitData.stepGeneratingDistr=_walk.stepGeneratingDistr;
	_chisqData.bestFitData.model=_chisqData.model;
	_chisqData.bestFitData.model.setVerbosityLevel(msgs->getVerbosity());
	_chisqData.bestFitData.temperatureAtStoring=_cooling.temperatureHistory.last();
	
	if (_chisqData.covarianceMatrixDiag.size()>0 and _chisqData.covDiagonal) {
		_chisqData.bestFitData.covDiag=_chisqData.covarianceMatrixDiag;
	}	
}
/***************************************************************************************/
MscsPDF1D cpedsMCMC::get1Dposterior(int paramID, long pdfPoints, string interpolationType) {
	
	if (paramID<_chisqData.bestFitData.posteriors1D.size()) { // if it was calculated then return it right away
		return _chisqData.bestFitData.posteriors1D[paramID];
	}

	msgs->say("Calculating 1-D PDF for parameter: %li",paramID,Medium);

	mscsFunction states=getParamValues(paramID);
	mscsFunction statesInv=states;
	statesInv.invert();
	statesInv.sortFunctionArgAscending();
	statesInv.invert();
	double increment=fabs(statesInv.getX(0)-statesInv.getX(1));
	
	states.checkRanges();
	states/=states.getMaxValue();
#ifdef DEBUG_MCMC_PDF
	states.save("raw-states-paramID"+msgs->toStr(paramID,"%i")+".L");
#endif
	states.removeSmaller(1e-10);
	states.sortFunctionArgAscending();
	states.checkRanges();
#ifdef DEBUG_MCMC_PDF
	cout << "Number of states after cutting off at L<1e-10: "<< states.pointsCount() <<"\n";
#endif
	
	MscsPDF1D pdf;
	
	
	cpedsList<long> bins;
	double pdfRange;//=sqrt(params.variance());
	if (states.pointsCount()<2) pdfRange=increment; else pdfRange=states.getMaxArg()-states.getMinArg();
	long pdfPoints_internal=long(sqrt(states.pointsCount()));
	double dx=pdfRange/pdfPoints_internal;
	
	
	bins.clear();
	pdf=states.binFunction(dx,bins,"bin_center","max");
	pdf.checkRanges();
#ifdef DEBUG_MCMC_PDF
	pdf.save("rusty-paramID"+msgs->toStr(paramID,"%i")+".pdf");
	states.save("states-paramID"+msgs->toStr(paramID,"%i")+".pdf");
#endif
	
	if (pdf.pointsCount()<6) { // we do not have enough points to probe the pdf, so we need to use some approximations

		
		// we will use gaussian PDF with fwhm = dx
		double x0=pdf.getX(pdf.getMaxValueIdx());
		printf("the number of PDF points for parameter %li is < 6, assuming gaussian PDF around %lE with fwhm=%lE\n",paramID,x0,dx);
		pdf.clearFunction();
		pdf.mkGauss(x0-pdfPoints/2*dx,x0+pdfPoints/2*dx,dx,1.0,x0,cpeds_fwhm2sigma(dx),0.0);
#ifdef DEBUG_MCMC_PDF
		pdf.save("assumed-paramID"+msgs->toStr(paramID,"%i")+".pdf");
//		exit(0);
#endif
		
//		assert(pdf.pointsCount()>5);
	}
	else {
		//	pdf.interpolate(pdf.getMinArg(),pdf.getMaxArg(),pdfPoints,true,"akima");
			pdf.interpolate(pdf.getMinArg(),pdf.getMaxArg(),pdfPoints,true,interpolationType);
			pdf/=pdf.getMaxValue();

#ifdef DEBUG_MCMC_PDF
			pdf.save("rusty-interpolated-paramID"+msgs->toStr(paramID,"%i")+".pdf");
#endif

			cpedsList<double> CR=pdf.getCR(0.9973).sort(12);
#ifdef DEBUG_MCMC_PDF
			CR.save("CR-paramID"+msgs->toStr(paramID,"%i")+".tmp");
#endif
			assert(CR.size()>1);
			if (CR.size()>2) { cout << "WARNING 0.9973 confidence interval suspicious. Will return a rusty pdf.\n"; }
			else {
				double margin=CR[CR.size()-1]-CR[0];
				pdf=pdf.cut(CR[0]-margin/2,CR[CR.size()-1]+margin/2);
#ifdef DEBUG_MCMC_PDF
				pdf.save("rusty-interpolated-cut-paramID"+msgs->toStr(paramID,"%i")+".pdf");
#endif
			}
	}
	return pdf;	

	
/*
	mscsFunction pdf;
	
	states.sortFunctionArgAscending();
	states.checkRanges();
	
	cpedsList<long> bins;
	double bin_factor=1;
	double pdfRange;//=sqrt(params.variance());
//	pdfRange=states.getX(states.pointsCount()*5/8)-states.getX(states.pointsCount()*3/8);
	pdfRange=states.getMaxArg()-states.getMinArg();
	
	printf("pdfRange: %f\n",pdfRange);
	long pdfPoints_internal=1000;
	double dx=pdfRange/pdfPoints_internal;
	printf("dx: %f\n",dx);
	
	
//	do {
		bins.clear();
		pdf=states.binFunction(dx,bins,"bin_center","max");
//		pdf=states.binFunction(dx*bin_factor,bins,"bin_center","max");
//		pdf.save("pdf");
//		bins.save("bins");
//		pdf/=bins;
		pdf.checkRanges();
		pdf/=pdf.getMaxValue();
		pdf.removeSmaller(1e-10);
		pdf.checkRanges();
		pdfRange=pdf.getMaxArg()-pdf.getMinArg();
		dx=pdfRange/pdfPoints_internal;
//		if (pdf.pointsCount()<pdfPoints_internal) {
//			msgs->say("adjusting PDF resolution (pdf points: %.0lf, bin_factor: %lf)",double(pdf.pointsCount()),bin_factor,Medium);
//			bin_factor/=1.04; 
//		}
//		else bin_factor*=1.04;
//	} while (pdf.pointsCount()<pdfPoints_internal);

	pdf=states.binFunction(dx,bins,"bin_center","max");

	pdf.removeSmaller(1e-10);
	pdf.checkRanges();
	pdf.interpolate(pdf.getMinArg(),pdf.getMaxArg(),pdfPoints,true,"cspline");
	pdf/=pdf.getMaxValue();
	
	return pdf;	
*/
}
/***************************************************************************************/
MscsPDF1D cpedsMCMC::get1Dposterior(string paramName, long pdfPoints, string interpolationType) {
	return get1Dposterior( getParameterByName(paramName),pdfPoints, interpolationType );	
}
/***************************************************************************************/
mscsFunction3dregc cpedsMCMC::get2Dposterior(string paramName1, string paramName2, long pdfPoints) {
	return get2Dposterior( getParameterByName(paramName1), getParameterByName(paramName2),pdfPoints );	
}
/***************************************************************************************/
mscsFunction3dregc cpedsMCMC::get2Dposterior(int paramID1, int paramID2, long pdfPoints) {
	
	if (paramij2idx(paramID1,paramID2)<_chisqData.bestFitData.posteriors2D.size()) { // if it was calculated then return it right away
		return _chisqData.bestFitData.posteriors2D[paramij2idx(paramID1,paramID2)];
	}
	
	
/*
	mscsFunction states1=getParamValues(paramID1);
	mscsFunction states2=getParamValues(paramID2);
	mscsFunction3dregc pdf, PDF;
	if (_walk.Nparam<2) return pdf;
	
	double p1Min,p1Max,p2Min,p2Max;
	if (_bestFitData.posteriors1D.size()!=dims()) calculate1Dposteriors();
	p1Min=posteriors()[paramID1].getMinArg();
	p1Max=posteriors()[paramID1].getMaxArg();
	p2Min=posteriors()[paramID2].getMinArg();
	p2Max=posteriors()[paramID2].getMaxArg();
	long pdfPoints_internal=500;
	long Nx=pdfPoints_internal;
	long Ny=pdfPoints_internal;
	cpedsPointSet3D ps(states1.pointsCount(),states1.toXList().toCarray(),states2.toXList().toCarray(),states1.toYList().toCarray(),states1.toYList().toCarray(),true);
	pdf.setSizeRange(Nx,Ny,1,p1Min,p2Min,0,p1Max,p2Max,0);
	pdf.allocFunctionSpace();
	pdf.populateField(ps,true,0,0,false,true,"max");
//	
	Nx=pdfPoints;
	Ny=pdfPoints;
	PDF.setSizeRange(Nx,Ny,1,p1Min,p2Min,0,p1Max,p2Max,0);
	PDF.allocFunctionSpace();
	for (long i = 0; i < Nx; i++) {
		for (long j = 0; j < Ny; j++) {
			PDF.setf(i,j,0, pdf.fxy(PDF.getX(i),PDF.getY(j)),0);
		}
	}

	PDF/=PDF.getMaxValue();
	
	return PDF;		
*/

	

	
	

	
	mscsFunction states1=getParamValues(paramID1);
	mscsFunction states2=getParamValues(paramID2);
	mscsFunction3dregc pdf;
	if (getNparam()<2) return pdf;

	
#ifdef DEBUG_MCMC_PDF2
	states1.save(getOutputDir()+"/"+getParameterDependentFileName("debug-states1",paramID1));
	states2.save(getOutputDir()+"/"+getParameterDependentFileName("debug-states2",paramID2));
#endif

	double p1Min,p1Max,p2Min,p2Max;
	if (_chisqData.bestFitData.posteriors1D.size()!=dims()) calculate1Dposteriors();
	cpedsList<double> CR;
	double margin;
	CR=posteriors()[paramID1].getCR(0.9973);
//	CR.print();
	CR.sort(12);
//	CR.print();
//	printf("p%i CR: %lf %lf\n",paramID1, CR[0],CR[CR.size()-1]);
	margin=fabs(CR[0]-CR[CR.size()-1]);
	p1Min=posteriors()[paramID1].getMinArg();
	p1Max=posteriors()[paramID1].getMaxArg();
//	posteriors()[paramID1].save("p1.1d");
	p1Min=CR[0]-margin/2;
	p1Max=CR[CR.size()-1]+margin/2;
//	printf("p1min: %lf p1max: %lf\n",p1Min,p1Max);
	CR=posteriors()[paramID2].getCR(0.9973);
//	CR.print();
	CR.sort(12);
//	CR.print();
//	posteriors()[paramID2].save("p2.1d");
//	printf("p%i CR: %lf %lf\n",paramID2, CR[0],CR[CR.size()-1]);
	margin=fabs(CR[0]-CR[CR.size()-1]);
	p2Min=posteriors()[paramID2].getMinArg();
	p2Max=posteriors()[paramID2].getMaxArg();
//	printf("p2min: %lf p2max: %lf\n",p2Min,p2Max);
	p2Min=CR[0]-margin/2;
	p2Max=CR[CR.size()-1]+margin/2;
#ifdef DEBUG_MCMC_PDF2
	printf("p1min: %lf p1max: %lf\n",p1Min,p1Max);
	printf("p2min: %lf p2max: %lf\n",p2Min,p2Max);
#endif
	long Nx=pdfPoints;
	long Ny=pdfPoints;
	cpedsPointSet3D ps(states1.pointsCount(),states1.toXList().toCarray(),states2.toXList().toCarray(),states1.toYList().toCarray(),states1.toYList().toCarray(),true);

	
	// maximize through parameter space (this pdf may have holes in it if not enough MCchains were run)
	pdf.setSizeRange(Nx,Ny,1,p1Min,p2Min,0,p1Max,p2Max,0);
	pdf.allocFunctionSpace();
	pdf.populateField(ps,true,0,0,false,true,"max");
	pdf/=pdf.getMaxValue();
#ifdef DEBUG_MCMC_PDF
	pdf.saveSlice(2,0,"2dposterior.slice");
	ps.save("2dposterior.ps");	
	printf("p1Min: %lE, p1Max: %lE, p2Min: %lE, p2Max: %lE\n",p1Min,p1Max,p2Min,p2Max);
#endif
	
	// create a new set of points for interpolation
	ps.clear();
	mscsVector<double> vals;
	for (long i = 0; i < Nx; i++) {
		for (long j = 0; j < Ny; j++) {
			if (pdf.fRe(i,j,0)>1e-10) {
				cpedsPoint3D p(pdf.getX(i),pdf.getY(j),pdf.fRe(i,j,0));
				p.setZ(0); // this is for the case when 2d PDF is interpolated using SP-interpolation
				ps.append(p);
				vals.push_back(pdf.fRe(i,j,0));
			}
		}
	}
	assert(ps.size()>0);
	
	subDomain_region_t r,t;
	r.subx=pdf.Nx();
	r.suby=pdf.Ny();
	r.subz=1;
	r.xmin=pdf.getMinX();
	r.xmax=pdf.getMaxX();
	r.ymin=pdf.getMinY();
	r.ymax=pdf.getMaxY();
	r.zmin=0;
	r.zmax=0;
	t=r;
	t.subx=2;
	t.suby=2;
	t.subz=1;

//	printf("r xmin: %lE, xmax: %lE \n",r.xmin,r.xmax);
//	printf("r ymin: %lE, ymax: %lE \n",r.ymin,r.ymax);
//	printf("r zmin: %lE, zmax: %lE \n",r.zmin,r.zmax);
	
	mscsVector<cpedsPoint3D> v=ps.exportAsVector();

#ifdef DEBUG_MCMC_PDF2
	pdf.setVerbosityLevel(Top);
#endif
	pdf.mkInterpolatedFieldScatter(r,t,v,vals,"gadget2");
//	ps.save("ps");
//	pdf.saveSlice(2,0,"pdf2d");
//	printf("%li\n",ps.size());
//	pdf.mkInterpolatedFieldTriangLinear2D(ps);

#ifdef DEBUG_MCMC_PDF2
	pdf/=pdf.getMaxValue();
	pdf.saveSlice(2,0,"pdfint2d");
	pdf.printInfo();
//	exit(0);
#endif


	
	
	
	
	
	
	
	
/*
	
	mscsFunction states1=getParamValues(paramID1);
	mscsFunction states2=getParamValues(paramID2);
	mscsFunction3dregc pdf;
	if (_walk.Nparam<2) return pdf;


	double p1Min,p1Max,p2Min,p2Max;
	if (_bestFitData.posteriors1D.size()!=dims()) calculate1Dposteriors();
	p1Min=posteriors()[paramID1].getMinArg();
	p1Max=posteriors()[paramID1].getMaxArg();
	p2Min=posteriors()[paramID2].getMinArg();
	p2Max=posteriors()[paramID2].getMaxArg();
	long Nx=pdfPoints;
	long Ny=pdfPoints;
	cpedsPointSet3D ps(states1.pointsCount(),states1.toXList().toCarray(),states2.toXList().toCarray(),states1.toYList().toCarray(),states1.toYList().toCarray(),true);
	pdf.setSizeRange(Nx,Ny,1,p1Min,p2Min,0,p1Max,p2Max,0);
	pdf.allocFunctionSpace();
	pdf.populateField(ps,true,0,0,false,true,"max");
	pdf/=pdf.getMaxValue();
	
	ps.clear();
	for (long i = 0; i < Nx; i++) {
		for (long j = 0; j < Ny; j++) {
			if (pdf.fRe(i,j,0)>1e-10) {
				cpedsPoint3D p(pdf.getX(i),pdf.getY(j),0);
				ps.append(p,pdf.fRe(i,j,0));
			}
		}
	}

	subDomain_region_t r,t;
	r.xmin=p1Min; 		r.xmax=p1Max;
	r.ymin=p2Min; 		r.ymax=p2Max;
	r.zmin=0; 			r.zmax=0;
	r.subx=Nx;
	r.suby=Ny;
	r.subz=1;
	t=r;
	t.subx=2;
	t.suby=2;
	t.subz=1;
	
	pdf.mkInterpolatedFieldScatter(r,t,ps.exportAsVector(),ps.values(),"gadget2",200,250);
	
	pdf/=pdf.getMaxValue();
*/
	
	return pdf;		
}
/***************************************************************************************/
void cpedsMCMC::append(MClink& link) {
	_MCaccepted.append(link);
}
/***************************************************************************************/
void cpedsMCMC::append(cpedsMCMC& chain) {
	setData(chain.getData());
	_MCaccepted.append(chain.acceptedStates());
	_MCrej.append(chain.rejectedStates());
	_MCbestFitChain.append(chain._MCbestFitChain);
	if (chain.bestFitLink().chisq() < bestFitLink().chisq() or bestFitLink().chisq()==-1) {
		bestFitLink()=chain.bestFitLink();
		_chisqData.bestFitData=chain.getBestFitData();
	}
	
	// append parts of structures 
	for (long i = 0; i < getNparam(); i++) {
		_walk.steps[i].append(chain.getMCsteps(i));
	}
}
/***************************************************************************************/
void cpedsMCMC::setup(cpedsMCMC& chain) {
	_cooling=chain.getCooling();
	_parameterSpace=chain.getParameterSpace();
	_walk=chain.getWalk();
	_chisqData=chain.getChisqData();
	_IOcontrol=chain.getIOcontrol();
	_convergence=chain.getConvergence();
}
/***************************************************************************************/
mscsFunction cpedsMCMC::getParamValues(int j) {
	mscsFunction p;
	
/*
	for (long i = 0; i < chain().size(); i++) {
		p.newPoint(chain(i).getParam(j),chain(i).L());
	}
	
	for (long i = 0; i < rejectedStates().size(); i++) {
		p.newPoint(rejectedStates()[i].getParam(j),rejectedStates()[i].L());
	}
*/

	for (long i = 0; i < acceptedStates().size(); i++) {
		p.newPoint(chain(i).getParam(j),chain(i).chisq());
	}
	
	for (long i = 0; i < rejectedStates().size(); i++) {
		p.newPoint(rejectedStates()[i].getParam(j),rejectedStates()[i].chisq());
	}

	p.checkRanges();
	p-=p.getMinValue(); // this effectively scales the PDF by a constant factor
	
#ifdef DEBUG_MCMC_PDF
	p.save("raw-states-paramID"+msgs->toStr(j,"%i")+".chisq");
#endif

	for (long i = 0; i < p.pointsCount(); i++) {
		p.setf(i,exp(-p.f(i)/2));
	}

	return p;
}
/***************************************************************************************/
int cpedsMCMC::paramij2idx(int i,int j) {
	assert(i<j);
	//TODO: make this smarter
	int idx=0;
	for (long ii = 0; ii < getNparam(); ii++) {
		for (long jj = ii+1; jj < getNparam(); jj++) {
			if (ii==i and jj==j) return idx;
			idx++;
		}
	}
	idx=-1;
	assert(idx>=0);
	return idx;
}
/***************************************************************************************/
void cpedsMCMC::setVerbocity(cpeds_VerbosityLevel verb) {
	msgs->setVerbosity(verb);
}
/***************************************************************************************/
/***************************************************************************************/
cpedsMCMC& cpedsMCMC::operator=(cpedsMCMC& rhs) {
	_MCaccepted=rhs.acceptedStates();
	_MCrej=rhs.rejectedStates();
	_bestFit=rhs.bestFitLink();
	_MCbestFitChain=rhs.BFchain();
	_startingLink=rhs.getStartingLink();
	_cooling=rhs.getCooling();
	_parameterSpace=rhs.getParameterSpace();
	_walk=rhs.getWalk();
	_chisqData=rhs.getChisqData();
	_IOcontrol=rhs.getIOcontrol();
	_convergence=rhs.getConvergence();

	return *this;	
}
/***************************************************************************************/
void cpedsMCMC::setID(long id) {
//	msgs->setSender(msgs->getSender()+"."+msgs->toStr(id));
//	id+=_cooling.coolingRNG.seed();
	
	_cooling.coolingRNG.seed(id);
	for (long i = 0; i < _walk.stepGeneratingRNG.size(); i++) { _walk.stepGeneratingRNG[i].seed(id); }
	for (long i = 0; i < _walk.rns.size(); i++) { _walk.rns[i].seed(id); }
	_walk.rnSgn.seed(id); 
	_chisqData.rnCov.seed(id);
	
	setRunIdx(id);
}
/***************************************************************************************/
void cpedsMCMC::coolForcibly() {
	if (_cooling.forcedCoolingTemperatureDecrement==-1) _cooling.forcedCoolingTemperatureDecrement=(getTemperature()-getFinalTemperature())/_cooling.maximalNumberOfForcedCoolings; 
//	if (_cooling.forcedCoolingTemperatureDecrement< 0.001*(_cooling.initialTemperature-_cooling.finalTemperature)) _cooling.forcedCoolingTemperatureDecrement=0.001*(_cooling.initialTemperature-_cooling.finalTemperature);
//	_cooling.temperature-=_cooling.forcedCoolingTemperatureDecrement;
	_cooling.forcedCoolingTemperatureDecrement=coolDownGeometric();
	msgs->say("Forcibly cooling by %f K",_cooling.forcedCoolingTemperatureDecrement,Medium);

	
	if (getTemperature()<getFinalTemperature()) _cooling.temperature=getFinalTemperature();
	_cooling.forcedCoolingTemperatureDecrement=-1;
	
	updateCoolingPDF();
	updateStepGeneratingPDFs();
	_walk.rejectionsCount=0;
}
/* ******************************************************************************************** */
void cpedsMCMC::appendLinkForConvergenceTest(MClink& link) {
	_convergence.chisq.append(link.chisq());
//	printf("appending link for convergence test: %lE\n",link.chisq());
	if (_convergence.chisq.size()>_convergence.statesToConvergenceCheck) _convergence.chisq.removeFirst();
}
/* ******************************************************************************************** */
void cpedsMCMC::addNewBestFitLink(MClink& l) {
	_bestFit=l; 
	_MCbestFitChain.append(_bestFit); 
	_IOcontrol.storeBestFitNow=true; 
	storeBestFit();
	if (_walk.uphillClimbing) _walk.keepDirectionNow=true; 

}
/* ******************************************************************************************** */
void cpedsMCMC::burnOut() {
	msgs->say("Starting burn-out states",Medium);
	setUpHillClimbing(false); // we move randomly from now on around the best-fit solution
	
	// update walk PDFs with expected step sizes fixed to well cover the high-likelihood region about the best fit solution
	calculate1Dposteriors();

	setInitialWalkStepSize(1); // we're not using this during burn-out, but I keep it here to remember of its existence

	for (long i = 0; i < dims(); i++) { // we need to find out the best step size in each direction for the burn-out
		double oneSigma=0;
		MscsPDF1D pdf=get1Dposterior(i);
//		pdf.save("burnout-pdf-param."+msgs->toStr(i));
		cpedsList<double> ran=pdf.getCR(0.68);
//		printf("ran param: %li\n",i);
//		ran.print();
		if (ran.size()>=2) {
			oneSigma=ran.last()-ran[0]; // 68% CR
		}
		else { 
			if (ran.size()==1) {
				long maxi=pdf.getMaxValueIdx();
				long st=maxi;
				oneSigma=fabs(pdf.getX(st)-ran[0]); // 68% CR				
			}
			else {
				// pdf is poorly probed after the walk, we will take guesses now.
				msgs->say("PDF for parameter %li is poorly probed after the walk, this is bad",i,Medium);
				
				long maxi=pdf.getMaxValueIdx();
				long st=maxi-1;
				long en=maxi+1;
				if (st<0) maxi=0;
				if (en>pdf.pointsCount()-1) en=pdf.pointsCount()-1;
				assert(en-st>0);
				oneSigma=pdf.getX(en)-pdf.getX(st); // 68% CR
			}
		}
		
		// expected stepPDF size
		double EXstep=oneSigma;
		msgs->say("Setting step PDF for parameter %.0lf to %lE",double(i),EXstep,Medium);
		
		updateStepGeneratingPDF(i,EXstep);
		
	}
	
	saveStepPDFs(getPartialFilesDirFull()+"/stepPDFburnOut");
	walkNstates(getBurnOutLength(),bestFitLink(),1);

}
/* ******************************************************************************************** */
void cpedsMCMC::updateStepGeneratingPDF(long param, double expectedStepSize) {
	msgs->say("Making step generating PDF for parameter: %li", param,Low);
	mscsFunction f("step generating PDF",Zero);
	double frac=1;

/*
	frac=_walk.initialPDFstepCRfraction;
	if (length()>=getBurnInLength()) {
		if (_walk.initialPDFstepCRfractionAfterBurnIn>0) frac=_walk.initialPDFstepCRfractionAfterBurnIn;
	}
*/

//	long i=param;
	
	f.setName(_parameterSpace[param].getName()+"PDF");
	f.mkBoltzmann(0,_cooling.initialMaximalEnergy*expectedStepSize,
			_cooling.initialMaximalEnergy*expectedStepSize/_cooling.stepGeneratingPDFPoints,
			expectedStepSize/CPEDS_kB); // we must provide temperature parameter in k_B units

	
	f.normalize(); // normalize
	_walk.stepGeneratingDistr[param]=f;
	
	
	// define step generating RNG
	_walk.stepGeneratingRNG[param].setPDF(f.pointsCount(),f.extractArguments(),f.extractValues());	


	
}
/* ******************************************************************************************** */
void cpedsMCMC::saveResiduals(int inputDataColumn) {
	string outfile=getOutputDir()+"/inputData.fit";
//		outfile=getRunDependentFileName(outfile);
	
	if (_chisqData.data2dCol>-1) { // this means that the data2d were set
		long bfCol,residCol;
		bfCol=_chisqData.data2d.Nx()-2;
		residCol=_chisqData.data2d.Nx()-1;
		
//		msgs->say("Saving residuals to test file: "+outfile+".test",Medium);
//		_chisqData.data2d.saveSlice(2,0,outfile,0);
//		_chisqData.data2d.setVerbosityLevel(Top);
//		_chisqData.data2d.printInfo();
//		printf("bf size: %li\n",_chisqData.bestFitData.model.pointsCount());
//		printf("data2dCol: %li\n",_chisqData.data2dCol);
//		exit(0);
		// calculate residuals
		for (long j = 0; j < _chisqData.data2d.Ny(); j++) {
			_chisqData.data2d(bfCol,j)=_chisqData.bestFitData.model.f(j);
			_chisqData.data2d(residCol,j)=_chisqData.data2d(_chisqData.data2dCol,j)-_chisqData.data2d(bfCol,j);
		}
		msgs->say("Saving residuals to file: "+outfile,Medium);
		_chisqData.data2d.saveSlice(2,0,outfile,0);
	}
	else {
		msgs->say("Residuals will not be saved.",Medium);
	}
}

/* ******************************************************************************************** */
MClink cpedsMCMC::getNextPointUphillGradient(const MClink& current) {
	msgs->say("generating the next point", Zero);
#ifdef DEBUGMCMC
	msgs->say("current state is:", Low);
	//	if (msgs->getVerbosity()>
	printLink(current);
	msgs->say("Best fit link is:", Low);
	printLink(_bestFit);	
#endif
	double grad[dims()];
	double delta[dims()];
	double sign;
	double tmp;
	long failures=0;
	MClink newLink(current),testLink;
//	printf("chisq cur: %lE\n",current.chisq());
	for (long i = 0; i < dims(); i++) {
		failures++;
//		delta[i]=_walk.initialPDFstepCRfractionAfterBurnIn*_parameterSpace.resolution(i)*getTemperature()/getInitialTemperature();
		delta[i]=_walk.initialPDFstepCRfractionAfterBurnIn*_parameterSpace.length(i)*getTemperature()/getInitialTemperature();
		testLink=current;
		testLink[i]+=delta[i];
		assert(delta[i]>0);
		// calculate gradient
		double X2=chisq(testLink);
		grad[i]=(X2-current.chisq())/current.chisq();
//		printf("i: %li, delta: %lE, chisq: %lE, grad[i]: %lE\n", i, delta[i],X2,grad[i]);
	}
	_walk.keepDirectionNow=false;

	// set new link
	for (long i = 0; i < dims(); i++) {
		newLink[i]-=grad[i]*delta[i];
		MCsteps(i).append(grad[i]*delta[i]); // store the step
	}
//	printf("new link\n");
//	newLink.printLink();
//	printf("\n\n");
	return newLink;		
}
