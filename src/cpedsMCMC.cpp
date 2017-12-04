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
	_IOcontrol.partialFilesDir=getOutputDir()+"/partial";
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
	
	_IOcontrol.PDF1d_pointsCount=500;
	_IOcontrol.PDF2d_pointsCount=500;
	_IOcontrol.calculate2Dpdf=false;
	
	//
	// walk parameters for step PDFs and RNGs
	// 
	_cooling.stepGeneratingPDFPoints=1000;
	_walk.runRNGseed=runSeed; // all other rngs will use this seed to initialize themselves (with different offsets of course)
	_walk.rnSgn.seed(_walk.runRNGseed); 
	_walk.rnSgn.setRNsType(_walk.rnSgn.uniform);
	_walk.rnSgn.setMinMax(-1,1);
	_walk.rnSgn.seedOffset(100+_IOcontrol.runOffset);
	_walk.Nparam=0;
	_cooling.convergenceTestMinimalLength=300;
	_cooling.maximalRejectionsCount=50;
	_cooling.maximalNumberOfForcedCoolings=10;
	_walk.maximalChainLength=200000;
	_walk.statesTot=0;
	_walk.initialPDFstepCRfraction=0.1;
	_walk.initialPDFstepCRfractionAfterBurnIn=0.05;
//	current()Direction=new double[Npar];
	_walk.currentDirection.makeLength(Npar);
	_walk.keepDirectionNow=false;
	_walk.keepDirectionStepGenerationFailuresNumber=50;
	_walk.rejectedLinkNumber=0;
	_walk.uphillClimbing=true;
	_walk.userOutputFrequency=1000;
	
	//
	// set the cooling PDF and corresponding RNG for "down-hill" walk
	//
	_cooling.initialTemperature=1000;
	_cooling.temperature=getInitialTemperature();
	_cooling.forcedCoolingTemperatureDecrement=-1;
	//	_temperatureHistory.append(_temperature);
	_cooling.finalTemperature=10;
	_cooling.initialMaximalEnergy=10; // i.e. 10 kT
	
	//
	// cooling variables
	//
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
	_cooling.coolingRNGseedOffset=101+_IOcontrol.runOffset;
	_cooling.coolingRNG.seedOffset(_cooling.coolingRNGseedOffset);
	_cooling.coolingRNG.seed(_walk.runRNGseed);
	_cooling.acceptWorseStateFromUniformDistr=true;
	setAcceptWorsePvalues(0.3,0.05);
	setAcceptWorseBoltzmannSchemeFactors(1,1);
	_cooling.coolingRate=0.5;
	
	_chisqData.bestFitNcovSimNum=5000;
	_chisqData.rnCov.setRNsType("gaussian_circle");
	_chisqData.rnCov.seedOffset(102+_IOcontrol.runOffset);
	_chisqData.rnCov.seed(_walk.runRNGseed);
	
	
	//
	// convergence variables
	//
	_convergence.relChisqChangeThres=0.001;
	_convergence.last_chisq=-1;
	_convergence.minimalLengthToTestConvergence=_cooling.convergenceTestMinimalLength;
	_convergence.nextConvergenceCheck=0;
	_convergence.statesToConvergenceCheck=1000;
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
void cpedsMCMC::addParameter(const mscsFunction& p) { 
	_parameterSpace.append(p); 
	_parameterSpace.last().checkRanges();

	_walk.Nparam++;
	// prepare RNGs for drawing starting point for the chain
	_walk.rns.append(cpedsRNG());
	_walk.rns.last().setRNsType(_walk.rns.last().from_array_invCDF);
	_walk.rns.last().setPDF(p.pointsCount(),p.extractArguments(),p.extractValues());
	_walk.rns.last().seedOffset(_walk.rns.size()+1+_IOcontrol.runOffset);
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
	_walk.stepGeneratingRNG.last().seedOffset(-_walk.stepGeneratingRNG.size()-1+_IOcontrol.runOffset);
	_walk.stepGeneratingRNG.last().seed(_walk.runRNGseed);
	
	
}
/***************************************************************************************/
void cpedsMCMC::addParameter(string param_name, double from, double to, long Npts) {
	mscsFunction param(param_name);
	double dx=(to-from)/Npts;
	param.mkConst(from,to,dx,1.0); 
	param.normalize();
	addParameter(param);
}

/***************************************************************************************/
void cpedsMCMC::startChain(bool followPriors) {
	
	MClink next(dims());
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
	msgs->say("Drawing random starting point in the prameter space",Medium);
	current().set(dims(),getStartingPoint());
	current().setAccepted(true);
	current().setIdx(0);
	_startingLink=current();
	// calculate the chisq there
	current().setChisq(chisq(current())); // get first chisq and assign chisq signature
	_convergence.last_chisq=current().chisq();
//	printInfo(); // print and store setup here to get the chisq assigned correctly before storing
	
	_MCaccepted.append(current());
	_bestFit=current();
	_cooling.temperatureHistory.append(getTemperature());
	
	bool walk=true;
	
	msgs->say("STARTING THE MC CHAIN",High);
	while (walk) {
		//
		// generate next step
		//
		_walk.statesTot++;
		if (length()==_cooling.convergenceTestMinimalLength) { 
			deltaX=current().chisq()-_bestFit.chisq();
			next=_bestFit;  // go back to the best fit link in the chain after the burn-in period and continue from there
			next.setIdx(_walk.statesTot);
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
			next.set(dims(),getNextPoint(current()));
			next.setIdx(_walk.statesTot);
			msgs->say("Current temperature is: "+msgs->toStr(getTemperature())+", final temperature is:"+msgs->toStr(getFinalTemperature()),Zero);
			next.setChisq(chisq(next));
		}
		
		
		
		//
		// decide what to do with the new state
		//
		
		// define the current accept worse threhold
		if (length()<=_cooling.convergenceTestMinimalLength) { acceptWorseThreshold=_cooling.acceptWorseThreshold; acceptWorseBoltzmannFactor=_cooling.acceptWorseBoltzmannFact; }
		else  { acceptWorseThreshold=_cooling.acceptWorseThresholdAfterBurnIn; acceptWorseBoltzmannFactor=_cooling.acceptWorseAfterBurnInBoltzmannFact; }
		
		
		if (next.chisq()<current().chisq()) {
			msgs->say("the new chisq is < than the current one: ACCEPTING the step",Zero);
			// take the step
			next.setAccepted(true);
			_MCaccepted.append(next);
			if (next.chisq()<_bestFit.chisq()) {
				_bestFit=next; 
				_MCbestFitChain.append(_bestFit); 
				_IOcontrol.storeBestFitNow=true; 
				storeBestFit();
				if (_walk.uphillClimbing) _walk.keepDirectionNow=true; 
			}
			_walk.rejectionsCount=0;
			// cool down
			msgs->say("Cooling down",Low);
			if (length() >= _cooling.convergenceTestMinimalLength) { // do not cool at all until the chain reaches minimal length suitable for testing the convergence
				/*
				if (length()==_convergenceTestMinimalLength) { 
					deltaX=_current.chisq()-_bestFit.chisq();
					_current=_bestFit;  // go back to the best fit link in the chain after the burn-in period and continue from there
					if (deltaX > 0) { // cool down at the same time to reduce chance of going away from the best fit solution
						coolDownLogLin(deltaX);
					}
				}
				else {
				 */
				if (getTemperature() > getFinalTemperature()) {
					//						deltaX=deltaChisq(length()-1); // commented out on Mar 17, 2011, 3:50:41 PM since we have accepted and rejected on separate lists.
					deltaX=current().chisq()-next.chisq();
					//						msgs->say("Delta chisq is: "+msgs->toStr(deltaX),Medium);
					//
					// here is the cooling scheme
					//
					if (deltaX > 0) {
						coolDownLogLin(deltaX);
						_acceptWorseRNs.clear();
					}
				}
				else { msgs->say("Cannot cool more. Reached the minimal temperature already Tmin= "+msgs->toStr(_cooling.finalTemperature),Low); }
				//				}
			}
			else { msgs->say("MC chain too short to cool down (accepted links: "+msgs->toStr(_MCaccepted.length())+", rejected links: "+msgs->toStr(_MCrej.length())+", together: "+msgs->toStr(length())+"). The minimal length to start cooling is: "+msgs->toStr(_cooling.convergenceTestMinimalLength),Low); }
			current()=next;
		}
		else {
			msgs->say("the new chisq is > than the current one",Zero);
			// take the step with probability according to the cooling PDF encoded in coolingRNG
#ifdef DEBUGMCMC
			// DEBUG STUFF - BEGIN
			if (length()>_convergenceTestMinimalLength) {
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
					acceptWorseState(true,next,_walk.rejectionsCount);
				}
				else acceptWorseState(false,next,_walk.rejectionsCount);
			}
			else {
				deltaX2perDOF=(next.chisq()-current().chisq())/dataSize();
				if (acceptWorseBoltzmannFactor*acceptWorseRN > deltaX2perDOF) { // take the step
					acceptWorseState(true,next,_walk.rejectionsCount);
				}
				else acceptWorseState(false,next,_walk.rejectionsCount);
			}
		}
		
		_cooling.temperatureHistory.append(getTemperature());
		//		printInfo();
		
		//
		// saving section
		//
		dumpAll();
		
		//
		// test convergence section NOT IMPLEMENTED YET
		//
		if (length()==_cooling.convergenceTestMinimalLength) { 
			msgs->say("Starting testing convergence at level: "+msgs->toStr(_convergence.relChisqChangeThres),Medium);
			_convergence.last_chisq=_bestFit.chisq(); 
			_convergence.nextConvergenceCheck=length()+_convergence.statesToConvergenceCheck;
		}
		if (length()==_convergence.nextConvergenceCheck) {
//			printf("bfX2: %lE bfX2last: %lE\n",_bestFit.chisq(),_convergence.last_chisq);
			if ((fabs(_convergence.last_chisq-_bestFit.chisq()))/_convergence.last_chisq < _convergence.relChisqChangeThres) {
				if (coolingPossible()) { 
					/* if the system is too hot then we can be here because the generated 
						steps are too big and MCMC departs too far from the bestFit solution and 
						the system won't have a chance to cool down. Cooling down here will force smaller steps 
						and will possibly improve the fit */
					msgs->say("Convergence (deltaX2<%lf) reached at T>Tfinal",_convergence.relChisqChangeThres,Medium);
					coolForcibly();
				}
				else {
					walk=false;
					msgs->say("Breaking the chain. Convergence (deltaX2<%lf) reached at T=Tfinal",_convergence.relChisqChangeThres,Medium);
				}
			}
			_convergence.last_chisq=_bestFit.chisq();
			_convergence.nextConvergenceCheck+=_convergence.statesToConvergenceCheck;
		}

		// test for maximal rejections count
		if (_walk.rejectionsCount>=_cooling.maximalRejectionsCount and length()>_cooling.convergenceTestMinimalLength) {
			if (forcedCoolingPossible()) { 
				msgs->say("Maximal rejections count (%li) reached reached at T>Tfinal",_walk.rejectionsCount,Medium);
				/* if the system is too hot then we can be here because the generated 
					steps are too big and MCMC departs too far from the bestFit solution and 
					the system won't have a chance to cool down. Cooling down here will force smaller steps 
					and will possibly improve the fit */
				coolForcibly();
			}
			else {
				walk=false;
			}
		}
		
		if (length()==_walk.maximalChainLength) walk=false;
		
		//
		// output
		//
		if (length() % getWalkInfoOutputFrequency() == 0 or walk==false or msgs->getVerbosity()>=High)
			msgs->say("MC length: %.0f, X2: %f, bf. X2: %f, cur. temp. %f, min. temp. %f ",double(length()),current().chisq(),_bestFit.chisq(),getTemperature(), getFinalTemperature(),Medium);
	}
	if (_walk.rejectionsCount>=_cooling.maximalRejectionsCount) { msgs->say("Breaking the chain. _maximalRejectionsCount reached: "+msgs->toStr(_walk.rejectionsCount),Medium); }
	if (length()==_walk.maximalChainLength) { msgs->say("Breaking the chain. The chain length reached: "+msgs->toStr(_walk.maximalChainLength),Medium); }
	
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
	if (length()>=_cooling.convergenceTestMinimalLength) {
		if (_walk.initialPDFstepCRfractionAfterBurnIn>0) frac=_walk.initialPDFstepCRfractionAfterBurnIn;
	}
	
	for (long i = 0; i < dims(); i++) {
		f.setName(_parameterSpace[i].getName()+"PDF");
		// define step size generating distribution
		//		f.mkBoltzmann(0,_initialMaximalEnergy*CPEDS_kB*_temperature,_initialMaximalEnergy*CPEDS_kB*_temperature/_stepGeneratingPDFPoints,_temperature);
		f.mkBoltzmann(0,_cooling.initialMaximalEnergy*CPEDS_kB*_cooling.initialTemperature,_cooling.initialMaximalEnergy*CPEDS_kB*getInitialTemperature()/_cooling.stepGeneratingPDFPoints,getTemperature());
		//		f.scaleX(0.5*(_parameterSpace[i].getMaxArg()-_parameterSpace[i].getMinArg())/f.getMaxArg());
		f.scaleX(frac*(_parameterSpace[i].getMaxArg()-_parameterSpace[i].getMinArg())/f.getMaxArg()); // make the maximal step size equal to half of parameter space prior - this is because the step can be positive or negative
		f.normalize(); // normalize
		_walk.stepGeneratingDistr[i]=f;
		
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
	for (long i = 0; i < dims(); i++) {		p[i]=_walk.rns[i].getRN();	}
	return p;
}
/***************************************************************************************/
double* cpedsMCMC::getNextPoint(const MClink& current) {
	msgs->say("generating the next point", Zero);
#ifdef DEBUGMCMC
	msgs->say("current state is:", Low);
	//	if (msgs->getVerbosity()>
	printLink(current);
	msgs->say("Best fit link is:", Low);
	printLink(_bestFit);	
#endif
	double* p=new double[dims()];
	//	double f; // step shrink factor
	double delta;
	double sign;
	double tmp;
	long failures;
	for (long i = 0; i < dims(); i++) {
		//		f=_temperature/_initialTemperature;
		failures=-1;
		do {
			failures++;
			delta=_walk.stepGeneratingRNG[i].getRN();
			
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
			//			msgs->say("new candidate: "+msgs->toStr(tmp)+", boundL: "+msgs->toStr(_parameterSpace[i].getX(0))+" boundU: "+msgs->toStr(_parameterSpace[i].getX(_parameterSpace[i].pointsCount()-1)), Low);
		} while (tmp < _parameterSpace[i].getX(0) || tmp > _parameterSpace[i].getX(_parameterSpace[i].pointsCount()-1));
		p[i]=tmp;
		//		msgs->say("accepted", Low);
	}
	_walk.keepDirectionNow=false;
	
	return p;	
}



/***************************************************************************************/
double cpedsMCMC::getSign() {
	double s=_walk.rnSgn.getRN();
	if (s<0) return -1;
	return 1;	
}
/***************************************************************************************/
void cpedsMCMC::setData(const mscsFunction& data) {	
	msgs->say("Setting data vector of length: "+msgs->toStr(data.pointsCount()),Medium);
	_chisqData.data=data; 	
	_chisqData.model=data; // this is to set the size and arguments
}
/***************************************************************************************/
const mscsFunction& cpedsMCMC::getData() const {
	return _chisqData.data;
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
	msgs->say("Burn-in period length (number of steps to start cooling down): "+msgs->toStr(_cooling.convergenceTestMinimalLength),Medium);
	msgs->say("maximal number of times to force cooling: "+msgs->toStr(_cooling.maximalNumberOfForcedCoolings),Medium);
	msgs->say("maximal number of rejections to force cooling or break chain: "+msgs->toStr(_cooling.maximalRejectionsCount),Medium);
	msgs->say("maximal chain length: "+msgs->toStr(_walk.maximalChainLength),Medium);
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
	msgs->say("partial files directory: "+_IOcontrol.partialFilesDir,Medium);
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
		_parameterSpace[i].save(getOutputDir()+"/parameterPrior.param."+msgs->toStr(i));
	}
}

/***************************************************************************************/
void cpedsMCMC::saveTemperature(string fname) {
	_cooling.temperatureHistory.save(fname);
}

/***************************************************************************************/
void cpedsMCMC::saveStepPDFs(string fname) {
	for (long j = 0; j < dims(); j++) {
		_walk.stepGeneratingDistr[j].save(fname+"-param_"+msgs->toStr(j)+".mc");
	}
}
/***************************************************************************************/
void cpedsMCMC::saveBestFitStepPDFs(string fname) {
	for (long j = 0; j < dims(); j++) {
		_chisqData.bestFitData.stepGeneratingDistr[j].save(fname+"-param_"+msgs->toStr(j)+".mc");
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

/***************************************************************************************/
/* NOTE: The indexing of partial files will be smaller by 1 wrt states number (statesTotal variable)
 * (for the case of dumping at every state), because the 
 * zero'th state is not saved, and the indexing variables used in this method work 
 * independently */

void cpedsMCMC::dumpAll() {
	char tmpch[10];
	string tmps;
	if (_IOcontrol.saveCoolingPDFEveryNthState!=0)
		if (_MCaccepted.size() % _IOcontrol.saveCoolingPDFEveryNthState == 0) {
			cpedsList<double> tmpl;
			sprintf(tmpch,"%05li",_IOcontrol.saveCoolingPDF_lastSavedNumber); tmps=tmpch;
			_cooling.coolingDistr.save(_IOcontrol.partialFilesDir+"/"+tmps+"-coolingPDF.step_"+msgs->toStr(_MCaccepted.size())); 
			_cooling.coolingDistrInvCDF.save(_IOcontrol.partialFilesDir+"/"+tmps+"-coolingInvCDF.step_"+msgs->toStr(_MCaccepted.size())); 
			_acceptWorseRNs.save(_IOcontrol.partialFilesDir+"/"+tmps+"-acceptWorseRNs.step_"+msgs->toStr(_MCaccepted.size())); 
			tmpl.append(_cooling.acceptWorseThreshold);
			tmpl.append(_cooling.acceptWorseThresholdAfterBurnIn);
			tmpl.save(_IOcontrol.partialFilesDir+"/"+tmps+"-acceptWorseThresholds.step_"+msgs->toStr(_MCaccepted.size())); 
			
			_IOcontrol.saveCoolingPDF_lastSavedNumber++;
		}
	if (_IOcontrol.saveStepPDFsEveryNthState!=0)
		if (_MCaccepted.size() % _IOcontrol.saveStepPDFsEveryNthState== 0) { 
			sprintf(tmpch,"%05li",_IOcontrol.saveStepPDFs_lastSavedNumber); tmps=tmpch;
			saveStepPDFs(_IOcontrol.partialFilesDir+"/"+tmps+"-stepPDFs.step_"+msgs->toStr(_MCaccepted.size())); 
			_IOcontrol.saveStepPDFs_lastSavedNumber++;
		}
	if (_IOcontrol.saveModelEveryNthState!=0)
		if (_MCaccepted.size() % _IOcontrol.saveModelEveryNthState== 0) { 
			sprintf(tmpch,"%05li",_IOcontrol.saveModel_lastSavedNumber); tmps=tmpch;
			_chisqData.model.save(_IOcontrol.partialFilesDir+"/"+tmps+"-model.step_"+msgs->toStr(_MCaccepted.size())); 
			_chisqData.chisqPerDOF.save(_IOcontrol.partialFilesDir+"/"+tmps+"-chisqPerDOF.step_"+msgs->toStr(_MCaccepted.size()));
			_IOcontrol.saveModel_lastSavedNumber++;				
		}
	if (_IOcontrol.dumpEveryMCMCstate) {
		saveAcceptedChisq(getOutputDir()+"/chisq.mc"); 
		saveParams(getOutputDir()+"/params.mc");
		saveChain(_MCrej,getOutputDir()+"/rejected.mc"); 
		saveChain(_MCaccepted,getOutputDir()+"/accepted.mc"); 
		QList<MClink> tmp=_MCrej; tmp.append(_MCaccepted);
		saveChain(tmp,getOutputDir()+"/all.mc"); 
		saveTemperature(getOutputDir()+"/temperature.mc");
		saveStepPDFs(getOutputDir()+"/stepPDFend");				
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
void cpedsMCMC::coolDownLogLin(double deltaX) {
	msgs->say("COOLING: current temperature: "+msgs->toStr(getTemperature()),Low);
	//	double deltaX2perDOF=deltaX/_data.pointsCount();
	
	if (log(deltaX) > 0) {
		_cooling.temperature-= getCoolingRate() * log(deltaX); // the cooling is proportional to log delta chisq
	}
	else {
		_cooling.temperature-= getCoolingRate() * deltaX; // the cooling is proportional to delta chisq							
	}
	
	
	if (getTemperature()<getFinalTemperature()) _cooling.temperature=_cooling.finalTemperature;
	
	msgs->say("COOLING: new temperature: "+msgs->toStr(getTemperature()),Low);
	
	updateStepGeneratingPDFs();		
	updateCoolingPDF();
	
}
/***************************************************************************************/
void cpedsMCMC::coolDownProp(double deltaX) {
	msgs->say("COOLING: current temperature: "+msgs->toStr(getTemperature()),High);
	double deltaX2perDOF=deltaX/dataSize();
	double chisqIni=_MCaccepted[_walk.rejectedLinkNumber].chisq()/dataSize();
	//	double dT=_finalTemperature/_initialTemperature*deltaX2perDOF*chisqIni;
	double dT=(_cooling.finalTemperature-_cooling.initialTemperature)*deltaX2perDOF/chisqIni; // this is under construction - what goes here ?
	
	_cooling.temperature-=_cooling.coolingRate*dT;
	
	
	if (getTemperature()<getFinalTemperature()) _cooling.temperature=getFinalTemperature();
	
	msgs->say("COOLING: new temperature: "+msgs->toStr(getTemperature()),Medium);
	
	updateStepGeneratingPDFs();		
	updateCoolingPDF();
	
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
	
	_startingLink.save(getOutputDir()+"/startingLink");
	saveData(getOutputDir()+"/inputData.txt");
	savePriors();
	_cooling.coolingDistr.save(getOutputDir()+"/coolingPDFini.mc");

}
/***************************************************************************************/
void cpedsMCMC::saveChain(QList<MClink>& chain, string fname) {
	if (chain.size()>0) {
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
		
		cpeds_matrix_save(m,fname);
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
	string cmd="if [ ! -d "+_IOcontrol.partialFilesDir+" ]; then mkdir -p "+_IOcontrol.partialFilesDir+"; fi";
	system(cmd.c_str());
}
/***************************************************************************************/
void cpedsMCMC::setOutputDir(string dirName) { 
	if (_IOcontrol.mcmcRunIdx==-1) {
		_IOcontrol.runFilesDir=dirName;
	}
	else _IOcontrol.runFilesDir=dirName+"."+msgs->toStr(_IOcontrol.mcmcRunIdx+_IOcontrol.dirStartIndex);
	_IOcontrol.partialFilesDir=_IOcontrol.runFilesDir+"/partial";
}
/***************************************************************************************/
void cpedsMCMC::saveFirstLinkLocation() {
	cpedsList<long> l;
	l.append(_walk.rejectedLinkNumber);
	l.save(_IOcontrol.partialFilesDir+"/first-link-index",false,"long");
}

/***************************************************************************************/
void cpedsMCMC::saveCurrentLength() {
	cpedsList<long> l;
	l.append(_MCaccepted.length());
	l.save(_IOcontrol.partialFilesDir+"/current-length",false,"long");
}
/***************************************************************************************/
void cpedsMCMC::saveCurrentTemperature() {
	cpedsList<double> ct;
	ct.append(_cooling.temperatureHistory.last());
	ct.save(_IOcontrol.partialFilesDir+"/current-temperature"); 			
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
	_cooling.convergenceTestMinimalLength=l; _convergence.minimalLengthToTestConvergence=l; 
	if (l>_walk.maximalChainLength) setMaximalChainLength(l);
}
/***************************************************************************************/
void cpedsMCMC::saveBestFitCovSimulations() {
	msgs->warning("you are calling an abstract saveBestFitCovSimulations() method. I hope this is what you want.", High);
}
/***************************************************************************************/
void cpedsMCMC::saveBestFit() {
	_chisqData.bestFitData.coolingDistr.save(_IOcontrol.partialFilesDir+"/bestFit-coolingPDF"); 
	saveBestFitStepPDFs(_IOcontrol.partialFilesDir+"/bestFit-stepPDFs"); 
	_chisqData.bestFitData.model.save(_IOcontrol.partialFilesDir+"/bestFit-model"); 			
	//		_temperatureHistory.save(_partialFilesDir+"/temperature-history"); 			
	_bestFit.save(_IOcontrol.partialFilesDir+"/bestFit-link");
	saveChain(_MCbestFitChain,getOutputDir()+"/bestFitChain");
	if (_chisqData.covarianceMatrixDiag.size()>0 and _chisqData.covDiagonal) {
		_chisqData.bestFitData.covDiag.save(_IOcontrol.partialFilesDir+"/bestFit-covarianceMatrixDiagonal",false,"double");
	}
	
	// save length
	cpedsList<long> l;
	l.append(_chisqData.bestFitData.lengthAtStoring);
	l.save(_IOcontrol.partialFilesDir+"/current-length",false,"long");
	
	//save temperature
	cpedsList<double> ct;
	ct.append(_chisqData.bestFitData.temperatureAtStoring);
	ct.save(_IOcontrol.partialFilesDir+"/current-temperature"); 			
	
	
}
/***************************************************************************************/
void cpedsMCMC::saveResults() {
	saveSetup();
	setDumpEveryMCMCstate(true);
	dumpAll();
	saveBestFit();
	save1Dposteriors();
	if (_IOcontrol.calculate2Dpdf) {
		save2Dposteriors();
		save2DCR(0.68);
		save2DCR(0.95);
		save2DCR(0.99);
	}

	saveBestFitCovSimulations();
}
/***************************************************************************************/
void cpedsMCMC::calculate1Dposteriors() {
	_chisqData.bestFitData.posteriors1D.clear();
	for (long i = 0; i < getNparam(); i++) {
		msgs->say("Calculating 1-D PDF for parameter: %li",i,Medium);
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
	if (_chisqData.bestFitData.posteriors1D.size()!=dims()) calculate1Dposteriors();
	for (long i = 0; i < getNparam(); i++) {
		_chisqData.bestFitData.posteriors1D[i].save(getOutputDir()+"/parameterPosterior.param"+msgs->toStr(i));
	}
	
}
/***************************************************************************************/
void cpedsMCMC::save2Dposteriors() {
	if (_chisqData.bestFitData.posteriors2D.size()!=double(dims()-1)/2*dims()) calculate2Dposteriors();
//	long k=0;
	for (long i = 0; i < getNparam(); i++) {
		for (long j = i+1; j < getNparam(); j++) {
			_chisqData.bestFitData.posteriors2D[paramij2idx(i,j)].saveHDF5(getOutputDir()+"/parameterPosterior.params"+msgs->toStr(i)+"_"+msgs->toStr(j)+".hdf5","L");
		}
	}
	
	//
	// 
	//
	
}
/***************************************************************************************/
void cpedsMCMC::save2DCR(double CL) {
	for (long i = 0; i < getNparam(); i++) {
		for (long j = i+1; j < getNparam(); j++) {
			double lvl;
			cpedsList<double> lvls;
			get2DCR(i,j,CL,&lvl).save(getOutputDir()+"/CR"+msgs->toStrf(CL,2)+".params"+msgs->toStr(i)+"_"+msgs->toStr(j));
			lvls.append(lvl);
			lvls.save(getOutputDir()+"/PDF-contour"+msgs->toStrf(CL,2)+".params"+msgs->toStr(i)+"_"+msgs->toStr(j));
		}
	}

}
/***************************************************************************************/
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
MscsPDF1D cpedsMCMC::get1Dposterior(int paramID, long pdfPoints) {
	mscsFunction states=getParamValues(paramID);

	states.checkRanges();
	states/=states.getMaxValue();
	states.removeSmaller(1e-10);
	states.checkRanges();
#ifdef DEBUG_MCMC_PDF
	cout << "Number of states after cutting off at L<1e-10: "<< states.pointsCount() <<"\n";
#endif
	states.sortFunctionArgAscending();
	
	MscsPDF1D pdf;
	
	
	cpedsList<long> bins;
	double pdfRange;//=sqrt(params.variance());
	pdfRange=states.getMaxArg()-states.getMinArg();
	long pdfPoints_internal=long(sqrt(states.pointsCount()));
	double dx=pdfRange/pdfPoints_internal;
	
	
	bins.clear();
	pdf=states.binFunction(dx,bins,"bin_center","max");
	pdf.checkRanges();
#ifdef DEBUG_MCMC_PDF
	pdf.save("rusty.pdf");
#endif

	pdf.interpolate(pdf.getMinArg(),pdf.getMaxArg(),pdfPoints,true,"cspline");
	pdf/=pdf.getMaxValue();

#ifdef DEBUG_MCMC_PDF
	pdf.save("rusty-interpolated.pdf");
#endif

	cpedsList<double> CR=pdf.getCR(0.9973).sort(12);
#ifdef DEBUG_MCMC_PDF
	CR.save("CR.tmp");
#endif
	if (CR.size()!=2) { cout << "WARNING 0.9973 confidence interval suspecious. Will return the rusty pdf.\n"; }
	else {
		double margin=CR[CR.size()-1]-CR[0];
		pdf=pdf.cut(CR[0]-margin/2,CR[CR.size()-1]+margin/2);
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
MscsPDF1D cpedsMCMC::get1Dposterior(string paramName, long pdfPoints) {
	return get1Dposterior( getParameterByName(paramName),pdfPoints );	
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
//	printf("p1min: %lf p1max: %lf\n",p1Min,p1Max);
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
//	printf("p2min: %lf p2max: %lf\n",p2Min,p2Max);
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
			if (pdf.fRe(i,j,0)>1e-5) {
				cpedsPoint3D p(pdf.getX(i),pdf.getY(j),pdf.fRe(i,j,0));
				ps.append(p);
			}
		}
	}
	assert(ps.size()>0);
//	ps.save("ps");
//	pdf.saveSlice(2,0,"pdf2d");
//	printf("%li\n",ps.size());
	pdf.mkInterpolatedFieldTriangLinear2D(ps);
	
	pdf/=pdf.getMaxValue();
//	pdf.saveSlice(2,0,"pdfint2d");
//	exit(0);


	
	
	
	
	
	
	
	
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
	_MCaccepted.append(chain.chain());
	_MCrej.append(chain.rejectedStates());
	_MCbestFitChain.append(chain._MCbestFitChain);
	if (chain.bestFitLink().chisq() < bestFitLink().chisq() or bestFitLink().chisq()==-1) {
		bestFitLink()=chain.bestFitLink();
		_chisqData.bestFitData=chain.getBestFitData();
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
	
	for (long i = 0; i < chain().size(); i++) {
		p.newPoint(chain(i).getParam(j),chain(i).L());
	}
	
	for (long i = 0; i < rejectedStates().size(); i++) {
		p.newPoint(rejectedStates()[i].getParam(j),rejectedStates()[i].L());
	}
	
	return p;
}
/***************************************************************************************/
int cpedsMCMC::paramij2idx(int i,int j) {
	int idx=0;
	for (long i = 0; i < dims(); i++) {
		for (long j = i+1; j < dims(); j++) {
			idx++;
		}
	}
	return idx-1;
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
	id+=_cooling.coolingRNG.seed();
	_cooling.coolingRNG.seed(id);
	for (long i = 0; i < _walk.stepGeneratingRNG.size(); i++) { _walk.stepGeneratingRNG[i].seed(id); }
	for (long i = 0; i < _walk.rns.size(); i++) { _walk.rns[i].seed(id); }
	_walk.rnSgn.seed(id); 
	_chisqData.rnCov.seed(id);
}
/***************************************************************************************/
void cpedsMCMC::coolForcibly() {
	if (_cooling.forcedCoolingTemperatureDecrement==-1) _cooling.forcedCoolingTemperatureDecrement=(getTemperature()-getFinalTemperature())/_cooling.maximalNumberOfForcedCoolings; 
	if (_cooling.forcedCoolingTemperatureDecrement< 0.01*(_cooling.initialTemperature-_cooling.finalTemperature)) _cooling.forcedCoolingTemperatureDecrement=0.01*(_cooling.initialTemperature-_cooling.finalTemperature);
	_cooling.temperature-=_cooling.forcedCoolingTemperatureDecrement;
	msgs->say("Forcibly cooling by %f K",_cooling.forcedCoolingTemperatureDecrement,Medium);
	if (getTemperature()<getFinalTemperature()) _cooling.temperature=getFinalTemperature();
	_cooling.forcedCoolingTemperatureDecrement=-1;
	
	updateCoolingPDF();
	updateStepGeneratingPDFs();
	_walk.rejectionsCount=0;
}
