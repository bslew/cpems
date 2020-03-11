/*!
 * \file cpedsMCMC.h
 *
 *  Created on: Nov 19, 2010
 *      Author: blew
 *      
 *      TODO: add saving the state of RNGs (low priority)
 *      TODO: add saving the whole MCMC (params, chisq, L) in one file
 *      TODO: store rejected links in a separate list (and store in the end) done
 *      
 *      TODO: what I still don't like is that the cooling is proportional to improvements in chisq.
 *      This is nice generally, but if the worse fit is accepted and then a better one is found again this
 *      also leads to cooling down even though the overall fit might have not improved.
 *      A workaround is to practically implement the option setAcceptWorseUniformScheme(false) which
 *      enables using boltzmann distribution to decide how bad the new state is and whether to take it or not.
 *      
 *      IDEA: Its possible to implement a communication between chains so that eg:
 *      - all chains after burn-in period switch to explore the region of the parameter space that best fits the data
 *      - the best fit chisqperDOF value was dynamically used between all chains to control the AcceptWorseBoltzmannSchemeFactors
 *      especially after the burn-in period. The best fit chain chisqperDOF value gives an estimate on how well the data can be fit
 *      with this model. Changing dynamically AcceptWorseBoltzmannSchemeFactors could shorten the total run time by stopping chain before
 *      the maximalChainLength condition is met simply by setting an appropriate setAcceptWorseBoltzmannSchemeFactors value.
 *      
 */

#ifndef CPEDSMCMC_H_
#define CPEDSMCMC_H_

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */

/* STANDALONE HEADERS */
#include "Mscs-function.h"
#include "Mscs-function3dregc.h"
#include "MClink.h"
#include "cpedsMC.h"
#include "cpeds-msgs.h"
#include "matrix.h"
#include "MscsPDF1D.h"
#include "MCparameterSpace.h"

/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */

/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
 \class 
 \brief Encapsulates the basic functionality of the MCMC for any number of dimensions.
 \details 
 
	 This is an abstract class and must be inherited by the class that will actually solve a problem.
	 In the derived class the method chisq must be implemented so that the chisq value is correctly 
	 calculated for a particular problem.

	The important part is how various MC RNGs are initialized with seeds. Basically they are all initialized with the seed based on time
	but they are differenciated based on different seed offsets.
	
	RNs for generating starting point for MCMC are started with seedOffsets proportional to the parameter 5*index+1
	RNs for generating MCMC steps are started with seedOffsets proportional to the parameter -5*index-1
	Cooling rng is started with seedOffset 100 - which guarantees fine initialization for 20 parameters
	The sign generagin rng is started with seed offset 0.
	The RNG for deriving covariance matrix from simualted data based on some model is started with seed offset 200
	
	 CHANGELOG:\n
	 * added support for uphill (directional) climbing depending on acceptance/rejection result from the previous state (Feb 3, 2011, 2:49:03 PM )
	 
	 * added support for cooling after the burn in period when returning to the best fit solution


 
 \date Nov 19, 2010, 2:35:34 PM 
 \author Bartosz Lew
 */

class cpedsMCMC {
	/***************************************************************************************/		
	/***************************************************************************************/		
	/***************************************************************************************/		
	/***************************************************************************************/		
	/***************************************************************************************/		
	/***************************************************************************************/		
	/***************************************************************************************/		
	/***************************************************************************************/		
	/***************************************************************************************/		
	public:
		typedef	enum { awSchemeUniform, awSchemeBoltzmann } acceptWorseScheme;
		typedef	enum { coolingScheme_none, coolingScheme_chisqPropLinLog, coolingScheme_chisqPerDoFprop, coolingScheme_geometric } coolingSchemeTraits;
		
		typedef struct {
				long stepGeneratingPDFPoints; //!< number of points that probe the cooling distributions
				long maximalRejectionsCount; //!< the maximal number of rejections of step proposals in the row after which the system will be forced to cool or the chain will be broken
				long maximalNumberOfForcedCoolings; //! the number of times when the system will be forced to cool down each time when _maximalRejectionsCount has been reached and _temperature > _finalTemperature
				
				double temperature; //!< defines the current temperature for the simulated annealing and its history
				mscsFunction temperatureHistory; //!< records the cooling history
				double initialTemperature; //! temperature to start the walk with
				double finalTemperature; //!< temperature at which we stop the parameter space walk
				double coolingRate; //!< temperature decrease multiplier for one Monte-Carlo step
				coolingSchemeTraits coolingScheme; //!< defines show the system cools upon chisq improvements
				 
				double forcedCoolingTemperatureDecrement;
				double initialMaximalEnergy; //!< energy in _initialTemperature * k_B units that define the initial range of PDF covered for step size calculations
				
				bool acceptWorseStateFromUniformDistr; //!< flag defining if the worse state should be accepted based on the fixed thresholds defined below and uniform distribution or based on Boltzmann distribution and delta chisq
				double acceptWorseThreshold; //!< threshold defining whether or not to accept the worse fit as a next state for a given system temperature during burn-in period.
				double acceptWorsePvalue; //!< corresponding P-value defining whether or not to accept the worse fit as a next state for a given system temperature during burn-in period.
				double acceptWorsePvalueAfterBurnIn; //!< P-value defining whether or not to accept the worse fit as a next state for a given system temperature after burn-in period.
				double acceptWorseThresholdAfterBurnIn; //!< threshold defining whether or not to accept the worse fit as a next state for a given system temperature after burn-in period.
				double acceptWorseAfterBurnInBoltzmannFact; //!< factor by which a RN is multiplied and compared against reduced chisq value to decide whether to accept the worse state than the current one after the burn-in period
				double acceptWorseBoltzmannFact; //!< factor by which a RN is multiplied and compared against reduced chisq value to decide whether to accept the worse state than the current one, during the burn-in period
				
				cpedsRNG coolingRNG;
				long coolingRNGseedOffset;
				mscsFunction coolingDistr; //!< 
				mscsFunction coolingDistrInvCDF;
		} cooling_t;
		
		typedef struct {
				int Nparam; //!< number of parameters
				QList<double> delta0; //!< defines the grid separations based on _parameter space give
						
				QList<mscsFunction> stepGeneratingDistr; //!< cooling distribution functions - one for each dimension
				QList<cpedsRNG> stepGeneratingRNG; //!< RNGs - one for each dimension

				QList<cpedsRNG> rns; //!< RNGs - one for each dimension: these are used to start the chain from random location

				cpedsRNG rnSgn; //!< This generator only generates the sign + or -
				MClink current; //!< The last accepted MClink
				MClink next; //!< currently not used
				MClink avgStep; //!< average step calculated using NavgSteps
				double likelihoodRejectionThreshold; //!< likelihood threshold below which links are rejected for posterior calculation
				long maximalChainLength; //!< the maximal length of the chain after which the chain will be broken
				long maximalAcceptedChainLength; //!< the maximal number of accepted links (default: 0 = not used)
				long statesTot; //!< total number of states (accepted and rejected) considered during the walk
				long statesTotBeforeBurnOut; //!< chain length after the walk is done and before burning-out
				double initialPDFstepCRfraction; //!< fraction of the parameter domain used for generation of initial step size PDFs
				double initialPDFstepCRfractionAfterBurnIn; //!< fraction of the parameter domain used for generation of initial step size PDFs after burn in period is finished (by default: 0 - ie. not used)
				long runRNGseed;
				bool poorX2improvement; //!< this flag becomes true whenever poor X2 improvement is detected; this may influence cooling
				bool uphillGradient; //!< this flag modifies the way of taking the next step; it calculates grad(chisq) and takes the steepest negative direction with step size proportional to the gradient and the initial step size (default: false)
				bool uphillClimbing; //!< this defines whether to keep the direction of walk if the previous state was accepted and deltaX2>0 (true) or otherwise whether the next step direction should be random (default: true)
				bool keepDirectionNow; //!< defines whether to keep the current direction when _uphillClimbing is on
				cpedsList<double> currentDirection; //!< holds the current climb direction

				cpedsMC steps; //!< holds history of MC steps in each direction (regardless of whether the step is accepted or not)
//				cpedsMC stepsBuff; //!< holds history of MC steps in each direction (regardless of whether the step is accepted or not)
//				vector< cpedsList<double> > steps; //!< holds history of MC steps in each direction (regardless of whether the step is accepted or not)
				vector< cpedsList<double> > stepsBuff; //!< holds history of NavgSteps accepted MC steps 
				long keepDirectionStepGenerationFailuresNumber; //!< maximal number of failures to generate the next step in given direction before dropping that direction (this prevents chain from trying to keep direction when staying on wall)
				long rejectedLinkNumber; //!< number of links that were prepended before the beginning of the chain to store the location of the true first link in the chain. -- not sure if this is accurate anymore...
				long rejectionsCount; //!< current number of rejected states
				long userOutputFrequency; //!< output info about MCMC walk every this number of MCMC states
				long convergenceTestMinimalLength; //!< the minimal length of the chain for testing the convergence; the cooling is not done for chains shorter than this value
				long burnOutLength; //!< number of steps to perform after convergence is reached.
				long NavgSteps; /*!< number of last steps used to calculate average step.
				 	 	 	 	 	 Used only for UphillGradient method to speed up
				 	 	 	 	 	 the convergence. */
				long momentum; /*!< number of accumulated average steps.
				 	 	 	 	 	 Used only for UphillGradient method to speed up
				 	 	 	 	 	 the convergence. */
				
				mscsFunction momentumFactorHistory;
				cpedsMC momentumHistory;
				cpedsMC etagradHistory;
				mscsFunction etaHistory;
				double d2X2rel;
				int d2X2rel_points; //!< number of last points used for d(deltaX2rel)/d step calculation
				mscsFunction dX2rel; //!< delta(chisq)/chisq used to control gradient walk parameters
				
		} walk_t;
		
		typedef struct {
				long minimalLengthToTestConvergence; //!< number of states to test the convergence; by default the length of the burn in period
				long nextConvergenceCheck; //!< counter
				long statesToConvergenceCheck; //!< number of states every which the convergence is checked
				double relChisqChangeThres; //!< fraction by which the variance must change at the least to continue the walk
				double convLastVariance; // not used
				double last_chisq; //!< chisq ast last convergence check
				cpedsList<double> chisq; //!< chisq values (both rejected and accepted) used to estimate the convergence
				 
		} convergence_t;
		
		typedef struct {
				mscsFunction coolingDistr; //!< 
				QList<mscsFunction> stepGeneratingDistr; //!< cooling distribution functions - one for each dimension
				mscsFunction model; //!< contains the best fit model 
				cpedsList<double> covDiag; //!< diagonal of the covariance matrix
				long lengthAtStoring;
				double temperatureAtStoring;
				QList<MscsPDF1D> posteriors1D; //!< one for each parameter
				QList<mscsFunction3dregc> posteriors2D; //!< one for each combination of parameters
//				map<int, int> params2idx;
				cpedsList<double> confidenceLevels; //!< CLs for which CR are calculated and stored. By default 68, 95 and 99 %
				QList<mscsFunction3dregc> confidenceRanges1D; /*!< 1D CRs, one per parameter. Each element is a 3xconfidenceLevels.size() matrix
					Each row of the matrix contains data for each CL. 
					1st column is the CL,\n
					2nd column is CRlow, \n
					3rd column is CRhigh.\n
					4th column is the PDF value that corresponds to the [CRlow,CRhigh] range
				*/
		} bestFitData_t;
		
		
		typedef struct {
				mscsFunction3dregc data2d; /*!<	contains the full input data structure; the column data2dcol is 
				 * used to evaluate chisq; the last 2 columns will contain best fit model and its residuals respectively
				 * The data structure below is initialized with the column data2dCol from data2d but it 
				 * may becore depreciated.
				 * If data2d is not initialized, then residuals will not be calculated.
				 * The residuals could be calculated using data below but often it is useful to combine the residuals
				 * with some of the other columns from the input data file for plotting against other columns which
				 * are missing in the data structure below. So it is convenient to provide the full input data and indicate
				 * which column should be used for fitting, then store the full input data along with two extra columns
				 * containing the best fit model and the resuduals.
				 */
				long data2dCol;
				mscsFunction data; //!< contains data used to evaluate chisq for various parameters of the model
				mscsFunction model; //!< contains the last derived model according to the _current set of parameters
				cpedsRNG rnCov; //!< This generator is used to generate the simulated data to derive the covariance matrix at every state
				cpedsList<double> chisqPerDOF; //!< holds the current chisq per DOF values
				cpedsList<double> dataMask; //!< the data mask list containing weights (0 or 1) for the data to be masked out. (Basically other weights are also possible)
				bool covDiagonal; //!< if true - _covDiag is assumed to hold the current diagonal covariance matrix. This might be is useful as a flag for saving and writting your reimplemented chisq() functions in derived class
		//		double* _covDiag; //!< diagonal of the covariance matrix -- THIS IS DEPRECIATED. Will use now _covarianceMatrixDiag
				cpedsList<double> covarianceMatrixDiag; //!< diagonal covariance matrix 
				double variance; //!< diagonal covariance matrix same for all data points
				mscsFunction3dregc covarianceMatrix; // full covariance matrix

				long bestFitNcovSimNum; //!< number of simulations to generate and save after the walk has finished for the best fit model. Used in derived class
				bestFitData_t bestFitData; //!< this is updated every time a better fit is found, and saved at the end of the walk
				
		} chisqData_t;
		
		
		typedef struct {
				int mcmcRunIdx; //!< identifier of the run used for uniquely naming the run directory
				long runOffset; //!< the global MCMC run offset used to initiate the RNGs with different seeds across different MCMC runs
				string partialFilesDir;
				string runFilesDir;
				string chisqSignature;
				long dirStartIndex; //!< with which to start indexing the output directories
				long saveCoolingPDFEveryNthState, saveCoolingPDF_lastSavedNumber;
				long saveStepPDFsEveryNthState, saveStepPDFs_lastSavedNumber;
				long saveModelEveryNthState, saveModel_lastSavedNumber;
				bool dumpEveryMCMCstate; //!< toggles if dump chisq, likelichood etc at model saving rate: useful for tracking the run on-fly
				bool storeBestFitNow; //!< controls storing best fit models
				QList<string> implementationNote; //!< this list will be stored in the setup file as an extra block. You can describe here the features of your particular implementation of the derived class
				long PDF1d_pointsCount; //!< number of points used to store estimated 1d PDFs
				long PDF2d_pointsCount; //!< number of points used to store estimated 2d PDFs
//				bool calculate1Dpdf;				
				bool calculate2Dpdf;	
				int saveFilePrecision; //!< parameter sent to setprecision in ostream when saving
		} IOcontrol_t;
		
//		cpedsMCMC(int runIdx=0);
		
		/*!
			\brief constructor of the MCMC class
			\details 
			@param Npar - number of parameters that will be used
			@param runIdx - unique MCMC run index. This number is used to uniquely name the directory storing all files generated by the run and also to
			provide appropriate run offset for RNGs when many parallel MCMC are started. If not given or -1 then it is not used for defining the output directory 
			@param runOffset - offset of the run used to uniquely initialize the RNGs with some seed. This number is multiplied by runIdx to obtain the final run offset that will be used.
			In practice this class will use a range of seed offsets randing from -100 to 100 regardless of the seed value. Using runOffset say 1000 will ensure that every MCMC run will be 
			started with possibly the same seed value (taked from the start time) will be offsetted by the unique number of the runIdx multiplied by the runOffset.
			@return
		
			\date Feb 3, 2011, 11:32:13 AM
			\author Bartosz Lew
		*/
		cpedsMCMC(long Npar=1,int runIdx=-1, long runOffset=1000, long runSeed=0, cpeds_VerbosityLevel verb=Medium);
		virtual ~cpedsMCMC();
		cpedsMCMC(cpedsMCMC &parent);
		
		//! initialize the basic parameters
		void initialize(int runIdx, long runOffset=0,long Npar=1, long runSeed=0);

		void setRunIdx(int idx) { _IOcontrol.mcmcRunIdx=idx; }
		int getRunIdx() const { return _IOcontrol.mcmcRunIdx; }
		string getRunDependentFileName(string fname);
		string getParameterDependentFileName(string fname, long paramID);
		string get2ParameterDependentFileName(string fname, long paramID1, long paramID2);
		void setVerbocity(cpeds_VerbosityLevel verb);
		
		/*!
			\brief chisq evaluation - abstract method to be implemented in the derived class fit for the particular problem.
			\details 
			@param data - data vector
			@param th - set of model parameters used to derive the expected values for every measured data point
			@param cov - covariance matrix for the input data; the reference is passed so that the object can be modified if
			the covariance matrix changes from one state to another which is generally true
			@return
		
			\date Nov 19, 2010, 10:04:48 PM
			\author Bartosz Lew
		*/
//		virtual double chisq(const matrix<double>& data, const MClink &th, const matrix<double>& cov);
		virtual double chisq(const mscsFunction& data, const MClink &th, matrix<double>& cov);
		virtual double chisq(const MClink &th);

		double getRelX2Improvement();
		double getX2Convergence();
		
		/*!
			\brief set the data for chisq calculation to be used throughout the chain walk
			\details 
			@param data - mscsFunction defining the data vector 

			NOTE: this method is depreciated in favour of setData(const mscsFunction3dregc& data);
		
			\date Dec 13, 2010, 7:23:59 PM
			\author Bartosz Lew
		*/
		void setData(const mscsFunction& data);
		/*!
			\brief provide 2D data array
			\details 
			@param data - 2d data array that will be used during chisq calculations
			@param colNo - column number in that array which holds the values that should be fitted

			The reason for storing the full 2D array rather than just a 1D vector of the values
			is that (i) sometimes the chisq depend on a number of other columns as well and (ii)
			it is convenient to store the residuals for the best fit model along with the 
			all other columns to ease plotting against other columns in the data. So this saves the need
			of reprocessing the outputs after the fitting. However, the 2d data need not be set in which
			case the residuals will not be calculated this way.
		
			\date Feb 14, 2018, 2:43:30 PM
		*/
		void setData2D(const mscsFunction3dregc& data, int colNo);
		const mscsFunction& getData() const;
		const mscsFunction& getModel() const;
		const mscsFunction3dregc& getData2D() const;
		double getData2D(long i, long j) const;

		void setDataMask(const cpedsList<double>& dataMask) { _chisqData.dataMask=dataMask; }
		long dataSize() const { return _chisqData.data.pointsCount(); }

		/*!
			\brief set the covariance matrix for all the states in the chain walk
			\details 
			@param cov - covaraince matrix that is independent from the values of the parameters 
		
			\date Dec 13, 2010, 7:22:21 PM
			\author Bartosz Lew
		*/
//		void setCovarianceMatrix(const matrix<double>& cov);
		void setCovarianceMatrix(const mscsFunction3dregc& cov);
		void setDiagonalCovarianceMatrix(const cpedsList<double>& cov);
		void setDiagonalCovarianceMatrix(double var);
		double getCov() { return _chisqData.variance; }
		double getCov(int i) { return _chisqData.covarianceMatrixDiag[i]; }
		double getCov(int i, int j) { return _chisqData.covarianceMatrix(i,j); }
		cpedsList<double>& covDiag() { return _chisqData.covarianceMatrixDiag; }
		mscsFunction3dregc& cov() { return _chisqData.covarianceMatrix; }
		double getConvergenceThres() const {	return _convergence.relChisqChangeThres;	}
		void setConvergenceThres(double relChisqChangeThres) { this->_convergence.relChisqChangeThres = relChisqChangeThres; }

		
		/*!
			\brief add a new parameter to the parameter space
			\details 
			@param p - function defining the parameter prior. The viable parameter ranges are defined by the arguments of the function.
			This function also defines the grid spacing in case of grid calculations
			@param parameter_full_name - a string defining a latex stype parameter name and units in which this parameter likelihood will be stored.
			It is useful to provide this string here because the string will be stored in the output likelihood hdf5 file
			and used latter for plotting with the provided plot_likelihood.py python script. (e.g. '$v_x$ [km/h]')
			
			The name of the parameter will be the same as the function name.
		
			\date Nov 19, 2010, 10:02:43 PM
			\author Bartosz Lew
		*/
		void addParameter(const mscsFunction& p, string parameter_full_name="");
		
		/*!
			\brief convenience function. Adds flat parameter prior defined on [from,to] and resolved with Npts
			\details 
			@param Npts - number of points that resolve flat prior function
		
			\date May 29, 2017, 10:02:09 PM
		*/
		void addParameter(string param_name, double from, double to, long Npts=1000, string param_full_name="");
		void addParameter(string param_name, double from, double to, double delta, string param_full_name);
		
		/*!
			\brief get parameter data by name
			\details 
			@param Pname - name of the parameter to get info about
			@return const reference to the parameter in formation
			
			If no parameter by the provided name was found then the function will return the first parameter information
			and the control flag ctrl will be set of -1 otherwise it will be set to 0;
			If ctrl is NULL then no control information will be returned.
		
			\date Jan 28, 2011, 2:08:43 PM
			\author Bartosz Lew
		*/
		const mscsFunction& getParameter(string Pname, int* ctrl=NULL) const;

		int getParameterByName(string Pname) const;

		int dims() const { return _walk.Nparam; }
		long sizeRejected() { return _MCrej.size(); }
		long sizeAccepted() { return _MCaccepted.size(); }
		long size() { return sizeAccepted()+sizeRejected(); }
		long length() { return _MCaccepted.size()+_MCrej.size(); }
		long lengthAccepted() { return _MCaccepted.size(); }
		cpedsMC& BFchain() { return _MCbestFitChain; }
		MClink& chain(long i) { return _MCaccepted[i]; }
		cpedsMC& rejectedStates() { return _MCrej; }
		cpedsMC& testStates() { return _MCtest; }
		cpedsMC& chain() { return _MCaccepted; }
		cpedsMC& acceptedStates() { return _MCaccepted; }
		MClink getBestFitLink() { return _bestFit; }
		MClink getStartingLink() { return _startingLink; }
		MClink generateStartingLink() { MClink l(getNparam()); l.setLink(getNparam(),getStartingPoint(),-1); return l; }
		mscsFunction getBestFitModel() { return _chisqData.bestFitData.model; }
		MClink& bestFitLink() { return _bestFit; }
		void addNewBestFitLink(MClink& l);

		/*!
			\brief starts a chain walk from a random place in the parameter space
			\details 
			@param followPriors - if true then the random point will be drawn according to the priors defined in the parameter space; 
			if false then flat priors will be used for all parameters (currently this is ON by default and the only one that works)
			@param startingLink - if given then it is used as starting link instead of the randomly generated one
				from the distribution of priors. Useful for debugging.
		
			\date Nov 19, 2010, 10:07:40 PM
			\author Bartosz Lew
		*/
		void startChain(MClink* startingLink=0, bool followPriors=true);
		/*!
			\brief run MCMC for Nstates another states around the starting link
			\details 
			@param how - 0 - walk freely; 1 - take always one step from the best fit link;
			2 -  take always one step from the best startingLink;
			@return
		
			\date Jan 25, 2018, 6:56:54 PM
		*/
		void walkNstates(long Nstates, MClink startingLink, int how=1);
		
		/*!
			\brief run MCMC for another N states in parallel assuming that steps are independent
			\details 
			@param N - number of steps that should be taken

			This method is useful for burn-in and burn-out phases where the next step
			does not depend on the current step.
			The next step is taken from the current PDFs for each dimension.
			
			\date Mar 6, 2020, 12:13:02 PM
		*/
		void walkNindependent_states(long N);
		void setMaximalChainLength(long l) { _walk.maximalChainLength=l; }
		void setMaximalAcceptedChainLength(long l) { _walk.maximalAcceptedChainLength=l; }
		/*!
			\brief set the maximal number of rejected states required for forced cooling
			\details 
			@param n - a reasonable practise is about 30*sqrt(number_parameters)
		
			\date Jan 9, 2018, 1:20:33 PM
		*/
		void setMaximalRejectionsCount(long n) { _cooling.maximalRejectionsCount=n; }
//		void setProbeRandomlyAtBurnIn(bool tf) { _walk; }
		void setBurnInLength(long l);
		long getBurnInLength() const { return _walk.convergenceTestMinimalLength; }
		void setBurnOutLength(long l);
		long getBurnOutLength() const { return _walk.burnOutLength; }
		void burnOut();
		/*!
			\brief deinfes the number that starts the indexing of the output directories for mcmc data
			\details 
			@param idx - starting index 
		
			\date Mar 4, 2011, 1:02:37 AM
			\author Bartosz Lew
		*/
		void setDirStIdx(long idx) { _IOcontrol.dirStartIndex=idx; }
		void updateStepGeneratingPDFs();
		void updateStepGeneratingPDF(long param, double expectedStepSize);
		void updateCoolingPDF();
		
		
		void printInfo() const;
		void printLastStep() const;

		QList<MscsPDF1D>& posteriors() { return _chisqData.bestFitData.posteriors1D; }
		MscsPDF1D get1Dposterior(string paramName, long pdfPoints=50, string interpolationType="steffen");
		MscsPDF1D get1Dposterior(int paramID, long pdfPoints=50, string interpolationType="steffen", bool recalculate=true);
		mscsFunction3dregc get2Dposterior(string paramName1, string paramName2, long pdfPoints=50);
		mscsFunction3dregc get2Dposterior(int paramID1, int paramID2, long pdfPoints=50);

		/*!
			\brief define the frequency at which the MCMC run states information is dumped defined in number of states between dumps
			\details 
			@param saveCoolingPDFEvery - how often to save the cooling pdf
			@param saveStepPDFsEvery - how often to save the step pdfs
			@param saveModelEvery - how often to save the fitted model
			@return
		
			\date Jan 28, 2011, 1:54:09 PM
			\author Bartosz Lew
		*/
		void setSaveParameters(long saveCoolingPDFEvery, long saveStepPDFsEvery, long saveModelEvery);
		
		void saveAcceptedChisq(string fname);
		cpedsMC loadMC(string fnameName);
		void loadMCrun(string dirName);
		void saveAcceptedParams(string fname);
		void saveMomentumHistory(string fname);
		void saveLearningRateHistory(string fname);
		void saveData(string fname);
		/*!
			\brief calculate and save residuals calculated wrt the provided input data column
			\details  
			@param inputDataColumn - column in the input data used to calculate residuals; 
			by default - the one indicated through setData2d()

			
		
			\date Feb 7, 2018, 3:31:44 PM
		*/
		void saveResiduals(int inputDataColumn=-1);

		void savePriors();
		void save1Dposteriors(string output_file_prefix="posterior1D");
		void save2Dposteriors();
//		void setCalculate1Dposteriors(bool tf) {_IOcontrol.calculate1Dpdf=tf; }
		void setCalculate2Dposteriors(bool tf) {_IOcontrol.calculate2Dpdf=tf; }
		/*!
			\brief re-calculate 1-D likelihoods based on all currently available links
			\details This method forces re-calculating the likelihoods.
		
			\date Mar 11, 2020, 3:27:03 PM
		*/
		void calculate1Dposteriors();
		/*!
			\brief calculate 1D CRs
			\details 
		
			\date Feb 24, 2018, 3:32:22 PM
		*/
		void calculate1DCRs();
		void calculate2Dposteriors();
		/*!
			\brief save vertices of a patch enclosing requested confidence level (CL)
			\details 
			@param CL - confidence level
		
			\date May 31, 2017, 10:28:34 AM
		*/
		
		void recalculateLikelihoods();
		long statesAboveRejectionThreshold();
		void save2DCR(double CL);
//		void save1DCR(double CL, string fname, string dset);
		cpedsList<double> get1DCR(int paramID1, double CL, double* LVL=NULL);
		mscsFunction get2DCR(int paramID1, int paramID2, double CL, double* LVL=NULL);
	
		void saveTemperature(string fname);
		void saveBestFitStepPDFs(string fname);
		void loadBestFitStepPDFs(string dirName);
		void saveStepPDFs(string fname);
//		mscsFunction getMCstepsFn(long param) const;
		cpedsMC getMCsteps() const { return _walk.steps; }
		cpedsMC& MCsteps();
		void saveMCsteps(string fname);
		void saveParameterNames(string fname="");
		/*!
			\brief sets a flag indicating whether or not the covariance matrix of the model in the problem is to be considered as diagonal
			\details 
			@param tf - true or false
			If true, the best fit diagonal matrix will be read from a linear structure and stored as a vector.
			By default the flag is set to false.
		
			\date Feb 2, 2011, 11:51:07 AM
			\author Bartosz Lew
		*/
		void setCovarianceMatrixDiagonal(bool tf) { _chisqData.covDiagonal=tf; }
		
		string getPartialFilesDirFull() const { return _IOcontrol.partialFilesDir; }
		string getPartialFilesDir() const { return "partial"; }
		
		void setCoolingRNGseedOffset(long seed) { _cooling.coolingRNGseedOffset=seed; }
		void setCoolingRate(double coolingRate) { _cooling.coolingRate=coolingRate; }
		void setCoolingScheme(coolingSchemeTraits coolingScheme) { _cooling.coolingScheme=coolingScheme; }
		coolingSchemeTraits getCoolingScheme() { return _cooling.coolingScheme; }
		void setDumpEveryMCMCstate(bool tf) { _IOcontrol.dumpEveryMCMCstate=tf; }
		void setOutputDir(string dirName);		
		void setUpHillClimbing(bool tf) { _walk.uphillClimbing=tf; }
		/*!
			\brief specifies how the next step is chosen
			\details 
			@param tf - if true then the next step is chosen by probing the d(chisq)/dth_i in each dimension and setting the
			next step in each direction that is proportional to the gradient in that direction. 

			The step size is controlled by the value set by setInitialWalkStepSize() method and the size of the parameter space
			along that direction.
			
			Setting this option to true, makes setUpHillClimbing() settings irrelevant,
			and also resets the cooling scheme to cpedsMCMC::coolingScheme_none (the cooling will be done abruptly only when convergence is 
			reached e.g. in cases when the steps are too big). In this mode the temperatures moderate the step sizes taken along the gradient.
			The lower temperatures the smaller steps.


			In this mode the rejections count is not incremented when poor X2 improvement is detected. ? is this useful
			and also resets the current, initial and final temperatures to 1 since the temperatures and cooling are irrelevant in this mode.? is this useful
		
			\date May 5, 2018, 9:04:25 PM
		*/
		void setUpHillGradient(bool tf) { 
			_walk.uphillGradient=tf; 
			if (_walk.uphillGradient) {
//				setTemperatures(1,1);
//				setTemperature(getInitialTemperature());
				setCoolingScheme(cpedsMCMC::coolingScheme_none);
			}
		}
		bool getUphillGradient() const { return _walk.uphillGradient; }
		void setAvgStepsCount(long N);
		void resetAvgStep();
		void resetMomentum() { _walk.momentum=0; }
		void setChisqSignature(string sig) { _IOcontrol.chisqSignature=sig; }
		const string getChisqSignature() const { return _IOcontrol.chisqSignature; }
		/*!
			\brief set the initial step size factor for all parameters
			\details 
			@param size - initial maximal value of the step size in units of the size of the 
			parameter space defined by the parameter prior.

			default value: 0.5

			This controls how large the step size will be for every parameter during the burn-in stage.

			By default size=0.5 corresponds to the half of the domain size over which the parameter has been defined.
			The MC step is drawn from Boltzmann distribution.
			This means that the expected initial step size for that parameter will be moderated by the system temperature.
			1/_initialMaximalEnergy (kT) of the domain size.
			
			The full domain size is covered by a scaled boltzmann PDF and sampled within 
			range 0 to _initialMaximalEnergy=10 by default 
			(which is numerically doable for the initial default temperature of 1000 K).
			This distribution is then scaled to the range [0,size]
			Since the expectation value for the boltzmann PDF is kT - the expected value of the step size will be 
			size/_initialMaximalEnergy=1/20 of the full domain size for any given parameter.

			By lowering the value of this parameter the movement of MCMC will more and more resemble bug rowing motion 
			rather than	jumping motion as in Gibbs sampling. 
			
			Generally it is a good idea to have this parameter small <<1 because
			once the burn-in period is finished the chain automatically switches to the found best fit solution, 
			cools down appropriately and continues from there. If the size parameter is large it is possible 
			that the system will not cool down sufficiently in that single switch to the best fit solution 
			(to amount of cooling is generally proportional function of deltaChisq) and it may be that the chain 
			will jump away from the best fit solution. It will try to find it again, but in practice it may not 
			return to the best fit solution and in the result the likelihood surface will be very poorly probed 
			around the best fit solution. By making this parameter small, the chain will not depart too much
			from the best fit solution after the burn-in period is finished. In fact, with this parameter small, 
			the having large burn-in periods (lengths) probably	doesn't make much sense either. 
			On the other hand having this parameter small, worsens the mixing of the chain, so it is advisable 
			to use large number	of chains in such case.
			
			Remember to set this parameter BEFORE adding a new parameter, as the parameter step generating PDFs 
			are generated using this information when the parameter	is added.
			
			An alternative solution is to use large size parameter during burn-in period, 
			and abruptly switch to smaller values after the burn-in period regenerating
			the step generating PDFs. Such behaviour is also possible and it is triggered by 
			setting setInitialWalkStepSize(size) to a non-zero value.
			By default this is not done.
			
			If you want to have randomly chosen parameter values during the burn-in state, i.e. 
			chosen from within the prior volume	and according to the prior PDF you should set this value to zero (0).
		
		
			\date Mar 11, 2011, 12:32:10 PM
			
			revision Feb 23, 2020, 5:19:09 PM
			
			If you use gradient descent then setting this to 1 should be OK
			
			\author Bartosz Lew
		*/
		void setInitialStepSize(double size) { 
			_walk.initialPDFstepCRfraction=size; 
			if (getNparam()>0 and size>0) {
				string s="Setting initial step size but "+msgs->toStr(getNparam())+"parameters have already been added. These settings will not alter the parameters that were already added";
				msgs->warning(s,High);
				addImplementationNote("");
				addImplementationNote("WARNING");
				addImplementationNote(s);
			}
		}
		
		double getInitialStepSize() const { return _walk.initialPDFstepCRfraction; }
		
		/*!
			\brief set initial and final system temperatures
			\details 
			@param initial - typical value 1000
			@param final - typical value 10; A value of 0 indicates that thae final temperature is not to be used
			but rather other criteria will decide whether the system has cooled enough to stop the chain.

			\date Jan 9, 2018, 12:34:52 PM
		*/
		void setTemperatures(double initial, double final) { _cooling.initialTemperature=initial; _cooling.finalTemperature=final; updateCoolingPDF();  }

		/*!
			\brief defines the step size during the MCMC walk
			\details 
			@param size
			
			If you use gradient descent then setting this to 1 should be OK. In this mode this variable
			effectively defines the step at which initially the derivatives are taken.
			During descent the actual sizes are adjusted by falling temperature.

			See also setInitialStepSize()
		
			\date Feb 23, 2020, 5:20:15 PM
		*/
		void setInitialWalkStepSize(double size) { _walk.initialPDFstepCRfractionAfterBurnIn=size; }
		
		void setWalkInfoOutputFrequency(long everyNstates);
		long getWalkInfoOutputFrequency() const { return _walk.userOutputFrequency;	}
		
		/*!
			\brief add a string about particular implementation 
			\details 
			@param implementation note
			
			These notes will be stored in the setup file during run for 
			a complete description of the run.
			You can paste here values of any extra technical parameters that you invented and use in the 
			the derived class that are needed to generate eg. covariance matrix for the model or
			a model for a given MCMC state.
			
			\see See also string defining the actual implementation of the chisq value
			
			\date Mar 4, 2011, 12:45:15 PM
			\author Bartosz Lew
		*/
		void addImplementationNote(string s) { _IOcontrol.implementationNote.append(s); }

		/*!
			\brief set the p-values before and after the burn-in period for accepting the worse fit as the next MCMC state
			\details 
			This function sets the p-values and derives the associated acceptance thresholds for the cooling distribution for the 
			burn-in stage and after burn-in states.
		
			\date Mar 14, 2011, 4:40:23 PM
			\author Bartosz Lew
		*/
		void setAcceptWorsePvalues(double bipv, double abipv);
		void setAcceptWorseUniformScheme(bool uniform);
		acceptWorseScheme getAcceptWorseScheme() const { if (_cooling.acceptWorseStateFromUniformDistr) return awSchemeUniform; return awSchemeBoltzmann; }
		
		/*!
			\brief factor by which a RN is multiplied and compared against reduced chisq value to decide whether to accept the worse state than the current one
			\details 
			@param factor for the burn-in period
			@param factor after burn-in period; the default -1 value means that the value will be set automatically to the best fit chisqPerDOF value after the burn-in 
			period is over.

			To decide if a new worse state is to be accepted for the setting setAcceptWorseUniformScheme(false) a boltzmann distribution is used
			defined upon range from 0 to _initialMaximalEnergy (=10 by default). A random number RN from this distribution is drawn and compared against
			chisqperDOF. if RN > chisqperDOF then the worse state will be accepted. These factors alter this condition to:
			if RN* factor > chisqperDOF then the worse state will be accepted. By default both factors for the burn-in period and after burn-in period are one,
			but this might be too small is some cases hence it is possible to alter them manually.
			These settings are only important for the setAcceptWorseUniformScheme(false) setting: i.e. for the boltzmann scheme of accepting/rejecting worse states than the current one.
			In such case settings setAcceptWorsePvalues(p1,p2) is ignored.
			
			Note that these settings have direct impact (combined with finalTemeperature) on how long the chain will last.
			The better fit you reach (or the better your data fit to the model) the higher probability that the accept-worse
			condition will be met and after reaching the finalTemperature and when the maximalRejectionsCount is reached the chain will continue
			after accepting the worse state, possibly until the chain reaches the maximal chain length.
			On the other hand when the data cannot be fitted well with any model, then RN*factor will more often be < chisqperDOF and 
			the worse state will more often be rejected and the maximal rejections count limit will be reached sooner.
		
			TODO: implement appropriate convergence condition !
		
			\date Mar 21, 2011, 12:23:02 PM
			\author Bartosz Lew
		*/
		void setAcceptWorseBoltzmannSchemeFactors(double bi, double abi=-1) { _cooling.acceptWorseBoltzmannFact=bi; _cooling.acceptWorseAfterBurnInBoltzmannFact=abi;  }
		
		/*!
			\brief abstract class called when walk finishes to save the the simulations for the best fit model
			\details 
			@param simNum - number of simulations to be generated for saving

			This method should be reimplemented in the derived class in order to save the
			simulations for the best fit model of the covariance matrix. These can be used
			eg for plotting the confidence ranges around the best fit model for a given 
			input dataset.
			
			\date Mar 16, 2011, 12:54:01 PM
			\author Bartosz Lew
		*/
		virtual void saveBestFitCovSimulations();
		void saveBestFit();
		/*!
			\brief save all MCMC run results
			\details 

			Save last states, 
			save best fit
			save 1d and 2d posteriors
			save best fit cov simulations (if the function was overloaded)
		
			\date May 29, 2017, 10:12:46 PM
		*/
		void saveResults();
		void setBestFitCovSimNum(long n) { _chisqData.bestFitNcovSimNum=n; }
	
		/*!
			\brief sets the maximal number of times the system can be cooled down.
			\details 
			@param n - number of forced coolings

			The system will be forcibly cooled down every time when the maximalRejectionsCount was reached if
			the system temperature is larger than the final temperature. 
			Additional cooling is possible every time the system finds a better fit of the model to the data
			than the last accepted state.
		
			\date Mar 18, 2011, 12:16:52 PM
			\author Bartosz Lew
		*/
		void setMaximalNumberOfForcedCoolings(long n) { _cooling.maximalNumberOfForcedCoolings=n; }


		/*!
			\brief append a new link to the chain
			\details 
			@param link - new link reference
			
			This operation does not update any of the object's internal variables.
			It only appends the link to the structure containing all MC links.
		
			\date May 24, 2017, 11:10:42 AM
		*/
		void append(MClink& link);

		
		void appendLinkForConvergenceTest(MClink& link);
		
		/*!
			\brief append the whole chain to this chain
			\details 
			@param chain - 
		
			This operation appends the accepted and rejected links to this chain and updates (if necessary) 
			the best fit link.

			\date May 24, 2017, 11:15:34 AM
		*/
		void append(cpedsMCMC& chain);
		/*!
			\brief configure this chain using data from chain
			\details 
			@param chain - cpedsMCMC from which all parameter space and configuration data are copied.
		
			\date Jun 13, 2017, 2:50:14 PM
		*/
		void setup(cpedsMCMC& chain);

		const chisqData_t& getChisqData() const {
			return _chisqData;
		}
		chisqData_t& chisqData() { return _chisqData; }
		
		const convergence_t& getConvergence() const {
			return _convergence;
		}
		
		const cooling_t& getCooling() const {
			return _cooling;
		}
		
		const IOcontrol_t& getIOcontrol() const {
			return _IOcontrol;
		}
		IOcontrol_t& IOcontrol() {
			return _IOcontrol;
		}
		
		const QList<MClink>& getCbestFitChain() const {
			return _MCbestFitChain;
		}
		
		const MCparameterSpace& getParameterSpace() const {
			return _parameterSpace;
		}
		
		MCparameterSpace& parameterSpace() {
			return _parameterSpace;
		}
		
		const walk_t& getWalk() const {
			return _walk;
		}
		
		const bestFitData_t& getBestFitData() const {
			return _chisqData.bestFitData;
		}
		
		double getCoolingRate() const {
			return _cooling.coolingRate;
		}
		
		double getFinalTemperature() const {
			return _cooling.finalTemperature;
		}
		
		double getInitialTemperature() const {
			return _cooling.initialTemperature;
		}
		
		/*!
			\brief Return the number of parameters currently in the parameter space.
			\details 
			@return
			
			This is different from dims() which returns the number of dimensions with which the object has been
			initialized
		
			\date Feb 1, 2018, 12:54:12 PM
		*/
		int getNparam() const {
			return _parameterSpace.size();
		}
		
		string getOutputDir() const { return _IOcontrol.runFilesDir; }
		
		double getTemperature() const {
			return _cooling.temperature;
		}
		
		void setTemperature(double temperature) {
			this->_cooling.temperature = temperature;
		}
		
		
		void setID(long id);
		
		cpedsMCMC& operator=(cpedsMCMC& rhs);
		
		//
		// info variables
		//
		cpedsMsgs* msgs; //!< standard object messenger. It is made public to provide a direct access to the object.

	
	
		/***************************************************************************************/		
		/***************************************************************************************/		
		/***************************************************************************************/		
		/***************************************************************************************/		
		/***************************************************************************************/		
		/***************************************************************************************/		
		/***************************************************************************************/		
		/***************************************************************************************/		
		/***************************************************************************************/		
	protected:
		/***************************************************************************************/
		// protected methods
		/***************************************************************************************/
		
		int paramij2idx(int i,int j);
		MClink& current() { return _walk.current; }
		MClink& nextCandidate() { return _walk.next; }
		MClink& avgStep() { return _walk.avgStep; }

		
		void saveCooling();
		void saveModel();
		void dumpAll(bool dumpNow=false);

		/*!
			\brief this function defines the actual cooling scheme
			\details 
			@param deltaX - is the curent chisq value decrease

			This function defines how the system should cool down at given decrease in chisq related
			with the model fitting the data better than the previous one.
			
			Generically we want the system to cool down from the _initialTemperature to _finalTemperature
			in steps that are proportional to the increases in the goodness of fit of the model to the data.
			
			The cooling function uses the current setting of the cooling speed factor called _coolingRate,
			which can assume values from 0 to 1 for slow and fast cooling respectively.
			
			\date Feb 1, 2011, 12:38:11 PM
			\author Bartosz Lew
		*/
		void coolDownLogLin(double deltaX);
		/*!
			\brief cool the system in a geometric sequence
			@param factor - factor by which the temperature will shrink
			\details 
			
			
			The current temperature is divided by 2^coolingRage
			For the typical coolingRate=1 this will cool system by 2.
			
			@return returns the temperature decrement by which the system cooled down.
		
			\date Jan 9, 2018, 12:57:49 PM
		*/
		double coolDownGeometric(double factor=2.0);
		void coolForcibly();
		void coolDownProp(double deltaX);

		/*!
			\brief saves the setup of the MCMC simulated annealling 
			\details 
		
			\date Feb 1, 2011, 12:51:43 PM
			\author Bartosz Lew
		*/
		void saveSetup();
		
		/*!
			\brief saves a MC chain into the file
			\details 
			@param
			@return
			Every row consists of dims()+2 entries. The entries from 0..dims()-1 correspond to the parameter values
			Column dims() gives the corresponding chisq and dims()+1 gives the derived likelihood
		
			\date Feb 1, 2011, 1:09:31 PM
			\author Bartosz Lew
		*/
		void saveChain(QList<MClink>& chain, string fname);
		
		/*!
			\brief save the current index of the starting link to file
			\details 
			Although the chain length should basically increase only due to the walk: i.e. at the end
			we also store all states that were rejected for some reason as they also probe the likelihood surface.
			These states are stored at the beginning of the chain, hence the location of the statring link changes.
			This function saves the current index of the starting link.
			
			I believe that this states recycling improves the performance.
			\date Mar 4, 2011, 8:56:42 AM
			\author Bartosz Lew
		*/
		void saveFirstLinkLocation();
		void saveCurrentLength();
		void saveCurrentTemperature();
		void printLink(const MClink& link, string comment="") const;

		void storeBestFit();
		
		void mkOutputDir();

		/*!
			\brief recalculate the accept-worse thresholds for the new system temperature based on the current p-values
			\details 
			This function should be called AFTER the new inverse CDF cooling distribution was calculated for the new temperature
			\date Mar 14, 2011, 5:54:07 PM
			\author Bartosz Lew
		*/
		void updateAcceptWorseThresholds();
		
		
		/*!
			\brief update the internal list of MC step sizes. Useful for debugging purposes
			\details 

			@return returns delta w.r.t previous step
			The returned link contains parameter differences from the last step to step-2.
			The index of the returned link is set to the index of the last step.
			
			This function assumes that the new link has been added, so it should be called
			after the link is added.
		
			\date Feb 12, 2020, 4:12:51 PM
		*/
		MClink update_deltas();
		
		/*!
			\brief extract parameter likelihoods from all accepted, rejected and test states
			\details 
			@param j - parameter index
			@return an unsorted function of parameter values (X) and their likelihoods (Y)
		
			\date May 24, 2017, 12:04:57 PM
		*/
		mscsFunction combineParamLikelihoods(int j);
		
		bool forcedCoolingPossible() { return _walk.rejectionsCount>=_cooling.maximalRejectionsCount and getTemperature()>getFinalTemperature(); }
		bool coolingPossible() { return getTemperature()>getFinalTemperature(); }
		
		/***************************************************************************************/
		// protected structures
		/***************************************************************************************/

		//
		// Markov Chain structures
		//
		
		cpedsMC _MCaccepted; //!< chain of all accepted links (including those that yield a higher chisq)
		cpedsMC _MCaccepted_deltas; /*!< TODO: chain of steps for accepted links 
								(including those that yield a higher chisq).
								There is a functionaliy overlap between this structure
								and _walk.steps structure which holds exactly the same 
								information */
		cpedsMC _MCrej; //!< list of rejected states
		cpedsMC _MCtest; //!< list of tested states
		MClink _bestFit; //!< stores the current best fit state
		cpedsMC _MCbestFitChain; //!< chain to hold the history of best fit links 
		MClink _startingLink; //!< the initial link
		
		//
		// variables controlling the annealing and cooling
		//
		cooling_t _cooling;
		
		//
		// variables controlling the walk in parameter space
		//
		
		
		MCparameterSpace _parameterSpace; //!< stores information on parameter space - ranges, and or grid, and parameter priors
		walk_t _walk;
		
		//
		// data used to calculate the chisq
		//
		chisqData_t _chisqData;
		
		//
		// run/debug/save control stuff
		//
		IOcontrol_t _IOcontrol;

		
		//
		// variables controlling the convergence -- NOT USED AT THE MOMENT
		//
		convergence_t _convergence;
		
		/***************************************************************************************/		
		/***************************************************************************************/		
		/***************************************************************************************/		
		/***************************************************************************************/		
		/***************************************************************************************/		
		/***************************************************************************************/		
		/***************************************************************************************/		
		/***************************************************************************************/		
		/***************************************************************************************/		
	private:
		/***************************************************************************************/
		// private methods
		/***************************************************************************************/
		
		/*!
			\brief define the starting point in the parameter space according to the given priors defined in the parameter space
			\details 
			@return N-element array of values corresponding to the parameters values to start the walk with, where N is the number of parameters
		
			\date Nov 22, 2010, 4:26:48 PM
			\author Bartosz Lew
		*/
		double* getStartingPoint();
		
		MClink getNextPoint(const MClink& current);
		MClink getNextPointUphillGradient(const MClink& current);
		
		double getSign();
		void generateInitialStepSizePDFsAfterBurnIn();
		
		//! returns the chisq(i-1)-chisq(i); i should be >=1; The positive value means the fit improves
		double deltaChisq(long i) { //printf("%lE %lE --- %lE %lE\n",chain(i-1).chisq(), chain(i).chisq(),   chain(i-1)[0],  chain(i)[0]); 
			if (i>=1)
				return chain(i-1).chisq()-chain(i).chisq(); 
			return chain(i).chisq(); 
		}
//		double deltaChisq(long i) { printf("%lE %lE\n",_MCchain[i-1].chisq(), _MCchain[i].chisq()); return _MCchain[i-1].chisq()-_MCchain[i].chisq(); }

		void acceptWorseState(bool acceptWorse, MClink& next, long& rejectionsCount);
		
		/***************************************************************************************/
		// private structures
		/***************************************************************************************/
		cpedsList<double> _acceptWorseRNs;
};

#endif /* CPEDSMCMC_H_ */










