#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <cpgplot.h>
#include <math.h>
#include <cfitsio/fitsio.h>
#include <string.h>
#include <time.h>
#include <tclap/CmdLine.h>
#include "cpeds-consts.h"
#include "cpeds-math.h"
#include "Mscs-map.h"
#include "Mscs-global-defs.h"
#include "mscsMapINC.h"
#include "MscsWMAPdata.h"
#include "Mscs-WMAPsimulation.h"
#include "Mscs-angular_power_spectrum.h"

#ifdef GOMPI
#include "mpi.h"
#endif

//void print_program_usage();
#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
using namespace TCLAP;
#define STD std
#else
#define STD
#endif

//declaration of the global variables
void parseOptions(int argc, char** argv);
long _sim_st,_sim_en, _ns, _lmax, _procNum, _centralLimitNumbers;
double _ClcalibVal, _precalibrate, _calibrate;
string  _Clcalib, _mask_file, _Cl, _DA, _yr, _outFilePref, _outFileSuff,_outDir, _almsFormat, _rng;
bool  _usePixtf, _savePartial,_precalib,_calib, _masked, _saveAlms, _Clconv;

//string*_input_file=NULL,_mask_file,_ft,_ord,_multi_mask,_op,_outfile,_sim_dir,_sim_file_pref,_sim_file_post;
//long _idnum,_lreg_num=0,_breg_num=0,_mmnside,_cmd_files_num,_files_num,_files_load_num,_sim_files_num;
//bool masked=false,multi_masked=false, _histogram, _two_files, _one_file,_mask_after,_dont_mask_before,_numeric,_invert,_sim,_mask_from_here,_dont_merge,_perform_under_mask,_predivide_by_var,_predivide_by_rms,_abs;
//double _value,_value2,_lowR,_hiR;
//-----------------------------------


int main(int argc, char **argv) {

	cpedsMsgs msgs("mk_gauss_map_simWMAP");
	msgs.setSaveRunWriteMode('a');
	msgs.saveThisRun(argc,argv);
	//----------------------------------------------------------------------------------------------------
	//PARSE THE COMMAND LINE ARGUMENTS
	//----------------------------------------------------------------------------------------------------
	parseOptions(argc,argv);
	Mscs_initiate_global_variables();
	//----------------------------------------------------------------------------------------------------

	cpedsRNG *rns;
	mscsAngularPowerSpectrum Cl;
	mscsMap *mask;
	mscsWMAPsimulation *sim;
	mscsWMAPspecifications wmapSpecs;
	mscsWMAPspecifications::WMAPversion wmapVersion=wmapSpecs.getVersion(_yr);
	double ClCalibrationFactor=1, ClunitConversionFactor=_ClcalibVal;
	int node_num=1, world_id=0;
	long Nsim, sim_st,sim_en;
	long Nsimpn;


#ifdef GOMPI
//	MPI_Status mpi_stat;
	int ierr;
//	long Ndirpn;
//	double from_slave, *from_slave1,*from_slave3,*from_slave5,*from_slave6,*from_slave7,*from_slave9;
#endif


	//------------------------------------------------------------------------------
	//
	// prepare to fork simulation range
	//

	sim_st=_sim_st;
	sim_en=_sim_en;
#ifdef GOMPI
	ierr = MPI_Init ( &argc, &argv ); // initiate MPI
	ierr = MPI_Comm_size ( MPI_COMM_WORLD, &node_num ); // Get the number of processes
	ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &world_id ); // Get the individual process ID

  // initiate the runs
	Nsim=_sim_en-_sim_st+1;
	Nsimpn=Nsim/node_num;
	sim_st=world_id*Nsimpn+_sim_st;
	sim_en=(world_id+1)*Nsimpn+_sim_st-1;
	if ( world_id == node_num-1 ) { sim_en=Nsim; } // make sure all simulations will be processed
#endif


	//
	// prepare output directories
	//
	string confFileName=_outDir+"/"+"gaussianSimulation.conf";
	string cmd="mkdir -p "+_outDir;
	system(cmd.c_str());


	//
	// initialize the shared memory space for these runs and other watching processes
	//

	// create shared memory block
	int fd;
	fd = shm_open(MSCS_GLOBAL__RUN_SIMULATION_SHM_KEY, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
	if (fd == -1) { msgs.criticalError("shared memory: shm_open failed",High); }
	if (ftruncate(fd, sizeof(Mscs_run_control_t)) == -1) { msgs.criticalError("shared memory: adjusting shared memory size failed",High); }


	/* Map shared memory object */
	MSCS_GLOBAL__RUN_SIMULATION_INFO= (Mscs_run_control_t*)mmap(NULL, sizeof(Mscs_run_control_t), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (MSCS_GLOBAL__RUN_SIMULATION_INFO == MAP_FAILED) { msgs.criticalError("shared memory: mapping failed",High); }


    // set the common run data
    string tmps;
#ifdef GOMPI
    if (world_id==0) {
    	MSCS_GLOBAL__RUN_SIMULATION_INFO->procNum=node_num;
    	MSCS_GLOBAL__RUN_SIMULATION_INFO->st_thread=0;
    	MSCS_GLOBAL__RUN_SIMULATION_INFO->en_thread=node_num-1;
    	tmps=confFileName+"__thread_";
    	strcpy(MSCS_GLOBAL__RUN_SIMULATION_INFO->sim_run_file_pref,tmps.c_str());
    	strcpy(MSCS_GLOBAL__RUN_SIMULATION_INFO->sim_run_file_suff,".run");
    }
#else
	MSCS_GLOBAL__RUN_SIMULATION_INFO->procNum=1;
	MSCS_GLOBAL__RUN_SIMULATION_INFO->st_thread=-1; // this is to distinguish from MPI run. for a single program run the file name doesn't contain any id numbers
	MSCS_GLOBAL__RUN_SIMULATION_INFO->en_thread=-1;
	tmps=confFileName;
	strcpy(MSCS_GLOBAL__RUN_SIMULATION_INFO->sim_run_file_pref,tmps.c_str());
	strcpy(MSCS_GLOBAL__RUN_SIMULATION_INFO->sim_run_file_suff,".run");
#endif

	//
	// define power spectrum
	//
	if (Cl.load(_Cl,true)!=cpedsSuccess) { msgs.error("no Cl file",High); return -1; }
//	Cl.print();
	if (Cl.get_lmax() < _lmax) msgs.warning("The requested lmax is larger than the largest l in the power spectrum file. Will assume zero power beyond maximal loaded multipole.",High);
	if (_Clconv) Cl.divide_llpotwoPI();
	if (_Clcalib != "nocalib") {
		if (_Clcalib=="Tcmb") { ClunitConversionFactor=double(T_CBR0*T_CBR0); }
		if (_Clcalib=="TcmbWMAP") {
			if (_yr=="WMAP1yr" or _yr=="WMAP3yrs" or _yr=="WMAP5yrs" or _yr=="WMAP7yrs") {
				ClCalibrationFactor=wmapSpecs.get_Cl_calibration(wmapVersion);
				ClunitConversionFactor=double(T_CBR0*T_CBR0);
			}
		}
	}
	Cl*=double(ClCalibrationFactor*ClunitConversionFactor);


	// TODO the error about cannot load pixtf should not exist there since we don't want to load pixtf from the cmd line.

	//
	// define the pixel transfer function
	//
	string pixtfType;
	if (_usePixtf) { pixtfType="healpix"; } else { pixtfType=""; }


	//
	// start random number generators
	//
	if (_rng=="gaussian")
		rns = new cpedsRNG("gaussian"); // TODO why the gaussian numbers via central limit theorem do a worse job than this method
	if (_rng=="gaussian_invcdf") {
		rns = new cpedsRNG("gaussian_invcdf");
		rns->setCentralLimitNumbers(_centralLimitNumbers);
	}
	else {
		rns = new cpedsRNG(_rng);
//		rns->setCentralLimitNumbers(_centralLimitNumbers);		
	}
	rns->seed(0);
	rns->seedOffset(world_id);


	//
	// load mask
	//
	if (_masked) { mask = new mscsMap("mask"); mask->loadfits(_mask_file,MSCS_WMAP_MASK_FITS_FILE_COLUMN_NAME,"m"); }
	else { mask=NULL; }


	//
	// gather run configuration information
	//
#ifdef GOMPI
	msgs.setLogFileName(confFileName+"__thread_"+msgs.toStr(long(world_id)));
#else
	msgs.setLogFileName(confFileName);
#endif
	msgs.loggingOn();
	msgs.clearLogFile();
	msgs.say("SIMULATION CONFIGURATION: ",Top);
	msgs.say("Simulation type: "+_DA,High);
	msgs.say("Simulated data release: "+_yr,High);

	msgs.say("Starting simulation: "+msgs.toStr(_sim_st),High);
	msgs.say("Ending simulation: "+msgs.toStr(_sim_en),High);
	msgs.say("Simulation nside: "+msgs.toStr(_ns),High);
	msgs.say("Simulation lmax: "+msgs.toStr(_lmax),High);
	msgs.say("Central limit numbers: "+msgs.toStr(_centralLimitNumbers),High);
	msgs.say("use pixel tf: "+msgs.toStr(_usePixtf),High);
	msgs.say("pixel tf: "+pixtfType,High);
	msgs.say("Cl file name: "+_Cl,High);
	msgs.say("Cl backup (calibrated) file name: angular_power_spectrum.cl",High);
#ifdef GOMPI
	if (world_id==0) Cl.save(_outDir+"/angular_power_spectrum.cl");
#endif
#ifndef GOMPI
	Cl.save(_outDir+"/angular_power_spectrum.cl");
#endif
	msgs.say("Cl calibration factor: "+msgs.toStr(ClCalibrationFactor),High);
	msgs.say("Cl unit conversion factor to [K]: "+msgs.toStr(ClunitConversionFactor),High);
	msgs.say("Using mask (for mean temperature calibrations): "+msgs.toStr(_masked)+" ("+_mask_file+")",High);
	msgs.say("Precalibration of mean is: "+msgs.toStr(_precalib)+ ", precalibration value ("+msgs.toStr(_precalibrate)+")",High);
	msgs.say("Calibration of mean is: "+msgs.toStr(_calib)+ ", calibration value ("+msgs.toStr(_calibrate)+")",High);

	msgs.say("FILES CONFIGURATION: ",Top);
	msgs.say("Simulation files prefix: "+_outFilePref,High);
	msgs.say("Simulation files suffix: "+_outFileSuff,High);
	msgs.say("Saving partial simulations: "+msgs.toStr(_savePartial),High);
	msgs.say("Simulation run working directory: "+_outDir,High);
	msgs.say("Signal alms files save: "+msgs.toStr(_saveAlms),High);
	msgs.say("Signal alms files format: "+_almsFormat,High);

	msgs.say("THREADS CONFIGURATION: ",Top);
	msgs.say("Number of threads: "+msgs.toStr(node_num),High);
	msgs.say("This thread id: "+msgs.toStr(world_id),High);
	msgs.say("This thread starting simulation: "+msgs.toStr(sim_st),High);
	msgs.say("This thread ending simulation: "+msgs.toStr(sim_en),High);
	msgs.say("This thread random seed: "+msgs.toStr(rns->seed()),High);
	msgs.say("This thread random seed offset: "+msgs.toStr(rns->seedOffset()),High);


	msgs.loggingOff();








	//
	// start simulations
	//

	for (int i=sim_st;i<=sim_en;i++) {
		// record current status
#ifdef GOMPI
		msgs.setLogFileName(confFileName+"__thread_"+msgs.toStr(long(world_id))+".run");
#else
		msgs.setLogFileName(confFileName+".run");
#endif
		msgs.loggingOn();
		msgs.clearLogFile();
		msgs.print(msgs.toStr(long(sim_en-sim_st+1))+" "+msgs.toStr(sim_st)+" "+msgs.toStr(sim_en)+" "+msgs.toStr(i),High);
		msgs.loggingOff();

		// generate simulation
		sim = new mscsWMAPsimulation(wmapSpecs.getDA(_DA),wmapSpecs.getVersion(_yr),_ns,_lmax,Cl,pixtfType,rns,_precalib, _precalibrate,_calib, _calibrate,mask);
		sim->setNames(_outDir,_outFilePref,i,_outFileSuff);
		sim->setSavePartialSimulations(_savePartial);
		sim->makeWMAPsimulation();

		// save simulation(s)
		sim->saveSimulation(_outDir,_outFilePref,i,_outFileSuff,_saveAlms);
		delete sim;
	}
	delete rns;
	if (mask!=NULL) delete mask;

#ifdef GOMPI
	ierr = MPI_Finalize ( ); // Shut down MPI.
#endif

	shm_unlink(MSCS_GLOBAL__RUN_SIMULATION_SHM_KEY);
	return 0;
}

/***************************************************************************************/
void parseOptions(int argc, char** argv) {
//	long i;
	vector<string> allowedStr;
	ValuesConstraint<string>* allowedClCal;
	ValuesConstraint<string>* allowedWMAPda;
	ValuesConstraint<string>* allowedWMAPyr;
	ValuesConstraint<string>* allowedRNG;
	ValuesConstraint<string>* allowedAmlsFmt;

	/*   try { */

	CmdLine cmd("mk_gauss_map_simWMAP\n generates gaussian simulations of the WMAP data.", ' ', Mscs_version.c_str() );

	//	UnlabeledMultiArg<string> infiles("infiles", "input file names (<100)","", "string");     	cmd.add( infiles );
	ValueArg<long> sim_st("","st","numer of simulation to start from (default: 1)",false,1,"long"); cmd.add(sim_st);
	ValueArg<long> sim_en("","en","numer of simulation to end at (default: 1000)",false,1000,"long"); cmd.add(sim_en);
	ValueArg<string> mask("m","mask"," mask file (prefix)",false,"nomask","string"); cmd.add(mask);

	ValueArg<double> precalibrate("","precalib","whether or not to shift mean of the maps to requested value before making INC (default: false)",false,0,"double"); cmd.add(precalibrate);
	ValueArg<double> calibrate("","calib","whether or not to shift mean of the maps to requested value after making INC (default: false)",false,0,"double"); cmd.add(calibrate);

	ValueArg<string> Cl("","cl","Theoretical CMB angular power spectrum to use for map genertion",true,"","string"); cmd.add(Cl);
	allowedStr.clear();
	allowedStr.push_back("nocalib");
	allowedStr.push_back("Tcmb");
	allowedStr.push_back("TcmbWMAP");
	allowedStr.push_back("other");
	allowedClCal = new ValuesConstraint<string>(allowedStr);
	ValueArg<string> Clcalib("","clcal","multiply theoretical cl by a calibration factor",false,"nocalib",allowedClCal); cmd.add(Clcalib);
	ValueArg<double> ClcalibVal("","cal","calibration factor. Cl will be multiplied by this value of clcal \"other\" was used",false,0,"double"); cmd.add(ClcalibVal);
	SwitchArg Clconv("","clconv", "whether or not to divide the input Cl by l(l+1)/twoPI factor before generating alms, default: false)", false);	cmd.add(Clconv);

	allowedStr.clear();
	allowedStr.push_back("K1");
	allowedStr.push_back("K2");
	allowedStr.push_back("Q1");
	allowedStr.push_back("Q2");
	allowedStr.push_back("V1");
	allowedStr.push_back("V2");
	allowedStr.push_back("W1");
	allowedStr.push_back("W2");
	allowedStr.push_back("W3");
	allowedStr.push_back("W4");
	allowedStr.push_back("Qinc");
	allowedStr.push_back("Vinc");
	allowedStr.push_back("Winc");
	allowedStr.push_back("QVWinc");
	allowedStr.push_back("Q_V_W_QVW_inc");
	allowedStr.push_back("VWinc");
	allowedStr.push_back("V_W_VW_inc");
	allowedWMAPda = new ValuesConstraint<string>(allowedStr);	
	ValueArg<string> DA("","DA","defines the type of map that will be generated",false,"VWinc",allowedWMAPda); cmd.add(DA);
	allowedStr.clear();
	allowedStr.push_back("WMAP1yr");
	allowedStr.push_back("WMAP3yrs");
	allowedStr.push_back("WMAP5yrs");
	allowedStr.push_back("WMAP7yrs");
	allowedWMAPyr = new ValuesConstraint<string>(allowedStr);	
	ValueArg<string> yr("","yr","defines the release of the WMAP data",false,"7yrs",allowedWMAPyr); cmd.add(yr);
	ValueArg<long> ns("","ns","Healpix resolution parameter (default: 512)",false,512,"long"); cmd.add(ns);
	ValueArg<long> lmax("","lmax","maximal multipole number (default: 1024)",false,1024,"long"); cmd.add(lmax);
	SwitchArg usePixtf("","usePixtf", "whether or not to use pixel transfer functions, default: false)", false);	cmd.add(usePixtf);

	allowedStr.clear();
	allowedStr.push_back("gaussian");
	allowedStr.push_back("gaussian_invcdf");
	allowedStr.push_back("gaussian_circle");
	allowedRNG= new ValuesConstraint<string>(allowedStr);	
	ValueArg<string> rng("","rng"," random number generator; default is 'gaussian' which is the generator based on the central limit theorem where the number of the uniformly distributed numbers is controlled by centralLimitNumbers option. 'gaussian_invcdf' is based on the inverting gaussian cdf method where cdf is generated within +-5.5 sigma; this value can be controlled in cpeds library (see cpeds_generate_gaussian_distribuant_function function in cpeds-math.c). Both methods use by default the uniform rng called 'gsl_rng_mt19937'  (default: gaussian_invcdf)",false,"gaussian_invcdf",allowedRNG); cmd.add(rng);
	ValueArg<long> centralLimitNumbers("","cln","number of numbers to use for generation of gaussian random numbers using central limit theorem (default: 20)",false,20,"long"); cmd.add(centralLimitNumbers);


	ValueArg<long> procNum("","proc","numer of parallel threads to start (default: 1)",false,1,"long"); cmd.add(procNum);

	SwitchArg savePartial("","save_partial", "whether or not to save the partial maps (eg. Q1, Q2 in case of Qinc request), default: false)", false);	cmd.add(savePartial);
	SwitchArg saveAlms("","save_alms", "whether or not to save the sigmal cmb alms, default: false)", false);	cmd.add(saveAlms);
	allowedStr.clear();
	allowedStr.push_back("RI");
	// allowedStr.push_back("lmRI");
	// allowedStr.push_back("lmRIMP");
	// allowedStr.push_back("lmMP");
	allowedAmlsFmt= new ValuesConstraint<string>(allowedStr);
	ValueArg<string> almsFormat("","almsFormat","defines the format for the saved sigmal alms (default: RI)",false,"RI",allowedAmlsFmt); cmd.add(almsFormat);

	ValueArg<string> outFilePref("","outPref","output file prefix (default: WMAPsim)",false,"WMAPsim","string"); cmd.add(outFilePref);
	ValueArg<string> outFileSuff("","outSuff","output file suffix",false,"","string"); cmd.add(outFileSuff);
	ValueArg<string> outDir("","outDir","output main directory (default: current dir)",false,"","string"); cmd.add(outDir);






	//
	// Parse the command line.
	//
	cmd.parse(argc,argv);
	//
	// Set variables
	//
	_sim_st = sim_st.getValue();
	_sim_en = sim_en.getValue();
	_mask_file = mask.getValue(); if (!mask.isSet()) { _masked = false; } else {_masked = true; }

	_precalibrate = precalibrate.getValue(); _precalib=precalibrate.isSet();
	_calibrate = calibrate.getValue(); _calib=calibrate.isSet();

	_Cl = Cl.getValue();
	_Clconv = Clconv.getValue();
	_Clcalib = Clcalib.getValue();
	_ClcalibVal = ClcalibVal.getValue();

	_DA = DA.getValue();
	_yr = yr.getValue();

	_ns = ns.getValue();
	_lmax = lmax.getValue();

	_usePixtf = usePixtf.getValue();
	_rng = rng.getValue();
	_centralLimitNumbers =  centralLimitNumbers.getValue();

	_procNum = procNum.getValue();

	_savePartial=savePartial.getValue();
	_saveAlms=saveAlms.getValue();
	_almsFormat=almsFormat.getValue();
	_outFilePref= outFilePref.getValue();
	_outFileSuff= outFileSuff.getValue();
	_outDir= outDir.getValue();


	delete allowedAmlsFmt;
	delete allowedRNG;
	delete allowedWMAPda;
	delete allowedWMAPyr;
	delete allowedClCal;

	/*   } catch ( ArgException& e ) */
	/*       { cout << "ERROR: " << e.error() << " " << e.argId() << endl; } */
}




