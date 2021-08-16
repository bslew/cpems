/*!
 * \file cpeds_calculate_cov.cpp -
 * 	calculates covariance matrix from a file - suitable for
 * arbitrary data sizes
 *
 *  Project: cpems
 *  Created on: Aug 13, 2021, 2:50:18 PM
 *  Author: blew
 */

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <cpems/Mscs-function.h>
#include <cpems/Mscs-function3dregc.h>
//#include <yaml-cpp/yaml.h>


#define PARAMETER_FILE "input.par"

using namespace std;


spdlog::logger getLogger(int verbosity);
string getCmdString(int argc, char** argv);
string getProgramVersionString();
boost::program_options::variables_map  parseOptions(int argc, char** argv);



int test_direction(boost::program_options::variables_map &opt,spdlog::logger& logger) {
	logger.info("dt [UTC]: {}",opt["dtSt"].as<string>());
	return 0;
	
}

/* ******************************************************************************************** */
/* ******************************************************************************************** */
/* ******************************************************************************************** */
int main(int argc, char **argv) {
	
	// initialize parser
	boost::program_options::variables_map opt=parseOptions(argc,argv);
	
	spdlog::logger logger=getLogger(opt["verbosity"].as<int>());
	logger.info("cpeds_calculate_cov - NEW RUN");
	logger.info(getCmdString(argc,argv));
	
	std::string infile;
//	std::cout << opt["input_file"].as< vector<string> >() << std::endl;
	if (opt.count("input-file")) {
		vector<string> input_files=opt["input-file"].as< vector<string> >();
		if (input_files.size()>=1) infile=input_files[0];
	}
	
	std::ifstream ifs(infile);
	//
	// save yaml config file for the dataset
	//
	//	YAML::Node config; 
	//	config["orig_ns"]=opt["ns"].as<long>();
	//	config["reg_ns"]=opt["nsc"].as<long>();
	//	config["hmin"]=opt["hmin"].as<double>();
	//	for (auto ch : channels) {
	//		config["channels"].push_back(ch);
	//	}
	//	std::ofstream fout(opt["odir"].as<string>()+"/config.yaml");
	//	fout << config;
	//	fout.close();
	
	
	if (opt["test"].as<bool>()) {
		auto ret=test_direction(opt,logger);
		exit(ret);
	}
	
	//
	// iterate over all input files in parallel
	//
	//	vector<string>::const_iterator fIt;
	long i=0;
	//#pragma omp parallel for private(i) firstprivate(channels) 
	//	for (auto infile : input_files) {
	//	for (fIt=input_files.begin(); fIt!=input_files.end();++fIt) {
	//	for (i=0;i<input_files.size();i++) {
	//		string infile=input_files[i];
	//    	logger.info("Processing file: {}", infile);
	//    	logger.debug("loading data...");
	
	//	}
	
	
	// prepare data for sasving
	
//	mscsFunction3dregc sun_data;
//	sun_data.setSize(4,sun_elev_avg.pointsCount());
//	sun_data.allocFunctionSpace();
	
//	for (long i=0;i<sun_elev_avg.pointsCount();i++) {
//		sun_data(0,i)=sun_elev_avg.getx(i);
//		sun_data(1,i)=sun_elev_avg.f(i);
//		sun_data(2,i)=sin(sun_elev_avg.f(i)*PI180);
//		sun_data(3,i)=sun_dist.f(i);
//	}
//	sun_data.setSavePrecision(8);
//	sun_data.saveSlice(2,0,opt["ofile"].as<string>());
	
//	sun_elev_avg.save(opt["ofile"].as<string>());
	//	DirectionAh ahSunCmd(opt["sunA"].as<double>(),opt["sunh"].as<double>());
	
	return 0;
}

/* ******************************************************************************************** */
boost::program_options::variables_map parseOptions(int argc, char** argv) {
	namespace po = boost::program_options;
	po::variables_map vm;
	
	
	try {
		int opt=1;
		bool boolopt;
		long longopt;
		double dbl;
		string stropt;
		vector<string> svec;
		stringstream ss;
		ss << std::getenv("HOME") << "/.config/cpems/config.txt";
		string config_file=ss.str();
		boost::filesystem::path cwd=boost::filesystem::current_path();
		//        string parameter_file=cwd.string()+"/input.par";
		string parameter_file="input.par";
		
		// Declare a group of options that will be 
		// allowed only on command line
		po::options_description generic("Generic options");
		generic.add_options()
            		("version,V", "print version string")
					//			("ifile,i", po::value<string>(&stropt)->multitoken()->required(), "input file name")
					//		    ("version,V", po::bool_switch()->default_value(false), "print program version and exit")
					("test", po::bool_switch()->default_value(false), "test positions")
					("help", "produce help message")
					("config,c", po::value<string>(&config_file)->default_value(config_file),
							"name of a global ACClass configuration file.")
							("param,p", po::value<string>(&parameter_file)->default_value(parameter_file),
									"parameter file")
									;
		
		// Declare a group of options that will be 
		// allowed both on command line and in
		// config file
		po::options_description config("Configuration");
		config.add_options()
            		("verbosity,v", po::value<int>(&opt)->default_value(2), 
            				"verbosity level")
					//            ("include-path,I", 
					//                 po::value< vector<string> >()->composing(), 
					//                 "include path")
					("ofile,o", po::value<string>()->default_value("cov.txt"), "output file name")
					("MJD", po::value<bool>()->default_value(false), "output MJD instead of JD")
					("dtEn", po::value<string>()->default_value("2021-01-01T10:00:00"), "Ending UTC date and time for generating "
							"list of sun positions. (eg. 2000-01-01T10:00:00)")
					("Navg", po::value<int>()->default_value(144), ""
							"Averaging scale - number of calculations"
							"that should be averaged for saving "
							"(eg. dt=10 and avg=144 means that position"
							"is calculated every 10 minutes and every 144"
							"positions aver averaged (1-day) before saving)")
					;
		
		// Hidden options, will be allowed both on command line and
		// in config file, but will not be shown to the user.
		po::options_description hidden("Hidden options");
		hidden.add_options()
		            ("input-file", po::value< vector<string> >(), "input files")
            		;
		
		
		po::options_description cmdline_options;
		cmdline_options.add(generic).add(config).add(hidden);
		
		po::options_description config_file_options;
		config_file_options.add(config).add(hidden);
		
		po::options_description parameter_file_options;
		parameter_file_options.add(config).add(hidden);
		
		po::options_description visible("cpeds_calculate_cov\n\n "
				"This program calculates covariance matrix from input data file."
				"\n\n. "
				"Each row in the input file should hold one full observation of all variables."
				"Each row should contain the same number of columns."
				"\n"
				"example usage: ${file_base} [input_files]"
				"\n"
				"\n\n"
				"Allowed options");
		visible.add(generic).add(config);
		
		po::positional_options_description p;
		p.add("input-file", -1);
		
		store(po::command_line_parser(argc, argv).
				options(cmdline_options).positional(p).run(), vm);
		notify(vm);
		
		// process config file
		ifstream ifs(config_file.c_str());
		if (!ifs) {
			cout << "cannot open config file: " << config_file << "\n";
			//            return 0;
		}
		else {
			store(parse_config_file(ifs, config_file_options), vm);
			notify(vm);
		}
		
		
		// process parameter file 
		ifs.open(parameter_file.c_str());
		if (!ifs.is_open()) {
			cout << "cannot open parameter file: " << parameter_file << "\n";
			//            return 0;
		}
		else {
			store(parse_config_file(ifs, parameter_file_options), vm);
			notify(vm);
		}
		
		
		if (vm.count("help")) {
			cout << visible << "\n";
			exit(0);
		}
		
		//    	if (vm["version"].as<bool>()) {
		//    		std::cout << getProgramVersionString() << std::endl;
		//    		exit(0);
		//    	}
		
		if (vm.count("version")) {
			std::cout << getProgramVersionString() << std::endl;
			exit(0);
		}
		
		/*
        if (vm.count("input-file")) {
            cout << "Input files are: \n" ;
            for ( auto& s : vm["input-file"].as< vector<string> >() ) {
            	cout << s << " ";
            }
            cout << "\n";
        }
        else {
        	cout << "No input files provided" << std::endl;
        	exit(0);
        }
		 */
		
		cout << "Verbosity level is " << opt << "\n";                
	}
	catch(exception& e)
	{
		cout << e.what() << "\n";
		exit(1);
	}    
	
	return vm;
	
}
/* ******************************************************************************************** */
string getProgramVersionString() {
	string rev;
	rev="0.1.1";
	
#ifdef GOMPI
	rev+=" (MPI)";
#endif
	
	return rev;
}
/* ******************************************************************************************** */
spdlog::logger getLogger(int verbosity) {
	auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
	console_sink->set_level(spdlog::level::info);
	//    console_sink->set_level(spdlog::level::debug);
	//    console_sink->set_pattern("[sun_elevation] [%^%l%$] %v");
	
	auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("cpeds_calculate_cov.log", true);
	file_sink->set_level(spdlog::level::debug);
	
	spdlog::logger logger("cpeds_calculate_cov", {console_sink, file_sink});
	logger.set_level(spdlog::level::info); 
	
	if (verbosity>2) {
		//		cout << "setting verbosity " << opt["verbosity"].as<int>() << "\n";
		logger.sinks()[0]->set_level(spdlog::level::debug);
		logger.sinks()[1]->set_level(spdlog::level::debug);
		logger.set_level(spdlog::level::debug);
	}
	
	return logger;
}
/* ******************************************************************************************** */
string getCmdString(int argc, char** argv) {
	stringstream ss;
	for (long i = 0; i < argc; i++) {
		ss << string(argv[i]) << " ";
	}
	return ss.str();
}
/* ******************************************************************************************** */
