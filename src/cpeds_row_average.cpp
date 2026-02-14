/*!
 * \file cpeds_row_average.cpp - reads stdin and calculates line by line averages and outputs to stdout
 *
 *  Project: cpems
 *  Created on: Aug 25, 2021 11:36:59 AM
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
#include "StringParser.h"
//#include 
//#include <yaml-cpp/yaml.h>


#define PARAMETER_FILE "input.par"

using namespace std;


spdlog::logger getLogger(int verbosity);
string getCmdString(int argc, char** argv);
string getProgramVersionString();
boost::program_options::variables_map  parseOptions(int argc, char** argv);

template<typename T> T vec_mean(std::vector<T> &v,long st,long en) {
	T m=0;
//	std::vector<T>::iterator it;
//	for (auto it=v.begin();it<v.end(); ++it) {
//		m+=*it;
//	}
//	return m/v.size();
	
	for (auto it=st;it<en; ++it) {
		m+=v[it];
	}
	return m/(en-st);
}

/* ******************************************************************************************** */
/* ******************************************************************************************** */
/* ******************************************************************************************** */
int main(int argc, char **argv) {
	
	// initialize parser
	boost::program_options::variables_map opt=parseOptions(argc,argv);

	spdlog::logger logger=getLogger(opt["verbosity"].as<int>());
//	logger.info("cpeds_row_average - NEW RUN");
//	logger.info(getCmdString(argc,argv));

//	std::cout << opt["input_file"].as< vector<string> >() << std::endl;
//	if (opt.count("input-file")) {
//	vector<string> input_files=opt["input-file"].as< vector<string> >();
		
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
	
	
	long n=opt["n"].as<long>();
	for (std::string line; std::getline(std::cin, line);) {
			cpeds::StringParser p;
			auto v=p.splitAs<double>(line, ' ');
			long ist=0;
			long ien=ist+n;
			if (ien>v.size()) ien=v.size();
			while (ist<v.size() and ist<ien-1) {
				std::cout << vec_mean(v,ist,ien) << " ";
				ist+=n;
				ien+=n;
				if (ien>v.size()) {
					std::cout << std::endl;
					break; 
					//ien=v.size();
				}
			}
			
	}
	
	return 0;
}

/* ******************************************************************************************** */
boost::program_options::variables_map parseOptions(int argc, char** argv) {
	namespace po = boost::program_options;
    po::variables_map vm;
	
	
    try {
    	int verbo=1;
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
//		    ("switch,s", po::bool_switch()->default_value(false), "switch option")
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
            ("verbosity,v", po::value<int>(&verbo)->default_value(0), 
                  "verbosity level")
//            ("include-path,I", 
//                 po::value< vector<string> >()->composing(), 
//                 "include path")
//			("ifile,i", po::value<string>(&stropt)->default_value(""), "Input file name.")
//			("odir", po::value<string>(&stropt)->default_value(""), "output file name")
//		    ("mask", po::value<string>(&stropt)->default_value(""), "specify mask as one of: \n"
//					"c,A,Z,r -- eg. c,-50,50,0.5 - will create a circular mask on sphere in deg."
//					"Azimuth from south westwards")
//		    ("show", po::value<bool>(&boolopt)->default_value(false), "shows the loaded image")
			("n", po::value<long>()->default_value(1), "Bin size. Number of samples per bin. "
					"Incomplete bins are ignored.")
			
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
        
        po::options_description visible("cpeds_row_average\n\n "
    			"This program ... "
    			"\n\n. "
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
            cerr << "cannot open config file: " << config_file << "\n";
//            return 0;
        }
        else {
            store(parse_config_file(ifs, config_file_options), vm);
            notify(vm);
        }
    

        // process parameter file 
        ifs.open(parameter_file.c_str());
        if (!ifs.is_open()) {
            cerr << "cannot open parameter file: " << parameter_file << "\n";
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

        if (vm.count("input-file")) {
			if (verbo>0) {
	            cout << "Input files are: \n" ;
	            for ( auto& s : vm["input-file"].as< vector<string> >() ) {
    	        	cout << s << " ";
        	    }
            	cout << "\n";
			}
        }
//        else {
//			if (verbo>0) {
//				cout << "No input files provided" << std::endl;
//				exit(0);
//			}
//        }

		if (verbo>0) cout << "Verbosity level is " << verbo << "\n";
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
//    console_sink->set_pattern("[cpeds_row_average] [%^%l%$] %v");

    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("cpeds_row_average.log", true);
    file_sink->set_level(spdlog::level::debug);

    spdlog::logger logger("cpeds_row_average", {console_sink, file_sink});
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
