/*!
 * \file test-cpedsMCMC-2Dcontours.cpp
 *
 *  Project: cpems
 *  Created on: Mar 6, 2020 2:39:28 PM
 *  Author: blew
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <tuple>
#include <boost/program_options.hpp>
#include <spdlog/spdlog.h>
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"


#include "cpedsMCMC.h"
#include "MscsPDF2D.h"

#define CONFIG_FILE "test-cpedsMCMC-2Dcontours.cfg"


spdlog::logger getLogger();
string getCmdString(int argc, char** argv);
string getProgramVersionString();
boost::program_options::variables_map  parseOptions(int argc, char** argv);



/* ******************************************************************************************** */
/* ******************************************************************************************** */
/* ******************************************************************************************** */
int main(int argc, char **argv) {
	
//	cpedsMsgs msgs("test-cpedsMCMC-2Dcontours",false,"",High);
//	msgs.setSaveRunWriteMode('a');
//	msgs.saveThisRun(argc,argv);
//	string s;
	// initialize parser
	boost::program_options::variables_map opt=parseOptions(argc,argv);

//	std::cout << opt["compareModels"].as<bool>();
	
	spdlog::logger logger=getLogger();
	logger.info("test-cpedsMCMC-2Dcontours - NEW RUN");
	logger.info(getCmdString(argc,argv));

	vector<string> input_file=opt["input-file"].as< vector<string> >();
    if (input_file.size()>0) {
    	logger.info("Loading data from: {}", input_file[0]);
    }

    MscsPDF2D pdf;
    pdf.setVerbosityLevel(High);
    pdf.loadHDF5(input_file[0],opt["dset"].as<string>());
    pdf.printInfo();
    pdf.saveSlice(2,0,input_file[0]+"-dset_"+opt["dset"].as<string>());
    double LVL;
	mscsFunction CR=pdf.getContour(opt["CL"].as<double>(),&LVL);

	CR.save(opt["ofile"].as<string>());
	cout << opt["CL"].as<double>() << "% CL contour is: " << LVL << "\n";
	
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
        string config_file;
    
        // Declare a group of options that will be 
        // allowed only on command line
        po::options_description generic("Generic options");
        generic.add_options()
            ("version,V", "print version string")
//			("ifile,i", po::value<string>(&stropt)->multitoken()->required(), "input file name")
//		    ("version,V", po::bool_switch()->default_value(false), "print program version and exit")
		    ("switch,s", po::bool_switch()->default_value(false), "switch option")
            ("help", "produce help message")
            ("config,c", po::value<string>(&config_file)->default_value(CONFIG_FILE),
                  "name of a file of a configuration.")
            ;
    
        // Declare a group of options that will be 
        // allowed both on command line and in
        // config file
        po::options_description config("Configuration");
        config.add_options()
            ("verbosity,v", po::value<int>(&opt)->default_value(10), 
                  "verbosity level")
//            ("include-path,I", 
//                 po::value< vector<string> >()->composing(), 
//                 "include path")
			("ofile,o", po::value<string>(&stropt)->default_value(""), "output file name")
			("dset", po::value<string>(&stropt)->default_value(""), "hdf5 dataset name")
			("CL", po::value<double>(&dbl)->default_value(0.68), "confidence level")
            ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden("Hidden options");
        hidden.add_options()
            ("input-file", po::value< vector<string> >(), "input .hdf5 file")
            ;

        
        po::options_description cmdline_options;
        cmdline_options.add(generic).add(config).add(hidden);

        po::options_description config_file_options;
        config_file_options.add(config).add(hidden);

        po::options_description visible("test-cpedsMCMC-2Dcontours\nInstitude of Astronomy, UMK, POLAND.\n "
    			"This program ...\n\n. "
    			"\n"
    			"example usage: test-cpedsMCMC-2Dcontours"
    			"\n"
    			"\n\n"
    			"Allowed options");
        visible.add(generic).add(config);
        
        po::positional_options_description p;
        p.add("input-file", -1);
        
        store(po::command_line_parser(argc, argv).
              options(cmdline_options).positional(p).run(), vm);
        notify(vm);
        
        ifstream ifs(config_file.c_str());
        if (!ifs) {
            cout << "can not open config file: " << config_file << "\n";
//            return 0;
        }
        else {
            store(parse_config_file(ifs, config_file_options), vm);
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
            cout << "Input files are: \n" ;
            for ( auto& s : vm["input-file"].as< vector<string> >() ) {
            	cout << s << " ";
            }
            cout << "\n";
        }

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
spdlog::logger getLogger() {
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_level(spdlog::level::info);
//    console_sink->set_pattern("[test-cpedsMCMC-2Dcontours] [%^%l%$] %v");

    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("test-cpedsMCMC-2Dcontours.log", false);
    file_sink->set_level(spdlog::level::debug);

    spdlog::logger logger("test-cpedsMCMC-2Dcontours", {console_sink, file_sink});
    logger.set_level(spdlog::level::info);

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
