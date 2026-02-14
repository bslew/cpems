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
#include "cpeds-math.h"
#include "cpeds-list.h"
//#include <cpems/Mscs-function3dregc.h>
//#include <Eigen/Geometry>
//#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/SVD>
//#include <Eigen/Eigen>
//#include <eigen3/Eigen/Eigen>
//#include <yaml-cpp/yaml.h>


#define PARAMETER_FILE "input.par"

using namespace std;
using namespace Eigen;

spdlog::logger getLogger(int verbosity);
string getCmdString(int argc, char** argv);
string getProgramVersionString();
boost::program_options::variables_map  parseOptions(int argc, char** argv);



//int test_direction(boost::program_options::variables_map &opt,spdlog::logger& logger) {
//	logger.info("dt [UTC]: {}",opt["dtSt"].as<string>());
//	return 0;
//	
//}


int save_matrix(Eigen::MatrixXd M, std::string fname) {
	int rows,cols;
	std::ofstream ofs;
	ofs.open(fname,std::ofstream::out);
		
	for (long i = 0; i < M.rows(); i++) {
		for (long j = 0; j < M.cols(); j++) {
			ofs << M(i,j) << " ";
		}
		ofs << "\n";
	}
//	ofs << "\n";
	ofs.close();
	return 0;
}

Eigen::MatrixXd load_matrix(std::string fname, 
		boost::program_options::variables_map& opt,
		spdlog::logger* logger=0) {
	int rows,cols,max_rows=-1,max_cols=-1;
//	char ch='\n';

	if (opt["Nrows"].as<long>()>0) {
		rows=opt["Nrows"].as<long>();
		if (logger!=0) logger->info("Ignoring counting file rows");
	}
	else {
		max_rows=cpeds_get_txt_file_non_empty_lines_count(fname);
		rows=max_rows;
	}
	
	if (opt["Ncols"].as<long>()>0) {
		cols=opt["Ncols"].as<long>();
		if (logger!=0) logger->info("Ignoring counting file cols");
	}
	else {
		max_cols=cpeds_get_file_cols_num_first_ln(fname);
		cols=max_cols;
	}
	if (logger!=0) logger->info("File {} has (rows,cols)=({},{})\n",fname, max_rows,max_cols);

	if (rows>max_rows and max_rows!=-1) rows=max_rows;
	if (cols>max_cols and max_cols!=-1) cols=max_cols;
	
	if (logger!=0) logger->info("Will load (rows,cols)=({},{})",rows,cols);
	
	Eigen::MatrixXd M(1,1);
	M.resize(rows, cols);
	
	long promptEvery=long(pow(10,long(log10(rows)-1)));
	double val;
	long i=0,r=0,c=0;
	std::ifstream ifs(fname);
	while (ifs >> val) {
		M(r,c)=val;
		if (opt["verbosity"].as<int>()>2)
			cout << "read " << val <<" into (r,c)=(" << r << "," << c << ")\n";

		i++;

		if (c+1==cols) ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

		c=i%cols;
		r=i/cols;
		
		if (r==rows) break;
		if (opt["verbosity"].as<int>()>1) {
			if (logger!=0 and c==0 and (r % promptEvery==0)) logger->info("Loaded {}/{} rows",r,rows);
		}
	}
	ifs.close();

	if (logger!=0) logger->info("Loaded matrix of size (rows,cols)=({},{})",M.rows(),M.cols());

	return M;
}

/* ******************************************************************************************** */
/* ******************************************************************************************** */
/* ******************************************************************************************** */
int main(int argc, char **argv) {
	
	// initialize parser
	boost::program_options::variables_map opt=parseOptions(argc,argv);
	
	spdlog::logger logger=getLogger(opt["verbosity"].as<int>());
	if (opt["verbosity"].as<int>()>0) {
		logger.info("cpeds_calculate_cov - NEW RUN");
		logger.info(getCmdString(argc,argv));
	}
	
	std::string infile;
	std::string outfile=opt["ofile"].as<string>();
//	std::cout << opt["input_file"].as< vector<string> >() << std::endl;
	if (opt.count("input-file")) {
		vector<string> input_files=opt["input-file"].as< vector<string> >();
		if (input_files.size()>=1) infile=input_files[0];
//		if (input_files.size()>=2) outfile=input_files[1]; 
	}
	
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
	Eigen::initParallel();
	
	if (opt["sliceCrossDiag"].as<bool>()) {
		Eigen::MatrixXd M;
		M=load_matrix(infile,opt);
		long cd=opt["i"].as<long>();
		long j=0;
		for (long i=M.rows()-cd-1;i>=0;i--) {
			if (j>0 and i>0) std::cout << j << " " << M(i,j) << "\n";
			j++;
		}
		return 0;
	}

	if (opt["diagStats"].as<bool>()) {
		Eigen::MatrixXd M;
		M=load_matrix(infile,opt);
		std::ofstream ofs(outfile);
		
		for (long i=0;i<M.rows();i++) {
			cpedsList<double> d;
			long j=0;
			long ii=i;
			while (j<M.cols() and ii<M.rows()) {
				d.append(M(ii++,j++));
			}
			ofs << i << " " << d.median() << " " << d.std() << "\n";
		}
		ofs.close();
		return 0;
		
	}
	
	
	if (opt["sliceDiag"].as<bool>()) {
		Eigen::MatrixXd M;
		M=load_matrix(infile,opt);
		long cd=opt["i"].as<long>();
		long j=0;
		for (long i=cd;i<M.rows();i++) {
			if (j>0 and i>0) std::cout << j << " " << M(i,j) << "\n";
			j++;
		}
		return 0;
	}
	
	if (opt["cov"].as<bool>()) {
		Eigen::MatrixXd dcov,C;
		dcov=load_matrix(infile,opt,&logger);
		long N=dcov.size();
		double* D = new double[N];
		for (long i = 0; i < dcov.rows(); i++) {
			for (long j = 0; j < dcov.cols(); j++) {
				D[dcov.cols()*i+j]=dcov(i,j);
			}
		}
//		double* Cov = cpeds_calculate_covariance_matrix_para(D,dcov.cols(),dcov.rows(),false,
//				opt["max_diag"].as<long>());
		logger.info("Calculating {} diagonals",opt["max_diag"].as<long>());
		double* Cov = cpeds_calculate_covariance_matrix_para(D,dcov.cols(),dcov.rows(),false,
				opt["max_diag"].as<long>());

		delete [] D;
		
		C.resize(dcov.cols(),dcov.cols());
		for (long i = 0; i < dcov.cols(); i++) {
			for (long j = 0; j < dcov.cols(); j++) {
				C(i,j)=Cov[dcov.cols()*i+j];
			}
		}
		delete [] Cov;
		
		logger.info("saving covariance matrix to: {}",outfile);
		save_matrix(C, outfile);
		return 0;
	}
	
	if (opt["invert"].as<bool>()) {
		Eigen::MatrixXd C,Cinv;

		C=load_matrix(infile,opt,&logger);
		Cinv=C.inverse();
		save_matrix(Cinv, outfile);
		logger.info("saving inverted matrix to: {}",outfile);
		return 0;
	}	
	
	if (opt["cond"].as<bool>()) {
		Eigen::MatrixXd C=load_matrix(infile,opt,&logger);
		
		BDCSVD<MatrixXd> svd(C);
		double cond = svd.singularValues()(0) 
		    / svd.singularValues()(svd.singularValues().size()-1);		
		logger.info("matrix condition number: {}",cond);
		return 0;
	}

	if (opt["det"].as<bool>()) {
		Eigen::MatrixXd C=load_matrix(infile,opt,&logger);
		
//		BDCSVD<MatrixXd> svd(C);
		double det = C.determinant();
		logger.info("det: {}",det);
		logger.info("ln det: {}",log(det));
		return 0;
	}

	// prepare data for saving
	
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
            		("verbosity,v", po::value<int>(&opt)->default_value(1), 
            				"verbosity level")
					//            ("include-path,I", 
					//                 po::value< vector<string> >()->composing(), 
					//                 "include path")
					("test", po::value<bool>()->default_value(false), "generate matrix and invert")
					("sliceCrossDiag", po::value<bool>()->default_value(false), "slice matrix along cross"
							"diagonal. Use option -i to specify distance from the main diagonal.")
					("sliceDiag", po::value<bool>()->default_value(false), "slice matrix along a"
							"diagonal. Use option -i to specify distance from the main diagonal.")
					("diagStats", po::value<bool>()->default_value(false), "calculate diagonal "
							"statistics (median value and std) and save to output file")
					("cond", po::value<bool>()->default_value(false), "condition number of input matrix")
					("det", po::value<bool>()->default_value(false), "determinant of input matrix")
					("invert", po::value<bool>()->default_value(false), "invert input matrix")
					("cov", po::value<bool>()->default_value(false), "calculate covariance matrix."
							"Rows in input file should contain complete observations of all variables.")
					("ofile,o", po::value<string>()->default_value("cov.txt"), "output file name")
					("max_diag", po::value<long>()->default_value(-1), "--cov suboption. "
							"Number of sub-diagonals to calculate."
							"-1 - calculates the full covariance matrix. "
							"0 - calculates only the diagonal, "
							"1 - calculates up to 1st sub-diagonal (inclusive)")
					("Nrows", po::value<long>()->default_value(-1), ""
							"Load only Nrows from input file. -1 loads the entire file (default)")
					("Ncols", po::value<long>()->default_value(-1), ""
							"Load only Ncols from input file. -1 loads all columns (default)")
					("i,id", po::value<long>()->default_value(0), ""
							"slicing parameter (default)")
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
		
		if (opt>0) cout << "Verbosity level is " << opt << "\n";                
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
