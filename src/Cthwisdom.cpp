/*!
 * \file Cthwisdom.cpp
 *
 *  Created on: Apr 29, 2019
 *      Author: blew
 */
#include <Cthwisdom.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "cpeds-list.h"

#define DEBUG_LOAD 1

#ifdef DEBUG_LOAD
typedef std::chrono::high_resolution_clock Clock;
#endif


Cthwisdom::Cthwisdom(double theta_min, double theta_max, double step, double bin) {
//	std::vector<long> j;
	Cthwisdom::element_t el;
	std::vector<element_t> dummy;
//	dummy.resize(Npix,el);

	long size=get_size(theta_min,theta_max,step);
	wisdom.idx.resize(size,el);
	wisdom.ang.resize(size);
	wisdom.hits.resize(size);
	
	wisdom.res=step;
	wisdom.bin=bin;
	wisdom.halfbin=bin/2;

	wisdom.min=theta_min; //-wisdom.halfbin;
//	wisdom.max=theta_max;
	wisdom.max=theta_min+(size-1)*step; //+wisdom.halfbin;
	wisdom.Npix=0;

	for (long i = 0; i < size; i++) {
//		wisdom.idx[i].resize(Npix-1,el);
		wisdom.ang[i]=theta_min+i*step;
		wisdom.hits[i]=0;
	}
}

Cthwisdom::~Cthwisdom() {
}
/* ******************************************************************************************** */
long Cthwisdom::get_size(double theta_min, double theta_max, double step) {
	return static_cast<long>((theta_max-theta_min)/step)+1;
}

/* ******************************************************************************************** */
void Cthwisdom::saveHDF(std::string fname) {
	if (fname=="") {	fname=get_wisdom_file_name();	}

}
/* ******************************************************************************************** */
void Cthwisdom::save(std::string fname) {
	std::ofstream ofile;
	
	if (fname=="") {	fname=get_wisdom_file_name();	}
	ofile.open(fname);
	
	// in the 1st row we save Cth parameters and Healpix map parameters
	ofile << wisdom.min << " " 
			<< wisdom.max << " " 
			<< wisdom.res << " " 
			<< wisdom.bin << " "  
			<< wisdom.Npix << " "  
			<< "\n";
	
	// in the 2nd row we save the angles that the Cth will have
	for (auto& ang : wisdom.ang) {
		ofile << ang << " ";
	}
	ofile << "\n";

	// in the 3rf row we save the number of pairs in bin
	for (auto& hits : wisdom.hits) {
		ofile << hits << " ";
	}
	ofile << "\n";

	// in each next row save all pairs of indexes that contribute to a single angle
	// (one angle per row)
	for (auto& idx : wisdom.idx) {
		pairs_t::iterator iit=idx.pairs.begin(); // i iterator
		while (iit!=idx.pairs.end()) {
			ofile << iit->first << " ";
			
			std::vector<long>::iterator jit=iit->second.begin();
//			for (auto& j : iit->second) {
			while (jit!=iit->second.end()) {
				ofile << (*jit);
				jit++;
				if (jit!=iit->second.end()) ofile << " ";
			}
			
			iit++;
			if (iit!=idx.pairs.end()) ofile << ",";
		}
		ofile << std::endl;
	}
	ofile.close();	
}
/* ******************************************************************************************** */
int Cthwisdom::load(std::string fname) {
#ifdef DEBUG_LOAD
	auto t0=Clock::now();
#endif
	
	
	std::ifstream ifile(fname);
	std::string line;

	
	if (fname=="") {	fname=get_wisdom_file_name();	}
	
	if (not ifile.is_open()) {
		throw "Cannot open file: "+fname+"\n";
	}
	// read 1st row
//	long dummyl;
	getline(ifile,line);
	std::istringstream fst(line);
//	fst >> dummyl;
	fst >> wisdom.min;
	fst >> wisdom.max;
	fst >> wisdom.res;
	fst >> wisdom.bin;
	fst >> wisdom.Npix;

	
	// read 2nd row
	getline(ifile,line);
	std::istringstream full_line(line);
	double ang;
	std::string ang_str;
	while (getline(full_line,ang_str,' ')) {
		std::stringstream iss(ang_str);
		iss >> ang;
		wisdom.ang.push_back(ang);
	}

	
	// read 3rd row
	getline(ifile,line);
	full_line.str(line);
	double hits;
	std::string hits_str;
	while (getline(full_line,hits_str,' ')) {
		std::stringstream iss(hits_str);
		iss >> hits;
		wisdom.hits.push_back(hits);
	}

	
	// remaining rows
#ifdef DEBUG_LOAD
	auto t1=Clock::now();
	cpedsList<double> load_time_rows;
	cpedsList<double> process_time_rows;
	cpedsList<double> mkwis_time_rows;
#endif
	
	while (getline(ifile,line)) {
		std::istringstream line_stream(line);

#ifdef DEBUG_LOAD
		auto t2=Clock::now();
		double load_time=double(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count())/1.e6;
//		std::cout << "Load line time [s]: " << load_time << std::endl;
		load_time_rows.append(load_time);
#endif
		

#ifdef DEBUG_LOAD
		auto proc_row_t1=Clock::now();
		double wis_time=0;
#endif
		std::string pairs;
		element_t el;
		while (getline(line_stream,pairs,',')) {
			std::istringstream pair_ss(pairs);
//			std::string jidx;
			long i=0;
			pair_ss >> i;
			long j;
			
#ifdef DEBUG_LOAD
			auto mkwis_t1=Clock::now();
#endif

			while (pair_ss >> j) {
				el.pairs[i].push_back(j);
			}
#ifdef DEBUG_LOAD
			auto mkwis_t2=Clock::now();
			wis_time+=double(std::chrono::duration_cast<std::chrono::microseconds>(mkwis_t2 - mkwis_t1).count())/1.e6;
#endif
			
			
		}
		wisdom.idx.push_back(el);
#ifdef DEBUG_LOAD
		mkwis_time_rows.append(wis_time);
//		std::cout << "make wisdom time [s]: " << wis_time << std::endl;

		auto proc_row_t2=Clock::now();
		double process_time=double(std::chrono::duration_cast<std::chrono::microseconds>(proc_row_t2 - proc_row_t1).count())/1.e6;
//		std::cout << "Process line time [s]: " << process_time << std::endl;
		process_time_rows.append(process_time);
#endif
		
#ifdef DEBUG_LOAD
		t1=Clock::now();
//		std::cout << "\n\n";
#endif
	}
	ifile.close();	
	

#ifdef DEBUG_LOAD
	auto t_fin=Clock::now();
	std::cout << "Load time total [s]: " 
			<< std::chrono::duration_cast<std::chrono::microseconds>(t_fin-t0).count()/1.e6
			<< std::endl;
	std::cout << "Load time split:\n";
	std::cout << "	- IO time [s]: " 
			<< load_time_rows.sum()
			<< std::endl;
	std::cout << "	- rows process time [s]: " 
			<< process_time_rows.sum()
			<< std::endl;
	std::cout << "		- making wisdom structure time [s]: " 
			<< mkwis_time_rows.sum()
			<< std::endl;
	
#endif
	return 0;
	
}

std::string Cthwisdom::get_wisdom_file_name() {
	std::stringstream wisss;
//	wisss << "wisdom-" << wisdom.Npix << "_" << std::setprecision(2) << wisdom.min << "_" << wisdom.max << "_" << wisdom.res << ".txt";
	wisss << "wisdom-" 
			<< std::setprecision(2) << std::fixed 
			<< wisdom.min << "_" 
			<< wisdom.max << "_" 
			<< wisdom.res << "_" 
			<< wisdom.bin << ".txt";
	return wisss.str();

}

void Cthwisdom::compile() {
	// calculate mean angles in bins
	for (long i=0; i<wisdom.ang.size(); i++) {
		if (wisdom.hits[i]>0) wisdom.ang[i]/=wisdom.hits[i];
	}
}

void Cthwisdom::add_wisdom(double angle, long i, long j) {
	long idx=get_ang_idx(angle);
	if (idx>=0) {
		wisdom.ang[idx]+=angle;
		wisdom.hits[idx]++;
		wisdom.idx[idx].pairs[i].push_back(j);
	}
}

long Cthwisdom::get_ang_idx(double angle) {
	if (angle<wisdom.min) return -1;
	if (angle>wisdom.max) return -1;
//	long n=static_cast<long>((angle-wisdom.min+wisdom.halfbin)/wisdom.res);
	long n=(angle-wisdom.min+wisdom.halfbin)/wisdom.res;
	double offset=angle-(wisdom.min+n*wisdom.res);
	if (fabs(offset)>wisdom.halfbin) return -1;
//	if (offset>0)
	return n;
}
