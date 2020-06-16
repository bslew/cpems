/*!
 * \file cpeds_get_healpix_nested_pixels.cpp
 *
 *  Created on: Jun 5, 2020
 *      Author: blew
 */

#include <tuple>
#include <iostream>
#include <string>
#include "cpeds-math.h"

int main(int argc, char** argv) {
	long ns=8;
	long cns=2;
	long pixID=2;
	
/*
	if (argc<4) {
		std::cout << "USAGE: cpeds_get_healpix_nested_pixels [ns cns pixID]" << std::endl;
		std::cout << std::endl;
		std::cout << "eg. " << "./cpeds_get_healpix_nested_pixels 16 8 0" << std::endl;
		std::cout << "assuming 8 2 1" << std::endl;
		ns=8;
		cns=2;
		pixID=1;
	}
	else {
		ns=std::stol(argv[1]);
		cns=std::stol(argv[2]);
		pixID=std::stol(argv[3]);
	}
*/
	
	
	auto p = cpeds_get_healpix_nested_pixels(ns,cns,pixID);

	std::cout << "ns: " << ns << std::endl;
	std::cout << "coarse ns: " << cns << std::endl;
	std::cout << "coarse pixID: " << pixID << std::endl;
	std::cout << std::endl;
	std::cout << "st: " << p.first << std::endl;
	std::cout << "count: " << p.second << std::endl;
	
	if (p.first!=32 or p.second!=16) return 1;
	
	return 0;
}
