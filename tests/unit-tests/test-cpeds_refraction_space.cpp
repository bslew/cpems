/*!
 * \file test-cpeds_refraction_space.cpp
 *
 *  Created on: Apr 29, 2020
 *      Author: blew
 */

#include <stdlib.h>
#include <string>
#include "cpeds-math.h"
#include <vector>


int main(int argc, char** argv) {
//	std::cout << argc;
	if (argc==1) { 
		std::cout << "USAGE: test-cpeds_refraction_sapce ZDobs" << std::endl
		<< "where ZDobs - observed zenith distance [deg]" << std::endl; 
	return 0; 
	}
	
	double ZDobs=std::stod(argv[1]);
	std::vector<double> Alt={0,300,1000};
	std::vector<double> Temp={-20,0,20};
	std::vector<double> Pres={980,1020};
	std::vector<double> RelHum={50,90,100};
	double lambda=0.0001;
	double lat=54;
	double lapse=0.0065;
	double acc=1e-8;
	
	double err, maxErr=0;
	
	for (auto alt : Alt) {
		for (auto T : Temp ) {
			for (auto P : Pres) {
				for (auto RH : RelHum) {

					double ZDspace=cpeds_refraction(ZDobs,alt,T,P,RH,lambda,lat,lapse,acc);
					double ZDobs_check=cpeds_refraction_space(ZDspace,alt,T,P,RH,lambda,lat,lapse,acc);
					err=ZDobs-ZDobs_check;
					if (err>maxErr) maxErr=err;
				}
			}
		}
	}
	
//	std::cout << "Input ZDobs [deg]: " << ZDobs << std::endl;
//	std::cout << "refraction at ZDobs [deg]: " << ref << std::endl;
//	std::cout << "ZDobs check [deg]: " << ZDobs_check << std::endl;
	std::cout << "Max. Error for ZDobs "<< ZDobs << " deg is [deg]: " << maxErr << std::endl;
	return 0;
}
