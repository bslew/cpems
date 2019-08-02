/*!
 * \file Mscs-map-correlation_function.cpp
 *
 *  Created on: May 14, 2019
 *      Author: blew
 */



#include "Mscs-map.h"

mscsCorrelationFunction mscsMap::calculate_angular_correlation_fn_wis(Cthwisdom& wisdom) {
	
	long point_num_C_th = wisdom.size();
	
	cpedsList<double> separation_number;
	//	cpedsList<double> separation_vals;
	
	mscsFunction Cth;
	Cth.setPointsNum(point_num_C_th);
	separation_number.makeLength(point_num_C_th);
	
	long i,j,Ndone=0;
	long pix_num = pixNum();
#pragma omp parallel 
	{
		double ang=0;
		long corr_i=0;
		
#pragma omp for schedule(dynamic)
		for (long ang_i=0;ang_i<wisdom.size();ang_i++) {
			
			Cthwisdom::pairs_t::iterator pair_it=wisdom[ang_i].pairs.begin();
			long separation_number=0;
			while (pair_it != wisdom[ang_i].pairs.end()) {
				if (not isMasked(pair_it->first)) {
					
					std::vector<long>::iterator j_it=pair_it->second.begin();
					double contrib=0;
					
					long step=1./wisdom.use_wisdom_fraction;
					while (j_it < pair_it->second.end()) {
						if (not isMasked((*j_it))) {
							separation_number++; // this stores the number of given separations on a sky ( for normalization purposes)
//							Cth[ang_i].rx() += ang;
							contrib+=T((*j_it));
						}
						
//						j_it++;
						std::advance(j_it,step);
					}
					Cth[ang_i].ry() += T(pair_it->first) *  contrib;
				}
#pragma omp critical
				{
					Ndone=Ndone+1;
					if (Ndone%10==0) {
						printf("calculating: %li of %li\r", Ndone, pix_num);
					}
				}
				
				pair_it++;
			}
			
			Cth[ang_i].rx()=wisdom(ang_i);
			Cth[ang_i].ry()/=separation_number;
		}
	}
	
	
	printf("cleaning up\n");
	Cth.removeNans();
	printf("Cth done\n");
	
	
	
	return Cth;
}
