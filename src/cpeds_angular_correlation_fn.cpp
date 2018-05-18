/*!
 * \file cpeds_angular_correlation_fn.cpp
 *
 *  Created on: Apr 19, 2018
 *      Author: blew
 */

#include "cpeds_angular_correlation_fn.h"

mscsFunction cpeds_calculate_angular_correlation_fn(cpedsDirectionSet ds, double theta_min, double theta_max, double resolution) {
	mscsFunction Cth;
	double ang;
	long i,j,corr_i;

	if (theta_min<resolution) printf("INFO: The requested minimal angle in the correlation function"
			"is smaller than the bin resolution.\n");
	
	long point_num_C_th = (int) (ceil((theta_max - theta_min) / resolution));
	
	cpedsList<double> separation_number;
	
	Cth.setPointsNum(point_num_C_th);
	separation_number.makeLength(point_num_C_th);

	
	long pix_num = ds.size();
	for (i = 0; i < pix_num; i++) {
		for (j = i; j < pix_num; j++) {
			
			ang = ds[i].angle(ds[j]); 
			
			if ((ang >= theta_min) && (ang <= theta_max)) {
				corr_i = (long) round((ang - theta_min) / resolution);
				separation_number[corr_i]++; // this stores the number of given separations on a sky ( for normalization purposes)
				Cth[corr_i].rx() += ang;
				Cth[corr_i].ry() += ds[i].val() * ds[j].val(); 
			}
		}
		
		printf("calculating: %li of %li\r", i, pix_num);
	}
	
	// normalization and averaging
	
	for (i = 0; i < point_num_C_th; i++) {
		Cth[i].rx() /= separation_number[i]; // th  -- this gives the average angle over all angles that fall into this range (bin) limited by the resolution parameter
		Cth[i].ry() /= separation_number[i]; // normalization of C(th)
	}
	
	
	return Cth;
}
