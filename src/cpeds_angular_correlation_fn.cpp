/*!
 * \file cpeds_angular_correlation_fn.cpp
 *
 *  Created on: Apr 19, 2018
 *      Author: blew
 */

#include "cpeds_angular_correlation_fn.h"
#include <omp.h>

mscsFunction cpeds_calculate_angular_correlation_fn(cpedsDirectionSet ds, double theta_min, double theta_max, double resolution) {

	if (theta_min<resolution) printf("INFO: The requested minimal angle in the correlation function"
			"is smaller than the bin resolution.\n");
	
	long point_num_C_th = (int) (ceil((theta_max - theta_min) / resolution));
	
	cpedsList<double> separation_number;
//	cpedsList<double> separation_vals;
	
	mscsFunction Cth;
	Cth.setPointsNum(point_num_C_th);
	separation_number.makeLength(point_num_C_th);
//	separation_vals.makeLength(point_num_C_th);

	long i,j,Ndone=0;
	long pix_num = ds.size();
#pragma omp parallel 
	{
		double ang=0;
		long corr_i=0;
		mscsFunction Cth_priv=Cth;
		cpedsList<double> separation_number_priv=separation_number;
		
#pragma omp for schedule(dynamic)
		for (long i = 0; i < pix_num-1; i++) {
			cpedsDirection d1=ds[i];
			for (long j = i+1; j < pix_num; j++) {
				cpedsDirection d2=ds[j];
				ang=d1.angle(d2);
				
				if ((ang >= theta_min) && (ang <= theta_max)) {
					corr_i = (long) round((ang - theta_min) / resolution);
					separation_number_priv[corr_i]++; // this stores the number of given separations on a sky ( for normalization purposes)
					//				separation_vals[corr_i].append()++; // this stores the number of given separations on a sky ( for normalization purposes)
					Cth_priv[corr_i].rx() += ang;
					Cth_priv[corr_i].ry() += d1.val() * d2.val(); 
//					CthX[corr_i] += ang;
//					CthY[corr_i] += ds[i].val() * ds[j].val(); 
				}
			}
#pragma omp critical
			{
			Ndone=Ndone+1;
			if (Ndone%10==0) {
				printf("calculating: %li of %li\r", Ndone, pix_num);
			}
			}
		}
#pragma omp single
		printf("\ndone\ncombining results\n");
		
#pragma omp critical
		{
			for (long k=0; k<point_num_C_th; k++) {
				Cth[k].rx()+=Cth_priv[k].rx();
				Cth[k].ry()+=Cth_priv[k].ry();
				separation_number[k]+=separation_number_priv[k];
			}
		}
	}
/*
	for (i = 0; i < pix_num; i++) {
		for (j = i; j < pix_num; j++) {
			
			ang = ds[i].angle(ds[j]); 
			
			if ((ang >= theta_min) && (ang <= theta_max)) {
				corr_i = (long) round((ang - theta_min) / resolution);
				separation_number[corr_i]++; // this stores the number of given separations on a sky ( for normalization purposes)
//				separation_vals[corr_i].append()++; // this stores the number of given separations on a sky ( for normalization purposes)
				Cth[corr_i].rx() += ang;
				Cth[corr_i].ry() += ds[i].val() * ds[j].val(); 
			}
		}
		
		printf("calculating: %li of %li\r", i, pix_num);
	}
*/

	// normalization and averaging
	printf("normalization\n");
	
	for (long i = 0; i < point_num_C_th; i++) {
//		if (separation_number[i]>0) {
			Cth[i].rx() *= PI180inv/separation_number[i]; // th  -- this gives the average angle over all angles that fall into this range (bin) limited by the resolution parameter
			Cth[i].ry() /= separation_number[i]; // normalization of C(th)
//		}
	}
	
	printf("cleaning up\n");
	Cth.removeNans();
	printf("Cth done\n");
	
	return Cth;
}


/* ******************************************************************************************** */
mscsFunction cpeds_calculate_angular_S_correlation_statistic(cpedsDirectionSet ds, double theta_min, double theta_max, double resolution) {
/*
	if (theta_min<resolution) printf("INFO: The requested minimal angle in the correlation function"
			"is smaller than the bin resolution.\n");
	
	long point_num_C_th = (int) (ceil((theta_max - theta_min) / resolution));
	
	cpedsList<double> separation_number;
	
	mscsFunction Cth;
	Cth.setPointsNum(point_num_C_th);
	separation_number.makeLength(point_num_C_th);

	long i,j,Ndone=0;
	long pix_num = ds.size();
#pragma omp parallel 
	{
		double ang=0;
		long corr_i=0;
		mscsFunction cross_term=Cth;
		mscsFunction auto_term1=Cth;
		mscsFunction auto_term2=Cth;
		cpedsList<double> separation_number_priv=separation_number;
		
#pragma omp for schedule(dynamic)
		for (long i = 0; i < pix_num-1; i++) {
			cpedsDirection d1=ds[i];
			auto_term1[corr_i].rx()++;
			auto_term1[corr_i].ry()+=d1.val()*d1.val();
			for (long j = i+1; j < pix_num; j++) {
				cpedsDirection d2=ds[j];
				auto_term2[corr_j].rx()++;
				auto_term2[corr_j].ry()+=d2.val()*d2.val();

				
				ang=d1.angle(d2);
				
				if ((ang >= theta_min) && (ang <= theta_max)) {
					corr_i = (long) round((ang - theta_min) / resolution);
					separation_number_priv[corr_i]++; // this stores the number of given separations on a sky ( for normalization purposes)
					cross_term[corr_i].rx() += ang;
					cross_term[corr_i].ry() += d1.val() * d2.val(); 
				}
			}
#pragma omp critical
			{
			Ndone=Ndone+1;
			if (Ndone%10==0) {
				printf("calculating: %li of %li\r", Ndone, pix_num);
			}
			}
		}
#pragma omp single
		printf("\ndone\ncombining results\n");
		
#pragma omp critical
		{
			for (long k=0; k<point_num_C_th; k++) {
				Cth[k].rx()+=cross_term[k].rx();
				Cth[k].ry()+=cross_term[k].ry();
				separation_number[k]+=separation_number_priv[k];
			}
		}
	}

	// normalization and averaging
	printf("normalization\n");
	
	for (long i = 0; i < point_num_C_th; i++) {
//		if (separation_number[i]>0) {
			Cth[i].rx() *= PI180inv/separation_number[i]; // th  -- this gives the average angle over all angles that fall into this range (bin) limited by the resolution parameter
			Cth[i].ry() /= separation_number[i]; // normalization of C(th)
//		}
	}
	
	printf("cleaning up\n");
	Cth.removeNans();
	printf("Cth done\n");
	
	return Cth;
	
*/
}
