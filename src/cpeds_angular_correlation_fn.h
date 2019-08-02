/*!
  \file cpeds-angular-correlation-function.h - 
*/


#ifndef SRC_CPEDS_ANGULAR_CORRELATION_FN_H_
#define SRC_CPEDS_ANGULAR_CORRELATION_FN_H_

#include "cpeds-direction_set.h"
#include "Mscs-function.h"
#include "Cthwisdom.h"

mscsFunction cpeds_calculate_angular_correlation_fn(cpedsDirectionSet ds, double theta_min, double theta_max, double resolution);
/*!
	\brief calculates 2-pt angular correclation function and wisdom (if given)
	\details 
	@param ds - direction set object that holds directions and values used to calculate the correlation function
	@param hp_idx - vector containing healpix pixel indexes (nested ordering). This is used by wisdom, since the same 
	wisdom is indented to be used with different maps and different masks, hence the only constant between different ds
	(cpedsDirectionSets) will be the angles between pairs of pixels with given index numbers.
	@param theta_min - defines angle calculation range
	@param theta_max - defines angle calculation range
	@param resolution - defines angle calculation step
	@param wisdom - pointer to wisdom object. This parameter is requred. If you don't want to calculate wisdom then use
		cpeds_calculate_angular_correlation_fn(cpedsDirectionSet ds, double theta_min, double theta_max, double resolution)
	@return

	\date May 7, 2019, 5:10:19 PM
*/
void cpeds_calculate_angular_correlation_fn_wisdom(cpedsDirectionSet ds, std::vector<long> hp_idx, double theta_min, double theta_max, double resolution, Cthwisdom* wisdom);
mscsFunction cpeds_calculate_angular_correlation_fn(cpedsDirectionSet ds, Cthwisdom* wisdom);

mscsFunction cpeds_calculate_angular_S_correlation_statistic(cpedsDirectionSet ds, double theta_min, double theta_max, double resolution);

#endif /* SRC_CPEDS_ANGULAR_CORRELATION_FN_H_ */ 



