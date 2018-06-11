/*!
  \file cpeds-angular-correlation-function.h - 
*/


#ifndef SRC_CPEDS_ANGULAR_CORRELATION_FN_H_
#define SRC_CPEDS_ANGULAR_CORRELATION_FN_H_

#include "cpeds-direction_set.h"
#include "Mscs-function.h"

mscsFunction cpeds_calculate_angular_correlation_fn(cpedsDirectionSet ds, double theta_min, double theta_max, double resolution);

mscsFunction cpeds_calculate_angular_S_correlation_statistic(cpedsDirectionSet ds, double theta_min, double theta_max, double resolution);

#endif /* SRC_CPEDS_ANGULAR_CORRELATION_FN_H_ */ 



