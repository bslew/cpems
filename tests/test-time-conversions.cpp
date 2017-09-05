/*!
 * \file test-time-conversions.cpp
 *
 *  Created on: Mar 5, 2013
 *      Author: blew
 */

#include "cpeds-math.h"
#include <sys/time.h>

int main() {
	double timestamp=1362494324.352013;
	printf("timestamp: %.10lf\n",timestamp);
	double JD=cpeds_timeSec_to_julian_time(timestamp);
	
	struct timeval  tv;
	struct timezone tz;
	struct tm      *t;
	for (long i = 0; i < 100; i++) {
		gettimeofday(&tv, &tz);

		printf("stamp: %.15lf JD: %.15lf\n",double(tv.tv_sec)+double(tv.tv_usec)*1e-6, cpeds_julian_time());
	}
//	printf("")
}
