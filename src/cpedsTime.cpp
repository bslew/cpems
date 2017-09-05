/*!
 * \file cpedsTime.cpp
 *
 *  Created on: Apr 27, 2012
 *      Author: blew
 */


#include "cpedsTime.h"
#include <stdio.h>


/***************************************************************************************/
cpedsTime::cpedsTime() {
	getrusage(RUSAGE_SELF, &ru);
	_init=ru.ru_utime;
}
/***************************************************************************************/
cpedsTime::~cpedsTime() {
}
/***************************************************************************************/
void cpedsTime::start(int timer) {
		_timer[timer]=_timeStart.count();
		getrusage(RUSAGE_SELF, &ru);
		_timeStart.append(ru.ru_utime);
		_timeStop.append(_init);
}


/***************************************************************************************/
void cpedsTime::stop(int timer) {
	getrusage(RUSAGE_SELF, &ru);
	_timeStop[_timer[timer]]=ru.ru_utime;
}


/***************************************************************************************/
double cpedsTime::timeDiff(int timer) {
	double tS = double(_timeStart[_timer[timer]].tv_sec)*1000000 + double(_timeStart[_timer[timer]].tv_usec);
	double tE = double(_timeStop[ _timer[timer]].tv_sec)*1000000 + double(_timeStop[ _timer[timer]].tv_usec);
	printf("tS: %lf\n",tS);
	printf("tE: %lf\n",tE);
	printf("tE: %lf\n",double(tE-tS)*1e-6);
	return double(tE-tS)*1e-6;
}


