/*!
 * \file cpedsTime.h
 * 
 * measures time between events in microseconds
 *
 *  Created on: Apr 27, 2012
 *      Author: blew
 */

#ifndef CPEDSTIME_H_
#define CPEDSTIME_H_

/* **************************************************************************************************** */
/* INCLUDES */
/* **************************************************************************************************** */

/* STANDALONE HEADERS */
#include <sys/time.h>
#include <sys/resource.h>
#include <QtCore/QList>
#include <QtCore/QMap>
/* INTERDEPENDENT HEADERS */

/* FORWARD DECLARATIONS */

/* USING NAMESPACES (only for inside-header files implementations) */

/* **************************************************************************************************** */
/* CLASS DECLARATION */
/* **************************************************************************************************** */
/*!
 \class 
 \brief Encapsulates 
 \details 
 
 \date Apr 27, 2012, 10:54:47 PM 
 \author Bartosz Lew
 */

class cpedsTime {
	public:
		typedef struct timeval timeStruct_t;
		cpedsTime();
		virtual ~cpedsTime();
		
		void start(int timer);
		void stop(int timer);
		double timeDiff(int timer);
		
	protected:
		QList<timeStruct_t> _timeStart; 
		QList<timeStruct_t> _timeStop; 
		QMap<int,int> _timer; 
		
		struct rusage ru;
		timeStruct_t _init;
};

#endif /* CPEDSTIME_H_ */
