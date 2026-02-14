/*
 * test-cpeds-templates.cpp
 *
 *  Created on: Jan 8, 2021
 *      Author: blew
 */




#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Hello
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include "cpeds-templates.h"

int export_array() {
	
	cpeds_queue<double> q;
	q << 1 << 2 << 3 << 3.1415926535;
	auto t=q.export_array();
	std::cout << std::setprecision(15) << t[3] << "\n" << t[3]- 3.1415926535 << "\n";
	if (t[3]==3.1415926535) return 0;
	if (fabs(t[3]- 3.1415926535)<1e-10) return 0;
	return 1;
}

int iter() {
	cpeds_queue<int> q;
	q << 1 << 2 << 3 << 4;

	for (cpeds_queue<int>::iterator it=q.begin();it!=q.end();++it) {
		std::cout << "iterating: " << *it << "\n";
	}
	return 0;
}

int get_size() {
	cpeds_queue<int> q;
	q << 1 << 2 << 3;
	return q.get_size();
}

BOOST_AUTO_TEST_CASE( testsize ) {   BOOST_CHECK( get_size() == 3 ); }
BOOST_AUTO_TEST_CASE( testiterator ) {   BOOST_CHECK( iter() == 0 ); }
BOOST_AUTO_TEST_CASE( testexport_array ) {   BOOST_CHECK( export_array() ==0); }



