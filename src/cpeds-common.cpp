/*!
 * \file cpeds-common.cpp
 *
 *  Created on: Feb 25, 2020
 *      Author: blew
 */

#include "cpeds-common.h"
#include <sstream>

std::vector<std::string> cpeds_strsplit(std::string cmd,char delim) {
	return cpeds_strToVec(cmd,delim);
}

std::vector<std::string> cpeds_strToVec(std::string cmd,char delim) {
    std::stringstream ss(cmd);
    std::string token;
    std::vector<std::string> val;
    while (std::getline(ss, token, delim)) {
        val.push_back(token);
    }
	return val;
}
