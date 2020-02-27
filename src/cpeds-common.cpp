/*!
 * \file cpeds-common.cpp
 *
 *  Created on: Feb 25, 2020
 *      Author: blew
 */

#include "cpeds-common.h"
#include <sstream>

std::vector<std::string> cpeds_strToVec(std::string cmd,char delim) {
    std::stringstream ss(cmd);
    std::string token;
    std::vector<std::string> val;
    while (std::getline(ss, token, delim)) {
        val.push_back(token);
    }
	return val;
}
