/*!
 * \file StringParser.cpp
 *
 *  Created on: Jun 17, 2020
 *      Author: blew
 */
#include "StringParser.h"
#include <iostream>

namespace cpeds {
	
	
	
	StringParser::StringParser() {
		// TODO Auto-generated constructor stub
		
	}
	
	StringParser::~StringParser() {
		// TODO Auto-generated destructor stub
	}
	
	bool StringParser::hasKey(std::string key) {
		std::map<std::string, std::string>::iterator it = dict.find(key);
		if(it != dict.end()) {		   return true;		}
		return false;
	}
	
	int StringParser::getValsCount(std::string key) {
		std::string val=dict[key];
		trim(val);
		return count_spaces(val)+1;
		
	}

	std::vector<std::string> StringParser::split(std::string cmd,char delim) {
	    std::stringstream ss(cmd);
	    std::string token;
	    std::vector<std::string> val;
	    while (std::getline(ss, token, delim)) {
	        val.push_back(token);
	    }
		return val;
	}
	
	std::map<std::string, std::string> StringParser::parse(std::string s,
			char fieldDelim, char kvDelim) {
		auto fields=split(s,fieldDelim);

		for (auto field : fields) {
			auto kv=split(field,kvDelim);
			if (kv.size()==2) {
				trim(kv[0]);
				trim(kv[1]);
				dict[kv[0]]=kv[1];
//				std::cout << kv[0] << " = " << kv[1] << "\n";
			}
		}
		return dict;
	}

//	std::string StringParser::remove_extension(std::string s) {
//		
//	}
	
	
} /* namespace acclass */

