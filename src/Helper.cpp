//
//	George O'Neill @ University of York, 2020/03/12
//
//	This file holds helper functions used in various other parts of the code
//
#include "../include/Helper.h"

namespace fnt {	//	create unique working area

void Helper::countdown() {	//	print progress

	if( --countdownN == 0) {	//	check progress

		percent++;	//	increment total bars
		std::cout << "\r" + std::string(percent, 'X') + std::string(100-percent, '-') + "\t" + std::to_string(percent) + "%";	//	create bar
		countdownN = increment;	//	reset countdown
		std::cout.flush();	//	print bar

	}	//	end progress check

}	//	end countdown

string Helper::sanitiser( string s ) {	//	sanitise histogram names

	s.erase(std::remove_if(s.begin(), s.end(), []( char const& c ) -> bool { return !std::isalnum(c); } ), s.end());	//	strip invalid name characters
	return s;	//	return edited string

}	//	end sanitiser

}	//	end namespace fnt


void Helper(){};	//	allow root compilation