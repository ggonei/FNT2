//
//	George O'Neill @ University of York, 2020/03/12
//
//	This file holds helper functions used in various other parts of the code
//
#ifndef HELPER_H
#define HELPER_H

//std libraries
#include <algorithm>
#include <iostream>

namespace fnt {	//	create unique working area

class Helper {	//	object for channels

	public:
		Helper(ULong64_t n = 0, double advamount = 0.01) { entries = n; increment = advamount * n, countdownN = increment; }	//	default constructor
		~Helper() {}	//	destructor
		
		void resetcountdown() { percent = 0; countdownN = increment; }	//	reset countdown
		void countdown();	//	print progress
		string sanitiser(string s);	//	sanitise histogram names


	private:
		Int_t percent = 0;	//	set countdown percentage
		ULong64_t entries, increment, countdownN;	//	number of entries for this helper, increment length, counter

	};

}	//	end namespace fnt
#endif