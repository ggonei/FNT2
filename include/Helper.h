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

//ROOT libraries
#include "TFitResult.h"
#include "TH1.h"
#include "TSpectrum.h"

namespace fnt {	//	create unique working area

class Helper {	//	object for channels

	public:
		Helper(ULong64_t n = 0, double advamount = 0.01) {	//	default constructor
			entries = n;	//	number of entries
			increment = advamount * n;	//	increment level
			countdownN = 1;	//	set countdown to 0 to print initial marker
		}	//	default constructor
		~Helper() {}	//	destructor
		
		void resetCountdown() { percent = -1; countdownN = 1; }	//	reset countdown
		void countdown();	//	print progress
		void peakf( TH1F* h, std::string s );	//	find peaks
		std::string sanitiser(std::string s);	//	sanitise histogram names


	private:
		Int_t percent = -1;	//	initialise countdown percentage marker
		ULong64_t entries, increment, countdownN;	//	number of entries for this helper, increment length, counter

	};

}	//	end namespace fnt
#endif