//
//	George O'Neill @ University of York, 2020
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
#include "TH1D.h"
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
		TFitResultPtr peakf( TH1D* h, std::string s );	//	find peaks
//		Int_t getPixelNE(Int_t n = 1) { return ( getPixelE(n) != nonPixel ? getPixelN(n) + n : nonPixel); }	//	get pixel number n rows above and right
		Int_t getPixelN(Int_t pixel, Int_t n = 1) { return pixel + n*8; }	//	get pixel number n rows above
		Int_t getPixelE(Int_t pixel, Int_t n = 1) { return ( ( 8 - ( pixel - 1 ) % 8 ) > n ? pixel + n : nonPixel ); }	//	get pixel number n after
		Int_t getPixelS(Int_t pixel, Int_t n = 1) { return pixel - n*8; }	//	get pixel number n rows below
		Int_t getPixelW(Int_t pixel, Int_t n = 1) { return ( ( ( pixel - 1 ) % 8 ) > ( n - 1 ) ? pixel - n : nonPixel ); }	//	get pixel number n before
		bool neighbour(Int_t pixel, Int_t other, Int_t n = 1) { return (getPixelN(pixel, n) == other) + (getPixelE(pixel, n) == other) + (getPixelS(pixel, n) == other) + (getPixelW(pixel, n) == other) + (( getPixelE(n) != nonPixel ? getPixelN(n) + n : nonPixel ) == other) + (( getPixelE(n) != nonPixel ? getPixelS(n) + n : nonPixel ) == other) + (( getPixelW(n) != nonPixel ? getPixelN(n) - n : nonPixel ) == other) + (( getPixelW(n) != nonPixel ? getPixelS(n) - n : nonPixel )); }	//	get neighbours
		std::string sanitiser(std::string s);	//	sanitise histogram names


	private:
		static const Int_t nonPixel = -2147483648;	//	initialise false pixel value
		Int_t percent = -1;	//	initialise countdown percentage marker
		ULong64_t entries, increment, countdownN;	//	number of entries for this helper, increment length, counter

	};

}	//	end namespace fnt
#endif