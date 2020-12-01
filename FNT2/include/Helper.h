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
#include <TFitResult.h>
#include <TH1D.h>
#include <TSpectrum.h>

namespace fnt {	//	create unique working area

class Helper {	//	object for channels

	public:
		Helper(unsigned long long n = 0, double advamount = 0.01) {	//	default constructor
			entries = n;	//	number of entries
			increment = n / (ULong64_t)(1 / advamount);	//	increment level
			countdownN = 1;	//	set countdown to 0 to print initial marker
		}	//	default constructor
		~Helper() {}	//	destructor
		
		void resetCountdown() { percent = -1; countdownN = 1; }	//	reset countdown
		void countdown();	//	print progress
		TFitResultPtr peakf( TH1D* h, std::string s );	//	find peaks
//		int getPixelNE(int n = 1) { return ( getPixelE(n) != nonPixel ? getPixelN(n) + n : nonPixel); }	//	get pixel number n rows above and right
		int getPixelN(int pixel, int n = 1) { return pixel + n*8; }	//	get pixel number n rows above
		int getPixelE(int pixel, int n = 1) { return ( ( 8 - ( pixel - 1 ) % 8 ) > n ? pixel + n : nonPixel ); }	//	get pixel number n after
		int getPixelS(int pixel, int n = 1) { return pixel - n*8; }	//	get pixel number n rows below
		int getPixelW(int pixel, int n = 1) { return ( ( ( pixel - 1 ) % 8 ) > ( n - 1 ) ? pixel - n : nonPixel ); }	//	get pixel number n before
		bool neighbour(int pixel, int other, int n = 1) { return (getPixelN(pixel, n) == other) + (getPixelE(pixel, n) == other) + (getPixelS(pixel, n) == other) + (getPixelW(pixel, n) == other) + (( getPixelE(n) != nonPixel ? getPixelN(n) + n : nonPixel ) == other) + (( getPixelE(n) != nonPixel ? getPixelS(n) + n : nonPixel ) == other) + (( getPixelW(n) != nonPixel ? getPixelN(n) - n : nonPixel ) == other) + (( getPixelW(n) != nonPixel ? getPixelS(n) - n : nonPixel )); }	//	get neighbours
		std::string sanitiser(std::string s);	//	sanitise histogram names


	private:
		static const int nonPixel = -2147483648;	//	initialise false pixel value
		int percent = -1;	//	initialise countdown percentage marker
		unsigned long long entries, increment, countdownN;	//	number of entries for this helper, increment length, counter

	};

}	//	end namespace fnt
#endif
