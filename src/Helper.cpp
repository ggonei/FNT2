//
//	George O'Neill @ University of York, 2020/03/12
//
//	This file holds helper functions used in various other parts of the code:
//		countdown
//		peakfinder
//		sanitiser
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


void Helper::peakf( TH1F* h, std::string s ) {	//	find main peaks of spectrum

	ULong64_t binMultiplier;	//	multiplier bins
	std::string channel = std::string(h->GetName()).substr(s.length());	//	channel name
	TSpectrum* finder = new TSpectrum();	//	peak fitting spectrum
	TFitResultPtr fitResult;	//	store fit result parameters
	TH1F* hFitted = h;	//	fitted histogram initialised as current
	const Int_t nbins = h->GetNbinsX();	//	number of bins in energy spectra
	Int_t tWidth = 50000, xmax = h->GetXaxis()->GetXmax(), xmin = h->GetXaxis()->GetXmin();	//	store time offset, x maximum, x minimum
	Double_t fitMean, fitSigma, source[nbins], dest[nbins], *xpeaks;	//	store fit parameters, histogram arrays

	for( Int_t i = 0; i < nbins; i++ )	source[i] = h->GetBinContent(i + 1);	//	set source array from histogram

	finder->SearchHighRes(source, dest, nbins, 10, 10, kFALSE, 2, kTRUE, 10);	//	find peaks

	if( finder->GetNPeaks() > 0 )	{	//	if there are peaks

		for (Int_t i = 0; i < nbins; i++ )	hFitted->SetBinContent(i + 1, dest[i]);	//	fill histogram to fit

		binMultiplier = (abs(xmin) + abs(xmax)) / h->GetNbinsX();	//	convert bin number to time
		xpeaks = finder->GetPositionX();	//	get peaks
		fitMean = xpeaks[0] * binMultiplier;	//	mean of fit
		fitResult = hFitted->Fit("gaus", "QNCS", "", fitMean - tWidth, fitMean + tWidth);	//	get fit information
		fitMean = fitResult->Parameter(1);	//	mean of fit
		fitSigma = 2*fitResult->Parameter(2);	//	sigma of fit
		
		if( fitSigma < 0.075 * (xmax - xmin) )	{	//	check fwtm is less than an amount of range

			std::cout << "Neutron " << s << " for channel " << channel << ":\t" << round(fitMean) << "\tgate low:\t" << round(fitMean - fitSigma) << "\tgate high:\t" << round(fitMean + fitSigma);	//	print neutron offset, gate
			if( sizeof(xpeaks)/sizeof(Double_t*) > 1)	std::cout << "\tgamma " + s + " offset:\t" << xpeaks[1] * binMultiplier;	//	print gamma offset

		}	else std::cout << "Channel " << channel << " is too noisy";	//	print failure
		std::cout << std::endl;	//	end line

	}	//	end peaks check

}	//	end peakfinder( TH1F )


std::string Helper::sanitiser( std::string s ) {	//	sanitise histogram names

	s.erase(std::remove_if(s.begin(), s.end(), []( char const& c ) -> bool { return !std::isalnum(c); } ), s.end());	//	strip invalid name characters
	return s;	//	return edited string

}	//	end sanitiser( string )


}	//	end namespace fnt


void Helper(){};	//	allow root compilation