//
//	George O'Neill @ University of York, 2020
//
//	This file holds helper functions used in various other parts of the code:
//		countdown
//		peakf
//		sanitiser
//
#include "../include/Helper.h"

namespace fnt{	//	create unique working area


void Helper::countdown(){	//	print progress

	if( !(--countdownN) ){	//	check progress

		percent++;	//	increment total bars
		std::cout << "\r" + std::string( percent, 'X' ) + std::string( 100 - percent, '-' ) + " " + std::to_string( percent ) + "%";	//	create bar
		countdownN = increment;	//	reset countdown
		std::cout.flush();	//	print bar

	}	//	end progress check

}	//	end countdown


TFitResultPtr Helper::peakf( TH1F *h, std::string s ){	//	find main peaks of spectrum

	std::string channel = std::string( h->GetName() ).substr( s.length() );	//	channel name
	TSpectrum *finder = new TSpectrum();	//	peak fitting spectrum
	TFitResultPtr fitResult;	//	store fit result parameters
	TH1F *hFitted = (TH1F *)h->Clone();	//	fitted histogram initialised as current
	const int nbins = h->GetNbinsX();	//	number of bins in energy spectra
	const double xmax = h->GetXaxis()->GetXmax(), xmin = h->GetXaxis()->GetXmin(), pCut = 0.9, binMultiplier = ( std::abs( xmin ) + std::abs( xmax ) ) / nbins;	//	store x maximum, x minimum, peak clip percentage, bin to value number
	double pWidth, fitMean, fitSigma, *source = new double[nbins], *dest = new double[nbins], *dest2 = new double[nbins], *xpeaks = 0;	//	store fit parameters, histogram arrays

	for( int i = 0; i < nbins; i++ )	//	for each bin
		source[i] = h->GetBinContent( i + 1 );	//	set array element from histogram

	if( s.find( "nrj" ) + 1 ){	//	if this is an energy spectrum

		finder->SearchHighRes( source, dest2, nbins, 1, 10, kFALSE, 3, kFALSE, 50 );	//	smooth peaks without removing background

	}	//	energy peaks
	else{	//	time peaks

		finder->SearchHighRes( source, dest, nbins, 6, 10, kFALSE, 3, kTRUE, 50 );	//	smooth peaks without removing background
		finder->SearchHighRes( dest, dest2, nbins, 6, 10, kTRUE, 3, kFALSE, 50 );	//	remove background without smoothing

	}	//	end peak filler

	if( finder->GetNPeaks() > 0 ){	//	if there are peaks

		for( int i = 0; i < nbins; i++ )	//	for each bin
			hFitted->SetBinContent( i + 1, dest2[i] );	//	fill histogram to fit

		xpeaks = finder->GetPositionX();	//	get peaks
		fitMean = xpeaks[0] * binMultiplier;	//	mean of fit
		pWidth = hFitted->GetBinContent( (int)xpeaks[0] ) * pCut;	//	store 90% of peak height
		fitSigma = ( h->FindLastBinAbove( pWidth ) - h->FindFirstBinAbove( pWidth ) ) / sqrt( -log10( pCut ) );	//	work out FWTM
		fitResult = hFitted->Fit( "gaus", "QNCS", "", ( xpeaks[0] - fitSigma ) * binMultiplier, ( xpeaks[0] + fitSigma ) * binMultiplier );	//	get fit information
		fitMean = fitResult->Parameter( 1 );	//	mean of fit
		fitSigma = 2 * fitResult->Parameter( 2 );	//	sigma of fit

		if( fitSigma > 0.075 * ( xmax - xmin ) )	//	check fwtm is less than an amount of range
			std::cout << "*** WARNING *** Channel " << channel << " is very noisy, with sigma of fit = " << fitSigma << std::endl;	//	print failure

		std::cout << s << " for channel " << channel << ":\t" << round( fitMean ) << "\tgate low:\t" << round( fitMean - fitSigma ) << "\tgate high:\t" << round( fitMean + fitSigma );	//	print neutron offset, gate

		if( sizeof( xpeaks ) / sizeof( double * ) > 1 )
			std::cout << "\tgamma " + s + " offset:\t" << xpeaks[1] * binMultiplier;	//	print gamma offset

		std::cout << std::endl;	//	end line

	}	//	end peaks check

	delete xpeaks;	//	remove peaks
	delete[] dest2;	//	remove second histogram array
	delete[] dest;	//	remove first histogram array
	delete[] source;	//	remove source histogram array
	delete hFitted;	//	remove histogram from memory and file
	delete finder;	//	remove spectrum
	return fitResult.Get();	//	send something useful back to user

}	//	end peakfinder( TH1F, string )


std::string Helper::sanitiser( std::string s ){	//	sanitise ROOT names

	s.erase( std::remove_if( s.begin(), s.end(), []( char const &c ) -> bool{ return !std::isalnum( c ); } ), s.end() );	//	strip invalid name characters
	return s;	//	return edited string

}	//	end sanitiser( string )


}	//	end namespace fnt


void Helper(){}	//	allow root compilation
