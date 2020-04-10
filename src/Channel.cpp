//
//	George O'Neill @ University of York, 2020
//
//	This file compartmentalises all channels into one and applies a calibration
//
#include "../include/Channel.h"

namespace fnt {	//	create unique working area


void Channel::addCalibration(Int_t t, char c, vector<Double_t> v) {	//	add calibration( offset, type, terms )

	timeOffset = t;	//	set time offset
	calType = c;	//	set calibration type
	terms = v;	//	set calibration terms

}	//	end addCalibration(short c, char v, vector<Double_t> t)



void Channel::addGate(char c, ULong64_t low, ULong64_t high) {	//	add gate( type, low, high )

	Gate* g = new Gate(c, low, high);	//	create new gate

	if( c == 'e' )	gates_e.push_back( g );	//	add gate to time collection
	else if( c == 't' )	gates_t.push_back( g );	//	add gate to energy collection

}	//	end addGate( char c, ULong64_t low, ULong64_t high )



Double_t Channel::adjE(Double_t e) {	//	adjust energy with a calibration( energy )

	Double_t x = e;	//	value to return defaulted if none passed to referenced energy

	if( calType == 'P' ) {	//	polynomial Ax^N + Bx^(N-1) ... + C

		Double_t hold = 0;	//	hold all values to do the calculation

		for( UInt_t i = 0; i < terms.size(); i++ ) {	//	loop over vector of terms
			hold = terms[i] * ( pow(e, i) );	//	A*x^N
			x += hold;	//	add calculated term to total
		}	//	close loop over vector of terms

	}	//	end check of polynomial
	else if( calType == 'e' ) {	//	exponential Ae^Bx + C

		x = terms[0] * exp( terms[1] * e ) + terms[2];	//	calculate exponential correction

	}	//	end check of exponential
	else if( calType == 'l' ) {	//	logarithmic AlogB(Cx) + D

		x = terms[0] * log(terms[2] * e) / log(terms[1]) + terms[3];	//	calculate logarithmic correction

	}	//	end check of logarithm

	return x;	//	return corrected energy

}	//	end adjE( Double_t e )



bool Channel::passes(Double_t e, char t/* = 'e'*/) {	//	check if energy is inside a gate( energy, type = energy)

	vector<Gate*>* v = &(t == 'e' ? gates_e : gates_t);	//	use time or energy gate collection
	
	if( v->size() == 0 )	return true;	//	if there are no gates then it must pass
	
	for( Gate* g : *v )	if( !g->passes(e) )	return false;	//	check if energy contained in each gate
	
	return true;	//	passed for loop so must be true

}	//	end passes( Double_t e, char t )


}	//	end namespace fnt



void Channel(){};	//	allow root compilation