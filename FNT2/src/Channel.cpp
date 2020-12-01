//
//	George O'Neill @ University of York, 2020
//
//	This file compartmentalises all channels into one and applies a calibration
//
#include "../include/Channel.h"

namespace fnt {	//	create unique working area


	void Channel::addCalibration( int t, char c, std::vector<double> vc, std::vector<double> ve ) {	//	add calibration( offset, type, energy, efficiency )

		timeOffset = t;	//	set time offset
		calType = c;	//	set calibration type
		calTerms = vc;	//	set calibration terms
		effTerms = ve;	//	set efficiency terms

	}	//	end addCalibration(short c, char v, std::vector<double> t)



	void Channel::addGate( char c, unsigned long long low, unsigned long long high ) {	//	add gate( type, low, high )

		Gate *g = new Gate( c, low, high );	//	create new gate

		if ( c == 'e' )	gates_e.push_back( g );	//	add gate to energy collection
		else if ( c == 't' )	gates_t.push_back( g );	//	add gate to time collection
		else if ( c == 'c' || c >= 'A' || c <= 'Z' )	gates_c.push_back( { g, c } );	//	add gate to coincidence collection

	}	//	end addGate( char c, unsigned long long low, unsigned long long high )



	double Channel::adjE( double e ) {	//	adjust energy with a calibration( energy )

		double x = 0;	//	value to return defaulted if none passed to referenced energy

		if ( calType == 'P' ) {	//	polynomial Ax^N + Bx^(N-1) ... + C

			double hold = 0;	//	hold all values to do the calculation

			for ( unsigned int i = 0; i < calTerms.size(); i++ ) {	//	loop over vector of calibration terms
				hold = calTerms[i] * ( pow( e, i ) );	//	A*x^N
				x += hold;	//	add calculated term to total
			}	//	close loop over vector of calibration terms

		}	//	end check of polynomial
		else if ( calType == 'e' ) {	//	exponential Ae^Bx + C

			x = calTerms[0] * exp( calTerms[1] * e ) + calTerms[2];	//	calculate exponential correction

		}	//	end check of exponential
		else if ( calType == 'l' ) {	//	logarithmic AlogB(Cx) + D

			x = calTerms[0] * log( calTerms[2] * e ) / log( calTerms[1] ) + calTerms[3];	//	calculate logarithmic correction

		}	//	end check of logarithm

		return x;	//	return corrected energy

	}	//	end adjE( double e )



	int Channel::passes( double e, char t/* = 'e'*/, bool a/* = false*/ ) {	//	check if energy is inside a gate( energy, type = energy, absolute = false)

		long unsigned int s;	//	gate size holder

		if ( t != 'e' && t != 't' ) {	//	coincidence gate

			int j = 0;	//	counter
			s = gates_c.size();	//	get number of gates

			if ( !s )	//	if there are no gates
				return -1;	//	then it must pass
			else	//	we must have gates
				for ( long unsigned int i = 0; i < s; i++ )	//	so for each gate
					if ( ( gates_c[i].second == t || ( gates_c[i].second == 'N' && ( t == 'C' || t == 'M' ) ) || gates_c[i].second == 'c' ) && ++j && gates_c[i].first->passes( e ) )	//	if the type matches and gate passes
						return (int)( a ? ++i : j );	//	then tell user which one it passed at

		}	//	end coincidence check
		else {	//	energy or time gate

			std::vector<Gate *> *v = &( t == 'e' ? gates_e : gates_t );	//	use correct gate collection
			s = v->size();	//	get number of gates

			if ( !s )	//	if there are no gates
				return -1;	//	then it must pass
			else	//	we must have gates
				for ( long unsigned int i = 0; i < s; i++ )	//	so for each gate
					if ( ( *v )[i]->passes( e ) )	//	if any pass
						return (int)i + 1;	//	then tell user which one it passed at

		}	//	end gate type check

		return false;	//	if we are here we did not succeed in the for loop

	}	//	end passes( double e, char t, bool a )


}	//	end namespace fnt


void Channel() {}	//	allow root compilation