//
//	George O'Neill @ University of York, 2020
//
//	This file compartmentalises all channels into one and applies a calibration
//
#ifndef CHANNELS_H
#define CHANNELS_H

//std libraries
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

//FNT2 libraries
#include "../include/Gate.h"

namespace fnt{	//	create unique working area

class Channel{	//	object for channels

public:
	Channel( char c, int n, double p = 0 ){	//	default constructor
		channelNumber = n;	//	set channel number
		type = c;	//	set channel type
		pixelNumber = p;	//	set pixel number
	}	//	end default constructor
	~Channel(){}	//	destructor

	void addCalibration( int t, char c, std::vector<double> vc, std::vector<double> ve );	//	add calibration( offset, type, energy, efficiency )
	void addGate( char c, unsigned long long low, unsigned long long high );	//	add gate with ( type, low, high )
	double adjE( double e );	//	adjust energy with a calibration( energy )
	double getPixelNumber(){
		return pixelNumber;
	}	//	get pixel number
	char getType(){
		return type;
	}	//	get detector type
	int getChannelNumber(){
		return channelNumber;
	}	//	get channel number
	int getTOffset(){
		return timeOffset;
	}	//	get time offset
	int passes( double e, char t = 'e', bool a = false );	//	check if energy is inside gate( energy, type = energy, absolute = false )
	double gatesTcover( int bpt ){
		double s = 0; for( unsigned int i = 0; i < gates_t.size(); i++ ) s += gates_t[i]->percentage( bpt ); return s;
	};
	Gate *getGate( unsigned long int g, char t = 'e' ){
		return ( t == 'e' ? gates_e[g >= gates_e.size() ? 0 : g] : gates_t[g >= gates_t.size() ? 0 : g] );
	}	//	get specific gate


private:
	char type, calType;	//	detector type, calibration type
	double pixelNumber;	//	pixel number
	int channelNumber, timeOffset;	//	channel number, offset of channel to reference time
	std::vector<double> calTerms = {}, effTerms = {};	//	calibration, efficiency terms
	std::vector<std::pair<Gate *, char>> gates_c = {};	//	std::vector of coincidence gates applied on detector
	std::vector<Gate *> gates_e = {};	//	std::vector of energy gates applied on detector
	std::vector<Gate *> gates_t = {};	//	std::vector of time gates applied on detector

};

}	//	end namespace fnt
#endif