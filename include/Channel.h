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

namespace fnt {	//	create unique working area

class Channel {	//	object for channels

	public:
		Channel( char c, Int_t n, Int_t p = -1 ) {	//	default constructor
			channelNumber = n;	//	set channel number
			type = c;	//	set channel type
			pixelNumber = p;	//	set pixel number
		}	//	end default constructor
		~Channel() {}	//	destructor
		
		void addCalibration(Int_t t, char c, vector<double> v);	//	add calibration( offset, type, terms )
		void addGate(char c, ULong64_t low, ULong64_t high);	//	add gate with ( type, low, high )
		Double_t adjE(double e);	//	adjust energy with a calibration( energy )
		Int_t getChannelNumber() { return channelNumber; }	//	get channel number
		Int_t getPixelNumber() { return pixelNumber; }	//	get pixel number
		Int_t getTOffset() { return timeOffset; };	//	get time offset
		bool passes(double e, char t = 'e');	//	check if energy is inside gate( energy, type = energy )


	private:
		char type, calType;	//	detector type, calibration type
		Int_t channelNumber, pixelNumber, timeOffset;	//	channel number, pixel number, offset of channel to reference time
		std::vector<double> terms = {};	//	calibration terms
		std::vector<Gate*> gates_e = {};	//	vector of energy gates applied on detector
		std::vector<Gate*> gates_t = {};	//	vector of time gates applied on detector

	};

}	//	end namespace fnt
#endif