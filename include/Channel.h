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
		Channel( char c, Int_t n, Double_t p = 0 ) {	//	default constructor
			channelNumber = n;	//	set channel number
			type = c;	//	set channel type
			pixelNumber = p;	//	set pixel number
		}	//	end default constructor
		~Channel() {}	//	destructor
		
		void addCalibration(Int_t t, char c, vector<double> vc, vector<double> ve);	//	add calibration( offset, type, energy, efficiency )
		void addGate(char c, ULong64_t low, ULong64_t high);	//	add gate with ( type, low, high )
		Double_t adjE(double e);	//	adjust energy with a calibration( energy )
		Double_t getPixelNumber() { return pixelNumber; }	//	get pixel number
		char getType() { return type; }	//	get detector type
		Int_t getChannelNumber() { return channelNumber; }	//	get channel number
		Int_t getTOffset() { return timeOffset; };	//	get time offset
		Int_t passes(double e, char t = 'e', bool a = false);	//	check if energy is inside gate( energy, type = energy, absolute = false )
		Double_t gatesTcover(Int_t bpt) { Double_t s = 0; for( unsigned int i =0; i < gates_t.size(); i++) s += gates_t[i]->percentage(bpt); return s; };


	private:
		char type, calType;	//	detector type, calibration type
		Double_t pixelNumber;	//	pixel number
		Int_t channelNumber, timeOffset;	//	channel number, offset of channel to reference time
		std::vector<double> calTerms = {}, effTerms = {};	//	calibration, efficiency terms
		std::vector<std::pair<Gate*, char>> gates_c = {};	//	vector of coincidence gates applied on detector
		std::vector<Gate*> gates_e = {};	//	vector of energy gates applied on detector
		std::vector<Gate*> gates_t = {};	//	vector of time gates applied on detector

	};

}	//	end namespace fnt
#endif