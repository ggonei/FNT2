//
//	George O'Neill @ University of York, 2020
//
//	This file is by design a simple separator for analysis
//
#ifndef ANALYSIS_H
#define ANALYSIS_H
#include "../include/FNT.h"

//	C++ libraries
#include <deque>

//	ROOT libraries
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"

namespace fnt {	//	create unique working area

class analysis {	//	main analysis object

	public:
		analysis(fnt::FNT* f)	{	//	default constructor
			std::cout << "Analysis object " << f << " has compiled successfully, continuing bespoke analysis from analysis.cpp..." << std::endl;	//	informative message
			newHists = f->histFile();	//	get histogram root file
			bigTree = f->getTree();	//	get big tree
			bigTree->SetBranchAddress("label", &label);	//	set address to store label
			bigTree->SetBranchAddress("time", &timeb);	//	set address to store time
			bigTree->SetBranchAddress("nrj", &nrj);	//	set address to store energy
			bigTree->SetBranchAddress("nrj2", &nrj2);	//	set address to store energy_slow
			bigTree->SetBranchAddress("timeLastB", &timeLastB);	//	set address to store last beam time recorded
			bigTree->SetBranchAddress("xPosition", &xPos);	//	set address to store last known x
			bigTree->SetBranchAddress("rPosition", &rPos);	//	set address to store last known theta
			bigTree->SetBranchAddress("timeOffset", &timeOffset);	//	set address to store time
			n = f->getEntries();	//	set limit for for loop
			bpt = f->getBPT();	//	get beam pulse time
			coincWind = f->getCoincWind();	//	get coincidence window time
			cR = f->getChanR();	//	get theta channel
			cX = f->getChanX();	//	get xpos channel
			chanR = cR->getChannelNumber();	//	get theta channel number
			chanX = cX->getChannelNumber();	//	get xpos channel number
			gmc = f->getMaxChannels();	//	get maximum channel
			do_analysis(f);	//	do analysis
			std::cout << std::endl;	//	output a clean line after the countdown
			histogram_operations(f);	//	do histogram operations
			newHists->Write("", TObject::kOverwrite);	//	write histograms to file
			newHists->Close();	//	close working file
			TFile::Open(f->getHFilename());	//	reopen histograms file
			std::cout << "Analysis complete!" << std::endl;	//	inform user
		}	//	end default constructor
		~analysis() {}	//	destructor

		void do_analysis(fnt::FNT* f);	//	analysis function
		void histogram_operations(fnt::FNT* f);	//	histogram manipulation
		void histogram_pretty(TTree* bigTree);	//	colourful hit pattern


	private:
		TChain* bigTree;	//	main tree
		TFile* newHists;	//	histogram root file
		Channel *cR, *cX;	//	position channels
		Int_t bpt, chanR, chanX, coincWind, il;	//	integer conversion of beam time, theta channel, x channel, label
		UChar_t label;	//	initialise variable to store labels
		UInt_t gmc, nrj, nrj2, xPos, rPos;	//	initialise variables to store maximum channel, nrj, nrj2, x position, theta position
		ULong64_t timeb, timeLastB, timeOffset, n;	//	initialise time branch variable, progress counter, time offset, and number of entries


};	//	end class analysis

}	//	end namespace fnt
#endif