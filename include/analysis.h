//
//	George O'Neill @ University of York, 2020
//
//	This file is by design a simple separator for analysis
//
#ifndef ANALYSIS_H
#define ANALYSIS_H
#include "../include/FNT.h"

//	ROOT libraries needed for histogram_pretty
//#include <TCanvas.h>
//#include <THStack.h>
//#include <TLegend.h>

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
			bigTree->SetBranchAddress("timeOffset", &timeOffset);	//	set address to store time
			bigTree->SetBranchAddress("movedFlag", &movedFlag);	//	set address to store flag for moved
			bigTree->SetBranchAddress("movedEntry", &movedEntry);	//	set address to store moved entries
			n = f->getEntries();	//	set limit for for loop
			chanB = f->getChanB()->getChannelNumber();	//	get beam channel number
			chanR = f->getChanR()->getChannelNumber();	//	get theta channel number
			cX = f->getChanX();	//	get x channel
			chanX = cX->getChannelNumber();	//	get xpos channel number
			gmc = f->getMaxChannels();	//	get maximum channel
			TStopwatch TS;	//	start profiler
				doAnalysis(f);	//	do analysis
			std::cout << std::endl << std::endl << "Analysis completed in " << TS.CpuTime() << "s (which is " << TS.RealTime() << "s in real time)" << std::endl;	//	profiler output
			std::cout << std::endl;	//	output a clean line after the countdown
			histogramOperations(f);	//	do histogram operations
			newHists->Write("", TObject::kOverwrite);	//	write histograms to file
			newHists->Close();	//	close working file
			TFile::Open(f->getHFilename());	//	reopen histograms file
			std::cout << "Analysis complete!" << std::endl;	//	inform user
		}	//	end default constructor
		~analysis() {}	//	destructor

		void doAnalysis(fnt::FNT* f);	//	analysis function( FNT* tree )
		void histogramOperations(fnt::FNT* f);	//	histogram manipulation( FNT* tree )
		void objectImg(TH2D* h, Double_t w = 1);	//	object image maker( TH2D* histogram, Double_t weight = 1 )
//		void histogramPretty(TTree* bigTree);	//	colourful hit pattern
		void setR(Int_t n);	//	set r positions( Int_t nrj )
		void setX(Int_t n);	//	set x positions( Int_t nrj )
		void setTrig();	//	set trigonometric variables


	private:
		const Double_t tableSize = 225.;	//	squared diameter of table in centimetres
		const Double_t objStep = 0.1;	//	object step size for smearing
		static const Int_t bpt = 400000;	//	time between beam pulses in ps ticks
		TChain* bigTree;	//	main tree
		TFile* newHists;	//	histogram root file
		Channel* cX;	//	allocate position channel pointer
		Bool_t movedFlag, cxFlag;	//	moved flag, in the interesting object region flag
		Double_t rcal, xcal, sinrpos, cosrpos, XsinR, XcosR, tableSizeDiff;	//	trigonometric holders
		Int_t chanB, chanR, chanX, il, xposition = -1, rposition = -1;	//	integer conversion of beam time, theta channel, x channel, label
		UChar_t label;	//	initialise variable to store labels
		UInt_t gmc, nrj, nrj2;	//	initialise variables to store maximum channel, nrj, nrj2
		ULong64_t timeb, timeOffset, n;	//	initialise time branch variable, time offset, number of entries
		std::vector<ULong64_t> *movedEntry = 0;	//	list of moved entries

};	//	end class analysis

}	//	end namespace fnt
#endif