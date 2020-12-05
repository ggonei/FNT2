//
//	George O'Neill @ University of York, 2020
//
//	This file is by design a simple separator for analysis
//
#ifndef ANALYSIS_H
#define ANALYSIS_H
#include "../include/FNT.h"

//	ROOT libraries needed for histogramPretty
//#include <TCanvas.h>
//#include <THStack.h>
//#include <TLegend.h>

namespace fnt{	//	create unique working area

class analysis{	//	main analysis object

public:
	analysis( FNT *f ){	//	default constructor
		std::cout << "Analysis object " << f << " has compiled successfully, continuing bespoke analysis from analysis.cpp..." << std::endl;	//	informative message
		newHists = f->histFile();	//	get histogram root file
		bigTree = f->getTree();	//	get big tree
		bigTree->SetBranchAddress( "label", &label );	//	set address to store label
		bigTree->SetBranchAddress( "time", &timeb );	//	set address to store time
		bigTree->SetBranchAddress( "nrj", &nrj );	//	set address to store energy
		bigTree->SetBranchAddress( "nrj2", &nrj2 );	//	set address to store energy_slow
		bigTree->SetBranchAddress( "timeOffset", &timeOffset );	//	set address to store time
		bigTree->SetBranchAddress( "movedFlag", &movedFlag );	//	set address to store flag for moved
		bigTree->SetBranchAddress( "movedEntry", &movedEntry );	//	set address to store moved entries
		n = f->getEntries();	//	set limit for for loop
		cR = f->getChanR();	//	get r channel
		cX = f->getChanX();	//	get x channel
		chanB = f->getChanB()->getChannelNumber();	//	get beam channel number
		chanR = cR->getChannelNumber();	//	get theta channel number
		chanX = cX->getChannelNumber();	//	get xpos channel number
		rmin = (double)cR->getGate( 0 )->min();	//	minimum part of gate for rotation calibration
		rdiv = ( (double)cR->getGate( 0 )->max() - rmin ) / 6.28318530718;	//	division part of calibration for rotation
		xmin = (double)cX->getGate( 1 )->min();	//	minimum part of gate for position calibration
		xdiv = ( (double)cX->getGate( 1 )->max() - xmin ) / ( 2 * xmoving );	//	division part of calibration for position
		gmc = f->getMaxChannels();	//	get maximum channel
		gammaBrowser = f->getH3( "Objectgamma0" );
		TStopwatch TS;	//	start profiler
		doAnalysis( f );	//	do analysis
		std::cout << std::endl << std::endl << "Analysis completed in " << TS.CpuTime() << "s (which is " << TS.RealTime() << "s in real time)" << std::endl;	//	profiler output
		std::cout << std::endl;	//	output a clean line after the countdown
		histogramOperations( f );	//	do histogram operations
		newHists->Write( "", TObject::kOverwrite );	//	write histograms to file
		newHists->Close();	//	close working file
		TFile::Open( f->getHFilename() );	//	reopen histograms file
		std::cout << "Analysis complete!" << std::endl;	//	inform user
	}	//	end default constructor
	~analysis(){}	//	destructor

	void doAnalysis( FNT *f );	//	analysis function( FNT* tree )
	void histogramOperations( fnt::FNT *f );	//	histogram manipulation( FNT* tree )
//		void histogramPretty(TTree* bigTree);	//	colourful hit pattern
	void objectImg( TH2F *h, double w = 1 );	//	object image maker( TH2F* histogram, double weight = 1 )
	void objectImg( double e, double w = 1 );	//	object image maker( double energy, double weight = 1 )
	void setR( int n );	//	set r positions( int nrj )
	void setX( int n );	//	set x positions( int nrj )
	void setTrig();	//	set trigonometric variables


private:
	const double tableSize = 225., xmoving = 14.5;	//	squared diameter of table in centimetres, actual table moving distance
	const double objStep = 0.1;	//	object step size for smearing
	static const int bpt = 399999;	//	time between beam pulses in ps ticks
	TChain *bigTree;	//	main tree
	TFile *newHists;	//	histogram root file
	TH3F *gammaBrowser;	//	fetch gamma browser
	Channel *cX, *cR;	//	allocate position channel pointer, rotation channel pointer
	Bool_t movedFlag, cxFlag, crFlag;	//	moved flag, in the interesting object region flag for x and r
	double rcal, rmin, rdiv, xcal, xmin, xdiv, sinrpos, cosrpos, XsinR, XcosR, tableSizeDiff;	//	trigonometric holders
	int chanB, chanR, chanX, il, xposition = -1, rposition = -1;	//	beam channel, theta channel, x channel, label
	UChar_t label;	//	initialise variable to store labels
	unsigned int gmc, nrj, nrj2;	//	initialise variables to store maximum channel, nrj, nrj2
	unsigned long long timeb, timeOffset, n;	//	initialise time branch variable, time offset, number of entries
	std::vector<unsigned long long> *movedEntry = 0;	//	list of moved entries

};	//	end class analysis

}	//	end namespace fnt
#endif
