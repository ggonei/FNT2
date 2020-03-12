//
//	George O'Neill @ University of York, 2020/02/18
//
//	This file creates an analysis object containing important constants which are used throughout analysis
//
#ifndef ANALYSER_H
#define ANALYSER_H

//std libraries
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

//	ROOT libraries
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TSystem.h>

//	FNT2 libraries
#include "../include/Channel.h"
#include "../include/Helper.h"

namespace fnt {	//	create unique working area

class FNT {	//	main analysis object

	public:
		FNT(const char* f = "out.root", const char* h = "histograms.root");	//	default constructor
		~FNT() {}	//	destructor
		
		Channel* addChannel( char d, int n, int p = -1 );	//	add new channel
		Channel* getChannel(int i) { if( channels.count(i) ) return channels.at(i); else return NULL; }	//	check for existence of and return channel( channel )
		Helper* helper;	//	required for countdown, sanitiser
		
		TChain* chainer();	//	combine tree files
		TChain* getTree() { return tree; }	//	get xpos channel
		TTree* newTree() { return addedTree; }	//	get new tree
		TFile* histFile() { return newHists; }	//	get histogram root file
		TH1F* getH1(string s) { return h1s.at(s); }	//	return histogram
		TH2F* getH2(string s) { return h2s.at(s); }	//	return histogram

		Int_t getBPT() { return beamPulseTime; }	//	get time between beam pulses
		Int_t getChanB() { return beamChannel; }	//	get beam channel
		Int_t getChanR() { return thetChannel; }	//	get tpos channel
		Int_t getChanX() { return xposChannel; }	//	get xpos channel
		const char* getFilename() { return filename; }	//	get root file name
		const char* getHFilename() { return fileRootH; }	//	get histograms file name
		std::string getFileCalibration() { return fileCalibration; }	//	get calibration file path
		std::string getFileChannels() { return fileChannel; }	//	get channel map path
		std::string getFileGates() { return fileGate; }	//	get channel map path
		std::string getFileFiles() { return fileFiles; }	//	get calibration file path
		std::string getFileHistos() { return fileHisto; }	//	get histogram file path
		std::vector<std::string> getFolders() { return folders; }	//	get folder names
		short getMaxChannels() { return totalChannels; }	//	get highest useful channel
		short getNumChans() { return channels.size(); }	//	get xpos channel

		void addFriend(TTree* t) { tree->AddFriend(t); };	//	add friend to our tree
		void setChanB(int n) { beamChannel = n; }	//	set beam clock( channel )
		void setChanR(int n) { thetChannel = n; }	//	set theta position( channel )
		void setChanX(int n) { xposChannel = n; }	//	set table x position( channel )
		void setRPosition(UInt_t nrj) { rPosition = nrj; }	//	set time of last theta position
		void setTimeLastB(ULong64_t timeb) { timeLastB = timeb; }	//	set time of last beam pulse
		void setXPosition(UInt_t nrj) { xPosition = nrj; }	//	set time of last x position

		bool addChannels();	//	add channels
		bool addGates();	//	add gates
		bool getHists();	//	get histograms
		bool histo(Int_t t, std::string n, Int_t xbins = 1000, Double_t xmin = 0, Double_t xmax = 10000, Int_t ybins = 1000, Double_t ymin = 0, Double_t ymax = 10000 );	//	add histogram
		ULong64_t getEntries() { return n; }	//	total entries in tree


	private:
		const char* filename;	//	path to tree root file
		const char* fileRootH;	//	path to histograms root file
		static const std::string filePrefix;	//	path to files
		static const std::string fileCalibration;	//	path to calibration file
		static const std::string fileChannel;	//	path to channel mapping
		static const std::string fileFiles;	//	path to calibration file
		static const std::string fileGate;	//	path to channel mapping
		static const std::string fileHisto;	//	path to histograms
		static const int beamPulseTime = 400000;	//	time between beam pulses in picoseconds
		static const int totalChannels = 98;	//	highest useful channel
		TChain* tree;	//	chain of files used for analysis;
		TFile* file;	//	tree file
		TFile* newHists;	//	histograms file
		TTree* addedTree;	//	new tree to add branches to
		UInt_t rPosition = 0, xPosition = 0;	//	variables to store theta and x
		ULong64_t timeLastB = 0, n = 0;	//	value for most recent beam pulse, number of entries
		TBranch *bp, *br, *bx;	//	new branches for last beam pulse, theta position, x position
		Int_t beamChannel = -1, thetChannel = -1, xposChannel = -1;	//	essential channels
		std::unordered_map<int, Channel*> channels;	//	all channels
		std::unordered_map<std::string, TH1F*> h1s;	//	1-D histograms
		std::unordered_map<std::string, TH2F*> h2s;	//	2-D histograms
		std::vector<Channel*> neutrons, germaniums, labr;	//	channels by detector group
		std::vector<std::string> folders;	//	folders


};	//	end class FNT

}	//	end namespace fnt
#endif