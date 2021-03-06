//
//	George O'Neill @ University of York, 2020
//
//	This file creates an analysis object containing important constants which are used throughout analysis
//
#ifndef ANALYSER_H
#define ANALYSER_H

//std libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

//	ROOT libraries
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TSystem.h>
#include <TStopwatch.h>

//	FNT2 libraries
#include "../include/Channel.h"
#include "../include/Helper.h"

namespace fnt{	//	create unique working area

class FNT{	//	main analysis object

public:
	FNT( std::string f = "out.root", std::string h = "histograms.root" );	//	default constructor
	~FNT(){}	//	destructor

	Channel *addChannel( char d, int n, double p = 0 );	//	add new channel
	Channel *getChannel( int i ){
		if( channels.count( i ) ) return channels.at( i ); else return NULL;
	}	//	check for existence of and return channel( channel )
	Channel *getChanB(){
		return getChannel( beamChannel );
	}	//	get beam channel
	Channel *getChanR(){
		return getChannel( thetChannel );
	}	//	get tpos channel
	Channel *getChanX(){
		return getChannel( xposChannel );
	}	//	get xpos channel
	Helper *helper;	//	required for countdown, sanitiser		
	TChain *chainer();	//	combine tree files
	TChain *getTree(){
		return tree;
	}	//	get xpos channel
	TFile *histFile(){
		return newHists;
	}	//	get histogram root file
	TH1F *getH1( std::string s ){
		return h1s.at( s );/*hs(&h1s,s);*/
	}	//	get 1-D histogram
	TH2F *getH2( std::string s ){
		return h2s.at( s );/*hs(&h2s,s);*/
	}	//	get 2-D histogram
	TH3F *getH3( std::string s ){
		return h3s.at( s );/*hs(&h3s,s);*/
	}	//	get 3-D histogram
	unsigned long long getEntries(){
		return nEntries;
	}	//	total entries in tree
	const char *getFilename(){
		return filename;
	}	//	get root file name
	const char *getHFilename(){
		return fileRootH;
	}	//	get histograms file name
	std::string getFileChannels(){
		return fileChannel;
	}	//	get channel map path
	std::string getFileGates(){
		return fileGate;
	}	//	get channel map path
	std::string getFileFiles(){
		return fileFiles;
	}	//	get path to DataTree root file list
	std::string getFileHistos(){
		return fileHisto;
	}	//	get histogram file path
	std::vector<std::string> getFolders(){
		return folders;
	}	//	get folder names
	short getMaxChannels(){
		return totalChannels;
	}	//	get highest useful channel
	short getNumChans(){
		return (short)channels.size();
	}	//	get xpos channel
	void addFriend( TTree *t ){
		tree->AddFriend( t );
	};	//	add friend to our tree
	void setMovedFlag( bool m ){
		movedFlag = m;
	}	//	set movedFlag flag
	void setTimeOffset( unsigned long long t ){
		timeOffset = t;
	}	//	set time offset
	bool addChannels();	//	add channels
	bool addGates();	//	add gates
	bool getHists();	//	get histograms
	bool histo( Float_t t, std::string n, int xbins = 1000, double xmin = 0, double xmax = 10000, int ybins = 1000, double ymin = 0, double ymax = 10000, int zbins = 1000, double zmin = 0, double zmax = 10000 );	//	add histogram


private:
	static const int totalChannels = 98;	//	highest useful channel
	static const std::string filePrefix;	//	path to files
	static const std::string fileOutPath;	//	path to output files
	static const std::string fileChannel;	//	path to channel mapping
	static const std::string fileFiles;	//	path to DataTree root file list
	static const std::string fileGate;	//	path to channel mapping
	static const std::string fileHisto;	//	path to histograms
	std::string filenames;	//	string version of path to tree
	std::string fileRootHs;	//	string version of histogram file location
	const char *filename;	//	path to tree root file
	const char *fileRootH;	//	path to histograms root file
	TChain *tree;	//	chain of files used for analysis;
	TFile *file = nullptr;	//	tree file
	TFile *newHists;	//	histograms file
	TTree *addedTree;	//	new tree to add first pass branches
	TTree *indexTree;	//	index tree for moved entries
	Bool_t movedFlag = false;	//	moved flag
	unsigned long long timeOffset = 0, nEntries = 0;	//	value for time offset, number of entries
	TBranch *bi, *bm, *bt, *be;	//	new branches for actual index, moved flag, time offset, moved entry list
	int beamChannel = -1, thetChannel = -1, xposChannel = -1;	//	essential channels
	std::unordered_map<int, Channel *> channels;	//	all channels
	std::unordered_map<std::string, TH1F *> h1s;	//	1-D histograms
	std::unordered_map<std::string, TH2F *> h2s;	//	2-D histograms
	std::unordered_map<std::string, TH3F *> h3s;	//	3-D histograms
	std::vector<Channel *> neutrons, germaniums, labr;	//	channels by detector group
	std::vector<std::string> folders;	//	folders
	std::vector<unsigned long long> movedEntry;	//	list of moved entries
/*		//	this template would be great, but it is really, unreasonably slow.  Use it for debugging histograms, otherwise stay away...
		template<class T> T hs(unordered_map<std::string, T> *h, std::string s) {	//	histogram lookup function

			try {	//	try...

				return h->at(s);	//	...to lookup histogram

			}	//	end attempt
			catch(const std::out_of_range& e) {	//	failing lookup

				std::cerr << std::endl << std::endl << "Failed at " << s << std::endl << std::endl;	//	tell user which histogram we are on
				exit(-1);	//	kill program as clean as possible

			}	//	end histogram fail

		}	//	histogram lookup
*/

};	//	end class FNT

}	//	end namespace fnt
#endif