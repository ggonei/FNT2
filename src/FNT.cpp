//
//	George O'Neill @ University of York, 2020/02/18
//
//	This file creates an analysis object containing important constants which are used throughout analysis
//
#include "../include/FNT.h"

namespace fnt {	//	create unique working area

const std::string FNT::filePrefix = "/mnt/d/fnt/macros/git/FNT2/inputs/";
const std::string FNT::fileChannel = filePrefix + "channels.txt";	//	path to channel mapping
const std::string FNT::fileFiles = filePrefix + "files.txt";	//	path to calibration file
const std::string FNT::fileGate = filePrefix + "gates.txt";	//	path to channel mapping
const std::string FNT::fileHisto = filePrefix + "histos.txt";	//	path to channel mapping

FNT::FNT(const char* f, const char* h)	{	//	default constructor

	ULong64_t timeb, timerunstart = -1, toffset = 0, tprevious = 0;	//	initialise time branch variable and offsets

	std::cout << "Tree output is in " << f << " and histograms will be saved in " << h << std::endl;	//	tell user some basic info
	filename = f;	//	set tree filename
	fileRootH = h;	//	set histogram filename

	addChannels();	//	add all channels list
	addGates();	//	add gates

	if( gSystem->AccessPathName(filename) ) {	//	if file does not exist

		std::cout << "No tree exists!  Creating..." << std::endl;	//	inform user we are making a new tree
		file = new TFile(filename, "RECREATE");	//	save output to new file
		addedTree = new TTree("newTree", "newTree");	//	editable tree
		bp = addedTree->Branch("timeLastB", &(timeLastB), "timeLastB/l");	//	create new branch for last beam pulse
		bx = addedTree->Branch("xPosition", &(xPosition), "xPosition/i");	//	create new branch for x position
		br = addedTree->Branch("rPosition", &(rPosition), "rPosition/i");	//	create new branch for theta position
		tree = chainer();	//	create tree
		tree->AddFriend(addedTree);	//	add expanded branches to existing data
		Int_t il;	//	integer conversion of label
		UChar_t label;	//	initialise variable to store labels
		UInt_t nrj;	//	initialise variable to store nrj
		tree->SetBranchAddress("label", &label);	//	set address to store label
		tree->SetBranchAddress("nrj", &nrj);	//	set address to store energy
		tree->SetBranchAddress("time", &timeb);	//	set address to store time
		n = tree->GetEntries();	//	initialise number of entries
		Int_t chanb = getChanB(), chanr = getChanR(), chanx = getChanX();	//	get channel positions once
		helper = new Helper(n);	//	required for countdown

		for( ULong64_t i = 0; i < n; i++ ) {	//	loop over all entries

			helper->countdown();	//	print progress
			tree->GetEntry(i);	//	grab label
			il = (Int_t) label;	//	convert label to integer
			if( timeb < timerunstart ) {	//	look for clock reset
				
				toffset = tprevious;	//	hold offset time
				timerunstart = timeb;	//	reset t0
				std::cout << "Clock has been reset at time " << timerunstart << " on entry number " << i << ", changing offset from " << toffset << " to " << tprevious << std::endl;	//	inform user
				
			}	//	end clock reset check

			if( il == chanb ) setTimeLastB(toffset + timeb);	//	update beam time if there is one

			else if( il == chanr ) setRPosition(nrj);	//	update theta if there is one

			else if( il == chanx ) setXPosition(nrj);	//	update xpos if there is one

			addedTree->Fill();	//	fill new branches
			tprevious = timeb;	//	store current event time for next event clock comparison
			
		}	//	end loop over all entries

		std::cout << std::endl;	//	output a clean line after the countdown

		addedTree->Write("", TObject::kOverwrite);	//	Write new tree to file
		tree->Write("", TObject::kOverwrite);	//	Write main tree to file

	}	//	file exists

	else {	//	file does exist

		std::cout << "Tree exists, using: " << filename << std::endl;	//	inform user we are using existing tree
		file = TFile::Open(filename);	//	load previous file
		tree = (TChain*)file->Get("DataTree");	//	load tree
		addedTree = tree->GetFriend("newTree");	//	load new tree
		tree->SetBranchAddress("time", &timeb);	//	set address to store time
		n = tree->GetEntries();	//	initialise number of entries
		helper = new Helper(n);	//	required for countdown
		
	}	//	end tree load

	tree->GetEntry(0);	//	get first entry
	std::cout << "There are " << n << " entries, with a minimum timestamp of " << timeb << " (" << timeb/(ULong64_t)pow(10,floor(log10(timeb))) << "e" << floor(log10(timeb))+1 << ") and a maximum timestamp of ";	//	output minimum timestamp
	tree->GetEntry(n-1);	//	get last entry
	std::cout << timeb << " (" << timeb/(ULong64_t)pow(10,floor(log10(timeb))) << "e" << floor(log10(timeb))+1 << ")" << std::endl;	//	output maximum timestamp
	getHists();	//	create histograms

}	//	end default constructor



Channel* FNT::addChannel( char d, int n, int p /* = -1 */ ) {	//	add channel( type, channel, pixel )

	Channel* x = new Channel(d, n, p);	//	create new channel
	channels.insert({n, x});	//	add a new channel to the map

	if( d == 'n' ) neutrons.push_back(x);	//	add channel to neutrons
	else if( d == 'g' ) germaniums.push_back(x);	//	add channel to germaniums
	else if( d == 'l' ) labr.push_back(x);	//	add channel to LaBr3s
	
	return x;

}



bool FNT::addChannels() {	//	add channels
	
	char d, c;	//	detector type, calibration type
	double n;	//	calibration term
	std::ifstream s( getFileChannels() );	//	open file
	int t;	//	time offset
	short i, p;	//	channel, pixel number
	std::string x = "";	//	hold lines
	std::stringstream h;	//	hold remainer of line in stream
	std::vector<double>	v;	//	hold calibration terms

	s.ignore(std::numeric_limits<std::streamsize>::max(), '\n');	//	skip opening comment line

	while( getline(s, x) ) {	//	while we have data...

		v.clear();	//	reset calibration term vector
		h.str(x);	//	fill holding stream
		h.clear();	//	reset stream pointer
		h >> i >> d >> p >> t >> c;	//	hold channel number, detector type, pixel number, time offset, and calibration type
		while( h >> n )	v.push_back(n);	//	rebuild double string
		Channel* newc = addChannel(d, i, p);	//	add channel
		newc->addCalibration(t, c, v);	//	add calibration to channel with ( offset, type, terms )

		if( d == 't' ) setChanB( i );	//	set brfp channel
		else if( d == 'r' ) setChanR( i );	//	set thet channel
		else if( d == 'x' ) setChanX( i );	//	set xpos channel

	}	//	end data read in

	return true;	//	success

}	//	end addChannels



bool FNT::addGates() {	//	add gates
	
	char c;	//	gate type
	std::ifstream s( getFileGates() );	//	open file
	ULong64_t low, high;	//	low, high values
	short i;	//	channel number
	std::string x = "";	//	hold lines
	
	s.ignore(std::numeric_limits<std::streamsize>::max(), '\n');	//	skip opening comment line
	
	while( getline(s, x) ) {	//	while we have data...
	
		s >> i >> c >> low >> high;	//	fill channel, type, low, and high
		if( getChannel(i) ) getChannel(i)->addGate(c, low, high);	//	add gate to channel with ( type, low, high )
		else std::cout << "Gate for channel " << i << " exists, but channel " << i << " doesn't exist!" << std::endl;	//	no channel to add calibration to
	
	}	//	end data read in
	
	return true;	//	success
	
}	//	end addGates



TChain* FNT::chainer() {	//	create new trees

	TChain* c = new TChain("DataTree", "DataTree");	//	create chain
	std::ifstream s( getFileFiles() );	//	open file
	char x[4096];	//	hold lines
	
	s.ignore(std::numeric_limits<std::streamsize>::max(), '\n');	//	skip opening comment line
	
	while( s.getline(x, sizeof(x), '\n') ) {	//	while we have data...
	
		if( x[strlen(x)-1] == '\r')	x[strlen(x)-1] = '\0';	//	remove carriage returns

		std::cout << "Adding " << x << " to tree..." << std::endl;	//	inform user of file
		c->Add(x);	//	add file to our chain
	
	}	//	end data read in
	
	return c;	//	return pointer to chain

}	//	end chainer



bool FNT::getHists() {	//	get histograms

	newHists = new TFile(fileRootH, "RECREATE");	//	create histogram file
	std::ifstream s( getFileHistos() );	//	open file
	Int_t t;	//	histogram type
	Int_t xbins = -1, ybins = -1;	//	histogram bin counts
	Double_t xmin = 0, xmax = 0, ymin = 0, ymax = 0;	//	histogram axis minmax
	std::string n = "", x = "";	//	hold name, lines
	
	s.ignore(std::numeric_limits<std::streamsize>::max(), '\n');	//	skip opening comment line
	
//	histo( 0, "test", 10, 0, 10 );
	
	while( getline(s, x) ) {	//	while we have data...
	
		s >> t >> n >> xbins >> xmin >> xmax;	//	fill values
		
		if( t < 1 )	{	//	if 1-D histograms
			
			if( t < 0 ) {	//	if multiple histograms
				
				folders.push_back(n);	//	use prefix to add to folder names
				
			}	//	end multiple histogram check
			
			t = abs(t);
			
			for( Int_t i = 0; i <= (t < getMaxChannels() ? t : getMaxChannels()); i++ )	{	//	loop over all histograms

				histo( 0, n + std::to_string(i), xbins, xmin, xmax );	//	create 1-d histogram

			}
			
		}
		else if( t == 1 ) {	// if a 2-D histogram
			
			s >> ybins >> ymin >> ymax;	//	fill remaining values
			histo( t, n, xbins, xmin, xmax, ybins, ymin, ymax );	//	create 2-d histogram
			
		}
	
	}	//	end data read in
	
	newHists->Write("", TObject::kOverwrite);	//	write histograms to file
	return true;	//	success

}	//	end getHists()



bool FNT::histo( Int_t t, std::string n, Int_t xbins, Double_t xmin, Double_t xmax, Int_t ybins, Double_t ymin, Double_t ymax ) {	//	add histogram

	std::string s = helper->sanitiser(n);	//	copy string for use for name
	if( xbins < 1 )	xbins = std::abs(xmin + xmax)/2;	//	set bin count to one per channel

	if( t < 1 )	h1s.insert({s, new TH1F( s.c_str(), n.c_str(), xbins, xmin, xmax )});	//	add 1-d histogram to vector
	
	else {	// if a 2-D histogram
		
		if( ybins < 1 )	ybins = std::abs(ymin + ymax)/2;	//	set bin count to one per channel
		
		h2s.insert({s, new TH2F( s.c_str(), n.c_str(), xbins, xmin, xmax, ybins, ymin, ymax )});	//	add 2-d histogram to vector
		
	}
	
	return true;	//	success

}	//	end adding histogram


}	//	end namespace fnt



void FNT(){};	//	allow root compilation