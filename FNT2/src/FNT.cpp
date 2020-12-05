//
//	George O'Neill @ University of York, 2020
//
//	This file creates an analysis object containing important constants which are used throughout analysis
//
#include "../include/FNT.h"

namespace fnt{	//	create unique working area

const std::string FNT::filePrefix = "/mnt/c/Users/ggon0/source/repos/FNT2/inputs/";	//	file path for file inputs
const std::string FNT::fileOutPath = "/mnt/d/fnt/sorted/";	//	file path for saving histogram names
const std::string FNT::fileChannel = filePrefix + "channels.dat";	//	path to channel mapping
const std::string FNT::fileFiles = filePrefix + "files.dat";	//	path to DataTree root file list
const std::string FNT::fileGate = filePrefix + "gates.dat";	//	path to channel mapping
const std::string FNT::fileHisto = filePrefix + "histos.dat";	//	path to channel mapping

FNT::FNT( std::string f, std::string h ){	//	default constructor

	unsigned long long jEntry, timeb, timerunstart = -1, tprevious = 0, tpreviousnoO = 0, thist = 0, tOffset = 0;	//	initialise entry store, time branch variable, run starting time, previous entry time, previous time without offset, this combined time time offset
	filenames = fileOutPath + f;	//	create path to ROOT file
	fileRootHs = fileOutPath + h;	//	create path to histograms
	filename = filenames.c_str();	//	set tree filename
	fileRootH = fileRootHs.c_str();	//	set histogram filename
	std::cout << "Tree output is in " << filename << " and histograms will be saved in " << fileRootH << std::endl;	//	tell user some basic info
	addChannels();	//	add all channels list
	addGates();	//	add gates
	TStopwatch TS;	//	start profiler

	if( gSystem->AccessPathName( filename ) ){	//	if file does not exist

		unsigned long long movedCtr = 0, movedSize;	//	moved entry counter, size of entry std::vector
		std::unordered_map<unsigned long long, std::vector<std::pair<unsigned long long, unsigned int>>> eindex;	//	index holder
		std::cout << "No tree exists!  Creating at " << filename << std::endl;	//	inform user we are making a new tree
		file = new TFile( filename, "RECREATE" );	//	save output to new file
		addedTree = new TTree( "newTree", "newTree" );	//	editable tree
		indexTree = new TTree( "iTree", "iTree" );	//	index tree
		bm = addedTree->Branch( "movedFlag", &( movedFlag ), "movedFlag/O" );	//	create new branch for moved flag
		bt = addedTree->Branch( "timeOffset", &( timeOffset ), "timeOffset/l" );	//	create new branch for time offset
		be = indexTree->Branch( "movedEntry", &( movedEntry ) );	//	create new branch from ordered std::vector of eindex values
		tree = chainer();	//	create tree
		tree->AddFriend( addedTree );	//	add expanded branches to existing data
		tree->AddFriend( indexTree );	//	add index branch to existing data
		tree->SetBranchAddress( "time", &timeb );	//	set address to store time
		nEntries = tree->GetEntries();	//	initialise number of entries
		std::cout << "Number of entries:	" << nEntries << std::endl;	//	inform user of number of entries
		helper = new Helper( nEntries );	//	required for countdown

		for( unsigned long long i = 0; i < nEntries; i++, jEntry = i ){	//	loop over all entries

			helper->countdown();	//	print progress
			tree->GetEntry( i );	//	grab entry

			if( ( timeb < 1000000000000 && tpreviousnoO  > 1000000000000 ) || timeb < timerunstart ){	//	look for clock reset

				if( i ){	//	we are past the first entry

					std::cout << std::endl << "Clock has been reset at time " << timeb << " (old run time " << timerunstart << ", old timeb " << tpreviousnoO << ") on entry number " << i << ", changing offset from " << tOffset << " to " << tOffset + tpreviousnoO - timeb << std::endl;	//	inform user and write out a new line to step out of counter
					tOffset += tpreviousnoO - timeb;	//	get time offset for use in setting beam pulse time
					setTimeOffset( tOffset );	//	set time offset

				}	//	end offset update

				timerunstart = timeb;	//	reset t0

			}	//	end clock reset check

			if( ( timeb + tOffset ) < tprevious ){	//	if our current time is less than our previous time

				thist = timeb + tOffset;	//	change previous time
				timeb++;	//	increment current time to make comparison easier

				while( movedFlag || thist < ( timeb + tOffset ) )	//	while our time is earlier than before
					tree->GetEntry( --jEntry );	//	get previous entry

				eindex[jEntry].push_back( { thist, i - jEntry } );	//	add entry number to the entry proceeding in original list
				setMovedFlag( true );	//	set indicator that entry has MovedFlag

			}	//	end out of sync entry check
			else{	//	if entry is time ordered

				tprevious = timeb + tOffset;	//	store current event time for next event clock comparison
				setMovedFlag( false );	//	initialise out of sync flag

			}	//	set values up for next time check

			tpreviousnoO = timeb;	//	store time without offset
			addedTree->Fill();	//	fill new branches

		}	//	end loop over all entries

		std::cout << std::endl << "Analysis completed in " << TS.CpuTime() << "s (which is " << TS.RealTime() << "s in real time), starting index creation..." << std::endl;	//	profiler output
		TS.Reset();	//	Reset profile time
		helper->resetCountdown();	//	reset countdown to first entry
		TS.Start();	//	Start profiling

		for( unsigned long long i = 0; i < nEntries; i++ ){	//	add shifted entry key

			helper->countdown();	//	print progress
			movedEntry.clear();	//	reset std::vector
			movedEntry.push_back( i );	//	add this entry to itself

			if( eindex.find( i ) != eindex.end() ){	//	if i is in the map

				sort( eindex[i].begin(), eindex[i].end() );	//	sort std::vector from lowest to highest
				movedSize = eindex[i].size();	//	get size of std::vector to use
				movedCtr += movedSize;	//	add to running counter

				for( unsigned long long j = 0; j < movedSize; j++ )	//	get each pair in std::vector
					movedEntry.push_back( eindex[i][j].second + i );	//	add entry number to std::vector

			}

			indexTree->Fill();	//	fill new branches			

		}	//	end loop to add shift entry key

		eindex.clear();	//	free up memory with no delete method available
		std::cout << std::endl << "Added index for " << movedCtr << " moved entries after " << eindex.size() << " entries in " << TS.CpuTime() << "s (which is " << TS.RealTime() << "s in real time), now our tree should exist" << std::endl;	//	profiler output
		addedTree->Write( "", TObject::kOverwrite );	//	Write new tree to file
		indexTree->Write( "", TObject::kOverwrite );	//	Write index to file
		tree->Write( "", TObject::kOverwrite );	//	Write main tree to file
		file->Close();	//	close properly
		delete file;	//	clean up completely

	}	//	end file does not exist

	std::cout << "Tree exists, using " << filename << std::endl;	//	inform user we are using existing tree
	file = new TFile( filename, "READ" );	//	read output file
	tree = (TChain *)file->Get( "DataTree" );	//	load tree
	addedTree = tree->GetFriend( "newTree" );	//	load new tree
	indexTree = tree->GetFriend( "iTree" );	//	load new tree
	tree->SetBranchAddress( "time", &timeb );	//	set address to store time
	nEntries = tree->GetEntries();	//	initialise number of entries
	helper = new Helper( nEntries );	//	required for countdown
	std::cout << "Tree loaded in " << TS.CpuTime() << "s (which is " << TS.RealTime() << "s in real time)" << std::endl;	//	profiler output
	jEntry = 1;	//	start check from the second entry as first entry is often strange
	tree->GetEntry( jEntry );	//	reset entry position

	while( timeb < ( unsigned long long ) 1 )	//	make sure we have a valid entry for time
		tree->GetEntry( jEntry++ );	//	get next entry

	std::cout << "There are " << nEntries << " entries, with a minimum timestamp of " << timeb << " (" << timeb / ( unsigned long long )pow( 10, floor( log10( (double)timeb ) ) ) << "e" << floor( log10( (double)timeb ) ) << ") and a maximum timestamp of ";	//	output minimum timestamp
	tree->GetEntry( nEntries - 1 );	//	get penultimate entry as last entry can be strange oh and do not worry about the timeb conversion it is only an estimate for the user
	std::cout << timeb << " (" << 1 + timeb / ( unsigned long long )pow( 10, floor( log10( (double)timeb ) ) ) << "e" << floor( log10( (double)timeb ) ) << ")" << std::endl;	//	output maximum timestamp
	getHists();	//	create histograms

}	//	end default constructor



Channel *FNT::addChannel( char d, int n, double p /* = 0 */ ){	//	add channel( type, channel, pixel )

	Channel *x = new Channel( d, n, p );	//	create new channel
	channels.insert( { n, x } );	//	add a new channel to the map

	if( d == 'C' || d == 'M' ) neutrons.push_back( x );	//	add channel to neutrons
	else if( d == 'G' ) germaniums.push_back( x );	//	add channel to germaniums
	else if( d == 'L' ) labr.push_back( x );	//	add channel to LaBr3s

	return x;

}


bool FNT::addChannels(){	//	add channels

	char d, c;	//	detector type, calibration type
	double n;	//	calibration term
	std::ifstream s( getFileChannels() );	//	open file
	double p;	//	pixel number
	int i, t;	//	channel, time offset
	std::string x = "";	//	hold lines
	std::stringstream h;	//	hold remainer of line in stream
	std::vector<double> vc, ve;	//	hold calibration terms, efficiency terms

	s.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );	//	skip opening comment line

	while( getline( s, x ) ){	//	while we have data...

		vc.clear();	//	reset calibration term std::vector
		ve.clear();	//	reset efficiency term std::vector
		h.str( x );	//	fill holding stream
		h.clear();	//	reset stream pointer
		h >> i >> d >> p >> t >> c;	//	hold channel number, detector type, pixel number, time offset, and calibration type

		while( h >> std::ws && h.peek() != 101 && h >> n )	//	while we have text which is not white space or 'e'
			vc.push_back( n );	//	rebuild calibration terms

		h.get();	//	skip over the efficiency holder

		while( h >> n )	//	while we have text
			ve.push_back( n );	//	rebuild efficiency terms

		Channel *newc = addChannel( d, i, p );	//	add channel
		std::cout << "Adding channel " << i << std::endl;	//	inform user
		newc->addCalibration( t, c, vc, ve );	//	add calibration to channel with ( offset, type, terms, efficiency )

		if( d == 'T' ) beamChannel = i;	//	set brfp channel
		else if( d == 'R' ) thetChannel = i;	//	set thet channel
		else if( d == 'X' ) xposChannel = i;	//	set xpos channel

	}	//	end data read in

	if( beamChannel < 0 || thetChannel < 0 || xposChannel < 0 )	//	essential channels not set
		std::cout << "WARNING: Essential channels not added!  Beam = " << beamChannel << ", r = " << thetChannel << ", x = " << xposChannel << std::endl;	//	warn user

	return true;	//	success

}	//	end addChannels



bool FNT::addGates(){	//	add gates

	char c;	//	gate type
	std::ifstream s( getFileGates() );	//	open file
	unsigned long long low, high;	//	low, high values
	int i;	//	channel number
	std::string x = "";	//	hold lines

//	s.ignore(std::numeric_limits<std::streamsize>::max(), '\n');	//	skip opening comment line

	while( getline( s, x ) ){	//	while we have data...

		s >> i >> c >> low >> high;	//	fill channel, type, low, and high

		if( getChannel( i ) ){	//	check for channel

			getChannel( i )->addGate( c, low, high );	//	add gate to channel with ( type, low, high )
			std::cout << ( c == 'e' ? "Energy" : ( c == 't' ? "Time" : "Coincidence" ) ) << " gate applied to channel " << i << " with minimum " << low << " and high " << high << std::endl;	//	inform user

		} else std::cout << "Gate for channel " << i << " exists, but channel " << i << " doesn't exist!" << std::endl;	//	no channel to add gate to

	}	//	end data read in

	return true;	//	success

}	//	end addGates



TChain *FNT::chainer(){	//	create new trees

	TChain *c = new TChain( "DataTree", "DataTree" );	//	create chain
	std::ifstream s( getFileFiles() );	//	open file
	char x[4096];	//	hold lines

	s.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );	//	skip opening comment line

	while( s.getline( x, sizeof( x ), '\n' ) ){	//	while we have data...

		if( x[strlen( x ) - 1] == '\r' )	x[strlen( x ) - 1] = '\0';	//	remove carriage returns

		std::cout << "Adding " << x << " to tree..." << std::endl;	//	inform user of file
		c->Add( x );	//	add file to our chain

	}	//	end data read in

	return c;	//	return pointer to chain

}	//	end chainer



bool FNT::getHists(){	//	get histograms

	newHists = new TFile( fileRootH, "RECREATE" );	//	create histogram file
	std::ifstream s( getFileHistos() );	//	open file
	int xbins = -1, ybins = -1, zbins = -1;	//	histogram bin counts
	double t, xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0;	//	type, histogram axis minmax
	std::string n = "", x = "";	//	hold name, lines

	s.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );	//	skip opening comment line

	while( getline( s, x ) ){	//	while we have data...

		s >> t >> n >> xbins >> xmin >> xmax;	//	fill values

		if( t < 1 ){	//	if 1-d histograms

			if( std::signbit( t ) ){	//	if 1-d histograms

				if( t < 0 ){	//	if multiple histograms

					folders.push_back( n );	//	use prefix to add to folder names

				}	//	end multiple histogram check

				t = std::abs( t );	//	convert the number of histograms

				for( int i = 0; i <= ( t < getMaxChannels() ? t : getMaxChannels() ); i++ ){	//	loop over all histograms

					histo( -0.0, n + std::to_string( i ), xbins, xmin, xmax );	//	create 1-d histogram

				}	//	end loop over histograms

			}	//	end 1-d histogram check
			else{	// else must be a 3-d histogram

				s >> ybins >> ymin >> ymax >> zbins >> zmin >> zmax;	//	fill values for 3-d
				histo( 0, n + "0", xbins, xmin, xmax, ybins, ymin, ymax, zbins, zmin, zmax );	//	create 3-d histogram

			}	//	end 3-d histogram check

		}	//	end 1-d histogram check
		else if( t > 0 ){	// if a 2-d histogram

			s >> ybins >> ymin >> ymax;	//	fill remaining values

			if( t > 1 ){	//	if multiple histograms

				folders.push_back( n );	//	use prefix to add to folder names

			}	//	end multiple histogram check

			for( int i = 0; i <= ( t < getMaxChannels() ? ( 1 == (int)t ? 0 : t ) : getMaxChannels() ); i++ ){	//	loop over all histograms

				histo( 1, n + std::to_string( i ), xbins, xmin, xmax, ybins, ymin, ymax );	//	create 2-d histogram

			}	//	end 2-d histogram loop

		}	//	end 2-d histogram check

		std::cout << "Added histogram collection with title:\t" << n << std::endl;	//	tell user where we are

	}	//	end data read in

	newHists->Write( "", TObject::kOverwrite );	//	write histograms to file
	return true;	//	success

}	//	end getHists()



bool FNT::histo( Float_t t, std::string n, int xbins, double xmin, double xmax, int ybins, double ymin, double ymax, int zbins, double zmin, double zmax ){	//	add histogram

	std::string s = helper->sanitiser( n );	//	copy string for use for name

	if( xbins < 1 ){	//	if no xbins have been set

		xbins = (int)std::abs( xmin + xmax ) / 2;	//	set bin count to one per channel

	}	//	end xbins fix

	if( t < 1 ){	//	if a 1-d histogram

		if( zbins < 1 ){	//	if no zbins have been set

			zbins = (int)std::abs( zmin + zmax ) / 2;	//	set bin count to one per channel

		}	//	end zbins fix

		if( std::signbit( t ) ){	//	check for signed zero

			h1s.insert( { s, new TH1F( s.c_str(), n.c_str(), xbins, xmin, xmax ) } );	//	fill 1-d histogram

		}	//	check for 1-d histogram
		else{	// else must be a 3-d histogram

			h3s.insert( { s, new TH3F( s.c_str(), n.c_str(), xbins, xmin, xmax, ybins, ymin, ymax, zbins, zmin, zmax ) } );	//	fill 3-d histogram

		}	//	end check of signed zero

	} else if( t > 0 ){	// if a 2-d histogram

		if( ybins < 1 ){	//	if no ybins have been set

			ybins = (int)std::abs( ymin + ymax ) / 2;	//	set bin count to one per channel

		}	//	end ybins fix

		h2s.insert( { s, new TH2F( s.c_str(), n.c_str(), xbins, xmin, xmax, ybins, ymin, ymax ) } );	//	add 2-d histogram to std::vector

	}	//	end histogram filling

	return true;	//	success

}	//	end adding histogram


}	//	end namespace fnt


void FNT(){}	//	allow root compilation