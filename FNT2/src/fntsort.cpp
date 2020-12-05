//
//	George O'Neill @ University of York, 2020/02/05
//
//	This is a simple program for adapting ROOT files created by 'ungroup2tree' at Orsay for FNT2, and:
//	-	Adds new branches to help time order
//	-	Stores a chained root file and checks whether it exists so should not be rewritten
//	-	Sets up channels, calibrations, gates, and histograms in FNT.h
//	-	Individually adds calibrations, gates, and adjustments in channels.h
//	-	Provides a simple pass check of all numeric conditions in gates.h
//	-	Allows analysis to be fully contained within analysis.cpp
//

//	FNT2 libraries
#include "../include/analysis.h"

int main(){	//	main method will run when compiled

//	fnt::FNT* f = runfnt();	//	method that actually does stuff
	new fnt::analysis( new fnt::FNT() );	//	perform analysis
	return 1;	//	exit gracefully

}	//	end main method for compilation


void fntsort(){}	//	allow root compilation