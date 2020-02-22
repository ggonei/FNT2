//
//	George O'Neill @ University of York, 2020/02/05
//
//	This is a simple program for adapting ROOT files created by 'ungroup2tree' at Orsay for FNT2, and:
//	-	Adds a new branch containing the time of the last beam pulse, last known x, and last known theta for each event
//	-	Sets up channels, calibrations, gates, and histograms in analyser.h
//	-	Individually adds calibrations, gates, and adjustments in channels.h
//	-	Provides a simple pass check of all numeric conditions in gates.h
//	-	Allows analysis to be fully contained within analysis.cpp
//	-	Stores a chained root file and checks whether it exists so should not be rewritten
//

//	FNT2 libraries
#include "include/analysis.h"

void fntsort() {	//	main method will run when compiled

//	fnt::FNT* f = runfnt();	//	method that actually does stuff
	new fnt::analysis(new fnt::FNT());	//	perform analysis

}	//	end main method for compilation