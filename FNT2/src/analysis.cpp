//
//	George O'Neill @ University of York, 2020
//
//	This file is by design a simple separator for analysis
//
#include "../include/analysis.h"
//#define TIMEFULLCAL true	//	uncomment to just generate 3-D time map against all pixels

namespace fnt{	//	create unique working area

void analysis::doAnalysis( FNT *f ){	//	start analysis

	Bool_t addOnce;	//	flag for coincidence gate for gammas
	Channel *c;	//	allocate channel pointer
	char dtypek;	//	detector type
	std::vector<std::pair<int, std::pair<Channel *, double>>> coincidences;	//	hold events occurring in same time pulse
	double nrjCal, npixel, w, fluxNorm = 0;	//	calibrated energy, pixel number, weight, neutron flux normalisation constant
	int timed, timeSinceB, npixelcj, npixelck, gateNumber, gateSpecif;	//	time difference, calibrated time since beam, x, r, coincident pixel number j, coincident pixel number k, gate number, specific gate type
	std::string s, bgHist;	//	string for names, background or normal
	long unsigned int coincSize = 0;	//	coincidence size
	unsigned int k = 0, neighbour = 0, lFlux = 0, rFlux = 0;	//	counter, neighbours, neutron counters
	std::unordered_map<int, bool> gammaFlag;	//	flag if gamma has been used
	std::map<std::pair<bool, int>, bool> gateFlag;	//	flag if gate has been used
	unsigned long long movedEntrySize, timeLastB = 0, timeLastBp = 0, timebo, tFluxIn = 0, tFluxT = 0;	//	entry size, time of last beam pulse, previous beam pulse time, offsetted time, flux timer start, flux timer total
	f->helper->resetCountdown();	//	reset countdown in case it has been used prior

	for( unsigned long long i = 0; i < n; i++ ){	//	loop over all entries

		f->helper->countdown();	//	print progress
		bigTree->GetEntry( i );	//	grab entry

		if( !movedFlag ){	//	construct events based on those in the correct position

			movedEntrySize = movedEntry->size();	//	store size for loop

			for( unsigned long long jEntry = 0; jEntry < movedEntrySize; jEntry++ ){	//	get entry by map

				bigTree->GetEntry( ( *movedEntry )[jEntry] );	//	grab ordered entry
				il = (int)label;	//	convert label to integer
				c = f->getChannel( il );	//	set channel
				f->getH1( "HitPatternAll0" )->Fill( il );	//	hit pattern

				if( c ){	//	skip unknown channels

					timebo = timeb + timeOffset;	//	offsetted time
					timed = (int)( timebo - timeLastB );	//	adjust timeb to have most of modulo removed

					if( il == chanB )	//	if we have a beam channel
						timeLastB = timebo;	//	update beam time
					else if( timed > bpt && timeLastB ){	//	else if we are due for a beam channel

						timeLastB = timebo - timed % bpt;	//	insert a fake beam time
						timed = (int)( timebo - timeLastB );	//	readjust timeb to have most of modulo removed

					}	//	end simulated beam pulse

					timeSinceB = ( bpt + timed + c->getTOffset() ) % bpt;	//	store aligned time shifted to always lie between 0 and bpt
					s = std::to_string( label );	//	get channel number stored as string instead of char array
					npixel = c->getPixelNumber();	//	store pixel
					nrjCal = c->adjE( nrj );	//	calibrate nrj
					f->getH1( "nrjRawC" + s )->Fill( nrj );	//	add to uncalibrated energy histogram
					f->getH1( "nrj2RawC" + s )->Fill( nrj2 );	//	add to slow uncalibrated energy histogram
					f->getH1( "timeRawC" + s )->Fill( timed );	//	add to uncalibrated time histogram
					f->getH1( "rateRawC" + s )->Fill( (double)timebo );	//	add to rate histogram
#ifndef NOBEAM
					f->getH3( "nrjCalVtimeCal0" )->Fill( nrjCal, timeSinceB, npixel );	//	corrected energy V corrected time

					if( cxFlag ){	//	if we are in useful parts of the table

						if( nrjCal > 0 && crFlag && timeLastB ){	//	if a real beam pulse has been seen

							if( (int)npixel ){	//	detector gate
#ifdef TIMEFULLCAL
								if( !tarr ) unsigned long long tarr[65];	//	initialise variable
								tarr[(int)npixel] = timeb + timeOffset;	//	add time of hit to array

								for( int k = 1; k < 65; k++ )	//	for each pixel
									if( tarr[k] > 0 )	//	if we have registered a hit previously
										f->getH3( "neutronTimeByPixel0" )->Fill( k, npixel, timeb + timeOffset - tarr[k] );	//	add to projection histogram
#else
								if( c->passes( timeSinceB, 't' ) && c->passes( nrjCal, 'e' ) ){	//	if passing channel gate

									f->getH3( "nrjGtdVtimeGtd0" )->Fill( nrjCal, timeSinceB, npixel );	//	gated corrected energy V gated corrected time
									w = fluxNorm;	//	set weight to one to add real count
									bgHist = "";	//	ensure normal histograms are filled

								}	//	end inside check
								else{	//	background event

									w = -1 * fluxNorm * c->gatesTcover( bpt );	//	set weight for background subtraction
									bgHist = "Bg";	//	fill background histograms

								}	//	end gate logic

								if( npixel > 0 ){	//	neutron gate alternatively il>10&&il<90&&(il%10>0&&il%10<9)

									f->getH2( "pixels" + bgHist + "0" )->Fill( ( (int)npixel - 1 ) % 8, floor( ( npixel - 1 ) / 8. ), w );	//	fill hit pattern
									f->getH1( "HitPatternNeutronGtd" + bgHist + "0" )->Fill( npixel, w );	//	neutrons hit pattern
									f->getH1( "FluxRneutron" + bgHist + "0" )->Fill( rposition, w );	//	intensity of neutrons at r
									f->getH1( "RealFluxRneutron" + bgHist + "0" )->Fill( rposition );	//	real intensity of neutrons at r
									f->getH1( "FluxXneutron" + bgHist + "0" )->Fill( xposition, w );	//	intensity of neutrons at x
									f->getH2( "FluxXRneutron" + bgHist + "0" )->Fill( xposition, rposition, w );	//	intensity of neutrons over positions
									objectImg( f->getH2( "Objectneutron" + bgHist + "0" ), w );	//	make object

								}	//	end neutron check
								else{	//	must be a germanium channel

									f->getH1( "FluxRgamma" + bgHist + "0" )->Fill( rposition, w );	//	intensity of gammas at r
									f->getH1( "FluxXgamma" + bgHist + "0" )->Fill( xposition, w );	//	intensity of gammas at x
									f->getH2( "FluxXRgamma" + bgHist + "0" )->Fill( xposition, rposition, w );	//	intensity of gammas over positions
									f->getH3( "XRgamma" + bgHist + "0" )->Fill( xposition, rposition, nrjCal, w );	//	XY-plane against energy
									objectImg( nrjCal, w );	//	Fill object picture
									gateSpecif = c->passes( nrjCal, 'G', 0 );	//	get the gate pass as specific gate

									if( gateSpecif )	//	if we have an interesting gamma
										objectImg( f->getH2( "Objectgamma" + bgHist + "G" + std::to_string( gateSpecif ) ), w );	//	make object

								}	//	end germanium fills

								if( timeLastB == timeLastBp ){	//	if this entry belongs to current beam pulse

									if( timeLastB > 0 )	//	if the entry is after the first beam pulse
										coincidences.push_back( { timeSinceB, {c, nrjCal} } );	//	add time since beam, channel, and energy to std::vector

								}	//	end check if beam pulse is the same
								else{	//	a new beam pulse

									sort( coincidences.begin(), coincidences.end() );	//	sort coincidences into time order
									coincSize = coincidences.size();	//	store coincidence size

	//								if( timeLastB < timeLastBp )	std::cout << std::endl << "Warning: Your coincidence std::vector is " << coincSize << " and possibly your beam pulses are out of order at " << i << " where timeLastB = " << timeLastB << " and previous timeLastB = " << timeLastBp << std::endl;	//	order check
									for( unsigned int j = 0; j < coincSize; j++ ){	//	for each coincidence

										k = j + 1;	//	set k to next item at the start of the loop
										neighbour = 0;	//	reset neighbour count
										addOnce = true;	//	set flag for adding to germanium one time only
										npixelcj = (int)coincidences[j].second.first->getPixelNumber();	//	get pixel number for first pixel

										while( k < coincSize ){	//	check k is within std::vector limits and within time constant

											npixelck = (int)coincidences[k].second.first->getPixelNumber();	//	get pixel number for second pixel
											dtypek = coincidences[k].second.first->getType();	//	fetch detector type for coincident pixel
											gateNumber = coincidences[j].second.first->passes( coincidences[k].first, dtypek, 1 );	//	get the gate pass as gate
											gateSpecif = coincidences[j].second.first->passes( coincidences[k].first, dtypek, 0 );	//	get the gate pass as specific gate
											f->getH3( "coincTimeMap0" )->Fill( npixelcj, npixelck, coincidences[k].first - coincidences[j].first );	//	fill coincidence pixel map for each time
											f->getH3( "pixelsfiring0" )->Fill( npixelcj, npixelck, gateNumber );	//	create maps of firing pixels in coincidence

											if( gateNumber ){	//	only proceed if we passed a gate

												if( npixelck > 0 ){	//	if it is a neutron channel do further checks

													if( !gateFlag[{std::signbit( npixelcj ), gateNumber}] /*&& ( (npixelcj > 0 && gateSpecif == gateNumber) || (npixelcj < 0 && gateSpecif) )*/ ){	//	if we have a clean event here coupled with a suitable parent event

														gateFlag[{std::signbit( npixelcj ), gateNumber}] = f->getH2( "nrjCoiVtimeCoiGate" + std::to_string( gateNumber ) + "C" + std::to_string( coincidences[j].second.first->getChannelNumber() ) )->Fill( coincidences[j].second.second, coincidences[j].first );	//	fill neutron if followed by neutron or gamma if followed by neutron and do not fire on this gate again

														if( npixelcj > 0 )	//	if we have a neutron
															neighbour += f->helper->neighbour( npixelcj, npixelck );	//	increment neighbour if there is one
														else if( addOnce && gateSpecif )	//	 if we have not added a gamma yet
															addOnce = 0 && f->getH2( "nrjCoiVtimeCoiGateG0" )->Fill( coincidences[j].second.second, coincidences[j].first );	//	add to a map which will not repeat

													}	//	end neutron checks

												}	//	end neutron channel check
												else if( !gammaFlag[k] && npixelcj > 0 && gateSpecif )	//	result is a neutron followed by gamma
													gammaFlag[k] = f->getH2( "nrjCoiVtimeCoiGate" + std::to_string( gateSpecif ) + "C" + std::to_string( std::abs( npixelck ) ) )->Fill( coincidences[k].second.second, coincidences[k].first );	//	fill gamma following neutron and do not fire on this flag again

											}	//	end gate logic

											k++;	//	increment k for next

										}	//	end coincidence time check

										if( npixelcj > 0 ){	//	check if we have a neutron pixel which is quicker than the methods below filling below zero

											if( neighbour )	//	if we see a neighbour
												objectImg( f->getH2( "ObjectneutronCoi0" ) );	//	make neutron object from this coincident guy

											if( neighbour < 9 )	//	if a 'real' number of neighbours was seen
												f->getH2( "pixelsNeighbours" + std::to_string( neighbour ) )->Fill( ( npixelcj - 1 ) % 8, floor( ( npixelcj - 1 ) / 8. ) );	//	add to appropriate hit patten
											else	//	fill overflow plot
												f->getH2( "pixelsNeighbours9" )->Fill( ( npixelcj - 1 ) % 8, floor( ( npixelcj - 1 ) / 8. ) );	//	overflow hit pattern

										}	//	end neutron pixel check

									}	//	end coincidence check for loop

									coincidences.clear();	//	reset coincidences std::vector
									gateFlag.clear();	//	reset flags for gates
									gammaFlag.clear();	//	reset flags for gammas

								}	//	end check of new beam pulse

								timeLastBp = timeLastB;	//	store this beam pulse time for a valid neutron for comparisons

							}	//	end pixel checks
							else if( std::signbit( npixel ) ){	//	LaBr3 as -0 fails first check

								f->getH1( "FluxRlabr0" )->Fill( rposition );	//	intensity of gammas at r
								f->getH1( "FluxXlabr0" )->Fill( xposition );	//	intensity of gammas at x
								f->getH2( "FluxXRlabr0" )->Fill( xposition, rposition );	//	intensity of gammas over positions
								f->getH3( "XRlabr0" )->Fill( xposition, rposition, nrjCal );	//	XY-plane against energy
#endif
							}	//	end doing stuff with LaBr3

						}	//	detector operations

					}	//	within table
					else{	//	outside table

						if( npixel > 0 && ( c->passes( timeSinceB, 't' ) && c->passes( nrjCal, 'e' ) ) ){	//	if we have a neutron

							if( ( xcal > 0 && lFlux != 0 ) || ( xcal < 0 && rFlux != 0 ) )	//	if the other side has a flux then edge is now entered
								tFluxIn = timebo;	//	store entry time to flux collector
							else tFluxT = timebo - tFluxIn;	//	otherwise store time in case we leave the flux collector

							if( xcal > 0 ){	//	if on the positive side of the table
								lFlux = 0;	//	reset the other side of the table flux factor to zero
								rFlux++;	//	increment the flux for the positive side
							} else	//	reset the positive side of flux to sero
							{	//	on the negative side of the table
								rFlux = 0;	//	set other side of the table flux factor to zero
								lFlux++;	//	increase the negative side
							}	//	end side check

							fluxNorm = (double)tFluxT / ( ( lFlux + rFlux ) * 10000000000. );	//	adjust flux adjuster by adjusting collected time into seconds and then dividing by 100 but taking the inverse

						}	//	end pixel check

					}	//	end table flag

					if( il == chanX )	//	if we have an x channel
						setX( nrj );	//	set x position
					else if( il == chanR )	//	if we have a radius channel
						setR( nrj );	//	set r position
#endif
				}	//	end valid channel check

			}	//	end reordered index loop

		}	//	check if entry matches real index

	}	//	end loop over all entries

}	//	end doAnalysis( FNT* f )


void analysis::histogramOperations( fnt::FNT *f ){	//	start histogram alterations

	TDirectory *d;	//	hold value for directories
	TH1F *h;	//	hold value for histogram to move
	std::string s;	//	string for histograms

	for( std::string folder : f->getFolders() ){	//	for each folder

		std::cout << "Moving histograms to hists_" << folder << " folder..." << std::endl;	//	inform user
		d = newHists->mkdir( ( "hists_" + folder ).c_str() );	//	make a folder

		for( unsigned int i = 0; i <= gmc; i++ ){	//	for each channel number

			s = f->helper->sanitiser( folder + std::to_string( i ) );	//	create name
			h = (TH1F *)newHists->FindObjectAny( s.c_str() );	//	histogram to move

			if( h ){	//	if histogram exists

				if( ( folder == "timeCalC" || folder == "nrjCalC" ) && h->GetEntries() > 0 ){	//	if time or energy spectrum

					f->helper->peakf( h, folder );	//	get peak positions

				}	//	end time spectrum specialisations

				h->SetDirectory( d );	//	move histogram
				newHists->Delete( ( s + ";1" ).c_str() );	//	delete from original folder

			}	//	end move of histogram

		}	//	end for loop over each histogram

	}	//	end for loop over each folder

}	//	end histogramOperations( FNT* f )


void analysis::objectImg( TH2F *h, double w /* = 1 */ ){	//	image maker

	for( double yd = -0.; yd < tableSizeDiff; yd += objStep ){	//	for a small step in our circle

		h->Fill( XcosR + ( yd * sinrpos ), XsinR - ( yd * cosrpos ), w );	//	object histogram
		yd = -1 * yd - !std::signbit( yd ) * objStep;	//	repeat for the negative value of yd

	}	//	end object maker loop

}	//	end objectImg( TH2F* h, double w /* = 1 */ )

void analysis::objectImg( double e, double w /* = 1 */ ){	//	image maker

	for( double yd = -0.; yd < tableSizeDiff; yd += objStep ){	//	for a small step in our circle

		gammaBrowser->Fill( XcosR + ( yd * sinrpos ), XsinR - ( yd * cosrpos ), e, w );	//	object histogram
		yd = -1 * yd - !std::signbit( yd ) * objStep;	//	repeat for the negative value of yd

	}	//	end object maker loop

}	//	end objectImg( double e, double w /* = 1 */ )

void analysis::setR( int np ){	//	set r positions( nrj )

	rposition = np;	//	set rposition variable
	rcal = ( np - rmin ) / rdiv;	//	calibrate r
	sinrpos = sin( rcal );	//	set sin of rcal
	cosrpos = cos( rcal );	//	set cos of rcal
	setTrig();	//	set numbers used in functions
	crFlag = cR->passes( rposition, 'e' );	//	set gate flag of position

}	//	set r positions


void analysis::setX( int np ){	//	set x positions( nrj )

	xposition = np;	//	set xposition variable
	xcal = ( ( np - xmin ) / xdiv ) - xmoving;	//	calibrate x
	tableSizeDiff = sqrt( tableSize - xcal * xcal );	//	set loop maximum
	setTrig();	//	set numbers used in functions
	cxFlag = ( cX->passes( xposition, 'e' ) == 1 );	//	set gate flag of position

}	//	set x positions( int )


void analysis::setTrig(){	//	set sin and cos

	XsinR = xcal * sinrpos;	//	for use in relevant loops
	XcosR = xcal * cosrpos;	//	for use in relevant loops

}	//	end setting of trigonometric variables

/*void analysis::histogramPretty( TTree *bigTree ){	//	output a colourful hit pattern

	std::cout << "Starting colourful hit pattern now...it's not efficient so it will take some time, I think it is with nlogn not n!" << std::endl;
	int bins[3] = { 256, 0, 255 };	//	histogram {bin count, X minimum, X maximum}
	int ic[7] = { kGray, kBlue, kViolet, kTeal, kRed, kPink, kGreen };	//	histogram colours
	TCanvas *cA = new TCanvas();	//	drawing canvas for histogram
	TH1F *hist_arr7[7];	//	histogram array
	TLegend *l = new TLegend( 0.7, 0.7, 1, 1 );	//	initialise legend
	bool tf[256] = { false };	//	initialise an array for all channels
	string sh[7] = { "All Channels", "Neutrons", "Germaniums", "LaBr3", "Nothing", "RF", "Table" };	//	histogram titles
	string sc[7] = { "", "((label%10>0&&label%10<9)&&label>10&&label<90)||label==2", "label==90||label==91", "label==3",  "label==0||label==1||label==4||label==7||label==8||label==10||label==19||label==20||label==29||label==30||label==39||label==40||label==49||label==50||label==59||label==60||label==234||label==235||label==243||label==244||label==245||label==246||label==247||label==248||label==249||label==250||label==253||label==254||label==255", "label==97", "label==95||label==96" };	//	histogram conditions GO

	for( unsigned int i = 0; i < ( sizeof( hist_arr7 ) / sizeof( TH1F * ) ); i++ ){	//	loop over histogram array

		hist_arr7[i] = new TH1F( ( "hist_arr[" + std::to_string( i ) + "]" ).c_str(), sh[i].c_str(), bins[0], bins[1], bins[2] );	//	create new histogram
		hist_arr7[i]->SetLineColor( ic[i] );	//	set histogram colour
		bigTree->Draw( ( "label>>hist_arr[" + std::to_string( i ) + "]" ).c_str(), sc[i].c_str() );	//	grab data from tree and add to histogram
		l->AddEntry( hist_arr7[i], sh[i].c_str() );	//	add histogram entry to legend with title

	}	//	end loop over histogram array

	cA->SetLogy();	//	logarithmic scale
	std::cout << "Colourful hit pattern is done!" << std::endl;	//	inform user

}	//	end colourful hit pattern
*/

}	//	end namespace fnt


void analysis(){}	//	allow root compilation