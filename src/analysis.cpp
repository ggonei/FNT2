//
//	George O'Neill @ University of York, 2020
//
//	This file is by design a simple separator for analysis
//
#include "../include/analysis.h"

namespace fnt {	//	create unique working area

void analysis::do_analysis(fnt::FNT* f) {	//	start analysis

	Channel* c;	//	allocate channel pointer
	std::vector<std::pair<Long64_t,std::array<UInt_t,3>>> coincidences;	//	hold events occurring in same time pulse
	Double_t flux = 0;	//	flux
	Int_t npixel, timeSinceB;	//	get pixel, time since beam, coincidence size
	std::string s;	//	string for names
	UInt_t coincSize = 0, k = 0, neighbour = 0, xPrev = 20001;	//	counters, neighbours, previous x position
	ULong64_t timeLastBp = 0, timein = 0, tarr[65];	//	previous beam pulse time, time in to beam pulse, coincidence array
	f->helper->resetCountdown();	//	reset countdown in case it has been used prior

	for( ULong64_t i = 0; i < n; i++ ) {	//	loop over all entries

		f->helper->countdown();	//	print progress
		bigTree->GetEntry(i);	//	grab energy info
		il = (Int_t) label;	//	convert label to integer
		c = f->getChannel(il);	//	set channel
		f->getH1("HitPattern0")->Fill(il);	//	hit pattern

		if( c ) {	//	skip unknown channels

			timeb = timeb + timeOffset;	//	adjust timeb to include offset
			s = to_string(il);	//	get channel number stored as string
			timeSinceB = (timeb-timeLastB-c->getTOffset())%bpt;	//	store aligned time
			f->getH1("timeC" + s)->Fill((timeb-timeLastB)%bpt);	//	add to time histogram
			npixel = c->getPixelNumber();	//	store pixel

			if( npixel > 0 )	{	//	neutron gate alternatively il>10&&il<90&&(il%10>0&&il%10<9)

				f->getH2("timeVneutrons0")->Fill((timeb-timeLastB)%bpt, npixel);	//	neutron channels V time
				f->getH2("timeadjVneutrons0")->Fill(timeSinceB, npixel);	//	time adjusted V channel

				if( c->passes(timeSinceB, 't') )	{	//	if passing channel gate

					if( xPos <= 20000 || xPos >= 202000 )	{	//	if we are at the edge of the table

						if( xPrev > 20000 && xPrev < 202000 )	{	//	if the last entry was on the table

							flux = 0;	//	reset flux
							timein = timeb;	//	set time into 

						}	//	end edge entry check

						flux++;	//	increment flux

					}	//	end table edge check
					else if( xPrev <= 20000 || xPrev >= 202000 )	{	//	else if the last entry was at the edge of the table

						flux = flux / ((timeb - timein) / 1e12);	//	edit flux from raw counts to per second
						f->getH1("NeutronFlux0")->Fill(timeb, flux);	//	plot flux

					}	//	end table edge check

					if( flux > 0 )	{	//	if flux is observed

						f->getH2("NeutronRatePerPixel0")->Fill(timeb, npixel);	//	add hit to rate for pixel

					}	//	end flux observed check

					f->getH2("pixels0")->Fill((npixel-1)%8, floor((npixel-1)/8));	//	fill hit pattern
					f->getH1("HitPatternNeutrons0")->Fill(il);	//	neutrons hit pattern
					f->getH2("neutronfluxatrx0")->Fill(xPos, rPos, 1/(1+flux));	//	intensity of neutrons over positions
					f->getH1("neutronfluxatx0")->Fill(xPos);	//	intensity of neutrons at x
					f->getH1("neutronfluxatr0")->Fill(rPos);	//	intensity of neutrons at r
					f->getH1("nrjadjC" + s)->Fill(c->adjE(nrj));	//	add to nrj histogram
					f->getH1("nrj2adjC" + s)->Fill(c->adjE(nrj2));	//	add to nrj2 histogram
					f->getH2("nrjadjVchan0")->Fill(c->adjE(nrj), il);	//	energy V channel
					f->getH2("nrj2adjVchan0")->Fill(c->adjE(nrj2), il);	//	nrj2 V channel
					xPrev = xPos;	//	store previous x position

					if( timeLastB == timeLastBp )	{	//	if this entry belongs to current beam pulse

						if( timeLastB > 0 )	{	//	if the entry is after the first beam pulse
							
							coincidences.push_back(std::make_pair((Long64_t) timeSinceB, std::array<UInt_t,3>{(UInt_t)npixel, xPos, rPos}));	//	add time, pixel and position to vector

						}	//	end first beam pulse check

					}	//	end check if beam pulse is the same
					else	{	//	a new beam pulse

						sort(coincidences.begin(), coincidences.end());	//	sort coincidences into time order
						coincSize = coincidences.size();	//	store coincidence size

						for( UInt_t j = 0; j < coincSize; j++ )	{	//	for each coincidence

							k = j + 1;	//	set k to next item start loop
							neighbour = 0;	//	reset neighbour count

							while( k < coincSize && coincidences[j].first > coincidences[k].first - coincWind )	{	//	check k is within vector limits and within time constant

								f->getH3("coincTimeMap0")->Fill( coincidences[j].second[0], coincidences[k].second[0], coincidences[k].first - coincidences[j].first );	//	fill coincidence pixel map for each time
								neighbour += f->helper->neighbour(coincidences[j].second[0], coincidences[k].second[0]);	//	increment neighbour if there is one
								k++;	//	increment k for next

							}	//	end coincidence time check

							f->getH2("pixelsCleaned" + to_string(neighbour))->Fill((coincidences[j].second[0]-1)%8, floor((coincidences[j].second[0]-1)/8));	//	fill hit pattern

						}

						coincidences.clear();	//	reset coincidences vector

					}	//	end check of new beam pulse

					timeLastBp = timeLastB;	//	store this beam pulse time for a valid neutron for comparisons

				}	//	end gate check

			}	//	end neutron match


			if(	il == chanX )	{	//	do stuff with the x position
				
				f->getH2("xVtime0")->Fill(timeb, nrj);	//	x position V time
				
			}	//	end doing stuff with x

			if(	il == chanR )	{	//	do stuff with the theta position
				
				f->getH2("rVtime0")->Fill(timeb, nrj);	//	theta position V time
				
			}	//	end doing stuff with theta

			f->getH1("nrjC" + s)->Fill(nrj);	//	add to nrj histogram
			f->getH1("nrj2C" + s)->Fill(nrj2);	//	add to nrj2 histogram
			f->getH2("nrjVchan0")->Fill(nrj, il);	//	nrj V time
			f->getH2("nrj2Vchan0")->Fill(nrj2, il);	//	nrj2 V time
			f->getH2("timeVchan0")->Fill((timeb-timeLastB)%bpt, il);	//	channel V time
			f->getH1("timeadjC" + s)->Fill(timeSinceB);	//	add to adjusted time histogram
			f->getH2("timeadjVchan0")->Fill(timeSinceB, il);	//	channel V channel
			f->getH2("nrjVtime" + s)->Fill(timeSinceB, nrj);	//	energy V time
			f->getH2("nrj2Vtime" + s)->Fill(timeSinceB, nrj2);	//	nrj2 V time

		}	//	end valid channel check

	}	//	end loop over all entries

}	//	end do_analysis( FNT* f )


void analysis::histogram_operations(fnt::FNT* f) {	//	start histogram alterations

	TDirectory* d;	//	hold value for directories
	TH1D* h;	//	hold value for histogram to move
	std::string s;	//	string for histograms

	for( string folder : f->getFolders() )	{	//	for each folder

		std::cout << "Moving histograms to hists_" << folder << " folder..." << std::endl;	//	inform user
		d = newHists->mkdir(("hists_" + folder).c_str());	//	make a folder

		for( UInt_t i = 0; i <= gmc; i++ ) {	//	for each channel number

			s = f->helper->sanitiser(folder + to_string(i));	//	create name
			h = (TH1D*)newHists->FindObjectAny(s.c_str());	//	histogram to move

			if( h ) {	//	if histogram exists

				if( folder == "timeC" || folder == "nrjC" )	{	//	if time or energy spectrum

					f->helper->peakf(h, folder);	//	get peak positions

				}	//	end time spectrum specialisations

				h->SetDirectory(d);	//	move histogram
				newHists->Delete((s + ";1").c_str());	//	delete from original folder

			}	//	end move of histogram

		}	//	end for loop over each histogram

	}	//	end for loop over each folder

}	//	end histogram_operations( FNT* f )


void analysis::histogram_pretty(TTree* bigTree) {	//	output a colourful hit pattern

	std::cout << "Starting colourful hit pattern now...it's not efficient so it will take some time, I think it is with nlogn not n!" << std::endl;
	Int_t bins[3] = {256, 0, 255};	//	histogram {bin count, X minimum, X maximum}
	Int_t ic[7] = {kGray, kBlue, kViolet, kTeal, kRed, kPink, kGreen};	//	histogram colours
	TCanvas* cA = new TCanvas();	//	drawing canvas for histogram
	TH1D* hist_arr7[7];	//	histogram array
	TLegend* l = new TLegend(0.7,0.7,1,1);	//	initialise legend
	bool tf[256] = { false };	//	initialise an array for all channels
	string sh[7] = {"All Channels", "Neutrons", "Germaniums", "LaBr3", "Nothing", "RF", "Table"};	//	histogram titles
	string sc[7] = {"", "((label%10>0&&label%10<9)&&label>10&&label<90)||label==2", "label==90||label==91", "label==3",  "label==0||label==1||label==4||label==7||label==8||label==10||label==19||label==20||label==29||label==30||label==39||label==40||label==49||label==50||label==59||label==60||label==234||label==235||label==243||label==244||label==245||label==246||label==247||label==248||label==249||label==250||label==253||label==254||label==255", "label==97", "label==95||label==96"};	//	histogram conditions GO

	for (UInt_t i = 0; i < (sizeof(hist_arr7)/sizeof(TH1D*)); i++) {	//	loop over histogram array

		hist_arr7[i] = new TH1D(("hist_arr[" + std::to_string(i) + "]").c_str(), sh[i].c_str(), bins[0], bins[1], bins[2]);	//	create new histogram
		hist_arr7[i]->SetLineColor(ic[i]);	//	set histogram colour
		bigTree->Draw(("label>>hist_arr[" + std::to_string(i) + "]").c_str(), sc[i].c_str());	//	grab data from tree and add to histogram
		l->AddEntry(hist_arr7[i], sh[i].c_str());	//	add histogram entry to legend with title

	}	//	end loop over histogram array

	cA->SetLogy();	//	logarithmic scale
	std::cout << "Colourful hit pattern is done!" << std::endl;	//	inform user

}	//	end colourful hit pattern


}	//	end namespace fnt


void analysis(){};	//	allow root compilation