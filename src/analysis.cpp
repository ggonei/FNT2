//
//	George O'Neill @ University of York, 2020/02/18
//
//	This file is by design a simple separator for analysis
//
#include "../include/analysis.h"

namespace fnt {	//	create unique working area

void analysis::do_analysis(fnt::FNT* f) {	//	start analysis

	Channel* c;	//	allocate channel pointer
	deque<pair<Long64_t,pair<Int_t,pair<Int_t,Int_t>>>> coincCheck;	//	hold coincident pixels
	vector<pair<Long64_t,pair<Int_t,pair<Int_t,Int_t>>>> coincidences;	//	hold events occurring in same time pulse
	Double_t flux = 0;	//	flux
	Int_t npixel, timeSinceB;	//	get pixel, time since beam
	std::string s;	//	string for names
	UInt_t xPrev = 20001;	//	previous x position
	ULong64_t timeLastBp = 0, timein = 0, tarr[65];	//	previous beam pulse time, time in to beam pulse, coincidence array
	f->helper->resetCountdown();	//	reset countdown in case it has been used prior
 TH3D* h3d = new TH3D("h","h",64,0,64,64,0,64,500,0,500);
 
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
				f->getH2("NeutronRate0")->Fill(timeb, npixel);	//	add hit to rate

				if( c->passes(timeSinceB, 't') )	{	//	if passing channel gate

					if( xPos <= 20000 || xPos >= 202000 )	{	//	if we are at the edge of the table

						if( xPrev > 20000 && xPrev < 202000 )	{	//	if the last entry was on the table

							flux = 0;	//	reset flux
							timein = timeb;	//	set time into 

						}

						flux++;	//	add one to flux

					}	//	end table edge check
					else if( xPrev <= 20000 || xPrev >= 202000 )	{	//	else if the last entry was at the edge of the table

						flux = flux / ((timeb - timein) / 1e12);	//	edit flux from raw counts to per second
						f->getH1("NeutronFlux0")->Fill(timeb, flux);	//	plot flux

					}

					f->getH2("neutronfluxatrx0")->Fill(xPos,rPos,1/(1+flux));	//	intensity of neutrons over positions
					if( flux > 0 )	f->getH2("NeutronRateGated0")->Fill(timeb, npixel);	//	add hit to rate
					f->getH2("pixels0")->Fill((npixel-1)%8, floor((npixel-1)/8));	//	fill hit pattern
					f->getH1("HitPatternNeutrons0")->Fill(il);	//	neutrons hit pattern
					f->getH1("neutronfluxatx0")->Fill(xPos);	//	intensity of neutrons at x
					f->getH1("neutronfluxatr0")->Fill(rPos);	//	intensity of neutrons at r
					f->getH1("nrjC" + s)->Fill(nrj);	//	add to nrj histogram
					f->getH1("nrj2C" + s)->Fill(nrj2);	//	add to nrj2 histogram
					xPrev = xPos;	//	store previous x position

					if( timeLastB == timeLastBp )	{	//	if this entry belongs to current beam pulse

						if( timeLastB > 0 )	coincidences.push_back(make_pair((Long64_t) timeSinceB, make_pair(npixel,make_pair(xPos,rPos))));	//	add time and pixel to vector

					}	//	end check if beam pulse is the same
					else	{	//	a new beam pulse

						sort(coincidences.begin(),coincidences.end());	//	sort coincidences into time order

						for( UInt_t j = 0; j < coincidences.size(); j++ )	{	//	for each entry in this beam pulse

							tarr[coincidences[j].second.first] = coincidences[j].first;	//	store the time for this pixel in an array							
							coincCheck.push_back(coincidences[j]);	//	add this pixel to coincidence time check

							for( ULong64_t tpiece : {(coincidences[j].second.first + 8 <65?tarr[coincidences[j].second.first + 8]:coincidences[j].first), (coincidences[j].second.first - 8 >0?tarr[coincidences[j].second.first - 8]:coincidences[j].first), ( ( 8 - ( coincidences[j].second.first - 1 ) % 8 ) > 1 ? tarr[coincidences[j].second.first + 1] : coincidences[j].first ), ( ( ( coincidences[j].second.first - 1 ) % 8 ) > 0 ? tarr[coincidences[j].second.first - 1] : coincidences[j].first )} )	f->getH2("multiplicityt0")->Fill(coincidences[j].first - tpiece, coincidences[j].second.first);	//	fill time map

						}	//	end for loop over coincidences
						
						UInt_t j=0, k = 0, neighbour = 0;

						for( pair<Long64_t,pair<Int_t,pair<Int_t,Int_t>>> p : coincCheck )	{

							k = j;
							neighbour = 0;
							while( k < coincCheck.size() && p.first > coincCheck[k].first - 40 ){
		//						if( !(p.second == coincCheck[k].second && p.first == coincCheck[k].first) )	f->getH2("coincidences0")->Fill(p.second, coincCheck[k].second);	//	fill coincidence map
/*								if( p.first > coincCheck[k].first - 1000 && p.first < coincCheck[k].first + 1000 ){f->getH2("coincidences10000")->Fill(p.second, coincCheck[k].second);	f->getH1("pixelsfired10C" + to_string(coincCheck[k].second))->Fill( p.second );}
								else if( p.first > coincCheck[k].first - 2000 && p.first < coincCheck[k].first + 2000 ){f->getH2("coincidences20000")->Fill(p.second, coincCheck[k].second);	f->getH1("pixelsfired20C" + to_string(coincCheck[k].second))->Fill( p.second );}
								else if( p.first > coincCheck[k].first - 5000 && p.first < coincCheck[k].first + 5000 ){f->getH2("coincidences50000")->Fill(p.second, coincCheck[k].second);	f->getH1("pixelsfired50C" + to_string(coincCheck[k].second))->Fill( p.second );}
								else if( p.first > coincCheck[k].first - 10000 && p.first < coincCheck[k].first + 10000 ){f->getH2("coincidences100000")->Fill(p.second, coincCheck[k].second);	f->getH1("pixelsfired100C" + to_string(coincCheck[k].second))->Fill( p.second );}
								else{std::cout << p.first << " doesn't pass with " << coincCheck[k].first << std::endl;}*/
							if( !(p.second.first == coincCheck[k].second.first && p.first == coincCheck[k].first) )								h3d->Fill(p.second.first,coincCheck[k].second.first, coincCheck[k].first - p.first);
								if( f->helper->neighbour(p.second.first, coincCheck[k].second.first) )	neighbour++;
								k++;
							}

//if(neighbour > 1)								std::cout << coincCheck.size() << "A" << k << "N:" << neighbour << std::endl;
							if( neighbour == 0 )	f->getH2("pixelsCleaned0")->Fill((p.second.first-1)%8, floor((p.second.first-1)/8));	//	fill hit pattern
							if( neighbour > 0 )	f->getH2("pixelsCleaned1")->Fill((p.second.first-1)%8, floor((p.second.first-1)/8));	//	fill hit pattern
							if( neighbour > 1 )	{
								f->getH2("pixelsCleaned2")->Fill((p.second.first-1)%8, floor((p.second.first-1)/8));	//	fill hit pattern
								f->getH2("neutronMfluxatrx0")->Fill(p.second.second.first,p.second.second.second);	//	intensity of neutrons over positions
							}
							if( neighbour > 2 )	f->getH2("pixelsCleaned3")->Fill((p.second.first-1)%8, floor((p.second.first-1)/8));	//	fill hit pattern
							if( neighbour > 3 )	f->getH2("pixelsCleaned4")->Fill((p.second.first-1)%8, floor((p.second.first-1)/8));	//	fill hit pattern
							if( neighbour > 4 )	f->getH2("pixelsCleaned5")->Fill((p.second.first-1)%8, floor((p.second.first-1)/8));	//	fill hit pattern
							if( neighbour > 5 )	f->getH2("pixelsCleaned6")->Fill((p.second.first-1)%8, floor((p.second.first-1)/8));	//	fill hit pattern
							if( neighbour > 6 )	f->getH2("pixelsCleaned7")->Fill((p.second.first-1)%8, floor((p.second.first-1)/8));	//	fill hit pattern
							if( neighbour > 7 )	f->getH2("pixelsCleaned8")->Fill((p.second.first-1)%8, floor((p.second.first-1)/8));	//	fill hit pattern
							j++;

						}

						coincidences.clear();	//	reset coincidences vector
						coincCheck.clear();	//	reset pixels in coincidence

					}	//	end check of new beam pulse

					timeLastBp = timeLastB;	//	store this beam pulse time for a valid neutron for comparisons

				}	//	end gate check

			}	//	end neutron match
			else
			{	//	do not apply extra conditions

				f->getH1("nrjC" + s)->Fill(nrj);	//	add to nrj histogram
				f->getH1("nrj2C" + s)->Fill(nrj2);	//	add to nrj2 histogram

			}	//	end no neutron match

			f->getH2("timeVchan0")->Fill((timeb-timeLastB)%bpt, il);	//	channel V time
			f->getH2("nrjVchan0")->Fill(nrj, il);	//	nrj V time
			f->getH2("nrj2Vchan0")->Fill(nrj2, il);	//	nrj2 V time

			if(	il == chanX )	{	//	do stuff with the x position
				
				f->getH2("xVtime0")->Fill(timeb, nrj);	//	x position V time
				
			}	//	end doing stuff with x

			if(	il == chanR )	{	//	do stuff with the theta position
				
				f->getH2("rVtime0")->Fill(timeb, nrj);	//	theta position V time
				
			}	//	end doing stuff with theta

			f->getH1("timeadjC" + s)->Fill(timeSinceB);	//	add to adjusted time histogram
			f->getH1("nrjadjC" + s)->Fill( c->adjE(nrj) );	//	add to adjusted nrj histogram
			f->getH1("nrj2adjC" + s)->Fill( c->adjE(nrj2) );	//	add to adjusted nrj2 histogram
			f->getH2("timeadjVchan0")->Fill(timeSinceB, il);	//	channel V channel
			f->getH2("nrjadjVchan0")->Fill(c->adjE(nrj), il);	//	energy V channel
			f->getH2("nrj2adjVchan0")->Fill(c->adjE(nrj2), il);	//	nrj2 V channel
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

		for( Int_t i = 0; i <= gmc; i++ ) {	//	for each channel number

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
	std::cout << "Colourful hit pattern is done!" << std::endl;

}	//	end colourful hit pattern


}	//	end namespace fnt


void analysis(){};	//	allow root compilation