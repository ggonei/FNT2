//
//	George O'Neill @ University of York, 2020/02/18
//
//	This file is by design a simple separator for analysis
//
#include "../include/analysis.h"

namespace fnt {	//	create unique working area


void analysis::do_analysis(fnt::FNT* f) {	//	start analysis

	Channel* c;	//	allocate channel pointer
	Int_t npixel, timeSinceB, k;	//	get pixel, reset value for progress counter, counter
	std::string s;	//	string for names
	f->helper->resetCountdown();	//	reset countdown in case it has been used prior
	ULong64_t timeLastBp = 0, tarr[gmc];	//	previous beam pulse time
	vector<pair<Long64_t,Int_t>> coincidences;	//	hold coincident timestamps
	deque<pair<Long64_t,Int_t>> coincCheck;

	for( ULong64_t i = 0; i < n; i++ ) {	//	loop over all entries

		f->helper->countdown();	//	print progress
		bigTree->GetEntry(i);	//	grab energy info
		il = (Int_t) label;	//	convert label to integer
		c = f->getChannel(il);	//	set channel
		f->getH1("HitPattern0")->Fill(il);	//	hit pattern

		if( c ) {	//	skip unknown channels

			timeb = timeb + timeOffset;	//	adjust timeb to include offset
			s = to_string(il);	//	get channel number stored as string
			f->getH1("timeC" + s)->Fill((timeb-timeLastB)%bpt);	//	add to time histogram
			timeSinceB = (timeb-timeLastB-c->getTOffset())%bpt;

			npixel = c->getPixelNumber();	//	store pixel

			if( npixel > 0 )	{	//	neutron gate alternatively il>10&&il<90&&(il%10>0&&il%10<9)

				f->getH2("timeVneutrons0")->Fill((timeb-timeLastB)%bpt, il);	//	neutron channels V time

				if( c->passes(timeSinceB, 't') )	{	//	if passing channel gate

					f->getH2("pixels0")->Fill(floor(npixel/8), npixel%8);	//	fill hit pattern
					f->getH1("HitPatternNeutrons0")->Fill(il);	//	neutrons hit pattern
					f->getH1("neutronfluxatx0")->Fill(xPos);	//	intensity of neutrons at x
					f->getH1("neutronfluxatr0")->Fill(rPos);	//	intensity of neutrons at r
					f->getH2("neutronfluxatrx0")->Fill(xPos,rPos);	//	intensity of neutrons over positions
					f->getH2("timeadjneutronsVchan0")->Fill(timeSinceB, il);	//	time adjusted V channel
					f->getH1("nrjC" + s)->Fill(nrj);	//	add to nrj histogram
					f->getH1("nrj2C" + s)->Fill(nrj2);	//	add to nrj2 histogram

					if( timeLastB == timeLastBp )	{	//	if this entry belongs to current beam pulse

						if( timeLastB > 0 && npixel < 61 )	coincidences.push_back(make_pair((Long64_t) timeSinceB, npixel));	//	add time and pixel to vector

					}	//	end check if beam pulse is the same
					else	{	//	a new beam pulse

						sort(coincidences.begin(),coincidences.end());	//	sort vector

						for( UInt_t j = 0; j < coincidences.size(); j++ )	{	//	for each entry in this beam pulse

							tarr[coincidences[j].second] = coincidences[j].first;	//	store the time for this pixel in an array							
							k = 1;	//	
							coincCheck.push_back(coincidences[j]);

							while( coincCheck[0].first < coincidences[j].first - 100000 )	coincCheck.pop_front();

							for( pair<Long64_t,Int_t> p : coincCheck )	{

								f->getH2("coincidences0")->Fill(p.second,coincidences[j].second);
								if( p.first > coincidences[j].first - 10000 && p.first < coincidences[j].first + 10000 )	f->getH1("pixelsfired10C" + to_string(coincidences[j].second))->Fill( p.second );
								else if( p.first > coincidences[j].first - 20000 && p.first < coincidences[j].first + 20000 )	f->getH1("pixelsfired20C" + to_string(coincidences[j].second))->Fill( p.second );
								else if( p.first > coincidences[j].first - 50000 && p.first < coincidences[j].first + 50000 )	f->getH1("pixelsfired50C" + to_string(coincidences[j].second))->Fill( p.second );
								else if( p.first > coincidences[j].first - 100000 && p.first < coincidences[j].first + 100000 )	f->getH1("pixelsfired100C" + to_string(coincidences[j].second))->Fill( p.second );

							}

							for( ULong64_t tpiece : tarr )	{	//	for each time

								f->getH2("multiplicityt0")->Fill(coincidences[j].first - tpiece, coincidences[j].second);	//	fill time map

							}	//	end for loop over time array

						}	//	end for loop over coincidences
						
						coincidences.clear();	//	reset coincidences vector
						coincCheck.clear();

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

	std::cout << std::endl;	//	output a clean line after the countdown

	TDirectory* d;	//	hold value for directories
	TH1F* h;	//	hold value for histogram to move
	const Int_t tnbins = f->getH1("timeC0")->GetNbinsX();	//	number of bins in time spectra
	Double_t dest[tnbins], source[tnbins], *xpeaks;	//	set histogram arrays
	TSpectrum* peakfinder = new TSpectrum();	//	peak fitting spectrum
	ULong64_t peakcount = 0, binMultiplier;	//	number of peaks

	for( string folder : f->getFolders() )	{	//	for each folder

		std::cout << "Moving histograms to hists_" << folder << " folder..." << std::endl;	//	inform user
		d = newHists->mkdir(("hists_" + folder).c_str());	//	make a folder
		
		for( Int_t i = 0; i <= gmc; i++ ) {	//	for each channel number

			s = f->helper->sanitiser(folder + to_string(i));	//	create name
			h = (TH1F*)newHists->FindObjectAny(s.c_str());	//	histogram to move

			if( h ) {	//	if histogram exists

				if( folder == "timeC" )	{	//	if time spectrum

						for( Int_t j = 0; j < tnbins; j++ )	source[j] = h->GetBinContent(j + 1);	//	set source array from histogram
						peakfinder->SearchHighRes(source, dest, tnbins, 10, 10, kFALSE, 2, kTRUE, 10);	//	find peaks
						xpeaks = peakfinder->GetPositionX();	//	get peak positions

						if( peakfinder->GetNPeaks() > 0 )	{	//	if there are peaks

							binMultiplier = (abs(h->GetXaxis()->GetXmin()) + abs(h->GetXaxis()->GetXmax())) / tnbins;
							std::cout << "Neutron time offset for channel " << i << ":\t" << xpeaks[0] * binMultiplier;	//	print neutron offset
							if( sizeof(xpeaks)/sizeof(Double_t*) > 1)	std::cout << ", gamma time offset" << xpeaks[1] * binMultiplier;	//	print gamma offset
							std::cout << std::endl;	//	end line
							xpeaks = {};	//	empty xpeaks
							delete xpeaks;	//	reset xpeaks
							delete peakfinder;	//	release peakfinder memory
							peakfinder = new TSpectrum();	//	add blank spectrum

						}	//	end peaks check

				}	//	end time spectrum specialisations

				h->SetDirectory(d);	//	move histogram
				newHists->Delete((s + ";1").c_str());	//	delete from original folder

			}	//	end move of histogram
			
		}	//	end for loop over each histogram

	}	//	end for loop over each folder

	newHists->Write("", TObject::kOverwrite);	//	write histograms to file
	
}	//	end do_analysis



void analysis::colourful_hp(TTree* bigTree) {	//	output a colourful hit pattern

	std::cout << "Starting colourful hit pattern now...it's not efficient so it will take some time, I think it is with nlogn not n!" << std::endl;
	Int_t bins[3] = {256, 0, 255};	//	histogram {bin count, X minimum, X maximum}
	Int_t ic[7] = {kGray, kBlue, kViolet, kTeal, kRed, kPink, kGreen};	//	histogram colours
	TCanvas* cA = new TCanvas();	//	drawing canvas for histogram
	TH1F* hist_arr7[7];	//	histogram array
	TLegend* l = new TLegend(0.7,0.7,1,1);	//	initialise legend
	bool tf[256] = { false };	//	initialise an array for all channels
	string sh[7] = {"All Channels", "Neutrons", "Germaniums", "LaBr3", "Nothing", "RF", "Table"};	//	histogram titles
	string sc[7] = {"", "((label%10>0&&label%10<9)&&label>10&&label<90)||label==2", "label==90||label==91", "label==3",  "label==0||label==1||label==4||label==7||label==8||label==10||label==19||label==20||label==29||label==30||label==39||label==40||label==49||label==50||label==59||label==60||label==234||label==235||label==243||label==244||label==245||label==246||label==247||label==248||label==249||label==250||label==253||label==254||label==255", "label==97", "label==95||label==96"};	//	histogram conditions GO


	for (UInt_t i = 0; i < (sizeof(hist_arr7)/sizeof(TH1F*)); i++) {	//	loop over histogram array

		hist_arr7[i] = new TH1F(("hist_arr[" + to_string(i) + "]").c_str(), sh[i].c_str(), bins[0], bins[1], bins[2]);	//	create new histogram
		hist_arr7[i]->SetLineColor(ic[i]);	//	set histogram colour
		bigTree->Draw(("label>>hist_arr[" + to_string(i) + "]").c_str(), sc[i].c_str());	//	grab data from tree and add to histogram
		l->AddEntry(hist_arr7[i], sh[i].c_str());	//	add histogram entry to legend with title

	}	//	end loop over histogram array


	cA->SetLogy();	//	logarithmic scale
	std::cout << "Colourful hit pattern is done!" << std::endl;

}	//	end colourful hit pattern


}	//	end namespace fnt


void analysis(){};	//	allow root compilation