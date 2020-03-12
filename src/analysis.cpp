//
//	George O'Neill @ University of York, 2020/02/18
//
//	This file is by design a simple separator for analysis
//
#include "../include/analysis.h"

namespace fnt {	//	create unique working area


void analysis::do_analysis(fnt::FNT* f) {	//	start analysis

	Channel* c;	//	allocate channel pointer
	Int_t npixel;	//	get reset value for progress counter

	for( ULong64_t i = 0; i < n; i++ ) {	//	loop over all entries

		f->helper->countdown();	//	print progress
		bigTree->GetEntry(i);	//	grab energy info
		il = (Int_t) label;	//	convert label to integer
		c = f->getChannel(il);	//	set channel
		f->getH1("HitPattern0")->Fill(il);	//	hit pattern

		if( c ) {	//	skip unknown channels

			f->getH1("timeC" + to_string(il))->Fill((timeb-timeLastB)%bpt);	//	add to time histogram
			f->getH1("nrjC" + to_string(il))->Fill(nrj);	//	add to nrj histogram
			f->getH1("nrj2C" + to_string(il))->Fill(nrj2);	//	add to nrj2 histogram

			npixel = c->getPixelNumber();	//	store pixel

			if( npixel > -1 )	{	//	neutron gate alternatively il>10&&il<90&&(il%10>0&&il%10<9)

				f->getH2("timeVneutrons")->Fill((timeb-timeLastB)%bpt, il);	//	neutron channels V time
				f->getH2("pixels")->Fill(floor(npixel/8), npixel%8);	//	fill hit pattern
if(il==55)	{
				f->getH1("HitPatternNeutrons0")->Fill(il);	//	neutrons hit pattern
				f->getH1("neutronfluxatx0")->Fill(xPos);	//	intensity of neutrons at x
				f->getH1("neutronfluxatr0")->Fill(rPos);	//	intensity of neutrons at r
				f->getH2("neutronfluxatrx")->Fill(xPos,rPos);	//	intensity of neutrons over positions
}

			}	//	end neutron match

			f->getH2("timeVchan")->Fill((timeb-timeLastB)%bpt, il);	//	channel V time
			f->getH2("nrjVchan")->Fill(nrj, il);	//	nrj V time
			f->getH2("nrj2Vchan")->Fill(nrj2, il);	//	nrj2 V time

			if(	il == chanX )	{	//	do stuff with the x position
				
				f->getH2("xVtime")->Fill(timeb, nrj);	//	x position V time
				
			}	//	end doing stuff with x

			if(	il == chanR )	{	//	do stuff with the theta position
				
				f->getH2("rVtime")->Fill(timeb, nrj);	//	theta position V time
				
			}	//	end doing stuff with theta

			f->getH1("timeadjC" + to_string(il))->Fill((timeb-timeLastB-(c->getTOffset()))%bpt);	//	add to adjusted time histogram
			f->getH1("nrjadjC" + to_string(il))->Fill( c->adjE(nrj) );	//	add to adjusted nrj histogram
			f->getH1("nrj2adjC" + to_string(il))->Fill( c->adjE(nrj2) );	//	add to adjusted nrj2 histogram
			f->getH2("timeadjVchan")->Fill((timeb-timeLastB-(c->getTOffset()))%bpt, il);	//	channel V time
			f->getH2("nrjadjVchan")->Fill(c->adjE(nrj), il);	//	nrj V time
			f->getH2("nrj2adjVchan")->Fill(c->adjE(nrj2), il);	//	nrj2 V time

		}	//	end valid channel check

	}	//	end loop over all entries

	std::cout << std::endl;	//	output a clean line after the countdown

	TDirectory* d;	//	hold value for directories
	TH1F* h;	//	hold value for histogram to move
	std::string s;	//	string for histogram name

	for( string folder : f->getFolders() )	{	//	for each folder

		std::cout << "Moving histograms to " << folder << "folder ..." << std::endl;	//	inform user
		d = newHists->mkdir(("hists_" + folder).c_str());	//	make a folder
		
		for( Int_t i = 0; i <= gmc; i++ ) {	//	for each channel number

			s = f->helper->sanitiser(folder + to_string(i));	//	create name
			h = (TH1F*)newHists->FindObjectAny(s.c_str());	//	histogram to move

			if( h ) {	//	if histogram exists

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