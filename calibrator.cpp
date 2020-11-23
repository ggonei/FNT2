/*	Calibrations for FNT2:
152Eu:
	90:
		Energy Calibration:
		0:	-0.0547329, 0.00102646	1:	0.00886368, 2.02635e-08
		Efficiency Calibration:
		0:	-137.338, 0.0962817			1:	77.3644, 0.0206774			2:	-9.80473, 0.00323811	3:	-0.868683, 0.000469529	4:	0.262858, 6.46264e-05   	5:	-0.0144665, 8.14819e-06 	6:	2.91573e+06
	91:
		Energy Calibration:
		0:	0.572201, 0.00126262		1:	0.0220622, 5.73668e-08
		Efficiency Calibration:
		0:	-421.135, 0.0890989			1:	219.271, 0.0189043			2:	-27.5265, 0.00294783	3:	-2.445, 0.000426333			4:	0.735513, 5.85643e-05   	5:	-0.0403394, 7.38059e-06 	6:	5.09471e+06
	3:
		Energy Calibration:
		0:	25.8179, 2.3254					1:	0.0249227, 0.000603929
		Efficiency Calibration:
		0:	-421.588, 0.592202			1:	219.177, 0.123318				2:	-27.5462, 0.0256794		3:	-2.4491, 0.00534738			4:	0.734661, 0.00111352			5:	-0.0405168, 0.000231876 	6:	176645

60Co:
	90:
		Energy Calibration:
		0:	0.430319, 0.0510943			1:	0.00885946, 3.65704e-07
		Efficiency Calibration:
		0:	7.91006, 0.0769945			1:	0.795861, 0.0110989			2:	0.0676043, 0.00157713 3:	0.0033063, 0.000221151	4:	-0.000402029, 3.06014e-05	5:	-0.000180786, 4.17418e-06	6:	1.24287e+06
	91:
		Energy Calibration:
		0:	-2.19462, 0.0515678			1:	0.0221174, 9.19492e-07
		Efficiency Calibration:
		0:	7.94413, 0.0829718			1:	0.797105, 0.0119612			2:	0.0672829, 0.00169977 3:	0.00319166, 0.000238362	4:	-0.000427853, 3.29848e-05	5:	-0.000185772, 4.49956e-06	6:	1.11757e+06
	3:
		Energy Calibration:
		0:	1323.32, 28.4498				1:	0.0294838, 0.0913926
		Efficiency Calibration:
		0:	7.61486, 0.335289				1:	0.751341, 0.0466015			2:	0.0609222, 0.0064771	3:	0.00230759, 0.000900247	4:	-0.000550729, 0.000125124	5:	-0.00020285, 1.73909e-05	6:	275228

252Cf:
	Looks like 152Eu
	15, 25, 34, 35, 41, 42, 43, 44, 53, 54, 55, 56, 57, 58, 64, 65, 74, 84 threshold is too high
	45, 46 are bad but not as bad, maybe different gains
	73, 75, 76, 77 have some noise in 'tails'
	85, 86, 87, 88 show basically nothing
	nrj2 are all lower by at least 10^3 with no patterns on above except 85-88
	range is correct
*/
template <typename T> double sgn(T v) {return (0. < v) - (v < 0.);}	//	typesafe sgn for real numbers for quartic


vector<std::pair<TGraphErrors*,TFitResultPtr>> calibrator() {

	const double interval = 0.005, threshold = 0.025, sigma = 5;	//	interval to look for similar peaks, threshold, sigma for peak finder
	const string prefix = "hists_nrjC", channels[] = {"90","91","3"};	//	histogram folder, histogram - prefix list
	const int order = 1, width = 500;	//	channel numbers, calibration polynomial order, half width of peak to fit over
	const int maxi = sizeof(channels)/sizeof(channels[0]);	//	store loop maximum
	double *npeaks, entryi, calibkey, sce, scd, tpx, frp[5], maxp;	//	peaks found in spectrum, first entry, key value, position error, integral error, log log efficiency maximum, quartic constants, efficiency maxmimum
	std::complex<double> quart;	//	complex quartic solutions
	int maxidx;	//	maximum count in vector
	short peaksTotal, adv, dotpos, declen, gpos;	//	number of peaks, determined peak, decimal placement, length of decimal string, point to change
	std::string scv, sci;	//	string value, integral
	TH1D *hFind;	//	working histogram
	TF1 *myfnl	= new TF1("myfnl", "[0]+[1]*TMath::Power(x,1)+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)+[5]*TMath::Power(x,5)");	//	efficiency function
	TF1 *myfnlp = new TF1("myfnlp","TMath::Power(TMath::Exp(1),[0]+[1]*TMath::Power(TMath::Log(x),1)+[2]*TMath::Power(TMath::Log(x),2)+[3]*TMath::Power(TMath::Log(x),3)+[4]*TMath::Power(TMath::Log(x),4)+[5]*TMath::Power(TMath::Log(x),5))/[6]",0,2000);	//	Efficiency curve
	TSpectrum *s;	//	workig spectrum for peak finding
	std::vector<pair<double,double>> means[maxi], integrals[maxi];	//	peak information
	TFitResultPtr fr, fm;	//	temporary fit result for efficiency, energy
	std::vector<TFitResultPtr> fitresult;	//	fitted peak parameters
	TF1 *gausfit;	//	fitted peak
	TGraphErrors *ge, *gm, *gi;	//	efficiency plot, energy calibration plot, relative intensity plot
	std::pair<std::map<double,std::pair<int,int>>::iterator,bool> ret;	//	duplicated entry
	std::map<double,std::pair<int,int>> calibstdpmap;	//	calibration percentages
	std::map<std::pair<std::pair<double,double>,std::pair<double,double>>,std::pair<std::pair<string,int>,std::pair<string,int>>> calib[maxi];	//	calibration term array
	std::map<double,int> meansl[maxi];	//	list of means and errors
	std::map<int,std::pair<TFitResultPtr,TFitResultPtr>> fitresults;	//	list of fits
	std::map<double,int>::iterator itub, itlb;	//	iterators for upper and lower bound for assignment
	std::vector<double> meansp[maxi];	//	list of means as percentages
	std::vector<std::pair<std::pair<int,int>,double>> meanspp[maxi];	//	list of means as percentages for relevant peaks
	std::vector<bool> t[maxi];	//	truth vector to determine cut off in matching
	std::vector<int> idx;	//	most probable peak matcher
	vector<std::pair<TGraphErrors*,TFitResultPtr>> hrv;	//	return plots

	struct sdcomp:public std::binary_function<std::pair<string,int>,std::pair<string,int>,bool> {	//	struct for map
		const bool operator()(const std::pair<string,int>& l,const std::pair<string,int>& r) {	//	comparator for struct
			return std::pair<float,int>(std::stod(l.first),l.second)<std::pair<float,int>(std::stod(r.first),r.second);	//	use efficiency without error to assign
		}	//	end comparator
	};	//	end sdcomp{}

	//	calibstd = {	{frequency/%, error}, {reference energy/keV, error}	}
	//	60Co
	/*const std::map<std::pair<string,int>,std::pair<string,int>,sdcomp> calibstd = {
		{{"0.0076",	5},	{"346.93",	7}},
		{{"0.0076",	8},	{"826.06",	3}},
		{{"99.9736",	7},	{"1173.237",	4}},
		{{"99.9856",	4},	{"1332.501",	5}},
		{{"0.00111",	18},	{"2158.57",	10}},
		{{"0.000002",	4},	{"2505",	0}},
	};*/
	//	152Eu
	const std::map<std::pair<string,int>,std::pair<string,int>,sdcomp> calibstd = {
		{{"28.58",	6},	{"121.7817",	3}},
		{{"26.5",	4},	{"344.2785",	12}},
		{{"21.005",	24},	{"1408.006",	3}},
		{{"14.605",	21},	{"964.079",	18}},
		{{"13.644",	21},	{"1112.074",	4}},
		{{"12.942",	19},	{"778.9040",	18}},
		{{"10.207",	21},	{"1085.869",	24}},
		{{"7.583",	19},	{"244.6975",	8}},
		{{"4.245",	19},	{"867.378",	4}},
		{{"2.821",	19},	{"443.965",	3}},
		{{"2.234",	4},	{"411.1163",	11}},
		{{"1.727",	6},	{"1089.737",	5}},
		{{"1.623",	8},	{"1299.140",	10}},
		{{"1.422",	6},	{"1212.948",	11}},
		{{"0.861",	5},	{"367.7887",	16}},
		{{"0.857",	8},	{"688.670",	5}},
		{{"0.646",	5},	{"1005.272",	17}},
		{{"0.502",	5},	{"1457.643",	11}},
		{{"0.489",	6},	{"563.990",	7}},
		{{"0.471",	4},	{"678.623",	5}},
		{{"0.459",	5},	{"586.2648",	25}},
		{{"0.447",	5},	{"295.9392",	17}},
		{{"0.427",	6},	{"919.330",	3}},
		{{"0.419",	3},	{"488.6792",	20}},
		{{"0.327",	19},	{"443.96",	4}},
		{{"0.320",	3},	{"810.451",	5}},
		{{"0.281",	5},	{"1528.103",	18}},
		{{"0.278",	8},	{"719.349",	4}},
		{{"0.278",	5},	{"926.317",	15}},
		{{"0.246",	8},	{"1084",	1}},
		{{"0.215",	24},	{"764.900",	9}},
		{{"0.188",	4},	{"1249.938",	13}},
		{{"0.186",	8},	{"1109.174",	12}},
		{{"0.172",	4},	{"674.675",	3}},
		{{"0.1660",	24},	{"841.570",	5}},
		{{"0.149",	8},	{"503.474",	5}},
		{{"0.1448",	19},	{"656.487",	5}},
		{{"0.135",	3},	{"963.390",	12}},
		{{"0.1290",	19},	{"566.430",	5}},
		{{"0.128",	8},	{"329.425",	21}},
		{{"0.1100",	19},	{"416.048",	8}},
		{{"0.105",	5},	{"1292.778",	19}},
		{{"0.094",	4},	{"768.944",	9}},
		{{"0.093",	8},	{"712.843",	6}},
		{{"0.086",	4},	{"901.181",	11}},
		{{"0.0729",	21},	{"271.131",	8}},
		{{"0.0729",	19},	{"930.580",	15}},
		{{"0.072",	3},	{"324.83",	3}},
		{{"0.072",	3},	{"251.630",	7}},
		{{"0.0669",	21},	{"896.58",	9}},
		{{"0.059",	8},	{"719.349",	4}},
		{{"0.052",	4},	{"520.227",	5}},
		{{"0.0507",	13},	{"315.174",	17}},
		{{"0.0427",	11},	{"534.245",	7}},
		{{"0.037",	4},	{"148.010",	17}},
		{{"0.037",	3},	{"1170.93",	11}},
		{{"0.036",	3},	{"340.40",	14}},
		{{"0.0335",	21},	{"275.449",	15}},
		{{"0.0334",	13},	{"1261.343",	23}},
		{{"0.032",	11},	{"595.61",	12}},
		{{"0.0313",	13},	{"990.19",	3}},
		{{"0.0294",	21},	{"493.508",	20}},
		{{"0.0292",	24},	{"482.31",	3}},
		{{"0.0263",	21},	{"794.81",	3}},
		{{"0.0257",	11},	{"1363.77",	5}},
		{{"0.0217",	24},	{"958.63",	5}},
		{{"0.0198",	5},	{"212.568",	15}},
		{{"0.0195",	16},	{"671.151",	13}},
		{{"0.0193",	16},	{"686.61",	5}},
		{{"0.0185",	13},	{"556.56",	3}},
		{{"0.0178",	11},	{"1348.10",	7}},
		{{"0.0170",	16},	{"675.1",	2}},
		{{"0.0166",	24},	{"664.78",	5}},
		{{"0.0163",	16},	{"905.9",	1}},
		{{"0.0161",	11},	{"839.36",	4}},
		{{"0.016",	8},	{"696.87",	19}},
		{{"0.016",	5},	{"125.69",	13}},
		{{"0.0150",	16},	{"523.13",	5}},
		{{"0.0141",	8},	{"974.09",	4}},
		{{"0.0141",	11},	{"1206.11",	15}},
		{{"0.0134",	16},	{"440.86",	10}},
		{{"0.0133",	16},	{"441.1",	4}},
		{{"0.0131",	6},	{"526.881",	20}},
		{{"0.0123",	11},	{"805.70",	7}},
		{{"0.0112",	8},	{"727.99",	14}},
		{{"0.0099",	8},	{"285.98",	3}},
		{{"0.0098",	11},	{"494.0",	3}},
		{{"0.0096",	4},	{"1769.09",	5}},
		{{"0.0094",	13},	{"237.31",	5}},
		{{"0.0093",	13},	{"351.66",	4}},
		{{"0.0091",	8},	{"616.05",	3}},
		{{"0.0083",	8},	{"269.86",	6}},
		{{"0.0082",	5},	{"1605.61",	7}},
		{{"0.008",	3},	{"395.75",	19}},
		{{"0.008",	3},	{"173.17",	15}},
		{{"0.0072",	16},	{"1315.31",	23}},
		{{"0.00679",	21},	{"192.60",	4}},
		{{"0.0063",	6},	{"1647.41",	14}},
		{{"0.0062",	8},	{"644.37",	5}},
		{{"0.0062",	8},	{"1674.30",	6}},
		{{"0.0061",	13},	{"195.05",	24}},
		{{"0.0059",	7},	{"207.6",	3}},
		{{"0.0059",	19},	{"1698.1",	4}},
		{{"0.0059",	16},	{"330.54",	10}},
		{{"0.0059",	11},	{"735.40",	10}},
		{{"0.0056",	24},	{"1485.9",	3}},
		{{"0.0055",	5},	{"209.41",	13}},
		{{"0.0053",	4},	{"1608.36",	8}},
		{{"0.0051",	8},	{"496.39",	3}},
		{{"0.0051",	8},	{"385.69",	20}},
		{{"0.0051",	11},	{"202.74",	13}},
		{{"0.0050",	11},	{"1314.67",	1}},
		{{"0.0048",	8},	{"571.6",	5}},
		{{"0.0048",	8},	{"1390.36",	16}},
		{{"0.0046",	8},	{"1001.1",	3}},
		{{"0.0046",	19},	{"239.42",	17}},
		{{"0.0045",	19},	{"557.91",	17}},
		{{"0.0042",	4},	{"496.39",	3}},
		{{"0.0040",	8},	{"357.26",	5}},
		{{"0.0040",	11},	{"175.18",	0}},
		{{"0.0039",	6},	{"538.29",	6}},
		{{"0.0037",	24},	{"968.00",	0}},
		{{"0.0035",	13},	{"389.07",	11}},
		{{"0.0034",	8},	{"703.25",	6}},
		{{"0.0034",	13},	{"937.05",	15}},
		{{"0.0032",	8},	{"683.32",	11}},
		{{"0.0032",	6},	{"423.45",	4}},
		{{"0.00294",	21},	{"387.90",	8}},
		{{"0.00292",	21},	{"387.90",	8}},
		{{"0.0027",	13},	{"562.93",	2}},
		{{"0.0021",	13},	{"316.3",	1}},
		{{"0.0016",	5},	{"703.54",	5}},
		{{"0.0016",	5},	{"320.03",	15}},
		{{"0.0014",	6},	{"482.3",	4}},
		{{"0.00126",	8},	{"1139",	1}},
		{{"0.00126",	21},	{"391.32",	14}},
		{{"0.00104",	21},	{"561.2",	5}},
		{{"0.00083",	21},	{"406.74",	15}},
		{{"0.00083",	21},	{"379.37",	6}},
		{{"0.00016",	5},	{"1635.2",	5}},
	};

	std::map<std::pair<string,int>,std::pair<string,int>>::const_reverse_iterator calibstdi = calibstd.rbegin();	//	iterator for standard

	if( prefix != "" ) {	//	if we need to dive into folder
		((TDirectory*) gROOT->FindObject(prefix.c_str()))->cd();	//	move to that folder
	}	//	close folder check

	for(int i = 0; i < maxi; i++) {	//	loop over channels to calibrate
		s = new TSpectrum();	//	peak fitting spectrum
		hFind = (TH1D*) (gROOT->FindObject(((prefix != "" ? prefix.substr(6) : "") + channels[i]).c_str()))->Clone();	//	make a copy of our original spectrum
		hFind->Add(s->Background(hFind, 20, "Compton"), -1);	//	subtract background manually
		peaksTotal = s->Search(hFind, sigma, "goff nodraw nobackground noMarkov", threshold);	//	find peaks
		npeaks = s->GetPositionX();	//	Get all peak positions

		for(int k = 0; k < peaksTotal; k++) {	//	loop over peaks
			fr = hFind->Fit("gaus", "QEMNS", "", npeaks[k] - width, npeaks[k] + width);	//	fit found peak
			means[i].push_back({fr->Parameter(1),fr->Error(1)});	//	store peak mean
			gausfit = (TF1*)gROOT->FindObject("gaus");	//	get fitted peak function
			integrals[i].push_back({gausfit->Integral(npeaks[k] - width,npeaks[k] + width),gausfit->IntegralError(npeaks[k] - width,npeaks[k] + width)});	//	determine area
		}	//	end loop over peaks

	}	//	end loop over channels

	for(int j = 0; j < calibstd.size(); j++) {	//	loop over all calibration standard entries

		calibstdi = calibstd.rbegin();	//	set iterator to most frequent peak in calibration standard
		std::advance(calibstdi,j);	//	set iterator to next entry in standard
		entryi = std::stod(calibstdi->second.first);	//	get numeric value of standard entry
		calibstdi = calibstd.rbegin();	//	move iterator back to most frequent peak in standard

		for(int k = 0; k < calibstd.size(); k++) {	//	loop over all calibration standard entries
		//	duplicates will exist but they are all at high n so do not concern us

			calibkey = std::stod(calibstdi->second.first) / entryi;	//	determine percentage at this level
			ret = calibstdpmap.insert({calibkey, {j, k}});	//	add percentage along with the peak that is highest and this peak position

			if(ret.second == false)	//	if our entry fails
				std::cout << "element already existed with a value of " << calibkey << " (" << j << "," << k << ")" << endl;	//	inform the user

			std::advance(calibstdi,1);	//	move to next entry in standard

		}

	}

	//	we make all our means percentages
	//	then with that set look for all matching percentages over a certain range
	//	then scan that for the most populous peak
	//	with this vector we then look for peaks which we are certain about because we have enough information to form a continuous set
	//	and then delete all others
	//	then iterate over the set to match known percentages with our continuous set that matches possible percentages
	//	and use the fact the first result has least error and is the one closest to the real percentage we ensure another peak cannot match

	for(int i = 0; i < maxi; i++) {	//	loop over all channels

		for(int j = 0; j < means[i].size(); j++)	//	loop over all found peaks
			meansp[i].push_back(means[i][j].first / means[i][0].first);	//	create a vector of percentages

		for(int j = 0; j < meansp[i].size(); j++)	//	loop over all found peak percentages
			for(std::map<double,std::pair<int,int>>::iterator it = calibstdpmap.lower_bound(meansp[i][j]*(1-interval)); it != calibstdpmap.lower_bound(meansp[i][j]*(1+interval)); ++it)	//	loop over bounded map
				meanspp[i].push_back({{it->second.first, it->second.second}, it->first});	//	add possible standard peak

		idx = vector<int>(means[i].size());	//	create empty vector based on number of found peaks

		for(int j = 0; j < meanspp[i].size(); j++)	//	loop over all standard percentages bounded by found peak percentages
			if(meanspp[i][j].first.first < means[i].size() && meanspp[i][j].first.first != meanspp[i][j].first.second)	//	increment vector of peak which corresponds to first
				idx[meanspp[i][j].first.first]++;	//	incrementation at position corresponding to highest peak before efficiency correction

		maxidx = std::max_element(idx.begin(),idx.end()) - idx.begin();	//	determine most frequent peak occurring first
		meansl[i].insert({1,maxidx});	//	add a value of 1 so vector matches itself since 0,0 will always take this position
		
		for(int j = 0; j < meanspp[i].size(); j++) {	//	loop over all standard percentages bounded by found peak percentages
			if(meanspp[i][j].first.first == maxidx) {	//	if state matches what we are looking for
				meansl[i].insert({meanspp[i][j].second, meanspp[i][j].first.second});	//	insert value to be used for matching

		t[i] = vector<bool>(meansl[i].size());	//	initialise a truth vector to false
		
		for(int j=0; j < meansl[i].size(); j++) {	//	loop over matching values

			for(std::map<double,int>::const_iterator it = meansl[i].begin(); it != meansl[i].end(); it++)	//	loop over all entries in map
				if((*it).second == j) {	//	if we find a consecutive value

					t[i][j] = true;	//	change value in truth vector
					break;	//	exit to loop at next j
					
				}	//	end consecutive value check

			if(!t[i][j])	//	if we exit without triggering a change in truth
				break;	//	skip checking higher numbers

		}	//	end loop over matching values

		for(std::map<double,int>::const_iterator it = meansl[i].begin(); it != meansl[i].end();)	//	loop over matching values
			if(!t[i][(*it).second])	//	if our truth entry is false
				meansl[i].erase(it++);	//	delete value to be matched
			else ++it;	//	if our truth entry is true then just increment iterator

		for(int j = 0; j < meansp[i].size(); j++) {	//	for found peak percentages

			calibstdi = calibstd.rbegin();	//	set calibration standard to start
			itub = meansl[i].lower_bound(meansp[i][j]);	//	determine value to be matched from logically ordered map
			itlb = itub;	//	set lower bound equal to upper bound

			if(itub != meansl[i].begin())	//	if we are not at the start of the map
				--itlb;	//	decrease lower bound to previous entry

			if(itub == meansl[i].end())	//	if we are past the end of the map
				--itub;	//	go to end of the map

			//	the next line is what should be changed for fits deviating by more than 25% from a linear fit
			adv = ((t[i][(*itub).second] && (*itub).second <= (*itlb).second) || !t[i][(*itlb).second] ? (*itub).second : (*itlb).second);	//	determine closest percentage

			if(t[i][adv]) {	//	if this entry is true

				t[i][adv] = false;	//	make sure it cannot be matched again

				advance(calibstdi, adv);	//	move iterator to appropriate position
			//debug
			//	std::cout << j << "\t" << means[i][j].first << "\t" << adv << "\t(" << (*itub).second << ":" << (*itlb).second << ")\t" << calibstdi->second.first << "\t" << meansp[i][j] << "\t" << (*itub).first << "\t" << (*itlb).first << "\t" << abs((*itub).first-meansp[i][j]) << "\t" << abs((*itlb).first-meansp[i][j]) << std::endl;
			//	std::cout << j << "\t" << integrals[i][j].first << "\t" << adv << "\t(" << (*itub).second << ":" << (*itlb).second << ")\t" << calibstdi->first.first; << std::endl;
				calib[i].insert({{means[i][j], integrals[i][j]}, {(*calibstdi).second, (*calibstdi).first}});	//	insert found peak with actual peak to a calibration map

			}	//	end truth check

		}	//	end found peak percentage loop

	}	//	end channel loop for calibration map creation

	for(int i = 0; i < maxi; i++) {	//	for each channel

		gm = new TGraphErrors();	//	reset means graph
		gi = new TGraphErrors();	//	reset integrals graph

		for(std::map<std::pair<std::pair<double,double>,std::pair<double,double>>,std::pair<std::pair<string,int>,std::pair<string,int>>>::const_iterator it = calib[i].begin(); it != calib[i].end(); ++it) {	//	for each entry in calibration

			gpos = gm->GetN();	//	get number of entries made so far
			scv = (*it).second.first.first;	//	fetch string for calibration position
			sci = (*it).second.second.first;	//	fetch string for integral amount
			dotpos = scv.find(".") + 1;	//	find position of decimal in mean
			declen = scv.substr(dotpos).length();	//	find precision of mean
			sce = (*it).second.first.second / (dotpos > 0 && declen > 0 ? TMath::Power(10,declen) : 1);	//	calculate actual mean error from given integer
			dotpos = sci.find(".") + 1;	//	find position of decimal in integral
			declen = sci.substr(dotpos).length();	//	find precision of integral
			scd = (*it).second.second.second/(dotpos > 0 && declen > 0 ? TMath::Power(10,declen) : 1);	//	calculate actual integral error from given integer
			gm->SetPoint(gpos, (*it).first.first.first, std::stod(scv));	//	plot energy point
			gm->SetPointError(gpos, (*it).first.first.second, sce);	//	plot energy point error
			gi->SetPoint(gpos, TMath::Log(std::stod(scv)), TMath::Log((*it).first.second.first / std::stod(sci)));	//	plot efficiency point
			gi->SetPointError(gpos, gi->GetPointX(gpos) * sce / std::stod(scv), gi->GetPointY(gpos) * TMath::Sqrt(TMath::Power((*it).first.second.second / (*it).first.second.first, 2) + TMath::Power(scd / std::stod(sci), 2)));	//	plot efficiency point error

		}	//	end plot of each calibration point

		for(fr = 0; fr == 0 || fr->MinFcnValue() != fr->Chi2();)	//	loop whilst fit is awful
			fr = gi->Fit(myfnl,"QEMNS");	//	fit

		fm = gm->Fit(((string)"pol"+char(order+'0')).c_str(), "QEMNS");	//	store fit from energy calibration
		fitresults.insert({i, {fm, fr}});	//	store fits for use in calculation step
		hrv.push_back({gm, fm});	//	store energy fit for returning to user
		hrv.push_back({gi, fr});	//	store efficiency fit for returning to user

	}	//	end loop over all channels

	for(std::map<int,std::pair<TFitResultPtr,TFitResultPtr>>::const_iterator it = fitresults.begin(); it != fitresults.end(); ++it) {
		for(int i = 0; i < 5; i++)	frp[i] = (i + 1)*(*it).second.second->Parameter(i + 1);
		std::cout << "Quartic roots (if only one principle (= real, positive) root it is used automatically):" << std::endl;
		tpx = -1;
		//	this is an amazing equation for calculating all quartic roots
		//	it could be functions, or structs, or a class making use of cubics and quadratics etc
		//	but IT IS ONE LINE TO OUTPUT QUARTIC ROOTS, which as far as I can find has not been written down generally
		//	For quintic roots it obviously cannot be done, "I have discovered a truly marvelous proof of this, which these comments are too narrow to contain."
		for(double s = -1; s <= 1; s = s + 2)	for(double t = -1; t <= 1; t = t + 2) {	//	create loop for each root

			//	this equation could be made completely complex safe but pow(x,1./3) is less precise than cbrt (see x=343 with std::cout.precision(std::numeric_limits<double>::max_digits10);)
			quart = (-3.*frp[3]+s*(std::sqrt(complex<double>(
				3.*(3.*std::pow(frp[3],2.)-8.*frp[4]*frp[2]+2.*frp[4]*std::cbrt(4.*(2.*std::pow(frp[2],3.)-9.*frp[3]*frp[2]*frp[1]+27.*frp[4]*std::pow(frp[1],2.)+27.*std::pow(frp[3],2.)*frp[0]-72.*frp[4]*frp[2]*frp[0]+std::sqrt(complex<double>(
					std::pow((2.*std::pow(frp[2],3.)-9.*frp[3]*frp[2]*frp[1]+27.*frp[4]*std::pow(frp[1],2.)+27.*std::pow(frp[3],2.)*frp[0]-72.*frp[4]*frp[2]*frp[0]),2.)-4.*std::pow((std::pow(frp[2],2.)-3.*frp[3]*frp[1]+12.*frp[4]*frp[0]),3.)
				))).real())+2.*frp[4]*std::cbrt(4.*(2.*std::pow(frp[2],3.)-9.*frp[3]*frp[2]*frp[1]+27.*frp[4]*std::pow(frp[1],2.)+27.*std::pow(frp[3],2.)*frp[0]-72.*frp[4]*frp[2]*frp[0]-std::sqrt(complex<double>(
					std::pow((2.*std::pow(frp[2],3.)-9.*frp[3]*frp[2]*frp[1]+27.*frp[4]*std::pow(frp[1],2.)+27.*std::pow(frp[3],2.)*frp[0]-72.*frp[4]*frp[2]*frp[0]),2.)-4.*std::pow((std::pow(frp[2],2.)-3.*frp[3]*frp[1]+12.*frp[4]*frp[0]),3.)
				))).real()))
			))
			+t*std::sqrt(
				3.*(3.*std::pow(frp[3],2.)-
				
				8.*frp[4]*frp[2]+2.*frp[4]*(-1.+std::sqrt(complex<double>(-3.)))/2.*std::cbrt(4.*(2.*std::pow(frp[2],3.)-9.*frp[3]*frp[2]*frp[1]+27.*frp[4]*std::pow(frp[1],2.)+27.*std::pow(frp[3],2.)*frp[0]-72.*frp[4]*frp[2]*frp[0]+std::sqrt(complex<double>(
					std::pow((2.*std::pow(frp[2],3.)-9.*frp[3]*frp[2]*frp[1]+27.*frp[4]*std::pow(frp[1],2.)+27.*std::pow(frp[3],2.)*frp[0]-72.*frp[4]*frp[2]*frp[0]),2.)-4.*std::pow((std::pow(frp[2],2.)-3.*frp[3]*frp[1]+12.*frp[4]*frp[0]),3.)
				))).real())
				
				+
				
				2.*frp[4]*(-1.+-1.*std::sqrt(complex<double>(-3.)))/2.*std::cbrt(4.*(2.*std::pow(frp[2],3.)-9.*frp[3]*frp[2]*frp[1]+27.*frp[4]*std::pow(frp[1],2.)+27.*std::pow(frp[3],2.)*frp[0]-72.*frp[4]*frp[2]*frp[0]-std::sqrt(complex<double>(
					std::pow((2.*std::pow(frp[2],3.)-9.*frp[3]*frp[2]*frp[1]+27.*frp[4]*std::pow(frp[1],2.)+27.*std::pow(frp[3],2.)*frp[0]-72.*frp[4]*frp[2]*frp[0]),2.)-4.*std::pow((std::pow(frp[2],2.)-3.*frp[3]*frp[1]+12.*frp[4]*frp[0]),3.)
				))).real())
				
				)
			))+t*sgn((sgn(-std::pow(frp[3],3.)+4.*frp[4]*frp[3]*frp[2]-8.*std::pow(frp[4],2.)*frp[1])-1./2.)*(sgn(std::max(std::pow((2.*std::pow(frp[2],3.)-9.*frp[3]*frp[2]*frp[1]+27.*frp[4]*std::pow(frp[1],2.)+27.*std::pow(frp[3],2.)*frp[0]-72.*frp[4]*frp[2]*frp[0]),2.)-4.*std::pow((std::pow(frp[2],2.)-3.*frp[3]*frp[1]+12.*frp[4]*frp[0]),3.),std::min(3.*std::pow(frp[3],2.)-8.*frp[4]*frp[2],3.*std::pow(frp[3],4.)+16.*std::pow(frp[4],2.)*std::pow(frp[2],2.)+16.*std::pow(frp[4],2.)*frp[3]*frp[1]-16.*frp[4]*std::pow(frp[3],2.)*frp[2]-64.*std::pow(frp[4],3.)*frp[0])))-1./2.))*std::sqrt(
				3.*(3.*std::pow(frp[3],2.)-8.*frp[4]*frp[2]+
				
				2.*frp[4]*(-1.+-1.*std::sqrt(complex<double>(-3.)))/2.*std::cbrt(4.*(2.*std::pow(frp[2],3.)-9.*frp[3]*frp[2]*frp[1]+27.*frp[4]*std::pow(frp[1],2.)+27.*std::pow(frp[3],2.)*frp[0]-72.*frp[4]*frp[2]*frp[0]+std::sqrt(complex<double>(
					std::pow((2.*std::pow(frp[2],3.)-9.*frp[3]*frp[2]*frp[1]+27.*frp[4]*std::pow(frp[1],2.)+27.*std::pow(frp[3],2.)*frp[0]-72.*frp[4]*frp[2]*frp[0]),2.)-4.*std::pow((std::pow(frp[2],2.)-3.*frp[3]*frp[1]+12.*frp[4]*frp[0]),3.)
				))).real())
				
				+
				
				2.*frp[4]*(-1.+std::sqrt(complex<double>(-3.)))/2.*std::cbrt(4.*(2.*std::pow(frp[2],3.)-9.*frp[3]*frp[2]*frp[1]+27.*frp[4]*std::pow(frp[1],2.)+27.*std::pow(frp[3],2.)*frp[0]-72.*frp[4]*frp[2]*frp[0]-std::sqrt(complex<double>(
					std::pow((2.*std::pow(frp[2],3.)-9.*frp[3]*frp[2]*frp[1]+27.*frp[4]*std::pow(frp[1],2.)+27.*std::pow(frp[3],2.)*frp[0]-72.*frp[4]*frp[2]*frp[0]),2.)-4.*std::pow((std::pow(frp[2],2.)-3.*frp[3]*frp[1]+12.*frp[4]*frp[0]),3.)
				))).real())
				
				)
			))/(12.*frp[4]);
			//	end my favourite, most ridiculous line of code ever
			std::cout << "s=" << s << ", t=" << t << ":\t" << quart.real() << "\t+ " << quart.imag() << "i" << std::endl;	//	output roots to user
			if(quart.real() >= 0) if(quart.imag() == 0 || (s == 1 && t == 1 && tpx == -1))	//	if we have a real, positive root, or we are at the last root
				tpx = quart.real();	//	we need to keep the real part

		}	//	end root finder

		std::cout << "Using root with real component " << tpx << std::endl;	//	tell user we have a root

		for(int i = 0; i < 6; i++)	//	for each parameter
			myfnlp->SetParameter(i, (*it).second.second->Parameter(i));	//	get parameter for the efficiency function

		myfnlp->SetParameter(6, 1);	//	assume we are not normalised
		maxp = myfnlp->Eval(TMath::Exp(tpx));	//	determine turning point
		myfnlp->SetParameter(6, maxp);	//	normalise
		std::cout << channels[(*it).first] << ":" << std::endl << "\tEnergy Calibration:" << std::endl;	//	inform user of channel and calibration type

		for(int j = 0; j < (*it).second.first->NTotalParameters(); j++)	//	for each parameter
			std::cout << "\t" << j << ":\t" << (*it).second.first->Parameter(j) << ", " << (*it).second.first->Error(j);	//	inform user

		std::cout << std::endl << "\tEfficiency Calibration:" << std::endl;	//	inform user of calibration type

		for(int j = 0; j < (*it).second.second->NTotalParameters(); j++)	//	for each parameter
			std::cout << "\t" << j << ":\t" << (*it).second.second->Parameter(j) << ", " << (*it).second.second->Error(j);	//	inform user

		std::cout << "\t6:\t" << maxp << std::endl;	//	inform user of maximum of efficiency curve
		ge = new TGraphErrors();	//	generate blank efficiency plot

		for(std::map<std::pair<std::pair<double,double>,std::pair<double,double>>,std::pair<std::pair<string,int>,std::pair<string,int>>>::const_iterator jt = calib[(*it).first].begin(); jt != calib[(*it).first].end(); ++jt) {	//	for each calibration point

			gpos = ge->GetN();	//	get number of entries made so far
			scv = (*jt).second.first.first;	//	fetch string for calibration position
			sci = (*jt).second.second.first;	//	fetch string for integral amount
			dotpos = scv.find(".") + 1;	//	find position of decimal in mean
			declen = scv.substr(dotpos).length();	//	find precision of mean
			sce = (*jt).second.first.second / (dotpos > 0 && declen > 0 ? TMath::Power(10,declen) : 1);	//	calculate actual mean error from given integer
			dotpos = sci.find(".") + 1;	//	find position of decimal in integral
			declen = sci.substr(dotpos).length();	//	find precision of integral
			scd = (*jt).second.second.second / (dotpos > 0 && declen > 0 ? TMath::Power(10,declen) : 1);	//	calculate actual integral error from given integer
			ge->SetPoint(gpos, std::stod(scv), ((*jt).first.second.first / std::stod(sci)) / myfnlp->GetParameter(6));	//	plot efficiency point
			ge->SetPointError(gpos, sce, TMath::Sqrt(TMath::Power((*jt).first.second.second / std::stod(sci), 2) + TMath::Power((*jt).first.second.first * scd / TMath::Power(std::stod(sci), 2), 2)) / myfnlp->GetParameter(6));	//	plot efficiency point error

		}	//	end loop over all calibration points

		ge->SetMinimum(0);	//	ensure plot starts at sensible point
		ge->SetMaximum(1.01);	//	ensure plot covers curve area
		ge->GetXaxis()->SetLimits(0, 2000);	//	set limits based on reasonable accuracy of plot
		ge->Draw("A*");	//	draw efficiency plot
		(new TGraph(myfnlp))->Draw("same");	//	add efficiency curve
		hrv.push_back({ge,(*it).second.second});	//	add efficiency fit for user

	}

	return hrv;	//	return fits to user for further manipulation

}	//	end calibrator()
/*	This is a useful list of commands to copy and paste for making some plots:
hrv = calibrator()
hrv[0].first->SetLineColor(kAzure-1)
hrv[0].first->SetMarkerColor(kAzure-1)
hrv[2].first->SetLineColor(kRed+1)
hrv[2].first->SetMarkerColor(kRed+1)
hrv[4].first->SetLineColor(kSpring-1)
hrv[4].first->SetMarkerColor(kSpring-1)
hrv[0].first->SetMinimum(0)
hrv[2].first->SetMinimum(0)
hrv[4].first->SetMinimum(0)
hrv[0].first->Draw("A*")
hrv[2].first->Draw("*same")
hrv[4].first->Draw("*same")
TF1* p = (TF1*)gROOT->FindObject("pol1")
p->SetRange(0,200000)
p->SetFitResult(*(hrv[0].second.Get()))
g = (new TGraph(p))
g->SetLineColor(kAzure)
g->Draw("same")
p->SetFitResult(*(hrv[2].second.Get()))
g = (new TGraph(p))
g->SetLineColor(kRed)
g->Draw("same")
p->SetFitResult(*(hrv[4].second.Get()))
g = (new TGraph(p))
g->SetLineColor(kGreen)
g->Draw("same")

hrv[6].first->SetLineColor(kAzure-1)
hrv[7].first->SetLineColor(kRed+1)
hrv[8].first->SetLineColor(kSpring-1)
hrv[6].first->SetMarkerColor(kAzure-1)
hrv[7].first->SetMarkerColor(kRed+1)
hrv[8].first->SetMarkerColor(kSpring-1)
hrv[6].first->Draw("A*")
hrv[7].first->Draw("*same")
hrv[8].first->Draw("*same")

myfnl->SetFitResult(*(hrv[8].second.Get()));
for(int i = 0; i < 6; i++) myfnlp->SetParameter(i,myfnl->GetParameter(i));
myfnlp->SetParameter(6,1)
tpx = 5.08371
maxp = myfnlp->Eval(TMath::Exp(tpx));
myfnlp->SetParameter(6, maxp);
g = (new TGraph(myfnlp))
g->SetLineColor(kRed)
g->Draw("same")
*/