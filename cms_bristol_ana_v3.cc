//For changes see ChangeLog

//c++ includes
#include <iostream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <time.h>
using namespace std;

//ROOT includes
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChainElement.h"
#include "TGraphAsymmErrors.h"
#include "TROOT.h"
#include "TMath.h"

// RooFit includes, needed by m3 fit
#include "TRandom.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooConstVar.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#define SHOW(a) std::cout << #a << ": " << (a) << std::endl
using namespace RooFit;

// defines ana class, including branches and leaves for tree
#include "cms_bristol_ana_v3.hh"

void ana::SetInputFile(const char* fname) {
	// check the cut values before reading input
	static bool first_time = true;
	if (first_time)
		PrintCuts();

	nfile += chain->Add(fname);
	if (GetTrigger())
		chain2->Add(fname);
	first_time = false;
}

/** Output file name: <first name>_<second name>.root/log
 * User set in script.C the front bit of output file name ("first name"), the hist and text
 * output files will be set automatically according to type of input files:
 *   data:  data
 *   MC:    mc_[ttbar,wjets, bce1,2,3, enri1,2,3] (individual sample) or mc_mixture (mixture)
 */
void ana::SetOutputFirstName(const string name) {

	cout << "Total number in chain: " << chain->GetEntries() << " events\n\n";

	string secondname = "noname";

	if (IsData()) {

		secondname = "data";

	} else {//FIXME: this can ba A LOT shorter
		mc_names.resize(NUMBER_OF_MC_TYPES);
		for (ushort x = kdata + 1; x < (short) NUMBER_OF_MC_TYPES; ++x) {
			mc_names[x] = mcfiles[x];
		}
		// TL 22-1-09
		//		mc_names.push_back("ttbar");
		//		mc_names.push_back("ttjet");
		//		mc_names.push_back("wjet");
		//		mc_names.push_back("wenu");//old
		//		mc_names.push_back("zjet");
		//		mc_names.push_back("zee");//old
		//		mc_names.push_back("enri1");
		//		mc_names.push_back("enri2");
		//		mc_names.push_back("enri3");
		//		mc_names.push_back("bce1");
		//		mc_names.push_back("bce2");
		//		mc_names.push_back("bce3");
		//		mc_names.push_back("vqq");
		//		mc_names.push_back("tW");
		//		mc_names.push_back("tchan");
		//		mc_names.push_back("schan");
		//		//FIXME: 04th step to new sample - add to vector
		//		mc_names.push_back("Zprime_M500GeV_W5GeV");
		//		mc_names.push_back("Zprime_M500GeV_W50GeV");
		//		mc_names.push_back("Zprime_M750GeV_W7500MeV");
		//		mc_names.push_back("Zprime_M1TeV_W10GeV");
		//		mc_names.push_back("Zprime_M1TeV_W100GeV");
		//		mc_names.push_back("Zprime_M1250GeV_W12500MeV");
		//		mc_names.push_back("Zprime_M1500GeV_W15GeV");
		//		mc_names.push_back("Zprime_M1500GeV_W150GeV");
		//		mc_names.push_back("Zprime_M2TeV_W20GeV");
		//		mc_names.push_back("Zprime_M2TeV_W200GeV");
		//		mc_names.push_back("Zprime_M3TeV_W30GeV");
		//		mc_names.push_back("Zprime_M3TeV_W300GeV");
		//		mc_names.push_back("Zprime_M4TeV_W40GeV");
		//		mc_names.push_back("Zprime_M4TeV_W400GeV");

		vector<int> nfiles(NUMBER_OF_MC_TYPES);
		vector<long> nevents(NUMBER_OF_MC_TYPES);

		cout << "List of input files:" << endl;

		// loop over all files in chain
		TObjArray *fileElements = chain->GetListOfFiles();
		TIter next(fileElements);
		TChainElement *chEl = 0;

		//		map<string, bool> mc_seen;
		vector<bool> mcseen;
		mcseen.resize(NUMBER_OF_MC_TYPES);

		while ((chEl = (TChainElement*) next())) {

			const string fname(chEl->GetTitle());
			cout << fname << ":  " << chEl->GetEntries() << " events";

			// try to find mc type i
			size_t i = kdata + 1;
			while (i < (short) NUMBER_OF_MC_TYPES) {
				const string mc = mc_names.at(i);
				if (fname.find(mc) != string::npos) { //match
					cout << "  " << mc << endl;
					nfiles[i]++;
					nevents[i] += chEl->GetEntries();
					// if first time seeing this mc, then add this mc to map.
					if (!mcseen[i])
						mcseen[i] = true;
					break;
				}
				++i;
			}
		}
		//store the event count in private variable with key=mc
		for (size_t i = kdata + 1; i < (short) NUMBER_OF_MC_TYPES; ++i) {
			nInitialEventMC_[i] = nevents.at(i);
			//			nInitialEventMC[mc_names.at(i)] = nevents.at(i);
		}

		// Set event weight map, and MC flag
		DefineCrossSection();
		SetWeights();
		//		SetMCFlag();

		// print out
		cout << "\n--------------------------------------------------------------------------------------------\n";
		cout << setw(10) << "mc" << setw(14) << "N(files)" << setw(18) << "N(events)" << setw(12) << "weight" << setw(15)
				<< "xsec" << setw(14) << "Nexp (L=" << intlumi << "/pb)";
		cout << "\n--------------------------------------------------------------------------------------------\n";
		for (size_t i = 0; i < mc_names.size(); ++i) {

			cout << setw(10) << mc_names.at(i) << setw(8) << nfiles.at(i) << " files" << setw(11) << nevents.at(i) << " events";
			cout << setw(11);

			if (nfiles.at(i) == 0)
				cout << "-";
			else
				cout << setprecision(4) << GetWeight(i);

			// print out cross-sections
			cout << setw(14) << cross_section_[i] << " pb";
			cout << setw(14) << cross_section_[i] * intlumi << endl;
			cout << setprecision(6);
		}
		cout << "--------------------------------------------------------------------------------------------\n";

		cout << "\n MC: Listing types of MC found in input\n";
		//FIXME: isMCPresent
		ushort NofMCseen = 0;
		ushort seenID = 0;
		for (ushort i = 0; i < (short) NUMBER_OF_MC_TYPES; ++i) {
			if (mcseen[i]) {
				cout << "  " << mc_names[i] << endl;
				++NofMCseen;
				seenID = i;
			}

		}
		//		for (map<string, bool>::iterator i = mc_seen.begin(); i != mc_seen.end(); ++i) {
		//			cout << "  " << i->first << endl;
		//		}

		if (NofMCseen == 1)
			secondname = "mc_" + mc_names[seenID]; //one MC only
		else if (NofMCseen > 1)
			secondname = "mc_mixture"; //there are more than 1 kind of MC in input
	}// data or MC

	const string outname1 = name + "_" + secondname + ".root";
	const string outname2 = name + "_" + secondname + ".log";
	cout << "\n Analysis outputs are: " << outname1 << endl << "                       " << outname2 << endl << endl;
	SetOutputHistFile(outname1);
	SetOutputTextFile(outname2);

}//end SetOutputFirstName

void ana::SetOutputHistFile(const string name, const string mode) {

	histf = new TFile(name.c_str(), mode.c_str(), "test file");
	outputHistFileName = name;
}

void ana::SetOutputTextFile(const string name, const string mode) {

	outfile = fopen(name.c_str(), mode.c_str());
	outputTextFileName = name;
}

void ana::CheckTrigger(bool val) {
	checkTrig = val;
}

void ana::CheckTrigger(bool val, const string hlt) {
	checkTrig = val;
	HLTBit = hlt;
}

void ana::SetLimit(int val) {
	if (val > 0) {
		numlimit = val;
	} else {
		numlimit = -1;
	}
}

//set kinematics cuts
void ana::SetEleETcut(float cut) {
	ELE_ETCUT = cut;
	nCutSetInScript++;
}
void ana::SetMuonPTcut(float cut) {
	MU_PTCUT = cut;
	nCutSetInScript++;
}
void ana::SetJetETcut(float cut) {
	JET_PTCUT = cut;
	nCutSetInScript++;
}
void ana::SetMETcut(float cut) {
	METCUT = cut;
	nCutSetInScript++;
}
void ana::SetHTcut(float cut) {
	HTCUT = cut;
	nCutSetInScript++;
}

void ana::PrintCuts() const {

	if (IsData())
		cout << " This run uses   REAL   Data\n\n";
	else
		cout << " This run uses   MONTE CARLO   Data\n\n";
	cout << "***********************************************\n\n";
	cout << "  LHC c.o.m Energy :  " << m_LHCEnergyInTeV << " TeV\n" << endl;
	cout << "***********************************************\n" << endl;
	if (checkTrig)
		cout << "  Trigger required :   " << HLTBit << endl;
	else
		cout << "  Trigger not required" << endl;
	cout << "\n***********************************************" << endl;
	if (nCutSetInScript < 4) {
		cout << "\n  WARNING! WARNING! YOU HAVE NOT SET ALL 4 CUTS!!!" << endl;
	} else {
		cout << "\n Okay, all 4 cuts are set." << endl;
	}
	cout << "\n This run will use the following cut values:\n" << endl;
	cout << " ELECTRON et cut =  " << ELE_ETCUT << "  GeV" << endl;
	cout << "     MUON pt cut =  " << MU_PTCUT << "  GeV" << endl;
	cout << "      JET pt cut =  " << JET_PTCUT << "  GeV" << endl;
	cout << "      MET    cut =  " << METCUT << "  GeV" << endl;
	cout << " ELECTRON ID requirement  =  " << printEleID() << endl;
	cout << "\n***********************************************" << endl;
	cout << "\n  HT cut in AES:  " << AES_HT_cut << " GeV";
	cout << "\n MET cut in AES:  " << AES_MET_cut << " GeV" << endl;
	if (AES_useSimpleZveto) {
		cout << "\n  In AES, exclude events with 2 reco-electrons" << endl;
	}
	cout << "\n***********************************************" << endl;
	cout << "\n Jet collection:  " << m_jetAlgo << endl;
	cout << "\n MET collection:  " << m_metAlgo << endl;
	cout << "\n***********************************************" << endl;
	cout << "\n use_old_Z_veto : " << use_old_Z_veto << endl;
	cout << "\n***********************************************" << endl;
	if (!m_applyMETcut && m_rejectEndcapEle) {
		cout << "\n  Run with Plan B:  MET is not used, add Eta(e) <1.442 as last cut" << endl;
		cout << "\n***********************************************" << endl;
	}
	if (m_applyMETcut && m_rejectEndcapEle) {
		cout << "\n  Run with Plan C:  Both MET and Eta(e) <1.442 cuts are used" << endl;
		cout << "\n***********************************************" << endl;
	}

}//PrintCuts

//TODO: replace with enum
long ana::GetNinit(const ushort mc) const {
	return nInitialEventMC_.at(mc);
	// if found this MC
	//	map<string, long>::const_iterator iter = nInitialEventMC.find(mc);
	//	if (iter != nInitialEventMC.end()) // this MC found in input
	//		return iter->second;
	//	else
	//		//return 1;
	//		return 0;
}

bool ana::hasSampleMC(const ushort mc) const {
	if (mc == kdata)//data is alway there
		return true;
	if (mc == kttbar || mc == kttjet)
		return (GetNinit(kttbar) > 0 || GetNinit(kttjet) > 0);
	else if (mc == kQCD)
		return GetNinit(kQCD) > 0 || GetNinit(kenri1) > 0 || GetNinit(kenri2) > 0 || GetNinit(kenri3) > 0 || GetNinit(kbce1) > 0
				|| GetNinit(kbce2) > 0 || GetNinit(kbce3) > 0;
	else if (mc == ksingleTop)
		return GetNinit(ksingleTop) > 0 || GetNinit(ktW) > 0 || GetNinit(ktchan) > 0 || GetNinit(kschan) > 0;
	else
		return GetNinit(mc) > 0;
}

bool ana::isMCType(const ushort mc) const {
	if (mc == kttbar || mc == kttjet)
		return fastmctype_ == kttbar || fastmctype_ == kttjet;
	else if (mc == kQCD)
		return fastmctype_ == kQCD || fastmctype_ == kenri1 || fastmctype_ == kenri2 || fastmctype_ == kenri3 || fastmctype_
				== kbce1 || fastmctype_ == kbce2 || fastmctype_ == kbce3;
	else if (mc == ksingleTop)
		return fastmctype_ == ksingleTop || fastmctype_ == ktW || fastmctype_ == ktchan || fastmctype_ == kschan;
	else
		return mc == fastmctype_;
}

void ana::DefineCrossSection() {
	cross_section_.resize(NUMBER_OF_MC_TYPES);
	if (m_LHCEnergyInTeV == 10) { //10TeV = default

		cout << "\nNOTE: 10 TeV cross section not to date!!!" << endl;
		cout << "\n check https://twiki.cern.ch/twiki/bin/viewauth/CMS/CrossSections_3XSeries\n" << endl;

		cross_section_[kttbar] = 414.; //xs 414 (NLO+NLL at 10TeV)   241.7 pb (LO)
		cross_section_[kttjet] = 414.;

		cross_section_[kWjets] = (40e3 * 1.14); //xs 40 nb  (K-factor=1.14)
		cross_section_[kZjets] = (3.7e3 * 1.14); //xs 3.7 nb (K-factor=1.14)

		//		cross_section["wenu"] = (11840 * 0.738 * 1.14); //xs 11840 pb (eff=0.738) (K-factor=1.14) // TO-BE-CHECKED
		//		cross_section["zee"] = (1944 * 1.14); //xs 1944 pb (K-factor=1.14)             // TO-BE-CHECKED

		cross_section_[kenri1] = 0.40e9 * 0.008; //xs 0.4 mb (*filter efficiency=0.008)
		cross_section_[kenri2] = 0.10e9 * 0.047; //xs 0.1 mb
		cross_section_[kenri3] = 1.9e6 * 0.15; //xs 1.9e-3 mb

		cross_section_[kbce1] = 0.40e9 * 4.8e-4; //xs 0.4 mb
		cross_section_[kbce2] = 0.10e9 * 2.4e-3; //xs 0.1 mb
		cross_section_[kbce3] = 1.9e6 * 0.012; //xs 1.9e-3 mb

		cross_section_[kVQQ] = 290.; //xs 290 pb (W->lep)(LO)
		cross_section_[ktW] = 29.; //xs  29 pb (NLO) inclusive t,W decay
		cross_section_[ktchan] = 130.0 * 0.324; //xs 130 pb (NLO) * BR(t->blnu)
		cross_section_[kschan] = 5.0 * 0.324; //xs   5 pb (NLO) * BR(t->blnu)
	} else if (m_LHCEnergyInTeV == 7) { //7TeV

		cross_section_[kttbar] = 165.; //was 187 pb
		cross_section_[kttjet] = 165.;

		cross_section_[kWjets] = 28000.; //was 24e3 * 1.14
		cross_section_[kZjets] = 2800.; //was 2.3e3 * 1.14

		//		cross_section["wenu"] = 7899 * 0.779 * 1.14; //xs 7899 pb (pythia LO) (eff=0.779) (K-fac=1.14)
		//		cross_section["zee"] = 1300 * 1.14; //xs 1300 pb (pythia LO) (K-fac=1.14)

		cross_section_[kenri1] = 0.2355e9 * 0.0073; //xs 0.2355 mb (filter efficiency=0.0073)
		cross_section_[kenri2] = 0.0593e9 * 0.059; //xs 0.0593 mb
		cross_section_[kenri3] = 0.906e6 * 0.148; //xs 0.906e-3 mb

		cross_section_[kbce1] = 0.2355e9 * 0.00046; //xs 0.2355 mb (filter efficiency=0.00046)
		cross_section_[kbce2] = 0.0593e9 * 0.00234; //xs 0.0593 mb
		cross_section_[kbce3] = 0.906e6 * 0.0104; //xs 0.906e-3 mb

		cross_section_[ktW] = 10.6; //xs  11 pb (NLO MCFM) inclusive t,W decay
		cross_section_[ktchan] = 63. * 0.324; //xs  63 pb (NLO MCFM) * 0.324 (Br(t->blnu)) = 20.412
		cross_section_[kschan] = 4.6 * 0.324; //4.6 pb x 0.324 (4.6pb is NNNLO) = 1.4904
		//FIXME: 06th step to new sample - set cross section
		//change to Z x-section (simple approach): sigma(Z) ~ 1/M_Z^2 => replace masses
		//dependence on width is neglected
		double zcross_normed = cross_section_[kZjets] * ZMASS * ZMASS;
		cross_section_[Zprime_M500GeV_W5GeV] = zcross_normed / 500 / 500;// or =50
		cross_section_[Zprime_M500GeV_W50GeV] = zcross_normed / 500 / 500;
		cross_section_[Zprime_M750GeV_W7500MeV] = zcross_normed / 750 / 750;//or = 15
		cross_section_[Zprime_M1TeV_W10GeV] = zcross_normed / 1000 / 1000;//or = 7
		cross_section_[Zprime_M1TeV_W100GeV] = zcross_normed / 1000 / 1000;
		cross_section_[Zprime_M1250GeV_W12500MeV] = zcross_normed / 1250 / 1250;
		cross_section_[Zprime_M1500GeV_W15GeV] = zcross_normed / 1500 / 1500;//or = 3
		cross_section_[Zprime_M1500GeV_W150GeV] = zcross_normed / 1500 / 1500;
		cross_section_[Zprime_M2TeV_W20GeV] = zcross_normed / 2000 / 2000;
		cross_section_[Zprime_M2TeV_W200GeV] = zcross_normed / 2000 / 2000;
		cross_section_[Zprime_M3TeV_W30GeV] = zcross_normed / 3000 / 3000;
		cross_section_[Zprime_M3TeV_W300GeV] = zcross_normed / 3000 / 3000;
		cross_section_[Zprime_M4TeV_W40GeV] = zcross_normed / 4000 / 4000;
		cross_section_[Zprime_M4TeV_W400GeV] = zcross_normed / 4000 / 4000;

	} else {
		if (!IsData())
			cout << "WARNING: Cross section values are not defined!" << endl;
	}

}//end DefineCrossSection
//-------------------------------------------------------------------------------------------

//TODO: change map to enum and vector -- here
double ana::GetCrossSection(const ushort mctype) const {
	if (mctype >= NUMBER_OF_MC_TYPES) {
		cout << "WARNING: Trying to access MC type out of range: " << mctype << " out of " << NUMBER_OF_MC_TYPES << endl;
	}
	double xsection = cross_section_.at(mctype);
	if (xsection == 0)
		cout << "WARNING: cross section not defined for MC Type No. " << mctype << " ('" << mclabel[mctype] << "')" << endl;
	return xsection;

}
//-------------------------------------------------------------------------------------------
// Declare the event weights
void ana::SetWeights() { //only if run on MC
	fastWeight_.resize(NUMBER_OF_MC_TYPES);
	fastWeight_[kdata] = 1.0;
	//	weightMap["data"] = 1.0;

	for (ushort i = 1; i < (short) NUMBER_OF_MC_TYPES; ++i) {
		long n = GetNinit(i);
		if (n > 0)
			fastWeight_[i] = GetCrossSection(i) * intlumi / n;
		//			weightMap[mc_names.at(i)] = GetCrossSection(i) * intlumi / n;
		else
			fastWeight_[i] = 0;
		//			weightMap[mc_names.at(i)] = 0;
	}

	//OctX for SD
	if (m_runOnSD) {
		if (m_LHCEnergyInTeV == 10) {
			cout << "\nSummer09 10 TeV, OctX SD skim efficiencies x prescale factor:" << endl;
			cout << "(note: nexp = Ninit * w / skim eff)" << endl;
			// rescale weight by filter efficiency (N_SD/N_ori), and prescale factor
			//                                  N_SD / N_ori     * pres
			const double skimEff_ttbar = 3270 / 529750. * 100;
			const double skimEff_wenu = 87514 / 2142960. * 20;
			const double skimEff_zee = 15092 / 2682355. * 100;
			const double skimEff_enri1 = 6242601 / 33638282.;
			const double skimEff_enri2 = 10792435 / 38360886.;
			const double skimEff_enri3 = 2538711 / 5729547.;
			const double skimEff_bce1 = 393019 / 2383833.;
			const double skimEff_bce2 = 769808 / 2035108.;
			const double skimEff_bce3 = 641522 / 1038080.;
			fastWeight_[kttbar] = fastWeight_[kttbar] * skimEff_ttbar;
			fastWeight_[kWjets] = fastWeight_[kWjets] * skimEff_wenu;
			fastWeight_[kZjets] = fastWeight_[kZjets] * skimEff_zee;
			fastWeight_[kenri1] = fastWeight_[kenri1] * skimEff_enri1;
			fastWeight_[kenri2] = fastWeight_[kenri2] * skimEff_enri2;
			fastWeight_[kenri3] = fastWeight_[kenri3] * skimEff_enri3;
			fastWeight_[kbce1] = fastWeight_[kbce1] * skimEff_bce1;
			fastWeight_[kbce2] = fastWeight_[kbce2] * skimEff_bce1;
			fastWeight_[kbce3] = fastWeight_[kbce3] * skimEff_bce1;
			//			weightMap["ttbar"] = weightMap["ttbar"] * skimEff_ttbar;
			//			weightMap["wenu"] = weightMap["wenu"] * skimEff_wenu;
			//			weightMap["zee"] = weightMap["zee"] * skimEff_zee;
			//			weightMap["enri1"] = weightMap["enri1"] * skimEff_enri1;
			//			weightMap["enri2"] = weightMap["enri2"] * skimEff_enri2;
			//			weightMap["enri3"] = weightMap["enri3"] * skimEff_enri3;
			//			weightMap["bce1"] = weightMap["bce1"] * skimEff_bce1;
			//			weightMap["bce2"] = weightMap["bce2"] * skimEff_bce2;
			//			weightMap["bce3"] = weightMap["bce3"] * skimEff_bce3;

			cout << "  skim eff   ttbar    " << skimEff_ttbar << endl;
			cout << "  skim eff   wenu     " << skimEff_wenu << endl;
			cout << "  skim eff   zee      " << skimEff_zee << endl;
			cout << "  skim eff   enri1    " << skimEff_enri1 << endl;
			cout << "  skim eff   enri2    " << skimEff_enri2 << endl;
			cout << "  skim eff   enri3    " << skimEff_enri3 << endl;
			cout << "  skim eff   bce1     " << skimEff_bce1 << endl;
			cout << "  skim eff   bce2     " << skimEff_bce2 << endl;
			cout << "  skim eff   bce3     " << skimEff_bce3 << endl;
			cout << endl;
		} else if (m_LHCEnergyInTeV == 7) {
			cout << "\nSummer09 7 TeV, OctX SD skim efficiencies x prescale factor:" << endl;
			cout << "(note: nexp = Ninit * w / skim eff)" << endl;
			// rescale weight by filter efficiency (N_SD/N_ori), and prescale factor
			//                                 N_SD / N_ori     * pres
			const double skimEff_ttbar = 3841 / 626610. * 100;
			const double skimEff_wenu = 85505 / 2078361. * 20;
			const double skimEff_zee = 14683 / 2538855. * 100;
			const double skimEff_enri1 = 6169999 / 33505929.;
			const double skimEff_enri2 = 9054696 / 32168675.;
			const double skimEff_enri3 = 2463429 / 5551386.;
			const double skimEff_bce1 = 432380 / 2752942.;
			const double skimEff_bce2 = 840100 / 2261916.;
			const double skimEff_bce3 = 682720 / 1097829.;

			fastWeight_[kttbar] = fastWeight_[kttbar] * skimEff_ttbar;
			fastWeight_[kWjets] = fastWeight_[kWjets] * skimEff_wenu;//FIXME wenu no included
			fastWeight_[kZjets] = fastWeight_[kZjets] * skimEff_zee;
			fastWeight_[kenri1] = fastWeight_[kenri1] * skimEff_enri1;
			fastWeight_[kenri2] = fastWeight_[kenri2] * skimEff_enri2;
			fastWeight_[kenri3] = fastWeight_[kenri3] * skimEff_enri3;
			fastWeight_[kbce1] = fastWeight_[kbce1] * skimEff_bce1;
			fastWeight_[kbce2] = fastWeight_[kbce2] * skimEff_bce1;
			fastWeight_[kbce3] = fastWeight_[kbce3] * skimEff_bce1;
			//			weightMap["ttbar"] = weightMap["ttbar"] * skimEff_ttbar;
			//			weightMap["wenu"] = weightMap["wenu"] * skimEff_wenu;
			//			weightMap["zee"] = weightMap["zee"] * skimEff_zee;
			//			weightMap["enri1"] = weightMap["enri1"] * skimEff_enri1;
			//			weightMap["enri2"] = weightMap["enri2"] * skimEff_enri2;
			//			weightMap["enri3"] = weightMap["enri3"] * skimEff_enri3;
			//			weightMap["bce1"] = weightMap["bce1"] * skimEff_bce1;
			//			weightMap["bce2"] = weightMap["bce2"] * skimEff_bce2;
			//			weightMap["bce3"] = weightMap["bce3"] * skimEff_bce3;

			cout << "  skim eff   ttbar    " << skimEff_ttbar << endl;
			cout << "  skim eff   wenu     " << skimEff_wenu << endl;
			cout << "  skim eff   zee      " << skimEff_zee << endl;
			cout << "  skim eff   enri1    " << skimEff_enri1 << endl;
			cout << "  skim eff   enri2    " << skimEff_enri2 << endl;
			cout << "  skim eff   enri3    " << skimEff_enri3 << endl;
			cout << "  skim eff   bce1     " << skimEff_bce1 << endl;
			cout << "  skim eff   bce2     " << skimEff_bce2 << endl;
			cout << "  skim eff   bce3     " << skimEff_bce3 << endl;
			cout << endl;
		}
	}
	// for skim
	if (run_on_octX_skim) {
		///                                          N_skim / N_ori    * pres
		fastWeight_[kbce1] = fastWeight_[kbce1] * 144791 / 2383833;
		//		weightMap["bce1"] = weightMap["bce1"] * 144791 / 2383833;
	}

	// for summer09 7TeV madgraph HLTskim
	if (m_runOnMyHLTskim) {
		cout << "\nSummer09 7 TeV, Madgraph my HLT skim efficiency:" << endl;
		cout << "(note: nexp = Ninit * w / skim eff)" << endl;
		///                         N_skim / N_ori    * pres
		const double skimEff_ttj = 610804 / 983964.;
		const double skimEff_wj = 2081537 / 8109289.;
		const double skimEff_zj = 381190 / 1068735.;
		fastWeight_[kttjet] = fastWeight_[kttjet] * skimEff_ttj;
		fastWeight_[kWjets] = fastWeight_[kWjets] * skimEff_wj;
		fastWeight_[kZjets] = fastWeight_[kZjets] * skimEff_zj;
		//		weightMap["ttjet"] = weightMap["ttjet"] * skimEff_ttj;
		//		weightMap["wjet"] = weightMap["wjet"] * skimEff_wj;
		//		weightMap["zjet"] = weightMap["zjet"] * skimEff_zj;

		cout << "  skim eff    ttjet      " << skimEff_ttj << endl;
		cout << "  skim eff    wjet       " << skimEff_wj << endl;
		cout << "  skim eff    zjet       " << skimEff_zj << endl;
		cout << endl;
	}

}//End SetWeights

//-------------------------------------------------------------------------------------------
//double ana::GetWeight(string mc) const {
//	return weightMap.find(mc)->second;
//}

double ana::GetWeight(short mc) const {
	return fastWeight_.at(mc);
}
//-------------------------------------------------------------------------------------------

ana::ana() {

	// Print out code compilation time and current time
	cout << "\nCompilation Date/Time = " << __DATE__ << "  " << __TIME__ << endl;
	TDatime now;
	cout << "Current local ";
	now.Print();
	cout << endl;

	cout << "\n***********************************************************************";
	cout << "\n*                                                                     *";
	cout << "\n*        C M S     B R I S T O L     T O P     A N A L Y S I S        *";
	cout << "\n*                                                                     *";
	cout << "\n***********************************************************************";
	cout << endl << endl;

	chain = new TChain("configurableAnalysis/eventB");
	chain2 = new TChain("configurableAnalysis/eventV");

	// Initialize private variables
	datafile = true;
	keepgood = true;
	checkTrig = true;
	HLTBit = "HLT_Ele15_LW_L1R";
	numlimit = -1;
	this_weight = 1.0;
	m_nGoodJet = 0;
	m_QCDest_reliso_bin_width = 0.01;
	m_doValidation = false;
	m_plotRelisoNES = true;
	m_debug = false;
	m_ConversionStudies = false;
	m_jetAlgo = "Default";
	m_metAlgo = "calojet_mujes";
	m_LHCEnergyInTeV = 7.0; //Default is 7 TeV
	m_runOnSD = false;
	m_runOnMyHLTskim = false;
	m_useMisslayers = false;
	m_ntoy = 2;

	// Default values of kinematic cuts
	ELE_ETCUT = 30.0;
	MU_PTCUT = 20.0;
	JET_PTCUT = 30.0;
	METCUT = 30.0;
	HTCUT = 0.0;
	m_applyMETcut = true;
	m_rejectEndcapEle = true;
	nCutSetInScript = 0;
	AES_HT_cut = 200.0;
	AES_MET_cut = 15.0;
	AES_useSimpleZveto = true;
	m_eID = robustTight; //enum

	// integrated luminosity
	intlumi = 20.0; //pb-1
	useNewReliso = true;
	doSystematics = "";
	sysSample = "";
	m_studyPDFunc = false;
	ScientificNotation = true;

	// initialize private variables (MC study)
	//	isMCTypePresent_.resize(NUMBER_OF_MC_TYPES);
	//	for (ushort type = 0; type < (short) NUMBER_OF_MC_TYPES; ++type) {
	//		isMCTypePresent_[type] = false;
	//	}
	nInitialEventMC_.resize(NUMBER_OF_MC_TYPES);

	//	mc_sample_has_ttbar = false;
	//	mc_sample_has_Wjet = false;
	//	mc_sample_has_Zjet = false;
	//	mc_sample_has_QCD = false;
	//	hasSampleMC(kenri1) = false;
	//	hasSampleMC(kenri2) = false;
	//	hasSampleMC(kenri3) = false;
	//	hasSampleMC(kbce1) = false;
	//	hasSampleMC(kbce2) = false;
	//	hasSampleMC(kbce3) = false;
	//	hasSampleMC(kVQQ) = false;
	//	mc_sample_has_singleTop = false;
	//	mc_sample_has_tW = false;
	//	mc_sample_has_tchan = false;
	//	mc_sample_has_schan = false;
	//	//FIXME: 07th step to new samples - set mc_sample_has*
	//	mc_sample_has_Zprime_M500GeV_W5GeV = false;
	//	mc_sample_has_Zprime_M500GeV_W50GeV = false;
	//	mc_sample_has_Zprime_M750GeV_W7500MeV = false;
	//	mc_sample_has_Zprime_M1TeV_W10GeV = false;
	//	mc_sample_has_Zprime_M1TeV_W100GeV = false;
	//	mc_sample_has_Zprime_M1250GeV_W12500MeV = false;
	//	mc_sample_has_Zprime_M1500GeV_W15GeV = false;
	//	mc_sample_has_Zprime_M1500GeV_W150GeV = false;
	//	mc_sample_has_Zprime_M2TeV_W20GeV = false;
	//	mc_sample_has_Zprime_M2TeV_W200GeV = false;
	//	mc_sample_has_Zprime_M3TeV_W30GeV = false;
	//	mc_sample_has_Zprime_M3TeV_W300GeV = false;
	//	mc_sample_has_Zprime_M4TeV_W40GeV = false;
	//	mc_sample_has_Zprime_M4TeV_W400GeV = false;

	//	isMCType(kttbar) = false;
	//	isMCType(kWjets) = false;
	//	isMCType(kZjets) = false;
	//	isMCType(kQCD) = false;
	//	isMCType(kenri1) = false;
	//	isMCType(kenri2) = false;
	//	isMCType(kenri3) = false;
	//	isMCType(kbce1) = false;
	//	isMCType(kbce2) = false;
	//	isMCType(kbce3) = false;
	//	isMCType(kVQQ) = false;
	//	isMCType(ksingleTop) = false;
	//	isMCType(ktW) = false;
	//	isMCType(ktchan) = false;
	//	isMCType(kschan) = false;
	//	//FIXME: 08th step to new samples - set is*
	//	isMCType(Zprime_M500GeV_W5GeV) = false;
	//	isMCType(Zprime_M500GeV_W50GeV) = false;
	//	isMCType(Zprime_M750GeV_W7500MeV) = false;
	//	isMCType(Zprime_M1TeV_W10GeV) = false;
	//	isMCType(Zprime_M1TeV_W100GeV) = false;
	//	isMCType(Zprime_M1250GeV_W12500MeV) = false;
	//	isMCType(Zprime_M1500GeV_W15GeV) = false;
	//	isMCType(Zprime_M1500GeV_W150GeV) = false;
	//	isMCType(Zprime_M2TeV_W20GeV) = false;
	//	isMCType(Zprime_M2TeV_W200GeV) = false;
	//	isMCType(Zprime_M3TeV_W30GeV) = false;
	//	isMCType(Zprime_M3TeV_W300GeV) = false;
	//	isMCType(Zprime_M4TeV_W40GeV) = false;
	//	isMCType(Zprime_M4TeV_W400GeV) = false;

	//856
	ConversionCounter = 0;
	for (int k = 0; k < 23; ++k) {
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 6; ++j) {
				ConversionArray[k][i][j] = 0;
			}
		}
	}
	// end 856

	mycounter = 0;
	SetNumberOfJetsUsedInReco(8);
}

ana::~ana() {
}

void ana::End() {
	cout << "End of EventLoop" << endl;
}

void ana::CheckAvailableJetMET() {

	vector<string> availableJET;
	vector<string> availableMET;

	// check if various jet/met is in the ntuple
	TObjArray *list = chain->GetListOfBranches();

	// -- Loop over all, and draw their variables into TCanvas c1
	for (int i = 0; i < list->GetEntries(); ++i) {
		string name = ((TBranch*) (*list)[i])->GetName();
		// jet
		if (name == "NjetsSC5")
			availableJET.push_back("SC5");
		else if (name == "NjetsSC7")
			availableJET.push_back("SC7");
		else if (name == "NjetsKT4")
			availableJET.push_back("KT4");
		else if (name == "NjetsKT6")
			availableJET.push_back("KT6");
		else if (name == "NPFJets")
			availableJET.push_back("PF(AK5)");
		else if (name == "NjetsJPTAK5")
			availableJET.push_back("JPT(AK5)");
		// met
		if (name == "NmetsSC5")
			availableMET.push_back("SC5");
		else if (name == "NmetsSC7")
			availableMET.push_back("SC7");
		else if (name == "NmetsKT4")
			availableMET.push_back("KT4");
		else if (name == "NmetsKT6")
			availableMET.push_back("KT6");
		else if (name == "NPFMets")
			availableMET.push_back("PF(AK5)");
		else if (name == "Ntcmets")
			availableMET.push_back("tcmet");
	}
	cout << " Available jet: ";
	for (unsigned int i = 0; i < availableJET.size(); ++i) {
		cout << " " << availableJET[i];
	}
	cout << endl;
	cout << " Available met: ";
	for (unsigned int i = 0; i < availableMET.size(); ++i) {
		cout << " " << availableMET[i];
	}
	cout << endl;
	availableJET.clear();
	availableMET.clear();

}// end CheckAvailableJetMET()
//------------------------------------------------------------------------------------


bool ana::EventLoop() {

	if (nfile == 0) {
		cout << "No input file found, stop." << endl;
		return false;
	}
	bookHistograms();
	CheckAvailableJetMET();

	Init(); //initialize branch

	// Get the number of events/entries in the file chain
	Long64_t nEvents = chain->GetEntries();
	Long64_t nEventsAvail = nEvents;
	if (nEvents == 0) {
		cout << "No input event found, stop." << endl;
		return false;
	}

	if (GetTrigger()) {
		chain->AddFriend(chain2);
	}

	// Read only selected branches to reduce cpu time
	ReadSelectedBranches();

	//ttbar decay code from PDG codes
	TH1F *all_mctype = new TH1F("all_mctype", "all_mctype (unweighted) ttbar decay modes", 10, 1, 11);

	//Electron Plots
	TH1F *h_hadOverEm_barrel = new TH1F("h_hadOverEm_barrel", "H/E - Barrel", 100, 0, 0.5);
	TH1F *h_EOverPIn_barrel = new TH1F("h_EOverPIn_barrel", "E/P (In) - Barrel", 100, 0, 10.0);
	TH1F *h_dEtaIn_barrel = new TH1F("h_dEtaIn_barrel", "Delta Eta (In) - Barrel", 100, -0.05, 0.05);
	TH1F *h_dPhiIn_barrel = new TH1F("h_dPhiIn_barrel", "Delta Phi (In) - Barrel", 100, -0.2, 0.2);
	TH1F *h_hadOverEm_endcap = new TH1F("h_hadOverEm_endcap", "H/E - Endcap", 100, 0, 0.5);
	TH1F *h_EOverPIn_endcap = new TH1F("h_EOverPIn_endcap", "E/P (In) - Endcap", 100, 0, 10.0);
	TH1F *h_dEtaIn_endcap = new TH1F("h_dEtaIn_endcap", "Delta Eta (In) - Endcap", 100, -0.05, 0.05);
	TH1F *h_dPhiIn_endcap = new TH1F("h_dPhiIn_endcap", "Delta Phi (In) - Endcap", 100, -0.2, 0.2);
	TH1F *h_cIso_barrel = new TH1F("h_cIso_barrel", "Cal Iso - Barrel", 20, 0, 10);
	TH1F *h_tIso_barrel = new TH1F("h_tIso_barrel", "Trk Iso - Barrel", 20, 0, 10);
	TH1F *h_cIso_endcap = new TH1F("h_cIso_endcap", "Cal Iso - Endcap", 20, 0, 10);
	TH1F *h_tIso_endcap = new TH1F("h_tIso_endcap", "Trk Iso - Endcap", 20, 0, 10);
	h_hadOverEm_barrel->Sumw2();
	h_EOverPIn_barrel ->Sumw2();
	h_dEtaIn_barrel ->Sumw2();
	h_dPhiIn_barrel ->Sumw2();
	h_hadOverEm_endcap->Sumw2();
	h_EOverPIn_endcap ->Sumw2();
	h_dEtaIn_endcap ->Sumw2();
	h_dPhiIn_endcap ->Sumw2();
	h_cIso_barrel->Sumw2();
	h_tIso_barrel->Sumw2();
	h_cIso_endcap->Sumw2();
	h_tIso_endcap->Sumw2();

	//Kinematic Plots
	TH1F *pass_ht = new TH1F("pass_ht", "HT (after all cuts)", 100, 0.0, 1000.0);
	TH1F *pass_met = new TH1F("pass_met", "MET (after all cuts)", 200, 0.0, 200.0);
	pass_ht->Sumw2();
	pass_met->Sumw2();

	TDirectory *dir_conv = histf->mkdir("conversion", "conversion");
	dir_conv->cd();
	//81FB
	Conv_Opti[0] = new TH2D("ConvOptimization_ttbar", "ConvOptimization_ttbar", 100, -1, 1, 100, -1, 1);
	Conv_OptiL[0] = new TH2D("ConvOptimization_ttbar_l", "ConvOptimization_ttbar_l", 100, -10, 10, 100, -10, 10);
	Conv_Optis[0] = new TH2D("ConvOptimization_ttbar_s", "ConvOptimization_ttbar_s", 100, -0.1, 0.1, 100, -0.1, 0.1);
	Conv_Opti[1] = new TH2D("ConvOptimization_enri3", "ConvOptimization_enri3", 100, -1, 1, 100, -1, 1);
	Conv_OptiL[1] = new TH2D("ConvOptimization_enri3_l", "ConvOptimization_enri3_l", 100, -10, 10, 100, -10, 10);
	Conv_Optis[1] = new TH2D("ConvOptimization_enri3_s", "ConvOptimization_enri3_s", 100, -0.1, 0.1, 100, -0.1, 0.1);
	Conv_Opti_extragran[0] = new TH2D("ConvOptimization_ttbar_eg", "ConvOptimization_ttbar_eg", 2000, -0.2, 0.2, 2000, -0.2, 0.2);
	Conv_Opti_extragran[1] = new TH2D("ConvOptimization_enri3_eg", "ConvOptimization_enri3_eg", 2000, -0.2, 0.2, 2000, -0.2, 0.2);

	//Conv_CheckDelR_GSFTk_ctfTk = new TH1D("Conv_CheckDelR_GSFTk_ctfTk","Conv_CheckDelR_GSFTk_ctfTk",1000,0,10);

	TDirectory *dir_Valid;
	TH1F *valid_HT[9][7];
	TH1F *valid_jetsEt[9][7];
	TH1F *valid_jetsEta[9][7];
	TH1F *valid_jetsPhi[9][7];
	TH1F *valid_jets1stEt[9][7];
	TH1F *valid_jets2ndEt[9][7];
	TH1F *valid_jets3rdEt[9][7];
	TH1F *valid_jets4thEt[9][7];
	TH1F *valid_eleEt[9][7];
	TH1F *valid_eleEta[9][7];
	TH1F *valid_elePhi[9][7];
	TH1F *valid_eleCalIso[9][7];
	TH1F *valid_eleTrkIso[9][7];
	TH1F *valid_eleRelIso[9][7];
	TH1F *valid_eled0[9][7];
	TH1F *valid_metEt[9][7];
	TH1F *valid_metPhi[9][7];
	TH1F *valid_genTT_pt[9][7];
	TH1F *valid_genT_pt[9][7];
	TH1F *valid_recoM3[9][7];
	TH1F *valid_mass_ee[9][7];
	TH1F *valid_recoM3_PTMax[9][7];
	TH1F *valid_numberTracks[9][7];
	TH1F *valid_trackPt[9][7];

	if (doValidation()) {
		dir_Valid = histf->mkdir("Valid", "Validation plots");
		dir_Valid->cd();
		valid_mkHisto_cut_njet(valid_HT, "ht", "ht", 1000, 0, 1000);
		valid_mkHisto_cut_njet(valid_jetsEt, "jetsEt", "all_jetsEt", 1000, 0, 1000);
		valid_mkHisto_cut_njet(valid_jetsEta, "jetsEta", "all_jetsEta", 1000, -5, 5);
		valid_mkHisto_cut_njet(valid_jetsPhi, "jetsPhi", "all_jetsPhi", 1000, -5, 5);
		valid_mkHisto_cut_njet(valid_jets1stEt, "jets1st", "clean_jets1st", 1000, 0, 1000);
		valid_mkHisto_cut_njet(valid_jets2ndEt, "jets2nd", "clean_jets2nd", 1000, 0, 1000);
		valid_mkHisto_cut_njet(valid_jets3rdEt, "jets3rd", "clean_jets3rd", 1000, 0, 1000);
		valid_mkHisto_cut_njet(valid_jets4thEt, "jets4th", "clean_jets4th", 1000, 0, 1000);
		valid_mkHisto_cut_njet(valid_eleEt, "eleEt", "eleEt", 1000, 0, 1000);
		valid_mkHisto_cut_njet(valid_eleEta, "eleEta", "eleEta", 1000, -5, 5);
		valid_mkHisto_cut_njet(valid_elePhi, "elePhi", "elePhi", 1000, -5, 5);
		valid_mkHisto_cut_njet(valid_eleCalIso, "eleCalIso", "eleCalIso", 1000, 0, 1);
		valid_mkHisto_cut_njet(valid_eleTrkIso, "eleTrkIso", "eleTrkIso", 1000, 0, 1);
		valid_mkHisto_cut_njet(valid_eleRelIso, "eleRelIso", "eleRelIso", 1000, 0, 1);
		valid_mkHisto_cut_njet(valid_eled0, "eled0", "eled0", 2000, 0, 2);
		valid_mkHisto_cut_njet(valid_metEt, "metEt", "metEt", 1000, 0, 1000);
		valid_mkHisto_cut_njet(valid_metPhi, "metPhi", "metPhi", 1000, -5, 5);
		valid_mkHisto_cut_njet(valid_genTT_pt, "genTT_pt", "genTT_pt", 1000, 0, 1000);
		valid_mkHisto_cut_njet(valid_genT_pt, "genT_pt", "genT_pt", 1000, 0, 1000);
		valid_mkHisto_cut_njet(valid_recoM3, "recoM3", "recoM3", 1000, 0, 1000);
		valid_mkHisto_cut_njet(valid_mass_ee, "mass_ee", "mass_ee", 1000, 0, 500);
		valid_mkHisto_cut_njet(valid_recoM3_PTMax, "recoM3_PTMax", "recoM3_PTMax", 1000, 0, 1000);
		valid_mkHisto_cut_njet(valid_numberTracks, "numberTracks", "numberTracks", 500, 0, 500);
		valid_mkHisto_cut_njet(valid_trackPt, "trackPt", "trackPt", 500, 0, 500);

	}
	//   valid_mkHisto_cut_njet(valid_ ,"","", , ,);
	//81FB end

	histf->cd();

	//TL: W+jets estimation: m3 plots (12-1-09)
	TDirectory *dir_wjets = histf->mkdir("Wjets_estimation", "m3 plots");
	dir_wjets->cd();
	// m3 for signal region (isolated ele) ---> TODO: how many bins to use, and what range for M3 fit
	h_hadTop_maxPT_mass_4j = new TH1D("hadTop_maxPT_mass_4j", "had top mass (m3) (>=4j)", nbm3, 0, nbm3); //<----
	h_hadTop_maxPT_pt_4j = new TH1D("hadTop_maxPT_pt_4j", "had top highest PT (>=4j)", 100, 0, 1000);
	// m3 for background region (non-isolated ele)
	h_hadTop_maxPT_mass_nonIso_4j = new TH1D("hadTop_maxPT_mass_nonIso_4j", "had top mass (m3) (>=4j) nonIsoE", nbm3, 0, nbm3);
	h_hadTop_maxPT_pt_nonIso_4j = new TH1D("hadTop_maxPT_pt_nonIso_4j", "had top highest PT (>=4j) nonIsoE", 100, 0, 1000);

	h_hadTop_maxPT_mass_4j->Sumw2();
	h_hadTop_maxPT_pt_4j ->Sumw2();
	h_hadTop_maxPT_mass_nonIso_4j->Sumw2();
	h_hadTop_maxPT_pt_nonIso_4j ->Sumw2();

	TDirectory *dir_wj2 = histf->mkdir("Wjets_estimation_nB1000", "m3 plots (1000 bins, 0-1000)");
	dir_wj2->cd();
	h_hadTop_maxPT_mass_4j_1000 = new TH1D("hadTop_maxPT_mass_4j", "had top mass (m3) (>=4j)", 1000, 0, 1000); //<----
	h_hadTop_maxPT_mass_nonIso_4j_1000 = new TH1D("hadTop_maxPT_mass_nonIso_4j", "had top mass (m3) (>=4j) nonIsoE", 1000, 0,
			1000);
	h_hadTop_maxPT_mass_4j_1000->Sumw2();
	h_hadTop_maxPT_mass_nonIso_4j_1000->Sumw2();

	// m3 for each type of MC (if running on MC)
	if (IsData() == false) {
		// (i) signal region
		//
		// TODO: To-be-reviewed: Number of bins used for M3
		//
		dir_wjets->cd();
		h_m3_tt = new TH1D("m3_tt", "m3 (ttbar MC)", nbm3, 0, nbm3);
		h_m3_wj = new TH1D("m3_wj", "m3 (W+jets MC)", nbm3, 0, nbm3);
		h_m3_zj = new TH1D("m3_zj", "m3 (Z+jets MC)", nbm3, 0, nbm3);
		h_m3_qcd = new TH1D("m3_qcd", "m3 (QCD MC)", nbm3, 0, nbm3);
		h_m3_vqq = new TH1D("m3_vqq", "m3 (VQQ MC)", nbm3, 0, nbm3);
		h_m3_singletop = new TH1D("m3_singletop", "m3 (single top MC)", nbm3, 0, nbm3);
		h_m3_bce[0] = new TH1D("m3_bce1", "m3 (QCD bce1 MC)", nbm3, 0, nbm3);
		h_m3_bce[1] = new TH1D("m3_bce2", "m3 (QCD bce2 MC)", nbm3, 0, nbm3);
		h_m3_bce[2] = new TH1D("m3_bce3", "m3 (QCD bce3 MC)", nbm3, 0, nbm3);
		h_m3_enri[0] = new TH1D("m3_enri1", "m3 (QCD enri1 MC)", nbm3, 0, nbm3);
		h_m3_enri[1] = new TH1D("m3_enri2", "m3 (QCD enri2 MC)", nbm3, 0, nbm3);
		h_m3_enri[2] = new TH1D("m3_enri3", "m3 (QCD enri3 MC)", nbm3, 0, nbm3);
		// (ii) control region
		h_m3_tt_control = new TH1D("m3_tt_control", "m3 (ttbar MC, control reg) ", nbm3, 0, nbm3);
		h_m3_wj_control = new TH1D("m3_wj_control", "m3 (W+jets MC, control reg)", nbm3, 0, nbm3);
		h_m3_zj_control = new TH1D("m3_zj_control", "m3 (Z+jets MC, control reg)", nbm3, 0, nbm3);
		h_m3_qcd_control = new TH1D("m3_qcd_control", "m3 (QCD MC, control reg)", nbm3, 0, nbm3);
		h_m3_vqq_control = new TH1D("m3_vqq_control", "m3 (VQQ MC, control reg)", nbm3, 0, nbm3);
		h_m3_singletop_control = new TH1D("m3_singletop_control", "m3 (single top MC, control reg)", nbm3, 0, nbm3);
		h_m3_bce_control[0] = new TH1D("m3_bce1_control", "m3 (QCD bce1 MC, control reg)", nbm3, 0, nbm3);
		h_m3_bce_control[1] = new TH1D("m3_bce2_control", "m3 (QCD bce2 MC, control reg)", nbm3, 0, nbm3);
		h_m3_bce_control[2] = new TH1D("m3_bce3_control", "m3 (QCD bce3 MC, control reg)", nbm3, 0, nbm3);
		h_m3_enri_control[0] = new TH1D("m3_enri1_control", "m3 (QCD enri1 MC, control reg)", nbm3, 0, nbm3);
		h_m3_enri_control[1] = new TH1D("m3_enri2_control", "m3 (QCD enri2 MC, control reg)", nbm3, 0, nbm3);
		h_m3_enri_control[2] = new TH1D("m3_enri3_control", "m3 (QCD enri3 MC, control reg)", nbm3, 0, nbm3);

		// 1000 bins
		dir_wj2->cd();
		// (i) signal region
		h_m3_tt_1000 = new TH1D("m3_tt", "m3 (ttbar MC)", 1000, 0, 1000);
		h_m3_wj_1000 = new TH1D("m3_wj", "m3 (W+jets MC)", 1000, 0, 1000);
		h_m3_zj_1000 = new TH1D("m3_zj", "m3 (Z+jets MC)", 1000, 0, 1000);
		h_m3_qcd_1000 = new TH1D("m3_qcd", "m3 (QCD MC)", 1000, 0, 1000);
		h_m3_vqq_1000 = new TH1D("m3_vqq", "m3 (VQQ MC)", 1000, 0, 1000);
		h_m3_singletop_1000 = new TH1D("m3_singletop", "m3 (single top MC)", 1000, 0, 1000);
		h_m3_bce_1000[0] = new TH1D("m3_bce1", "m3 (QCD bce1 MC)", 1000, 0, 1000);
		h_m3_bce_1000[1] = new TH1D("m3_bce2", "m3 (QCD bce2 MC)", 1000, 0, 1000);
		h_m3_bce_1000[2] = new TH1D("m3_bce3", "m3 (QCD bce3 MC)", 1000, 0, 1000);
		h_m3_enri_1000[0] = new TH1D("m3_enri1", "m3 (QCD enri1 MC)", 1000, 0, 1000);
		h_m3_enri_1000[1] = new TH1D("m3_enri2", "m3 (QCD enri2 MC)", 1000, 0, 1000);
		h_m3_enri_1000[2] = new TH1D("m3_enri3", "m3 (QCD enri3 MC)", 1000, 0, 1000);
		// (ii) control region
		h_m3_tt_control_1000 = new TH1D("m3_tt_control", "m3 (ttbar MC, control reg) ", 1000, 0, 1000);
		h_m3_wj_control_1000 = new TH1D("m3_wj_control", "m3 (W+jets MC, control reg)", 1000, 0, 1000);
		h_m3_zj_control_1000 = new TH1D("m3_zj_control", "m3 (Z+jets MC, control reg)", 1000, 0, 1000);
		h_m3_qcd_control_1000 = new TH1D("m3_qcd_control", "m3 (QCD MC, control reg)", 1000, 0, 1000);
		h_m3_vqq_control_1000 = new TH1D("m3_vqq_control", "m3 (VQQ MC, control reg)", 1000, 0, 1000);
		h_m3_singletop_control_1000 = new TH1D("m3_singletop_control", "m3 (single top MC, control reg)", 1000, 0, 1000);
		h_m3_bce_control_1000[0] = new TH1D("m3_bce1_control", "m3 (QCD bce1 MC, control reg)", 1000, 0, 1000);
		h_m3_bce_control_1000[1] = new TH1D("m3_bce2_control", "m3 (QCD bce2 MC, control reg)", 1000, 0, 1000);
		h_m3_bce_control_1000[2] = new TH1D("m3_bce3_control", "m3 (QCD bce3 MC, control reg)", 1000, 0, 1000);
		h_m3_enri_control_1000[0] = new TH1D("m3_enri1_control", "m3 (QCD enri1 MC, control reg)", 1000, 0, 1000);
		h_m3_enri_control_1000[1] = new TH1D("m3_enri2_control", "m3 (QCD enri2 MC, control reg)", 1000, 0, 1000);
		h_m3_enri_control_1000[2] = new TH1D("m3_enri3_control", "m3 (QCD enri3 MC, control reg)", 1000, 0, 1000);

		h_m3_tt ->Sumw2();
		h_m3_wj ->Sumw2();
		h_m3_zj ->Sumw2();
		h_m3_qcd->Sumw2();
		h_m3_vqq->Sumw2();
		h_m3_singletop->Sumw2();
		h_m3_tt_control ->Sumw2();
		h_m3_wj_control ->Sumw2();
		h_m3_zj_control ->Sumw2();
		h_m3_qcd_control->Sumw2();
		h_m3_vqq_control->Sumw2();
		h_m3_singletop_control->Sumw2();
		h_m3_tt_1000 ->Sumw2(); //1000 bins
		h_m3_wj_1000 ->Sumw2();
		h_m3_zj_1000 ->Sumw2();
		h_m3_qcd_1000->Sumw2();
		h_m3_vqq_1000->Sumw2();
		h_m3_singletop_1000->Sumw2();
		h_m3_tt_control_1000 ->Sumw2();
		h_m3_wj_control_1000 ->Sumw2();
		h_m3_zj_control_1000 ->Sumw2();
		h_m3_qcd_control_1000->Sumw2();
		h_m3_vqq_control_1000->Sumw2();
		h_m3_singletop_control_1000->Sumw2();

		for (short i = 0; i < 3; ++i) {
			h_m3_bce[i]->Sumw2();
			h_m3_enri[i]->Sumw2();
			h_m3_bce_control[i]->Sumw2();
			h_m3_enri_control[i]->Sumw2();
			h_m3_bce_1000[i]->Sumw2();
			h_m3_enri_1000[i]->Sumw2();
			h_m3_bce_control_1000[i]->Sumw2();
			h_m3_enri_control_1000[i]->Sumw2();
		}

	}//if MC


	histf->cd();

	//// const int ntjet(5);    //made as global

	// Collect event count after selection
	TH2D *Signal_njetsVcuts = new TH2D("Signal_njetsVcuts", "Events V Cuts (signal)", ntjet + 1, 0, ntjet + 1, nstage, 0, nstage);
	TH2D *QCD_njetsVcuts = new TH2D("QCD_njetsVcuts", "Events V Cuts (QCD)", ntjet + 1, 0, ntjet + 1, nstage, 0, nstage);
	TH2D *Wjets_njetsVcuts = new TH2D("Wjets_njetsVcuts", "Events V Cuts (W+jets)", ntjet + 1, 0, ntjet + 1, nstage, 0, nstage);
	TH2D *Zjets_njetsVcuts = new TH2D("Zjets_njetsVcuts", "Events V Cuts (Z+jets)", ntjet + 1, 0, ntjet + 1, nstage, 0, nstage);
	TH2D *VQQ_njetsVcuts = new TH2D("VQQ_njetsVcuts", "Events V Cuts (VQQ)", ntjet + 1, 0, ntjet + 1, nstage, 0, nstage);
	TH2D *SingleTop_njetsVcuts = new TH2D("SingleTop_njetsVcuts", "Events V Cuts (Single top)", ntjet + 1, 0, ntjet + 1, nstage,
			0, nstage);

	// Histogram to store signal eff.acc
	TH1D *SignalVars = new TH1D("SignalVars", "Signal acceptance and efficiency", 4, 0, 4);

	//--------------------------------------------------
	// [MC] PDF uncertainties - signal accep.eff (9 Feb 2010)
	//--------------------------------------------------
	TDirectory *dir_pdfunc(0);
	TH1F *h_pdf_total(0);
	TH1F *h_pdf_pass(0);
	TH1F *h_pdf_eff(0);
	if (!IsData() && m_studyPDFunc) { //for signal only
		dir_pdfunc = histf->mkdir("pdf_unc", "PDF uncertainty (signal accep.eff)");
		dir_pdfunc->cd();
		h_pdf_total = new TH1F("pdf_total_sig", "N(Signal event) Total", 46, -1, 45);
		h_pdf_pass = new TH1F("pdf_pass_sig", "N(Signal event) Pass", 46, -1, 45);
		h_pdf_eff = new TH1F("pdf_eff_sig", "Signal Accep.eff", 46, -1, 45);
		h_pdf_total->GetXaxis()->SetTitle("PDF set");
		h_pdf_pass->GetXaxis()->SetTitle("PDF set");
		h_pdf_eff->GetXaxis()->SetTitle("PDF set");
		h_pdf_total->SetFillColor(kGray);
		h_pdf_pass->SetFillColor(kGray);
		h_pdf_eff->SetMarkerStyle(20);
		h_pdf_total->Sumw2();
		h_pdf_pass->Sumw2();
		h_pdf_eff->Sumw2();
	}

	//-----------
	// Cut names
	//-----------
	vector<string> ve;
	ve.push_back("INITIAL   ");
	ve.push_back("GOODRUN   ");
	ve.push_back("TRIGGER   ");
	ve.push_back("$\\ge$1TEle"); //NEW
	ve.push_back("$\\ge$1TISO");
	ve.push_back("=1TISO    ");
	ve.push_back("!MUON");
	ve.push_back("MET       ");
	ve.push_back("!Z        ");
	ve.push_back("!CONV     ");
	if (!m_rejectEndcapEle)
		ve.push_back("!DIFFZ    ");
	else
		ve.push_back("BARREL    ");
	ve.push_back("HT        ");
	ve.push_back("TAGGABLE  ");
	ve.push_back("$\\ge$1+BTAG");
	ve.push_back("$\\ge$2+BTAG");
	ve.push_back("$\\ge$1-BTAG");

	//For MC, added index for ttbar event type
	double e_plus_jet[nstage][ntjet][nmctype]; //, mu_plus_jet[nstage][ntjet][nmctype];
	double e_plus_jet_weighted[nstage][ntjet][nmctype]; //TL 15-1-09: create a new array for weighted-sum?

	//Ignore mu_plus_jet, for a mu_plus_jet analysis only
	for (short i = 0; i < nstage; ++i) {
		for (short j = 0; j < ntjet; ++j) {
			for (short k = 0; k < nmctype; ++k) {
				e_plus_jet[i][j][k] = 0;
				e_plus_jet_weighted[i][j][k] = 0;
			}
		}
	}

	//Vectors to store 4-vectors for selected objects
	std::vector<TLorentzVector> electrons;
	std::vector<unsigned int> ii_electrons; //Added 28-9-09 (index of GoodEle)
	std::vector<TLorentzVector> electrons_barrel;//SEARCH: vector definitions
	std::vector<TLorentzVector> electrons_endcap;
	std::vector<TLorentzVector> iso_electrons;
	std::vector<TLorentzVector> iso_electrons_barrel;
	std::vector<TLorentzVector> iso_electrons_endcap;
	std::vector<TLorentzVector> muons;
	std::vector<TLorentzVector> iso_muons;
	std::vector<TLorentzVector> goodJets;
	std::vector<std::vector<double> > btag_information;
	btag_information.resize(NUMBER_OF_BTAGS);
	std::vector<TLorentzVector> met;
	std::vector<float> electrons_isoval;
	std::vector<float> electrons_isoval2;

	unsigned int index_selected_ele = 0; //index of the selected isolated electron


	int nZ_TEST = 0; //TEST
	int nZ_ORIG = 0; //TEST


	int counter = 0;
	int counter_pass = 0;
	int printLevel = 0;

	cout << endl << " Starting! nEvents is " << nEvents << " counter is " << counter << endl;
	if (GetLimit() > 0)
		nEvents = GetLimit();
	cout << " Limited nEvents to " << nEvents << endl;

	// MC information (declared before event loop)
	//	string this_mc = "data"; //default
	int mctype = 0;

	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//
	//                         Loop over the events
	//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	fprintf(outfile, "%7s %6s %10s %10s %8s %8s %6s %13s %6s %8s %8s   %s\n", "no", "run", "event", "lumiSec", "tree", "entry",
			"ntj", "leadingJetET", "met", "HT", "mc", "file");
	fprintf(outfile, "---------------------------------------------------------");
	fprintf(outfile, "---------------------------------------------------------\n");
	if (m_debug)
		cout << "\n Will start main Event Loop" << endl;

	// Z stuff
	int nEvent_Z = 0;
	int nEvent_Zee = 0;
	int nEvent_Zmm = 0;
	int nEvent_Ztt = 0;
	int nEvent_2orMoreEle = 0;
	int nEvent_EleMatch = 0;
	int nEvent_Z_et = 0;
	int nEvent_Z_eta = 0;
	int nEvent_Z_eteta = 0;
	int nEvent_Z_eteta_d0 = 0;
	int nEvent_Z_eteta_d0_L = 0;
	int nEvent_Z_eteta_d0_T = 0;
	int nEvent_Z_eteta_d0_RL = 0;
	int nEvent_Z_eteta_d0_RT = 0;
	int nEvent_Zee_pass = 0;
	int nEvent_Zee_highMET = 0;
	int nEvent_Zee_lowMET = 0;
	int nEvent_Zmm_pass = 0;
	int nEvent_Zmm_highMET = 0;
	int nEvent_Zmm_lowMET = 0;
	int nEvent_Ztt_pass = 0;
	int nEvent_Ztt_highMET = 0;
	int nEvent_Ztt_lowMET = 0;

	// Monitor number of events with DR(e,mu) < 0.1
	int nEvent_DR_ele_muo_less_than_01 = 0;

	// record list of "interesting" files that contain selected events
	set<string> interestingFiles;

	//----------------------------------------------------------------------------
	//                   Beginning  of  Event  Loop
	//----------------------------------------------------------------------------

	for (Long64_t ev = 0; ev < nEvents; ++ev) {

		// NB: LoadTree will complain (harmlessly) about unknown branch if SetBranchAddress
		//     was called for a branch not present in the ntuple
		Long64_t lflag = chain->LoadTree(ev);
		if (lflag < 0)
			break;

		++counter;
		//if(counter<900) continue;


		int nbytes = chain->GetEntry(ev);
		//     if(ev>100) m_debug = false; // turn off m_debug message after 100 events


		// 1- Check for good run (only applicable to data)
		bool goodrun = true;

		if (IsData() && KeepGoodRuns()) {
			//goodrun = mytool->good_run(run);
		}

		if (m_debug || ev < 10 || ev % 5000 == 0) {
			if (m_debug || ev < 10)
				cout << "\nBegin processing event no. " << ev + 1;
			else
				cout << "\nBegin processing event no. " << (ev + 1) / 1000 << " k";
			cout << ". GoodRun=" << goodrun;
			cout << "  << Run " << run << ", Event " << event << ", LumiSection " << lumiBlock << " >>  ";
			cout << printTimeNow() << endl;
			cout << " entry " << lflag << ", tree # " << chain->GetTreeNumber() << ", file "
					<< chain->GetCurrentFile()->GetName() << ", EvtSize " << nbytes << endl;
		}

		// 2- Check Trigger
		if (m_debug)
			cout << " Checking Trigger" << endl;

		bool fired_single_em = true;

		if (GetTrigger()) {

			fired_single_em = false; //reset default value to "fail" if checking trigger

			if (ev == 0)
				cout << "Check HLT trigger: " << HLTBit << endl; //nonIso electron trigger
			// fired_single_em = (bool)HLT_Ele15_LW_L1R;
			fired_single_em = passHLT();

			if (ev < 10) { //print first 10 events
				if (fired_single_em == true)
					cout << " PASS Single Electron Trigger" << endl;
				else
					cout << " FAIL Single Electron Trigger" << endl;
			}
		}

		// If running on MC, check type of MC from input filename
		// note: it is sufficient to check once per tree
		// However, if run on soup, then need to check every event
		bool RunOnMixedMC = false; //29-9-09

		if (IsData() == false) {

			if (RunOnMixedMC || //check once per event
					!RunOnMixedMC && lflag == 0) { //check once per input file

				// Reset flags
				fastmctype_ = kdata;

				if (!RunOnMixedMC) {
					if (m_debug)
						cout << " [MC] Checking MC type from input filename" << endl;

					//cout << "\nThis is the 1st event in the tree" << endl;
					//cout << "Global event no. "<< ev <<". Entry number in this tree: " << lflag << endl; //TL: 14-1-09
					const string fname(chain->GetCurrentFile()->GetName());
					//cout << "file: " << fname << endl;

					// figure out what type of MC from the file name

					int nmatch = 0;

					for (size_t i = kdata + 1; i < mc_names.size(); ++i) {
						if (fname.find(mc_names[i]) != string::npos) {
							//							this_mc = mc_names[i];
							fastmctype_ = (MCType) i;
							cout << mc_names[i] << " matched to " << fastmctype_ << " : " << mclabel[fastmctype_] << endl;
							nmatch++;
							break; //exit for loop when a match is found
						}
					}
					if (nmatch == 0)
						cout << "ERROR. File name does not match to any type of MC" << endl;

				} else {

					//29-9-09
					if (m_debug)
						cout << " [MC] Checking Event type from MC truth" << endl;
					//this_mc = CheckEventTypeFromMcTruth(); //TEMPORARY
					//cout << "this_mc = " << this_mc << endl;
					cout << "ERROR: Should not use this piece of code!" << endl;
				}

				// set event weight of current event
				//  this_weight = weightMap[this_mc];
				//				this_weight = GetWeight(this_mc);
				this_weight = GetWeight(fastmctype_);
				if (!isMCType(kttbar))
					mctype = fastmctype_ + NUMBER_OF_DECAY_BRANCHES;//equivalent to the lines below
				//				if (m_debug)
				cout << "FAST MC Type is: " << fastmctype_ << endl;
				cout << "mctype is: " << mctype << endl;

				// set nmctype code
				//				if (this_mc == "ttbar" || this_mc == "ttjet") {
				//					isMCType(kttbar) = true;
				//					fastmctype_ = kttbar;
				//				} //mctype 1-10
				//				else if (this_mc == "wjet" || this_mc == "wenu") {
				//					isMCType(kWjets) = true;
				//					mctype = 11;
				//					fastmctype_ = kWjets;
				//				} else if (this_mc == "zjet" || this_mc == "zee") {
				//					isMCType(kZjets) = true;
				//					mctype = 12;
				//					fastmctype_ = kZjets;
				//				} else if (this_mc == "enri1") {
				//					isMCType(kQCD) = true;
				//					isMCType(kenri1) = true;
				//					mctype = 13;
				//					fastmctype_ = kenri1;
				//				} else if (this_mc == "enri2") {
				//					isMCType(kQCD) = true;
				//					isMCType(kenri2) = true;
				//					mctype = 14;
				//					fastmctype_ = kenri2;
				//				} else if (this_mc == "enri3") {
				//					isMCType(kQCD) = true;
				//					isMCType(kenri3) = true;
				//					mctype = 15;
				//					fastmctype_ = kenri3;
				//				} else if (this_mc == "bce1") {
				//					isMCType(kQCD) = true;
				//					isMCType(kbce1) = true;
				//					mctype = 16;
				//					fastmctype_ = kbce1;
				//				} else if (this_mc == "bce2") {
				//					isMCType(kQCD) = true;
				//					isMCType(kbce2) = true;
				//					mctype = 17;
				//					fastmctype_ = kbce2;
				//				} else if (this_mc == "bce3") {
				//					isMCType(kQCD) = true;
				//					isMCType(kbce3) = true;
				//					mctype = 18;
				//					fastmctype_ = kbce3;
				//				} else if (this_mc == "vqq") {
				//					isMCType(kVQQ) = true;
				//					mctype = 19;
				//					fastmctype_ = kVQQ;
				//				} else if (this_mc == "tW") {
				//					isMCType(ksingleTop) = true;
				//					isMCType(ktW) = true;
				//					mctype = 20;
				//					fastmctype_ = ktW;
				//				} else if (this_mc == "tchan") {
				//					isMCType(ksingleTop) = true;
				//					isMCType(ktchan) = true;
				//					mctype = 21;
				//					fastmctype_ = ktchan;
				//				} else if (this_mc == "schan") {
				//					isMCType(ksingleTop) = true;
				//					isMCType(kschan) = true;
				//					mctype = 22;
				//					fastmctype_ = kschan;
				//				}
				//				//FIXME: 10th step to new sample - set mctype and is*
				//				else if (this_mc == "Zprime_M500GeV_W5GeV") {
				//					isMCType(Zprime_M500GeV_W5GeV) = true;
				//					mctype = 23;
				//					fastmctype_ = Zprime_M500GeV_W5GeV;
				//				} else if (this_mc == "Zprime_M500GeV_W50GeV") {
				//					isMCType(Zprime_M500GeV_W50GeV) = true;
				//					mctype = 24;
				//					fastmctype_ = Zprime_M500GeV_W50GeV;
				//				}
				//
				//				else if (this_mc == "Zprime_M750GeV_W7500MeV") {
				//					isMCType(Zprime_M750GeV_W7500MeV) = true;
				//					mctype = 25;
				//					fastmctype_ = Zprime_M750GeV_W7500MeV;
				//				}
				//
				//				else if (this_mc == "Zprime_M1TeV_W10GeV") {
				//					isMCType(Zprime_M1TeV_W10GeV) = true;
				//					mctype = 26;
				//					fastmctype_ = Zprime_M1TeV_W10GeV;
				//				} else if (this_mc == "Zprime_M1TeV_W100GeV") {
				//					isMCType(Zprime_M1TeV_W100GeV) = true;
				//					mctype = 27;
				//					fastmctype_ = Zprime_M1TeV_W100GeV;
				//				}
				//
				//				else if (this_mc == "Zprime_M1250GeV_W12500MeV") {
				//					isMCType(Zprime_M1250GeV_W12500MeV) = true;
				//					mctype = 28;
				//					fastmctype_ = Zprime_M1250GeV_W12500MeV;
				//				}
				//
				//				else if (this_mc == "Zprime_M1500GeV_W15GeV") {
				//					isMCType(Zprime_M1500GeV_W15GeV) = true;
				//					mctype = 29;
				//					fastmctype_ = Zprime_M1500GeV_W15GeV;
				//				} else if (this_mc == "Zprime_M1500GeV_W150GeV") {
				//					isMCType(Zprime_M1500GeV_W150GeV) = true;
				//					mctype = 30;
				//					fastmctype_ = Zprime_M1500GeV_W150GeV;
				//				} else if (this_mc == "Zprime_M2TeV_W20GeV") {
				//					isMCType(Zprime_M2TeV_W20GeV) = true;
				//					mctype = 31;
				//					fastmctype_ = Zprime_M2TeV_W20GeV;
				//				} else if (this_mc == "Zprime_M2TeV_W200GeV") {
				//					isMCType(Zprime_M2TeV_W200GeV) = true;
				//					mctype = 32;
				//					fastmctype_ = Zprime_M2TeV_W200GeV;
				//				} else if (this_mc == "Zprime_M3TeV_W30GeV") {
				//					isMCType(Zprime_M3TeV_W30GeV) = true;
				//					mctype = 33;
				//					fastmctype_ = Zprime_M3TeV_W30GeV;
				//				} else if (this_mc == "Zprime_M3TeV_W300GeV") {
				//					isMCType(Zprime_M3TeV_W300GeV) = true;
				//					mctype = 34;
				//					fastmctype_ = Zprime_M3TeV_W300GeV;
				//				} else if (this_mc == "Zprime_M4TeV_W40GeV") {
				//					isMCType(Zprime_M4TeV_W40GeV) = true;
				//					mctype = 35;
				//					fastmctype_ = Zprime_M4TeV_W40GeV;
				//				} else if (this_mc == "Zprime_M4TeV_W400GeV") {
				//					isMCType(Zprime_M4TeV_W400GeV) = true;
				//					mctype = 36;
				//					fastmctype_ = Zprime_M4TeV_W400GeV;
				//				}

				if (m_debug) {
					//					cout << "+++> Now starts running on  " << this_mc << "  events" << endl;
					cout << "+++> Now starts running on  " << mc_names[fastmctype_] << "  events" << endl;
					cout << " current entry: num " << lflag << ", tree # " << chain->GetTreeNumber() << "\n filename: "
							<< chain->GetCurrentFile()->GetName() << endl << endl;
				}

			}//end if check MC
		}//end if run on MC


		//--------------------
		//  For TTbar MC
		//--------------------
		if (!IsData() && isMCType(kttbar)) { //If it isn't data figure out what kind of MC (for ttbar)

			if (m_debug)
				cout << " [MC] signal, checking type of ttbar decay." << endl;

			if (Nmc_doc == 0) {
				mctype = 1;
				static bool first_time = true;
				if (first_time) {
					cout
							<< " [MC] Nmc_doc = 0. No MC info in ntuple! Set mctype to 1 so that all ttbar events appear in ev column"
							<< endl;
					first_time = false;
				}
			} else {

				//Check ttbar decay
				int nlep = 0;
				int nt = 0;
				int nW = 0;
				int nb = 0;
				int nWfromt = 0;
				int nbfromt = 0;

				if (printLevel > 0)
					PrintGenParticles();

				//Loop over mc documentation lines
				for (unsigned int i = 0; i < Nmc_doc; ++i) {

					if (fabs(mc_doc_id->at(i)) == 6)
						nt++;
					if (fabs(mc_doc_id->at(i)) == 24)
						nW++;
					if (fabs(mc_doc_id->at(i)) == 24 && fabs(mc_doc_mother_id->at(i)) == 6)
						nWfromt++;
					if (fabs(mc_doc_id->at(i)) == 5)
						nb++;
					if (fabs(mc_doc_id->at(i)) == 5 && fabs(mc_doc_mother_id->at(i)) == 6)
						nbfromt++;

					if (fabs(mc_doc_mother_id->at(i)) == 24) { //Look at W daughters
						if (fabs(mc_doc_id->at(i)) == 11)
							nlep += 1; //add 1 for an electron
						if (fabs(mc_doc_id->at(i)) == 13)
							nlep += 3; //add 3 for a muon
						if (fabs(mc_doc_id->at(i)) == 15)
							nlep += 7; //add 7 for a tau
					}
				}

				//some reassignment for easy print-out in table form
				if (nlep == 1)
					mctype = 1; //evqq
				else if (nlep == 3)
					mctype = 2; //mvqq
				else if (nlep == 7)
					mctype = 3; //tvqq
				else if (nlep == 2)
					mctype = 4; //evev
				else if (nlep == 6)
					mctype = 5; //mvmv
				else if (nlep == 14)
					mctype = 6; //tvtv
				else if (nlep == 4)
					mctype = 7; //evmv
				else if (nlep == 8)
					mctype = 8; //evtv
				else if (nlep == 10)
					mctype = 9; //mvtv
				else if (nlep == 0)
					mctype = 10; //qqqq
				else
					std::cout << " ERROR:  nlep not defined " << nlep << std::endl;

				if (mctype > 10) {
					mctype = 10;
					std::cout << " ERROR: mctype is greater than 10! " << std::endl;
				}

				all_mctype->Fill(mctype); //TL note: unweighted. fill with weight?
				if (printLevel > 0)
					std::cout << " ttbar code " << nlep << " " << mctype << std::endl;
			}//if Nmc_doc>0
		}//if ttbar MC
		//		cout << "MC Type is: " << mctype << endl;
		//-------------------------------
		// Inspecting Z+jets M C events
		//-------------------------------
		if (m_debug)
			cout << " Starting Inspection of Z+jets MC" << endl;
		bool isZee = false;
		bool isZmm = false;
		bool isZtt = false;
		vector<TLorentzVector> Zele;

		if (!IsData() && isMCType(kZjets) && Nmc_doc > 0) {

			//if(m_debug) cout << "checking Z+jets MC"<< endl;
			//Loop over mc documentation lines, look for Z->ee
			for (unsigned int i = 0; i < Nmc_doc; ++i) {

				if (mc_doc_status->at(i) != 3)
					continue; //documentation line only

				//if( fabs(mc_doc_id->at(i))==11 && fabs(mc_doc_mother_id->at(i))==23 ) {

				// Look for electron
				if (fabs(mc_doc_id->at(i)) == 11) { //don't check Z parent as in some events, the Z is missing

					isZee = true;
					Zele.push_back(TLorentzVector(mc_doc_px->at(i), mc_doc_py->at(i), mc_doc_pz->at(i), mc_doc_energy->at(i)));
					//plot eta/ET of Z electrons
					//special fill Zee only
					fasthist_[fastmctype_][h_Zee_eta]->Fill(mc_doc_eta->at(i));
					fasthist_[fastmctype_][h_Zee_pt]->Fill(mc_doc_pt->at(i));

					if (ev < 10) {
						cout << "\nZj MC:  event #" << ev << endl;
						cout << "generator ele:  id = " << mc_doc_id->at(i) << "  pt = " << mc_doc_pt->at(i) << "  eta = "
								<< mc_doc_eta->at(i) << endl;
					}

				}//is Z->ee
				//if( fabs(mc_doc_id->at(i))==13 && fabs(mc_doc_mother_id->at(i))==23 ) isZmm=true;
				//if( fabs(mc_doc_id->at(i))==15 && fabs(mc_doc_mother_id->at(i))==23 ) isZtt=true;
				if (fabs(mc_doc_id->at(i)) == 13)
					isZmm = true;
				if (fabs(mc_doc_id->at(i)) == 15)
					isZtt = true;
			}//loop over genPar

			// loop over reco-photons in Z+jets events
			for (unsigned int i = 0; i < Nphotons; ++i) {
				fasthist_[fastmctype_][h_Z_photon_eta]->Fill(photons_eta->at(i));
				fasthist_[fastmctype_][h_Z_photon_et]->Fill(photons_et->at(i));
				if (isZee) {
					fasthist_[fastmctype_][h_Zee_photon_eta]->Fill(photons_eta->at(i));
					fasthist_[fastmctype_][h_Zee_photon_et]->Fill(photons_et->at(i));
					fasthist_[fastmctype_][h_Zee_photon_eteta_2D]->Fill(photons_et->at(i), photons_eta->at(i));
				}
			}
			fasthist_[fastmctype_][h_Z_Nphotons]->Fill(Nphotons);
			if (isZee)
				fasthist_[fastmctype_][h_Zee_Nphotons]->Fill(Nphotons); //no cut

			nEvent_Z++;
			if (isZee)
				nEvent_Zee++;
			else if (isZmm)
				nEvent_Zmm++;
			else if (isZtt)
				nEvent_Ztt++;
			else {
				cout << "\nThis Z decay is not e,mu,tau! What is it?" << endl;
				cout << "Print all MC info" << endl;
				PrintGenParticles();
			}
			if (isZee) {
				int nGenBasicEle_Zee = 0;
				for (unsigned int i = 0; i < Zele.size(); ++i) {
					if (Zele[i].Et() > 30.0 && fabs(Zele[i].Eta()) < 2.5)
						nGenBasicEle_Zee++;
				}
				if (fastmctype_ == kZjets)
					fasthist_[fastmctype_][h_nGenBasicEle_Zee_allj]->Fill(nGenBasicEle_Zee);
			}
		}// Zjets MC


		// 2b - Get Beam Spot
		/*
		 if(m_debug) cout << " Get Beam Spot" << endl;
		 //     cout << "number of beam spot = " << NbeamSpot << endl;
		 const float bsx = beamSpot_x->at(0);
		 const float bsy = beamSpot_y->at(0);
		 //     cout << "beam spot (x,y): "<< bsx << ", " << bsy << endl;
		 */

		// 3- Find number of of electrons, muons, jets - after energy correctionss

		// 3A - Find number of electrons
		if (m_debug)
			cout << " Starting ELECTRON section" << endl;

		int nGoodEle = 0;
		int nGoodEle_barrel = 0;
		int nGoodEle_endcap = 0;
		int nGoodIsoEle = 0;
		int nGoodIsoEle_barrel = 0;
		int nGoodIsoEle_endcap = 0;
		int nGoodNonIsoEle = 0;
		int nGoodNonIsoEle_barrel = 0;
		int nGoodNonIsoEle_endcap = 0;
		float this_isoele_d0 = 0;

		int nHighETEle = 0; //pass et
		int nBasicEle = 0; //pass et,eta
		int nBasicD0Ele = 0; //pass et,eta,d0
		int nLooseEle = 0; //pass et,eta,d0,id
		int nTightEle = 0;
		int nRobustLooseEle = 0;

		//cout << "\n==> There are  " << Nels << "  reco-electrons in this event."<<endl;

		// basic
		fillHistoDataAndMC(h_nele, Nels, this_weight);

		for (unsigned int i = 0; i < Nels; ++i) {

			//cout << "\n--> electron no." << i+1 << endl;
			//cout << "  this electron ET/eta = " << els_et->at(i) << " / " << els_eta->at(i) << endl;

			//Compute Relative isolations
			//float RelCalIso = els_cIso->at(i)/els_et->at(i);
			//float RelTrkIso = els_tIso->at(i)/els_et->at(i);

			//Compute Combined RelIso for electrons
			float CombRelIso = getRelIso(i); // 11-11-09
			float CombRelIso2 = els_et->at(i) / (els_et->at(i) + els_dr04EcalRecHitSumEt->at(i) + els_dr04HcalTowerSumEt->at(i)
					+ els_tIso->at(i)); //norm


			if (m_debug)
				cout << "-> fill electron id histogram" << endl;
			fillHistoDataAndMC(h_eid, 0.0, this_weight); //count all electrons
			fillHistoDataAndMC(knele_cuts, 0.0, this_weight);

			if (els_et->at(i) > ELE_ETCUT)
				nHighETEle++;

			// Basic kinematics (14-8-09)
			// - all electrons
			fillHistoDataAndMC(h_ele_ET_all, els_et->at(i), this_weight);
			fillHistoDataAndMC(kele_eta_all, els_eta->at(i), this_weight);
			fillHistoDataAndMC(kele_phi_all, els_phi->at(i), this_weight);
			fillHistoDataAndMC(kele_iso_all, CombRelIso, this_weight);
			if (i < 3) { //first 3 electrons
				fillHistoDataAndMC(h_ele_ET_all + i + 1, els_et->at(i), this_weight);
				fillHistoDataAndMC(kele_eta_all + i + 1, els_eta->at(i), this_weight);
				fillHistoDataAndMC(kele_phi_all + i + 1, els_phi->at(i), this_weight);
				fillHistoDataAndMC(kele_iso_all + i + 1, CombRelIso, this_weight);
			}
			//FIXME: make this with function bool PassCut(EnumCuts cut)
			//Make Et and Eta cuts
			if (els_et->at(i) > ELE_ETCUT && fabs(els_eta->at(i)) < 2.5) {

				if (m_debug)
					cout << "-> this electron has passed ET and ETA cut" << endl;

				nBasicEle++;

				// Count how many electron passing each type of electron ID
				if (m_debug)
					cout << "-> filling electron id histogram (pass Et,eta cuts)" << endl;
				//				fillHistoDataAndMC(h_eid, 1., this_weight); // count all ele passing et,eta cut
				fillHistoDataAndMC(knele_cuts, 1., this_weight);

				// Calculate d0 w.r.t beam spot position
				float d0_corrected = compute_d0("electron", i);

				// Store electron d0 information
				if (m_debug)
					cout << "-> filling electron d0 histograms" << endl;
				//				fillHistoDataAndMC(h_ed0_unCor, fabs(els_d0dum->at(i)), this_weight);
				//				fillHistoDataAndMC(h_ed0,d0_corrected, this_weight);
				fillHistoDataAndMC(h_ed0_unCor, fabs(els_d0dum->at(i)), this_weight);
				fillHistoDataAndMC(h_ed0, fabs(d0_corrected), this_weight);

				// (6 Mar 09) Apply d0 cut on electron, 200 micron = 0.02cm
				if (fabs(d0_corrected) < 0.02) {

					nBasicD0Ele++;
					// Count how many electron passing each type of electron ID
					if (m_debug)
						cout << "-> filling electron id histograms (pass d0 cut)" << endl;
					fillHistoDataAndMC(knele_cuts, 2., this_weight);
					if (els_looseId->at(i) > 0)
						fillHistoDataAndMC(knele_cuts, 3., this_weight);
					if (els_tightId->at(i) > 0)
						fillHistoDataAndMC(knele_cuts, 4., this_weight);
					if (els_robustLooseId->at(i) > 0)
						fillHistoDataAndMC(knele_cuts, 5., this_weight);
					if (els_robustTightId->at(i) > 0)
						fillHistoDataAndMC(knele_cuts, 6., this_weight);
					if (els_robustHighEnergyId->at(i) > 0)
						fillHistoDataAndMC(knele_cuts, 7., this_weight);

					// Check isolation of electron
					bool isIsolated = false;
					if (useNewReliso)
						isIsolated = CombRelIso < 0.1; //new reliso
					else
						isIsolated = CombRelIso2 > 0.85; //old reliso
					if (m_debug) {
						cout << " ele  Old RelIso = " << CombRelIso2;
						cout << "  New RelIso = " << CombRelIso;
						cout << "  isIsolated = " << isIsolated << endl;
					}

					// Count number of good electron passing each eid (ignore those in eta gap: 1.442-1.56)
					if (fabs(els_eta->at(i)) < 1.442 || fabs(els_eta->at(i)) > 1.560) {
						if (els_looseId->at(i) > 0)
							nLooseEle++;
						if (els_tightId->at(i) > 0)
							nTightEle++;
						if (els_robustLooseId->at(i) > 0)
							nRobustLooseEle++;
					}

					//Separate into endcap and barrel

					if (fabs(els_eta->at(i)) < 1.442) { //Barrel

						//Find out what isolation looks like (out of the box)
						h_cIso_barrel->Fill(els_dr04EcalRecHitSumEt->at(i) + els_dr04HcalTowerSumEt->at(i), this_weight); //<-- add weight
						h_tIso_barrel->Fill(els_tIso->at(i), this_weight);

						//Plot electron ID quantities
						h_hadOverEm_barrel->Fill(els_hadOverEm->at(i), this_weight);
						h_EOverPIn_barrel->Fill(els_eOverPIn->at(i), this_weight);
						h_dEtaIn_barrel->Fill(els_dEtaIn->at(i), this_weight);
						h_dPhiIn_barrel->Fill(els_dPhiIn->at(i), this_weight);

						//Apply "Robust Tight" Electron ID

						//	     if (els_robustTightId->at(i) == true){
						if (passEleID(i)) {

							//Store 4 vector for "good" electron and increment counters
							TLorentzVector ele(els_px->at(i), els_py->at(i), els_pz->at(i), els_energy->at(i));
							nGoodEle_barrel++;
							nGoodEle++;
							electrons.push_back(ele);
							ii_electrons.push_back(i); //index
							electrons_barrel.push_back(ele);
							electrons_isoval.push_back(CombRelIso); //needed by QCD estimation
							electrons_isoval2.push_back(CombRelIso2); //needed by QCD estimation

							if (isIsolated) { //Isolated

								//Store 4 vector for isolated electron and increment counters
								nGoodIsoEle_barrel++;
								nGoodIsoEle++;
								iso_electrons.push_back(ele);
								iso_electrons_barrel.push_back(ele);
								this_isoele_d0 = fabs(d0_corrected);
								index_selected_ele = i;
							} else { //Non-Isolated

								nGoodNonIsoEle_barrel++;
								nGoodNonIsoEle++;
							}

						}
					}// end barrel
					else if (fabs(els_eta->at(i)) > 1.560 && fabs(els_eta->at(i)) < 2.5) { //Endcap

						//Find out what isolation looks like (out of the box)
						h_cIso_endcap->Fill(els_dr04EcalRecHitSumEt->at(i) + els_dr04HcalTowerSumEt->at(i), this_weight);
						h_tIso_endcap->Fill(els_tIso->at(i), this_weight);

						//Plot electron ID quantities
						h_hadOverEm_endcap->Fill(els_hadOverEm->at(i), this_weight);
						h_EOverPIn_endcap->Fill(els_eOverPIn->at(i), this_weight);
						h_dEtaIn_endcap->Fill(els_dEtaIn->at(i), this_weight);
						h_dPhiIn_endcap->Fill(els_dPhiIn->at(i), this_weight);

						//Apply "Robust Tight" Electron ID
						// if (els_robustTightId->at(i) > 0 ) {
						if (passEleID(i)) {

							TLorentzVector ele(els_px->at(i), els_py->at(i), els_pz->at(i), els_energy->at(i));
							nGoodEle_endcap++;
							nGoodEle++;
							electrons.push_back(ele);
							ii_electrons.push_back(i); //index
							electrons_endcap.push_back(ele);
							electrons_isoval.push_back(CombRelIso); //needed by QCD estimation
							electrons_isoval2.push_back(CombRelIso2); //needed by QCD estimation

							//if (CombRelIso2 > 0.85) { //Isolated (endcap)
							if (isIsolated) { //Isolated (endcap)

								//Store 4 vector for isolated electron and increment counters
								nGoodIsoEle_endcap++;
								nGoodIsoEle++;
								iso_electrons.push_back(ele);
								iso_electrons_endcap.push_back(ele);
								this_isoele_d0 = fabs(d0_corrected);
								index_selected_ele = i;
							} else { //Non-Isolated

								nGoodNonIsoEle_endcap++;
								nGoodNonIsoEle++;
							}

						}//end if pass eidRobustTight

					}// end endcap
					else { //Electrons in Gap
					}
				}// end d0 cut
			}// end et,eta cuts
		}// end electrons loop

		if (m_debug) {
			cout << "\n-- summary -- " << endl;
			cout << "Nels   = " << Nels << endl;
			cout << "nHighETEle = " << nHighETEle << endl;
			cout << "nBasicEle  = " << nBasicEle << endl;
			cout << "nBasicD0Ele = " << nBasicD0Ele << endl;
			cout << "nLooseEle  = " << nLooseEle << endl;
			cout << "nTightEle  = " << nTightEle << endl;
			cout << "nRobustLooseEle = " << nRobustLooseEle << endl;
			cout << "nGoodEle = " << nGoodEle << endl;
		}

		// (Z->ee only) for event with exactly 2 reco-ele, match them to gen Ele from Z decay
		if (isMCType(kZjets) && isZee) {
			int nEleMatch = 0;
			vector<int> recoElectronMatchToZee_index;
			if (Nels >= 2) {

				nEvent_2orMoreEle++;
				for (unsigned int i = 0; i < Nels; ++i) { //loop over reco-ele
					//cout << " reco-ele no." << i <<  endl;

					TLorentzVector this_reco_ele(els_px->at(i), els_py->at(i), els_pz->at(i), els_energy->at(i));

					for (unsigned int j = 0; j < Zele.size(); ++j) { //loop over gen-ele from Z
						//cout << " gen-ele no." << j <<  endl;

						float delR_recoE_genE = calcDeltaR(this_reco_ele, Zele.at(j));
						//cout << "DR(reco " << i << ",gen " << j << ") = "<< delR_recoE_genE << endl;
						if (delR_recoE_genE < 0.1) {
							nEleMatch++;
							recoElectronMatchToZee_index.push_back(i);
						}
					}
				}// loop over ele
			}// 2 or more ele

			int nEle_Z_et = 0;
			int nEle_Z_eta = 0;
			int nEle_Z_eteta = 0;
			int nEle_Z_eteta_d0 = 0;
			int nEle_Z_eteta_d0_RL = 0;
			int nEle_Z_eteta_d0_L = 0;
			int nEle_Z_eteta_d0_T = 0;
			int nEle_Z_eteta_d0_RT = 0;

			if (nEleMatch == 2) {
				nEvent_EleMatch++;
				//cout << "info: this event has 2 reo-electrons matched to Z->ee"<< endl;
				for (short i = 0; i < 2; ++i) {
					int k = recoElectronMatchToZee_index.at(i);

					if (els_et->at(k) > 30.0)
						nEle_Z_et++; //et cut only
					if (fabs(els_eta->at(k)) < 2.5)
						nEle_Z_eta++; //eta cut only

					if (els_et->at(k) > 30.0 && fabs(els_eta->at(k)) < 2.5) { //et and eta cut
						nEle_Z_eteta++;

						// Calculate d0 w.r.t beam spot position
						float d0_corrected = compute_d0("electron", k);

						if (fabs(d0_corrected) < 0.02) {
							nEle_Z_eteta_d0++;

							if (els_robustLooseId->at(k))
								nEle_Z_eteta_d0_RL++;
							if (els_robustTightId->at(k))
								nEle_Z_eteta_d0_RT++;
							if (els_tightId->at(k))
								nEle_Z_eteta_d0_T++;
							if (els_looseId->at(k))
								nEle_Z_eteta_d0_L++;
						}//d0
					}//et,eta
				}//ele loop
			}//if nEleMatch==2

			// count how many events with 2 or more electrons (pass each cut)
			if (nEle_Z_et >= 2)
				nEvent_Z_et++;
			if (nEle_Z_eta >= 2)
				nEvent_Z_eta++;
			if (nEle_Z_eteta >= 2)
				nEvent_Z_eteta++;
			if (nEle_Z_eteta_d0 >= 2)
				nEvent_Z_eteta_d0++;
			if (nEle_Z_eteta_d0_L >= 2)
				nEvent_Z_eteta_d0_L++;
			if (nEle_Z_eteta_d0_T >= 2)
				nEvent_Z_eteta_d0_T++;
			if (nEle_Z_eteta_d0_RT >= 2)
				nEvent_Z_eteta_d0_RT++;
			if (nEle_Z_eteta_d0_RL >= 2)
				nEvent_Z_eteta_d0_RL++;
		}//Z->ee MC


		// 3B - Find number of muons
		if (m_debug)
			cout << " Starting MUON section" << endl;

		int nGoodMu = 0;
		int nGoodIsoMu = 0;
		int nGoodNonIsoMu = 0;

		for (unsigned int i = 0; i < Nmus; ++i) {

			//Make Pt and Eta cuts (consider global muons only)
			if (mus_id_AllGlobalMuons->at(i) > 0 && mus_cm_pt->at(i) > MU_PTCUT && fabs(mus_cm_eta->at(i)) < 2.1) {

				// Correct muon d0 using BeamSpot
				float mu_d0_corrected = compute_d0("muon", i);

				//Plot Muon ID quantities
				fillHistoDataAndMC(h_muon_chi2, mus_cm_chi2->at(i) / mus_cm_ndof->at(i), this_weight); // "cm" means combinedMuon
				fillHistoDataAndMC(h_muon_d0_unCor, mus_cm_d0dum->at(i), this_weight);
				fillHistoDataAndMC(h_muon_d0, mu_d0_corrected, this_weight);
				fillHistoDataAndMC(h_muon_hits, mus_tkHits->at(i), this_weight);

				//Apply V+jets Muon ID (consider only Global Muons)
				if ((mus_cm_chi2->at(i) / mus_cm_ndof->at(i)) < 10 && fabs(mu_d0_corrected) < 0.02 && mus_tkHits->at(i) >= 11
						&& mus_hcalvetoDep->at(i) < 6.0 && // veto cone energy
						mus_ecalvetoDep->at(i) < 4.0) { //what else?

					//Store 4 vector for "good" muon and increment counters
					TLorentzVector muo(mus_cm_px->at(i), mus_cm_py->at(i), mus_cm_pz->at(i), mus_energy->at(i)); //global muon

					nGoodMu++;
					muons.push_back(muo);

					//Compute Combined RelIso for muons
					float CombRelIso = (mus_tIso->at(i) + mus_cIso->at(i)) / mus_cm_pt->at(i); //new reliso
					float CombRelIso2 = 1 / (1 + CombRelIso); //norm (old)

					// Apply isolation cut
					bool isIsolated = false;
					if (useNewReliso)
						isIsolated = CombRelIso < 0.05; //new reliso
					else
						isIsolated = CombRelIso2 > 0.85; //old reliso
					if (m_debug) {
						cout << " muon  Old RelIso = " << CombRelIso2;
						cout << "  New RelIso = " << CombRelIso;
						cout << "  isIsolated = " << isIsolated << endl;
					}

					if (isIsolated) { //Isolated

						//Store 4 vector for isolated muon and increment counters
						nGoodIsoMu++;
						iso_muons.push_back(muo);
					} else { //Non-Isolated
						nGoodNonIsoMu++;
					}

				}// end muon ID cut
			}// end pt,eta cut
		}// end muons loop


		// 3C - Set muon, Z veto, and conversion flags
		if (m_debug)
			cout << " Starting MUON Veto section" << endl;
		bool isMuon = false;
		if (nGoodIsoMu > 0)
			isMuon = true;

		if (m_debug)
			cout << " Starting CONVERSION Veto section" << endl;
		bool isConversion = false;
		if (nGoodIsoEle == 1) { // If 1 electron only, check if this electron is a conversion ... if more than one ignore (will be rejected anyawy)
			isConversion = ConversionFinder(iso_electrons.at(0), mctype, index_selected_ele);
			if (DoConversionStudies()) {
				ConversionMCMatching(iso_electrons.at(0), mctype, isConversion);
			}
			//OptimiseConversionFinder(iso_electrons.at(0), mctype);
		}

		// Is the fist isolated electron in Barrel?
		bool isBarrel = false;
		if (iso_electrons.size() > 0 && fabs(iso_electrons.at(0).Eta()) < 1.442)
			isBarrel = true;

		//-----------------------------------------------------------------------------------
		//  Z veto (normal selection)
		//-----------------------------------------------------------------------------------
		// Two components:
		//  (1) If event contain 2 or more isolated electrons, flag as Z.
		//  (2) If event has just 1 isolated electron, plus a robustLoose electron with 20 GeV,
		//      and if the di-electron invariant mass falls in window 75-106 GeV, flag as Z.
		//-----------------------------------------------------------------------------------
		if (m_debug)
			cout << " Starting Z Veto section" << endl;

		// Component (1)
		bool isZ_twoE = false;
		if (nGoodIsoEle >= 2) {
			isZ_twoE = true;
			if (nGoodIsoEle >= 3)
				cout << "*** Multilepton (>=3) Event!!! ***" << endl;
		}

		// Component (2)
		bool isZ_mee_NEW = false;

		if (nGoodIsoEle == 1 && Nels > 1) {

			//cout << "TEST: N(e) = " << Nels << endl;
			for (unsigned int j = 0; j < Nels; ++j) { //e loop

				if (j == index_selected_ele)
					continue;

				//consider only 2nd electron with ET>20 GeV, eta<2.5, d0<200um, RobustLoose
				if (els_et->at(j) < 20.0)
					continue;
				if (fabs(els_eta->at(j)) > 2.5)
					continue;
				if (fabs(compute_d0("electron", j)) > 0.02)
					continue;
				if (els_robustLooseId->at(j) < 1)
					continue;

				TLorentzVector loose(els_px->at(j), els_py->at(j), els_pz->at(j), els_energy->at(j));
				float mass = (iso_electrons.at(0) + loose).M();
				//cout << "   m(e,e) = " << mass << endl;
				fillHistoDataAndMC(h_mass_diele_new, mass, this_weight);
				if (mass > 76 && mass < 106)
					isZ_mee_NEW = true; //within window, flag

			}//e loop

		}//if more than 2 electron

		if (isZ_mee_NEW == true)
			nZ_TEST++;

		//----------------
		///1310 TEST  13-10-09
		//----------------
		/*
		 int  test_nRL = 0;
		 for ( unsigned int i=0; i<Nels; ++i) {
		 if ( els_et->at(i) > 20.0 &&
		 fabs(els_eta->at(i)) < 2.5 &&
		 els_robustLooseId->at(i) > 0 ) ++test_nRL;
		 }
		 bool isZ_twoRL = (test_nRL >= 2) ;
		 */
		//bool isZ = isZ_twoE;
		//if(test_nRL>=2) isZ = true;

		//float mass_ee = -1;

		//1310 ----------


		//-----------------------------------------------------------------------------------
		//  Z veto (normal selection)   [ Original ]
		//-----------------------------------------------------------------------------------
		// 13 Mar 09: ask for 1 robustTight electron ("good"), plus a reco-electron (not cut)
		// ie ask for 2 reco-el, and then require one to be tight (also pass et,eta,d0 cuts)
		//-----------------------------------------------------------------------------------

		if (m_debug)
			cout << " Apply Z veto" << endl;

		bool isZ_mee_ORI = false; //mee = mass(e,e)
		float mass_ee = -1;

		if (nGoodIsoEle == 1 && Nels > 1) { //better (because there are Z events with 3 ele, one is low pt (243evts/32995

			// (a) if there are only 1 RT electron, take another electron (most energetic one if there
			//     are more than one.
			if (nGoodEle == 1) {
				// check index of the RT electron
				TLorentzVector e1(electrons.at(0)); //RT
				int index = -1;
				for (unsigned int i = 0; i < Nels; ++i) { //e loop
					if (fabs(e1.Px() - els_px->at(i)) < 1e-6)
						index = i;//match
				}
				if (index < 0) {
					cout << "WARNING: cannot find index of the RobustTight electron!" << endl;
				}
				if (m_debug) {
					cout << "electrons.size = " << electrons.size() << endl;
					cout << "electrons.at(0) px py pz E = " << electrons.at(0).Px() << ", " << electrons.at(0).Py() << ", "
							<< electrons.at(0).Pz() << ", " << electrons.at(0).E() << endl;
					cout << "e1.px py pz E = " << e1.Px() << ", " << e1.Py() << ", " << e1.Pz() << ", " << e1.E() << endl;
					for (unsigned int i = 0; i < Nels; ++i) { //e loop
						cout << "   els_px(py,pz,E)->at(" << i << ") = " << els_px->at(i);
						cout << ", " << els_py->at(i) << ", " << els_pz->at(i) << ", " << els_energy->at(i) << endl;
					}
				}
				vector<TLorentzVector> loose_electrons;
				if (index >= 0) { //bug-fix Z veto (Aug09)
					//find 2nd electron passing ET>20, |eta|<2.5, d0, RobustLoose
					for (unsigned int i = 0; i < Nels; ++i) {
						if ((i - index) == 0)
							continue; //skip the Good Electron
						if (els_et->at(i) > 20 && fabs(els_eta->at(i)) < 2.5 && fabs(compute_d0("electron", i)) < 0.02
								&& els_robustLooseId->at(i) == true) {
							loose_electrons.push_back(TLorentzVector(els_px->at(i), els_py->at(i), els_pz->at(i), els_energy->at(
									i)));
						}
					}
				}
				if (loose_electrons.size() > 0) {
					float mass = (e1 + loose_electrons.at(0)).M();
					mass_ee = mass;
				}
			}

			// (b) if there are >=2 RT electrons, use them to calculate inv mass(e,e)
			//     - if there are more than 2, use the first 2.
			else if (nGoodEle >= 2) {
				float mass = (electrons.at(0) + electrons.at(1)).M();
				mass_ee = mass;
			}
		}// 2 reco-e, of which 1 is RT
		if (mass_ee > 0)
			fillHistoDataAndMC(h_mass_diele, mass_ee, this_weight);
		//			fillHistoDataAndMC(h_mass_diele, mass_ee, this_weight);
		if (mass_ee > 76 && mass_ee < 106)
			isZ_mee_ORI = true;
		//if(isZ_mee) isZ = true;

		if (isZ_mee_ORI)
			nZ_ORIG++;

		//--------------------------------------
		// 2-12-09: Revisit Z veto (alternative)
		// NB: not clear if this is better, need more MC stats to test
		//--------------------------
		/*
		 // Loop over all els and find "RL" e
		 bool isZ_mee_2RL = false;
		 vector<TLorentzVector> myRLe;

		 for (unsigned int j=0; j<Nels; ++j){ //e loop

		 //consider only ele with ET>20 GeV, eta<2.5, d0<200um, RobustLoose
		 if( els_et->at(j) < 20.0 ) continue;
		 if( fabs( els_eta->at(j) ) > 2.5 ) continue;
		 if( fabs(compute_d0("electron",j)) > 0.02 ) continue;
		 if( els_robustLooseId->at(j) < 1 ) continue;

		 TLorentzVector loose(els_px->at(j),els_py->at(j),els_pz->at(j),els_energy->at(j));
		 myRLe.push_back(loose);
		 }
		 // find all myRLe pair
		 for (unsigned int i=0; i<myRLe.size(); ++i) {
		 for (unsigned int j=0; j<myRLe.size(); ++j) {
		 if(i==j) continue; //skip pair of ele with itself
		 float mass = ( myRLe.at(i) + myRLe.at(j) ).M();
		 //cout << "   m(e,e) = " << mass << endl;
		 //fillHistoDataAndMC( h_mass_diele_new, mass, this_weight );
		 if ( mass > 76.  &&  mass < 106. ) isZ_mee_2RL = true; //within window, flag
		 }
		 }//if more than 2 electron
		 */

		// Choose which Z-veto to use
		//----------------------------
		bool isZ_mee;
		if (use_old_Z_veto)
			isZ_mee = isZ_mee_ORI; //old
		else
			isZ_mee = isZ_mee_NEW; //new

		bool isZ = isZ_twoE || isZ_mee;
		//bool isZ = isZ_twoE || isZ_mee_2RL;


		// 3D - Find number of good Primary Vertices
		// Ignore for now, PV currently not in ntuple?

		// ---->  Add PV here?  [ Aug 09 ]


		//std::cout << " Quality " << zv->quality << " Vertex z = " << zv->z_pos << " +- " << zv->z_err
		//    << " Ntrk " << zv->n_trk << " Sum_pt " << zv->sum_pt << std::endl;


		// 3E - Find number of jets
		if (m_debug)
			cout << " Starting JET section" << endl;
		int nGoodJet = 0;
		int ntj = 0;
		if (m_jetAlgo == "Default") {
			for (unsigned int i = 0; i < Njets; ++i) {
				if (jets_pt->at(i) > JET_PTCUT && fabs(jets_eta->at(i)) < 2.4 && jets_emf->at(i) > 0.01) { //<-- NEW: jets_emf
					TLorentzVector jt(jets_px->at(i), jets_py->at(i), jets_pz->at(i), jets_energy->at(i));
					nGoodJet++;//FIXME: isn't jets.size() enough?
					goodJets.push_back(jt);
					btag_information[btag_TC_highEff].push_back(jets_btag_TC_highEff->at(i));
					btag_information[btag_TC_highPur].push_back(jets_btag_TC_highPur->at(i));
					btag_information[btag_JetBProb].push_back(jets_btag_jetBProb->at(i));
					btag_information[btag_JetProb].push_back(jets_btag_jetProb->at(i));
					btag_information[btag_softEle].push_back(jets_btag_softEle->at(i));
					btag_information[btag_softMuon].push_back(jets_btag_softMuon->at(i));
					btag_information[btag_softMuonNoIP].push_back(jets_btag_softMuonNoIP->at(i));
					btag_information[btag_secondaryVertex].push_back(jets_btag_secVertex->at(i));
					btag_information[btag_fake].push_back(jets_parton_Id->at(i));
				}
			}
		} else if (m_jetAlgo == "SC5") {
			for (unsigned int i = 0; i < NjetsSC5; ++i) {
				if (jetsSC5_pt->at(i) > JET_PTCUT && fabs(jetsSC5_eta->at(i)) < 2.4 && jetsSC5_emf->at(i) > 0.01) {
					TLorentzVector jt(jetsSC5_px->at(i), jetsSC5_py->at(i), jetsSC5_pz->at(i), jetsSC5_energy->at(i));
					nGoodJet++;//FIXME: isn't jets.size() enough?
					goodJets.push_back(jt);
					btag_information[btag_TC_highEff].push_back(jetsSC5_btag_TC_highEff->at(i));
					btag_information[btag_TC_highPur].push_back(jetsSC5_btag_TC_highPur->at(i));
					btag_information[btag_JetBProb].push_back(jetsSC5_btag_jetBProb->at(i));
					btag_information[btag_JetProb].push_back(jetsSC5_btag_jetProb->at(i));
					btag_information[btag_softEle].push_back(jetsSC5_btag_softEle->at(i));
					btag_information[btag_softMuon].push_back(jetsSC5_btag_softMuon->at(i));
					btag_information[btag_softMuonNoIP].push_back(jetsSC5_btag_softMuonNoIP->at(i));
					btag_information[btag_secondaryVertex].push_back(jetsSC5_btag_secVertex->at(i));
				}
			}
		} else if (m_jetAlgo == "KT4") {
			for (unsigned int i = 0; i < NjetsKT4; ++i) {
				if (jetsKT4_pt->at(i) > JET_PTCUT && fabs(jetsKT4_eta->at(i)) < 2.4 && jetsKT4_emf->at(i) > 0.01) {
					TLorentzVector jt(jetsKT4_px->at(i), jetsKT4_py->at(i), jetsKT4_pz->at(i), jetsKT4_energy->at(i));
					nGoodJet++;//FIXME: isn't jets.size() enough?
					goodJets.push_back(jt);
					btag_information[btag_TC_highEff].push_back(jetsKT4_btag_TC_highEff->at(i));
					btag_information[btag_TC_highPur].push_back(jetsKT4_btag_TC_highPur->at(i));
					btag_information[btag_JetBProb].push_back(jetsKT4_btag_jetBProb->at(i));
					btag_information[btag_JetProb].push_back(jetsKT4_btag_jetProb->at(i));
					btag_information[btag_softEle].push_back(jetsKT4_btag_softEle->at(i));
					btag_information[btag_softMuon].push_back(jetsKT4_btag_softMuon->at(i));
					btag_information[btag_softMuonNoIP].push_back(jetsKT4_btag_softMuonNoIP->at(i));
					btag_information[btag_secondaryVertex].push_back(jetsKT4_btag_secVertex->at(i));
				}
			}
		} else if (m_jetAlgo == "pfjet") {
			for (unsigned int i = 0; i < NPFJets; ++i) {
				if (PFJets_pt->at(i) > JET_PTCUT && fabs(PFJets_eta->at(i)) < 2.4) {// NB: EMF cannot be calculated for PFJets
					TLorentzVector jt(PFJets_px->at(i), PFJets_py->at(i), PFJets_pz->at(i), PFJets_energy->at(i));
					nGoodJet++;//FIXME: isn't jets.size() enough?
					goodJets.push_back(jt);
					btag_information[btag_TC_highEff].push_back(PFJets_btag_TC_highEff->at(i));
					btag_information[btag_TC_highPur].push_back(PFJets_btag_TC_highPur->at(i));
					btag_information[btag_JetBProb].push_back(PFJets_btag_jetBProb->at(i));
					btag_information[btag_JetProb].push_back(PFJets_btag_jetProb->at(i));
					btag_information[btag_softEle].push_back(PFJets_btag_softEle->at(i));
					btag_information[btag_softMuon].push_back(PFJets_btag_softMuon->at(i));
					btag_information[btag_softMuonNoIP].push_back(PFJets_btag_softMuonNoIP->at(i));
					btag_information[btag_secondaryVertex].push_back(PFJets_btag_secVertex->at(i));
					btag_information[btag_fake].push_back(PFJets_parton_Id->at(i));
				}
			}

		} else {
			cout << "ERROR: wrong jet collection." << endl;
		}

		// 3F - Do "jet cleaning"
		if (m_debug)
			cout << " Do jet cleaning" << endl;

		float delR_jet_ele = 9999;

		//Loop over selected jets
		for (unsigned int i = 0; i < goodJets.size(); ++i) {

			//Loop over selected (isolated) electrons
			//for(unsigned int j = 0 ; j < iso_electrons.size(); ++j) {
			for (unsigned int j = 0; j < electrons.size(); ++j) {

				//For each jet, compute dR to *all* selected electrons
				delR_jet_ele = calcDeltaR(goodJets.at(i), electrons.at(j));
				//delR_jet_ele = calcDeltaR(jets.at(i),iso_electrons.at(j));

				//If delta R is within a cone of 0.3, remove this jet from selected jet list, since it's an electron
				if (delR_jet_ele < 0.3) {
					nGoodJet--;
					goodJets.erase(goodJets.begin() + i);
					i--;
					break;
				}
			}
		}// end jets loop


		// (TL) NEW: 14 Aug 09
		// basic plot: jet collection after cleaning
		//		fillHistoDataAndMC(h_njet, float(nGoodJet), this_weight);
		fillHistoDataAndMC(knjets, double(nGoodJet), this_weight);
		m_nGoodJet = nGoodJet; //set private var for this event//FIXME: why no use this from the very beginning?


		if (Njet() < 0) {
			cout << "ERROR!!! negative njet: " << Njet() << endl;
		}

		for (unsigned int i = 0; i < goodJets.size(); ++i) {
			fillHistoDataAndMC(kjet_pt_all, goodJets.at(i).Pt(), this_weight);
			fillHistoDataAndMC(kjet_eta_all, goodJets.at(i).Eta(), this_weight);
			fillHistoDataAndMC(kjet_phi_all, goodJets.at(i).Phi(), this_weight);
			if (i < 4) {// fill for first 4 jets
				fillHistoDataAndMC(kjet_pt_all + i + 1, goodJets.at(i).Pt(), this_weight);
				fillHistoDataAndMC(kjet_eta_all + i + 1, goodJets.at(i).Eta(), this_weight);
				fillHistoDataAndMC(kjet_phi_all + i + 1, goodJets.at(i).Phi(), this_weight);
			}
		}

		// fill electron-counting histo
		fillHisto_Njet_DataAndMC(h_nEle_all, Nels, 1.0);
		fillHisto_Njet_DataAndMC(h_nEle_all, Nels, 1.0);
		fillHisto_Njet_DataAndMC(h_nEle_s1, nBasicEle, 1.0); //pass et,eta
		fillHisto_Njet_DataAndMC(h_nEle_s2, nBasicD0Ele, 1.0); //also pass d0
		fillHisto_Njet_DataAndMC(h_nEle_s3_idLoose, nLooseEle, 1.0); //also pass eidL
		fillHisto_Njet_DataAndMC(h_nEle_s3_idTight, nTightEle, 1.0); //also pass eidT
		fillHisto_Njet_DataAndMC(h_nEle_s3_idRL, nRobustLooseEle, 1.0); //also pass eidRL
		fillHisto_Njet_DataAndMC(h_nEle_s3_idRT, nGoodEle, 1.0); //also pass eidRT

		// 3G - Order jets in Et

		// IGNORE FOR NOW - not necessary

		// 3H - Check z of lepton within to z of PV (within 2cm)?

		// IGNORE FOR NOW - no PV info in ntuples

		bool isDifferentInteraction = false;

		// 4 - Find numbers of taggable and tagged jets

		// IGNORE FOR NOW - no tagging info in ntuples

		int ntaggable = 0;
		int nbtagP = 0;
		int nbtagN = 0;

		// 5 Get MET

		if (m_debug)
			cout << " Starting MET section" << endl;

		if (Nmets > 1)
			cout << "WARNING: more than 1 met object, using index = 0 one" << endl;

		// NOTE, filling met4v (MET_x,MET_y,0, sqrt(MET_x^2 + MET_y^2))
		// e.g. to get |met| = met4v.Pt();
		// e.g to get met_phi, do arccos(met4v.px()/met4v.E()))

		//TLorentzVector met4v(mets_ex->at(0),mets_ey->at(0),0,mets_et->at(0));
		//met.push_back(met4v);
		//const float this_met = met.at(0).Pt();

		// Use muon-corrected MET
		//const float this_met     = mets_et_muonCor->at(0);
		//const float this_met_phi = mets_phi_muonCor->at(0);
		float this_met = 0;
		float this_mex = 0;
		float this_mey = 0;
		float this_met_phi = 0;

		if (m_metAlgo == "Default") {
			this_met = mets_et_muonCor->at(0);
			this_met_phi = mets_phi_muonCor->at(0);
			this_mex = mets_ex->at(0);
			this_mey = mets_ey->at(0);

		} else if (m_metAlgo == "calomet_mujes") {
			this_met = mets_et->at(0);
			this_met_phi = mets_phi->at(0);
			this_mex = mets_ex->at(0);
			this_mey = mets_ey->at(0);

		} else if (m_metAlgo == "SC5") {
			this_met = metsSC5_et->at(0);
			this_met_phi = metsSC5_phi->at(0);
			this_mex = metsSC5_ex->at(0);
			this_mey = metsSC5_ey->at(0);

		} else if (m_metAlgo == "KT4") {
			this_met = metsKT4_et->at(0);
			this_met_phi = metsKT4_phi->at(0);
			this_mex = metsKT4_ex->at(0);
			this_mey = metsKT4_ey->at(0);

		} else if (m_metAlgo == "tcmet") {
			this_met = tcmets_et->at(0);
			this_met_phi = tcmets_phi->at(0);
			this_mex = tcmets_ex->at(0);
			this_mey = tcmets_ey->at(0);

		} else if (m_metAlgo == "pfmet") {
			this_met = PFMets_et->at(0);
			this_met_phi = PFMets_phi->at(0);
			this_mex = PFMets_ex->at(0);
			this_mey = PFMets_ey->at(0);
		}

		//FIXME: SET new met:
		TVector2 met2v_mu;
		TVector2 met2v_t1;
		TLorentzVector tmpmet;
		//Corrected filling FIXED MET!
		tmpmet. SetPxPyPzE(this_mex, this_mey, 0, this_met);
		met.push_back(tmpmet);
		met2v_mu.SetMagPhi(mets_et_muonCor->at(0), mets_phi_muonCor->at(0));
		met2v_t1.SetMagPhi(mets_et->at(0), mets_phi->at(0));

		fillHistoDataAndMC(h_metAlone, this_met, this_weight);
		fillHistoDataAndMC(h_metAlone_phi, this_met_phi, this_weight);

		// 6 - Calculate some kinematic variables

		if (m_debug)
			cout << " Starting HT calculation section" << endl;

		// HT =scalar sum of MET, electron PT, muon PT and jet PT
		float ht = this_met;
		for (unsigned int i = 0; i < iso_electrons.size(); ++i) {
			ht += iso_electrons.at(i).Pt();
		}
		for (unsigned int i = 0; i < iso_muons.size(); ++i) {
			ht += iso_muons.at(i).Pt();
		}
		for (unsigned int i = 0; i < goodJets.size(); ++i) {
			ht += goodJets.at(i).Pt();
		}

		// compute mT(W)
		//     float this_mtw = -1;
		float this_mu_mtw = -1;
		float this_t1_mtw = -1;
		float this_mu_DPhiEmet = -1;
		float this_t1_DPhiEmet = -1;

		if (iso_electrons.size() > 0) {

			TVector2 e2v(iso_electrons.at(0).Px(), iso_electrons.at(0).Py());

			this_mu_mtw = compute_mtw(e2v, met2v_mu); // muon-corr. caloMET
			this_t1_mtw = compute_mtw(e2v, met2v_t1); // type 1 caloMET
			fillHistoDataAndMC(h_mtw_mu_incl, this_mu_mtw, this_weight);
			fillHistoDataAndMC(h_mtw_t1_incl, this_t1_mtw, this_weight);

			/*
			 cout << "\n\nele:     px = " << e2v.X()  << ",  py = " << e2v.Y() << ",  phi = " << e2v.Phi()<< endl;
			 cout << "miss (mu): px = " << met2v_mu.X() << ",  py = " << met2v_mu.Y() << ",  phi = " << met2v_mu.Phi() << endl;
			 float dphi = e2v.Phi() - met2v_mu.Phi() ;
			 cout << "phi(ele) - phi(miss) = " << dphi << endl;
			 if(fabs(dphi) > TMath::Pi() ) dphi = 2*TMath::Pi() - fabs(dphi);
			 cout << " correct = " << dphi << endl;
			 */
			//       TVector2 eee(iso_electrons.at(0).P() ;

			this_mu_DPhiEmet = (e2v).DeltaPhi(met2v_mu);
			this_t1_DPhiEmet = (e2v).DeltaPhi(met2v_t1);
			/*
			 cout  << "\nthis_mu_DPhiEmet = " << this_mu_DPhiEmet << endl;
			 cout  << " ele phi = " << iso_electrons.at(0).Phi() << endl;
			 cout  << " miss phi (mu) = "<< met2v_mu.Phi() << endl;
			 cout  << " miss phi (t1) = "<< met2v_t1.Phi() << endl;
			 */
			fillHistoDataAndMC(h_DPhiEmet_mu_incl, this_mu_DPhiEmet, this_weight);
			fillHistoDataAndMC(h_DPhiEmet_t1_incl, this_t1_DPhiEmet, this_weight);

		}

		//*******************************************************************************************
		//*******************************************************************************************
		// Now select events based on quantities filled above
		//*******************************************************************************************
		//*******************************************************************************************


		//Event counts
		//stage 0 - initial count from files
		//stage 1 - good run checked
		//stage 2 - trigger checked
		//stage 3 - >= 1 lepton passing cuts
		//stage 4 - == 1 lepton passing cuts
		//stage 5 - pass muon veto
		//stage 6 - pass MET cut
		//stage 7 - pass Z veto
		//stage 8 - not a conversion
		//stage 9 - letpton close to PV
		//stage 10 - pass HT cut
		//stage 11 - >= 1 taggable jet
		//stage 12 - >= 1+ tagged jet
		//stage 13 - >= 2+ tagged jet
		//stage 14 - >= 1- tagged jet

		//Fill counters
		ntj = nGoodJet;
		if (nGoodJet > 4)
			ntj = 4;

		bool e_plus_jet_pass = false;

		if (m_debug) {
			cout << "this: weight: " << this_weight;
			//cout << "wmap[this_mc]: "<< weightMap[this_mc]
			cout << " mc: " << mclabel[fastmctype_];
			//cout << " type: " << mctype << endl;
			cout << endl;
		}

		//------------------------------------------------------------------------------------
		//81FB  produce validation plots
		//-------------------------------
		if (doValidation()) {
			if (m_debug)
				cout << " Produce validation plots" << endl;
			const bool Boolcuts[9] = { 1, fired_single_em, (nGoodEle > 0), (nGoodIsoEle > 0), (nGoodIsoEle == 1), !isMuon,
					(this_met > METCUT), !isZ, !isConversion };

			valid_fillHisto(valid_HT, Boolcuts, nGoodJet, ht);

			for (unsigned int i = 0; i < goodJets.size(); ++i) {
				valid_fillHisto(valid_jetsEt, Boolcuts, nGoodJet, goodJets.at(i).Et());
				valid_fillHisto(valid_jetsEta, Boolcuts, nGoodJet, goodJets.at(i).Phi());
				valid_fillHisto(valid_jetsPhi, Boolcuts, nGoodJet, goodJets.at(i).Eta());
			}

			if (ntj > 0)
				valid_fillHisto(valid_jets1stEt, Boolcuts, nGoodJet, goodJets.at(0).Et());
			if (ntj > 1)
				valid_fillHisto(valid_jets2ndEt, Boolcuts, nGoodJet, goodJets.at(1).Et());
			if (ntj > 2)
				valid_fillHisto(valid_jets3rdEt, Boolcuts, nGoodJet, goodJets.at(2).Et());
			if (ntj > 3)
				valid_fillHisto(valid_jets4thEt, Boolcuts, nGoodJet, goodJets.at(3).Et());

			vector<TLorentzVector> valid_eles;
			for (unsigned int i = 0; i < Nels; ++i) {
				bool eleBoolcuts[9] = { 1, fired_single_em, (nGoodEle > 0), (nGoodIsoEle > 0), (nGoodIsoEle == 1), !isMuon,
						(this_met > METCUT), !isZ, !isConversion };

				float d0_corrected = fabs(compute_d0("electron", i)); //abs
				float RelCalIso = (els_dr04EcalRecHitSumEt->at(i) + els_dr04HcalTowerSumEt->at(i)) / els_et->at(i);
				float RelTrkIso = els_tIso->at(i) / els_et->at(i);

				if ((RelCalIso + RelTrkIso) > 0.1) {
					eleBoolcuts[3] = 0;
				}
				if (els_et->at(i) < ELE_ETCUT || d0_corrected > 0.02) {
					eleBoolcuts[2] = 0;
				}
				//if(els_robustTightId->at(i) == false || fabs(els_eta->at(i))>2.5 ) { eleBoolcuts[2] = 0; }
				if (passEleID(i) == false || fabs(els_eta->at(i)) > 2.5) {
					eleBoolcuts[2] = 0;
				}

				valid_fillHisto(valid_eleEt, eleBoolcuts, nGoodJet, els_et->at(i));
				valid_fillHisto(valid_eleEta, eleBoolcuts, nGoodJet, els_eta->at(i));
				valid_fillHisto(valid_elePhi, eleBoolcuts, nGoodJet, els_phi->at(i));
				valid_fillHisto(valid_eled0, eleBoolcuts, nGoodJet, d0_corrected);
				valid_fillHisto(valid_eleCalIso, eleBoolcuts, nGoodJet, RelCalIso);
				valid_fillHisto(valid_eleTrkIso, eleBoolcuts, nGoodJet, RelTrkIso);
				valid_fillHisto(valid_eleRelIso, eleBoolcuts, nGoodJet, (RelCalIso + RelTrkIso));
				TLorentzVector eles(els_px->at(i), els_py->at(i), els_pz->at(i), els_energy->at(i));
				valid_eles.push_back(eles);
			}//end of Nels loop

			if (nGoodEle > 0) {
				for (size_t i = 0; i < valid_eles.size(); ++i) {
					double mass_mee = (electrons.at(0) + valid_eles.at(i)).M();
					valid_fillHisto(valid_mass_ee, Boolcuts, nGoodJet, mass_mee);
				}
			}

			valid_fillHisto(valid_metEt, Boolcuts, nGoodJet, this_met);
			valid_fillHisto(valid_metPhi, Boolcuts, nGoodJet, this_met_phi);

			valid_fillHisto(valid_numberTracks, Boolcuts, nGoodJet, Ntracks);
			for (unsigned int i = 0; i < Ntracks; ++i) {
				valid_fillHisto(valid_trackPt, Boolcuts, nGoodJet, tracks_pt->at(i));
			}

			if (nGoodJet > 2) {
				pair<double, double> res = compute_M3(goodJets);
				valid_fillHisto(valid_recoM3, Boolcuts, nGoodJet, res.first);
				valid_fillHisto(valid_recoM3_PTMax, Boolcuts, nGoodJet, res.second);
			}

			if (!IsData() && isMCType(kttbar)) { //check some gen. quantities for top quark
				double genttbarpt = 0;
				double genttbarpx = 0;
				double genttbarpy = 0;
				for (unsigned int g = 0; g < Nmc_doc; g++) {
					if (fabs(mc_doc_id->at(g)) == 6) {
						valid_fillHisto(valid_genT_pt, Boolcuts, nGoodJet, mc_doc_pt->at(g));
						genttbarpx += mc_doc_px->at(g);
						genttbarpy += mc_doc_py->at(g);
					}
				}
				genttbarpt = sqrt(genttbarpx * genttbarpx + genttbarpy * genttbarpy);
				valid_fillHisto(valid_genTT_pt, Boolcuts, nGoodJet, genttbarpt);
			}
		}//end doValidation
		//81FB end
		//------------------------------- End validation ---------------------------


		e_plus_jet[0][ntj][mctype]++;
		e_plus_jet_weighted[0][ntj][mctype] += this_weight; //1 event x weight


		if (m_debug)
			cout << " Applying event selection" << endl;

		if (goodrun) { //Event in GOOD run list

			e_plus_jet[1][ntj][mctype]++;
			e_plus_jet_weighted[1][ntj][mctype] += this_weight;

			if (fired_single_em) { //Trigger
				e_plus_jet[2][ntj][mctype]++;
				e_plus_jet_weighted[2][ntj][mctype] += this_weight;

				// Electron checks
				if (nGoodEle > 0) {
					e_plus_jet[3][ntj][mctype]++;
					e_plus_jet_weighted[3][ntj][mctype] += this_weight;

					// Isolated electron
					if (nGoodIsoEle > 0) {
						e_plus_jet[4][ntj][mctype]++;
						e_plus_jet_weighted[4][ntj][mctype] += this_weight;

						if (nGoodIsoEle == 1) {
							e_plus_jet[5][ntj][mctype]++;
							e_plus_jet_weighted[5][ntj][mctype] += this_weight;

							if (!isMuon) { // Muon Veto
								e_plus_jet[6][ntj][mctype]++;
								e_plus_jet_weighted[6][ntj][mctype] += this_weight;

								if (this_met > METCUT) { // MET
									e_plus_jet[7][ntj][mctype]++;
									e_plus_jet_weighted[7][ntj][mctype] += this_weight;

									if (!isZ) { // Z Veto
										e_plus_jet[8][ntj][mctype]++;
										e_plus_jet_weighted[8][ntj][mctype] += this_weight;

										if (!isConversion) { //Conversion Veto
											e_plus_jet[9][ntj][mctype]++;
											e_plus_jet_weighted[9][ntj][mctype] += this_weight;

											if ((m_rejectEndcapEle == false && !isDifferentInteraction) || // PV check (DIFFZ)
													(m_rejectEndcapEle == true && isBarrel)) { // ele eta cut

												e_plus_jet[10][ntj][mctype]++;
												e_plus_jet_weighted[10][ntj][mctype] += this_weight;

												if (ht >= HTCUT) { // HT cut (30-1-09)
													e_plus_jet[11][ntj][mctype]++;
													e_plus_jet_weighted[11][ntj][mctype] += this_weight;

													e_plus_jet_pass = true;

													if (ntaggable > 0) { // taggable
														e_plus_jet[12][ntj][mctype]++;
														e_plus_jet_weighted[12][ntj][mctype] += this_weight;

														if (nbtagP >= 1) { //at least one +tag
															e_plus_jet[13][ntj][mctype]++;
															e_plus_jet_weighted[13][ntj][mctype] += this_weight;

															if (nbtagP >= 2) { // at least two +tags
																e_plus_jet[14][ntj][mctype]++;
																e_plus_jet_weighted[14][ntj][mctype] += this_weight;
															}
														}
														if (nbtagN >= 1) { //at least one -tag
															e_plus_jet[15][ntj][mctype]++;
															e_plus_jet_weighted[15][ntj][mctype] += this_weight;
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//+++++++++++++++++++++++++  Analysis!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		if (m_debug)
			cout << " Starting << ANALYSIS >>" << endl;

		if (goodrun && fired_single_em && nGoodEle > 0) {
			//			fasthist_[fastmctype_][h_met_ante_ISO]->Fill(this_met, this_weight); //user-chosen MET
			//			fasthist_[fastmctype_][h_met_ante_ISO_mu]->Fill(met2v_mu.Mod(), this_weight);
			//			fasthist_[fastmctype_][h_met_ante_ISO_t1]->Fill(met2v_t1.Mod(), this_weight);
			fillHistoDataAndMC(h_met_ante_ISO, this_met, this_weight); //user-chosen MET
			fillHistoDataAndMC(h_met_ante_ISO_mu, met2v_mu.Mod(), this_weight);
			fillHistoDataAndMC(h_met_ante_ISO_t1, met2v_t1.Mod(), this_weight);
		}

		//E+jets Analysis
		if (e_plus_jet_pass && nGoodJet >= 4) {

			++counter_pass;

			//Sample hist
			pass_ht->Fill(ht, this_weight);
			pass_met->Fill(this_met, this_weight);
			fillHistoDataAndMC(h_ed0_pass, this_isoele_d0, this_weight);
			//			fillHistoDataAndMC(h_ed0_pass, this_isoele_d0, this_weight);

			//Print out for each selected event
			fprintf(outfile, "%7d %6d %10d %10d %8d %8lld %5dj %11.2f %8.2f %8.2f %8s   %s\n", counter_pass, run, event,
					lumiBlock, chain->GetTreeNumber(), lflag, nGoodJet, goodJets.at(0).Pt(), this_met, ht,
					mclabel[fastmctype_].c_str(), chain->GetCurrentFile()->GetName());
			interestingFiles.insert(chain->GetCurrentFile()->GetName());
		}

		// (21 Feb 09) make some kinematics plots for events passing N-1 cuts (HT,MET)
		if (m_debug)
			cout << "[M_DEBUG] Filling N-1 histograms" << endl;
		if (goodrun && fired_single_em && nGoodIsoEle == 1 && !isMuon && !isZ && !isConversion && !isDifferentInteraction) {

			if (this_met > METCUT) //after all but HT cut
				fillHisto_Njet_DataAndMC(h_ht, ht, this_weight);

			if (ht >= HTCUT) { //after all but MET cut
				fillHisto_Njet_DataAndMC(h_met, this_met, this_weight); //user-chosen MET
				fillHisto_Njet_DataAndMC(h_met_mu, met2v_mu.Mod(), this_weight);
				fillHisto_Njet_DataAndMC(h_met_t1, met2v_t1.Mod(), this_weight);
				if (isBarrel) { //26-1-10
					fillHisto_Njet_DataAndMC(h_met_BA, this_met, this_weight);
					fillHisto_Njet_DataAndMC(h_met_mu_BA, met2v_mu.Mod(), this_weight);
					fillHisto_Njet_DataAndMC(h_met_t1_BA, met2v_t1.Mod(), this_weight);
				}
			}
			// after all but MET/HT/Njet cuts
			fillHisto_Njet_DataAndMC(h_mtw_mu, this_mu_mtw, this_weight);
			fillHisto_Njet_DataAndMC(h_mtw_t1, this_t1_mtw, this_weight);
			fillHisto_Njet_DataAndMC(h_DPhiEmet_mu, this_mu_DPhiEmet, this_weight);
			fillHisto_Njet_DataAndMC(h_DPhiEmet_t1, this_t1_DPhiEmet, this_weight);

			// inspect distributions of the selected events
			if (Njet() >= 4) {
				//				fasthist_[fastmctype_][h_exp_ele_et]->Fill(iso_electrons.at(0).Et(), this_weight);
				//				fasthist_[fastmctype_][h_exp_ele_eta]->Fill(iso_electrons.at(0).Eta(), this_weight);
				//				fasthist_[fastmctype_][h_exp_j0_pt]->Fill(jets.at(0).Pt(), this_weight);
				//				fasthist_[fastmctype_][h_exp_j1_pt]->Fill(jets.at(1).Pt(), this_weight);
				//				fasthist_[fastmctype_][h_exp_DRej]->Fill(iso_electrons.at(0).DeltaR(jets.at(0)), this_weight);
				//				fasthist_[fastmctype_][h_exp_DPhiej]->Fill(iso_electrons.at(0).DeltaPhi(jets.at(0)), this_weight);
				//				fasthist_[fastmctype_][h_exp_DRjj]->Fill(jets.at(0).DeltaR(jets.at(1)), this_weight);
				//				fasthist_[fastmctype_][h_exp_DPhijj]->Fill(jets.at(0).DeltaPhi(jets.at(1)), this_weight);
				//				fasthist2D_[fastmctype_][h_exp_met_v_eeta]->Fill(this_met, iso_electrons.at(0).Eta(), this_weight);
				fillHistoDataAndMC(h_exp_ele_et, iso_electrons.at(0).Et(), this_weight);
				fillHistoDataAndMC(h_exp_ele_eta, iso_electrons.at(0).Eta(), this_weight);
				fillHistoDataAndMC(h_exp_j0_pt, goodJets.at(0).Pt(), this_weight);
				fillHistoDataAndMC(h_exp_j1_pt, goodJets.at(1).Pt(), this_weight);
				fillHistoDataAndMC(h_exp_DRej, iso_electrons.at(0).DeltaR(goodJets.at(0)), this_weight);
				fillHistoDataAndMC(h_exp_DPhiej, iso_electrons.at(0).DeltaPhi(goodJets.at(0)), this_weight);
				fillHistoDataAndMC(h_exp_DRjj, goodJets.at(0).DeltaR(goodJets.at(1)), this_weight);
				fillHistoDataAndMC(h_exp_DPhijj, goodJets.at(0).DeltaPhi(goodJets.at(1)), this_weight);
				fillHistoDataAndMC(h_exp_met_v_eeta, this_met, iso_electrons.at(0).Eta(), this_weight);
			}
		}
		if (m_debug)
			cout << "[M_DEBUG] After N-1" << endl;

		//----------------------------------------------------
		//  Make isolation plots needed by QCD Estimation
		// TL: 1D method: plot Reliso distribution for events passing all cuts except electron isolation and nj
		// NOTE: at the moment still using the old definition of RelIso = et/(et+SumIso), to be changed.
		// notes: all cuts except ISO and Njet
		//cout << "200) isZ = " << isZ << endl;

		if (m_debug)
			cout << " Starting << Making isolation plot >>" << endl;

		float CombRelIso = -1.0;
		int ii_GoodEle_mostIso = -1;

		//cout << "(B4) ii_GoodEle_mostIso = " << ii_GoodEle_mostIso << endl;

		if (nGoodEle > 0) {
			if (nGoodEle == 1) {
				// need isolation values of the GoodEle but not available!!
				// can resolve this by creating a companion vector to store the CombRelIso values of the GoodEles
				//  <electrons> ++ <electrons_isoval>
				CombRelIso = electrons_isoval.at(0);
				//NormCombRelIso = electrons_isoval2.at(0);
				ii_GoodEle_mostIso = ii_electrons.at(0);
			} else {
				// pick the most isolated one, loop over GoodEle, because if we happen to have 2 GoodEle,
				// and one pass isolation cut and the other fails, then if we choose randomly we may choose the
				// non-isolated one, and contradict our event selection.
				// bug fix (17 Feb 09)
				float mostIso = 9999999.9; //most-isolated = smallest Reliso (not normalized)
				float mostIso2 = -9999999.9; //most-isolated = largest Reliso (normalized)
				for (int i = 0; i < nGoodEle; ++i) {
					float thisIso = electrons_isoval.at(i);
					float thisIso2 = electrons_isoval2.at(i);
					if (thisIso < mostIso) { //RelIso
						mostIso = thisIso;
						ii_GoodEle_mostIso = ii_electrons.at(i);
					}
					if (thisIso2 > mostIso2) {
						mostIso2 = thisIso2;
					} //NormReliso
				}
				CombRelIso = mostIso;
				//NormCombRelIso = mostIso2;
			}

			// To identify the index of the most isolated GoodEle in the Nels collection
			//
			//        for(unsigned int i=0; i<Nels; ++i){
			// 	 float tempComIso = (els_tIso->at(i) + els_dr04EcalRecHitSumEt->at(i) + els_dr04HcalTowerSumEt->at(i))/els_et->at(i);
			// 	 if( (tempComIso - CombRelIso)/CombRelIso  < 1e-5 ){ ii_GoodEle_mostIso = i; break;}
			// 	 //changed 280909 to make relative comparison due to diff of 1.01e-6 for one event.e-5 = 0.001%
			//        }

			if (ii_GoodEle_mostIso < 0) { //bugfix: 9-9-09
				cout << "\nERROR: PROBLEM, COULDNT FIND MOST_ISO_ELECTRON AGAIN" << endl;
				cout << "event number is" << ev << endl;
				cout << " index of most isolated GoodEle (ii_GoodEle_mostIso): " << ii_GoodEle_mostIso << endl;
				cout << " Nels: " << Nels << ",  nGoodEle: " << nGoodEle << endl;
				printf(" smallest isolation: %12.10f\n", CombRelIso);
				cout << " Isolation of all electrons: " << endl;
				for (unsigned int j = 0; j < Nels; ++j) {
					//float tmpIso = (els_tIso->at(j) + els_dr04EcalRecHitSumEt->at(j) + els_dr04HcalTowerSumEt->at(j))/els_et->at(j);
					float tmpIso = getRelIso(j);
					printf("ele %d : %12.10f\n", j, tmpIso);
				}
			}

		}//nGoodEle > 0

		//if electron is not isolated, apply conversion algo on it so that all electrons in the distribution
		//are treated the same wrt conversion algo
		//987

		if (m_debug)
			cout << " Apply << Conversion algo >> on all 'most-iso' good electrons" << endl;

		if (CombRelIso > 0.1 && nGoodEle > 0) {
			//cout <<"ii_GoodEle_mostIso = "<< ii_GoodEle_mostIso  << endl;
			TLorentzVector eles_temp(els_px->at(ii_GoodEle_mostIso), els_py->at(ii_GoodEle_mostIso), els_pz->at(
					ii_GoodEle_mostIso), els_energy->at(ii_GoodEle_mostIso));
			isConversion = ConversionFinder(eles_temp, mctype, ii_GoodEle_mostIso);
		}

		//-------------------------------------------------------------------
		//
		//  QCD estimation: Normal & Anti Event Selection
		//
		//-------------------------------------------------------------------
		if (m_debug)
			cout << " Starting << QCD estimation NES AES >>" << endl;

		/// bool isZ_mee_AES = false;
		/// bool isZ_mep_AES = false;
		float mass_ep = -1; //NEW

		if (goodrun && fired_single_em && nGoodEle > 0 && !isMuon && !isZ && !isConversion && !isDifferentInteraction && ht
				>= HTCUT) { // <-- HT cut


			// Apply missing ET cut (Normal Selection)
			bool passALL = true;
			if (m_applyMETcut && this_met < METCUT)
				passALL = false;

			// Apply eta cut if specified (optional)
			if (m_rejectEndcapEle && fabs(els_eta->at(ii_GoodEle_mostIso)) > 1.442)
				passALL = false;

			if (passALL) {
				if (m_debug)
					cout << "QCDest: event pass all cuts (except isol and nj), make isolation plots" << endl;
				fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso, CombRelIso, this_weight);
			}

			//=====================================
			// Tighter Z veto in Anti Selection
			//================= 8-6-09 ==========
			// move the following block to outside MET cut
			//=====================================
			// - M(e,e): require 1 GoodEle (no isolation requirement), and 1 "robustloose" electron
			// - M(e,e) veto window: 65-106 GeV
			// - M(e,pho): take the 1st GoodEle, and the 1st photon (no cut)
			// - M(e,pho) veto window: 40-110 GeV
			//=====================================


			//------------------------
			// 8-6-09: Define flags for AES
			bool flag_AES_pass_metcut = this_met < AES_MET_cut;
			bool flag_AES_pass_HTcut = false;
			bool flag_AES_pass_tighterZveto_mee = true;
			bool flag_AES_pass_tighterZveto_mep = true;
			//------------------------


			//1510  Simple Z veto: exclude event with 2 Nels from AES
			//       if(AES_useSimpleZveto_TwoRL) flag_AES_pass_tighterZveto_mee = !isZ_twoRL;
			if (AES_useSimpleZveto) {
				if (Nels >= 2)
					flag_AES_pass_tighterZveto_mee = false; //exclude
			}

			//======================================
			//  M(e,e): di-electron inv mass [AES]
			//======================================

			if (AES_useSimpleZveto == false && nGoodEle > 0 && Nels > 1) {

				flag_AES_pass_tighterZveto_mee = false;

				//cout << "TEST: N(e) = " << Nels << endl;
				for (unsigned int i = 0; i < Nels; ++i) { //1st e loop

					// Require 1 GoodEle
					if (els_et->at(i) < 30.0)
						continue;
					if (fabs(els_eta->at(i)) > 2.5)
						continue;
					if (fabs(els_eta->at(i)) > 1.442 && fabs(els_eta->at(i)) < 1.56)
						continue; //ignore gap
					if (fabs(compute_d0("electron", i)) > 0.02)
						continue;
					if (passEleID(i) == false)
						continue;

					for (unsigned int j = 0; j < Nels; ++j) { //2nd e loop
						if (j == i)
							continue;

						// Require "loose" electron: consider only electrons with ET>20 GeV, eta<2.5, d0<200um, RobustLoose
						if (els_et->at(j) < 20.0)
							continue;
						if (fabs(els_eta->at(j)) > 2.5)
							continue;
						if (fabs(compute_d0("electron", j)) > 0.02)
							continue; //d0 cut
						if (els_robustLooseId->at(j) < 1)
							continue;

						TLorentzVector tight(els_px->at(i), els_py->at(i), els_pz->at(i), els_energy->at(i));
						TLorentzVector loose(els_px->at(j), els_py->at(j), els_pz->at(j), els_energy->at(j));

						float mass = (tight + loose).M();
						// if fall in the window, flag as Z
						if (mass < 65. || mass > 106.)
							flag_AES_pass_tighterZveto_mee = true; // not Z, keep event

					}//2nd e loop
				}//1st e loop
				//cout << "isZ_mee_TEST = " << isZ_mee_TEST << endl;
			}//if more than 2 electron


			//=============================================
			//  M(e,pho): ele-photon invariant mass [AES]
			//=============================================
			// revise this...
			// at the moment, it only consider the inv mass of the first goodEle and the first photon

			if (nGoodEle >= 1 && Nphotons >= 1) {

				flag_AES_pass_tighterZveto_mep = false;
				TLorentzVector photon1(photons_px->at(0), photons_py->at(0), photons_pz->at(0), photons_energy->at(0));
				mass_ep = (electrons.at(0) + photon1).M(); //GoodEle + photon
				if (mass_ep < 40. || mass_ep > 110.)
					flag_AES_pass_tighterZveto_mep = true; //outside window, keep event

			}//require 1 GoodEle, 1 photon

			// tighter Z veto for AES
			//// if(!isZ_mee_AES) flag_AES_pass_tighterZveto_mee = true;
			//// if(!isZ_mep_AES) flag_AES_pass_tighterZveto_mep = true;
			//================= 8-6-09 ==========


			// Fill histograms for Z veto (AES), lowMET, 1j (already pass the normal Z veto)

			if (flag_AES_pass_metcut && Njet() == 1) {

				if (mass_ee > 0)
					//					fasthist_[fastmctype_][h_mass_diele_lowMet_1j]->Fill(mass_ee, this_weight); //plot m(ee) before cut
					fillHistoDataAndMC(h_mass_diele_lowMet_1j, mass_ee, this_weight); //plot m(ee) before cut
				if (flag_AES_pass_tighterZveto_mee) { //event survive mee cut
					if (mass_ep > 0)
						//						fasthist_[fastmctype_][h_mass_ephoton_lowMet_1j]->Fill(mass_ep, this_weight); //plot m(ep) before cut
						fillHistoDataAndMC(h_mass_ephoton_lowMet_1j, mass_ep, this_weight); //plot m(ep) before cut
				}
				// record number of ele/photon in this event
				//				fasthist_[fastmctype_][h_Nele_lowMet_1j]->Fill(Nels, this_weight);
				//				fasthist_[fastmctype_][h_Nphoton_lowMet_1j]->Fill(Nphotons, this_weight);
				fillHistoDataAndMC(h_Nele_lowMet_1j, Nels, this_weight);
				fillHistoDataAndMC(h_Nphoton_lowMet_1j, Nphotons, this_weight);
				for (unsigned int i = 0; i < Nphotons; ++i) {
					//					fasthist_[fastmctype_][h_photon_eta_lowMet_1j]->Fill(photons_eta->at(i), this_weight);
					//					fasthist_[fastmctype_][h_photon_et_lowMet_1j]->Fill(photons_et->at(i), this_weight);
					fillHistoDataAndMC(h_photon_eta_lowMet_1j, photons_eta->at(i), this_weight);
					fillHistoDataAndMC(h_photon_et_lowMet_1j, photons_et->at(i), this_weight);
					if (i == 0) {
						//						fasthist_[fastmctype_][h_photon1_eta_lowMet_1j]->Fill(photons_eta->at(i), this_weight);
						//						fasthist_[fastmctype_][h_photon1_et_lowMet_1j]->Fill(photons_et->at(i), this_weight);
						fillHistoDataAndMC(h_photon1_eta_lowMet_1j, photons_eta->at(i), this_weight);
						fillHistoDataAndMC(h_photon1_et_lowMet_1j, photons_et->at(i), this_weight);
					}
				}
			}

			//===============================
			// AES:  HT cut (HT<100)
			//===============================

			int ii_AES = ii_GoodEle_mostIso; // 15-6-09, move to front

			float ht_AES = this_met;
			ht_AES += els_et->at(ii_AES);
			for (unsigned int i = 0; i < goodJets.size(); ++i) {
				ht_AES += goodJets.at(i).Pt();
			}

			if (ht_AES < AES_HT_cut)
				flag_AES_pass_HTcut = true;
			//===============================


			//987
			//----------------------------------------------
			// 8-6-09: fill histo for n-1 AES reliso plots
			//----------------------------------------------
			// Before AES cut
			fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_before, CombRelIso, this_weight);

			// AES: all cut
			if (flag_AES_pass_metcut && flag_AES_pass_HTcut && flag_AES_pass_tighterZveto_mee && flag_AES_pass_tighterZveto_mep) {

				fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES, CombRelIso, this_weight); //All AES cut
				//cout << "AES pass all" << endl;
			}

			// AES: N-1 (met)
			if (flag_AES_pass_HTcut) {

				fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_justHT, CombRelIso, this_weight);

				if (flag_AES_pass_tighterZveto_mee && flag_AES_pass_tighterZveto_mep) {
					//cout << "amtb AES - met" << endl;
					fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_minusMET, CombRelIso, this_weight);
				}
			}
			// AES: N-1 (HT)
			if (flag_AES_pass_metcut) {

				fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_justMET, CombRelIso, this_weight);

				if (flag_AES_pass_tighterZveto_mee && flag_AES_pass_tighterZveto_mep) {
					//cout << "amtb AES - ht" << endl;
					fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_minusHT, CombRelIso, this_weight);
				}

				// AES: N-1 (tighter Z, both mee, mep)
				if (flag_AES_pass_HTcut) {
					//cout << "amtb AES - Z" << endl;
					fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_minusTighterZ, CombRelIso, this_weight);
				}
			}//pass AES_MET
			else { //event has high MET
				fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_justHighMET, CombRelIso, this_weight);
			}
			// AES: just tighter Z
			if (flag_AES_pass_tighterZveto_mee && flag_AES_pass_tighterZveto_mep)
				fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_justZ, CombRelIso, this_weight);

		}//pass some cuts
		// TL: end ------------------------------


		//=====================================================================
		//
		//
		//    QCD Control Region (AES)        plan A      4-12-09
		//
		//--------------------------------------------------------
		// Trial AES definition for Plan A: MET>20.
		//
		// Common: met<15, HT<150,  exactly 1 GoodEle (RT).
		//
		// Definitions: A1 : leave out !conv cut
		//              A2 : invert !conv cut, ie plot only conversion ele.
		//              A3 : invert d0 cut = d0 > 200um.
		//--------------------------------------------------------
		// First make sure selected events do not enter control sample
		//     if(m_debug) cout << << endl;
		if (e_plus_jet_pass == false) {

			if (goodrun && fired_single_em && nGoodEle == 1 && !isMuon && !isZ && !isDifferentInteraction) {

				// the Good Electron
				unsigned int theGE = ii_electrons.at(0); //first good ele

				float ht_AES = this_met + els_et->at(theGE);
				for (unsigned int i = 0; i < goodJets.size(); ++i) {
					ht_AES += goodJets.at(i).Pt();
				}

				if (this_met < AES_MET_cut && ht_AES < AES_HT_cut) {
					// satisfy criteria for plan A control sample

					const double this_et = els_et->at(theGE);
					const double this_iso = getRelIso(theGE);

					// A1: no conv cuts
					//cout << "[M_DEBUG] A1" << endl;
					if (this_et > 20.) {
						fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planA1_e20, this_iso, this_weight);
						if (this_et > 30.)
							fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planA1_e30, this_iso, this_weight);
					}

					// A2: conv only
					//cout << "[M_DEBUG] A2" << endl;
					//	   TLorentzVector eles_temp(els_px->at(0),els_py->at(0),els_pz->at(0),els_energy->at(0));
					bool this_is_conv = ConversionFinder(electrons.at(0), mctype, 0); //1st good ele
					if (this_is_conv) {
						if (this_et > 20.) {
							fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planA2_e20, this_iso, this_weight);
							if (this_et > 30.)
								fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planA2_e30, this_iso, this_weight);
						}
					}

					// A3: invert d0 (>200 um)
					//cout << "[M_DEBUG] A3" << endl;
					if (fabs(compute_d0("electron", theGE)) > 0.02) {
						if (this_et > 20.) {
							fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planA3_e20, this_iso, this_weight);
							if (this_et > 30.)
								fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planA3_e30, this_iso, this_weight);
						}
					}
				}
			}
		}
		// END: plan A control sample
		//=====================================================================


		//=====================================================================
		//
		//    QCD Control Region (AES)        plan B      10-11-09
		//
		//----------------------------------------------------
		// Trial AES definition for Plan B: noMET, eta(e)<1.5
		//
		// Definitions: B1: - 1 electron with ET>20, |eta|<1.442, d0<200um, RT, not conversion.
		//                  - 1 jet
		//                  - simple Z veto: no 2 reco ele
		//                  - std cuts: HLT, !mu
		//
		// Definitions: B2: - 1 electron with ET<30, |eta|<1.442, d0<200um, RT, not conversion.
		//                  - 1 jet
		//                  - simple Z veto: no 2 reco ele
		//                  - std cuts: HLT, !mu
		// Definitions: B3: - 1 electron withh ET>20/30, |eta|<1.442, d0<200um, fail RT, not conversion.
		//                  - 1 jet
		//                  - simple Z veto: no 2 reco ele
		//                  - std cuts: HTL, !mu
		// Definitions: B4: same as B3 but failing RL ID.
		// Definitions: B5: same as B3 but failing category-based Loose ID.
		// Definitions: B6: same as B3 but failing category-based Tight ID.
		// Definitions: B7: same as B1 but inverting d0 cut, req d0>200um..
		// Definitions: B8: same as B7 but fail RT.
		//-----------------------------------------------------
		// First make sure selected events do not enter control sample
		if (e_plus_jet_pass == false) {
			if (goodrun && fired_single_em && Nels == 1 && !isMuon && !isDifferentInteraction) {//AAAAA

				// common criteria for B1 & 2
				//if(  fabs(els_eta->at(0))           < 1.442  &&
				if (fabs(els_eta->at(0)) < 2.5 && fabs(compute_d0("electron", 0)) < 0.02) {

					//FIXME: if true? obsolete if
					if (true) {//FIXME: no check for the variable which creates histograms

						const double this_iso = getRelIso(0);

						if (els_robustTightId->at(0) > 0) { //pass RT ID

							// Definition B1:  lowered EleET cut
							//------------------------------------
							if (els_et->at(0) > 20.) {
								fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planB1_e20, this_iso, this_weight);
								if (els_et->at(0) > 30.)
									fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planB1_e30, this_iso, this_weight);
							}
							// Definition B2:  inverted EleET cut
							//------------------------------------
							if (els_et->at(0) < 30.) {
								fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planB2_e30, this_iso, this_weight);
								if (els_et->at(0) < 20.)
									fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planB2_e20, this_iso, this_weight);
							}

						} else {
							// Definition B3: fail RT ID
							//---------------------------
							if (els_et->at(0) > 20.) {
								fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planB3_e20, this_iso, this_weight);
								if (els_et->at(0) > 30.)
									fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planB3_e30, this_iso, this_weight);
							}
							if (fabs(els_eta->at(0)) < 1.442) { //BARREL
								if (els_et->at(0) > 20.) {
									fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planB3b_e20, this_iso, this_weight);
									if (els_et->at(0) > 30.)
										fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planB3b_e30, this_iso, this_weight);
								}
							}
						}
						// 2-12-09
						if (els_et->at(0) > 20.) {
							if (els_robustLooseId->at(0) == 0)
								fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planB4_e20, this_iso, this_weight);
							if (els_looseId->at(0) == 0)
								fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planB5_e20, this_iso, this_weight);
							if (els_tightId->at(0) == 0)
								fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planB6_e20, this_iso, this_weight);

							if (els_et->at(0) > 30.) {
								if (els_robustLooseId->at(0) == 0)
									fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planB4_e30, this_iso, this_weight);
								if (els_looseId->at(0) == 0)
									fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planB5_e30, this_iso, this_weight);
								if (els_tightId->at(0) == 0)
									fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planB6_e30, this_iso, this_weight);
							}
						}

					}
				}
				// B7: invert d0 cut
				if (els_et->at(0) > 20. && fabs(els_eta->at(0)) < 2.5 && fabs(compute_d0("electron", 0)) > 0.02) { //<---

					const double this_iso = getRelIso(0);
					if (els_robustTightId->at(0) > 0) { //pass RT ID
						fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planB7_e20, this_iso, this_weight);
						if (els_et->at(0) > 30.)
							fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planB7_e30, this_iso, this_weight);
					} else { //fail RT ID
						fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planB8_e20, this_iso, this_weight);
						if (els_et->at(0) > 30.)
							fillHisto_Njet_DataAndMC(h_QCDest_CombRelIso_AES_planB8_e30, this_iso, this_weight);
					}//RT ID
				}
			}
		}// plan B control sample
		//=====================================================================


		//---------------------------------------------------------------
		// 9 Jun 09: inspect isolation/met at each key stage of cuts
		//           to study correlation of iso:met for QCD events
		// level 1: HLT
		//       2: nGoodEle>0
		//       3: mu, Z veto
		//       4: conversion cut
		//---------------------------------------------------------------
		if (m_plotRelisoNES) {

			if (goodrun && fired_single_em) {

				string metside = "loMET"; //low met
				ushort mettype = 0;//lomet
				if (this_met > METCUT) {
					metside = "hiMET"; //high met
					mettype = 3;
				}

				// Pass L1
				if (m_debug)
					cout << "-> Filling Reliso NES histograms, L1" << endl;
				// fill for all electrons
				for (unsigned int ie = 0; ie < Nels; ie++) {

					float tmpRelIso = getRelIso(ie);

					// barrel or endcap?
					string etaside = "barrel";
					ushort etatype = 0; //barrel
					if (fabs(els_eta->at(ie)) > 1.5) {
						etaside = "endcap";
						etatype = 1;
					}
					iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES, Level_1, tmpRelIso, this_met, this_weight);
					iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES_uw, Level_1, tmpRelIso, this_met, 1);

					iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES, Level_1, tmpRelIso, this_weight);//no met cut
					iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_barrel + etatype, Level_1, tmpRelIso, this_weight);//no met cut
					iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET + mettype, Level_1, tmpRelIso, this_weight);
					iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET_barrel + etatype + mettype, Level_1, tmpRelIso,
							this_weight);

					// fill only for electrons with at least 30 GeV ET, and pass eta cut
					if (els_et->at(ie) > ELE_ETCUT && fabs(els_eta->at(ie)) < 2.5 && (fabs(els_eta->at(ie)) < 1.442 || fabs(
							els_eta->at(ie)) > 1.56)) {
						iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES, Level_1b, tmpRelIso, this_met, this_weight);
						iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES_uw, Level_1b, tmpRelIso, this_met, 1);

						iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES, Level_1b, tmpRelIso, this_weight);//no met cut
						iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_barrel + etatype, Level_1b, tmpRelIso, this_weight);//no met cut
						iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET + mettype, Level_1b, tmpRelIso, this_weight);
						iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET_barrel + etatype + mettype, Level_1b,
								tmpRelIso, this_weight);

						// d0 cut
						// Calculate d0 w.r.t beam spot
						float d0_corrected = compute_d0("electron", ie);

						if (fabs(d0_corrected) < 0.02) {
							iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES, Level_1c, tmpRelIso, this_met, this_weight);
							iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES_uw, Level_1c, tmpRelIso, this_met, 1);

							iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES, Level_1c, tmpRelIso, this_weight);//no met cut
							iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_barrel + etatype, Level_1c, tmpRelIso,
									this_weight);//no met cut
							iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET + mettype, Level_1c, tmpRelIso, this_weight);
							iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET_barrel + etatype + mettype, Level_1c,
									tmpRelIso, this_weight);

							// eID (barrel)
							bool pass_eid_c0 = passEleID(ie); //out-of-box eID variable
							bool pass_eid_c1 = false;//Question: What are they?
							bool pass_eid_c2 = false;
							bool pass_eid_c3 = false;
							bool pass_eid_c4 = false;
							bool pass_eid = false;

							// electron ID variables
							//
							// 31X: eID cuts have changed!!! values need to be updated!  TODO
							//

							if (fabs(els_eta->at(ie)) < 1.442) { // barrel
								// Updated eID cut values: tight Fixed Threshold "RobustTight"
								// ref: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCutBasedElectronID
								if (els_hadOverEm->at(ie) < 0.01)
									pass_eid_c1 = true;
								if (fabs(els_dEtaIn->at(ie)) < 0.0040)
									pass_eid_c2 = true;
								if (fabs(els_dPhiIn->at(ie)) < 0.025)
									pass_eid_c3 = true;
								if (els_sigmaIEtaIEta->at(ie) < 0.0099)
									pass_eid_c4 = true;

							} else if (fabs(els_eta->at(ie)) > 1.56) { // endcap

								if (els_hadOverEm->at(ie) < 0.01)
									pass_eid_c1 = true;
								if (fabs(els_dEtaIn->at(ie)) < 0.0066)
									pass_eid_c2 = true;
								if (fabs(els_dPhiIn->at(ie)) < 0.020)
									pass_eid_c3 = true;
								if (els_sigmaIEtaIEta->at(ie) < 0.028)
									pass_eid_c4 = true;

							}

							pass_eid = pass_eid_c0 && pass_eid_c1 && pass_eid_c2 && pass_eid_c3 && pass_eid_c4; //all 5

							if (pass_eid_c1) {
								iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES, Level_1d1, tmpRelIso, this_met, this_weight);
								iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES_uw, Level_1d1, tmpRelIso, this_met, 1);

								iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES, Level_1d1, tmpRelIso, this_weight);//no met cut
								iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_barrel + etatype, Level_1d1, tmpRelIso,
										this_weight);//no met cut
								iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET + mettype, Level_1d1, tmpRelIso,
										this_weight);
								iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET_barrel + etatype + mettype, Level_1d1,
										tmpRelIso, this_weight);

								if (pass_eid_c2) {
									iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES, Level_1d2, tmpRelIso, this_met, this_weight);
									iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES_uw, Level_1d2, tmpRelIso, this_met, 1);

									iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES, Level_1d2, tmpRelIso, this_weight);//no met cut
									iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_barrel + etatype, Level_1d2, tmpRelIso,
											this_weight);//no met cut
									iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET + mettype, Level_1d2, tmpRelIso,
											this_weight);
									iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET_barrel + etatype + mettype,
											Level_1d2, tmpRelIso, this_weight);

									if (pass_eid_c3) {
										iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES, Level_1d3, tmpRelIso, this_met,
												this_weight);
										iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES_uw, Level_1d3, tmpRelIso, this_met, 1);

										iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES, Level_1d3, tmpRelIso, this_weight);//no met cut
										iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_barrel + etatype, Level_1d3,
												tmpRelIso, this_weight);//no met cut
										iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET + mettype, Level_1d3,
												tmpRelIso, this_weight);
										iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET_barrel + etatype + mettype,
												Level_1d3, tmpRelIso, this_weight);

										if (pass_eid_c4) {
											iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES, Level_1d4, tmpRelIso, this_met,
													this_weight);
											iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES_uw, Level_1d4, tmpRelIso, this_met,
													1);

											iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES, Level_1d4, tmpRelIso,
													this_weight);//no met cut
											iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_barrel + etatype, Level_1d4,
													tmpRelIso, this_weight);//no met cut
											iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET + mettype, Level_1d4,
													tmpRelIso, this_weight);
											iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET_barrel + etatype + mettype,
													Level_1d4, tmpRelIso, this_weight);

											if (pass_eid_c0) {
												iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES, Level_1d5, tmpRelIso, this_met,
														this_weight);
												iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES_uw, Level_1d5, tmpRelIso,
														this_met, 1);

												iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES, Level_1d5, tmpRelIso,
														this_weight);//no met cut
												iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_barrel + etatype, Level_1d5,
														tmpRelIso, this_weight);//no met cut
												iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET + mettype, Level_1d5,
														tmpRelIso, this_weight);
												iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET_barrel + etatype
														+ mettype, Level_1d5, tmpRelIso, this_weight);

											}//pass c0
										}//pass c4
									}//pass c3
								}//pass c2
							}//pass c1
						}//pass d0
					}//pass et,eta
				}// loop over Nels


				if (nGoodEle > 0) {

					// barrel or endcap?
					if (ii_GoodEle_mostIso < 0)
						cout << "error/warning: could not find most isolated GoodEle." << endl;
					string etaside = "barrel";
					ushort etatype = 0;
					if (fabs(els_eta->at(ii_GoodEle_mostIso)) > 1.56) {
						etaside = "endcap";
						etatype = 1;
					}

					if (CombRelIso < 0) {
						static short mm = 0;
						if (mm == 0)
							cout << "(Printing the first 20 occurances.)" << endl;
						if (mm < 20) {
							cout << "** attention: negative CombRelIso (" << CombRelIso << ")" << endl;
							++mm;
						}
					}

					if (m_debug)
						cout << "-> Filling Reliso NES histograms, L2" << endl;

					iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES, Level_2, CombRelIso, this_met, this_weight);
					iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES_uw, Level_2, CombRelIso, this_met, 1);

					iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES, Level_2, CombRelIso, this_weight);//no met cut
					iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_barrel + etatype, Level_2, CombRelIso, this_weight);//no met cut
					iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET + mettype, Level_2, CombRelIso, this_weight);
					iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET_barrel + etatype + mettype, Level_2, CombRelIso,
							this_weight);

					if (!isMuon && !isZ) {

						if (m_debug)
							cout << "-> Filling Reliso NES histograms, L3" << endl;
						iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES, Level_3, CombRelIso, this_met, this_weight);
						iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES_uw, Level_3, CombRelIso, this_met, 1);

						iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES, Level_3, CombRelIso, this_weight);//no met cut
						iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_barrel + etatype, Level_3, CombRelIso, this_weight);//no met cut
						iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET + mettype, Level_3, CombRelIso, this_weight);
						iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET_barrel + etatype + mettype, Level_3,
								CombRelIso, this_weight);

						if (!isConversion && !isDifferentInteraction) { //NB: no HT cut

							if (m_debug)
								cout << "-> Filling Reliso NES histograms, L4" << endl;
							iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES, Level_4, CombRelIso, this_met, this_weight);
							iso_fillHisto_nlevel_nj_nmc(h_QCDest_isoVmet_NES_uw, Level_4, CombRelIso, this_met, 1);

							iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES, Level_4, CombRelIso, this_weight);//no met cut
							iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_barrel + etatype, Level_4, CombRelIso,
									this_weight);//no met cut
							iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET + mettype, Level_4, CombRelIso, this_weight);
							iso_fillHisto_nlevel_nj_nmc(h_QCDest_CombRelIso_NES_loMET_barrel + etatype + mettype, Level_4,
									CombRelIso, this_weight);
						}
					}
				}//nGoodEle > 0
			}//goodrun, trigger
		}//m_plotRelisoNES switch


		// study Z veto ---------
		if (goodrun && fired_single_em && nGoodEle > 0 && !isMuon && !isZ && !isConversion && !isDifferentInteraction && ht
				>= HTCUT) {

			if (isZee) {
				nEvent_Zee_pass++;
				if (this_met > METCUT)
					nEvent_Zee_highMET++;
				else
					nEvent_Zee_lowMET++;
			} else if (isZmm) {
				nEvent_Zmm_pass++;
				if (this_met > METCUT)
					nEvent_Zmm_highMET++;
				else
					nEvent_Zmm_lowMET++;
			} else if (isZtt) {
				nEvent_Ztt_pass++;
				if (this_met > METCUT)
					nEvent_Ztt_highMET++;
				else
					nEvent_Ztt_lowMET++;
			}
		}

		// Add m3 plots for W+jets estimation
		//---------------------------------------
		// Plot for event with 4 or more jets only
		// Require at least 3 jest to reconstruct t->jjj

		// notes: all cuts except ISO
		if (goodrun && fired_single_em && nGoodEle > 0 && this_met > METCUT && !isMuon && !isZ && !isConversion
				&& !isDifferentInteraction && ht >= HTCUT && nGoodJet >= 4) {
			reco_hadronicTop_highestTopPT(goodJets, nGoodIsoEle);
			//			reco_mttbar(jets, electrons, met, nGoodIsoEle);
			//new version with Chi2 fit
			reco_Mttbar(goodJets, electrons[0], met[0]);
			fillBtagHistograms(goodJets, btag_information);
		}

		// Add Delta R(e,mu)
		//-------------------
		// Get isolated electrons
		for (size_t e = 0; e < iso_electrons.size(); e++) {

			if (iso_electrons.at(e).Pt() < 1)
				cout << "iso ele pt < 1, ev:" << ev << endl;
			// Get all GoodMuon
			for (size_t m = 0; m < muons.size(); m++) {
				if (muons.at(m).Pt() < 1)
					cout << "mu pt < 1, ev:" << ev << endl;
				float DRemu = (iso_electrons.at(e)).DeltaR(muons.at(m));
				if (ev < 100)
					cout << "ev: " << ev << ", DRemu: " << DRemu << endl;
				//				fasthist_[fastmctype_][h_DRemu_selE_GoodMu]->Fill(DRemu, this_weight);
				fillHistoDataAndMC(h_DRemu_selE_GoodMu, DRemu, this_weight);
				if (e_plus_jet_pass && Njet() >= 4) { // pass ALL cuts
					//					fasthist_[fastmctype_][h_DRemu_selE_GoodMu_pass]->Fill(DRemu, this_weight);
					fillHistoDataAndMC(h_DRemu_selE_GoodMu_pass, DRemu, this_weight);
					if (DRemu < 0.1)
						nEvent_DR_ele_muo_less_than_01++;
				}
			}//mu
		}//e


		if (m_debug) {
			cout << "\n---------------------------";
			cout << "\n End: Clearing vectors";
			cout << "\n---------------------------\n";
		}
		electrons.clear();
		ii_electrons.clear();
		electrons_barrel.clear();
		electrons_endcap.clear();
		iso_electrons.clear();
		iso_electrons_barrel.clear();
		iso_electrons_endcap.clear();
		muons.clear();
		iso_muons.clear();
		goodJets.clear();
		for (ushort btag = 0; btag < NUMBER_OF_BTAGS; ++btag) {
			btag_information[btag].clear();
		}

		met.clear();
		electrons_isoval.clear();
		electrons_isoval2.clear();

		//---------------------------
		// PDF uncertainties (begin)
		//---------------------------
		if (!IsData() && isMCType(kttbar) && m_studyPDFunc) {
			if (ev == 1) {
				cout << "PDF uncertainties" << endl;
				cout << "PDFWcteq66_0   " << PDFWcteq66_0 << endl;
				cout << "PDFWcteq66_1   " << PDFWcteq66_1 << endl;
				cout << "PDFWcteq66_2   " << PDFWcteq66_2 << endl;
				cout << "PDFWcteq66_3   " << PDFWcteq66_3 << endl;
				cout << "PDFWcteq66_4   " << PDFWcteq66_4 << endl;
				cout << "PDFWcteq66_10  " << PDFWcteq66_10 << endl;
				cout << "PDFWcteq66_20  " << PDFWcteq66_20 << endl;
				cout << "PDFWcteq66_30  " << PDFWcteq66_30 << endl;
				cout << "PDFWcteq66_40  " << PDFWcteq66_40 << endl;
				cout << "PDFWcteq66_44  " << PDFWcteq66_44 << endl;
			}
			fillHisto_PDF_weights(h_pdf_total);
			if (e_plus_jet_pass)
				fillHisto_PDF_weights(h_pdf_pass);
		}
		//-------------------------
		// PDF uncertainties (end)
		//-------------------------


	}// loop over all the events

	//----------------------------------------------------------------------------
	//                         End  of  Event  Loop
	//----------------------------------------------------------------------------


	if (!IsData()) {
		cout << "\n (num of event flagged as Z, di-e inv mass)" << endl;
		cout << " nZ_ORIG (old) = " << nZ_ORIG << endl;
		cout << " nZ_TEST (new) = " << nZ_TEST << endl << endl;
	}

	const bool study_z_veto = false;
	if (study_z_veto) {
		cout << endl;
		cout << "nEvent looped over (Z+jets)            = " << nEvent_Z << endl;
		cout << "|->nEvent which is really Z->ee        = " << nEvent_Zee << endl;
		cout << "|->nEvent which is really Z->mumu      = " << nEvent_Zmm << endl;
		cout << "|->nEvent which is really Z->tautau    = " << nEvent_Ztt << endl;
		cout << endl;
		cout << "|->nEvent which is really Z->ee        = " << nEvent_Zee << endl;
		cout << "  |-> 2 or more reco-ele (no cut)      = " << nEvent_2orMoreEle << endl;
		cout << "    |-> subset: match to genEle Z->ee  = " << nEvent_EleMatch << " (" << (float) nEvent_EleMatch
				/ nEvent_2orMoreEle * 100.0 << " %)" << endl;
		cout << "       |-> pass cut (2 or more):" << endl;
		cout << "       |-> subset: et only             = " << nEvent_Z_et << endl;
		cout << "       |-> subset: eta only            = " << nEvent_Z_eta << endl;
		cout << "       |-> subset: et,eta              = " << nEvent_Z_eteta << endl;
		cout << "       |-> subset: et,eta,d0           = " << nEvent_Z_eteta_d0 << endl;
		cout << "       |-> subset: et,eta,d0,idL       = " << nEvent_Z_eteta_d0_L << endl;
		cout << "       |-> subset: et,eta,d0,idT       = " << nEvent_Z_eteta_d0_T << endl;
		cout << "       |-> subset: et,eta,d0,idRL      = " << nEvent_Z_eteta_d0_RL << endl;
		cout << "       |-> subset: et,eta,d0,idRT      = " << nEvent_Z_eteta_d0_RT << endl;
		cout << "|-> true Z, reco-ele, pass et,eta:" << endl;
		cout << "    nEvent     |-> 0 reco-e            = " << fasthist_jetbinned_[kZjets][AllJets][h_nEle_s1]->GetBinContent(1)
				<< endl; //allj
		cout << "               |-> 1 reco-e            = " << fasthist_jetbinned_[kZjets][AllJets][h_nEle_s1]->GetBinContent(2)
				<< endl;
		cout << "               |-> 2 reco-e            = " << fasthist_jetbinned_[kZjets][AllJets][h_nEle_s1]->GetBinContent(3)
				<< endl;
		cout << "               |-> 3 reco-e            = " << fasthist_jetbinned_[kZjets][AllJets][h_nEle_s1]->GetBinContent(4)
				<< endl;
		cout << "               |-> 4 reco-e            = " << fasthist_jetbinned_[kZjets][AllJets][h_nEle_s1]->GetBinContent(5)
				<< endl;
		cout << "|-> true Z->ee, GenEle pass et,eta:" << endl;
		cout << "    nEvent     |-> 0 gen-e             = " << fasthist_[fastmctype_][h_nGenBasicEle_Zee_allj]->GetBinContent(1)
				<< endl;
		cout << "               |-> 1 gen-e             = " << fasthist_[fastmctype_][h_nGenBasicEle_Zee_allj]->GetBinContent(2)
				<< endl;
		cout << "               |-> 2 gen-e             = " << fasthist_[fastmctype_][h_nGenBasicEle_Zee_allj]->GetBinContent(3)
				<< endl;
		cout << "               |-> 3 gen-e             = " << fasthist_[fastmctype_][h_nGenBasicEle_Zee_allj]->GetBinContent(4)
				<< endl;
		cout << "               |-> 4 gen-e             = " << fasthist_[fastmctype_][h_nGenBasicEle_Zee_allj]->GetBinContent(5)
				<< endl;
		cout << "Zee --> pass QCDest scope (reliso)     = " << nEvent_Zee_pass << endl;
		cout << "    |-> high met                       = " << nEvent_Zee_highMET << endl;
		cout << "    |-> low met                        = " << nEvent_Zee_lowMET << endl;
		cout << "Zmm --> pass QCDest scope (reliso)     = " << nEvent_Zmm_pass << endl;
		cout << "    |-> high met                       = " << nEvent_Zmm_highMET << endl;
		cout << "    |-> low met                        = " << nEvent_Zmm_lowMET << endl;
		cout << "Ztt --> pass QCDest scope (reliso)     = " << nEvent_Ztt_pass << endl;
		cout << "    |-> high met                       = " << nEvent_Ztt_highMET << endl;
		cout << "    |-> low met                        = " << nEvent_Ztt_lowMET << endl;
		cout << endl;
	}
	//---------
	fprintf(outfile, "\nList of files that contain selected events:\n\n");
	for (set<string>::const_iterator it = interestingFiles.begin(); it != interestingFiles.end(); ++it) {
		fprintf(outfile, "%s\n", it->c_str());
	}
	fclose(outfile);

	//MC - make bin zero the total bin for data-like print-out of numbers of events
	for (int i = 0; i < nstage; ++i) {
		for (int j = 0; j < ntjet; ++j) {
			for (int k = 1; k < nmctype; ++k) {
				e_plus_jet[i][j][0] += e_plus_jet[i][j][k];
				e_plus_jet_weighted[i][j][0] += e_plus_jet_weighted[i][j][k];
			}
		}
	}

	// Fill histograms (event count tables) (for MC only)
	if (!IsData()) {
		for (short i = 0; i < nstage; ++i) {
			for (short j = 0; j < ntjet; ++j) {
				if (hasSampleMC(kttbar)) {//MC contains Signal
					for (short k = 1; k < 11; ++k) {
						Signal_njetsVcuts->Fill(j, nstage - i - 1, e_plus_jet_weighted[i][j][k]);
						Signal_njetsVcuts->Fill(5, nstage - i - 1, e_plus_jet_weighted[i][j][k]);//sum for all jets
					}
				}
				if (hasSampleMC(kQCD)) {//MC contains QCD
					for (short k = 13; k < 19; ++k) { // mctype is 13 to 18 for QCD
						QCD_njetsVcuts->Fill(j, nstage - i - 1, e_plus_jet_weighted[i][j][k]);
						QCD_njetsVcuts->Fill(5, nstage - i - 1, e_plus_jet_weighted[i][j][k]);//sum for all jets
					}
				}
				if (hasSampleMC(kWjets)) {//MC contains Wjets (mctype=11)
					Wjets_njetsVcuts->Fill(j, nstage - i - 1, e_plus_jet_weighted[i][j][11]);
					Wjets_njetsVcuts->Fill(5, nstage - i - 1, e_plus_jet_weighted[i][j][11]);//sum for all jets
				}
				if (hasSampleMC(kZjets)) {//MC contains Zjets (mctype=12)
					Zjets_njetsVcuts->Fill(j, nstage - i - 1, e_plus_jet_weighted[i][j][12]);
					Zjets_njetsVcuts->Fill(5, nstage - i - 1, e_plus_jet_weighted[i][j][12]);//sum for all jets
				}
				if (hasSampleMC(kVQQ)) {//MC contains VQQ (mctype=19)
					VQQ_njetsVcuts->Fill(j, nstage - i - 1, e_plus_jet_weighted[i][j][19]);
					VQQ_njetsVcuts->Fill(5, nstage - i - 1, e_plus_jet_weighted[i][j][19]);//sum for all jets
				}
				if (hasSampleMC(ksingleTop)) {//MC contains single top (mctypes are 20-22)
					for (int k = 20; k < 23; ++k) {
						SingleTop_njetsVcuts->Fill(j, nstage - i - 1, e_plus_jet_weighted[i][j][k]);
						SingleTop_njetsVcuts->Fill(5, nstage - i - 1, e_plus_jet_weighted[i][j][k]);//sum for all jets
					}
				}
			}
		}
	}//end MC

	// Set labels of the histograms (cuts and Njet)
	SetHistoLabelCutNjet(Signal_njetsVcuts, ve);
	SetHistoLabelCutNjet(QCD_njetsVcuts, ve);
	SetHistoLabelCutNjet(Wjets_njetsVcuts, ve);
	SetHistoLabelCutNjet(Zjets_njetsVcuts, ve);
	SetHistoLabelCutNjet(VQQ_njetsVcuts, ve);
	SetHistoLabelCutNjet(SingleTop_njetsVcuts, ve);

	// Fill Signal accpetance and efficiency histograms
	if (hasSampleMC(kttbar)) {
		float total = 0;
		float totalb = 0;
		for (int k = 1; k < 11; ++k) {
			totalb += e_plus_jet[9][4][k];
			for (int j = 0; j < ntjet; ++j) {
				total += e_plus_jet[0][j][k];
			}
		}

		const double SignalAcc = totalb / total;
		const double SignalAccUnc = sqrt(SignalAcc * (1 - SignalAcc) / totalb);

		SignalVars->SetBinContent(1, SignalAcc);
		SignalVars->SetBinContent(2, SignalAccUnc);

		SignalVars->SetBinContent(3, totalb); //signal passed event selection
		SignalVars->SetBinContent(4, e_plus_jet[9][4][1]);//evqq channel passed event selection
	}

	// Histogram for PDF uncertainties
	if (m_studyPDFunc) {
		h_pdf_eff->Divide(h_pdf_pass, h_pdf_total, 1., 1., "B");
	}

	//Summary Info
	cout << "\n\\begin{verbatim}";
	cout << "\n%------------------------------------------------------------------------------------\n";
	cout << "                             Summary Information";
	cout << "\n%------------------------------------------------------------------------------------\n";

	cout << "Integrated luminosity assumed = " << intlumi << "/pb" << endl;
	if (GetTrigger())
		cout << "Applied HLT trigger: " << HLTBit << endl;
	cout << " Electron ET cut =  " << ELE_ETCUT << "  GeV" << endl;
	cout << "     Muon PT cut =  " << MU_PTCUT << "  GeV" << endl;
	cout << "      Jet PT cut =  " << JET_PTCUT << "  GeV" << endl;
	cout << "      MET    cut =  " << METCUT << "  GeV" << endl;
	cout << "       HT    cut =  " << HTCUT << "  GeV" << endl;
	cout << " Electron ID = " << printEleID() << endl;
	cout << " Electron RelIso formula = ";
	if (useNewReliso)
		cout << "new";
	else
		cout << "old";
	if (!m_applyMETcut && m_rejectEndcapEle)
		cout << "\n Using plan B, does not use MET, use barrel ele only (eta<1.442)" << endl;
	if (m_applyMETcut && m_rejectEndcapEle)
		cout << "\n Using plan C, both MET & ele eta cuts are used (eta<1.442)" << endl;
	cout << "\n\\end{verbatim}\n" << endl;

	cout << "\n Monitor:" << endl;
	cout << " Number of event with DR(e,mu) less than 0.1  =  " << nEvent_DR_ele_muo_less_than_01 << endl;
	cout << " e: isolated electron" << endl;
	cout << " mu: all GoodMuon\n" << endl;

	cout << endl << "%------------------------------------------------------------------";
	cout << endl << "%------------------------------------------------------------------";
	cout << endl << "           Event counts per jet multiplicity";
	cout << endl << "%------------------------------------------------------------------";
	cout << endl << "%------------------------------------------------------------------";
	cout << endl << "\\\\\n";
	cout << endl << "\\begin{tabular}{|l|rrrrr|r|}\n\\hline";

	cout << endl << "Unweighted numbers" << endl;
	cout.precision(0);
	DrawEventPerNjetTable(e_plus_jet, ve); //unweighted number

	if (IsData() == false) {
		cout << "Normalized for " << intlumi << " pb$^{-1}$" << endl;
		cout.precision(myprec);
		//cout.fixed;
		DrawEventPerNjetTable(e_plus_jet_weighted, ve); //weighted number
	}
	cout << "\\end{tabular}\n" << endl;

	// FB+TL (21-1-09)
	// When running on MC, produce tables of event count according to type of MC

	if (IsData() == false) {

		cout << endl;
		DrawSignalBGTable(e_plus_jet_weighted, ve); //weighted table

		DrawMCTypeTable(e_plus_jet, "Events Table (per MC type; unweighted)", ve); //unweighted table
		DrawMCTypeTable(e_plus_jet_weighted, "Events Table (per MC type; weighted)  ", ve); //weighted table

		// Break down table for QCD
		if (hasSampleMC(kQCD)) {
			DrawQCDTable(e_plus_jet, "QCD Table (unweighted)", ve);
			DrawQCDTable(e_plus_jet_weighted, "QCD Table (weighted)  ", ve);
		}

		// Break down table for single top
		if (hasSampleMC(ksingleTop)) {
			DrawSingleTopTable(e_plus_jet, "Single Top Table (unweighted)", ve);
			DrawSingleTopTable(e_plus_jet_weighted, "Single Top Table (weighted)  ", ve);
		}

		// Acceptance note printout (ttbar)
		// >= 3 tight jets, >= 4 tight jets

		//cout.precision(0); //reset precision

		// TL (21-1-09):  only run this if we've ran over some ttbar MC in input
		if (hasSampleMC(kttbar)) {

			DrawSignalAcceptanceTable(e_plus_jet, ve);

			// Set labels for all_mctype histogram (ttbar decay modes)
			all_mctype->GetXaxis()->SetBinLabel(1, "evqq"); //semilep
			all_mctype->GetXaxis()->SetBinLabel(2, "mvqq");
			all_mctype->GetXaxis()->SetBinLabel(3, "tvqq");
			all_mctype->GetXaxis()->SetBinLabel(4, "evev"); //dilepton
			all_mctype->GetXaxis()->SetBinLabel(5, "mvmv");
			all_mctype->GetXaxis()->SetBinLabel(6, "tvtv");
			all_mctype->GetXaxis()->SetBinLabel(7, "evmu");
			all_mctype->GetXaxis()->SetBinLabel(8, "evtv");
			all_mctype->GetXaxis()->SetBinLabel(9, "mvtv");
			all_mctype->GetXaxis()->SetBinLabel(10, "qqqq"); //fully hadronic

		}//end mc_sample_has_ttbar


		// (19 Feb 09) print how many QCD events we have after all cuts except reliso nad njet
		if (hasSampleMC(kQCD)) {
			const int QCD_bc = kQCD;
			cout << "\n New RelIso  mc" << setw(10) << intlumi << "/pb" << endl;
			cout << "   1j" << setw(10) << fasthist_jetbinned_[QCD_bc][OneJet][h_QCDest_CombRelIso]->GetEntries() << setw(10)
					<< fasthist_jetbinned_[QCD_bc][OneJet][h_QCDest_CombRelIso]->Integral() << endl;
			cout << "   2j" << setw(10) << fasthist_jetbinned_[QCD_bc][TwoJets][h_QCDest_CombRelIso]->GetEntries() << setw(10)
					<< fasthist_jetbinned_[QCD_bc][TwoJets][h_QCDest_CombRelIso]->Integral() << endl;
			cout << "   3j" << setw(10) << fasthist_jetbinned_[QCD_bc][ThreeJets][h_QCDest_CombRelIso]->GetEntries() << setw(10)
					<< fasthist_jetbinned_[QCD_bc][ThreeJets][h_QCDest_CombRelIso]->Integral() << endl;
			cout << " >=4j" << setw(10) << fasthist_jetbinned_[QCD_bc][FourOrMoreJets][h_QCDest_CombRelIso]->GetEntries()
					<< setw(10) << fasthist_jetbinned_[QCD_bc][FourOrMoreJets][h_QCDest_CombRelIso]->Integral() << endl;
			//			cout << "   2j" << setw(10) << h_QCDest_CombRelIso[2][QCD_bc]->GetEntries() << setw(10)
			//					<< h_QCDest_CombRelIso[2][QCD_bc]->Integral() << endl;
			//			cout << "   3j" << setw(10) << h_QCDest_CombRelIso[3][QCD_bc]->GetEntries() << setw(10)
			//					<< h_QCDest_CombRelIso[3][QCD_bc]->Integral() << endl;
			//			cout << " >=4j" << setw(10) << h_QCDest_CombRelIso[5][QCD_bc]->GetEntries() << setw(10)
			//					<< h_QCDest_CombRelIso[5][QCD_bc]->Integral() << endl;
		}

	}//end MC


	cout.precision(myprec); //reset precision

	//Print event tables with errors
	PrintErrorTables(e_plus_jet, e_plus_jet_weighted, ve); //TEMPORARY Turn off

	//Print event tables with errors, but without applying scientific notation
	ScientificNotation = false;
	PrintErrorTables(e_plus_jet, e_plus_jet_weighted, ve); //TEMPORARY Turn off
	End();

	// Close the histogram file
	histf->Write();
	histf->Close();

	//856
	if (DoConversionStudies()) {
		PrintConversionTable();
	}
	//end 856

	cout << "\n***********************************************************************";
	cout << "\n*                                                                     *";
	cout << "\n*           A N A L Y S I S       C O M P L E T E D                   *";
	cout << "\n*                                                                     *";
	cout << "\n***********************************************************************";
	cout << "\n*   Available events:  " << left << setw(47) << nEventsAvail << "*";
	cout << "\n*   Limit:             " << left << setw(47) << nEvents << "*";
	cout << "\n*   Processed events:  " << left << setw(47) << counter << "*";
	cout << "\n*   Passed events:     " << left << setw(47) << counter_pass << "*";
	cout << "\n***********************************************************************\n";

	// (TL) Print current time
	TDatime now;
	cout << endl << "Current local ";
	now.Print();
	cout << endl;

	return true;

}// end EventLoop()


void ana::PrintErrorTables(const double e_plus_jet[][5][nmctype], const double e_plus_jet_weighted[][5][nmctype],
		vector<string> ve) const {

	const int mynstage = 11; //up to !DIFFZ

	double e_plus_jet_errors[mynstage][5][nstage];
	double e_plus_jet_effic[mynstage][5][nmctype];
	double e_plus_jet_effic_unc[mynstage][5][nmctype];

	for (short i = 0; i < mynstage; ++i) {
		for (short j = 0; j < 5; ++j) {
			for (short k = 0; k < 24; ++k) {
				e_plus_jet_errors[i][j][k] = 0;
				e_plus_jet_effic[i][j][k] = 0;
				e_plus_jet_effic_unc[i][j][k] = 0;
			}
		}
	}

	string ttsample = "ttbar";
	if (GetNinit(kttjet) > 0)
		ttsample = "ttjet";
	//	if (GetNinit("TTJet") > 0)
	//		ttsample = "TTJet";

	string wjetSample = "wjet"; //changed from char*
	//  if(nInitialEventMC["WJET"]>0) wjetSample = "WJET";//or whatever the other wjet sample is


	//TODO: is this needed?
	//FIXME: 11th step to new sample - set kindexmcNames
	//	string kIndexmcNames[] = { "", ttsample, ttsample, ttsample, ttsample, ttsample, ttsample, ttsample, ttsample, ttsample,
	//			ttsample, wjetSample, "zjet", "enri1", "enri2", "enri3", "bce1", "bce2", "bce3", "vqq", "tW", "tchan", "schan",
	//			"Zprime_M500GeV_W5GeV", "Zprime_M500GeV_W50GeV", "Zprime_M750GeV_W7500MeV", "Zprime_M1TeV_W10GeV",
	//			"Zprime_M1TeV_W100GeV", "Zprime_M1250GeV_W12500MeV", "Zprime_M1500GeV_W15GeV", "Zprime_M1500GeV_W150GeV",
	//			"Zprime_M2TeV_W20GeV", "Zprime_M2TeV_W200GeV", "Zprime_M3TeV_W30GeV", "Zprime_M3TeV_W300GeV", "Zprime_M4TeV_W40GeV",
	//			"Zprime_M4TeV_W400GeV" };

	for (short i = 0; i < mynstage; ++i) {
		for (short j = 0; j < 5; ++j) {
			for (short k = 1; k < 23; ++k) {

				long ni = 0;
				if (i < 11 && ttsample == "ttbar")
					ni = GetNinit(kttbar);
				else if (i < 11 && ttsample == "ttjet")
					ni = GetNinit(kttjet);
				else
					ni = GetNinit(k - 11 + 3);

				//----------------------
				// Ask FB to check....
				if (ni == 0)
					continue;

				if (e_plus_jet[i][j][k] > 0) {
					e_plus_jet_effic[i][j][k] = e_plus_jet[i][j][k] / ni;
					e_plus_jet_effic_unc[i][j][k] = sqrt(e_plus_jet_effic[i][j][k] * (1 - e_plus_jet_effic[i][j][k]) / ni);
				} else {
					// 0 event pass
					e_plus_jet_effic_unc[i][j][k] = GetBayesUncertainty(ni);
				}
				e_plus_jet_errors[i][j][k] = e_plus_jet_effic_unc[i][j][k] * GetCrossSection(k) * intlumi; //TL 21-8-09 (ask FB to check)
				//----------------------

			}//end of k loop

			e_plus_jet_effic[i][j][23] = 0;
			for (int k = 1; k < 11; ++k) {
				e_plus_jet_effic[i][j][23] += e_plus_jet_effic[i][j][k];
			}
			if (i < 11 && ttsample == "ttbar")
				e_plus_jet_effic_unc[i][j][23] = sqrt(e_plus_jet_effic[i][j][23] * (1 - e_plus_jet_effic[i][j][23]) / GetNinit(
						kttbar));
			else if (i < 11 && ttsample == "ttjet")
				e_plus_jet_effic_unc[i][j][23] = sqrt(e_plus_jet_effic[i][j][23] * (1 - e_plus_jet_effic[i][j][23]) / GetNinit(
						kttjet));

			if (i < 11 && ttsample == "ttbar")
				e_plus_jet_errors[i][j][23] = e_plus_jet_effic_unc[i][j][23] * GetCrossSection(kttbar) * intlumi; //TL 21-8-09 (ask FB to check)
			else if (i < 11 && ttsample == "ttjet")
				e_plus_jet_errors[i][j][23] = e_plus_jet_effic_unc[i][j][23] * GetCrossSection(kttjet) * intlumi; //TL 21-8-09 (ask FB to check)
			//			e_plus_jet_errors[i][j][23] = e_plus_jet_effic_unc[i][j][23] * GetCrossSection(ttsample) * intlumi; //TL 21-8-09 (ask FB to check)

		}
	}

	// QCD sum
	double Sum_Effic2[mynstage][ntjet];
	double Sum_Effic_unc_pos[mynstage][ntjet];
	double Sum_Effic_unc_neg[mynstage][ntjet];
	double Sum_pass_weighted2[mynstage][ntjet];
	// single top sum
	double STopSum_Effic2[mynstage][ntjet];
	double STopSum_Effic_unc_pos[mynstage][ntjet];
	double STopSum_Effic_unc_neg[mynstage][ntjet];
	double STopSum_pass_weighted2[mynstage][ntjet];

	for (short i = 0; i < mynstage; ++i) {
		for (short j = 0; j < 5; ++j) {
			Sum_Effic2[i][j] = 0.0;
			Sum_Effic_unc_pos[i][j] = 0.0;
			Sum_Effic_unc_neg[i][j] = 0.0;
			Sum_pass_weighted2[i][j] = 0.0;

			STopSum_Effic2[i][j] = 0.0;
			STopSum_Effic_unc_pos[i][j] = 0.0;
			STopSum_Effic_unc_neg[i][j] = 0.0;
			STopSum_pass_weighted2[i][j] = 0.0;
		}
	}

	// QCD
	for (short i = 0; i < mynstage; ++i) {
		for (short j = 0; j < 5; ++j) {
			for (short k = 13; k < 19; ++k) {
				Sum_pass_weighted2[i][j] += e_plus_jet_weighted[i][j][k];
				double cer = e_plus_jet_errors[i][j][k];
				Sum_Effic_unc_pos[i][j] += (cer * cer);
				if (e_plus_jet_effic[i][j][k] != 0) {
					Sum_Effic_unc_neg[i][j] += (cer * cer);
				}
			}
			Sum_Effic_unc_pos[i][j] = sqrt(Sum_Effic_unc_pos[i][j]);
			Sum_Effic_unc_neg[i][j] = sqrt(Sum_Effic_unc_neg[i][j]);
		}
	}

	// Single top
	for (int i = 0; i < mynstage; ++i) {
		for (int j = 0; j < 5; ++j) {
			for (int k = 20; k < 23; ++k) {
				STopSum_pass_weighted2[i][j] += e_plus_jet_weighted[i][j][k];
				double cer = e_plus_jet_errors[i][j][k];
				STopSum_Effic_unc_pos[i][j] += (cer * cer);
				if (e_plus_jet_effic[i][j][k] != 0) {
					STopSum_Effic_unc_neg[i][j] += (cer * cer);
				}
			}
			STopSum_Effic_unc_pos[i][j] = sqrt(STopSum_Effic_unc_pos[i][j]);
			STopSum_Effic_unc_neg[i][j] = sqrt(STopSum_Effic_unc_neg[i][j]);
		}
	}

	ofstream myfile;
	myfile.open("Analyzeroutput.txt", ios::app);

	myfile.setf(ios::fixed, ios::floatfield);

	myfile << endl;
	myfile << endl;
	myfile << endl;
	myfile << endl;
	myfile << endl;

	double Allevents[12]; //mynstage
	double AlleventsUncPos[12];
	double AlleventsUncNeg[12];
	double JustSignal[12];
	double JustSignalUnc[12];
	double JustBG[12];
	double JustBGUncPos[12];
	double JustBGUncNeg[12];

	for (short i = 0; i < 12; ++i) {
		Allevents[i] = 0;
		AlleventsUncPos[i] = 0;
		AlleventsUncNeg[i] = 0;
		JustBG[i] = 0;
		JustBGUncPos[i] = 0;
		JustBGUncNeg[i] = 0;
		JustSignal[i] = 0;
		JustSignalUnc[i] = 0;
	}

	myfile.precision(1);

	myfile << "\\begin{tabular}{|l|cccccc|c|}" << endl << "\\hline" << endl;

	myfile << setw(23) << " Cuts  &" << setw(23) << " \\ttbar{}  &" << setw(23) << " W+jets   &" << setw(23) << " Z+jets   &"
			<< setw(23) << " QCD   &" << setw(23) << " VQQ   &" << setw(23) << " Single Top   &" << setw(20) << " Total"
			<< " \\\\ \n\\hline" << endl;

	// First 7 cuts (up to mu-veto)
	for (short i = 0; i < 7; ++i) {
		double TotalEvents = 0;
		double TotalErrorPos = 0;
		double TotalErrorNeg = 0;

		printCutStage(myfile, i, ve.at(i));

		double total = 0;
		double totalerr = 0;

		// signal
		for (short j = 0; j < 5; ++j) {
			for (short k = 1; k < 11; ++k) {
				total += e_plus_jet_weighted[i][j][k];
			}
			totalerr += e_plus_jet_errors[i][j][23] * e_plus_jet_errors[i][j][23];
		}
		printLine(myfile, total, sqrt(totalerr));
		TotalEvents += total;
		TotalErrorPos += totalerr;
		TotalErrorNeg += totalerr;
		JustSignal[i] += total;
		JustSignalUnc[i] += totalerr;

		// Backgrounds
		// W+jets
		total = 0;
		totalerr = 0;
		for (short j = 0; j < 5; ++j) {
			total += e_plus_jet_weighted[i][j][11];
			totalerr += e_plus_jet_errors[i][j][11] * e_plus_jet_errors[i][j][11];
		}
		printLine(myfile, total, sqrt(totalerr));
		TotalEvents += total;
		TotalErrorPos += totalerr;
		TotalErrorNeg += totalerr;
		JustBG[i] += total;
		JustBGUncPos[i] += totalerr;
		JustBGUncNeg[i] += totalerr;

		// Z+jets
		total = 0;
		totalerr = 0;
		for (short j = 0; j < 5; ++j) {
			total += e_plus_jet_weighted[i][j][12];
			totalerr += e_plus_jet_errors[i][j][12] * e_plus_jet_errors[i][j][12];
		}
		printLine(myfile, total, sqrt(totalerr));
		TotalEvents += total;
		TotalErrorPos += totalerr;
		TotalErrorNeg += totalerr;
		JustBG[i] += total;
		JustBGUncPos[i] += totalerr;
		JustBGUncNeg[i] += totalerr;

		total = 0;
		totalerr = 0;
		double totalNerr = 0;
		for (short j = 0; j < 5; ++j) {
			total += Sum_pass_weighted2[i][j];
			totalerr += Sum_Effic_unc_pos[i][j] * Sum_Effic_unc_pos[i][j];
			totalNerr += Sum_Effic_unc_neg[i][j] * Sum_Effic_unc_neg[i][j];
		}
		printLine(myfile, total, sqrt(totalerr));
		TotalEvents += total;
		TotalErrorPos += totalerr;
		TotalErrorNeg += totalNerr;
		JustBG[i] += total;
		JustBGUncPos[i] += totalerr;
		JustBGUncNeg[i] += totalNerr;

		// VQQ
		total = 0;
		totalerr = 0;
		for (short j = 0; j < 5; ++j) {
			total += e_plus_jet_weighted[i][j][19];
			totalerr += e_plus_jet_errors[i][j][19] * e_plus_jet_errors[i][j][19];
		}
		printLine(myfile, total, sqrt(totalerr));

		// Single top
		total = 0;
		totalerr = 0;
		for (short j = 0; j < 5; ++j) {
			total += STopSum_pass_weighted2[i][j];
			totalerr += STopSum_Effic_unc_pos[i][j] * STopSum_Effic_unc_pos[i][j];
		}
		printLine(myfile, total, sqrt(totalerr));
		TotalEvents += total;
		TotalErrorPos += totalerr;
		TotalErrorNeg += totalerr;
		JustBG[i] += total;
		JustBGUncPos[i] += totalerr;
		JustBGUncNeg[i] += totalerr;

		// ALL
		printLine(myfile, TotalEvents, sqrt(TotalErrorPos));
		myfile << " \\\\" << endl;
		Allevents[i] = TotalEvents;
		AlleventsUncPos[i] = TotalErrorPos;
		AlleventsUncNeg[i] = TotalErrorNeg;
	}

	// For cut 4mj to DIFFZ
	for (short i = 6; i < mynstage; ++i) { //cut (up to !DIFFZ)

		double TotalEvents = 0;
		double TotalErrorPos = 0;
		double TotalErrorNeg = 0;

		if (i == 6) {
			printCutStage(myfile, i, "$\\ge$4 jets");
		} else {
			printCutStage(myfile, i, ve.at(i));
		}

		double total = 0;
		for (short k = 1; k < 11; ++k) {
			total += e_plus_jet_weighted[i][4][k];
		}

		printLine(myfile, total, e_plus_jet_errors[i][4][23]);
		TotalErrorPos += e_plus_jet_errors[i][4][23] * e_plus_jet_errors[i][4][23];
		TotalErrorNeg += e_plus_jet_errors[i][4][23] * e_plus_jet_errors[i][4][23];
		TotalEvents += total;
		JustSignal[i + 1] += total;
		JustSignalUnc[i + 1] += e_plus_jet_errors[i][4][23] * e_plus_jet_errors[i][4][23];

		// W+jets
		printLine(myfile, e_plus_jet_weighted[i][4][11], e_plus_jet_errors[i][4][11]);
		TotalErrorPos += e_plus_jet_errors[i][4][11] * e_plus_jet_errors[i][4][11];
		TotalErrorNeg += e_plus_jet_errors[i][4][11] * e_plus_jet_errors[i][4][11];
		TotalEvents += e_plus_jet_weighted[i][4][11];
		JustBG[i + 1] += e_plus_jet_weighted[i][4][11];
		JustBGUncPos[i + 1] += e_plus_jet_errors[i][4][11] * e_plus_jet_errors[i][4][11];
		JustBGUncNeg[i + 1] += e_plus_jet_errors[i][4][11] * e_plus_jet_errors[i][4][11];

		// Z+jets
		printLine(myfile, e_plus_jet_weighted[i][4][12], e_plus_jet_errors[i][4][12]);
		TotalErrorPos += e_plus_jet_errors[i][4][12] * e_plus_jet_errors[i][4][12];
		TotalErrorNeg += e_plus_jet_errors[i][4][12] * e_plus_jet_errors[i][4][12];
		TotalEvents += e_plus_jet_weighted[i][4][12];
		JustBG[i + 1] += e_plus_jet_weighted[i][4][12];
		JustBGUncPos[i + 1] += e_plus_jet_errors[i][4][12] * e_plus_jet_errors[i][4][12];
		JustBGUncNeg[i + 1] += e_plus_jet_errors[i][4][12] * e_plus_jet_errors[i][4][12];

		// QCD
		printLine(myfile, Sum_pass_weighted2[i][4], Sum_Effic_unc_pos[i][4]);
		TotalErrorPos += Sum_Effic_unc_pos[i][4] * Sum_Effic_unc_pos[i][4];
		TotalErrorNeg += Sum_Effic_unc_neg[i][4] * Sum_Effic_unc_neg[i][4];
		TotalEvents += Sum_pass_weighted2[i][4];
		JustBG[i + 1] += Sum_pass_weighted2[i][4];
		JustBGUncPos[i + 1] += Sum_Effic_unc_pos[i][4] * Sum_Effic_unc_pos[i][4];
		JustBGUncNeg[i + 1] += Sum_Effic_unc_neg[i][4] * Sum_Effic_unc_neg[i][4];

		// VQQ
		printLine(myfile, e_plus_jet_weighted[i][4][19], e_plus_jet_errors[i][4][19]);

		// Single top
		printLine(myfile, STopSum_pass_weighted2[i][4], STopSum_Effic_unc_pos[i][4]);
		TotalErrorPos += STopSum_Effic_unc_pos[i][4];
		TotalErrorNeg += STopSum_Effic_unc_pos[i][4];
		TotalEvents += STopSum_pass_weighted2[i][4];
		JustBG[i + 1] += STopSum_pass_weighted2[i][4];
		JustBGUncPos[i + 1] += STopSum_Effic_unc_pos[i][4];
		JustBGUncNeg[i + 1] += STopSum_Effic_unc_neg[i][4];

		// All
		Allevents[i + 1] += TotalEvents;
		AlleventsUncPos[i + 1] += TotalErrorPos;
		AlleventsUncNeg[i + 1] += TotalErrorNeg;
		printLine(myfile, TotalEvents, sqrt(TotalErrorPos));
		myfile << " \\\\" << endl;
	}
	myfile << "  \\hline" << endl << "\\end{tabular}" << endl;
	myfile << endl << endl << endl;

	// Break Down of Single Top (with errors)
	//----------------------------------------
	myfile << "\\begin{tabular}{|l|rrr|r|}" << endl;
	myfile << "\\hline" << endl;
	myfile << "\\multicolumn{5}{|l|}";
	myfile << "{Break down of expected single top events for " << (int) intlumi << "/pb}";
	myfile << "\\\\\n\\hline" << endl;

	myfile << "            Cut      ";
	myfile << " &" << setw(21) << "tW-chan  ";
	myfile << " &" << setw(21) << "t-chan  ";
	myfile << " &" << setw(21) << "s-chan  ";
	myfile << " &" << setw(29) << "AllSingleTop \\\\\\hline" << endl;

	short njbegin = 0;

	for (short i = 0; i < 11; ++i) { //cut stage (up to HT)
		double allSingleTop = 0;
		double allSingleTopEr = 0;
		printCutStage(myfile, i, ve.at(i));
		for (short k = 20; k < 23; ++k) { //mctype (single top): 20-22
			double totalT = 0;
			double totalEr = 0;
			for (short j = njbegin; j < ntjet; ++j) { //njet
				totalT += e_plus_jet_weighted[i][j][k];
				totalEr += e_plus_jet_errors[i][j][k] * e_plus_jet_errors[i][j][k];
			}
			totalEr = sqrt(totalEr);
			allSingleTop += totalT;

			printLine(myfile, totalT, totalEr);
		}
		for (short j = njbegin; j < ntjet; ++j) {
			allSingleTopEr += STopSum_Effic_unc_pos[i][j] * STopSum_Effic_unc_pos[i][j];
		}
		allSingleTopEr = sqrt(allSingleTopEr);

		myfile << " & " << setw(10) << ScrNum(allSingleTop) << "$\\pm$" << setw(4) << left << ScrNum(allSingleTopEr) << right
				<< " \\\\" << endl;

		// insert >=4j cut
		if (ve.at(i) == "!MUON") {
			njbegin = 4;
			ve.at(i) = "$\\ge$4 jets";
			i--;
		}
	}

	myfile << "\\hline" << endl;
	myfile << "\\end{tabular}\\\\[5mm]" << endl;
	myfile << endl << endl;
	ve.at(6) = "!MUON"; //restore


	//BEGIN QCD BREAKDOWN TABLE
	//-------------------------
	myfile << "\\begin{tabular}{|l|cccccc|c|}" << endl;
	myfile << "\\hline" << endl;
	myfile << "\\multicolumn{8}{|l|}";
	myfile << "{Break down of expected QCD events for " << intlumi << "/pb}";
	myfile << "\\\\\n\\hline" << endl;

	myfile << "            Cut      ";
	myfile << " &" << setw(24) << "enri1 ";
	myfile << " &" << setw(24) << "enri2 ";
	myfile << " &" << setw(24) << "enri3 ";
	myfile << " &" << setw(24) << "bce1 ";
	myfile << " &" << setw(24) << "bce2 ";
	myfile << " &" << setw(24) << "bce3 ";
	myfile << " &" << setw(27) << "AllQCD \\\\\\hline" << endl;

	//ntjet = 5;
	njbegin = 0;

	for (short i = 0; i < 11; ++i) { //cut stage (up to HT)
		double totalAllQCD = 0;
		double totalAllQCDErPos = 0;
		double totalAllQCDErNeg = 0;
		//myfile << "" << setw(2) << i << " " << setw(11) << left << ve.at(i) << right;
		printCutStage(myfile, i, ve.at(i));

		for (short k = 13; k < 19; ++k) { //mctype (QCD): 13-18
			double totalT = 0;
			double totalEr = 0;
			for (short j = njbegin; j < ntjet; ++j) { //njet
				totalT += e_plus_jet_weighted[i][j][k];
				totalEr += e_plus_jet_errors[i][j][k] * e_plus_jet_errors[i][j][k];
			}
			totalEr = sqrt(totalEr);
			totalAllQCD += totalT;

			myfile << " &" << setw(12) << ScrNum(totalT);
			if (totalT != 0)
				myfile << "$\\pm$" << setw(7) << left << ScrNum(totalEr) << right;
			else {
				myfile << "$<$" << setw(9) << left << ScrNum(totalEr) << right;
			}
		}
		for (short j = njbegin; j < ntjet; ++j) {
			totalAllQCDErPos += Sum_Effic_unc_pos[i][j] * Sum_Effic_unc_pos[i][j];
			totalAllQCDErNeg += Sum_Effic_unc_neg[i][j] * Sum_Effic_unc_neg[i][j];
		}

		totalAllQCDErPos = sqrt(totalAllQCDErPos);
		totalAllQCDErNeg = sqrt(totalAllQCDErNeg);
		myfile << " & " << setw(13) << ScrNum(totalAllQCD) << "$\\pm$" << setw(8) << left << ScrNum(totalAllQCDErPos) << right
				<< " \\\\" << endl;

		if (ve.at(i) == "!MUON") {
			njbegin = 4;
			ve.at(i) = "$\\ge$4 jets";
			i--;
		}
	}
	myfile << "\\hline" << endl;
	myfile << "\\end{tabular}\\\\[5mm]\n" << endl << endl;

	ve.at(6) = "!MUON";//restore


	myfile << "\n%---------------------------------------------------------------------\n";
	myfile << "       Expected Signal and Background for " << intlumi << "/pb";
	myfile << "\n%---------------------------------------------------------------------\n\n";
	myfile << "\\begin{tabular}{|l|r|rr|r|}\\hline" << endl;
	myfile << "          Cut        " << " &  " << setw(30) << "Total Events " << " &  " << setw(24) << "Total Signal (S) "
			<< " &  " << setw(26) << "Total Background (B) " << " &  " << setw(20) << "S/B   \\\\\n\\hline" << endl;

	// Insert >=4j cut after muon-veto
	const short muVeto_pos = 7;
	ve.insert(ve.begin() + muVeto_pos, "$\\ge$4 jets");

	for (short i = 0; i <= 11; ++i) { //loop over cuts (up to DIFFZ)

		printCutStage(myfile, i, ve.at(i));

		myfile << " & " << setw(12) << right << ScrNum(Allevents[i]) << setw(24) << left << Form("$^{+%s}_{-%s}$", ScrNum(sqrt(
				AlleventsUncPos[i])).c_str(), ScrNum(sqrt(AlleventsUncNeg[i])).c_str()) << left;
		myfile << " & " << setw(10) << right << ScrNum(JustSignal[i]) << setw(10) << left << "$\\pm$" + ScrNum(sqrt(
				JustSignalUnc[i]));
		myfile << " & " << setw(12) << right << ScrNum(JustBG[i]) << setw(24) << left << Form("$^{+%s}_{-%s}$", ScrNum(sqrt(
				JustBGUncPos[i])).c_str(), ScrNum(sqrt(JustBGUncNeg[i])).c_str()) << right;

		//Print Signal-to-background ratio (S/B)
		if (JustBG[i] > 0) {
			double errSBGPos = sqrt(JustSignalUnc[i] / (JustSignal[i] * JustSignal[i]) + JustBGUncPos[i]
					/ (JustBG[i] * JustBG[i]));
			myfile << " & " << setw(11) << ScrNum(JustSignal[i] / JustBG[i]) << "$\\pm$" << ScrNum(errSBGPos * JustSignal[i]
					/ JustBG[i]);
		} else {
			myfile << " & " << setw(12) << "-";
		}
		myfile << " \\\\" << endl;

	}// end loop over cut (nstage)
	myfile << "\\hline\n\\end{tabular}\n" << endl << endl << endl;

	// restore (take out 4mj)
	ve.erase(ve.begin() + muVeto_pos);

	myfile << "%---------------------------------------" << endl;
	myfile << "            NjetVcut table " << endl;
	myfile << "%---------------------------------------\n" << endl;
	PrintError_NjetVcut(myfile, e_plus_jet_weighted, e_plus_jet_errors, ve);

	myfile.close();

}//end PrintErrorTables

void ana::printLine(ofstream &myfile, double num, double err) const {
	myfile << " &" << setw(9) << ScrNum(num) << "$\\pm$" << setw(7) << left << ScrNum(err) << right;
}

string ana::ScrNum(double num) const {

	string Number = "";
	if (ScientificNotation) {
		if (num / 10000.0 > 1) {
			Number = Form("%.1E", num);
		} else {
			Number = Form("%.1f", num);
		}
	} else {
		Number = Form("%.1f", num);
	}
	return Number;
}

//---------------------------------------------------------------------
double ana::GetBayesUncertainty(int Ninitial) const {

	if (Ninitial == 0)
		return 0; //NEW, added by TL (Ask FB to check)

	TH1D *h_total = new TH1D("total", "events that are incident", 1, 0, 1);
	h_total->SetBinContent(1, Ninitial);
	TH1D *h_pass = new TH1D("pass", "events that pass", 1, 0, 1);
	h_pass->SetBinContent(1, 0);

	TGraphAsymmErrors *tge = new TGraphAsymmErrors();
	tge->BayesDivide(h_pass, h_total);
	delete h_total;
	delete h_pass;
	double error1 = tge->GetErrorYhigh(0);
	delete tge;
	return error1;

}
//---------------------------------------------------------------------


//----------------------------------------------------------------------------------
void ana::PrintError_NjetVcut(ofstream& myfile, const double e_plus_jet_weighted[][5][nmctype],
		const double e_plus_jet_errors[][5][nstage], vector<string>& ve) const {

	myfile << "\\begin{tabular}{|c|ccccc|c|}" << endl;
	myfile << "\\hline" << endl;

	myfile << setw(21) << "Cut  " << " &" << setw(26) << "0-jet " << " &" << setw(26) << "1-jet " << " &" << setw(26)
			<< "2-jets " << " &" << setw(26) << "3-jets " << " &" << setw(26) << "$\\ge$4-jets " << " &" << setw(36)
			<< "Total   \\\\\n\\hline\n";

	for (short i = 0; i < 11; ++i) { //nstage

		printCutStage(myfile, i, ve.at(i));

		double total = 0;
		double PosError = 0;
		double NegError = 0;
		for (short j = 0; j < ntjet; ++j) {
			double myTotal = 0;
			double myPosErr = 0;
			double myNegError = 0;
			for (short k = 1; k < 23; ++k) {
				if (k == 19)
					continue;//skip vqq
				myTotal += e_plus_jet_weighted[i][j][k];
				myPosErr += e_plus_jet_errors[i][j][k] * e_plus_jet_errors[i][j][k];
				if (e_plus_jet_weighted[i][j][k] != 0) {
					myNegError += e_plus_jet_errors[i][j][k] * e_plus_jet_errors[i][j][k];
				}
			}
			myfile << " & " << setw(12) << ScrNum(myTotal) << "$\\pm$" << setw(8) << left << ScrNum(sqrt(myPosErr)) << right;
			total += myTotal;
			PosError += myPosErr;
			NegError += myNegError;
		}//loop jets
		myfile << " & " << setw(13) << fixed << ScrNum(total) << "$\\pm$" << setw(8) << left << ScrNum(sqrt(PosError))
				<< " \\\\ \n" << right;
	}//loop cuts
	myfile << "\n\\hline\n\\end{tabular}" << endl;

}//end PrintError_NjetVcut
//----------------------------------------------------------------------------------


float ana::calcDeltaPhi(const TLorentzVector& p1, const TLorentzVector& p2) const {
	return p1.DeltaPhi(p2);
}

float ana::calcDeltaEta(const TLorentzVector& p1, const TLorentzVector& p2) const {
	float delEta(999.0);
	delEta = p1.PseudoRapidity() - p2.PseudoRapidity();
	return delEta;
}

// Calculate Delta R
float ana::calcDeltaR(const TLorentzVector& p1, const TLorentzVector& p2) const {
	return p1.DeltaR(p2);
}

//Calculate Delta R
float ana::calcDeltaR(const float phi1, const float eta1, const float phi2, const float eta2) const {

	// Calculate delta phi and delta eta
	float delPhi = fabs(phi1 - phi2);
	if (delPhi > TMath::Pi())
		delPhi = 2.0 * TMath::Pi() - delPhi;
	float delEta = eta1 - eta2;

	// Calculate the delta R
	float delR = sqrt(delPhi * delPhi + delEta * delEta);

	return delR;
}

//---------------------------------------------------------------------------------------------
//                                   QCD Estimation
//---------------------------------------------------------------------------------------------
// global var
//const float QCDest_reliso_bin_width = 0.1;

void ana::EstimateQCD() {
	EstimateQCD(outputHistFileName);
}

bool ana::EstimateQCD(const string inputFile) {

	cout << "--------------------------------------------------------" << endl;
	cout << "                    QCD Estimation" << endl;
	cout << "--------------------------------------------------------" << endl;

	string isoVariable = "QCDest_CombRelIso"; //new reliso
	if (!useNewReliso)
		isoVariable = "QCDest_NormCombRelIso"; //old reliso

	// Open the histogram file to add histogram
	TFile f(inputFile.c_str(), "UPDATE");
	outfile = fopen(outputTextFileName.c_str(), "a"); //append

	TCanvas c1("c1", "Electron isolation", 980, 700);
	c1.Divide(2, 2);

	//f.Print();
	const int nrange = 1;// number of trial fit ranges

	const int nj = 4;

	double sig_from = 0; //new reliso
	double sig_upto = 0.1;
	if (!useNewReliso) { //old reliso
		sig_from = 0.85;
		sig_upto = 1.0;
	}
	if (useNewReliso)
		cout << "Using New Reliso formulation" << endl;
	else
		cout << "Using Old Reliso formulation" << endl;
	cout << "Signal region: " << sig_from << " to " << sig_upto << endl;

	// use single function with fixed range
	const double fit_range_low_oldReliso[nj] = { 0.4, 0.4, 0.4, 0.4 };
	const double fit_range_up_oldReliso[nj] = { 0.8, 0.8, 0.8, 0.8 };

	const double fit_range_low_newReliso[nj] = { 0.2, 0.2, 0.2, 0.2 }; //<-----
	const double fit_range_up_newReliso[nj] = { 1.0, 1.0, 1.0, 1.0 }; //<-----

	double fit_range_low[nj];
	double fit_range_up[nj];

	int rebin_oldReliso[nj] = { 1, 2, 5, 10 };
	int rebin_newReliso[nj] = { 10, 10, 10, 10 }; //<----
	int rebin[nj];

	for (int i = 0; i < nj; ++i) { //set the fit range according to reliso choice
		if (useNewReliso) {
			fit_range_low[i] = fit_range_low_newReliso[i];
			fit_range_up[i] = fit_range_up_newReliso[i];
			rebin[i] = rebin_newReliso[i];
		} else {
			fit_range_low[i] = fit_range_low_oldReliso[i];
			fit_range_up[i] = fit_range_up_oldReliso[i];
			rebin[i] = rebin_oldReliso[i];
		}
	}

	string func = "gaus";
	if (!useNewReliso)
		func = "gaus"; //old reliso
	cout << "Fit function: " << func << endl;

	//all=s+b, for 0j,1j,2j,3j,4j,>=4j
	double nall_actual_ctr[nrange][nj];
	double nqcd_actual_sig[nrange][nj]; //qcd only
	double nqcd_actual_ctr[nrange][nj];
	double nqcd_actual_all[nj];
	double nqcd_actual_all_unweighted[nj];
	double n_extrap[nrange][nj];
	double n_extrap_err_plus[nrange][nj];
	double n_extrap_err_minus[nrange][nj];
	double fit_chi2[nrange][nj];
	int fit_ndf[nrange][nj];

	// initialize arrays to zero
	for (int i = 0; i < nrange; ++i) {
		for (int j = 0; j < nj; ++j) {
			nall_actual_ctr[i][j] = 0;
			nqcd_actual_sig[i][j] = 0;
			nqcd_actual_ctr[i][j] = 0;
			nqcd_actual_all[j] = 0;
			nqcd_actual_all_unweighted[j] = 0;
			n_extrap[i][j] = 0;
			n_extrap_err_plus[i][j] = 0;
			n_extrap_err_minus[i][j] = 0;
			fit_chi2[i][j] = 0;
			fit_ndf[i][j] = 0;
		}
	}

	// take CombRelIso distributions, fit in range 0-0.85, extrapolate to 0.85-1.0
	const string myjetbin[4] = { "1 jet", "2 jets", "3 jets", "#geq4 jets" }; //root label
	const string jetlabel[4] = { "1j", "2j", "3j", "$\\ge$4j" }; //latex label
	const string jettext[4] = { "1j", "2j", "3j", ">=4j" }; //text label

	for (int i = 0; i < nrange; ++i) {

		// 27 Feb 09
		vector<double> mpv;

		for (int j = 0; j < nj; ++j) { //nj: 0,1,2,3,4,>=4j -> 1,2,3,>=4j


			const double fit_from = fit_range_low[j];
			const double fit_upto = fit_range_up[j];
			cout << "\n\n-----------------------" << endl;
			cout << " " << jettext[j] << "  fit range: " << fit_from << " to " << fit_upto << endl;

			char *hname = Form("%s_%dj", isoVariable.c_str(), j + 1);
			if (j + 1 == 4)
				hname = Form("%s_4mj", isoVariable.c_str());
			TH1F *allEvent = (TH1F*) f.Get(Form("data/QCD_estimation/%s", hname)); //all events
			TH1F *qcdEvent = (TH1F*) f.Get(Form("QCD/QCD_estimation/%s", hname)); //QCD only
			if (allEvent == 0 || qcdEvent == 0) {
				cout << "Warning: zero pointer to reliso histograms, exit EstimateQCD()" << endl;
				return false;
			}

			const double binw = 0.01; //original binw = 0.01 = 1.1/110 (old) = 10/1000 (new)
			// signal bin (86, 100) for old reliso
			// signal bin ( 1, 100) for new reliso
			const int ctr_bin_low = (int) (fit_from / binw + 1);
			const int ctr_bin_up = (int) (fit_upto / binw);
			const int sig_bin_low = (int) (sig_from / binw + 1);
			const int sig_bin_up = (int) (sig_upto / binw);
			cout << "     signal bin:  " << sig_bin_low << " to " << sig_bin_up;
			cout << "  control bin:  " << ctr_bin_low << " to " << ctr_bin_up << endl;

			nall_actual_ctr[i][j] = allEvent->Integral(ctr_bin_low, ctr_bin_up); //actual number of _all_ event in control region

			nqcd_actual_ctr[i][j] = qcdEvent->Integral(ctr_bin_low, ctr_bin_up); //actual number of _QCD_ event in control region
			nqcd_actual_sig[i][j] = qcdEvent->Integral(sig_bin_low, sig_bin_up); //actual number of _QCD_ event in signal region

			nqcd_actual_all_unweighted[j] = qcdEvent->GetEntries();
			nqcd_actual_all[j] = qcdEvent->Integral();
			nqcd_actual_ctr[i][j] = qcdEvent->Integral(ctr_bin_low, ctr_bin_up);

			const int this_rebin = rebin[j];

			// set value of global variable
			//QCDest_reliso_bin_width = binw*this_rebin; //binw = original bin width = 0.01
			Set_Reliso_bin_width(binw * this_rebin); //binw = original bin width = 0.01
			cout << "     rebin = " << this_rebin << ", binw = " << Get_Reliso_bin_width() << endl;

			TH1F* all = (TH1F*) allEvent->Clone();
			TH1F* qcd = (TH1F*) qcdEvent->Clone();

			all->Rebin(this_rebin); //all events
			qcd->Rebin(this_rebin); //qcd only


			// Draw
			c1.cd(j + 1);
			qcd->SetFillColor(kAzure - 9); //blue
			qcd->SetLineColor(kAzure + 2); //blue (darker)
			all->Draw("ahist");
			if (useNewReliso)
				all->GetXaxis()->SetRangeUser(0, 1.6);
			all->SetMinimum(0);
			all->SetTitle(Form("RelIso distribution (%s)", myjetbin[j].c_str()));
			all->GetXaxis()->SetTitle("RelIso");

			// Do the fit if there are non-0 event in control region

			if (nall_actual_ctr[i][j] > 0) {

				// For 3j, >=4j, set limit on MPV (for Landau only)
				//---------------------------------
				if (useNewReliso && j > 1 && func == "landau") { //j=2 is 3j

					// MPV-Constrained Landau Fit
					TF1 *fitf = new TF1("landau", "landau", 0, 2);

					// Set limit on MPV parameter for 3,>=4 jet bins
					if (mpv.size() >= 2) {

						cout << "mpv fitted: 1j= " << mpv.at(0) << ", 2j=" << mpv.at(1) << endl;
						if (mpv.at(0) <= mpv.at(1)) {
							fitf->SetParLimits(1, mpv.at(0), mpv.at(1));
						} else {
							fitf->SetParLimits(1, mpv.at(1), mpv.at(0));
						}
					} else {
						cout << "\n could not find 2 mpv values from 1j,2j fits.\n" << endl;
					}

					all->Fit(fitf, "V", "ah", fit_from, fit_upto); //"Q": Quiet, "0": do not draw, "ah": (no axis)
					delete fitf;//test
				} else {
					// Free Fit
					all->Fit(func.c_str(), "V", "ah", fit_from, fit_upto); //"Q": Quiet, "0": do not draw, "ah": (no axis)
				}

				qcd->Draw("ahist same"); //qcd part, no axis
				qcd->Draw("ae same");
				all->Draw("ae same");
				all->Draw("axis same");

				TF1 *myf = all->GetFunction(func.c_str());
				myf->SetLineColor(kRed);
				myf->Draw("same");

				TF1 *myf2 = (TF1*) myf->Clone();
				myf2->SetLineColor(kBlue);
				myf2->SetRange(sig_from, sig_upto); //new
				myf2->Draw("same");

				// resize the stat box
				gPad->Update();
				TPaveStats *st = (TPaveStats*) all->FindObject("stats");
				st->SetX1NDC(0.6);
				st->SetY1NDC(0.5);
				st->SetX2NDC(0.97);
				st->SetY2NDC(0.97);
				st->Draw();

				// Get fit results
				//------------------
				fit_chi2[i][j] = myf->GetChisquare();
				fit_ndf[i][j] = myf->GetNDF();
				// Extrapolate fitted function

				n_extrap[i][j] = myf->Integral(sig_from, sig_upto) / Get_Reliso_bin_width();

				// calculate uncertainty on n_extrap by varying the fit function
				if (func == "gaus" || func == "landau") {
					pair<double, double> estimate_error = estimateQCD_computeFitResultUncertainty(n_extrap[i][j], myf);
					n_extrap_err_plus[i][j] = estimate_error.first;
					n_extrap_err_minus[i][j] = estimate_error.second;
				}

				// Get fit parameters
				if (func == "landau") {
					mpv.push_back(myf->GetParameter(1));
				}

				// Add legend
				//-------------
				TLegend *leg = new TLegend(0.22, 0.55, 0.6, 0.88);
				leg->SetFillStyle(0);
				leg->SetBorderSize(0);
				leg->AddEntry(all, "All events (s+b)");
				leg->AddEntry(qcd, "QCD events", "f");
				leg->AddEntry(myf, "Fit to all events");
				leg->AddEntry(myf2, "Extrapolation");
				leg->Draw();

			}// end fit
		}// end nj loop: 0 to 4 jets

		cout << "-----------------------------------------------" << endl;
		cout << "  " << intlumi << "/pb     rb      True    Estimate     Diff" << endl;
		printf("   1 jet:   %2d %10.1f  %10.1f  %6.1f %%\n", rebin[0], nqcd_actual_sig[i][0], n_extrap[i][0], (n_extrap[i][0]
				/ nqcd_actual_sig[i][0] - 1) * 100);
		printf("   2 jet:   %2d %10.1f  %10.1f  %6.1f %%\n", rebin[1], nqcd_actual_sig[i][1], n_extrap[i][1], (n_extrap[i][1]
				/ nqcd_actual_sig[i][1] - 1) * 100);
		printf("   3 jet:   %2d %10.1f  %10.1f  %6.1f %%\n", rebin[2], nqcd_actual_sig[i][2], n_extrap[i][2], (n_extrap[i][2]
				/ nqcd_actual_sig[i][2] - 1) * 100);
		printf(" >=4 jet:   %2d %10.1f  %10.1f  %6.1f %%\n", rebin[3], nqcd_actual_sig[i][3], n_extrap[i][3], (n_extrap[i][3]
				/ nqcd_actual_sig[i][3] - 1) * 100);
		cout << "-----------------------------------------------" << endl;

		for (short j = 0; j < 4; ++j) {
			cout << "Unc of QCD estimate (" << jettext[j] << "):  +" << n_extrap_err_plus[i][j] << " / -"
					<< n_extrap_err_minus[i][j] << " events" << endl;
		}

		cout << "\nEstimated number (>=4 jets): " << n_extrap[0][3] << " events";
		cout << " (+" << n_extrap_err_plus[i][3] << "/-" << n_extrap_err_minus[i][3] << " events)" << endl;

		// save canvas (for each fit range)
		c1.SaveAs(Form("c1_reliso_fit_%d.ps", i + 1)); // hard-coded filename for now, can easily be changed if we want

	}// end nrange loop
	c1.Close();

	// Print result table
	//---------------------
	fprintf(outfile, "\n%%- - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
	fprintf(outfile, "%%    QCD Estimation Results\n");
	fprintf(outfile, "%%- - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
	fprintf(outfile, "%%\\scriptsize\n"); //small font size
	fprintf(outfile, "\\addtolength{\\tabcolsep}{-0.9mm}\n");
	fprintf(outfile, "\\begin{tabular}{|r|rr|rr|rr|rcccc|}\\hline\n");
	fprintf(outfile, " & \\multicolumn{2}{c|}{Total QCD}\n");
	fprintf(outfile, " & \\multicolumn{2}{c|}{Control region}\n");
	fprintf(outfile, " & \\multicolumn{2}{c|}{Signal region (QCD)}\n");
	fprintf(outfile, " & \\multicolumn{5}{c|}{Fit results} \\\\\n");
	fprintf(outfile, "         &   True &     True &  True (All) &  True (QCD) & True");
	fprintf(outfile, " & Estimate & Diff & $\\chi^2$/dof & Function & Range & Bin\\\\\n");
	fprintf(outfile,
			"         & (MC) & (%2.0f/pb) & (%2.0f/pb) & (%2.0f/pb) & (%2.0f/pb) & (%2.0f/pb) &&&&&width \\\\ \\hline\n",
			intlumi, intlumi, intlumi, intlumi, intlumi);

	for (int i = 0; i < nrange; ++i) { //range

		for (int j = 0; j < nj; ++j) {
			double diff = -1;
			if (nqcd_actual_sig[i][j] > 0)
				diff = (n_extrap[i][j] / nqcd_actual_sig[i][j] - 1) * 100;

			fprintf(outfile, "%8s & ", jetlabel[j].c_str()); // %10.4g = 4 sig figures
			fprintf(outfile, "%6.0f & %9.1f & ", nqcd_actual_all_unweighted[j], nqcd_actual_all[j]); //total QCD
			fprintf(outfile, "%8.1f & %8.1f & ", nall_actual_ctr[i][j], nqcd_actual_ctr[i][j]); //control region: all, QCD
			fprintf(outfile, "%8.1f & ", nqcd_actual_sig[i][j]); //signal region: actual QCD

			if (nall_actual_ctr[i][j] > 0) {
				fprintf(outfile, "%8.1f & ", n_extrap[i][j]); //QCD estimate
				fprintf(outfile, "%5.1f \\%% & %6.2f & ", diff, fit_chi2[i][j] / fit_ndf[i][j]);
				fprintf(outfile, "%6s & %2.1f--%2.2f & %5.2f\\\\\n", func.c_str(), fit_range_low[j], fit_range_up[j], rebin[j]
						/ (110 / 1.1));

			} else {
				//if zero actual event, print "-" in the extrapolation result column
				fprintf(outfile, "%8s & %8s & %6s & %6s & %9s & -\\\\\n", "-", "-", "-", "-", "-");
			}

		}// end nj loop
		fprintf(outfile, "\\hline\n");
	}// end fit-ranges loop

	fprintf(outfile, "\\end{tabular}\n");
	fprintf(outfile, "%%- - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
	fclose(outfile);

	cout << "\n  Actual number (>=4 jets): " << nqcd_actual_sig[0][3] << " +/- " << sqrt(nqcd_actual_sig[0][3])
			<< " events.  (unc=sqrt(N))" << endl;

	//add QCD estimate results to root file
	TH1D *QCDEstimate = new TH1D("QCDEstimate", "QCD estimate results (>=4j)", 8, 0, 8);
	// (i) range 1: 0 to 0.85
	QCDEstimate->SetBinContent(1, n_extrap[0][3]);
	QCDEstimate->SetBinContent(2, n_extrap_err_plus[0][3]);
	QCDEstimate->SetBinContent(3, n_extrap_err_minus[0][3]);
	//average error = (plus + minus)/2
	QCDEstimate->SetBinContent(4, (n_extrap_err_plus[0][3] + n_extrap_err_minus[0][3]) / 2); //average
	QCDEstimate->Write();
	f.Close();

	cout << "--------------------------------------------------------" << endl;
	cout << "                    QCD Estimation (end)" << endl;
	cout << "--------------------------------------------------------" << endl;
	return true;

}// end of EstimateQCD()

//-------------------------------------------------------------------------------------------
// Method to calculate QCD estimate uncertainty, by varying each of the 3 parameters
// in the Gaussian fit (ie constant, mean, sigma), by 1sigma plus and minus.
// The asymmetric uncertainty is calculated by adding in quadrature individual deviations.
// -> adapt for the case when we use landau to fit (also 3 parameters: constant, MPV, sigma).
//-------------------------------------------------------------------------------------------
pair<double, double> ana::estimateQCD_computeFitResultUncertainty(const double est, TF1 *myf) const {

	// get fit results
	double constant = myf->GetParameter(0); //1st parameter: constant
	double mean = myf->GetParameter(1); //2nd parameter: mean (gaus) or MPV (landau)
	double sigma = myf->GetParameter(2); //3rd parameter: sigma
	double constant_err = myf->GetParError(0);
	double mean_err = myf->GetParError(1);
	double sigma_err = myf->GetParError(2);
	double est_const_pos;
	double est_const_neg;
	double est_mean_pos;
	double est_mean_neg;
	double est_sigma_pos;
	double est_sigma_neg;
	double new_est[2];

	// 1) vary constant, get new estimate
	new_est[0] = estimateQCD_newEstimateByVaryingFitFunction(constant + constant_err, mean, sigma);
	new_est[1] = estimateQCD_newEstimateByVaryingFitFunction(constant - constant_err, mean, sigma);
	pair<double, double> est_mean_pair = estimateQCD_assign_pos_neg_estimate(est, new_est);
	est_const_pos = est_mean_pair.first;
	est_const_neg = est_mean_pair.second;

	// 2) vary mean(or MPV), get new estimate
	new_est[0] = estimateQCD_newEstimateByVaryingFitFunction(constant, mean + mean_err, sigma);
	new_est[1] = estimateQCD_newEstimateByVaryingFitFunction(constant, mean - mean_err, sigma);
	est_mean_pair = estimateQCD_assign_pos_neg_estimate(est, new_est);
	est_mean_pos = est_mean_pair.first;
	est_mean_neg = est_mean_pair.second;

	// 3) vary sigma, get new estimate
	new_est[0] = estimateQCD_newEstimateByVaryingFitFunction(constant, mean, sigma + sigma_err);
	new_est[1] = estimateQCD_newEstimateByVaryingFitFunction(constant, mean, sigma - sigma_err);
	est_mean_pair = estimateQCD_assign_pos_neg_estimate(est, new_est);
	est_sigma_pos = est_mean_pair.first;
	est_sigma_neg = est_mean_pair.second;

	// 4) Compute total uncertainty (sum in quadrature)
	double unc_pos = sqrt(TMath::Power((est_const_pos - est), 2) + TMath::Power((est_mean_pos - est), 2) + TMath::Power(
			(est_sigma_pos - est), 2));
	double unc_neg = sqrt(TMath::Power((est_const_neg - est), 2) + TMath::Power((est_mean_neg - est), 2) + TMath::Power(
			(est_sigma_neg - est), 2));
	pair<double, double> unc(unc_pos, unc_neg);
	return unc;
}

double ana::estimateQCD_newEstimateByVaryingFitFunction(const double this_constant, const double this_mean,
		const double this_sigma) const {

	// 1) create new function using each error of fit parameters
	TF1 *f1a;
	if (useNewReliso)
		f1a = new TF1("f1a", "landau", 0, 1.1);
	else
		f1a = new TF1("f1a", "gaus", 0, 1);

	f1a->SetParameter("Constant", this_constant);
	if (useNewReliso)
		f1a->SetParameter("MPV", this_mean); //landau
	else
		f1a->SetParameter("Mean", this_mean); //gaus
	f1a->SetParameter("Sigma", this_sigma);

	double new_est = 0;
	if (useNewReliso)
		new_est = f1a->Integral(0, 0.1) / Get_Reliso_bin_width();
	else
		new_est = f1a->Integral(0.85, 1.0) / Get_Reliso_bin_width();

	delete f1a;
	return new_est;
}
// Method to check whether the new estimate is larger or smaller than the nominal estimate
pair<double, double> ana::estimateQCD_assign_pos_neg_estimate(const double est, const double new_est[2]) const {

	double pos = est;
	double neg = est;

	if (new_est[0] > est && new_est[1] < est) {
		pos = new_est[0];
		neg = new_est[1];
	} else if (new_est[0] < est && new_est[1] > est) {
		pos = new_est[1];
		neg = new_est[0];
	} else {
		if (new_est[0] > est && new_est[1] > est) { //+ve side
			if (new_est[0] >= new_est[1]) {
				pos = new_est[0];
			} else {
				pos = new_est[1];
			}
		}
		if (new_est[0] < est && new_est[1] < est) { //-ve side
			if (new_est[0] <= new_est[1]) {
				neg = new_est[0];
			} else {
				neg = new_est[1];
			}
		}
	}
	pair<double, double> retval(pos, neg);
	return retval;
}
//----------------------------------------------------------------------------------------
//                                  QCD Estimation (end)
//----------------------------------------------------------------------------------------


//81FB
void ana::valid_mkHisto_cut_njet(TH1F* h[][7], const string name, const string title, const int nbin, const float xlow,
		const float xup) {

	const string jetname[7] = { "0j", "1j", "2j", "3j", "4j", "4mj", "allj" };
	const string jetlabel[7] = { "0j", "1j", "2j", "3j", "4j", ">=4j", "allj" };
	const string cutname[9] = { "Initial", "Trigger", ">=1T ele", ">=1T Iso", "==1T Iso", "noMuon", "MET", "Z veto", "noConv" };
	const string cutlabel[9] = { "Initial", "Trigger", ">=1Tele", ">=1TIso", "==1TIso", "noMuon", "MET", "Zveto", "noConv" };

	char hname[70];
	char htitle[100];
	for (short i = 0; i < 9; ++i) { //cut
		for (short j = 0; j < 7; ++j) { //nj
			sprintf(hname, "%s_%s_%s", name.c_str(), cutlabel[i].c_str(), jetname[j].c_str());
			sprintf(htitle, "%s %s (%s)", title.c_str(), cutname[i].c_str(), jetlabel[j].c_str());
			h[i][j] = new TH1F(hname, htitle, nbin, xlow, xup);
		}
	}
}

void ana::valid_fillHisto(TH1F* h[][7], const bool cuts[8], int nj, double value) const {

	for (short i = 0; i < 9; ++i) { //cut
		if (cuts[i]) {
			h[i][6]->Fill(value);//all jets
			if (nj < 5) {
				h[i][nj]->Fill(value);
			}//fill 0-4j
			if (nj > 3) {
				h[i][5]->Fill(value);
			}//fill >=4j
		} else {
			break;
		}
	}
}
//81FB end


//----- my method to add a set of 1D histograms (for each type of MC when running on MC) ------
void ana::addHistoDataAndMC(ushort type, ushort hist, const string name, const string title, const int nbins, const float xmin,
		const float xmax) {
	if (hist >= NUMBER_OF_1D_HISTOGRAMS) {
		cout << "WARNING: 1D Histogram No " << hist << "for MC type " << type << "out of bounds" << endl;
		return;
	} else if (type >= NUMBER_OF_MC_TYPES) {
		cout << "WARNING: MC type " << type << "out of bounds" << endl;
		return;
	}
	if (fasthist_[type][hist] > 0)
		cout << "WARNING: Overwriting existing histogram: Hist '" << fasthist_[type][hist]->GetName() << "'" << endl;

	string htitle = title + "(" + mclabel[type] + ")";
	fasthist_[type][hist] = new TH1F(name.c_str(), htitle.c_str(), nbins, xmin, xmax);
	fasthist_[type][hist]->Sumw2();

}
////----- my method to add a set of 2D histograms (for each type of MC when running on MC) ------
void ana::addHistoDataAndMC(ushort type, ushort hist, const string name, const string title, const int nbinsX, const float xmin,
		const float xmax, const int nbinsY, const float ymin, const float ymax) {
	if (hist >= NUMBER_OF_2D_HISTOGRAMS) {
		cout << "WARNING: 2D Histogram No " << hist << "for MC type " << type << "out of bounds" << endl;
		return;
	} else if (type >= NUMBER_OF_MC_TYPES) {
		cout << "WARNING: MC type " << type << "out of bounds" << endl;
		return;
	}
	if (fasthist2D_[type][hist] > 0)
		cout << "WARNING: Overwriting existing histogram: Hist '" << fasthist2D_[type][hist]->GetName() << "'" << endl;

	string histogram_title = title + "(" + mclabel[type] + ")";
	fasthist2D_[type][hist] = new TH2F(name.c_str(), histogram_title.c_str(), nbinsX, xmin, xmax, nbinsY, ymin, ymax);

}
////-------- my method to fill 1D histograms (acc to MC type when running on MC) ---------
void ana::fillHistoDataAndMC(ushort hist, const float value, const double weight) {
	if (hist >= NUMBER_OF_1D_HISTOGRAMS) {
		cout << "WARNING: Trying to fill 1D histogram out of range: Hist no " << hist + 1 << "out of " << NUMBER_OF_1D_HISTOGRAMS
				<< endl;
		return;
	}
	if (fasthist_[fastmctype_][hist] != NULL && fasthist_[kdata][hist] != NULL) {
		fasthist_[fastmctype_][hist]->Fill(value, weight);

		if (fastmctype_ != kdata)
			fasthist_[kdata][hist]->Fill(value, weight);
		if (fastmctype_ >= kenri1 && fastmctype_ <= kbce3)
			fasthist_[kQCD][hist]->Fill(value, weight);
		if (fastmctype_ >= ktW && fastmctype_ <= kschan)
			fasthist_[ksingleTop][hist]->Fill(value, weight);
	} else
		cout << "WARNING: Trying to fill not initialized 1D histogram: " << hist + 1 << "out of " << NUMBER_OF_1D_HISTOGRAMS
				<< endl;
}
////--------------------------------------------------------------------------------------
////-------- my method to fill 2D histograms (acc to MC type when running on MC) ---------
void ana::fillHistoDataAndMC(ushort hist, const float Xvalue, const float Yvalue, const double weight) {
	if (hist >= NUMBER_OF_2D_HISTOGRAMS) {
		cout << "WARNING: Trying to fill 2D histogram out of range: Hist no(+1) " << hist + 1 << "out of "
				<< NUMBER_OF_2D_HISTOGRAMS << endl;
		return;
	}
	if (fasthist2D_[fastmctype_][hist] != NULL && fasthist2D_[kdata][hist] != NULL) {
		fasthist2D_[fastmctype_][hist]->Fill(Xvalue, Yvalue, weight);
		if (fastmctype_ != kdata)
			fasthist2D_[kdata][hist]->Fill(Xvalue, Yvalue, weight);
		if (fastmctype_ >= kenri1 && fastmctype_ <= kbce3)
			fasthist2D_[kQCD][hist]->Fill(Xvalue, Yvalue, weight);
		if (fastmctype_ >= ktW && fastmctype_ <= kschan)
			fasthist2D_[ksingleTop][hist]->Fill(Xvalue, Yvalue, weight);
	} else
		cout << "WARNING: Trying to fill not initialized 2D histogram: " << hist + 1 << "out of " << NUMBER_OF_2D_HISTOGRAMS
				<< endl;

}
////--------------------------------------------------------------------------------------
//////-------------- my method to fill 2D histograms acc to njet & mctype -----------------------------

////----- add a set of 1D histograms (for each type of MC when running on MC) ------
void ana::addHisto_Njet_DataAndMC(ushort type, ushort hist, const string name, const string title, const int nbin,
		const float xlow, const float xup) {

	const string jetname[NUMBER_OF_JETBINS] = { "0j", "1j", "2j", "3j", "4j", "4mj", "allj" };
	const string jetlabel[NUMBER_OF_JETBINS] = { "0j", "1j", "2j", "3j", "4j", ">=4j", "allj" };
	//	cout << "histname" << name << " number " << hist << endl;
	for (short jbin = 0; jbin < NUMBER_OF_JETBINS; ++jbin) {
		//		cout << jbin << " ";
		string hname = name + "_" + jetname[jbin];
		string htitle = title + " (" + jetlabel[jbin] + ", " + mclabel[type] + ")";
		if (fasthist_jetbinned_[type][jbin][hist] > 0)
			cout << "WARNING: Overwriting existing histogram: Hist '" << fasthist_jetbinned_[type][jbin][hist]->GetName() << "'"
					<< endl;
		fasthist_jetbinned_[type][jbin][hist] = new TH1F(hname.c_str(), htitle.c_str(), nbin, xlow, xup);
		fasthist_jetbinned_[type][jbin][hist]->Sumw2();
	}
}

void ana::iso_fillHisto_nlevel_nj_nmc(const ushort hist, ushort level, const float value, const float weight) {
	if (hist < NUMBER_OF_NLEVEL_JBINNED_1D_HISTOGRAMS && level < NUMBER_OF_LEVELS) {
		//		if (fasthistNlevel1D_[fastmctype_][level][jetbin][hist] != NULL) {
		if (Njet() < 5) {
			fasthistNlevel1D_[fastmctype_][level][Njet()][hist]->Fill(value, weight);

			if (fastmctype_ != kdata)
				fasthistNlevel1D_[kdata][level][Njet()][hist]->Fill(value, weight);
			if (fastmctype_ >= kenri1 && fastmctype_ <= kbce3)
				fasthistNlevel1D_[kQCD][level][Njet()][hist]->Fill(value, weight);
			if (fastmctype_ >= ktW && fastmctype_ <= kschan)
				fasthistNlevel1D_[ksingleTop][level][Njet()][hist]->Fill(value, weight);
		}
		if (Njet() >= 4) {
			fasthistNlevel1D_[fastmctype_][level][FourOrMoreJets][hist]->Fill(value, weight);

			if (fastmctype_ != kdata)
				fasthistNlevel1D_[kdata][level][FourOrMoreJets][hist]->Fill(value, weight);
			if (fastmctype_ >= kenri1 && fastmctype_ <= kbce3)
				fasthistNlevel1D_[kQCD][level][FourOrMoreJets][hist]->Fill(value, weight);
			if (fastmctype_ >= ktW && fastmctype_ <= kschan)
				fasthistNlevel1D_[ksingleTop][level][FourOrMoreJets][hist]->Fill(value, weight);
		}

		fasthistNlevel1D_[fastmctype_][level][AllJets][hist]->Fill(value, weight);
		if (fastmctype_ != kdata)
			fasthistNlevel1D_[kdata][level][AllJets][hist]->Fill(value, weight);
		if (fastmctype_ >= kenri1 && fastmctype_ <= kbce3)
			fasthistNlevel1D_[kQCD][level][AllJets][hist]->Fill(value, weight);
		if (fastmctype_ >= ktW && fastmctype_ <= kschan)
			fasthistNlevel1D_[ksingleTop][level][AllJets][hist]->Fill(value, weight);
		//		} else
		//			cout << "Trying to fill not initiated histogram" << endl;
	}
}
void ana::iso_fillHisto_nlevel_nj_nmc(const ushort hist, ushort level, const float valueX, const float valueY, const float weight) {
	if (hist < NUMBER_OF_NLEVEL_JBINNED_2D_HISTOGRAMS && level < NUMBER_OF_LEVELS)
		//		if (fasthistNlevel2D_[fastmctype_][level][jetbin][hist] != NULL) {
		if (Njet() < 5) {
			fasthistNlevel2D_[fastmctype_][level][Njet()][hist]->Fill(valueX, valueY, weight);

			if (fastmctype_ != kdata)
				fasthistNlevel2D_[kdata][level][Njet()][hist]->Fill(valueX, valueY, weight);
			if (fastmctype_ >= kenri1 && fastmctype_ <= kbce3)
				fasthistNlevel2D_[kQCD][level][Njet()][hist]->Fill(valueX, valueY, weight);
			if (fastmctype_ >= ktW && fastmctype_ <= kschan)
				fasthistNlevel2D_[ksingleTop][level][Njet()][hist]->Fill(valueX, valueY, weight);
		}
	if (Njet() >= 4) {
		fasthistNlevel2D_[fastmctype_][level][FourOrMoreJets][hist]->Fill(valueX, valueY, weight);
		if (fastmctype_ != kdata)
			fasthistNlevel2D_[kdata][level][FourOrMoreJets][hist]->Fill(valueX, valueY, weight);
		if (fastmctype_ >= kenri1 && fastmctype_ <= kbce3)
			fasthistNlevel2D_[kQCD][level][FourOrMoreJets][hist]->Fill(valueX, valueY, weight);
		if (fastmctype_ >= ktW && fastmctype_ <= kschan)
			fasthistNlevel2D_[ksingleTop][level][FourOrMoreJets][hist]->Fill(valueX, valueY, weight);
	}

	fasthistNlevel2D_[fastmctype_][level][AllJets][hist]->Fill(valueX, valueY, weight);
	if (fastmctype_ != kdata)
		fasthistNlevel2D_[kdata][level][AllJets][hist]->Fill(valueX, valueY, weight);
	if (fastmctype_ >= kenri1 && fastmctype_ <= kbce3)
		fasthistNlevel2D_[kQCD][level][AllJets][hist]->Fill(valueX, valueY, weight);
	if (fastmctype_ >= ktW && fastmctype_ <= kschan)
		fasthistNlevel2D_[ksingleTop][level][AllJets][hist]->Fill(valueX, valueY, weight);
	//		} else
	//			cout << "Trying to fill not initiated histogram" << endl;
}
void ana::iso_addHisto_nlevel_nj_nmc(ushort type, ushort hist, const string name, const string title, const int nbins,
		const float xmin, const float xmax) {
	if (m_debug)
		cout << "Starting << iso_addisto_nlevel_nj_nmc >>: " << name << endl;
	const string level_name[NUMBER_OF_LEVELS] = { "L1", "L1b", "L1c", "L1d1", "L1d2", "L1d3", "L1d4", "L1d5", "L2", "L3", "L4" };
	const string level_label[NUMBER_OF_LEVELS] = { "L1 HLT", "L1b E_{T},#eta", "L1c d0", "L1d1 eID H/E",
			"L1d2 eID |#Delta#sigma_{in}|", "L1d3 eID |#Delta#phi_{in}|", "L1d4 eID #sigma_{i#eta i#eta}", "L1d5 eID",
			"L2 GoodEle", "L3 MuZVeto", "L4 Conv" };
	const string jetname[NUMBER_OF_JETBINS] = { "0j", "1j", "2j", "3j", "4j", "4mj", "allj" };
	const string jetlabel[NUMBER_OF_JETBINS] = { "0j", "1j", "2j", "3j", "4j", ">=4j", "allj" };

	for (ushort level = 0; level < NUMBER_OF_LEVELS; ++level) {
		for (ushort jbin = 0; jbin < NUMBER_OF_JETBINS; ++jbin) {
			string hname = name + "_" + level_name[level] + "_" + jetname[jbin];
			string htitle = title + " (" + level_label[level] + "," + jetlabel[jbin] + ", " + mclabel[type] + ")";
			if (fasthistNlevel1D_[type][level][jbin][hist] > 0)
				cout << "WARNING: Overwriting existing histogram: Hist '"
						<< fasthistNlevel1D_[type][level][jbin][hist]->GetName() << "'" << endl;
			fasthistNlevel1D_[type][level][jbin][hist] = new TH1D(hname.c_str(), htitle.c_str(), nbins, xmin, xmax);
			fasthistNlevel1D_[type][level][jbin][hist]->Sumw2();
		}
	}
}

void ana::iso_addHisto_nlevel_nj_nmc(ushort type, ushort hist, const string name, const string title, const int nbinsX,
		const float xmin, const float xmax, const int nbinsY, const float ymin, const float ymax) {
	if (m_debug)
		cout << "Starting << iso_addisto_nlevel_nj_nmc >>: " << name << endl;
	const string level_name[NUMBER_OF_LEVELS] = { "L1", "L1b", "L1c", "L1d1", "L1d2", "L1d3", "L1d4", "L1d5", "L2", "L3", "L4" };
	const string level_label[NUMBER_OF_LEVELS] = { "L1 HLT", "L1b E_{T},#eta", "L1c d0", "L1d1 eID H/E",
			"L1d2 eID |#Delta#sigma_{in}|", "L1d3 eID |#Delta#phi_{in}|", "L1d4 eID #sigma_{i#eta i#eta}", "L1d5 eID",
			"L2 GoodEle", "L3 MuZVeto", "L4 Conv" };
	const string jetname[NUMBER_OF_JETBINS] = { "0j", "1j", "2j", "3j", "4j", "4mj", "allj" };
	const string jetlabel[NUMBER_OF_JETBINS] = { "0j", "1j", "2j", "3j", "4j", ">=4j", "allj" };

	for (ushort level = 0; level < NUMBER_OF_LEVELS; ++level) {
		for (ushort jbin = 0; jbin < NUMBER_OF_JETBINS; ++jbin) {
			string hname = name + "_" + level_name[level] + "_" + jetname[jbin];
			string htitle = title + " (" + level_label[level] + "," + jetlabel[jbin] + ", " + mclabel[type] + ")";

			if (fasthistNlevel2D_[type][level][jbin][hist] > 0)
				cout << "WARNING: Overwriting existing histogram: Hist '"
						<< fasthistNlevel2D_[type][level][jbin][hist]->GetName() << "'" << endl;
			fasthistNlevel2D_[type][level][jbin][hist] = new TH2D(hname.c_str(), htitle.c_str(), nbinsX, xmin, xmax, nbinsY,
					ymin, ymax);
			fasthistNlevel2D_[type][level][jbin][hist]->Sumw2();
		}
	}
}
////-------------- my method to fill  1D  histograms acc to njet & mctype ------------------
void ana::fillHisto_Njet_DataAndMC(ushort hist, const float value, const double w) {

	if (m_debug)
		cout << "\nStart of << fillHisto_Njet_DataAndMC >>: " << fasthist_jetbinned_[fastmctype_][0][hist]->GetName() << endl;
	if (Njet() < 5) {
		fasthist_jetbinned_[fastmctype_][Njet()][hist]->Fill(value, w);

		if (fastmctype_ != kdata)
			fasthist_jetbinned_[kdata][Njet()][hist]->Fill(value, w);
		if (fastmctype_ >= kenri1 && fastmctype_ <= kbce3)
			fasthist_jetbinned_[kQCD][Njet()][hist]->Fill(value, w);
		if (fastmctype_ >= ktW && fastmctype_ <= kschan)
			fasthist_jetbinned_[ksingleTop][Njet()][hist]->Fill(value, w);
	}
	if (Njet() >= 4) {
		fasthist_jetbinned_[fastmctype_][FourOrMoreJets][h_nEle_all]->Fill(value, w);

		if (fastmctype_ != kdata)
			fasthist_jetbinned_[kdata][FourOrMoreJets][h_nEle_all]->Fill(value, w);
		if (fastmctype_ >= kenri1 && fastmctype_ <= kbce3)
			fasthist_jetbinned_[kQCD][FourOrMoreJets][hist]->Fill(value, w);
		if (fastmctype_ >= ktW && fastmctype_ <= kschan)
			fasthist_jetbinned_[ksingleTop][FourOrMoreJets][hist]->Fill(value, w);
	}
	fasthist_jetbinned_[fastmctype_][AllJets][h_nEle_all]->Fill(value, w);

	if (fastmctype_ != kdata)
		fasthist_jetbinned_[kdata][AllJets][h_nEle_all]->Fill(value, w);
	if (fastmctype_ >= kenri1 && fastmctype_ <= kbce3)
		fasthist_jetbinned_[kQCD][AllJets][hist]->Fill(value, w);
	if (fastmctype_ >= ktW && fastmctype_ <= kschan)
		fasthist_jetbinned_[ksingleTop][AllJets][hist]->Fill(value, w);

	if (m_debug)
		cout << "End of << fillHisto_Njet_DataAndMC >>" << endl;
}
//----------------------------------------------------------------------------------------
void ana::fillHisto_PDF_weights(TH1F* h) {
	// 45 bins
	h->Fill(-1, 1);

	h->Fill(0., PDFWcteq66_0);
	h->Fill(1., PDFWcteq66_1);
	h->Fill(2., PDFWcteq66_2);
	h->Fill(3., PDFWcteq66_3);
	h->Fill(4., PDFWcteq66_4);
	h->Fill(5., PDFWcteq66_5);
	h->Fill(6., PDFWcteq66_6);
	h->Fill(7., PDFWcteq66_8);
	h->Fill(8., PDFWcteq66_8);
	h->Fill(9., PDFWcteq66_9);
	h->Fill(10., PDFWcteq66_10);

	h->Fill(11., PDFWcteq66_11);
	h->Fill(12., PDFWcteq66_12);
	h->Fill(13., PDFWcteq66_13);
	h->Fill(14., PDFWcteq66_14);
	h->Fill(15., PDFWcteq66_15);
	h->Fill(16., PDFWcteq66_16);
	h->Fill(17., PDFWcteq66_18);
	h->Fill(18., PDFWcteq66_18);
	h->Fill(19., PDFWcteq66_19);
	h->Fill(20., PDFWcteq66_20);

	h->Fill(21., PDFWcteq66_21);
	h->Fill(22., PDFWcteq66_22);
	h->Fill(23., PDFWcteq66_23);
	h->Fill(24., PDFWcteq66_24);
	h->Fill(25., PDFWcteq66_25);
	h->Fill(26., PDFWcteq66_26);
	h->Fill(27., PDFWcteq66_28);
	h->Fill(28., PDFWcteq66_28);
	h->Fill(29., PDFWcteq66_29);
	h->Fill(30., PDFWcteq66_30);

	h->Fill(31., PDFWcteq66_31);
	h->Fill(32., PDFWcteq66_32);
	h->Fill(33., PDFWcteq66_33);
	h->Fill(34., PDFWcteq66_34);
	h->Fill(35., PDFWcteq66_35);
	h->Fill(36., PDFWcteq66_36);
	h->Fill(37., PDFWcteq66_38);
	h->Fill(38., PDFWcteq66_38);
	h->Fill(39., PDFWcteq66_39);
	h->Fill(40., PDFWcteq66_40);

	h->Fill(41., PDFWcteq66_41);
	h->Fill(42., PDFWcteq66_42);
	h->Fill(43., PDFWcteq66_43);
	h->Fill(44., PDFWcteq66_44);

}

//----------------------------------------------------------------------------------------
// Return code of current event
//   0   data
//   1   ttbar
//   2   QCD
//   3     enri1
//   4     enri2
//   5     enri3
//   6     bce1
//   7     bce2
//   8     bce3
//   9   wj
//  10   zj
//  11   vqq
//  12   singleTop
//  13    tW
//  14    tchan
//  15    schan


//----------------------------------------------------------------------------------------
//                                  Conversion Finder
//----------------------------------------------------------------------------------------
// This was originally ConversionFinder but has been replaced (07-01-10). It is left here until the new routine has been thoroughly tested.
bool ana::ConversionFinder2(const TLorentzVector& e1, int mctype, int index_selected_ele) {

	if (m_debug)
		cout << "Starting << ConversionFinder >>" << endl;
	bool isthisConversion = false;

	const float bfield = 3.8;

	float phi_ie = e1.Phi();
	float eta_ie = e1.PseudoRapidity();

	int mytrackref = static_cast<int> (els_closestCtfTrackRef->at(index_selected_ele));

	//declare track variables
	float phi1, eta1, phi2, eta2;
	float tk1Curvature;
	float tk2Curvature;
	float tk1r;
	float tk2r;
	float d0;
	float d02;
	float tk1x;
	float tk1y;
	float tk2x;
	float tk2y;
	float distmag;
	float dist;
	float dcot;

	//first loop over tracks, trying to match to the electron
	for (unsigned int i = 0; i < Ntracks; ++i) {
		if (mytrackref > -1) {
			i = mytrackref;
		}

		phi1 = tracks_phi->at(i);
		eta1 = tracks_eta->at(i);

		float EleTrackDelR = calcDeltaR(phi1, eta1, phi_ie, eta_ie);

		//consider only tracks that fall within a small cone around electron
		if (EleTrackDelR > 0.3 && mytrackref < 0)
			continue;

		//calculate curvature and radius of track
		tk1Curvature = -0.3 * bfield * tracks_chg->at(i) / (100 * (tracks_pt->at(i)));
		tk1r = fabs(1 / tk1Curvature);

		// NEW (Aug 09)
		d0 = compute_d0("track", i);

		//calculate coordinates of centre of track "circle"
		tk1x = ((1 / tk1Curvature) - d0) * cos(phi1);
		tk1y = ((1 / tk1Curvature) - d0) * sin(phi1);
		//loop over the tracks again, try to find other track in a conversion
		for (unsigned int j = 0; j < Ntracks; ++j) {
			//avoid using the same track:
			if (i == j)
				continue;
			phi2 = tracks_phi->at(j);
			eta2 = tracks_eta->at(j);

			//For a restriction on the position of the second track with respect to the
			//electron, uncomment the line following. The main advantage is a slight
			//reduction in signal loss with very little change in the conversion
			//rejection rate
			//if( calcDeltaR( phi2, eta2, phi_ie, eta_ie) > 0.5 ) continue;

			//check the tracks have opposite charges
			if (tracks_chg->at(i) * tracks_chg->at(j) > 0.0)
				continue;

			tk2Curvature = -0.3 * bfield * tracks_chg->at(j) / (100 * (tracks_pt->at(j)));
			tk2r = fabs(1 / tk2Curvature);

			d02 = compute_d0("track", j);

			tk2x = ((1 / tk2Curvature) - d02) * cos(phi2);
			tk2y = ((1 / tk2Curvature) - d02) * sin(phi2);

			//calculate distance between curves
			distmag = (tk2x - tk1x) * (tk2x - tk1x) + (tk2y - tk1y) * (tk2y - tk1y);
			distmag = sqrt(distmag);
			dist = distmag - (tk1r + tk2r);
			dcot = 1 / tan(tracks_theta->at(i)) - 1 / tan(tracks_theta->at(j));

			// FB/TL (13 Feb 09): update cut values for better performace
			// studies shown that using the new criteria, acceptance of conversion
			// events increase from 33% to ~50%, with littler effect on the signal
			// rejection (1.4% to 2%)
			// if( dist > -0.02 && dist < 0.01 && fabs(dcot) < 0.02)
			if (dist > -0.04 && dist < 0.04 && fabs(dcot) < 0.03) {
				//std::cout <<"Found a Conversion!"<<std::endl;
				isthisConversion = true;
				break;
			}

		}//end of second track loop

		if (isthisConversion)
			break;
		if (mytrackref > -1)
			break;
	}//end of first track loop

	//856
	if (DoConversionStudies() && Nmc_doc == 0) {
		// print only first 20 occurances
		static short mm = 0;
		if (mm == 0)
			cout << "(Printing first 20 occurances.)" << endl;
		if (mm < 20) {
			cout << "PROBLEM: Nmc_doc = 0. Cannot do conversion studies (involve MC matching)." << endl;
			++mm;
		}
	}

	if (DoConversionStudies() && Nmc_doc > 0) {

		int didConv = 0;
		if (isthisConversion) {
			didConv = 1;
		}

		//  int mctype = 1;
		ConversionArray[mctype][didConv][0]++;

		float mc_phi;
		float mc_eta;
		int ii = 0;
		float tempDelR = 100;

		for (unsigned int i = 0; i < Nmc_doc; ++i) {
			mc_phi = mc_doc_phi->at(i);
			mc_eta = mc_doc_eta->at(i);
			float empDelR = calcDeltaR(mc_phi, mc_eta, phi_ie, eta_ie);
			if (empDelR > 0.3)
				continue;

			if (empDelR < tempDelR) {
				tempDelR = empDelR;
				ii = i;
			}
		}//end mc particle loop

		if (fabs(mc_doc_id->at(ii)) == 11) {
			if (m_debug)
				cout << "It matches closest to a real electron" << endl;
			ConversionArray[mctype][didConv][1]++;
		}

		else if (fabs(mc_doc_id->at(ii)) == 22) {
			if (m_debug)
				cout << "It matches closest to a photon" << endl;
			ConversionArray[mctype][didConv][2]++;
		}

		else if (fabs(mc_doc_id->at(ii)) == 111 || fabs(mc_doc_id->at(ii)) == 221) {
			if (m_debug)
				cout << "It matches closest to a pi zero or eta" << endl;
			ConversionArray[mctype][didConv][3]++;
		} else {
			if (m_debug)
				cout << "It matches closest to something else" << endl;
			ConversionArray[mctype][didConv][4]++;
		}

		//end 856
	}

	return isthisConversion;

}

//return generated particle match to reco object based on delR
float ana::MCTruthMatch(float eta, float phi) {
	float mc_phi, mc_eta, tempDelR, finalDelR = 100;
	int ii = 0;
	for (unsigned int i = 0; i < Nmc_doc; ++i) {
		mc_phi = mc_doc_phi->at(i);
		mc_eta = mc_doc_eta->at(i);
		tempDelR = calcDeltaR(mc_phi, mc_eta, phi, eta);
		if (tempDelR > 0.3)
			continue;

		if (tempDelR < finalDelR) {
			finalDelR = tempDelR;
			ii = i;
		}
	}//end mc particle loop
	if (ii == -1)
		cout << "ERROR: index for mc_doc_id is -1" << endl;
	return mc_doc_id->at(ii);
}

/**
 * Match a reconstructed object to a collecton of MC partons
 * using a DeltaR and Pt requirement.
 * Cuts used in CMS AN 2009/095:
 * DeltaR < 0.4
 * relative Pt difference < 3
 *
 * @param recoObj the reconstructed object
 * @param eventpartons the MC partons
 * @return ID of matched parton
 */
short ana::match(TLorentzVector recoObj, vector<TLorentzVector> eventpartons) {
	vector<short> results;
	short id = -1;
	//find closest parton to the reconstructed object
	short closest = findClosest(eventpartons, recoObj);
	//deltaR requirement
	bool dr = recoObj.DeltaR(eventpartons[closest]) < 0.4;
	//Pt requirement
	bool pt = fabs(recoObj.Pt() - eventpartons[closest].Pt()) / eventpartons[closest].Pt() < 3.0;
	if (dr && pt)
		id = closest;
	return id;
}

bool ana::ConversionFinder(const TLorentzVector& e1, int mctype, int index_selected_ele) {

	//Should do no more than return whether the electron passes the conversion algo or not
	//This currently (07-01-10) doesn't actually use mctype or e1. For now leave in case they might be useful for
	//later alterations


	if (m_debug)
		cout << "Starting << ConversionFinder >>" << endl;
	bool isthisConversion = false;

	//if the electron has 3 or more missing layers, declare it a conversion
	if (m_useMisslayers) {
		if (els_innerLayerMissingHits->at(index_selected_ele) > 2)
			return true;
	}
	//else, continue the routine. Note, we may want to apply this separately to the routine. For now, leave it here (07-01-10)


	const float bfield = 3.8;

	int mytrackref = static_cast<int> (els_closestCtfTrackRef->at(index_selected_ele));

	float phi1, eta1, phi2, eta2;
	float tk1Curvature, tk2Curvature;
	float tk1r, tk2r;
	float d0, d02;
	float tk1x, tk1y, tk2x, tk2y;
	float distmag, dist, dcot;
	float EleTrackDelR, theta1, charge1;
	float NumMissLayers = -1;

	//list of cut values subject to possible optimisation:
	// 1) EleTrackDelR > 0.3
	// 2) Dist
	// 3) Dcot
	// 4) Bool opposite track charges - required or not.
	// 5) EleTrackDelR between both tracks < 0.3


	//  mytrackref = 0;
	//if track ref is valid, AND the ctf track and gsf track have a minimum required number of matching hits
	// then use this CTF track (remember to skip in loop over track collection). Otherwise we use the GSF track
	if (mytrackref >= 0 && els_shFracInnerHits->at(index_selected_ele) > 0.45) {//use CTF params
		phi1 = tracks_phi->at(mytrackref);
		eta1 = tracks_eta->at(mytrackref);
		charge1 = tracks_chg->at(mytrackref);
		d0 = compute_d0("track", mytrackref);
		tk1Curvature = -0.3 * bfield * charge1 / (100 * (tracks_pt->at(mytrackref)));
		theta1 = tracks_theta->at(mytrackref);
		if (m_useMisslayers)
			NumMissLayers = tracks_innerLayerMissingHits->at(mytrackref);
	} else {//use GSF params
		phi1 = els_tk_phi->at(index_selected_ele);
		eta1 = els_tk_eta->at(index_selected_ele);
		charge1 = els_tk_charge->at(index_selected_ele);
		d0 = compute_d0("electron", index_selected_ele);
		tk1Curvature = -0.3 * bfield * charge1 / (100 * (els_tk_pt->at(index_selected_ele)));
		theta1 = els_tk_theta->at(index_selected_ele);
		if (m_useMisslayers)
			NumMissLayers = els_innerLayerMissingHits->at(index_selected_ele);
	}

	tk1r = fabs(1 / tk1Curvature);
	tk1x = ((1 / tk1Curvature) - d0) * cos(phi1);
	tk1y = ((1 / tk1Curvature) - d0) * sin(phi1);

	for (unsigned int j = 0; j < Ntracks; ++j) {//loop over the track collection

		if ((abs(mytrackref) - j) == 0)
			continue;//works if trackref is valid or not (invalid it will be less than zero, which j is not)
		phi2 = tracks_phi->at(j);
		eta2 = tracks_eta->at(j);

		//require difference in gsf and ctf charges to be different
		if (charge1 * tracks_chg->at(j) > 0.0)
			continue;

		EleTrackDelR = calcDeltaR(phi1, eta1, phi2, eta2);
		if (EleTrackDelR > 0.3)
			continue;

		tk2Curvature = -0.3 * bfield * tracks_chg->at(j) / (100 * (tracks_pt->at(j)));
		tk2r = fabs(1 / tk2Curvature);

		d02 = compute_d0("track", j);

		tk2x = ((1 / tk2Curvature) - d02) * cos(phi2);
		tk2y = ((1 / tk2Curvature) - d02) * sin(phi2);

		distmag = (tk2x - tk1x) * (tk2x - tk1x) + (tk2y - tk1y) * (tk2y - tk1y);
		distmag = sqrt(distmag);
		dist = distmag - (tk1r + tk2r);
		dcot = 1 / tan(theta1) - 1 / tan(tracks_theta->at(j));

		if (dist > -0.04 && dist < 0.04 && fabs(dcot) < 0.03) {
			bool mlPass = true;
			//require similar number of missing layers for each track
			if (m_useMisslayers) {
				if (fabs(NumMissLayers - tracks_innerLayerMissingHits->at(j)) > 2) {
					mlPass = false;
				}
			}
			if (mlPass) {
				isthisConversion = true;
				break;
			}
		}

		if (isthisConversion)
			break;
	}//end of second track loop

	return isthisConversion;
}

void ana::ConversionMCMatching(const TLorentzVector& e1, int mctype, bool isthisConversion) {

	float phi_ie = e1.Phi();
	float eta_ie = e1.PseudoRapidity();

	int didConv = 0;
	if (isthisConversion) {
		didConv = 1;
	}

	//totals the number of electrons flagged or not flagged as conversions
	ConversionArray[mctype][didConv][0]++;

	float mc_phi, mc_eta, empDelR;
	int ii = 0;
	float tempDelR = 100;

	if (Nmc_doc > 0) {//only run if there are MC particles in the file

		for (unsigned int i = 0; i < Nmc_doc; ++i) {
			mc_phi = mc_doc_phi->at(i);
			mc_eta = mc_doc_eta->at(i);
			empDelR = calcDeltaR(mc_phi, mc_eta, phi_ie, eta_ie);
			if (empDelR > 0.3)
				continue;

			if (empDelR < tempDelR) {
				tempDelR = empDelR;
				ii = i;
			}
		}//end mc particle loop

		if (fabs(mc_doc_id->at(ii)) == 11) {
			if (m_debug)
				cout << "It matches closest to a real electron" << endl;
			ConversionArray[mctype][didConv][1]++;
		}

		else if (fabs(mc_doc_id->at(ii)) == 22) {
			if (m_debug)
				cout << "It matches closest to a photon" << endl;
			ConversionArray[mctype][didConv][2]++;
		}

		else if (fabs(mc_doc_id->at(ii)) == 111 || fabs(mc_doc_id->at(ii)) == 221) {
			if (m_debug)
				cout << "It matches closest to a pi zero or eta" << endl;
			ConversionArray[mctype][didConv][3]++;
		} else if (fabs(mc_doc_id->at(ii)) == 211) {
			if (m_debug)
				cout << "It matches closest to charged pion" << endl;
			ConversionArray[mctype][didConv][4]++;
		} else {
			if (m_debug)
				cout << "It matches closest to something else" << endl;
			ConversionArray[mctype][didConv][5]++;
		}

	}

	return;

}

// 06-01-10: Should remove this or at least change it - the algo is obsolete.
//856
void ana::OptimiseConversionFinder(const TLorentzVector& e1, int mctype) {

	//  bool isthisConversion = false;

	const float bfield = 3.8;

	float phi_ie = e1.Phi();
	float eta_ie = e1.PseudoRapidity();

	float mc_phi;
	float mc_eta;
	int ii = 0;
	float tempDelR = 100;

	for (unsigned int i = 0; i < Nmc_doc; ++i) {
		mc_phi = mc_doc_phi->at(i);
		mc_eta = mc_doc_eta->at(i);
		float empDelR = calcDeltaR(mc_phi, mc_eta, phi_ie, eta_ie);
		if (empDelR > 0.3)
			continue;

		if (empDelR < tempDelR) {
			tempDelR = empDelR;
			ii = i;
		}
	}//end mc particle loop

	//if particle doesn't match a photon and isn't ttbar, return. For ttbar, want all
	if (fabs(mc_doc_id->at(ii)) != 22 && fabs(mc_doc_id->at(ii)) != 111 && mctype > 10) {
		return;
	}

	//declare track variables
	float phi1;
	float eta1;
	float phi2;
	float eta2;
	float tk1Curvature;
	float tk2Curvature;
	float tk1r;
	float tk2r;
	float d0;
	float d02;
	float tk1x;
	float tk1y;
	float tk2x;
	float tk2y;
	float distmag;
	float dist;
	float dcot;

	float Mindist = 100;
	float Mindcot = 100;
	float min_didc = 1000;

	//first loop over tracks, trying to match to the electron
	for (unsigned int i = 0; i < Ntracks; ++i) {

		phi1 = tracks_phi->at(i);
		eta1 = tracks_eta->at(i);

		float EleTrackDelR = calcDeltaR(phi1, eta1, phi_ie, eta_ie);

		//consider only tracks that fall within a small cone around electron
		if (EleTrackDelR > 0.3)
			continue;

		//calculate curvature and radius of track
		tk1Curvature = -0.3 * bfield * tracks_chg->at(i) / (100 * (tracks_pt->at(i)));
		tk1r = fabs(1 / tk1Curvature);

		d0 = compute_d0("track", i);

		//calculate coordinates of centre of track "circle"
		tk1x = ((1 / tk1Curvature) - d0) * cos(phi1);
		tk1y = ((1 / tk1Curvature) - d0) * sin(phi1);

		//loop over the tracks again, try to find other track in a conversion
		for (unsigned int j = 0; j < Ntracks; ++j) {
			//avoid using the same track:
			if (i == j)
				continue;

			phi2 = tracks_phi->at(j);
			eta2 = tracks_eta->at(j);

			//check the tracks have opposite charges
			if (tracks_chg->at(i) * tracks_chg->at(j) > 0.0)
				continue;

			tk2Curvature = -0.3 * bfield * tracks_chg->at(j) / (100 * (tracks_pt->at(j)));
			tk2r = fabs(1 / tk2Curvature);

			d02 = compute_d0("track", j);

			tk2x = ((1 / tk2Curvature) - d02) * cos(phi2);
			tk2y = ((1 / tk2Curvature) - d02) * sin(phi2);

			//calculate distance between curves
			distmag = (tk2x - tk1x) * (tk2x - tk1x) + (tk2y - tk1y) * (tk2y - tk1y);
			distmag = sqrt(distmag);
			dist = distmag - (tk1r + tk2r);
			dcot = 1 / tan(tracks_theta->at(i)) - 1 / tan(tracks_theta->at(j));

			//float temp_mindidc = sqrt(dist*dist/(3.2*3.2) + dcot*dcot/(2.7*2.7));
			float temp_mindidc = sqrt(dist * dist + dcot * dcot);

			if (temp_mindidc < min_didc) {
				min_didc = temp_mindidc;
				Mindist = dist;
				Mindcot = dcot;
			}
		}//end of second track loop
	}//end of first track loop

	//  return isthisConversion;
	if (mctype < 11) {
		Conv_Opti[0]->Fill(Mindist, Mindcot);
		Conv_Optis[0]->Fill(Mindist, Mindcot);
		Conv_Opti_extragran[0]->Fill(Mindist, Mindcot);
	}
	if (mctype == 15) {
		Conv_Opti[1]->Fill(Mindist, Mindcot);
		Conv_Optis[1]->Fill(Mindist, Mindcot);
		Conv_Opti_extragran[1]->Fill(Mindist, Mindcot);
	}

}

void ana::PrintConversionTable() {

	TString MySamples[14] = { "ttbar", "W+jet", "Z+Jet", "Enri1", "Enri2", "Enri3", "bce1", "bce2", "bce3", "vqq", "tW", "tchan",
			"schan", "data" };
	TString ConvNames[13] = { "&  Iso Ele ", "&  Convers ", "& Electron ", "&   Photon ", "&   PiZero ", "&  ChrPion ",
			"&    Other ", "& Non Conv ", "& Electron ", "&   Photon ", "&   PiZero ", "&  ChrPion ", "&    Other " };

	std::cout << std::endl << std::endl << "-------------------------------------------------------------------";
	std::cout << std::endl << "-------------------------------------------------------------------";
	std::cout << std::endl << "----------------Conversion Breakdown-------------------------------";
	std::cout << std::endl << "-------------------------------------------------------------------";
	std::cout << std::endl;
	std::cout
			<< "Shows Conversions with matching particle, and non conversions with matching particle, for all isolated reco electrons"
			<< std::endl << std::endl;

	cout << "mctype ";
	for (int i = 0; i < 12; ++i)
		std::cout << ConvNames[i];
	std::cout << ConvNames[12] << " \\\\" << std::endl;

	int ff;

	for (int k = 10; k != 1; ++k) {
		if (k == 11 || k == 12)
			continue;
		if (k == 10) {
			for (int kkc = 0; kkc < 2; ++kkc) {
				for (int j = 0; j < 6; ++j) {
					for (int kk = 1; kk < 10; ++kk) {
						ConversionArray[10][kkc][j] += ConversionArray[kk][kkc][j];
					}
				}
			}
		}
		if (k == 0) {
			ff = 13;
		} else {
			ff = k - 10;
		}
		std::cout << setw(7) << MySamples[ff] << "&  " << std::setw(7) << (ConversionArray[k][0][0] + ConversionArray[k][1][0]);

		for (int i = 1; i > -1; --i) {
			for (int j = 0; j < 6; ++j) {
				std::cout << " & " << std::setw(8) << ConversionArray[k][i][j];
			}
		}
		std::cout << "  \\\\" << std::endl;
		if (k == 22)
			k = -1;
	}
	cout << endl;

}
// end 856


/**
 * ----------- my method to compute hadronic top mass (m3) --------------
 * m3 = invariant mass of 3 jets which have highest vector sum PT
 * @param jetColl
 * @param nGoodIsoEle
 */
void ana::reco_hadronicTop_highestTopPT(const vector<TLorentzVector>& jetColl, const int nGoodIsoEle) {

	pair<double, double> res = compute_M3(jetColl);
	const double max_top_PT = res.second;
	const double m3 = res.first; //renamed from "had_top_mass"

	// fill histo
	if (nGoodIsoEle > 0) { // this is signal region

		if (jetColl.size() >= 4) {
			h_hadTop_maxPT_pt_4j->Fill(max_top_PT, this_weight);
			h_hadTop_maxPT_mass_4j->Fill(m3, this_weight);
			h_hadTop_maxPT_mass_4j_1000->Fill(m3, this_weight);//15Jan10

			if (IsData() == false) { // fill according to type of MC
				if (isMCType(kttbar)) {
					h_m3_tt->Fill(m3, this_weight);
					h_m3_tt_1000->Fill(m3, this_weight);
				} else if (isMCType(kWjets)) {
					h_m3_wj->Fill(m3, this_weight);
					h_m3_wj_1000->Fill(m3, this_weight);
				} else if (isMCType(kZjets)) {
					h_m3_zj->Fill(m3, this_weight);
					h_m3_zj_1000->Fill(m3, this_weight);
				} else if (isMCType(kQCD)) {
					h_m3_qcd->Fill(m3, this_weight);
					h_m3_qcd_1000->Fill(m3, this_weight);
					if (isMCType(kenri3)) {
						h_m3_enri[2]->Fill(m3, this_weight);
						h_m3_enri_1000[2]->Fill(m3, this_weight);
					} else if (isMCType(kenri2)) {
						h_m3_enri[1]->Fill(m3, this_weight);
						h_m3_enri_1000[1]->Fill(m3, this_weight);
					} else if (isMCType(kenri1)) {
						h_m3_enri[0]->Fill(m3, this_weight);
						h_m3_enri_1000[0]->Fill(m3, this_weight);
					} else if (isMCType(kbce3)) {
						h_m3_bce[2]->Fill(m3, this_weight);
						h_m3_bce_1000[2]->Fill(m3, this_weight);
					} else if (isMCType(kbce2)) {
						h_m3_bce[1]->Fill(m3, this_weight);
						h_m3_bce_1000[1]->Fill(m3, this_weight);
					} else if (isMCType(kbce1)) {
						h_m3_bce[0]->Fill(m3, this_weight);
						h_m3_bce_1000[0]->Fill(m3, this_weight);
					}
				} else if (isMCType(ksingleTop)) {
					h_m3_singletop->Fill(m3, this_weight);
					h_m3_singletop_1000->Fill(m3, this_weight);
				} else if (isMCType(kVQQ)) {
					h_m3_vqq->Fill(m3, this_weight);
					h_m3_vqq_1000->Fill(m3, this_weight);
				}
			}//MC
		}
	}// end pass_iso_ele (signal region)
	else { // fail pass_iso_ele, this is bg region
		if (jetColl.size() >= 4) {
			h_hadTop_maxPT_pt_nonIso_4j->Fill(max_top_PT, this_weight);
			h_hadTop_maxPT_mass_nonIso_4j->Fill(m3, this_weight);
			h_hadTop_maxPT_mass_nonIso_4j_1000->Fill(m3, this_weight);
			if (IsData() == false) { // fill according to type of MC
				if (isMCType(kttbar)) {
					h_m3_tt_control->Fill(m3, this_weight);
					h_m3_tt_control_1000->Fill(m3, this_weight);
				} else if (isMCType(kWjets)) {
					h_m3_wj_control->Fill(m3, this_weight);
					h_m3_wj_control_1000->Fill(m3, this_weight);
				} else if (isMCType(kZjets)) {
					h_m3_zj_control->Fill(m3, this_weight);
					h_m3_zj_control_1000->Fill(m3, this_weight);
				} else if (isMCType(kQCD)) {
					h_m3_qcd_control->Fill(m3, this_weight);
					h_m3_qcd_control_1000->Fill(m3, this_weight);
					if (isMCType(kenri3)) {
						h_m3_enri_control[2]->Fill(m3, this_weight);
						h_m3_enri_control_1000[2]->Fill(m3, this_weight);
					} else if (isMCType(kenri2)) {
						h_m3_enri_control[1]->Fill(m3, this_weight);
						h_m3_enri_control_1000[1]->Fill(m3, this_weight);
					} else if (isMCType(kenri1)) {
						h_m3_enri_control[0]->Fill(m3, this_weight);
						h_m3_enri_control_1000[0]->Fill(m3, this_weight);
					} else if (isMCType(kbce3)) {
						h_m3_bce_control[2]->Fill(m3, this_weight);
						h_m3_bce_control_1000[2]->Fill(m3, this_weight);
					} else if (isMCType(kbce2)) {
						h_m3_bce_control[1]->Fill(m3, this_weight);
						h_m3_bce_control_1000[1]->Fill(m3, this_weight);
					} else if (isMCType(kbce1)) {
						h_m3_bce_control[0]->Fill(m3, this_weight);
						h_m3_bce_control_1000[0]->Fill(m3, this_weight);
					}
				} else if (isMCType(ksingleTop)) {
					h_m3_singletop_control->Fill(m3, this_weight);
					h_m3_singletop_control->Fill(m3, this_weight);
				} else if (isMCType(kVQQ)) {
					h_m3_vqq_control->Fill(m3, this_weight);
					h_m3_vqq_control->Fill(m3, this_weight);
				}
			}
		}
	}

}// end reco_hadronicTop_highestTopPT
//-------------------------------------------------------------------------------------
// NEW
/**
 *
 * @param jetColl
 * @return <M3,3j_vecsum_pt>
 */
pair<double, double> ana::compute_M3(const std::vector<TLorentzVector>& jetColl) const {

	double max_top_PT = 0;
	double had_top_mass = 0; //ie m3
	vector<unsigned int> which_3jet(3); //3 places for int

	// Permute over all 3-jet combination, compute vector sum PT, find 3 jets which give the highest value
	for (unsigned int a = 0; a < (jetColl.size() - 2); ++a) {

		TLorentzVector j1(jetColl[a].Px(), jetColl[a].Py(), jetColl[a].Pz(), jetColl[a].Energy());

		for (unsigned int b = a + 1; b < (jetColl.size() - 1); ++b) {

			TLorentzVector j2(jetColl[b].Px(), jetColl[b].Py(), jetColl[b].Pz(), jetColl[b].Energy());

			for (unsigned int c = b + 1; c < jetColl.size(); ++c) {

				TLorentzVector j3(jetColl[c].Px(), jetColl[c].Py(), jetColl[c].Pz(), jetColl[c].Energy());

				TLorentzVector top(j1 + j2 + j3); //construct a hadronic top candidate
				double this_top_PT = top.Pt();

				// find top quark candidate with the highest PT (ie 3-jet combination with highest PT)
				if (this_top_PT > max_top_PT) {
					max_top_PT = this_top_PT;
					had_top_mass = top.M(); //invariant mass of top
					which_3jet[0] = a;
					which_3jet[1] = b;
					which_3jet[2] = c;
				}
			}//3rd jet
		}//2nd jet
	}//1st jet
	//cout << "BBB selected Top (max pt):  pt=" << max_top_PT << "  m=" << had_top_mass << endl;

	return pair<double, double> (had_top_mass, max_top_PT); //M3,pt
}

// return <M3,3j_vecsum_pt, selected_jets, rejected_jets>
//pair<vector<double> , vector<unsigned int> > ana::compute_M3_modified(const vector<TLorentzVector>& jetColl) const {
//
//	double max_top_PT = 0;
//	double had_top_mass = 0; //ie m3
//	vector<unsigned int> which_3jet(3); //3 places for int
//	vector<double> mass_and_pt(2);
//
//	// Permute over all 3-jet combination, compute vector sum PT, find 3 jets which give the highest value
//	for (unsigned int a = 0; a < (jetColl.size() - 2); ++a) {
//
//		TLorentzVector j1(jetColl[a].Px(), jetColl[a].Py(), jetColl[a].Pz(), jetColl[a].Energy());
//
//		for (unsigned int b = a + 1; b < (jetColl.size() - 1); ++b) {
//
//			TLorentzVector j2(jetColl[b].Px(), jetColl[b].Py(), jetColl[b].Pz(), jetColl[b].Energy());
//
//			for (unsigned int c = b + 1; c < jetColl.size(); ++c) {
//
//				TLorentzVector j3(jetColl[c].Px(), jetColl[c].Py(), jetColl[c].Pz(), jetColl[c].Energy());
//
//				TLorentzVector top(j1 + j2 + j3); //construct a hadronic top candidate
//				double this_top_PT = top.Pt();
//
//				// find top quark candidate with the highest PT (ie 3-jet combination with highest PT)
//				if (this_top_PT > max_top_PT) {
//					max_top_PT = this_top_PT;
//					had_top_mass = top.M(); //invariant mass of top
//					which_3jet[0] = a;
//					which_3jet[1] = b;
//					which_3jet[2] = c;
//				}
//			}//3rd jet
//		}//2nd jet
//	}//1st jet
//
//	mass_and_pt.push_back(had_top_mass);
//	mass_and_pt.push_back(max_top_PT);
//	return pair<vector<double> , vector<unsigned int> > (mass_and_pt, which_3jet); //M3,pt, selected jets, rejected jets
//}

//----------------------------------------------------------------------------------------
//                                 W+jets Estimation
//----------------------------------------------------------------------------------------
// Descriptions:
//   o  The code are designed to run on data _or_ MC.
//   o  When running on data, it takes the m3 template shapes for ttbar and W+jets from MC.
//      The shapes for QCD is taken from data in the control region (nonIsoE), that is 
//      we assume QCD is dominating.
//   o  When running on MC only, it takes the m3 shapes from MC for ttbar, W+jets and QCD.
//      The shapes in the signal region are used to generate the pseudo-data. 
//      The number of event generated is the sum of Poisson-fluctuated n(tt), n(wj) and n(QCD).
//      To do the fit, we use template shape from MC (signal region) for ttbar and W+jets;
//      for QCD, we can choose to use the m3 distribution from control region for: 
//      (a) all events (tt+wj+QCD), or (b) only QCD events. At the moment the code assume
//      approach (b). (To use (a), uncomment lines that says "Use all events for QCD template")
//   o  In both cases, the size of QCD int the fit is constrained to be a Gaussian, with
//      mean equals to the QCD estimate, and width equals to the uncertainty of the estimate.
//   o  The code assumes that the various m3 template histograms from MC are read in from
//      a root file which will be made using a separate macro.
//   o  The fit produces the number of ttbar, W+jets and QCD. It also produces a root file
//      called m3_out.root to store histograms.
//
// IMPORTANT NOTICE: the fit code requires ROOT v5.22.00.
//

void ana::EstimateWjets() {

	EstimateWjets(outputHistFileName, ""); //run on MC
}

bool ana::EstimateWjets(const string inputFile_data, const string inputFile_mc) {

	cout << "--------------------------------------------------------" << endl;
	cout << "                    W+jets Estimation" << endl;
	cout << "--------------------------------------------------------" << endl;

	vector<double> nFitResults;

	//--------------------------
	// import m3 distributions
	//--------------------------
	// aim: to import m3 from one root file

	TFile *fdata(0); //data: real or simulated
	TFile *fmc(0); //template from MC
	fdata = new TFile(inputFile_data.c_str(), "UPDATE");
	if (IsData()) {
		fmc = new TFile(inputFile_mc.c_str()); //separate file
	} else {
		fmc = fdata; //point to same file
	}

	cout << "inputFile_data = " << inputFile_data << endl;
	cout << "inputFile_mc = " << inputFile_mc << endl;

	// used as template to fit, also to generate pseudo-data
	// m3 in signal region (iso)
	string dirW = "Wjets_estimation"; //960 bins
	if (m3_use_1000_bins)
		dirW = "Wjets_estimation_1000"; //1000 bins
	cout << "--------------------------------------------" << endl;
	cout << "Taking m3 histograms from dir: " << dirW << endl;
	cout << "--------------------------------------------" << endl;
	TH1D *h_m3_tt = (TH1D*) fmc->Get(Form("%s/m3_tt", dirW.c_str()));
	TH1D *h_m3_wj = (TH1D*) fmc->Get(Form("%s/m3_wj", dirW.c_str()));

	// 21-5-09: Take wj template from Fastsim normal MC
	const bool take_wj_template_from_fastsim = false;
	if (take_wj_template_from_fastsim == true) {
		cout << "NOTE: Taking W+jets shape from wj Fastsim as the default template." << endl;
		TFile *fwj_fastsim = new TFile("test_mc_wjet_fastsim_nom.root");
		h_m3_wj = (TH1D*) fwj_fastsim->Get(Form("%s/m3_wj", dirW.c_str()));
	}

	TH1D *h_m3_zj(0);
	TH1D *h_m3_stop(0);
	//NB: for QCD template, we use data in control region (was called "control_all")
	TH1D *h_m3_qcd_control(0); //only used when doing "QCD template unc"

	if (IsData() == false) { //MC
		h_m3_qcd = (TH1D*) fmc->Get(Form("%s/m3_qcd", dirW.c_str()));
		h_m3_zj = (TH1D*) fmc->Get(Form("%s/m3_zj", dirW.c_str()));
		h_m3_stop = (TH1D*) fmc->Get(Form("%s/m3_singletop", dirW.c_str()));
		h_m3_qcd_control = (TH1D*) fmc->Get(Form("%s/m3_qcd_control", dirW.c_str()));
	}

	TH1D *h_m3_MCdata_control = (TH1D*) fdata->Get(Form("%s/hadTop_maxPT_mass_nonIso_4j", dirW.c_str()));
	TH1D *control_all = (TH1D*) h_m3_MCdata_control->Clone("control_all");

	//-----------------------------------------------------
	// May 22: adapt code to use Thorsten's M3 histograms
	//-----------------------------------------------------
	const bool use_Thorsten_histo = false;
	const bool use_Thorsten_eventNo = false;

	if (use_Thorsten_histo) {
		cout << "\n *** Using Thorsten's M3 histograms to fit ***" << endl;

		TFile *ftt = new TFile("Thorsten/ttbar.root");
		TFile *fwj = new TFile("Thorsten/wjets_Fastsim.root"); //from FASTsim
		TFile *fqcd = new TFile("Thorsten/QCD_template.root");
		TFile *fstop = new TFile("Thorsten/stop.root");

		h_m3_tt = (TH1D*) ftt->Get("h_m3_100");
		h_m3_wj = (TH1D*) fwj->Get("h_m3_100");
		control_all = (TH1D*) fqcd->Get("h_m3_100");
		h_m3_stop = (TH1D*) fstop->Get("h_m3_100");
	}

	// 1 June 09: to do ttbar ISR/FSR systematics, use fastsim ttbar nominal as reference
	if (sysSample.find("ttbar_fastsim_ISRFSR") != string::npos) { //match "ttbar_fastsim_ISRFSR"
		TFile *ftt = new TFile("test_mc_ttbar_fastsim_ISRFSR_nom.root");
		h_m3_tt = (TH1D*) ftt->Get(Form("%s/m3_tt", dirW.c_str()));
		cout << "\n*** NOTE: For signal, use Fastsim ttbar pythia (nominal ISR/FSR) shape. ***" << endl;
	}

	// 18 June 09: to do ttjet scale/threshold systematics, use fastsim ttjet nominal as reference
	if (sysSample.find("ttjet_fastsim") != string::npos) { //match "ttjet_fastsim"
		TFile *ftt = new TFile("test_mc_ttjet_fastsim_nom.root");
		h_m3_tt = (TH1D*) ftt->Get(Form("%s/m3_tt", dirW.c_str()));
		cout << "\n*** NOTE: To do ttjet scale/threshold systematics, use normal Fastsim ttjet shape as template. ***\n" << endl;
	}

	cout << "m3 distributions:" << endl;
	cout << "pointer tt: " << h_m3_tt << endl;
	cout << "pointer wj: " << h_m3_wj << endl;
	cout << "pointer zj: " << h_m3_zj << endl;
	if (h_m3_tt > 0)
		h_m3_tt->ls();
	if (h_m3_wj > 0)
		h_m3_wj->ls();
	if (h_m3_tt == 0 || h_m3_wj == 0 || control_all == 0) {
		cout << "Warning: zero pointer to either h_m3_tt or h_m3_wj or control_all. Exit EstimateWjets()." << endl;
		return false;
	}
	// If we don't find the m3 histograms, stop execution of code.
	if (IsData() == false) {
		if (h_m3_zj == 0 || h_m3_stop == 0) {
			cout << "Warning: could not find some m3 histogram(s). Exit EstimateWjets()" << endl;
			return false;
		}
	}

	// Define (default) fit templates
	//--------------------------------
	TH1D *temp_tt = (TH1D*) h_m3_tt->Clone("temp_tt");
	TH1D *temp_wj = (TH1D*) h_m3_wj->Clone("temp_wj");
	TH1D *temp_stop = (TH1D*) h_m3_stop->Clone("temp_stop");
	TH1D *temp_qcd = (TH1D*) control_all->Clone("temp_qcd");

	// Define (default) pseudodata models
	//------------------------------------
	TH1D *model_tt = (TH1D*) h_m3_tt->Clone("model_tt");
	TH1D *model_wj = (TH1D*) h_m3_wj->Clone("model_wj");
	TH1D *model_zj = (TH1D*) h_m3_zj->Clone("model_zj");
	TH1D *model_stop = (TH1D*) h_m3_stop->Clone("model_stop");
	TH1D *model_qcd = (TH1D*) control_all->Clone("model_qcd");

	//---------------------------------------------
	// 25 Apr 09: Systematic Shape Uncertainties
	//---------------------------------------------
	// Swith: Use a different shape to evaluate systematic unc

	if (doSystematics == "shape" || doSystematics == "accepShape") {

		cout << "\nDoing *shape* uncertainty\n" << endl;

		TFile *fmc2(0);

		// (1) if vary  ttbar
		if (sysSample.find("tt") != string::npos) { //match "tt"

			if (!use_Thorsten_histo) { //use our histo

				if (sysSample == "ttbar")
					fmc2 = new TFile("test_mc_ttbar.root");
				else if (sysSample == "ttjet_largerISR")
					fmc2 = new TFile("test_mc_ttjet_largerISR.root"); //old
				else if (sysSample == "ttjet_smallerISR")
					fmc2 = new TFile("test_mc_ttjet_smallerISR.root"); //old
				else if (sysSample == "ttbar_fastsim_ISRFSR_large")
					fmc2 = new TFile("test_mc_ttbar_fastsim_ISRFSR_large.root");
				else if (sysSample == "ttbar_fastsim_ISRFSR_small")
					fmc2 = new TFile("test_mc_ttbar_fastsim_ISRFSR_small.root");
				// new sample 15-Jun-09
				else if (sysSample == "ttjet_fastsim_nom")
					fmc2 = new TFile("test_mc_ttjet_fastsim_nom.root");
				else if (sysSample == "ttjet_fastsim_T40")
					fmc2 = new TFile("test_mc_ttjet_fastsim_T40.root");
				else if (sysSample == "ttjet_fastsim_T10")
					fmc2 = new TFile("test_mc_ttjet_fastsim_T10.root");
				else if (sysSample == "ttjet_fastsim_scaleUp")
					fmc2 = new TFile("test_mc_ttjet_fastsim_scaleUp.root");
				else if (sysSample == "ttjet_fastsim_scaleDown")
					fmc2 = new TFile("test_mc_ttjet_fastsim_scaleDown.root");

				if (fmc2 > 0)
					model_tt = (TH1D*) fmc2->Get(Form("%s/m3_tt", dirW.c_str()))->Clone("model_tt");

			} else { //use Thorsten's histo
				cout << "NOTE: Use Thorsten's histo for ttbar shape" << endl;
				if (sysSample == "ttbar")
					fmc2 = new TFile("Thorsten/ttbar_Pythia.root");
				model_tt = (TH1D*) fmc2->Get("h_m3_100")->Clone("model_tt");
			}
		}

		// (2) if vary  wj
		else if (sysSample.find("wj") != string::npos) { //match "wj"

			if (!use_Thorsten_histo) { //use our histo

				if (sysSample == "wjet_fastsim_nom")
					fmc2 = new TFile("test_mc_wjet_fastsim_nom.root");
				else if (sysSample == "wjet_thres20")
					fmc2 = new TFile("test_mc_wjet_thres20.root");
				else if (sysSample == "wjet_thres5")
					fmc2 = new TFile("test_mc_wjet_thres5.root");
				else if (sysSample == "wjet_scaleUp")
					fmc2 = new TFile("test_mc_wjet_scaleUp.root");
				else if (sysSample == "wjet_scaleDown")
					fmc2 = new TFile("test_mc_wjet_scaleDown.root");
				model_wj = (TH1D*) fmc2->Get(Form("%s/m3_wj", dirW.c_str()))->Clone("model_wj");

			} else { //use Thorsten's histo
				cout << "NOTE: Use Thorsten's histo for wj systematics shape" << endl;
				if (sysSample == "wjet_fastsim_nom")
					fmc2 = new TFile("Thorsten/wjets_FastSim.root");
				else if (sysSample == "wjet_thres20")
					fmc2 = new TFile("Thorsten/wjets_FastSim_T20.root");
				else if (sysSample == "wjet_thres5")
					fmc2 = new TFile("Thorsten/wjets_FastSim_T5.root");
				else if (sysSample == "wjet_scaleUp")
					fmc2 = new TFile("Thorsten/wjets_FastSim_scaleup.root");
				else if (sysSample == "wjet_scaleDown")
					fmc2 = new TFile("Thorsten/wjets_FastSim_scaledown.root");
				model_wj = (TH1D*) fmc2->Get("h_m3_100")->Clone("model_wj");
			}
		}
		cout << ",  file: " << fmc2->GetName() << endl;
	}//end if vary shape

	//---------------------------
	// SYS: Unc due to QCD template
	//---------------------------
	if (doSystematics == "QCDtemplate") {
		cout << "For 8.1.2, Uncertainty due to QCD template" << endl;
		cout << "Use only QCD events in control region as data model" << endl;
		model_qcd = (TH1D*) h_m3_qcd_control->Clone("model_qcd");
	}

	// 1-2-09: use more appropriate binning (1000/12.5 bins = 80GeV/bin)
	int rB = 1;
	if (temp_tt->GetNbinsX() == 100)
		rB = 8; //12.5 bins
	else if (temp_tt->GetNbinsX() == 960)
		rB = 80; //12 bins
	else if (temp_tt->GetNbinsX() == 1000)
		rB = 80; //12.5 bins

	temp_tt->Rebin(rB);
	temp_wj->Rebin(rB);
	temp_qcd->Rebin(rB);
	temp_stop->Rebin(rB);

	model_tt->Rebin(rB);
	model_wj->Rebin(rB);
	model_zj->Rebin(rB);
	model_qcd->Rebin(rB);
	model_stop->Rebin(rB);

	if (temp_tt->GetNbinsX() == 75)
		temp_tt->Rebin(6);
	if (model_tt->GetNbinsX() == 75)
		model_tt->Rebin(6);

	//----------------------------------------


	//----------------------------------------
	//  If running on data, nfit is set automatically to 1.
	//  If running on MC, can specify how many toy MC to generate and fit
	//
	int nfit = m_ntoy;
	if (IsData())
		nfit = 1;
	//----------------------------------------

	//-----------------------------------------------------
	// true value: expected signal and bgnd events from MC

	TH2D *SigEvents = (TH2D*) fmc->Get("Signal_njetsVcuts");
	TH2D *WjetsEvents = (TH2D*) fmc->Get("Wjets_njetsVcuts");
	TH2D *ZjetsEvents = (TH2D*) fmc->Get("Zjets_njetsVcuts");
	TH2D *QCDEvents = (TH2D*) fmc->Get("QCD_njetsVcuts");
	TH2D *SingleTopEvents = (TH2D*) fmc->Get("SingleTop_njetsVcuts");

	// new event yields for wj fastsim sys samples (to do accept unc)
	double ntrue_wjet_T20 = 108.9; //125.06;
	double ntrue_wjet_T5 = 82.7; //95.01;
	double ntrue_wjet_scaleUp = 54.6; //62.73;
	double ntrue_wjet_scaleDown = 160.3; //184.09;
	double ntrue_ttbar = 180.7;
	const double ntrue_ttbar_fastsim_ISRFSR_large = 181.95; // unscaled = 159.7;
	const double ntrue_ttbar_fastsim_ISRFSR_small = 169.78; // unscaled = 149.0;
	// ttjet fastsim (scaled to Fullsim)
	const double ntrue_ttjet_fastsim_T40 = 188.4;
	const double ntrue_ttjet_fastsim_T10 = 190.1;
	const double ntrue_ttjet_fastsim_scaleUp = 187.3;
	const double ntrue_ttjet_fastsim_scaleDown = 191.8;

	if (use_Thorsten_eventNo) {
		ntrue_wjet_T20 = 74;
		ntrue_wjet_T5 = 55;
		ntrue_wjet_scaleUp = 36;
		ntrue_wjet_scaleDown = 99;
		ntrue_ttbar = 172;
	}

	double ntt_true = SigEvents->GetBinContent(5, 5); //180.70
	double nwj_true = WjetsEvents->GetBinContent(5, 5);
	double nzj_true = ZjetsEvents->GetBinContent(5, 5);
	double nqcd_true = QCDEvents->GetBinContent(5, 5);
	double nstop_true = SingleTopEvents->GetBinContent(5, 5);

	if (use_Thorsten_eventNo) {
		// Use Thorsten's event yields
		cout << "NOTE: Using Thorsten's event numbers." << endl;
		ntt_true = 177;
		nwj_true = 60;
		nzj_true = 12;
		nqcd_true = 30;
		nstop_true = 8;
	}

	const double nstop_exp = nstop_true; //used in Gaus-constrain


	//----------------------------------
	// QCD normalization uncertainty
	//----------------------------------

	// Default
	double nqcd_est = nqcd_true; // = mean of Gaussian constrain
	double nqcd_est_err = nqcd_true * 0.5; // = width of Gaussian constrain

	// Do as Thorsten, vary the "true" mean in number of QCD evnt in pseudodata

	// Our numbers
	if (doSystematics == "QCDnorm_genest") {
		nqcd_true = 17.0; // QCd est
		nqcd_est = 17.0;
		nqcd_est_err = 8.5;
	} else if (doSystematics == "QCDnorm_genest15") {
		nqcd_true = 25.5; // QCD est +50%
		nqcd_est = 17.0;
		nqcd_est_err = 8.5;
	} else if (doSystematics == "QCDnorm_genest05") {
		nqcd_true = 8.5; // QCD est -50%
		nqcd_est = 17.0;
		nqcd_est_err = 8.5;
	}

	if (use_Thorsten_eventNo) {
		// Thorsten's numbers
		if (doSystematics == "QCDnorm_genest15") {
			nqcd_true = 45; // QCD est +50%
			nqcd_est = 30;
			nqcd_est_err = 15;
		} else if (doSystematics == "QCDnorm_genest05") {
			nqcd_true = 15; // QCD est -50%
			nqcd_est = 30;
			nqcd_est_err = 15;
		}
	}
	cout << "QCD mean of pseudodata (true) = " << nqcd_true << endl;
	cout << "QCD Gaussian-constraint, mean = " << nqcd_est;
	cout << "  width = " << nqcd_est_err << endl;

	//--------------------------------------
	// Single Top normalization uncertainty
	//--------------------------------------
	const double nstop_true_plus = nstop_true * 1.46;
	const double nstop_true_minus = nstop_true * (1 - 0.46);
	if (doSystematics == "stop_plus") { //+46%
		cout << "NOTE: more single top in pseudodata: + 46%, " << nstop_true_plus << endl;
		nstop_true = nstop_true_plus;
	} else if (doSystematics == "stop_minus") { //-46%
		cout << "NOTE: less single top in pseudodata: - 46%, " << nstop_true_minus << endl;
		nstop_true = nstop_true_minus;
	}

	//------------------------------------
	// Acceptance (rate) uncertainty
	//------------------------------------
	// 19-5-09
	if (doSystematics == "accep" || doSystematics == "accepShape") {
		cout << "\n Doing *acceptance* uncertainty\n" << endl;
		// wj
		if (sysSample == "wjet_fastsim_nom") {
			cout << "use the same nwj_true from fullsim" << endl;
			//nwj_true = ntrue_wjet_fastsim_nom;
		} else if (sysSample == "wjet_thres20")
			nwj_true = ntrue_wjet_T20;
		else if (sysSample == "wjet_thres5")
			nwj_true = ntrue_wjet_T5;
		else if (sysSample == "wjet_scaleUp")
			nwj_true = ntrue_wjet_scaleUp;
		else if (sysSample == "wjet_scaleDown")
			nwj_true = ntrue_wjet_scaleDown;
		// tt
		else if (sysSample == "ttbar")
			ntt_true = ntrue_ttbar;
		else if (sysSample == "ttbar_fastsim_ISRFSR_nom") {
			cout << "use the same ntt_true from fullsim" << endl;
			//ntt_true = ntrue_ttbar_fastsim_ISRFSR_nom;
		} else if (sysSample == "ttbar_fastsim_ISRFSR_large")
			ntt_true = ntrue_ttbar_fastsim_ISRFSR_large;
		else if (sysSample == "ttbar_fastsim_ISRFSR_small")
			ntt_true = ntrue_ttbar_fastsim_ISRFSR_small;
		else if (sysSample == "ttjet_fastsim_T40")
			ntt_true = ntrue_ttjet_fastsim_T40;
		else if (sysSample == "ttjet_fastsim_T10")
			ntt_true = ntrue_ttjet_fastsim_T10;
		else if (sysSample == "ttjet_fastsim_scaleUp")
			ntt_true = ntrue_ttjet_fastsim_scaleUp;
		else if (sysSample == "ttjet_fastsim_scaleDown")
			ntt_true = ntrue_ttjet_fastsim_scaleDown;
	}

	cout << "MC true expected number (signal region, isol)" << endl;
	cout << "True n(tt):  " << ntt_true << endl;
	cout << "True n(wj):  " << nwj_true << endl;
	cout << "True n(zj):  " << nzj_true << endl;
	cout << "True n(qcd): " << nqcd_true << endl;
	cout << "True n(singletop): " << nstop_true << endl;

	cout << "QCD Estimate used: n(qcd):  " << nqcd_est << " +/- " << nqcd_est_err << endl;
	//-----------------------------------------------------


	//---------------
	// add  singlet top as Gaussian constraint like qcd
	const double nstop_exp_err = sqrt(nstop_exp + 0.3 * 0.3 * nstop_exp * nstop_exp);
	//---------------

	// Book histograms to store fit results
	TH1D *h_ntt_fit = new TH1D("h_ntt_fit", "n(tt) fit", 100, 0, 500 * intlumi / 20);
	TH1D *h_nwj_fit = new TH1D("h_nwj_fit", "n(wj) fit", 100, -100, 450 * intlumi / 20);
	TH1D *h_nqcd_fit = new TH1D("h_nqcd_fit", "n(qcd) fit", 100, 0 * intlumi / 20, 80 * intlumi / 20); //rescale (20)
	TH1D *h_nstop_fit = new TH1D("h_nstop_fit", "n(singletop) fit", 1000, 0, 20 * intlumi / 20);
	// good fit only
	TH1D *h_ntt_goodfit = new TH1D("h_ntt_goodfit", "n(tt) goodfit", 100, 0, 500 * intlumi / 20);
	TH1D *h_nwj_goodfit = new TH1D("h_nwj_goodfit", "n(wj) goodfit", 100, -100, 450 * intlumi / 20);
	TH1D *h_nqcd_goodfit = new TH1D("h_nqcd_goodfit", "n(qcd) goodfit", 100, 0 * intlumi / 20, 80 * intlumi / 20);
	TH1D *h_nstop_goodfit = new TH1D("h_nstop_goodfit", "n(singletop) goodfit", 1000, 0, 20 * intlumi / 20);

	TH1D *h_ntt_fiterr = new TH1D("h_ntt_fiterr", "error of n(ntt) fit", 1000, 0, 200);
	TH1D *h_nwj_fiterr = new TH1D("h_nwj_fiterr", "error of n(wj) fit", 1000, 0, 200);
	TH1D *h_nqcd_fiterr = new TH1D("h_nqcd_fiterr", "error of n(qcd) fit", 1000, 0, 200);//mod
	TH1D *h_nstop_fiterr = new TH1D("h_nstop_fiterr", "error of n(singletop) fit", 1000, 0, 200);

	// Declare histograms to store number of generated (Poisson-fluctuated) events
	// (when running on MC).
	TH1D *h_ntt_fluc(0);
	TH1D *h_nwj_fluc(0);
	TH1D *h_nzj_fluc(0);
	TH1D *h_nqcd_fluc(0);
	TH1D *h_nstop_fluc(0);
	TH1D *h_ntt_pull(0); //pull of fit results = (fit-true)/error
	TH1D *h_nwj_pull(0);
	TH1D *h_nqcd_pull(0);
	TH1D *h_nstop_pull(0);
	TH1D *h_ntt_goodpull(0); //pull of fit results = (fit-true)/error
	TH1D *h_nwj_goodpull(0);
	TH1D *h_nqcd_goodpull(0);
	TH1D *h_nstop_goodpull(0);

	if (IsData() == false) {
		h_ntt_fluc = new TH1D("h_ntt_fluc", "generated number of tt events (Poisson)", 100, 0, 1000);
		h_nwj_fluc = new TH1D("h_nwj_fluc", "generated number of wj events (Poisson)", 100, 0, 500);
		h_nzj_fluc = new TH1D("h_nzj_fluc", "generated number of zj events (Poisson)", 100, 0, 500);
		h_nqcd_fluc = new TH1D("h_nqcd_fluc", "generated number of qcd events (Poisson)", 100, 0, 200); //rescale (20)
		h_nstop_fluc = new TH1D("h_nstop_fluc", "generated number of singletop events (Poisson)", 100, 0, 100);

		// pull
		h_ntt_pull = new TH1D("h_ntt_pull", "pull of n(tt)", 100, -5, 5);
		h_nwj_pull = new TH1D("h_nwj_pull", "pull of n(wj)", 100, -5, 5);
		h_nqcd_pull = new TH1D("h_nqcd_pull", "pull of n(qcd)", 1200, -7, 5);
		h_nstop_pull = new TH1D("h_nstop_pull", "pull of n(singletop)", 1000, -5, 5);

		h_ntt_goodpull = new TH1D("h_ntt_goodpull", "pull of n(tt) good", 100, -5, 5);
		h_nwj_goodpull = new TH1D("h_nwj_goodpull", "pull of n(wj) good", 100, -5, 5);
		h_nqcd_goodpull = new TH1D("h_nqcd_goodpull", "pull of n(qcd) good", 1200, -7, 5);
		h_nstop_goodpull = new TH1D("h_nstop_goodpull", "pull of n(singletop) good", 1000, -5, 5);

	}
	// keep a few set of pseudo-data
	TH1D *h_m3_sum_toy_keep[5];

	for (short i = 0; i < 5; ++i) {
		h_m3_sum_toy_keep[i] = new TH1D(Form("h_m3_sum_toy_%d", i + 1), Form("m3 sum toy mc %d", i + 1), 100, 0, 1000);
	}

	// keep data (real/sim) to write to root file
	TH1D *h_m3_data_keep(0);
	TH1D *h_m3_data_control_keep(0);

	// keep a few qcd pdf (use data in control region)
	RooPlot* m3_qcd_pdf[5];

	bool printText = true;

	RooPlot* m3frame1(0);
	RooPlot* m3frame2(0);
	RooPlot* m3frame3(0);
	RooPlot* m3frame4(0);
	RooPlot* m3frame5(0);
	RooPlot* m3frame6(0);
	RooPlot* m3frame7(0);
	RooPlot* m3frame8(0);
	RooPlot* m3frame9(0);
	RooPlot* m3frame10(0);
	RooPlot* m3frame11(0);

	//---------------
	// perform fit
	//---------------


	for (int n = 0; n < nfit; ++n) {

		if (n > 20)
			printText = false;
		cout << "Fits carried out: " << n + 1 << " of " << nfit << endl;
		if ((n + 1) % 100 == 1 || (n + 1) == nfit) {
			printText = true;
		}

		//--------------------------------------
		// 1) Get the m3 distribution
		//    (can be real data or pseudo-data)
		//--------------------------------------
		TH1D *h_m3_data(0); // in signal region (Iso)
		TH1D *h_m3_data_control(0); // in control region (nonIso) - take this to be QCD shape


		int ntt_fluc = 0;
		int nwj_fluc = 0;
		int nzj_fluc = 0;
		int nqcd_fluc = 0;
		int nstop_fluc = 0;

		if (IsData() == true) { // Real LHC data

			TH1D *h_m3_LHCdata = (TH1D*) fdata->Get(Form("%s/hadTop_maxPT_mass_4j", dirW.c_str()));
			TH1D *h_m3_LHCdata_control = (TH1D*) fdata->Get(Form("%s/hadTop_maxPT_mass_nonIso_4j", dirW.c_str()));
			h_m3_data = (TH1D*) h_m3_LHCdata->Clone("h_m3_data"); //rename
			h_m3_data_control = (TH1D*) h_m3_LHCdata_control->Clone("h_m3_data_control");

		} else { // Monte Carlo

			// a) Poission fluctuated number of events
			ntt_fluc = gRandom->Poisson(ntt_true);
			nwj_fluc = gRandom->Poisson(nwj_true);
			nzj_fluc = gRandom->Poisson(nzj_true);
			nqcd_fluc = gRandom->Poisson(nqcd_true);
			nstop_fluc = gRandom->Poisson(nstop_true);
			h_ntt_fluc ->Fill(ntt_fluc);
			h_nwj_fluc ->Fill(nwj_fluc);
			h_nzj_fluc ->Fill(nzj_fluc);
			h_nqcd_fluc ->Fill(nqcd_fluc);
			h_nstop_fluc->Fill(nstop_fluc);

			if (printText) {
				cout << "\n>>> toy experiment #" << n + 1 << " of " << nfit << endl;
				cout << "-------------------------------------" << endl;
				cout << "Poisson-fluc n(tt):  " << ntt_fluc << endl;
				cout << "Poisson-fluc n(wj):  " << nwj_fluc << endl;
				cout << "Poisson-fluc n(zj):  " << nzj_fluc << endl;
				cout << "Poisson-fluc n(qcd): " << nqcd_fluc << endl;
				cout << "Poisson-fluc n(singletop): " << nstop_fluc << endl;
			}

			// b) Book histos to store pseudo-data (can use smaller binning if needed)
			// TEST 15-Jan-10
			const int nbin = temp_tt->GetNbinsX();
			const float max = temp_tt->GetXaxis()->GetXmax();
			if (printText)
				cout << "\n m3 Pseudo-data: " << nbin << " bins, range 0-" << max << endl << endl;
			// signal region
			TH1D *h_m3_tt_toy = new TH1D("h_m3_tt_toy", "h_m3_tt_toy", nbin, 0, max);
			TH1D *h_m3_wj_toy = new TH1D("h_m3_wj_toy", "h_m3_wj_toy", nbin, 0, max);
			TH1D *h_m3_zj_toy = new TH1D("h_m3_zj_toy", "h_m3_zj_toy", nbin, 0, max);
			TH1D *h_m3_qcd_toy = new TH1D("h_m3_qcd_toy", "h_m3_qcd_toy", nbin, 0, max);
			TH1D *h_m3_stop_toy = new TH1D("h_m3_stop_toy", "h_m3_stop_toy", nbin, 0, max);
			TH1D *h_m3_sum_toy = new TH1D("h_m3_sum_toy", "h_m3_sum_toy", nbin, 0, max);

			//-------------------------------------------------
			// c) Generate pseudo-data for this toy experiment
			//-------------------------------------------------
			for (int i = 0; i < ntt_fluc; ++i) {
				h_m3_tt_toy->Fill(model_tt->GetRandom());
			}
			for (int i = 0; i < nwj_fluc; ++i) {
				h_m3_wj_toy->Fill(model_wj->GetRandom());
			}
			for (int i = 0; i < nzj_fluc; ++i) {
				h_m3_zj_toy->Fill(model_zj->GetRandom());
			}
			for (int i = 0; i < nstop_fluc; ++i) {
				h_m3_stop_toy->Fill(model_stop->GetRandom());
			}
			for (int i = 0; i < nqcd_fluc; ++i) {
				h_m3_qcd_toy->Fill(model_qcd->GetRandom());
			}

			// d) prepares the pseudo-data in signal region
			h_m3_sum_toy->Add(h_m3_tt_toy);
			h_m3_sum_toy->Add(h_m3_wj_toy);
			if (doSystematics == "addZj_fitWithWjShape") {
				h_m3_sum_toy->Add(h_m3_zj_toy);
			}
			h_m3_sum_toy->Add(h_m3_stop_toy);
			h_m3_sum_toy->Add(h_m3_qcd_toy);

			if (n < 5) {
				h_m3_sum_toy_keep[n] = (TH1D*) h_m3_sum_toy->Clone(Form("h_m3_sum_toy_%d", n + 1));
				h_m3_sum_toy_keep[n]->SetTitle(Form("m3 sum toy mc %d", n + 1));
			}

			const double n_gen_event = h_m3_sum_toy->GetEntries();
			if (printText)
				cout << "Total gen event (Poisson-fluc): " << n_gen_event << endl;

			// prepare fake data distributions
			h_m3_data = (TH1D*) h_m3_sum_toy->Clone("h_m3_data");
			h_m3_data_control = (TH1D*) control_all->Clone("h_m3_data_control");

			delete h_m3_tt_toy;
			delete h_m3_wj_toy;
			delete h_m3_zj_toy;
			delete h_m3_qcd_toy;
			delete h_m3_stop_toy;
			delete h_m3_sum_toy;

		}// end MC

		if (n == 0) {
			h_m3_data_keep = (TH1D*) h_m3_data->Clone("h_m3_data_keep");
			h_m3_data_control_keep = (TH1D*) h_m3_data_control->Clone("h_m3_data_control_keep");
		}

		// number of event observed in data (either real or fake)
		const double n_event_obs = h_m3_data->GetEntries();

		// export generated m3 histogram (pseudo-data) to RooDataHist object
		RooRealVar m3("m3", "m3", 0, 960); //changed from 1000 to 960
		RooDataHist data("data", "dataset with m3", m3, h_m3_data); //real or MC data


		//------------------------
		//  2)      f i t
		//------------------------
		// signal and background histograms to make PDFs
		RooDataHist rh_tt("rh_tt", "tt", m3, temp_tt);
		RooDataHist rh_wj("rh_wj", "wj", m3, temp_wj);
		RooDataHist rh_stop("rh_stop", "singletop", m3, temp_stop);
		RooDataHist rh_qcd("rh_qcd", "qcd", m3, temp_qcd);

		// use 0-th order interpolation
		RooHistPdf pdf_tt("pdf_tt", "Signal pdf", m3, rh_tt, 0);
		RooHistPdf pdf_wj("pdf_wj", "W+jets pdf", m3, rh_wj, 0);
		RooHistPdf pdf_qcd("pdf_qcd", "QCD pdf ", m3, rh_qcd, 0);
		RooHistPdf pdf_singletop("pdf_singletop", "single top pdf", m3, rh_stop, 0);

		// create fit model and parameters                       min    max   unit
		// limit on each parameter is 0 -- total number of observed events
		RooRealVar ntt("ntt", "number of tt signal events", -150, 2 * n_event_obs, "event");
		RooRealVar nwj("nwj", "number of W+jets bgnd events", -150, 2 * n_event_obs, "event");

		//RooRealVar ntt ("ntt",  "number of tt signal events",   -150, n_event_obs, "event");
		//RooRealVar nwj ("nwj",  "number of W+jets bgnd events", -150, n_event_obs, "event");

		RooRealVar nqcd("nqcd", "number of QCD bgnd events", nqcd_est, 0, n_event_obs, "event");
		RooRealVar nstop("nstop", "number of single top bgnd events", nstop_exp, 0, n_event_obs, "event");

		// add constrain to QCD PDF
		RooConstVar qcd_est("qcd_est", "qcd est", (nqcd_est)); //mean of gaus
		RooConstVar qcd_est_unc("qcd_est_unc", "uncertainty of qcd est", nqcd_est_err); //sigma of gaus
		RooGaussian nqcd_constraint("nqcd_constraint", "nqcd_constraint", nqcd, qcd_est, qcd_est_unc);

		// multiply constraint with QCD PDF
		RooProdPdf pdf_qcd_constraint("pdf_qcd_constraint", "constrained QCD pdf", RooArgSet(pdf_qcd, nqcd_constraint));

		// Mar 09: add constraint to single top PDF using expected value from MC
		RooConstVar singletop_exp("singletop_exp", "singletop exp", (nstop_exp)); //mean of gauss
		RooConstVar singletop_exp_unc("singletop_exp_unc", "uncertainty of singletop exp", nstop_exp_err); //sigma of gaus
		RooGaussian nstop_constraint("nstop_constraint", "nstop_constraint", nstop, singletop_exp, singletop_exp_unc);
		// Multiply constraint with single top PDF
		RooProdPdf pdf_singletop_constraint("pdf_singletop_constraint", "constrained single top pdf", RooArgSet(pdf_singletop,
				nstop_constraint));

		// Create a model with constrained QCD and single top
		//----------------
		if (printText)
			cout << " ---> Create data model" << endl;
		RooAddPdf model2("model2", "sig+wj+qcd+stop", RooArgList(pdf_tt, pdf_wj, pdf_qcd_constraint, pdf_singletop_constraint),
				RooArgList(ntt, nwj, nqcd, nstop));

		// Fit model to data  (with internal constraint on nqcd and nstop)
		//--------------------
		if (printText)
			cout << " ---> Perform RooFit" << endl;
		model2.fitTo(data, Extended(kTRUE), Constrain(RooArgList(nqcd, nstop)), PrintLevel(-1));

		if (n == 0) {
			m3frame1 = m3.frame();
			m3frame2 = m3.frame();
			m3frame3 = m3.frame();
			m3frame4 = m3.frame();
			m3frame5 = m3.frame();
			m3frame6 = m3.frame();
			m3frame7 = m3.frame();
			m3frame8 = m3.frame();
			m3frame9 = m3.frame();
			m3frame10 = m3.frame();
			m3frame11 = m3.frame();
			//m3frame21 = m3.frame();
			//m3frame22 = m3.frame();


			//-----------------
			// 15-12-09
			//-----------------
			// error type (Def=Poission for integer histogram bin content)
			// ref: http://root.cern.ch/root/html/RooAbsData.html
			//  RooAbsData::Poission
			//  RooAbsData::SumW2
			//  RooAbsData::None
			//  RooAbsData::Auto
			RooAbsData::ErrorType etype = RooAbsData::Poisson; //Def
			if (!IsData())
				etype = RooAbsData::SumW2;
			//-----------------


			rh_tt.plotOn(m3frame1, MarkerSize(1), DataError(etype));
			rh_wj.plotOn(m3frame2, MarkerSize(1), DataError(etype));
			rh_qcd.plotOn(m3frame3, MarkerSize(1), DataError(etype));

			pdf_tt.plotOn(m3frame4);
			pdf_wj.plotOn(m3frame5);
			pdf_qcd.plotOn(m3frame6);
			pdf_qcd_constraint.plotOn(m3frame7);

			data.plotOn(m3frame8, MarkerSize(1)); //plot pseudo-data
			model2.plotOn(m3frame8); //plot composite pdf (s+b model)
			model2.plotOn(m3frame8, Components(pdf_tt), LineStyle(kDashed), LineColor(kBlue + 1));
			model2.plotOn(m3frame8, Components(pdf_wj), LineStyle(kDashed), LineColor(kRed + 1));
			model2.plotOn(m3frame8, Components(pdf_qcd_constraint), LineStyle(kDashed), LineColor(kOrange - 6));

			// single top plots
			rh_stop.plotOn(m3frame9, MarkerSize(1), DataError(etype));
			pdf_singletop.plotOn(m3frame10);
			pdf_singletop_constraint.plotOn(m3frame11);
		}

		// QCD pdf
		if (n < 5) {
			m3_qcd_pdf[n] = m3.frame();
			m3_qcd_pdf[n]->SetNameTitle(Form("m3_qcd_pdf_%d", n + 1), Form("m3 pdf: qcd (control region) toy %d", n + 1));
			pdf_qcd.plotOn(m3_qcd_pdf[n]);
		}

		if (printText) {
			cout << "\n" << n + 1 << ") Fit results (with constraint):\n";
			ntt.Print();
			nwj.Print();
			nqcd.Print();
			nstop.Print();
		}
		// record fit results
		double ntt_fit = ntt.getVal();
		double nwj_fit = nwj.getVal();
		double nqcd_fit = nqcd.getVal();
		double nstop_fit = nstop.getVal();

		double ntt_fiterr = ntt.getError();
		double nwj_fiterr = nwj.getError();
		double nqcd_fiterr = nqcd.getError();
		double nstop_fiterr = nstop.getError();

		h_ntt_fiterr ->Fill(ntt_fiterr);
		h_nwj_fiterr ->Fill(nwj_fiterr);
		h_nqcd_fiterr->Fill(nqcd_fiterr);
		h_nstop_fiterr->Fill(nstop_fiterr);

		h_ntt_fit ->Fill(ntt_fit);
		h_nwj_fit ->Fill(nwj_fit);
		h_nqcd_fit->Fill(nqcd_fit);
		h_nstop_fit->Fill(nstop_fit);

		// calculate the pull when running on MC
		if (IsData() == false) {
			double ntt_pull = (ntt_fit - ntt_true) / ntt_fiterr;
			double nwj_pull = (nwj_fit - nwj_true) / nwj_fiterr;
			if (doSystematics == "addZj_fitWithWjShape") {
				nwj_pull = (nwj_fit - nwj_true - nzj_true) / nwj_fiterr; //include also z+jets
			}
			double nqcd_pull = (nqcd_fit - nqcd_true) / nqcd_fiterr;
			double nstop_pull = (nstop_fit - nstop_true) / nstop_fiterr;

			h_ntt_pull ->Fill(ntt_pull);
			h_nwj_pull ->Fill(nwj_pull);
			h_nqcd_pull->Fill(nqcd_pull);
			h_nstop_pull->Fill(nstop_pull);

			// to exclude problematic fit
			if (ntt_fiterr < 100 && nwj_fiterr > 20) {
				h_ntt_goodfit ->Fill(ntt_fit);
				h_nwj_goodfit ->Fill(nwj_fit);
				h_nqcd_goodfit->Fill(nqcd_fit);
				h_nstop_goodfit->Fill(nstop_fit);
				h_ntt_goodpull ->Fill(ntt_pull);
				h_nwj_goodpull ->Fill(nwj_pull);
				h_nqcd_goodpull->Fill(nqcd_pull);
				h_nstop_goodpull->Fill(nstop_pull);
			}
			if (printText) {
				cout << "fit ntt:  " << ntt_fit << " pm " << ntt_fiterr;
				cout << " (pull = " << ntt_pull << ")" << endl;
				cout << "fit nwj:  " << nwj_fit << " pm " << nwj_fiterr;
				cout << " (pull = " << nwj_pull << ")" << endl;
				cout << "fit nqcd: " << nqcd_fit << " pm " << nqcd_fiterr;
				cout << " (pull = " << nqcd_pull << ")" << endl;
				cout << "fit nstop: " << nstop_fit << " pm " << nstop_fiterr;
				cout << " (pull = " << nstop_pull << ")" << endl;
				cout << endl;

				//record results for the 1st toy MC only
				if (n == 0) {
					nFitResults.push_back(ntt_fit);
					nFitResults.push_back(ntt_fiterr);
					nFitResults.push_back(nwj_fit);
					nFitResults.push_back(nwj_fiterr);
					nFitResults.push_back(nqcd_fit);
					nFitResults.push_back(nqcd_fiterr);
					nFitResults.push_back(nstop_fit);
					nFitResults.push_back(nstop_fiterr);
				}
				// 1 Feb 09: Print result table
				printf("\n               %8s  %6s %18s  %12s\n", "true", "fluc", "Gaus constraint", "fit");
				printf(" ttbar           %8.2f %6d  %18s  %8.2f +/- %4.2f\n", ntt_true, ntt_fluc, "-", ntt_fit, ntt_fiterr);
				// W only
				printf(" W+jets          %8.2f %6d  %18s  %8.2f +/- %4.2f\n", nwj_true, nwj_fluc, "-", nwj_fit, nwj_fiterr);
				// W+Z
				if (doSystematics == "addZj_fitWithWjShape") {
					printf(" W/Z+jets      %8.2f %6d  %18s  %8.2f +/- %4.2f\n", nwj_true + nzj_true, nwj_fluc + nzj_fluc, "-",
							nwj_fit, nwj_fiterr);
				}
				printf(" QCD             %8.2f %6d  %8.2f +/- %4.2f  %8.2f +/- %4.2f\n", nqcd_true, nqcd_fluc, nqcd_est,
						nqcd_est_err, nqcd_fit, nqcd_fiterr);
				printf(" single top      %8.2f %6d  %8.2f +/- %4.2f  %8.2f +/- %4.2f\n", nstop_true, nstop_fluc, nstop_exp,
						nstop_exp_err, nstop_fit, nstop_fiterr);

				printf("\n");
				printf(" W+jets          %8.2f %6d\n", nwj_true, nwj_fluc);
				printf(" Z+jets          %8.2f %6d\n", nzj_true, nzj_fluc);
				printf(" single top      %8.2f %6d\n", nstop_true, nstop_fluc);
			}
		}

		delete h_m3_data;
		delete h_m3_data_control;

	}// end nfit loop

	//add M3 method data to root file

	TH1D *Signal_fitted = new TH1D("Signal_fitted", "Signal fitted", 2, 0, 2);
	TH1D *Wjet_fitted = new TH1D("Wjet_fitted", "Wjet fitted", 2, 0, 2);
	TH1D *QCD_fitted = new TH1D("QCD_fitted", "QCD fitted", 2, 0, 2);
	TH1D *SingleTop_fitted = new TH1D("SingleTop_fitted", "Single Top fitted", 2, 0, 2);

	Signal_fitted->SetBinContent(1, nFitResults.at(0));
	Signal_fitted->SetBinContent(2, nFitResults.at(1));

	Wjet_fitted->SetBinContent(1, nFitResults.at(2));
	Wjet_fitted->SetBinContent(2, nFitResults.at(3));

	QCD_fitted->SetBinContent(1, nFitResults.at(4));
	QCD_fitted->SetBinContent(2, nFitResults.at(5));

	SingleTop_fitted->SetBinContent(1, nFitResults.at(6));
	SingleTop_fitted->SetBinContent(2, nFitResults.at(7));

	// remove pointer to histograms in input root file
	delete SigEvents;
	delete QCDEvents;
	delete WjetsEvents;
	delete ZjetsEvents;
	delete SingleTopEvents;

	delete h_m3_tt; //template
	delete h_m3_wj; //template
	delete h_m3_MCdata_control;
	delete control_all;
	if (IsData() == false) {
		delete h_m3_qcd;
		delete h_m3_zj;
		delete h_m3_stop;
	}

	if (fdata->GetDirectory("m3_fit") > 0) {
		fdata->rmdir("m3_fit");
	}

	TDirectory *m3_fit = fdata->mkdir("m3_fit");
	m3_fit->cd();

	// write histo to root file
	h_m3_data_keep->Write();
	h_m3_data_control_keep->SetTitle("m3 (data in control reg, used as QCD template)");
	h_m3_data_control_keep->Write();

	m3frame1->SetName("m3frame1");
	m3frame1->SetTitle("m3 hist: tt (mc)");
	m3frame1->Write(); //rh_tt
	m3frame2->SetName("m3frame2");
	m3frame2->SetTitle("m3 hist: wj (mc)");
	m3frame2->Write(); //rh_wj
	m3frame3->SetName("m3frame3");
	m3frame3->SetTitle("m3 hist: qcd (data control)");
	m3frame3->Write(); //rh_data_control
	m3frame4->SetName("m3frame4");
	m3frame4->SetTitle("m3 pdf: tt (mc)");
	m3frame4->Write(); //pdf_tt
	m3frame5->SetName("m3frame5");
	m3frame5->SetTitle("m3 pdf: wj (mc)");
	m3frame5->Write(); //pdf_wj
	m3frame6->SetName("m3frame6");
	m3frame6->SetTitle("m3 pdf: qcd (control)");
	m3frame6->Write();//pdf_qcd_control
	m3frame7->SetName("m3frame7");
	m3frame7->SetTitle("m3 pdf: qcd constraint");
	m3frame7->Write(); //pdf_qcd_constraint
	m3frame8->SetName("m3frame8");
	m3frame8->SetTitle("m3 pseudo-data and fit");
	m3frame8->Write();
	m3frame9->SetName("m3frame9");
	m3frame9->SetTitle("m3 hist: single top (mc)");
	m3frame9->Write();
	m3frame10->SetName("m3frame10");
	m3frame10->SetTitle("m3 pdf: single top (mc)");
	m3frame10->Write();
	m3frame11->SetName("m3frame11");
	m3frame11->SetTitle("m3 pdf: single top constraint (mc)");
	m3frame11->Write();

	h_ntt_fit->Write();
	h_nwj_fit->Write();
	h_nqcd_fit->Write();
	h_nstop_fit->Write();

	h_ntt_goodfit->Write();
	h_nwj_goodfit->Write();
	h_nqcd_goodfit->Write();
	h_nstop_goodfit->Write();

	h_ntt_fiterr->Write();
	h_nwj_fiterr->Write();
	h_nqcd_fiterr->Write();
	h_nstop_fiterr->Write();

	h_ntt_pull->Write();
	h_nwj_pull->Write();
	h_nqcd_pull->Write();
	h_nstop_pull->Write();

	h_ntt_goodpull->Write();
	h_nwj_goodpull->Write();
	h_nqcd_goodpull->Write();
	h_nstop_goodpull->Write();

	h_ntt_fluc->Write();
	h_nwj_fluc->Write();
	h_nzj_fluc->Write();
	h_nqcd_fluc->Write();
	h_nstop_fluc->Write();

	const int nn = min(nfit, 5);
	for (short i = 0; i < nn; ++i)
		h_m3_sum_toy_keep[i]->Write();
	for (short i = 0; i < nn; ++i)
		m3_qcd_pdf[i]->Write();

	Signal_fitted->Write();
	Wjet_fitted->Write();
	QCD_fitted->Write();
	SingleTop_fitted->Write();

	model_tt->Write();
	model_wj->Write();
	model_zj->Write();
	model_qcd->Write();
	model_stop->Write();

	temp_tt->Write();
	temp_wj->Write();
	temp_qcd->Write();
	temp_stop->Write();

	fdata->Close();
	fmc->Close();

	cout << "qcd est: " << nqcd_est << ", unc: " << nqcd_est_err << endl;
	cout << "--------------------------------------------------------" << endl;
	cout << "                    W+jets Estimation (end)" << endl;
	cout << "--------------------------------------------------------\n" << endl;

	return true;
}
//----------------------------------------------------------------------------------------
//                                 W+jets Estimation (end)
//----------------------------------------------------------------------------------------

// Add 19-5-09
void ana::StudySystematics(const string name, const string name2) {
	doSystematics = name;
	sysSample = name2;
	cout << "Study Systematics" << endl;
	cout << "doSystematics = " << doSystematics << endl;
	cout << "sysSample = " << sysSample << endl;
}

//----------------------------------------------------------------------------------------
void ana::DrawEventPerNjetTable(const double nevent[][5][nmctype], const vector<string>& ve) const {

	cout << setw(23) << " &" << setw(13) << "0-jet" << " &" << setw(13) << "1-jet" << " &" << setw(13) << "2-jets" << " &"
			<< setw(13) << "3-jets" << " &" << setw(13) << "$\\ge$4-jets" << " &" << setw(23) << "Total \\\\\n\\hline";
	cout << endl;

	for (short i = 0; i < 11; ++i) { //nstage
		printCutStage(i, ve.at(i));
		double total = 0;
		for (short j = 0; j < ntjet; ++j) {
			cout << " & " << setw(12) << fixed << nevent[i][j][0];
			total += nevent[i][j][0];
		}
		cout << " & " << setw(13) << fixed << total << " \\\\" << endl;
	}
	cout << "\n\\hline" << endl;

}//end DrawEventPerNjetTable


//----------------------------------------------------------------------------------------
void ana::DrawMCTypeTable(const double nevent[14][5][nmctype], const string title, vector<string> ve) const {
	vector<string> ve_initial = ve;
	static bool first_time = true;
	if (first_time)
		cout.precision(0); //unweighted table
	else
		cout.precision(myprec); //weighted table

	if (first_time)
		cout << "\\newpage\n" << endl;

	cout << "\n%------------------------------------------------------------------------\n";
	cout << "%       " << title;
	cout << "\n%------------------------------------------------------------------------\n";

	if (first_time) {
		cout << "\\begin{tabular}{|l|rrrrrrr|r|}" << endl;
		cout << "\\hline" << endl;
	}
	cout << "\\multicolumn{9}{|l|}";

	if (first_time)
		cout << "{Actual number of MC events passing selection}";
	else
		cout << "{Expected number of events for " << intlumi << "/pb}";

	cout << "\\\\\\hline" << endl;

	cout << "     Cut       " << " &" << setw(13) << "\\ttbar{} " << " &" << setw(13) << "W+jets " << " &" << setw(13)
			<< "Z+jets " << " &" << setw(13) << "QCD " << " &" << setw(13) << "VQQ " << " &" << setw(13) << "Single Top " << " &"
			<< setw(13) << "Z'(500, 5) " << " &" << setw(25) << "Total  \\\\\n\\hline" << endl;

	double totalT; //total for a particular mc type
	double totalA; //total for all event types at each cut stage

	short njbegin = 0;

	for (short i = 0; i < 11; ++i) {//loop over cuts (up to HT)

		totalA = 0; //reset to zero for each cut
		printCutStage(i, ve.at(i));

		//Signal
		totalT = getTotalEvents(nevent, i, kttbar, kttbar + NUMBER_OF_DECAY_BRANCHES, njbegin);//TODO: check this
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;//add signal to total

		//Wjets
		totalT = getTotalEvents(nevent, i, kWjets, kWjets, njbegin);
		;
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;//add wjets to total

		//Zjets
		totalT = getTotalEvents(nevent, i, kZjets, kZjets, njbegin);
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;//add zjets to total

		//QCD//loop over QCD mc types
		totalT = getTotalEvents(nevent, i, kenri1, kbce3, njbegin);
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;//add QCD to total

		//VQQ
		totalT = getTotalEvents(nevent, i, kVQQ, kVQQ, njbegin);
		;
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;//add VQQ to total

		//single top
		totalT = getTotalEvents(nevent, i, ktW, kschan, njbegin);//loop over single top mc types (20-22)
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT; //add single top to total

		//FIXME 17th step to new sample (optional)
		totalT = getTotalEvents(nevent, i, Zprime_M500GeV_W5GeV, Zprime_M500GeV_W5GeV, njbegin);
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT; //add Z' M500GeV_W5GeV to total

		//print total column:
		cout << " & " << setw(13) << fixed << totalA << " \\\\" << endl;

		//after muon-veto, place nj>=4 cut
		if (ve.at(i) == "!MUON") {
			njbegin = 4;
			ve.at(i) = "$\\ge$4 jets";
			i--;
		}

	}// end loop over cut (nstage)
	cout << "\\hline" << endl;
	if (!first_time) {
		cout << "\\end{tabular}\\\\[5mm]" << endl;
	}
	//restore ve
	ve = ve_initial;
	//Zprime table 1% width. TODO: make an own function for this

	if (first_time) {
		//cout << "\\\\[1em]\n" << endl;
		cout << "\\begin{tabular}{|l|rrrrrrrr|r|}" << endl;
		cout << "\\hline" << endl;
	}
	cout << "\\multicolumn{10}{|l|}";
	if (first_time)
		cout << "{Actual number of MC events passing selection}";
	else
		cout << "{Expected number of events for " << intlumi << "/pb}";
	cout << "\\\\\\hline" << endl;

	cout << "           Cut       " << " &" << setw(13) << "Z'(500,5) " << " &" << setw(13) << "Z'(750,7.5) " << " &" << setw(13)
			<< "Z'(1000,10)" << " &" << setw(13) << "Z'(1250,12.5) " << " &" << setw(13) << "Z'(1500,15) " << " &" << setw(13)
			<< "Z'(2000,20)" << " &" << setw(13) << "Z'(3000,30) " << " &" << setw(13) << "Z'(4000,40) " << " &" << setw(25)
			<< "Total  \\\\\n\\hline" << endl;
	njbegin = 0;

	for (short i = 0; i < 11; ++i) {//loop over cuts (up to HT)
		//		if(i == 0){
		//		for(short xxx = 0; xxx < NUMBER_OF_MC_TYPES + 7; ++xxx){
		//			cout << "showing number " << xxx << " : ";
		//			cout << getTotalEvents(nevent, i, xxx,xxx+1, 0) << endl;
		//		}}
		totalA = 0; //reset to zero for each cut
		printCutStage(i, ve.at(i));
		totalT = getTotalEvents(nevent, i, Zprime_M500GeV_W5GeV, Zprime_M500GeV_W5GeV, njbegin);//Z' M500GeV_W5GeV
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;//add Z' M500GeV_W5GeV to total

		totalT = getTotalEvents(nevent, i, Zprime_M750GeV_W7500MeV, Zprime_M750GeV_W7500MeV, njbegin);//750, 7.5
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;

		totalT = getTotalEvents(nevent, i, Zprime_M1TeV_W10GeV, Zprime_M1TeV_W10GeV, njbegin);//1000, 10
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;

		totalT = getTotalEvents(nevent, i, Zprime_M1250GeV_W12500MeV, Zprime_M1250GeV_W12500MeV, njbegin);//1250, 12.5
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;

		totalT = getTotalEvents(nevent, i, Zprime_M1500GeV_W15GeV, Zprime_M1500GeV_W15GeV, njbegin);//1500, 15
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;

		totalT = getTotalEvents(nevent, i, Zprime_M2TeV_W20GeV, Zprime_M2TeV_W20GeV, njbegin);//2000, 20
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;

		totalT = getTotalEvents(nevent, i, Zprime_M3TeV_W30GeV, Zprime_M3TeV_W30GeV, njbegin);//3000, 30
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;

		totalT = getTotalEvents(nevent, i, Zprime_M4TeV_W40GeV, Zprime_M4TeV_W40GeV, njbegin);//4000, 40
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;

		//print total column:
		cout << " & " << setw(13) << fixed << totalA << " \\\\" << endl;

		//after muon-veto, place nj>=4 cut
		if (ve.at(i) == "!MUON") {
			njbegin = 4;
			ve.at(i) = "$\\ge$4 jets";
			i--;
		}

	}// end loop over cut (nstage)
	//restore ve
	ve = ve_initial;
	cout << "\\hline" << endl;
	if (!first_time) {
		cout << "\\end{tabular}\\\\[5mm]" << endl;
	}
	njbegin = 0;
	//Zprime table 10% width. TODO: make an own function for this
	if (first_time) {
		//cout << "\\\\[1em]\n" << endl;
		cout << "\\begin{tabular}{|l|rrrrrr|r|}" << endl;
		cout << "\\hline" << endl;
	}
	cout << "\\multicolumn{8}{|l|}";
	if (first_time)
		cout << "{Actual number of MC events passing selection}";
	else
		cout << "{Expected number of events for " << intlumi << "/pb}";
	cout << "\\\\\\hline" << endl;

	cout << "           Cut       " << " &" << setw(13) << "Z'(500,50) " << " &" << setw(13) << "Z'(1000,100)" << " &"
			<< setw(13) << "Z'(1500,150)" << " &" << setw(13) << "Z'(2000,200)" << " &" << setw(13) << "Z'(3000,300) " << " &"
			<< setw(13) << "Z'(4000,400) " << " &" << setw(25) << "Total  \\\\\n\\hline" << endl;

	for (short i = 0; i < 11; ++i) {//loop over cuts (up to HT)
		totalA = 0; //reset to zero for each cut
		printCutStage(i, ve.at(i));
		totalT = getTotalEvents(nevent, i, Zprime_M500GeV_W50GeV, Zprime_M500GeV_W50GeV, njbegin);//Z' M500GeV_W50GeV
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;//add Z' M500GeV_W5GeV to total

		totalT = getTotalEvents(nevent, i, Zprime_M1TeV_W100GeV, Zprime_M1TeV_W100GeV, njbegin);//1000, 100
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;

		totalT = getTotalEvents(nevent, i, Zprime_M1500GeV_W150GeV, Zprime_M1500GeV_W150GeV, njbegin);//1500, 150
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;

		totalT = getTotalEvents(nevent, i, Zprime_M2TeV_W200GeV, Zprime_M2TeV_W200GeV, njbegin);//2000, 200
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;

		totalT = getTotalEvents(nevent, i, Zprime_M3TeV_W300GeV, Zprime_M3TeV_W300GeV, njbegin);//3000, 300
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;

		totalT = getTotalEvents(nevent, i, Zprime_M4TeV_W400GeV, Zprime_M4TeV_W400GeV, njbegin);//4000, 400
		cout << " & " << setw(12) << fixed << totalT;
		totalA += totalT;

		//print total column:
		cout << " & " << setw(13) << fixed << totalA << " \\\\" << endl;
		if (ve.at(i) == "!MUON") {
			njbegin = 4;
			ve.at(i) = "$\\ge$4 jets";
			i--;
		}
	}// end loop over cut (nstage)
	cout << "\\hline" << endl;
	if (!first_time) {
		cout << "\\end{tabular}\\\\[5mm]" << endl;
	}
	first_time = false;
}
//end DrawMCTypeTable
//---------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
void ana::DrawSignalBGTable(const double nevent[][5][nmctype], vector<string> ve) const {

	cout << "\n%---------------------------------------------------------------------\n";
	cout << "       Expected Signal and Background for " << intlumi << "/pb";
	cout << "\n%---------------------------------------------------------------------\n\n";

	cout << "\\begin{tabular}{|l|r|rr|r|}\\hline" << endl;
	cout << "           Cut       " << " &  " << setw(11) << "Total Events" << " &  " << setw(14) << "Total Signal (S)" << " &  "
			<< setw(16) << "Total Background (B)" << " &  " << setw(10) << "S/B \\\\\n\\hline" << endl;

	short njbegin = 0;
	//short ntjet = 5;

	for (short i = 0; i < 11; ++i) { //loop over cuts

		double total_sig = 0;
		double total_bkg = 0;

		//Total Signal
		for (short k = 1; k < 11; ++k) { //loop over ttbar mc types (code 1 to 10)
			for (short j = njbegin; j < ntjet; ++j) {
				total_sig += nevent[i][j][k];
			} //sum up jet bins
		}

		//Total BG
		for (short k = 11; k < 23; ++k) { //loop over all bg mc types (code 11 to 22)
			for (short j = njbegin; j < ntjet; ++j) {
				total_bkg += nevent[i][j][k];
			} //sum up jet bins
		}

		//Print cut, S+B, S, B
		printCutStage(i, ve.at(i));
		cout << " & " << setw(13) << fixed << total_sig + total_bkg;
		cout << " & " << setw(13) << fixed << total_sig;
		cout << " & " << setw(13) << fixed << total_bkg;

		//Print Signal-to-background ratio (S/B)
		if (total_bkg > 0)
			cout << " & " << setw(11) << fixed << total_sig / total_bkg;
		else
			cout << " & " << setw(12) << "-";
		cout << " \\\\" << endl;

		// after mu-veto, place nj>=4 cut
		if (ve.at(i) == "!MUON") {
			njbegin = 4;
			ve.at(i) = "$\\ge$4 jets";
			i--;
		}

	}// end loop over cut (nstage)

	cout << "\\hline\n\\end{tabular}\n" << endl;
}
//end DrawSignalBGTable
//---------------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------------
// Draw event count table for break down of QCD
void ana::DrawQCDTable(const double nevent[14][5][nmctype], const string QCDtitle, vector<string> ve) const {

	static bool first_time = true;
	if (first_time)
		cout.precision(0); //unweighted table
	else
		cout.precision(myprec); //weighted table


	if (first_time)
		cout << "\\newpage\n" << endl;
	cout << "%------------------------------------------------------------------------" << endl;
	cout << "%                   " << QCDtitle << endl;
	cout << "%------------------------------------------------------------------------" << endl;

	if (first_time) {
		cout << "\\begin{tabular}{|l|rrrrrr|r|}" << endl;
		cout << "\\hline" << endl;
	}
	cout << "\\multicolumn{8}{|l|}";
	if (first_time)
		cout << "{Break down of actual number of QCD events passing selection}";
	else
		cout << "{Break down of expected QCD events for " << intlumi << "/pb}";
	cout << "\\\\\n\\hline" << endl;

	cout << "           Cut       ";
	cout << " &" << setw(13) << "enri1";
	cout << " &" << setw(13) << "enri2";
	cout << " &" << setw(13) << "enri3";
	cout << " &" << setw(13) << "bce1";
	cout << " &" << setw(13) << "bce2";
	cout << " &" << setw(13) << "bce3";
	cout << " &" << setw(23) << "AllQCD \\\\\\hline" << endl;

	//short ntjet = 5;
	short njbegin = 0;

	for (short i = 0; i < 11; ++i) { //cut stage
		double totalAllQCD = 0;
		printCutStage(i, ve.at(i));
		for (short k = 13; k < 19; ++k) { //mctype (QCD): 13-18
			double totalT = 0;
			for (short j = njbegin; j < ntjet; ++j) { //njet
				totalT += nevent[i][j][k];
			}
			totalAllQCD += totalT;
			cout << " & " << setw(12) << totalT;
		}
		cout << " & " << setw(13) << totalAllQCD << " \\\\" << endl;

		// insert >=4j cut
		if (ve.at(i) == "!MUON") {
			njbegin = 4;
			ve.at(i) = "$\\ge$4 jets";
			i--;
		}
	}
	cout << "\\hline" << endl;
	if (!first_time)
		cout << "\\end{tabular}\\\\[5mm]\n" << endl;
	first_time = false;

}//end DrawQCDTable
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// Draw event count table for break down of Single Top
void ana::DrawSingleTopTable(const double nevent[14][5][nmctype], const string title, vector<string> ve) const {

	static bool first_time = true;
	if (first_time)
		cout.precision(0); //unweighted table
	else
		cout.precision(myprec); //weighted table


	if (first_time)
		cout << "\\newpage\n" << endl;
	cout << "%------------------------------------------------------------------------" << endl;
	cout << "%                      " << title << endl;
	cout << "%------------------------------------------------------------------------" << endl;

	if (first_time) {
		cout << "\\begin{tabular}{|l|rrr|r|}" << endl;
		cout << "\\hline" << endl;
	}
	cout << "\\multicolumn{5}{|l|}";
	if (first_time)
		cout << "{Break down of actual number of single top events passing selection}";
	else
		cout << "{Break down of expected single top events for " << intlumi << "/pb}";
	cout << "\\\\\n\\hline" << endl;

	cout << "          Cut        ";
	cout << " &" << setw(13) << "tW-chan";
	cout << " &" << setw(13) << "t-chan";
	cout << " &" << setw(13) << "s-chan";
	cout << " &" << setw(23) << "AllSingleTop \\\\\\hline" << endl;

	//  short ntjet = 5;
	short njbegin = 0;

	for (short i = 0; i < 11; ++i) { //cut stage (up to HT)
		double allSingleTop = 0;
		printCutStage(i, ve.at(i));
		for (short k = 20; k < 23; ++k) { //mctype (single top): 20-22
			double totalT = 0;
			for (short j = njbegin; j < ntjet; ++j) { //njet
				totalT += nevent[i][j][k];
			}
			allSingleTop += totalT;
			cout << " & " << setw(12) << totalT;
		}
		cout << " & " << setw(13) << allSingleTop << " \\\\" << endl;
		//FIXME: reducing the iterating variable is not a good idea
		//insert >=4j cut
		if (ve.at(i) == "!MUON") {
			njbegin = 4;
			ve.at(i) = "$\\ge$4 jets";
			i--;
		}

	}
	cout << "\\hline" << endl;
	if (!first_time)
		cout << "\\end{tabular}\\\\[5mm]" << endl;
	first_time = false;

}//end DrawSingleTopTable
//---------------------------------------------------------------------------------------------


void ana::DrawSignalAcceptanceTable(const double nevent[][5][nmctype], vector<string> ve) const {

	cout.precision(0); //reset precision

	vector<string> tt;
	tt.push_back("&   Total ");
	tt.push_back("&    evqq ");
	tt.push_back("&    mvqq ");
	tt.push_back("&    tvqq ");
	tt.push_back("&    evev ");
	tt.push_back("&    mvmv ");
	tt.push_back("&    tvtv ");
	tt.push_back("&    evmv ");
	tt.push_back("&    evtv ");
	tt.push_back("&    mvtv ");
	tt.push_back("&    qqqq ");

	for (short h = 3; h < 5; h++) {

		TString njets_text = "$\\ge$"; //"$>=$";
		njets_text += h;
		njets_text += " jets";

		if (h == 3)
			cout << endl << "\\newpage";
		cout << endl << "%------------------------------------------------------------------";
		cout << endl << "%------------------------------------------------------------------";
		cout << endl << "      Acceptance for ttbar MC with " << njets_text << "  (unweighted)";
		cout << endl << "%------------------------------------------------------------------";
		cout << endl << "%------------------------------------------------------------------";
		cout << endl << "\\begin{tabular}{|l|r|rrrrrrrrrr|}\\hline";

		cout << endl << "           Cut        ";
		for (unsigned int k = 0; k < tt.size(); ++k)
			cout << tt.at(k);
		cout << "\\\\\n\\hline" << endl;

		int startnj = 0;
		for (short i = 0; i < 11; ++i) { //nstage

			printCutStage(i, ve.at(i));

			// calculate the sum of all ttbar first
			double total_tt = 0;
			for (short k = 1; k < 11; ++k) {
				for (short j = startnj; j < ntjet; ++j)
					total_tt += nevent[i][j][k];
			}
			// now print
			cout << " &  " << setw(6) << total_tt;
			for (short k = 1; k < 11; ++k) {
				double total = 0;
				for (short j = startnj; j < ntjet; ++j)
					total += nevent[i][j][k];
				cout << " &  " << setw(6) << total;
			}
			cout << " \\\\ ";

			if (ve.at(i) == "!MUON") { // extra step for >= n jet cut
				startnj = h;
				ve.at(i) = njets_text;
				i--;
			}
			cout << endl;
		}
		cout << "\\hline\n\\end{tabular}" << endl;
		if (h == 3) {
			cout << "\\\\[5mm]" << endl;
			ve.at(6) = "!MUON"; //reset
		}
	}

}//end DrawSignalAcceptanceTable
//------------------------------------------------------------------------------------


void ana::SetHistoLabelCutNjet(TH2D *this_njetsVcuts, vector<string>& ve) const {

	// ve.size() = nstage
	for (size_t i = 0; i < ve.size(); ++i) {
		this_njetsVcuts->GetYaxis()->SetBinLabel(i + 1, ve.at(ve.size() - i - 1).c_str()); //i=0, bin=1, stage=13, cut=-1btag
	}
	this_njetsVcuts->GetXaxis()->SetBinLabel(1, "0-jet");
	this_njetsVcuts->GetXaxis()->SetBinLabel(2, "1-jet");
	this_njetsVcuts->GetXaxis()->SetBinLabel(3, "2-jets");
	this_njetsVcuts->GetXaxis()->SetBinLabel(4, "3-jets");
	this_njetsVcuts->GetXaxis()->SetBinLabel(5, "#geq4-jets");
	this_njetsVcuts->GetXaxis()->SetBinLabel(6, "Total");
	this_njetsVcuts->SetLabelSize(0.07, "x");
	this_njetsVcuts->GetYaxis()->SetLabelOffset(0.0);
	this_njetsVcuts->SetOption("text");
	this_njetsVcuts->SetMarkerSize(2.0);
	this_njetsVcuts->SetStats(0);

}//end SetHistoLabelCutNjet
//---------------------------------------------------------------------------------------------

void ana::SetHistoLabelEleID(TH1F *eid[]) const {

	short nhisto = 1; //data
	if (!IsData())
		nhisto = mcsize;
	for (short i = 0; i < nhisto; ++i) {
		if (i > 0 && hasSampleMC(i) == false)
			continue;
		eid[i]->GetXaxis()->SetBinLabel(1, "allEle");
		eid[i]->GetXaxis()->SetBinLabel(2, "ET/Eta");
		eid[i]->GetXaxis()->SetBinLabel(3, "d0");
		eid[i]->GetXaxis()->SetBinLabel(4, "eidLoose");
		eid[i]->GetXaxis()->SetBinLabel(5, "eidTight");
		eid[i]->GetXaxis()->SetBinLabel(6, "eidRobustLoose");
		eid[i]->GetXaxis()->SetBinLabel(7, "eidRobustTight");
		eid[i]->SetLabelSize(0.04, "x");
		eid[i]->SetOption("hist text0");
		eid[i]->SetMarkerSize(1.5);
		eid[i]->SetStats(0);
	}
}//end SetHistoLabelEleID
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------

//--------- compute (signed) d0 corrected w.r.t Beam Spot -------------------------------------
float ana::compute_d0(const string lepton, const int i) const {

	float vx = 0;
	float vy = 0;
	float px = 0;
	float py = 0;

	if (lepton == "electron") {
		vx = els_vx->at(i);
		vy = els_vy->at(i);
		px = els_vpx->at(i);
		py = els_vpy->at(i);
	} else if (lepton == "muon") {
		vx = mus_tk_vx->at(i);
		vy = mus_tk_vy->at(i);
		px = mus_tk_px->at(i);
		py = mus_tk_py->at(i);
	} else { //track
		vx = tracks_vx->at(i);
		vy = tracks_vy->at(i);
		px = tracks_px->at(i);
		py = tracks_py->at(i);
	}

	float refx = beamSpot_x->at(0);
	float refy = beamSpot_y->at(0);

	float d0_corrected = -(-(vx - refx) * py + (vy - refy) * px) / TMath::Hypot(px, py); //d0 = -dxy

	return d0_corrected;
}//end compute_d0
//---------------------------------------------------------------------------------------------

//--------- compute mT(W) given electron 4-vector and MET -------------------------------------
float ana::compute_mtw(const TVector2& e, const TVector2& miss) const {

	float leppx = e.Px();
	float leppy = e.Py();
	float lepet = e.Mod();
	float mspx = miss.Px();
	float mspy = miss.Py();
	float mset = miss.Mod();
	float z1 = pow((lepet + mset), 2);
	float z2 = pow((leppx + mspx), 2) + pow((leppy + mspy), 2);

	float z3 = z1 - z2;
	float mtlm = -1.0;

	if (z3 > 0) {
		mtlm = sqrt(z3);
	}
	return mtlm;

}//end compute_mtw
//---------------------------------------------------------------------------------------------

void ana::PrintGenParticles() const {
	cout << setfill('-') << setw(105) << "" << setfill(' ') << endl;
	cout << setw(3) << "" << setw(8) << "PDG id" << setw(8) << "Status" << setw(10) << "Parent id" << setw(15) << "Mass" << setw(
			12) << "Px" << setw(12) << "Py" << setw(12) << "Pz" << setw(12) << "Pt" << setw(12) << "Energy" << endl;
	cout << setfill('-') << setw(105) << "" << setfill(' ') << endl;
	for (unsigned int i = 0; i < Nmc_doc; ++i) {
		cout << setw(3) << i << setw(8) << mc_doc_id->at(i) << setw(8) << mc_doc_status->at(i) << setw(10)
				<< mc_doc_mother_id->at(i) << setw(15) << mc_doc_mass->at(i) << setw(12) << mc_doc_px->at(i) << setw(12)
				<< mc_doc_py->at(i) << setw(12) << mc_doc_pz->at(i) << setw(12) << mc_doc_pt->at(i) << setw(12)
				<< mc_doc_energy->at(i) << endl;
	}
	cout << setfill('-') << setw(105) << "" << setfill(' ') << endl;
}//end PrintGenParticles
//---------------------------------------------------------------------------------------------

void ana::printCutStage(int i, string cut) const {
	cout << " Stage " << setw(2) << i << " " << setw(11) << left << cut << right;
}
//---------------------------------------------------------------------------------------------

void ana::printCutStage(ofstream& os, int i, string cut) const {
	os << " Stage " << setw(2) << i << " " << setw(11) << left << cut << right;
}
//---------------------------------------------------------------------------------------------

// *** NOT READY FOR USED ***
string ana::CheckEventTypeFromMcTruth() const {
	bool found_top = false;
	bool found_W = false;
	bool found_Z = false;
	// Look only at "documentation lines" (in Pythia), ie status=3 particles
	for (unsigned int i = 0; i < Nmc_doc; ++i) {
		if (mc_doc_status->at(i) != 3)
			break;
		if (fabs(mc_doc_id->at(i)) == 6)
			found_top = true;
		else if (fabs(mc_doc_id->at(i)) == 24)
			found_W = true;
		else if (fabs(mc_doc_id->at(i)) == 22)
			found_Z = true;
	}
	if (found_top)
		return "ttbar";
	if (found_W)
		return "wenu";
	if (found_Z)
		return "zee";
	return "none";
}
//---------------------------------------------------------------------------------------------

float ana::getRelIso(int i) const {
	return (els_dr04EcalRecHitSumEt->at(i) + els_dr04HcalTowerSumEt->at(i) + els_tIso->at(i)) / els_et->at(i);
}
//---------------------------------------------------------------------------------------------

bool ana::passEleID(unsigned int i) const {
	if (EleID() == robustTight)
		return (els_robustTightId->at(i) > 0);
	if (EleID() == robustLoose)
		return (els_robustLooseId->at(i) > 0);
	if (EleID() == tight)
		return (els_tightId->at(i) > 0);
	if (EleID() == loose)
		return (els_looseId->at(i) > 0);
	if (EleID() == robustHighenergy)
		return (els_robustHighEnergyId->at(i) > 0);
	if (EleID() == none)
		return true;
	return false;
}
//---------------------------------------------------------------------------------------------

string ana::printEleID() const {
	switch (EleID()) {
	case robustTight:
		return "robustTight";
	case robustLoose:
		return "robustLoose";
	case tight:
		return "tight";
	case loose:
		return "loose";
	case none:
		return "none";
	default:
		return "robustTight (Default)";
	}
}
//---------------------------------------------------------------------------------------------

bool ana::passHLT() const {
	if (HLTBit == "HLT_Ele15_LW_L1R")
		return (bool) HLT_Ele15_LW_L1R;
	if (HLTBit == "HLT_Ele15_SW_L1R")
		return (bool) HLT_Ele15_SW_L1R;
	return false;
}
//---------------------------------------------------------------------------------------------

string ana::printTimeNow() const {
	TDatime now;
	return now.AsSQLString();
}

/**
 * Creating histogram names for 1D histograms
 */

void ana::bookHistograms() {
	fasthist_.resize(NUMBER_OF_MC_TYPES);
	fasthist2D_.resize(NUMBER_OF_MC_TYPES);
	fasthist_jetbinned_.resize(NUMBER_OF_MC_TYPES);
	fasthistNlevel1D_.resize(NUMBER_OF_MC_TYPES);
	fasthistNlevel2D_.resize(NUMBER_OF_MC_TYPES);
	for (ushort type = 0; type < (short) NUMBER_OF_MC_TYPES; ++type) {
		if (!hasSampleMC(type))
			continue;
		fasthist_[type].resize(NUMBER_OF_1D_HISTOGRAMS);
		fasthist2D_[type].resize(NUMBER_OF_2D_HISTOGRAMS);
		fasthist_jetbinned_[type].resize(NUMBER_OF_JETBINS);
		fasthist_jetbinned_[type].resize(NUMBER_OF_JETBINS);
		for (ushort jbin = 0; jbin < NUMBER_OF_JETBINS; ++jbin) {
			fasthist_jetbinned_[type][jbin].resize(NUMBER_OF_JETBINNED_HISTOGRAMS);
		}
		fasthistNlevel1D_[type].resize(NUMBER_OF_LEVELS);
		fasthistNlevel2D_[type].resize(NUMBER_OF_LEVELS);
		for (ushort level = 0; level < NUMBER_OF_LEVELS; ++level) {
			fasthistNlevel1D_[type][level].resize(NUMBER_OF_JETBINS);
			fasthistNlevel2D_[type][level].resize(NUMBER_OF_JETBINS);
			for (ushort jbin = 0; jbin < NUMBER_OF_JETBINS; ++jbin) {
				fasthistNlevel1D_[type][level][jbin].resize(NUMBER_OF_NLEVEL_JBINNED_1D_HISTOGRAMS);
				fasthistNlevel2D_[type][level][jbin].resize(NUMBER_OF_NLEVEL_JBINNED_2D_HISTOGRAMS);
			}
		}
		gROOT->cd();
		TDirectory *dir = histf->mkdir(mclabel[type].c_str());
		dir->cd();
		bookBtagHistograms(type, dir);
		TDirectory *exo = dir->mkdir("exotic_top", "Exotic Top searches");
		exo->cd();
		bookExoticTopMatchedHistograms(type);
		bookExoticTopMCHistograms(type);
		bookExoticTopHistograms(type);
		bookExoticTopBtaggedHistograms(type);
		dir->cd();
		//
		TDirectory *basic = dir->mkdir("basic", "Basic kinematics");
		basic->cd();
		bookBasicHistograms(type);
		dir->cd();

		TDirectory *dir_DRemu = dir->mkdir("DRemu", "Delta R(e,mu)");
		dir_DRemu->cd();
		bookDRemuHistograms(type);
		dir->cd();

		TDirectory *dir_explore = dir->mkdir("explore", "explore new var after all but MET cuts (>=4j)");
		dir_explore->cd();
		bookExploreHistograms(type);
		dir->cd();

		//-------------------
		//  Study of Z veto
		//-------------------
		bookZvetoHistograms(type, dir);

		// Electron d0 plots (bin: 0.5cm for 5000 bins, ie 1 micron per bin)
		TDirectory *dir_ed0_unCor = dir->mkdir("electron_d0_unCor", "electron d0 (uncorrected)");
		dir_ed0_unCor->cd();
		addHistoDataAndMC(type, h_ed0_unCor, "ed0_unCor", "ele |d0| (pass ET,eta cut)", 3000, 0, 0.3);

		TDirectory *dir_ed0 = dir->mkdir("electron_d0", "electron d0 (corrected using BeamSpot)");
		dir_ed0->cd();
		addHistoDataAndMC(type, h_ed0, "ed0", "ele |d0| (pass ET,eta cut)", 3000, 0, 0.3);
		addHistoDataAndMC(type, h_ed0_pass, "ed0_pass", "ele |d0| (pass ET,eta cut)", 3000, 0, 0.3);
		dir->cd();
		addHistoDataAndMC(type, h_muon_chi2, "muon_chi2", "muon chi2 (normalized) (pass pt,eta cut)", 100, 0, 20);
		addHistoDataAndMC(type, h_muon_d0_unCor, "muon_d0_unCor", "muon d0 (uncorrected )(pass pt,eta cut) ", 1000, -0.5, 0.5);
		addHistoDataAndMC(type, h_muon_d0, "muon_d0", "muon d0 (corrected using BeamSpot) (pass pt,eta cut)", 1000, -0.5, 0.5);
		addHistoDataAndMC(type, h_muon_hits, "muon_hits", "muon hits (pass pt,eta cut)", 40, 0, 40);

		TDirectory *dir_MET = dir->mkdir("MET", "MET after all but MET cut");
		dir_MET->cd();
		addHistoDataAndMC(type, h_met_ante_ISO, "met_ante_ISO", "Missing ET (1 GooEle, before ISO cut)", 200, 0, 200);
		addHistoDataAndMC(type, h_met_ante_ISO_mu, "met_ante_ISO_mu", "muMET (1 GooEle, before ISO cut)", 200, 0, 200);
		addHistoDataAndMC(type, h_met_ante_ISO_t1, "met_ante_ISO_t1", "Type1MET (1 GooEle, before ISO cut)", 200, 0, 200);
		dir->cd();
		TDirectory *dir_mtw = dir->mkdir("mtw", "mT(W) W transverse mass");
		dir_mtw->cd();
		addHistoDataAndMC(type, h_mtw_mu_incl, "mtw_mu_incl", "mT(W) #mu-cor caloMet (incl, leading isolated ele)", 100, 0, 200);
		addHistoDataAndMC(type, h_mtw_t1_incl, "mtw_t1_incl", "mT(W) type1 caloMet (incl, leading isolated ele)", 100, 0, 200);
		dir->cd();

		TDirectory *dir_DPhiEmet = dir->mkdir("DPhi_ele_met", "Delta phi(isoele,met)");
		dir_DPhiEmet->cd();
		addHistoDataAndMC(type, h_DPhiEmet_mu_incl, "DPhiEmet_mu_incl",
				"#Delta#Phi(e,met) #mu-cor caloMet (incl, leading isolated ele)", 64, -3.2, 3.2);
		addHistoDataAndMC(type, h_DPhiEmet_t1_incl, "DPhiEmet_t1_incl",
				"#Delta#Phi(e,met) type 1 caloMet (incl, leading isolated ele)", 64, -3.2, 3.2);
		dir->cd();

		//jet binned histograms
		TDirectory *dir_nEle = dir->mkdir("electron_count");
		TDirectory *dir_HT = dir->mkdir("HT", "HT after all but HT cut");

		dir_nEle->cd();
		bookElectronCountHistograms(type);
		dir_HT->cd();
		addHisto_Njet_DataAndMC(type, h_ht, "HT", "HT (after all but HT cut)", 100, 0, 1000);
		dir_MET->cd();
		bookNewKinematicHistograms(type);
		dir_mtw->cd();
		addHisto_Njet_DataAndMC(type, h_mtw_mu, "mtw_mu", "mT(W) #mu-cor caloMet (after all but MET cut)", 100, 0, 200);
		addHisto_Njet_DataAndMC(type, h_mtw_t1, "mtw_t1", "mT(W) type1 caloMet (after all but MET cut)", 100, 0, 200);
		dir_DPhiEmet->cd();
		addHisto_Njet_DataAndMC(type, h_DPhiEmet_mu, "DPhiEmet_mu", "#Delta#Phi(e,met) #mu-cor caloMet (after all but MET cut)",
				64, -3.2, 3.2);
		addHisto_Njet_DataAndMC(type, h_DPhiEmet_t1, "DPhiEmet_t1", "#Delta#Phi(e,met) type 1 caloMet (after all but MET cut)",
				64, -3.2, 3.2);

		bookQCDEstimationHistograms(type, dir);

		for (unsigned short ihist = 0; ihist < NUMBER_OF_1D_HISTOGRAMS; ++ihist) {
			if (ihist >= h_nGenBasicEle_Zee_allj && ihist <= h_Zee_Nphotons && type != kZjets)
				continue;
			if (fasthist_[type][ihist] > 0) {

			} else {
				cout << "ERROR: 1D Histogram (" << ihist << "/" << NUMBER_OF_1D_HISTOGRAMS << ") '" /*<< histnames_[ihist]*/
				<< "' not initiated" << endl;
				if (ihist > 0)
					cout << "last good histogram was '" << fasthist_[type][ihist - 1]->GetName() << "'" << endl;
				exit(-1);
			}
		}
		for (ushort jetbin = 0; jetbin < NUMBER_OF_JETBINS; ++jetbin) {
			for (unsigned short ihist = 0; ihist < NUMBER_OF_JETBINNED_HISTOGRAMS; ++ihist) {
				if (fasthist_jetbinned_[type][jetbin][ihist] == NULL) {
					cout << "ERROR: Histogram (jetbin " << jetbin << ") (" << ihist << "/" << NUMBER_OF_JETBINNED_HISTOGRAMS
							<< ") '" << ihist << "' not initiated" << endl;
					exit(-1);
				}
			}
		}

		for (ushort level = 0; level < NUMBER_OF_LEVELS; ++level) {
			for (ushort jetbin = 0; jetbin < NUMBER_OF_JETBINS; ++jetbin) {
				for (unsigned short ihist = 0; ihist < NUMBER_OF_NLEVEL_JBINNED_1D_HISTOGRAMS; ++ihist) {
					if (fasthistNlevel1D_[type][level][jetbin][ihist] == NULL) {
						cout << "ERROR: 1D Histogram (level " << level << ", jetbin " << jetbin << ") (" << ihist << "/"
								<< NUMBER_OF_NLEVEL_JBINNED_1D_HISTOGRAMS << ") '" << ihist << "' not initiated" << endl;
						if (ihist > 0)
							cout << "last good histogram was '" << fasthistNlevel1D_[type][level][jetbin][ihist - 1]->GetName()
									<< "'" << endl;
						exit(-1);
					}
				}
				for (unsigned short ihist = 0; ihist < NUMBER_OF_NLEVEL_JBINNED_2D_HISTOGRAMS; ++ihist) {
					if (fasthistNlevel2D_[type][level][jetbin][ihist] == NULL) {
						cout << "ERROR: 2D Histogram (level " << level << ", jetbin " << jetbin << ") (" << ihist << "/"
								<< NUMBER_OF_NLEVEL_JBINNED_2D_HISTOGRAMS << ") '" << ihist << "' not initiated" << endl;
						if (ihist > 0)
							cout << "last good histogram was '" << fasthistNlevel2D_[type][level][jetbin][ihist - 1]->GetName()
									<< "'" << endl;
						exit(-1);
					}
				}
			}
		}
		//		}

		//2D histograms

		for (unsigned short ihist = 0; ihist < NUMBER_OF_2D_HISTOGRAMS; ++ihist) {
			if (ihist == h_Zee_photon_eteta_2D && type != kZjets)
				continue;
			if (fasthist2D_[type][ihist] > 0) {

			} else {
				cout << "WARNING: 2D Histogram (" << ihist << "/" << NUMBER_OF_2D_HISTOGRAMS << ") '" /*<< histnames2D_[ihist]*/
				<< "' not initiated" << endl;
				exit(-1);
			}
		}

	}
	gROOT->cd();
}

void ana::bookQCDEstimationHistograms(ushort type, TDirectory *parent) {
	TDirectory *dir_QCD = parent->mkdir("QCD_estimation");
	dir_QCD->cd();
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso, "QCDest_CombRelIso", "RelIso", 1000, 0, 10);
	//--------------------------------
	// AES (Anti Event Selection)
	//--------------------------------
	TDirectory *dir_AES = dir_QCD->mkdir("AES", "anti event selection");
	dir_AES->cd();
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES, "QCDest_CombRelIso_AES", "RelIso (AES)", 1000, 0, 10);//ALL AES
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_minusMET, "QCDest_CombRelIso_AES_minusMET", "RelIso (AES-met)", 1000,
			0, 10);// AES except MET cut
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_minusHT, "QCDest_CombRelIso_AES_minusHT", "RelIso (AES-HT)", 1000, 0,
			10);// AES except HT cut
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_minusTighterZ, "QCDest_CombRelIso_AES_minusTighterZ",
			"RelIso (AES-TighterZ)", 1000, 0, 10);// AES except tighter Z veto cut (both mee and mep)
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_before, "QCDest_CombRelIso_AES_before", "RelIso (Before AES)", 1000, 0,
			10);// Before AES
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_justMET, "QCDest_CombRelIso_AES_justMET", "RelIso (AES: just MET)",
			1000, 0, 10);// AES: just MET (<x)
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_justHighMET, "QCDest_CombRelIso_AES_justHighMET",
			"RelIso (AES: just highMET)", 1000, 0, 10);// AES: just high MET (>x)
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_justHT, "QCDest_CombRelIso_AES_justHT", "RelIso (AES: just HT)", 1000,
			0, 10); // AES: just HT
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_justZ, "QCDest_CombRelIso_AES_justZ", "RelIso (AES: just tightZ)",
			1000, 0, 10);// AES: just Tighter Z (both)
	// 4-12-09: Trial definitions of AES for Plan B
	TDirectory *dir_QCD_planA = dir_QCD->mkdir("QCD_planA");
	dir_QCD_planA->cd();
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planA1_e20, "QCDest_CombRelIso_AES_planA1_e20",
			"RelIso (AES A1, no conv cut, E_{T}>20)", 1000, 0, 10);
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planA1_e30, "QCDest_CombRelIso_AES_planA1_e30",
			"RelIso (AES A1, no conv cut, E_{T}>30)", 1000, 0, 10);
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planA2_e20, "QCDest_CombRelIso_AES_planA2_e20",
			"RelIso (AES A2, conv, E_{T}>20)", 1000, 0, 10);
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planA2_e30, "QCDest_CombRelIso_AES_planA2_e30",
			"RelIso (AES A2, conv, E_{T}>30)", 1000, 0, 10);
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planA3_e20, "QCDest_CombRelIso_AES_planA3_e20",
			"RelIso (AES A3 |d_{0}|>200um, E_{T}>20)", 1000, 0, 10);
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planA3_e30, "QCDest_CombRelIso_AES_planA3_e30",
			"RelIso (AES A3 |d_{0}|>200um, E_{T}>30)", 1000, 0, 10);

	// 10-11-09: Trial definitions of AES for Plan B
	TDirectory *dir_QCD_planB = dir_QCD->mkdir("QCD_planB");
	dir_QCD_planB->cd();
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planB1_e20, "QCDest_CombRelIso_AES_planB1_e20",
			"RelIso (AES B1 E_{T}>20)", 1000, 0, 10);//EleET >x
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planB1_e30, "QCDest_CombRelIso_AES_planB1_e30",
			"RelIso (AES B1 E_{T}>30)", 1000, 0, 10);
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planB2_e20, "QCDest_CombRelIso_AES_planB2_e20",
			"RelIso (AES B2 E_{T}<20)", 1000, 0, 10);//EleET < x
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planB2_e30, "QCDest_CombRelIso_AES_planB2_e30",
			"RelIso (AES B2 E_{T}<30)", 1000, 0, 10);
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planB3_e20, "QCDest_CombRelIso_AES_planB3_e20",
			"RelIso (AES B3 RT=0 E_{T}>20)", 1000, 0, 10);//fail RT ID
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planB3_e30, "QCDest_CombRelIso_AES_planB3_e30",
			"RelIso (AES B3 RT=0 E_{T}>30)", 1000, 0, 10);
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planB3b_e20, "QCDest_CombRelIso_AES_planB3b_e20",
			"RelIso (AES B3 RT=0 EB E_{T}>20)", 1000, 0, 10);//new//fail RT ID (BARREL)
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planB3b_e30, "QCDest_CombRelIso_AES_planB3b_e30",
			"RelIso (AES B3 RT=0 EB E_{T}>30)", 1000, 0, 10);//new
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planB4_e20, "QCDest_CombRelIso_AES_planB4_e20",
			"RelIso (AES B4 RL=0 E_{T}>20)", 1000, 0, 10);//fail RL ID
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planB4_e30, "QCDest_CombRelIso_AES_planB4_e30",
			"RelIso (AES B4 RL=0 E_{T}>30)", 1000, 0, 10);
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planB5_e20, "QCDest_CombRelIso_AES_planB5_e20",
			"RelIso (AES B5 cL=0 E_{T}>20)", 1000, 0, 10);//fail Loose ID
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planB5_e30, "QCDest_CombRelIso_AES_planB5_e30",
			"RelIso (AES B5 cL=0 E_{T}>30)", 1000, 0, 10);
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planB6_e20, "QCDest_CombRelIso_AES_planB6_e20",
			"RelIso (AES B6 cT=0 E_{T}>20)", 1000, 0, 10);//fail Tight ID
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planB6_e30, "QCDest_CombRelIso_AES_planB6_e30",
			"RelIso (AES B6 cT=0 E_{T}>30)", 1000, 0, 10);
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planB7_e20, "QCDest_CombRelIso_AES_planB7_e20",
			"RelIso (AES B7 |d_{0}|>200um RT=1 E_{T}>20)", 1000, 0, 10);//d0 > 200um, pass RT
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planB7_e30, "QCDest_CombRelIso_AES_planB7_e30",
			"RelIso (AES B7 |d_{0}|>200um RT=1 E_{T}>30)", 1000, 0, 10);
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planB8_e20, "QCDest_CombRelIso_AES_planB8_e20",
			"RelIso (AES B8 |d_{0}|>200um RT=0 E_{T}>20)", 1000, 0, 10);//d0 > 200um, fail RT
	addHisto_Njet_DataAndMC(type, h_QCDest_CombRelIso_AES_planB8_e30, "QCDest_CombRelIso_AES_planB8_e30",
			"RelIso (AES B8 |d_{0}|>200um RT=0 E_{T}>30)", 1000, 0, 10);

	//-------------------------------
	// NES (Normal Event Selection)
	//-------------------------------
	// 4 levels: after HLT(all electrons); nGoodEle>0; muZveto; conv
	// notes: [4]:levels, [7]:njet
	// scatter plot, iso:met
	if (m_plotRelisoNES) {
		TDirectory *dir_NES = dir_QCD->mkdir("NES", "normal event selection");
		dir_NES->cd();

		iso_addHisto_nlevel_nj_nmc(type, h_QCDest_isoVmet_NES, "QCDest_isoVmet_NES", "RelIso v MET (NES)", 70, 0, 1.4, 50, 0, 50);//iso,met
		iso_addHisto_nlevel_nj_nmc(type, h_QCDest_isoVmet_NES_barrel, "QCDest_isoVmet_NES_barrel", "RelIso v MET (NES barrel)",
				70, 0, 1.4, 50, 0, 50);
		iso_addHisto_nlevel_nj_nmc(type, h_QCDest_isoVmet_NES_endcap, "QCDest_isoVmet_NES_endcap", "RelIso v MET (NES endcap)",
				70, 0, 1.4, 50, 0, 50);
		// unweighted
		iso_addHisto_nlevel_nj_nmc(type, h_QCDest_isoVmet_NES_uw, "QCDest_isoVmet_NES_uw", "RelIso v MET (NES uw)", 70, 0, 1.4,
				50, 0, 50);
		iso_addHisto_nlevel_nj_nmc(type, h_QCDest_isoVmet_NES_uw_barrel, "QCDest_isoVmet_NES_uw_barrel",
				"RelIso v MET (NES uw barrel)", 70, 0, 1.4, 50, 0, 50);
		iso_addHisto_nlevel_nj_nmc(type, h_QCDest_isoVmet_NES_uw_endcap, "QCDest_isoVmet_NES_uw_endcap",
				"RelIso v MET (NES uw endcap)", 70, 0, 1.4, 50, 0, 50);

		iso_addHisto_nlevel_nj_nmc(type, h_QCDest_CombRelIso_NES, "QCDest_CombRelIso_NES", "RelIso (NES)", 200, 0, 2.0);
		iso_addHisto_nlevel_nj_nmc(type, h_QCDest_CombRelIso_NES_barrel, "QCDest_CombRelIso_NES_barrel", "RelIso (NES barrel)",
				200, 0, 2.0);
		iso_addHisto_nlevel_nj_nmc(type, h_QCDest_CombRelIso_NES_endcap, "QCDest_CombRelIso_NES_endcap", "RelIso (NES endcap)",
				200, 0, 2.0);

		iso_addHisto_nlevel_nj_nmc(type, h_QCDest_CombRelIso_NES_loMET, "QCDest_CombRelIso_NES_loMET", "RelIso (NES loMET)", 200,
				0, 2.0);
		iso_addHisto_nlevel_nj_nmc(type, h_QCDest_CombRelIso_NES_loMET_barrel, "QCDest_CombRelIso_NES_loMET_barrel",
				"RelIso (NES loMET barrel)", 200, 0, 2.0);
		iso_addHisto_nlevel_nj_nmc(type, h_QCDest_CombRelIso_NES_loMET_endcap, "QCDest_CombRelIso_NES_loMET_endcap",
				"RelIso (NES loMET endcap)", 200, 0, 2.0);

		iso_addHisto_nlevel_nj_nmc(type, h_QCDest_CombRelIso_NES_hiMET, "QCDest_CombRelIso_NES_hiMET", "RelIso (NES hiMET)", 200,
				0, 2.0);
		iso_addHisto_nlevel_nj_nmc(type, h_QCDest_CombRelIso_NES_hiMET_barrel, "QCDest_CombRelIso_NES_hiMET_barrel",
				"RelIso (NES hiMET barrel)", 200, 0, 2.0);
		iso_addHisto_nlevel_nj_nmc(type, h_QCDest_CombRelIso_NES_hiMET_endcap, "QCDest_CombRelIso_NES_hiMET_endcap",
				"RelIso (NES hiMET endcap)", 200, 0, 2.0);
	}
	parent->cd();
}//end bookQCDEstimationHistograms
void ana::bookNewKinematicHistograms(ushort type) {
	addHisto_Njet_DataAndMC(type, h_met, "met", "MET (n-1)", 200, 0, 200);
	addHisto_Njet_DataAndMC(type, h_met_mu, "met_mu", "#muMET (n-1)", 200, 0, 200);
	addHisto_Njet_DataAndMC(type, h_met_t1, "met_t1", "Type1MET (n-1)", 200, 0, 200);
	addHisto_Njet_DataAndMC(type, h_met_BA, "met_BA", "MET (n-1 barrel)", 200, 0, 200);
	addHisto_Njet_DataAndMC(type, h_met_mu_BA, "met_mu_BA", "#muMET (n-1 barrel)", 200, 0, 200);
	addHisto_Njet_DataAndMC(type, h_met_t1_BA, "met_t1_BA", "Type1MET (n-1 barrel)", 200, 0, 200);
}
void ana::bookElectronCountHistograms(ushort type) {
	// for each njet and mctype
	addHisto_Njet_DataAndMC(type, h_nEle_all, "nEle_all", "N(e) all", 6, 0, 6);
	addHisto_Njet_DataAndMC(type, h_nEle_s1, "nEle_s1_EtEta", "N(e) pass ET,eta cuts", 6, 0, 6);
	addHisto_Njet_DataAndMC(type, h_nEle_s2, "nEle_s2_d0", "N(e) pass ET,eta,d0 cuts", 6, 0, 6);
	addHisto_Njet_DataAndMC(type, h_nEle_s3_idLoose, "nEle_s3_idLoose", "N(e) pass ET,eta,d0,looseId cuts", 6, 0, 6);
	addHisto_Njet_DataAndMC(type, h_nEle_s3_idTight, "nEle_s3_idTight", "N(e) pass ET,eta,d0,tightId cuts", 6, 0, 6);
	addHisto_Njet_DataAndMC(type, h_nEle_s3_idRL, "nEle_s3_idRobustLoose", "N(e) pass ET,eta,d0,robustLooseId cuts", 6, 0, 6);
	addHisto_Njet_DataAndMC(type, h_nEle_s3_idRT, "nEle_s3_idRobustTight", "N(e) pass ET,eta,d0,robustTightId cuts", 6, 0, 6);

	addHistoDataAndMC(type, h_eid, "eid", "Total N(ele) pass each cut", 7, 0, 7);//all, ET+eta, d0, loose, tight, RL, RT
	fasthist_[type][h_eid]->GetXaxis()->SetBinLabel(1, "allEle");
	fasthist_[type][h_eid]->GetXaxis()->SetBinLabel(2, "ET/Eta");
	fasthist_[type][h_eid]->GetXaxis()->SetBinLabel(3, "d0");
	fasthist_[type][h_eid]->GetXaxis()->SetBinLabel(4, "eidLoose");
	fasthist_[type][h_eid]->GetXaxis()->SetBinLabel(5, "eidTight");
	fasthist_[type][h_eid]->GetXaxis()->SetBinLabel(6, "eidRobustLoose");
	fasthist_[type][h_eid]->GetXaxis()->SetBinLabel(7, "eidRobustTight");
	fasthist_[type][h_eid]->SetLabelSize(0.04, "x");
	fasthist_[type][h_eid]->SetOption("hist text0");
	fasthist_[type][h_eid]->SetMarkerSize(1.5);
	fasthist_[type][h_eid]->SetStats(0);
}

void ana::bookExoticTopHistograms(ushort type) {
	addHistoDataAndMC(type, h_neutrino_pz, "neutrino_pz", "reconstructed p_z(#nu)", 50, -500, 500);
	addHistoDataAndMC(type, h_mttbar, "mttbar", "reconstructed M_{ttbar}", 10000, 0, 10000);
	addHistoDataAndMC(type, h_mWlep, "mWlep", "reconstructed leptonic W mass", 500, 0, 500);
	addHistoDataAndMC(type, h_mWhad, "mWhad", "reconstructed hadronic W mass", 500, 0, 500);
	addHistoDataAndMC(type, h_ptRel_ele_jet, "ptRel_ele_jet", "relative p_{T} electron - closest jet", 2000, 0, 2000);
	addHistoDataAndMC(type, h_minDeltaR_ele_Jet, "minDeltaR_ele_Jet", "min. #Delta R(electron, jet)", 200, 0, 2);
	addHistoDataAndMC(type, h_mtlep, "mtlep", "reconstructed leptonic top mass", 1000, 0, 1000);
	addHistoDataAndMC(type, h_mthad, "mthad", "reconstructed hadronic top mass", 1000, 0, 1000);
	addHistoDataAndMC(type, h_tlep_pt, "thad_pt", "reconstructed leptonic top p_{T}", 2000, 0, 2000);
	addHistoDataAndMC(type, h_thad_pt, "tlep_pt", "reconstructed hadronic top p_{T}", 2000, 0, 2000);
	addHistoDataAndMC(type, h_angle_b_ele, "angle_b_ele", "angle between b-jet and electron", 400, 0, 4);
	addHistoDataAndMC(type, h_ptratio, "ptratio", "Pt ratio (top/W)", 600, -6, 6);
	addHistoDataAndMC(type, h_pttbar, "pttbar", "Pttbar/HT", 200, 0, 2);
	addHistoDataAndMC(type, h_htsystem, "htsystem", "HTsystem", 240, 0, 1.2);
	addHistoDataAndMC(type, h_angle_b_e_vs_mttbar, "_angle_b_e_vs_mttbar",
			"angle between b-jet and electron vs reconstructed mttbar", 400, 0, 4, 1000, 0, 5000);
	addHistoDataAndMC(type, h_ptratio_vs_mttbar, "ptratio_vs_mttbar", "Pt ratio (top/W) vs reconstructed mttbar", 600, -6, 6,
			1000, 0, 5000);
	addHistoDataAndMC(type, h_pttbar_vs_mttbar, "pttbar_vs_mttbar", "Pttbar/HT vs reconstructed mttbar", 200, 0, 2, 1000, 0, 5000);
	addHistoDataAndMC(type, h_htsystem_vs_mttbar, "htsystem_vs_mttbar", "HTsystem vs reconstructed mttbar", 240, 0, 1.2, 1000, 0,
			5000);
	addHistoDataAndMC(type, h_Chi2Leptonic, "Chi2Leptonic", "#Chi^{2}_{lep}", 1000, 0, 100);
	addHistoDataAndMC(type, h_Chi2Hadronic, "Chi2Hadronic", "#Chi^{2}_{had}", 1000, 0, 100);
	addHistoDataAndMC(type, h_Chi2Global, "Chi2Global", "#Chi^{2}_{global}", 1000, 0, 100);
	addHistoDataAndMC(type, h_Chi2Total, "Chi2Total", "#Chi^{2}_{total}", 3000, 0, 300);
	addHistoDataAndMC(type, h_ptRel_vs_deltaRmin, "kptRel_vs_deltaRmin", "ptrel vs min. #Delta R", 600, 0, 6, 2000, 0, 2000);
	addHistoDataAndMC(type, h_mWhad_vs_mthad, "mWhad_vs_mthad", "M(W_{had}) vs M(t_{had})", 1000, 0, 1000, 500, 0, 500);
	addHistoDataAndMC(type, k2D_mWlep_vs_mtlep, "mWlep_vs_mtlep", "M(W_{lep}) vs M(t_{lep})", 1000, 0, 1000, 500, 0, 500);

}

/**
 * Book histograms which are filled using parton matched objects
 * @param type MCType
 */
void ana::bookExoticTopMatchedHistograms(ushort type) {
	addHistoDataAndMC(type, h_mttbar_matched, "mttbar_matched", "reconstructed M_{ttbar} (matched)", 10000, 0, 10000);
	addHistoDataAndMC(type, h_mtlep_matched, "mtlep_matched", "reconstructed leptonic top mass (matched)", 1000, 0, 1000);
	addHistoDataAndMC(type, h_mthad_matched, "mthad_matched", "reconstructed leptonic top mass (matched)", 1000, 0, 1000);
	addHistoDataAndMC(type, h_thad_pt_matched, "thad_pt_matched", "reconstructed hadronic top p_{T} (matched)", 2000, 0, 2000);
	addHistoDataAndMC(type, h_tlep_pt_matched, "tlep_pt_matched", "reconstructed leptonic top p_{T} (matched)", 2000, 0, 2000);

	addHistoDataAndMC(type, h_mWhad_matched, "mWhad_matched", "reconstructed hadronic W mass(matched)", 500, 0, 500);
	addHistoDataAndMC(type, h_mWlep_matched, "mWlep_matched", "reconstructed leptonic W mass(matched)", 500, 0, 500);

	addHistoDataAndMC(type, h_angle_b_ele_matched, "angle_b_ele_matched", "angle between b-jet and electron(matched)", 400, 0, 4);
	addHistoDataAndMC(type, h_ptratio_matched, "ptratio_matched", "Pt ratio (top/W)(matched)", 600, -6, 6);
	addHistoDataAndMC(type, h_ptratio2_matched, "ptratio2_matched", "Pt ratio (b/W) (matched)", 600, -6, 6);
	addHistoDataAndMC(type, h_pttbar_matched, "pttbar_matched", "Pttbar/HT(matched)", 200, 0, 2);
	addHistoDataAndMC(type, h_htsystem_matched, "htsystem_matched", "HTsystem (matched)", 240, 0, 1.2);

	addHistoDataAndMC(type, h_Chi2Leptonic_matched, "Chi2Leptonic_matched", "#Chi^{2}_{lep} (matched)", 1000, 0, 100);
	addHistoDataAndMC(type, h_Chi2Hadronic_matched, "Chi2Hadronic_matched", "#Chi^{2}_{had} (matched)", 1000, 0, 100);
	addHistoDataAndMC(type, h_Chi2Global_matched, "Chi2Global_matched", "#Chi^{2}_{global} (matched)", 1000, 0, 100);
	addHistoDataAndMC(type, h_Chi2Total_matched, "Chi2Total_matched", "#Chi^{2}_{total} (matched)", 3000, 0, 300);
}

void ana::bookExoticTopBtaggedHistograms(ushort type) {
	addHistoDataAndMC(type, h_mttbar_btag_fake, "mttbar_btag_fake", "reconstructed M_{ttbar} (fake btag)", 10000, 0, 10000);
	addHistoDataAndMC(type, h_mWhad_btag_fake, "mWhad_btag_fake", "reconstructed hadronic W mass (fake btag)", 500, 0, 500);
	addHistoDataAndMC(type, h_mWlep_btag_fake, "mWlep_btag_fake", "reconstructed leptonic W mass (fake btag)", 500, 0, 500);
	addHistoDataAndMC(type, h_mtlep_btag_fake, "mtlep_btag_fake", "reconstructed leptonic top mass (fake btag)", 1000, 0, 1000);
	addHistoDataAndMC(type, h_mthad_btag_fake, "mthad_btag_fake", "reconstructed hadronic top mass (fake btag)", 1000, 0, 1000);
	addHistoDataAndMC(type, k2D_mWlep_vs_mtlep_btag_fake, "mWlep_vs_mtlep_btag_fake", "M(W_{lep}) vs M(t_{lep}) (fake btag)",
			1000, 0, 1000, 500, 0, 500);
	addHistoDataAndMC(type, h_mWhad_vs_mthad_btag_fake, "mWhad_vs_mthad_btag_fake", "M(W_{had}) vs M(t_{had}) (fake btag)", 1000,
			0, 1000, 500, 0, 500);
}

void ana::bookExoticTopMCHistograms(ushort type) {

	addHistoDataAndMC(type, h_mttbar_mc_smeared, "mttbar_mc_smeared", "reconstructed M_{ttbar} (smeared mc)", 10000, 0, 10000);
	addHistoDataAndMC(type, h_mttbar_diff_reco_and_mc, "mttbar_diff_reco_mc", "reconstructed - generated M_{ttbar}", 1000, -500,
			500);
	addHistoDataAndMC(type, h_neutrino_pz_mc, "neutrino_pz_mc", "generated p_z(#nu)", 50, -500, 500);
	addHistoDataAndMC(type, h_mttbar_mc, "mttbar_mc", "generated M_{ttbar}", 10000, 0, 10000);
	//TODO: decide if to keep
	addHistoDataAndMC(type, h_mZprime_mc, "mZprime_mc", "generated M_{Z'}", 10000, 0, 10000);

	addHistoDataAndMC(type, h_mWlep_mc, "mWlep_mc", "generated leptonic W mass", 500, 0, 500);
	addHistoDataAndMC(type, h_mWhad_mc, "mWhad_mc", "generated hadronic W mass", 500, 0, 500);

	addHistoDataAndMC(type, h_mtlep_mc, "mtlep_mc", "generated leptonic top mass", 1000, 0, 1000);
	addHistoDataAndMC(type, h_mthad_mc, "mthad_mc", "generated hadronic top mass", 1000, 0, 1000);
	addHistoDataAndMC(type, h_tlep_pt_mc, "tlep_pt_mc", "generated leptonic top p_{T}", 2000, 0, 2000);
	addHistoDataAndMC(type, h_thad_pt_mc, "thad_pt_mc", "generated hadronic top p_{T}", 2000, 0, 2000);

	addHistoDataAndMC(type, h_angle_b_ele_mc, "angle_b_ele_mc", "angle between b-jet and electron (MC)", 400, 0, 4);
	addHistoDataAndMC(type, h_ptratio_mc, "ptratio_mc", "Pt ratio (top/W) (MC)", 600, -6, 6);
	addHistoDataAndMC(type, h_ptratio2_mc, "ptratio2_mc", "Pt ratio (b/W) (MC)", 600, -6, 6);
	addHistoDataAndMC(type, h_pttbar_mc, "pttbar_mc", "Pttbar/HT (MC)", 200, 0, 2);
	addHistoDataAndMC(type, h_htsystem_mc, "htsystem_mc", "HTsystem (MC)", 240, 0, 1.2);

	addHistoDataAndMC(type, h_Chi2Leptonic_mc, "Chi2Leptonic_mc", "#Chi^{2}_{lep} (mc)", 1000, 0, 100);
	addHistoDataAndMC(type, h_Chi2Hadronic_mc, "Chi2Hadronic_mc", "#Chi^{2}_{had} (mc)", 1000, 0, 100);
	addHistoDataAndMC(type, h_Chi2Global_mc, "Chi2Global_mc", "#Chi^{2}_{global} (mc)", 1000, 0, 100);
	addHistoDataAndMC(type, h_Chi2Total_mc, "Chi2Total_mc", "#Chi^{2}_{total} (mc)", 3000, 0, 300);
}

/**
 * Basic Kinematic Plots (Aug 09)
 * @param type MC type (or data)
 */
void ana::bookBasicHistograms(ushort type) {
	//-------------
	//  Electrons
	//-------------
	addHistoDataAndMC(type, h_nele, "nele", "N(e) (no cut)", 6, 0, 6);
	addHistoDataAndMC(type, h_ele_ET_all, "ele_ET", "E_{T}(all e) (nocut)", 50, 0, 100);
	addHistoDataAndMC(type, kele_ET_1, "ele1_ET", "E_{T}(leading e: no cut)", 50, 0, 100);
	addHistoDataAndMC(type, kele_ET_2, "ele2_ET", "E_{T}(2^{nd} leading e: no cut)", 50, 0, 100);
	addHistoDataAndMC(type, kele_ET_3, "ele3_ET", "E_{T}(3^{rd} leading e: no cut)", 50, 0, 100);

	addHistoDataAndMC(type, kele_eta_all, "ele_eta", "#eta(all e) (nocut)", 50, -2.5, 2.5);
	addHistoDataAndMC(type, kele_eta_1, "ele1_eta", "#eta(leading e: no cut)", 50, -2.5, 2.5);
	addHistoDataAndMC(type, kele_eta_2, "ele2_eta", "#eta(2^{nd} leading e: no cut)", 50, -2.5, 2.5);
	addHistoDataAndMC(type, kele_eta_3, "ele3_eta", "#eta(3^{rd} leading e: no cut)", 50, -2.5, 2.5);

	addHistoDataAndMC(type, kele_phi_all, "ele_phi", "#phi(all e) (nocut)", 50, -2.5, 2.5);
	addHistoDataAndMC(type, kele_phi_1, "ele1_phi", "#phi(leading e: no cut)", 50, -2.5, 2.5);
	addHistoDataAndMC(type, kele_phi_2, "ele2_phi", "#phi(2^{nd} leading e: no cut)", 50, -2.5, 2.5);
	addHistoDataAndMC(type, kele_phi_3, "ele3_phi", "#phi(3^{rd} leading e: no cut)", 50, -2.5, 2.5);

	addHistoDataAndMC(type, kele_iso_all, "ele_iso", "RelIso(all e) (no cut", 1000, 0, 2);
	addHistoDataAndMC(type, kele_iso_1, "ele1_iso", "RelIso(leading e: no cut", 1000, 0, 2);
	addHistoDataAndMC(type, kele_iso_2, "ele2_iso", "RelIso(2^{nd} leading e: no cut", 1000, 0, 2);
	addHistoDataAndMC(type, kele_iso_3, "ele3_iso", "RelIso(3^{rd} leading e: no cut", 1000, 0, 2);

	addHistoDataAndMC(type, knele_cuts, "nele_cuts", "N(ele) pass each cut", 8, 0, 8);//FIXME: change nbins to kNumCuts
	fasthist_[type][knele_cuts]->GetXaxis()->SetBinLabel(1, "allEle");
	fasthist_[type][knele_cuts]->GetXaxis()->SetBinLabel(2, "ET/Eta");
	fasthist_[type][knele_cuts]->GetXaxis()->SetBinLabel(3, "d0");
	fasthist_[type][knele_cuts]->GetXaxis()->SetBinLabel(4, "eidLoose");
	fasthist_[type][knele_cuts]->GetXaxis()->SetBinLabel(5, "eidTight");
	fasthist_[type][knele_cuts]->GetXaxis()->SetBinLabel(6, "eidRobustLoose");
	fasthist_[type][knele_cuts]->GetXaxis()->SetBinLabel(7, "eidRobustTight");
	fasthist_[type][knele_cuts]->SetLabelSize(0.04, "x");
	fasthist_[type][knele_cuts]->SetOption("hist text0");
	fasthist_[type][knele_cuts]->SetMarkerSize(1.5);
	fasthist_[type][knele_cuts]->SetStats(0);

	//-----------
	//  J E T S
	//-----------
	// 14-8-09: jets after jet-cleaning
	addHistoDataAndMC(type, knjets, "njet", "jet multiplicity (ET/#eta)", 10, 0, 10);
	addHistoDataAndMC(type, kjet_pt_all, "jet_PT", "P_{T}(all jets) (pass P_{T}/#eta)", 100, 0, 200);
	addHistoDataAndMC(type, kjet_pt_1, "jet1_PT", "P_{T}(leading jet) (pass P_{T}/#eta)", 100, 0, 200);
	addHistoDataAndMC(type, kjet_pt_2, "jet2_PT", "P_{T}(2^{nd} leading jet) (pass P_{T}/#eta)", 100, 0, 200);
	addHistoDataAndMC(type, kjet_pt_3, "jet3_PT", "P_{T}(3^{rd} leading jet) (pass P_{T}/#eta)", 100, 0, 200);
	addHistoDataAndMC(type, kjet_pt_4, "jet4_PT", "P_{T}(4^{th} leading jet) (pass P_{T}/#eta)", 100, 0, 200);

	addHistoDataAndMC(type, kjet_eta_all, "jet_eta", "#eta(all jets) (pass P_{T}/#eta)", 50, -2.5, 2.5);
	addHistoDataAndMC(type, kjet_eta_1, "jet1_eta", "#eta(leading jet) (pass P_{T}/#eta)", 50, -2.5, 2.5);
	addHistoDataAndMC(type, kjet_eta_2, "jet2_eta", "#eta(2^{nd} leading jet) (pass P_{T}/#eta)", 50, -2.5, 2.5);
	addHistoDataAndMC(type, kjet_eta_3, "jet3_eta", "#eta(3^{rd} leading jet) (pass P_{T}/#eta)", 50, -2.5, 2.5);
	addHistoDataAndMC(type, kjet_eta_4, "jet4_eta", "#eta(4^{th} leading jet) (pass P_{T}/#eta)", 50, -2.5, 2.5);

	addHistoDataAndMC(type, kjet_phi_all, "jet_phi", "#phi(all jets) (pass P_{T}/#eta)", 50, -3.2, 3.2);
	addHistoDataAndMC(type, kjet_phi_1, "jet1_phi", "#phi(leading jet) (pass P_{T}/#eta)", 50, -3.2, 3.2);
	addHistoDataAndMC(type, kjet_phi_2, "jet2_phi", "#phi(2^{nd} leading jet) (pass P_{T}/#eta)", 50, -3.2, 3.2);
	addHistoDataAndMC(type, h_jet_phi_3, "jet3_phi", "#phi(3^{rd} leading jet) (pass P_{T}/#eta)", 50, -3.2, 3.2);
	addHistoDataAndMC(type, h_jet_phi_4, "jet4_phi", "#phi(4^{th} leading jet) (pass P_{T}/#eta)", 50, -3.2, 3.2);

	//-----
	// MET
	//-----
	addHistoDataAndMC(type, h_metAlone, "metAlone", "Missing ET (no cut)", 200, 0, 200);
	addHistoDataAndMC(type, h_metAlone_phi, "metAlone_phi", "MET #phi (no cut)", 50, -3.2, 3.2);
}

void ana::bookDRemuHistograms(ushort type) {
	addHistoDataAndMC(type, h_DRemu_selE_GoodMu, "DRemu_selE_GoodMu", "#DeltaR(e,#mu) (selE,GoodMu)", 60, 0, 6);
	addHistoDataAndMC(type, h_DRemu_selE_GoodMu_pass, "DRemu_selE_GoodMu_pass", "#DeltaR(e,#mu) (selE,GoodMu) (passAllCut)", 60,
			0, 6);
}

void ana::bookExploreHistograms(ushort type) {
	addHistoDataAndMC(type, h_exp_ele_et, "ele_et", "ele E_{T}", 50, 20, 120);
	addHistoDataAndMC(type, h_exp_ele_eta, "ele_eta", "ele #eta", 50, -2.5, 2.5);
	addHistoDataAndMC(type, h_exp_j0_pt, "j0_pt", "leading jet p_{T}", 100, 0, 200);
	addHistoDataAndMC(type, h_exp_j1_pt, "j1_pt", "2nd leading jet p_{T}", 100, 0, 200);

	addHistoDataAndMC(type, h_exp_DRej, "DRej", "#DeltaR(e,j0)", 60, 0, 6);
	addHistoDataAndMC(type, h_exp_DPhiej, "DPhiej", "#Delta#Phi(e,j0)", 64, -3.2, 3.2);
	addHistoDataAndMC(type, h_exp_DRjj, "DRjj", "#DeltaR(j0,j1)", 60, 0, 6);
	addHistoDataAndMC(type, h_exp_DPhijj, "DPhijj", "#Delta#Phi(j0,j1)", 64, -3.2, 3.2);

	addHistoDataAndMC(type, h_exp_met_v_eeta, "met_v_eeta", "met v #eta(e)", 50, 0, 100, 50, -2.5, 2.5);
}

void ana::bookZvetoHistograms(ushort type, TDirectory *parent) {
	TDirectory *dir_zveto = parent->mkdir("z_veto", "Z veto");
	dir_zveto->cd();
	if (type == kZjets) {//This only for Z->ee MC
		addHistoDataAndMC(type, h_nGenBasicEle_Zee_allj, "nGenBasicEle_Zee_allj", "N(GenEle) (Z->ee, allj) pass ET,eta cut", 6,
				0, 6);
		fasthist_[type][h_nGenBasicEle_Zee_allj]->SetDrawOption("hist text0");
		fasthist_[type][h_nGenBasicEle_Zee_allj]->SetMarkerSize(3);

		addHistoDataAndMC(type, h_Zee_eta, "Zee_eta", "#eta of gen-ele from Z", 100, -5, 5);
		addHistoDataAndMC(type, h_Zee_pt, "Zee_pt", "pt of gen-ele from Z", 100, 0, 100);
		addHistoDataAndMC(type, h_Z_photon_eta, "Z_photon_eta", "#eta of reco-photons in Z events", 100, -4, 4);
		addHistoDataAndMC(type, h_Z_photon_et, "Z_photon_et", "ET of reco-photons in Z->ee events", 100, 0, 100);

		addHistoDataAndMC(type, h_Zee_photon_eta, "Zee_photon_eta", "#eta of reco-photons in Z->ee events", 100, -4, 4);
		addHistoDataAndMC(type, h_Zee_photon_et, "Zee_photon_et", "ET of reco-photons in Z->ee events", 100, 0, 100);

		addHistoDataAndMC(type, h_Zee_photon_eteta_2D, "Zee_photon_eteta_2D", "ET vs eta, of reco-photons in Z->ee events", 100,
				0, 100, 100, -4, 4);

		addHistoDataAndMC(type, h_Z_Nphotons, "Z_Nphotons", "N(reco-photons) in Z events (e,mu,tau)", 6, 0, 6);
		addHistoDataAndMC(type, h_Zee_Nphotons, "Zee_Nphotons", "N(reco-photons) in Z->ee events", 6, 0, 6);
	}
	TDirectory *dir_z1 = parent->mkdir("NES", "Normal Event Selection");
	dir_z1->cd();
	addHistoDataAndMC(type, h_mass_diele, "mass_diele", "M(e,e) (sel,rloose)", 100, 0, 180);
	addHistoDataAndMC(type, h_mass_diele_new, "mass_diele_new", "M(e,e) (sel,rloose) NEW", 100, 0, 180);
	TDirectory *dir_z2 = parent->mkdir("AES", "Anti Event Selection");
	dir_z2->cd();
	addHistoDataAndMC(type, h_mass_diele_lowMet_1j, "mass_diele_lowMet_1j", "M(e,e) (1 rtight, 1 rloose) lowMET 1j", 100, 0, 180);
	addHistoDataAndMC(type, h_mass_ephoton_lowMet_1j, "mass_ephoton_lowMet_1j",
			"M(e,#gamma) (1 rtight, 1 #gamma(nocut)) lowMET 1j", 100, 0, 180);
	addHistoDataAndMC(type, h_Nele_lowMet_1j, "Nele_lowMet_1j", "N(reco-ele) lowMET 1j", 6, 0, 6);
	addHistoDataAndMC(type, h_Nphoton_lowMet_1j, "Nphoton_lowMet_1j", "N(reco-photons) lowMET 1j", 6, 0, 6);
	addHistoDataAndMC(type, h_photon_eta_lowMet_1j, "photon_eta_lowMet_1j", "#eta(reco-photons) lowMET 1j", 100, -4, 4);
	addHistoDataAndMC(type, h_photon_et_lowMet_1j, "photon_et_lowMet_1j", "E_{T}(reco-photons) lowMET 1j", 100, 0, 100);
	addHistoDataAndMC(type, h_photon1_eta_lowMet_1j, "photon1_eta_lowMet_1j", "#eta(leading reco-photon) lowMET 1j", 100, -4, 4);
	addHistoDataAndMC(type, h_photon1_et_lowMet_1j, "photon1_et_lowMet_1j", "E_{T}(leading reco-photon) lowMET 1j", 100, 0, 100);
	parent->cd();
}

void ana::bookBtagHistograms(ushort type, TDirectory *parent) {
	TDirectory *btag_dir = parent->mkdir("btags", "btag algorithms");
	btag_dir->cd();
	addHistoDataAndMC(type, h_btag_TC_highEff_b, "btag_TC_highEff_b", "btag track counting high efficiency (b-quarks)", 100,
			-110, 110);
	addHistoDataAndMC(type, h_btag_TC_highEff_c, "btag_TC_highEff_c", "btag track counting high efficiency (c-quarks)", 100,
			-110, 110);
	addHistoDataAndMC(type, h_btag_TC_highEff_g, "btag_TC_highEff_g", "btag track counting high efficiency (gluons)", 100, -110,
			110);
	addHistoDataAndMC(type, h_btag_TC_highEff_uds, "btag_TC_highEff_uds", "btag track counting high efficiency (uds-quarks)",
			100, -110, 110);

	addHistoDataAndMC(type, h_btag_TC_highPur_b, "btag_TC_highPur_b", "btag track counting high purity (b-quarks)", 80, -80, 80);
	addHistoDataAndMC(type, h_btag_TC_highPur_c, "btag_TC_highPur_c", "btag track counting high purity (c-quarks)", 80, -80, 80);
	addHistoDataAndMC(type, h_btag_TC_highPur_g, "btag_TC_highPur_g", "btag track counting high purity (gluons)", 80, -80, 80);
	addHistoDataAndMC(type, h_btag_TC_highPur_uds, "btag_TC_highPur_uds", "btag track counting high purity (uds-quarks)", 80,
			-80, 80);

	addHistoDataAndMC(type, h_btag_JetBProb_b, "btag_JetBProb_b", "btag JetB Propability (b-quarks)", 100, 0, 12);
	addHistoDataAndMC(type, h_btag_JetBProb_c, "btag_JetBProb_c", "btag JetB Propability (c-quarks)", 100, 0, 12);
	addHistoDataAndMC(type, h_btag_JetBProb_g, "btag_JetBProb_g", "btag JetB Propability (gluons)", 100, 0, 12);
	addHistoDataAndMC(type, h_btag_JetBProb_uds, "btag_JetBProb_uds", "btag JetB Propability (uds-quarks)", 100, 0, 12);

	addHistoDataAndMC(type, h_btag_JetProb_b, "btag_JetProb_b", "btag Jet Propability (b-quarks)", 30, 0, 3.6);
	addHistoDataAndMC(type, h_btag_JetProb_c, "btag_JetProb_c", "btag Jet Propability (c-quarks)", 30, 0, 3.6);
	addHistoDataAndMC(type, h_btag_JetProb_g, "btag_JetProb_g", "btag Jet Propability (gluons)", 30, 0, 3.6);
	addHistoDataAndMC(type, h_btag_JetProb_uds, "btag_JetProb_uds", "btag Jet Propability (uds-quarks)", 30, 0, 3.6);

	addHistoDataAndMC(type, h_btag_secondaryVertex_b, "btag_secondaryVertex_b", "btag secondary vertex (b-quarks)", 80, 0, 6.4);
	addHistoDataAndMC(type, h_btag_secondaryVertex_c, "btag_secondaryVertex_c", "btag secondary vertex (c-quarks)", 80, 0, 6.4);
	addHistoDataAndMC(type, h_btag_secondaryVertex_g, "btag_secondaryVertex_g", "btag secondary vertex (gluons)", 80, 0, 6.4);
	addHistoDataAndMC(type, h_btag_secondaryVertex_uds, "btag_secondaryVertex_uds", "btag secondary vertex (uds-quarks)", 30, 0,
			6.2);

	//FIXME: put proper bins
	addHistoDataAndMC(type, h_btag_softEle_b, "btag_softEle_b", "btag soft electron (b-quarks)", 100, 0, 12);
	addHistoDataAndMC(type, h_btag_softEle_c, "btag_softEle_c", "btag soft electron (c-quarks)", 100, 0, 12);
	addHistoDataAndMC(type, h_btag_softEle_g, "btag_softEle_g", "btag soft electron (gluons)", 100, 0, 12);
	addHistoDataAndMC(type, h_btag_softEle_uds, "btag_softEle_uds", "btag soft electron (uds-quarks)", 30, 0, 3.6);
	addHistoDataAndMC(type, h_btag_softMuon_b, "btag_softMuon_b", "btag soft muon (b-quarks)", 100, 0, 12);
	addHistoDataAndMC(type, h_btag_softMuon_c, "btag_softMuon_c", "btag soft muon (c-quarks)", 100, 0, 12);
	addHistoDataAndMC(type, h_btag_softMuon_g, "btag_softMuon_g", "btag soft muon (gluons)", 100, 0, 12);
	addHistoDataAndMC(type, h_btag_softMuon_uds, "btag_softMuon_uds", "btag soft muon (uds-quarks)", 30, 0, 3.6);
	addHistoDataAndMC(type, h_btag_softMuonNoIP_b, "btag_softMuonNoIP_b", "btag soft muon (b-quarks) no IP", 100, 0, 12);
	addHistoDataAndMC(type, h_btag_softMuonNoIP_c, "btag_softMuonNoIP_c", "btag soft muon (c-quarks) no IP", 100, 0, 12);
	addHistoDataAndMC(type, h_btag_softMuonNoIP_g, "btag_softMuonNoIP_g", "btag soft muon (gluons) no IP", 100, 0, 12);
	addHistoDataAndMC(type, h_btag_softMuonNoIP_uds, "btag_softMuonNoIP_uds", "btag soft muon (uds-quarks) no IP", 30, 0, 3.6);

	parent->cd();
}

/**
 * Reconstruct the invariant mass of the ttbar system
 * @param jets reconstructed jets
 * @param electron one electron candidate
 * @param met missing transverse energy
 */
void ana::reco_Mttbar(const vector<TLorentzVector>& jets, const TLorentzVector& electron, const TLorentzVector& met) {
	TLorentzVector blepjet, bhadjet, q1, q2, Wlep1, Wlep2, tlep1, tlep2, thad, Whad, tlep, neutrino;
	pair<TLorentzVector, TLorentzVector> neutrinos = reconstruct_neutrinos(electron, met);
	Wlep1 = TLorentzVector(neutrinos.first + electron);
	Wlep2 = TLorentzVector(neutrinos.second + electron);
	ushort neutrinoID(0);
	short bhadjetID(-1), blepjetID(-1), q1ID(-1), q2ID(-1);
	ushort njets = number_of_jets_to_use_for_reco;

	double chi2Tot = 999999;
	double ht = GetHT(jets, 8);
	for (ushort blep_id = 0; blep_id < njets; blep_id++) {//get leptonic b-jet
		for (ushort bhad_id = 0; bhad_id < njets; bhad_id++) {//get hadronic b-jet
			if (bhad_id == blep_id)
				continue;
			for (ushort quark1_id = 0; quark1_id < njets; quark1_id++) {//get quark from W decay
				if (quark1_id == blep_id || quark1_id == bhad_id)
					continue;
				for (ushort quark2_id = 0; quark2_id < njets; quark2_id++) {//get quark' from W decay
					if (quark2_id == quark1_id || quark2_id == bhad_id || quark2_id == blep_id)
						continue;
					blepjet = jets[blep_id];
					bhadjet = jets[bhad_id];
					q1 = jets[quark1_id];
					q2 = jets[quark2_id];
					tlep1 = TLorentzVector(blepjet + Wlep1);
					tlep2 = TLorentzVector(blepjet + Wlep2);
					Whad = TLorentzVector(q1 + q2);
					thad = TLorentzVector(bhadjet + Whad);

					double ptratio = TMath::Log(thad.Pt() / Whad.Pt());
					double angle = blepjet.Angle(electron.Vect());
					double chi2lep_1 = GetChi2Leptonic(tlep1.M(), angle);
					double chi2lep_2 = GetChi2Leptonic(tlep2.M(), angle);

					double chi2had = GetChi2Hadronic(Whad.M(), thad.M(), ptratio);

					TLorentzVector tmp1 = TLorentzVector(tlep1 + thad);
					TLorentzVector tmp2 = TLorentzVector(tlep2 + thad);

					double htsystem = (blepjet.Pt() + bhadjet.Pt() + q1.Pt() + q2.Pt()) / ht;

					double chi2global_1 = GetChi2Global(tmp1.Pt() / ht, htsystem);
					double chi2global_2 = GetChi2Global(tmp2.Pt() / ht, htsystem);

					double chi2total_1 = chi2lep_1 + chi2had + chi2global_1;
					double chi2total_2 = chi2lep_2 + chi2had + chi2global_2;
					if (chi2total_1 < chi2total_2 && chi2total_1 < chi2Tot) {
						chi2Tot = chi2total_1;
						neutrinoID = 1;
					} else if (chi2total_1 > chi2total_2 && chi2total_2 < chi2Tot) {
						chi2Tot = chi2total_2;
						neutrinoID = 2;
					}
					bhadjetID = bhad_id;
					blepjetID = blep_id;
					q1ID = quark1_id;
					q2ID = quark2_id;
				}
			}
		}
	}
	//just in case nothing was filled:
	if (neutrinoID == 0 || bhadjetID < 0 || blepjetID < 0 || q1ID < 0 || q2ID < 0)
		return;
	//create best combination:
	TLorentzVector Wlep;
	if (neutrinoID == 1) {
		Wlep = TLorentzVector(neutrinos.first + electron);
		neutrino = neutrinos.first;
		if (Wlep.M() <= 0) {
			Wlep = TLorentzVector(neutrinos.second + electron);
			neutrino = neutrinos.second;
		}
		if (Wlep.M() <= 0)
			cout << "WARNING: No valid neutrino solution found or wrong electron used to reconstruct W (W mass is negative)"
					<< Wlep.M() << endl;
	} else {
		Wlep = TLorentzVector(neutrinos.second + electron);
		neutrino = neutrinos.second;
		if (Wlep.M() <= 0) {
			Wlep = TLorentzVector(neutrinos.first + electron);
			neutrino = neutrinos.first;
		}
		if (Wlep.M() <= 0)
			cout << "WARNING: No valid neutrino solution found or wrong electron used to reconstruct W (W mass is negative)"
					<< Wlep.M() << endl;
	}
	if (m_debug)
		cout << "Wlep M: " << Wlep.M() << ", P: " << Wlep.P() << endl;
	//	Wlep = TLorentzVector(neutrino + electron);
	blepjet = jets[blepjetID];
	bhadjet = jets[bhadjetID];
	q1 = jets[q1ID];
	q2 = jets[q2ID];
	Whad = TLorentzVector(q1 + q2);
	tlep = TLorentzVector(blepjet + Wlep);
	thad = TLorentzVector(bhadjet + Whad);
	TLorentzVector ttbar = TLorentzVector(thad + tlep);

	double htsystem = (blepjet.Pt() + bhadjet.Pt() + q1.Pt() + q2.Pt()) / ht;
	double minDeltaR = blepjet.DeltaR(electron);
	double costheta = blepjet.Dot(electron) / electron.P() / blepjet.P();
	double sintheta = TMath::Sqrt(1 - costheta * costheta);
	double ptrel = sintheta * electron.P();
	double chi2lep = GetChi2Leptonic(tlep.M(), blepjet.Angle(electron.Vect()));
	double chi2had = GetChi2Hadronic(Whad.M(), thad.M(), TMath::Log(thad.Pt() / Whad.Pt()));
	double chi2global = GetChi2Global(ttbar.Pt() / ht, htsystem);
	double angle = blepjet.Angle(electron.Vect());
	double ptratio = TMath::Log(thad.Pt() / Whad.Pt());
	TLorentzVector tmp1 = TLorentzVector(tlep + thad);
	//fill histograms with best choice
	fillHistoDataAndMC(h_mthad, thad.M(), this_weight);
	fillHistoDataAndMC(h_mtlep, tlep.M(), this_weight);
	fillHistoDataAndMC(h_mWhad, Whad.M(), this_weight);
	fillHistoDataAndMC(h_mWlep, Wlep.M(), this_weight);
	fillHistoDataAndMC(h_mttbar, ttbar.M(), this_weight);
	fillHistoDataAndMC(h_neutrino_pz, neutrino.Pz(), this_weight);
	fillHistoDataAndMC(h_pttbar, ttbar.Pt() / ht, this_weight);
	fillHistoDataAndMC(h_htsystem, htsystem, this_weight);
	fillHistoDataAndMC(h_ptRel_ele_jet, ptrel, this_weight);
	fillHistoDataAndMC(h_minDeltaR_ele_Jet, minDeltaR, this_weight);
	fillHistoDataAndMC(h_ptRel_vs_deltaRmin, minDeltaR, ptrel, this_weight);
	fillHistoDataAndMC(h_Chi2Global, chi2global, this_weight);
	fillHistoDataAndMC(h_Chi2Leptonic, chi2lep, this_weight);
	fillHistoDataAndMC(h_Chi2Hadronic, chi2had, this_weight);
	fillHistoDataAndMC(h_Chi2Total, chi2Tot, this_weight);

	fillHistoDataAndMC(h_pttbar_vs_mttbar, tmp1.Pt() / ht, ttbar.M(), this_weight);
	fillHistoDataAndMC(h_ptratio_vs_mttbar, ptratio, ttbar.M(), this_weight);
	fillHistoDataAndMC(h_htsystem_vs_mttbar, htsystem, ttbar.M(), this_weight);
	fillHistoDataAndMC(h_angle_b_e_vs_mttbar, angle, ttbar.M(), this_weight);

	fillHistoDataAndMC(k2D_mWlep_vs_mtlep, tlep.M(), Wlep.M(), this_weight);
	fillHistoDataAndMC(h_mWhad_vs_mthad, thad.M(), Whad.M(), this_weight);

	if (fastmctype_ == kttbar || fastmctype_ == kttjet || fastmctype_ >= Zprime_M500GeV_W5GeV && fastmctype_
			<= Zprime_M4TeV_W400GeV) {
		//only for ttbar and Zprime
		fillMCTopEventHists(ttbar.M());//fill MC information
		reco_Mttbar_matched(jets, electron, met);
		reco_Mttbar_btagged(jets, electron, met);
	}

}

/**
 * Reconstruct the invariant mass of the ttbar system
 * @param jets reconstructed jets
 * @param electron one electron candidate
 * @param met missing transverse energy
 */
void ana::reco_Mttbar_btagged(const vector<TLorentzVector>& jets, const TLorentzVector& electron, const TLorentzVector& met) {
	TLorentzVector blepjet, bhadjet, q1, q2, neutrino, Wlep1, Wlep2, tlep1, tlep2, thad, Whad, tlep;
	pair<TLorentzVector, TLorentzVector> neutrinos = reconstruct_neutrinos(electron, met);
	Wlep1 = TLorentzVector(neutrinos.first + electron);
	Wlep2 = TLorentzVector(neutrinos.second + electron);
	ushort neutrinoID(0);
	short bhadjetID(-1), blepjetID(-1), q1ID(-1), q2ID(-1);
	ushort njets = number_of_jets_to_use_for_reco;

	double chi2Tot = 999999;
	double ht = GetHT(jets, 8);
	//	for (ushort a = 0; a < njets; a++) {//get leptonic b-jet
	//		for (ushort b = 0; b < njets; b++) {//get hadronic b-jet
	//			if (b == a)
	//				continue;
	//fix b jets:
	vector<TLorentzVector> MCEvent = GetMCTopEvent();
	if (MCEvent.size() == 0)
		return;
	//	vector<TLorentzVector> tmp;
	//	tmp.push_back(MCEvent[kblep]);
	short a = match(MCEvent[kblep], jets);
	//	for (ushort x = 0; x < njets; ++x) {
	//		short id = match(jets[x], tmp);
	//		if (id >= 0)
	//			a = id;
	//	}
	if (a < 0)
		return;
	blepjet = jets[a];
	//	tmp.clear();
	//	tmp.push_back(MCEvent[kbhad]);
	short b = match(MCEvent[kbhad], jets);
	//	for (ushort x = 0; x < njets; ++x) {
	//		short id = match(jets[x], tmp);
	//		if (id >= 0)
	//			b = id;
	//	}
	if (b < 0)
		return;
	bhadjet = jets[b];
	for (ushort c = 0; c < njets; c++) {//get quark from W decay
		if (c == a || c == b)
			continue;
		for (ushort d = 0; d < njets; d++) {//get quark' from W decay
			if (d == c || d == b || d == a)
				continue;

			q1 = jets[c];
			q2 = jets[d];
			tlep1 = TLorentzVector(blepjet + Wlep1);
			tlep2 = TLorentzVector(blepjet + Wlep2);
			Whad = TLorentzVector(q1 + q2);
			thad = TLorentzVector(bhadjet + Whad);

			double ptratio = TMath::Log(thad.Pt() / Whad.Pt());
			double angle = blepjet.Angle(electron.Vect());

			double chi2lep_1 = GetChi2Leptonic(tlep1.M(), angle);
			double chi2lep_2 = GetChi2Leptonic(tlep2.M(), angle);

			double chi2had = GetChi2Hadronic(Whad.M(), thad.M(), ptratio);

			TLorentzVector tmp1 = TLorentzVector(tlep1 + thad);
			TLorentzVector tmp2 = TLorentzVector(tlep2 + thad);

			double htsystem = (blepjet.Pt() + bhadjet.Pt() + q1.Pt() + q2.Pt()) / ht;

			double chi2global_1 = GetChi2Global(tmp1.Pt() / ht, htsystem);
			double chi2global_2 = GetChi2Global(tmp2.Pt() / ht, htsystem);

			double chi2total_1 = chi2lep_1 + chi2had + chi2global_1;
			double chi2total_2 = chi2lep_2 + chi2had + chi2global_2;
			if (chi2total_1 < chi2total_2 && chi2total_1 < chi2Tot) {
				chi2Tot = chi2total_1;
				neutrinoID = 1;
			} else if (chi2total_1 > chi2total_2 && chi2total_2 < chi2Tot) {
				chi2Tot = chi2total_2;
				neutrinoID = 2;
			}
			bhadjetID = b;
			blepjetID = a;
			q1ID = c;
			q2ID = d;
		}
		//	}
		//}
	}
	//just in case nothing was filled:
	if (neutrinoID == 0 || bhadjetID < 0 || blepjetID < 0 || q1ID < 0 || q2ID < 0)
		return;
	//create best combination:
	TLorentzVector Wlep;
	if (neutrinoID == 1)
		neutrino = neutrinos.first;
	else
		neutrino = neutrinos.second;

	Wlep = TLorentzVector(neutrino + electron);
	blepjet = jets[blepjetID];
	bhadjet = jets[bhadjetID];
	q1 = jets[q1ID];
	q2 = jets[q2ID];
	Whad = TLorentzVector(q1 + q2);
	tlep = TLorentzVector(blepjet + Wlep);
	thad = TLorentzVector(bhadjet + Whad);
	TLorentzVector resonance = TLorentzVector(thad + tlep);

	//fill histograms with best choice
	fillHistoDataAndMC(h_mthad_btag_fake, thad.M(), this_weight);
	fillHistoDataAndMC(h_mtlep_btag_fake, tlep.M(), this_weight);
	fillHistoDataAndMC(h_mWhad_btag_fake, Whad.M(), this_weight);
	fillHistoDataAndMC(h_mWlep_btag_fake, Wlep.M(), this_weight);
	fillHistoDataAndMC(h_mttbar_btag_fake, resonance.M(), this_weight);
	fillHistoDataAndMC(h_mWhad_vs_mthad_btag_fake, thad.M(), Whad.M(), this_weight);
	fillHistoDataAndMC(k2D_mWlep_vs_mtlep_btag_fake, tlep.M(), Wlep.M(), this_weight);
}
/**
 * Reconstruct the invariant mass of the ttbar system using parton matched jets only
 * @param jets reconstructed jets
 * @param electron one electron candidate
 * @param met missing transverse energy
 */
void ana::reco_Mttbar_matched(const vector<TLorentzVector>& jets, const TLorentzVector& electron, const TLorentzVector& met) {
	vector<TLorentzVector> MCEvent = GetMCTopEvent();
	if (MCEvent.size() == 0)
		return;
	vector<TLorentzVector> mcpartons;
	mcpartons.push_back(MCEvent[kq1]);
	mcpartons.push_back(MCEvent[kq2]);
	mcpartons.push_back(MCEvent[kbhad]);
	mcpartons.push_back(MCEvent[kblep]);
	vector<pair<ushort, ushort> > ids;
	short id;

	for (ushort x = 0; x < jets.size(); x++) {
		id = match(jets[x], mcpartons);
		if (id >= 0)
			ids.push_back(pair<ushort, ushort> (x, id));
	}

	//remove duplicates
	if (ids.size() > 4) {
		cout << "WARNING: duplicate matches with MC partons" << endl;
		//		ids.erase(unique( ids.begin(), ids.end()), ids.end());
	}

	if (ids.size() == 4) {
		//fill histos
		TLorentzVector blepjet, bhadjet, q1jet, q2jet;
		for (ushort x = 0; x < ids.size(); x++) {
			ushort mc = ids[x].second;
			if (mc == 0)
				q1jet = jets[ids[x].first];
			else if (mc == 1)
				q2jet = jets[ids[x].first];
			else if (mc == 2)
				bhadjet = jets[ids[x].first];
			else if (mc == 3)
				blepjet = jets[ids[x].first];
		}
		TLorentzVector thad, tlep, Whad, Wlep;
		Whad = TLorentzVector(q1jet + q2jet);
		TLorentzVector neutrino;
		neutrino = reconstruct_neutrino(electron, met);
		Wlep = TLorentzVector(neutrino + electron);
		thad = TLorentzVector(Whad + bhadjet);
		tlep = TLorentzVector(Wlep + blepjet);

		double ht = GetHT(jets, 8);
		double htsystem = (bhadjet.Pt() + blepjet.Pt() + q1jet.Pt() + q2jet.Pt()) / ht;

		TLorentzVector temp = TLorentzVector(tlep + thad);
		TLorentzVector temp2 = TLorentzVector(blepjet + bhadjet + q1jet + q2jet + neutrino + electron);
		double pttbar = temp.Pt() / ht;
		pttbar = temp2.Pt() / ht;

		double chi2lep = GetChi2Leptonic(tlep.M(), blepjet.Angle(electron.Vect()));
		double chi2had = GetChi2Hadronic(Whad.M(), thad.M(), TMath::Log(thad.Pt() / Whad.Pt()));
		double chi2global = GetChi2Global(temp.Pt() / ht, htsystem);
		double chi2Tot = chi2lep + chi2had + chi2global;
		TLorentzVector resonance = TLorentzVector(thad + tlep);

		fillHistoDataAndMC(h_thad_pt_matched, thad.Pt(), this_weight);
		fillHistoDataAndMC(h_mthad_matched, thad.M(), this_weight);
		fillHistoDataAndMC(h_tlep_pt_matched, tlep.Pt(), this_weight);
		fillHistoDataAndMC(h_mtlep_matched, tlep.M(), this_weight);
		fillHistoDataAndMC(h_mWhad_matched, Whad.M(), this_weight);
		fillHistoDataAndMC(h_mWlep_matched, Wlep.M(), this_weight);
		fillHistoDataAndMC(h_angle_b_ele_matched, blepjet.Angle(electron.Vect()), this_weight);
		fillHistoDataAndMC(h_ptratio_matched, TMath::Log(thad.Pt() / Whad.Pt()), this_weight);
		fillHistoDataAndMC(h_ptratio2_matched, TMath::Log(bhadjet.Pt() / Whad.Pt()), this_weight);
		fillHistoDataAndMC(h_htsystem_matched, htsystem, this_weight);
		fillHistoDataAndMC(h_pttbar_matched, pttbar, this_weight);
		fillHistoDataAndMC(h_Chi2Global_matched, chi2global, this_weight);
		fillHistoDataAndMC(h_Chi2Leptonic_matched, chi2lep, this_weight);
		fillHistoDataAndMC(h_Chi2Hadronic_matched, chi2had, this_weight);
		fillHistoDataAndMC(h_Chi2Total_matched, chi2Tot, this_weight);
		fillHistoDataAndMC(h_mttbar_matched, resonance.M(), this_weight);
	}

}

TLorentzVector ana::reconstruct_neutrino(const TLorentzVector& electron, const TLorentzVector& met) {
	TLorentzVector neutrino1, neutrino2, neutrino;
	pair<double, double> pz = compute_neutrino_momentum_z(met, electron);
	double energy1 = TMath::Sqrt(met.Pt() * met.Pt() + pz.first * pz.first);
	double energy2 = TMath::Sqrt(met.Pt() * met.Pt() + pz.second * pz.second);
	neutrino1.SetPxPyPzE(met.Px(), met.Py(), pz.first, energy1);
	neutrino2.SetPxPyPzE(met.Px(), met.Py(), pz.second, energy2);
	//choose right solution
	bool stage1 = fabs(neutrino1.Pz() - electron.Pz()) < fabs(neutrino2.Pz() - electron.Pz());
	stage1 ? neutrino = neutrino1 : neutrino = neutrino2;
	if (fabs(neutrino.Pz()) > 300)
		fabs(neutrino1.Pz()) < fabs(neutrino2.Pz()) ? neutrino = neutrino1 : neutrino = neutrino2;

	return neutrino;
}

pair<TLorentzVector, TLorentzVector> ana::reconstruct_neutrinos(const TLorentzVector& electron, const TLorentzVector& met) {
	TLorentzVector neutrino1, neutrino2;
	pair<double, double> pz = compute_neutrino_momentum_z(met, electron);
	double energy1 = TMath::Sqrt(met.Et() * met.Et() + pz.first * pz.first);
	double energy2 = TMath::Sqrt(met.Et() * met.Et() + pz.second * pz.second);
	neutrino1.SetPxPyPzE(met.Px(), met.Py(), pz.first, energy1);
	neutrino2.SetPxPyPzE(met.Px(), met.Py(), pz.second, energy2);
	return pair<TLorentzVector, TLorentzVector> (neutrino1, neutrino2);
}

pair<unsigned short, unsigned short> ana::reco_hadronic_W(const std::vector<TLorentzVector>& jets, short ommit_blep_id) const {
	double tempmass = 0;
	unsigned short jet1(0), jet2(0);
	for (unsigned short i = 0; i < jets.size(); i++) {
		if (i == ommit_blep_id)
			continue;

		for (unsigned short j = 0; j < jets.size(); j++) {
			if (j == ommit_blep_id)
				continue;
			TLorentzVector tempW = TLorentzVector(jets[i] + jets[j]);
			//if this is closer than last one:
			if (fabs(tempW.M() - WMASS) < fabs(tempmass - WMASS)) {
				tempmass = tempW.M();
				jet1 = i;
				jet2 = j;
			}

		}
	}
	return pair<unsigned short, unsigned short> (jet1, jet2);
}

/**
 * Find the closest(in DeltaR) object from a group to another object
 * @param group group of reconstructed objects
 * @param obj reconstructed object
 * @return ID of closest object from the group
 */
short ana::findClosest(const vector<TLorentzVector>& group, const TLorentzVector& obj) {
	short cId = -1;
	float deltaR = 99999;
	for (unsigned int i = 0; i < group.size(); i++) {
		float dr = calcDeltaR(group[i], obj);
		if (dr < deltaR) {
			deltaR = dr;
			cId = i;
		}
	}
	return cId;
}

/**
 * Used to calculate the neutrino longitudinal momentum (P_z)
 * @param met missing transverse energy
 * @param electron the electron candidate
 * @return the two solutions for the longitudinal momentum
 */
pair<double, double> ana::compute_neutrino_momentum_z(const TLorentzVector& met, const TLorentzVector& electron) {
	double pte, pze, ptnu, pznu1, pznu2, Ee, pxe, pye, pxnu, pynu;
	pznu1 = 0;
	pznu2 = 0;
	pxe = electron.Px();
	pye = electron.Py();
	pte = electron.Pt();
	pze = electron.Pz();
	pxnu = met.Px();
	pynu = met.Py();
	ptnu = met.Et();
	double wmass = smearValueWithGaus(WMASS, WMASSERROR);
	Ee = electron.Energy();
	//	cout << "--------------------" << endl;
	//	cout << "E(e) = " << Ee << ", Px(e) = " << pxe<< ", Py(e) = " << pye<< ", Pt(e) = " << pte << ", Pz(e)= " << pze << endl;
	//	cout << "Px(nu) = " << pxnu<< ", Py(nu) = " << pynu<< ", Pt(nu) = " << met.Et() <<endl;
	//	cout << "E(e) = " << Ee << ", Pt(e) = " << pte << ", Pz(e)= " << pze << endl;
	//	cout << "MET = " << met.Et() << endl;
	//	double k = wmass * wmass * 0.5 + ptnu * pte;//checked// all roots are real
	double k = wmass * wmass * 0.5 + (pxe * pxnu + pye * pynu);//checked// some complex roots
	double l = Ee * Ee - pze * pze;//checked
	double b = pze * k / l;
	double q = (k * k - Ee * Ee * ptnu * ptnu) / l;
	double r = b * b + q;
	//	cout << "k(" << k << ") l(" << l << ") q(" << q << ")" << endl;
	//	cout << "b(" << b << ") r(" << r << ")" << endl;
	if (r >= 0) {//positive root
		pznu1 = b + sqrt(r);
		pznu2 = b - sqrt(r);
	} else if (r < 0) {//complex root
		pznu1 = b;//take only real part
		pznu2 = b;//take only real part
		//		cout << "complex root " << endl;
	}
	return pair<double, double> (pznu1, pznu2);
}

/**
 * Calculate the total number of events for a given MC type, cut and startpoint for Njet bin
 *
 * @param nevent the counter
 * @param cut the cut to sum for
 * @param k_start starting point for MC
 * @param k_end end point for MC
 * @param njbegin start point for NJet bin
 * @return number of events for the given requirements
 */
double ana::getTotalEvents(const double nevent[][5][nmctype], short cut, short k_start, short k_end, short njbegin) const {
	double totalT = 0; //total for a particular mc type
	ushort offset = 0;
	if (k_start > kttjet)
		offset = NUMBER_OF_DECAY_BRANCHES;
	for (short k = k_start + offset; k <= k_end + offset; ++k) { //loop over mc types
		for (short j = njbegin; j < ntjet; ++j) {
			totalT += nevent[cut][j][k];
		} //sum up jet bins
	}
	return totalT;
}

/**
 * Read MC truth to reconstruct the semileptonic (electron) ttbar event:
 * both b-jets, Ws, top quarks as well as the neutrino, electron and both quarks from
 * the hadronically decaying W
 *
 * @return vector of event particles
 */
vector<TLorentzVector> ana::GetMCTopEvent() {
	TLorentzVector blep_mc, thad_mc, tlep_mc, bhad_mc, Whad_mc, Wlep_mc, neutrino_mc, electron_mc;
	vector<TLorentzVector> QsFromW, MCEvent;
	//use float to prevent overflow
	float tlepPID(0), thadPID(0), WlepPID(0), WhadPID(0), blepPID(0), bhadPID(0), electronPID(0);
	//use neutrino from W-decay from top-decay as reference particle
	for (unsigned int x = 0; x < Nmc_doc; x++) {
		bool isneutrino = fabs(mc_doc_id->at(x)) == 12//electron neutrino
				&& fabs(mc_doc_mother_id->at(x)) == 24//W-boson
				&& fabs(mc_doc_grandmother_id->at(x)) == 6;//top quark
		if (isneutrino) {
			tlepPID = mc_doc_grandmother_id->at(x);
			thadPID = -tlepPID;
			WlepPID = mc_doc_mother_id->at(x);
			WhadPID = -WlepPID;
			blepPID = 5 * (WlepPID / 24);
			bhadPID = -blepPID;
			electronPID = -11 * (mc_doc_id->at(x) / 12); //opposite sign of the electron neutrino
			neutrino_mc = TLorentzVector(mc_doc_px->at(x), mc_doc_py->at(x), mc_doc_pz->at(x), mc_doc_energy->at(x));
		}
	}
	if (tlepPID == 0 || thadPID == 0 || WlepPID == 0 || WhadPID == 0 || blepPID == 0 || bhadPID == 0 || electronPID == 0)
		return MCEvent;
	for (unsigned int x = 0; x < Nmc_doc; x++) {
		bool isBlep = mc_doc_id->at(x) == blepPID && mc_doc_grandmother_id->at(x) == tlepPID;
		bool isBhad = mc_doc_id->at(x) == bhadPID && mc_doc_grandmother_id->at(x) == thadPID;
		bool isThad = mc_doc_id->at(x) == thadPID;
		bool isTlep = mc_doc_id->at(x) == tlepPID;
		bool isWhad = mc_doc_id->at(x) == WhadPID;
		bool isWlep = mc_doc_id->at(x) == WlepPID;
		bool isElectron = mc_doc_id->at(x) == electronPID && mc_doc_mother_id->at(x) == WlepPID;
		bool isQfromW = mc_doc_mother_id->at(x) == WhadPID;

		if (isBlep)
			blep_mc = TLorentzVector(mc_doc_px->at(x), mc_doc_py->at(x), mc_doc_pz->at(x), mc_doc_energy->at(x));
		else if (isBhad)
			bhad_mc = TLorentzVector(mc_doc_px->at(x), mc_doc_py->at(x), mc_doc_pz->at(x), mc_doc_energy->at(x));
		else if (isThad) {
			thad_mc = TLorentzVector(mc_doc_px->at(x), mc_doc_py->at(x), mc_doc_pz->at(x), mc_doc_energy->at(x));
		} else if (isTlep) {
			tlep_mc = TLorentzVector(mc_doc_px->at(x), mc_doc_py->at(x), mc_doc_pz->at(x), mc_doc_energy->at(x));
		} else if (isWhad)
			Whad_mc = TLorentzVector(mc_doc_px->at(x), mc_doc_py->at(x), mc_doc_pz->at(x), mc_doc_energy->at(x));
		else if (isWlep)
			Wlep_mc = TLorentzVector(mc_doc_px->at(x), mc_doc_py->at(x), mc_doc_pz->at(x), mc_doc_energy->at(x));
		else if (isElectron)
			electron_mc = TLorentzVector(mc_doc_px->at(x), mc_doc_py->at(x), mc_doc_pz->at(x), mc_doc_energy->at(x));
		else if (isQfromW)
			QsFromW.push_back(TLorentzVector(mc_doc_px->at(x), mc_doc_py->at(x), mc_doc_pz->at(x), mc_doc_energy->at(x)));
	}

	MCEvent.resize(NUMBER_OF_EVENT_PARTONS);
	MCEvent[ktlep] = tlep_mc;
	MCEvent[kthad] = thad_mc;
	MCEvent[kblep] = blep_mc;
	MCEvent[kbhad] = bhad_mc;
	MCEvent[kWhad] = Whad_mc;
	MCEvent[kWlep] = Wlep_mc;
	MCEvent[kq1] = QsFromW[0];
	MCEvent[kq2] = QsFromW[1];
	MCEvent[kelectron] = electron_mc;
	MCEvent[kneutrino] = neutrino_mc;
	if (QsFromW.size() > 2) {
		cout << "WARNING: Too many quarks from W-decay. " << QsFromW.size() << " partons found" << endl;
	}
	return MCEvent;

}

/**
 * Fills histograms with MC truth information
 */
void ana::fillMCTopEventHists(double reco_mttbar) {
	vector<TLorentzVector> MCEvent = GetMCTopEvent();
	if (MCEvent.size() == 0)
		return;
	TLorentzVector ttbar = TLorentzVector(MCEvent[ktlep] + MCEvent[kthad]);

	fasthist_[fastmctype_][h_neutrino_pz_mc]->Fill(MCEvent[kneutrino].Pz());
	fasthist_[fastmctype_][h_mtlep_mc]->Fill(MCEvent[ktlep].M());
	fasthist_[fastmctype_][h_mthad_mc]->Fill(MCEvent[kthad].M());
	fasthist_[fastmctype_][h_mWlep_mc]->Fill(MCEvent[kWlep].M());
	fasthist_[fastmctype_][h_mWhad_mc]->Fill(MCEvent[kWhad].M());
	fasthist_[fastmctype_][h_mttbar_mc]->Fill(ttbar.M());
	fasthist_[fastmctype_][h_mttbar_diff_reco_and_mc]->Fill(reco_mttbar - ttbar.M());
	fasthist_[fastmctype_][h_mttbar_mc_smeared]->Fill(smearValueWithGaus(ttbar.M(), 90)); //smear 90 GeV width
	fasthist_[fastmctype_][h_thad_pt_mc]->Fill(MCEvent[kthad].Pt());
	fasthist_[fastmctype_][h_tlep_pt_mc]->Fill(MCEvent[ktlep].Pt());
	double ptratio = TMath::Log(MCEvent[kthad].Pt() / MCEvent[kWhad].Pt());
	fasthist_[fastmctype_][h_ptratio_mc]->Fill(ptratio);
	ptratio = TMath::Log(MCEvent[kbhad].Pt() / MCEvent[kWhad].Pt());
	fasthist_[fastmctype_][h_ptratio2_mc]->Fill(ptratio);
	if (fastmctype_ >= Zprime_M500GeV_W5GeV && fastmctype_ <= Zprime_M4TeV_W400GeV) {
		for (unsigned int x = 0; x < Nmc_doc; x++) {
			bool zprime = fabs(mc_doc_id->at(x)) == 32;
			if (zprime)
				fasthist_[fastmctype_][h_mZprime_mc]->Fill(mc_doc_mass->at(x));
		}

	}

}

/**
 * Calculate the Chi2 for the leptonic site
 *
 * @param Wmass leptonic W mass
 * @param tmass leptonic top quark mass
 * @param angle angle between
 * @return
 */
double ana::GetChi2Leptonic(double tmass, double angle) {
	const double angle_mc = 1.0972;
	const double angle_sigma = 0.86786;
	const double tmass_mc = 166.134;
	const double tmass_sigma = 52.743;

	double an = TMath::Power(angle - angle_mc, 2) / (2 * angle_sigma * angle_sigma);
	double tm = TMath::Power(tmass - tmass_mc, 2) / (2 * tmass_sigma * tmass_sigma);
	return 1 / TMath::Sqrt(2) * (an + tm);
}

double ana::GetChi2Hadronic(double Wmass, double tmass, double ptratio) {
	const double Wmass_mc = 90.38;
	const double Wmass_sigma = 17.47;
	const double tmass_mc = 182.476;
	const double tmass_sigma = 21.858;
	const double ptratio_mc = 0.1259;
	const double ptratio_sigma = 0.4013;
	double pt = TMath::Power(ptratio - ptratio_mc, 2) / (2 * ptratio_sigma * ptratio_sigma);
	double wm = TMath::Power(Wmass - Wmass_mc, 2) / (2 * Wmass_sigma * Wmass_sigma);
	double tm = TMath::Power(tmass - tmass_mc, 2) / (2 * tmass_sigma * tmass_sigma);
	return 1 / TMath::Sqrt(3) * (wm + pt + tm);

}

/**
 * Scalar sum of PT (HT)
 * @param jets: jet collection
 * @param N: the maximal number of jets which should be used to calculate HT
 * @return HT
 */
double ana::GetHT(const std::vector<TLorentzVector>& jets, ushort N) {
	double HT(0);
	ushort limit = jets.size();
	if (limit > N + 1)//set to maximal number of jets
		limit = N + 1;
	for (ushort a = 0; a < limit; a++) {
		HT += jets[a].Pt();
	}
	return HT;

}

double ana::GetChi2Global(double pttbar, double htsystem) {
	//	const double pttbar_mc = 0.156677;
	const double pttbar_sigma = 0.0677653;
	//not using values from the gaus fit, since HTsystem is almost always == 1;
	const double htsystem_mc = 1.0;
	//taking sigma to be 10% so the mean of the bump is included.
	const double htsystem_sigma = 0.1;
	//not using pttbar_mc since pttbar->0 for higher Z mass
	double pt = TMath::Power(pttbar, 2) / (pttbar_sigma * pttbar_sigma * 2);
	double ht = TMath::Power(htsystem - htsystem_mc, 2) / (2 * htsystem_sigma * htsystem_sigma);
	return 1 / TMath::Sqrt(2) * (pt + ht);

}

double ana::smearValueWithGaus(double value, double gaus_width) {
	return value + gRandom->Gaus(0, gaus_width);
}

//double ana::btag(ushort jetid, ushort btag) {
//	return btag_information[btag][jetid];
//}

void ana::fillBtagHistograms(std::vector<TLorentzVector> jets, std::vector<std::vector<double> > btag_information) {
	for (ushort jetID = 0; jetID < jets.size(); ++jetID) {
		double pdgID = 0;
		if (m_jetAlgo == "Default" || m_jetAlgo == "pfjet")
			pdgID = fabs(btag_information[btag_fake][jetID]);
		else
			pdgID = fabs(MCTruthMatch(jets[jetID].Eta(), jets[jetID].Phi()));

		if (pdgID == 5) {
			//fill B jets
			fillHistoDataAndMC(h_btag_JetBProb_b, btag_information[btag_JetBProb][jetID], this_weight);
			fillHistoDataAndMC(h_btag_JetProb_b, btag_information[btag_JetProb][jetID], this_weight);
			fillHistoDataAndMC(h_btag_secondaryVertex_b, btag_information[btag_secondaryVertex][jetID], this_weight);
			fillHistoDataAndMC(h_btag_TC_highEff_b, btag_information[btag_TC_highEff][jetID], this_weight);
			fillHistoDataAndMC(h_btag_TC_highPur_b, btag_information[btag_TC_highPur][jetID], this_weight);
			fillHistoDataAndMC(h_btag_softEle_b, btag_information[btag_softEle][jetID], this_weight);
			fillHistoDataAndMC(h_btag_softMuon_b, btag_information[btag_softMuon][jetID], this_weight);
			fillHistoDataAndMC(h_btag_softMuonNoIP_b, btag_information[btag_softMuonNoIP][jetID], this_weight);
		} else if (pdgID == 4) {
			fillHistoDataAndMC(h_btag_JetBProb_c, btag_information[btag_JetBProb][jetID], this_weight);
			fillHistoDataAndMC(h_btag_JetProb_c, btag_information[btag_JetProb][jetID], this_weight);
			fillHistoDataAndMC(h_btag_secondaryVertex_c, btag_information[btag_secondaryVertex][jetID], this_weight);
			fillHistoDataAndMC(h_btag_TC_highEff_c, btag_information[btag_TC_highEff][jetID], this_weight);
			fillHistoDataAndMC(h_btag_TC_highPur_c, btag_information[btag_TC_highPur][jetID], this_weight);
			fillHistoDataAndMC(h_btag_softEle_c, btag_information[btag_softEle][jetID], this_weight);
			fillHistoDataAndMC(h_btag_softMuon_c, btag_information[btag_softMuon][jetID], this_weight);
			fillHistoDataAndMC(h_btag_softMuonNoIP_c, btag_information[btag_softMuonNoIP][jetID], this_weight);
		} else if (pdgID >= 1 && pdgID <= 3) {
			//fill uds
			fillHistoDataAndMC(h_btag_JetBProb_uds, btag_information[btag_JetBProb][jetID], this_weight);
			fillHistoDataAndMC(h_btag_JetProb_uds, btag_information[btag_JetProb][jetID], this_weight);
			fillHistoDataAndMC(h_btag_secondaryVertex_uds, btag_information[btag_secondaryVertex][jetID], this_weight);
			fillHistoDataAndMC(h_btag_TC_highEff_uds, btag_information[btag_TC_highEff][jetID], this_weight);
			fillHistoDataAndMC(h_btag_TC_highPur_uds, btag_information[btag_TC_highPur][jetID], this_weight);
			fillHistoDataAndMC(h_btag_softEle_uds, btag_information[btag_softEle][jetID], this_weight);
			fillHistoDataAndMC(h_btag_softMuon_uds, btag_information[btag_softMuon][jetID], this_weight);
			fillHistoDataAndMC(h_btag_softMuonNoIP_uds, btag_information[btag_softMuonNoIP][jetID], this_weight);
		} else if (pdgID == 21) {
			//fill gluons
			fillHistoDataAndMC(h_btag_JetBProb_g, btag_information[btag_JetBProb][jetID], this_weight);
			fillHistoDataAndMC(h_btag_JetProb_g, btag_information[btag_JetProb][jetID], this_weight);
			fillHistoDataAndMC(h_btag_secondaryVertex_g, btag_information[btag_secondaryVertex][jetID], this_weight);
			fillHistoDataAndMC(h_btag_TC_highEff_g, btag_information[btag_TC_highEff][jetID], this_weight);
			fillHistoDataAndMC(h_btag_TC_highPur_g, btag_information[btag_TC_highPur][jetID], this_weight);
			fillHistoDataAndMC(h_btag_softEle_g, btag_information[btag_softEle][jetID], this_weight);
			fillHistoDataAndMC(h_btag_softMuon_g, btag_information[btag_softMuon][jetID], this_weight);
			fillHistoDataAndMC(h_btag_softMuonNoIP_g, btag_information[btag_softMuonNoIP][jetID], this_weight);
		}

	}
}
//-- eof ------------------------------------------------------------------------------------------
