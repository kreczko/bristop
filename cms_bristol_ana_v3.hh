#include "TSelector.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "NTupleReader.h"

#include <vector>
#include <iostream>
#include <map>
//FIXME: remove using namespace from header file
using namespace std;
typedef unsigned short ushort;

#ifndef CMS_BRISTOL_ANA_HH
#define CMS_BRISTOL_ANA_HH
// Global variables/constants
#define WMASS 80.398
#define WMASSERROR 0.025
#define ZMASS 91.1876
const int ntjet(5);
const int myprec(1); //no of decimal point for weighted nEvent
const int nbm3 = 960;
const bool m3_use_1000_bins = false;
const bool run_on_octX_skim = 0; // <---- set temporary swith here
const bool use_old_Z_veto = false; //TEMPORARY
//TODO: remove one of the two. Reason: they have the same information
//FIXME: 01st step to new sample - add mcname and mclabel
const string mcname[] = { "data", "ttbar", "QCD", "enri1", "enri2", "enri3", "bce1", "bce2", "bce3", "wj", "zj", "vqq",
		"singleTop", "tW", "tchan", "schan", "Zprime_M500GeV_W5GeV", "Zprime_M500GeV_W50GeV", "Zprime_M750GeV_W7500MeV",
		"Zprime_M1TeV_W10GeV", "Zprime_M1TeV_W100GeV", "Zprime_M1250GeV_W12500MeV", "Zprime_M1500GeV_W15GeV",
		"Zprime_M1500GeV_W150GeV", "Zprime_M2TeV_W20GeV", "Zprime_M2TeV_W200GeV", "Zprime_M3TeV_W30GeV", "Zprime_M3TeV_W300GeV",
		"Zprime_M4TeV_W40GeV", "Zprime_M4TeV_W400GeV" };
const string mclabel[] = { "data", "signal", "QCD", "enri1", "enri2", "enri3", "bce1", "bce2", "bce3", "W+jets", "Z+jets", "VQQ",
		"singleTop", "tW", "t-chan", "s-chan", "Zprime_M500GeV_W5GeV", "Zprime_M500GeV_W50GeV", "Zprime_M750GeV_W7500MeV",
		"Zprime_M1TeV_W10GeV", "Zprime_M1TeV_W100GeV", "Zprime_M1250GeV_W12500MeV", "Zprime_M1500GeV_W15GeV",
		"Zprime_M1500GeV_W150GeV", "Zprime_M2TeV_W20GeV", "Zprime_M2TeV_W200GeV", "Zprime_M3TeV_W30GeV", "Zprime_M3TeV_W300GeV",
		"Zprime_M4TeV_W40GeV", "Zprime_M4TeV_W400GeV" };
//make static in order to avoid ld duplication
//static const char* histnames[] = { "neutrino_pz", "neutrino_pz_mc", "mttbar", "mttbar_mc", "mZprime_mc", "mWlep", "mWlep_mc",
//		"mWhad", "mWhad_mc", "minDeltaR_ele_Jet", "ptRel_ele_jet", "mtlep", "mtlep_mc", "mthad", "mthad_mc", "thad_pt",
//		"thad_pt_mc", "tlep_pt", "tlep_pt_mc", "angle_b_ele", "ptratio", "pttbar", "htsystem", "angle_b_ele_matched",
//		"mtlep_matched", "mthad_matched", "mWlep_matched", "mWhad_matched", "ptratio_matched", "pttbar_matched",
//		"htsystem_matched", "WhadPartons", "Chi2Leptonic", "Chi2Leptonic_matched", "Chi2Hadronic", "Chi2Hadronic_matched",
//		"Chi2Global", "Chi2Global_matched", "Chi2Total", "Chi2Total_matched" };
//static const char* histnames2D[] = { "kptRel_vs_deltaRmin" };
const short int mcsize = sizeof(mcname) / sizeof(mcname[0]);
const int nmctype(mcsize + 7); //extend to include wj, zj, QCD, VQQ, single top
const int nstage(mcsize + 8); //add >=1Tele


class ana: public NTupleReader {

public:
	// The output file with histograms
	TFile* histf;

	//The text file with the candidates
	FILE* outfile;

	ana();
	~ana();

	void End();

	bool EventLoop();// the main analysis

	//Methods to call from anascript.C
	void SetInputFile(const char* fname);
	void SetOutputFirstName(const string name);
	void SetOutputHistFile(const string name, const string mode = "RECREATE");
	void SetOutputTextFile(const string name, const string mode = "w");

	void SetGoodRuns(bool f) {
		keepgood = f;
	}
	; // use all runs or just good runs?
	void SetData(bool f) {
		datafile = f;
	}
	; // is this a data file
	void CheckTrigger(bool); // check trigger fired?  time-consuming...
	void CheckTrigger(bool, const string hlt);
	void SetLimit(int n); // testing first few events
	void EstimateQCD(); // run QCD estimation
	bool EstimateQCD(const string file); // run QCD estimation
	void EstimateWjets(); // run W+jets estimation
	bool EstimateWjets(const string data, const string mc = ""); // run W+jets estimation
	void SetNtoyForM3Fit(int val) {
		m_ntoy = val;
	}
	; //number of toy exp to run for m3

	void SetEleETcut(float);
	void SetMuonPTcut(float);
	void SetJetETcut(float);
	void SetMETcut(float);
	void SetHTcut(float);

	// electron ID
	enum eID {
		robustTight, robustLoose, loose, tight, robustHighenergy, none
	};
	/**
	 * Enum for MC types
	 */
	enum MCType {
		kdata, ksignal, kQCD, kenri1, kenri2, kenri3, kbce1, kbce2, kbce3, kWjets, kZjets, kVQQ, ksingleTop, ktW, ktchan, kschan,
		Zprime_M500GeV_W5GeV, Zprime_M500GeV_W50GeV, Zprime_M750GeV_W7500MeV, Zprime_M1TeV_W10GeV, Zprime_M1TeV_W100GeV,
		Zprime_M1250GeV_W12500MeV, Zprime_M1500GeV_W15GeV, Zprime_M1500GeV_W150GeV, Zprime_M2TeV_W20GeV, Zprime_M2TeV_W200GeV,
		Zprime_M3TeV_W30GeV, Zprime_M3TeV_W300GeV, Zprime_M4TeV_W40GeV, Zprime_M4TeV_W400GeV, kNumMCTypes
	};

	/**
	 * Enum for 1-D histograms
	 */
	enum EHist {
		kneutrino_pz, kneutrino_pz_mc, kMttbar, kMttbar_mc, kMZprime_mc, kmWlep, kmWlep_mc, kmWhad, kmWhad_mc,
		kminDeltaR_ele_Jet, kptRel_ele_jet, kmtlep, kmtlep_mc, kmthad, kmthad_mc, kthad_pt, kthad_pt_mc, ktlep_pt, ktlep_pt_mc,
		kangle_b_ele, kptratio, kpttbar, khtsystem, kangle_b_ele_matched, kmtlep_matched, kmthad_matched, kmWhad_matched,
		kmWlep_matched, kptratio_matched, kpttbar_matched, khtsystem_matched, kWhadPartons, kChi2Leptonic, kChi2Leptonic_matched,
		kChi2Hadronic, kChi2Hadronic_matched, kChi2Global, kChi2Global_matched, kChi2Total, kChi2Total_matched, ktlep_pt_matched,
		kthad_pt_matched, kangle_b_ele_mc, kptratio_mc, kpttbar_mc, khtsystem_mc, kChi2Leptonic_mc, kChi2Hadronic_mc,
		kChi2Global_mc, kChi2Total_mc, kNumHists
	};

	/**
	 * Enum for 2-D histograms
	 */
	enum E2DHist {
		k2D_ptRel_vs_deltaRmin, k2D_NumHists
	};

	/**
	 * Enum for MC ttbar event partons
	 */
	enum MCEventPartons {
		kthad, ktlep, kelectron, kneutrino, kbhad, kblep, kWhad, kWlep, kq1, kq2, kNumMCEventPartons
	};

	/**
	 * vector for 1D histograms
	 */
	vector<vector<TH1F*> > fasthist_;
	/**
	 * vector for 2D histograms
	 */
	vector<vector<TH2F*> > fasthist2D_;
	/**
	 * 1D histogram names
	 */
	vector<string> histnames_;
	/**
	 * 1D histogram names (const char*) for histogram creation
	 */
	vector<const char*> histnames_c_;
	/**
	 * 2D histograms names
	 */
	vector<string> histnames2D_;
	/**
	 * 2D histograms names (const char*) for histogram creation
	 */
	vector<const char*> histnames2D_c;
	/**
	 * MC type flag
	 */
	MCType fastmctype_;
	/**
	 * Flag if MC type is present in dataset
	 */
	vector<bool> isMCTypePresent_;
	/**
	 * Set electron ID
	 * @param val new electron ID
	 */
	void SetEleID(eID val) {
		m_eID = val;
	}
	;

	void SetAESHTcut(float cut) {
		AES_HT_cut = cut;
	}
	;
	void SetAESMETcut(float cut) {
		AES_MET_cut = cut;
	}
	;
	void SetAESZveto_TwoEle(bool f) {
		AES_useSimpleZveto = f;
	}
	;

	void ApplyMETcut(bool f) {
		m_applyMETcut = f;
	}
	;
	void RejectEndcapEle(bool f) {
		m_rejectEndcapEle = f;
	}
	;

	void StudySystematics(const string, const string);
	void StudyPDFunc(bool f) {
		m_studyPDFunc = f;
	}
	;

	// Switches
	void Validation(bool val) {
		m_doValidation = val;
	}
	;
	void ConversionStudySwitch(bool val) {
		m_ConversionStudies = val;
	}
	;
	void PlotRelisoNES(bool val) {
		m_plotRelisoNES = val;
	}
	;
	void SetJetAlgo(string val) {
		m_jetAlgo = val;
	}
	;
	void SetMETAlgo(string val) {
		m_metAlgo = val;
	}
	;
	void SetLHCEnergyInTeV(double val) {
		m_LHCEnergyInTeV = val;
	}
	;
	void SetRunOnSD(bool val) {
		m_runOnSD = val;
	}
	;
	void SetRunOnMyHLTskim(bool val) {
		m_runOnMyHLTskim = val;
	}
	;

	//switch to preclue missing layers as these are not in 314 data samples
	void UseMissLayers(bool val) {
		m_useMisslayers = val;
	}
	;

private:

	//-----------------
	// Intenal methods
	//-----------------
	bool KeepGoodRuns() const {
		return keepgood;
	}
	bool IsData() const {
		return datafile;
	}
	;
	//	;
	int GetLimit() const {
		return numlimit;
	}
	;
	eID EleID() const {
		return m_eID;
	}
	;
	void PrintCuts() const; // print kinematic cuts
	string printEleID() const;
	void CheckAvailableJetMET();
	//	void Init(); //initialize tree branches
	//	void ReadSelectedBranches() const;

	//geometry
	float calcDeltaR(const float phi1, const float eta1, const float phi2, const float eta2) const;
	float calcDeltaR(const TLorentzVector& p1, const TLorentzVector& p2) const;

	float calcDeltaEta(const TLorentzVector& p1, const TLorentzVector& p2) const;
	float calcDeltaPhi(const TLorentzVector& p1, const TLorentzVector& p2) const;

	bool ConversionFinder(const TLorentzVector& e1, int mctype, int index_selected_ele);
	bool ConversionFinder2(const TLorentzVector& e1, int mctype, int index_selected_ele);
	void ConversionMCMatching(const TLorentzVector& e1, int mctype, bool isthisConversion);
	float MCTruthMatch(float eta, float phi);
	short match(TLorentzVector, vector<TLorentzVector> );
	void CreateHistogramNames();
	void Create2DHistogramNames();

	// QCD estimation
	pair<double, double> estimateQCD_computeFitResultUncertainty(const double est, TF1* f) const;
	double estimateQCD_newEstimateByVaryingFitFunction(const double p1, const double p2, const double p3) const;
	pair<double, double> estimateQCD_assign_pos_neg_estimate(const double est, const double new_est[2]) const;
	void Set_Reliso_bin_width(float bw) {
		m_QCDest_reliso_bin_width = bw;
	}
	;
	float Get_Reliso_bin_width() const {
		return m_QCDest_reliso_bin_width;
	}
	;

	// Histograms
	void addHistoNjet(TH1F* h[], const string, const string, const string, const int, const float, const float);
	void addHistoNjet(TH2F* h[], const string, const string, const string, const int, const float, const float, const int,
			const float, const float);
	void fillHistoNjet(TH1F* h[], const float value, const double weight);
	void fillHistoNjet(TH2F* h[], const float value, const float value2, const double weight);

	void fillHistoNjet_DataAndMC(const string hname, const float value, const double weight);
	void fillHistoNjet_DataAndMC(const string hname, const float v1, const float v2, const double weight);

	void fillHistoNjet2D(TH1F* h[][mcsize], const int ec, const float value, const double weight);

	// Helper methods to fill histograms for data and each type of MC
	void addHistoDataAndMC(TH1F* h[], const string, const string, const int, const float, const float) const;
	void addHistoDataAndMC(TH2F* h[], const string, const string, const int, const float, const float, const int, const float,
			const float) const;
	void fillHistoDataAndMC(TH1F* h[], const float value, const double weight) const;
	void fillHistoDataAndMC(TH2F* h[], const float v1, const float v2, const double weight) const;

	void addHisto_Njet_DataAndMC(TH1F* h[7][mcsize], const string, const string, const int, const float, const float);
	void fillHisto_Njet_DataAndMC(TH1F* h[7][mcsize], const float value, const double w);

	void fillHisto_PDF_weights(TH1F* h);
	void fillMCTopEventHists();

	// W+jets estimation
	void reco_hadronicTop_highestTopPT(const std::vector<TLorentzVector>&, const int nGoodIsoEle);
	pair<double, double> compute_M3(const std::vector<TLorentzVector>&) const;
	void reco_Mttbar(const std::vector<TLorentzVector>& jets, const TLorentzVector& electron, const TLorentzVector& met);
	void reco_Mttbar_matched(const std::vector<TLorentzVector>& jets, const TLorentzVector& electron, const TLorentzVector& met);
	TLorentzVector reconstruct_neutrino(const TLorentzVector&, const TLorentzVector&);
	pair<TLorentzVector, TLorentzVector> reconstruct_neutrinos(const TLorentzVector&, const TLorentzVector&);
	pair<ushort, ushort> reco_hadronic_W(const std::vector<TLorentzVector>& jets, short ommit_blep_id) const;
	short findClosest(const std::vector<TLorentzVector>&, const TLorentzVector&);

	vector<TLorentzVector> GetMCTopEvent();
	void SetHistoLabelCutNjet(TH2D *this_njetVcuts, vector<string>& ve) const;
	void SetHistoLabelEleID(TH1F *eid[]) const;

	void DefineCrossSection();
	double GetCrossSection(const string) const;
	void SetEventWeightMap();
	void SetWeights();
	double GetWeight(const string) const;
	double GetWeight(const short) const;
	long GetNinit(const string) const;

	// print event-count tables
	void DrawEventPerNjetTable(const double nevt[][5][nmctype], const vector<string>& ve) const;
	void DrawSignalBGTable(const double nevt[][5][nmctype], const vector<string> ve) const;
	void DrawMCTypeTable(const double nevt[][5][nmctype], const string title, vector<string> ve) const;
	//  void DrawMCTypeTable(const double nevt[][5][nmctype], const string title, vector<string> cuts, bool weighted) const;
	void DrawQCDTable(const double nevt[][5][nmctype], const string title, vector<string> ve) const;
	void DrawSingleTopTable(const double nevt[][5][nmctype], const string title, vector<string> ve) const;
	void DrawSignalAcceptanceTable(const double nevent[][5][nmctype], vector<string> ve) const;
	void printCutStage(int, string) const;
	void printCutStage(ofstream& os, int, string) const;
	double getTotalEvents(const double nevt[][5][nmctype], short cut, short k_start, short k_end, short njbegin) const;

	// print event-count tables (with errors)
	void PrintErrorTables(const double e_plus_jet[][5][nmctype], const double e_plus_jet_weighted[][5][nmctype],
			vector<string> ve) const;
	void PrintError_NjetVcut(ofstream&, const double[][5][nmctype], const double[][5][nstage], vector<string>&) const;
	double GetBayesUncertainty(int Ninitial) const;
	double GetChi2Leptonic(double Wmass, double tmass, double angle);
	double GetChi2Hadronic(double Wmass, double tmass, double ptratio);
	double GetChi2Global(double pttbar, double htsystem);
	double GetHT(const std::vector<TLorentzVector>& jets, ushort N);
	string ScrNum(double num) const;
	void printLine(ofstream &myfile, const double, const double) const;
	bool ScientificNotation;

	bool is_mc_present(const int) const;
	float compute_d0(const string, const int) const;
	float compute_mtw(const TVector2&, const TVector2&) const;
	pair<double, double> compute_neutrino_momentum_z(double met, const TLorentzVector& electron);
	void PrintGenParticles() const;

	// validation plots
	void valid_mkHisto_cut_njet(TH1F* h[][7], const string, const string, const int, const float, const float);
	void valid_fillHisto(TH1F* h[][7], const bool cuts[8], int nj, double value) const;

	int Njet() const {
		return m_nGoodJet;
	}
	; //num of cleaned jet in event
	bool doValidation() const {
		return m_doValidation;
	}
	;
	bool DoConversionStudies() {
		return m_ConversionStudies;
	}
	;
	// conversion
	void PrintConversionTable();
	int ConversionCounter;
	int ConversionArray[nmctype][2][6];
	void OptimiseConversionFinder(const TLorentzVector& e1, int mctype);

	TH2D *Conv_Opti[2];
	TH2D *Conv_Optis[2];
	TH2D *Conv_OptiL[2];
	TH2D *Conv_Opti_extragran[2];
	TH1D *Conv_CheckDelR_GSFTk_ctfTk;
	int mycounter;

	// new (to run on mixed MC) (not currently needed)
	string CheckEventTypeFromMcTruth() const;

	float getRelIso(int) const;
	bool passEleID(unsigned int) const;
	bool passHLT() const;
	string printTimeNow() const;
	void bookHistograms();
	void bookChi2MatchedHists(ushort type);
	void bookMCHists(ushort type);

	//--------------------
	// private variables
	//--------------------
	//	bool datafile;
	bool keepgood;
	//	bool checkTrig;
	//	string HLTBit;
	int numlimit;
	int nfile;
	string outputHistFileName;
	string outputTextFileName;
	double this_weight; // current event weight
	int m_nGoodJet;
	float m_QCDest_reliso_bin_width;
	bool m_doValidation;
	bool m_plotRelisoNES;
	//	bool m_debug;
	bool m_ConversionStudies;
	//	string m_jetAlgo;
	//	string m_metAlgo;
	double m_LHCEnergyInTeV;
	bool m_runOnSD;
	bool m_runOnMyHLTskim;
	//	bool m_useMisslayers;

	// cuts
	float ELE_ETCUT;
	float MU_PTCUT;
	float JET_PTCUT;
	float METCUT;
	float HTCUT;
	int nCutSetInScript;
	bool m_applyMETcut;
	bool m_rejectEndcapEle;
	eID m_eID;
	float AES_HT_cut;
	float AES_MET_cut;
	bool AES_useSimpleZveto;
	float intlumi; // integrated luminosity assumed
	bool useNewReliso;

	string doSystematics;
	string sysSample;
	bool m_studyPDFunc;

	map<string, double> cross_section;
	map<string, long> nInitialEventMC; //initial number of event, used to compute event weight
	map<string, double> weightMap;
	vector<double> fastWeight;
	// Wjet estimation
	TH1D *h_hadTop_maxPT_mass_4j; // 960 bins (0-960)
	TH1D *h_hadTop_maxPT_pt_4j;
	TH1D *h_hadTop_maxPT_mass_nonIso_4j;
	TH1D *h_hadTop_maxPT_pt_nonIso_4j;
	TH1D *h_m3_tt;
	TH1D *h_m3_wj;
	TH1D *h_m3_zj;
	TH1D *h_m3_qcd;
	TH1D *h_m3_vqq;
	TH1D *h_m3_singletop;
	TH1D *h_m3_bce[3];
	TH1D *h_m3_enri[3];
	TH1D *h_m3_tt_control;
	TH1D *h_m3_wj_control;
	TH1D *h_m3_zj_control;
	TH1D *h_m3_qcd_control;
	TH1D *h_m3_vqq_control;
	TH1D *h_m3_singletop_control;
	TH1D *h_m3_bce_control[3];
	TH1D *h_m3_enri_control[3];
	// 1000 bins (0-1000)
	TH1D *h_hadTop_maxPT_mass_4j_1000;
	TH1D *h_hadTop_maxPT_mass_nonIso_4j_1000;
	TH1D *h_m3_tt_1000;
	TH1D *h_m3_wj_1000;
	TH1D *h_m3_zj_1000;
	TH1D *h_m3_qcd_1000;
	TH1D *h_m3_vqq_1000;
	TH1D *h_m3_singletop_1000;
	TH1D *h_m3_bce_1000[3];
	TH1D *h_m3_enri_1000[3];
	TH1D *h_m3_tt_control_1000;
	TH1D *h_m3_wj_control_1000;
	TH1D *h_m3_zj_control_1000;
	TH1D *h_m3_qcd_control_1000;
	TH1D *h_m3_vqq_control_1000;
	TH1D *h_m3_singletop_control_1000;
	TH1D *h_m3_bce_control_1000[3];
	TH1D *h_m3_enri_control_1000[3];
	int m_ntoy;

	// MC flag
	void SetMCFlag();
	vector<string> mc_names;
	bool mc_sample_has_ttbar;
	bool mc_sample_has_Wjet;
	bool mc_sample_has_Zjet;
	bool mc_sample_has_QCD;
	bool mc_sample_has_enri1;
	bool mc_sample_has_enri2;
	bool mc_sample_has_enri3;
	bool mc_sample_has_bce1;
	bool mc_sample_has_bce2;
	bool mc_sample_has_bce3;
	bool mc_sample_has_VQQ;
	bool mc_sample_has_singleTop;
	bool mc_sample_has_tW;
	bool mc_sample_has_tchan;
	bool mc_sample_has_schan;
	//FIXME: 02nd step to new sample - declare mc_sample_has*
	bool mc_sample_has_Zprime_M500GeV_W5GeV, mc_sample_has_Zprime_M500GeV_W50GeV, mc_sample_has_Zprime_M750GeV_W7500MeV,
			mc_sample_has_Zprime_M1TeV_W10GeV, mc_sample_has_Zprime_M1TeV_W100GeV, mc_sample_has_Zprime_M1250GeV_W12500MeV,
			mc_sample_has_Zprime_M1500GeV_W15GeV, mc_sample_has_Zprime_M1500GeV_W150GeV, mc_sample_has_Zprime_M2TeV_W20GeV,
			mc_sample_has_Zprime_M2TeV_W200GeV, mc_sample_has_Zprime_M3TeV_W30GeV, mc_sample_has_Zprime_M3TeV_W300GeV,
			mc_sample_has_Zprime_M4TeV_W40GeV, mc_sample_has_Zprime_M4TeV_W400GeV;

	// MC type of this event
	bool isTTbar;
	bool isWjets;
	bool isZjets;
	bool isQCD;
	bool isEnri1;
	bool isEnri2;
	bool isEnri3;
	bool isBce1;
	bool isBce2;
	bool isBce3;
	bool isVQQ;
	bool isSingleTop;
	bool isTW;
	bool isTchan;
	bool isSchan;
	//FIXME: 03rd step to new sample - declare isMCtype
	bool isZprime_M500GeV_W5GeV, isZprime_M500GeV_W50GeV, isZprime_M750GeV_W7500MeV, isZprime_M1TeV_W10GeV,
			isZprime_M1TeV_W100GeV, isZprime_M1250GeV_W12500MeV, isZprime_M1500GeV_W15GeV, isZprime_M1500GeV_W150GeV,
			isZprime_M2TeV_W20GeV, isZprime_M2TeV_W200GeV, isZprime_M3TeV_W30GeV, isZprime_M3TeV_W300GeV, isZprime_M4TeV_W40GeV,
			isZprime_M4TeV_W400GeV;

};

#endif
