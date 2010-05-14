#include "TSelector.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "NTupleReader.h"
#include "Tools/Logger.hh"
#include <vector>
#include <iostream>
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
const string mcname[] = { "data", "ttbar", "ttjet", "wj", "zj", "enri1", "enri2", "enri3", "bce1", "bce2", "bce3", "vqq", "tW",
		"tchan", "schan", "Zprime_M500GeV_W5GeV", "Zprime_M500GeV_W50GeV", "Zprime_M750GeV_W7500MeV", "Zprime_M1TeV_W10GeV",
		"Zprime_M1TeV_W100GeV", "Zprime_M1250GeV_W12500MeV", "Zprime_M1500GeV_W15GeV", "Zprime_M1500GeV_W150GeV",
		"Zprime_M2TeV_W20GeV", "Zprime_M2TeV_W200GeV", "Zprime_M3TeV_W30GeV", "Zprime_M3TeV_W300GeV", "Zprime_M4TeV_W40GeV",
		"Zprime_M4TeV_W400GeV", "QCD", "singleTop" };

//TODO: remove them and put them in where the mcnames are filled
const string mcfiles[] = { "data", "ttbar", "ttjet", "wjet", "zjet", "enri1", "enri2", "enri3", "bce1", "bce2", "bce3", "vqq",
		"tW", "tchan", "schan", "Zprime_M500GeV_W5GeV", "Zprime_M500GeV_W50GeV", "Zprime_M750GeV_W7500MeV",
		"Zprime_M1TeV_W10GeV", "Zprime_M1TeV_W100GeV", "Zprime_M1250GeV_W12500MeV", "Zprime_M1500GeV_W15GeV",
		"Zprime_M1500GeV_W150GeV", "Zprime_M2TeV_W20GeV", "Zprime_M2TeV_W200GeV", "Zprime_M3TeV_W30GeV", "Zprime_M3TeV_W300GeV",
		"Zprime_M4TeV_W40GeV", "Zprime_M4TeV_W400GeV", "QCD", "singleTop" };

const string mclabel[] = { "data", "signal", "signalJ", "W+jets", "Z+jets", "enri1", "enri2", "enri3", "bce1", "bce2", "bce3",
		"VQQ", "tW", "t-chan", "s-chan", "Zprime_M500GeV_W5GeV", "Zprime_M500GeV_W50GeV", "Zprime_M750GeV_W7500MeV",
		"Zprime_M1TeV_W10GeV", "Zprime_M1TeV_W100GeV", "Zprime_M1250GeV_W12500MeV", "Zprime_M1500GeV_W15GeV",
		"Zprime_M1500GeV_W150GeV", "Zprime_M2TeV_W20GeV", "Zprime_M2TeV_W200GeV", "Zprime_M3TeV_W30GeV", "Zprime_M3TeV_W300GeV",
		"Zprime_M4TeV_W40GeV", "Zprime_M4TeV_W400GeV", "QCD", "singleTop" };
const short int mcsize = sizeof(mcname) / sizeof(mcname[0]);
const int nmctype(mcsize + 8); //extend to include wj, zj, QCD, VQQ, single top
const int nstage(16); //add >=1Tele


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
	}//number of toy exp to run for m3

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
		kdata, kttbar, kttjet, kWjets, kZjets, kenri1, kenri2, kenri3, kbce1, kbce2, kbce3, kVQQ, ktW, ktchan, kschan,
		Zprime_M500GeV_W5GeV, Zprime_M500GeV_W50GeV, Zprime_M750GeV_W7500MeV, Zprime_M1TeV_W10GeV, Zprime_M1TeV_W100GeV,
		Zprime_M1250GeV_W12500MeV, Zprime_M1500GeV_W15GeV, Zprime_M1500GeV_W150GeV, Zprime_M2TeV_W20GeV, Zprime_M2TeV_W200GeV,
		Zprime_M3TeV_W30GeV, Zprime_M3TeV_W300GeV, Zprime_M4TeV_W40GeV, Zprime_M4TeV_W400GeV, kQCD, ksingleTop/*30*/,
		NUMBER_OF_MC_TYPES
	};

	/**
	 * Enum for 1-D histograms
	 */
	enum EHist {
		h_neutrino_pz, h_neutrino_pz_mc, h_mttbar_matched, h_mttbar_mc, h_mttbar_mc_smeared, h_mttbar_diff_reco_and_mc,
		h_mZprime_mc, h_mWlep, h_mWlep_mc, h_mWhad_mc, h_minDeltaR_ele_Jet, h_ptRel_ele_jet, h_mtlep_mc, h_mthad_mc, h_thad_pt,
		h_thad_pt_mc, h_tlep_pt, h_tlep_pt_mc, h_angle_b_ele, h_ptratio, h_pttbar, h_htsystem, h_angle_b_ele_matched,
		h_mtlep_matched, h_mthad_matched, h_mWhad_matched, h_mWlep_matched, h_ptratio_matched, h_ptratio2_matched, h_ptratio2_mc,
		h_pttbar_matched, h_htsystem_matched, h_Chi2Leptonic, h_Chi2Leptonic_matched, h_Chi2Hadronic, h_Chi2Hadronic_matched,
		h_Chi2Global, h_Chi2Global_matched, h_Chi2Total, h_Chi2Total_matched, h_tlep_pt_matched, h_thad_pt_matched,
		h_angle_b_ele_mc, h_ptratio_mc, h_pttbar_mc, h_htsystem_mc, h_Chi2Leptonic_mc, h_Chi2Hadronic_mc, h_Chi2Global_mc,
		h_Chi2Total_mc,/*end of exotic top*/
		h_nele, h_ele_ET_all, h_ele_ET_1, h_ele_ET_2, h_ele_ET_3, h_ele_eta_all, h_ele_eta_1, h_ele_eta_2, h_ele_eta_3, h_ele_phi_all,
		h_ele_phi_1, h_ele_phi_2, h_ele_phi_3, h_ele_iso_all, h_ele_iso_1, h_ele_iso_2, h_ele_iso_3, h_nele_cuts, h_eid,/*electrons*/
		h_njets, h_jet_pt_all, h_jet_pt_1, h_jet_pt_2, h_jet_pt_3, h_jet_pt_4, h_jet_eta_all, h_jet_eta_1, h_jet_eta_2, h_jet_eta_3,
		h_jet_eta_4, h_jet_phi_all, h_jet_phi_1, h_jet_phi_2, h_jet_phi_3, h_jet_phi_4,/*jets*/h_metAlone, h_metAlone_phi,
		h_DRemu_selE_GoodMu, h_DRemu_selE_GoodMu_pass, h_exp_ele_et, h_exp_ele_eta, h_exp_j0_pt, h_exp_j1_pt, h_exp_DRej,
		h_exp_DPhiej, h_exp_DRjj, h_exp_DPhijj, h_nGenBasicEle_Zee_allj, h_Zee_eta, h_Zee_pt, h_Z_photon_eta, h_Z_photon_et,
		h_Zee_photon_eta, h_Zee_photon_et, h_Z_Nphotons, h_Zee_Nphotons, h_mass_diele, h_mass_diele_new, h_mass_diele_lowMet_1j,
		h_mass_ephoton_lowMet_1j, h_Nele_lowMet_1j, h_Nphoton_lowMet_1j, h_photon_eta_lowMet_1j, h_photon_et_lowMet_1j,
		h_photon1_eta_lowMet_1j, h_photon1_et_lowMet_1j, h_ed0_unCor, h_ed0, h_ed0_pass, h_muon_chi2, h_muon_d0_unCor, h_muon_d0,
		h_muon_hits, h_met_ante_ISO, h_met_ante_ISO_mu, h_met_ante_ISO_t1, h_mtw_mu_incl, h_mtw_t1_incl, h_DPhiEmet_mu_incl,
		h_DPhiEmet_t1_incl, h_btag_TC_highEff_b, h_btag_TC_highEff_c, h_btag_TC_highEff_uds, h_btag_TC_highEff_g,
		h_btag_TC_highPur_b, h_btag_TC_highPur_c, h_btag_TC_highPur_uds, h_btag_TC_highPur_g, h_btag_JetBProb_b,
		h_btag_JetBProb_c, h_btag_JetBProb_uds, h_btag_JetBProb_g, h_btag_JetProb_b, h_btag_JetProb_c, h_btag_JetProb_uds,
		h_btag_JetProb_g, h_btag_secondaryVertex_b, h_btag_secondaryVertex_c, h_btag_secondaryVertex_uds,
		h_btag_secondaryVertex_g, h_btag_softEle_b, h_btag_softEle_c, h_btag_softEle_uds, h_btag_softEle_g, h_btag_softMuon_b,
		h_btag_softMuon_c, h_btag_softMuon_uds, h_btag_softMuon_g, h_btag_softMuonNoIP_b, h_btag_softMuonNoIP_c,
		h_btag_softMuonNoIP_uds, h_btag_softMuonNoIP_g, h_numberOfBtags, h_hadTop_maxPT_mass_4j, h_hadTop_maxPT_pt_4j,
		h_hadTop_maxPT_mass_nonIso_4j, h_hadTop_maxPT_pt_nonIso_4j, h_m3, h_m3_control, h_hadTop_maxPT_mass_4j_1000,
		h_hadTop_maxPT_mass_nonIso_4j_1000, h_m3_1000, h_m3_control_1000, NUMBER_OF_1D_HISTOGRAMS
	};

	/**
	 * Enum for 1-D njet binned histograms
	 */
	enum EJbinHist {
		h_nEle_all, h_nEle_s1, h_nEle_s2, h_nEle_s3_idLoose, h_nEle_s3_idTight, h_nEle_s3_idRL, h_nEle_s3_idRT, h_ht, h_met,
		h_met_mu, h_met_t1, h_met_BA, h_met_mu_BA, h_met_t1_BA, h_mtw_mu, h_mtw_t1, h_DPhiEmet_mu, h_DPhiEmet_t1,
		h_QCDest_CombRelIso, h_QCDest_CombRelIso_AES, h_QCDest_CombRelIso_AES_minusMET, h_QCDest_CombRelIso_AES_minusHT,
		h_QCDest_CombRelIso_AES_minusTighterZ, h_QCDest_CombRelIso_AES_before, h_QCDest_CombRelIso_AES_justMET,
		h_QCDest_CombRelIso_AES_justHighMET, h_QCDest_CombRelIso_AES_justHT, h_QCDest_CombRelIso_AES_justZ,
		h_QCDest_CombRelIso_AES_planA1_e20, h_QCDest_CombRelIso_AES_planA1_e30, h_QCDest_CombRelIso_AES_planA2_e20,
		h_QCDest_CombRelIso_AES_planA2_e30, h_QCDest_CombRelIso_AES_planA3_e20, h_QCDest_CombRelIso_AES_planA3_e30,
		h_QCDest_CombRelIso_AES_planB1_e20, h_QCDest_CombRelIso_AES_planB1_e30, h_QCDest_CombRelIso_AES_planB2_e20,
		h_QCDest_CombRelIso_AES_planB2_e30, h_QCDest_CombRelIso_AES_planB3_e20, h_QCDest_CombRelIso_AES_planB3_e30,
		h_QCDest_CombRelIso_AES_planB3b_e20, h_QCDest_CombRelIso_AES_planB3b_e30, h_QCDest_CombRelIso_AES_planB4_e20,
		h_QCDest_CombRelIso_AES_planB4_e30, h_QCDest_CombRelIso_AES_planB5_e20, h_QCDest_CombRelIso_AES_planB5_e30,
		h_QCDest_CombRelIso_AES_planB6_e20, h_QCDest_CombRelIso_AES_planB6_e30, h_QCDest_CombRelIso_AES_planB7_e20,
		h_QCDest_CombRelIso_AES_planB7_e30, h_QCDest_CombRelIso_AES_planB8_e20, h_QCDest_CombRelIso_AES_planB8_e30,
		NUMBER_OF_JETBINNED_HISTOGRAMS
	};

	enum ENlevelJbinned2DHists {
		h_QCDest_isoVmet_NES, h_QCDest_isoVmet_NES_barrel, h_QCDest_isoVmet_NES_endcap, h_QCDest_isoVmet_NES_uw,
		h_QCDest_isoVmet_NES_uw_barrel, h_QCDest_isoVmet_NES_uw_endcap, NUMBER_OF_NLEVEL_JBINNED_2D_HISTOGRAMS
	};

	enum ENlevelJbinned1DHists {
		h_QCDest_CombRelIso_NES, h_QCDest_CombRelIso_NES_barrel, h_QCDest_CombRelIso_NES_endcap, h_QCDest_CombRelIso_NES_loMET,
		h_QCDest_CombRelIso_NES_loMET_barrel, h_QCDest_CombRelIso_NES_loMET_endcap, h_QCDest_CombRelIso_NES_hiMET,
		h_QCDest_CombRelIso_NES_hiMET_barrel, h_QCDest_CombRelIso_NES_hiMET_endcap, NUMBER_OF_NLEVEL_JBINNED_1D_HISTOGRAMS
	};

	enum Ebtag1DHists {
		h_mttbar, h_mWhad, h_mtlep, h_mthad, NUMBER_OF_BTAG_1D_HISTOGRAMS
	};
	/**
	 * Enum for 2-D histograms
	 */
	enum E2DHist {
		h_ptRel_vs_deltaRmin, h_mWhad_vs_mthad, h_mWhad_vs_mthad_btag_fake, k2D_mWlep_vs_mtlep, k2D_mWlep_vs_mtlep_btag_fake,
		h_exp_met_v_eeta, h_Zee_photon_eteta_2D, h_ptratio_vs_mttbar, h_htsystem_vs_mttbar, h_pttbar_vs_mttbar,
		h_angle_b_e_vs_mttbar, NUMBER_OF_2D_HISTOGRAMS
	};

	enum Ejetbins {
		NoJet, OneJet, TwoJets, ThreeJets, FourJets, FourOrMoreJets, AllJets, NUMBER_OF_JETBINS
	};

	enum ECuts {
		kNumCuts
	};

	/**
	 * Enum for MC ttbar event partons
	 */
	enum MCEventPartons {
		kthad, ktlep, kelectron, kneutrino, kbhad, kblep, kWhad, kWlep, kq1, kq2, NUMBER_OF_EVENT_PARTONS
	};

	enum MCDecayBranch {
		evqq, mvqq, tvqq, evev, mvmv, tvtv, evmv, evtv, mvtv, qqqq, NUMBER_OF_DECAY_BRANCHES
	};

	enum ELevel {
		Level_1, Level_1b, Level_1c, Level_1d1, Level_1d2, Level_1d3, Level_1d4, Level_1d5, Level_2, Level_3, Level_4,
		NUMBER_OF_LEVELS
	};

	enum Enum_BTAGS {
		btag_type_none, btag_type_fake, btag_type_TrackCount_highEff, btag_type_TrackCount_highPur, btag_type_JetBProb,
		btag_type_JetProb, btag_type_secondaryVertex, btag_softEle, btag_type_softMuon, btag_type_softMuonNoIP, NUMBER_OF_BTAGS
	};
private:
	Logger *logger;
	vector<string> mc_names;
	/**
	 * vector for 1D histograms
	 */
	vector<vector<TH1F*> > fasthist_;
	vector<vector<vector<TH1F*> > > btag_hist_;
	/**
	 * vector for 1D histograms
	 */
	vector<vector<vector<TH1F*> > > fasthist_jetbinned_;
	/**
	 * vector for 2D histograms
	 */
	vector<vector<TH2F*> > fasthist2D_;
	//mctype, level, jetbin, hist
	vector<vector<vector<vector<TH1D*> > > > fasthistNlevel1D_;
	vector<vector<vector<vector<TH2D*> > > > fasthistNlevel2D_;
	/**
	 * MC type flag
	 */
	MCType fastmctype_;
	/**
	 * Flag if MC type is present in dataset
	 * @deprecated use isMCPresent(MCType) instead
	 */
	vector<bool> isMCTypePresent_;

	/**
	 * Initial number of events per MCType
	 */
	vector<long> nInitialEventMC_;
	/**
	 * Cross sections for different MCTypes
	 */
	vector<double> cross_section_;

	/**
	 * Events weights for different MCTypes
	 */
	vector<double> fastWeight_;
	short number_of_jets_to_use_for_reco;
	bool useIsoElectronForReco;
	std::vector<std::vector<double> > btag_information;

	ushort N_positive_btags, N_negative_btags;
public:

	void ResetBtagCount(const std::vector<TLorentzVector> & jets, const ushort btag_type);
	void SetNumberOfJetsUsedInReco(short number) {
		number_of_jets_to_use_for_reco = number;
	}

	void UseIsoElectronForReco(bool flag) {
		useIsoElectronForReco = flag;
	}
	/**
	 * Set electron ID
	 * @param val new electron ID
	 */
	void SetEleID(eID val) {
		m_eID = val;
	}

	void SetAESHTcut(float cut) {
		AES_HT_cut = cut;
	}
	void SetAESMETcut(float cut) {
		AES_MET_cut = cut;
	}
	void SetAESZveto_TwoEle(bool f) {
		AES_useSimpleZveto = f;
	}

	void ApplyMETcut(bool f) {
		m_applyMETcut = f;
	}
	void RejectEndcapEle(bool f) {
		m_rejectEndcapEle = f;
	}

	void StudySystematics(const string, const string);
	void StudyPDFunc(bool f) {
		m_studyPDFunc = f;
	}

	// Switches
	void Validation(bool val) {
		m_doValidation = val;
	}
	void ConversionStudySwitch(bool val) {
		m_ConversionStudies = val;
	}
	void PlotRelisoNES(bool val) {
		m_plotRelisoNES = val;
	}
	void SetJetAlgo(string val) {
		m_jetAlgo = val;
	}
	void SetMETAlgo(string val) {
		m_metAlgo = val;
	}
	void SetLHCEnergyInTeV(double val) {
		m_LHCEnergyInTeV = val;
	}
	void SetRunOnSD(bool val) {
		m_runOnSD = val;
	}
	void SetRunOnMyHLTskim(bool val) {
		m_runOnMyHLTskim = val;
	}

	//switch to preclue missing layers as these are not in 314 data samples
	void UseMissLayers(bool val) {
		m_useMisslayers = val;
	}

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

	double smearValueWithGaus(double value, double gaus_width);
	//	;
	int GetLimit() const {
		return numlimit;
	}
	eID EleID() const {
		return m_eID;
	}

	//	double btag(ushort jetid, ushort btag_type);
	void PrintCuts() const; // print kinematic cuts
	string printEleID() const;
	void CheckAvailableJetMET();
	//	void Init(); //initialize tree branches
	//	void ReadSelectedBranches() const;

	//geometry
	float calcDeltaR(const float phi1, const float eta1, const float phi2, const float eta2) const;
	/**
	 * @deprecated
	 * @param p1
	 * @param p2
	 * @return
	 */
	float calcDeltaR(const TLorentzVector& p1, const TLorentzVector& p2) const;

	float calcDeltaEta(const TLorentzVector& p1, const TLorentzVector& p2) const;
	/**
	 * @deprecated
	 * @param p1
	 * @param p2
	 * @return
	 */
	float calcDeltaPhi(const TLorentzVector& p1, const TLorentzVector& p2) const;

	bool ConversionFinder(const TLorentzVector& e1, int mctype, int index_selected_ele);
	bool ConversionFinder2(const TLorentzVector& e1, int mctype, int index_selected_ele);
	void ConversionMCMatching(const TLorentzVector& e1, int mctype, bool isthisConversion);
	float MCTruthMatch(float eta, float phi);
	short match(TLorentzVector, vector<TLorentzVector> );

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
	//	void addHistoNjet(TH1F* h[], const string, const string, const string, const int, const float, const float);
	//	void addHistoNjet(TH2F* h[], const string, const string, const string, const int, const float, const float, const int,
	//			const float, const float);
	//	void addHisto_Njet_DataAndMC(TH1F* h[7][mcsize], const string, const string, const int, const float, const float);
	void addHisto_Njet_DataAndMC(ushort, ushort, const string, const string, const int, const float, const float);
	void addHisto_btag_DataAndMC(ushort, ushort, const string, const string, const int, const float, const float);
	void addHistoDataAndMC(ushort, ushort, const string, const string, const int, const float, const float);
	void addHistoDataAndMC(ushort, ushort, const string, const string, const int, const float, const float, const int,
			const float, const float);
	void iso_addHisto_nlevel_nj_nmc(ushort, ushort, const string, const string, const int, const float, const float);
	void iso_addHisto_nlevel_nj_nmc(ushort, ushort, const string, const string, const int, const float, const float, const int,
			const float, const float);

	void bookHistograms();
	void bookExoticTopHistograms(ushort type);
	void bookExoticTopMatchedHistograms(ushort type);
	void bookExoticTopBtaggedHistograms(ushort type);
	void bookExoticTopMCHistograms(ushort type);
	void bookBasicHistograms(ushort type);
	void bookDRemuHistograms(ushort type);
	void bookExploreHistograms(ushort type);
	void bookZvetoHistograms(ushort type, TDirectory *parent);
	void bookElectronCountHistograms(ushort type);
	void bookNewKinematicHistograms(ushort type);
	void bookQCDEstimationHistograms(ushort type, TDirectory *parent);
	void bookWjetsEstimationHistograms(ushort type, TDirectory *parent);
	void bookBtagHistograms(ushort type, TDirectory *parent);

	//	const char* dynTitle(string title, ushort type);
	//	const char* jbinHName(string name, ushort jbin);
	//	const char* jbinHTitle(string title, ushort jbin, ushort type);

	// Helper methods to fill histograms for data and each type of MC
	//	void fillHistoDataAndMC(TH1F* h[], const float value, const double weight) const;
	void fillHistoDataAndMC(ushort hist, const float value, const double weight);
	void fillHistoDataAndMC(ushort hist, const float Xvalue, const float Yvalue, const double weight);
	void iso_fillHisto_nlevel_nj_nmc(const ushort hist, ushort level, const float value, const float weight);
	void iso_fillHisto_nlevel_nj_nmc(const ushort hist, ushort level, const float valueX, const float valueY, const float weight);
	//	void fillHistoDataAndMC(TH2F* h[], const float v1, const float v2, const double weight) const;
	//	void fillHisto_Njet_DataAndMC(TH1F* h[7][mcsize], const float value, const double w);
	void fillHisto_Njet_DataAndMC(ushort hist, const float value, const double w);
	void fillHisto_btag_DataAndMC(ushort btag, ushort hist, const float value, const double w);
	void fillHisto_PDF_weights(TH1F* h);
	void fillMCTopEventHists(double reco_mttbar);
	void fillBtagHistograms(std::vector<TLorentzVector> jets);//TODO: pass size instead of whole vector
	//	void fillHistoNjet(TH1F* h[], const float value, const double weight);
	//	void fillHistoNjet(TH2F* h[], const float value, const float value2, const double weight);
	//	void fillHistoNjet_DataAndMC(const string hname, const float value, const double weight);
	//	void fillHistoNjet_DataAndMC(const string hname, const float v1, const float v2, const double weight);
	//	void fillHistoNjet2D(TH1F* h[][mcsize], const int ec, const float value, const double weight);

	// W+jets estimation
	void reco_hadronicTop_highestTopPT(const std::vector<TLorentzVector>&, const int nGoodIsoEle);
	pair<double, double> compute_M3(const std::vector<TLorentzVector>&) const;
	void reco_Mttbar(const std::vector<TLorentzVector>& jets, const TLorentzVector& electron, const TLorentzVector& met);
	//	void reco_Mttbar_btagged(const std::vector<TLorentzVector>& jets, const TLorentzVector& electron, const TLorentzVector& met);
	void reco_Mttbar_matched(const std::vector<TLorentzVector>& jets, const TLorentzVector& electron, const TLorentzVector& met);
	TLorentzVector reconstruct_neutrino(const TLorentzVector&, const TLorentzVector&);
	pair<TLorentzVector, TLorentzVector> reconstruct_neutrinos(const TLorentzVector&, const TLorentzVector&);
	pair<ushort, ushort> reco_hadronic_W(const std::vector<TLorentzVector>& jets, short ommit_blep_id) const;
	short findClosest(const std::vector<TLorentzVector>&, const TLorentzVector&);

	vector<TLorentzVector> GetMCTopEvent();
	void SetHistoLabelCutNjet(TH2D *this_njetVcuts, vector<string>& ve) const;
	void SetHistoLabelEleID(TH1F *eid[]) const;

	void DefineCrossSection();
	double GetCrossSection(const ushort) const;
	void SetWeights();
	//	double GetWeight(const string) const;
	double GetWeight(const short) const;
	long GetNinit(const ushort) const;

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
	double GetChi2Leptonic(double tmass, double angle);
	double GetChi2Hadronic(double Wmass, double tmass, double ptratio);
	double GetChi2Global(double pttbar, double htsystem);
	double GetHT(const std::vector<TLorentzVector>& jets, ushort N);
	bool GetBtagFlag(ushort jetID, ushort btag);
	string ScrNum(double num) const;
	void printLine(ofstream &myfile, const double, const double) const;
	bool ScientificNotation;

	//	bool is_mc_present(const int) const;
	float compute_d0(const string, const int) const;
	float compute_mtw(const TVector2&, const TVector2&) const;
	pair<double, double> compute_neutrino_momentum_z(const TLorentzVector& met, const TLorentzVector& electron);
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
	double weight; // current event weight
	int m_nGoodJet;
	float m_QCDest_reliso_bin_width;
	bool m_doValidation;
	bool m_plotRelisoNES;
	//	bool debug_flag;
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

	// Wjet estimation
	//	TH1D *h_hadTop_maxPT_mass_4j; // 960 bins (0-960)
	//	TH1D *h_hadTop_maxPT_pt_4j;
	//	TH1D *h_hadTop_maxPT_mass_nonIso_4j;
	//	TH1D *h_hadTop_maxPT_pt_nonIso_4j;
	//	TH1D *h_m3_tt;
	//	TH1D *h_m3_wj;
	//	TH1D *h_m3_zj;
	//	TH1D *h_m3_qcd;
	//	TH1D *h_m3_vqq;
	//	TH1D *h_m3_singletop;
	//	TH1D *h_m3_bce[3];
	//	TH1D *h_m3_enri[3];
	//	TH1D *h_m3_tt_control;
	//	TH1D *h_m3_wj_control;
	//	TH1D *h_m3_zj_control;
	//	TH1D *h_m3_qcd_control;
	//	TH1D *h_m3_vqq_control;
	//	TH1D *h_m3_singletop_control;
	//	TH1D *h_m3_bce_control[3];
	//	TH1D *h_m3_enri_control[3];
	// 1000 bins (0-1000)
	//	TH1D *h_hadTop_maxPT_mass_4j_1000;
	//	TH1D *h_hadTop_maxPT_mass_nonIso_4j_1000;
	//	TH1D *h_m3_tt_1000;
	//	TH1D *h_m3_wj_1000;
	//	TH1D *h_m3_zj_1000;
	//	TH1D *h_m3_qcd_1000;
	//	TH1D *h_m3_vqq_1000;
	//	TH1D *h_m3_singletop_1000;
	//	TH1D *h_m3_bce_1000[3];
	//	TH1D *h_m3_enri_1000[3];
	//	TH1D *h_m3_tt_control_1000;
	//	TH1D *h_m3_wj_control_1000;
	//	TH1D *h_m3_zj_control_1000;
	//	TH1D *h_m3_qcd_control_1000;
	//	TH1D *h_m3_vqq_control_1000;
	//	TH1D *h_m3_singletop_control_1000;
	//	TH1D *h_m3_bce_control_1000[3];
	//	TH1D *h_m3_enri_control_1000[3];
	int m_ntoy;

	// MC flag
	bool hasSampleMC(const ushort mc) const;
	bool isMCType(const ushort mc) const;

};

#endif
