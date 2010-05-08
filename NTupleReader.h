/*
 * NTupleReader.h
 *
 *  Created on: Mar 22, 2010
 *      Author: lkreczko
 */

#ifndef NTUPLEREADER_H_
#define NTUPLEREADER_H_
#include "TChain.h"
#include "TBranch.h"

#include <vector>
#include <iostream>

//#ifdef __MAKECINT__
//#pragma link C++ class std::vector<float>+;
//#endif

using namespace std;
/*
 * Class to read the Ntuples from the ROOT file.
 * Includes several flags to minimize the amount of data loaded
 */
class NTupleReader {
public:
	NTupleReader();
	~NTupleReader() {
	}
	;

	/**
	 * @return the debug flag
	 */
	bool GetDebug() const;

	/**
	 * Set the debug flag
	 * @param val new value
	 */
	void SetDebug(bool);

	/**
	 * Get the trigger flag
	 * @return the trigger flag
	 */
	bool GetTrigger() const;

	/**
	 * Set the trigger flag
	 * @param val new value
	 */
	void SetTrigger(bool val);
	void SetTrigger(bool, const string);

	bool GetMissLayersFlag();
	void SetMissLayersFlag(bool);

	bool GetStudyPDFunc();
	void SetStudyPDFunc(bool);

	string GetHLTBit();
	void SetHLTBit(string);

protected:
	// chain = list of root files containing the same tree
	/**
	 * Main TChain holding all variables
	 */
	TChain* chain;
	/**
	 * Secondary TChain holding HLT information.
	 */
	TChain* chain2;

	// Declaration of leaf types
	UInt_t NPFJets;
	vector<float> *PFJets_energy;
	vector<float> *PFJets_et;
	vector<float> *PFJets_eta;
	vector<float> *PFJets_phi;
	vector<float> *PFJets_pt;
	vector<float> *PFJets_px;
	vector<float> *PFJets_py;
	vector<float> *PFJets_pz;
	vector<float> *PFJets_status;
	vector<float> *PFJets_theta;
	vector<float> *PFJets_chgEmE;
	vector<float> *PFJets_chgHadE;
	vector<float> *PFJets_chgMuE;
	vector<float> *PFJets_chg_Mult;
	vector<float> *PFJets_neutralEmE;
	vector<float> *PFJets_neutralHadE;
	vector<float> *PFJets_neutral_Mult;
	vector<float> *PFJets_mu_Mult;
	vector<float> *PFJets_mass;
	vector<float> *PFJets_btag_TC_highEff;
	vector<float> *PFJets_btag_TC_highPur;
	vector<float> *PFJets_btag_jetBProb;
	vector<float> *PFJets_btag_jetProb;
	vector<float> *PFJets_btag_softEle;
	vector<float> *PFJets_btag_softMuon;
	vector<float> *PFJets_btag_softMuonNoIP;
	vector<float> *PFJets_btag_secVertex;
	vector<float> *PFJets_parton_Id;
	UInt_t NPFMets;
	vector<float> *PFMets_et;
	vector<float> *PFMets_phi;
	vector<float> *PFMets_ex;
	vector<float> *PFMets_ey;
	vector<float> *PFMets_sumEt;
	UInt_t NbeamSpot;
	vector<float> *beamSpot_x;
	vector<float> *beamSpot_y;
	vector<float> *beamSpot_z;
	vector<float> *beamSpot_x0Error;
	vector<float> *beamSpot_y0Error;
	vector<float> *beamSpot_z0Error;
	vector<float> *beamSpot_sigmaZ;
	vector<float> *beamSpot_sigmaZ0Error;
	vector<float> *beamSpot_dxdz;
	vector<float> *beamSpot_dxdzError;
	vector<float> *beamSpot_dydz;
	vector<float> *beamSpot_dydzError;
	vector<float> *beamSpot_beamWidthX;
	vector<float> *beamSpot_beamWidthY;
	vector<float> *beamSpot_beamWidthXError;
	vector<float> *beamSpot_beamWidthYError;
	UInt_t Nels;
	vector<float> *els_energy;
	vector<float> *els_et;
	vector<float> *els_eta;
	vector<float> *els_phi;
	vector<float> *els_pt;
	vector<float> *els_px;
	vector<float> *els_py;
	vector<float> *els_pz;
	vector<float> *els_status;
	vector<float> *els_theta;
	vector<float> *els_closestCtfTrackRef;
	vector<float> *els_isEcalDriven;
	vector<float> *els_isTrackerDriven;
	vector<float> *els_dr03EcalRecHitSumEt;
	vector<float> *els_dr04EcalRecHitSumEt;
	vector<float> *els_dr03HcalTowerSumEt;
	vector<float> *els_dr04HcalTowerSumEt;
	vector<float> *els_gen_id;
	vector<float> *els_gen_phi;
	vector<float> *els_gen_pt;
	vector<float> *els_gen_pz;
	vector<float> *els_gen_px;
	vector<float> *els_gen_py;
	vector<float> *els_gen_eta;
	vector<float> *els_gen_theta;
	vector<float> *els_gen_et;
	vector<float> *els_gen_mother_id;
	vector<float> *els_gen_mother_phi;
	vector<float> *els_gen_mother_pt;
	vector<float> *els_gen_mother_pz;
	vector<float> *els_gen_mother_px;
	vector<float> *els_gen_mother_py;
	vector<float> *els_gen_mother_eta;
	vector<float> *els_gen_mother_theta;
	vector<float> *els_gen_mother_et;
	vector<float> *els_tightId;
	vector<float> *els_looseId;
	vector<float> *els_robustTightId;
	vector<float> *els_robustLooseId;
	vector<float> *els_robustHighEnergyId;
	vector<float> *els_cIso;
	vector<float> *els_tIso;
	vector<float> *els_ecalIso;
	vector<float> *els_hcalIso;
	vector<float> *els_chi2;
	vector<float> *els_charge;
	vector<float> *els_caloEnergy;
	vector<float> *els_hadOverEm;
	vector<float> *els_eOverPIn;
	vector<float> *els_eSeedOverPOut;
	vector<float> *els_eSCraw;
	vector<float> *els_eSeed;
	vector<float> *els_sigmaEtaEta;
	vector<float> *els_sigmaIEtaIEta;
	vector<float> *els_scE1x5;
	vector<float> *els_scE2x5Max;
	vector<float> *els_scE5x5;
	vector<float> *els_dEtaIn;
	vector<float> *els_dPhiIn;
	vector<float> *els_dEtaOut;
	vector<float> *els_dPhiOut;
	vector<float> *els_numvalhits;
	vector<float> *els_numlosthits;
	vector<float> *els_basicClustersSize;
	vector<float> *els_tk_pt;
	vector<float> *els_tk_phi;
	vector<float> *els_tk_eta;
	vector<float> *els_tk_theta;
	vector<float> *els_tk_charge;
	vector<float> *els_shFracInnerHits;
	vector<float> *els_d0dum;
	vector<float> *els_dz;
	vector<float> *els_vx;
	vector<float> *els_vy;
	vector<float> *els_vz;
	vector<float> *els_ndof;
	vector<float> *els_ptError;
	vector<float> *els_d0dumError;
	vector<float> *els_dzError;
	vector<float> *els_etaError;
	vector<float> *els_phiError;
	vector<float> *els_cpx;
	vector<float> *els_cpy;
	vector<float> *els_cpz;
	vector<float> *els_vpx;
	vector<float> *els_vpy;
	vector<float> *els_vpz;
	vector<float> *els_cx;
	vector<float> *els_cy;
	vector<float> *els_cz;
	vector<float> *els_innerLayerMissingHits;
	vector<float> *els_isEE;
	vector<float> *els_isEEGap;
	vector<float> *els_isEB;
	vector<float> *els_isEBGap;
	vector<float> *els_isConvertedPhoton;
	UInt_t Njets;
	vector<float> *jets_energy;
	vector<float> *jets_et;
	vector<float> *jets_eta;
	vector<float> *jets_phi;
	vector<float> *jets_pt;
	vector<float> *jets_px;
	vector<float> *jets_py;
	vector<float> *jets_pz;
	vector<float> *jets_status;
	vector<float> *jets_theta;
	vector<float> *jets_parton_Id;
	vector<float> *jets_parton_motherId;
	vector<float> *jets_parton_pt;
	vector<float> *jets_parton_phi;
	vector<float> *jets_parton_eta;
	vector<float> *jets_parton_Energy;
	vector<float> *jets_parton_mass;
	vector<float> *jets_parton_motherID;
	vector<float> *jets_gen_et;
	vector<float> *jets_gen_pt;
	vector<float> *jets_gen_eta;
	vector<float> *jets_gen_phi;
	vector<float> *jets_gen_mass;
	vector<float> *jets_gen_Energy;
	vector<float> *jets_gen_Id;
	vector<float> *jets_gen_motherID;
	vector<float> *jets_gen_threeCharge;
	vector<float> *jets_partonFlavour;
	vector<float> *jets_btag_TC_highPur;
	vector<float> *jets_btag_TC_highEff;
	vector<float> *jets_btag_jetProb;
	vector<float> *jets_btag_jetBProb;
	vector<float> *jets_btag_softEle;
	vector<float> *jets_btag_softMuon;
	vector<float> *jets_btag_softMuonNoIP;
	vector<float> *jets_btag_secVertex;
	vector<float> *jets_chgEmE;
	vector<float> *jets_chgHadE;
	vector<float> *jets_chgMuE;
	vector<float> *jets_chg_Mult;
	vector<float> *jets_neutralEmE;
	vector<float> *jets_neutralHadE;
	vector<float> *jets_neutral_Mult;
	vector<float> *jets_mu_Mult;
	vector<float> *jets_emf;
	vector<float> *jets_ehf;
	vector<float> *jets_n60;
	vector<float> *jets_n90;
	vector<float> *jets_area;
	vector<float> *jets_mass;

	UInt_t NjetsKT4;
	vector<float> *jetsKT4_energy;
	vector<float> *jetsKT4_et;
	vector<float> *jetsKT4_eta;
	vector<float> *jetsKT4_phi;
	vector<float> *jetsKT4_pt;
	vector<float> *jetsKT4_px;
	vector<float> *jetsKT4_py;
	vector<float> *jetsKT4_pz;
	vector<float> *jetsKT4_status;
	vector<float> *jetsKT4_theta;
	vector<float> *jetsKT4_btag_TC_highPur;
	vector<float> *jetsKT4_btag_TC_highEff;
	vector<float> *jetsKT4_btag_jetProb;
	vector<float> *jetsKT4_btag_jetBProb;
	vector<float> *jetsKT4_btag_softEle;
	vector<float> *jetsKT4_btag_softMuon;
	vector<float> *jetsKT4_btag_softMuonNoIP;
	vector<float> *jetsKT4_btag_secVertex;
	vector<float> *jetsKT4_chgEmE;
	vector<float> *jetsKT4_chgHadE;
	vector<float> *jetsKT4_chgMuE;
	vector<float> *jetsKT4_chg_Mult;
	vector<float> *jetsKT4_neutralEmE;
	vector<float> *jetsKT4_neutralHadE;
	vector<float> *jetsKT4_neutral_Mult;
	vector<float> *jetsKT4_mu_Mult;
	vector<float> *jetsKT4_emf;
	vector<float> *jetsKT4_ehf;
	vector<float> *jetsKT4_n60;
	vector<float> *jetsKT4_n90;
	vector<float> *jetsKT4_area;
	vector<float> *jetsKT4_mass;

	UInt_t NjetsSC5;
	vector<float> *jetsSC5_energy;
	vector<float> *jetsSC5_et;
	vector<float> *jetsSC5_eta;
	vector<float> *jetsSC5_phi;
	vector<float> *jetsSC5_pt;
	vector<float> *jetsSC5_px;
	vector<float> *jetsSC5_py;
	vector<float> *jetsSC5_pz;
	vector<float> *jetsSC5_status;
	vector<float> *jetsSC5_theta;
	vector<float> *jetsSC5_btag_TC_highPur;
	vector<float> *jetsSC5_btag_TC_highEff;
	vector<float> *jetsSC5_btag_jetProb;
	vector<float> *jetsSC5_btag_jetBProb;
	vector<float> *jetsSC5_btag_softEle;
	vector<float> *jetsSC5_btag_softMuon;
	vector<float> *jetsSC5_btag_softMuonNoIP;
	vector<float> *jetsSC5_btag_secVertex;
	vector<float> *jetsSC5_chgEmE;
	vector<float> *jetsSC5_chgHadE;
	vector<float> *jetsSC5_chgMuE;
	vector<float> *jetsSC5_chg_Mult;
	vector<float> *jetsSC5_neutralEmE;
	vector<float> *jetsSC5_neutralHadE;
	vector<float> *jetsSC5_neutral_Mult;
	vector<float> *jetsSC5_mu_Mult;
	vector<float> *jetsSC5_emf;
	vector<float> *jetsSC5_ehf;
	vector<float> *jetsSC5_n60;
	vector<float> *jetsSC5_n90;
	vector<float> *jetsSC5_area;
	vector<float> *jetsSC5_mass;
	UInt_t Nmc_doc;
	vector<float> *mc_doc_id;
	vector<float> *mc_doc_pt;
	vector<float> *mc_doc_px;
	vector<float> *mc_doc_py;
	vector<float> *mc_doc_pz;
	vector<float> *mc_doc_eta;
	vector<float> *mc_doc_phi;
	vector<float> *mc_doc_theta;
	vector<float> *mc_doc_energy;
	vector<float> *mc_doc_status;
	vector<float> *mc_doc_charge;
	vector<float> *mc_doc_mother_id;
	vector<float> *mc_doc_grandmother_id;
	vector<float> *mc_doc_ggrandmother_id;
	vector<float> *mc_doc_mother_pt;
	vector<float> *mc_doc_vertex_x;
	vector<float> *mc_doc_vertex_y;
	vector<float> *mc_doc_vertex_z;
	vector<float> *mc_doc_mass;
	vector<float> *mc_doc_numOfDaughters;
	vector<float> *mc_doc_numOfMothers;
	UInt_t Nmets;
	vector<float> *mets_et;
	vector<float> *mets_phi;
	vector<float> *mets_ex;
	vector<float> *mets_ey;
	vector<float> *mets_gen_et;
	vector<float> *mets_gen_phi;
	vector<float> *mets_sign;
	vector<float> *mets_sumEt;
	vector<float> *mets_unCPhi;
	vector<float> *mets_unCPt;
	vector<float> *mets_et_muonCor;
	vector<float> *mets_phi_muonCor;
	vector<float> *mets_et_JESCor;
	vector<float> *mets_phi_JESCor;
	UInt_t NmetsKT4;
	vector<float> *metsKT4_et;
	vector<float> *metsKT4_phi;
	vector<float> *metsKT4_ex;
	vector<float> *metsKT4_ey;
	vector<float> *metsKT4_sumEt;
	vector<float> *metsKT4_et_JESCor;
	vector<float> *metsKT4_phi_JESCor;
	UInt_t NmetsSC5;
	vector<float> *metsSC5_et;
	vector<float> *metsSC5_phi;
	vector<float> *metsSC5_ex;
	vector<float> *metsSC5_ey;
	vector<float> *metsSC5_sumEt;
	vector<float> *metsSC5_et_JESCor;
	vector<float> *metsSC5_phi_JESCor;
	UInt_t Nmus;
	vector<float> *mus_energy;
	vector<float> *mus_et;
	vector<float> *mus_eta;
	vector<float> *mus_phi;
	vector<float> *mus_pt;
	vector<float> *mus_px;
	vector<float> *mus_py;
	vector<float> *mus_pz;
	vector<float> *mus_status;
	vector<float> *mus_theta;
	vector<float> *mus_gen_id;
	vector<float> *mus_gen_phi;
	vector<float> *mus_gen_pt;
	vector<float> *mus_gen_pz;
	vector<float> *mus_gen_px;
	vector<float> *mus_gen_py;
	vector<float> *mus_gen_eta;
	vector<float> *mus_gen_theta;
	vector<float> *mus_gen_et;
	vector<float> *mus_gen_mother_id;
	vector<float> *mus_gen_mother_phi;
	vector<float> *mus_gen_mother_pt;
	vector<float> *mus_gen_mother_pz;
	vector<float> *mus_gen_mother_px;
	vector<float> *mus_gen_mother_py;
	vector<float> *mus_gen_mother_eta;
	vector<float> *mus_gen_mother_theta;
	vector<float> *mus_gen_mother_et;
	vector<float> *mus_tkHits;
	vector<float> *mus_cIso;
	vector<float> *mus_tIso;
	vector<float> *mus_ecalIso;
	vector<float> *mus_hcalIso;
	vector<float> *mus_ecalvetoDep;
	vector<float> *mus_hcalvetoDep;
	vector<float> *mus_calEnergyEm;
	vector<float> *mus_calEnergyHad;
	vector<float> *mus_calEnergyHo;
	vector<float> *mus_calEnergyEmS9;
	vector<float> *mus_calEnergyHadS9;
	vector<float> *mus_calEnergyHoS9;
	vector<float> *mus_iso03_sumPt;
	vector<float> *mus_iso03_emEt;
	vector<float> *mus_iso03_hadEt;
	vector<float> *mus_iso03_hoEt;
	vector<float> *mus_iso03_nTracks;
	vector<float> *mus_iso05_sumPt;
	vector<float> *mus_iso05_emEt;
	vector<float> *mus_iso05_hadEt;
	vector<float> *mus_iso05_hoEt;
	vector<float> *mus_iso05_nTracks;
	vector<float> *mus_charge;
	vector<float> *mus_cm_chi2;
	vector<float> *mus_cm_ndof;
	vector<float> *mus_cm_chg;
	vector<float> *mus_cm_pt;
	vector<float> *mus_cm_px;
	vector<float> *mus_cm_py;
	vector<float> *mus_cm_pz;
	vector<float> *mus_cm_eta;
	vector<float> *mus_cm_phi;
	vector<float> *mus_cm_theta;
	vector<float> *mus_cm_d0dum;
	vector<float> *mus_cm_dz;
	vector<float> *mus_cm_vx;
	vector<float> *mus_cm_vy;
	vector<float> *mus_cm_vz;
	vector<float> *mus_cm_numvalhits;
	vector<float> *mus_cm_numlosthits;
	vector<float> *mus_cm_d0dumErr;
	vector<float> *mus_cm_dzErr;
	vector<float> *mus_cm_ptErr;
	vector<float> *mus_cm_etaErr;
	vector<float> *mus_cm_phiErr;
	vector<float> *mus_tk_chi2;
	vector<float> *mus_tk_ndof;
	vector<float> *mus_tk_chg;
	vector<float> *mus_tk_pt;
	vector<float> *mus_tk_px;
	vector<float> *mus_tk_py;
	vector<float> *mus_tk_pz;
	vector<float> *mus_tk_eta;
	vector<float> *mus_tk_phi;
	vector<float> *mus_tk_theta;
	vector<float> *mus_tk_d0dum;
	vector<float> *mus_tk_dz;
	vector<float> *mus_tk_vx;
	vector<float> *mus_tk_vy;
	vector<float> *mus_tk_vz;
	vector<float> *mus_tk_numvalhits;
	vector<float> *mus_tk_numlosthits;
	vector<float> *mus_tk_d0dumErr;
	vector<float> *mus_tk_dzErr;
	vector<float> *mus_tk_ptErr;
	vector<float> *mus_tk_etaErr;
	vector<float> *mus_tk_phiErr;
	vector<float> *mus_stamu_chi2;
	vector<float> *mus_stamu_ndof;
	vector<float> *mus_stamu_chg;
	vector<float> *mus_stamu_pt;
	vector<float> *mus_stamu_px;
	vector<float> *mus_stamu_py;
	vector<float> *mus_stamu_pz;
	vector<float> *mus_stamu_eta;
	vector<float> *mus_stamu_phi;
	vector<float> *mus_stamu_theta;
	vector<float> *mus_stamu_d0dum;
	vector<float> *mus_stamu_dz;
	vector<float> *mus_stamu_vx;
	vector<float> *mus_stamu_vy;
	vector<float> *mus_stamu_vz;
	vector<float> *mus_stamu_numvalhits;
	vector<float> *mus_stamu_numlosthits;
	vector<float> *mus_stamu_d0dumErr;
	vector<float> *mus_stamu_dzErr;
	vector<float> *mus_stamu_ptErr;
	vector<float> *mus_stamu_etaErr;
	vector<float> *mus_stamu_phiErr;
	vector<float> *mus_num_matches;
	vector<float> *mus_id_All;
	vector<float> *mus_id_AllGlobalMuons;
	vector<float> *mus_id_AllStandAloneMuons;
	vector<float> *mus_id_AllTrackerMuons;
	vector<float> *mus_id_TrackerMuonArbitrated;
	vector<float> *mus_id_AllArbitrated;
	vector<float> *mus_id_GlobalMuonPromptTight;
	vector<float> *mus_id_TMLastStationLoose;
	vector<float> *mus_id_TMLastStationTight;
	vector<float> *mus_id_TM2DCompatibilityLoose;
	vector<float> *mus_id_TM2DCompatibilityTight;
	vector<float> *mus_id_TMOneStationLoose;
	vector<float> *mus_id_TMOneStationTight;
	vector<float> *mus_id_TMLastStationOptimizedLowPtLoose;
	vector<float> *mus_id_TMLastStationOptimizedLowPtTight;
	UInt_t Nphotons;
	vector<float> *photons_energy;
	vector<float> *photons_et;
	vector<float> *photons_eta;
	vector<float> *photons_phi;
	vector<float> *photons_pt;
	vector<float> *photons_px;
	vector<float> *photons_py;
	vector<float> *photons_pz;
	vector<float> *photons_status;
	vector<float> *photons_theta;
	vector<float> *photons_hadOverEM;
	vector<float> *photons_scEnergy;
	vector<float> *photons_scRawEnergy;
	vector<float> *photons_scEta;
	vector<float> *photons_scPhi;
	vector<float> *photons_scEtaWidth;
	vector<float> *photons_scPhiWidth;
	vector<float> *photons_tIso;
	vector<float> *photons_ecalIso;
	vector<float> *photons_hcalIso;
	vector<float> *photons_isoEcalRecHitDR04;
	vector<float> *photons_isoHcalRecHitDR04;
	vector<float> *photons_isoSolidTrkConeDR04;
	vector<float> *photons_isoHollowTrkConeDR04;
	vector<float> *photons_nTrkSolidConeDR04;
	vector<float> *photons_nTrkHollowConeDR04;
	vector<float> *photons_isoEcalRecHitDR03;
	vector<float> *photons_isoHcalRecHitDR03;
	vector<float> *photons_isoSolidTrkConeDR03;
	vector<float> *photons_isoHollowTrkConeDR03;
	vector<float> *photons_nTrkSolidConeDR03;
	vector<float> *photons_nTrkHollowConeDR03;
	vector<float> *photons_isAlsoElectron;
	vector<float> *photons_hasPixelSeed;
	vector<float> *photons_isConverted;
	vector<float> *photons_isEBGap;
	vector<float> *photons_isEEGap;
	vector<float> *photons_isEBEEGap;
	vector<float> *photons_isEBPho;
	vector<float> *photons_isEEPho;
	vector<float> *photons_isLoosePhoton;
	vector<float> *photons_isTightPhoton;
	vector<float> *photons_r9;
	vector<float> *photons_gen_et;
	vector<float> *photons_gen_eta;
	vector<float> *photons_gen_phi;
	vector<float> *photons_gen_id;
	UInt_t Npv;
	vector<float> *pv_x;
	vector<float> *pv_y;
	vector<float> *pv_z;
	vector<float> *pv_xErr;
	vector<float> *pv_yErr;
	vector<float> *pv_zErr;
	UInt_t Ntcmets;
	vector<float> *tcmets_et;
	vector<float> *tcmets_phi;
	vector<float> *tcmets_ex;
	vector<float> *tcmets_ey;
	vector<float> *tcmets_sumEt;
	vector<float> *tcmets_et_muonCor;
	vector<float> *tcmets_phi_muonCor;
	UInt_t Ntracks;
	vector<float> *tracks_chi2;
	vector<float> *tracks_ndof;
	vector<float> *tracks_chg;
	vector<float> *tracks_pt;
	vector<float> *tracks_px;
	vector<float> *tracks_py;
	vector<float> *tracks_pz;
	vector<float> *tracks_eta;
	vector<float> *tracks_phi;
	vector<float> *tracks_theta;
	vector<float> *tracks_d0dum;
	vector<float> *tracks_dz;
	vector<float> *tracks_vx;
	vector<float> *tracks_vy;
	vector<float> *tracks_vz;
	vector<float> *tracks_numvalhits;
	vector<float> *tracks_numlosthits;
	vector<float> *tracks_d0dumErr;
	vector<float> *tracks_dzErr;
	vector<float> *tracks_ptErr;
	vector<float> *tracks_etaErr;
	vector<float> *tracks_phiErr;
	vector<float> *tracks_Nrechits;
	vector<float> *tracks_innerHitX;
	vector<float> *tracks_innerHitY;
	vector<float> *tracks_innerHitZ;
	vector<float> *tracks_outerHitX;
	vector<float> *tracks_outerHitY;
	vector<float> *tracks_outerHitZ;
	vector<float> *tracks_highPurity;
	vector<float> *tracks_innerLayerMissingHits;
	UInt_t run;
	UInt_t event;
	UInt_t lumiBlock;

	// List of branches
	//TBranch        *b_NTEST;   //! TEST
	//TBranch        *b_TEST_x;   //! TEST

	TBranch *b_NPFJets; //!
	TBranch *b_PFJets_energy; //!
	TBranch *b_PFJets_et; //!
	TBranch *b_PFJets_eta; //!
	TBranch *b_PFJets_phi; //!
	TBranch *b_PFJets_pt; //!
	TBranch *b_PFJets_px; //!
	TBranch *b_PFJets_py; //!
	TBranch *b_PFJets_pz; //!
	TBranch *b_PFJets_status; //!
	TBranch *b_PFJets_theta; //!
	TBranch *b_PFJets_chgEmE; //!
	TBranch *b_PFJets_chgHadE; //!
	TBranch *b_PFJets_chgMuE; //!
	TBranch *b_PFJets_chg_Mult; //!
	TBranch *b_PFJets_neutralEmE; //!
	TBranch *b_PFJets_neutralHadE; //!
	TBranch *b_PFJets_neutral_Mult; //!
	TBranch *b_PFJets_mu_Mult; //!
	TBranch *b_PFJets_mass; //!
	TBranch *b_PFJets_btag_TC_highEff;//!
	TBranch *b_PFJets_btag_TC_highPur;//!
	TBranch *b_PFJets_btag_jetBProb;//!
	TBranch *b_PFJets_btag_jetProb;//!
	TBranch *b_PFJets_btag_softEle;//!
	TBranch *b_PFJets_btag_softMuon;//!
	TBranch *b_PFJets_btag_softMuonNoIP;//!
	TBranch *b_PFJets_btag_secVertex;//!
	TBranch *b_PFJets_parton_Id;//!
	TBranch *b_NPFMets; //!
	TBranch *b_PFMets_et; //!
	TBranch *b_PFMets_phi; //!
	TBranch *b_PFMets_ex; //!
	TBranch *b_PFMets_ey; //!
	TBranch *b_PFMets_sumEt; //!
	TBranch *b_NbeamSpot; //!
	TBranch *b_beamSpot_x; //!
	TBranch *b_beamSpot_y; //!
	TBranch *b_beamSpot_z; //!
	TBranch *b_beamSpot_x0Error; //!
	TBranch *b_beamSpot_y0Error; //!
	TBranch *b_beamSpot_z0Error; //!
	TBranch *b_beamSpot_sigmaZ; //!
	TBranch *b_beamSpot_sigmaZ0Error; //!
	TBranch *b_beamSpot_dxdz; //!
	TBranch *b_beamSpot_dxdzError; //!
	TBranch *b_beamSpot_dydz; //!
	TBranch *b_beamSpot_dydzError; //!
	TBranch *b_beamSpot_beamWidthX; //!
	TBranch *b_beamSpot_beamWidthY; //!
	TBranch *b_beamSpot_beamWidthXError; //!
	TBranch *b_beamSpot_beamWidthYError; //!
	TBranch *b_Nels; //!
	TBranch *b_els_energy; //!
	TBranch *b_els_et; //!
	TBranch *b_els_eta; //!
	TBranch *b_els_phi; //!
	TBranch *b_els_pt; //!
	TBranch *b_els_px; //!
	TBranch *b_els_py; //!
	TBranch *b_els_pz; //!
	TBranch *b_els_status; //!
	TBranch *b_els_theta; //!
	TBranch *b_els_closestCtfTrackRef; //!
	TBranch *b_els_isEcalDriven; //!
	TBranch *b_els_isTrackerDriven; //!
	TBranch *b_els_dr03EcalRecHitSumEt; //!
	TBranch *b_els_dr04EcalRecHitSumEt; //!
	TBranch *b_els_dr03HcalTowerSumEt; //!
	TBranch *b_els_dr04HcalTowerSumEt; //!
	TBranch *b_els_gen_id; //!
	TBranch *b_els_gen_phi; //!
	TBranch *b_els_gen_pt; //!
	TBranch *b_els_gen_pz; //!
	TBranch *b_els_gen_px; //!
	TBranch *b_els_gen_py; //!
	TBranch *b_els_gen_eta; //!
	TBranch *b_els_gen_theta; //!
	TBranch *b_els_gen_et; //!
	TBranch *b_els_gen_mother_id; //!
	TBranch *b_els_gen_mother_phi; //!
	TBranch *b_els_gen_mother_pt; //!
	TBranch *b_els_gen_mother_pz; //!
	TBranch *b_els_gen_mother_px; //!
	TBranch *b_els_gen_mother_py; //!
	TBranch *b_els_gen_mother_eta; //!
	TBranch *b_els_gen_mother_theta; //!
	TBranch *b_els_gen_mother_et; //!
	TBranch *b_els_tightId; //!
	TBranch *b_els_looseId; //!
	TBranch *b_els_robustTightId; //!
	TBranch *b_els_robustLooseId; //!
	TBranch *b_els_robustHighEnergyId; //!
	TBranch *b_els_cIso; //!
	TBranch *b_els_tIso; //!
	TBranch *b_els_ecalIso; //!
	TBranch *b_els_hcalIso; //!
	TBranch *b_els_chi2; //!
	TBranch *b_els_charge; //!
	TBranch *b_els_caloEnergy; //!
	TBranch *b_els_hadOverEm; //!
	TBranch *b_els_eOverPIn; //!
	TBranch *b_els_eSeedOverPOut; //!
	TBranch *b_els_eSCraw; //!
	TBranch *b_els_eSeed; //!
	TBranch *b_els_sigmaEtaEta; //!
	TBranch *b_els_sigmaIEtaIEta; //!
	TBranch *b_els_scE1x5; //!
	TBranch *b_els_scE2x5Max; //!
	TBranch *b_els_scE5x5; //!
	TBranch *b_els_dEtaIn; //!
	TBranch *b_els_dPhiIn; //!
	TBranch *b_els_dEtaOut; //!
	TBranch *b_els_dPhiOut; //!
	TBranch *b_els_numvalhits; //!
	TBranch *b_els_numlosthits; //!
	TBranch *b_els_basicClustersSize; //!
	TBranch *b_els_tk_pt; //!
	TBranch *b_els_tk_phi; //!
	TBranch *b_els_tk_eta; //!
	TBranch *b_els_tk_charge; //!
	TBranch *b_els_tk_theta; //!
	TBranch *b_els_shFracInnerHits; //!
	TBranch *b_els_d0dum; //!
	TBranch *b_els_dz; //!
	TBranch *b_els_vx; //!
	TBranch *b_els_vy; //!
	TBranch *b_els_vz; //!
	TBranch *b_els_ndof; //!
	TBranch *b_els_ptError; //!
	TBranch *b_els_d0dumError; //!
	TBranch *b_els_dzError; //!
	TBranch *b_els_etaError; //!
	TBranch *b_els_phiError; //!
	TBranch *b_els_cpx; //!
	TBranch *b_els_cpy; //!
	TBranch *b_els_cpz; //!
	TBranch *b_els_vpx; //!
	TBranch *b_els_vpy; //!
	TBranch *b_els_vpz; //!
	TBranch *b_els_cx; //!
	TBranch *b_els_cy; //!
	TBranch *b_els_cz; //!
	TBranch *b_els_innerLayerMissingHits; //!
	TBranch *b_els_isEE; //!
	TBranch *b_els_isEEGap; //!
	TBranch *b_els_isEB; //!
	TBranch *b_els_isEBGap; //!
	TBranch *b_els_isConvertedPhoton; //!
	TBranch *b_Njets; //!
	TBranch *b_jets_energy; //!
	TBranch *b_jets_et; //!
	TBranch *b_jets_eta; //!
	TBranch *b_jets_phi; //!
	TBranch *b_jets_pt; //!
	TBranch *b_jets_px; //!
	TBranch *b_jets_py; //!
	TBranch *b_jets_pz; //!
	TBranch *b_jets_status; //!
	TBranch *b_jets_theta; //!
	TBranch *b_jets_parton_Id; //!
	TBranch *b_jets_parton_motherId; //!
	TBranch *b_jets_parton_pt; //!
	TBranch *b_jets_parton_phi; //!
	TBranch *b_jets_parton_eta; //!
	TBranch *b_jets_parton_Energy; //!
	TBranch *b_jets_parton_mass; //!
	TBranch *b_jets_parton_motherID; //!
	TBranch *b_jets_gen_et; //!
	TBranch *b_jets_gen_pt; //!
	TBranch *b_jets_gen_eta; //!
	TBranch *b_jets_gen_phi; //!
	TBranch *b_jets_gen_mass; //!
	TBranch *b_jets_gen_Energy; //!
	TBranch *b_jets_gen_Id; //!
	TBranch *b_jets_gen_motherID; //!
	TBranch *b_jets_gen_threeCharge; //!
	TBranch *b_jets_partonFlavour; //!
	TBranch *b_jets_btag_TC_highPur; //!
	TBranch *b_jets_btag_TC_highEff; //!
	TBranch *b_jets_btag_jetProb; //!
	TBranch *b_jets_btag_jetBProb; //!
	TBranch *b_jets_btag_softEle; //!
	TBranch *b_jets_btag_softMuon; //!
	TBranch *b_jets_btag_softMuonNoIP; //!
	TBranch *b_jets_btag_secVertex; //!
	TBranch *b_jets_chgEmE; //!
	TBranch *b_jets_chgHadE; //!
	TBranch *b_jets_chgMuE; //!
	TBranch *b_jets_chg_Mult; //!
	TBranch *b_jets_neutralEmE; //!
	TBranch *b_jets_neutralHadE; //!
	TBranch *b_jets_neutral_Mult; //!
	TBranch *b_jets_mu_Mult; //!
	TBranch *b_jets_emf; //!
	TBranch *b_jets_ehf; //!
	TBranch *b_jets_n60; //!
	TBranch *b_jets_n90; //!
	TBranch *b_jets_area; //!
	TBranch *b_jets_mass; //!

	TBranch *b_NjetsKT4; //!
	TBranch *b_jetsKT4_energy; //!
	TBranch *b_jetsKT4_et; //!
	TBranch *b_jetsKT4_eta; //!
	TBranch *b_jetsKT4_phi; //!
	TBranch *b_jetsKT4_pt; //!
	TBranch *b_jetsKT4_px; //!
	TBranch *b_jetsKT4_py; //!
	TBranch *b_jetsKT4_pz; //!
	TBranch *b_jetsKT4_status; //!
	TBranch *b_jetsKT4_theta; //!
	TBranch *b_jetsKT4_btag_TC_highPur; //!
	TBranch *b_jetsKT4_btag_TC_highEff; //!
	TBranch *b_jetsKT4_btag_jetProb; //!
	TBranch *b_jetsKT4_btag_jetBProb; //!
	TBranch *b_jetsKT4_btag_softEle; //!
	TBranch *b_jetsKT4_btag_softMuon; //!
	TBranch *b_jetsKT4_btag_softMuonNoIP; //!
	TBranch *b_jetsKT4_btag_secVertex; //!
	TBranch *b_jetsKT4_chgEmE; //!
	TBranch *b_jetsKT4_chgHadE; //!
	TBranch *b_jetsKT4_chgMuE; //!
	TBranch *b_jetsKT4_chg_Mult; //!
	TBranch *b_jetsKT4_neutralEmE; //!
	TBranch *b_jetsKT4_neutralHadE; //!
	TBranch *b_jetsKT4_neutral_Mult; //!
	TBranch *b_jetsKT4_mu_Mult; //!
	TBranch *b_jetsKT4_emf; //!
	TBranch *b_jetsKT4_ehf; //!
	TBranch *b_jetsKT4_n60; //!
	TBranch *b_jetsKT4_n90; //!
	TBranch *b_jetsKT4_area; //!
	TBranch *b_jetsKT4_mass; //!

	///3099
	TBranch *b_NjetsSC5; //!
	TBranch *b_jetsSC5_energy; //!
	TBranch *b_jetsSC5_et; //!
	TBranch *b_jetsSC5_eta; //!
	TBranch *b_jetsSC5_phi; //!
	TBranch *b_jetsSC5_pt; //!
	TBranch *b_jetsSC5_px; //!
	TBranch *b_jetsSC5_py; //!
	TBranch *b_jetsSC5_pz; //!
	TBranch *b_jetsSC5_status; //!
	TBranch *b_jetsSC5_theta; //!
	TBranch *b_jetsSC5_btag_TC_highPur; //!
	TBranch *b_jetsSC5_btag_TC_highEff; //!
	TBranch *b_jetsSC5_btag_jetProb; //!
	TBranch *b_jetsSC5_btag_jetBProb; //!
	TBranch *b_jetsSC5_btag_softEle; //!
	TBranch *b_jetsSC5_btag_softMuon; //!
	TBranch *b_jetsSC5_btag_softMuonNoIP; //!
	TBranch *b_jetsSC5_btag_secVertex; //!
	TBranch *b_jetsSC5_chgEmE; //!
	TBranch *b_jetsSC5_chgHadE; //!
	TBranch *b_jetsSC5_chgMuE; //!
	TBranch *b_jetsSC5_chg_Mult; //!
	TBranch *b_jetsSC5_neutralEmE; //!
	TBranch *b_jetsSC5_neutralHadE; //!
	TBranch *b_jetsSC5_neutral_Mult; //!
	TBranch *b_jetsSC5_mu_Mult; //!
	TBranch *b_jetsSC5_emf; //!
	TBranch *b_jetsSC5_ehf; //!
	TBranch *b_jetsSC5_n60; //!
	TBranch *b_jetsSC5_n90; //!
	TBranch *b_jetsSC5_area; //!
	TBranch *b_jetsSC5_mass; //!

	TBranch *b_Nmc_doc; //!
	TBranch *b_mc_doc_id; //!
	TBranch *b_mc_doc_pt; //!
	TBranch *b_mc_doc_px; //!
	TBranch *b_mc_doc_py; //!
	TBranch *b_mc_doc_pz; //!
	TBranch *b_mc_doc_eta; //!
	TBranch *b_mc_doc_phi; //!
	TBranch *b_mc_doc_theta; //!
	TBranch *b_mc_doc_energy; //!
	TBranch *b_mc_doc_status; //!
	TBranch *b_mc_doc_charge; //!
	TBranch *b_mc_doc_mother_id; //!
	TBranch *b_mc_doc_grandmother_id; //!
	TBranch *b_mc_doc_ggrandmother_id; //!
	TBranch *b_mc_doc_mother_pt; //!
	TBranch *b_mc_doc_vertex_x; //!
	TBranch *b_mc_doc_vertex_y; //!
	TBranch *b_mc_doc_vertex_z; //!
	TBranch *b_mc_doc_mass; //!
	TBranch *b_mc_doc_numOfDaughters; //!
	TBranch *b_mc_doc_numOfMothers; //!
	TBranch *b_Nmets; //!
	TBranch *b_mets_et; //!
	TBranch *b_mets_phi; //!
	TBranch *b_mets_ex; //!
	TBranch *b_mets_ey; //!
	TBranch *b_mets_gen_et; //!
	TBranch *b_mets_gen_phi; //!
	TBranch *b_mets_sign; //!
	TBranch *b_mets_sumEt; //!
	TBranch *b_mets_unCPhi; //!
	TBranch *b_mets_unCPt; //!
	TBranch *b_mets_et_muonCor; //!
	TBranch *b_mets_phi_muonCor; //!
	TBranch *b_mets_et_JESCor; //!
	TBranch *b_mets_phi_JESCor; //!

	//3099
	TBranch *b_NmetsKT4; //!
	TBranch *b_metsKT4_et; //!
	TBranch *b_metsKT4_phi; //!
	TBranch *b_metsKT4_ex; //!
	TBranch *b_metsKT4_ey; //!
	TBranch *b_metsKT4_sumEt; //!
	TBranch *b_metsKT4_et_JESCor; //!
	TBranch *b_metsKT4_phi_JESCor; //!
	TBranch *b_NmetsSC5; //!
	TBranch *b_metsSC5_et; //!
	TBranch *b_metsSC5_phi; //!
	TBranch *b_metsSC5_ex; //!
	TBranch *b_metsSC5_ey; //!
	TBranch *b_metsSC5_sign; //!
	TBranch *b_metsSC5_sumEt; //!
	TBranch *b_metsSC5_et_JESCor; //!
	TBranch *b_metsSC5_phi_JESCor; //!

	TBranch *b_Nmus; //!
	TBranch *b_mus_energy; //!
	TBranch *b_mus_et; //!
	TBranch *b_mus_eta; //!
	TBranch *b_mus_phi; //!
	TBranch *b_mus_pt; //!
	TBranch *b_mus_px; //!
	TBranch *b_mus_py; //!
	TBranch *b_mus_pz; //!
	TBranch *b_mus_status; //!
	TBranch *b_mus_theta; //!
	TBranch *b_mus_gen_id; //!
	TBranch *b_mus_gen_phi; //!
	TBranch *b_mus_gen_pt; //!
	TBranch *b_mus_gen_pz; //!
	TBranch *b_mus_gen_px; //!
	TBranch *b_mus_gen_py; //!
	TBranch *b_mus_gen_eta; //!
	TBranch *b_mus_gen_theta; //!
	TBranch *b_mus_gen_et; //!
	TBranch *b_mus_gen_mother_id; //!
	TBranch *b_mus_gen_mother_phi; //!
	TBranch *b_mus_gen_mother_pt; //!
	TBranch *b_mus_gen_mother_pz; //!
	TBranch *b_mus_gen_mother_px; //!
	TBranch *b_mus_gen_mother_py; //!
	TBranch *b_mus_gen_mother_eta; //!
	TBranch *b_mus_gen_mother_theta; //!
	TBranch *b_mus_gen_mother_et; //!
	TBranch *b_mus_tkHits; //!
	TBranch *b_mus_cIso; //!
	TBranch *b_mus_tIso; //!
	TBranch *b_mus_ecalIso; //!
	TBranch *b_mus_hcalIso; //!
	TBranch *b_mus_ecalvetoDep; //!
	TBranch *b_mus_hcalvetoDep; //!
	TBranch *b_mus_calEnergyEm; //!
	TBranch *b_mus_calEnergyHad; //!
	TBranch *b_mus_calEnergyHo; //!
	TBranch *b_mus_calEnergyEmS9; //!
	TBranch *b_mus_calEnergyHadS9; //!
	TBranch *b_mus_calEnergyHoS9; //!
	TBranch *b_mus_iso03_sumPt; //!
	TBranch *b_mus_iso03_emEt; //!
	TBranch *b_mus_iso03_hadEt; //!
	TBranch *b_mus_iso03_hoEt; //!
	TBranch *b_mus_iso03_nTracks; //!
	TBranch *b_mus_iso05_sumPt; //!
	TBranch *b_mus_iso05_emEt; //!
	TBranch *b_mus_iso05_hadEt; //!
	TBranch *b_mus_iso05_hoEt; //!
	TBranch *b_mus_iso05_nTracks; //!
	TBranch *b_mus_charge; //!
	TBranch *b_mus_cm_chi2; //!
	TBranch *b_mus_cm_ndof; //!
	TBranch *b_mus_cm_chg; //!
	TBranch *b_mus_cm_pt; //!
	TBranch *b_mus_cm_px; //!
	TBranch *b_mus_cm_py; //!
	TBranch *b_mus_cm_pz; //!
	TBranch *b_mus_cm_eta; //!
	TBranch *b_mus_cm_phi; //!
	TBranch *b_mus_cm_theta; //!
	TBranch *b_mus_cm_d0dum; //!
	TBranch *b_mus_cm_dz; //!
	TBranch *b_mus_cm_vx; //!
	TBranch *b_mus_cm_vy; //!
	TBranch *b_mus_cm_vz; //!
	TBranch *b_mus_cm_numvalhits; //!
	TBranch *b_mus_cm_numlosthits; //!
	TBranch *b_mus_cm_d0dumErr; //!
	TBranch *b_mus_cm_dzErr; //!
	TBranch *b_mus_cm_ptErr; //!
	TBranch *b_mus_cm_etaErr; //!
	TBranch *b_mus_cm_phiErr; //!
	TBranch *b_mus_tk_chi2; //!
	TBranch *b_mus_tk_ndof; //!
	TBranch *b_mus_tk_chg; //!
	TBranch *b_mus_tk_pt; //!
	TBranch *b_mus_tk_px; //!
	TBranch *b_mus_tk_py; //!
	TBranch *b_mus_tk_pz; //!
	TBranch *b_mus_tk_eta; //!
	TBranch *b_mus_tk_phi; //!
	TBranch *b_mus_tk_theta; //!
	TBranch *b_mus_tk_d0dum; //!
	TBranch *b_mus_tk_dz; //!
	TBranch *b_mus_tk_vx; //!
	TBranch *b_mus_tk_vy; //!
	TBranch *b_mus_tk_vz; //!
	TBranch *b_mus_tk_numvalhits; //!
	TBranch *b_mus_tk_numlosthits; //!
	TBranch *b_mus_tk_d0dumErr; //!
	TBranch *b_mus_tk_dzErr; //!
	TBranch *b_mus_tk_ptErr; //!
	TBranch *b_mus_tk_etaErr; //!
	TBranch *b_mus_tk_phiErr; //!
	TBranch *b_mus_stamu_chi2; //!
	TBranch *b_mus_stamu_ndof; //!
	TBranch *b_mus_stamu_chg; //!
	TBranch *b_mus_stamu_pt; //!
	TBranch *b_mus_stamu_px; //!
	TBranch *b_mus_stamu_py; //!
	TBranch *b_mus_stamu_pz; //!
	TBranch *b_mus_stamu_eta; //!
	TBranch *b_mus_stamu_phi; //!
	TBranch *b_mus_stamu_theta; //!
	TBranch *b_mus_stamu_d0dum; //!
	TBranch *b_mus_stamu_dz; //!
	TBranch *b_mus_stamu_vx; //!
	TBranch *b_mus_stamu_vy; //!
	TBranch *b_mus_stamu_vz; //!
	TBranch *b_mus_stamu_numvalhits; //!
	TBranch *b_mus_stamu_numlosthits; //!
	TBranch *b_mus_stamu_d0dumErr; //!
	TBranch *b_mus_stamu_dzErr; //!
	TBranch *b_mus_stamu_ptErr; //!
	TBranch *b_mus_stamu_etaErr; //!
	TBranch *b_mus_stamu_phiErr; //!
	TBranch *b_mus_num_matches; //!
	TBranch *b_mus_id_All; //!
	TBranch *b_mus_id_AllGlobalMuons; //!
	TBranch *b_mus_id_AllStandAloneMuons; //!
	TBranch *b_mus_id_AllTrackerMuons; //!
	TBranch *b_mus_id_TrackerMuonArbitrated; //!
	TBranch *b_mus_id_AllArbitrated; //!
	TBranch *b_mus_id_GlobalMuonPromptTight; //!
	TBranch *b_mus_id_TMLastStationLoose; //!
	TBranch *b_mus_id_TMLastStationTight; //!
	TBranch *b_mus_id_TM2DCompatibilityLoose; //!
	TBranch *b_mus_id_TM2DCompatibilityTight; //!
	TBranch *b_mus_id_TMOneStationLoose; //!
	TBranch *b_mus_id_TMOneStationTight; //!
	TBranch *b_mus_id_TMLastStationOptimizedLowPtLoose; //!
	TBranch *b_mus_id_TMLastStationOptimizedLowPtTight; //!
	TBranch *b_Nphotons; //!
	TBranch *b_photons_energy; //!
	TBranch *b_photons_et; //!
	TBranch *b_photons_eta; //!
	TBranch *b_photons_phi; //!
	TBranch *b_photons_pt; //!
	TBranch *b_photons_px; //!
	TBranch *b_photons_py; //!
	TBranch *b_photons_pz; //!
	TBranch *b_photons_status; //!
	TBranch *b_photons_theta; //!
	TBranch *b_photons_hadOverEM; //!
	TBranch *b_photons_scEnergy; //!
	TBranch *b_photons_scRawEnergy; //!
	TBranch *b_photons_scEta; //!
	TBranch *b_photons_scPhi; //!
	TBranch *b_photons_scEtaWidth; //!
	TBranch *b_photons_scPhiWidth; //!
	TBranch *b_photons_tIso; //!
	TBranch *b_photons_ecalIso; //!
	TBranch *b_photons_hcalIso; //!
	TBranch *b_photons_isoEcalRecHitDR04; //!
	TBranch *b_photons_isoHcalRecHitDR04; //!
	TBranch *b_photons_isoSolidTrkConeDR04; //!
	TBranch *b_photons_isoHollowTrkConeDR04; //!
	TBranch *b_photons_nTrkSolidConeDR04; //!
	TBranch *b_photons_nTrkHollowConeDR04; //!
	TBranch *b_photons_isoEcalRecHitDR03; //!
	TBranch *b_photons_isoHcalRecHitDR03; //!
	TBranch *b_photons_isoSolidTrkConeDR03; //!
	TBranch *b_photons_isoHollowTrkConeDR03; //!
	TBranch *b_photons_nTrkSolidConeDR03; //!
	TBranch *b_photons_nTrkHollowConeDR03; //!
	TBranch *b_photons_isAlsoElectron; //!
	TBranch *b_photons_hasPixelSeed; //!
	TBranch *b_photons_isConverted; //!
	TBranch *b_photons_isEBGap; //!
	TBranch *b_photons_isEEGap; //!
	TBranch *b_photons_isEBEEGap; //!
	TBranch *b_photons_isEBPho; //!
	TBranch *b_photons_isEEPho; //!
	TBranch *b_photons_isLoosePhoton; //!
	TBranch *b_photons_isTightPhoton; //!
	TBranch *b_photons_r9; //!
	TBranch *b_photons_gen_et; //!
	TBranch *b_photons_gen_eta; //!
	TBranch *b_photons_gen_phi; //!
	TBranch *b_photons_gen_id; //!
	TBranch *b_Npv; //!
	TBranch *b_pv_x; //!
	TBranch *b_pv_y; //!
	TBranch *b_pv_z; //!
	TBranch *b_pv_xErr; //!
	TBranch *b_pv_yErr; //!
	TBranch *b_pv_zErr; //!
	TBranch *b_Ntcmets; //!
	TBranch *b_tcmets_et; //!
	TBranch *b_tcmets_phi; //!
	TBranch *b_tcmets_ex; //!
	TBranch *b_tcmets_ey; //!
	TBranch *b_tcmets_sumEt; //!
	TBranch *b_tcmets_et_muonCor; //!
	TBranch *b_tcmets_phi_muonCor; //!
	TBranch *b_Ntracks; //!
	TBranch *b_tracks_chi2; //!
	TBranch *b_tracks_ndof; //!
	TBranch *b_tracks_chg; //!
	TBranch *b_tracks_pt; //!
	TBranch *b_tracks_px; //!
	TBranch *b_tracks_py; //!
	TBranch *b_tracks_pz; //!
	TBranch *b_tracks_eta; //!
	TBranch *b_tracks_phi; //!
	TBranch *b_tracks_theta; //!
	TBranch *b_tracks_d0dum; //!
	TBranch *b_tracks_dz; //!
	TBranch *b_tracks_vx; //!
	TBranch *b_tracks_vy; //!
	TBranch *b_tracks_vz; //!
	TBranch *b_tracks_numvalhits; //!
	TBranch *b_tracks_numlosthits; //!
	TBranch *b_tracks_d0dumErr; //!
	TBranch *b_tracks_dzErr; //!
	TBranch *b_tracks_ptErr; //!
	TBranch *b_tracks_etaErr; //!
	TBranch *b_tracks_phiErr; //!
	TBranch *b_tracks_Nrechits; //!
	TBranch *b_tracks_innerHitX; //!
	TBranch *b_tracks_innerHitY; //!
	TBranch *b_tracks_innerHitZ; //!
	TBranch *b_tracks_outerHitX; //!
	TBranch *b_tracks_outerHitY; //!
	TBranch *b_tracks_outerHitZ; //!
	TBranch *b_tracks_highPurity; //!
	TBranch *b_tracks_innerLayerMissingHits; //!
	TBranch *b_run; //!
	TBranch *b_event; //!
	TBranch *b_lumiBlock; //!


	//-----------------------
	// HLT electron trigger
	//-----------------------
	// For 8E29 (ref: https://twiki.cern.ch/twiki/bin/view/CMS/TSG_13_VIII_09_8E29)
	// For 1E31 (ref: https://twiki.cern.ch/twiki/bin/view/CMS/TSG_13_VIII_09_1E31)

	Double_t HLT_Ele10_LW_EleId_L1R;
	Double_t HLT_Ele10_SW_L1R;
	Double_t HLT_Ele15_LW_L1R;
	Double_t HLT_Ele15_SC10_LW_L1R;
	Double_t HLT_Ele15_SW_L1R;
	Double_t HLT_Ele15_SiStrip_L1R;
	Double_t HLT_Ele20_LW_L1R;

	TBranch *b_HLT_Ele10_LW_EleId_L1R; //!
	TBranch *b_HLT_Ele10_SW_L1R; //!
	TBranch *b_HLT_Ele15_LW_L1R; //!
	TBranch *b_HLT_Ele15_SC10_LW_L1R; //!
	TBranch *b_HLT_Ele15_SW_L1R; //!
	TBranch *b_HLT_Ele15_SiStrip_L1R; //!
	TBranch *b_HLT_Ele20_LW_L1R; //!

	// PDF Weights
	//------------
	Double_t PDFWcteq66_0;
	Double_t PDFWcteq66_1;
	Double_t PDFWcteq66_2;
	Double_t PDFWcteq66_3;
	Double_t PDFWcteq66_4;
	Double_t PDFWcteq66_5;
	Double_t PDFWcteq66_6;
	Double_t PDFWcteq66_7;
	Double_t PDFWcteq66_8;
	Double_t PDFWcteq66_9;
	Double_t PDFWcteq66_10;
	Double_t PDFWcteq66_11;
	Double_t PDFWcteq66_12;
	Double_t PDFWcteq66_13;
	Double_t PDFWcteq66_14;
	Double_t PDFWcteq66_15;
	Double_t PDFWcteq66_16;
	Double_t PDFWcteq66_17;
	Double_t PDFWcteq66_18;
	Double_t PDFWcteq66_19;
	Double_t PDFWcteq66_20;
	Double_t PDFWcteq66_21;
	Double_t PDFWcteq66_22;
	Double_t PDFWcteq66_23;
	Double_t PDFWcteq66_24;
	Double_t PDFWcteq66_25;
	Double_t PDFWcteq66_26;
	Double_t PDFWcteq66_27;
	Double_t PDFWcteq66_28;
	Double_t PDFWcteq66_29;
	Double_t PDFWcteq66_30;
	Double_t PDFWcteq66_31;
	Double_t PDFWcteq66_32;
	Double_t PDFWcteq66_33;
	Double_t PDFWcteq66_34;
	Double_t PDFWcteq66_35;
	Double_t PDFWcteq66_36;
	Double_t PDFWcteq66_37;
	Double_t PDFWcteq66_38;
	Double_t PDFWcteq66_39;
	Double_t PDFWcteq66_40;
	Double_t PDFWcteq66_41;
	Double_t PDFWcteq66_42;
	Double_t PDFWcteq66_43;
	Double_t PDFWcteq66_44;

	TBranch *b_PDFWcteq66_0; //!
	TBranch *b_PDFWcteq66_1; //!
	TBranch *b_PDFWcteq66_2; //!
	TBranch *b_PDFWcteq66_3; //!
	TBranch *b_PDFWcteq66_4; //!
	TBranch *b_PDFWcteq66_5; //!
	TBranch *b_PDFWcteq66_6; //!
	TBranch *b_PDFWcteq66_7; //!
	TBranch *b_PDFWcteq66_8; //!
	TBranch *b_PDFWcteq66_9; //!
	TBranch *b_PDFWcteq66_10; //!
	TBranch *b_PDFWcteq66_11; //!
	TBranch *b_PDFWcteq66_12; //!
	TBranch *b_PDFWcteq66_13; //!
	TBranch *b_PDFWcteq66_14; //!
	TBranch *b_PDFWcteq66_15; //!
	TBranch *b_PDFWcteq66_16; //!
	TBranch *b_PDFWcteq66_17; //!
	TBranch *b_PDFWcteq66_18; //!
	TBranch *b_PDFWcteq66_19; //!
	TBranch *b_PDFWcteq66_20; //!
	TBranch *b_PDFWcteq66_21; //!
	TBranch *b_PDFWcteq66_22; //!
	TBranch *b_PDFWcteq66_23; //!
	TBranch *b_PDFWcteq66_24; //!
	TBranch *b_PDFWcteq66_25; //!
	TBranch *b_PDFWcteq66_26; //!
	TBranch *b_PDFWcteq66_27; //!
	TBranch *b_PDFWcteq66_28; //!
	TBranch *b_PDFWcteq66_29; //!
	TBranch *b_PDFWcteq66_30; //!
	TBranch *b_PDFWcteq66_31; //!
	TBranch *b_PDFWcteq66_32; //!
	TBranch *b_PDFWcteq66_33; //!
	TBranch *b_PDFWcteq66_34; //!
	TBranch *b_PDFWcteq66_35; //!
	TBranch *b_PDFWcteq66_36; //!
	TBranch *b_PDFWcteq66_37; //!
	TBranch *b_PDFWcteq66_38; //!
	TBranch *b_PDFWcteq66_39; //!
	TBranch *b_PDFWcteq66_40; //!
	TBranch *b_PDFWcteq66_41; //!
	TBranch *b_PDFWcteq66_42; //!
	TBranch *b_PDFWcteq66_43; //!
	TBranch *b_PDFWcteq66_44; //!

	bool datafile, m_useMisslayers, checkTrig, m_studyPDFunc, debug_flag;
	string m_jetAlgo, m_metAlgo, HLTBit;

	void Init(); //initialize tree branches
	void ReadSelectedBranches() const;
	bool IsData() const;
};
#endif /* NTUPLEREADER_H_ */

/*
 * NTupleReader.cc
 *
 *  Created on: Mar 22, 2010
 *      Author: lkreczko
 */

