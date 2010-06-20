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

#ifdef __MAKECINT__
#pragma link C++ class std::vector<float>+;
#endif

/*
 * Class to read the Ntuples from the ROOT file.
 * Includes several flags to minimize the amount of data loaded
 */
class NTupleReader {
public:
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
	std::vector<float> *PFJets_energy;
	std::vector<float> *PFJets_et;
	std::vector<float> *PFJets_eta;
	std::vector<float> *PFJets_phi;
	std::vector<float> *PFJets_pt;
	std::vector<float> *PFJets_px;
	std::vector<float> *PFJets_py;
	std::vector<float> *PFJets_pz;
	std::vector<float> *PFJets_status;
	std::vector<float> *PFJets_theta;
	std::vector<float> *PFJets_chgEmE;
	std::vector<float> *PFJets_chgHadE;
	std::vector<float> *PFJets_chgMuE;
	std::vector<float> *PFJets_chg_Mult;
	std::vector<float> *PFJets_neutralEmE;
	std::vector<float> *PFJets_neutralHadE;
	std::vector<float> *PFJets_neutral_Mult;
	std::vector<float> *PFJets_mu_Mult;
	std::vector<float> *PFJets_mass;
	std::vector<float> *PFJets_btag_TC_highEff;
	std::vector<float> *PFJets_btag_TC_highPur;
	std::vector<float> *PFJets_btag_jetBProb;
	std::vector<float> *PFJets_btag_jetProb;
//	std::vector<float> *PFJets_btag_softEle;
//	std::vector<float> *PFJets_btag_softMuon;
//	std::vector<float> *PFJets_btag_softMuonNoIP;
	std::vector<float> *PFJets_btag_secVertex;
	std::vector<float> *PFJets_parton_Id;
	UInt_t NPFMets;
	std::vector<float> *PFMets_et;
	std::vector<float> *PFMets_phi;
	std::vector<float> *PFMets_ex;
	std::vector<float> *PFMets_ey;
	std::vector<float> *PFMets_sumEt;
	UInt_t NbeamSpot;
	std::vector<float> *beamSpot_x;
	std::vector<float> *beamSpot_y;
	std::vector<float> *beamSpot_z;
	std::vector<float> *beamSpot_x0Error;
	std::vector<float> *beamSpot_y0Error;
	std::vector<float> *beamSpot_z0Error;
	std::vector<float> *beamSpot_sigmaZ;
	std::vector<float> *beamSpot_sigmaZ0Error;
	std::vector<float> *beamSpot_dxdz;
	std::vector<float> *beamSpot_dxdzError;
	std::vector<float> *beamSpot_dydz;
	std::vector<float> *beamSpot_dydzError;
	std::vector<float> *beamSpot_beamWidthX;
	std::vector<float> *beamSpot_beamWidthY;
	std::vector<float> *beamSpot_beamWidthXError;
	std::vector<float> *beamSpot_beamWidthYError;
	UInt_t Nels;
	std::vector<float> *els_energy;
	std::vector<float> *els_et;
	std::vector<float> *els_eta;
	std::vector<float> *els_phi;
	std::vector<float> *els_pt;
	std::vector<float> *els_px;
	std::vector<float> *els_py;
	std::vector<float> *els_pz;
	std::vector<float> *els_status;
	std::vector<float> *els_theta;
	std::vector<float> *els_closestCtfTrackRef;
	std::vector<float> *els_isEcalDriven;
	std::vector<float> *els_isTrackerDriven;
	std::vector<float> *els_dr03EcalRecHitSumEt;
	std::vector<float> *els_dr04EcalRecHitSumEt;
	std::vector<float> *els_dr03HcalTowerSumEt;
	std::vector<float> *els_dr04HcalTowerSumEt;
	std::vector<float> *els_gen_id;
	std::vector<float> *els_gen_phi;
	std::vector<float> *els_gen_pt;
	std::vector<float> *els_gen_pz;
	std::vector<float> *els_gen_px;
	std::vector<float> *els_gen_py;
	std::vector<float> *els_gen_eta;
	std::vector<float> *els_gen_theta;
	std::vector<float> *els_gen_et;
	std::vector<float> *els_gen_mother_id;
	std::vector<float> *els_gen_mother_phi;
	std::vector<float> *els_gen_mother_pt;
	std::vector<float> *els_gen_mother_pz;
	std::vector<float> *els_gen_mother_px;
	std::vector<float> *els_gen_mother_py;
	std::vector<float> *els_gen_mother_eta;
	std::vector<float> *els_gen_mother_theta;
	std::vector<float> *els_gen_mother_et;
	std::vector<float> *els_tightId;
	std::vector<float> *els_looseId;
	std::vector<float> *els_robustTightId;
	std::vector<float> *els_robustLooseId;
	std::vector<float> *els_robustHighEnergyId;
	std::vector<float> *els_cIso;
	std::vector<float> *els_tIso;
	std::vector<float> *els_ecalIso;
	std::vector<float> *els_hcalIso;
	std::vector<float> *els_chi2;
	std::vector<float> *els_charge;
	std::vector<float> *els_caloEnergy;
	std::vector<float> *els_hadOverEm;
	std::vector<float> *els_eOverPIn;
	std::vector<float> *els_eSeedOverPOut;
	std::vector<float> *els_eSCraw;
	std::vector<float> *els_eSeed;
	std::vector<float> *els_sigmaEtaEta;
	std::vector<float> *els_sigmaIEtaIEta;
	std::vector<float> *els_scE1x5;
	std::vector<float> *els_scE2x5Max;
	std::vector<float> *els_scE5x5;
	std::vector<float> *els_dEtaIn;
	std::vector<float> *els_dPhiIn;
	std::vector<float> *els_dEtaOut;
	std::vector<float> *els_dPhiOut;
	std::vector<float> *els_numvalhits;
	std::vector<float> *els_numlosthits;
	std::vector<float> *els_basicClustersSize;
	std::vector<float> *els_tk_pt;
	std::vector<float> *els_tk_phi;
	std::vector<float> *els_tk_eta;
	std::vector<float> *els_tk_theta;
	std::vector<float> *els_tk_charge;
	std::vector<float> *els_shFracInnerHits;
	std::vector<float> *els_d0dum;
	std::vector<float> *els_dz;
	std::vector<float> *els_vx;
	std::vector<float> *els_vy;
	std::vector<float> *els_vz;
	std::vector<float> *els_ndof;
	std::vector<float> *els_ptError;
	std::vector<float> *els_d0dumError;
	std::vector<float> *els_dzError;
	std::vector<float> *els_etaError;
	std::vector<float> *els_phiError;
	std::vector<float> *els_cpx;
	std::vector<float> *els_cpy;
	std::vector<float> *els_cpz;
	std::vector<float> *els_vpx;
	std::vector<float> *els_vpy;
	std::vector<float> *els_vpz;
	std::vector<float> *els_cx;
	std::vector<float> *els_cy;
	std::vector<float> *els_cz;
	std::vector<float> *els_innerLayerMissingHits;
	std::vector<float> *els_isEE;
	std::vector<float> *els_isEEGap;
	std::vector<float> *els_isEB;
	std::vector<float> *els_isEBGap;
	std::vector<float> *els_isConvertedPhoton;
	UInt_t Njets;
	std::vector<float> *jets_energy;
	std::vector<float> *jets_et;
	std::vector<float> *jets_eta;
	std::vector<float> *jets_phi;
	std::vector<float> *jets_pt;
	std::vector<float> *jets_px;
	std::vector<float> *jets_py;
	std::vector<float> *jets_pz;
	std::vector<float> *jets_status;
	std::vector<float> *jets_theta;
	std::vector<float> *jets_parton_Id;
	std::vector<float> *jets_parton_motherId;
	std::vector<float> *jets_parton_pt;
	std::vector<float> *jets_parton_phi;
	std::vector<float> *jets_parton_eta;
	std::vector<float> *jets_parton_Energy;
	std::vector<float> *jets_parton_mass;
	std::vector<float> *jets_parton_motherID;
	std::vector<float> *jets_gen_et;
	std::vector<float> *jets_gen_pt;
	std::vector<float> *jets_gen_eta;
	std::vector<float> *jets_gen_phi;
	std::vector<float> *jets_gen_mass;
	std::vector<float> *jets_gen_Energy;
	std::vector<float> *jets_gen_Id;
	std::vector<float> *jets_gen_motherID;
	std::vector<float> *jets_gen_threeCharge;
	std::vector<float> *jets_partonFlavour;
	std::vector<float> *jets_btag_TC_highPur;
	std::vector<float> *jets_btag_TC_highEff;
	std::vector<float> *jets_btag_jetProb;
	std::vector<float> *jets_btag_jetBProb;
//	std::vector<float> *jets_btag_softEle;
//	std::vector<float> *jets_btag_softMuon;
//	std::vector<float> *jets_btag_softMuonNoIP;
	std::vector<float> *jets_btag_secVertex;
	std::vector<float> *jets_chgEmE;
	std::vector<float> *jets_chgHadE;
	std::vector<float> *jets_chgMuE;
	std::vector<float> *jets_chg_Mult;
	std::vector<float> *jets_neutralEmE;
	std::vector<float> *jets_neutralHadE;
	std::vector<float> *jets_neutral_Mult;
	std::vector<float> *jets_mu_Mult;
	std::vector<float> *jets_emf;
	std::vector<float> *jets_ehf;
	std::vector<float> *jets_n60;
	std::vector<float> *jets_n90;
	std::vector<float> *jets_area;
	std::vector<float> *jets_mass;

	UInt_t NjetsKT4;
	std::vector<float> *jetsKT4_energy;
	std::vector<float> *jetsKT4_et;
	std::vector<float> *jetsKT4_eta;
	std::vector<float> *jetsKT4_phi;
	std::vector<float> *jetsKT4_pt;
	std::vector<float> *jetsKT4_px;
	std::vector<float> *jetsKT4_py;
	std::vector<float> *jetsKT4_pz;
	std::vector<float> *jetsKT4_status;
	std::vector<float> *jetsKT4_theta;
	std::vector<float> *jetsKT4_btag_TC_highPur;
	std::vector<float> *jetsKT4_btag_TC_highEff;
	std::vector<float> *jetsKT4_btag_jetProb;
	std::vector<float> *jetsKT4_btag_jetBProb;
//	std::vector<float> *jetsKT4_btag_softEle;
//	std::vector<float> *jetsKT4_btag_softMuon;
//	std::vector<float> *jetsKT4_btag_softMuonNoIP;
	std::vector<float> *jetsKT4_btag_secVertex;
	std::vector<float> *jetsKT4_chgEmE;
	std::vector<float> *jetsKT4_chgHadE;
	std::vector<float> *jetsKT4_chgMuE;
	std::vector<float> *jetsKT4_chg_Mult;
	std::vector<float> *jetsKT4_neutralEmE;
	std::vector<float> *jetsKT4_neutralHadE;
	std::vector<float> *jetsKT4_neutral_Mult;
	std::vector<float> *jetsKT4_mu_Mult;
	std::vector<float> *jetsKT4_emf;
	std::vector<float> *jetsKT4_ehf;
	std::vector<float> *jetsKT4_n60;
	std::vector<float> *jetsKT4_n90;
	std::vector<float> *jetsKT4_area;
	std::vector<float> *jetsKT4_mass;

	UInt_t NjetsSC5;
	std::vector<float> *jetsSC5_energy;
	std::vector<float> *jetsSC5_et;
	std::vector<float> *jetsSC5_eta;
	std::vector<float> *jetsSC5_phi;
	std::vector<float> *jetsSC5_pt;
	std::vector<float> *jetsSC5_px;
	std::vector<float> *jetsSC5_py;
	std::vector<float> *jetsSC5_pz;
	std::vector<float> *jetsSC5_status;
	std::vector<float> *jetsSC5_theta;
	std::vector<float> *jetsSC5_btag_TC_highPur;
	std::vector<float> *jetsSC5_btag_TC_highEff;
	std::vector<float> *jetsSC5_btag_jetProb;
	std::vector<float> *jetsSC5_btag_jetBProb;
//	std::vector<float> *jetsSC5_btag_softEle;
//	std::vector<float> *jetsSC5_btag_softMuon;
//	std::vector<float> *jetsSC5_btag_softMuonNoIP;
	std::vector<float> *jetsSC5_btag_secVertex;
	std::vector<float> *jetsSC5_chgEmE;
	std::vector<float> *jetsSC5_chgHadE;
	std::vector<float> *jetsSC5_chgMuE;
	std::vector<float> *jetsSC5_chg_Mult;
	std::vector<float> *jetsSC5_neutralEmE;
	std::vector<float> *jetsSC5_neutralHadE;
	std::vector<float> *jetsSC5_neutral_Mult;
	std::vector<float> *jetsSC5_mu_Mult;
	std::vector<float> *jetsSC5_emf;
	std::vector<float> *jetsSC5_ehf;
	std::vector<float> *jetsSC5_n60;
	std::vector<float> *jetsSC5_n90;
	std::vector<float> *jetsSC5_area;
	std::vector<float> *jetsSC5_mass;
	UInt_t Nmc_doc;
	std::vector<float> *mc_doc_id;
	std::vector<float> *mc_doc_pt;
	std::vector<float> *mc_doc_px;
	std::vector<float> *mc_doc_py;
	std::vector<float> *mc_doc_pz;
	std::vector<float> *mc_doc_eta;
	std::vector<float> *mc_doc_phi;
	std::vector<float> *mc_doc_theta;
	std::vector<float> *mc_doc_energy;
	std::vector<float> *mc_doc_status;
	std::vector<float> *mc_doc_charge;
	std::vector<float> *mc_doc_mother_id;
	std::vector<float> *mc_doc_grandmother_id;
	std::vector<float> *mc_doc_ggrandmother_id;
	std::vector<float> *mc_doc_mother_pt;
	std::vector<float> *mc_doc_vertex_x;
	std::vector<float> *mc_doc_vertex_y;
	std::vector<float> *mc_doc_vertex_z;
	std::vector<float> *mc_doc_mass;
	std::vector<float> *mc_doc_numOfDaughters;
	std::vector<float> *mc_doc_numOfMothers;
	UInt_t Nmets;
	std::vector<float> *mets_et;
	std::vector<float> *mets_phi;
	std::vector<float> *mets_ex;
	std::vector<float> *mets_ey;
	std::vector<float> *mets_gen_et;
	std::vector<float> *mets_gen_phi;
	std::vector<float> *mets_sign;
	std::vector<float> *mets_sumEt;
	std::vector<float> *mets_unCPhi;
	std::vector<float> *mets_unCPt;
	std::vector<float> *mets_et_muonCor;
	std::vector<float> *mets_phi_muonCor;
	std::vector<float> *mets_et_JESCor;
	std::vector<float> *mets_phi_JESCor;
	UInt_t NmetsKT4;
	std::vector<float> *metsKT4_et;
	std::vector<float> *metsKT4_phi;
	std::vector<float> *metsKT4_ex;
	std::vector<float> *metsKT4_ey;
	std::vector<float> *metsKT4_sumEt;
	std::vector<float> *metsKT4_et_JESCor;
	std::vector<float> *metsKT4_phi_JESCor;
	UInt_t NmetsSC5;
	std::vector<float> *metsSC5_et;
	std::vector<float> *metsSC5_phi;
	std::vector<float> *metsSC5_ex;
	std::vector<float> *metsSC5_ey;
	std::vector<float> *metsSC5_sumEt;
	std::vector<float> *metsSC5_et_JESCor;
	std::vector<float> *metsSC5_phi_JESCor;
	UInt_t Nmus;
	std::vector<float> *mus_energy;
	std::vector<float> *mus_et;
	std::vector<float> *mus_eta;
	std::vector<float> *mus_phi;
	std::vector<float> *mus_pt;
	std::vector<float> *mus_px;
	std::vector<float> *mus_py;
	std::vector<float> *mus_pz;
	std::vector<float> *mus_status;
	std::vector<float> *mus_theta;
	std::vector<float> *mus_gen_id;
	std::vector<float> *mus_gen_phi;
	std::vector<float> *mus_gen_pt;
	std::vector<float> *mus_gen_pz;
	std::vector<float> *mus_gen_px;
	std::vector<float> *mus_gen_py;
	std::vector<float> *mus_gen_eta;
	std::vector<float> *mus_gen_theta;
	std::vector<float> *mus_gen_et;
	std::vector<float> *mus_gen_mother_id;
	std::vector<float> *mus_gen_mother_phi;
	std::vector<float> *mus_gen_mother_pt;
	std::vector<float> *mus_gen_mother_pz;
	std::vector<float> *mus_gen_mother_px;
	std::vector<float> *mus_gen_mother_py;
	std::vector<float> *mus_gen_mother_eta;
	std::vector<float> *mus_gen_mother_theta;
	std::vector<float> *mus_gen_mother_et;
	std::vector<float> *mus_tkHits;
	std::vector<float> *mus_cIso;
	std::vector<float> *mus_tIso;
	std::vector<float> *mus_ecalIso;
	std::vector<float> *mus_hcalIso;
	std::vector<float> *mus_ecalvetoDep;
	std::vector<float> *mus_hcalvetoDep;
	std::vector<float> *mus_calEnergyEm;
	std::vector<float> *mus_calEnergyHad;
	std::vector<float> *mus_calEnergyHo;
	std::vector<float> *mus_calEnergyEmS9;
	std::vector<float> *mus_calEnergyHadS9;
	std::vector<float> *mus_calEnergyHoS9;
	std::vector<float> *mus_iso03_sumPt;
	std::vector<float> *mus_iso03_emEt;
	std::vector<float> *mus_iso03_hadEt;
	std::vector<float> *mus_iso03_hoEt;
	std::vector<float> *mus_iso03_nTracks;
	std::vector<float> *mus_iso05_sumPt;
	std::vector<float> *mus_iso05_emEt;
	std::vector<float> *mus_iso05_hadEt;
	std::vector<float> *mus_iso05_hoEt;
	std::vector<float> *mus_iso05_nTracks;
	std::vector<float> *mus_charge;
	std::vector<float> *mus_cm_chi2;
	std::vector<float> *mus_cm_ndof;
	std::vector<float> *mus_cm_chg;
	std::vector<float> *mus_cm_pt;
	std::vector<float> *mus_cm_px;
	std::vector<float> *mus_cm_py;
	std::vector<float> *mus_cm_pz;
	std::vector<float> *mus_cm_eta;
	std::vector<float> *mus_cm_phi;
	std::vector<float> *mus_cm_theta;
	std::vector<float> *mus_cm_d0dum;
	std::vector<float> *mus_cm_dz;
	std::vector<float> *mus_cm_vx;
	std::vector<float> *mus_cm_vy;
	std::vector<float> *mus_cm_vz;
	std::vector<float> *mus_cm_numvalhits;
	std::vector<float> *mus_cm_numlosthits;
	std::vector<float> *mus_cm_d0dumErr;
	std::vector<float> *mus_cm_dzErr;
	std::vector<float> *mus_cm_ptErr;
	std::vector<float> *mus_cm_etaErr;
	std::vector<float> *mus_cm_phiErr;
	std::vector<float> *mus_tk_chi2;
	std::vector<float> *mus_tk_ndof;
	std::vector<float> *mus_tk_chg;
	std::vector<float> *mus_tk_pt;
	std::vector<float> *mus_tk_px;
	std::vector<float> *mus_tk_py;
	std::vector<float> *mus_tk_pz;
	std::vector<float> *mus_tk_eta;
	std::vector<float> *mus_tk_phi;
	std::vector<float> *mus_tk_theta;
	std::vector<float> *mus_tk_d0dum;
	std::vector<float> *mus_tk_dz;
	std::vector<float> *mus_tk_vx;
	std::vector<float> *mus_tk_vy;
	std::vector<float> *mus_tk_vz;
	std::vector<float> *mus_tk_numvalhits;
	std::vector<float> *mus_tk_numlosthits;
	std::vector<float> *mus_tk_d0dumErr;
	std::vector<float> *mus_tk_dzErr;
	std::vector<float> *mus_tk_ptErr;
	std::vector<float> *mus_tk_etaErr;
	std::vector<float> *mus_tk_phiErr;
	std::vector<float> *mus_stamu_chi2;
	std::vector<float> *mus_stamu_ndof;
	std::vector<float> *mus_stamu_chg;
	std::vector<float> *mus_stamu_pt;
	std::vector<float> *mus_stamu_px;
	std::vector<float> *mus_stamu_py;
	std::vector<float> *mus_stamu_pz;
	std::vector<float> *mus_stamu_eta;
	std::vector<float> *mus_stamu_phi;
	std::vector<float> *mus_stamu_theta;
	std::vector<float> *mus_stamu_d0dum;
	std::vector<float> *mus_stamu_dz;
	std::vector<float> *mus_stamu_vx;
	std::vector<float> *mus_stamu_vy;
	std::vector<float> *mus_stamu_vz;
	std::vector<float> *mus_stamu_numvalhits;
	std::vector<float> *mus_stamu_numlosthits;
	std::vector<float> *mus_stamu_d0dumErr;
	std::vector<float> *mus_stamu_dzErr;
	std::vector<float> *mus_stamu_ptErr;
	std::vector<float> *mus_stamu_etaErr;
	std::vector<float> *mus_stamu_phiErr;
	std::vector<float> *mus_num_matches;
	std::vector<float> *mus_id_All;
	std::vector<float> *mus_id_AllGlobalMuons;
	std::vector<float> *mus_id_AllStandAloneMuons;
	std::vector<float> *mus_id_AllTrackerMuons;
	std::vector<float> *mus_id_TrackerMuonArbitrated;
	std::vector<float> *mus_id_AllArbitrated;
	std::vector<float> *mus_id_GlobalMuonPromptTight;
	std::vector<float> *mus_id_TMLastStationLoose;
	std::vector<float> *mus_id_TMLastStationTight;
	std::vector<float> *mus_id_TM2DCompatibilityLoose;
	std::vector<float> *mus_id_TM2DCompatibilityTight;
	std::vector<float> *mus_id_TMOneStationLoose;
	std::vector<float> *mus_id_TMOneStationTight;
	std::vector<float> *mus_id_TMLastStationOptimizedLowPtLoose;
	std::vector<float> *mus_id_TMLastStationOptimizedLowPtTight;
	UInt_t Nphotons;
	std::vector<float> *photons_energy;
	std::vector<float> *photons_et;
	std::vector<float> *photons_eta;
	std::vector<float> *photons_phi;
	std::vector<float> *photons_pt;
	std::vector<float> *photons_px;
	std::vector<float> *photons_py;
	std::vector<float> *photons_pz;
	std::vector<float> *photons_status;
	std::vector<float> *photons_theta;
	std::vector<float> *photons_hadOverEM;
	std::vector<float> *photons_scEnergy;
	std::vector<float> *photons_scRawEnergy;
	std::vector<float> *photons_scEta;
	std::vector<float> *photons_scPhi;
	std::vector<float> *photons_scEtaWidth;
	std::vector<float> *photons_scPhiWidth;
	std::vector<float> *photons_tIso;
	std::vector<float> *photons_ecalIso;
	std::vector<float> *photons_hcalIso;
	std::vector<float> *photons_isoEcalRecHitDR04;
	std::vector<float> *photons_isoHcalRecHitDR04;
	std::vector<float> *photons_isoSolidTrkConeDR04;
	std::vector<float> *photons_isoHollowTrkConeDR04;
	std::vector<float> *photons_nTrkSolidConeDR04;
	std::vector<float> *photons_nTrkHollowConeDR04;
	std::vector<float> *photons_isoEcalRecHitDR03;
	std::vector<float> *photons_isoHcalRecHitDR03;
	std::vector<float> *photons_isoSolidTrkConeDR03;
	std::vector<float> *photons_isoHollowTrkConeDR03;
	std::vector<float> *photons_nTrkSolidConeDR03;
	std::vector<float> *photons_nTrkHollowConeDR03;
	std::vector<float> *photons_isAlsoElectron;
	std::vector<float> *photons_hasPixelSeed;
	std::vector<float> *photons_isConverted;
	std::vector<float> *photons_isEBGap;
	std::vector<float> *photons_isEEGap;
	std::vector<float> *photons_isEBEEGap;
	std::vector<float> *photons_isEBPho;
	std::vector<float> *photons_isEEPho;
	std::vector<float> *photons_isLoosePhoton;
	std::vector<float> *photons_isTightPhoton;
	std::vector<float> *photons_r9;
	std::vector<float> *photons_gen_et;
	std::vector<float> *photons_gen_eta;
	std::vector<float> *photons_gen_phi;
	std::vector<float> *photons_gen_id;
	UInt_t Npv;
	std::vector<float> *pv_x;
	std::vector<float> *pv_y;
	std::vector<float> *pv_z;
	std::vector<float> *pv_xErr;
	std::vector<float> *pv_yErr;
	std::vector<float> *pv_zErr;
	UInt_t Ntcmets;
	std::vector<float> *tcmets_et;
	std::vector<float> *tcmets_phi;
	std::vector<float> *tcmets_ex;
	std::vector<float> *tcmets_ey;
	std::vector<float> *tcmets_sumEt;
	std::vector<float> *tcmets_et_muonCor;
	std::vector<float> *tcmets_phi_muonCor;
	UInt_t Ntracks;
	std::vector<float> *tracks_chi2;
	std::vector<float> *tracks_ndof;
	std::vector<float> *tracks_chg;
	std::vector<float> *tracks_pt;
	std::vector<float> *tracks_px;
	std::vector<float> *tracks_py;
	std::vector<float> *tracks_pz;
	std::vector<float> *tracks_eta;
	std::vector<float> *tracks_phi;
	std::vector<float> *tracks_theta;
	std::vector<float> *tracks_d0dum;
	std::vector<float> *tracks_dz;
	std::vector<float> *tracks_vx;
	std::vector<float> *tracks_vy;
	std::vector<float> *tracks_vz;
	std::vector<float> *tracks_numvalhits;
	std::vector<float> *tracks_numlosthits;
	std::vector<float> *tracks_d0dumErr;
	std::vector<float> *tracks_dzErr;
	std::vector<float> *tracks_ptErr;
	std::vector<float> *tracks_etaErr;
	std::vector<float> *tracks_phiErr;
	std::vector<float> *tracks_Nrechits;
	std::vector<float> *tracks_innerHitX;
	std::vector<float> *tracks_innerHitY;
	std::vector<float> *tracks_innerHitZ;
	std::vector<float> *tracks_outerHitX;
	std::vector<float> *tracks_outerHitY;
	std::vector<float> *tracks_outerHitZ;
	std::vector<float> *tracks_highPurity;
	std::vector<float> *tracks_innerLayerMissingHits;
	UInt_t run_number;
	UInt_t event_number;
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
//	TBranch *b_PFJets_btag_softEle;//!
//	TBranch *b_PFJets_btag_softMuon;//!
//	TBranch *b_PFJets_btag_softMuonNoIP;//!
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
//	TBranch *b_jets_btag_softEle; //!
//	TBranch *b_jets_btag_softMuon; //!
//	TBranch *b_jets_btag_softMuonNoIP; //!
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
//	TBranch *b_jetsKT4_btag_softEle; //!
//	TBranch *b_jetsKT4_btag_softMuon; //!
//	TBranch *b_jetsKT4_btag_softMuonNoIP; //!
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
//	TBranch *b_jetsSC5_btag_softEle; //!
//	TBranch *b_jetsSC5_btag_softMuon; //!
//	TBranch *b_jetsSC5_btag_softMuonNoIP; //!
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
private:
	bool datafile, m_useMisslayers, checkTrig, m_studyPDFunc, debug_flag;
	std::string m_jetAlgo, m_metAlgo, HLTBit;

public:
	NTupleReader();
	~NTupleReader() {
	}
	;

	bool GetDebug() const;
	std::string GetHLTBit() const;
	std::string GetJetAlgorithm() const;
	std::string GetMETAlgorithm() const;
	bool GetMissLayersFlag() const;
	bool GetStudyPDFunc() const;
	bool GetTrigger() const;

	bool IsData() const;

	void SetData(bool f);
	void SetDebug(bool);
	void SetHLTBit(std::string);
	void SetJetAlgo(std::string);
	void SetMETAlgo(std::string);
	void SetStudyPDFunc(bool);
	void SetTrigger(bool val);
	void SetTrigger(bool, const std::string);

	void UseMissLayers(bool val);

	void Init(); //initialize tree branches
	void SelectBranches();

private:
	void EnablePFJets();
	void LoadPFJets();


};
#endif /* NTUPLEREADER_H_ */

/*
 * NTupleReader.cc
 *
 *  Created on: Mar 22, 2010
 *      Author: lkreczko
 */

