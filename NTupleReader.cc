/*
 * NTupleReader.cc
 *
 *  Created on: Mar 23, 2010
 *      Author: lkreczko
 */

#include "NTupleReader.h"


NTupleReader::NTupleReader(){

}

bool NTupleReader::GetDebug() const{
		return m_debug;
	}

void NTupleReader::SetDebug(bool val) {
		m_debug = val;
	}


bool NTupleReader::IsData() const {
		return datafile;
	}

bool NTupleReader::GetTrigger() const {
	return checkTrig;
}

void NTupleReader::SetTrigger(bool val) {
	checkTrig = val;
}

void NTupleReader::SetTrigger(bool val, const string hlt) {
	checkTrig = val;
	HLTBit = hlt;
}

bool NTupleReader::GetMissLayersFlag(){
	return m_useMisslayers;
}

void NTupleReader::SetMissLayersFlag(bool val) {
 		m_useMisslayers = val;
}

bool NTupleReader::GetStudyPDFunc() {
		return m_studyPDFunc;
}
void NTupleReader::SetStudyPDFunc(bool f) {
		m_studyPDFunc = f;
}

string NTupleReader::GetHLTBit(){
	return HLTBit;
}

void NTupleReader::SetHLTBit(string bit){
	HLTBit = bit;
}

void NTupleReader::Init() {

	cout << "\n Initializing tree branches\n" << endl;

	// Set object pointer (v3)
	PFJets_energy = 0;
	PFJets_et = 0;
	PFJets_eta = 0;
	PFJets_phi = 0;
	PFJets_pt = 0;
	PFJets_px = 0;
	PFJets_py = 0;
	PFJets_pz = 0;
	PFJets_status = 0;
	PFJets_theta = 0;
	PFJets_chgEmE = 0;
	PFJets_chgHadE = 0;
	PFJets_chgMuE = 0;
	PFJets_chg_Mult = 0;
	PFJets_neutralEmE = 0;
	PFJets_neutralHadE = 0;
	PFJets_neutral_Mult = 0;
	PFJets_mu_Mult = 0;
	PFJets_mass = 0;
	PFMets_et = 0;
	PFMets_phi = 0;
	PFMets_ex = 0;
	PFMets_ey = 0;
	PFMets_sumEt = 0;
	beamSpot_x = 0;
	beamSpot_y = 0;
	beamSpot_z = 0;
	beamSpot_x0Error = 0;
	beamSpot_y0Error = 0;
	beamSpot_z0Error = 0;
	beamSpot_sigmaZ = 0;
	beamSpot_sigmaZ0Error = 0;
	beamSpot_dxdz = 0;
	beamSpot_dxdzError = 0;
	beamSpot_dydz = 0;
	beamSpot_dydzError = 0;
	beamSpot_beamWidthX = 0;
	beamSpot_beamWidthY = 0;
	beamSpot_beamWidthXError = 0;
	beamSpot_beamWidthYError = 0;
	els_energy = 0;
	els_et = 0;
	els_eta = 0;
	els_phi = 0;
	els_pt = 0;
	els_px = 0;
	els_py = 0;
	els_pz = 0;
	els_status = 0;
	els_theta = 0;
	els_closestCtfTrackRef = 0;
	els_isEcalDriven = 0;
	els_isTrackerDriven = 0;
	els_dr03EcalRecHitSumEt = 0;
	els_dr04EcalRecHitSumEt = 0;
	els_dr03HcalTowerSumEt = 0;
	els_dr04HcalTowerSumEt = 0;
	els_gen_id = 0;
	els_gen_phi = 0;
	els_gen_pt = 0;
	els_gen_pz = 0;
	els_gen_px = 0;
	els_gen_py = 0;
	els_gen_eta = 0;
	els_gen_theta = 0;
	els_gen_et = 0;
	els_gen_mother_id = 0;
	els_gen_mother_phi = 0;
	els_gen_mother_pt = 0;
	els_gen_mother_pz = 0;
	els_gen_mother_px = 0;
	els_gen_mother_py = 0;
	els_gen_mother_eta = 0;
	els_gen_mother_theta = 0;
	els_gen_mother_et = 0;
	els_tightId = 0;
	els_looseId = 0;
	els_robustTightId = 0;
	els_robustLooseId = 0;
	els_robustHighEnergyId = 0;
	els_cIso = 0;
	els_tIso = 0;
	els_ecalIso = 0;
	els_hcalIso = 0;
	els_chi2 = 0;
	els_charge = 0;
	els_caloEnergy = 0;
	els_hadOverEm = 0;
	els_eOverPIn = 0;
	els_eSeedOverPOut = 0;
	els_eSCraw = 0;
	els_eSeed = 0;
	els_sigmaEtaEta = 0;
	els_sigmaIEtaIEta = 0;
	els_scE1x5 = 0;
	els_scE2x5Max = 0;
	els_scE5x5 = 0;
	els_dEtaIn = 0;
	els_dPhiIn = 0;
	els_dEtaOut = 0;
	els_dPhiOut = 0;
	els_numvalhits = 0;
	els_numlosthits = 0;
	els_basicClustersSize = 0;
	els_tk_pt = 0;
	els_tk_phi = 0;
	els_tk_eta = 0;
	els_tk_charge = 0;
	els_tk_theta = 0;
	els_shFracInnerHits = 0;
	els_d0dum = 0;
	els_dz = 0;
	els_vx = 0;
	els_vy = 0;
	els_vz = 0;
	els_ndof = 0;
	els_ptError = 0;
	els_d0dumError = 0;
	els_dzError = 0;
	els_etaError = 0;
	els_phiError = 0;
	els_cpx = 0;
	els_cpy = 0;
	els_cpz = 0;
	els_vpx = 0;
	els_vpy = 0;
	els_vpz = 0;
	els_cx = 0;
	els_cy = 0;
	els_cz = 0;
	els_isEE = 0;
	els_isEEGap = 0;
	els_isEB = 0;
	els_isEBGap = 0;
	els_isConvertedPhoton = 0;
	els_innerLayerMissingHits = 0;
	jets_energy = 0;
	jets_et = 0;
	jets_eta = 0;
	jets_phi = 0;
	jets_pt = 0;
	jets_px = 0;
	jets_py = 0;
	jets_pz = 0;
	jets_status = 0;
	jets_theta = 0;
	jets_parton_Id = 0;
	jets_parton_motherId = 0;
	jets_parton_pt = 0;
	jets_parton_phi = 0;
	jets_parton_eta = 0;
	jets_parton_Energy = 0;
	jets_parton_mass = 0;
	jets_parton_motherID = 0;
	jets_gen_et = 0;
	jets_gen_pt = 0;
	jets_gen_eta = 0;
	jets_gen_phi = 0;
	jets_gen_mass = 0;
	jets_gen_Energy = 0;
	jets_gen_Id = 0;
	jets_gen_motherID = 0;
	jets_gen_threeCharge = 0;
	jets_partonFlavour = 0;
	jets_btag_TC_highPur = 0;
	jets_btag_TC_highEff = 0;
	jets_btag_jetProb = 0;
	jets_btag_jetBProb = 0;
	jets_btag_softEle = 0;
	jets_btag_softMuon = 0;
	jets_btag_softMuonNoIP = 0;
	jets_btag_secVertex = 0;
	jets_chgEmE = 0;
	jets_chgHadE = 0;
	jets_chgMuE = 0;
	jets_chg_Mult = 0;
	jets_neutralEmE = 0;
	jets_neutralHadE = 0;
	jets_neutral_Mult = 0;
	jets_mu_Mult = 0;
	jets_emf = 0;
	jets_ehf = 0;
	jets_n60 = 0;
	jets_n90 = 0;
	jets_area = 0;
	jets_mass = 0;

	jetsKT4_energy = 0;
	jetsKT4_et = 0;
	jetsKT4_eta = 0;
	jetsKT4_phi = 0;
	jetsKT4_pt = 0;
	jetsKT4_px = 0;
	jetsKT4_py = 0;
	jetsKT4_pz = 0;
	jetsKT4_status = 0;
	jetsKT4_theta = 0;
	jetsKT4_btag_TC_highPur = 0;
	jetsKT4_btag_TC_highEff = 0;
	jetsKT4_btag_jetProb = 0;
	jetsKT4_btag_jetBProb = 0;
	jetsKT4_btag_softEle = 0;
	jetsKT4_btag_softMuon = 0;
	jetsKT4_btag_softMuonNoIP = 0;
	jetsKT4_btag_secVertex = 0;
	jetsKT4_chgEmE = 0;
	jetsKT4_chgHadE = 0;
	jetsKT4_chgMuE = 0;
	jetsKT4_chg_Mult = 0;
	jetsKT4_neutralEmE = 0;
	jetsKT4_neutralHadE = 0;
	jetsKT4_neutral_Mult = 0;
	jetsKT4_mu_Mult = 0;
	jetsKT4_emf = 0;
	jetsKT4_ehf = 0;
	jetsKT4_n60 = 0;
	jetsKT4_n90 = 0;
	jetsKT4_area = 0;
	jetsKT4_mass = 0;

	jetsSC5_energy = 0;
	jetsSC5_et = 0;
	jetsSC5_eta = 0;
	jetsSC5_phi = 0;
	jetsSC5_pt = 0;
	jetsSC5_px = 0;
	jetsSC5_py = 0;
	jetsSC5_pz = 0;
	jetsSC5_status = 0;
	jetsSC5_theta = 0;
	jetsSC5_btag_TC_highPur = 0;
	jetsSC5_btag_TC_highEff = 0;
	jetsSC5_btag_jetProb = 0;
	jetsSC5_btag_jetBProb = 0;
	jetsSC5_btag_softEle = 0;
	jetsSC5_btag_softMuon = 0;
	jetsSC5_btag_softMuonNoIP = 0;
	jetsSC5_btag_secVertex = 0;
	jetsSC5_chgEmE = 0;
	jetsSC5_chgHadE = 0;
	jetsSC5_chgMuE = 0;
	jetsSC5_chg_Mult = 0;
	jetsSC5_neutralEmE = 0;
	jetsSC5_neutralHadE = 0;
	jetsSC5_neutral_Mult = 0;
	jetsSC5_mu_Mult = 0;
	jetsSC5_emf = 0;
	jetsSC5_ehf = 0;
	jetsSC5_n60 = 0;
	jetsSC5_n90 = 0;
	jetsSC5_area = 0;
	jetsSC5_mass = 0;
	mc_doc_id = 0;
	mc_doc_pt = 0;
	mc_doc_px = 0;
	mc_doc_py = 0;
	mc_doc_pz = 0;
	mc_doc_eta = 0;
	mc_doc_phi = 0;
	mc_doc_theta = 0;
	mc_doc_energy = 0;
	mc_doc_status = 0;
	mc_doc_charge = 0;
	mc_doc_mother_id = 0;
	mc_doc_grandmother_id = 0;
	mc_doc_ggrandmother_id = 0;
	mc_doc_mother_pt = 0;
	mc_doc_vertex_x = 0;
	mc_doc_vertex_y = 0;
	mc_doc_vertex_z = 0;
	mc_doc_mass = 0;
	mc_doc_numOfDaughters = 0;
	mc_doc_numOfMothers = 0;
	mets_et = 0;
	mets_phi = 0;
	mets_ex = 0;
	mets_ey = 0;
	mets_gen_et = 0;
	mets_gen_phi = 0;
	mets_sign = 0;
	mets_sumEt = 0;
	mets_unCPhi = 0;
	mets_unCPt = 0;
	mets_et_muonCor = 0;
	mets_phi_muonCor = 0;
	mets_et_JESCor = 0;
	mets_phi_JESCor = 0;
	metsKT4_et = 0;
	metsKT4_phi = 0;
	metsKT4_ex = 0;
	metsKT4_ey = 0;
	metsKT4_sumEt = 0;
	metsKT4_et_JESCor = 0;
	metsKT4_phi_JESCor = 0;
	metsSC5_et = 0;
	metsSC5_phi = 0;
	metsSC5_ex = 0;
	metsSC5_ey = 0;
	metsSC5_sumEt = 0;
	metsSC5_et_JESCor = 0;
	metsSC5_phi_JESCor = 0;
	mus_energy = 0;
	mus_et = 0;
	mus_eta = 0;
	mus_phi = 0;
	mus_pt = 0;
	mus_px = 0;
	mus_py = 0;
	mus_pz = 0;
	mus_status = 0;
	mus_theta = 0;
	mus_gen_id = 0;
	mus_gen_phi = 0;
	mus_gen_pt = 0;
	mus_gen_pz = 0;
	mus_gen_px = 0;
	mus_gen_py = 0;
	mus_gen_eta = 0;
	mus_gen_theta = 0;
	mus_gen_et = 0;
	mus_gen_mother_id = 0;
	mus_gen_mother_phi = 0;
	mus_gen_mother_pt = 0;
	mus_gen_mother_pz = 0;
	mus_gen_mother_px = 0;
	mus_gen_mother_py = 0;
	mus_gen_mother_eta = 0;
	mus_gen_mother_theta = 0;
	mus_gen_mother_et = 0;
	mus_tkHits = 0;
	mus_cIso = 0;
	mus_tIso = 0;
	mus_ecalIso = 0;
	mus_hcalIso = 0;
	mus_ecalvetoDep = 0;
	mus_hcalvetoDep = 0;
	mus_calEnergyEm = 0;
	mus_calEnergyHad = 0;
	mus_calEnergyHo = 0;
	mus_calEnergyEmS9 = 0;
	mus_calEnergyHadS9 = 0;
	mus_calEnergyHoS9 = 0;
	mus_iso03_sumPt = 0;
	mus_iso03_emEt = 0;
	mus_iso03_hadEt = 0;
	mus_iso03_hoEt = 0;
	mus_iso03_nTracks = 0;
	mus_iso05_sumPt = 0;
	mus_iso05_emEt = 0;
	mus_iso05_hadEt = 0;
	mus_iso05_hoEt = 0;
	mus_iso05_nTracks = 0;
	mus_charge = 0;
	mus_cm_chi2 = 0;
	mus_cm_ndof = 0;
	mus_cm_chg = 0;
	mus_cm_pt = 0;
	mus_cm_px = 0;
	mus_cm_py = 0;
	mus_cm_pz = 0;
	mus_cm_eta = 0;
	mus_cm_phi = 0;
	mus_cm_theta = 0;
	mus_cm_d0dum = 0;
	mus_cm_dz = 0;
	mus_cm_vx = 0;
	mus_cm_vy = 0;
	mus_cm_vz = 0;
	mus_cm_numvalhits = 0;
	mus_cm_numlosthits = 0;
	mus_cm_d0dumErr = 0;
	mus_cm_dzErr = 0;
	mus_cm_ptErr = 0;
	mus_cm_etaErr = 0;
	mus_cm_phiErr = 0;
	mus_tk_chi2 = 0;
	mus_tk_ndof = 0;
	mus_tk_chg = 0;
	mus_tk_pt = 0;
	mus_tk_px = 0;
	mus_tk_py = 0;
	mus_tk_pz = 0;
	mus_tk_eta = 0;
	mus_tk_phi = 0;
	mus_tk_theta = 0;
	mus_tk_d0dum = 0;
	mus_tk_dz = 0;
	mus_tk_vx = 0;
	mus_tk_vy = 0;
	mus_tk_vz = 0;
	mus_tk_numvalhits = 0;
	mus_tk_numlosthits = 0;
	mus_tk_d0dumErr = 0;
	mus_tk_dzErr = 0;
	mus_tk_ptErr = 0;
	mus_tk_etaErr = 0;
	mus_tk_phiErr = 0;
	mus_stamu_chi2 = 0;
	mus_stamu_ndof = 0;
	mus_stamu_chg = 0;
	mus_stamu_pt = 0;
	mus_stamu_px = 0;
	mus_stamu_py = 0;
	mus_stamu_pz = 0;
	mus_stamu_eta = 0;
	mus_stamu_phi = 0;
	mus_stamu_theta = 0;
	mus_stamu_d0dum = 0;
	mus_stamu_dz = 0;
	mus_stamu_vx = 0;
	mus_stamu_vy = 0;
	mus_stamu_vz = 0;
	mus_stamu_numvalhits = 0;
	mus_stamu_numlosthits = 0;
	mus_stamu_d0dumErr = 0;
	mus_stamu_dzErr = 0;
	mus_stamu_ptErr = 0;
	mus_stamu_etaErr = 0;
	mus_stamu_phiErr = 0;
	mus_num_matches = 0;
	mus_id_All = 0;
	mus_id_AllGlobalMuons = 0;
	mus_id_AllStandAloneMuons = 0;
	mus_id_AllTrackerMuons = 0;
	mus_id_TrackerMuonArbitrated = 0;
	mus_id_AllArbitrated = 0;
	mus_id_GlobalMuonPromptTight = 0;
	mus_id_TMLastStationLoose = 0;
	mus_id_TMLastStationTight = 0;
	mus_id_TM2DCompatibilityLoose = 0;
	mus_id_TM2DCompatibilityTight = 0;
	mus_id_TMOneStationLoose = 0;
	mus_id_TMOneStationTight = 0;
	mus_id_TMLastStationOptimizedLowPtLoose = 0;
	mus_id_TMLastStationOptimizedLowPtTight = 0;
	photons_energy = 0;
	photons_et = 0;
	photons_eta = 0;
	photons_phi = 0;
	photons_pt = 0;
	photons_px = 0;
	photons_py = 0;
	photons_pz = 0;
	photons_status = 0;
	photons_theta = 0;
	photons_hadOverEM = 0;
	photons_scEnergy = 0;
	photons_scRawEnergy = 0;
	photons_scEta = 0;
	photons_scPhi = 0;
	photons_scEtaWidth = 0;
	photons_scPhiWidth = 0;
	photons_tIso = 0;
	photons_ecalIso = 0;
	photons_hcalIso = 0;
	photons_isoEcalRecHitDR04 = 0;
	photons_isoHcalRecHitDR04 = 0;
	photons_isoSolidTrkConeDR04 = 0;
	photons_isoHollowTrkConeDR04 = 0;
	photons_nTrkSolidConeDR04 = 0;
	photons_nTrkHollowConeDR04 = 0;
	photons_isoEcalRecHitDR03 = 0;
	photons_isoHcalRecHitDR03 = 0;
	photons_isoSolidTrkConeDR03 = 0;
	photons_isoHollowTrkConeDR03 = 0;
	photons_nTrkSolidConeDR03 = 0;
	photons_nTrkHollowConeDR03 = 0;
	photons_isAlsoElectron = 0;
	photons_hasPixelSeed = 0;
	photons_isConverted = 0;
	photons_isEBGap = 0;
	photons_isEEGap = 0;
	photons_isEBEEGap = 0;
	photons_isEBPho = 0;
	photons_isEEPho = 0;
	photons_isLoosePhoton = 0;
	photons_isTightPhoton = 0;
	photons_r9 = 0;
	photons_gen_et = 0;
	photons_gen_eta = 0;
	photons_gen_phi = 0;
	photons_gen_id = 0;
	pv_x = 0;
	pv_y = 0;
	pv_z = 0;
	pv_xErr = 0;
	pv_yErr = 0;
	pv_zErr = 0;
	tcmets_et = 0;
	tcmets_phi = 0;
	tcmets_ex = 0;
	tcmets_ey = 0;
	tcmets_sumEt = 0;
	tcmets_et_muonCor = 0;
	tcmets_phi_muonCor = 0;
	tracks_chi2 = 0;
	tracks_ndof = 0;
	tracks_chg = 0;
	tracks_pt = 0;
	tracks_px = 0;
	tracks_py = 0;
	tracks_pz = 0;
	tracks_eta = 0;
	tracks_phi = 0;
	tracks_theta = 0;
	tracks_d0dum = 0;
	tracks_dz = 0;
	tracks_vx = 0;
	tracks_vy = 0;
	tracks_vz = 0;
	tracks_numvalhits = 0;
	tracks_numlosthits = 0;
	tracks_d0dumErr = 0;
	tracks_dzErr = 0;
	tracks_ptErr = 0;
	tracks_etaErr = 0;
	tracks_phiErr = 0;
	tracks_Nrechits = 0;
	tracks_innerHitX = 0;
	tracks_innerHitY = 0;
	tracks_innerHitZ = 0;
	tracks_outerHitX = 0;
	tracks_outerHitY = 0;
	tracks_outerHitZ = 0;
	tracks_highPurity = 0;
	tracks_innerLayerMissingHits = 0;

	// Set branch addresses and branch pointers
	if (m_jetAlgo == "pfjet") {
		chain->SetBranchAddress("NPFJets", &NPFJets, &b_NPFJets);
		chain->SetBranchAddress("PFJets_energy", &PFJets_energy, &b_PFJets_energy);
		chain->SetBranchAddress("PFJets_et", &PFJets_et, &b_PFJets_et);
		chain->SetBranchAddress("PFJets_eta", &PFJets_eta, &b_PFJets_eta);
		chain->SetBranchAddress("PFJets_phi", &PFJets_phi, &b_PFJets_phi);
		chain->SetBranchAddress("PFJets_pt", &PFJets_pt, &b_PFJets_pt);
		chain->SetBranchAddress("PFJets_px", &PFJets_px, &b_PFJets_px);
		chain->SetBranchAddress("PFJets_py", &PFJets_py, &b_PFJets_py);
		chain->SetBranchAddress("PFJets_pz", &PFJets_pz, &b_PFJets_pz);
		chain->SetBranchAddress("PFJets_status", &PFJets_status, &b_PFJets_status);
		chain->SetBranchAddress("PFJets_theta", &PFJets_theta, &b_PFJets_theta);
		chain->SetBranchAddress("PFJets_chgEmE", &PFJets_chgEmE, &b_PFJets_chgEmE);
		chain->SetBranchAddress("PFJets_chgHadE", &PFJets_chgHadE, &b_PFJets_chgHadE);
		chain->SetBranchAddress("PFJets_chgMuE", &PFJets_chgMuE, &b_PFJets_chgMuE);
		chain->SetBranchAddress("PFJets_chg_Mult", &PFJets_chg_Mult, &b_PFJets_chg_Mult);
		chain->SetBranchAddress("PFJets_neutralEmE", &PFJets_neutralEmE, &b_PFJets_neutralEmE);
		chain->SetBranchAddress("PFJets_neutralHadE", &PFJets_neutralHadE, &b_PFJets_neutralHadE);
		chain->SetBranchAddress("PFJets_neutral_Mult", &PFJets_neutral_Mult, &b_PFJets_neutral_Mult);
		chain->SetBranchAddress("PFJets_mu_Mult", &PFJets_mu_Mult, &b_PFJets_mu_Mult);
		chain->SetBranchAddress("PFJets_mass", &PFJets_mass, &b_PFJets_mass);
	}
	if (m_metAlgo == "pfmet") {
		chain->SetBranchAddress("NPFMets", &NPFMets, &b_NPFMets);
		chain->SetBranchAddress("PFMets_et", &PFMets_et, &b_PFMets_et);
		chain->SetBranchAddress("PFMets_phi", &PFMets_phi, &b_PFMets_phi);
		chain->SetBranchAddress("PFMets_ex", &PFMets_ex, &b_PFMets_ex);
		chain->SetBranchAddress("PFMets_ey", &PFMets_ey, &b_PFMets_ey);
		chain->SetBranchAddress("PFMets_sumEt", &PFMets_sumEt, &b_PFMets_sumEt);
	}
	chain->SetBranchAddress("NbeamSpot", &NbeamSpot, &b_NbeamSpot);
	chain->SetBranchAddress("beamSpot_x", &beamSpot_x, &b_beamSpot_x);
	chain->SetBranchAddress("beamSpot_y", &beamSpot_y, &b_beamSpot_y);
	chain->SetBranchAddress("beamSpot_z", &beamSpot_z, &b_beamSpot_z);
	chain->SetBranchAddress("beamSpot_x0Error", &beamSpot_x0Error, &b_beamSpot_x0Error);
	chain->SetBranchAddress("beamSpot_y0Error", &beamSpot_y0Error, &b_beamSpot_y0Error);
	chain->SetBranchAddress("beamSpot_z0Error", &beamSpot_z0Error, &b_beamSpot_z0Error);
	chain->SetBranchAddress("beamSpot_sigmaZ", &beamSpot_sigmaZ, &b_beamSpot_sigmaZ);
	chain->SetBranchAddress("beamSpot_sigmaZ0Error", &beamSpot_sigmaZ0Error, &b_beamSpot_sigmaZ0Error);
	chain->SetBranchAddress("beamSpot_dxdz", &beamSpot_dxdz, &b_beamSpot_dxdz);
	chain->SetBranchAddress("beamSpot_dxdzError", &beamSpot_dxdzError, &b_beamSpot_dxdzError);
	chain->SetBranchAddress("beamSpot_dydz", &beamSpot_dydz, &b_beamSpot_dydz);
	chain->SetBranchAddress("beamSpot_dydzError", &beamSpot_dydzError, &b_beamSpot_dydzError);
	chain->SetBranchAddress("beamSpot_beamWidthX", &beamSpot_beamWidthX, &b_beamSpot_beamWidthX);
	chain->SetBranchAddress("beamSpot_beamWidthY", &beamSpot_beamWidthY, &b_beamSpot_beamWidthY);
	chain->SetBranchAddress("beamSpot_beamWidthXError", &beamSpot_beamWidthXError, &b_beamSpot_beamWidthXError);
	chain->SetBranchAddress("beamSpot_beamWidthYError", &beamSpot_beamWidthYError, &b_beamSpot_beamWidthYError);
	chain->SetBranchAddress("Nels", &Nels, &b_Nels);
	chain->SetBranchAddress("els_energy", &els_energy, &b_els_energy);
	chain->SetBranchAddress("els_et", &els_et, &b_els_et);
	chain->SetBranchAddress("els_eta", &els_eta, &b_els_eta);
	chain->SetBranchAddress("els_phi", &els_phi, &b_els_phi);
	chain->SetBranchAddress("els_pt", &els_pt, &b_els_pt);
	chain->SetBranchAddress("els_px", &els_px, &b_els_px);
	chain->SetBranchAddress("els_py", &els_py, &b_els_py);
	chain->SetBranchAddress("els_pz", &els_pz, &b_els_pz);
	chain->SetBranchAddress("els_status", &els_status, &b_els_status);
	chain->SetBranchAddress("els_theta", &els_theta, &b_els_theta);
	chain->SetBranchAddress("els_closestCtfTrackRef", &els_closestCtfTrackRef, &b_els_closestCtfTrackRef);
	chain->SetBranchAddress("els_isEcalDriven", &els_isEcalDriven, &b_els_isEcalDriven);
	chain->SetBranchAddress("els_isTrackerDriven", &els_isTrackerDriven, &b_els_isTrackerDriven);
	chain->SetBranchAddress("els_dr03EcalRecHitSumEt", &els_dr03EcalRecHitSumEt, &b_els_dr03EcalRecHitSumEt);
	chain->SetBranchAddress("els_dr04EcalRecHitSumEt", &els_dr04EcalRecHitSumEt, &b_els_dr04EcalRecHitSumEt);
	chain->SetBranchAddress("els_dr03HcalTowerSumEt", &els_dr03HcalTowerSumEt, &b_els_dr03HcalTowerSumEt);
	chain->SetBranchAddress("els_dr04HcalTowerSumEt", &els_dr04HcalTowerSumEt, &b_els_dr04HcalTowerSumEt);
	chain->SetBranchAddress("els_gen_id", &els_gen_id, &b_els_gen_id);
	chain->SetBranchAddress("els_gen_phi", &els_gen_phi, &b_els_gen_phi);
	chain->SetBranchAddress("els_gen_pt", &els_gen_pt, &b_els_gen_pt);
	chain->SetBranchAddress("els_gen_pz", &els_gen_pz, &b_els_gen_pz);
	chain->SetBranchAddress("els_gen_px", &els_gen_px, &b_els_gen_px);
	chain->SetBranchAddress("els_gen_py", &els_gen_py, &b_els_gen_py);
	chain->SetBranchAddress("els_gen_eta", &els_gen_eta, &b_els_gen_eta);
	chain->SetBranchAddress("els_gen_theta", &els_gen_theta, &b_els_gen_theta);
	chain->SetBranchAddress("els_gen_et", &els_gen_et, &b_els_gen_et);
	chain->SetBranchAddress("els_gen_mother_id", &els_gen_mother_id, &b_els_gen_mother_id);
	chain->SetBranchAddress("els_gen_mother_phi", &els_gen_mother_phi, &b_els_gen_mother_phi);
	chain->SetBranchAddress("els_gen_mother_pt", &els_gen_mother_pt, &b_els_gen_mother_pt);
	chain->SetBranchAddress("els_gen_mother_pz", &els_gen_mother_pz, &b_els_gen_mother_pz);
	chain->SetBranchAddress("els_gen_mother_px", &els_gen_mother_px, &b_els_gen_mother_px);
	chain->SetBranchAddress("els_gen_mother_py", &els_gen_mother_py, &b_els_gen_mother_py);
	chain->SetBranchAddress("els_gen_mother_eta", &els_gen_mother_eta, &b_els_gen_mother_eta);
	chain->SetBranchAddress("els_gen_mother_theta", &els_gen_mother_theta, &b_els_gen_mother_theta);
	chain->SetBranchAddress("els_gen_mother_et", &els_gen_mother_et, &b_els_gen_mother_et);
	chain->SetBranchAddress("els_tightId", &els_tightId, &b_els_tightId);
	chain->SetBranchAddress("els_looseId", &els_looseId, &b_els_looseId);
	chain->SetBranchAddress("els_robustTightId", &els_robustTightId, &b_els_robustTightId);
	chain->SetBranchAddress("els_robustLooseId", &els_robustLooseId, &b_els_robustLooseId);
	chain->SetBranchAddress("els_robustHighEnergyId", &els_robustHighEnergyId, &b_els_robustHighEnergyId);
	chain->SetBranchAddress("els_cIso", &els_cIso, &b_els_cIso);
	chain->SetBranchAddress("els_tIso", &els_tIso, &b_els_tIso);
	chain->SetBranchAddress("els_ecalIso", &els_ecalIso, &b_els_ecalIso);
	chain->SetBranchAddress("els_hcalIso", &els_hcalIso, &b_els_hcalIso);
	chain->SetBranchAddress("els_chi2", &els_chi2, &b_els_chi2);
	chain->SetBranchAddress("els_charge", &els_charge, &b_els_charge);
	chain->SetBranchAddress("els_caloEnergy", &els_caloEnergy, &b_els_caloEnergy);
	chain->SetBranchAddress("els_hadOverEm", &els_hadOverEm, &b_els_hadOverEm);
	chain->SetBranchAddress("els_eOverPIn", &els_eOverPIn, &b_els_eOverPIn);
	chain->SetBranchAddress("els_eSeedOverPOut", &els_eSeedOverPOut, &b_els_eSeedOverPOut);
	chain->SetBranchAddress("els_eSCraw", &els_eSCraw, &b_els_eSCraw);
	chain->SetBranchAddress("els_eSeed", &els_eSeed, &b_els_eSeed);
	chain->SetBranchAddress("els_sigmaEtaEta", &els_sigmaEtaEta, &b_els_sigmaEtaEta);
	chain->SetBranchAddress("els_sigmaIEtaIEta", &els_sigmaIEtaIEta, &b_els_sigmaIEtaIEta);
	chain->SetBranchAddress("els_scE1x5", &els_scE1x5, &b_els_scE1x5);
	chain->SetBranchAddress("els_scE2x5Max", &els_scE2x5Max, &b_els_scE2x5Max);
	chain->SetBranchAddress("els_scE5x5", &els_scE5x5, &b_els_scE5x5);
	chain->SetBranchAddress("els_dEtaIn", &els_dEtaIn, &b_els_dEtaIn);
	chain->SetBranchAddress("els_dPhiIn", &els_dPhiIn, &b_els_dPhiIn);
	chain->SetBranchAddress("els_dEtaOut", &els_dEtaOut, &b_els_dEtaOut);
	chain->SetBranchAddress("els_dPhiOut", &els_dPhiOut, &b_els_dPhiOut);
	chain->SetBranchAddress("els_numvalhits", &els_numvalhits, &b_els_numvalhits);
	chain->SetBranchAddress("els_numlosthits", &els_numlosthits, &b_els_numlosthits);
	chain->SetBranchAddress("els_basicClustersSize", &els_basicClustersSize, &b_els_basicClustersSize);
	chain->SetBranchAddress("els_tk_pt", &els_tk_pt, &b_els_tk_pt);
	chain->SetBranchAddress("els_tk_phi", &els_tk_phi, &b_els_tk_phi);
	chain->SetBranchAddress("els_tk_eta", &els_tk_eta, &b_els_tk_eta);
	chain->SetBranchAddress("els_tk_charge", &els_tk_charge, &b_els_tk_charge);
	chain->SetBranchAddress("els_tk_theta", &els_tk_theta, &b_els_tk_theta);
	chain->SetBranchAddress("els_shFracInnerHits", &els_shFracInnerHits, &b_els_shFracInnerHits);
	chain->SetBranchAddress("els_d0dum", &els_d0dum, &b_els_d0dum);
	chain->SetBranchAddress("els_dz", &els_dz, &b_els_dz);
	chain->SetBranchAddress("els_vx", &els_vx, &b_els_vx);
	chain->SetBranchAddress("els_vy", &els_vy, &b_els_vy);
	chain->SetBranchAddress("els_vz", &els_vz, &b_els_vz);
	chain->SetBranchAddress("els_ndof", &els_ndof, &b_els_ndof);
	chain->SetBranchAddress("els_ptError", &els_ptError, &b_els_ptError);
	chain->SetBranchAddress("els_d0dumError", &els_d0dumError, &b_els_d0dumError);
	chain->SetBranchAddress("els_dzError", &els_dzError, &b_els_dzError);
	chain->SetBranchAddress("els_etaError", &els_etaError, &b_els_etaError);
	chain->SetBranchAddress("els_phiError", &els_phiError, &b_els_phiError);
	chain->SetBranchAddress("els_cpx", &els_cpx, &b_els_cpx);
	chain->SetBranchAddress("els_cpy", &els_cpy, &b_els_cpy);
	chain->SetBranchAddress("els_cpz", &els_cpz, &b_els_cpz);
	chain->SetBranchAddress("els_vpx", &els_vpx, &b_els_vpx);
	chain->SetBranchAddress("els_vpy", &els_vpy, &b_els_vpy);
	chain->SetBranchAddress("els_vpz", &els_vpz, &b_els_vpz);
	chain->SetBranchAddress("els_cx", &els_cx, &b_els_cx);
	chain->SetBranchAddress("els_cy", &els_cy, &b_els_cy);
	chain->SetBranchAddress("els_cz", &els_cz, &b_els_cz);
	chain->SetBranchAddress("els_isEE", &els_isEE, &b_els_isEE);
	chain->SetBranchAddress("els_isEEGap", &els_isEEGap, &b_els_isEEGap);
	chain->SetBranchAddress("els_isEB", &els_isEB, &b_els_isEB);
	chain->SetBranchAddress("els_isEBGap", &els_isEBGap, &b_els_isEBGap);
	chain->SetBranchAddress("els_isConvertedPhoton", &els_isConvertedPhoton, &b_els_isConvertedPhoton);
	if (m_useMisslayers)
		chain->SetBranchAddress("els_innerLayerMissingHits", &els_innerLayerMissingHits, &b_els_innerLayerMissingHits);

	chain->SetBranchAddress("Njets", &Njets, &b_Njets);
	chain->SetBranchAddress("jets_energy", &jets_energy, &b_jets_energy);
	chain->SetBranchAddress("jets_et", &jets_et, &b_jets_et);
	chain->SetBranchAddress("jets_eta", &jets_eta, &b_jets_eta);
	chain->SetBranchAddress("jets_phi", &jets_phi, &b_jets_phi);
	chain->SetBranchAddress("jets_pt", &jets_pt, &b_jets_pt);
	chain->SetBranchAddress("jets_px", &jets_px, &b_jets_px);
	chain->SetBranchAddress("jets_py", &jets_py, &b_jets_py);
	chain->SetBranchAddress("jets_pz", &jets_pz, &b_jets_pz);
	chain->SetBranchAddress("jets_status", &jets_status, &b_jets_status);
	chain->SetBranchAddress("jets_theta", &jets_theta, &b_jets_theta);
	chain->SetBranchAddress("jets_parton_Id", &jets_parton_Id, &b_jets_parton_Id);
	chain->SetBranchAddress("jets_parton_motherId", &jets_parton_motherId, &b_jets_parton_motherId);
	chain->SetBranchAddress("jets_parton_pt", &jets_parton_pt, &b_jets_parton_pt);
	chain->SetBranchAddress("jets_parton_phi", &jets_parton_phi, &b_jets_parton_phi);
	chain->SetBranchAddress("jets_parton_eta", &jets_parton_eta, &b_jets_parton_eta);
	chain->SetBranchAddress("jets_parton_Energy", &jets_parton_Energy, &b_jets_parton_Energy);
	chain->SetBranchAddress("jets_parton_mass", &jets_parton_mass, &b_jets_parton_mass);
	chain->SetBranchAddress("jets_parton_motherID", &jets_parton_motherID, &b_jets_parton_motherID);
	chain->SetBranchAddress("jets_gen_et", &jets_gen_et, &b_jets_gen_et);
	chain->SetBranchAddress("jets_gen_pt", &jets_gen_pt, &b_jets_gen_pt);
	chain->SetBranchAddress("jets_gen_eta", &jets_gen_eta, &b_jets_gen_eta);
	chain->SetBranchAddress("jets_gen_phi", &jets_gen_phi, &b_jets_gen_phi);
	chain->SetBranchAddress("jets_gen_mass", &jets_gen_mass, &b_jets_gen_mass);
	chain->SetBranchAddress("jets_gen_Energy", &jets_gen_Energy, &b_jets_gen_Energy);
	chain->SetBranchAddress("jets_gen_Id", &jets_gen_Id, &b_jets_gen_Id);
	chain->SetBranchAddress("jets_gen_motherID", &jets_gen_motherID, &b_jets_gen_motherID);
	chain->SetBranchAddress("jets_gen_threeCharge", &jets_gen_threeCharge, &b_jets_gen_threeCharge);
	chain->SetBranchAddress("jets_partonFlavour", &jets_partonFlavour, &b_jets_partonFlavour);
	chain->SetBranchAddress("jets_btag_TC_highPur", &jets_btag_TC_highPur, &b_jets_btag_TC_highPur);
	chain->SetBranchAddress("jets_btag_TC_highEff", &jets_btag_TC_highEff, &b_jets_btag_TC_highEff);
	chain->SetBranchAddress("jets_btag_jetProb", &jets_btag_jetProb, &b_jets_btag_jetProb);
	chain->SetBranchAddress("jets_btag_jetBProb", &jets_btag_jetBProb, &b_jets_btag_jetBProb);
	chain->SetBranchAddress("jets_btag_softEle", &jets_btag_softEle, &b_jets_btag_softEle);
	chain->SetBranchAddress("jets_btag_softMuon", &jets_btag_softMuon, &b_jets_btag_softMuon);
	chain->SetBranchAddress("jets_btag_softMuonNoIP", &jets_btag_softMuonNoIP, &b_jets_btag_softMuonNoIP);
	chain->SetBranchAddress("jets_btag_secVertex", &jets_btag_secVertex, &b_jets_btag_secVertex);
	chain->SetBranchAddress("jets_chgEmE", &jets_chgEmE, &b_jets_chgEmE);
	chain->SetBranchAddress("jets_chgHadE", &jets_chgHadE, &b_jets_chgHadE);
	chain->SetBranchAddress("jets_chgMuE", &jets_chgMuE, &b_jets_chgMuE);
	chain->SetBranchAddress("jets_chg_Mult", &jets_chg_Mult, &b_jets_chg_Mult);
	chain->SetBranchAddress("jets_neutralEmE", &jets_neutralEmE, &b_jets_neutralEmE);
	chain->SetBranchAddress("jets_neutralHadE", &jets_neutralHadE, &b_jets_neutralHadE);
	chain->SetBranchAddress("jets_neutral_Mult", &jets_neutral_Mult, &b_jets_neutral_Mult);
	chain->SetBranchAddress("jets_mu_Mult", &jets_mu_Mult, &b_jets_mu_Mult);
	chain->SetBranchAddress("jets_emf", &jets_emf, &b_jets_emf);
	chain->SetBranchAddress("jets_ehf", &jets_ehf, &b_jets_ehf);
	chain->SetBranchAddress("jets_n60", &jets_n60, &b_jets_n60);
	chain->SetBranchAddress("jets_n90", &jets_n90, &b_jets_n90);
	chain->SetBranchAddress("jets_area", &jets_area, &b_jets_area);
	chain->SetBranchAddress("jets_mass", &jets_mass, &b_jets_mass);

	if (m_jetAlgo == "KT4") {
		chain->SetBranchAddress("jetsKT4_energy", &jetsKT4_energy, &b_jetsKT4_energy);
		chain->SetBranchAddress("jetsKT4_et", &jetsKT4_et, &b_jetsKT4_et);
		chain->SetBranchAddress("jetsKT4_eta", &jetsKT4_eta, &b_jetsKT4_eta);
		chain->SetBranchAddress("jetsKT4_phi", &jetsKT4_phi, &b_jetsKT4_phi);
		chain->SetBranchAddress("jetsKT4_pt", &jetsKT4_pt, &b_jetsKT4_pt);
		chain->SetBranchAddress("jetsKT4_px", &jetsKT4_px, &b_jetsKT4_px);
		chain->SetBranchAddress("jetsKT4_py", &jetsKT4_py, &b_jetsKT4_py);
		chain->SetBranchAddress("jetsKT4_pz", &jetsKT4_pz, &b_jetsKT4_pz);
		chain->SetBranchAddress("jetsKT4_status", &jetsKT4_status, &b_jetsKT4_status);
		chain->SetBranchAddress("jetsKT4_theta", &jetsKT4_theta, &b_jetsKT4_theta);
		chain->SetBranchAddress("jetsKT4_btag_TC_highPur", &jetsKT4_btag_TC_highPur, &b_jetsKT4_btag_TC_highPur);
		chain->SetBranchAddress("jetsKT4_btag_TC_highEff", &jetsKT4_btag_TC_highEff, &b_jetsKT4_btag_TC_highEff);
		chain->SetBranchAddress("jetsKT4_btag_jetProb", &jetsKT4_btag_jetProb, &b_jetsKT4_btag_jetProb);
		chain->SetBranchAddress("jetsKT4_btag_jetBProb", &jetsKT4_btag_jetBProb, &b_jetsKT4_btag_jetBProb);
		chain->SetBranchAddress("jetsKT4_btag_softEle", &jetsKT4_btag_softEle, &b_jetsKT4_btag_softEle);
		chain->SetBranchAddress("jetsKT4_btag_softMuon", &jetsKT4_btag_softMuon, &b_jetsKT4_btag_softMuon);
		chain->SetBranchAddress("jetsKT4_btag_softMuonNoIP", &jetsKT4_btag_softMuonNoIP, &b_jetsKT4_btag_softMuonNoIP);
		chain->SetBranchAddress("jetsKT4_btag_secVertex", &jetsKT4_btag_secVertex, &b_jetsKT4_btag_secVertex);
		chain->SetBranchAddress("jetsKT4_chgEmE", &jetsKT4_chgEmE, &b_jetsKT4_chgEmE);
		chain->SetBranchAddress("jetsKT4_chgHadE", &jetsKT4_chgHadE, &b_jetsKT4_chgHadE);
		chain->SetBranchAddress("jetsKT4_chgMuE", &jetsKT4_chgMuE, &b_jetsKT4_chgMuE);
		chain->SetBranchAddress("jetsKT4_chg_Mult", &jetsKT4_chg_Mult, &b_jetsKT4_chg_Mult);
		chain->SetBranchAddress("jetsKT4_neutralEmE", &jetsKT4_neutralEmE, &b_jetsKT4_neutralEmE);
		chain->SetBranchAddress("jetsKT4_neutralHadE", &jetsKT4_neutralHadE, &b_jetsKT4_neutralHadE);
		chain->SetBranchAddress("jetsKT4_neutral_Mult", &jetsKT4_neutral_Mult, &b_jetsKT4_neutral_Mult);
		chain->SetBranchAddress("jetsKT4_mu_Mult", &jetsKT4_mu_Mult, &b_jetsKT4_mu_Mult);
		chain->SetBranchAddress("jetsKT4_emf", &jetsKT4_emf, &b_jetsKT4_emf);
		chain->SetBranchAddress("jetsKT4_ehf", &jetsKT4_ehf, &b_jetsKT4_ehf);
		chain->SetBranchAddress("jetsKT4_n60", &jetsKT4_n60, &b_jetsKT4_n60);
		chain->SetBranchAddress("jetsKT4_n90", &jetsKT4_n90, &b_jetsKT4_n90);
		chain->SetBranchAddress("jetsKT4_area", &jetsKT4_area, &b_jetsKT4_area);
		chain->SetBranchAddress("jetsKT4_mass", &jetsKT4_mass, &b_jetsKT4_mass);
	}

	if (m_jetAlgo == "SC5") {
		chain->SetBranchAddress("NjetsSC5", &NjetsSC5, &b_NjetsSC5);
		chain->SetBranchAddress("jetsSC5_energy", &jetsSC5_energy, &b_jetsSC5_energy);
		chain->SetBranchAddress("jetsSC5_et", &jetsSC5_et, &b_jetsSC5_et);
		chain->SetBranchAddress("jetsSC5_eta", &jetsSC5_eta, &b_jetsSC5_eta);
		chain->SetBranchAddress("jetsSC5_phi", &jetsSC5_phi, &b_jetsSC5_phi);
		chain->SetBranchAddress("jetsSC5_pt", &jetsSC5_pt, &b_jetsSC5_pt);
		chain->SetBranchAddress("jetsSC5_px", &jetsSC5_px, &b_jetsSC5_px);
		chain->SetBranchAddress("jetsSC5_py", &jetsSC5_py, &b_jetsSC5_py);
		chain->SetBranchAddress("jetsSC5_pz", &jetsSC5_pz, &b_jetsSC5_pz);
		chain->SetBranchAddress("jetsSC5_status", &jetsSC5_status, &b_jetsSC5_status);
		chain->SetBranchAddress("jetsSC5_theta", &jetsSC5_theta, &b_jetsSC5_theta);
		chain->SetBranchAddress("jetsSC5_btag_TC_highPur", &jetsSC5_btag_TC_highPur, &b_jetsSC5_btag_TC_highPur);
		chain->SetBranchAddress("jetsSC5_btag_TC_highEff", &jetsSC5_btag_TC_highEff, &b_jetsSC5_btag_TC_highEff);
		chain->SetBranchAddress("jetsSC5_btag_jetProb", &jetsSC5_btag_jetProb, &b_jetsSC5_btag_jetProb);
		chain->SetBranchAddress("jetsSC5_btag_jetBProb", &jetsSC5_btag_jetBProb, &b_jetsSC5_btag_jetBProb);
		chain->SetBranchAddress("jetsSC5_btag_softEle", &jetsSC5_btag_softEle, &b_jetsSC5_btag_softEle);
		chain->SetBranchAddress("jetsSC5_btag_softMuon", &jetsSC5_btag_softMuon, &b_jetsSC5_btag_softMuon);
		chain->SetBranchAddress("jetsSC5_btag_softMuonNoIP", &jetsSC5_btag_softMuonNoIP, &b_jetsSC5_btag_softMuonNoIP);
		chain->SetBranchAddress("jetsSC5_btag_secVertex", &jetsSC5_btag_secVertex, &b_jetsSC5_btag_secVertex);
		chain->SetBranchAddress("jetsSC5_chgEmE", &jetsSC5_chgEmE, &b_jetsSC5_chgEmE);
		chain->SetBranchAddress("jetsSC5_chgHadE", &jetsSC5_chgHadE, &b_jetsSC5_chgHadE);
		chain->SetBranchAddress("jetsSC5_chgMuE", &jetsSC5_chgMuE, &b_jetsSC5_chgMuE);
		chain->SetBranchAddress("jetsSC5_chg_Mult", &jetsSC5_chg_Mult, &b_jetsSC5_chg_Mult);
		chain->SetBranchAddress("jetsSC5_neutralEmE", &jetsSC5_neutralEmE, &b_jetsSC5_neutralEmE);
		chain->SetBranchAddress("jetsSC5_neutralHadE", &jetsSC5_neutralHadE, &b_jetsSC5_neutralHadE);
		chain->SetBranchAddress("jetsSC5_neutral_Mult", &jetsSC5_neutral_Mult, &b_jetsSC5_neutral_Mult);
		chain->SetBranchAddress("jetsSC5_mu_Mult", &jetsSC5_mu_Mult, &b_jetsSC5_mu_Mult);
		chain->SetBranchAddress("jetsSC5_emf", &jetsSC5_emf, &b_jetsSC5_emf);
		chain->SetBranchAddress("jetsSC5_ehf", &jetsSC5_ehf, &b_jetsSC5_ehf);
		chain->SetBranchAddress("jetsSC5_n60", &jetsSC5_n60, &b_jetsSC5_n60);
		chain->SetBranchAddress("jetsSC5_n90", &jetsSC5_n90, &b_jetsSC5_n90);
		chain->SetBranchAddress("jetsSC5_area", &jetsSC5_area, &b_jetsSC5_area);
		chain->SetBranchAddress("jetsSC5_mass", &jetsSC5_mass, &b_jetsSC5_mass);
	}
	chain->SetBranchAddress("Nmets", &Nmets, &b_Nmets);
	chain->SetBranchAddress("mets_et", &mets_et, &b_mets_et);
	chain->SetBranchAddress("mets_phi", &mets_phi, &b_mets_phi);
	chain->SetBranchAddress("mets_ex", &mets_ex, &b_mets_ex);
	chain->SetBranchAddress("mets_ey", &mets_ey, &b_mets_ey);
	chain->SetBranchAddress("mets_gen_et", &mets_gen_et, &b_mets_gen_et);
	chain->SetBranchAddress("mets_gen_phi", &mets_gen_phi, &b_mets_gen_phi);
	chain->SetBranchAddress("mets_sign", &mets_sign, &b_mets_sign);
	chain->SetBranchAddress("mets_sumEt", &mets_sumEt, &b_mets_sumEt);
	chain->SetBranchAddress("mets_unCPhi", &mets_unCPhi, &b_mets_unCPhi);
	chain->SetBranchAddress("mets_unCPt", &mets_unCPt, &b_mets_unCPt);
	chain->SetBranchAddress("mets_et_muonCor", &mets_et_muonCor, &b_mets_et_muonCor);
	chain->SetBranchAddress("mets_phi_muonCor", &mets_phi_muonCor, &b_mets_phi_muonCor);
	chain->SetBranchAddress("mets_et_JESCor", &mets_et_JESCor, &b_mets_et_JESCor);
	chain->SetBranchAddress("mets_phi_JESCor", &mets_phi_JESCor, &b_mets_phi_JESCor);

	if (m_metAlgo == "KT4") {
		chain->SetBranchAddress("NmetsKT4", &NmetsKT4, &b_NmetsKT4);
		chain->SetBranchAddress("metsKT4_et", &metsKT4_et, &b_metsKT4_et);
		chain->SetBranchAddress("metsKT4_phi", &metsKT4_phi, &b_metsKT4_phi);
		chain->SetBranchAddress("metsKT4_ex", &metsKT4_ex, &b_metsKT4_ex);
		chain->SetBranchAddress("metsKT4_ey", &metsKT4_ey, &b_metsKT4_ey);
		chain->SetBranchAddress("metsKT4_sumEt", &metsKT4_sumEt, &b_metsKT4_sumEt);
		chain->SetBranchAddress("metsKT4_et_JESCor", &metsKT4_et_JESCor, &b_metsKT4_et_JESCor);
		chain->SetBranchAddress("metsKT4_phi_JESCor", &metsKT4_phi_JESCor, &b_metsKT4_phi_JESCor);
	}
	if (m_metAlgo == "SC5") {
		chain->SetBranchAddress("NmetsSC5", &NmetsSC5, &b_NmetsSC5);
		chain->SetBranchAddress("metsSC5_et", &metsSC5_et, &b_metsSC5_et);
		chain->SetBranchAddress("metsSC5_phi", &metsSC5_phi, &b_metsSC5_phi);
		chain->SetBranchAddress("metsSC5_ex", &metsSC5_ex, &b_metsSC5_ex);
		chain->SetBranchAddress("metsSC5_ey", &metsSC5_ey, &b_metsSC5_ey);
		chain->SetBranchAddress("metsSC5_sumEt", &metsSC5_sumEt, &b_metsSC5_sumEt);
		chain->SetBranchAddress("metsSC5_et_JESCor", &metsSC5_et_JESCor, &b_metsSC5_et_JESCor);
		chain->SetBranchAddress("metsSC5_phi_JESCor", &metsSC5_phi_JESCor, &b_metsSC5_phi_JESCor);
	}
	if (m_metAlgo == "tcmet") {
		chain->SetBranchAddress("Ntcmets", &Ntcmets, &b_Ntcmets);
		chain->SetBranchAddress("tcmets_et", &tcmets_et, &b_tcmets_et);
		chain->SetBranchAddress("tcmets_phi", &tcmets_phi, &b_tcmets_phi);
		chain->SetBranchAddress("tcmets_ex", &tcmets_ex, &b_tcmets_ex);
		chain->SetBranchAddress("tcmets_ey", &tcmets_ey, &b_tcmets_ey);
		chain->SetBranchAddress("tcmets_sumEt", &tcmets_sumEt, &b_tcmets_sumEt);
		chain->SetBranchAddress("tcmets_et_muonCor", &tcmets_et_muonCor, &b_tcmets_et_muonCor);
		chain->SetBranchAddress("tcmets_phi_muonCor", &tcmets_phi_muonCor, &b_tcmets_phi_muonCor);
	}

	chain->SetBranchAddress("Nmus", &Nmus, &b_Nmus);
	chain->SetBranchAddress("mus_energy", &mus_energy, &b_mus_energy);
	chain->SetBranchAddress("mus_et", &mus_et, &b_mus_et);
	chain->SetBranchAddress("mus_eta", &mus_eta, &b_mus_eta);
	chain->SetBranchAddress("mus_phi", &mus_phi, &b_mus_phi);
	chain->SetBranchAddress("mus_pt", &mus_pt, &b_mus_pt);
	chain->SetBranchAddress("mus_px", &mus_px, &b_mus_px);
	chain->SetBranchAddress("mus_py", &mus_py, &b_mus_py);
	chain->SetBranchAddress("mus_pz", &mus_pz, &b_mus_pz);
	chain->SetBranchAddress("mus_status", &mus_status, &b_mus_status);
	chain->SetBranchAddress("mus_theta", &mus_theta, &b_mus_theta);
	chain->SetBranchAddress("mus_gen_id", &mus_gen_id, &b_mus_gen_id);
	chain->SetBranchAddress("mus_gen_phi", &mus_gen_phi, &b_mus_gen_phi);
	chain->SetBranchAddress("mus_gen_pt", &mus_gen_pt, &b_mus_gen_pt);
	chain->SetBranchAddress("mus_gen_pz", &mus_gen_pz, &b_mus_gen_pz);
	chain->SetBranchAddress("mus_gen_px", &mus_gen_px, &b_mus_gen_px);
	chain->SetBranchAddress("mus_gen_py", &mus_gen_py, &b_mus_gen_py);
	chain->SetBranchAddress("mus_gen_eta", &mus_gen_eta, &b_mus_gen_eta);
	chain->SetBranchAddress("mus_gen_theta", &mus_gen_theta, &b_mus_gen_theta);
	chain->SetBranchAddress("mus_gen_et", &mus_gen_et, &b_mus_gen_et);
	chain->SetBranchAddress("mus_gen_mother_id", &mus_gen_mother_id, &b_mus_gen_mother_id);
	chain->SetBranchAddress("mus_gen_mother_phi", &mus_gen_mother_phi, &b_mus_gen_mother_phi);
	chain->SetBranchAddress("mus_gen_mother_pt", &mus_gen_mother_pt, &b_mus_gen_mother_pt);
	chain->SetBranchAddress("mus_gen_mother_pz", &mus_gen_mother_pz, &b_mus_gen_mother_pz);
	chain->SetBranchAddress("mus_gen_mother_px", &mus_gen_mother_px, &b_mus_gen_mother_px);
	chain->SetBranchAddress("mus_gen_mother_py", &mus_gen_mother_py, &b_mus_gen_mother_py);
	chain->SetBranchAddress("mus_gen_mother_eta", &mus_gen_mother_eta, &b_mus_gen_mother_eta);
	chain->SetBranchAddress("mus_gen_mother_theta", &mus_gen_mother_theta, &b_mus_gen_mother_theta);
	chain->SetBranchAddress("mus_gen_mother_et", &mus_gen_mother_et, &b_mus_gen_mother_et);
	chain->SetBranchAddress("mus_tkHits", &mus_tkHits, &b_mus_tkHits);
	chain->SetBranchAddress("mus_cIso", &mus_cIso, &b_mus_cIso);
	chain->SetBranchAddress("mus_tIso", &mus_tIso, &b_mus_tIso);
	chain->SetBranchAddress("mus_ecalIso", &mus_ecalIso, &b_mus_ecalIso);
	chain->SetBranchAddress("mus_hcalIso", &mus_hcalIso, &b_mus_hcalIso);
	chain->SetBranchAddress("mus_ecalvetoDep", &mus_ecalvetoDep, &b_mus_ecalvetoDep);
	chain->SetBranchAddress("mus_hcalvetoDep", &mus_hcalvetoDep, &b_mus_hcalvetoDep);
	chain->SetBranchAddress("mus_calEnergyEm", &mus_calEnergyEm, &b_mus_calEnergyEm);
	chain->SetBranchAddress("mus_calEnergyHad", &mus_calEnergyHad, &b_mus_calEnergyHad);
	chain->SetBranchAddress("mus_calEnergyHo", &mus_calEnergyHo, &b_mus_calEnergyHo);
	chain->SetBranchAddress("mus_calEnergyEmS9", &mus_calEnergyEmS9, &b_mus_calEnergyEmS9);
	chain->SetBranchAddress("mus_calEnergyHadS9", &mus_calEnergyHadS9, &b_mus_calEnergyHadS9);
	chain->SetBranchAddress("mus_calEnergyHoS9", &mus_calEnergyHoS9, &b_mus_calEnergyHoS9);
	chain->SetBranchAddress("mus_iso03_sumPt", &mus_iso03_sumPt, &b_mus_iso03_sumPt);
	chain->SetBranchAddress("mus_iso03_emEt", &mus_iso03_emEt, &b_mus_iso03_emEt);
	chain->SetBranchAddress("mus_iso03_hadEt", &mus_iso03_hadEt, &b_mus_iso03_hadEt);
	chain->SetBranchAddress("mus_iso03_hoEt", &mus_iso03_hoEt, &b_mus_iso03_hoEt);
	chain->SetBranchAddress("mus_iso03_nTracks", &mus_iso03_nTracks, &b_mus_iso03_nTracks);
	chain->SetBranchAddress("mus_iso05_sumPt", &mus_iso05_sumPt, &b_mus_iso05_sumPt);
	chain->SetBranchAddress("mus_iso05_emEt", &mus_iso05_emEt, &b_mus_iso05_emEt);
	chain->SetBranchAddress("mus_iso05_hadEt", &mus_iso05_hadEt, &b_mus_iso05_hadEt);
	chain->SetBranchAddress("mus_iso05_hoEt", &mus_iso05_hoEt, &b_mus_iso05_hoEt);
	chain->SetBranchAddress("mus_iso05_nTracks", &mus_iso05_nTracks, &b_mus_iso05_nTracks);
	chain->SetBranchAddress("mus_charge", &mus_charge, &b_mus_charge);
	chain->SetBranchAddress("mus_cm_chi2", &mus_cm_chi2, &b_mus_cm_chi2);
	chain->SetBranchAddress("mus_cm_ndof", &mus_cm_ndof, &b_mus_cm_ndof);
	chain->SetBranchAddress("mus_cm_chg", &mus_cm_chg, &b_mus_cm_chg);
	chain->SetBranchAddress("mus_cm_pt", &mus_cm_pt, &b_mus_cm_pt);
	chain->SetBranchAddress("mus_cm_px", &mus_cm_px, &b_mus_cm_px);
	chain->SetBranchAddress("mus_cm_py", &mus_cm_py, &b_mus_cm_py);
	chain->SetBranchAddress("mus_cm_pz", &mus_cm_pz, &b_mus_cm_pz);
	chain->SetBranchAddress("mus_cm_eta", &mus_cm_eta, &b_mus_cm_eta);
	chain->SetBranchAddress("mus_cm_phi", &mus_cm_phi, &b_mus_cm_phi);
	chain->SetBranchAddress("mus_cm_theta", &mus_cm_theta, &b_mus_cm_theta);
	chain->SetBranchAddress("mus_cm_d0dum", &mus_cm_d0dum, &b_mus_cm_d0dum);
	chain->SetBranchAddress("mus_cm_dz", &mus_cm_dz, &b_mus_cm_dz);
	chain->SetBranchAddress("mus_cm_vx", &mus_cm_vx, &b_mus_cm_vx);
	chain->SetBranchAddress("mus_cm_vy", &mus_cm_vy, &b_mus_cm_vy);
	chain->SetBranchAddress("mus_cm_vz", &mus_cm_vz, &b_mus_cm_vz);
	chain->SetBranchAddress("mus_cm_numvalhits", &mus_cm_numvalhits, &b_mus_cm_numvalhits);
	chain->SetBranchAddress("mus_cm_numlosthits", &mus_cm_numlosthits, &b_mus_cm_numlosthits);
	chain->SetBranchAddress("mus_cm_d0dumErr", &mus_cm_d0dumErr, &b_mus_cm_d0dumErr);
	chain->SetBranchAddress("mus_cm_dzErr", &mus_cm_dzErr, &b_mus_cm_dzErr);
	chain->SetBranchAddress("mus_cm_ptErr", &mus_cm_ptErr, &b_mus_cm_ptErr);
	chain->SetBranchAddress("mus_cm_etaErr", &mus_cm_etaErr, &b_mus_cm_etaErr);
	chain->SetBranchAddress("mus_cm_phiErr", &mus_cm_phiErr, &b_mus_cm_phiErr);
	chain->SetBranchAddress("mus_tk_chi2", &mus_tk_chi2, &b_mus_tk_chi2);
	chain->SetBranchAddress("mus_tk_ndof", &mus_tk_ndof, &b_mus_tk_ndof);
	chain->SetBranchAddress("mus_tk_chg", &mus_tk_chg, &b_mus_tk_chg);
	chain->SetBranchAddress("mus_tk_pt", &mus_tk_pt, &b_mus_tk_pt);
	chain->SetBranchAddress("mus_tk_px", &mus_tk_px, &b_mus_tk_px);
	chain->SetBranchAddress("mus_tk_py", &mus_tk_py, &b_mus_tk_py);
	chain->SetBranchAddress("mus_tk_pz", &mus_tk_pz, &b_mus_tk_pz);
	chain->SetBranchAddress("mus_tk_eta", &mus_tk_eta, &b_mus_tk_eta);
	chain->SetBranchAddress("mus_tk_phi", &mus_tk_phi, &b_mus_tk_phi);
	chain->SetBranchAddress("mus_tk_theta", &mus_tk_theta, &b_mus_tk_theta);
	chain->SetBranchAddress("mus_tk_d0dum", &mus_tk_d0dum, &b_mus_tk_d0dum);
	chain->SetBranchAddress("mus_tk_dz", &mus_tk_dz, &b_mus_tk_dz);
	chain->SetBranchAddress("mus_tk_vx", &mus_tk_vx, &b_mus_tk_vx);
	chain->SetBranchAddress("mus_tk_vy", &mus_tk_vy, &b_mus_tk_vy);
	chain->SetBranchAddress("mus_tk_vz", &mus_tk_vz, &b_mus_tk_vz);
	chain->SetBranchAddress("mus_tk_numvalhits", &mus_tk_numvalhits, &b_mus_tk_numvalhits);
	chain->SetBranchAddress("mus_tk_numlosthits", &mus_tk_numlosthits, &b_mus_tk_numlosthits);
	chain->SetBranchAddress("mus_tk_d0dumErr", &mus_tk_d0dumErr, &b_mus_tk_d0dumErr);
	chain->SetBranchAddress("mus_tk_dzErr", &mus_tk_dzErr, &b_mus_tk_dzErr);
	chain->SetBranchAddress("mus_tk_ptErr", &mus_tk_ptErr, &b_mus_tk_ptErr);
	chain->SetBranchAddress("mus_tk_etaErr", &mus_tk_etaErr, &b_mus_tk_etaErr);
	chain->SetBranchAddress("mus_tk_phiErr", &mus_tk_phiErr, &b_mus_tk_phiErr);
	chain->SetBranchAddress("mus_stamu_chi2", &mus_stamu_chi2, &b_mus_stamu_chi2);
	chain->SetBranchAddress("mus_stamu_ndof", &mus_stamu_ndof, &b_mus_stamu_ndof);
	chain->SetBranchAddress("mus_stamu_chg", &mus_stamu_chg, &b_mus_stamu_chg);
	chain->SetBranchAddress("mus_stamu_pt", &mus_stamu_pt, &b_mus_stamu_pt);
	chain->SetBranchAddress("mus_stamu_px", &mus_stamu_px, &b_mus_stamu_px);
	chain->SetBranchAddress("mus_stamu_py", &mus_stamu_py, &b_mus_stamu_py);
	chain->SetBranchAddress("mus_stamu_pz", &mus_stamu_pz, &b_mus_stamu_pz);
	chain->SetBranchAddress("mus_stamu_eta", &mus_stamu_eta, &b_mus_stamu_eta);
	chain->SetBranchAddress("mus_stamu_phi", &mus_stamu_phi, &b_mus_stamu_phi);
	chain->SetBranchAddress("mus_stamu_theta", &mus_stamu_theta, &b_mus_stamu_theta);
	chain->SetBranchAddress("mus_stamu_d0dum", &mus_stamu_d0dum, &b_mus_stamu_d0dum);
	chain->SetBranchAddress("mus_stamu_dz", &mus_stamu_dz, &b_mus_stamu_dz);
	chain->SetBranchAddress("mus_stamu_vx", &mus_stamu_vx, &b_mus_stamu_vx);
	chain->SetBranchAddress("mus_stamu_vy", &mus_stamu_vy, &b_mus_stamu_vy);
	chain->SetBranchAddress("mus_stamu_vz", &mus_stamu_vz, &b_mus_stamu_vz);
	chain->SetBranchAddress("mus_stamu_numvalhits", &mus_stamu_numvalhits, &b_mus_stamu_numvalhits);
	chain->SetBranchAddress("mus_stamu_numlosthits", &mus_stamu_numlosthits, &b_mus_stamu_numlosthits);
	chain->SetBranchAddress("mus_stamu_d0dumErr", &mus_stamu_d0dumErr, &b_mus_stamu_d0dumErr);
	chain->SetBranchAddress("mus_stamu_dzErr", &mus_stamu_dzErr, &b_mus_stamu_dzErr);
	chain->SetBranchAddress("mus_stamu_ptErr", &mus_stamu_ptErr, &b_mus_stamu_ptErr);
	chain->SetBranchAddress("mus_stamu_etaErr", &mus_stamu_etaErr, &b_mus_stamu_etaErr);
	chain->SetBranchAddress("mus_stamu_phiErr", &mus_stamu_phiErr, &b_mus_stamu_phiErr);
	chain->SetBranchAddress("mus_num_matches", &mus_num_matches, &b_mus_num_matches);
	chain->SetBranchAddress("mus_id_All", &mus_id_All, &b_mus_id_All);
	chain->SetBranchAddress("mus_id_AllGlobalMuons", &mus_id_AllGlobalMuons, &b_mus_id_AllGlobalMuons);
	chain->SetBranchAddress("mus_id_AllStandAloneMuons", &mus_id_AllStandAloneMuons, &b_mus_id_AllStandAloneMuons);
	chain->SetBranchAddress("mus_id_AllTrackerMuons", &mus_id_AllTrackerMuons, &b_mus_id_AllTrackerMuons);
	chain->SetBranchAddress("mus_id_TrackerMuonArbitrated", &mus_id_TrackerMuonArbitrated, &b_mus_id_TrackerMuonArbitrated);
	chain->SetBranchAddress("mus_id_AllArbitrated", &mus_id_AllArbitrated, &b_mus_id_AllArbitrated);
	chain->SetBranchAddress("mus_id_GlobalMuonPromptTight", &mus_id_GlobalMuonPromptTight, &b_mus_id_GlobalMuonPromptTight);
	chain->SetBranchAddress("mus_id_TMLastStationLoose", &mus_id_TMLastStationLoose, &b_mus_id_TMLastStationLoose);
	chain->SetBranchAddress("mus_id_TMLastStationTight", &mus_id_TMLastStationTight, &b_mus_id_TMLastStationTight);
	chain->SetBranchAddress("mus_id_TM2DCompatibilityLoose", &mus_id_TM2DCompatibilityLoose, &b_mus_id_TM2DCompatibilityLoose);
	chain->SetBranchAddress("mus_id_TM2DCompatibilityTight", &mus_id_TM2DCompatibilityTight, &b_mus_id_TM2DCompatibilityTight);
	chain->SetBranchAddress("mus_id_TMOneStationLoose", &mus_id_TMOneStationLoose, &b_mus_id_TMOneStationLoose);
	chain->SetBranchAddress("mus_id_TMOneStationTight", &mus_id_TMOneStationTight, &b_mus_id_TMOneStationTight);
	chain->SetBranchAddress("mus_id_TMLastStationOptimizedLowPtLoose", &mus_id_TMLastStationOptimizedLowPtLoose,
			&b_mus_id_TMLastStationOptimizedLowPtLoose);
	chain->SetBranchAddress("mus_id_TMLastStationOptimizedLowPtTight", &mus_id_TMLastStationOptimizedLowPtTight,
			&b_mus_id_TMLastStationOptimizedLowPtTight);
	chain->SetBranchAddress("Nphotons", &Nphotons, &b_Nphotons);
	chain->SetBranchAddress("photons_energy", &photons_energy, &b_photons_energy);
	chain->SetBranchAddress("photons_et", &photons_et, &b_photons_et);
	chain->SetBranchAddress("photons_eta", &photons_eta, &b_photons_eta);
	chain->SetBranchAddress("photons_phi", &photons_phi, &b_photons_phi);
	chain->SetBranchAddress("photons_pt", &photons_pt, &b_photons_pt);
	chain->SetBranchAddress("photons_px", &photons_px, &b_photons_px);
	chain->SetBranchAddress("photons_py", &photons_py, &b_photons_py);
	chain->SetBranchAddress("photons_pz", &photons_pz, &b_photons_pz);
	chain->SetBranchAddress("photons_status", &photons_status, &b_photons_status);
	chain->SetBranchAddress("photons_theta", &photons_theta, &b_photons_theta);
	chain->SetBranchAddress("photons_hadOverEM", &photons_hadOverEM, &b_photons_hadOverEM);
	chain->SetBranchAddress("photons_scEnergy", &photons_scEnergy, &b_photons_scEnergy);
	chain->SetBranchAddress("photons_scRawEnergy", &photons_scRawEnergy, &b_photons_scRawEnergy);
	chain->SetBranchAddress("photons_scEta", &photons_scEta, &b_photons_scEta);
	chain->SetBranchAddress("photons_scPhi", &photons_scPhi, &b_photons_scPhi);
	chain->SetBranchAddress("photons_scEtaWidth", &photons_scEtaWidth, &b_photons_scEtaWidth);
	chain->SetBranchAddress("photons_scPhiWidth", &photons_scPhiWidth, &b_photons_scPhiWidth);
	chain->SetBranchAddress("photons_tIso", &photons_tIso, &b_photons_tIso);
	chain->SetBranchAddress("photons_ecalIso", &photons_ecalIso, &b_photons_ecalIso);
	chain->SetBranchAddress("photons_hcalIso", &photons_hcalIso, &b_photons_hcalIso);
	chain->SetBranchAddress("photons_isoEcalRecHitDR04", &photons_isoEcalRecHitDR04, &b_photons_isoEcalRecHitDR04);
	chain->SetBranchAddress("photons_isoHcalRecHitDR04", &photons_isoHcalRecHitDR04, &b_photons_isoHcalRecHitDR04);
	chain->SetBranchAddress("photons_isoSolidTrkConeDR04", &photons_isoSolidTrkConeDR04, &b_photons_isoSolidTrkConeDR04);
	chain->SetBranchAddress("photons_isoHollowTrkConeDR04", &photons_isoHollowTrkConeDR04, &b_photons_isoHollowTrkConeDR04);
	chain->SetBranchAddress("photons_nTrkSolidConeDR04", &photons_nTrkSolidConeDR04, &b_photons_nTrkSolidConeDR04);
	chain->SetBranchAddress("photons_nTrkHollowConeDR04", &photons_nTrkHollowConeDR04, &b_photons_nTrkHollowConeDR04);
	chain->SetBranchAddress("photons_isoEcalRecHitDR03", &photons_isoEcalRecHitDR03, &b_photons_isoEcalRecHitDR03);
	chain->SetBranchAddress("photons_isoHcalRecHitDR03", &photons_isoHcalRecHitDR03, &b_photons_isoHcalRecHitDR03);
	chain->SetBranchAddress("photons_isoSolidTrkConeDR03", &photons_isoSolidTrkConeDR03, &b_photons_isoSolidTrkConeDR03);
	chain->SetBranchAddress("photons_isoHollowTrkConeDR03", &photons_isoHollowTrkConeDR03, &b_photons_isoHollowTrkConeDR03);
	chain->SetBranchAddress("photons_nTrkSolidConeDR03", &photons_nTrkSolidConeDR03, &b_photons_nTrkSolidConeDR03);
	chain->SetBranchAddress("photons_nTrkHollowConeDR03", &photons_nTrkHollowConeDR03, &b_photons_nTrkHollowConeDR03);
	chain->SetBranchAddress("photons_isAlsoElectron", &photons_isAlsoElectron, &b_photons_isAlsoElectron);
	chain->SetBranchAddress("photons_hasPixelSeed", &photons_hasPixelSeed, &b_photons_hasPixelSeed);
	chain->SetBranchAddress("photons_isConverted", &photons_isConverted, &b_photons_isConverted);
	chain->SetBranchAddress("photons_isEBGap", &photons_isEBGap, &b_photons_isEBGap);
	chain->SetBranchAddress("photons_isEEGap", &photons_isEEGap, &b_photons_isEEGap);
	chain->SetBranchAddress("photons_isEBEEGap", &photons_isEBEEGap, &b_photons_isEBEEGap);
	chain->SetBranchAddress("photons_isEBPho", &photons_isEBPho, &b_photons_isEBPho);
	chain->SetBranchAddress("photons_isEEPho", &photons_isEEPho, &b_photons_isEEPho);
	chain->SetBranchAddress("photons_isLoosePhoton", &photons_isLoosePhoton, &b_photons_isLoosePhoton);
	chain->SetBranchAddress("photons_isTightPhoton", &photons_isTightPhoton, &b_photons_isTightPhoton);
	chain->SetBranchAddress("photons_r9", &photons_r9, &b_photons_r9);
	chain->SetBranchAddress("photons_gen_et", &photons_gen_et, &b_photons_gen_et);
	chain->SetBranchAddress("photons_gen_eta", &photons_gen_eta, &b_photons_gen_eta);
	chain->SetBranchAddress("photons_gen_phi", &photons_gen_phi, &b_photons_gen_phi);
	chain->SetBranchAddress("photons_gen_id", &photons_gen_id, &b_photons_gen_id);
	chain->SetBranchAddress("Npv", &Npv, &b_Npv);
	chain->SetBranchAddress("pv_x", &pv_x, &b_pv_x);
	chain->SetBranchAddress("pv_y", &pv_y, &b_pv_y);
	chain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
	chain->SetBranchAddress("pv_xErr", &pv_xErr, &b_pv_xErr);
	chain->SetBranchAddress("pv_yErr", &pv_yErr, &b_pv_yErr);
	chain->SetBranchAddress("pv_zErr", &pv_zErr, &b_pv_zErr);
	chain->SetBranchAddress("Ntracks", &Ntracks, &b_Ntracks);
	chain->SetBranchAddress("tracks_chi2", &tracks_chi2, &b_tracks_chi2);
	chain->SetBranchAddress("tracks_ndof", &tracks_ndof, &b_tracks_ndof);
	chain->SetBranchAddress("tracks_chg", &tracks_chg, &b_tracks_chg);
	chain->SetBranchAddress("tracks_pt", &tracks_pt, &b_tracks_pt);
	chain->SetBranchAddress("tracks_px", &tracks_px, &b_tracks_px);
	chain->SetBranchAddress("tracks_py", &tracks_py, &b_tracks_py);
	chain->SetBranchAddress("tracks_pz", &tracks_pz, &b_tracks_pz);
	chain->SetBranchAddress("tracks_eta", &tracks_eta, &b_tracks_eta);
	chain->SetBranchAddress("tracks_phi", &tracks_phi, &b_tracks_phi);
	chain->SetBranchAddress("tracks_theta", &tracks_theta, &b_tracks_theta);
	chain->SetBranchAddress("tracks_d0dum", &tracks_d0dum, &b_tracks_d0dum);
	chain->SetBranchAddress("tracks_dz", &tracks_dz, &b_tracks_dz);
	chain->SetBranchAddress("tracks_vx", &tracks_vx, &b_tracks_vx);
	chain->SetBranchAddress("tracks_vy", &tracks_vy, &b_tracks_vy);
	chain->SetBranchAddress("tracks_vz", &tracks_vz, &b_tracks_vz);
	chain->SetBranchAddress("tracks_numvalhits", &tracks_numvalhits, &b_tracks_numvalhits);
	chain->SetBranchAddress("tracks_numlosthits", &tracks_numlosthits, &b_tracks_numlosthits);
	chain->SetBranchAddress("tracks_d0dumErr", &tracks_d0dumErr, &b_tracks_d0dumErr);
	chain->SetBranchAddress("tracks_dzErr", &tracks_dzErr, &b_tracks_dzErr);
	chain->SetBranchAddress("tracks_ptErr", &tracks_ptErr, &b_tracks_ptErr);
	chain->SetBranchAddress("tracks_etaErr", &tracks_etaErr, &b_tracks_etaErr);
	chain->SetBranchAddress("tracks_phiErr", &tracks_phiErr, &b_tracks_phiErr);
	chain->SetBranchAddress("tracks_Nrechits", &tracks_Nrechits, &b_tracks_Nrechits);
	chain->SetBranchAddress("tracks_innerHitX", &tracks_innerHitX, &b_tracks_innerHitX);
	chain->SetBranchAddress("tracks_innerHitY", &tracks_innerHitY, &b_tracks_innerHitY);
	chain->SetBranchAddress("tracks_innerHitZ", &tracks_innerHitZ, &b_tracks_innerHitZ);
	chain->SetBranchAddress("tracks_outerHitX", &tracks_outerHitX, &b_tracks_outerHitX);
	chain->SetBranchAddress("tracks_outerHitY", &tracks_outerHitY, &b_tracks_outerHitY);
	chain->SetBranchAddress("tracks_outerHitZ", &tracks_outerHitZ, &b_tracks_outerHitZ);
	chain->SetBranchAddress("tracks_highPurity", &tracks_highPurity, &b_tracks_highPurity);
	if (m_useMisslayers)
		chain->SetBranchAddress("tracks_innerLayerMissingHits", &tracks_innerLayerMissingHits, &b_tracks_innerLayerMissingHits);
	chain->SetBranchAddress("run", &run, &b_run);
	chain->SetBranchAddress("event", &event, &b_event);
	chain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);

	// HLT tree (nonIso ele trigger)
	//chain2->SetBranchAddress("HLT_Ele10_SW_L1R", &HLT_Ele10_SW_L1R, &b_HLT_Ele10_SW_L1R); //(8E29, startup)
	//cout << "GetTrigger()" << GetTrigger()<< endl;
	if (GetTrigger()) {
		if (HLTBit == "HLT_Ele15_LW_L1R")
			chain2->SetBranchAddress("HLT_Ele15_LW_L1R", &HLT_Ele15_LW_L1R, &b_HLT_Ele15_LW_L1R); //v2 (1E31, ideal)
		if (HLTBit == "HLT_Ele15_SW_L1R")
			chain2->SetBranchAddress("HLT_Ele15_SW_L1R", &HLT_Ele15_SW_L1R, &b_HLT_Ele15_SW_L1R);
	}
	///------------------------  MC Truth info  ------------------------------------
	if (!IsData()) {
		cout << " Set MC branch address" << endl;
		chain->SetBranchAddress("Nmc_doc", &Nmc_doc, &b_Nmc_doc);
		chain->SetBranchAddress("mc_doc_id", &mc_doc_id, &b_mc_doc_id);
		chain->SetBranchAddress("mc_doc_pt", &mc_doc_pt, &b_mc_doc_pt);
		chain->SetBranchAddress("mc_doc_px", &mc_doc_px, &b_mc_doc_px);
		chain->SetBranchAddress("mc_doc_py", &mc_doc_py, &b_mc_doc_py);
		chain->SetBranchAddress("mc_doc_pz", &mc_doc_pz, &b_mc_doc_pz);
		chain->SetBranchAddress("mc_doc_eta", &mc_doc_eta, &b_mc_doc_eta);
		chain->SetBranchAddress("mc_doc_phi", &mc_doc_phi, &b_mc_doc_phi);
		chain->SetBranchAddress("mc_doc_theta", &mc_doc_theta, &b_mc_doc_theta);
		chain->SetBranchAddress("mc_doc_energy", &mc_doc_energy, &b_mc_doc_energy);
		chain->SetBranchAddress("mc_doc_status", &mc_doc_status, &b_mc_doc_status);
		chain->SetBranchAddress("mc_doc_charge", &mc_doc_charge, &b_mc_doc_charge);
		chain->SetBranchAddress("mc_doc_mother_id", &mc_doc_mother_id, &b_mc_doc_mother_id);
		chain->SetBranchAddress("mc_doc_grandmother_id", &mc_doc_grandmother_id, &b_mc_doc_grandmother_id);
		chain->SetBranchAddress("mc_doc_ggrandmother_id", &mc_doc_ggrandmother_id, &b_mc_doc_ggrandmother_id);
		chain->SetBranchAddress("mc_doc_mother_pt", &mc_doc_mother_pt, &b_mc_doc_mother_pt);
		chain->SetBranchAddress("mc_doc_vertex_x", &mc_doc_vertex_x, &b_mc_doc_vertex_x);
		chain->SetBranchAddress("mc_doc_vertex_y", &mc_doc_vertex_y, &b_mc_doc_vertex_y);
		chain->SetBranchAddress("mc_doc_vertex_z", &mc_doc_vertex_z, &b_mc_doc_vertex_z);
		chain->SetBranchAddress("mc_doc_mass", &mc_doc_mass, &b_mc_doc_mass);
		chain->SetBranchAddress("mc_doc_numOfDaughters", &mc_doc_numOfDaughters, &b_mc_doc_numOfDaughters);
		chain->SetBranchAddress("mc_doc_numOfMothers", &mc_doc_numOfMothers, &b_mc_doc_numOfMothers);

		if (m_studyPDFunc) {
			cout << "Set PDF branch address" << endl;
			chain2->SetBranchAddress("PDFWcteq66_0", &PDFWcteq66_0, &b_PDFWcteq66_0);
			chain2->SetBranchAddress("PDFWcteq66_1", &PDFWcteq66_1, &b_PDFWcteq66_1);
			chain2->SetBranchAddress("PDFWcteq66_10", &PDFWcteq66_10, &b_PDFWcteq66_10);
			chain2->SetBranchAddress("PDFWcteq66_11", &PDFWcteq66_11, &b_PDFWcteq66_11);
			chain2->SetBranchAddress("PDFWcteq66_12", &PDFWcteq66_12, &b_PDFWcteq66_12);
			chain2->SetBranchAddress("PDFWcteq66_13", &PDFWcteq66_13, &b_PDFWcteq66_13);
			chain2->SetBranchAddress("PDFWcteq66_14", &PDFWcteq66_14, &b_PDFWcteq66_14);
			chain2->SetBranchAddress("PDFWcteq66_15", &PDFWcteq66_15, &b_PDFWcteq66_15);
			chain2->SetBranchAddress("PDFWcteq66_16", &PDFWcteq66_16, &b_PDFWcteq66_16);
			chain2->SetBranchAddress("PDFWcteq66_17", &PDFWcteq66_17, &b_PDFWcteq66_17);
			chain2->SetBranchAddress("PDFWcteq66_18", &PDFWcteq66_18, &b_PDFWcteq66_18);
			chain2->SetBranchAddress("PDFWcteq66_19", &PDFWcteq66_19, &b_PDFWcteq66_19);
			chain2->SetBranchAddress("PDFWcteq66_2", &PDFWcteq66_2, &b_PDFWcteq66_2);
			chain2->SetBranchAddress("PDFWcteq66_20", &PDFWcteq66_20, &b_PDFWcteq66_20);
			chain2->SetBranchAddress("PDFWcteq66_21", &PDFWcteq66_21, &b_PDFWcteq66_21);
			chain2->SetBranchAddress("PDFWcteq66_22", &PDFWcteq66_22, &b_PDFWcteq66_22);
			chain2->SetBranchAddress("PDFWcteq66_23", &PDFWcteq66_23, &b_PDFWcteq66_23);
			chain2->SetBranchAddress("PDFWcteq66_24", &PDFWcteq66_24, &b_PDFWcteq66_24);
			chain2->SetBranchAddress("PDFWcteq66_25", &PDFWcteq66_25, &b_PDFWcteq66_25);
			chain2->SetBranchAddress("PDFWcteq66_26", &PDFWcteq66_26, &b_PDFWcteq66_26);
			chain2->SetBranchAddress("PDFWcteq66_27", &PDFWcteq66_27, &b_PDFWcteq66_27);
			chain2->SetBranchAddress("PDFWcteq66_28", &PDFWcteq66_28, &b_PDFWcteq66_28);
			chain2->SetBranchAddress("PDFWcteq66_29", &PDFWcteq66_29, &b_PDFWcteq66_29);
			chain2->SetBranchAddress("PDFWcteq66_3", &PDFWcteq66_3, &b_PDFWcteq66_3);
			chain2->SetBranchAddress("PDFWcteq66_30", &PDFWcteq66_30, &b_PDFWcteq66_30);
			chain2->SetBranchAddress("PDFWcteq66_31", &PDFWcteq66_31, &b_PDFWcteq66_31);
			chain2->SetBranchAddress("PDFWcteq66_32", &PDFWcteq66_32, &b_PDFWcteq66_32);
			chain2->SetBranchAddress("PDFWcteq66_33", &PDFWcteq66_33, &b_PDFWcteq66_33);
			chain2->SetBranchAddress("PDFWcteq66_34", &PDFWcteq66_34, &b_PDFWcteq66_34);
			chain2->SetBranchAddress("PDFWcteq66_35", &PDFWcteq66_35, &b_PDFWcteq66_35);
			chain2->SetBranchAddress("PDFWcteq66_36", &PDFWcteq66_36, &b_PDFWcteq66_36);
			chain2->SetBranchAddress("PDFWcteq66_37", &PDFWcteq66_37, &b_PDFWcteq66_37);
			chain2->SetBranchAddress("PDFWcteq66_38", &PDFWcteq66_38, &b_PDFWcteq66_38);
			chain2->SetBranchAddress("PDFWcteq66_39", &PDFWcteq66_39, &b_PDFWcteq66_39);
			chain2->SetBranchAddress("PDFWcteq66_4", &PDFWcteq66_4, &b_PDFWcteq66_4);
			chain2->SetBranchAddress("PDFWcteq66_40", &PDFWcteq66_40, &b_PDFWcteq66_40);
			chain2->SetBranchAddress("PDFWcteq66_41", &PDFWcteq66_41, &b_PDFWcteq66_41);
			chain2->SetBranchAddress("PDFWcteq66_42", &PDFWcteq66_42, &b_PDFWcteq66_42);
			chain2->SetBranchAddress("PDFWcteq66_43", &PDFWcteq66_43, &b_PDFWcteq66_43);
			chain2->SetBranchAddress("PDFWcteq66_44", &PDFWcteq66_44, &b_PDFWcteq66_44);
			chain2->SetBranchAddress("PDFWcteq66_5", &PDFWcteq66_5, &b_PDFWcteq66_5);
			chain2->SetBranchAddress("PDFWcteq66_6", &PDFWcteq66_6, &b_PDFWcteq66_6);
			chain2->SetBranchAddress("PDFWcteq66_7", &PDFWcteq66_7, &b_PDFWcteq66_7);
			chain2->SetBranchAddress("PDFWcteq66_8", &PDFWcteq66_8, &b_PDFWcteq66_8);
			chain2->SetBranchAddress("PDFWcteq66_9", &PDFWcteq66_9, &b_PDFWcteq66_9);
		}
	}
	///------------------------  MC Truth info (END) ------------------------------------

}//End Init()

void NTupleReader::ReadSelectedBranches() const {

		chain->SetBranchStatus("*", 0); //disable all branches
		chain->SetBranchStatus("run", 1);
		chain->SetBranchStatus("event", 1);
		chain->SetBranchStatus("lumiBlock", 1);
		chain->SetBranchStatus("beamSpot_x", 1); //beam spot
		chain->SetBranchStatus("beamSpot_y", 1);
		chain->SetBranchStatus("Nels", 1); //electrons
		chain->SetBranchStatus("els_px", 1);
		chain->SetBranchStatus("els_py", 1);
		chain->SetBranchStatus("els_pz", 1);
		chain->SetBranchStatus("els_energy", 1);
		chain->SetBranchStatus("els_et", 1);
		chain->SetBranchStatus("els_eta", 1);
		chain->SetBranchStatus("els_phi", 1);//z study
		chain->SetBranchStatus("els_looseId", 1);
		chain->SetBranchStatus("els_tightId", 1);
		chain->SetBranchStatus("els_robustLooseId", 1);
		chain->SetBranchStatus("els_robustTightId", 1);
		chain->SetBranchStatus("els_robustHighEnergyId", 1);
		chain->SetBranchStatus("els_sigmaIEtaIEta", 1); //sigma i eta i eta
		chain->SetBranchStatus("els_hadOverEm", 1);
		chain->SetBranchStatus("els_dEtaIn", 1);
		chain->SetBranchStatus("els_dPhiIn", 1);
		chain->SetBranchStatus("els_eOverPIn", 1);
		chain->SetBranchStatus("els_tIso", 1);
		chain->SetBranchStatus("els_dr04EcalRecHitSumEt", 1);
		chain->SetBranchStatus("els_dr04HcalTowerSumEt", 1);
		chain->SetBranchStatus("els_d0dum", 1);
		chain->SetBranchStatus("els_vx", 1);
		chain->SetBranchStatus("els_vy", 1);
		chain->SetBranchStatus("els_vpx", 1);
		chain->SetBranchStatus("els_vpy", 1);
		chain->SetBranchStatus("els_closestCtfTrackRef", 1);
		chain->SetBranchStatus("els_tk_pt", 1);
		chain->SetBranchStatus("els_tk_phi", 1);
		chain->SetBranchStatus("els_tk_eta", 1);
		chain->SetBranchStatus("els_tk_charge", 1);
		chain->SetBranchStatus("els_tk_theta", 1);
		chain->SetBranchStatus("els_shFracInnerHits", 1);
		if (m_useMisslayers)
			chain->SetBranchStatus("els_innerLayerMissingHits", 1);
		chain->SetBranchStatus("Nmus", 1); //muons
		chain->SetBranchStatus("mus_cm_px", 1); //global muon
		chain->SetBranchStatus("mus_cm_py", 1);
		chain->SetBranchStatus("mus_cm_pz", 1);
		chain->SetBranchStatus("mus_energy", 1);
		chain->SetBranchStatus("mus_cm_pt", 1);
		chain->SetBranchStatus("mus_cm_eta", 1);
		chain->SetBranchStatus("mus_cm_chi2", 1);
		chain->SetBranchStatus("mus_cm_ndof", 1);
		chain->SetBranchStatus("mus_cm_d0dum", 1);
		chain->SetBranchStatus("mus_tk_vx", 1);
		chain->SetBranchStatus("mus_tk_vy", 1);
		chain->SetBranchStatus("mus_tk_px", 1);
		chain->SetBranchStatus("mus_tk_py", 1);
		chain->SetBranchStatus("mus_tkHits", 1);
		chain->SetBranchStatus("mus_tIso", 1);
		chain->SetBranchStatus("mus_cIso", 1);
		chain->SetBranchStatus("mus_id_AllGlobalMuons", 1); //new
		chain->SetBranchStatus("mus_ecalvetoDep", 1); //new
		chain->SetBranchStatus("mus_hcalvetoDep", 1); //new

		if (m_jetAlgo == "Default") { //using default jet-met
			chain->SetBranchStatus("Njets", 1); //jets
			chain->SetBranchStatus("jets_px", 1);
			chain->SetBranchStatus("jets_py", 1);
			chain->SetBranchStatus("jets_pz", 1);
			chain->SetBranchStatus("jets_energy", 1);
			chain->SetBranchStatus("jets_eta", 1);
			chain->SetBranchStatus("jets_pt", 1);
			chain->SetBranchStatus("jets_emf", 1);

		} else if (m_jetAlgo == "pfjet") { //PFJet
			chain->SetBranchStatus(Form("NPFJets"), 1); //jets
			chain->SetBranchStatus("PFJets_px", 1);
			chain->SetBranchStatus("PFJets_py", 1);
			chain->SetBranchStatus("PFJets_pz", 1);
			chain->SetBranchStatus("PFJets_energy", 1);
			chain->SetBranchStatus("PFJets_eta", 1);
			chain->SetBranchStatus("PFJets_pt", 1);

		} else { //if not using default jets
			chain->SetBranchStatus(Form("Njets%s", m_jetAlgo.c_str()), 1); //jets
			chain->SetBranchStatus(Form("jets%s_px", m_jetAlgo.c_str()), 1);
			chain->SetBranchStatus(Form("jets%s_py", m_jetAlgo.c_str()), 1);
			chain->SetBranchStatus(Form("jets%s_pz", m_jetAlgo.c_str()), 1);
			chain->SetBranchStatus(Form("jets%s_energy", m_jetAlgo.c_str()), 1);
			chain->SetBranchStatus(Form("jets%s_eta", m_jetAlgo.c_str()), 1);
			chain->SetBranchStatus(Form("jets%s_pt", m_jetAlgo.c_str()), 1);
			chain->SetBranchStatus(Form("jets%s_emf", m_jetAlgo.c_str()), 1);
		}

		// 30-10-09
		chain->SetBranchStatus("Nmets", 1);
		chain->SetBranchStatus("mets_et_muonCor", 1);
		chain->SetBranchStatus("mets_phi_muonCor", 1);
		chain->SetBranchStatus("mets_et", 1);
		chain->SetBranchStatus("mets_phi", 1);
		chain->SetBranchStatus("mets_ex", 1);
		chain->SetBranchStatus("mets_ey", 1);

		if (m_metAlgo == "tcmet") { //tcMET
			chain->SetBranchStatus("Ntcmets", 1);
			chain->SetBranchStatus("tcmets_et", 1);
			chain->SetBranchStatus("tcmets_phi", 1);
			chain->SetBranchStatus("tcmets_ex", 1);
			chain->SetBranchStatus("tcmets_ey", 1);
		} else if (m_metAlgo == "pfmet") { //PFMET
			chain->SetBranchStatus("NPFMets", 1);
			chain->SetBranchStatus("PFMets_et", 1);
			chain->SetBranchStatus("PFMets_phi", 1);
			chain->SetBranchStatus("PFMets_ex", 1);
			chain->SetBranchStatus("PFMets_ey", 1);

		} else if (m_metAlgo == "SC5" || m_metAlgo == "SC7" || m_metAlgo == "KT4" || m_metAlgo == "KT6") { //CaloMET: SC5/7, KT4/6
			chain->SetBranchStatus(Form("Nmets%s", m_metAlgo.c_str()), 1);
			chain->SetBranchStatus(Form("mets%s_et", m_metAlgo.c_str()), 1);
			chain->SetBranchStatus(Form("mets%s_phi", m_metAlgo.c_str()), 1);
			chain->SetBranchStatus(Form("mets%s_ex", m_metAlgo.c_str()), 1);
			chain->SetBranchStatus(Form("mets%s_ey", m_metAlgo.c_str()), 1);
		}

		chain->SetBranchStatus("Ntracks", 1); //tracks
		chain->SetBranchStatus("tracks_pt", 1);
		chain->SetBranchStatus("tracks_phi", 1);
		chain->SetBranchStatus("tracks_eta", 1);
		chain->SetBranchStatus("tracks_chg", 1);
		chain->SetBranchStatus("tracks_theta", 1);
		chain->SetBranchStatus("tracks_vx", 1);
		chain->SetBranchStatus("tracks_vy", 1);
		chain->SetBranchStatus("tracks_px", 1);
		chain->SetBranchStatus("tracks_py", 1);
		if (m_useMisslayers)
			chain->SetBranchStatus("tracks_innerLayerMissingHits", 1);

		chain->SetBranchStatus("Nphotons", 1); //z study
		chain->SetBranchStatus("photons_eta", 1); //z study
		chain->SetBranchStatus("photons_et", 1); //z study
		chain->SetBranchStatus("photons_phi", 1); //z study
		chain->SetBranchStatus("photons_px", 1); //z study
		chain->SetBranchStatus("photons_py", 1); //z study
		chain->SetBranchStatus("photons_pz", 1); //z study
		chain->SetBranchStatus("photons_energy", 1); //z study
		///-------------------------- MC Truth info ------------------------------
		if (!IsData()) {
			chain->SetBranchStatus("Nmc_doc", 1); // generator particles
			chain->SetBranchStatus("mc_doc_id", 1);
			chain->SetBranchStatus("mc_doc_status", 1);
			chain->SetBranchStatus("mc_doc_mother_id", 1);
			chain->SetBranchStatus("mc_doc_grandmother_id", 1);
			//		chain->SetBranchStatus("mc_doc_mother_id", 1);
			chain->SetBranchStatus("mc_doc_mass", 1);
			chain->SetBranchStatus("mc_doc_px", 1);
			chain->SetBranchStatus("mc_doc_py", 1);
			chain->SetBranchStatus("mc_doc_pz", 1);
			chain->SetBranchStatus("mc_doc_pt", 1);
			chain->SetBranchStatus("mc_doc_eta", 1);
			chain->SetBranchStatus("mc_doc_phi", 1); //z study
			chain->SetBranchStatus("mc_doc_energy", 1); //z study

			if (m_studyPDFunc) {
				for (short i = 0; i <= 44; i++) {
					string pdf = Form("PDFWcteq66_%u", i);//unsigned %u
					if (m_debug)
						cout << " pdf  " << pdf << endl;
					chain->SetBranchStatus(pdf.c_str(), 1);
				}
			}
		}
		///------------------------  MC Truth info (END)  ---------------------------

		if (GetTrigger()) {
			//chain->SetBranchStatus("HLT_Ele15_LW_L1R",1); //trigger (8e29)
			//chain->SetBranchStatus("HLT_Ele10_SW_L1R",1); //trigger (1e30)
			chain->SetBranchStatus(HLTBit.c_str(), 1);
		}
	}//end ReadSelectedBranches
//	------------------------------------------------------------------------------------