/*
 * NTupleReader.cc
 *
 *  Created on: Mar 23, 2010
 *      Author: lkreczko
 */

#include "NTupleReader.hh"

NTupleReader::NTupleReader() :
	ntupleChain(0), hltChain(0),

	NPFJets(0), PFJets_energy(0), PFJets_et(0), PFJets_eta(0), PFJets_phi(0), PFJets_pt(0), PFJets_px(0), PFJets_py(0),
			PFJets_pz(0), PFJets_status(0), PFJets_theta(0), PFJets_chgEmE(0), PFJets_chgHadE(0), PFJets_chgMuE(0),
			PFJets_chg_Mult(0), PFJets_neutralEmE(0), PFJets_neutralHadE(0), PFJets_neutral_Mult(0), PFJets_mu_Mult(0),
			PFJets_mass(0), PFJets_btag_TC_highEff(0), PFJets_btag_TC_highPur(0), PFJets_btag_jetBProb(0),
			PFJets_btag_jetProb(0), /*PFJets_btag_softEle(0), PFJets_btag_softMuon(0), PFJets_btag_softMuonNoIP(0),*/
			PFJets_btag_secVertex(0), PFJets_parton_Id(0),

			NPFMets(0), PFMets_et(0), PFMets_phi(0), PFMets_ex(0), PFMets_ey(0), PFMets_sumEt(0),

			NbeamSpot(0), beamSpot_x(0), beamSpot_y(0), beamSpot_z(0), beamSpot_x0Error(0), beamSpot_y0Error(0),
			beamSpot_z0Error(0), beamSpot_sigmaZ(0), beamSpot_sigmaZ0Error(0), beamSpot_dxdz(0), beamSpot_dxdzError(0),
			beamSpot_dydz(0), beamSpot_dydzError(0), beamSpot_beamWidthX(0), beamSpot_beamWidthY(0),
			beamSpot_beamWidthXError(0), beamSpot_beamWidthYError(0),

			Nels(0), els_energy(0), els_et(0), els_eta(0), els_phi(0), els_pt(0), els_px(0), els_py(0), els_pz(0),
			els_status(0), els_theta(0), els_closestCtfTrackRef(0), els_isEcalDriven(0), els_isTrackerDriven(0),
			els_dr03EcalRecHitSumEt(0), els_dr04EcalRecHitSumEt(0), els_dr03HcalTowerSumEt(0),
			els_dr04HcalTowerSumEt(0), els_gen_id(0), els_gen_phi(0), els_gen_pt(0), els_gen_pz(0), els_gen_px(0),
			els_gen_py(0), els_gen_eta(0), els_gen_theta(0), els_gen_et(0), els_gen_mother_id(0),
			els_gen_mother_phi(0), els_gen_mother_pt(0), els_gen_mother_pz(0), els_gen_mother_px(0), els_gen_mother_py(
					0), els_gen_mother_eta(0), els_gen_mother_theta(0), els_gen_mother_et(0), els_tightId(0),
			els_looseId(0), els_robustTightId(0), els_robustLooseId(0), els_robustHighEnergyId(0), els_cIso(0),
			els_tIso(0), els_ecalIso(0), els_hcalIso(0), els_chi2(0), els_charge(0), els_caloEnergy(0),
			els_hadOverEm(0), els_eOverPIn(0), els_eSeedOverPOut(0), els_eSCraw(0), els_eSeed(0), els_sigmaEtaEta(0),
			els_sigmaIEtaIEta(0), els_scE1x5(0), els_scE2x5Max(0), els_scE5x5(0), els_dEtaIn(0), els_dPhiIn(0),
			els_dEtaOut(0), els_dPhiOut(0), els_numvalhits(0), els_numlosthits(0), els_basicClustersSize(0), els_tk_pt(
					0), els_tk_phi(0), els_tk_eta(0), els_tk_theta(0), els_tk_charge(0), els_shFracInnerHits(0),
			els_d0dum(0), els_dz(0), els_vx(0), els_vy(0), els_vz(0), els_ndof(0), els_ptError(0), els_d0dumError(0),
			els_dzError(0), els_etaError(0), els_phiError(0), els_cpx(0), els_cpy(0), els_cpz(0), els_vpx(0),
			els_vpy(0), els_vpz(0), els_cx(0), els_cy(0), els_cz(0), els_innerLayerMissingHits(0), els_isEE(0),
			els_isEEGap(0), els_isEB(0), els_isEBGap(0), els_isConvertedPhoton(0),

			Njets(0), jets_energy(0), jets_et(0), jets_eta(0), jets_phi(0), jets_pt(0), jets_px(0), jets_py(0),
			jets_pz(0), jets_status(0), jets_theta(0), jets_parton_Id(0), jets_parton_motherId(0), jets_parton_pt(0),
			jets_parton_phi(0), jets_parton_eta(0), jets_parton_Energy(0), jets_parton_mass(0),
			jets_parton_motherID(0), jets_gen_et(0), jets_gen_pt(0), jets_gen_eta(0), jets_gen_phi(0),
			jets_gen_mass(0), jets_gen_Energy(0), jets_gen_Id(0), jets_gen_motherID(0), jets_gen_threeCharge(0),
			jets_partonFlavour(0), jets_btag_TC_highPur(0), jets_btag_TC_highEff(0), jets_btag_jetProb(0),
			/*jets_btag_jetBProb(0), jets_btag_softEle(0), jets_btag_softMuon(0), jets_btag_softMuonNoIP(0),*/
			jets_btag_secVertex(0), jets_chgEmE(0), jets_chgHadE(0), jets_chgMuE(0), jets_chg_Mult(0), jets_neutralEmE(
					0), jets_neutralHadE(0), jets_neutral_Mult(0), jets_mu_Mult(0), jets_emf(0), jets_ehf(0), jets_n60(
					0), jets_n90(0), jets_area(0), jets_mass(0),

			NjetsKT4(0), jetsKT4_energy(0), jetsKT4_et(0), jetsKT4_eta(0), jetsKT4_phi(0), jetsKT4_pt(0),
			jetsKT4_px(0), jetsKT4_py(0), jetsKT4_pz(0), jetsKT4_status(0), jetsKT4_theta(0),
			jetsKT4_btag_TC_highPur(0), jetsKT4_btag_TC_highEff(0), jetsKT4_btag_jetProb(0), jetsKT4_btag_jetBProb(0),
			/*jetsKT4_btag_softEle(0), jetsKT4_btag_softMuon(0), jetsKT4_btag_softMuonNoIP(0),*/jetsKT4_btag_secVertex(
					0), jetsKT4_chgEmE(0), jetsKT4_chgHadE(0), jetsKT4_chgMuE(0), jetsKT4_chg_Mult(0),
			jetsKT4_neutralEmE(0), jetsKT4_neutralHadE(0), jetsKT4_neutral_Mult(0), jetsKT4_mu_Mult(0), jetsKT4_emf(0),
			jetsKT4_ehf(0), jetsKT4_n60(0), jetsKT4_n90(0), jetsKT4_area(0), jetsKT4_mass(0),

			NjetsSC5(0), jetsSC5_energy(0), jetsSC5_et(0), jetsSC5_eta(0), jetsSC5_phi(0), jetsSC5_pt(0),
			jetsSC5_px(0), jetsSC5_py(0), jetsSC5_pz(0), jetsSC5_status(0), jetsSC5_theta(0),
			jetsSC5_btag_TC_highPur(0), jetsSC5_btag_TC_highEff(0), jetsSC5_btag_jetProb(0), jetsSC5_btag_jetBProb(0),
			/*jetsSC5_btag_softEle(0), jetsSC5_btag_softMuon(0), jetsSC5_btag_softMuonNoIP(0),*/jetsSC5_btag_secVertex(
					0), jetsSC5_chgEmE(0), jetsSC5_chgHadE(0), jetsSC5_chgMuE(0), jetsSC5_chg_Mult(0),
			jetsSC5_neutralEmE(0), jetsSC5_neutralHadE(0), jetsSC5_neutral_Mult(0), jetsSC5_mu_Mult(0), jetsSC5_emf(0),
			jetsSC5_ehf(0), jetsSC5_n60(0), jetsSC5_n90(0), jetsSC5_area(0), jetsSC5_mass(0),

			Nmc_doc(0), mc_doc_id(0), mc_doc_pt(0), mc_doc_px(0), mc_doc_py(0), mc_doc_pz(0), mc_doc_eta(0),
			mc_doc_phi(0), mc_doc_theta(0), mc_doc_energy(0), mc_doc_status(0), mc_doc_charge(0), mc_doc_mother_id(0),
			mc_doc_grandmother_id(0), mc_doc_ggrandmother_id(0), mc_doc_mother_pt(0), mc_doc_vertex_x(0),
			mc_doc_vertex_y(0), mc_doc_vertex_z(0), mc_doc_mass(0), mc_doc_numOfDaughters(0), mc_doc_numOfMothers(0),

			Nmets(0), mets_et(0), mets_phi(0), mets_ex(0), mets_ey(0), mets_gen_et(0), mets_gen_phi(0), mets_sign(0),
			mets_sumEt(0), mets_unCPhi(0), mets_unCPt(0), mets_et_muonCor(0), mets_phi_muonCor(0), mets_et_JESCor(0),
			mets_phi_JESCor(0),

			NmetsKT4(0), metsKT4_et(0), metsKT4_phi(0), metsKT4_ex(0), metsKT4_ey(0), metsKT4_sumEt(0),
			metsKT4_et_JESCor(0), metsKT4_phi_JESCor(0),

			NmetsSC5(0), metsSC5_et(0), metsSC5_phi(0), metsSC5_ex(0), metsSC5_ey(0), metsSC5_sumEt(0),
			metsSC5_et_JESCor(0), metsSC5_phi_JESCor(0),

			Nmus(0), mus_energy(0), mus_et(0), mus_eta(0), mus_phi(0), mus_pt(0), mus_px(0), mus_py(0), mus_pz(0),
			mus_status(0), mus_theta(0), mus_gen_id(0), mus_gen_phi(0), mus_gen_pt(0), mus_gen_pz(0), mus_gen_px(0),
			mus_gen_py(0), mus_gen_eta(0), mus_gen_theta(0), mus_gen_et(0), mus_gen_mother_id(0),
			mus_gen_mother_phi(0), mus_gen_mother_pt(0), mus_gen_mother_pz(0), mus_gen_mother_px(0), mus_gen_mother_py(
					0), mus_gen_mother_eta(0), mus_gen_mother_theta(0), mus_gen_mother_et(0), mus_tkHits(0),
			mus_cIso(0), mus_tIso(0), mus_ecalIso(0), mus_hcalIso(0), mus_ecalvetoDep(0), mus_hcalvetoDep(0),
			mus_calEnergyEm(0), mus_calEnergyHad(0), mus_calEnergyHo(0), mus_calEnergyEmS9(0), mus_calEnergyHadS9(0),
			mus_calEnergyHoS9(0), mus_iso03_sumPt(0), mus_iso03_emEt(0), mus_iso03_hadEt(0), mus_iso03_hoEt(0),
			mus_iso03_nTracks(0), mus_iso05_sumPt(0), mus_iso05_emEt(0), mus_iso05_hadEt(0), mus_iso05_hoEt(0),
			mus_iso05_nTracks(0), mus_charge(0), mus_cm_chi2(0), mus_cm_ndof(0), mus_cm_chg(0), mus_cm_pt(0),
			mus_cm_px(0), mus_cm_py(0), mus_cm_pz(0), mus_cm_eta(0), mus_cm_phi(0), mus_cm_theta(0), mus_cm_d0dum(0),
			mus_cm_dz(0), mus_cm_vx(0), mus_cm_vy(0), mus_cm_vz(0), mus_cm_numvalhits(0), mus_cm_numlosthits(0),
			mus_cm_d0dumErr(0), mus_cm_dzErr(0), mus_cm_ptErr(0), mus_cm_etaErr(0), mus_cm_phiErr(0), mus_tk_chi2(0),
			mus_tk_ndof(0), mus_tk_chg(0), mus_tk_pt(0), mus_tk_px(0), mus_tk_py(0), mus_tk_pz(0), mus_tk_eta(0),
			mus_tk_phi(0), mus_tk_theta(0), mus_tk_d0dum(0), mus_tk_dz(0), mus_tk_vx(0), mus_tk_vy(0), mus_tk_vz(0),
			mus_tk_numvalhits(0), mus_tk_numlosthits(0), mus_tk_d0dumErr(0), mus_tk_dzErr(0), mus_tk_ptErr(0),
			mus_tk_etaErr(0), mus_tk_phiErr(0), mus_stamu_chi2(0), mus_stamu_ndof(0), mus_stamu_chg(0),
			mus_stamu_pt(0), mus_stamu_px(0), mus_stamu_py(0), mus_stamu_pz(0), mus_stamu_eta(0), mus_stamu_phi(0),
			mus_stamu_theta(0), mus_stamu_d0dum(0), mus_stamu_dz(0), mus_stamu_vx(0), mus_stamu_vy(0), mus_stamu_vz(0),
			mus_stamu_numvalhits(0), mus_stamu_numlosthits(0), mus_stamu_d0dumErr(0), mus_stamu_dzErr(0),
			mus_stamu_ptErr(0), mus_stamu_etaErr(0), mus_stamu_phiErr(0), mus_num_matches(0), mus_id_All(0),
			mus_id_AllGlobalMuons(0), mus_id_AllStandAloneMuons(0), mus_id_AllTrackerMuons(0),
			mus_id_TrackerMuonArbitrated(0), mus_id_AllArbitrated(0), mus_id_GlobalMuonPromptTight(0),
			mus_id_TMLastStationLoose(0), mus_id_TMLastStationTight(0), mus_id_TM2DCompatibilityLoose(0),
			mus_id_TM2DCompatibilityTight(0), mus_id_TMOneStationLoose(0), mus_id_TMOneStationTight(0),
			mus_id_TMLastStationOptimizedLowPtLoose(0), mus_id_TMLastStationOptimizedLowPtTight(0),

			Nphotons(0), photons_energy(0), photons_et(0), photons_eta(0), photons_phi(0), photons_pt(0),
			photons_px(0), photons_py(0), photons_pz(0), photons_status(0), photons_theta(0), photons_hadOverEM(0),
			photons_scEnergy(0), photons_scRawEnergy(0), photons_scEta(0), photons_scPhi(0), photons_scEtaWidth(0),
			photons_scPhiWidth(0), photons_tIso(0), photons_ecalIso(0), photons_hcalIso(0),
			photons_isoEcalRecHitDR04(0), photons_isoHcalRecHitDR04(0), photons_isoSolidTrkConeDR04(0),
			photons_isoHollowTrkConeDR04(0), photons_nTrkSolidConeDR04(0), photons_nTrkHollowConeDR04(0),
			photons_isoEcalRecHitDR03(0), photons_isoHcalRecHitDR03(0), photons_isoSolidTrkConeDR03(0),
			photons_isoHollowTrkConeDR03(0), photons_nTrkSolidConeDR03(0), photons_nTrkHollowConeDR03(0),
			photons_isAlsoElectron(0), photons_hasPixelSeed(0), photons_isConverted(0), photons_isEBGap(0),
			photons_isEEGap(0), photons_isEBEEGap(0), photons_isEBPho(0), photons_isEEPho(0), photons_isLoosePhoton(0),
			photons_isTightPhoton(0), photons_r9(0), photons_gen_et(0), photons_gen_eta(0), photons_gen_phi(0),
			photons_gen_id(0),

			Npv(0), pv_x(0), pv_y(0), pv_z(0), pv_xErr(0), pv_yErr(0), pv_zErr(0),

			Ntcmets(0), tcmets_et(0), tcmets_phi(0), tcmets_ex(0), tcmets_ey(0), tcmets_sumEt(0), tcmets_et_muonCor(0),
			tcmets_phi_muonCor(0),

			Ntracks(0), tracks_chi2(0), tracks_ndof(0), tracks_chg(0), tracks_pt(0), tracks_px(0), tracks_py(0),
			tracks_pz(0), tracks_eta(0), tracks_phi(0), tracks_theta(0), tracks_d0dum(0), tracks_dz(0), tracks_vx(0),
			tracks_vy(0), tracks_vz(0), tracks_numvalhits(0), tracks_numlosthits(0), tracks_d0dumErr(0),
			tracks_dzErr(0), tracks_ptErr(0), tracks_etaErr(0), tracks_phiErr(0), tracks_Nrechits(0), tracks_innerHitX(
					0), tracks_innerHitY(0), tracks_innerHitZ(0), tracks_outerHitX(0), tracks_outerHitY(0),
			tracks_outerHitZ(0), tracks_highPurity(0), tracks_innerLayerMissingHits(0),

			run_number(0), event_number(0), lumiBlock(0),

			HLT_Ele10_LW_EleId_L1R(0), HLT_Ele10_SW_L1R(0), HLT_Ele15_LW_L1R(0), HLT_Ele15_SC10_LW_L1R(0),
			HLT_Ele15_SW_L1R(0), HLT_Ele15_SiStrip_L1R(0), HLT_Ele20_LW_L1R(0),

			b_HLT_Ele10_LW_EleId_L1R(0), b_HLT_Ele10_SW_L1R(0), b_HLT_Ele15_LW_L1R(0), b_HLT_Ele15_SC10_LW_L1R(0),
			b_HLT_Ele15_SW_L1R(0), b_HLT_Ele15_SiStrip_L1R(0), b_HLT_Ele20_LW_L1R(0),

			PDFWcteq66_0(0), PDFWcteq66_1(0), PDFWcteq66_2(0), PDFWcteq66_3(0), PDFWcteq66_4(0), PDFWcteq66_5(0),
			PDFWcteq66_6(0), PDFWcteq66_7(0), PDFWcteq66_8(0), PDFWcteq66_9(0), PDFWcteq66_10(0), PDFWcteq66_11(0),
			PDFWcteq66_12(0), PDFWcteq66_13(0), PDFWcteq66_14(0), PDFWcteq66_15(0), PDFWcteq66_16(0), PDFWcteq66_17(0),
			PDFWcteq66_18(0), PDFWcteq66_19(0), PDFWcteq66_20(0), PDFWcteq66_21(0), PDFWcteq66_22(0), PDFWcteq66_23(0),
			PDFWcteq66_24(0), PDFWcteq66_25(0), PDFWcteq66_26(0), PDFWcteq66_27(0), PDFWcteq66_28(0), PDFWcteq66_29(0),
			PDFWcteq66_30(0), PDFWcteq66_31(0), PDFWcteq66_32(0), PDFWcteq66_33(0), PDFWcteq66_34(0), PDFWcteq66_35(0),
			PDFWcteq66_36(0), PDFWcteq66_37(0), PDFWcteq66_38(0), PDFWcteq66_39(0), PDFWcteq66_40(0), PDFWcteq66_41(0),
			PDFWcteq66_42(0), PDFWcteq66_43(0), PDFWcteq66_44(0),

			datafile(true), m_useMisslayers(false), checkTrig(true), m_studyPDFunc(false), debug_flag(false),
			m_jetAlgo("Default") {
}

bool NTupleReader::IsData() const {
	return datafile;
}

bool NTupleReader::GetDebug() const {
	return this->debug_flag;
}

std::string NTupleReader::GetHLTBit() const {
	return HLTBit;
}

std::string NTupleReader::GetJetAlgorithm() const {
	return this->m_jetAlgo;
}

std::string NTupleReader::GetMETAlgorithm() const {
	return this->m_metAlgo;
}

bool NTupleReader::GetMissLayersFlag() const {
	return m_useMisslayers;
}

bool NTupleReader::GetStudyPDFunc() const {
	return m_studyPDFunc;
}
bool NTupleReader::UseTrigger() const {
	return checkTrig;
}

void NTupleReader::SetData(bool f) {
	this->datafile = f;
}
void NTupleReader::SetDebug(bool val) {
	this->debug_flag = val;
}

void NTupleReader::SetHLTBit(std::string bit) {
	HLTBit = bit;
}

void NTupleReader::SetJetAlgo(std::string val) {
	m_jetAlgo = val;
}
void NTupleReader::SetMETAlgo(std::string val) {
	m_metAlgo = val;
}
void NTupleReader::SetStudyPDFunc(bool f) {
	m_studyPDFunc = f;
}

void NTupleReader::SetTrigger(bool val) {
	checkTrig = val;
}

void NTupleReader::SetTrigger(bool val, const std::string hlt) {
	checkTrig = val;
	HLTBit = hlt;
}

//switch to preclue missing layers as these are not in 314 data samples
void NTupleReader::UseMissLayers(bool val) {
	m_useMisslayers = val;
}
void NTupleReader::Init() {

	std::cout << "\n Initializing tree branches\n" << std::endl;

	// Set branch addresses and branch pointers
	if (m_jetAlgo == "pfjet") {
		LoadPFJets();
	}
	if (m_metAlgo == "pfmet") {
		ntupleChain->SetBranchAddress("NPFMets", &NPFMets);
		ntupleChain->SetBranchAddress("PFMets_et", &PFMets_et);
		ntupleChain->SetBranchAddress("PFMets_phi", &PFMets_phi);
		ntupleChain->SetBranchAddress("PFMets_ex", &PFMets_ex);
		ntupleChain->SetBranchAddress("PFMets_ey", &PFMets_ey);
		ntupleChain->SetBranchAddress("PFMets_sumEt", &PFMets_sumEt);
	}
	ntupleChain->SetBranchAddress("NbeamSpot", &NbeamSpot);
	ntupleChain->SetBranchAddress("beamSpot_x", &beamSpot_x);
	ntupleChain->SetBranchAddress("beamSpot_y", &beamSpot_y);
	ntupleChain->SetBranchAddress("beamSpot_z", &beamSpot_z);
	ntupleChain->SetBranchAddress("beamSpot_x0Error", &beamSpot_x0Error);
	ntupleChain->SetBranchAddress("beamSpot_y0Error", &beamSpot_y0Error);
	ntupleChain->SetBranchAddress("beamSpot_z0Error", &beamSpot_z0Error);
	ntupleChain->SetBranchAddress("beamSpot_sigmaZ", &beamSpot_sigmaZ);
	ntupleChain->SetBranchAddress("beamSpot_sigmaZ0Error", &beamSpot_sigmaZ0Error);
	ntupleChain->SetBranchAddress("beamSpot_dxdz", &beamSpot_dxdz);
	ntupleChain->SetBranchAddress("beamSpot_dxdzError", &beamSpot_dxdzError);
	ntupleChain->SetBranchAddress("beamSpot_dydz", &beamSpot_dydz);
	ntupleChain->SetBranchAddress("beamSpot_dydzError", &beamSpot_dydzError);
	ntupleChain->SetBranchAddress("beamSpot_beamWidthX", &beamSpot_beamWidthX);
	ntupleChain->SetBranchAddress("beamSpot_beamWidthY", &beamSpot_beamWidthY);
	ntupleChain->SetBranchAddress("beamSpot_beamWidthXError", &beamSpot_beamWidthXError);
	ntupleChain->SetBranchAddress("beamSpot_beamWidthYError", &beamSpot_beamWidthYError);
	ntupleChain->SetBranchAddress("Nels", &Nels);
	ntupleChain->SetBranchAddress("els_energy", &els_energy);
	ntupleChain->SetBranchAddress("els_et", &els_et);
	ntupleChain->SetBranchAddress("els_eta", &els_eta);
	ntupleChain->SetBranchAddress("els_phi", &els_phi);
	ntupleChain->SetBranchAddress("els_pt", &els_pt);
	ntupleChain->SetBranchAddress("els_px", &els_px);
	ntupleChain->SetBranchAddress("els_py", &els_py);
	ntupleChain->SetBranchAddress("els_pz", &els_pz);
	ntupleChain->SetBranchAddress("els_status", &els_status);
	ntupleChain->SetBranchAddress("els_theta", &els_theta);
	ntupleChain->SetBranchAddress("els_closestCtfTrackRef", &els_closestCtfTrackRef);
	ntupleChain->SetBranchAddress("els_isEcalDriven", &els_isEcalDriven);
	ntupleChain->SetBranchAddress("els_isTrackerDriven", &els_isTrackerDriven);
	ntupleChain->SetBranchAddress("els_dr03EcalRecHitSumEt", &els_dr03EcalRecHitSumEt);
	ntupleChain->SetBranchAddress("els_dr04EcalRecHitSumEt", &els_dr04EcalRecHitSumEt);
	ntupleChain->SetBranchAddress("els_dr03HcalTowerSumEt", &els_dr03HcalTowerSumEt);
	ntupleChain->SetBranchAddress("els_dr04HcalTowerSumEt", &els_dr04HcalTowerSumEt);
	ntupleChain->SetBranchAddress("els_gen_id", &els_gen_id);
	ntupleChain->SetBranchAddress("els_gen_phi", &els_gen_phi);
	ntupleChain->SetBranchAddress("els_gen_pt", &els_gen_pt);
	ntupleChain->SetBranchAddress("els_gen_pz", &els_gen_pz);
	ntupleChain->SetBranchAddress("els_gen_px", &els_gen_px);
	ntupleChain->SetBranchAddress("els_gen_py", &els_gen_py);
	ntupleChain->SetBranchAddress("els_gen_eta", &els_gen_eta);
	ntupleChain->SetBranchAddress("els_gen_theta", &els_gen_theta);
	ntupleChain->SetBranchAddress("els_gen_et", &els_gen_et);
	ntupleChain->SetBranchAddress("els_gen_mother_id", &els_gen_mother_id);
	ntupleChain->SetBranchAddress("els_gen_mother_phi", &els_gen_mother_phi);
	ntupleChain->SetBranchAddress("els_gen_mother_pt", &els_gen_mother_pt);
	ntupleChain->SetBranchAddress("els_gen_mother_pz", &els_gen_mother_pz);
	ntupleChain->SetBranchAddress("els_gen_mother_px", &els_gen_mother_px);
	ntupleChain->SetBranchAddress("els_gen_mother_py", &els_gen_mother_py);
	ntupleChain->SetBranchAddress("els_gen_mother_eta", &els_gen_mother_eta);
	ntupleChain->SetBranchAddress("els_gen_mother_theta", &els_gen_mother_theta);
	ntupleChain->SetBranchAddress("els_gen_mother_et", &els_gen_mother_et);
	ntupleChain->SetBranchAddress("els_tightId", &els_tightId);
	ntupleChain->SetBranchAddress("els_looseId", &els_looseId);
	ntupleChain->SetBranchAddress("els_robustTightId", &els_robustTightId);
	ntupleChain->SetBranchAddress("els_robustLooseId", &els_robustLooseId);
	ntupleChain->SetBranchAddress("els_robustHighEnergyId", &els_robustHighEnergyId);
	ntupleChain->SetBranchAddress("els_cIso", &els_cIso);
	ntupleChain->SetBranchAddress("els_tIso", &els_tIso);
	ntupleChain->SetBranchAddress("els_ecalIso", &els_ecalIso);
	ntupleChain->SetBranchAddress("els_hcalIso", &els_hcalIso);
	ntupleChain->SetBranchAddress("els_chi2", &els_chi2);
	ntupleChain->SetBranchAddress("els_charge", &els_charge);
	ntupleChain->SetBranchAddress("els_caloEnergy", &els_caloEnergy);
	ntupleChain->SetBranchAddress("els_hadOverEm", &els_hadOverEm);
	ntupleChain->SetBranchAddress("els_eOverPIn", &els_eOverPIn);
	ntupleChain->SetBranchAddress("els_eSeedOverPOut", &els_eSeedOverPOut);
	ntupleChain->SetBranchAddress("els_eSCraw", &els_eSCraw);
	ntupleChain->SetBranchAddress("els_eSeed", &els_eSeed);
	ntupleChain->SetBranchAddress("els_sigmaEtaEta", &els_sigmaEtaEta);
	ntupleChain->SetBranchAddress("els_sigmaIEtaIEta", &els_sigmaIEtaIEta);
	ntupleChain->SetBranchAddress("els_scE1x5", &els_scE1x5);
	ntupleChain->SetBranchAddress("els_scE2x5Max", &els_scE2x5Max);
	ntupleChain->SetBranchAddress("els_scE5x5", &els_scE5x5);
	ntupleChain->SetBranchAddress("els_dEtaIn", &els_dEtaIn);
	ntupleChain->SetBranchAddress("els_dPhiIn", &els_dPhiIn);
	ntupleChain->SetBranchAddress("els_dEtaOut", &els_dEtaOut);
	ntupleChain->SetBranchAddress("els_dPhiOut", &els_dPhiOut);
	ntupleChain->SetBranchAddress("els_numvalhits", &els_numvalhits);
	ntupleChain->SetBranchAddress("els_numlosthits", &els_numlosthits);
	ntupleChain->SetBranchAddress("els_basicClustersSize", &els_basicClustersSize);
	ntupleChain->SetBranchAddress("els_tk_pt", &els_tk_pt);
	ntupleChain->SetBranchAddress("els_tk_phi", &els_tk_phi);
	ntupleChain->SetBranchAddress("els_tk_eta", &els_tk_eta);
	ntupleChain->SetBranchAddress("els_tk_charge", &els_tk_charge);
	ntupleChain->SetBranchAddress("els_tk_theta", &els_tk_theta);
	ntupleChain->SetBranchAddress("els_shFracInnerHits", &els_shFracInnerHits);
	ntupleChain->SetBranchAddress("els_d0dum", &els_d0dum);
	ntupleChain->SetBranchAddress("els_dz", &els_dz);
	ntupleChain->SetBranchAddress("els_vx", &els_vx);
	ntupleChain->SetBranchAddress("els_vy", &els_vy);
	ntupleChain->SetBranchAddress("els_vz", &els_vz);
	ntupleChain->SetBranchAddress("els_ndof", &els_ndof);
	ntupleChain->SetBranchAddress("els_ptError", &els_ptError);
	ntupleChain->SetBranchAddress("els_d0dumError", &els_d0dumError);
	ntupleChain->SetBranchAddress("els_dzError", &els_dzError);
	ntupleChain->SetBranchAddress("els_etaError", &els_etaError);
	ntupleChain->SetBranchAddress("els_phiError", &els_phiError);
	ntupleChain->SetBranchAddress("els_cpx", &els_cpx);
	ntupleChain->SetBranchAddress("els_cpy", &els_cpy);
	ntupleChain->SetBranchAddress("els_cpz", &els_cpz);
	ntupleChain->SetBranchAddress("els_vpx", &els_vpx);
	ntupleChain->SetBranchAddress("els_vpy", &els_vpy);
	ntupleChain->SetBranchAddress("els_vpz", &els_vpz);
	ntupleChain->SetBranchAddress("els_cx", &els_cx);
	ntupleChain->SetBranchAddress("els_cy", &els_cy);
	ntupleChain->SetBranchAddress("els_cz", &els_cz);
	ntupleChain->SetBranchAddress("els_isEE", &els_isEE);
	ntupleChain->SetBranchAddress("els_isEEGap", &els_isEEGap);
	ntupleChain->SetBranchAddress("els_isEB", &els_isEB);
	ntupleChain->SetBranchAddress("els_isEBGap", &els_isEBGap);
	ntupleChain->SetBranchAddress("els_isConvertedPhoton", &els_isConvertedPhoton);
	if (m_useMisslayers)
		ntupleChain->SetBranchAddress("els_innerLayerMissingHits", &els_innerLayerMissingHits);

	ntupleChain->SetBranchAddress("Njets", &Njets);
	ntupleChain->SetBranchAddress("jets_energy", &jets_energy);
	ntupleChain->SetBranchAddress("jets_et", &jets_et);
	ntupleChain->SetBranchAddress("jets_eta", &jets_eta);
	ntupleChain->SetBranchAddress("jets_phi", &jets_phi);
	ntupleChain->SetBranchAddress("jets_pt", &jets_pt);
	ntupleChain->SetBranchAddress("jets_px", &jets_px);
	ntupleChain->SetBranchAddress("jets_py", &jets_py);
	ntupleChain->SetBranchAddress("jets_pz", &jets_pz);
	ntupleChain->SetBranchAddress("jets_status", &jets_status);
	ntupleChain->SetBranchAddress("jets_theta", &jets_theta);
	ntupleChain->SetBranchAddress("jets_parton_Id", &jets_parton_Id);
	ntupleChain->SetBranchAddress("jets_parton_motherId", &jets_parton_motherId);
	ntupleChain->SetBranchAddress("jets_parton_pt", &jets_parton_pt);
	ntupleChain->SetBranchAddress("jets_parton_phi", &jets_parton_phi);
	ntupleChain->SetBranchAddress("jets_parton_eta", &jets_parton_eta);
	ntupleChain->SetBranchAddress("jets_parton_Energy", &jets_parton_Energy);
	ntupleChain->SetBranchAddress("jets_parton_mass", &jets_parton_mass);
	ntupleChain->SetBranchAddress("jets_parton_motherID", &jets_parton_motherID);
	ntupleChain->SetBranchAddress("jets_gen_et", &jets_gen_et);
	ntupleChain->SetBranchAddress("jets_gen_pt", &jets_gen_pt);
	ntupleChain->SetBranchAddress("jets_gen_eta", &jets_gen_eta);
	ntupleChain->SetBranchAddress("jets_gen_phi", &jets_gen_phi);
	ntupleChain->SetBranchAddress("jets_gen_mass", &jets_gen_mass);
	ntupleChain->SetBranchAddress("jets_gen_Energy", &jets_gen_Energy);
	ntupleChain->SetBranchAddress("jets_gen_Id", &jets_gen_Id);
	ntupleChain->SetBranchAddress("jets_gen_motherID", &jets_gen_motherID);
	ntupleChain->SetBranchAddress("jets_gen_threeCharge", &jets_gen_threeCharge);
	ntupleChain->SetBranchAddress("jets_partonFlavour", &jets_partonFlavour);
	ntupleChain->SetBranchAddress("jets_btag_TC_highPur", &jets_btag_TC_highPur);
	ntupleChain->SetBranchAddress("jets_btag_TC_highEff", &jets_btag_TC_highEff);
	ntupleChain->SetBranchAddress("jets_btag_jetProb", &jets_btag_jetProb);
	ntupleChain->SetBranchAddress("jets_btag_jetBProb", &jets_btag_jetBProb);
	//	chain->SetBranchAddress("jets_btag_softEle", &jets_btag_softEle);
	//	chain->SetBranchAddress("jets_btag_softMuon", &jets_btag_softMuon);
	//	chain->SetBranchAddress("jets_btag_softMuonNoIP", &jets_btag_softMuonNoIP);
	ntupleChain->SetBranchAddress("jets_btag_secVertex", &jets_btag_secVertex);
	ntupleChain->SetBranchAddress("jets_chgEmE", &jets_chgEmE);
	ntupleChain->SetBranchAddress("jets_chgHadE", &jets_chgHadE);
	ntupleChain->SetBranchAddress("jets_chgMuE", &jets_chgMuE);
	ntupleChain->SetBranchAddress("jets_chg_Mult", &jets_chg_Mult);
	ntupleChain->SetBranchAddress("jets_neutralEmE", &jets_neutralEmE);
	ntupleChain->SetBranchAddress("jets_neutralHadE", &jets_neutralHadE);
	ntupleChain->SetBranchAddress("jets_neutral_Mult", &jets_neutral_Mult);
	ntupleChain->SetBranchAddress("jets_mu_Mult", &jets_mu_Mult);
	ntupleChain->SetBranchAddress("jets_emf", &jets_emf);
	ntupleChain->SetBranchAddress("jets_ehf", &jets_ehf);
	ntupleChain->SetBranchAddress("jets_n60", &jets_n60);
	ntupleChain->SetBranchAddress("jets_n90", &jets_n90);
	ntupleChain->SetBranchAddress("jets_area", &jets_area);
	ntupleChain->SetBranchAddress("jets_mass", &jets_mass);

	if (m_jetAlgo == "KT4") {
		ntupleChain->SetBranchAddress("jetsKT4_energy", &jetsKT4_energy);
		ntupleChain->SetBranchAddress("jetsKT4_et", &jetsKT4_et);
		ntupleChain->SetBranchAddress("jetsKT4_eta", &jetsKT4_eta);
		ntupleChain->SetBranchAddress("jetsKT4_phi", &jetsKT4_phi);
		ntupleChain->SetBranchAddress("jetsKT4_pt", &jetsKT4_pt);
		ntupleChain->SetBranchAddress("jetsKT4_px", &jetsKT4_px);
		ntupleChain->SetBranchAddress("jetsKT4_py", &jetsKT4_py);
		ntupleChain->SetBranchAddress("jetsKT4_pz", &jetsKT4_pz);
		ntupleChain->SetBranchAddress("jetsKT4_status", &jetsKT4_status);
		ntupleChain->SetBranchAddress("jetsKT4_theta", &jetsKT4_theta);
		ntupleChain->SetBranchAddress("jetsKT4_btag_TC_highPur", &jetsKT4_btag_TC_highPur);
		ntupleChain->SetBranchAddress("jetsKT4_btag_TC_highEff", &jetsKT4_btag_TC_highEff);
		ntupleChain->SetBranchAddress("jetsKT4_btag_jetProb", &jetsKT4_btag_jetProb);
		ntupleChain->SetBranchAddress("jetsKT4_btag_jetBProb", &jetsKT4_btag_jetBProb);
		//		chain->SetBranchAddress("jetsKT4_btag_softEle", &jetsKT4_btag_softEle);
		//		chain->SetBranchAddress("jetsKT4_btag_softMuon", &jetsKT4_btag_softMuon);
		//		chain->SetBranchAddress("jetsKT4_btag_softMuonNoIP", &jetsKT4_btag_softMuonNoIP);
		ntupleChain->SetBranchAddress("jetsKT4_btag_secVertex", &jetsKT4_btag_secVertex);
		ntupleChain->SetBranchAddress("jetsKT4_chgEmE", &jetsKT4_chgEmE);
		ntupleChain->SetBranchAddress("jetsKT4_chgHadE", &jetsKT4_chgHadE);
		ntupleChain->SetBranchAddress("jetsKT4_chgMuE", &jetsKT4_chgMuE);
		ntupleChain->SetBranchAddress("jetsKT4_chg_Mult", &jetsKT4_chg_Mult);
		ntupleChain->SetBranchAddress("jetsKT4_neutralEmE", &jetsKT4_neutralEmE);
		ntupleChain->SetBranchAddress("jetsKT4_neutralHadE", &jetsKT4_neutralHadE);
		ntupleChain->SetBranchAddress("jetsKT4_neutral_Mult", &jetsKT4_neutral_Mult);
		ntupleChain->SetBranchAddress("jetsKT4_mu_Mult", &jetsKT4_mu_Mult);
		ntupleChain->SetBranchAddress("jetsKT4_emf", &jetsKT4_emf);
		ntupleChain->SetBranchAddress("jetsKT4_ehf", &jetsKT4_ehf);
		ntupleChain->SetBranchAddress("jetsKT4_n60", &jetsKT4_n60);
		ntupleChain->SetBranchAddress("jetsKT4_n90", &jetsKT4_n90);
		ntupleChain->SetBranchAddress("jetsKT4_area", &jetsKT4_area);
		ntupleChain->SetBranchAddress("jetsKT4_mass", &jetsKT4_mass);
	}

	if (m_jetAlgo == "SC5") {
		ntupleChain->SetBranchAddress("NjetsSC5", &NjetsSC5);
		ntupleChain->SetBranchAddress("jetsSC5_energy", &jetsSC5_energy);
		ntupleChain->SetBranchAddress("jetsSC5_et", &jetsSC5_et);
		ntupleChain->SetBranchAddress("jetsSC5_eta", &jetsSC5_eta);
		ntupleChain->SetBranchAddress("jetsSC5_phi", &jetsSC5_phi);
		ntupleChain->SetBranchAddress("jetsSC5_pt", &jetsSC5_pt);
		ntupleChain->SetBranchAddress("jetsSC5_px", &jetsSC5_px);
		ntupleChain->SetBranchAddress("jetsSC5_py", &jetsSC5_py);
		ntupleChain->SetBranchAddress("jetsSC5_pz", &jetsSC5_pz);
		ntupleChain->SetBranchAddress("jetsSC5_status", &jetsSC5_status);
		ntupleChain->SetBranchAddress("jetsSC5_theta", &jetsSC5_theta);
		ntupleChain->SetBranchAddress("jetsSC5_btag_TC_highPur", &jetsSC5_btag_TC_highPur);
		ntupleChain->SetBranchAddress("jetsSC5_btag_TC_highEff", &jetsSC5_btag_TC_highEff);
		ntupleChain->SetBranchAddress("jetsSC5_btag_jetProb", &jetsSC5_btag_jetProb);
		ntupleChain->SetBranchAddress("jetsSC5_btag_jetBProb", &jetsSC5_btag_jetBProb);
		//		chain->SetBranchAddress("jetsSC5_btag_softEle", &jetsSC5_btag_softEle);
		//		chain->SetBranchAddress("jetsSC5_btag_softMuon", &jetsSC5_btag_softMuon);
		//		chain->SetBranchAddress("jetsSC5_btag_softMuonNoIP", &jetsSC5_btag_softMuonNoIP);
		ntupleChain->SetBranchAddress("jetsSC5_btag_secVertex", &jetsSC5_btag_secVertex);
		ntupleChain->SetBranchAddress("jetsSC5_chgEmE", &jetsSC5_chgEmE);
		ntupleChain->SetBranchAddress("jetsSC5_chgHadE", &jetsSC5_chgHadE);
		ntupleChain->SetBranchAddress("jetsSC5_chgMuE", &jetsSC5_chgMuE);
		ntupleChain->SetBranchAddress("jetsSC5_chg_Mult", &jetsSC5_chg_Mult);
		ntupleChain->SetBranchAddress("jetsSC5_neutralEmE", &jetsSC5_neutralEmE);
		ntupleChain->SetBranchAddress("jetsSC5_neutralHadE", &jetsSC5_neutralHadE);
		ntupleChain->SetBranchAddress("jetsSC5_neutral_Mult", &jetsSC5_neutral_Mult);
		ntupleChain->SetBranchAddress("jetsSC5_mu_Mult", &jetsSC5_mu_Mult);
		ntupleChain->SetBranchAddress("jetsSC5_emf", &jetsSC5_emf);
		ntupleChain->SetBranchAddress("jetsSC5_ehf", &jetsSC5_ehf);
		ntupleChain->SetBranchAddress("jetsSC5_n60", &jetsSC5_n60);
		ntupleChain->SetBranchAddress("jetsSC5_n90", &jetsSC5_n90);
		ntupleChain->SetBranchAddress("jetsSC5_area", &jetsSC5_area);
		ntupleChain->SetBranchAddress("jetsSC5_mass", &jetsSC5_mass);
	}
	ntupleChain->SetBranchAddress("Nmets", &Nmets);
	ntupleChain->SetBranchAddress("mets_et", &mets_et);
	ntupleChain->SetBranchAddress("mets_phi", &mets_phi);
	ntupleChain->SetBranchAddress("mets_ex", &mets_ex);
	ntupleChain->SetBranchAddress("mets_ey", &mets_ey);
	ntupleChain->SetBranchAddress("mets_gen_et", &mets_gen_et);
	ntupleChain->SetBranchAddress("mets_gen_phi", &mets_gen_phi);
	ntupleChain->SetBranchAddress("mets_sign", &mets_sign);
	ntupleChain->SetBranchAddress("mets_sumEt", &mets_sumEt);
	ntupleChain->SetBranchAddress("mets_unCPhi", &mets_unCPhi);
	ntupleChain->SetBranchAddress("mets_unCPt", &mets_unCPt);
	ntupleChain->SetBranchAddress("mets_et_muonCor", &mets_et_muonCor);
	ntupleChain->SetBranchAddress("mets_phi_muonCor", &mets_phi_muonCor);
	ntupleChain->SetBranchAddress("mets_et_JESCor", &mets_et_JESCor);
	ntupleChain->SetBranchAddress("mets_phi_JESCor", &mets_phi_JESCor);

	if (m_metAlgo == "KT4") {
		ntupleChain->SetBranchAddress("NmetsKT4", &NmetsKT4);
		ntupleChain->SetBranchAddress("metsKT4_et", &metsKT4_et);
		ntupleChain->SetBranchAddress("metsKT4_phi", &metsKT4_phi);
		ntupleChain->SetBranchAddress("metsKT4_ex", &metsKT4_ex);
		ntupleChain->SetBranchAddress("metsKT4_ey", &metsKT4_ey);
		ntupleChain->SetBranchAddress("metsKT4_sumEt", &metsKT4_sumEt);
		ntupleChain->SetBranchAddress("metsKT4_et_JESCor", &metsKT4_et_JESCor);
		ntupleChain->SetBranchAddress("metsKT4_phi_JESCor", &metsKT4_phi_JESCor);
	}
	if (m_metAlgo == "SC5") {
		ntupleChain->SetBranchAddress("NmetsSC5", &NmetsSC5);
		ntupleChain->SetBranchAddress("metsSC5_et", &metsSC5_et);
		ntupleChain->SetBranchAddress("metsSC5_phi", &metsSC5_phi);
		ntupleChain->SetBranchAddress("metsSC5_ex", &metsSC5_ex);
		ntupleChain->SetBranchAddress("metsSC5_ey", &metsSC5_ey);
		ntupleChain->SetBranchAddress("metsSC5_sumEt", &metsSC5_sumEt);
		ntupleChain->SetBranchAddress("metsSC5_et_JESCor", &metsSC5_et_JESCor);
		ntupleChain->SetBranchAddress("metsSC5_phi_JESCor", &metsSC5_phi_JESCor);
	}
	if (m_metAlgo == "tcmet") {
		ntupleChain->SetBranchAddress("Ntcmets", &Ntcmets);
		ntupleChain->SetBranchAddress("tcmets_et", &tcmets_et);
		ntupleChain->SetBranchAddress("tcmets_phi", &tcmets_phi);
		ntupleChain->SetBranchAddress("tcmets_ex", &tcmets_ex);
		ntupleChain->SetBranchAddress("tcmets_ey", &tcmets_ey);
		ntupleChain->SetBranchAddress("tcmets_sumEt", &tcmets_sumEt);
		ntupleChain->SetBranchAddress("tcmets_et_muonCor", &tcmets_et_muonCor);
		ntupleChain->SetBranchAddress("tcmets_phi_muonCor", &tcmets_phi_muonCor);
	}

	ntupleChain->SetBranchAddress("Nmus", &Nmus);
	ntupleChain->SetBranchAddress("mus_energy", &mus_energy);
	ntupleChain->SetBranchAddress("mus_et", &mus_et);
	ntupleChain->SetBranchAddress("mus_eta", &mus_eta);
	ntupleChain->SetBranchAddress("mus_phi", &mus_phi);
	ntupleChain->SetBranchAddress("mus_pt", &mus_pt);
	ntupleChain->SetBranchAddress("mus_px", &mus_px);
	ntupleChain->SetBranchAddress("mus_py", &mus_py);
	ntupleChain->SetBranchAddress("mus_pz", &mus_pz);
	ntupleChain->SetBranchAddress("mus_status", &mus_status);
	ntupleChain->SetBranchAddress("mus_theta", &mus_theta);
	ntupleChain->SetBranchAddress("mus_gen_id", &mus_gen_id);
	ntupleChain->SetBranchAddress("mus_gen_phi", &mus_gen_phi);
	ntupleChain->SetBranchAddress("mus_gen_pt", &mus_gen_pt);
	ntupleChain->SetBranchAddress("mus_gen_pz", &mus_gen_pz);
	ntupleChain->SetBranchAddress("mus_gen_px", &mus_gen_px);
	ntupleChain->SetBranchAddress("mus_gen_py", &mus_gen_py);
	ntupleChain->SetBranchAddress("mus_gen_eta", &mus_gen_eta);
	ntupleChain->SetBranchAddress("mus_gen_theta", &mus_gen_theta);
	ntupleChain->SetBranchAddress("mus_gen_et", &mus_gen_et);
	ntupleChain->SetBranchAddress("mus_gen_mother_id", &mus_gen_mother_id);
	ntupleChain->SetBranchAddress("mus_gen_mother_phi", &mus_gen_mother_phi);
	ntupleChain->SetBranchAddress("mus_gen_mother_pt", &mus_gen_mother_pt);
	ntupleChain->SetBranchAddress("mus_gen_mother_pz", &mus_gen_mother_pz);
	ntupleChain->SetBranchAddress("mus_gen_mother_px", &mus_gen_mother_px);
	ntupleChain->SetBranchAddress("mus_gen_mother_py", &mus_gen_mother_py);
	ntupleChain->SetBranchAddress("mus_gen_mother_eta", &mus_gen_mother_eta);
	ntupleChain->SetBranchAddress("mus_gen_mother_theta", &mus_gen_mother_theta);
	ntupleChain->SetBranchAddress("mus_gen_mother_et", &mus_gen_mother_et);
	ntupleChain->SetBranchAddress("mus_tkHits", &mus_tkHits);
	ntupleChain->SetBranchAddress("mus_cIso", &mus_cIso);
	ntupleChain->SetBranchAddress("mus_tIso", &mus_tIso);
	ntupleChain->SetBranchAddress("mus_ecalIso", &mus_ecalIso);
	ntupleChain->SetBranchAddress("mus_hcalIso", &mus_hcalIso);
	ntupleChain->SetBranchAddress("mus_ecalvetoDep", &mus_ecalvetoDep);
	ntupleChain->SetBranchAddress("mus_hcalvetoDep", &mus_hcalvetoDep);
	ntupleChain->SetBranchAddress("mus_calEnergyEm", &mus_calEnergyEm);
	ntupleChain->SetBranchAddress("mus_calEnergyHad", &mus_calEnergyHad);
	ntupleChain->SetBranchAddress("mus_calEnergyHo", &mus_calEnergyHo);
	ntupleChain->SetBranchAddress("mus_calEnergyEmS9", &mus_calEnergyEmS9);
	ntupleChain->SetBranchAddress("mus_calEnergyHadS9", &mus_calEnergyHadS9);
	ntupleChain->SetBranchAddress("mus_calEnergyHoS9", &mus_calEnergyHoS9);
	ntupleChain->SetBranchAddress("mus_iso03_sumPt", &mus_iso03_sumPt);
	ntupleChain->SetBranchAddress("mus_iso03_emEt", &mus_iso03_emEt);
	ntupleChain->SetBranchAddress("mus_iso03_hadEt", &mus_iso03_hadEt);
	ntupleChain->SetBranchAddress("mus_iso03_hoEt", &mus_iso03_hoEt);
	ntupleChain->SetBranchAddress("mus_iso03_nTracks", &mus_iso03_nTracks);
	ntupleChain->SetBranchAddress("mus_iso05_sumPt", &mus_iso05_sumPt);
	ntupleChain->SetBranchAddress("mus_iso05_emEt", &mus_iso05_emEt);
	ntupleChain->SetBranchAddress("mus_iso05_hadEt", &mus_iso05_hadEt);
	ntupleChain->SetBranchAddress("mus_iso05_hoEt", &mus_iso05_hoEt);
	ntupleChain->SetBranchAddress("mus_iso05_nTracks", &mus_iso05_nTracks);
	ntupleChain->SetBranchAddress("mus_charge", &mus_charge);
	ntupleChain->SetBranchAddress("mus_cm_chi2", &mus_cm_chi2);
	ntupleChain->SetBranchAddress("mus_cm_ndof", &mus_cm_ndof);
	ntupleChain->SetBranchAddress("mus_cm_chg", &mus_cm_chg);
	ntupleChain->SetBranchAddress("mus_cm_pt", &mus_cm_pt);
	ntupleChain->SetBranchAddress("mus_cm_px", &mus_cm_px);
	ntupleChain->SetBranchAddress("mus_cm_py", &mus_cm_py);
	ntupleChain->SetBranchAddress("mus_cm_pz", &mus_cm_pz);
	ntupleChain->SetBranchAddress("mus_cm_eta", &mus_cm_eta);
	ntupleChain->SetBranchAddress("mus_cm_phi", &mus_cm_phi);
	ntupleChain->SetBranchAddress("mus_cm_theta", &mus_cm_theta);
	ntupleChain->SetBranchAddress("mus_cm_d0dum", &mus_cm_d0dum);
	ntupleChain->SetBranchAddress("mus_cm_dz", &mus_cm_dz);
	ntupleChain->SetBranchAddress("mus_cm_vx", &mus_cm_vx);
	ntupleChain->SetBranchAddress("mus_cm_vy", &mus_cm_vy);
	ntupleChain->SetBranchAddress("mus_cm_vz", &mus_cm_vz);
	ntupleChain->SetBranchAddress("mus_cm_numvalhits", &mus_cm_numvalhits);
	ntupleChain->SetBranchAddress("mus_cm_numlosthits", &mus_cm_numlosthits);
	ntupleChain->SetBranchAddress("mus_cm_d0dumErr", &mus_cm_d0dumErr);
	ntupleChain->SetBranchAddress("mus_cm_dzErr", &mus_cm_dzErr);
	ntupleChain->SetBranchAddress("mus_cm_ptErr", &mus_cm_ptErr);
	ntupleChain->SetBranchAddress("mus_cm_etaErr", &mus_cm_etaErr);
	ntupleChain->SetBranchAddress("mus_cm_phiErr", &mus_cm_phiErr);
	ntupleChain->SetBranchAddress("mus_tk_chi2", &mus_tk_chi2);
	ntupleChain->SetBranchAddress("mus_tk_ndof", &mus_tk_ndof);
	ntupleChain->SetBranchAddress("mus_tk_chg", &mus_tk_chg);
	ntupleChain->SetBranchAddress("mus_tk_pt", &mus_tk_pt);
	ntupleChain->SetBranchAddress("mus_tk_px", &mus_tk_px);
	ntupleChain->SetBranchAddress("mus_tk_py", &mus_tk_py);
	ntupleChain->SetBranchAddress("mus_tk_pz", &mus_tk_pz);
	ntupleChain->SetBranchAddress("mus_tk_eta", &mus_tk_eta);
	ntupleChain->SetBranchAddress("mus_tk_phi", &mus_tk_phi);
	ntupleChain->SetBranchAddress("mus_tk_theta", &mus_tk_theta);
	ntupleChain->SetBranchAddress("mus_tk_d0dum", &mus_tk_d0dum);
	ntupleChain->SetBranchAddress("mus_tk_dz", &mus_tk_dz);
	ntupleChain->SetBranchAddress("mus_tk_vx", &mus_tk_vx);
	ntupleChain->SetBranchAddress("mus_tk_vy", &mus_tk_vy);
	ntupleChain->SetBranchAddress("mus_tk_vz", &mus_tk_vz);
	ntupleChain->SetBranchAddress("mus_tk_numvalhits", &mus_tk_numvalhits);
	ntupleChain->SetBranchAddress("mus_tk_numlosthits", &mus_tk_numlosthits);
	ntupleChain->SetBranchAddress("mus_tk_d0dumErr", &mus_tk_d0dumErr);
	ntupleChain->SetBranchAddress("mus_tk_dzErr", &mus_tk_dzErr);
	ntupleChain->SetBranchAddress("mus_tk_ptErr", &mus_tk_ptErr);
	ntupleChain->SetBranchAddress("mus_tk_etaErr", &mus_tk_etaErr);
	ntupleChain->SetBranchAddress("mus_tk_phiErr", &mus_tk_phiErr);
	ntupleChain->SetBranchAddress("mus_stamu_chi2", &mus_stamu_chi2);
	ntupleChain->SetBranchAddress("mus_stamu_ndof", &mus_stamu_ndof);
	ntupleChain->SetBranchAddress("mus_stamu_chg", &mus_stamu_chg);
	ntupleChain->SetBranchAddress("mus_stamu_pt", &mus_stamu_pt);
	ntupleChain->SetBranchAddress("mus_stamu_px", &mus_stamu_px);
	ntupleChain->SetBranchAddress("mus_stamu_py", &mus_stamu_py);
	ntupleChain->SetBranchAddress("mus_stamu_pz", &mus_stamu_pz);
	ntupleChain->SetBranchAddress("mus_stamu_eta", &mus_stamu_eta);
	ntupleChain->SetBranchAddress("mus_stamu_phi", &mus_stamu_phi);
	ntupleChain->SetBranchAddress("mus_stamu_theta", &mus_stamu_theta);
	ntupleChain->SetBranchAddress("mus_stamu_d0dum", &mus_stamu_d0dum);
	ntupleChain->SetBranchAddress("mus_stamu_dz", &mus_stamu_dz);
	ntupleChain->SetBranchAddress("mus_stamu_vx", &mus_stamu_vx);
	ntupleChain->SetBranchAddress("mus_stamu_vy", &mus_stamu_vy);
	ntupleChain->SetBranchAddress("mus_stamu_vz", &mus_stamu_vz);
	ntupleChain->SetBranchAddress("mus_stamu_numvalhits", &mus_stamu_numvalhits);
	ntupleChain->SetBranchAddress("mus_stamu_numlosthits", &mus_stamu_numlosthits);
	ntupleChain->SetBranchAddress("mus_stamu_d0dumErr", &mus_stamu_d0dumErr);
	ntupleChain->SetBranchAddress("mus_stamu_dzErr", &mus_stamu_dzErr);
	ntupleChain->SetBranchAddress("mus_stamu_ptErr", &mus_stamu_ptErr);
	ntupleChain->SetBranchAddress("mus_stamu_etaErr", &mus_stamu_etaErr);
	ntupleChain->SetBranchAddress("mus_stamu_phiErr", &mus_stamu_phiErr);
	ntupleChain->SetBranchAddress("mus_num_matches", &mus_num_matches);
	ntupleChain->SetBranchAddress("mus_id_All", &mus_id_All);
	ntupleChain->SetBranchAddress("mus_id_AllGlobalMuons", &mus_id_AllGlobalMuons);
	ntupleChain->SetBranchAddress("mus_id_AllStandAloneMuons", &mus_id_AllStandAloneMuons);
	ntupleChain->SetBranchAddress("mus_id_AllTrackerMuons", &mus_id_AllTrackerMuons);
	ntupleChain->SetBranchAddress("mus_id_TrackerMuonArbitrated", &mus_id_TrackerMuonArbitrated);
	ntupleChain->SetBranchAddress("mus_id_AllArbitrated", &mus_id_AllArbitrated);
	ntupleChain->SetBranchAddress("mus_id_GlobalMuonPromptTight", &mus_id_GlobalMuonPromptTight);
	ntupleChain->SetBranchAddress("mus_id_TMLastStationLoose", &mus_id_TMLastStationLoose);
	ntupleChain->SetBranchAddress("mus_id_TMLastStationTight", &mus_id_TMLastStationTight);
	ntupleChain->SetBranchAddress("mus_id_TM2DCompatibilityLoose", &mus_id_TM2DCompatibilityLoose);
	ntupleChain->SetBranchAddress("mus_id_TM2DCompatibilityTight", &mus_id_TM2DCompatibilityTight);
	ntupleChain->SetBranchAddress("mus_id_TMOneStationLoose", &mus_id_TMOneStationLoose);
	ntupleChain->SetBranchAddress("mus_id_TMOneStationTight", &mus_id_TMOneStationTight);
	ntupleChain->SetBranchAddress("mus_id_TMLastStationOptimizedLowPtLoose", &mus_id_TMLastStationOptimizedLowPtLoose);
	ntupleChain->SetBranchAddress("mus_id_TMLastStationOptimizedLowPtTight", &mus_id_TMLastStationOptimizedLowPtTight);
	ntupleChain->SetBranchAddress("Npv", &Npv);
	ntupleChain->SetBranchAddress("pv_x", &pv_x);
	ntupleChain->SetBranchAddress("pv_y", &pv_y);
	ntupleChain->SetBranchAddress("pv_z", &pv_z);
	ntupleChain->SetBranchAddress("pv_xErr", &pv_xErr);
	ntupleChain->SetBranchAddress("pv_yErr", &pv_yErr);
	ntupleChain->SetBranchAddress("pv_zErr", &pv_zErr);
	ntupleChain->SetBranchAddress("run", &run_number);
	ntupleChain->SetBranchAddress("event", &event_number);
	ntupleChain->SetBranchAddress("lumiblock", &lumiBlock);

	// HLT tree (nonIso ele trigger)
	//chain2->SetBranchAddress("HLT_Ele10_SW_L1R", &HLT_Ele10_SW_L1R);
	//cout << "GetTrigger()" << GetTrigger()<< endl;
	if (UseTrigger()) {
		if (HLTBit == "HLT_Ele15_LW_L1R")
			hltChain->SetBranchAddress("HLT_Ele15_LW_L1R", &HLT_Ele15_LW_L1R);
		if (HLTBit == "HLT_Ele15_SW_L1R")
			hltChain->SetBranchAddress("HLT_Ele15_SW_L1R", &HLT_Ele15_SW_L1R);
	}
	///------------------------  MC Truth info  ------------------------------------
	if (!IsData()) {
		LoadMCInformation();
	}
	///------------------------  MC Truth info (END) ------------------------------------

}//End Init()

void NTupleReader::SelectBranches() {

	ntupleChain->SetBranchStatus("*", 0); //disable all branches
	ntupleChain->SetBranchStatus("run", 1);
	ntupleChain->SetBranchStatus("event", 1);
	ntupleChain->SetBranchStatus("lumiblock", 1);
	ntupleChain->SetBranchStatus("beamSpot_x", 1); //beam spot
	ntupleChain->SetBranchStatus("beamSpot_y", 1);
	ntupleChain->SetBranchStatus("Nels", 1); //electrons
	ntupleChain->SetBranchStatus("els_px", 1);
	ntupleChain->SetBranchStatus("els_py", 1);
	ntupleChain->SetBranchStatus("els_pz", 1);
	ntupleChain->SetBranchStatus("els_energy", 1);
	ntupleChain->SetBranchStatus("els_et", 1);
	ntupleChain->SetBranchStatus("els_eta", 1);
	ntupleChain->SetBranchStatus("els_phi", 1);//z study
	ntupleChain->SetBranchStatus("els_looseId", 1);
	ntupleChain->SetBranchStatus("els_tightId", 1);
	ntupleChain->SetBranchStatus("els_robustLooseId", 1);
	ntupleChain->SetBranchStatus("els_robustTightId", 1);
	ntupleChain->SetBranchStatus("els_robustHighEnergyId", 1);
	ntupleChain->SetBranchStatus("els_sigmaIEtaIEta", 1); //sigma i eta i eta
	ntupleChain->SetBranchStatus("els_hadOverEm", 1);
	ntupleChain->SetBranchStatus("els_dEtaIn", 1);
	ntupleChain->SetBranchStatus("els_dPhiIn", 1);
	ntupleChain->SetBranchStatus("els_eOverPIn", 1);
	ntupleChain->SetBranchStatus("els_tIso", 1);
	ntupleChain->SetBranchStatus("els_dr04EcalRecHitSumEt", 1);
	ntupleChain->SetBranchStatus("els_dr04HcalTowerSumEt", 1);
	ntupleChain->SetBranchStatus("els_d0dum", 1);
	ntupleChain->SetBranchStatus("els_vx", 1);
	ntupleChain->SetBranchStatus("els_vy", 1);
	ntupleChain->SetBranchStatus("els_vpx", 1);
	ntupleChain->SetBranchStatus("els_vpy", 1);
	ntupleChain->SetBranchStatus("els_closestCtfTrackRef", 1);
	ntupleChain->SetBranchStatus("els_tk_pt", 1);
	ntupleChain->SetBranchStatus("els_tk_phi", 1);
	ntupleChain->SetBranchStatus("els_tk_eta", 1);
	ntupleChain->SetBranchStatus("els_tk_charge", 1);
	ntupleChain->SetBranchStatus("els_tk_theta", 1);
	ntupleChain->SetBranchStatus("els_shFracInnerHits", 1);
	ntupleChain->SetBranchStatus("els_innerLayerMissingHits", 1);
	ntupleChain->SetBranchStatus("Nmus", 1); //muons
	ntupleChain->SetBranchStatus("mus_cm_px", 1); //global muon
	ntupleChain->SetBranchStatus("mus_cm_py", 1);
	ntupleChain->SetBranchStatus("mus_cm_pz", 1);
	ntupleChain->SetBranchStatus("mus_energy", 1);
	ntupleChain->SetBranchStatus("mus_cm_pt", 1);
	ntupleChain->SetBranchStatus("mus_cm_eta", 1);
	ntupleChain->SetBranchStatus("mus_cm_chi2", 1);
	ntupleChain->SetBranchStatus("mus_cm_ndof", 1);
	ntupleChain->SetBranchStatus("mus_cm_d0dum", 1);
	ntupleChain->SetBranchStatus("mus_tk_vx", 1);
	ntupleChain->SetBranchStatus("mus_tk_vy", 1);
	ntupleChain->SetBranchStatus("mus_tk_px", 1);
	ntupleChain->SetBranchStatus("mus_tk_py", 1);
	ntupleChain->SetBranchStatus("mus_tkHits", 1);
	ntupleChain->SetBranchStatus("mus_tIso", 1);
	ntupleChain->SetBranchStatus("mus_cIso", 1);
	ntupleChain->SetBranchStatus("mus_id_AllGlobalMuons", 1); //new
	ntupleChain->SetBranchStatus("mus_ecalvetoDep", 1); //new
	ntupleChain->SetBranchStatus("mus_hcalvetoDep", 1); //new

	if (m_jetAlgo == "Default") { //using default jet-met
		ntupleChain->SetBranchStatus("Njets", 1); //jets
		ntupleChain->SetBranchStatus("jets_px", 1);
		ntupleChain->SetBranchStatus("jets_py", 1);
		ntupleChain->SetBranchStatus("jets_pz", 1);
		ntupleChain->SetBranchStatus("jets_energy", 1);
		ntupleChain->SetBranchStatus("jets_eta", 1);
		ntupleChain->SetBranchStatus("jets_pt", 1);
		ntupleChain->SetBranchStatus("jets_emf", 1);
		ntupleChain->SetBranchStatus("jets_emf", 1);
		//btagging
		ntupleChain->SetBranchStatus("jets_btag_TC_highEff", 1);
		ntupleChain->SetBranchStatus("jets_btag_TC_highPur", 1);
		ntupleChain->SetBranchStatus("jets_btag_jetBProb", 1);
		ntupleChain->SetBranchStatus("jets_btag_jetProb", 1);
		//		chain->SetBranchStatus("jets_btag_softEle", 1);
		//		chain->SetBranchStatus("jets_btag_softMuon", 1);
		//		chain->SetBranchStatus("jets_btag_softMuonNoIP", 1);
		ntupleChain->SetBranchStatus("jets_btag_secVertex", 1);
		ntupleChain->SetBranchStatus("jets_parton_Id", 1);

	} else if (m_jetAlgo == "pfjet") { //PFJet
		EnablePFJets();
	} else { //if not using default jets
		ntupleChain->SetBranchStatus(Form("Njets%s", m_jetAlgo.c_str()), 1); //jets
		ntupleChain->SetBranchStatus(Form("jets%s_px", m_jetAlgo.c_str()), 1);
		ntupleChain->SetBranchStatus(Form("jets%s_py", m_jetAlgo.c_str()), 1);
		ntupleChain->SetBranchStatus(Form("jets%s_pz", m_jetAlgo.c_str()), 1);
		ntupleChain->SetBranchStatus(Form("jets%s_energy", m_jetAlgo.c_str()), 1);
		ntupleChain->SetBranchStatus(Form("jets%s_eta", m_jetAlgo.c_str()), 1);
		ntupleChain->SetBranchStatus(Form("jets%s_pt", m_jetAlgo.c_str()), 1);
		ntupleChain->SetBranchStatus(Form("jets%s_emf", m_jetAlgo.c_str()), 1);
		//btagging
		ntupleChain->SetBranchStatus(Form("jets%s_btag_TC_highEff", m_jetAlgo.c_str()), 1);
		ntupleChain->SetBranchStatus(Form("jets%s_btag_TC_highPur", m_jetAlgo.c_str()), 1);
		ntupleChain->SetBranchStatus(Form("jets%s_btag_jetBProb", m_jetAlgo.c_str()), 1);
		ntupleChain->SetBranchStatus(Form("jets%s_btag_jetProb", m_jetAlgo.c_str()), 1);
		//		chain->SetBranchStatus(Form("jets%s_btag_softEle", m_jetAlgo.c_str()), 1);
		//		chain->SetBranchStatus(Form("jets%s_btag_softMuon", m_jetAlgo.c_str()), 1);
		//		chain->SetBranchStatus(Form("jets%s_btag_softMuonNoIP", m_jetAlgo.c_str()), 1);
		ntupleChain->SetBranchStatus(Form("jets%s_btag_secVertex", m_jetAlgo.c_str()), 1);
	}

	// 30-10-09
	ntupleChain->SetBranchStatus("Nmets", 1);
	ntupleChain->SetBranchStatus("mets_et_muonCor", 1);
	ntupleChain->SetBranchStatus("mets_phi_muonCor", 1);
	ntupleChain->SetBranchStatus("mets_et", 1);
	ntupleChain->SetBranchStatus("mets_phi", 1);
	ntupleChain->SetBranchStatus("mets_ex", 1);
	ntupleChain->SetBranchStatus("mets_ey", 1);

	if (m_metAlgo == "tcmet") { //tcMET
		ntupleChain->SetBranchStatus("Ntcmets", 1);
		ntupleChain->SetBranchStatus("tcmets_et", 1);
		ntupleChain->SetBranchStatus("tcmets_phi", 1);
		ntupleChain->SetBranchStatus("tcmets_ex", 1);
		ntupleChain->SetBranchStatus("tcmets_ey", 1);
	} else if (m_metAlgo == "pfmet") { //PFMET
		ntupleChain->SetBranchStatus("NPFMets", 1);
		ntupleChain->SetBranchStatus("PFMets_et", 1);
		ntupleChain->SetBranchStatus("PFMets_phi", 1);
		ntupleChain->SetBranchStatus("PFMets_ex", 1);
		ntupleChain->SetBranchStatus("PFMets_ey", 1);

	} else if (m_metAlgo == "SC5" || m_metAlgo == "SC7" || m_metAlgo == "KT4" || m_metAlgo == "KT6") { //CaloMET: SC5/7, KT4/6
		ntupleChain->SetBranchStatus(Form("Nmets%s", m_metAlgo.c_str()), 1);
		ntupleChain->SetBranchStatus(Form("mets%s_et", m_metAlgo.c_str()), 1);
		ntupleChain->SetBranchStatus(Form("mets%s_phi", m_metAlgo.c_str()), 1);
		ntupleChain->SetBranchStatus(Form("mets%s_ex", m_metAlgo.c_str()), 1);
		ntupleChain->SetBranchStatus(Form("mets%s_ey", m_metAlgo.c_str()), 1);
	}

	///-------------------------- MC Truth info ------------------------------
	if (!IsData()) {
		EnableMCInformation();
	}
	///------------------------  MC Truth info (END)  ---------------------------

	if (UseTrigger()) {
		//chain->SetBranchStatus("HLT_Ele15_LW_L1R",1); //trigger (8e29)
		//chain->SetBranchStatus("HLT_Ele10_SW_L1R",1); //trigger (1e30)
		ntupleChain->SetBranchStatus(HLTBit.c_str(), 1);
	}
}//end ReadSelectedBranches

void NTupleReader::EnablePFJets() {
	ntupleChain->SetBranchStatus(Form("NPFJets"), 1); //jets
	ntupleChain->SetBranchStatus("PFJets_px", 1);
	ntupleChain->SetBranchStatus("PFJets_py", 1);
	ntupleChain->SetBranchStatus("PFJets_pz", 1);
	ntupleChain->SetBranchStatus("PFJets_energy", 1);
	ntupleChain->SetBranchStatus("PFJets_eta", 1);
	ntupleChain->SetBranchStatus("PFJets_pt", 1);
	//btagging
	ntupleChain->SetBranchStatus("PFJets_btag_TC_highEff", 1);
	ntupleChain->SetBranchStatus("PFJets_btag_TC_highPur", 1);
	ntupleChain->SetBranchStatus("PFJets_btag_jetBProb", 1);
	ntupleChain->SetBranchStatus("PFJets_btag_jetProb", 1);
	ntupleChain->SetBranchStatus("PFJets_btag_softEle", 1);
	ntupleChain->SetBranchStatus("PFJets_btag_softMuon", 1);
	ntupleChain->SetBranchStatus("PFJets_btag_softMuonNoIP", 1);
	ntupleChain->SetBranchStatus("PFJets_btag_secVertex", 1);
	ntupleChain->SetBranchStatus("PFJets_parton_Id", 1);
}

void NTupleReader::LoadPFJets() {
	ntupleChain->SetBranchAddress("NPFJets", &NPFJets);
	ntupleChain->SetBranchAddress("PFJets_energy", &PFJets_energy);
	ntupleChain->SetBranchAddress("PFJets_et", &PFJets_et);
	ntupleChain->SetBranchAddress("PFJets_eta", &PFJets_eta);
	ntupleChain->SetBranchAddress("PFJets_phi", &PFJets_phi);
	ntupleChain->SetBranchAddress("PFJets_pt", &PFJets_pt);
	ntupleChain->SetBranchAddress("PFJets_px", &PFJets_px);
	ntupleChain->SetBranchAddress("PFJets_py", &PFJets_py);
	ntupleChain->SetBranchAddress("PFJets_pz", &PFJets_pz);
	ntupleChain->SetBranchAddress("PFJets_status", &PFJets_status);
	ntupleChain->SetBranchAddress("PFJets_theta", &PFJets_theta);
	ntupleChain->SetBranchAddress("PFJets_chgEmE", &PFJets_chgEmE);
	ntupleChain->SetBranchAddress("PFJets_chgHadE", &PFJets_chgHadE);
	ntupleChain->SetBranchAddress("PFJets_chgMuE", &PFJets_chgMuE);
	ntupleChain->SetBranchAddress("PFJets_chg_Mult", &PFJets_chg_Mult);
	ntupleChain->SetBranchAddress("PFJets_neutralEmE", &PFJets_neutralEmE);
	ntupleChain->SetBranchAddress("PFJets_neutralHadE", &PFJets_neutralHadE);
	ntupleChain->SetBranchAddress("PFJets_neutral_Mult", &PFJets_neutral_Mult);
	ntupleChain->SetBranchAddress("PFJets_mu_Mult", &PFJets_mu_Mult);
	ntupleChain->SetBranchAddress("PFJets_mass", &PFJets_mass);
	ntupleChain->SetBranchAddress("PFJets_parton_Id", &PFJets_parton_Id);
	ntupleChain->SetBranchAddress("PFJets_btag_TC_highPur", &PFJets_btag_TC_highPur);
	ntupleChain->SetBranchAddress("PFJets_btag_TC_highEff", &PFJets_btag_TC_highEff);
	ntupleChain->SetBranchAddress("PFJets_btag_jetProb", &PFJets_btag_jetProb);
	ntupleChain->SetBranchAddress("PFJets_btag_jetBProb", &PFJets_btag_jetBProb);
	//	chain->SetBranchAddress("PFJets_btag_softEle", &PFJets_btag_softEle);
	//	chain->SetBranchAddress("PFJets_btag_softMuon", &PFJets_btag_softMuon);
	//	chain->SetBranchAddress("PFJets_btag_softMuonNoIP", &PFJets_btag_softMuonNoIP);
	ntupleChain->SetBranchAddress("PFJets_btag_secVertex", &PFJets_btag_secVertex);
}

void NTupleReader::EnableMCInformation() {
	ntupleChain->SetBranchStatus("Nmc_doc", 1); // generator particles
	ntupleChain->SetBranchStatus("mc_doc_id", 1);
	ntupleChain->SetBranchStatus("mc_doc_status", 1);
	ntupleChain->SetBranchStatus("mc_doc_mother_id", 1);
	ntupleChain->SetBranchStatus("mc_doc_grandmother_id", 1);
	//		chain->SetBranchStatus("mc_doc_mother_id", 1);
	ntupleChain->SetBranchStatus("mc_doc_mass", 1);
	ntupleChain->SetBranchStatus("mc_doc_px", 1);
	ntupleChain->SetBranchStatus("mc_doc_py", 1);
	ntupleChain->SetBranchStatus("mc_doc_pz", 1);
	ntupleChain->SetBranchStatus("mc_doc_pt", 1);
	ntupleChain->SetBranchStatus("mc_doc_eta", 1);
	ntupleChain->SetBranchStatus("mc_doc_phi", 1); //z study
	ntupleChain->SetBranchStatus("mc_doc_energy", 1); //z study

	if (m_studyPDFunc) {
		for (short i = 0; i <= 44; i++) {
			std::string pdf = Form("PDFWcteq66_%u", i);//unsigned %u
			if (debug_flag)
				std::cout << " pdf  " << pdf << std::endl;
			ntupleChain->SetBranchStatus(pdf.c_str(), 1);
		}
	}
}

void NTupleReader::LoadMCInformation() {
	std::cout << " Set MC branch address" << std::endl;
	ntupleChain->SetBranchAddress("Nmc_doc", &Nmc_doc);
	ntupleChain->SetBranchAddress("mc_doc_id", &mc_doc_id);
	ntupleChain->SetBranchAddress("mc_doc_pt", &mc_doc_pt);
	ntupleChain->SetBranchAddress("mc_doc_px", &mc_doc_px);
	ntupleChain->SetBranchAddress("mc_doc_py", &mc_doc_py);
	ntupleChain->SetBranchAddress("mc_doc_pz", &mc_doc_pz);
	ntupleChain->SetBranchAddress("mc_doc_eta", &mc_doc_eta);
	ntupleChain->SetBranchAddress("mc_doc_phi", &mc_doc_phi);
	ntupleChain->SetBranchAddress("mc_doc_theta", &mc_doc_theta);
	ntupleChain->SetBranchAddress("mc_doc_energy", &mc_doc_energy);
	ntupleChain->SetBranchAddress("mc_doc_status", &mc_doc_status);
	ntupleChain->SetBranchAddress("mc_doc_charge", &mc_doc_charge);
	ntupleChain->SetBranchAddress("mc_doc_mother_id", &mc_doc_mother_id);
	ntupleChain->SetBranchAddress("mc_doc_grandmother_id", &mc_doc_grandmother_id);
	ntupleChain->SetBranchAddress("mc_doc_ggrandmother_id", &mc_doc_ggrandmother_id);
	ntupleChain->SetBranchAddress("mc_doc_mother_pt", &mc_doc_mother_pt);
	ntupleChain->SetBranchAddress("mc_doc_vertex_x", &mc_doc_vertex_x);
	ntupleChain->SetBranchAddress("mc_doc_vertex_y", &mc_doc_vertex_y);
	ntupleChain->SetBranchAddress("mc_doc_vertex_z", &mc_doc_vertex_z);
	ntupleChain->SetBranchAddress("mc_doc_mass", &mc_doc_mass);
	ntupleChain->SetBranchAddress("mc_doc_numOfDaughters", &mc_doc_numOfDaughters);
	ntupleChain->SetBranchAddress("mc_doc_numOfMothers", &mc_doc_numOfMothers);

	if (m_studyPDFunc) {
		std::cout << "Set PDF branch address" << std::endl;
		hltChain->SetBranchAddress("PDFWcteq66_0", &PDFWcteq66_0);
		hltChain->SetBranchAddress("PDFWcteq66_1", &PDFWcteq66_1);
		hltChain->SetBranchAddress("PDFWcteq66_10", &PDFWcteq66_10);
		hltChain->SetBranchAddress("PDFWcteq66_11", &PDFWcteq66_11);
		hltChain->SetBranchAddress("PDFWcteq66_12", &PDFWcteq66_12);
		hltChain->SetBranchAddress("PDFWcteq66_13", &PDFWcteq66_13);
		hltChain->SetBranchAddress("PDFWcteq66_14", &PDFWcteq66_14);
		hltChain->SetBranchAddress("PDFWcteq66_15", &PDFWcteq66_15);
		hltChain->SetBranchAddress("PDFWcteq66_16", &PDFWcteq66_16);
		hltChain->SetBranchAddress("PDFWcteq66_17", &PDFWcteq66_17);
		hltChain->SetBranchAddress("PDFWcteq66_18", &PDFWcteq66_18);
		hltChain->SetBranchAddress("PDFWcteq66_19", &PDFWcteq66_19);
		hltChain->SetBranchAddress("PDFWcteq66_2", &PDFWcteq66_2);
		hltChain->SetBranchAddress("PDFWcteq66_20", &PDFWcteq66_20);
		hltChain->SetBranchAddress("PDFWcteq66_21", &PDFWcteq66_21);
		hltChain->SetBranchAddress("PDFWcteq66_22", &PDFWcteq66_22);
		hltChain->SetBranchAddress("PDFWcteq66_23", &PDFWcteq66_23);
		hltChain->SetBranchAddress("PDFWcteq66_24", &PDFWcteq66_24);
		hltChain->SetBranchAddress("PDFWcteq66_25", &PDFWcteq66_25);
		hltChain->SetBranchAddress("PDFWcteq66_26", &PDFWcteq66_26);
		hltChain->SetBranchAddress("PDFWcteq66_27", &PDFWcteq66_27);
		hltChain->SetBranchAddress("PDFWcteq66_28", &PDFWcteq66_28);
		hltChain->SetBranchAddress("PDFWcteq66_29", &PDFWcteq66_29);
		hltChain->SetBranchAddress("PDFWcteq66_3", &PDFWcteq66_3);
		hltChain->SetBranchAddress("PDFWcteq66_30", &PDFWcteq66_30);
		hltChain->SetBranchAddress("PDFWcteq66_31", &PDFWcteq66_31);
		hltChain->SetBranchAddress("PDFWcteq66_32", &PDFWcteq66_32);
		hltChain->SetBranchAddress("PDFWcteq66_33", &PDFWcteq66_33);
		hltChain->SetBranchAddress("PDFWcteq66_34", &PDFWcteq66_34);
		hltChain->SetBranchAddress("PDFWcteq66_35", &PDFWcteq66_35);
		hltChain->SetBranchAddress("PDFWcteq66_36", &PDFWcteq66_36);
		hltChain->SetBranchAddress("PDFWcteq66_37", &PDFWcteq66_37);
		hltChain->SetBranchAddress("PDFWcteq66_38", &PDFWcteq66_38);
		hltChain->SetBranchAddress("PDFWcteq66_39", &PDFWcteq66_39);
		hltChain->SetBranchAddress("PDFWcteq66_4", &PDFWcteq66_4);
		hltChain->SetBranchAddress("PDFWcteq66_40", &PDFWcteq66_40);
		hltChain->SetBranchAddress("PDFWcteq66_41", &PDFWcteq66_41);
		hltChain->SetBranchAddress("PDFWcteq66_42", &PDFWcteq66_42);
		hltChain->SetBranchAddress("PDFWcteq66_43", &PDFWcteq66_43);
		hltChain->SetBranchAddress("PDFWcteq66_44", &PDFWcteq66_44);
		hltChain->SetBranchAddress("PDFWcteq66_5", &PDFWcteq66_5);
		hltChain->SetBranchAddress("PDFWcteq66_6", &PDFWcteq66_6);
		hltChain->SetBranchAddress("PDFWcteq66_7", &PDFWcteq66_7);
		hltChain->SetBranchAddress("PDFWcteq66_8", &PDFWcteq66_8);
		hltChain->SetBranchAddress("PDFWcteq66_9", &PDFWcteq66_9);
	}
}

void NTupleReader::EnablePhotons() {
	ntupleChain->SetBranchStatus("Nphotons", 1); //z study
	ntupleChain->SetBranchStatus("photons_eta", 1); //z study
	ntupleChain->SetBranchStatus("photons_et", 1); //z study
	ntupleChain->SetBranchStatus("photons_phi", 1); //z study
	ntupleChain->SetBranchStatus("photons_px", 1); //z study
	ntupleChain->SetBranchStatus("photons_py", 1); //z study
	ntupleChain->SetBranchStatus("photons_pz", 1); //z study
	ntupleChain->SetBranchStatus("photons_energy", 1); //z study
}

void NTupleReader::LoadPhotons() {
	ntupleChain->SetBranchAddress("Nphotons", &Nphotons);
	ntupleChain->SetBranchAddress("photons_energy", &photons_energy);
	ntupleChain->SetBranchAddress("photons_et", &photons_et);
	ntupleChain->SetBranchAddress("photons_eta", &photons_eta);
	ntupleChain->SetBranchAddress("photons_phi", &photons_phi);
	ntupleChain->SetBranchAddress("photons_pt", &photons_pt);
	ntupleChain->SetBranchAddress("photons_px", &photons_px);
	ntupleChain->SetBranchAddress("photons_py", &photons_py);
	ntupleChain->SetBranchAddress("photons_pz", &photons_pz);
	ntupleChain->SetBranchAddress("photons_status", &photons_status);
	ntupleChain->SetBranchAddress("photons_theta", &photons_theta);
	ntupleChain->SetBranchAddress("photons_hadOverEM", &photons_hadOverEM);
	ntupleChain->SetBranchAddress("photons_scEnergy", &photons_scEnergy);
	ntupleChain->SetBranchAddress("photons_scRawEnergy", &photons_scRawEnergy);
	ntupleChain->SetBranchAddress("photons_scEta", &photons_scEta);
	ntupleChain->SetBranchAddress("photons_scPhi", &photons_scPhi);
	ntupleChain->SetBranchAddress("photons_scEtaWidth", &photons_scEtaWidth);
	ntupleChain->SetBranchAddress("photons_scPhiWidth", &photons_scPhiWidth);
	ntupleChain->SetBranchAddress("photons_tIso", &photons_tIso);
	ntupleChain->SetBranchAddress("photons_ecalIso", &photons_ecalIso);
	ntupleChain->SetBranchAddress("photons_hcalIso", &photons_hcalIso);
	ntupleChain->SetBranchAddress("photons_isoEcalRecHitDR04", &photons_isoEcalRecHitDR04);
	ntupleChain->SetBranchAddress("photons_isoHcalRecHitDR04", &photons_isoHcalRecHitDR04);
	ntupleChain->SetBranchAddress("photons_isoSolidTrkConeDR04", &photons_isoSolidTrkConeDR04);
	ntupleChain->SetBranchAddress("photons_isoHollowTrkConeDR04", &photons_isoHollowTrkConeDR04);
	ntupleChain->SetBranchAddress("photons_nTrkSolidConeDR04", &photons_nTrkSolidConeDR04);
	ntupleChain->SetBranchAddress("photons_nTrkHollowConeDR04", &photons_nTrkHollowConeDR04);
	ntupleChain->SetBranchAddress("photons_isoEcalRecHitDR03", &photons_isoEcalRecHitDR03);
	ntupleChain->SetBranchAddress("photons_isoHcalRecHitDR03", &photons_isoHcalRecHitDR03);
	ntupleChain->SetBranchAddress("photons_isoSolidTrkConeDR03", &photons_isoSolidTrkConeDR03);
	ntupleChain->SetBranchAddress("photons_isoHollowTrkConeDR03", &photons_isoHollowTrkConeDR03);
	ntupleChain->SetBranchAddress("photons_nTrkSolidConeDR03", &photons_nTrkSolidConeDR03);
	ntupleChain->SetBranchAddress("photons_nTrkHollowConeDR03", &photons_nTrkHollowConeDR03);
	ntupleChain->SetBranchAddress("photons_isAlsoElectron", &photons_isAlsoElectron);
	ntupleChain->SetBranchAddress("photons_hasPixelSeed", &photons_hasPixelSeed);
	ntupleChain->SetBranchAddress("photons_isConverted", &photons_isConverted);
	ntupleChain->SetBranchAddress("photons_isEBGap", &photons_isEBGap);
	ntupleChain->SetBranchAddress("photons_isEEGap", &photons_isEEGap);
	ntupleChain->SetBranchAddress("photons_isEBEEGap", &photons_isEBEEGap);
	ntupleChain->SetBranchAddress("photons_isEBPho", &photons_isEBPho);
	ntupleChain->SetBranchAddress("photons_isEEPho", &photons_isEEPho);
	ntupleChain->SetBranchAddress("photons_isLoosePhoton", &photons_isLoosePhoton);
	ntupleChain->SetBranchAddress("photons_isTightPhoton", &photons_isTightPhoton);
	ntupleChain->SetBranchAddress("photons_r9", &photons_r9);
	ntupleChain->SetBranchAddress("photons_gen_et", &photons_gen_et);
	ntupleChain->SetBranchAddress("photons_gen_eta", &photons_gen_eta);
	ntupleChain->SetBranchAddress("photons_gen_phi", &photons_gen_phi);
	ntupleChain->SetBranchAddress("photons_gen_id", &photons_gen_id);
}

void NTupleReader::EnableTracks() {
	ntupleChain->SetBranchStatus("Ntracks", 1); //tracks
	ntupleChain->SetBranchStatus("tracks_pt", 1);
	ntupleChain->SetBranchStatus("tracks_phi", 1);
	ntupleChain->SetBranchStatus("tracks_eta", 1);
	ntupleChain->SetBranchStatus("tracks_chg", 1);
	ntupleChain->SetBranchStatus("tracks_theta", 1);
	ntupleChain->SetBranchStatus("tracks_vx", 1);
	ntupleChain->SetBranchStatus("tracks_vy", 1);
	ntupleChain->SetBranchStatus("tracks_px", 1);
	ntupleChain->SetBranchStatus("tracks_py", 1);
	ntupleChain->SetBranchStatus("tracks_innerLayerMissingHits", 1);
}

void NTupleReader::LoadTracks() {
	ntupleChain->SetBranchAddress("Ntracks", &Ntracks);
	ntupleChain->SetBranchAddress("tracks_chi2", &tracks_chi2);
	ntupleChain->SetBranchAddress("tracks_ndof", &tracks_ndof);
	ntupleChain->SetBranchAddress("tracks_chg", &tracks_chg);
	ntupleChain->SetBranchAddress("tracks_pt", &tracks_pt);
	ntupleChain->SetBranchAddress("tracks_px", &tracks_px);
	ntupleChain->SetBranchAddress("tracks_py", &tracks_py);
	ntupleChain->SetBranchAddress("tracks_pz", &tracks_pz);
	ntupleChain->SetBranchAddress("tracks_eta", &tracks_eta);
	ntupleChain->SetBranchAddress("tracks_phi", &tracks_phi);
	ntupleChain->SetBranchAddress("tracks_theta", &tracks_theta);
	ntupleChain->SetBranchAddress("tracks_d0dum", &tracks_d0dum);
	ntupleChain->SetBranchAddress("tracks_dz", &tracks_dz);
	ntupleChain->SetBranchAddress("tracks_vx", &tracks_vx);
	ntupleChain->SetBranchAddress("tracks_vy", &tracks_vy);
	ntupleChain->SetBranchAddress("tracks_vz", &tracks_vz);
	ntupleChain->SetBranchAddress("tracks_numvalhits", &tracks_numvalhits);
	ntupleChain->SetBranchAddress("tracks_numlosthits", &tracks_numlosthits);
	ntupleChain->SetBranchAddress("tracks_d0dumErr", &tracks_d0dumErr);
	ntupleChain->SetBranchAddress("tracks_dzErr", &tracks_dzErr);
	ntupleChain->SetBranchAddress("tracks_ptErr", &tracks_ptErr);
	ntupleChain->SetBranchAddress("tracks_etaErr", &tracks_etaErr);
	ntupleChain->SetBranchAddress("tracks_phiErr", &tracks_phiErr);
	ntupleChain->SetBranchAddress("tracks_Nrechits", &tracks_Nrechits);
	ntupleChain->SetBranchAddress("tracks_innerHitX", &tracks_innerHitX);
	ntupleChain->SetBranchAddress("tracks_innerHitY", &tracks_innerHitY);
	ntupleChain->SetBranchAddress("tracks_innerHitZ", &tracks_innerHitZ);
	ntupleChain->SetBranchAddress("tracks_outerHitX", &tracks_outerHitX);
	ntupleChain->SetBranchAddress("tracks_outerHitY", &tracks_outerHitY);
	ntupleChain->SetBranchAddress("tracks_outerHitZ", &tracks_outerHitZ);
	ntupleChain->SetBranchAddress("tracks_highPurity", &tracks_highPurity);
	ntupleChain->SetBranchAddress("tracks_innerLayerMissingHits", &tracks_innerLayerMissingHits);
}
