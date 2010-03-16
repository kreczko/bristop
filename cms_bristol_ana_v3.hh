class ana {

public:
    
  // chain = list of root files containing the same tree
  TChain* chain;
  TChain* chain2;
    
  // The output file with histograms
  TFile* histf;

  //The text file with the candidates
  FILE* outfile;
          
   // Declaration of leaf types
   UInt_t          NPFJets;
   vector<float>   *PFJets_energy;
   vector<float>   *PFJets_et;
   vector<float>   *PFJets_eta;
   vector<float>   *PFJets_phi;
   vector<float>   *PFJets_pt;
   vector<float>   *PFJets_px;
   vector<float>   *PFJets_py;
   vector<float>   *PFJets_pz;
   vector<float>   *PFJets_status;
   vector<float>   *PFJets_theta;
   vector<float>   *PFJets_chgEmE;
   vector<float>   *PFJets_chgHadE;
   vector<float>   *PFJets_chgMuE;
   vector<float>   *PFJets_chg_Mult;
   vector<float>   *PFJets_neutralEmE;
   vector<float>   *PFJets_neutralHadE;
   vector<float>   *PFJets_neutral_Mult;
   vector<float>   *PFJets_mu_Mult;
   vector<float>   *PFJets_mass;
   UInt_t          NPFMets;
   vector<float>   *PFMets_et;
   vector<float>   *PFMets_phi;
   vector<float>   *PFMets_ex;
   vector<float>   *PFMets_ey;
   vector<float>   *PFMets_sumEt;
   UInt_t          NbeamSpot;
   vector<float>   *beamSpot_x;
   vector<float>   *beamSpot_y;
   vector<float>   *beamSpot_z;
   vector<float>   *beamSpot_x0Error;
   vector<float>   *beamSpot_y0Error;
   vector<float>   *beamSpot_z0Error;
   vector<float>   *beamSpot_sigmaZ;
   vector<float>   *beamSpot_sigmaZ0Error;
   vector<float>   *beamSpot_dxdz;
   vector<float>   *beamSpot_dxdzError;
   vector<float>   *beamSpot_dydz;
   vector<float>   *beamSpot_dydzError;
   vector<float>   *beamSpot_beamWidthX;
   vector<float>   *beamSpot_beamWidthY;
   vector<float>   *beamSpot_beamWidthXError;
   vector<float>   *beamSpot_beamWidthYError;
   UInt_t          Nels;
   vector<float>   *els_energy;
   vector<float>   *els_et;
   vector<float>   *els_eta;
   vector<float>   *els_phi;
   vector<float>   *els_pt;
   vector<float>   *els_px;
   vector<float>   *els_py;
   vector<float>   *els_pz;
   vector<float>   *els_status;
   vector<float>   *els_theta;
   vector<float>   *els_closestCtfTrackRef;
   vector<float>   *els_isEcalDriven;
   vector<float>   *els_isTrackerDriven;
   vector<float>   *els_dr03EcalRecHitSumEt;
   vector<float>   *els_dr04EcalRecHitSumEt;
   vector<float>   *els_dr03HcalTowerSumEt;
   vector<float>   *els_dr04HcalTowerSumEt;
   vector<float>   *els_gen_id;
   vector<float>   *els_gen_phi;
   vector<float>   *els_gen_pt;
   vector<float>   *els_gen_pz;
   vector<float>   *els_gen_px;
   vector<float>   *els_gen_py;
   vector<float>   *els_gen_eta;
   vector<float>   *els_gen_theta;
   vector<float>   *els_gen_et;
   vector<float>   *els_gen_mother_id;
   vector<float>   *els_gen_mother_phi;
   vector<float>   *els_gen_mother_pt;
   vector<float>   *els_gen_mother_pz;
   vector<float>   *els_gen_mother_px;
   vector<float>   *els_gen_mother_py;
   vector<float>   *els_gen_mother_eta;
   vector<float>   *els_gen_mother_theta;
   vector<float>   *els_gen_mother_et;
   vector<float>   *els_tightId;
   vector<float>   *els_looseId;
   vector<float>   *els_robustTightId;
   vector<float>   *els_robustLooseId;
   vector<float>   *els_robustHighEnergyId;
   vector<float>   *els_cIso;
   vector<float>   *els_tIso;
   vector<float>   *els_ecalIso;
   vector<float>   *els_hcalIso;
   vector<float>   *els_chi2;
   vector<float>   *els_charge;
   vector<float>   *els_caloEnergy;
   vector<float>   *els_hadOverEm;
   vector<float>   *els_eOverPIn;
   vector<float>   *els_eSeedOverPOut;
   vector<float>   *els_eSCraw;
   vector<float>   *els_eSeed;
   vector<float>   *els_sigmaEtaEta;
   vector<float>   *els_sigmaIEtaIEta;
   vector<float>   *els_scE1x5;
   vector<float>   *els_scE2x5Max;
   vector<float>   *els_scE5x5;
   vector<float>   *els_dEtaIn;
   vector<float>   *els_dPhiIn;
   vector<float>   *els_dEtaOut;
   vector<float>   *els_dPhiOut;
   vector<float>   *els_numvalhits;
   vector<float>   *els_numlosthits;
   vector<float>   *els_basicClustersSize;
   vector<float>   *els_tk_pt;
   vector<float>   *els_tk_phi;
   vector<float>   *els_tk_eta;
   vector<float>   *els_tk_theta;
   vector<float>   *els_tk_charge;
   vector<float>   *els_shFracInnerHits;
   vector<float>   *els_d0dum;
   vector<float>   *els_dz;
   vector<float>   *els_vx;
   vector<float>   *els_vy;
   vector<float>   *els_vz;
   vector<float>   *els_ndof;
   vector<float>   *els_ptError;
   vector<float>   *els_d0dumError;
   vector<float>   *els_dzError;
   vector<float>   *els_etaError;
   vector<float>   *els_phiError;
   vector<float>   *els_cpx;
   vector<float>   *els_cpy;
   vector<float>   *els_cpz;
   vector<float>   *els_vpx;
   vector<float>   *els_vpy;
   vector<float>   *els_vpz;
   vector<float>   *els_cx;
   vector<float>   *els_cy;
   vector<float>   *els_cz;
   vector<float>   *els_innerLayerMissingHits;
   vector<float>   *els_isEE;
   vector<float>   *els_isEEGap;
   vector<float>   *els_isEB;
   vector<float>   *els_isEBGap;
   vector<float>   *els_isConvertedPhoton;
   UInt_t          Njets;
   vector<float>   *jets_energy;
   vector<float>   *jets_et;
   vector<float>   *jets_eta;
   vector<float>   *jets_phi;
   vector<float>   *jets_pt;
   vector<float>   *jets_px;
   vector<float>   *jets_py;
   vector<float>   *jets_pz;
   vector<float>   *jets_status;
   vector<float>   *jets_theta;
   vector<float>   *jets_parton_Id;
   vector<float>   *jets_parton_motherId;
   vector<float>   *jets_parton_pt;
   vector<float>   *jets_parton_phi;
   vector<float>   *jets_parton_eta;
   vector<float>   *jets_parton_Energy;
   vector<float>   *jets_parton_mass;
   vector<float>   *jets_parton_motherID;
   vector<float>   *jets_gen_et;
   vector<float>   *jets_gen_pt;
   vector<float>   *jets_gen_eta;
   vector<float>   *jets_gen_phi;
   vector<float>   *jets_gen_mass;
   vector<float>   *jets_gen_Energy;
   vector<float>   *jets_gen_Id;
   vector<float>   *jets_gen_motherID;
   vector<float>   *jets_gen_threeCharge;
   vector<float>   *jets_partonFlavour;
   vector<float>   *jets_btag_TC_highPur;
   vector<float>   *jets_btag_TC_highEff;
   vector<float>   *jets_btag_jetProb;
   vector<float>   *jets_btag_jetBProb;
   vector<float>   *jets_btag_softEle;
   vector<float>   *jets_btag_softMuon;
   vector<float>   *jets_btag_softMuonNoIP;
   vector<float>   *jets_btag_secVertex;
   vector<float>   *jets_chgEmE;
   vector<float>   *jets_chgHadE;
   vector<float>   *jets_chgMuE;
   vector<float>   *jets_chg_Mult;
   vector<float>   *jets_neutralEmE;
   vector<float>   *jets_neutralHadE;
   vector<float>   *jets_neutral_Mult;
   vector<float>   *jets_mu_Mult;
   vector<float>   *jets_emf;
   vector<float>   *jets_ehf;
   vector<float>   *jets_n60;
   vector<float>   *jets_n90;
   vector<float>   *jets_area;
   vector<float>   *jets_mass;

  UInt_t          NjetsKT4;
  vector<float>   *jetsKT4_energy;
  vector<float>   *jetsKT4_et;
  vector<float>   *jetsKT4_eta;
  vector<float>   *jetsKT4_phi;
  vector<float>   *jetsKT4_pt;
  vector<float>   *jetsKT4_px;
  vector<float>   *jetsKT4_py;
  vector<float>   *jetsKT4_pz;
  vector<float>   *jetsKT4_status;
  vector<float>   *jetsKT4_theta;
  vector<float>   *jetsKT4_btag_TC_highPur;
  vector<float>   *jetsKT4_btag_TC_highEff;
  vector<float>   *jetsKT4_btag_jetProb;
  vector<float>   *jetsKT4_btag_jetBProb;
  vector<float>   *jetsKT4_btag_softEle;
  vector<float>   *jetsKT4_btag_softMuon;
  vector<float>   *jetsKT4_btag_softMuonNoIP;
  vector<float>   *jetsKT4_btag_secVertex;
  vector<float>   *jetsKT4_chgEmE;
  vector<float>   *jetsKT4_chgHadE;
  vector<float>   *jetsKT4_chgMuE;
  vector<float>   *jetsKT4_chg_Mult;
  vector<float>   *jetsKT4_neutralEmE;
  vector<float>   *jetsKT4_neutralHadE;
  vector<float>   *jetsKT4_neutral_Mult;
  vector<float>   *jetsKT4_mu_Mult;
  vector<float>   *jetsKT4_emf;
  vector<float>   *jetsKT4_ehf;
  vector<float>   *jetsKT4_n60;
  vector<float>   *jetsKT4_n90;
  vector<float>   *jetsKT4_area;
  vector<float>   *jetsKT4_mass;

  UInt_t          NjetsKT6;
  vector<float>   *jetsKT6_energy;
  vector<float>   *jetsKT6_et;
  vector<float>   *jetsKT6_eta;
  vector<float>   *jetsKT6_phi;
  vector<float>   *jetsKT6_pt;
  vector<float>   *jetsKT6_px;
  vector<float>   *jetsKT6_py;
  vector<float>   *jetsKT6_pz;
  vector<float>   *jetsKT6_status;
  vector<float>   *jetsKT6_theta;
  vector<float>   *jetsKT6_btag_TC_highPur;
  vector<float>   *jetsKT6_btag_TC_highEff;
  vector<float>   *jetsKT6_btag_jetProb;
  vector<float>   *jetsKT6_btag_jetBProb;
  vector<float>   *jetsKT6_btag_softEle;
  vector<float>   *jetsKT6_btag_softMuon;
  vector<float>   *jetsKT6_btag_softMuonNoIP;
  vector<float>   *jetsKT6_btag_secVertex;
  vector<float>   *jetsKT6_chgEmE;
  vector<float>   *jetsKT6_chgHadE;
  vector<float>   *jetsKT6_chgMuE;
  vector<float>   *jetsKT6_chg_Mult;
  vector<float>   *jetsKT6_neutralEmE;
  vector<float>   *jetsKT6_neutralHadE;
  vector<float>   *jetsKT6_neutral_Mult;
  vector<float>   *jetsKT6_mu_Mult;
  vector<float>   *jetsKT6_emf;
  vector<float>   *jetsKT6_ehf;
  vector<float>   *jetsKT6_n60;
  vector<float>   *jetsKT6_n90;
  vector<float>   *jetsKT6_area;
  vector<float>   *jetsKT6_mass;

   UInt_t          NjetsSC5;
   vector<float>   *jetsSC5_energy;
   vector<float>   *jetsSC5_et;
   vector<float>   *jetsSC5_eta;
   vector<float>   *jetsSC5_phi;
   vector<float>   *jetsSC5_pt;
   vector<float>   *jetsSC5_px;
   vector<float>   *jetsSC5_py;
   vector<float>   *jetsSC5_pz;
   vector<float>   *jetsSC5_status;
   vector<float>   *jetsSC5_theta;
   vector<float>   *jetsSC5_btag_TC_highPur;
   vector<float>   *jetsSC5_btag_TC_highEff;
   vector<float>   *jetsSC5_btag_jetProb;
   vector<float>   *jetsSC5_btag_jetBProb;
   vector<float>   *jetsSC5_btag_softEle;
   vector<float>   *jetsSC5_btag_softMuon;
   vector<float>   *jetsSC5_btag_softMuonNoIP;
   vector<float>   *jetsSC5_btag_secVertex;
   vector<float>   *jetsSC5_chgEmE;
   vector<float>   *jetsSC5_chgHadE;
   vector<float>   *jetsSC5_chgMuE;
   vector<float>   *jetsSC5_chg_Mult;
   vector<float>   *jetsSC5_neutralEmE;
   vector<float>   *jetsSC5_neutralHadE;
   vector<float>   *jetsSC5_neutral_Mult;
   vector<float>   *jetsSC5_mu_Mult;
   vector<float>   *jetsSC5_emf;
   vector<float>   *jetsSC5_ehf;
   vector<float>   *jetsSC5_n60;
   vector<float>   *jetsSC5_n90;
   vector<float>   *jetsSC5_area;
   vector<float>   *jetsSC5_mass;

   UInt_t          NjetsSC7;
   vector<float>   *jetsSC7_energy;
   vector<float>   *jetsSC7_et;
   vector<float>   *jetsSC7_eta;
   vector<float>   *jetsSC7_phi;
   vector<float>   *jetsSC7_pt;
   vector<float>   *jetsSC7_px;
   vector<float>   *jetsSC7_py;
   vector<float>   *jetsSC7_pz;
   vector<float>   *jetsSC7_status;
   vector<float>   *jetsSC7_theta;
   vector<float>   *jetsSC7_btag_TC_highPur;
   vector<float>   *jetsSC7_btag_TC_highEff;
   vector<float>   *jetsSC7_btag_jetProb;
   vector<float>   *jetsSC7_btag_jetBProb;
   vector<float>   *jetsSC7_btag_softEle;
   vector<float>   *jetsSC7_btag_softMuon;
   vector<float>   *jetsSC7_btag_softMuonNoIP;
   vector<float>   *jetsSC7_btag_secVertex;
   vector<float>   *jetsSC7_chgEmE;
   vector<float>   *jetsSC7_chgHadE;
   vector<float>   *jetsSC7_chgMuE;
   vector<float>   *jetsSC7_chg_Mult;
   vector<float>   *jetsSC7_neutralEmE;
   vector<float>   *jetsSC7_neutralHadE;
   vector<float>   *jetsSC7_neutral_Mult;
   vector<float>   *jetsSC7_mu_Mult;
   vector<float>   *jetsSC7_emf;
   vector<float>   *jetsSC7_ehf;
   vector<float>   *jetsSC7_n60;
   vector<float>   *jetsSC7_n90;
   vector<float>   *jetsSC7_area;
   vector<float>   *jetsSC7_mass;

   UInt_t          NjetsJPTAK5;
   vector<float>   *jetsJPTAK5_energy;
   vector<float>   *jetsJPTAK5_et;
   vector<float>   *jetsJPTAK5_eta;
   vector<float>   *jetsJPTAK5_phi;
   vector<float>   *jetsJPTAK5_pt;
   vector<float>   *jetsJPTAK5_px;
   vector<float>   *jetsJPTAK5_py;
   vector<float>   *jetsJPTAK5_pz;
   vector<float>   *jetsJPTAK5_status;
   vector<float>   *jetsJPTAK5_theta;
   vector<float>   *jetsJPTAK5_btag_TC_highPur;
   vector<float>   *jetsJPTAK5_btag_TC_highEff;
   vector<float>   *jetsJPTAK5_btag_jetProb;
   vector<float>   *jetsJPTAK5_btag_jetBProb;
   vector<float>   *jetsJPTAK5_btag_softEle;
   vector<float>   *jetsJPTAK5_btag_softMuon;
   vector<float>   *jetsJPTAK5_btag_softMuonNoIP;
   vector<float>   *jetsJPTAK5_btag_secVertex;
   vector<float>   *jetsJPTAK5_chgEmE;
   vector<float>   *jetsJPTAK5_chgHadE;
   vector<float>   *jetsJPTAK5_chgMuE;
   vector<float>   *jetsJPTAK5_chg_Mult;
   vector<float>   *jetsJPTAK5_neutralEmE;
   vector<float>   *jetsJPTAK5_neutralHadE;
   vector<float>   *jetsJPTAK5_neutral_Mult;
   vector<float>   *jetsJPTAK5_mu_Mult;
   vector<float>   *jetsJPTAK5_emf;
   vector<float>   *jetsJPTAK5_ehf;
   vector<float>   *jetsJPTAK5_n60;
   vector<float>   *jetsJPTAK5_n90;
   vector<float>   *jetsJPTAK5_area;
   vector<float>   *jetsJPTAK5_mass;

   UInt_t          Nmc_doc;
   vector<float>   *mc_doc_id;
   vector<float>   *mc_doc_pt;
   vector<float>   *mc_doc_px;
   vector<float>   *mc_doc_py;
   vector<float>   *mc_doc_pz;
   vector<float>   *mc_doc_eta;
   vector<float>   *mc_doc_phi;
   vector<float>   *mc_doc_theta;
   vector<float>   *mc_doc_energy;
   vector<float>   *mc_doc_status;
   vector<float>   *mc_doc_charge;
   vector<float>   *mc_doc_mother_id;
   vector<float>   *mc_doc_grandmother_id;
   vector<float>   *mc_doc_ggrandmother_id;
   vector<float>   *mc_doc_mother_pt;
   vector<float>   *mc_doc_vertex_x;
   vector<float>   *mc_doc_vertex_y;
   vector<float>   *mc_doc_vertex_z;
   vector<float>   *mc_doc_mass;
   vector<float>   *mc_doc_numOfDaughters;
   vector<float>   *mc_doc_numOfMothers;
   UInt_t          Nmets;
   vector<float>   *mets_et;
   vector<float>   *mets_phi;
   vector<float>   *mets_ex;
   vector<float>   *mets_ey;
   vector<float>   *mets_gen_et;
   vector<float>   *mets_gen_phi;
   vector<float>   *mets_sign;
   vector<float>   *mets_sumEt;
   vector<float>   *mets_unCPhi;
   vector<float>   *mets_unCPt;
   vector<float>   *mets_et_muonCor;
   vector<float>   *mets_phi_muonCor;
   vector<float>   *mets_et_JESCor;
   vector<float>   *mets_phi_JESCor;
   UInt_t          NmetsKT4;
   vector<float>   *metsKT4_et;
   vector<float>   *metsKT4_phi;
   vector<float>   *metsKT4_ex;
   vector<float>   *metsKT4_ey;
   vector<float>   *metsKT4_sumEt;
   vector<float>   *metsKT4_et_JESCor;
   vector<float>   *metsKT4_phi_JESCor;
   UInt_t          NmetsKT6;
   vector<float>   *metsKT6_et;
   vector<float>   *metsKT6_phi;
   vector<float>   *metsKT6_ex;
   vector<float>   *metsKT6_ey;
   vector<float>   *metsKT6_sumEt;
   vector<float>   *metsKT6_et_JESCor;
   vector<float>   *metsKT6_phi_JESCor;
   UInt_t          NmetsSC5;
   vector<float>   *metsSC5_et;
   vector<float>   *metsSC5_phi;
   vector<float>   *metsSC5_ex;
   vector<float>   *metsSC5_ey;
   vector<float>   *metsSC5_sumEt;
   vector<float>   *metsSC5_et_JESCor;
   vector<float>   *metsSC5_phi_JESCor;
   UInt_t          NmetsSC7;
   vector<float>   *metsSC7_et;
   vector<float>   *metsSC7_phi;
   vector<float>   *metsSC7_ex;
   vector<float>   *metsSC7_ey;
   vector<float>   *metsSC7_sumEt;
   vector<float>   *metsSC7_et_JESCor;
   vector<float>   *metsSC7_phi_JESCor;
   UInt_t          Nmus;
   vector<float>   *mus_energy;
   vector<float>   *mus_et;
   vector<float>   *mus_eta;
   vector<float>   *mus_phi;
   vector<float>   *mus_pt;
   vector<float>   *mus_px;
   vector<float>   *mus_py;
   vector<float>   *mus_pz;
   vector<float>   *mus_status;
   vector<float>   *mus_theta;
   vector<float>   *mus_gen_id;
   vector<float>   *mus_gen_phi;
   vector<float>   *mus_gen_pt;
   vector<float>   *mus_gen_pz;
   vector<float>   *mus_gen_px;
   vector<float>   *mus_gen_py;
   vector<float>   *mus_gen_eta;
   vector<float>   *mus_gen_theta;
   vector<float>   *mus_gen_et;
   vector<float>   *mus_gen_mother_id;
   vector<float>   *mus_gen_mother_phi;
   vector<float>   *mus_gen_mother_pt;
   vector<float>   *mus_gen_mother_pz;
   vector<float>   *mus_gen_mother_px;
   vector<float>   *mus_gen_mother_py;
   vector<float>   *mus_gen_mother_eta;
   vector<float>   *mus_gen_mother_theta;
   vector<float>   *mus_gen_mother_et;
   vector<float>   *mus_tkHits;
   vector<float>   *mus_cIso;
   vector<float>   *mus_tIso;
   vector<float>   *mus_ecalIso;
   vector<float>   *mus_hcalIso;
   vector<float>   *mus_ecalvetoDep;
   vector<float>   *mus_hcalvetoDep;
   vector<float>   *mus_calEnergyEm;
   vector<float>   *mus_calEnergyHad;
   vector<float>   *mus_calEnergyHo;
   vector<float>   *mus_calEnergyEmS9;
   vector<float>   *mus_calEnergyHadS9;
   vector<float>   *mus_calEnergyHoS9;
   vector<float>   *mus_iso03_sumPt;
   vector<float>   *mus_iso03_emEt;
   vector<float>   *mus_iso03_hadEt;
   vector<float>   *mus_iso03_hoEt;
   vector<float>   *mus_iso03_nTracks;
   vector<float>   *mus_iso05_sumPt;
   vector<float>   *mus_iso05_emEt;
   vector<float>   *mus_iso05_hadEt;
   vector<float>   *mus_iso05_hoEt;
   vector<float>   *mus_iso05_nTracks;
   vector<float>   *mus_charge;
   vector<float>   *mus_cm_chi2;
   vector<float>   *mus_cm_ndof;
   vector<float>   *mus_cm_chg;
   vector<float>   *mus_cm_pt;
   vector<float>   *mus_cm_px;
   vector<float>   *mus_cm_py;
   vector<float>   *mus_cm_pz;
   vector<float>   *mus_cm_eta;
   vector<float>   *mus_cm_phi;
   vector<float>   *mus_cm_theta;
   vector<float>   *mus_cm_d0dum;
   vector<float>   *mus_cm_dz;
   vector<float>   *mus_cm_vx;
   vector<float>   *mus_cm_vy;
   vector<float>   *mus_cm_vz;
   vector<float>   *mus_cm_numvalhits;
   vector<float>   *mus_cm_numlosthits;
   vector<float>   *mus_cm_d0dumErr;
   vector<float>   *mus_cm_dzErr;
   vector<float>   *mus_cm_ptErr;
   vector<float>   *mus_cm_etaErr;
   vector<float>   *mus_cm_phiErr;
   vector<float>   *mus_tk_chi2;
   vector<float>   *mus_tk_ndof;
   vector<float>   *mus_tk_chg;
   vector<float>   *mus_tk_pt;
   vector<float>   *mus_tk_px;
   vector<float>   *mus_tk_py;
   vector<float>   *mus_tk_pz;
   vector<float>   *mus_tk_eta;
   vector<float>   *mus_tk_phi;
   vector<float>   *mus_tk_theta;
   vector<float>   *mus_tk_d0dum;
   vector<float>   *mus_tk_dz;
   vector<float>   *mus_tk_vx;
   vector<float>   *mus_tk_vy;
   vector<float>   *mus_tk_vz;
   vector<float>   *mus_tk_numvalhits;
   vector<float>   *mus_tk_numlosthits;
   vector<float>   *mus_tk_d0dumErr;
   vector<float>   *mus_tk_dzErr;
   vector<float>   *mus_tk_ptErr;
   vector<float>   *mus_tk_etaErr;
   vector<float>   *mus_tk_phiErr;
   vector<float>   *mus_stamu_chi2;
   vector<float>   *mus_stamu_ndof;
   vector<float>   *mus_stamu_chg;
   vector<float>   *mus_stamu_pt;
   vector<float>   *mus_stamu_px;
   vector<float>   *mus_stamu_py;
   vector<float>   *mus_stamu_pz;
   vector<float>   *mus_stamu_eta;
   vector<float>   *mus_stamu_phi;
   vector<float>   *mus_stamu_theta;
   vector<float>   *mus_stamu_d0dum;
   vector<float>   *mus_stamu_dz;
   vector<float>   *mus_stamu_vx;
   vector<float>   *mus_stamu_vy;
   vector<float>   *mus_stamu_vz;
   vector<float>   *mus_stamu_numvalhits;
   vector<float>   *mus_stamu_numlosthits;
   vector<float>   *mus_stamu_d0dumErr;
   vector<float>   *mus_stamu_dzErr;
   vector<float>   *mus_stamu_ptErr;
   vector<float>   *mus_stamu_etaErr;
   vector<float>   *mus_stamu_phiErr;
   vector<float>   *mus_num_matches;
   vector<float>   *mus_id_All;
   vector<float>   *mus_id_AllGlobalMuons;
   vector<float>   *mus_id_AllStandAloneMuons;
   vector<float>   *mus_id_AllTrackerMuons;
   vector<float>   *mus_id_TrackerMuonArbitrated;
   vector<float>   *mus_id_AllArbitrated;
   vector<float>   *mus_id_GlobalMuonPromptTight;
   vector<float>   *mus_id_TMLastStationLoose;
   vector<float>   *mus_id_TMLastStationTight;
   vector<float>   *mus_id_TM2DCompatibilityLoose;
   vector<float>   *mus_id_TM2DCompatibilityTight;
   vector<float>   *mus_id_TMOneStationLoose;
   vector<float>   *mus_id_TMOneStationTight;
   vector<float>   *mus_id_TMLastStationOptimizedLowPtLoose;
   vector<float>   *mus_id_TMLastStationOptimizedLowPtTight;
   UInt_t          Nphotons;
   vector<float>   *photons_energy;
   vector<float>   *photons_et;
   vector<float>   *photons_eta;
   vector<float>   *photons_phi;
   vector<float>   *photons_pt;
   vector<float>   *photons_px;
   vector<float>   *photons_py;
   vector<float>   *photons_pz;
   vector<float>   *photons_status;
   vector<float>   *photons_theta;
   vector<float>   *photons_hadOverEM;
   vector<float>   *photons_scEnergy;
   vector<float>   *photons_scRawEnergy;
   vector<float>   *photons_scEta;
   vector<float>   *photons_scPhi;
   vector<float>   *photons_scEtaWidth;
   vector<float>   *photons_scPhiWidth;
   vector<float>   *photons_tIso;
   vector<float>   *photons_ecalIso;
   vector<float>   *photons_hcalIso;
   vector<float>   *photons_isoEcalRecHitDR04;
   vector<float>   *photons_isoHcalRecHitDR04;
   vector<float>   *photons_isoSolidTrkConeDR04;
   vector<float>   *photons_isoHollowTrkConeDR04;
   vector<float>   *photons_nTrkSolidConeDR04;
   vector<float>   *photons_nTrkHollowConeDR04;
   vector<float>   *photons_isoEcalRecHitDR03;
   vector<float>   *photons_isoHcalRecHitDR03;
   vector<float>   *photons_isoSolidTrkConeDR03;
   vector<float>   *photons_isoHollowTrkConeDR03;
   vector<float>   *photons_nTrkSolidConeDR03;
   vector<float>   *photons_nTrkHollowConeDR03;
   vector<float>   *photons_isAlsoElectron;
   vector<float>   *photons_hasPixelSeed;
   vector<float>   *photons_isConverted;
   vector<float>   *photons_isEBGap;
   vector<float>   *photons_isEEGap;
   vector<float>   *photons_isEBEEGap;
   vector<float>   *photons_isEBPho;
   vector<float>   *photons_isEEPho;
   vector<float>   *photons_isLoosePhoton;
   vector<float>   *photons_isTightPhoton;
   vector<float>   *photons_r9;
   vector<float>   *photons_gen_et;
   vector<float>   *photons_gen_eta;
   vector<float>   *photons_gen_phi;
   vector<float>   *photons_gen_id;
   UInt_t          Npv;
   vector<float>   *pv_x;
   vector<float>   *pv_y;
   vector<float>   *pv_z;
   vector<float>   *pv_xErr;
   vector<float>   *pv_yErr;
   vector<float>   *pv_zErr;
   UInt_t          Ntcmets;
   vector<float>   *tcmets_et;
   vector<float>   *tcmets_phi;
   vector<float>   *tcmets_ex;
   vector<float>   *tcmets_ey;
   vector<float>   *tcmets_sumEt;
   vector<float>   *tcmets_et_JESCor;
   vector<float>   *tcmets_phi_JESCor;
   vector<float>   *tcmets_et_muonCor;
   vector<float>   *tcmets_phi_muonCor;
   UInt_t          Ntracks;
   vector<float>   *tracks_chi2;
   vector<float>   *tracks_ndof;
   vector<float>   *tracks_chg;
   vector<float>   *tracks_pt;
   vector<float>   *tracks_px;
   vector<float>   *tracks_py;
   vector<float>   *tracks_pz;
   vector<float>   *tracks_eta;
   vector<float>   *tracks_phi;
   vector<float>   *tracks_theta;
   vector<float>   *tracks_d0dum;
   vector<float>   *tracks_dz;
   vector<float>   *tracks_vx;
   vector<float>   *tracks_vy;
   vector<float>   *tracks_vz;
   vector<float>   *tracks_numvalhits;
   vector<float>   *tracks_numlosthits;
   vector<float>   *tracks_d0dumErr;
   vector<float>   *tracks_dzErr;
   vector<float>   *tracks_ptErr;
   vector<float>   *tracks_etaErr;
   vector<float>   *tracks_phiErr;
   vector<float>   *tracks_Nrechits;
   vector<float>   *tracks_innerHitX;
   vector<float>   *tracks_innerHitY;
   vector<float>   *tracks_innerHitZ;
   vector<float>   *tracks_outerHitX;
   vector<float>   *tracks_outerHitY;
   vector<float>   *tracks_outerHitZ;
   vector<float>   *tracks_highPurity;
   vector<float>    *tracks_innerLayerMissingHits;
   UInt_t          run;
   UInt_t          event;
   UInt_t          lumiBlock;

   // List of branches
   //TBranch        *b_NTEST;   //! TEST
   //TBranch        *b_TEST_x;   //! TEST

   TBranch        *b_NPFJets;   //!
   TBranch        *b_PFJets_energy;   //!
   TBranch        *b_PFJets_et;   //!
   TBranch        *b_PFJets_eta;   //!
   TBranch        *b_PFJets_phi;   //!
   TBranch        *b_PFJets_pt;   //!
   TBranch        *b_PFJets_px;   //!
   TBranch        *b_PFJets_py;   //!
   TBranch        *b_PFJets_pz;   //!
   TBranch        *b_PFJets_status;   //!
   TBranch        *b_PFJets_theta;   //!
   TBranch        *b_PFJets_chgEmE;   //!
   TBranch        *b_PFJets_chgHadE;   //!
   TBranch        *b_PFJets_chgMuE;   //!
   TBranch        *b_PFJets_chg_Mult;   //!
   TBranch        *b_PFJets_neutralEmE;   //!
   TBranch        *b_PFJets_neutralHadE;   //!
   TBranch        *b_PFJets_neutral_Mult;   //!
   TBranch        *b_PFJets_mu_Mult;   //!
   TBranch        *b_PFJets_mass;   //!
   TBranch        *b_NPFMets;   //!
   TBranch        *b_PFMets_et;   //!
   TBranch        *b_PFMets_phi;   //!
   TBranch        *b_PFMets_ex;   //!
   TBranch        *b_PFMets_ey;   //!
   TBranch        *b_PFMets_sumEt;   //!
   TBranch        *b_NbeamSpot;   //!
   TBranch        *b_beamSpot_x;   //!
   TBranch        *b_beamSpot_y;   //!
   TBranch        *b_beamSpot_z;   //!
   TBranch        *b_beamSpot_x0Error;   //!
   TBranch        *b_beamSpot_y0Error;   //!
   TBranch        *b_beamSpot_z0Error;   //!
   TBranch        *b_beamSpot_sigmaZ;   //!
   TBranch        *b_beamSpot_sigmaZ0Error;   //!
   TBranch        *b_beamSpot_dxdz;   //!
   TBranch        *b_beamSpot_dxdzError;   //!
   TBranch        *b_beamSpot_dydz;   //!
   TBranch        *b_beamSpot_dydzError;   //!
   TBranch        *b_beamSpot_beamWidthX;   //!
   TBranch        *b_beamSpot_beamWidthY;   //!
   TBranch        *b_beamSpot_beamWidthXError;   //!
   TBranch        *b_beamSpot_beamWidthYError;   //!
   TBranch        *b_Nels;   //!
   TBranch        *b_els_energy;   //!
   TBranch        *b_els_et;   //!
   TBranch        *b_els_eta;   //!
   TBranch        *b_els_phi;   //!
   TBranch        *b_els_pt;   //!
   TBranch        *b_els_px;   //!
   TBranch        *b_els_py;   //!
   TBranch        *b_els_pz;   //!
   TBranch        *b_els_status;   //!
   TBranch        *b_els_theta;   //!
   TBranch        *b_els_closestCtfTrackRef;   //!
   TBranch        *b_els_isEcalDriven;   //!
   TBranch        *b_els_isTrackerDriven;   //!
   TBranch        *b_els_dr03EcalRecHitSumEt;   //!
   TBranch        *b_els_dr04EcalRecHitSumEt;   //!
   TBranch        *b_els_dr03HcalTowerSumEt;   //!
   TBranch        *b_els_dr04HcalTowerSumEt;   //!
   TBranch        *b_els_gen_id;   //!
   TBranch        *b_els_gen_phi;   //!
   TBranch        *b_els_gen_pt;   //!
   TBranch        *b_els_gen_pz;   //!
   TBranch        *b_els_gen_px;   //!
   TBranch        *b_els_gen_py;   //!
   TBranch        *b_els_gen_eta;   //!
   TBranch        *b_els_gen_theta;   //!
   TBranch        *b_els_gen_et;   //!
   TBranch        *b_els_gen_mother_id;   //!
   TBranch        *b_els_gen_mother_phi;   //!
   TBranch        *b_els_gen_mother_pt;   //!
   TBranch        *b_els_gen_mother_pz;   //!
   TBranch        *b_els_gen_mother_px;   //!
   TBranch        *b_els_gen_mother_py;   //!
   TBranch        *b_els_gen_mother_eta;   //!
   TBranch        *b_els_gen_mother_theta;   //!
   TBranch        *b_els_gen_mother_et;   //!
   TBranch        *b_els_tightId;   //!
   TBranch        *b_els_looseId;   //!
   TBranch        *b_els_robustTightId;   //!
   TBranch        *b_els_robustLooseId;   //!
   TBranch        *b_els_robustHighEnergyId;   //!
   TBranch        *b_els_cIso;   //!
   TBranch        *b_els_tIso;   //!
   TBranch        *b_els_ecalIso;   //!
   TBranch        *b_els_hcalIso;   //!
   TBranch        *b_els_chi2;   //!
   TBranch        *b_els_charge;   //!
   TBranch        *b_els_caloEnergy;   //!
   TBranch        *b_els_hadOverEm;   //!
   TBranch        *b_els_eOverPIn;   //!
   TBranch        *b_els_eSeedOverPOut;   //!
   TBranch        *b_els_eSCraw;   //!
   TBranch        *b_els_eSeed;   //!
   TBranch        *b_els_sigmaEtaEta;   //!
   TBranch        *b_els_sigmaIEtaIEta;   //!
   TBranch        *b_els_scE1x5;   //!
   TBranch        *b_els_scE2x5Max;   //!
   TBranch        *b_els_scE5x5;   //!
   TBranch        *b_els_dEtaIn;   //!
   TBranch        *b_els_dPhiIn;   //!
   TBranch        *b_els_dEtaOut;   //!
   TBranch        *b_els_dPhiOut;   //!
   TBranch        *b_els_numvalhits;   //!
   TBranch        *b_els_numlosthits;   //!
   TBranch        *b_els_basicClustersSize;   //!
   TBranch        *b_els_tk_pt;   //!
   TBranch        *b_els_tk_phi;   //!
   TBranch        *b_els_tk_eta;   //!
   TBranch        *b_els_tk_charge;   //!
   TBranch        *b_els_tk_theta;   //!     
   TBranch        *b_els_shFracInnerHits; //!     
   TBranch        *b_els_d0dum;   //!
   TBranch        *b_els_dz;   //!
   TBranch        *b_els_vx;   //!
   TBranch        *b_els_vy;   //!
   TBranch        *b_els_vz;   //!
   TBranch        *b_els_ndof;   //!
   TBranch        *b_els_ptError;   //!
   TBranch        *b_els_d0dumError;   //!
   TBranch        *b_els_dzError;   //!
   TBranch        *b_els_etaError;   //!
   TBranch        *b_els_phiError;   //!
   TBranch        *b_els_cpx;   //!
   TBranch        *b_els_cpy;   //!
   TBranch        *b_els_cpz;   //!
   TBranch        *b_els_vpx;   //!
   TBranch        *b_els_vpy;   //!
   TBranch        *b_els_vpz;   //!
   TBranch        *b_els_cx;   //!
   TBranch        *b_els_cy;   //!
   TBranch        *b_els_cz;   //!
   TBranch        *b_els_innerLayerMissingHits;  //! 
   TBranch        *b_els_isEE;   //!
   TBranch        *b_els_isEEGap;   //!
   TBranch        *b_els_isEB;   //!
   TBranch        *b_els_isEBGap;   //!
   TBranch        *b_els_isConvertedPhoton;   //!
   TBranch        *b_Njets;   //!
   TBranch        *b_jets_energy;   //!
   TBranch        *b_jets_et;   //!
   TBranch        *b_jets_eta;   //!
   TBranch        *b_jets_phi;   //!
   TBranch        *b_jets_pt;   //!
   TBranch        *b_jets_px;   //!
   TBranch        *b_jets_py;   //!
   TBranch        *b_jets_pz;   //!
   TBranch        *b_jets_status;   //!
   TBranch        *b_jets_theta;   //!
   TBranch        *b_jets_parton_Id;   //!
   TBranch        *b_jets_parton_motherId;   //!
   TBranch        *b_jets_parton_pt;   //!
   TBranch        *b_jets_parton_phi;   //!
   TBranch        *b_jets_parton_eta;   //!
   TBranch        *b_jets_parton_Energy;   //!
   TBranch        *b_jets_parton_mass;   //!
   TBranch        *b_jets_parton_motherID;   //!
   TBranch        *b_jets_gen_et;   //!
   TBranch        *b_jets_gen_pt;   //!
   TBranch        *b_jets_gen_eta;   //!
   TBranch        *b_jets_gen_phi;   //!
   TBranch        *b_jets_gen_mass;   //!
   TBranch        *b_jets_gen_Energy;   //!
   TBranch        *b_jets_gen_Id;   //!
   TBranch        *b_jets_gen_motherID;   //!
   TBranch        *b_jets_gen_threeCharge;   //!
   TBranch        *b_jets_partonFlavour;   //!
   TBranch        *b_jets_btag_TC_highPur;   //!
   TBranch        *b_jets_btag_TC_highEff;   //!
   TBranch        *b_jets_btag_jetProb;   //!
   TBranch        *b_jets_btag_jetBProb;   //!
   TBranch        *b_jets_btag_softEle;   //!
   TBranch        *b_jets_btag_softMuon;   //!
   TBranch        *b_jets_btag_softMuonNoIP;   //!
   TBranch        *b_jets_btag_secVertex;   //!
   TBranch        *b_jets_chgEmE;   //!
   TBranch        *b_jets_chgHadE;   //!
   TBranch        *b_jets_chgMuE;   //!
   TBranch        *b_jets_chg_Mult;   //!
   TBranch        *b_jets_neutralEmE;   //!
   TBranch        *b_jets_neutralHadE;   //!
   TBranch        *b_jets_neutral_Mult;   //!
   TBranch        *b_jets_mu_Mult;   //!
   TBranch        *b_jets_emf;   //!
   TBranch        *b_jets_ehf;   //!
   TBranch        *b_jets_n60;   //!
   TBranch        *b_jets_n90;   //!
   TBranch        *b_jets_area;   //!
   TBranch        *b_jets_mass;   //!

   TBranch        *b_NjetsKT4;   //!
   TBranch        *b_jetsKT4_energy;   //!
   TBranch        *b_jetsKT4_et;   //!
   TBranch        *b_jetsKT4_eta;   //!
   TBranch        *b_jetsKT4_phi;   //!
   TBranch        *b_jetsKT4_pt;   //!
   TBranch        *b_jetsKT4_px;   //!
   TBranch        *b_jetsKT4_py;   //!
   TBranch        *b_jetsKT4_pz;   //!
   TBranch        *b_jetsKT4_status;   //!
   TBranch        *b_jetsKT4_theta;   //!
   TBranch        *b_jetsKT4_btag_TC_highPur;   //!
   TBranch        *b_jetsKT4_btag_TC_highEff;   //!
   TBranch        *b_jetsKT4_btag_jetProb;   //!
   TBranch        *b_jetsKT4_btag_jetBProb;   //!
   TBranch        *b_jetsKT4_btag_softEle;   //!
   TBranch        *b_jetsKT4_btag_softMuon;   //!
   TBranch        *b_jetsKT4_btag_softMuonNoIP;   //!
   TBranch        *b_jetsKT4_btag_secVertex;   //!
   TBranch        *b_jetsKT4_chgEmE;   //!
   TBranch        *b_jetsKT4_chgHadE;   //!
   TBranch        *b_jetsKT4_chgMuE;   //!
   TBranch        *b_jetsKT4_chg_Mult;   //!
   TBranch        *b_jetsKT4_neutralEmE;   //!
   TBranch        *b_jetsKT4_neutralHadE;   //!
   TBranch        *b_jetsKT4_neutral_Mult;   //!
   TBranch        *b_jetsKT4_mu_Mult;   //!
   TBranch        *b_jetsKT4_emf;   //!
   TBranch        *b_jetsKT4_ehf;   //!
   TBranch        *b_jetsKT4_n60;   //!
   TBranch        *b_jetsKT4_n90;   //!
   TBranch        *b_jetsKT4_area;   //!
   TBranch        *b_jetsKT4_mass;   //!

   TBranch        *b_NjetsKT6;   //!
   TBranch        *b_jetsKT6_energy;   //!
   TBranch        *b_jetsKT6_et;   //!
   TBranch        *b_jetsKT6_eta;   //!
   TBranch        *b_jetsKT6_phi;   //!
   TBranch        *b_jetsKT6_pt;   //!
   TBranch        *b_jetsKT6_px;   //!
   TBranch        *b_jetsKT6_py;   //!
   TBranch        *b_jetsKT6_pz;   //!
   TBranch        *b_jetsKT6_status;   //!
   TBranch        *b_jetsKT6_theta;   //!
   TBranch        *b_jetsKT6_btag_TC_highPur;   //!
   TBranch        *b_jetsKT6_btag_TC_highEff;   //!
   TBranch        *b_jetsKT6_btag_jetProb;   //!
   TBranch        *b_jetsKT6_btag_jetBProb;   //!
   TBranch        *b_jetsKT6_btag_softEle;   //!
   TBranch        *b_jetsKT6_btag_softMuon;   //!
   TBranch        *b_jetsKT6_btag_softMuonNoIP;   //!
   TBranch        *b_jetsKT6_btag_secVertex;   //!
   TBranch        *b_jetsKT6_chgEmE;   //!
   TBranch        *b_jetsKT6_chgHadE;   //!
   TBranch        *b_jetsKT6_chgMuE;   //!
   TBranch        *b_jetsKT6_chg_Mult;   //!
   TBranch        *b_jetsKT6_neutralEmE;   //!
   TBranch        *b_jetsKT6_neutralHadE;   //!
   TBranch        *b_jetsKT6_neutral_Mult;   //!
   TBranch        *b_jetsKT6_mu_Mult;   //!
   TBranch        *b_jetsKT6_emf;   //!
   TBranch        *b_jetsKT6_ehf;   //!
   TBranch        *b_jetsKT6_n60;   //!
   TBranch        *b_jetsKT6_n90;   //!
   TBranch        *b_jetsKT6_area;   //!
   TBranch        *b_jetsKT6_mass;   //!

   ///3099
   TBranch        *b_NjetsSC5;   //!
   TBranch        *b_jetsSC5_energy;   //!
   TBranch        *b_jetsSC5_et;   //!
   TBranch        *b_jetsSC5_eta;   //!
   TBranch        *b_jetsSC5_phi;   //!
   TBranch        *b_jetsSC5_pt;   //!
   TBranch        *b_jetsSC5_px;   //!
   TBranch        *b_jetsSC5_py;   //!
   TBranch        *b_jetsSC5_pz;   //!
   TBranch        *b_jetsSC5_status;   //!
   TBranch        *b_jetsSC5_theta;   //!
   TBranch        *b_jetsSC5_btag_TC_highPur;   //!
   TBranch        *b_jetsSC5_btag_TC_highEff;   //!
   TBranch        *b_jetsSC5_btag_jetProb;   //!
   TBranch        *b_jetsSC5_btag_jetBProb;   //!
   TBranch        *b_jetsSC5_btag_softEle;   //!
   TBranch        *b_jetsSC5_btag_softMuon;   //!
   TBranch        *b_jetsSC5_btag_softMuonNoIP;   //!
   TBranch        *b_jetsSC5_btag_secVertex;   //!
   TBranch        *b_jetsSC5_chgEmE;   //!
   TBranch        *b_jetsSC5_chgHadE;   //!
   TBranch        *b_jetsSC5_chgMuE;   //!
   TBranch        *b_jetsSC5_chg_Mult;   //!
   TBranch        *b_jetsSC5_neutralEmE;   //!
   TBranch        *b_jetsSC5_neutralHadE;   //!
   TBranch        *b_jetsSC5_neutral_Mult;   //!
   TBranch        *b_jetsSC5_mu_Mult;   //!
   TBranch        *b_jetsSC5_emf;   //!
   TBranch        *b_jetsSC5_ehf;   //!
   TBranch        *b_jetsSC5_n60;   //!
   TBranch        *b_jetsSC5_n90;   //!
   TBranch        *b_jetsSC5_area;   //!
   TBranch        *b_jetsSC5_mass;   //!

   TBranch        *b_NjetsSC7;   //!
   TBranch        *b_jetsSC7_energy;   //!
   TBranch        *b_jetsSC7_et;   //!
   TBranch        *b_jetsSC7_eta;   //!
   TBranch        *b_jetsSC7_phi;   //!
   TBranch        *b_jetsSC7_pt;   //!
   TBranch        *b_jetsSC7_px;   //!
   TBranch        *b_jetsSC7_py;   //!
   TBranch        *b_jetsSC7_pz;   //!
   TBranch        *b_jetsSC7_status;   //!
   TBranch        *b_jetsSC7_theta;   //!
   TBranch        *b_jetsSC7_btag_TC_highPur;   //!
   TBranch        *b_jetsSC7_btag_TC_highEff;   //!
   TBranch        *b_jetsSC7_btag_jetProb;   //!
   TBranch        *b_jetsSC7_btag_jetBProb;   //!
   TBranch        *b_jetsSC7_btag_softEle;   //!
   TBranch        *b_jetsSC7_btag_softMuon;   //!
   TBranch        *b_jetsSC7_btag_softMuonNoIP;   //!
   TBranch        *b_jetsSC7_btag_secVertex;   //!
   TBranch        *b_jetsSC7_chgEmE;   //!
   TBranch        *b_jetsSC7_chgHadE;   //!
   TBranch        *b_jetsSC7_chgMuE;   //!
   TBranch        *b_jetsSC7_chg_Mult;   //!
   TBranch        *b_jetsSC7_neutralEmE;   //!
   TBranch        *b_jetsSC7_neutralHadE;   //!
   TBranch        *b_jetsSC7_neutral_Mult;   //!
   TBranch        *b_jetsSC7_mu_Mult;   //!
   TBranch        *b_jetsSC7_emf;   //!
   TBranch        *b_jetsSC7_ehf;   //!
   TBranch        *b_jetsSC7_n60;   //!
   TBranch        *b_jetsSC7_n90;   //!
   TBranch        *b_jetsSC7_area;   //!
   TBranch        *b_jetsSC7_mass;   //!

   TBranch        *b_NjetsJPTAK5;   //!
   TBranch        *b_jetsJPTAK5_energy;   //!
   TBranch        *b_jetsJPTAK5_et;   //!
   TBranch        *b_jetsJPTAK5_eta;   //!
   TBranch        *b_jetsJPTAK5_phi;   //!
   TBranch        *b_jetsJPTAK5_pt;   //!
   TBranch        *b_jetsJPTAK5_px;   //!
   TBranch        *b_jetsJPTAK5_py;   //!
   TBranch        *b_jetsJPTAK5_pz;   //!
   TBranch        *b_jetsJPTAK5_status;   //!
   TBranch        *b_jetsJPTAK5_theta;   //!
   TBranch        *b_jetsJPTAK5_btag_TC_highPur;   //!
   TBranch        *b_jetsJPTAK5_btag_TC_highEff;   //!
   TBranch        *b_jetsJPTAK5_btag_jetProb;   //!
   TBranch        *b_jetsJPTAK5_btag_jetBProb;   //!
   TBranch        *b_jetsJPTAK5_btag_softEle;   //!
   TBranch        *b_jetsJPTAK5_btag_softMuon;   //!
   TBranch        *b_jetsJPTAK5_btag_softMuonNoIP;   //!
   TBranch        *b_jetsJPTAK5_btag_secVertex;   //!
   TBranch        *b_jetsJPTAK5_chgEmE;   //!
   TBranch        *b_jetsJPTAK5_chgHadE;   //!
   TBranch        *b_jetsJPTAK5_chgMuE;   //!
   TBranch        *b_jetsJPTAK5_chg_Mult;   //!
   TBranch        *b_jetsJPTAK5_neutralEmE;   //!
   TBranch        *b_jetsJPTAK5_neutralHadE;   //!
   TBranch        *b_jetsJPTAK5_neutral_Mult;   //!
   TBranch        *b_jetsJPTAK5_mu_Mult;   //!
   TBranch        *b_jetsJPTAK5_emf;   //!
   TBranch        *b_jetsJPTAK5_ehf;   //!
   TBranch        *b_jetsJPTAK5_n60;   //!
   TBranch        *b_jetsJPTAK5_n90;   //!
   TBranch        *b_jetsJPTAK5_area;   //!
   TBranch        *b_jetsJPTAK5_mass;   //!
 
   TBranch        *b_Nmc_doc;   //!
   TBranch        *b_mc_doc_id;   //!
   TBranch        *b_mc_doc_pt;   //!
   TBranch        *b_mc_doc_px;   //!
   TBranch        *b_mc_doc_py;   //!
   TBranch        *b_mc_doc_pz;   //!
   TBranch        *b_mc_doc_eta;   //!
   TBranch        *b_mc_doc_phi;   //!
   TBranch        *b_mc_doc_theta;   //!
   TBranch        *b_mc_doc_energy;   //!
   TBranch        *b_mc_doc_status;   //!
   TBranch        *b_mc_doc_charge;   //!
   TBranch        *b_mc_doc_mother_id;   //!
   TBranch        *b_mc_doc_grandmother_id;   //!
   TBranch        *b_mc_doc_ggrandmother_id;   //!
   TBranch        *b_mc_doc_mother_pt;   //!
   TBranch        *b_mc_doc_vertex_x;   //!
   TBranch        *b_mc_doc_vertex_y;   //!
   TBranch        *b_mc_doc_vertex_z;   //!
   TBranch        *b_mc_doc_mass;   //!
   TBranch        *b_mc_doc_numOfDaughters;   //!
   TBranch        *b_mc_doc_numOfMothers;   //!
   TBranch        *b_Nmets;   //!
   TBranch        *b_mets_et;   //!
   TBranch        *b_mets_phi;   //!
   TBranch        *b_mets_ex;   //!
   TBranch        *b_mets_ey;   //!
   TBranch        *b_mets_gen_et;   //!
   TBranch        *b_mets_gen_phi;   //!
   TBranch        *b_mets_sign;   //!
   TBranch        *b_mets_sumEt;   //!
   TBranch        *b_mets_unCPhi;   //!
   TBranch        *b_mets_unCPt;   //!
   TBranch        *b_mets_et_muonCor;   //!
   TBranch        *b_mets_phi_muonCor;   //!
   TBranch        *b_mets_et_JESCor;   //!
   TBranch        *b_mets_phi_JESCor;   //!

   //3099
   TBranch        *b_NmetsKT4;   //!
   TBranch        *b_metsKT4_et;   //!
   TBranch        *b_metsKT4_phi;   //!
   TBranch        *b_metsKT4_ex;   //!
   TBranch        *b_metsKT4_ey;   //!
   TBranch        *b_metsKT4_sumEt;   //!
   TBranch        *b_metsKT4_et_JESCor;   //!
   TBranch        *b_metsKT4_phi_JESCor;   //!
   TBranch        *b_NmetsKT6;   //!
   TBranch        *b_metsKT6_et;   //!
   TBranch        *b_metsKT6_phi;   //!
   TBranch        *b_metsKT6_ex;   //!
   TBranch        *b_metsKT6_ey;   //!
   TBranch        *b_metsKT6_sumEt;   //!
   TBranch        *b_metsKT6_et_JESCor;   //!
   TBranch        *b_metsKT6_phi_JESCor;   //!
   TBranch        *b_NmetsSC5;   //!
   TBranch        *b_metsSC5_et;   //!
   TBranch        *b_metsSC5_phi;   //!
   TBranch        *b_metsSC5_ex;   //!
   TBranch        *b_metsSC5_ey;   //!
   TBranch        *b_metsSC5_sign;   //!
   TBranch        *b_metsSC5_sumEt;   //!
   TBranch        *b_metsSC5_et_JESCor;   //!
   TBranch        *b_metsSC5_phi_JESCor;   //!
   TBranch        *b_NmetsSC7;   //!
   TBranch        *b_metsSC7_et;   //!
   TBranch        *b_metsSC7_phi;   //!
   TBranch        *b_metsSC7_ex;   //!
   TBranch        *b_metsSC7_ey;   //!
   TBranch        *b_metsSC7_sign;   //!
   TBranch        *b_metsSC7_sumEt;   //!
   TBranch        *b_metsSC7_et_JESCor;   //!
   TBranch        *b_metsSC7_phi_JESCor;   //!
  
   TBranch        *b_Nmus;   //!
   TBranch        *b_mus_energy;   //!
   TBranch        *b_mus_et;   //!
   TBranch        *b_mus_eta;   //!
   TBranch        *b_mus_phi;   //!
   TBranch        *b_mus_pt;   //!
   TBranch        *b_mus_px;   //!
   TBranch        *b_mus_py;   //!
   TBranch        *b_mus_pz;   //!
   TBranch        *b_mus_status;   //!
   TBranch        *b_mus_theta;   //!
   TBranch        *b_mus_gen_id;   //!
   TBranch        *b_mus_gen_phi;   //!
   TBranch        *b_mus_gen_pt;   //!
   TBranch        *b_mus_gen_pz;   //!
   TBranch        *b_mus_gen_px;   //!
   TBranch        *b_mus_gen_py;   //!
   TBranch        *b_mus_gen_eta;   //!
   TBranch        *b_mus_gen_theta;   //!
   TBranch        *b_mus_gen_et;   //!
   TBranch        *b_mus_gen_mother_id;   //!
   TBranch        *b_mus_gen_mother_phi;   //!
   TBranch        *b_mus_gen_mother_pt;   //!
   TBranch        *b_mus_gen_mother_pz;   //!
   TBranch        *b_mus_gen_mother_px;   //!
   TBranch        *b_mus_gen_mother_py;   //!
   TBranch        *b_mus_gen_mother_eta;   //!
   TBranch        *b_mus_gen_mother_theta;   //!
   TBranch        *b_mus_gen_mother_et;   //!
   TBranch        *b_mus_tkHits;   //!
   TBranch        *b_mus_cIso;   //!
   TBranch        *b_mus_tIso;   //!
   TBranch        *b_mus_ecalIso;   //!
   TBranch        *b_mus_hcalIso;   //!
   TBranch        *b_mus_ecalvetoDep;   //!
   TBranch        *b_mus_hcalvetoDep;   //!
   TBranch        *b_mus_calEnergyEm;   //!
   TBranch        *b_mus_calEnergyHad;   //!
   TBranch        *b_mus_calEnergyHo;   //!
   TBranch        *b_mus_calEnergyEmS9;   //!
   TBranch        *b_mus_calEnergyHadS9;   //!
   TBranch        *b_mus_calEnergyHoS9;   //!
   TBranch        *b_mus_iso03_sumPt;   //!
   TBranch        *b_mus_iso03_emEt;   //!
   TBranch        *b_mus_iso03_hadEt;   //!
   TBranch        *b_mus_iso03_hoEt;   //!
   TBranch        *b_mus_iso03_nTracks;   //!
   TBranch        *b_mus_iso05_sumPt;   //!
   TBranch        *b_mus_iso05_emEt;   //!
   TBranch        *b_mus_iso05_hadEt;   //!
   TBranch        *b_mus_iso05_hoEt;   //!
   TBranch        *b_mus_iso05_nTracks;   //!
   TBranch        *b_mus_charge;   //!
   TBranch        *b_mus_cm_chi2;   //!
   TBranch        *b_mus_cm_ndof;   //!
   TBranch        *b_mus_cm_chg;   //!
   TBranch        *b_mus_cm_pt;   //!
   TBranch        *b_mus_cm_px;   //!
   TBranch        *b_mus_cm_py;   //!
   TBranch        *b_mus_cm_pz;   //!
   TBranch        *b_mus_cm_eta;   //!
   TBranch        *b_mus_cm_phi;   //!
   TBranch        *b_mus_cm_theta;   //!
   TBranch        *b_mus_cm_d0dum;   //!
   TBranch        *b_mus_cm_dz;   //!
   TBranch        *b_mus_cm_vx;   //!
   TBranch        *b_mus_cm_vy;   //!
   TBranch        *b_mus_cm_vz;   //!
   TBranch        *b_mus_cm_numvalhits;   //!
   TBranch        *b_mus_cm_numlosthits;   //!
   TBranch        *b_mus_cm_d0dumErr;   //!
   TBranch        *b_mus_cm_dzErr;   //!
   TBranch        *b_mus_cm_ptErr;   //!
   TBranch        *b_mus_cm_etaErr;   //!
   TBranch        *b_mus_cm_phiErr;   //!
   TBranch        *b_mus_tk_chi2;   //!
   TBranch        *b_mus_tk_ndof;   //!
   TBranch        *b_mus_tk_chg;   //!
   TBranch        *b_mus_tk_pt;   //!
   TBranch        *b_mus_tk_px;   //!
   TBranch        *b_mus_tk_py;   //!
   TBranch        *b_mus_tk_pz;   //!
   TBranch        *b_mus_tk_eta;   //!
   TBranch        *b_mus_tk_phi;   //!
   TBranch        *b_mus_tk_theta;   //!
   TBranch        *b_mus_tk_d0dum;   //!
   TBranch        *b_mus_tk_dz;   //!
   TBranch        *b_mus_tk_vx;   //!
   TBranch        *b_mus_tk_vy;   //!
   TBranch        *b_mus_tk_vz;   //!
   TBranch        *b_mus_tk_numvalhits;   //!
   TBranch        *b_mus_tk_numlosthits;   //!
   TBranch        *b_mus_tk_d0dumErr;   //!
   TBranch        *b_mus_tk_dzErr;   //!
   TBranch        *b_mus_tk_ptErr;   //!
   TBranch        *b_mus_tk_etaErr;   //!
   TBranch        *b_mus_tk_phiErr;   //!
   TBranch        *b_mus_stamu_chi2;   //!
   TBranch        *b_mus_stamu_ndof;   //!
   TBranch        *b_mus_stamu_chg;   //!
   TBranch        *b_mus_stamu_pt;   //!
   TBranch        *b_mus_stamu_px;   //!
   TBranch        *b_mus_stamu_py;   //!
   TBranch        *b_mus_stamu_pz;   //!
   TBranch        *b_mus_stamu_eta;   //!
   TBranch        *b_mus_stamu_phi;   //!
   TBranch        *b_mus_stamu_theta;   //!
   TBranch        *b_mus_stamu_d0dum;   //!
   TBranch        *b_mus_stamu_dz;   //!
   TBranch        *b_mus_stamu_vx;   //!
   TBranch        *b_mus_stamu_vy;   //!
   TBranch        *b_mus_stamu_vz;   //!
   TBranch        *b_mus_stamu_numvalhits;   //!
   TBranch        *b_mus_stamu_numlosthits;   //!
   TBranch        *b_mus_stamu_d0dumErr;   //!
   TBranch        *b_mus_stamu_dzErr;   //!
   TBranch        *b_mus_stamu_ptErr;   //!
   TBranch        *b_mus_stamu_etaErr;   //!
   TBranch        *b_mus_stamu_phiErr;   //!
   TBranch        *b_mus_num_matches;   //!
   TBranch        *b_mus_id_All;   //!
   TBranch        *b_mus_id_AllGlobalMuons;   //!
   TBranch        *b_mus_id_AllStandAloneMuons;   //!
   TBranch        *b_mus_id_AllTrackerMuons;   //!
   TBranch        *b_mus_id_TrackerMuonArbitrated;   //!
   TBranch        *b_mus_id_AllArbitrated;   //!
   TBranch        *b_mus_id_GlobalMuonPromptTight;   //!
   TBranch        *b_mus_id_TMLastStationLoose;   //!
   TBranch        *b_mus_id_TMLastStationTight;   //!
   TBranch        *b_mus_id_TM2DCompatibilityLoose;   //!
   TBranch        *b_mus_id_TM2DCompatibilityTight;   //!
   TBranch        *b_mus_id_TMOneStationLoose;   //!
   TBranch        *b_mus_id_TMOneStationTight;   //!
   TBranch        *b_mus_id_TMLastStationOptimizedLowPtLoose;   //!
   TBranch        *b_mus_id_TMLastStationOptimizedLowPtTight;   //!
   TBranch        *b_Nphotons;   //!
   TBranch        *b_photons_energy;   //!
   TBranch        *b_photons_et;   //!
   TBranch        *b_photons_eta;   //!
   TBranch        *b_photons_phi;   //!
   TBranch        *b_photons_pt;   //!
   TBranch        *b_photons_px;   //!
   TBranch        *b_photons_py;   //!
   TBranch        *b_photons_pz;   //!
   TBranch        *b_photons_status;   //!
   TBranch        *b_photons_theta;   //!
   TBranch        *b_photons_hadOverEM;   //!
   TBranch        *b_photons_scEnergy;   //!
   TBranch        *b_photons_scRawEnergy;   //!
   TBranch        *b_photons_scEta;   //!
   TBranch        *b_photons_scPhi;   //!
   TBranch        *b_photons_scEtaWidth;   //!
   TBranch        *b_photons_scPhiWidth;   //!
   TBranch        *b_photons_tIso;   //!
   TBranch        *b_photons_ecalIso;   //!
   TBranch        *b_photons_hcalIso;   //!
   TBranch        *b_photons_isoEcalRecHitDR04;   //!
   TBranch        *b_photons_isoHcalRecHitDR04;   //!
   TBranch        *b_photons_isoSolidTrkConeDR04;   //!
   TBranch        *b_photons_isoHollowTrkConeDR04;   //!
   TBranch        *b_photons_nTrkSolidConeDR04;   //!
   TBranch        *b_photons_nTrkHollowConeDR04;   //!
   TBranch        *b_photons_isoEcalRecHitDR03;   //!
   TBranch        *b_photons_isoHcalRecHitDR03;   //!
   TBranch        *b_photons_isoSolidTrkConeDR03;   //!
   TBranch        *b_photons_isoHollowTrkConeDR03;   //!
   TBranch        *b_photons_nTrkSolidConeDR03;   //!
   TBranch        *b_photons_nTrkHollowConeDR03;   //!
   TBranch        *b_photons_isAlsoElectron;   //!
   TBranch        *b_photons_hasPixelSeed;   //!
   TBranch        *b_photons_isConverted;   //!
   TBranch        *b_photons_isEBGap;   //!
   TBranch        *b_photons_isEEGap;   //!
   TBranch        *b_photons_isEBEEGap;   //!
   TBranch        *b_photons_isEBPho;   //!
   TBranch        *b_photons_isEEPho;   //!
   TBranch        *b_photons_isLoosePhoton;   //!
   TBranch        *b_photons_isTightPhoton;   //!
   TBranch        *b_photons_r9;   //!
   TBranch        *b_photons_gen_et;   //!
   TBranch        *b_photons_gen_eta;   //!
   TBranch        *b_photons_gen_phi;   //!
   TBranch        *b_photons_gen_id;   //!
   TBranch        *b_Npv;   //!
   TBranch        *b_pv_x;   //!
   TBranch        *b_pv_y;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_pv_xErr;   //!
   TBranch        *b_pv_yErr;   //!
   TBranch        *b_pv_zErr;   //!
   TBranch        *b_Ntcmets;   //!
   TBranch        *b_tcmets_et;   //!
   TBranch        *b_tcmets_phi;   //!
   TBranch        *b_tcmets_ex;   //!
   TBranch        *b_tcmets_ey;   //!
   TBranch        *b_tcmets_sumEt;   //!
   TBranch        *b_tcmets_et_muonCor;   //!
   TBranch        *b_tcmets_phi_muonCor;   //!
   TBranch        *b_tcmets_et_JESCor;   //!
   TBranch        *b_tcmets_phi_JESCor;   //!
   TBranch        *b_Ntracks;   //!
   TBranch        *b_tracks_chi2;   //!
   TBranch        *b_tracks_ndof;   //!
   TBranch        *b_tracks_chg;   //!
   TBranch        *b_tracks_pt;   //!
   TBranch        *b_tracks_px;   //!
   TBranch        *b_tracks_py;   //!
   TBranch        *b_tracks_pz;   //!
   TBranch        *b_tracks_eta;   //!
   TBranch        *b_tracks_phi;   //!
   TBranch        *b_tracks_theta;   //!
   TBranch        *b_tracks_d0dum;   //!
   TBranch        *b_tracks_dz;   //!
   TBranch        *b_tracks_vx;   //!
   TBranch        *b_tracks_vy;   //!
   TBranch        *b_tracks_vz;   //!
   TBranch        *b_tracks_numvalhits;   //!
   TBranch        *b_tracks_numlosthits;   //!
   TBranch        *b_tracks_d0dumErr;   //!
   TBranch        *b_tracks_dzErr;   //!
   TBranch        *b_tracks_ptErr;   //!
   TBranch        *b_tracks_etaErr;   //!
   TBranch        *b_tracks_phiErr;   //!
   TBranch        *b_tracks_Nrechits;   //!
   TBranch        *b_tracks_innerHitX;   //!
   TBranch        *b_tracks_innerHitY;   //!
   TBranch        *b_tracks_innerHitZ;   //!
   TBranch        *b_tracks_outerHitX;   //!
   TBranch        *b_tracks_outerHitY;   //!
   TBranch        *b_tracks_outerHitZ;   //!
   TBranch        *b_tracks_highPurity;   //!
   TBranch        *b_tracks_innerLayerMissingHits;  //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumiBlock;   //!


  //-----------------------
   // HLT electron trigger
  //-----------------------
  // For 8E29 (ref: https://twiki.cern.ch/twiki/bin/view/CMS/TSG_13_VIII_09_8E29)
  // For 1E31 (ref: https://twiki.cern.ch/twiki/bin/view/CMS/TSG_13_VIII_09_1E31)

  Double_t        HLT_Ele10_LW_EleId_L1R;
  Double_t        HLT_Ele10_SW_L1R;
  Double_t        HLT_Ele15_LW_L1R;
  Double_t        HLT_Ele15_SC10_LW_L1R;
  Double_t        HLT_Ele15_SW_L1R;
  Double_t        HLT_Ele15_SiStrip_L1R;
  Double_t        HLT_Ele20_LW_L1R;

  TBranch        *b_HLT_Ele10_LW_EleId_L1R;   //!
  TBranch        *b_HLT_Ele10_SW_L1R;   //!
  TBranch        *b_HLT_Ele15_LW_L1R;   //!
  TBranch        *b_HLT_Ele15_SC10_LW_L1R;   //!
  TBranch        *b_HLT_Ele15_SW_L1R;   //!
  TBranch        *b_HLT_Ele15_SiStrip_L1R;   //!
  TBranch        *b_HLT_Ele20_LW_L1R;   //!

  // PDF Weights
  //------------
  Double_t        PDFWcteq66_0;
  Double_t        PDFWcteq66_1;
  Double_t        PDFWcteq66_2;
  Double_t        PDFWcteq66_3;
  Double_t        PDFWcteq66_4;
  Double_t        PDFWcteq66_5;
  Double_t        PDFWcteq66_6;
  Double_t        PDFWcteq66_7;
  Double_t        PDFWcteq66_8;
  Double_t        PDFWcteq66_9;
  Double_t        PDFWcteq66_10;
  Double_t        PDFWcteq66_11;
  Double_t        PDFWcteq66_12;
  Double_t        PDFWcteq66_13;
  Double_t        PDFWcteq66_14;
  Double_t        PDFWcteq66_15;
  Double_t        PDFWcteq66_16;
  Double_t        PDFWcteq66_17;
  Double_t        PDFWcteq66_18;
  Double_t        PDFWcteq66_19;
  Double_t        PDFWcteq66_20;
  Double_t        PDFWcteq66_21;
  Double_t        PDFWcteq66_22;
  Double_t        PDFWcteq66_23;
  Double_t        PDFWcteq66_24;
  Double_t        PDFWcteq66_25;
  Double_t        PDFWcteq66_26;
  Double_t        PDFWcteq66_27;
  Double_t        PDFWcteq66_28;
  Double_t        PDFWcteq66_29;
  Double_t        PDFWcteq66_30;
  Double_t        PDFWcteq66_31;
  Double_t        PDFWcteq66_32;
  Double_t        PDFWcteq66_33;
  Double_t        PDFWcteq66_34;
  Double_t        PDFWcteq66_35;
  Double_t        PDFWcteq66_36;
  Double_t        PDFWcteq66_37;
  Double_t        PDFWcteq66_38;
  Double_t        PDFWcteq66_39;
  Double_t        PDFWcteq66_40;
  Double_t        PDFWcteq66_41;
  Double_t        PDFWcteq66_42;
  Double_t        PDFWcteq66_43;
  Double_t        PDFWcteq66_44;

  TBranch        *b_PDFWcteq66_0;   //!
  TBranch        *b_PDFWcteq66_1;   //!
  TBranch        *b_PDFWcteq66_2;   //!
  TBranch        *b_PDFWcteq66_3;   //!
  TBranch        *b_PDFWcteq66_4;   //!
  TBranch        *b_PDFWcteq66_5;   //!
  TBranch        *b_PDFWcteq66_6;   //!
  TBranch        *b_PDFWcteq66_7;   //!
  TBranch        *b_PDFWcteq66_8;   //!
  TBranch        *b_PDFWcteq66_9;   //!
  TBranch        *b_PDFWcteq66_10;   //!
  TBranch        *b_PDFWcteq66_11;   //!
  TBranch        *b_PDFWcteq66_12;   //!
  TBranch        *b_PDFWcteq66_13;   //!
  TBranch        *b_PDFWcteq66_14;   //!
  TBranch        *b_PDFWcteq66_15;   //!
  TBranch        *b_PDFWcteq66_16;   //!
  TBranch        *b_PDFWcteq66_17;   //!
  TBranch        *b_PDFWcteq66_18;   //!
  TBranch        *b_PDFWcteq66_19;   //!
  TBranch        *b_PDFWcteq66_20;   //!
  TBranch        *b_PDFWcteq66_21;   //!
  TBranch        *b_PDFWcteq66_22;   //!
  TBranch        *b_PDFWcteq66_23;   //!
  TBranch        *b_PDFWcteq66_24;   //!
  TBranch        *b_PDFWcteq66_25;   //!
  TBranch        *b_PDFWcteq66_26;   //!
  TBranch        *b_PDFWcteq66_27;   //!
  TBranch        *b_PDFWcteq66_28;   //!
  TBranch        *b_PDFWcteq66_29;   //!
  TBranch        *b_PDFWcteq66_30;   //!
  TBranch        *b_PDFWcteq66_31;   //!
  TBranch        *b_PDFWcteq66_32;   //!
  TBranch        *b_PDFWcteq66_33;   //!
  TBranch        *b_PDFWcteq66_34;   //!
  TBranch        *b_PDFWcteq66_35;   //!
  TBranch        *b_PDFWcteq66_36;   //!
  TBranch        *b_PDFWcteq66_37;   //!
  TBranch        *b_PDFWcteq66_38;   //!
  TBranch        *b_PDFWcteq66_39;   //!
  TBranch        *b_PDFWcteq66_40;   //!
  TBranch        *b_PDFWcteq66_41;   //!
  TBranch        *b_PDFWcteq66_42;   //!
  TBranch        *b_PDFWcteq66_43;   //!
  TBranch        *b_PDFWcteq66_44;   //!


  ana();
  ~ana(){};
    

  bool EventLoop();// the main analysis 

  //Methods to call from anascript.C
  void	  SetInputFile(const char* fname);
  void	  SetOutputFirstName(const string name);
  void	  SetOutputHistFile(const string name, const string mode="RECREATE");
  void	  SetOutputTextFile(const string name, const string mode="w");

  void	  SetGoodRuns(bool f) { keepgood = f; };   // use all runs or just good runs?
  void	  SetData(bool f)     { datafile = f; };   // is this a data file
  void    CheckTrigger(bool);  // check trigger fired?  time-consuming...
  void    CheckTrigger(bool, const string hlt);
  void    SetLimit(int n);     // testing first few events
  void    EstimateQCD();         // run QCD estimation
  bool    EstimateQCD(const string file);  // run QCD estimation
  void    EstimateWjets();       // run W+jets estimation
  bool    EstimateWjets(const string data, const string mc=""); // run W+jets estimation
  void    SetNtoyForM3Fit( int val ) { m_ntoy = val; } ; //number of toy exp to run for m3

  void    SetEleETcut(float);
  void    SetMuonPTcut(float);
  void    SetJetETcut(float);
  void    SetMETcut(float);
  void    SetHTcut(float);

  // electron ID
  enum  eID { robustTight, robustLoose, loose, tight, none };
  void    SetEleID(eID val) { m_eID = val;};

  void    SetAESHTcut(float cut)     { AES_HT_cut        = cut; };
  void    SetAESMETcut(float cut)    { AES_MET_cut       = cut; };
  void    SetAESZveto_TwoEle(bool f) { AES_useSimpleZveto = f; };

  void    ApplyMETcut(bool f)        { m_applyMETcut     = f; };
  void    RejectEndcapEle(bool f)    { m_rejectEndcapEle = f; };

  void    StudySystematics(const string&, const string&);
  void    StudyPDFunc(bool f)        { m_studyPDFunc     = f; };

  // Switches
  void Validation(bool val)            { m_doValidation      = val; };
  void ConversionStudySwitch(bool val) { m_ConversionStudies = val; };
  void StudyZveto(bool val)            { m_studyZveto        = val; };
  void PlotRelisoNES(bool val)         { m_plotRelisoNES     = val; };
  void SetDebug(bool val)              { m_debug             = val; };
  void SetJetAlgo(string val)          { m_jetAlgo           = val; };
  void SetMETAlgo(string val)          { m_metAlgo           = val; };
  void SetLHCEnergyInTeV(double val)   { m_LHCEnergyInTeV    = val; };
  void SetRunOnSD(bool val)            { m_runOnSD           = val; };
  void SetRunOnMyHLTskim(bool val)     { m_runOnMyHLTskim    = val; };


  //switch to preclue missing layers as these are not in 314 data samples
  void UseMissLayers(bool val)         { m_useMisslayers     = val; };

  void SignalIsAlpgen(bool val)        { signal_is_Alpgen    = val; };
  void SignalAlpgenThreshold(unsigned int val) { signal_Alpgen_matching_threshold = val; };

private:

  //-----------------
  // Intenal methods
  //-----------------
  bool	 KeepGoodRuns()       const { return keepgood; }
  bool	 IsData()             const { return datafile; };
  bool   GetTrigger()         const { return checkTrig; };
  int    GetLimit()           const { return numlimit; };
  eID    EleID()              const { return m_eID; };
  void   PrintCuts() const;   // print kinematic cuts
  string printEleID() const;
  void   CheckAvailableJetMET();
  void   Init(); //initialize tree branches
  void   ReadSelectedBranches() const;

  //geometry
  float calcDeltaR(const float phi1, const float eta1, const float phi2, const float eta2) const;
  float calcDeltaR(const TLorentzVector& p1,const TLorentzVector& p2) const;

  float calcDeltaEta(const TLorentzVector& p1,const TLorentzVector& p2) const;
  float calcDeltaPhi(const TLorentzVector& p1,const TLorentzVector& p2) const;

  bool  ConversionFinder(const TLorentzVector& e1, int mctype, int index_selected_ele);
  bool  ConversionFinder2(const TLorentzVector& e1, int mctype, int index_selected_ele);
  void  ConversionMCMatching(const TLorentzVector& e1, int mctype, bool isthisConversion);
  float MCTruthMatch(float eta, float phi);

  // QCD estimation
  pair<double,double> estimateQCD_computeFitResultUncertainty( const double est, TF1* f ) const;
  double estimateQCD_newEstimateByVaryingFitFunction( const double p1, const double p2, const double p3 ) const;
  pair<double,double> estimateQCD_assign_pos_neg_estimate( const double est, const double new_est[2] ) const;
  void Set_Reliso_bin_width(float bw) { m_QCDest_reliso_bin_width = bw; };
  float Get_Reliso_bin_width() const { return m_QCDest_reliso_bin_width; };

  // Histograms
  void BookHistograms();
  void BookHistograms_valid();
  void BookHistograms_basicKin();
  void BookHistograms_explore();
  void BookHistograms_nEle();
  void BookHistograms_eid();
  void BookHistograms_ed0();
  void BookHistograms_zVeto();
  void BookHistograms_met();
  void BookHistograms_HT();
  void BookHistograms_mtw();
  void BookHistograms_DRemu();
  void BookHistograms_DPhiEmet();
  void BookHistograms_DPhimetjet();
  void BookHistograms_conv();
  void BookHistograms_QCD();
  void BookHistograms_QCD_planA();
  void BookHistograms_QCD_planB();
  void BookHistograms_wj();
  void BookHistograms_event_table();
  void BookHistograms_btag();
  void BookHistograms_PDFunc();


  // Helper methods to make/fill histograms (vectors)
  typedef vector<vector<TH1*> > v2D_TH1; //TH1[][]
  typedef vector<vector<TH2*> > v2D_TH2; //TH2[][]
  typedef vector<v2D_TH1>       v3D_TH1; //TH1[][][]
  typedef vector<v2D_TH2>       v3D_TH2; //TH2[][][]

  void addHistoDataAndMC( vector<TH1*>& h, const string&, const string&, const int&, const float&, const float& );
  void addHistoDataAndMC( vector<TH2*>& h, const string&, const string&, const int&, const float&, const float&, 
  			  const int&, const float&, const float& );
  void fillHistoDataAndMC( const vector<TH1*>& h, const float& v ) const; //take out weight
  void fillHistoDataAndMC( const vector<TH2*>& h, const float& v1, const float& v2, const double& w ) const;

  void addHisto_Njet_DataAndMC( v2D_TH1& h, const string&, const string&, const int&, const float&, const float&);
  void fillHisto_Njet_DataAndMC( v2D_TH1& h, const float& v, const double& w ) const;
  void fillHistoNjet2D( v2D_TH1& h, const int& ec, const float& v, const double& w ) const;
  // validation
  void valid_mkHisto_cut_njet(v2D_TH1& h, const string&, const string&, const int&, const float&, const float& );
  void valid_fillHisto(v2D_TH1& h, const bool cuts[8], const double& value) const;
  // Reliso NES plots
  void iso_addHisto_nlevel_nj_nmc( v3D_TH1& h, const string&, const string&, const int&, const float&, const float& );
  void iso_addHisto_nlevel_nj_nmc( v3D_TH2& h, const string&, const string&, const int&, const float&, const float&,
  				   const int&, const float&, const float& );
  void iso_fillHisto_NES( const int& ilevel, const float& iso, const float& met, const bool& inBarrel ) const;
  void iso_fillHisto_nlevel_nj_nmc( const v2D_TH1& h, const float& ) const ;
  void iso_fillHisto_nlevel_nj_nmc( const v2D_TH2& h, const float&, const float&, const double& w ) const ;

  void fillHisto_event_tables();
  void fillHisto_PDF_weights( TH1F* h );


  // W+jets estimation
  void reco_hadronicTop_highestTopPT( const std::vector<TLorentzVector>&, const int nGoodIsoEle );
  pair<double,double> compute_M3(const std::vector<TLorentzVector>&) const;

  void SetHistoLabelCutNjet( TH2D *this_njetVcuts, const vector<string>& ve ) const;
  void SetHistoLabelEleID( const vector<TH1*>& h ) const;

  void   DefineCrossSection();
  void   DefineCrossSectionAlpgen7TeV();
  double GetCrossSection(const string) const;
  void   SetEventWeightMap();
  double GetWeight(const string) const;
  long   GetNinit(const string) const;

  // event-count tables
  void FillEventCounter(const int, const int&, const int&);
  void DrawEventPerNjetTable() const;
  void DrawSignalBGTable() const;
  void DrawMCTypeTable(    const string title ) const;
  void DrawQCDTable(       const string title ) const;
  void DrawSingleTopTable( const string title ) const;
  void DrawTTnjetTable(    const string title ) const;
  void DrawSignalAcceptanceTable( vector<string> ve) const;
  void printCutStage(const int&,const string&) const;
  void printCutStage(ofstream& os,const int&,const string&) const;

  // print event-count tables (with errors)
  void PrintErrorTables( vector<string> ve ) const;
  void PrintError_NjetVcut(ofstream&, const double [][5][24] ) const;

  double GetBayesUncertainty(int Ninitial) const;
  string ScrNum(double num) const;
  void printLine(ofstream &myfile, const double, const double) const;
  bool ScientificNotation;

  bool  is_mc_present(const int&) const;
  float compute_d0(const string&, const int&) const;
  float compute_mtw(const TVector2&, const TVector2&) const;

  void  PrintGenParticles() const;

  // conversion
  bool   DoConversionStudies()       { return m_ConversionStudies; };
  void PrintConversionTable();
  int ConversionCounter;
  int ConversionArray[23][2][6];
  void OptimiseConversionFinder(const TLorentzVector& e1, int mctype);  

//   TH2D *Conv_Opti[2];
//   TH2D *Conv_Optis[2];
//   TH2D *Conv_OptiL[2];
//   TH2D *Conv_Opti_extragran[2];
  vector<TH2D*> Conv_Opti;  //[2];
  vector<TH2D*> Conv_Optis; //[2];
  vector<TH2D*> Conv_OptiL; //[2];
  vector<TH2D*> Conv_Opti_extragran;//[2];
  TH1D *Conv_CheckDelR_GSFTk_ctfTk;
  int mycounter;

  // new (to run on mixed MC) (not currently needed)
  string CheckEventTypeFromMcTruth() const;

  float  getRelIso(const unsigned int&) const; 
  bool   passEleID(const unsigned int&) const;
  bool   passHLT() const;
  string printTimeNow() const;
  void   DoBTagging(const vector<TLorentzVector>&);
  bool   jetNotNearElectron(const TLorentzVector& j, const vector<TLorentzVector>& e) const;

  //--------------------
  // private variables
  //--------------------
  bool   datafile;
  bool   keepgood;
  bool   checkTrig;
  string HLTBit;
  int    numlimit;
  int    nfile;
  string outputHistFileName;
  string outputTextFileName;
  double this_weight;    // current event weight  
  int    m_nGoodJet;  //number of cleaned, good jets
  float  m_QCDest_reliso_bin_width;
  bool   m_doValidation;
  bool   m_studyZveto;
  bool   m_plotRelisoNES;
  bool   m_debug;
  bool   m_ConversionStudies;
  string m_jetAlgo;
  string m_metAlgo;
  double m_LHCEnergyInTeV;
  bool   m_runOnSD;
  bool   m_runOnMyHLTskim;
  bool   m_useMisslayers;
  short  m_muonCutNum;
  int    m_ntoy;
  int    m_nbtag_TCHE;
  int    m_nbtag_TCHP;
  int    m_nbtag_SSV;
  // pass flag
  bool   pass_met;
  // event counts
  //double e_plus_jet[nstage][ntjet][nmctype];
  vector<vector<vector<int> > >    e_plus_jet; //3D
  vector<vector<vector<double> > > e_plus_jet_weighted; //3D

  // cuts
  vector<string> ve; //list of cuts (excl nj)
  vector<string> ve2; //list of cuts (incl nj)
  float ELE_ETCUT;
  float MU_PTCUT;
  float JET_PTCUT;
  float METCUT;
  float HTCUT;
  int   nCutSetInScript;
  bool  m_applyMETcut;
  bool  m_rejectEndcapEle;
  eID   m_eID;
  float AES_HT_cut;
  float AES_MET_cut;
  bool  AES_useSimpleZveto;
  float intlumi;   // integrated luminosity assumed
  bool  useNewReliso;

  string doSystematics;
  string sysSample;
  bool   m_studyPDFunc;

  map<string,double> cross_section;
  map<string,long>   nInitialEventMC; //initial number of event, used to compute event weight
  map<string,double> weightMap;

  // Histograms
  //------------
  // validation
  v2D_TH1     valid_HT;         //[9][7]
  v2D_TH1     valid_jetsEt;     //[9][7]
  v2D_TH1     valid_jetsEta;    //[9][7]
  v2D_TH1     valid_jetsPhi;    //[9][7]
  v2D_TH1     valid_jets1stEt;  //[9][7]
  v2D_TH1     valid_jets2ndEt;  //[9][7]
  v2D_TH1     valid_jets3rdEt;  //[9][7]
  v2D_TH1     valid_jets4thEt;  //[9][7]
  v2D_TH1     valid_eleEt;      //[9][7]
  v2D_TH1     valid_eleEta;     //[9][7]
  v2D_TH1     valid_elePhi;     //[9][7]
  v2D_TH1     valid_eleCalIso;  //[9][7]
  v2D_TH1     valid_eleTrkIso;  //[9][7]
  v2D_TH1     valid_eleRelIso;  //[9][7]
  v2D_TH1     valid_eled0;      //[9][7]
  v2D_TH1     valid_metEt;      //[9][7]
  v2D_TH1     valid_metPhi;     //[9][7]
  v2D_TH1     valid_genTT_pt;   //[9][7]
  v2D_TH1     valid_genT_pt;    //[9][7]
  v2D_TH1     valid_recoM3;     //[9][7]
  v2D_TH1     valid_mass_ee;    //[9][7]
  v2D_TH1     valid_recoM3_PTMax; //[9][7]
  v2D_TH1     valid_numberTracks; //[9][7]
  v2D_TH1     valid_trackPt;      //[9][7]

  // basic
  //v2D_TH1F h_ele_ET2;//[4][nclass];
  // - ele
  vector<TH1*>  h_nele;    //[nclass];
  v2D_TH1       h_ele_ET;  //[4][nclass];
  v2D_TH1       h_ele_eta; //[4][nclass];
  v2D_TH1       h_ele_phi; //[4][nclass];
  v2D_TH1       h_ele_iso; //[4][nclass];
  // - jets
  vector<TH1*>  h_njet;    //[nclass]; //per MC type
  v2D_TH1       h_jet_PT;  //[5][nclass];
  v2D_TH1       h_jet_eta; //[5][nclass];
  v2D_TH1       h_jet_phi; //[5][nclass];
  // - met
  vector<TH1*>  h_metAlone;     //[nclass]; //per MC type
  vector<TH1*>  h_metAlone_phi; //[nclass];

  // explore
  vector<TH1*>  h_exp_ele_et;     //[nclass];  // selected ele et
  vector<TH1*>  h_exp_ele_eta;    //[nclass];  // selected ele eta
  vector<TH1*>  h_exp_j0_pt;      //[nclass];  // leading jet pt
  vector<TH1*>  h_exp_j1_pt;      //[nclass];  // 2n-leading jet pt
  vector<TH1*>  h_exp_DRej;       //[nclass];  // DR(e,j0)
  vector<TH1*>  h_exp_DPhiej;     //[nclass];  // DPhi(e,j0)
  vector<TH1*>  h_exp_DRjj;       //[nclass];  // DR(j0,j1)
  vector<TH1*>  h_exp_DPhijj;     //[nclass];  // DPhi(j0,j1)
  vector<TH2*>  h_exp_met_v_eeta; //[nclass];  // met:ele_eta

  // ele count
  v2D_TH1       h_nEle_all; //[7][nclass];
  v2D_TH1       h_nEle_s1;
  v2D_TH1       h_nEle_s2;
  v2D_TH1       h_nEle_s3_idLoose;
  v2D_TH1       h_nEle_s3_idTight;
  v2D_TH1       h_nEle_s3_idRL;
  v2D_TH1       h_nEle_s3_idRT;

  // ele id
  vector<TH1*>  h_eid; //[nclass];
  // ele d0
  vector<TH1*>  h_ed0_unCor;//[nclass]
  vector<TH1*>  h_ed0;      //[nclass] w.r.t beam spot
  vector<TH1*>  h_ed0_pass; //[nclass] w.r.t beam spot

  // z veto
  TH1F         *h_nGenBasicEle_Zee_allj;
  TH1F         *h_Zee_eta;
  TH1F         *h_Zee_pt;
  TH1F         *h_Z_photon_eta;
  TH1F         *h_Z_photon_et;
  TH1F         *h_Zee_photon_eta;
  TH1F         *h_Zee_photon_et;
  TH2F         *h_Zee_photon_eteta_2D;			
  TH1F         *h_Z_Nphotons;
  TH1F         *h_Zee_Nphotons;
  vector<TH1*>  h_mass_diele;            //[nclass];
  vector<TH1*>  h_mass_diele_new;        //[nclass];
  vector<TH1*>  h_mass_diele_lowMet_1j;  //[nclass];
  vector<TH1*>  h_mass_ephoton_lowMet_1j;//[nclass];
  vector<TH1*>  h_Nele_lowMet_1j;        //[nclass];
  vector<TH1*>  h_Nphoton_lowMet_1j;     //[nclass];
  vector<TH1*>  h_photon_eta_lowMet_1j;  //[nclass];
  vector<TH1*>  h_photon_et_lowMet_1j;   //[nclass];
  vector<TH1*>  h_photon1_eta_lowMet_1j; //[nclass];
  vector<TH1*>  h_photon1_et_lowMet_1j;  //[nclass];

  // MET
  v2D_TH1       h_met;                  //[7][nclass]; //user-chosen MET
  v2D_TH1       h_met_mu;               //[7][nclass]; //muon-MET
  v2D_TH1       h_met_t1;               //[7][nclass]; //type1-MET
  v2D_TH1       h_met_BA;               //[7][nclass]; //user-chosen MET (Barrel)
  v2D_TH1       h_met_mu_BA;            //[7][nclass]; //muon-MET (Barrel)
  v2D_TH1       h_met_t1_BA;            //[7][nclass]; //type1-MET (Barrel)
  vector<TH1*>  h_met_ante_ISO;         //[nclass]; //user-chosen MET
  vector<TH1*>  h_met_ante_ISO_mu;      //[nclass];
  vector<TH1*>  h_met_ante_ISO_t1;      //[nclass];
  v2D_TH1       h_met_gen;              //[7][nclass]; //BA+EN
  v2D_TH1       h_met_gen_diff_t1;      //[7][nclass]; //delta et (t1-gen)
  v2D_TH1       h_met_gen_diff_mu;      //[7][nclass]; //delta et (mu-gen)
  v2D_TH1       h_met_gen_dphi_t1;      //[7][nclass]; //delta phi (t1-gen)
  v2D_TH1       h_met_gen_dphi_mu;      //[7][nclass]; //delta phi (mu-gen)
  v2D_TH1       h_met_gen_BA;           //[7][nclass]; //Barrel
  v2D_TH1       h_met_gen_diff_t1_BA;   //[7][nclass];
  v2D_TH1       h_met_gen_diff_mu_BA;   //[7][nclass];
  v2D_TH1       h_met_gen_dphi_t1_BA;   //[7][nclass];
  v2D_TH1       h_met_gen_dphi_mu_BA;   //[7][nclass];

  // HT
  v2D_TH1       h_HT;     //[7][nclass];

  // mtw
  vector<TH1*>  h_mtw_mu_incl; //[nclass];   //inclusive
  vector<TH1*>  h_mtw_t1_incl; //[nclass];   //type 1 calomet
  v2D_TH1       h_mtw_mu;      //[7][nclass];  // after all but MET cut
  v2D_TH1       h_mtw_t1;      //[7][nclass];

  // DR(ele,mu)
  vector<TH1*>  h_DRemu_selE_GoodMu;      // [nclass];
  vector<TH1*>  h_DRemu_selE_GoodMu_pass; // [nclass];

  // Dphi(ele,met)
  vector<TH1*>  h_DPhiEmet_mu_incl;   //[nclass];
  vector<TH1*>  h_DPhiEmet_t1_incl;   //[class];
  v2D_TH1       h_DPhiEmet_mu;        //[7][nclass];
  v2D_TH1       h_DPhiEmet_t1;        //[7][nclass];

  // Dphi(met,jet)
  v2D_TH1       h_DPhiMetJet_mu_goodE;  //[7][nclass]; //pass HLT, e30, eta<2.5
  v2D_TH1       h_DPhiMetJet_t1_goodE;  //[7][nclass];
  v2D_TH1       h_DPhiMetJet_gen_goodE; //[7][nclass]; //gen
  v2D_TH1       h_DPhiMetJet_mu;        //[7][nclass]; //pass all cut except met
  v2D_TH1       h_DPhiMetJet_t1;        //[7][nclass];
  v2D_TH1       h_DPhiMetJet_gen;       //[7][nclass];

  // QCD estimation
  v2D_TH1       h_QCDest_CombRelIso;    //[7][nclass];     //"new" formulation (0-infinity)
  // dir /AES/
  v2D_TH1       h_QCDest_CombRelIso_AES;    //[7][nclass];  //ALL AES
  v2D_TH1       h_QCDest_CombRelIso_AES_minusMET;    //[7][nclass];  // AES except MET cut
  v2D_TH1       h_QCDest_CombRelIso_AES_minusHT;    //[7][nclass];   // AES except HT cut
  v2D_TH1       h_QCDest_CombRelIso_AES_minusTighterZ;    //[7][nclass];  // AES except tighter Z veto cut (both mee and mep)
  v2D_TH1       h_QCDest_CombRelIso_AES_before;    //[7][nclass];      // Before AES
  v2D_TH1       h_QCDest_CombRelIso_AES_justMET;    //[7][nclass];     // AES: just MET (<x)
  v2D_TH1       h_QCDest_CombRelIso_AES_justHighMET;    //[7][nclass]; // AES: just high MET (>x)
  v2D_TH1       h_QCDest_CombRelIso_AES_justHT;    //[7][nclass];      // AES: just HT
  v2D_TH1       h_QCDest_CombRelIso_AES_justZ;    //[7][nclass];       // AES: just Tighter Z (both)
  // dir /NES/
  v3D_TH2       h_QCDest_isoVmet_NES;     //[nLevel][7][nmc] weighted
  v3D_TH2       h_QCDest_isoVmet_NES_barrel;
  v3D_TH2       h_QCDest_isoVmet_NES_endcap;
  v3D_TH2       h_QCDest_isoVmet_NES_uw;  //[nLevel][7][nmc] unweighted
  v3D_TH2       h_QCDest_isoVmet_NES_uw_barrel;
  v3D_TH2       h_QCDest_isoVmet_NES_uw_endcap;  

  // NB: actually the following plots can be derived from the isoVmet scatter plot using projection,
  //     but right now keeping them so that we can just plot it without extra macro.
  // a) reliso: no MET cut
  v3D_TH1       h_QCDest_iso_NES; //[nLevel=11][nj=7][nmc=16] weighted
  v3D_TH1       h_QCDest_iso_NES_barrel;
  v3D_TH1       h_QCDest_iso_NES_endcap;
  // b) reliso: fail MET cut
  v3D_TH1       h_QCDest_iso_NES_loMET;
  v3D_TH1       h_QCDest_iso_NES_loMET_barrel;
  v3D_TH1       h_QCDest_iso_NES_loMET_endcap;
  // c) reliso: pass MET cut
  v3D_TH1       h_QCDest_iso_NES_hiMET;
  v3D_TH1       h_QCDest_iso_NES_hiMET_barrel;
  v3D_TH1       h_QCDest_iso_NES_hiMET_endcap;

  // AES plan A
  v2D_TH1       h_QCDest_CombRelIso_AES_planA1_e20; //[7][nclass]//no !conv cut
  v2D_TH1       h_QCDest_CombRelIso_AES_planA1_e30; //[7][nclass]
  v2D_TH1       h_QCDest_CombRelIso_AES_planA2_e20; //[7][nclass]//invert !conv cut
  v2D_TH1       h_QCDest_CombRelIso_AES_planA2_e30; //[7][nclass]
  v2D_TH1       h_QCDest_CombRelIso_AES_planA3_e20; //[7][nclass]//invert d0 cut
  v2D_TH1       h_QCDest_CombRelIso_AES_planA3_e30; //[7][nclass]

  // AES plan B
  v2D_TH1       h_QCDest_CombRelIso_AES_planB1_e20; //[7][nclass]//EleET >x
  v2D_TH1       h_QCDest_CombRelIso_AES_planB1_e30; //[7][nclass]
  v2D_TH1       h_QCDest_CombRelIso_AES_planB2_e20; //[7][nclass]//EleET < x
  v2D_TH1       h_QCDest_CombRelIso_AES_planB2_e30; //[7][nclass]
  v2D_TH1       h_QCDest_CombRelIso_AES_planB3_e20; //[7][nclass]//fail RT ID
  v2D_TH1       h_QCDest_CombRelIso_AES_planB3_e30; //[7][nclass]
  v2D_TH1       h_QCDest_CombRelIso_AES_planB3b_e20; //[7][nclass]//fail RT ID (BARREL)
  v2D_TH1       h_QCDest_CombRelIso_AES_planB3b_e30; //[7][nclass]
  v2D_TH1       h_QCDest_CombRelIso_AES_planB4_e20; //[7][nclass]//fail RL ID
  v2D_TH1       h_QCDest_CombRelIso_AES_planB4_e30; //[7][nclass]
  v2D_TH1       h_QCDest_CombRelIso_AES_planB5_e20; //[7][nclass]//fail Loose ID
  v2D_TH1       h_QCDest_CombRelIso_AES_planB5_e30; //[7][nclass]
  v2D_TH1       h_QCDest_CombRelIso_AES_planB6_e20; //[7][nclass]//fail Tight ID
  v2D_TH1       h_QCDest_CombRelIso_AES_planB6_e30; //[7][nclass]
  v2D_TH1       h_QCDest_CombRelIso_AES_planB7_e20; //[7][nclass]//d0 > 200um, pass RT
  v2D_TH1       h_QCDest_CombRelIso_AES_planB7_e30; //[7][nclass]
  v2D_TH1       h_QCDest_CombRelIso_AES_planB8_e20; //[7][nclass]//d0 > 200um, fail RT
  v2D_TH1       h_QCDest_CombRelIso_AES_planB8_e30; //[7][nclass]


  // Wjet estimation
  TH1D *h_hadTop_maxPT_mass_4j;   // 960 bins (0-960)
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

  // event table
  TH2D *Signal_njetsVcuts;
  TH2D *QCD_njetsVcuts;
  TH2D *Wjets_njetsVcuts;
  TH2D *Zjets_njetsVcuts;
  TH2D *VQQ_njetsVcuts; 
  TH2D *SingleTop_njetsVcuts;
  TH2D *Data_njetsVcuts;

  // btag
  vector<TH1*> h_nbtag_TCHE; //[nclass]
  vector<TH1*> h_nbtag_TCHP; //[nclass]
  vector<TH1*> h_nbtag_SSV;  //[nclass]

  // signal PDF unc
  TH1F *h_pdf_total;
  TH1F *h_pdf_pass;
  TH1F *h_pdf_eff;



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

  // z studies
  bool isZee, isZmm, isZtt;
  vector<TLorentzVector> Zele;

  bool signal_is_Alpgen;
  unsigned int signal_Alpgen_matching_threshold;

};


void ana::Init(){
   cout << "\n Initializing tree branches\n"  << endl;

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

   jetsKT6_energy = 0;
   jetsKT6_et = 0;
   jetsKT6_eta = 0;
   jetsKT6_phi = 0;
   jetsKT6_pt = 0;
   jetsKT6_px = 0;
   jetsKT6_py = 0;
   jetsKT6_pz = 0;
   jetsKT6_status = 0;
   jetsKT6_theta = 0;
   jetsKT6_btag_TC_highPur = 0;
   jetsKT6_btag_TC_highEff = 0;
   jetsKT6_btag_jetProb = 0;
   jetsKT6_btag_jetBProb = 0;
   jetsKT6_btag_softEle = 0;
   jetsKT6_btag_softMuon = 0;
   jetsKT6_btag_softMuonNoIP = 0;
   jetsKT6_btag_secVertex = 0;
   jetsKT6_chgEmE = 0;
   jetsKT6_chgHadE = 0;
   jetsKT6_chgMuE = 0;
   jetsKT6_chg_Mult = 0;
   jetsKT6_neutralEmE = 0;
   jetsKT6_neutralHadE = 0;
   jetsKT6_neutral_Mult = 0;
   jetsKT6_mu_Mult = 0;
   jetsKT6_emf = 0;
   jetsKT6_ehf = 0;
   jetsKT6_n60 = 0;
   jetsKT6_n90 = 0;
   jetsKT6_area = 0;
   jetsKT6_mass = 0;

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

   jetsSC7_energy = 0;
   jetsSC7_et = 0;
   jetsSC7_eta = 0;
   jetsSC7_phi = 0;
   jetsSC7_pt = 0;
   jetsSC7_px = 0;
   jetsSC7_py = 0;
   jetsSC7_pz = 0;
   jetsSC7_status = 0;
   jetsSC7_theta = 0;
   jetsSC7_btag_TC_highPur = 0;
   jetsSC7_btag_TC_highEff = 0;
   jetsSC7_btag_jetProb = 0;
   jetsSC7_btag_jetBProb = 0;
   jetsSC7_btag_softEle = 0;
   jetsSC7_btag_softMuon = 0;
   jetsSC7_btag_softMuonNoIP = 0;
   jetsSC7_btag_secVertex = 0;
   jetsSC7_chgEmE = 0;
   jetsSC7_chgHadE = 0;
   jetsSC7_chgMuE = 0;
   jetsSC7_chg_Mult = 0;
   jetsSC7_neutralEmE = 0;
   jetsSC7_neutralHadE = 0;
   jetsSC7_neutral_Mult = 0;
   jetsSC7_mu_Mult = 0;
   jetsSC7_emf = 0;
   jetsSC7_ehf = 0;
   jetsSC7_n60 = 0;
   jetsSC7_n90 = 0;
   jetsSC7_area = 0;
   jetsSC7_mass = 0;

   jetsJPTAK5_energy = 0;
   jetsJPTAK5_et = 0;
   jetsJPTAK5_eta = 0;
   jetsJPTAK5_phi = 0;
   jetsJPTAK5_pt = 0;
   jetsJPTAK5_px = 0;
   jetsJPTAK5_py = 0;
   jetsJPTAK5_pz = 0;
   jetsJPTAK5_status = 0;
   jetsJPTAK5_theta = 0;
   jetsJPTAK5_btag_TC_highPur = 0;
   jetsJPTAK5_btag_TC_highEff = 0;
   jetsJPTAK5_btag_jetProb = 0;
   jetsJPTAK5_btag_jetBProb = 0;
   jetsJPTAK5_btag_softEle = 0;
   jetsJPTAK5_btag_softMuon = 0;
   jetsJPTAK5_btag_softMuonNoIP = 0;
   jetsJPTAK5_btag_secVertex = 0;
   jetsJPTAK5_chgEmE = 0;
   jetsJPTAK5_chgHadE = 0;
   jetsJPTAK5_chgMuE = 0;
   jetsJPTAK5_chg_Mult = 0;
   jetsJPTAK5_neutralEmE = 0;
   jetsJPTAK5_neutralHadE = 0;
   jetsJPTAK5_neutral_Mult = 0;
   jetsJPTAK5_mu_Mult = 0;
   jetsJPTAK5_emf = 0;
   jetsJPTAK5_ehf = 0;
   jetsJPTAK5_n60 = 0;
   jetsJPTAK5_n90 = 0;
   jetsJPTAK5_area = 0;
   jetsJPTAK5_mass = 0;

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
   metsKT6_et = 0;
   metsKT6_phi = 0;
   metsKT6_ex = 0;
   metsKT6_ey = 0;
   metsKT6_sumEt = 0;
   metsKT6_et_JESCor = 0;
   metsKT6_phi_JESCor = 0;
   metsSC5_et = 0;
   metsSC5_phi = 0;
   metsSC5_ex = 0;
   metsSC5_ey = 0;
   metsSC5_sumEt = 0;
   metsSC5_et_JESCor = 0;
   metsSC5_phi_JESCor = 0;
   metsSC7_et = 0;
   metsSC7_phi = 0;
   metsSC7_ex = 0;
   metsSC7_ey = 0;
   metsSC7_sumEt = 0;
   metsSC7_et_JESCor = 0;
   metsSC7_phi_JESCor = 0;
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
   tcmets_et_JESCor = 0;
   tcmets_phi_JESCor = 0;
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
   if( m_jetAlgo=="pfjet" ){
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
   if( m_metAlgo=="pfmet" ){
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
   if(m_useMisslayers) chain->SetBranchAddress("els_innerLayerMissingHits", &els_innerLayerMissingHits, &b_els_innerLayerMissingHits);

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

   if( m_jetAlgo=="KT4" ) {     
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
   else if( m_jetAlgo=="KT6" ) {
     chain->SetBranchAddress("jetsKT6_energy", &jetsKT6_energy, &b_jetsKT6_energy);
     chain->SetBranchAddress("jetsKT6_et", &jetsKT6_et, &b_jetsKT6_et);
     chain->SetBranchAddress("jetsKT6_eta", &jetsKT6_eta, &b_jetsKT6_eta);
     chain->SetBranchAddress("jetsKT6_phi", &jetsKT6_phi, &b_jetsKT6_phi);
     chain->SetBranchAddress("jetsKT6_pt", &jetsKT6_pt, &b_jetsKT6_pt);
     chain->SetBranchAddress("jetsKT6_px", &jetsKT6_px, &b_jetsKT6_px);
     chain->SetBranchAddress("jetsKT6_py", &jetsKT6_py, &b_jetsKT6_py);
     chain->SetBranchAddress("jetsKT6_pz", &jetsKT6_pz, &b_jetsKT6_pz);
     chain->SetBranchAddress("jetsKT6_status", &jetsKT6_status, &b_jetsKT6_status);
     chain->SetBranchAddress("jetsKT6_theta", &jetsKT6_theta, &b_jetsKT6_theta);
     chain->SetBranchAddress("jetsKT6_btag_TC_highPur", &jetsKT6_btag_TC_highPur, &b_jetsKT6_btag_TC_highPur);
     chain->SetBranchAddress("jetsKT6_btag_TC_highEff", &jetsKT6_btag_TC_highEff, &b_jetsKT6_btag_TC_highEff);
     chain->SetBranchAddress("jetsKT6_btag_jetProb", &jetsKT6_btag_jetProb, &b_jetsKT6_btag_jetProb);
     chain->SetBranchAddress("jetsKT6_btag_jetBProb", &jetsKT6_btag_jetBProb, &b_jetsKT6_btag_jetBProb);
     chain->SetBranchAddress("jetsKT6_btag_softEle", &jetsKT6_btag_softEle, &b_jetsKT6_btag_softEle);
     chain->SetBranchAddress("jetsKT6_btag_softMuon", &jetsKT6_btag_softMuon, &b_jetsKT6_btag_softMuon);
     chain->SetBranchAddress("jetsKT6_btag_softMuonNoIP", &jetsKT6_btag_softMuonNoIP, &b_jetsKT6_btag_softMuonNoIP);
     chain->SetBranchAddress("jetsKT6_btag_secVertex", &jetsKT6_btag_secVertex, &b_jetsKT6_btag_secVertex);
     chain->SetBranchAddress("jetsKT6_chgEmE", &jetsKT6_chgEmE, &b_jetsKT6_chgEmE);
     chain->SetBranchAddress("jetsKT6_chgHadE", &jetsKT6_chgHadE, &b_jetsKT6_chgHadE);
     chain->SetBranchAddress("jetsKT6_chgMuE", &jetsKT6_chgMuE, &b_jetsKT6_chgMuE);
     chain->SetBranchAddress("jetsKT6_chg_Mult", &jetsKT6_chg_Mult, &b_jetsKT6_chg_Mult);
     chain->SetBranchAddress("jetsKT6_neutralEmE", &jetsKT6_neutralEmE, &b_jetsKT6_neutralEmE);
     chain->SetBranchAddress("jetsKT6_neutralHadE", &jetsKT6_neutralHadE, &b_jetsKT6_neutralHadE);
     chain->SetBranchAddress("jetsKT6_neutral_Mult", &jetsKT6_neutral_Mult, &b_jetsKT6_neutral_Mult);
     chain->SetBranchAddress("jetsKT6_mu_Mult", &jetsKT6_mu_Mult, &b_jetsKT6_mu_Mult);
     chain->SetBranchAddress("jetsKT6_emf", &jetsKT6_emf, &b_jetsKT6_emf);
     chain->SetBranchAddress("jetsKT6_ehf", &jetsKT6_ehf, &b_jetsKT6_ehf);
     chain->SetBranchAddress("jetsKT6_n60", &jetsKT6_n60, &b_jetsKT6_n60);
     chain->SetBranchAddress("jetsKT6_n90", &jetsKT6_n90, &b_jetsKT6_n90);
     chain->SetBranchAddress("jetsKT6_area", &jetsKT6_area, &b_jetsKT6_area);
     chain->SetBranchAddress("jetsKT6_mass", &jetsKT6_mass, &b_jetsKT6_mass);
   }
   else if( m_jetAlgo=="SC5" ) {     
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
   else if( m_jetAlgo=="SC7" ) {     
     chain->SetBranchAddress("NjetsSC7", &NjetsSC7, &b_NjetsSC7);
     chain->SetBranchAddress("jetsSC7_energy", &jetsSC7_energy, &b_jetsSC7_energy);
     chain->SetBranchAddress("jetsSC7_et", &jetsSC7_et, &b_jetsSC7_et);
     chain->SetBranchAddress("jetsSC7_eta", &jetsSC7_eta, &b_jetsSC7_eta);
     chain->SetBranchAddress("jetsSC7_phi", &jetsSC7_phi, &b_jetsSC7_phi);
     chain->SetBranchAddress("jetsSC7_pt", &jetsSC7_pt, &b_jetsSC7_pt);
     chain->SetBranchAddress("jetsSC7_px", &jetsSC7_px, &b_jetsSC7_px);
     chain->SetBranchAddress("jetsSC7_py", &jetsSC7_py, &b_jetsSC7_py);
     chain->SetBranchAddress("jetsSC7_pz", &jetsSC7_pz, &b_jetsSC7_pz);
     chain->SetBranchAddress("jetsSC7_status", &jetsSC7_status, &b_jetsSC7_status);
     chain->SetBranchAddress("jetsSC7_theta", &jetsSC7_theta, &b_jetsSC7_theta);
     chain->SetBranchAddress("jetsSC7_btag_TC_highPur", &jetsSC7_btag_TC_highPur, &b_jetsSC7_btag_TC_highPur);
     chain->SetBranchAddress("jetsSC7_btag_TC_highEff", &jetsSC7_btag_TC_highEff, &b_jetsSC7_btag_TC_highEff);
     chain->SetBranchAddress("jetsSC7_btag_jetProb", &jetsSC7_btag_jetProb, &b_jetsSC7_btag_jetProb);
     chain->SetBranchAddress("jetsSC7_btag_jetBProb", &jetsSC7_btag_jetBProb, &b_jetsSC7_btag_jetBProb);
     chain->SetBranchAddress("jetsSC7_btag_softEle", &jetsSC7_btag_softEle, &b_jetsSC7_btag_softEle);
     chain->SetBranchAddress("jetsSC7_btag_softMuon", &jetsSC7_btag_softMuon, &b_jetsSC7_btag_softMuon);
     chain->SetBranchAddress("jetsSC7_btag_softMuonNoIP", &jetsSC7_btag_softMuonNoIP, &b_jetsSC7_btag_softMuonNoIP);
     chain->SetBranchAddress("jetsSC7_btag_secVertex", &jetsSC7_btag_secVertex, &b_jetsSC7_btag_secVertex);
     chain->SetBranchAddress("jetsSC7_chgEmE", &jetsSC7_chgEmE, &b_jetsSC7_chgEmE);
     chain->SetBranchAddress("jetsSC7_chgHadE", &jetsSC7_chgHadE, &b_jetsSC7_chgHadE);
     chain->SetBranchAddress("jetsSC7_chgMuE", &jetsSC7_chgMuE, &b_jetsSC7_chgMuE);
     chain->SetBranchAddress("jetsSC7_chg_Mult", &jetsSC7_chg_Mult, &b_jetsSC7_chg_Mult);
     chain->SetBranchAddress("jetsSC7_neutralEmE", &jetsSC7_neutralEmE, &b_jetsSC7_neutralEmE);
     chain->SetBranchAddress("jetsSC7_neutralHadE", &jetsSC7_neutralHadE, &b_jetsSC7_neutralHadE);
     chain->SetBranchAddress("jetsSC7_neutral_Mult", &jetsSC7_neutral_Mult, &b_jetsSC7_neutral_Mult);
     chain->SetBranchAddress("jetsSC7_mu_Mult", &jetsSC7_mu_Mult, &b_jetsSC7_mu_Mult);
     chain->SetBranchAddress("jetsSC7_emf", &jetsSC7_emf, &b_jetsSC7_emf);
     chain->SetBranchAddress("jetsSC7_ehf", &jetsSC7_ehf, &b_jetsSC7_ehf);
     chain->SetBranchAddress("jetsSC7_n60", &jetsSC7_n60, &b_jetsSC7_n60);
     chain->SetBranchAddress("jetsSC7_n90", &jetsSC7_n90, &b_jetsSC7_n90);
     chain->SetBranchAddress("jetsSC7_area", &jetsSC7_area, &b_jetsSC7_area);
     chain->SetBranchAddress("jetsSC7_mass", &jetsSC7_mass, &b_jetsSC7_mass);
   }
   else if( m_jetAlgo=="JPTAK5" ) {     
     chain->SetBranchAddress("NjetsJPTAK5", &NjetsJPTAK5, &b_NjetsJPTAK5);
     chain->SetBranchAddress("jetsJPTAK5_energy", &jetsJPTAK5_energy, &b_jetsJPTAK5_energy);
     chain->SetBranchAddress("jetsJPTAK5_et", &jetsJPTAK5_et, &b_jetsJPTAK5_et);
     chain->SetBranchAddress("jetsJPTAK5_eta", &jetsJPTAK5_eta, &b_jetsJPTAK5_eta);
     chain->SetBranchAddress("jetsJPTAK5_phi", &jetsJPTAK5_phi, &b_jetsJPTAK5_phi);
     chain->SetBranchAddress("jetsJPTAK5_pt", &jetsJPTAK5_pt, &b_jetsJPTAK5_pt);
     chain->SetBranchAddress("jetsJPTAK5_px", &jetsJPTAK5_px, &b_jetsJPTAK5_px);
     chain->SetBranchAddress("jetsJPTAK5_py", &jetsJPTAK5_py, &b_jetsJPTAK5_py);
     chain->SetBranchAddress("jetsJPTAK5_pz", &jetsJPTAK5_pz, &b_jetsJPTAK5_pz);
     chain->SetBranchAddress("jetsJPTAK5_status", &jetsJPTAK5_status, &b_jetsJPTAK5_status);
     chain->SetBranchAddress("jetsJPTAK5_theta", &jetsJPTAK5_theta, &b_jetsJPTAK5_theta);
     chain->SetBranchAddress("jetsJPTAK5_btag_TC_highPur", &jetsJPTAK5_btag_TC_highPur, &b_jetsJPTAK5_btag_TC_highPur);
     chain->SetBranchAddress("jetsJPTAK5_btag_TC_highEff", &jetsJPTAK5_btag_TC_highEff, &b_jetsJPTAK5_btag_TC_highEff);
     chain->SetBranchAddress("jetsJPTAK5_btag_jetProb", &jetsJPTAK5_btag_jetProb, &b_jetsJPTAK5_btag_jetProb);
     chain->SetBranchAddress("jetsJPTAK5_btag_jetBProb", &jetsJPTAK5_btag_jetBProb, &b_jetsJPTAK5_btag_jetBProb);
     chain->SetBranchAddress("jetsJPTAK5_btag_softEle", &jetsJPTAK5_btag_softEle, &b_jetsJPTAK5_btag_softEle);
     chain->SetBranchAddress("jetsJPTAK5_btag_softMuon", &jetsJPTAK5_btag_softMuon, &b_jetsJPTAK5_btag_softMuon);
     chain->SetBranchAddress("jetsJPTAK5_btag_softMuonNoIP", &jetsJPTAK5_btag_softMuonNoIP, &b_jetsJPTAK5_btag_softMuonNoIP);
     chain->SetBranchAddress("jetsJPTAK5_btag_secVertex", &jetsJPTAK5_btag_secVertex, &b_jetsJPTAK5_btag_secVertex);
     chain->SetBranchAddress("jetsJPTAK5_chgEmE", &jetsJPTAK5_chgEmE, &b_jetsJPTAK5_chgEmE);
     chain->SetBranchAddress("jetsJPTAK5_chgHadE", &jetsJPTAK5_chgHadE, &b_jetsJPTAK5_chgHadE);
     chain->SetBranchAddress("jetsJPTAK5_chgMuE", &jetsJPTAK5_chgMuE, &b_jetsJPTAK5_chgMuE);
     chain->SetBranchAddress("jetsJPTAK5_chg_Mult", &jetsJPTAK5_chg_Mult, &b_jetsJPTAK5_chg_Mult);
     chain->SetBranchAddress("jetsJPTAK5_neutralEmE", &jetsJPTAK5_neutralEmE, &b_jetsJPTAK5_neutralEmE);
     chain->SetBranchAddress("jetsJPTAK5_neutralHadE", &jetsJPTAK5_neutralHadE, &b_jetsJPTAK5_neutralHadE);
     chain->SetBranchAddress("jetsJPTAK5_neutral_Mult", &jetsJPTAK5_neutral_Mult, &b_jetsJPTAK5_neutral_Mult);
     chain->SetBranchAddress("jetsJPTAK5_mu_Mult", &jetsJPTAK5_mu_Mult, &b_jetsJPTAK5_mu_Mult);
     chain->SetBranchAddress("jetsJPTAK5_emf", &jetsJPTAK5_emf, &b_jetsJPTAK5_emf);
     chain->SetBranchAddress("jetsJPTAK5_ehf", &jetsJPTAK5_ehf, &b_jetsJPTAK5_ehf);
     chain->SetBranchAddress("jetsJPTAK5_n60", &jetsJPTAK5_n60, &b_jetsJPTAK5_n60);
     chain->SetBranchAddress("jetsJPTAK5_n90", &jetsJPTAK5_n90, &b_jetsJPTAK5_n90);
     chain->SetBranchAddress("jetsJPTAK5_area", &jetsJPTAK5_area, &b_jetsJPTAK5_area);
     chain->SetBranchAddress("jetsJPTAK5_mass", &jetsJPTAK5_mass, &b_jetsJPTAK5_mass);
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
   
   if( m_metAlgo=="KT4" ) {
     chain->SetBranchAddress("NmetsKT4", &NmetsKT4, &b_NmetsKT4);
     chain->SetBranchAddress("metsKT4_et", &metsKT4_et, &b_metsKT4_et);
     chain->SetBranchAddress("metsKT4_phi", &metsKT4_phi, &b_metsKT4_phi);
     chain->SetBranchAddress("metsKT4_ex", &metsKT4_ex, &b_metsKT4_ex);
     chain->SetBranchAddress("metsKT4_ey", &metsKT4_ey, &b_metsKT4_ey);
     chain->SetBranchAddress("metsKT4_sumEt", &metsKT4_sumEt, &b_metsKT4_sumEt);
     chain->SetBranchAddress("metsKT4_et_JESCor", &metsKT4_et_JESCor, &b_metsKT4_et_JESCor);
     chain->SetBranchAddress("metsKT4_phi_JESCor", &metsKT4_phi_JESCor, &b_metsKT4_phi_JESCor);
   }
   else if( m_metAlgo=="KT6" ) {
     chain->SetBranchAddress("NmetsKT6", &NmetsKT6, &b_NmetsKT6);
     chain->SetBranchAddress("metsKT6_et", &metsKT6_et, &b_metsKT6_et);
     chain->SetBranchAddress("metsKT6_phi", &metsKT6_phi, &b_metsKT6_phi);
     chain->SetBranchAddress("metsKT6_ex", &metsKT6_ex, &b_metsKT6_ex);
     chain->SetBranchAddress("metsKT6_ey", &metsKT6_ey, &b_metsKT6_ey);
     chain->SetBranchAddress("metsKT6_sumEt", &metsKT6_sumEt, &b_metsKT6_sumEt);
     chain->SetBranchAddress("metsKT6_et_JESCor", &metsKT6_et_JESCor, &b_metsKT6_et_JESCor);
     chain->SetBranchAddress("metsKT6_phi_JESCor", &metsKT6_phi_JESCor, &b_metsKT6_phi_JESCor);
   }
   else if( m_metAlgo=="SC5" ) {
     chain->SetBranchAddress("NmetsSC5", &NmetsSC5, &b_NmetsSC5);
     chain->SetBranchAddress("metsSC5_et", &metsSC5_et, &b_metsSC5_et);
     chain->SetBranchAddress("metsSC5_phi", &metsSC5_phi, &b_metsSC5_phi);
     chain->SetBranchAddress("metsSC5_ex", &metsSC5_ex, &b_metsSC5_ex);
     chain->SetBranchAddress("metsSC5_ey", &metsSC5_ey, &b_metsSC5_ey);
     chain->SetBranchAddress("metsSC5_sumEt", &metsSC5_sumEt, &b_metsSC5_sumEt);
     chain->SetBranchAddress("metsSC5_et_JESCor", &metsSC5_et_JESCor, &b_metsSC5_et_JESCor);
     chain->SetBranchAddress("metsSC5_phi_JESCor", &metsSC5_phi_JESCor, &b_metsSC5_phi_JESCor); 
   } 
   else if( m_metAlgo=="SC7" ) {
     chain->SetBranchAddress("NmetsSC7", &NmetsSC7, &b_NmetsSC7);
     chain->SetBranchAddress("metsSC7_et", &metsSC7_et, &b_metsSC7_et);
     chain->SetBranchAddress("metsSC7_phi", &metsSC7_phi, &b_metsSC7_phi);
     chain->SetBranchAddress("metsSC7_ex", &metsSC7_ex, &b_metsSC7_ex);
     chain->SetBranchAddress("metsSC7_ey", &metsSC7_ey, &b_metsSC7_ey);
     chain->SetBranchAddress("metsSC7_sumEt", &metsSC7_sumEt, &b_metsSC7_sumEt);
     chain->SetBranchAddress("metsSC7_et_JESCor", &metsSC7_et_JESCor, &b_metsSC7_et_JESCor);
     chain->SetBranchAddress("metsSC7_phi_JESCor", &metsSC7_phi_JESCor, &b_metsSC7_phi_JESCor); 
   }
   if( m_metAlgo=="tcmet" ) {
     chain->SetBranchAddress("Ntcmets", &Ntcmets, &b_Ntcmets);
     chain->SetBranchAddress("tcmets_et", &tcmets_et, &b_tcmets_et);
     chain->SetBranchAddress("tcmets_phi", &tcmets_phi, &b_tcmets_phi);
     chain->SetBranchAddress("tcmets_ex", &tcmets_ex, &b_tcmets_ex);
     chain->SetBranchAddress("tcmets_ey", &tcmets_ey, &b_tcmets_ey);
     chain->SetBranchAddress("tcmets_sumEt", &tcmets_sumEt, &b_tcmets_sumEt);
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
   chain->SetBranchAddress("mus_id_TMLastStationOptimizedLowPtLoose", &mus_id_TMLastStationOptimizedLowPtLoose, &b_mus_id_TMLastStationOptimizedLowPtLoose);
   chain->SetBranchAddress("mus_id_TMLastStationOptimizedLowPtTight", &mus_id_TMLastStationOptimizedLowPtTight, &b_mus_id_TMLastStationOptimizedLowPtTight);
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
   if(m_useMisslayers) chain->SetBranchAddress("tracks_innerLayerMissingHits",&tracks_innerLayerMissingHits, &b_tracks_innerLayerMissingHits);
   chain->SetBranchAddress("run", &run, &b_run);
   chain->SetBranchAddress("event", &event, &b_event);
   chain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);

   // HLT tree (nonIso ele trigger)
   //chain2->SetBranchAddress("HLT_Ele10_SW_L1R", &HLT_Ele10_SW_L1R, &b_HLT_Ele10_SW_L1R); //(8E29, startup)
   //cout << "GetTrigger()" << GetTrigger()<< endl;
   if( GetTrigger() ) {
     if(HLTBit=="HLT_Ele15_LW_L1R") chain2->SetBranchAddress("HLT_Ele15_LW_L1R", &HLT_Ele15_LW_L1R, &b_HLT_Ele15_LW_L1R); //v2 (1E31, ideal)
     if(HLTBit=="HLT_Ele15_SW_L1R") chain2->SetBranchAddress("HLT_Ele15_SW_L1R", &HLT_Ele15_SW_L1R, &b_HLT_Ele15_SW_L1R);
   }
   
   ///------------------------  MC Truth info  ------------------------------------
   if( !IsData() ) {
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

     if( m_studyPDFunc ){
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
