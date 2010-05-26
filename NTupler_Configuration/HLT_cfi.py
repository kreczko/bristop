## HLT paths (label=HLT)

import FWCore.ParameterSet.Config as cms

TriggerPSets = cms.PSet(


      HLT_Activity_L1A = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Activity_PixelClusters = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Activity_DT = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Activity_DT_Tuned = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Activity_Ecal = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Activity_EcalREM = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_SelectEcalSpikes_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_SelectEcalSpikesHighEt_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Jet6U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Jet6U_NoBPTX = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Jet10U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Jet10U_NoBPTX = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Jet15U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Jet30U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Jet50U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1SingleForJet = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1SingleForJet_NoBPTX = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1SingleCenJet = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1SingleCenJet_NoBPTX = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1SingleTauJet = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1SingleTauJet_NoBPTX = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_FwdJet20U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DiJetAve15U_8E29 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DiJetAve30U_8E29 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoubleJet15U_ForwardBackward = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_QuadJet15U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1MET20 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_MET45 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_MET100 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_HT100U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1MuOpen = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1MuOpen_NoBPTX = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1MuOpen_AntiBPTX = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Mu = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Mu20 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L2Mu0 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L2Mu3 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L2Mu9 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L2Mu11 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L2DoubleMu0 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_IsoMu3 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu3 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu5 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu9 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1DoubleMuOpen = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoubleMu0 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoubleMu3 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu0_L1MuOpen = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu3_L1MuOpen = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu5_L1MuOpen = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu0_L2Mu0 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu3_L2Mu0 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu5_L2Mu0 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu0_Track0_Jpsi = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu3_Track0_Jpsi = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu5_Track0_Jpsi = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1SingleEG2 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1SingleEG2_NoBPTX = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1SingleEG5 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1SingleEG5_NoBPTX = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1SingleEG8 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1SingleEG20_NoBPTX = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1DoubleEG5 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_EgammaSuperClusterOnly_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Ele10_LW_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Ele10_LW_EleId_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Ele15_LW_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Ele15_SC10_LW_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Ele15_SiStrip_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Ele20_LW_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoubleEle5_SW_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Photon10_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Photon15_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Photon15_TrackIso_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Photon15_LooseEcalIso_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Photon20_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Photon30_L1R_8E29 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoublePhoton4_eeRes_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoublePhoton4_Jpsi_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoublePhoton4_Upsilon_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoublePhoton5_Jpsi_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoublePhoton5_Upsilon_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoublePhoton5_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoublePhoton10_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_SingleLooseIsoTau20 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoubleLooseIsoTau15 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_BTagIP_Jet50U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_BTagMu_Jet10U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_StoppedHSCP_8E29 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Mu14_L1SingleEG10 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Mu14_L1SingleJet6U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Mu14_L1ETM30 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_ZeroBias = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_MinBiasBSC = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_MinBiasBSC_NoBPTX = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_MinBiasBSC_OR = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_MinBiasHcal = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_MinBiasEcal = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_ZeroBiasPixel_SingleTrack = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_MinBiasPixel_SingleTrack = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_MinBiasPixel_DoubleTrack = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_MinBiasPixel_DoubleIsoTrack5 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_CSCBeamHalo = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_CSCBeamHaloOverlapRing1 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_CSCBeamHaloOverlapRing2 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_CSCBeamHaloRing2or3 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_BackwardBSC = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_ForwardBSC = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_HighMultiplicityBSC = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_SplashBSC = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1_BscMinBiasOR_BptxPlusORMinus = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1_BscMinBiasOR_BptxPlusORMinus_NoBPTX = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_RPCBarrelCosmics = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_TrackerCosmics = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Tech_BSC_halo = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Tech_BSC_halo_forPhysicsBackground = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Tech_RPC_TTU_RBst1_collisions = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_IsoTrackHE_8E29 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_IsoTrackHB_8E29 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_HcalPhiSym = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_HcalNZS_8E29 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DTErrors = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1_BscMinBiasOR_BeamGas = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_HighMult40 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Calibration = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_EcalCalibration = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_HcalCalibration = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Random = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1_HFtech = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Tech_HCAL_HF_coincidence_PM = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_GlobalRunHPDNoise = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_TechTrigHCALNoise = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1_BPTX = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1_BPTX_MinusOnly = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1_BPTX_PlusOnly = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L2Mu0_NoVertex = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_TkMu3_NoVertex = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_PhysicsDeclared = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_LogMonitor = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLTriggerFinalPath = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')
      )

)
