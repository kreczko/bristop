## HLT paths (label=HLT)

import FWCore.ParameterSet.Config as cms

TriggerPSets = cms.PSet(


      HLT_Activity_CSC = cms.PSet(
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

      HLT_Activity_Ecal_SC7 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Activity_Ecal_SC17 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Jet6U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Jet10U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Jet15U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Jet15U_HcalNoiseFiltered = cms.PSet(
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

      HLT_Jet70U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Jet100U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DiJetAve15U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DiJetAve30U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DiJetAve50U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DiJetAve70U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoubleJet15U_ForwardBackward = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoubleJet25U_ForwardBackward = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_QuadJet15U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_QuadJet20U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_QuadJet25U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1ETT100 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_EcalOnly_SumEt160 = cms.PSet(
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

      HLT_MET65 = cms.PSet(
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

      HLT_HT120U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_HT140U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1MuOpen = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1MuOpen_DT = cms.PSet(
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

      HLT_L2Mu0_NoVertex = cms.PSet(
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

      HLT_L2Mu25 = cms.PSet(
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

      HLT_Mu7 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu9 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu11 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_IsoMu9 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu20_NoVertex = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1DoubleMuOpen = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L2DoubleMu0 = cms.PSet(
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

      HLT_Mu5_L2Mu0 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu0_Track0_Jpsi = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu0_TkMu0_OST_Jpsi = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu3_Track3_Jpsi = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu3_TkMu0_OST_Jpsi = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu5_Track0_Jpsi = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu5_TkMu0_OST_Jpsi = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1SingleEG2 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1SingleEG5 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1SingleEG8 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1DoubleEG5 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Ele10_SW_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Ele12_SW_TightEleId_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Ele12_SW_TightEleIdIsol_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Ele12_SW_TightEleIdIsol_NoDEtaInEE_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Ele17_SW_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Ele17_SW_CaloEleId_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),


      
      HLT_Ele17_SW_TightEleId_L1R = cms.PSet(
    src = cms.InputTag("TriggerResults","","HLT"),
    method = cms.string('HLTBitVariable')
    ),
      
      

      HLT_Ele17_SW_LooseEleId_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Ele17_SW_EleId_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Ele22_SW_CaloEleId_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Ele22_SW_TighterEleId_L1R_v2 = cms.PSet(
    src = cms.InputTag("TriggerResults","","HLT"),
    method = cms.string('HLTBitVariable')
    ),
      
      
      HLT_Ele40_SW_L1R = cms.PSet(
    src = cms.InputTag("TriggerResults","","HLT"),
    method = cms.string('HLTBitVariable')
    
    ),

      HLT_DoubleEle4_SW_eeRes_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoubleEle10_SW_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Photon10_Cleaned_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Photon15_Cleaned_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Photon20_NoHE_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Photon20_Cleaned_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Photon30_Cleaned_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Photon50_NoHE_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Photon50_NoHE_Cleaned_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoublePhoton5_CEP_L1R = cms.PSet(
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

      HLT_DoublePhoton15_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoublePhoton17_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_SingleIsoTau20_Trk5_MET20 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_SingleIsoTau20_Trk15_MET20 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_SingleIsoTau30_Trk5_MET20 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_SingleIsoTau30_Trk5_L120or30 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoubleIsoTau15_OneLeg_Trk5 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_DoubleIsoTau15_Trk5 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_BTagMu_Jet10U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_BTagMu_Jet20U = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_StoppedHSCP = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L2Mu5_Photon9_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_Mu5_Photon9_Cleaned_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_ZeroBias = cms.PSet(
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

      HLT_MultiVertex6 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_MultiVertex8_L1ETT60 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1_BptxXOR_BscMinBiasOR = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Tech_BSC_minBias_OR = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Tech_BSC_minBias = cms.PSet(
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

      HLT_L1Tech_BSC_HighMultiplicity = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Tech_RPC_TTU_RBst1_collisions = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_L1Tech_HCAL_HF = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_TrackerCosmics = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_RPCBarrelCosmics = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_IsoTrackHE = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_IsoTrackHB = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_HcalPhiSym = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_HcalNZS = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_PixelTracks_Multiplicity70 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_PixelTracks_Multiplicity85 = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_PixelTracks_Multiplicity100 = cms.PSet(
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

      HLT_DTErrors = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      ),

      HLT_LogMonitor = cms.PSet(
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

      ##extra triggers for earlier runs
      HLT_Ele10_LW_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')
      ),


      HLT_Ele15_LW_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')
      ),

      HLT_Ele15_SW_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')
      ),

      HLT_Ele15_SW_CaloEleId_L1R = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')
      ),


      HLT_HLTriggerFinalPath = cms.PSet(
         src = cms.InputTag("TriggerResults","","HLT"),
         method = cms.string('HLTBitVariable')

      )
      
)
