import FWCore.ParameterSet.Config as cms

##
## TL added:
##    - pv_rho; jets: CHF,NHF, CEF, NEF
##    - els: scEta, scPhi
##
## Changes for 3_6_0 (to run on 35X-reco input):
##    - keep 'btag_secVertex:bDiscriminator("simpleSecondaryVertexBJetTags")'
##    - removed SC5, SC7 jetMET
##    - review btag var
##


myHLT = "REDIGI"

### put the usual cff fragment here
InputTagDistributorService = cms.Service("InputTagDistributorService")

VariableHelperService = cms.Service("VariableHelperService")

UpdaterService = cms.Service("UpdaterService")

TFileService = cms.Service("TFileService",
                           fileName = cms.string('ntuple.root')
                           )

#########################################################################
### below are the usual ntupler cfi
#########################################################################
basicKinematicLeaves = cms.PSet(
    status = cms.string('status'),
    phi = cms.string('phi'),
    pt = cms.string('pt'),
    pz = cms.string('pz'),
    px = cms.string('px'),
    py = cms.string('py'),
    eta = cms.string('eta'),
    theta = cms.string('theta'),
    et = cms.string('et'),
    energy = cms.string('energy')
)



if myHLT=="HLT":
    from HLT_cfi import TriggerPSets
else:
    from HLT_redigi_cfi import TriggerPSets



## TL added 7-1-10 to simplify code
met_commonVarList = cms.vstring(
    'et:et',
    'phi:phi',
    'ex:px',
    'ey:py',
    'sign:metSignificance',
    'sumEt:sumEt',
    'unCPt:uncorrectedPt',
    'unCPhi:uncorrectedPhi',
    'et_muonCor:uncorrectedPt("uncorrJES")',
    'phi_muonCor:uncorrectedPhi("uncorrJES")',
    'et_JESCor:uncorrectedPt("uncorrMUON")',
    'phi_JESCor:uncorrectedPhi("uncorrMUON")'
    )

## Ref: http://cmslxr.fnal.gov/lxr/source/PhysicsTools/PatAlgos/python/tools/jetTools.py?v=CMSSW_3_5_7
jet_btagInfo = cms.vstring(
    'btag_TC_highPur:bDiscriminator("trackCountingHighPurBJetTags")',
    'btag_TC_highEff:bDiscriminator("trackCountingHighEffBJetTags")',
    'btag_jetProb:bDiscriminator("jetProbabilityBJetTags")',
    'btag_jetBProb:bDiscriminator("jetBProbabilityBJetTags")',
     ## soft ele
    #'btag_softEleByPt:bDiscriminator("softElectronByPtBJetTags")',
    #'btag_softEleByIP3d:bDiscriminator("softElectronByIP3dBJetTags")',
     ## soft muon
    #'btag_softMuon:bDiscriminator("softMuonBJetTags")',
    #'btag_softMuonByPt:bDiscriminator("softMuonByPtBJetTags")',
    #'btag_softMuonByIP3d:bDiscriminator("softMuonByIP3dBJetTags")',
    ## SSV
    'btag_secVertex:bDiscriminator("simpleSecondaryVertexBJetTags")', #RECO 3_5_X = SSVHE
    'btag_ssvHE:bDiscriminator("simpleSecondaryVertexHighEffBJetTags")', #RECO 3_6_X
    'btag_ssvHP:bDiscriminator("simpleSecondaryVertexHighPurBJetTags")', #RECO 3_6_X
    ## CSV
    #'btag_csv:bDiscriminator("combinedSecondaryVertexBJetTags")',
    #'btag_csvMVA:bDiscriminator("combinedSecondaryVertexMVABJetTags")'
)

jet_commonVarList = cms.vstring(    
    'chgEmE:chargedEmEnergy',
    'chgHadE:chargedHadronEnergy',
    'chgMuE:chargedMuEnergy',
    'chg_Mult:chargedMultiplicity',#NCH
    'neutralEmE:neutralEmEnergy',
    'neutralHadE:neutralHadronEnergy',
    'neutral_Mult:neutralMultiplicity',
    'mu_Mult:muonMultiplicity',
    'emf:emEnergyFraction',
    'ehf:energyFractionHadronic',
    'n60:n60',
    'n90:n90',
    'nConstituents:nConstituents', ##TL
    'area:towersArea',
    'mass:mass',
    'id_fHPD:jetID.fHPD',  ## jetID
    'id_fRBX:jetID.fRBX',
    'id_n90hits:jetID.n90Hits',
    'id_fSubDetector1:jetID.fSubDetector1',
    'id_fSubDetector2:jetID.fSubDetector2',
    'id_fSubDetector3:jetID.fSubDetector3',
    'id_fSubDetector4:jetID.fSubDetector4',
    'id_restrictedEMF:jetID.restrictedEMF',
    'id_nHCALTowers:jetID.nHCALTowers',
    'id_nECALTowers:jetID.nECALTowers',
    'id_approximatefHPD:jetID.approximatefHPD',
    'id_approximatefRBX:jetID.approximatefRBX',
    'id_hitsInN90:jetID.hitsInN90',
    'charge:charge',##
    'etaetaMoment:etaetaMoment',
    'etaphiMoment:etaphiMoment',
    'phiphiMoment:phiphiMoment',
    'pileup:pileup',
    'towersArea:towersArea',
    'sv_nVertices:tagInfoSecondaryVertex("secondaryVertex").nVertices', ## SV num
    'sv_nSelectedTracks:tagInfoSecondaryVertex("secondaryVertex").nSelectedTracks',
    'sv_nVertexTracks:tagInfoSecondaryVertex("secondaryVertex").nVertexTracks',
    'sv_hasTracks:tagInfoSecondaryVertex("secondaryVertex").hasTracks',
    'sv0_ntk:tagInfoSecondaryVertex("secondaryVertex").secondaryVertex(0).tracksSize()', ##dist look ok when nVertices>0.
    'sv0_chi2:tagInfoSecondaryVertex("secondaryVertex").secondaryVertex(0).chi2()',
    'sv0_ndof:tagInfoSecondaryVertex("secondaryVertex").secondaryVertex(0).ndof()',
    'numberOfDaughters:numberOfDaughters',
    'RawJetEnergy:correctedJet("pat::JetCorrFactors::Raw").energy',
    'CHF:correctedJet("raw","").chargedHadronEnergyFraction',
    'NHF:correctedJet("raw","").neutralHadronEnergyFraction',
    'NEF:correctedJet("raw","").neutralEmEnergyFraction',
    'CEF:correctedJet("raw","").chargedEmEnergyFraction'
    )

jet_mcInfo = cms.vstring(
    'parton_Id:genParton.pdgId',
    'parton_motherId:genParton.mother.pdgId', 
    'parton_pt:genParton.pt', 
    'parton_phi:genParton.phi', 
    'parton_eta:genParton.eta', 
    'parton_Energy:genParton.energy', 
    'parton_mass:genParton.mass',  
    'parton_motherID:genParton.mother.pdgId',
    #'parton_grandmotherID:genParton.mother.mother.pdgId',
    'gen_et:genJet.et', 
    'gen_pt:genJet.pt', 
    'gen_eta:genJet.eta', 
    'gen_phi:genJet.phi', 
    'gen_mass:genJet.mass', 
    'gen_Energy:genJet.energy', 
    'gen_Id:genJet.pdgId', 
    'gen_motherID:genJet.mother.pdgId', 
    'gen_threeCharge:genJet.threeCharge',
    'partonFlavour:partonFlavour'
    )
default_jet_var = cms.vstring()
default_jet_var.extend( jet_commonVarList )
default_jet_var.extend( jet_btagInfo )
default_jet_var.extend( jet_mcInfo )



# note: no btag or mc info for extra jet collections (exception: PF)
extra_jet_var = jet_commonVarList


ele_vars = cms.vstring(
    ## generator
    'gen_id:genLepton.pdgId',
    'gen_phi:genLepton.phi',
    'gen_pt:genLepton.pt',
    'gen_pz:genLepton.pz',
    'gen_px:genLepton.px',
    'gen_py:genLepton.py',
    'gen_eta:genLepton.eta',
    'gen_theta:genLepton.theta',
    'gen_et:genLepton.et',                       
    'gen_mother_id:genLepton.mother.pdgId',
    'gen_mother_phi:genLepton.mother.phi',
    'gen_mother_pt:genLepton.mother.pt',
    'gen_mother_pz:genLepton.mother.pz',
    'gen_mother_px:genLepton.mother.px',
    'gen_mother_py:genLepton.mother.py',
    'gen_mother_eta:genLepton.mother.eta',
    'gen_mother_theta:genLepton.mother.theta',
    'gen_mother_et:genLepton.mother.et',                       
    ## isolation
    'cIso:caloIso', 
    'tIso:trackIso', 
    'ecalIso:ecalIso',
    'hcalIso:hcalIso',
    'dr03EcalRecHitSumEt:dr03EcalRecHitSumEt',
    'dr03HcalTowerSumEt:dr03HcalTowerSumEt',
    'dr03TkSumPt:dr03TkSumPt',
    'dr03HcalDepth1TowerSumEt:dr03HcalDepth1TowerSumEt',
    'dr03HcalDepth2TowerSumEt:dr03HcalDepth2TowerSumEt',                        
    'dr04EcalRecHitSumEt:dr04EcalRecHitSumEt',                        
    'dr04HcalTowerSumEt:dr04HcalTowerSumEt',
    'dr04HcalDepth1TowerSumEt:dr04HcalDepth1TowerSumEt',
    'dr04HcalDepth2TowerSumEt:dr04HcalDepth2TowerSumEt',                        
    'dr04TkSumPt:dr04TkSumPt',
    ## ID
    'tightId:electronID("eidTight")',
    'looseId:electronID("eidLoose")',
    'robustTightId:electronID("eidRobustTight")',
    'robustLooseId:electronID("eidRobustLoose")',
    'robustHighEnergyId:electronID("eidRobustHighEnergy")',
    'closestCtfTrackRef:closestCtfTrackRef.key',
    'shFracInnerHits:shFracInnerHits',   
    ## flag
    'isEcalDriven:ecalDrivenSeed',#isEcalDriven',
    'isTrackerDriven:trackerDrivenSeed',#isTrackerDriven',
    'isEE:isEE',
    'isEEGap:isEEGap',
    'isEB:isEB',
    'isEBGap:isEEGap',
    'isConvertedPhoton:isConvertedPhoton',   
    'chi2:gsfTrack.chi2', 
    #'class:classification', 
    'charge:charge', 
    'caloEnergy:caloEnergy', 
    'hadOverEm:hadronicOverEm', 
    'eOverPIn:eSuperClusterOverP', 
    'eSeedOverPOut:eSeedClusterOverPout', 
    'eSCraw:superCluster.rawEnergy', 
    'eSeed:superCluster.seed.energy', 
    'sigmaEtaEta:scSigmaEtaEta',
    'sigmaIEtaIEta:scSigmaIEtaIEta',
    'scE1x5:scE1x5',
    'scE2x5Max:scE2x5Max',
    'scE5x5:scE5x5',
    'dEtaIn:deltaEtaSuperClusterTrackAtVtx', 
    'dPhiIn:deltaPhiSuperClusterTrackAtVtx', 
    'dEtaOut:deltaEtaSeedClusterTrackAtCalo', 
    'dPhiOut:deltaPhiSeedClusterTrackAtCalo', 
    'numvalhits:gsfTrack.numberOfValidHits', 
    'numlosthits:gsfTrack.numberOfLostHits', 
    #'numCluster:numberOfClusters', 
    'basicClustersSize:basicClustersSize',    
    'tk_theta:gsfTrack.theta',
    'tk_charge:gsfTrack.charge',
    'tk_pt:gsfTrack.pt', 
    'tk_phi:gsfTrack.phi', 
    'tk_eta:gsfTrack.eta', 
    'd0dum:gsfTrack.d0',
    'dz:gsfTrack.dz', 
    'vx:gsfTrack.vx', 
    'vy:gsfTrack.vy', 
    'vz:gsfTrack.vz', 
    'ndof:gsfTrack.ndof', 
    'ptError:gsfTrack.ptError', 
    'd0dumError:gsfTrack.d0Error', 
    'dzError:gsfTrack.dzError', 
    'etaError:gsfTrack.etaError', 
    'phiError:gsfTrack.phiError', 
    'cpx:trackMomentumAtCalo.x', 
    'cpy:trackMomentumAtCalo.y', 
    'cpz:trackMomentumAtCalo.z', 
    'vpx:trackMomentumAtVtx.x', 
    'vpy:trackMomentumAtVtx.y', 
    'vpz:trackMomentumAtVtx.z', 
    'cx:TrackPositionAtCalo.x', 
    'cy:TrackPositionAtCalo.y', 
    'cz:TrackPositionAtCalo.z', 
    # new addition (18-12-09): hit pattern 
    'innerLayerMissingHits:gsfTrack.trackerExpectedHitsInner().numberOfHits()',
    'outerLayerMissingHits:gsfTrack.trackerExpectedHitsOuter().numberOfHits()',
    # added 01-04-10, further innerHit pattern info
    'innerLayerMissingHitsPixels:gsfTrack.trackerExpectedHitsInner().numberOfValidPixelHits()',
    'innerLayerMissingHitsStrips:gsfTrack.trackerExpectedHitsInner().numberOfValidStripHits()',
    'numPixelHits:gsfTrack.hitPattern.numberOfValidPixelHits()',
    'numStripHits:gsfTrack.hitPattern.numberOfValidStripHits()',
    'fbrem:fbrem',
    'dB:dB',    #d0 w.r.t beam spot/PV (PAT)
    'edB:edB',  #d0 error w.r.t beam spot/PV (PAT)
    'scEta:superCluster.eta', #supercluster eta
    'scPhi:superCluster.phi'  #supercluster phi
    )





muon_vars = cms.vstring(
    'gen_id:genLepton.pdgId',
    'gen_phi:genLepton.phi',
    'gen_pt:genLepton.pt',
    'gen_pz:genLepton.pz',
    'gen_px:genLepton.px',
    'gen_py:genLepton.py',
    'gen_eta:genLepton.eta',
    'gen_theta:genLepton.theta',
    'gen_et:genLepton.et',                       
    'gen_mother_id:genLepton.mother.pdgId',
    'gen_mother_phi:genLepton.mother.phi',
    'gen_mother_pt:genLepton.mother.pt',
    'gen_mother_pz:genLepton.mother.pz',
    'gen_mother_px:genLepton.mother.px',
    'gen_mother_py:genLepton.mother.py',
    'gen_mother_eta:genLepton.mother.eta',
    'gen_mother_theta:genLepton.mother.theta',
    'gen_mother_et:genLepton.mother.et',                       
    'tkHits:track.hitPattern.numberOfValidHits', 
    'cIso:caloIso', 
    'tIso:trackIso',
    'ecalIso:ecalIso',
    'hcalIso:hcalIso',
    'ecalvetoDep:ecalIsoDeposit.candEnergy',
    'hcalvetoDep:hcalIsoDeposit.candEnergy',    
    'calEnergyEm:calEnergy.em',
    'calEnergyHad:calEnergy.had',
    'calEnergyHo:calEnergy.ho',
    'calEnergyEmS9:calEnergy.emS9',
    'calEnergyHadS9:calEnergy.hadS9',
    'calEnergyHoS9:calEnergy.hoS9',
    'iso03_sumPt:isolationR03.sumPt',
    'iso03_emEt:isolationR03.emEt',
    'iso03_hadEt:isolationR03.hadEt',
    'iso03_hoEt:isolationR03.hoEt',
    'iso03_nTracks:isolationR03.nTracks',
    'iso05_sumPt:isolationR05.sumPt',
    'iso05_emEt:isolationR05.emEt',
    'iso05_hadEt:isolationR05.hadEt',
    'iso05_hoEt:isolationR05.hoEt',
    'iso05_nTracks:isolationR05.nTracks',
    'charge:charge', 
    'cm_chi2:combinedMuon.chi2', 
    'cm_ndof:combinedMuon.ndof', 
    'cm_chg:combinedMuon.charge', 
    'cm_pt:combinedMuon.pt', 
    'cm_px:combinedMuon.px', 
    'cm_py:combinedMuon.py', 
    'cm_pz:combinedMuon.pz', 
    'cm_eta:combinedMuon.eta', 
    'cm_phi:combinedMuon.phi', 
    'cm_theta:combinedMuon.theta', 
    'cm_d0dum:combinedMuon.d0', 
    'cm_dz:combinedMuon.dz', 
    'cm_vx:combinedMuon.vx', 
    'cm_vy:combinedMuon.vy', 
    'cm_vz:combinedMuon.vz', 
    'cm_numvalhits:combinedMuon.numberOfValidHits', 
    'cm_numlosthits:combinedMuon.numberOfLostHits', 
    'cm_d0dumErr:combinedMuon.d0Error', 
    'cm_dzErr:combinedMuon.dzError', 
    'cm_ptErr:combinedMuon.ptError', 
    'cm_etaErr:combinedMuon.etaError', 
    'cm_phiErr:combinedMuon.phiError', 
    'tk_chi2:track.chi2',
    'tk_ndof:track.ndof', 
    'tk_chg:track.charge', 
    'tk_pt:track.pt', 
    'tk_px:track.px', 
    'tk_py:track.py', 
    'tk_pz:track.pz', 
    'tk_eta:track.eta', 
    'tk_phi:track.phi', 
    'tk_theta:track.theta', 
    'tk_d0dum:track.d0', 
    'tk_dz:track.dz', 
    'tk_vx:track.vx', 
    'tk_vy:track.vy', 
    'tk_vz:track.vz', 
    'tk_numvalhits:track.numberOfValidHits', 
    'tk_numlosthits:track.numberOfLostHits', 
    'tk_d0dumErr:track.d0Error', 
    'tk_dzErr:track.dzError', 
    'tk_ptErr:track.ptError', 
    'tk_etaErr:track.etaError', 
    'tk_phiErr:track.phiError', 
    'stamu_chi2:standAloneMuon.chi2', 
    'stamu_ndof:standAloneMuon.ndof', 
    'stamu_chg:standAloneMuon.charge', 
    'stamu_pt:standAloneMuon.pt', 
    'stamu_px:standAloneMuon.px', 
    'stamu_py:standAloneMuon.py', 
    'stamu_pz:standAloneMuon.pz', 
    'stamu_eta:standAloneMuon.eta', 
    'stamu_phi:standAloneMuon.phi', 
    'stamu_theta:standAloneMuon.theta', 
    'stamu_d0dum:standAloneMuon.d0',
    'stamu_dz:standAloneMuon.dz', 
    'stamu_vx:standAloneMuon.vx', 
    'stamu_vy:standAloneMuon.vy', 
    'stamu_vz:standAloneMuon.vz', 
    'stamu_numvalhits:standAloneMuon.numberOfValidHits', 
    'stamu_numlosthits:standAloneMuon.numberOfLostHits', 
    'stamu_d0dumErr:standAloneMuon.d0Error', 
    'stamu_dzErr:standAloneMuon.dzError', 
    'stamu_ptErr:standAloneMuon.ptError', 
    'stamu_etaErr:standAloneMuon.etaError', 
    'stamu_phiErr:standAloneMuon.phiError', 
    'num_matches:numberOfMatches',
    'id_All:isGood("All")',
    'id_AllGlobalMuons:isGood("AllGlobalMuons")',
    'id_AllStandAloneMuons:isGood("AllStandAloneMuons")',
    'id_AllTrackerMuons:isGood("AllTrackerMuons")',
    'id_TrackerMuonArbitrated:isGood("TrackerMuonArbitrated")',
    'id_AllArbitrated:isGood("AllArbitrated")',
    'id_GlobalMuonPromptTight:isGood("GlobalMuonPromptTight")',
    'id_TMLastStationLoose:isGood("TMLastStationLoose")',                        
    'id_TMLastStationTight:isGood("TMLastStationTight")',
    'id_TM2DCompatibilityLoose:isGood("TM2DCompatibilityLoose")',
    'id_TM2DCompatibilityTight:isGood("TM2DCompatibilityTight")',                        
    'id_TMOneStationLoose:isGood("TMOneStationLoose")',
    'id_TMOneStationTight:isGood("TMOneStationTight")',
    'id_TMLastStationOptimizedLowPtLoose:isGood("TMLastStationOptimizedLowPtLoose")',
    'id_TMLastStationOptimizedLowPtTight:isGood("TMLastStationOptimizedLowPtTight")',
    'globalTk_normChi2:globalTrack.normalizedChi2()',
    'globalTk_Chi2:globalTrack.chi2()',
    'dB:dB',
    'edB:edB'
    )


# MC only
from PDFw_cfi import PDFPSets

# MC genParticles
mc_doc_vars = cms.vstring(
    'id:pdgId', 
    'pt:pt', 
    'px:px', 
    'py:py', 
    'pz:pz',
    'eta:eta',
    'phi:phi',
    'theta:theta', 
    'energy:energy',
    'status:status', 
    'charge:charge', 
    'mother_id:mother.pdgId',
    'grandmother_id:mother.mother.pdgId', 
    'ggrandmother_id:mother.mother.mother.pdgId', 
    'mother_pt:mother.pt',           
    'vertex_x:vertex.x', 
    'vertex_y:vertex.y', 
    'vertex_z:vertex.z', 
    'mass:mass', 
    'numOfDaughters:numberOfDaughters',
    'numOfMothers:numberOfMothers'
    )



pv_vars = cms.vstring(
    'x:x',
    'y:y',
    'z:z',
    'xErr:xError',
    'yErr:yError',
    'zErr:zError',
    'chi2:chi2',
    'ndof:ndof',
    'normalizedChi2:normalizedChi2',#
    'tracksSize:tracksSize',
    'isValid:isValid',#
    'isFake:isFake', #
    'rho:position.rho' #TL 16 Apr    
)



#Main module to create the two trees and histograms (these histograms are not from the trees)
configurableAnalysis = cms.EDFilter("ConfigurableAnalysis",
    Selections = cms.PSet(
        filters = cms.PSet(
              none = cms.PSet(),
#             leadingMuon = cms.PSet(
#                 src = cms.string('muons'),
# #                cut = cms.vstring('pt > 20.0 & abs(eta) < 3.0'),
#                 cut = cms.vstring(''),
#                 selector = cms.string('patMuonSEventSelector')
# #                selector = cms.string('')
#             )#,
# #            leadingElectron = cms.PSet(
# #                src = cms.string('electrons'),
# #                cut = cms.vstring('pt > 20.0 & abs(eta) < 3.0'),
# #                selector = cms.string('patElectronSEventSelector')
# #            )
        ),
        selections = cms.PSet(
            minSelection = cms.PSet(
#                 #vstring filterOrder = { "leadingElectron" }
                  filterOrder = cms.vstring('none'),
                  makeFinalPlots = cms.bool(False),
                  makeSummaryTable = cms.bool(True),
                  makeContentPlots = cms.bool(False),
                  makeAllButOnePlots = cms.bool(False),
                  ntuplize = cms.bool(True),
                  nMonitor = cms.uint32(1000),
                  makeCumulativePlots = cms.bool(False)
             ),
        ),
     ),
      
     Plotter = cms.PSet(),
#         TH1s = cms.PSet(

#         ),
#         ComponentName = cms.string('VariablePlotter'),
#         TProfiles = cms.PSet(

#         ),
#         TH2s = cms.PSet(

#         )
#     ),
                                    
    Variables = cms.PSet(
        TriggerPSets,

#         L1Bit = cms.PSet(
#             src = cms.InputTag("gtDigis"),
#             method = cms.string('ComputedVariable'),
#             computer = cms.string('L1BitComputer')
#             ),

        ##############################
        #   PDF Weights
        ##############################
        # TL 8 Feb 2010: PDF weights
        PDFPSets
                
            

     #   csaWeight = cms.PSet(
     #       src = cms.InputTag("csaweightproducer","weight"),
     #       method = cms.string('DoubleVar')
     #   ),
    #    procIDSplit = cms.PSet(
    #        weightLabel = cms.string('csaweightproducer'),
    #        maxID = cms.uint32(70),
    #        method = cms.string('ProcessIdSplitter'),
    #        lumi = cms.double(1000.0)
    #    )
    ),
    workAsASelector = cms.bool(False), #True
    flows = cms.vstring('minSelection'),
    InputTags = cms.PSet(

        # Uncomment the following 3 lines if running on MC
        genParticles = cms.InputTag("genParticles"),
        mets      = cms.InputTag("patMETs"),
        electrons = cms.InputTag("cleanPatElectrons"),
        ele2      = cms.InputTag("cleanPatEle2"), # 2nd ele for d0(offlinePrimaryVerticesWithBS)
        ele3      = cms.InputTag("cleanPatEle3"), # 3rd ele for d0(offlineBeamSpot)
        muons     = cms.InputTag("cleanPatMuons"),
        mu2       = cms.InputTag("cleanPatMu2"),
        mu3       = cms.InputTag("cleanPatMu3"),
        jets      = cms.InputTag("cleanPatJets"),

        ## Extra jet/met collections
        #metsKT4 = cms.InputTag("patMETsKT4Calo"),
        #metsKT6 = cms.InputTag("patMETsKT6Calo"),
        tcmets  = cms.InputTag("patMETsTC"),        
        #pfMets = cms.InputTag("patMETsPF"),       

        #jetsKT4    = cms.InputTag("cleanPatJetsKT4Calo"),
        #jetsKT6    = cms.InputTag("cleanPatJetsKT6Calo"),
        #jetsAK5JPT = cms.InputTag("cleanPatJetsAK5JPT"), #jet-plus-tracks
        jetsAK5JPT = cms.InputTag("patJetsAK5JPT"), #jet-plus-tracks
        #jetsPFAK5 = cms.InputTag("cleanPatJetsAK5PF"), ##note addJetCollection does not make MET for PF

        ## Particle Flow (PF2PAT) Kachanon
        pfMets         = cms.InputTag("patMETsPF"),
        pfJets         = cms.InputTag("selectedPatJetsPF"),
        pfElectronsIso = cms.InputTag("selectedPatElectronsPF"),
		pfElectronsAll = cms.InputTag("selectedPatElectronsToto"),
        pfMuons        = cms.InputTag("selectedPatMuonsPF"),
        mypfJets       = cms.InputTag("cleanPatJetsAK5PF"),
        
        # Others
        #conversions = cms.InputTag("conversions"),
        # syntax: cms.InputTag(label,instance,process)
        #Kshort = cms.InputTag("generalV0Candidates","Kshort"),
        #Lambda = cms.InputTag("generalV0Candidates","Lambda"),

        primaryVertex   = cms.InputTag("offlinePrimaryVertices"),
        primaryVertexBS = cms.InputTag("offlinePrimaryVerticesWithBS"),
        beamspot        = cms.InputTag("offlineBeamSpot")
        
    ),
                                    
    Ntupler = cms.PSet(
        branchesPSet = cms.PSet(
           treeName = cms.string('eventB'),

           pv = cms.PSet(
               src = cms.string('primaryVertex'),
                  leaves = cms.PSet(
                     vars = pv_vars                      
                  ),
               Class = cms.string('reco::Vertex')
           ),

           pvB = cms.PSet(
               src = cms.string('primaryVertexBS'),
                  leaves = cms.PSet(
                     vars = pv_vars                      
                  ),
               Class = cms.string('reco::Vertex')
           ),


           beamSpot = cms.PSet(  
           src = cms.string('beamspot'),
            leaves = cms.PSet(
                vars = cms.vstring(
                     'x:position.x',
                     'y:position.y',
                     'z:position.z',
                     'x0Error:x0Error',
                     'y0Error:y0Error',
                     'z0Error:z0Error',
                     'sigmaZ:sigmaZ',
                     'sigmaZ0Error:sigmaZ0Error',
                     'dxdz:dxdz',
                     'dxdzError:dxdzError',
                     'dydz:dydz',
                     'dydzError:dydzError',
                     'beamWidthX:BeamWidthX',#addedFB
                     'beamWidthY:BeamWidthY',#addedFB
                     'beamWidthXError:BeamWidthXError',#addedFB
                     'beamWidthYError:BeamWidthYError'#addedFB
                       )
                ),
            Class = cms.string('reco::BeamSpot')
            ),
                                                                              
            
            ## TL added 18 May
            hcalNoiseSummary = cms.PSet(
            src = cms.InputTag("hcalnoise"),
            leaves = cms.PSet(
                vars = cms.vstring(
                'passLooseNoiseFilter:passLooseNoiseFilter',
                'passTightNoiseFilter:passTightNoiseFilter',
                'passHighLevelNoiseFilter:passHighLevelNoiseFilter',
                'noiseFilterStatus:noiseFilterStatus',
                'noiseType:noiseType',
                'eventEMEnergy:eventEMEnergy',
                'eventHadEnergy:eventHadEnergy',
                'eventTrackEnergy:eventTrackEnergy',
                'eventEMFraction:eventEMFraction',
                'eventChargeFraction:eventChargeFraction',
                'min10GeVHitTime:min10GeVHitTime',
                'max10GeVHitTime:max10GeVHitTime',
                'rms10GeVHitTime:rms10GeVHitTime',
                'min25GeVHitTime:min25GeVHitTime',
                'max25GeVHitTime:max25GeVHitTime',
                'rms25GeVHitTime:rms25GeVHitTime',
                'num10GeVHits:num10GeVHits',
                'num25GeVHits:num25GeVHits',
                'minE2TS:minE2TS',
                'minE10TS:minE10TS',
                'minE2Over10TS:minE2Over10TS',
                'maxE2Over10TS:maxE2Over10TS',
                'maxZeros:maxZeros',
                'maxHPDHits:maxHPDHits',
                'maxRBXHits:maxRBXHits',
                'minHPDEMF:minHPDEMF',
                'minRBXEMF:minRBXEMF',
                'numProblematicRBXs:numProblematicRBXs',
                'maxHPDNoOtherHits:maxHPDNoOtherHits'
                )
                ),
            Class = cms.string('HcalNoiseSummary')
            ),
            
            hcalNoiseRBX = cms.PSet(
            src = cms.InputTag("hcalnoise"),
            leaves = cms.PSet(
                vars = cms.vstring(
                'idnumber:idnumber',
                'allChargeTotal:allChargeTotal',
                'allChargeHighest2TS:allChargeHighest2TS',
                'allChargeHighest3TS:allChargeHighest3TS',
                'totalZeros:totalZeros',
                'maxZeros:maxZeros',
                'recHitEnergy:recHitEnergy(1.5)',
                'minRecHitTime:minRecHitTime(20.0)',
                'maxRecHitTime:maxRecHitTime(20.0)',
                'numRecHits:numRecHits(1.5)',
                'caloTowerHadE:caloTowerHadE',
                'caloTowerEmE:caloTowerEmE',
                'caloTowerTotalE:caloTowerTotalE',
                'caloTowerEmFraction:caloTowerEmFraction'
                )
                ),
            Class = cms.string('reco::HcalNoiseRBX')
            ),
            


#            hemi = cms.PSet(
#                src = cms.InputTag("allLayer1Hemispheres"),
#                leaves = cms.PSet(
#                    basicKinematicLeaves
#                ),
#                Class = cms.string('pat::Hemisphere')
#            ),


             mus = cms.PSet(
             src = cms.string('muons'),
             leaves = cms.PSet(
                 basicKinematicLeaves,
                 vars = muon_vars
                 ),             
             Class = cms.string('pat::Muon')
             ),
             

            mus2 = cms.PSet(
            src = cms.string('mu2'),
            leaves = cms.PSet(
                vars = cms.vstring('px', 'py',
                                   'd0_pvbs:dB',
                                   'ed0_pvbs:edB')
                ),
            Class = cms.string('pat::Muon')
            ),

            mus3 = cms.PSet(
            src = cms.string('mu3'),
            leaves = cms.PSet(
                vars = cms.vstring('px', 'py',
                                   'd0_bs:dB',
                                   'ed0_bs:edB')
                ),
            Class = cms.string('pat::Muon')
            ),



            mets = cms.PSet(
            src = cms.InputTag("patMETs"),
            leaves = cms.PSet(
                vars = cms.vstring('et:et', 
                                   'phi:phi', 
                                   'ex:px', 
                                   'ey:py', 
                                   'gen_et:genMET.et', 
                                   'gen_phi:genMET.phi',
                                   'sign:metSignificance',
                                   'sumEt:sumEt', 
                                   'unCPhi:uncorrectedPhi', 
                                   'unCPt:uncorrectedPt',
                                   'et_muonCor:uncorrectedPt("uncorrJES")',
                                   'phi_muonCor:uncorrectedPhi("uncorrJES")',
                                   'et_JESCor:uncorrectedPt("uncorrMUON")',
                                   'phi_JESCor:uncorrectedPhi("uncorrMUON")'
                                   )
                ),
            Class = cms.string('pat::MET')
            ),
            
            
            #metsKT4 = cms.PSet(
            #src = cms.InputTag("patMETsKT4Calo"),
            #leaves = cms.PSet(
            #    vars = met_commonVarList
            #    ),                
            #Class = cms.string('pat::MET')
            #),

            #metsKT6 = cms.PSet(
            #src = cms.InputTag("patMETsKT6Calo"),
            #leaves = cms.PSet(
            #    vars = met_commonVarList
            #    ),
            #Class = cms.string('pat::MET')
            #),

            tcmets = cms.PSet(
            src = cms.InputTag("patMETsTC"),
            leaves = cms.PSet(
                vars = met_commonVarList
                ),
            Class = cms.string('pat::MET')
            ),         


            photons = cms.PSet(
            src = cms.InputTag("cleanPatPhotons"),#clean<=>all
                leaves = cms.PSet(
                    basicKinematicLeaves,
                    vars= cms.vstring('hadOverEM:hadronicOverEm',
                                      'scEnergy:superCluster.energy',
                                      'scRawEnergy:superCluster.rawEnergy',
                                      'scEta:superCluster.position.eta',
                                      'scPhi:superCluster.position.phi',
                                      'scEtaWidth:superCluster.etaWidth',
                                      'scPhiWidth:superCluster.phiWidth',
                                      'tIso:trackIso',
                                      'ecalIso:ecalIso',
                                      'hcalIso:hcalIso',
                                      
                                      #'isoEcalRecHit:isolationEcalRecHit',
                                      #'isoHcalRecHit:isolationHcalRecHit',
                                      #'isoSolidTrkCone:isolationSolidTrkCone',
                                      #'isoHollowTrkCone:isolationHollowTrkCone',
                                      #'nTrkSolidCone:nTrkSolidCone',
                                      #'nTrkHollowCone:nTrkHollowCone',
                                      
                                      'isoEcalRecHitDR04:ecalRecHitSumEtConeDR04',
                                      'isoHcalRecHitDR04:hcalTowerSumEtConeDR04',
                                      'isoSolidTrkConeDR04:trkSumPtSolidConeDR04',
                                      'isoHollowTrkConeDR04:trkSumPtHollowConeDR04',
                                      'nTrkSolidConeDR04:nTrkSolidConeDR04',
                                      'nTrkHollowConeDR04:nTrkSolidConeDR04',
                                      'isoEcalRecHitDR03:ecalRecHitSumEtConeDR03',
                                      'isoHcalRecHitDR03:hcalTowerSumEtConeDR03',
                                      'isoSolidTrkConeDR03:trkSumPtSolidConeDR03',
                                      'isoHollowTrkConeDR03:trkSumPtHollowConeDR03',
                                      'nTrkSolidConeDR03:nTrkSolidConeDR03',
                                      'nTrkHollowConeDR03:nTrkSolidConeDR03',
                                      

                                      #'isAlsoElectron:isAlsoElectron',
                                      'isAlsoElectron:isElectron',
                                      
                                      'hasPixelSeed:hasPixelSeed',
                                      #'isConverted:isConverted',
                                      'isConverted:isConvertedPhoton',
                                      'isEBGap:isEBGap',
                                      'isEEGap:isEEGap',
                                      'isEBEEGap:isEBEEGap',
                                      'isEBPho:isEB',#changed from isEBPho
                                      'isEEPho:isEE',#changed from isEEPho
                                      #'isLooseEM:isLooseEM',
                                      #'isLoosePhoton:isLoosePhoton',
                                      #'isTightPhoton:isTightPhoton',
                                      'isLoosePhoton:photonID("PhotonCutBasedIDLoose")',
                                      'isTightPhoton:photonID("PhotonCutBasedIDLoose")',
                                      #Currently (28-07) all photons are defined as LooseEM, but this will be added again as
                                      # they are going to relax pre-selection cuts on the photon 
                                      #'isLooseEM:photonID()',
                                      
                                      'r9:r9',
                                      'gen_et:genPhoton.et',
                                      'gen_eta:genPhoton.eta',
                                      'gen_phi:genPhoton.phi',
                                      'gen_id:genPhoton.pdgId'
                                      )
                                      ),
                    
                Class = cms.string('pat::Photon')
            ),




            tracks = cms.PSet(
                src = cms.InputTag("generalTracks"),
                leaves = cms.PSet(
                    vars = cms.vstring('chi2:chi2', 
                        'ndof:ndof', 
                        'chg:charge', 
                        'pt:pt', 
                        'px:px', 
                        'py:py', 
                        'pz:pz', 
                        'eta:eta', 
                        'phi:phi', 
                        'theta:theta', 
                        'd0dum:d0', 
                        'dz:dz', 
                        'vx:vx', 
                        'vy:vy', 
                        'vz:vz', 
                        'numvalhits:numberOfValidHits', 
                        'numlosthits:numberOfLostHits', 
                        'd0dumErr:d0Error', 
                        'dzErr:dzError', 
                        'ptErr:ptError', 
                        'etaErr:etaError', 
                        'phiErr:phiError', 
                        'Nrechits:recHitsSize', 
                        'innerHitX:innerPosition.x', 
                        'innerHitY:innerPosition.y', 
                        'innerHitZ:innerPosition.z', 
                        'outerHitX:outerPosition.x', 
                        'outerHitY:outerPosition.y', 
                        'outerHitZ:outerPosition.z',
                        'highPurity:quality("highPurity")',
                        # new addition (18-12-09): hit pattern
                        'innerLayerMissingHits:trackerExpectedHitsInner().numberOfHits()',
                        #added 1-04-10
                        'innerLayerMissingHitsPixels:trackerExpectedHitsInner().numberOfValidPixelHits()',
                        # new addition (21-1-10)
                        'numPixelHits:hitPattern.numberOfValidPixelHits()',
                        'numStripHits:hitPattern.numberOfValidStripHits()',
                        'hasValidHitInFirstPixelBarrel:hitPattern.hasValidHitInFirstPixelBarrel()'
                                       )
                ),
                Class = cms.string('reco::Track')
            ),


             els = cms.PSet(
             src = cms.string('electrons'),
             leaves = cms.PSet(
                 basicKinematicLeaves,
                 vars = ele_vars
                 ),
             Class = cms.string('pat::Electron')
             ),


            ## NEW 10 MAY 10 TL
            els2 = cms.PSet(
            src = cms.string('ele2'),
            leaves = cms.PSet(
                vars = cms.vstring('px', 'py',
                                   'd0_pvbs:dB',
                                   'ed0_pvbs:edB')
                ),
            Class = cms.string('pat::Electron')                
            ),
            
            els3 = cms.PSet(
            src = cms.string('ele3'),
            leaves = cms.PSet(
                vars = cms.vstring('px', 'py',
                                   'd0_bs:dB',
                                   'ed0_bs:edB')
                ),
            Class = cms.string('pat::Electron')                
            ),


            jets = cms.PSet(
            src = cms.string('jets'),
            leaves = cms.PSet(
                basicKinematicLeaves,
                vars = default_jet_var  ## contain also gen-info
                ),
            Class = cms.string('pat::Jet')
            ),

            ###############################################################
            ## Note: at the moment the extra jet collections do not contain
            ##       generator information
            ###############################################################
            
            #jetsKT4 = cms.PSet(
            #src = cms.string('jetsKT4'),
            #leaves = cms.PSet(
            #    basicKinematicLeaves,
            #    vars = extra_jet_var
            #    ),
            #Class = cms.string('pat::Jet')
            #),
                
            #jetsKT6 = cms.PSet(
            #src = cms.string('jetsKT6'),
            #leaves = cms.PSet(
            #    basicKinematicLeaves,
            #    vars = extra_jet_var
            #    ),
            #Class = cms.string('pat::Jet')
            #),
            
            jetsJPTAK5 = cms.PSet(
            src = cms.string('jetsAK5JPT'),
            leaves = cms.PSet(
                basicKinematicLeaves,                
                vars = extra_jet_var
                ),            
            Class = cms.string('pat::Jet')
            ),                


            PFJets2 = cms.PSet(
            src = cms.string('mypfJets'), #jetsPFAK5
            leaves = cms.PSet(
                basicKinematicLeaves,
                vars = default_jet_var   ## var as in default
                ),
            Class = cms.string('pat::Jet')
            ),
            
                                                    
            #########################
            ## Particle Flow (PF2PAT) Kachanon
            #########################
            PFJets = cms.PSet(
            src = cms.string('pfJets'), #jetsPFAK5
            leaves = cms.PSet(
                basicKinematicLeaves,
                vars = default_jet_var  ## var as in default
                ),
            Class = cms.string('pat::Jet')
            ),

            PFMets = cms.PSet(
            src = cms.string('pfMets'),
            leaves = cms.PSet(                
                vars = met_commonVarList                
                ),
            Class = cms.string('pat::MET')
            ),

            PFElsAll = cms.PSet(
            src = cms.string('pfElectronsAll'),
            leaves = cms.PSet(
                basicKinematicLeaves,
                vars = ele_vars
                ),
            Class = cms.string('pat::Electron')
            ),
            
            PFElsIso = cms.PSet(
            src = cms.string('pfElectronsIso'),
            leaves = cms.PSet(
                basicKinematicLeaves,
                vars = ele_vars
                ),
            Class = cms.string('pat::Electron')
            ),

			
			PFMuons = cms.PSet(
            src = cms.string('pfMuons'),
            leaves = cms.PSet(
                basicKinematicLeaves,
                vars = muon_vars
                ),
            Class = cms.string('pat::Muon')
            ),


            

######################################  MC Truth (begin)  ########################################

            genEv = cms.PSet(
            src = cms.InputTag("generator"),
            leaves = cms.PSet(
                vars = cms.vstring('signalProcessID',
                                   'weight',
                                   'qScale',
                                   'x1:pdf.x.first',
                                   'x2:pdf.x.second',
                                   'xpdf1:pdf.xPDF.first',
                                   'xpdf2:pdf.xPDF.second'
#                                   'id1:pdf.id.first',#get error as below:
#Expression parser error:data member named "first" for type "pair<int,int>" is not publically accessible. (char 7)
                                   )
                ),
            Class = cms.string("GenEventInfoProduct")
            ),

            mc_doc = cms.PSet(
            src = cms.InputTag("genParticles"),
            leaves = cms.PSet(
                vars = mc_doc_vars
                ),
            selection = cms.string(''),
            Class = cms.string('reco::GenParticle')
            )

############################# end of MC Truth ##########################################################


            
        ),
        
        ComponentName = cms.string('CompleteNTupler'),
        useTFileService = cms.bool(True), ## false for EDM; true for non EDM

        variablesPSet = cms.PSet(
            #use all the variables from the PSet above
            allVariables = cms.bool(True),
            treeName = cms.string('eventV')
        ),

        # TL add 16 Jun to test adhoc ntupler
        AdHocNPSet = cms.PSet(
            treeName = cms.string('eventA'),
            ecalBarrelRecHit = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
            HLTProcessName = cms.string( myHLT ) #cms.string('REDIGI')
        )
        
    )
)
