# This is an example PAT configuration showing the usage of PAT on minbias data
#
# Things to do:
#   - Add pdf weights
#   - Check using correct global tag
#   - Add JPT collection
#   - Use correct JEC set
#   - Eventually use combined PAT+PATPF configuration once it's finished
#      (if it does what we want)
#   - Add jetIDs
#   - Add hemisphere info
#
#
#
# TL notes:
#  - change GT to V5
#  - change label fo PFJet from PFAK5 to AK5
#  - added PhysDecl filter
#  - adapted for 3_6_0
#  - added switch for running on minbias SD EG data
#  - 361: added 2ndmus, mus_pvB, mus2_d0_pvbs
#  - add PF2PAT from Kachanon
#  - add HLT report
#


# Select input type ("Data", "Data_SD", "MinBias_MC" or "Spring10_Physics_MC")
selection = "Data_SD"


# Main switch (only one should be ON)
runOn_Data                = bool(0)
runOn_Data_SD             = bool(0)
runOn_MinBias_MC          = bool(0)
runOn_Spring10_Physics_MC = bool(0)

if selection=="Data":
    runOn_Data = True
elif selection=="Data_SD":
    runOn_Data_SD = True
elif selection=="MinBias_MC":
    runOn_MinBias_MC = True
else:
    runOn_Spring10_Physics_MC = True


# switch for input: data=1, MC=0
# default value for runnning on MinBias Data
#global realData
#global myHLT
# import global variables


#from globals_data import *
realData         = True
myHLT            = "HLT"
JECSetName       = "Spring10"
wantPatTuple     = bool(0)
doPDFweights     = bool(0)
applyBSCTrigOnMC = bool(1)
runOn35xInput    = bool(0)  #<--- EG_v2 (CMSSW 361p2)
outname          = "nTuple_data.root"
patname          = "pat_data.root"


if runOn_Data_SD:
    realData         = True
    myHLT            = "HLT"
    wantPatTuple     = False
    outname          = "nTuple_data_SD.root"
    patname          = "pat_data_SD.root"

elif runOn_MinBias_MC:
    realData         = False
    wantPatTuple     = False
    doPDFweights     = False
    applyBSCTrigOnMC = True
    runOn35xInput    = True #check
    outname          = "nTuple_minbias.root"
    patname          = "pat_minbias.root"

elif runOn_Spring10_Physics_MC:
    realData         = False
    myHLT            = "REDIGI"
    wantPatTuple     = False
    doPDFweights     = False  ## No PDF weights
    #doPDFweights    = True  ##<--- Write PDF weights
    applyBSCTrigOnMC = False
    runOn35xInput    = True
    outname          = "nTuple_mc.root"
    patname          = "pat_mc.root"


print "***********************************"
if runOn_Data:
    print " Running on Data"
if runOn_Data_SD:
    print " Running on SD Data"
if runOn_MinBias_MC:
    print " Running on MinBias MC"
if runOn_Spring10_Physics_MC:
    print " Running on Spring10 physics MC"
print "***********************************"
print "Switches:"
print "  realData:         ", realData
print "  myHLT:            ", myHLT
print "  JECSetName:       ", JECSetName
print "  wantPatTuple:     ", wantPatTuple
print "  doPDFweights:     ", doPDFweights
print "  applyBSCTrigOnMC: ", applyBSCTrigOnMC
print "  runOn35xInput:    ", runOn35xInput
print "***********************************"

if not realData and applyBSCTrigOnMC:
    print "***********************************"
    print " Running on MC, take note that following L1 Tech Trg BSC bits are required!!!"
    print "   (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)"
    print "***********************************"





# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *



## PF2PAT Was Here.

######################### Kachanon adds PF2PAT #########################
from PhysicsTools.PatAlgos.tools.pfTools import *
postfix = "PF"

if realData:
    usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=False, postfix=postfix)
    # turn to false when running on data
    getattr(process, "patElectrons"+postfix).embedGenMatch = False
    getattr(process, "patMuons"+postfix).embedGenMatch = False
    getattr(process, "patJets"+postfix).embedGenPartonMatch = False
    getattr(process, "patJets"+postfix).embedGenJetMatch = False
else:
    usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=True, postfix=postfix)
    getattr(process, "patElectrons"+postfix).embedGenMatch = True
    getattr(process, "patMuons"+postfix).embedGenMatch = True

# To tell PF2PAT to produce pf all electrons in stead of pf isolated electrons
#process.pfIsolatedElectrons.combinedIsolationCut = cms.double(10000)
getattr(process, "pfIsolatedElectrons"+postfix).combinedIsolationCut = cms.double(10000)
########################################################################





###################
#  Add JPT jets
###################
# Ref: https://twiki.cern.ch/twiki/bin/view/CMS/JetPlusTracksCorrections
process.load("RecoJets.Configuration.RecoJPTJets_cff")

# Jet Correction (in 3_6_0)
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

# to run on 35X input sample
# ref: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePATRecipes#CMSSW_3_6_X
from PhysicsTools.PatAlgos.tools.cmsswVersionTools import *
if runOn35xInput:
    run36xOn35xInput( process )




#process.GlobalTag.globaltag =   cms.string('MC_3XY_V18::All')
#process.GlobalTag.globaltag =   cms.string('MC_3XY_V25::All')
#process.GlobalTag.globaltag =   cms.string('START3X_V25B::All')
#process.GlobalTag.globaltag =   cms.string('GR10_P_V4::All')
process.GlobalTag.globaltag =   cms.string('GR10_P_V5::All') #TL, 16 Apr
#START3X_V20 



from PhysicsTools.PatAlgos.tools.coreTools import *

#switch off new tau features introduced in 33X to restore 31X defaults
# new feaures: - shrinkingConeTaus instead of fixedCone ones
#              - TaNC discriminants attached for shrinkingConeTaus
#              - default preselection on cleaningLayer1
from PhysicsTools.PatAlgos.tools.tauTools import *
#switchTo31Xdefaults(process)


# turn off MC matching for the process
if realData:
    print "*********************"
    print "Turn off MC matching"
    print "*********************"
    removeMCMatching(process, ["All"])




# set jet corrections
#print "*******************************************************"
#print "Calling switchJECSet() to set jet energy corrections: ", JECSetName
#print "*******************************************************"
#from PhysicsTools.PatAlgos.tools.jetTools import switchJECSet
#switchJECSet( process, JECSetName )

from PhysicsTools.PatAlgos.tools.jetTools import *

# remove the tag infos
# need to keep tag infos to access SV info via pat::jet->tagInfoSecondaryVertex()
process.patJets.addTagInfos = True  #D=False
# require jet pt > 10 (L2+L3 corrected)
process.selectedPatJets.cut = cms.string('pt > 10')
# look for at least one jet
process.countPatJets.minNumber = 0












### 6-1-10 ###
## TL added 6-1-10: add new jet collections
## Refer to PhysicsTools/PatAlgos/test/patLayer1_fromAOD_jetSuite_full_cfg.py
## The following block needs to be after 'removeMCMatching', otherwise will get error of 'genMetCalo not found'.
# uncomment the following line to add tcMET to the event content
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
##addPfMET(process, 'PF')  ##<-- gives error of genMetCalo not found if put before removeMCMatching

# produce jpt corrected calo jets, which are not on AOD per default
#process.load("PhysicsTools.PatAlgos.recoLayer0.jetPlusTrack_cff")



print "*******************************"
print "Calling addJetCollection()"
print " Note: b-tagging and SSV info are not stored for the extra jet collections"
print "*******************************"





# add ak5GenJets as they are not in the Spring10 physics MC
# ref: https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/899/1/1/1/1.html

# ak5GenJets are NOT there: First load the needed modules
if runOn_Spring10_Physics_MC:
    process.load("RecoJets.Configuration.GenJetParticles_cff")
    process.load("RecoJets.JetProducers.ak5GenJets_cfi")
# process.p = cms.Path(process.genParticlesForJets *
#                      process.ak5GenJets )
                        

####################################
## JEC for the defaul jet collection
####################################
process.patJetCorrFactors.corrSample = JECSetName


####################################
##
##    Add extra jet collections
##
####################################
# JPT jets
# when using 36X input, doesn't seem to make any difference with addJetCollection
# or addJetCollection35X
#
addJetCollection(process,cms.InputTag('JetPlusTrackZSPCorJetAntiKt5'),
#addJetCollection35X(process,cms.InputTag('JetPlusTrackZSPCorJetAntiKt5'),
                 'AK5', 'JPT',
                 doJTA        = True,
                 doBTagging   = False, #off
                 jetCorrLabel = ('AK5','JPT'), #None,
                 doType1MET   = False,
                 doL1Cleaning = False,
                 doL1Counters = False,
                 genJetCollection = cms.InputTag("ak5GenJets"),
                 doJetID      = True,
                 jetIdLabel   = "ak5"
                 )

# kt4Calo jets
addJetCollection(process,cms.InputTag('kt4CaloJets'),
#addJetCollection35X(process,cms.InputTag('kt4CaloJets'),
                 'KT4','Calo',
                 doJTA        = True,
                 doBTagging   = False, #off
                 jetCorrLabel = ('KT4','Calo'),
                 doType1MET   = True,
                 doL1Cleaning = True,
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("kt4GenJets"),
                 doJetID      = True,
                 jetIdLabel   = "kt4"
                 )
# kt6Calo jets
addJetCollection(process,cms.InputTag('kt6CaloJets'),
#addJetCollection35X(process,cms.InputTag('kt6CaloJets'),
                 'KT6','Calo',
                 doJTA        = True,
                 doBTagging   = False, #off
                 jetCorrLabel = ('KT6','Calo'),
                 doType1MET   = True,
                 doL1Cleaning = True,
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("kt6GenJets"),
                 doJetID      = True,
                 jetIdLabel   = "kt6"
                 )

# ak5 Paticle Flow jets
#addJetCollection(process,cms.InputTag('ak5PFJets'),
#addJetCollection35X(process,cms.InputTag('ak5PFJets'),
#                 'AK5','PF',
#                 doJTA        = True,
#                 doBTagging   = True, #ON
#                 jetCorrLabel =  ('AK5','PF'), #test
#                 doType1MET   = False,
#                 doL1Cleaning = True,
#                 doL1Counters = False,
#                 genJetCollection=cms.InputTag("ak5GenJets"),
#                 doJetID      = False
#                 )

# no btag for extra jets
process.patJetsKT4Calo.addTagInfos = False
process.patJetsKT6Calo.addTagInfos = False
process.patJetsAK5JPT.addTagInfos  = False
#process.patJetsAK5PF.addTagInfos   = True
process.patJetsPF.addTagInfos   = True


################################
##  Jet Energy Correction
################################

process.patJetCorrFactorsKT4Calo.corrSample = JECSetName
process.patJetCorrFactorsKT6Calo.corrSample = JECSetName
process.patJetCorrFactorsAK5JPT.corrSample  = JECSetName
#process.patJetCorrFactorsAK5PF.corrSample   = JECSetName
process.patJetCorrFactorsPF.corrSample   = JECSetName



# Input files
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()



# test real data (copied from PatExample)
readFiles.extend( [
    #'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/F4C92A98-163C-DF11-9788-0030487C7392.root'
    # one of the /MinimumBias/Commissioning10-SD_EG-v9/RECO
    #LFN: /store/data/Commissioning10/MinimumBias/RECO/v9/000/133/887/1CD67BCA-EA50-DF11-A8EC-00E0817918A7.root 
    #'rfio:/castor/cern.ch/user/t/tcheng/1CD67BCA-EA50-DF11-A8EC-00E0817918A7.root'
    # test using 361p2 data file (/EG/Run2010A-PromptReco-v2/RECO)
    'rfio:/castor/cern.ch/cms/store/data/Run2010A/EG/RECO/v2/000/136/089/BE05CF8D-6E67-DF11-B5D9-000423D99996.root'
    ] );


if realData:
    process.source.fileNames = readFiles
else:
    if runOn_MinBias_MC:
        ## test using Spring10 MinBias MC
        process.source.fileNames = cms.untracked.vstring(
            '/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V26A_356ReReco-v1/0007/00536B4F-E23D-DF11-A2CF-003048678B26.root'
            )
    else:
        ## test using Spring10 ttjet: /TTbarJets-madgraph/Spring10_START3X_V26_S09-v1/GEN-SIM-RECO (need run33xOnReReco)
        process.source.fileNames = cms.untracked.vstring(
            'rfio:/castor/cern.ch/user/t/tcheng/A4121AB4-0747-DF11-8984-0030487F171B.root'
        )

            

# configure HLT
if realData:
    if not runOn_Data_SD:
        process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
        process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
        process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
        process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')
        
else:
    #    if applyBSCTrigOnMC:
    process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
    process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
    process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
    #process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')
    process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')
    




## TL 16 Apr
## - Added physics declared, scraping filter, pv filter
## - following the example in PhysicsTools/PatExamples/test/patLayer1_fromRECO_7TeV_firstdata_cfg.py
# require physics declared
#process.hltPhysicsDeclared = cms.EDFilter("PhysDecl",
#   applyfilter = cms.untracked.bool(True),
#   HLTriggerResults = cms.InputTag("TriggerResults","","HLT")
#)
# Update for 357
if runOn_Data and not runOn_Data_SD:
    process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
    process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'

    
# HLT Trigger Report
# process.hlTrigReport
#myHLT="HLT" #"REDIGI"
process.hlTrigReport = cms.EDAnalyzer("HLTrigReport",
    HLTriggerResults = cms.InputTag("TriggerResults","",myHLT)
)




# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
#process.MessageLogger.cerr.INFO.limit = 10
#process.MessageLogger.cerr.threshold = "DEBUG"
process.MessageLogger.categories.append("HLTrigReport")
#process.MessageLogger.hlTrigReport.limit = 1000




# process all the events
process.maxEvents.input = 10 #20000
process.options.wantSummary = True

#process.out.outputCommands += (['keep *_*_*_*'
#                               ])


#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('pat_dataPATLayer1_Output.fromAOD_full.root'),
    outputCommands = cms.untracked.vstring('drop *',
                                           #*patEventContent)
#                                           'keep *_cleanPatElectrons_*_*',
                                           'keep *_cleanPat*_*_*',
#                                           'keep *_cleanPatEle2_*_*',
#                                           'keep *_cleanPatEle3_*_*',
                                           'keep *_cleanPatMuons_*_*',
                                           'keep *_cleanPatMu2_*_*',
                                           'keep *_cleanPatMu3_*_*',
                                           'keep *_cleanPatJets_*_*',
                                           'keep *_cleanPatJetsAK5PF_*_*',
                                           'keep *_cleanPatJetsAK5JPT_*_*',
                                           'keep *_patMETs_*_*',
                                           'keep *_patMETsPF_*_*',
                                           'keep *_patMETsTC_*_*',
                                           'keep *_cleanPatPhotons_*_*',
                                           'keep *_cleanPatTaus_*_*',
                                           #'keep recoSuperCluster_*_*_*',
                                           'keep *_towerMaker_*_*',
                                           'keep recoCaloTower_*_*_*',
                                           'keep recoTracks_generalTracks_*_*',
                                           'keep edmTriggerResults_TriggerResults_*_*',
                                           'keep *_offlinePrimaryVertices_*_*',
                                           'keep *_offlinePrimaryVerticesWithBS_*_*',
                                           'keep *_offlineBeamSpot_*_*'
                                           #,'keep *'
                                           )
    )

process.out.fileName = patname
################# Above are examples from PhysicsTools ##############################




# add PDFWeightProducer (8-2-2010)

if doPDFweights:
    process.pdfWeights = cms.EDProducer(
        "PdfWeightProducer",
        PdfInfoTag = cms.untracked.InputTag("generator"),
        PdfSetNames = cms.untracked.vstring("cteq66.LHgrid"
                                            #, "MRST2006nnlo.LHgrid"
                                            #, "MRST2007lomod.LHgrid"
                                            )
        )



#####################################################################
#
# Add second electron collections for d0 w.r.t beam-constrained PV
# TL: 10 May 2010
#
#####################################################################
#process.patElectrons.pvSrc = "offlinePrimaryVerticesWithBS"

import PhysicsTools.PatAlgos.producersLayer1.electronProducer_cfi as eleProducer
process.patEle2 = eleProducer.patElectrons.clone(
    pvSrc = "offlinePrimaryVerticesWithBS"
    )

process.patEle3 = eleProducer.patElectrons.clone(
    usePV = False,
    beamLineSrc = "offlineBeamSpot"
    )

    
import PhysicsTools.PatAlgos.cleaningLayer1.electronCleaner_cfi as eleCleaner
process.cleanPatEle2 = eleCleaner.cleanPatElectrons.clone(
    src = "patEle2"
    )
process.cleanPatEle3 = eleCleaner.cleanPatElectrons.clone(
    src = "patEle3"
    )

process.patElectronsPF.pvSrc = "offlinePrimaryVertices"


################
#  muons d0
################
import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi as muProducer
process.patMu2 = muProducer.patMuons.clone(
    pvSrc = "offlinePrimaryVerticesWithBS"
    )
process.patMu3 = muProducer.patMuons.clone(
    usePV = False,
    beamLineSrc = "offlineBeamSpot"
    )

    
import PhysicsTools.PatAlgos.cleaningLayer1.muonCleaner_cfi as muCleaner
process.cleanPatMu2 = muCleaner.cleanPatMuons.clone(
    src = "patMu2"
    )
process.cleanPatMu3 = muCleaner.cleanPatMuons.clone(
    src = "patMu3"
    )


process.myExtraLepton = cms.Sequence(
    process.patEle2 *
    process.patEle3 *
    process.patMu2 *
    process.patMu3 *    
    process.cleanPatEle2 *
    process.cleanPatEle3 *
    process.cleanPatMu2 *
    process.cleanPatMu3    
    )




# Add ntupler
#print "----------------------------------------------------------------------"
#print "NOTE: The HLT list in the ntuple may need to be reviewed when running"
#print "      on data samples later than Summer09 (312)."
#print "      Use HLTrigReport to find out available HLT paths."
#print "----------------------------------------------------------------------"


# note: need to take out genParticles in cfi as it throws a lot of error messages when running on data.
if realData:
    process.load("Workspace.ConfigurableAnalysis.361.configurableAnalysis_data_cfi")
else:
    if runOn_Spring10_Physics_MC:
        process.load("Workspace.ConfigurableAnalysis.361.configurableAnalysis_mc_redigi_cfi")
    else:
        process.load("Workspace.ConfigurableAnalysis.361.configurableAnalysis_mc_cfi")

#process.MessageLogger.destinations += cms.string("configurableAnalysis")
#process.MessageLogger.configurableAnalysis.threshold = 'INFO'

# extraContent = [
#     'keep *_generalTracks_*_*',
#     'keep recoPhotonCores_*_*_*',
#     'keep recoCaloClusters_*_*_*',
#     'keep *_iterativeCone5PFJets_*_*',
#     'keep *_pfMet_*_*',
#     'keep *_layer1METsTC_*_*',
#     "keep *_cleanLayer1Jets*_*_*",
#     "keep *_selectedLayer1Jets*_*_*",
#     "keep *_layer1METs*_*_*"
# ]

process.TFileService.fileName = outname

#to prevent patTuple output, delete the outpath:
if not wantPatTuple:
    del process.outpath



process.myJPT = cms.Sequence(
    process.recoJPTJets *
    process.ak5JPTJetsL2L3 
    )
# If running on 36X-reco, then no need to do recoJPTJets
if not runOn35xInput:
    process.myJPT.remove(process.recoJPTJets)


    
# let it run
if runOn_Data: 
    process.p = cms.Path(
        process.hltLevel1GTSeed*
        process.hltPhysicsDeclared*
        process.hlTrigReport *
        process.recoJPTJets * #jpt in 3_6
        #process.jptCaloJets* #jpt in 3_5
        process.patDefaultSequence *
        process.myExtraLepton
        + getattr(process,"patPF2PATSequence"+postfix) *
        process.configurableAnalysis
        )
    
elif runOn_Data_SD: 
    process.p = cms.Path(
        process.hlTrigReport *
        process.myJPT *
        process.patDefaultSequence *
        process.myExtraLepton
        + getattr(process,"patPF2PATSequence"+postfix) *
        process.configurableAnalysis
        )

elif runOn_MinBias_MC:
    
    process.p = cms.Path(
        process.hltLevel1GTSeed *
        process.hlTrigReport *
        process.myJPT *
        process.patDefaultSequence *
        process.myExtraLepton
        + getattr(process,"patPF2PATSequence"+postfix)
        )

    if not applyBSCTrigOnMC:
        process.p.remove( process.hltLevel1GTSeed )
        
    process.p *= process.configurableAnalysis

        
elif runOn_Spring10_Physics_MC:
    
    process.p = cms.Path(
        process.hlTrigReport *
        process.genParticlesForJets *
        process.ak5GenJets *
        process.myJPT *
        process.patDefaultSequence *
        process.myExtraLepton
        + getattr(process,"patPF2PATSequence"+postfix)
        )
        
    ## add PDF weighs as required
    if doPDFweights:
        process.p *= process.pdfWeights
        
    process.p *= process.configurableAnalysis

