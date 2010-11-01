# Things that have changed from 06 Oct '10:
# 
#   process.patJets.addTagInfos = False and process.patJets.addBTagInfo = true gives the ssvHE values, for rereco AND aod
# - NOTE: PFmet seems to be the same added via Kachanon's lines, or using addPfMET(process, 'PF')
#
# This is an example PAT configuration showing the usage of PAT on minbias data
#
# More things to do
# - Residual JEC not yet availble in 38X (04-10-10). Available later in the week perhaps (Kostas)
# - Jet id for PF, and electron d0
# - apply appropriate jec sets 
# - add triggers
# - in 38X, the ecal EE corrections need not be applied as it is done in reco,
#   so the following sequences are no longer needed:
#     process.gsfElectrons *
#     process.eIdSequence * 
#
# For 38X, we need / don't need:
#   - Do NOT need spike cleaning (can remove spike cleaning from adhocnTuple)
#   - Do NOT need to apply ecal EE corrections
#   - DO need to apply HBHE noise filter
#   - DO need to apply scraping filter
#
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
# Change Data_SD to Data_EG
#

import sys
print sys.argv

# Select input type ("Data", "Data_EG", "MinBias_MC" or "Spring10_Physics_MC")
selection = "Data_EG"
#selection = "Spring10_Physics_MC"


# Main switch (only one should be ON)
runOn_Data                = bool(0)
runOn_Data_EG             = bool(0)
runOn_MinBias_MC          = bool(0)
runOn_Spring10_Physics_MC = bool(0)

if selection=="Data":
    runOn_Data = True
elif selection=="Data_EG":
    runOn_Data_EG = True
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
runOn35xInput    = bool(0)
outname          = "nTuple_data.root"
patname          = "pat_data.root"


if runOn_Data_EG:
    realData         = True
    myHLT            = "HLT"
    wantPatTuple     = False
    outname          = "nTuple_data_EG.root"
    patname          = "pat_data_EG.root"

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
if runOn_Data_EG:
    print " Running on Data EG skim"
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






######################### Kachanon adds PF2PAT #########################

#from PhysicsTools.PatAlgos.tools.pfTools import *
from Workspace.ConfigurableAnalysis.pfTools import *

postfix = "PF"
postfix_Toto = "Toto"

cloneProcessingSnippet(process, process.makePatElectrons, postfix_Toto)
getattr(process, "patElectrons"+postfix_Toto).pfElectronSource = cms.InputTag("pfIsolatedElectronsTotoPF")

# for MC only
if not realData: 
    getattr(process, "patElectrons"+postfix_Toto).genParticleMatch = cms.InputTag("electronMatch"+postfix_Toto)

process.selectedPatElectronsToto = process.selectedPatElectrons.clone()
getattr(process, "selectedPatElectrons"+postfix_Toto).src = cms.InputTag("patElectrons"+postfix_Toto)

getattr(process, "patElectrons"+postfix_Toto).useParticleFlow = True
getattr(process, "patElectrons"+postfix_Toto).userIsolation   = cms.PSet()
getattr(process, "patElectrons"+postfix_Toto).isoDeposits = cms.PSet(
    pfChargedHadrons = cms.InputTag("isoDepElectronWithCharged" + postfix),
    pfNeutralHadrons = cms.InputTag("isoDepElectronWithNeutral" + postfix),
    pfPhotons = cms.InputTag("isoDepElectronWithPhotons" + postfix)
    )
getattr(process, "patElectrons"+postfix_Toto).isolationValues = cms.PSet(
    pfChargedHadrons = cms.InputTag("isoValElectronWithCharged" + postfix),
    pfNeutralHadrons = cms.InputTag("isoValElectronWithNeutral" + postfix),
    pfPhotons = cms.InputTag("isoValElectronWithPhotons" + postfix)
    )

#For this version of PatAlgos, patElectronIsolation is commented out
#getattr(process, "makePatElectrons", postfix_Toto).remove(getattr(process, "patElectronIsolation"+postfix_Toto))

if realData:
    usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=False, postfix=postfix)
    # turn to false when running on data
    getattr(process, "patElectrons"+postfix).embedGenMatch = False
    getattr(process, "patElectrons"+postfix_Toto).addGenMatch = False #TL add 9 Jun
    getattr(process, "patElectrons"+postfix_Toto).embedGenMatch = False #TL add 9 Jun   
    getattr(process, "patMuons"+postfix).embedGenMatch = False
    getattr(process, "patJets"+postfix).embedGenPartonMatch = False
    getattr(process, "patJets"+postfix).embedGenJetMatch = False
else:
    usePF2PAT(process,runPF2PAT=True, jetAlgo='AK5', runOnMC=True, postfix=postfix)
    getattr(process, "patElectrons"+postfix).embedGenMatch = True
    getattr(process, "patMuons"+postfix).embedGenMatch = True


process.morePFElectron = cms.Sequence(
		#getattr(process, "pfMET"+postfix) *
		getattr(process, "pfNoPileUpSequence"+postfix) *
		getattr(process, "pfAllNeutralHadrons"+postfix) *
		getattr(process, "pfAllChargedHadrons"+postfix) *
		getattr(process, "pfAllPhotons"+postfix) *
		getattr(process, "pfMuonSequence"+postfix) *
		getattr(process, "pfMuonSequence"+postfix) *
		getattr(process, "pfNoMuon"+postfix) *
		getattr(process, "pfElectronSequence"+postfix) *
		getattr(process, "makePatElectrons"+postfix_Toto) *
		getattr(process, "selectedPatElectrons"+postfix_Toto)		
)

# remove electronMatchToto
if realData:
    process.makePatElectronsToto.remove(process.electronMatchToto)
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
#if runOn35xInput:
#    run36xOn35xInput( process )




#process.GlobalTag.globaltag =   cms.string('MC_3XY_V18::All')
#process.GlobalTag.globaltag =   cms.string('MC_3XY_V25::All')
#process.GlobalTag.globaltag =   cms.string('START3X_V25B::All')
#process.GlobalTag.globaltag =   cms.string('GR10_P_V4::All')
#process.GlobalTag.globaltag =   cms.string('GR10_P_V5::All') #TL, 16 Apr
#process.GlobalTag.globaltag =   cms.string('GR10_P_V10::All') #TL, 16 Apr
process.GlobalTag.globaltag =   cms.string('GR_R_38X_V13::All') #

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
process.patJets.addTagInfos = False
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
#addPfMET(process, 'PF')  ##<-- gives error of genMetCalo not found if put before removeMCMatching

# produce jpt corrected calo jets, which are not on AOD per default
#process.load("PhysicsTools.PatAlgos.recoLayer0.jetPlusTrack_cff")



print "*******************************"
print "Calling addJetCollection()"
print " Note: b-tagging and SSV info are not stored for the extra jet collections"
print "*******************************"

print process.patJetCorrFactors.corrSample



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
#process.patJetCorrFactors.corrSample = JECSetName
#The default in 382 appears to be spring10 

####################################
##
##    Add extra jet collections
##
####################################
# JPT jets
#addJetCollection(process,cms.InputTag('JetPlusTrackZSPCorJetAntiKt5'),
addJetCollection(process,cms.InputTag('JetPlusTrackZSPCorJetAntiKt5'),
                 'AK5', 'JPT',
                 doJTA        = True,
                 doBTagging   = True, #off
                 jetCorrLabel = ('AK5','JPT'), #None,
                 doType1MET   = False,
                 doL1Cleaning = False,
                 doL1Counters = False,
                 genJetCollection = cms.InputTag("ak5GenJets"),
                 doJetID      = True,
                 jetIdLabel   = "ak5"
                 )

print process.patJetCorrFactors.corrSample

# kt4Calo jets
#addJetCollection(process,cms.InputTag('kt4CaloJets'),
#addJetCollection(process,cms.InputTag('kt4CaloJets'),
#                 'KT4','Calo',
#                 doJTA        = True,
#                 doBTagging   = False, #off
#                 jetCorrLabel = ('KT4','Calo'),
#                 doType1MET   = True,
#                 doL1Cleaning = True,
#                 doL1Counters = False,
#                 genJetCollection=cms.InputTag("kt4GenJets"),
#                 doJetID      = True,
#                 jetIdLabel   = "kt4"
#                 )
# kt6Calo jets
#addJetCollection(process,cms.InputTag('kt6CaloJets'),
#addJetCollection(process,cms.InputTag('kt6CaloJets'),
#                 'KT6','Calo',
#                 doJTA        = True,
#                 doBTagging   = False, #off
#                 jetCorrLabel = ('KT6','Calo'),
#                 doType1MET   = True,
#                 doL1Cleaning = True,
#                 doL1Counters = False,
#                 genJetCollection=cms.InputTag("kt6GenJets"),
#                 doJetID      = True,
#                 jetIdLabel   = "kt6"
#                 )

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
#process.patJetsKT4Calo.addTagInfos = False
#process.patJetsKT6Calo.addTagInfos = False
process.patJetsAK5JPT.addTagInfos  =  False
#process.patJetsAK5PF.addTagInfos   = True
process.patJetsPF.addTagInfos   = False#True

#
#process.patJets.addBTagInfo = False

################################
##  Jet Energy Correction
################################

#process.patJetCorrFactorsKT4Calo.corrSample = JECSetName
#process.patJetCorrFactorsKT6Calo.corrSample = JECSetName
#process.patJetCorrFactorsAK5JPT.corrSample  = JECSetName
###process.patJetCorrFactorsAK5PF.corrSample   = JECSetName
#process.patJetCorrFactorsPF.corrSample   = JECSetName



# Input files
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

#This produces the collection cleanPatJetsAK5PF
process.load("RecoJets.JetProducers.ak5PFJets_cfi")
addJetCollection(process, cms.InputTag('ak5PFJets::PAT'), 'AK5', 'PF', jetCorrLabel=('AK5','PF'), doType1MET=False, doJetID = False)
# make sure to run process.ak5PFJets before PAT, for example:
process.patDefaultSequence = cms.Sequence(process.ak5PFJets * process.patDefaultSequence)


# test real data (copied from PatExample)
readFiles.extend( [
    #'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/440/F4C92A98-163C-DF11-9788-0030487C7392.root'
    # one of the /MinimumBias/Commissioning10-SD_EG-v9/RECO
    #LFN: /store/data/Commissioning10/MinimumBias/RECO/v9/000/133/887/1CD67BCA-EA50-DF11-A8EC-00E0817918A7.root 
    #'rfio:/castor/cern.ch/user/t/tcheng/1CD67BCA-EA50-DF11-A8EC-00E0817918A7.root'
    # sample file from: /EG/Run2010A-May27thReReco_v1/RECO 
    #    /store/data/Run2010A/EG/RECO/May27thReReco_v1/0174/FCA740BC-956A-DF11-A982-003048F0E186.root
    ##'rfio:/castor/cern.ch/user/b/bostock/24947446-5B87-DF11-A504-0030487CAF5E.root'
    #'rfio:/castor/cern.ch/user/b/bostock/reco_run_139195_lumi_77_ev_69244083.root'
    #'rfio:/castor/cern.ch/cms/store/data/Run2010A/EG/RECO/v2/000/136/100/64D2B03E-0468-DF11-94F9-000423D99996.root'
    #'rfio:/castor/cern.ch/user/t/tcheng/FCA740BC-956A-DF11-A982-003048F0E186.root'
    #'file:/storage/top/FED5D572-8DBA-DF11-B636-0030487FA609.root' ##re-reco Sep3 sample
    #'file:/storage/top/6625EA96-45C7-DF11-860F-001D09F24F65.root' #2010B AOD sample
    'file:/storage/top/mc/045D512E-FFCF-DF11-BDC3-003048F117B6.root'#2010B sample
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
            #'rfio:/castor/cern.ch/user/t/tcheng/A4121AB4-0747-DF11-8984-0030487F171B.root'
        )

            

# configure HLT
if realData:
    process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
    process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')
    process.hltLevel1GTSeed.L1TechTriggerSeeding = cms.bool(True)
    process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0')
    #process.hltLevel1GTSeed.L1SeedsLogicalExpression = cms.string('0 AND (40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)')
    
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
if runOn_Data and not runOn_Data_EG:
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
    outputCommands = cms.untracked.vstring('keep *',#'drop *',
                                           #*patEventContent)
                                           'keep *_patElectrons_*_*',
                                           'keep *_patElectronsToto_*_*',
                                           'keep *_patElectronsPF_*_*',
                                           'keep *_selectedPatElectronsToto_*_*',
                                           'keep *_selectedPatElectronsPF_*_*',
                                           'keep *_selectedPatMuonsPF_*_*',
                                           'keep *_selectedPatJetsPF_*_*',
                                           
                                           'keep *_cleanPatElectrons_*_*',
                                           'keep *_cleanPatEle2_*_*',
                                           'keep *_cleanPatEle3_*_*',
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
                                           'keep recoSuperCluster_*_*_*',
                                           'keep recoTracks_generalTracks_*_*',
                                           'keep edmTriggerResults_TriggerResults_*_*',
                                           'keep *_offlinePrimaryVertices_*_*',
                                           'keep *_offlinePrimaryVerticesWithBS_*_*',
                                           'keep *_offlineBeamSpot_*_*'
										   
                                           )
    )

process.out.fileName = patname
################# Above are examples from PhysicsTools ##############################

#commented next line out as not available (yet) in 38X 
#process.load("RecoEgamma.EgammaTools.correctedElectronsProducer_cfi")
print "WARNING"
print "Electrons have not been corrected for alignment"

#from RecoEgamma.EgammaTools.correctedElectronsProducer_cfi import *
#process.patElectrons.electronSource = "gsfElectrons::PAT"



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
		process.cleanPatMu3 *
		process.morePFElectron

)




# Add ntupler
#print "----------------------------------------------------------------------"
#print "NOTE: The HLT list in the ntuple may need to be reviewed when running"
#print "      on data samples later than Summer09 (312)."
#print "      Use HLTrigReport to find out available HLT paths."
#print "----------------------------------------------------------------------"


# note: need to take out genParticles in cfi as it throws a lot of error messages when running on data.
print "loading CFA config"
if realData:
    process.load("Workspace.ConfigurableAnalysis.362.configurableAnalysis_data_cfi")
else:
    if runOn_Spring10_Physics_MC:
        process.load("Workspace.ConfigurableAnalysis.362.configurableAnalysis_mc_redigi_cfi")
    else:
        process.load("Workspace.ConfigurableAnalysis.362.configurableAnalysis_mc_cfi")

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


process.load("CommonTools.RecoAlgos.HBHENoiseFilter_cfi")
#process.p = cms.Path(process.HBHENoiseFilter)

process.TFileService.fileName = 'nTuple_data.root'
#process.TFileService.fileName = outname

#to prevent patTuple output, delete the outpath:
if not wantPatTuple:
    del process.outpath



#process.myJPT = cms.Sequence(
#    process.recoJPTJets *
#    process.ak5JPTJetsL2L3 
#    )
# If running on 36X-reco, then no need to do recoJPTJets
#if not runOn35xInput:
#    process.myJPT.remove(process.recoJPTJets)


    
# let it run
if runOn_Data: 
    process.p = cms.Path(
        process.hltLevel1GTSeed*
        process.hltPhysicsDeclared*
        process.hlTrigReport *
        process.recoJPTJets * #jpt in 3_6
        #process.jptCaloJets* #jpt in 3_5
        process.patDefaultSequence *
        getattr(process,"patPF2PATSequence"+postfix) *
        process.myExtraLepton *
		process.configurableAnalysis
        )
    
elif runOn_Data_EG:
    process.p = cms.Path(
        process.HBHENoiseFilter *
        process.hlTrigReport *
        #process.myJPT *
        #process.gsfElectrons *
        #process.eIdSequence *
        process.patDefaultSequence *
        getattr(process,"patPF2PATSequence"+postfix) *
        process.myExtraLepton 
        *process.configurableAnalysis
        )


elif runOn_MinBias_MC:
    
    process.p = cms.Path(
        process.hltLevel1GTSeed *
        process.hlTrigReport *
        #process.myJPT *
        process.patDefaultSequence *
        process.myExtraLepton *
        getattr(process,"patPF2PATSequence"+postfix)
        )

    if not applyBSCTrigOnMC:
        process.p.remove( process.hltLevel1GTSeed )
        
    process.p *= process.configurableAnalysis

        
elif runOn_Spring10_Physics_MC:
    
    process.p = cms.Path(
        process.hlTrigReport *
        process.genParticlesForJets *
        process.ak5GenJets *
        #process.myJPT *
        process.patDefaultSequence *
        process.myExtraLepton *
        getattr(process,"patPF2PATSequence"+postfix)
        )
        
    ## add PDF weighs as required
    if doPDFweights:
        process.p *= process.pdfWeights
        
    process.p *= process.configurableAnalysis

