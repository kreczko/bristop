import FWCore.ParameterSet.Config as cms

# import the PatTuple sequence
from Workspace.ConfigurableAnalysis.tlkn_361_testFromRECO_mc_spring10_cfg import *

# use GT for MC as in the RECO
process.GlobalTag.globaltag = 'START3X_V26::All'

# reduce stdout
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# rename output file
process.TFileService.fileName = 'nTuple_ttjet.root'
