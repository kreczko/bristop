import FWCore.ParameterSet.Config as cms

# import the PatTuple sequence
from Workspace.ConfigurableAnalysis.tlkn_361_testFromRECO_data_SD_cfg import *

#process.GlobalTag.globaltag = 'GR10_P_V5::All'
# Data reprocessing global tags >= 361
process.GlobalTag.globaltag = 'GR_R_36X_V10::All'

# reduce stdout
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
#process.MessageLogger.cerr.INFO.limit = 10
#process.MessageLogger.cerr.threshold = 'ERROR'

# rename output file
process.TFileService.fileName = 'nTuple_data.root'
process.out.fileName          = 'pat_data.root'
