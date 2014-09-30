# Auto generated configuration file
# using: 
# Revision: 1.303 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: STEP2 --step RAW2DIGI,L1Reco,RECO,VALIDATION:validation_prod --conditions START42_V11::All --datamix NODATAMIXER --customise Configuration/GlobalRuns/reco_TLR_42X.customisePPMC --eventcontent RECOSIM --datatier GEN-SIM-RECO
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('Configuration.StandardSequences.ReMixingSeeds_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/pompyt_raw.root')
    #fileNames = cms.untracked.vstring('file:pompyt_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT.root')
fileNames = cms.untracked.vstring(
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_86_1_aMd.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_87_1_RFp.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_88_1_IZn.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_89_1_spY.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_90_1_iwH.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_91_1_pBV.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_92_1_RbK.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_93_1_PVD.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_94_1_8tU.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_95_1_wQO.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_96_1_C87.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_98_1_Rlf.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_99_1_B5F.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_9_1_Z7T.root'
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.303 $'),
    annotation = cms.untracked.string('STEP2 nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('QCD_pompyt_START42_13_RAW2DIGI_L1Reco_RECO_E.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
)

# Additional output definition

# Other statements
process.mix.playback = True
process.GlobalTag.globaltag = 'START42_V13::All'

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.validation_step = cms.EndPath(process.validation_prod)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.validation_step,process.endjob_step,process.RECOSIMoutput_step)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.GlobalRuns.reco_TLR_42X
from Configuration.GlobalRuns.reco_TLR_42X import customisePPMC 

#call to customisation function customisePPMC imported from Configuration.GlobalRuns.reco_TLR_42X
process = customisePPMC(process)

# End of customisation functions
