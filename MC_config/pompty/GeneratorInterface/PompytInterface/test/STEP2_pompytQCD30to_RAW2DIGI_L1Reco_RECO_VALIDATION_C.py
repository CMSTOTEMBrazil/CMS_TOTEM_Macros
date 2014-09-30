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
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_48_1_mjQ.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_49_2_Jga.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_4_1_6nT.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_50_3_31q.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_51_3_yGM.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_52_1_MGx.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_53_1_P24.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_54_1_1m5.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_55_1_05l.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_56_1_h7I.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_58_1_aDN.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_59_1_QIq.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_5_1_Ud6.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_60_1_TD2.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_61_1_5u3.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_62_1_lGY.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_63_1_IQV.root', 
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_63_1_IQV.root', 
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_64_1_aIh.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_65_1_7rM.root'
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
    fileName = cms.untracked.string('QCD_pompyt_START42_13_RAW2DIGI_L1Reco_RECOi_C.root'),
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
