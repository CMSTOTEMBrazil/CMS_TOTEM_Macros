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
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_66_1_wqv.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_67_1_ORO.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_68_1_eYS.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_69_1_SpA.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_70_1_x6C.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_71_1_GL5.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_72_1_QDF.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_73_1_YyC.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_74_1_CLo.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_75_1_BcZ.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_76_1_Zl3.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_77_1_OaR.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_78_1_QWL.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_79_1_4d0.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_7_3_ctD.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_80_1_Mzl.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_81_1_DhP.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_82_1_YYy.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_83_1_Pk8.root',
'file:/storage2/dilson/diffractivebJets/pompytprod/CMSSW_4_2_8/src/GeneratorInterface/PompytInterface/test/crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_85_1_oKS.root'
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
    fileName = cms.untracked.string('QCD_pompyt_START42_13_RAW2DIGI_L1Reco_RECO_D.root'),
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
