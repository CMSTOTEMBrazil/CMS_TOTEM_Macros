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
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_100_1_bBz.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_10_1_Sgj.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_11_1_4mj.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_12_1_nu2.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_13_1_qvl.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_14_1_cL3.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_15_1_5Hg.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_16_1_ajj.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_17_1_0yG.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_18_1_5Tc.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_19_1_6tC.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_1_1_hgI.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_20_1_Yli.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_21_1_Flm.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_22_1_txX.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_23_1_ix7.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_24_1_Kfe.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_27_1_aGO.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_28_1_68j.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_29_1_zwz.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_2_1_QD6.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_30_1_yzZ.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_31_1_QfA.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_32_1_p6D.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_33_1_7pZ.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_34_1_c2p.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_35_1_xeC.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_36_1_Bq3.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_37_1_NBZ.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_38_1_aZ8.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_39_1_uVZ.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_3_3_pkH.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_40_1_ABz.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_41_1_9Zk.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_42_1_KLS.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_43_1_NmD.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_44_1_oSa.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_45_1_YjS.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_46_1_e5k.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_47_1_yJW.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_48_1_mjQ.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_4_1_6nT.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_50_3_31q.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_51_3_yGM.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_52_1_MGx.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_53_1_P24.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_54_1_1m5.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_55_1_05l.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_56_1_h7I.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_58_1_aDN.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_59_1_QIq.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_5_1_Ud6.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_60_1_TD2.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_61_1_5u3.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_62_1_lGY.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_63_1_IQV.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_64_1_aIh.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_65_1_7rM.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_66_1_wqv.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_67_1_ORO.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_68_1_eYS.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_69_1_SpA.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_70_1_x6C.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_71_1_GL5.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_72_1_QDF.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_73_1_YyC.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_74_1_CLo.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_75_1_BcZ.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_76_1_Zl3.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_77_1_OaR.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_78_1_QWL.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_79_1_4d0.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_7_3_ctD.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_80_1_Mzl.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_81_1_DhP.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_82_1_YYy.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_83_1_Pk8.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_85_1_oKS.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_86_1_aMd.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_87_1_RFp.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_88_1_IZn.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_89_1_spY.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_90_1_iwH.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_91_1_pBV.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_92_1_RbK.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_93_1_PVD.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_94_1_8tU.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_95_1_wQO.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_96_1_C87.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_98_1_Rlf.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_99_1_B5F.root',
'file:crab_0_111007_093928/res/pompyt_QCD30to_START42_V13_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_9_1_Z7T.root'

)

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
    fileName = cms.untracked.string('QCD_pompyt_START42_13_RAW2DIGI_L1Reco_RECO.root'),
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
