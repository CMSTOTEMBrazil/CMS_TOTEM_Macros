#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms

process = cms.Process("GEN")

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")

from Configuration.Generator.PythiaUESettings_cfi import *

process.source = cms.Source("EmptySource")

process.generator = cms.EDProducer("PompytProducer",
    pythiaPylistVerbosity = cms.untracked.int32(1),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(1),
    comEnergy = cms.untracked.double(7000.0),
##  pomeron direction
    pomZ = cms.untracked.double(1.0),
##  pomPdiss = 1 (pp->Xp),  2( pp->XY)
    pomPdiss = cms.untracked.int32(1),
##  xi max
    xiMax = cms.untracked.double(0.2),
####
    pdiss_B0 = cms.untracked.double(1.5),
    pdiss_AP = cms.untracked.double(0.016),
    pdiss_GLU = cms.untracked.double(0.5),
    pdiss_MMAX = cms.untracked.double(2000.0),
    PythiaParameters = cms.PSet(
        pythiaUESettingsBlock,
        pytParameters = cms.vstring('MSEL=1', 
            'CKIN(3)=5', 
            'CKIN(4)=15'),
        parameterSets = cms.vstring('pythiaUESettings', 
##        parameterSets = cms.vstring(
            'pytParameters')
    )
)

process.RandomNumberGeneratorService.generator.initialSeed = cms.untracked.uint32(123456788)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")

process.caloJetMetGen=cms.Sequence(process.genParticles + (process.genJetParticles+ process.recoGenJets))

process.genParticles.abortOnUnknownPDGCode = False

##
process.load("FastSimulation.Configuration.EventContent_cff")
process.output = cms.OutputModule("PoolOutputModule",
     process.FEVTSIMEventContent,
     fileName = cms.untracked.string('pompyt_gen.root')
)
##process.outpath = cms.EndPath(process.FEVT)

process.p1 = cms.Path(process.generator)
process.p2 = cms.Path(process.caloJetMetGen)
process.out_step = cms.EndPath(process.output)

process.schedule = cms.Schedule(process.p1,process.p2,process.out_step)
