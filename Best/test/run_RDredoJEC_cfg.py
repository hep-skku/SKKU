import FWCore.ParameterSet.Config as cms
import os

MC_flag = False
#MC_flag = True

process = cms.Process("Ana")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")
#process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
#process.GlobalTag.globaltag = "START53_V27::All"
process.GlobalTag.globaltag = "FT53_V21A_AN6::All"  #Data

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
       # 'file:./patRefSel_muJets.root',    
    ),
)
for line in open('../samples/Run2012B.txt').readlines():

    line = line.strip("'\", \n")
    if '.root' not in line: continue
    process.source.fileNames.append(line)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("result_Run2012.root"),
)

process.genParticleCount = cms.EDFilter("GenParticleCountFilter",
    src = cms.InputTag("genParticles"),
    cut = cms.string("status == 3"),
    ids = cms.vint32(13, -13, 11, -11),
    minNumber = cms.uint32(1), 
    maxNumber = cms.uint32(1),
)

process.genParticleTauVeto = cms.EDFilter("GenParticleCountFilter",
    src = cms.InputTag("genParticles"),
    cut = cms.string("status == 3"),
    ids = cms.vint32(15, -15),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(0),
)

process.genParticleMuVeto = cms.EDFilter("GenParticleCountFilter",
    src = cms.InputTag("genParticles"),
    cut = cms.string("status == 3"),
    ids = cms.vint32(13,-13),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(0),
)

process.printDecay = cms.EDAnalyzer("ParticleDecayDrawer",
    src = cms.InputTag("genParticles"),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(False)
)

process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25),
)

process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter", 
    src = cms.InputTag('offlinePrimaryVertices'),
    filterParams =  cms.PSet(
        minNdof = cms.double(4.),
        maxZ    = cms.double(24.), 
        maxRho  = cms.double(2.)
    ),
    filter = cms.bool(True),
)

process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
HLTPaths = {
    "SingleMu_7E33":["HLT_IsoMu24_eta2p1_v*",],
    "Mu_53X_GTV7":["HLT_IsoMu24_eta2p1_v*",],
}
if MC_flag:
    process.hltHighLevel.HLTPaths = HLTPaths["Mu_53X_GTV7"]
else:
    process.hltHighLevel.HLTPaths = HLTPaths["SingleMu_7E33"]

#process.load("SKKU.Best.TopCleanJetSelectorMC_cfi")
process.load("SKKU.Best.TopCleanJetSelectorRD_cfi")
process.load("SKKU.Best.EventWeightProducer_cfi")
#process.load("SKKU.Best.pdfWeight_cff")

process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")

process.cleanJets.redoJES = cms.bool(True)
process.cleanJets.jecFileNames = [
    "SKKU/Best/data/Summer13_V4/Summer13_V4_DATA_txts_fromDB/Summer13_V4_DATA_L1FastJet_AK5PFchs.txt",
    "SKKU/Best/data/Summer13_V4/Summer13_V4_DATA_txts_fromDB/Summer13_V4_DATA_L2Relative_AK5PFchs.txt",
    "SKKU/Best/data/Summer13_V4/Summer13_V4_DATA_txts_fromDB/Summer13_V4_DATA_L3Absolute_AK5PFchs.txt",
    #"SKKU/Best/data/Summer13_V4/Summer13_V4_DATA_txts_fromDB/Summer13_V4_DATA_L2L3Residual_AK5PFchs.txt",
    "SKKU/Best/data/Summer13_V4/PTFIXV2_FT_53_V21_AN5_private_L2L3Residual_AK5PFchs.txt",
]
process.cleanJets.uncFilename = cms.string("SKKU/Best/data/Summer13_V5_DATA_UncertaintySources_AK5PFchs.txt")
#process.cleanJets.uncFilename = cms.string("SKKU/Best/data/Summer13_V5_MC_Uncertainty_AK5PFchs.txt")
process.cleanJets.uncSource   = cms.string("Total")

process.decaySubset.fillMode = cms.string("kME")

process.event = cms.EDAnalyzer("EventTupleProducerMuon",
    doMCMatch = cms.bool(True),
    ttGenEvent = cms.InputTag("genEvt"),
    #pdfWeights = cms.InputTag("pdfWeight"),
    gen = cms.InputTag("genParticles"),
    jet = cms.InputTag("cleanJets"),
    met = cms.InputTag("patMETsPF"),
    useEventCounter = cms.bool( True ),
    vertex = cms.InputTag("goodOfflinePrimaryVertices"),
    lepton = cms.InputTag("goodPatMuonsPF"),
    #lepton = cms.InputTag("goodElectronsPF"),
    leptonCut = cms.string(
        "abs(eta) < 2.1 && pt > 30 && abs(dB) < 0.2"
        " && isPFMuon && isGlobalMuon && normChi2 < 10"
        " && track.hitPattern.trackerLayersWithMeasurement > 5"
        " && globalTrack.hitPattern.numberOfValidMuonHits > 0"
        " && innerTrack.hitPattern.numberOfValidPixelHits > 0"
        " && numberOfMatchedStations > 1"
        " && (chargedHadronIso+max(0.,neutralHadronIso+photonIso-0.5*puChargedHadronIso)) < 0.12*pt" # relative isolation w/ Delta beta corrections (factor 0.5)
    ),
    jetCut = cms.string(
      #" pt > 20 "
      " abs(eta) < 2.5 && pt > 30 && numberOfDaughters() > 1"
      " && neutralHadronEnergyFraction() < 0.99 && neutralEmEnergyFraction() < 0.99"
      " && (abs(eta) >= 2.4 || chargedEmEnergyFraction() < 0.99)"
      " && (abs(eta) >= 2.4 || chargedHadronEnergyFraction() > 0)"
      " && (abs(eta) >= 2.4 || chargedMultiplicity() > 0)"
    ),
    bTagType = cms.string("combinedSecondaryVertexBJetTags"),
)

process.eventUp = process.event.clone(
    jet = cms.InputTag("cleanJets", "up"),
    met = cms.InputTag("cleanJets", "up"),
    )
process.eventDn = process.event.clone(
    jet = cms.InputTag("cleanJets", "dn"),
    met = cms.InputTag("cleanJets", "dn"),
    )
process.eventJERUp = process.event.clone(
    jet = cms.InputTag("cleanJets", "resUp"),
    met = cms.InputTag("cleanJets", "resUp"),
    )
process.eventJERDn = process.event.clone(
    jet = cms.InputTag("cleanJets", "resDn"),
    met = cms.InputTag("cleanJets", "resDn"),
    )

if MC_flag:
    process.hltHighLevel.HLTPaths = HLTPaths["Mu_53X_GTV7"]
    process.p = cms.Path(
        process.genParticleCount + process.genParticleTauVeto
        #+ process.genParticleMuVeto
        #  process.printDecay
        + process.goodOfflinePrimaryVertices
        + process.hltHighLevel
        * process.makeGenEvt
        #* process.pdfWeight
        * process.PUweight
        * process.cleanJets
        * process.event
        * process.eventUp
        * process.eventDn
        * process.eventJERUp
        * process.eventJERDn
        )
else:
    process.hltHighLevel.HLTPaths = HLTPaths["SingleMu_7E33"]
    process.p = cms.Path(
        process.goodOfflinePrimaryVertices
        + process.hltHighLevel
        * process.cleanJets
        * process.event
        * process.eventUp
        * process.eventDn
        )

