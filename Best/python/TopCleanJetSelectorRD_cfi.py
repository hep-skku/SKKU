import FWCore.ParameterSet.Config as cms

cleanJets = cms.EDFilter("TopJetSelector",
    doFilter = cms.bool(False),
    debug = cms.untracked.bool(False),

    isMC = cms.bool(False),

    jet = cms.InputTag("veryLoosePatJetsPF"),
    met = cms.InputTag("patMETsPF"),

    redoJES = cms.bool(True),
    rho = cms.InputTag("ak5PFJets", "rho"),
    vertex = cms.InputTag("goodOfflinePrimaryVertices"),
    jecFileNames = cms.vstring(
        "SKKU/Best/data/Summer13_V4/Summer13_V4_Data_txts_fromDB/Summer13_V4_Data_L1FastJet_AK5PFchs.txt",
        "SKKU/Best/data/Summer13_V4/Summer13_V4_Data_txts_fromDB/Summer13_V4_Data_L2Relative_AK5PFchs.txt",
        "SKKU/Best/data/Summer13_V4/Summer13_V4_Data_txts_fromDB/Summer13_V4_Data_L3Absolute_AK5PFchs.txt",
        "SKKU/Best/data/Summer13_V4/Summer13_V4_Data_txts_fromDB/Summer13_V4_Data_L2L3Residual_AK5PFchs.txt",
    ),

    selection = cms.PSet(
        cut = cms.string(""),
        minPt = cms.double(30),
        maxEta = cms.double(2.5),
    ),
    cleaning = cms.PSet(
        overlapCands = cms.VInputTag(),
        overlapDeltaR = cms.double(0.5),
        #cleanMethod = cms.untracked.string("subtract"),
        cleanMethod = cms.string(""),
        #cleanMethod = cms.untracked.string(""),
    ),
    #jecFileRD = cms.string("KrAFT/RecoSelectorTools/data/JEC/Summer13_V4/Summer13_V4_DATA_Uncertainty_AK5PFchs.txt"),
    #jecFileMC = cms.string("KrAFT/RecoSelectorTools/data/JEC/Summer13_V4/Summer13_V4_MC_Uncertainty_AK5PFchs.txt"),
    #jecSource = cms.string("Total"),
    uncFilename = cms.string("SKKU/Best/data/Summer13_V5_MC_Uncertainty_AK5PFchs.txt"),
    uncSource = cms.string(""),

    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999),
)

