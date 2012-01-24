import FWCore.ParameterSet.Config as cms
process = cms.Process("ZMuSkim")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "START42_V12::All"
process.prefer("GlobalTag")

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #'/../user/h/hkseo/RPC/ZMuSkim/zmumuJet_384_10_1_Ibq.root'
    #'file:/tmp/hkseo/FE67038D-A57C-E011-9E83-0030487F16BF.root'
    'file:/korserv001/mskim/rpc/PromptReco-v4_163592_0A68E6C1-377B-E011-8BC8-00237DA28240.root'
    )
)

# WZMu Skim
process.load("DPGAnalysis.Skims.ZMuSkim_cff")
process.Zskim = cms.Sequence(
    process.ZMuHLTFilter *
    process.looseMuonsForZ *
    process.tightMuonsForZ *
    process.dimuons *
    process.dimuonsFilter
    )

# Good vertex requirement
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24),
                                           maxd0 = cms.double(2)
                                           )


process.p1 = cms.Path(process.primaryVertexFilter*process.Zskim)
#process.p1 = cms.Path(process.diMuonSelSeq*process.demo)
#process.p2 = cms.Path(process.pfMetWMuNuSeq*process.demo)

# Output module configuration
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('ZMuSkim.root'),
                               # save only events passing the full path
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p1') )
                               )
process.outpath = cms.EndPath(process.out)