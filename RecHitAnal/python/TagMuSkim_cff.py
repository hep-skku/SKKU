import FWCore.ParameterSet.Config as cms

### HLT filter
import copy
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
ZMuHLTFilter = copy.deepcopy(hltHighLevel)
#ZMuHLTFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT") #default for Data
ZMuHLTFilter.throw = cms.bool(False)
ZMuHLTFilter.HLTPaths = ["HLT_Mu*","HLT_IsoMu*","HLT_DoubleMu*"]

### Z -> MuMu candidates

# Get muons of needed quality for Zs
looseMuonsForZ = cms.EDFilter("MuonSelector",
                             src = cms.InputTag("muons"),
                             cut = cms.string('pt > 10 && abs(eta)<2.4 && isGlobalMuon = 1 && isTrackerMuon = 1 && abs(innerTrack().dxy)<2.0'),
                             filter = cms.bool(True)                                
                             )

tightMuonsForZ = cms.EDFilter("MuonSelector",
                             src = cms.InputTag("looseMuonsForZ"),
                             cut = cms.string('pt > 15 && globalTrack().normalizedChi2<10.0 && isolationR03().sumPt<3.0 && (isolationR03().emEt+isolationR03().hadEt+isolationR03().sumPt)<0.2*pt && track().hitPattern().numberOfValidTrackerHits()>10'),
                             filter = cms.bool(True)                                
                             )

generalTracksFilter = cms.EDFilter("TrackCountFilter",
                                   src = cms.InputTag('generalTracks'),
                                   cut = cms.string('pt > 5 && abs(eta)<2.4'),
                                   minNumber = cms.uint32(2)
                                   )

# TagMu Skim sequence
TagMuSelSeq = cms.Sequence(ZMuHLTFilter *
                           looseMuonsForZ *
                           tightMuonsForZ *
                           generalTracksFilter
                           )

