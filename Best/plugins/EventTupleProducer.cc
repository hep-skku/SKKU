#ifndef BESTAnalysis_TTbarLeptonJet_EventTupleProducer_H
#define BESTAnalysis_TTbarLeptonJet_EventTupleProducer_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
//#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
//#include "DataFormats/Common/interface/MergeableCounter.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"

#include "CommonTools/Utils/interface/PtComparator.h"

#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TTree.h"
#include "TH1F.h"

#include <memory>
#include <vector>
#include <string>

template<typename Lepton>
class EventTupleProducer : public edm::EDAnalyzer
{
public:
  EventTupleProducer(const edm::ParameterSet& pset);
  ~EventTupleProducer();

  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

private:
  bool hasMother(const reco::Candidate* p, const int pdgId);

private:
  bool doMCMatch_;

  // Input objects
  edm::InputTag genLabel_;
  edm::InputTag leptonLabel_;
  edm::InputTag jetLabel_;
  edm::InputTag metLabel_;
  edm::InputTag vertexLabel_;
  std::string bTagType_;
  edm::InputTag ttGenEvent_;

  // Cuts
  StringCutObjectSelector<Lepton, true>* isGoodLepton_;
  StringCutObjectSelector<pat::Jet, true>* isGoodJet_;

  // Output tree
  TTree* tree_;
  int run_, lumi_, event_;
  math::XYZTLorentzVector lepton_, met_;
  int charge_;
  double eventWeight_, eventWeightUp_, eventWeightDn_;
  double ptWeight_, ptWeightUp_, ptWeightDn_;
  int nPV_;
  std::vector<math::XYZTLorentzVector>* jets_;
  std::vector<double>* bTags_;
  std::vector<int>* jetsFlavor_;
  std::vector<int>* jetMCBits_;
//  std::vector<double>* pTGen_; 
  std::vector<int>* pdgGen_;
//  std::vector<double>* pdfWeight_;
};

template<typename Lepton>
EventTupleProducer<Lepton>::EventTupleProducer(const edm::ParameterSet& pset)
{
  doMCMatch_ = pset.getParameter<bool>("doMCMatch");

  // Input labels
  genLabel_ = pset.getParameter<edm::InputTag>("gen");
  jetLabel_ = pset.getParameter<edm::InputTag>("jet");
  metLabel_ = pset.getParameter<edm::InputTag>("met");
  leptonLabel_ = pset.getParameter<edm::InputTag>("lepton");
  vertexLabel_ = pset.getParameter<edm::InputTag>("vertex");
  ttGenEvent_ = pset.getParameter<edm::InputTag>("ttGenEvent");

  std::string leptonCut = pset.getParameter<std::string>("leptonCut");
  isGoodLepton_ = new StringCutObjectSelector<Lepton, true>(leptonCut);
  std::string jetCut = pset.getParameter<std::string>("jetCut");
  isGoodJet_ = new StringCutObjectSelector<pat::Jet, true>(jetCut);

  bTagType_ = pset.getParameter<std::string>("bTagType");

  // Output histograms and tree
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "Mixed event tree");
  tree_->Branch("run", &run_, "run/I");
  tree_->Branch("lumi", &lumi_, "lumi/I");
  tree_->Branch("event", &event_, "event/I");

  tree_->Branch("lepton", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>", &lepton_);
  tree_->Branch("charge", &charge_, "charge/I");
  tree_->Branch("met", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>", &met_);
  tree_->Branch("vertex", &nPV_, "nPV_/I");
  
  jets_ = new std::vector<math::XYZTLorentzVector>();
  bTags_ = new std::vector<double>();
  jetsFlavor_ = new std::vector<int>();
  jetMCBits_ = new std::vector<int>();
//  pTGen_ = new std::vector<double>();
  pdgGen_ = new std::vector<int>();
//  pdfWeight_ = new std::vector<double>();
  //eventWeight_ = new std::vector<double>();
  tree_->Branch("jets", "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &jets_);
  tree_->Branch("bTags", "std::vector<double>", &bTags_);
  tree_->Branch("jetsFlavor", "std::vector<int>", &jetsFlavor_);
  tree_->Branch("jetMCBits", "std::vector<int>", &jetMCBits_);
//  tree_->Branch("pTGen", "std::vector<double>", &pTGen_);
  tree_->Branch("pdgGen", "std::vector<int>", &pdgGen_);
//  tree_->Branch("pdfWeight", "std::vector<double>", &pdfWeight_);
  tree_->Branch("PUweight", &eventWeight_, "PUweight/d");
  tree_->Branch("PUweightup", &eventWeightUp_, "PUweightup/d");
  tree_->Branch("PUweightdn", &eventWeightDn_, "PUweightdn/d");
  tree_->Branch("Ptweight", &ptWeight_, "Ptweight/d");
  tree_->Branch("Ptweightup", &ptWeightUp_, "Ptweightup/d");
  tree_->Branch("Ptweightdn", &ptWeightDn_, "Ptweightdn/d");
}

template<typename Lepton>
EventTupleProducer<Lepton>::~EventTupleProducer()
{
  if ( jets_ ) delete jets_;
  if ( bTags_ ) delete bTags_;
  if ( jetsFlavor_ ) delete jetsFlavor_;
  if ( jetMCBits_ ) delete jetMCBits_;
//  if ( pTGen_ ) delete pTGen_;
//  if ( pdfWeight_ ) delete pdfWeight_;
  if ( pdgGen_ ) delete pdgGen_;
}

template<typename Lepton>
void EventTupleProducer<Lepton>::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  using namespace std;

  jets_->clear();
  bTags_->clear();
  jetsFlavor_->clear();
  jetMCBits_->clear();
//  pTGen_->clear();
//  pdfWeight_->clear();
  pdgGen_->clear();

  edm::Handle<edm::View<reco::Vertex> > vertexHandle;
  event.getByLabel(vertexLabel_, vertexHandle);
  nPV_ = vertexHandle->size();

//  eventWeight_ = eventWeightUp_ = eventWeightDn_ = 1;
 
  edm::Handle<double> eventWeightHandle;
  event.getByLabel(edm::InputTag("PUweight", "weight"), eventWeightHandle);
  eventWeight_= *eventWeightHandle;

  edm::Handle<double> eventWeightUpHandle;
  event.getByLabel(edm::InputTag("PUweight", "weightplus"), eventWeightUpHandle);
  eventWeightUp_= *eventWeightUpHandle;

  edm::Handle<double> eventWeightDnHandle;
  event.getByLabel(edm::InputTag("PUweight", "weightminus"), eventWeightDnHandle);
  eventWeightDn_= *eventWeightDnHandle;
/*
  edm::Handle<std::vector<double> > pdfWeightHandle;
  event.getByLabel(edm::InputTag("pdfWeight"), pdfWeightHandle);
  for ( int i=0, n=pdfWeightHandle->size(); i<n; i++ )
  {
    pdfWeight_->push_back(pdfWeightHandle->at(i));
  }
*/
  edm::Handle<edm::View<Lepton> > leptonHandle;
  event.getByLabel(leptonLabel_, leptonHandle);
  int nPassingLepton = 0;
  int leptonIdx = -1;
  for ( int i=0, n=leptonHandle->size(); i<n; ++i )
  {
    if ( !(*isGoodLepton_)(leptonHandle->at(i)) ) continue;

    ++nPassingLepton;
    leptonIdx = i;
  }
  if ( nPassingLepton != 1 ) return;
  lepton_ = leptonHandle->at(leptonIdx).p4();
  charge_ = leptonHandle->at(leptonIdx).charge();

  //--https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
  edm::Handle<TtGenEvent> genEvt;
  event.getByLabel(ttGenEvent_, genEvt);

  //--https://cmssdt.cern.ch/SDT/lxr/source/AnalysisDataFormats/TopObjects/src/TtGenEvent.cc
  //--TtGenEvent::leptonicDecayTop(bool excludeTauLeptons) const //default: excludeTauLeptons=false

  //--pT reweight lepton+jets and dilepton channel 
  double topPt = 0, topBarPt = 0;
  double SF_top = 1, SF_topBar = 1;
  ptWeight_ = ptWeightUp_ = ptWeightDn_ = 1;

  //bool fromWBoson = true;
  int nLeptons = genEvt->numberOfLeptons();
  if ( nLeptons > 0 ) {
    topPt    = genEvt->top()->pt();
    topBarPt = genEvt->topBar()->pt();
    if ( nLeptons == 1 ) {
      //topPt    = genEvt->leptonicDecayTop()->pt(); topBarPt = genEvt->hadronicDecayTop()->pt();
      SF_top = TMath::Exp(0.159+((-0.00141)*topPt));
      SF_topBar = TMath::Exp(0.159+((-0.00141)*topBarPt));
    } else if ( nLeptons == 2 ) {
      //genEvt->leptonicDecayTop()->pdgId()>0 ? topPt = genEvt->leptonicDecayTop()->pt() : topBarPt = genEvt->leptonicDecayTop()->pt();
      SF_top = TMath::Exp(0.148+((-0.00129)*topPt));
      SF_topBar = TMath::Exp(0.148+((-0.00129)*topBarPt));
    }
    ptWeight_ = sqrt(SF_top*SF_topBar);
    ptWeightUp_ = ptWeight_*ptWeight_; 
    ptWeightDn_ = 1;
  }

  edm::Handle<edm::View<pat::MET> > metHandle;
  event.getByLabel(metLabel_, metHandle);
  met_ = metHandle->at(0).p4();

  const reco::Candidate* genLepB, * genHadB, * genHadJ1, * genHadJ2;
  genLepB = genHadB = genHadJ1 = genHadJ2 = 0;
  if ( doMCMatch_ )
  {
    edm::Handle<reco::GenParticleCollection> genHandle;
    event.getByLabel(genLabel_, genHandle);
    if ( !genHandle.isValid() ) 
    {
      doMCMatch_ = false;
    }
    else
    {
      // Find top quark from the genParticles
      const reco::GenParticle* t1 = 0, * t2 = 0;
      const reco::GenParticle* w1 = 0, * w2 = 0;
      const reco::GenParticle* b1 = 0, * b2 = 0;
      const reco::GenParticle* l1 = 0, * l2 = 0;

      std::vector<const reco::GenParticle*> w1Decay, w2Decay;
      for ( int i=0, n=genHandle->size(); i<n; ++i )
      {
        const reco::GenParticle& p = genHandle->at(i);
        if ( p.status() != 3 ) continue;

        switch(p.pdgId())
        {
          case   6: t1 = &p; break;
          case  -6: t2 = &p; break;
          case  24: w1 = &p; break;
          case -24: w2 = &p; break;
          case   5: b1 = &p; break;
          case  -5: b2 = &p; break;
          case  15: break;
          case -15: break;
          case -11: case -13: l1 = &p; break;
          case  11: case  13: l2 = &p; break;
          //case -11: case -13: case -15: l1 = &p; break;
          //case  11: case  13: case  15: l2 = &p; break;
          case  12: case  14: case  16:
          case -12: case -14: case -16: break;
          default:
            if ( hasMother(&p, 24) ) w1Decay.push_back(&p);
            else if ( hasMother(&p, -24) ) w2Decay.push_back(&p);
        }
      //pTGen_->push_back(p.pt());
      pdgGen_->push_back(p.pdgId());

      }

      if ( !t1 or !t2 or !w1 or !w2 or !b1 or !b2 ) return;
      if ( (l1 and l2) or (!l1 and !l2) ) return;

      if ( l1 and w2Decay.size() > 1 )
      {
        genLepB = b1;
        genHadB = b2;

        genHadJ1 = w2Decay[0];
        genHadJ2 = w2Decay[1];
      }
      else if ( l2 and w1Decay.size() > 1 )
      {
        genLepB = b2;
        genHadB = b1;

        genHadJ1 = w1Decay[0];
        genHadJ2 = w1Decay[1];
      }
      else
      {
        cout << "FATAL: This should not happen\n";
        return;
      }
    }
  }

  edm::Handle<edm::View<pat::Jet> > jetHandle;
  event.getByLabel(jetLabel_, jetHandle);
  for ( int i=0, n=jetHandle->size(); i<n; ++i )
  {
    const pat::Jet& jet = jetHandle->at(i);

    if ( !(*isGoodJet_)(jet) ) continue;
    //cout<<jet.partonFlavour() <<endl;
    jets_->push_back(jet.p4());
    bTags_->push_back(jet.bDiscriminator(bTagType_));
    jetsFlavor_->push_back(jet.partonFlavour());

    int jetMCBit = 0;
    const double dRLepB = genLepB ? deltaR(jet, *genLepB) : 1e9;
    const double dRHadB = genHadB ? deltaR(jet, *genHadB) : 1e9;
    const double dRHadJ1 = genHadJ1 ? deltaR(jet, *genHadJ1) : 1e9;
    const double dRHadJ2 = genHadJ2 ? deltaR(jet, *genHadJ2) : 1e9;

    if ( dRLepB < 0.5 ) jetMCBit |= 1;
    if ( dRHadB < 0.5 ) jetMCBit |= 2;
    if ( dRHadJ1 < 0.5 or dRHadJ2 < 0.5 ) jetMCBit |= 4;

    jetMCBits_->push_back(jetMCBit);
  }
  if ( jets_->size() < 3 ) return;

  std::sort(jets_->begin(), jets_->end(), GreaterByPt<math::XYZTLorentzVector>());

  // Now put jets in current event to the event cache
  run_ = event.run();
  lumi_ = event.luminosityBlock();
  event_ = event.id().event();

  tree_->Fill();

}

template<typename Lepton>
bool EventTupleProducer<Lepton>::hasMother(const reco::Candidate* p, const int pdgId)
{
  for ( int i=0, n=p->numberOfMothers(); i<n; ++i )
  {
    const reco::Candidate* m = p->mother(i);
    if ( !m ) return false;
    else if ( m->pdgId() == pdgId ) return true;
    else if ( hasMother(m, pdgId) ) return true;
  }

  return false;
}

typedef EventTupleProducer<pat::Electron> EventTupleProducerElectron;
typedef EventTupleProducer<pat::Muon> EventTupleProducerMuon;

DEFINE_FWK_MODULE(EventTupleProducerElectron);
DEFINE_FWK_MODULE(EventTupleProducerMuon);

#endif
