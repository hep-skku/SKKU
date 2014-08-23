#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "Math/GenVector/LorentzVector.h"
#include "Math/VectorUtil_Cint.h"

using namespace std;
using namespace ROOT::Math::VectorUtil;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

const int cut_minNJet = 4;
//const double cut_minJetPt = 35;
//const double cut_maxJetEta = 2.5;
const double cut_maxJetEta = 2.4; //since Dec 24, 2013 (twiki for SFb)
const double cut_minLeadJetPt = 45; //considered to be disabled from Dec 30, 2013

//const int cut_minNBjet = 2;
const double cut_bTag = 0.679; //medium

const double cut_drSEvt = 0.4; //for SameEvent
const double cut_drBEvt = 0.5; //for BiEvent (used to be 0.4)

struct EventTopology
{
  int run, lumi, event;
//  double PUweight; //PUweightup, PUweightdn
  LorentzVector* lepton;
  LorentzVector* met;
  std::vector<LorentzVector>* jets;
  int charge;
  std::vector<double>* bTags;
  std::vector<int>* jetMCBits;
//  std::vector<double>* pTGen;

  EventTopology()
  {
    lepton = new LorentzVector;
    met = new LorentzVector;
    jets = new std::vector<LorentzVector>();
    bTags = new std::vector<double>();
    jetMCBits = new std::vector<int>();
//    pTGen = new std::vector<double>();
  };
};

void setBranch(TTree* tree, EventTopology& e)
{
  tree->SetBranchAddress("run"  , &e.run  );
  tree->SetBranchAddress("lumi" , &e.lumi );
  tree->SetBranchAddress("event", &e.event);
//  tree->SetBranchAddress("PUweight", &e.PUweight);

  tree->SetBranchAddress("lepton", &e.lepton);
  tree->SetBranchAddress("met", &e.met);
  tree->SetBranchAddress("jets", &e.jets);

  tree->SetBranchAddress("charge", &e.charge);
  tree->SetBranchAddress("bTags", &e.bTags);
  tree->SetBranchAddress("jetMCBits", &e.jetMCBits);
//  tree->SetBranchAddress("pTGen", &e.pTGen);
}

void mixData(TString input = "Run2012", int cutPt = 40, int cutNb = 2)
{
  cout << endl << " ************************************************* " << endl
    //<< endl << " WARNING! PUweight = 1, so need to update PUweight " << endl
       << endl << " maxJetEta = " <<cut_maxJetEta<<" (used to be 2.5) " << endl
       << endl << " ************************************************* " << endl << endl;
  const char* sample = input;
  const int cut_minJetPt = cutPt;
  const int cut_minNBjet = cutNb;

  TFile* file1 = TFile::Open(Form("../2_runNtuple/ntuples/result_%s.root", sample));
  TFile* file2 = TFile::Open(Form("../2_runNtuple/ntuples/result_%s.root", sample));
  TTree* tree1 = (TTree*)file1->Get("event/tree");//eventUp, eventDn, eventJERUp, eventJERDn
  TTree* tree2 = (TTree*)file2->Get("event/tree");

  EventTopology event1, event2;

  setBranch(tree1, event1);
  setBranch(tree2, event2);

  TFile* outFile = TFile::Open(Form("hist/hist_%s_pt%d_nj%d_nb%d.root", sample, cut_minJetPt, cut_minNJet, cut_minNBjet), "RECREATE");
  TDirectory* dirSEvt = outFile->mkdir("SEvt");
  TDirectory* dirBEvt = outFile->mkdir("BEvt");

  dirSEvt->cd();
  const double binShift = 2.5; //0;
  const double mJJMax = 500+binShift;
  const double mJJBMax = 1000+ 2*binShift;
  const int nbin = 100;
  const TString axisTitleMw = ";Dijet mass (GeV/c^{2});Events per 5GeV/c^{2}";
  const TString axisTitleMt = ";Three jet mass (GeV/c^{2});Events per 10GeV/c^{2}";

  TH1F* hSEvt_Mw = new TH1F("hMw", "Dijet mass"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hSEvt_Mw_JJ = new TH1F("hMw_JJ", "J+J"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hSEvt_Mw_JK = new TH1F("hMw_JK", "J+K"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hSEvt_Mw_KK = new TH1F("hMw_KK", "K+K"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hSEvt_Mw_HBJ = new TH1F("hMw_HBJ", "HadB+J"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hSEvt_Mw_LBJ = new TH1F("hMw_LBJ", "LepB+J"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hSEvt_Mw_BX = new TH1F("hMw_BX", "B+X"+axisTitleMw, nbin, binShift, mJJMax);

  TH1F* hSEvt_Mt = new TH1F("hMt", "Three jet mass"+axisTitleMt, nbin, 2*binShift, mJJBMax);
  TH1F* hSEvt_Mt_JJHB = new TH1F("hMt_JJHB", "J+J+HadB"+axisTitleMt, nbin, 2*binShift, mJJBMax);
  TH1F* hSEvt_Mt_JJLB = new TH1F("hMt_JJLB", "J+J+LepB"+axisTitleMt, nbin, 2*binShift, mJJBMax);
  //TH1F* hSEvt_Mt_KXY = new TH1F("hMt_KXY", "K+XY"+axisTitleMt, nbin, 2*binShift, mJJBMax);
  //TH1F* hSEvt_Mt_KKX = new TH1F("hMt_KKX", "KK+X"+axisTitleMt, nbin, 2*binShift, mJJBMax);
  TH1F* hSEvt_Mt_XYZ = new TH1F("hMt_XYZ", "Others"+axisTitleMt, nbin, 2*binShift, mJJBMax);
  TH1F* hSEvt_Mt_HBLBJ = new TH1F("hMt_HBLBJ", "J+HadB+LepB"+axisTitleMt, nbin, 2*binShift, mJJBMax);
  TH1F* hSEvt_Mt_SFb = new TH1F("hMt_SFb", "SFb", 100, 0.9, 1);
  TH1F* hSEvt_nJets = new TH1F("hnjets","hnjets",10,0,10);

  dirBEvt->cd();
  TH1F* hBEvt_Mw = new TH1F("hMw", "Dijet mass"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hBEvt_Mw_JJ = new TH1F("hMw_JJ", "J+J"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hBEvt_Mw_JK = new TH1F("hMw_JK", "J+K"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hBEvt_Mw_KK = new TH1F("hMw_KK", "K+K"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hBEvt_Mw_HBJ = new TH1F("hMw_HBJ", "HadB+J"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hBEvt_Mw_LBJ = new TH1F("hMw_LBJ", "LepB+J"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hBEvt_Mw_BX = new TH1F("hMw_BX", "B+X"+axisTitleMw, nbin, binShift, mJJMax);

  TH1F* hBEvt_Mt = new TH1F("hMt", "Three jet mass"+axisTitleMt, nbin, 2*binShift, mJJBMax);
  TH1F* hBEvt_Mt_JJHB = new TH1F("hMt_JJHB", "J+J+HadB"+axisTitleMt, nbin, 2*binShift, mJJBMax);
  TH1F* hBEvt_Mt_JJLB = new TH1F("hMt_JJLB", "J+J+LepB"+axisTitleMt, nbin, 2*binShift, mJJBMax);
  //TH1F* hBEvt_Mt_KXY = new TH1F("hMt_KXY", "K+XY"+axisTitleMt, nbin, 2*binShift, mJJBMax);
  //TH1F* hBEvt_Mt_KKX = new TH1F("hMt_KKX", "KK+X"+axisTitleMt, nbin, 2*binShift, mJJBMax);
  TH1F* hBEvt_Mt_XYZ = new TH1F("hMt_XYZ", "Others"+axisTitleMt, nbin, 2*binShift, mJJBMax);
  TH1F* hBEvt_Mt_HBLBJ = new TH1F("hMt_HBLBJ", "J+HadB+LepB"+axisTitleMt, nbin, 2*binShift, mJJBMax);
  TH1F* hBEvt_Mt_SFb = new TH1F("hMt_SFb", "SFb", 100, 0.9, 1);
  TH1F* hBEvt_nJets = new TH1F("hnjets","hnjets",10,0,10);

  //--https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#2012_Data_and_MC_EPS13_prescript (with ttbar measurements)
  float ptmin[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600};
  float ptmax[] = {30, 40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 600, 800};

  double SFb_error[] = {
  0.0415707,
  0.0204209,
  0.0223227,
  0.0206655,
  0.0199325,
  0.0174121,
  0.0202332,
  0.0182446,
  0.0159777,
  0.0218531,
  0.0204688,
  0.0265191,
  0.0313175,
  0.0415417,
  0.0740446,
  0.0596716 };

  double SFb = 0;
  double SFb_up = 0;
  double SFb_dn = 0;

  double SFb2 = 0;
  double SFb_up2 = 0;
  double SFb_dn2 = 0;

  float etamin[] = {0.0, 0.5, 1.1, 1.7, 2.3};
  float etamax[] = {0.5, 1.1, 1.7, 2.3, 5.0};
  
  double JER[] = {1.052, 1.057, 1.096, 1.134, 1.288};
  double JERup[] = {0.990, 1.001, 1.032, 1.042, 1.089};
  double JERdn[] = {1.115, 1.114, 1.161, 1.228, 1.488};
 
  double SFJER = 0;
  double SFJER_up = 0;
  double SFJER_dn = 0;

  int nLeadJetsInSameEvt = 0, nLeadJetsInBiEvt = 0;
  int nBiEvent = 0, nPassedBiEvent = 0; // Variables to restore overlap removal scale
  for ( int iEvent1=0, nEvent=tree1->GetEntries(); iEvent1<nEvent; ++iEvent1 )
  {
    tree1->GetEntry(iEvent1);
    const LorentzVector& lepton1 = *event1.lepton;

    //if ( event1.jets->size() < cut_minNJet ) continue;
    //if ( event1.jets->at(0).pt() < cut_minLeadJetPt or event1.jets->at(1).pt() < cut_minLeadJetPt ) continue;
    int nJets1 = 0, nBjets1 = 0;
    for ( int j1=0, nj=event1.jets->size(); j1<nj; ++j1 )
    {
      const LorentzVector jet1 = event1.jets->at(j1);
      const double bTag1 = event1.bTags->at(j1);
      //if ( event1.jets->at(0).pt() < cut_minLeadJetPt or event1.jets->at(1).pt() < cut_minLeadJetPt ) continue;
      if ( jet1.pt() > cut_minJetPt and abs(jet1.eta()) < cut_maxJetEta ) ++nJets1;
      if ( bTag1 > cut_bTag ) ++nBjets1;
    }
    if ( nJets1 < cut_minNJet || nBjets1 < cut_minNBjet ) continue;
    //cout << "iEvent1 = " << iEvent1 << ": " << nJets1 << " " << nBjets1 << endl;
    //int nBjet1 = 0;
    //for ( int j1=0, nj=event1.bTags->size(); j1<nj; ++j1 )
    //{
    //  if ( event1.bTags->at(j1) > cut_bTag ) ++nBjet1;
    //}
    //if ( nBjet1 < cut_minNBjet ) continue;

    //int njet1 = 0, njet2 = 0, nbjet1 = 0, nbjet2 = 0;
    // Make SEvt event jet combinations
    for ( int j1=0, nj=event1.jets->size(); j1<nj; ++j1 )
    {
      const LorentzVector jet1 = event1.jets->at(j1);
      const double bTag1 = event1.bTags->at(j1);
      const int mcBit1 = event1.jetMCBits->at(j1);
      //if ( DeltaR(jet1, lepton1) < 0.3 ) continue;

      if ( jet1.pt() < cut_minJetPt or abs(jet1.eta()) > cut_maxJetEta ) continue;
      //if ( jet1.pt() < cut_minLeadJetPt or abs(jet1.eta()) > cut_maxJetEta ) continue;
      //if ( mcBit1 == 1 or mcBit1 == 2 ) continue;
      //if ( mcBit1&3 ) continue;
      if ( bTag1 > cut_bTag ) continue;
      if (jet1.pt() > cut_minLeadJetPt) ++nLeadJetsInSameEvt;
      
      for ( int j2=j1+1; j2<nj; ++j2 )
      {
        const LorentzVector jet2 = event1.jets->at(j2);
        const double bTag2 = event1.bTags->at(j2);
        const int mcBit2 = event1.jetMCBits->at(j2);

        //if ( DeltaR(jet2, lepton1) < 0.3 ) continue;
        if ( jet2.pt() < cut_minJetPt or abs(jet2.eta()) > cut_maxJetEta ) continue;
        //if ( jet2.pt() < cut_minLeadJetPt or abs(jet2.eta()) > cut_maxJetEta ) continue;
        //if ( mcBit2 == 1 or mcBit2 == 2 ) continue;
        if ( bTag2 > cut_bTag ) continue;
        //if ( mcBit2&3 ) continue;

        if ( DeltaR(jet1, jet2) < cut_drSEvt ) continue;

        LorentzVector jj = event1.jets->at(j1)+event1.jets->at(j2);
        const double mJJ = jj.mass(); //min(mJJMax-1e-3, jj.mass());

        hSEvt_Mw->Fill(mJJ);
  if(input.Contains("Run2012") && jet1.pt() > cut_minLeadJetPt && jet2.pt() > cut_minLeadJetPt) hSEvt_Mw_JJ->Fill(mJJ);

        if ( mcBit1&4 and mcBit2&4 ) hSEvt_Mw_JJ->Fill(mJJ);
        else if ( (mcBit1&2 and mcBit2&4) or (mcBit1&4 and mcBit2&2) ) hSEvt_Mw_HBJ->Fill(mJJ);
        else if ( (mcBit1&1 and mcBit2&4) or (mcBit1&4 and mcBit2&1) ) hSEvt_Mw_LBJ->Fill(mJJ);
        else if ( mcBit1&3 or mcBit2&3 ) hSEvt_Mw_BX->Fill(mJJ);
        else if ( mcBit1 == 0 and mcBit2 == 0 ) hSEvt_Mw_KK->Fill(mJJ);
        else hSEvt_Mw_JK->Fill(mJJ);

        for ( int j3=0; j3<nj; ++j3 )
        {
          const LorentzVector jet3 = event1.jets->at(j3);
          const double bTag3 = event1.bTags->at(j3);
          const int mcBit3 = event1.jetMCBits->at(j3);

          //if ( DeltaR(jet3, lepton1) < 0.3 ) continue;
 
          if ( jet3.pt() < cut_minJetPt or abs(jet3.eta()) > cut_maxJetEta ) continue;
          if ( bTag3 <= cut_bTag ) continue;
          if ( DeltaR(jet1, jet3) < cut_drSEvt ) continue;
          if ( DeltaR(jet2, jet3) < cut_drSEvt ) continue;
          
    for ( int k=0; k<16; ++k ) {
      if ( jet3.pt() > ptmin[k] and jet3.pt() < ptmax[k] ) {
        SFb = 0.938887+(0.00017124* (jet3.pt()))+(-2.76366e-07*(jet3.pt())*(jet3.pt()));
        SFb_up = SFb + SFb_error[k];
        SFb_dn = SFb - SFb_error[k];
        //cout << "iEvent = " << iEvent1 << ": b-jet pT = " << jet3.pt() << " (k = " << k << ")" << endl;
        break;
      }
    }
    if( jet3.pt() > ptmax[15] ) {
      SFb = 0.938887+(0.00017124*ptmax[15])+(-2.76366e-07*(ptmax[15]*ptmax[15]));
      SFb_up = SFb + 2*SFb_error[15];
      SFb_dn = SFb - 2*SFb_error[15];
    }
          LorentzVector jjb = jj + event1.jets->at(j3);
          const double mJJB = jjb.mass();
          //double PUweight1 = event1.PUweight; //PUweight1 = 1;
    double SFbPUw1 = SFb;//*PUweight1; if(input.Contains("Run2012")) SFb = SFbPUw1 = 1;
          hSEvt_Mt->Fill(mJJB,SFbPUw1);
    if(input.Contains("Run2012") && jet1.pt() > cut_minLeadJetPt && jet2.pt() > cut_minLeadJetPt) hSEvt_Mt_JJHB->Fill(mJJB,SFbPUw1);

          if ( (mcBit1&2 and mcBit2&4 and mcBit3&4) or 
               (mcBit1&4 and mcBit2&2 and mcBit3&4) or 
               (mcBit1&4 and mcBit2&4 and mcBit3&2) ) hSEvt_Mt_JJHB->Fill(mJJB,SFbPUw1);
          else if ( (mcBit1&1 and mcBit2&4 and mcBit3&4) or
                    (mcBit1&4 and mcBit2&1 and mcBit3&4) or 
                    (mcBit1&4 and mcBit2&4 and mcBit3&1) ) hSEvt_Mt_JJLB->Fill(mJJB,SFbPUw1);
          else if ( (mcBit1&1 and mcBit2&2 and mcBit3&4) or
                    (mcBit1&1 and mcBit2&4 and mcBit3&2) or
                    (mcBit1&2 and mcBit2&1 and mcBit3&4) or
                    (mcBit1&2 and mcBit2&4 and mcBit3&1) or
                    (mcBit1&4 and mcBit2&1 and mcBit3&2) or
                    (mcBit1&4 and mcBit2&2 and mcBit3&1) ) hSEvt_Mt_HBLBJ->Fill(mJJB,SFbPUw1);
          //else if ( mcBit1+mcBit2 == 0 or mcBit2+mcBit3 == 0 or mcBit3+mcBit1 == 0 ) hSEvt_Mt_KKX->Fill(mJJB);
          //else if ( mcBit1 ==0 or mcBit2 == 0 or mcBit3 == 0 ) hSEvt_Mt_KXY->Fill(mJJB);
          else hSEvt_Mt_XYZ->Fill(mJJB,SFbPUw1);
          //if(abs(jet3.eta())>2.4) cout << "iEvent = " << iEvent1 << ": eta & pT = " << jet3.eta() << " " << jet3.pt() << " mJJB = " << mJJB << " (SFb = " << SFb << ": up = " << SFb_up << " dn = " << SFb_dn << ")" << endl;
    hSEvt_Mt_SFb->Fill(SFb);
        }
      }
      //break;

    }
    //cout << "iEvent = " << iEvent1 << ": nLeadJets In SameEvt = " << nLeadJetsInSameEvt << endl;
    hSEvt_nJets->Fill(nLeadJetsInSameEvt); nLeadJetsInSameEvt = 0;

    // Make Bi event jet combinations
    for ( int iEvent2 = iEvent1+1; ;++iEvent2 )
    {
      if ( iEvent2 == nEvent ) iEvent2 = 0;
      tree2->GetEntry(iEvent2);

      const LorentzVector& lepton2 = *event2.lepton;

      // ** DONOT use "continue" HERE because it causes event-mixing happen again (i.e. artifitial effect)
      //if ( event2.jets->size() < cut_minNJet ) continue;
      //if ( event2.jets->at(0).pt() < cut_minLeadJetPt or event2.jets->at(1).pt() < cut_minLeadJetPt ) continue;
      //cout << "-> " << iEvent2 << ": " << event2.jets->at(0).pt() << " " << event2.jets->at(1).pt() << endl;
      int nJets2 = 0, nBjets2 = 0;
      for ( int j2=0, nj=event2.jets->size(); j2<nj; ++j2 )
      {
  const LorentzVector jet2 = event2.jets->at(j2);
  const double bTag2 = event2.bTags->at(j2);
  //if ( event2.jets->at(0).pt() < cut_minLeadJetPt or event2.jets->at(1).pt() < cut_minLeadJetPt ) continue;
  if ( jet2.pt() > cut_minJetPt and abs(jet2.eta()) < cut_maxJetEta ) ++nJets2;
  if ( bTag2 > cut_bTag ) ++nBjets2;
      }
      if ( nJets2 < cut_minNJet || nBjets2 < cut_minNBjet ) continue;
      //int nBjet2 = 0;
      //for ( int j2=0, nj=event2.bTags->size(); j2<nj; ++j2 )
      //{
      //  if ( event2.bTags->at(j2) > cut_bTag ) ++nBjet2;
      //}
      //if ( nBjet2 < cut_minNBjet ) continue;

      //cout << "iEvent2 = " << iEvent2 << ": " << nJets2 << " " << nBjets2 << endl;
      break;

    }

    for ( int j1=0, nj1=event1.jets->size(); j1<nj1; ++j1 )
    {
      const LorentzVector jet1 = event1.jets->at(j1);
      const double bTag1 = event1.bTags->at(j1);
      const int mcBit1 = event1.jetMCBits->at(j1);

      //if ( DeltaR(jet1, lepton1) < 0.3 ) continue;
      if ( jet1.pt() < cut_minJetPt or abs(jet1.eta()) > cut_maxJetEta ) continue;
      //if ( jet1.pt() < cut_minLeadJetPt or abs(jet1.eta()) > cut_maxJetEta ) continue;
      if ( bTag1 > cut_bTag ) continue;
      //if ( mcBit1&3 ) continue;
      //if ( mcBit1 == 1 or mcBit1 == 2 ) continue;
      if (jet1.pt() > cut_minLeadJetPt) ++nLeadJetsInBiEvt;

      for ( int j2=0, nj2=event2.jets->size(); j2<nj2; ++j2 )
      {
        const LorentzVector jet2 = event2.jets->at(j2);
        const double bTag2 = event2.bTags->at(j2);
        const int mcBit2 = event2.jetMCBits->at(j2);

        //if ( DeltaR(jet2, lepton1) < 0.3 ) continue;
        if ( jet2.pt() < cut_minJetPt or abs(jet2.eta()) > cut_maxJetEta ) continue;
        //if ( jet2.pt() < cut_minLeadJetPt or abs(jet2.eta()) > cut_maxJetEta ) continue;
        if ( bTag2 > cut_bTag ) continue;
        //if ( mcBit2&3 ) continue;
        //if ( mcBit2 == 1 or mcBit2 == 2 ) continue;

        ++nBiEvent;
        if ( DeltaR(jet1, jet2) < cut_drBEvt ) continue;
        ++nPassedBiEvent;

  //if (jet2.pt() > cut_minLeadJetPt) ++nLeadJetsInBiEvt;

        LorentzVector jj = event1.jets->at(j1)+event2.jets->at(j2);
        const double mJJ = jj.mass(); //min(mJJMax-1e-3, jj.mass());

        hBEvt_Mw->Fill(mJJ);
  if(input.Contains("Run2012") && jet1.pt() > cut_minLeadJetPt && jet2.pt() > cut_minLeadJetPt) hBEvt_Mw_JJ->Fill(mJJ);

        if ( mcBit1&4 and mcBit2&4 ) hBEvt_Mw_JJ->Fill(mJJ);
        else if ( (mcBit1&2 and mcBit2&4) or (mcBit1&4 and mcBit2&2) ) hBEvt_Mw_HBJ->Fill(mJJ);
        else if ( (mcBit1&1 and mcBit2&4) or (mcBit1&4 and mcBit2&1) ) hBEvt_Mw_LBJ->Fill(mJJ);
        else if ( mcBit1&3 or mcBit2&3 ) hBEvt_Mw_BX->Fill(mJJ);
        else if ( mcBit1 == 0 and mcBit2 == 0 ) hBEvt_Mw_KK->Fill(mJJ);
        else hBEvt_Mw_JK->Fill(mJJ);

        for ( int j3=0; j3<nj2; ++j3 )
        {
          const LorentzVector jet3 = event2.jets->at(j3);
          const double bTag3 = event2.bTags->at(j3);
          const int mcBit3 = event2.jetMCBits->at(j3);

          //if ( DeltaR(jet3, lepton1) < 0.3 ) continue;
          if ( jet3.pt() < cut_minJetPt or abs(jet3.eta()) > cut_maxJetEta ) continue;
          if ( bTag3 <= cut_bTag ) continue;
          if ( DeltaR(jet1, jet3) < cut_drBEvt ) continue;
          if ( DeltaR(jet2, jet3) < cut_drBEvt ) continue;

    for ( int k=0; k<16; ++k ) {
      if ( jet3.pt() > ptmin[k] and jet3.pt() < ptmax[k] ) {
        SFb2 = 0.938887+(0.00017124* (jet3.pt()))+(-2.76366e-07*(jet3.pt())*(jet3.pt()));
        SFb_up2 = SFb2 + SFb_error[k];
        SFb_dn2 = SFb2 - SFb_error[k];
        break;
      }
    }
    if( jet3.pt() > ptmax[15] ) {
      SFb2 = 0.938887+(0.00017124*ptmax[15])+(-2.76366e-07*(ptmax[15]*ptmax[15]));
      SFb_up2 = SFb2 + 2*SFb_error[15];
      SFb_dn2 = SFb2 - 2*SFb_error[15];
    }

          LorentzVector jjb = jj + event2.jets->at(j3);
          const double mJJB = jjb.mass();
          //double PUweight1 = event1.PUweight;
          //double PUweight2 = event2.PUweight;
          //double PUweight3 = PUweight1*PUweight2; //PUweight3 = 1;
    double SFbPUw3 = SFb2;//*PUweight3; if(input.Contains("Run2012")) SFb2 = SFbPUw3 = 1;
          hBEvt_Mt->Fill(mJJB,SFbPUw3);

    if(input.Contains("Run2012") && jet1.pt() > cut_minLeadJetPt && jet2.pt() > cut_minLeadJetPt) hBEvt_Mt_JJHB->Fill(mJJB,SFbPUw3);

          if ( mcBit1&4 and mcBit2&4 ) hBEvt_Mt_JJHB->Fill(mJJB,SFbPUw3);
/*
          if ( (mcBit1&2 and mcBit2&4 and mcBit3&4) or 
               (mcBit1&4 and mcBit2&2 and mcBit3&4) or 
               (mcBit1&4 and mcBit2&4 and mcBit3&2) ) hBEvt_Mt_JJHB->Fill(mJJB,SFbPUw3);
          else if ( (mcBit1&1 and mcBit2&4 and mcBit3&4) or
                    (mcBit1&4 and mcBit2&1 and mcBit3&4) or 
                    (mcBit1&4 and mcBit2&4 and mcBit3&1) ) hBEvt_Mt_JJLB->Fill(mJJB,SFbPUw3);
          else if ( (mcBit1&1 and mcBit2&2 and mcBit3&4) or
                    (mcBit1&2 and mcBit2&1 and mcBit3&4) or
                    (mcBit1&4 and mcBit2&1 and mcBit3&2) ) hBEvt_Mt_HBLBJ->Fill(mJJB,SFbPUw3);
*/
          //else if ( mcBit1+mcBit2 == 0 or mcBit2+mcBit3 == 0 or mcBit3+mcBit1 == 0 ) hBEvt_Mt_KKX->Fill(mJJB);
          //`else if ( mcBit1 ==0 or mcBit2 == 0 or mcBit3 == 0 ) hBEvt_Mt_KXY->Fill(mJJB);
          else hBEvt_Mt_XYZ->Fill(mJJB,SFbPUw3);
    hBEvt_Mt_SFb->Fill(SFb2);
        }
      }
      //break;

    }
    hBEvt_nJets->Fill(nLeadJetsInBiEvt); nLeadJetsInBiEvt = 0;

    //cout << event1.lepton->pt() << ' ' << event2.lepton->pt() << endl;
    //cout << event1.charge << ' ' << event2.charge << endl;
    //cout << event1.bTags->size() << ' ' << event2.bTags->size() << endl;
  }

  // Restore scale factor by skipping overlapping jets
  //hBEvt_Mw->Scale(0.5*nBiEvent/nPassedBiEvent);
  //hBEvt_Mw->Scale(0.5);
  //hBEvt_Mt->Scale(1./3);

  outFile->Write();

}
