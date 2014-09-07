#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/RotationZ.h"
#include "Math/VectorUtil_Cint.h"

using namespace std;
using namespace ROOT::Math::VectorUtil;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

const double mW_world = 80.387;
const double mW_error = 0.019;

//const int cut_minNJet = 4;
//const double cut_minJetPt = 35;
//const double cut_maxJetEta = 2.5;
const double cut_maxJetEta = 2.4; //since Dec 24, 2013 (twiki for SFb)
const double cut_minLeadJetPt = 35; //considered to be disabled from Dec 30, 2013
const double cut_maxMuonEta = 0.8;

//const int cut_minNBjet = 2;
const double cut_bTag = 0.679; //medium
//const double cut_bTag = 0.898; //Tight
//const double cut_bTag = 0.244; //Low

const double cut_drSEvt = 0.4; //for SameEvent
const double cut_drBEvt = 0.5; //for BiEvent (used to be 0.4)

struct EventTopology
{
  int run, lumi, event;
  double PUweight;
  double Ptweight;
  LorentzVector* lepton;
  LorentzVector* met;
  std::vector<LorentzVector>* jets;
  int charge;
  int vertex;
  std::vector<double>* bTags;
  std::vector<int>* jetMCBits;
//  std::vector<double>* pdgGen;
  //  std::vector<double>* pTGen;

  EventTopology()
  {
    lepton = new LorentzVector;
    met = new LorentzVector;
    jets = new std::vector<LorentzVector>();
    bTags = new std::vector<double>();
    jetMCBits = new std::vector<int>();
//    pdgGen = new std::vector<double>();
    //    pTGen = new std::vector<double>();
  };
};

void setBranch(TTree* tree, EventTopology& e)
{
  tree->SetBranchAddress("run"  , &e.run  );
  tree->SetBranchAddress("lumi" , &e.lumi );
  tree->SetBranchAddress("event", &e.event);
  tree->SetBranchAddress("PUweight", &e.PUweight);
  tree->SetBranchAddress("Ptweight", &e.Ptweight);

  tree->SetBranchAddress("lepton", &e.lepton);
  tree->SetBranchAddress("met", &e.met);
  tree->SetBranchAddress("jets", &e.jets);

  tree->SetBranchAddress("charge", &e.charge);
  tree->SetBranchAddress("bTags", &e.bTags);
  tree->SetBranchAddress("jetMCBits", &e.jetMCBits);
  tree->SetBranchAddress("vertex", &e.vertex);
//  tree->SetBranchAddress("pdgGen", &e.pdgGen);
  //  tree->SetBranchAddress("pTGen", &e.pTGen);
}

void mix(TString input = "MSDecays_172", int cutNj = 4, int cutPt = 35, int cutNb = 2, int cutQ = 0)
{
  cout << endl << " ************************************************* " << endl
    //<< endl << " WARNING! PUweight = 1, so need to update PUweight " << endl
    << endl << " charge = " << cutQ << endl
    << endl << " maxJetEta = " <<cut_maxJetEta<<" (used to be 2.5) " << endl
    << endl << " cut_bTag = " << cut_bTag << ", cutNb = " << cutNb <<endl
    << endl << " ************************************************* " << endl << endl; 
  const char* sample = input;
  const int cut_minJetPt = cutPt;
  const int cut_minNBjet = cutNb;
  const int cut_minNJet = cutNj;
  //const int cut_minNLightJet = cutNj - cutNb; if(cut_minNLightJet<=0) { cout << "Error! checking nLightJets requirement." << endl; return; }

  TFile* file1 = TFile::Open(Form("../ntuples/result_%s.root", sample));
  TFile* file2 = TFile::Open(Form("../ntuples/result_%s.root", sample));
  //TFile* file1 = TFile::Open(Form("ehkwon_ntuples/result_%s.root", sample));
  //TFile* file2 = TFile::Open(Form("ehkwon_ntuples/result_%s.root", sample));
  TTree* tree1 = (TTree*)file1->Get("event/tree");//eventUp, eventDn,eventJERUp,eventJERDn
  TTree* tree2 = (TTree*)file2->Get("event/tree");

  EventTopology event1, event2;

  setBranch(tree1, event1);
  setBranch(tree2, event2);

  //TFile* outFile = TFile::Open(Form("hgood/hgood_%s_pt%d_nj%d_nb%d_lq%d.root", sample, cut_minJetPt, cut_minNJet, cut_minNBjet, cutQ), "RECREATE");
  //TFile* outFile = TFile::Open(Form("ehkwon_hist/hist_%s_pt%d_nj%d_nb%d.root", sample, cut_minJetPt, cut_minNJet, cut_minNBjet), "RECREATE");
  //TFile* outFile = TFile::Open(Form("hist_%s_pt%d_nj%d_nb%d_lq%d_nlj3_lessbin.root", sample, cut_minJetPt, cut_minNJet, cut_minNBjet, cutQ), "RECREATE");
  //TFile* outFile = TFile::Open(Form("hist_%s_pt%d_nj%d_nb%d_lq%d_nlj3_CSVT.root", sample, cut_minJetPt, cut_minNJet, cut_minNBjet, cutQ), "RECREATE");
  TFile* outFile = TFile::Open(Form("hist_%s_pt%d_nj%d_nb%d_lq%d_nlj3.root", sample, cut_minJetPt, cut_minNJet, cut_minNBjet, cutQ), "RECREATE");
  TDirectory* dirSEvt = outFile->mkdir("SEvt");
  TDirectory* dirBEvt = outFile->mkdir("BEvt");

  dirSEvt->cd();
  const double binShift = 0;
  const double mJJMax = 500+binShift;
  //const double binShift2 = 2*binShift;
  const double binShift2 = binShift;
  const double mJJBMax = 1000+ binShift2;
  //const double mJJBMax = mJJMax;

  //const double binShift = 2.5;
  //const double mJJMax = 500+binShift;
  //const double binShift2 = 1;
  const double rMax = 6;

  const int nbin = 100;
  const int nbinR = 50;
  cout << "nBin: " << nbin << ", nBinR: " << nbinR << endl;
  const TString axisTitleMw = ";Dijet mass (GeV/c^{2});Events per 5GeV/c^{2}";
  const TString axisTitleMt = ";Three jet mass (GeV/c^{2});Events per 10GeV/c^{2}";
  //const TString axisTitleMt = ";R32;Events";

  TH1F* hnvertex = new TH1F("hnvertex", "vertex", 50, binShift, 50+binShift);
  TH1F* hSEvt_Mw = new TH1F("hMw", "Dijet mass"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hSEvt_Mw_JJ = new TH1F("hMw_JJ", "J+J"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hSEvt_Mw_JK = new TH1F("hMw_JK", "J+K"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hSEvt_Mw_KK = new TH1F("hMw_KK", "K+K"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hSEvt_Mw_HBJ = new TH1F("hMw_HBJ", "HadB+J"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hSEvt_Mw_LBJ = new TH1F("hMw_LBJ", "LepB+J"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hSEvt_Mw_BX = new TH1F("hMw_BX", "B+X"+axisTitleMw, nbin, binShift, mJJMax);

  TH1F* nhSEvt_Mt = new TH1F("nhMt", "Three jet mass"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* nhSEvt_Mt_JJHB = new TH1F("nhMt_JJHB", "J+J+HadB"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* nhSEvt_Mt_XYZ = new TH1F("nhMt_XYZ", "Others"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* rhSEvt_Mt = new TH1F("rhMt", "Three jet mass;R32;Events", nbinR, 1, rMax);
  TH1F* rhSEvt_Mt_JJHB = new TH1F("rhMt_JJHB", "J+J+HadB;R32;Events", nbinR, 1, rMax);
  TH1F* rhSEvt_Mt_XYZ = new TH1F("rhMt_XYZ", "Others;R32;Events", nbinR, 1, rMax);

  TH2F* rhSEvt_Mw_vsR = new TH2F("rhMw_vsR", "Mw vs R32", nbin, 1, rMax, nbinR, binShift, mJJMax);
  TH2F* rhSEvt_Mt_vsR = new TH2F("rhMt_vsR", "Mt vs R32", nbin, 1, rMax, nbinR, binShift2, mJJBMax);

  TH1F* nhSEvt_Mt_A = new TH1F("nhMt_A", "JJ+X"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* nhSEvt_Mt_B = new TH1F("nhMt_B", "Others"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* rhSEvt_Mt_A = new TH1F("rhMt_A", "JJ+X;R32;Events", nbin, binShift2, mJJBMax);
  TH1F* rhSEvt_Mt_B = new TH1F("rhMt_B", "Others;R32;Events", nbin, binShift2, mJJBMax);

  TH1F* hSEvt_Mt = new TH1F("hMt", "Three jet mass"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hSEvt_Mt_JJHB = new TH1F("hMt_JJHB", "J+J+HadB"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hSEvt_Mt_JJLB = new TH1F("hMt_JJLB", "J+J+LepB"+axisTitleMt, nbin, binShift2, mJJBMax);
  //TH1F* hSEvt_Mt_KXY = new TH1F("hMt_KXY", "K+XY"+axisTitleMt, nbin, binShift2, mJJBMax);
  //TH1F* hSEvt_Mt_KKX = new TH1F("hMt_KKX", "KK+X"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hSEvt_Mt_XYZ = new TH1F("hMt_XYZ", "Others"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hSEvt_Mt_HBLBJ = new TH1F("hMt_HBLBJ", "J+HadB+LepB"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hSEvt_Mt_NEW = new TH1F("hMt_NEW", "JJ+B"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hSEvt_Mt_NEW1 = new TH1F("hMt_NEW1", "JJ+B"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hSEvt_Mt_NEW2 = new TH1F("hMt_NEW2", "JJ+B"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hSEvt_Mt_NEW3 = new TH1F("hMt_NEW3", "JJ+B"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hSEvt_Mt_Z = new TH1F("hMt_Z", "Others"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hSEvt_Mt_A = new TH1F("hMt_A", "JJ+X"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hSEvt_Mt_B = new TH1F("hMt_B", "Others"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hSEvt_Mt_SFb = new TH1F("hMt_SFb", "SFb", 100, 0.9, 1);
  TH1F* hSEvt_nLeadJets = new TH1F("hnleadjets","hnleadjets",10,0,10);
  TH1F* hSEvt_nJets = new TH1F("hnjets","hnjets",15,0,15);
  TH1F* hSEvt_nBJets = new TH1F("hnbjets","hnbjets",10,0,10);
  TH1F* hSEvt_nLightJets = new TH1F("hnlightjets","hnlightjets",10,0,10);
  TH1F* hSEvt_B_eta = new TH1F("hbeta","hbeta",100,-3,3);
  TH1F* hSEvt_B_phi = new TH1F("hbphi","hbphi",100,-TMath::Pi(),TMath::Pi());
  TH1F* hSEvt_B_pt = new TH1F("hbpt","hbpt",100,30,330);

  double jet2_ptbins[48]={0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,260,280,300};
  TH1F* hSEvt_Mu_pt = new TH1F("Mu_pt","Mu_pt",100,0,500);
  TH1F* hSEvt_Mu_eta = new TH1F("Mu_eta","Mu_eta",25, -2.5, 2.5);
  TH1F* hSEvt_LJets_pt = new TH1F("hSEvt_LJets_pt","hSEvt_LJets_pt",100,0,500);
  TH1F* hSEvt_LJets_eta = new TH1F("hSEvt_LJets_eta","hSEvt_LJets_eta",24,-2.4,2.4);
  TH1F* hSEvt_L2Jets_pt = new TH1F("hSEvt_L2Jets_pt","hSEvt_L2Jets_pt",100,0,500);
  TH1F* hSEvt_L2Jets_pt1 = new TH1F("hSEvt_L2Jets_pt1","hSEvt_L2Jets_pt1",47,jet2_ptbins);
  TH1F* hSEvt_L2Jets_eta = new TH1F("hSEvt_L2Jets_eta","hSEvt_L2Jets_eta",24,-2.4,2.4);
  TH1F* hSEvt_L3Jets_pt = new TH1F("hSEvt_L3Jets_pt","hSEvt_L3Jets_pt",100,0,500);
  TH1F* hSEvt_L3Jets_eta = new TH1F("hSEvt_L3Jets_eta","hSEvt_L3Jets_eta",24,-2.4,2.4);
  TH1F* hSEvt_L4Jets_pt = new TH1F("hSEvt_L4Jets_pt","hSEvt_L4Jets_pt",100,0,500);
  TH1F* hSEvt_L4Jets_eta = new TH1F("hSEvt_L4Jets_eta","hSEvt_L4Jets_eta",24,-2.4,2.4);


  TH1F* hSEvt_pt_Mw = new TH1F("hMw_pt_Mw", "J+J",100,0,500);
  TH1F* hSEvt_eta_Mw = new TH1F("hMw_eta_Mw", "J+J",100,-3,3);
  TH1F* hSEvt_phi_Mw = new TH1F("hMw_phi_Mw", "J+J",100,-TMath::Pi(),TMath::Pi());
  TH1F* hSEvt_dR_Mw = new TH1F("hMw_dR_Mw", "J+J",60,0,6);
  TH1F* hSEvt_dPhi_Mw = new TH1F("hMw_dPhi_Mw", "J+J",100,-3,3);
  TH1F* hSEvt_dEta_Mw = new TH1F("hMw_dEta_Mw", "J+J",100,-6,6);
  TH1F* hSEvt_pt_JJ = new TH1F("hMw_pt_JJ", "J+J",100,0,500);
  TH1F* hSEvt_eta_JJ = new TH1F("hMw_eta_JJ", "J+J",100,-3,3);
  TH1F* hSEvt_phi_JJ = new TH1F("hMw_phi_JJ", "J+J",100,-TMath::Pi(),TMath::Pi());
  TH1F* hSEvt_dR_JJ = new TH1F("hMw_dR_JJ", "J+J",60,0,6);
  TH1F* hSEvt_dPhi_JJ = new TH1F("hMw_dPhi_JJ", "J+J",100,-3,3);
  TH1F* hSEvt_dEta_JJ = new TH1F("hMw_dEta_JJ", "J+J",100,-6,6);
  TH1F* hSEvt_pt_JK = new TH1F("hMw_pt_JK", "J+J",100,0,500);
  TH1F* hSEvt_eta_JK = new TH1F("hMw_eta_JK", "J+J",100,-3,3);
  TH1F* hSEvt_phi_JK = new TH1F("hMw_phi_JK", "J+J",100,-TMath::Pi(),TMath::Pi());
  TH1F* hSEvt_dR_JK = new TH1F("hMw_dR_JK", "J+J",60,0,6);
  TH1F* hSEvt_dPhi_JK = new TH1F("hMw_dPhi_JK", "J+J",100,-3,3);
  TH1F* hSEvt_dEta_JK = new TH1F("hMw_dEta_JK", "J+J",100,-6,6);
  TH1F* hSEvt_pt_KK = new TH1F("hMw_pt_KK", "J+J",100,0,500);
  TH1F* hSEvt_eta_KK = new TH1F("hMw_eta_KK", "J+J",100,-3,3);
  TH1F* hSEvt_phi_KK = new TH1F("hMw_phi_KK", "J+J",100,-TMath::Pi(),TMath::Pi());
  TH1F* hSEvt_dR_KK = new TH1F("hMw_dR_KK", "J+J",60,0,6);
  TH1F* hSEvt_dPhi_KK = new TH1F("hMw_dPhi_KK", "J+J",100,-3,3);
  TH1F* hSEvt_dEta_KK = new TH1F("hMw_dEta_KK", "J+J",100,-6,6);
  TH1F* hSEvt_pt_HBJ = new TH1F("hMw_pt_HBJ", "J+J",100,0,500);
  TH1F* hSEvt_eta_HBJ = new TH1F("hMw_eta_HBJ", "J+J",100,-3,3);
  TH1F* hSEvt_phi_HBJ = new TH1F("hMw_phi_HBJ", "J+J",100,-TMath::Pi(),TMath::Pi());
  TH1F* hSEvt_dR_HBJ = new TH1F("hMw_dR_HBJ", "J+J",60,0,6);
  TH1F* hSEvt_dPhi_HBJ = new TH1F("hMw_dPhi_HBJ", "J+J",100,-3,3);
  TH1F* hSEvt_dEta_HBJ = new TH1F("hMw_dEta_HBJ", "J+J",100,-6,6);
  TH1F* hSEvt_pt_LBJ = new TH1F("hMw_pt_LBJ", "J+J",100,0,500);
  TH1F* hSEvt_eta_LBJ = new TH1F("hMw_eta_LBJ", "J+J",100,-3,3);
  TH1F* hSEvt_phi_LBJ = new TH1F("hMw_phi_LBJ", "J+J",100,-TMath::Pi(),TMath::Pi());
  TH1F* hSEvt_dR_LBJ = new TH1F("hMw_dR_LBJ", "J+J",60,0,6);
  TH1F* hSEvt_dPhi_LBJ = new TH1F("hMw_dPhi_LBJ", "J+J",100,-3,3);
  TH1F* hSEvt_dEta_LBJ = new TH1F("hMw_dEta_LBJ", "J+J",100,-6,6);
  TH1F* hSEvt_pt_BX = new TH1F("hMw_pt_BX", "J+J",100,0,500);
  TH1F* hSEvt_eta_BX = new TH1F("hMw_eta_BX", "J+J",100,-3,3);
  TH1F* hSEvt_phi_BX = new TH1F("hMw_phi_BX", "J+J",100,-TMath::Pi(),TMath::Pi());
  TH1F* hSEvt_dR_BX = new TH1F("hMw_dR_BX", "J+J",60,0,6);
  TH1F* hSEvt_dPhi_BX = new TH1F("hMw_dPhi_BX", "J+J",100,-3,3);
  TH1F* hSEvt_dEta_BX = new TH1F("hMw_dEta_BX", "J+J",100,-6,6);

  TH1F* hSEvt_dR_LK = new TH1F("hMw_dR_LK", "J+J",60,0,6);
  TH1F* hSEvt_dPhi_LK = new TH1F("hMw_dPhi_LK", "J+J",100,-3,3);

  TH1F* hSEvt_pt_Mt = new TH1F("hMt_pt_Mt", "J+J+HadB",100,0,500);
  TH1F* hSEvt_eta_Mt = new TH1F("hMt_eta_Mt", "J+J+HadB",100,-3,3);
  TH1F* hSEvt_phi_Mt = new TH1F("hMt_phi_Mt", "J+J+HadB",100,-TMath::Pi(),TMath::Pi());
  TH1F* hSEvt_dR_Mt = new TH1F("hMt_dR_Mt", "J+J+HadB",60,0,6);
  TH1F* hSEvt_dPhi_Mt = new TH1F("hMt_dPhi_Mt", "J+J+HadB",100,-3,3);
  TH1F* hSEvt_dEta_Mt = new TH1F("hMt_dEta_Mt", "J+J+HadB",100,-6,6);
  TH1F* hSEvt_pt_JJHB = new TH1F("hMt_pt_JJHB", "J+J+HadB",100,0,500);
  TH1F* hSEvt_eta_JJHB = new TH1F("hMt_eta_JJHB", "J+J+HadB",100,-3,3);
  TH1F* hSEvt_phi_JJHB = new TH1F("hMt_phi_JJHB", "J+J+HadB",100,-TMath::Pi(),TMath::Pi());
  TH1F* hSEvt_dR_JJHB = new TH1F("hMt_dR_JJHB", "J+J+HadB",60,0,6);
  TH1F* hSEvt_dPhi_JJHB = new TH1F("hMt_dPhi_JJHB", "J+J+HadB",100,-3,3);
  TH1F* hSEvt_dEta_JJHB = new TH1F("hMt_dEta_JJHB", "J+J+HadB",100,-6,6);
  TH1F* hSEvt_pt_JJLB = new TH1F("hMt_pt_JJLB", "J+J+LepB",100,0,500);
  TH1F* hSEvt_eta_JJLB = new TH1F("hMt_eta_JJLB", "J+J+LepB",100,-3,3);
  TH1F* hSEvt_phi_JJLB = new TH1F("hMt_phi_JJLB", "J+J+LepB",100,-TMath::Pi(),TMath::Pi());
  TH1F* hSEvt_dR_JJLB = new TH1F("hMt_dR_JJLB", "J+J+LepB",60,0,6);
  TH1F* hSEvt_dPhi_JJLB = new TH1F("hMt_dPhi_JJLB", "J+J+LepB",100,-3,3);
  TH1F* hSEvt_dEta_JJLB = new TH1F("hMt_dEta_JJLB", "J+J+LepB",100,-6,6);
  TH1F* hSEvt_pt_HBLBJ = new TH1F("hMt_pt_HBLBJ", "J+HadB+LepB",100,0,500);
  TH1F* hSEvt_eta_HBLBJ = new TH1F("hMt_eta_HBLBJ", "J+HadB+LepB",100,-3,3);
  TH1F* hSEvt_phi_HBLBJ = new TH1F("hMt_phi_HBLBJ", "J+HadB+LepB",100,-TMath::Pi(),TMath::Pi());
  TH1F* hSEvt_dR_HBLBJ = new TH1F("hMt_dR_HBLBJ", "J+HadB+LepB",60,0,6);
  TH1F* hSEvt_dPhi_HBLBJ = new TH1F("hMt_dPhi_HBLBJ", "J+HadB+LepB",100,-3,3);
  TH1F* hSEvt_dEta_HBLBJ = new TH1F("hMt_dEta_HBLBJ", "J+HadB+LepB",100,-6,6);
  TH1F* hSEvt_pt_XYZ = new TH1F("hMt_pt_XYZ", "J+HadB+LepB",100,0,500);
  TH1F* hSEvt_eta_XYZ = new TH1F("hMt_eta_XYZ", "J+HadB+LepB",100,-3,3);
  TH1F* hSEvt_phi_XYZ = new TH1F("hMt_phi_XYZ", "J+HadB+LepB",100,-TMath::Pi(),TMath::Pi());
  TH1F* hSEvt_dR_XYZ = new TH1F("hMt_dR_XYZ", "J+HadB+LepB",60,0,6);
  TH1F* hSEvt_dPhi_XYZ = new TH1F("hMt_dPhi_XYZ", "J+HadB+LepB",100,-3,3);
  TH1F* hSEvt_dEta_XYZ = new TH1F("hMt_dEta_XYZ", "J+HadB+LepB",100,-6,6);

  TH2F* hDPhiDPhi_JKLK = new TH2F("hDPhiDPhi_JKLK", "#DeltaPhi vs #DeltaPhi;#Delta#phi(J,K);#Delta#phi(L,K)", 100, -TMath::Pi(), TMath::Pi(), 100, TMath::Pi(), TMath::Pi());

  dirBEvt->cd();
  TH1F* hBEvt_Mw = new TH1F("hMw", "Dijet mass"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hBEvt_Mw_JJ = new TH1F("hMw_JJ", "J+J"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hBEvt_Mw_JK = new TH1F("hMw_JK", "J+K"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hBEvt_Mw_KK = new TH1F("hMw_KK", "K+K"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hBEvt_Mw_HBJ = new TH1F("hMw_HBJ", "HadB+J"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hBEvt_Mw_LBJ = new TH1F("hMw_LBJ", "LepB+J"+axisTitleMw, nbin, binShift, mJJMax);
  TH1F* hBEvt_Mw_BX = new TH1F("hMw_BX", "B+X"+axisTitleMw, nbin, binShift, mJJMax);

  TH1F* nhBEvt_Mt = new TH1F("nhMt", "Three jet mass"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* nhBEvt_Mt_JJHB = new TH1F("nhMt_JJHB", "J+J+HadB"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* nhBEvt_Mt_XYZ = new TH1F("nhMt_XYZ", "Others"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* rhBEvt_Mt = new TH1F("rhMt", "Three jet mass;R32;Events", nbinR, 1, rMax);
  TH1F* rhBEvt_Mt_JJHB = new TH1F("rhMt_JJHB", "J+J+HadB;R32;Events", nbinR, 1, rMax);
  TH1F* rhBEvt_Mt_XYZ = new TH1F("rhMt_XYZ", "Others;R32;Events", nbinR, 1, rMax);

  TH2F* rhBEvt_Mw_vsR = new TH2F("rhMw_vsR", "Mw vs R32", nbin, 1, rMax, nbinR, binShift, mJJMax);
  TH2F* rhBEvt_Mt_vsR = new TH2F("rhMt_vsR", "Mt vs R32", nbin, 1, rMax, nbinR, binShift2, mJJBMax);

  TH1F* hBEvt_Mt_A = new TH1F("hMt_A", "JJ+X"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hBEvt_Mt_B = new TH1F("hMt_B", "Others"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* nhBEvt_Mt_A = new TH1F("nhMt_A", "JJ+X"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* nhBEvt_Mt_B = new TH1F("nhMt_B", "Others"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* rhBEvt_Mt_A = new TH1F("rhMt_A", "JJ+X;R32;Events", nbin, binShift2, mJJBMax);
  TH1F* rhBEvt_Mt_B = new TH1F("rhMt_B", "Others;R32;Events", nbin, binShift2, mJJBMax);

  TH1F* hBEvt_Mt = new TH1F("hMt", "Three jet mass"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hBEvt_Mt_JJHB = new TH1F("hMt_JJHB", "J+J+HadB"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hBEvt_Mt_JJLB = new TH1F("hMt_JJLB", "J+J+LepB"+axisTitleMt, nbin, binShift2, mJJBMax);
  //TH1F* hBEvt_Mt_KXY = new TH1F("hMt_KXY", "K+XY"+axisTitleMt, nbin, binShift2, mJJBMax);
  //TH1F* hBEvt_Mt_KKX = new TH1F("hMt_KKX", "KK+X"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hBEvt_Mt_XYZ = new TH1F("hMt_XYZ", "Others"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hBEvt_Mt_HBLBJ = new TH1F("hMt_HBLBJ", "J+HadB+LepB"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hBEvt_Mt_NEW = new TH1F("hMt_NEW", "JJ+B"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hBEvt_Mt_OLD = new TH1F("hMt_OLD", "J+HadB+LepB"+axisTitleMt, nbin, binShift2, mJJBMax);
  TH1F* hBEvt_Mt_SFb = new TH1F("hMt_SFb", "SFb", 100, 0.9, 1);
  TH1F* hBEvt_nLeadJets = new TH1F("hnleadjets","hnleadjets",10,0,10);
  TH1F* hBEvt_nJets = new TH1F("hnjets","hnjets",15,0,15);
  TH1F* hBEvt_nBJets = new TH1F("hnbjets","hnbjets",10,0,10);

  TH1F* hBEvt_pt_Mw = new TH1F("hMw_pt_Mw", "J+J",100,0,500);
  TH1F* hBEvt_eta_Mw = new TH1F("hMw_eta_Mw", "J+J",100,-3,3);
  TH1F* hBEvt_phi_Mw = new TH1F("hMw_phi_Mw", "J+J",100,-TMath::Pi(),TMath::Pi());
  TH1F* hBEvt_dR_Mw = new TH1F("hMw_dR_Mw", "J+J",60,0,6);
  TH1F* hBEvt_dPhi_Mw = new TH1F("hMw_dPhi_Mw", "J+J",100,-3,3);
  TH1F* hBEvt_dEta_Mw = new TH1F("hMw_dEta_Mw", "J+J",100,-6,6);
  TH1F* hBEvt_pt_JJ = new TH1F("hMw_pt_JJ", "J+J",100,0,500);
  TH1F* hBEvt_eta_JJ = new TH1F("hMw_eta_JJ", "J+J",100,-3,3);
  TH1F* hBEvt_phi_JJ = new TH1F("hMw_phi_JJ", "J+J",100,-TMath::Pi(),TMath::Pi());
  TH1F* hBEvt_dR_JJ = new TH1F("hMw_dR_JJ", "J+J",60,0,6);
  TH1F* hBEvt_dPhi_JJ = new TH1F("hMw_dPhi_JJ", "J+J",100,-3,3);
  TH1F* hBEvt_dEta_JJ = new TH1F("hMw_dEta_JJ", "J+J",100,-6,6);
  TH1F* hBEvt_pt_JK = new TH1F("hMw_pt_JK", "J+J",100,0,500);
  TH1F* hBEvt_eta_JK = new TH1F("hMw_eta_JK", "J+J",100,-3,3);
  TH1F* hBEvt_phi_JK = new TH1F("hMw_phi_JK", "J+J",100,-TMath::Pi(),TMath::Pi());
  TH1F* hBEvt_dR_JK = new TH1F("hMw_dR_JK", "J+J",60,0,6);
  TH1F* hBEvt_dPhi_JK = new TH1F("hMw_dPhi_JK", "J+J",100,-3,3);
  TH1F* hBEvt_dEta_JK = new TH1F("hMw_dEta_JK", "J+J",100,-6,6);
  TH1F* hBEvt_pt_KK = new TH1F("hMw_pt_KK", "J+J",100,0,500);
  TH1F* hBEvt_eta_KK = new TH1F("hMw_eta_KK", "J+J",100,-3,3);
  TH1F* hBEvt_phi_KK = new TH1F("hMw_phi_KK", "J+J",100,-TMath::Pi(),TMath::Pi());
  TH1F* hBEvt_dR_KK = new TH1F("hMw_dR_KK", "J+J",60,0,6);
  TH1F* hBEvt_dPhi_KK = new TH1F("hMw_dPhi_KK", "J+J",100,-3,3);
  TH1F* hBEvt_dEta_KK = new TH1F("hMw_dEta_KK", "J+J",100,-6,6);
  TH1F* hBEvt_pt_HBJ = new TH1F("hMw_pt_HBJ", "J+J",100,0,500);
  TH1F* hBEvt_eta_HBJ = new TH1F("hMw_eta_HBJ", "J+J",100,-3,3);
  TH1F* hBEvt_phi_HBJ = new TH1F("hMw_phi_HBJ", "J+J",100,-TMath::Pi(),TMath::Pi());
  TH1F* hBEvt_dR_HBJ = new TH1F("hMw_dR_HBJ", "J+J",60,0,6);
  TH1F* hBEvt_dPhi_HBJ = new TH1F("hMw_dPhi_HBJ", "J+J",100,-3,3);
  TH1F* hBEvt_dEta_HBJ = new TH1F("hMw_dEta_HBJ", "J+J",100,-6,6);
  TH1F* hBEvt_pt_LBJ = new TH1F("hMw_pt_LBJ", "J+J",100,0,500);
  TH1F* hBEvt_eta_LBJ = new TH1F("hMw_eta_LBJ", "J+J",100,-3,3);
  TH1F* hBEvt_phi_LBJ = new TH1F("hMw_phi_LBJ", "J+J",100,-TMath::Pi(),TMath::Pi());
  TH1F* hBEvt_dR_LBJ = new TH1F("hMw_dR_LBJ", "J+J",60,0,6);
  TH1F* hBEvt_dPhi_LBJ = new TH1F("hMw_dPhi_LBJ", "J+J",100,-3,3);
  TH1F* hBEvt_dEta_LBJ = new TH1F("hMw_dEta_LBJ", "J+J",100,-6,6);
  TH1F* hBEvt_pt_BX = new TH1F("hMw_pt_BX", "J+J",100,0,500);
  TH1F* hBEvt_eta_BX = new TH1F("hMw_eta_BX", "J+J",100,-3,3);
  TH1F* hBEvt_phi_BX = new TH1F("hMw_phi_BX", "J+J",100,-TMath::Pi(),TMath::Pi());
  TH1F* hBEvt_dR_BX = new TH1F("hMw_dR_BX", "J+J",60,0,6);
  TH1F* hBEvt_dPhi_BX = new TH1F("hMw_dPhi_BX", "J+J",100,-3,3);
  TH1F* hBEvt_dEta_BX = new TH1F("hMw_dEta_BX", "J+J",100,-6,6);

  TH1F* hBEvt_pt_Mt = new TH1F("hMt_pt_Mt", "J+J+HadB",100,0,500);
  TH1F* hBEvt_eta_Mt = new TH1F("hMt_eta_Mt", "J+J+HadB",100,-3,3);
  TH1F* hBEvt_phi_Mt = new TH1F("hMt_phi_Mt", "J+J+HadB",100,-TMath::Pi(),TMath::Pi());
  TH1F* hBEvt_dR_Mt = new TH1F("hMt_dR_Mt", "J+J+HadB",60,0,6);
  TH1F* hBEvt_dEta_Mt = new TH1F("hMt_dEta_Mt", "J+J+HadB",100,-6,6);
  TH1F* hBEvt_dPhi_Mt = new TH1F("hMt_dPhi_Mt", "J+J+HadB",100,-3,3);
  TH1F* hBEvt_pt_JJHB = new TH1F("hMt_pt_JJHB", "J+J+HadB",100,0,500);
  TH1F* hBEvt_eta_JJHB = new TH1F("hMt_eta_JJHB", "J+J+HadB",100,-3,3);
  TH1F* hBEvt_phi_JJHB = new TH1F("hMt_phi_JJHB", "J+J+HadB",100,-TMath::Pi(),TMath::Pi());
  TH1F* hBEvt_dR_JJHB = new TH1F("hMt_dR_JJHB", "J+J+HadB",60,0,6);
  TH1F* hBEvt_dPhi_JJHB = new TH1F("hMt_dPhi_JJHB", "J+J+HadB",100,-3,3);
  TH1F* hBEvt_dEta_JJHB = new TH1F("hMt_dEta_JJHB", "J+J+HadB",100,-6,6);
  TH1F* hBEvt_pt_JJLB = new TH1F("hMt_pt_JJLB", "J+J+LepB",100,0,500);
  TH1F* hBEvt_eta_JJLB = new TH1F("hMt_eta_JJLB", "J+J+LepB",100,-3,3);
  TH1F* hBEvt_phi_JJLB = new TH1F("hMt_phi_JJLB", "J+J+LepB",100,-TMath::Pi(),TMath::Pi());
  TH1F* hBEvt_dR_JJLB = new TH1F("hMt_dR_JJLB", "J+J+LepB",60,0,6);
  TH1F* hBEvt_dPhi_JJLB = new TH1F("hMt_dPhi_JJLB", "J+J+LepB",100,-3,3);
  TH1F* hBEvt_dEta_JJLB = new TH1F("hMt_dEta_JJLB", "J+J+LepB",100,-6,6);
  TH1F* hBEvt_pt_HBLBJ = new TH1F("hMt_pt_HBLBJ", "J+HadB+LepB",100,0,500);
  TH1F* hBEvt_eta_HBLBJ = new TH1F("hMt_eta_HBLBJ", "J+HadB+LepB",100,-3,3);
  TH1F* hBEvt_phi_HBLBJ = new TH1F("hMt_phi_HBLBJ", "J+HadB+LepB",100,-TMath::Pi(),TMath::Pi());
  TH1F* hBEvt_dR_HBLBJ = new TH1F("hMt_dR_HBLBJ", "J+HadB+LepB",60,0,6);
  TH1F* hBEvt_dPhi_HBLBJ = new TH1F("hMt_dPhi_HBLBJ", "J+HadB+LepB",100,-3,3);
  TH1F* hBEvt_dEta_HBLBJ = new TH1F("hMt_dEta_HBLBJ", "J+HadB+LepB",100,-6,6);
  TH1F* hBEvt_pt_XYZ = new TH1F("hMt_pt_XYZ", "J+HadB+LepB",100,0,500);
  TH1F* hBEvt_eta_XYZ = new TH1F("hMt_eta_XYZ", "J+HadB+LepB",100,-3,3);
  TH1F* hBEvt_phi_XYZ = new TH1F("hMt_phi_XYZ", "J+HadB+LepB",100,-TMath::Pi(),TMath::Pi());
  TH1F* hBEvt_dR_XYZ = new TH1F("hMt_dR_XYZ", "J+HadB+LepB",60,0,6);
  TH1F* hBEvt_dPhi_XYZ = new TH1F("hMt_dPhi_XYZ", "J+HadB+LepB",100,-3,3);
  TH1F* hBEvt_dEta_XYZ = new TH1F("hMt_dEta_XYZ", "J+HadB+LepB",100,-6,6);


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

  float muEtaMin[] = {-2.1, -1.6, -1.2, -0.9, -0.6, -0.3, -0.2, 0.2, 0.3, 0.6, 0.9, 1.2, 1.6};
  float muEtaMax[] = {-1.6, -1.2, -0.9, -0.6, -0.3, -0.2,  0.2, 0.3, 0.6, 0.9, 1.2, 1.6, 2.1}; 

  double mutrigEff[] = {
    0.998523,
    0.993233,
    0.964378,
    0.981626,
    0.986639,
    0.946682,
    0.983689,
    0.961673,
    0.985527,
    0.985813,
    0.965818,
    0.971026,
    1.02413};

  double muidEff[] = {
    0.989381,
    0.998003,
    0.993071,
    0.994778,
    0.995176,
    0.984401,
    0.992150,
    0.984043,
    0.995343,
    0.993261,
    0.992932,
    0.999929,
    0.997083};

  double muisoEff[] = {
    1.005630,
    1.000320,
    0.999505,
    0.997581,
    0.995127,
    0.992673,
    0.992707,
    0.994813,
    0.996828,
    0.999479,
    1.001430,
    0.999707,
    1.004670};

  double muidiso1 = 0;
  double muidiso2 = 0;
  double muScale1 = 0;
  double muScale2 = 0;

  for(int s= 0; s<13; s++)
  { 
    cout << muidEff[s]*muisoEff[s]*mutrigEff[s] <<endl;
  }

  int nLeadJetsInSameEvt = 0, nLeadJetsInBiEvt = 0;
  int nBiEvent = 0, nPassedBiEvent = 0; // Variables to restore overlap removal scale
  for ( int iEvent1=0, nEvent=tree1->GetEntries(); iEvent1<nEvent; ++iEvent1 )
  {
    tree1->GetEntry(iEvent1);
    const LorentzVector& lepton1 = *event1.lepton;
    if( cutQ != 0 && event1.charge != cutQ ) continue;

    //if ( event1.jets->size() < cut_minNJet ) continue;
    //if ( event1.jets->at(0).pt() < cut_minLeadJetPt or event1.jets->at(1).pt() < cut_minLeadJetPt ) continue;
    int nJets1 = 0, nBjets1 = 0, nLightJets1 = 0;
    for ( int j1=0, nj=event1.jets->size(); j1<nj; ++j1 )
    {
      const LorentzVector jet1 = event1.jets->at(j1);
      const double bTag1 = event1.bTags->at(j1);
      //if ( event1.jets->at(0).pt() < cut_minLeadJetPt or event1.jets->at(1).pt() < cut_minLeadJetPt ) continue;
      ////if ( jet1.pt() > cut_minJetPt or abs(jet1.eta()) < cut_maxJetEta ) ++nJets1; //a bug found at Jan 28, 2014
      if ( jet1.pt() > cut_minJetPt and abs(jet1.eta()) < cut_maxJetEta ) 
      {
        ++nJets1;
        if ( bTag1 > cut_bTag ) { ++nBjets1; hSEvt_B_eta->Fill(jet1.eta()); hSEvt_B_phi->Fill(jet1.phi()); hSEvt_B_pt->Fill(jet1.pt()); } 
        else ++nLightJets1;
      }
    }
    //if ( nJets1 < cut_minNJet || nBjets1 < cut_minNBjet || nLightJets1 < cut_minNLightJet ) continue;
    if ( nJets1 < cut_minNJet || nBjets1 < cut_minNBjet ) continue;
    //if ( nLightJets1 < 2 || nBjets1 < 2 ) continue;
    //    if ( nLightJets1 != 2 || nBjets1 != 2 ) continue;   
    //    if ( abs(lepton1.eta()) >= cut_maxMuonEta ) continue;
    //    if ( abs(lepton1.eta()) < cut_maxMuonEta ) continue;
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
      double jet1Eta = jet1.eta();
      double PUweight1 = event1.PUweight; //PUweight1 = 1;
      double Ptweight1 = event1.Ptweight;

      for ( int l1=0; l1<13; ++l1 ) 
      {
        if ( lepton1.eta() >= muEtaMin[l1] and lepton1.eta() < muEtaMax[l1] ) 
        {
          muScale1 = muidEff[l1]*muisoEff[l1]*mutrigEff[l1];
        }
      }    

      double weight1 = PUweight1*muScale1;

      //if ( DeltaR(jet1, lepton1) < 0.3 ) continue;

      if ( jet1.pt() < cut_minJetPt or abs(jet1.eta()) > cut_maxJetEta ) continue;
      //if ( jet1.pt() < cut_minLeadJetPt or abs(jet1.eta()) > cut_maxJetEta ) continue;
      //if ( mcBit1 == 1 or mcBit1 == 2 ) continue;
      //if ( mcBit1&3 ) continue;
      if ( bTag1 > cut_bTag ) continue;
      hSEvt_LJets_pt->Fill(event1.jets->at(0).pt(),weight1);  
      hSEvt_LJets_eta->Fill(event1.jets->at(0).eta(),weight1);
      hSEvt_L2Jets_pt->Fill(event1.jets->at(1).pt(),weight1);hSEvt_L2Jets_pt1->Fill(event1.jets->at(1).pt(),weight1);
      hSEvt_L2Jets_eta->Fill(event1.jets->at(1).eta(),weight1);
      hSEvt_L3Jets_pt->Fill(event1.jets->at(2).pt(),weight1);
      hSEvt_L3Jets_eta->Fill(event1.jets->at(2).eta(),weight1);
      hSEvt_L4Jets_pt->Fill(event1.jets->at(3).pt(),weight1);
      hSEvt_L4Jets_eta->Fill(event1.jets->at(3).eta(),weight1); 
      //      if (jet1.pt() > cut_minLeadJetPt) ++nLeadJetsInSameEvt;
      //      if (jet1.pt() > cut_minLeadJetPt) { hSEvt_LJets_pt->Fill(jet1.pt(),weight1); hSEvt_LJets_eta->Fill(jet1.eta(),weight1);} 


      for ( int j2=j1+1; j2<nj; ++j2 )
      {
        const LorentzVector jet2 = event1.jets->at(j2);
        const double bTag2 = event1.bTags->at(j2);
        const int mcBit2 = event1.jetMCBits->at(j2);
        double jet2Eta = jet2.eta();

        //if ( DeltaR(jet2, lepton1) < 0.3 ) continue;
        if ( jet2.pt() < cut_minJetPt or abs(jet2.eta()) > cut_maxJetEta ) continue;
        //if ( jet2.pt() < cut_minLeadJetPt or abs(jet2.eta()) > cut_maxJetEta ) continue;
        //if ( mcBit2 == 1 or mcBit2 == 2 ) continue;
        if ( bTag2 > cut_bTag ) continue;
        //if ( mcBit2&3 ) continue;

        if ( DeltaR(jet1, jet2) < cut_drSEvt ) continue;

        LorentzVector jj = event1.jets->at(j1)+event1.jets->at(j2);
        const double mJJ = jj.mass(); //min(mJJMax-1e-3, jj.mass());
        const double etaJJ = jj.eta();
        const double ptJJ = jj.pt();
        const double phiJJ = jj.phi();
        double dRJJ = ROOT::Math::VectorUtil::DeltaR(jet1,jet2);
        double dPhiJJ = ROOT::Math::VectorUtil::DeltaPhi(jet1,jet2);
        double dEtaJJ = jet1Eta-jet2Eta;
        double dRLK = ROOT::Math::VectorUtil::DeltaR(lepton1,jet1);
        double dPhiLK = ROOT::Math::VectorUtil::DeltaPhi(lepton1,jet2); 

        hSEvt_Mw->Fill(mJJ,weight1);
        hSEvt_eta_Mw->Fill(etaJJ,weight1);
        hSEvt_pt_Mw->Fill(ptJJ,weight1);
        hSEvt_phi_Mw->Fill(phiJJ,weight1);
        hSEvt_dR_Mw->Fill(dRJJ,weight1);
        hSEvt_dPhi_Mw->Fill(dPhiJJ,weight1);
        hSEvt_dEta_Mw->Fill(dEtaJJ,weight1);

        if(input.Contains("Run2012") && jet1.pt() > cut_minLeadJetPt && jet2.pt() > cut_minLeadJetPt) 
        { 
          hSEvt_Mw_JJ->Fill(mJJ,weight1);
          hSEvt_eta_JJ->Fill(etaJJ,weight1);
          hSEvt_pt_JJ->Fill(ptJJ,weight1);
          hSEvt_phi_JJ->Fill(phiJJ,weight1);
          hSEvt_dR_JJ->Fill(dRJJ,weight1);
          hSEvt_dPhi_JJ->Fill(dPhiJJ,weight1);
          hSEvt_dEta_JJ->Fill(dEtaJJ,weight1);
        }
        if ( mcBit1&4 and mcBit2&4 ) 
        { 
          hSEvt_Mw_JJ->Fill(mJJ,weight1); 
          hSEvt_eta_JJ->Fill(etaJJ,weight1); 
          hSEvt_pt_JJ->Fill(ptJJ,weight1); 
          hSEvt_phi_JJ->Fill(phiJJ,weight1);
          hSEvt_dR_JJ->Fill(dRJJ,weight1);
          hSEvt_dPhi_JJ->Fill(dPhiJJ,weight1);
          hSEvt_dEta_JJ->Fill(dEtaJJ,weight1);
          hSEvt_dR_LK->Fill(dRLK,weight1);
          hSEvt_dPhi_LK->Fill(dPhiLK,weight1);

        }
        else if ( (mcBit1&2 and mcBit2&4) or (mcBit1&4 and mcBit2&2) ) 
        {
          hSEvt_Mw_HBJ->Fill(mJJ,weight1);
          hSEvt_eta_HBJ->Fill(etaJJ,weight1);
          hSEvt_pt_HBJ->Fill(ptJJ,weight1);
          hSEvt_phi_HBJ->Fill(phiJJ,weight1);
          hSEvt_dR_HBJ->Fill(dRJJ,weight1);
          hSEvt_dPhi_HBJ->Fill(dPhiJJ,weight1);
          hSEvt_dEta_HBJ->Fill(dEtaJJ,weight1);
        }
        else if ( (mcBit1&1 and mcBit2&4) or (mcBit1&4 and mcBit2&1) ) 
        {
          hSEvt_Mw_LBJ->Fill(mJJ,weight1);
          hSEvt_eta_LBJ->Fill(etaJJ,weight1);
          hSEvt_pt_LBJ->Fill(ptJJ,weight1);
          hSEvt_phi_LBJ->Fill(phiJJ,weight1);
          hSEvt_dR_LBJ->Fill(dRJJ,weight1);
          hSEvt_dPhi_LBJ->Fill(dPhiJJ,weight1);
          hSEvt_dEta_LBJ->Fill(dEtaJJ,weight1);
        }
        else if ( mcBit1&3 or mcBit2&3 ) 
        {
          hSEvt_Mw_BX->Fill(mJJ,weight1);
          hSEvt_eta_BX->Fill(etaJJ,weight1);
          hSEvt_pt_BX->Fill(ptJJ,weight1);
          hSEvt_phi_BX->Fill(phiJJ,weight1);
          hSEvt_dR_BX->Fill(dRJJ,weight1);
          hSEvt_dPhi_BX->Fill(dPhiJJ,weight1);
          hSEvt_dEta_BX->Fill(dEtaJJ,weight1);
        }
        else if ( mcBit1 == 0 and mcBit2 == 0 )
        {
          hSEvt_Mw_KK->Fill(mJJ,weight1);
          hSEvt_eta_KK->Fill(etaJJ,weight1);
          hSEvt_pt_KK->Fill(ptJJ,weight1);
          hSEvt_phi_KK->Fill(phiJJ,weight1);
          hSEvt_dR_KK->Fill(dRJJ,weight1);
          hSEvt_dPhi_KK->Fill(dPhiJJ,weight1);
          hSEvt_dEta_KK->Fill(dEtaJJ,weight1);
        }
        else if ( ( mcBit1&4 and mcBit2 == 0) )// or ( mcBit1 == 0 and mcBit2&4) ) 
        { 
          hSEvt_Mw_JK->Fill(mJJ,weight1);
          hSEvt_eta_JK->Fill(etaJJ,weight1);
          hSEvt_pt_JK->Fill(ptJJ,weight1);
          hSEvt_phi_JK->Fill(phiJJ,weight1);
          hSEvt_dR_JK->Fill(dRJJ,weight1);
          hSEvt_dPhi_JK->Fill(dPhiJJ,weight1);
          hSEvt_dEta_JK->Fill(dEtaJJ,weight1);

          hDPhiDPhi_JKLK->Fill(dPhiJJ,dPhiLK);
        }


        for ( int j3=0; j3<nj; ++j3 )
        {
          const LorentzVector jet3 = event1.jets->at(j3);
          const double bTag3 = event1.bTags->at(j3);
          const int mcBit3 = event1.jetMCBits->at(j3);
          double jet3Eta = jet3.eta();

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
          const double etaJJB = jjb.eta();
          const double ptJJB = jjb.pt();
          const double phiJJB = jjb.phi();
          double dRJJB = ROOT::Math::VectorUtil::DeltaR(jj,jet3);
          double dPhiJJB = ROOT::Math::VectorUtil::DeltaPhi(jj,jet3);
          double dEtaJJB = etaJJ-jet3Eta;
          const double mJJB_new = (mJJB/mJJ)*mW_world;
          double SFbPUw1 = SFb*weight1; if(input.Contains("Run2012")) SFb = SFbPUw1 = 1;
          hSEvt_Mt->Fill(mJJB,SFbPUw1); 
          nhSEvt_Mt->Fill(mJJB_new,SFbPUw1); 
          rhSEvt_Mt->Fill(mJJB/mJJ,SFbPUw1);
          hSEvt_eta_Mt->Fill(etaJJB,SFbPUw1);
          hSEvt_pt_Mt->Fill(ptJJB,SFbPUw1);
          hSEvt_phi_Mt->Fill(phiJJB,SFbPUw1);
          hSEvt_dR_Mt->Fill(dRJJB,SFbPUw1);
          hSEvt_dPhi_Mt->Fill(dPhiJJB,SFbPUw1);
          hSEvt_dEta_Mt->Fill(dEtaJJB,SFbPUw1);
          if(input.Contains("Run2012") && jet1.pt() > cut_minLeadJetPt && jet2.pt() > cut_minLeadJetPt) 
          { 
            hSEvt_Mt_JJHB->Fill(mJJB,SFbPUw1); 
            hSEvt_eta_JJHB->Fill(etaJJB,SFbPUw1); 
            hSEvt_pt_JJHB->Fill(ptJJB,SFbPUw1); 
            hSEvt_phi_JJHB->Fill(phiJJB,SFbPUw1); 
            hSEvt_dR_JJHB->Fill(dRJJB,SFbPUw1);
            hSEvt_dPhi_JJHB->Fill(dPhiJJB,SFbPUw1);
            hSEvt_dEta_JJHB->Fill(dEtaJJB,SFbPUw1);
          }

          //          if ( (mcBit1&4 and mcBit2&2 and mcBit3&4) or (mcBit1&4 and mcBit2&4 and mcBit3&2) ) hSEvt_Mt_JJHB->Fill(mJJB,SFbPUw1);
          //          else if ( (mcBit1&4 and mcBit2&1 and mcBit3&4) or (mcBit1&4 and mcBit2&4 and mcBit3&1) ) hSEvt_Mt_JJLB->Fill(mJJB,SFbPUw1);
          //          else if ( (mcBit1&4 and mcBit2 == 0 and mcBit3&4) or (mcBit1&4 and mcBit2&4 and mcBit3 == 0) ) hSEvt_Mt_HBLBJ->Fill(mJJB,SFbPUw1);
          //          else if ( (mcBit1&2 and mcBit2&4 and mcBit3&4) or (mcBit1&1 and mcBit2&4 and mcBit3&4) ) hSEvt_Mt_NEW->Fill(mJJB,SFbPUw1);

          if ( (mcBit1==2 and mcBit2&4 and mcBit3&4) or
              (mcBit1&4 and mcBit2==2 and mcBit3&4) or
              (mcBit1&4 and mcBit2&4 and mcBit3==2) ) hSEvt_Mt_NEW->Fill(mJJB,SFbPUw1);
          else if ( (mcBit1==1 and mcBit2&4 and mcBit3&4) or
              (mcBit1&4 and mcBit2==1 and mcBit3&4) or
              (mcBit1&4 and mcBit2&4 and mcBit3==1) ) hSEvt_Mt_NEW2->Fill(mJJB,SFbPUw1);
          else hSEvt_Mt_Z->Fill(mJJB,SFbPUw1);

          if ( (mcBit1==1 and mcBit2&4 and mcBit3&4) or
              (mcBit1&4 and mcBit2==1 and mcBit3&4) or
              (mcBit1&4 and mcBit2&4 and mcBit3==1) ) hSEvt_Mt_NEW1->Fill(mJJB,SFbPUw1);
          else if ( (mcBit1==2 and mcBit2&4 and mcBit3&4) or
              (mcBit1&4 and mcBit2==2 and mcBit3&4) or 
              (mcBit1&4 and mcBit2&4 and mcBit3==2) ) hSEvt_Mt_NEW3->Fill(mJJB,SFbPUw1);

          if ( (mcBit1&2 and mcBit2&4 and mcBit3&4) or
              (mcBit1&4 and mcBit2&2 and mcBit3&4) or 
              (mcBit1&4 and mcBit2&4 and mcBit3&2) ) { hSEvt_Mt_A->Fill(mJJB,SFbPUw1); nhSEvt_Mt_A->Fill(mJJB_new,SFbPUw1); rhSEvt_Mt_A->Fill(mJJB/mJJ,SFbPUw1); }
          else { hSEvt_Mt_B->Fill(mJJB,SFbPUw1); nhSEvt_Mt_B->Fill(mJJB_new,SFbPUw1); rhSEvt_Mt_B->Fill(mJJB/mJJ,SFbPUw1); }


          if ( (mcBit1&4 and mcBit2&4 and mcBit3&2) ) 
          { 
            hSEvt_Mt_JJHB->Fill(mJJB,SFbPUw1); 
            nhSEvt_Mt_JJHB->Fill(mJJB_new,SFbPUw1); 
            rhSEvt_Mt_JJHB->Fill(mJJB/mJJ,SFbPUw1); 
            hSEvt_eta_JJHB->Fill(etaJJB,SFbPUw1); 
            hSEvt_pt_JJHB->Fill(ptJJB,SFbPUw1); 
            hSEvt_phi_JJHB->Fill(phiJJB,SFbPUw1); 
            hSEvt_dR_JJHB->Fill(dRJJB,SFbPUw1);
            hSEvt_dPhi_JJHB->Fill(dPhiJJB,SFbPUw1);
            hSEvt_dEta_JJHB->Fill(dEtaJJB,SFbPUw1);
          }
          /*
             if ( (mcBit1&2 and mcBit2&4 and mcBit3&4) or 
             (mcBit1&4 and mcBit2&2 and mcBit3&4) or 
             (mcBit1&4 and mcBit2&4 and mcBit3&2) ) { hSEvt_Mt_JJHB->Fill(mJJB,SFbPUw1); nhSEvt_Mt_JJHB->Fill(mJJB_new,SFbPUw1); rhSEvt_Mt_JJHB->Fill(mJJB/mJJ,SFbPUw1); }
           */
             else if ( (mcBit1&1 and mcBit2&4 and mcBit3&4) or
                 (mcBit1&4 and mcBit2&1 and mcBit3&4) or 
                 (mcBit1&4 and mcBit2&4 and mcBit3&1) )
             {
               hSEvt_Mt_JJLB->Fill(mJJB,SFbPUw1);
               hSEvt_eta_JJLB->Fill(mJJB,SFbPUw1);
               hSEvt_pt_JJLB->Fill(mJJB,SFbPUw1);
               hSEvt_phi_JJLB->Fill(phiJJB,SFbPUw1);
               hSEvt_dR_JJLB->Fill(dRJJB,SFbPUw1);
               hSEvt_dPhi_JJLB->Fill(dPhiJJB,SFbPUw1);
               hSEvt_dEta_JJLB->Fill(dEtaJJB,SFbPUw1);
             }
             else if ( (mcBit1&1 and mcBit2&2 and mcBit3&4) or
                 (mcBit1&1 and mcBit2&4 and mcBit3&2) or
                 (mcBit1&2 and mcBit2&1 and mcBit3&4) or
                 (mcBit1&2 and mcBit2&4 and mcBit3&1) or
                 (mcBit1&4 and mcBit2&1 and mcBit3&2) or
                 (mcBit1&4 and mcBit2&2 and mcBit3&1) )
             {
               hSEvt_Mt_HBLBJ->Fill(mJJB,SFbPUw1);
               hSEvt_eta_HBLBJ->Fill(etaJJB,SFbPUw1);
               hSEvt_pt_HBLBJ->Fill(ptJJB,SFbPUw1);
               hSEvt_phi_HBLBJ->Fill(phiJJB,SFbPUw1);
               hSEvt_dR_HBLBJ->Fill(dRJJB,SFbPUw1);
               hSEvt_dPhi_HBLBJ->Fill(dPhiJJB,SFbPUw1);
               hSEvt_dEta_HBLBJ->Fill(dEtaJJB,SFbPUw1);
             }
             //else if ( mcBit1+mcBit2 == 0 or mcBit2+mcBit3 == 0 or mcBit3+mcBit1 == 0 ) hSEvt_Mt_KKX->Fill(mJJB);
             //else if ( mcBit1 ==0 or mcBit2 == 0 or mcBit3 == 0 ) hSEvt_Mt_KXY->Fill(mJJB);
             else
             { 
               hSEvt_Mt_XYZ->Fill(mJJB,SFbPUw1); 
               nhSEvt_Mt_XYZ->Fill(mJJB_new,SFbPUw1); 
               rhSEvt_Mt_XYZ->Fill(mJJB/mJJ,SFbPUw1); 
               rhSEvt_Mw_vsR->Fill(mJJB/mJJ,mJJ,SFbPUw1); 
               rhSEvt_Mt_vsR->Fill(mJJB/mJJ,mJJB,SFbPUw1); 
               hSEvt_Mt_XYZ->Fill(mJJB,SFbPUw1);
               hSEvt_eta_XYZ->Fill(etaJJB,SFbPUw1);
               hSEvt_pt_XYZ->Fill(ptJJB,SFbPUw1);
               hSEvt_phi_XYZ->Fill(phiJJB,SFbPUw1);
               hSEvt_dR_XYZ->Fill(dRJJB,SFbPUw1);
               hSEvt_dPhi_XYZ->Fill(dPhiJJB,SFbPUw1);
               hSEvt_dEta_XYZ->Fill(dEtaJJB,SFbPUw1);
             }
             //if(abs(jet3.eta())>2.4) cout << "iEvent = " << iEvent1 << ": eta & pT = " << jet3.eta() << " " << jet3.pt() << " mJJB = " << mJJB << " (SFb = " << SFb << ": up = " << SFb_up << " dn = " << SFb_dn << ")" << endl;
             hSEvt_Mt_SFb->Fill(SFb);
        }
      }
      //break;
    }

    //cout << "iEvent = " << iEvent1 << ": nLeadJets In SameEvt = " << nLeadJetsInSameEvt << endl;
    for ( int l1=0; l1<13; ++l1 ) 
    {
      if ( lepton1.eta() >= muEtaMin[l1] and lepton1.eta() < muEtaMax[l1] ) 
      {
        muScale1 = muidEff[l1]*muisoEff[l1]*mutrigEff[l1];

        hnvertex->Fill(event1.vertex,event1.PUweight*muScale1);

        hSEvt_Mu_pt->Fill(lepton1.pt(),event1.PUweight*muScale1);
        hSEvt_Mu_eta->Fill(lepton1.eta(),event1.PUweight*muScale1);
      }
    }    

    //hSEvt_nLeadJets->Fill(nLeadJetsInSameEvt,event1.PUweight); nLeadJetsInSameEvt = 0;
    hSEvt_nJets->Fill(nJets1,event1.PUweight); hSEvt_nBJets->Fill(nBjets1,event1.PUweight);
    hSEvt_nLightJets->Fill(nLightJets1,event1.PUweight);

    // Make Bi event jet combinations
    for ( int iEvent2 = iEvent1+1; ;++iEvent2 )
    {
      if ( iEvent2 == nEvent ) iEvent2 = 0;
      tree2->GetEntry(iEvent2);

      const LorentzVector& lepton2 = *event2.lepton;

      if( cutQ != 0 && event2.charge != cutQ ) continue;

      // ** DONOT use "continue" HERE because it causes event-mixing happen again (i.e. artifitial effect)
      //if ( event2.jets->size() < cut_minNJet ) continue;
      //if ( event2.jets->at(0).pt() < cut_minLeadJetPt or event2.jets->at(1).pt() < cut_minLeadJetPt ) continue;
      //cout << "-> " << iEvent2 << ": " << event2.jets->at(0).pt() << " " << event2.jets->at(1).pt() << endl;
      int nJets2 = 0, nBjets2 = 0, nLightJets2 = 0;
      for ( int j2=0, nj=event2.jets->size(); j2<nj; ++j2 )
      {
        const LorentzVector jet2 = event2.jets->at(j2);
        const double bTag2 = event2.bTags->at(j2);
        //if ( event2.jets->at(0).pt() < cut_minLeadJetPt or event2.jets->at(1).pt() < cut_minLeadJetPt ) continue;
        ////if ( jet2.pt() > cut_minJetPt or abs(jet2.eta()) < cut_maxJetEta ) ++nJets2; //a bug found at Jan 28, 2014
        if ( jet2.pt() > cut_minJetPt and abs(jet2.eta()) < cut_maxJetEta ) 
        {
          ++nJets2;
          if ( bTag2 > cut_bTag ) ++nBjets2;
          else ++nLightJets2;
        }
      }
      //if ( nJets2 < cut_minNJet || nBjets2 < cut_minNBjet || nLightJets2 < cut_minNLightJet ) continue;
      if ( nJets2 < cut_minNJet || nBjets2 < cut_minNBjet ) continue;
      //if ( nLightJets2 < 2 || nBjets2 < 2 ) continue;
      //      if ( nLightJets2 != 2 || nBjets2 != 2 ) continue;
      //      if ( abs(lepton2.eta()) >= cut_maxMuonEta ) continue;
      //      if ( abs(lepton2.eta()) < cut_maxMuonEta ) continue;
      else { hBEvt_nJets->Fill(nJets2,event2.PUweight); hBEvt_nBJets->Fill(nBjets2,event2.PUweight); }

      //int nBjet2 = 0;
      //for ( int j2=0, nj=event2.bTags->size(); j2<nj; ++j2 )
      //{
      //  if ( event2.bTags->at(j2) > cut_bTag ) ++nBjet2;
      //}
      //if ( nBjet2 < cut_minNBjet ) continue;

      //cout << "iEvent2 = " << iEvent2 << ": " << nJets2 << " " << nBjets2 << endl;
      break;

    }
    const LorentzVector& lepton2 = *event2.lepton;
    double dPhiLepton = lepton1.phi() - lepton2.phi();//ROOT::Math::VectorUtil::DeltaPhi(lepton1,lepton2); 
    ROOT::Math::RotationZ leptonRotation(dPhiLepton);
    for ( int j1=0, nj1=event1.jets->size(); j1<nj1; ++j1 )
    {
      const LorentzVector jet1 = event1.jets->at(j1);
      const double bTag1 = event1.bTags->at(j1);
      const int mcBit1 = event1.jetMCBits->at(j1);
      double jet1Eta = jet1.eta();

      //if ( DeltaR(jet1, lepton1) < 0.3 ) continue;
      if ( jet1.pt() < cut_minJetPt or abs(jet1.eta()) > cut_maxJetEta ) continue;
      //if ( jet1.pt() < cut_minLeadJetPt or abs(jet1.eta()) > cut_maxJetEta ) continue;
      if ( bTag1 > cut_bTag ) continue;
      //if ( mcBit1&3 ) continue;
      //if ( mcBit1 == 1 or mcBit1 == 2 ) continue;
      //if (jet1.pt() > cut_minLeadJetPt) ++nLeadJetsInBiEvt;

      for ( int j2=0, nj2=event2.jets->size(); j2<nj2; ++j2 )
      {
        //const LorentzVector jet2 = leptonRotation*(event2.jets->at(j2));
        const LorentzVector jet2 =event2.jets->at(j2);
        const double bTag2 = event2.bTags->at(j2);
        const int mcBit2 = event2.jetMCBits->at(j2);
        double jet2Eta = jet2.eta();

        //if ( DeltaR(jet2, lepton1) < 0.3 ) continue;
        if ( jet2.pt() < cut_minJetPt or abs(jet2.eta()) > cut_maxJetEta ) continue;
        //if ( jet2.pt() < cut_minLeadJetPt or abs(jet2.eta()) > cut_maxJetEta ) continue;
        //if ( j1 == j2 ) continue;
        if ( bTag2 > cut_bTag ) continue;
        //if ( mcBit2&3 ) continue;
        //if ( mcBit2 == 1 or mcBit2 == 2 ) continue;

        ++nBiEvent;
        if ( DeltaR(jet1, jet2) < cut_drBEvt ) continue;
        ++nPassedBiEvent;

        LorentzVector jj = jet1 + jet2; //event1.jets->at(j1)+event2.jets->at(j2);
        const double mJJ = jj.mass(); //min(mJJMax-1e-3, jj.mass());
        const double ptJJ = jj.pt();
        const double etaJJ = jj.eta();
        const double phiJJ = jj.phi();
        double dRJJ = ROOT::Math::VectorUtil::DeltaR(jet1,jet2);
        double dPhiJJ = ROOT::Math::VectorUtil::DeltaPhi(jet1,jet2);
        double dEtaJJ = jet1Eta-jet2Eta;
        double PUweight1 = event1.PUweight;
        double PUweight2 = event2.PUweight;
        double PUweight3 = PUweight1*PUweight2; //PUweight3 = 1;
        double Ptweight1 = event1.Ptweight;
        double Ptweight2 = event2.Ptweight;
        double Ptweight3 = Ptweight1*Ptweight2;

        for ( int l1=0; l1<13; ++l1 )
        {
          if ( lepton1.eta() >= muEtaMin[l1] and lepton1.eta() < muEtaMax[l1] )
          {
            muScale1 = muidEff[l1]*muisoEff[l1]*mutrigEff[l1];
          }
        }
        for ( int l2=0; l2<13; ++l2 )
        {
          if ( lepton2.eta() >= muEtaMin[l2] and lepton2.eta() < muEtaMax[l2] )
          {
            muScale2 = muidEff[l2]*muisoEff[l2]*mutrigEff[l2];
          }
        }

        double weight2 = PUweight3*muScale1*muScale2;

        hBEvt_Mw->Fill(mJJ,weight2);
        hBEvt_pt_Mw->Fill(ptJJ,weight2);
        hBEvt_eta_Mw->Fill(etaJJ,weight2);
        hBEvt_phi_Mw->Fill(phiJJ,weight2);
        hBEvt_dR_Mw->Fill(dRJJ,weight2);
        hBEvt_dPhi_Mw->Fill(dPhiJJ,weight2);
        hBEvt_dEta_Mw->Fill(dEtaJJ,weight2);
        if(input.Contains("Run2012") && jet1.pt() > cut_minLeadJetPt && jet2.pt() > cut_minLeadJetPt) 
        { 
          hBEvt_Mw_JJ->Fill(mJJ,weight2); 
          hBEvt_pt_JJ->Fill(ptJJ,weight2); 
          hBEvt_eta_JJ->Fill(etaJJ,weight2); 
          hBEvt_phi_JJ->Fill(phiJJ,weight2); 
          hBEvt_dR_JJ->Fill(dRJJ,weight2); 
          hBEvt_dPhi_JJ->Fill(dPhiJJ,weight2);
          hBEvt_dEta_JJ->Fill(dEtaJJ,weight2);
        }

        if ( mcBit1&4 and mcBit2&4 ) 
        { 
          hBEvt_Mw_JJ->Fill(mJJ,weight2); 
          hBEvt_pt_JJ->Fill(ptJJ,weight2); 
          hBEvt_eta_JJ->Fill(etaJJ,weight2); 
          hBEvt_phi_JJ->Fill(phiJJ,weight2);
          hBEvt_dR_JJ->Fill(dRJJ,weight2);
          hBEvt_dPhi_JJ->Fill(dPhiJJ,weight2);
          hBEvt_dEta_JJ->Fill(dEtaJJ,weight2);
        }
        else if ( (mcBit1&2 and mcBit2&4) or (mcBit1&4 and mcBit2&2) ) 
        {
          hBEvt_Mw_HBJ->Fill(mJJ,weight2);
          hBEvt_Mw_HBJ->Fill(mJJ,weight2);
          hBEvt_pt_HBJ->Fill(ptJJ,weight2);
          hBEvt_eta_HBJ->Fill(etaJJ,weight2);
          hBEvt_phi_HBJ->Fill(phiJJ,weight2);
          hBEvt_dR_HBJ->Fill(dRJJ,weight2);
          hBEvt_dPhi_HBJ->Fill(dPhiJJ,weight2);
          hBEvt_dEta_HBJ->Fill(dEtaJJ,weight2);
        }
        else if ( (mcBit1&1 and mcBit2&4) or (mcBit1&4 and mcBit2&1) ) 
        {
          hBEvt_Mw_LBJ->Fill(mJJ,weight2);
          hBEvt_Mw_LBJ->Fill(mJJ,weight2);
          hBEvt_pt_LBJ->Fill(ptJJ,weight2);
          hBEvt_eta_LBJ->Fill(etaJJ,weight2);
          hBEvt_phi_LBJ->Fill(phiJJ,weight2);
          hBEvt_dR_LBJ->Fill(dRJJ,weight2);
          hBEvt_dPhi_LBJ->Fill(dPhiJJ,weight2);
          hBEvt_dEta_LBJ->Fill(dEtaJJ,weight2);
        }
        else if ( mcBit1&3 or mcBit2&3 ) 
        {
          hBEvt_Mw_BX->Fill(mJJ,weight2);
          hBEvt_eta_BX->Fill(etaJJ,weight2);
          hBEvt_pt_BX->Fill(ptJJ,weight2);
          hBEvt_phi_BX->Fill(phiJJ,weight2);
          hBEvt_dR_BX->Fill(dRJJ,weight2);
          hBEvt_dPhi_BX->Fill(dPhiJJ,weight2);
          hBEvt_dEta_BX->Fill(dEtaJJ,weight2);
        }
        else if ( mcBit1 == 0 and mcBit2 == 0 ) 
        {
          hBEvt_Mw_KK->Fill(mJJ,weight2);
          hBEvt_Mw_KK->Fill(mJJ,weight2);
          hBEvt_pt_KK->Fill(ptJJ,weight2);
          hBEvt_eta_KK->Fill(etaJJ,weight2);
          hBEvt_phi_KK->Fill(phiJJ,weight2);
          hBEvt_dR_KK->Fill(dRJJ,weight2);
          hBEvt_dPhi_KK->Fill(dPhiJJ,weight2);
          hBEvt_dEta_KK->Fill(dEtaJJ,weight2);
        }
        else 
        {
          hBEvt_Mw_JK->Fill(mJJ,weight2);
          hBEvt_Mw_JK->Fill(mJJ,weight2);
          hBEvt_pt_JK->Fill(ptJJ,weight2);
          hBEvt_eta_JK->Fill(etaJJ,weight2);
          hBEvt_phi_JK->Fill(phiJJ,weight2);
          hBEvt_dR_JK->Fill(dRJJ,weight2);
          hBEvt_dPhi_JK->Fill(dPhiJJ,weight2);
          hBEvt_dEta_JK->Fill(dEtaJJ,weight2);
        }	
        for ( int j3=0; j3<nj2; ++j3 )
        {
          const LorentzVector jet3 = event2.jets->at(j3);
          const double bTag3 = event2.bTags->at(j3);
          const int mcBit3 = event2.jetMCBits->at(j3);
          double jet3Eta = jet3.eta();

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
          const double ptJJB = jjb.pt();
          const double etaJJB = jjb.eta();
          const double phiJJB = jjb.phi();
          double dRJJB = ROOT::Math::VectorUtil::DeltaR(jj,jet3);
          double dPhiJJB = ROOT::Math::VectorUtil::DeltaPhi(jj,jet3);
          double dEtaJJB = etaJJ-jet3Eta;
          const double mJJB_new = (mJJB/mJJ)*mW_world;
          double SFbPUw3 = SFb2*weight2; if(input.Contains("Run2012")) SFb2 = SFbPUw3 = 1;
          hBEvt_Mt->Fill(mJJB,SFbPUw3); 
          nhBEvt_Mt->Fill(mJJB_new,SFbPUw3); 
          rhBEvt_Mt->Fill(mJJB/mJJ,SFbPUw3);
          hBEvt_pt_Mt->Fill(ptJJB,SFbPUw3);
          hBEvt_eta_Mt->Fill(etaJJB,SFbPUw3);
          hBEvt_phi_Mt->Fill(phiJJB,SFbPUw3);
          hBEvt_dR_Mt->Fill(dRJJB,SFbPUw3);
          hBEvt_dPhi_Mt->Fill(dPhiJJB,SFbPUw3);
          hBEvt_dEta_Mt->Fill(dEtaJJB,SFbPUw3);
          if(input.Contains("Run2012") && jet1.pt() > cut_minLeadJetPt && jet2.pt() > cut_minLeadJetPt) 
          { 
            hBEvt_Mt_JJHB->Fill(mJJB,SFbPUw3); 
            hBEvt_pt_JJHB->Fill(ptJJB,SFbPUw3); 
            hBEvt_eta_JJHB->Fill(etaJJB,SFbPUw3); 
            hBEvt_phi_JJHB->Fill(phiJJB,SFbPUw3);
            hBEvt_dR_JJHB->Fill(dRJJB,SFbPUw3);
            hBEvt_dPhi_JJHB->Fill(dPhiJJB,SFbPUw3);
            hBEvt_dEta_JJHB->Fill(dEtaJJB,SFbPUw3);
          }

          if ( (mcBit1&4 and mcBit2&4 and mcBit3&3) ) hBEvt_Mt_NEW->Fill(mJJB,SFbPUw3);

          if ( (mcBit1&2 and mcBit2&4 and mcBit3&4) or
              (mcBit1&4 and mcBit2&2 and mcBit3&4) or
              (mcBit1&4 and mcBit2&4 and mcBit3&2) ) { hBEvt_Mt_A->Fill(mJJB,SFbPUw3); nhBEvt_Mt_A->Fill(mJJB_new,SFbPUw3); rhBEvt_Mt_A->Fill(mJJB/mJJ,SFbPUw3); }
          else { hBEvt_Mt_B->Fill(mJJB,SFbPUw3); nhBEvt_Mt_B->Fill(mJJB_new,SFbPUw3); rhBEvt_Mt_B->Fill(mJJB/mJJ,SFbPUw3); }

          if ( (mcBit1&4 and mcBit2&4 and mcBit3&2) ) 
          { 
            hBEvt_Mt_JJHB->Fill(mJJB,SFbPUw3); 
            nhBEvt_Mt_JJHB->Fill(mJJB_new,SFbPUw3); 
            rhBEvt_Mt_JJHB->Fill(mJJB/mJJ,SFbPUw3); 
            hBEvt_pt_JJHB->Fill(ptJJB,SFbPUw3); 
            hBEvt_eta_JJHB->Fill(etaJJB,SFbPUw3); 
            hBEvt_phi_JJHB->Fill(phiJJB,SFbPUw3);
            hBEvt_dR_JJHB->Fill(dRJJB,SFbPUw3);
            hBEvt_dPhi_JJHB->Fill(dPhiJJB,SFbPUw3);
            hBEvt_dEta_JJHB->Fill(dEtaJJB,SFbPUw3);
          }
          //          if ( (mcBit1&4 and mcBit2&2 and mcBit3&4) or (mcBit1&4 and mcBit2&4 and mcBit3&2) ) hBEvt_Mt_JJHB->Fill(mJJB,SFbPUw3);
          //          else if ( (mcBit1&4 and mcBit2&1 and mcBit3&4) or (mcBit1&4 and mcBit2&4 and mcBit3&1) ) hBEvt_Mt_JJLB->Fill(mJJB,SFbPUw3);
          //          else if ( (mcBit1&4 and mcBit2 == 0 and mcBit3&4) or (mcBit1&4 and mcBit2&4 and mcBit3 == 0) ) hBEvt_Mt_HBLBJ->Fill(mJJB,SFbPUw3);
          //else if ( (mcBit1&2 and mcBit2&4 and mcBit3&4) ) hBEvt_Mt_NEW->Fill(mJJB,SFbPUw3);
          ////else if ( (mcBit1&2 and mcBit2&4 and mcBit3&4) or (mcBit1&1 and mcBit2&4 and mcBit3&4) ) hBEvt_Mt_NEW->Fill(mJJB,SFbPUw3);

          //          if ( (mcBit1&4 and mcBit2&4 and mcBit3&2) ) hBEvt_Mt_JJHB->Fill(mJJB,SFbPUw3);
          //          else if ( (mcBit1&4 and mcBit2&4 and mcBit3&1) ) hBEvt_Mt_JJLB->Fill(mJJB,SFbPUw3);
          //          else if ( (mcBit1&4 and mcBit2&2 and mcBit3&4) or (mcBit1&4 and mcBit2&1 and mcBit3&4) ) hBEvt_Mt_HBLBJ->Fill(mJJB,SFbPUw3);

          //          if ( (mcBit1&4 and mcBit2&2 and mcBit3&4) or (mcBit1&4 and mcBit2&4 and mcBit3&2) ) hBEvt_Mt_JJHB->Fill(mJJB,SFbPUw3);
          //          else if ( (mcBit1&4 and mcBit2&1 and mcBit3&4) or (mcBit1&4 and mcBit2&4 and mcBit3&1) ) hBEvt_Mt_JJLB->Fill(mJJB,SFbPUw3);

          //          if ( (mcBit1&2 and mcBit2&4 and mcBit3&4) or 
          //               (mcBit1&4 and mcBit2&2 and mcBit3&4) or 
          //               (mcBit1&4 and mcBit2&4 and mcBit3&2) ) hBEvt_Mt_JJHB->Fill(mJJB,SFbPUw3);
          else if ( (mcBit1&1 and mcBit2&4 and mcBit3&4) or
              (mcBit1&4 and mcBit2&1 and mcBit3&4) or 
              (mcBit1&4 and mcBit2&4 and mcBit3&1) )
          {
            hBEvt_Mt_JJLB->Fill(mJJB,SFbPUw3);
            hBEvt_pt_JJLB->Fill(mJJB,SFbPUw3);
            hBEvt_eta_JJLB->Fill(etaJJB,SFbPUw3);
            hBEvt_phi_JJLB->Fill(phiJJB,SFbPUw3);
            hBEvt_dR_JJLB->Fill(dRJJB,SFbPUw3);
            hBEvt_dPhi_JJLB->Fill(dPhiJJB,SFbPUw3);
            hBEvt_dEta_JJLB->Fill(dEtaJJB,SFbPUw3);
          }
          else if ( (mcBit1&1 and mcBit2&2 and mcBit3&4) or
              (mcBit1&2 and mcBit2&1 and mcBit3&4) or
              (mcBit1&4 and mcBit2&1 and mcBit3&2) or
              (mcBit1&2 and mcBit2&4 and mcBit3&1) or
              (mcBit1&4 and mcBit2&1 and mcBit3&2) or
              (mcBit1&4 and mcBit2&2 and mcBit3&1) ) 
          {
            hBEvt_Mt_HBLBJ->Fill(mJJB,SFbPUw3); //added on Jan 23, 2014
            hBEvt_pt_HBLBJ->Fill(mJJB,SFbPUw3);
            hBEvt_eta_HBLBJ->Fill(etaJJB,SFbPUw3);
            hBEvt_phi_HBLBJ->Fill(phiJJB,SFbPUw3);
            hBEvt_dR_HBLBJ->Fill(dRJJB,SFbPUw3);
            hBEvt_dPhi_HBLBJ->Fill(dPhiJJB,SFbPUw3);
            hBEvt_dEta_HBLBJ->Fill(dEtaJJB,SFbPUw3);
          }

          //else if ( mcBit1+mcBit2 == 0 or mcBit2+mcBit3 == 0 or mcBit3+mcBit1 == 0 ) hBEvt_Mt_KKX->Fill(mJJB);
          //`else if ( mcBit1 ==0 or mcBit2 == 0 or mcBit3 == 0 ) hBEvt_Mt_KXY->Fill(mJJB);
          else 
          { 
            hBEvt_Mt_XYZ->Fill(mJJB,SFbPUw3); 
            nhBEvt_Mt_XYZ->Fill(mJJB_new,SFbPUw3); 
            rhBEvt_Mt_XYZ->Fill(mJJB/mJJ,SFbPUw3); 
            rhBEvt_Mw_vsR->Fill(mJJB/mJJ,mJJ,SFbPUw3); 
            rhBEvt_Mt_vsR->Fill(mJJB/mJJ,mJJB,SFbPUw3); 
            hBEvt_Mt_XYZ->Fill(mJJB,SFbPUw3);
            hBEvt_eta_XYZ->Fill(etaJJB,SFbPUw3);
            hBEvt_pt_XYZ->Fill(ptJJB,SFbPUw3);
            hBEvt_phi_XYZ->Fill(phiJJB,SFbPUw3);
            hBEvt_dR_XYZ->Fill(dRJJB,SFbPUw3);
            hBEvt_dPhi_XYZ->Fill(dPhiJJB,SFbPUw3);
            hBEvt_dEta_XYZ->Fill(dEtaJJB,SFbPUw3);
          }
          hBEvt_Mt_SFb->Fill(SFb2);
        }
      }
      //break;

    }
    //hBEvt_nLeadJets->Fill(nLeadJetsInBiEvt,event2.PUweight); nLeadJetsInBiEvt = 0;

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
