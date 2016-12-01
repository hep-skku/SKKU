#define MisID_cxx

#include "MisID.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>

MisID::MisID(TTree* tree)
{
  ptbins_ = {4, 5, 6, 8, 10, 15, 20, 30, 200};
  aetabins_ = {0, 0.9, 1.2, 1.6, 2.1, 2.4};

  std::vector<const char*> vboolVarNames = {
    "mu_isTight", "mu_isMedium", "mu_isLoose",
    "mu_isSoft", "mu_isHighPt",
    "mu_isGLB", "mu_isTRK", "mu_isSTA", "mu_isRPC",
    "mu_isGLBPT",
    "mu_isTMLastLoose", "mu_isTMLastTight",
    "mu_isTM2DLoose", "mu_isTM2DTight",
    "mu_isOneLoose", "mu_isOneTight",
    "mu_isLastLowPtLoose", "mu_isLastLowPtTight",
    "mu_isGMTkChi2Compat", "mu_isGMStaChi2Compat", "mu_isGMTkKinkTight", 
    "mu_isTMLastAngLoose", "mu_isTMLastAngTight", "mu_isTMOneAngLoose", "mu_isTMOneAngTight",
    "mu_isRPCMuLoose", "mu_RPCLoose", "mu_TStLoose",
  };

  std::vector<const char*> vintVarNames = {
    "trk_pdgId", "mu_q"
  };

  std::vector<const char*> floatVarNames = {
    "vtx_mass", "vtx_pt", "vtx_lxy", "vtx_vz",
    "gen_dR"
  };

  std::vector<const char*> vfloatVarNames = {
    "trk_pt", "trk_eta", "trk_phi",
    "mu_dR", "mu_pt",
  };

  for ( auto name : vboolVarNames ) {
    vboolVars_[name] = new TTreeReaderValue<std::vector<bool>>(fReader, name);
  }
  for ( auto name : vintVarNames ) {
    vintVars_[name] = new TTreeReaderValue<std::vector<int>>(fReader, name);
  }
  for ( auto name : floatVarNames ) {
    floatVars_[name] = new TTreeReaderValue<float>(fReader, name);
  }
  for ( auto name : vfloatVarNames ) {
    vfloatVars_[name] = new TTreeReaderValue<std::vector<float>>(fReader, name);
  }

}

void MisID::Begin(TTree * tree)
{
  TString option = GetOption();
  auto options = option.Tokenize(",");
  for ( int i=0; i<options->GetEntries(); ++i ) {
    TString opt = options->At(i)->GetName();
    if ( opt.Length() > 5 and opt.BeginsWith("mode=") ) mode_ = opt(5,opt.Length());
    cout << "Option" << i << " = \"" << opt << "\"" << endl;
  }
  cout << "MODE=" << mode_ << endl;
}

void MisID::SlaveBegin(TTree * /*tree*/)
{
  TString option = GetOption();
  auto options = option.Tokenize(",");
  for ( int i=0; i<options->GetEntries(); ++i ) {
    TString opt = options->At(i)->GetName();
    if ( opt.Length() > 5 and opt.BeginsWith("mode=") ) mode_ = opt(5,opt.Length());
  }

  const int nbinMass = 50;
  const double minMass = mode_ == "ks" ? 0.45 : mode_ == "phi" ? 1.00 : mode_ == "lamb" ? 1.10 : 3.00;
  const double maxMass = mode_ == "ks" ? 0.55 : mode_ == "phi" ? 1.04 : mode_ == "lamb" ? 1.13 : 3.18;

  // book histograms

  cout << "NVars=" << vboolVars_.size() << endl;
  for ( auto key = vboolVars_.begin(); key != vboolVars_.end(); ++key ) {
    const char* idName = key->first.c_str();
    for ( int leg=1; leg<=2; ++leg ) {
      book(Form("%s_leg%d/pt/hFrame", idName, leg), "pt;p_{T} (GeV)", ptbins_);
      for ( unsigned int i=0; i<ptbins_.size(); ++i ) {
        book(Form("%s_leg%d/pt/bin%d/hPass", idName, leg, i+1), "pass;Mass (GeV)", nbinMass, minMass, maxMass);
        book(Form("%s_leg%d/pt/bin%d/hFail", idName, leg, i+1), "fail;Mass (GeV)", nbinMass, minMass, maxMass);
      }

      book(Form("%s_leg%d/abseta/hFrame", idName, leg), "abseta;|#eta|", aetabins_);
      for ( unsigned int i=0; i<aetabins_.size(); ++i ) {
        book(Form("%s_leg%d/abseta/bin%d/hPass", idName, leg, i+1), "pass;Mass (GeV)", nbinMass, minMass, maxMass);
        book(Form("%s_leg%d/abseta/bin%d/hFail", idName, leg, i+1), "fail;Mass (GeV)", nbinMass, minMass, maxMass);
      }
    }
  }
}

Bool_t MisID::Process(Long64_t entry)
{
  fReader.SetEntry(entry);

  const double mass = **floatVars_["vtx_mass"];
  const double lxy = **floatVars_["vtx_lxy"];
  if ( lxy > 4 ) return false;

  for ( auto key = vboolVars_.begin(); key != vboolVars_.end(); ++key ) {
    const char* idName = key->first.c_str();
    const double maxAEta = string(idName).find("RPC") != string::npos ? 2.1 : 2.4;
    const double minPt = 4;

    for ( int leg=1; leg<=2; ++leg ) {
      const bool idRes = ( (*vfloatVars_["mu_dR"])->at(leg-1) < 0.01 and 
                           (*vboolVars_[idName])->at(leg-1) );
      const double pt = (*vfloatVars_["trk_pt"])->at(leg-1);
      const double aeta = std::abs((*vfloatVars_["trk_eta"])->at(leg-1));
      // Apply common acceptance cuts
      if ( aeta > maxAEta or pt < minPt ) continue;

      auto hFramePt = (TH1D*)h_[Form("%s_leg%d/pt/hFrame", idName, leg)];
      auto hFrameAEta = (TH1D*)h_[Form("%s_leg%d/abseta/hFrame", idName, leg)];
      const int ptbin = hFramePt->GetXaxis()->FindBin(pt);
      const int aetabin = hFrameAEta->GetXaxis()->FindBin(aeta);

      if ( idRes ) {
        fill(Form("%s_leg%d/pt/bin%d/hPass", idName, leg, ptbin), mass);
        fill(Form("%s_leg%d/abseta/bin%d/hPass", idName, leg, aetabin), mass);
      }
      else {
        fill(Form("%s_leg%d/pt/bin%d/hFail", idName, leg, ptbin), mass);
        fill(Form("%s_leg%d/abseta/bin%d/hFail", idName, leg, aetabin), mass);
      }
    }
  }

  return kTRUE;
}

