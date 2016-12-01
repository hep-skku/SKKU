//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov 28 22:53:09 2016 by ROOT version 6.06/01
// from TTree tree/tree
// found on file: /xrootd/store/user/jhgoh/MuonMisID/20161125_1/JetHT/JetHT_Run2016B-23Sep2016-v3/161125_042627/0000/ntuple_1.root
//////////////////////////////////////////////////////////

#ifndef MisID_h
#define MisID_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>
#include <string>
#include <map>
#include <TH1D.h>
#include <TH2D.h>

class MisID : public TSelector {
public :
  TTreeReader     fReader;  //!the tree reader
  TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

  // Readers to access the data (delete the ones you do not need).
  TTreeReaderValue<UChar_t> run = {fReader, "run"};
  TTreeReaderValue<UChar_t> lumi = {fReader, "lumi"};
  TTreeReaderValue<ULong64_t> event = {fReader, "event"};

  TTreeReaderValue<UChar_t> nPV = {fReader, "nPV"};
  TTreeReaderValue<UChar_t> nSV = {fReader, "nSV"};
  TTreeReaderValue<UChar_t> nGen = {fReader, "nGen"};

  std::vector<double> ptbins_, aetabins_;
  std::map<std::string, TTreeReaderValue<float>*> floatVars_;
  std::map<std::string, TTreeReaderValue<float>> floatVars2_;
  std::map<std::string, TTreeReaderValue<std::vector<float>>*> vfloatVars_;
  std::map<std::string, TTreeReaderValue<std::vector<int>>*> vintVars_;
  std::map<std::string, TTreeReaderValue<std::vector<bool>>*> vboolVars_;

  MisID(TTree * /*tree*/ =0);
  virtual ~MisID() {
    for ( auto key = floatVars_.begin(); key != floatVars_.end(); ++key ) delete key->second;
    for ( auto key = vfloatVars_.begin(); key != vfloatVars_.end(); ++key ) delete key->second;
    for ( auto key = vintVars_.begin(); key != vintVars_.end(); ++key ) delete key->second;
    for ( auto key = vboolVars_.begin(); key != vboolVars_.end(); ++key ) delete key->second;
  }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }

  TString mode_;
  std::map<std::string, TH1*> h_;

  void book(const char* name, const char* title,
            const std::vector<double>& xbins, const std::vector<double>& ybins = std::vector<double>())
  {
    if ( ybins.empty() ) h_[name] = new TH1D(name, title, xbins.size()-1, &xbins[0]);
    else h_[name] = new TH2D(name, title, xbins.size()-1, &xbins[0], ybins.size()-1, &ybins[0]);
  }

  void book(const char* name, const char* title,
            const int nbinx, const double xmin, const double xmax,
            const int nbiny=0, const double ymin=0, const double ymax=0)
  {
    if ( nbiny == 0 ) h_[name] = new TH1D(name, title, nbinx, xmin, xmax);
    else h_[name] = new TH2D(name, title, nbinx, xmin, xmax, nbiny, ymin, ymax);
  }

  void fill(const char* name, const double x, const double y=1)
  {
    auto h = h_.find(name);
    if ( h == h_.end() ) return;
    h->second->Fill(x, y);
  }

  virtual void    SlaveTerminate()
  {
    for ( auto x = h_.begin(); x != h_.end(); ++x ) {
      fOutput->Add(x->second);
    }
  }

  virtual void    Terminate()
  {
    cout << "MisID::Terminate() " << mode_ << endl;
    TFile fout(Form("hist_%s.root", mode_.Data()), "recreate");
    fout.cd();
    for ( int i=0, n=fOutput->GetEntries(); i<n; ++i ) {
      auto obj = fOutput->At(i);
      TH1* h1 = dynamic_cast<TH2D*>(obj);
      if ( !h1 ) h1 = dynamic_cast<TH1D*>(obj);
      if ( !h1 ) continue;

      TDirectory* dir = &fout;
      auto hPath = (mode_+"/"+h1->GetName()).Tokenize("/");
      TString hName = hPath->At(hPath->GetEntries()-1)->GetName();
      for ( int i=0; i<hPath->GetEntries()-1; ++i ) {
        TString dirName = hPath->At(i)->GetName();
        TDirectory* tdir = dir->GetDirectory(dirName);
        if ( tdir ) dir = tdir;
        else dir = dir->mkdir(dirName);
      }

      dir->cd();
      h1->SetName(hName);
      h1->Write();

/*      cout << "==== " << h1->GetName() << " ====\n";
      cout << "* nEntries = " << h1->GetEntries() << endl;
      cout << "* Integral = " << h1->Integral() << endl;
      cout << "* Mean     = " << h1->GetMean() << endl;
      cout << "* RMS      = " << h1->GetRMS() << endl;*/
    }
    fout.Close();
  }

  ClassDef(MisID,0);

};

#endif

#ifdef MisID_cxx
void MisID::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t MisID::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef MisID_cxx
