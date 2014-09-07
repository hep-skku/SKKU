//root -l hrandom.C\(\"MSDecays_central_non\",\"MSDecays_central\"\)
//root -l hrandom.C\(\"MSDecays_central_non\",\"Ptweight\"\)
//root -l hrandom.C\(\"MSDecays_central_non\",\"Up\"\)
//root -l -b -q hrandom.C\(\"MSDecays_central_non\",\"Dn\"\)
//root -l -b -q hrandom.C\(\"MSDecays_central_non_FlavorPureBottom\",\"Up\"\); root -l -b -q hrandom.C\(\"MSDecays_central_non_FlavorPureBottom\",\"Dn\"\)
//root -l -b -q hrandom.C\(\"MSDecays_central_non_FlavorPureGluon\",\"Up\"\); root -l -b -q hrandom.C\(\"MSDecays_central_non_FlavorPureGluon\",\"Dn\"\)
//root -l -b -q hrandom.C\(\"MSDecays_central_non_FlavorPureQuark\",\"Up\"\); root -l -b -q hrandom.C\(\"MSDecays_central_non_FlavorPureQuark\",\"Dn\"\)
//root -l -b -q hrandom.C\(\"MSDecays_central_non_FlavorPureCharm\",\"Up\"\); root -l -b -q hrandom.C\(\"MSDecays_central_non_FlavorPureCharm\",\"Dn\"\)

//root -l -b -q hrandom.C\(\"MSDecays_169\",\"\"\)
//root -l -b -q hrandom.C\(\"MSDecays_169_non\",\"\"\)
//root -l -b -q hrandom.C\(\"MSDecays_169_part1_non\",\"\"\)

//root -l -b -q hrandom.C\(\"MSDecays_central\"\)
//root -l -b -q hrandom.C\(\"MSDecays_central_non\"\)
//root -l hrandom.C\(\"MSDecays_central_non\",\"\"\,10000,\"30\"\)
//root -l -b -q hrandom.C\(\"ALL_processes\"\)
//hadd hist_ALL_processes_pt35_nj4_nb2_lq0.root hist_MSDecays_central_non_pt35_nj4_nb2_lq0.root hist_T_schannel_pt35_nj4_nb2_lq0.root hist_T_tchannel_v1_pt35_nj4_nb2_lq0.root hist_T_tWchannel_pt35_nj4_nb2_lq0.root hist_Tbar_schannel_pt35_nj4_nb2_lq0.root hist_Tbar_tchannel_pt35_nj4_nb2_lq0.root hist_Tbar_tWchannel_pt35_nj4_nb2_lq0.root hist_WJetsToLNu_v2_pt35_nj4_nb2_lq0.root hist_DYJets_M50_pt35_nj4_nb2_lq0.root hist_WW_pt35_nj4_nb2_lq0.root hist_WZ_pt35_nj4_nb2_lq0.root hist_ZZ_pt35_nj4_nb2_lq0.root

//root -l -b -q hrandom.C\(\"MSDecays_central_non_CorrelationGroupIntercalibration\",\"Up\"\); root -l -b -q hrandom.C\(\"MSDecays_central_non_CorrelationGroupIntercalibration\",\"Dn\"\)
//root -l -b -q hrandom.C\(\"MSDecays_central_non_CorrelationGroupUncorrelated\",\"Up\"\); root -l -b -q hrandom.C\(\"MSDecays_central_non_CorrelationGroupUncorrelated\",\"Dn\"\)
//root -l -b -q hrandom.C\(\"MSDecays_central_non_CorrelationGroupMPFInSitu\",\"Up\"\); root -l -b -q hrandom.C\(\"MSDecays_central_non_CorrelationGroupMPFInSitu\",\"Dn\"\)

//root -l hrandom.C\(\"MSDecays_central_non_PDF\",\"PDFWeight0\"\);

void hrandom(TString mc,/* TString jes = "",*/ int numToys = 1000, TString pt="35", TString nb="2")
{

  gROOT->ProcessLine(".x ~/rootlogon.C");
  tdrStyle->cd();

  gStyle->SetStatY(0.87);
  gStyle->SetStatX(0.93);
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.15);

  gStyle->SetOptStat(1111111);

  TH1::AddDirectory(kFALSE); //sets a global switch disabling the reference

  //TString pt = "35";
  //TString nb = "2";
  TString lq = "0";

  int ndata = 139859;

  //const char* fileName = "hpseudo/hnewPU_"+mc+"_pt"+pt+"_nj4_nb"+nb+"_lq"+lq;
  //const char* name = "hpas/hnewPU";

  //const char* fileName = "hpseudo/hpseudo"+jes+"_"+mc+"_pt"+pt+"_nj4_nb"+nb+"_lq"+lq;
  //const char* name = "hpas/hist"+jes+"_"+mc+"_pt"+pt+"_nj4_nb"+nb+"_lq"+lq;
  //const char* name = "hpas/hnew"+jes+"_"+mc+"_pt"+pt+"_nj4_nb"+nb+"_lq"+lq; //b-JES by flavor and JER at 8 TeV
  //const char* name = "hpas/hist"+jes+"_"+mc+"_pt"+pt+"_nj4_nb"+nb+"_lq"+lq; //b-JES by flavor and JER at 8 TeV
  const char* name = "hist/hist_"+mc+"_pt"+pt+"_nj4_nb"+nb+"_lq"+lq+"_nlj3";
  const char* fileName = "hpseudo/hpseudo_"+mc+"_pt"+pt+"_nj4_nb"+nb+"_lq"+lq+"_nlj3";

  //printf("%s.root\n",name); break;

  const int m = 1;
  const int n = numToys; //number of toys

  TFile *f[m];

  f[0] = new TFile(Form("%s.root",name ));
  //f[1] = new TFile(Form("%s_%s_pt%s_nj4_nb%s_lq%s.root",name,signal.Data(),pt.Data(),nb.Data(),lq.Data())); //for biW

  string sel;

  TH1F *hsameTop[m], *hbiTop[m], *rhsameTop[m], *rhbiTop[m], *hbiTop_JJHB[m], *rhbiTop_JJHB[m]; 
  TH1F *hMt[n],  *rhMt[n],  *hMt_JJHB[n], *rhMt_JJHB[n], *hMt_XYZ[n], *rhMt_XYZ[n];
  TH1F *hMt_[n], *rhMt_[n];
  //TH1F *hMt_JJHB_[n], *rhMt_JJHB_[n];

  for(int k=0; k<m; ++k) {
    
    sel = Form("SEvt");

    hsameTop[k] = (TH1F*)f[k]->Get(Form("%s/%s",sel.c_str(),"hMt")); hsameTop[k]->SetName(Form("hMtop%d",k));
    rhsameTop[k] = (TH1F*)f[k]->Get(Form("%s/%s",sel.c_str(),"rhMt")); rhsameTop[k]->SetName(Form("rhMtop%d",k));
    hsameTop[k]->GetXaxis()->SetTitle("m_{jjb} (GeV)"); rhsameTop[k]->GetXaxis()->SetTitle("m_{jjb}/m_{jj}");

    hsameTop_JJHB[k] = (TH1F*)f[k]->Get(Form("%s/%s",sel.c_str(),"hMt_JJHB")); hsameTop_JJHB[k]->SetName(Form("hMtop_JJHB%d",k));
    rhsameTop_JJHB[k] = (TH1F*)f[k]->Get(Form("%s/%s",sel.c_str(),"rhMt_JJHB")); rhsameTop_JJHB[k]->SetName(Form("rhMtop_JJHB%d",k));
    hsameTop_JJHB[k]->GetXaxis()->SetTitle("m_{jjb} (GeV)"); rhsameTop_JJHB[k]->GetXaxis()->SetTitle("m_{jjb}/m_{jj}");

    hsameTop_XYZ[k] = (TH1F*)f[k]->Get(Form("%s/%s",sel.c_str(),"hMt_XYZ")); hsameTop_XYZ[k]->SetName(Form("hMtop_XYZ%d",k));
    rhsameTop_XYZ[k] = (TH1F*)f[k]->Get(Form("%s/%s",sel.c_str(),"rhMt_XYZ")); rhsameTop_XYZ[k]->SetName(Form("rhMtop_XYZ%d",k));
    hsameTop_XYZ[k]->GetXaxis()->SetTitle("m_{jjb} (GeV)"); rhsameTop_XYZ[k]->GetXaxis()->SetTitle("m_{jjb}/m_{jj}");

    hsameTop[k]->SetMarkerStyle(20); rhsameTop[k]->SetMarkerStyle(20);

    sel = Form("BEvt");

    hbiTop[k] = (TH1F*)f[k]->Get(Form("%s/%s",sel.c_str(),"hMt")); hbiTop[k]->SetName(Form("hMtop%d",k));
    rhbiTop[k] = (TH1F*)f[k]->Get(Form("%s/%s",sel.c_str(),"rhMt")); rhbiTop[k]->SetName(Form("rhMtop%d",k));
    hbiTop[k]->GetXaxis()->SetTitle("m_{jjb} (GeV)"); rhbiTop[k]->GetXaxis()->SetTitle("m_{jjb}/m_{jj}");

    hbiTop_JJHB[k] = (TH1F*)f[k]->Get(Form("%s/%s",sel.c_str(),"hMt_JJHB")); hbiTop_JJHB[k]->SetName(Form("hMtop_JJHB%d",k));
    rhbiTop_JJHB[k] = (TH1F*)f[k]->Get(Form("%s/%s",sel.c_str(),"rhMt_JJHB")); rhbiTop_JJHB[k]->SetName(Form("rhMtop_JJHB%d",k));
    hbiTop_JJHB[k]->GetXaxis()->SetTitle("m_{jjb} (GeV)"); rhbiTop_JJHB[k]->GetXaxis()->SetTitle("m_{jjb}/m_{jj}");

    hbiTop[k]->SetMarkerStyle(24); rhbiTop[k]->SetMarkerStyle(24);

  }

  for(int i=0; i<n; ++i) {

    hMt[i]  = (TH1F*)hsameTop[0]->Clone(Form("hMt%d",i)); hMt[i]->Reset();
    rhMt[i]  = (TH1F*)rhsameTop[0]->Clone(Form("rhMt%d",i)); rhMt[i]->Reset();

    hMt_JJHB[i] = (TH1F*)hbiTop_JJHB[0]->Clone(Form("hMt_JJHB%d",i)); hMt_JJHB[i]->Reset();
    rhMt_JJHB[i] = (TH1F*)rhbiTop_JJHB[0]->Clone(Form("rhMt_JJHB%d",i)); rhMt_JJHB[i]->Reset();

    hMt_[i] = (TH1F*)hbiTop[0]->Clone(Form("hMt%d",i)); hMt_[i]->Reset();
    rhMt_[i] = (TH1F*)rhbiTop[0]->Clone(Form("rhMt%d",i)); rhMt_[i]->Reset();

    //hMt_JJHB_[i] = (TH1F*)hbiTop_JJHB[0]->Clone(Form("hMt_JJHB%d",i)); hMt_JJHB_[i]->Reset();
    //rhMt_JJHB_[i] = (TH1F*)rhbiTop_JJHB[0]->Clone(Form("rhMt_JJHB%d",i)); rhMt_JJHB_[i]->Reset();

    hMt_XYZ[i] = (TH1F*)hbiTop_JJHB[0]->Clone(Form("hMt_XYZ%d",i)); hMt_XYZ[i]->Reset();
    rhMt_XYZ[i] = (TH1F*)rhbiTop_JJHB[0]->Clone(Form("rhMt_XYZ%d",i)); rhMt_XYZ[i]->Reset();

  }

  /*
  //--too slow
  for(int i=0; i<n; ++i) {
    for(int j=0; j<ndata; ++j) {
      TRandom3 r(0);
      double sameMt = hsameTop[0]->GetRandom(),           sameR = rhsameTop[0]->GetRandom();
      double sameMt_JJHB = hsameTop_JJHB[0]->GetRandom(), sameR_JJHB = rhsameTop_JJHB[0]->GetRandom();
      double biMt = hbiTop[0]->GetRandom(),               biR  = rhbiTop[0]->GetRandom();
      double biMt_JJHB = hbiTop_JJHB[0]->GetRandom(),     biR_JJHB = rhbiTop_JJHB[0]->GetRandom();
      hMt[i]->Fill(sameMt); rhMt[i]->Fill(sameR); hMt_JJHB[i]->Fill(sameMt_JJHB); rhMt_JJHB[i]->Fill(sameR_JJHB);
      hMt_[i]->Fill(biMt); rhMt_[i]->Fill(biR); hMt_JJHB_[i]->Fill(biMt_JJHB); rhMt_JJHB_[i]->Fill(biR_JJHB);
    }
  }
  */

  //TH1F *hrnd[m]; hrnd[0] = (TH1F*)hsameTop[0]->Clone(Form("hrnd%d",0)); hrnd[0]->Reset();
  //hrnd[0]->FillRandom(hsameTop[0],ndata);

  for(int i=0; i<n; ++i) {
    hMt[i]->FillRandom(hsameTop[0],ndata); rhMt[i]->FillRandom(rhsameTop[0],ndata);
    hMt_JJHB[i]->FillRandom(hsameTop_JJHB[0],ndata); rhMt_JJHB[i]->FillRandom(rhsameTop_JJHB[0],ndata);
    hMt_[i]->FillRandom(hbiTop[0],ndata); rhMt_[i]->FillRandom(rhbiTop[0],ndata);
    //hMt_JJHB_[i]->FillRandom(hbiTop_JJHB[0],ndata); rhMt_JJHB_[i]->FillRandom(rhbiTop_JJHB[0],ndata);
    hMt_XYZ[i]->FillRandom(hsameTop_XYZ[0],ndata); rhMt_XYZ[i]->FillRandom(rhsameTop_XYZ[0],ndata);
  }

  TCanvas *ccv = new TCanvas("ccv","ccv",0,0,500,500);
  hsameTop[0]->Draw();

  TCanvas *ccv2 = new TCanvas("ccv2","ccv2",50,50,500,700);
  ccv2->Divide(1,2);
  ccv2->cd(1); hMt[0]->Draw();
  hsameTop[0]->DrawNormalized("esames",ndata);
  ccv2->cd(2); hMt_[0]->Draw();
  hbiTop[0]->DrawNormalized("esames",ndata);

  TCanvas *ccv3 = new TCanvas("ccv3","ccv3",100,100,600,600);
  ccv3->Divide(3,3);
  for(int i=0; i<9; ++i) {
    ccv3->cd(i+1); hMt[i]->Draw();
  }

  TCanvas *ccv4 = new TCanvas("ccv4","ccv4",150,150,600,600);
  ccv4->Divide(3,3);
  for(int i=0; i<9; ++i) {
    ccv4->cd(i+1); hMt_[i]->Draw();
  }

  TString fileout = Form("%s_%d.root",fileName,numToys);  // name of the output of the code
  TFile * outFile = new TFile(fileout, "RECREATE");
  TDirectory* dirSEvt = outFile->mkdir("SEvt");
  TDirectory* dirBEvt = outFile->mkdir("BEvt");
  dirSEvt->cd(); hsameTop[0]->Write(); rhsameTop[0]->Write(); hsameTop_JJHB[0]->Write(); rhsameTop_JJHB[0]->Write();
  //for(int i=0; i<n; ++i) { hMt[i]->Write(); rhMt[i]->Write(); }
  for(int i=0; i<n; ++i) { hMt[i]->Write(); rhMt[i]->Write(); hMt_JJHB[i]->Write(); rhMt_JJHB[i]->Write(); hMt_XYZ[i]->Write(); rhMt_XYZ[i]->Write(); }
  dirBEvt->cd(); hbiTop[0]->Write(); rhbiTop[0]->Write(); hbiTop_JJHB[0]->Write(); rhbiTop_JJHB[0]->Write();
  //for(int i=0; i<n; ++i) { hMt_[i]->Write(); rhMt_[i]->Write(); hMt_JJHB_[i]->Write(); rhMt_JJHB_[i]->Write(); }
  outFile->Close();


}
