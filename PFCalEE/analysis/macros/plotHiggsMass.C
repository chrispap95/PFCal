#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TGaxis.h"

#include "effSigmaMacro.C"

int plotHiggsMass(){//main

  const unsigned nP = 3;//2;
  unsigned puVal[nP] = {0,140,200};

  const unsigned nV = 1;//7;

  std::string label[nV] = {
    //"True E, Shower pos",
    //"True E, Vtx smear",
    //"True E, Shower angle",
    "Reco E, True pos",
    //"Reco E, Shower pos",
    //"Reco E, Vtx smear",
    //"Reco E, Shower angle",
  };

  double val[nV];
  double valerr[nV];
  for (unsigned iV(0);iV<nV;++iV){
    val[iV] = iV+0.5;
    valerr[iV] = 0;
  }

  //fill arrays
  double mean[nP][nV];
  double meanerr[nP][nV];
  double sigma[nP][nV];
  double sigmaerr[nP][nV];
  
  TFile *f[nP];
  f[0] = TFile::Open("/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalTDR//gittestV8/version63/model2/HggLarge/pu0.root");
  f[1] = TFile::Open("/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalTDR//gittestV8/version63/model2/HggLarge/pu140.root");
  f[2] = TFile::Open("/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalTDR//gittestV8/version63/model2/HggLarge/pu200.root");

  TCanvas *mycfit[3];
  mycfit[0] = new TCanvas("myc1","myc1",1500,1000);
  mycfit[1] = new TCanvas("myc2","myc2",1500,1000);
  mycfit[2] = new TCanvas("myc3","myc3",1500,1000);
  if (nV==3){
    mycfit[0]->Divide(2,2);
    mycfit[1]->Divide(2,2);
    mycfit[2]->Divide(2,2);
  }
  else {
    mycfit[0]->Divide(2,1);
    mycfit[1]->Divide(2,1);
    mycfit[2]->Divide(2,1);
  }
  TCanvas *myc = new TCanvas("myc","myc",1);
  gStyle->SetOptStat("eMRuo");
  gStyle->SetOptFit(1111);
  //gStyle->SetStatH(0.2);
  //gStyle->SetStatW(0.4);
  for (unsigned ipu(0);ipu<nP;++ipu){//loop on pu
    if (!f[ipu]) {
      std::cout << " -- Input file for pu " << puVal[ipu] << " not found ! " << std::endl;
      return 1;
    }
    f[ipu]->cd("HiggsMass");
    TH1F *hMass[nV];
    //hMass[0] = (TH1F*)gDirectory->Get("p_position_trueE");
    //hMass[1] = (TH1F*)gDirectory->Get("p_position_vtxsmear_trueE");
    //hMass[2] = (TH1F*)gDirectory->Get("p_angle_trueE");
    hMass[0] = (TH1F*)gDirectory->Get("p_trueDir_recoE");
    if (!hMass[0]) return 1;
    //hMass[3] = (TH1F*)gDirectory->Get("p_position_recoE");
    //hMass[4] = (TH1F*)gDirectory->Get("p_position_vtxsmear_recoE");
    //hMass[6] = (TH1F*)gDirectory->Get("p_angle_recoE");
   
    TH1F *trueM = (TH1F*)gDirectory->Get("p_trueDir_trueE");
    if (!trueM) return 1;
    mycfit[ipu]->cd(1);
    trueM->Draw();

    for (unsigned iV(0);iV<nV;++iV){//loop on masses
      mycfit[ipu]->cd(2+iV);
      hMass[iV]->GetXaxis()->SetRangeUser(105,145);
      hMass[iV]->Draw();
      hMass[iV]->Fit("gaus","+","same",118,132);
      TF1 *fit = (TF1*)hMass[iV]->GetFunction("gaus");
      if (!fit) return 1;
      mean[ipu][iV] = fit->GetParameter(1);
      meanerr[ipu][iV] = fit->GetParError(1);
      sigma[ipu][iV] = fit->GetParameter(2);
      sigmaerr[ipu][iV] = fit->GetParError(2);
      //mean[ipu][iV] = hMass[iV]->GetMean();
      //meanerr[ipu][iV] = hMass[iV]->GetMeanError();
      //sigma[ipu][iV] = hMass[iV]->GetRMS();
      //sigmaerr[ipu][iV] = hMass[iV]->GetRMSError();
    }
    std::ostringstream lsave;
    lsave << "TDRPLOTS/HiggsMasses_pu" << puVal[ipu] << ".pdf";
    mycfit[ipu]->Print(lsave.str().c_str());

  }

  return 1;

  myc->cd();
  gPad->SetGridy(1);

  gStyle->SetOptStat(0);
  TGraphErrors *grMass[nP];
  TGraphErrors *grSigma[nP];
 
  for (unsigned iH(0); iH<2; ++iH){
    TGraphErrors *gr[nP];

    for (unsigned iP(0); iP<nP; ++iP){
      
      grMass[iP] = new TGraphErrors(nV,val,mean[iP],valerr,meanerr[iP]);
      grSigma[iP] = new TGraphErrors(nV,val,sigma[iP],valerr,sigmaerr[iP]);
      
      grMass[iP]->SetTitle(";;Mass (GeV)");
      grSigma[iP]->SetTitle(";;#sigma_{M} (GeV)");
      
      gr[iP] = iH==0? grMass[iP] : grSigma[iP];
      gr[iP]->SetMarkerStyle(21+iP);
      gr[iP]->SetMarkerColor(iP+1);
      if (iH==1) {
	gr[iP]->SetMinimum(0);
	gr[iP]->SetMaximum(4.5);
      }
      else {
	gr[iP]->SetMinimum(120);
	gr[iP]->SetMaximum(126);
      }
      TAxis *ax = gr[iP]->GetHistogram()->GetXaxis();
      Double_t x1 = ax->GetBinLowEdge(1);
      Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
      gr[iP]->GetHistogram()->GetXaxis()->Set(nV,x1,x2);
      
      for(Int_t k=0;k<nV;k++){
	gr[iP]->GetHistogram()->GetXaxis()->SetBinLabel(k+1,label[k].c_str());
      }
      
      if (iP==0) gr[iP]->Draw("AP");
      else gr[iP]->Draw("P");

      TLatex lat;
      lat.SetTextColor(iP+1);
      char buf[100];
      sprintf(buf,"PU = %d",puVal[iP]);
      lat.DrawLatexNDC(0.2+0.2*iP,0.2,buf);
    }//loop on PU

    TLatex lat;
    lat.DrawLatexNDC(0.1,0.92,"Pythia gg#rightarrow Higgs, H#rightarrow #gamma#gamma");
    lat.DrawLatexNDC(0.2,0.85,"1.6 < #eta_{#gamma^{1}} < 2.8, p_{T}>40 GeV");
    lat.DrawLatexNDC(0.2,0.8,"1.6 < #eta_{#gamma^{2}} < 2.8, p_{T}>40 GeV");
    myc->Update();
    if (iH==0) myc->Print("TDRPLOTS/SummaryHiggsMass.pdf");
    else {
      myc->Print("TDRPLOTS/SummaryHiggsReso.pdf");
      myc->Print("TDRPLOTS/SummaryHiggsReso.C");
    }
  }//loop on histos


  return 0;


}//main
