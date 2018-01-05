#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<map>
  
#include "TH1.h"
#include "TH1F.h"

#include "vector"



void ScaleRes() 
{ 

  gStyle->SetOptStat(0);
  gStyle->SetOptFit();
 
  //set up canvas
  int W = 800;
  int H = 600;
  TString canvName = "resplot";
  TCanvas* canv = new TCanvas(canvName,canvName,50,50,W,H);
  // references for T, B, L, R
  float T = 0.08*H;
  float B = 0.12*H; 
  float L = 0.12*W;
  float R = 0.04*W;
  
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );
  canv->SetTickx(0);
  canv->SetTicky(0);
  

  //set up legend  
  float x1_l = 0.75;
  float y1_l = 0.60;
  
  float dx_l = 0.60;
  float dy_l = 0.1;
  float x0_l = x1_l-dx_l;
  float y0_l = y1_l-dy_l;
  TLegend *lgd = new TLegend(x0_l,y0_l,x1_l, y1_l); 
  lgd->SetBorderSize(0); lgd->SetTextSize(0.04); lgd->SetTextFont(62); lgd->SetFillColor(0);


  //get data
  const unsigned n_f=5;
  unsigned ets[n_f]={10.,20.,30.,50.,100.};
  double scemean[n_f]={0.,0.,0.,0.,0.};
  double scemeane[n_f]={0.,0.,0.,0.,0.};
  double sceres[n_f]={0.,0.,0.,0.,0.};
  double scerese[n_f]={0.,0.,0.,0.,0.};
  double energies[n_f]={0.,0.,0.,0.,0.};
  double energiese[n_f]={0.,0.,0.,0.,0.};

  double sceeta=1.7;
  double scetheta = 2.*atan(exp(-sceeta));

  std::cout<<"will robinson"<<std::endl;
  for(int i(0);i<n_f;i++) {
    energies[i]=ets[i]/sin(scetheta);
    std::cout<<" energy is "<<energies[i]<<std::endl;
    std::ostringstream sceName;
    sceName << "et" << ets[i] << "_eta1.7.root";
    std::cout<<"opening "<<sceName.str().c_str()<<std::endl;
    TFile *f1 = new TFile(sceName.str().c_str());
    TH1F *AAA = static_cast<TH1F*>(f1->Get("h_egenrecopphi")->Clone());
    TFitResultPtr r = AAA->Fit("gaus","S");
    double a0 = r->Parameter(0);
    double a1 = r->Parameter(1);
    double a2 = r->Parameter(2);
    double a0e = r->ParError(0);
    double a1e = r->ParError(1);
    double a2e = r->ParError(2);

    scemean[i]=a1;
    scemeane[i] = a1e;
    sceres[i]=a2/a1;
    scerese[i] = sceres[i]*sqrt((a2e/a2)**2+(a1e/a1)**2);

    std::cout<<" fit results for "<<energies[i]<<" mean "<<scemean[i]<<"+-"<<scemeane[i]<<" res "<<sceres[i]<<"+-"<<scerese[i]<<std::endl;

  }

  TGraphErrors *gr = new TGraphErrors(n_f,energies,sceres,energiese,scerese);
  gr->SetTitle("pion resolution versus energy");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  TF1  *f2 = new TF1("f2","sqrt(([0]/sqrt(x))**2+([1]/x)**2+([2])**2)");
  gr->Fit("f2","q");
  gr->Draw("AP");
  
  TLatex* t = new TLatex();
  t->SetTextFont(62);
  t->SetTextColor(36);
  t->SetTextSize(0.04);
  t->DrawLatex(50.,0.042," energy (GeV)");
  t->SetTextAngle(90);
  t->DrawLatex(-20.,0.08," width/mean");

  canv->Update();
  canv->Print(canvName+".pdf",".pdf");

  return;
}


