#include "TMath.h"
#include <math.h>
#include "TF1.h"
#include "TF2.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveLabel.h"
#include "TAxis.h"
#include "TH1.h"
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<map>

#include "TH1.h"
#include "TH1F.h"

#include "vector"


//gStyle->SetOptStat(0);
//gStyle->SetOptFit();

/*
Double_t gamma(Double_t *x,Double_t *par) {
  Double_t gammaVal = TMath::Gamma(Double_t par[0], Double_t x[0]);
  return fitval;
}
*/

void Gamma1() {
  //gSystem->Load("libMathCore");
  TFile *f = new TFile("et100_eta1.7_deadfrac0p5_test.root");
  TH2F *h_lay_energy2D = (TH2F*)f->Get("h_lay_energy2D");
  TProfile *prof = h_lay_energy2D->ProfileX();
  //TH1D *oneD = h_lay_energy2D->ProjectionX();
//  TF1 *gamma = new TF1("gamma", "[0]*TMath::GammaDist(x,[1],[2],[3])", 1, 25);
  TF1 *gamma = new TF1("gamma", "TMath::GammaDist(x,[0],[1],[2])", 0, 25);
  TF1 *gamma_pdf = new TF1("gamma_pdf","ROOT::Math::gamma_pdf([0],[1],[2])",0,25);

  Double_t norm = 1;
  Double_t scale = norm/(prof->Integral());
  prof->Scale(scale);

  gamma->SetParameters(3, 0, 0.5);
//  gamma->SetParameters(8., 0, 2.4, 8.5);
  //gamma_pdf->SetParameters(1,0.5,0);
  gamma->SetParNames("gamma","mu","beta","amplitude");
  prof->Fit("gamma","S");
  //Double_t chi2 = r->Chi2();
//  Double_t NDF = r->NDF();
  //r->Print("V");
//  prof->Fit("gamma","S");
}
