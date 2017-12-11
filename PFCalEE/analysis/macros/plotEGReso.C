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
#include "TProfile.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TGaxis.h"

#include "TDRStyle.h"

struct FitResult{
  double chi2;
  unsigned ndf;
  double mean;
  double sigma;
  double meanerr;
  double sigmaerr;
};

double E(const unsigned pT, const unsigned eta){
  return pT*cosh(eta/10.);
};

double pT(const unsigned E, const unsigned eta){
  return E/cosh(eta/10.);
};

void drawChi2(TCanvas *myc,TH1F ** p_chi2ndf){
  
  gStyle->SetOptStat("eMRuo");
  gStyle->SetStatH(0.4);
  gStyle->SetStatW(0.4);
  for (unsigned iSR(0); iSR<8;++iSR){
    myc->cd(iSR+1);
    if (p_chi2ndf[iSR]) p_chi2ndf[iSR]->Draw();
  }
  
  myc->Update();
  std::ostringstream lsave;
  lsave << "PLOTS/EnergyFitQuality.pdf";
  myc->Print(lsave.str().c_str());
}

TPad* plot_ratio(TPad *canv, bool up){
  canv->SetFillColor      (0);
  canv->SetBorderMode     (0);
  canv->SetBorderSize     (10);
  // Set margins to reasonable defaults
  canv->SetLeftMargin     (0.18);
  canv->SetLeftMargin     (0.17);
  canv->SetRightMargin    (0.05);
  canv->SetTopMargin      (0.12);
  canv->SetBottomMargin   (0.18);
  // Setup a frame which makes sense
  canv->SetFrameFillStyle (0);
  canv->SetFrameLineStyle (0);
  canv->SetFrameBorderMode(0);
  canv->SetFrameBorderSize(10);
  canv->SetFrameFillStyle (0);
  canv->SetFrameLineStyle (0);
  canv->SetFrameBorderMode(0);
  canv->SetFrameBorderSize(10);      

  canv->cd();
  TPad *pad = 0;
  if (up){
    pad = new TPad("upper","pad",0, 0.26 ,1 ,1);
    pad->SetBottomMargin(0.05);
    pad->SetTopMargin(0.09);
    pad->Draw();
    pad->cd();
    return pad;
  }
  else {
    pad = new TPad("lower","pad",0, 0   ,1 ,0.26);  
    pad->SetTopMargin(0.05);
    pad->SetBottomMargin(0.24);
    pad->Draw();
    return pad;
  }

};


unsigned fitEnergy(TH1F *hist,
		   TPad *pad,
		   std::string unitStr,
		   FitResult & lres,
		   unsigned isr){
  
  pad->cd();
  //double eMin = hist->GetMean()-5*hist->GetRMS();
  //double eMax = hist->GetMean()+5*hist->GetRMS();
  //hist->GetXaxis()->SetRangeUser(eMin,eMax);
  hist->Draw("PE");


  double nRMSm = isr<1? 1 : 2;
  double nRMSp = 2;
  
  TF1 *fitResult = new TF1("fitResult","[0]*TMath::Gaus(x,[1],[2],0)",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  double tmpmean = hist->GetXaxis()->GetBinCenter(hist->GetMaximumBin());
  /*  if (tmpmean<3) {
    //find next maximum
    std::cout << " !!ERROR!! mismatch of events leading to peak at 0 !! " << std::endl;
    double max = 0;
    for (unsigned ix(4); ix<hist->GetNbinsX()+1; ++ix){
      if (hist->GetBinContent(ix)>max) {
	max = hist->GetBinContent(ix);
	tmpmean = hist->GetXaxis()->GetBinCenter(ix);
      }
    }
    }*/
  fitResult->SetParameters(hist->Integral(),
			   tmpmean,
			   hist->GetRMS());

  std::cout << " Initial params: "  << fitResult->GetParameter(0) << " "<< fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
  << std::endl;


  int status = hist->Fit("fitResult","L0QEMI","",
			 fitResult->GetParameter(1)-nRMSm*fitResult->GetParameter(2),
			 fitResult->GetParameter(1)+nRMSp*fitResult->GetParameter(2));
  
  
  //std::cout << " First fit: " << status << " " << fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
  //<< std::endl;
  
  if ((status != 0 && status != 4000) || fitResult->GetChisquare()/fitResult->GetNDF()>20){
    std::cout << " -- Bad fit ! Try again..." << std::endl;
    status = hist->Fit("fitResult","L0QEMI","",
   		       fitResult->GetParameter(1)-1*fitResult->GetParameter(2),
		       fitResult->GetParameter(1)+2*fitResult->GetParameter(2));
    
    //   std::cout << " Second fit: " << status << " " << fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
    // 	      << std::endl;
  }
  
  // std::cout << " Final fit: " << fitResult->GetParameter(1) << " " << fitResult->GetParameter(2)
  // 	    << std::endl;
  
  fitResult->SetLineColor(2);
  fitResult->Draw("same");
  
  if (status != 0 && status != 4000) {
    std::cout << " Warning! Fit failed with status " << status << "! Please have a look at the verbose output below...." << std::endl;
    hist->Fit("fitResult","L0EMI","",
	      fitResult->GetParameter(1)-nRMSm*fitResult->GetParameter(2),
	      fitResult->GetParameter(1)+nRMSp*fitResult->GetParameter(2));
    //totalE for pu140 is expected to be pathological :/
    //if (isr!=7) return 1;
  }

  char buf[500];
  TLatex lat;
  double latx = hist->GetXaxis()->GetXmin()+(hist->GetXaxis()->GetXmax()-hist->GetXaxis()->GetXmin())/20.;
  double laty = hist->GetMaximum();
  sprintf(buf,"<Efit> = %3.3f +/- %3.3f %s",fitResult->GetParameter(1),fitResult->GetParError(1),unitStr.c_str());
  lat.DrawLatex(latx,laty*0.9,buf);
  sprintf(buf,"RMSfit = %3.3f +/- %3.3f %s",fitResult->GetParameter(2),fitResult->GetParError(2),unitStr.c_str());
  lat.DrawLatex(latx,laty*0.8,buf);
  sprintf(buf,"RMS/meanfit = %3.3f",fitResult->GetParameter(2)/fitResult->GetParameter(1));
  lat.DrawLatex(latx,laty*0.7,buf);
  
  sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fitResult->GetChisquare(),fitResult->GetNDF(),fitResult->GetChisquare()/fitResult->GetNDF());
  lat.DrawLatex(latx,laty*0.6,buf);
  
  lres.chi2 = fitResult->GetChisquare();
  lres.ndf = fitResult->GetNDF();
  lres.mean = fitResult->GetParameter(1);
  lres.meanerr = fitResult->GetParError(1);
  lres.sigma = fitResult->GetParameter(2);
  lres.sigmaerr = fitResult->GetParError(2);

  return 0;
};

TPad* plotCalibration(TGraphErrors *gr,TPad *pad,bool doRatio, TGraphErrors *grDelta,std::string unit, double & calib,double & calibErr, double & offset, double & offsetErr,const unsigned eta, const bool dovsE){

  TPad *upper = 0;
  TPad *lower = 0;

  if (!doRatio) pad->cd();
  else {
    pad->Clear();
    upper = plot_ratio(pad, true);
    lower = plot_ratio(pad, false);
    upper->cd();
  }
  gr->GetXaxis()->SetLabelSize(0.0);
  gr->GetXaxis()->SetTitleSize(0.0);
  gr->GetYaxis()->SetLabelSize(0.06);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetTitleOffset(7);
  gr->GetYaxis()->SetTitleOffset(0.9);
  
  gr->Draw("ap");

  if (!doRatio){
    if (!dovsE) gr->GetXaxis()->SetTitle("p_{T} (GeV)");
    else gr->GetXaxis()->SetTitle("E (GeV)");
  }
  else gr->GetXaxis()->SetTitle("");
  
  gr->GetYaxis()->SetTitle(("Average energy deposited ("+unit+")").c_str()); 
  char buf[500];
  TF1 *fitFunc=new TF1("calib","[0]+[1]*x",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());
  fitFunc->SetLineColor(6);

  //if (dovsE) gr->Fit(fitFunc,"RIME","same",0,200);
  //else gr->Fit(fitFunc,"RIME","same",0,pT(200,eta));
  if (dovsE) gr->Fit(fitFunc,"IME","same");
  else gr->Fit(fitFunc,"IME","same");
  TLatex lat;
  lat.SetTextColor(6);
  lat.SetTextSize(0.1);
  if (!dovsE) sprintf(buf,"<E> #propto a + b #times p_{T} ");
  else sprintf(buf,"<E> #propto a + b #times E ");
  lat.DrawLatexNDC(0.2,0.85,buf);
  sprintf(buf,"a = %3.3f #pm %3.3f %s",fitFunc->GetParameter(0),fitFunc->GetParError(0),unit.c_str());
  lat.DrawLatexNDC(0.2,0.7,buf);
  sprintf(buf,"b = %3.3f #pm %3.3f %s/GeV",fitFunc->GetParameter(1),fitFunc->GetParError(1),unit.c_str());
  lat.DrawLatexNDC(0.2,0.55,buf);
  sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
  lat.DrawLatexNDC(0.2,0.4,buf);

  calib = fitFunc->GetParameter(1);
  offset = fitFunc->GetParameter(0);
  calibErr = fitFunc->GetParError(1);
  offsetErr = fitFunc->GetParError(0);

  if (doRatio){
    //draw deltaE/E vs E
    lower->cd();
    gPad->SetLogx(0);
    gPad->SetGridx(1);
    gPad->SetGridy(1);

    double loffset = fitFunc->GetParameter(0);
    double lslope = fitFunc->GetParameter(1);
    double range = 0.05;
    if (unit=="GeV") {
      loffset=0;
      lslope=1;
      range = 0.1;
    }

    //fill delta
    for (int ip(0);ip<gr->GetN();++ip){
      double x=0;
      double y=0;
      gr->GetPoint(ip,x,y);
      grDelta->SetPoint(ip,x,((y-loffset)/lslope-x)/x);
      double err = gr->GetErrorY(ip)/lslope*1./x;
      grDelta->SetPointError(ip,0,err);
      std::cout << "Calib " << ip << " Egen=" << x << " Erec=" << y << " delta=" << ((y-loffset)/lslope-x)/x << std::endl;
    }
    grDelta->SetTitle("");
    grDelta->SetMinimum(-1.*range);
    grDelta->SetMaximum(range);
    grDelta->GetXaxis()->SetLabelSize(0.15);
    grDelta->GetXaxis()->SetTitleSize(0.15);
    grDelta->GetYaxis()->SetLabelSize(0.12);
    grDelta->GetYaxis()->SetTitleSize(0.15);
    grDelta->GetXaxis()->SetTitleOffset(0.5);
    grDelta->GetYaxis()->SetTitleOffset(0.3);
    
    grDelta->Draw("ap");
    //grDelta->GetYaxis()->SetRangeUser(0,grDelta->GetYaxis()->GetXmax());
    if (!dovsE) grDelta->GetXaxis()->SetTitle("p_{T} (GeV)");
    else grDelta->GetXaxis()->SetTitle("E (GeV)");
    grDelta->GetYaxis()->SetTitle("(#Delta E)/E");

    TLine *line = new TLine(grDelta->GetXaxis()->GetXmin(),0,grDelta->GetXaxis()->GetXmax(),0);
    line->SetLineColor(2);//kYellow+4);
    line->Draw();
    

    lower->Update();

  }

  return upper;
};

bool plotResolution(TGraphErrors *gr,TPad *pad,
		    const unsigned ipu,
		    const unsigned eta,
		    const double & stoch0,
		    const double & const0,
		    const double & noise0,
		    double & stoch,double & stochErr, 
		    double & constant, double & constErr,
		    double & noise,double & noiseErr,
		    const bool dovsE){

  pad->cd();
  gr->GetXaxis()->SetLabelSize(0.06);
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetLabelSize(0.06);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetTitleOffset(0.7);
  gr->GetYaxis()->SetTitleOffset(0.8);
  gr->SetMinimum(0);
  gr->SetMaximum(0.3);
  gr->Draw("ap");
  if (!dovsE) gr->GetXaxis()->SetTitle("p_{T} (GeV)");
  else gr->GetXaxis()->SetTitle("E (GeV)");
  gr->GetYaxis()->SetTitle("#sigma/E");

  TF1 *fitFunc =new TF1("reso","sqrt([0]*[0]/x+[1]*[1]+[2]*[2]/(x*x))",gr->GetXaxis()->GetXmin(),gr->GetXaxis()->GetXmax());

  fitFunc->SetParameter(0,stoch0);
  fitFunc->SetParLimits(0,0,1);
  fitFunc->SetParameter(1,const0);
  fitFunc->SetParLimits(1,0,1);
  fitFunc->SetParameter(2,noise0/2.);
  //fitFunc->SetParLimits(2,0,noise0);

  if (ipu<2) 
    fitFunc->FixParameter(2,noise0);
  if (ipu>=2) fitFunc->FixParameter(1,const0);
  
  std::cout << " Initial params: "  << fitFunc->GetParameter(0) << " "<< fitFunc->GetParameter(1) << " " << fitFunc->GetParameter(2)
	    << std::endl;

  int status = gr->Fit(fitFunc,"BIME0");

  if (fitFunc->GetChisquare()/fitFunc->GetNDF()>50){
    fitFunc->ReleaseParameter(1);
    fitFunc->SetParameters(stoch0,const0,noise0);
    if (!dovsE) status = gr->Fit(fitFunc,"BIME0");//,"",5,110);
    else status = gr->Fit(fitFunc,"BIME0");//,"",10,400);
  }

  fitFunc->SetLineColor(6);
  fitFunc->Draw("same");

  stoch = fitFunc->GetParameter(0);
  stochErr = fitFunc->GetParError(0);
  constant = fitFunc->GetParameter(1);
  constErr = fitFunc->GetParError(1);
  noise = fitFunc->GetParameter(2);
  noiseErr = fitFunc->GetParError(2);

  char buf[500];
  TLatex lat;
  lat.SetTextSize(0.1);
  lat.SetTextColor(6);
  if (!dovsE) sprintf(buf,"#frac{#sigma}{E} #propto #frac{s}{#sqrt{p_{T}}} #oplus c #oplus #frac{n}{p_{T}}");
  else sprintf(buf,"#frac{#sigma}{E} #propto #frac{s}{#sqrt{E}} #oplus c #oplus #frac{n}{E}");

  lat.DrawLatexNDC(0.5,0.85,buf);
  sprintf(buf,"s=%3.3f #pm %3.3f",stoch,stochErr);
  lat.DrawLatexNDC(0.5,0.7,buf);
  sprintf(buf,"c=%3.3f #pm %3.3f",constant,constErr);
  lat.DrawLatexNDC(0.5,0.55,buf);
  sprintf(buf,"n=%3.3f #pm %3.3f",noise,noiseErr);
  lat.DrawLatexNDC(0.5,0.4,buf);
  //sprintf(buf,"status = %d, #chi^{2}/N = %3.1f/%d = %3.1f",status,fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
  sprintf(buf,"#chi^{2}/N = %3.1f/%d = %3.1f",fitFunc->GetChisquare(),fitFunc->GetNDF(),fitFunc->GetChisquare()/fitFunc->GetNDF());
  lat.DrawLatexNDC(0.5,0.25,buf);
  
  if (status != 0 && status != 4000) {
    std::cout << " -- Fit failed with status " << status << std::endl;
    return false;
  }

  return true;
};

bool retrievePuSigma(TTree *atree, TTree *atreePu, 
		     TH1F **p_sigma,
		     TH2F **p_subtrEvspuE,
		     double * calib,
		     double etaval){

  unsigned nLayers = 28;
  unsigned nSR=5;
  if (!atree || !atreePu) {
    std::cout << " -- Info, both trees were not found:" << atree << " " << atreePu << std::endl;
    return true;
  }
  else std::cout << " Trees found." << std::endl;

  unsigned evtIndex = 0;
  unsigned evtIndexPu = 0;
  atree->SetBranchAddress("eventIndex",&evtIndex);
  atreePu->SetBranchAddress("eventIndex",&evtIndexPu);

  if (atree->GetBranchStatus("eventIndex") != 1 || 
      atreePu->GetBranchStatus("eventIndex") != 1) {
    std::cout << " -- Error! Branch eventIndex not found: " << atree->GetBranchStatus("eventIndex") << " " << atreePu->GetBranchStatus("eventIndex") << std::endl;
    return false;
  }
  double totalE = 0;
  double totalEpu = 0;
  atree->SetBranchAddress("wgtEtotal",&totalE);
  atreePu->SetBranchAddress("wgtEtotal",&totalEpu);

  std::vector<std::vector<double> > energySR[2];
  std::vector<std::vector<double> > subtractedenergySR;
  std::vector<double> absweight;

  std::vector<double> emptyvec;
  emptyvec.resize(nSR,0);
  energySR[0].resize(nLayers,emptyvec);
  energySR[1].resize(nLayers,emptyvec);
  subtractedenergySR.resize(nLayers,emptyvec);
  absweight.resize(nLayers,1);

  std::ostringstream label;
  for (unsigned iL(0); iL<nLayers;++iL){
    absweight[iL] = 10.0166;
    //label.str("");
    //label << "absweight_" << iL;
    //atree->SetBranchAddress(label.str().c_str(),&absweight[iL]);
    //atreePu->SetBranchAddress(label.str().c_str(),&absweight[iL]);

    for (unsigned iSR(0);iSR<nSR;++iSR){
      label.str("");
      label << "energy_" << iL << "_SR" << iSR;
      atree->SetBranchAddress(label.str().c_str(),&energySR[0][iL][iSR]);
      atreePu->SetBranchAddress(label.str().c_str(),&energySR[1][iL][iSR]);
      label.str("");
      label << "subtractedenergy_" << iL << "_SR" << iSR;
      atreePu->SetBranchAddress(label.str().c_str(),&subtractedenergySR[iL][iSR]);
    }
  }
  absweight[0] = 20.3628;
  absweight[nLayers-1] = 13.0629;
  int nEvtsPu = atreePu->GetEntries();
  int nEvts = atree->GetEntries();
  if (nEvts!=nEvtsPu) {
    std::cout << " -- Warning, not all events found! nEvts = " << nEvts << " nEvtsPu=" << nEvtsPu << std::endl;
    //return;
  }

  int ievtpu(0);

  unsigned nSkipped = 0;

  for (int ievt(0); ievt<nEvts && ievtpu<nEvtsPu; ++ievt,++ievtpu){//loop on entries
    atree->GetEntry(ievt);
    atreePu->GetEntry(ievtpu);
    
    if (ievt%100 == 0) {
      //if (nSkipped < 10) 
      //std::cout << "... Processing entry: evt " << ievt << " idx=" << evtIndex << " / pu evt " << ievtpu << " idx=" << evtIndexPu << std::endl;
    }

    if (evtIndexPu!=evtIndex) {
      //std::cout << " -- Different events found: " << evtIndex << " " << evtIndexPu;
      if (evtIndexPu<evtIndex) {
	//std::cout << ". Skipping in 140pu tree." << std::endl;
	ievt--;
      }
      else {
	//std::cout << ". Skipping in 0 pu tree." << std::endl;
	ievtpu--;
      }
      nSkipped++;
      continue;
    }
    //if (energySR[1][14][2]<energySR[0][14][2] && evtIndex<10) {
    //std::cout << "nopu " << ievt << " idx=" << evtIndex << " pu " << ievtpu << " idx=" << evtIndexPu << std::endl
    //		<< "SR2 E[14] = " << energySR[0][14][2] << " " << energySR[1][14][2] << std::endl;
    //	}

    double puE[nSR];
    double E[nSR];
    double subtrE[nSR];
    for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
      puE[iSR]=0;
      E[iSR]=0;
      subtrE[iSR]=0;
      for (unsigned iL(0);iL<nLayers;++iL){
	E[iSR]   += absweight[iL]*energySR[0][iL][iSR];
	subtrE[iSR]   += absweight[iL]*subtractedenergySR[iL][iSR];
	puE[iSR] += absweight[iL]*energySR[1][iL][iSR];
      }//loop on layers
      p_sigma[iSR]->Fill((puE[iSR]-E[iSR])/calib[iSR]);
      p_subtrEvspuE[iSR]->Fill((puE[iSR]-E[iSR])/calib[iSR],(puE[iSR]-subtrE[iSR])/calib[iSR]);
    }//loop on SR

    //p_sigma[nSR+2]->Fill((totalEpu-totalE)/tanh(etaval)*1./calib[nSR+2]);
    if (p_sigma[nSR]) p_sigma[nSR]->Fill((totalEpu-totalE)/calib[nSR]);

  }//loop on events

  for (unsigned iSR(0); iSR<nSR+1;++iSR){
    if (p_sigma[iSR]) std::cout << " -- SR " << iSR << " - Found " << p_sigma[iSR]->GetEntries() << " events for pu subtraction sigma! nEvts = " << nEvts << " nEvtsPu=" << nEvtsPu << " nskipped = " << nSkipped << std::endl;
  }

  return true;

};

int plotEGReso(){//main

  SetTdrStyle();

  bool dovsE = true;
  bool processNoFitFiles = false;

  bool doBackLeakCor = false;

  const unsigned nIC = 1;
  const unsigned ICval[nIC] = {3};//0,1,2,3,4,5,10,15,20,50};

  const unsigned nPu = 4;//4;
  unsigned pu[nPu] = {0,0,140,200};//,140,200};

  const unsigned nS = 1;
  std::string scenario[nS] = {
    "model2/gamma/"
  };

  std::string foutname = "PLOTS/PuSubtraction.root";
  TFile *fout = TFile::Open(foutname.c_str(),"RECREATE");

  const unsigned neta = 3;
  unsigned eta[neta]={17,20,24};
  //const unsigned neta = 7;
  //unsigned eta[neta]={17,19,21,23,25,27,29};

  int mycolor[7] = {1,2,4,6,7,8,9};

  double etaval[neta];
  double etaerr[neta];

  const unsigned nEvtMin = 150;

  TString pSuffix = doBackLeakCor?"backLeakCor":"";

  
  const unsigned nV = 1;
  TString version[nV] = {"63"};//,"0"};
  
  const unsigned nLayers = 28;
  const unsigned nSR = 6;
  const double radius[nSR] = {13,15,20,23,26,53};
  double noise100[nSR];
  double noise200[nSR];
  double noise300[nSR];
  unsigned ncellsSmall[nSR] = {7,13,19,31,37,151};
  unsigned ncellsLarge[nSR] = {7,7,13,19,19,85};
  double fitQual[nSR];
  double srval[nSR];
  double srerr[nSR];
  fitQual[0] = 50;
  for (unsigned iSR(0); iSR<nSR;++iSR){
    srval[iSR] = iSR*1.;
    srerr[iSR] = 0.;
    if (iSR>0) fitQual[iSR] = 30;
    noise100[iSR] = sqrt(pow(sqrt(10*ncellsSmall[iSR])*0.00192,2)+pow(sqrt(10*ncellsSmall[iSR])*0.00241,2)+pow(sqrt(8*ncellsSmall[iSR])*0.00325,2));
    noise200[iSR] = sqrt(pow(sqrt(10*ncellsLarge[iSR])*0.00097,2)+pow(sqrt(10*ncellsLarge[iSR])*0.00121,2)+pow(sqrt(8*ncellsLarge[iSR])*0.00164,2));
    noise300[iSR] = sqrt(pow(sqrt(10*ncellsLarge[iSR])*0.00049,2)+pow(sqrt(10*ncellsLarge[iSR])*0.00062,2)+pow(sqrt(8*ncellsLarge[iSR])*0.00083,2));
    std::cout << " Noise vals for iSR " << iSR << " = " << noise100[iSR] << " " << noise200[iSR] << " " << noise300[iSR] << " GeV." << std::endl;
  }


  double calib[nPu][neta][nSR];
  double calibErr[nPu][neta][nSR];
  double offset[nPu][neta][nSR];
  double offsetErr[nPu][neta][nSR];

  double sigmaStoch[nPu][neta][nSR];
  double sigmaStochErr[nPu][neta][nSR];
  double sigmaConst[nPu][neta][nSR];
  double sigmaConstErr[nPu][neta][nSR];
  double sigmaNoise[nPu][neta][nSR];
  double sigmaNoiseErr[nPu][neta][nSR];
  
  std::ostringstream saveName;

  //unsigned genEnAll[]={3,5,7,10,20,30,40,50,60,70,80,90,100,125,150,175,200};
  unsigned genEnAll[]={5,10,20,30,40,60,80,100,150,200};
  //unsigned genEnAll[]={7,10,20,30,40};
  const unsigned nGenEnAll=sizeof(genEnAll)/sizeof(unsigned);


  //canvas so they are created only once
  TCanvas *mycE[nGenEnAll];
  for (unsigned iE(0); iE<nGenEnAll;++iE){
    std::ostringstream lName;
    lName << "mycE" << genEnAll[iE];
    mycE[iE] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    mycE[iE]->Divide(3,2);
  }
  TCanvas *mycE2D[nGenEnAll];
  if (doBackLeakCor){
    for (unsigned iE(0); iE<nGenEnAll;++iE){
      std::ostringstream lName;
      lName << "mycE2D" << genEnAll[iE];
      mycE2D[iE] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
      mycE2D[iE]->Divide(3,2);
    }
  }

  TCanvas *mycSig = new TCanvas("mycSig","puSigma",1);
  //mycSig->Divide(2,2);
  TCanvas *mycCalibEta = new TCanvas("mycCalibEta","mycCalibEta",1500,1000);
  TCanvas *mycOffsetEta = new TCanvas("mycOffsetEta","mycOffsetEta",1500,1000);
  
  const unsigned nCanvas = 4;  
  TCanvas *myc[nCanvas];
  for (unsigned iC(0);iC<nCanvas;++iC){
    std::ostringstream lName;
    lName << "myc" << iC;
    myc[iC] = new TCanvas(lName.str().c_str(),lName.str().c_str(),1500,1000);
    myc[iC]->Divide(2,3);
  }

  TCanvas *mycReso = new TCanvas("mycReso","mycReso",750,750);
  TCanvas *mycR = new TCanvas("mycR","Sampling",1500,1000);
  TCanvas *mycC = new TCanvas("mycC","Constant",1500,1000);
  TCanvas *mycN = new TCanvas("mycN","Noise",1500,1000);

  TH1F *p_chi2ndf[nSR];
  for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
    std::ostringstream label;
    label << "p_chi2ndf_SR" << iSR;
    p_chi2ndf[iSR] = new TH1F(label.str().c_str(),";#chi^{2}/N;entries",500,0,50);
    p_chi2ndf[iSR]->StatOverflows();
  }


  for (unsigned ic(0);ic<nIC;++ic){//loop on intercalib
    
    for (unsigned iV(0); iV<nV;++iV){//loop on versions
      for (unsigned iS(0); iS<nS;++iS){//loop on scenarios
	TString plotDir = "/afs/cern.ch/work/a/amagnan/PFCalEEAna/HGCalTDR/gittestV8/version"+version[iV]+"/"+scenario[iS]+"/";

	TFile *fcalib;
	std::ostringstream label;
	label << plotDir ;
	label << "/CalibReso";
	if (dovsE) label << "_vsE";
	label << "_IC" << ICval[ic];
	label << pSuffix;
	label << ".root";
	fcalib = TFile::Open(label.str().c_str(),"RECREATE");
	
	TTree *ltree[neta][nPu][nGenEnAll];
	TGraphErrors *resoRecoFit[nPu][neta][nSR];
	
	for (unsigned ieta(0);ieta<neta;++ieta){//loop on eta
	  
	  etaval[ieta] = eta[ieta]/10.;
	  etaerr[ieta] = 0;
	  double backLeakCor[nGenEnAll][nSR];

	  for (unsigned ipu(0); ipu<nPu; ++ipu){//loop on pu
	    //if (ieta<2 && ipu==2) continue;
	    unsigned puOption = pu[ipu];

	    std::string unit = "MIPS";
	    if (ipu>0) unit = "GeV";
	    
	    //identify valid energy values
	    bool skip[nGenEnAll];
	    unsigned nValid = 0;
	    for (unsigned iE(0); iE<nGenEnAll; ++iE){
	      ltree[ieta][ipu][iE] = 0;
	      skip[iE] = false;
	      TFile *inputFile = 0;
	      std::ostringstream linputStr;
	      linputStr << plotDir ;
	      if (ipu>=2 && eta[ieta]==17) linputStr << "300u/";
	      else if (ipu>=2 && eta[ieta]==20) linputStr << "200u/";
	      linputStr << "eta" << eta[ieta] << "_et" << genEnAll[iE];// << "_pu" << pu[ipu];
	      if (ipu>=2) linputStr << "_Pu" << pu[ipu];
	      linputStr << "_IC" << ICval[ic];// << "_Si2";
	      if (!processNoFitFiles) linputStr << ".root";
	      else linputStr << "_nofit.root";
	      inputFile = TFile::Open(linputStr.str().c_str());
	      if (!inputFile) {
		std::cout << " -- Error, input file " << linputStr.str() << " cannot be opened. Skipping..." << std::endl;
		//	    return 1;
		skip[iE] = true;
	      }
	      else {
		inputFile->cd("Energies");
		ltree[ieta][ipu][iE] = (TTree*)gDirectory->Get("Ereso");
		
		if (!ltree[ieta][ipu][iE]){
		  std::cout << " -- File " << inputFile->GetName() << " sucessfully opened but tree Ereso not found! Skipping." << std::endl;
		  skip[iE] = true;
		} else { 
		  std::cout << " -- File " << inputFile->GetName() << " sucessfully opened and tree found." << std::endl;
		  nValid++;
		}
	      }
	    }
	    
	    unsigned newidx = 0;
	    const unsigned nGenEn = nValid;
	    unsigned genEn[nGenEn];
	    unsigned oldIdx[nGenEn];
	    for (unsigned iE(0); iE<nGenEnAll; ++iE){	  
	      if (!skip[iE]) {
		genEn[newidx]=genEnAll[iE];
		oldIdx[newidx]=iE;
		newidx++;
	      }
	    }
	    
	    TH1F *p_Ereco[nGenEn][nSR];
	    TH2F *p_ErecovsEback[nGenEn][nSR];
	    TGraphErrors *calibRecoFit[nSR];
	    TGraphErrors *calibRecoDelta[nSR];
	    TGraphErrors *corrBackLeakFit[nSR];

	    //draw calibration curves
	    for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
	      TString srStr = "";
	      srStr += iSR;
	      srStr += "eta";
	      srStr += eta[ieta];
	      srStr += "pu";
	      srStr += ipu;
	      calibRecoFit[iSR] = new TGraphErrors();
	      calibRecoFit[iSR]->SetName("calibRecoFit"+srStr);
	      calibRecoFit[iSR]->SetTitle("");
	      calibRecoFit[iSR]->SetMarkerStyle(20);
	      calibRecoFit[iSR]->SetMarkerColor(1);
	      calibRecoFit[iSR]->SetLineColor(1);
	      calibRecoDelta[iSR] = (TGraphErrors *) calibRecoFit[iSR]->Clone("calibRecoDelta"+srStr);
	      resoRecoFit[ipu][ieta][iSR] = (TGraphErrors *) calibRecoFit[iSR]->Clone("resoRecoFit"+srStr);

	      if (doBackLeakCor) corrBackLeakFit[iSR] = (TGraphErrors *) calibRecoFit[iSR]->Clone("corrBackLeakFit"+srStr);

	    }
	    
	    //get calib and offset from 0 pu file for each SR.
	    
	    std::ostringstream lsave;
	    lsave << "eta"<< eta[ieta]  << "_pu" << puOption;
	    fout->mkdir(lsave.str().c_str());
	    fout->cd(lsave.str().c_str());
	    gStyle->SetOptStat("eMRuo");
	    TH1F *p_sigma[nSR];
	    TH2F *p_subtrEvspuE[nSR];
	    if (pu[ipu]>0){
	      for (unsigned iSR(0); iSR<nSR;++iSR){
		label.str("");
		label << "p_sigma_" << iSR << "_" << lsave.str();
		//if (iSR<(nSR-1) ) 
		p_sigma[iSR] = new TH1F(label.str().c_str(),";PuE-E (GeV)",500,-20,20);
		//else p_sigma[iSR] = new TH1F(label.str().c_str(),";PuE-E (GeV)",5000,0,10000);
		label.str("");
		label << "p_subtrEvspuE_" << iSR << "_" << lsave.str();
		p_subtrEvspuE[iSR] = new TH2F(label.str().c_str(),";PuE-E (GeV);PuE-PuEsubtr (GeV)",
					      500,-20,20,500,-20,20);
		//p_sigma[iSR]->StatOverflows();
	      }
	    }
	    
	    for (unsigned iE(0); iE<nGenEn; ++iE){
	      
	      std::cout << "- Processing energy : " << genEn[iE] 
			<< std::endl;
	      
	      TFile *inputFile = 0;
	      std::ostringstream linputStr;
	      linputStr << plotDir ;
	      if (ipu>=2 && eta[ieta]==17) linputStr << "300u/";
	      else if (ipu>=2 && eta[ieta]==20) linputStr << "200u/";
	      linputStr << "eta" << eta[ieta] << "_et" << genEn[iE];
	      // << "_pu" << pu[ipu];
	      if (ipu>=2) linputStr << "_Pu" << pu[ipu];
	      linputStr << "_IC" << ICval[ic];// << "_Si2";
	      if (!processNoFitFiles) linputStr << ".root" ;
	      else linputStr << "_nofit.root";
	      inputFile = TFile::Open(linputStr.str().c_str());
	      inputFile->cd("Energies");
	      
	      std::cout << " -- Tree entries for eta=" << eta[ieta] << " pu=" << pu[ipu] << " : " << ltree[ieta][ipu][oldIdx[iE]]->GetEntries() << std::endl;
	      
	      for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
		std::cout << " --Processing signal region: " << iSR << std::endl;
		
		mycE[iE]->cd(iSR+1);
		gStyle->SetOptStat(0);
		gStyle->SetOptFit(0);
		
		std::ostringstream lName,lNameTot,lNameBack;
		lName.str("");
		lNameTot.str("");
		lNameBack.str("");
		lName << std::setprecision(6);
		if (iSR<6){
		  for (unsigned iL(0);iL<nLayers;++iL){	      
		    //if (iL==0) lName << "absweight_" << iL  << "*subtractedenergy_" << iL << "_SR" << iSR ;
		    //else lName << "+" << "absweight_" << iL  << "*subtractedenergy_" << iL << "_SR" << iSR ;
		    //if (iL==0) lName << "absweight_" << iL  << "*energy_" << iL << "_SR" << iSR ;
		    //else lName << "+" << "absweight_" << iL  << "*energy_" << iL << "_SR" << iSR ;
		    if (iL==0) lNameTot << "20.3628*energy_" << iL << "_SR" << iSR ;
		    else if (iL==nLayers-1) lNameTot << "+13.0629*energy_" << iL << "_SR" << iSR ;
		    else lNameTot << "+10.0166*energy_" << iL << "_SR" << iSR ;
		    //use always largest area for correction
		    if (iL==nLayers-4) lNameBack << "10.0166*energy_" << iL << "_SR" << iSR ;
		    else if (iL==nLayers-1) lNameBack << "+13.0629*energy_" << iL << "_SR" << iSR ;
		    else if (iL>nLayers-4) lNameBack << "+10.0166*energy_" << iL << "_SR" << iSR ;
		  }
		}
		else lNameTot << "wgtEtotal";///" << tanh(etaval[ieta]);
		if (iSR==4){
		  //std::cout << lName.str() << std::endl;
		  //return 1;
		}

		if (doBackLeakCor && ipu==1) {
		  mycE2D[iE]->cd(iSR+1);
		  lName.str("");
		  lName << "(";
		  lName << lNameTot.str();
		  lName << " - " << offset[0][ieta][iSR] << ")/" << calib[0][ieta][iSR];
		  lName << ":(" << lNameBack.str() << ")/(" << lNameTot.str() << ")";
		  std::cout << lName.str().c_str() << std::endl;
		  ltree[ieta][ipu][oldIdx[iE]]->Draw(lName.str().c_str(),"","colz");
		  lName.str("");
		  lName << "energy" << genEn[iE] << "_SR" << iSR << "_vsBackFraction";
		  p_ErecovsEback[iE][iSR] = (TH2F*)(gPad->GetPrimitive("htemp"))->Clone(lName.str().c_str()); // 2D
		  //if (iSR==7) p_Ereco[iE][iSR] = (TH1F*)gDirectory->Get("p_wgtEtotal");
		  if (!p_ErecovsEback[iE][iSR]){
		    std::cout << " -- ERROR, pointer for histogram " << lName.str() << " is null." << std::endl;
		    return 1;
		  }

		  p_ErecovsEback[iE][iSR]->SetTitle(";E_{24-27}/E_{tot};E_{tot} (GeV)");
		  lName << "_pfx";
		  TProfile *tmpProf = p_ErecovsEback[iE][iSR]->ProfileX(lName.str().c_str());
		  tmpProf->SetMarkerStyle(20);
		  tmpProf->SetMarkerColor(1);
		  tmpProf->SetLineColor(1);
		  p_ErecovsEback[iE][iSR]->Draw("colz");
		  tmpProf->Draw("PEsame");
		  tmpProf->Fit("pol1","","same");

		  //tmpProf->SetStats(1);
		  //gStyle->SetOptFit(1111);
		  TF1 *fitcor = (TF1*)tmpProf->GetFunction("pol1");
		  if (!fitcor) {
		    std::cout << " Fit failed for back leakage correction" << std::endl;
		    return 1;
		  }

		  backLeakCor[oldIdx[iE]][iSR] = fitcor->GetParameter(1);

		  char buf[500];
		  TLatex lat;
		  sprintf(buf,"E=%3.3f #times f_{back} + %3.3f",backLeakCor[oldIdx[iE]][iSR],fitcor->GetParameter(0));
		  lat.DrawLatexNDC(0.25,0.85,buf);

		  Int_t np=corrBackLeakFit[iSR]->GetN();
		  //if (!dovsE) corrBackLeakFit[iSR]->SetPoint(np,genEn[iE],backLeakCor[oldIdx[iE]][iSR]);
		  corrBackLeakFit[iSR]->SetPoint(np,E(genEn[iE],eta[ieta]),backLeakCor[oldIdx[iE]][iSR]);
		  corrBackLeakFit[iSR]->SetPointError(np,0.0,fitcor->GetParError(1));

		  //if (iSR==4) return 1;
		}

		mycE[iE]->cd(iSR+1);

		lName.str("");
		if (ipu>0) lName << "(";
		lName << lNameTot.str();
		if (ipu>0) lName << " - " << offset[0][ieta][iSR] << ")/" << calib[0][ieta][iSR];
		if (doBackLeakCor && ipu>0 && E(genEn[iE],eta[ieta])>300) {
		  lName << " - " << backLeakCor[oldIdx[iE]][iSR] << "*(" << lNameBack.str() << ")/(" << lNameTot.str() << ")";
		}

		std::ostringstream lcut;
		lcut << lName.str() << ">0.5*";
		lcut << E(genEn[iE],eta[ieta]);

		ltree[ieta][ipu][oldIdx[iE]]->Draw(lName.str().c_str(),lcut.str().c_str(),"");
		lName.str("");
		lName << "energy" << genEn[iE] << "_SR" << iSR ;
		p_Ereco[iE][iSR] = (TH1F*)(gPad->GetPrimitive("htemp"))->Clone(lName.str().c_str()); // 1D
		//if (iSR==7) p_Ereco[iE][iSR] = (TH1F*)gDirectory->Get("p_wgtEtotal");
		if (!p_Ereco[iE][iSR]){
		  std::cout << " -- ERROR, pointer for histogram " << lName.str() << " is null." << std::endl;
		  return 1;
		}
		std::cout << " --- Reco E = entries " << p_Ereco[iE][iSR]->GetEntries() 
			  << " mean " << p_Ereco[iE][iSR]->GetMean() 
			  << " rms " << p_Ereco[iE][iSR]->GetRMS() 
			  << " overflows " << p_Ereco[iE][iSR]->GetBinContent(p_Ereco[iE][iSR]->GetNbinsX()+1)
			  << std::endl;
		
		
		p_Ereco[iE][iSR]->SetTitle((";E ("+unit+");events").c_str());
		
		//take min 20 bins
		//if(p_Ereco[iE][iSR]->GetNbinsX()>40) p_Ereco[iE][iSR]->Rebin(2);
		
		//skip data with too little stat
		if (p_Ereco[iE][iSR]->GetEntries()<nEvtMin) {
		  gPad->Clear();
		  continue;
		}


		TPad *lpad = (TPad*)(mycE[iE]->cd(iSR+1));
		FitResult lres;
		if (fitEnergy(p_Ereco[iE][iSR],lpad,unit,lres,iSR)!=0) return 1;
		lpad->cd();
		char buf[500];
		sprintf(buf,"#gamma p_{T}=%d GeV + PU %d",genEn[iE],pu[ipu]);
		TLatex lat;
		lat.SetTextSize(0.05);
		lat.DrawLatexNDC(0.25,0.965,buf);
		sprintf(buf,"#eta=%3.1f, r = %3.0f mm",etaval[ieta],radius[iSR]);
		lat.SetTextSize(0.06);
		lat.DrawLatexNDC(0.15,0.87,buf);
		if (iSR==4) lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
		p_chi2ndf[iSR]->Fill(lres.chi2/lres.ndf);
		
		//filter out bad points
		if (lres.chi2/lres.ndf > fitQual[iSR]) {
		  std::cout << " --- INFO! Point Egen=" 
			    << genEn[iE] 
			    << " eta=" << etaval[ieta]
			    << " pu=" << pu[ipu]
			    << " skipped, chi2/ndf = "
			    << lres.chi2/lres.ndf
			    << std::endl;
		  continue;
		}
		
		Int_t np=calibRecoFit[iSR]->GetN();
		//if (!dovsE) calibRecoFit[iSR]->SetPoint(np,genEn[iE],lres.mean);
		//else 
		calibRecoFit[iSR]->SetPoint(np,E(genEn[iE],eta[ieta]),lres.mean);
		calibRecoFit[iSR]->SetPointError(np,0.0,lres.meanerr);
		
		//use truth: residual calib won't affect resolution....
		double reso = fabs(lres.sigma/lres.mean);
		//if (ipu>1) reso = fabs(lres.sigma/(dovsE?E(genEn[iE],eta[ieta]):genEn[iE]));
		if (ipu>1) reso = fabs(lres.sigma/(E(genEn[iE],eta[ieta])));
		if (!dovsE) resoRecoFit[ipu][ieta][iSR]->SetPoint(np,genEn[iE],reso);
		else resoRecoFit[ipu][ieta][iSR]->SetPoint(np,E(genEn[iE],eta[ieta]),reso);
		double errFit = reso*sqrt(pow(lres.sigmaerr/lres.sigma,2)+pow(lres.meanerr/lres.mean,2));
		resoRecoFit[ipu][ieta][iSR]->SetPointError(np,0,errFit);
		
	      }//loop on SR
	      
	      if (doBackLeakCor && ipu==1){
		saveName.str("");
		saveName << plotDir << "/ErecovsbackFraction_eta" << eta[ieta] << "_pu" << puOption;
		saveName << "_E" << genEn[iE] << pSuffix;
		mycE2D[iE]->Update();
		mycE2D[iE]->Print((saveName.str().c_str()+pSuffix)+".pdf");
		mycE2D[iE]->Print((saveName.str().c_str()+pSuffix)+".C");
	      }

	      saveName.str("");
	      saveName << plotDir << "/Ereco_eta" << eta[ieta] << "_pu" << puOption;
	      if (ipu==0) saveName << "raw";
	      saveName << "_E" << genEn[iE] << pSuffix;
	      mycE[iE]->Update();
	      mycE[iE]->Print((saveName.str().c_str()+pSuffix)+".pdf");
	      mycE[iE]->Print((saveName.str().c_str()+pSuffix)+".C");
	      
	      //fill sigma PU for lower energies
	      if (ipu>1){// && genEn[iE]>5 && genEn[iE]<40) {
		bool success = retrievePuSigma(ltree[ieta][1][oldIdx[iE]], ltree[ieta][ipu][oldIdx[iE]],
					       p_sigma,
					       p_subtrEvspuE,
					       calib[0][ieta],
					       etaval[ieta]);
		if (!success) return 1;
		mycSig->cd();
		p_subtrEvspuE[2]->Draw("colz");
		
	      }
	      
	      
	    }//loop on energies
	    
	    drawChi2(myc[2],p_chi2ndf);


	    //plot and fit calib
	    for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
	      //if (pu[ipu]!=0 && iSR==nSR-1) continue;
	      TPad *lpad = (TPad*)(myc[0]->cd(iSR+1));
	      TPad *upper = plotCalibration(calibRecoFit[iSR],lpad,
					    true,calibRecoDelta[iSR],
					    unit,
					    calib[ipu][ieta][iSR],
					    calibErr[ipu][ieta][iSR],
					    offset[ipu][ieta][iSR],
					    offsetErr[ipu][ieta][iSR],
					    eta[ieta],true);
	      upper->cd();
	      char buf[500];
	      sprintf(buf,"#gamma #eta=%3.1f + PU %d",etaval[ieta],pu[ipu]);
	      TLatex lat;
	      lat.SetTextSize(0.1);
	      lat.DrawLatexNDC(0.7,0.15,buf);
	      sprintf(buf,"r = %3.0f mm",radius[iSR]);
	      lat.DrawLatexNDC(0.7,0.3,buf);
	      lat.SetTextSize(0.06);
	      if (iSR==0) lat.DrawLatexNDC(0.01,0.94,"HGCAL G4 standalone");
	    
	    }

	    myc[0]->Update();
	    lsave.str("");
	    lsave << plotDir << "/";
	    if (ipu==0) lsave << "CalibMipToGeV";
	    else lsave << "Calib";
	    lsave << "_eta" << eta[ieta] << "_pu" << puOption << pSuffix;
	    //if (dovsE) 
	    lsave << "_vsE";
	    myc[0]->Print((lsave.str()+".pdf").c_str());
	    myc[0]->Print((lsave.str()+".C").c_str());

	    //plot reso
	    for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
	      //if (pu[ipu]!=0 && iSR==nSR-1) continue;
	      TPad *lpad = (TPad*)(myc[1]->cd(iSR+1));

	      double stoch0 = pu[ipu]==0? (dovsE?0.25 : 0.14) : sigmaStoch[1][ieta][iSR];
	      double const0 = pu[ipu]==0? 0.01 : sigmaConst[1][ieta][iSR];

	      //limit range to get more realistic RMS ?
	      //if (pu[ipu]!=0 && iSR<(nSR-1)) p_sigma[iSR]->GetXaxis()->SetRangeUser(-5,15);
	      double noise0 = pu[ipu]==0? (eta[ieta]==24?noise100[iSR]:eta[ieta]==20?noise200[iSR]:eta[ieta]==17?noise300[iSR]:0) : p_sigma[iSR]->GetRMS();
	   

	      bool success = plotResolution(resoRecoFit[ipu][ieta][iSR],lpad,
					    ipu,eta[ieta],
					    stoch0,const0,noise0,
					    sigmaStoch[ipu][ieta][iSR],
					    sigmaStochErr[ipu][ieta][iSR],
					    sigmaConst[ipu][ieta][iSR],
					    sigmaConstErr[ipu][ieta][iSR],
					    sigmaNoise[ipu][ieta][iSR],
					    sigmaNoiseErr[ipu][ieta][iSR],
					    dovsE);
	      lpad->cd();
	      char buf[500];
	      sprintf(buf,"#gamma #eta=%3.1f + PU %d",etaval[ieta],pu[ipu]);
	      TLatex lat;
	      lat.SetTextSize(0.1);
	      lat.DrawLatexNDC(0.15,0.85,buf);
	      sprintf(buf,"r = %3.0f mm",radius[iSR]);
	      lat.DrawLatexNDC(0.15,0.75,buf);
	      if (iSR==6) lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");

	      if (!success) {
		continue;
		//return 1;	    
	      }


	      lsave.str("");
	      lsave << "SR" << iSR;
	      fcalib->mkdir(lsave.str().c_str());
	      fcalib->cd(lsave.str().c_str());
	      calibRecoFit[iSR]->Write();
	      calibRecoDelta[iSR]->Write();
	      resoRecoFit[ipu][ieta][iSR]->Write();
	      if (doBackLeakCor && ipu==1) corrBackLeakFit[iSR]->Write();
	    }//loop on SR

	    myc[1]->Update();
	    lsave.str("");
	    lsave << plotDir << "/";
	    if (ipu==0) lsave << "ResoRaw";
	    else lsave << "Reso";
	    lsave << "_eta" << eta[ieta] << "_pu" << puOption << pSuffix;
	    if (dovsE) lsave << "_vsE";
	    myc[1]->Print((lsave.str()+".pdf").c_str());
	    myc[1]->Print((lsave.str()+".C").c_str());

	    if (doBackLeakCor && ipu==1){
	      for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
		//plot back correction slope
		if (doBackLeakCor && ipu==1){
		  myc[3]->cd(iSR+1);
		  //if (!dovsE) corrBackLeakFit[iSR]->SetTitle(";p_{T} (GeV);back cor");
		  corrBackLeakFit[iSR]->SetTitle(";E (GeV);back cor");
		  corrBackLeakFit[iSR]->Draw("APE");
		}
	      }
	      myc[3]->Update();
	      lsave.str("");
	      lsave << plotDir << "/";
	      lsave << "BackLeakCor";
	      lsave << "_eta" << eta[ieta] << "_pu" << puOption << pSuffix;
	      //if (dovsE) 
	      lsave << "_vsE";
	      myc[3]->Print((lsave.str()+".pdf").c_str());
	      myc[3]->Print((lsave.str()+".C").c_str());
	      //return 1;
	    }

	  }//loop on pu
	}//loop on eta

	fout->Write();
	if (nGenEnAll==1) continue;

	TLatex lat;
	char buf[500];

      
	//plot and fit calib parameters vs eta
	TGraphErrors *grSlope[nPu][nSR];
	TGraphErrors *grOffset[nPu][nSR];

	gStyle->SetOptStat(0);
	//gStyle->SetOptFit(1111);
	gStyle->SetStatX(0.64);//top-right corner
	gStyle->SetStatY(0.85);
	gStyle->SetStatW(0.25);
	gStyle->SetStatH(0.25);
	gStyle->SetStatColor(0);

	for (unsigned ipu(0); ipu<nPu; ++ipu){//loop on pu
	  mycCalibEta->Clear();
	  mycOffsetEta->Clear();
	  mycCalibEta->Divide(3,2);
	  mycOffsetEta->Divide(3,2);

	  for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
	    //if (pu[ipu]!=0 && iSR==nSR-1) continue;
	    grSlope[ipu][iSR] = new TGraphErrors();
	    std::ostringstream llabel;
	    llabel << "grSlope_" << ipu << "_SR" << iSR;
	    grSlope[ipu][iSR]->SetName(llabel.str().c_str());
	    llabel.str("");
	    llabel << ";#eta;Calib slope (" ;
	    if (ipu==0) llabel << "MIPs/GeV)";
	    else llabel << "GeV/GeV)";
	    grSlope[ipu][iSR]->SetTitle(llabel.str().c_str());
	    grSlope[ipu][iSR]->SetMarkerStyle(20+iSR);
	    grSlope[ipu][iSR]->SetMarkerColor(1+iSR);
	    grSlope[ipu][iSR]->SetLineColor(1+iSR);
	    llabel.str("");
	    llabel << "grOffset_" << ipu << "_SR" << iSR;
	    grOffset[ipu][iSR] = (TGraphErrors *)grSlope[ipu][iSR]->Clone(llabel.str().c_str());
	    llabel.str("");
	    llabel << ";#eta;Calib offset (" ;
	    if (ipu==0) llabel << "MIPs)";
	    else llabel << "GeV)";
	    grOffset[ipu][iSR]->SetTitle(llabel.str().c_str());

	    for (unsigned ieta(0);ieta<neta;++ieta){//loop on eta
	      grSlope[ipu][iSR]->SetPoint(ieta,etaval[ieta],calib[ipu][ieta][iSR]);
	      grSlope[ipu][iSR]->SetPointError(ieta,0.0,calibErr[ipu][ieta][iSR]);

	      grOffset[ipu][iSR]->SetPoint(ieta,etaval[ieta],offset[ipu][ieta][iSR]);
	      grOffset[ipu][iSR]->SetPointError(ieta,0.0,offsetErr[ipu][ieta][iSR]);
	    }//loop on eta
	    mycCalibEta->cd(iSR+1);
	    grSlope[ipu][iSR]->Draw();
	    grSlope[ipu][iSR]->Fit("pol2","","same");
	     
	    lat.SetTextSize(0.06);
	    sprintf(buf,"pu %d, r = %3.0f mm",pu[ipu],radius[iSR]);
	    lat.DrawLatexNDC(0.2,0.87,buf);
	    if (iSR==4) lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");

	    mycOffsetEta->cd(iSR+1);
	    grOffset[ipu][iSR]->Draw();
	    grOffset[ipu][iSR]->Fit("pol2","","same");

	    lat.DrawLatexNDC(0.2,0.87,buf);
	    if (iSR==4) lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");

	    std::ostringstream lsave;
	    lsave.str("");
	    lsave << "SR" << iSR;
	    fcalib->cd(lsave.str().c_str());
	    grOffset[ipu][iSR]->Write();
	    grSlope[ipu][iSR]->Write();
	  }//loop on SR
	  mycCalibEta->Update();
	  std::ostringstream lsave;
	  lsave << plotDir << "/CalibSlopevsEta_pu" << ipu;
	  if (ipu==0) lsave << "raw";
	  if (dovsE) lsave << "_vsE";
	  lsave << pSuffix << ".pdf";
	  mycCalibEta->Print(lsave.str().c_str());
	  lsave.str("");
	  lsave << plotDir << "/CalibOffsetvsEta_pu" << ipu;
	  if (ipu==0) lsave << "raw";
	  if (dovsE) lsave << "_vsE";
	  lsave << pSuffix << ".pdf";
	  mycOffsetEta->Update();
	  mycOffsetEta->Print(lsave.str().c_str());

	}//loop on pu

	TGraph *grDummy = new TGraph();
	grDummy->SetName("grDummy");
	if (dovsE){
	  grDummy->SetPoint(0,E(genEnAll[0],eta[0]),0.1);
	  grDummy->SetPoint(1,E(genEnAll[nGenEnAll-1],eta[neta-1]),0.1);
	} else {
	  grDummy->SetPoint(0,genEnAll[0],0.1);
	  grDummy->SetPoint(1,genEnAll[nGenEnAll-1],0.1);
	}
	grDummy->SetLineColor(10);
	grDummy->SetMarkerColor(10);
	grDummy->SetMinimum(0);
	if (!dovsE) grDummy->GetXaxis()->SetTitle("p_{T} (GeV)");
	else grDummy->GetXaxis()->SetTitle("E (GeV)");
	grDummy->GetYaxis()->SetTitle("#sigma/E");
	grDummy->GetYaxis()->SetTitleOffset(1.2);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	for (unsigned ipu(1); ipu<nPu; ++ipu){//loop on pu
	  //if (pu[ipu]==0) 
	  grDummy->SetMaximum(0.12);
	  if (!dovsE) {
	    grDummy->SetMaximum(0.09);
	  }
	  for (unsigned iSR(0); iSR<nSR;++iSR){//loop on signal region
	    //if (pu[ipu]!=0 && iSR==nSR-1) continue;
	    //if (!dovsE && iSR==nSR-1) grDummy->SetMaximum(0.15);
	    if (!dovsE) grDummy->SetMaximum(0.15);
	  
	    mycReso->cd();
	    if (dovsE) gPad->SetLogx(1);
	    gPad->SetGridx(1);
	    gPad->SetGridy(1);
	    TLegend *legeta =  new TLegend(0.7,0.7-0.07*neta,0.94,0.94);
	    legeta->SetFillColor(10);
	    grDummy->Draw("AP");
	    for (unsigned ieta(0);ieta<neta;++ieta){//loop on eta
	      resoRecoFit[ipu][ieta][iSR]->SetLineColor(mycolor[ieta]);
	      resoRecoFit[ipu][ieta][iSR]->SetMarkerColor(mycolor[ieta]);
	      resoRecoFit[ipu][ieta][iSR]->SetMarkerStyle(ieta+20);
	      resoRecoFit[ipu][ieta][iSR]->Draw("PLsame");
	      label.str("");
	      label << "#eta = " << etaval[ieta];
	      legeta->AddEntry(resoRecoFit[ipu][ieta][iSR],label.str().c_str(),"P");
	      //add ref pu=0 curve
	      if (pu[ipu]!=0) {
		resoRecoFit[1][ieta][iSR]->SetLineWidth(2);
		resoRecoFit[1][ieta][iSR]->SetLineStyle(2);
		resoRecoFit[1][ieta][iSR]->SetLineColor(mycolor[ieta]);
		resoRecoFit[1][ieta][iSR]->Draw("Lsame");
		legeta->AddEntry(resoRecoFit[1][ieta][iSR],"Ref pu=0","L");
	      }
	    }//loop on eta
	    legeta->Draw("same");
	    sprintf(buf,"#gamma + PU %d",pu[ipu]);
	    //lat.SetTextSize(0.1);
	    lat.DrawLatexNDC(0.45,0.85,buf);
	    sprintf(buf,"r = %3.0f mm",radius[iSR]);
	    lat.DrawLatexNDC(0.45,0.75,buf);
	    lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
	    mycReso->Update();
	    std::ostringstream lsave;
	    lsave << plotDir << "/Resolution_pu" << pu[ipu] << "_SR" << iSR;
	    if (dovsE) lsave << "_vsE";
	    lsave << pSuffix;
	    mycReso->Print((lsave.str()+".pdf").c_str());
	    mycReso->Print((lsave.str()+".C").c_str());
	  
	  }
	}


	TLegend *leg = new TLegend(0.75,0.8,0.94,0.94);
	leg->SetFillColor(10);

	TPad *mypad[3][neta];

	if (neta>1){
	  mycR->Divide(2,1);
	  mycC->Divide(2,1);
	  mycN->Divide(2,1);
	  TPad *left[3];
	  TPad *right[3];
	  left[0] = (TPad*)mycR->cd(1);
	  right[0] = (TPad*)mycR->cd(2);
	  left[1] = (TPad*)mycC->cd(1);
	  right[1] = (TPad*)mycC->cd(2);
	  left[2] = (TPad*)mycN->cd(1);
	  right[2] = (TPad*)mycN->cd(2);
	  for (unsigned iC(0);iC<3;++iC){
	    left[iC]->Divide(1,4);
	    right[iC]->Divide(1,3);
	    for (unsigned ieta=0; ieta<neta;++ieta){//loop on pt values
	      if (ieta<4) mypad[iC][ieta] = (TPad*)left[iC]->GetPad(ieta+1);
	      else mypad[iC][ieta] = (TPad*)right[iC]->GetPad((ieta-4)+1);
	    }
	  }
	}
	else {
	  mypad[0][0] = (TPad*)mycR->cd();
	  mypad[1][0] = (TPad*)mycC->cd();
	  mypad[2][0] = (TPad*)mycN->cd();
	}

	TGraphErrors *grStoch[nPu-1][neta];
	TGraphErrors *grConst[nPu-1][neta];
	TGraphErrors *grNoise[nPu-1][neta];
	for (unsigned ieta(0); ieta<neta;++ieta){
	
	  for (unsigned ipu(0); ipu<( (nPu>1)?(nPu-1):nPu ); ++ipu){//loop on pu
	  
	    grStoch[ipu][ieta] = new TGraphErrors(nSR,srval,sigmaStoch[ipu+1][ieta],srerr,sigmaStochErr[ipu+1][ieta]);
	    grConst[ipu][ieta] = new TGraphErrors(nSR,srval,sigmaConst[ipu+1][ieta],srerr,sigmaConstErr[ipu+1][ieta]);
	    grNoise[ipu][ieta] = new TGraphErrors(nSR,srval,sigmaNoise[ipu+1][ieta],srerr,sigmaNoiseErr[ipu+1][ieta]);

	    TGraphErrors *gr=0;
	    for (unsigned iP(0);iP<3;++iP){
	      gr = (iP==0) ? grStoch[ipu][ieta] : (iP==1) ? grConst[ipu][ieta] : grNoise[ipu][ieta];
	      if (!gr) continue;
	      if (iP==0) mypad[0][ieta]->cd();
	      else if (iP==1) mypad[1][ieta]->cd();
	      else mypad[2][ieta]->cd();
	      gPad->SetGridy(1);
	      gr->SetLineColor(ipu+1);
	      gr->SetMarkerColor(ipu+1);
	      gr->SetMarkerStyle(ipu+21);
	      gr->SetMinimum(0);
	      if (iP<2) gr->SetMaximum(iP==0?0.5:0.06);
	      else gr->SetMaximum(5);
	      if (iP==0) gr->SetTitle(";SR;sampling term (GeV^{#frac{1}{2}})");
	      else if (iP==1) gr->SetTitle(";SR;constant term");
	      else gr->SetTitle(";SR;noise term (GeV)");
	      gr->Draw( (ipu==0) ? "APEL" : "PEL");
	      label.str("");
	      if (ieta==0 && iP==0){
		label.str("");
		label << "PU " << pu[ipu+1];
		leg->AddEntry(gr,label.str().c_str(),"P");
	      }
	      if (ipu==nPu-2){
		label.str("");
		label << "#eta=" << etaval[ieta];
		lat.SetTextSize(0.06);
		lat.DrawLatexNDC(0.2,0.85,label.str().c_str());
		leg->Draw("same");
	      }
	      if (ieta==3 && ipu==0) lat.DrawLatexNDC(0.01,0.01,"HGCAL G4 standalone");
	    }
	  }
	}
      
	mycR->Update();
	std::ostringstream lsave;
	lsave << plotDir << "/SamplingTerm_vsSR" << pSuffix;
	if (dovsE) lsave << "_vsE";
	mycR->Print((lsave.str()+".pdf").c_str());
	mycR->Print((lsave.str()+".png").c_str());

	mycC->Update();
	lsave.str("");
	lsave << plotDir << "/ConstantTerm_vsSR" << pSuffix;
	if (dovsE) lsave << "_vsE";
	mycC->Print((lsave.str()+".pdf").c_str());
	mycC->Print((lsave.str()+".png").c_str());
      
	mycN->Update();
	lsave.str("");
	lsave << plotDir << "/NoiseTerm_vsSR" << pSuffix;
	if (dovsE) lsave << "_vsE";
	mycN->Print((lsave.str()+".pdf").c_str());
	mycN->Print((lsave.str()+".png").c_str());
      
      
	fcalib->Write();
      }//loop on scenarios
      
    }//loop on versions
    
  }//loop on IC vals
  
  
  return 0;
  
}//main
