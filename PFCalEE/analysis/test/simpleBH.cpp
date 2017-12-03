#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<map>
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/function.hpp"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSParameters.hh"
#include "HGCSSCalibration.hh"
#include "HGCSSDigitisation.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"
#include "HGCSSPUenergy.hh"
#include "HGCSSSimpleHit.hh"

#include "PositionFit.hh"
#include "SignalRegion.hh"

#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"

#include "utilities.h"

using boost::lexical_cast;
namespace po=boost::program_options;

int main(int argc, char** argv){//main  


  //Input output and config options
  std::string cfg;
  unsigned pNevts;
  std::string inFilePath;
  std::string outFilePath;
  unsigned nRuns;
  unsigned debug;
  double etamean;
  double deta;
  po::options_description preconfig("Configuration"); 
  preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
  po::notify(vm);
  po::options_description config("Configuration");
  config.add_options()
    //Input output and config options //->required()
    ("pNevts,n",       po::value<unsigned>(&pNevts)->default_value(0))
    ("inFilePath,i",   po::value<std::string>(&inFilePath)->required())
    ("outFilePath,o",  po::value<std::string>(&outFilePath)->required())
    ("nRuns",        po::value<unsigned>(&nRuns)->default_value(0))
    ("debug,d",        po::value<unsigned>(&debug)->default_value(0))
    ("etamean,e",      po::value<double>(&etamean)->default_value(2.8))
    ("deta",      po::value<double>(&deta)->default_value(0.05))
    ;
  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);

  std::cout << " -- Input parameters: " << std::endl
	    << " -- Input file path: " << inFilePath << std::endl
	    << " -- Output file path: " << outFilePath << std::endl
	    << " -- mean eta: " << etamean 
	    << std::endl
	    << " -- Processing ";
  if (pNevts == 0) std::cout << "all events." << std::endl;
  else std::cout << pNevts << " events per run." << std::endl;

  TRandom3 lRndm(0);
  std::cout << " -- Random number seed: " << lRndm.GetSeed() << std::endl;


  //
  // hardcoded
  //


  //global threshold to reduce size of noise hits
  const double threshMin = 0.5;

  std::cout << " ---- Selection settings: ---- " << std::endl
	    << " -------threshMin " << threshMin << std::endl
	    << " ------------------------------------------" << std::endl;



  /////////////////////////////////////////////////////////////
  //input
  /////////////////////////////////////////////////////////////
  TChain *lTree = new TChain("RecoTree");
  TFile * recFile = 0;
  if (nRuns == 0){
    if (!testInputFile(inFilePath,recFile)) return 1;
    lTree->AddFile(inFilePath.c_str());
  }
  else {
    for (unsigned i(0);i<nRuns;++i){
      std::ostringstream lstrrec;
      lstrrec << inFilePath << "_run" << i << ".root";
      if (!testInputFile(lstrrec.str(),recFile)) continue;
      lTree->AddFile(lstrrec.str().c_str());
    }
  }
  if (!lTree){
    std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  std::cout << " Trees added." << std::endl;


  HGCSSInfo * info=(HGCSSInfo*)recFile->Get("Info");

  assert(info);

  const unsigned versionNumber = info->version();
  const unsigned model = info->model();
  const unsigned shape = info->shape();
  const double cellSize = info->cellSize();
  const double calorSizeXY = info->calorSizeXY();

  bool isTBsetup = (model != 2);
  bool bypassR = false;
  if (isTBsetup) bypassR = true;

  HGCSSDetector & myDetector = theDetector();
  myDetector.buildDetector(versionNumber,true,false,bypassR);

  //corrected for Si-Scint overlap
  const unsigned nLayers = 52;//etamean<2.3? myDetector.nLayers(): 52;

  std::cout << " -- Calor size XY = " << calorSizeXY
	    << ", version number = " << versionNumber 
	    << ", model = " << model << std::endl
	    << " -- cellSize = " << cellSize
	    << ", shape = " << shape
	    << ", nLayers = " << nLayers
	    << std::endl;
  HGCSSGeometryConversion geomConv(model,cellSize,bypassR,3);

  geomConv.setXYwidth(calorSizeXY);
  geomConv.setVersion(versionNumber);
  
  if (shape==2) geomConv.initialiseDiamondMap(calorSizeXY,10.);
  else if (shape==3) geomConv.initialiseTriangleMap(calorSizeXY,10.*sqrt(2.));
  else if (shape==1) geomConv.initialiseHoneyComb(calorSizeXY,cellSize);
  else if (shape==4) geomConv.initialiseSquareMap(calorSizeXY,10.);

  //square map for BHCAL
  geomConv.initialiseSquareMap1(1.4,3.0,0,2*TMath::Pi(),0.01745);//eta phi segmentation
  geomConv.initialiseSquareMap2(1.4,3.0,0,2*TMath::Pi(),0.02182);//eta phi segmentation
  std::vector<unsigned> granularity;
  granularity.resize(myDetector.nLayers(),1);
  geomConv.setGranularity(granularity);
  geomConv.initialiseHistos();

  //////////////////////////////////////////////////
  //////////////////////////////////////////////////
  ///////// Output File // /////////////////////////
  //////////////////////////////////////////////////
  //////////////////////////////////////////////////

  TFile *outputFile = TFile::Open(outFilePath.c_str(),"RECREATE");
  
  if (!outputFile) {
    std::cout << " -- Error, output file " << outFilePath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
  }
  outputFile->cd();

  TTree *miptree=new TTree("MipTree","HGC MipAnalysis tree");
  HGCSSSimpleHitVec miphitvec;
  const unsigned nHitsInit = 20000;
  miphitvec.resize(nHitsInit);//,dummy);

  std::ostringstream label;
  //root doesn't like . in branch names.....
  label << "hits";
  miptree->Branch(label.str().c_str(),"std::vector<HGCSSSimpleHit>",&miphitvec);

  TH1F* h_energy = new TH1F("h_energy","hit energy",1000,0.,5.);
  TH1F* h_z = new TH1F("h_z","z of hit",1000,3000.,5200);
  TH2F* h_xy = new TH2F("h_xy","xy of hit",1000,-1200,1200,1000,-1200,1200);

  TH2F* h_xy36 = new TH2F("h_xy36","xy of hit",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_xy37 = new TH2F("h_xy37","xy of hit",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_xy38 = new TH2F("h_xy38","xy of hit",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_xy39 = new TH2F("h_xy39","xy of hit",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_xy40 = new TH2F("h_xy40","xy of hit",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_xy41 = new TH2F("h_xy41","xy of hit",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_xy42 = new TH2F("h_xy42","xy of hit",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_xy43 = new TH2F("h_xy43","xy of hit",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_xy44 = new TH2F("h_xy44","xy of hit",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_xy45 = new TH2F("h_xy45","xy of hit",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_xy46 = new TH2F("h_xy46","xy of hit",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_xy47 = new TH2F("h_xy47","xy of hit",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_xy48 = new TH2F("h_xy48","xy of hit",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_xy49 = new TH2F("h_xy49","xy of hit",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_xy50 = new TH2F("h_xy50","xy of hit",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_xy51 = new TH2F("h_xy51","xy of hit",1000,-1200,1200,1000,-1200,1200);

  
  ///////////////////////////////////////////////////////
  //////////////////  start event loop
  //////////////////////////////////////////////////////


  const unsigned nEvts = ((pNevts > lTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lTree->GetEntries()) : pNevts) ;
  
  std::cout << " -- Processing " << nEvts << " events out of " << lTree->GetEntries() << std::endl;


  //loop on events
  HGCSSEvent * event = 0;
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  unsigned nPuVtx = 0;

  lTree->SetBranchAddress("HGCSSEvent",&event);
  lTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  if (lTree->GetBranch("nPuVtx")) lTree->SetBranchAddress("nPuVtx",&nPuVtx);

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    lTree->GetEntry(ievt);

    miphitvec.clear();


    bool isScint = false;
    if (debug) std::cout << " - Event contains " << (*rechitvec).size() << " rechits." << std::endl;
    unsigned notInEtaRange = 0;
    unsigned notInLayerRange = 0;
    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
	HGCSSRecoHit lHit = (*rechitvec)[iH];
	double leta = lHit.eta();
	if (debug>1) std::cout << " -- hit " << iH << " eta " << leta << std::endl; 
	//clean up rechit collection
	if (fabs(leta-etamean)>= deta){
	  rechitvec->erase(rechitvec->begin()+iH);
	  --iH;
	  notInEtaRange++;
	  continue;
	}
	unsigned layer = lHit.layer();
	if (layer>=myDetector.nLayers()) {
	  rechitvec->erase(rechitvec->begin()+iH);
	  --iH;
	  notInLayerRange++;
	  continue;
	}
	const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
	isScint = subdet.isScint;
	TH2Poly *map = isScint?(subdet.type==DetectorEnum::BHCAL1?geomConv.squareMap1():geomConv.squareMap2()): shape==4?geomConv.squareMap() : shape==2?geomConv.diamondMap() : shape==3? geomConv.triangleMap(): geomConv.hexagonMap();

	unsigned cellid = map->FindBin(lHit.get_x(),lHit.get_y());
	geomConv.fill(lHit.layer(),lHit.energy(),0,cellid,lHit.get_z());

	HGCSSSimpleHit myHit;
	myHit.setE(lHit.energy());
	myHit.setx(lHit.get_x());
	myHit.sety(lHit.get_y());
	myHit.setz(lHit.get_z());
	myHit.setLayer(lHit.layer());
        miphitvec.push_back(myHit);

	h_energy->Fill(lHit.energy());
	h_z->Fill(lHit.get_z());
	h_xy->Fill(lHit.get_x(),lHit.get_y());
	int ilayer = lHit.layer();
	if(ilayer==36) h_xy36->Fill(lHit.get_x(),lHit.get_y());
	if(ilayer==37) h_xy37->Fill(lHit.get_x(),lHit.get_y());
	if(ilayer==38) h_xy38->Fill(lHit.get_x(),lHit.get_y());
	if(ilayer==39) h_xy39->Fill(lHit.get_x(),lHit.get_y());
	if(ilayer==40) h_xy40->Fill(lHit.get_x(),lHit.get_y());
	if(ilayer==41) h_xy41->Fill(lHit.get_x(),lHit.get_y());
	if(ilayer==42) h_xy42->Fill(lHit.get_x(),lHit.get_y());
	if(ilayer==43) h_xy43->Fill(lHit.get_x(),lHit.get_y());
	if(ilayer==44) h_xy44->Fill(lHit.get_x(),lHit.get_y());
	if(ilayer==45) h_xy45->Fill(lHit.get_x(),lHit.get_y());
	if(ilayer==46) h_xy46->Fill(lHit.get_x(),lHit.get_y());
	if(ilayer==47) h_xy47->Fill(lHit.get_x(),lHit.get_y());
	if(ilayer==48) h_xy48->Fill(lHit.get_x(),lHit.get_y());
	if(ilayer==49) h_xy49->Fill(lHit.get_x(),lHit.get_y());
	if(ilayer==50) h_xy50->Fill(lHit.get_x(),lHit.get_y());
	if(ilayer==51) h_xy51->Fill(lHit.get_x(),lHit.get_y());

      }//loop on hits

      if (debug) std::cout << " - In eta range, event contains " << (*rechitvec).size() << " rechits." << std::endl;

    
    miptree->Fill();

    geomConv.initialiseHistos();

  }//loop on entries


  outputFile->cd();
  outputFile->Write();
  outputFile->Close();

  return 0;
  
}//main
