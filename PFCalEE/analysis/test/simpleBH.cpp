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
  //std::string inFilePath;
  std::string outFilePath;
  std::string filePath;
  std::string digifilePath;
  unsigned nRuns;
  std::string simFileName;
  std::string recoFileName;
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
    ("outFilePath,o",  po::value<std::string>(&outFilePath)->required())
    ("filePath,i",     po::value<std::string>(&filePath)->required())
    ("digifilePath", po::value<std::string>(&digifilePath)->default_value(""))
    ("simFileName,s",  po::value<std::string>(&simFileName)->required())
    ("recoFileName,r", po::value<std::string>(&recoFileName)->required())
    ("nRuns",        po::value<unsigned>(&nRuns)->default_value(0))
    ("debug,d",        po::value<unsigned>(&debug)->default_value(0))
    ("etamean,e",      po::value<double>(&etamean)->default_value(2.8))
    ("deta",      po::value<double>(&deta)->default_value(0.05))
    ;
  po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
  po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
  po::notify(vm);


  std::string inFilePath = filePath+simFileName;

  std::cout << " -- Input parameters: " << std::endl
            << " -- Input file path: " << filePath << std::endl
            << " -- Digi Input file path: " << digifilePath << std::endl
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


  std::ostringstream inputsim;
  inputsim << filePath << "/" << simFileName;
  std::ostringstream inputrec;
  if (digifilePath.size()==0)
    inputrec << filePath << "/" << recoFileName;
  else
    inputrec << digifilePath << "/" << recoFileName;

  std::cout << inputsim.str() << " " << inputrec.str() << std::endl;

  HGCSSInfo * info;

  TChain *lSimTree = new TChain("HGCSSTree");
  TChain *lRecTree = 0;

  TFile * simFile = 0;
  TFile * recFile = 0;

  if (recoFileName.find("Digi") != recoFileName.npos)
    lRecTree = new TChain("RecoTree");
  else lRecTree = new TChain("PUTree");

  if (nRuns == 0){
    if (!testInputFile(inputsim.str(),simFile)) return 1;
    lSimTree->AddFile(inputsim.str().c_str());
    if (simFile) info =(HGCSSInfo*)simFile->Get("Info");
    else {
      std::cout << " -- Error in getting information from simfile!" << std::endl;
      return 1;
    }
    if (!testInputFile(inputrec.str(),recFile)) return 1;
    lRecTree->AddFile(inputrec.str().c_str());
  }
  else {
    for (unsigned i(0);i<nRuns;++i){
      std::ostringstream lstrsim;
      std::ostringstream lstrrec;
      lstrsim << inputsim.str() << "_run" << i << ".root";
      if (testInputFile(lstrsim.str(),simFile)){
        if (simFile) info =(HGCSSInfo*)simFile->Get("Info");
        else {
	  std::cout << " -- Error in getting information from simfile!" << std::endl;
          return 1;
        }
      }
      else continue;
      lstrrec << inputrec.str() << "_run" << i << ".root";
      if (!testInputFile(lstrrec.str(),recFile)) continue;
      lSimTree->AddFile(lstrsim.str().c_str());
      lRecTree->AddFile(lstrrec.str().c_str());
    }
  }

  if (!lSimTree){
    std::cout << " -- Error, tree HGCSSTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  if (!lRecTree){
    std::cout << " -- Error, tree RecoTree cannot be opened. Exiting..." << std::endl;
    return 1;
  }





  //assert(info);

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




  std::ostringstream label;
  //root doesn't like . in branch names.....
  //label << "hits";
  //miptree->Branch(label.str().c_str(),"std::vector<HGCSSSimpleHit>",&miphitvec);

  TH1F* h_energy = new TH1F("h_energy","hit energy",1000,0.,5.);
  TH1F* h_z = new TH1F("h_z","z of hit",5000,3100.,5200);
  TH1F* h_z1 = new TH1F("h_z1","z of hit",5000,3150.,3550);
  TH1F* h_z2 = new TH1F("h_z2","z of hit",5000,3550.,5200);
  TH2F* h_xy = new TH2F("h_xy","xy of hit",1000,-1200,1200,1000,-1200,1200);
  TH1F* h_l = new TH1F("h_l","layer of hit",80,0.,80.);
  TH1F* h_l2 = new TH1F("h_l2","layer of hit",30,50,80.);
  TH2F* h_zl = new TH2F("h_zl","z vs l of hit",5000,4300.,5200,25,30.,55.);

  TH2F* h_nsxy36 = new TH2F("h_nsxy36","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy37 = new TH2F("h_nsxy37","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy38 = new TH2F("h_nsxy38","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy39 = new TH2F("h_nsxy39","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy40 = new TH2F("h_nsxy40","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy41 = new TH2F("h_nsxy41","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy42 = new TH2F("h_nsxy42","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy43 = new TH2F("h_nsxy43","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy44 = new TH2F("h_nsxy44","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy45 = new TH2F("h_nsxy45","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy46 = new TH2F("h_nsxy46","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy47 = new TH2F("h_nsxy47","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy48 = new TH2F("h_nsxy48","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy49 = new TH2F("h_nsxy49","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy50 = new TH2F("h_nsxy50","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_nsxy51 = new TH2F("h_nsxy51","xy of hit not scint",1000,-1200,1200,1000,-1200,1200);

  TH2F* h_sxy36 = new TH2F("h_sxy36","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy37 = new TH2F("h_sxy37","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy38 = new TH2F("h_sxy38","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy39 = new TH2F("h_sxy39","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy40 = new TH2F("h_sxy40","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy41 = new TH2F("h_sxy41","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy42 = new TH2F("h_sxy42","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy43 = new TH2F("h_sxy43","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy44 = new TH2F("h_sxy44","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy45 = new TH2F("h_sxy45","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy46 = new TH2F("h_sxy46","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy47 = new TH2F("h_sxy47","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy48 = new TH2F("h_sxy48","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy49 = new TH2F("h_sxy49","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy50 = new TH2F("h_sxy50","xy of hit scint",1000,-1200,1200,1000,-1200,1200);
  TH2F* h_sxy51 = new TH2F("h_sxy51","xy of hit scint",1000,-1200,1200,1000,-1200,1200);

  TH2F* h_Egenreco = new TH2F("h_Egenreco","E reco sum versus gen",100,0.,50.,100,0.,1.1);
  TH1F* h_egenreco = new TH1F("h_egenreco","E reco sum over gen",100,0.,1.1);

  
  ///////////////////////////////////////////////////////
  //////////////////  start event loop
  //////////////////////////////////////////////////////


  //  const unsigned nEvts = ((pNevts > lRecTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lRecTree->GetEntries()) : pNevts) ;
  
  //std::cout << " -- Processing " << nEvts << " events out of " << 
  //   lRecTree->GetEntries()<< std::endl;


  const unsigned nEvts = ((pNevts > lSimTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(lSimTree->GetEntries()) : pNevts) ;

  std::cout << " -- Processing " << nEvts << " events out of " << lSimTree->GetEntries() << " " << lRecTree->GetEntries() << std::endl;


  //loop on events
  HGCSSEvent * event = 0;
  HGCSSEvent * eventRec = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  std::vector<HGCSSGenParticle> * genvec = 0;
  unsigned nPuVtx = 0;



  lSimTree->SetBranchAddress("HGCSSEvent",&event);
  lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);

  lRecTree->SetBranchAddress("HGCSSEvent",&eventRec);
  lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  if (lRecTree->GetBranch("nPuVtx")) lRecTree->SetBranchAddress("nPuVtx",&nPuVtx);



  unsigned ievtRec = 0;
  unsigned nSkipped = 0;
  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (ievtRec>=lRecTree->GetEntries()) continue;

    if (debug) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;


    lSimTree->GetEntry(ievt);
    lRecTree->GetEntry(ievtRec);
    if (nPuVtx>0 && eventRec->eventNumber()==0 && event->eventNumber()!=0) {
      std::cout << " skip !" << ievt << " " << ievtRec << std::endl;
      nSkipped++;
      continue;
    }

    double ptgen=-1.;
    double Egen=-1.;
    double ptgenpx=-1.;
    double ptgenpy=-1.;
    double ptgenpz=-1.;
    int pidgen=-1;
    if((*genvec).size()>0) {
      pidgen=(*genvec)[0].pdgid();
      ptgenpx=(*genvec)[0].px()/1000.;
      ptgenpy=(*genvec)[0].py()/1000.;
      ptgenpz=(*genvec)[0].pz()/1000.;
      ptgen=sqrt(ptgenpx*ptgenpx+ptgenpy*ptgenpy);
      Egen=sqrt(ptgenpx*ptgenpx+ptgenpy*ptgenpy+ptgenpz*ptgenpz);
    }
    if(debug) {
      std::cout<<" gen vec size is "<<(*genvec).size()<<std::endl;
      std::cout<<" first gen "<<ptgen<<" "<<Egen<<" "<<pidgen<<std::endl;
      for (unsigned iP(0); iP<(*genvec).size(); ++iP){
        std::cout<<" gen particle "<<iP<<" is "<<(*genvec)[iP].pdgid()<<std::endl;
      }
    }


    bool isScint = false;
    if (debug) std::cout << " - Event contains " << (*rechitvec).size() << " rechits." << std::endl;
    unsigned notInEtaRange = 0;
    unsigned notInLayerRange = 0;
    double rechitsumE=0.;
    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
	HGCSSRecoHit lHit = (*rechitvec)[iH];
	double leta = lHit.eta();
	if (debug>20) std::cout << " -- hit " << iH << " eta " << leta << std::endl; 
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


	rechitsumE+=lHit.energy();
	h_energy->Fill(lHit.energy());
	h_z->Fill(lHit.get_z());
	h_z1->Fill(lHit.get_z());
	h_z2->Fill(lHit.get_z());
	h_l->Fill(lHit.layer()+0.5);
	int ixx=lHit.layer();
	if(ixx>52) ixx=ixx-17;
	h_zl->Fill(lHit.get_z(),ixx);
	h_l2->Fill(lHit.layer()+0.5);
	h_xy->Fill(lHit.get_x(),lHit.get_y());
	int ilayer = ixx;

	if(debug>20) {
	  if(isScint) std::cout<<"Will Robinson: found one "<<lHit.energy()<<" "<<lHit.layer()<<" "<<ilayer<<std::endl;
	}



	if(ilayer==36) {
	  if(isScint) h_sxy36->Fill(lHit.get_x(),lHit.get_y());
	  else h_nsxy36->Fill(lHit.get_x(),lHit.get_y());
	};
	if(ilayer==37) {
	  if(isScint) h_sxy37->Fill(lHit.get_x(),lHit.get_y());
	  else h_nsxy37->Fill(lHit.get_x(),lHit.get_y());
	};
	if(ilayer==38) {
	  if(isScint) h_sxy38->Fill(lHit.get_x(),lHit.get_y());
	  else h_nsxy38->Fill(lHit.get_x(),lHit.get_y());
	};
	if(ilayer==39) {
	  if(isScint) h_sxy39->Fill(lHit.get_x(),lHit.get_y());
	  else h_nsxy39->Fill(lHit.get_x(),lHit.get_y());
	};
	if(ilayer==40) {
	  if(isScint) h_sxy40->Fill(lHit.get_x(),lHit.get_y());
	  else h_nsxy40->Fill(lHit.get_x(),lHit.get_y());
	};
	if(ilayer==41) {
	  if(isScint) h_sxy41->Fill(lHit.get_x(),lHit.get_y());
	  else h_nsxy41->Fill(lHit.get_x(),lHit.get_y());
	};
	if(ilayer==42) {
	  if(isScint) h_sxy42->Fill(lHit.get_x(),lHit.get_y());
	  else h_nsxy42->Fill(lHit.get_x(),lHit.get_y());
	};
	if(ilayer==43) {
	  if(isScint) h_sxy43->Fill(lHit.get_x(),lHit.get_y());
	  else h_nsxy43->Fill(lHit.get_x(),lHit.get_y());
	};
	if(ilayer==44) {
	  if(isScint) h_sxy44->Fill(lHit.get_x(),lHit.get_y());
	  else h_nsxy44->Fill(lHit.get_x(),lHit.get_y());
	};
	if(ilayer==45) {
	  if(isScint) h_sxy45->Fill(lHit.get_x(),lHit.get_y());
	  else h_nsxy45->Fill(lHit.get_x(),lHit.get_y());
	};
	if(ilayer==46) {
	  if(isScint) h_sxy46->Fill(lHit.get_x(),lHit.get_y());
	  else h_nsxy46->Fill(lHit.get_x(),lHit.get_y());
	};
	if(ilayer==47) {
	  if(isScint) h_sxy47->Fill(lHit.get_x(),lHit.get_y());
	  else h_nsxy47->Fill(lHit.get_x(),lHit.get_y());
	};
	if(ilayer==48) {
	  if(isScint) h_sxy48->Fill(lHit.get_x(),lHit.get_y());
	  else h_nsxy48->Fill(lHit.get_x(),lHit.get_y());
	};
	if(ilayer==49) {
	  if(isScint) h_sxy49->Fill(lHit.get_x(),lHit.get_y());
	  else h_nsxy49->Fill(lHit.get_x(),lHit.get_y());
	};
	if(ilayer==50) {
	  if(isScint) h_sxy50->Fill(lHit.get_x(),lHit.get_y());
	  else h_nsxy50->Fill(lHit.get_x(),lHit.get_y());
	};
	if(ilayer==51) {
	  if(isScint) h_sxy51->Fill(lHit.get_x(),lHit.get_y());
	  else h_nsxy51->Fill(lHit.get_x(),lHit.get_y());
	};



      }//loop on hits
    if(debug>1) {
      std::cout<<" reco gen are "<<rechitsumE<<" "<<Egen<<std::endl;
    }
    h_Egenreco->Fill(Egen,rechitsumE/Egen);
    h_egenreco->Fill(rechitsumE/Egen);


      if (debug) std::cout << " - In eta range, event contains " << (*rechitvec).size() << " rechits." << std::endl;

    
      //miptree->Fill();

    geomConv.initialiseHistos();
    ievtRec++;
  }//loop on entries

  if(debug) std::cout<<"writing files"<<std::endl;

  outputFile->cd();
  outputFile->Write();
  outputFile->Close();

  return 0;
  
}//main
