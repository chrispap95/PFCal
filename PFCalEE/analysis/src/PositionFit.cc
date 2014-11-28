#include <iomanip>

#include "PositionFit.hh"
#include "Clusterizer.hh"
#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSDetector.hh"

double PositionFit::DeltaPhi(const double & phi1, const double & phi2){
  double dphi = phi1 - phi2;
  if (dphi< (-1.*TMath::Pi())) dphi += 2*TMath::Pi();
  if (dphi>TMath::Pi()) dphi -= 2*TMath::Pi();
  return dphi;
}

std::pair<unsigned, std::pair<double,double> > PositionFit::findMajorityValue(std::vector<std::pair<double,double> > & values) const
{
  
  unsigned lTot = values.size();
  if (!lTot) return std::pair<unsigned, std::pair<double,double> >(0,std::pair<double,double>(0,0));
  
  std::sort(values.begin(),values.end());

  unsigned lMajorityCounter = 0;
  std::pair<double,double> lMaj = std::pair<double,double>(0,0);
  
  std::vector<std::pair<double,double> >::iterator lIter = values.begin();
  for ( ; lIter != values.end(); ) {
    //std::cout << lIter->first << " " << lIter->second << std::endl;
    unsigned lCounter = std::count(lIter,values.end(),*lIter);
    if (lCounter > lMajorityCounter) {
      lMajorityCounter = lCounter;
      lMaj = *lIter;
    }
    lIter += lCounter;
  }
    
  //std::cout << " -- Found majority value " << lMaj.first << " " << lMaj.second << " for " << lMajorityCounter << " elements out of " << values.size() << "." << std::endl;

  return std::pair<unsigned, std::pair<double,double> >(lMajorityCounter,lMaj);
  
}

PositionFit::PositionFit(const unsigned nSR,
			 const double & residualMax, 
			 const unsigned nLayers, 
			 const unsigned nSiLayers,
			 const bool applyPuMixFix,
			 const unsigned debug,
			 const bool doMatrix
			 ){
  nSR_ = nSR;
  residualMax_ = residualMax;
  chi2ndfmax_ = 20;
  seedMipThreshold_ = 10;
  maxdR_ = 0.1;
  nLayers_ = nLayers;
  nSiLayers_ = nSiLayers;
  debug_ = debug;
  useMeanPU_ = true;
  fixForPuMixBug_ = applyPuMixFix;
  doMatrix_ = doMatrix;

  p_nGenParticles = 0;
  //p_numberOfMaxTried = 0;
  //p_dRMaxTruth = 0;
  p_hitMeanPuContrib = 0;
  p_hitEventPuContrib = 0;
  p_diffPuContrib = 0;

  p_etavsphi = 0;
  p_etavsphi_max = 0;
  p_etavsphi_truth = 0;
  p_yvsx_max = 0;
  p_yvsx_truth = 0;

  p_residuals_x = 0;
  p_residuals_y = 0;
  p_errorMatrix = 0;
  p_corrMatrix = 0;
  p_chi2[0] = 0;
  p_chi2[1] = 0;

  p_chi2overNDF[0] = 0;
  p_impactX[0] = 0;
  p_impactY[0] = 0;
  p_impactX_residual = 0;
  p_impactY_residual = 0;
  p_tanAngleX[0] = 0;
  p_tanAngleY[0] = 0;
  p_tanAngleX_residual = 0;
  p_tanAngleY_residual = 0;
  p_positionReso = 0;
  p_angularReso = 0;
  p_chi2overNDF[1] = 0;
  p_impactX[1] = 0;
  p_impactY[1] = 0;
  p_tanAngleX[1] = 0;
  p_tanAngleY[1] = 0;

}

void PositionFit::initialise(TFile *outputFile,
			     const std::string outputDir,
			     const std::string outFolder, 
			     const HGCSSGeometryConversion & geomConv, 
			     const HGCSSPUenergy & puDensity){

  outputDir_ = outputDir;
  setOutputFile(outputFile);

  outFolder_ = outFolder;
  matrixFolder_ = outFolder_;

  geomConv_ = geomConv;
  puDensity_ = puDensity;

  nL_mean_.resize(nLayers_,0);
  mean_[0].resize(nLayers_,0);
  mean_[1].resize(nLayers_,0);

  nL_sigma_.resize(nLayers_,nL_mean_);
  sigma_[0].resize(nLayers_,mean_[0]);
  sigma_[1].resize(nLayers_,mean_[0]);
  for (unsigned iL(0);iL<nLayers_;++iL){
    nL_sigma_[iL].resize(nLayers_,0);
    sigma_[0][iL].resize(nLayers_,0);
    sigma_[1][iL].resize(nLayers_,0);
  }

}

void PositionFit::initialiseClusterHistograms(){

  outputFile_->cd(outputDir_.c_str());

  p_nClusters = new TH1F("p_nClusters",";n_{clusters};n_{events}",1000,0,3000);

  p_clusnHits_all = new TH1F("p_clusnHits_all",";n_{hits};n_{clusters}",500,0,1000);
  p_seedEoverE_all = new TH1F("p_seedEoverE_all",";seedE/E;n_{clusters}",100,0,1);
  p_clusLayer_all = new TH1F("p_clusLayer_all",";cluster layer;n_{clusters}",nLayers_,0,nLayers_);
  p_clusWidth_all = new TH1F("p_clusWidth_all",";cluster width (layers);n_{clusters}",nLayers_,0,nLayers_);
  p_seeddeta_all = new TH1F("p_seeddeta_all",";#Delta#eta(seed,cluster);n_{clusters}",100,-0.5,0.5);
  p_seeddphi_all = new TH1F("p_seeddphi_all",";#Delta#phi(seed,cluster);n_{clusters}",100,-0.5,0.5);

  p_clusnHits_sel = new TH1F("p_clusnHits_sel",";n_{hits};n_{events}",500,0,1000);
  p_seedEoverE_sel = new TH1F("p_seedEoverE_sel",";seedE/E;n_{events}",100,0,1);
  p_clusLayer_sel = new TH1F("p_clusLayer_sel",";cluster layer;n_{events}",nLayers_,0,nLayers_);
  p_clusWidth_sel = new TH1F("p_clusWidth_sel",";cluster width (layers);n_{events}",nLayers_,0,nLayers_);
  p_seeddeta_sel = new TH1F("p_seeddeta_sel",";#Delta#eta(seed,cluster);n_{events}",100,-0.1,0.1);
  p_seeddphi_sel = new TH1F("p_seeddphi_sel",";#Delta#phi(seed,cluster);n_{events}",100,-0.1,0.1);

  p_mindRtruth = new TH1F("p_mindRtruth",";mindR(cluster,truth);n_{events}",100,0,0.5);

}

void PositionFit::initialisePositionHistograms(){
  //for yvsx plots
  unsigned nX=2*1700/10-2;
  double minX=-1.*nX*5-5,maxX=nX*5+5;
  double minY=minX,maxY=maxX;
  nX += 1;
  unsigned nY = nX;

  outputFile_->cd(outputDir_.c_str());

  p_nGenParticles = new TH1F("p_nGenParticles",";nGenParticles",10,0,10);
  //p_numberOfMaxTried = new TH1F("p_numberOfMaxTried",";max cells tried",250,0,250);
  //p_dRMaxTruth = new TH1F("p_dRMaxTruth",";#Delta R(max,truth)",200,0,0.1);
  //p_dRMaxTruth->StatOverflows();

  p_genxy.resize(nLayers_,0);
  p_recoxy.resize(nLayers_,0);
  std::ostringstream lName;
  for (unsigned iL(0); iL<nLayers_; ++iL){
    lName.str("");
    lName << "p_genxy_" << iL;
    p_genxy[iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",
			   nX*10,minX,maxX,
			   nY*10,minY,maxY);
    lName.str("");
    lName << "p_recoxy_" << iL;
    p_recoxy[iL] = new TH2F(lName.str().c_str(),";x(mm);y(mm)",
   			    nX,minX,maxX,
   			    nY,minY,maxY);
  }

  p_hitMeanPuContrib = new TH2F("p_hitMeanPuContrib",";layer;E_{PU} (MIPs) from mean;hits",nLayers_,0,nLayers_,1000,0,50);
  p_hitMeanPuContrib->StatOverflows();

  p_hitEventPuContrib = new TH2F("p_hitEventPuContrib",";layer;E_{PU} (MIPs) from RC;hits",nLayers_,0,nLayers_,1000,0,50);
  p_hitEventPuContrib->StatOverflows();

  p_diffPuContrib = new TH1F("p_diffPuContrib",";E_{PU}^{RC}-E_{PU}^{avg} (MIPs) in SR2;events",1000,-2000,2000);
  p_diffPuContrib->StatOverflows();

  p_etavsphi = new TH2F("p_etavsphi",";#phi_{hit};#eta_{hit};n_{hits}",
			900,-3.1416,3.1416,
			250,1.4,3.0);

  p_etavsphi_max = new TH2F("p_etavsphi_max",";#phi_{max};#eta_{max};n_{events}",
			900,-3.1416,3.1416,
			250,1.4,3.0);
  p_etavsphi_truth = new TH2F("p_etavsphi_truth",";#phi_{gen};#eta_{gen};n_{events}",
			900,-3.1416,3.1416,
			250,1.4,3.0);

  p_yvsx_max = new TH2F("p_yvsx_max",";x_{max};y_{max};n_{events}",
			nX,minX,maxX,nY,minY,maxY);
  p_yvsx_truth = new TH2F("p_yvsx_truth",";x_{gen};y_{gen};n_{events}",
			  nX,minX,maxX,nY,minY,maxY);


  p_residuals_x = new TH1F("p_residuals_x",";xreco-xtruth (mm)",1000,-50,50);
  p_residuals_y = new TH1F("p_residuals_y",";yreco-ytruth (mm)",1000,-50,50);
  p_residuals_x->StatOverflows();
  p_residuals_y->StatOverflows();

}

void PositionFit::initialiseFitHistograms(){

  outputFile_->cd(outputDir_.c_str());

  //check if already defined
  if (!p_chi2[0]){

    p_recoXvsLayer = new TH2F("p_recoXvsLayer",";layer;weighted x (mm);n_{events}",nLayers_,0,nLayers_,200,-100,100);
    p_recoYvsLayer = new TH2F("p_recoYvsLayer",";layer;weighted y (mm);n_{events}",nLayers_,0,nLayers_,1200,300,1500);
    p_recoZvsLayer = new TH2F("p_recoZvsLayer",";layer;avg z (mm);n_{events}",nLayers_,0,nLayers_,3000,3170,3470);
    p_truthXvsLayer = new TH2F("p_truthXvsLayer",";layer;weighted x (mm);n_{events}",nLayers_,0,nLayers_,200,-100,100);
    p_truthYvsLayer = new TH2F("p_truthYvsLayer",";layer;weighted x (mm);n_{events}",nLayers_,0,nLayers_,1200,300,1500);
    p_fitXvsLayer = new TH2F("p_fitXvsLayer",";layer;fit x (mm);n_{events}",nLayers_,0,nLayers_,200,-100,100);
    p_fitYvsLayer = new TH2F("p_fitYvsLayer",";layer;fit y (mm);n_{events}",nLayers_,0,nLayers_,1200,300,1500);
    p_nLayersFit = new TH1F("p_nLayersFit",";#layers in fit;n_{events}",31,-0.5,30.5);

    p_chi2[0] = new TH1F("p_chi2",";#chi^{2};n_{events}",1000,0,5000);
    p_chi2[1] = new TH1F("p_chi2_truth",";#chi^{2};n_{events}",1000,0,5000);
    p_chi2overNDF[0] = new TH1F("p_chi2overNDF",";#chi^{2}/NDF;n_{events}",1000,0,500);
    p_chi2overNDF[1] = new TH1F("p_chi2overNDF_truth",";#chi^{2}/NDF;n_{events}",1000,0,500);
    for (unsigned rt(0); rt<2;++rt){
      p_chi2[rt]->StatOverflows();
      p_chi2overNDF[rt]->StatOverflows();
    }
    
    p_impactX[0] = new TH1F("p_impactX",";x front face impact (mm);n_{events}",500,-100,100);
    p_impactX[1] = new TH1F("p_impactX_truth",";x front face impact (mm);n_{events}",500,-100,100);
    p_impactY[0] = new TH1F("p_impactY",";y front face impact (mm);n_{events}",1200,300,1500);
    p_impactY[1] = new TH1F("p_impactY_truth",";y front face impact (mm);n_{events}",1200,300,1500);
    p_tanAngleX[0] = new TH1F("p_tanAngleX",";x direction tanAngle (rad);n_{events}",500,-1,1);
    p_tanAngleX[1] = new TH1F("p_tanAngleX_truth",";x direction tanAngle (rad);n_{events}",500,-1,1);
    p_tanAngleY[0] = new TH1F("p_tanAngleY",";y direction tanAngle (rad);n_{events}",500,-1,1);
    p_tanAngleY[1] = new TH1F("p_tanAngleY_truth",";y direction tanAngle (rad);n_{events}",500,-1,1);

    p_impactX_residual = new TH1F("p_impactX_residual",";residual x front face impact (mm);n_{events}",200,-10,10);
    p_tanAngleX_residual = new TH1F("p_tanAngleX_residual",";residual x direction tanAngle (rad);n_{events}",200,-0.1,0.1);
    p_impactY_residual = new TH1F("p_impactY_residual",";residual y front face impact (mm);n_{events}",200,-10,10);
    p_tanAngleY_residual = new TH1F("p_tanAngleY_residual",";residual y direction tanAngle (rad);n_{events}",200,-0.1,0.1);

    
    p_positionReso = new TH1F("p_positionReso",";#sigma_{x,y} (mm);n_{events}",500,0,50);
    p_angularReso = new TH1F("p_angularReso",";#sigma_{#theta} (rad);n_{events}",500,0,1);
  }
}

bool PositionFit::getGlobalMaximum(const unsigned ievt, 
				   const unsigned nVtx, 
				   std::vector<HGCSSRecoHit> *rechitvec,
				   const ROOT::Math::XYZVector & truthPos0, 
				   double & aPhimax,double & aEtamax){

  bool oneresult = true;
  TH2F *letavsphi = new TH2F("letavsphi",";#phi;#eta;hits",
			     900,-3.1416,3.1416,
			     250,1.4,3.0);
  
  /*
  TH2F *letavsphi_layer[nLayers_];
  if (nVtx>0){
    std::ostringstream lName;
    for (unsigned iL(0); iL<nLayers_; ++iL){
      lName.str("");
      lName << "letavsphi_layer" << iL;
      letavsphi_layer[iL] = new TH2F(lName.str().c_str(),";#phi;#eta;hits",
				     900,-3.1416,3.1416,
				     250,1.4,3.0);
    }
  }
  */

  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
    const HGCSSRecoHit & lHit = (*rechitvec)[iH];
    
    double posx = lHit.get_x();
    if (fixForPuMixBug_) posx-=1.25;
    double posy = lHit.get_y();
    if (fixForPuMixBug_) posy-=1.25;
    double posz = lHit.get_z();
    //double radius = sqrt(posx*posx+posy*posy);
    double energy = lHit.energy();
    //unsigned layer = lHit.layer();
    //if (energy>1) std::cout << "Hit " << layer << " " << posx << " " << posy << " " << posz << " " << energy << std::endl;
    
    if (debug_>1) {
      std::cout << " --  RecHit " << iH << "/" << (*rechitvec).size() << " --" << std::endl
		<< " --  position x,y " << posx << "," << posy << std::endl;
      lHit.Print(std::cout);
    }
    
    ROOT::Math::XYZVector pos(posx,posy,posz);
    letavsphi->Fill(pos.phi(),pos.eta(),energy);

    //if (nVtx>0){
    //letavsphi_layer[layer]->Fill(pos.phi(),pos.eta(),energy);
    //}
    p_etavsphi->Fill(pos.phi(),pos.eta(),energy);
  }//loop on hits

  if (debug_)  std::cout << std::endl;

  //add histograms for all layers but iL
  //for (unsigned iL(0); iL<nLayers_; ++iL){
  //letavsphi_layer[iL]->Add(letavsphi,letavsphi_layer[iL],1,-1);
  //}
  
  //get position of maximum E tower
  int maxbin = -1;
  int binx=-1,biny=-1,binz=-1;
  //compare with truth
  double etatruth = truthPos0.eta();
  double phitruth = truthPos0.phi();
  double dRmin = 10;
  int binxmin=-1,binymin=-1;
  unsigned counter = 0;
  while (dRmin>maxdR_) {
    letavsphi->SetBinContent(maxbin,0);
    maxbin = letavsphi->GetMaximumBin();
    letavsphi->GetBinXYZ(maxbin,binx,biny,binz);
    double leta = letavsphi->GetYaxis()->GetBinCenter(biny);
    double lphi = letavsphi->GetXaxis()->GetBinCenter(binx);
    double deta = leta-etatruth;
    double dphi = DeltaPhi(lphi,phitruth);
    double dR = sqrt(pow(deta,2)+pow(dphi,2));
    if (dR<dRmin){
      dRmin=dR;
      binxmin=binx;
      binymin=biny;
    }
    if (counter>200) {
      break;
    }
    counter++;
  }

  if (debug_ || (counter>200)) std::cout << " -- Maximum found after " << counter << " trials, dRmin = " << dRmin;
  if (counter>200) std::cout << " --> away from truth pos for evt " << ievt ;
  if (debug_ || (counter>200)) std::cout << std::endl;  

  if (counter>1) oneresult = false;
  aPhimax =letavsphi->GetXaxis()->GetBinCenter(binxmin); 
  aEtamax =letavsphi->GetYaxis()->GetBinCenter(binymin); 

  //p_numberOfMaxTried->Fill(counter);
  //p_dRMaxTruth->Fill(dRmin);
  /*
  if (nVtx>0){
    double binxsize =  letavsphi->GetXaxis()->GetBinWidth(binx);
    double binysize =  letavsphi->GetYaxis()->GetBinWidth(biny);
    
    //allow for +/- 1 bin in each direction
    double dRmax = 2*sqrt(pow(binxsize,2)+pow(binysize,2));
    
    //unsigned counter = 0;
    //std::vector<double> layervec;
    //std::vector<double> etavec;
    //std::vector<double> phivec;
    std::vector<std::pair<double,double> > etaphivec;

    for (unsigned iL(0); iL<nLayers_; ++iL){
      int maxbin_layer = letavsphi_layer[iL]->GetMaximumBin();
      int binxl,binyl,binzl;
      letavsphi_layer[iL]->GetBinXYZ(maxbin_layer,binxl,binyl,binzl);
      double leta = letavsphi_layer[iL]->GetYaxis()->GetBinCenter(binyl);
      double lphi = letavsphi_layer[iL]->GetXaxis()->GetBinCenter(binxl);
      etaphivec.push_back(std::pair<double,double>(leta,lphi));
      //double deta = leta-aEtamax;
      //double dphi = DeltaPhi(lphi,aPhimax);
      //double dR = sqrt(pow(deta,2)+pow(dphi,2));
      //if (dR>dRmax) {
      //layervec.push_back(iL);
      //etavec.push_back(leta);
      //phivec.push_back(lphi);
      //counter++;
      //}
    }
    
    std::pair<unsigned, std::pair<double,double> > lmaj = findMajorityValue(etaphivec);
    double leta = lmaj.second.first;
    double lphi = lmaj.second.second;
    double deta = leta-aEtamax;
    double dphi = DeltaPhi(lphi,aPhimax);
    double dR = sqrt(pow(deta,2)+pow(dphi,2));
    
    if (dR>dRmax) {
      std::cout << " -- Warning ! Event " << ievt << " with " << lmaj.first << " majority layers dRmax=" << dRmax << " away, probably from PU." << std::endl
		<< " All layers 0-29 eta-phi = " << aEtamax << " " << aPhimax << std::endl
		<< " Maj value = " << leta << " " << lphi << std::endl;
      oneresult = false;
    }
    // if (counter > 0){
    //std::cout << " -- Warning ! Event " << ievt << " with " << counter << " layers dRmax=" << dRmax << " away, probably from PU." << std::endl
    // 		<< " All layers 0-29 eta-phi = " << aEtamax << " " << aPhimax << std::endl;
    //   for (unsigned iC(0); iC<counter;++iC){
    // 	std::cout << " Removing layer " << layervec[iC] << " eta-phi = " << etavec[iC] << " " << phivec[iC] << std::endl;
    //   }
  //}
  }//if nvtx>0
  */

  if (debug_) std::cout << " MaxE cell eta,phi = " << aEtamax << " " << aPhimax << std::endl;

  letavsphi->Delete();
  /*if (nVtx>0) {
    for (unsigned iL(0); iL<nLayers_; ++iL){
      letavsphi_layer[iL]->Delete();
    }
    }*/

  return oneresult;
}

void PositionFit::findSeeds(std::vector<HGCSSRecoHit> *rechitvec,
			    std::vector<bool> & seedable){
  
  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
    const HGCSSRecoHit & lHit = (*rechitvec)[iH];
    if (lHit.energy()>seedMipThreshold_) seedable[iH] = true;
  }
  
}

unsigned PositionFit::getClusters(std::vector<HGCSSRecoHit> *rechitvec,
				  HGCSSClusterVec & output){

  const unsigned nHits = (*rechitvec).size();
  Clusterizer lClusterizer(debug_);
  std::vector<bool> rechitMask;
  rechitMask.resize(nHits,true);
  std::vector<bool> seedable;
  seedable.resize(nHits,false);

  findSeeds(rechitvec,seedable);

  lClusterizer.buildClusters(rechitvec,rechitMask,seedable,output);

  return output.size();

}

bool PositionFit::setTruthInfo(std::vector<HGCSSGenParticle> *genvec, const int G4TrackID){

  bool found = false;
  for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles    
    //if ((*genvec).size()!= 1) (*genvec)[iP].Print(std::cout);
    if ((*genvec)[iP].trackID()==G4TrackID){
      found = true;
      //set direction
      truthDir_ = Direction((*genvec)[iP].px()/(*genvec)[iP].pz(),(*genvec)[iP].py()/(*genvec)[iP].pz());
      double zvtx = (*genvec)[iP].z()-((*genvec)[iP].y()/truthDir_.tanangle_y);
      double xvtx = (*genvec)[iP].x()-(((*genvec)[iP].z()-zvtx)*truthDir_.tanangle_x);
      truthVtx_ = ROOT::Math::XYZPoint(xvtx,0,zvtx);
      //in GeV
      truthE_ = (*genvec)[iP].E()/1000.;

      if (debug_) std::cout << " Found truth vertex pos at: x=" << xvtx << " z=" << zvtx << " xpos was: " << (*genvec)[iP].x() << std::endl;

    }

  }

  return found;

}

bool PositionFit::getTruthPosition(std::vector<HGCSSGenParticle> *genvec,std::vector<ROOT::Math::XYPoint> & truthPos, const int G4TrackID){

  bool found = false;
  for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles    
    //if ((*genvec).size()!= 1) (*genvec)[iP].Print(std::cout);
    if ((*genvec)[iP].trackID()==G4TrackID){
      if ((*genvec)[iP].pdgid() != 22) {
	std::cout << " -- Error ! Particle trackID " << G4TrackID << " is not a photon ! pdgid = " << (*genvec)[iP].pdgid() << std::endl;
	break;
      }
      found = true;
      double x0 = (*genvec)[iP].x();
      double y0 = (*genvec)[iP].y();
      double z0 = (*genvec)[iP].z();
      //double p = sqrt(pow((*genvec)[iP].px(),2)+pow((*genvec)[iP].py(),2)+pow((*genvec)[iP].pz(),2));
      //std::cout << "init : " << x0 << " " << y0 << " " << z0 << std::endl;
      //fill layers by propagating with momentum
      //ROOT::Math::XYZVector unit((*genvec)[iP].px()/p,(*genvec)[iP].py()/p,(*genvec)[iP].pz()/p);
      
      //std::cout << " Gen particle eta,phi = " << unit.eta() << " " << unit.phi() << std::endl;
      p_etavsphi_truth->Fill((*genvec)[iP].phi(),(*genvec)[iP].eta());
      //std::cout << " Truth pos eta-phi = " << (*genvec)[iP].eta() << " " << (*genvec)[iP].phi() << std::endl;
      
      for (unsigned iL(0); iL<nLayers_; ++iL){
	//double xy = (avgZ_[iL]-z0)/sinh(unit.eta());
	//double x = xy*cos(unit.phi())+x0;
	//double y = xy*sin(unit.phi())+y0;
	double x = x0+(avgZ_[iL]-z0)*(*genvec)[iP].px()/(*genvec)[iP].pz();
	double y = y0+(avgZ_[iL]-z0)*(*genvec)[iP].py()/(*genvec)[iP].pz();
	//std::cout << "Lay " << iL << ": " << x << " " << y << " " << avgZ_[iL] << std::endl;
	p_genxy[iL]->Fill(x,y,1);
	truthPos[iL] = ROOT::Math::XYPoint(x,y);
      }
      
    }
  }//loop on gen particles

  if (!found){
    std::cout << " - Info: no photon G4trackID=" << G4TrackID << " found, already converted or outside of acceptance ..." << std::endl;
    //std::cout << " -- Number of genparticles: " << (*genvec).size() << std::endl;
    //for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles 
    //std::cout << " --- particle " << iP << std::endl;
    //(*genvec)[iP].Print(std::cout);
    //}
  }
  else {
    p_nGenParticles->Fill((*genvec).size());
  }
  return found;

}

bool PositionFit::getZpositions(){
  std::ifstream fin;
  std::ostringstream finname;
  finname << outFolder_ << "/zPositions.dat";
  fin.open(finname.str());
  if (!fin.is_open()){
    std::cout << " Cannot open input file " << finname.str() << "! Refilling now..." << std::endl;
    return false;
  }
  
  std::cout << " Reading z position per layer from input file " << finname.str() << std::endl;

  std::vector<unsigned> layerId;
  std::vector<double> posz;
  layerId.reserve(nLayers_);
  posz.reserve(nLayers_);

  while (!fin.eof()){
    unsigned l=nLayers_;
    double z=0;
    fin>>l>>z;
    if (l<nLayers_){
      avgZ_.push_back(z);
      std::cout << " Layer " << l << ", z = " << z << std::endl;
    }
  }
  
  if (avgZ_.size() != nLayers_) {
    std::cout << " -- Warning! Problem in extracting z positions, did not find one value per layer. Please check input file: " << finname.str() << std::endl
	      << " Proceeding to refilling them from simtree." << std::endl;
    return false;
  }

  fin.close();
  return true;
 
}

void PositionFit::getZpositions(TTree *aSimTree,
				const unsigned nEvts){
  

  std::vector<HGCSSSimHit> * simhitvec = 0;
  aSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);

  std::ofstream fout;
  std::ostringstream foutname;
  foutname << outFolder_ << "/zPositions.dat";
  fout.open(foutname.str());
  if (!fout.is_open()){
    std::cout << " Cannot open outfile " << foutname.str() << " for writing ! Exiting..." << std::endl;
    exit(1);
  }
  
  std::cout << "--- Filling z positions:" << std::endl
	    << "- Processing = " << nEvts  << " events out of " << aSimTree->GetEntries() << std::endl;

  avgZ_.resize(nLayers_,0);


  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug_) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    aSimTree->GetEntry(ievt);
    for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on rechits
      const HGCSSSimHit & lHit = (*simhitvec)[iH];
      unsigned layer = lHit.layer();
      if (layer >= nLayers_) {
	continue;
      }

      //discard some si layers...
      if (lHit.silayer() >= nSiLayers_) continue; 

      double posz = lHit.get_z();
      //get z position of hits
      if (avgZ_[layer]<posz) avgZ_[layer]=posz;
    }
  }

  std::cout << " --- Z positions of layers: " << std::endl;
  for (unsigned iL(0); iL<nLayers_;++iL){
    std::cout << " Layer " << iL << ", z = " << avgZ_[iL] << std::endl;
    fout << iL << " " << avgZ_[iL] << std::endl;
  }

  fout.close();
  
}

void PositionFit::getInitialPositions(TTree *aSimTree, 
				      TTree *aRecTree,
				      const unsigned nEvts,
				      const unsigned G4TrackID){

  initialiseClusterHistograms();
  initialisePositionHistograms();
  //HGCSSDetector & myDetector = theDetector();

  //////////////////////////////////////////////////
  ///////// Event loop /////////////////////////////
  //////////////////////////////////////////////////

  HGCSSEvent * event = 0;
  std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSSimHit> * simhitvec = 0;
  std::vector<HGCSSRecoHit> * rechitvec = 0;
  std::vector<HGCSSGenParticle> * genvec = 0;
  unsigned nPuVtx = 0;

  aSimTree->SetBranchAddress("HGCSSEvent",&event);
  aSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  aSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  aSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);
  
  aRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  if (aRecTree->GetBranch("nPuVtx")) aRecTree->SetBranchAddress("nPuVtx",&nPuVtx);

  std::cout << "- Processing = " << nEvts  << " events out of " << aSimTree->GetEntries() << std::endl;

  //bool firstEvent = true;

  unsigned nConvertedPhotons = 0;
  unsigned nTooFar = 0;

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    if (debug_) std::cout << "... Processing entry: " << ievt << std::endl;
    else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
    
    aSimTree->GetEntry(ievt);

    aRecTree->GetEntry(ievt);

    if (debug_) std::cout << " nPuVtx = " << nPuVtx << std::endl;

    if (debug_){
      std::cout << "... Size of hit vectors: sim = " <<  (*simhitvec).size() << ", reco = " << (*rechitvec).size()<< " gen " << (*genvec).size() << std::endl;
    }

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    //////// output files to save position for chi2 fit //////////
    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    std::vector<ROOT::Math::XYPoint> truthPos;
    truthPos.resize(nLayers_,ROOT::Math::XYPoint(0,0));
    bool found = getTruthPosition(genvec,truthPos,G4TrackID);
    if (!found) {
      nConvertedPhotons++;
      continue;
    }
    
    p_yvsx_truth->Fill(truthPos[10].X(),truthPos[10].Y());


    //get position of maximum to start with

    //from global max -- not working in PU
    /*
    double phimax = 0;
    double etamax = 0;
    ROOT::Math::XYZVector truthPos0(truthPos[0].X(),truthPos[0].Y(),avgZ_[0]);
    bool oneresult = getGlobalMaximum(ievt,nPuVtx,rechitvec,truthPos0,phimax,etamax);
    if (!oneresult) nMultipleMax++;
    */
    //from clusters -- using Lindsey's clustering
    HGCSSClusterVec lClusVec;
    unsigned nClusters = getClusters(rechitvec,lClusVec);

    p_nClusters->Fill(nClusters);

    if (nClusters == 0) {
      std::cout << " -- Error, no cluster found !" << std::endl;
      exit(1);
      //continue;
    }
    else if (debug_) std::cout << " -- Found " << nClusters << " clusters." << std::endl;

    //loop over clusters to find closest to the truth
    //get eta-phi of the seed ...
    ROOT::Math::XYZVector truthPos0(truthPos[0].X(),truthPos[0].Y(),avgZ_[0]);
    double etatruth = truthPos0.eta();
    double phitruth = truthPos0.phi();
    double dRmin = 10;
    unsigned clusIdx = 0;
    for (unsigned iClus(0); iClus<nClusters;++iClus){
      const HGCSSCluster & lCluster = lClusVec[iClus];
      double leta = lCluster.eta();
      double lphi = lCluster.phi();
      p_clusnHits_all->Fill(lCluster.nRecHits());
      if (lCluster.energy()>0) p_seedEoverE_all->Fill(lCluster.getSeedE()/lCluster.energy());
      p_clusLayer_all->Fill(lCluster.layer());
      p_clusWidth_all->Fill(lCluster.width());
      p_etavsphi->Fill(lphi,leta,lCluster.energy());
      double deta = leta-etatruth;
      double dphi = DeltaPhi(lphi,phitruth);
      double dR = sqrt(pow(deta,2)+pow(dphi,2));
      if (dR<dRmin){
	dRmin = dR;
	clusIdx = iClus;
      }
      double detas = leta-lCluster.getSeedEta();
      double dphis = DeltaPhi(lphi,lCluster.getSeedPhi());
      p_seeddeta_all->Fill(detas);
      p_seeddphi_all->Fill(dphis);
      if (debug_) {
	std::cout << " Cluster " << iClus << ":" ;
	lCluster.Print(std::cout);
      }
    }

    const HGCSSCluster & lCluster = lClusVec[clusIdx];
    double phimax = lCluster.phi();//getSeedPhi();
    double etamax = lCluster.eta();//getSeedEta();
    
    if (debug_){
      std::cout << " -- evt " << ievt << " found cluster mindR= " << dRmin;
      lCluster.Print(std::cout);
      std::cout << " Truth pos eta-phi = " << truthPos0.eta() << " " << truthPos0.phi() << std::endl;
    }
    
    p_mindRtruth->Fill(dRmin);

    if (dRmin > 0.1) {
      nTooFar++;
      continue;
    }

    //fill selected cluster histos
    p_clusnHits_sel->Fill(lCluster.nRecHits());
    if (lCluster.energy()>0) p_seedEoverE_sel->Fill(lCluster.getSeedE()/lCluster.energy());
    p_clusLayer_sel->Fill(lCluster.layer());
    p_clusWidth_sel->Fill(lCluster.width());
    double detas = lCluster.eta()-lCluster.getSeedEta();
    double dphis = DeltaPhi(lCluster.phi(),lCluster.getSeedPhi());
    p_seeddeta_sel->Fill(detas);
    p_seeddphi_sel->Fill(dphis);
    p_etavsphi_max->Fill(phimax,etamax);

    std::vector<double> xmax;
    xmax.resize(nLayers_,0);
    std::vector<double> ymax;
    ymax.resize(nLayers_,0);
    getMaximumCell(rechitvec,phimax,etamax,xmax,ymax);
    p_yvsx_max->Fill(xmax[10],ymax[10]);

    //get PU contrib from elsewhere in the event
    //loop over phi with same etamax
    //take average per layer: not all 9 cells of 3*3 area have hits...
    std::vector<double> puE;
    puE.resize(nLayers_,0);
    if (nPuVtx>0){
      unsigned nRandomCones = 50;
      double phistep = TMath::Pi()/nRandomCones;
      if (debug_) std::cout << "--- etamax = " << etamax << " phimax=" << phimax << " phistep = " << phistep << std::endl;
      for (unsigned ipm(0);ipm<nRandomCones;++ipm){
	std::vector<double> xmaxrc;
	xmaxrc.resize(nLayers_,0);
	std::vector<double> ymaxrc;
	ymaxrc.resize(nLayers_,0);
	double phirc = phimax-TMath::Pi();
	if (phirc < -1.*TMath::Pi()) phirc+=2.*TMath::Pi();
	if (phirc > TMath::Pi()) phirc-=2.*TMath::Pi();
	if (ipm%2==0) phirc += ipm/2*phistep+phistep/2.;
	else  phirc = phirc - ipm/2*phistep-phistep/2.;
	if (phirc < -1.*TMath::Pi()) phirc+=2.*TMath::Pi();
	if (phirc > TMath::Pi()) phirc-=2.*TMath::Pi();
	//take from geom to not be biased by hit having PU, because
	//not from geom means find cell with a hit closest to maxpos...
	getMaximumCellFromGeom(phirc,etamax,xmaxrc,ymaxrc);
	if (debug_>1) std::cout << "rc #" << ipm << " phirc=" << phirc << " xmax[10]=" << xmaxrc[10] << " ymax[10]=" << ymaxrc[10] << " r=" << sqrt(pow(xmaxrc[10],2)+pow(ymaxrc[10],2)) << std::endl;
	getPuContribution(rechitvec,xmaxrc,ymaxrc,puE);
      }
      
      //normalise to one cell: must count cells with 0 hit !
      //use cell size at etamax...
      for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
	if (debug_) std::cout << "layer " << iL ;
	unsigned nCells = nRandomCones*nSR_*geomConv_.cellSize()/geomConv_.cellSize(iL,etamax);
	puE[iL] = puE[iL]/(nCells*nCells);
	if (debug_) std::cout << " Epu=" << puE[iL] << std::endl;	
      }

    }//if PU

    std::vector<ROOT::Math::XYPoint> recoPos;
    std::vector<double> recoE;
    recoPos.resize(nLayers_,ROOT::Math::XYPoint(0,0));
    recoE.resize(nLayers_,0);
    std::vector<unsigned> nHits;
    nHits.resize(nLayers_,0);

    //get energy-weighted position and energy around maximum
    getEnergyWeightedPosition(rechitvec,nPuVtx,xmax,ymax,recoPos,recoE,nHits,puE);

    std::ofstream fout;
    std::ostringstream foutname;
    foutname << outFolder_ << "/initialPos_evt" << ievt << ".dat";
    fout.open(foutname.str());
    if (!fout.is_open()){
      std::cout << " Cannot open outfile " << foutname.str() << " for writing ! Exiting..." << std::endl;
      exit(1);
    }

    if (debug_) std::cout << " Summary of reco and truth positions:" << std::endl;
    for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
      if (debug_) std::cout << iL << " nHits=" << nHits[iL] << " Max=(" << xmax[iL] << "," << ymax[iL] << ")\t Reco=(" << recoPos[iL].X() << "," << recoPos[iL].Y() << ")\t Truth=(" << truthPos[iL].X() << "," << truthPos[iL].Y() << ")" << std::endl;
      fout << iL << " " << recoPos[iL].X() << " " << recoPos[iL].Y() << " " << truthPos[iL].X() << " " << truthPos[iL].Y() ;
      if (!doMatrix_) fout << " " << recoE[iL];
      fout << std::endl;
    }

    fout.close();

    if (doMatrix_) fillErrorMatrix(recoPos,truthPos,nHits);

    //firstEvent = false;
  }//loop on entries

  std::cout << " -- Number of converted photons: " << nConvertedPhotons << std::endl;
  std::cout << " -- Number of events with closest cluster away from truth within dR " << maxdR_ << " : " << nTooFar << std::endl;

}

void PositionFit::getMaximumCellFromGeom(const double & phimax,const double & etamax,std::vector<double> & xmax,std::vector<double> & ymax){

  for (unsigned iL(0); iL<nLayers_;++iL){
    double theta = 2*atan(exp(-etamax));
    double rho = avgZ_[iL]/cos(theta);
    xmax[iL] = rho*tan(theta)*cos(phimax);
    ymax[iL] = rho*tan(theta)*sin(phimax);
  }//loop on layers

}

void PositionFit::getMaximumCell(std::vector<HGCSSRecoHit> *rechitvec,const double & phimax,const double & etamax,std::vector<double> & xmax,std::vector<double> & ymax){
  
  std::vector<double> dRmin;
  dRmin.resize(nLayers_,10);
  
  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
    const HGCSSRecoHit & lHit = (*rechitvec)[iH];
    
    unsigned layer = lHit.layer();
    if (layer >= nLayers_) {
      continue;
    }
    
    double posx = lHit.get_x();
    if (fixForPuMixBug_) posx-=1.25;
    double posy = lHit.get_y();
    if (fixForPuMixBug_) posy-=1.25;
    double posz = lHit.get_z();
    if (debug_>1) {
      std::cout << " --  RecoHit " << iH << "/" << (*rechitvec).size() << " --" << std::endl
		<< " --  position x,y " << posx << "," << posy << std::endl;
      lHit.Print(std::cout);
    }


    ROOT::Math::XYZVector pos(posx,posy,posz);
    double deta = fabs(pos.eta()-etamax);
    double dphi = DeltaPhi(pos.phi(),phimax);
    
    double dR = sqrt(pow(deta,2)+pow(dphi,2));
    if (dR<dRmin[layer]) {
      dRmin[layer] = dR;
      xmax[layer] = posx;
      ymax[layer] = posy;
    }
    
    double energy = lHit.energy();
    p_recoxy[layer]->Fill(posx,posy,energy);
    
  }//loop on rechits
    
}

void PositionFit::getEnergyWeightedPosition(std::vector<HGCSSRecoHit> *rechitvec,const unsigned nPU, const std::vector<double> & xmax,const std::vector<double> & ymax,std::vector<ROOT::Math::XYPoint> & recoPos,std::vector<double> & recoE,std::vector<unsigned> & nHits,std::vector<double> & puE,const bool puSubtracted){
  
  double eSum_puMean = 0;
  double eSum_puEvt = 0;

  double step = geomConv_.cellSize()*nSR_/2.+0.1;//+0.1 to accomodate double precision
  if (debug_) std::cout << "step = " << step << std::endl;

  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
    const HGCSSRecoHit & lHit = (*rechitvec)[iH];
    double energy = lHit.energy();//in MIP already...
    unsigned layer = lHit.layer();
    double posx = lHit.get_x();
    if (fixForPuMixBug_) posx-=1.25;
    double posy = lHit.get_y();
    if (fixForPuMixBug_) posy-=1.25;


    if (fabs(posx-xmax[layer]) < step && 
	fabs(posy-ymax[layer]) < step){
      if (puSubtracted) {
	double leta = lHit.eta();
	if (debug_>1) std::cout << " -- Hit " << iH << ", eta=" << leta << ", energy before PU subtraction: " << energy << " after: " ;
	double lCorMean =  puDensity_.getDensity(leta,layer,geomConv_.cellSizeInCm(layer,leta),nPU);
	eSum_puMean += lCorMean;
	p_hitMeanPuContrib->Fill(layer,lCorMean);
	double lCorEvent = puE[layer];
	eSum_puEvt += lCorEvent;
	p_hitEventPuContrib->Fill(layer,lCorEvent);
	double lCor = 0;
	if (useMeanPU_) lCor = lCorMean;
	else lCor = lCorEvent;
	energy = std::max(0.,energy - lCor);
	if (debug_>1) std::cout << energy << std::endl;
      }
      recoPos[layer].SetX(recoPos[layer].X() + posx*energy);
      recoPos[layer].SetY(recoPos[layer].Y() + posy*energy);
      recoE[layer] += energy;
      if (energy>0) nHits[layer]++;
    }
    
  }//loop on rechits

  p_diffPuContrib->Fill(eSum_puEvt-eSum_puMean);

  for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
    if (nHits[iL]==0) continue;
    recoPos[iL].SetX(recoPos[iL].X()/recoE[iL]);
    recoPos[iL].SetY(recoPos[iL].Y()/recoE[iL]);
  }
  
}

void PositionFit::getPuContribution(std::vector<HGCSSRecoHit> *rechitvec, const std::vector<double> & xmax,const std::vector<double> & ymax,std::vector<double> & puE){

  double step = geomConv_.cellSize()*nSR_/2.+0.1;

  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on rechits
    const HGCSSRecoHit & lHit = (*rechitvec)[iH];
    double energy = lHit.energy();//in MIP already...
    unsigned layer = lHit.layer();
    double posx = lHit.get_x();
    if (fixForPuMixBug_) posx-=1.25;
    double posy = lHit.get_y();
    if (fixForPuMixBug_) posy-=1.25;

    //std::cout << "- iH " << iH << " x=" << posx << " xmax=" << xmax[layer] << " y=" << posy << " ymax=" << ymax[layer] << " step " << step << std::endl;
    if (fabs(posx-xmax[layer]) < step && 
	fabs(posy-ymax[layer]) < step){
      if (debug_>1) std::cout << "- iH " << iH 
			     << " x=" << posx 
			     << " xmax=" << xmax[layer] 
			     << " y=" << posy 
			     << " ymax=" << ymax[layer] 
			     << " step " << step 
			     << " --- Pass, layer" << layer 
			     << " E="<< energy << std::endl;
      puE[layer] += energy;
    }

  }//loop on rechits
  //exit(1);
}

void PositionFit::fillErrorMatrix(const std::vector<ROOT::Math::XYPoint> & recoPos,const std::vector<ROOT::Math::XYPoint> & truthPos,const std::vector<unsigned> & nHits){

  for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
    if (nHits[iL]==0) continue;
    double residual_xi = recoPos[iL].X()-truthPos[iL].X();
    double residual_yi = recoPos[iL].Y()-truthPos[iL].Y();
    p_residuals_x->Fill(residual_xi);
    p_residuals_y->Fill(residual_yi);
    if (fabs(residual_xi)>residualMax_ || fabs(residual_yi)>residualMax_) continue;
    mean_[0][iL] += residual_xi;
    mean_[1][iL] += residual_yi;
    ++nL_mean_[iL];
    for (unsigned jL(0);jL<nLayers_;++jL){//loop on layers
      if (nHits[jL]==0) continue;
      double residual_xj = recoPos[jL].X()-truthPos[jL].X();
      double residual_yj = recoPos[jL].Y()-truthPos[jL].Y();
      if (fabs(residual_xj)>residualMax_ || fabs(residual_yj)>residualMax_) continue;
      double sigma_x = residual_xi*residual_xj;
      double sigma_y = residual_yi*residual_yj;
      sigma_[0][iL][jL] += sigma_x;
      sigma_[1][iL][jL] += sigma_y;
      ++nL_sigma_[iL][jL];
    }//loop on layers
  }//loop on layers

}


void PositionFit::finaliseErrorMatrix(){
  //finalise error matrix

  std::ofstream fmatrix;
  std::ostringstream fmatrixname;
  fmatrixname << matrixFolder_ << "/errorMatrix.dat";
  fmatrix.open(fmatrixname.str());
  if (!fmatrix.is_open()){
    std::cout << " Cannot open outfile " << fmatrixname.str() << " for writing ! Exiting..." << std::endl;
    exit(1);
  }

  matrix_.ResizeTo(nLayers_,nLayers_);

  //set mean values first
  for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
    mean_[0][iL] = mean_[0][iL]/nL_mean_[iL];
    mean_[1][iL] = mean_[1][iL]/nL_mean_[iL];
  }
  //set sigmas and fill matrix
  for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
    for (unsigned jL(0);jL<nLayers_;++jL){//loop on layers
      sigma_[0][iL][jL] = sigma_[0][iL][jL]/nL_sigma_[iL][jL];
      sigma_[1][iL][jL] = sigma_[1][iL][jL]/nL_sigma_[iL][jL];
      //consider average of both x and y in one matrix
      matrix_[iL][jL] = 0.5*(sigma_[0][iL][jL]-mean_[0][iL]*mean_[0][jL]+
			     sigma_[1][iL][jL]-mean_[1][iL]*mean_[1][jL]);
      //matrix_[jL][iL] = matrix_[iL][jL];
      
      fmatrix << iL << " " << jL << " " << std::setprecision(15) << matrix_[iL][jL] << std::endl;

      //if (iL!=jL){
      //p_matrix->Fill(jL,iL,matrix_[iL][jL]);
      //}
    }
  }

  fmatrix.close();

}

void PositionFit::fillCorrelationMatrix(){
  std::cout << " -- Filling correlation matrix" << std::endl;
  outputFile_->cd(outputDir_.c_str());

  p_errorMatrix = new TH2D("p_errorMatrix",";i;j;M_{ij}",
			   nLayers_,0,nLayers_,
			   nLayers_,0,nLayers_);
  p_corrMatrix = new TH2D("p_corrMatrix",";i;j;M_{ij}",
			  nLayers_,0,nLayers_,
			  nLayers_,0,nLayers_);

  corrMatrix_.ResizeTo(nLayers_,nLayers_);
  for (unsigned iL(0);iL<nLayers_;++iL){//loop on layers
    for (unsigned jL(0);jL<nLayers_;++jL){//loop on layers
      p_errorMatrix->Fill(iL,jL,matrix_[iL][jL]);
      if (matrix_[iL][iL]!=0 && matrix_[jL][jL]!= 0){
	corrMatrix_[iL][jL] =matrix_[iL][jL]/sqrt(matrix_[iL][iL]*matrix_[jL][jL]); 
	p_corrMatrix->Fill(iL,jL,corrMatrix_[iL][jL]);
      }
    }
  }
}


bool PositionFit::fillMatrixFromFile(){

  std::ifstream fmatrix;
  std::ostringstream fmatrixname;
  fmatrixname << matrixFolder_ << "/errorMatrix.dat";
  fmatrix.open(fmatrixname.str());
  if (!fmatrix.is_open()){
    std::cout << " -- Cannot open outfile " << fmatrixname.str() << "! Refilling the matrix..." << std::endl;
    return false;
  }

  matrix_.ResizeTo(nLayers_,nLayers_);
  if (debug_) std::cout << " -- Error matrix: " << std::endl;
  while (!fmatrix.eof()){
    unsigned iL=nLayers_;
    unsigned jL=nLayers_;
    double m=0;
    fmatrix>>iL>>jL>>m;
    if (iL<nLayers_ && jL<nLayers_){
      if (debug_) std::cout << std::setprecision(15) << iL << " " << jL << " " << m << std::endl;
      matrix_[iL][jL] = m;
    }
  }
  
  fmatrix.close();

  return true;
}

bool PositionFit::initialiseLeastSquareFit(){
  initialiseFitHistograms();

  //try reading matrix from file, if fail
  //return false and refill the matrix.
  if (doMatrix_ && !fillMatrixFromFile()) return false;
  
  //fill matrices
  if (doMatrix_) fillCorrelationMatrix();

  //get back data for each event and perform chi2 fit:
  std::cout << " -- Performing chi2 fit for each event" << std::endl;

  nInvalidFits_=0;
  nFailedFitsAfterCut_=0;

  //open new file to save accurate positions
  std::ostringstream foutname;
  foutname << outFolder_ << "/accuratePos.dat";
  fout_.open(foutname.str());
  if (!fout_.is_open()){
    std::cout << " Cannot open outfile " << foutname.str() << " for writing ! Exiting..." << std::endl;
    exit(1);
  }

  return true;

}

unsigned PositionFit::performLeastSquareFit(const unsigned ievt,
					    FitResult & fit){

    //cut outliers
  unsigned fitres = fitEvent(ievt,fit,true);
    if (fitres==1){
      // std::cout << " -- Number of genparticles: " << (*genvec).size() << std::endl;
      // for (unsigned iP(0); iP<(*genvec).size(); ++iP){//loop on gen particles 
      // 	std::cout << " --- particle " << iP << std::endl;
      // 	(*genvec)[iP].Print(std::cout);
      // }
      std::cout << " -- Event " << ievt << " skipped." << std::endl;
    }
    if (fitres>1) {
      //std::cout << " ---- Fit failed ! Reprocess event " << ievt << " cutting outliers" << std::endl;
      //nFailedFits++;
      //if (fitEvent(ievt,fit,true)>1)
      std::cout << " -- Event " << ievt << " failed fit." << std::endl;
      nFailedFitsAfterCut_++;
    }
  
    return fitres;
}

void PositionFit::finaliseFit(){
  fout_.close();    
  outputFile_->Flush();
  std::cout << " -- Number of invalid fits: " << nInvalidFits_ << std::endl;
  std::cout << " -- Number of fits failed after cutting outliers: " << nFailedFitsAfterCut_ << std::endl;

}

unsigned PositionFit::fitEvent(const unsigned ievt, 
			       FitResult & fit,
			       const bool cutOutliers){

  std::vector<unsigned> layerId;
  std::vector<double> posx;
  std::vector<double> posy;
  std::vector<double> posz;
  std::vector<double> posxtruth;
  std::vector<double> posytruth;
  std::vector<double> Ereco;
  layerId.reserve(nLayers_);
  posx.reserve(nLayers_);
  posy.reserve(nLayers_);
  posz.reserve(nLayers_);
  posxtruth.reserve(nLayers_);
  posytruth.reserve(nLayers_);
  Ereco.reserve(nLayers_);
  if (!getPositionFromFile(ievt,
			   layerId,posx,posy,posz,
			   posxtruth,posytruth,
			   Ereco,
			   cutOutliers)){
    nInvalidFits_++;
    return 1;
  }

  const unsigned nL = layerId.size();
  
  
  //fill some control histograms
  //only once
  if (cutOutliers){
    for (unsigned iL(0); iL<nL;++iL){
      //std::cout << layerId[iL] << " " << posx[iL] << " " << posy[iL] << " " << posz[iL] << std::endl;
      p_recoXvsLayer->Fill(layerId[iL],posx[iL]);
      p_recoYvsLayer->Fill(layerId[iL],posy[iL]);
      p_recoZvsLayer->Fill(layerId[iL],posz[iL]);
      p_truthXvsLayer->Fill(layerId[iL],posxtruth[iL]);
      p_truthYvsLayer->Fill(layerId[iL],posytruth[iL]);
    }
  }

  p_nLayersFit->Fill(nL);

  //if less than 3 valid layers: no point doing a fit !!
  if (nL<3){
    nInvalidFits_++;
    return 1;
  }
  
  //number of points: x and y per layer minus number of parameters: 2 for x + 2 for y.
  double ndf = 2*nL-4;
  
  //Get error matrix removing lines with zero hits
  TMatrixDSym e(nL);
  TVectorD u(nL),z(nL),x(nL),y(nL);
  
  for(unsigned i(0);i<nL;++i) {
    u(i)=1.0;
    z(i)=posz[i];
    //std::cout << "fit() z(" << i << ") = " << z(i) << std::endl;
    
    for(unsigned j(i);j<nL;++j) {
      e(i,j)=matrix_(layerId[i],layerId[j]);
      e(j,i)=matrix_(layerId[j],layerId[i]);
    }
  }
  
  e.Invert();
  
  //do fit for reco and truth
  double positionFF[2][2];
  double TanAngle[2][2];
 
  for (unsigned rt(0); rt<2;++rt){
    if (debug_) {
      std::cout << "... Processing ";
      if (rt==0) std::cout << " fit to reco position.";
      else std::cout << " fit to truth position.";
      std::cout << std::endl;
    }
    double chiSq(0.0);
    double position[2];
    
    TMatrixD fitMatrix(4,4);
    for (unsigned ii(0);ii<4;++ii){
      for (unsigned ij(0);ij<4;++ij){
	fitMatrix[ii][ij]=0;
      }
    }

    //resolve equation for x and y separately
    for(unsigned xy(0);xy<2;xy++) {//loop on x or y
      if (debug_) {
	std::cout << "... Processing ";
	if (xy==0) std::cout << " fit to x position.";
	else std::cout << " fit to y position.";
	std::cout << std::endl;
      }
      for(unsigned i(0);i<nL;i++) {
	x(i)= rt==0 ? ((xy==0) ? posx[i] : posy[i]) : ((xy==0) ? posxtruth[i] : posytruth[i]);
	//std::cout << "fit() x(" << i << ") = " << x(i) << std::endl;
      }
      
      TMatrixD w(2,2);
      TVectorD v(2),p(2);
      
      w(0,0)=u*(e*u);
      w(0,1)=u*(e*z);
      w(1,0)=z*(e*u);
      w(1,1)=z*(e*z);
      
      v(0)=u*(e*x);
      v(1)=z*(e*x);
      
      w.Invert();
      
      p=w*v;
      if (debug_) {
	std::cout << "fit() w(0,0) = " << w(0,0) << std::endl;
	std::cout << "fit() w(0,1) = " << w(0,1) << std::endl;
	std::cout << "fit() w(1,0) = " << w(1,0) << std::endl;
	std::cout << "fit() w(1,1) = " << w(1,1) << std::endl;	
	std::cout << "fit() p(0) = " << p(0) << std::endl;
	std::cout << "fit() p(1) = " << p(1) << std::endl;
      }
      
      position[xy] = p(0);
      positionFF[rt][xy] = p(0)+p(1)*posz[0];
      TanAngle[rt][xy] = p(1);

      //sanity check for nan values
      if (w(0,0)==w(0,0)) fitMatrix[2*xy][2*xy]=w(0,0);
      if (w(0,1)==w(0,1)) fitMatrix[2*xy][2*xy+1]=w(0,1);
      if (w(1,0)==w(1,0)) fitMatrix[2*xy+1][2*xy]=w(1,0);
      if (w(1,1)==w(1,1)) fitMatrix[2*xy+1][2*xy+1]=w(1,1);
      
      
      TVectorD dp(nL);
      for(unsigned i(0);i<nL;i++) {
	dp(i)=x(i)-p(0)-p(1)*z(i);
      }
      
      chiSq+=dp*(e*dp);
    }//loop on x or y
    
    //chi2 test
    if (chiSq/ndf>chi2ndfmax_) {
      std::cout << " ---- Fit failed for event " << ievt << std::endl;
      std::cout << "Chi2/ndf = " << chiSq << "/" << ndf << "=" << chiSq/ndf << std::endl;
      //std::cout << "fitw(0,0) = " << fitMatrix[0][0] << std::endl;
      //std::cout << "fitw(1,1) = " << fitMatrix[1][1] << std::endl;
      //std::cout << "fitw(2,2) = " << fitMatrix[2][2] << std::endl;
      //std::cout << "fitw(3,3) = " << fitMatrix[3][3] << std::endl;	
      std::cout << "ecal frontface position = " << positionFF[0][0] << " " << positionFF[0][1] << std::endl;
      std::cout << "ecal frontface tanAngle = " << TanAngle[0][0] << " " << TanAngle[0][1] << std::endl;
      return 2;
    }

    p_chi2[rt]->Fill(chiSq);
    p_chi2overNDF[rt]->Fill(chiSq/ndf);
    p_impactX[rt]->Fill(positionFF[rt][0]);
    p_tanAngleX[rt]->Fill(TanAngle[rt][0]);
    p_impactY[rt]->Fill(positionFF[rt][1]);
    p_tanAngleY[rt]->Fill(TanAngle[rt][1]);

    if (rt==0) {
      p_positionReso->Fill(sqrt(fitMatrix[0][0]));
      p_angularReso->Fill(sqrt(fitMatrix[1][1]));
      
      fout_ << ievt << " " 
	   << position[0] << " " 
	   << sqrt(fitMatrix[0][0]) << " " 
	   << TanAngle[0][0] << " " 
	   << sqrt(fitMatrix[1][1]) << " "
	   << position[1] << " " 
	   << sqrt(fitMatrix[2][2]) << " "
	   << TanAngle[0][1] << " "
	   << sqrt(fitMatrix[3][3])
	   << std::endl;
    

      for (unsigned iL(0); iL<nL;++iL){
	double x = position[0]+TanAngle[0][0]*posz[iL];
	double y = position[1]+TanAngle[0][1]*posz[iL];
	p_fitXvsLayer->Fill(layerId[iL],x);
	p_fitYvsLayer->Fill(layerId[iL],y);
	//eventPos[layerId[iL]] = ROOT::Math::XYZVector(x,y,posz[iL]);
      }
      fit.found = true;
      fit.pos_x = position[0];
      fit.tanangle_x = TanAngle[0][0];
      fit.pos_y = position[1];
      fit.tanangle_y = TanAngle[0][1];
    }

  }//reco or truth

  recoDir_ = Direction(TanAngle[0][0],TanAngle[0][1]);

  p_impactX_residual->Fill(positionFF[0][0]-positionFF[1][0]);
  p_tanAngleX_residual->Fill(TanAngle[0][0]-TanAngle[1][0]);
  p_impactY_residual->Fill(positionFF[0][1]-positionFF[1][1]);
  p_tanAngleY_residual->Fill(TanAngle[0][1]-TanAngle[1][1]);
      
  //std::cout << " -- Size of eventPos=" << eventPos.size() << std::endl;
  return 0;
}

bool PositionFit::getPositionFromFile(const unsigned ievt,
				      std::vector<unsigned> & layerId,
				      std::vector<double> & posx,
				      std::vector<double> & posy,
				      std::vector<double> & posz,
				      std::vector<double> & posxtruth,
				      std::vector<double> & posytruth,
				      std::vector<double> & E,
				      bool cutOutliers,
				      bool print){


  std::ifstream fin;
  std::ostringstream finname;
  finname << outFolder_ << "/initialPos_evt" << ievt << ".dat";
  fin.open(finname.str());
  if (!fin.is_open()){
    if (print) std::cout << " Cannot open input file " << finname.str() << "!" << std::endl;
      return false;
  }
  
  while (!fin.eof()){
    unsigned l=nLayers_;
    double xr=0,yr=0,xt=0,yt=0,e=0;
    fin>>l>>xr>>yr>>xt>>yt;
    if (!doMatrix_) fin>>e;
    if (l<nLayers_){
      bool pass=fabs(xr-xt)<residualMax_ && fabs(yr-yt)<residualMax_;
      if (!cutOutliers || (cutOutliers && pass)){
	layerId.push_back(l);
	posx.push_back(xr);
	posy.push_back(yr);
	posz.push_back(avgZ_[l]);
	posxtruth.push_back(xt);
	posytruth.push_back(yt);
      }
      //use all for energy estimate
      if (!doMatrix_) E.push_back(e);
    }
  }
  
  fin.close();
  /*
  //@TODO to use something else than truth info :/
  if (cutOutliers){
    //
    for (unsigned i(0);i<layerId.size();++i){
      double xt=?;
      double yt=?;
      double xr=posx[i];
      double yr=posy[i];
      bool pass=fabs(xr-xt)<residualMax_ && fabs(yr-yt)<residualMax_;
      std::cout << i << " l=" << layerId[i] << " pass=" << pass << " xr=" << xr << " yr=" << yr << std::endl;
      if (!pass) {
	std::cout << " ---- erase layer " << *(layerId.begin()+i) << std::endl;
	layerId.erase(layerId.begin()+i);
	posx.erase(posx.begin()+i);
	posy.erase(posy.begin()+i);
	posxtruth.erase(posxtruth.begin()+i);
	posytruth.erase(posytruth.begin()+i);
	i--;
      }
    }
  }
  */

  return true;
}




