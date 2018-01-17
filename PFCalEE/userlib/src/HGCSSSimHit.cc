#include "HGCSSSimHit.hh"

#include <iomanip>
#include <cmath>
#include <stdlib.h>

HGCSSSimHit::HGCSSSimHit(const G4SiHit & aSiHit, 
			 const unsigned & asilayer, 
			 TH2Poly* map, 
			 const float ,
			 const bool etaphimap){
  energy_ = aSiHit.energy;
  //energy weighted time
  //PS: need to call calculateTime() after all hits 
  //have been added to have divided by totalE!!
  time_ = aSiHit.time*aSiHit.energy;
  zpos_ = aSiHit.hit_z;
  setLayer(aSiHit.layer,asilayer);

  //coordinates in mm
  //double z = aSiHit.hit_x;
  double x = aSiHit.hit_x;
  double y = aSiHit.hit_y;
  double eee = aSiHit.energy;
  //cellid encoding:
  //map->Reset("");
  //map->Fill(x,y);
  //GetMaximumBin doesn't work :(

  assert(map);
  if (etaphimap){
    ROOT::Math::XYZPoint pos = ROOT::Math::XYZPoint(x,y,zpos_);
    cellid_ = map->FindBin(pos.eta(),pos.phi());
  }
  else {
    cellid_ = map->FindBin(x,y);
      double dbgphi=atan2(y,x);
      double dbgr=sqrt(x*x+y*y);
      double dbgtantheta=dbgr/zpos_;
      double dbgtheta=atan(dbgtantheta);
      double dbgeta=-log(tan(dbgtheta/2.));

    if(cellid_>1000000) {
      std::cout<<"danger danger id energy x y r z phi eta "<<cellid_<<" "<<eee<<" "<<x<<","<<y<<" "<<dbgr<<" "<<zpos_<<" "<<dbgphi<<" "<<dbgeta<<std::endl;
    } else {
      std::cout<<"good          id energy x y r z phi eta "<<cellid_<<" "<<eee<<" "<<x<<","<<y<<" "<<dbgr<<" "<<zpos_<<" "<<dbgphi<<" "<<dbgeta<<std::endl;

    }
  }
  
  //for (int ix(1);ix<map->GetNumberOfBins()+1; ++ix){
  //if (map->GetBinContent(ix)!=0)
    //cellid_ = ix;
    //std::cout << ix << " " << map->GetBinContent(ix) << std::endl;
  //}

  /*
  TIter next(map->GetBins());
  TObject *obj=0; 
  TH2PolyBin *polyBin = 0;

  while ((obj=next())){
    polyBin=(TH2PolyBin*)obj;
    int id = polyBin->GetBinNumber();
    if (id==cellid_) break; 
  }

  std::cout << " - Sanity check: x,y = " << x << " " << y 
	    << " cellid=" << cellid_ 
	    << " polybin# " << polyBin->GetBinNumber() << " area " << polyBin->GetArea() << std::endl
	    << " x bin = " << polyBin->GetXMin() << "-" << polyBin->GetXMax() << std::endl
	    << " y bin = " << polyBin->GetYMin() << "-" << polyBin->GetYMax() << std::endl
	    << " middle = " << (polyBin->GetXMax()+polyBin->GetXMin())/2 
	    << " " << (polyBin->GetYMax()+polyBin->GetYMin())/2 
	    << std::endl;
*/


  //bool x_side = x>0 ? true : false;
  //bool y_side = y>0 ? true : false;
  //unsigned x_cell = static_cast<unsigned>(fabs(x)/(cellSize*getGranularity()));
  //unsigned y_cell = static_cast<unsigned>(fabs(y)/(cellSize*getGranularity()));

  //encodeCellId(x_side,y_side,x_cell,y_cell);

  nGammas_= 0;
  nElectrons_ = 0;
  nMuons_ = 0;
  nNeutrons_ = 0;
  nProtons_ = 0;
  nHadrons_ = 0;
  if(abs(aSiHit.pdgId)==22) nGammas_++;
  else if(abs(aSiHit.pdgId)==11) nElectrons_++;
  else if(abs(aSiHit.pdgId)==13) nMuons_++;
  else if(abs(aSiHit.pdgId)==2112) nNeutrons_++;
  else if(abs(aSiHit.pdgId)==2212) nProtons_++;
  else nHadrons_++;

  trackIDMainParent_ = aSiHit.parentId;
  energyMainParent_ = aSiHit.energy;

}

/*void HGCSSSimHit::encodeCellId(const bool x_side,const bool y_side,const unsigned x_cell,const unsigned y_cell){
  cellid_ = 
    x_side | (x_cell<<1) |
    (y_side<<16) | (y_cell<<17);

  // std::cout << " Cross-check of encoding: cellid=" << cellid_ << std::endl
  // 	    << " x_side " << x_side << " " << get_x_side() << std::endl
  // 	    << " y_side " << y_side << " " << get_y_side() << std::endl
  // 	    << " x_cell " << x_cell << " " << get_x_cell() << std::endl
  // 	    << " y_cell " << y_cell << " " << get_y_cell() << std::endl
  //   ;
  }*/

void HGCSSSimHit::Add(const G4SiHit & aSiHit){

  time_ = time_ + aSiHit.time*aSiHit.energy;
  //PS: need to call calculateTime() after all hits 
  //have been added to have divided by totalE!!

  if(abs(aSiHit.pdgId)==22) nGammas_++;
  else if(abs(aSiHit.pdgId)==11) nElectrons_++;
  else if(abs(aSiHit.pdgId)==13) nMuons_++;
  else if(abs(aSiHit.pdgId)==2112) nNeutrons_++;
  else if(abs(aSiHit.pdgId)==2212) nProtons_++;
  else nHadrons_++;

  energy_ += aSiHit.energy;
  if (aSiHit.energy > energyMainParent_){
    trackIDMainParent_ = aSiHit.parentId;
    energyMainParent_ = aSiHit.energy;
  }

}

/*double HGCSSSimHit::eta() const {
  double x = get_x();
  double y = get_y();
  double theta = acos(fabs(zpos_)/sqrt(zpos_*zpos_+x*x+y*y));
  double leta = -log(tan(theta/2.));
  if (zpos_>0) return leta;
  else return -leta;
  }*/



std::pair<double,double> HGCSSSimHit::get_xy(const HGCSSSubDetector & subdet,
					     const HGCSSGeometryConversion & aGeom,
					     const unsigned shape) const {
  if (subdet.isScint){
    std::pair<double,double> etaphi = subdet.type==DetectorEnum::BHCAL1?aGeom.squareGeom1.find(cellid_)->second : aGeom.squareGeom2.find(cellid_)->second;
    //convert back to x-y
    double theta = 2*atan(exp(-1.*etaphi.first));
    double r = zpos_/cos(theta);
    double x = r*sin(theta)*cos(etaphi.second);
    double y = r*sin(theta)*sin(etaphi.second);
    return std::pair<double,double>(x,y);
  }
  else if (shape==4) return aGeom.squareGeom.find(cellid_)->second;
  else {
    if (shape==2) return aGeom.diamGeom.find(cellid_)->second;
    else if (shape==3) return aGeom.triangleGeom.find(cellid_)->second;
    else return aGeom.hexaGeom.find(cellid_)->second;
  }

}

ROOT::Math::XYZPoint HGCSSSimHit::position(const HGCSSSubDetector & subdet,
					   const HGCSSGeometryConversion & aGeom,
					   const unsigned shape) const{
  std::pair<double,double> xy = get_xy(subdet,aGeom,shape);
  return ROOT::Math::XYZPoint(xy.first/10.,xy.second/10.,zpos_/10.);
}

double HGCSSSimHit::theta(const HGCSSSubDetector & subdet,
			  const HGCSSGeometryConversion & aGeom,
			  const unsigned shape) const {
  return 2*atan(exp(-1.*eta(subdet,aGeom,shape)));
}

double HGCSSSimHit::eta(const HGCSSSubDetector & subdet,
			const HGCSSGeometryConversion & aGeom,
			const unsigned shape) const {
  return position(subdet,aGeom,shape).eta();
}

double HGCSSSimHit::phi(const HGCSSSubDetector & subdet,
			const HGCSSGeometryConversion & aGeom,
			const unsigned shape) const {
  return position(subdet,aGeom,shape).phi();
}






void HGCSSSimHit::Print(std::ostream & aOs) const{
  aOs << "====================================" << std::endl
      << " = Layer " << layer() << " siLayer " << silayer() << " cellid " << cellid_ << std::endl
      << " = Energy " << energy_ << " time " << time_ << std::endl
      << " = g " << nGammas_ 
      << " e " << nElectrons_ 
      << " mu " << nMuons_ 
      << " neutron " << nNeutrons_ 
      << " proton " << nProtons_ 
      << " had " << nHadrons_ 
      << std::endl
      << " = main parent: trackID " << trackIDMainParent_ << " efrac " << mainParentEfrac()
      << std::endl
      << "====================================" << std::endl;

}

