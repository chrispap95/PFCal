//
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
#ifdef _DEBUG
#include "debug_new.h"
#endif
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
//#include "HGCSSSimpleHit.hh"

#include "PositionFit.hh"
#include "SignalRegion.hh"

#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"

#include "utilities.h"

using boost::lexical_cast;
namespace po=boost::program_options;

bool domap = true;


void SNAPrec_rl(TH2F* h_1,std::vector<HGCSSRecoHit> *rechitvec) {
    std::cout<<"SNAP"<<std::endl;
    double ymax = h_1->GetYaxis()->GetXmax();
    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
        HGCSSRecoHit lHit = (*rechitvec)[iH];
        double xh=lHit.get_x();
        double yh=lHit.get_y();
        unsigned ilayer=lHit.layer();
        //unsigned ixx=ilayer;
        //if(ixx>52) ixx=ixx-17;
        double rh=sqrt(xh*xh+yh*yh);
        double Eh=lHit.energy();
        h_1->Fill(ilayer+0.5,std::min(rh,ymax-200.),Eh);
    }

    return;
}

void SNAPrec_rz(TH2F* h_1,std::vector<HGCSSRecoHit> *rechitvec) {
    std::cout<<"SNAP"<<std::endl;
    double xmax = h_1->GetXaxis()->GetXmax();
    double ymax = h_1->GetYaxis()->GetXmax();
    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
        HGCSSRecoHit lHit = (*rechitvec)[iH];
        double xh=lHit.get_x();
        double yh=lHit.get_y();
        double zh=lHit.get_z();
        double rh=sqrt(xh*xh+yh*yh);
        double Eh=lHit.energy();
        h_1->Fill(std::min(zh,xmax-0.001),std::min(rh,ymax-200.),Eh);
    }

    return;
}

void SNAPsim(TH2F* h_1,std::vector<HGCSSSimHit> *simhitvec,
    HGCSSDetector & myDetector,
    HGCSSGeometryConversion & aGeom,
    unsigned shape
) {
    std::cout<<"SNAP"<<std::endl;
    double xmax = h_1->GetXaxis()->GetXmax();
    double ymax = h_1->GetYaxis()->GetXmax();
    for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
        HGCSSSimHit lHit = (*simhitvec)[iH];
        unsigned layer = lHit.layer();
        const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
        ROOT::Math::XYZPoint pp = lHit.position(subdet,aGeom,shape);
        double Eh=lHit.energy();
        //std::cout<<Eh<<" "<<pp.z()<<" "<<pp.r()<<std::endl;
        h_1->Fill(std::min(pp.z(),xmax-0.001),std::min(pp.r(),ymax-0.001),Eh);
    }

    return;
}



double DeltaR(double eta1,double phi1,double eta2,double phi2){
    double dr=99999.;
    double deta=fabs(eta1-eta2);
    double dphi=fabs(phi1-phi2);
    if(dphi>TMath::Pi()) dphi=2.*TMath::Pi()-dphi;
    dr=sqrt(deta*deta+dphi*dphi);
    return dr;
}


int main(int argc, char** argv){//main



    ///////////////////////////////////////////////////////////////
    // initialize some variables
    ///////////////////////////////////////////////////////////////


    const unsigned nscintlayer=16;
    const unsigned scintoffset=36;
    const unsigned nsilayer=52;
    //  const unsigned nDeadSiliconCell=620100;//DPG presentation slide number 3

    bool ipfirst[nscintlayer];
    for(int i(0);i<nscintlayer;i++) {ipfirst[i]=true;}

    unsigned scintminid[nscintlayer]={145,179,219,206,159,205,274,284,28,223,47,24,219,1,188,211};
    unsigned scintmaxid[nscintlayer]={32545,32579,32619,32606,20823,20869,20938,20948,20692,20887,20711,20688,20883,20665,20852,20875};

    unsigned siminid[nsilayer]={34987,34497,33997,33510,33007,33001,32505,32011,31515,31024,30531,30518,30025,29528,29035,28548,28045,28036,27542,27049,26552,26062,25562,25556,25057,24566,24066,23576,22080,20105,18607,17114,15136,13643,12150,10177,34987,34497,33507,33007,33001,32505,32011,31515,31024,30531,30518,30025,29528,29035,28548};
    unsigned simaxid[nsilayer]={250502,250992,251492,251979,2522479,252488,252982,253475,253972,254465,254955,254971,255464,255958,256454,256941,257441,257451,257947,258437,258937,259427,259927,259930,260430,260920,261420,261910,263406,265382,266882,268375,270347,271844,273337,275316,250502,250992,251979,252479,252488,252982,253475,253972,254465,254955,254971,255464,255958,256454,256941};


    //Input output and config options
    std::string cfg;
    unsigned pNevts;
    std::string outFilePath;
    std::string filePath;
    std::string digifilePath;
    unsigned nRuns;
    std::string simFileName;
    std::string recoFileName;
    unsigned debug;
    double deadfrac;
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
    ("deadfrac",    po::value<double>(&deadfrac)->default_value(0))
    ;
    po::store(po::command_line_parser(argc, argv).options(config).allow_unregistered().run(), vm);
    po::store(po::parse_config_file<char>(cfg.c_str(), config), vm);
    po::notify(vm);


    std::string inFilePath = filePath+simFileName;

    std::cout << " -- Input parameters: " << std::endl
    << " -- Input file path: " << filePath << std::endl
    << " -- Digi Input file path: " << digifilePath << std::endl
    << " -- Output file path: " << outFilePath << std::endl
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

    HGCSSCalibration mycalib(inFilePath);

    //corrected for Si-Scint overlap
    const unsigned nLayers = 52;//


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
    geomConv.initialiseSquareMap1(1.4,3.0,-1.*TMath::Pi(),TMath::Pi(),0.01745);//eta phi segmentation
    geomConv.initialiseSquareMap2(1.4,3.0,-1.*TMath::Pi(),TMath::Pi(),0.02182);//eta phi segmentation
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
    const int nlay=28;

    TH1F* h_rechitsum = new TH1F("h_rechitsum","Rechitsum Energy",250,100,400.);
    ///////////////////////////////
    // for missing channel study //
    ///////////////////////////////

    //####### SCINTILLATOR
    std::set<std::pair<unsigned, unsigned>> deadlist;
    unsigned nchan=0;
    for(int i(0);i<nscintlayer;i++) {
        nchan+=(scintmaxid[i]-scintminid[i]);
    }
    std::cout<<"total number scintillator channels is "<<nchan<<std::endl;

    //***** Define Scintillator randomly 1% dead channels parameters
    unsigned ndead=deadfrac*nchan;
    unsigned ld;
    unsigned cd;
    unsigned range;

    for(int i(0);i<ndead;i++) {
        ld=lRndm.Integer(nscintlayer);
        range=scintmaxid[ld]-scintminid[ld];
        cd=scintminid[ld]+(lRndm.Integer(range));
        //std::cout<<ld<<" "<<cd<<std::endl;
        deadlist.insert(std::make_pair(ld,cd));
    }

    //########## SILICON
    std::set<std::pair<unsigned, unsigned>> deadlistsi;
    unsigned nsichan=0;
    for(int i(0);i<nsilayer;i++) {
        nsichan+=(simaxid[i]-siminid[i]);
    }
    //***** Define Silicon randomly 1% dead channels parameters
    std::cout<<"total number silicon channels is "<<nsichan<<std::endl;
    unsigned nsidead=deadfrac*nsichan;
    std::cout<<"total number of dead silicon channels is "<<nsidead<<std::endl;

    const unsigned nDeadSiliconCell= 69917;
    //const unsigned nDeadSiliconCell= 139836;

    std::pair<unsigned, unsigned> adj_to_dead[nDeadSiliconCell*2];//to define average energy in layers plus and minus 1
    std::pair<unsigned, unsigned> adj_to_dead_inlay[nDeadSiliconCell*6];//to define average energy in layer in cells plus and minus 1
    unsigned ld_si;
    unsigned cd_si;
    unsigned range_si;

    unsigned array_deadcells[nDeadSiliconCell];
    unsigned array_deadlayers[nDeadSiliconCell];
    double gausampl[nDeadSiliconCell];

    for(int i(0);i<nsidead;i++) {
        ld_si=lRndm.Integer(nsilayer);
        range_si=simaxid[ld_si]-siminid[ld_si];
        cd_si=siminid[ld_si]+(lRndm.Integer(range_si));
        deadlistsi.insert(std::make_pair(ld_si,cd_si));
    }
    unsigned n = 0;
    for(auto itr=deadlistsi.begin();itr!=deadlistsi.end();itr++ ) {
        //    std::cout<<(*itr).first<<"     "<<(*itr).second<<std::endl;
        array_deadcells[n]=(*itr).second;
        array_deadlayers[n]=(*itr).first;
        adj_to_dead[2*n] = {(*itr).first-1, (*itr).second};
        adj_to_dead[2*n+1] = {(*itr).first+1, (*itr).second};
        //    adj_to_deadinlay[2*n] = {(*itr).first, (*itr).second-1};
        //    adj_to_deadinlay[2*n+1] = {(*itr).first, (*itr).second+1};

        adj_to_dead_inlay[2*n] = {(*itr).first, (*itr).second-497};
        adj_to_dead_inlay[2*n+1] = {(*itr).first, (*itr).second-496};
        adj_to_dead_inlay[2*n+2] = {(*itr).first, (*itr).second-1};
        adj_to_dead_inlay[2*n+3] = {(*itr).first, (*itr).second+1};
        adj_to_dead_inlay[2*n+4] = {(*itr).first, (*itr).second+496};
        adj_to_dead_inlay[2*n+5] = {(*itr).first, (*itr).second+497};

        n+=1;
        //    std::cout<<"n  "<<n<<" array_deadcells  "<<array_deadcells[nDeadSiliconCell]<<std::endl;
    }



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
    std::vector<HGCSSSimHit> * alusimhitvec = 0;
    std::vector<HGCSSRecoHit> * rechitvec = 0;
    std::vector<HGCSSGenParticle> * genvec = 0;
    unsigned nPuVtx = 0;



    lSimTree->SetBranchAddress("HGCSSEvent",&event);
    lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
    lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
    //  lSimTree->SetBranchAddress("HGCSSAluSimHitVec",&alusimhitvec);
    lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);

    lRecTree->SetBranchAddress("HGCSSEvent",&eventRec);
    lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
    if (lRecTree->GetBranch("nPuVtx")) lRecTree->SetBranchAddress("nPuVtx",&nPuVtx);



    unsigned ievtRec = 0;
    unsigned nSkipped = 0;
    std::vector<double> absW;
    bool firstEvent = true;
    bool firstEvent2 = true;
    unsigned ibinScintmin[nscintlayer];
    for(int i(0);i<nscintlayer;i++) {ibinScintmin[i]=4294967295;}
    unsigned ibinScintmax[nscintlayer];
    for(int i(0);i<nscintlayer;i++) {ibinScintmax[i]=0;}
    unsigned badlaymin=10000;

    double rmaxs[nscintlayer] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double rmins[nscintlayer]={5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000};
    double rmaxns[nscintlayer] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double rminns[nscintlayer]={5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000};


    double xmax[nscintlayer] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double xmin[nscintlayer]={5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000};
    double ymax[nscintlayer] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double ymin[nscintlayer]={5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000};

    int isnap=-1;
    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
        if (ievtRec>=lRecTree->GetEntries()) continue;

        mycalib.setVertex(0.,0.,0.);

        if (debug) std::cout << std::endl<<std::endl<<"... Processing entry: " << ievt << std::endl;
        else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;
        mycalib.setVertex(0.,0.,0.);
        lSimTree->GetEntry(ievt);
        lRecTree->GetEntry(ievtRec);
        if (nPuVtx>0 && eventRec->eventNumber()==0 && event->eventNumber()!=0) {
            std::cout << " skip !" << ievt << " " << ievtRec << std::endl;
            nSkipped++;
            continue;
        }

        if(firstEvent) {
            firstEvent=false;
            double absweight=0;
            for (unsigned iL(0); iL<(*ssvec).size();++iL){
                if(iL<((*ssvec).size()-1)) {
                    unsigned next=iL+1;
                    absweight=(((*ssvec)[iL].voldEdx())+((*ssvec)[next].voldEdx()))/2. ;
                } else{
                    absweight+=(*ssvec)[iL].voldEdx()  ;
                }
                absW.push_back(absweight);
                absweight=0;
            }
        }

        double ptgen=-1.;
        double Egen=-1.;
        double ptgenpx=-1.;
        double ptgenpy=-1.;
        double ptgenpz=-1.;
        double etagen=99999.;
        double phigen=99999.;
        double thetagen=-1.;
        int pidgen=-1;
        if((*genvec).size()>0) {
            pidgen=(*genvec)[0].pdgid();
            etagen=(*genvec)[0].eta();
            phigen=(*genvec)[0].phi();
            thetagen=(*genvec)[0].theta();
            ptgenpx=(*genvec)[0].px()/1000.;
            ptgenpy=(*genvec)[0].py()/1000.;
            ptgenpz=(*genvec)[0].pz()/1000.;
            ptgen=sqrt(ptgenpx*ptgenpx+ptgenpy*ptgenpy);
            Egen=sqrt(ptgenpx*ptgenpx+ptgenpy*ptgenpy+ptgenpz*ptgenpz);
        }
        //  h_getaphi->Fill(etagen,phigen);
        bool isScint = false;
        double coneSize = 0.3;
        double rechitsum = 0;
        double rechitsumdead = 0;
        double rechitBHsum = 0;
        double rechitsumdead_Si = 0;
        double rechitsumdead_Silay8 = 0;
        double rechitsumlaypn = 0;
        double rechitsum_six_inlay = 0;
        double rechitsumgaus = 0;

        const unsigned nlay=27;
        double rechitsumnodeadinlay[nlay] = {};
        double rechitsumnodeadinlayind[nlay] = {};
        double rechitsuminlay[nlay] = {};
        double rechitsumlay[nlay] = {};
        double rechsumdead[nlay] = {};
        for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
            HGCSSRecoHit lHit = (*rechitvec)[iH];
            unsigned layer = lHit.layer();
            unsigned ixx=layer;
            double xh=lHit.get_x();
            double yh=lHit.get_y();
            double zh=lHit.get_z();
            //      double rh=sqrt(xh*xh+yh*yh);
            double lenergy=lHit.energy()*absW[layer]/1000.;
            double leta = lHit.eta();
            double lphi = lHit.phi();
            double dR=DeltaR(etagen,phigen,leta,lphi);
            //*** compute the 53 mm cut on the gen hits (Chris email)
            double ltheta = 2.*atan(exp(-leta));
            double lr=zh*tan(ltheta);
            double rgen = zh*tan(thetagen);//thetagen defined in line 847
            double xgen = rgen*cos(phigen);
            double ygen = rgen*sin(phigen);
            double Rgen = sqrt(xgen*xgen+ygen*ygen);
            double dR1 = fabs(sqrt((xgen-xh)*(xgen-xh)+(ygen-yh)*(ygen-yh)));
            //std::cout<<"dR1  "<<dR1<<std::endl;
            //***
            if(ixx>nscintlayer+scintoffset) ixx=ixx-nscintlayer-1;
            int ip=ixx-scintoffset;
            const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
            isScint = subdet.isScint;
            TH2Poly *map = isScint?(subdet.type==DetectorEnum::BHCAL1?geomConv.squareMap1():geomConv.squareMap2()): shape==4?geomConv.squareMap() : shape==2?geomConv.diamondMap() : shape==3? geomConv.triangleMap(): geomConv.hexagonMap();
            unsigned cellid = 0;
            ROOT::Math::XYZPoint pos = ROOT::Math::XYZPoint(lHit.get_x(),lHit.get_y(),lHit.get_z());
            if (isScint){//**********************
                cellid = map->FindBin(pos.eta(),pos.phi());
            } else {
                cellid = map->FindBin(lHit.get_x(),lHit.get_y());

            }
            if(dR<coneSize && dR1<53 && layer<28) {
                //        double b = 8.72;
                //	double b = 14.6304;
                //     double c = 0.0246166544214353467;
                //	double c = 0.06462;
                double b = 10.027;//13.65;//10.027;
                double c = 0.05025;//0.061294596;//0.05025;
                double gausfun = -1.0;
                rechitsum+=lenergy;
                if(!isScint) {
                    for (unsigned i=1; i<27;i++){
                        if (layer==i){
                            rechitsumlay[i]+=lenergy;
                            rechitsuminlay[i]=lenergy;
                            //h_inlay_energy[i]->Fill(rechitsuminlay[i]);
                        }
                    }
                    std::pair<unsigned,unsigned> tempsi(layer,cellid);
                    std::set<std::pair<unsigned,unsigned>>::iterator ibc=deadlistsi.find(tempsi);
                    if(ibc==deadlistsi.end()) {
                        rechitsumdead_Si+=lenergy;
                        if(layer == 8){ rechitsumdead_Silay8+=lenergy;}
                        for (unsigned i=1; i<27;i++){
                            if (layer==i){
                                rechitsumnodeadinlay[i]+=lenergy;
                                rechitsumnodeadinlayind[i]=lenergy;
                            }
                        }
                    } else {
                        for (unsigned i=1; i<27;i++){
                            if(layer==i){
                                rechsumdead[i]+=lenergy;
                            }
                        }
                    }
                    for(int i(0); i < nDeadSiliconCell; i++){
                        if(cellid == array_deadcells[i]){
                            gausfun = exp(-0.5*c*(layer-b));
                            if(layer == array_deadlayers[i]+1){
                                gausampl[i]= lenergy/gausfun;
                                //		std::cout<<"amplitude: "<<gausampl[i]<<std::endl;
                            }
                            if(layer == array_deadlayers[i]){
                                rechitsumgaus+= gausampl[i] * gausfun;
                                //		std::cout<<"one dead cell is found: "<<std::endl;
                            }
                        }
                    }
                    for(int i(0); i <2*nDeadSiliconCell; i++){
                        if(cellid == adj_to_dead[i].second && layer == adj_to_dead[i].first){
                            rechitsumlaypn += lenergy/2;
                        }
                    }
                    for(int i(0); i <6*nDeadSiliconCell; i++){
                        if(cellid == adj_to_dead_inlay[i].second && layer == adj_to_dead[i].first){
                            rechitsum_six_inlay += lenergy/6;
                            //	      std::cout<<"inlay rechitsumave: "<<rechitsum_six_inlay<<std::endl;
                        }
                    }
                }
            }
        }  // end loop over hits
        h_rechitsum->Fill(rechitsum);
        ievtRec++;
    }//loop on entries

    outputFile->cd();
    outputFile->Write();
    outputFile->Close();

    return 0;

}//main
