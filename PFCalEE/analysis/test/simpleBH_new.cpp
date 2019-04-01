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


int main(int argc, char** argv){ //main
    ///////////////////////////////////////////////////////////////
    // initialize some variables
    ///////////////////////////////////////////////////////////////

    const unsigned nscintlayer=16;    //! scintillator layers
    const unsigned scintoffset=36;    //! ??

    //! Create an array of booleans and initialize to true
    bool ipfirst[nscintlayer];
    for(int i(0);i<nscintlayer;i++) {ipfirst[i]=true;}

    //! Define the the id ranges for each scintillator layer
    unsigned scintminid[nscintlayer]={145,179,219,206,159,205,274,284,28,223,47,24,219,1,188,211};
    unsigned scintmaxid[nscintlayer]={32545,32579,32619,32606,20823,20869,20938,20948,20692,20887,20711,20688,20883,20665,20852,20875};


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
    po::options_description preconfig("Configuration");  //! use the boost library to create program options
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
    const double threshMin = 0.5;     //! not sure

    std::cout << " ---- Selection settings: ---- " << std::endl
    << " -------threshMin " << threshMin << std::endl
    << " ------------------------------------------" << std::endl;



    /////////////////////////////////////////////////////////////
    //input
    /////////////////////////////////////////////////////////////


    std::ostringstream inputsim;
    inputsim << filePath << "/" << simFileName; //! pass dir+filename for simulation file in this ostringstream
    std::ostringstream inputrec;
    if (digifilePath.size()==0)
    inputrec << filePath << "/" << recoFileName;
    else
    inputrec << digifilePath << "/" << recoFileName; //! pass dir+filename for simulation file in this ostringstream

    std::cout << inputsim.str() << " " << inputrec.str() << std::endl;

    HGCSSInfo * info; //! Can store a few variables describing version, model, cellSize and shape for HGC

    TChain *lSimTree = new TChain("HGCSSTree"); //! Create a TChain to append multiple TTrees
    TChain *lRecTree = 0;

    TFile * simFile = 0;
    TFile * recFile = 0;

    if (recoFileName.find("Digi") != recoFileName.npos) //! is this correct?
    lRecTree = new TChain("RecoTree");
    else lRecTree = new TChain("PUTree");

    if (nRuns == 0){ //! Execute only during first run - Get all simulation files and append to TChains
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

    const unsigned versionNumber = info->version(); //! Retrieve values from info object
    const unsigned model = info->model();
    const unsigned shape = info->shape();
    const double cellSize = info->cellSize();
    const double calorSizeXY = info->calorSizeXY();

    bool isTBsetup = (model != 2);
    bool bypassR = false;
    if (isTBsetup) bypassR = true;

    HGCSSDetector & myDetector = theDetector(); //! define an object which is friend to HGCSSDetector class. It can access the detector geometry parameter arrays
    myDetector.buildDetector(versionNumber,true,false,bypassR); //! Initializes the detector with all subdetectors inside

    HGCSSCalibration mycalib(inFilePath); //! Define a calibration object

    //corrected for Si-Scint overlap
    const unsigned nLayers = 52;//


    std::cout << " -- Calor size XY = " << calorSizeXY
    << ", version number = " << versionNumber
    << ", model = " << model << std::endl
    << " -- cellSize = " << cellSize
    << ", shape = " << shape
    << ", nLayers = " << nLayers
    << std::endl;


    HGCSSGeometryConversion geomConv(model,cellSize,bypassR,3); //! not sure what's the functionality
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
    const int nsnap=10;

    TH1F* h_ECone = new TH1F("h_ECone","rechit Energy;E_{rechit}(GeV);Events",1500,0,1500);

    ///////////////////////////////
    // for missing channel study //
    ///////////////////////////////

    //const double deadfrac=0.00003;
    //const double deadfrac=0.5;
    std::set<std::pair<unsigned, unsigned>> deadlist;
    unsigned nchan=0;
    for(int i(0);i<nscintlayer;i++) {
        nchan+=(scintmaxid[i]-scintminid[i]);
    }
    std::cout<<"total number scintillator channels is "<<nchan<<std::endl;
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

    //std::cout<<" dead list is "<<std::endl;
    //for(auto itr=deadlist.begin();itr!=deadlist.end();itr++ ) {
    //std::cout<<(*itr).first<<" "<<(*itr).second<<std::endl;
    //}


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
    //lSimTree->SetBranchAddress("HGCSSAluSimHitVec",&alusimhitvec);
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

        lSimTree->GetEntry(ievt);
        lRecTree->GetEntry(ievtRec);
        if (nPuVtx>0 && eventRec->eventNumber()==0 && event->eventNumber()!=0) {
            std::cout << " skip !" << ievt << " " << ievtRec << std::endl;
            nSkipped++;
            continue;
        }

        if(firstEvent) {
            firstEvent=false;
            std::cout<<" size of ssvec of weights is "<<(*ssvec).size()<<std::endl;
            double absweight=0;
            for (unsigned iL(0); iL<(*ssvec).size();++iL){
                if(iL<((*ssvec).size()-1)) {
                    unsigned next=iL+1;
                    absweight=(((*ssvec)[iL].voldEdx())+((*ssvec)[next].voldEdx()))/2.;
                } else{
                    absweight+=(*ssvec)[iL].voldEdx();
                }
                absW.push_back(absweight);
                absweight=0;
            }
            std::cout << " -- AbsWeight size: " << absW.size() << std::endl;
            std::cout<<" values are ";
            for (unsigned iL(0); iL<(*ssvec).size();++iL){
                std::cout<<" "<<absW[iL];
            }
            std::cout<<std::endl;
        }

        //
        //!information about generator level particles

        double ptgen=-1.;
        double Egen=-1.;
        double ptgenpx=-1.;
        double ptgenpy=-1.;
        double ptgenpz=-1.;
        double etagen=99999.;
        double phigen=99999.;
        double thetagen=-1.;
        int pidgen=-1;
        if((*genvec).size()>0) { //!check if particle entry is valid and fill particle variables
            pidgen=(*genvec)[0].pdgid();
            ptgenpx=(*genvec)[0].px()/1000.;
            ptgenpy=(*genvec)[0].py()/1000.;
            ptgenpz=(*genvec)[0].pz()/1000.;
            ptgen=sqrt(ptgenpx*ptgenpx+ptgenpy*ptgenpy);
            etagen=(*genvec)[0].eta();
            Egen=sqrt(ptgenpx*ptgenpx+ptgenpy*ptgenpy+ptgenpz*ptgenpz);
            phigen=(*genvec)[0].phi();
            thetagen=(*genvec)[0].theta();
        }
        double etaaxis=etagen;
        double phiaxis=phigen;


        bool isScint = false;

        // make some simple plots about all the raw rechits, without the weights
        unsigned iMax=-1;
        double MaxE=-1.;
        bool snap=false;
        //double aistep=20000;
        double aistep=400;
        int istep=aistep;
        double stepsize=6000/aistep;
        for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
            HGCSSRecoHit lHit = (*rechitvec)[iH];
            double xh=lHit.get_x();
            double yh=lHit.get_y();
            double zh=lHit.get_z();
            double rh=sqrt(xh*xh+yh*yh);
            double Eh=lHit.energy();
            double leta = lHit.eta();
            double lphi = lHit.phi();
            unsigned layer = lHit.layer();
            unsigned ixx=layer;

            if(ixx>nscintlayer+scintoffset) ixx=ixx-nscintlayer-1;
            int ip=ixx-scintoffset;
            const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
            isScint = subdet.isScint;



            TH2Poly *map = isScint?(subdet.type==DetectorEnum::BHCAL1?geomConv.squareMap1():geomConv.squareMap2()): shape==4?geomConv.squareMap() : shape==2?geomConv.diamondMap() : shape==3? geomConv.triangleMap(): geomConv.hexagonMap();

            unsigned cellid = 0;
            ROOT::Math::XYZPoint pos = ROOT::Math::XYZPoint(lHit.get_x(),lHit.get_y(),lHit.get_z());
            if (isScint){
                double aaaphi = pos.phi();
                ///if(aaaphi<0) aaaphi+=2.*TMath::Pi();
                cellid = map->FindBin(pos.eta(),aaaphi);
                //cellid = map->FindBin(pos.eta(),pos.phi());
            } else {
                cellid = map->FindBin(lHit.get_x(),lHit.get_y());
            }
            geomConv.fill(layer,Eh,0,cellid,zh);

            if(isScint) {
                if(firstEvent2) {
                    firstEvent2=false;
                    std::cout<<std::endl<<std::endl;
                    std::cout<<"check check layer is "<<layer<<std::endl;
                    double abceta=1.7;
                    double abcphi;
                    unsigned acellid=0;
                    for(int i(0);i<100;i++) {
                        abcphi=(6.2/100.)*i;
                        acellid = map->FindBin(abceta,abcphi);
                        std::cout<<abceta<<" "<<abcphi<<" "<<acellid<<std::endl;
                    }
                    abceta=2.0;
                    for(int i(0);i<100;i++) {
                        abcphi=(6.2/100.)*i;
                        acellid = map->FindBin(abceta,abcphi);
                        std::cout<<abceta<<" "<<abcphi<<" "<<acellid<<std::endl;
                    }
                }
            }

            if(Eh>MaxE) {
                MaxE=Eh;
                iMax=iH;
            }

            if(ip>=0) { // look at layers with scintillator
                if(xh>5000) {  // weird hits
                    if(layer<badlaymin) badlaymin=layer;
                    if(!snap) {
                        if(isnap<nsnap-1) {
                            snap=true;
                            isnap+=1;
                            std::cout<<" isnap is "<<isnap<<std::endl;
                        }
                    }
                    double rr=zh*tan(thetagen);
                    double xx=rr*cos(phigen);
                    double yy=rr*sin(phigen);
                } else {  //normal his
                }
                if(cellid<4000000000) {  // temp fix until Anne-Marie fixes cell geo problem
                    if(xh<xmin[ip]) {xmin[ip]=xh;}
                    if(xh>xmax[ip]) {xmax[ip]=xh;}
                    if(yh<ymin[ip]) {ymin[ip]=yh;}
                    if(yh>ymax[ip]) {ymax[ip]=yh;}

                    if(isScint) {
                        //	  std::cout<<"ip rh is "<<ip<<" "<<rh<<std::endl;
                        //std::cout<<rmins[ip]<<" "<<rmaxs[ip]<<std::endl;

                        if(domap) {
                            if(ipfirst[ip]) {
                                ipfirst[ip]=false;
                                std::cout<<"doing ip "<<ip<<std::endl;
                                for(int iia(0);iia<istep;iia++) {
                                    for(int iib(0);iib<istep;iib++) {
                                        //double xxa=-3000.+stepsize*iia;
                                        //double yya=-3000.+stepsize*iib;
                                        double etaa=1.2+((4.-1.2)/aistep)*iia;
                                        double phia=-3.14+((6.4)/aistep)*iib;
                                        double thetaa=2.*atan(exp(-etaa));
                                        double ra=lHit.get_z()*tan(thetaa);
                                        double xxa=ra*cos(phia);
                                        double yya=ra*sin(phia);
                                        ROOT::Math::XYZPoint apos = ROOT::Math::XYZPoint(xxa,yya,lHit.get_z());
                                        //double bbbphi = phia;
                                        double bbbphi=apos.phi();
                                        //if(bbbphi<0) bbbphi+=2.*TMath::Pi();
                                        //unsigned acellid = map->FindBin(etaa,bbbphi);
                                        unsigned acellid = map->FindBin(apos.eta(),bbbphi);
                                        //unsigned acellid = map->FindBin(apos.eta(),apos.phi());
                                        if(acellid<4000000000) {
                                            if((apos.x()<-500)&&(apos.y()<-500) ) {
                                                if((apos.x()>-600)&&(apos.y()>-600) ) {
                                                    //std::cout<<"michael "<<ip<<" "<<apos.x()<<" "<<apos.y()<<" "<<apos.eta()<<" "<<apos.phi()<<" "<<acellid<<std::endl;
                                                }
                                            }
                                            if((apos.x()>500)&&(apos.y()>500) ) {
                                                if((apos.x()<600)&&(apos.y()<600) ) {
                                                    //std::cout<<"michael "<<ip<<" "<<apos.x()<<" "<<apos.y()<<" "<<apos.eta()<<" "<<apos.phi()<<" "<<acellid<<std::endl;
                                                }
                                            }

                                            if(acellid>ibinScintmax[ip]) ibinScintmax[ip]=acellid;
                                            if(acellid<ibinScintmin[ip]) ibinScintmin[ip]=acellid;
                                        }
                                    }
                                }
                            }
                        }

                        if(rh<rmins[ip]) {rmins[ip]=rh;}
                        if(rh>rmaxs[ip]) {rmaxs[ip]=rh;}
                        //std::cout<<rmins[ip]<<" "<<rmaxs[ip]<<std::endl;

                    } else {  // end if scint
                        if(rh<rminns[ip]) {rminns[ip]=rh;}
                        if(rh>rmaxns[ip]) {rmaxns[ip]=rh;}
                    }
                }  // end temp fix
            }  // end loooking at layers with scipt (ip>0)
        }  // end loop over hits



        HGCSSRecoHit lHit = (*rechitvec)[iMax];
        double maxeta = lHit.eta();
        double maxphi=lHit.phi();
        double maxE=lHit.energy();


        // make e/p plots for various cones around gen particle, using weights this time
        // for now, scale up by 10.  don't know why.  asking Anne-Marie
        // also change to GeV for this section
        const unsigned isize=5;
        double coneSize[isize]={0.1,0.2,0.3,0.4,0.5};
        double rechitsum[isize]={0.,0.,0.,0.,0.};
        double rechitsumdead[isize]={0.,0.,0.,0.,0.};
        double rechitBHsum[isize]={0.,0.,0.,0.};

        double etaW=0.;
        double phiW=0.;
        double norm=0.;
        for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
            HGCSSRecoHit lHit = (*rechitvec)[iH];
            unsigned layer = lHit.layer();
            int ip=layer-scintoffset-nscintlayer-1;
            const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
            isScint = subdet.isScint;
            double leta = lHit.eta();
            double lphi = lHit.phi();
            double lenergy=lHit.energy()*absW[layer]/1000.;
            TH2Poly *map = isScint?(subdet.type==DetectorEnum::BHCAL1?geomConv.squareMap1():geomConv.squareMap2()): shape==4?geomConv.squareMap() : shape==2?geomConv.diamondMap() : shape==3? geomConv.triangleMap(): geomConv.hexagonMap();
            unsigned cellid = 0;
            ROOT::Math::XYZPoint pos = ROOT::Math::XYZPoint(lHit.get_x(),lHit.get_y(),lHit.get_z());
            if (isScint){
                double aaaphi = pos.phi();
                //if(aaaphi<0) aaaphi+=2.*TMath::Pi();
                cellid = map->FindBin(pos.eta(),aaaphi);
                //cellid = map->FindBin(pos.eta(),pos.phi());
            } else {
                cellid = map->FindBin(lHit.get_x(),lHit.get_y());
            }

            //clean up rechit collection
            norm+=lenergy;
            etaW+=leta*lenergy;
            phiW+=lphi*lenergy;

            double dR=DeltaR(etaaxis,phiaxis,leta,lphi);
            //double dR=fabs(etagen-leta);
            double xh=lHit.get_x();
            double yh=lHit.get_y();
            double zh=lHit.get_z();
            double rgen = zh*tan(thetagen);//thetagen defined in line 847
            double xgen = rgen*cos(phigen);
            double ygen = rgen*sin(phigen);
            double Rgen = sqrt(xgen*xgen+ygen*ygen);
            double dR1 = fabs(sqrt((xgen-xh)*(xgen-xh)+(ygen-yh)*(ygen-yh)));
            for(unsigned ii(0);ii<isize;ii++) {
                if(dR<coneSize[ii] && dR1<53 && layer<28) {
                    rechitsum[ii]+=lenergy;

                    if(isScint) {
                        rechitBHsum[ii]+=lenergy;
                        std::pair<unsigned,unsigned> temp(ip,cellid);
                        std::set<std::pair<unsigned,unsigned>>::iterator ibc=deadlist.find(temp);
                        //std::cout<<"michael "<<temp.first<<" "<<temp.second<<std::endl;
                        //std::cout<<" I am here"<<std::endl;
                        //if(ibc!=deadlist.end()) std::cout<<" found on dead list"<<std::endl;
                        if(ibc==deadlist.end()) {
                            rechitsumdead[ii]+=lenergy;
                        }
                    } else {
                        rechitsumdead[ii]+=lenergy;
                    }
                }
            }
        }//loop on hits

        double frac=-0.05;
        double notBH=rechitsum[3]-rechitBHsum[3];
        if(rechitsum[3]>0) frac=rechitBHsum[3]/rechitsum[3];

        h_ECone->Fill(rechitsum[3]);

        geomConv.initialiseHistos();
        ievtRec++;
    }//loop on entries

    outputFile->cd();
    outputFile->Write();
    outputFile->Close();

    return 0;
}
