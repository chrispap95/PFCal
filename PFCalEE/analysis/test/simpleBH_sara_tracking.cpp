#include<string>
#include<iostream>
#include<iomanip>
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
#include "PositionFit.hh"
#include "SignalRegion.hh"
#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"
#include "utilities.h"

using boost::lexical_cast;
namespace po=boost::program_options;

double zlayer[27] = {3207.14, 3222.53, 3231.62, 3246.91, 3256.01, 3271.31, 3280.41, 3295.71,
    3304.81, 3320.11, 3329.2, 3344.51, 3353.6, 3368.9, 3378, 3393.31, 3402.4, 3417.71,
    3426.81, 3442.11, 3451.2, 3466.51, 3475.61, 3490.91, 3500.02, 3515.32, 3524.43
};

double DeltaR(double eta1,double phi1,double eta2,double phi2){
    double dr=99999.;
    double deta=fabs(eta1-eta2);
    double dphi=fabs(phi1-phi2);
    if(dphi>TMath::Pi()) dphi=2.*TMath::Pi()-dphi;
    dr=sqrt(deta*deta+dphi*dphi);
    return dr;
}

std::pair<double, double> TrackId(int layer, double theta, double phi){
    double z = zlayer[layer-1];
    double x = z*cos(phi)*tan(theta);
    double y = z*sin(phi)*tan(theta);
    return std::make_pair(x,y);
}

int main(int argc, char** argv){
    /**********************************
    ** initialize some variables
    **********************************/
    const unsigned nscintlayer=16;
    const unsigned scintoffset=36;

    //Input output and config options
    std::string cfg;
    unsigned pNevts;
    std::string outFilePath;
    std::string filePath;
    std::string digifilePath;
    unsigned nRuns;
    std::string simFileName;
    std::string recoFileName;
    std::string TrFilePath;
    unsigned debug;
    bool Trsample;
    po::options_description preconfig("Configuration");
    preconfig.add_options()("cfg,c",po::value<std::string>(&cfg)->required());
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(preconfig).allow_unregistered().run(), vm);
    po::notify(vm);
    po::options_description config("Configuration");
    config.add_options()
    //Input output and config options //->required()
    ("pNevts,n",        po::value<unsigned>(&pNevts)->default_value(0))
    ("outFilePath,o",   po::value<std::string>(&outFilePath)->required())
    ("filePath,i",      po::value<std::string>(&filePath)->required())
    ("digifilePath",    po::value<std::string>(&digifilePath)->default_value(""))
    ("simFileName,s",   po::value<std::string>(&simFileName)->required())
    ("recoFileName,r",  po::value<std::string>(&recoFileName)->required())
    ("nRuns",           po::value<unsigned>(&nRuns)->default_value(0))
    ("debug,d",         po::value<unsigned>(&debug)->default_value(0))
    ("Trsample",        po::value<bool>(&Trsample)->default_value(1)) //Generate Tr study training sample
    ("TrFilePath",      po::value<std::string>(&TrFilePath)->default_value("training_sample.root")) //File to export data for Tr
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

    /**********************************
    ** Input
    **********************************/
    // Get Sim and Digi files
    std::ostringstream inputsim;
    inputsim << filePath << "/" << simFileName;
    std::ostringstream inputrec;
    if (digifilePath.size()==0) inputrec << filePath << "/" << recoFileName;
    else inputrec << digifilePath << "/" << recoFileName;

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
            std::cout << " -- \e[31mError\e[0m in getting information from simfile!" << std::endl;
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
                    std::cout << " -- \e[31mError\e[0m in getting information from simfile!" << std::endl;
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
        std::cout << " -- \e[31mError\e[0m, tree HGCSSTree cannot be opened. Exiting..." << std::endl;
        return 1;
    }
    if (!lRecTree){
        std::cout << " -- \e[31mError\e[0m, tree RecoTree cannot be opened. Exiting..." << std::endl;
        return 1;
    }

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

    /**********************************
    ** Output
    **********************************/
    // Define output file and the histograms contained
    TFile *outputFile = TFile::Open(outFilePath.c_str(),"RECREATE");

    if (!outputFile) {
        std::cout << " -- \e[31mError\e[0m, output file " << outFilePath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
        return 1;
    }
    else std::cout << " -- opening output file: " << outputFile->GetName() << " \e[32msuccess\e[0m" << std::endl;

    outputFile->cd();

    /**********************************
    ** Track Study output section
    **     - Trlayer is track cell layer
    **     - Treta is gen eta
    **     - Trphi is gen phi
    **     - Trni is ith track cell neighbor
    **     - Trtrack is track cell rechit
    ** We also need to create a ?set? container to store the values before writing to TTree
    **********************************/
    TFile* fout = new TFile(TrFilePath.c_str(),"RECREATE");
    float Trlayer,Trcellid, Treta, Trphi, Trn1, Trn2, Trn3, Trn4, Trn5, Trn6, Trtrack, Trnup, Trndown;
    TTree* t1 = new TTree("t1","sample");
    t1->Branch("Trlayer",&Trlayer,"Trlayer/F");
    t1->Branch("Trcellid",&Trcellid,"Trcellid/F");
    t1->Branch("Treta",&Treta,"Treta/F");
    t1->Branch("Trphi",&Trphi,"Trphi/F");
    t1->Branch("Trn1",&Trn1,"Trn1/F");
    t1->Branch("Trn2",&Trn2,"Trn2/F");
    t1->Branch("Trn3",&Trn3,"Trn3/F");
    t1->Branch("Trn4",&Trn4,"Trn4/F");
    t1->Branch("Trn5",&Trn5,"Trn5/F");
    t1->Branch("Trn6",&Trn6,"Trn6/F");
    t1->Branch("Trtrack",&Trtrack,"Trtrack/F");
    t1->Branch("Trnup",&Trnup,"Trnup/F");
    t1->Branch("Trndown",&Trndown,"Trndown/F");

    /*
    ** Define a vector of the array:
    ** {track cell layer, tr id, tr eta, tr phi, Trn1, Trn2, Trn3, Trn4, Trn5, Trn6, Trnup,
    ** Trndown, track cell rechit}
    */
    std::set<std::pair<unsigned, unsigned>> tracklist;
    std::vector<std::array<float, 13>> Trvectorev;

    /**********************************
    **  Start event loop
    **********************************/
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
    std::vector<double> absW;
    bool firstEvent = true;

    for (unsigned ievt(0); ievt<1; ++ievt){//loop on entries (events)
        if (ievtRec>=lRecTree->GetEntries()) continue;

        mycalib.setVertex(0.,0.,0.);

        int digits = 0;
        int ievt_temp = ievt;
        while (ievt_temp) {
            ievt_temp /= 10;
            digits++;
        }
        std::string mv_dl = "\e[2D\e[0K";
        if(firstEvent) {
            std::cout <<"\n\n... Processing entry: " << ievt;
        }else {
            std::cout << mv_dl << ievt;
        }

        mycalib.setVertex(0.,0.,0.);
        lSimTree->GetEntry(ievt);
        lRecTree->GetEntry(ievtRec);
        if (nPuVtx>0 && eventRec->eventNumber()==0 && event->eventNumber()!=0) {
            std::cout << " \e[31mskip!\e[0m" << ievt << " " << ievtRec << std::endl;
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
                } else absweight+=(*ssvec)[iL].voldEdx();
                absW.push_back(absweight);
                absweight=0;
            }
        }


        double etagen=99999.;
        double phigen=99999.;
        double thetagen=-1.;
        if((*genvec).size()>0) {
            etagen=(*genvec)[0].eta();
            phigen=(*genvec)[0].phi();
            thetagen=(*genvec)[0].theta();
        }

        bool isScint = false;
        double coneSize = 0.3;
        double rechitsum = 0;

        const unsigned nlay=27;

        /*
        ** Now we need to populate a list with the predicted cellids for each
        ** layer. This accomplished through the use of TrackId() for each layer.
        */
        TH2Poly* map2 = geomConv.hexagonMap();
        std::array<float, 13> Trarr;
        for(unsigned k(0); k < 13; ++k) Trarr[k] = 0;
        for (unsigned iL(0); iL < nlay; ++iL){
            Trarr[0] = iL;
            std::pair<double, double> xypair = TrackId(iL, thetagen, phigen);
            Trarr[1] = map2->FindBin(xypair.first, xypair.second);
            Trvectorev.push_back(Trarr);
            tracklist.insert(std::make_pair(Trarr[0],Trarr[1]));
        }

        // Populate vectors for neighboring cells of track cells
        std::vector<std::pair<unsigned, unsigned>> neighbors;
        std::vector<std::pair<unsigned, unsigned>> neighbors_inlay;
        for(auto itr=tracklist.begin();itr!=tracklist.end();itr++ ) {
            neighbors.push_back({(*itr).first-1, (*itr).second});
            neighbors.push_back({(*itr).first+1, (*itr).second});
            neighbors_inlay.push_back({(*itr).first, (*itr).second-497});
            neighbors_inlay.push_back({(*itr).first, (*itr).second-496});
            neighbors_inlay.push_back({(*itr).first, (*itr).second-1});
            neighbors_inlay.push_back({(*itr).first, (*itr).second+1});
            neighbors_inlay.push_back({(*itr).first, (*itr).second+496});
            neighbors_inlay.push_back({(*itr).first, (*itr).second+497});
        }

        for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
            HGCSSRecoHit lHit = (*rechitvec)[iH];
            unsigned layer = lHit.layer();
            unsigned ixx=layer;
            double xh=lHit.get_x();
            double yh=lHit.get_y();
            double zh=lHit.get_z();
            double lenergy=lHit.energy()*absW[layer]/1000.;
            double leta = lHit.eta();
            double lphi = lHit.phi();
            double dR=DeltaR(etagen,phigen,leta,lphi);
            //*** compute the 53 mm cut on the gen hits (Chris email)
            double rgen = zh*tan(thetagen);//thetagen defined in line 847
            double xgen = rgen*cos(phigen);
            double ygen = rgen*sin(phigen);
            double dR1 = fabs(sqrt((xgen-xh)*(xgen-xh)+(ygen-yh)*(ygen-yh)));

            if(ixx>nscintlayer+scintoffset) ixx=ixx-nscintlayer-1;
            const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
            isScint = subdet.isScint;
            TH2Poly *map = isScint?(subdet.type==DetectorEnum::BHCAL1?geomConv.squareMap1():geomConv.squareMap2()): shape==4?geomConv.squareMap() : shape==2?geomConv.diamondMap() : shape==3? geomConv.triangleMap(): geomConv.hexagonMap();
            unsigned cellid = 0;
            ROOT::Math::XYZPoint pos = ROOT::Math::XYZPoint(lHit.get_x(),lHit.get_y(),lHit.get_z());
            if (isScint) cellid = map->FindBin(pos.eta(),pos.phi());
            else cellid = map->FindBin(lHit.get_x(),lHit.get_y());

            if(dR<coneSize && dR1<53 && layer<28) {
                rechitsum+=lenergy;
                if(!isScint) {
                    /* Track code
                    ** Input track cells eta, phi and rechits
                    */
                    for(auto itr = Trvectorev.begin(); itr != Trvectorev.end(); itr++) {
                        if((*itr)[0] == layer && (*itr)[1] == cellid){
                            (*itr)[2] = etagen;
                            (*itr)[3] = phigen;
                            (*itr)[10] = lenergy;
                        }
                    }

                    for(unsigned i(0); i < 6*27-5; i++){
                        //Loop over track cells
                        //Different layer averaging
                        if(i < 2*27-1) {
                            if(cellid == neighbors[i].second && layer == neighbors[i].first) {
                                /* Track code
                                ** Find different layers neighbors
                                ** Hint: Incorporate to Loop below...
                                */
                                int n_track = (int)((float)i/2.);
                                int j;
                                if (i-2*n_track == 0) j = -1;
                                else if (i-2*n_track == 1) j = 1;
                                else {
                                    std::cout << "check 1 error" << std::endl;
                                    return 1;
                                }
                                // Write the data in Trvectorev vector
                                auto itr = Trvectorev.begin();
                                while (itr != Trvectorev.end()) {
                                    if((*itr)[0]+j == layer && (*itr)[1] == cellid) {
                                        (*itr)[i-2*n_track+11] = lenergy;
                                        break;
                                    }
                                    itr++;
                                }
                            }
                        }

                        //In layer averaging
                        if(cellid == neighbors_inlay[i].second && layer == neighbors_inlay[i].first){
                            /* Track code
                            ** First find which is the dead cell and what neighbor we are at
                            */
                            int n_track = (int)((float)i/6.);
                            int j;
                            if (i-6*n_track == 0) j = -497;
                            else if (i-6*n_track == 1) j = -496;
                            else if (i-6*n_track == 2) j = -1;
                            else if (i-6*n_track == 3) j = 1;
                            else if (i-6*n_track == 4) j = 496;
                            else j = 497;
                            // Write the data in Trvectorev vector
                            auto itr = Trvectorev.begin();
                            while (itr != Trvectorev.end()) {
                                if((*itr)[0] == layer && (*itr)[1]+j == cellid) {
                                    (*itr)[i-6*n_track+4] = lenergy;
                                    break;
                                }
                                itr++;
                            }
                        }
                    }
                }
            }
        }// end loop over hits
        ievtRec++;
    }//loop on entries

    std::cout << "\nCreating Track sample ..." << std::endl;
    for(auto itr = Trvectorev.begin(); itr != Trvectorev.end(); ++itr) {
        if ((*itr)[2]>0) {
            Trlayer = (*itr)[0];
            Trcellid = (*itr)[1];
            Treta = (*itr)[2];
            Trphi = (*itr)[3];
            Trn1 = (*itr)[4];
            Trn2 = (*itr)[5];
            Trn3 = (*itr)[6];
            Trn4 = (*itr)[7];
            Trn5 = (*itr)[8];
            Trn6 = (*itr)[9];
            Trtrack = (*itr)[10];
            Trnup = (*itr)[11];
            Trndown = (*itr)[12];
            t1->Fill();
        }
    }
    fout->cd();
    fout->Write();
    fout->Close();

    std::cout << "Writing files ..." << std::endl;
    outputFile->cd();
    outputFile->Write();
    outputFile->Close();

    return 0;
}
