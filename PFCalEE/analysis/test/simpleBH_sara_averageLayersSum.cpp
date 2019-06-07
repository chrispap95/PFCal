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

int main(int argc, char** argv){

    /**********************************
    ** initialize some variables
    **********************************/

    const unsigned nsilayer=52;
    //const unsigned nDeadSiliconCell=620100;//DPG presentation slide number 3

    // Define cell id boundaries

    unsigned siminid[nsilayer]={34987,34497,33997,33510,33007,33001,32505,32011,31515,31024,30531,30518,30025,29528,29035,28548,28045,
        28036,27542,27049,26552,26062,25562,25556,25057,24566,24066,23576,22080,20105,18607,17114,15136,13643,12150,10177,34987,34497,
        33507,33007,33001,32505,32011,31515,31024,30531,30518,30025,29528,29035,28548};
    unsigned simaxid[nsilayer]={250502,250992,251492,251979,2522479,252488,252982,253475,253972,254465,254955,254971,255464,255958,
        256454,256941,257441,257451,257947,258437,258937,259427,259927,259930,260430,260920,261420,261910,263406,265382,266882,268375,
        270347,271844,273337,275316,250502,250992,251979,252479,252488,252982,253475,253972,254465,254955,254971,255464,255958,256454,
        256941};

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
    bool adjacent;
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
    ("deadfrac",        po::value<double>(&deadfrac)->default_value(0))
    ("adjacent",        po::value<bool>(&adjacent)->default_value(0)) //Conduct study for restricted number of adjacent dead cells
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

    //global threshold to reduce size of noise hits
    const double threshMin = 0.5;

    std::cout << " ---- Selection settings: ---- " << std::endl
    << " -------threshMin " << threshMin << std::endl
    << " ------------------------------------------" << std::endl;


    /**********************************
    ** Input
    **********************************/
    // Get Sim and Digi files
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
        std::cout << " -- Error, output file " << outFilePath << " cannot be opened. Please create output directory. Exiting..." << std::endl;
        return 1;
    }
    else {
        std::cout << " -- output file " << outputFile->GetName() << " successfully opened." << std::endl;
    }
    outputFile->cd();

    TH1F* h_rechitsum = new TH1F("h_rechitsum","Rechitsum silicon",700,0,700.);
    TH1F* h_rechitsumdead_Si = new TH1F("h_rechitsumdead_Si","Rechitsum dead silicon",700,0,700.);
    TH1F* h_rechitsumave = new TH1F("h_rechitsumave","Sum energy average method",700,0,700.);
    TH1F* h_rechitsumave_layerSum = new TH1F("h_rechitsumave_layerSum","Sum energy average layer sums method",700,0,700.);


    /**********************************
    ** for missing channel study
    **********************************/

    // SILICON
    std::set<std::pair<unsigned, unsigned>> deadlistsi;
    unsigned nsichan=0;
    for(unsigned i(0);i<nsilayer;i++) {
        nsichan+=(simaxid[i]-siminid[i]);
    }
    //***** Define Silicon randomly 1% dead channels parameters
    std::cout<<"total number silicon channels is "<<nsichan<<std::endl;
    unsigned nsidead=deadfrac*nsichan;
    std::cout<<"total number of dead silicon channels is "<<nsidead<<std::endl;

    unsigned ld_si;
    unsigned cd_si;
    unsigned range_si;

    // Kill cells and calculate statistics on adjacent dead cells
    unsigned N_try = 0; // Number of trials to kill cells
    float N_cluster2 = 0; // Number of dead cells clusters (n_dead = 2)
    float N_clusters = 0; // Number of dead cells clusters (n_dead > 2)
    for(unsigned i(0);i<nsidead;i++) {
        N_try++;
        ld_si=lRndm.Integer(nsilayer);
        range_si=simaxid[ld_si]-siminid[ld_si];
        cd_si=siminid[ld_si]+(lRndm.Integer(range_si));
        // Enforce that any dead cell has no more than one adjacent dead cell
        unsigned adj_ok = 0;
        if (adjacent) {
            if (deadlistsi.find(std::make_pair(ld_si, cd_si-497)) != deadlistsi.end()) adj_ok++;
            if (deadlistsi.find(std::make_pair(ld_si, cd_si-496)) != deadlistsi.end()) adj_ok++;
            if (deadlistsi.find(std::make_pair(ld_si, cd_si-1)) != deadlistsi.end()) adj_ok++;
            if (deadlistsi.find(std::make_pair(ld_si, cd_si+1)) != deadlistsi.end()) adj_ok++;
            if (deadlistsi.find(std::make_pair(ld_si, cd_si+496)) != deadlistsi.end()) adj_ok++;
            if (deadlistsi.find(std::make_pair(ld_si, cd_si+497)) != deadlistsi.end()) adj_ok++;
            if (adj_ok > 1) {
                --i;
                N_cluster2++;
                continue;
            }
            if (adj_ok < 2) {
                deadlistsi.insert(std::make_pair(ld_si,cd_si));
            }
            if (adj_ok == 1) N_clusters++;
        }
        else {
            deadlistsi.insert(std::make_pair(ld_si,cd_si));
        }
    }

    // Print statistics on adjacent dead cells
    if (adjacent) {
        std::cout << std::string(120,'-') << std::endl
        << std::string(49,'-') << " Dead cells statistics "
        << std::string(48,'-') << std::endl
        << std::string(120,'-') << std::endl
        << "Number of dead cells clusters: " << N_clusters << std::endl
        << "Fraction of dead cluster cells: "
        << N_clusters*2./13983553. << std::endl
        << "Fraction of dead cells having a dead neighbor: "
        << N_clusters*2./nsidead << std::endl
        << "Dead fraction: " << deadfrac << std::endl
        << "Times the code tried to create clusters with more than 2 dead cells: "
        << N_cluster2 << std::endl
        << std::string(120,'-') << std::endl;
    }

    std::vector<std::pair<unsigned, unsigned>> adj_to_dead;//to define average energy in layers plus and minus 1
    unsigned n = 0;
    for(auto itr=deadlistsi.begin();itr!=deadlistsi.end();itr++ ) {
        adj_to_dead.push_back({(*itr).first-1, (*itr).second});
        adj_to_dead.push_back({(*itr).first+1, (*itr).second});
        n+=1;
    }

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

    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries (events)
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
            if(debug>5) std::cout<<" size of ssvec of weights is "<<(*ssvec).size()<<std::endl;
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
            if(debug>5) std::cout << " -- AbsWeight size: " << absW.size() << std::endl;
            if(debug>5) std::cout<<" values are ";
            for (unsigned iL(0); iL<(*ssvec).size();++iL){
                if(debug>5) std::cout<<" "<<absW[iL];
            }
            if(debug>5) std::cout<<std::endl;
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
        if (debug) std::cout << " - Event contains " << (*rechitvec).size() << " rechits." << std::endl;
        double coneSize = 0.3;
        double rechitsum = 0;
        double rechitsumdead_Si = 0;
        double rechitsumlaypn = 0;
        double rechitsumlayerSum = 0;
        int nlay = 27;
        double layersum[nlay]; //Sum of all rechits in one layers
        for (int iL = 0; iL < nlay; ++iL) {
            layersum[iL] = 0;
        }

        for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
            HGCSSRecoHit lHit = (*rechitvec)[iH];
            unsigned layer = lHit.layer();
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
                rechitsum+=lenergy;
                if(!isScint) {
                    std::pair<unsigned,unsigned> tempsi(layer,cellid);
                    std::set<std::pair<unsigned,unsigned>>::iterator ibc=deadlistsi.find(tempsi);
                    if(ibc==deadlistsi.end()) {
                        rechitsumdead_Si+=lenergy;
                        layersum[layer-1]+=lenergy;
                    }
                    for(unsigned i(0); i < adj_to_dead.size(); i++){
                        if(cellid == adj_to_dead[i].second && layer == adj_to_dead[i].first){
                            rechitsumlaypn += lenergy/2;
                        }
                    }
                }
            }
        }// end loop over hits
        double rechitsumave=rechitsumlaypn+rechitsumdead_Si;
        h_rechitsumave->Fill(rechitsumave);
        h_rechitsum->Fill(rechitsum);
        h_rechitsumdead_Si->Fill(rechitsumdead_Si);
        ievtRec++;
    }//loop on entries

    if(debug) std::cout << "Writing files ..." << std::endl;
    outputFile->cd();
    outputFile->Write();
    outputFile->Close();

    return 0;
}
