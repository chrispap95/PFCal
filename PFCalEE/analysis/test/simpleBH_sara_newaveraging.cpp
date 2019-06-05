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

std::set<int> neighborhood_list(int cellint, int n) {
    /* Gives a set of the cell IDs of a neighborhood.
    ** Parameter n gives the radius that is included eg. n=1 includes the nearest
    ** neighbors, n=2 is nearest neighbors + neighbors' neighbors, etc.
    */

    std::set<int> v;
    v.insert(cellint);

    for (int i(1); i <= n; ++i) {
        for (auto itr = v.begin(); itr != v.end(); itr++) {
            v.insert((*itr)-497);
            v.insert((*itr)-496);
            v.insert((*itr)-1);
            v.insert((*itr)+1);
            v.insert((*itr)+496);
            v.insert((*itr)+497);
        }
    }
    return v;
}

int main(int argc, char** argv){
    /**********************************
    ** initialize some variables
    **********************************/
    const unsigned nscintlayer=16;
    const unsigned scintoffset=36;
    const unsigned nsilayer=52;

    unsigned siminid[nsilayer]={34987,34497,33997,33510,33007,33001,32505,32011,31515,31024,30531,30518,30025,29528,29035,28548,28045,
        28036,27542,27049,26552,26062,25562,25556,25057,24566,24066,23576,22080,20105,18607,17114,15136,13643,12150,10177,34987,34497,
        33507,33007,33001,32505,32011,31515,31024,30531,30518,30025,29528,29035,28548
    };
    unsigned simaxid[nsilayer]={250502,250992,251492,251979,2522479,252488,252982,253475,253972,254465,254955,254971,255464,255958,
        256454,256941,257441,257451,257947,258437,258937,259427,259927,259930,260430,260920,261420,261910,263406,265382,266882,268375,
        270347,271844,273337,275316,250502,250992,251979,252479,252488,252982,253475,253972,254465,254955,254971,255464,255958,256454,
        256941
    };

    //Input output and config options
    std::string cfg;
    unsigned pNevts;
    std::string outFilePath;
    std::string filePath;
    std::string digifilePath;
    unsigned nRuns;
    std::string simFileName;
    std::string recoFileName;
    std::string AvFilePath;
    unsigned debug;
    double deadfrac;
    bool adjacent, Avsample;
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
    ("Avsample",        po::value<bool>(&Avsample)->default_value(1)) //Generate ML study training sample
    ("AvFilePath",      po::value<std::string>(&AvFilePath)->default_value("av_sample.root")) //File to export data for ML
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
        std::cout << " -- \e[31mError\e[0m, output file " << outFilePath
        << " cannot be opened. Please create output directory. Exiting..." << std::endl;
        return 1;
    }
    else std::cout << " -- opening output file: " << outputFile->GetName() << " \e[32msuccess\e[0m" << std::endl;

    outputFile->cd();

    TH1F* h_test = new TH1F("h_test","h_test",100,0,40);
    TH1F* h_rechitsum = new TH1F("h_rechitsum","Rechitsum silicon;Enegry [GeV]",700,0,700.);
    TH1F* h_rechitsumdead_Si = new TH1F("h_rechitsumdead_Si","Rechitsum dead silicon; Energy [GeV]",700,0,700.);
    TH1F* h_rechitsumave = new TH1F("h_rechitsumave","Sum energy averaging method 2.0; Energy [GeV]",700,0,700.);
    TH2F* h_avbias = new TH2F("h_avbias","Averaging 2.0 method bias;Rechit_{true} [GeV];(Rechit_{true}-Rechit_{av})",
    100,0,40,100,-40,40);

    /*
    ** Define a vector of the array:
    ** {dead cell(dc) layer, dc id, dc eta, dc phi, MLn1, MLn2, MLn3, MLn4, MLn5, MLn6, dc rechit}
    */
    std::array<float, 5> avRechits;
    std::vector<std::array<float, 5>> avRechitsVector;

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

    //const unsigned nDeadSiliconCell = 13983553*deadfrac;

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
        if (deadlistsi.find(std::make_pair(ld_si, cd_si)) != deadlistsi.end()) {
            --i;
            continue;
        }
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
            if (adj_ok < 2) deadlistsi.insert(std::make_pair(ld_si,cd_si));
            if (adj_ok == 1) N_clusters++;
        }
        else {
            deadlistsi.insert(std::make_pair(ld_si,cd_si));
            avRechits[0] = ld_si;
            avRechits[1] = cd_si;
            avRechitsVector.push_back(avRechits);
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

    /**********************************
    **  start event loop
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

    std::array<float, 27> Trarr;
    std::array<float, 27> Trdiff;

    float avdevquad = 0;
    int ndev = 0;

    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries (events)
        for(std::vector<std::array<float, 5>>::iterator itr = avRechitsVector.begin(); itr != avRechitsVector.end(); ++itr){
            (*itr)[2]=0;
            (*itr)[3]=0;
            (*itr)[4]=0;
        }
        if (ievtRec>=lRecTree->GetEntries()) continue;

        mycalib.setVertex(0.,0.,0.);

        if (debug) std::cout << std::endl<<std::endl<<"... Processing entry: " << ievt << std::endl;
        else if (ievt%50 == 0) std::cout << "... Processing entry: " << ievt << std::endl;

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
        double rechitsumdead_Si = 0;
        double rechitsumlaypn = 0;

        const unsigned nlay=27;

        /*
        ** Now we need to populate a list with the predicted cellids for each
        ** layer. This accomplished through the use of TrackId() for each layer.
        */
        TH2Poly* map2 = geomConv.hexagonMap();
        std::set<std::pair<unsigned, unsigned>> tracklist;

        for(unsigned k(0); k < 27; ++k) Trarr[k] = 0;
        for (unsigned iL(1); iL <= nlay; ++iL){
            std::pair<double, double> xypair = TrackId(iL, thetagen, phigen);
            tracklist.insert(std::make_pair(iL,map2->FindBin(xypair.first, xypair.second)));
            Trarr[iL-1] = map2->FindBin(xypair.first, xypair.second);
            if(iL == 1) Trdiff[iL-1] = 0;
            else Trdiff[iL-1] = Trarr[iL-1]-Trarr[iL-2];
        }

        for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop on hits
            HGCSSRecoHit lHit = (*rechitvec)[iH];
            unsigned layer = lHit.layer();
            unsigned ixx=layer;
            double xh=lHit.get_x();
            double yh=lHit.get_y();
            double zh=lHit.get_z();
            double lenergy=lHit.energy()*absW[layer]/1000;
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
                    std::pair<unsigned,unsigned> tempsi(layer,cellid);
                    std::set<std::pair<unsigned,unsigned>>::iterator ibc=deadlistsi.find(tempsi);
                    if(ibc==deadlistsi.end()) {
                        rechitsumdead_Si+=lenergy;
                    }else {
                        for(std::vector<std::array<float, 5>>::iterator itr = avRechitsVector.begin(); itr != avRechitsVector.end(); ++itr){
                            if(cellid == (*itr)[1] && layer == (*itr)[0]) {
                                (*itr)[3] = lenergy;
                            }
                        }
                    }
                    for(std::vector<std::array<float, 5>>::iterator itr = avRechitsVector.begin(); itr != avRechitsVector.end(); ++itr){
                        //Loop over dead cells
                        //Different layer averaging
                        if(cellid == (*itr)[1]+Trdiff[layer-1] && layer == (*itr)[0]+1) {
                            rechitsumlaypn += lenergy/2;
                            (*itr)[4] = lenergy;
                        }else if(cellid == (*itr)[1]+Trdiff[layer-1] && layer == (*itr)[0]-1) {
                            rechitsumlaypn += lenergy/2;
                            (*itr)[2] = lenergy;
                        }
                    }
                }
            }
        }// end loop over hits
        for(std::vector<std::array<float, 5>>::iterator itr = avRechitsVector.begin(); itr != avRechitsVector.end(); ++itr) {
            float avPrediction = (*itr)[2]/2+(*itr)[4]/2;
            if((*itr)[3] && avPrediction){
            h_avbias->Fill((*itr)[3],(*itr)[3]-avPrediction);
            ndev++;
            avdevquad+=pow((*itr)[3]-avPrediction,2);
            h_test->Fill((*itr)[3]);
        }
        }

        double rechitsumave = rechitsumlaypn+rechitsumdead_Si;
        h_rechitsumave->Fill(rechitsumave);
        h_rechitsum->Fill(rechitsum);
        h_rechitsumdead_Si->Fill(rechitsumdead_Si);
        ievtRec++;
    }//loop on entries

    std::cout << "average quadratic deviation = " << sqrt(avdevquad/(float)ndev) << std::endl;
    std::cout << "ndev = " << ndev << std::endl;
    std::cout << "avdevquad = " << avdevquad << std::endl;
    std::cout << "Writing files ..." << std::endl;
    outputFile->cd();
    outputFile->Write();
    outputFile->Close();

    return 0;
}
