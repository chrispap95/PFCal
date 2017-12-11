#include "HiggsMass.hh"

void HiggsMass::setRecoInfo(const TLorentzVector & g1,
			    const TLorentzVector & g2,
			    const ROOT::Math::XYZPoint & posFF1,
			    const ROOT::Math::XYZPoint & posFF2){
  g1_ = g1;
  g2_ = g2;
  posFF1_ = posFF1;
  posFF2_ = posFF2;
}

void HiggsMass::setTruthInfo(const TLorentzVector & g1,
			     const TLorentzVector & g2,
			     const ROOT::Math::XYZPoint & vtx1,
			     const ROOT::Math::XYZPoint & vtx2){

  tg1_ = g1;
  tg2_ = g2;
  tvtx1_ = vtx1;
  tvtx2_ = vtx2;
}


void HiggsMass::initialiseHistograms(TFile *fout, 
				     const std::string folder){

  fout->mkdir(folder.c_str());
  fout->cd(folder.c_str());

  tree_ = new TTree("Photons","higgs photon 4 vectors");
  tree_->Branch("eventIndex",&evtIdx_);
  //Higgs
  tree_->Branch("MH",&MH_);
  tree_->Branch("pTH",&pTH_);
  tree_->Branch("etaH",&etaH_);
  tree_->Branch("truthMH",&truthMH_);
  tree_->Branch("truthpTH",&truthpTH_);
  tree_->Branch("truthetaH",&truthetaH_);

  //photon1
  tree_->Branch("E1",&E1_);
  tree_->Branch("eta1",&eta1_);
  tree_->Branch("phi1",&phi1_);
  tree_->Branch("truthE1",&truthE1_);
  tree_->Branch("trutheta1",&trutheta1_);
  tree_->Branch("truthphi1",&truthphi1_);

  //photon2
  tree_->Branch("E2",&E2_);
  tree_->Branch("eta2",&eta2_);
  tree_->Branch("phi2",&phi2_);
  tree_->Branch("truthE2",&truthE2_);
  tree_->Branch("trutheta2",&trutheta2_);
  tree_->Branch("truthphi2",&truthphi2_);

  p_trueDir_trueE = new TH1F("p_trueDir_trueE",";true M_{H} (GeV)",400,0,200);
  p_trueDir_recoE = new TH1F("p_trueDir_recoE",";M_{H} (GeV) [P_{true},E_{reco}]",400,0,200);
  p_position_trueE = new TH1F("p_position_trueE",";M_{H} (GeV) [P_{pos},E_{true}]",400,0,200);
  p_position_recoE = new TH1F("p_position_recoE",";M_{H} (GeV) [P_{pos},E_{reco}]",400,0,200);
  p_angle_trueE = new TH1F("p_angle_trueE",";M_{H} (GeV) [P_{angle},E_{true}]",400,0,200);
  p_angle_recoE = new TH1F("p_angle_recoE",";M_{H} (GeV) [P_{angle},E_{reco}]",400,0,200);
  p_position_vtxsmear_trueE = new TH1F("p_position_vtxsmear_trueE",";M_{H} (GeV) [P_{pos,vtx smeared},E_{true}]",400,0,200);
  p_position_vtxsmear_recoE = new TH1F("p_position_vtxsmear_recoE",";M_{H} (GeV) [P_{pos,vtx smeared},E_{reco}]",400,0,200);

  p_ErecooverEtrue = new TH1F("p_ErecooverEtrue",";E_{#gamma}/E_{truth};showers",400,0,1.5);

  p_vtx_x = new TH1F("p_vtx_x",";vtx x (mm)",100,-1,1);
  p_vtx_y = new TH1F("p_vtx_y",";vtx y (mm)",100,-1,1);
  p_vtx_z = new TH1F("p_vtx_z",";vtx z (mm)",100,-50,50);

  p_dvtx_x = new TH1F("p_dvtx_x",";dvtx x (mm)",100,-1,1);
  p_dvtx_y = new TH1F("p_dvtx_y",";dvtx y (mm)",100,-1,1);
  p_dvtx_z = new TH1F("p_dvtx_z",";dvtx z (mm)",100,-50,50);

  p_pT_Higgs = new TH1F("p_pT_Higgs",";p_{T}^{H} (GeV)",200,0,400);
  p_eta_Higgs = new TH1F("p_eta_Higgs",";#eta^{H}",150,1.3,6);

  p_pTvseta[0] = new TH2F("p_pTvseta",";#eta^{H};p_{T}^{H} (GeV)",
			  150,1.3,6,
			  200,0,400);
  p_MvspT[0] = new TH2F("p_MvspT",";p_{T}^{H} (GeV);M_{H} (GeV)",
			200,0,400,
			100,105,145);

  p_pT_gamma1[0] = new TH1F("p_pT_gamma1",";p_{T}^{#gamma1}",200,0,300);
  p_eta_gamma1[0] = new TH1F("p_eta_gamma1",";#eta^{#gamma1}",120,1.3,3.5);
  p_phi_gamma1[0] = new TH1F("p_phi_gamma1",";#phi^{#gamma1}",120,-3.1416,3.1416);

  p_pT_gamma2[0] = new TH1F("p_pT_gamma2",";p_{T}^{#gamma2}",200,0,300);
  p_eta_gamma2[0] = new TH1F("p_eta_gamma2",";#eta^{#gamma2}",120,1.3,3.5);
  p_phi_gamma2[0] = new TH1F("p_phi_gamma2",";#phi^{#gamma2}",120,-3.1416,3.1416);

  p_pTvseta[1] = new TH2F("p_pTvseta_truth",";#eta^{H};p_{T}^{H} (GeV)",
			  150,1.3,6,
			  200,0,400);
  p_MvspT[1] = new TH2F("p_MvspT_truth",";p_{T}^{H} (GeV);M_{H} (GeV)",
		     200,0,400,
		     100,105,145);

  p_pT_gamma1[1] = new TH1F("p_pT_gamma1_truth",";p_{T}^{#gamma1}",200,0,300);
  p_eta_gamma1[1] = new TH1F("p_eta_gamma1_truth",";#eta^{#gamma1}",120,1.3,3.5);
  p_phi_gamma1[1] = new TH1F("p_phi_gamma1_truth",";#phi^{#gamma1}",120,-3.1416,3.1416);

  p_pT_gamma2[1] = new TH1F("p_pT_gamma2_truth",";p_{T}^{#gamma2}",200,0,300);
  p_eta_gamma2[1] = new TH1F("p_eta_gamma2_truth",";#eta^{#gamma2}",120,1.3,3.5);
  p_phi_gamma2[1] = new TH1F("p_phi_gamma2_truth",";#phi^{#gamma2}",120,-3.1416,3.1416);

}

void HiggsMass::fillHistograms(){

  TLorentzVector lg1 = g1_.Pt()>g2_.Pt() ? g1_ : g2_;
  TLorentzVector lg2 = g1_.Pt()>g2_.Pt() ? g2_ : g1_;

  p_pT_gamma1[0]->Fill(lg1.Pt());

  p_eta_gamma1[0]->Fill(lg1.Eta());
  p_phi_gamma1[0]->Fill(lg1.Phi());

  p_pT_gamma2[0]->Fill(lg2.Pt());
  p_eta_gamma2[0]->Fill(lg2.Eta());
  p_phi_gamma2[0]->Fill(lg2.Phi());

  TLorentzVector rh = g1_+g2_;
  p_pTvseta[0]->Fill(rh.Eta(),rh.Pt());
  p_MvspT[0]->Fill(rh.Pt(),rh.M());


  p_ErecooverEtrue->Fill(g1_.E()/tg1_.E());
  p_ErecooverEtrue->Fill(g2_.E()/tg2_.E());
  //std::cout << "DD reco1 " << g1_.E() << " t1 " << tg1_.E() << " reco2 " << g2_.E() << " t2 " << tg2_.E() << std::endl;



  TLorentzVector tg1 = tg1_.Pt()>tg2_.Pt() ? tg1_ : tg2_;
  TLorentzVector tg2 = tg1_.Pt()>tg2_.Pt() ? tg2_ : tg1_;

  p_pT_gamma1[1]->Fill(tg1.Pt());
  p_eta_gamma1[1]->Fill(tg1.Eta());
  p_phi_gamma1[1]->Fill(tg1.Phi());

  p_pT_gamma2[1]->Fill(tg2.Pt());
  p_eta_gamma2[1]->Fill(tg2.Eta());
  p_phi_gamma2[1]->Fill(tg2.Phi());

  TLorentzVector th = tg1_+tg2_;

  MH_ = rh.M();
  pTH_ = rh.Pt();
  etaH_ = rh.Eta();
  truthMH_ = th.M();
  truthpTH_ = th.Pt();
  truthetaH_ = th.Eta();

  //keep truth associated to reco!!!
  E1_ = g1_.E();
  eta1_ = g1_.Eta();
  phi1_ = g1_.Phi();
  truthE1_ = tg1_.E();
  trutheta1_ = tg1_.Eta();
  truthphi1_ = tg1_.Phi();

  E2_ = g2_.E();
  eta2_ = g2_.Eta();
  phi2_ = g2_.Phi();
  truthE2_ = tg2_.E();
  trutheta2_ = tg2_.Eta();
  truthphi2_ = tg2_.Phi();

  tree_->Fill();


  p_pTvseta[1]->Fill(th.Eta(),th.Pt());
  p_MvspT[1]->Fill(th.Pt(),th.M());
  p_pT_Higgs->Fill(th.Pt());
  p_eta_Higgs->Fill(th.Eta());

  p_vtx_x->Fill(tvtx1_.x());
  p_vtx_y->Fill(tvtx1_.y());
  p_vtx_z->Fill(tvtx1_.z());

  p_dvtx_x->Fill(tvtx2_.x()-tvtx1_.x());
  p_dvtx_y->Fill(tvtx2_.y()-tvtx1_.y());
  p_dvtx_z->Fill(tvtx2_.z()-tvtx1_.z());

  //true
  p_trueDir_trueE->Fill(th.M());

  //true dir, reco E
  lg1.SetPtEtaPhiE(tg1_.Pt()*g1_.E()/tg1_.E(),tg1_.Eta(),tg1_.Phi(),g1_.E());
  lg2.SetPtEtaPhiE(tg2_.Pt()*g2_.E()/tg2_.E(),tg2_.Eta(),tg2_.Phi(),g2_.E());  
  TLorentzVector lh = lg1+lg2;
  p_trueDir_recoE->Fill(lh.M());

  //pos+true vtx, true E
  lg1 = tg1_;
  double dx = posFF1_.x()-tvtx1_.x();
  double dy = posFF1_.y()-tvtx1_.y();
  double dz = posFF1_.z()-tvtx1_.z();
  double norm = sqrt(dx*dx+dy*dy+dz*dz);

  TVector3 vtxFF1 = TVector3(dx/norm,dy/norm,dz/norm)*tg1_.E();
  lg1.SetVect(vtxFF1);
  lg2 = tg2_;
  dx = posFF2_.x()-tvtx2_.x();
  dy = posFF2_.y()-tvtx2_.y();
  dz = posFF2_.z()-tvtx2_.z();
  norm = sqrt(dx*dx+dy*dy+dz*dz);
  TVector3 vtxFF2 = TVector3(dx/norm,dy/norm,dz/norm)*tg2_.E();
  lg2.SetVect(vtxFF2);
  lh = lg1+lg2;
  p_position_trueE->Fill(lh.M());

  //pos+true vtx, reco E
  lg1.SetE(g1_.E());
  lg1.SetVect(vtxFF1*(g1_.E()/tg1_.E()));
  lg2.SetE(g2_.E());
  lg2.SetVect(vtxFF2*(g2_.E()/tg2_.E()));
  lh = lg1+lg2;
  p_position_recoE->Fill(lh.M());

  //reco angle recoE
  lg1 = g1_;
  lg2 = g2_;
  lh = lg1+lg2;
  p_angle_recoE->Fill(lh.M());

  //reco angle true E
  lg1.SetPtEtaPhiE(g1_.Pt()*tg1_.E()/g1_.E(),g1_.Eta(),g1_.Phi(),tg1_.E());
  lg2.SetPtEtaPhiE(g2_.Pt()*tg2_.E()/g2_.E(),g2_.Eta(),g2_.Phi(),tg2_.E());
  lh = lg1+lg2;
  p_angle_trueE->Fill(lh.M());

  //pos+smeared vtx, true E
  lg1 = tg1_;
  dx = posFF1_.x()-tvtx1_.x();
  dy = posFF1_.y()-tvtx1_.y();
  double vtxz = rand_.Gaus(tvtx1_.z(),50);
  dz = posFF1_.z()- vtxz;
  norm = sqrt(dx*dx+dy*dy+dz*dz);

  vtxFF1 = TVector3(dx/norm,dy/norm,dz/norm)*tg1_.E();
  lg1.SetVect(vtxFF1);
  lg2 = tg2_;
  dx = posFF2_.x()-tvtx2_.x();
  dy = posFF2_.y()-tvtx2_.y();
  dz = posFF2_.z()- vtxz;
  norm = sqrt(dx*dx+dy*dy+dz*dz);
  vtxFF2 = TVector3(dx/norm,dy/norm,dz/norm)*tg2_.E();
  lg2.SetVect(vtxFF2);
  lh = lg1+lg2;
  p_position_vtxsmear_trueE->Fill(lh.M());

  //pos+smeared vtx, reco E
  lg1.SetE(g1_.E());
  lg1.SetVect(vtxFF1*(g1_.E()/tg1_.E()));
  lg2.SetE(g2_.E());
  lg2.SetVect(vtxFF2*(g2_.E()/tg2_.E()));
  lh = lg1+lg2;
  p_position_vtxsmear_recoE->Fill(lh.M());


}
