{
    // Import files
    TFile* f2 = TFile::Open("et60_eta1.7_test_deadfrac010.root");
    TFile* f3 = TFile::Open("et60_eta1.7_test_deadfrac020.root");
    TFile* f4 = TFile::Open("et60_eta1.7_test_deadfrac030.root");
    TFile* f6 = TFile::Open("et60_eta1.7_test_deadfrac050.root");
    TFile* f8 = TFile::Open("et60_eta1.7_test_deadfrac070.root");

    // Import histograms
    TH1F* h1 = (TH1F*)f2->Get("h_rechitsum");
    TH1F* h2 = (TH1F*)f2->Get("h_rechitsumdead_Si");
    TH1F* h3 = (TH1F*)f3->Get("h_rechitsumdead_Si");
    TH1F* h4 = (TH1F*)f4->Get("h_rechitsumdead_Si");
    TH1F* h6 = (TH1F*)f6->Get("h_rechitsumdead_Si");
    TH1F* h8 = (TH1F*)f8->Get("h_rechitsumdead_Si");

    // Define fit functions
    TF1* g1 = new TF1("g1","gaus");
    TF1* g2 = new TF1("g2","gaus");
    TF1* g3 = new TF1("g3","gaus");
    TF1* g4 = new TF1("g4","gaus");
    TF1* g5 = new TF1("g5","gaus");
    TF1* g6 = new TF1("g6","gaus");
    TF1* g7 = new TF1("g7","gaus");
    TF1* g8 = new TF1("g8","gaus");
    TF1* g9 = new TF1("g9","gaus");

    g1->SetNpx(200);
    g2->SetNpx(200);
    g3->SetNpx(200);
    g4->SetNpx(200);
    g5->SetNpx(200);
    g6->SetNpx(200);
    g7->SetNpx(200);
    g8->SetNpx(200);
    g9->SetNpx(200);

    h1->Rebin();
    h2->Rebin();
    h3->Rebin();
    h4->Rebin();
    h8->Rebin();
    h6->Rebin();


    gStyle->SetOptStat(00000000);
    gStyle->SetOptFit(1111);
    TCanvas* c = new TCanvas("c","c",1200,600);
    Bool_t setlog = false;
    c->Divide(3,2);
    if (setlog) {
        c->cd(1)->SetLogy();
        c->cd(2)->SetLogy();
        c->cd(3)->SetLogy();
        c->cd(4)->SetLogy();
        c->cd(5)->SetLogy();
        c->cd(6)->SetLogy();
        c->cd(7)->SetLogy();
        c->cd(8)->SetLogy();
        c->cd(9)->SetLogy();
    }
    c->cd(1);
    h1->GetXaxis()->SetRangeUser(120,220);
    h1->Draw(); h1->SetTitle("rechit Energy 60 GeV for 0%;Energy [GeV];Hits");
    h1->Fit("g1","");
    h1->Fit("g1","","",g1->GetParameter(1)-2.5*g1->GetParameter(2),g1->GetParameter(1)+3*g1->GetParameter(2));
    c->cd(2);
    h2->GetXaxis()->SetRangeUser(120,220);
    h2->Draw(); h2->SetTitle("rechit Energy 60 GeV for 1% ;Energy [GeV];Hits");
    h2->Fit("g2","");
    h2->Fit("g2","","",g2->GetParameter(1)-2*g2->GetParameter(2),g2->GetParameter(1)+3*g2->GetParameter(2));
    c->cd(3);
    h3->GetXaxis()->SetRangeUser(120,220);
    h3->Draw(); h3->SetTitle("rechit Energy 60 GeV for 2%;Energy [GeV];Hits");
    h3->Fit("g3","");
    h3->Fit("g3","","",g3->GetParameter(1)-1.3*g3->GetParameter(2),g3->GetParameter(1)+3*g3->GetParameter(2));
    c->cd(4);
    h4->GetXaxis()->SetRangeUser(120,220);
    h4->Draw(); h4->SetTitle("rechit Energy 60 GeV for 3%;Hits;Energy [GeV]");
    h4->Fit("g4","");
    h4->Fit("g4","","",g4->GetParameter(1)-1.1*g4->GetParameter(2),g4->GetParameter(1)+3*g4->GetParameter(2));
    c->cd(5);
    h6->GetXaxis()->SetRangeUser(120,220);
    h6->Draw(); h6->SetTitle("rechit Energy 60 GeV for 5%;Energy [GeV];Hits");
    h6->Fit("g6","");
    h6->Fit("g6","","",g6->GetParameter(1)-0.7*g6->GetParameter(2),g6->GetParameter(1)+3*g6->GetParameter(2));
    c->cd(6);
    h8->GetXaxis()->SetRangeUser(120,220);
    h8->Draw(); h8->SetTitle("rechit Energy 60 GeV for 7%;Energy [GeV];Hits");
    h8->Fit("g8","");
    h8->Fit("g8","","",g8->GetParameter(1)-0.7*g8->GetParameter(2),g8->GetParameter(1)+3*g8->GetParameter(2));
    if (setlog) c->Print("rechitsum_log.eps");
    //else c->Print("presentation/rechitsum_deadfrac050_av.eps");
    double gy[6] = {g1->GetParameter(2)/g1->GetParameter(1),g2->GetParameter(2)/g2->GetParameter(1),g3->GetParameter(2)/g3->GetParameter(1),
        g4->GetParameter(2)/g4->GetParameter(1),g6->GetParameter(2)/g6->GetParameter(1),g8->GetParameter(2)/g8->GetParameter(1)
    };
    double gx[6] = {0,0.01,0.02,0.03,0.05,0.07};
    double egx[6] = {0,0,0,0,0,0};
    double egy[6] = {gy[0]*sqrt(pow(g1->GetParError(1)/g1->GetParameter(1),2)+pow(g1->GetParError(2)/g1->GetParameter(2),2)),
        gy[1]*sqrt(pow(g2->GetParError(1)/g2->GetParameter(1),2)+pow(g2->GetParError(2)/g2->GetParameter(2),2)),
        gy[2]*sqrt(pow(g3->GetParError(1)/g3->GetParameter(1),2)+pow(g3->GetParError(2)/g3->GetParameter(2),2)),
        gy[3]*sqrt(pow(g4->GetParError(1)/g4->GetParameter(1),2)+pow(g4->GetParError(2)/g4->GetParameter(2),2)),
        gy[4]*sqrt(pow(g6->GetParError(1)/g6->GetParameter(1),2)+pow(g6->GetParError(2)/g6->GetParameter(2),2)),
        gy[5]*sqrt(pow(g8->GetParError(1)/g8->GetParameter(1),2)+pow(g8->GetParError(2)/g8->GetParameter(2),2))
    };

    TGraphErrors *gr = new TGraphErrors(6,gx,gy,egx,egy);
    gr->SetTitle("Resolution vs deadfrac for 60 GeV");
    gr->GetYaxis()->SetTitle("#frac{width}{mean}");
    gr->GetXaxis()->SetTitle("deadfrac");
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    TCanvas* c1 = new TCanvas("c1","c1",1);
    gr->Draw("AP");
    TFile* out = new TFile("resvsdeadfrac.root","RECREATE");
    gr->Write();
    out->Close();
}
