{
    // Import files
    //TFile* f1 = TFile::Open("rootFiles/et5_eta1.7_test_deadfrac070.root");
    TFile* f2 = TFile::Open("rootFiles/et10_eta1.7_test_deadfrac050_corr_noadj.root");
    TFile* f3 = TFile::Open("rootFiles/et20_eta1.7_test_deadfrac050_corr_noadj.root");
    TFile* f4 = TFile::Open("rootFiles/et30_eta1.7_test_deadfrac050_corr_noadj.root");
    TFile* f5 = TFile::Open("rootFiles/et40_eta1.7_test_deadfrac050_corr_noadj.root");
    TFile* f6 = TFile::Open("rootFiles/et60_eta1.7_test_deadfrac050_corr_noadj.root");
    //TFile* f7 = TFile::Open("rootFiles/et80_eta1.7_test_deadfrac070.root");
    TFile* f8 = TFile::Open("rootFiles/et100_eta1.7_test_deadfrac050_corr_noadj.root");
    //TFile* f9 = TFile::Open("rootFiles/et200_eta1.7_test_deadfrac070.root");

    // Import histograms
    //TH1F* h1 = (TH1F*)f1->Get("h_rechitsumave");
    TH1F* h2 = (TH1F*)f2->Get("h_rechitsumdead_Si");
    TH1F* h3 = (TH1F*)f3->Get("h_rechitsumdead_Si");
    TH1F* h4 = (TH1F*)f4->Get("h_rechitsumdead_Si");
    TH1F* h5 = (TH1F*)f5->Get("h_rechitsumdead_Si");
    TH1F* h6 = (TH1F*)f6->Get("h_rechitsumdead_Si");
    //TH1F* h7 = (TH1F*)f7->Get("h_rechitsumdead_Siave");
    TH1F* h8 = (TH1F*)f8->Get("h_rechitsumdead_Si");
    //TH1F* h9 = (TH1F*)f9->Get("h_rechitsumave");

    // Define fit functions
    TF1* tstud = new TF1("tstud","[0]*sqrt([3])*TMath::Gamma(([2]+1)/2)*(1+[3]*(x-[1])^2/[2])**(-0.5*([2]+1))/(sqrt(TMath::Pi()*[2])*TMath::Gamma(0.5*[2]))",20,55);
    tstud->SetParNames("Scale","#mu","#nu","#lambda");
    tstud->SetParameters(100,30,1,1);
    tstud->SetParLimits(0,0.1,15000);
    tstud->SetParLimits(1,2,50);
    tstud->SetParLimits(2,0.1,10);
    // Define fit functions
    TF1* tstud2 = new TF1("tstud2","[0]*sqrt([3])*TMath::Gamma(([2]+1)/2)*(1+([3]*x-[1])^2/[2])**(-0.5*([2]+1))/(sqrt(TMath::Pi()*[2])*TMath::Gamma(0.5*[2]))",40,85);
    tstud2->SetParNames("Scale","#mu","#nu","#lambda");
    tstud2->SetParameters(100,55,1,1);
    tstud2->SetParLimits(0,0.1,15000);
    tstud2->SetParLimits(1,10,60);
    tstud2->SetParLimits(2,0.1,10);
    // Define fit functions
    TF1* tstud3 = new TF1("tstud3","[0]*sqrt([3])*TMath::Gamma(([2]+1)/2)*(1+[3]*(x-[1])^2/[2])**(-0.5*([2]+1))/(sqrt(TMath::Pi()*[2])*TMath::Gamma(0.5*[2]))",75,120);
    tstud3->SetParNames("Scale","#mu","#nu","#lambda");
    tstud3->SetParameters(100,80,1,1);
    tstud3->SetParLimits(0,0.1,15000);
    tstud3->SetParLimits(1,50,90);
    tstud3->SetParLimits(2,0.1,10);
    // Define fit functions
    TF1* tstud4 = new TF1("tstud4","[0]*sqrt([3])*TMath::Gamma(([2]+1)/2)*(1+[3]*(x-[1])^2/[2])**(-0.5*([2]+1))/(sqrt(TMath::Pi()*[2])*TMath::Gamma(0.5*[2]))",155,220);
    tstud4->SetParNames("Scale","#mu","#nu","#lambda");
    tstud4->SetParameters(100,160,1,1);
    tstud4->SetParLimits(0,0.1,15000);
    tstud4->SetParLimits(1,140,180);
    tstud4->SetParLimits(2,0.1,10);
    // Define fit functions
    TF1* tstud5 = new TF1("tstud5","[0]*sqrt([3])*TMath::Gamma(([2]+1)/2)*(1+[3]*(x-[1])^2/[2])**(-0.5*([2]+1))/(sqrt(TMath::Pi()*[2])*TMath::Gamma(0.5*[2]))",260,400);
    tstud5->SetParNames("Scale","#mu","#nu","#lambda");
    tstud5->SetParameters(100,280,1,1);
    tstud5->SetParLimits(0,0.1,15000);
    tstud5->SetParLimits(1,200,300);
    tstud5->SetParLimits(2,0.1,10);


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
    tstud->SetNpx(700);
    tstud2->SetNpx(700);
    tstud3->SetNpx(700);
    tstud4->SetNpx(700);
    tstud5->SetNpx(700);

    //h2->Rebin();
    //h3->Rebin();
    //h4->Rebin();
    //h5->Rebin();
    //h6->Rebin();
    h6->Rebin();
    h8->Rebin(4);

    gStyle->SetOptStat(00000000);
    gStyle->SetOptFit(1111);
    TCanvas* c = new TCanvas("c","c",1000,600);
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
    //h1->GetXaxis()->SetRangeUser(0,30);
    //h1->Draw(); //h1->SetTitle("rechit Energy 5 GeV");
    //h1->Fit("g1","");
    //h1->Fit("g1","","",g1->GetParameter(1)-3*g1->GetParameter(2),g1->GetParameter(1)+3*g1->GetParameter(2));
    //c->cd(2);
    h2->GetXaxis()->SetRangeUser(10,55);
    h2->Draw(); h2->SetTitle("rechit Energy 10 GeV ;Energy [GeV];Hits");
    h2->Fit("tstud","","",20,35);
    //h2->Fit("g2","","",g2->GetParameter(1)-3*g2->GetParameter(2),g2->GetParameter(1)+3*g2->GetParameter(2));
    c->cd(2);
    h3->GetXaxis()->SetRangeUser(30,85);
    h3->Draw(); h3->SetTitle("rechit Energy 20 GeV;Energy [GeV];Hits");
    h3->Fit("tstud2","","",40,70);
    //h3->Fit("g3","","",g3->GetParameter(1)-2*g3->GetParameter(2),g3->GetParameter(1)+3*g3->GetParameter(2));
    c->cd(3);
    h4->GetXaxis()->SetRangeUser(50,120);
    h4->Draw(); h4->SetTitle("rechit Energy 30 GeV;Hits;Energy [GeV]");
    h4->Fit("tstud3","","",77,100);
    //h4->Fit("g4","","",g4->GetParameter(1)-1.5*g4->GetParameter(2),g4->GetParameter(1)+3*g4->GetParameter(2));
    c->cd(4);
    //h5->GetXaxis()->SetRangeUser(60,160);
    //h5->Draw(); h5->SetTitle("rechit Energy 40 GeV");
    //h5->Fit("g5","");
    //h5->Fit("g5","","",g5->GetParameter(1)-1.5*g5->GetParameter(2),g5->GetParameter(1)+3*g5->GetParameter(2));
    //c->cd(5);
    h6->GetXaxis()->SetRangeUser(120,220);
    h6->Draw(); h6->SetTitle("rechit Energy 60 GeV;Energy [GeV];Hits");
    h6->Fit("tstud4","","",157,190);
    //h6->Fit("g6","","",g6->GetParameter(1)-1*g6->GetParameter(2),g6->GetParameter(1)+3*g6->GetParameter(2));
    //c->cd(5);
    //h7->GetXaxis()->SetRangeUser(180,280);
    //h7->Draw(); //h7->SetTitle("rechit Energy 80 GeV");
    //h7->Fit("g7","");
    //h7->Fit("g7","","",g7->GetParameter(1)-1.5*g7->GetParameter(2),g7->GetParameter(1)+1.5*g7->GetParameter(2));
    c->cd(5);
    h8->GetXaxis()->SetRangeUser(150,400);
    h8->Draw(); h8->SetTitle("rechit Energy 100 GeV;Energy [GeV];Hits");
    h8->Fit("tstud5","","",265,310);
    //h8->Fit("g8","","",g8->GetParameter(1)-1*g8->GetParameter(2),g8->GetParameter(1)+3*g8->GetParameter(2));
    //c->cd(9);
    //h9->GetXaxis()->SetRangeUser(480,660);
    //h9->Draw(); //h9->SetTitle("rechit Energy 200 GeV");
    //h9->Fit("g9","");
    //h9->Fit("g9","","",g9->GetParameter(1)-2*g9->GetParameter(2),g9->GetParameter(1)+2*g9->GetParameter(2));
    if (setlog) c->Print("plots/rechitsum_log.eps");
    //else c->Print("plots/rechitsum_deadfrac050_av.eps");
}
