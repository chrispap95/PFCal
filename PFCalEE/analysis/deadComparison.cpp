{
    TCanvas* c = new TCanvas("c","canvas",1);
    gStyle->SetOptFit(0);
    TFile* f1 = TFile::Open("out00.root");
    TFile* f2 = TFile::Open("out50.root");
    TFile* f3 = TFile::Open("out50noadj.root");
    TFile* f4 = TFile::Open("out50corr_av.root");
    TFile* f5 = TFile::Open("out50corr_noadj_av.root");
    //TFile* f6 = TFile::Open("resvsdeadfrac_ave.root");
    //TFile* f7 = TFile::Open("out50_noadj.root");

    TGraphErrors* gr1 = (TGraphErrors*)f1->Get("Graph");
    TGraphErrors* gr2 = (TGraphErrors*)f2->Get("Graph");
    TGraphErrors* gr3 = (TGraphErrors*)f3->Get("Graph");
    TGraphErrors* gr4 = (TGraphErrors*)f4->Get("Graph");
    TGraphErrors* gr5 = (TGraphErrors*)f5->Get("Graph");
    //TGraphErrors* gr6 = (TGraphErrors*)f6->Get("Graph");
    //TGraphErrors* gr7 = (TGraphErrors*)f7->Get("Graph");
    TMultiGraph* g = new TMultiGraph();

    gr1->SetTitle("deadfrac = 0%");
    gr2->SetTitle("deadfrac = 5%");
    gr3->SetTitle("deadfrac = 5% no adjacent");
    gr4->SetTitle("deadfrac = 5% averaging corrected");
    gr5->SetTitle("deadfrac = 5% noadj aver. corrected");
    //gr6->SetTitle("E_{t}=60GeV averaging");
    //gr7->SetTitle("deadfrac = 5% no adjacent");

    gr1->SetMarkerColor(kBlack);
    gr2->SetMarkerColor(kBlue);
    gr3->SetMarkerColor(kMagenta);
    gr4->SetMarkerColor(kRed);
    gr5->SetMarkerColor(kYellow);
    //gr6->SetMarkerColor(kGreen);
    //gr7->SetMarkerColor(kCyan);

    gr1->GetFunction("f2")->SetLineColor(kBlack);
    gr2->GetFunction("f2")->SetLineColor(kBlue);
    gr3->GetFunction("f2")->SetLineColor(kMagenta);
    gr4->GetFunction("f2")->SetLineColor(kRed);
    gr5->GetFunction("f2")->SetLineColor(kYellow);
    //gr6->GetFunction("f2")->SetLineColor(kGreen);
    //gr7->GetFunction("f2")->SetLineColor(kCyan);

    g->Add(gr1,"P");
    g->Add(gr2,"P");
    g->Add(gr3,"P");
    g->Add(gr4,"P");
    g->Add(gr5,"P");
    //g->Add(gr6,"P");
    //g->Add(gr2,"P");
    //g->Add(gr7,"P");
    g->SetTitle("Resolution vs Energy");
    g->Draw("AP");
    c->BuildLegend();
    g->GetXaxis()->SetTitle("E [GeV]");
    g->GetYaxis()->SetTitle("#frac{width}{mean}");
    gPad->Update();
    //c->Print("presentation/deadfrac_comparison05_noadj.eps")
}
