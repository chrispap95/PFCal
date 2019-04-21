{
    TRandom3 r(0);
    TH1F* h1 = new TH1F("h1","toy mc", 100, 0, 40);

    for (int i = 0; i < 1000000; ++i){
        double x = 40*r.Rndm();
        double u = r.Rndm();
        double w = r.Rndm();
        double a0 = 176.853;
        double a1 = 6025.48;
        double a2 = -177001;
        double a3 = -8.7379;
        double a4 = -5361;
        double a5 = 3481.65;
        double y = a0/x+a1/pow(x,2)+a4/pow(x,3)+a5/pow(x,4)+a3;

        if (y/97068 > u) {
            if (1) h1->Fill(x);
        }
    }
    h1->Draw();
}//28752
