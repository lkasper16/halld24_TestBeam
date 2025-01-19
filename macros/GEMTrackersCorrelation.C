
vector<Double_t> GetLinearFitParameters(TH2F *h2, Double_t xmin, Double_t xmax, Double_t step);

void GEMTrackersCorrelation(TString rootFile){
    TFile *f = TFile::Open(rootFile);
    TString h2name_gemtrkr_x = "hgemtrkr_double_x";

    TH2F *h2_gemtrkr_x = (TH2F*)f->Get(h2name_gemtrkr_x);

    // Get the fit parameters
    vector<Double_t> fitParams = GetLinearFitParameters(h2_gemtrkr_x, 100, 300, 1);
    cout << "Fit parameters: " << fitParams[0] << " " << fitParams[1] << endl;
}

vector<Double_t> GetLinearFitParameters(TH2F *h2, Double_t xmin, Double_t xmax, Double_t step) {
    TGraph *g = new TGraph();
    vector<Double_t> fitParams;
    for (Double_t x = xmin; x < xmax; x += step) {
        TH1D *h1 = h2->ProjectionY("h1", x, x+step);
        Double_t X2 = h1->GetBinCenter(h1->GetMaximumBin());
        Double_t X1 = x + step/2;
        g->SetPoint(g->GetN(), X1, X2);
    }

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    g->Draw("AP");

    // Fit g with a linear function
    TF1 *f1 = new TF1("f1", "[0] + [1]*x", xmin, xmax);
    g->Fit(f1, "R");

    // Get the fit parameters
    fitParams.push_back(f1->GetParameter(0));
    fitParams.push_back(f1->GetParameter(1));

    // g->SaveAs("GEMTrackersCorrelation.pdf");
    return fitParams;   
}