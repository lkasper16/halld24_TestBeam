
const TString tree_name = "gem_trackers_hits";
const TString RootOutputPath = "/work/halld2/home/nseptian/halld24_TestBeam/RootOutput/";
// no radiator min,max = 0, 150
// hermes radiator min,max = 150, 267
// VU foil radiator min,max = 267, 390
// Fleece radiator min,max = 390, 511

const vector<TString> radiator_names = {"No Radiator", "Hermes Radiator", "VU Foil Radiator", "Fleece Radiator"};

const vector<pair<Int_t,Int_t>> ychan_limits = {
    {0, 150},
    {150, 267},
    {267, 390},
    {390, 527}
};

void PrintRadiator1DHistogram(vector<ROOT::RDF::RResultPtr<TH1D>> v_h1, TString pdfName){
    TCanvas *c = new TCanvas("c", "c", 1600, 1200);
    c->Divide(2, 2);
    for (size_t i = 0; i < v_h1.size(); ++i) {
        c->cd(i + 1);
        v_h1[i]->Draw();
    }
    c->SaveAs(pdfName);
    delete c;
};

void SetYChAxisRadiator(vector<ROOT::RDF::RResultPtr<TH1D>> v_h1){
    for (size_t i = 0; i < v_h1.size(); ++i) {
        v_h1[i]->GetXaxis()->SetRangeUser(ychan_limits[i].first, ychan_limits[i].second);
        v_h1[i]->GetXaxis()->SetTitle("GEM-TRD Channel (X)");
    }
};

void SetDeltaXAxisRadiator(vector<ROOT::RDF::RResultPtr<TH1D>> v_h1){
    for (size_t i = 0; i < v_h1.size(); ++i) {
        v_h1[i]->GetXaxis()->SetTitle("GEM Trackers #DeltaX (mm)");
        v_h1[i]->GetYaxis()->SetTitle("Counts");
    }
};

void PrintRadiator2DHistogram(vector<ROOT::RDF::RResultPtr<TH2D>> v_h2, TString pdfName){
    TCanvas *c = new TCanvas("c", "c", 1600, 1200);
    c->Divide(2, 2);
    for (size_t i = 0; i < v_h2.size(); ++i) {
        c->cd(i + 1);
        v_h2[i]->Draw("COLZ");
    }
    c->SaveAs(pdfName);
    delete c;
};

void PrintEfficiencyHistogram(vector<ROOT::RDF::RResultPtr<TH1D>> v_h1, vector<ROOT::RDF::RResultPtr<TH1D>> v_h1_true, TString pdfName){
    TCanvas *c = new TCanvas("c", "c", 1600, 1200);
    TH1D *h1[4];
    c->Divide(2, 2);
    for (size_t i = 0; i < v_h1.size(); ++i) {
        c->cd(i + 1);
        h1[i] = (TH1D*)v_h1_true[i]->Clone(Form("%s_clone", v_h1[i]->GetName()));
        h1[i]->Divide(v_h1[i].GetPtr());
        h1[i]->Draw();
    }
    c->SaveAs(pdfName);
    delete c;
    for (size_t i = 0; i < v_h1.size(); ++i) delete h1[i];
};

Bool_t GetIsMatchHit(Double_t extrp_y, vector<Double_t> v_trd_y, vector<bool> v_trk_ytime_coincidence){
    const Double_t trd_y_diff_threshold = 20; // same as in trd class, can be changed to study the effect of this
    for (size_t i = 0; i < v_trd_y.size(); ++i) {
        if ((abs(v_trd_y[i] - extrp_y) < trd_y_diff_threshold) && v_trk_ytime_coincidence[i]) return kTRUE;
    }
    return kFALSE;
}

void study_efficiency(TString evt_root_file=RootOutputPath+"trd_singleTrackHits_Run_004349_001.root", TString hist_root_file=""){

    ROOT::EnableImplicitMT();
    
    auto df = ROOT::RDataFrame(tree_name.Data(), evt_root_file.Data());
    auto h_ych = df.Histo1D({"h_ych", "Y Channel", 528, -0.5,527.5}, "ych");

    // TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    // h_ych->Draw("HIST");
    // c1->SaveAs("ych.pdf");
    
    vector<ROOT::RDF::RResultPtr<TH1D>> h_ych_filtered, h_is_match_hit, h_is_match_hit_true, h_ych_filtered_true, h_gem_trk_delta_x;
    vector<ROOT::RDF::RResultPtr<TH1D>> h_gem_trk_delta_x_cut, h_ych_filtered_true_deltaxcut;
    vector<ROOT::RDF::RResultPtr<TH2D>> h_ytime_vs_ych_filtered, h_trkrx1_vs_trkrx2_filtered, h_trkrx1_vs_trdx_filtered;
    vector<ROOT::RDF::RResultPtr<TH2D>> h_trkrx1_vs_trkrx2_filtered_deltaxcut;
    
    vector<TString> h_radiator_names_var = radiator_names;
    for (auto &name : h_radiator_names_var) name.ReplaceAll(" ", "_");
    
    for (size_t i = 0; i < ychan_limits.size(); ++i) {
        auto df_filtered = df.Filter(Form("ych > %d && ych < %d", ychan_limits[i].first, ychan_limits[i].second));

        df_filtered = df_filtered.Define("is_match_hit", GetIsMatchHit, {"extrp_y", "v_trd_y", "trk_ytime_coincidence"});
        h_is_match_hit.push_back(df_filtered.Histo1D({Form("h_is_match_hit_%s", h_radiator_names_var[i].Data()), radiator_names[i], 2, -0.5, 1.5}, "is_match_hit"));
        h_ych_filtered.push_back(df_filtered.Histo1D({Form("h_ych_%s", h_radiator_names_var[i].Data()), radiator_names[i], 528, -0.5, 527.5}, "ych"));

        df_filtered = df_filtered.Filter("is_match_hit == 1");
        h_is_match_hit_true.push_back(df_filtered.Histo1D({Form("h_is_match_hit_true_%s", h_radiator_names_var[i].Data()), radiator_names[i], 2, -0.5, 1.5}, "is_match_hit"));
        h_ych_filtered_true.push_back(df_filtered.Histo1D({Form("h_ych_%s_true", h_radiator_names_var[i].Data()), radiator_names[i], 528, -0.5, 527.5}, "ych"));

        h_gem_trk_delta_x.push_back(df_filtered.Histo1D({Form("h_gem_trk_delta_x_%s", h_radiator_names_var[i].Data()), radiator_names[i], 40, 0.0, 5}, "GEMTrkrsDeltaX"));
        h_trkrx1_vs_trkrx2_filtered.push_back(df_filtered.Histo2D({Form("h_trkrx1_vs_trkrx2_%s", h_radiator_names_var[i].Data()), radiator_names[i], 512, -0.5, 511.5, 512, -0.5, 511.5}, "xchtrkr1", "xchtrkr2"));

        df_filtered = df_filtered.Filter("GEMTrkrsDeltaX < 0.5");
        h_gem_trk_delta_x_cut.push_back(df_filtered.Histo1D({Form("h_gem_trk_delta_x_cut_%s", h_radiator_names_var[i].Data()), radiator_names[i], 40, 0.0, 5}, "GEMTrkrsDeltaX"));
        h_trkrx1_vs_trkrx2_filtered_deltaxcut.push_back(df_filtered.Histo2D({Form("h_trkrx1_vs_trkrx2_%s_deltaxcut", h_radiator_names_var[i].Data()), radiator_names[i], 512, -0.5, 511.5, 512, -0.5, 511.5}, "xchtrkr1", "xchtrkr2"));
        h_ych_filtered_true_deltaxcut.push_back(df_filtered.Histo1D({Form("h_ych_%s_true_deltaxcut", h_radiator_names_var[i].Data()), radiator_names[i], 528, -0.5, 527.5}, "ych"));

        // df_filtered = df_filtered.Filter("GENTrkrsDeltaX < ");
        // h_ytime_vs_ych_filtered.push_back(df_filtered.Histo2D({Form("h_ytime_vs_ych_%s", h_radiator_names_var[i].Data()), radiator_names[i], 251, -0.5, 250, 528, -0.5, 527.5}, "ytime", "ych"));
        // h_trkrx1_vs_trdx_filtered.push_back(df_filtered.Histo2D({Form("h_trkrx1_vs_trdx_%s", h_radiator_names_var[i].Data()), radiator_names[i], 201, -1, 401, 528, -0.5, 527.5}, "xchtrkr1", "ych"));
        // h_trkrx1_vs_trdx_filtered.back()->GetXaxis()->SetTitle("GEM tracker 1 hit x channel");
        // h_trkrx1_vs_trdx_filtered.back()->GetYaxis()->SetTitle("GEM-TRD x channel");

    }

    SetYChAxisRadiator(h_ych_filtered);
    SetYChAxisRadiator(h_ych_filtered_true);
    SetYChAxisRadiator(h_ych_filtered_true_deltaxcut);
    SetDeltaXAxisRadiator(h_gem_trk_delta_x);

    // PrintRadiator1DHistogram(h_ych_filtered, "ych_filtered.pdf");
    PrintRadiator1DHistogram(h_is_match_hit, "is_match_hit.pdf");
    PrintRadiator1DHistogram(h_is_match_hit_true, "is_match_hit_true.pdf");
    PrintRadiator1DHistogram(h_ych_filtered, "ych_filtered.pdf");
    PrintRadiator1DHistogram(h_ych_filtered_true, "ych_filtered_true.pdf");
    PrintEfficiencyHistogram(h_ych_filtered, h_ych_filtered_true, "ych_efficiency.pdf");
    PrintRadiator1DHistogram(h_gem_trk_delta_x, "gem_trk_delta_x.pdf");

    PrintRadiator2DHistogram(h_trkrx1_vs_trkrx2_filtered, "trkrx1_vs_trkrx2.pdf");

    PrintRadiator1DHistogram(h_gem_trk_delta_x_cut, "gem_trk_delta_x_cut.pdf");
    PrintEfficiencyHistogram(h_ych_filtered, h_ych_filtered_true_deltaxcut, "ych_efficiency_deltaxcut.pdf");

}
