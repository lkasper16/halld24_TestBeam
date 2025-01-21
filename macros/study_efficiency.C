
const TString tree_name = "gem_trackers_hits";
const TString RootOutputPath = "/work/halld2/home/nseptian/halld24_TestBeam/RootOutput/halld24/";

// no radiator min,max = 0, 150
// hermes radiator min,max = 150, 267
// VU foil radiator min,max = 267, 390
// Fleece radiator min,max = 390, 511

// Detectors coordinate from Survey

Double_t Y0=80.59990;
Double_t v0=104.69986;
Double_t x0=386.76755;

Double_t Y1=(81.18291+80.72291)/2.;
Double_t v1=(104.76201+104.76228)/2.;
Double_t x1=(388.82626+388.83477)/2.;

Double_t y2=(81.25059+80.79095)/2.;
Double_t v2=(104.76775+104.77306)/2.;
Double_t x2=(389.41947+389.43687)/2.;

Double_t y3=(81.32000+80.77612+81.32446+80.78202)/4.;
Double_t v3=(105.11588+105.11209+104.38180+104.37675)/4.;
Double_t x3=(391.54396+391.53969+391.54030+391.53604)/4.;

Double_t Bl = 1.6617 * 36*0.0254; // 1.6T x 36" field integral

Double_t dx1=(x1-x0)*1000.;
Double_t dY1=(Y1-Y0)*1000.-256.*0.8-31.75; 
Double_t dx2=(x2-x0)*1000.;
Double_t dy2=(y2-Y0)*1000.-256.*0.8-31.75-2.4;
Double_t dx3=(x3-x0)*1000.+31.35;
Double_t dy3=(y3-Y0)*1000.-264.;  

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

void PrintDividedHistogram(vector<ROOT::RDF::RResultPtr<TH1D>> v_h1, vector<ROOT::RDF::RResultPtr<TH1D>> v_h2, TString pdfName){
    TCanvas *c = new TCanvas("c", "c", 1600, 1200);
    TH1D *h1[4];
    c->Divide(2, 2);
    for (size_t i = 0; i < v_h1.size(); ++i) {
        c->cd(i + 1);
        h1[i] = (TH1D*)v_h2[i]->Clone(Form("%s_clone", v_h1[i]->GetName()));
        h1[i]->Divide(v_h1[i].GetPtr());
        h1[i]->Draw();
    }
    c->SaveAs(pdfName);
    delete c;
    for (size_t i = 0; i < v_h1.size(); ++i) delete h1[i];
};

Double_t GetTRDHitsExtrpYCh(Double_t extrp_y) {
    return (extrp_y - dy3);
}

Bool_t GetIsMatchHit(Double_t extrp_y, vector<Double_t> v_trd_y, vector<bool> v_trk_ytime_coincidence){
    const Double_t trd_y_diff_threshold = 20; // same as in trd class, can be changed to study the effect of this
    cout << "Number of TRD hits: " << v_trd_y.size() << endl;
    for (size_t i = 0; i < v_trd_y.size(); ++i) {
        if ((abs(v_trd_y[i] - extrp_y) < trd_y_diff_threshold) && v_trk_ytime_coincidence[i]) return kTRUE;
    }
    return kFALSE;
}

void study_efficiency(TString evt_root_file=RootOutputPath+"trd_singleTrackHits_Run_004349.root", TString hist_root_file=""){

    ROOT::EnableImplicitMT();
    
    auto df = ROOT::RDataFrame(tree_name.Data(), evt_root_file.Data());
    auto h_ych = df.Histo1D({"h_ych", "Y Channel", 528, -0.5,527.5}, "ych");

    // TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    // h_ych->Draw("HIST");
    // c1->SaveAs("ych.pdf");
    
    vector<ROOT::RDF::RResultPtr<TH1D>> h_ych_filtered, h_extrp_ych_filtered;
    vector<ROOT::RDF::RResultPtr<TH1D>> h_is_match_hit, h_is_match_hit_MatchHit, h_ych_filtered_MatchHit, h_extrp_ych_filtered_MatchHit, h_gem_trk_delta_x;
    vector<ROOT::RDF::RResultPtr<TH1D>> h_gem_trk_delta_x_cut, h_ych_filtered_MatchHit_deltaxcut, h_extrp_ych_filtered_MatchHit_deltaxcut;
    vector<ROOT::RDF::RResultPtr<TH2D>> h_ytime_vs_ych_filtered, h_trkrx1_vs_trkrx2_filtered, h_trkrx1_vs_trdx_filtered;
    vector<ROOT::RDF::RResultPtr<TH2D>> h_trkrx1_vs_trkrx2_filtered_deltaxcut;
    
    vector<TString> h_radiator_names_var = radiator_names;
    for (auto &name : h_radiator_names_var) name.ReplaceAll(" ", "_");
    
    for (size_t i = 0; i < ychan_limits.size(); ++i) {
        auto df_filtered = df.Filter(Form("ych > %d && ych < %d", ychan_limits[i].first, ychan_limits[i].second));

        df_filtered = df_filtered.Define("is_match_hit", GetIsMatchHit, {"extrp_y", "v_trd_y", "trk_ytime_coincidence"});
        df_filtered = df_filtered.Define("extrp_ych", GetTRDHitsExtrpYCh, {"extrp_y"});

        h_is_match_hit.push_back(df_filtered.Histo1D({Form("h_is_match_hit_%s", h_radiator_names_var[i].Data()), radiator_names[i], 2, -0.5, 1.5}, "is_match_hit"));
        h_ych_filtered.push_back(df_filtered.Histo1D({Form("h_ych_%s", h_radiator_names_var[i].Data()), radiator_names[i], 528, -0.5, 527.5}, "ych"));
        h_extrp_ych_filtered.push_back(df_filtered.Histo1D({Form("h_extrp_ych_%s", h_radiator_names_var[i].Data()), radiator_names[i], 528, -0.5, 527.5}, "extrp_ych"));

        df_filtered = df_filtered.Filter("is_match_hit == 1");
        h_is_match_hit_MatchHit.push_back(df_filtered.Histo1D({Form("h_is_match_hit_MatchHit_%s", h_radiator_names_var[i].Data()), radiator_names[i], 2, -0.5, 1.5}, "is_match_hit"));
        h_ych_filtered_MatchHit.push_back(df_filtered.Histo1D({Form("h_ych_%s_MatchHit", h_radiator_names_var[i].Data()), radiator_names[i], 528, -0.5, 527.5}, "ych"));
        h_extrp_ych_filtered_MatchHit.push_back(df_filtered.Histo1D({Form("h_extrp_ych_%s_MatchHit", h_radiator_names_var[i].Data()), radiator_names[i], 528, -0.5, 527.5}, "extrp_ych"));

        h_gem_trk_delta_x.push_back(df_filtered.Histo1D({Form("h_gem_trk_delta_x_%s", h_radiator_names_var[i].Data()), radiator_names[i], 40, 0.0, 5}, "GEMTrkrsDeltaX"));
        h_trkrx1_vs_trkrx2_filtered.push_back(df_filtered.Histo2D({Form("h_trkrx1_vs_trkrx2_%s", h_radiator_names_var[i].Data()), radiator_names[i], 512, -0.5, 511.5, 512, -0.5, 511.5}, "xchtrkr1", "xchtrkr2"));

        df_filtered = df_filtered.Filter("GEMTrkrsDeltaX < 5");
        h_gem_trk_delta_x_cut.push_back(df_filtered.Histo1D({Form("h_gem_trk_delta_x_cut_%s", h_radiator_names_var[i].Data()), radiator_names[i], 40, 0.0, 5}, "GEMTrkrsDeltaX"));
        h_trkrx1_vs_trkrx2_filtered_deltaxcut.push_back(df_filtered.Histo2D({Form("h_trkrx1_vs_trkrx2_%s_deltaxcut", h_radiator_names_var[i].Data()), radiator_names[i], 512, -0.5, 511.5, 512, -0.5, 511.5}, "xchtrkr1", "xchtrkr2"));
        h_ych_filtered_MatchHit_deltaxcut.push_back(df_filtered.Histo1D({Form("h_ych_%s_MatchHit_deltaxcut", h_radiator_names_var[i].Data()), radiator_names[i], 528, -0.5, 527.5}, "ych"));
        h_extrp_ych_filtered_MatchHit_deltaxcut.push_back(df_filtered.Histo1D({Form("h_extrp_ych_%s_MatchHit_deltaxcut", h_radiator_names_var[i].Data()), radiator_names[i], 528, -0.5, 527.5}, "extrp_ych"));

        // h_ytime_vs_ych_filtered.push_back(df_filtered.Histo2D({Form("h_ytime_vs_ych_%s", h_radiator_names_var[i].Data()), radiator_names[i], 251, -0.5, 250, 528, -0.5, 527.5}, "ytime", "ych"));
        // h_trkrx1_vs_trdx_filtered.push_back(df_filtered.Histo2D({Form("h_trkrx1_vs_trdx_%s", h_radiator_names_var[i].Data()), radiator_names[i], 201, -1, 401, 528, -0.5, 527.5}, "xchtrkr1", "ych"));
        // h_trkrx1_vs_trdx_filtered.back()->GetXaxis()->SetTitle("GEM tracker 1 hit x channel");
        // h_trkrx1_vs_trdx_filtered.back()->GetYaxis()->SetTitle("GEM-TRD x channel");

    }

    SetYChAxisRadiator(h_ych_filtered);
    SetYChAxisRadiator(h_ych_filtered_MatchHit);
    SetYChAxisRadiator(h_ych_filtered_MatchHit_deltaxcut);
    SetDeltaXAxisRadiator(h_gem_trk_delta_x);

    // PrintRadiator1DHistogram(h_ych_filtered, "ych_filtered.pdf");
    PrintRadiator1DHistogram(h_is_match_hit, "is_match_hit.pdf");
    PrintRadiator1DHistogram(h_is_match_hit_MatchHit, "is_match_hit_MatchHit.pdf");
    PrintRadiator1DHistogram(h_ych_filtered, "ych_filtered.pdf");
    PrintRadiator1DHistogram(h_extrp_ych_filtered, "extrp_ych_filtered.pdf");

    PrintRadiator1DHistogram(h_ych_filtered_MatchHit, "ych_filtered_MatchHit.pdf");
    PrintDividedHistogram(h_ych_filtered, h_ych_filtered_MatchHit, "ych_matchhit_ratio.pdf");
    PrintRadiator1DHistogram(h_gem_trk_delta_x, "gem_trk_delta_x.pdf");
    PrintRadiator2DHistogram(h_trkrx1_vs_trkrx2_filtered, "trkrx1_vs_trkrx2.pdf");
    PrintRadiator1DHistogram(h_extrp_ych_filtered_MatchHit, "extrp_ych_filtered_MatchHit.pdf");

    PrintRadiator1DHistogram(h_gem_trk_delta_x_cut, "gem_trk_delta_x_cut.pdf");
    PrintRadiator2DHistogram(h_trkrx1_vs_trkrx2_filtered_deltaxcut, "trkrx1_vs_trkrx2_deltaxcut.pdf");
    PrintDividedHistogram(h_ych_filtered, h_ych_filtered_MatchHit_deltaxcut, "ych_deltaxcut_ratio.pdf");
    PrintRadiator1DHistogram(h_extrp_ych_filtered_MatchHit_deltaxcut, "extrp_ych_filtered_MatchHit_deltaxcut.pdf");

    PrintDividedHistogram(h_extrp_ych_filtered, h_extrp_ych_filtered_MatchHit, "extrp_ych_matchhit_efficiency.pdf");
    PrintDividedHistogram(h_extrp_ych_filtered, h_extrp_ych_filtered_MatchHit_deltaxcut, "extrp_ych_matchhit_efficiency_deltaxcut.pdf");

}
