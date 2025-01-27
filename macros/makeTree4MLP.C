const TString tree_name = "gem_trackers_hits";
const TString RootOutputPath = "/work/halld2/home/nseptian/halld24_TestBeam/RootOutput/halld24/";

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

Double_t GetTRDHitsExtrpYCh(Double_t extrp_y) {
    // cout << "extrp_y: " << extrp_y << endl;
    return (extrp_y - dy3);
}

Bool_t GetIsMatchHit(Double_t extrp_y, vector<Double_t> v_trd_y, vector<bool> v_trk_ytime_coincidence){
    const Double_t trd_y_diff_threshold = 20; // same as in trd class, can be changed to study the effect of this
    // cout << "Number of TRD hits: " << v_trd_y.size() << endl;
    if (v_trd_y.size() == 0) return kFALSE;
    // need to evaluate if we only want to compare the max of v_trd_y with the extrp_y
    // for (size_t i = 0; i < v_trd_y.size(); ++i) {
    //     if ((abs(v_trd_y[i] - extrp_y) < trd_y_diff_threshold) && v_trk_ytime_coincidence[i]) return kTRUE;
    // }
    // return kFALSE;

    Double_t max_trd_y = *std::max_element(v_trd_y.begin(), v_trd_y.end());

    if (((max_trd_y - extrp_y) < trd_y_diff_threshold)) return kTRUE;
    else return kFALSE;

}

void PlotDiagnostics(ROOT::RDF::RNode df, TString tag){
    auto h_extrp_ych = df.Histo1D({"h_extrp_ych", "extrp_ych", 528, -0.5, 527.5}, "extrp_ych");
    auto h_ych = df.Histo1D({"h_ych", "ych", 528, -0.5, 527.5}, "v_trd_y");
    auto h_is_match_hit = df.Histo1D({"h_is_match_hit", "is_match_hit", 2, -0.5, 1.5}, "is_match_hit");
    auto h_gem_trk_delta_x = df.Histo1D({"h_gem_trk_delta_x", "GEMTrkrsDeltaX", 40, 0.0, 5}, "GEMTrkrsDeltaX");
    auto h_trkrx1_vs_trkrx2 = df.Histo2D({"h_trkrx1_vs_trkrx2", "trkrx1_vs_trkrx2", 512, -0.5, 511.5, 512, -0.5, 511.5}, "xchtrkr1", "xchtrkr2");
    auto h_trkrx1_vs_trdx = df.Histo2D({"h_trkrx1_vs_trdx", "trkrx1_vs_trdx", 512, -0.5, 511.5, 512, -0.5, 511.5}, "xchtrkr1", "extrp_ych");

    TCanvas *c = new TCanvas("c", "c", 1600, 1200);
    c->Divide(2, 3);
    c->cd(1);
    h_extrp_ych->Draw();
    c->cd(2);
    h_ych->Draw();
    c->cd(3);
    h_is_match_hit->Draw();
    c->cd(4);
    h_gem_trk_delta_x->Draw();
    c->cd(5);
    h_trkrx1_vs_trkrx2->Draw("COLZ");
    c->cd(6);
    h_trkrx1_vs_trdx->Draw("COLZ");
    c->SaveAs("diagnostics_"+tag+".pdf");
    delete c;
}

void PlotDiagnostics2(ROOT::RDF::RNode df, TString tag){
    auto h2_zpos_vs_xpos = df.Histo2D({"h2_zpos_vs_xpos", "Z Position (time) vs X Position", 120, 0, 120, 200, 0, 200}, "xpos", "zpos", "dedx");
    // auto h2_zpos_vs_ypos = df.Histo2D({"h2_zpos_vs_ypos", "Z Position vs Y Position", 200, 0, 200, 200, 0, 200}, "ypos", "zpos");
    auto h2_ytime_vs_ych = df.Histo2D({"h2_ytime_vs_ych", "Y Max Time vs Y Max Channel (log scale)", 528, -0.5, 527.5, 200, 0, 200}, "ych", "ytime", "yamp");
    auto h2_zpos_vs_dedx = df.Histo2D({"h2_zpos_vs_dedx", "Z Position vs dEdx", 100, 150, 2000, 200, 0, 200}, "dedx", "zpos");
    auto h2_ytime_vs_yamp = df.Histo2D({"h2_ytime_vs_yamp", "Y Max Time vs Y Max Amplitude", 100, 150, 2000, 200, 0, 200}, "yamp", "ytime");

    TCanvas *c = new TCanvas("c", "c", 1600, 800);
    c->Divide(2, 2);
    c->cd(1);
    h2_zpos_vs_xpos->Draw("COLZ");
    c->cd(2)->SetLogz();
    h2_ytime_vs_ych->Draw("COLZ");
        
    // c->cd(2);
    // h2_zpos_vs_ypos->Draw("COLZ");
    c->cd(3)->SetLogz();
    h2_zpos_vs_dedx->Draw("COLZ");

    c->cd(4);
    h2_ytime_vs_yamp->Draw("COLZ");

    c->SaveAs("diagnostics2_"+tag+".pdf");

    delete c;
}

void makeTree4MLP(TString evt_root_file_wRadiator=RootOutputPath+"trd_singleTrackHits_Run_004349.root", TString evt_root_file_noRadiator=RootOutputPath+"trd_singleTrackHits_Run_004349.root"){

    ROOT::EnableImplicitMT();
    
    auto df_wRadiators = ROOT::RDataFrame(tree_name.Data(), evt_root_file_wRadiator.Data());

    auto df_wRadiators1 = df_wRadiators.Define("extrp_ych", GetTRDHitsExtrpYCh, {"extrp_y"});
    auto df_wRadiators2 = df_wRadiators1.Define("is_match_hit", GetIsMatchHit, {"extrp_y", "v_trd_y", "trk_ytime_coincidence"});
    PlotDiagnostics(df_wRadiators2, "wRadiators");
    PlotDiagnostics2(df_wRadiators2, "wRadiators");
    auto df_wRadiators3 = df_wRadiators2.Filter("is_match_hit == 1");
    PlotDiagnostics(df_wRadiators3, "wRadiatorsMatchHit");
    PlotDiagnostics2(df_wRadiators3, "wRadiatorsMatchHit");
    auto df_wRadiators4 = df_wRadiators3.Filter("GEMTrkrsDeltaX < 5");
    PlotDiagnostics(df_wRadiators4, "wRadiatorsMatchHitDeltaXCut");
    PlotDiagnostics2(df_wRadiators4, "wRadiatorsMatchHitDeltaXCut");

    auto df_wRadiators5 = df_wRadiators4.Filter("ych > 160");
    PlotDiagnostics(df_wRadiators5, "wRadiatorsMatchHitDeltaXCutYChCut");
    PlotDiagnostics2(df_wRadiators5, "wRadiatorsMatchHitDeltaXCutYChCut");

    auto df_wRadiators6 = df_wRadiators5.Filter("dedx > 600");
    PlotDiagnostics(df_wRadiators6, "wRadiatorsMatchHitDeltaXCutYChCutYampCut");
    PlotDiagnostics2(df_wRadiators6, "wRadiatorsMatchHitDeltaXCutYChCutYampCut");


}