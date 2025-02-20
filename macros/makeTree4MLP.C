const TString tree_name = "gem_trackers_hits";
const TString RootOutputPath = "/work/halld2/home/nseptian/halld24_TestBeam/RootOutput/";

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

Double_t GETYChExtrpYChDiff(Double_t extrp_ych, Int_t ych) {
    return (extrp_ych - ych);
}

Bool_t GetIsMatchHit(Double_t extrp_y, vector<Double_t> v_trd_y, vector<bool> v_trk_ytime_coincidence){
    const Double_t trd_y_diff_threshold = 20; // same as in trd class, can be changed to study the effect of this
    // cout << "Number of TRD hits: " << v_trd_y.size() << endl;
    if (v_trd_y.size() == 0) return kFALSE;
    // need to evaluate if we only want to compare the max of v_trd_y with the extrp_y
    for (size_t i = 0; i < v_trd_y.size(); ++i) {
        if ((abs(v_trd_y[i] - extrp_y) < trd_y_diff_threshold) && v_trk_ytime_coincidence[i]) return kTRUE;
    }
    return kFALSE;
}

void PlotDiagnostics(ROOT::RDF::RNode df, TString tag){
    auto h_ych = df.Histo1D({"h_ych", "GEM-TRD Y Channel Max Amp", 528, -0.5, 527.5}, "ych");
    auto h_extrp_ych = df.Histo1D({"h_extrp_ych", "GEM-TRD Y Channel Extrapolated from Trackers", 528, -0.5, 527.5}, "extrp_ych");
    auto h_ych_diff = df.Histo1D({"h_ych_diff", "Y Channel Difference", 201, -100.5, 100.5}, "ydiff");
    auto h_trkrx1_vs_trkrx2 = df.Histo2D({"h_trkrx1_vs_trkrx2", "GEM Trackers Correlation", 512, -0.5, 511.5, 512, -0.5, 511.5}, "xchtrkr1", "xchtrkr2");
    auto h_gem_trk_delta_x = df.Histo1D({"h_gem_trk_delta_x", "#Delta X_{trks}", 40, 0.0, 8}, "GEMTrkrsDeltaX");
    auto h_trkrx1_vs_trdx = df.Histo2D({"h_trkrx1_vs_trdx", "GEM Tracker 1 vs GEM-TRD Y Ch Max", 512, -0.5, 511.5, 512, -0.5, 511.5}, "xchtrkr1", "ych");

    TCanvas *c = new TCanvas("c", "c", 1600, 800);
    c->Divide(3, 2);
    c->cd(1);
    h_ych->Draw();
    c->cd(2);
    h_extrp_ych->Draw();
    c->cd(3);
    h_ych_diff->Draw();
    c->cd(4);
    h_trkrx1_vs_trkrx2->Draw("COLZ");
    c->cd(5);
    h_gem_trk_delta_x->Draw();
    c->cd(6);
    h_trkrx1_vs_trdx->Draw("COLZ");
    c->SaveAs("diagnostics_"+tag+".pdf");
    delete c;
}

void PlotDiagnostics2(ROOT::RDF::RNode df, TString tag){
    auto h2_zpos_vs_xpos = df.Histo2D({"h2_zpos_vs_xpos", "Z Position (time) vs X Position", 120, 0, 120, 200, 0, 200}, "xpos", "zpos", "dedx");
    // auto h2_zpos_vs_ypos = df.Histo2D({"h2_zpos_vs_ypos", "Z Position vs Y Position", 200, 0, 200, 200, 0, 200}, "ypos", "zpos");
    auto h2_ytime_vs_ych = df.Histo2D({"h2_ytime_vs_ych", "Y Max Time vs Y Max Channel", 528, -0.5, 527.5, 200, 0, 200}, "ych", "ytime", "yamp");
    auto h2_zpos_vs_dedx = df.Histo2D({"h2_zpos_vs_dedx", "Z Position vs dEdx", 100, 150, 5000, 200, 0, 200}, "dedx", "zpos");
    auto h2_ytime_vs_yamp = df.Histo2D({"h2_ytime_vs_yamp", "Y Max Time vs Y Max Amplitude", 100, 150, 5000, 200, 0, 200}, "yamp", "ytime");

    TCanvas *c = new TCanvas("c", "c", 1600, 800);
    c->Divide(2, 2);
    c->cd(1);
    h2_zpos_vs_xpos->Draw("COLZ");
    c->cd(2);
    h2_ytime_vs_ych->Draw("COLZ");
        
    // c->cd(2);
    // h2_zpos_vs_ypos->Draw("COLZ");
    c->cd(3);
    h2_zpos_vs_dedx->Draw("COLZ");

    c->cd(4);
    h2_ytime_vs_yamp->Draw("COLZ");

    c->SaveAs("diagnostics2_"+tag+".pdf");

    delete c;
}

void PlotDiagnostics3(ROOT::RDF::RNode df, TString tag){
    auto h2_zpos_vs_xpos = df.Histo2D({"h2_zpos_vs_xpos", "Z Position (time) vs X Position", 120, 0, 120, 200, 0, 200}, "xpos", "zpos", "dedx");
    auto h2_ytime_vs_ych = df.Histo2D({"h2_ytime_vs_ych", "Y Max Time vs Y Max Channel", 528, -0.5, 527.5, 200, 0, 200}, "ych", "ytime", "yamp");

    auto h1_zpos = df.Histo1D({"h1_zpos", "Z Position (time)", 100, 0, 200}, "zpos", "dedx");
    auto h1_ytime = df.Histo1D({"h1_ytime", "Y Max Time (Z position)", 100, 0, 200}, "ytime", "yamp");

    TCanvas *c = new TCanvas("c", "c", 1600, 800);
    c->Divide(2, 2);
    c->cd(1);
    h2_zpos_vs_xpos->Draw("COLZ");
    c->cd(2);
    h2_ytime_vs_ych->Draw("COLZ");
    c->cd(3);
    h1_zpos->Draw();
    c->cd(4);
    h1_ytime->Draw();
    c->SaveAs("diagnostics3_"+tag+".pdf");
    delete c;
}

void PlotDiagnostics4(ROOT::RDF::RNode df, TString tag){
    auto h_xytimediff = df.Histo1D({"h_xytimediff", "X - Y Time Difference", 100, -50, 50}, "xytimediff");
    // auto h_xpos_filtered = df.Histo1D({"h_xpos_filtered", "X Position Filtered by X,Y Time Difference", 100, 0, 200}, "xpos_filtered_time");
    // auto h_ypos_filtered = df.Histo1D({"h_ypos_filtered", "Y Position Filtered by X,Y Time Difference", 100, 0, 200}, "ypos_filtered_time");

    TCanvas *c = new TCanvas("c", "c", 800, 800);
    // c->Divide(2,1);
    // c->cd(1);
    h_xytimediff->Draw();
    // c->cd(2);
    // h_xpos_filtered->Draw();
    // c->cd(3);
    // h_ypos_filtered->Draw();
    c->SaveAs("diagnostics4_"+tag+".pdf");
    delete c;
}

void PlotDiagnosticTrackers(ROOT::RDF::RNode df, TString tag){
    auto h_ych = df.Histo1D({"h_ych", "GEM-TRD Y Channel Max Amp", 528, -0.5, 527.5}, "ych");
    auto h_extrp_ych = df.Histo1D({"h_extrp_ych", "GEM-TRD Y Channel Extrapolated from Trackers", 528, -0.5, 527.5}, "extrp_ych");
    auto h_trkrx1_vs_trkrx2 = df.Histo2D({"h_trkrx1_vs_trkrx2", "GEM Trackers Correlation", 512, -0.5, 511.5, 512, -0.5, 511.5}, "xchtrkr1", "xchtrkr2");
    auto h_trkrx1_vs_ych = df.Histo2D({"h_trkrx1_vs_ych", "GEM Tracker 1 vs GEM-TRD Correlation", 512, -0.5, 511.5, 528, -0.5, 527.5}, "xchtrkr1", "ych");

    TCanvas *c = new TCanvas("c", "c", 1600, 800);
    c->Divide(2, 2);
    c->cd(1);
    h_ych->Draw();
    c->cd(2);
    h_extrp_ych->Draw();
    c->cd(3);
    h_trkrx1_vs_trkrx2->Draw("COLZ");
    c->cd(4);
    h_trkrx1_vs_ych->Draw("COLZ");
    c->SaveAs("diagnostics_trackers_"+tag+".pdf");
    delete c;
}

void PlotDiagnosticTrackers2(ROOT::RDF::RNode df, TString tag){
    auto h_ych = df.Histo1D({"h_ych", "GEM-TRD Y Channel Max Amp", 528, -0.5, 527.5}, "ych");
    auto h_extrp_ych = df.Histo1D({"h_extrp_ych", "GEM-TRD Y Channel Extrapolated from Trackers", 528, -0.5, 527.5}, "extrp_ych");
    auto h_ydiff = df.Histo1D({"h_ydiff", "Y Channel Extrapolated - Y Channel", 528, -0.5, 527.5}, "ydiff");

    TCanvas *c = new TCanvas("c", "c", 1600, 800);
    c->Divide(3,1);
    c->cd(1);
    h_ych->Draw();
    c->cd(2);
    h_extrp_ych->Draw();
    c->cd(3);
    h_ydiff->Draw();
    c->SaveAs("diagnostics_trackers2_"+tag+".pdf");
    delete c;
}

void PlotDiagnosticsXYTCoin(ROOT::RDF::RNode df, TString tag){
    auto h_xpos = df.Histo2D({"h_xpos", "no time cut", 121, -0.5, 120.5, 200, 0, 200}, "xpos", "zpos", "dedx");
    auto h_ypos = df.Histo2D({"h_ypos", "no time cut", 528, -0.5, 527.5, 200, 0, 200}, "ypos", "ypos_time" ,"ypos_amp");
    auto h_xpos_tcoin = df.Histo2D({"h_xpos_tcoin", "|t_{X}-t_{Y}| < 2", 121, -0.5, 120.5, 200, 0, 200}, "v_tcoin_xpos", "v_tcoin_xpos_time", "v_tcoin_xpos_amp");
    auto h_ypos_tcoin = df.Histo2D({"h_ypos_tcoin", "|t_{X}-t_{Y}| < 2", 528, -0.5, 527.5, 200, 0, 200}, "v_tcoin_ypos", "v_tcoin_ypos_time", "v_tcoin_ypos_amp");

    TCanvas *c = new TCanvas("c", "c", 1600, 800);
    c->Divide(2, 2);
    c->cd(1);
    h_xpos->GetXaxis()->SetTitle("X channel");
    h_xpos->GetYaxis()->SetTitle("time (sample)");
    h_xpos->Draw("COLZ");
    c->cd(2);
    h_ypos->GetXaxis()->SetTitle("Y channel");
    h_ypos->GetYaxis()->SetTitle("time (sample)");
    h_ypos->Draw("COLZ");
    c->cd(3);
    h_xpos_tcoin->GetXaxis()->SetTitle("X channel");
    h_xpos_tcoin->GetYaxis()->SetTitle("time (sample)");
    h_xpos_tcoin->Draw("COLZ");
    c->cd(4);
    h_ypos_tcoin->GetXaxis()->SetTitle("Y channel");
    h_ypos_tcoin->GetYaxis()->SetTitle("time (sample)");
    h_ypos_tcoin->Draw("COLZ");
    c->SaveAs("diagnostics_xy_tcoin_"+tag+".pdf");
    delete c;
}

void PlotDiagnosticsXYTCoin2(ROOT::RDF::RNode df,TString tag){
    auto h_xpos_amp = df.Histo2D({"h_xpos_amp", "no time cut", 3000, 0, 3000, 121, -0.5, 120.5}, "dedx", "xpos");
    auto h_ypos_amp = df.Histo2D({"h_ypos_amp", "no time cut", 3000, 0, 3000, 528, -0.5, 527.5}, "ypos_amp", "ypos");
    auto h_xpos_tcoin_amp = df.Histo2D({"h_xpos_tcoin_amp", "|t_{X}-t_{Y}| < 2", 3000, 0, 3000, 121, -0.5, 120.5}, "v_tcoin_xpos_amp", "v_tcoin_xpos");
    auto h_ypos_tcoin_amp = df.Histo2D({"h_ypos_tcoin_amp", "|t_{X}-t_{Y}| < 2", 3000, 0, 3000, 528, -0.5, 527.5}, "v_tcoin_ypos_amp", "v_tcoin_ypos");

    TCanvas *c = new TCanvas("c", "c", 1600, 800);
    c->Divide(2, 2);
    c->cd(1);
    h_xpos_amp->GetXaxis()->SetTitle("X channel amplitude");
    h_xpos_amp->GetYaxis()->SetTitle("X channel");
    h_xpos_amp->Draw("COLZ");
    c->cd(2);
    h_ypos_amp->GetXaxis()->SetTitle("Y channel amplitude");
    h_ypos_amp->GetYaxis()->SetTitle("Y channel");
    h_ypos_amp->Draw("COLZ");
    c->cd(3);
    h_xpos_tcoin_amp->GetXaxis()->SetTitle("X channel amplitude");
    h_xpos_tcoin_amp->GetYaxis()->SetTitle("X channel");
    h_xpos_tcoin_amp->Draw("COLZ");
    c->cd(4);
    h_ypos_tcoin_amp->GetXaxis()->SetTitle("Y channel amplitude");
    h_ypos_tcoin_amp->GetYaxis()->SetTitle("Y channel");
    h_ypos_tcoin_amp->Draw("COLZ");
    c->SaveAs("diagnostics_xy_tcoin_amp_"+tag+".pdf");
    delete c;


}

void PlotDiagnosticsXYTCoin3(ROOT::RDF::RNode df,TString tag){
    auto h_xpos_tcoin = df.Histo2D({"h_xpos_tcoin", "|t_{X}-t_{Y}| < 2", 121, -0.5, 120.5, 200, 0, 200}, "v_tcoin_xpos", "v_tcoin_xpos_time", "v_tcoin_xpos_amp");
    auto h_ypos_tcoin = df.Histo2D({"h_ypos_tcoin", "|t_{X}-t_{Y}| < 2", 528, -0.5, 527.5, 200, 0, 200}, "v_tcoin_ypos", "v_tcoin_ypos_time", "v_tcoin_ypos_amp");
    auto h_xpos_tcoin_amp = df.Histo2D({"h_xpos_tcoin_amp", "|t_{X}-t_{Y}| < 2", 3000, 0, 3000, 121, -0.5, 120.5}, "v_tcoin_xpos_amp", "v_tcoin_xpos");
    auto h_ypos_tcoin_amp = df.Histo2D({"h_ypos_tcoin_amp", "|t_{X}-t_{Y}| < 2", 3000, 0, 3000, 528, -0.5, 527.5}, "v_tcoin_ypos_amp", "v_tcoin_ypos");

    TCanvas *c = new TCanvas("c", "c", 1600, 800);
    c->Divide(2, 2);
    c->cd(1);
    h_xpos_tcoin->GetXaxis()->SetTitle("X channel");
    h_xpos_tcoin->GetYaxis()->SetTitle("time (sample)");
    h_xpos_tcoin->Draw("COLZ");
    c->cd(2);
    h_ypos_tcoin->GetXaxis()->SetTitle("Y channel");
    h_ypos_tcoin->GetYaxis()->SetTitle("time (sample)");
    h_ypos_tcoin->Draw("COLZ");
    c->cd(3);
    h_xpos_tcoin_amp->GetXaxis()->SetTitle("X channel amplitude");
    h_xpos_tcoin_amp->GetYaxis()->SetTitle("X channel");
    h_xpos_tcoin_amp->Draw("COLZ");
    c->cd(4);
    h_ypos_tcoin_amp->GetXaxis()->SetTitle("Y channel amplitude");
    h_ypos_tcoin_amp->GetYaxis()->SetTitle("Y channel");
    h_ypos_tcoin_amp->Draw("COLZ");
    c->SaveAs("diagnostics3_xy_tcoin_amp_"+tag+".pdf");
    delete c;
}

void PlotDiagnosticsXYTCoin4(ROOT::RDF::RNode df,TString tag){
    auto h_xpos_pulsecut = df.Histo2D({"h_xpos_pulsecut", "", 121, -0.5, 120.5, 200, 0, 200}, "v_pulsecut_xpos", "v_pulsecut_xpos_time", "v_pulsecut_xpos_amp");
    auto h_ypos_pulsecut = df.Histo2D({"h_ypos_pulsecut", "", 528, -0.5, 527.5, 200, 0, 200}, "v_pulsecut_ypos", "v_pulsecut_ypos_time", "v_pulsecut_ypos_amp");
    auto h_xpos_pulsecut_time = h_xpos_pulsecut->ProjectionY();
    auto h_ypos_pulsecut_time = h_ypos_pulsecut->ProjectionY();

    TCanvas *c = new TCanvas("c", "c", 1600, 800);
    c->Divide(2, 2);
    c->cd(1);
    h_xpos_pulsecut->GetXaxis()->SetTitle("X channel");
    h_xpos_pulsecut->GetYaxis()->SetTitle("time (sample)");
    h_xpos_pulsecut->Draw("COLZ");
    c->cd(2);
    h_ypos_pulsecut->GetXaxis()->SetTitle("Y channel");
    h_ypos_pulsecut->GetYaxis()->SetTitle("time (sample)");
    h_ypos_pulsecut->Draw("COLZ");
    c->cd(3);
    h_xpos_pulsecut_time->GetXaxis()->SetTitle("time (sample)");
    h_xpos_pulsecut_time->GetYaxis()->SetTitle("Counts");
    h_xpos_pulsecut_time->Draw();
    c->cd(4);
    h_ypos_pulsecut_time->GetXaxis()->SetTitle("time (sample)");
    h_ypos_pulsecut_time->GetYaxis()->SetTitle("Counts");
    h_ypos_pulsecut_time->Draw();
    c->SaveAs("diagnostics4_xy_pulsecut_"+tag+".pdf");
    delete c;
}
    


void PlotDiagnosticsXYMaxAmp(ROOT::RDF::RNode df, TString tag){
    auto h_xch_vs_time = df.Histo2D({"h_xch_vs_time", "", 121, -0.5, 120.5, 200, 0, 200}, "xch", "xtime", "xamp");
    auto h_ych_vs_time = df.Histo2D({"h_ych_vs_time", "", 528, -0.5, 527.5, 200, 0, 200}, "ych", "ytime", "yamp");
    auto h_xytime_diff = df.Histo1D({"h_xytime_diff", "", 100, -50, 50}, "xymaxamptimediff");

    TCanvas *c = new TCanvas("c", "c", 1200, 600);
    c->Divide(3,1);
    c->cd(1);
    h_xch_vs_time->GetXaxis()->SetTitle("X channel");
    h_xch_vs_time->GetYaxis()->SetTitle("time (sample)");
    h_xch_vs_time->Draw("COLZ");
    c->cd(2);
    h_ych_vs_time->GetXaxis()->SetTitle("Y channel");
    h_ych_vs_time->GetYaxis()->SetTitle("time (sample)");
    h_ych_vs_time->Draw("COLZ");
    c->cd(3);
    h_xytime_diff->GetXaxis()->SetTitle("t_{X} - t_{Y} (sample)");
    h_xytime_diff->GetYaxis()->SetTitle("Counts");
    h_xytime_diff->Draw();
    c->SaveAs("diagnostics_xy_maxamp_"+tag+".pdf");
    delete c;
}

std::vector<float> GetXYTimeDiff(const std::vector<float>& xtime, const std::vector<float>& ytime) {
    std::vector<float> result;
    for (size_t i = 0; i < xtime.size(); ++i) {
        for (size_t j = 0; j < ytime.size(); ++j) {
            // cout << "xtime: " << xtime[i] << " ytime: " << ytime[j] << endl;
            result.push_back(xtime[i] - ytime[i]);
        }
    }
    return result;
}

int GetXYMaxAmpTimeDiff(int xtime, int ytime) {
    return xtime - ytime;
}

std::vector<int> GetXFilteredTime(const std::vector<int>& xpos, const std::vector<float>& xtime, const std::vector<float>& ytime) {
    std::vector<int> result;
    const float cutoff = 10.0;
    for (size_t i = 0; i < xpos.size(); ++i) {
        for (size_t j = 0; j < ytime.size(); ++j) {
            if (TMath::Abs(xtime[i] - ytime[j]) < cutoff) {
                // cout << "xpos: " << xpos[i] << endl;
                result.push_back(xpos[i]);
                break;
            }
        }
    }
    return result;
}

std::vector<int> GetYFilteredTime(const std::vector<int>& ypos, const std::vector<float>& xtime, const std::vector<float>& ytime) {
    std::vector<int> result;
    const float cutoff = 10.0;
    for (size_t i = 0; i < ypos.size(); ++i) {
        for (size_t j = 0; j < xtime.size(); ++j) {
            if (TMath::Abs(xtime[j] - ytime[i]) < cutoff) {
                cout << "ypos: " << ypos[i] << endl;
                result.push_back(ypos[i]);
                break;
            }
        }
    }
    return result;
}


std::vector<int> GetXPosXChDiff(const std::vector<int>& xpos, int xch) {
    // xpos is all the pulse above the threshold, xch is the channel with the max amplitude

    std::vector<int> result;
    for (size_t i = 0; i < xpos.size(); ++i) {
        result.push_back(xpos[i] - xch);
    }
    return result;
}

std::vector<int> GetYPosYChDiff(const std::vector<int>& ypos, int ych) {
    // ypos is all the pulse above the threshold, ych is the channel with the max amplitude

    std::vector<int> result;
    for (size_t i = 0; i < ypos.size(); ++i) {
        result.push_back(ypos[i] - ych);
    }
    return result;
}

void PlotDiagnosticsPulseMaxDiff(ROOT::RDF::RNode df, TString tag){
    auto h_xpos_xch_diff = df.Histo1D({"h_xpos_xch_diff", "", 100, -50, 50}, "xpos_xch_diff");
    auto h_ypos_ych_diff = df.Histo1D({"h_ypos_ych_diff", "", 100, -50, 50}, "ypos_ych_diff");

    TCanvas *c = new TCanvas("c", "c", 1600, 800);
    c->Divide(2, 1);
    c->cd(1);
    h_xpos_xch_diff->GetXaxis()->SetTitle("X_{ch}^{pulse} - X_{ch}^{max}");
    h_xpos_xch_diff->GetYaxis()->SetTitle("Counts");
    h_xpos_xch_diff->Draw();
    c->cd(2);
    h_ypos_ych_diff->GetXaxis()->SetTitle("Y_{ch}^{pulse} - Y_{ch}^{max}");
    h_ypos_ych_diff->GetYaxis()->SetTitle("Counts");
    h_ypos_ych_diff->Draw();
    c->SaveAs("diagnostics_pulse_max_diff_"+tag+".pdf");
    delete c;
}

bool FilterXPulseMaxDiff(const std::vector<int>& xpos_xch_diff, const std::vector<int>& ypos_ych_diff) {
    for (size_t i = 0; i < xpos_xch_diff.size(); ++i) {
        if (TMath::Abs(xpos_xch_diff[i]) > 2) return false;
    }
    return true;
}

bool FilterYPulseMaxDiff(const std::vector<int>& ypos_ych_diff, const std::vector<int>& xpos_xch_diff) {
    for (size_t i = 0; i < ypos_ych_diff.size(); ++i) {
        if (TMath::Abs(ypos_ych_diff[i]) > 2) return false;
    }
    return true;
}

bool FilterYPos(const std::vector<int>& ypos) {
    for (size_t i = 0; i < ypos.size(); ++i) {
        if (ypos[i] < 150 || ypos[i] > 260) return false;
    }
    return true;
}

void makeTree4MLP(TString evt_root_file_wRadiator=RootOutputPath+"out_hits_Run_004349_thr150.root", TString evt_root_file_noRadiator=RootOutputPath+"trd_singleTrackHits_Run_004375.root", TString tag = "test_4349_thr150") {
    gStyle->SetOptStat(0);
    ROOT::EnableImplicitMT();
    
    auto df_wRadiators = ROOT::RDataFrame(tree_name.Data(), evt_root_file_wRadiator.Data());

    auto df_wRadiators1 = df_wRadiators.Define("extrp_ych", GetTRDHitsExtrpYCh, {"extrp_y"});
    // auto df_wRadiators2 = df_wRadiators1.Define("is_match_hit", GetIsMatchHit, {"extrp_y", "v_trd_y", "trk_ytime_coincidence"});
    auto df_wRadiators2 = df_wRadiators1.Define("ydiff", GETYChExtrpYChDiff, {"extrp_ych", "ych"});
    // PlotDiagnosticTrackers(df_wRadiators2, tag+"wRadiators");
    // PlotDiagnosticTrackers2(df_wRadiators2, tag+"wRadiators");
    // PlotDiagnostics(df_wRadiators2, tag+"wRadiators");
    // PlotDiagnostics2(df_wRadiators2, tag+"wRadiators");

    // PlotDiagnosticsXYTCoin(df_wRadiators2, tag+"wRadiators");
    PlotDiagnosticsXYTCoin2(df_wRadiators2, tag+"wRadiators");

    auto df_wRadiators3 = df_wRadiators2.Define("xymaxamptimediff", GetXYMaxAmpTimeDiff, {"xtime", "ytime"});
    PlotDiagnosticsXYMaxAmp(df_wRadiators3, tag+"wRadiators");
    auto df_wRadiators4 = df_wRadiators3.Filter("TMath::Abs(xymaxamptimediff) < 2");
    PlotDiagnosticsXYMaxAmp(df_wRadiators4, tag+"wRadiatorsXYMaxAmpTimeDiffCut");

    auto df_wRadiators5 = df_wRadiators4.Filter("GEMTrkrsDeltaX < 5");
    // PlotDiagnostics(df_wRadiators5, tag+"wRadiatorsDeltaXCut");
    // PlotDiagnostics2(df_wRadiators5, tag+"wRadiatorsDeltaXCut");

    auto df_wRadiator6 = df_wRadiators5.Filter("TMath::Abs(ydiff) < 5");
    PlotDiagnosticsXYMaxAmp(df_wRadiator6, tag+"wRadiatorsXYMaxAmpTimeDiffCutDeltaXCutYChDiffCut");
    PlotDiagnosticsXYTCoin(df_wRadiator6, tag+"wRadiatorsXYMaxAmpTimeDiffCutDeltaXCutYChDiffCut");
    PlotDiagnosticsXYTCoin2(df_wRadiator6, tag+"wRadiatorsXYMaxAmpTimeDiffCutDeltaXCutYChDiffCut");

    auto df_wRadiator7 = df_wRadiator6.Define("xpos_xch_diff", GetXPosXChDiff, {"xpos", "xch"});
    auto df_wRadiator8 = df_wRadiator7.Define("ypos_ych_diff", GetYPosYChDiff, {"ypos", "ych"});

    PlotDiagnosticsPulseMaxDiff(df_wRadiator8, tag+"wRadiatorsXYMaxAmpTimeDiffCutDeltaXCutYChDiffCut");

    // auto df_wRadiator9 = df_wRadiator8.Filter(FilterPulseMaxDiff, {"xpos_xch_diff", "ypos_ych_diff"});
    // PlotDiagnosticsPulseMaxDiff(df_wRadiator9, tag+"wRadiatorsXYMaxAmpTimeDiffCutDeltaXCutYChDiffCutPulseMaxDiffCut");
    // PlotDiagnosticsXYTCoin(df_wRadiator9, tag+"wRadiatorsXYMaxAmpTimeDiffCutDeltaXCutYChDiffCutPulseMaxDiffCut");
    PlotDiagnosticsXYTCoin4(df_wRadiator8, tag+"wRadiatorsXYMaxAmpTimeDiffCutDeltaXCutYChDiffCutPulseMaxDiffCut");

    df_wRadiator8.Snapshot(tree_name.Data(), (tag+"_wRadiators_clustering.root").Data(), {"v_pulsecut_xpos", "v_pulsecut_xpos_time", "v_pulsecut_xpos_amp", "v_pulsecut_ypos", "v_pulsecut_ypos_time", "v_pulsecut_ypos_amp"});

    // PlotDiagnosticsXYTCoin3(df_wRadiators3, tag+"wRadiatorsDeltaXCut");


    // PlotDiagnosticTrackers(df_wRadiators3, tag+"wRadiatorsDeltaXCut");
    // PlotDiagnosticTrackers2(df_wRadiators3, tag+"wRadiatorsDeltaXCut");
    // PlotDiagnostics(df_wRadiators3, tag+"wRadiatorsDeltaXCut");
    // PlotDiagnostics2(df_wRadiators3, tag+"wRadiatorsDeltaXCut");

    // auto df_wRadiator4 = df_wRadiators3.Filter("TMath::Abs(ydiff) < 20");
    // PlotDiagnosticTrackers(df_wRadiator4, tag+"wRadiatorsDeltaXCutYChDiffCut");
    // PlotDiagnosticTrackers2(df_wRadiator4, tag+"wRadiatorsDeltaXCutYChDiffCut");
    // PlotDiagnostics(df_wRadiator4, tag+"wRadiatorsDeltaXCutYChDiffCut");
    // PlotDiagnostics2(df_wRadiator4, tag+"wRadiatorsDeltaXCutYChDiffCut");
    // PlotDiagnostics3(df_wRadiator4, tag+"wRadiatorsDeltaXCutYChDiffCut");

    // // define new collumn vector xtime - ytime
    // auto df_wRadiator5 = df_wRadiator4.Define("xytimediff", GetXYTimeDiff, {"zpos", "ypos_time"});
    // df_wRadiator5 = df_wRadiator5.Define("xpos_filtered_time", GetXFilteredTime, {"xpos", "zpos", "ypos_time"});
    // df_wRadiator5 = df_wRadiator5.Define("ypos_filtered_time", GetYFilteredTime, {"ypos", "zpos", "ypos_time"});

    // PlotDiagnostics4(df_wRadiator5, tag+"wRadiatorsDeltaXCutYChDiffCut");
    

    // auto df_wRadiators3 = df_wRadiators2.Filter("is_match_hit == 1");
    // PlotDiagnostics(df_wRadiators3, "wRadiatorsMatchHit");
    // PlotDiagnostics2(df_wRadiators3, "wRadiatorsMatchHit");

    // auto df_wRadiators5 = df_wRadiators4.Filter("ych > 160");
    // PlotDiagnostics(df_wRadiators5, "wRadiatorsMatchHitDeltaXCutYChCut");
    // PlotDiagnostics2(df_wRadiators5, "wRadiatorsMatchHitDeltaXCutYChCut");

    // auto df_wRadiators6 = df_wRadiators5.Filter("dedx > 600");
    // PlotDiagnostics(df_wRadiators6, "wRadiatorsMatchHitDeltaXCutYChCutYampCut");
    // PlotDiagnostics2(df_wRadiators6, "wRadiatorsMatchHitDeltaXCutYChCutYampCut");

    // auto df_noRadiators = ROOT::RDataFrame(tree_name.Data(), evt_root_file_noRadiator.Data());
    // auto df_noRadiators1 = df_noRadiators.Define("extrp_ych", GetTRDHitsExtrpYCh, {"extrp_y"});
    // auto df_noRadiators2 = df_noRadiators1.Define("ydiff", GETYChExtrpYChDiff, {"extrp_ych", "ych"});
    // PlotDiagnosticTrackers(df_noRadiators2, "noRadiators");
    // PlotDiagnosticTrackers2(df_noRadiators2, "noRadiators");
    // PlotDiagnostics(df_noRadiators2, "noRadiators");
    // PlotDiagnostics2(df_noRadiators2, "noRadiators");

    // auto df_noRadiators3 = df_noRadiators2.Filter("GEMTrkrsDeltaX < 5");
    // PlotDiagnosticTrackers(df_noRadiators3, "noRadiatorsDeltaXCut");
    // PlotDiagnosticTrackers2(df_noRadiators3, "noRadiatorsDeltaXCut");
    // PlotDiagnostics(df_noRadiators3, "noRadiatorsDeltaXCut");
    // PlotDiagnostics2(df_noRadiators3, "noRadiatorsDeltaXCut");

    // auto df_noRadiator4 = df_noRadiators3.Filter("TMath::Abs(ydiff) < 20");
    // PlotDiagnosticTrackers(df_noRadiator4, "noRadiatorsDeltaXCutYChDiffCut");
    // PlotDiagnosticTrackers2(df_noRadiator4, "noRadiatorsDeltaXCutYChDiffCut");
    // PlotDiagnostics(df_noRadiator4, "noRadiatorsDeltaXCutYChDiffCut");
    // PlotDiagnostics2(df_noRadiator4, "noRadiatorsDeltaXCutYChDiffCut");
    // PlotDiagnostics3(df_noRadiator4, "noRadiatorsDeltaXCutYChDiffCut");

}