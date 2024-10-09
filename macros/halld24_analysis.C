
TString rootfiles_dir = "/work/halld2/home/nseptian/halld24_TestBeam/RootOutput/halld24/";

TFile* LoadRootFile(TString filename);
int halld24_analysis(){
    
    TString run_number = "004324";
    TString filename = "Run_" + run_number + "_Output.root";
    TFile* file = LoadRootFile(filename);
    if (file == NULL) return 1;

    
    TCanvas *c1 = new TCanvas("c1","c1",800,600);

    TString h1_name = "f125_el_Ymax";

    TH1F* h1 = (TH1F*)file->Get(h1_name);

    h1->Draw();
    TString pdf_name = h1_name + "_Run_" + run_number + ".pdf";
    c1->SaveAs(pdf_name);

    h1_name = "f125_el_Xmax";
    h1 = (TH1F*)file->Get(h1_name);
    h1->Draw();
    pdf_name = h1_name + "_Run_" + run_number + ".pdf";
    c1->SaveAs(pdf_name);

    TString h2_name = "f125_yamp2d";
    TH2F* h2 = (TH2F*)file->Get(h2_name);
    h2->Draw("colz");
    pdf_name = h2_name + "_Run_" + run_number + ".pdf";
    c1->SaveAs(pdf_name);

    // draw projection Y of 2D f125_yamp2d
    h2->ProjectionY()->Draw();
    pdf_name = h2_name + "_ProjectionY_Run_" + run_number + ".pdf";
    c1->SaveAs(pdf_name);

    // draw projection of Ymax in a range
    Double_t dt_range[2] = {100.0, 150};
    h2->ProjectionY("f125_yamp2d_Ymax", h2->GetXaxis()->FindBin(dt_range[0]), h2->GetXaxis()->FindBin(dt_range[1]))->Draw();
    pdf_name = h2_name + "_ProjectionY_DT_"+Form("%.0f_%.0f", dt_range[0], dt_range[1])+"_Run_" + run_number + ".pdf";
    c1->SaveAs(pdf_name);

    return 0;
}

TFile* LoadRootFile(TString filename){
    TFile* file = new TFile(rootfiles_dir+filename);
    if (!file->IsOpen()){
        std::cout << "Error: file " << filename << " not found" << std::endl;
        return NULL;
    }
    return file;
}