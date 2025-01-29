#define trdclass_halld24_cxx
#include "trdclass_halld24.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TSocket.h"
#include "TMarker.h"
#include "TMultiGraph.h"
#include "PlotLib.C"
#include <TError.h>
#include <TStopwatch.h>
#include <iostream>
#include <fstream>
#include "GNN/gnn_model.h"
#include "GNN/gnn_model.cpp"
#include "GNN/toGraph.cpp"

#define NPRT 1000
#define USE_TRK
#define MAX_PRINT 1
//#define USE_250_PULSE
#define USE_125_RAW
#define BUFSIZE 128000
#define DEBUG 0
#define MAX_CLUST 500
#define MAX_NODES 100
#define USE_GNN    1
#define USE_CLUST  1
#define USE_FIT    1
#define SAVE_TRACK_HITS
//#define WRITE_CSV
#define SAVE_PDF

//-- For single evt clustering display, uncomment BOTH:
// #define SHOW_EVT_DISPLAY
// #define SHOW_EVTbyEVT

void WriteToCSV(std::ofstream &csvFile, float v1, float v2, float v3, float v4, float v5) {
  csvFile<<v1<<","<<v2<<","<<v3<<","<<v4<<","<<v5<<std::endl;
}

//=================================================
//            Detector Channel Mapping
//=================================================

// -- Large GEM-TRD fADC Mapping--
#define first_slot 3
#define gem_x_slot 0
#define gem_x_ch0 0
#define gem_y_slot 1
#define gem_y_ch0 48
#define last_slot 8
#define last_ch 71

int GetGEMXChan(int fch, int slot) {
  int sl=slot-first_slot;
  if((sl>=gem_x_slot&&sl<gem_y_slot)||(sl==gem_y_slot&&fch<gem_y_ch0)) {
    int ch=sl*72+fch-gem_x_ch0;
    return ch;
  } else { return -1; }
}

int GetGEMYChan(int fch, int slot) {
  int sl=slot-first_slot;
  if(sl>7)sl=sl-2;
  if((sl>gem_y_slot&&sl<last_slot)||(sl==gem_y_slot&&fch>=gem_y_ch0)||(sl==last_slot&&fch<=last_ch)) {
    int ch=(sl-gem_y_slot)*72+fch-gem_y_ch0;
    return ch;
  } else { return -1; }
}

//-- GEM-Tracker 1 & 2 SRS Mapping
int GetTrackerChan(int chn) {
  int apvn=chn/128;
  int apvch=127-(chn-apvn*128);
  int ch = apvch+128*apvn;
  return ch;
}

double GetTrackersDeltaX(double GEMTrkX1, double GEMTrkX2) {
  double p0 = -40.8951;
  double p1 = 1.28738;

  double deltaX1 = GEMTrkX1 - (GEMTrkX2-p0)/p1;
  double deltaX2 = GEMTrkX2 - (p0 + p1*GEMTrkX1);

  return std::sqrt(deltaX1*deltaX1 + deltaX2*deltaX2);
}

//=========== End Detector Channel Mapping =============

void trdclass_halld24::Loop() {
  //   In a ROOT session, you can do:
  //      root> .L trdclass_halld24.C
  //      root> trdclass_halld24 t(RunNum)
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries
  
  if (fChain == 0) return;
  
  //-- Create output .root file to store all histograms
  TFile* fOut;
  char rootFileName[256]; sprintf(rootFileName, "RootOutput/halld24/Run_%06d_Output.root", RunNum);
  fOut = new TFile(rootFileName, "RECREATE");
  fOut->cd();
  #ifdef WRITE_CSV
    char csvTitle[100];
    sprintf(csvTitle,"EventByEvent_Run00%d.csv",RunNum);
    std::ofstream csvFile(csvTitle);
  #endif
  TList *HistList = new TList();
  int THRESH=600; //-- GEM-TRD ADC threshold
  
  //========================================================================
  //            B o o k    H i s t o g r a m s
  //========================================================================
  
  hcount= new TH1D("hcount","Count",3,0,3);
  hcount->SetStats(0);  hcount->SetFillColor(38); hcount->SetMinimum(1.);
  #if ROOT_VERSION_CODE > ROOT_VERSION(6,0,0)
    hcount->SetCanExtend(TH1::kXaxis);
  #else
    hcount->SetBit(TH1::kCanRebin);
  #endif
  //-- Large GEM-TRD fADC Amplitude Distributions
  f125_el_x = new TH1F("f125_el_x","GEM-TRD f125 Pulse Peak Amp in X; ADC Amplitude ; Counts ",102,0.5,4096.5);
  f125_el_y = new TH1F("f125_el_y","GEM-TRD f125 Pulse Peak Amp in Y; ADC Amplitude ; Counts ",102,0.5,4096.5);
  f125_el_Xmax = new TH1F("f125_el_Xmax","GEM-TRD f125 Max Pulse Amp in X; ADC Amplitude ; Counts ",102,0.5,4096.5);
  f125_el_Ymax = new TH1F("f125_el_Ymax","GEM-TRD f125 Max Pulse Amp in Y; ADC Amplitude ; Counts ",102,0.5,4096.5);
  f125_el_xch_max = new TH1F("f125_el_xch_max","GEM-TRD f125 X Channel of Max Pulse; GEMTRD Channel (X) ; Counts ",128,-0.5,127.5);
	f125_el_ych_max = new TH1F("f125_el_ych_max","GEM-TRD f125 Y Channel of Max Pulse; GEMTRD Channel (Y) ; Counts ",528,-0.5,527.5);
	
  //-- GEM-TRD ADC Amplitude info in Time (2D)
  TH2F *f125_ydiff2d = new TH2F("f125_ydiff2d","GEM-TRD Y Pulse Amp in Time ; Time Response (8ns) ; GEMTRD Channel (Y) ",250,0.5,250.5,400,-50.,50.);  
  TH2F *f125_yamp = new TH2F("f125_yamp","GEM-TRD Y Pulse Amp in Time (w Hit Match Cond.); Time Response (8ns) ; GEMTRD Channel (Y) ",250,0.5,250.5,528,-0.5,527.5);  
  f125_yamp2d = new TH2F("f125_yamp2d","GEM-TRD Amp. in Time (w Time Coinc. Cond.); Time Response (8ns) ; GEMTRD Channel (Y) ",250,0.5,250.5,528,-0.5,527.5);
  f125_xamp2d = new TH2F("f125_xamp2d","GEM-TRD Amp. in Time ; Time Response (8ns) ; GEMTRD Channel (X) ",250,0.5,250.5,128,-0.5,127.5);
  
  //-- GEM Trackers (SRS)
  hgemtrkr1_peak_xy = new TH2F("hgemtrkr1_peak_xy","GEM-TRKR1 Pulse Peak X-Y Correlation; Pulse Peak X [mm]; Pulse Peak Y [mm] ",512,-0.5,511.5,128,-0.5,127.5);
	hgemtrkr2_peak_xy = new TH2F("hgemtrkr2_peak_xy","GEM-TRKR2 Pulse Peak X-Y Correlation; Pulse Peak X [mm]; Pulse Peak Y [mm] ",512,-0.5,511.5,128,-0.5,127.5);
	hgemtrkr1_max_x = new TH1F("hgemtrkr1_max_x","GEM-TRKR1 Max Pulse Amp in X; Max ADC Amplitude ",102,0.5,4096.5);
	hgemtrkr2_max_x = new TH1F("hgemtrkr2_max_x","GEM-TRKR2 Max Pulse Amp in X; Max ADC Amplitude ",102,0.5,4096.5);
	hgemtrkr_double_y = new TH2F("hgemtrkr_double_y","GEM-TRKR1 & GEM-TRKR2 Hit Correlation in Y; GEMTRKR1 Pulse Peak Y [mm]; GEMTRKR2 Pulse Peak Y [mm] ",128,-0.5,127.5,128,-0.5,127.5);
	hgemtrkr_double_x = new TH2F("hgemtrkr_double_x","GEM-TRKR1 & GEM-TRKR2 Hit Correlation in X; GEMTRKR1 Pulse Peak X [mm]; GEMTRKR2 Pulse Peak X [mm] ",512,-0.5,511.5,512,-0.5,511.5);
  hgemtrkr1_peak_x = new TH1F("hgemtrkr1_peak_x"," GEM-TRKR1 Pulse Distribution in X ; X [mm] ",512,-0.5,511.5);
  hgemtrkr1_peak_y = new TH1F("hgemtrkr1_peak_y"," GEM-TRKR1 Pulse Distribution in Y ; Y [mm] ",128,-0.5,127.5);
	hgemtrkr2_peak_x = new TH1F("hgemtrkr2_peak_x"," GEM-TRKR2 Pulse Distribution in X ; X [mm] ",512,-0.5,511.5);
  hgemtrkr2_peak_y = new TH1F("hgemtrkr2_peak_y"," GEM-TRKR2 Pulse Distribution in Y ; Y [mm] ",128,-0.5,127.5);
  hgemtrkr_peak_delta_x = new TH1F("hgemtrkr_peak_delta_x"," GEM-TRKR1 & GEM-TRKR2 Delta X; Delta X [mm] ",41,-0.5,40.5);
  
  //-- LUBOMIR'S Track Matching Histograms
  TH1F *thits = new TH1F("thits","",528.,-0.5,527.5);
  //TH1F *tnohits = new TH1F("tnohits","",528.,-0.5,527.5);
  TH1F *teff = new TH1F("teff","",528.,-0.5,527.5);
  TH1F *teff2d = new TH1F("teff2d","",528.,-0.5,527.5);
  //TH1F *thits1 = new TH1F("thits1","",512.,-0.5,511.5);
  //TH1F *tnohits1 = new TH1F("tnohits1","",512.,-0.5,511.5);
  //TH1F *teff1 = new TH1F("teff1","",512.,-0.5,511.5);
  //TH1F *thits2 = new TH1F("thits2","",512.,-0.5,511.5);
  //TH1F *tnohits2 = new TH1F("tnohits2","",512.,-0.5,511.5);
  //TH1F *teff2 = new TH1F("teff2","",512.,-0.5,511.5);
  // float x1=2057.4+10.;
  // float x2=2654.3+10.;
  // float x3=4749.8+10.;
  // float dy1=104.8+20.;
  // float dy2=172.5+20.;
  // float dy3=202.;

  
  // Detectors coordinate from Survey

  float y0=80.59990;
  float v0=104.69986;
  float x0=386.76755;

  float y1=(81.18291+80.72291)/2.;
  float v1=(104.76201+104.76228)/2.;
  float x1=(388.82626+388.83477)/2.;

  float y2=(81.25059+80.79095)/2.;
  float v2=(104.76775+104.77306)/2.;
  float x2=(389.41947+389.43687)/2.;

  float y3=(81.32000+80.77612+81.32446+80.78202)/4.;
  float v3=(105.11588+105.11209+104.38180+104.37675)/4.;
  float x3=(391.54396+391.53969+391.54030+391.53604)/4.;

  float Bl = 1.6617 * 36*0.0254; // 1.6T x 36" field integral


  float dx1=(x1-x0)*1000.;
  float dy1=(y1-y0)*1000.-256.*0.8-31.75; 
  float dx2=(x2-x0)*1000.;
  float dy2=(y2-y0)*1000.-256.*0.8-31.75-2.4;
  float dx3=(x3-x0)*1000.+31.35;
  float dy3=(y3-y0)*1000.-264.;  

  
  //-- Clustering and track finding histos
  gErrorIgnoreLevel = kBreak; // Suppress warning messages from empty chi^2 fit data
  TF1 fx1("fx1","pol1",40,150);
  f125_fit = new TH2F("f125_fit","GEM-TRD Track Fit in Time; Time Response (8ns) ; GEMTRD Channel (Y)",250,0.5,250.5,528,-0.5,527.5);
  //----------------- TCanvas c2 for Sergey's FPGA Display ----------
  #ifdef SHOW_EVT_DISPLAY
    char c2Title[256]; sprintf(c2Title,"Event_Display_Run=%d",RunNum);
    TCanvas *c2 = new TCanvas("FPGA",c2Title,1000,100,1500,1300);
    c2->Divide(3,2); c2->cd(1);
  #endif
  int nx0=100;       int ny0=528;
  double Ymin=0.;    double Ymax=528.;
  double Xmin=0.;    double Xmax=30.;
  hevt  = new TH2F("hevt","GEM-TRD Event display; z pos [mm]; y pos [mm]",nx0,Xmin,Xmax,ny0,Ymin,Ymax); hevt->SetStats(0); hevt->SetMaximum(10.);
  hevtc = new TH2F("hevtc","Clustering; FADC bins; GEM strips",nx0,-0.5,nx0-0.5,ny0,-0.5,ny0-0.5); hevtc->SetStats(0); hevtc->SetMinimum(0.07); hevtc->SetMaximum(40.);
  hevtf = new TH2F("hevtf","Clusters for FPGA; z pos [mm]; y pos [mm]",nx0,Xmin,Xmax,ny0,Ymin,Ymax);  hevtf->SetStats(0); hevtf->SetMaximum(10.);
  //f125_el_Yraw = new TH2F("f125_el_Yraw","GEM-TRD Raw Y Response; Time Response (8ns); Channel (Y)",250,0.5,250.5,528,-0.5,527.5);
  gem_zHist = new TH1F("gem_zHist","gem_zHist",20,80.,200.);
  
  //============= End Book Histograms =================
  
  int gem_xch_max;
  int gem_ych_max;
  int gemtrkr1_xch_max;
  int gemtrkr1_ych_max;
  int gemtrkr2_xch_max;
  int gemtrkr2_ych_max;
  double gem_xamp_max;
  double gem_yamp_max;
  double gemtrkr1_xamp_max;
  double gemtrkr1_yamp_max;
  double gemtrkr2_xamp_max;
  double gemtrkr2_yamp_max;
  double chi2cc_gem;
  double gem_integral;
  double a0, a1;
  double a0nn, a1nn;
  double chi2nn_gem;
  int NTRACKS;
  int gem_xtime_max;
  int gem_ytime_max;
  ULong64_t gt1_idx_x, gt1_idx_y;
  ULong64_t gt2_idx_x, gt2_idx_y;
  double GEMTrkrsDeltaX;
  double extrp_y;
  vector<double> v_trd_y;
  vector<bool> v_trk_xtime_coincidence;
  vector<bool> v_trk_ytime_coincidence;

  //=================================================================
  //  Create GEM-TRD TTree of Hit Info for NN Rej Factor Calculation
  //=================================================================

  TFile* fHits; //-- Create output .root file to store GEM-TRD TTree for passing to NN
  #ifdef SAVE_TRACK_HITS
    char hitsFileName[256]; sprintf(hitsFileName, "RootOutput/halld24/trd_singleTrackHits_Run_%06d.root", RunNum);
    fHits = new TFile(hitsFileName, "RECREATE");
    EVENT_VECT_GEM = new TTree("gem_hits","GEM TTree with single track hit info");
    EVENT_VECT_GEM->Branch("event_num",&event_num,"event_num/I");
    EVENT_VECT_GEM->Branch("nhit",&gem_nhit,"gem_nhit/I");
    EVENT_VECT_GEM->Branch("xpos",&gem_xpos);
    EVENT_VECT_GEM->Branch("ypos",&gem_ypos);
    EVENT_VECT_GEM->Branch("ypos_amp",&gem_ypos_amp);
    EVENT_VECT_GEM->Branch("ypos_time",&gem_ypos_time);
    EVENT_VECT_GEM->Branch("zpos",&gem_zpos);
    EVENT_VECT_GEM->Branch("dedx",&gem_dedx);
    EVENT_VECT_GEM->Branch("zHist",&gem_zHist_vect);
    EVENT_VECT_GEM->Branch("xposc",&clu_xpos);
    EVENT_VECT_GEM->Branch("yposc",&clu_ypos);
    EVENT_VECT_GEM->Branch("zposc",&clu_zpos);
    EVENT_VECT_GEM->Branch("dedxc",&clu_dedx);
    EVENT_VECT_GEM->Branch("widthc",&clu_width);
    EVENT_VECT_GEM->Branch("xch",&gem_xch_max);
    EVENT_VECT_GEM->Branch("ych",&gem_ych_max);
    EVENT_VECT_GEM->Branch("xamp",&gem_xamp_max);
    EVENT_VECT_GEM->Branch("yamp",&gem_yamp_max);
    EVENT_VECT_GEM->Branch("xtime",&gem_xtime_max);
    EVENT_VECT_GEM->Branch("ytime",&gem_ytime_max);
    EVENT_VECT_GEM->Branch("xchtrkr1",&gemtrkr1_xch_max);
    EVENT_VECT_GEM->Branch("ychtrkr1",&gemtrkr1_ych_max);
    EVENT_VECT_GEM->Branch("xchtrkr2",&gemtrkr2_xch_max);
    EVENT_VECT_GEM->Branch("ychtrkr2",&gemtrkr2_ych_max);
    EVENT_VECT_GEM->Branch("xamptrkr1",&gemtrkr1_xamp_max);
    EVENT_VECT_GEM->Branch("yamptrkr1",&gemtrkr1_yamp_max);
    EVENT_VECT_GEM->Branch("xamptrkr2",&gemtrkr2_xamp_max);
    EVENT_VECT_GEM->Branch("yamptrkr2",&gemtrkr2_yamp_max);
    EVENT_VECT_GEM->Branch("nxtrkr1",&gt1_idx_x);
    EVENT_VECT_GEM->Branch("nytrkr1",&gt1_idx_y);
    EVENT_VECT_GEM->Branch("nxtrkr2",&gt2_idx_x);
    EVENT_VECT_GEM->Branch("nytrkr2",&gt2_idx_y);
    EVENT_VECT_GEM->Branch("chi2nn",&chi2nn_gem);
    EVENT_VECT_GEM->Branch("chi2",&chi2cc_gem);
    EVENT_VECT_GEM->Branch("Fint",&gem_integral);
    EVENT_VECT_GEM->Branch("a0",&a0);
    EVENT_VECT_GEM->Branch("a1",&a1);
    EVENT_VECT_GEM->Branch("a0nn",&a0nn);
    EVENT_VECT_GEM->Branch("a1nn",&a1nn);
    EVENT_VECT_GEM->Branch("ntr",&NTRACKS,"NTRACKS/I");
    EVENT_VECT_GEM->Branch("GEMTrkrsDeltaX",&GEMTrkrsDeltaX);

    // clone the GEM-TRD TTree for the GEM Trackers with new name
    EVENT_VECT_GEM_TRACKERS = (TTree*)EVENT_VECT_GEM->CloneTree(0);
    EVENT_VECT_GEM_TRACKERS->SetName("gem_trackers_hits");
    EVENT_VECT_GEM_TRACKERS->SetTitle("GEM Trackers TTree with single track hit info");
    EVENT_VECT_GEM_TRACKERS->SetDirectory(fHits);
    EVENT_VECT_GEM_TRACKERS->Branch("GEMTrkrNHitMatch",&gem_trk_hit);
    // TRD hit position extrappolated from GEM Trackers
    EVENT_VECT_GEM_TRACKERS->Branch("extrp_y", &extrp_y);
    // TRD hit position relative to the beam line
    EVENT_VECT_GEM_TRACKERS->Branch("v_trd_y", &v_trd_y);
    EVENT_VECT_GEM_TRACKERS->Branch("trk_xtime_coincidence", &v_trk_xtime_coincidence);
    EVENT_VECT_GEM_TRACKERS->Branch("trk_ytime_coincidence", &v_trk_ytime_coincidence);
  #endif
  
  //==============================================================================
  //                      E v e n t    L o o p
  //==============================================================================
  
  TStopwatch timer;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  if (MaxEvt>0) nentries=MaxEvt;  //-- limit number of events for test
  Long64_t jentry=0;
  printf("***>>>  Begin Event Loop - 1st evt=%lld, Last evt=%lld \n",FirstEvt,MaxEvt);
  timer.Start();
  
  for (jentry=FirstEvt; jentry<nentries; jentry++) { //-- Event Loop --
    event_num = jentry;
    Count("EVT");
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (jentry<MAX_PRINT || !(jentry%1000)) {
      printf("------- evt=%llu  f125_wraw_count=%llu f125_pulse_count=%llu f250_wraw_count=%llu, gem_peak_count=%llu \n", jentry, f125_wraw_count, f125_pulse_count, f250_wraw_count, gem_peak_count);
    }
    
    //-- Reset TTree Hit Info for each event
    gem_nhit=0;
    clu_nhit=0;
    trkr_hit=0;
    gem_xpos.clear();
    gem_ypos.clear();
    gem_ypos_amp.clear();
    gem_ypos_time.clear();
    gem_zpos.clear();
    gem_dedx.clear();
    gem_zHist->Reset();
    gem_zHist_vect.clear();
    clu_xpos.clear();
    clu_ypos.clear();
    clu_zpos.clear();
    clu_dedx.clear();
    clu_width.clear();
    f125_fit->Reset();
    chi2cc_gem=-999.;
    chi2nn_gem=-999.;
    gem_integral=0.;
    gem_yamp_max =-1;
    gem_ych_max=-1;
    gem_ytime_max=-1;
    gem_xamp_max =-1;
    gem_xch_max=-1;
    gem_xtime_max=-1;
    GEMTrkrsDeltaX = 99999;
    gem_trk_hit=0;
    extrp_y=0;
    v_trd_y.clear();
    v_trk_xtime_coincidence.clear();
    v_trk_ytime_coincidence.clear();
    
    int gem_xtime[f125_pulse_count];
    int gem_ytime[f125_pulse_count];
    int nxpulse=0;
    int nypulse=0;
    
    //======================================================================
    //        Process SRS Information
    //======================================================================
    
    //===========================================================
    //    GEMTracker (SRS) Distributions
    //===========================================================
    
    double gemtrkr1_peak_pos_y[gem_peak_count];
    double gemtrkr1_peak_pos_x[gem_peak_count];
    double gemtrkr1_peak_height_y[gem_peak_count];
    double gemtrkr1_peak_height_x[gem_peak_count];
    int gemtrkr1_peak_ch_y[gem_peak_count];
    int gemtrkr1_peak_ch_x[gem_peak_count];
    double gemtrkr2_peak_pos_y[gem_peak_count];
    double gemtrkr2_peak_pos_x[gem_peak_count];
    double gemtrkr2_peak_height_y[gem_peak_count];
    double gemtrkr2_peak_height_x[gem_peak_count];
    int gemtrkr2_peak_ch_y[gem_peak_count];
    int gemtrkr2_peak_ch_x[gem_peak_count];
    gt1_idx_x=0.;
    gt1_idx_y=0.;
    gt2_idx_x=0.;
    gt2_idx_y=0.;
    gemtrkr1_xch_max=-1.;
    gemtrkr1_ych_max=-1.;
    gemtrkr1_xamp_max=0.;
    gemtrkr1_yamp_max=0.;
    gemtrkr2_xch_max=-1.;
    gemtrkr2_ych_max=-1.;
    gemtrkr2_xamp_max=0.;
    gemtrkr2_yamp_max=0.;
    
    for (ULong64_t i=0; i<gem_peak_count; i++) { //-- SRS Pulse Peak Info Loop
      //-- Match SRS info to correct detectors and detector planes
      //-- GEM Tracker 1
      if (gem_peak_plane_name->at(i) == "GEMTR1Y") {
        gemtrkr1_peak_pos_y[gt1_idx_y] = GetTrackerChan(gem_peak_index->at(i));
        gemtrkr1_peak_height_y[gt1_idx_y] = gem_peak_height->at(i);
        gt1_idx_y++;
      } if (gem_peak_plane_name->at(i) == "GEMTR1X") {
        gemtrkr1_peak_pos_x[gt1_idx_x] = GetTrackerChan(gem_peak_index->at(i));
        gemtrkr1_peak_height_x[gt1_idx_x] = gem_peak_height->at(i);
        gt1_idx_x++;
      }
      //-- GEM Tracker 2
      if (gem_peak_plane_name->at(i) == "GEMTR2Y") {
				gemtrkr2_peak_pos_y[gt2_idx_y] = GetTrackerChan(gem_peak_index->at(i));
        gemtrkr2_peak_height_y[gt2_idx_y] = gem_peak_height->at(i);
        gt2_idx_y++;
			} if (gem_peak_plane_name->at(i) == "GEMTR2X") {
        gemtrkr2_peak_pos_x[gt2_idx_x] = GetTrackerChan(gem_peak_index->at(i));
        gemtrkr2_peak_height_x[gt2_idx_x] = gem_peak_height->at(i);
        gt2_idx_x++;
			}
    } //-- End SRS Peak Loop

    // reject noise using trackers correlation
    if (gt1_idx_x==1 && gt2_idx_x==1) {
      trkr_hit=1;
      GEMTrkrsDeltaX = GetTrackersDeltaX(gemtrkr1_peak_pos_x[0], gemtrkr2_peak_pos_x[0]);
      }
    else {
      gt1_idx_x = 0;
      gt2_idx_x = 0;
    }
    
    //-- Fill GEMTracker Histos
    if (gt1_idx_x>0) {
      for (ULong64_t j=0; j<gt1_idx_x; j++) {
        hgemtrkr1_peak_x->Fill(gemtrkr1_peak_pos_x[j]);
        if(gemtrkr1_peak_height_x[j]>gemtrkr1_xamp_max){
           gemtrkr1_xamp_max=gemtrkr1_peak_height_x[j];
           gemtrkr1_xch_max=gemtrkr1_peak_pos_x[j];
        }
        if (gt1_idx_y>0) {
        for (ULong64_t k=0; k<gt1_idx_y; k++) {
          hgemtrkr1_peak_xy->Fill(gemtrkr1_peak_pos_x[k], gemtrkr1_peak_pos_y[j]);
        	}
        }
				if (gt2_idx_x>0) {
        	for (ULong64_t i=0; i<gt2_idx_x; i++) {
            hgemtrkr_double_x->Fill(gemtrkr1_peak_pos_x[j], gemtrkr2_peak_pos_x[i]);
            hgemtrkr_peak_delta_x->Fill(GEMTrkrsDeltaX);
        	}
        }
      }
   	}
    
    if (gt2_idx_x>0) {
      for (ULong64_t j=0; j<gt2_idx_x; j++) {
        hgemtrkr2_peak_x->Fill(gemtrkr2_peak_pos_x[j]);
        if(gemtrkr2_peak_height_x[j]>gemtrkr2_xamp_max){
           gemtrkr2_xamp_max=gemtrkr2_peak_height_x[j];
           gemtrkr2_xch_max=gemtrkr2_peak_pos_x[j];
        }
				if (gt2_idx_y>0) {
					for (ULong64_t k=0; k<gt2_idx_y; k++) {
						hgemtrkr_double_y->Fill(gemtrkr1_peak_pos_y[j], gemtrkr2_peak_pos_y[k]);
					}
				}
      }
    }
		
		if (gt2_idx_x>0) {
      for (ULong64_t j=0; j<gt2_idx_x; j++) {
        hgemtrkr2_peak_x->Fill(gemtrkr2_peak_pos_x[j]);
        if(gemtrkr2_peak_height_x[j]>gemtrkr2_xamp_max) {
           gemtrkr2_xamp_max=gemtrkr2_peak_height_x[j];
           gemtrkr2_xch_max=gemtrkr2_peak_pos_x[j];
        }
        if (gt2_idx_y>0) {
        for (ULong64_t k=0; k<gt2_idx_y; k++) {
          hgemtrkr2_peak_xy->Fill(gemtrkr2_peak_pos_x[k], gemtrkr2_peak_pos_y[j]);
        	}
        }
      }
   	}
		
    if (gt2_idx_y>0) {
      for (ULong64_t j=0; j<gt2_idx_y; j++) {
        hgemtrkr2_peak_y->Fill(gemtrkr2_peak_pos_y[j]);
        if(gemtrkr2_peak_height_y[j]>gemtrkr2_yamp_max) {
           gemtrkr2_yamp_max=gemtrkr2_peak_height_y[j];
           gemtrkr2_ych_max=gemtrkr2_peak_pos_y[j];
	      }
			}
		}
    
		if(gemtrkr1_xamp_max>0) hgemtrkr1_max_x->Fill(gemtrkr1_xamp_max);
		if(gemtrkr2_xamp_max>0) hgemtrkr2_max_x->Fill(gemtrkr2_xamp_max);
    
    //============== End Process SRS Information ===================
    
    //=================================================================
    //            Process fADC125 Pulse Information (GEM-TRD)
    //=================================================================
    
for (ULong64_t i=0; i<f125_pulse_count; i++) {
      
      float peak_amp = f125_pulse_peak_amp->at(i);
      float ped = f125_pulse_pedestal->at(i);
      if (0 > ped || ped > 200 ) ped = 100;
      float amp = peak_amp-ped;
      if (amp<0) amp=0;
      float time = f125_pulse_peak_time->at(i);
      int fADCSlot = f125_pulse_slot->at(i);
      int fADCChan = f125_pulse_channel->at(i);
      int gemChanX = GetGEMXChan(fADCChan, fADCSlot);
      int gemChanY = GetGEMYChan(fADCChan, fADCSlot);
      
      // store gem-trd timing and track information
      if (gemChanY>-1 && amp>THRESH) {
        gem_ypos.push_back(gemChanY);
        gem_ypos_amp.push_back(amp);
        gem_ypos_time.push_back(time);
        if (gem_yamp_max<amp) {
          gem_yamp_max=amp;
          gem_ych_max=gemChanY;
          gem_ytime_max=time;
        }
        gem_ytime[nypulse]=time;
        nypulse++;
        //-- Fill Histo to Perform Chi^2 Track Fitting On
        f125_fit->Fill(time,gemChanY,amp); //only has two entries maybe should adjust the threshold
      }
      
      if (gemChanX>-1 && amp>THRESH) {
        gem_xpos.push_back(gemChanX);
        gem_dedx.push_back(amp);
        gem_zpos.push_back(time);
        gem_zHist->Fill(time,amp);
        gem_nhit++;
        if (gem_xamp_max<amp) {
          gem_xamp_max=amp;
          gem_xch_max=gemChanX;
          gem_xtime_max=time;
        }
        gem_xtime[nxpulse]=time;
        nxpulse++;
      }
    } //====== End First Fadc125 Pulse Loop ======
    
    //-- Fill GEM-TRD Max Distribution Histos
    if(gem_xamp_max>0) f125_el_Xmax->Fill(gem_xamp_max);
    if(gem_yamp_max>0) f125_el_Ymax->Fill(gem_yamp_max);
		if(gem_xch_max>0) f125_el_xch_max->Fill(gem_xch_max);
		if(gem_ych_max>0) f125_el_ych_max->Fill(gem_ych_max);
    for (int i=1; i<21; i++) {
      gem_zHist_vect.push_back(gem_zHist->GetBinContent(i));
    }
    //-- Perform Chi^2 Fit on Track Left in GEM-TRD
    if (f125_fit->GetEntries()!=0) {
      std::pair<Double_t, Double_t> fitResult = TrkFit(f125_fit,fx1,"fx1",1);
      chi2cc_gem = fitResult.first;
      gem_integral = fitResult.second;
    }
    a0 = fx1.GetParameter(0);
    a1 = fx1.GetParameter(1);
    
    //---- Perform Another Fadc125 Pulse Loop with Hit-Matching Condition with GEM-Trackers Applied ----
    
    if (gt1_idx_x==1&&gt2_idx_x==1) {
      Count("trk_hit");
      //-- Convert GEMTracker channels to distance (mm)
      float y1=dy1+gemtrkr1_xch_max*0.8;
      float y2=dy2+gemtrkr2_xch_max*0.8;
      float a=(y2-y1)/(x2-x1);
      float b=((y1)*x2-(y2)*x1)/(x2-x1);
      float extr = a*x3+b; //-- Track extrapolation line
      extrp_y = extr;
      thits->Fill(extr-dy3);
      bool match = false;
      bool match2d = false;
      
      for (ULong64_t i=0; i<f125_pulse_count; i++) {
      
      float peak_amp = f125_pulse_peak_amp->at(i);
      float ped = f125_pulse_pedestal->at(i);
      if (0 > ped || ped > 200 ) ped = 100;
      float amp = peak_amp-ped;
      if (amp<0) amp=0;
      float time = f125_pulse_peak_time->at(i);
      int fADCSlot = f125_pulse_slot->at(i);
      int fADCChan = f125_pulse_channel->at(i);
      int gemChanX = GetGEMXChan(fADCChan, fADCSlot);
      int gemChanY = GetGEMYChan(fADCChan, fADCSlot);
      
      if (gemChanY>-1 && amp>THRESH ) {
        gem_trk_hit=0;
        float y3=dy3+gemChanY;
        v_trd_y.push_back(y3);
        f125_ydiff2d->Fill(time,y3-extr+a*(time-45)*0.3);
        if (abs(y3-extr)<20) { //-- if actual hit location & track extrapolation difference is within 20 mm
          Count("gem_trk_hit");
          bool tcoin=false;
          for (int it=0; it<nxpulse; it++){
            if (abs(time-gem_xtime[it])<10) tcoin=true; //-- if time between hits is less than 10 time samples, time coincidence is TRUE
          }
          v_trk_ytime_coincidence.push_back(tcoin);
          if (!match) {
            teff->Fill(extr-dy3);
            match = true;
          }
          gem_trk_hit++;
          f125_yamp->Fill(time,gemChanY,amp);
          if (tcoin) {
            f125_yamp2d->Fill(time,gemChanY,amp);
            f125_el_Ymax->Fill(gem_yamp_max);
            f125_el_y->Fill(amp);
            if(!match2d){
              teff2d->Fill(extr-dy3);
              match2d = true;
            }
          }
        }
      }
      if (gemChanX>-1 && amp>THRESH) {
        bool tcoin=false;
        for (int it=0; it<nypulse; it++){
          if (abs(time-gem_ytime[it])<5) tcoin=true;
        }
        v_trk_xtime_coincidence.push_back(tcoin);
        if (tcoin) {
          f125_xamp2d->Fill(time,gemChanX,amp);
          f125_el_Xmax->Fill(gem_xamp_max);
          f125_el_x->Fill(amp);
        }
      }
    } //-- End second fADC125 pulse loop
  } //-- End GEMTracker hit condition 
   
  //============== End Process fADC125 Pulse Information ==================
   
    //=========================================================================
    //            Process Fa125 RAW Data (Sergey's Clustering)
    //=========================================================================
    
    #ifdef USE_125_RAW
      
      hevt->Reset();
      hevtc->Reset();
      hevtf->Reset();
      int hentries=0;
      
      for (ULong64_t i=0; i<f125_wraw_count; i++) { //-- fadc125 raw hit loop
        
        int fadc_window = f125_wraw_samples_count->at(i);
        int fADCSlot = f125_wraw_slot->at(i);
        int fADCChan = f125_wraw_channel->at(i);
				int gemChan = GetGEMYChan(fADCChan, fADCSlot);
				if (gemChan<0) continue;
				
    		double DEDX_THR = THRESH;
    		int TimeWindowStart = 45;
    		int TimeMin = 0;
    		int TimeMax = 100;
		    
        //-- pedestal calculation for f125 raw data
        int nped = 0, ped = 100;
        double ped_sum = 0.;
        for (int si=TimeWindowStart-15; si<TimeWindowStart; si++) {
          int ped_samp = f125_wraw_samples->at(f125_wraw_samples_index->at(i)+si);
          ped_sum += ped_samp;
          nped++;
        }
        ped = ped_sum / nped;
        if (0. > ped || ped > 200 ) ped = 100;
        
	      for (int si=0; si<fadc_window; si++) {
      	  int time=si;
      	  int adc = f125_wraw_samples->at(f125_wraw_samples_index->at(i)+si);
      	  adc = adc-ped;
          //if (adc>4090) printf("!!!!!!!!!!!!!!!!!!!!!! ADC 125 overflow: %d \n",adc);
      	  if (adc>DEDX_THR) {
      	    double adc_fill=adc;
      	    //if (gemChan>-1) { f125_el_Yraw->Fill(time,gemChan,adc); }
      	    time-=TimeWindowStart;
      	    if ( TimeMin > time || time > TimeMax ) continue; // --- drop early and late hits from clustering ---
      	    hevtc->SetBinContent(time,gemChan,adc/100.);
      	    hevt->SetBinContent(time,gemChan,adc/100.);
      	    hentries++;
      	  }
      	} // --  end of timing samples loop
      } // -- end of fadc125 raw hit loop
      
      #ifdef SHOW_EVT_DISPLAY
        if (jentry<NPRT) {
      	  c2->cd(1);   hevt->Draw("colz");
      	  c2->cd(2);   hevtf->Draw("text");
      	  c2->Modified(); c2->Update();
        }
      #endif
      
      //=============================================================================
      //          Begin Sergey's NN Clustering & Track Fitting for GEM-TRD
      //=============================================================================
      
      if (hentries<100) continue; //-- skip event !!!!!
      #if (USE_CLUST>0)
        // ---------------- histogram hit distribution clustering -----------------
        float clust_Xpos[MAX_CLUST];
        float clust_Ypos[MAX_CLUST];
        float clust_Zpos[MAX_CLUST];
        float clust_dEdx[MAX_CLUST];
        float clust_Size[MAX_CLUST];
        float clust_Width[MAX_CLUST][3];  // y1, y2, dy ; strips
        float clust_Length[MAX_CLUST][3]; // x1, x2, dx ; time
        
        float hits_Xpos[500];
        float hits_Ypos[500];
        float hits_Zpos[500];
        float hits_dEdx[500];
        float hits_Size[MAX_CLUST];
        float hits_Width[MAX_CLUST];  // y1, y2, dy ; strips
        float hits_Length[MAX_CLUST]; // x1, x2, dx ; time
        
        for (int k=0; k<MAX_CLUST; k++) {
        	clust_Xpos[k]=0; clust_Ypos[k]=0; clust_Zpos[k]=0; clust_dEdx[k]=0;  clust_Size[k]=0;
        	clust_Width[k][0]=999999;   	clust_Width[k][1]=-999999;   	clust_Width[k][2]=0;
        	clust_Length[k][0]=999999;  	clust_Length[k][1]=-999999;  	clust_Length[k][2]=0;
        }
        float CL_DIST=2.7; // mm
        int nclust=0;
        
        TH2F* hp = hevt; // -- hevt and hevtc should be same bin size
        TH2F* hpc = hevtc;
        int nx=hp->GetNbinsX();    int ny=hp->GetNbinsY();
        double xmi=hp->GetXaxis()->GetBinLowEdge(1);     double xma=hp->GetXaxis()->GetBinUpEdge(nx);
        double ymi=hp->GetYaxis()->GetBinLowEdge(1);     double yma=hp->GetYaxis()->GetBinUpEdge(ny);
        double binx = (xma-xmi)/nx;      double biny = (yma-ymi)/ny;
        double THR2 = 1.2;
        
        // scan BinContent of hp = hevt
        for (int ix=0; ix<nx; ix++) {  //-------------------- clustering loop ------------------------------------
  	      for (int iy=0; iy<ny; iy++) {
      	    double c1 = hpc->GetBinContent(ix,iy);                    // energy
      	    double x1=double(ix)/double(nx)*(xma-xmi)+xmi-binx/2.;    // drift time
      	    double y1=double(iy)/double(ny)*(yma-ymi)+ymi-biny/2.;    // X strip
            
      	    if (c1<THR2) continue;

            // first iteration if c1 > threshold, set clust_Xpos and clust_Zpos
      	    if (nclust==0) {
      	      clust_Xpos[nclust]=y1; clust_Ypos[nclust]=0; clust_Zpos[nclust]=x1;  clust_dEdx[nclust]=c1;  clust_Size[nclust]=1;
      	      clust_Width[nclust][0]=y1;   	clust_Width[nclust][1]=y1;   	clust_Width[nclust][2]=0;
      	      clust_Length[nclust][0]=x1;  	clust_Length[nclust][1]=x1;  	clust_Length[nclust][2]=0;
      	      nclust++; continue;
      	    }
            
        	  int added=0;
        	  for (int k=0; k<nclust; k++) {
        	    double dist=sqrt(pow((y1-clust_Xpos[k]),2.)+pow((x1-clust_Zpos[k]),2.)); //--- dist hit to clusters
        	    // check the distance from the x1,y1 to the center of the cluster based on the radius (2.7 mm)
        	    if (dist<CL_DIST) {
                // if it's within the radius set this as a new position using weighted average to approximate the new central position
        	      clust_Xpos[k]=(y1*c1+clust_Xpos[k]*clust_dEdx[k])/(c1+clust_dEdx[k]);  //--  new X pos
        	      clust_Zpos[k]=(x1*c1+clust_Zpos[k]*clust_dEdx[k])/(c1+clust_dEdx[k]);  //--  new Z pos
                // new dedx is the sum of the two weighted averaged amplitude
        	      clust_dEdx[k]=c1+clust_dEdx[k];  // new dEdx
                // add cluster size
        	      clust_Size[k]=1+clust_Size[k];  // clust size in pixels

                // update cluster width in x and y
        	      if (y1<clust_Width[k][0]) clust_Width[k][0]=y1; if (y1>clust_Width[k][1]) clust_Width[k][1]=y1; clust_Width[k][2]=clust_Width[k][1]-clust_Width[k][0];
        	      if (x1<clust_Length[k][0]) clust_Length[k][0]=x1;if (x1>clust_Length[k][1]) clust_Length[k][1]=x1;clust_Length[k][2]=clust_Length[k][1]-clust_Length[k][0];
        	      hpc->SetBinContent(ix,iy,k+1.);
        	      added=1; break;
        	    }
        	  }

            // if it's outside the radius, set this as a new center
        	  if (added==0) {
        	    if (nclust+1>=MAX_CLUST) continue;
        	    clust_Xpos[nclust]=y1; clust_Ypos[nclust]=0; clust_Zpos[nclust]=x1;  clust_dEdx[nclust]=c1;  clust_Size[nclust]=1;
        	    clust_Width[nclust][0]=y1;   	clust_Width[nclust][1]=y1;   	clust_Width[nclust][2]=0;
        	    clust_Length[nclust][0]=x1;  	clust_Length[nclust][1]=x1;  	clust_Length[nclust][2]=0;
        	    nclust++;
        	  }
        	}
        } //----------- end  clustering loop ---------------
        
        int MinClustSize=10;
        double MinClustWidth=0.001;
        double MinClustLength=0.01;
        double zStart =  0.; // mm
        double zEnd   = 40.; // mm
        int ii=0;
        
        for (int k=0; k<nclust; k++) {
      	  //-------------  Cluster Filter -----------------
      	  if (clust_Size[k] >= MinClustSize && zStart < clust_Zpos[k] && clust_Zpos[k] < zEnd && clust_Width[k][2]>MinClustWidth ) {
      	    hits_Xpos[ii]=clust_Xpos[k];
      	    hits_Ypos[ii]=clust_Ypos[k];
      	    hits_Zpos[ii]=clust_Zpos[k];
      	    hits_dEdx[ii]=clust_dEdx[k];
      	    hits_Width[ii]=clust_Width[k][2];
            hits_Length[ii]=clust_Length[k][2];
      	    ii++;
      	  } else {
      	    //printf(" <--- skip \n");
      	  }
        }
        int nhits=ii;
        //------------------  end histogram hit distribution clustering  ------------------------
        
        //================ Draw HITS and CLUST  ======================
        #ifdef SHOW_EVTbyEVT
          char hevtTitle[80]; sprintf(hevtTitle,"GEM-TRD: Event=%lld Run=%d #Tracks=%d; z pos [mm]; y pos [mm]",jentry,RunNum,NTRACKS);
          hevt->SetTitle(hevtTitle);
          if (jentry<NPRT) {
        	  c2->cd(1); gPad->Modified(); gPad->Update();
        	  int COLMAP[]={1,2,3,4,6,5};
        	  int pmt=22 ,pmt0 = 20; // PM type
        	  for (int i = 0; i < nclust; i++) {
        	    TMarker m = TMarker(clust_Zpos[i],clust_Xpos[i],pmt);
        	    int tcol=2; //min(tracks[i],6);
        	    if (clust_Size[i]<MinClustSize) pmt=22; else pmt=pmt0;
        	    int mcol = COLMAP[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(pmt);
        	    m.SetMarkerSize(0.7+clust_dEdx[i]/300);
        	    m.DrawClone();
        	  }
        	  gPad->Modified(); gPad->Update();
          }
        #endif
        
        #if (USE_GNN==1)   // GNN MC
          //----------------------------------------------------------
          //--   Send to Model simulation
          //----------------------------------------------------------
          
          //printf("**> Start Model simulation nclust=%d nhits=%d \n",nclust,nhits);
          std::vector<int> tracks(nhits, 0);
          std::vector<float> Xcl;
          std::vector<float> Zcl;
          Xcl.clear();
          Zcl.clear();
          for (int n=0; n<nhits; n++) {
        	  Xcl.push_back(hits_Xpos[n]);
        	  Zcl.push_back(hits_Zpos[n]);
          }
          doPattern(Xcl, Zcl, tracks);  //---- call GNN ---
          //printf("**> End Model simulation \n"); //===================================================
          #ifdef SHOW_EVTbyEVT
            if (jentry<NPRT) {
              c2->cd(2); gPad->Modified(); gPad->Update();
              int COLMAP[]={1,2,3,4,6,5};
              for(int i = 0; i < (int)tracks.size(); i++){
            	  TMarker m = TMarker(hits_Zpos[i],hits_Xpos[i],24);
            	  int tcol=min(tracks[i],6);
            	  int mcol = COLMAP[tcol];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(41);     m.SetMarkerSize(1.5);
            	  m.DrawClone();  gPad->Modified(); gPad->Update();
              }
              //////////////////////////////////////////gPad->Modified(); gPad->Update();
              //printf("**> End Cluster Plot \n");
            }
          #endif
          //--------------------------------------------------
          //----           Track fitting                 -----
          //--------------------------------------------------
          
          std::vector<std::vector<float>> TRACKS;  // -- [Nnodes][N_params]
          TRACKS.resize(nhits);
          std::vector<float> hit_coord(2,0);
          std::vector<int>  TRACKS_N(nhits, 0); // [Nnodes]
          for (int i = 0; i < nhits; i++)  { TRACKS_N[i] = 0;  }
          std::vector<float> xz(2,0);
          
          for (int i2 = 0; i2 < nhits; i2++) {
        	  int num =  tracks[i2];
        	  int num2 = std::max(0, std::min(num, nhits - 1));
        	  xz[0]=Xcl[i2]; 	xz[1]=Zcl[i2];
        	  TRACKS[num2].push_back(Xcl[i2]); TRACKS[num2].push_back(Zcl[i2]);
        	  TRACKS_N[num2]++;
          }
          
          #if (DEBUG > 1)
            for (int i2 = 0; i2 < nhits; i2++) {
    	        printf(" trdID=%d n_hits=%d v_size=%d \n",i2,TRACKS_N[i2],TRACKS[i2].size());
    	        for (int i3 = 0; i3 < TRACKS[i2].size(); i3+=2) {
    	          printf(" i2=%d  i3=%d nhits-on-track=%d sz2=%d\n",i2,i3,TRACKS_N[i2], TRACKS[i2].size());
    	          printf(" trkID=%d  hit=%d x=%f z=%f \n",i2,i3/2,TRACKS[i2].at(i3),TRACKS[i2].at(i3+1));
    	        }
    	        if ( TRACKS_N[i2]>0) printf("\n");
            }
          #endif
          //--- end tracks sorting ---
          
          #if (USE_FIT==1)
            
            //-----------------------------------
            //---       linear fitting          ---
            //-----------------------------------
            
            static TMultiGraph *mg;
            if (mg != NULL ) delete mg;
            mg = new TMultiGraph();
            mg->SetTitle(" ML-FPGA response; z pos,mm; y pos,mm ");
            NTRACKS=0;
            Double_t p0, p1;
            // look for tracks
            for (int i2 = 1; i2 < nhits; i2++) {  // tracks loop; zero track -> noise
              
      	      if (TRACKS_N[i2]<2) continue;   //---- select 2 (x,z) and more hits on track ----
          	  std::vector<Double_t> x;
          	  std::vector<Double_t> y;
          	  for (int i3 = 0; i3 < (int)TRACKS[i2].size(); i3+=2) {
          	    x.push_back(TRACKS[i2].at(i3+1));
          	    y.push_back(TRACKS[i2].at(i3));
          	  }
          	  #ifdef SHOW_EVTbyEVT
            	  TGraph *g = new TGraph(TRACKS_N[i2], &x[0], &y[0]);  g->SetMarkerStyle(21); g->SetMarkerColor(i2);
            	  TF1 *f = new TF1("f", "[1] * x + [0]");
            	  g->Fit(f,"q0");
            	  //  --- get fit parameters ---
            	  TF1 *ffunc=g->GetFunction("f");
            	  p0=ffunc->GetParameter(0);
            	  p1=ffunc->GetParameter(1);
                Double_t chi2x_nn = ffunc->GetChisquare();
                Double_t Ndfx_nn = ffunc->GetNDF();
                chi2nn_gem=chi2x_nn/Ndfx_nn;

                if (jentry<NPRT) {
    	            mg->Add(g,"p");
                }
              #endif
  	          NTRACKS++;
            }  //  end tracks loop
            
            if(NTRACKS==1) {
              Count("sglTRK");
              a0nn=p0;
              a1nn=p1;
            }
            if (NTRACKS>1) Count("multTRK");
            #ifdef SHOW_EVTbyEVT
              if (NTRACKS<1) continue; //-- # of tracks seen inside GEM-TRD
              if (gem_trk_hit<1) continue; //-- # of tracks correlated between the GEM-TRD & the GEM Trackers
              //if (jentry<NPRT) {
                c2->cd(3); mg->Draw("AP");
                mg->GetXaxis()->SetLimits(Xmin,Xmax);
                mg->SetMinimum(Ymin);
                mg->SetMaximum(Ymax);
                gPad->Modified(); gPad->Update();
              //}
            #endif
          #endif // USE_FIT
        #endif // USE_GNN MC
        
        //******************************************************************************
        #ifdef SHOW_EVTbyEVT
          //if (NTRACKS<1) continue;
          cout<<"Event#="<<event_num<<" #ofTracks="<<NTRACKS<<endl;
          #ifdef WRITE_CSV
            WriteToCSV(csvFile,event_num,NTRACKS,a0,a1,chi2cc_gem);
          #endif
          c2->cd(4);  f125_fit->Draw("box");    fx1.Draw("same");   gPad->Modified(); gPad->Update();
          c2->cd(5);  hevt->Draw("colz");       gPad->Modified(); gPad->Update();
          printf(" all done, click middle pad ... a0=%f a1=%f (%f deg)  fx1(150)=%f chi2cc_gem=%f  \n",a0,a1,a1/3.1415*180.,fx1.Eval(150.),chi2cc_gem);
          c2->cd(5); gPad->WaitPrimitive();
        #endif
        
      #endif  // (USE_CLUST)
      //=======================  End Fa125 RAW Data Rrocessing Loop  =======================================
    #endif   //  USE 125 RAW
    
    //-- Fill GEM-TRD TTree of Hit Info
    #ifdef SAVE_TRACK_HITS
      if (gem_nhit>0) EVENT_VECT_GEM->Fill();
      if (trkr_hit>0) EVENT_VECT_GEM_TRACKERS->Fill();
    #endif
  } //===================== End of Event Loop ===============================
  timer.Stop();
  cout<<"***>>>  End Event Loop, Elapsed Time:"<<endl; timer.Print();
  
  #ifdef WRITE_CSV
    csvFile.close();
  #endif
  
  //=====================================================================================
  //===                 S A V E   G E M - T R D   H I T   T T R E E                  ====
  //=====================================================================================
  
  #ifdef SAVE_TRACK_HITS
    printf("Writing TTree Hit Info File... \n");
    fHits->cd();
    EVENT_VECT_GEM->Write();
    EVENT_VECT_GEM_TRACKERS->Write();
    fHits->Close();
    printf("TTree file written & closed OK \n");
  #endif
  
  //=====================================================================================
  //===                 P L O T    H I S T O G R A M S                                ===
  //=====================================================================================
  const char *OutputDir="RootOutput/halld24";
  #ifdef SAVE_PDF
    char ctit[120];
    sprintf(G_DIR,"%s/Run_%06d",OutputDir,RunNum);
    sprintf(ctit,"File=%s",G_DIR);
    bool COMPACT=false;
    TCanvas *cc;
    int nxd=3;
    int nyd=5;
    char pdfname[120];  sprintf(pdfname,"%s_evdisp.pdf",G_DIR);  //c0->Print(pdfname);
    
    //---------------------  page 1 --------------------
    htitle(" Count ");
    cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hcount->Draw();
    
    //---------------------  page 2 --------------------
    htitle(" (SRS) GEM-TRACKERS  ");   if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=4;
    cc=NextPlot(nxd,nyd); hgemtrkr1_peak_x->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr2_peak_x->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr1_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr2_peak_y->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr1_max_x->Draw();
    cc=NextPlot(nxd,nyd); hgemtrkr2_max_x->Draw();
    
    htitle(" (SRS) GEM-TRACKERS  ");   if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=4;
    cc=NextPlot(nxd,nyd); hgemtrkr1_peak_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr2_peak_xy->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_double_y->Draw("colz");
    cc=NextPlot(nxd,nyd); hgemtrkr_double_x->Draw("colz");
    
    //---------------------  page 3 --------------------
    htitle(" GEMTRD (fadc125) Amp ");    if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=4;
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  f125_el_x->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  f125_el_Xmax->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  f125_el_y->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  f125_el_Ymax->Draw();
    
    //---------------------  page 3 --------------------
    htitle(" GEMTRD (fadc125) Amp ");    if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=4;
    cc=NextPlot(nxd,nyd);   f125_el_xch_max->Draw();
    cc=NextPlot(nxd,nyd);   f125_el_ych_max->Draw();
    cc=NextPlot(nxd,nyd);   f125_xamp2d->Draw("colz");
    cc=NextPlot(nxd,nyd);   f125_yamp2d->Draw("colz");
    
    //---------------------  page 3 --------------------
    htitle(" GEMTRD (fadc125) Amp (Track Matching Condition)");    if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=4;
    cc=NextPlot(nxd,nyd);   f125_ydiff2d->Draw("colz");
    cc=NextPlot(nxd,nyd);   f125_yamp->Draw("colz");
    cc=NextPlot(nxd,nyd);   teff->Draw();
    cc=NextPlot(nxd,nyd);   teff2d->Draw();
    cc=NextPlot(nxd,nyd);   thits->Draw();
    
    //--- close PDF file ----
    cc=NextPlot(-1,-1);
  #endif
  
  //=====================================================================================
  //===                 S A V E   H I S T O G R A M S                                ====
  //=====================================================================================
  
  fOut->cd();
  cout<<"Writing Output File: "<<rootFileName<<endl;
  fOut->Write();
  fOut->Close();
  delete fOut;
  
  cout<<"================= END OF RUN "<<RunNum<<" ================"<<endl;
}
//=========== The End ===================
