//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 29 20:15:03 2023 by ROOT version 6.28/00
// from TTree events/jana4ml4fpga_tree_v1
// found on file: ROOT/Run_003264.root
//////////////////////////////////////////////////////////

#ifndef trdclass_halld24_h
#define trdclass_halld24_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <sstream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include <iostream>
#include "TRandom3.h"
#include "TCanvas.h"
#include "TMultiLayerPerceptron.h"
#include "TMLPAnalyzer.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TString.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLinearFitter.h"
#include "TGraphErrors.h"
#include <TError.h>
#include "TF1.h"
#include "TPaveStats.h"
#include "TCutG.h"
#include "TProfile.h"
#include "TBox.h"
#include "TSystem.h"
#include "TSystem.h"
#include "stdio.h"


// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class trdclass_halld24 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   ULong64_t       srs_raw_count;
   vector<unsigned int> *srs_raw_roc;
   vector<unsigned int> *srs_raw_slot;
   vector<unsigned int> *srs_raw_channel;
   vector<unsigned int> *srs_raw_apv_id;
   vector<unsigned int> *srs_raw_channel_apv;
   vector<unsigned short> *srs_raw_best_sample;
   vector<unsigned short> *srs_raw_samples_index;
   vector<unsigned short> *srs_raw_samples_count;
   vector<unsigned short> *srs_raw_samples;
   ULong64_t       f125_wraw_count;
   vector<unsigned int> *f125_wraw_roc;
   vector<unsigned int> *f125_wraw_slot;
   vector<unsigned int> *f125_wraw_channel;
   vector<bool>    *f125_wraw_invalid_samples;
   vector<bool>    *f125_wraw_overflow;
   vector<unsigned int> *f125_wraw_itrigger;
   vector<unsigned short> *f125_wraw_samples_index;
   vector<unsigned short> *f125_wraw_samples_count;
   vector<unsigned short> *f125_wraw_samples;
   ULong64_t       f250_wraw_count;
   vector<unsigned int> *f250_wraw_roc;
   vector<unsigned int> *f250_wraw_slot;
   vector<unsigned int> *f250_wraw_channel;
   vector<bool>    *f250_wraw_invalid_samples;
   vector<bool>    *f250_wraw_overflow;
   vector<unsigned int> *f250_wraw_itrigger;
   vector<unsigned short> *f250_wraw_samples_index;
   vector<unsigned short> *f250_wraw_samples_count;
   vector<unsigned short> *f250_wraw_samples;
   ULong64_t       f125_pulse_count;
   vector<unsigned int> *f125_pulse_roc;
   vector<unsigned int> *f125_pulse_slot;
   vector<unsigned int> *f125_pulse_channel;
   vector<unsigned int> *f125_pulse_npk;
   vector<unsigned int> *f125_pulse_le_time;
   vector<unsigned int> *f125_pulse_time_quality_bit;
   vector<unsigned int> *f125_pulse_overflow_count;
   vector<unsigned int> *f125_pulse_pedestal;
   vector<unsigned int> *f125_pulse_integral;
   vector<unsigned int> *f125_pulse_peak_amp;
   vector<unsigned int> *f125_pulse_peak_time;
   vector<unsigned int> *f125_pulse_word1;
   vector<unsigned int> *f125_pulse_word2;
   vector<unsigned int> *f125_pulse_nsamples_pedestal;
   vector<unsigned int> *f125_pulse_nsamples_integral;
   vector<bool>    *f125_pulse_emulated;
   vector<unsigned int> *f125_pulse_le_time_emulated;
   vector<unsigned int> *f125_pulse_time_quality_bit_emulated;
   vector<unsigned int> *f125_pulse_overflow_count_emulated;
   vector<unsigned int> *f125_pulse_pedestal_emulated;
   vector<unsigned int> *f125_pulse_integral_emulated;
   vector<unsigned int> *f125_pulse_peak_amp_emulated;
   vector<unsigned int> *f125_pulse_peak_time_emulated;
   ULong64_t       f250_pulse_count;
   vector<unsigned int> *f250_pulse_roc;
   vector<unsigned int> *f250_pulse_slot;
   vector<unsigned int> *f250_pulse_channel;
   vector<unsigned int> *f250_pulse_event_within_block;
   vector<bool>    *f250_pulse_qf_pedestal;
   vector<unsigned int> *f250_pulse_pedestal;
   vector<unsigned int> *f250_pulse_integral;
   vector<bool>    *f250_pulse_qf_nsa_beyond_ptw;
   vector<bool>    *f250_pulse_qf_overflow;
   vector<bool>    *f250_pulse_qf_underflow;
   vector<unsigned int> *f250_pulse_nsamples_over_threshold;
   vector<unsigned int> *f250_pulse_course_time;
   vector<unsigned int> *f250_pulse_fine_time;
   vector<unsigned int> *f250_pulse_pulse_peak;
   vector<bool>    *f250_pulse_qf_vpeak_beyond_nsa;
   vector<bool>    *f250_pulse_qf_vpeak_not_found;
   vector<bool>    *f250_pulse_qf_bad_pedestal;
   vector<unsigned int> *f250_pulse_pulse_number;
   vector<unsigned int> *f250_pulse_nsamples_integral;
   vector<unsigned int> *f250_pulse_nsamples_pedestal;
   vector<bool>    *f250_pulse_emulated;
   vector<unsigned int> *f250_pulse_integral_emulated;
   vector<unsigned int> *f250_pulse_pedestal_emulated;
   vector<unsigned int> *f250_pulse_time_emulated;
   vector<unsigned int> *f250_pulse_course_time_emulated;
   vector<unsigned int> *f250_pulse_fine_time_emulated;
   vector<unsigned int> *f250_pulse_pulse_peak_emulated;
   vector<unsigned int> *f250_pulse_qf_emulated;
   ULong64_t       gem_scluster_count;
   vector<double>  *gem_scluster_x;
   vector<double>  *gem_scluster_y;
   vector<double>  *gem_scluster_energy;
   vector<double>  *gem_scluster_adc;
   ULong64_t       srs_prerecon_count;
   vector<double>  *srs_prerecon_y;
   vector<double>  *srs_prerecon_x;
   ULong64_t       gem_peak_count;
   vector<unsigned int> *gem_peak_plane_id;
   vector<string>  *gem_peak_plane_name;
   vector<unsigned int> *gem_peak_index;
   vector<unsigned int> *gem_peak_apv_id;
   vector<double>  *gem_peak_height;
   vector<double>  *gem_peak_width;
   vector<double>  *gem_peak_area;
   vector<double>  *gem_peak_real_pos;

   // List of branches
   //============== SRS Raw Data Variables ===============
   TBranch        *b_srs_raw_count;   //!
   TBranch        *b_srs_raw_roc;   //!
   TBranch        *b_srs_raw_slot;   //!
   TBranch        *b_srs_raw_channel;   //!
   TBranch        *b_srs_raw_apv_id;   //!
   TBranch        *b_srs_raw_channel_apv;   //!
   TBranch        *b_srs_raw_best_sample;   //!
   TBranch        *b_srs_raw_samples_index;   //!
   TBranch        *b_srs_raw_samples_count;   //!
   TBranch        *b_srs_raw_samples;   //!
   //============== fADC125 Raw Data Variables ===============
   TBranch        *b_f125_wraw_count;   //!
   TBranch        *b_f125_wraw_roc;   //!
   TBranch        *b_f125_wraw_slot;   //!
   TBranch        *b_f125_wraw_channel;   //!
   TBranch        *b_f125_wraw_invalid_samples;   //!
   TBranch        *b_f125_wraw_overflow;   //!
   TBranch        *b_f125_wraw_itrigger;   //!
   TBranch        *b_f125_wraw_samples_index;   //!
   TBranch        *b_f125_wraw_samples_count;   //!
   TBranch        *b_f125_wraw_samples;   //!
   //============== fADC250 Raw Data Variables ===============
   TBranch        *b_f250_wraw_count;   //!
   TBranch        *b_f250_wraw_roc;   //!
   TBranch        *b_f250_wraw_slot;   //!
   TBranch        *b_f250_wraw_channel;   //!
   TBranch        *b_f250_wraw_invalid_samples;   //!
   TBranch        *b_f250_wraw_overflow;   //!
   TBranch        *b_f250_wraw_itrigger;   //!
   TBranch        *b_f250_wraw_samples_index;   //!
   TBranch        *b_f250_wraw_samples_count;   //!
   TBranch        *b_f250_wraw_samples;   //!
   //============== fADC125 Pulse Data Variables ===============
   TBranch        *b_f125_pulse_count;   //!
   TBranch        *b_f125_pulse_roc;   //!
   TBranch        *b_f125_pulse_slot;   //!
   TBranch        *b_f125_pulse_channel;   //!
   TBranch        *b_f125_pulse_npk;   //!
   TBranch        *b_f125_pulse_le_time;   //!
   TBranch        *b_f125_pulse_time_quality_bit;   //!
   TBranch        *b_f125_pulse_overflow_count;   //!
   TBranch        *b_f125_pulse_pedestal;   //!
   TBranch        *b_f125_pulse_integral;   //!
   TBranch        *b_f125_pulse_peak_amp;   //!
   TBranch        *b_f125_pulse_peak_time;   //!
   TBranch        *b_f125_pulse_word1;   //!
   TBranch        *b_f125_pulse_word2;   //!
   TBranch        *b_f125_pulse_nsamples_pedestal;   //!
   TBranch        *b_f125_pulse_nsamples_integral;   //!
   TBranch        *b_f125_pulse_emulated;   //!
   TBranch        *b_f125_pulse_le_time_emulated;   //!
   TBranch        *b_f125_pulse_time_quality_bit_emulated;   //!
   TBranch        *b_f125_pulse_overflow_count_emulated;   //!
   TBranch        *b_f125_pulse_pedestal_emulated;   //!
   TBranch        *b_f125_pulse_integral_emulated;   //!
   TBranch        *b_f125_pulse_peak_amp_emulated;   //!
   TBranch        *b_f125_pulse_peak_time_emulated;   //!
   //============== fADC250 Pulse Data Variables ===============
   TBranch        *b_f250_pulse_count;   //!
   TBranch        *b_f250_pulse_roc;   //!
   TBranch        *b_f250_pulse_slot;   //!
   TBranch        *b_f250_pulse_channel;   //!
   TBranch        *b_f250_pulse_event_within_block;   //!
   TBranch        *b_f250_pulse_qf_pedestal;   //!
   TBranch        *b_f250_pulse_pedestal;   //!
   TBranch        *b_f250_pulse_integral;   //!
   TBranch        *b_f250_pulse_qf_nsa_beyond_ptw;   //!
   TBranch        *b_f250_pulse_qf_overflow;   //!
   TBranch        *b_f250_pulse_qf_underflow;   //!
   TBranch        *b_f250_pulse_nsamples_over_threshold;   //!
   TBranch        *b_f250_pulse_course_time;   //!
   TBranch        *b_f250_pulse_fine_time;   //!
   TBranch        *b_f250_pulse_pulse_peak;   //!
   TBranch        *b_f250_pulse_qf_vpeak_beyond_nsa;   //!
   TBranch        *b_f250_pulse_qf_vpeak_not_found;   //!
   TBranch        *b_f250_pulse_qf_bad_pedestal;   //!
   TBranch        *b_f250_pulse_pulse_number;   //!
   TBranch        *b_f250_pulse_nsamples_integral;   //!
   TBranch        *b_f250_pulse_nsamples_pedestal;   //!
   TBranch        *b_f250_pulse_emulated;   //!
   TBranch        *b_f250_pulse_integral_emulated;   //!
   TBranch        *b_f250_pulse_pedestal_emulated;   //!
   TBranch        *b_f250_pulse_time_emulated;   //!
   TBranch        *b_f250_pulse_course_time_emulated;   //!
   TBranch        *b_f250_pulse_fine_time_emulated;   //!
   TBranch        *b_f250_pulse_pulse_peak_emulated;   //!
   TBranch        *b_f250_pulse_qf_emulated;   //!
   //============== SRS Cluster (Pulse?) Data Variables ===============
   TBranch        *b_gem_scluster_count;   //!
   TBranch        *b_gem_scluster_x;   //!
   TBranch        *b_gem_scluster_y;   //!
   TBranch        *b_gem_scluster_energy;   //!
   TBranch        *b_gem_scluster_adc;   //!
   TBranch        *b_srs_prerecon_count;   //!
   TBranch        *b_srs_prerecon_y;   //!
   TBranch        *b_srs_prerecon_x;   //!
   TBranch        *b_gem_peak_count;   //!
   TBranch        *b_gem_peak_plane_id;   //!
   TBranch        *b_gem_peak_plane_name;   //!
   TBranch        *b_gem_peak_index;   //!
   TBranch        *b_gem_peak_apv_id;   //!
   TBranch        *b_gem_peak_height;   //!
   TBranch        *b_gem_peak_width;   //!
   TBranch        *b_gem_peak_area;   //!
   TBranch        *b_gem_peak_real_pos;   //!

   trdclass_halld24(int RunNum, int MaxEvt, int FirstEvt, int FileNum);
   virtual ~trdclass_halld24();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   std::pair<Double_t, Double_t>   TrkFit(TH2F *h2, TF1 &fx, const char *cfx, int rob);
   void Count(const char *tit);
   void Count(const char *tit, double cut1);
   void Count(const char *tit, double cut1, double cut2);
   
   //==================  histograms ========================
  
   int RunNum;
   int FileNum;
   Long64_t MaxEvt;
   Long64_t FirstEvt;
   TH1F *h250_size;
   TH1D *hcount;
   TH1D *hNTracks, *hNTracks_e;
    TH1F *gem_trk_fit_integral;

   TH1F *hgemtrkr1_x;
   TH1F *hgemtrkr1_y;
   TH1F *hgemtrkr1_max_x;
   TH1F *hgemtrkr2_max_x;
   TH1F *hgemtrkr1_peak_x, *hgemtrkr2_peak_x, *hgemtrkr_peak_delta_x;
   TH1F *hgemtrkr1_peak_y, *hgemtrkr2_peak_y;
   TH2F *hgemtrkr_double_x, *hgemtrkr_double_y, *hgemtrkr1_peak_xy, *hgemtrkr2_peak_xy, *hgemtrkr1_ch_xy, *hgemtrkr1_peak_xy_chi2, *srs_gemtrd_el, *srs_etrd_beam, *srs_etrd_ratio;
   TH1F *f125_el_x, *f125_el_y, *f125_el_Xmax, *f125_el_Ymax, *f125_el_xch_max, *f125_el_ych_max, *f125_el_chi2, *f125_el_fita;
   TH2F *f125_xamp2d, *f125_yamp2d, *f125_el_Yraw, *f125_el_evt_display, *f125_el_Xraw, *f125_el_fit, *f125_fit, *f125_el_clu2d;
   
   TH2F *hevt, *hevtc, *hevti, *hevtf, *hevtL;
   
   //----- EVENT STRUCTURE -----
   TTree *EVENT_VECT_GEM;
   TTree *EVENT_VECT_GEM_TRACKERS;

   //---------------------------
   int event_num;
   int gem_nhit;
   int trkr_hit;
   int clu_nhit;
   int gem_trk_hit;

   std::vector <int> gem_xpos;
   std::vector <int> gem_ypos;
   std::vector <float> gem_zpos;
   std::vector <float> gem_dedx;
   //std::vector <int> gem_trackID;
   //std::vector <bool> gem_parID;
   std::vector <float> gem_zHist_vect;
   std::vector <float> clu_xpos;
   std::vector <float> clu_ypos;
   std::vector <float> clu_zpos;
   std::vector <float> clu_dedx;
   std::vector <float> clu_width;
   TH1F *gem_zHist;
   std::vector <int> gem_xch_max;
   std::vector <int> gem_ych_max;
   std::vector <float> gem_xamp_max;
   std::vector <float> gem_yamp_max;
   std::vector <float> gem_xtime_max;
   std::vector <float> gem_ytime_max;
   std::vector <float> chi2cc_gem;
   std::vector <float> chi2nn_gem;
   std::vector <float> gem_integral;
   std::vector <float> a0;
   std::vector <float> a1;
   std::vector <float> a0nn;
   std::vector <float> a1nn;
   std::vector <int> gemtrkr_xch_max;
   std::vector <int> gemtrkr_ych_max;
   std::vector <int> gemtrkr_2_xch_max;
   std::vector <int> gemtrkr_2_ych_max;
   std::vector <float> gemtrkr_xamp_max;
   std::vector <float> gemtrkr_yamp_max;
   std::vector <float> gemtrkr_2_xamp_max;
   std::vector <float> gemtrkr_2_yamp_max;
   std::vector <int> nxtrkr1; 
   std::vector <int> nytrkr1;
   std::vector <int> nxtrkr2;
   std::vector <int> nytrkr2;
   
   //=============================================
};

#endif

#ifdef trdclass_halld24_cxx
trdclass_halld24::trdclass_halld24(int RunNum_in, int MaxEvt_in=0, int FirstEvt_in=0, int FileNum_in=-1) : fChain(0)
{
  RunNum=RunNum_in;
  MaxEvt=MaxEvt_in;
  FirstEvt=FirstEvt_in;
  FileNum=FileNum_in;
  TTree *tree=NULL;
  
  printf("====== trdclass constructor Run=%d File=%d  ============\n",RunNum,FileNum);

  TChain* chain;
  char chainFiles[128];
  
  if (FileNum<0.) {
    sprintf(chainFiles,"ROOT/Run_%06d_*.root",RunNum);
    chain = new TChain("events");
    printf("Chaining data files into TTree chain for Run %d \n",RunNum);
    chain->Add(chainFiles);
   //  chain->Print();
    tree=chain;
  }
  
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if (tree == 0) {
     char FileName[128];
     sprintf(FileName,"ROOT/Run_%06d_%03d.root",RunNum,FileNum);
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(FileName);
      if (!f || !f->IsOpen()) {
         f = new TFile(FileName);
       printf("---> Open file = %s \n",FileName);
     } else {
       sprintf(FileName,"ROOT/Run_%06d.root",RunNum);
       f = (TFile*)gROOT->GetListOfFiles()->FindObject(FileName);
       printf("---> Open file = %s \n",FileName);
      }   
      f->GetObject("events",tree);
   }
   Init(tree);
}

trdclass_halld24::~trdclass_halld24()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t trdclass_halld24::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t trdclass_halld24::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void trdclass_halld24::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   srs_raw_roc = 0;
   srs_raw_slot = 0;
   srs_raw_channel = 0;
   srs_raw_apv_id = 0;
   srs_raw_channel_apv = 0;
   srs_raw_best_sample = 0;
   srs_raw_samples_index = 0;
   srs_raw_samples_count = 0;
   srs_raw_samples = 0;
   f125_wraw_roc = 0;
   f125_wraw_slot = 0;
   f125_wraw_channel = 0;
   f125_wraw_invalid_samples = 0;
   f125_wraw_overflow = 0;
   f125_wraw_itrigger = 0;
   f125_wraw_samples_index = 0;
   f125_wraw_samples_count = 0;
   f125_wraw_samples = 0;
   f250_wraw_roc = 0;
   f250_wraw_slot = 0;
   f250_wraw_channel = 0;
   f250_wraw_invalid_samples = 0;
   f250_wraw_overflow = 0;
   f250_wraw_itrigger = 0;
   f250_wraw_samples_index = 0;
   f250_wraw_samples_count = 0;
   f250_wraw_samples = 0;
   f125_pulse_roc = 0;
   f125_pulse_slot = 0;
   f125_pulse_channel = 0;
   f125_pulse_npk = 0;
   f125_pulse_le_time = 0;
   f125_pulse_time_quality_bit = 0;
   f125_pulse_overflow_count = 0;
   f125_pulse_pedestal = 0;
   f125_pulse_integral = 0;
   f125_pulse_peak_amp = 0;
   f125_pulse_peak_time = 0;
   f125_pulse_word1 = 0;
   f125_pulse_word2 = 0;
   f125_pulse_nsamples_pedestal = 0;
   f125_pulse_nsamples_integral = 0;
   f125_pulse_emulated = 0;
   f125_pulse_le_time_emulated = 0;
   f125_pulse_time_quality_bit_emulated = 0;
   f125_pulse_overflow_count_emulated = 0;
   f125_pulse_pedestal_emulated = 0;
   f125_pulse_integral_emulated = 0;
   f125_pulse_peak_amp_emulated = 0;
   f125_pulse_peak_time_emulated = 0;
   f250_pulse_roc = 0;
   f250_pulse_slot = 0;
   f250_pulse_channel = 0;
   f250_pulse_event_within_block = 0;
   f250_pulse_qf_pedestal = 0;
   f250_pulse_pedestal = 0;
   f250_pulse_integral = 0;
   f250_pulse_qf_nsa_beyond_ptw = 0;
   f250_pulse_qf_overflow = 0;
   f250_pulse_qf_underflow = 0;
   f250_pulse_nsamples_over_threshold = 0;
   f250_pulse_course_time = 0;
   f250_pulse_fine_time = 0;
   f250_pulse_pulse_peak = 0;
   f250_pulse_qf_vpeak_beyond_nsa = 0;
   f250_pulse_qf_vpeak_not_found = 0;
   f250_pulse_qf_bad_pedestal = 0;
   f250_pulse_pulse_number = 0;
   f250_pulse_nsamples_integral = 0;
   f250_pulse_nsamples_pedestal = 0;
   f250_pulse_emulated = 0;
   f250_pulse_integral_emulated = 0;
   f250_pulse_pedestal_emulated = 0;
   f250_pulse_time_emulated = 0;
   f250_pulse_course_time_emulated = 0;
   f250_pulse_fine_time_emulated = 0;
   f250_pulse_pulse_peak_emulated = 0;
   f250_pulse_qf_emulated = 0;
   gem_scluster_x = 0;
   gem_scluster_y = 0;
   gem_scluster_energy = 0;
   gem_scluster_adc = 0;
   srs_prerecon_y = 0;
   srs_prerecon_x = 0;
   gem_peak_plane_id = 0;
   gem_peak_plane_name = 0;
   gem_peak_index = 0;
   gem_peak_apv_id = 0;
   gem_peak_height = 0;
   gem_peak_width = 0;
   gem_peak_area = 0;
   gem_peak_real_pos = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("srs_raw_count", &srs_raw_count, &b_srs_raw_count);
   fChain->SetBranchAddress("srs_raw_roc", &srs_raw_roc, &b_srs_raw_roc);
   fChain->SetBranchAddress("srs_raw_slot", &srs_raw_slot, &b_srs_raw_slot);
   fChain->SetBranchAddress("srs_raw_channel", &srs_raw_channel, &b_srs_raw_channel);
   fChain->SetBranchAddress("srs_raw_apv_id", &srs_raw_apv_id, &b_srs_raw_apv_id);
   fChain->SetBranchAddress("srs_raw_channel_apv", &srs_raw_channel_apv, &b_srs_raw_channel_apv);
   fChain->SetBranchAddress("srs_raw_best_sample", &srs_raw_best_sample, &b_srs_raw_best_sample);
   fChain->SetBranchAddress("srs_raw_samples_index", &srs_raw_samples_index, &b_srs_raw_samples_index);
   fChain->SetBranchAddress("srs_raw_samples_count", &srs_raw_samples_count, &b_srs_raw_samples_count);
   fChain->SetBranchAddress("srs_raw_samples", &srs_raw_samples, &b_srs_raw_samples);
   fChain->SetBranchAddress("f125_wraw_count", &f125_wraw_count, &b_f125_wraw_count);
   fChain->SetBranchAddress("f125_wraw_roc", &f125_wraw_roc, &b_f125_wraw_roc);
   fChain->SetBranchAddress("f125_wraw_slot", &f125_wraw_slot, &b_f125_wraw_slot);
   fChain->SetBranchAddress("f125_wraw_channel", &f125_wraw_channel, &b_f125_wraw_channel);
   fChain->SetBranchAddress("f125_wraw_invalid_samples", &f125_wraw_invalid_samples, &b_f125_wraw_invalid_samples);
   fChain->SetBranchAddress("f125_wraw_overflow", &f125_wraw_overflow, &b_f125_wraw_overflow);
   fChain->SetBranchAddress("f125_wraw_itrigger", &f125_wraw_itrigger, &b_f125_wraw_itrigger);
   fChain->SetBranchAddress("f125_wraw_samples_index", &f125_wraw_samples_index, &b_f125_wraw_samples_index);
   fChain->SetBranchAddress("f125_wraw_samples_count", &f125_wraw_samples_count, &b_f125_wraw_samples_count);
   fChain->SetBranchAddress("f125_wraw_samples", &f125_wraw_samples, &b_f125_wraw_samples);
   fChain->SetBranchAddress("f250_wraw_count", &f250_wraw_count, &b_f250_wraw_count);
   fChain->SetBranchAddress("f250_pulse_roc", &f250_pulse_roc, &b_f250_pulse_roc);
   fChain->SetBranchAddress("f250_pulse_slot", &f250_pulse_slot, &b_f250_pulse_slot);
   fChain->SetBranchAddress("f250_pulse_channel", &f250_pulse_channel, &b_f250_pulse_channel);
   fChain->SetBranchAddress("f250_wraw_roc", &f250_wraw_roc, &b_f250_wraw_roc);
   fChain->SetBranchAddress("f250_wraw_slot", &f250_wraw_slot, &b_f250_wraw_slot);
   fChain->SetBranchAddress("f250_wraw_channel", &f250_wraw_channel, &b_f250_wraw_channel);
   fChain->SetBranchAddress("f250_wraw_invalid_samples", &f250_wraw_invalid_samples, &b_f250_wraw_invalid_samples);
   fChain->SetBranchAddress("f250_wraw_overflow", &f250_wraw_overflow, &b_f250_wraw_overflow);
   fChain->SetBranchAddress("f250_wraw_itrigger", &f250_wraw_itrigger, &b_f250_wraw_itrigger);
   fChain->SetBranchAddress("f250_wraw_samples_index", &f250_wraw_samples_index, &b_f250_wraw_samples_index);
   fChain->SetBranchAddress("f250_wraw_samples_count", &f250_wraw_samples_count, &b_f250_wraw_samples_count);
   fChain->SetBranchAddress("f250_wraw_samples", &f250_wraw_samples, &b_f250_wraw_samples);
   fChain->SetBranchAddress("f125_pulse_count", &f125_pulse_count, &b_f125_pulse_count);
   fChain->SetBranchAddress("f125_pulse_roc", &f125_pulse_roc, &b_f125_pulse_roc);
   fChain->SetBranchAddress("f125_pulse_slot", &f125_pulse_slot, &b_f125_pulse_slot);
   fChain->SetBranchAddress("f125_pulse_channel", &f125_pulse_channel, &b_f125_pulse_channel);
   fChain->SetBranchAddress("f125_pulse_npk", &f125_pulse_npk, &b_f125_pulse_npk);
   fChain->SetBranchAddress("f125_pulse_le_time", &f125_pulse_le_time, &b_f125_pulse_le_time);
   fChain->SetBranchAddress("f125_pulse_time_quality_bit", &f125_pulse_time_quality_bit, &b_f125_pulse_time_quality_bit);
   fChain->SetBranchAddress("f125_pulse_overflow_count", &f125_pulse_overflow_count, &b_f125_pulse_overflow_count);
   fChain->SetBranchAddress("f125_pulse_pedestal", &f125_pulse_pedestal, &b_f125_pulse_pedestal);
   fChain->SetBranchAddress("f125_pulse_integral", &f125_pulse_integral, &b_f125_pulse_integral);
   fChain->SetBranchAddress("f125_pulse_peak_amp", &f125_pulse_peak_amp, &b_f125_pulse_peak_amp);
   fChain->SetBranchAddress("f125_pulse_peak_time", &f125_pulse_peak_time, &b_f125_pulse_peak_time);
   fChain->SetBranchAddress("f125_pulse_word1", &f125_pulse_word1, &b_f125_pulse_word1);
   fChain->SetBranchAddress("f125_pulse_word2", &f125_pulse_word2, &b_f125_pulse_word2);
   fChain->SetBranchAddress("f125_pulse_nsamples_pedestal", &f125_pulse_nsamples_pedestal, &b_f125_pulse_nsamples_pedestal);
   fChain->SetBranchAddress("f125_pulse_nsamples_integral", &f125_pulse_nsamples_integral, &b_f125_pulse_nsamples_integral);
   fChain->SetBranchAddress("f125_pulse_emulated", &f125_pulse_emulated, &b_f125_pulse_emulated);
   fChain->SetBranchAddress("f125_pulse_le_time_emulated", &f125_pulse_le_time_emulated, &b_f125_pulse_le_time_emulated);
   fChain->SetBranchAddress("f125_pulse_time_quality_bit_emulated", &f125_pulse_time_quality_bit_emulated, &b_f125_pulse_time_quality_bit_emulated);
   fChain->SetBranchAddress("f125_pulse_overflow_count_emulated", &f125_pulse_overflow_count_emulated, &b_f125_pulse_overflow_count_emulated);
   fChain->SetBranchAddress("f125_pulse_pedestal_emulated", &f125_pulse_pedestal_emulated, &b_f125_pulse_pedestal_emulated);
   fChain->SetBranchAddress("f125_pulse_integral_emulated", &f125_pulse_integral_emulated, &b_f125_pulse_integral_emulated);
   fChain->SetBranchAddress("f125_pulse_peak_amp_emulated", &f125_pulse_peak_amp_emulated, &b_f125_pulse_peak_amp_emulated);
   fChain->SetBranchAddress("f125_pulse_peak_time_emulated", &f125_pulse_peak_time_emulated, &b_f125_pulse_peak_time_emulated);
   fChain->SetBranchAddress("f250_pulse_count", &f250_pulse_count, &b_f250_pulse_count);
   fChain->SetBranchAddress("f250_pulse_event_within_block", &f250_pulse_event_within_block, &b_f250_pulse_event_within_block);
   fChain->SetBranchAddress("f250_pulse_qf_pedestal", &f250_pulse_qf_pedestal, &b_f250_pulse_qf_pedestal);
   fChain->SetBranchAddress("f250_pulse_pedestal", &f250_pulse_pedestal, &b_f250_pulse_pedestal);
   fChain->SetBranchAddress("f250_pulse_integral", &f250_pulse_integral, &b_f250_pulse_integral);
   fChain->SetBranchAddress("f250_pulse_qf_nsa_beyond_ptw", &f250_pulse_qf_nsa_beyond_ptw, &b_f250_pulse_qf_nsa_beyond_ptw);
   fChain->SetBranchAddress("f250_pulse_qf_overflow", &f250_pulse_qf_overflow, &b_f250_pulse_qf_overflow);
   fChain->SetBranchAddress("f250_pulse_qf_underflow", &f250_pulse_qf_underflow, &b_f250_pulse_qf_underflow);
   fChain->SetBranchAddress("f250_pulse_nsamples_over_threshold", &f250_pulse_nsamples_over_threshold, &b_f250_pulse_nsamples_over_threshold);
   fChain->SetBranchAddress("f250_pulse_course_time", &f250_pulse_course_time, &b_f250_pulse_course_time);
   fChain->SetBranchAddress("f250_pulse_fine_time", &f250_pulse_fine_time, &b_f250_pulse_fine_time);
   fChain->SetBranchAddress("f250_pulse_pulse_peak", &f250_pulse_pulse_peak, &b_f250_pulse_pulse_peak);
   fChain->SetBranchAddress("f250_pulse_qf_vpeak_beyond_nsa", &f250_pulse_qf_vpeak_beyond_nsa, &b_f250_pulse_qf_vpeak_beyond_nsa);
   fChain->SetBranchAddress("f250_pulse_qf_vpeak_not_found", &f250_pulse_qf_vpeak_not_found, &b_f250_pulse_qf_vpeak_not_found);
   fChain->SetBranchAddress("f250_pulse_qf_bad_pedestal", &f250_pulse_qf_bad_pedestal, &b_f250_pulse_qf_bad_pedestal);
   fChain->SetBranchAddress("f250_pulse_pulse_number", &f250_pulse_pulse_number, &b_f250_pulse_pulse_number);
   fChain->SetBranchAddress("f250_pulse_nsamples_integral", &f250_pulse_nsamples_integral, &b_f250_pulse_nsamples_integral);
   fChain->SetBranchAddress("f250_pulse_nsamples_pedestal", &f250_pulse_nsamples_pedestal, &b_f250_pulse_nsamples_pedestal);
   fChain->SetBranchAddress("f250_pulse_emulated", &f250_pulse_emulated, &b_f250_pulse_emulated);
   fChain->SetBranchAddress("f250_pulse_integral_emulated", &f250_pulse_integral_emulated, &b_f250_pulse_integral_emulated);
   fChain->SetBranchAddress("f250_pulse_pedestal_emulated", &f250_pulse_pedestal_emulated, &b_f250_pulse_pedestal_emulated);
   fChain->SetBranchAddress("f250_pulse_time_emulated", &f250_pulse_time_emulated, &b_f250_pulse_time_emulated);
   fChain->SetBranchAddress("f250_pulse_course_time_emulated", &f250_pulse_course_time_emulated, &b_f250_pulse_course_time_emulated);
   fChain->SetBranchAddress("f250_pulse_fine_time_emulated", &f250_pulse_fine_time_emulated, &b_f250_pulse_fine_time_emulated);
   fChain->SetBranchAddress("f250_pulse_pulse_peak_emulated", &f250_pulse_pulse_peak_emulated, &b_f250_pulse_pulse_peak_emulated);
   fChain->SetBranchAddress("f250_pulse_qf_emulated", &f250_pulse_qf_emulated, &b_f250_pulse_qf_emulated);
   fChain->SetBranchAddress("gem_scluster_count", &gem_scluster_count, &b_gem_scluster_count);
   fChain->SetBranchAddress("gem_scluster_x", &gem_scluster_x, &b_gem_scluster_x);
   fChain->SetBranchAddress("gem_scluster_y", &gem_scluster_y, &b_gem_scluster_y);
   fChain->SetBranchAddress("gem_scluster_energy", &gem_scluster_energy, &b_gem_scluster_energy);
   fChain->SetBranchAddress("gem_scluster_adc", &gem_scluster_adc, &b_gem_scluster_adc);
   fChain->SetBranchAddress("srs_prerecon_count", &srs_prerecon_count, &b_srs_prerecon_count);
   fChain->SetBranchAddress("srs_prerecon_y", &srs_prerecon_y, &b_srs_prerecon_y);
   fChain->SetBranchAddress("srs_prerecon_x", &srs_prerecon_x, &b_srs_prerecon_x);
   fChain->SetBranchAddress("gem_peak_count", &gem_peak_count, &b_gem_peak_count);
   fChain->SetBranchAddress("gem_peak_plane_id", &gem_peak_plane_id, &b_gem_peak_plane_id);
   fChain->SetBranchAddress("gem_peak_plane_name", &gem_peak_plane_name, &b_gem_peak_plane_name);
   fChain->SetBranchAddress("gem_peak_index", &gem_peak_index, &b_gem_peak_index);
   fChain->SetBranchAddress("gem_peak_apv_id", &gem_peak_apv_id, &b_gem_peak_apv_id);
   fChain->SetBranchAddress("gem_peak_height", &gem_peak_height, &b_gem_peak_height);
   fChain->SetBranchAddress("gem_peak_width", &gem_peak_width, &b_gem_peak_width);
   fChain->SetBranchAddress("gem_peak_area", &gem_peak_area, &b_gem_peak_area);
   fChain->SetBranchAddress("gem_peak_real_pos", &gem_peak_real_pos, &b_gem_peak_real_pos);
   Notify();
}

Bool_t trdclass_halld24::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void trdclass_halld24::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t trdclass_halld24::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

std::pair<Double_t, Double_t> trdclass_halld24::TrkFit(TH2F *h2_evt, TF1 &fx, const char *cfx, int rob )
{
  //----------  SF fit ---------------------------
  /*
    h1f->Fit("gaus")
    TF1 * f = h1f->GetFunction("gaus")
    f->GetNDF()
    f->GetChisquare()
    f->GetProb()
    Int_t bin = h3->GetBin(binx,biny,binz);
    Float_t y = h3->GetBinContent(bin);
    virtual Double_t TH2::GetBinContent     (       Int_t   binx,           Int_t   biny    )
  */
  gErrorIgnoreLevel = kBreak; // Suppress warning messages from empty fit data
  
  TCutG *cutgx = new TCutG("cutgx",5);
  cutgx->SetPoint(0,40,0);
  cutgx->SetPoint(1,150,0);
  cutgx->SetPoint(2,150,528);
  cutgx->SetPoint(3,40,528);
  cutgx->SetPoint(4,40,0);
  
  TProfile *profx = h2_evt->ProfileX("profx", 0, 528, "[cutgx]");
  if (rob>0) {
    profx->Fit(cfx,"QNR+rob=0.75"); //  "+rob=0.75"
  } else {
    profx->Fit(cfx,"QNR");
  }
  Double_t chi2x = fx.GetChisquare();
  Double_t Ndfx = fx.GetNDF();
  Double_t integral = profx->Integral(profx->GetXaxis()->FindBin(40.), profx->GetXaxis()->FindBin(150.));
  //Double_t p0x = fx.GetParameter(0);
  //Double_t p1x = fx.GetParameter(1);
  
  //int kfit = 0;
  ////if (chi2x/Ndfx<100 && chi2y/Ndfy<10 && Ndfx>10 && Ndfy>10) {
  //if (chi2x/Ndfx<100 && Ndfx>10) {
  //  kfit=1;
  //}
  double chi2=chi2x/Ndfx; if (Ndfx<3) chi2=-chi2;
  //return chi2;
  return std::make_pair(chi2, integral);
}


//==================================================================
void trdclass_halld24::Count(const char *tit) {
  hcount->Fill(tit,1);
}
void trdclass_halld24::Count(const char *tit, double cut1) {
  char clab[20];
  sprintf(clab,"%s_%.1f",tit,cut1);
  hcount->Fill(clab,1);
}
void  trdclass_halld24::Count(const char *tit, double cut1, double cut2) {
  char clab[20];
  sprintf(clab,"%s_%.1f_%.1f",tit,cut1,cut2);
  hcount->Fill(clab,1);
}
//------------------------------------------------------------------




#endif // #ifdef trdclass_halld24_cxx
