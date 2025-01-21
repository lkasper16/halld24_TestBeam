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
#include "TF1.h"
#include "TPaveStats.h"
#include "TMultiGraph.h"
#include "TSystem.h"
#include "stdio.h"
#include "trd_mlp_halld24.h"
#include "PlotLib.C"

//--  0-dEdx 1-Params 2-dEdx+Params ; change MAXpar to 17 !
#define NN_MODE 3
//#define VERBOSE
#define ANALYZE_MERGED 1
#define MMG_RUN
//#define USE_CLUSTERS
#define FERMI_NN
#define NO_RAD_COMPARE 1

void Count(const char *tit);
void Count(const char *tit, double cut1);
void Count(const char *tit, double cut1, double cut2);
TH1D *hcount;

const int MaxNDEslices =10;
const int NFixed = 8;

#if  NN_MODE == 3
  const int NDEslices = 10;
  const int MAXpar = NFixed+NDEslices;
#elif  NN_MODE == 4
  const int NDEslices = 4;
  const int MAXpar = NFixed+NDEslices;
#else 
ne rabotaet !!!!
#endif

Float_t Par[MAXpar];
Int_t type, ievent, channel;
TCanvas *c2=NULL;
TH2F *dispe, *disppi;

//----------------------------------------------------------------------
int hgive1(TH1 *hp, int *nx, double *xmi, double *xma) {
  *nx=hp->GetNbinsX();
  *xmi=hp->GetBinLowEdge(1);
  *xma=hp->GetBinLowEdge(*nx)+hp->GetBinWidth(*nx);
  
  return hp->GetEntries();
}

//-------------------------------------------------------------------------------
int hscale(TH1 *he, TH1 *hpi, double scale, int NORM, int DRAW) {
  int ret=0;
  int noent_e= he->GetEntries(); 
  int noent_pi= hpi->GetEntries(); 
  double escale = 1; 
  
  if (scale>0) escale=scale;
  else if (noent_e>0) escale=(double)noent_pi/(double)noent_e; 
  else    ret=1; // -- 1 = error
  string name = hpi->GetName();
  #ifdef VERBOSE
    printf("NORM :: hist=%s<<< noent e=%d pi=%d e-scale=%f \n",name.c_str(),noent_e,noent_pi,escale);
  #endif
  
  if (NORM) he->Scale(escale);
  if (DRAW>0) {
    double maxe = he->GetMaximum(); 
    double maxpi = hpi->GetMaximum(); 
    if (maxe>maxpi) {
      he->Draw("hist");  hpi->Draw("histsames");
    } else {
      hpi->Draw("hist");  he->Draw("histsames");
    }
    he->SetLineColor(2);  hpi->SetLineColor(4);
  }
  
  if (DRAW>1) {
    gPad->Update();  
    TH1 *h1 = he;
    TH1 *h0 = hpi; 
    TPaveStats *ps1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
    ps1->SetY1NDC(0.57);  ps1->SetY2NDC(0.75); ps1->SetTextColor(kRed);
    TPaveStats *ps0 = (TPaveStats*)h0->GetListOfFunctions()->FindObject("stats");
    ps0->SetTextColor(kBlue);
    gPad->Modified(); gPad->Update(); 
  }
  return ret;
}

//----------------------- Rejection Factor Calculation ---------------------------------
std::pair<double,double> Reject(TH1 *hp, TH1 *he, double thr) {
  double s[202][2], r1, w, Rej, XMI, XMA, es, ps, e1=0., e2=0., p1=0., p2=0.;
  int NX, i0, noente=0, noentp=0;
  Rej=-1; r1=1-thr; i0=thr; w=thr-i0;
  //............ electron 10% of integral level ..........................
  noente=hgive1(he,&NX,&XMI,&XMA);
  TString etitle = he->GetTitle();
  #ifdef VERBOSE
    cout << "===> " << etitle << " <===  NX e = " << NX <<  " w = " << w <<  " i0 = " << i0 << endl;
  #endif
  if (NX>200) { cout << "error NX e = " << NX << endl; exit(1); }
  if (noente>0) {
    es=0;
    for (int i=0; i<=NX+1; i++) {
      es=es+he->GetBinContent(i);
      s[i][1]=es;
    }
    e2=0;
    #ifdef VERBOSE
      cout << " es = " << es << " noente= " << noente << endl;
    #endif
    if (es>0) {
      for (int i=0; i<=NX+1; i++) { 
	      e1=e2;
        e2=s[i][1]/es;  
	      if(e2>r1) {
          i0=i-1;
          break;
        } 
      }
      w=(r1-e1)/(e2-e1);
      Rej = i0+w;
    }
  }
  //............ now pions level at this threshold .......................
  noentp=hgive1(hp,&NX,&XMI,&XMA);
  #ifdef VERBOSE
    cout << "NX pi = " << NX <<  " i0=" << i0 << " thr=" << he->GetBinLowEdge(i0) << endl;
  #endif
  double relRejError = 0.;
  if (noentp>0) {
    ps=0;
    for (int i=0; i<=NX+1; i++) {
      ps=ps+hp->GetBinContent(i);
      s[i][0]=ps;
    }
    p1=s[i0][0]/ps;
    p2=s[i0+1][0]/ps;
    Rej = 1-(p1+(p2-p1)*w);
    double rejError = sqrt(ps-s[i0][0]);
    relRejError = 1./rejError;
    cout << "===> " << etitle << " <=== Relative Error = " << relRejError << endl;
  }
  cout << "===> " << etitle << " <=== Rejection = " << 1./Rej << endl;
  return std::make_pair(Rej, relRejError);
}

//---------------------------------------------------------------------
#if ANALYZE_MERGED
int fill_trees(TTree *ttree_hits, TTree *signal, TTree *background, TTree *sig_tst, TTree *bg_tst, int RunNum, int *rtw1, int *rtw3, int nEntries) {
#else
int fill_trees(TTree *ttree_hits, TTree *signal, TTree *background, TTree *sig_tst, TTree *bg_tst, int RunNum, int *rtw1, int *rtw3) {
#endif
  // Set object pointer
  ecal_energy=0;
  presh_energy=0;
  mult_energy=0;
  gem_nhit = 0;
  gem_nclu = 0;
  mmg1_nhit = 0;
  mmg1_nclu = 0;
  xpos = 0;
  zpos = 0;
  dedx = 0;
  parID = NULL;
  zHist = 0;
  xposc = 0;
  zposc = 0;
  dedxc = 0;
  widthc = 0;
  xposc_max = 0;
  zposc_max = 0;
  dedxc_max = 0;
  widthc_max = 0;
  
  // Set branch addresses and branch pointers
  if (!ttree_hits) return -1;
  TTree *fChain;
  fChain = ttree_hits;
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress("event_num", &event_num, &b_event_num);
  fChain->SetBranchAddress("ecal_energy", &ecal_energy, &b_ecal_energy);
  fChain->SetBranchAddress("presh_energy", &presh_energy, &b_presh_energy);
  fChain->SetBranchAddress("mult_energy", &mult_energy, &b_mult_energy);
  fChain->SetBranchAddress("parID", &parID, &b_parID);
  fChain->SetBranchAddress("xpos", &xpos, &b_xpos);
  fChain->SetBranchAddress("zpos", &zpos, &b_zpos);
  fChain->SetBranchAddress("dedx", &dedx, &b_dedx);
  fChain->SetBranchAddress("zHist", &zHist, &b_zHist);
  fChain->SetBranchAddress("xposc", &xposc, &b_xposc);
  fChain->SetBranchAddress("zposc", &zposc, &b_zposc);
  fChain->SetBranchAddress("dedxc", &dedxc, &b_dedxc);
  fChain->SetBranchAddress("widthc", &widthc, &b_dedxc);
  fChain->SetBranchAddress("xposc_max", &xposc_max, &b_xposc_max);
  fChain->SetBranchAddress("zposc_max", &zposc_max, &b_zposc_max);
  fChain->SetBranchAddress("dedxc_max", &dedxc_max, &b_dedxc_max);
  fChain->SetBranchAddress("widthc_max", &widthc_max, &b_widthc_max);
  #ifdef MMG_RUN
    fChain->SetBranchAddress("nhit", &mmg1_nhit, &b_mmg1_nhit);
    fChain->SetBranchAddress("nclu", &mmg1_nclu, &b_mmg1_nclu);
  #else
    fChain->SetBranchAddress("nhit", &gem_nhit, &b_gem_nhit);
    fChain->SetBranchAddress("nclu", &gem_nclu, &b_gem_nclu);
  #endif
  //========================================
  //TList *HistList = new TList();
  TH2F *hits2d_e = new TH2F("hits2d_e","hits2d_e",125,0.5,250.5,240,-0.5,239.5);
  TH2F *hits2d_pi = new TH2F("hits2d_pi","hits2d_pi",125,0.5,250.5,240,-0.5,239.5);
  
  TH2F *aver2d_e = new TH2F("aver2d_e","aver-rms electrons",120,0.,240.,100,0.,10.);
  TH2F *aver2d_pi = new TH2F("aver2d_pi","aver-rms pions",120,0.,240.,100,0.,10.);
  
  TH1F *hNhits  = new TH1F("hNhits","",70,-0.5,69.5);  hNhits->SetStats(0);
  TH1F *hNclu  = new TH1F("hNclu","",70,-0.5,69.5);  hNclu->SetStats(0);
  
  TH1F *helectron_maxamp = new TH1F("helectron_maxamp","Electron max amp",200,0.5,4000.5);
  TH1F *hpion_maxamp = new TH1F("hpion_maxamp","Pion max amp",200,0.5,4000.5);
  TH1F *helectron_dedxtotal = new TH1F("helectron_dedxtotal","Electron e_total",600,0.5,150000.5);
  TH1F *hpion_dedxtotal = new TH1F("hpion_dedxtotal","Pion e_total",600,0.5,150000.5);
  
  TH1F *time_e  = new TH1F("time_e","Electron Amplitude in Time",250,0.5,250.5);
  TH1F *time_pi = new TH1F("time_pi","Pion Amplitude in Time",250,0.5,250.5);
  
  TH1F *par_e[MAXpar]; 
  TH1F *par_pi[MAXpar];
  
  for (int ip=0; ip<MAXpar; ip++) {
    char hname[80];
    sprintf(hname,"par_e_%d ",ip);
    par_e[ip] = new TH1F(hname,hname,100,-0.5,99.5);
    sprintf(hname,"par_pi_%d ",ip);
    par_pi[ip] = new TH1F(hname,hname,100,-0.5,99.5);
  }
  
  //========================================
  TRandom3 *rndm = new TRandom3();
  int nx;
  double xmi, xma;
  int noent = hgive1(par_e[0], &nx, &xmi, &xma);
  #ifdef VERBOSE
    cout << " hist: nx=" << nx << " xmi=" << xmi << " xma=" << xma << " noent= " << noent << endl;
  #endif
  Long64_t nentries = ttree_hits->GetEntries();
  Long64_t nbytes = 0; //////
  int ntest=nentries*0.25; ////////
  #ifdef VERBOSE
    //cout << " ntest=" << ntest << endl;
  #endif
  
  //==============================================================================
  const int NDE=MaxNDEslices;
  double dEdx[NDE];
  int NPF=MAXpar; // number of parameters filled
  int ntrk_e=0, ntrk_pi=0;
  int e_chan1=0;    //-- first TR channel
  int e_chan2=0;    //-- last  TR channel
  int pi_chan1=0;   //-- first pion (no-rad) channlel
  int pi_chan2=0;   //-- last  pion (no-rad) channel
  
  cout << "============== BEGIN RUN " << RunNum << "===============" <<endl;
  //--------------------------------------------------------------------------------
  //                   E V E N T    L O O P
  //--------------------------------------------------------------------------------
  for (Long64_t iev = 0; iev < nentries; iev++) {  
    nbytes += ttree_hits->GetEntry(iev); ///////
    ievent=iev;
    Count("EVT");
    if (!(iev % 1000)) cout<<" ===> Event = "<<iev<<" out of "<<nentries<<endl;
    for (int ip=0; ip<MAXpar; ip++) Par[ip]=0;
    for (int i=0; i<NDE; i++) dEdx[i]=0;
    channel = 22;
    //------------------------------------------------------------------------------
    //                 G E M 
    //------------------------------------------------------------------------------
    
    int USE_TRACK = 0;
    /*
    if ( USE_TRACK ) {
      if (abs(13-channel/27-trkch)<3.) { // tracking cut 
	      h2xdiff->Fill(channel,trkch);
	      hbeamX->Fill(channel);
      }
    } else {
      h2xdiff->Fill(channel,trkch);
      hbeamX->Fill(channel);
    }
    */
    float amax2=-1.;
    float emax2=-1.;
    
    int khit=0;
    int NTR=0;
    
    float THR1=100.;
    float THR2=500.;
    float Escale=400.; //1.
    float Ascale=40.; //1.
       
    float atot=0.;
    float etot=0.;
    float etrzon=0.;
    
    int tw1=0;
    int tw2=tw1+1;
    int tw3=tw2+1;
    
    //---- run specific ---
    if ( 5000 < RunNum && RunNum <= 5999 ) { //------------------  CERN   2024  ---------------------
      
      Escale=10.;
      Ascale=10.;       
      e_chan1=100;        //-- first TR channel
      e_chan2=130;        //-- last  TR channel
      pi_chan1=e_chan1;   //-- first pion (no-rad) channel
      pi_chan2=e_chan2;   //-- last  pion (no-rad) channel
      tw1 = 60;
      tw2 = 95;
      tw3 = 140;
      
      switch (RunNum) {
        
        #ifdef MMG_RUN
        case 5264:   tw1=46; tw2=95;  tw3=193; e_chan1=102;   e_chan2=127;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- 15cm Fleece
        case 5268:   tw1=46; tw2=95;  tw3=193; e_chan1=102;   e_chan2=127;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- No Rad
        case 5278:   tw1=46; tw2=95;  tw3=193; e_chan1=104;   e_chan2=127;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- 15cm Fleece
        case 5283:   tw1=46; tw2=95;  tw3=193; e_chan1=102;   e_chan2=127;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- 23cm Fleece
        case 5284:   tw1=46; tw2=95;  tw3=193; e_chan1=102;   e_chan2=127;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- 23cm Fleece
        
        //-- For Second Xe Bottle !!
        case 5301:   tw1=49; tw2=105; tw3=193; e_chan1=102; e_chan2=128;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- No Rad
        case 5303:   tw1=49; tw2=105; tw3=185; e_chan1=101; e_chan2=128;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- Foil
        case 5304:   tw1=49; tw2=105; tw3=185; e_chan1=101; e_chan2=128;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- Foil
        case 5306:   tw1=55; tw2=105; tw3=193; e_chan1=102; e_chan2=128;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- Foil

        #else
        case 5252:   tw1=65; tw2=115; tw3=135; e_chan1=104;   e_chan2=127;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- 15cm Fleece
        case 5254:   tw1=65; tw2=118; tw3=134; e_chan1=112;   e_chan2=127;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- 15cm Fleece
        case 5256:   tw1=65; tw2=118; tw3=134; e_chan1=112;   e_chan2=127;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- 15cm Fleece
        case 5264:   tw1=65; tw2=90;  tw3=141; e_chan1=109;   e_chan2=130;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- 15cm Fleece
        case 5268:   tw1=67; tw2=95;  tw3=141; e_chan1=109;   e_chan2=129;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- No Rad
        case 5278:   tw1=61; tw2=90;  tw3=140; e_chan1=106;   e_chan2=133;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- 15cm Fleece
        case 5281:   tw1=67; tw2=104; tw3=141; e_chan1=104;   e_chan2=127;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- 23cm Fleece
        case 5282:   tw1=67; tw2=104; tw3=141; e_chan1=104;   e_chan2=129;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- 23cm Fleece
        case 5283:   tw1=61; tw2=90;  tw3=140; e_chan1=106;   e_chan2=133;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- 23cm Fleece
        case 5284:   tw1=61; tw2=90;  tw3=140; e_chan1=106;   e_chan2=133;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- 23cm Fleece
        
        //-- For Second Xe Bottle !!
        case 5301:   tw1=68; tw2=105; tw3=170; e_chan1=106; e_chan2=133;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- No Rad
        case 5302:   tw1=70; tw2=135; tw3=170; e_chan1=104; e_chan2=127;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- Foil
        case 5303:   tw1=68; tw2=105; tw3=170; e_chan1=106; e_chan2=133;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- Foil
        case 5304:   tw1=70; tw2=135; tw3=170; e_chan1=104; e_chan2=127;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- Foil
        case 5306:   tw1=68; tw2=105; tw3=170; e_chan1=106; e_chan2=133;  pi_chan1=e_chan1; pi_chan2=e_chan2;   break; //-- Foil
        #endif
        
        default:     tw1=60; tw2=95;  tw3=140; e_chan1=100; e_chan2=130;  pi_chan1=e_chan1; pi_chan2=e_chan2;
      }
    } else {
      printf("Run %d not found in the CERN 2024 range - exiting \n",RunNum);
      exit(1);
    }
    
    //----------------------------------------------------------------
    //             P I D   S E L E C T I O N
    //----------------------------------------------------------------
    
    double dt=(tw3-tw1)/NDE;
    
    #if 1
      type=-1;
      #if NO_RAD_COMPARE
        if (parID->size()>0) {
          if(parID->at(0)==1 && ecal_energy > 4000. && presh_energy > 1000. ) {
            type=1; ntrk_e++; Count("ntrk_e");
          } else  if(parID->at(0)==0) {
            type=0; ntrk_pi++; Count("ntrk_pi");
          }
        }
      #else
        if (parID->size()>0) {
          if(parID->at(0)==1 && ecal_energy > 4000. && presh_energy > 1000. ) {
	          type=1; ntrk_e++; Count("ntrk_e");
          } else  if(parID->at(0)==0 && ecal_energy < 400. && presh_energy < 400. ) { 
	          type=0; ntrk_pi++; Count("ntrk_pi");
          }
        }
      #endif
    #else
      type=-1;
      if(ecal_energy > 4000.) {
        type=1; ntrk_e++;
      } else if(ecal_energy < 400.) {
        type=0; ntrk_pi++;
      }
    #endif
    
    //----------------------------------------------------------------
    //             C L U S T E R S 
    //----------------------------------------------------------------
    /*
    int clu_nhits=xposc->size();
    float max_widthc=0, max_dedxc=0;
    for (int icl=0; icl<clu_nhits; icl++) {
      if (widthc->at(icl)>max_widthc) max_widthc=widthc->at(icl);
      if (dedxc->at(icl)>max_dedxc) max_dedxc=dedxc->at(icl);
    }
    */
    #ifdef MMG_RUN
      int clu_nhits=mmg1_nclu;
      hNhits->Fill(mmg1_nhit);
      hNclu->Fill(mmg1_nclu);
    #else
      int clu_nhits=gem_nclu;
      hNhits->Fill(gem_nhit);
      hNclu->Fill(gem_nclu);
    #endif
    float max_widthc=widthc_max;
    float max_dedxc=dedxc_max;
    
    //=================== Loop to calculate average & RMS beam position ========
    double xaver=0, xaver2=0;
    int naver=0;
    #ifdef MMG_RUN
      #ifdef USE_CLUSTERS
        for (int i=0;i<mmg1_nhit;i++){ // count fixed parameters 
          if (tw1 > zposc->at(i) || zposc->at(i) > tw3) continue;
          xaver+=xposc->at(i); naver++;
        }
      #else
        for (int i=0;i<mmg1_nhit;i++){ // count fixed parameters 
          if (tw1 > zpos->at(i) || zpos->at(i) > tw3) continue;
          xaver+=xpos->at(i); naver++; 
        }
      #endif
    #else
      #ifdef USE_CLUSTERS
        for (int i=0;i<gem_nhit;i++){ // count fixed parameters 
          if (tw1 > zposc->at(i) || zposc->at(i) > tw3) continue;
          xaver+=xposc->at(i); naver++;
        }
      #else
        for (int i=0;i<gem_nhit;i++){ // count fixed parameters 
          if (tw1 > zpos->at(i) || zpos->at(i) > tw3) continue;
          xaver+=xpos->at(i); naver++;
        }
      #endif
    #endif
    xaver=xaver/naver;
    
    //====================================================================
    #ifdef MMG_RUN
    for (int i=0;i<mmg1_nhit;i++) { // count fixed parameters
      Count("mmg1_Hits");
      #ifdef USE_CLUSTERS
        if (tw1 > zposc->at(i) || zposc->at(i) > tw3) continue;
      #else
        if (tw1 > zpos->at(i) || zpos->at(i) > tw3) continue;
      #endif
      Count("mmg1_zHits");
    #else
    for (int i=0;i<gem_nhit;i++) { // count fixed parameters 
      Count("Gem_Hits");
      #ifdef USE_CLUSTERS
        if (tw1 > zposc->at(i) || zposc->at(i) > tw3) continue;
      #else
        if (tw1 > zpos->at(i) || zpos->at(i) > tw3) continue;
      #endif
      Count("Gem_zHits");
    #endif
      
    #ifdef USE_CLUSTERS
      xaver2+=((xposc->at(i)-xaver)*(xposc->at(i)-xaver));
      if (dedxc->at(i)>THR1) {
  	    if (type==1) hits2d_e->Fill(zposc->at(i),xposc->at(i),dedxc->at(i));
	      else if (type==0) hits2d_pi->Fill(zposc->at(i),xposc->at(i),dedxc->at(i));
      }
      if(dedxc->at(i)>amax2 && tw1<zposc->at(i) && zposc->at(i)<tw3) { //-- tw2 or tw1 ????
	      amax2=dedxc->at(i);
      }
      if(dedxc->at(i)>emax2) {  // currently same as amp
	      emax2=dedxc->at(i);
      }
      
      if (type>=0) {  //-- rad.pos. window OK
	      if (tw1<zposc->at(i) && zposc->at(i)<tw3 && dedxc->at(i)>THR1) {
	        khit++; etot+=dedxc->at(i); atot+=dedxc->at(i);
	        int ibin=(zposc->at(i)-tw1)/dt; ibin=min(max(0,ibin),(NDE-1)); dEdx[ibin]+=dedxc->at(i)/10;
	      }
	      if (tw2 < zposc->at(i) && zposc->at(i) < tw3) { etrzon+=dedxc->at(i); }
	      if (dedxc->at(i)>THR2) { NTR++; Count("NTR"); }
      }
      
      if (dedxc->at(i)>THR1) {
	      if (type==1) time_e->Fill(zposc->at(i),dedxc->at(i));
	      if (type==0) time_pi->Fill(zposc->at(i),dedxc->at(i));
      }
    #else
    xaver2+=((xpos->at(i)-xaver)*(xpos->at(i)-xaver));
      if (dedx->at(i)>THR1) {
        if (type==1) hits2d_e->Fill(zpos->at(i),xpos->at(i),dedx->at(i));
        else if (type==0) hits2d_pi->Fill(zpos->at(i),xpos->at(i),dedx->at(i));
      }
      if(dedx->at(i)>amax2 && tw1<zpos->at(i) && zpos->at(i)<tw3) { //-- tw2 or tw1 ????
        amax2=dedx->at(i);
      }
      if(dedx->at(i)>emax2) {  // currently same as amp
        emax2=dedx->at(i);
      }
      
      if (type>=0 ) {  //-- rad.pos. window OK
        if (tw1<zpos->at(i) && zpos->at(i)<tw3 && dedx->at(i)>THR1) {
          khit++; etot+=dedx->at(i); atot+=dedx->at(i);
          int ibin=(zpos->at(i)-tw1)/dt;  ibin=min(max(0,ibin),(NDE-1)); dEdx[ibin]+=dedx->at(i)/10;
        }
        if (tw2 < zpos->at(i) && zpos->at(i) < tw3)  {  etrzon+=dedx->at(i); }
        if (dedx->at(i)>THR2) { NTR++;  Count("NTR"); }
      }

      if (dedx->at(i)>THR1)  {
        if (type==1) time_e->Fill(zpos->at(i),dedx->at(i));
        if (type==0) time_pi->Fill(zpos->at(i),dedx->at(i));
      }
    #endif
    } //--- End Loop over detector nhits
    
    //=================================================================
    
    xaver2=sqrt(xaver2/naver); 
    if (type==1 ) { 
      Count("el");
      aver2d_e->Fill(xaver,xaver2); 
    } else {
      Count("pi");
      aver2d_pi->Fill(xaver,xaver2);
    }
    if (e_chan1 > xaver || xaver > e_chan2) continue; //--- radiator  area  ; for Y - need a track
    Count("radAreaHits");
    
    //--------------------------------------------------------------------------------
    //                    electron case
    //--------------------------------------------------------------------------------
    if (type==1) {
      helectron_maxamp->Fill(amax2);
      helectron_dedxtotal->Fill(etot);
    }
    //--------------------------------------------------------------------------------
    //                    pion case
    //--------------------------------------------------------------------------------
    if (type==0) {
      hpion_maxamp->Fill(amax2);
      hpion_dedxtotal->Fill(etot);
    }
    //-----------------------------------------------
    if (type<0) continue;
    Count("type_ok");
    
    if (NN_MODE==3) {  //-- dEdx(amp) + Par
      if (MAXpar<(NFixed+NDE)) { printf("ERROR :: MAXpar array too small =%d \n",MAXpar); exit(1); }
      Par[0]=amax2/Ascale/5.;
      Par[1]=clu_nhits*5;
      Par[2]=xaver2*8;
      Par[3]=atot/Escale/50.; //atot or etot ???
      Par[4]=etrzon/Escale/50.;
      Par[5]=NTR;
      Par[6]=max_widthc*5.;
      Par[7]=max_dedxc/Ascale/5.;
      int np=NDE;
      double coef=Ascale/2.;
      for (int ip=0; ip<np; ip++) { 
	      Par[ip+NFixed]=dEdx[ip]/coef; 
      }
    } else if (NN_MODE==4) {  //-- dEdx(amp) + Par
      Par[0]=amax2/Ascale/5.;
      Par[1]=khit;  // -- hits > THR1 = 100
      Par[2]=xaver2*8;
      Par[3]=atot/Escale/50.;
      Par[4]=etrzon/Escale/50.;
      Par[5]=NTR;
      Par[6]=max_widthc;
      Par[7]=max_dedxc/Ascale/5.;
      int np=NDEslices; // NDE;
      double coef=Ascale/2.;
      for (int ip=0; ip<np; ip++) { 
	      Par[ip+NFixed]=dEdx[NDE-1-ip]/coef; 
      }
    } else {   //-- dEdx only
      int np=min(MAXpar,NDE);
      double coef=Ascale*3.;
      for (int ip=0; ip<np; ip++) { 
	      Par[ip]=dEdx[ip]/coef; 
      }
    }
    
    *rtw1=tw1;
    *rtw3=tw3;
    
    if (type==1) {    // electron
      if (rndm->Rndm()<0.1)   sig_tst->Fill();
      else                    signal->Fill();
      for (int ip=0; ip<NPF; ip++) {
        par_e[ip]->Fill(Par[ip]);
      }
    }
    if (type==0) {    // pion
      if (rndm->Rndm()<0.1)   bg_tst->Fill();
      else                    background->Fill();
      for (int ip=0; ip<NPF; ip++) {
        par_pi[ip]->Fill(Par[ip]);
      }
    }
    
  }//-------------------- End event loop --------------------
  
  //------------------ Plot Histograms -----------------------------
  
  double escale_trk=1;   if (ntrk_e>0) escale_trk=(double)ntrk_pi/(double)ntrk_e;
  #ifdef VERBOSE
    printf("escale_trk=%f\n",escale_trk);
  #endif
  int NORM=1; //-- scale hist using num entries
  int DRAW=2; // 1=draw 2= two stat boxes
  int nxd=3;
  int nyd=5;
  int COMPACT=0;
  TCanvas *c1,*c0; 
  char ctit[120];
  #if ANALYZE_MERGED
    #ifdef MMG_RUN
      sprintf(G_DIR,"mlpOutput/halld24/merged/hd_rawdata_mmg_%06d_%06dEntries.root",RunNum,nEntries);
    #else
      sprintf(G_DIR,"mlpOutput/halld24/merged/hd_rawdata_%06d_%06dEntries.root",RunNum,nEntries);
    #endif
  #else
    #ifdef MMG_RUN
      sprintf(G_DIR,"mlpOutput/halld24/hd_rawdata_mmg_%06d.root",RunNum);
    #else
      sprintf(G_DIR,"mlpOutput/halld24/hd_rawdata_%06d.root",RunNum);
    #endif
  #endif
  sprintf(ctit,"File=%s",G_DIR);
  htitle(ctit);
  
  c1=NextPlot(nxd,nyd);   hits2d_e->Draw("colz"); 
  TLine *lin1 = new TLine(0.,e_chan1,250.,e_chan1);    TLine *lin2 = new TLine(0.,e_chan2,250.,e_chan2);
  lin1->SetLineColor(kRed);   lin2->SetLineColor(kRed);   lin1->Draw();  lin2->Draw();
  gPad->Modified(); gPad->Update();
  
  c1=NextPlot(nxd,nyd);   hits2d_pi->Draw("colz");
  TLine *lin1p = new TLine(0.,pi_chan1-1,250.,pi_chan1); TLine *lin2p = new TLine(0.,pi_chan2-1,250.,pi_chan2);
  lin1p->SetLineColor(kCyan); lin2p->SetLineColor(kCyan); lin1p->Draw(); lin2p->Draw();
  gPad->Modified(); gPad->Update();
  
  #ifdef VERBOSE
    printf(" Draw Lines :: %d %d %d %d \n",e_chan1,e_chan2, pi_chan1, pi_chan2);
  #endif
  
  c1=NextPlot(nxd,nyd);  aver2d_e->Draw("colz");
  c1=NextPlot(nxd,nyd);  aver2d_pi->Draw("colz");
  c1=NextPlot(nxd,nyd);  hNclu->SetLineColor(6);  hNclu->Draw();  hNhits->SetLineColor(209);  hNhits->Draw("same");
  TLegend *l1 = new TLegend(.7, .7, .9, .9);
  #ifdef MMG_RUN
  l1->AddEntry(hNhits, "mmg1_nhits");
  l1->AddEntry(hNclu, "mmg1_nclu");
  #else
  l1->AddEntry(hNhits, "gem_nhits");
  l1->AddEntry(hNclu, "gem_nclu");
  #endif
  l1->Draw();
  
  //---------------------------------------------------------------------
  cout << " +++++++  Rejection From Amplitude Distributions  +++++++" << endl;
  std::pair<double,double> rej70 = Reject(hpion_maxamp, helectron_maxamp, 0.7);
  cout << " Ampl: e = 70% , Eff pi = " << rej70.first*100. << "% , Rejection = " << 1./rej70.first << endl;
  std::pair<double,double> rej90 = Reject(hpion_maxamp, helectron_maxamp, 0.9);
  cout << " Ampl: e = 90% , Eff pi = " << rej90.first*100. << "% , Rejection = " << 1./rej90.first << endl;
  cout << "+++++++++++++++++++++++++++++++++++\n" << endl;
  
  c1=NextPlot(nxd,nyd);   hscale(helectron_maxamp,hpion_maxamp,0.,NORM,2);
  c1=NextPlot(nxd,nyd);   hscale(helectron_dedxtotal,hpion_dedxtotal,0.,NORM,2);
  c1=NextPlot(nxd,nyd);   hscale(time_e,time_pi,escale_trk,NORM,2); //--- scale time hist here ---
  
  for (int ip=0; ip<NPF; ip++) {
    c1=NextPlot(nxd,nyd);
    if (NN_MODE==0 || (NN_MODE > 1 && ip>=NFixed))  gPad->SetLogy();
    hscale(par_e[ip],par_pi[ip],escale_trk,1,2);
  }
  c1=NextPlot(-1,-1);
  
  //------------------------------------------------------------
  //char pngname[120];
  //sprintf(pngname,"%s_dqm_m%d.png",G_DIR,NN_MODE);
  //c1->Print(pngname);
  char pdfname[120];
  sprintf(pdfname,"%s_dqm_m%d.pdf",G_DIR,NN_MODE);
  c1->Print(pdfname);
  
  c0 = new TCanvas("Plot","Plot",400,200,1200,900);     c0->Divide(2,2);
  c0->cd(1);  hscale(time_e,time_pi,1.,0,2); // --- already scaled before 
  
  c0->cd(3);  hits2d_e->Draw("colz"); // --- already scaled before 
  lin1->Draw();  lin2->Draw();  
  gPad->Modified(); gPad->Update();
  
  c0->cd(4);  hits2d_pi->Draw("colz"); // --- already scaled before 
  lin1p->Draw(); lin2p->Draw();
  gPad->Modified(); gPad->Update();
  
  #ifdef VERBOSE
    printf("Draw Lines on c0: e1=%d e2=%d \n",e_chan1,e_chan2);
  #endif
  c0->cd();
  c0->Modified(); c0->Update();
  char pdfname0[120];
  sprintf(pdfname0,"%s_dqm_m%d_time.pdf",G_DIR,NN_MODE);
  c0->Print(pdfname0);
  
  //------------------------------------------------------------
  return NN_MODE;
} //--- End fill _trees() ---

//==================================================================

#if ANALYZE_MERGED
void trd_mlp_halld24(int RunNum, int nEntries) {
#else
void trd_mlp_halld24(int RunNum) {
#endif
  
  hcount= new TH1D("hcount","Count",3,0,3); 
  hcount->SetStats(0);   hcount->SetFillColor(38);   hcount->SetMinimum(1.);
  #if ROOT_VERSION_CODE > ROOT_VERSION(6,0,0)
    hcount->SetCanExtend(TH1::kXaxis);
  #else
    hcount->SetBit(TH1::kCanRebin);
  #endif
  
  gStyle->SetTitleFontSize(.06);
  gStyle->SetLabelSize(.05, "XY");
  gStyle->SetTitleSize(0.05,"XY");
  
  char rootfile[256];
  #if ANALYZE_MERGED
    #if NO_RAD_COMPARE
      sprintf(rootfile,"RootOutput/halld24/merged/nrc_trd_singleTrackHits_Run_%06d_%06dEntries.root",RunNum,nEntries);
    #else
      sprintf(rootfile,"RootOutput/halld24/merged/trd_singleTrackHits_Run_%06d_%06dEntries.root",RunNum,nEntries);
    #endif
  #else
    sprintf(rootfile,"RootOutput/halld24/trd_singleTrackHits_Run_%06d.root",RunNum);
  #endif
  char basename[200];
  char *hd = strstr(rootfile,"/"); 
  strncpy(basename,&hd[1],200-1);   char *dot= strstr(basename,"."); *dot=0;
  printf("Base of input file name = %s \n",basename);
  #ifdef FERMI_NN
    Int_t ntrain=50;    // epoch //-- FERMI
  #else
    Int_t ntrain=101;   // epoch //-- CERN
  #endif
  int Nmod=3;
  // Prepare inputs
  // The 2 trees are merged into one, and a "type" branch, 
  // equal to 1 for the signal and 0 for the background is added.
  TFile *inputFile = 0;
  inputFile = new TFile(rootfile);
  if (!inputFile) {
    cout<<"Input file "<<rootfile<<" does not exist - exiting..."<<endl;
    return;
  }
  #ifdef MMG_RUN
    TTree *ttree_hits = (TTree *) inputFile->Get("mmg1_hits");
  #else
    TTree *ttree_hits = (TTree *) inputFile->Get("gem_hits");
  #endif
  char mlpname[128];
  #if ANALYZE_MERGED
    #ifdef MMG_RUN
      sprintf(mlpname,"mlpOutput/halld24/merged/mlp_mmg_run%06d_%06dEntries.root",RunNum,nEntries);
    #else
      sprintf(mlpname,"mlpOutput/halld24/merged/mlp_run%06d_%06dEntries.root",RunNum,nEntries);
    #endif
  #else
    #ifdef MMG_RUN
      sprintf(mlpname,"mlpOutput/halld24/mlp_mmg_run%06d.root",RunNum);
    #else
      sprintf(mlpname,"mlpOutput/halld24/mlp_run%06d.root",RunNum);
    #endif
  #endif
  TFile* outputFile = new TFile(mlpname,"RECREATE");
  
  TTree *sig_tst = new TTree("sig_tst", "Filtered Events"); 
  TTree *bg_tst = new TTree("bg_tst", "Filtered Events"); 
  TTree *signal = new TTree("signal", "Filtered Events"); 
  TTree *background = new TTree("background", "Filtered Events"); 
  TTree *simu = new TTree("simu", "Filtered Events");
  #if ANALYZE_MERGED
    #ifdef MMG_RUN
      char ctit[120];   sprintf(ctit,"hd_rawdata_mmg_%06d_%06d",RunNum,nEntries);
    #else
      char ctit[120];   sprintf(ctit,"hd_rawdata_%06d_%06d",RunNum,nEntries);
    #endif
  #else
    #ifdef MMG_RUN
      char ctit[120];   sprintf(ctit,"hd_rawdata_mmg_%06d_000",RunNum); 
    #else
      char ctit[120];   sprintf(ctit,"hd_rawdata_%06d_000",RunNum);
    #endif
  #endif
  dispe =  new TH2F("dispe","disp e; time ; X strip",250,-0.5,249.5,150,50.5,200.5);
  disppi = new TH2F("disppi","disp pi; time ; X strip",250,-0.5,249.5,150,50.5,200.5);
  dispe->SetStats(0);
  disppi->SetStats(0);
  double fsiz=0.02;
  dispe->GetYaxis()->SetTitleSize(fsiz);   dispe->GetXaxis()->SetTitleSize(fsiz);  dispe->GetZaxis()->SetTitleSize(fsiz);
  dispe->GetYaxis()->SetLabelSize(fsiz);   dispe->GetXaxis()->SetLabelSize(fsiz);  dispe->GetZaxis()->SetLabelSize(fsiz);
  disppi->GetYaxis()->SetTitleSize(fsiz);  disppi->GetXaxis()->SetTitleSize(fsiz); disppi->GetZaxis()->SetTitleSize(fsiz);
  disppi->GetYaxis()->SetLabelSize(fsiz);  disppi->GetXaxis()->SetLabelSize(fsiz); disppi->GetZaxis()->SetLabelSize(fsiz);
  
  for (int ip=0; ip<MAXpar; ip++) {
    char bname[80], tname[80];
    sprintf(bname,"par%d",ip);     sprintf(tname,"par%d/F",ip);
    sig_tst->Branch(bname, &Par[ip], tname);
    bg_tst->Branch(bname, &Par[ip], tname);
    signal->Branch(bname, &Par[ip], tname);
    background->Branch(bname, &Par[ip], tname);
    simu->Branch(bname, &Par[ip], tname);
  }
  simu->Branch("type", &type, "type/I");
  simu->Branch("ievent", &ievent, "ievent/I");
  sig_tst->Branch("ievent", &ievent, "ievent/I");
  bg_tst->Branch("ievent", &ievent, "ievent/I");
  sig_tst->Branch("channel", &channel, "channel/I");
  bg_tst->Branch("channel", &channel, "channel/I");
  signal->Branch("ievent", &ievent, "ievent/I");
  background->Branch("ievent", &ievent, "ievent/I");
  
  #ifdef VERBOSE
    cout << " Get trees signal=" << signal << endl;
  #endif
  int rtw1, rtw3;
  #if ANALYZE_MERGED
    int nn_mode = fill_trees(ttree_hits, signal, background, sig_tst, bg_tst, RunNum, &rtw1, &rtw3, nEntries);
  #else
    int nn_mode = fill_trees(ttree_hits, signal, background, sig_tst, bg_tst, RunNum, &rtw1, &rtw3);
  #endif
  
  #ifdef VERBOSE
    signal->Print();
    background->Print();
    sig_tst->Print();
    bg_tst->Print();
  #endif
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  //-----------------------------------------------
  type=1;
  Int_t i;
  for (i=0; i<signal->GetEntries(); i++) {
    signal->GetEntry(i);
    simu->Fill();
  }
  //-----------------------------------------------
  type=0;
  for (i=0; i<background->GetEntries(); i++) {
    background->GetEntry(i);
    simu->Fill();
  }
  //-----------------------------------------------
  // Build and train the NN par1 is used as a weight since we are primarly 
  // interested  by high pt events.
  // The datasets used here are the same as the default ones.
  //----------------------------
  const int NPAR=MAXpar;
  string INL, NNcfg;
  for (int il=0; il<NPAR; il++) {
    stringstream ss;  ss << il;  string si = ss.str();
    INL=INL+"@par"+si;   if (il<(NPAR-1)) INL=INL+",";
  }
  #ifdef FERMI_NN
    NNcfg=INL+":25:8:type"; //-- FERMI
  #else
    NNcfg=INL+":35:5:type"; //-- CERN
  #endif
  #ifdef VERBOSE 
    cout<<" INL="<<INL<<endl;
    cout<<" NNcfg="<<NNcfg<<endl;
  #endif
  TMultiLayerPerceptron *mlp = new TMultiLayerPerceptron(NNcfg.data(),simu,"Entry$%2","(Entry$+1)%2");
    //new TMultiLayerPerceptron(NNcfg.data(),simu,"Entry$%2","(Entry$+1)%2", TNeuron::kTanh);
    //new TMultiLayerPerceptron("@par1,@par2,@par3,@par4,@par5,@par6:16:8:type",
    //new TMultiLayerPerceptron("par1,par2,par3,par4,par5,par6:16:8:type",
    //new TMultiLayerPerceptron("@par1,@par3:8:3:type",
    //                          simu,"Entry$%2","(Entry$+1)%2");
  //----------------------------
  
  // The learning method is defined using the TMultiLayerPerceptron::SetLearningMethod() . Learning methods are :
  //
  // TMultiLayerPerceptron::kStochastic,
  // TMultiLayerPerceptron::kBatch,
  // TMultiLayerPerceptron::kSteepestDescent,
  // TMultiLayerPerceptron::kRibierePolak,
  // TMultiLayerPerceptron::kFletcherReeves,
  // TMultiLayerPerceptron::kBFGS
  
  //===========================================================================
  // Use TMLPAnalyzer to see what it looks for
  TCanvas* mlpa_canvas = new TCanvas("mlpa_canvas","Network Analysis",1600,300,800,900);
  mlpa_canvas->Divide(3,3);
  int ipad=1;
  mlpa_canvas->cd(ipad++);
  TLatex latex;
  char text[200];
  latex.SetTextSize(0.05);
  latex.SetTextAlign(13);  //align at top
  double ystep=0.07, ypos=0.95;
  #ifdef MMG_RUN
    sprintf(text,"MMG1-TRD Mode=%d",nn_mode);
  #else
    sprintf(text,"GEM-TRD Mode=%d",nn_mode);
  #endif
  latex.DrawLatex(0.05,ypos-=ystep,text);
  latex.DrawLatex(0.05,ypos-=ystep,rootfile);
  latex.DrawLatex(0.05,ypos-=ystep,NNcfg.data());
  
  //===========================================================================
  //=====            Train                                             ========
  //===========================================================================
  mlpa_canvas->cd(ipad++);
  #ifdef VERBOSE
    mlp->Train(ntrain, "text,graph,update=10,current");
  #else
    mlp->Train(ntrain, "graph,update=10,current");
  #endif
  mlp->Export("trd_ann","C++");
  
  //-------------------------------
  TMLPAnalyzer ana(mlp);
  ana.GatherInformations(); // Initialisation
  #ifdef VERBOSE
    ana.CheckNetwork(); // output to the terminal
  #endif
  mlpa_canvas->cd(ipad++);
  ana.DrawDInputs(); // shows how each variable influences the network
  //HistList->Add(ana.DrawDInputs());
  mlpa_canvas->cd(ipad++);
  mlp->Draw(); // shows the network structure
  //HistList->Add(mlp);
  mlpa_canvas->cd(ipad++);
  //gPad->WaitPrimitive(); /////////////////////////////////////////////
  ana.DrawNetwork(0,"type==1","type==0"); // draws the resulting network
  //HistList->Add(ana.DrawNetwork(0,"type==1","type==0"));
  // Use the NN to plot the results for each sample
  // This will give approx. the same result as DrawNetwork.
  // All entries are used, while DrawNetwork focuses on 
  // the test sample. Also the xaxis range is manually set.
  TH1F *bg = new TH1F("bgh", "NN bg output, single mod",115, -0.05, 1.1);
  TH1F *sig = new TH1F("sigh", "NN sig output, single mod",115, -0.05, 1.1);
  
  char htit[120]; sprintf(htit,"NN output, Nmod=%d",Nmod);
  TH1F *bgm = new TH1F("bgm",htit, 115, -.05, 1.1);
  TH1F *sigm = new TH1F("sigm",htit,115, -.05, 1.1);
  TH1F *hrejection_errors = new TH1F("hrejection_errors","Rej Factor Relative Error; Electron Purity %; Rejection Factor",6,62.5,92.5); hrejection_errors->SetStats(0);
  
  //---------------------------------------------------------------------
  //---------           test net                              -----------
  //---------------------------------------------------------------------
  
  double tout1=0.44, tout2=0.46;
  int GAUSS=0; double g_mean1=0.2, g_sigma1=0.3, g_mean2=0.8, g_sigma2=0.3;
  int DISP=0;
  if (DISP>0) {
    c2 = new TCanvas("Event Display","Event Display",1200,10,2000,1000); c2->SetRightMargin(0.15);
    c2->Divide(2,1);
  }
  
  TRandom3 *rndm = new TRandom3();
  Double_t params[NPAR];
  double sum=0;
  for (i=0; i<bg_tst->GetEntries(); i++) {
    bg_tst->GetEntry(i);
    for (int ip=0; ip<NPAR; ip++ ) params[ip] = Par[ip];
    double out = mlp->Evaluate(0, params);
    if (GAUSS) out=rndm->Gaus(g_mean1,g_sigma1);
    sum+=out;
    if (!((i+1)%Nmod)) {
      bgm->Fill(sum/Nmod);
      sum=0;
    }
    bg->Fill(out);
    
    if (tout1<out && out<tout2) { 
      #ifdef VERBOSE
        for (int ip=0; ip<NPAR; ip++) printf(" %f ",params[ip]); printf(" type = %d  out=%f iev=%d \n",type,out,ievent);
      #endif
    }
    if (out>0.7 && (DISP==1 || DISP==3 )) { 
      #ifdef VERBOSE
        printf(" pi high : out=%f type=%d iev=%d par=%5.1f %5.0f  %5.1f %5.1f %5.1f  \n",out,type,ievent,Par[0],Par[1],Par[2],Par[3],Par[4]);
      #endif
      ttree_hits->GetEntry(ievent);
      c2->cd(); disppi->Reset();
      char htit[128]; sprintf(htit,"#pi event %d; time 8 ns/bin ; x-strip ",ievent);   disppi->SetTitle(htit);
      #ifdef MMG_RUN
      for (int i=0;i<mmg1_nhit;i++) {
        #ifdef USE_CLUSTERS
          if (dedxc->at(i)>0) disppi->Fill(xposc->at(i),zposc->at(i),dedxc->at(i));
        #else
          if (dedx->at(i)>0) disppi->Fill(xpos->at(i),zpos->at(i),dedx->at(i));
        #endif
      }
      #else
      for (int i=0;i<gem_nhit;i++) {
        #ifdef USE_CLUSTERS
          if (dedxc->at(i)>0) disppi->Fill(xposc->at(i),zposc->at(i),dedxc->at(i));
        #else
	        if (dedx->at(i)>0) disppi->Fill(xpos->at(i),zpos->at(i),dedx->at(i));
        #endif
      }
      #endif
      if (disppi->GetEntries()!=0) disppi->Draw("colz");
      TLine lin1(0.,rtw1,300.,rtw1);   TLine lin2(0.,rtw3,300.,rtw3);  lin1.SetLineColor(kBlue); lin2.SetLineColor(kBlue); lin1.Draw(); lin2.Draw();  
      c2->Modified(); c2->Update();
    }
  } //-- End bg_tst entry loop
  
  sum=0;
  for (i=0; i<sig_tst->GetEntries(); i++) {
    sig_tst->GetEntry(i);
    for (int ip=0; ip<NPAR; ip++ ) params[ip] = Par[ip];
    double out = mlp->Evaluate(0, params);
    if (GAUSS) out=rndm->Gaus(g_mean2,g_sigma2);
    sum+=out;
    if (!((i+1)%Nmod)) {
      sigm->Fill(sum/Nmod);
      sum=0;
    }
    sig->Fill(out);
    
    if (tout1<out && out<tout2) {
      #ifdef VERBOSE
        for (int ip=0; ip<NPAR; ip++) printf(" %f ",params[ip]); printf(" type = %d out=%f iev=%d \n",type,out,ievent);
      #endif
    }
    
    //-------------------------------------------------------------
    if (out<0.1 && DISP>1 ) {
      #ifdef VERBOSE
        printf(" e low : out=%f type=%d iev=%d par=%5.1f %5.0f  %5.1f %5.1f %5.1f  \n",out,type,ievent,Par[0],Par[1],Par[2],Par[3],Par[4]);
      #endif
      ttree_hits->GetEntry(ievent);
      c2->cd(1);  dispe->Reset();
      char htit[128]; sprintf(htit,"electrons event %d ; time 8 ns/bin ; x-strip",ievent);  dispe->SetTitle(htit);
      int ii=0; double x[1000], y[1000], ey[1000];
      double yaver=0;
      #ifdef MMG_RUN
      for (int i=0;i<mmg1_nhit;i++) {
      #else
      for (int i=0;i<gem_nhit;i++) {
	    #endif
      #ifdef USE_CLUSTERS
        if (dedxc->at(i)>0) {
          dispe->Fill(zposc->at(i),xposc->at(i),dedxc->at(i));
          if (rtw1 < zposc->at(i) && zposc->at(i) < rtw3) {
            x[ii]=zposc->at(i);  y[ii]=xposc->at(i); if (dedxc->at(i)==0) ey[ii]=100; else ey[ii]=1000./dedxc->at(i); ii++; yaver+=y[ii];
          }
        }
      #else
        if (dedx->at(i)>0) {
	        dispe->Fill(zpos->at(i),xpos->at(i),dedx->at(i));
	        if (rtw1 < zpos->at(i) && zpos->at(i) < rtw3) {
	          x[ii]=zpos->at(i);  y[ii]=xpos->at(i); if (dedx->at(i)==0) ey[ii]=100; else ey[ii]=1000./dedx->at(i); ii++; yaver+=y[ii];
  	      }
	      }
      #endif
      }
      yaver/=ii;
      TF1 ffit1("ffit1", "pol1", rtw1, rtw3);  ffit1.SetLineColor(kBlue);
      TF1 ffit2("ffit2", "pol1", rtw1, rtw3);  ffit2.SetLineColor(kRed);
      
      TGraphErrors grr(ii, x, y, 0, ey);   grr.SetMinimum(150);   grr.SetMaximum(250);
      Double_t p1=0, p0=0;
      grr.Draw("ap");      
      if (ii>2) {
	      grr.Fit(&ffit1);
        grr.Fit(&ffit2, "+rob=0.75");
  	    p1=ffit2.GetParameter(1);
	      p0=ffit2.GetParameter(0);
	      #ifdef VERBOSE
          printf("p0=%f yaver=%f \n",p0,yaver);
        #endif
      }
      c2->cd(2);
      //gPad->WitPrimitive(); ///////////////
      if (dispe->GetEntries()!=0) dispe->Draw("colz");
      TLine lin1(rtw1,40.,rtw1,240.);   TLine lin2(rtw3,40.,rtw3,240.);  lin1.SetLineColor(kRed); lin2.SetLineColor(kRed); lin1.Draw(); lin2.Draw();
      ffit1.Draw("same");
      ffit2.Draw("same");
      c2->Modified(); c2->Update();
      c2->WaitPrimitive();
    } //-- End DISP>1
  } //-- End sig_tst Entry Loop
  
  //---------------------------------------------------------------------
  //---------               plot NN output                   ------------
  //---------------------------------------------------------------------
  mlpa_canvas->cd(ipad++);
  bg->SetLineColor(kBlue);
  bg->SetFillStyle(3345);   bg->SetFillColor(kBlue);
  sig->SetLineColor(kRed);
  sig->SetFillStyle(3354); sig->SetFillColor(kRed);
  bg->SetStats(0);
  sig->SetStats(0);
  bg->Draw();
  std::pair<double,double> rejLine90 = Reject(bg,sig,0.9);
  std::pair<double,double> rejLine70 = Reject(bg, sig, 0.7);
  TLine *line90 = new TLine(rejLine90.first,0,rejLine90.first,bg->GetMaximum()+1);
  TLine *line70 = new TLine(rejLine70.first,0,rejLine70.first,bg->GetMaximum()+1);
  line90->SetLineStyle(kDashed);   line70->SetLineStyle(kDashed);
  line90->Draw();  line70->Draw();
  sig->Draw("same");
  gPad->Modified(); gPad->Update();
  TLegend *legend = new TLegend(.7, .7, .9, .9);
  legend->AddEntry(bg, "Background (#pi)");
  legend->AddEntry(sig, "Signal (el)");
  legend->Draw();
  
  //---------  plot NN mod sum  ------------
  mlpa_canvas->cd(ipad++);
  bgm->SetLineColor(kBlue);
  bgm->SetFillStyle(3345);   bg->SetFillColor(kBlue);
  sigm->SetLineColor(kRed);
  sigm->SetFillStyle(3354); sig->SetFillColor(kRed);
  bgm->SetStats(0);
  sigm->SetStats(0);
  bgm->Draw();
  sigm->Draw("same");
  TLegend *legend2 = new TLegend(.7, .7, .9, .9);
  legend2->AddEntry(bgm, "Background (#pi)");
  legend2->AddEntry(sigm, "Signal (el)");
  legend2->Draw();
  //---------------------------------------------------------------------
  
  //---------- GRAPH 1 mod ------------------------
  const Int_t ngr = 9;
  const Int_t kNMAX = 9;
  Double_t *Xgr = new Double_t[kNMAX];
  Double_t *Ygr = new Double_t[kNMAX];
  cout << "\n=====================================================" << endl;
  cout << "======== BEGIN bg,sig Rejection Calculation =========\n" << endl;
  for ( int i=0; i<ngr; i++) {
    double eeff=(i+1)*0.1;
    std::pair<double,double> rej = Reject(bg, sig, eeff);
    Xgr[i] = eeff;
    Ygr[i] = rej.first;
    cout << " i=" << i << " x=" << Xgr[i] << " y=" << Ygr[i] << endl;
    cout << "---------------------------------------- \n" << endl;
  }
  cout << "\n======== END bg,sig Rejection Calculation =========" << endl;
  cout << "===================================================\n" << endl;
  TGraph *gr = new TGraph(ngr,Xgr,Ygr); gr->SetName("e #pi efficiency");gr->SetTitle("Efficiency single module");
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr,"lp");
  mg->SetTitle("Efficiency single module");
  mg->GetXaxis()->SetTitle("El Efficiency"); mg->GetYaxis()->SetTitle("#pi Efficiency");
  mlpa_canvas->cd(ipad++);
  mg->Draw("a");
  
  //---------------------------------------------------------------------
  mlpa_canvas->cd(ipad++);   gPad->SetLogy(); hcount->Draw();
  //gPad->WaitPrimitive();
  //---------------------------------------------------------------------
  std::pair<double,double> rej70 = Reject(bg, sig, 0.7);
  cout << "---------------------------------------- \n" << endl;
  std::pair<double,double> rej80 = Reject(bg, sig, 0.8);
  cout << "---------------------------------------- \n" << endl;
  std::pair<double,double> rej85 = Reject(bg, sig, 0.85);
  cout << "---------------------------------------- \n" << endl;
  std::pair<double,double> rej90 = Reject(bg, sig, 0.9);
  cout << "---------------------------------------- \n" << endl;
  
  double rejectionValues[4] = {rej70.first, rej80.first, rej85.first, rej90.first};
  double rejErrors[4] = {rej70.second, rej80.second, rej85.second, rej90.second};
  for (int i=1; i<=6; i++) {
    if (i==2) {
      hrejection_errors->SetBinContent(i, 1./rejectionValues[i-2]);
      hrejection_errors->SetBinError(i, rejErrors[i-2]*hrejection_errors->GetBinContent(i));
    }
    if (i>3) {
      hrejection_errors->SetBinContent(i, 1./rejectionValues[i-3]);
      hrejection_errors->SetBinError(i, rejErrors[i-3]*hrejection_errors->GetBinContent(i));
    }
  }
  
  #ifdef VERBOSE
    cout << " 1-Mod e=70% , Eff pi = " << rej70.first*100. << "% ,  Rejection =" << 1./rej70.first << endl;
    cout << " 1-Mod pi=90%, Eff pi = " << rej90.first*100. << "% ,  Rejection =" << 1./rej90.first << endl;
  #endif
  //---------------------------------------------
  #ifdef VERBOSE
  std::pair<double,double> rej70m = Reject(bgm, sigm, 0.7);
  cout << "---------------------------------------- \n" << endl;
  std::pair<double,double> rej90m = Reject(bgm, sigm, 0.9);
  cout << "---------------------------------------- \n" << endl;
    cout << "Nmod=" << Nmod << " Nmod=" << Nmod << " e=70% , Eff pi = " << rej70m.first*100. << "% ,  Rejection =" << 1./rej70m.first << endl;
    cout << " Nmod=" << Nmod << " pi=90%, Eff pi = " << rej90m.first*100. << "% ,  Rejection =" << 1./rej90m.first << endl;
  #endif
  //---------------------------------------------
  mlpa_canvas->cd(1);
  stringstream ss;
  ss << "Nmod=" << 1 << ", e=70%, #pi=" << std::fixed << std::setprecision(2) << rej70.first*100. << "%, Rej=" << 1./rej70.first ;  string str = ss.str();
  latex.DrawLatex(0.05,ypos-=ystep,str.data());
  //--
  ss.str("");  ss.clear();
  ss << "Nmod=" << 1 << ", e=80%, Eff #pi=" << rej80.first*100. << "%, Rej=" << 1./rej80.first ;  string str0 = ss.str();
  latex.DrawLatex(0.05,ypos-=ystep,str0.data());
  //--
  ss.str("");  ss.clear(); 
  ss << "Nmod=" << 1 << ", e=85%, Eff #pi=" << rej85.first*100. << "%, Rej=" << 1./rej85.first ;  string str1 = ss.str();
  latex.DrawLatex(0.05,ypos-=ystep,str1.data());
  //--
  ss.str("");  ss.clear(); 
  ss << "Nmod=" << 1 << ", e=90%, Eff #pi=" << rej90.first*100. << "%, Rej=" << 1./rej90.first ;  string str2 = ss.str();
  latex.DrawLatex(0.05,ypos-=ystep,str2.data());
  latex.DrawLatex(0.05,ypos-=ystep,"--------------");
  //--
  #ifdef VERBOSE
  ss.str("");  ss.clear(); 
  ss << " Nmod=" << Nmod << " e=70% , Eff #pi = " << rej70m.first*100. << "% ,  Rej =" << 1./rej70m.first ;  str2 = ss.str();
  latex.DrawLatex(0.05,ypos-=ystep,str2.data());
  //--
  ss.str("");  ss.clear(); 
  ss << " Nmod=" << Nmod << " e=90% , Eff #pi = " << rej90m.first*100. << "% ,  Rej =" << 1./rej90m.first ;  str2 = ss.str();
  latex.DrawLatex(0.05,ypos-=ystep,str2.data());
  #endif
  //---------------------------------------------
  #ifdef MMG_RUN
    sprintf(text,"mlpOutput/%s_mmg_m%d.png",basename,nn_mode);
  #else
    sprintf(text,"mlpOutput/%s_m%d.png",basename,nn_mode);
  #endif
  //mlpa_canvas->WaitPrimitive();
  mlpa_canvas->Print(text);
  mlpa_canvas->cd(0);
  //mlpa_canvas->WaitPrimitive();
  
  //HistList->Write("HistDQM", TObject::kSingleKey);
  outputFile->Write();
  delete inputFile;
  cout << "============== END RUN " << RunNum << "===============" <<endl;
}


//==================================================================
void Count(const char *tit) {
  hcount->Fill(tit,1);
}
void Count(const char *tit, double cut1) {
  char clab[20];
  sprintf(clab,"%s_%.1f",tit,cut1);
  hcount->Fill(clab,1);
}
void Count(const char *tit, double cut1, double cut2) {
  char clab[20];
  sprintf(clab,"%s_%.1f_%.1f",tit,cut1,cut2);
  hcount->Fill(clab,1);
}

//------------------------------------------------------------------
