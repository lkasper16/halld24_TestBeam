//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jun 17 22:57:15 2023 by ROOT version 6.14/04
// from TTree gem_hits/GEM TTree with single track hit info
// found on file: trd_singleTrackHits_Run_003200.root
//////////////////////////////////////////////////////////


#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"


// Declaration of leaf types
Int_t           event_num;
Int_t           gem_nhit;
Int_t           gem_nclu;

//vector<bool>    *parID;
vector<int>     *xpos;
vector<float>   *zpos;
vector<float>   *dedx;
vector<float>   *zHist;
vector<float>   *xposc;
vector<float>   *zposc;
vector<float>   *dedxc;
vector<float>   *widthc;

float     xposc_max;
float     zposc_max;
float     dedxc_max;
float     widthc_max;

// List of branches
TBranch        *b_event_num;   //!
TBranch        *b_gem_nhit;   //!
TBranch        *b_gem_nclu;
TBranch        *b_xpos;   //!
TBranch        *b_zpos;   //!
TBranch        *b_dedx;   //!
//TBranch        *b_parID;   //!
TBranch        *b_zHist;   //!
TBranch        *b_xposc;   //!
TBranch        *b_zposc;   //!
TBranch        *b_dedxc;   //!
TBranch        *b_widthc;   //!
TBranch        *b_xposc_max;
TBranch        *b_zposc_max;
TBranch        *b_dedxc_max;
TBranch        *b_widthc_max;
