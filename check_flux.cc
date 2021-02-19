#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

//#include<sys/stat.h>

#include<map>
#include<set>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"

#include "TVector3.h"
#include "TRotation.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"

#include "TRandom3.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"

#include "TPrincipal.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TPad.h"
#include "TLegend.h"
#include "TString.h"
#include "TColor.h"

#include "draw.icc"

void check_flux()
{
  //////////////////////////////////////////////////////////////////////////////////////// Draw style

  gStyle->SetOptStat(0);
  //gStyle->SetPalette(kBird);

  double snWidth = 2;

  // use medium bold lines and thick markers
  gStyle->SetLineWidth(snWidth);
  gStyle->SetFrameLineWidth(snWidth);
  gStyle->SetHistLineWidth(snWidth);
  gStyle->SetFuncWidth(snWidth);
  gStyle->SetGridWidth(snWidth);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.0);
  gStyle->SetEndErrorSize(4);
  gStyle->SetEndErrorSize(0);

  ////////////////////////////////////////////////////////////////////////////////////////
  
  TString roostr = "";

  TFile *roofile = new TFile("./new_checkout/checkout_prodgenie_bnb_nu_overlay_run1.root", "read");
  TTree *tree = (TTree*)roofile->Get("wcpselection/T_PFeval");
  int entries = tree->GetEntries();
  cout<<" ---> entries "<<entries<<endl;
  
  Float_t         mcflux_genx;
  Float_t         mcflux_geny;
  Float_t         mcflux_genz;
  Float_t         truth_nu_pos[4];
  Float_t         truth_nu_momentum[4];
  TBranch        *b_mcflux_genx;   //!
  TBranch        *b_mcflux_geny;   //!
  TBranch        *b_mcflux_genz;   //!
  TBranch        *b_truth_nu_pos;   //!
  TBranch        *b_truth_nu_momentum;   //!

  tree->SetBranchAddress("mcflux_genx", &mcflux_genx, &b_mcflux_genx);
  tree->SetBranchAddress("mcflux_geny", &mcflux_geny, &b_mcflux_geny);
  tree->SetBranchAddress("mcflux_genz", &mcflux_genz, &b_mcflux_genz);  
  tree->SetBranchAddress("truth_nu_pos", truth_nu_pos, &b_truth_nu_pos);
  tree->SetBranchAddress("truth_nu_momentum", truth_nu_momentum, &b_truth_nu_momentum);

  TH1D *h1_diff_direction = new TH1D("h1_diff_direction", "", 100, -0.01, 0.01);
  
  
  for(int ientry=0; ientry<entries; ientry++) {
    tree->GetEntry( ientry );

    TVector3 vtx_window( mcflux_genx*100, mcflux_geny*100, mcflux_genz*100 );
    TVector3 vtx_det( truth_nu_pos[0], truth_nu_pos[1], truth_nu_pos[2] );
    TVector3 vtx_direction = vtx_det - vtx_window;
    TVector3 vtx_direction_unit = vtx_direction * (1./vtx_direction.Mag());

    TVector3 momentum_( truth_nu_momentum[0], truth_nu_momentum[1], truth_nu_momentum[2] );
    TVector3 momentum_direction_unit = momentum_ * (1./truth_nu_momentum[3]);

    TVector3 vc_diff = vtx_direction_unit - momentum_direction_unit;       
    
    cout<<TString::Format(" ---> %4d, x/y/z %13.10f %13.10f %13.10f",
			  ientry,
			  vc_diff.X(),
			  vc_diff.Y(),
			  vc_diff.Z()
			  )<<endl;

    h1_diff_direction->Fill( vc_diff.X() );
    h1_diff_direction->Fill( vc_diff.Y() );
    h1_diff_direction->Fill( vc_diff.Z() );
  }

  cout<<" ---> Integral "<<h1_diff_direction->Integral()<<endl;
  
  TCanvas *canv_h1_diff_direction = new TCanvas("canv_h1_diff_direction", "canv_h1_diff_direction", 900, 650);
  func_canv_margin(canv_h1_diff_direction, 0.15, 0.1, 0.1, 0.15);
  h1_diff_direction->Draw("hist");
  h1_diff_direction->SetLineColor(kBlue);
  func_title_size(h1_diff_direction, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h1_diff_direction, "Difference in x-dir, y-dir, z-dir", "Entries");
  h1_diff_direction->GetXaxis()->CenterTitle();
  h1_diff_direction->GetYaxis()->CenterTitle();
  h1_diff_direction->GetXaxis()->SetNdivisions(506);
  h1_diff_direction->GetXaxis()->SetTitleOffset(1.2);
  h1_diff_direction->Draw("same axis");
  canv_h1_diff_direction->SaveAs("canv_h1_diff_direction.png");
}
