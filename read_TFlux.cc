#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<sys/stat.h>

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

////////////////////////////////////////////////////////////////////////////////////////// TFlux

class TFlux {
public:

  ////////////////////////// constructor
  
  TFlux() {
    vc_fluxfile_list.clear();
    val_POT = 0;

    // DocDB-13804, 15522
    vc_uB_LArTPC_halfXYZ.SetX( (254.8+1.55)/2 );// cm
    vc_uB_LArTPC_halfXYZ.SetY( (117.47+115.53)/2 );
    vc_uB_LArTPC_halfXYZ.SetZ( (1036.9-0.1)/2 );
    
    vc_uB_LArTPC_center.SetX( -1.55   + vc_uB_LArTPC_halfXYZ.X() );
    vc_uB_LArTPC_center.SetY( -115.53 + vc_uB_LArTPC_halfXYZ.Y() );
    vc_uB_LArTPC_center.SetZ( -0.1    + vc_uB_LArTPC_halfXYZ.Z() );
  }
  
  ////////////////////////// function member

  void Exe_fluxfile(int bgn, int end);
  bool Flag_file_exist(TString filename);

  void Exe_read_tree(int out_index);

  bool Flag_pass_XYplane_Zfixed(double Zfixed, TVector3 point_at_window, TVector3 direction_unit);  
  bool Flag_pass_Circle_fixed(TVector3 point_fixed, double radius, TVector3 point_at_window, TVector3 direction_unit);  
  bool Flag_pass_Square_fixed(TVector3 point_fixed, double radius, TVector3 point_at_window, TVector3 direction_unit);
  bool Flag_pass_Polygon_Zfixed(vector<double> vtx_XY_fixed, TVector3 point_at_window, TVector3 direction_unit);
  
  bool Flag_pass_Cuboid_fixed(TVector3 center_fixed, TVector3 halfXYZ, TVector3 point_at_window, TVector3 direction_unit);

  
  ////////////////////////// data member

  vector<TString>vc_fluxfile_list;

  double val_POT;

  TVector3 vc_uB_LArTPC_center;
  TVector3 vc_uB_LArTPC_halfXYZ;
 
};

//////////////////////////// ccc

//     Y            X
//     |           -
//     |         -
//     |       -
//     |     -
//     |   -
//     | -
//     |------------Z


bool TFlux::Flag_pass_Cuboid_fixed(TVector3 center_fixed, TVector3 halfXYZ, TVector3 point_at_window, TVector3 direction_unit)
{
  bool flag = false;
  int check_case = 0;
  
  // principle: if there is >=one intersection point, then the ray (outside the cuboid) passes througth the cuboid

  double plane_XY_low = center_fixed.Z() - halfXYZ.Z();
  double plane_XY_hgh = center_fixed.Z() + halfXYZ.Z();

  double plane_XZ_low = center_fixed.Y() - halfXYZ.Y();
  double plane_XZ_hgh = center_fixed.Y() + halfXYZ.Y();

  double plane_YZ_low = center_fixed.X() - halfXYZ.X();
  double plane_YZ_hgh = center_fixed.X() + halfXYZ.X();

  /// case plane_XY_low
  {
    double step_length = ( plane_XY_low - point_at_window.Z() )/direction_unit.Z();
    if( step_length>0 ) {
      double intersection_X = point_at_window.X() + step_length*direction_unit.X();
      double intersection_Y = point_at_window.Y() + step_length*direction_unit.Y();
      double intersection_Z = point_at_window.Z() + step_length*direction_unit.Z();
      TVector3 point_intersection(intersection_X, intersection_Y, intersection_Z);
      //TVector3 point_center(center_fixed.X(), center_fixed.Y(), center_fixed.Z());
      TVector3 point_center(center_fixed.X(), center_fixed.Y(), plane_XY_low);
      TVector3 vc_diff = point_intersection - point_center;
      if( fabs(vc_diff.X())<=halfXYZ.X() && fabs(vc_diff.Y())<=halfXYZ.Y() ) {
	flag = true;
	check_case = 1;
	goto comehereAA;
      }                  
    }
  }
    
  /// case plane_XZ_low
  {
    double step_length = ( plane_XZ_low - point_at_window.Y() )/direction_unit.Y();
    
    double intersection_X = point_at_window.X() + step_length*direction_unit.X();
    double intersection_Y = point_at_window.Y() + step_length*direction_unit.Y();
    double intersection_Z = point_at_window.Z() + step_length*direction_unit.Z();
    TVector3 point_intersection(intersection_X, intersection_Y, intersection_Z);
    //TVector3 point_center(center_fixed.X(), center_fixed.Y(), center_fixed.Z());
    TVector3 point_center(center_fixed.X(), plane_XZ_low, center_fixed.Z());
    TVector3 vc_diff = point_intersection - point_center;
    if( fabs(vc_diff.X())<=halfXYZ.X() && fabs(vc_diff.Z())<=halfXYZ.Z() ) {
      flag = true;
      check_case = 3;
      goto comehereAA;
    }                  
  }
 
  /// case plane_XZ_hgh
  {
    double step_length = ( plane_XZ_hgh - point_at_window.Y() )/direction_unit.Y();
    
    double intersection_X = point_at_window.X() + step_length*direction_unit.X();
    double intersection_Y = point_at_window.Y() + step_length*direction_unit.Y();
    double intersection_Z = point_at_window.Z() + step_length*direction_unit.Z();
    TVector3 point_intersection(intersection_X, intersection_Y, intersection_Z);
    //TVector3 point_center(center_fixed.X(), center_fixed.Y(), center_fixed.Z());
    TVector3 point_center(center_fixed.X(), plane_XZ_hgh, center_fixed.Z());
    TVector3 vc_diff = point_intersection - point_center;
    if( fabs(vc_diff.X())<=halfXYZ.X() && fabs(vc_diff.Z())<=halfXYZ.Z() ) {
      flag = true;
      check_case = 4;
      goto comehereAA;
    }                  
  }
    
  /// case plane_YZ_low
  {
    double step_length = ( plane_YZ_low - point_at_window.X() )/direction_unit.X();
    
    double intersection_X = point_at_window.X() + step_length*direction_unit.X();
    double intersection_Y = point_at_window.Y() + step_length*direction_unit.Y();
    double intersection_Z = point_at_window.Z() + step_length*direction_unit.Z();
    TVector3 point_intersection(intersection_X, intersection_Y, intersection_Z);
    //TVector3 point_center(center_fixed.X(), center_fixed.Y(), center_fixed.Z());
    TVector3 point_center(plane_YZ_low, center_fixed.Y(), center_fixed.Z());
    TVector3 vc_diff = point_intersection - point_center;
    if( fabs(vc_diff.Y())<=halfXYZ.Y() && fabs(vc_diff.Z())<=halfXYZ.Z() ) {
      flag = true;
      check_case = 5;
      goto comehereAA;
    }                  
  }
        
  /// case plane_YZ_hgh
  {
    double step_length = ( plane_YZ_hgh - point_at_window.X() )/direction_unit.X();
    
    double intersection_X = point_at_window.X() + step_length*direction_unit.X();
    double intersection_Y = point_at_window.Y() + step_length*direction_unit.Y();
    double intersection_Z = point_at_window.Z() + step_length*direction_unit.Z();
    TVector3 point_intersection(intersection_X, intersection_Y, intersection_Z);
    //TVector3 point_center(center_fixed.X(), center_fixed.Y(), center_fixed.Z());
    TVector3 point_center(plane_YZ_hgh, center_fixed.Y(), center_fixed.Z());
    TVector3 vc_diff = point_intersection - point_center;
    if( fabs(vc_diff.Y())<=halfXYZ.Y() && fabs(vc_diff.Z())<=halfXYZ.Z() ) {
      flag = true;
      check_case = 6;
      goto comehereAA;
    }                  
  }
    
 comehereAA:
  //if( check_case!=0 ) cout<<" - "<<check_case<<endl;
  //if( !((check_case==0) || (check_case==1)) ) cout<<" - "<<check_case<<endl;  
  
  return flag;
}

bool TFlux::Flag_pass_Square_fixed(TVector3 point_fixed, double radius, TVector3 point_at_window, TVector3 direction_unit)
{
  bool flag = false;
  
  double step_length = ( point_fixed.Z() - point_at_window.Z() )/direction_unit.Z();
  
  double intersection_X = point_at_window.X() + step_length*direction_unit.X();
  double intersection_Y = point_at_window.Y() + step_length*direction_unit.Y();
  double intersection_Z = point_fixed.Z();

  TVector3 point_intersection(intersection_X, intersection_Y, intersection_Z);
  TVector3 vc_diff = point_intersection - point_fixed;
  
  if(  fabs(vc_diff.X())<=radius && fabs(vc_diff.Y())<=radius  ) {
    flag = true;
  }
  
  return flag;
}

bool TFlux::Flag_pass_Circle_fixed(TVector3 point_fixed, double radius, TVector3 point_at_window, TVector3 direction_unit)
{
  bool flag = false;
  
  double step_length = ( point_fixed.Z() - point_at_window.Z() )/direction_unit.Z();
  
  double intersection_X = point_at_window.X() + step_length*direction_unit.X();
  double intersection_Y = point_at_window.Y() + step_length*direction_unit.Y();
  double intersection_Z = point_fixed.Z();
  
  TVector3 point_intersection(intersection_X, intersection_Y, intersection_Z);
  
  if( (point_intersection-point_fixed).Mag() <= radius ) {
    flag = true;
  }
  
  return flag;
}

bool TFlux::Flag_pass_XYplane_Zfixed(double Zfixed, TVector3 point_at_window, TVector3 direction_unit)
{
  bool flag = false;

  double step_length = ( Zfixed - point_at_window.Z() )/direction_unit.Z();
  // direction: step_length >0 or <0
  
  double intersection_X = point_at_window.X() + step_length*direction_unit.X();
  double intersection_Y = point_at_window.Y() + step_length*direction_unit.Y();
  
  double diff_X = fabs( intersection_X - vc_uB_LArTPC_center.X() );
  double diff_Y = fabs( intersection_Y - vc_uB_LArTPC_center.Y() );
  
  if( diff_X<=vc_uB_LArTPC_halfXYZ.X() && diff_Y<=vc_uB_LArTPC_halfXYZ.Y() ) {
    flag = true;      
  }    
  
  return flag;
}



//////////////////////////// ccc

void TFlux::Exe_fluxfile(int bgn, int end)
{
  cout<<endl<<" ---> read flux file(s)"<<endl;

  for(int idx=bgn; idx<=end; idx++) {
    TString roostr = TString::Format("/home/xji/data0/work/101_flux_Xs/bnb_gsimple_fluxes_01.09.2019_463/converted_beammc_wincorr_%04d.root", idx);
    //TString roostr = TString::Format("/pnfs/uboone/persistent/uboonebeam/bnb_gsimple/bnb_gsimple_fluxes_01.09.2019_463/converted_beammc_wincorr_%04d.root", idx);
    
    if( Flag_file_exist(roostr) ) {
      vc_fluxfile_list.push_back( roostr );
      cout<<" - "<<roostr<<endl;
    }
  }

  cout<<TString::Format("      input file(s) %d", (int)(vc_fluxfile_list.size()) )<<endl<<endl;
}

bool TFlux::Flag_file_exist(TString filename)
{
  bool flag = false;
  struct stat st;
  if( stat(filename, &st)==0 ) flag = true;
  return flag;    
}

//////////////////////////// ccc

void TFlux::Exe_read_tree(int out_index)
{
  cout<<" ---> read tree "<<endl<<endl;
  
  cout<<" uB LArTPC centerX and halfX, "<<vc_uB_LArTPC_center.X()<<"\t"<<vc_uB_LArTPC_halfXYZ.X()<<endl;
  cout<<" uB LArTPC centerY and halfY, "<<vc_uB_LArTPC_center.Y()<<"\t"<<vc_uB_LArTPC_halfXYZ.Y()<<endl;
  cout<<" uB LArTPC centerZ and halfZ, "<<vc_uB_LArTPC_center.Z()<<"\t"<<vc_uB_LArTPC_halfXYZ.Z()<<endl<<endl;
  
  TChain *tree_flux = new TChain("flux");
  TChain *tree_meta = new TChain("meta");

  for(vector<TString>::iterator it_vc = vc_fluxfile_list.begin(); it_vc!=vc_fluxfile_list.end(); ++it_vc) {
    TString filename = *it_vc;
    tree_flux->Add(filename);
    tree_meta->Add(filename);
  }

  //////////////////////////////////////
  
  Double_t        protons;
  TBranch        *b_meta_protons;   //!
  tree_meta->SetBranchAddress("protons", &protons, &b_meta_protons);
  
  double total_pot = 0;
  int entries_meta = tree_meta->GetEntries();
  cout<<endl<<" ---> entries_meta "<<entries_meta<<endl;
  for(int ientry=0; ientry<entries_meta; ientry++) {
    tree_meta->GetEntry( ientry );
    protons = b_meta_protons->GetLeaf("protons")->GetValue();
    total_pot += protons;
  }
  cout<<" ---> total POT "<<total_pot<<endl<<endl;

  val_POT = total_pot;
  
  //////////////////////////////////////

  const int numu = 14;
  
  // Decay mode that produced neutrino: 
  // 1  K0L -> nue pi- e+
  // 2  K0L -> nuebar pi+ e-
  // 3  K0L -> numu pi- mu+
  // 4  K0L -> numubar pi+ mu-
  // 5  K+  -> numu mu+
  // 6  K+  -> nue pi0 e+
  // 7  K+  -> numu pi0 mu+
  // 8  K-  -> numubar mu-
  // 9  K-  -> nuebar pi0 e-
  // 10  K-  -> numubar pi0 mu-
  // 11  mu+ -> numubar nue e+
  // 12  mu- -> numu nuebar e-
  // 13  pi+ -> numu mu+
  // 14  pi- -> numubar mu-
	
  // Declaration of leaf types
  // genie::flux::GSimpleNtpEntry *entry;
  Double_t        wgt;
  Double_t        vtxx;// vertex
  Double_t        vtxy;
  Double_t        vtxz;
  Double_t        dist;
  Double_t        px;
  Double_t        py;
  Double_t        pz;
  Double_t        E;
  Int_t           pdg;// nue 12, numu 14
  UInt_t          metakey;
  // genie::flux::GSimpleNtpNuMI *numi;
  Double_t        tpx;// parent momentum exiting the target
  Double_t        tpy;
  Double_t        tpz;
  Double_t        vx;
  Double_t        vy;
  Double_t        vz;
  Double_t        pdpx;// parent momentum at decay point
  Double_t        pdpy;
  Double_t        pdpz;
  Double_t        pppx;
  Double_t        pppy;
  Double_t        pppz;
  Int_t           ndecay;
  Int_t           ptype;// Parent GEANT code particle ID 
  Int_t           ppmedium;
  Int_t           tptype;// Parent particle ID exiting the target (GEANT code) 
  Int_t           run;
  Int_t           evtno;
  Int_t           entryno;
  
  // List of branches
  TBranch        *b_entry_wgt;   //!
  TBranch        *b_entry_vtxx;   //!
  TBranch        *b_entry_vtxy;   //!
  TBranch        *b_entry_vtxz;   //!
  TBranch        *b_entry_dist;   //!
  TBranch        *b_entry_px;   //!
  TBranch        *b_entry_py;   //!
  TBranch        *b_entry_pz;   //!
  TBranch        *b_entry_E;   //!
  TBranch        *b_entry_pdg;   //!
  TBranch        *b_entry_metakey;   //!
  TBranch        *b_numi_tpx;   //!
  TBranch        *b_numi_tpy;   //!
  TBranch        *b_numi_tpz;   //!
  TBranch        *b_numi_vx;   //!
  TBranch        *b_numi_vy;   //!
  TBranch        *b_numi_vz;   //!
  TBranch        *b_numi_pdpx;   //!
  TBranch        *b_numi_pdpy;   //!
  TBranch        *b_numi_pdpz;   //!
  TBranch        *b_numi_pppx;   //!
  TBranch        *b_numi_pppy;   //!
  TBranch        *b_numi_pppz;   //!
  TBranch        *b_numi_ndecay;   //!
  TBranch        *b_numi_ptype;   //!
  TBranch        *b_numi_ppmedium;   //!
  TBranch        *b_numi_tptype;   //!
  TBranch        *b_numi_run;   //!
  TBranch        *b_numi_evtno;   //!
  TBranch        *b_numi_entryno;   //!
  
  // Set branch addresses and branch pointers
  tree_flux->SetBranchAddress("wgt", &wgt, &b_entry_wgt);
  tree_flux->SetBranchAddress("vtxx", &vtxx, &b_entry_vtxx);
  tree_flux->SetBranchAddress("vtxy", &vtxy, &b_entry_vtxy);
  tree_flux->SetBranchAddress("vtxz", &vtxz, &b_entry_vtxz);
  tree_flux->SetBranchAddress("dist", &dist, &b_entry_dist);
  tree_flux->SetBranchAddress("px", &px, &b_entry_px);
  tree_flux->SetBranchAddress("py", &py, &b_entry_py);
  tree_flux->SetBranchAddress("pz", &pz, &b_entry_pz);
  tree_flux->SetBranchAddress("E", &E, &b_entry_E);
  tree_flux->SetBranchAddress("pdg", &pdg, &b_entry_pdg);
  tree_flux->SetBranchAddress("metakey", &metakey, &b_entry_metakey);
  tree_flux->SetBranchAddress("tpx", &tpx, &b_numi_tpx);
  tree_flux->SetBranchAddress("tpy", &tpy, &b_numi_tpy);
  tree_flux->SetBranchAddress("tpz", &tpz, &b_numi_tpz);
  tree_flux->SetBranchAddress("vx", &vx, &b_numi_vx);
  tree_flux->SetBranchAddress("vy", &vy, &b_numi_vy);
  tree_flux->SetBranchAddress("vz", &vz, &b_numi_vz);
  tree_flux->SetBranchAddress("pdpx", &pdpx, &b_numi_pdpx);
  tree_flux->SetBranchAddress("pdpy", &pdpy, &b_numi_pdpy);
  tree_flux->SetBranchAddress("pdpz", &pdpz, &b_numi_pdpz);
  tree_flux->SetBranchAddress("pppx", &pppx, &b_numi_pppx);
  tree_flux->SetBranchAddress("pppy", &pppy, &b_numi_pppy);
  tree_flux->SetBranchAddress("pppz", &pppz, &b_numi_pppz);
  tree_flux->SetBranchAddress("ndecay", &ndecay, &b_numi_ndecay);
  tree_flux->SetBranchAddress("ptype", &ptype, &b_numi_ptype);
  tree_flux->SetBranchAddress("ppmedium", &ppmedium, &b_numi_ppmedium);
  tree_flux->SetBranchAddress("tptype", &tptype, &b_numi_tptype);
  tree_flux->SetBranchAddress("run", &run, &b_numi_run);
  tree_flux->SetBranchAddress("evtno", &evtno, &b_numi_evtno);
  tree_flux->SetBranchAddress("entryno", &entryno, &b_numi_entryno);

  int entries = tree_flux->GetEntries();

  ////////////////

  cout<<Form(" ---> entries_flux_event %d", entries)<<endl<<endl;

  ////// from Zarko
  double beam2det_x0 = 1.24325;
  double beam2det_y0 = -0.0093;
  double beam2det_z0 = -463.363525;
  TVector3 beam2det(beam2det_x0, beam2det_y0, beam2det_z0);
  TRotation rr;
  rr.RotateX(0.016/10.712);
  rr.RotateY(0.036/10.712);

  //////
  double x_min = 1e6; double x_max = -1e6;
  double y_min = 1e6; double y_max = -1e6;

  ////////////////////////////////////////

  TH1D *h1_POT = new TH1D("h1_POT", "", 1, 0, 1); h1_POT->SetBinContent(1, val_POT);

  TH1D *h1_numu_window = new TH1D("h1_numu_window", "", 100, 0, 5);
  TH1D *h1_numu_XYplane = new TH1D("h1_numu_XYplane", "", 100, 0, 5);
  TH1D *h1_numu_Circle = new TH1D("h1_numu_Circle", "", 100, 0, 5);
  TH1D *h1_numu_Square = new TH1D("h1_numu_Square", "", 100, 0, 5);
  TH1D *h1_numu_Cuboid = new TH1D("h1_numu_Cuboid", "", 100, 0, 5);
  
  ////////////////////////////////////////
  
  for(int ientry=0; ientry<entries; ientry++) {
    if( ientry%max(1, entries/10)==0 ) cout<<Form(" ---> processing %6.3f, %12d", ientry*1./entries, ientry)<<endl;
    else if( ientry==entries-1 ) cout<<Form(" ---> processing %6.3f, %12d", ientry*1./entries, ientry)<<endl;

    tree_flux->GetEntry(ientry);

    //////////////////
    
    vtxx = b_entry_vtxx->GetLeaf("vtxx")->GetValue() *100;// m ---> cm, neutrino at window (at detector coordinate)
    vtxy = b_entry_vtxy->GetLeaf("vtxy")->GetValue() *100;
    vtxz = b_entry_vtxz->GetLeaf("vtxz")->GetValue() *100;

    px = b_entry_px->GetLeaf("px")->GetValue();
    py = b_entry_py->GetLeaf("py")->GetValue();
    pz = b_entry_pz->GetLeaf("pz")->GetValue();

    dist = b_entry_dist->GetLeaf("dist")->GetValue();

    pdg = b_entry_pdg->GetLeaf("pdg")->GetValue();
    E = b_entry_E->GetLeaf("E")->GetValue();// GeV
    
    vx = b_numi_vx->GetLeaf("vx")->GetValue();// cm, neutrino at generated position (at beam coordinate)
    vy = b_numi_vy->GetLeaf("vy")->GetValue();
    vz = b_numi_vz->GetLeaf("vz")->GetValue();
        
    TVector3 vtx_nu_at_window(vtxx, vtxy, vtxz);

    TVector3 direction_unit(px/E, py/E, pz/E);
    double val_phi = direction_unit.Phi() * 180./3.1415926;
    double val_theta = direction_unit.Theta() * 180./3.1415926;
    
    //////////////////

    TVector3 point_fixed = vc_uB_LArTPC_center;
    double radius_fixed = 50;
    TVector3 halfXYZ_fixed = vc_uB_LArTPC_halfXYZ;
    
    bool flag_pass_XYplane = Flag_pass_XYplane_Zfixed(point_fixed.Z(), vtx_nu_at_window, direction_unit);
    bool flag_pass_Circle  = Flag_pass_Circle_fixed(point_fixed, radius_fixed, vtx_nu_at_window, direction_unit);
    bool flag_pass_Square  = Flag_pass_Square_fixed(point_fixed, radius_fixed, vtx_nu_at_window, direction_unit);
    bool flag_pass_Cuboid  = Flag_pass_Cuboid_fixed(vc_uB_LArTPC_center, vc_uB_LArTPC_halfXYZ, vtx_nu_at_window, direction_unit);

    if( pdg==numu ) {

      h1_numu_window->Fill( E );
      
      if( flag_pass_XYplane ) h1_numu_XYplane->Fill( E );
      if( flag_pass_Circle )  h1_numu_Circle->Fill( E );
      if( flag_pass_Square )  h1_numu_Square->Fill( E );
      if( flag_pass_Cuboid )  h1_numu_Cuboid->Fill( E );
    }
    
  }  
  cout<<endl;

  TFile *outfile = new TFile(TString::Format("subfile_%06d.root", out_index), "recreate");
  h1_POT->Write();
  h1_numu_window->Write();
  h1_numu_XYplane->Write();
  h1_numu_Circle->Write();
  h1_numu_Square->Write();
  h1_numu_Cuboid->Write();
  outfile->Close();
  
  ///////////////////////////////// plotting
  
  double area_XYplane = vc_uB_LArTPC_halfXYZ.X() * vc_uB_LArTPC_halfXYZ.Y() * 4;
  double area_Circle = 3.1415926 * 50*50;
  double area_Square = 50*50 *4;
  double area_Cuboid = area_XYplane;
  double bin_width = h1_numu_XYplane->GetBinWidth(1);

}

//////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{  

  int infile_bgn = 0;
  int infile_end = 0;
  int outfile_idx = 1;

  for(int i=1; i<argc; i++) {    
    if( strcmp(argv[i],"-ia")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>infile_bgn ) ) { cerr<<" ---> Error infile_bgn !"<<endl; exit(1); }
    }
    
    if( strcmp(argv[i],"-ib")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>infile_end ) ) { cerr<<" ---> Error infile_end !"<<endl; exit(1); }
    }
    
    if( strcmp(argv[i],"-o")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>outfile_idx ) ) { cerr<<" ---> Error outfile_idx !"<<endl; exit(1); }
    }   
  }

  cout<<endl<<TString::Format(" ---> infile bgn/end %4d %4d, outfile %4d", infile_bgn, infile_end, outfile_idx)<<endl<<endl;
  
  ////////////////////////////////////////////////////////////////
  
  TString roostr = "";
  
  TFlux *flux_test = new TFlux();
  flux_test->Exe_fluxfile(infile_bgn, infile_end);
  flux_test->Exe_read_tree(outfile_idx);

  return 0;
}
