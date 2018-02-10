
///////////////////////////////////// Macro description:








////////////////////////////////////////////////////////////

#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "TFile.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TH2.h"
#include "TColor.h"
#include "TMath.h"
#include "TString.h"
#include "TLatex.h"
#include "TText.h"
#include "TVirtualPad.h"
#include "TPaletteAxis.h"
#include "TPad.h"
#include "TCanvas.h"

TVirtualPad* draw_single_map(TH2* h_to_draw,Bool_t logz);
void draw_grid(TVirtualPad* pad,Double_t canv_x_size,Double_t canv_y_size,Bool_t draw_label=false,Int_t n_coo=6,Int_t lat_label=12);
void draw_pull(TH2* h_ratio,TH2* h_num,TH2* h_den,TH1D* h_pull[],Int_t idx_pool);

using namespace std;

void draw_maps(const TString results_alias) {

  Bool_t logz = false;
  
  ////////////////////////////// Set input and output path
  
  TString input_results_path = "../results/";
  TString pull_results_path = "../results/pull/";
  TString aitoff_maps_results_path = "../results/aitoff_maps/";
  
  input_results_path+=results_alias;
  input_results_path+="_maps_result.root";

  pull_results_path+=results_alias;
  pull_results_path+="_pull_result.root";

  aitoff_maps_results_path+=results_alias;
  aitoff_maps_results_path+="_aitoff_maps_result.root";
  
  ///////////////////////////////////////////////////////

  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(8);

  TFile *input_results_file = new TFile(input_results_path.Data());
  if(input_results_file->IsZombie()) {
    cout<<"\n\nError opening file. Macro finished\n\n";
    exit(-1);
  }
  
  ///////////////////////////////////////////// Satellite's maps
  
  ////////// Draw satellite pointing
  
  TH2* pointing = (TH2D*)(input_results_file->Get("pointing"));
  pointing->SetMinimum(1);
  draw_single_map(pointing,true);

  ////////// Draw satellite orbit                                                                                                                                                                               

  TH2* orbit = (TH2D*)(input_results_file->Get("orbit"));
  pointing->SetMinimum(1);
  draw_single_map(orbit,true);

  ////////// Draw satellite orbit
  
  TH2* nevents = (TH2D*)(input_results_file->Get("nevents"));
  pointing->SetMinimum(1);
  draw_single_map(nevents,true);
  
  ////////// Draw satellite rate
  
  TH2* rate = (TH2D*)(input_results_file->Get("rate"));
  pointing->SetMinimum(1);
  draw_single_map(rate,true);



  /////////////////////////////////////////////////////////////// Isotropic and anisotropic maps


  
  ///////// Draw isotropic maps
  
  TH2* iso_map_inf_stat = (TH2D*)(input_results_file->Get("iso_map_inf_stat"));
  draw_single_map(iso_map_inf_stat,logz);

  TH2* iso_map_real_stat_1 = (TH2D*)(input_results_file->Get("iso_map_real_stat_1"));
  draw_single_map(iso_map_real_stat_1,logz);

  TH2* iso_map_real_stat_2 = (TH2D*)(input_results_file->Get("iso_map_real_stat_2"));
  draw_single_map(iso_map_real_stat_2,logz);

  TH2* iso_map_shuf_1 = (TH2D*)(input_results_file->Get("iso_map_shuf_1"));
  draw_single_map(iso_map_shuf_1,logz);

  TH2* iso_map_shuf_2 = (TH2D*)(input_results_file->Get("iso_map_shuf_2"));
  draw_single_map(iso_map_shuf_2,logz);


  ///////// Draw anisotropic maps  

  TH2* ani_map_inf_stat = (TH2D*)(input_results_file->Get("ani_map_inf_stat"));
  draw_single_map(ani_map_inf_stat,logz);

  TH2* ani_map_real_stat_1 = (TH2D*)(input_results_file->Get("ani_map_real_stat_1"));
  draw_single_map(ani_map_real_stat_1,logz);

  TH2* ani_map_real_stat_2 = (TH2D*)(input_results_file->Get("ani_map_real_stat_2"));
  draw_single_map(ani_map_real_stat_2,logz);

  TH2* ani_map_shuf_1 = (TH2D*)(input_results_file->Get("ani_map_shuf_1"));
  draw_single_map(ani_map_shuf_1,logz);

  TH2* ani_map_shuf_2 = (TH2D*)(input_results_file->Get("ani_map_shuf_2"));
  draw_single_map(ani_map_shuf_2,logz);


  ///////// Draw cloned maps

  TH2* C_ani_map_inf_stat = (TH2D*)(input_results_file->Get("C_ani_map_inf_stat"));  
  TH2* C_ani_map_real_stat_1 = (TH2D*)(input_results_file->Get("C_ani_map_real_stat_1"));
  TH2* C_ani_map_real_stat_2 = (TH2D*)(input_results_file->Get("C_ani_map_real_stat_2"));
  TH2* C_ani_map_shuf_1 = (TH2D*)(input_results_file->Get("C_ani_map_shuf_1"));
  TH2* C_ani_map_shuf_2 = (TH2D*)(input_results_file->Get("C_ani_map_shuf_2"));
  
  TH2* C_iso_map_inf_stat = (TH2D*)(input_results_file->Get("C_iso_map_inf_stat"));  
  TH2* C_iso_map_real_stat_1 = (TH2D*)(input_results_file->Get("C_iso_map_real_stat_1"));
  TH2* C_iso_map_real_stat_2 = (TH2D*)(input_results_file->Get("C_iso_map_real_stat_2"));
  TH2* C_iso_map_shuf_1 = (TH2D*)(input_results_file->Get("C_iso_map_shuf_1"));
  TH2* C_iso_map_shuf_2 = (TH2D*)(input_results_file->Get("C_iso_map_shuf_2"));

  
  /////////////////////////////////////////////////////////////////////////////// MAPS COMPARISON AND PULL BUILDING

  Int_t n_pool=10;
  TH1D *h_pull[n_pool];

  /////////// Isotropic Vs Isotropic with real statistic 1/2

  TH2* iso_map_vs_iso_map_real_stat = (TH2*)(iso_map_real_stat_1->Clone("iso_map_vs_iso_map_real_stat"));
  iso_map_vs_iso_map_real_stat->Divide(iso_map_real_stat_2);
  iso_map_vs_iso_map_real_stat->SetMinimum(0);
  iso_map_vs_iso_map_real_stat->SetMaximum(2);
  iso_map_vs_iso_map_real_stat->SetTitle("Ratio Isotropic/Isotropic (real statistic)");
  iso_map_vs_iso_map_real_stat->GetZaxis()->SetTitle("Ratio");
  draw_single_map(iso_map_vs_iso_map_real_stat,logz);
  draw_pull(iso_map_vs_iso_map_real_stat,iso_map_real_stat_1,iso_map_real_stat_2,h_pull,0);

  
  /////////// Isotropic Vs Isotropic shuffled with real statistic 1
  
  TH2* iso_map_vs_iso_map_shuf_1 = (TH2*)(iso_map_real_stat_1->Clone("iso_map_vs_iso_map_shuf_1"));
  iso_map_vs_iso_map_shuf_1->Divide(iso_map_shuf_1);
  iso_map_vs_iso_map_shuf_1->SetMinimum(0);
  iso_map_vs_iso_map_shuf_1->SetMaximum(2);
  iso_map_vs_iso_map_shuf_1->SetTitle("Ratio Isotropic/Isotropic Shuffled (real statistic 1)");
  iso_map_vs_iso_map_shuf_1->GetZaxis()->SetTitle("Ratio");
  draw_single_map(iso_map_vs_iso_map_shuf_1,logz);
  draw_pull(iso_map_vs_iso_map_shuf_1,iso_map_real_stat_1,iso_map_shuf_1,h_pull,1);

  
  /////////// Anisotropic Vs Isotropic with infinite statistic
  
  TH2* ani_map_vs_iso_map_inf_stat = (TH2*)(ani_map_inf_stat->Clone("ani_map_vs_iso_map_inf_stat"));
  ani_map_vs_iso_map_inf_stat->Divide(iso_map_inf_stat);
  ani_map_vs_iso_map_inf_stat->SetMinimum(0);
  ani_map_vs_iso_map_inf_stat->SetMaximum(2);
  ani_map_vs_iso_map_inf_stat->SetTitle("Ratio Anisotropic/Isotropic (infinite statistic)");
  ani_map_vs_iso_map_inf_stat->GetZaxis()->SetTitle("Ratio");
  draw_single_map(ani_map_vs_iso_map_inf_stat,logz);
  draw_pull(ani_map_vs_iso_map_inf_stat,ani_map_inf_stat,iso_map_inf_stat,h_pull,2);


  /////////// Anisotropic Vs Isotropic with real statistic 1
  
  TH2* ani_map_vs_iso_map_real_stat_1 = (TH2*)(ani_map_real_stat_1->Clone("ani_map_vs_iso_map_real_stat_1"));
  ani_map_vs_iso_map_real_stat_1->Divide(iso_map_real_stat_1);
  ani_map_vs_iso_map_real_stat_1->SetMinimum(0);
  ani_map_vs_iso_map_real_stat_1->SetMaximum(2);
  ani_map_vs_iso_map_real_stat_1->SetTitle("Ratio Anisotropic/Isotropic (real statistic 1)");
  ani_map_vs_iso_map_real_stat_1->GetZaxis()->SetTitle("Ratio");
  draw_single_map(ani_map_vs_iso_map_real_stat_1,logz);
  draw_pull(ani_map_vs_iso_map_real_stat_1,ani_map_real_stat_1,iso_map_real_stat_1,h_pull,3);
  

  /////////// Anisotropic Vs Anisotropic Shuffled with real statistic 1
  
  TH2* ani_map_vs_ani_map_shuf_1 = (TH2*)(ani_map_real_stat_1->Clone("ani_map_vs_ani_map_shuf_1"));
  ani_map_vs_ani_map_shuf_1->Divide(ani_map_shuf_1);
  ani_map_vs_ani_map_shuf_1->SetMinimum(0);
  ani_map_vs_ani_map_shuf_1->SetMaximum(2);
  ani_map_vs_ani_map_shuf_1->SetTitle("Ratio Anisotropic/Anisotropic Shuffled (real statistic 1)");
  ani_map_vs_ani_map_shuf_1->GetZaxis()->SetTitle("Ratio");
  draw_single_map(ani_map_vs_ani_map_shuf_1,logz);
  draw_pull(ani_map_vs_ani_map_shuf_1,ani_map_real_stat_1,ani_map_shuf_1,h_pull,4);                   


  //********** Cloned maps **************

  /////////// Cloned Isotropic Vs Isotropic with real statistic 1/2

  TH2* C_iso_map_vs_iso_map_real_stat = (TH2*)(C_iso_map_real_stat_1->Clone("C_iso_map_vs_iso_map_real_stat"));
  C_iso_map_vs_iso_map_real_stat->Divide(C_iso_map_real_stat_2);
  C_iso_map_vs_iso_map_real_stat->SetMinimum(0);
  C_iso_map_vs_iso_map_real_stat->SetMaximum(2);
  C_iso_map_vs_iso_map_real_stat->SetTitle("Ratio Isotropic/Isotropic - Cloned Maps - (real statistic)");
  C_iso_map_vs_iso_map_real_stat->GetZaxis()->SetTitle("Ratio");
  draw_single_map(C_iso_map_vs_iso_map_real_stat,logz);
  draw_pull(C_iso_map_vs_iso_map_real_stat,C_iso_map_real_stat_1,C_iso_map_real_stat_2,h_pull,5);

  /////////// Cloned Isotropic Vs Isotropic shuffled with real statistic 1
  
  TH2* C_iso_map_vs_iso_map_shuf_1 = (TH2*)(C_iso_map_real_stat_1->Clone("C_iso_map_vs_iso_map_shuf_1"));
  C_iso_map_vs_iso_map_shuf_1->Divide(C_iso_map_shuf_1);
  C_iso_map_vs_iso_map_shuf_1->SetMinimum(0);
  C_iso_map_vs_iso_map_shuf_1->SetMaximum(2);
  C_iso_map_vs_iso_map_shuf_1->SetTitle("Ratio Isotropic/Isotropic Shuffled (real statistic 1)");
  C_iso_map_vs_iso_map_shuf_1->GetZaxis()->SetTitle("Ratio");
  draw_single_map(C_iso_map_vs_iso_map_shuf_1,logz);
  draw_pull(C_iso_map_vs_iso_map_shuf_1,C_iso_map_real_stat_1,C_iso_map_shuf_1,h_pull,6);

  /////////// Anisotropic Vs Isotropic with infinite statistic
  
  TH2* C_ani_map_vs_iso_map_inf_stat = (TH2*)(C_ani_map_inf_stat->Clone("C_ani_map_vs_iso_map_inf_stat"));
  C_ani_map_vs_iso_map_inf_stat->Divide(C_iso_map_inf_stat);
  C_ani_map_vs_iso_map_inf_stat->SetMinimum(0);
  C_ani_map_vs_iso_map_inf_stat->SetMaximum(2);
  C_ani_map_vs_iso_map_inf_stat->SetTitle("Ratio Anisotropic/Isotropic (infinite statistic)");
  C_ani_map_vs_iso_map_inf_stat->GetZaxis()->SetTitle("Ratio");
  draw_single_map(C_ani_map_vs_iso_map_inf_stat,logz);
  draw_pull(C_ani_map_vs_iso_map_inf_stat,C_ani_map_inf_stat,C_iso_map_inf_stat,h_pull,7);

  /////////// Anisotropic Vs Isotropic with real statistic 1
  
  TH2* C_ani_map_vs_iso_map_real_stat_1 = (TH2*)(C_ani_map_real_stat_1->Clone("C_ani_map_vs_iso_map_real_stat_1"));
  C_ani_map_vs_iso_map_real_stat_1->Divide(C_iso_map_real_stat_1);
  C_ani_map_vs_iso_map_real_stat_1->SetMinimum(0);
  C_ani_map_vs_iso_map_real_stat_1->SetMaximum(2);
  C_ani_map_vs_iso_map_real_stat_1->SetTitle("Ratio Anisotropic/Isotropic (real statistic 1)");
  C_ani_map_vs_iso_map_real_stat_1->GetZaxis()->SetTitle("Ratio");
  draw_single_map(C_ani_map_vs_iso_map_real_stat_1,logz);
  draw_pull(C_ani_map_vs_iso_map_real_stat_1,C_ani_map_real_stat_1,C_iso_map_real_stat_1,h_pull,8);
  
  /////////// Anisotropic Vs Anisotropic Shuffled with real statistic 1
  
  TH2* C_ani_map_vs_ani_map_shuf_1 = (TH2*)(C_ani_map_real_stat_1->Clone("C_ani_map_vs_ani_map_shuf_1"));
  C_ani_map_vs_ani_map_shuf_1->Divide(ani_map_shuf_1);
  C_ani_map_vs_ani_map_shuf_1->SetMinimum(0);
  C_ani_map_vs_ani_map_shuf_1->SetMaximum(2);
  C_ani_map_vs_ani_map_shuf_1->SetTitle("Ratio Anisotropic/Anisotropic Shuffled (real statistic 1)");
  C_ani_map_vs_ani_map_shuf_1->GetZaxis()->SetTitle("Ratio");
  draw_single_map(C_ani_map_vs_ani_map_shuf_1,logz);
  draw_pull(C_ani_map_vs_ani_map_shuf_1,C_ani_map_real_stat_1,C_ani_map_shuf_1,h_pull,9);
  
  //////////////////////////// Writing Aitoff maps


  // TFile *aitoff_maps = new TFile(aitoff_maps_results_path.Data(),"RECREATE");
  

  ///////////////////////////////////////////////////////////////////////////////////


  //////////////////////////// Writing Pools
  
  TFile *pools_results = new TFile(pull_results_path.Data(),"RECREATE");
  if(pools_results->IsZombie()) {
    cout<<"\n\nError writing pools results file! Macro finished\n\n";
    exit(-1);
  }
  pools_results->cd();
  gStyle->SetOptStat(1111);
  
  for(Int_t p_idx=0; p_idx<n_pool; p_idx++) {
    if(p_idx==0)
      cout<<"\n\nStandard hstos gaussian fits\n\n";
    if(p_idx==5)
      cout<<"\n\nCloned histos gaussian fits\n\n";

    h_pull[p_idx]->Fit("gaus");
    h_pull[p_idx]->Write();
  }

  pools_results->Write();
  
}


//////////////////////////////////////////////////////////////////////////// Drawing functions.....



TVirtualPad* draw_single_map(TH2* h_to_draw,Bool_t logz) {

  Double_t canv_x_size=1024;
  Double_t canv_y_size=768;
  
  // Add a little style
  
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(8);

  TH2* h_to_draw_clone = (TH2D*)(h_to_draw->Clone(Form("%s_clone", h_to_draw->GetName())));
  
  h_to_draw_clone->GetXaxis()->SetLabelSize(0.0);
  h_to_draw_clone->GetXaxis()->SetTickSize(0.0);
  h_to_draw_clone->GetXaxis()->SetTitleOffset(0.5);
  h_to_draw_clone->GetYaxis()->SetLabelSize(0.0);
  h_to_draw_clone->GetYaxis()->SetTickSize(0.0);
  h_to_draw_clone->GetYaxis()->SetTitleOffset(0.5);

  TCanvas* c_aitoff = new TCanvas(Form("Aitoff_%s",h_to_draw->GetName()),h_to_draw->GetTitle(),1.1*canv_x_size,canv_y_size);
  c_aitoff->cd();

  TPad* pad = new TPad("pad","", 0, 0, 1.0/1.1-0.075, 1.0, 0, 0);
  pad->SetNumber(1);
  pad->Draw();
  c_aitoff->cd(1);
  gPad->SetRightMargin(0.01);

  h_to_draw_clone->Draw("aitoff");                  //Draw the histo with the aitoff projection !! You can also use the mercatore projection, if you want ! 
  
  c_aitoff->cd();
  TPaletteAxis* pal = new TPaletteAxis(1.0/1.1-0.075,0.1,1.0/1.1-0.075+0.05,0.90,h_to_draw_clone);
  pal->Draw();

  draw_grid(c_aitoff->cd(1),canv_x_size,canv_y_size);

  TVirtualPad* pad_2_return = c_aitoff->cd(1);
  pad_2_return->SetLogz(logz);

  c_aitoff->cd();
  c_aitoff->SetLogz(logz);

  return pad_2_return;  
  
}

void draw_grid(TVirtualPad* pad,Double_t canv_x_size,Double_t canv_y_size,Bool_t draw_label,Int_t n_coo,Int_t lat_label) {

  Double_t convx=((1.0/(canv_x_size*1.02))*(10.0/3.0));
  Double_t convy=((1.0/(canv_y_size*0.85*0.8))*(10.0/3.0));
  Bool_t d_grid = true;
  
  Float_t la, lo, _x, _y, z;

  const Int_t Nl = 19;    // Number of drawn latitudes
  const Int_t NL = 19;    // Number of drawn longitudes
  Int_t       M  = 90;

  Double_t _M;
  
  TGraph  *latitudes[Nl];
  TGraph  *longitudes[NL];

  for (int j=0;j<Nl;++j) {
    latitudes[j]=new TGraph();
    latitudes[j]->SetMarkerSize(0.3);
    latitudes[j]->SetMarkerColor(0);      // -> Set grid latitude marker color
    la = -90+180/(Nl-1)*j;
    _M = M;

    // The number of point do draw shold be set according to the latitude. Not all the latitudes need the same number of points

    if (fabs(la)>50.0) {
      _M = M*3.0/5.0;
    }
    if (fabs(la)>60.0) {
      _M = M*2.0/5.0;
    }
    if (fabs(la)>70.0) {
      _M = M*1.0/5.0;
    }
    for (int i=0;i<_M+1;++i) {
      lo = -180+360/_M*i;
      z  = sqrt(1+cos(la*TMath::DegToRad())*cos(lo*TMath::DegToRad()/2));
      _x  = 180*cos(la*TMath::DegToRad())*sin(lo*TMath::DegToRad()/2)/z;
      _y  = 90*sin(la*TMath::DegToRad())/z;
      latitudes[j]->SetPoint(i,_x*convx,_y*convy);
    }
  }


  for (int j=0;j<NL;++j) {
    longitudes[j]=new TGraph();
    longitudes[j]->SetMarkerSize(0.3);
    longitudes[j]->SetMarkerColor(0);    //Set grid longitude marker color
    lo = -180+360/(NL-1)*j;
    for (int i=0;i<M+1;++i) {
      la = -90+180/M*i;
      if (fabs(la)<85.0) {
        z  = sqrt(1+cos(la*TMath::DegToRad())*cos(lo*TMath::DegToRad()/2));
        _x  = 180*cos(la*TMath::DegToRad())*sin(lo*TMath::DegToRad()/2)/z;
        _y  = 90*sin(la*TMath::DegToRad())/z;
        longitudes[j]->SetPoint(i,_x*convx,_y*convy);
      }
    }
  }

  if (d_grid) {
    
    // Draw the grid. That is done drawing each single latitude and longitude point.
    
    pad->cd();
    for (int j=0;j<Nl;++j) {
      if (j!=0 && j!=(Nl-1)) latitudes[j]->Draw("P");
    }
    for (int j=0;j<NL;++j) {
      longitudes[j]->Draw("P");
    }
  }

  // draw label on map

  if (draw_label==true){
    char coo_name[255];
    double* vec_xx = (double*)latitudes[lat_label]->GetX();
    int Npoints = (int)latitudes[lat_label]->GetN();
    double delta_coo = (double)Npoints/(double)n_coo;
    double delta_coo_for_label = 360./(double)n_coo;
    double coo_for_label = -180.;
    int num_coo = 0;
    for (int nn=0; nn<=n_coo; nn++){
      double y = latitudes[lat_label]->Eval(vec_xx[num_coo]);
      double x = vec_xx[num_coo];
      sprintf(coo_name,"%.0f.",coo_for_label);
      TText* tex = new TLatex(x,y,coo_name);
      tex->SetTextFont(40);
      tex->SetTextSize(0.03);
      tex->Draw();
      num_coo = num_coo+delta_coo;
      coo_for_label = coo_for_label+delta_coo_for_label;
    }
  }
}


void draw_pull(TH2* h_ratio,TH2* h_num,TH2* h_den,TH1D* h_pull[],Int_t idx_pool) {

  Double_t pull_min = 0;
  Double_t pull_max = 0;

  Double_t y_rec,y_exp,sigma;
  Double_t pull;

  /*
    
    This first two loops are just needed because I wanto to obtain max and min value for the pull histo.
    These two values will be needed to proprly set the pool's histos max and min value for the X axis.
    The second loop just fills poll's histo ! So obviously I need a couple of loops.

   */

  for(Int_t xx=1; xx<=h_ratio->GetNbinsX(); xx++) {
    for(Int_t yy=1; yy<=h_ratio->GetNbinsY(); yy++) {
      y_rec = h_num->GetBinContent(xx,yy);
      y_exp = h_den->GetBinContent(xx,yy);
      sigma = sqrt(y_exp+y_rec);
      if (sigma>0) {
        pull = (y_rec-y_exp)/sigma;
        if (pull<pull_min)
	  pull_min=pull;
        if (pull>pull_max)
	  pull_max=pull;
      } //end if
    } //end for yy
  } //end for xx

  //h_pull[idx_pool] = new TH1D(Form("%s_pull", h_ratio->GetName()), Form("Pull of %s", h_ratio->GetTitle()), 4*3*5*7, 1.2*pull_min, 1.2*pull_max);
  h_pull[idx_pool] = new TH1D(Form("%s_pull", h_ratio->GetName()), Form("Pull of %s", h_ratio->GetTitle()),150, 1.2*pull_min, 1.2*pull_max);
  h_pull[idx_pool]->GetXaxis()->SetTitle("(y_{rec}-y_{ref})/(#sqrt{#sigma^{2}_{rec}+#sigma^{2}_{ref}})");
  h_pull[idx_pool]->GetYaxis()->SetTitle("Entries");

  for(Int_t xx=1; xx<=h_ratio->GetNbinsX(); xx++) {
    for(Int_t yy=1; yy<=h_ratio->GetNbinsY(); yy++) {
      y_rec = h_num->GetBinContent(xx, yy);
      y_exp = h_den->GetBinContent(xx, yy);
      sigma = sqrt(y_exp+y_rec);
      if (y_rec<0 || y_exp<0)
	printf("%f %f\n", y_rec, y_exp);
      if (sigma>0) {
        pull = (y_rec-y_exp)/sigma;
        h_pull[idx_pool]->Fill(pull);
      }
    }
  }
}
