
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
#include "TRandom3.h"

TVirtualPad* draw_single_map(TH2* h_to_draw,Bool_t logz);
void draw_grid(TVirtualPad* pad,Double_t canv_x_size,Double_t canv_y_size,Bool_t draw_label=false,Int_t n_coo=6,Int_t lat_label=12);
void draw_pull(TH2* h_ratio,TH2* h_num,TH2* h_den,TH1D* h_pull[],Int_t idx_pool);
void filling_entries(TH1I* h_entries,TH2* h_map);
void pull_correlation_MC(TString out_path);
void obtain_correlation_map(TH2* h_map1,TH2* h_map2,TH2I* h_correl);

using namespace std;

void draw_maps(const TString results_alias) {
    Bool_t logz = false;
  
    ////////////////////////////// Set input and output path
  
    TString input_results_path = "../results/";
    TString pull_results_path = "../results/pull/";
    TString aitoff_maps_results_path = "../results/aitoff_maps/";
    TString pull_correlation_path = "../results/pull/correlation_study/";
    
    input_results_path+=results_alias;
    input_results_path+="_maps_result.root";

    pull_results_path+=results_alias;
    pull_results_path+="_pull_result.root";

    aitoff_maps_results_path+=results_alias;
    aitoff_maps_results_path+="_aitoff_maps_result.root";

    pull_correlation_path+=results_alias;
    pull_correlation_path+="_pool_correlation_study.root";
    
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


    ///////// Load cloned maps - I'll use that ones just for the pools !!
    
    TH2* C_iso_map_real_stat_1 = (TH2D*)(input_results_file->Get("C_iso_map_real_stat_1"));
    TH2* C_iso_map_real_stat_2 = (TH2D*)(input_results_file->Get("C_iso_map_real_stat_2"));
    TH2* C_ani_map_real_stat_1 = (TH2D*)(input_results_file->Get("C_ani_map_real_stat_1"));
    TH2* C_ani_map_real_stat_2 = (TH2D*)(input_results_file->Get("C_ani_map_real_stat_2"));
    TH2* C_iso_map_shuf_1 = (TH2D*)(input_results_file->Get("C_iso_map_shuf_1"));
    TH2* C_iso_map_shuf_2 = (TH2D*)(input_results_file->Get("C_iso_map_shuf_2"));
    TH2* C_ani_map_shuf_1 = (TH2D*)(input_results_file->Get("C_ani_map_shuf_1"));
    TH2* C_ani_map_shuf_2 = (TH2D*)(input_results_file->Get("C_ani_map_shuf_2"));
    TH2* C_ani_map_inf_stat = (TH2D*)(input_results_file->Get("C_ani_map_inf_stat"));
    TH2* C_iso_map_inf_stat = (TH2D*)(input_results_file->Get("C_iso_map_inf_stat"));
      
    /////////////////////////////////////////////////////////////////////////////// MAPS COMPARISON AND PULL BUILDING

    Int_t n_pool=10;
    TH1D *h_pull[n_pool];
    
    gStyle->SetOptStat(1);
  
    TH1I *iso_map_real_stat_1_entries = new TH1I("iso_map_real_stat_1_entries","Entries Isotropic Map Real Statistic 1",100,0,800);
    TH1I *iso_map_real_stat_2_entries = new TH1I("iso_map_real_stat_2_entries","Entries Isotropic Map Real Statistic 2",100,0,800);
    TH1I *iso_map_shuf_1_entries = new TH1I("iso_map_shuf_1_entries","Entries Isotropic Map Shffled RS 1",100,0,800);
    TH1I *iso_map_inf_stat_entries = new TH1I("iso_map_inf_stat_entries","Entries Isotropic Map Infinite Statistic",100,0,13000);
    TH1I *ani_map_inf_stat_entries = new TH1I("ani_map_inf_stat_entries","Entries Anisotropic Map Infinite Statistic",100,0,22000);
    TH1I *ani_map_real_stat_1_entries = new TH1I("ani_map_real_stat_1_entries","Entries Anisotropic Map Real Statistic 1",100,0,1300);
    TH1I *ani_map_shuf_1_entries = new TH1I("ani_map_shuf_1_entries","Entries Anisotropic Map Shuffled RS 1",100,0,1300);
  
    /////////// Filling entries histo

    filling_entries(iso_map_real_stat_1_entries,iso_map_real_stat_1);
    filling_entries(iso_map_real_stat_2_entries,iso_map_real_stat_2);
    filling_entries(iso_map_shuf_1_entries,iso_map_shuf_1);
    filling_entries(iso_map_inf_stat_entries,iso_map_inf_stat);
    filling_entries(ani_map_inf_stat_entries,ani_map_inf_stat);
    filling_entries(ani_map_real_stat_1_entries,ani_map_real_stat_1);
    filling_entries(ani_map_shuf_1_entries,ani_map_shuf_1);
    
    /////////// Obtaining maps correlation histos !
    
    TH2I *iso_real_stat_correl = new TH2I("iso_real_stat_correl","Correlation Isotropic Maps Real Statistic",1000,0,800,1000,0,800); // isotropic real stat 1 / real stat 2
    TH2I *iso_real_shuf_correl = new TH2I("iso_real_shuf_correl","Correlation Isotropic Maps Real Statistic and Shuffled",1000,0,800,1000,0,800); // isotropic real stat 1 / shuf 1
    TH2I *ani_iso_inf_stat_correl = new TH2I("ani_iso_inf_stat_correl","Correlation Anisotropic/Isotropic Map Infinite Statistic",1000,0,22000,1000,0,22000); // anisotropic/isotropic infinite statistic
    TH2I *ani_iso_real_stat_1_correl = new TH2I("ani_iso_real_stat_1_correl","Correlation Anisotropic/Isotropic Map Real Statistic 1",1000,0,1300,1000,0,1300); // anisotropic/isotropic real statistic 1
    TH2I *ani_real_shuf_correl = new TH2I("ani_real_shuf_correl","Correlation Anisotropic Maps Real Statistic and Shuffled",1000,0,1300,1000,0,1300); // anisotropic real statistic/anisotropic shuffled real statistic 1
    
    obtain_correlation_map(iso_map_real_stat_1,iso_map_real_stat_2,iso_real_stat_correl);
    obtain_correlation_map(iso_map_real_stat_1,iso_map_shuf_1,iso_real_shuf_correl);
    obtain_correlation_map(ani_map_inf_stat,iso_map_inf_stat,ani_iso_inf_stat_correl);
    obtain_correlation_map(ani_map_real_stat_1,iso_map_real_stat_1,ani_iso_real_stat_1_correl);
    obtain_correlation_map(ani_map_real_stat_1,ani_map_shuf_1,ani_real_shuf_correl);
    
    /////////// Isotropic Vs Isotropic with real statistic 1/2

    TH2* iso_map_vs_iso_map_real_stat = (TH2*)(iso_map_real_stat_1->Clone("iso_map_vs_iso_map_real_stat"));
    iso_map_vs_iso_map_real_stat->Divide(iso_map_real_stat_2);
    iso_map_vs_iso_map_real_stat->SetMinimum(0);
    iso_map_vs_iso_map_real_stat->SetMaximum(2);
    iso_map_vs_iso_map_real_stat->SetTitle("Ratio Isotropic/Isotropic (real statistic)");
    iso_map_vs_iso_map_real_stat->GetZaxis()->SetTitle("Ratio");
    draw_single_map(iso_map_vs_iso_map_real_stat,logz);
    gStyle->SetOptStat(1);
    draw_pull(iso_map_vs_iso_map_real_stat,iso_map_real_stat_1,iso_map_real_stat_2,h_pull,0);
    
  
    /////////// Isotropic Vs Isotropic shuffled with real statistic 1
  
    TH2* iso_map_vs_iso_map_shuf_1 = (TH2*)(iso_map_real_stat_1->Clone("iso_map_vs_iso_map_shuf_1"));
    iso_map_vs_iso_map_shuf_1->Divide(iso_map_shuf_1);
    iso_map_vs_iso_map_shuf_1->SetMinimum(0);
    iso_map_vs_iso_map_shuf_1->SetMaximum(2);
    iso_map_vs_iso_map_shuf_1->SetTitle("Ratio Isotropic/Isotropic Shuffled (real statistic 1)");
    iso_map_vs_iso_map_shuf_1->GetZaxis()->SetTitle("Ratio");
    draw_single_map(iso_map_vs_iso_map_shuf_1,logz);
    gStyle->SetOptStat(1);
    draw_pull(iso_map_vs_iso_map_shuf_1,iso_map_real_stat_1,iso_map_shuf_1,h_pull,1);

      
    /////////// Anisotropic Vs Isotropic with infinite statistic
      
    TH2* ani_map_vs_iso_map_inf_stat = (TH2*)(ani_map_inf_stat->Clone("ani_map_vs_iso_map_inf_stat"));
    ani_map_vs_iso_map_inf_stat->Divide(iso_map_inf_stat);
    ani_map_vs_iso_map_inf_stat->SetMinimum(0);
    ani_map_vs_iso_map_inf_stat->SetMaximum(2);
    ani_map_vs_iso_map_inf_stat->SetTitle("Ratio Anisotropic/Isotropic (infinite statistic)");
    ani_map_vs_iso_map_inf_stat->GetZaxis()->SetTitle("Ratio");
    draw_single_map(ani_map_vs_iso_map_inf_stat,logz);
    gStyle->SetOptStat(1);
    draw_pull(ani_map_vs_iso_map_inf_stat,ani_map_inf_stat,iso_map_inf_stat,h_pull,2);


    /////////// Anisotropic Vs Isotropic with real statistic 1
      
    TH2* ani_map_vs_iso_map_real_stat_1 = (TH2*)(ani_map_real_stat_1->Clone("ani_map_vs_iso_map_real_stat_1"));
    ani_map_vs_iso_map_real_stat_1->Divide(iso_map_real_stat_1);
    ani_map_vs_iso_map_real_stat_1->SetMinimum(0);
    ani_map_vs_iso_map_real_stat_1->SetMaximum(2);
    ani_map_vs_iso_map_real_stat_1->SetTitle("Ratio Anisotropic/Isotropic (real statistic 1)");
    ani_map_vs_iso_map_real_stat_1->GetZaxis()->SetTitle("Ratio");
    draw_single_map(ani_map_vs_iso_map_real_stat_1,logz);
    gStyle->SetOptStat(1);
    draw_pull(ani_map_vs_iso_map_real_stat_1,ani_map_real_stat_1,iso_map_real_stat_1,h_pull,3);
      

    /////////// Anisotropic Vs Anisotropic Shuffled with real statistic 1
      
    TH2* ani_map_vs_ani_map_shuf_1 = (TH2*)(ani_map_real_stat_1->Clone("ani_map_vs_ani_map_shuf_1"));
    ani_map_vs_ani_map_shuf_1->Divide(ani_map_shuf_1);
    ani_map_vs_ani_map_shuf_1->SetMinimum(0);
    ani_map_vs_ani_map_shuf_1->SetMaximum(2);
    ani_map_vs_ani_map_shuf_1->SetTitle("Ratio Anisotropic/Anisotropic Shuffled (real statistic 1)");
    ani_map_vs_ani_map_shuf_1->GetZaxis()->SetTitle("Ratio");
    draw_single_map(ani_map_vs_ani_map_shuf_1,logz);
    gStyle->SetOptStat(1);
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
    gStyle->SetOptStat(1);
    draw_pull(C_iso_map_vs_iso_map_real_stat,C_iso_map_real_stat_1,C_iso_map_real_stat_2,h_pull,5);

    /////////// Cloned Isotropic Vs Isotropic shuffled with real statistic 1
    
    TH2* C_iso_map_vs_iso_map_shuf_1 = (TH2*)(C_iso_map_real_stat_1->Clone("C_iso_map_vs_iso_map_shuf_1"));
    C_iso_map_vs_iso_map_shuf_1->Divide(C_iso_map_shuf_1);
    C_iso_map_vs_iso_map_shuf_1->SetMinimum(0);
    C_iso_map_vs_iso_map_shuf_1->SetMaximum(2);
    C_iso_map_vs_iso_map_shuf_1->SetTitle("Ratio Isotropic/Isotropic Shuffled (real statistic 1)");
    C_iso_map_vs_iso_map_shuf_1->GetZaxis()->SetTitle("Ratio");
    draw_single_map(C_iso_map_vs_iso_map_shuf_1,logz);
    gStyle->SetOptStat(1);
    draw_pull(C_iso_map_vs_iso_map_shuf_1,C_iso_map_real_stat_1,C_iso_map_shuf_1,h_pull,6);

    /////////// Anisotropic Vs Isotropic with infinite statistic
      
    TH2* C_ani_map_vs_iso_map_inf_stat = (TH2*)(C_ani_map_inf_stat->Clone("C_ani_map_vs_iso_map_inf_stat"));
    C_ani_map_vs_iso_map_inf_stat->Divide(C_iso_map_inf_stat);
    C_ani_map_vs_iso_map_inf_stat->SetMinimum(0);
    C_ani_map_vs_iso_map_inf_stat->SetMaximum(2);
    C_ani_map_vs_iso_map_inf_stat->SetTitle("Ratio Anisotropic/Isotropic (infinite statistic)");
    C_ani_map_vs_iso_map_inf_stat->GetZaxis()->SetTitle("Ratio");
    draw_single_map(C_ani_map_vs_iso_map_inf_stat,logz);
    gStyle->SetOptStat(1);
    draw_pull(C_ani_map_vs_iso_map_inf_stat,C_ani_map_inf_stat,C_iso_map_inf_stat,h_pull,7);

    /////////// Anisotropic Vs Isotropic with real statistic 1
      
    TH2* C_ani_map_vs_iso_map_real_stat_1 = (TH2*)(C_ani_map_real_stat_1->Clone("C_ani_map_vs_iso_map_real_stat_1"));
    C_ani_map_vs_iso_map_real_stat_1->Divide(C_iso_map_real_stat_1);
    C_ani_map_vs_iso_map_real_stat_1->SetMinimum(0);
    C_ani_map_vs_iso_map_real_stat_1->SetMaximum(2);
    C_ani_map_vs_iso_map_real_stat_1->SetTitle("Ratio Anisotropic/Isotropic (real statistic 1)");
    C_ani_map_vs_iso_map_real_stat_1->GetZaxis()->SetTitle("Ratio");
    draw_single_map(C_ani_map_vs_iso_map_real_stat_1,logz);
    gStyle->SetOptStat(1);
    draw_pull(C_ani_map_vs_iso_map_real_stat_1,C_ani_map_real_stat_1,C_iso_map_real_stat_1,h_pull,8);
      
    /////////// Anisotropic Vs Anisotropic Shuffled with real statistic 1
      
    TH2* C_ani_map_vs_ani_map_shuf_1 = (TH2*)(C_ani_map_real_stat_1->Clone("C_ani_map_vs_ani_map_shuf_1"));
    C_ani_map_vs_ani_map_shuf_1->Divide(C_ani_map_shuf_1);
    C_ani_map_vs_ani_map_shuf_1->SetMinimum(0);
    C_ani_map_vs_ani_map_shuf_1->SetMaximum(2);
    C_ani_map_vs_ani_map_shuf_1->SetTitle("Ratio Anisotropic/Anisotropic Shuffled (real statistic 1)");
    C_ani_map_vs_ani_map_shuf_1->GetZaxis()->SetTitle("Ratio");
    draw_single_map(C_ani_map_vs_ani_map_shuf_1,logz);
    gStyle->SetOptStat(1);
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
    gStyle->SetOptStat(1);
    gStyle->SetOptStat(1111);
    
    for(Int_t p_idx=0; p_idx<n_pool; p_idx++) {
    if(p_idx==0)
        cout<<"\n\nStandard hstos gaussian fits\n\n";
        if(p_idx==5)
            cout<<"\n\nCloned histos gaussian fits\n\n";

        h_pull[p_idx]->Fit("gaus");
        h_pull[p_idx]->Write();
    }
    
    

    iso_map_real_stat_1_entries->Write();
    iso_map_real_stat_2_entries->Write();
    iso_map_shuf_1_entries->Write();
    iso_map_inf_stat_entries->Write();
    ani_map_inf_stat_entries->Write();
    ani_map_real_stat_1_entries->Write();
    ani_map_shuf_1_entries->Write();
    
    iso_real_stat_correl->Write();
    iso_real_shuf_correl->Write();
    ani_iso_inf_stat_correl->Write();
    ani_iso_real_stat_1_correl->Write();
    ani_real_shuf_correl->Write();
    
    pools_results->Write();
    pools_results->Close();
    
    pull_correlation_MC(pull_correlation_path);
      
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
          if(y_rec<10 || y_exp<10) {
        //cout<<"\nLow single bin statistic !!";
        continue;
          }
          //sigma = sqrt(y_exp+y_rec);
          sigma = sqrt(y_exp)+sqrt(y_rec);
          //sigma = sqrt(y_exp);
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
          if(y_rec<10 || y_exp<10) {
        //cout<<"\nLow single bin statistic !!";
        continue;
          }
          sigma = sqrt(y_exp)+sqrt(y_rec);
          //sigma = sqrt(y_exp);
          if (y_rec<0 || y_exp<0)
        printf("%f %f\n", y_rec, y_exp);
          if (sigma>0) {
            pull = (y_rec-y_exp)/sigma;
            h_pull[idx_pool]->Fill(pull);
          }
        }
      }
    }

void filling_entries(TH1I* h_entries,TH2* h_map) {
    for(Int_t idx_bx=1; idx_bx<=h_map->GetNbinsX(); idx_bx++)
        for(Int_t idx_by=1; idx_by<=h_map->GetNbinsY(); idx_by++)
            if(h_map->GetBinContent(idx_bx,idx_by)>0)
                h_entries->Fill(h_map->GetBinContent(idx_bx,idx_by));
    
}

void pull_correlation_MC(TString out_path) {
    UInt_t seed=22;
    Double_t n_events=1e+9;
    Int_t ent_1,ent_2;
    Double_t sigma,pull;
    TRandom3 *rand = new TRandom3(seed);
    Double_t tmp_bin1,tmp_bin2;
    
    TFile *results = new TFile(out_path.Data(),"RECREATE");
    if(results->IsZombie()) {
        cout<<"\n\nError writing pool correlation study output file\n\n";
        exit(-1);
    }
    
    // Two independent uniform distributions
    
    gStyle->SetOptStat(1);
    TH1D *h_uniform1 = new TH1D("h_uniform1","Unifrom Distribution",1e+5,0,1);
    TH1D *h_uniform2 = new TH1D("h_uniform2","Uniform Distribution",1e+5,0,1);
    TH1D *h_pull_uniform = new TH1D("h_pool_uniform","Pool from Uniform Distributions",1000,-5,5);
    TH1I *h_uniform1_entries = new TH1I("h_uniform1_entries","Entries Distribution Uniform 1",1000,9e+3,11e+3);
    TH1I *h_uniform2_entries = new TH1I("h_uniform2_entries","Entries Distribution Uniform 2",1000,9e+3,11e+3);
    TH2D *uniform_correl = new TH2D("uniform_correl","Correlation Uniform Distribution",1000,9e+3,11e+3,1000,9e+3,11e+3);
    
    //TH2I *uniform_correl = new TH2I("uniform_correl","Uniform Maps Correlation",1000)
    
    for(Int_t idx_l=0; idx_l<n_events; idx_l++) {
        h_uniform1->Fill(rand->Uniform(0,1));
        h_uniform2->Fill(rand->Uniform(0,1));
    }
    
    for(Int_t idx_bx=1; idx_bx<=h_uniform1->GetNbinsX(); idx_bx++) {
        h_uniform1_entries->Fill(h_uniform1->GetBinContent(idx_bx));
        h_uniform2_entries->Fill(h_uniform2->GetBinContent(idx_bx));
    }
    
    for(Int_t idx_bx=1; idx_bx<=h_uniform1->GetNbinsX(); idx_bx++) {
        tmp_bin1=h_uniform1->GetBinContent(idx_bx);
        tmp_bin2=h_uniform2->GetBinContent(idx_bx);
        uniform_correl->Fill(tmp_bin1,tmp_bin2);
    }
        
    for(Int_t idx_b=1; idx_b<=h_uniform1->GetNbinsX(); idx_b++) {
        ent_1=h_uniform1->GetBinContent(idx_b);
        ent_2=h_uniform2->GetBinContent(idx_b);
        if(ent_1<10 || ent_2<10)
            continue;
        //sigma=sqrt(ent_1+ent_2);
        sigma=sqrt(ent_1+ent_2+2*(uniform_correl->GetCorrelationFactor())*sqrt(ent_1)*sqrt(ent_2));
        if(sigma>0) {
            pull=(ent_1-ent_2)/sigma;
            h_pull_uniform->Fill(pull);
        }
    }
    
    // Two independent gaussian distributions
    
    gStyle->SetOptStat(1);
    TH1D *h_gaus1 = new TH1D("h_gaus1","Gaussian Distribution",1e+5,-5,5);
    TH1D *h_gaus2 = new TH1D("h_gaus2","Gaussian Distribution",1e+5,-5,5);
    TH1D *h_pull_gaus = new TH1D("h_pool_gaus","Pool from Gaussian Distributions",1000,-5,5);
    TH1I *h_gaus1_entries = new TH1I("h_gaus1_entries","Entries Distribution Uniform 1",1000,0,45e+3);
    TH1I *h_gaus2_entries = new TH1I("h_gaus2_entries","Entries Distribution Uniform 2",1000,0,45e+3);
    TH2D *gaus_correl = new TH2D("gaus_correl","Correlation Gaussian Distribution",1000,0,45e+3,1000,0,45e+3);
    
    for(Int_t idx_l=0; idx_l<n_events; idx_l++) {
        h_gaus1->Fill(rand->Gaus(0,1));
        h_gaus2->Fill(rand->Gaus(0,1));
    }
    
    for(Int_t idx_bx=1; idx_bx<=h_gaus1->GetNbinsX(); idx_bx++) {
        h_gaus1_entries->Fill(h_gaus1->GetBinContent(idx_bx));
        h_gaus2_entries->Fill(h_gaus2->GetBinContent(idx_bx));
    }
    
    for(Int_t idx_bx=1; idx_bx<=h_gaus1->GetNbinsX(); idx_bx++) {
        tmp_bin1=h_gaus1->GetBinContent(idx_bx);
        tmp_bin2=h_gaus2->GetBinContent(idx_bx);
        gaus_correl->Fill(tmp_bin1,tmp_bin2);
    }
    
    for(Int_t idx_b=1; idx_b<=h_gaus1->GetNbinsX(); idx_b++) {
        ent_1=h_gaus1->GetBinContent(idx_b);
        ent_2=h_gaus2->GetBinContent(idx_b);
        if(ent_1<10 || ent_2<10)
            continue;
        sigma=sqrt(ent_1+ent_2+2*(gaus_correl->GetCorrelationFactor())*sqrt(ent_1)*sqrt(ent_2));
        //sigma=sqrt(ent_1)+sqrt(ent_2);
        //sigma=sqrt(ent_1+ent_2);
        if(sigma>0) {
            pull=(ent_1-ent_2)/sigma;
            h_pull_gaus->Fill(pull);
        }
    }
    
    
    results->Write();
    results->Close();
    
}

void obtain_correlation_map(TH2* h_map1,TH2* h_map2,TH2I* h_correl) {
    Int_t tmp_bin1,tmp_bin2;
    for(Int_t idx_bx=1; idx_bx<=h_map1->GetNbinsX(); idx_bx++)
        for(Int_t idx_by=1; idx_by<=h_map1->GetNbinsY(); idx_by++) {
            tmp_bin1=h_map1->GetBinContent(idx_bx,idx_by);
            tmp_bin2=h_map2->GetBinContent(idx_bx,idx_by);
            h_correl->Fill(tmp_bin1,tmp_bin2);
        }
}
