//////////////////////////// Software description:                                                                                                                                                                 




/////////////////////////////////////////////////////////////

#include "MyHead.h"

int main(int argc,char *argv[]) {

  Float_t g_lat=0,g_lon=0,geo_lat=0,geo_lon=0,sat_ra[3],sat_dec[3];
  UInt_t sec=0;
  UShort_t n_ev=0;
  bool good=true;
  string log_path = output_path_creator(0),root_out_path = output_path_creator(1);
  Int_t chain_entries;                  // -> Variable that stores the whole number of tchain entries
  Double_t perc=0;                      // -> Just used to store the percentage values
  Int_t n_poisson_events=0;             // -> Number of events, generated in the function "fill_maps_for_one_second", following a Poisson distribution. Details in the comments of that function
  
  //////////// Costheta flat binning variables

  Int_t n_bin_lon=360;                  // -> Number of bins along longitude axis
  Double_t lon_bin_min=-180.0;          // -> Set max and min for longitude binning
  Double_t lon_bin_max=180.0;

  Int_t n_bin_lat=180;                  // -> Number of bins along latitude axis
  Double_t lat_bin_min=-90.0;           // -> Set max and min for latitude binning
  Double_t lat_bin_max=90.0;

  Double_t* binning;                    // -> Array used to store the custom binning intervals !!!

  /////////////////////////////////////////////////////////////

  
  //////////////////////// Variable description
  /*

    sat_ra[3]     -> Is the array that stores the right ascension values of the satellite
    sat_dec[3]    -> Is the array that stores the celestial declination values of the satellite
    
    geo_lat       -> Is the geographic latitude
    geo_lon       -> Is the geographic longitude

    g_lat         -> Is the galactic latitude
    g_lon         -> Is the galactic longitude
                     These two are the "absolute coordinates" used to plot the final maps !!!!

    n_ev          -> Is the nuber of triggered events in a second
    sec           -> Is DAMPE's acquisition second number

    good          -> Is the status of the SBI

   */


  //////////////////////////////////////    

  ///////// Set initial arrays' variables to 0

  for(Int_t idx=0; idx<3; idx++) {
    sat_ra[idx]=0;
    sat_dec[idx]=0;
  }
  
  ofstream output_log_file(log_path);     //log file creation !                                                                                                                                                 
  if(!output_log_file.is_open()) {
    cout<<"\n\nCannot create output file! Program finished !"<<endl;
    exit(-1);
  }

  log_file_init(output_log_file);
   
  ////////////////////////////////////   

  TChain *tree= new TChain("SBItree");      //Defining a TTree to read SBI data file
  read_SBI_data(g_lat,g_lon,geo_lat,geo_lon,sat_ra,sat_dec,sec,n_ev,good,tree,sbi_path,output_log_file);    //That function fills the tree reading from the SBI data files

  gRandom->SetSeed(random_seed);  //Set the random seed

  create_binning(n_bin_lat,lat_bin_min,lat_bin_max,binning,true);  //Create binning

  TFile *results_file = new TFile(root_out_path.c_str(),"RECREATE");

  /////////////////////////// Create histos....

  // ---------- Anisotropic maps with infinite and realistic statistics
  
  TH2D* ani_map_inf_stat = new TH2D("ani_map_inf_stat","Anisotropic Map (infinite statistic); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
  TH2D* ani_map_real_stat_1 = new TH2D("ani_map_real_stat_1","Anisotropic Map (real statistic 1); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
  TH2D* ani_map_real_stat_2 = new TH2D("ani_map_real_stat_2","Anisotropic Map (real statistic 2); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);

  TH2D* ani_map_shuf_1 = new TH2D("ani_map_shuf_1","Shuffled Anisotropic Map (real statistic 1); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
  TH2D* ani_map_shuf_2 = new TH2D("ani_map_shuf_2","Shuffled Anisotropic Map (real statistic 2; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries", n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);


  // ---------- Isotropic maps with infinite and realistic statistics 
  
  TH2D* iso_map_inf_stat = new TH2D("iso_map_inf_stat","Isotropic Map (infinite statistic); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
  TH2D* iso_map_real_stat_1 = new TH2D("iso_map_real_stat_1","Isotropic Map (real statistic 1); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
  TH2D* iso_map_real_stat_2 = new TH2D("iso_map_real_stat_2","Isotropic Map (real statistic 2); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
  
  TH2D* iso_map_shuf_1 = new TH2D("iso_map_shuf_1","Shuffled Isotropic Map (real statistic 1); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
  TH2D* iso_map_shuf_2 = new TH2D("iso_map_shuf_2","Shuffled Isotropic Map (real statistic 2); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);


  // ---------- Maps of Satellite Pointing, Exposure, Events and Rate
  
  TH2D* pointing = new TH2D("pointing","Satellite Pointing; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",n_bin_lon,lon_bin_min,lon_bin_max,n_bin_lat,binning);
  TH2D* orbit = new TH2D("orbit","Exposure; Geographic Longitude (#circ);  Geographic Latitude (#circ); Exposure (s)", 180, -180, 180, 90, -90.0, 90.0);
  TH2D* nevents = new TH2D("nevents","Events; Geographic Longitude (#circ);  Geographic Latitude (#circ); Entries", 180, -180, 180, 90, -90.0, 90.0);
  TH2D* rate = new TH2D("rate","Rate; Geographic Longitude (#circ);  Geographic Latitude (#circ); Rate (hz)", 180, -180, 180, 90, -90.0, 90.0);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  chain_entries = tree->GetEntries();
 
  for(Int_t tree_idx=0; tree_idx<chain_entries; tree_idx++) {
    tree->GetEntry(tree_idx);
    if((sec%100000)==0) continue;               //there was a bug in the SBI production
    if(!good) continue;                         //good second

    if (((Double_t)tree_idx/chain_entries)>(perc*0.01)) {
      cout<<"[ "<<perc<<" % ]"<<endl;
      output_log_file<<"[ "<<perc<<" % ]"<<endl;
      perc++;
    }

    // I change the longitudes interval, from [0,360] to [-180,180]
    
    if (g_lon>180) g_lon-=360;
    if (geo_lon>180) geo_lon-=360;

    pointing->Fill(g_lon,g_lat);
    orbit->Fill(geo_lon,geo_lat);           //once per second
    nevents->Fill(geo_lon,geo_lat,n_ev);    //once per event
    rate->Fill(geo_lon,geo_lat,n_ev);       //once per event and eventually divided by the exposure          

    //change to radians

    for (int i=0; i<3; i++) {
      sat_ra[i]*=TMath::DegToRad();
      sat_dec[i]*=TMath::DegToRad();
    }

    
    // Filling isotropic maps -------
    
    fill_maps_for_one_second(output_log_file,n_ev,n_poisson_events,iso_map_inf_stat,iso_map_real_stat_1,iso_map_real_stat_2,sat_ra,sat_dec,true,false);
    fill_maps_for_one_second(output_log_file,n_poisson_events,n_poisson_events,iso_map_inf_stat,iso_map_shuf_1,iso_map_shuf_2,sat_ra,sat_dec,false,false);

    // Filling anisotropic maps -------
    
    fill_maps_for_one_second(output_log_file,n_ev,n_poisson_events,ani_map_inf_stat,ani_map_real_stat_1,ani_map_real_stat_2,sat_ra,sat_dec,true,true);
    fill_maps_for_one_second(output_log_file,n_poisson_events,n_poisson_events,ani_map_inf_stat,ani_map_shuf_1,ani_map_shuf_2,sat_ra,sat_dec,false,true);

  }

  rate->Divide(orbit);
  
  cout<<"\n\nSimulation Completed !\n\n";
  output_log_file << "\n\nSimulation Completed !\n\n";
  
  //////////////////////////////////// Write final histos and closing output files
  
  output_log_file.close();
  results_file->Write();
  results_file->Close();
  
  return 0;
  
}
