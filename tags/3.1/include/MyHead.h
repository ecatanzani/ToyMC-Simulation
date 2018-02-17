///// Head file ///////                                                                                                                                                                                            
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <ctime>
#include <vector>
#include <string>

///// ROOT libraries

#include "TRandom3.h"
#include "TMath.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TTree.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TNamed.h"
#include "TFile.h"
#include "TChain.h"
#include "TColor.h"

///////////////// Other dependancies

#include "orbitStruct.h"

///////////////////////////////////////////////////////

using namespace std;

#define EPS 1.e-12

const static UInt_t random_seed = 9;
const static Int_t number_SBI_files = 3;
static Int_t multiplication_factor = 10;

const static time_t time_stamp=time(0);     //Setting timestamp for the out files 

//////////////////////////// Input path for SBI data files: 

const static TString sbi_path = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/Stuff/toymc-anisotropy-svn/trunk/SBI_data";      
const static TString sbi_subsample = "010";
const static string string_sbi_subsample = "010";      //This is usefull into the function that writes log files and output ROOT files

//////////////////////////// Outh paths for logs and ROOT files:             !!!!!!!!!!!!!!! ATTENCTION !!! Here pat is written for TRUNK directory ! Modify it if necessary !!!!

const static string output_log = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/Stuff/toymc-anisotropy-svn/trunk/logs/";
const static string output_root = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/Stuff/toymc-anisotropy-svn/trunk/results/"; 

//////////////////////////////////////////////////////////////////////////////////                         

//////////////////////////// Input path for the acceptance final histo:

const static string acceptance_final_plot = "/Users/enrico/Documents/Università/Magistrale/Tesi/MyCode/dampe-gacceptance-svn/trunk/results/1515941993_acceptance_result.root";  //  !!!!!!!!!!!!!!! ATTENCTION !!! Remember to wtite the whole path with the correct final acceptance ploth path

////////////////////////////////////////////////////////////////////////////////// 

extern string output_path_creator(const Int_t out_choose);
extern void log_file_init(ofstream &out_file);
extern void read_SBI_data(Float_t &galactic_lat,Float_t &galactic_lon,Float_t &geographic_lat,Float_t &geographic_lon,Float_t sat_ra[],Float_t sat_dec[],UInt_t &sec,UShort_t &n_events,bool &good,TChain* tree,TString sbi_data_path,ofstream &out_file);
extern Bool_t check_sbi_loading(Float_t galactic_lat,Float_t galactic_lon,Float_t geographic_lat,Float_t geographic_lon,Float_t sat_ra[],Float_t sat_dec[],UInt_t sec,UShort_t n_events);
extern void reinitialize_all_variables(Float_t &galactic_lat,Float_t &galactic_lon,Float_t &geographic_lat,Float_t &geographic_lon,Float_t sat_ra[],Float_t sat_dec[],UInt_t &sec,UShort_t &n_events);
extern Bool_t chech_if_null_variable(Float_t in_variable);
extern void create_binning(Int_t n_bin_lat,Double_t lat_bin_min,Double_t lat_bin_max,Double_t* &binning,Bool_t cos_scale);
extern void fill_maps_for_one_second(ofstream &out_file,UShort_t n_ev,Int_t &n_poisson_events,TH2D* map_inf_stat,TH2D* map_real_stat_1,TH2D* map_real_stat_2,Float_t sat_ra[],Float_t sat_dec[],Bool_t poissonize,Bool_t anisotropize);
extern void from_satellite_to_celestial(Float_t ra[],Float_t dec[],double vectorin[],AtPolarVect &vector_out);
extern void AtVect_To_AtPolarVect(double in_vector[],AtPolarVect &vector_out);
extern void invert_AtPolarVect_direction(AtPolarVect vector_out,AtPolarVect &vector_out_inv);
extern void AtPolarVect_to_vector(AtPolarVect &input_polar,double out_array[]);
extern void from_celestial_to_galactic(Double_t ra,Double_t dec,Double_t &l,Double_t &b);
extern void filling_h_cloned(TH2D* h_cloned,TH2D* h_original);
