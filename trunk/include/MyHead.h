///// Head file	///////

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include <ctime>

///// ROOT libraries

#include "TFile.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

using namespace std;

//#define trials_number 1000
//#define events_per_trial 1000

#define trials_number 10
#define events_per_trial 1e+9

const static Int_t n_checked_detectors = 2;                                   //Number of detectors to check with this Toy... 
const static Double_t Lgen = 3.90;                                            //Side of the generation plane in meters 
const static Double_t Ldet = 0.026*22;                                        //Side of the square detector in meters  
const static Double_t Hdet = 0.026*14;                                        //Height of the square detector in meters
const static Double_t h_detectors[n_checked_detectors]={Hdet/2,-Hdet/2};      // !!!!!!!!!!!!!!!  -> -> REMBER TO MODIFY THIS ARRAY WHEN MORE LAYES ARE NEEDEED !!!!

const static Double_t phi_min = 0;
const static Double_t phi_max = 2*TMath::Pi();
const static Double_t cos2_theta_min = 0;
const static Double_t cos2_theta_max = 1;

const static UInt_t trandom3_seed = 9;

////// REMEMBER TO CHAGE output_path VARIABLE IN logs_witer.cpp WHEN CHANGE COMPUTER /////////////

extern string output_path_creator(time_t &time_stamp,const Int_t out_choose);
extern void log_file_init(ofstream &out_file,time_t &time_stamp);
extern void generate_coordinate(Double_t X[],TRandom3 *r_gen);
extern void generate_theta_phi(Double_t &theta,Double_t &phi,TRandom3 *r_gen);
extern void obtain_direction(Double_t theta,Double_t phi,Double_t dir[]);
extern void check_acceptance(Double_t X[],Double_t dir[],Double_t theta,Bool_t &accepted_event);
extern void propagate(Double_t X[],Double_t dir[],Double_t theta,Double_t X_detectors[][2],const Double_t h_detectors[],Int_t n_layer);
extern Double_t Get_Telescope_Analysical_Acceptance(Double_t _l, Double_t _a1, Double_t _b1, Double_t _a2, Double_t _b2);
extern Double_t Get_Acceptance(Double_t ratio);
