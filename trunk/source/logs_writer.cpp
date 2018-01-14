
#include "MyHead.h"

string output_path_creator(time_t &time_stamp,const Int_t out_choose) {
  // out_choose == 0 means we are creating the path for a log file
  // out_choose == 0 means we are creating the path for a ROOT output file

  string output_path;
  time_t t_stamp;
  
  switch(out_choose) {

  case 0: 
    output_path = "/Users/enrico/Documents/Università/Magistrale/Tesi/Anisotropy/MyCode/ToyMC_Anisotrypy_Study/toymc-anisotropy-study-code/logs/";
    t_stamp = time(0);
    time_stamp = t_stamp;
    output_path+=to_string((long long)t_stamp);
    output_path+=".txt";
    cout<<"\nWritten log file: -> \t "<<output_path<<endl;
    break;

  case 1:
    output_path = "/Users/enrico/Documents/Università/Magistrale/Tesi/Anisotropy/MyCode/ToyMC_Anisotrypy_Study/toymc-anisotropy-study-code/results/";
    t_stamp = time(0);
    time_stamp = t_stamp;
    output_path+=to_string((long long)t_stamp);
    output_path+="_acceptance_result.root";
    cout<<"\nWritten ROOT file: -> \t "<<output_path<<endl;
    break;

  }
  
  return output_path;
}

void log_file_init(ofstream &out_file,time_t &time_stamp) {
  out_file << "********************* Automatic Log File Generator *******************"<<endl<<endl;
  
  out_file << "////////////////////////// Simulation Parameters //////////////////////////"<<endl<<endl;
  out_file << "Simulation timestamp: "<<time_stamp<<endl;
  out_file << "Simulation TRandom3 seed: "<<trandom3_seed<<endl;
  out_file << "Number of trials: "<<trials_number<<endl;
  out_file << "Number of events per trial: "<<events_per_trial<<endl<<endl;
  
  out_file << "////////////////////////// Detector Parameters //////////////////////////"<<endl<<endl;
  out_file << "Generation plane (in meters): "<<Lgen<<endl;
  out_file << "Detector lenght (in meters): "<<Ldet<<endl;
  out_file << "Detector high (in meters): "<<Hdet<<endl;
  out_file << "Number of checked layers for geometrical acceptace computing: "<<n_checked_detectors<<endl;

  out_file << "*\n*\n*\n"<<endl<<endl;
}
