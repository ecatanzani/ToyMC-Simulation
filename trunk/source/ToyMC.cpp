//////////////////////////// Software description:




/////////////////////////////////////////////////////////////
#include "MyHead.h"

int main(int argc,char *argv[]) 
{
  Double_t X[3],dir[3],theta,phi,acceptance[trials_number],acceptance_rel_diff[trials_number],analytical_acceptance,costheta,bin_content;
  Int_t perc,accepted_events[trials_number];
  time_t time_stamp;
  string log_path = output_path_creator(time_stamp,0),root_out_path = output_path_creator(time_stamp,1);
  Bool_t accepted_event;
  
  //////////////////////////////////////
  
  ofstream output_log_file(log_path);     //log file creation !
  if(!output_log_file.is_open()) {
    cout<<"\n\nCannot create output file! Program finished !"<<endl;
    exit(-1);
  }

  log_file_init(output_log_file,time_stamp);
  
  ///////////////////////// Defining histos....

  TH1F *hacc = new TH1F("hacc","Acceptance",1000,-0.5,0.5);
  hacc->GetXaxis()->SetTitle("Acceptance (m^{2}sr)");
  hacc->GetYaxis()->SetTitle("Trials");
  hacc->SetLineColor(kBlack);
  hacc->SetFillColor(kOrange-9);

  TH1F *haccreldiff = new TH1F("haccreldiff","Acceptance Relative Error",100, -0.5, 0.5);
  haccreldiff->GetXaxis()->SetTitle("(Acc_{Meas} - Acc_{True}) / Acc_{True}");
  haccreldiff->GetYaxis()->SetTitle("Trials");
  haccreldiff->SetLineColor(kBlack);
  haccreldiff->SetFillColor(kSpring+1);

  TH2D* theta_phi_all = new TH2D("Generated_Theta_Phi","Generated Theta/Phi", 1000, 0, 1, 1000, 0, 2.0*TMath::Pi());
  theta_phi_all->SetXTitle("cos(#theta)");
  theta_phi_all->SetYTitle("#phi (rad)");
  
  TH2D* theta_phi_acceptance = new TH2D("Acceptance","Theta/Phi Acceptance", 1000, 0, 1, 1000, 0, 2.0*TMath::Pi());
  theta_phi_acceptance->SetXTitle("cos(#theta)");
  theta_phi_acceptance->SetYTitle("#phi (rad)");
  theta_phi_acceptance->SetZTitle("Acceptance (m^{2})");


  //////////////////////////////////////////////////////////////////
  
  
  TRandom3* rand_gen = new TRandom3(trandom3_seed);

  TFile *acceptance_results = new TFile(root_out_path.c_str(),"RECREATE");
  
  for(Int_t idx=0; idx<trials_number; idx++)
    accepted_events[idx]=0;
  
  for(Int_t idx_trial=0; idx_trial<trials_number; idx_trial++) {
    printf("\n Try[%i]\n",idx_trial);
    //cout << "\n Try["<<idx_trial<<"]\n";
    output_log_file << "\n Try["<<idx_trial<<"]\n";
    perc=0;
    for(Int_t idx_event=0; idx_event<events_per_trial; idx_event++) {
      if(((Double_t)idx_event/events_per_trial)*100>=perc) {
	perc++;
	printf("\t-> %i%%",perc);
	//cout<<"\t-> "<<perc<<"%";
	output_log_file<<"\t-> "<<perc<<"%";
      }
      generate_coordinate(X,rand_gen);
      generate_theta_phi(theta,phi,rand_gen);
      theta_phi_all->Fill(fabs(cos(theta)),phi);
      obtain_direction(theta,phi,dir);
      check_acceptance(X,dir,theta,accepted_event);
      if(accepted_event) {
	accepted_events[idx_trial]++;
	theta_phi_acceptance->Fill(fabs(cos(theta)),phi);
      }
    } //end for on events for each trial
  } //end for on trials

  ///////////////////// Computing acceptance /////////////////////

  analytical_acceptance = Get_Telescope_Analysical_Acceptance(Hdet,Ldet,Ldet,Ldet,Ldet);

  for(Int_t idx_trial=0; idx_trial<trials_number; idx_trial++) {
    acceptance[idx_trial] = Get_Acceptance( (Double_t)accepted_events[idx_trial]/events_per_trial );
    acceptance_rel_diff[idx_trial] = (acceptance[idx_trial] - analytical_acceptance) / analytical_acceptance;
    cout<<"\n\nAccepted "<<accepted_events[idx_trial]<<" events over "<<events_per_trial<<endl;
    cout<<"Acceptance: "<<acceptance[idx_trial]<<" m^2sr"<<endl;
    cout<<"Acceptance relatice difference: "<<acceptance_rel_diff[idx_trial]<<endl;

    output_log_file << "\n\nAccepted "<<accepted_events[idx_trial]<<" events over "<<events_per_trial<<endl;
    output_log_file << "Acceptance: "<<acceptance[idx_trial]<<" m^2sr"<<endl;
    output_log_file << "Acceptance relatice difference: "<<acceptance_rel_diff[idx_trial]<<endl;
    
    hacc->Fill(acceptance[idx_trial]);
    haccreldiff->Fill(acceptance[idx_trial]);
  }


  ///////////////////// Computing final histos ///////////////////// 
  
  theta_phi_acceptance->Divide(theta_phi_all);
  theta_phi_acceptance->Scale(TMath::Power(Lgen,2));

  for (Int_t xx=1; xx<=theta_phi_acceptance->GetNbinsX(); xx++) {
    costheta = theta_phi_acceptance->GetXaxis()->GetBinCenter(xx);
    for (Int_t yy=1; yy<=theta_phi_acceptance->GetNbinsY(); yy++) {
      bin_content = theta_phi_acceptance->GetBinContent(xx, yy);
      theta_phi_acceptance->SetBinContent(xx, yy, costheta*bin_content);
    }
  }

  
  /////////////////////////// Closing files and writing objects ////////////////////////

  hacc->Fit("gaus");
  hacc->Write();

  haccreldiff->Fit("gaus");
  haccreldiff->Write();

  theta_phi_all->Write();
  theta_phi_acceptance->Write();

  
  output_log_file.close();
  acceptance_results->Close();

  cout<<"\n\nSimulation Completed !\n\n";
  output_log_file << "\n\nSimulation Completed !\n\n";  
  
  return 0;
  
}
