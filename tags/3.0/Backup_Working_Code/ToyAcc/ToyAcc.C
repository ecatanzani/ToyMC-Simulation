#include "TRandom3.h"
#include "TMath.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TFile.h"

const static Double_t Lgen = 3.90; //Side of the generation plane in meters
const static Double_t Ldet = 0.026*22; //Side of the square detector in meters
const static Double_t Hdet = 0.026*14; //Height of the square detector in meters

void GenXY(Double_t Xgen[3]);
void GenDir(Double_t &theta, Double_t &phi, Double_t dir[3]);
void GenThetaPhi(Double_t &thetagen, Double_t &phigen);
void GetDir(Double_t theta, Double_t phi, Double_t dir[3]);
bool IsInAcceptance( const Double_t Xgen[3], const Double_t dir[3], bool accup, bool accdown );
void Propagate( const Double_t Xgen[3], const Double_t dir[3], const Double_t z, Double_t X[3]);
Double_t GetAcc(Double_t ratio){ return Lgen*Lgen*TMath::Pi()*ratio; } 
Double_t TelescopeAnaliticalAcceptance(Double_t _l, Double_t _a1, Double_t _b1, Double_t _a2, Double_t _b2);

TRandom3 *rgen;

#define NTRY 1

int ToyAcc(int trialspertime=((int)(1e+9))){

  rgen = new TRandom3(9);

  int accepted[NTRY], total[NTRY];
  for(int itry=0; itry<NTRY; itry++) {
    accepted[itry]=0; total[itry]=0;
  }
  
  TH1F *hacc = new TH1F("hacc","Acceptance",100,7.5e+3,7.7e+3);
  hacc->GetXaxis()->SetTitle("Acceptance (m^{2}sr)");
  hacc->GetYaxis()->SetTitle("Trials");
  hacc->SetLineColor(kBlack);
  hacc->SetFillColor(kOrange-9);
  
  TH1F *haccreldiff = new TH1F("haccreldiff","Acceptance Relative Error", 100, -1e-2, +1e-2);
  haccreldiff->GetXaxis()->SetTitle("(Acc_{Meas} - Acc_{True}) / Acc_{True}");
  haccreldiff->GetYaxis()->SetTitle("Trials");
  haccreldiff->SetLineColor(kBlack);
  haccreldiff->SetFillColor(kSpring+1);
  
  TH2D* gen= new TH2D("gen","GeneratedThetaPhi",         1000, 0, 1, 1000, 0, 2.0*TMath::Pi());
  TH2D* accettanza = new TH2D("Acceptance","Acceptance", 1000, 0, 1, 1000, 0, 2.0*TMath::Pi());
  
  for(int itry=0; itry<NTRY; itry++){
    printf("%s::Try[%d]\n",__FUNCTION__,itry);
    int perc=0;
    for(int i=0; i<trialspertime; i++){
      if (((1.0*i)/trialspertime)>(perc*0.01)) {
	perc++;
	printf("%d%%\n", perc);
      }
      total[itry]++;
      Double_t Xgen[3];
      GenXY(Xgen);
      Double_t theta, phi, dir[3];
      GenDir(theta,phi,dir);
      bool accup, accdown, inacc;
      gen->Fill(fabs(TMath::Cos(theta)),phi);
      inacc = IsInAcceptance(Xgen,dir,accup,accdown);
      if (inacc){
	accepted[itry]++;
	accettanza->Fill(fabs(TMath::Cos(theta)),phi);
      }
    }
  }
  
  Double_t acceptance[NTRY];
  Double_t accreldiff[NTRY];
  Double_t analacc = TelescopeAnaliticalAcceptance(Hdet, Ldet, Ldet, Ldet, Ldet);
  for(int itry=0; itry<NTRY; itry++) {
    acceptance[itry] = GetAcc( (Double_t)accepted[itry]/total[itry] );
    accreldiff[itry] = (acceptance[itry]-analacc)/analacc;
    printf("%s::Acceptance[%d]:%d/%d\n",__FUNCTION__,itry,accepted[itry],total[itry]);
    printf("%s::Acceptance[%d]:%fm2sr\n",__FUNCTION__,itry,acceptance[itry]);
    printf("%s::Acceptance[%d]:%.2f\n",__FUNCTION__,itry,accreldiff[itry]);
    
    hacc->Fill(acceptance[itry]);
    haccreldiff->Fill(accreldiff[itry]);
  }
  
  /*
    TCanvas *cacc = new TCanvas("cacc","cacc",1000,500);
     cacc->Divide(2,1);
     cacc->cd(1);
     hacc->Draw("");
     hacc->Fit("gaus");
     cacc->cd(2);
     haccreldiff->Draw("");
     haccreldiff->Fit("gaus");
  */

  //  printf("%f/%f = %f -> %f m2sr\n", accettanza->Integral(), gen->Integral(), accettanza->Integral()/gen->Integral(), Lgen*Lgen*TMath::Pi()*accettanza->Integral()/gen->Integral());
  accettanza->Divide(gen);
  accettanza->Scale(Lgen*Lgen);//we don't need the Pi() term since is A(theta, phi). If it was be A(theta), we would need something like costheta*2*Pi().
  //  For A(theta, phi) the term is L^2*cos(theta). This "cos(theta)" is the scalar product costheta
  for (int xx=1; xx<=accettanza->GetNbinsX(); xx++) {
    double costheta = accettanza->GetXaxis()->GetBinCenter(xx);
    for (int yy=1; yy<=accettanza->GetNbinsY(); yy++) {
      double bc = accettanza->GetBinContent(xx, yy);
      accettanza->SetBinContent(xx, yy, costheta*bc);
    }
  }  
  //  printf("    %f m2sr\n", accettanza->Integral("width"));
  
  TCanvas *c = new TCanvas();
  c->cd();
  accettanza->SetXTitle("cos (#theta)");
  accettanza->SetYTitle("#phi (rad)");
  accettanza->SetZTitle("Acceptance (m^{2})");
  accettanza->DrawCopy("colz");
  
  TFile *fout= new TFile ("acceptance.root","recreate");
  fout->cd();
  accettanza->Write();
  fout->Close();
  
  /*  
      TCanvas *d = new TCanvas();
      d->cd();
      gen->Draw("colz");*/
  
  return 0;
}

  
void GenXY(Double_t Xgen[3]){
  Xgen[0] = rgen->Uniform( -Lgen/2., Lgen/2.);
  Xgen[1] = rgen->Uniform( -Lgen/2., Lgen/2.);
  Xgen[2] = -Lgen/2;
  return;
}

void GenDir(Double_t &theta, Double_t &phi, Double_t dir[3]){
  GenThetaPhi(theta,phi);
  GetDir(theta,phi,dir);
  return;
  }

void GenThetaPhi(Double_t &thetagen, Double_t &phigen){
  static const Double_t ctheta2min = 0;
  static const Double_t ctheta2max = pow(cos(TMath::Pi()),2);
  static const Double_t phimin = 0;
  static const Double_t phimax = 2.0*TMath::Pi();
  Double_t ctheta2 = rgen->Uniform( ctheta2min, ctheta2max );
  phigen = rgen->Uniform( phimin, phimax );
  thetagen = TMath::Pi() - acos( sqrt(ctheta2) );
  return;
}


void GetDir(Double_t theta, Double_t phi, Double_t dir[3]){
  dir[0] = sin(theta)*cos(phi);
  dir[1] = sin(theta)*sin(phi);
  dir[2] = cos(theta);
  return;
}

void Propagate( const Double_t Xgen[3], const Double_t dir[3], const Double_t z, Double_t X[3]){
  Double_t dz = fabs(Xgen[2]-z); 
  Double_t r =  dz / dir[2];
  X[0] = Xgen[0] + r*dir[0];
  X[1] = Xgen[1] + r*dir[1];
  X[2] = z;
  return;
}

bool IsInAcceptance( const Double_t Xgen[3], const Double_t dir[3], bool accup, bool accdown ){

  accup =false;
  accdown=false;
  
  const static Double_t zup   = +Hdet/2;
  const static Double_t zdown = -Hdet/2;
  Double_t Xup[3],Xdown[3];
  Propagate( Xgen, dir, zup,   Xup   );
  Propagate( Xgen, dir, zdown, Xdown );
  if( fabs(Xup[0]  -0)<Ldet/2 && fabs(Xup[1]  -0)<Ldet/2 ) { accup=true; } 
  if( fabs(Xdown[0]-0)<Ldet/2 && fabs(Xdown[1]-0)<Ldet/2 ) { accdown=true; }
  
  if(accup&&accdown) return true;
  else return false;
}

Double_t TelescopeAnaliticalAcceptance(Double_t _l, Double_t _a1, Double_t _b1, Double_t _a2, Double_t _b2){
  //l: distance from planes
  //a1,b1: sides of upper squares
  //a2,b2: sides of lower squares
  Double_t a = 0.5*(_a1+_a2);
  Double_t b = 0.5*(_b1+_b2);
  Double_t c = 0.5*(_a1-_a2);
  Double_t d = 0.5*(_b1-_b2);
  Double_t a2 = pow(a,2);
  Double_t b2 = pow(b,2);
  Double_t c2 = pow(c,2);
  Double_t d2 = pow(d,2);
  Double_t l2 = pow(_l,2);

  Double_t G=0;
  G += l2 * log( (l2+a2+d2)/(l2+a2+b2) * (l2+c2+b2)/(l2+c2+d2) );
  G += 2*a*sqrt(l2+b2)*atan( a/sqrt(l2+a2) );
  G += 2*b*sqrt(l2+a2)*atan( b/sqrt(l2+b2) );
  G += 2*c*sqrt(l2+d2)*atan( c/sqrt(l2+d2) );
  G += 2*d*sqrt(l2+c2)*atan( d/sqrt(l2+c2) );
  G -= 2*a*sqrt(l2+d2)*atan( a/sqrt(l2+d2) );
  G -= 2*b*sqrt(l2+c2)*atan( b/sqrt(l2+c2) );
  G -= 2*c*sqrt(l2+b2)*atan( c/sqrt(l2+b2) );
  G -= 2*d*sqrt(l2+a2)*atan( d/sqrt(l2+a2) );
  
  return G;
}
