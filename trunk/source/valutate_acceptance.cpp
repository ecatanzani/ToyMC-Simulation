
#include "MyHead.h"

void check_acceptance(Double_t X[],Double_t dir[],Double_t theta,Bool_t &accepted_event) {
  Double_t X_detectors[3][n_checked_detectors];
  Bool_t checked_layer[n_checked_detectors]={false,false};
  
  for(Int_t idx_layer=0; idx_layer<n_checked_detectors; idx_layer++) {
    propagate(X,dir,theta,X_detectors,h_detectors,idx_layer);
    if( fabs(X_detectors[0][idx_layer])<Ldet/2 && fabs(X_detectors[1][idx_layer])<Ldet/2)
      checked_layer[idx_layer]=true;
  }

  accepted_event = true;   //Set initial status of this variable
  for(Int_t idx_layer=0; idx_layer<n_checked_detectors; idx_layer++) {
    accepted_event *= checked_layer[idx_layer];
    if(!accepted_event)
      break;
  }
}

void propagate(Double_t X[],Double_t dir[],Double_t theta,Double_t X_detectors[][2],const Double_t h_detectors[],Int_t n_layer) {
  Double_t dz = fabs(X[2]-h_detectors[n_layer]);
  Double_t r = dz/dir[2];

  X_detectors[0][n_layer] = X[0] + r*dir[0];
  X_detectors[1][n_layer] = X[1] + r*dir[1];
  X_detectors[2][n_layer] = h_detectors[n_layer];
}

Double_t Get_Telescope_Analysical_Acceptance(Double_t _l, Double_t _a1, Double_t _b1, Double_t _a2, Double_t _b2) {

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

  Double_t G = 0;
  
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

Double_t Get_Acceptance(Double_t ratio) {
  return TMath::Pi()*TMath::Power(Lgen,2)*ratio;
}
