
#include "MyHead.h"

void generate_coordinate(Double_t X[],TRandom3 *r_gen) {
  X[0] = r_gen->Uniform(-Lgen/2,Lgen/2);
  X[1] = r_gen->Uniform(-Lgen/2,Lgen/2);
  X[2] = -Lgen/2;
}

void generate_theta_phi(Double_t &theta,Double_t &phi,TRandom3 *r_gen) {
  phi = r_gen->Uniform(phi_min,phi_max);
  theta = acos(sqrt(r_gen->Uniform(cos2_theta_min,cos2_theta_max)));
  // theta = TMath::Pi() - acos(sqrt(r_gen->Uniform(cos2_theta_min,cos2_theta_max)));
}

void obtain_direction(Double_t theta,Double_t phi,Double_t dir[]) {
  //Obviously I cannot include also 'r' in this directions, because that will depend on the position of the layers
  dir[0] = sin(theta)*cos(phi);
  dir[1] = sin(theta)*sin(phi);
  dir[2] = cos(theta);
}
