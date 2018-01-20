#include "TMath.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TTree.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TNamed.h"
#include "TFile.h"
#include "TPaletteAxis.h"


void rapp(){
  TFile *i= new TFile("Isotropy_map232.root");
  TFile *ai= new TFile("Anisotropy_map232.root");
  TFile *s= new TFile("Shuffled_map232.root");
  TH2D *iso = (TH2D*)i->Get("mapi");
  TH2D *aniso = (TH2D*)ai->Get("map");
  TH2D *shuffle = (TH2D*)s->Get("Shuf");
  TH2D *cloneaiso = (TH2D*)(aniso->Clone("cloneaiso"));
  aniso->Divide(shuffle);
  cloneaiso->Divide(iso);

  aniso->SetXTitle("Glon");
  aniso->SetYTitle("Glat");
  cloneaiso->SetXTitle("Glon");
  cloneaiso->SetYTitle("Glat");

  TFile *fout1= new TFile ("232Anivs232iso.root","recreate");
  fout1->cd();
  cloneaiso->Write();
  fout1->Close();

  TFile *fout2= new TFile ("232Anivs232Shuf.root","recreate");
  fout2->cd();
  aniso->Write();
  fout2->Close();

  TCanvas *d= new TCanvas();
  d->cd();
  aniso->Draw("aitoff");
  TCanvas *c =new TCanvas();
  c->cd();
  c->SetName("pino");
  cloneaiso->Draw("aitoff");
}
