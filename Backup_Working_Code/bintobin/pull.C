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

void pull(){
  TFile *i= new TFile("Isotropy_map232.root");
  TFile *ai= new TFile("Anisotropy_map232.root");
  TFile *s= new TFile("Shuffled_map232.root");
  TH2D *iso = (TH2D*)i->Get("mapi");
  TH2D *aniso = (TH2D*)ai->Get("map");
  TH2D *shuffle = (TH2D*)s->Get("Shuf");
  TH2D *pullmapiso = (TH2D*)(aniso->Clone("pullmapiso"));
  TH2D *pullmapshuf = (TH2D*)(aniso->Clone("pullmapshuf"));
  pullmapiso->Add(iso,-1);
  pullmapshuf->Add(shuffle,-1);
  TH1D *pulliso = new TH1D("pulliso","pulliso",1000,-100,100);
  TH1D *pullshuf= new TH1D("pullshuf","pullshuf",1000,-100,100);
  Double_t deltaiso=1,deltashuf=1;
  for(int ii=1;ii<=pullmapiso->GetNbinsX();ii++){
    for(int jj=1;jj<=pullmapiso->GetNbinsY();jj++){
      if (iso->GetBinContent(ii,jj)>=5){
	deltaiso=pullmapiso->GetBinContent(ii,jj)/sqrt(iso->GetBinContent(ii,jj));
	pulliso->Fill(deltaiso);
	}
      if(shuffle->GetBinContent(ii,jj)>=5){
      deltashuf=pullmapshuf->GetBinContent(ii,jj)/sqrt(shuffle->GetBinContent(ii,jj));
      pullshuf->Fill(deltashuf);
       }
    }
  }
  TCanvas *d= new TCanvas();
  d->cd();
  pulliso->Draw();
  TCanvas *c =new TCanvas();
  c->cd();
  pullshuf->Draw();

}
