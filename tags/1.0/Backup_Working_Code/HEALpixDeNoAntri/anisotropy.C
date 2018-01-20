#include "TVirtualPad.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TPaletteAxis.h"
#include "TGraph.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TText.h"
#include "TLatex.h"

void drawgrid(TVirtualPad* pad, double canvxsize, double canvysize, bool draw_label=false, int n_coo=6, int lat_label=12);
void createbinning(int nbinlat, double latbinmin, double latbinmax, double*& binning, bool cosscale=false);

void anisotropy(){
  
  // Add a little style
  gStyle->SetOptStat(0);

  gStyle->SetNumberContours(8);

  //------------------ histos --------------------
  int nbinlon=180;
  double lonbinmin=-180.0;
  double lonbinmax=180.0;
  
  int nbinlat=180;
  double latbinmin=-90.0;
  double latbinmax=90.0;

  double* binning;

  bool cosscale=true;
  createbinning(nbinlat, latbinmin,latbinmax, binning, cosscale);

  TH2F* ha = new TH2F("ha","Aitoff; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries", nbinlon, lonbinmin, lonbinmax, nbinlat, binning);
  //--------------------------------------------------

  //------------------ fill the histo ----------------
  // 0=isotropic
  // 1=east-west dipole
  // 2=forward-backward dipole
  // 3=north-south dipole
  int anisotropytype=2;
  
  double x;
  double y, cosy;
  double w;
  for (int ii=0; ii<10000000; ii++) {
    x = gRandom->Uniform(-180.0, 180.0);
    if (cosscale) {
      cosy = gRandom->Uniform(-1.0, 1.0);
      y = TMath::RadToDeg()*TMath::ACos(cosy) - 90.0;
    }
    else {//this is not really isotropical, but appears flat in a flat-theta histo
      y = gRandom->Uniform(-90.0, 90.0);
    }
    if (anisotropytype==0) w=1;
    else if (anisotropytype==1) w = 1.0 + TMath::Sin(x*TMath::DegToRad())*TMath::Cos(y*TMath::DegToRad());//E-W dipole
    else if (anisotropytype==2) w = 1.0 + TMath::Cos(x*TMath::DegToRad())*TMath::Cos(y*TMath::DegToRad());//F-B dipole
    else if (anisotropytype==3) w = 1.0 + TMath::Cos(y*TMath::DegToRad());//N-S dipole

    ha->Fill(x, y, w);
  }
  //--------------------------------------------------
  
  //----------------- draw the histo -----------------

  double canvxsize=1024;
  double canvysize=768;
  
  TCanvas* c_colz = new TCanvas("ColZ", "ColZ", canvxsize, canvysize);
  c_colz->cd();
  ha->DrawCopy("colz");

  ha->GetXaxis()->SetLabelSize(0.0);
  ha->GetXaxis()->SetTickSize(0.0);
  ha->GetXaxis()->SetTitleOffset(0.5);
  ha->GetYaxis()->SetLabelSize(0.0);
  ha->GetYaxis()->SetTickSize(0.0);
  ha->GetYaxis()->SetTitleOffset(0.5);
  TCanvas* c_aitoff = new TCanvas("Aitoff", "Aitoff", 1.1*canvxsize, canvysize);
  c_aitoff->cd();
  TPad* pad = new TPad("pad","", 0, 0, 1.0/1.1-0.075, 1.0, 0, 0);
  pad->SetNumber(1);
  pad->Draw();
  c_aitoff->cd(1);
  gPad->Draw();
  // gPad->SetTopMargin(0);
  // gPad->SetLeftMargin(0);
  gPad->SetRightMargin(0.01);
  //  gPad->SetBottomMargin(0);
  ha->Draw("Aitoff");
  //  ha->Draw("Mercator");
  //  ha->Draw("Sinusoidal");
  //  ha->Draw("Parabolic");
  c_aitoff->cd();
  TPaletteAxis* pal = new TPaletteAxis(1.0/1.1-0.075,0.1,1.0/1.1-0.075+0.05,0.90,ha);
  //  pal->SetLabelSize(0.02);
  pal->Draw();
  //--------------------------------------------------

  drawgrid(c_aitoff->cd(1), canvxsize, canvysize);
  
  c_aitoff->cd();
  
  return;
}

void drawgrid(TVirtualPad* pad, double canvxsize, double canvysize, bool draw_label, int n_coo, int lat_label){

  //  pad->Dump();
  double convx=((1.0/(canvxsize*1.02))*(10.0/3.0));
  double convy=((1.0/(canvysize*0.85*0.8))*(10.0/3.0));
  
  float la, lo, _x, _y, z;
  
  const int Nl = 19; // Number of drawn latitudes
  const int NL = 19; // Number of drawn longitudes
  int       M  = 90;
  
  TGraph  *latitudes[Nl];
  TGraph  *longitudes[NL];
  
  for (int j=0;j<Nl;++j) {
    latitudes[j]=new TGraph();
    latitudes[j]->SetMarkerSize(0.3);
    la = -90+180/(Nl-1)*j;
    double _M = M;
    if (fabs(la)>50.0) {
      _M = M*3.0/5.0;
    }
    if (fabs(la)>60.0) {
      _M = M*2.0/5.0;
    }
    if (fabs(la)>70.0) {
      _M = M*1.0/5.0;
    }
    for (int i=0;i<_M+1;++i) {
      lo = -180+360/_M*i;
      z  = sqrt(1+cos(la*TMath::DegToRad())*cos(lo*TMath::DegToRad()/2));
      _x  = 180*cos(la*TMath::DegToRad())*sin(lo*TMath::DegToRad()/2)/z;
      _y  = 90*sin(la*TMath::DegToRad())/z;
      latitudes[j]->SetPoint(i,_x*convx,_y*convy);
    }
  } 
  
  for (int j=0;j<NL;++j) {
    longitudes[j]=new TGraph();
    longitudes[j]->SetMarkerSize(0.3);
    lo = -180+360/(NL-1)*j;
    for (int i=0;i<M+1;++i) {
      la = -90+180/M*i;
      if (fabs(la)<85.0) {
	z  = sqrt(1+cos(la*TMath::DegToRad())*cos(lo*TMath::DegToRad()/2));
	_x  = 180*cos(la*TMath::DegToRad())*sin(lo*TMath::DegToRad()/2)/z;
	_y  = 90*sin(la*TMath::DegToRad())/z;
	longitudes[j]->SetPoint(i,_x*convx,_y*convy);
      }
    }
  }
  
  bool drawgrid=true;
  if (drawgrid) {
    // Draw the grid
    pad->cd();
    for (int j=0;j<Nl;++j) {
      if (j!=0 && j!=(Nl-1)) latitudes[j]->Draw("P");
    }
    for (int j=0;j<NL;++j) {
      longitudes[j]->Draw("P");
    }
  }


  // draw label on map

  if (draw_label==true){
    char coo_name[255];
    double* vec_xx = (double*)latitudes[lat_label]->GetX();
    int Npoints = (int)latitudes[lat_label]->GetN();
    double delta_coo = (double)Npoints/(double)n_coo; 
    double delta_coo_for_label = 360./(double)n_coo; 
    double coo_for_label = -180.;
    int num_coo = 0;
    for (int nn=0; nn<=n_coo; nn++){
      double y = latitudes[lat_label]->Eval(vec_xx[num_coo]);
      double x = vec_xx[num_coo];
      sprintf(coo_name,"%.0f.",coo_for_label);
      TText* tex = new TLatex(x,y,coo_name);
      tex->SetTextFont(40);
      tex->SetTextSize(0.03);
      tex->Draw();
      num_coo = num_coo+delta_coo;
      coo_for_label = coo_for_label+delta_coo_for_label;
    }
  }

  return;
}

void createbinning(int nbinlat, double latbinmin, double latbinmax, double*& binning, bool cosscale){
  
  binning = new double[nbinlat+1];
  
  if (cosscale==false) {
    binning[0]=latbinmin;
    //    printf("%d) %f\n", 0, binning[0]);
    for (int ii=1; ii<=nbinlat; ii++) {
      binning[ii] = binning[ii-1] + (latbinmax-latbinmin)/nbinlat;
      //      printf("%d) %f\n", ii, binning[ii]);
    }
  }
  else {
    double binningshift = (TMath::Cos(TMath::DegToRad()*(latbinmax-90.0))-TMath::Cos(TMath::DegToRad()*(latbinmin-90.0)))/nbinlat;
    binning[0]=latbinmin;
    //    printf("%d) %f\n", 0, binning[0]);
    for (int ii=1; ii<=nbinlat; ii++) {
      double binningmenouno = TMath::DegToRad()*(binning[ii-1]-90.0);
      binning[ii] = 90.0-(TMath::RadToDeg()*TMath::ACos(TMath::Cos(binningmenouno)+binningshift));
      //      printf("%f = TMath::ACos(TMath::Cos(%f)+%f) = TMath::ACos(%f+%f)\n", binning[ii], TMath::RadToDeg()*binningmenouno, binningshift, TMath::Cos(binningmenouno), binningshift);
      //      printf("%d) %f\n", ii, binning[ii]);
      //      gSystem->Sleep(1000);
    }
  }

  return;
}
