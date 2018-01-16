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

#include <fstream>

#include "various/orbitStruct.h"
#include "HEALpixDeNoAntri/anisotropy.C"

using namespace std;

/// Flag to denote that the function has sucessfully completed its job
#define NORMAL_END 0
#define NULL_VECTOR -1
#define EPS 1.e-12

int InvertDirection(AtPolarVect pvectorin, AtPolarVect *pvectorout);
//int atVectToPol(AtVect x, AtPolarVect *y);
int atVectToPol(double x[3], AtPolarVect *y);
//int FromSatelliteToCelestial(float ra[3], float  dec[3], AtVect vectorin, AtPolarVect *pvectorout );
int FromSatelliteToCelestial(float ra[3], float  dec[3], double vectorin[3], AtPolarVect *pvectorout );
int fromCelestialtoGalactic(double ra, double dec, double &l, double &b);
//int atPolToVect(AtPolarVect *x, AtVect y);
int atPolToVect(AtPolarVect *x, double y[3]);


int FillMapsForOneSecond(Int_t nev,
			 TH2D* mapstatinf,
			 TH2D* mapstatreal1,
			 TH2D* mapstatreal2,
			 Float_t satra[], Float_t satdec[],
			 bool poissonize=true,
			 bool anisotropize=false,
			 int multfact=10
			 );

void drawmaps(TString filename="./foutput_0103.root");
TVirtualPad* drawmap(TH2* htodraw, Bool_t klogz=false);
void draworbit(TString filename="./foutput_0103.root");
TVirtualPad* drawpull(TH2* hratio, TH2* hnum, TH2* hden);

void mappe(){

  //-------------- setting a deterministic seed -------

  gRandom->SetSeed(233445);

  //-------------- reading the SBI ----------
  
  TChain *tree= new TChain("SBItree");

  TString sbisubsample = "010";
  tree->Add(Form("SBI/%s*_SBI.root", sbisubsample.Data()));
            
  Float_t glat=0, glon=0;
  Float_t geola=0,geolo=0;
  Float_t satra[3], satdec[3];
  //Float_t zenra=0,zendec=0;
  UInt_t sec=0;
  UShort_t nev=0;
  bool good=true;
  tree->SetBranchAddress("second",&sec);
  tree->SetBranchAddress("goodsbi",&good);
  tree->SetBranchAddress("glon",&glon);
  tree->SetBranchAddress("glat",&glat);
  tree->SetBranchAddress("nev",&nev);
  //  tree->SetBranchAddress("ra_zenith",&zenra);
  //  tree->SetBranchAddress("dec_zenith",&zendec);
  tree->SetBranchAddress("ra_scx",&satra[0]);
  tree->SetBranchAddress("ra_scy",&satra[1]);
  tree->SetBranchAddress("ra_scz",&satra[2]);
  tree->SetBranchAddress("dec_scx",&satdec[0]);
  tree->SetBranchAddress("dec_scy",&satdec[1]);
  tree->SetBranchAddress("dec_scz",&satdec[2]);
  tree->SetBranchAddress("lat_geo",&geola);
  tree->SetBranchAddress("lon_geo",&geolo);
  tree->GetEntry(0);
  
  //----- costheta-flat binning ---------------
  int nbinlon=360;
  double lonbinmin=-180.0;
  double lonbinmax=180.0;
  
  int nbinlat=180;
  double latbinmin=-90.0;
  double latbinmax=90.0;

  double* binning;

  createbinning(nbinlat, latbinmin,latbinmax, binning, true);
  //----- histos ------------------------------

  TFile* foutput = new TFile(Form("foutput_%s.root", sbisubsample.Data()), "RECREATE");
  foutput->cd();
  
  TH2D* mapanistatinf = new TH2D("mapanistatinf","Mappa Anisotropa (statistica infinita); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries", nbinlon, lonbinmin, lonbinmax, nbinlat, binning);
  TH2D* mapanistatreal1 = new TH2D("mapanistatreal1","Mappa Anisotropa (statistica realistica 1); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries", nbinlon, lonbinmin, lonbinmax, nbinlat, binning);
  TH2D* mapanistatreal2 = new TH2D("mapanistatreal2","Mappa Anisotropa statistica realistica 2); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries", nbinlon, lonbinmin, lonbinmax, nbinlat, binning);
  TH2D* mapanishuf1 = new TH2D("mapanishuf1","Mappa Anisotropa Shufflata; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries", nbinlon, lonbinmin, lonbinmax, nbinlat, binning);
  TH2D* mapanishuf2 = new TH2D("mapanishuf2","Mappa Anisotropa Shufflata; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries", nbinlon, lonbinmin, lonbinmax, nbinlat, binning);
  
  TH2D* mapisostatinf = new TH2D("mapisostatinf","Mappa Isotropa (statistica infinita); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries", nbinlon, lonbinmin, lonbinmax, nbinlat, binning);
  TH2D* mapisostatreal1 = new TH2D("mapisostatreal1","Mappa Isotropa (statistica realistica 1); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries", nbinlon, lonbinmin, lonbinmax, nbinlat, binning);
  TH2D* mapisostatreal2 = new TH2D("mapisostatreal2","Mappa Isotropa statistica realistica 2); Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries",nbinlon, lonbinmin, lonbinmax, nbinlat, binning);
  TH2D* mapisoshuf1 = new TH2D("mapisoshuf1","Mappa Isotropa Shufflata; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries", nbinlon, lonbinmin, lonbinmax, nbinlat, binning);
  TH2D* mapisoshuf2 = new TH2D("mapisoshuf2","Mappa Isotropa Shufflata; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries", nbinlon, lonbinmin, lonbinmax, nbinlat, binning);
  
  TH2D* point = new TH2D("point","Pointing del Satellite; Galactic Longitude (#circ);  Galactic Latitude (#circ); Entries", nbinlon, lonbinmin, lonbinmax, nbinlat, binning);
  
  TH2D* orbit = new TH2D("orbit","Exposure; Geographic Longitude (#circ);  Geographic Latitude (#circ); Exposure (s)", 180, -180, 180, 90, -90.0, 90.0);
  TH2D* nevents = new TH2D("nevents","Events; Geographic Longitude (#circ);  Geographic Latitude (#circ); Entries", 180, -180, 180, 90, -90.0, 90.0);
  TH2D* rate = new TH2D("rate","Rate; Geographic Longitude (#circ);  Geographic Latitude (#circ); Rate (hz)", 180, -180, 180, 90, -90.0, 90.0);

  int entries = (int)tree->GetEntries();
  int perc=0;
  int treenumber=-1;
  for(int ientry=0; ientry<entries; ientry++){
    
    if (((1.0*ientry)/entries)>(perc*0.01)) {
      printf("%d%%\n", perc);
      perc++;
    }
    
    tree->GetEntry(ientry);   
    if (tree->GetTreeNumber()!=treenumber) {
      treenumber=tree->GetTreeNumber();
      printf("Tree %i opened\n", treenumber);
    }
    
    if((sec%100000)==0) continue;//there was a bug in the SBI production
    if(!good) continue;//good second

    // from [0,360] to [-180,180]
    if (glon>180) glon-=360;
    if (geolo>180) geolo-=360;
    
    point->Fill(glon, glat);
    orbit->Fill(geolo, geola);//once per second
    nevents->Fill(geolo, geola, nev);//once per event
    rate->Fill(geolo, geola, nev);//once per event and eventuallt divided by the exposure

    //change to radians
    for (int i=0; i<3; i++) {
      satra[i]*=TMath::DegToRad();
      satdec[i]*=TMath::DegToRad();
    }

    int kev=0;

    kev = FillMapsForOneSecond(nev, mapisostatinf, mapisostatreal1, mapisostatreal2, satra, satdec, true, false);
    FillMapsForOneSecond(kev, mapisostatinf, mapisoshuf1, mapisoshuf2, satra, satdec, false, false);

    kev = FillMapsForOneSecond(nev, mapanistatinf, mapanistatreal1, mapanistatreal2, satra, satdec, true, true);
    FillMapsForOneSecond(kev, mapanistatinf, mapanishuf1, mapanishuf2, satra, satdec, false, false);
  }

  rate->Divide(orbit);
  
  foutput->Write();
  foutput->Close();

  return;

}

int FillMapsForOneSecond(Int_t nev,
			 TH2D* mapstatinf,
			 TH2D* mapstatreal1,
			 TH2D* mapstatreal2,
			 Float_t satra[], Float_t satdec[],
			 bool poissonize,
			 bool anisotropize,
			 int multfact
			 ){

  /*
    EXPLANATIONS REGARDING THIS FUNCTION:

    The name of this function is FillMapsForOneSecond; the reason for that is that each entry of the TTree represents a second of DAMPE data taking.
    We have choosen to divide DAMPE acquisition data in single files of 100000 seconds, that means maybe 1.15 days of data.
    The reason for that is purely practical: in that way you can better ricognise at which file a single avent comes from; for example you take an electron triggered at 110653 seconds (from a reference time, obviously) and so you know that you have to take that file containing info from second 100000 to 200000, and so on.
    No physical reasons there are under this choose, just to avoid to make "retro-engeneering" to calculate the correct day of acquisition (needed, for example, if you create a file each single day).
    So we get information for each second, and all these informations were stored inside the TTree. Some informations, like teh position, is not precise at the second, so you can interpolate ! That's easy.
    
    nev is the number of events triggered in that second. We have an histo just for that into the SBI ROOT file.
    
    What I want to do in this function is just to fill all the different maps I will need.
    To do that I have just selected a second from the data TTree (in the function where FillMapsForOneSecond is called) and I now make a loop on the events triggered in that second (parameter passed to the function).
    Now I use a variable called multifact just to increment my statistic; that means for example that if multifact is 10 for each single event triggered by DAMPE I obtain 10 more events.
    If I want to build a realistic statistics map I just need one of all these 10 events (maybe the first or the second or what I want), but I take all of them if I want to build an infinite statistics map.
    Infact when we build the infinite statistic map we use the Fill method for each iteration of the loop on multifact, but I take just iteration number 0 and 1 for the real statistic maps.
    Why I build 2 real statistic maps ? The answer is easy, just to be able to confront them ! So I cam make the ratio beetween to maps with realistic statistic (maybe both isotrope or anisotrope or a mixed ratio).
    I could need, as I will, also a shuffled map. But I have just built it ! Let me explain: if I build 2 maps with realistic statistics, and if for each "event" (or for each sets of sub events (remember that the loop on multifact variable is just used to increment the statistic of the event) caming from the same physically triggered event) the local satellite variable costheta and phi are choosen randomly, each map is the shuffled map of the other ! That's simply geius !

    So I just have two shuffled maps (respectively) with just that function.

    When for each event (or sub-event) I choose costheta and phi, I can recognise a precise point of the sky (obviously knowing the satellite position in that moment, or better in that second).
    At that point (and just at that point) I can aky myself if the eletrons of positrons I am receiving are more or less from the mean, or in other words if an anisotropy is present.
    So I can "build" that anisotropy with the waith w (building the dipole that I prefer).
    But now w is a real number, and I cannot use that as the number of received particles. So I create a Poisson variable randomly froma poisson distribution with "w" as expected value. 
    Now that is a natural number and I can use that as the number or received particles.

    !!!!!!!!!!!!!!!!! Here the software has a problem: let me explain. 
    What we do in the real life is taking data and, after some time, we get the map of the observed particles from the satellite's surved sky.
    From that map we start working asking ourself if any kind of anisotropy is present.
    So we can found some sky's zones where an excess of particles are coming from (or maybe the countrary). But we don't know is that access was generated at the same "second", at the same istant; at big probability it was not ! 
    What we have done in our simulation, at the countray, is to simulate a some kind of anisotropy (excess or not, because we generate a poisson random variable were we only precise the expected value, so it could fluctuate upper or not the expected value itself) and we inject ALL of that particles AT THE SAME TIME and AT THE SAME POSITION into the map. That create a STATISTICAL CORRELATION beetween such events, and we can note that when we create the POOL histo !!!!

    However, let that bug as it was and go on with the simulation.

    Now I just do all that things with that simple function call:

    kev = FillMapsForOneSecond(nev, mapisostatinf, mapisostatreal1, mapisostatreal2, satra, satdec, true, false);
    
    kev is the number of poisson generated particles from the anisotropy weight w ! 
    Why I need to save it ? Just because if I want to compare two maps (with realistic statistic or not) I need two maps with the same number of particles !
    So now I call the same function (as before) but changing the anysotropy value (from true to false, or viceversa) and passing as number of events to generate the number of poisson particles generated before !

    You can simply see that from that function calls:

    kev = FillMapsForOneSecond(nev, mapisostatinf, mapisostatreal1, mapisostatreal2, satra, satdec, true, false);
    FillMapsForOneSecond(kev, mapisostatinf, mapisoshuf1, mapisoshuf2, satra, satdec, false, false);

    kev = FillMapsForOneSecond(nev, mapanistatinf, mapanistatreal1, mapanistatreal2, satra, satdec, true, true);
    FillMapsForOneSecond(kev, mapanistatinf, mapanishuf1, mapanishuf2, satra, satdec, false, false);

    We have talked about the bool just few moments ago; that's a statistical method that permits to understand how good is your simulation method to produce sky maps.

    Imagine you have a map of the real DAMPE data flight and one created by the simulation method; that maps are, at the and (as we just said) bidimensional histograms, so just bins filled with the number of events for each one of them.
    We have to compare the number of events of the real map with the simulated one; for example take one bin on the real map that counts 1034 particles, versus the corrisponding one on the simulated one with 1000 particles (the simulated one has 1000 particles for each bin). For another one we have 913 on the real one versus 1000 on the simulated map.
    What we build is the histogram of the difference between the really observed particles and the simulated particle (related to the same bin) divided the square root of the simulated particles for bin (that is, at the end, the standard deviation of the number of simulated events for each bin). Adding more events to both of the maps, we expect that the distribution shrinks, at limit becoming a delta distribution (a spike !!)
    !!!!!!!!!!!!!!!!!!!!!! That ratio behaves as a chi^2 !!!!!!!!!!!!!!!!!!
    What we expect on that distribution is that it respects a normal distribution, with mean 0 (if it has not mean 0 it means we have 2 maps with 2 different numbers of total events and that is not correct) and sigma 1.
    That's the shuffling method. If the fit of the obtained histo with the normal distribution is good we are ok, otherwise we have some problem with the simulated map.
    
    Remembering the bag mentioned before, what we obtain is a distribution that not fits well with the normal one, and so probably the cause is the correlation mentioned before in the simlation of the anisotropy with the w weight.


    Maybe we have to think about that and find a better simulation method.


   */
  
  static TFile* fileacc= new TFile("./ToyAcc/acceptance_final.root");
  static TH2D* acc = (TH2D*)fileacc->Get("Acceptance");
  
  Double_t b=0, l=0, w=0; 
  int p=0;
  
  AtVect vectorin;
  AtPolarVect vectorout, vectoroutinv;
  
  Double_t costheta=-999, phi=-999, ra=-999, dec=-999;

  if (multfact<3) {
    printf("Multiplication factor less than 3 has no meaning. Set to 3.\n");
    multfact=3;
  }
  
  int kev=0;
  
  for (int i=0; i<nev;i++){
    for (int mm=0; mm<multfact; mm++){
      
      acc->GetRandom2(costheta,phi);
      
      vectorin[0]=sin(acos(costheta))*cos(phi);
      vectorin[1]=sin(acos(costheta))*sin(phi);
      vectorin[2]=costheta;
      
      //from coordinates in the satellite frame (cosine directors) to RA and DEC in the celestial frame
      FromSatelliteToCelestial(satra, satdec, vectorin, &vectorout);
      
      //for some reasons we've to invert the direction (exiting <-> entering)
      InvertDirection(vectorout, &vectoroutinv);    
      ra=vectoroutinv.lon;
      dec=vectoroutinv.lat;
      //printf("%f %f\n",ra,dec);
      
      //change to degrees
      ra=ra*TMath::RadToDeg();
      dec=dec*TMath::RadToDeg();
      
      //from celestial co'os to galactic ones
      fromCelestialtoGalactic(ra,dec,l,b);

      if (isnan(l)) {
	mm--;//to not to "loose" this event
	continue;
      }
      
      // ----- weight to create anisotropic maps --------
      if (anisotropize) {
	w = 1.0 + TMath::Sin(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());//E-W dipole
	//	w = 1.0 + TMath::Cos(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());//F-B dipole
	//	w = 1.0 + TMath::Cos(b*TMath::DegToRad());//N-S dipole
	if (w<0 || w>2) printf("Problem: w=%f\n", w);
	
      }
      else {
	w=1;
      }

      p=-1;      
      if (poissonize) {
	while (p<0 || p>999) {
	  p = gRandom->Poisson(w);
	  if (p<0 || p>999) {
	    printf ("Problem p=%d (w=%f, %E) (l,b = %f,%f)\n", p, w, w, l, b);
	  }
	}
      }
      else {
	p=1;
      }
      
      // from [0,360] to [-180,180]
      if (l>180.0) l-=360.0;
      
      mapstatinf->Fill(l, b, p);
      if (mm==0) {
	mapstatreal1->Fill(l, b, p);
	kev+=p;
      }
      if (mm==1) mapstatreal2->Fill(l, b, p);      
    } //end loop on multfact !!!
  } //end loop on events !! 
  
  return kev;
}

//------------------- drawing functions ---------------

void drawmaps(TString filename){

  // Double_t __r__[7]    = {0.0, 0.0, 0.0, 0.1, 1.0, 0.8, 0.5};
  // Double_t __g__[7]    = {0.0, 0.2, 0.6, 0.7, 0.6, 0.2, 0.0};
  // Double_t __b__[7]    = {0.1, 0.8, 1.0, 0.1, 0.0, 0.0, 0.0};

  // const int nRGBs=4;
  // Double_t __r__[nRGBs]    = {1.0, 0.5,  1.0,  0.84};
  // Double_t __g__[nRGBs]    = {1.0, 0.0,  0.6,  1.0};
  // Double_t __b__[nRGBs]    = {1.0, 0.0,  0.0,  0.2};
  // Double_t __stop__[nRGBs] = {0.0, 0.001, 0.5, 1.0};
  // TColor::CreateGradientColorTable(nRGBs, __stop__, __r__, __g__, __b__, 255);
  
  // Add a little style
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(8);

  TH2* htodrawclone = NULL;
  
  TFile* finput = new TFile(filename.Data());

  // KEY: TH2D	mapanistatinf;1	Mappa Anisotropa (statistica infinita)
  // KEY: TH2D	mapanistatreal1;1	Mappa Anisotropa (statistica realistica 1)
  // KEY: TH2D	mapanistatreal2;1	Mappa Anisotropa statistica realistica 2)
  // KEY: TH2D	mapanishuf1;1	Mappa Anisotropa Shufflata
  // KEY: TH2D	mapanishuf2;1	Mappa Anisotropa Shufflata
  // KEY: TH2D	mapisostatinf;1	Mappa Isotropa (statistica infinita)
  // KEY: TH2D	mapisostatreal1;1	Mappa Isotropa (statistica realistica 1)
  // KEY: TH2D	mapisostatreal2;1	Mappa Isotropa statistica realistica 2)
  // KEY: TH2D	mapisoshuf1;1	Mappa Isotropa Shufflata
  // KEY: TH2D	mapisoshuf2;1	Mappa Isotropa Shufflata

  Bool_t klogz = false;

  TH2* point = (TH2D*)(finput->Get("point"));
  point->SetMinimum(1);
  drawmap(point, true);
  
  TH2* mapisostatinf = (TH2D*)(finput->Get("mapisostatinf"));
  drawmap(mapisostatinf, klogz);
  
  TH2*  mapisostatreal1 = (TH2D*)(finput->Get("mapisostatreal1"));
  drawmap(mapisostatreal1, klogz);

  TH2* mapisostatreal2 = (TH2D*)(finput->Get("mapisostatreal2"));
  drawmap(mapisostatreal2, klogz);

  TH2* mapisoshuf1 = (TH2D*)(finput->Get("mapisoshuf1"));
  drawmap(mapisoshuf1, klogz);
  
  TH2* mapanistatinf = (TH2D*)(finput->Get("mapanistatinf"));
  drawmap(mapanistatinf, klogz);
  
  TH2* mapanistatreal1 = (TH2D*)(finput->Get("mapanistatreal1"));
  drawmap(mapanistatreal1, klogz);

  TH2* mapanishuf1 = (TH2D*)(finput->Get("mapanishuf1"));
  drawmap(mapanishuf1, klogz);

  TH2* mapisooverisostatreal1 = (TH2*)(mapisostatreal1->Clone("mapisooverisostatreal1"));
  mapisooverisostatreal1->Divide(mapisostatreal2);
  mapisooverisostatreal1->SetMinimum(0);
  mapisooverisostatreal1->SetMaximum(2);
  mapisooverisostatreal1->SetTitle("Ratio Isotropa/Isotropa (statistica realistica 1/2)");
  mapisooverisostatreal1->GetZaxis()->SetTitle("Ratio");
  drawmap(mapisooverisostatreal1, klogz);
  drawpull(mapisooverisostatreal1, mapisostatreal1, mapisostatreal2);
  
  TH2* mapiso1overshuf1 = (TH2*)(mapisostatreal1->Clone("mapiso1overshuf1"));
  mapiso1overshuf1->Divide(mapisoshuf1);
  mapiso1overshuf1->SetMinimum(0);
  mapiso1overshuf1->SetMaximum(2);
  mapiso1overshuf1->SetTitle("Ratio Isotropa/Shuffled (statistica realistica 1)");
  mapiso1overshuf1->GetZaxis()->SetTitle("Ratio");
  drawmap(mapiso1overshuf1, klogz);
  drawpull(mapiso1overshuf1, mapisostatreal1, mapisoshuf1);
  
  TH2* mapanioverisostatinf = (TH2*)(mapanistatinf->Clone("mapanioverisostatinf"));
  mapanioverisostatinf->Divide(mapisostatinf);
  mapanioverisostatinf->SetMinimum(0);
  mapanioverisostatinf->SetMaximum(2);
  mapanioverisostatinf->SetTitle("Ratio Anisotropa/Isotropa (statistica infinita)");
  mapanioverisostatinf->GetZaxis()->SetTitle("Ratio");
  drawmap(mapanioverisostatinf, klogz);
  drawpull(mapanioverisostatinf, mapanistatinf, mapisostatinf);

  TH2* mapanioverisostatreal1 = (TH2*)(mapanistatreal1->Clone("mapanioverisostatreal1"));
  mapanioverisostatreal1->Divide(mapisostatreal1);
  mapanioverisostatreal1->SetMinimum(0);
  mapanioverisostatreal1->SetMaximum(2);
  mapanioverisostatreal1->SetTitle("Ratio Anisotropa/Isotropa (statistica realistica 1)");
  mapanioverisostatreal1->GetZaxis()->SetTitle("Ratio");
  drawmap(mapanioverisostatreal1, klogz);
  drawpull(mapanioverisostatreal1, mapanistatreal1, mapisostatreal1);

  TH2* mapani1overshuf1 = (TH2*)(mapanistatreal1->Clone("mapani1overshuf1"));
  mapani1overshuf1->Divide(mapanishuf1);
  mapani1overshuf1->SetMinimum(0);
  mapani1overshuf1->SetMaximum(2);
  mapani1overshuf1->SetTitle("Ratio Anisotropa/Shuffled (statistica realistica 1)");
  mapani1overshuf1->GetZaxis()->SetTitle("Ratio");
  drawmap(mapani1overshuf1, klogz);
  drawpull(mapani1overshuf1, mapanistatreal1, mapanishuf1);

  /*
  TH2D* hnum = new TH2D("hnum", "hnum", 100, 0, 100, 100, 0, 100);
  TH2D* hden = new TH2D("hden", "hden", 100, 0, 100, 100, 0, 100);
  // for (int ii=0; ii<10000000; ii++) {
  //   double x = gRandom->Uniform(0, 100);
  //   double y = gRandom->Uniform(0, 100);
  //   hnum->Fill(x, y, gRandom->Poisson(64));
  //   x = gRandom->Uniform(0, 100);
  //   y = gRandom->Uniform(0, 100);
  //   hden->Fill(x, y, gRandom->Poisson(64));
  // }
  for (int xx=0; xx<100; xx++) {
    for (int yy=0; yy<100; yy++) {
      for (int ii=0; ii<10000; ii++) {
	hnum->Fill(xx, yy, gRandom->Poisson(64));
	hden->Fill(xx, yy, gRandom->Poisson(64));
      }
    }
  }

  TH2* hratio = (TH2*)(hnum->Clone("hratio"));
  hratio->Divide(hden);
  drawpull(hratio, hnum, hden);
  */
  
  return;
}

TVirtualPad* drawpull(TH2* hratio, TH2* hnum, TH2* hden){

  double pullmin = 0;
  double pullmax = 0;

  ////////////////
  /*
    This forst two loops are just needed because I wanto to obtain max and min value for the pull histo.
    These two values will be needed to proprly set the pool's histos max and min value for the X axis.


    The second loop just fills poll's histo ! So obviously I need a couple of loops.
  */
  ///////////////
  
  for (int xx=1; xx<=hratio->GetNbinsX(); xx++) {
    for (int yy=1; yy<=hratio->GetNbinsY(); yy++) {
      double yrec = hnum->GetBinContent(xx, yy);
      double yexp = hden->GetBinContent(xx, yy);
      double sigma = sqrt(yexp+yrec);
      if (sigma>0) {
       	double pull = (yrec-yexp)/sigma;
	if (pull<pullmin) pullmin=pull;
	if (pull>pullmax) pullmax=pull;
      }
    }
  }
  
  TH1D* hpull = new TH1D(Form("%s_pull", hratio->GetName()), Form("Pull di %s", hratio->GetTitle()), 4*3*5*7, 1.2*pullmin, 1.2*pullmax);
  hpull->GetXaxis()->SetTitle("(y_{rec}-y_{ref})/(#sqrt{#sigma^{2}_{rec}+#sigma^{2}_{ref}})");
  hpull->GetYaxis()->SetTitle("Entries");

  //  TGraph* gr = new TGraph();
  
  for (int xx=1; xx<=hratio->GetNbinsX(); xx++) {
    for (int yy=1; yy<=hratio->GetNbinsY(); yy++) {
      double yrec = hnum->GetBinContent(xx, yy);
      double yexp = hden->GetBinContent(xx, yy);
      double sigma = sqrt(yexp+yrec);
      if (yrec<0 || yexp<0) printf("%f %f\n", yrec, yexp);
      if (sigma>0) {
	double pull = (yrec-yexp)/sigma;
	hpull->Fill(pull);
	// if ((fabs(pull))<0.1*2.0*20.0/1000.0) {
	//   //	  printf("%f %E \n", fabs(pull), fabs(pull));
	//   gr->SetPoint(gr->GetN(), hratio->GetXaxis()->GetBinCenter(xx), hratio->GetYaxis()->GetBinCenter(yy));
	// }
      }
    }
  }

  //  gr->Draw("P");
  
  double canvxsize=1024;
  double canvysize=768;
  
  TCanvas* c_pull = new TCanvas(Form("Pull_%s", hratio->GetName()), hratio->GetTitle(), 1.1*canvxsize, canvysize);
  c_pull->cd();
  hpull->Draw();
  
  return c_pull;  
}

TVirtualPad* drawmap(TH2* htodraw, Bool_t klogz){

  // Add a little style
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(8);

  TH2* htodrawclone = (TH2D*)(htodraw->Clone(Form("%s_clone", htodraw->GetName())));
  //  printf("%p\n", htodrawclone);
  
  double canvxsize=1024;
  double canvysize=768;
  
  htodrawclone->GetXaxis()->SetLabelSize(0.0);
  htodrawclone->GetXaxis()->SetTickSize(0.0);
  htodrawclone->GetXaxis()->SetTitleOffset(0.5);
  htodrawclone->GetYaxis()->SetLabelSize(0.0);
  htodrawclone->GetYaxis()->SetTickSize(0.0);
  htodrawclone->GetYaxis()->SetTitleOffset(0.5);
  
  TCanvas* c_aitoff = new TCanvas(Form("Aitoff_%s", htodraw->GetName()), htodraw->GetTitle(), 1.1*canvxsize, canvysize);
  c_aitoff->cd();
  TPad* pad = new TPad("pad","", 0, 0, 1.0/1.1-0.075, 1.0, 0, 0);
  pad->SetNumber(1);
  pad->Draw();
  c_aitoff->cd(1);
  //  gPad->Draw();
  // gPad->SetTopMargin(0);
  // gPad->SetLeftMargin(0);
  gPad->SetRightMargin(0.01);
  //  gPad->SetBottomMargin(0);
  htodrawclone->Draw("aitoff");
  c_aitoff->cd();
  TPaletteAxis* pal = new TPaletteAxis(1.0/1.1-0.075, 0.1, 1.0/1.1-0.075+0.05, 0.90, htodrawclone);
  //  pal->SetLabelSize(0.02);
  pal->Draw();
  //--------------------------------------------------
  
  drawgrid(c_aitoff->cd(1), canvxsize, canvysize);
  
  TVirtualPad* pad2return = c_aitoff->cd(1);
  pad2return->SetLogz(klogz);

  c_aitoff->cd();
  c_aitoff->SetLogz(klogz);
  
  return pad2return;
}

void draworbit(TString filename) {

  // Add a little style
  gStyle->SetOptStat(0);
  
  TH2F* earth = new TH2F("earth", "Earth", 180, -180, 180, 180, -90, 90);
  earth->GetXaxis()->SetTitle("Geographic Longitude (#circ)");
  earth->GetYaxis()->SetTitle("Geographic Latitude (#circ)");
  
  TString dat="various/earth.dat";
  ifstream in;
  in.open(dat.Data());
  Float_t x,y;
  while (1) {
    in >> x >> y;
    if (!in.good()) break;
    earth->Fill(x,y, 1);
  }
  in.close();

  TFile* finput = new TFile(filename.Data());

  TH2D* orbit = (TH2D*)finput->Get("orbit");
  TH2D* nevents = (TH2D*)finput->Get("nevents");
  TH2D* rate = (TH2D*)finput->Get("rate");

  orbit->SetMarkerStyle(1);
  orbit->SetMarkerColor(kRed+2);
  //  orbit->SetMarkerSize(0.2);
  
  TCanvas* c_orbit = new TCanvas("C_Orbit", "Orbit");
  earth->DrawCopy();
  orbit->Draw("same");//this is beautiful only if done with few orbits

  earth->SetMarkerStyle(1);
  
  TCanvas* c_rate = new TCanvas("C_Rate", "Rate");
  rate->Draw("colz");
  earth->DrawCopy("same");

  return;
}

//------ from DAMPESW/OrbitTools-----------

int fromCelestialtoGalactic(double ra, double dec, double &l, double &b){
  
  double ragc=192.85948,decgc=27.12825,lcp=122.932;
 
  double ragcr=ragc*TMath::DegToRad();
  double decgcr=decgc*TMath::DegToRad();
  double lcpr=lcp*TMath::DegToRad();
  double decr=dec*TMath::DegToRad();
  double rar=ra*TMath::DegToRad();
  
  double br=asin(sin(decgcr)*sin(decr)+cos(decgcr)*cos(decr)*(cos(rar-ragcr)));
 
  //Calculate sin(theta) and cos(theta)

  double sin_t=(cos(decr)*sin(rar-ragcr)/cos(br));
  double cos_t=((cos(decgcr)*sin(decr)-sin(decgcr)*cos(decr)*cos(rar-ragcr))/cos(br));

  //Use signs of sin_t and cos_t to determine which quadrant t lies in

  double t=0;
  if((sin_t>=0)&&(cos_t>=0)){
    t=asin(sin_t);
  }
  else if((sin_t>=0)&&(cos_t<0)){
    t=TMath::Pi()-asin(sin_t);
  }
  else if((sin_t<0)&&(cos_t<0)){
    t=TMath::Pi()-asin(sin_t);
  }
  else if((sin_t<0)&&(cos_t>=0)){
    // t=2.*TMath::Pi()-asin(sin_t);
    t=asin(sin_t);
  }
  //     Convert t to l
  double lr = lcpr-t;
  //     Ensure l lies in the range 0-2pi 

  if(lr<0) {lr=lr+(2.*TMath::Pi());}
  if(lr>(2.*TMath::Pi())){lr=lr-(2.*TMath::Pi());}

  b = br*TMath::RadToDeg();
  l = lr*TMath::RadToDeg();
  return 0;
}

//int FromSatelliteToCelestial(float ra[3], float dec[3], AtVect vectorin, AtPolarVect *pvectorout ){
int FromSatelliteToCelestial(float ra[3], float dec[3], double vectorin[3], AtPolarVect *pvectorout ){
  
  // ra and dec in radiants
  // ra and dec of the three satellite axies
  
  //double pmom1[3];
  //double pmom[3];
 
  double ux1[3];
  double uy1[3];
  double uz1[3];
  double rax = ra[0];
  double ray = ra[1];
  double raz = ra[2];
  double decx = dec[0];
  double decy = dec[1];
  double decz = dec[2];
  //double pvpart_px = vectorin[0];
  //double pvpart_py = vectorin[1];
  //double pvpart_pz = vectorin[2];
  //double gst,theta0;
  AtVect vectorout;
  // (1) Take direction cosines of instrument axes in the celestial reference frame

  ux1[0] = cos(decx)*cos(rax);
  ux1[1] = cos(decx)*sin(rax);
  ux1[2] = sin(decx);
  uy1[0] = cos(decy)*cos(ray);
  uy1[1] = cos(decy)*sin(ray);
  uy1[2] = sin(decy);
  uz1[0] = cos(decz)*cos(raz);
  uz1[1] = cos(decz)*sin(raz);
  uz1[2] = sin(decz);

  // (2) Rotate the particle momentum in the celestial reference frame
  for (int i=0; i<3; i++){
    vectorout[i] = vectorin[0] * ux1[i] + vectorin[1] * uy1[i] + vectorin[2] * uz1[i];
  }
  atVectToPol(vectorout, pvectorout);
  
  return 0;
}

//int atVectToPol(AtVect x, AtPolarVect *y) {
int atVectToPol(double x[3], AtPolarVect *y) {
  double  norm01, c, s;
  norm01 = x[0]*x[0] + x[1]*x[1];
  if ( (y->r = sqrt( norm01 + x[2]*x[2] )) == 0.0 ) { 
    y->lon = y->lat = 0.0;
    return (NULL_VECTOR); 
  }
  norm01 = sqrt( norm01 );
  y->lat = asin( x[2]/y->r );
  c = x[0]/norm01;  s = x[1]/norm01;
  if (norm01 < EPS) { y->lon = 0.0; }
  else if ( fabs( s ) < EPS ) { y->lon = (1.0 - c/fabs(c))*TMath::Pi()/2.0; }
  else { y->lon = atan((1.0-c)/s)*2.0; }
  while( y->lon >= 2.*TMath::Pi() ) { y->lon -= 2.*TMath::Pi(); }
  while( y->lon < 0.0 )    { y->lon += 2.*TMath::Pi(); }
  return (NORMAL_END);
}

int InvertDirection(AtPolarVect pvectorin, AtPolarVect *pvectorout){
  AtVect vectorin;
  AtVect vectorout;
  atPolToVect(&pvectorin, vectorin);
  for(int i=0; i<3;i++){
    vectorout[i]=-vectorin[i];
  }
  atVectToPol(vectorout,pvectorout);
  return 0;
}

// Converting the coordinate from Polar to Cartesian (POLREC)
// x: input
// y: result (output)
//int atPolToVect(AtPolarVect *x, AtVect y){
int atPolToVect(AtPolarVect *x, double y[3]){
   y[0] = (x->r)*cos(x->lat)*cos(x->lon);
   y[1] = (x->r)*cos(x->lat)*sin(x->lon);
   y[2] = (x->r)*sin(x->lat);
  return (NORMAL_END);
}


