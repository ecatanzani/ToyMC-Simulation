
#include "MyHead.h"

void fill_maps_for_one_second(ofstream &out_file,UShort_t n_ev,Int_t &n_poisson_events,TH2D* map_inf_stat,TH2D* map_real_stat_1,TH2D* map_real_stat_2,Float_t sat_ra[],Float_t sat_dec[],Bool_t poissonize,Bool_t anisotropize) {
  
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
    
    fill_maps_for_one_second(n_ev,n_poisson_events,iso_map_inf_stat,iso_map_real_stat_1,iso_map_real_stat_1,sat_ra,sat_dec,true,false);

    "n_poisson_events" is the number of poisson generated particles from the anisotropy weight w ! 
    Why I need to save it ? Just because if I want to compare two maps (with realistic statistic or not) I need two maps with the same number of particles !
    So now I call the same function (as before) but changing the anysotropy value (from true to false, or viceversa) and passing as number of events to generate the number of poisson particles generated before !
    You can simply see that from that function calls:
    
    // Filling isotropic maps -------                                                                                                                                                                             

    fill_maps_for_one_second(n_ev,n_poisson_events,iso_map_inf_stat,iso_map_real_stat_1,iso_map_real_stat_1,sat_ra,sat_dec,true,false);
    fill_maps_for_one_second(n_poisson_events,n_poisson_events,iso_map_inf_stat,iso_map_shuf_1,iso_map_shuf_2,sat_ra,sat_dec,false,false);

    // Filling anisotropic maps -------                                                                                                                                                                            
    fill_maps_for_one_second(n_ev,n_poisson_events,ani_map_stat_inf,ani_map_real_stat_1,ani_map_real_stat_2,sat_ra,sat_dec,true,true);
    fill_maps_for_one_second(n_poisson_events,n_poisson_events,ani_map_stat_inf,ani_map_shuf_1,ani_map_shuf_2,sat_ra,sat_dec,false,false);

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

  
  Double_t b=0,l=0,anisotropy_eccess=0;
  Double_t costheta=-999, phi=-999, ra=-999, dec=-999;
  Int_t p_events=0;

  AtVect vector_in;                                                          //This is just a common array with 3 pomponents
  AtPolarVect vector_out, vector_out_inv;                                    //This, at the countrary, is a structure for the 3-dimensional vector in polar coordinates.
                                                                             //That structure includes the radial component, latitude and longitude

  /*
    ATTENCTION: explanation of the variable: physical meaning

    vector_in                -> This is the cosine directors of the coming events in the local satellite frame
    vector_out               -> Is the polar vector that describes the position in the sky (radial component, right ascension and celestial declination)
    vector_out_inv           -> Is the output vector, but inverted !!!

   */ 
  
  if(multiplication_factor<3) {
    cout<<"\nMultiplication factor less than 3 has no meaning. Set to 3.\n";
    out_file<<"\nMultiplication factor less than 3 has no meaning. Set to 3"<<endl;
    multiplication_factor=3;
  }
  
  static TFile* acc_file= new TFile(acceptance_final_plot.c_str());          //Opening the final acceptance file
  if(acc_file->IsZombie()) {
    cout<<"\n\nError opening DAMPE acceptance TFile. Prorgram finished \n\n";
    out_file<<"\n\nError opening DAMPE acceptance TFile. Prorgram finished \n\n";
    exit(-1);
  }
  
  static TH2D* acc = (TH2D*)acc_file->Get("Acceptance");                     //Get the final acceptance histo

  if(n_ev==n_poisson_events) {
    n_ev=n_poisson_events;
    n_poisson_events = 0;
  }
  else
    n_poisson_events=0;

  for(Int_t tr_ev=0; tr_ev<n_ev; tr_ev++)  {                                 //loop on the triggered events
    for(Int_t m_fact=0; m_fact<multiplication_factor; m_fact++) {            //loop on the multiplication factor

      acc->GetRandom2(costheta,phi);

      //---------------- Filling the array in the reference frame of the satellite
      
      vector_in[0]=sin(acos(costheta))*cos(phi);
      vector_in[1]=sin(acos(costheta))*sin(phi);
      vector_in[2]=costheta;

      //---------------- From coordinates in the satellite frame (cosine directors) to RA and DEC in the celestial frame
      // Infact what is needed is the position of the satellite (stored in the array sat_ra and sat_dec) in right ascension and celestial declination and the local coorinate of the entering events (as cosine directors). In taht way we'll obtain the right ascension and celestial declination of the events coming from the "sources" in the sky map !!!

      from_satellite_to_celestial(sat_ra,sat_dec,vector_in,vector_out);

      //for some reasons we've to invert the direction (exiting <-> entering)
      
      invert_AtPolarVect_direction(vector_out,vector_out_inv);
      ra=vector_out_inv.lon;
      dec=vector_out_inv.lat;

      //change to degrees
      
      ra=ra*TMath::RadToDeg();
      dec=dec*TMath::RadToDeg();

      //From celestial to galactic coordinate
      
      from_celestial_to_galactic(ra,dec,l,b);
      
      if (isnan(l)) {
	m_fact--;         //In that way I don't "loose" this event
	continue;
      }


      // ----- weight to create anisotropic maps --------    CHOOSE THE DIPOLE YOUO WANT......
      
      if (anisotropize) {
        anisotropy_eccess = 1.0 + TMath::Sin(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());                     //E-W dipole
	// anisotropy_eccess = 1.0 + TMath::Cos(l*TMath::DegToRad())*TMath::Cos(b*TMath::DegToRad());                  //F-B dipole
        // anisotropy_eccess = 1.0 + TMath::Cos(b*TMath::DegToRad());                                                  //N-S dipole
	if (anisotropy_eccess<0 || anisotropy_eccess>2)
	  cout<<"Problem: anisotropy_eccess="<<anisotropy_eccess<<"\n";
      }
      else
        anisotropy_eccess=1;


      // ----- choose if you want to poissonize or not --------  
      
      p_events=-1;
      if (poissonize) {
        while (p_events<0 || p_events>999) {
          p_events = gRandom->Poisson(anisotropy_eccess);
          if (p_events<0 || p_events>999) {
            cout<<"\nProblem number poisson events = "<<p_events<<" (anisotropy_eccess = "<<anisotropy_eccess<<" ) (galactic latitude = "<<l<<" - galactic longitude = "<<b<<endl;
	    out_file<<"\nProblem number poisson events = "<<p_events<<" (anisotropy_eccess = "<<anisotropy_eccess<<" ) (galactic latitude = "<<l<<" - galactic longitude = "<<b<<endl;
	  }
        }
      }
      else {
        p_events = 1;
      }

      
      // transform latitude's interval from [0,360] to [-180,180]                                                                                                                                                                                
      if (l>180.0) l-=360.0;

      map_inf_stat->Fill(l,b,p_events);
      if (m_fact==0) {
        map_real_stat_1->Fill(l,b,p_events);
        n_poisson_events+=p_events;
      }
      if (m_fact==1)
	map_real_stat_2->Fill(l,b,p_events);
      
    } //end loop on multiplication factor variable
  } //end loop on n_ev variable
  
}

void from_satellite_to_celestial(Float_t ra[],Float_t dec[],double vectorin[],AtPolarVect &vector_out) {

  // ra (right ascension) and dec (celestial declination) are in radiants
  // ra and dec of the three satellite axies

  Float_t ux1[3];
  Float_t uy1[3];
  Float_t uz1[3];
  Float_t rax = ra[0];
  Float_t ray = ra[1];
  Float_t raz = ra[2];
  Float_t decx = dec[0];
  Float_t decy = dec[1];
  Float_t decz = dec[2];

  AtVect tmp_vector_out;
  
  // (1) Take director cosine of instrument axes in the celestial reference frame
  
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
  
  for (Int_t idx=0; idx<3; idx++)
    tmp_vector_out[idx] = vectorin[0] * ux1[idx] + vectorin[1] * uy1[idx] + vectorin[2] * uz1[idx];
  
  AtVect_To_AtPolarVect(tmp_vector_out,vector_out);
  
}

void AtVect_To_AtPolarVect(double in_vector[],AtPolarVect &vector_out) {
  double  norm01, c, s;

  norm01 = TMath::Power(in_vector[0],2) + TMath::Power(in_vector[1],2);

  if ( ( vector_out.r = sqrt( norm01 + TMath::Power(in_vector[2],2) ) ) == 0.0 )
    vector_out.lon = vector_out.lat = 0.0;
  
  norm01 = sqrt(norm01);
  vector_out.lat = asin(in_vector[2]/vector_out.r);
  c = in_vector[0]/norm01;
  s = in_vector[1]/norm01;

  if(norm01 < EPS)
    vector_out.lon = 0.0;
  else if (fabs(s) < EPS) {
    vector_out.lon = (1.0 - c/fabs(c))*TMath::Pi()/2.0;
  }
  else
    vector_out.lon = atan((1.0-c)/s)*2.0;

  while(vector_out.lon >= 2.*TMath::Pi())
    vector_out.lon -= 2.*TMath::Pi();
  
  while(vector_out.lon < 0.0)
    vector_out.lon += 2.*TMath::Pi();
  
}

void invert_AtPolarVect_direction(AtPolarVect vector_out,AtPolarVect &vector_out_inv) {

  AtVect in_array;
  AtVect out_array;

  AtPolarVect_to_vector(vector_out,in_array);
  for(Int_t idx=0; idx<3; idx++)
    out_array[idx]=-in_array[idx];
  AtVect_To_AtPolarVect(out_array,vector_out_inv);

}

void AtPolarVect_to_vector(AtPolarVect &input_polar,double out_array[]) {
  // Converting the coordinate from Polar to Cartesian (POLREC)

  out_array[0] = (input_polar.r)*cos(input_polar.lat)*cos(input_polar.lon);
  out_array[1] = (input_polar.r)*cos(input_polar.lat)*sin(input_polar.lon);
  out_array[2] = (input_polar.r)*sin(input_polar.lat);

}

void from_celestial_to_galactic(Double_t ra,Double_t dec,Double_t &l,Double_t &b) {

  Double_t ragc=192.85948,decgc=27.12825,lcp=122.932;

  Double_t ragcr=ragc*TMath::DegToRad();
  Double_t decgcr=decgc*TMath::DegToRad();
  Double_t lcpr=lcp*TMath::DegToRad();
  Double_t decr=dec*TMath::DegToRad();
  Double_t rar=ra*TMath::DegToRad();

  Double_t br=asin(sin(decgcr)*sin(decr)+cos(decgcr)*cos(decr)*(cos(rar-ragcr)));

  //Calculate sin(theta) and cos(theta)
  
  Double_t sin_t=(cos(decr)*sin(rar-ragcr)/cos(br));
  Double_t cos_t=((cos(decgcr)*sin(decr)-sin(decgcr)*cos(decr)*cos(rar-ragcr))/cos(br));

  //Use signs of sin_t and cos_t to determine which quadrant t lies in
  
  Double_t t=0;
  if((sin_t>=0)&&(cos_t>=0))
    t=asin(sin_t);
  else if((sin_t>=0)&&(cos_t<0))
    t=TMath::Pi()-asin(sin_t);
  else if((sin_t<0)&&(cos_t<0))
    t=TMath::Pi()-asin(sin_t);
  else if((sin_t<0)&&(cos_t>=0))
    // t=2.*TMath::Pi()-asin(sin_t);                                                                                                                                                                              
    t=asin(sin_t);
  
  // Convert t to l                                                                                                                                                                                            

  Double_t lr = lcpr-t;

  // Ensure l lies in the range 0-2pi                                                                                                                                                                          

  if(lr<0) {lr=lr+(2.*TMath::Pi());}
  if(lr>(2.*TMath::Pi())){lr=lr-(2.*TMath::Pi());}

  b = br*TMath::RadToDeg();
  l = lr*TMath::RadToDeg();

}
