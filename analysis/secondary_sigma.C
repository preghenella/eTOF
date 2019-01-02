// define PDG codes
int pdgcode[5] = {11, 13, 211, 321, 2212};

// function prototypes
std::vector<std::vector<Track *>> vertexer(TClonesArray *branchTrack);
std::map<Track *, int> TOFpid(std::vector<Track *> &tracks);
TLorentzVector invariantMass(TH2 *h, 
		   std::vector<std::pair<Track *, double>> trackmasspairs,
		   double ymin, double ymax);
void invariantMass_Dplus(TH2 *h, 
			 std::vector<Track *> &kaonsO, 
			 std::vector<Track *> &pionsS,
			 double ymin, double ymax);
void invariantMass_Ds(TH2 *h, 
		      std::vector<Track *> &kaonsO, 
		      std::vector<Track *> &pions,
		      std::vector<Track *> &kaonsS, 
		      double ymin, double ymax);
void invariantMass_Lambdac(TH2 *h, 
			   std::vector<Track *> &kaonsO, 
			   std::vector<Track *> &pionsS,
			   std::vector<Track *> &protonsS, 
			   double ymin, double ymax);
void invariantMass_Sigmac(TH2 *h, 
			  std::vector<Track *> &kaonsO, 
			  std::vector<Track *> &pionsS,
			  std::vector<Track *> &protonsS, 
			  std::vector<Track *> &pionsP,
			  double ymin, double ymax);

// main function
void secondary_sigma(const char *inputFile = "delphes.root")
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");

  // Book histograms
  TH1 *hNevents = new TH1F("hNevents", "", 1, 0, 1);
  TH1 *hNruns = new TH1F("hNruns", "", 1, 0, 1);
  TH1 *hXsec = new TH1F("hXsec", "", 1, 0, 1); 
  //
  TH2 *hMass_Dplus = new TH2F("hMass_Dplus", "", 1000, 1., 3.5, 100, 0., 10.);
  TH2 *hMass_Ds = new TH2F("hMass_Ds", "", 1000, 1., 3.5, 100, 0., 10.);
  TH2 *hMass_Lambdac = new TH2F("hMass_Lambdac", "", 1000, 1., 3.5, 100, 0., 10.);
  TH2 *hMass_Sigmac = new TH2F("hMass_Sigmac", "", 1000, 2.4, 2.7, 100, 0., 10.);
  TH2 *hMass_Sigmac0 = new TH2F("hMass_Sigmac0", "", 1000, 2.4, 2.7, 100, 0., 10.);
  TH2 *hMassT_Dplus = new TH2F("hMassT_Dplus", "", 1000, 1., 3.5, 100, 0., 10.);
  TH2 *hMassT_Ds = new TH2F("hMassT_Ds", "", 1000, 1., 3.5, 100, 0., 10.);
  TH2 *hMassT_Lambdac = new TH2F("hMassT_Lambdac", "", 1000, 1., 3.5, 100, 0., 10.);
  TH2 *hMassT_Sigmac = new TH2F("hMassT_Sigmac", "", 1000, 2.4, 2.7, 100, 0., 10.);
  TH2 *hMassT_Sigmac0 = new TH2F("hMassT_Sigmac0", "", 1000, 2.4, 2.7, 100, 0., 10.);
  TH2 *hMassI_Dplus = new TH2F("hMassI_Dplus", "", 1000, 1., 3.5, 100, 0., 10.);
  TH2 *hMassI_Ds = new TH2F("hMassI_Ds", "", 1000, 1., 3.5, 100, 0., 10.);
  TH2 *hMassI_Lambdac = new TH2F("hMassI_Lambdac", "", 1000, 1., 3.5, 100, 0., 10.);
  TH2 *hMassI_Sigmac = new TH2F("hMassI_Sigmac", "", 1000, 2.4, 2.7, 100, 0., 10.);
  TH2 *hMassI_Sigmac0 = new TH2F("hMassI_Sigmac0", "", 1000, 2.4, 2.7, 100, 0., 10.);
  
  // Loop over all events
  for (Int_t ientry = 0; ientry < numberOfEntries; ++ientry) {
    
    // Load selected branches with data from specified event
    treeReader->ReadEntry(ientry);

    // select primaries
    std::vector<Track *> primaries;
    // loop over tracks
    for (Int_t itrack = 0; itrack < branchTrack->GetEntries(); ++itrack) {
      // get tracks
      Track *track = (Track *)branchTrack->At(itrack);
      // skip non-primary tracks
      if (track->X != 0. || track->Y != 0. || track->Z != 0.) continue; 
      primaries.push_back(track);
    } // loop over tracks
    
    // perform TOF pid
    auto pid = TOFpid(primaries);

    std::vector<Track *> Pidentified[5][2], PidentifiedT[5][2], PidentifiedI[5][2];
    for (auto &track : primaries) {
      for (int ipart = 0; ipart < 5; ++ipart) {
	Pidentified[ipart][track->Charge > 0 ? 0 : 1].push_back(track);
	if (!pid.count(track)) // no PID available
	  PidentifiedT[ipart][track->Charge > 0 ? 0 : 1].push_back(track);
	else if (pid[track] & (1 << ipart)) // PID compatible
	  PidentifiedT[ipart][track->Charge > 0 ? 0 : 1].push_back(track);
	if (TMath::Abs(track->PID) == pdgcode[ipart]) // ideal PID
	  PidentifiedI[ipart][track->Charge > 0 ? 0 : 1].push_back(track);
      }
    }

    // run secondary vertexer
    auto vertices = vertexer(branchTrack);

    // loop over vertices
    for (auto &vertex : vertices) {
      
      // check vertex size
      if (vertex.size() != 3) continue;

      // perform TOF pid
      auto pid = TOFpid(vertex);
            
      // define vertex charge
      int vertexCharge = 0;
      for (auto &track : vertex) vertexCharge += track->Charge;

      
      std::vector<Track *> identified[5][2], identifiedT[5][2], identifiedI[5][2];
      for (auto &track : vertex) {
	for (int ipart = 0; ipart < 5; ++ipart) {
	  identified[ipart][track->Charge > 0 ? 0 : 1].push_back(track);
	  if (!pid.count(track)) // no PID available
	    identifiedT[ipart][track->Charge > 0 ? 0 : 1].push_back(track);
	  else if (pid[track] & (1 << ipart)) // PID compatible
	    identifiedT[ipart][track->Charge > 0 ? 0 : 1].push_back(track);
	  if (TMath::Abs(track->PID) == pdgcode[ipart]) // ideal PID
	    identifiedI[ipart][track->Charge > 0 ? 0 : 1].push_back(track);
	}
      }

      // no PID
      invariantMass_Dplus(hMass_Dplus, 
			  identified[3][vertexCharge > 0 ? 1 : 0], 
			  identified[2][vertexCharge > 0 ? 0 : 1], 
			  -2., 2.);
      invariantMass_Ds(hMass_Ds, 
		       identified[3][vertexCharge > 0 ? 1 : 0], 
		       identified[2][vertexCharge > 0 ? 0 : 1],
		       identified[3][vertexCharge > 0 ? 0 : 1],
		       -2., 2.);
      invariantMass_Lambdac(hMass_Lambdac,
			    identified[3][vertexCharge > 0 ? 1 : 0],
			    identified[2][vertexCharge > 0 ? 0 : 1],
			    identified[4][vertexCharge > 0 ? 0 : 1], 
			    -2., 2.);
      invariantMass_Sigmac(hMass_Sigmac,
			   identified[3][vertexCharge > 0 ? 1 : 0],
			   identified[2][vertexCharge > 0 ? 0 : 1],
			   identified[4][vertexCharge > 0 ? 0 : 1], 
			   Pidentified[2][vertexCharge > 0 ? 0 : 1], 
			    -2., 2.);
      invariantMass_Sigmac(hMass_Sigmac0,
			   identified[3][vertexCharge > 0 ? 1 : 0],
			   identified[2][vertexCharge > 0 ? 0 : 1],
			   identified[4][vertexCharge > 0 ? 0 : 1], 
			   Pidentified[2][vertexCharge > 0 ? 1 : 0], 
			    -2., 2.);

      // TOF PID if available
      invariantMass_Dplus(hMassT_Dplus, 
			  identifiedT[3][vertexCharge > 0 ? 1 : 0], 
			  identifiedT[2][vertexCharge > 0 ? 0 : 1], 
			  -2., 2.);
      invariantMass_Ds(hMassT_Ds, 
		       identifiedT[3][vertexCharge > 0 ? 1 : 0], 
		       identifiedT[2][vertexCharge > 0 ? 0 : 1],
		       identifiedT[3][vertexCharge > 0 ? 0 : 1],
		       -2., 2.);
      invariantMass_Lambdac(hMassT_Lambdac,
			    identifiedT[3][vertexCharge > 0 ? 1 : 0],
			    identifiedT[2][vertexCharge > 0 ? 0 : 1],
			    identifiedT[4][vertexCharge > 0 ? 0 : 1], 
			    -2., 2.);
      invariantMass_Sigmac(hMassT_Sigmac,
			   identifiedT[3][vertexCharge > 0 ? 1 : 0],
			   identifiedT[2][vertexCharge > 0 ? 0 : 1],
			   identifiedT[4][vertexCharge > 0 ? 0 : 1], 
			   PidentifiedT[2][vertexCharge > 0 ? 0 : 1], 
			    -2., 2.);
      invariantMass_Sigmac(hMassT_Sigmac0,
			   identifiedT[3][vertexCharge > 0 ? 1 : 0],
			   identifiedT[2][vertexCharge > 0 ? 0 : 1],
			   identifiedT[4][vertexCharge > 0 ? 0 : 1], 
			   PidentifiedT[2][vertexCharge > 0 ? 1 : 0], 
			    -2., 2.);

      // ideal PID
      invariantMass_Dplus(hMassI_Dplus, 
			  identifiedI[3][vertexCharge > 0 ? 1 : 0], 
			  identifiedI[2][vertexCharge > 0 ? 0 : 1], 
			  -2., 2.);
      invariantMass_Ds(hMassI_Ds, 
		       identifiedI[3][vertexCharge > 0 ? 1 : 0], 
		       identifiedI[2][vertexCharge > 0 ? 0 : 1],
		       identifiedI[3][vertexCharge > 0 ? 0 : 1],
		       -2., 2.);
      invariantMass_Lambdac(hMassI_Lambdac,
			    identifiedI[3][vertexCharge > 0 ? 1 : 0],
			    identifiedI[2][vertexCharge > 0 ? 0 : 1],
			    identifiedI[4][vertexCharge > 0 ? 0 : 1], 
			    -2., 2.);
      invariantMass_Sigmac(hMassI_Sigmac,
			   identifiedI[3][vertexCharge > 0 ? 1 : 0],
			   identifiedI[2][vertexCharge > 0 ? 0 : 1],
			   identifiedI[4][vertexCharge > 0 ? 0 : 1], 
			   PidentifiedI[2][vertexCharge > 0 ? 0 : 1], 
			    -2., 2.);
      invariantMass_Sigmac(hMassI_Sigmac0,
			   identifiedI[3][vertexCharge > 0 ? 1 : 0],
			   identifiedI[2][vertexCharge > 0 ? 0 : 1],
			   identifiedI[4][vertexCharge > 0 ? 0 : 1], 
			   PidentifiedI[2][vertexCharge > 0 ? 1 : 0], 
			    -2., 2.);

    }

  } // loop over events

  // fill run information
  HepMCEvent *event = (HepMCEvent *)branchEvent->At(0);
  hNruns->Fill(0);
  hNevents->SetBinContent(1, numberOfEntries);
  hXsec->SetBinContent(1, event->CrossSection);
  hXsec->SetBinError(1, event->CrossSectionError);

  TFile *fout = TFile::Open("secondary_sigma.root", "RECREATE");
  //
  hNevents->Write();
  hNruns->Write();
  hXsec->Write();
  //
  hMass_Dplus->Write();
  hMass_Ds->Write();
  hMass_Lambdac->Write();
  hMass_Sigmac->Write();
  hMass_Sigmac0->Write();
  hMassT_Dplus->Write();
  hMassT_Ds->Write();
  hMassT_Lambdac->Write();
  hMassT_Sigmac->Write();
  hMassT_Sigmac0->Write();
  hMassI_Dplus->Write();
  hMassI_Ds->Write();
  hMassI_Lambdac->Write();
  hMassI_Sigmac->Write();
  hMassI_Sigmac0->Write();
  fout->Close();

}

/************************************************************
 * vertexer
 ************************************************************/

std::vector<std::vector<Track *>>
vertexer(TClonesArray *branchTrack)
{
  
  // define seeds
  std::vector<std::pair<Track *, bool>> seeds;
  // loop over tracks
  for (Int_t itrack = 0; itrack < branchTrack->GetEntries(); ++itrack) {
    // get tracks
    Track *track = (Track *)branchTrack->At(itrack);
    // skip primary tracks
    if (track->X == 0. && track->Y == 0. && track->Z == 0.) continue; 
    // create and add unused seed
    std::pair<Track *, bool> seed = {track, false};
    seeds.push_back(seed);
  } // loop over tracks
  
  // define secondary vertices
  std::vector<std::vector<Track *>> vertices;
  // loop over seeds
  for (auto &seed : seeds) {
    // skip if already used
    if (seed.second) continue;
    double seedX = seed.first->X;
    double seedY = seed.first->Y;
    double seedZ = seed.first->Z;
    // create candidate vertex
    std::vector<Track *> vertex;
    vertex.push_back(seed.first);
    // mark seed as used
    seed.second = true;
    
    // loop over seeds
    for (auto &track : seeds) {
      // skip if already used
      if (track.second) continue;
      // check position
      if (TMath::Abs(track.first->X - seedX) != 0.) continue;
      if (TMath::Abs(track.first->Y - seedY) != 0.) continue;
      if (TMath::Abs(track.first->Z - seedZ) != 0.) continue;
      // add track to candidate vertex
      vertex.push_back(track.first);
      // mark seed as used
      track.second = true;
    }
    // check vertex multiplicity
    if (vertex.size() < 2) {
      // release seed
      seed.second = false;
      continue;
    }
    // add vertex
    vertices.push_back(vertex);
  }
  
  return vertices;
}

/************************************************************
 * TOF pid
 ************************************************************/
			  
std::map<Track *, int>
TOFpid(std::vector<Track *> &tracks)
{

  // define PID map
  std::map<Track *, int> pid;
  // loop over tracks
  for (auto &track : tracks) {

    /** eTOF cuts **/
    double deltat[5], nsigma[5];
    int etof = eTOF(track, deltat, nsigma);
    if (etof < 0) continue;

    // loop over mass hypoteses
    pid[track] = 0;
    for (Int_t ipart = 0; ipart < 5; ++ipart)
      if (TMath::Abs(nsigma[ipart]) < 3.) 
	pid[track] |= 1 << ipart;
    
  } 

  return pid;
}

/************************************************************
 * invariant mass
 ************************************************************/

TLorentzVector
invariantMass(TH2 *h, 
	      std::vector<std::pair<Track *, double>> trackmasspairs,
	      double ymin, double ymax)
{
  TLorentzVector LV, lv;
  for (auto &trackmasspair : trackmasspairs) {
    auto &track = trackmasspair.first;
    auto mass = trackmasspair.second;
    lv.SetPtEtaPhiM(track->PT, track->Eta, track->Phi, mass);
    LV += lv;
  }
  if (LV.Rapidity() < ymin || LV.Rapidity() > ymax) return LV;
  if (h) h->Fill(LV.Mag(), LV.Pt());
  return LV;
}

void
invariantMass_Dplus(TH2 *h, 
		    std::vector<Track *> &kaons, 
		    std::vector<Track *> &pions,
		    double ymin, double ymax)
{
  if (kaons.size() == 0 || pions.size() < 2) return;
  for (int ikaon = 0; ikaon < kaons.size(); ++ikaon) {
    auto &kaon = kaons[ikaon];
    for (int ipion1 = 0; ipion1 < pions.size(); ++ipion1) {
      auto &pion1 = pions[ipion1];
      if (pion1 == kaon) continue;
      for (int ipion2 = ipion1 + 1; ipion2 < pions.size(); ++ipion2) {
	auto &pion2 = pions[ipion2];
	if (pion2 == kaon || pion2 == pion1) continue;
	invariantMass(h, { {pion1, 0.139570},  {pion2, 0.139570}, {kaon, 0.493600} }, ymin, ymax);
      }
    }
  }
}

void
invariantMass_Ds(TH2 *h, 
		 std::vector<Track *> &kaonsO, 
		 std::vector<Track *> &pions,
		 std::vector<Track *> &kaonsS, 
		 double ymin, double ymax)
{
  if (kaonsO.size() == 0 || pions.size() == 0 || kaonsS.size() == 0) return;
  for (int ikaonO = 0; ikaonO < kaonsO.size(); ++ikaonO) {
    auto &kaonO = kaonsO[ikaonO];
    for (int ipion = 0; ipion < pions.size(); ++ipion) {
      auto &pion = pions[ipion];
      if (pion == kaonO) continue;
      for (int ikaonS = 0; ikaonS < kaonsS.size(); ++ikaonS) {
	auto &kaonS = kaonsS[ikaonS];
	if (kaonS == kaonO || kaonS == pion) continue;
	invariantMass(h, { {pion, 0.139570},  {kaonO, 0.493600}, {kaonS, 0.493600} }, ymin, ymax);
      }
    }
  }
}

void
invariantMass_Lambdac(TH2 *h, 
		      std::vector<Track *> &kaons, 
		      std::vector<Track *> &pions,
		      std::vector<Track *> &protons, 
		      double ymin, double ymax)
{
  if (kaons.size() == 0 || pions.size() == 0 || protons.size() == 0) return;
  for (auto &kaon : kaons) {
    for (auto &proton : protons) {
      if (proton == kaon) continue;
      for (auto &pion : pions) {
	if (pion == proton || pion == kaon) continue;
	// proton candidate has higher momentum than pion
	//	if (proton->P < pion->P) continue;
	invariantMass(h, { {pion, 0.139570}, {kaon, 0.493600}, {proton, 0.938270} }, ymin, ymax);
     }
    }
  }
}

void
invariantMass_Sigmac(TH2 *h, 
		      std::vector<Track *> &kaons, 
		      std::vector<Track *> &pions,
		      std::vector<Track *> &protons, 
		      std::vector<Track *> &pionsP,
		      double ymin, double ymax)
{
  if (kaons.size() == 0 || pions.size() == 0 || protons.size() == 0 || pionsP.size() == 0) return;
  for (auto &kaon : kaons) {
    for (auto &proton : protons) {
      if (proton == kaon) continue;
      for (auto &pion : pions) {
	if (pion == proton || pion == kaon) continue;
	// proton candidate has higher momentum than pion
	//	if (proton->P < pion->P) continue;
	auto LV = invariantMass(nullptr, { {pion, 0.139570}, {kaon, 0.493600}, {proton, 0.938270} }, ymin, ymax);
	if (LV.M() < 2.26 || LV.M() > 2.32) continue;
	Track lambdac;
	lambdac.PT = LV.Pt();
	lambdac.Eta = LV.Eta();
	lambdac.Phi = LV.Phi();
	for (auto &pionP : pionsP) {
	  if (pionP == pion || pionP == kaon || pionP == proton) continue;
	  invariantMass(h, { {&lambdac, 2.28646}, {pionP, 0.139570} }, ymin, ymax);
	}
      }
    }
  }
}

