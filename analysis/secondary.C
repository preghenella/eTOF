/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, and plot simple quantities such as the jet pt and the di-electron invariant
mass.

root -l examples/Example1.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

//------------------------------------------------------------------------------

// define TOF (barrel/forward)
double bTOFr    = eTOF_radius; // [cm]
double bTOFz    = eTOF_length; // [cm]
double fTOFrin  = 50.;         // [cm]
double fTOFrout = eTOF_radius; // [cm]
double fTOFz    = eTOF_length; // [cm]
double TOFreso  = eTOF_sigmat; // [ns]

// define masses
double mass[5] = {0.00051099891, 
		  0.10565800, 
		  0.13957000, 
		  0.49367700, 
		  0.93827200};
// define PDG codes
int pdgcode[5] = {11, 13, 211, 321, 2212};

// function prototypes
std::vector<std::vector<Track *>> vertexer(TClonesArray *branchTrack);
std::map<Track *, int> TOFpid(std::vector<Track *> &tracks);
void invariantMass(TH2 *h, 
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

// main function
void secondary(const char *inputFile = "delphes.root")
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");

  // Book histograms
  TH2 *hMass_Dplus = new TH2F("hMass_Dplus", "", 500, 0., 5., 100, 0., 10.);
  TH2 *hMass_Ds = new TH2F("hMass_Ds", "", 500, 0., 5., 100, 0., 10.);
  TH2 *hMass_Lambdac = new TH2F("hMass_Lambdac", "", 500, 0., 5., 100, 0., 10.);
  TH2 *hMassT_Dplus = new TH2F("hMassT_Dplus", "", 500, 0., 5., 100, 0., 10.);
  TH2 *hMassT_Ds = new TH2F("hMassT_Ds", "", 500, 0., 5., 100, 0., 10.);
  TH2 *hMassT_Lambdac = new TH2F("hMassT_Lambdac", "", 500, 0., 5., 100, 0., 10.);
  TH2 *hMassI_Dplus = new TH2F("hMassI_Dplus", "", 500, 0., 5., 100, 0., 10.);
  TH2 *hMassI_Ds = new TH2F("hMassI_Ds", "", 500, 0., 5., 100, 0., 10.);
  TH2 *hMassI_Lambdac = new TH2F("hMassI_Lambdac", "", 500, 0., 5., 100, 0., 10.);
    
  // Loop over all events
  for (Int_t ientry = 0; ientry < numberOfEntries; ++ientry) {
    
    // Load selected branches with data from specified event
    treeReader->ReadEntry(ientry);

    // run vertexer
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

#if 0
      // define pion, kaon and proton candidates
      std::vector<Track *> kaonsS, protonsS, pionsS;
      std::vector<Track *> kaonsO, protonsO, Opions;
      std::vector<Track *> kaonsT, protonsT, pionsT;
      std::vector<Track *> kaonsI, protonsI, pionsI;
      for (auto &track : vertex) {
	
	// kaon candidate, has opposite charge than vertex
	if (track->Charge != vertexCharge) {
	  kaonsO.push_back(track);
	  if (!pid.count(track)) // no PID available
	    kaonsT.push_back(track);
	  else if (pid[track] & (1 << 3)) // PID available (kaon)
	    kaonsT.push_back(track);
	  if (TMath::Abs(track->PID) == 321) // ideal PID
	    kaonsI.push_back(track);
	}
	// proton and pion candidates, same charge as vertex
	else {
	  pionsS.push_back(track);
	  kaonsS.push_back(track);
	  protonsS.push_back(track);
	  if (!pid.count(track)) { // no PID available
	    protonsT.push_back(track);
	    pionsT.push_back(track);
	  }
	  else { // PID available
	    if (pid[track] & (1 << 4)) // proton
	      protonsT.push_back(track);
	    if (pid[track] & (1 << 2)) // pion
	      pionsT.push_back(track);
	  }
	  if (TMath::Abs(track->PID) == 2212) // ideal PID
	    protonsI.push_back(track);
	  if (TMath::Abs(track->PID) == 211) // ideal PID
	    pionsI.push_back(track);
	}
      }    
      
      invariantMass_Dplus(hMass_Dplus, kaonsO, pionsS, -2., 2.);
      invariantMass_Ds(hMass_Ds, kaonsO, pionsS, kaonsS, -2., 2.);
      invariantMass_Lambdac(hMass_Lambdac, kaonsO, pionsS, protonsS, -2., 2.);
      //      invariantMass(hMass_tof, pionsT, kaonsT, protonsT, 1.);
      //      invariantMass(hMass_ideal, pionsI, kaonsI, protonsI, 1.);
#endif

    }
  } // loop over events

  TFile *fout = TFile::Open("secondary.root", "RECREATE");
  hMass_Dplus->Write();
  hMass_Ds->Write();
  hMass_Lambdac->Write();
  hMassT_Dplus->Write();
  hMassT_Ds->Write();
  hMassT_Lambdac->Write();
  hMassI_Dplus->Write();
  hMassI_Ds->Write();
  hMassI_Lambdac->Write();
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

    // TOF matching
    Double_t R = TMath::Sqrt(track->XOuter * track->XOuter +
			     track->YOuter * track->YOuter) * 0.1; // [cm]
    Double_t Z = track->ZOuter * 0.1; // [cm]
    bool hasTOF = false;
    if (TMath::Abs(R - bTOFr) < 0.001 && TMath::Abs(Z) < bTOFr) hasTOF = true;
    else if (R < fTOFrout && R > fTOFrin && TMath::Abs(TMath::Abs(Z) - fTOFz) < 0.001) hasTOF = true;
    else continue;

    Double_t tof = track->TOuter * 1.e9; // [ns]
    Double_t L = track->L * 0.1;         // [cm]
    Double_t p = track->P;               // [GeV/c]
    Double_t ep = track->ErrorP;         // [GeV/c]
    Double_t c = 29.9792458;             // [cm/ns]

    // R+HACK (to avoid momentum smearing and texp resolution)
    GenParticle *particle = (GenParticle *) track->Particle.GetObject();
    p = particle->P;

    // loop over mass hypoteses
    pid[track] = 0;
    for (Int_t ipart = 0; ipart < 5; ++ipart) {
      Double_t texp = L / c / p * TMath::Sqrt(mass[ipart] * mass[ipart] + p * p);
      Double_t etexp = L / c * mass[ipart] * mass[ipart] / p / p / TMath::Sqrt(mass[ipart] * mass[ipart] + p * p) * ep;
      
      Double_t sigma = TMath::Sqrt(etexp * etexp + TOFreso * TOFreso);
      Double_t delta = tof - texp;
      Double_t nsigma = delta / sigma;
      if (TMath::Abs(nsigma) < 3.) 
	pid[track] |= 1 << ipart;
    }
  } 

  return pid;
}

/************************************************************
 * invariant mass
 ************************************************************/

void
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
  if (LV.Rapidity() < ymin || LV.Rapidity() > ymax) return;
  h->Fill(LV.Mag(), LV.Pt());

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
	invariantMass(h, { {pion1, 0.139570},  {pion2, 0.139570}, {kaon, 0.493677} }, ymin, ymax);
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
	invariantMass(h, { {pion, 0.139570},  {kaonO, 0.493677}, {kaonS, 0.493677} }, ymin, ymax);
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
	invariantMass(h, { {pion, 0.139570}, {kaon, 0.493677}, {proton, 0.93827200} }, ymin, ymax);
     }
    }
  }
}

