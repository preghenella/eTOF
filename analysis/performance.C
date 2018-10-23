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

TFormula ptError("ptError", "(4. / 2.) * \
    ( \
     (abs(y) < 2.66 && abs(y) >= 1.74 ) * 2 * sqrt( 8.56036e-05^2 * x^2 +0.0148987^2    ) + \
     (abs(y) < 1.74 && abs(y) >= 1.32 ) * sqrt( 8.56036e-05^2 * x^2 +0.0148987^2    ) + \
     (abs(y) < 1.32 && abs(y) >= 0.88 ) * sqrt( 1.12382e-05^2 * x^2 +0.00391722^2   ) + \
     (abs(y) < 0.88 && abs(y) >= 0.45 ) * sqrt( 1.16768e-05^2 * x^2 +0.00255204^2    ) + \
     (abs(y) < 0.45 && abs(y) >= 0.18 ) * sqrt( 1.28327e-05^2 * x^2 +0.00220587^2   ) + \
     (abs(y) < 0.18)                      * sqrt( 1.32845e-05^2 * x^2 +0.00209325^2   ) \
     )", 2); 

void performance(const char *inputFile = "delphes.root")
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
  Int_t pbins = 2000;
  Double_t pmin = 0., pmax = 100.; // log10(GeV/c)
  
  TH2 *hBetaP[2], *hDeltaT[2][5], *hNsigma[2][5], *hAll[5], *hAcceptance[2][5];
  TH2 *hCandidate[2][5], *hWrong[2][5];
  const char *tofname[2] = {"barrel", "forward"};
  const char *pname[5] = {"electron", "muon", "pion", "kaon", "proton"};
  const char *platex[5] = {"e", "#mu", "#pi", "K", "p"};

  for (Int_t ipart = 0; ipart < 5; ++ipart) {
    hAll[ipart] = new TH2F(Form("hAll_%s", pname[ipart]), ";#it{y};#it{p} (GeV/#it{c})", 500, -5., 5., pbins, pmin, pmax);
  }

  for (Int_t iTOF = 0; iTOF < 2; ++iTOF) {
    
    hBetaP[iTOF] = new TH2F(Form("hBetaP_%s", tofname[iTOF]), ";#it{p} (GeV/#it{c});v/#it{c}", pbins, pmin, pmax, 1000, 0.2, 1.2);
    
    for (Int_t ipart = 0; ipart < 5; ++ipart) {

      hAcceptance[iTOF][ipart] = new TH2F(Form("hAcceptance_%s_%s", tofname[iTOF], pname[ipart]), ";#it{y};#it{p} (GeV/#it{c})", 500, -5., 5., pbins, pmin, pmax);

      hDeltaT[iTOF][ipart] = new TH2F(Form("hDeltaT_%s_%s", tofname[iTOF], pname[ipart]), Form(";#it{p} (GeV/#it{c});#Deltat_{%s} (ps)", platex[ipart]), pbins, pmin, pmax, 200, -2.44, 2.44);
      hNsigma[iTOF][ipart] = new TH2F(Form("hNsigma_%s_%s", tofname[iTOF], pname[ipart]), Form(";#it{p} (GeV/#it{c});n#sigma_{%s}", platex[ipart]), pbins, pmin, pmax, 500, -25., 25.);

      
      hCandidate[iTOF][ipart] = new TH2F(Form("hCandidate_%s_%s", tofname[iTOF], pname[ipart]), ";#it{y};#it{p} (GeV/#it{c});", 500, -5., 5., pbins, pmin, pmax);
      hWrong[iTOF][ipart] = new TH2F(Form("hWrong_%s_%s", tofname[iTOF], pname[ipart]), ";#it{y};#it{p} (GeV/#it{c});", 500., -5., 5., pbins, pmin, pmax);
    }
      
  }
    
  // Loop over all events
  for (Int_t ientry = 0; ientry < numberOfEntries; ++ientry) {
    
    // Load selected branches with data from specified event
    treeReader->ReadEntry(ientry);
    
    // If event contains at least 1 track
    Long64_t numberOfTracks = branchTrack->GetEntries();
    if (numberOfTracks <= 0) continue;
    
    // loop over tracks
    for (Int_t itrack = 0; itrack < numberOfTracks; ++itrack) {
      Track *track = (Track *)branchTrack->At(itrack);
      GenParticle *particle = (GenParticle *) track->Particle.GetObject();

      Int_t m1 = particle->M1;
      GenParticle *mother = (GenParticle *)branchParticle->At(m1);
      //      if (TMath::Abs(mother->PID) != 4122) continue;

      Int_t pid[5] = {11, 13, 211, 321, 2212};

      Double_t y_true = particle->Rapidity;
      Double_t p_true = particle->P;
      Double_t pt_true = particle->PT;
      for (Int_t ipart = 0; ipart < 5; ++ipart)
	if (TMath::Abs(particle->PID) == pid[ipart])
	  hAll[ipart]->Fill(y_true, p_true);

      Double_t R = TMath::Sqrt(track->XOuter * track->XOuter +
			       track->YOuter * track->YOuter) * 0.1;
      Double_t Z = track->ZOuter * 0.1;

      /** define TOF hit (barrel/forward) **/
      Double_t bTOFr    = eTOF_radius; // [cm]
      Double_t bTOFz    = eTOF_length; // [cm]
      Double_t fTOFrin  = 50.;         // [cm]
      Double_t fTOFrout = eTOF_radius; // [cm]
      Double_t fTOFz    = eTOF_length; // [cm]
      Double_t TOFreso  = eTOF_sigmat;  // [ns]
      Int_t iTOF = -1;
      if (TMath::Abs(R - bTOFr) < 0.001 && TMath::Abs(Z) < bTOFr) iTOF = 0;
      else if (R < fTOFrout && R > fTOFrin && TMath::Abs(TMath::Abs(Z) - fTOFz) < 0.001) iTOF = 1;
      else continue;

      /** primary DCA cuts **/
      if (TMath::Abs(track->D0) > 0.1) continue;
      if (TMath::Abs(track->DZ) > 1.) continue;

      /** get info **/
      Double_t tof = track->TOuter * 1.e9; // [ns]
      Double_t L = track->L * 0.1; // [cm]
      Double_t eta = track->Eta;
      Double_t pt = track->PT; // [GeV/c]
      Double_t ept = pt * ptError.Eval(pt, eta); // [GeV/c]
      Double_t p = track->P; // [GeV/c]
      Double_t ep = ept * cosh(eta); // [GeV/c]
      Double_t c = 29.9792458; // [cm/ns]
      Double_t beta = L / tof / c;
      hBetaP[iTOF]->Fill(p, beta);

      /** compute expected TOF **/
      Double_t mass[5] = {0.00051099891, 0.10565800, 0.13957000, 0.49367700, 0.93827200};
      for (Int_t ipart = 0; ipart < 5; ++ipart) {

	if (TMath::Abs(track->PID) == pid[ipart]) {
	  hAcceptance[iTOF][ipart]->Fill(y_true, p_true);
	}

	Double_t texp = L / c / p * TMath::Sqrt(mass[ipart] * mass[ipart] + p * p);
	Double_t etexp = L / c * mass[ipart] * mass[ipart] / p / p / TMath::Sqrt(mass[ipart] * mass[ipart] + p * p) * ep;

	Double_t sigma = TMath::Sqrt(etexp * etexp + TOFreso * TOFreso);
	Double_t delta = tof - texp;
	Double_t nsigma = delta / sigma;
	hDeltaT[iTOF][ipart]->Fill(p, delta);
	hNsigma[iTOF][ipart]->Fill(p, nsigma);

	if (TMath::Abs(nsigma) < 3.) {
	  hCandidate[iTOF][ipart]->Fill(y_true, p_true);
	  if (TMath::Abs(track->PID) != pid[ipart])
	    hWrong[iTOF][ipart]->Fill(y_true, p_true);
	}
      }
      
    } // loop over tracks
  } // loop over events
  
  // Show resulting histograms
  TFile *fout = TFile::Open("performance.root", "RECREATE");
  for (Int_t ipart = 0; ipart < 5; ++ipart) {
    hAll[ipart]->Write();
  }
  for (Int_t iTOF = 0; iTOF < 2; ++iTOF) {
    hBetaP[iTOF]->Write();
    for (Int_t ipart = 0; ipart < 5; ++ipart) {
      hAcceptance[iTOF][ipart]->Write();
      hDeltaT[iTOF][ipart]->Write();
      hNsigma[iTOF][ipart]->Write();
      hCandidate[iTOF][ipart]->Write();
      hWrong[iTOF][ipart]->Write();
    }
  }
  fout->Close();
    
}

void post()
{

  // Show resulting histograms
  TFile *f = TFile::Open("analysis.root", "UPDATE");
  TH2 *hCandidate[2][5], *hWrong[2][5], *hContamination[2][5];
  const char *tofname[2] = {"barrel", "forward"};
  const char *pname[5] = {"electron", "muon", "pion", "kaon", "proton"};
  for (Int_t iTOF = 0; iTOF < 2; ++iTOF) {
    for (Int_t ipart = 0; ipart < 5; ++ipart) {
      hCandidate[iTOF][ipart] = (TH2 *)f->Get(Form("hCandidate_%s_%s", tofname[iTOF], pname[ipart]));
      hWrong[iTOF][ipart] = (TH2 *)f->Get(Form("hWrong_%s_%s", tofname[iTOF], pname[ipart]));
      hContamination[iTOF][ipart] = (TH2 *)hWrong[iTOF][ipart]->Clone(Form("hContamination_%s_%s", tofname[iTOF], pname[ipart]));
      hContamination[iTOF][ipart]->Sumw2();
      hContamination[iTOF][ipart]->Divide(hContamination[iTOF][ipart], hCandidate[iTOF][ipart], 1., 1., "B");
      hContamination[iTOF][ipart]->SetTitle(";#it{y};#it{p} (GeV/#it{c});contamination");
      hContamination[iTOF][ipart]->Write();
    }
  }
  f->Close();
  
}
