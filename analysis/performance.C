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
  Int_t pbins = 200;
  double pmin = -2., pmax = 2.; // log10(GeV/c)
  
  TH1 *hNhits[2];
  TH2 *hBetaP[2], *hDeltaT[2][5], *hNsigma[2][5], *hAll[5], *hAcceptance[2][5];
  TH2 *hCandidate[2][5], *hWrong[2][5];
  const char *tofname[2] = {"barrel", "forward"};
  const char *pname[5] = {"electron", "muon", "pion", "kaon", "proton"};
  const char *platex[5] = {"e", "#mu", "#pi", "K", "p"};

  for (Int_t ipart = 0; ipart < 5; ++ipart) {
    hAll[ipart] = new TH2F(Form("hAll_%s", pname[ipart]), ";#it{y};#it{p} (GeV/#it{c})", 500, -5., 5., pbins, pmin, pmax);
  }

  for (Int_t iTOF = 0; iTOF < 2; ++iTOF) {
    
    hNhits[iTOF] = new TH1F(Form("hNhits_%s", tofname[iTOF]), "", 100, 0., 100.);

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
    int nhits[2] = {0, 0};
    for (Int_t itrack = 0; itrack < numberOfTracks; ++itrack) {
      Track *track = (Track *)branchTrack->At(itrack);
      GenParticle *particle = (GenParticle *) track->Particle.GetObject();

      Int_t m1 = particle->M1;
      GenParticle *mother = (GenParticle *)branchParticle->At(m1);
      //      if (TMath::Abs(mother->PID) != 4122) continue;

      Int_t pid[5] = {11, 13, 211, 321, 2212};

      double y_true = particle->Rapidity;
      double p_true = particle->P;
      double pt_true = particle->PT;
      for (Int_t ipart = 0; ipart < 5; ++ipart)
	if (TMath::Abs(particle->PID) == pid[ipart])
	  hAll[ipart]->Fill(y_true, log10(p_true));

      /** primary DCA cuts **/
      if (TMath::Abs(track->D0) > 0.1) continue;
      if (TMath::Abs(track->DZ) > 1.) continue;

      /** eTOF cuts **/
      double deltat[5], nsigma[5];
      int etof = eTOF(track, deltat, nsigma);
      if (etof < 0) continue;
      nhits[etof]++;

      /** get info **/
      double tof = track->TOuter * 1.e9; // [ns]
      double L   = track->L * 0.1; // [cm]
      double p   = track->P; // [GeV/c]
      double c   = 29.9792458; // [cm/ns]
      double beta = L / tof / c;
      hBetaP[etof]->Fill(log10(p), beta);

      for (Int_t ipart = 0; ipart < 5; ++ipart) {
	
	if (TMath::Abs(track->PID) == pid[ipart])
	  hAcceptance[etof][ipart]->Fill(y_true, log10(p_true));
	
	hDeltaT[etof][ipart]->Fill(log10(p), deltat[ipart]);
	hNsigma[etof][ipart]->Fill(log10(p), nsigma[ipart]);
	
	if (TMath::Abs(nsigma[ipart]) < 3.) {
	  hCandidate[etof][ipart]->Fill(y_true, log10(p_true));
	  if (TMath::Abs(track->PID) != pid[ipart])
	    hWrong[etof][ipart]->Fill(y_true, log10(p_true));
	}
      }
      
    } // loop over tracks

    hNhits[0]->Fill(nhits[0]);
    hNhits[1]->Fill(nhits[1]);

  } // loop over events
  
  // Show resulting histograms
  TFile *fout = TFile::Open("performance.root", "RECREATE");
  for (Int_t ipart = 0; ipart < 5; ++ipart) {
    hAll[ipart]->Write();
  }
  for (Int_t etof = 0; etof < 2; ++etof) {
    hBetaP[etof]->Write();
    hNhits[etof]->Write();
    for (Int_t ipart = 0; ipart < 5; ++ipart) {
      hAcceptance[etof][ipart]->Write();
      hDeltaT[etof][ipart]->Write();
      hNsigma[etof][ipart]->Write();
      hCandidate[etof][ipart]->Write();
      hWrong[etof][ipart]->Write();
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
  for (Int_t etof = 0; etof < 2; ++etof) {
    for (Int_t ipart = 0; ipart < 5; ++ipart) {
      hCandidate[etof][ipart] = (TH2 *)f->Get(Form("hCandidate_%s_%s", tofname[etof], pname[ipart]));
      hWrong[etof][ipart] = (TH2 *)f->Get(Form("hWrong_%s_%s", tofname[etof], pname[ipart]));
      hContamination[etof][ipart] = (TH2 *)hWrong[etof][ipart]->Clone(Form("hContamination_%s_%s", tofname[etof], pname[ipart]));
      hContamination[etof][ipart]->Sumw2();
      hContamination[etof][ipart]->Divide(hContamination[etof][ipart], hCandidate[etof][ipart], 1., 1., "B");
      hContamination[etof][ipart]->SetTitle(";#it{y};#it{p} (GeV/#it{c});contamination");
      hContamination[etof][ipart]->Write();
    }
  }
  f->Close();
  
}
