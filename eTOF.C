#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

int eTOF(Track *track, double *deltat, double *nsigma)
{
  /** define TOF hit (barrel/forward) **/
  int etof = -1;
  double R = TMath::Sqrt(track->XOuter * track->XOuter +
			   track->YOuter * track->YOuter) * 0.1;
  double Z = track->ZOuter * 0.1;

  double Ro = fabs(R - eTOF_radius);
  double Za = fabs(fabs(Z) - eTOF_length);
  if (Ro == 0.) etof = 0;
  else if (Za == 0. && R > 50.) etof = 1;
  else return -1;

  /** get info **/
  double tof  = track->TOuter * 1.e9; // [ns]
  double L    = track->L * 0.1; // [cm]
  double eta  = track->Eta;
  double pt   = track->PT; // [GeV/c]
  double ept  = pt * eTOF_ptresolution.Eval(pt, eta); // [GeV/c]
  double p    = track->P; // [GeV/c]
  double p2   = p * p;
  double ep   = ept * cosh(eta); // [GeV/c]
  double c    = 29.9792458; // [cm/ns]
  double Lc   = L / c;

  /** perform PID **/
  double mass[5] = {0.00051099891, 0.10565800, 0.139570, 0.493600, 0.938270};
  for (Int_t ipart = 0; ipart < 5; ++ipart) {
    double mass2 = mass[ipart] * mass[ipart];
    double texp = Lc / p * TMath::Sqrt(mass2 + p2);
    double etexp = Lc * mass2 / p2 / TMath::Sqrt(mass2 + p2) * ep;    
    double sigma = TMath::Sqrt(etexp * etexp + eTOF_sigmat * eTOF_sigmat);
    deltat[ipart] = tof - texp;
    nsigma[ipart] = deltat[ipart] / sigma;
  }

  return etof;
}
