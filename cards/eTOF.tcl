#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {

    ParticlePropagator

    TrackMerger
    TrackSmearing
    TimeSmearing

    TreeWriter

}

#################################
# Propagate particles in cylinder
#################################

module ParticlePropagator ParticlePropagator {

    set InputArray Delphes/stableParticles
    set OutputArray stableParticles
    set ChargedHadronOutputArray chargedHadrons
    set ElectronOutputArray electrons
    set MuonOutputArray muons

    # radius of the magnetic field coverage, in m
    set Radius $eTOF_Radius
    # half-length of the magnetic field coverage, in m
    set HalfLength $eTOF_HalfLength

    # magnetic field
    set Bz $eMAG_Bz
    
}

##############
# Track merger
##############

module Merger TrackMerger {
    # add InputArray InputArray
    add InputArray ParticlePropagator/chargedHadrons
    add InputArray ParticlePropagator/electrons
    add InputArray ParticlePropagator/muons
    set OutputArray tracks
}

########################################
# Track smearing
########################################

module TrackSmearing TrackSmearing {
    set InputArray TrackMerger/tracks
    set OutputArray tracks
    # rescaled CLIC resolution formula 
    set PResolutionFormula {
	(4. / 2.) * 
	(
	 (abs(eta) < 2.66 && abs(eta) >= 1.74 ) * 2 * sqrt( 8.56036e-05^2 * pt^2 +0.0148987^2    ) +
	 (abs(eta) < 1.74 && abs(eta) >= 1.32 ) * sqrt( 8.56036e-05^2 * pt^2 +0.0148987^2    ) +
	 (abs(eta) < 1.32 && abs(eta) >= 0.88 ) * sqrt( 1.12382e-05^2 * pt^2 +0.00391722^2   ) +
	 (abs(eta) < 0.88 && abs(eta) >= 0.45 ) * sqrt( 1.16768e-05^2 * pt^2 +0.00255204^2    ) +
	 (abs(eta) < 0.45 && abs(eta) >= 0.18 ) * sqrt( 1.28327e-05^2 * pt^2 +0.00220587^2   ) +
	 (abs(eta) < 0.18)                      * sqrt( 1.32845e-05^2 * pt^2 +0.00209325^2   )
	 )
    }
    ###
    set CtgThetaResolutionFormula { 0.0 }
    set PhiResolutionFormula { 0.0 }
    set D0ResolutionFormula { 0.0 }
    set DZResolutionFormula { 0.0 }
}

####################################
# Time smearing
####################################

module TimeSmearing TimeSmearing {
    set InputArray TrackSmearing/tracks
    set OutputArray tracks
    set TimeResolution $eTOF_TimeResolution
}

##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
    # add Branch InputArray BranchName BranchClass
    add Branch Delphes/allParticles Particle GenParticle
    add Branch TimeSmearing/tracks Track Track
}
