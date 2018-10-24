#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {

    ParticlePropagator

    TrackMerger
    TrackingEfficiency
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

    set Radius eTOF_Radius
    set HalfLength eTOF_HalfLength
    set Bz eMAG_Bz    
}

##############
# Track merger
##############

module Merger TrackMerger {
    add InputArray ParticlePropagator/chargedHadrons
    add InputArray ParticlePropagator/electrons
    add InputArray ParticlePropagator/muons
    set OutputArray tracks
}

####################################
# Tracking efficiency
####################################

module Efficiency TrackingEfficiency {
    set InputArray TrackMerger/tracks
    set OutputArray tracks
    set EfficiencyFormula { eTOF_EfficiencyFormula }       
}

########################################
# Track smearing
########################################

module TrackSmearing TrackSmearing {
    set InputArray TrackingEfficiency/tracks
    set OutputArray tracks
    set PResolutionFormula { eTOF_PResolutionFormula }
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
    set TimeResolution eTOF_TimeResolution
}

##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
    # add Branch InputArray BranchName BranchClass
    add Branch Delphes/allParticles Particle GenParticle
    add Branch TimeSmearing/tracks Track Track
}
