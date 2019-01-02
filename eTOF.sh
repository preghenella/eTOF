#! /bin/bash

nruns=1
workdir=""
analysis=""
params=""
### pythia variables
nevents=10000
sseed=123456789

### detector variables
bz=2.       # [T]
radius=100. # [cm]
length=150. # [cm]
sigmat=0.03 # [ns]
tracker=EIC

### EIC tracker definition
efficiency_EIC="\
                    (pt <= 0.1) * (0.00) + \
(abs(eta) <= 3.0) * (pt > 0.1)  * (1.00) + \
(abs(eta) >  3.0)               * (0.00)   \
"
ptresolution_EIC="\
(3.5 \/ $bz) * ( \
 (abs(eta) <= 3.0)                     * sqrt(0.001^2 + pt^2*1.e-5^2) + \
 (abs(eta) > 1.0 \&\& abs(eta) <= 3.0) * sqrt(0.01^2  + pt^2*1.e-4^2)   \
)"

### CLIC tracker definition
efficiency_CLIC="\
(abs(eta) > 2.54) * (0.000) + \
(energy >= 80) * (abs(eta) < 2.54)  * (1.000) + \
(energy < 80 \&\& energy >= 3) * (abs(eta) <=2.54 \&\& abs(eta) > 2.34)  * (0.994) + \
(energy < 80 \&\& energy >= 3) * (abs(eta) <= 2.34) * (1.000) + \
(energy < 3) * (abs(eta) <= 2.54 \&\& abs(eta) > 0.55 ) * (0.000) + \
(energy < 3) * (abs(eta) <= 0.55 ) * (1.000) \
"
ptresolution_CLIC="\
(4. \/ $bz) * ( \
 (abs(eta) < 2.66 \&\& abs(eta) >= 1.74 ) * 2 * sqrt(8.56036e-05^2 * pt^2 +0.0148987^2)  + \
 (abs(eta) < 1.74 \&\& abs(eta) >= 1.32 )     * sqrt(8.56036e-05^2 * pt^2 +0.0148987^2)  + \
 (abs(eta) < 1.32 \&\& abs(eta) >= 0.88 )     * sqrt(1.12382e-05^2 * pt^2 +0.00391722^2) + \
 (abs(eta) < 0.88 \&\& abs(eta) >= 0.45 )     * sqrt(1.16768e-05^2 * pt^2 +0.00255204^2) + \
 (abs(eta) < 0.45 \&\& abs(eta) >= 0.18 )     * sqrt(1.28327e-05^2 * pt^2 +0.00220587^2) + \
 (abs(eta) < 0.18)                            * sqrt(1.32845e-05^2 * pt^2 +0.00209325^2)   \
)"

### CMS tracker definition
efficiency_CMS="\
                                        (pt <= 0.1)               * (0.00) + \
(abs(eta) <= 1.5)                     * (pt > 0.1 \&\& pt <= 1.0) * (0.70) + \
(abs(eta) <= 1.5)                     * (pt > 1.0)                * (0.95) + \
(abs(eta) > 1.5 \&\& abs(eta) <= 2.5) * (pt > 0.1 \&\& pt <= 1.0) * (0.60) + \
(abs(eta) > 1.5 \&\& abs(eta) <= 2.5) * (pt > 1.0)                * (0.85) + \
(abs(eta) > 2.5)                                                  * (0.00)   \
"
ptresolution_CMS="\
(3.8 \/ $bz) * 0.1 * ( \
 (abs(eta) <= 0.5)                     * (pt > 0.1) * sqrt(0.06^2 + pt^2*1.3e-3^2) + \
 (abs(eta) > 0.5 \&\& abs(eta) <= 1.5) * (pt > 0.1) * sqrt(0.10^2 + pt^2*1.7e-3^2) + \
 (abs(eta) > 1.5 \&\& abs(eta) <= 2.5) * (pt > 0.1) * sqrt(0.25^2 + pt^2*3.1e-3^2)   \
)"

### ATLAS tracker definition
efficiency_ATLAS="\
                                        (pt <= 0.1)               * (0.00) + \
(abs(eta) <= 1.5)                     * (pt > 0.1 \&\& pt <= 1.0) * (0.70) + \
(abs(eta) <= 1.5)                     * (pt > 1.0)                * (0.95) + \
(abs(eta) > 1.5 \&\& abs(eta) <= 2.5) * (pt > 0.1 \&\& pt <= 1.0) * (0.60) + \
(abs(eta) > 1.5 \&\& abs(eta) <= 2.5) * (pt > 1.0)                * (0.85) + \
(abs(eta) > 2.5)                                                  * (0.00)   \
"
ptresolution_ATLAS="\
(2.0 \/ $bz) * 0.1 * ( \
 (abs(eta) <= 0.5)                     * (pt > 0.1) * sqrt(0.06^2 + pt^2*1.3e-3^2) + \
 (abs(eta) > 0.5 \&\& abs(eta) <= 1.5) * (pt > 0.1) * sqrt(0.10^2 + pt^2*1.7e-3^2) + \
 (abs(eta) > 1.5 \&\& abs(eta) <= 2.5) * (pt > 0.1) * sqrt(0.25^2 + pt^2*3.1e-3^2)   \
)"

### program options
while [ ! -z "$1" ]; do
    option="$1"
    shift
    if [ "$option" = "--help" ]; then
        echo "usage: ./eTOF.sh"
        echo "options:"
        echo "          --nruns     number of runs [$nruns]"
        echo "          --workdir   running directory"
        echo "          --nevents   number of events [$nevents]"
        echo "          --sseed     start seed [$sseed]"
        echo "          --params    Pythia6 parameters filename"
        echo "          --analysis  ROOT analysis macro filename"
        echo "          --bz        magnetic field [$bz] T"
        echo "          --radius    radius of TOF layer [$radius] cm"
        echo "          --length    length of TOF layer [$length] cm"
        echo "          --sigmat    time resolution of TOF layer [$sigmat] ns"
        echo "          --tracker   tracker definition [$tracker]"
        exit
    elif [ "$option" = "--running" ]; then
	nProcesses=$(ls .running.* 2> /dev/null | wc -l)
	ls .running.*
	echo "[---] currently running $nProcesses runs"
	exit
    elif [ "$option" = "--nruns" ]; then
        nruns="$1"
        shift
    elif [ "$option" = "--workdir" ]; then
        workdir="$1"
        shift
    elif [ "$option" = "--nevents" ]; then
        nevents="$1"
        shift
    elif [ "$option" = "--sseed" ]; then
        sseed="$1"
        shift
    elif [ "$option" = "--params" ]; then
        params="$1"
        shift
    elif [ "$option" = "--analysis" ]; then
        analysis="$1"
        shift
    elif [ "$option" = "--bz" ]; then
        bz="$1"
        shift
    elif [ "$option" = "--radius" ]; then
        radius=$1
	length=$radius
        shift
    elif [ "$option" = "--length" ]; then
        length=$1
        shift
    elif [ "$option" = "--sigmat" ]; then
        sigmat=$1
        shift
    elif [ "$option" = "--tracker" ]; then
        tracker=$1
        shift
    fi
done

main() {
    ### init
    splash
    check
    rm -rf $workdir
    ### run
    for runid in $(seq 0 $(($nruns-1)) ); do
	### pause if too many runs running
	nProcesses=$(ls .running.* 2> /dev/null | wc -l)
	if [ $nProcesses -ge 16 ]; then
	    echo "[---] sleeping: currently running $nProcesses runs"
	fi
	while [ $nProcesses -ge 16 ]; do
            sleep 5;
	    nProcesses=$(ls .running.* 2> /dev/null | wc -l)
	done
	run $runid &
	sleep 0.5
    done
    nProcesses=$(ls .running.* 2> /dev/null | wc -l)
    echo "[---] waiting: currently running $nProcesses runs"
    wait
    for ianalysis in $analysis; do
	analysisName=$(basename -- "$ianalysis")
	analysisName=${analysisName%.*}
	echo "[---] merging analysis: hadd -k -f $workdir/$analysisName.root $workdir/*/$analysisName.root"
	hadd -k -f $workdir/$analysisName.root $workdir/*/$analysisName.root > $workdir/hadd.$analysisName.log
    done
    grep "All included subprocesses" $workdir/*/pythia.log | awk '{ gsub("D", "E", $11); NEVENTS += $8; XSECTION += $11 } END { print NEVENTS " events"; print XSECTION/NR " (mb)" }' > $workdir/pythia.summary
    echo "[---] cleaning: find $workdir -type d -not -name $workdir -not -name 000 -exec rm -rf {} +"
    find $workdir -type d -not -name $workdir -not -name 000 -exec rm -rf {} +
    rm -rf $workdir/hadd.*log
    echo "[---] done"
}

check() {

    case $tracker in
	"ideal")
	    efficiency="1.0"
	    ptresolution="0.0"
	    ;;
	"EIC")
	    efficiency=$efficiency_EIC
	    ptresolution=$ptresolution_EIC
	    ;;
	"CMS")
	    efficiency=$efficiency_CMS
	    ptresolution=$ptresolution_CMS
	    ;;
	"ATLAS")
	    efficiency=$efficiency_ATLAS
	    ptresolution=$ptresolution_ATLAS
	    ;;
	"CLIC")
	    efficiency=$efficiency_CLIC
	    ptresolution=$ptresolution_CLIC
	    ;;
	"")
	    echo "[-E-] tracker definition not defined"
	    ./eTOF.sh --help
	    exit
	    ;;
	"*")
	    echo "[-E-] tracker definition '$tracker' does not exist"
	    exit
    esac

    if [ "$params" == "" ]; then
	echo "[-E-] Pythia6 parameters file not defined"
	./eTOF.sh --help
	exit
    elif [ ! -f $params ]; then
	echo "[-E-] Pythia6 parameters file '$params$' not found"
	exit
    fi

    if [ "$analysis" == "" ]; then
	echo "[-E-] ROOT analysis macro file not defined"
	./eTOF.sh --help
	exit
    else
	for ianalysis in $analysis; do
	    if [ ! -f $ianalysis ]; then
		echo "[-E-] ROOT analysis macro file '$ianalysis' not found"
		exit
	    fi
	done
    fi

    if [ "$workdir" == "" ]; then
	echo "error: working directory not defined"
	./eTOF.sh --help
	exit
    elif [ -d $workdir ]; then
	echo "[-W-] working directory '$workdir' already exists: ^C to interrupt"
	sleep 5
	echo "[-W-] working directory '$workdir' will be deleted"
    fi
}

splash() {

    echo "##### eTOF #####################"
    echo "bz           [T]  $bz"
    echo "radius      [cm]  $radius"
    echo "length      [cm]  $length"
    echo "sigmat      [ns]  $sigmat"
    echo "tracker           $tracker"
    echo "################################"
    echo "nruns             $nruns"
    echo "workdir           $workdir"
    echo "nevents           $nevents"
    echo "sseed             $sseed"
    echo "params            $params"
    echo "analysis          $analysis"
    echo "################################"

}

run() {

    seed=$(($sseed + $1 * 2))
    runid=$(printf "%03d" $1)
    touch .running.$workdir.$runid
    if [ $? -ne 0 ]; then
	echo "[$runid] error while initialising the run"
	rm -rf .running.$workdir.$runid
	exit
    fi
    rundir="$workdir/$runid"
    mkdir -p $rundir
    if [ $? -ne 0 ]; then
	echo "[$runid] error while creating run directory"
	rm -rf .running.$workdir.$runid
	exit
    fi
    mkfifo $rundir/pythia.hepmc
    if [ $? -ne 0 ]; then
	echo "[$runid] error while creating FIFO"
	rm -rf .running.$workdir.$runid
	exit
    fi
    generate_card $1 > $rundir/card.tcl
    if [ $? -ne 0 ]; then
	echo "[$runid] error while generating card"
	rm -rf .running.$workdir.$runid
	exit
    fi
    cp $params $rundir/pythia.params
    if [ $? -ne 0 ]; then
	echo "[$runid] error copying params"
	rm -rf .running.$workdir.$runid
	exit
    fi
    for ianalysis in $analysis; do
	analysisName=$(basename -- "$ianalysis")
	analysisName=${analysisName%.*}
	generate_analysis > $rundir/$analysisName.C
	if [ $? -ne 0 ]; then
	    echo "[$runid] error while generating analysis"
	    rm -rf .running.$workdir.$runid
	    exit
	fi
    done
    cd $rundir
    echo "[$runid] working directory: `pwd`" 
    ### run PYTHIA
    echo "[$runid] running PYTHIA: agile-runmc Pythia6:HEAD -s $seed -n $nevents -P pythia.params"
    agile-runmc Pythia6:HEAD -s $seed -n $nevents -P pythia.params -o pythia.hepmc >& pythia.log &
    ### run DELPHES
    echo "[$runid] running DELPHES: DelphesHepMC card.tcl delphes.root" 
    cat pythia.hepmc | DelphesHepMC card.tcl delphes.root - >& delphes.log
    ### run ROOT
    for ianalysis in $analysis; do
	analysisName=$(basename -- "$ianalysis")
	analysisName=${analysisName%.*}
	echo "[$runid] running ANALYSIS: root -l -b -q $analysisName.C"
	root -b -q $analysisName.C >& $analysisName.log
    done
    ### clean
    rm pythia.hepmc
    rm delphes.root
    cd - >& /dev/null
    rm .running.$workdir.$runid
    echo "[$runid] processing complete"

}

function join_by { local IFS="$1"; shift; echo "$*"; }

generate_analysis() {
    formula=$ptresolution
    formula=${formula//\\/}
    formula=${formula//pt/x}
    formula=${formula//eta/y}
    echo "#define eTOF_radius $radius // [cm]"
    echo "#define eTOF_length $length // [cm]"
    echo "#define eTOF_sigmat $sigmat // [ns]"
    echo "#define magneticBz  $bz     // [T]"
    echo
    echo "TFormula eTOF_ptresolution(\"eTOF_ptresolution\", \"$formula\", 2);"
    echo
    cat eTOF.C
    echo
    cat $analysis
}

generate_card() {

    echo "set RandomSeed $seed"
    echo
    cat cards/eTOF.tcl \
	| sed -e "s/eTOF_Radius/${radius}e-2/g" \
	| sed -e "s/eTOF_HalfLength/${length}e-2/g" \
	| sed -e "s/eTOF_TimeResolution/${sigmat}e-9/g" \
	| sed -e "s/eTOF_PResolutionFormula/${ptresolution}/g" \
	| sed -e "s/eTOF_EfficiencyFormula/${efficiency}/g" \
	| sed -e "s/eMAG_Bz/${bz}/g"
}

main