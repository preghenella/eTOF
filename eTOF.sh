#! /bin/bash

nruns=1
workdir=""
analysis=""
params=""
### pythia variables
nevents=10000

### detector variables
bz=2.          # [T]
radius=200.    # [cm]
length=$radius # [cm]
sigmat=0.03    # [ns]

### program options
while [ ! -z "$1" ]; do
    option="$1"
    shift
    if [ "$option" = "--help" ]; then
        echo "usage: ./eTOF.sh"
        echo "options:"
        echo "          --nruns     number of runs [$nruns]"
        echo "          --workdir   running directory [$workdir]"
        echo "          --nevents   number of events [$nevents]"
        echo "          --seed      random seed [$seed]"
        echo "          --params    Pythia6 parameters filename"
        echo "          --analysis  ROOT analysis macro filename"
        echo "          --bz        magnetic field in T [$bz]"
        echo "          --radius    position of TOF layer $radius"
        echo "          --length    length of TOF layer $length"
        echo "          --sigmat    time resolution of TOF layer $length"
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
    elif [ "$option" = "--seed" ]; then
        seed="$1"
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
    fi
done

main() {
    ### check
    check
    ### init
    splash
    ### run
    for runid in $(seq 0 $(($nruns-1)) ); do
	### pause if too many runs running
	nProcesses=$(ls .running.* 2> /dev/null | wc -l)
#	echo "[...] currently running: $nAgileProcesses agile-runmc processes"
#	echo "[...] currently running: $nROOTProcesses root.exe processes"
	echo "[...] currently running: $nProcesses runs"
	while [ $nProcesses -ge 10 ]; do
            sleep 5s;
	    nProcesses=$(ls .running.* 2> /dev/null | wc -l)
	    echo "[...] currently running: $nProcesses runs"
	done
	run $runid &
	sleep 1s
    done
    wait
}

check() {

    if [ "$params" == "" ]; then
	echo "error: Pythia6 parameters filename not defined"
	./eTOF.sh --help
	exit
    elif [ ! -f $params ]; then
	echo "error: Pythia6 parameters filename not found ($params)"
	exit
    fi

    if [ "$analysis" == "" ]; then
	echo "error: ROOT analysis macro filename not defined"
	./eTOF.sh --help
	exit
    else
	for ianalysis in $analysis; do
	    if [ ! -f $ianalysis ]; then
		echo "error: ROOT analysis macro filename not found ($ianalysis)"
		exit
	    fi
	done
    fi
}

splash() {

    echo "##### eTOF #####################"
    echo "bz           [T]  $bz"
    echo "radius      [cm]  $radius"
    echo "length      [cm]  $length"
    echo "sigmat      [ns]  $sigmat"
    echo "################################"
    echo "nruns             $nruns"
    echo "workdir           $workdir"
    echo "nevents           $nevents"
    echo "params            $params"
    echo "analysis          $analysis"
    echo "################################"

}

run() {

    seed=$((123456789 + $1 * 2))
    runid=$(printf "%03d" $1)
    touch .running.$runid
    rundir="$workdir/$runid"
    mkdir -p $rundir
    rm -rf $rundir/*
    mkfifo $rundir/pythia.hepmc
    generate_card $1 > $rundir/card.tcl
    cp $params $rundir/pythia.params
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
	generate_analysis > $rundir/$analysisName.C
	echo "[$runid] running ANALYSIS: root -l -b -q $analysisName.C"
	root -b -q $analysisName.C >& $analysisName.log
    done
    ### clean
    rm pythia.hepmc
    cd - >& /dev/null
    rm .running.$runid
    echo "[$runid] processing complete"

}

function join_by { local IFS="$1"; shift; echo "$*"; }

generate_analysis() {
    echo "#define eTOF_radius $radius // [cm]"
    echo "#define eTOF_length $length // [cm]"
    echo "#define eTOF_sigmat $sigmat // [ns]"
    echo "#define magneticBz  $bz     // [T]"
    echo
    cat $analysis
}

generate_card() {

    seed=$((123456789 + $1 * 2))
    echo "#######################################"
    echo "# External variables"
    echo "#######################################"
    echo
    echo "set eTOF_Radius "$radius"e-2"
    echo "set eTOF_HalfLength "$length"e-2"
    echo "set eTOF_TimeResolution "$sigmat"e-9"
    echo "set eMAG_Bz $bz"
    echo "set RandomSeed $seed"
    echo
    cat eTOF_base.tcl
}

main