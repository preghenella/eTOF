fchunk=0
lchunk=0

nruns=100
nevents=100000
tracker=ATLAS
#analysis="analysis/secondary_sigma.C"
analysis="analysis/performance.C"
process="all"
#q2min=(0 1 2 4  8 16 32  64) 
#q2max=(1 2 4 8 16 32 64 128)
q2min=(1)
q2max=(-1)

for iprocess in $process; do
for itracker in $tracker; do
for iq2 in $(seq 0 $((${#q2min[@]} - 1))); do
for ichunk in $(seq $fchunk $lchunk); do

    sseed=$((123456789 + $ichunk * $nruns * 10))
    workdir="${iprocess}_processes__${q2min[$iq2]}_q2_${q2max[$iq2]}__chunk_${ichunk}"
    params/mkpypar.sh --process $iprocess --Q2min ${q2min[$iq2]} --Q2max ${q2max[$iq2]} --output $workdir.pythia.params
    echo "[$workdir] started: logging to '$workdir.log'"
    time -p ./eTOF.sh --workdir $workdir  --tracker $itracker --params $workdir.pythia.params  --analysis "$analysis" --nruns $nruns --nevents $nevents --sseed $sseed &> $workdir.log
    rm -rf $workdir.pythia.params
    echo "[$workdir] done"

done
done
done
done