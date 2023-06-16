#!/bin/bash

EMAIL=faith.ocitti@tufts.edu
MIN_C=3
MAX_C=100
CLUSTER_OUTPUT=../results/simple_inverse/clusters
LOGFILE_OUTPUT=../results/simple_inverse/logs
echo "MAX: $MAX_C , MIN: $MIN_C"
if [ ! -d $CLUSTER_OUTPUT ]; then mkdir -p $CLUSTER_OUTPUT; fi
if [ ! -d $LOGFILE_OUTPUT ]; then mkdir -p $LOGFILE_OUTPUT; fi

for DREAM in 1 2 3
do
    for IC in 3 10 50 100
    do
        OFILE=../results/simple_inverse/clusters/d${DREAM}_ic${IC}
        NETWORK=../data/networks/DREAM_files/dream_${DREAM}.txt
        LOGFILE=../results/simple_inverse/logs/d${DREAM}_ic${IC}
        SBATCH_OPTS="--job-name=d${DREAM}-ic${IC}-si --nodes=1 --ntasks=2 --cpus-per-task=5 --mem=32GB --time=07-00:00:00 --output=${LOGFILE}.%j.out --error=${LOGFILE}.%j.err --mail-type=ALL --mail-user=${EMAIL} --partition=preempt "
        echo "DREAM: $DREAM, IC: $IC, OUTPUT_FILE: $OFILE"
        echo $SBATCH_OPTS
        sbatch $SBATCH_OPTS ../src/clustering_main.py --network $NETWORK --min-c $MIN_C --max-c $MAX_C --output $OFILE --spectral $IC --simple-inverse 1
    done
done
