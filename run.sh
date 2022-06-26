#!/bin/bash

SBATCH_OPTS="--mem=64GB --time=2-00:00:00 --partition=preempt"

NETWORK=data/networks/DREAM_files/dream_3.txt
#NETWORK=data/networks/test_net.txt
MIN_C=5
MAX_C=100
OUTPUT=data/results/DREAM-3
echo "MAX: $MAX_C , MIN: $MIN_C"
OFILE=${OUTPUT}/d3_${MIN_C}_${MAX_C}.json
if [ ! -d $OUTPUT ]; then mkdir -p $OUTPUT; fi
LOGFILE=data/results/DREAM-3/RESULTS.log


sbatch $SBATCH_OPTS -o $LOGFILE ./src/clustering_main.py --network $NETWORK --min-c $MIN_C --max-c $MAX_C --output $OFILE
