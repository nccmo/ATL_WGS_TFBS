#!/bin/bash
#$ -S /bin/bash
#$ -o log/runall/
#$ -e log/log/runall/
#$ -cwd
#$ -l os7

DATETIME=`date +%Y%m%d_%H%M%S%3N` # YMD+HMS+millisec

qsub -N tfbs01_$DATETIME x01_make_bin_for_tfbs_combined.sh
qsub -N tfbs02_$DATETIME -hold_jid tfbs01_$DATETIME x02_count_overlapping_mutation_combined.sh
qsub -N tfbs03_$DATETIME -hold_jid tfbs02_$DATETIME x03_count_peripheral_mutation.sh
qsub -N tfbs04_$DATETIME -hold_jid tfbs03_$DATETIME x04b_run_R_analysis_parallel.sh

