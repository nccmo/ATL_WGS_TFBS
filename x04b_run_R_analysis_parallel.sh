#!/bin/bash
#$ -S /bin/bash
#$ -o log/04b_analysis_parallel/launch/
#$ -e log/04b_analysis_parallel/launch/
#$ -cwd
#$ -l os7

set -eux

tst0=2021-04-24

num_samples=150

########################
tstENC=${tst0}_ENC
tstOTHER=${tst0}_othermotifs
########################


# qsub ENCODE result ------------------------------------------------

for slop in 0 ; do
  for window in 30 ;do
    for rep in 4 ; do
      for bin in 1000 ; do
        qsub -N slop${slop}_win${window}_rep${rep}b_bin${bin} ./R/04b_qsub_b.sh ${slop} ${window} ${rep} ${bin} ${tstENC} ${num_samples}
      done
    done
  done
done


# qsub other HBZ-motifs ---------------------------------------------

for slop in 0; do
  for window in 30; do
    for rep in 4; do
      for bin in 1000 ; do
        for cell in KK1; do
          for qval in 5; do
            for motif in 1 2; do
              qsub -N slop${slop}_win${window}_rep${rep}b_bin${bin}_${cell}_q${qval}_motif${motif} ./R/04b_qsub_hbz_chip_other_motifs.sh ${slop} ${window} ${rep} ${bin} ${cell} ${qval} ${motif} ${tstOTHER} ${num_samples}
            done
          done
        done
      done
    done
  done
done
