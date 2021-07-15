#!/bin/bash
#$ -S /bin/bash
#$ -o log/04b_analysis_parallel/
#$ -e log/04b_analysis_parallel/
#$ -cwd
#$ -l os7


# ${slop} ${window} ${rep} ${bin} ${tstENC} ${n_sample}

# $1 $2 $3 $4 will be recieved as below:
# -------------------------------------
# slop <- commandArgs(trailingOnly = TRUE)[1]
# window <- commandArgs(trailingOnly = TRUE)[2]
# rep <- commandArgs(trailingOnly = TRUE)[3]
# bin <- commandArgs(trailingOnly = TRUE)[4]
# -------------------------------------

Rscript --vanilla --slave /home/kogu/analysis/ner_analysis/R/04b_analysis_parallel_b.R $1 $2 $3 $4 $5 $6
                          
