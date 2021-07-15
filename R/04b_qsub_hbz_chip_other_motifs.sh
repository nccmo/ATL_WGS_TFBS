#!/bin/bash
#$ -S /bin/bash
#$ -o log/04b_analysis_parallel/
#$ -e log/04b_analysis_parallel/
#$ -cwd
#$ -l os7

# -------------------------------------
# slop <- commandArgs(trailingOnly = TRUE)[1] %>% as.numeric()
# window <- commandArgs(trailingOnly = TRUE)[2] %>% as.numeric()
# rep <- commandArgs(trailingOnly = TRUE)[3] %>% as.numeric()
# bin <- commandArgs(trailingOnly = TRUE)[4] %>% as.numeric()
# cell <- commandArgs(trailingOnly = TRUE)[5] %>% as.character()
# qval <- commandArgs(trailingOnly = TRUE)[6] %>% as.numeric()
# -------------------------------------



Rscript --vanilla --slave /home/kogu/analysis/ner_analysis/R/04b_analysis_parallel_hbz_chip_other_motifs.R $1 $2 $3 $4 $5 $6 $7 $8 $9

