#!/bin/bash
#$ -S /bin/bash
#$ -o log/01_make_bin_for_tfbs/
#$ -e log/01_make_bin_for_tfbs/
#$ -cwd
#$ -l os7

set -eu
set -x

cd /home/kogu/analysis/ner_analysis

# --------------------------------------------------------------------------------------------
# Encode motifs under ChIP-seq motifs.

##########################
# define variables
##########################
tst0=2021-04-24 # timestamp or any distinguishable name

paths_dnase=/path/to/dnase

path_tfbs=/path/to/merged_ENCODE.tf.bound.union_sort.bed
tfbs_name="merged_ENCODE.tf.bound.union"

genome_size=/path/to/genome_size.tsv


##########################
tst=${tst0}_ENC
outdir=result_bin/$tst
slop_size="0"

##########################
# main funcition

for path_dnase in `eval echo ${paths_dnase}`; do
  dnase=`basename $path_dnase | perl -pe 's/\..*$//'`
  tfbs=$tfbs_name

  resource_d=${outdir}/resource/${dnase}/

  mkdir -p $resource_d


  # make slopped file (-b option: extend both side by x bp)
  for x in `echo $slop_size`; do
    # extend dnase window x bases in both side
    bedtools slop \
    -i ${path_dnase} \
    -g ${genome_size} \
    -b ${x} \
    > ${resource_d}/${dnase}.slop${x}.bed

    # TFBS windows related to DNAse window
    bedtools intersect \
    -wa \
    -a ${path_tfbs} \
    -b ${resource_d}/${dnase}.slop${x}.bed \
    -sorted \
    -u \
    > ${outdir}/${tfbs}__${dnase}.slop${x}.bed

    # TFBS windows NOT-related to DNAse window (-v option)
    bedtools intersect \
    -wa \
    -a ${path_tfbs} \
    -b ${resource_d}/${dnase}.slop${x}.bed \
    -sorted \
    -v \
    > ${outdir}/${tfbs}__${dnase}.slop${x}_v.bed

  done

done

#--------------------------------------------------------------------------------

# HBZ motifs


##########################
# define variables
##########################
tst=${tst0}_othermotifs

paths_dnase=/path/to/dnase
path_tfbs=/path/to/HBZ_motifs.tsv
genome_size=/path/to/genome_size.tsv


##########################
outdir=result_bin/$tst
slop_size="0"

##########################
# main funcition

for path_tfbs in `eval echo ${paths_tfbs}`
do
  for path_dnase in `eval echo ${paths_dnase}`
  do
    dnase=`basename $path_dnase .bed`
    tfbs=`basename $path_tfbs .bed`
    
    resource_d=${outdir}/resource/${dnase}/

    mkdir -p $resource_d
    # mkdir -p $stats_d
    
    # echo $dnase
    
    # make slopped file (-b option: extend both side by x bp)
    for x in `echo $slop_size`; do
      # extend dnase window x bases in both side
      bedtools slop \
      -i ${path_dnase} \
      -g ${genome_size} \
      -b ${x} \
      > ${resource_d}/${dnase}.slop${x}.bed
      
      # TFBS windows related to DNAse window
      bedtools intersect \
      -wa \
      -a ${path_tfbs} \
      -b ${resource_d}/${dnase}.slop${x}.bed \
      -sorted \
      -u \
      > ${outdir}/${tfbs}__${dnase}.slop${x}.bed
      
      # TFBS windows NOT-related to DNAse window (-v option)
      bedtools intersect \
      -wa \
      -a ${path_tfbs} \
      -b ${resource_d}/${dnase}.slop${x}.bed \
      -sorted \
      -v \
      > ${outdir}/${tfbs}__${dnase}.slop${x}_v.bed
      
    done
  done
done

