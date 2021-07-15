#!/bin/bash
#$ -S /bin/bash
#$ -o log/02_count_overlapping_mutation/
#$ -e log/02_count_overlapping_mutation/
#$ -cwd
#$ -l os7

set -eu
# set -x
set -o pipefail

cd /path/to/analysis_dir

tst0=2021-04-24
path_mutation=/path/to/mutation.bed
num_samples=150


# variables usually not changed ---------------------------------------------------------
path_callable=/path/to/callable.bed.gz
path_cds=/path/to/PCAWG_test_genomic_elements_pc_only.bed6



window="30"


# -------------------------------------------------------------------------------------------------------
echo ""
echo Encode
echo ""

##########################
# define variables
##########################
tst=${tst0}_ENC


##########################
bin_dir=result_bin/$tst
out_dir=result_count/$tst


##########################
# ENCODE 1/3 step

for bed in $bin_dir/*.bed ; do
  echo "-------------------"

  # echo $bed
  prefix=`basename $bed .bed`
  echo $prefix

  for x in $(echo $window)
  do
    resource_dir=${out_dir}/resource/wlen${x}
    mkdir -p ${resource_dir}

    cat $bed \
    | awk -v "window=$x" 'BEGIN{OFS = "\t"} NR>1{m = int(($2+$3)/2); print $1, m-1-window, m+window, $4}' \
    | bedtools intersect \
        -a stdin \
        -b ${path_cds} \
        -v \
    | bedtools intersect \
        -a stdin \
        -b ${path_callable} \
        -wa \
        >  ${resource_dir}/${prefix}__${x}bp_noCDS_callable__pre.bed
  #  >  ${resource_dir}/${prefix}__${x}bp_noCDS_callable.bed


   # merge bin --------------------------------------------
    rm -rf ${resource_dir}/${prefix}_tmp/
    mkdir -p ${resource_dir}/${prefix}_tmp/

    cat ${resource_dir}/${prefix}__${x}bp_noCDS_callable__pre.bed \
    | awk -v d="${resource_dir}/${prefix}_tmp/" '{fh=d $4 ".bed"; print $0 >> fh}'
    for f in ${resource_dir}/${prefix}_tmp/*.bed
    do 
      dir_n=`dirname $f`
      base_n=`basename $f .bed`
      sort -k 1,1 -k2,2n $f > $dir_n/$base_n.sort.bed
      bedtools merge -i $dir_n/$base_n.sort.bed | awk -v tf=$base_n '{print $0 "\t" tf}'  > $dir_n/$base_n.merged.bed
    done
    # summarize
    cat ${resource_dir}/${prefix}_tmp/*.merged.bed | sort -k 1,1 -k2,2n > ${resource_dir}/${prefix}__${x}bp_noCDS_callable.bed
   #  --------------------------------------------

    mkdir -p ${out_dir}/wlen${x}/

    bedtools intersect \
        -a ${resource_dir}/${prefix}__${x}bp_noCDS_callable.bed \
        -b ${path_mutation} \
        -c \
        > ${out_dir}/wlen${x}/${prefix}__${x}bp_noCDS_callable__count.txt

    echo "${prefix}, window: ${x}bp"
    cat ${out_dir}/wlen${x}/${prefix}__${x}bp_noCDS_callable__count.txt \
    | awk -v n_samples="$num_samples" '{len = $3-$2 + len; sum = sum + $5} END{print len, sum, len/NR, sum/len/n_samples*1e6}'

  done
done


# --------------------------------------------------------------------------
echo ""
echo HBZ motifs
echo ""
##########################
# define variables
##########################
tst=${tst0}_othermotifs

#########################
bin_dir=result_bin/$tst
out_dir=result_count/$tst


##########################
# main funcition 3/3 step

for bed in $bin_dir/*.bed ; do
  echo "-------------------"

  # echo $bed
  prefix=`basename $bed .bed`
  # echo $prefix

  for x in $(echo $window)
    do
    resource_dir=${out_dir}/resource/wlen${x}
    mkdir -p ${resource_dir}

    cat $bed \
    | awk -v "window=$x" 'BEGIN{OFS = "\t"} NR>1{m = int(($2+$3)/2); print $1, m-1-window, m+window, "HBZ"}' \
    | bedtools intersect \
        -a stdin \
        -b ${path_cds} \
        -v \
    | bedtools intersect \
        -a stdin \
        -b ${path_callable} \
        -wa \
          >  ${resource_dir}/${prefix}__${x}bp_noCDS_callable__pre.bed
  #  >  ${resource_dir}/${prefix}__${x}bp_noCDS_callable.bed


   # merge bin --------------------------------------------
    rm -rf ${resource_dir}/${prefix}_tmp/
    mkdir -p ${resource_dir}/${prefix}_tmp/

    cp ${resource_dir}/${prefix}__${x}bp_noCDS_callable__pre.bed ${resource_dir}/${prefix}_tmp/

    for f in ${resource_dir}/${prefix}_tmp/*.bed
    do 
      dir_n=`dirname $f`
      base_n=`basename $f .bed`
      sort -k 1,1 -k2,2n $f > $dir_n/$base_n.sort.bed
      bedtools merge -i $dir_n/$base_n.sort.bed | awk -v tf="HBZ" '{print $0 "\t" tf}'  > $dir_n/$base_n.merged.bed
    done
    
    # summarize
    cat ${resource_dir}/${prefix}_tmp/*.merged.bed | sort -k 1,1 -k2,2n > ${resource_dir}/${prefix}__${x}bp_noCDS_callable.bed
   #  --------------------------------------------

    mkdir -p ${out_dir}/wlen${x}/

    bedtools intersect \
        -a ${resource_dir}/${prefix}__${x}bp_noCDS_callable.bed \
        -b ${path_mutation} \
        -c \
        > ${out_dir}/wlen${x}/${prefix}__${x}bp_noCDS_callable__count.txt

    echo "${prefix}, window: ${x}bp"
    cat ${out_dir}/wlen${x}/${prefix}__${x}bp_noCDS_callable__count.txt \
    | awk -v n_samples="$num_samples" '{len = $3-$2 + len; sum = sum + $5} END{if(NR>0)print len, sum, len/NR, sum/len/n_samples*1e6}'

  done
done


exit 0
