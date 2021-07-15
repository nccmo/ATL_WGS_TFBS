#!/bin/bash
#$ -S /bin/bash
#$ -o log/03_count_peripheral_mutation/
#$ -e log/03_count_peripheral_mutation/
#$ -cwd
#$ -l os7

set -eu
set -x

cd /home/kogu/analysis/ner_analysis

# ----------------------------------------------------------------------------------------------------
# I first named the "Flanking" region "Periphearal (peri*)" and used this term hereafter.
# This naming may be confusiing, but has not been corrected.
# ----------------------------------------------------------------------------------------------------


tst0=2021-04-24
path_mutation=/path/to/mutation.bed
num_samples=150
# variables usually not changed ---------------------------------------------------------
path_callable=/path/to/callable.bed.gz
path_cds=/path/to/PCAWG_test_genomic_elements_pc_only.bed6

window="30"


# -------------------------------------------------------------------------------------------------------
# Encode

##########################
# define variables
##########################
tst=${tst0}_ENC

##########################
# directories: same as file-02

bin_dir=result_bin/$tst
out_dir=result_count/$tst


##########################
# main funcition

for bed in $bin_dir/*.bed
do
  echo "-------------------"
  # echo $bed
  prefix=`basename $bed .bed`
  echo $prefix
  
  for x in $(echo $window)
  do
    resource_dir=${out_dir}/resource/wlen${x}
    mkdir -p ${resource_dir}
    # peripheral (2x+1)bp bin
    cat $bed \
    | awk -v "halfw=$x" 'BEGIN{OFS = "\t"; width = 2*halfw + 1 } NR>1{m = int(($2+$3)/2); print $1, m-1000-width, m-1000, $4; print $1, m+1000, m+1000+width, $4}' \
    | bedtools intersect \
    -a stdin \
    -b ${path_cds} \
    -v \
    | bedtools intersect \
    -a stdin \
    -b ${path_callable} \
    -wa \
   >  ${resource_dir}/${prefix}_peri__${x}bp_noCDS_callable__pre.bed

   # merge bin --------------------------------------------
    rm -rf ${resource_dir}/${prefix}_tmp_peri/
    mkdir -p ${resource_dir}/${prefix}_tmp_peri/

    cat ${resource_dir}/${prefix}_peri__${x}bp_noCDS_callable__pre.bed \
    | awk -v d="${resource_dir}/${prefix}_tmp_peri/" '{fh=d $4 ".bed"; print $0 >> fh}'
    for f in ${resource_dir}/${prefix}_tmp_peri/*.bed
    do 
      dir_n=`dirname $f`
      base_n=`basename $f .bed`
      sort -k 1,1 -k2,2n $f > $dir_n/$base_n.sort.bed
      bedtools merge -i $dir_n/$base_n.sort.bed | awk -v tf=$base_n '{print $0 "\t" tf}'  > $dir_n/$base_n.merged.bed
    done
    # summarize
    cat ${resource_dir}/${prefix}_tmp_peri/*.merged.bed | sort -k 1,1 -k2,2n > ${resource_dir}/${prefix}_peri__${x}bp_noCDS_callable.bed
   #  --------------------------------------------


    mkdir -p ${out_dir}/wlen${x}/
    
    bedtools intersect \
    -a ${resource_dir}/${prefix}_peri__${x}bp_noCDS_callable.bed \
    -b ${path_mutation} \
    -c \
    > ${out_dir}/wlen${x}/${prefix}_peri__${x}bp_noCDS_callable__count.txt
    
    echo "${prefix} flanking region (peri), window: ${x}bp"
    cat ${out_dir}/wlen${x}/${prefix}_peri__${x}bp_noCDS_callable__count.txt \
    | awk -v n_samples="$num_samples" '{len = $3-$2 + len; sum = sum + $5} END{print len, sum, len/NR, sum/len/n_samples*1e6}'
    
  done
done


# --------------------------------------------------------------------------
# HBZ motifs

##########################
# define variables
##########################
tst=${tst0}_othermotifs


##########################
# directories: same as file-02

bin_dir=result_bin/$tst
out_dir=result_count/$tst


##########################
# main funcition

for bed in $bin_dir/*.bed
do
  echo "-------------------"
  # echo $bed
  prefix=`basename $bed .bed`
  echo $prefix
  
  for x in $(echo $window)
  do
    resource_dir=${out_dir}/resource/wlen${x}
    mkdir -p ${resource_dir}
    # peripheral (2x+1)bp bin
    cat $bed \
    | awk -v "halfw=$x" 'BEGIN{OFS = "\t"; width = 2*halfw + 1 } NR>1{m = int(($2+$3)/2); print $1, m-1000-width, m-1000, "HBZ"; print $1, m+1000, m+1000+width, "HBZ"}' \
    | bedtools intersect \
    -a stdin \
    -b ${path_cds} \
    -v \
    | bedtools intersect \
    -a stdin \
    -b ${path_callable} \
    -wa \
   >  ${resource_dir}/${prefix}_peri__${x}bp_noCDS_callable__pre.bed


   # merge bin --------------------------------------------
    rm -rf ${resource_dir}/${prefix}_tmp_peri/
    mkdir -p ${resource_dir}/${prefix}_tmp_peri/

    cp ${resource_dir}/${prefix}_peri__${x}bp_noCDS_callable__pre.bed  ${resource_dir}/${prefix}_tmp_peri/

    for f in ${resource_dir}/${prefix}_tmp_peri/*.bed
    do 
      dir_n=`dirname $f`
      base_n=`basename $f .bed`
      sort -k 1,1 -k2,2n $f > $dir_n/$base_n.sort.bed
      bedtools merge -i $dir_n/$base_n.sort.bed | awk -v tf="HBZ" '{print $0 "\t" tf}'  > $dir_n/$base_n.merged.bed
    done
    # summarize
    cat ${resource_dir}/${prefix}_tmp_peri/*.merged.bed | sort -k 1,1 -k2,2n > ${resource_dir}/${prefix}_peri__${x}bp_noCDS_callable.bed
 
   #  --------------------------------------------
    
    mkdir -p ${out_dir}/wlen${x}/
    
    bedtools intersect \
    -a ${resource_dir}/${prefix}_peri__${x}bp_noCDS_callable.bed \
    -b ${path_mutation} \
    -c \
    > ${out_dir}/wlen${x}/${prefix}_peri__${x}bp_noCDS_callable__count.txt
    
    echo "${prefix} flanking region (peri), window: ${x}bp"
    cat ${out_dir}/wlen${x}/${prefix}_peri__${x}bp_noCDS_callable__count.txt \
    | awk '{len = $3-$2 + len; sum = sum + $5} END{if(NR>0) print len, sum, len/NR, sum/len/72*1e6}'
    
  done
done

