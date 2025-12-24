#!/bin/bash
# script to combine STAR gene feature counts - simplified version (not run or tested)

# to run
# cd /users/ramess/Code/nxf_workflows/ramess/OAC_Trial_RNASeq
# bash combine_featureCounts_v2.sh

# for ASPPs_in_mouse_pancreatic_cancer
#fcDir="/mnt/lustre/processed/Lu_lab/ASPPs_in_mouse_pancreatic_cancer/RNAseq/Results/featureCounts"
#fcDir="/mnt/lustre/users/ramess/ASPPs_in_mouse_pancreatic_cancer/Results_eYFP/featureCounts"
#fcPattern="STAR_gene.featureCounts.txt"
#resultFile="ASPPs_in_mouse_pancreatic_cancer_STAR_featureCounts_eYFP.tsv"
#exclHeaderPatterns=(".bam")

# for Fetal_Oesophagus_Development/RNASeq
#fcDir="/mnt/lustre/processed/Lu_lab/Fetal_Oesophagus_Development/RNASeq/Results/featureCounts"
#fcPattern="STAR_gene.featureCounts.txt"
#resultDir="${fcDir}"
#resultFile="${fcDir}/STAR_featureCounts.tsv"
#exclHeaderPatterns=()

# for OAC_Immuno_Trial/RNAseq
fcDir="/mnt/lustre/processed/Lu_lab/OAC_Immuno_Trial/RNAseq/Results/featureCounts"
fcPattern="STAR_gene.featureCounts.txt"
resultDir="${fcDir}"
resultFile="${fcDir}/STAR_featureCounts.tsv"
exclHeaderPatterns=(".bam")

# create array with list of files
IFS=' ' read -ra myarr <<< `ls ${fcDir}/*${fcPattern}`

# copy Geneid
cut -f 1 ${myarr[0]} | tail -n +2 > "${resultFile}"

# loop through files
for fc in "${myarr[@]}"
do
  # cut column with the number of reads assigned to the gene
  cut -f 7 ${fc} | tail -n +2 > fc.tsv
  
  mv "${resultFile}" tmp.tsv
  
  # combine read counts
  paste tmp.tsv fc.tsv > "${resultFile}"
  
  rm tmp.tsv
done

# clean up column headers
if [[ ! -z "${exclHeaderPatterns}" ]]
then
  for p in "${exclHeaderPatterns[@]}"
  do
  	echo "p: ${p}"
    sed -i "s@${p}@@g" "${resultFile}"
  done
fi
