#!/bin/bash

# Mon 14 Feb 2022 05:13:14 PM PST
# script to convert a copy number .bed file (assumes chr/start/end/cn) into a comma separated fasta file (needed for cnp2cnp). Prints to stdout
# example fasta format is:
# > sampleName
# 2,1,2,3

p="$1" # "paramsA"
k="$2" # "10"
c="$3" # "20"
bedType="$4" # "mean", "mode", "median". Note cnp2cnp expects integers (which mean cannot guarantee)
fileKey="$5" # filekey to find output files
rev="$6" # T/F. should the order of the bed file be reversed?
nearest="$7" # "nearest10_"

dataset=$(echo "$p" | sed 's/\//_/')
bedFiles=$(find "$p" -maxdepth 1 -name "output*""$fileKey""_""$nearest""$dataset""_k""$k""_c""$c""*""$bedType""*.bed")

if [[ -z "$bedFiles" ]] ; then
  exit 1
fi

if [[ "$rev" == "T" ]] ; then
  bedFiles=$(echo "$bedFiles" | tac)
fi

echo "$bedFiles" | while read bed ; do
  line1=">""$bed"
  line2=$(awk '{printf("%.0f\n", $4)}' "$bed" | tr '\n' ','i | sed 's/,$//')
  
  echo "$line1"
  echo "$line2"
done

