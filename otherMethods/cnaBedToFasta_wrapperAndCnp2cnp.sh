#!/bin/bash

# Mon 14 Feb 2022 09:16:12 PM PST
# wrapper to call cnaBedToFasta.sh and cnp2cnp on all parameter sets
# assumes lists of simulated tumor files are available

# cnp2cnp: git clone git@github.com:AEVO-lab/cnp2cnp.git
# Note: zzs is the adaptation of MEDICC implemented in cnp2cnp

cnp2cnpDir="/space/s2/sandra/methodComparisons/cnp2cnp"
k=10
for p in "paramsA" "paramsB" "paramsC" "paramsD" ; do
  dataset=$(echo "$p" | sed 's/\//_/')
  for c in "20" "40" "60" "80" "100" "120" ; do
    # first do true bed files if they don't already exist
    trueBase="$p""/tumor_depths_""$c""_true"
    if [[ ! -f "$trueBase""Rev.zzs.cnp2cnp" ]] ; then
      echo -n > "$trueBase"".fasta"
      while read simu ; do
        bed=$(echo $simu | sed -e 's/simu/true/' -e 's/depth/bed/')
        line1=">""$bed"
        line2=$(awk '{printf("%.0f\n", $4)}' "$bed" | tr '\n' ','i | sed 's/,$//')
        echo -e "$line1""\n""$line2" >> "$trueBase"".fasta"
      done < "$p""/tumor_depths_""$c"
      python "$cnp2cnpDir""/cnp2cnp.py" -m "matrix" -d any -i "$trueBase"".fasta" -o "$trueBase"".cnp2cnp"
      python "$cnp2cnpDir""/cnp2cnp.py" -m "matrix" -d zzs -i "$trueBase"".fasta" -o "$trueBase"".zzs.cnp2cnp"

      echo -n > "$trueBase""Rev.fasta"
      tac "$p""/tumor_depths_""$c" | while read simu ; do
        bed=$(echo $simu | sed -e 's/simu/true/' -e 's/depth/bed/')
        line1=">""$bed"
        line2=$(awk '{printf("%.0f\n", $4)}' "$bed" | tr '\n' ','i | sed 's/,$//')
        echo -e "$line1""\n""$line2" >> "$trueBase""Rev.fasta"
      done 
      python "$cnp2cnpDir""/cnp2cnp.py" -m "matrix" -d any -i "$trueBase""Rev.fasta" -o "$trueBase""Rev.cnp2cnp"
      python "$cnp2cnpDir""/cnp2cnp.py" -m "matrix" -d zzs -i "$trueBase"".fasta" -o "$trueBase""Rev.zzs.cnp2cnp"
    fi

    for fileKey in "sconce2" "sconce2_nearest10" ; do # note: these should be changed to be unique identifiers per run of sconce2 (ie needed to find output files)
      # check if last cell in this tumor_depths_$c set has a mode file (ie is this run done)
      lastCell=$(basename "$(tail -n 1 "$p""/tumor_depths_""$c")")
      lastCellMode="$p""/output_""$fileKey""_""$dataset""_k""$k""_c""$c""__""$lastCell""__k""$k""__mode.bed"
      if [[ ! -f "$lastCellMode" ]] ; then
        echo "could not find last cell's mode bed file: $lastCellMode. skipping"
        continue
      fi
      for bedType in "sconce" "mean" "median" "mode" ; do
        outBase="$p""/output_""$fileKey""_""$dataset""_k""$k""_c""$c""_rounded""$bedType"
        echo "going to run ""$outBase"

        # create fastas
        ./cnaBedToFasta.sh "$p" "$k" "$c" "$bedType" "$fileKey" "F" > "$outBase"".fasta"
        ./cnaBedToFasta.sh "$p" "$k" "$c" "$bedType" "$fileKey" "T" > "$outBase""Rev.fasta"

        # forward
        python "$cnp2cnpDir""/cnp2cnp.py" -m "matrix" -d any -i "$outBase"".fasta" -o "$outBase"".cnp2cnp"
        python "$cnp2cnpDir""/cnp2cnp.py" -m "matrix" -d zzs -i "$outBase"".fasta" -o "$outBase"".zzs.cnp2cnp"

        # reversed. calc both since cnp2cnp is not symmetric
        python "$cnp2cnpDir""/cnp2cnp.py" -m "matrix" -d any -i "$outBase""Rev.fasta" -o "$outBase""Rev.cnp2cnp"
        python "$cnp2cnpDir""/cnp2cnp.py" -m "matrix" -d zzs -i "$outBase""Rev.fasta" -o "$outBase""Rev.zzs.cnp2cnp"

      done # output bed file type loop
    done # fileKey loop
  done # cell subset loop
done # parameter set loop

