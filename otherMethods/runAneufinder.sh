#!/bin/bash

# Tue 22 Feb 2022 06:08:28 PM PST

# example usage: echo -e "paramsA\nparamsB\nparamsC\nparamsD" | parallel -j 2 --delay 2 ./runAneufinder.sh {}

export p="$1"

export scripts="../scripts"
export ref="../../hg19/hg19_lite.fa"
export windows="../../hg19/hg19_lite.250kb_windows"
export OMP_NUM_THREADS=1

# aneufinder
aneufinderDir=$p/aneu
mkdir -p $aneufinderDir

{ time Rscript --vanilla "$scripts"/runAneufinder_noGCmap.R $p $aneufinderDir &> $aneufinderDir/aneufinderTiming.log ; } 2> $aneufinderDir/aneufinderTiming.total.log

# aneufinder needs to be split from segments to bins
while read cellID ; do
  aneuSegs=$(find $aneufinderDir -name "$cellID.aneu.bed") # simulations
  aneuBins="$aneufinderDir""/""$cellID"".aneufinderBins.bed"
  bedtools intersect -a "$aneuSegs" -b "$windows" | sort -k1,1V -k2,2n > "$aneuBins"
done < $p/tumor_ids_hg19

