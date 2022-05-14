# SCONCE2

This directory is for the program SCONCE2.

Preprint at https://www.biorxiv.org/content/10.1101/2022.05.12.491742v1

## Dependencies
SCONCE2 is written in C++11, and requires
- GNU make (tested on v4.1)
- g++ (tested on 7.5.0)
- BOOST (tested on v1.66.0)
- GSL (tested on v2.4)

Additional [R scripts](scripts/) require
- R (tested on v4.2.0)
- ape (tested on v5.6-2)
- cowplot (tested on v1.1.1)
- ggplot2 (tested on v3.3.5)
- ggtree (tested on v3.4.0)
- grid (tested on v0.5-1)
- gtools (tested on v3.9.2)
- phangorn (tested on v2.8.1)
- plyr (tested on v1.8.7)
- reshape2 (tested on v1.4.4)
- scales (tested on v1.2.0)
- stringr (tested on v1.4.0)

SCONCE2 was developed and tested on Ubuntu 18.04.

## Installation instructions
1. Clone this repo:
```
git clone git@github.com:NielsenBerkeleyLab/sconce2.git
```
2. Run `make`. This will build intermediates into the `build/` directory and create an executable named `sconce2`.


## Brief parameter descriptions for SCONCE2
- `--diploid` This file should be the averaged read depth across all diploid cells. It can be generated using [scripts/avgDiploid.R](scripts/avgDiploid.R).
- `--tumorFileList` This file should list filenames of the per window read depth of tumor cells to analyze.
- `--meanVarCoefFile` This file should define the coefficients for the relationship between the mean and variance of the negative binomial distribution used for emission probabilities. It should be generated using [scripts/fitMeanVarRlnshp.R](scripts/fitMeanVarRlnshp.R).
- `--outputBase` This gives the basename for all [output files](#output-files).
- `--maxKploid` This gives the maximum allowed ploidy (recommended `k=10`).
- `-j` This gives how many threads sconce2 should use (default 1)
- `--sconceEstimatesPath` This is useful if one needs to stop and resume an analysis, or if one has previously run sconce and would like to skip the first sconce step. This should be the `--outputBase` value of the previous run. If previous runs cannot be found, sconce2 will be run.
- `--pairedEstimatesPath` as above, but for the paired estimates. Should be the `--outputBase` value of the previous run.
- Run `./sconce2 -h` for the full list of options


## Input files
Averaged diploid read depth files should be tab separated, with columns `<chr>\t<start>\t<end>\t<meanReadDepth>\t<varianceOfReadDepth>`. They should be generated using [scripts/avgDiploid.R](scripts/avgDiploid.R), given a file providing a list of paths to the observed diploid read depths. For example:
```
Rscript scripts/avgDiploid.R test/diploidFileList test/test_healthy_avg.bed
```
produces the following output:
```
$ head test/ref_healthy_avg.bed
chr1	0	250000	314.96	2413.43272727273
chr1	250000	500000	318.9	2269.84848484848
chr1	500000	750000	322.03	2715.32232323232
chr1	750000	1000000	321.58	1781.64
chr1	1000000	1250000	318.57	2367.21727272727
chr1	1250000	1500000	318.14	2574.24282828283
chr1	1500000	1750000	327.14	2439.11151515152
chr1	1750000	2000000	313.13	1945.95262626263
chr1	2000000	2250000	329.12	2260.00565656566
chr1	2250000	2500000	321.15	2663.78535353535
```

`tumorFileList` should list tumor read depth files to analyze, one per line. Tumor read depth files should be tab separated, with columns `<chr>\t<start>\t<end>\t<readDepth>`. See [simulations/README.md](simulations/README.md) for how to generate simulations with this format. For real data, a tool like [bedtools coverage](https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html) can be used to create this file from a bam file. For example:
```
$ cat test/tumorFileList
test/simu_cancer_cell_0.hg19_lite.bed
test/simu_cancer_cell_1.hg19_lite.bed
test/simu_cancer_cell_2.hg19_lite.bed
test/simu_cancer_cell_3.hg19_lite.bed

$ head test/simu_cancer_cell_0.hg19_lite.bed
chr1	0	250000	348
chr1	250000	500000	401
chr1	500000	750000	307
chr1	750000	1000000	349
chr1	1000000	1250000	259
chr1	1250000	1500000	337
chr1	1500000	1750000	300
chr1	1750000	2000000	295
chr1	2000000	2250000	362
chr1	2250000	2500000	350
```

Mean and variance coefficient files should have one parameter and value pair per line. They should be generated using [scripts/fitMeanVarRlnshp.R](scripts/fitMeanVarRlnshp.R). For example:
```
Rscript scripts/fitMeanVarRlnshp.R test/diploidFileList test/ref.meanVar
```
produces the following output:
```
$ cat test/ref.meanVar
intercept=10.79950512968
slope= 1.18499408815
poly2= 0.01910756218
```

## Example test run
To ensure `sconce2` was built and the above scripts were run correctly, we include some test files. Run the following:
```
mkdir -p test && cd test
export windows="../../hg19/hg19_lite.250kb_windows"
cp ../simulations/inputFiles/infileA.txt . && cp ../simulations/inputFiles/paramfile.txt .
../simulations/sconce_sim infileA.txt paramfile.txt 0 > test_sim.log 2>&1
mkdir -p diploid
mkdir -p cancer
mv *healthy* diploid
find . -maxdepth 1 -name "*_cancer_cell_[0-9]*" -exec mv {} cancer \;
find . -name "*_cell_*[0-9]" | while read simu ; do
  depth="$simu"".hg19_lite.bed"
  paste "$windows" <(awk '{print $4}' $simu) | awk 'BEGIN{OFS="\t"} {if($4 == "") {$4=0} print}' > $depth
done

find diploid -name "simu_healthy*hg19_lite.bed" | sort -V > diploidFileList
find cancer -name "simu_cancer*hg19_lite.bed" | sort -V | head -n 4 > tumorFileList

Rscript --vanilla ../scripts/avgDiploid.R diploidFileList test_healthy_avg.bed
Rscript --vanilla ../scripts/fitMeanVarRlnshp.R diploidFileList test.meanVar

../sconce2 -d test_healthy_avg.bed -t tumorFileList -k 5 -v --meanVarCoefFile test.meanVar -o test_k5_c4 --saveSconce --summarizeAll -j 6 --sconceEstimatesPath test_k5_c4 --pairedEstimatesPath test_k5_c4 > test_k5_c4.log 2> test_k5_c4.err
```
Your output (with the exception of timing information and randomized ordering due to multithreading) should match the provided `test/ref*` files.


## Output files
SCONCE2 will create the following files automatically:
- `<output>.hmm` This file contains the state of the HMM after the Baum Welch step and the state of the HMM after the BFGS step.
- `<output>__<cellName>__k<maxKploid>.sconceParams` This file contains the model estimates from running sconce on this individual cell. Useful for resuming an analysis.
- `<output>__<cell0Name>__<cell1Name>__k<maxKploid>.pairedParams` This file contains the model estimates from running cell0 and cell1 as a pair. Useful for resuming an analysis.
- `<output>__sconce__<cellName>__k<maxKploid>.bed` This file contains the copy number calls in tab separated bed format for this individual cell from running sconce. Only saved if the `--saveSconce` option is passed.
- `<output>__pair_<cell0Name>_<cell1Name>__<cellXName>__k<maxKploid>.bed` This file contains the copy number calls in tab separated bed format for cellX (out of the pair cell0 and cell1), created from the joint/paired analysis of cell0 and cell1.
test_k5_c4__pair_simu_cancer_cell_2.hg19_lite.bed_simu_cancer_cell_3.hg19_lite.bed__cell_simu_cancer_cell_3.hg19_lite.bed__k5.bed
- `<output>__<cellName>__k<maxKploid>__<summaryMethod>.bed` This file contains the summary copy number calls in tab separated bed format for this indivdual cell, created using <summaryMethod>. Only mean is saved by default, median and mode can be enabled by passing `--summarizeAll` or `--summarizeMedian`/`--summarizeMode`, respectively.


SCONCE2 also prints log messages to stdout and error messages to stderr.
If the `--verbose` flag is used, debugging statements will be printed to stderr.

## Simulations
To compile and run the simulation program, see [simulations/README.md](simulations/README.md).

