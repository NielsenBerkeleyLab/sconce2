library(ggplot2)
library(stringr)

# script to calculate the median euclidean distances between each SCONCE profile for each dataset, for the 1st, 10th, and 20th nearest cells

k <- 10
numCellsList <- c(20, 40, 60, 80, 100, 120)
filekeys <- c("sconce2_nearest10")
paramSets <- c("paramsA", "paramsB", "paramsC", "paramsD")
forceRecalc <- T

dataDir <- "."
outputDir <- paste0(dataDir, "plots/")
if(!dir.exists(outputDir)) {
  dir.create(outputDir)
}

plotWidth <- 8
plotHeight <- 5.5



calcPairwiseDistMatSummary <- function(logFile, tumorDepths, summaryFile, forceRecalc=F) {
  if(!forceRecalc && file.exists(summaryFile)) {
    summaryDat <- read.table(summaryFile, header=T)
    return(summaryDat)
  }
  distMat <- read.table(text=system(paste0("grep -A 20 \"this->stage1PairwiseDistMat\" ", logFile, " | tail -n+2 | cut --complement -f 1 | sed 's/\t$//'"), intern=T), sep="\t")

  summaryDat <- do.call(rbind, lapply(colnames(distMat), FUN=function(col) {
    currDat <- distMat[,col]
    currDat <- currDat[currDat != 0]
    sorted <- sort(currDat)
    data.frame(first=sorted[1], tenth=sorted[10], last=sorted[length(sorted)])
  }))
  depthFiles <- read.table(tumorDepths)
  cellNames <- sapply(depthFiles, FUN=function(x) {str_extract(x, "cancer_cell_[0-9]*.")})
  rownames(summaryDat) <- cellNames
  write.table(summaryDat, file=summaryFile, quote=F, sep="\t", row.names=T, col.names=T)
  summaryDat
}

summaryDistMatDatList <- list()
combinedSummaryDistDatList <- list()
for(paramSet in paramSets) {
  for(key in filekeys) {
    for(numCells in numCellsList) {
      currOutputDir <- paste0(dataDir, "/", paramSet, "/plots/")
      if(!dir.exists(currOutputDir)) {
        dir.create(currOutputDir)
      }
      datasetName <- paste0(paramSet, ", ", key, ", ", numCells)

      # read pairwise distance matrices
      print(paste0("reading dist mat ", datasetName))
      logFile <- system(paste0("find ", dataDir, "/", paramSet, " -maxdepth 1 -name \"output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, ".log\""), intern=T) # based on scAllP_*sh outBase variable
      tumorDepths <- system(paste0("find ", dataDir, "/", paramSet, " -maxdepth 1 -name \"tumor_depths_", numCells, "\""), intern=T)
      summaryFile <- paste0(currOutputDir, "output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, "_pairwiseDistMatSummary") # based on scAllP_*sh outBase variable
      summaryDistMatDat <- calcPairwiseDistMatSummary(logFile, tumorDepths, summaryFile, forceRecalc)
      summaryDistMatDatList[[datasetName]] <- summaryDistMatDat
    }
    # summarize across cells in this dataset
    combinedSummaryDistDatList[[paramSet]] <- do.call(rbind, summaryDistMatDatList[grepl(paramSet, names(summaryDistMatDatList))])
  }
}

# calc median of each paramSet's first, tenth, and last cols
medianDists <- data.frame(do.call(rbind, lapply(paramSets, FUN=function(paramSet) {sapply(combinedSummaryDistDatList[[paramSet]], median)})))
rownames(medianDists) <- paramSets

# write to text and latex files
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_combinedPairwiseDistMatSummary")
write.table(medianDists, file=paste0(outputDir, "/", outputFile, ".txt"), quote=F, sep="\t", row.names=T, col.names=NA)

texColnames <- "& \\textbf{nearest/most similar} & \\textbf{10th nearest} & \\textbf{furthest/least similar}  \\\\ \\hline"
medianDistsTex <- medianDists
medianDistsTex$paramSet <- paste0(" \\\\ \\hline % ", rownames(medianDistsTex))
write.table(c(texColnames, paste0(paste0(c("A", "B", "C", "D"), " & ", sep=""), apply(medianDistsTex[,-ncol(medianDistsTex)], 1, FUN=function(x) {paste0(x, collapse=" & ")}),  medianDistsTex[,ncol(medianDistsTex)], sep="")), file=paste0(outputDir, "/", outputFile, ".tex"), quote=F, row.names=F, col.names=F)
paste0(paste0(c("A", "B", "C", "D"), " & ", sep=""), apply(medianDistsTex[,-ncol(medianDistsTex)], 1, FUN=function(x) {paste0(x, collapse=" & ")}),  medianDistsTex[,ncol(medianDistsTex)], sep="")

