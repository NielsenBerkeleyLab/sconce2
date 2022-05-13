library(ape)
library(plyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(cowplot)

source("readTreeBranches.R")

# script to make plots of correlations of true branch lengths (from newick trees) and estimated tree branch lengths. Using correlation bc tree branches are scaled differently

dataDir <- "."
outputDir <- paste0(dataDir, "plots/")
if(!dir.exists(outputDir)) {
  dir.create(outputDir)
}
k <- 10
numCellsList <- c(20, 40, 60, 80, 100, 120)
filekeys <- c(
              "sconce2")
              #"sconce2_nearest10")
paramSets <- c("paramsA", "paramsB", "paramsC", "paramsD")
forceRecalc <- F

treeBranchList <- list()
corrPlotList <- list()

for(numCells in numCellsList) {
  for(paramSet in paramSets) {
    for(key in filekeys) {
      shortName <- paste0(gsub("params", "p", paramSet), "_", key, "_c", numCells)
      print(paste0("reading ", shortName))
      newickString <- newickStrings[[paramSet]]
      branchLength <- branchLengths[[paramSet]]
      hmmFile <- system(paste0("find ", dataDir, "/", paramSet, " -name \"output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, ".hmm\""), intern=T) # based on scAllP_*sh outBase variable
      if(length(hmmFile) == 0) {
        next
      }
      mergedTreeBranches <- getTreeBranches(newickString, branchLength, hmmFile, forceRecalc)
      if(is.null(mergedTreeBranches)) {
        next
      }
      treeBranchList[[shortName]] <- mergedTreeBranches

      print(paste0("plotting ", shortName))
      p <- makeCorrPlot(mergedTreeBranches, shortName)
      corrPlotList[[shortName]] <- p
    }
  }
}

# plot combining cell subsets across each paramSet
combinedTreeBranchList <- list()
combinedCorrPlotList <- list()
for(paramSet in paramSets) {
  combinedTreeBranchList[[paramSet]] <- do.call(rbind, treeBranchList[grepl(gsub("params", "p", paramSet), names(treeBranchList))])
  combinedCorrPlotList[[paramSet]] <- makeCorrPlot(combinedTreeBranchList[[paramSet]], NA)
}
legend <- get_legend(combinedCorrPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom"))
pGrid <- plot_grid(plotlist=lapply(combinedCorrPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=2)
toSave <- pGrid
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_combinedtreeBranchCorr")

png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth, height=plotHeight, res=600, units="in"); plot(toSave); dev.off()
pdf(paste0(outputDir, "/", outputFile, ".pdf"), width=plotWidth, height=plotHeight); plot(toSave); dev.off()

