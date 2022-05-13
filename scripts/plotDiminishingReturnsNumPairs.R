library(ggplot2)
library(gtools)
library(stringr)
library(cowplot)

# Sun 27 Feb 2022 06:33:02 PM PST
# makes plots of how change in SSE (relative to SCONCE) drops as you add cells, but there's a point of diminishing returns

source("readBedFilesPairs.R")
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
selectionMethods <- c("nearest", "random", "furthest")

set.seed(0)
dimReturnsDatList <- list()
dimReturnsPlotList <- list()
for(numCells in numCellsList) {
  for(paramSet in paramSets) {
    for(key in filekeys) {
      currOutputDir <- paste0(dataDir, "/", paramSet, "/plots/")
      if(!dir.exists(currOutputDir)) {
        dir.create(currOutputDir)
      }
      shortName <- paste0(gsub("params", "p", paramSet), "_", key, "_c", numCells)
      print(paste0("reading ", shortName))

      dimReturnsFile <- paste0(currOutputDir, "output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, "_diminishingReturnsChangeSSE") # based on scAllP_*sh outBase variable
      dimReturns <- calcDiminishingReturns(paramSet, key, dimReturnsFile, numCells, forceRecalc)
      if(is.null(dimReturns)) {
        next
      }
      dimReturnsDatList[[shortName]] <- dimReturns

      p <- makeDiminishingReturnsPlot(dimReturns, shortName)
      dimReturnsPlotList[[shortName]] <- p
    }
  }
}

# plot combining each subset of cells within a paramset
combinedDimReturnsDatList <- list()
combinedDimReturnsPlotList <- list()
summarizedImprList <- list()
for(paramSet in paramSets) {
  combinedDimReturnsDatList[[paramSet]] <- do.call(rbind, dimReturnsDatList[grepl(gsub("params", "p", paramSet), names(dimReturnsDatList))])
  combinedDimReturnsPlotList[[paramSet]] <- makeDiminishingReturnsPlot(combinedDimReturnsDatList[[paramSet]], NA)
  summarizedImprList[[paramSet]] <- summarizeDiminishingReturns(combinedDimReturnsDatList[[paramSet]])
  summarizedImprList[[paramSet]]$paramSet <- paramSet
}
legend <- get_legend(combinedDimReturnsPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom"))
pGrid <- plot_grid(plotlist=lapply(combinedDimReturnsPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=2)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_combinedDiminishingReturnChangeSSE")

png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth, height=plotHeight, res=600, units="in"); plot(toSave); dev.off()
pdf(paste0(outputDir, "/", outputFile, ".pdf"), width=plotWidth, height=plotHeight); plot(toSave); dev.off()

# write a table of the mean change in SSE
summarizedImpr <- do.call(rbind, summarizedImprList)

# write a table of minimum/biggest drop in SSE, drop at nearest10, and drop when using all 20 cells
write.table(summarizedImpr, paste0(outputDir, "/", outputFile, ".txt"), sep="\t", col.names=T, row.names=F, quote=F)
nearestExtremes <- do.call(rbind, lapply(paramSets, FUN=function(pSet) {
  currDat <- subset(summarizedImpr, paramSet == pSet & selectionMethod == "nearest")
  data.frame(paramSet=pSet, impr=c(min(currDat$meanImpr), currDat$meanImpr[10], currDat$meanImpr[length(currDat$meanImpr)]), numCellsSummarized=c(which.min(currDat$meanImpr), 10, length(currDat$meanImpr)), type=c("min", "nearest10", "all"))
}))
write.table(nearestExtremes, paste0(outputDir, "/", outputFile, "_nearestExtremes.txt"), sep="\t", col.names=T, row.names=F, quote=F)

