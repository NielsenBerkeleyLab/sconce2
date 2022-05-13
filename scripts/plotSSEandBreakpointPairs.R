library(reshape2)
library(stringr)
library(ggplot2)
library(cowplot)

# script to plot SSE and breakpoint detection accuracy on paired data

# will create countedBreakpointsTable*txt, *sumSqSconceOnePairMeanMedianMode_indv.txt, and *sumSqSconceOnePairMeanMedianMode_pairs.txt if they don't exist, and will read from file if they do

source("readBedFilesPairs.R")

dataDir <- "."
k <- 10
numCellsList <- c(20, 40, 60, 80, 100, 120)
filekeys <- c(
              "sconce2")
              #"sconce2_nearest10")
paramSets <- c("paramsA", "paramsB", "paramsC", "paramsD")
forceRecalc <- F
inclAneu <- T

breakpointDatList <- list()
breakpointPlotList <- list()
sseDatList <- list()
ssePlotList <- list()

for(numCells in numCellsList) {
  for(paramSet in paramSets) {
    for(key in filekeys) {
      currOutputDir <- paste0(dataDir, "/", paramSet, "/plots/")
      if(!dir.exists(currOutputDir)) {
        dir.create(currOutputDir)
      }
      datasetName <- paste0(paramSet, ", ", key, ", ", numCells)
      print(paste0("reading SSE ", datasetName))
      sseFile <- paste0(currOutputDir, "output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, "_sumSqSconceOnePairMeanMedianMode") # based on scAllP_*sh outBase variable
      sseDat <- calcAllSSE(paramSet, key, sseFile, numCells, forceRecalc, inclAneu)
      sseDatList[[datasetName]] <- sseDat

      print(paste0("reading breakpoint ", datasetName))
      breakpointFile <- paste0(currOutputDir, "output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, "_breakpointDistSconceOnePairMeanMedianMode") # based on scAllP_*sh outBase variable
      breakpointDat <- calcAllBreakpoints(paramSet, key, breakpointFile, numCells, forceRecalc, inclAneu)
      breakpointDatList[[datasetName]] <- breakpointDat

      print(paste0("saving tex files ", datasetName))
      medianDistOmegaTexFile <- paste0(currOutputDir, "output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, "_medianDistOmegaSconceOnePairMeanMedianMode") # based on scAllP_*sh outBase variable
      breakpointDat$paramSet <- paramSet
      writeMedianDistOmegaTexFiles(breakpointDat, medianDistOmegaTexFile)
    }
  }
  for(paramSet in paramSets) {
    for(key in filekeys) {
      datasetName <- paste0(paramSet, ", ", key, ", ", numCells)
      shortName <- paste0(gsub("params", "p", paramSet), "_", key, "_c", numCells)
      if(datasetName %in% names(sseDatList)) {
        print(paste0("plotting SSE ", datasetName))
        ssePlot <- makeSSEPlot(sseDatList[[datasetName]], shortName)
        ssePlotList[[datasetName]] <- ssePlot
      }
      if(datasetName %in% names(breakpointDatList)) {
        print(paste0("plotting breakpoint ", datasetName))
        breakpointPlot <- makeBreakpointPlot(breakpointDatList[[datasetName]], shortName)
        breakpointPlotList[[datasetName]] <- breakpointPlot
      }
    }
  }
}

# plot combining each subset of cells within a paramset
combinedBreakpointDatList <- list()
combinedBreakpointPlotList <- list()
combinedSSEDatList <- list()
combinedSSEPlotList <- list()
print("combining all numCells within each paramSet")
for(paramSet in paramSets) {
  combinedBreakpointDatList[[paramSet]] <- do.call(rbind, breakpointDatList[grepl(paramSet, names(breakpointDatList))])
  combinedBreakpointDatList[[paramSet]]$paramSet <- paramSet
  combinedBreakpointPlotList[[paramSet]] <- makeBreakpointPlot(combinedBreakpointDatList[[paramSet]], NA)
  combinedSSEDatList[[paramSet]] <- do.call(rbind, sseDatList[grepl(paramSet, names(sseDatList))])
  combinedSSEPlotList[[paramSet]] <- makeSSEPlot(combinedSSEDatList[[paramSet]], NA)
}
print("plotting SSE and breakpoint combined numCells within each paramSet")
legend <- get_legend(combinedBreakpointPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom") + guides(colour=guide_legend(nrow=1)))
pGrid <- plot_grid(plotlist=lapply(c(combinedSSEPlotList, combinedBreakpointPlotList), FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=2)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_combinedSSEbreakpointDistSconceMeanMedianMode")

png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth, height=plotHeight, res=600, units="in"); plot(toSave); dev.off()
pdf(paste0(outputDir, "/", outputFile, ".pdf"), width=plotWidth, height=plotHeight); plot(toSave); dev.off()

print("plotting SSE combined numCells within each paramSet")
pGrid <- plot_grid(plotlist=lapply(combinedSSEPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=2)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_combinedSSESconceMeanMedianMode")

png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth, height=plotHeight, res=600, units="in"); plot(toSave); dev.off()
pdf(paste0(outputDir, "/", outputFile, ".pdf"), width=plotWidth, height=plotHeight); plot(toSave); dev.off()

print("plotting breakpoint combined numCells within each paramSet")
pGrid <- plot_grid(plotlist=lapply(combinedBreakpointPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=2)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_combinedBreakpointDistSconceMeanMedianMode")

png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth, height=plotHeight, res=600, units="in"); plot(toSave); dev.off()
pdf(paste0(outputDir, "/", outputFile, ".pdf"), width=plotWidth, height=plotHeight); plot(toSave); dev.off()

print("saving tex files combined numCells within each paramSet")
allBreakpointDat <- do.call(rbind, combinedBreakpointDatList)
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_combinedMedianDistOmegaSconceOnePairMeanMedianMode")
writeMedianDistOmegaTexFiles(allBreakpointDat, paste0(outputDir, "/", outputFile))

