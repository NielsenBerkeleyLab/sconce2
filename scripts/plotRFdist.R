library(ggtree)
library(phangorn)
library(ape)
library(reshape2)
library(stringr)
library(ggplot2)
library(cowplot)
library(scales)

source("readTreeBranches.R")

# script to calculate the similarity/distance between estimated trees
# used for comparing nj on t2+t3 distances and some other tree building metric
# plots Robinson-Foulds distances and neighbor joining trees

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
treeListList <- list()
treeDistList <- list()
njTreesPlotListList <- list()
rfPlotList <- list()

for(numCells in numCellsList) {
  for(paramSet in paramSets) {
    for(key in filekeys) {
      shortName <- paste0(gsub("params", "p", paramSet), "_", key, "_c", numCells)
      print(paste0("reading ", shortName))
      newickString <- newickStrings[[paramSet]]
      branchLength <- branchLengths[[paramSet]]
      hmmFile <- system(paste0("find ", dataDir, "/", paramSet, " -maxdepth 1 -name \"output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells, ".hmm\""), intern=T) # based on scAllP_*sh outBase variable
      if(length(hmmFile) == 0) {
        next
      }
      mergedTreeBranches <- getTreeBranches(newickString, branchLength, hmmFile, forceRecalc)
      if(is.null(mergedTreeBranches)) {
        next
      }
      treeBranchList[[shortName]] <- mergedTreeBranches

      outputFile <- paste0(dataDir, "/plots/njTrees_", gsub("/", "_", paramSet), "_", key, "_k", k, "_c", numCells, "_cnp.png")
      treeList <- createNJtrees(mergedTreeBranches, paramSet, numCells, key)
      if(is.null(treeList)) {
        next
      }
      treeListList[[shortName]] <- treeList

      njTreesPlotList <- makeNJtreePlots(treeList, outputFile)
      njTreesPlotListList[[shortName]] <- njTreesPlotList

      dists <- calcRFdists(treeList)
      # if nearest10, remove t2_t3 from plots (since there might be missing pairs)
      if(grepl("nearest10", key, ignore.case=T)) {
        dists <- subset(dists, variable != "t2_t3")
      }
      dists$paramSet <- paramSet
      treeDistList[[shortName]] <- dists

      p <- makeRFdistPlot(dists, shortName)
      rfPlotList[[shortName]] <- p
    }
  }
}

# plot combining each subset of cells within a paramset
print("combining cell subsets across paramsets")
combinedTreeDistList <- list()
combinedTreeDistPlotList <- list()
combinedTreeDistFiltPlotList <- list()
combinedTreeDistFacetPlotList <- list()
print("plotting cell subsets across paramsets")
for(paramSet in paramSets) {
  combinedTreeDistList[[paramSet]] <- do.call(rbind, treeDistList[grepl(gsub("params", "p", paramSet), names(treeDistList))])
  combinedTreeDistPlotList[[paramSet]] <- makeRFdistPlot(combinedTreeDistList[[paramSet]], NA)
}
legend <- get_legend(combinedTreeDistPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom") + guides(colour=guide_legend(nrow=1)))
pGrid <- plot_grid(plotlist=lapply(combinedTreeDistPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=2)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_combinedRfDistWithZZS")

# small plot
# euc, cnp2, zzs on sconce, t2+t3
print("plotting just euc/cnp/zzs on sconce, t2+t3")
for(paramSet in paramSets) {
  combinedTreeDistFiltPlotList[[paramSet]] <- makeRFdistPlot(subset(combinedTreeDistList[[paramSet]], variable %in% c("eucSconce", "cnpSconce", "zzsSconce", "t2_t3")), NA)
}
legend <- get_legend(combinedTreeDistFiltPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom") + guides(colour=guide_legend(nrow=1)))
pGrid <- plot_grid(plotlist=lapply(combinedTreeDistFiltPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=2)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_combinedRfDistWithZZS_filt")

png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth, height=plotHeight, res=600, units="in"); plot(toSave); dev.off()
pdf(paste0(outputDir, "/", outputFile, ".pdf"), width=plotWidth, height=plotHeight); plot(toSave); dev.off()

# faceted plot: sections/facets are distance metric, true/sconce/mean/median/mode[t2+t3] are box plots
print("plotting faceted plot")
for(paramSet in paramSets) {
  p <- makeFacetedRFdistPlot(combinedTreeDistList[[paramSet]], NA)
  p <- p + theme(strip.text.x = element_text(size = 6)) # shrink facet labels
  combinedTreeDistFacetPlotList[[paramSet]] <- p
}
legend <- get_legend(combinedTreeDistFacetPlotList[[1]] + theme(legend.box.margin=margin(0, 0, 0, 0), legend.position="bottom") + guides(colour=guide_legend(nrow=1)))
pGrid <- plot_grid(plotlist=lapply(combinedTreeDistFacetPlotList, FUN=function(x) {x + theme(legend.position="none")}), align='vh', labels="AUTO", nrow=2)
toSave <- plot_grid(pGrid, legend, ncol=1, rel_heights=c(1, 0.05))
outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_", filekeys[1], "-", filekeys[length(filekeys)], "_k", k, "_c", paste0(numCellsList, collapse="-c"), "_combinedRfDistWithZZS_facet")

png(paste0(outputDir, "/", outputFile, ".png"), width=plotWidth, height=plotHeight, res=600, units="in"); plot(toSave); dev.off()
pdf(paste0(outputDir, "/", outputFile, ".pdf"), width=plotWidth, height=plotHeight); plot(toSave); dev.off()
#save_plot(paste0(outputDir, "/", outputFile, ".eps"), toSave, device=cairo_ps, dpi=600, base_width=plotWidth, base_height=plotHeight)

