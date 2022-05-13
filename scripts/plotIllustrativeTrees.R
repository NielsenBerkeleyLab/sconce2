library(phangorn)
library(ape)
library(reshape2)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggtree)

source("readTreeBranches.R")

# script to plot tiny illustration trees of 8 cells per paramSet

numCellsList <- 8
paramSets <- c("paramsA", "paramsB", "paramsC", "paramsD")

for(numCells in numCellsList) {
treeList <- list()
treePlotList <- list()
for(paramSet in paramSets) {
  newickString <- newickStrings[paramSet]
  treeObj <- read.tree(text=newickString)
  smallTree <- keep.tip(treeObj, 1:numCells)
  treeList[[paramSet]] <- smallTree

  treePlotList[[paramSet]] <- ggtree(treeList[[paramSet]], ladderize=F) + geom_tiplab()
}

outputFile <- paste0(gsub("/", "_", paramSets[1]), "-", gsub("/", "_", paramSets[length(paramSets)]), "_c", paste0(numCellsList, collapse="-c"), "_illustrativeTrees")
toSave <- plot_grid(plotlist=treePlotList, align='vh', labels="AUTO", nrow=2)
png(paste0(outputDir, "/", outputFile, ".png"), width=4, height=4, res=600, units="in"); plot(toSave); dev.off()
pdf(paste0(outputDir, "/", outputFile, ".pdf"), width=4, height=4); plot(toSave); dev.off()
}

