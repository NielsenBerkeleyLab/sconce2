library(ggplot2)
library(reshape2)
library(scales)
library(grid)
library(stringr)

# Thu 03 Mar 2022 03:51:17 PM PST
# script to plot specific indices of a pair of cells to show that using a pair of cells has lower SSE than using sconce due to
#   1. summarizing across multiple cells/inference runs
#   2. even in one pair, get more info using just a pair

paramSet <- "paramsC"
numCells <- 100
cellNumA <- 111
minIdx <- 8280
maxIdx <- 8390

cellNumB <- 69
arrows <- data.frame(x=   c(8347,8324,8355,8317,8355,8317,8355,8317,8355,    8315,8324,8317,8317,8317),
                     xend=c(8347,8324,8355,8317,8355,8317,8355,8317,8355,    8315,8324,8317,8317,8317),
                     y=rep(8.4,14),
                     yend=rep(7.5,14),
                     source=c("sconce", "pair", "pair", "mean", "mean", "median", "median", "mode", "mode", "sconce", "pair", "mean", "median", "mode"),
                     cell=c(rep(paste0("cancer_cell_", cellNumA, "."),9), rep(paste0("cancer_cell_", cellNumB, "."),5)))

cellNumB <- 59
arrows <- data.frame(x=   c(8347,8316,8355,8317,8355,8317,8355,8317,8355,    8311,8316,8317,8317,8317),
                     xend=c(8347,8316,8355,8317,8355,8317,8355,8317,8355,    8311,8316,8317,8317,8317),
                     y=rep(9.5,14),
                     yend=rep(8,14),
                     source=c("sconce", "pair", "pair", "mean", "mean", "median", "median", "mode", "mode", "sconce", "pair", "mean", "median", "mode"),
                     cell=c(rep(paste0("cancer_cell_", cellNumA, "."),9), rep(paste0("cancer_cell_", cellNumB, "."),5)))



k <- 10
dataDir <- "."
outputDir <- paste0(dataDir, "plots/")
filekeys <- c("sconce2")
key <- filekeys

plotRoot <- paste0(dataDir, "/plots/betterBounds_cells_", cellNumA, "_", cellNumB, "_", gsub("/", "_", paramSet), "_", key, "_k", k, "_c", numCells)

header <- c("chr", "start", "end")
cellIDs <- paste0("cancer_cell_", c(cellNumA, cellNumB), ".")
names(cellIDs) <- c("cellNameA", "cellNameB")

simuDatA <- read.table(system(paste0("find ", dataDir, "/", paramSet, " -name \"*simu_cancer_cell_", cellNumA, ".hg19_lite.depth\""), intern=T), stringsAsFactors=F, header=F)
simuDatB <- read.table(system(paste0("find ", dataDir, "/", paramSet, " -name \"*simu_cancer_cell_", cellNumB, ".hg19_lite.depth\""), intern=T), stringsAsFactors=F, header=F)

trueDatA <- read.table(system(paste0("find ", dataDir, "/", paramSet, " -name \"true_cancer_cell_", cellNumA, ".hg19_lite.bed\""), intern=T), stringsAsFactors=F, header=F)
trueDatB <- read.table(system(paste0("find ", dataDir, "/", paramSet, " -name \"true_cancer_cell_", cellNumB, ".hg19_lite.bed\""), intern=T), stringsAsFactors=F, header=F)

outputFilePrefix <- paste0("output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells) # based on scAllP_*sh outBase variable
sconceDatA <- read.table(system(paste0("find ", dataDir, "/", paramSet, " -maxdepth 1 -name \"", outputFilePrefix, "__sconce__simu_cancer_cell_", cellNumA, ".hg19_lite.depth*.bed\""), intern=T), stringsAsFactors=F, header=F)
sconceDatB <- read.table(system(paste0("find ", dataDir, "/", paramSet, " -maxdepth 1 -name \"", outputFilePrefix, "__sconce__simu_cancer_cell_", cellNumB, ".hg19_lite.depth*.bed\""), intern=T), stringsAsFactors=F, header=F)

unfiltPairsFileList <- system(paste0("find ", dataDir, paramSet,  " -maxdepth 1 -name \"", outputFilePrefix, "__pair_*.bed\" | sort -V"), intern=T)
pairsFileList <- unfiltPairsFileList[sapply(lapply(str_extract_all(unfiltPairsFileList, "cancer_cell_[0-9]*."), unique), FUN=function(cellNames) {cellNames[1] %in% cellIDs && cellNames[2] %in% cellIDs})]

pairFileA <- pairsFileList[grepl(paste0("simu_cancer_cell_", cellNumA, ".hg19_lite.depth__k"), pairsFileList)]
pairFileB <- pairsFileList[grepl(paste0("simu_cancer_cell_", cellNumB, ".hg19_lite.depth__k"), pairsFileList)]
pairDatA <- read.table(pairFileA, stringsAsFactors=F, header=F)
pairDatB <- read.table(pairFileB, stringsAsFactors=F, header=F)

meanDatA <- read.table(system(paste0("find ", dataDir, paramSet,   " -maxdepth 1 -name \"", outputFilePrefix, "__simu_cancer_cell_", cellNumA, ".hg19*__mean.bed\""), intern=T), stringsAsFactors=F, header=F)
meanDatB <- read.table(system(paste0("find ", dataDir, paramSet,   " -maxdepth 1 -name \"", outputFilePrefix, "__simu_cancer_cell_", cellNumB, ".hg19*__mean.bed\""), intern=T), stringsAsFactors=F, header=F)

medianDatA <- read.table(system(paste0("find ", dataDir, paramSet,   " -maxdepth 1 -name \"", outputFilePrefix, "__simu_cancer_cell_", cellNumA, ".hg19*__median.bed\""), intern=T), stringsAsFactors=F, header=F)
medianDatB <- read.table(system(paste0("find ", dataDir, paramSet,   " -maxdepth 1 -name \"", outputFilePrefix, "__simu_cancer_cell_", cellNumB, ".hg19*__median.bed\""), intern=T), stringsAsFactors=F, header=F)

modeDatA <- read.table(system(paste0("find ", dataDir, paramSet,   " -maxdepth 1 -name \"", outputFilePrefix, "__simu_cancer_cell_", cellNumA, ".hg19*__mode.bed\""), intern=T), stringsAsFactors=F, header=F)
modeDatB <- read.table(system(paste0("find ", dataDir, paramSet,   " -maxdepth 1 -name \"", outputFilePrefix, "__simu_cancer_cell_", cellNumB, ".hg19*__mode.bed\""), intern=T), stringsAsFactors=F, header=F)

dipAvg <- read.table(system(paste0("find ", paramSet, " -name simu_healthy_avg.bed"), intern=T), stringsAsFactors=F, header=F)


colnames(simuDatA) <- colnames(simuDatB) <- c(header, "depth")
colnames(trueDatA) <- colnames(trueDatB) <- c(header, "true")
colnames(sconceDatA) <- colnames(sconceDatB) <- c(header, "copyNumber")
colnames(pairDatA) <- colnames(pairDatB) <- c(header, "copyNumber")
colnames(meanDatA) <- colnames(meanDatB) <- c(header, "copyNumber")
colnames(medianDatA) <- colnames(medianDatB) <- c(header, "copyNumber")
colnames(modeDatA) <- colnames(modeDatB) <- c(header, "copyNumber")
colnames(dipAvg) <- c(header, "mean", "var")

simuDatA$idx <- simuDatB$idx <- 1:nrow(simuDatA)
trueDatA$idx <- trueDatB$idx <- 1:nrow(trueDatA)
sconceDatA$idx <- sconceDatB$idx <- 1:nrow(sconceDatA)
pairDatA$idx <- pairDatB$idx <- 1:nrow(pairDatA)
meanDatA$idx <- meanDatB$idx <- 1:nrow(meanDatA)
medianDatA$idx <- medianDatB$idx <- 1:nrow(medianDatA)
modeDatA$idx <- modeDatB$idx <- 1:nrow(modeDatA)
dipAvg$idx <- 1:nrow(dipAvg)

sconceDatA$source <- sconceDatB$source <- "sconce"
pairDatA$source <- pairDatB$source <- "pair"
meanDatA$source <- meanDatB$source <- "mean"
medianDatA$source <- medianDatB$source <- "median"
modeDatA$source <- modeDatB$source <- "mode"

simuDatA$cell <- trueDatA$cell <- sconceDatA$cell <- pairDatA$cell <- meanDatA$cell <- medianDatA$cell <- modeDatA$cell <- cellIDs["cellNameA"]
simuDatB$cell <- trueDatB$cell <- sconceDatB$cell <- pairDatB$cell <- meanDatB$cell <- medianDatB$cell <- modeDatB$cell <- cellIDs["cellNameB"]

cnvScalingA <- sum(simuDatA$depth) / mean(trueDatA$true) / nrow(simuDatA)
cnvScalingB <- sum(simuDatB$depth) / mean(trueDatB$true) / nrow(simuDatB)
dipScaling <- mean(dipAvg$mean) / 2

inferredPloidiesA <- rbind(sconceDatA, pairDatA, meanDatA, medianDatA, modeDatA)
inferredPloidiesB <- rbind(sconceDatB, pairDatB, meanDatB, medianDatB, modeDatB)
inferredPloidies <- rbind(inferredPloidiesA, inferredPloidiesB)

simuDat <- rbind(simuDatA, simuDatB)
trueDat <- rbind(trueDatA, trueDatB)


sumSq <- do.call(rbind, lapply(cellIDs, FUN=function(currCell) {
  dat <- do.call(rbind, lapply(unique(inferredPloidies$source), FUN=function(currSource) {
    sse <- sum((subset(trueDat, cell == currCell)$true - subset(inferredPloidies, source == currSource & cell == currCell)$copyNumber)^2)
    data.frame(source=currSource, SSE=sprintf("SSE = %.3f", sse))
  }))
  dat$cell <- currCell
  dat
}))

orderedSources <- c("sconce", "pair", "mean", "median", "mode")
facetLabels <- c("SCONCE", "one pair", "mean", "median", "mode", "cell A", "cell B") # get all caps SCONCE and generic cell labels
names(facetLabels) <- c("sconce", "pair", "mean", "median", "mode", as.name(cellIDs["cellNameA"]), as.name(cellIDs["cellNameB"]))
inferredPloidies$source <- factor(inferredPloidies$source, levels=orderedSources)
sumSq$source <- factor(sumSq$source, levels=orderedSources)
orderedCells <- cellIDs
inferredPloidies$cell <- factor(inferredPloidies$cell, levels=orderedCells)
sumSq$cell <- factor(sumSq$cell, levels=orderedCells)
trueDat$cell <- factor(trueDat$cell, levels=orderedCells)
simuDatA$cell <- factor(simuDatA$cell, levels=orderedCells)
simuDatB$cell <- factor(simuDatB$cell, levels=orderedCells)

arrows$source <- factor(arrows$source, levels=orderedSources)
arrows$cell <- factor(arrows$cell, levels=orderedCells)

copyNumberAxisMax <- ceiling(max(max(subset(simuDatA, idx > minIdx & idx < maxIdx)$depth)/cnvScalingA, max(subset(simuDatB, idx > minIdx & idx < maxIdx)$depth)/cnvScalingB))

cellStripFills <- c("#1CADE4", "#FFC000")
cellTextCols <- c("gray10", "gray10")


changeSourceStripColors <- function(p, side, g=NULL) {
  if(is.null(g)) {
    g <- ggplot_gtable(ggplot_build(p))
  }
  stripr <- which(grepl(paste0("strip-", side), g$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- sourceStripFills[k]
    j <- which(grepl('title', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- sourceTextCols[k] # use str(g$grobs) to find elements to change
    k <- k+1
  }
  g
}

changeCellStripColors <- function(p, side, g=NULL) {
  if(is.null(g)) {
    g <- ggplot_gtable(ggplot_build(p))
  }
  stripr <- which(grepl(paste0("strip-", side), g$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- cellStripFills[k]
    j <- which(grepl('title', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$children[[1]]$gp$col <- cellTextCols[k] # use str(g$grobs) to find elements to change
    k <- k+1
  }
  g
}

# sconce | one pair | mean | median | mode
#                                          | cell A
#                                          | cell B
p <- ggplot(subset(inferredPloidies, idx > minIdx & idx < maxIdx), aes(x=idx, y=copyNumber, colour=cell)) +
  theme_bw() + guides(colour="none") +
  geom_line(data=subset(dipAvg, idx > minIdx & idx < maxIdx), aes(x=idx, y=mean / dipScaling), colour="lightblue", alpha=1) +
  geom_ribbon(data=subset(dipAvg, idx > minIdx & idx < maxIdx), aes(x=idx, ymin=mean/dipScaling - sqrt(var/(dipScaling^2)), ymax=mean/dipScaling + sqrt(var/(dipScaling^2))), fill="lightblue", alpha=0.3, inherit.aes=F) +
  geom_point(data=subset(simuDatA, idx > minIdx & idx < maxIdx), aes(x=idx, y=depth / cnvScalingA), alpha=0.5, size=0.85, colour="darkgray") +
  geom_point(data=subset(simuDatB, idx > minIdx & idx < maxIdx), aes(x=idx, y=depth / cnvScalingB), alpha=0.5, size=0.85, colour="darkgray") +
  geom_step(size=1.15) +
  facet_grid(cell ~ source, labeller=as_labeller(facetLabels)) +
  geom_step(data=subset(trueDat, idx > minIdx & idx < maxIdx), aes(x=idx, y=true), colour="red", linetype="dashed") +
  geom_segment(data=arrows, aes(x=x, xend=xend, y=y, yend=yend), colour="black", arrow=arrow(length=unit(0.5, "char"))) +
  scale_y_continuous(sec.axis=sec_axis(~ . * dipScaling, name="scaled read depth"), limits=c(0, copyNumberAxisMax), breaks=seq(0, copyNumberAxisMax, 2)) +
  labs(y="copy number", x="genomic index") +
  geom_text(data=sumSq, aes(x=-Inf, y=Inf, label=SSE, vjust=1.5, hjust=-0.1), colour="black") +
  scale_colour_manual(values=cellStripFills) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

plotWidth <- 8
plotHeight <- 4
outputFile <- paste0(plotRoot, "_all.pdf")
pdf(outputFile, width=plotWidth, height=plotHeight)
g <- changeCellStripColors(p, "r")
g <- changeSourceStripColors(p, "t",g)
grid::grid.draw(g)
dev.off()


#######################################################################
#######################################################################
#######################################################################
##### Real data
#######################################################################
#######################################################################
#######################################################################
orderedSources <- c("sconce", "pair", "mean", "median", "mode")
cellStripFills <- c("#1CADE4", "#FFC000")
cellTextCols <- c("gray10", "gray10")
sourceStripFills <- hue_pal()(length(orderedSources) + 1) # plus 1 for aneufinder
sourceTextCols <- rep("gray10", length(sourceStripFills))
header <- c("chr", "start", "end", "copyNumber")
k <- 10
filekeys <- c("sconce2")

############
# navin data
dataDir <- "Navin_Nature2011_hg19/output/"
currOutputDir <- paste0(dataDir, "/plots/")
paramSet <- "Navin_Nature2011_hg19"
cellRegex <- "SRR[0-9]*."
if(!dir.exists(currOutputDir)) {
  dir.create(currOutputDir)
}

# plotting "SRR054596.", "SRR054609.", 6500-6700
cellA <- "SRR054596."
cellB <- "SRR054609."
minIdx <- 6500
maxIdx <- 6700

cellIDs <- c(cellA, cellB)
names(cellIDs) <- cellIDs
orderedCells <- cellIDs

arrows <- data.frame(x=   c(6520,6520,6628,6628,6628,6628,    6520,6520,6628,6628,6628,6628),
                     xend=c(6520,6520,6628,6628,6628,6628,    6520,6520,6628,6628,6628,6628),
                     y=c(rep(9, 6), rep(8,6)),
                     yend=c(rep(7.5,6), rep(6.5,6)),
                     source=c("sconce", "pair", "pair", "mean", "median", "mode",       "sconce", "pair", "pair", "mean", "median", "mode"),
                     cell=c(rep(cellA,6), rep(cellB,6)))
arrows$source <- factor(arrows$source, levels=orderedSources)
arrows$cell <- factor(arrows$cell, levels=orderedCells)


plotRoot <- paste0(dataDir, "/plots/betterBounds_cells_", cellA, "_", cellB, "_", paramSet, "_", key, "_k", k, "_c", numCells, "_", minIdx, "-", maxIdx)

outputFilePrefix <- paste0("output_", key, "_", paramSet, "_k", k)
unfiltSconceFileList <- system(paste0("find ", dataDir, " -maxdepth 1 -name \"", outputFilePrefix, "*__sconce__*.bed\" | sort -V"), intern=T)
unfiltPairsFileList <- system(paste0("find ", dataDir,   " -maxdepth 1 -name \"", outputFilePrefix, "*__pair_*.bed\" | sort -V"), intern=T)
unfiltMeanFileList <- system(paste0("find ", dataDir,     " -maxdepth 1 -name \"", outputFilePrefix, "*__mean.bed\" | sort -V"), intern=T)
unfiltMedianFileList <- system(paste0("find ", dataDir, " -maxdepth 1 -name \"", outputFilePrefix, "*__median.bed\" | sort -V"), intern=T)
unfiltModeFileList <- system(paste0("find ", dataDir,     " -maxdepth 1 -name \"", outputFilePrefix, "*__mode.bed\" | sort -V"), intern=T)

sconceFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltSconceFileList[grepl(cellID, unfiltSconceFileList, fixed=T)]})
pairsFileList <- unfiltPairsFileList[sapply(lapply(str_extract_all(unfiltPairsFileList, cellRegex), unique), FUN=function(cellNames) {cellNames[1] %in% cellIDs && cellNames[2] %in% cellIDs})]
names(pairsFileList) <- sapply(str_extract_all(pairsFileList, cellRegex), FUN=function(cellNames){cellNames[3]})
meanFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltMeanFileList[grepl(cellID, unfiltMeanFileList, fixed=T)]})
medianFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltMedianFileList[grepl(cellID, unfiltMedianFileList, fixed=T)]})
modeFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltModeFileList[grepl(cellID, unfiltModeFileList, fixed=T)]})

readsDatA <- read.table(system(paste0("find ", dataDir, "/../ -name \"", cellA, "cov_unif_250kb\""), intern=T), stringsAsFactors=F, header=F)
sconceDatA <- read.table(sconceFileList[cellA])
pairDatA <- read.table(pairsFileList[cellA])
meanDatA <- read.table(meanFileList[cellA])
medianDatA <- read.table(medianFileList[cellA])
modeDatA <- read.table(modeFileList[cellA])

readsDatB <- read.table(system(paste0("find ", dataDir, "/../ -name \"", cellB, "cov_unif_250kb\""), intern=T), stringsAsFactors=F, header=F)
sconceDatB <- read.table(sconceFileList[cellB])
pairDatB <- read.table(pairsFileList[cellB])
meanDatB <- read.table(meanFileList[cellB])
medianDatB <- read.table(medianFileList[cellB])
modeDatB <- read.table(modeFileList[cellB])

dipAvg <- read.table(system(paste0("find ", dataDir, "/../diploid -name diploid_avg_cov_unif_250kb.bed"), intern=T), stringsAsFactors=F, header=F)

header <- header[1:3]
colnames(readsDatA) <- colnames(readsDatB) <- c(header, "depth")
colnames(sconceDatA) <- colnames(sconceDatB) <- c(header, "copyNumber")
colnames(pairDatA) <- colnames(pairDatB) <- c(header, "copyNumber")
colnames(meanDatA) <- colnames(meanDatB) <- c(header, "copyNumber")
colnames(medianDatA) <- colnames(medianDatB) <- c(header, "copyNumber")
colnames(modeDatA) <- colnames(modeDatB) <- c(header, "copyNumber")
colnames(dipAvg) <- c(header, "mean", "var")

readsDatA$idx <- readsDatB$idx <- 1:nrow(readsDatA)
sconceDatA$idx <- sconceDatB$idx <- 1:nrow(sconceDatA)
pairDatA$idx <- pairDatB$idx <- 1:nrow(pairDatA)
meanDatA$idx <- meanDatB$idx <- 1:nrow(meanDatA)
medianDatA$idx <- medianDatB$idx <- 1:nrow(medianDatA)
modeDatA$idx <- modeDatB$idx <- 1:nrow(modeDatA)
dipAvg$idx <- 1:nrow(dipAvg)

sconceDatA$source <- sconceDatB$source <- "sconce"
pairDatA$source <- pairDatB$source <- "pair"
meanDatA$source <- meanDatB$source <- "mean"
medianDatA$source <- medianDatB$source <- "median"
modeDatA$source <- modeDatB$source <- "mode"

readsDatA$cell <- sconceDatA$cell <- pairDatA$cell <- meanDatA$cell <- medianDatA$cell <- modeDatA$cell <- factor(cellA, levels=orderedCells)
readsDatB$cell <- sconceDatB$cell <- pairDatB$cell <- meanDatB$cell <- medianDatB$cell <- modeDatB$cell <- factor(cellB, levels=orderedCells)

inferredPloidiesA <- rbind(sconceDatA, pairDatA, meanDatA, medianDatA, modeDatA)
inferredPloidiesB <- rbind(sconceDatB, pairDatB, meanDatB, medianDatB, modeDatB)
inferredPloidies <- rbind(inferredPloidiesA, inferredPloidiesB)

sconceDat <- rbind(sconceDatA, sconceDatB)
sconceDat$source <- NULL

facetLabels <- c("SCONCE", "one pair", "mean", "median", "mode", "cell A", "cell B") # get all caps SCONCE and generic cell labels
names(facetLabels) <- c("sconce", "pair", "mean", "median", "mode", as.name(cellA), as.name(cellB))
inferredPloidies$source <- factor(inferredPloidies$source, levels=orderedSources)
inferredPloidies$cell <- factor(inferredPloidies$cell, levels=orderedCells)

scalingFactorsA <-  sapply(orderedSources, FUN=function(x) { sum(readsDatA$depth) / mean(inferredPloidies[inferredPloidies$source == x & inferredPloidies$cell == cellA, "copyNumber"], na.rm=T) / nrow(readsDatA)})
scalingFactorsB <-  sapply(orderedSources, FUN=function(x) { sum(readsDatB$depth) / mean(inferredPloidies[inferredPloidies$source == x & inferredPloidies$cell == cellB, "copyNumber"], na.rm=T) / nrow(readsDatB)})
cnvScalingA <- mean(scalingFactorsA)
cnvScalingB <- mean(scalingFactorsB)
dipScaling <- mean(dipAvg$mean) / 2

copyNumberAxisMax <- ceiling(max(max(subset(readsDatA, idx > minIdx & idx < maxIdx)$depth)/cnvScalingA, max(subset(readsDatB, idx > minIdx & idx < maxIdx)$depth)/cnvScalingB))

p <- ggplot(subset(inferredPloidies, idx > minIdx & idx < maxIdx), aes(x=idx, y=copyNumber, colour=cell)) +
  theme_bw() + guides(colour="none") +
  geom_line(data=subset(dipAvg, idx > minIdx & idx < maxIdx), aes(x=idx, y=mean / dipScaling), colour="lightblue", alpha=1) +
  geom_ribbon(data=subset(dipAvg, idx > minIdx & idx < maxIdx), aes(x=idx, ymin=mean/dipScaling - sqrt(var/(dipScaling^2)), ymax=mean/dipScaling + sqrt(var/(dipScaling^2))), fill="lightblue", alpha=0.3, inherit.aes=F) +
  geom_point(data=subset(readsDatA, idx > minIdx & idx < maxIdx), aes(x=idx, y=depth / cnvScalingA), alpha=0.5, size=0.85, colour="darkgray") +
  geom_point(data=subset(readsDatB, idx > minIdx & idx < maxIdx), aes(x=idx, y=depth / cnvScalingB), alpha=0.5, size=0.85, colour="darkgray") +
  geom_step(size=1.15) +
  facet_grid(cell ~ source, labeller=as_labeller(facetLabels)) +
  geom_segment(data=arrows, aes(x=x, xend=xend, y=y, yend=yend), colour="black", arrow=arrow(length=unit(0.5, "char"))) +
  scale_y_continuous(sec.axis=sec_axis(~ . * dipScaling, name="scaled read depth"), limits=c(0, copyNumberAxisMax), breaks=seq(0, copyNumberAxisMax, 2)) +
  labs(y="copy number", x="genomic index") +
  scale_colour_manual(values=cellStripFills) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

plotWidth <- 8
plotHeight <- 4
outputFile <- paste0(plotRoot, "_all.pdf")
pdf(outputFile, width=plotWidth, height=plotHeight)
g <- changeCellStripColors(p, "r")
g <- changeSourceStripColors(p, "t",g)
grid::grid.draw(g)
dev.off()
##########
# end navin data



############
# 10x data
dataDir <- "10x/output"
currOutputDir <- paste0(dataDir, "/plots/")
paramSet <- "10x_breast_tissue"
cellRegex <- "cell_[0-9]*_[ACGT]*-1."
if(!dir.exists(currOutputDir)) {
  dir.create(currOutputDir)
}

# plotting "cell_1538_TAGCCGGCAAGAACTA-1.", "cell_1360_GCATACAAGTAACCCT-1.", 3900-4100
cellA <- "cell_1538_TAGCCGGCAAGAACTA-1."
cellB <- "cell_1360_GCATACAAGTAACCCT-1."
minIdx <- 3900
maxIdx <- 4120

cellIDs <- c(cellA, cellB)
names(cellIDs) <- cellIDs
orderedCells <- cellIDs

arrows <- data.frame(x=   c(4024,4024, 3932, 4024,4024,4024,    4024,4024, 3932, 4024,4024,4024),
                     xend=c(4024,4024, 3932, 4024,4024,4024,    4024,4024, 3932, 4024,4024,4024),
                     y=c(rep(8, 6), rep(8.5,6)),
                     yend=c(rep(6.5,6), rep(7,6)),
                     source=c("sconce", "pair", "mean", "mean", "median", "mode",       "sconce", "pair", "mean", "mean", "median", "mode"),
                     cell=c(rep(cellA,6), rep(cellB,6)))
arrows$source <- factor(arrows$source, levels=orderedSources)
arrows$cell <- factor(arrows$cell, levels=orderedCells)


plotRoot <- paste0(dataDir, "/plots/betterBounds_cells_", cellA, "_", cellB, "_", paramSet, "_", key, "_k", k, "_c", numCells, "_", minIdx, "-", maxIdx)

outputFilePrefix <- paste0("output_", key, "_", paramSet, "_k", k)
unfiltSconceFileList <- system(paste0("find ", dataDir, " -maxdepth 1 -name \"", outputFilePrefix, "*__sconce__*.bed\" | sort -V"), intern=T)
unfiltPairsFileList <- system(paste0("find ", dataDir,   " -maxdepth 1 -name \"", outputFilePrefix, "*__pair_*.bed\" | sort -V"), intern=T)
unfiltMeanFileList <- system(paste0("find ", dataDir,     " -maxdepth 1 -name \"", outputFilePrefix, "*__mean.bed\" | sort -V"), intern=T)
unfiltMedianFileList <- system(paste0("find ", dataDir, " -maxdepth 1 -name \"", outputFilePrefix, "*__median.bed\" | sort -V"), intern=T)
unfiltModeFileList <- system(paste0("find ", dataDir,     " -maxdepth 1 -name \"", outputFilePrefix, "*__mode.bed\" | sort -V"), intern=T)

sconceFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltSconceFileList[grepl(cellID, unfiltSconceFileList, fixed=T)]})
pairsFileList <- unfiltPairsFileList[sapply(lapply(str_extract_all(unfiltPairsFileList, cellRegex), unique), FUN=function(cellNames) {cellNames[1] %in% cellIDs && cellNames[2] %in% cellIDs})]
names(pairsFileList) <- sapply(str_extract_all(pairsFileList, cellRegex), FUN=function(cellNames){cellNames[3]})
meanFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltMeanFileList[grepl(cellID, unfiltMeanFileList, fixed=T)]})
medianFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltMedianFileList[grepl(cellID, unfiltMedianFileList, fixed=T)]})
modeFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltModeFileList[grepl(cellID, unfiltModeFileList, fixed=T)]})

readsDatA <- read.table(system(paste0("find ", dataDir, "/../ -name \"", cellA, "q20.cov_unif_250kb\""), intern=T), stringsAsFactors=F, header=F)
sconceDatA <- read.table(sconceFileList[cellA])
pairDatA <- read.table(pairsFileList[cellA])
meanDatA <- read.table(meanFileList[cellA])
medianDatA <- read.table(medianFileList[cellA])
modeDatA <- read.table(modeFileList[cellA])

readsDatB <- read.table(system(paste0("find ", dataDir, "/../ -name \"", cellB, "q20.cov_unif_250kb\""), intern=T), stringsAsFactors=F, header=F)
sconceDatB <- read.table(sconceFileList[cellB])
pairDatB <- read.table(pairsFileList[cellB])
meanDatB <- read.table(meanFileList[cellB])
medianDatB <- read.table(medianFileList[cellB])
modeDatB <- read.table(modeFileList[cellB])

dipAvg <- read.table(system(paste0("find ", dataDir, "/../breast_tissue_A_2k -name breast_tissue_A_2k_avg_cov_unif_250kb.bed"), intern=T), stringsAsFactors=F, header=F)
dipAvg <- dipAvg[-as.numeric(rownames(tail(dipAvg[order(dipAvg[,4]),],5))),] # remove top 5 windows with biggest dip counts, ruins plotting

header <- header[1:3]
colnames(readsDatA) <- colnames(readsDatB) <- c(header, "depth")
colnames(sconceDatA) <- colnames(sconceDatB) <- c(header, "copyNumber")
colnames(pairDatA) <- colnames(pairDatB) <- c(header, "copyNumber")
colnames(meanDatA) <- colnames(meanDatB) <- c(header, "copyNumber")
colnames(medianDatA) <- colnames(medianDatB) <- c(header, "copyNumber")
colnames(modeDatA) <- colnames(modeDatB) <- c(header, "copyNumber")
colnames(dipAvg) <- c(header, "mean", "var")

readsDatA$idx <- readsDatB$idx <- 1:nrow(readsDatA)
sconceDatA$idx <- sconceDatB$idx <- 1:nrow(sconceDatA)
pairDatA$idx <- pairDatB$idx <- 1:nrow(pairDatA)
meanDatA$idx <- meanDatB$idx <- 1:nrow(meanDatA)
medianDatA$idx <- medianDatB$idx <- 1:nrow(medianDatA)
modeDatA$idx <- modeDatB$idx <- 1:nrow(modeDatA)
dipAvg$idx <- 1:nrow(dipAvg)

sconceDatA$source <- sconceDatB$source <- "sconce"
pairDatA$source <- pairDatB$source <- "pair"
meanDatA$source <- meanDatB$source <- "mean"
medianDatA$source <- medianDatB$source <- "median"
modeDatA$source <- modeDatB$source <- "mode"

readsDatA$cell <- sconceDatA$cell <- pairDatA$cell <- meanDatA$cell <- medianDatA$cell <- modeDatA$cell <- factor(cellA, levels=orderedCells)
readsDatB$cell <- sconceDatB$cell <- pairDatB$cell <- meanDatB$cell <- medianDatB$cell <- modeDatB$cell <- factor(cellB, levels=orderedCells)

inferredPloidiesA <- rbind(sconceDatA, pairDatA, meanDatA, medianDatA, modeDatA)
inferredPloidiesB <- rbind(sconceDatB, pairDatB, meanDatB, medianDatB, modeDatB)
inferredPloidies <- rbind(inferredPloidiesA, inferredPloidiesB)

sconceDat <- rbind(sconceDatA, sconceDatB)
sconceDat$source <- NULL

facetLabels <- c("SCONCE", "one pair", "mean", "median", "mode", "cell A", "cell B") # get all caps SCONCE and generic cell labels
names(facetLabels) <- c("sconce", "pair", "mean", "median", "mode", cellA, cellB)
inferredPloidies$source <- factor(inferredPloidies$source, levels=orderedSources)
inferredPloidies$cell <- factor(inferredPloidies$cell, levels=orderedCells)

scalingFactorsA <-  sapply(orderedSources, FUN=function(x) { sum(readsDatA$depth) / mean(inferredPloidies[inferredPloidies$source == x & inferredPloidies$cell == cellA, "copyNumber"], na.rm=T) / nrow(readsDatA)})
scalingFactorsB <-  sapply(orderedSources, FUN=function(x) { sum(readsDatB$depth) / mean(inferredPloidies[inferredPloidies$source == x & inferredPloidies$cell == cellB, "copyNumber"], na.rm=T) / nrow(readsDatB)})
cnvScalingA <- mean(scalingFactorsA)
cnvScalingB <- mean(scalingFactorsB)
dipScaling <- mean(dipAvg$mean) / 2

copyNumberAxisMax <- ceiling(max(max(subset(readsDatA, idx > minIdx & idx < maxIdx)$depth)/cnvScalingA, max(subset(readsDatB, idx > minIdx & idx < maxIdx)$depth)/cnvScalingB))

p <- ggplot(subset(inferredPloidies, idx > minIdx & idx < maxIdx), aes(x=idx, y=copyNumber, colour=cell)) +
  theme_bw() + guides(colour="none") +
  geom_line(data=subset(dipAvg, idx > minIdx & idx < maxIdx), aes(x=idx, y=mean / dipScaling), colour="lightblue", alpha=1) +
  geom_ribbon(data=subset(dipAvg, idx > minIdx & idx < maxIdx), aes(x=idx, ymin=mean/dipScaling - sqrt(var/(dipScaling^2)), ymax=mean/dipScaling + sqrt(var/(dipScaling^2))), fill="lightblue", alpha=0.3, inherit.aes=F) +
  geom_point(data=subset(readsDatA, idx > minIdx & idx < maxIdx), aes(x=idx, y=depth / cnvScalingA), alpha=0.5, size=0.85, colour="darkgray") +
  geom_point(data=subset(readsDatB, idx > minIdx & idx < maxIdx), aes(x=idx, y=depth / cnvScalingB), alpha=0.5, size=0.85, colour="darkgray") +
  geom_step(size=1.15) +
  facet_grid(cell ~ source, labeller=as_labeller(facetLabels)) +
  geom_segment(data=arrows, aes(x=x, xend=xend, y=y, yend=yend), colour="black", arrow=arrow(length=unit(0.5, "char"))) +
  scale_y_continuous(sec.axis=sec_axis(~ . * dipScaling, name="scaled read depth"), limits=c(0, copyNumberAxisMax), breaks=seq(0, copyNumberAxisMax, 2)) +
  labs(y="copy number", x="genomic index") +
  scale_colour_manual(values=cellStripFills) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

plotWidth <- 8
plotHeight <- 4
outputFile <- paste0(plotRoot, "_all.pdf")
pdf(outputFile, width=plotWidth, height=plotHeight)
g <- changeCellStripColors(p, "r")
g <- changeSourceStripColors(p, "t",g)
grid::grid.draw(g)
dev.off()
###########
## end 10x data

