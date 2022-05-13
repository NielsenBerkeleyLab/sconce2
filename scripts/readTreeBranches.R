# script of shared tree reading functions and constants (newick strings, branch lengths)

###########################################
# Tree manipulations for trees read from newick strings
###########################################
# node 127 == label 1 == cellnum 127
# node 128 == label 2 == cellnum 126
# node 254 == label 128 == cellnum 0
convertCellNumToNewickLabel <- function(cellFilename) { # ex simu_cancer_cell_0.hg19_lite.depth
  cellNum <- as.numeric(str_extract(str_extract(cellFilename, "cell_[0-9]+"), "[0-9]+"))
  label <- 128 - cellNum
  label
}

# given mrcaMat, distMat, nodeDepths, and branchLength, calculates all possible tree branch combinations
# left | right | t1 | t2 | t3
calcTreeBranches <- function(mrcaMat, distMat, nodeDepths, branchLength) {
  do.call(rbind, lapply(1:(ncol(mrcaMat)), FUN=function(left) {
    do.call(rbind, lapply(1:ncol(mrcaMat), FUN=function(right) {
      currMrca <- mrcaMat[left, right]
      t1 <- nodeDepths[currMrca] + branchLength # add for time before first split
      t2 <- distMat[left, currMrca]
      t3 <- distMat[right, currMrca]
      results <- data.frame(left=left, right=right, t1=t1, t2=t2, t3=t3, t2_t3=t2+t3)
      results
    }))
  }))
}

###########################################
# File reading functions
###########################################
# reads .hmm files to get tree branch estimates and generates true tree branches; returns a matrix with
# left (newick label) | right | variable (branch number) | cell0 (filename) | cell1 | inferred | true
getTreeBranches <- function(newickString, branchLength, hmmFile, forceRecalc=F) {
  outputFile <- paste0(hmmFile, ".treeBranches.txt")
  if(!forceRecalc && file.exists(outputFile)) {
    merged <- read.table(outputFile, sep="\t", header=T, stringsAsFactors=F)
    return(merged)
  }
  treeObj <- read.tree(text=newickString)
  mrcaMat <- mrca(treeObj)
  distMat <- dist.nodes(treeObj)
  nodeDepths <- node.depth.edgelength(treeObj)

  # left | right | t1 | t2 | t3
  # left/right correspond to node labels in the newick file
  treeBranches <- calcTreeBranches(mrcaMat, distMat, nodeDepths, branchLength)
  treeBranches_m <- melt(treeBranches, id.vars=c("left", "right"))

  # cell0 | cell1 | t1
  # cell0 | cell1 | t2
  # cell0 | cell1 | t3
  # ...
  # need to separate cell name reading from param reading in case doing nearest n
  hmmLines <- system(paste0("/bin/bash -c ", shQuote(sprintf("grep \"^HMM\" %s | grep \",\" | sed -e 's/HMM\\s\\+//' -e 's/(//' -e 's/)//' -e 's/,/\t/' -e 's/://' -e 's/\\s\\+/\t/'", hmmFile))), intern=T)
  if(length(hmmLines) == 0) {
    return(NULL)
  }
  analyzedHmmNames <- read.table(text=hmmLines)
  colnames(analyzedHmmNames) <- c("hmmIdx", "cell0", "cell1")
  paramIdxToExtract <- do.call(c, lapply(analyzedHmmNames$hmmIdx, FUN=function(i) {c(3*i, (3*i)+1, (3*i)+2)})) # don't extract hmms we skipped
  paramIdxToExtract <- paramIdxToExtract + 1 # paramsToEst is 0 indexed, R vector is 1 indexed
  unlabeledHmmParamsRaw <- read.table(text=system(paste0("/bin/bash -c ", shQuote(sprintf("sed '1,/FINAL HMM/d' %s | tac | sed '/paramsToEst/q' | tac | sed '/fixedParams/q' | tail -n +2 | head -n -2", hmmFile))), intern=T))
  hmmParamsRaw <- data.frame(cell0=do.call(c, lapply(analyzedHmmNames$cell0, rep, 3)), cell1=do.call(c, lapply(analyzedHmmNames$cell1, rep, 3)), value=unlabeledHmmParamsRaw[paramIdxToExtract,])
  
  hmmParamsRaw$variable <- c("t1", "t2", "t3")

  # cell0 | cell1 | t1 | t2 | t3
  # cell0/cell1 correspond to indv cell filenames
  hmmParams_wide <- dcast(hmmParamsRaw, cell0 + cell1 ~ variable)
  hmmParams_wide$t2_t3 <- hmmParams_wide$t2 + hmmParams_wide$t3
  hmmParams_wide_m <- melt(hmmParams_wide, id.vars=c("cell0", "cell1"))

  hmmParams_wide_m$left <- sapply(hmmParams_wide_m$cell0, convertCellNumToNewickLabel)
  hmmParams_wide_m$right <- sapply(hmmParams_wide_m$cell1, convertCellNumToNewickLabel)

  merged <- merge(hmmParams_wide_m, treeBranches_m, by=c("left", "right", "variable"))
  colnames(merged) <- c("left", "right", "variable", "cell0", "cell1", "inferred", "true")

  write.table(merged, file=outputFile, quote=F, sep="\t", col.names=T, row.names=F)

  merged
}

# returns matrix of true copy numbers from true_cancer_cell*, rows are cells, cols are genomic positions
readAllTumorTrueBeds <- function(dataDir, paramSet, numCells) {
  tumorListFile <- paste0(dataDir, paramSet, "/tumor_depths_", numCells)
  if(!file.exists(tumorListFile)) {
    next
  }
  tumorList <- read.table(tumorListFile)
  tumorTrueBeds <- sapply(tumorList, FUN=function(x) {gsub("depth", "bed", gsub("simu", "true", x))})
  
  tumorDepths <- lapply(tumorTrueBeds, FUN=function(x) {
    bed <- read.table(x, sep="\t", header=F, stringsAsFactors=F)
    colnames(bed) <- c("chr", "start", "end", "cn")
    bed$cellFile <- str_extract(x, "cancer_cell_[0-9]*.")
    bed$newickLab <- factor(128 - as.numeric(gsub(".$", "", gsub("cancer_cell_", "", bed$cellFile))))
    bed
  })
  
  newickLabs <- factor(128 - as.numeric(gsub(".$", "", gsub("cancer_cell_", "",str_extract(tumorList$V1, "cancer_cell_[0-9]*.")))))
  allDepths <- do.call(rbind, lapply(tumorDepths, FUN=function(x) {x$cn}))
  rownames(allDepths) <- newickLabs
  allDepths
}

# returns matrix of inferred copy numbers from outputType ("__sconce__*.bed", "*__mean.bed", "*__median.bed", "*__mode.bed"), rows are cells, cols are genomic positions
readAllOutputBedFiles <- function(dataDir, paramSet, key, numCells, outputType) {
  outputFilePrefix <- paste0("output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells) # based on scAllP_*sh outBase variable
  unfiltFileList <- system(paste0("find ", dataDir, paramSet, " -maxdepth 1 -name \"", outputFilePrefix, outputType, "\" | sort -V"), intern=T)
  tumorDepths <- read.table(paste0(dataDir, paramSet, "/tumor_depths_", numCells), stringsAsFactors=F)$V1
  cellIDs <- str_extract(tumorDepths, "cancer_cell_[0-9]*.")
  filtFileList <- sapply(cellIDs, FUN=function(cellID) {unfiltFileList[grepl(cellID, unfiltFileList, fixed=T)]})
  dat <- do.call(rbind, lapply(cellIDs, FUN=function(cellID) {
    bed <- read.table(filtFileList[grepl(cellID, filtFileList, fixed=T)], stringsAsFactors=F, sep="\t", header=F)
    bed$V4
  }))
  newickLabs <- factor(128 - as.numeric(gsub(".$", "", gsub("cancer_cell_", "",str_extract(tumorDepths, "cancer_cell_[0-9]*.")))))
  rownames(dat) <- newickLabs
  dat
}

readcnp2cnpFile <- function(cnp2cnpFile) {
  dat <- read.table(cnp2cnpFile, skip=1, row.names=1, stringsAsFactors=F)
  colnames(dat) <- rownames(dat) <- as.character(sapply(rownames(dat), convertCellNumToNewickLabel))
  as.matrix(dat)
}

###########################################
# Distance matrix functions
###########################################
# given mergedTreeBranches, creates a distance matrix based on t2+t3 values for uniqNodes
createInferred_t2_t3_distMat <- function(mergedTreeBranches, uniqNodes) {
  #inferredDistMat <- matrix(rep(0, length(uniqNodes)^2), ncol=length(uniqNodes)) # values from .hmm file
  inferredDistMat <- matrix(rep(NA, length(uniqNodes)^2), ncol=length(uniqNodes)) # values from .hmm file
  diag(inferredDistMat) <- 0
  dimnames(inferredDistMat) <- list(uniqNodes, uniqNodes)
  t2_t3 <- subset(mergedTreeBranches, variable == "t2_t3") # t2+t3 values from .hmm file
  for(i in 1:nrow(t2_t3)) {
    left <- as.character(t2_t3[i,"left"])
    right <- as.character(t2_t3[i,"right"])
    inferredDistMat[left, right] <- t2_t3[t2_t3$left == left & t2_t3$right == right, "inferred"]
    inferredDistMat[right, left] <- t2_t3[t2_t3$left == left & t2_t3$right == right, "inferred"]
  }
  inferredDistMat
}

# given mergedTreeBranches, creates a distance matrix based on true values (from the newick file) for uniqNodes
createTrueDistMat <- function(mergedTreeBranches, uniqNodes) {
  trueDistMat <- matrix(rep(0, length(uniqNodes)^2), ncol=length(uniqNodes)) # values from newick file
  dimnames(trueDistMat) <- list(uniqNodes, uniqNodes)
  t2_t3 <- subset(mergedTreeBranches, variable == "t2_t3") # t2+t3 values from .hmm file
  
  for(i in 1:nrow(t2_t3)) {
    left <- as.character(t2_t3[i,"left"])
    right <- as.character(t2_t3[i,"right"])
    trueDistMat[left, right] <- t2_t3[t2_t3$left == left & t2_t3$right == right, "true"]
    trueDistMat[right, left] <- t2_t3[t2_t3$left == left & t2_t3$right == right, "true"]
  }
  trueDistMat
}

# given a cnp2cnp output file name, reads it in, converts it to a dist mat
createCnpDistMat <- function(cnpFilename, cnpRevFilename) {
  if(!file.exists(cnpFilename)) {
    warning(paste0("could not find cnp2cnp files", cnpFilename, ", ", cnpRevFilename, ". Did you run scripts/cnaBedToFasta_wrapper.sh?"))
    return(NULL)
  }
  cnpMat <- readcnp2cnpFile(cnpFilename)
  cnpRevMat <- readcnp2cnpFile(cnpRevFilename)
  cnpMat <- cnpMat + cnpRevMat[rownames(cnpMat), colnames(cnpMat)]
  cnpMat
}

###########################################
# Neighbor joining and tree calculation functions
###########################################
createNJtrees <- function(mergedTreeBranches, paramSet, numCells, key) {
  uniqNodes <- as.character(sort(unique(c(mergedTreeBranches$left, mergedTreeBranches$right))))
  inferredDistMat <- createInferred_t2_t3_distMat(mergedTreeBranches, uniqNodes)
  trueDistMat <- createTrueDistMat(mergedTreeBranches, uniqNodes)

  # get cnp2cnp values
  trueCnpBase <- paste0(dataDir, "/", paramSet, "/tumor_depths_", numCells)
  cnpTrueMat <- createCnpDistMat(paste0(trueCnpBase, "_true.cnp2cnp"), paste0(trueCnpBase, "_trueRev.cnp2cnp"))

  cnpBase <- paste0(dataDir, "/", paramSet, "/output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells)
  cnpSconceMat <- createCnpDistMat(paste0(cnpBase, "_roundedsconce.cnp2cnp"), paste0(cnpBase, "_roundedsconceRev.cnp2cnp"))
  if(is.null(cnpSconceMat)) {
    return(NULL)
  }
  cnpMeanMat <- createCnpDistMat(paste0(cnpBase, "_roundedmean.cnp2cnp"), paste0(cnpBase, "_roundedmeanRev.cnp2cnp"))
  cnpMedianMat <- createCnpDistMat(paste0(cnpBase, "_roundedmedian.cnp2cnp"), paste0(cnpBase, "_roundedmedianRev.cnp2cnp"))
  cnpModeMat <- createCnpDistMat(paste0(cnpBase, "_roundedmode.cnp2cnp"), paste0(cnpBase, "_roundedmodeRev.cnp2cnp"))

  cnpTrueTree   <- nj(cnpTrueMat) # tree built on cnp2cnp metric btn true copy number profiles
  cnpSconceTree <- nj(cnpSconceMat) # tree built on cnp2cnp metric btn sconce copy number profiles
  cnpMeanTree   <- nj(cnpMeanMat) # tree built on cnp2cnp metric btn mean copy number profiles
  cnpMedianTree <- nj(cnpMedianMat) # tree built on cnp2cnp metric btn median copy number profiles
  cnpModeTree   <- nj(cnpModeMat) # tree built on cnp2cnp metric btn mode copy number profiles

  # zzs metric
  trueZzsBase <- paste0(dataDir, "/", paramSet, "/tumor_depths_", numCells)
  zzsTrueMat <- createCnpDistMat(paste0(trueZzsBase, "_true.zzs.cnp2cnp"), paste0(trueZzsBase, "_trueRev.zzs.cnp2cnp"))
  zzsBase <- paste0(dataDir, "/", paramSet, "/output_", key, "_", gsub("/", "_", paramSet), "_k", k, "_c", numCells)
  zzsSconceMat <- createCnpDistMat(paste0(zzsBase, "_roundedsconce.zzs.cnp2cnp"), paste0(zzsBase, "_roundedsconceRev.zzs.cnp2cnp"))
  zzsMeanMat <- createCnpDistMat(paste0(zzsBase, "_roundedmean.zzs.cnp2cnp"), paste0(zzsBase, "_roundedmeanRev.zzs.cnp2cnp"))
  zzsMedianMat <- createCnpDistMat(paste0(zzsBase, "_roundedmedian.zzs.cnp2cnp"), paste0(zzsBase, "_roundedmedianRev.zzs.cnp2cnp"))
  zzsModeMat <- createCnpDistMat(paste0(zzsBase, "_roundedmode.zzs.cnp2cnp"), paste0(zzsBase, "_roundedmodeRev.zzs.cnp2cnp"))

  zzsTrueTree   <- nj(zzsTrueMat) # tree built on zzs cnp2cnp metric btn true copy number profiles
  zzsSconceTree <- nj(zzsSconceMat) # tree built on zzs cnp2cnp metric btn sconce copy number profiles
  zzsMeanTree   <- nj(zzsMeanMat) # tree built on zzs cnp2cnp metric btn mean copy number profiles
  zzsMedianTree <- nj(zzsMedianMat) # tree built on zzs cnp2cnp metric btn median copy number profiles
  zzsModeTree   <- nj(zzsModeMat) # tree built on zzs cnp2cnp metric btn mode copy number profiles

  fullTrueNewickTree <- read.tree(text=newickStrings[paramSet]) # tree read directly from newick file
  trueDistTree <- keep.tip(fullTrueNewickTree, uniqNodes) # tree built directly from newick file, filtered for tips in this subset
  inferredDistTree <- njs(inferredDistMat) # tree built on inferred t2+t3 branch lengths, using njs for missing values if nearest10

  # manually calculate euclidean distance
  # get output from our bed files
  trueCNs <- readAllTumorTrueBeds(dataDir, paramSet, numCells) # values from true_cancer_cell* files
  sconceCNs <- readAllOutputBedFiles(dataDir, paramSet, key, numCells, "*__sconce__*.bed") # values from output*bed files
  meanCNs <- readAllOutputBedFiles(dataDir, paramSet, key, numCells, "*__mean.bed")
  medianCNs <- readAllOutputBedFiles(dataDir, paramSet, key, numCells, "*__median.bed")
  modeCNs <- readAllOutputBedFiles(dataDir, paramSet, key, numCells, "*__mode.bed")
  eucTrueTree   <- nj(dist(trueCNs, method="euclidean")) # tree built on euclidean distance btn true copy number profiles
  eucSconceTree <- nj(dist(sconceCNs, method="euclidean")) # tree built on euclidean distance btn sconce copy number profiles
  eucMeanTree   <- nj(dist(meanCNs, method="euclidean")) # tree built on euclidean distance btn mean copy number profiles
  eucMedianTree <- nj(dist(medianCNs, method="euclidean")) # tree built on euclidean distance btn median copy number profiles
  eucModeTree   <- nj(dist(modeCNs, method="euclidean")) # tree built on euclidean distance btn mode copy number profiles

  list(eucTrueTree=eucTrueTree, eucSconceTree=eucSconceTree, eucMeanTree=eucMeanTree, eucMedianTree=eucMedianTree, eucModeTree=eucModeTree, cnpTrueTree=cnpTrueTree, cnpSconceTree=cnpSconceTree, cnpMeanTree=cnpMeanTree, cnpMedianTree=cnpMedianTree, cnpModeTree=cnpModeTree, zzsTrueTree=zzsTrueTree, zzsSconceTree=zzsSconceTree, zzsMeanTree=zzsMeanTree, zzsMedianTree=zzsMedianTree, zzsModeTree=zzsModeTree, trueDistTree=trueDistTree, inferredDistTree=inferredDistTree)
}

calcRFdists <- function(treeList) {
  dists <- data.frame(eucTrue=  treedist(treeList[["trueDistTree"]], treeList[["eucTrueTree"]]),
                      eucSconce=treedist(treeList[["trueDistTree"]], treeList[["eucSconceTree"]]),
                      eucMean=  treedist(treeList[["trueDistTree"]], treeList[["eucMeanTree"]]),
                      eucMedian=treedist(treeList[["trueDistTree"]], treeList[["eucMedianTree"]]),
                      eucMode=  treedist(treeList[["trueDistTree"]], treeList[["eucModeTree"]]),
                      cnpTrue=  treedist(treeList[["trueDistTree"]], treeList[["cnpTrueTree"]]),
                      cnpSconce=treedist(treeList[["trueDistTree"]], treeList[["cnpSconceTree"]]),
                      cnpMean=  treedist(treeList[["trueDistTree"]], treeList[["cnpMeanTree"]]),
                      cnpMedian=treedist(treeList[["trueDistTree"]], treeList[["cnpMedianTree"]]),
                      cnpMode=  treedist(treeList[["trueDistTree"]], treeList[["cnpModeTree"]]),
                      zzsTrue=  treedist(treeList[["trueDistTree"]], treeList[["zzsTrueTree"]]),
                      zzsSconce=treedist(treeList[["trueDistTree"]], treeList[["zzsSconceTree"]]),
                      zzsMean=  treedist(treeList[["trueDistTree"]], treeList[["zzsMeanTree"]]),
                      zzsMedian=treedist(treeList[["trueDistTree"]], treeList[["zzsMedianTree"]]),
                      zzsMode=  treedist(treeList[["trueDistTree"]], treeList[["zzsModeTree"]]),
                      t2_t3=    treedist(treeList[["trueDistTree"]], treeList[["inferredDistTree"]]))
  dists$distMetric <- rownames(dists)
  dists_m <- melt(dists)
  dists_m$cleanTree <- factor(sapply(dists_m$variable, cleanTreeName), levels=orderedTreeNames)

  # reformat for better plotting labels
  dists_m$metric <- sapply(dists_m$cleanTree, FUN=function(x) {s <- unlist(str_split(x, " ")); paste0(s[1], "~", s[2])})
  dists_m$metric <- gsub("dist", "distance", dists_m$metric)
  dists_m$metric <- gsub("inferred~t2+t3", "t[2]+t[3]", dists_m$metric, fixed=T)
  dists_m$metric <- factor(dists_m$metric, levels=metricStrings)

  dists_m$bed <- sapply(dists_m$cleanTree, FUN=function(x) {s <- unlist(str_split(x, " "));s[4]})
  dists_m$bed[is.na(dists_m$bed)] <- "t[2]+t[3]"
  dists_m$bed <- factor(dists_m$bed, levels=bedStrings)

  dists_m
}

makeNJtreePlots <- function(treeList, outputFile) {
  njTreesPlotList <- list()
  for(treeName in names(treeList)) {
    title <- cleanTreeName(treeName) # TODO probably want to reorder and reduce, but not sure how yet
    p <- ggtree(treeList[[treeName]]) + geom_tiplab() + labs(title=title)
    njTreesPlotList[[treeName]] <- p
  }
  toSave <- plot_grid(plotlist=njTreesPlotList, align='vh', labels="AUTO")
  png(outputFile, width=9, height=6, units="in", res=600)
  plot(toSave); dev.off()
  njTreesPlotList
}

makeRFdistPlot <- function(dists, title) {
  metricColors <- hue_pal()(length(metricStrings)) # euc, cnp2cnp, medicc/zzs, t2+t3
  p <- ggplot(subset(dists, distMetric == "symmetric.difference"), aes(x=variable, colour=metric, y=value)) + geom_boxplot() + theme_bw() + theme(legend.title=element_blank(), axis.text.x=element_blank()) + scale_colour_manual(labels=scales::parse_format(), values=metricColors)
  if(!is.na(title)) {
    p <- p + labs(x="distance metric", y="RF distance", title=title)
  } else {
    p <- p + labs(x="distance metric", y="RF distance")
  }
  p
}

# faceted plot: sections/facets are distance metric, true/sconce/mean/median/mode[/t2+t3] are box plots
makeFacetedRFdistPlot <- function(dists, title) {
  bedColors <- hue_pal()(length(bedStrings)) # sconce, mean, median, mode, true, t2+t3
  p <- ggplot(subset(dists, distMetric == "symmetric.difference"), aes(x=bed, colour=bed, y=value)) + geom_boxplot() + theme_bw() + theme(legend.title=element_blank(), axis.text.x=element_blank()) + facet_grid(~metric, scales="free_x", space="free", labeller=label_parsed) + scale_colour_manual(labels=scales::parse_format(), values=bedColors)

  if(!is.na(title)) {
    p <- p  + labs(x="distance metric", y="RF distance", title=title)
  } else {
    p <- p  + labs(x="distance metric", y="RF distance")
  }
  p
}

# write median RF distances to text and tex files for easy copy/paste into latex tables
writeMedianRFdistTexFiles <- function(dists, medianRFdistFile) {
  programs <- levels(dists$variable)
  currParamSets <- unique(dists$paramSet)

  medianDists <- as.data.frame(t(sapply(currParamSets, FUN=function(set) {
    sapply(programs, FUN=function(program) {
      median(subset(dists, paramSet == set & variable == program & distMetric == "symmetric.difference")$value)
    })
  })))
  medianDists$paramSet <- rownames(medianDists)
  medianDists_m <- melt(medianDists)
  medianDists_m$cleanTree <- sapply(medianDists_m$variable, cleanTreeName)

  medianDists_m$metric <- sapply(medianDists_m$cleanTree, FUN=function(x) {s <- unlist(str_split(x, " ")); paste0(s[1], " ", s[2])})
  medianDists_m$metric <- gsub("dist", "distance", medianDists_m$metric)
  medianDists_m$metric <- gsub("inferred ", "", medianDists_m$metric)
  medianDists_m$metric <- factor(medianDists_m$metric, levels=c("Euclidean distance", "cnp2cnp distance", "MEDICC distance", "t2+t3"))

  medianDists_m$bed <- sapply(medianDists_m$cleanTree, FUN=function(x) {s <- unlist(str_split(x, " "));s[4]})
  medianDists_m$bed[is.na(medianDists_m$bed)] <- "t2+t3"
  medianDists_m$bed <- factor(medianDists_m$bed, levels=c("SCONCE", "mean", "median", "mode", "true", "t2+t3"))

  medianDists_w <- dcast(medianDists_m[,c("paramSet", "value", "metric", "bed")],paramSet + metric ~ bed)
  treeNames <- LETTERS[1:length(currParamSets)]
  names(treeNames) <- currParamSets
  medianDists_w$treeName <- treeNames[medianDists_w$paramSet]
  medianDists_w <- medianDists_w[,c("treeName", "metric", "SCONCE", "mean", "median", "mode", "true", "t2+t3", "paramSet")]

  write.table(medianDists_w, file=paste0(medianRFdistFile, "_medianDists.txt"), sep="\t", quote=F, row.names=T, col.names=NA)

  texColnames <- paste0("\\textbf{", colnames(medianDists_w), "}")
  names(texColnames) <- colnames(medianDists_w)

  medianDists_w_rounded <- as.data.frame(do.call(cbind, lapply(colnames(medianDists_w), FUN=function(col) {
    if(is.numeric(medianDists_w[,col])) {
      round(medianDists_w[,col], digits=4)
    } else {as.character(medianDists_w[,col])}
  })))
  colnames(medianDists_w_rounded) <- colnames(medianDists_w)

  medianDistsTex <- rbind(texColnames, medianDists_w_rounded)
  medianDistsTex$paramSet <- paste0(" \\\\ \\hline % ", medianDistsTex$paramSet)
  write.table(paste0(apply(medianDistsTex[,-ncol(medianDistsTex)], 1, FUN=function(x) {paste0(x, collapse=" & ")}),  medianDistsTex[,ncol(medianDistsTex)]), file=paste0(medianRFdistFile, "_medianDists.tex"), col.names=F, row.names=F, quote=F)
}


cleanTreeName <- function(f) {
  if(grepl("eucTrue", f)) {
    return("Euclidean dist on true CNPs")
  } else if(grepl("t2_t3", f) || grepl("inferred", f)) {
    return("inferred t2+t3")
  } else if(grepl("eucSconce", f)) {
    return("Euclidean dist on SCONCE CNPs")
  } else if(grepl("eucMean", f)) {
    return("Euclidean dist on mean CNPs")
  } else if(grepl("eucMedian", f)) {
    return("Euclidean dist on median CNPs")
  } else if(grepl("eucMode", f)) {
    return("Euclidean dist on mode CNPs")
  } else if(grepl("cnpTrue", f)) {
    return("cnp2cnp dist on true CNPs")
  } else if(grepl("cnpSconce", f)) {
    return("cnp2cnp dist on SCONCE CNPs")
  } else if(grepl("cnpMean", f)) {
    return("cnp2cnp dist on mean CNPs")
  } else if(grepl("cnpMedian", f)) {
    return("cnp2cnp dist on median CNPs")
  } else if(grepl("cnpMode", f)) {
    return("cnp2cnp dist on mode CNPs")
  } else if(grepl("zzsTrue", f)) {
    return("MEDICC dist on true CNPs")
  } else if(grepl("zzsSconce", f)) {
    return("MEDICC dist on SCONCE CNPs")
  } else if(grepl("zzsMean", f)) {
    return("MEDICC dist on mean CNPs")
  } else if(grepl("zzsMedian", f)) {
    return("MEDICC dist on median CNPs")
  } else if(grepl("zzsMode", f)) {
    return("MEDICC dist on mode CNPs")
  }
  return(f)
}

orderedTreeNames <- c("Euclidean dist on true CNPs",   
                      "Euclidean dist on SCONCE CNPs",
                      "Euclidean dist on mean CNPs",   
                      "Euclidean dist on median CNPs",
                      "Euclidean dist on mode CNPs",
                      "cnp2cnp dist on true CNPs",
                      "cnp2cnp dist on SCONCE CNPs",   
                      "cnp2cnp dist on mean CNPs",
                      "cnp2cnp dist on median CNPs",   
                      "cnp2cnp dist on mode CNPs",   
                      "MEDICC dist on true CNPs",
                      "MEDICC dist on SCONCE CNPs",   
                      "MEDICC dist on mean CNPs",
                      "MEDICC dist on median CNPs",   
                      "MEDICC dist on mode CNPs",   
                      "inferred t2+t3")

metricStrings <- c("Euclidean~distance", "cnp2cnp~distance", "MEDICC~distance", "t[2]+t[3]")
bedStrings <- c("SCONCE", "mean", "median", "mode", "true", "t[2]+t[3]")


###########################################
# Tree branch correlation plotting functions
###########################################
# expects results from getTreeBranches
makeCorrPlot <- function(merged_df, title) {
  # left | right | variable | cell0 | cell1 | inferred | true
  # based on https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
  lm_eqn <- function(df){
    m <- lm(inferred ~ true, df)
    rSq <- substitute(italic(r)^2~"="~r2, list(r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(rSq))
  }

  merged_df$variable <- factor(merged_df$variable, levels=c("t1", "t2", "t3", "t2_t3"))
  merged_df$labels <- factor(merged_df$variable)
  levels(merged_df$labels) <- c("t[1]", "t[2]", "t[3]", "t[2]+t[3]")

  lmStr <- ddply(merged_df,.(labels),lm_eqn)

  p <- ggplot(merged_df, aes(x=true, y=inferred, colour=labels)) + geom_point(alpha=0.5) + facet_grid(.~labels, labeller=label_parsed, scales="free_x") + scale_color_discrete(breaks=levels(merged_df$labels), labels=parse(text=levels(merged_df$labels))) + geom_smooth(method="lm", colour="gray35", show.legend=F) + geom_text(data=lmStr, aes(x=-Inf, y=Inf, label=V1, hjust=0, vjust=1), parse=TRUE, colour="black", show.legend=F) + theme_bw() + theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5)) + guides(colour=guide_legend(override.aes=list(alpha=1)))
  if(!is.na(title)) {
    p <- p + labs(colour="branch", x="true", y="inferred", title=title) 
  } else {
    p <- p + labs(colour="branch", x="true", y="inferred")
  }
  p
}

###########################################
# newick strings for simulations
###########################################
paramSets <- c("paramsA", "paramsB", "paramsC", "paramsD")
pAstr <- "(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((1:0.0078125,2:0.0078125):0.0078125,3:0.015625):0.0078125,4:0.0234375):0.0078125,5:0.03125):0.0078125,6:0.0390625):0.0078125,7:0.046875):0.0078125,8:0.0546875):0.0078125,9:0.0625):0.0078125,10:0.0703125):0.0078125,11:0.078125):0.0078125,12:0.0859375):0.0078125,13:0.09375):0.0078125,14:0.1015625):0.0078125,15:0.109375):0.0078125,16:0.1171875):0.0078125,17:0.125):0.0078125,18:0.1328125):0.0078125,19:0.140625):0.0078125,20:0.1484375):0.0078125,21:0.15625):0.0078125,22:0.1640625):0.0078125,23:0.171875):0.0078125,24:0.1796875):0.0078125,25:0.1875):0.0078125,26:0.1953125):0.0078125,27:0.203125):0.0078125,28:0.2109375):0.0078125,29:0.21875):0.0078125,30:0.2265625):0.0078125,31:0.234375):0.0078125,32:0.2421875):0.0078125,33:0.25):0.0078125,34:0.2578125):0.0078125,35:0.265625):0.0078125,36:0.2734375):0.0078125,37:0.28125):0.0078125,38:0.2890625):0.0078125,39:0.296875):0.0078125,40:0.3046875):0.0078125,41:0.3125):0.0078125,42:0.3203125):0.0078125,43:0.328125):0.0078125,44:0.3359375):0.0078125,45:0.34375):0.0078125,46:0.3515625):0.0078125,47:0.359375):0.0078125,48:0.3671875):0.0078125,49:0.375):0.0078125,50:0.3828125):0.0078125,51:0.390625):0.0078125,52:0.3984375):0.0078125,53:0.40625):0.0078125,54:0.4140625):0.0078125,55:0.421875):0.0078125,56:0.4296875):0.0078125,57:0.4375):0.0078125,58:0.4453125):0.0078125,59:0.453125):0.0078125,60:0.4609375):0.0078125,61:0.46875):0.0078125,62:0.4765625):0.0078125,63:0.484375):0.0078125,64:0.4921875):0.0078125,65:0.5):0.0078125,66:0.5078125):0.0078125,67:0.515625):0.0078125,68:0.5234375):0.0078125,69:0.53125):0.0078125,70:0.5390625):0.0078125,71:0.546875):0.0078125,72:0.5546875):0.0078125,73:0.5625):0.0078125,74:0.5703125):0.0078125,75:0.578125):0.0078125,76:0.5859375):0.0078125,77:0.59375):0.0078125,78:0.6015625):0.0078125,79:0.609375):0.0078125,80:0.6171875):0.0078125,81:0.625):0.0078125,82:0.6328125):0.0078125,83:0.640625):0.0078125,84:0.6484375):0.0078125,85:0.65625):0.0078125,86:0.6640625):0.0078125,87:0.671875):0.0078125,88:0.6796875):0.0078125,89:0.6875):0.0078125,90:0.6953125):0.0078125,91:0.703125):0.0078125,92:0.7109375):0.0078125,93:0.71875):0.0078125,94:0.7265625):0.0078125,95:0.734375):0.0078125,96:0.7421875):0.0078125,97:0.75):0.0078125,98:0.7578125):0.0078125,99:0.765625):0.0078125,100:0.7734375):0.0078125,101:0.78125):0.0078125,102:0.7890625):0.0078125,103:0.796875):0.0078125,104:0.8046875):0.0078125,105:0.8125):0.0078125,106:0.8203125):0.0078125,107:0.828125):0.0078125,108:0.8359375):0.0078125,109:0.84375):0.0078125,110:0.8515625):0.0078125,111:0.859375):0.0078125,112:0.8671875):0.0078125,113:0.875):0.0078125,114:0.8828125):0.0078125,115:0.890625):0.0078125,116:0.8984375):0.0078125,117:0.90625):0.0078125,118:0.9140625):0.0078125,119:0.921875):0.0078125,120:0.9296875):0.0078125,121:0.9375):0.0078125,122:0.9453125):0.0078125,123:0.953125):0.0078125,124:0.9609375):0.0078125,125:0.96875):0.0078125,126:0.9765625):0.0078125,127:0.984375):0.0078125,128:0.9921875);"
pBstr <- "(((((((1:0.125,2:0.125):0.125,(3:0.125,4:0.125):0.125):0.125,((5:0.125,6:0.125):0.125,(7:0.125,8:0.125):0.125):0.125):0.125,(((9:0.125,10:0.125):0.125,(11:0.125,12:0.125):0.125):0.125,((13:0.125,14:0.125):0.125,(15:0.125,16:0.125):0.125):0.125):0.125):0.125,((((17:0.125,18:0.125):0.125,(19:0.125,20:0.125):0.125):0.125,((21:0.125,22:0.125):0.125,(23:0.125,24:0.125):0.125):0.125):0.125,(((25:0.125,26:0.125):0.125,(27:0.125,28:0.125):0.125):0.125,((29:0.125,30:0.125):0.125,(31:0.125,32:0.125):0.125):0.125):0.125):0.125):0.125,(((((33:0.125,34:0.125):0.125,(35:0.125,36:0.125):0.125):0.125,((37:0.125,38:0.125):0.125,(39:0.125,40:0.125):0.125):0.125):0.125,(((41:0.125,42:0.125):0.125,(43:0.125,44:0.125):0.125):0.125,((45:0.125,46:0.125):0.125,(47:0.125,48:0.125):0.125):0.125):0.125):0.125,((((49:0.125,50:0.125):0.125,(51:0.125,52:0.125):0.125):0.125,((53:0.125,54:0.125):0.125,(55:0.125,56:0.125):0.125):0.125):0.125,(((57:0.125,58:0.125):0.125,(59:0.125,60:0.125):0.125):0.125,((61:0.125,62:0.125):0.125,(63:0.125,64:0.125):0.125):0.125):0.125):0.125):0.125):0.125,((((((65:0.125,66:0.125):0.125,(67:0.125,68:0.125):0.125):0.125,((69:0.125,70:0.125):0.125,(71:0.125,72:0.125):0.125):0.125):0.125,(((73:0.125,74:0.125):0.125,(75:0.125,76:0.125):0.125):0.125,((77:0.125,78:0.125):0.125,(79:0.125,80:0.125):0.125):0.125):0.125):0.125,((((81:0.125,82:0.125):0.125,(83:0.125,84:0.125):0.125):0.125,((85:0.125,86:0.125):0.125,(87:0.125,88:0.125):0.125):0.125):0.125,(((89:0.125,90:0.125):0.125,(91:0.125,92:0.125):0.125):0.125,((93:0.125,94:0.125):0.125,(95:0.125,96:0.125):0.125):0.125):0.125):0.125):0.125,(((((97:0.125,98:0.125):0.125,(99:0.125,100:0.125):0.125):0.125,((101:0.125,102:0.125):0.125,(103:0.125,104:0.125):0.125):0.125):0.125,(((105:0.125,106:0.125):0.125,(107:0.125,108:0.125):0.125):0.125,((109:0.125,110:0.125):0.125,(111:0.125,112:0.125):0.125):0.125):0.125):0.125,((((113:0.125,114:0.125):0.125,(115:0.125,116:0.125):0.125):0.125,((117:0.125,118:0.125):0.125,(119:0.125,120:0.125):0.125):0.125):0.125,(((121:0.125,122:0.125):0.125,(123:0.125,124:0.125):0.125):0.125,((125:0.125,126:0.125):0.125,(127:0.125,128:0.125):0.125):0.125):0.125):0.125):0.125):0.125);"
pCstr <- "(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((1:0.0078125,2:0.0078125):0.0078125,3:0.0078125):0.0078125,4:0.0078125):0.0078125,5:0.0078125):0.0078125,6:0.0078125):0.0078125,7:0.0078125):0.0078125,8:0.0078125):0.0078125,9:0.0078125):0.0078125,10:0.0078125):0.0078125,11:0.0078125):0.0078125,12:0.0078125):0.0078125,13:0.0078125):0.0078125,14:0.0078125):0.0078125,15:0.0078125):0.0078125,16:0.0078125):0.0078125,17:0.0078125):0.0078125,18:0.0078125):0.0078125,19:0.0078125):0.0078125,20:0.0078125):0.0078125,21:0.0078125):0.0078125,22:0.0078125):0.0078125,23:0.0078125):0.0078125,24:0.0078125):0.0078125,25:0.0078125):0.0078125,26:0.0078125):0.0078125,27:0.0078125):0.0078125,28:0.0078125):0.0078125,29:0.0078125):0.0078125,30:0.0078125):0.0078125,31:0.0078125):0.0078125,32:0.0078125):0.0078125,33:0.0078125):0.0078125,34:0.0078125):0.0078125,35:0.0078125):0.0078125,36:0.0078125):0.0078125,37:0.0078125):0.0078125,38:0.0078125):0.0078125,39:0.0078125):0.0078125,40:0.0078125):0.0078125,41:0.0078125):0.0078125,42:0.0078125):0.0078125,43:0.0078125):0.0078125,44:0.0078125):0.0078125,45:0.0078125):0.0078125,46:0.0078125):0.0078125,47:0.0078125):0.0078125,48:0.0078125):0.0078125,49:0.0078125):0.0078125,50:0.0078125):0.0078125,51:0.0078125):0.0078125,52:0.0078125):0.0078125,53:0.0078125):0.0078125,54:0.0078125):0.0078125,55:0.0078125):0.0078125,56:0.0078125):0.0078125,57:0.0078125):0.0078125,58:0.0078125):0.0078125,59:0.0078125):0.0078125,60:0.0078125):0.0078125,61:0.0078125):0.0078125,62:0.0078125):0.0078125,63:0.0078125):0.0078125,64:0.0078125):0.0078125,65:0.0078125):0.0078125,66:0.0078125):0.0078125,67:0.0078125):0.0078125,68:0.0078125):0.0078125,69:0.0078125):0.0078125,70:0.0078125):0.0078125,71:0.0078125):0.0078125,72:0.0078125):0.0078125,73:0.0078125):0.0078125,74:0.0078125):0.0078125,75:0.0078125):0.0078125,76:0.0078125):0.0078125,77:0.0078125):0.0078125,78:0.0078125):0.0078125,79:0.0078125):0.0078125,80:0.0078125):0.0078125,81:0.0078125):0.0078125,82:0.0078125):0.0078125,83:0.0078125):0.0078125,84:0.0078125):0.0078125,85:0.0078125):0.0078125,86:0.0078125):0.0078125,87:0.0078125):0.0078125,88:0.0078125):0.0078125,89:0.0078125):0.0078125,90:0.0078125):0.0078125,91:0.0078125):0.0078125,92:0.0078125):0.0078125,93:0.0078125):0.0078125,94:0.0078125):0.0078125,95:0.0078125):0.0078125,96:0.0078125):0.0078125,97:0.0078125):0.0078125,98:0.0078125):0.0078125,99:0.0078125):0.0078125,100:0.0078125):0.0078125,101:0.0078125):0.0078125,102:0.0078125):0.0078125,103:0.0078125):0.0078125,104:0.0078125):0.0078125,105:0.0078125):0.0078125,106:0.0078125):0.0078125,107:0.0078125):0.0078125,108:0.0078125):0.0078125,109:0.0078125):0.0078125,110:0.0078125):0.0078125,111:0.0078125):0.0078125,112:0.0078125):0.0078125,113:0.0078125):0.0078125,114:0.0078125):0.0078125,115:0.0078125):0.0078125,116:0.0078125):0.0078125,117:0.0078125):0.0078125,118:0.0078125):0.0078125,119:0.0078125):0.0078125,120:0.0078125):0.0078125,121:0.0078125):0.0078125,122:0.0078125):0.0078125,123:0.0078125):0.0078125,124:0.0078125):0.0078125,125:0.0078125):0.0078125,126:0.0078125):0.0078125,127:0.0078125):0.0078125,128:0.0078125);"
pDstr <- "(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((1:0.0078125,2:0.0078125):0.0078125,3:0.0116008332014688):0.0078125,4:0.0155916583871406):0.0078125,5:0.0187593545442356):0.0078125,6:0.0213880438928388):0.0078125,7:0.0236355600838684):0.0078125,8:0.0255989540548129):0.0078125,9:0.0273422630494552):0.0078125,10:0.0289100028577532):0.0078125,11:0.030334390020206):0.0078125,12:0.0316395000566371):0.0078125,13:0.0328438022427875):0.0078125,14:0.0339617775837302):0.0078125,15:0.0350049915886093):0.0078125,16:0.0359828286429374):0.0078125,17:0.036903008609003):0.0078125,18:0.0377719588833015):0.0078125,19:0.0385950879118784):0.0078125,20:0.0393769899318549):0.0078125,21:0.0401216007115761):0.0078125,22:0.0408323177288084):0.0078125,23:0.0415120941115399):0.0078125,24:0.0421635129313762):0.0078125,25:0.0427888465850756):0.0078125,26:0.0433901047189515):0.0078125,27:0.0439690732514026):0.0078125,28:0.0445273464075131):0.0078125,29:0.0450663532160049):0.0078125,30:0.0455873795792758):0.0078125,31:0.0460915867756631):0.0078125,32:0.0465800270645798):0.0078125,33:0.0470536569225207):0.0078125,34:0.0475133483289498):0.0078125,35:0.0479598984370744):0.0078125,36:0.0483940378992214):0.0078125,37:0.0488164380653954):0.0078125,38:0.0492277172332503):0.0078125,39:0.049628446095653):0.0078125,40:0.0500191525063924):0.0078125,41:0.0504003256639661):0.0078125,42:0.0507724197966948):0.0078125,43:0.0511358574188424):0.0078125,44:0.051491032216317):0.0078125,45:0.0518383116114016):0.0078125,46:0.0521780390484309):0.0078125,47:0.0525105360360778):0.0078125,48:0.0528361039767072):0.0078125,49:0.0531550258088966):0.0078125,50:0.0534675674855661):0.0078125,51:0.0537739793070762):0.0078125,52:0.0540744971260433):0.0078125,53:0.0543693434384093):0.0078125,54:0.054658728373412):0.0078125,55:0.054942850593497):0.0078125,56:0.0552218981138258):0.0078125,57:0.0554960490498514):0.0078125,58:0.0557654723004074):0.0078125,59:0.0560303281728737):0.0078125,60:0.0562907689562142):0.0078125,61:0.0565469394470188):0.0078125,62:0.0567989774330987):0.0078125,63:0.0570470141386813):0.0078125,64:0.057291174634807):0.0078125,65:0.057531578218141):0.0078125,66:0.0577683387610733):0.0078125,67:0.0580015650356803):0.0078125,68:0.0582313610138545):0.0078125,69:0.0584578261456767):0.0078125,70:0.0586810556178959):0.0078125,71:0.0589011405941983):0.0078125,72:0.0591181684387847):0.0078125,73:0.0593322229246294):0.0078125,74:0.0595433844276601):0.0078125,75:0.0597517301079894):0.0078125,76:0.0599573340792189):0.0078125,77:0.0601602675667444):0.0078125,78:0.0603605990559111):0.0078125,79:0.0605583944307883):0.0078125,80:0.0607537171042681):0.0078125,81:0.0609466281401304):0.0078125,82:0.0611371863676629):0.0078125,83:0.0613254484893743):0.0078125,84:0.061511469182294):0.0078125,85:0.0616953011933111):0.0078125,86:0.0618769954289694):0.0078125,87:0.0620566010400999):0.0078125,88:0.0622341655016433):0.0078125,89:0.0624097346879875):0.0078125,90:0.0625833529441181):0.0078125,91:0.0627550631528585):0.0078125,92:0.0629249067984554):0.0078125,93:0.0630929240267443):0.0078125,94:0.0632591537021137):0.0078125,95:0.0634236334614717):0.0078125,96:0.0635863997654003):0.0078125,97:0.0637474879466715):0.0078125,98:0.0639069322562904):0.0078125,99:0.0640647659072084):0.0078125,100:0.0642210211158543):0.0078125,101:0.0643757291416057):0.0078125,102:0.0645289203243274):0.0078125,103:0.064680624120086):0.0078125,104:0.0648308691351476):0.0078125,105:0.0649796831583561):0.0078125,106:0.0651270931919842):0.0078125,107:0.0652731254811416):0.0078125,108:0.0654178055418231):0.0078125,109:0.0655611581876684):0.0078125,110:0.0657032075555058):0.0078125,111:0.0658439771297455):0.0078125,112:0.0659834897656831):0.0078125,113:0.0661217677117721):0.0078125,114:0.0662588326309198):0.0078125,115:0.0663947056208564):0.0078125,116:0.066529407233626):0.0078125,117:0.0666629574942454):0.0078125,118:0.0667953759185718):0.0078125,119:0.0669266815304186):0.0078125,120:0.0670568928779594):0.0078125,121:0.0671860280494532):0.0078125,122:0.0673141046883257):0.0078125,123:0.0674411400076365):0.0078125,124:0.0675671508039636):0.0078125,125:0.0676921534707327):0.0078125,126:0.0678161640110162):0.0078125,127:0.0679391980498297):0.0078125,128:0.0680612708459476);"
newickStrings <- c(pAstr, pBstr, pCstr, pDstr)
names(newickStrings) <- paramSets

pAbl <- 1/128 # maximally imbalanced, ultrametric
pBbl <- 1/8   # perfectly balanced
pCbl <- 1/128 # maximally imbalanced, not ultrametric
pDbl <- 1/128 # maximally imbalanced, not ultrametric
branchLengths <- c(pAbl, pBbl, pCbl, pDbl)
names(branchLengths) <- paramSets

dataDir <- "."
outputDir <- paste0(dataDir, "plots/")
if(!dir.exists(outputDir)) {
  dir.create(outputDir)
}

plotWidth <- 8
plotHeight <- 5.5

