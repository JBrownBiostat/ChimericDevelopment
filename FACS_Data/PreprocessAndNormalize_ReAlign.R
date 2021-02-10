##################################################
## Preprocessing and normalization of FACS data ##
##################################################
WorkDir <- "~/Documents/"
setwd(WorkDir)


###############
## Libraries ##
###############
suppressMessages(library(ggplot2))
suppressMessages(library(scran))
suppressMessages(library(topGO))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(org.Mm.eg.db))


##########
## Data ##
##########
facsDat <- read.table("FACS_Data/RealignedData/genes.ec.tab", 
                      sep = "\t", header = 1)

themeObj <- theme(
  title = element_text(face = "bold", size = 8),
  axis.title = element_text(size = 7),
  axis.text = element_text(size = 6), 
  legend.title = element_text(size = 6),
  legend.text = element_text(size = 6)
)


#####################
## Pre-format data ##
#####################
allGenes <- facsDat$gene_id
speciesInd <- sapply(allGenes, function(x) {strsplit(x, "|", fixed = T)[[1]][2]})
allGenes <- sapply(allGenes, function(x) {strsplit(x, "|", fixed = T)[[1]][1]})
hInd <- speciesInd == "hg19"
mInd <- speciesInd == "mm10"

mouseDat <- list(
  EC = as.matrix(facsDat[, 2:36]),
  Day = c(0, rep(c(1:5, 8:12, 15, 17, 19, 23, 26, 29, 33), 2)),
  Mix = rep(c("M100", "Mix"), c(18, 17))
)
colnames(mouseDat$EC) <- paste0("d", mouseDat$Day)
rownames(mouseDat$EC) <- allGenes

humanDat <- list(
  EC = as.matrix(facsDat[, 37:71]),
  Day = c(0, rep(c(1:5, 8:12, 15, 17, 19, 23, 26, 29, 33), 2)),
  Mix = rep(c("H100", "Mix"), c(18, 17))
)
colnames(humanDat$EC) <- paste0("d", humanDat$Day)
rownames(humanDat$EC) <- allGenes


############################
## Check sequencing depth ##
############################
plot(mouseDat$Day, log10(colSums(mouseDat$EC)), pch = 19, cex = 0.5)

plot(humanDat$Day, log10(colSums(humanDat$EC)), pch = 19, cex = 0.5)
# remove mix day 29 (<1e3 expected counts; typical samples have >1e6)
remInd <- which(log10(colSums(humanDat$EC)) < 3)
humanDat$EC <- humanDat$EC[, -remInd]
humanDat$Day <- humanDat$Day[-remInd]
humanDat$Mix <- humanDat$Mix[-remInd]


##############################
## Check misalignment rates ##
##############################
plotDat <- data.frame(
  day = mouseDat$Day,
  errRate = colSums(mouseDat$EC[hInd, ]) / colSums(mouseDat$EC),
  expt = mouseDat$Mix
)
p.M <- ggplot(plotDat, aes(x = day, y = errRate, color = expt)) + 
  theme_classic() + 
  geom_point() + 
  ylim(c(0, 0.1)) +
  labs(x = "day", y = "total error", 
       title = "Mouse FACS empirical misalignment", 
       color = "Species\nmixture") +
  themeObj
p.M

plotDat <- data.frame(
  day = humanDat$Day,
  errRate = colSums(humanDat$EC[mInd, ]) / 
    colSums(humanDat$EC),
  expt = humanDat$Mix
)
p.H <- ggplot(plotDat, aes(x = day, y = errRate, color = expt)) + 
  theme_classic() + 
  geom_point() + 
  ylim(c(0, 0.1)) +
  labs(x = "day", y = "total error", 
       title = "Human FACS empirical misalignment", 
       color = "Species\nmixture") +
  themeObj
p.H

save(p.H, p.M, 
     file = "FACS_Data/MisalignFigs.RData", 
     compress = "xz")


#############################
## Misalignment enrichment ##
#############################
goEnrich <- function(geneVec, exprMat, mapping){
  statVec <- factor(as.integer(apply(exprMat, 1, quantile, probs = 0.8) >= 20))
  names(statVec) <- geneVec
  GOdata <- new("topGOdata", ontology = "BP", allGenes = statVec,
                annot = annFUN.org, mapping = mapping, ID = "Symbol", nodeSize = 100)
  allGO = usedGO(object = GOdata)
  if(length(allGO) == 0){return(NA)}
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  topResults <- GenTable(GOdata, classicFisher = resultFisher, topNodes = length(allGO))
  topResults$FDR <- p.adjust(topResults$classicFisher, method = "BH")
  if(any(is.na(topResults$FDR))) {
    topResults$FDR[is.na(topResults$FDR)] <- min(
      topResults$FDR[!is.na(topResults$FDR)]
    )
  }
  return(topResults)
}

hMix.misEnrich <- goEnrich(
  geneVec = rownames(humanDat$EC)[mInd],
  exprMat = humanDat$EC[mInd, humanDat$Mix == "Mix"], 
  mapping = "org.Mm.eg.db"
)

mMix.misEnrich <- goEnrich(
  geneVec = rownames(mouseDat$EC)[hInd],
  exprMat = mouseDat$EC[hInd, mouseDat$Mix == "Mix"], 
  mapping = "org.Hs.eg.db"
)

save(hMix.misEnrich, mMix.misEnrich,
     file = "FACS_Data/misalignEnrich.RData",
     compress = "xz")


############################
## Normalize observations ##
############################
mouseSF <- calculateSumFactors(mouseDat$EC[mInd, ])
mouseDat$SF <- mouseSF
mouseDat$normDat <- t(t(mouseDat$EC) / mouseSF)

humanSF <- calculateSumFactors(humanDat$EC[hInd, ])
humanDat$SF <- humanSF
humanDat$normDat <- t(t(humanDat$EC) / humanSF)


###########################
## Final format and save ##
###########################
mGeneOrd <- order(allGenes[mInd])
m100Ind <- mouseDat$Mix == "M100"
m100.facs33 <- list(
  EC = mouseDat$EC[mInd, m100Ind][mGeneOrd, ],
  normDat = mouseDat$normDat[mInd, m100Ind][mGeneOrd, ],
  Day = mouseDat$Day[m100Ind],
  SF = mouseDat$SF[m100Ind]
)

mMixInd <- mouseDat$Mix == "Mix" | mouseDat$Day == 0
mMix.facs33 <- list(
  EC = mouseDat$EC[mInd, mMixInd][mGeneOrd, ],
  normDat = mouseDat$normDat[mInd, mMixInd][mGeneOrd, ],
  Day = mouseDat$Day[mMixInd],
  SF = mouseDat$SF[mMixInd]
)

hGeneOrd <- order(allGenes[hInd])
h100Ind <- humanDat$Mix == "H100"
h100.facs33 <- list(
  EC = humanDat$EC[hInd, h100Ind][hGeneOrd, ],
  normDat = humanDat$normDat[hInd, h100Ind][hGeneOrd, ],
  Day = humanDat$Day[h100Ind],
  SF = humanDat$SF[h100Ind]
)

hMixInd <- humanDat$Mix == "Mix" | humanDat$Day == 0
hMix.facs33 <- list(
  EC = humanDat$EC[hInd, hMixInd][hGeneOrd, ],
  normDat = humanDat$normDat[hInd, hMixInd][hGeneOrd, ],
  Day = humanDat$Day[hMixInd],
  SF = humanDat$SF[hMixInd]
)

save(m100.facs33, mMix.facs33, h100.facs33, hMix.facs33,
     file = "FACS_Data/NormDat.RData",
     compress = "xz")


