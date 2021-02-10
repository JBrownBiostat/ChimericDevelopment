###################################
## Trendy Results: GO Enrichment ##
###################################
WorkDir <- "~/Documents/"
setwd(WorkDir)


###############
## Libraries ##
###############
suppressMessages(library(topGO))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(org.Mm.eg.db))


##########
## Data ##
##########
load("FACS_Data/NormDat.RData")
load("FACS_Data/trendyFit.RData")
load("FACS_Data/GeneClasses.RData")


##########################
## Enrichment Functions ##
##########################
goEnrich <- function(refGenes, selGenes, mapping){
  geneVec <- factor(as.integer(refGenes %in% selGenes))
  names(geneVec) <- refGenes
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneVec,
                annot = annFUN.org, mapping = mapping, ID = "Symbol", nodeSize = 100)
  allGO = usedGO(object = GOdata)
  if(length(allGO) == 0){return(NA)}
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  topResults <- GenTable(GOdata, classicFisher = resultFisher, topNodes = length(allGO))
  topResults$FDR <- p.adjust(topResults$classicFisher, method = "BH")
  return(list(topResults = topResults))
}

listEnrich <- function(refGenes, selList, mapping){
  namesVec <- names(selList)
  listGO <- list()
  for(i in 1:length(namesVec)){
    if(length(selList[[i]]) > 0){
      listGO[[i]] <- goEnrich(refGenes, selList[[i]], mapping)
      names(listGO)[i] <- namesVec[i]
    }else{
      listGO[[i]] <- NA
      names(listGO)[i] <- namesVec[i]
    }
  }
  return(listGO)
}


####################
## Run Enrichment ##
####################

## Human ##
refGenes <- names(h100FACS33.trendy)

H100_HMix_AccDec_GO <- list(accelerated = listEnrich(refGenes, H100_HMix_Facs_AccDec$accelerated, "org.Hs.eg.db"),
                            decelerated = listEnrich(refGenes, H100_HMix_Facs_AccDec$decelerated, "org.Hs.eg.db"))
H100_HMix_UpDown_GO <- list(deUp = listEnrich(refGenes, H100_HMix_Facs_UpDown$deUp, "org.Hs.eg.db"),
                            deDown = listEnrich(refGenes, H100_HMix_Facs_UpDown$deDown, "org.Hs.eg.db"))

## Mouse ##
refGenes <- names(m100FACS33.trendy)

M100_MMix_AccDec_GO <- list(accelerated = listEnrich(refGenes, M100_MMix_Facs_AccDec$accelerated, "org.Mm.eg.db"),
                            decelerated = listEnrich(refGenes, M100_MMix_Facs_AccDec$decelerated, "org.Mm.eg.db"))
M100_MMix_UpDown_GO <- list(deUp = listEnrich(refGenes, M100_MMix_Facs_UpDown$deUp, "org.Mm.eg.db"),
                            deDown = listEnrich(refGenes, M100_MMix_Facs_UpDown$deDown, "org.Mm.eg.db"))


##################
## Save Results ##
##################
save(H100_HMix_AccDec_GO, H100_HMix_UpDown_GO, 
     M100_MMix_AccDec_GO, M100_MMix_UpDown_GO, 
     file="FACS_Data/EnrichmentResults.RData", 
     compress = "xz")


