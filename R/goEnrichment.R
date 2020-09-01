###################################
## Trendy Results: GO Enrichment ##
###################################


###############
## Libraries ##
###############
suppressMessages(library(topGO))
suppressMessages(library(org.Hs.eg.db))


##########
## Data ##
##########
setwd("~/Documents/Research/Chris_Mixing_2019/Final_MixData_March25/")
load("PreprocessingAndNormalization/NormDat.RData")
load("Trendy/trendyFit.RData")
load("Trendy/GeneClasses.RData")


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
writeList <- function(listGO, saveDir) {
  system(paste0("mkdir ", saveDir))
  for(i in names(listGO)) {
    if(!is.na(listGO[[i]])) {
      write.csv(listGO[[i]]$topResults,
                file = paste0(saveDir, i, "_GOTable.csv"))
    }
  }
}


####################
## Run Enrichment ##
####################
saveDir <- paste0(getwd(), "/Trendy/EnrichmentTables/")
system(paste0("mkdir ", saveDir))

## Human ##
refGenes <- names(h100.trendy)

H100_H10_AccDec_GO <- list(accelerated=listEnrich(refGenes,H100_H10_AccDec$accelerated,"org.Hs.eg.db"),
                           decelerated=listEnrich(refGenes,H100_H10_AccDec$decelerated,"org.Hs.eg.db"))
system(paste0("mkdir ", saveDir, "/H100_H10/"))
writeList(H100_H10_AccDec_GO$accelerated, paste0(saveDir, "H100_H10/accelerated/"))
writeList(H100_H10_AccDec_GO$decelerated, paste0(saveDir, "H100_H10/decelerated/"))

H100_H85_AccDec_GO <- list(accelerated=listEnrich(refGenes,H100_H85_AccDec$accelerated,"org.Hs.eg.db"),
                           decelerated=listEnrich(refGenes,H100_H85_AccDec$decelerated,"org.Hs.eg.db"))
system(paste0("mkdir ", saveDir, "/H100_H85/"))
writeList(H100_H85_AccDec_GO$accelerated, paste0(saveDir, "H100_H85/accelerated/"))
writeList(H100_H85_AccDec_GO$decelerated, paste0(saveDir, "H100_H85/decelerated/"))

H100_H10_UpDown_GO <- list(deUp=listEnrich(refGenes,H100_H10_UpDown$deUp,"org.Hs.eg.db"),
                           deDown=listEnrich(refGenes,H100_H10_UpDown$deDown,"org.Hs.eg.db"))
writeList(H100_H10_UpDown_GO$deUp, paste0(saveDir, "H100_H10/deUp/"))
writeList(H100_H10_UpDown_GO$deDown, paste0(saveDir, "H100_H10/deDown/"))

H100_H85_UpDown_GO <- list(deUp=listEnrich(refGenes,H100_H85_UpDown$deUp,"org.Hs.eg.db"),
                           deDown=listEnrich(refGenes,H100_H85_UpDown$deDown,"org.Hs.eg.db"))
writeList(H100_H85_UpDown_GO$deUp, paste0(saveDir, "H100_H85/deUp/"))
writeList(H100_H85_UpDown_GO$deDown, paste0(saveDir, "H100_H85/deDown/"))

save(H100_H10_AccDec_GO,H100_H85_AccDec_GO,H100_H10_UpDown_GO,H100_H85_UpDown_GO,
     file="Trendy/EnrichmentResults.RData")

## Mouse ##
refGenes <- names(m100.trendy)

M100_M15_AccDec_GO <- list(accelerated=listEnrich(refGenes,M100_M15_AccDec$accelerated,"org.Mm.eg.db"),
                           decelerated=listEnrich(refGenes,M100_M15_AccDec$decelerated,"org.Mm.eg.db"))

M100_M90_AccDec_GO <- list(accelerated=listEnrich(refGenes,M100_M90_AccDec$accelerated,"org.Mm.eg.db"),
                           decelerated=listEnrich(refGenes,M100_M90_AccDec$decelerated,"org.Mm.eg.db"))

M100_M15_UpDown_GO <- list(deUp=listEnrich(refGenes,M100_M15_UpDown$deUp,"org.Mm.eg.db"),
                           deDown=listEnrich(refGenes,M100_M15_UpDown$deDown,"org.Mm.eg.db"))

M100_M90_UpDown_GO <- list(deUp=listEnrich(refGenes,M100_M90_UpDown$deUp,"org.Mm.eg.db"),
                           deDown=listEnrich(refGenes,M100_M90_UpDown$deDown,"org.Mm.eg.db"))


##################
## Save Results ##
##################
rm(list=setdiff(ls(),c("H100_H10_AccDec_GO","H100_H85_AccDec_GO","H100_H10_UpDown_GO","H100_H85_UpDown_GO",
                       "M100_M15_AccDec_GO","M100_M90_AccDec_GO","M100_M15_UpDown_GO","M100_M90_UpDown_GO")))
save(H100_H10_AccDec_GO,H100_H85_AccDec_GO,H100_H10_UpDown_GO,H100_H85_UpDown_GO,
     M100_M15_AccDec_GO,M100_M90_AccDec_GO,M100_M15_UpDown_GO,M100_M90_UpDown_GO,
     file="Trendy/EnrichmentResults.RData", 
     compress = "xz")




