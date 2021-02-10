####################################
## CoDEx deconvolution with MuSiC ##
####################################
WorkDir <- "~/Documents/"
setwd(paste0(WorkDir, "/CoDEx/"))


###############
## Libraries ##
###############
library(ggplot2)
library(gridExtra)
library(viridis)
library(MuSiC)
library(Matrix)
library(Biobase)
library(xbioc)
library(matrixStats)
library(DirichletReg)
library(splines)


##########
## Data ##
##########
load("../PreprocessingAndNormalization/NormDat.RData")
load("../Trendy/GeneClasses.RData")
load("../Trendy/trendyFit.RData")
load("../Trendy/EnrichmentResults.RData")
load("Data/raw_counts_mat.rdata")
metaDat <- read.csv("Data/cell_metadata.csv")
rownames(metaDat) <- metaDat$Cell
CoDEx <- raw_counts_mat[rowSums(raw_counts_mat > 0) >= 10, ]
CoDEx <- CoDEx[, metaDat$Cell]
rm(raw_counts_mat); gc()

set.seed(1)
subTypes <- c("vRG", "oRG", "PgS", "PgG2M", "IP", "ExN", "ExM", "ExM-U", "ExDp1", "ExDp2")
sampInd <- sample(which(metaDat$Cluster %in% subTypes), 1e4)
CoDEx.ES <- ExpressionSet(as.matrix(CoDEx[, sampInd]))

themeObj <- theme(
  title = element_text(face = "bold", size = 8),
  axis.title = element_text(size = 7),
  axis.text = element_text(size = 6), 
  legend.title = element_text(size = 6),
  legend.text = element_text(size = 6)
)


###############
## Functions ##
###############
plotProps.func <- function(DeconvDat, refDat, subTypes, title) {
  Y <- DR_data(DeconvDat$Est.prop.weighted)
  modDat <- data.frame(
    DeconvDat$Est.prop.weighted,
    day = refDat$days_jitter
  )
  modDat$Y <- Y
  mod <- DirichReg(Y ~ bs(day, df = 4), modDat)
  modPred <- predict(mod)
  
  plotDat <- data.frame(
    propVec = c(modPred),
    CellAnn = rep(colnames(DeconvDat$Est.prop.weighted), 
                  rep(nrow(DeconvDat$Est.prop.weighted), 
                      ncol(DeconvDat$Est.prop.weighted))),
    Day = rep(refDat$days_jitter, ncol(DeconvDat$Est.prop.weighted))
  )
  plotDat$CellAnn <- factor(plotDat$CellAnn, ordered = T,
                            levels = subTypes)
  p <- ggplot(plotDat, aes(x = Day, y = propVec, color = CellAnn, fill = CellAnn)) + 
    theme_classic() + 
    geom_area() + 
    labs(x = "Day", y = "Prop.", title = title, 
         color = "Neural\nstage", fill = "Neural\nstage") + 
    geom_vline(xintercept = 5, color = "red", lty = "dashed") + 
    geom_vline(xintercept = 10, color = "red", lty = "dashed") + 
    geom_vline(xintercept = 15, color = "red", lty = "dashed") + 
    geom_vline(xintercept = 20, color = "red", lty = "dashed") + 
    geom_vline(xintercept = 30, color = "red", lty = "dashed") + 
    themeObj
  p
  return(p)
}


###################
## Deconvolution ##
###################
# Human #
H100_Deconv <- music_prop(bulk.eset = ExpressionSet(as.matrix(dat.h100$allDay)), 
                          sc.eset = CoDEx.ES, 
                          clusters = metaDat$Cluster[sampInd], 
                          samples = metaDat$Cell[sampInd])

H85_Deconv <- music_prop(bulk.eset = ExpressionSet(as.matrix(dat.h85$allDay)), 
                         sc.eset = CoDEx.ES, 
                         clusters = metaDat$Cluster[sampInd], 
                         samples = metaDat$Cell[sampInd])

H10_Deconv <- music_prop(bulk.eset = ExpressionSet(as.matrix(dat.h10$allDay)), 
                         sc.eset = CoDEx.ES, 
                         clusters = metaDat$Cluster[sampInd], 
                         samples = metaDat$Cell[sampInd])
save(H100_Deconv, H85_Deconv, H10_Deconv,
     file = "Deconvolution_Human_All.RData", 
     compress = T)


p.h100 <- plotProps.func(H100_Deconv, dat.h100, subTypes, "H100")
p.h85 <- plotProps.func(H85_Deconv, dat.h85, subTypes, "H85")
p.h10 <- plotProps.func(H10_Deconv, dat.h10, subTypes, "H10")

jpeg("HumanProps.jpeg", 
     height = 9, width = 8, units = "in", res = 300, quality = 0.95)
grid.arrange(p.h100 + theme(legend.position = "none"), 
             p.h85 + theme(legend.position = "none"), 
             p.h10 + theme(legend.position = "bottom"), 
             ncol = 1, heights = c(0.75, 0.75, 1))
dev.off()









