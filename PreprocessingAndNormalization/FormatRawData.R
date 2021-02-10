#####################
## Format Raw Data ##
#####################
WorkDir <- "~/Documents/"


###############
## Libraries ##
###############
library(Matrix)


###############
## Load Data ##
###############
setwd(WorkDir)
Rep1.expr <- read.table("RawData/Sub_441_unpooled_b6187_human_and_mouse_07a9de1136075dba/genes.ec.tab", 
                        header = TRUE, sep = "\t", check.names = FALSE)
Rep1.fpkm <- read.table("RawData/Sub_441_unpooled_b6187_human_and_mouse_07a9de1136075dba/genes.fpkm.tab", 
                        header = TRUE, sep = "\t", check.names = FALSE)
Rep23.expr <- read.table("RawData/Sub_464_157df_human_and_mouse_c6dce33905b35919/genes.ec.tab", 
                         header = TRUE, sep = "\t", check.names = FALSE)
Rep23.fpkm <- read.table("RawData/Sub_464_157df_human_and_mouse_c6dce33905b35919/genes.fpkm.tab", 
                         header = TRUE, sep = "\t", check.names = FALSE)


###############
## Pool Runs ##
###############
pool1 <- Rep1.expr[, 1:2]
pool1.fpkm <- Rep1.fpkm[, 1:2]
rep1.names <- colnames(Rep1.expr)
iter <- 2; n <- ncol(Rep1.expr); run <- rep1.names[iter]; col.iter <- 2
while(iter < n) {
  iter <- iter + 1
  if(run != rep1.names[iter]) {
    col.iter <- col.iter + 1
    run <- rep1.names[iter]
    pool1 <- cbind(pool1, Rep1.expr[, iter])
    pool1.fpkm <- cbind(pool1.fpkm, Rep1.fpkm[, iter])
    colnames(pool1)[col.iter] <- run
    colnames(pool1.fpkm)[col.iter] <- run
    next
  }
  pool1[, col.iter] <- pool1[, col.iter] + Rep1.expr[, iter]
  pool1.fpkm[, col.iter] <- pool1.fpkm[, col.iter] + Rep1.fpkm[, iter]
}

pool23 <- Rep23.expr[, 1:2]
pool23.fpkm <- Rep23.fpkm[, 1:2]
rep23.names <- colnames(Rep23.expr)
iter <- 2; n <- ncol(Rep23.expr); run <- rep23.names[iter]; col.iter <- 2
while(iter < n) {
  iter <- iter + 1
  if(run != rep23.names[iter]) {
    col.iter <- col.iter + 1
    run <- rep23.names[iter]
    pool23 <- cbind(pool23, Rep23.expr[, iter])
    pool23.fpkm <- cbind(pool23.fpkm, Rep23.fpkm[, iter])
    colnames(pool23)[col.iter] <- run
    colnames(pool23.fpkm)[col.iter] <- run
    next
  }
  pool23[, col.iter] <- pool23[, col.iter] + Rep23.expr[, iter]
  pool23.fpkm[, col.iter] <- pool23.fpkm[, col.iter] + Rep23.fpkm[, iter]
}
pool2 <- pool23[, 1:101]
pool2.fpkm <- pool23.fpkm[, 1:101]
pool3 <- pool23[, c(1, 102:201)]
pool3.fpkm <- pool23.fpkm[, c(1, 102:201)]
rm(list = c("pool23", "pool23.fpkm"))


#######################
## Split Human/Mouse ##
#######################
hm.1 <- sapply(as.character(pool1$gene_id), function(x) {strsplit(x, "|", fixed = TRUE)[[1]][2]})
genes.1 <- sapply(as.character(pool1$gene_id), function(x) {strsplit(x, "|", fixed = TRUE)[[1]][1]})
H1 <- pool1[hm.1 == "hg19", -1]
H1.fpkm <- pool1.fpkm[hm.1 == "hg19", -1]
rownames(H1) <- genes.1[hm.1 == "hg19"]
rownames(H1.fpkm) <- genes.1[hm.1 == "hg19"]
M1 <- pool1[hm.1 == "mm10", -1]
M1.fpkm <- pool1.fpkm[hm.1 == "mm10", -1]
rownames(M1) <- genes.1[hm.1 == "mm10"]
rownames(M1.fpkm) <- genes.1[hm.1 == "mm10"]

hm.2 <- sapply(as.character(pool2$gene_id), function(x) {strsplit(x, "|", fixed = TRUE)[[1]][2]})
genes.2 <- sapply(as.character(pool2$gene_id), function(x) {strsplit(x, "|", fixed = TRUE)[[1]][1]})
H2 <- pool2[hm.2 == "hg19", -1]
H2.fpkm <- pool2.fpkm[hm.2 == "hg19", -1]
rownames(H2) <- genes.2[hm.2 == "hg19"]
rownames(H2.fpkm) <- genes.2[hm.2 == "hg19"]
M2 <- pool2[hm.2 == "mm10", -1]
M2.fpkm <- pool2.fpkm[hm.2 == "mm10", -1]
rownames(M2) <- genes.2[hm.2 == "mm10"]
rownames(M2.fpkm) <- genes.2[hm.2 == "mm10"]

hm.3 <- sapply(as.character(pool3$gene_id), function(x) {strsplit(x, "|", fixed = TRUE)[[1]][2]})
genes.3 <- sapply(as.character(pool3$gene_id), function(x) {strsplit(x, "|", fixed = TRUE)[[1]][1]})
H3 <- pool3[hm.3 == "hg19", -1]
H3.fpkm <- pool3.fpkm[hm.3 == "hg19", -1]
rownames(H3) <- genes.3[hm.3 == "hg19"]
rownames(H3.fpkm) <- genes.3[hm.3 == "hg19"]
M3 <- pool3[hm.3 == "mm10", -1]
M3.fpkm <- pool3.fpkm[hm.3 == "mm10", -1]
rownames(M3) <- genes.3[hm.3 == "mm10"]
rownames(M3.fpkm) <- genes.3[hm.3 == "mm10"]


#############################
## Format Percentage Mixes ##
#############################
## Replicate 1 ##
exptGroup <- sapply(colnames(H1), function(x) {strsplit(x, "_", fixed = TRUE)[[1]][1]})
H0.1 <- Matrix(as.matrix(H1[, exptGroup == "0"]))
H10.1 <- Matrix(as.matrix(H1[, exptGroup == "10"]))
H85.1 <- Matrix(as.matrix(H1[, exptGroup == "85"]))
H100.1 <- Matrix(as.matrix(H1[, exptGroup == "100"]))
HCont.1 <- Matrix(as.matrix(H1[, exptGroup == "control"]))

M100.1 <- Matrix(as.matrix(M1[, exptGroup == "0"]))
M90.1 <- Matrix(as.matrix(M1[, exptGroup == "10"]))
M15.1 <- Matrix(as.matrix(M1[, exptGroup == "85"]))
M0.1 <- Matrix(as.matrix(M1[, exptGroup == "100"]))
MCont.1 <- Matrix(as.matrix(M1[, exptGroup == "control"]))

H0.1.fpkm <- Matrix(as.matrix(H1.fpkm[, exptGroup == "0"]))
H10.1.fpkm <- Matrix(as.matrix(H1.fpkm[, exptGroup == "10"]))
H85.1.fpkm <- Matrix(as.matrix(H1.fpkm[, exptGroup == "85"]))
H100.1.fpkm <- Matrix(as.matrix(H1.fpkm[, exptGroup == "100"]))
HCont.1.fpkm <- Matrix(as.matrix(H1.fpkm[, exptGroup == "control"]))

M100.1.fpkm <- Matrix(as.matrix(M1.fpkm[, exptGroup == "0"]))
M90.1.fpkm <- Matrix(as.matrix(M1.fpkm[, exptGroup == "10"]))
M15.1.fpkm <- Matrix(as.matrix(M1.fpkm[, exptGroup == "85"]))
M0.1.fpkm <- Matrix(as.matrix(M1.fpkm[, exptGroup == "100"]))
MCont.1.fpkm <- Matrix(as.matrix(M1.fpkm[, exptGroup == "control"]))

## Replicate 2 ##
exptGroup <- sapply(colnames(H2), function(x) {strsplit(x, "_", fixed = TRUE)[[1]][1]})
H0.2 <- Matrix(as.matrix(H2[, exptGroup == "0"]))
H10.2 <- Matrix(as.matrix(H2[, exptGroup == "10"]))
H85.2 <- Matrix(as.matrix(H2[, exptGroup == "85"]))
H100.2 <- Matrix(as.matrix(H2[, exptGroup == "100"]))
HCont.2 <- Matrix(as.matrix(H2[, exptGroup == "control"]))

M100.2 <- Matrix(as.matrix(M2[, exptGroup == "0"]))
M90.2 <- Matrix(as.matrix(M2[, exptGroup == "10"]))
M15.2 <- Matrix(as.matrix(M2[, exptGroup == "85"]))
M0.2 <- Matrix(as.matrix(M2[, exptGroup == "100"]))
MCont.2 <- Matrix(as.matrix(M2[, exptGroup == "control"]))

H0.2.fpkm <- Matrix(as.matrix(H2.fpkm[, exptGroup == "0"]))
H10.2.fpkm <- Matrix(as.matrix(H2.fpkm[, exptGroup == "10"]))
H85.2.fpkm <- Matrix(as.matrix(H2.fpkm[, exptGroup == "85"]))
H100.2.fpkm <- Matrix(as.matrix(H2.fpkm[, exptGroup == "100"]))
HCont.2.fpkm <- Matrix(as.matrix(H2.fpkm[, exptGroup == "control"]))

M100.2.fpkm <- Matrix(as.matrix(M2.fpkm[, exptGroup == "0"]))
M90.2.fpkm <- Matrix(as.matrix(M2.fpkm[, exptGroup == "10"]))
M15.2.fpkm <- Matrix(as.matrix(M2.fpkm[, exptGroup == "85"]))
M0.2.fpkm <- Matrix(as.matrix(M2.fpkm[, exptGroup == "100"]))
MCont.2.fpkm <- Matrix(as.matrix(M2.fpkm[, exptGroup == "control"]))

## Replicate 3 ##
exptGroup <- sapply(colnames(H3), function(x) {strsplit(x, "_", fixed = TRUE)[[1]][1]})
H0.3 <- Matrix(as.matrix(H3[, exptGroup == "0"]))
H10.3 <- Matrix(as.matrix(H3[, exptGroup == "10"]))
H85.3 <- Matrix(as.matrix(H3[, exptGroup == "85"]))
H100.3 <- Matrix(as.matrix(H3[, exptGroup == "100"]))
HCont.3 <- Matrix(as.matrix(H3[, exptGroup == "control"]))

M100.3 <- Matrix(as.matrix(M3[, exptGroup == "0"]))
M90.3 <- Matrix(as.matrix(M3[, exptGroup == "10"]))
M15.3 <- Matrix(as.matrix(M3[, exptGroup == "85"]))
M0.3 <- Matrix(as.matrix(M3[, exptGroup == "100"]))
MCont.3 <- Matrix(as.matrix(M3[, exptGroup == "control"]))

H0.3.fpkm <- Matrix(as.matrix(H3.fpkm[, exptGroup == "0"]))
H10.3.fpkm <- Matrix(as.matrix(H3.fpkm[, exptGroup == "10"]))
H85.3.fpkm <- Matrix(as.matrix(H3.fpkm[, exptGroup == "85"]))
H100.3.fpkm <- Matrix(as.matrix(H3.fpkm[, exptGroup == "100"]))
HCont.3.fpkm <- Matrix(as.matrix(H3.fpkm[, exptGroup == "control"]))

M100.3.fpkm <- Matrix(as.matrix(M3.fpkm[, exptGroup == "0"]))
M90.3.fpkm <- Matrix(as.matrix(M3.fpkm[, exptGroup == "10"]))
M15.3.fpkm <- Matrix(as.matrix(M3.fpkm[, exptGroup == "85"]))
M0.3.fpkm <- Matrix(as.matrix(M3.fpkm[, exptGroup == "100"]))
MCont.3.fpkm <- Matrix(as.matrix(M3.fpkm[, exptGroup == "control"]))


####################################
## Match genes between replicates ##
####################################
hGenes <- sort(intersect(rownames(H100.1),rownames(H100.2)))
mGenes <- sort(intersect(rownames(M100.1),rownames(M100.2)))

H0.1 <- H0.1[hGenes,]
H10.1 <- H10.1[hGenes,]
H85.1 <- H85.1[hGenes,]
H100.1 <- H100.1[hGenes,]
HCont.1 <- HCont.1[hGenes,]

H0.1.fpkm <- H0.1.fpkm[hGenes,]
H10.1.fpkm <- H10.1.fpkm[hGenes,]
H85.1.fpkm <- H85.1.fpkm[hGenes,]
H100.1.fpkm <- H100.1.fpkm[hGenes,]
HCont.1.fpkm <- HCont.1.fpkm[hGenes,]

H0.2 <- H0.2[hGenes,]
H10.2 <- H10.2[hGenes,]
H85.2 <- H85.2[hGenes,]
H100.2 <- H100.2[hGenes,]

H0.2.fpkm <- H0.2.fpkm[hGenes,]
H10.2.fpkm <- H10.2.fpkm[hGenes,]
H85.2.fpkm <- H85.2.fpkm[hGenes,]
H100.2.fpkm <- H100.2.fpkm[hGenes,]

H0.3 <- H0.3[hGenes,]
H10.3 <- H10.3[hGenes,]
H85.3 <- H85.3[hGenes,]
H100.3 <- H100.3[hGenes,]

H0.3.fpkm <- H0.3.fpkm[hGenes,]
H10.3.fpkm <- H10.3.fpkm[hGenes,]
H85.3.fpkm <- H85.3.fpkm[hGenes,]
H100.3.fpkm <- H100.3.fpkm[hGenes,]

M100.1 <- M100.1[mGenes,]
M90.1 <- M90.1[mGenes,]
M15.1 <- M15.1[mGenes,]
M0.1 <- M0.1[mGenes,]
MCont.1 <- MCont.1[mGenes,]

M100.1.fpkm <- M100.1.fpkm[mGenes,]
M90.1.fpkm <- M90.1.fpkm[mGenes,]
M15.1.fpkm <- M15.1.fpkm[mGenes,]
M0.1.fpkm <- M0.1.fpkm[mGenes,]
MCont.1.fpkm <- MCont.1.fpkm[mGenes,]

M100.2 <- M100.2[mGenes,]
M90.2 <- M90.2[mGenes,]
M15.2 <- M15.2[mGenes,]
M0.2 <- M0.2[mGenes,]

M100.2.fpkm <- M100.2.fpkm[mGenes,]
M90.2.fpkm <- M90.2.fpkm[mGenes,]
M15.2.fpkm <- M15.2.fpkm[mGenes,]
M0.2.fpkm <- M0.2.fpkm[mGenes,]

M100.3 <- M100.3[mGenes,]
M90.3 <- M90.3[mGenes,]
M15.3 <- M15.3[mGenes,]
M0.3 <- M0.3[mGenes,]

M100.3.fpkm <- M100.3.fpkm[mGenes,]
M90.3.fpkm <- M90.3.fpkm[mGenes,]
M15.3.fpkm <- M15.3.fpkm[mGenes,]
M0.3.fpkm <- M0.3.fpkm[mGenes,]


##################
## Save results ##
##################
save(H100.1, H85.1, H10.1, H0.1, HCont.1,
     H100.1.fpkm, H85.1.fpkm, H10.1.fpkm, H0.1.fpkm, HCont.1.fpkm,
     H100.2, H85.2, H10.2, H0.2,
     H100.2.fpkm, H85.2.fpkm, H10.2.fpkm, H0.2.fpkm,
     H100.3, H85.3, H10.3, H0.3,
     H100.3.fpkm, H85.3.fpkm, H10.3.fpkm, H0.3.fpkm,
     M100.1, M90.1, M15.1, M0.1, MCont.1,
     M100.1.fpkm, M90.1.fpkm, M15.1.fpkm, M0.1.fpkm, MCont.1.fpkm,
     M100.2, M90.2, M15.2, M0.2,
     M100.2.fpkm, M90.2.fpkm, M15.2.fpkm, M0.2.fpkm,
     M100.3, M90.3, M15.3, M0.3,
     M100.3.fpkm, M90.3.fpkm, M15.3.fpkm, M0.3.fpkm,
     file="PreprocessingAndNormalization/FormatedRawData.RData",
     compress="xz")







