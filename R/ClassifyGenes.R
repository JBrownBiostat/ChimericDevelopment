####################################
## Trendy Results: Classify Genes ##
####################################


###############
## Libraries ##
###############
suppressMessages(library(Matrix))
suppressMessages(library(Trendy))


##########
## Data ##
##########
setwd("~/Documents/Research/Chris_Mixing_2019/Final_MixData_March25/")
load("PreprocessingAndNormalization/NormDat.RData")
load("Trendy/trendyFit.RData")


#################
## FilterGenes ##
#################
dat.h100$allDay <- dat.h100$allDay[hGenes,]
dat.h85$allDay <- dat.h85$allDay[hGenes,]
dat.h10$allDay <- dat.h10$allDay[hGenes,]

dat.m100$allDay <- dat.m100$allDay[mGenes,]
dat.m90$allDay <- dat.m90$allDay[mGenes,]
dat.m15$allDay <- dat.m15$allDay[mGenes,]


####################
## Classify Genes ##
####################
## Accelerated/Decelerated ##
# Functions #
checkPeak <- function(trendy.ref, trendy.test){
  ret <- list(gene = c(), refBreak = c(), testBreak = c())
  for(i in names(trendy.ref)){
    ref <- trendy.ref[[i]]
    test <- trendy.test[[i]]
    if(!is.na(ref$Breakpoints[1]) & !is.na(test$Breakpoints[1])){
      breakRef <- 0; breakTest <- 0
      for(j in 1:length(ref$Breakpoints)) {
        if(ref$Segment.Trends[j] == 1 & ref$Segment.Trends[j + 1] %in% c(0, -1)){
          breakRef <- ref$Breakpoints[j]
          break
        }
        if(j == length(ref$Breakpoints)) {next}
      }
      for(j in 1:length(test$Breakpoints)) {
        if(test$Segment.Trends[j] == 1 & test$Segment.Trends[j + 1] %in% c(0, -1)){
          breakTest <- test$Breakpoints[j]
          break
        }
        if(j == length(test$Breakpoints)) {next}
      }
      if(breakRef != 0 & breakTest != 0) {
        l <- length(ret$gene)
        ret$gene[l+1] <- i
        ret$refBreak[l+1] <- breakRef
        ret$testBreak[l+1] <- breakTest
      }
    }
  }
  return(ret)
}

checkUpTrend <- function(trendy.ref, trendy.test){
  ret <- list(gene = c(), refBreak = c(), testBreak = c(), 
              refSlope = c(), testSlope = c(),
              refSlopeAtTest = c(), testSlopeAtRef = c())
  for(i in names(trendy.ref)){
    ref <- trendy.ref[[i]]
    test <- trendy.test[[i]]
    if(any(ref$Segment.Trends == 1) & any(test$Segment.Trends == 1)){
      if(ref$Segment.Trends[1] != 0){
        upSeg <- which.max(ref$Segment.Trends == 1)
        if(any(ref$Segment.Trends[1:(upSeg-1)] == -1)){next}
      }
      if(test$Segment.Trends[1] != 0){
        upSeg <- which.max(test$Segment.Trends == 1)
        if(any(test$Segment.Trends[1:(upSeg-1)] == -1)){next}
      }
      l <- length(ret$gene)
      ret$gene[l+1] <- i
      ret$refBreak[l+1] <- c(0, ref$Breakpoint)[which.max(ref$Segment.Trends == 1)]
      ret$testBreak[l+1] <- c(0, test$Breakpoint)[which.max(test$Segment.Trends == 1)]
      ret$refSlope[l+1] <- ref$Segment.Slopes[which.max(ref$Segment.Trends == 1)]
      ret$testSlope[l+1] <- test$Segment.Slopes[which.max(test$Segment.Trends == 1)]
      ret$refSlopeAtTest[l+1] <- ifelse(is.na(ref$Breakpoints[1]),
                                        ref$Segment.Slopes[1],
                                        ref$Segment.Slopes[which.max(c(ref$Breakpoints, 42) >= ret$testBreak[l + 1])])
      ret$testSlopeAtRef[l+1] <- ifelse(is.na(test$Breakpoints[1]),
                                        test$Segment.Slopes[1],
                                        test$Segment.Slopes[which.max(c(test$Breakpoints, 42) >= ret$testBreak[l + 1])])
    }
  }
  return(ret)
}

findAccDec_Genes <- function(trendy.ref, trendy.test){
  # Find common peaking/up trending genes
  peaks <- checkPeak(trendy.ref, trendy.test)
  upTrends <- checkUpTrend(trendy.ref, trendy.test)
  # Initialize return
  accGenes <- list(peak = c(),
                   earlyUp = c())
  decGenes <- list(peak = c(),
                   lateUp = c())
  # Peak analysis
  accGenes$peak <- peaks$gene[peaks$refBreak - peaks$testBreak >= 2]
  decGenes$peak <- peaks$gene[peaks$testBreak - peaks$refBreak >= 2]
  
  # Early/Late up trend
  accGenes$earlyUp <- upTrends$gene[(upTrends$refBreak - upTrends$testBreak >= 2) & 
                                    (upTrends$testSlope >= 5 * upTrends$refSlopeAtTest)]
  decGenes$lateUp <- upTrends$gene[(upTrends$testBreak - upTrends$refBreak >= 2) & 
                                   (upTrends$refSlope >= 5 * upTrends$testSlopeAtRef)]
  
  return(list(accelerated=accGenes,decelerated=decGenes))
}

# Comparisons #
H100_H10_AccDec <- findAccDec_Genes(h100.trendy,h10.trendy)
H100_H85_AccDec <- findAccDec_Genes(h100.trendy,h85.trendy)
M100_M15_AccDec <- findAccDec_Genes(m100.trendy,m15.trendy)
M100_M90_AccDec <- findAccDec_Genes(m100.trendy,m90.trendy)

## Differentially Expressed: Up/Down ##
# Functions #
findUpDown_Genes <- function(trendy.ref, trendy.test, day0.ref, day0.test){
  # Initialize Return
  upGenes <- list(highMaxExpr_x3 = c(),
                  deOn = c())
  downGenes <- list(lowMaxExpr_x3 = c(),
                    deOff = c())
  
  # Check for difference in maximum fitted expression
  for(i in names(trendy.ref)){
    ref <- trendy.ref[[i]]
    test <- trendy.test[[i]]
    if((max(test$Fitted.Values[-day0.test]) + 1) >= 3 * (max(ref$Fitted.Values[-day0.ref]) + 1)){
      l <- length(upGenes$highMaxExpr_x3)
      upGenes$highMaxExpr_x3[l + 1] <- i
    } else if((max(ref$Fitted.Values[-day0.ref]) + 1) >= 3 * (max(test$Fitted.Values[-day0.test]) + 1)){
      l <- length(downGenes$lowMaxExpr_x3)
      downGenes$lowMaxExpr_x3[l + 1] <- i
    }
  }
  
  for(i in names(trendy.ref)){
    ref <- trendy.ref[[i]]
    test <- trendy.test[[i]]
    if(sum((test$Fitted.Values[-day0.test] + 1) >= 3 * (max(ref$Fitted.Values[-day0.ref]) + 1)) >= 3 & 
       max(ref$Fitted.Values[-day0.ref]) <= 10){
      l <- length(upGenes$deOn)
      upGenes$deOn[l + 1] <- i
    } else if(sum((ref$Fitted.Values[-day0.ref] + 1) >= 3 * (max(test$Fitted.Values[-day0.test]) + 1)) >= 3 & 
              max(test$Fitted.Values[-day0.test]) <= 10){
      l <- length(downGenes$deOff)
      downGenes$deOff[l + 1] <- i
    }
  }
  
  return(list(deUp=upGenes,deDown=downGenes))
}

# Comparisons #
H100_H10_UpDown <- findUpDown_Genes(h100.trendy,h10.trendy, which(dat.h100$days == 0), which(dat.h10$days == 0))
H100_H85_UpDown <- findUpDown_Genes(h100.trendy,h85.trendy, which(dat.h100$days == 0), which(dat.h85$days == 0))
M100_M15_UpDown <- findUpDown_Genes(m100.trendy,m15.trendy, which(dat.m100$days == 0), which(dat.m15$days == 0))
M100_M90_UpDown <- findUpDown_Genes(m100.trendy,m90.trendy, which(dat.m100$days == 0), which(dat.m90$days == 0))


##################
## Save Results ##
##################
save(H100_H10_AccDec,H100_H85_AccDec,M100_M15_AccDec,M100_M90_AccDec,
     H100_H10_UpDown,H100_H85_UpDown,M100_M15_UpDown,M100_M90_UpDown,
     file="Trendy/GeneClasses.RData",
     compress="xz")








