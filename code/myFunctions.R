# Alice Balard
# useful functions for data analysis infection experiments
# March 2019

# Load libraries
if(!require(multcomp)){install.packages("multcomp")}
listLib <- c("ggplot2", "gridExtra", "reshape2", "scales", "lme4", 
             "lmerTest", "plyr", "dplyr", "tidyr", "reshape", "multcomp")
lapply(listLib, require, character.only = TRUE)

#### Drop unused levels for a full dataframe
dropLevelsAllFactorsDF <- function(df){
  data.frame(lapply(df, function(x) if (is.factor(x)) droplevels(x) else x))
}

# calculate OPG
calculateOPG <- function(ExpeDF){
  ExpeDF$mean_Neubauer <- 
    (ExpeDF$Neubauer1 + ExpeDF$Neubauer2 + ExpeDF$Neubauer3 + ExpeDF$Neubauer4) / 4
  # NB! Limit of detection = 1 oocysts
  ExpeDF$mean_Neubauer[ExpeDF$Neubauer1 + ExpeDF$Neubauer2 + ExpeDF$Neubauer3 + ExpeDF$Neubauer4 == 1] <- 0
  ExpeDF$oocysts.per.tube <- ExpeDF$mean_Neubauer * 10000 * ExpeDF$dilution_ml
  ExpeDF$OPG <- ExpeDF$oocysts.per.tube / ExpeDF$fecweight
  ## If we don't have the fecal weight BUT we counted in Neubauer chamber 0, then OPG = 0
  ExpeDF$oocysts.per.tube[ExpeDF$fecweight == 0 & ExpeDF$mean_Neubauer == 0] <- 0
  ExpeDF$OPG[ExpeDF$fecweight == 0 & ExpeDF$mean_Neubauer == 0] <- 0
  ## If we don't have the fecal weight/have it at zero BUT we counted stuff: our bad, NAs, snif
  ExpeDF$OPG[ExpeDF$OPG %in% Inf] <- NA
  ## Round for later count data analysis. 0.5 oocysts makes no sense anyway!
  ExpeDF$OPG <- round(ExpeDF$OPG)
  return(ExpeDF)
}

# create a table with tolerance factor, max.loss, max.OPG and sum.oocysts concatenated
makeSummaryTable <- function(df){
  # start at dpi 4 to avoid the second day post-traumatic infection protocol 
  dfwindow <- df[df$dpi %in% 4:11,]
  # minimum weight in g / in ratio and associated original weight
  X <- as.data.frame(
    dfwindow %>% dplyr::group_by(EH_ID) %>%
      dplyr::slice(which.min(weight)) %>%
      dplyr::select(EH_ID, weight, relWL, HI, startingWeight, ageAtInfection, Sex,
                    Mouse_genotype, Eimeria_species, Mouse_subspecies,
                    infection_isolate, Exp_ID, Batch, dpi))
  # time to min host weight loss peak
  names(X)[names(X) %in% "dpi"] = "dpi_minWeight"
  names(X)[names(X) %in% "weight"] = "minweight"
  # maximum OPG
  Y <- as.data.frame(
    df %>% dplyr::group_by(EH_ID) %>%
      dplyr::slice(which.max(OPG)) %>%
      dplyr::select(EH_ID, OPG, oocysts.per.tube, fecweight, dpi))
  # time to parasite shedding peak
  names(Y)[names(Y) %in% "dpi"] = "dpi_max.OPG"
  names(Y)[names(Y) %in% "OPG"] = "max.OPG"
  names(Y)[names(Y) %in% "oocysts.per.tube"] = "max.oocysts.per.tube"
  #round
  Y$max.OPG <- round(Y$max.OPG)
  # sum oocysts shed along full infection
  Z <- as.data.frame(
    df %>% dplyr::group_by(EH_ID) %>%
      dplyr::summarise(sumoocysts.per.tube = sum(oocysts.per.tube, na.rm = T)))
  # NEW FOR SLOPE: oocysts at dpi0
  W <- df[df$dpi == 0, c("EH_ID", "OPG", "dpi")]
  names(W)[names(W) %in% "dpi"] = "dpi_start.OPG"
  names(W)[names(W) %in% "OPG"] = "start.OPG"
    # merge
  fullDF <- merge(X, Y); fullDF <- merge(fullDF, Z); fullDF <- merge(fullDF, W)
  return(fullDF)
}

findGoodDist <- function(x, distribs){ 
  fits = lapply(distribs, function(i){
    fitdistrplus::fitdist(data = x, distr = i)})
  names(fits) = distribs
  par(mfrow=c(2,2))
  denscomp(fits, addlegend = T, legendtext=distribs)
  cdfcomp(fits, addlegend = T, legendtext=distribs)
  qqcomp(fits, addlegend = T, legendtext=distribs)
  ppcomp(fits, addlegend = T, legendtext=distribs)
  par(mfrow=c(1,1))
  return(list(fits = fits, aic = lapply(fits, function(x) x$aic), loglik = lapply(fits, function(x) x$loglik)))
} 

# calculate weightloss
calculateWeightLoss <- function(x, startingDay = 0){
  # define weight at infection
  A = x[x$dpi == startingDay, c("weight", "EH_ID")]
  names(A)[1] = "startingWeight"
  x = merge(x, A)
  x$relWL <- (x$startingWeight - x$weight)/x$startingWeight
  x$relWL[x$relWL < 0] <- 0
  return(x)
}
# 
# # calculate maximum realtive weight loss
# getMaxLoss <- function(df){
#   max.loss <- do.call("rbind", by(df, df$EH_ID, function (x){
#     m.loss <- which(x$relativeWeight == min(x$relativeWeight, na.rm=TRUE))
#     x[m.loss,]
#   }))
#   max.loss <- max.loss[!duplicated(max.loss$EH_ID),]
#   names(max.loss)[names(max.loss) %in% "dpi"] <- "dpi_maxweightloss"
#   names(max.loss)[names(max.loss) %in% "relativeWeight"] <- "minRelativeWeight"
#   return(max.loss)
# }
# 
# # calculate day with higher shedding peak
# getMaxOPG <- function(df){
#   max.opg <- do.call("rbind", by(df, df$EH_ID, function (x){
#     m.opg <- which(x$OPG == max(x$OPG, na.rm=TRUE))
#     x[m.opg,]
#   }))
#   max.opg <- max.opg[!duplicated(max.opg$EH_ID),]
#   names(max.opg)[names(max.opg) %in% "dpi"] <- "dpi_maxOPG"
#   names(max.opg)[names(max.opg) %in% "OPG"] <- "maxOPG"
#   max.opg$maxOPG_inmillion = max.opg$maxOPG/1e6
#   return(max.opg)
# }
# 
# # to take with care if the animal is still infected when sacrificed...
# getSumOPG <- function(df){
#   all.sum.opg <- do.call("rbind", by(df, df$EH_ID, function (x){
#     x$sum.opg <- sum(x$OPG, na.rm = TRUE)
#     x
#   }))
#   sumOpg <- all.sum.opg[!duplicated(all.sum.opg$EH_ID),]
#   sumOpg$sum.oocysts_inmillion = sumOpg$sum.opg/1e6
#   return(sumOpg)
# }

#https://rdrr.io/github/ProjectMOSAIC/mosaic/src/R/Tukey.R
# Compute Tukey Honest Significant Differences
# Create a set of confidence intervals on the differences between the means of the levels of a factor with the specified family-wise probability of coverage. The intervals are based on the Studentized range statistic, Tukey's ‘Honest Significant Difference’ method.
TukeyHSD.lm <- function(x, which, ordered = FALSE, conf.level=0.95, ...) {
  stats::TukeyHSD( aov(x), which = which, ordered = ordered, conf.level = conf.level, ...)
}

mytukey <- function(m1){
  mytukey = TukeyHSD.lm(m1)
  mytukeyDF = lapply(mytukey, function(x) {
    df <- data.frame(x)
    data.frame(combi = rownames(df[df$p.adj < 0.05 & !is.na(df$p.adj),]),
               padj = df$p.adj[df$p.adj < 0.05 & !is.na(df$p.adj)],
               diff = df$diff[df$p.adj < 0.05 & !is.na(df$p.adj)])})
  return(mytukeyDF)
}

### Create matrix for posthoc test
#   getMatrixPostHoc <- function(postHoc){
#   pvalues <- postHoc$test$pvalues
#   pvalues <- unlist(lapply(pvalues, function(i){if(i < 0.001) {
#     i <- "< 0.001"
#   } else if (i < 0.01){
#     i <- "< 0.01"
#   } else {i <- round(i, 2)}}))
#   
#   upperTriangle <- paste0("est:", round(postHoc$test$coefficients, 2), ",\nStd.Error:",round(postHoc$test$sigma, 2))
#   lowerTriangle <- paste0("z value:", round(postHoc$test$tstat, 2), ",\nPr(>|z|):", pvalues)
#   getMatrix <- function(postHoc, toFillUpWith){
#     x <- strsplit(names(postHoc$test$coefficients), " - ")
#     rownames <- as.character(unlist(lapply(x, `[[`, 1)))
#     colnames <- as.character(unlist(lapply(x, `[[`, 2)))
#     myDf_PostHoc <- data.frame(colnames = colnames, rownames = rownames, toFillUpWith = toFillUpWith)
#     
#     #fill gaps to have 12 rows and 12 columns
#     myDf_PostHoc <- rbind(data.frame(colnames = myDf_PostHoc$colnames[1],
#                                      rownames = myDf_PostHoc$colnames[1], toFillUpWith = "diagonal"),
#                           myDf_PostHoc,
#                           data.frame(colnames = myDf_PostHoc$rownames[length(myDf_PostHoc$rownames)],
#                                      rownames = myDf_PostHoc$rownames[length(myDf_PostHoc$rownames)], toFillUpWith = "diagonal"))
#     
#     # Create matrix
#     mat <- reshape(myDf_PostHoc, idvar = "colnames", timevar = "rownames", direction = "wide")
#     rownames(mat) <- mat$colnames
#     mat <- mat[!names(mat) %in% "colnames"]
#     colnames(mat) <- gsub("toFillUpWith.", "", colnames(mat))
#     return(mat)
#   }
#   ### make matrix 1: estimates and standard errors
#   mat1 <- getMatrix(postHoc, upperTriangle)
#   ### make matrix 2: z value and P
#   mat2 <- getMatrix(postHoc, lowerTriangle)
#   mat3 <- t(mat2)
#   ## Combine
#   mat <- mat1
#   mat[] <- lapply(mat, as.character)
#   mat[lower.tri(mat)] <- mat3[lower.tri(mat3)]
#   mat <- as.matrix(mat)
#   diag(mat) <- ""
#   return(mat)
# }
