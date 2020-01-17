# Alice Balard
# useful functions for data analysis infection experiments
# March 2019

### Create matrix for posthoc test
getMatrixPostHoc <- function(postHoc){
  
  pvalues <- postHoc$test$pvalues
  pvalues <- unlist(lapply(pvalues, function(i){if(i < 0.001) {
    i <- "< 0.001"
  } else if (i < 0.01){
    i <- "< 0.01"
  } else {i <- round(i, 2)}}))
  
  upperTriangle <- paste0("est:", round(postHoc$test$coefficients, 2), ", Std.Error:",round(postHoc$test$sigma, 2))
  lowerTriangle <- paste0("z value:", round(postHoc$test$tstat, 2), ", Pr(>|z|):", pvalues)
  getMatrix <- function(postHoc, toFillUpWith){
    x <- strsplit(names(postHoc$test$coefficients), " - ")
    rownames <- as.character(unlist(lapply(x, `[[`, 1)))
    colnames <- as.character(unlist(lapply(x, `[[`, 2)))
    myDf_PostHoc <- data.frame(colnames = colnames, rownames = rownames, toFillUpWith = toFillUpWith)
    
    #fill gaps to have 12 rows and 12 columns
    myDf_PostHoc <- rbind(data.frame(colnames = myDf_PostHoc$colnames[1],
                                     rownames = myDf_PostHoc$colnames[1], toFillUpWith = "diagonal"),
                          myDf_PostHoc,
                          data.frame(colnames = myDf_PostHoc$rownames[length(myDf_PostHoc$rownames)],
                                     rownames = myDf_PostHoc$rownames[length(myDf_PostHoc$rownames)], toFillUpWith = "diagonal"))
    
    # Create matrix
    mat <- reshape(myDf_PostHoc, idvar = "colnames", timevar = "rownames", direction = "wide")
    rownames(mat) <- mat$colnames
    mat <- mat[!names(mat) %in% "colnames"]
    colnames(mat) <- gsub("toFillUpWith.", "", colnames(mat))
    return(mat)
  }
  ### make matrix 1: estimates and standard errors
  mat1 <- getMatrix(postHoc, upperTriangle)
  ### make matrix 2: z value and P
  mat2 <- getMatrix(postHoc, lowerTriangle)
  mat3 <- t(mat2)
  ## Combine
  mat <- mat1
  mat[] <- lapply(mat, as.character)
  mat[lower.tri(mat)] <- mat3[lower.tri(mat3)]
  mat <- as.matrix(mat)
  diag(mat) <- ""
  return(mat)
}


# Load libraries
if(!require(multcomp)){install.packages("multcomp")}
listLib <- c("ggplot2", "gridExtra", "reshape2", "scales", "lme4", 
             "lmerTest", "plyr", "dplyr", "tidyr", "reshape", "multcomp")
lapply(listLib, require, character.only = TRUE)

makeMiceGenotypeAndIsolate <- function(df){
  # NB. Let's not consider which parent is which, but make A_B mouse = B_A mouse
  # we don't have enough individuals to test this effect, and we are not interested in it anyway!
  df$Mouse_strain <- as.character(df$Mouse_strain)
  x <- strsplit(df$Mouse_strain, "_")
  y <- lapply(x, sort)
  z <- unlist(lapply(y, FUN = function(x){paste(x, collapse="-")}))
  df$Mouse_genotype <- z
  ## Order the levels to be more clear in later plots (parents will be low and down 
  ## on the legend, hybrids in between...)
  df$Mouse_genotype <- factor(df$Mouse_genotype,
                              levels = c("NMRI", 
                                         "WSB", "WP", "PWD1",
                                         "SCHUNT-SCHUNT", "STRA-STRA",
                                         "SCHUNT-STRA",
                                         "BUSNA-STRA","PWD-SCHUNT",
                                         "BUSNA-PWD",
                                         "BUSNA-BUSNA", "PWD-PWD"),
                              labels = c("NMRI", 
                                         "MMd_F0 (Ws-Ws)", "Mmm-Mmd_Hybrid (WP)", "MMm_F0 (Pw1-Pw1)",
                                         "MMd_F0 (Sc-Sc)", "MMd_F0 (St-St)", 
                                         "MMd_F1 (Sc-St)",
                                         "Mmm-Mmd_F1 (Bu-St)", "Mmm-Mmd_F1 (Pw-Sc)",
                                         "MMm_F1 (Bu-Pw)",
                                         "MMm_F0 (Bu-Bu)", "MMm_F0 (Pw-Pw)"))
  df$Mouse_subspecies <- NA
  df$Mouse_subspecies[df$Mouse_genotype %in% "NMRI"] <- "NMRI"
  df$Mouse_subspecies[df$Mouse_genotype %in% c("MMd_F0 (Ws-Ws)", "MMd_F0 (St-St)", "MMd_F0 (Sc-Sc)", "MMd_F1 (Sc-St)")] <- "M.m.dom"
  df$Mouse_subspecies[df$Mouse_genotype %in% c("MMm_F0 (Bu-Bu)", "MMm_F0 (Pw-Pw)", "MMm_F0 (Pw1-Pw1)", "MMm_F1 (Bu-Pw)")] <- "M.m.mus"
  df$Mouse_subspecies[df$Mouse_genotype %in% c("Mmm-Mmd_Hybrid (WP)", "Mmm-Mmd_F1 (Bu-St)", "Mmm-Mmd_F1 (Pw-Sc)")] <- "Hybrid_mus_dom"
  
  df$crossingLevel <- NA
  df$crossingLevel[df$Mouse_genotype %in% "NMRI"] <- "inbredNMRI"
  df$crossingLevel[grep("F0", df$Mouse_genotype)] <- "F0"
  df$crossingLevel[grep("F1", df$Mouse_genotype)] <- "F1"
  
  df$infection_isolate <- factor(df$infection_isolate,
                                 levels = c("E139", "E64", "E88", "EfLab"),
                                 labels = c("E.ferrisi (E139)",
                                            "E.ferrisi (E64)",
                                            "E.falciformis (E88)",
                                            "E.falciformis (EfLab)"))
  # erase useless level
  df$infection_isolate <- droplevels(df$infection_isolate)
  df$Mouse_genotype <- droplevels(df$Mouse_genotype)
  df$Mouse_subspecies <- factor(df$Mouse_subspecies, levels = c("M.m.dom", "Hybrid_mus_dom", "M.m.mus"))
  df$Mouse_subspecies <- droplevels(as.factor(df$Mouse_subspecies))
  df$crossingLevel <- droplevels(as.factor(df$crossingLevel))
  
  # make average HI
  df <- merge(df, forMap[c("Mouse_genotype", "HI")], all.x = T)
  df[grep(pattern = "Sc-St", df$Mouse_genotype), "HI"] <-
    sum(forMap[grep(pattern = "Sc", forMap$Mouse_genotype), "HI"],
        forMap[grep(pattern = "St", forMap$Mouse_genotype), "HI"]) / 2
  df[grep(pattern = "Bu-Pw", df$Mouse_genotype), "HI"] <-
    sum(forMap[grep(pattern = "Bu", forMap$Mouse_genotype), "HI"],
        forMap[grep(pattern = "Pw", forMap$Mouse_genotype), "HI"]) / 2
  df[grep(pattern = "Bu-St", df$Mouse_genotype), "HI"] <-
    sum(forMap[grep(pattern = "Bu", forMap$Mouse_genotype), "HI"],
        forMap[grep(pattern = "St", forMap$Mouse_genotype), "HI"]) / 2
  df[grep(pattern = "Pw-Sc", df$Mouse_genotype), "HI"] <-
    sum(forMap[grep(pattern = "Pw", forMap$Mouse_genotype), "HI"],
        forMap[grep(pattern = "Sc", forMap$Mouse_genotype), "HI"]) / 2
  
  # Manual correction
  df$fecweight[df$fecweight > 100 & !is.na(df$fecweight)] <- 
    df$fecweight[df$fecweight > 100 & !is.na(df$fecweight)]/1000
  df$fecweight[df$fecweight > 10 & !is.na(df$fecweight)] <- 
    df$fecweight[df$fecweight > 10 & !is.na(df$fecweight)] / 10
  # Calculate oocysts
  df$mean_Neubauer <- 
    (df$Neubauer1 + df$Neubauer2 + df$Neubauer3 + df$Neubauer4) / 4
  # NB! Limit of detection = 1 oocysts
  df$mean_Neubauer[df$Neubauer1 + df$Neubauer2 + df$Neubauer3 + df$Neubauer4 == 1] <- 0
  df$oocysts.per.tube <- df$mean_Neubauer * 10000 * df$dilution_ml
  
  return(df)
}

# create a table with tolerance factor, max.loss, max.OPG and sum.oocysts concatenated
makeSummaryTable <- function(df){
  # minimum weight in g / in ratio and associated original weight
  X <- as.data.frame(
    df %>% dplyr::group_by(EH_ID) %>%
      dplyr::slice(which.min(weight)) %>%
      dplyr::select(EH_ID, weight, HI, startingWeight, ageAtInfection, Sex,
                    Mouse_genotype, Eimeria_species, Mouse_subspecies,
                    infection_isolate, Exp_ID, dpi))
  # time to min host weight loss peak
  names(X)[names(X) %in% "dpi"] = "dpi_minWeight"
  names(X)[names(X) %in% "weight"] = "minweight"
  # maximum oocysts and associated fecweight
  Y <- as.data.frame(
    df %>% dplyr::group_by(EH_ID) %>%
      dplyr::slice(which.max(oocysts.per.tube)) %>%
      dplyr::select(EH_ID, oocysts.per.tube, fecweight, dpi, weight))
  # time to parasite shedding peak
  names(Y)[names(Y) %in% "dpi"] = "dpi_max.oocysts.per.tube"
  names(Y)[names(Y) %in% "weight"] = "weight_at_dpi_max.oocysts.per.tube"
  names(Y)[names(Y) %in% "oocysts.per.tube"] = "max.oocysts.per.tube"
  Y$max.OPG <- Y$max.oocysts.per.tube / Y$fecweight
  #round
  Y$max.OPG <- round(Y$max.OPG)
    # sum oocysts shed along full infection
  Z <- as.data.frame(
    df %>% dplyr::group_by(EH_ID) %>%
      dplyr::summarise(sumoocysts.per.tube = sum(oocysts.per.tube, na.rm = T)))
  # merge
  fullDF <- merge(X, Y); fullDF <- merge(fullDF, Z)
  
  # Calculate relative weight loss
  fullDF$relWL <- (fullDF$startingWeight - fullDF$minweight)/fullDF$startingWeight
  fullDF$relWL[fullDF$relWL < 0] <- 0

  return(fullDF)
}

# Define function to be used to test, get the log lik and aic
tryDistrib <- function(x, distrib){
  # deals with fitdistr error:
  fit <- tryCatch(MASS::fitdistr(x, distrib), error=function(err) "fit failed")
  return(list(fit = fit,
              loglik = tryCatch(fit$loglik, error=function(err) "no loglik computed"), 
              AIC = tryCatch(fit$aic, error=function(err) "no aic computed")))
}

findGoodDist <- function(x, distribs, distribs2){
  l =lapply(distribs, function(i) tryDistrib(x, i))
  names(l) <- distribs
  print(l)
  listDistr <- lapply(distribs2, function(i) fitdistrplus::fitdist(x,i))
  par(mfrow=c(2,2))
  denscomp(listDistr, legendtext=distribs2)
  cdfcomp(listDistr, legendtext=distribs2)
  qqcomp(listDistr, legendtext=distribs2)
  ppcomp(listDistr, legendtext=distribs2)
  par(mfrow=c(1,1))
}


# calculate weightloss
calculateWeightLoss <- function(x, startingDay = 0){
  # define weight at infection
  A = x[x$dpi == startingDay, c("weight", "EH_ID")]
  names(A)[1] = "startingWeight"
  x = merge(x, A)
  # Cumulative sum of our weight loss
  x = x %>% 
    group_by(EH_ID) %>% 
    dplyr::arrange(dpi, .by_group = TRUE) %>%
    dplyr::mutate(relativeWeight = weight / startingWeight * 100)
  x = data.frame(x)
  return(x)
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
  ExpeDF$OPG <- round(ExpeDF$OPG)
  return(ExpeDF)
}

# calculate maximum realtive weight loss
getMaxLoss <- function(df){
  max.loss <- do.call("rbind", by(df, df$EH_ID, function (x){
    m.loss <- which(x$relativeWeight == min(x$relativeWeight, na.rm=TRUE))
    x[m.loss,]
  }))
  max.loss <- max.loss[!duplicated(max.loss$EH_ID),]
  names(max.loss)[names(max.loss) %in% "dpi"] <- "dpi_maxweightloss"
  names(max.loss)[names(max.loss) %in% "relativeWeight"] <- "minRelativeWeight"
  return(max.loss)
}

# calculate day with higher shedding peak
getMaxOPG <- function(df){
  max.opg <- do.call("rbind", by(df, df$EH_ID, function (x){
    m.opg <- which(x$OPG == max(x$OPG, na.rm=TRUE))
    x[m.opg,]
  }))
  max.opg <- max.opg[!duplicated(max.opg$EH_ID),]
  names(max.opg)[names(max.opg) %in% "dpi"] <- "dpi_maxOPG"
  names(max.opg)[names(max.opg) %in% "OPG"] <- "maxOPG"
  max.opg$maxOPG_inmillion = max.opg$maxOPG/1e6
  return(max.opg)
}

# to take with care if the animal is still infected when sacrificed...
getSumOPG <- function(df){
  all.sum.opg <- do.call("rbind", by(df, df$EH_ID, function (x){
    x$sum.opg <- sum(x$OPG, na.rm = TRUE)
    x
  }))
  sumOpg <- all.sum.opg[!duplicated(all.sum.opg$EH_ID),]
  sumOpg$sum.oocysts_inmillion = sumOpg$sum.opg/1e6
  return(sumOpg)
}

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
