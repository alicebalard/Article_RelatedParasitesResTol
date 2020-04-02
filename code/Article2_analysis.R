### Code for data analysis of Article 2
### August 2019
### Alice Balard

## INFO
# Mouse AA_0088, HI = 0.2
# Mouse AA_0064, HI = 0.08
# Mouse AA_0139, HI = 0.85

#### Load data and functions within ####
source("dataPreparationArticle2.R")
library(cowplot)
library(ggplot2)
library(dplyr)

## Different datasets as follow:

# FULL = DSart2 / art2SummaryDF

# conservative 1 = remove mice without oocysts at peak day
# DSart2_conservative1 ; art2SummaryDF_conservative1 # 99 mice

# conservative 2 = remove mice with contamination or anthelminthic
# DSart2_conservative2 ; art2SummaryDF_conservative2 # 68 mice

###### Map of our samples FIGURE 1 (with design) ######
hmhzline <- read.csv("../data/HMHZ.csv")
# change not visible color
forMap$color <- as.character(forMap$color)
forMap$color[forMap$color %in% "green"] <- "green4"
forMap$color <- as.factor(forMap$color)

area <- get_stamenmap(bbox = c(8, 48, 18, 54), zoom = 6, maptype = "toner-lite") 
map <- ggmap(area) +
  geom_path(hmhzline, mapping =  aes(x = lon, y = lat), col = "purple", size = 1) +
  geom_label_repel(data = forMap,
                   aes(longitude, latitude, label = Name, fill = color),
                   box.padding = 2, size = 7, col = "white", segment.colour = "black") +
  geom_point(data = forMap, aes(longitude, latitude, col = color), size = 6) +
  scale_color_manual(values = as.character(levels(forMap$color))) +
  scale_fill_manual(values = as.character(levels(forMap$color))) +
  theme_bw() +
  theme(legend.position = 'none', axis.ticks=element_blank())
map 
pdf(file = "../figures/Fig1_temp.pdf", width = 8, height = 8)
map
dev.off()

###### what is the overall peak day for each parasite isolate? ######
aggregate(art2SummaryDF$dpi_max.OPG,
          list(art2SummaryDF$infection_isolate), 
          function(x) {paste(length(x), median(x), round(sd(x),2))})
# Brandenburg139 (E. ferrisi) 25 6 0.73
# Brandenburg64 (E. ferrisi) 56 6 1.01
# Brandenburg88 (E. falciformis) 27 8 1.78
aggregate(art2SummaryDF$dpi_minWeight,
          list(art2SummaryDF$infection_isolate), 
          function(x) {paste(length(x), median(x), round(sd(x),2))})
# Brandenburg139 (E. ferrisi) 25 5 2.14
# Brandenburg64 (E. ferrisi) 56 5 1.92
# Brandenburg88 (E. falciformis) 27 9 1.49

## Make table with batches (only batch 1 was treated with anthelminthics)
resume <- data.frame(table(art2SummaryDF$Batch, art2SummaryDF$Mouse_genotype, art2SummaryDF$infection_isolate))
resume <- resume[order(resume$Var1) & resume$Freq != 0,]
test <- data.frame(Batch = resume$Var1)
test$group <- paste(resume$Var2,resume$Var3)
test$freq <- resume$Freq
test <- test[order(test$Batch),]
test <- reshape(test, idvar = "Batch", v.names = "freq", timevar = "group", direction="wide")
# write.csv(test,
#           "../figures/TableAllBatches.csv", row.names = F) # NB done for FULL DS

## Age of mice
range(as.numeric(art2SummaryDF$ageAtInfection))

###### what is the overall prepatent period for each parasite isolate? ######
d <- as.data.frame(
  DSart2[!is.na(DSart2$OPG) & DSart2$OPG > 0,] %>% 
    dplyr::group_by(EH_ID) %>%
    dplyr::slice(which.min(dpi)) %>%
    dplyr::select(EH_ID, weight, HI, startingWeight, ageAtInfection, Sex,
                  Mouse_genotype, Eimeria_species, Mouse_subspecies,
                  infection_isolate, Exp_ID, dpi))
aggregate(d$dpi,
          list(d$infection_isolate), 
          function(x) {paste(length(x), median(x), round(sd(x),2))})
# 1    Brandenburg139 (E. ferrisi)  25 5 0.2
# 2     Brandenburg64 (E. ferrisi) 56 5 1.58
# 3 Brandenburg88 (E. falciformis) 23 7 2.57

###### Course of infection FIGURE 2 ######
forplot <- DSart2 %>%
  group_by(infection_isolate, dpi) %>%
  summarise(mean = mean(OPG*10e-6, na.rm = TRUE),
            sd = sd(OPG*10e-6, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)

F2.1 <- ggplot(forplot, aes(dpi, mean, group = infection_isolate, col = infection_isolate)) + 
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = .2)+
  ylab("million oocysts per gram of feces") +
  scale_x_continuous(breaks = 0:11, name = "days post infection") +
  scale_color_manual(values = c("darkgreen", "lightgreen", "orange"))+
  theme(legend.position = c(0.25, 0.8)) +
  labs(color = "Eimeria isolate") 

forplot2 <-  DSart2 %>%
  group_by(infection_isolate, dpi) %>%
  summarise(mean = mean(relativeWeight, na.rm = TRUE),
            sd = sd(relativeWeight, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)

F2.2 <- ggplot(forplot2, aes(dpi, mean, group = infection_isolate, 
                             col = infection_isolate)) + 
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = .2)+
  ylab("relative weight compared to day 0 (%)") +
  scale_x_continuous(breaks = 0:11, name = "days post infection") +
  scale_color_manual(values = c("darkgreen", "lightgreen", "orange"))+
  theme(legend.position = c(0.25, 0.2)) +
  labs(color = "Eimeria isolate") 

Fig2 <- cowplot::plot_grid(F2.1, F2.2,
                           labels=c("A", "B"), label_size = 20)

pdf(file = "../figures/Fig2_temp.pdf", width = 10, height = 5)
Fig2
dev.off()

## Correlation sum of oocysts / peak oocysts
ggplot(art2SummaryDF, aes(sumoocysts.per.tube, max.oocysts.per.tube)) + 
  geom_smooth(method = "lm")+ geom_point()

cor(art2SummaryDF$sumoocysts.per.tube , art2SummaryDF$max.oocysts.per.tube,
    method = "pearson")

###############################################################
########## Define our indexes and their distribution ##########
###############################################################
## RESISTANCE: inverse of OPG
xRes <- round(as.numeric(na.omit(art2SummaryDF$max.OPG)))
hist(xRes, breaks = 100)
findGoodDist(x = xRes, distribs = c("norm", "nbinom"))
### nbinom for resistance

################################
##### Statistical analyses #####
################################
## LRT test
homemadeGtest <- function(full, base){
  dLL = logLik(full) - logLik(base)
  dDF = base$df.residual - full$df.residual
  pvalue <- 1 - stats::pchisq(2*dLL, df=dDF)
  formatC(pvalue, format = "e", digits = 2)
  chisqvalue <- stats::qchisq(p = pvalue, df=dDF)
  return(paste0("G=",round(2*dLL, 1), " ,df=", dDF, " ,p=", signif(pvalue, digits=2)))
}
## NB. lrtest from pckage lmtest shows similar results ^^ 
## I'm reinventing the wheel again

## LRT significance for each factor
myLRTsignificanceFactors <- function(modFull, modPar, modMouse, modInt){
  # print("significance of parasite:")
  return(list(signifParasite = homemadeGtest(modFull, modPar),
              signifMouse = homemadeGtest(modFull, modMouse),
              signifInter = homemadeGtest(modFull, modInt)))
}

testSignif <- function(dataframe, which){
  if(which == "RES"){
    modFULL <- glm.nb(max.OPG ~ infection_isolate*Mouse_genotype, data = dataframe)
    modPara <- glm.nb(max.OPG ~ Mouse_genotype, data = dataframe)
    modMous <- glm.nb(max.OPG ~ infection_isolate, data = dataframe)
    modinter <- glm.nb(max.OPG ~ infection_isolate+Mouse_genotype, data = dataframe)
  } else if (which == "RES_ZI") { # for zero inflated
    modFULL <- zeroinfl(max.OPG ~ infection_isolate*Mouse_genotype, data = dataframe, dist = "negbin")
    modPara <- zeroinfl(max.OPG ~ Mouse_genotype, data = dataframe, dist = "negbin")
    modMous <- zeroinfl(max.OPG ~ infection_isolate, data = dataframe, dist = "negbin")
    modinter <- zeroinfl(max.OPG ~ infection_isolate+Mouse_genotype, data = dataframe, dist = "negbin")
  } else if (which == "IMP"){
    modFULL <- lm(relWL~infection_isolate*Mouse_genotype, data = dataframe)
    modPara <- lm(relWL~Mouse_genotype, data = dataframe)
    modMous <- lm(relWL~infection_isolate, data = dataframe)
    modinter <- lm(relWL~infection_isolate+Mouse_genotype, data = dataframe)
  } else if (which == "TOL"){
    modFULL <- lm(relWL ~ 0 + max.OPG : (infection_isolate * Mouse_genotype), data = dataframe)
    modPara <- lm(relWL ~ 0 + max.OPG : (Mouse_genotype), data = dataframe)
    modMous <- lm(relWL ~ 0 + max.OPG : (infection_isolate), data = dataframe)
    modinter <- lm(relWL ~ 0 + max.OPG : (infection_isolate + Mouse_genotype), data = dataframe)
  }
  return(list(modfull = modFULL, 
              LRT = myLRTsignificanceFactors(modFULL, modPara, modMous, modinter)))
}

testSignifWithinParas <- function(dataframe, which){
  if(which == "RES"){
    modFULL <- glm.nb(max.OPG ~ Mouse_genotype, data = dataframe)
    mod0 <- glm.nb(max.OPG ~ 1, data = dataframe)
  } else if (which == "RES_ZI"){
    modFULL <- zeroinfl(max.OPG ~ Mouse_genotype, data = dataframe, dist = "negbin")
    mod0 <- zeroinfl(max.OPG ~ 1, data = dataframe, dist = "negbin")
  } else if (which == "IMP"){
    modFULL <- lm(relWL ~ Mouse_genotype, data = dataframe)
    mod0 <- lm(relWL ~ 1, data = dataframe)
  } else if (which == "TOL"){
    modFULL <- lm(relWL ~ 0 + max.OPG : Mouse_genotype, data = dataframe)
    mod0 <- lm(relWL ~ 0 + max.OPG, data = dataframe)
  }
  G <- homemadeGtest(modFULL, mod0)
  return(list(modfull = modFULL, LRT = G))
}

######### STEP 1. Full model to see significance of all variables
######### STEP 2. If parasites significant, model within this infection group
######### STEP 3. If mouse significant, post-hoc test
######### Plot all

# to apply on our 3 DF:
MyListSumma <- list(full = art2SummaryDF, cons1 = art2SummaryDF_conservative1, 
                    cons2 = art2SummaryDF_conservative2)

######### STEP 1. Full model to see significance of all variables

### Res
lapply(MyListSumma, function(x){testSignif(x,"RES")$LRT}) # consistent
# Predicted values:
getPred <- function(x, which){
  pred <- ggpredict(testSignif(x, which)$modfull, terms = c("Mouse_genotype", "infection_isolate"))
  pred <- (data.frame(pred))
  names(pred)[names(pred) %in% c("x", "group")] <- c("Mouse_genotype", "infection_isolate")
  return(pred)}
predResList <- lapply(MyListSumma, function(x) getPred(x, "RES"))

### Imp
lapply(MyListSumma, function(x){testSignif(x,"IMP")$LRT}) # consistent
# Predicted values:
predImpList <- lapply(MyListSumma, function(x) getPred(x, "IMP"))

### Tol
lapply(MyListSumma, function(x){testSignif(x,"TOL")$LRT}) # consistent

# Predicted values of slopes:
getPredTol <- function(x){
  predTolSlopes <- ggpredict(testSignif(x, "TOL")$modfull, terms = c("Mouse_genotype", "infection_isolate"), 
                             condition = c(max.OPG = 1000000))  ## For a million OPG
  predTolSlopes <- data.frame(predTolSlopes)
  names(predTolSlopes)[names(predTolSlopes) %in% c("x", "group")] <- c("Mouse_genotype", "infection_isolate")
  predTolSlopes$group <- paste0(predTolSlopes$Mouse_genotype, predTolSlopes$infection_isolate)
  return(predTolSlopes)}

predTolList <- lapply(MyListSumma, getPredTol)

######### STEP 2. Model within each infection group
listPar <- list("Brandenburg139 (E. ferrisi)","Brandenburg64 (E. ferrisi)", "Brandenburg88 (E. falciformis)")
names(listPar) <- c("Brandenburg139", "Brandenburg64", "Brandenburg88")

### Res
lapply(MyListSumma, function(xlist){
  lapply(listPar, function(xpar){
    testSignifWithinParas(xlist[xlist$infection_isolate %in% xpar,], "RES")$LRT})
}) # consistent. 64 and 88

warnings()

test <- art2SummaryDF[art2SummaryDF$Eimeria_species %in% "E.falciformis",]
test <- art2SummaryDF
test <- art2SummaryDF[art2SummaryDF$max.OPG!=0,]
# test <- test[test$max.OPG !=0,]
modFULL <- glm.nb(max.OPG ~ Mouse_genotype, data = test)
mod0 <- glm.nb(max.OPG ~ 1, data = test)

art2SummaryDF[art2SummaryDF$max.OPG == 0, ]
  


library(pscl)
modFULL1 <- zeroinfl(max.OPG ~ Mouse_genotype, data = test, dist = "negbin")
# modFULL #EM=TRUE
mod01 <- zeroinfl(max.OPG ~ 1, data = test, dist = "negbin")
# mod0


library(lmtest)
homemadeGtest(modFULL, mod0)
lrtest(modFULL, mod0)

homemadeGtest(modFULL1, mod01)
lrtest(modFULL1, mod01)

modFULL
modFULL1

### Imp
lapply(MyListSumma, function(xlist){
  lapply(listPar, function(xpar){
    testSignifWithinParas(xlist[xlist$infection_isolate %in% xpar,], "IMP")$LRT})
}) # consistent. 64 and 88

### Tol
lapply(MyListSumma, function(xlist){
  lapply(listPar, function(xpar){
    testSignifWithinParas(xlist[xlist$infection_isolate %in% xpar,], "TOL")$LRT})
}) # consistent: only 88

######### STEP 3. If mouse significant, post-hoc test

### Res: Brandenburg64 (E. ferrisi) & Brandenburg88 (E. falciformis)
library(emmeans)

lapply(MyListSumma, function(xlist){
  lapply(listPar[c(2,3)], function(xPar){
    lsmeans(glm.nb(max.OPG ~ Mouse_genotype, 
                   data = xlist[xlist$infection_isolate %in% xPar,]),
            pairwise ~ Mouse_genotype, adjust = "tukey")
  })
})
# E64 - full + consistent conserv1 + consistent conserv2 (+STRA - PWD)
# SCHUNT - BUSNA   -0.882 0.272 Inf -3.242  0.0065 
# SCHUNT - PWD     -1.309 0.277 Inf -4.719  <.0001 

# E88 - full PB (zeros...)
# conserv1: STRA - PWD        1.632 0.384 Inf  4.254  0.0001 + consistent conserv2 + SCHUNT - STRA

### Imp: Brandenburg64 (E. ferrisi) & Brandenburg88 (E. falciformis)
lapply(MyListSumma, function(xlist){
  lapply(listPar[c(2,3)], function(xPar){
    lsmeans(lm(relWL ~ Mouse_genotype, 
               data = xlist[xlist$infection_isolate %in% xPar,]),
            pairwise ~ Mouse_genotype, adjust = "tukey")
  })
})
# E64 - full - consistent cons 1 - consistent cons 2
# SCHUNT - PWD    -0.0459 0.0162 52 -2.830  0.0324  
# STRA - PWD      -0.0579 0.0159 52 -3.632  0.0035 

# E88- full - consistent cons 1 & cons 2 (STRA - PWD only)
# STRA - BUSNA    -0.1154 0.0319 23 -3.623  0.0073 
# STRA - PWD      -0.1313 0.0319 23 -4.122  0.0022 

### Tol: Brandenburg88 (E. falciformis)
lapply(MyListSumma, function(xlist){
  mod88 <- lm(relWL ~ 0 + max.OPG : Mouse_genotype, 
              data = na.omit(xlist[xlist$infection_isolate %in% "Brandenburg88 (E. falciformis)",]))
  lsmeans(mod88, pairwise ~ max.OPG : Mouse_genotype, adjust = "tukey")
})
# E88 nothing full , cons1 & cons2 STRA - PWD

##########
## plot ##
##########
## To add Ns on top of bars
getNs <- function(proxy, df, groupMus = "Mouse_genotype", groupPar = "infection_isolate"){
  noNA = df[!is.na(df[[proxy]]),]
  noNA$groupMus = noNA[[groupMus]]
  noNA$groupPar = noNA[[groupPar]]
  tab = table(noNA$groupPar, noNA$groupMus)
  Ns = as.character(as.vector(t(tab)[as.vector(t(tab))!=0]))
  return(Ns)
}

############
## Resistance
# plot marginal effects of interaction terms by isolates & strains
posx.2 <- c(0.8+c(0,1/8,2/8,3/8),1.8+c(0,1/8,2/8,3/8),2.8+c(0,1/8,2/8,3/8))
get_plotR <- function(dataframe){
  plot_model(testSignif(dataframe, "RES")$modfull,
             type = "int", dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
    scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                       name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +
    scale_y_continuous("(predicted) maximum million oocysts per gram of feces", 
                       breaks = seq(0, 5e6, 0.5e6),
                       labels = as.character(seq(0, 5e6, 0.5e6)/1e6))+
    ggtitle("Maximum parasite load = (inverse of) resistance \n(mean and 95%CI)") +
    xlab("Eimeria isolate") +
    theme(axis.title.x = element_text(hjust=1), axis.text=element_text(size=13)) +
    geom_text(aes(x=posx.2,y=0,label=getNs("max.OPG", dataframe)),vjust=0)
} 

plotResList <- lapply(MyListSumma, function(x) get_plotR(x))

############
## Impact ##
get_plotI <- function(dataframe){
  plot_model(testSignif(dataframe, "IMP")$modfull,
             type = "int",dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
    scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                       name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +
    xlab("Eimeria isolate") +
    ggtitle("Maximum weight loss \n(mean and 95%CI)") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 5L), 
                       name = "(predicted) maximum weight loss compared to day of infection")+
    theme(axis.title.x = element_text(hjust=1), axis.text=element_text(size=13)) +
    geom_text(aes(x=posx.2,y=0,label=getNs("relWL", dataframe)),vjust=0)
}

plotImpList <- lapply(MyListSumma, function(x) get_plotI(x))

### Tolerance
get_plotT <- function(mypredTolSlopes, mydata){
  # make line up to 5e6 OPG for plot
  pts <- mypredTolSlopes
  names(pts) <- "name"
  pts <- data.frame(pts)
  colnames(pts) <- gsub("name.", "", colnames(pts))
  pts$predicted <- pts$predicted*5
  pts$relWL_OPGnull <- 0
  names(pts)[names(pts) %in% c("predicted")] <- "relWL_5MOPG"
  pts <- melt(pts, measure.vars = c("relWL_OPGnull", "relWL_5MOPG"))
  names(pts)[names(pts) %in% c("variable", "value")] <- c("max.OPG", "relWL")
  pts$max.OPG <- as.character(pts$max.OPG)
  pts$max.OPG[pts$max.OPG %in% "relWL_OPGnull"] <- "0"
  pts$max.OPG[pts$max.OPG %in% "relWL_5MOPG"] <- "5e6"
  pts$max.OPG <- as.numeric(pts$max.OPG)
  pts$group <- factor(paste0(pts$Mouse_genotype, pts$infection_isolate))
  
  ggplot(pts, aes(x = max.OPG, y = relWL, col = Mouse_genotype)) +
    geom_line(aes(group = group)) +
    facet_grid(.~infection_isolate) +
    scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                       name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +
    scale_x_continuous("maximum million oocysts per gram of feces",
                       breaks = seq(0, 5000000, 1000000),
                       labels = seq(0, 5000000, 1000000)/1000000) +
    scale_y_continuous(name = "maximum weight loss compared to day of infection",
                       breaks = seq(0,0.3, 0.05),
                       labels = scales::percent_format(accuracy = 5L)) +
    geom_point(data = mydata, size = 4, pch = 1)+
    coord_cartesian(ylim=c(0, 0.30)) +
    theme(legend.position = "top") +
    ggtitle("Tolerance \n(slope of B (max weight loss) on A (max parasite load), per genotype)")
}

plotTolList <- list(full = get_plotT(predTolList["full"], art2SummaryDF),
                    cons1 = get_plotT(predTolList["cons1"], art2SummaryDF_conservative1),
                    cons2 = get_plotT(predTolList["cons2"], art2SummaryDF_conservative2))

# Fig 3.
getFigRIT <- function(which){
  # cowplot::plot_grid(plotResList[[which]] + theme(legend.position = "none"),
  #                    plotImpList[[which]] + theme(legend.position = "none"),
  #                    plotTolList[[which]]+ theme(legend.position = "none"),
  #                    labels=c("A", "B", "C"), label_size = 20,nrow = 2,
  #                    rel_widths = c(1, 1, 2))
  cowplot::ggdraw() +
    draw_plot(plotResList[[which]] + theme(legend.position = "none"), 0, .5, .49, .49) +
    draw_plot(plotImpList[[which]] + theme(legend.position = "none"), .5, .5, .49, .49) +
    draw_plot(plotTolList[[which]]+ theme(legend.position = "none"), 0, 0, .7, .49) +
    draw_plot_label(c("A", "B", "C"), c(.01, .51, .01), c(1, 1, .5), size = 15)
}

listPlotsRIT <- lapply(c("full", "cons1", "cons2"), getFigRIT)
listPlotsRIT[1]
listPlotsRIT[2]
listPlotsRIT[3]

Fig3 <- listPlotsRIT[[1]]
Fig3
pdf(file = "../figures/Fig3_temp.pdf",
    width = 10, height = 10)
Fig3
dev.off()

#### Final: Res-Tol plots
getMergeRT <- function(x, y){
  names(x) <- "name"
  x <- data.frame(x)
  colnames(x) <- gsub("name.", "", colnames(x))
  names(x)[names(x) %in% c("predicted", "conf.low",  "std.error", "conf.high")] <- 
    paste(names(x)[names(x) %in% c("predicted", "conf.low",  "std.error", "conf.high")], "Res", sep = "_")
  names(y) <- "name"
  y <- data.frame(y)
  colnames(y) <- gsub("name.", "", colnames(y))
  names(y)[names(y) %in% c("predicted", "conf.low",  "std.error", "conf.high")] <- 
    paste(names(y)[names(y) %in% c("predicted", "conf.low",  "std.error", "conf.high")], "Tol", sep = "_")
  
  merge(x, y)
}

finalplotDF_full <- getMergeRT(predResList["full"], predTolList["full"])
finalplotDF_cons1 <- getMergeRT(predResList["cons1"], predTolList["cons1"])
finalplotDF_cons2 <- getMergeRT(predResList["cons2"], predTolList["cons2"])

list_finalplotDf <- list(finalplotDF_full = finalplotDF_full, 
                         finalplotDF_cons1 = finalplotDF_cons1,
                         finalplotDF_cons2 = finalplotDF_cons2)

## test correlations
lapply(list_finalplotDf, function(df){
  lapply(listPar, function(x){
    cor.test(df[df$infection_isolate %in% x, "predicted_Res"], 
             df[df$infection_isolate %in% x, "predicted_Tol"], 
             method="spearman")})
})

# Plot
getPlot <- function(df){
  finalplotDF_Efal <- df[grep("falciformis", df$group),]
  finalplotDF_Efer <- df[grep("ferrisi", df$group),]
  
  finalplot_Efer <- ggplot(finalplotDF_Efer, aes(x = predicted_Res, y = predicted_Tol)) +
    geom_errorbar(aes(ymin = conf.low_Tol, ymax = conf.high_Tol), color = "grey") +
    geom_errorbarh(aes(xmin = conf.low_Res, xmax = conf.high_Res), color = "grey") +
    geom_point(aes(col = Mouse_genotype, pch = infection_isolate), size = 7)+
    scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                       name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +
    scale_x_continuous("(predicted) maximum million oocysts per gram of feces",
                       breaks = seq(0, 3500000, 500000),
                       labels = seq(0, 3500000, 500000)/1000000)+
    scale_y_continuous("% weight loss by million OPG shed", labels = scales::percent_format(accuracy = 5L))
  
  finalplot_Efal <- ggplot(finalplotDF_Efal, aes(x = predicted_Res, y = predicted_Tol)) +
    geom_errorbar(aes(ymin = conf.low_Tol, ymax = conf.high_Tol), color = "grey") +
    geom_errorbarh(aes(xmin = conf.low_Res, xmax = conf.high_Res), color = "grey") +
    geom_point(aes(col = Mouse_genotype, pch = infection_isolate), size = 7)+
    scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                       name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +
    scale_x_continuous("(predicted) maximum million oocysts per gram of feces",
                       breaks = seq(0, 3500000, 500000),
                       labels = seq(0, 3500000, 500000)/1000000)+
    scale_y_continuous("% weight loss by million OPG shed", labels = scales::percent_format(accuracy = 5L))
  return(list(finalplot_Efer = finalplot_Efer, finalplot_Efal = finalplot_Efal))
}

listPlots <- lapply(list_finalplotDf, getPlot)
listPlots$finalplotDF_full$finalplot_Efer
listPlots$finalplotDF_cons1$finalplot_Efer
listPlots$finalplotDF_cons2$finalplot_Efer
listPlots$finalplotDF_full$finalplot_Efal # SCHUNT and BUSNA change places when change DF
listPlots$finalplotDF_cons1$finalplot_Efal
listPlots$finalplotDF_cons2$finalplot_Efal
  
pdf(file = "../figures/Fig4_Efer_temp.pdf",
    width = 10, height = 7)
listPlots$finalplotDF_full$finalplot_Efer
dev.off()

pdf(file = "../figures/Fig5_Efal_temp.pdf",
    width = 10, height = 7)
listPlots$finalplotDF_full$finalplot_Efal
dev.off()
