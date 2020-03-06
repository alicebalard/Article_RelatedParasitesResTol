### Code for data analysis of Article 2
### August 2019
### Alice Balard

## INFO
# Mouse AA_0088, HI = 0.2
# Mouse AA_0064, HI = 0.08
# Mouse AA_0139, HI = 0.85

#### Load data and functions within ####
source("1_dataPreparation.R")
library(cowplot)
library(ggplot2)
library(dplyr)

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
# pdf(file = "../figures/Fig1.pdf", width = 8, height = 8)
# map
# dev.off()

######################################
########## Read information ##########
######################################

## Read batches (Exp_003 treated by anthelminthcs only)
table(art2al_SUMdf$infection_isolate, art2al_SUMdf$Mouse_genotype, art2al_SUMdf$Exp_ID)

## Read anthelminthics
table(art2al_SUMdf$infection_isolate, art2al_SUMdf$Mouse_genotype, art2al_SUMdf$anth)

## Make subdata, removing coinfected (N=9) and anthelminthic trt mice (N=22)
contaAnimals <- art2al_RAWdf[art2al_RAWdf$oocysts.per.tube > 0 & !is.na(art2al_RAWdf$oocysts.per.tube) &
                               art2al_RAWdf$dpi == 0, "EH_ID"]
SUBsummaryDF77mice <- art2al_SUMdf[!art2al_SUMdf$EH_ID %in% contaAnimals &
                                         art2al_SUMdf$anth == FALSE,]

## which mice?
problemMice <- art2al_SUMdf[art2al_SUMdf$EH_ID %in% contaAnimals |
               art2al_SUMdf$anth == TRUE,]

table(problemMice$Mouse_genotype, problemMice$infection_isolate)
  # 2 SCHUNT & 1 PWD falciformis, rest ferrisi
table(art2al_SUMdf$Mouse_genotype, art2al_SUMdf$infection_isolate)

## Age of mice
range(as.numeric(art2al_RAWdf$ageAtInfection))

###### what is the overall peak day for each parasite isolate? ######
aggregate(art2al_SUMdf$dpi_max.OPG,
          list(art2al_SUMdf$infection_isolate), 
          function(x) {paste(length(x), median(x), round(sd(x),2))})
aggregate(art2al_SUMdf$dpi_minWeight,
          list(art2al_SUMdf$infection_isolate), 
          function(x) {paste(length(x), median(x), round(sd(x),2))})

###### what is the overall prepatent period for each parasite isolate? ######
d <- as.data.frame(
  art2al_RAWdf[!is.na(art2al_RAWdf$OPG) & art2al_RAWdf$OPG > 0,] %>% 
    dplyr::group_by(EH_ID) %>%
    dplyr::slice(which.min(dpi)) %>%
    dplyr::select(EH_ID, weight, HI, startingWeight, ageAtInfection, Sex,
                  Mouse_genotype, Eimeria_species, Mouse_subspecies,
                  infection_isolate, Exp_ID, dpi))
aggregate(d$dpi,
          list(d$infection_isolate), 
          function(x) {paste(length(x), median(x), round(sd(x),2))})

###### Course of infection FIGURE 2 ######
forplot <- art2al_RAWdf %>%
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
  ylab("OPG (x10^6)") +
  scale_x_continuous(breaks = 0:11, name = "days post infection") +
  theme(legend.position = c(0.25, 0.8)) +
  labs(color = "Eimeria isolate") 
F2.1

forplot2 <- art2al_RAWdf %>%
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
  theme(legend.position = c(0.25, 0.2)) +
  labs(color = "Eimeria isolate") 
F2.2

Fig2 <- cowplot::plot_grid(F2.1, F2.2,
                  labels=c("A", "B"), label_size = 20)
# pdf(file = "../figures/Fig2.pdf", width = 10, height = 5)
Fig2
# dev.off()

## Correlation sum of oocysts / peak oocysts
ggplot(art2al_SUMdf, aes(sumoocysts.per.tube, max.oocysts.per.tube)) + 
  geom_smooth(method = "lm")+ geom_point()

cor(art2al_SUMdf$sumoocysts.per.tube , art2al_SUMdf$max.oocysts.per.tube,
    method = "pearson")

###############################################################
########## Define our indexes and their distribution ##########
###############################################################

## RESISTANCE: inverse of OPG
## We round
xRes <- round(as.numeric(na.omit(art2al_SUMdf$max.OPG)))
hist(xRes, breaks = 100)
descdist(xRes)
#pdf("../figures/supfig1.1.pdf")
findGoodDist(x = xRes, distribs = c("normal", "negative binomial"), 
             distribs2 = c("norm", "nbinom"))
#dev.off()
### nbinom for resistance

## IMPACT ON HEALTH
xImp <- as.numeric(na.omit(art2al_SUMdf$relWL))
hist(xImp, breaks = 100)
descdist(xImp)
#pdf("../figures/supfig1.2.pdf")
findGoodDist(x = xImp+ 0.01, distribs = c("normal", "weibull"), 
             distribs2 = c("norm", "weibull"))
#dev.off()
### weibull for impact on health
art2al_SUMdf$impact <- art2al_SUMdf$relWL + 0.01
SUBsummaryDF77mice$impact <- SUBsummaryDF77mice$relWL + 0.01

# 9 mice died before peak

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
  print(paste0("G=",round(2*dLL, 1), " ,df=", dDF, " ,p=", pvalue))
}

## LRT significance for each factor
myLRTsignificanceFactors <- function(modFull, modPar, modMouse, modInt){
  print("significance of parasite:")
  homemadeGtest(modFull, modPar)
  print("significance of mouse:")
  homemadeGtest(modFull, modMouse)
  print("significance of interaction:")
  homemadeGtest(modFull, modInt)
}

testSignif <- function(dataframe, which){
  if(which == "RES"){
    modFULL <- glm.nb(max.OPG ~ infection_isolate*Mouse_genotype, data = dataframe)
    modPara <- glm.nb(max.OPG ~ Mouse_genotype, data = dataframe)
    modMous <- glm.nb(max.OPG ~ infection_isolate, data = dataframe)
    modinter <- glm.nb(max.OPG ~ infection_isolate+Mouse_genotype, data = dataframe)
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
  return(list(modfull = modFULL, LRT = myLRTsignificanceFactors(modFULL, modPara, modMous, modinter)))
}

testSignifWithinParas <- function(dataframe, which){
  if(which == "RES"){
    modFULL <- glm.nb(max.OPG ~ Mouse_genotype, data = dataframe)
    mod0 <- glm.nb(max.OPG ~ 1, data = dataframe)
  } else if (which == "IMP"){
    modFULL <- lm(relWL ~ Mouse_genotype, data = dataframe)
    mod0 <- lm(relWL ~ 1, data = dataframe)
  } else if (which == "TOL"){
    modFULL <- lm(relWL ~ 0 + max.OPG : Mouse_genotype, data = dataframe)
    mod0 <- lm(relWL ~ 0 + max.OPG, data = dataframe)
  }
  print("significance of mouse:")
  G <- homemadeGtest(modFULL, mod0)
  return(list(modfull = modFULL, LRT = G))
}

######### STEP 1. Full model to see significance of all variables
######### STEP 2. If parasites significant, model within this infection group
######### STEP 3. If mouse significant, post-hoc test
######### Plot all

######### STEP 1. Full model to see significance of all variables

### Res
testSignif(art2al_SUMdf, "RES")$LRT 
testSignif(SUBsummaryDF77mice, "RES")$LRT # consistent
# Predicted values:
predRes <- ggpredict(testSignif(art2al_SUMdf, "RES")$modfull, terms = c("Mouse_genotype", "infection_isolate"))
predRes <- (data.frame(predRes))
names(predRes)[names(predRes) %in% c("x", "group")] <- c("Mouse_genotype", "infection_isolate")

predRes77 <- ggpredict(testSignif(SUBsummaryDF77mice, "RES")$modfull, terms = c("Mouse_genotype", "infection_isolate"))
predRes77 <- (data.frame(predRes77))
predRes77

### Imp
testSignif(art2al_SUMdf, "IMP")$LRT
testSignif(SUBsummaryDF77mice, "IMP")$LRT # consistent
predImp <- ggpredict(testSignif(art2al_SUMdf, "IMP")$modfull, terms = c("Mouse_genotype", "infection_isolate"))
predImp <- data.frame(predImp)
predImp
# art2al_SUMdf %>%
#   group_by(Mouse_genotype, infection_isolate) %>% 
#   summarise(meanImp = mean(relWL, na.rm = T))

### Tol
testSignif(art2al_SUMdf, "TOL")$LRT
testSignif(SUBsummaryDF77mice, "TOL")$LRT # consistent

# Predicted values of slopes:
predTolSlopes <- ggpredict(testSignif(art2al_SUMdf, "TOL")$modfull, 
                           terms = c("Mouse_genotype", "infection_isolate"), 
                           condition = c(max.OPG = 1000000))  ## For a million OPG
predTolSlopes <- data.frame(predTolSlopes)
names(predTolSlopes)[names(predTolSlopes) %in% c("x", "group")] <- c("Mouse_genotype", "infection_isolate")
predTolSlopes$group <- paste0(predTolSlopes$Mouse_genotype, predTolSlopes$infection_isolate)
predTolSlopes

predTolSlopes77 <- ggpredict(testSignif(SUBsummaryDF77mice, "TOL")$modfull, 
                           terms = c("Mouse_genotype", "infection_isolate"), 
                           condition = c(max.OPG = 1000000))  ## For a million OPG
predTolSlopes77 <- data.frame(predTolSlopes77)
names(predTolSlopes77)[names(predTolSlopes77) %in% c("x", "group")] <- c("Mouse_genotype", "infection_isolate")
predTolSlopes77$group <- paste0(predTolSlopes77$Mouse_genotype, predTolSlopes77$infection_isolate)
predTolSlopes77

######### STEP 2. If parasites significant, model within this infection group
### Res
sapply(levels(art2al_SUMdf$infection_isolate), function(x){
  testSignifWithinParas(art2al_SUMdf[art2al_SUMdf$infection_isolate %in% x,], "RES")$LRT})
# Brandenburg139 (E. ferrisi)     Brandenburg64 (E. ferrisi) Brandenburg88 (E. falciformis) 
# "G=4.6 ,df=3 ,p=0.201444"       "G=19 ,df=3 ,p=0.000276"     "G=11.6 ,df=3 ,p=0.008946"
sapply(levels(SUBsummaryDF77mice$infection_isolate), function(x){
  testSignifWithinParas(SUBsummaryDF77mice[SUBsummaryDF77mice$infection_isolate %in% x,], "RES")$LRT})
# consistent

### Imp
sapply(levels(art2al_SUMdf$infection_isolate), function(x){
  testSignifWithinParas(art2al_SUMdf[art2al_SUMdf$infection_isolate %in% x,], "IMP")$LRT})
# Brandenburg139 (E. ferrisi)             Brandenburg64 (E. ferrisi)         Brandenburg88 (E. falciformis) 
# "G=0.6 ,df=3 ,p=0.9035"  "G=14.6 ,df=3 ,p=0.00217" "G=18.3 ,df=3 ,p=0.000379" 
sapply(levels(SUBsummaryDF77mice$infection_isolate), function(x){
  testSignifWithinParas(SUBsummaryDF77mice[SUBsummaryDF77mice$infection_isolate %in% x,], "IMP")$LRT})
# consistent

### Tol
sapply(levels(art2al_SUMdf$infection_isolate), function(x){
  testSignifWithinParas(art2al_SUMdf[art2al_SUMdf$infection_isolate %in% x,], "TOL")$LRT})
# Brandenburg139 (E. ferrisi)           Brandenburg64 (E. ferrisi)       Brandenburg88 (E. falciformis) 
# "G=1.2 ,df=3 ,p=0.74265"   "G=4.1 ,df=3 ,p=0.25548" "G=10.3 ,df=3 ,p=0.01601" 
sapply(levels(SUBsummaryDF77mice$infection_isolate), function(x){
  testSignifWithinParas(SUBsummaryDF77mice[SUBsummaryDF77mice$infection_isolate %in% x,], "TOL")$LRT})
# consistent

######### STEP 3. If mouse significant, post-hoc test

### Res: Brandenburg64 (E. ferrisi) & Brandenburg88 (E. falciformis)
library(emmeans)
lapply(c("Brandenburg64 (E. ferrisi)", "Brandenburg88 (E. falciformis)"), function(x){
  lsmeans(glm.nb(max.OPG ~ Mouse_genotype, 
                      data = art2al_SUMdf[art2al_SUMdf$infection_isolate %in% x,]),
               pairwise ~ Mouse_genotype, adjust = "tukey")
})
lapply(c("Brandenburg64 (E. ferrisi)", "Brandenburg88 (E. falciformis)"), function(x){
  lsmeans(glm.nb(max.OPG ~ Mouse_genotype, 
                 data = SUBsummaryDF77mice[SUBsummaryDF77mice$infection_isolate %in% x,]),
          pairwise ~ Mouse_genotype, adjust = "tukey")
})

### Imp: Brandenburg64 (E. ferrisi) & Brandenburg88 (E. falciformis)
lapply(c("Brandenburg64 (E. ferrisi)", "Brandenburg88 (E. falciformis)"), function(x){
  mod <- lm(relWL ~ Mouse_genotype, data = art2al_SUMdf[art2al_SUMdf$infection_isolate %in% x,])
  lsmeans(mod, pairwise ~ Mouse_genotype, adjust = "tukey")
})
lapply(c("Brandenburg64 (E. ferrisi)", "Brandenburg88 (E. falciformis)"), function(x){
  mod <- lm(relWL ~ Mouse_genotype, data = SUBsummaryDF77mice[SUBsummaryDF77mice$infection_isolate %in% x,])
  lsmeans(mod, pairwise ~ Mouse_genotype, adjust = "tukey")
})

### Tol: Brandenburg88 (E. falciformis)
mod88 <- lm(relWL ~ 0 + max.OPG : Mouse_genotype, 
            data = na.omit(art2al_SUMdf[art2al_SUMdf$infection_isolate %in% "Brandenburg88 (E. falciformis)",]))
lsmeans(mod88, pairwise ~ max.OPG : Mouse_genotype, adjust = "tukey")

mod88_77 <- lm(relWL ~ 0 + max.OPG : Mouse_genotype, 
            data = na.omit(SUBsummaryDF77mice[SUBsummaryDF77mice$infection_isolate %in% "Brandenburg88 (E. falciformis)",]))
lsmeans(mod88_77, pairwise ~ max.OPG : Mouse_genotype, adjust = "tukey") # consistent

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
get_plotR_STRAINS <- function(dataframe){
  plot_model(testSignif(dataframe, "RES")$modfull,
             type = "int", dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
    scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                       name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +
    scale_y_continuous("(predicted) maximum oocysts per gram of feces (x10e6)", 
                       breaks = seq(0, 3500000, 500000),
                       labels = seq(0, 3500000, 500000)/1000000)+
    ggtitle("Maximum parasite load \n(mean and 95%CI)") +
    xlab("Eimeria isolate") +
    theme(axis.title.x = element_text(hjust=1), axis.text=element_text(size=13)) +
    geom_text(aes(x=posx.2,y=120000,label=getNs("max.OPG", dataframe)),vjust=0)
} 

plotR_STRAINS <- get_plotR_STRAINS(art2al_SUMdf)
plotR_STRAINS
plotR_STRAINS_77mice <- get_plotR_STRAINS(SUBsummaryDF77mice)
plotR_STRAINS_77mice

############
## Impact ##

## NB: art2al_SUMdf$impact <- art2al_SUMdf$relWL + 0.01
art2al_SUMdf %>%
  group_by(Mouse_subspecies, Eimeria_species) %>% 
  summarise(meanImp = mean(impact, na.rm = T))

get_plotI_STRAINS <- function(dataframe){
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

plotI_STRAINS  <- get_plotI_STRAINS(art2al_SUMdf)
plotI_STRAINS
plotI_STRAINS_77mice <- get_plotI_STRAINS(SUBsummaryDF77mice)
plotI_STRAINS_77mice

# Fig 3.
Fig3 <- cowplot::plot_grid(plotR_STRAINS + theme(legend.position = "none"),
                           plotI_STRAINS + theme(legend.position = "none"),
                           plotR_STRAINS,
                           labels=c("A", "B", "C"), label_size = 20)

Fig3
pdf(file = "../figures/Fig3.pdf",
    width = 9, height = 9)
Fig3
dev.off()

### SUB df
pdf(file = "../figures/FigSTRAINS_77mice.pdf",
    width = 9, height = 9)
cowplot::plot_grid(
  plotR_STRAINS_77mice + theme(legend.position = "none"),
  plotI_STRAINS_77mice + theme(legend.position = "none"),
  plotR_STRAINS_77mice,
  labels=c("A", "B", "C"), label_size = 20)
dev.off()

### Tolerance

makeTolPlot <- function(mypredTolSlopes, mydata){
  # make line up to 5e6 OPG for plot
  pts <- mypredTolSlopes
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
    scale_x_continuous("maximum oocysts per gram of feces (x10e6)", 
                       breaks = seq(0, 5000000, 1000000),
                       labels = seq(0, 5000000, 1000000)/1000000) +
    scale_y_continuous(name = "maximum weight loss compared to day of infection",
                       breaks = seq(0,0.3, 0.05), 
                       labels = scales::percent_format(accuracy = 5L)) +
    geom_point(data = mydata, size = 4, pch = 1)+
    coord_cartesian(ylim=c(0, 0.30)) +
    theme(legend.position = "top")
}

T1 <- makeTolPlot(predTolSlopes, art2al_SUMdf)
pdf(file = "../figures/Fig4.pdf",
    width = 12, height = 6)
T1
dev.off()

T2 <- makeTolPlot(predTolSlopes77, SUBsummaryDF77mice)
pdf(file = "../figures/Figtolslopes77.pdf",
    width = 12, height = 6)
T2
dev.off()

#### Final: Res-Tol plot to illustrate

# 12 groups
names(predRes)[names(predRes) %in% c("predicted", "conf.low", "std.error", "conf.high")] <- 
  paste(names(predRes)[names(predRes) %in% c("predicted", "conf.low",  "std.error", "conf.high")], "Res", sep = "_")
names(predTolSlopes)[names(predTolSlopes) %in% c("predicted", "conf.low",  "std.error", "conf.high")] <- 
  paste(names(predTolSlopes)[names(predTolSlopes) %in% c("predicted", "conf.low",  "std.error", "conf.high")], "Tol", sep = "_")

finalplotDF <- merge(predRes, predTolSlopes)
finalplot <- ggplot(finalplotDF, aes(x = predicted_Res, y = predicted_Tol)) +
  geom_point( aes(col = Mouse_genotype, pch = infection_isolate), size = 5)+
  scale_fill_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                     name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +
  geom_smooth(method = "lm", se = F, aes(group = infection_isolate, linetype = infection_isolate)) +
  scale_color_manual(values = as.character(levels(forMap$color)))  +
  scale_x_continuous("(predicted) maximum oocysts per gram of feces (x10e6)", 
                     breaks = seq(0, 3500000, 500000),
                     labels = seq(0, 3500000, 500000)/1000000)+
  scale_y_continuous("% weight loss by million OPG shed", labels = scales::percent_format(accuracy = 5L))
finalplot
  