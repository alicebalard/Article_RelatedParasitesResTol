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
# Fig2
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

## TOLERANCE
hist(na.omit(art2al_SUMdf$relWL / art2al_SUMdf$max.OPG), breaks = 100)
range(art2al_SUMdf$relWL / art2al_SUMdf$max.OPG, na.rm = T)
range(art2al_SUMdf$relWL / art2al_SUMdf$max.OPG + 1e-8, na.rm = T)
range(log10(art2al_SUMdf$relWL / art2al_SUMdf$max.OPG + 1e-8), na.rm = T)
range(log10(art2al_SUMdf$relWL / art2al_SUMdf$max.OPG + 1e-8)/-1, na.rm = T)
range(log10(art2al_SUMdf$relWL / art2al_SUMdf$max.OPG + 1e-8)/-8, na.rm = T)
hist(na.omit(log10(art2al_SUMdf$relWL / art2al_SUMdf$max.OPG + 1e-8)/-8), breaks = 100)

art2al_SUMdf$ToleranceIndex <- log10(
  art2al_SUMdf$relWL / art2al_SUMdf$max.OPG + 1e-8) / (-8)
SUBsummaryDF77mice$ToleranceIndex <- log10(
  SUBsummaryDF77mice$relWL / SUBsummaryDF77mice$max.OPG + 1e-8) / (-8)

plotChoiceTolIndex <- ggplot(art2al_SUMdf, 
                             aes(x=max.OPG, y =relWL, fill = ToleranceIndex))+
  geom_point(size = 4, pch =21) +
  scale_fill_gradient(low = "green", high = "red") +
  xlab("Parasite density (oocysts per mouse gram) at peak day") +
  ylab("maximumrelative weight loss compared to day 0 (%)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  labs(fill = "Tolerance Index") +
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )

plotChoiceTolIndex

xTol <- as.numeric(na.omit(art2al_SUMdf$ToleranceIndex))

# for tolerance, we need to find a way to deal with the massive extreme values
hist(xTol, breaks = 1000)
descdist(xTol)

#pdf("../figures/supfig1.3.pdf")
findGoodDist(x = xTol, distribs = c("normal"), 
             distribs2 = c("norm"))
#dev.off()

art2al_SUMdf[is.na(art2al_SUMdf$ToleranceIndex), 
                 c("EH_ID", "Mouse_genotype", "relWL", "max.OPG")]
# 9 mice died before peak
#
# pdf(file = "../figures/choiceIndexes.pdf", width = 10, height = 5)
# cowplot::plot_grid(plotChoiceResIndex,
#                    plotChoiceTolIndex,
#                    labels=c("A", "B"), label_size = 15)
# dev.off()

################################
##### Statistical analyses #####
################################
## LRT test
homemadeGtest <- function(full, base){
  dLL = logLik(full) - logLik(base)
  dDF = base$df.residual - full$df.residual
  pvalue <- 1 - stats::pchisq(2*dLL, df=dDF)
  chisqvalue <- stats::qchisq(p = pvalue, df=dDF)
  print(paste0("G=",round(2*dLL, 1), " ,df=", dDF, " ,p=", round(pvalue, 6)))
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

# for posthoc tests
art2al_SUMdf$intFacSPECIES <- interaction(art2al_SUMdf$Eimeria_species, 
                                              art2al_SUMdf$Mouse_subspecies, drop=T)
art2al_SUMdf$intFacSTRAINS <- interaction(art2al_SUMdf$infection_isolate, 
                                              art2al_SUMdf$Mouse_genotype, drop=T)
SUBsummaryDF77mice$intFacSPECIES <- interaction(SUBsummaryDF77mice$Eimeria_species, 
                                                SUBsummaryDF77mice$Mouse_subspecies, drop=T)
SUBsummaryDF77mice$intFacSTRAINS <- interaction(SUBsummaryDF77mice$infection_isolate, 
                                                SUBsummaryDF77mice$Mouse_genotype, drop=T)

testSignif <- function(dataframe, which, level){
  if(which == "RES"){
    if(level == "SPECIES"){
      modFULL <- glm.nb(max.OPG ~ Eimeria_species*Mouse_subspecies, data = dataframe)
      modPara <- glm.nb(max.OPG ~ Mouse_subspecies, data = dataframe)
      modMous <- glm.nb(max.OPG ~ Eimeria_species, data = dataframe)
      modinter <- glm.nb(max.OPG ~ Eimeria_species+Mouse_subspecies, data = dataframe)
    } else if (level == "STRAINS"){
      modFULL <- glm.nb(max.OPG ~ infection_isolate*Mouse_genotype, data = dataframe)
      modPara <- glm.nb(max.OPG ~ Mouse_genotype, data = dataframe)
      modMous <- glm.nb(max.OPG ~ infection_isolate, data = dataframe)
      modinter <- glm.nb(max.OPG ~ infection_isolate+Mouse_genotype, data = dataframe)
    }
  } else if (which == "IMP"){
    if(level == "SPECIES"){
      modFULL <- survreg(Surv(impact)~Eimeria_species*Mouse_subspecies, data = dataframe, dist="weibull")
      modPara <- survreg(Surv(impact)~Mouse_subspecies, data = dataframe, dist="weibull")
      modMous <- survreg(Surv(impact)~Eimeria_species, data = dataframe, dist="weibull")
      modinter <- survreg(Surv(impact)~Eimeria_species+Mouse_subspecies, data = dataframe, dist="weibull")
    } else if (level == "STRAINS"){
      modFULL <- survreg(Surv(impact)~infection_isolate*Mouse_genotype, data = dataframe, dist="weibull")
      modPara <- survreg(Surv(impact)~Mouse_genotype, data = dataframe, dist="weibull")
      modMous <- survreg(Surv(impact)~infection_isolate, data = dataframe, dist="weibull")
      modinter <- survreg(Surv(impact)~infection_isolate+Mouse_genotype, data = dataframe, dist="weibull")
    }
  } else if (which == "TOL"){
    if(level == "SPECIES"){
      modFULL <- lm(ToleranceIndex ~ Eimeria_species*Mouse_subspecies, data = dataframe)
      modPara <- lm(ToleranceIndex ~ Mouse_subspecies, data = dataframe)
      modMous <- lm(ToleranceIndex ~ Eimeria_species, data = dataframe)
      modinter <- lm(ToleranceIndex ~ Eimeria_species+Mouse_subspecies, data = dataframe)
    } else if (level == "STRAINS"){
      modFULL <- lm(ToleranceIndex ~ infection_isolate*Mouse_genotype, data = dataframe)
      modPara <- lm(ToleranceIndex ~ Mouse_genotype, data = dataframe)
      modMous <- lm(ToleranceIndex ~ infection_isolate, data = dataframe)
      modinter <- lm(ToleranceIndex ~ infection_isolate+Mouse_genotype, data = dataframe)
    }
  }
  return(list(modfull = modFULL, LRT = myLRTsignificanceFactors(modFULL, modPara, modMous, modinter)))
}

testPostHoc <- function(dataframe, which, level){
  if(level == "SPECIES"){
    if(which == "RES"){
      mod_multicomp <- glm.nb(max.OPG ~ intFacSPECIES, data = dataframe)
    } else if(which == "IMP"){
      mod_multicomp <- survreg(Surv(impact)~intFacSPECIES, data = dataframe, dist="weibull")
    } else if(which == "TOL"){
      mod_multicomp <- lm(ToleranceIndex ~ intFacSPECIES, data = dataframe)
    }
    return(summary(glht(mod_multicomp, linfct=mcp(intFacSPECIES = "Tukey"))))
  }
  if(level == "STRAINS"){
    if(which == "RES"){
      mod_multicomp <- glm.nb(max.OPG ~ intFacSTRAINS, data = dataframe)
    } else if(which == "IMP"){
      mod_multicomp <- survreg(Surv(impact)~intFacSTRAINS, data = dataframe, dist="weibull")
    } else if(which == "TOL"){
      mod_multicomp <- lm(ToleranceIndex ~ intFacSTRAINS, data = dataframe)
    }
    return(summary(glht(mod_multicomp, linfct=mcp(intFacSTRAINS = "Tukey"))))
  }
}

###############################
## Test factors significance ##
###############################

# Resistance
testSignif(art2al_SUMdf, "RES", "SPECIES")
testSignif(SUBsummaryDF77mice, "RES", "SPECIES")
testSignif(art2al_SUMdf, "RES", "STRAINS")
testSignif(SUBsummaryDF77mice, "RES", "STRAINS")
# Impact
testSignif(art2al_SUMdf, "IMP", "SPECIES")
testSignif(SUBsummaryDF77mice, "IMP", "SPECIES")
testSignif(art2al_SUMdf, "IMP", "STRAINS")
testSignif(SUBsummaryDF77mice, "IMP", "STRAINS")
# Tolerance
testSignif(art2al_SUMdf, "TOL", "SPECIES")
testSignif(SUBsummaryDF77mice, "TOL", "SPECIES")
testSignif(art2al_SUMdf, "TOL", "STRAINS")
testSignif(SUBsummaryDF77mice, "TOL", "STRAINS")

####################
## Post-hoc tests ##
####################

# to avoid running these long test all the time
doYouRun = "foncebebe"

if (doYouRun == "foncebebe"){
  # Resistance
  testPostHoc(art2al_SUMdf, "RES", "SPECIES")
  testPostHoc(SUBsummaryDF77mice, "RES", "SPECIES")
  testPostHoc(art2al_SUMdf, "RES", "STRAINS")
  testPostHoc(SUBsummaryDF77mice, "RES", "STRAINS")
  # Impact
  testPostHoc(art2al_SUMdf, "IMP", "SPECIES")
  testPostHoc(SUBsummaryDF77mice, "IMP", "SPECIES")
  testPostHoc(art2al_SUMdf, "IMP", "STRAINS")
  testPostHoc(SUBsummaryDF77mice, "IMP", "STRAINS")
}

## Translation of 1% because Weibull doesn't support nul data
coef(testSignif(art2al_SUMdf, "IMP", "SPECIES")$modfull)
coefImp <- exp(coef(testSignif(art2al_SUMdf, "IMP", "SPECIES")$modfull))
coefImp[1] - 0.01# Efer-MmD: 6.1%
coefImp[1] * coefImp[2] - 0.01 # Efal-MmD: 9.3%
coefImp[1] * coefImp[3] - 0.01# Efer-Mmm: 8.3%
coefImp[1] * coefImp[2] * coefImp[3] * coefImp[4] -0.01 # Efal-Mmm: 18.7%

if (doYouRun == "foncebebe"){
  # Tolerance
  testPostHoc(art2al_SUMdf, "TOL", "SPECIES")
  testPostHoc(SUBsummaryDF77mice, "TOL", "SPECIES")
  testPostHoc(art2al_SUMdf, "TOL", "STRAINS")
  testPostHoc(SUBsummaryDF77mice, "TOL", "STRAINS")
}

#################
## save output ##
#################
doYouSave = "foncebebe"
if (doYouSave == "foncebebe"){
  # Resistance
  write.csv(getMatrixPostHoc(testPostHoc(art2al_SUMdf, "RES", "SPECIES")),
            "../figures/posthocResSPECIES.csv")
  write.csv(getMatrixPostHoc(testPostHoc(SUBsummaryDF77mice, "RES", "SPECIES")),
            "../figures/posthocResSPECIES_77mice.csv")
  write.csv(getMatrixPostHoc(testPostHoc(art2al_SUMdf, "RES", "STRAINS")),
            "../figures/posthocResSTRAINS.csv")
  write.csv(getMatrixPostHoc(testPostHoc(SUBsummaryDF77mice, "RES", "STRAINS")),
            "../figures/posthocResSTRAINS_77mice.csv")
  # Impact
  write.csv(getMatrixPostHoc(testPostHoc(art2al_SUMdf, "IMP", "SPECIES")),
            "../figures/posthocImpSPECIES.csv")
  write.csv(getMatrixPostHoc(testPostHoc(SUBsummaryDF77mice, "IMP", "SPECIES")),
            "../figures/posthocImpSPECIES_77mice.csv")
  write.csv(getMatrixPostHoc(testPostHoc(art2al_SUMdf, "IMP", "STRAINS")),
            "../figures/posthocImpSTRAINS.csv")
  write.csv(getMatrixPostHoc(testPostHoc(SUBsummaryDF77mice, "IMP", "STRAINS")),
            "../figures/posthocImpSTRAINS_77mice.csv")
  # Tolerance
  write.csv(getMatrixPostHoc(testPostHoc(art2al_SUMdf, "TOL", "SPECIES")),
            "../figures/posthocTolSPECIES.csv")
  write.csv(getMatrixPostHoc(testPostHoc(SUBsummaryDF77mice, "TOL", "SPECIES")),
            "../figures/posthocTolSPECIES_77mice.csv")
  write.csv(getMatrixPostHoc(testPostHoc(art2al_SUMdf, "TOL", "STRAINS")),
            "../figures/posthocTolSTRAINS.csv")
  write.csv(getMatrixPostHoc(testPostHoc(SUBsummaryDF77mice, "TOL", "STRAINS")),
            "../figures/posthocTolSTRAINS_77mice.csv")
}

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
# plot marginal effects of interaction terms
posx.1 <- c(0.9,1.1, 1.9,2.1)
get_plotR_SPECIES <- function(dataframe){
  plot_model(testSignif(dataframe, "RES", "SPECIES")$modfull,
             type = "int", dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
    scale_color_manual(values = c("blue", "red"),
                       name = "Mouse subspecies",labels = c("Mmd", "Mmm")) +
    ggtitle("Maximum parasite load \n(mean and 95%CI)") +
    scale_y_continuous("(predicted) maximum oocysts per gram of feces (x10e6)", 
                       breaks = seq(0, 2500000, 500000),
                       labels = seq(0, 2500000, 500000)/1000000)+
    xlab("Eimeria species") +
    theme(axis.title.x = element_text(hjust=1), axis.text=element_text(size=13)) +
    geom_text(aes(x=posx.1,y=90000,label=getNs("max.OPG", dataframe, 
                                               "Mouse_subspecies", "Eimeria_species")), vjust=0)
}

# plot marginal effects of interaction terms by isolates & strains
posx.2 <- c(0.8+c(0,1/8,2/8,3/8),1.8+c(0,1/8,2/8,3/8),2.8+c(0,1/8,2/8,3/8))
get_plotR_STRAINS <- function(dataframe){
  plot_model(testSignif(dataframe, "RES", "STRAINS")$modfull,
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

plotR_SPECIES <- get_plotR_SPECIES(art2al_SUMdf)
plotR_SPECIES
plotR_SPECIES_77mice <- get_plotR_SPECIES(SUBsummaryDF77mice)
plotR_SPECIES_77mice

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

## NB. translate back 0.01
transValuesImp <- seq(0.01,0.31, 0.05)
as.character(transValuesImp - 0.01)
realValuesImpLabels <- c("0%", "5%", "10%", "15%", "20%", "25%", "30%")

get_plotI_SPECIES <- function(dataframe){
  plot_model(testSignif(dataframe, "IMP", "SPECIES")$modfull,
             type = "int",dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
    scale_color_manual(values = c("blue","red"),
                       name = "Mouse subspecies",labels = c("Mmd", "Mmm")) +
    xlab("Eimeria species") +
    ggtitle("Maximum weight loss \n(mean and 95%CI)") +
    scale_y_continuous(breaks = transValuesImp, labels = realValuesImpLabels, 
                       name = "(predicted) maximum weight loss compared to day of infection")+
    theme(axis.title.x = element_text(hjust=1), axis.text=element_text(size=13)) +
    geom_text(aes(x=posx.1,y=0,label=getNs("relWL", dataframe,
                                           "Mouse_subspecies", "Eimeria_species")),vjust=0)
}

plotI_SPECIES <- get_plotI_SPECIES(art2al_SUMdf)
plotI_SPECIES
plotI_SPECIES_77mice <- get_plotI_SPECIES(SUBsummaryDF77mice)
plotI_SPECIES_77mice

get_plotI_STRAINS <- function(dataframe){
  plot_model(testSignif(dataframe, "IMP", "STRAINS")$modfull,
             type = "int",dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
    scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                       name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +
    xlab("Eimeria isolate") +
    ggtitle("Maximum weight loss \n(mean and 95%CI)") +
    scale_y_continuous(breaks = transValuesImp, labels = realValuesImpLabels, 
                       name = "(predicted) maximum weight loss compared to day of infection")+
    theme(axis.title.x = element_text(hjust=1), axis.text=element_text(size=13)) +
    geom_text(aes(x=posx.2,y=0,label=getNs("relWL", dataframe)),vjust=0)
}

plotI_STRAINS  <- get_plotI_STRAINS(art2al_SUMdf)
plotI_STRAINS
plotI_STRAINS_77mice <- get_plotI_STRAINS(SUBsummaryDF77mice)

###############
## Tolerance ##

get_plotT_SPECIES <- function(dataframe){
  plot_model(testSignif(dataframe, "TOL", "SPECIES")$modfull, 
             type = "int", dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
    scale_color_manual(values = c("blue", "red"),
                       name = "Mouse subspecies",labels = c("Mmd", "Mmm")) +
    xlab("Eimeria species") +
    ggtitle("Tolerance index \n(mean and 95%CI)") +
    scale_y_continuous("(predicted) tolerance index")+
    theme(axis.title.x = element_text(hjust=1), axis.text = element_text(size=13))+
    geom_text(aes(x=posx.1,y=0.7,label=getNs("ToleranceIndex", dataframe,
                                             "Mouse_subspecies", "Eimeria_species")),vjust=0)
}
plotT_SPECIES <- get_plotT_SPECIES(art2al_SUMdf)
plotT_SPECIES
plotT_SPECIES_77mice <- get_plotT_SPECIES(SUBsummaryDF77mice)

get_plotT_STRAINS <- function(dataframe){
  plot_model(testSignif(dataframe, "TOL", "STRAINS")$modfull,
             type = "int", dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
    scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                       name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +
    xlab("Eimeria isolate") +
    ggtitle("Tolerance index \n(mean and 95%CI)") +
    scale_y_continuous("(predicted) tolerance index")+
    theme(axis.title.x = element_text(hjust=1), axis.text = element_text(size=13))+
    geom_text(aes(x=posx.2,y=0.7,label=getNs("ToleranceIndex", dataframe)),vjust=0)
}
plotT_STRAINS <- get_plotT_STRAINS(art2al_SUMdf)
plotT_STRAINS
plotT_STRAINS_77mice <- get_plotT_STRAINS(SUBsummaryDF77mice)

# Fig 3.
Fig3 <- cowplot::plot_grid(plotR_SPECIES + theme(legend.position = "none"),
                           plotI_SPECIES + theme(legend.position = "none"),
                           plotT_SPECIES + theme(legend.position = "none"), 
                           plotT_SPECIES,
                           labels=c("A", "B", "C", "D"), label_size = 20)
  
Fig3
pdf(file = "../figures/Fig3.pdf",
    width = 9, height = 9)
Fig3
dev.off()

## Fig 4
Fig4 <-  cowplot::plot_grid(
  plotR_STRAINS + theme(legend.position = "none"),
  plotI_STRAINS + theme(legend.position = "none"),
  plotT_STRAINS+ theme(legend.position = "none"), 
  plotT_STRAINS,
  labels=c("A", "B", "C", "D"), label_size = 20)
Fig4

pdf(file = "../figures/Fig4.pdf",
    width = 9, height = 9)
Fig4
dev.off()

### SUB df
pdf(file = "../figures/FigSPECIES_77mice.pdf",
    width = 9, height = 9)
cowplot::plot_grid(
  plotR_SPECIES_77mice + theme(legend.position = "none"),
  plotI_SPECIES_77mice + theme(legend.position = "none"),
  plotT_SPECIES_77mice+ theme(legend.position = "none"), 
  plotT_SPECIES_77mice,
  labels=c("A", "B", "C", "D"), label_size = 20)
dev.off()

pdf(file = "../figures/FigSTRAINS_77mice.pdf",
    width = 9, height = 9)
cowplot::plot_grid(
  plotR_STRAINS_77mice + theme(legend.position = "none"),
  plotI_STRAINS_77mice + theme(legend.position = "none"),
  plotT_STRAINS_77mice+ theme(legend.position = "none"), 
  plotT_STRAINS_77mice,
  labels=c("A", "B", "C", "D"), label_size = 20)
dev.off()

### Second part: correlation resistance / tolerance

# take the predictions from before 
mydatx <- ggeffects::ggpredict(
  model = testSignif(art2al_SUMdf, "RES", "SPECIES")$modfull, 
  terms = c("Eimeria_species", "Mouse_subspecies"), ci.lvl = 0.95)
names(mydatx)[2:5] <- paste0(names(mydatx)[2:5], "_OPG")

mydaty <- ggeffects::ggpredict(
  model = testSignif(art2al_SUMdf, "IMP", "SPECIES")$modfull, 
  terms = c("Eimeria_species", "Mouse_subspecies"), ci.lvl = 0.95)
names(mydaty)[2:5] <- paste0(names(mydaty)[2:5], "_WL")

## NB. translate back 0.01 the y axis because we did 
# art2al_SUMdf$impact <- art2al_SUMdf$relWL + 0.01
mydaty$predicted_WL <- mydaty$predicted_WL - 0.01
mydaty$conf.low_WL <- mydaty$conf.low_WL - 0.01
mydaty$conf.high_WL <- mydaty$conf.high_WL - 0.01

mydat <- merge(data.frame(mydatx), data.frame(mydaty))

# Trade-off plot
Fig5_couplingplot <- ggplot(mydat, aes(x = predicted_OPG, y = predicted_WL, col = group)) +
  geom_errorbar(aes(ymin = conf.low_WL, ymax = conf.high_WL))+
  geom_errorbarh(aes(xmin = conf.low_OPG, xmax = conf.high_OPG))+
  geom_smooth(method = "lm", col = "black", alpha = .2, aes(linetype = x)) + 
  theme_bw()+
  scale_color_manual(values = c("blue", "red")) +
  scale_x_log10(name = "(predicted) maximum oocysts per gram of feces")+
  scale_y_continuous(name = "(predicted) maximum weight loss \ncompared to day of infection",
                     breaks = seq(0,0.3, 0.05), 
                     labels = c("0%", "5%", "10%", "15%", "20%", "25%", "30%")) +
  ggtitle("Maximum weight loss = f(maximum parasite load) \n(mean and 95%CI)")

pdf("../figures/Fig5_couplingplot.pdf", width = 6, height = 5)
Fig5_couplingplot  
dev.off()

# https://stats.idre.ucla.edu/r/seminars/interactions-r/
# Test the difference between slopes (simplified):
art2al_SUMdf$group <- paste0(art2al_SUMdf$Eimeria_species, "_", art2al_SUMdf$Mouse_subspecies)

modO_W <- lm(formula = relWL ~ max.OPG * Eimeria_species, data = art2al_SUMdf)
modO_W
summary(modO_W)
#p-value of the t-statistic for the interaction between ResistanceIndex and Eimeria_species: 1.39e-05 ***

# The interaction term is significant, which suggests
# that the relationship of relWL by maxOPG does vary by Eimeria spieces
modO_WA <- lm(formula = relWL ~ max.OPG * Eimeria_species, data = art2al_SUMdf)
modO_WB <- lm(formula = relWL ~ max.OPG + Eimeria_species, data = art2al_SUMdf)
modO_WC <- lm(formula = relWL ~ max.OPG, data = art2al_SUMdf)
modO_WD <- lm(formula = relWL ~ Eimeria_species, data = art2al_SUMdf)


olB <- lm(formula = ToleranceIndex ~ ResistanceIndex + Eimeria_species, data = summaryDF108mice)
modResTolC <- lm(formula = ToleranceIndex ~ ResistanceIndex, data = summaryDF108mice)
modResTolD <- lm(formula = ToleranceIndex ~ Eimeria_species, data = summaryDF108mice)

# signif interaction:
homemadeGtest(modResTolA, modResTolB)
# signif eimeria:
homemadeGtest(modResTolA, modResTolC)
# signif resistance index:
homemadeGtest(modResTolA, modResTolD)

# Since our goal is to obtain simple slopes of parasite:
library(lsmeans)
emtrends(modResTolA, ~ Eimeria_species, var="ResistanceIndex")

# The 95% confidence interval does not contain zero for E. falciformis but contains zero
# for E. ferrisi, so the simple slope is significant for E. falciformis but not for E. ferrisi

########### Remove the outlier:

summaryDF_rmoutlier <- summaryDF108mice[!summaryDF108mice$ResistanceIndex < 0.25 & !is.na(summaryDF108mice$ResistanceIndex),]

modResTolA2 <- lm(formula = ToleranceIndex ~ ResistanceIndex * Eimeria_species, data = summaryDF_rmoutlier)
modResTolB2 <- lm(formula = ToleranceIndex ~ ResistanceIndex + Eimeria_species, data = summaryDF_rmoutlier)
modResTolC2 <- lm(formula = ToleranceIndex ~ ResistanceIndex, data = summaryDF_rmoutlier)
modResTolD2 <- lm(formula = ToleranceIndex ~ Eimeria_species, data = summaryDF_rmoutlier)

# signif interaction:
homemadeGtest(modResTolA2, modResTolB2)
# signif eimeria:
homemadeGtest(modResTolA2, modResTolC2)
# signif resistance index:
homemadeGtest(modResTolA2, modResTolD2)

# Since our goal is to obtain simple slopes of parasite:
emtrends(modResTolA2, ~ Eimeria_species, var="ResistanceIndex")


# for supp data
SUBsummaryDF77mice$group <- paste(SUBsummaryDF77mice$Mouse_subspecies, SUBsummaryDF77mice$Eimeria_species, sep = "_")
gd2 <- SUBsummaryDF77mice %>%
  group_by(group, Mouse_subspecies, Eimeria_species) %>% 
  summarise(
    Wmean = mean(relWL, na.rm = T),
    Wsd = sd(relWL, na.rm = T),
    Wn = sum(!is.na(relWL)),
    Omean = mean(max.OPG*10e-6, na.rm = T),
    Osd = sd(max.OPG*10e-6, na.rm = T),
    On = sum(!is.na(max.OPG))) %>%
  mutate(Wse = Wsd / sqrt(Wn),
         Wlower.ci = Wmean - qt(1 - (0.05 / 2), Wn - 1) * Wse,
         Wupper.ci = Wmean + qt(1 - (0.05 / 2), Wn - 1) * Wse,
         Ose = Osd / sqrt(On),
         Olower.ci = Omean - qt(1 - (0.05 / 2), On - 1) * Ose,
         Oupper.ci = Omean + qt(1 - (0.05 / 2), On - 1) * Ose)

Fig5_couplingplot <- ggplot(gd2, aes(x = Omean, y = Wmean, col = Mouse_subspecies)) +
  geom_errorbar(aes(ymin = Wlower.ci, ymax = Wupper.ci))+
  geom_errorbarh(aes(xmin = Olower.ci, xmax = Oupper.ci))+
  geom_smooth(method = "lm", col = "black", alpha = .2, aes(linetype = Eimeria_species)) + 
  theme_bw()+
  scale_color_manual(values = c("blue", "red")) +
  ylab(label = "Minimum relative weight compared to day 0") +
  scale_x_continuous(name = "maximum oocysts per gram of feces (x10^6)")+
  scale_y_continuous(breaks = seq(0,0.20, 0.05), labels = c("0%", "5%", "10%", "15%", "20%"), 
                     name = "(predicted) maximum weight loss compared to day of infection") +
  ggtitle("Maximum weight loss = f(maximum parasite load) \n(mean and 95%CI)")

Fig5_couplingplot

pdf("../figures/suppleResTolPlot77mice.pdf", width = 12, height = 9)
restolplot77  
dev.off()



