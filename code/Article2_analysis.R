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
  }
  return(list(modfull = modFULL, LRT = myLRTsignificanceFactors(modFULL, modPara, modMous, modinter)))
}

testPostHoc <- function(dataframe, which, level){
  if(level == "SPECIES"){
    if(which == "RES"){
      mod_multicomp <- glm.nb(max.OPG ~ intFacSPECIES, data = dataframe)
    } else if(which == "IMP"){
      mod_multicomp <- survreg(Surv(impact)~intFacSPECIES, data = dataframe, dist="weibull")
    } 
    return(summary(glht(mod_multicomp, linfct=mcp(intFacSPECIES = "Tukey"))))
  }
  if(level == "STRAINS"){
    if(which == "RES"){
      mod_multicomp <- glm.nb(max.OPG ~ intFacSTRAINS, data = dataframe)
    } else if(which == "IMP"){
      mod_multicomp <- survreg(Surv(impact)~intFacSTRAINS, data = dataframe, dist="weibull")
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
testSignif(SUBsummaryDF77mice, "IMP", "STRAINS") ####### !!!
## Translation of 1% because Weibull doesn't support nul data
coef(testSignif(art2al_SUMdf, "IMP", "SPECIES")$modfull)
coefImp <- exp(coef(testSignif(art2al_SUMdf, "IMP", "SPECIES")$modfull))
coefImp[1] - 0.01# Efer-MmD: 6.1%
coefImp[1] * coefImp[2] - 0.01 # Efal-MmD: 9.3%
coefImp[1] * coefImp[3] - 0.01# Efer-Mmm: 8.3%
coefImp[1] * coefImp[2] * coefImp[3] * coefImp[4] -0.01 # Efal-Mmm: 18.7%

####################
## Post-hoc tests ##
####################

# to avoid running these long test all the time
doYouRun = "keepit"

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
plotI_STRAINS_77mice

# # Fig 3.
# Fig3 <- cowplot::plot_grid(plotR_SPECIES + theme(legend.position = "none"),
#                            plotI_SPECIES + theme(legend.position = "none"),
#                            plotR_SPECIES,
#                            labels=c("A", "B", "C"), label_size = 20)
#   
# Fig3
# pdf(file = "../figures/Fig3.pdf",
#     width = 9, height = 9)
# Fig3
# dev.off()

## Fig 4
Fig4 <-  cowplot::plot_grid(
  plotR_STRAINS + theme(legend.position = "none"),
  plotI_STRAINS + theme(legend.position = "none"),
  plotR_STRAINS,
  labels=c("A", "B", "C"), label_size = 20)
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
  plotR_SPECIES_77mice,
  labels=c("A", "B", "C"), label_size = 20)
dev.off()

pdf(file = "../figures/FigSTRAINS_77mice.pdf",
    width = 9, height = 9)
cowplot::plot_grid(
  plotR_STRAINS_77mice + theme(legend.position = "none"),
  plotI_STRAINS_77mice + theme(legend.position = "none"),
  plotR_STRAINS_77mice,
  labels=c("A", "B", "C"), label_size = 20)
dev.off()

########################################
### Second part: assessing tolerance ###
########################################

art2al_SUMdf$group <- factor(paste(art2al_SUMdf$Mouse_genotype, art2al_SUMdf$infection_isolate, sep = "_"))

modfull <- lm(relWL ~ max.OPG * infection_isolate * Mouse_genotype, data = art2al_SUMdf)
modNOmouse <- lm(relWL ~ max.OPG * infection_isolate, data = art2al_SUMdf)
modNOisolate <- lm(relWL ~ max.OPG * Mouse_genotype, data = art2al_SUMdf)
modNOint <- lm(relWL ~ max.OPG + infection_isolate + Mouse_genotype, data = art2al_SUMdf)

homemadeGtest(modfull, modNOisolate)
homemadeGtest(modfull, modNOmouse)
homemadeGtest(modfull, modNOint)

# Calculate slopes for each group:
df.allgroups <- data.frame(infection_isolate = rep(levels(art2al_SUMdf$infection_isolate),4),
                           Mouse_genotype = c(rep(levels(art2al_SUMdf$Mouse_genotype)[1], 3),
                                              rep(levels(art2al_SUMdf$Mouse_genotype)[2], 3),
                                              rep(levels(art2al_SUMdf$Mouse_genotype)[3], 3),
                                              rep(levels(art2al_SUMdf$Mouse_genotype)[4], 3)),
                           group = NA, slope = NA)
df.allgroups$group <- paste(df.allgroups$Mouse_genotype, df.allgroups$infection_isolate, sep = "_")

for (i in 1:length(df.allgroups$group)){
  mod.subgroup <- lm(relWL ~ 0 + max.OPG, 
                     data = art2al_SUMdf[art2al_SUMdf$group %in% df.allgroups$group[i],])
  df.allgroups$slope[i] <- coef(mod.subgroup)
  df.allgroups$down[i] <- confint(mod.subgroup)[1]
  df.allgroups$up[i] <- confint(mod.subgroup)[2]
  
  }

df.allgroups$origin.y <- 0
df.allgroups <- melt(df.allgroups, measure.vars = c("slope", "origin.y"))
df.allgroups$down[df.allgroups$variable %in% "origin.y"] <- 0
df.allgroups$up[df.allgroups$variable %in% "origin.y"] <- 0
df.allgroups$max.OPG <- as.numeric(as.character(
  plyr::mapvalues(df.allgroups$variable, 
                  from = c("slope", "origin.y"), to = c(1, 0))))
df.allgroups$max.OPG <- df.allgroups$max.OPG*4e6
df.allgroups$relWL <- df.allgroups$value * 4e6

T1 <- ggplot(data = df.allgroups, aes(x = max.OPG, y = relWL, group = group, 
                                      col = Mouse_genotype)) +
  geom_smooth(method = "lm", se = F) +
  geom_point(data = art2al_SUMdf, size = 4, pch = 3)+
  facet_grid(.~infection_isolate) +
  scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                     name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +  scale_x_continuous(name = "maximum oocysts per gram of feces") +
  scale_y_continuous(name = "maximum weight loss compared to day of infection",
                     breaks = seq(0,0.3, 0.05), 
                     labels = c("0%", "5%", "10%", "15%", "20%", "25%", "30%"))+
  coord_cartesian(ylim=c(0, 0.3)) +
  theme(legend.position = "top")
T1
  
## Means & 95% CI of each slope
T2 <- ggplot(data = df.allgroups[df.allgroups$variable %in% "slope",],
       aes(x = infection_isolate, y = value, col = Mouse_genotype)) +
  geom_point(size = 4, position = position_dodge(width = .5))+
  geom_errorbar(aes(ymin = down, ymax = up), width = .1,
                position = position_dodge(width = .5))+
  # facet_grid(.~infection_isolate) +
  scale_y_continuous(name = "slope of relative weight loss on OPG (inverse of tolerance)")+
  coord_cartesian(ylim=c(-0.3e-7,4.5e-7)) +
  scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                     name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +
  theme(legend.position = "top")
  
FigTOL <- cowplot::plot_grid(T1, T2,
                             labels=c("A", "B"), label_size = 20)

FigTOL

########## Coupling res tol
d <- df.allgroups[df.allgroups$variable %in% "slope",
                  c("infection_isolate", "Mouse_genotype", "down", "up", "value")]
d

# take the predictions from before 
mydatx <- ggeffects::ggpredict(
  model = testSignif(art2al_SUMdf, "RES", "STRAINS")$modfull, 
  terms = c("infection_isolate", "Mouse_genotype"), ci.lvl = 0.95)
names(mydatx)[2:5] <- paste0(names(mydatx)[2:5], "_OPG")
mydatx <- data.frame(mydatx)
names(mydatx)[names(mydatx)%in% "x"] <- "infection_isolate"
names(mydatx)[names(mydatx)%in% "group"] <- "Mouse_genotype"

d <- merge(d, mydatx)

figCoupl <- ggplot(d, aes(x = predicted_OPG, y = value,col = Mouse_genotype)) +
  geom_point(size = 4)+
  geom_errorbar(aes(ymin = down, ymax = up), width = .1) +
  geom_errorbarh(aes(xmin = conf.low_OPG, xmax = conf.high_OPG), width = .1,
                position = position_dodge(width = .5)) +
  coord_cartesian(ylim=c(-0.3e-7,4.5e-7))+
  scale_y_continuous(name = "slope of relative weight loss on OPG (inverse of tolerance)")+
  scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                     name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +
  facet_grid(.~infection_isolate)
figCoupl

## NB supplementary.