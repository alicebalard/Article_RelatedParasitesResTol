### Code for data analysis of Article 2
### August 2019
### Alice Balard

## INFO
# Mouse AA_0088, HI = 0.2
# Mouse AA_0064, HI = 0.08
# Mouse AA_0139, HI = 0.85

#### Load data and functions within ####
source("dataPreparation.R")

######################################
########## Read information ##########
######################################

DF_all # all 168 mice
length(na.omit(ALL_summary_F0$relWL))# 108 F0 mice
length(na.omit(ALL_summary_F1$relWL))# 60 F1 mice
length(na.omit(ALL_summary_F0$peak.oocysts.per.g.mouse))# 99 F0 mice
length(na.omit(ALL_summary_F1$peak.oocysts.per.g.mouse))# 59 F1 mice
length(na.omit(ALL_summary_F0$invtolerance))# 103 F0 mice
length(na.omit(ALL_summary_F1$invtolerance))# 60 F0 mice

ALL_summary$Exp_2groups <- "group2"
ALL_summary$Exp_2groups[ALL_summary$Exp_ID %in% c("Exp_003", "Exp_004")] <- "group1"
ALL_summary[ALL_summary$Exp_2groups %in% "group1",] %>%
  dplyr::select(Mouse_genotype, infection_isolate) %>%
  dplyr::group_by(Mouse_genotype, infection_isolate)%>%
  dplyr::summarize(sum=n())%>%
  tidyr::spread(infection_isolate, value = sum)
ALL_summary[ALL_summary$Exp_2groups %in% "group2",] %>%
  dplyr::select(Mouse_genotype, infection_isolate) %>%
  dplyr::group_by(Mouse_genotype, infection_isolate)%>%
  dplyr::summarize(sum=n())%>%
  tidyr::spread(infection_isolate, value = sum)

###### what is the overall peak day for each parasite isolate? ######
aggregate(ALL_summary$dpi_max.oocysts.per.tube,
          list(ALL_summary$infection_isolate), function(x) {paste(length(x), median(x), round(sd(x),2))})
aggregate(ALL_summary$dpi_minWeight,
          list(ALL_summary$infection_isolate), function(x) {paste(length(x), median(x), round(sd(x),2))})

###### what is the overall prepatent period for each parasite isolate? ######
d <- as.data.frame(
  DF_all[!is.na(DF_all$oocysts.per.tube) & DF_all$oocysts.per.tube > 0,] %>% 
    dplyr::group_by(EH_ID) %>%
    dplyr::slice(which.min(dpi)) %>%
    dplyr::select(EH_ID, weight, HI, startingWeight, ageAtInfection, Sex,
                  Mouse_genotype, Eimeria_species, Mouse_subspecies,
                  infection_isolate, Exp_ID, dpi))
aggregate(d$dpi,
          list(d$infection_isolate), function(x) {paste(length(x), median(x), round(sd(x),2))})

ggplot(DF_all, aes(dpi, oocysts.per.tube, group = EH_ID)) + geom_line() + facet_grid(.~infection_isolate)

#################################################
########## Chose correct distributions ##########
#################################################

xRes <- as.numeric(na.omit(ALL_summary$peak.oocysts.per.g.mouse))
xImp <- as.numeric(na.omit(ALL_summary$relWL))
xTol <- as.numeric(na.omit(ALL_summary$invtolerance))

hist(xRes, breaks = 100)
descdist(xRes)
findGoodDist(x = xRes, distribs = c("normal", "negative binomial", "poisson"), 
             distribs2 = c("norm", "nbinom", "pois"))
### nbinom for resistance

hist(xImp, breaks = 100)
descdist(xImp)
findGoodDist(x = xImp+ 0.01, distribs = c("normal", "beta", "gamma", "weibull"), 
             distribs2 = c("norm", "beta", "gamma", "weibull"))
### weibull for impact on health
ALL_summary$impact <- ALL_summary$relWL +0.01
ALL_summary_F0$impact <- ALL_summary_F0$relWL +0.01
ALL_summary_F1$impact <- ALL_summary_F1$relWL +0.01

# for tolerance, we need to find a way to deal with the massive extreme values
hist(xTol, breaks = 1000)
descdist(xTol)
## too skewed, so we log transform it
translation <- 1e-8
xTol2 <- log10(xTol+translation)
range(xTol2)
hist(xTol2, breaks = 1000)
descdist(xTol2)
findGoodDist(x = xTol2, distribs = c("normal"), 
             distribs2 = c("norm"))

### We translate, log transfor the data to deal with skewness
ALL_summary$tolerance <- log10(ALL_summary$invtolerance+translation) + 8
ALL_summary_F0$tolerance <- log10(ALL_summary_F0$invtolerance+translation) + 8
ALL_summary_F1$tolerance <- log10(ALL_summary_F1$invtolerance+translation) + 8

ALL_summary[is.na(ALL_summary$tolerance), c("EH_ID", "Mouse_genotype", "relWL", "peak.oocysts.per.g.mouse")]

##################################
########## F0 mice glms ##########
##################################

################
## Resistance ##
################

# Resistance on F0 (parental strains)
modRes <- glm.nb(peak.oocysts.per.g.mouse ~ infection_isolate*Mouse_genotype, data = ALL_summary_F0)
anova(modRes, test = "LRT")
# SIGNIF infection isolate (p-value = 0.01902) + interactions with mice (p-value = 8.432e-05)
plot_model(modRes, type = "int", dot.size = 4, dodge = .5)

modRes2 <- glm.nb(peak.oocysts.per.g.mouse ~ Eimeria_species*Mouse_subspecies, data = ALL_summary_F0)
anova(modRes2)
# SIGNIF infection isolate (p-value = 0.02236) + interactions with mice (p-value = 6.52e-07)
plot_model(modRes2, type = "int", dot.size = 4, dodge = .5)

############
## Impact ##
############

## Translation of 1% because Weibull doesn't support nul data
modImp <- survreg(Surv(impact)~infection_isolate*Mouse_genotype, data = ALL_summary_F0, dist="weibull")
anova(modImp)
length(ALL_summary_F0$relWL)
# Eimeria isolate significant

## Plot
plot_model(modImp, type = "int",dot.size = 4, dodge = .5)

modImp2 <- survreg(Surv(impact)~Eimeria_species*Mouse_subspecies, data = ALL_summary_F0, dist="weibull")
anova(modImp2) # Eimeria species AND mouse subspecies significant
plot_model(modImp2, type = "int",dot.size = 4, dodge = .5) 

###############
## Tolerance ##
###############

length(na.omit(ALL_summary_F0$tolerance))
modTol <- lm(tolerance ~ infection_isolate*Mouse_genotype, data = ALL_summary_F0)
anova(modTol)
# mouse genotype & interactions significant

## Plot
plot_model(modTol, type = "int", dot.size = 4, dodge = .5)

modTol2 <- lm(tolerance ~ Eimeria_species*Mouse_subspecies, data = ALL_summary_F0)
anova(modTol2)
plot_model(modTol2, type = "int", dot.size = 4, dodge = .5)

################################
########## F1 mice ML ##########
################################

ALL_summary_F1$hybridLevel <- "hybrid"
ALL_summary_F1$hybridLevel[ALL_summary_F1$HI<0.2] <- "outbred domesticus"
ALL_summary_F1$hybridLevel[ALL_summary_F1$HI>0.8] <- "outbred musculus" 
ALL_summary_F1$hybridLevel <- factor(ALL_summary_F1$hybridLevel, levels = c("outbred domesticus", "hybrid", "outbred musculus"))

################
## Resistance ##
################

ALL_summary_F1$Eimeria_species <- as.factor(ALL_summary_F1$Eimeria_species)

# to fix the scaling issue
hist(ALL_summary_F1$peak.oocysts.per.g.mouse, breaks = 100)
hist(round(ALL_summary_F1$peak.oocysts.per.g.mouse/1000), breaks = 100)

nrow(ALL_summary_F1[!is.na(ALL_summary_F1$peak.oocysts.per.g.mouse),])
nrow(ALL_summary_F1[!is.na(ALL_summary_F1$peak.oocysts.per.g.mouse) & 
                      ALL_summary_F1$Eimeria_species %in% "E.falciformis",])
nrow(ALL_summary_F1[!is.na(ALL_summary_F1$peak.oocysts.per.g.mouse) & 
                      ALL_summary_F1$Eimeria_species %in% "E.ferrisi",])

ALL_summary_F1$peak.oocysts.per.g.mouse.div100 <- round(ALL_summary_F1$peak.oocysts.per.g.mouse/100)

## 1. discrete lm : test factors P & H & interactions
modResF1 <- glm.nb(peak.oocysts.per.g.mouse ~ hybridLevel*Eimeria_species, data = ALL_summary_F1)
anova(modResF1, test = "LRT") # eimeria, host, inter signif

plot_model(modResF1, type = "int", dot.size = 4, dodge = .5,
                      colors = c("orange", "darkgreen"), title = "Resistance")

## 2. continuous ML : test HV
fitRes <- parasiteLoad::analyse(data = ALL_summary_F1,
                                response = "peak.oocysts.per.g.mouse.div100", 
                                group = "Eimeria_species", model = "negbin")

# [1] "Testing H0 no alpha vs alpha"
# dLL dDF     pvalue
# 1 2.58   1 0.02318184
# [1] "Testing H1 no alpha vs alpha"
# dLL dDF    pvalue
# 1 2.69   1 0.0202649
# [1] "Testing H2 groupA no alpha vs alpha"
# dLL dDF     pvalue
# 1 1.92   1 0.04995505
# [1] "Testing H2 groupB no alpha vs alpha"
# dLL dDF    pvalue
# 1 1.29   1 0.1082959
# [1] "Testing H3 groupA no alpha vs alpha"
# dLL dDF    pvalue
# 1 1.53   1 0.0801773
# [1] "Testing H3 groupB no alpha vs alpha"
# dLL dDF   pvalue
# 1 1.15   1 0.128581
# [1] "Testing H1 vs H0"
# dLL dDF    pvalue
# 1 0.57   2 0.5683046
# [1] "Testing H2 vs H0"
# dLL dDF     pvalue
# 1 5.21   4 0.03405811
# [1] "Testing H3 vs H1"
# dLL dDF     pvalue
# 1 8.67   6 0.00810378
# [1] "Testing H3 vs H2"
# dLL dDF     pvalue
# 1 4.03   4 0.08931408

coef(fitRes$H1)
# L1          L2          A1          A2       alpha           Z 
# 746.7278613 629.2621833   0.6075115   0.9480888   0.7319809  -0.8364519 

# to plot on axis
scaleOO <- c(50, 100, 200, 500, 1000, 3000)

parasiteLoad::bananaPlot(mod = fitRes$H1,
                                        data = ALL_summary_F1,
                                        response = "peak.oocysts.per.g.mouse.div100",
                                        islog10 = T, group = "Eimeria_species",
                                        cols = c("orange", "darkgreen")) + theme_bw() +
  ylab("Maximum parasite load (oocysts) per mouse gram") +
  xlab("Hybrid index") +
  scale_y_log10(breaks=scaleOO, labels = format(scaleOO*100, scientific = TRUE))+
  theme(axis.title.x = element_text(hjust=1), axis.text = element_text(size=13))

coef(fitRes$H3$groupA)
# L1          L2          A1          A2       alpha           Z 
# 888.4995253 235.4619325   0.8652611   0.6782603   0.8651792  -1.0900823 

coef(fitRes$H3$groupB)
# L1          L2          A1          A2       alpha           Z 
# 642.9575588 891.7035835   0.2169692   0.6284915   0.5955537  -0.1252396 

parasiteLoad::bananaPlot(mod = fitRes$H3,
                                        data = ALL_summary_F1,
                                        response = "peak.oocysts.per.g.mouse.div100",
                                        islog10 = T, group = "Eimeria_species",
                                        cols = c("orange", "darkgreen"))

######################
## Impact on health ##
######################

nrow(ALL_summary_F1[!is.na(ALL_summary_F1$impact),])
nrow(ALL_summary_F1[!is.na(ALL_summary_F1$impact) & 
                      ALL_summary_F1$Eimeria_species %in% "E.falciformis",])
nrow(ALL_summary_F1[!is.na(ALL_summary_F1$impact) & 
                      ALL_summary_F1$Eimeria_species %in% "E.ferrisi",])

## 1. discrete lm : test factors P & H & interactions
modImpF1 <- survreg(Surv(impact)~ hybridLevel*Eimeria_species, data = ALL_summary_F1, dist="weibull")
anova(modImpF1) # Eimeria signif. 0.022
plot_model(modImpF1, type = "int", dot.size = 4, dodge = .5, 
           colors = c("orange", "darkgreen"), title = "Impact on host health") 

## 2. discrete lm : test HV
modImpF1 <- survreg(Surv(impact)~ hybridLevel*Eimeria_species, data = ALL_summary_F1, dist="weibull")

modResF1_fal <- survreg(Surv(impact)~ hybridLevel*Eimeria_species,  dist="weibull",
                        data = ALL_summary_F1[ALL_summary_F1$Eimeria_species %in% "E.falciformis",])
anova(modResF1_fal) # no signif

modResF1_fer <- survreg(Surv(impact)~ hybridLevel*Eimeria_species,  dist="weibull",
                        data = ALL_summary_F1[ALL_summary_F1$Eimeria_species %in% "E.ferrisi",])
anova(modResF1_fer) # no signif

## 3. continuous ML : test HV
fitImp <- parasiteLoad::analyse(data = ALL_summary_F1,
                                response = "impact", 
                                group = "Eimeria_species", model = "weibull")

# [1] "Testing H0 no alpha vs alpha"
# dLL dDF     pvalue
# 1 1.94   1 0.04914124
# [1] "Testing H1 no alpha vs alpha"
# dLL dDF     pvalue
# 1 1.79   1 0.05840147
# [1] "Testing H2 groupA no alpha vs alpha"
# dLL dDF     pvalue
# 1 1.58   1 0.07533839
# [1] "Testing H2 groupB no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.26   1 0.4702023
# [1] "Testing H3 groupA no alpha vs alpha"
# dLL dDF    pvalue
# 1 1.62   1 0.0720455
# [1] "Testing H3 groupB no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.12   1 0.6177118
# [1] "Testing H1 vs H0"
# dLL dDF    pvalue
# 1 0.66   1 0.2509531
# [1] "Testing H2 vs H0"
# dLL dDF  pvalue
# 1 2.63   3 0.15408
# [1] "Testing H3 vs H1"
# dLL dDF    pvalue
# 1 3.14   4 0.1789046
# [1] "Testing H3 vs H2"
# dLL dDF    pvalue
# 1 1.17   2 0.3091788

coef(fitImp$H1)
# L1         L2      alpha    myshape 
# 0.07698095 0.10305627 0.58406414 1.46473120 

parasiteLoad::bananaPlot(mod = fitImp$H1,
                                        data = ALL_summary_F1,
                                        response = "impact",
                                        islog10 = F, group = "Eimeria_species",
                                        cols = c("orange", "darkgreen")) + theme_bw() +
  ylab("Impact") +
  xlab("Hybrid index") +
  theme(axis.title.x = element_text(hjust=1), axis.text = element_text(size=13))

coef(fitImp$H3$groupA)
# L1         L2      alpha    myshape 
# 0.09544317 0.13216564 0.73310201 1.54840752 

coef(fitImp$H3$groupB)
# L1         L2      alpha    myshape 
# 0.05095116 0.07761396 0.22188069 1.55068685 

parasiteLoad::bananaPlot(mod = fitImp$H3,
                                        data = ALL_summary_F1,
                                        response = "impact",
                                        islog10 = F, group = "Eimeria_species",
                                        cols = c("orange", "darkgreen")) 
###############
## Tolerance ##
###############
nrow(ALL_summary_F1[!is.na(ALL_summary_F1$tolerance),])
nrow(ALL_summary_F1[!is.na(ALL_summary_F1$tolerance) & 
                      ALL_summary_F1$Eimeria_species %in% "E.falciformis",])
nrow(ALL_summary_F1[!is.na(ALL_summary_F1$tolerance) & 
                      ALL_summary_F1$Eimeria_species %in% "E.ferrisi",])

## 1. discrete lm : test factors P & H & interactions
modTolF1 <- lm(tolerance ~  hybridLevel*Eimeria_species, data = ALL_summary_F1)
anova(modTolF1) # eimeria signif 0.039
plot_model(modTolF1, type = "int", dot.size = 4, dodge = .5, 
           colors = c("orange", "darkgreen"), title = "Tolerance") 

## 2. discrete lm : test HV
modTolF1_fal <- lm(tolerance ~  hybridLevel,
                   data = ALL_summary_F1[ALL_summary_F1$Eimeria_species %in% "E.falciformis",])
anova(modTolF1_fal) # no signif

modTolF1_fer <- lm(tolerance ~  hybridLevel,
                   data = ALL_summary_F1[ALL_summary_F1$Eimeria_species %in% "E.ferrisi",])
anova(modTolF1_fer) # no signif

## 3. continuous ML : test HV
fitTol <- parasiteLoad::analyse(data = ALL_summary_F1,
                                response = "tolerance", 
                                group = "Eimeria_species", model = "normal")

# [1] "Testing H0 no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.18   1 0.5476505
# [1] "Testing H1 no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.16   1 0.5765093
# [1] "Testing H2 groupA no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.3   1 0.4393456
# [1] "Testing H2 groupB no alpha vs alpha"
# dLL dDF    pvalue
# 1   0   1 0.9983804
# [1] "Testing H3 groupA no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.5   1 0.3174401
# [1] "Testing H3 groupB no alpha vs alpha"
# dLL dDF    pvalue
# 1 0.05   1 0.7523796
# [1] "Testing H1 vs H0"
# dLL dDF     pvalue
# 1 2.68   1 0.02067993
# [1] "Testing H2 vs H0"
# dLL dDF    pvalue
# 1 2.4   3 0.1867651
# [1] "Testing H3 vs H1"
# dLL dDF    pvalue
# 1 3.55   4 0.1303202
# [1] "Testing H3 vs H2"
# dLL dDF     pvalue
# 1 3.83   2 0.02173701

coef(fitTol$H1)
# L1          L2       alpha        mysd 
# 6.42333772  5.65680781 -0.04236129  0.86514631 

parasiteLoad::bananaPlot(mod = fitTol$H1,
                                        data = ALL_summary_F1,
                                        response = "tolerance",
                                        islog10 = F, group = "Eimeria_species",
                                        cols = c("orange", "darkgreen"))

coef(fitTol$H3$groupA)
# L1         L2      alpha       mysd 
# 6.1448290  5.1919438 -0.1252122  0.9120794

coef(fitTol$H3$groupB)
# L1         L2      alpha       mysd 
# 6.79413631 5.96694171 0.02684156 0.73549792


parasiteLoad::bananaPlot(mod = fitTol$H3,
                                        data = ALL_summary_F1,
                                        response = "tolerance",
                                        islog10 = F, group = "Eimeria_species",
                                        cols = c("orange", "darkgreen")) 

################ The end ################ 
