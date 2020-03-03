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
homemadeGtest <- function(full, base, round = "yes"){
  dLL = logLik(full) - logLik(base)
  dDF = base$df.residual - full$df.residual
  pvalue <- 1 - stats::pchisq(2*dLL, df=dDF)
  if(round == "yes"){
    pvalue <- round(pvalue, 6)
  }
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

# for posthoc tests
art2al_SUMdf$intFacSTRAINS <- interaction(art2al_SUMdf$infection_isolate, 
                                              art2al_SUMdf$Mouse_genotype, drop=T)
SUBsummaryDF77mice$intFacSTRAINS <- interaction(SUBsummaryDF77mice$infection_isolate, 
                                                SUBsummaryDF77mice$Mouse_genotype, drop=T)

testSignif <- function(dataframe, which){
  if(which == "RES"){
    modFULL <- glm.nb(max.OPG ~ infection_isolate*Mouse_genotype, data = dataframe)
    modPara <- glm.nb(max.OPG ~ Mouse_genotype, data = dataframe)
    modMous <- glm.nb(max.OPG ~ infection_isolate, data = dataframe)
    modinter <- glm.nb(max.OPG ~ infection_isolate+Mouse_genotype, data = dataframe)
  } else if (which == "IMP"){
    modFULL <- survreg(Surv(impact)~infection_isolate*Mouse_genotype, data = dataframe, dist="weibull")
    modPara <- survreg(Surv(impact)~Mouse_genotype, data = dataframe, dist="weibull")
    modMous <- survreg(Surv(impact)~infection_isolate, data = dataframe, dist="weibull")
    modinter <- survreg(Surv(impact)~infection_isolate+Mouse_genotype, data = dataframe, dist="weibull")
  }
  return(list(modfull = modFULL, LRT = myLRTsignificanceFactors(modFULL, modPara, modMous, modinter)))
}

testPostHoc <- function(dataframe, which){
  if(which == "RES"){
    mod_multicomp <- glm.nb(max.OPG ~ intFacSTRAINS, data = dataframe)
  } else if(which == "IMP"){
    mod_multicomp <- survreg(Surv(impact)~intFacSTRAINS, data = dataframe, dist="weibull")
  } 
  return(summary(glht(mod_multicomp, linfct=mcp(intFacSTRAINS = "Tukey"))))
}

###############################
## Test factors significance ##
###############################

# Resistance
coef(testSignif(art2al_SUMdf, "RES"))
# [1] "significance of parasite:"
# [1] "G=35.5 ,df=8 ,p=2.2e-05"
# [1] "significance of mouse:"
# [1] "G=36.3 ,df=9 ,p=3.5e-05"
# [1] "significance of interaction:"
# [1] "G=21.8 ,df=6 ,p=0.00131"
coef(testSignif(SUBsummaryDF77mice, "RES"))
# [1] "significance of parasite:"
# [1] "G=28.5 ,df=8 ,p=0.000393"
# [1] "significance of mouse:"
# [1] "G=31.8 ,df=9 ,p=0.000216"
# [1] "significance of interaction:"
# [1] "G=19.9 ,df=6 ,p=0.002856"

ggpredict(testSignif(art2al_SUMdf, "RES")$modfull, terms = c("Mouse_genotype", "infection_isolate"))
# Predicted values of max.OPG
# x = Mouse_genotype

  # # infection_isolate = Brandenburg139 (E. ferrisi)
  # x predicted std.error conf.low conf.high
  # MMd_F0 (Sc-Sc)  508998.8     0.269 300677.9  861652.3
  # MMd_F0 (St-St)  649984.3     0.269 383961.5 1100317.6
  # MMm_F0 (Bu-Bu)  503079.7     0.269 297181.3  851632.1
  # MMm_F0 (Pw-Pw)  896458.8     0.269 529560.2 1517558.3
  # 
  # # infection_isolate = Brandenburg64 (E. ferrisi)
  # x predicted std.error  conf.low conf.high
  # MMd_F0 (Sc-Sc)  456384.5     0.176  323345.7  644161.5
  # MMd_F0 (St-St)  878777.5     0.170  629926.8 1225935.9
  # MMm_F0 (Bu-Bu) 1102932.1     0.176  781421.1 1556726.8
  # MMm_F0 (Pw-Pw) 1689484.9     0.182 1181520.0 2415836.6
  # 
  # # infection_isolate = Brandenburg88 (E. falciformis)
  # x predicted std.error  conf.low conf.high
  # MMd_F0 (Sc-Sc) 1136492.8     0.269  671354.2 1923896.4
  # MMd_F0 (St-St) 2118816.9     0.249 1301479.0 3449448.5
  # MMm_F0 (Bu-Bu) 1392974.0     0.465  559719.4 3466695.5
  # MMm_F0 (Pw-Pw)  414254.5     0.329  217406.1  789337.5

# Impact
coef(testSignif(art2al_SUMdf, "IMP")) #  "G=10.3 ,df=6 ,p=0.114453"
# [1] "significance of parasite:"
# [1] "G=30.7 ,df=8 ,p=0.000159"
# [1] "significance of mouse:"
# [1] "G=23 ,df=9 ,p=0.006115"
# [1] "significance of interaction:"
# [1] "G=10.3 ,df=6 ,p=0.114453"
coef(testSignif(SUBsummaryDF77mice, "IMP"))
# [1] "significance of parasite:"
# [1] "G=38 ,df=8 ,p=8e-06"
# [1] "significance of mouse:"
# [1] "G=31.6 ,df=9 ,p=0.000235"
# [1] "significance of interaction:"
# [1] "G=23.8 ,df=6 ,p=0.000579"

## Translation of 1% because Weibull doesn't support nul data
coef(testSignif(art2al_SUMdf, "IMP")$modfull)
coefImp <- exp(coef(testSignif(art2al_SUMdf, "IMP")$modfull))
# marginal effect for each combination:
predImp <- ggpredict(testSignif(art2al_SUMdf, "IMP")$modfull, terms = c("Mouse_genotype", "infection_isolate"))
## NB substract 0.01 to each value, as it was added to model!
predImp <- data.frame(predImp)
predImp[c("predicted", "conf.low", "conf.high")]  <- predImp[c("predicted", "conf.low", "conf.high")] - 0.01
predImp
# x                 predicted  std.error conf.low   conf.high                          group
# 1  MMd_F0 (Sc-Sc) 0.08162452 0.2204047 0.04948437 0.13113037    Brandenburg139 (E. ferrisi)
# 2  MMd_F0 (Sc-Sc) 0.05407199 0.1574559 0.03705887 0.07723583     Brandenburg64 (E. ferrisi)
# 3  MMd_F0 (Sc-Sc) 0.12165590 0.2388303 0.07244183 0.20024856 Brandenburg88 (E. falciformis)
# 4  MMd_F0 (St-St) 0.10504376 0.2443933 0.06125824 0.17573384    Brandenburg139 (E. ferrisi)
# 5  MMd_F0 (St-St) 0.04110312 0.1528678 0.02787268 0.05895547     Brandenburg64 (E. ferrisi)
# 6  MMd_F0 (St-St) 0.06991219 0.2212570 0.04179391 0.11329553 Brandenburg88 (E. falciformis)
# 7  MMm_F0 (Bu-Bu) 0.06322701 0.2384227 0.03589079 0.10684687    Brandenburg139 (E. ferrisi)
# 8  MMm_F0 (Bu-Bu) 0.07587708 0.1564709 0.05319591 0.10669857     Brandenburg64 (E. ferrisi)
# 9  MMm_F0 (Bu-Bu) 0.18057140 0.2204249 0.11371766 0.28355112 Brandenburg88 (E. falciformis)
# 10 MMm_F0 (Pw-Pw) 0.08692038 0.2385478 0.05072435 0.14469181    Brandenburg139 (E. ferrisi)
# 11 MMm_F0 (Pw-Pw) 0.10082160 0.1624209 0.07060676 0.14236224     Brandenburg64 (E. ferrisi)
# 12 MMm_F0 (Pw-Pw) 0.19613163 0.2204155 0.12382171 0.30751385 Brandenburg88 (E. falciformis)

####################
## Post-hoc tests ##
####################

# to avoid running these long test all the time
doYouRun = "keepit"

if (doYouRun == "foncebebe"){
  # Resistance
  testPostHoc(art2al_SUMdf, "RES", "STRAINS")
  testPostHoc(SUBsummaryDF77mice, "RES", "STRAINS")
  # Impact
  testPostHoc(art2al_SUMdf, "IMP", "STRAINS")
  testPostHoc(SUBsummaryDF77mice, "IMP", "STRAINS")
}

#################
## save output ##
#################
doYouSave = "notthistime"
if (doYouSave == "foncebebe"){
  # Resistance
  write.csv(getMatrixPostHoc(testPostHoc(art2al_SUMdf, "RES")),
            "../figures/Tab_posthocResSTRAINS.csv")
  write.csv(getMatrixPostHoc(testPostHoc(SUBsummaryDF77mice, "RES")),
            "../figures/Tab_posthocResSTRAINS_77mice.csv")
  # Impact
  write.csv(getMatrixPostHoc(testPostHoc(art2al_SUMdf, "IMP")),
            "../figures/Tab_posthocImpSTRAINS.csv")
  write.csv(getMatrixPostHoc(testPostHoc(SUBsummaryDF77mice, "IMP")),
            "../figures/Tab_posthocImpSTRAINS_77mice.csv")
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

## NB. translate back 0.01
transValuesImp <- seq(0.01,0.31, 0.05)
as.character(transValuesImp - 0.01)
realValuesImpLabels <- c("0%", "5%", "10%", "15%", "20%", "25%", "30%")

get_plotI_STRAINS <- function(dataframe){
  plot_model(testSignif(dataframe, "IMP")$modfull,
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

########################################
### Second part: assessing tolerance ###
########################################
toleranceAnalysis <- function(dataset){
  modfull <- lm(relWL ~ 0 + max.OPG : (infection_isolate * Mouse_genotype), data = dataset)
  modNOmouse <- lm(relWL ~ 0 + max.OPG : (infection_isolate), data = dataset)
  modNOisolate <- lm(relWL ~ 0 + max.OPG : (Mouse_genotype), data = dataset)
  modNOint <- lm(relWL ~ 0 + max.OPG : (infection_isolate + Mouse_genotype), data = dataset)
  
  testIsolate <- homemadeGtest(modfull, modNOisolate)
  testMouse <- homemadeGtest(modfull, modNOmouse)
  testInteraction <- homemadeGtest(modfull, modNOint)

  # Calculate slopes for each group:
  predTolSlopes <- ggpredict(modfull, terms = c("Mouse_genotype", "infection_isolate"), condition = c(max.OPG = 1000000))  ## And plot
  predTolSlopes <- data.frame(predTolSlopes)
  names(predTolSlopes)[names(predTolSlopes) %in% c("x", "group")] <- c("Mouse_genotype", "infection_isolate")
  predTolSlopes$group <- paste0(predTolSlopes$Mouse_genotype, predTolSlopes$infection_isolate)
  
  # make line up to 5e6 OPG for plot
  pts <- predTolSlopes
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
  
  T1 <- ggplot(pts, aes(x = max.OPG, y = relWL, col = Mouse_genotype)) +
    geom_line(aes(group = group)) +
    facet_grid(.~infection_isolate) +
    scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                       name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +  
    scale_x_continuous("maximum oocysts per gram of feces (x10e6)", 
                       breaks = seq(0, 5000000, 1000000),
                       labels = seq(0, 5000000, 1000000)/1000000) +
    scale_y_continuous(name = "maximum weight loss compared to day of infection",
                       breaks = seq(0,0.3, 0.05), 
                       labels = c("0%", "5%", "10%", "15%", "20%", "25%", "30%")) +
    geom_point(data = dataset, size = 4, pch = 1)+
    coord_cartesian(ylim=c(0, 0.30)) +
    theme(legend.position = "top")
  predTolSlopes$predicted = round(predTolSlopes$predicted,2)
  return(list(testIsolate = testIsolate, testMouse = testMouse, testInteraction = testInteraction,
              predTolSlopes = predTolSlopes, T1 = T1))
}

t <- toleranceAnalysis(dataset = art2al_SUMdf)
# "G=30.2 ,df=8 ,p=0.000197"
# "G=30.6 ,df=9 ,p=0.000341"
# "G=24 ,df=6 ,p=0.000513"
pdf(file = "../figures/Fig4.pdf",
    width = 12, height = 6)
t$T1
dev.off()

##### Spearman
func <- function(dataset, what){
  # dataset = dataset[c("relWL", "Mouse_genotype", "infection_isolate", "max.OPG")]
  # add a null point for zero intercept
  dataset <- rbind(dataset[1,], dataset)
  dataset[1,c("relWL", "max.OPG")] <- c(0,0)
  c <- cor.test(dataset$relWL, dataset$max.OPG, method = "spearman")
  if (what == "pvalue"){
    round(c$p.value, 3)
  } else if (what == "N"){
    nrow(na.omit(data.frame(dataset$relWL, dataset$max.OPG))) -1
  } else if (what == "statistic"){
    round(c$statistic, 2)
  } else if (what == "rho"){
    round(c$estimate, 2)
    }
}

func2 <- function(dataset){
  dataset$group = paste0(dataset$Mouse_genotype, dataset$infection_isolate)  
  P <- plyr::ddply(dataset, .(group), func, what = "pvalue")
  names(P)[names(P) %in% "V1"] <- "pvalue"
  N <- plyr::ddply(dataset, .(group), func, what = "N")
  names(N)[names(N) %in% "V1"] <- "N"
  S <- plyr::ddply(dataset, .(group), func, what = "statistic")
  names(S)[names(S) %in% "V1"] <- "statistic"
  R <- plyr::ddply(dataset, .(group), func, what = "rho")
  names(R)[names(R) %in% "V1"] <- "rho"
  merge(merge(merge(P, N),S), R)
}

final <- merge(func2(art2al_SUMdf), t$predTolSlopes)
final <- final[order(final$infection_isolate),]
write.csv(final, file = "../figures/TabToleranceSpearman.csv", row.names = F)

### with conservative dataset
tol77mice <- toleranceAnalysis(dataset = SUBsummaryDF77mice)
tol77mice$T1

pdf(file = "../figures/tol77mice.pdf",
    width = 12, height = 6)
tol77mice$FigTOL
dev.off()

final77 <- merge(func2(SUBsummaryDF77mice), tol77mice$predTolSlopes)
final77 <- final77[order(final77$infection_isolate),]
final77

write.csv(final77, file = "../figures/TabToleranceSpearman_77mice.csv", row.names = F)
