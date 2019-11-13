### Code for data analysis of Article 2
### August 2019
### Alice Balard

## INFO
# Mouse AA_0088, HI = 0.2
# Mouse AA_0064, HI = 0.08
# Mouse AA_0139, HI = 0.85

#### Load data and functions within ####
source("dataPreparation.R")

###### Map of our samples FIGURE 1 (with design) ######
hmhzline <- read.csv("../data/HMHZ.csv")
area <- get_stamenmap(bbox = c(9, 49, 17, 53.7), zoom = 7,
                      maptype = "toner-lite")

map <- ggmap(area) +
  geom_path(hmhzline, mapping =  aes(x = lon, y = lat), col = "purple", size = 4) +
  geom_label_repel(data = forMap,
                   aes(longitude, latitude, label = Name),
                   box.padding = 2, size = 5) +
  geom_point(data = forMap, aes(longitude, latitude, col = color), size = 6) +
  scale_color_manual(values = as.character(levels(forMap$color))) +
  theme_bw() +
  theme(legend.position = 'none', axis.ticks=element_blank())
map 

pdf(file = "../figures/Fig1.pdf", width = 5, height =5)
map
dev.off()

######################################
########## Read information ##########
######################################

table(summaryDF108mice$infection_isolate, summaryDF108mice$Mouse_genotype, summaryDF108mice$anth)

## Age of mice
range(as.numeric(rawDF108mice$ageAtInfection))

###### what is the overall peak day for each parasite isolate? ######
aggregate(summaryDF108mice$dpi_max.oocysts.per.tube,
          list(summaryDF108mice$infection_isolate), 
          function(x) {paste(length(x), median(x), round(sd(x),2))})
aggregate(summaryDF108mice$dpi_minWeight,
          list(summaryDF108mice$infection_isolate), 
          function(x) {paste(length(x), median(x), round(sd(x),2))})

###### what is the overall prepatent period for each parasite isolate? ######
d <- as.data.frame(
  rawDF108mice[!is.na(rawDF108mice$oocysts.per.tube) & rawDF108mice$oocysts.per.tube > 0,] %>% 
    dplyr::group_by(EH_ID) %>%
    dplyr::slice(which.min(dpi)) %>%
    dplyr::select(EH_ID, weight, HI, startingWeight, ageAtInfection, Sex,
                  Mouse_genotype, Eimeria_species, Mouse_subspecies,
                  infection_isolate, Exp_ID, dpi))
aggregate(d$dpi,
          list(d$infection_isolate), 
          function(x) {paste(length(x), median(x), round(sd(x),2))})

###### Course of infection FIGURE 2 ######
rawDF108mice$parasiteDensity <- rawDF108mice$oocysts.per.tube / rawDF108mice$weight

forplot <- rawDF108mice %>%
  group_by(infection_isolate, dpi) %>%
  summarise(mean = mean(parasiteDensity*10e-6, na.rm = TRUE),
            sd = sd(parasiteDensity*10e-6, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)

F2.1 <- ggplot(forplot, aes(dpi, mean, group = infection_isolate, col = infection_isolate)) + 
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = .2)+
  ylab("parasite density (10e6 oocysts per mouse gram)") +
  scale_x_continuous(breaks = 0:11, name = "days post infection") +
  theme(legend.position = c(0.25, 0.8)) +
  labs(color = "Eimeria isolate") 

forplot2 <- rawDF108mice %>%
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

library(cowplot)
Fig2 <- plot_grid(F2.1, F2.2,
                  labels=c("A", "B"), label_size = 20)
pdf(file = "../figures/Fig2.pdf", width = 10, height = 5)
Fig2
dev.off()

## Correlation sum of oocysts / peak oocysts
ggplot(summaryDF108mice, aes(sumoocysts.per.tube, max.oocysts.per.tube)) + 
  geom_smooth(method = "lm")+ geom_point()

cor(summaryDF108mice$sumoocysts.per.tube , summaryDF108mice$max.oocysts.per.tube,
    method = "pearson")

###############################################################
########## Define our indexes and their distribution ##########
###############################################################

## RESISTANCE
# for models, inverse of parasite density:
summaryDF108mice$peak.oocysts.per.g.mouse
# for plots, invert + translation positive (300000) / 10000:
getResistanceIndex <- function(x){
  y = (- x + 300000)/300000
  return(y)}

summaryDF108mice$ResistanceIndex <- getResistanceIndex(summaryDF108mice$peak.oocysts.per.g.mouse)
plotChoiceResIndex <- ggplot(summaryDF108mice, aes(peak.oocysts.per.g.mouse, ResistanceIndex)) + 
  geom_point(size = 4, pch =21) + 
  ylab("Resistance Index") +
  xlab("Parasite density (oocysts per mouse gram) at peak day")
plotChoiceResIndex

xRes <- as.numeric(na.omit(summaryDF108mice$peak.oocysts.per.g.mouse))
hist(xRes, breaks = 100)
descdist(xRes)
pdf("../figures/supfig1.1.pdf")
findGoodDist(x = xRes, distribs = c("normal", "negative binomial"), 
             distribs2 = c("norm", "nbinom"))
dev.off()
### nbinom for resistance

## IMPACT ON HEALTH
xImp <- as.numeric(na.omit(summaryDF108mice$relWL))
hist(xImp, breaks = 100)
descdist(xImp)
pdf("../figures/supfig1.2.pdf")
findGoodDist(x = xImp+ 0.01, distribs = c("normal", "weibull"), 
             distribs2 = c("norm", "weibull"))
dev.off()
### weibull for impact on health
summaryDF108mice$impact <- summaryDF108mice$relWL +0.01

## TOLERANCE
summaryDF108mice$ToleranceIndex <- log10(
  summaryDF108mice$relWL / summaryDF108mice$peak.oocysts.per.g.mouse + 1e-8) / (-8)

plotChoiceTolIndex <- ggplot(summaryDF108mice, 
                             aes(x=peak.oocysts.per.g.mouse, y =relWL, fill = ToleranceIndex))+
  geom_point(size = 4, pch =21) +
  scale_fill_gradient(low = "white", high = "black") +
  xlab("Parasite density (oocysts per mouse gram) at peak day") +
  ylab("Maximum relative weight loss compared to day 0 (%)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  labs(fill = "Tolerance Index") +
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )
  
plotChoiceTolIndex

xTol <- as.numeric(na.omit(summaryDF108mice$ToleranceIndex))

# for tolerance, we need to find a way to deal with the massive extreme values
hist(xTol, breaks = 1000)
descdist(xTol)

pdf("../figures/supfig1.3.pdf")
findGoodDist(x = xTol, distribs = c("normal"), 
             distribs2 = c("norm"))
dev.off()

summaryDF108mice[is.na(summaryDF108mice$ToleranceIndex), c("EH_ID", "Mouse_genotype", "relWL", "peak.oocysts.per.g.mouse")]
# 9 mice died before peak

#
pdf(file = "../figures/choiceIndexes.pdf", width = 10, height = 5)
plot_grid(plotChoiceResIndex,
          plotChoiceTolIndex,
          labels=c("A", "B"), label_size = 15)
dev.off()

################################
##### Statistical analyses #####
################################

################
## Resistance ##
################
# Mmd vs Mmm
modResSubsp <- glm.nb(peak.oocysts.per.g.mouse ~ Eimeria_species*Mouse_subspecies, data = summaryDF108mice)

modResSubsp
anova(modResSubsp)
summary(modResSubsp)
a <- aov(modResSubsp)
summary.lm(a)
# SIGNIF infection isolate (p-value = 0.02236) + interactions with mice (p-value = 6.52e-07)

summary(modResSubsp)

# by strains
modResStrain <- glm.nb(peak.oocysts.per.g.mouse ~ infection_isolate*Mouse_genotype, 
                 data = summaryDF108mice)
anova(modResStrain, test = "LRT")
# SIGNIF infection isolate (p-value = 0.01902) + interactions with mice (p-value = 8.432e-05)

# postHoc test
summaryDF108mice$intFac <- interaction(summaryDF108mice$infection_isolate, 
                                       summaryDF108mice$Mouse_genotype, drop=T)

modResStrainformulticomp <- glm.nb(peak.oocysts.per.g.mouse ~ intFac, data = summaryDF108mice)
postHocRes <- summary(glht(modResStrainformulticomp, linfct=mcp(intFac = "Tukey")))

# Create a matrix to present the post-hoc tests 
getFullMatrix <- function(postHoc){
  
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
  mat1 <- getMatrix(postHocRes, upperTriangle)
  ### make matrix 2: z value and P
  mat2 <- getMatrix(postHocRes, lowerTriangle)
  mat3 <- t(mat2)
  ## Combine
  mat <- mat1
  mat[] <- lapply(mat, as.character)
  mat[lower.tri(mat)] <- mat3[lower.tri(mat3)]
  mat <- as.matrix(mat)
  diag(mat) <- ""
  return(mat)
}

myMatpostHocRes <- getFullMatrix(postHocRes)
write.csv(myMatpostHocRes, "../figures/supTablePostHocRes.csv")

## Results:
# Brandenburg88 (E. falciformis).MMm_F0 (Pw-Pw) - Brandenburg64 (E. ferrisi).MMd_F0 (Sc-Sc) == 0       0.0294 *  
# Brandenburg88 (E. falciformis).MMm_F0 (Pw-Pw) - Brandenburg88 (E. falciformis).MMd_F0 (Sc-Sc) == 0   0.0368 *  
# Brandenburg88 (E. falciformis).MMm_F0 (Pw-Pw) - Brandenburg64 (E. ferrisi).MMd_F0 (St-St) == 0        <0.01 ** 
# Brandenburg88 (E. falciformis).MMm_F0 (Pw-Pw) - Brandenburg88 (E. falciformis).MMd_F0 (St-St) == 0    <0.01 ***
# Brandenburg88 (E. falciformis).MMm_F0 (Pw-Pw) - Brandenburg139 (E. ferrisi).MMm_F0 (Bu-Bu) == 0      0.0318 *  
# Brandenburg88 (E. falciformis).MMm_F0 (Bu-Bu) - Brandenburg64 (E. ferrisi).MMm_F0 (Bu-Bu) == 0       0.0189 *  
# Brandenburg88 (E. falciformis).MMm_F0 (Pw-Pw) - Brandenburg64 (E. ferrisi).MMm_F0 (Bu-Bu) == 0        <0.01 ***
# Brandenburg64 (E. ferrisi).MMm_F0 (Pw-Pw) - Brandenburg88 (E. falciformis).MMm_F0 (Bu-Bu) == 0       0.0226 *  
# Brandenburg88 (E. falciformis).MMm_F0 (Pw-Pw) - Brandenburg139 (E. ferrisi).MMm_F0 (Pw-Pw) == 0       <0.01 ** 
# Brandenburg88 (E. falciformis).MMm_F0 (Pw-Pw) - Brandenburg64 (E. ferrisi).MMm_F0 (Pw-Pw) == 0        <0.01 ***

## And plot:
## To add Ns on top of bars
getNs <- function(proxy, df, groupMus = "Mouse_genotype", groupPar = "infection_isolate"){
  noNA = df[!is.na(df[[proxy]]),]
  noNA$groupMus = noNA[[groupMus]]
  noNA$groupPar = noNA[[groupPar]]
  tab = table(noNA$groupPar, noNA$groupMus)
  Ns = as.character(as.vector(t(tab)[as.vector(t(tab))!=0]))
  return(Ns)
}

## Fig 3
# plot marginal effects of interaction terms
posx.1 <- c(0.9,1.1, 1.9,2.1)

## FOR PLOT, use Resistance Index 
y = seq(0,1,0.05)
x = -y * 300000 + 300000

plotR_F0_subsp <- plot_model(modResSubsp, type = "int", dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
  scale_color_manual(values = c("blue", "red"),
                     name = "Mouse subspecies",labels = c("Mmd", "Mmm")) +
  ggtitle("Resistance") +
  scale_y_continuous("(predicted) Resistance Index",
                     trans = "reverse",
                     breaks = x,
                     labels = y) +
  xlab("Eimeria species") +
  theme(axis.title.x = element_text(hjust=1), axis.text=element_text(size=13)) +
  geom_text(aes(x=posx.1,y=90000,label=getNs("peak.oocysts.per.g.mouse", summaryDF108mice, 
                                            "Mouse_subspecies", "Eimeria_species")), vjust=0) 

plotR_F0_subsp

# plot marginal effects of interaction terms by isolates & strains
posx.2 <- c(0.8+c(0,1/8,2/8,3/8),1.8+c(0,1/8,2/8,3/8),2.8+c(0,1/8,2/8,3/8))
## and plot
plotR_F0_strains <- plot_model(modResStrain, type = "int", dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
  scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                     name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +
  ggtitle("Resistance") +
  scale_y_continuous("(predicted) Resistance Index",
                     trans = "reverse",
                     breaks = x,
                     labels = y) +
  xlab("Eimeria isolate") +
  theme(axis.title.x = element_text(hjust=1), axis.text=element_text(size=13)) +
  geom_text(aes(x=posx.2,y=120000,label=getNs("peak.oocysts.per.g.mouse", summaryDF108mice)),vjust=0) 
plotR_F0_strains

############
## Impact ##
############
modImpSubsp <- survreg(Surv(impact)~Eimeria_species*Mouse_subspecies, data = summaryDF108mice, dist="weibull")
anova(modImpSubsp) # Eimeria species AND mouse subspecies significant

## Translation of 1% because Weibull doesn't support nul data
modImpStrain <- survreg(Surv(impact)~infection_isolate*Mouse_genotype, data = summaryDF108mice, dist="weibull")
anova(modImpStrain)
length(summaryDF108mice$relWL)
# Eimeria isolate significant

## post-hoc Tukey test
modImpStrainformulticomp <- survreg(Surv(impact)~intFac, data = summaryDF108mice, dist="weibull")
postHocImp <- summary(glht(modImpStrainformulticomp, linfct=mcp(intFac = "Tukey")))

myMatpostHocImp <- getFullMatrix(postHocImp)
write.csv(myMatpostHocImp, "../figures/supTablePostHocImp.csv")

# Results (significant only)  
# Brandenburg88 (E. falciformis).MMm_F0 (Bu-Bu) - Brandenburg64 (E. ferrisi).MMd_F0 (Sc-Sc) == 0        <0.01 ** 
# Brandenburg88 (E. falciformis).MMm_F0 (Pw-Pw) - Brandenburg64 (E. ferrisi).MMd_F0 (Sc-Sc) == 0        <0.01 ***
# Brandenburg64 (E. ferrisi).MMd_F0 (St-St) - Brandenburg88 (E. falciformis).MMd_F0 (Sc-Sc) == 0       0.0360 *  
# Brandenburg88 (E. falciformis).MMm_F0 (Bu-Bu) - Brandenburg64 (E. ferrisi).MMd_F0 (St-St) == 0        <0.01 ***
# Brandenburg64 (E. ferrisi).MMm_F0 (Pw-Pw) - Brandenburg64 (E. ferrisi).MMd_F0 (St-St) == 0           0.0222 *  
# Brandenburg88 (E. falciformis).MMm_F0 (Pw-Pw) - Brandenburg64 (E. ferrisi).MMd_F0 (St-St) == 0        <0.01 ***

plotI_F0_subsp <- plot_model(modImpSubsp, type = "int",dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
  scale_color_manual(values = c("blue","red"),
                     name = "Mouse subspecies",labels = c("Mmd", "Mmm")) +
  xlab("Eimeria species") +
  ggtitle("Impact on host health") +
  ylab("(predicted) maximum weight loss") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme(axis.title.x = element_text(hjust=1), axis.text=element_text(size=13)) +
  geom_text(aes(x=posx.1,y=0,label=getNs("relWL", summaryDF108mice,
                                         "Mouse_subspecies", "Eimeria_species")),vjust=0)
plotI_F0_subsp

plotI_F0_strains <- plot_model(modImpStrain, type = "int",dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
  scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                     name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +
  xlab("Eimeria isolate") +
  ggtitle("Impact on host health") +
  ylab("(predicted) maximum weight loss") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme(axis.title.x = element_text(hjust=1), axis.text=element_text(size=13)) +
  geom_text(aes(x=posx.2,y=0,label=getNs("relWL", summaryDF108mice)),vjust=0)
plotI_F0_strains

###############
## Tolerance ##
###############
modTolSubspecies <- lm(ToleranceIndex ~ Eimeria_species*Mouse_subspecies, data = summaryDF108mice)
anova(modTolSubspecies)

length(na.omit(summaryDF108mice$ToleranceIndex))
modTolStrain <- lm(ToleranceIndex ~ infection_isolate*Mouse_genotype, data = summaryDF108mice)
anova(modTolStrain)
# mouse genotype & interactions significant
TukeyHSD(modTolStrain)
# Brandenburg88 (E. falciformis):MMm_F0 (Pw-Pw)-Brandenburg64 (E. ferrisi):MMm_F0 (Pw-Pw)     0.0276440
# Brandenburg88 (E. falciformis):MMm_F0 (Pw-Pw)-Brandenburg64 (E. ferrisi):MMm_F0 (Bu-Bu)     0.0187152
# Brandenburg88 (E. falciformis):MMm_F0 (Pw-Pw)-Brandenburg88 (E. falciformis):MMd_F0 (St-St) 0.0197022
# Brandenburg88 (E. falciformis):MMm_F0 (Pw-Pw)-Brandenburg139 (E. ferrisi):MMd_F0 (St-St)    0.0009261
# Brandenburg88 (E. falciformis):MMm_F0 (Pw-Pw)-Brandenburg88 (E. falciformis):MMd_F0 (Sc-Sc) 0.0385139
# Brandenburg88 (E. falciformis):MMm_F0 (Pw-Pw)-Brandenburg64 (E. ferrisi):MMd_F0 (Sc-Sc)     0.0068948
# Brandenburg88 (E. falciformis):MMm_F0 (Pw-Pw)-Brandenburg64 (E. ferrisi):MMd_F0 (St-St)     0.0008689

## post-hoc Tukey test -> for lm should be the same results with TukeyHSD and glht: indeed :) 
modTolStrainformulticomp <- lm(ToleranceIndex ~ intFac, data = summaryDF108mice)
  postHocTol <- summary(glht(modTolStrainformulticomp, linfct=mcp(intFac = "Tukey")))
# Results (significant only)
# Brandenburg88 (E. falciformis).MMm_F0 (Pw-Pw) - Brandenburg64 (E. ferrisi).MMd_F0 (Sc-Sc) == 0        <0.01 ** 
# Brandenburg88 (E. falciformis).MMm_F0 (Pw-Pw) - Brandenburg88 (E. falciformis).MMd_F0 (Sc-Sc) == 0   0.0356 *  
# Brandenburg88 (E. falciformis).MMm_F0 (Pw-Pw) - Brandenburg139 (E. ferrisi).MMd_F0 (St-St) == 0       <0.01 ***
# Brandenburg88 (E. falciformis).MMm_F0 (Pw-Pw) - Brandenburg64 (E. ferrisi).MMd_F0 (St-St) == 0        <0.01 ***
# Brandenburg88 (E. falciformis).MMm_F0 (Pw-Pw) - Brandenburg88 (E. falciformis).MMd_F0 (St-St) == 0   0.0185 *  
# Brandenburg88 (E. falciformis).MMm_F0 (Pw-Pw) - Brandenburg64 (E. ferrisi).MMm_F0 (Bu-Bu) == 0       0.0169 *  
# Brandenburg88 (E. falciformis).MMm_F0 (Pw-Pw) - Brandenburg64 (E. ferrisi).MMm_F0 (Pw-Pw) == 0       0.0256 *  

myMatpostHocTol<- getFullMatrix(postHocTol)
write.csv(myMatpostHocTol, "../figures/supTablePostHocTol.csv")

plotT_F0_subsp <- plot_model(modTolSubspecies, type = "int", dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
  scale_color_manual(values = c("blue", "red"),
                     name = "Mouse subspecies",labels = c("Mmd", "Mmm")) +
  xlab("Eimeria species") +
  ylab("(predicted) Tolerance index")+
  ggtitle("Tolerance") +
  theme(axis.title.x = element_text(hjust=1), axis.text = element_text(size=13))+
  geom_text(aes(x=posx.1,y=0.4,label=getNs("ToleranceIndex", summaryDF108mice,
                                           "Mouse_subspecies", "Eimeria_species")),vjust=0)
plotT_F0_subsp

plotT_F0_strain <- plot_model(modTolStrain, type = "int", dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
  scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                     name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +
  xlab("Eimeria isolate") +
  ylab("(predicted) Tolerance index")+
  ggtitle("Tolerance") +
  theme(axis.title.x = element_text(hjust=1), axis.text = element_text(size=13))+
  geom_text(aes(x=posx.2,y=0.4,label=getNs("ToleranceIndex", summaryDF108mice)),vjust=0)
plotT_F0_strain

# Fig 3.
Fig3 <- plot_grid(plotR_F0_subsp + theme(legend.position = "none"),
                  plotI_F0_subsp + theme(legend.position = "none"),
                  plotT_F0_subsp+ theme(legend.position = "none"), 
                  plotT_F0_subsp ,
                  labels=c("A", "B", "C", "D"), label_size = 20)

Fig3
pdf(file = "../figures/Fig3.pdf",
    width = 9, height = 9)
Fig3
dev.off()

## Fig 4
Fig4 <- plot_grid(plotR_F0_strains + theme(legend.position = "none"),
                  plotI_F0_strains + theme(legend.position = "none"),
                  plotT_F0_strain + theme(legend.position = "none"), 
                  plotT_F0_strain,
                  labels=c("A", "B", "C", "D"), label_size = 20)

Fig4

pdf(file = "../figures/Fig4.pdf",
    width = 9, height = 9)
Fig4
dev.off()

# Pretty table outputs
# library(stargazer)
# stargazer(modResSubsp, modImpSubsp, modTolSubspecies)#,  type = "html")
# stargazer(modResStrain, modImpStrain, modTolStrain,  type = "html")

### Second part: correlation resistance / tolerance

# Calculate mean per group
# Packages we need
library(ggplot2)
library(dplyr)

# Create a group-means data set
gd <- summaryDF108mice %>%
  group_by(Mouse_genotype, infection_isolate) %>% 
  summarise(
    ResistanceIndexMean = mean(ResistanceIndex, na.rm = T),
    ResistanceIndexSd = sd(ResistanceIndex, na.rm = T),
    Impact = mean(impact, na.rm=T),
    ToleranceIndexMean = mean(ToleranceIndex, na.rm = T),
    ToleranceIndexSd = sd(ToleranceIndex, na.rm = T)
  )

gd

# define the 8 groups
summaryDF108mice$group <- paste(summaryDF108mice$Mouse_genotype, summaryDF108mice$infection_isolate, sep = "_")

# Plot both data sets
restolplot <- ggplot(summaryDF108mice, aes(x = ResistanceIndex, y = ToleranceIndex)) +
  geom_smooth(method = "lm", col = "black", alpha = .2, aes(linetype = Eimeria_species)) +
  geom_point(alpha = .4, aes(col = Mouse_genotype, fill = Mouse_genotype, shape = infection_isolate), size = 4) +
  geom_point(data = gd, aes(x = ResistanceIndexMean, y = ToleranceIndexMean,
                            fill = Mouse_genotype, shape = infection_isolate), size = 10) +
  theme_bw()+
  scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1")) +
  scale_fill_manual(values = c("blue", "cornflowerblue", "red4", "indianred1")) +
  scale_shape_manual(values = c(24,22,21)) +
  ylab(label = "Tolerance index") +
  scale_x_continuous(name = "Resistance index")
restolplot

pdf("../figures/Fig5.pdf", width = 12, height = 9)
restolplot  
dev.off()

# https://stats.idre.ucla.edu/r/seminars/interactions-r/
# Test the difference between slopes:
modResTol <- lm(formula = ToleranceIndex ~ ResistanceIndex * Eimeria_species, data = summaryDF108mice)

summary(modResTol)
#p-value of the t-statistic for the interaction between ResistanceIndex and Eimeria_species: 1.39e-05 ***

# ResistanceIndex:Eimeria_speciesE.falciformis -1.17999    0.25744  -4.584 1.39e-05 ***
# The interaction ResistanceIndex * Eimeria_species is significant, which suggests
# that the relationship of ToleranceIndex by ResistanceIndex does vary by Eimeria spieces

# Since our goal is to obtain simple slopes of parasite:
library(lsmeans)
emtrends(modResTol, ~ Eimeria_species, var="ResistanceIndex")

# The 95% confidence interval does not contain zero for females but contains zero for males, so the simple slope is significant for females but not for males.

# plot with sd 
ggplot(summaryDF108mice, aes(x = ResistanceIndex, y = ToleranceIndex)) +
  geom_smooth(method = "lm", col = "black", alpha = .2, aes(linetype = Eimeria_species)) +
  geom_point(alpha = .4, aes(col = Mouse_genotype, fill = Mouse_genotype, shape = infection_isolate), size = 4) +
  geom_point(data = gd, aes(x = ResistanceIndexMean, y = ToleranceIndexMean,
                            fill = Mouse_genotype, shape = infection_isolate), size = 10) +
  geom_errorbar(data = gd, aes(x = ResistanceIndexMean, y = ToleranceIndexMean, col = Mouse_genotype,
                               ymin = ToleranceIndexMean - ToleranceIndexSd, 
                               ymax = ToleranceIndexMean + ToleranceIndexSd)) +
  geom_errorbarh(data = gd, aes(x = ResistanceIndexMean, y = ToleranceIndexMean, col = Mouse_genotype,
                               xmin = ResistanceIndexMean - ResistanceIndexSd, 
                               xmax = ResistanceIndexMean + ResistanceIndexSd)) +
  theme_bw()+
  scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1")) +
  scale_fill_manual(values = c("blue", "cornflowerblue", "red4", "indianred1")) +
  scale_shape_manual(values = c(24,22,21)) +
  ylab(label = "Tolerance index") +
  scale_x_continuous(name = "Resistance index")

##########
# Some toy
rawDF108mice$oocysts.per.g.mouse <- rawDF108mice$oocysts.per.tube / rawDF108mice$weight

ggplot(rawDF108mice, aes(x = oocysts.per.g.mouse, y =relativeWeight,
                         color = Mouse_genotype, shape = infection_isolate, fill = infection_isolate)) +
  geom_smooth(method = "lm")+
  # geom_point() +
  # scale_x_log10()+
  theme_bw()+
  scale_fill_manual(values = c(1,2,3)) +
  scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1")) +
  scale_shape_manual(values = c(15,16,10)) 
  
# Bonus: Disease trajectory?

# NB quite a bunch of animals died before end

infDF <- rawDF108mice %>%
  dplyr::group_by(Mouse_genotype, infection_isolate, dpi)%>%
  dplyr::summarise(meanOO = mean(oocysts.per.tube, na.rm = T),
                   meanWR = mean(weightRelativeToInfection, na.rm = T)) %>%
  as.data.frame()

infDF$group <- paste(infDF$Mouse_genotype, infDF$infection_isolate, sep = "_")

ggplot(infDF, aes(x = meanOO , y =meanWR, col = as.factor(dpi), group = group)) +
  geom_point() + 
  geom_path() +
  facet_grid(infection_isolate~Mouse_genotype) +
  scale_x_log10() +
  geom_hline(yintercept = 100) +
  geom_label(aes(label = dpi))


