### Code for data analysis of Article 2
### August 2019
### Alice Balard

## INFO
# Mouse AA_0088, HI = 0.2
# Mouse AA_0064, HI = 0.08
# Mouse AA_0139, HI = 0.85

#### Load data and functions within ####
source("dataPreparation_168mice.R")
library(cowplot)
library(ggplot2)
library(dplyr)
library(lmtest)
library(pscl)
library(scales)
mycolors <- scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=8))

## Different datasets as follow:

# FULL = DSart2 / art2SummaryDF
nrow(art2SummaryDF) # 168

# conservative = remove mice with contamination or anthelminthic
# DSart2_conservative2 ; art2SummaryDF_conservative2 # 118 mice
nrow(art2SummaryDF_conservative)

# for further plots:
art2SummaryDF$label <- art2SummaryDF$Mouse_genotype
levels(art2SummaryDF$label) <- as.character(c(1:8))
art2SummaryDF$Genotype <- paste0(art2SummaryDF$label, ". ", art2SummaryDF$Mouse_genotype)
art2SummaryDF_conservative$label <- art2SummaryDF_conservative$Mouse_genotype
levels(art2SummaryDF_conservative$label) <- as.character(c(1:8))
art2SummaryDF_conservative$Genotype <- paste0(art2SummaryDF_conservative$label, ". ", art2SummaryDF_conservative$Mouse_genotype)

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

### Table 1:Infection experiment design
tab1 <- table(art2SummaryDF$Mouse_genotype,art2SummaryDF$Mouse_subspecies,art2SummaryDF$Sex,art2SummaryDF$infection_isolate)
tab1 <- data.frame(tab1)
tab1 <- tab1[tab1$Freq != 0,]
tab1wide <- dcast(tab1, Var1 + Var2 ~ Var4 + Var3, value.var="Freq")
tab1wide
write.csv(tab1wide, "../figures/Table1_temp.csv", row.names = F)

###### how many mice died before shedding oocysts? ######
table(art2SummaryDF$sumoocysts.per.tube == 0)
# 5 mice did not shed oocysts, all infected by E88
deadNoOO <- art2SummaryDF[art2SummaryDF$sumoocysts.per.tube == 0,"EH_ID"]
dfNoOO <- DSart2[DSart2$EH_ID %in% deadNoOO,]
ggplot(dfNoOO, aes(x=dpi, y=relWL, group=EH_ID, col=EH_ID)) + geom_line() + geom_point()

###### what is the overall peak day for each parasite isolate? ######
aggregate(art2SummaryDF$dpi_max.OPG,
          list(art2SummaryDF$infection_isolate), 
          function(x) {paste(length(x), median(x), round(sd(x),2))})
# 1    Brandenburg139 (E. ferrisi) 25 6 0.73
# 2     Brandenburg64 (E. ferrisi) 87 6 0.86
# 3 Brandenburg88 (E. falciformis) 56 8 1.34
aggregate(art2SummaryDF$dpi_minWeight,
          list(art2SummaryDF$infection_isolate), 
          function(x) {paste(length(x), median(x), round(sd(x),2))})
# Brandenburg139 (E. ferrisi) 25 5 2.08
# Brandenburg64 (E. ferrisi) 87 5 1.69
# Brandenburg88 (E. falciformis)  56 9 1.58

## Make table with batches (only batch 1 was treated with anthelminthics)
tS1 <- data.frame(table(art2SummaryDF$Batch, art2SummaryDF$Mouse_genotype, art2SummaryDF$infection_isolate))
tS1 <- tS1[order(tS1$Var1) & tS1$Freq != 0,]
tS1 <- tS1[order(tS1$Var1),]
tS1 <- dcast(tS1, Var1 + Var3 ~ Var2, value.var="Freq")
tS1
write.csv(tS1, "../figures/TableS1_temp.csv", row.names = F) # NB done for FULL DS

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
# 2     Brandenburg64 (E. ferrisi) 87 5 1.89
# 3 Brandenburg88 (E. falciformis) 51 7 3.18

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
  geom_line(aes(linetype=infection_isolate)) +
  scale_linetype_manual(values = c(1,2,1)) +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = .2)+
  ylab("million oocysts per gram of feces") +
  scale_x_continuous(breaks = 0:11, name = "days post infection") +
  scale_color_manual(values = c("darkgreen", "#00c000ff", "orange"))+
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
  geom_line(aes(linetype=infection_isolate)) +
  scale_linetype_manual(values = c(1,2,1)) +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = .2)+
  ylab("relative weight compared to day 0 (%)") +
  scale_x_continuous(breaks = 0:11, name = "days post infection") +
  scale_color_manual(values = c("darkgreen", "#00c000ff", "orange"))+
  theme(legend.position = c(0.25, 0.2)) +
  labs(color = "Eimeria isolate") 

Fig2 <- cowplot::plot_grid(F2.1, F2.2,
                           labels=c("A", "B"), label_size = 20)

pdf(file = "../figures/Fig2_temp.pdf", width = 10, height = 5)
Fig2
dev.off()

## Correlation sum of oocysts / peak oocysts
corSumMax <- cor.test(art2SummaryDF$sumoocysts.per.tube, 
                      art2SummaryDF$max.oocysts.per.tube,
                      method = "spearman", exact=F)
ggplot(art2SummaryDF, aes(sumoocysts.per.tube, max.oocysts.per.tube)) + 
  geom_smooth(method = "lm")+ geom_point()+
  geom_label(aes(label=paste0("Spearman coefficient: ", as.character(round(corSumMax$estimate, 2)),
                              "\np-value= ",  as.character(signif(corSumMax$p.value, digits=2))),
                 x =6e6, y =2e6)) +
  xlab("Sum of oocyst") + ylab("Oocysts at peak day")

###############################################################
########## Define our indexes and their distribution ##########
###############################################################
## RESISTANCE: inverse of OPG
xRes <- round(as.numeric(na.omit(art2SummaryDF$max.OPG)))
hist(xRes, breaks = 100)
findGoodDist(x = xRes, distribs = c("norm", "nbinom"))
### nbinom for resistance

#####################
##### Functions #####
#####################

## LRT test
homemadeGtest <- function(full, base){
  dLL = logLik(full) - logLik(base)
  dDF = base$df.residual - full$df.residual
  pvalue <- 1 - stats::pchisq(2*dLL, df=dDF)
  formatC(pvalue, format = "e", digits = 2)
  chisqvalue <- stats::qchisq(p = pvalue, df=dDF)
  return(paste0("G=",round(2*dLL, 1), ", df=", dDF, ", p=", signif(pvalue, digits=2)))
}
## NB. lrtest from pckage lmtest shows similar results ^^ 
## I'm reinventing the wheel again

## LRT test
homemadeGtest <- function(full, base){
  dLL = logLik(full) - logLik(base)
  dDF = base$df.residual - full$df.residual
  pvalue <- 1 - stats::pchisq(2*dLL, df=dDF)
  formatC(pvalue, format = "e", digits = 2)
  chisqvalue <- stats::qchisq(p = pvalue, df=dDF)
  return(paste0("G=",round(2*dLL, 1), ", df=", dDF, ", p=", signif(pvalue, digits=2)))
}

## LRT significance for each factor
myLRTsignificanceFactors <- function(modFull, modPar, modMouse, modInt){
  # print("significance of parasite:")
  return(list(signifParasite = homemadeGtest(modFull, modPar),
              signifMouse = homemadeGtest(modFull, modMouse),
              signifInter = homemadeGtest(modFull, modInt)))
}

testSignif <- function(dataframe, which){
  if(which == "RES"){
    modFULL <- glm.nb(max.OPG ~ infection_isolate*Genotype, data = dataframe)
    modPara <- glm.nb(max.OPG ~ Genotype, data = dataframe)
    modMous <- glm.nb(max.OPG ~ infection_isolate, data = dataframe)
    modinter <- glm.nb(max.OPG ~ infection_isolate+Genotype, data = dataframe)
  } else if (which == "RES_ZI") { # for zero inflated
    modFULL <- zeroinfl(max.OPG ~ infection_isolate*Genotype, data = dataframe, dist = "negbin")
    modPara <- zeroinfl(max.OPG ~ Genotype, data = dataframe, dist = "negbin")
    modMous <- zeroinfl(max.OPG ~ infection_isolate, data = dataframe, dist = "negbin")
    modinter <- zeroinfl(max.OPG ~ infection_isolate+Genotype, data = dataframe, dist = "negbin")
  } else if (which == "IMP"){
    modFULL <- lm(relWL~infection_isolate*Genotype, data = dataframe)
    modPara <- lm(relWL~Genotype, data = dataframe)
    modMous <- lm(relWL~infection_isolate, data = dataframe)
    modinter <- lm(relWL~infection_isolate+Genotype, data = dataframe)
  } else if (which == "TOL"){
    modFULL <- lm(relWL ~ 0 + max.OPG : (infection_isolate * Genotype), data = dataframe)
    modPara <- lm(relWL ~ 0 + max.OPG : (Genotype), data = dataframe)
    modMous <- lm(relWL ~ 0 + max.OPG : (infection_isolate), data = dataframe)
    modinter <- lm(relWL ~ 0 + max.OPG : (infection_isolate + Genotype), data = dataframe)
  }
  return(list(modfull = modFULL, 
              LRT = myLRTsignificanceFactors(modFULL, modPara, modMous, modinter)))
}

testSignifWithinParas <- function(dataframe, which){
  if(which == "RES"){
    modFULL <- glm.nb(max.OPG ~ Genotype, data = dataframe)
    mod0 <- glm.nb(max.OPG ~ 1, data = dataframe)
  } else if (which == "RES_ZI"){
    modFULL <- zeroinfl(max.OPG ~ Genotype, data = dataframe, dist = "negbin")
    mod0 <- zeroinfl(max.OPG ~ 1, data = dataframe, dist = "negbin")
  } else if (which == "IMP"){
    modFULL <- lm(relWL ~ Genotype, data = dataframe)
    mod0 <- lm(relWL ~ 1, data = dataframe)
  } else if (which == "TOL"){
    modFULL <- lm(relWL ~ 0 + max.OPG : Genotype, data = dataframe)
    mod0 <- lm(relWL ~ 0 + max.OPG, data = dataframe)
  }
  G <- homemadeGtest(modFULL, mod0)
  return(list(modfull = modFULL, LRT = G))
}

# Get predicted values for resistance and impact:
getPred <- function(x, which){
  pred <- ggpredict(testSignif(x, which)$modfull, terms = c("Genotype", "infection_isolate"))
  pred <- (data.frame(pred))
  names(pred)[names(pred) %in% c("x", "group")] <- c("Genotype", "infection_isolate")
  # remove misleading predictions for factors with no value
  pred$group <- paste0(pred$Genotype,pred$infection_isolate)
  pred <- pred[pred$group %in% unique(paste0(x$Genotype,x$infection_isolate)),]
  return(pred)
}

# Predicted values of slopes:
getPredTol <- function(x){
  predTolSlopes <- ggpredict(testSignif(x, "TOL")$modfull, terms = c("Genotype", "infection_isolate"),
                             condition = c(max.OPG = 1000000))  ## For a million OPG
  predTolSlopes <- data.frame(predTolSlopes)
  names(predTolSlopes)[names(predTolSlopes) %in% c("x", "group")] <- c("Genotype", "infection_isolate")
  predTolSlopes$group <- paste0(predTolSlopes$Genotype, predTolSlopes$infection_isolate)
  # remove misleading predictions for factors with no value
  predTolSlopes <- predTolSlopes[predTolSlopes$group %in% unique(paste0(x$Genotype,x$infection_isolate)),]
  return(predTolSlopes)
}

## Plots
get_plotR <- function(df, model, cols= mycolors[c(1,2,7,8)]){
  plot_model(model, type = "int", dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
    scale_color_manual(values = cols)+
    scale_y_log10("(predicted) maximum million OPG \n(oocysts per gram of feces)", 
                  breaks = seq(0, 5e6, 0.5e6),
                  labels = as.character(seq(0, 5e6, 0.5e6)/1e6))+
    ggtitle("Maximum parasite load = (inverse of) resistance \n(mean and 95%CI)") +
    xlab("Eimeria isolate") + 
    theme(axis.title.x = element_text(hjust=1), axis.text=element_text(size=13))
} 

get_plotI <- function(df, model, cols=mycolors[c(1,2,7,8)]){
  plot_model(model, type = "int", dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
    scale_color_manual(values = cols)+
    scale_y_continuous(labels = scales::percent_format(accuracy = 5L), 
                       name = "(predicted) maximum weight loss \nrelative to day of infection")+
    ggtitle("Maximum weight loss \n(mean and 95%CI)") +
    xlab("Eimeria isolate") + 
    theme(axis.title.x = element_text(hjust=1), axis.text=element_text(size=13))
} 

get_plotT <- function(df, N, cols=mycolors[c(1,2,7,8)]){ # n is the maximum OPG for end of geom_line
  # make line up to 5e6 OPG for plot
  pts <- getPredTol(df)
  pts$predicted <- pts$predicted*N
  pts$relWL_OPGnull <- 0
  names(pts)[names(pts) %in% c("predicted")] <- "relWL_NMOPG"
  pts <- melt(pts, measure.vars = c("relWL_OPGnull", "relWL_NMOPG"))
  names(pts)[names(pts) %in% c("variable", "value")] <- c("max.OPG", "relWL")
  pts$max.OPG <- as.character(pts$max.OPG)
  pts$max.OPG[pts$max.OPG %in% "relWL_OPGnull"] <- "0"
  pts$max.OPG[pts$max.OPG %in% "relWL_NMOPG"] <- paste0(N, "e6")
  pts$max.OPG <- as.numeric(pts$max.OPG)
  pts$group <- factor(paste0(pts$Mouse_genotype, pts$infection_isolate))
  
  pts$label <- pts$Genotype
  pts$label[pts$max.OPG %in% 0] <- NA
  
  ggplot(pts, aes(x = max.OPG, y = relWL, col = Genotype)) +
    geom_line(aes(group = Genotype)) +
    facet_grid(.~infection_isolate) +
    geom_label(aes(label = label), nudge_x = 0.25e6, na.rm = T)+
    scale_color_manual(values = cols) +
    scale_x_continuous("maximum million oocysts per gram of feces",
                       breaks = seq(0, 4000000, 1000000),
                       labels = seq(0, 4000000, 1000000)/1000000) +
    scale_y_continuous(name = "maximum weight loss \nrelative to day of infection",
                       breaks = seq(0,0.3, 0.05),
                       labels = scales::percent_format(accuracy = 5L)) +
    geom_point(data = df, size = 4, alpha = .5)+
    coord_cartesian(xlim=c(0, N*1000000), ylim=c(0, 0.25)) +
    theme(legend.position = "none") +
    ggtitle("Tolerance \n(slope of B (max weight loss) on A (max parasite load), per genotype)")
}

################################
##### Statistical analyses #####
################################

# to apply on our 2 DF (one more conservative):
MyListDF <- list(full = art2SummaryDF, cons = art2SummaryDF_conservative)
listPar <- list("Brandenburg139 (E. ferrisi)","Brandenburg64 (E. ferrisi)", "Brandenburg88 (E. falciformis)")
names(listPar) <- c("Brandenburg139", "Brandenburg64", "Brandenburg88")

############### 1. Local adaptation of pure strains for E. ferrisi
MyListDF_locad <- MyListDF

MyListDF_locad$full$infection_isolate <- relevel(MyListDF_locad$full$infection_isolate, "Brandenburg64 (E. ferrisi)")
MyListDF_locad$full <- MyListDF_locad$full[!grepl("-", MyListDF_locad$full$Mouse_genotype),]
MyListDF_locad$full <- dropLevelsAllFactorsDF(
  MyListDF_locad$full[grep("ferrisi", MyListDF_locad$full$infection_isolate),])

MyListDF_locad$cons$infection_isolate <- relevel(MyListDF_locad$cons$infection_isolate, "Brandenburg64 (E. ferrisi)")
MyListDF_locad$cons <- MyListDF_locad$cons[!grepl("-", MyListDF_locad$cons$Mouse_genotype),]
MyListDF_locad$cons <- dropLevelsAllFactorsDF(
  MyListDF_locad$cons[grep("ferrisi", MyListDF_locad$cons$infection_isolate),])

### Res
lapply(MyListDF_locad, function(x){testSignif(x,"RES")$LRT}) # interaction factor not significant
### Imp
lapply(MyListDF_locad, function(x){testSignif(x,"IMP")$LRT}) # interaction factor not significant (apart conserv.)
### Tol
lapply(MyListDF_locad, function(x){testSignif(x,"TOL")$LRT}) # interaction factor not significant

## Plot Figure 3:

listPlotRes_LA <- lapply(MyListDF_locad, function(x){
  df = x
  mod = testSignif(x,"RES")$modfull
  get_plotR(df, mod)
})

listPlotImp_LA <- lapply(MyListDF_locad, function(x){
  df = x
  mod = testSignif(x,"IMP")$modfull
  get_plotI(df, mod)
})

listPlotTol_LA <- lapply(MyListDF_locad, function(x){
  df = x
  get_plotT(df, 4)
})

Fig3 <- cowplot::ggdraw() +
  draw_plot(listPlotRes_LA$full + theme(legend.position = "none"), 0, .5, .49, .49) +
  draw_plot(listPlotImp_LA$full + theme(legend.position = "none"), .5, .5, .49, .49) +
  draw_plot(listPlotTol_LA$full+ theme(legend.position = "none"), 0, 0, .7, .49) +
  draw_plot_label(c("A", "B", "C"), c(.01, .51, .01), c(1, 1, .5), size = 15)

pdf(file = "../figures/Fig3_temp.pdf",
    width = 8, height = 8)
Fig3
dev.off()

pdf(file = "../figures/SupplS2_part2_temp.pdf",
    width = 8, height = 8)
cowplot::ggdraw() +
  draw_plot(listPlotRes_LA$cons + theme(legend.position = "none"), 0, .5, .49, .49) +
  draw_plot(listPlotImp_LA$cons + theme(legend.position = "none"), .5, .5, .49, .49) +
  draw_plot(listPlotTol_LA$cons+ theme(legend.position = "none"), 0, 0, .7, .49) +
  draw_plot_label(c("A", "B", "C"), c(.01, .51, .01), c(1, 1, .5), size = 15)
dev.off()

############### 2. E. ferrisi64 and E.fal88: res, impact, tolerance
MyListDF_6488 <- MyListDF

MyListDF_6488$full <- dropLevelsAllFactorsDF(MyListDF_6488$full[!grepl("139", MyListDF_6488$full$infection_isolate),])
MyListDF_6488$cons <- dropLevelsAllFactorsDF(MyListDF_6488$cons[!grepl("139", MyListDF_6488$cons$infection_isolate),])

##### 2.1. within E64
##### 2.2. within E88
parasites <- list("Brandenburg64 (E. ferrisi)", "Brandenburg88 (E. falciformis)")

## RES
lapply(MyListDF_6488, function(xlist){
    testSignifWithinParas(xlist[xlist$infection_isolate %in% parasites[1],], "RES")$LRT})
lapply(MyListDF_6488, function(xlist){
  testSignifWithinParas(xlist[xlist$infection_isolate %in% parasites[2],], "RES_ZI")$LRT})
## diff res in E64 "G=26.6, df=7, p=4e-04" and in E88 "G=28.6, df=14, p=0.012"

#### posthoc:

## IMP
lapply(MyListDF_6488, function(xlist){
  testSignifWithinParas(xlist[xlist$infection_isolate %in% parasites[1],], "IMP")$LRT})
# "G=21.5, df=7, p=0.0031"
lapply(MyListDF_6488, function(xlist){
  testSignifWithinParas(xlist[xlist$infection_isolate %in% parasites[2],], "IMP")$LRT})
# "G=21, df=7, p=0.0038"

#### posthoc:

## TOL
lapply(MyListDF_6488, function(xlist){
  testSignifWithinParas(xlist[xlist$infection_isolate %in% parasites[1],], "TOL")$LRT})
# "G=6.8, df=7, p=0.45"
lapply(MyListDF_6488, function(xlist){
  testSignifWithinParas(xlist[xlist$infection_isolate %in% parasites[2],], "TOL")$LRT})
# "G=13.9, df=7, p=0.054"

#### posthoc:

##### 2.3. between E64 and E88

### Res
lapply(MyListDF_6488, function(x){testSignif(x,"RES")$LRT}) # interaction factor significant (P and M *)
### Imp
lapply(MyListDF_6488, function(x){testSignif(x,"IMP")$LRT}) # interaction factor not significant (P and M *)
### Tol
lapply(MyListDF_locad, function(x){testSignif(x,"TOL")$LRT}) # interaction factor not significant (nothing *)

## Plot Figure 4:
listPlotRes_6488 <- lapply(MyListDF_6488, function(x){
  df = x
  mod = testSignif(x,"RES")$modfull
  get_plotR(df, mod, cols = mycolors)
})

listPlotImp_6488 <- lapply(MyListDF_6488, function(x){
  df = x
  mod = testSignif(x,"IMP")$modfull
  get_plotI(df, mod, cols = mycolors)
})

listPlotTol_6488 <- lapply(MyListDF_6488, function(x){
  df = x
  get_plotT(df, 4.5, cols = mycolors)
})

Fig4 <- cowplot::ggdraw() +
  draw_plot(listPlotRes_6488$full + theme(legend.position = "none"), 0, .5, .49, .49) +
  draw_plot(listPlotImp_6488$full + theme(legend.position = "none"), .5, .5, .49, .49) +
  draw_plot(listPlotTol_6488$full+ theme(legend.position = "none"), 0, 0, .7, .49) +
  draw_plot_label(c("A", "B", "C"), c(.01, .51, .01), c(1, 1, .5), size = 15)
Fig4

pdf(file = "../figures/Fig4_temp.pdf",
    width = 8, height = 8)
Fig4
dev.off()

pdf(file = "../figures/SupplS2_part3_temp.pdf", width = 8, height = 8)
cowplot::ggdraw() +
  draw_plot(listPlotRes_6488$cons + theme(legend.position = "none"), 0, .5, .49, .49) +
  draw_plot(listPlotImp_6488$cons + theme(legend.position = "none"), .5, .5, .49, .49) +
  draw_plot(listPlotTol_6488$cons+ theme(legend.position = "none"), 0, 0, .7, .49) +
  draw_plot_label(c("A", "B", "C"), c(.01, .51, .01), c(1, 1, .5), size = 15)
dev.off()

############### 3. E. ferrisi64 and E.fal88: coupling

# Predicted values of resistance and tolerance per mouse genotype:
getMergeRIT <- function(x, y, z){
  colnames(x) <- gsub("name.", "", colnames(x))
  names(x)[names(x) %in% c("predicted", "conf.low",  "std.error", "conf.high")] <-
    paste(names(x)[names(x) %in% c("predicted", "conf.low",  "std.error", "conf.high")], "OPG", sep = "_")
  colnames(y) <- gsub("name.", "", colnames(y))
  names(y)[names(y) %in% c("predicted", "conf.low",  "std.error", "conf.high")] <-
    paste(names(y)[names(y) %in% c("predicted", "conf.low",  "std.error", "conf.high")], "relWL", sep = "_")
  colnames(z) <- gsub("name.", "", colnames(z))
  names(z)[names(z) %in% c("predicted", "conf.low",  "std.error", "conf.high")] <-
    paste(names(z)[names(z) %in% c("predicted", "conf.low",  "std.error", "conf.high")], "Tol", sep = "_")
  merge(merge(x, y), z)
}

## function to reverse and log10 resistance axis:
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

getBigPlot <- function(x){
  predRes <- getPred(x, "RES")
  predImp <- getPred(x, "IMP")
  predTol <- getPredTol(x)
  finalplotDF <- getMergeRIT(predRes, predImp, predTol)
  finalplotDF
  # test correlations:
  listPar <- list("Brandenburg64 (E. ferrisi)", "Brandenburg88 (E. falciformis)")

  l1 <- lapply(listPar, function(x){
    c <- cor.test(finalplotDF[finalplotDF$infection_isolate %in% x, "predicted_OPG"],
                  finalplotDF[finalplotDF$infection_isolate %in% x, "predicted_relWL"],
                  method="spearman")
    paste0(as.character(c$method), ": ", round(c$estimate, 2), ",\n p-value=",  signif(c$p.value, digits=2))
  })

  addCortext1 <- data.frame(infection_isolate = unlist(listPar),
                            testcor = unlist(l1))

  l2 <- lapply(listPar, function(x){
    c <- cor.test(finalplotDF[finalplotDF$infection_isolate %in% x, "predicted_OPG"],
                  finalplotDF[finalplotDF$infection_isolate %in% x, "predicted_Tol"],
                  method="spearman")
    paste0(as.character(c$method), ": ", round(c$estimate, 2), ",\n p-value=",  signif(c$p.value, digits=2))
  })

  addCortext2 <- data.frame(infection_isolate = unlist(listPar),
                            testcor = unlist(l2))

  ## Plot raw (WL vs OPG) and res-tol

  ## Plot1
  P1 <- ggplot(finalplotDF, aes(x=predicted_OPG, y=predicted_relWL)) +
    geom_errorbarh(aes(xmin = conf.low_OPG, xmax = conf.high_OPG), color = "grey") +
    geom_errorbar(aes(ymin = conf.low_relWL, ymax = conf.high_relWL), color = "grey") +
    geom_point(aes(col = Genotype), size = 7)+
    facet_grid(.~infection_isolate)+
    scale_x_log10(name = "Maximum million oocysts per gram of feces (OPG)",
                  breaks = c(100000, 300000, 500000, 1000000, 2000000, 3000000),
                  labels = c(100000, 300000, 500000, 1000000, 2000000, 3000000)/1000000) +
    scale_y_continuous(name = "Maximum relative weight loss",
                       labels = scales::percent_format(accuracy = 5L))+
    geom_text(aes(label=substring(Genotype, 1, 1)), col = "white")+
    scale_color_manual(values = mycolors) +
    geom_smooth(method = "lm", se = F, col = "black")+
    geom_text(data = addCortext1, aes(label = testcor, x = 1.8e6, y = 0.17))
  P1

  ## Plot2
  P2 <- ggplot(finalplotDF, aes(x = predicted_OPG, y = predicted_Tol)) +
    geom_errorbar(aes(ymin = conf.low_Tol, ymax = conf.high_Tol), color = "grey") +
    geom_errorbarh(aes(xmin = conf.low_OPG, xmax = conf.high_OPG), color = "grey") +
    geom_point(aes(col = Genotype), size = 7)+
    scale_x_continuous(trans=reverselog_trans(10), "RESISTANCE \n(inverse of) maximum million OPG", 
                       breaks = c(100000, 300000, 500000, 1000000, 2000000, 3000000),
                       labels = c(100000, 300000, 500000, 1000000, 2000000, 3000000)/1000000) +
    scale_y_continuous(trans=reverselog_trans(10), name = "TOLERANCE \n(inverse of) slope of maximum weight loss on maximum OPG")+
    facet_grid(.~infection_isolate)+
    geom_text(aes(label=substring(Genotype, 1, 1)), col = "white")+
    coord_cartesian(ylim=c(0.5,0.01))+
    scale_color_manual(values = mycolors) +
    geom_smooth(method = "lm", se = F, col = "black") +
    geom_text(data = addCortext2, aes(label = testcor, x = 1.5e6, y = 0.22))
  P2
  
  bigPlot <- cowplot::ggdraw() +
    draw_plot(P1, 0, .5, 1, .5) +
    draw_plot(P2, 0, 0, 1, .5)
  bigPlot
}

listBigPlot_6488 <- lapply(MyListDF_6488, function(x){
  getBigPlot(x)
})

pdf(file = "../figures/Fig5_temp.pdf",
    width = 8, height = 8)
listBigPlot_6488$full
dev.off()







lapply(MyListDF_F0, function(xlist){
  lapply(listPar[c(1,2)], function(xpar){
    testSignifWithinParas(xlist[xlist$infection_isolate %in% xpar,], "IMP")$LRT})
}) ## diff imp in E64, not in E138

lapply(MyListDF_F0, function(xlist){
  lapply(listPar[c(1,2)], function(xpar){
    testSignifWithinParas(xlist[xlist$infection_isolate %in% xpar,], "TOL")$LRT})
}) ## no diff tol



############### 2.2. E64 vs E88 for all strains

### Res
lapply(MyListDF, function(xlist){
  lapply(listPar, function(xpar){
    testSignifWithinParas(xlist[xlist$infection_isolate %in% xpar,], "RES")$LRT})
})  # consistent. 64 different. bad fit for E88, likely because zeros

reswithinpar <- lapply(MyListDF, function(xlist){
  lapply(listPar, function(xpar){
    testSignifWithinParas(xlist[xlist$infection_isolate %in% xpar,], "RES")})
})
## E.falciformis (4 individuals with zeros)

ZI88 <- testSignifWithinParas(art2SummaryDF[art2SummaryDF$infection_isolate %in% listPar[3],], "RES_ZI")
# best fit? zero inflated far better
lrtest(reswithinpar$full$Brandenburg88$modfull, ZI88$modfull)
ZI88$LRT #G=28.6 ,df=14 ,p=0.012
ZI88

## conservative
ZI88_cons <- testSignifWithinParas(art2SummaryDF_conservative[
  art2SummaryDF_conservative$infection_isolate %in% listPar[3],], "RES_ZI")
# best fit? zero inflated far better
lrtest(reswithinpar$cons$Brandenburg88$modfull, ZI88_cons$modfull)
ZI88_cons$LRT #G=24 ,df=14 ,p=0.046
ZI88_cons

### Imp
lapply(MyListDF, function(xlist){
  lapply(listPar, function(xpar){
    testSignifWithinParas(xlist[xlist$infection_isolate %in% xpar,], "IMP")$LRT})
}) # 139, 64 and 88

### Tol
lapply(MyListDF, function(xlist){
  lapply(listPar, function(xpar){
    testSignifWithinParas(xlist[xlist$infection_isolate %in% xpar,], "TOL")$LRT})
}) # none sauf in cons

##########
## plot ##
##########

############
## Resistance
# plot marginal effects of interaction terms by isolates & strains


testSignif(art2SummaryDF, "RES")$modfull

get_plotR <- function(dataframe){
  plot_model(testSignif(dataframe, "RES")$modfull,
             type = "int", dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
    scale_color_manual(values = mycolors)+
    scale_y_log10("(predicted) maximum million oocysts per gram of feces", 
                  breaks = seq(0, 5e6, 0.5e6),
                  labels = as.character(seq(0, 5e6, 0.5e6)/1e6))+
    ggtitle("Maximum parasite load = (inverse of) resistance \n(mean and 95%CI)") +
    xlab("Eimeria isolate") +
    theme(axis.title.x = element_text(hjust=1), axis.text=element_text(size=13))
} 

plotResList <- lapply(MyListDF, function(x) get_plotR(x))
plotResList[1]
############
## Impact ##
get_plotI <- function(dataframe){
  plot_model(testSignif(dataframe, "IMP")$modfull,
             type = "int",dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
    scale_color_manual(values = mycolors)+
    xlab("Eimeria isolate") +
    ggtitle("Maximum weight loss \n(mean and 95%CI)") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 5L), 
                       name = "(predicted) maximum weight loss compared to day of infection")+
    theme(axis.title.x = element_text(hjust=1), axis.text=element_text(size=13)) #+
}

plotImpList <- lapply(MyListDF, function(x) get_plotI(x))
plotImpList[1]

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
  
  pts$label <- pts$Mouse_genotype
  levels(pts$label) <- as.character(c(1:8))
  pts$Genotype <- paste0(pts$label, ". ", pts$Mouse_genotype)
  pts$label[pts$max.OPG == 0] <- NA_character_
  mydata$label <- mydata$Mouse_genotype
  levels(mydata$label) <- as.character(c(1:8))
  mydata$Genotype <- paste0(mydata$label, ". ", mydata$Mouse_genotype)
  
  ggplot(pts, aes(x = max.OPG, y = relWL, col = Genotype)) +
    geom_line(aes(group = group)) +
    facet_grid(.~infection_isolate) +
    geom_label(aes(label = label), nudge_x = -10, na.rm = TRUE)+
    scale_color_manual(values = mycolors) +
    scale_x_continuous("maximum million oocysts per gram of feces",
                       breaks = seq(0, 5000000, 1000000),
                       labels = seq(0, 5000000, 1000000)/1000000) +
    scale_y_continuous(name = "maximum weight loss compared to day of infection",
                       breaks = seq(0,0.3, 0.05),
                       labels = scales::percent_format(accuracy = 5L)) +
    geom_point(data = mydata, size = 4, alpha = .5)+
    coord_cartesian(ylim=c(0, 0.30)) +
    theme(legend.position = "top") +
    ggtitle("Tolerance \n(slope of B (max weight loss) on A (max parasite load), per genotype)")
}

plotTolList <- list(full = get_plotT(predTolList["full"], art2SummaryDF),
                    cons = get_plotT(predTolList["cons"], art2SummaryDF_conservative))

plotTolList[1]

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

listPlotsRIT <- lapply(c("full", "cons"), getFigRIT)
listPlotsRIT[1]
listPlotsRIT[2]

Fig3 <- listPlotsRIT[[1]]
Fig3
pdf(file = "../figures/Fig3_temp.pdf",
    width = 10, height = 10)
Fig3
dev.off()

pdf(file = "../figures/SupplS2_part2_temp.pdf",
    width = 10, height = 10)
listPlotsRIT[[2]]
dev.off()

bigPlot <- getBigPlot(art2SummaryDF)

pdf(file = "../figures/bigPlot.pdf",
    width = 17, height = 12)
bigPlot
dev.off()

bigPlot
## rm outlier
# bigPlot2 <- getBigPlot(art2SummaryDF[!art2SummaryDF$Mouse_genotype %in% c("PWD", "PWD-SCHUNT", "BUSNA-PWD"),])
# pdf(file = "../figures/bigPlot2.pdf",
#     width = 17, height = 12)
# bigPlot2
# dev.off()

getBigPlot(art2SummaryDF[!grepl("PWD",art2SummaryDF$Mouse_genotype),])
getBigPlot(art2SummaryDF[!grepl("BUSNA",art2SummaryDF$Mouse_genotype),])
getBigPlot(art2SummaryDF[!grepl("STRA",art2SummaryDF$Mouse_genotype),])
getBigPlot(art2SummaryDF[!grepl("SCHUNT",art2SummaryDF$Mouse_genotype),])

# mortality test
table(DSart2$infection_isolate,DSart2$dpi, is.na(DSart2$weight))

# , ,  = TRUE
#                                 0  1  2  3  4  5  6  7  8  9 10 11
# Brandenburg139 (E. ferrisi)     0  0  0  0  0  0  0  0  0  0  0  1
# Brandenburg64 (E. ferrisi)      0  0  0  0  0  0  0  0  0 37 37 37
# Brandenburg88 (E. falciformis)  0  0  0  0  0  0  0  0  0  8 12 13

# Efer: no collection until the end.
# Peak WL median: Efal=9, Efer=5
# Efal: 13 animals died before end
# Q°. is it mouse related?
dfEfal <- DSart2[grepl("falci", DSart2$infection_isolate),]
table(dfEfal$Mouse_genotype[dfEfal$dpi ==1])
# SCHUNT        STRA SCHUNT-STRA  BUSNA-STRA  PWD-SCHUNT   BUSNA-PWD       BUSNA         PWD 
# 6           7           8           8           6           7           7           7 

tab = table(dfEfal$Mouse_genotype[dfEfal$dpi == 11], is.na(dfEfal$weight[dfEfal$dpi == 11]))
tab
chisq.test(tab, simulate.p.value = TRUE)
#               FALSE TRUE
# SCHUNT          6    0
# STRA            7    0
# SCHUNT-STRA     8    0
# BUSNA-STRA      8    0
# PWD-SCHUNT      6    0
# BUSNA-PWD       4    3
# BUSNA           3    4
# PWD             1    6

# PWD and BUSNA smaller, weaker   
# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# X-squared = 31.957, df = NA, p-value = 0.0004998
