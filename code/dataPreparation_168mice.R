ExpeDF_003_4 <- readRDS(file = "../data/ExpeDF_003_4_Alice.rds")
ExpeDF_005 <- readRDS(file = "../data/ExpeDF_005_Alice.rds")
forMap <- read.csv("../data/geolocalisation.csv")
source("myFunctions.R")

# Define a theme
theme_set(theme_bw())

## Packages
list.of.packages <- c("parasiteLoad", "bbmle", "devtools", "optimx", # for bbmle it needs to be required(?)
                      "ggplot2", "VennDiagram","fitdistrplus", # evaluate distribution
                      "epiR", # Sterne's exact method
                      "ggmap", "gridExtra",# several plots in one panel
                      "wesanderson", # nice colors
                      "ggpubr", "tidyr", "stats", "ggrepel",
                      "lme4", "lmerTest", "reshape",
                      "ggeffects", "sjmisc", "sjPlot")

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(list.of.packages)

## Reinstall the package in case I updated it
# devtools::install_github("alicebalard/parasiteLoad")
# library(parasiteLoad)
# devtools::install_github("vqv/ggbiplot")
# library(ggbiplot)

###### data cleaning ######
DF_all<- merge(ExpeDF_003_4, ExpeDF_005, all = T)

# NB. Let's not consider which parent is which, but make A_B mouse = B_A mouse
# we don't have enough individuals to test this effect, and we are not interested in it anyway!
DF_all$Mouse_strain <- as.character(DF_all$Mouse_strain)
x <- strsplit(DF_all$Mouse_strain, "_")
y <- lapply(x, sort)
z <- unlist(lapply(y, FUN = function(x){paste(x, collapse="-")}))
DF_all$Mouse_genotype <- z
rm(x, y, z)

## Order the levels to be more clear in later plots (parents will be low and down 
## on the legend, hybrids in between...)
DF_all$Mouse_genotype <- factor(
  DF_all$Mouse_genotype,
  levels = c("NMRI", "WSB", "WP", "PWD1", "SCHUNT-SCHUNT", 
             "STRA-STRA", "SCHUNT-STRA", "BUSNA-STRA","PWD-SCHUNT",
             "BUSNA-PWD", "BUSNA-BUSNA", "PWD-PWD"),
  labels = c("NMRI", "MMd_F0 (Ws-Ws)", "Mmm-Mmd_Hybrid (WP)", "MMm_F0 (Pw1-Pw1)", "MMd_F0 (Sc-Sc)", 
             "MMd_F0 (St-St)", "MMd_F1 (Sc-St)", "Mmm-Mmd_F1 (Bu-St)", "Mmm-Mmd_F1 (Pw-Sc)",
             "MMm_F1 (Bu-Pw)", "MMm_F0 (Bu-Bu)", "MMm_F0 (Pw-Pw)"))

# add subspecies information based on genotype
DF_all$Mouse_subspecies <- NA
DF_all$Mouse_subspecies[DF_all$Mouse_genotype %in% "NMRI"] <- "NMRI"
DF_all$Mouse_subspecies[DF_all$Mouse_genotype %in% c("MMd_F0 (Ws-Ws)", "MMd_F0 (St-St)", "MMd_F0 (Sc-Sc)", "MMd_F1 (Sc-St)")] <- "M.m.dom"
DF_all$Mouse_subspecies[DF_all$Mouse_genotype %in% c("MMm_F0 (Bu-Bu)", "MMm_F0 (Pw-Pw)", "MMm_F0 (Pw1-Pw1)", "MMm_F1 (Bu-Pw)")] <- "M.m.mus"
DF_all$Mouse_subspecies[DF_all$Mouse_genotype %in% c("Mmm-Mmd_Hybrid (WP)", "Mmm-Mmd_F1 (Bu-St)", "Mmm-Mmd_F1 (Pw-Sc)")] <- "Hybrid_mus_dom"
DF_all$Mouse_subspecies <- factor(DF_all$Mouse_subspecies, levels = c("M.m.dom", "Hybrid_mus_dom", "M.m.mus"))

# add crossing level (F0, F1...) information based on genotype
DF_all$crossingLevel <- NA
DF_all$crossingLevel[DF_all$Mouse_genotype %in% "NMRI"] <- "inbredNMRI"
DF_all$crossingLevel[grep("F0", DF_all$Mouse_genotype)] <- "F0"
DF_all$crossingLevel[grep("F1", DF_all$Mouse_genotype)] <- "F1"

# rename properly infection isolates
DF_all$infection_isolate <- factor(
  DF_all$infection_isolate,
  levels = c("E139", "E64", "E88", "EfLab"),
  labels = c("Brandenburg139 (E. ferrisi)", "Brandenburg64 (E. ferrisi)", 
             "Brandenburg88 (E. falciformis)", "EfLab (E.falciformis)"))

# Put E.ferrisi first (prettier for plots)
DF_all$Eimeria_species <- relevel(as.factor(DF_all$Eimeria_species), "E.ferrisi")

# drop unused levels
DF_all <- data.frame(lapply(DF_all, function(x) if (is.factor(x)) droplevels(x) else x))

# add hybrid index for inbred strains, based on genotype
DF_all$HI <- 0
DF_all$HI[grep("MMm", DF_all$Mouse_genotype)] <- 1
DF_all$HI[grep("Mmm-Mmd", DF_all$Mouse_genotype)] <- 0.5

# Calculate OPG
DF_all <- calculateOPG(DF_all)

# Calculate relative weight loss in individuals
DF_all <- calculateWeightLoss(DF_all)

## First batch (Expe_003) had anthelminthics, did not kill worms, stopped after
DF_all$anth<- FALSE
DF_all$anth[DF_all$Exp_ID %in% "Exp_003"] <- TRUE

#####################################
#### Dataset for article 2 Alice #### 
#####################################
# DSart2 <- DF_all[grep("F0", DF_all$Mouse_genotype),] #reuse if needed
DSart2 <- DF_all
DSart2 <- dropLevelsAllFactorsDF(DSart2)

# rename batches
DSart2$Batch[DSart2$Exp_ID %in% "Exp_003"] <- "B1"
DSart2$Batch[DSart2$Exp_ID %in% "Exp_004"] <- "B2"
DSart2$Batch[grep("Exp_005_1", DSart2$Exp_ID)] <- "B3"
DSart2$Batch[grep("Exp_005_2", DSart2$Exp_ID)] <- "B4"

# rename (shorter) mouse strains
# levels(DSart2$Mouse_genotype) <- c("SCHUNT", "STRA", "BUSNA", "PWD")
levels(DSart2$Mouse_genotype) <- factor(levels(DSart2$Mouse_genotype), 
       levels = c("MMd_F0 (Sc-Sc)","MMd_F0 (St-St)", "MMd_F1 (Sc-St)", "Mmm-Mmd_F1 (Bu-St)", "Mmm-Mmd_F1 (Pw-Sc)",
                  "MMm_F1 (Bu-Pw)", "MMm_F0 (Bu-Bu)", "MMm_F0 (Pw-Pw)"),
       labels = c("SCHUNT", "STRA", "SCHUNT-STRA",  "STRA-BUSNA", "SCHUNT-PWD",
                  "PWD-BUSNA", "BUSNA", "PWD"))

# 9 animals lost more than 20%, because they died overnight. 
# 8 were between 18 and 20% but got better
pb <- DSart2$EH_ID[DSart2$relWL >0.18 & !is.na(DSart2$relWL)]
ggplot(DSart2[DSart2$EH_ID %in% pb,], aes(x=dpi, y=relWL, group = EH_ID, col=EH_ID)) + 
  geom_point() + geom_line() + geom_hline(yintercept = 0.18, col="red") +
  facet_grid(.~infection_isolate)
# For the 9, we will remove the last weight point, as it's the weight of a dehydrated cadaver
diedOvernight <- DSart2$EH_ID[DSart2$relWL >0.20 & !is.na(DSart2$relWL)]
DSart2$relWL[
  DSart2$EH_ID %in% diedOvernight & DSart2$relWL > 0.20 & !is.na(DSart2$relWL)] <-NA

# Summarize all
art2SummaryDF <- makeSummaryTable(DSart2)

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
# 1    Brandenburg139 (E. ferrisi) 25 5 2.14
# 2     Brandenburg64 (E. ferrisi) 87 5 1.69
# 3 Brandenburg88 (E. falciformis)  56 9 1.5

#########################################################
#### which died before end? or no oocysts collected? #### 

###### For OPG (resistance measure):
## Eferrisi: did everyone shed at dpi 6?
# mice shedding 0 OPG at dpi 6? 
DSart2[DSart2$Eimeria_species %in% "E.ferrisi" & DSart2$dpi ==6 & 
         (DSart2$OPG ==0 | is.na(DSart2$OPG)),  c("EH_ID", "infection_isolate", "Mouse_genotype")] 
# LM0164 Brandenburg139 (E. ferrisi) MMd_F0 (Sc-Sc) # 1 problematic individual, liquid feces
ggplot(DSart2[DSart2$EH_ID %in% "LM0164",], aes(dpi, OPG, group=EH_ID)) + geom_line() 

## Efalciformis: did everyone shed at dpi 8?
d <- DSart2[DSart2$Eimeria_species %in% "E.falciformis" & DSart2$dpi ==8 & 
              (DSart2$OPG ==0 | is.na(DSart2$OPG)),  c("EH_ID", "infection_isolate", "Mouse_genotype")] 
d
# EH_ID              infection_isolate Mouse_genotype
# LM0187 Brandenburg88 (E. falciformis) MMm_F0 (Pw-Pw) -> no sign of shedding, last dpi7 *
# LM0193 Brandenburg88 (E. falciformis) MMm_F0 (Bu-Bu) -> shed only at dpi7 *
# LM0202 Brandenburg88 (E. falciformis) MMm_F0 (Bu-Bu) -> shed at dpi 7, then 9, 10, 11
# LM0230 Brandenburg88 (E. falciformis) MMm_F0 (Bu-Bu) -> no sign of shedding, last dpi7 *
# LM0237 Brandenburg88 (E. falciformis) MMm_F0 (Pw-Pw) -> shed at dpi7, 0 at dpi8, nothing after
# LM0244 Brandenburg88 (E. falciformis) MMm_F0 (Bu-Bu) -> no sign of shedding, last dpi8 (0)
# LM0250 Brandenburg88 (E. falciformis)      BUSNA-PWD -> no sign of shedding, last dpi8 (0)
# LM0253 Brandenburg88 (E. falciformis) MMm_F0 (Pw-Pw) -> no sign of shedding, last dpi8 (0) *
# LM0255 Brandenburg88 (E. falciformis) MMm_F0 (Bu-Bu) -> shed at dpi 7, then 10, 11; nothing between
ggplot(DSart2[DSart2$EH_ID %in% d$EH_ID[9],], aes(dpi, OPG, group=EH_ID, col = EH_ID)) + 
  geom_point() + geom_line() 

###### For weight:
## Eferrisi: did everyone had weight measured at dpi 5?
DSart2[DSart2$Eimeria_species %in% "E.ferrisi" & DSart2$dpi == 5 & 
         (DSart2$weight ==0 | is.na(DSart2$weight)),  c("EH_ID", "infection_isolate", "Mouse_genotype")] 
# NO PB HERE

## Efalciformis: did everyone had weight measured at dpi 9?
d2 <- DSart2[DSart2$Eimeria_species %in% "E.falciformis" & DSart2$dpi == 9 & 
               (DSart2$weight ==0 | is.na(DSart2$weight)),  c("EH_ID", "infection_isolate", "Mouse_genotype")] 
d2 # died before peak
# EH_ID              infection_isolate Mouse_genotype
# LM0187 Brandenburg88 (E. falciformis)            PWD
# LM0193 Brandenburg88 (E. falciformis)          BUSNA
# LM0230 Brandenburg88 (E. falciformis)          BUSNA
# LM0243 Brandenburg88 (E. falciformis)      BUSNA-PWD
# LM0245 Brandenburg88 (E. falciformis)          BUSNA
# LM0250 Brandenburg88 (E. falciformis)      BUSNA-PWD
# LM0252 Brandenburg88 (E. falciformis)            PWD
# LM0253 Brandenburg88 (E. falciformis)            PWD
ggplot(DSart2[DSart2$EH_ID %in% d2$EH_ID[8],], aes(dpi, relWL, group=EH_ID, col = EH_ID)) + 
  geom_point() + geom_line() 

## CCL 
# 6 mice died at peak day (9) Efal, last weight considered
# 8 Efal + 1 Efer mice didn't shed at peak day, remove them in conservative dataset
miceProblematicNoOo <- c(as.character(d$EH_ID),"LM0164")

##### Third pb: anthelminthic & contamination
conta <- DSart2[DSart2$dpi %in% 0 & !DSart2$OPG %in% 0 & !is.na(DSart2$OPG), ]
table(conta$infection_isolate)

############ Make all datasets:
# FULL = DSart2 / art2SummaryDF

# conservative = remove mice with contamination or anthelminthic
DSart2_conservative <- DSart2[DSart2$anth == F & !DSart2$EH_ID %in% conta$EH_ID,]
art2SummaryDF_conservative <- makeSummaryTable(DSart2_conservative) # 118 mice

### The end ###