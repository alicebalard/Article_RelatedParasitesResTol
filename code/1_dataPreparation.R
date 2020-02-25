### data preparation
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

# Summarize all
ALL_summary <- makeSummaryTable(DF_all)

############################## 
#### PB1. dead before end #### 
##############################

###### what is the overall peak day for each parasite isolate? ######
aggregate(ALL_summary$dpi_max.OPG,
          list(ALL_summary$infection_isolate), 
          function(x) {paste(length(x), median(x), round(sd(x),2))})
aggregate(ALL_summary$dpi_minWeight,
          list(ALL_summary$infection_isolate), 
          function(x) {paste(length(x), median(x), round(sd(x),2))})

# Which animals died before median peak day?
notCompleteFal <- unique(DF_all[is.na(DF_all$weight) &  DF_all$dpi %in% 9 &
                                  DF_all$Eimeria_species %in% "E.falciformis" ,"EH_ID"]) 
notCompleteFer <- unique(DF_all[is.na(DF_all$weight) &  DF_all$dpi %in% 5 &
                                  DF_all$Eimeria_species %in% "E.ferrisi" ,"EH_ID"]) 
notComplete <- c(as.character(notCompleteFal), as.character(notCompleteFer))
# mice died before the end (E88)
notComplete
ALL_summary[ALL_summary$EH_ID %in% notComplete,]

d <- DF_all[DF_all$EH_ID %in% notComplete,]
ggplot(d, aes(dpi, weight, group=EH_ID, col = Eimeria_species)) + 
  geom_point() + geom_line() + facet_grid(.~Exp_ID)
ggplot(d, aes(dpi, oocysts.per.tube, group=EH_ID, col = Eimeria_species)) + 
  geom_point() + geom_line() + facet_grid(.~Exp_ID)
table(ALL_summary[ALL_summary$EH_ID %in% notComplete, "Mouse_genotype"])
# MMm_F1 (Bu-Pw)     MMm_F0 (Bu-Bu)     MMm_F0 (Pw-Pw) 
# 2                  3                  3 

############################################
#### PB2. no oocysts at median peak day #### 
############################################

# All problematic animals: no shedding at all, no oocysts at median OO peak

# How many mice did not shed at all?
noOO <- ALL_summary[ALL_summary$sumoocysts.per.tube == 0,"EH_ID"]

DF_all[DF_all$EH_ID %in% noOO & DF_all$dpi >7,c("EH_ID", "dpi", "weight", "oocysts.per.tube", "Eimeria", "Mouse_strain")] %>%
  dplyr::arrange(EH_ID, dpi)
# 5 mice didn't show Oo but they died before the end of the experiment: NA for oo
# 187,230,244,250,253. 244 died after peak WL

# which mice did not have oocysts collected on peak day?
m1 <- as.character(DF_all[DF_all$Eimeria_species %in% "E.falciformis" &
                      (DF_all$oocysts.per.tube %in% 0 |is.na(DF_all$oocysts.per.tube)) &
                      DF_all$dpi == 8, "EH_ID"])
m2 <- as.character(DF_all[DF_all$Eimeria_species %in% "E.ferrisi" &
                      (DF_all$oocysts.per.tube %in% 0 | is.na(DF_all$oocysts.per.tube)) &
                      DF_all$dpi == 6, "EH_ID"])
ALL_summary[ALL_summary$EH_ID %in% c(m1, m2), c("Mouse_genotype", "infection_isolate")] %>% 
  arrange(Mouse_genotype) %>% 
  dplyr::group_by(Mouse_genotype,infection_isolate) %>%
  dplyr::summarize(sum=n()) 
# 10 animals didn't have collection of oocysts at peak day (4 Bu-Bu E88 + 1 Pw-Pw E88) -> set to NA
ALL_summary$max.OPG[ALL_summary$EH_ID %in% c(m1,m2)] <- NA

###############################################
#### Finalize datasets for article 2 Alice #### 
###############################################
ALL_summary_F0 <- ALL_summary[grep("F0", ALL_summary$Mouse_genotype),]
ALL_summary_F1 <- ALL_summary[grep("F1", ALL_summary$Mouse_genotype),]

# Coinfection to adress somehow. Try at the end with / without second batch
conta <- DF_all[DF_all$dpi %in% 0 & !DF_all$oocysts.per.tube %in% 0 & !is.na(DF_all$oocysts.per.tube), ]
table(conta$infection_isolate)

# full data used for Article 2
art2al_RAWdf <- dropLevelsAllFactorsDF(DF_all[grep("F0", DF_all$Mouse_genotype),])
art2al_SUMdf <- dropLevelsAllFactorsDF(ALL_summary_F0)

length(na.omit(art2al_SUMdf$relWL))# 108 F0 mice
length(na.omit(art2al_SUMdf$max.OPG))# 99 F0 mice

## Set infection group (group 1 had anthelminthics, did not kill worms, stopped after: test effect)
art2al_SUMdf$batch <- art2al_SUMdf$Exp_ID
levels(art2al_SUMdf$batch) <- c("B1", "B2", "B3", "B4", "B5", "B6")

art2al_SUMdf$anth<- FALSE
art2al_SUMdf$anth[art2al_SUMdf$batch %in% "B1"] <- TRUE

### The end ###