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
                      "ggeffects")

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(list.of.packages)

library(sjmisc)
library(sjPlot)

## Reinstall the package in case I updated it
devtools::install_github("alicebalard/parasiteLoad")
library(parasiteLoad)
devtools::install_github("vqv/ggbiplot")
library(ggbiplot)

###### data preparation ######
DF_all<- merge(ExpeDF_003_4, ExpeDF_005, all = T)
DF_all <- makeMiceGenotypeAndIsolate(DF_all)
length(unique(DF_all$EH_ID)) # 168 mice
ALL_summary <- makeSummaryTable(DF_all[DF_all$dpi %in% 4:11,])

############################## 
#### PB1. dead before end #### 
##############################

# Which animals died before median peak day?
notCompleteFal <- unique(DF_all[is.na(DF_all$weight) &  DF_all$dpi %in% 9 &
                                  DF_all$Eimeria_species %in% "E.falciformis" ,"EH_ID"]) 
notCompleteFer <- unique(DF_all[is.na(DF_all$weight) &  DF_all$dpi %in% 5 &
                                  DF_all$Eimeria_species %in% "E.ferrisi" ,"EH_ID"]) 
notComplete <- c(as.character(notCompleteFal), as.character(notCompleteFer)) # N = 12 mice died before the end (E88)
notComplete

d <- DF_all[DF_all$EH_ID %in% notComplete,]
ggplot(d, aes(dpi, weight, group=EH_ID, col = Eimeria_species)) + 
  geom_point() + geom_line() + facet_grid(.~Exp_ID)
ggplot(d, aes(dpi, oocysts.per.tube, group=EH_ID, col = Eimeria_species)) + 
  geom_point() + geom_line() + facet_grid(.~Exp_ID)
table(ALL_summary[ALL_summary$EH_ID %in% notComplete, "Mouse_genotype"])
# MMm_F1 (Bu-Pw)     MMm_F0 (Bu-Bu)     MMm_F0 (Pw-Pw) 
# 2                  4                  6

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
ALL_summary$peak.oocysts.per.g.mouse[ALL_summary$EH_ID %in% c(m1,m2)] <- NA

##########################
#### Finalize dataset #### 
##########################

# Calculate tolerance (inverse)
ALL_summary$invtolerance <- ALL_summary$relWL / ALL_summary$peak.oocysts.per.g.mouse

table(is.na(ALL_summary$invtolerance))

ALL_summary_F0 <- ALL_summary[grep("F0", ALL_summary$Mouse_genotype),]
ALL_summary_F1 <- ALL_summary[grep("F1", ALL_summary$Mouse_genotype),]

# Coinfection to adress somehow. Try at the end with / without second batch
conta <- DF_all[DF_all$dpi %in% 0 & !DF_all$oocysts.per.tube %in% 0 & !is.na(DF_all$oocysts.per.tube), ]
table(conta$infection_isolate)

### The end ###
