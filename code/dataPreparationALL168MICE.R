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
