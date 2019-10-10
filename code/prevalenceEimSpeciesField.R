# Analysis of prevalence Eimeria 2014 to 2017
library(ggplot2)
prevDF <- readRDS("../data/speciesIdentField1417.rds")
hiDF <- readRDS("../data/MiceTable_fullEimeriaInfos_2014to2017.rds")

# manual correction
prevDF$Mouse_ID <- as.character(prevDF$Mouse_ID)
prevDF$Mouse_ID <- gsub(" ", "", prevDF$Mouse_ID)

# 3 double infection, all ferrisi + vermiformis
prevDF$Species <- as.character(prevDF$Species)
prevDF[prevDF$Mouse_ID %in% "AA_0100","Species"] <- c("Double (E_ferrisi, E_vermiformis)")
prevDF[grep("Double", prevDF$Species),"Species"] <- c("Double (E_ferrisi, E_vermiformis)")

## before AA_0047, it was a unique farm with a lot of animals caught, we remove to not bias 
toExclude <- c(paste0("A_000", 1:3), paste0("AA_000", 1:9), paste0("AA_00", 10:46))
prevDF <- prevDF[!prevDF$Mouse_ID %in% toExclude,]
hiDF <- hiDF[!hiDF$Mouse_ID %in% toExclude,]

# Merge all info
fieldDF <- merge(prevDF, hiDF[c("HI","Mouse_ID", "Longitude", "Latitude")])

# remove mice without HI (N = 136)
fieldDF <- fieldDF[!is.na(fieldDF$HI),]

# How many mice? N = 1046
length(unique(fieldDF$Mouse_ID))

# plot
ggplot(fieldDF, aes(x = Species, y = HI, col = Species)) +
  geom_violin() +
  geom_jitter(width = 0.05) + theme_bw()

# calculate prevalence (rough cut in the middle)
fieldDF$MusSubSp[fieldDF$HI < 0.5] <- "Mmd"
fieldDF$MusSubSp[fieldDF$HI >= 0.5] <- "Mmm"

# by halves
library(dplyr)
D1 <- fieldDF %>%
  dplyr::group_by(MusSubSp) %>%
  dplyr::summarize(n_halfHZ = n()) %>% 
  as.data.frame()

D2 <- fieldDF %>%
  dplyr::group_by(MusSubSp, Species) %>%
  dplyr::summarize(n = n()) %>%
  as.data.frame()

prevalenceDF <- merge(D1, D2)
prevalenceDF$prevalence <- round(prevalenceDF$n / prevalenceDF$n_halfHZ *100, 2)

prevalenceDF

library(tidyr)
spread(prevalenceDF[c("MusSubSp", "prevalence", "Species")], MusSubSp, prevalence)

# check y farm that there is no bias
fieldDF$locality <- paste(fieldDF$Latitude, fieldDF$Longitude, sep = "_")

# define by what is a locality infected

# library(plyr)
# fieldDFLOC <- plyr::ddply(fieldDF, "locality", head, 1)

D3 <- fieldDFLOC %>%
  dplyr::group_by(MusSubSp) %>%
  dplyr::summarize(n_halfHZ = n()) %>% 
  as.data.frame()

D4 <- fieldDFLOC %>%
  dplyr::group_by(MusSubSp, Species) %>%
  dplyr::summarize(n = n()) %>%
  as.data.frame()

prevalenceDFLOC <- merge(D3, D4)
prevalenceDFLOC$prevalence <- round(prevalenceDFLOC$n / prevalenceDFLOC$n_halfHZ *100, 2)

prevalenceDFLOC

library(tidyr)
spread(prevalenceDF[c("MusSubSp", "prevalence", "Species")], MusSubSp, prevalence)

