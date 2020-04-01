### data preparation all 168 mice
source("dataPreparationALL168MICE.R")

#####################################
#### Dataset for article 2 Alice #### 
#####################################
DSart2 <- DF_all[grep("F0", DF_all$Mouse_genotype),]
DSart2 <- dropLevelsAllFactorsDF(DSart2)

# rename batches
DSart2$Batch <- DSart2$Exp_ID
levels(DSart2$Batch) <- c("B1","B2","B3","B4","B5","B6")

# rename (shorter) mouse strains
levels(DSart2$Mouse_genotype) <- c("SCHUNT", "STRA", "BUSNA", "PWD")

# Summarize all
art2SummaryDF <- makeSummaryTable(DSart2)

###### what is the overall peak day for each parasite isolate? ######
aggregate(art2SummaryDF$dpi_max.OPG,
          list(art2SummaryDF$infection_isolate), 
          function(x) {paste(length(x), median(x), round(sd(x),2))})
# Brandenburg139 (E. ferrisi) 25 6 0.73
# Brandenburg64 (E. ferrisi) 56 6 1.01
# Brandenburg88 (E. falciformis) 27 8 1.78
aggregate(art2SummaryDF$dpi_minWeight,
          list(art2SummaryDF$infection_isolate), 
          function(x) {paste(length(x), median(x), round(sd(x),2))})
# Brandenburg139 (E. ferrisi) 25 5 2.14
# Brandenburg64 (E. ferrisi) 56 5 1.92
# Brandenburg88 (E. falciformis) 27 9 1.49

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
# LM0253 Brandenburg88 (E. falciformis) MMm_F0 (Pw-Pw) -> no sign of shedding, last dpi8 (0) *
# LM0255 Brandenburg88 (E. falciformis) MMm_F0 (Bu-Bu) -> shed at dpi 7, then 10, 11; nothing between
ggplot(DSart2[DSart2$EH_ID %in% d$EH_ID[8],], aes(dpi, OPG, group=EH_ID, col = EH_ID)) + 
  geom_point() + geom_line() 

###### For weight:
## Eferrisi: did everyone had weight measured at dpi 5?
DSart2[DSart2$Eimeria_species %in% "E.ferrisi" & DSart2$dpi == 5 & 
         (DSart2$weight ==0 | is.na(DSart2$weight)),  c("EH_ID", "infection_isolate", "Mouse_genotype")] 
# NO PB HERE

## Efalciformis: did everyone had weight measured at dpi 9?
d2 <- DSart2[DSart2$Eimeria_species %in% "E.falciformis" & DSart2$dpi == 9 & 
         (DSart2$weight ==0 | is.na(DSart2$weight)),  c("EH_ID", "infection_isolate", "Mouse_genotype")] 
# LM0187 Brandenburg88 (E. falciformis) MMm_F0 (Pw-Pw) <- >20% WL at dpi8, sacrificed or dead
# LM0193 Brandenburg88 (E. falciformis) MMm_F0 (Bu-Bu) <- >20% WL at dpi8, sacrificed or dead
# LM0230 Brandenburg88 (E. falciformis) MMm_F0 (Bu-Bu) <- >20% WL at dpi8, sacrificed or dead
# LM0245 Brandenburg88 (E. falciformis) MMm_F0 (Bu-Bu) <- >20% WL at dpi8, sacrificed or dead
# LM0252 Brandenburg88 (E. falciformis) MMm_F0 (Pw-Pw) <- >15% WL at dpi8, dead
# LM0253 Brandenburg88 (E. falciformis) MMm_F0 (Pw-Pw) <- >12% WL at dpi8, dead
ggplot(DSart2[DSart2$EH_ID %in% d2$EH_ID[6],], aes(dpi, relWL, group=EH_ID, col = EH_ID)) + 
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

# conservative 1 = remove mice without oocysts at peak day
DSart2_conservative1 <- DSart2[!DSart2$EH_ID %in% miceProblematicNoOo,] # 108 mice
art2SummaryDF_conservative1 <- makeSummaryTable(DSart2_conservative1) # 99 mice

# conservative 2 = remove mice with contamination or anthelminthic
DSart2_conservative2 <- DSart2_conservative1[DSart2_conservative1$anth == F & !DSart2_conservative1$EH_ID %in% conta$EH_ID,]
art2SummaryDF_conservative2 <- makeSummaryTable(DSart2_conservative2) # 68 mice

### The end ###