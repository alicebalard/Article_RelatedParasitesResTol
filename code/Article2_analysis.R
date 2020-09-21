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
# mycolors <- scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=8))
# New colors, scales too difficult to distinguish
mycolors <- c("blue", "blue", "blue","purple","purple", "red", "red", "red")

## Different datasets as follow:

# FULL = DSart2 / art2SummaryDF
nrow(art2SummaryDF) # 168
##### NEW remove E139
art2SummaryDF <- art2SummaryDF[!art2SummaryDF$infection_isolate %in% "Brandenburg139 (E. ferrisi)",]
nrow(art2SummaryDF)

# conservative = remove mice with contamination or anthelminthic
# DSart2_conservative2 ; art2SummaryDF_conservative2 # 118 mice
nrow(art2SummaryDF_conservative)
##### NEW remove E139
art2SummaryDF_conservative <- art2SummaryDF_conservative[!art2SummaryDF_conservative$infection_isolate %in% "Brandenburg139 (E. ferrisi)",]
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

# to pretty plot
d <- aggregate(as.numeric(art2SummaryDF$ageAtInfection),
               list(Mouse_genotype = art2SummaryDF$Mouse_genotype), mean)
names(d)[2] <- "mean_age_genotype"
df <- merge(art2SummaryDF, d)
df$Mouse_genotype <- reorder(df$Mouse_genotype, as.numeric(df$mean_age_genotype))

avg <- df %>%
    summarize(avg = mean(mean_age_genotype, na.rm = T)) %>%
    pull(avg)

unique(data.frame(df$Mouse_genotype, df$mean_age_genotype))

ggplot(df, aes(Mouse_genotype, as.numeric(ageAtInfection), col = Mouse_genotype)) +
    geom_segment(aes(x = Mouse_genotype, xend = Mouse_genotype,
                     y = avg, yend = mean_age_genotype),
                 size = 0.8) +
    geom_hline(aes(yintercept = avg), color = "gray70", size = 0.6) +
    
    coord_flip() +
    geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
    geom_point(aes(Mouse_genotype, as.numeric(mean_age_genotype)), size = 5) 
    
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

## How many mice died before the end?
table(DSart2$infection_isolate, DSart2$dpi, is.na(DSart2$weight))
#   , ,  = FALSE                  0  1  2  3  4  5  6  7  8  9 10 11
# Brandenburg139 (E. ferrisi)    25 25 25 25 25 25 25 25 25 25 24 24
# Brandenburg64 (E. ferrisi)     87 87 87 87 87 87 87 87 87 50 50 28
# Brandenburg88 (E. falciformis) 56 56 56 56 56 56 56 56 51 46 43 43
# , ,  = TRUE                     0  1  2  3  4  5  6  7  8  9 10 11
# Brandenburg139 (E. ferrisi)     0  0  0  0  0  0  0  0  0  0  1  1
# Brandenburg64 (E. ferrisi)      0  0  0  0  0  0  0  0  0 37 37 37
# Brandenburg88 (E. falciformis)  0  0  0  0  0  0  0  0  5 10 13 13

# Peak WL median: Efal=9, Efer=5
# Efer139: 2 died before the end (>3days after the peak)
# Efer64: no one died before the end / no collection until the end for one batch (Batch 4) stop at dpi8
dfE64 <- DSart2[grepl("64", DSart2$infection_isolate),]
table(dfE64$Batch[dfE64$dpi ==9], dfE64$Mouse_genotype[dfE64$dpi ==9], is.na(dfE64$weight[dfE64$dpi ==9]))

# Efal: 13 animals died before end
# Q°. is it mouse related?
dfEfal <- DSart2[grepl("falci", DSart2$infection_isolate),]
table(dfEfal$Mouse_genotype[dfEfal$dpi ==1])
# SCHUNT        STRA SCHUNT-STRA  BUSNA-STRA  PWD-SCHUNT   BUSNA-PWD       BUSNA         PWD 
# 6           7           8           8           6           7           7           7 
tab = table(dfEfal$Mouse_genotype[dfEfal$dpi == 11], is.na(dfEfal$weight[dfEfal$dpi == 11]))
tab
write.csv(tab, "../figures/Table2_temp.csv")

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

## LRT significance for each factor
myLRTsignificanceFactors <- function(modFull, modPar, modMouse, modInt){
        # print("significance of parasite:")
        return(list(signifParasite = lrtest(modFull, modPar),
                    signifMouse = lrtest(modFull, modMouse),
                    signifInter = lrtest(modFull, modInt)))
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

testSignifSUBSPECIES <- function(dataframe, which){
    if(which == "RES"){
        modFULL <- glm.nb(max.OPG ~ infection_isolate*Mouse_subspecies, data = dataframe)
        modPara <- glm.nb(max.OPG ~ Mouse_subspecies, data = dataframe)
        modMous <- glm.nb(max.OPG ~ infection_isolate, data = dataframe)
        modinter <- glm.nb(max.OPG ~ infection_isolate+Mouse_subspecies, data = dataframe)
    } else if (which == "RES_ZI") { # for zero inflated
        modFULL <- zeroinfl(max.OPG ~ infection_isolate*Mouse_subspecies, data = dataframe, dist = "negbin")
        modPara <- zeroinfl(max.OPG ~ Mouse_subspecies, data = dataframe, dist = "negbin")
        modMous <- zeroinfl(max.OPG ~ infection_isolate, data = dataframe, dist = "negbin")
        modinter <- zeroinfl(max.OPG ~ infection_isolate+Mouse_subspecies, data = dataframe, dist = "negbin")
    } else if (which == "IMP"){
        modFULL <- lm(relWL~infection_isolate*Mouse_subspecies, data = dataframe)
        modPara <- lm(relWL~Mouse_subspecies, data = dataframe)
        modMous <- lm(relWL~infection_isolate, data = dataframe)
        modinter <- lm(relWL~infection_isolate+Mouse_subspecies, data = dataframe)
    } else if (which == "TOL"){
        modFULL <- lm(relWL ~ 0 + max.OPG : (infection_isolate * Mouse_subspecies), data = dataframe)
        modPara <- lm(relWL ~ 0 + max.OPG : (Mouse_subspecies), data = dataframe)
        modMous <- lm(relWL ~ 0 + max.OPG : (infection_isolate), data = dataframe)
        modinter <- lm(relWL ~ 0 + max.OPG : (infection_isolate + Mouse_subspecies), data = dataframe)
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
        G <- lrtest(modFULL, mod0)
        return(list(modfull = modFULL, LRT = G))
}

# Get predicted values for resistance, impact or tolerance:
getPredictedValues <- function(x, which, oneortwo){ # one parasite at a time or 2?
        if(oneortwo == 1){
                ## check possible error
                if((length(unique(x$infection_isolate)) !=1) == TRUE){
                        print("Whoo WOOOO Whoo more than 1 parasites isolate!!!")
                }
                if (which == "TOL"){  
                        pred <- ggpredict(testSignifWithinParas(x, "TOL")$modfull, terms = c("Genotype"),
                                          condition = c(max.OPG = 1000000))  ## For a million OPG
                        pred <- (data.frame(pred))
                } else{ 
                        pred <- ggpredict(testSignifWithinParas(x, which)$modfull)
                        pred <- (data.frame(pred$Genotype))
                }
                names(pred)[names(pred) %in% "x"] <- "Genotype"
        }
        if(oneortwo == 2){
                ## check possible error
                if((length(unique(x$infection_isolate)) !=2) == TRUE){
                        print("Whoo WOOOO Whoo that's not just 2 parasites isolates!!!")
                }
                if (which == "TOL"){ 
                        pred <- ggpredict(testSignif(x, "TOL")$modfull, terms = c("Genotype", "infection_isolate"),
                                          condition = c(max.OPG = 1000000))  ## For a million OPG
                } else{ 
                        pred <- ggpredict(testSignif(x, which)$modfull, terms = c("Genotype", "infection_isolate"))
                }
                pred <- (data.frame(pred))
                names(pred)[names(pred) %in% c("x", "group")] <- c("Genotype", "infection_isolate")
                # remove misleading predictions for factors with no value
                pred$group <- paste0(pred$Genotype,pred$infection_isolate)
                pred <- pred[pred$group %in% unique(paste0(x$Genotype,x$infection_isolate)),]
        }
        return(pred)
}

## Plots function

get_plot_par <- function(df, npara, cols, plottype, N = 4, linetype, pointshape){ # n is the maximum OPG for end of geom_line
        plotdf <- getPredictedValues(df, plottype, oneortwo=npara)
        if(plottype != "TOL"){
                plot <- ggplot(plotdf, aes(x=Genotype, y=predicted, col=Genotype))+
                        geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width=.1) +
                        geom_point(size = 5, aes(pch = Genotype)) +
                        scale_color_manual(values = cols)+
                        scale_shape_manual(values = pointshape) +
                        theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
                if(plottype == "RES" | plottype == "RES_ZI"){
                        plot <- plot + scale_y_log10("(predicted) maximum million OPG \n(oocysts per gram of feces)",
                                                     breaks = seq(0, 5e6, 0.5e6),
                                                     labels = as.character(seq(0, 5e6, 0.5e6)/1e6))+
                                ggtitle("Maximum parasite load \n=(inverse of) resistance")
                } else if (plottype == "IMP"){
                        plot <- plot + scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                                          name = "(predicted) maximum weight loss \nrelative to day of infection")+
                                ggtitle("Maximum weight loss")
                }
        }
        if (plottype == "TOL"){
                pts <- getPredictedValues(df, which = "TOL", oneortwo = npara)
                pts$predicted <- pts$predicted*N
                pts$relWL_OPGnull <- 0
                names(pts)[names(pts) %in% c("predicted")] <- "relWL_NMOPG"
                pts <- melt(pts, measure.vars = c("relWL_OPGnull", "relWL_NMOPG"))
                names(pts)[names(pts) %in% c("variable", "value")] <- c("max.OPG", "relWL")
                pts$max.OPG <- as.character(pts$max.OPG)
                pts$max.OPG[pts$max.OPG %in% "relWL_OPGnull"] <- "0"
                pts$max.OPG[pts$max.OPG %in% "relWL_NMOPG"] <- paste0(N, "e6")
                pts$max.OPG <- as.numeric(pts$max.OPG)
                pts$group <- factor(paste0(pts$Genotype, pts$infection_isolate))
                
                pts$label <- pts$Genotype
                pts$label[pts$max.OPG %in% 0] <- NA
                
                plot <- ggplot(pts, aes(x = max.OPG, y = relWL, col = Genotype)) +
                        geom_line(aes(group = Genotype)) +
                        geom_label(aes(label = substring(label, 1, 1)), na.rm = T)+
                        scale_color_manual(values = cols) +
                        scale_x_continuous("maximum million OPG \n(oocysts per gram of feces)",
                                           breaks = seq(0, 4000000, 1000000),
                                           labels = seq(0, 4000000, 1000000)/1000000) +
                        scale_y_continuous(name = "maximum weight loss \nrelative to day of infection",
                                           breaks = seq(0,0.3, 0.05),
                                           labels = scales::percent_format(accuracy = 5L)) +
                        geom_point(data = df, size = 4, alpha = .5, aes(pch=Genotype))+
                        scale_shape_manual(values = pointshape) +
                        coord_cartesian(xlim=c(1, N*1000000), ylim=c(0, 0.25)) +
                        ggtitle("Tolerance \n(slope of B (max weight loss) on A (max parasite load), per genotype)")
        }
        plot <- plot + theme(text = element_text(size=15), plot.title = element_text(size=15))
        if (npara == 2){ plot <- plot+  facet_grid(.~infection_isolate) }
        return(plot)
}

getCouplingPlot <- function(x, isZINB = F, posx1, posy1, posx2, posy2, mycolors, myshapes){ ## NB. for RES Efalci, ZINB model
        # make dataframe for plot
        if(isZINB==TRUE){
                predRes <- getPredictedValues(x, "RES_ZI",1)
        } else {
                predRes <- getPredictedValues(x, "RES", 1)
        }
        predImp <- getPredictedValues(x, "IMP", 1)
        predTol <- getPredictedValues(x, "TOL", 1)
        names(predRes)[names(predRes) %in% c("predicted", "conf.low",  "std.error", "conf.high")] <-
                paste(names(predRes)[names(predRes) %in% c("predicted", "conf.low",  "std.error", "conf.high")], "OPG", sep = "_")
        names(predImp)[names(predImp) %in% c("predicted", "conf.low",  "std.error", "conf.high")] <-
                paste(names(predImp)[names(predImp) %in% c("predicted", "conf.low",  "std.error", "conf.high")], "relWL", sep = "_")
        names(predTol)[names(predTol) %in% c("predicted", "conf.low",  "std.error", "conf.high")] <-
                paste(names(predTol)[names(predTol) %in% c("predicted", "conf.low",  "std.error", "conf.high")], "Tol", sep = "_")
        finalplotDF <- merge(merge(predRes, predImp), predTol, by = "Genotype")
        # test correlations:
        c1 <- cor.test(finalplotDF$predicted_OPG, finalplotDF$predicted_relWL, method="spearman")
        testcor1 <- paste0(as.character(c1$method), ": ", signif(c1$estimate, digits=2), "\n p-value=", signif(c1$p.value, digits=2))
        c2 <- cor.test(finalplotDF$predicted_OPG, finalplotDF$predicted_Tol, method="spearman")
        testcor2 <- paste0(as.character(c2$method), ": ", signif(c2$estimate, digits=2), "\n p-value=", signif(c2$p.value, digits=2))
        # Plot
        ## Plot1
        P1 <- ggplot(finalplotDF, aes(x=predicted_OPG, y=predicted_relWL)) +
                geom_smooth(method = "lm", se = F, col = "black")+
                geom_errorbarh(aes(xmin = conf.low_OPG, xmax = conf.high_OPG), color = "grey") +
                geom_errorbar(aes(ymin = conf.low_relWL, ymax = conf.high_relWL), color = "grey") +
                geom_point(aes(col = Genotype), size = 7, pch = 22, fill = "white")+
                scale_x_continuous(name = "Maximum million OPG",
                                   breaks = seq(0.5,5,0.5)*1e6,
                                   labels = seq(0.5,5,0.5)) +
                scale_y_continuous(name = "Maximum relative weight loss",
                                   labels = scales::percent_format(accuracy = 1))+
                geom_text(aes(label=substring(Genotype, 1, 1), col = Genotype))+
                scale_color_manual(values = mycolors) +
                scale_shape_manual(values = myshapes) +
                annotate("text", x = posx1, y = posy1, label = testcor1)+
                theme(legend.position = "none")+ theme(text = element_text(size=15))
        
        P2 <- ggplot(finalplotDF, aes(x = predicted_OPG, y = -predicted_Tol)) +
                geom_smooth(method = "lm", se = F, col = "black") +
                geom_errorbar(aes(ymin = -conf.low_Tol, ymax = -conf.high_Tol), color = "grey") +
                geom_errorbarh(aes(xmin = conf.low_OPG, xmax = conf.high_OPG), color = "grey") +
                geom_point(aes(col = Genotype), size = 7, pch = 22, fill = "white")+
                scale_x_continuous("Maximum million OPG \n (inverse of) RESISTANCE ",
                                               breaks = seq(0.5,5,0.5)*1e6,
                                               labels = seq(0.5,5,0.5))+
                scale_y_continuous(name = "TOLERANCE (inverse of) slope of\n maximum weight loss on maximum OPG ")+
                geom_text(aes(label=substring(Genotype, 1, 1), col = Genotype))+
                scale_color_manual(values = mycolors) +
                scale_shape_manual(values = myshapes) +
                annotate("text", x = posx2, y = posy2, label = testcor2)+
                theme(legend.position = "none")+ theme(text = element_text(size=15))
        return(list(P1, P2))
}

################################
##### Statistical analyses #####
################################

# to apply on our 2 DF (one more conservative):
MyListDF <- list(full = art2SummaryDF, cons = art2SummaryDF_conservative)
listPar <- list("Brandenburg139 (E. ferrisi)","Brandenburg64 (E. ferrisi)", "Brandenburg88 (E. falciformis)")
names(listPar) <- c("Brandenburg139", "Brandenburg64", "Brandenburg88")

############### 1. E. ferrisi64 and E.fal88: res, impact, tolerance
MyListDF_64 <- MyListDF
MyListDF_64$full <- dropLevelsAllFactorsDF(MyListDF_64$full[grepl("64", MyListDF_64$full$infection_isolate),])
MyListDF_64$cons <- dropLevelsAllFactorsDF(MyListDF_64$cons[grepl("64", MyListDF_64$cons$infection_isolate),])

MyListDF_88 <- MyListDF
MyListDF_88$full <- dropLevelsAllFactorsDF(MyListDF_88$full[grepl("88", MyListDF_88$full$infection_isolate),])
MyListDF_88$cons <- dropLevelsAllFactorsDF(MyListDF_88$cons[grepl("88", MyListDF_88$cons$infection_isolate),])

MyListDF_64_88 <- MyListDF
MyListDF_64_88$full <- rbind(MyListDF_64$full, MyListDF_88$full)
MyListDF_64_88$cons <- rbind(MyListDF_64$cons, MyListDF_88$cons)

## test factor significance
# RES
# nbin or zeroinfl
mod1 <- glm.nb(max.OPG ~ infection_isolate*Genotype, data = MyListDF_64_88$full)
mod2 <- zeroinfl(max.OPG ~ infection_isolate*Genotype, data = MyListDF_64_88$full, dist = "negbin")

summary(mod1)
summary(mod2)

AIC(mod1)
AIC(mod2)

mod1$twologlik/2
mod2$loglik

# mod2 selected. 

# RES
lapply(MyListDF_64_88, function(xlist){testSignif(xlist, "RES")$LRT})# similar
lapply(MyListDF_64_88, function(xlist){testSignif(xlist, "RES_ZI")$LRT})

# IMP
lapply(MyListDF_64_88, function(xlist){testSignif(xlist, "IMP")$LRT})
## multiple comparison for impact on weight
MyListDF_64_88$full$MoPa <- as.factor(paste(MyListDF_64_88$full$Genotype, MyListDF_64_88$full$infection_isolate))
mimp <- lm(relWL~MoPa, data = MyListDF_64_88$full)
ANOVA <- aov(mimp)
# Tukey test to study each pair of treatment :
TUKEY <- TukeyHSD(x=ANOVA, 'MoPa', conf.level=0.95)
# Extract only significant differences:
tukeyDF <- data.frame(TUKEY$MoPa)
tukeyDF[tukeyDF$p.adj<0.05,]
# Keep only statistical differences WITHIN one parasite isolate or one mouse group

# TOL
lapply(MyListDF_64_88, function(xlist){testSignif(xlist, "TOL")$LRT})

##### 2.1. within E64
# RES
lapply(MyListDF_64, function(xlist){testSignifWithinParas(xlist, "RES")$LRT})
# IMP
lapply(MyListDF_64, function(xlist){testSignifWithinParas(xlist, "IMP")$LRT})
# TOL
lapply(MyListDF_64, function(xlist){testSignifWithinParas(xlist, "TOL")$LRT})

##### 2.2. within E88
# RES
lapply(MyListDF_88, function(xlist){testSignifWithinParas(xlist, "RES_ZI")$LRT})
# IMP
lapply(MyListDF_88, function(xlist){testSignifWithinParas(xlist, "IMP")$LRT})
# TOL
lapply(MyListDF_88, function(xlist){testSignifWithinParas(xlist, "TOL")$LRT})

##### 2.3. plot within E64 (figure 4)
## Plot Figure 4:
listPlotRes_64 <- lapply(MyListDF_64, function(x){get_plot_par(x, npara = 1,  plottype = "RES", cols = c("blue", "blue", "blue", "purple", "purple", "red", "red","red", "red"),
                                                               pointshape = c(15,16,17,15,16,15,16, 17))})
listPlotImp_64 <-lapply(MyListDF_64, function(x){get_plot_par(x, npara = 1,  plottype = "IMP", cols = c("blue", "blue", "blue", "purple", "purple", "red", "red","red", "red"),
                                                              pointshape = c(15,16,17,15,16,15,16, 17))})
listPlotTol_64 <- lapply(MyListDF_64, function(x){get_plot_par(x, npara = 1,  plottype = "TOL", cols = c("blue", "blue", "blue", "purple", "purple", "red", "red","red", "red"),
                                                               pointshape = c(15,16,17,15,16,15,16, 17),
                                                               linetype = c(1,2,3,1,2,1,2,3),)})
listBigPlot_64 <- lapply(MyListDF_64, function(x){
    getCouplingPlot(x, posx1 = 1.5e6, posy1 = 0.01, posx2 = 1.5e6, posy2 = -0.09, mycolors = c("blue", "blue", "blue", "purple", "purple", "red", "red","red", "red"),
                    myshapes = c(15,16,17,15,16,15,16, 17))
})
# to correct
make5panelsPlot <- function(res, imp, tol, bp1, bp2){
    p1 <- plot_grid(res + theme(legend.position = "none"), labels = c('A'), label_size = 20)
    p2 <- plot_grid(p1, tol, labels = c("", "C"), label_size = 20,
                    ncol =1, align = "v", rel_widths = c(1.5, .5))
    p3 <- plot_grid(bp1, bp2, ncol =1, align = "v", labels = c("D", "E"), label_size = 20)
    return(plot_grid(p2,p3,rel_widths = c(1, 1)))
} 
Fig4 <- make5panelsPlot(res = listPlotRes_64$full, imp = listPlotImp_64$full, tol = listPlotTol_64$full, 
                        bp1 = listBigPlot_64$full[[1]], bp2 = listBigPlot_64$full[[2]])

pdf(file ="../figures/Fig4_temp.pdf",  width = 15, height = 10)
Fig4
dev.off()

# Suppl
SupplFig4 <- make5panelsPlot(res = listPlotRes_64$cons, imp = listPlotImp_64$cons, tol = listPlotTol_64$cons, 
                             bp1 = listBigPlot_64$cons[[1]], bp2 = listBigPlot_64$cons[[2]])

pdf(file ="../figures/SupplFig4_temp.pdf",  width = 15, height = 10)
SupplFig4
dev.off()

##### 2.4. plot within E88 (figure 6)
listPlotRes_88 <- lapply(MyListDF_88, function(x){get_plot_par(x, npara = 1,  plottype = "RES_ZI", cols = c("blue", "blue", "blue", "purple", "purple", "red", "red","red", "red"),
                                                               pointshape = c(15,16,17,15,16,15,16, 17))})
listPlotImp_88 <-lapply(MyListDF_88, function(x){get_plot_par(x, npara = 1,  plottype = "IMP", cols = c("blue", "blue", "blue", "purple", "purple", "red", "red","red", "red"),
                                                              pointshape = c(15,16,17,15,16,15,16, 17))})
listPlotTol_88 <- lapply(MyListDF_88, function(x){get_plot_par(x, npara = 1,  plottype = "TOL", cols = c("blue", "blue", "blue", "purple", "purple", "red", "red","red", "red"),
                                                               pointshape = c(15,16,17,15,16,15,16, 17),
                                                               linetype = c(1,2,3,1,2,1,2,3),)})
listBigPlot_88 <- lapply(MyListDF_88, function(x){
    getCouplingPlot(x,  isZINB=TRUE, posx1 = 1.5e6, posy1 = 0.01, posx2 = 1.5e6, posy2 = -0.09, mycolors = c("blue", "blue", "blue", "purple", "purple", "red", "red","red", "red"),
                    myshapes = c(15,16,17,15,16,15,16, 17))
})

# test outlier
MyListDF_88_full_noPWD <- MyListDF_88$full[!MyListDF_88$full$Mouse_genotype %in% "PWD",]
listBigPlot_88_noPWD <- getCouplingPlot(
    MyListDF_88_full_noPWD, isZINB=TRUE, posx1 = 1.5e6, posy1 = 0.01, posx2 = 1.5e6, posy2 = -0.09, 
    mycolors = c("blue", "blue", "blue", "purple", "purple", "red", "red","red", "red"),
    myshapes = c(15,16,17,15,16,15,16, 17))
listBigPlot_88_noPWD
# test outlier

Fig5 <- make5panelsPlot(res = listPlotRes_88$full, imp = listPlotImp_88$full, tol = listPlotTol_88$full, 
                        bp1 = listBigPlot_88$full[[1]], bp2 = listBigPlot_88$full[[2]])

Fig5

pdf(file ="../figures/Fig5_temp.pdf",  width = 15, height = 10)
Fig5
dev.off()

# Suppl
SupplFig5 <- make5panelsPlot(res = listPlotRes_88$cons, imp = listPlotImp_88$cons, tol = listPlotTol_88$cons, 
                             bp1 = listBigPlot_88$cons[[1]], bp2 = listBigPlot_88$cons[[2]])

pdf(file ="../figures/SupplFig5_temp.pdf",  width = 15, height = 10)
SupplFig5
dev.off()

## NEW. test coupling TOGETHER of E64 and E88 for resistance and tolerance
# Get predicted values for resistance, impact or tolerance:
dfResPred64_88 <- getPredictedValues(MyListDF_64_88$full, "RES", 2)
dfImpPred64_88 <- getPredictedValues(MyListDF_64_88$full, "IMP", 2)
dfTolPred64_88 <- getPredictedValues(MyListDF_64_88$full, "TOL", 2)
# merge
names(dfResPred64_88)[names(dfResPred64_88) %in% c("predicted", "conf.low",  "std.error", "conf.high")] <-
    paste(names(dfResPred64_88)[names(dfResPred64_88) %in% c("predicted", "conf.low",  "std.error", "conf.high")], "OPG", sep = "_")
names(dfImpPred64_88)[names(dfImpPred64_88) %in% c("predicted", "conf.low",  "std.error", "conf.high")] <-
    paste(names(dfImpPred64_88)[names(dfImpPred64_88) %in% c("predicted", "conf.low",  "std.error", "conf.high")], "relWL", sep = "_")
names(dfTolPred64_88)[names(dfTolPred64_88) %in% c("predicted", "conf.low",  "std.error", "conf.high")] <-
    paste(names(dfTolPred64_88)[names(dfTolPred64_88) %in% c("predicted", "conf.low",  "std.error", "conf.high")], "Tol", sep = "_")

finalPredDF_64_88 <- merge(merge(dfResPred64_88, dfImpPred64_88), dfTolPred64_88)

# OPG-relWL
lm1 <- lm(predicted_OPG ~ predicted_relWL * infection_isolate, data = finalPredDF_64_88)
anova(lm1) # interaction significant
# Obtain slopes
lm1$coefficients # positive for Efer, ne
# plot 
ggplot(finalPredDF_64_88, aes(x= predicted_OPG, y = predicted_relWL))+
    geom_point(aes(col = Genotype), size = 4) +
    geom_smooth(method = "lm", aes(linetype=infection_isolate))

# remove outlier PWD
lm1 <- lm(predicted_OPG ~ predicted_relWL * infection_isolate, data = finalPredDF_64_88[finalPredDF_64_88$Genotype!= "8. PWD",])
anova(lm1) # interaction significant
# Obtain slopes
lm1$coefficients
lm1.lst <- lstrends(lm1, "infection_isolate", var="predicted_relWL")
# Compare slopes
pairs(lm1.lst) # difference of slope
# plot 
ggplot(finalPredDF_64_88[finalPredDF_64_88$Genotype!= "8. PWD",], aes(x= predicted_OPG, y = predicted_relWL))+
    geom_point(aes(col = Genotype)) +
    geom_smooth(method = "lm", aes(col=infection_isolate))

# RES-TOL
lm2 <- lm(predicted_OPG ~ predicted_Tol * infection_isolate, data = finalPredDF_64_88)
lm2
anova(lm2) # NO significant interactions, but tolerance different for both parasites
# Obtain slopes
lm2$coefficients
lm2.lst <- lstrends(lm2, "infection_isolate", var="predicted_Tol")
# Compare slopes
pairs(lm2.lst) # no difference of slope
# plot 
ggplot(finalPredDF_64_88, aes(x= predicted_OPG, y = -predicted_Tol))+
    geom_point(aes(col = Genotype)) +
    geom_smooth(method = "lm", aes(linetype=infection_isolate))
    
# remove outlier PWD
lm2 <- lm(predicted_OPG ~ -predicted_Tol * infection_isolate, data = finalPredDF_64_88[finalPredDF_64_88$Genotype!= "8. PWD",])
anova(lm2) # NO significant interactions, but tolerance different for both parasites
# Obtain slopes
lm2$coefficients
lm2.lst <- lstrends(lm2, "infection_isolate", var="predicted_Tol")
# Compare slopes
pairs(lm2.lst) # no difference of slope
# plot 
ggplot(finalPredDF_64_88[finalPredDF_64_88$Genotype!= "8. PWD",], aes(x= predicted_OPG, y = -predicted_Tol))+
    geom_point(aes(col = Genotype)) +
    geom_smooth(method = "lm", aes(col=infection_isolate))

############### Local adaptation of pure strains for E. ferrisi
MyListDF_locad <- MyListDF

MyListDF_locad$full$infection_isolate <- relevel(MyListDF_locad$full$infection_isolate, "Brandenburg64 (E. ferrisi)")
MyListDF_locad$full <- MyListDF_locad$full[!grepl("-", MyListDF_locad$full$Mouse_genotype),]
MyListDF_locad$full <- dropLevelsAllFactorsDF(
    MyListDF_locad$full[grep("ferrisi", MyListDF_locad$full$infection_isolate),])

MyListDF_locad$cons$infection_isolate <- relevel(MyListDF_locad$cons$infection_isolate, "Brandenburg64 (E. ferrisi)")
MyListDF_locad$cons <- MyListDF_locad$cons[!grepl("-", MyListDF_locad$cons$Mouse_genotype),]
MyListDF_locad$cons <- dropLevelsAllFactorsDF(
    MyListDF_locad$cons[grep("ferrisi", MyListDF_locad$cons$infection_isolate),])

## Matching vs unmatching
MyListDF_locad$full$combi <- "unmatching"
MyListDF_locad$full$combi[(MyListDF_locad$full$Mouse_subspecies %in% "M.m.dom" & MyListDF_locad$full$infection_isolate %in% "Brandenburg64 (E. ferrisi)")|
                        (MyListDF_locad$full$Mouse_subspecies %in% "M.m.mus" & MyListDF_locad$full$infection_isolate %in% "Brandenburg139 (E. ferrisi)")] <- "matching"
MyListDF_locad$full$combi <- factor(MyListDF_locad$full$combi)
MyListDF_locad$cons$combi <- "unmatching"
MyListDF_locad$cons$combi[(MyListDF_locad$cons$Mouse_subspecies %in% "M.m.dom" & MyListDF_locad$cons$infection_isolate %in% "Brandenburg64 (E. ferrisi)")|
                              (MyListDF_locad$cons$Mouse_subspecies %in% "M.m.mus" & MyListDF_locad$cons$infection_isolate %in% "Brandenburg139 (E. ferrisi)")] <- "matching"
MyListDF_locad$cons$combi <- factor(MyListDF_locad$cons$combi)

### Res
# no zero, so we use "classic" negbin:
modResLA <- glm.nb(max.OPG ~ combi, data = MyListDF_locad$full)
summary(modResLA) # p= 0.00323 ** 

plot_model(modResLA, type = "pred")#, show.data = TRUE)

ggplot(MyListDF_locad$full, aes(y = max.OPG, x = combi)) +
    geom_boxplot() +
    geom_jitter(aes(col= Genotype, shape = infection_isolate), size = 4) +
    scale_shape_manual(values = c(21,16)) +
    scale_color_manual(values = c("blue", "blue4", "red", "red4")) +
    facet_grid(.~Exp_ID)

#cons
summary(glm.nb(max.OPG ~ combi, data = MyListDF_locad$cons))
plot_model(glm.nb(max.OPG ~ combi, data = MyListDF_locad$cons), type = "pred")
ggplot(MyListDF_locad$cons, aes(y = max.OPG, x = combi)) +
    geom_boxplot() +
    geom_jitter(aes(col= Genotype, shape = infection_isolate), size = 4) +
    scale_shape_manual(values = c(21,16)) +
    scale_color_manual(values = c("blue", "blue4", "red", "red4")) +
    facet_grid(.~Exp_ID)

### Tol
modTolLA <- lm(relWL ~ 0 + max.OPG : combi, data = MyListDF_locad$full)
modTolLA
summary(modTolLA)
anova(modTolLA) # p=2.664e-15 ***

plot_model(modTolLA)
plot_model(modTolLA, type = "int", terms = "combi", show.data = TRUE)

ggplot(MyListDF_locad$full, aes(x=max.OPG, y=relWL, col = combi)) +
    geom_smooth(method = "lm", formula=y~x-1) +
    geom_point()

##LOCAL ADAPTATION OLD

MyListDF_locad$full$infection_isolate <- relevel(MyListDF_locad$full$infection_isolate, "Brandenburg64 (E. ferrisi)")
MyListDF_locad$full <- MyListDF_locad$full[!grepl("-", MyListDF_locad$full$Mouse_genotype),]
MyListDF_locad$full <- dropLevelsAllFactorsDF(
        MyListDF_locad$full[grep("ferrisi", MyListDF_locad$full$infection_isolate),])

MyListDF_locad$cons$infection_isolate <- relevel(MyListDF_locad$cons$infection_isolate, "Brandenburg64 (E. ferrisi)")
MyListDF_locad$cons <- MyListDF_locad$cons[!grepl("-", MyListDF_locad$cons$Mouse_genotype),]
MyListDF_locad$cons <- dropLevelsAllFactorsDF(
        MyListDF_locad$cons[grep("ferrisi", MyListDF_locad$cons$infection_isolate),])

### Res
lapply(MyListDF_locad, function(x){testSignif(x,"RES")$LRT}) # interaction factor significant
lapply(MyListDF_locad, function(x){testSignifSUBSPECIES(x,"RES")$LRT}) # interaction factor not significant
mR <- testSignifSUBSPECIES(MyListDF_locad$full,"RES")
plot_model(mR$)
### Imp
lapply(MyListDF_locad, function(x){testSignif(x,"IMP")$LRT}) # interaction factor not significant (apart conserv.)
### Tol
lapply(MyListDF_locad, function(x){testSignif(x,"TOL")$LRT}) # interaction factor not significant
lapply(MyListDF_locad, function(x){testSignifSUBSPECIES(x,"TOL")$LRT})

# ## Plot Figure 3:
# listPlotRes_LA <- lapply(MyListDF_locad, function(x){get_plot_par(x, npara = 2, plottype = "RES", cols = c("blue", "blue", "red", "red"), pointshape = c(15,16,15,16))})
# listPlotImp_LA <- lapply(MyListDF_locad, function(x){get_plot_par(x, npara = 2, plottype = "IMP", cols = c("blue", "blue", "red", "red"), pointshape = c(15,16,15,16))})
# listPlotTol_LA <- lapply(MyListDF_locad, function(x){get_plot_par(x, npara = 2, plottype = "TOL", linetype = c(1,2,1,2), cols = c("blue", "blue", "red", "red"), pointshape = c(15,16,15,16))})
# 
# Fig3 <- cowplot::ggdraw() +
#         draw_plot(listPlotRes_LA$full + theme(legend.position = "none"), 0, .5, .49, .49) +
#         draw_plot(listPlotImp_LA$full + theme(legend.position = "none"), .5, .5, .49, .49) +
#         draw_plot(listPlotTol_LA$full, 0, 0, .7, .49) +
#         draw_plot_label(c("A", "B", "C"), c(.01, .51, .01), c(1, 1, .5), size = 15)
# 
# pdf(file = "../figures/Fig3_temp.pdf",
#     width = 10, height = 10)
# Fig3
# dev.off()
# 
# pdf(file = "../figures/SupplS2_Fig3_temp.pdf",
#     width = 10, height = 10)
# cowplot::ggdraw() +
#         draw_plot(listPlotRes_LA$cons + theme(legend.position = "none"), 0, .5, .49, .49) +
#         draw_plot(listPlotImp_LA$cons + theme(legend.position = "none"), .5, .5, .49, .49) +
#         draw_plot(listPlotTol_LA$cons+ theme(legend.position = "none"), 0, 0, .7, .49) +
#         draw_plot_label(c("A", "B", "C"), c(.01, .51, .01), c(1, 1, .5), size = 15)
# dev.off()
