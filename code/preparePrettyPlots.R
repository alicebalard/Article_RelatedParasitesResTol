## PLOTS pretty for publication

source("Article2_analysis.R")

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

pdf(file = "../figures/mapIsolatesandStrains.pdf", width = 5, height =5)
map
dev.off()

  ## To add Ns on top of bars
getNs <- function(proxy, df){
  noNA = df[!is.na(df[[proxy]]),]
  tab = table(noNA$infection_isolate, noNA$Mouse_genotype)
  Ns = as.character(as.vector(t(tab)[as.vector(t(tab))!=0]))
  return(Ns)
}

## Fig 2

# plot marginal effects of interaction terms
posx <- c(0.8+c(0,1/8,2/8,3/8),1.8+c(0,1/8,2/8,3/8),2.8+c(0,1/8,2/8,3/8))

## and plot
plotR_F0 <- plot_model(modRes, type = "int", dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
  scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                     name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +
  ggtitle("Resistance", subtitle = "low parasite burden = high resistance") +
  scale_y_log10("(predicted) max parasite density", 
                breaks = c(2000,5000, 10000, 20000, 50000, 100000,200000), 
                labels = c(2,5, 10, 20, 50, 100,200))+
  scale_x_continuous("Eimeria isolate", breaks = 1:3, labels = c("Brandenburg139", "Brandenburg64", "Brandenburg88")) +
  theme(axis.title.x = element_text(hjust=1), axis.text=element_text(size=13)) +
  geom_text(aes(x=posx,y=3000,label=getNs("peak.oocysts.per.g.mouse", ALL_summary_F0)),vjust=0) 
plotR_F0

plotI_F0 <- plot_model(modImp, type = "int",dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
  scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                     name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +
  xlab("Eimeria isolate") +
  ggtitle("Impact on host health") +
  ylab("(predicted) max impact") +
  scale_x_continuous(breaks = 1:3, labels = c("Brandenburg139", "Brandenburg64", "Brandenburg88")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme(axis.title.x = element_text(hjust=1), axis.text=element_text(size=13)) +
  geom_text(aes(x=posx,y=0,label=getNs("relWL", ALL_summary_F0)),vjust=0)
plotI_F0

plotT_F0 <- plot_model(modTol, type = "int", dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
  scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1"),
                     name = "Mouse strain",labels = c("SCHUNT", "STRA", "BUSNA", "PWD")) +
  xlab("Eimeria isolate") +
  ylab("(predicted) tolerance index")+
  ggtitle("Tolerance", subtitle = "low values = high tolerance") +
  scale_x_continuous(breaks = 1:3, labels = c("Brandenburg139", "Brandenburg64", "Brandenburg88"))+
  theme(axis.title.x = element_text(hjust=1), axis.text = element_text(size=13))+
  geom_text(aes(x=posx,y=0.2,label=getNs("tolerance", ALL_summary_F0)),vjust=0)
plotT_F0

# Fig 2.
library(cowplot)
Fig2 <- plot_grid(plotR_F0 + theme(legend.position = "none"),
                  plotI_F0 + theme(legend.position = "none"),
                  plotT_F0 + theme(legend.position = "none"), 
                  labels=c("A", "B", "C"), label_size = 20)

Fig2

pdf(file = "../figures/Fig2.pdf",
    width = 9, height = 9)
Fig2
dev.off()

## Fig 3
## To add Ns on top of bars
getNsF1 <- function(proxy, df){
  noNA = df[!is.na(df[[proxy]]),]
  tab = table(noNA$hybridLevel, noNA$Eimeria_species)
  Ns = as.character(as.vector(t(tab)[as.vector(t(tab))!=0]))
  return(Ns)
}

# plot marginal effects of interaction terms
posxF1 <- c(0.8+c(0,3/8),1.8+c(0,3/8),2.8+c(0,3/8))

lmResF1plot <- plot_model(modResF1, type = "int", dot.size = 4, dodge = .5,
                          colors = c("orange", "darkgreen"))+
  theme_bw() +
  ggtitle(label = "Resistance", subtitle = "low parasite burden = high resistance")+
  scale_y_log10("(predicted) max parasite density", 
                breaks = c(20,50, 100, 200, 500, 1000,2000)*100, 
                labels = c(2,5, 10, 20, 50, 100,200))+
  xlab("Mouse genotype") +
  theme(axis.title.x = element_text(hjust=1), axis.text = element_text(size=13)) +
  geom_text(aes(x=c(0.9,1.1,1.9,2.1,2.9,3.1),
                y=13000,label=getNsF1("peak.oocysts.per.g.mouse", ALL_summary_F1)),vjust=0)

lmResF1plot

bananaResH3 <- parasiteLoad::bananaPlot(mod = fitRes$H3,
                                        data = ALL_summary_F1,
                                        response = "peak.oocysts.per.g.mouse.div100",
                                        islog10 = T, group = "Eimeria_species",
                                        cols = c("orange", "darkgreen")) + theme_bw() +
  ylab("Maximum parasite load (oocysts) per mouse gram") +
  xlab("Hybrid index") +
  ggtitle(label = "Resistance", subtitle = "test for hybrid vigor")+
  scale_y_log10("max parasite density", 
                breaks = c(20, 50, 100, 200, 500, 1000, 2000), 
                labels = c(2,5, 10, 20, 50, 100,200))+
  theme(axis.title.x = element_text(hjust=1), axis.text = element_text(size=13))

bananaResH3 

lmImpF1plot <- plot_model(modImpF1, type = "int", dot.size = 4, dodge = .5, 
                          colors = c("orange", "darkgreen"), title = "Impact on host health") + theme_bw() +
  ylab("(predicted) max impact") +
  xlab("Mouse genotype") +
  theme(axis.title.x = element_text(hjust=1), axis.text = element_text(size=13)) +
  geom_text(aes(x=c(0.9,1.1,1.9,2.1,2.9,3.1),
                y=0,label=getNsF1("relWL", ALL_summary_F1)),vjust=0) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))
lmImpF1plot

bananaImpH3 <- parasiteLoad::bananaPlot(mod = fitImp$H3,
                                        data = ALL_summary_F1,
                                        response = "impact",
                                        islog10 = F, group = "Eimeria_species",
                                        cols = c("orange", "darkgreen")) + theme_bw() +
  ylab("max impact") +
  xlab("Hybrid index") +
  ggtitle(label = "Impact of host health", subtitle = "test for hybrid vigor")+
  theme(axis.title.x = element_text(hjust=1), axis.text = element_text(size=13))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))
bananaImpH3 

lmTolF1plot <- plot_model(modTolF1, type = "int", dot.size = 4, dodge = .5, 
                          colors = c("orange", "darkgreen"), title = "Tolerance") + theme_bw() +
  ggtitle("Tolerance", subtitle = "low values = high tolerance") +
  xlab("Mouse genotype") +
  ylab("tolerance index") +
  theme(axis.title.x = element_text(hjust=1), axis.text = element_text(size=13))+
  geom_text(aes(x=c(0.9,1.1,1.9,2.1,2.9,3.1),
                y=0,label=getNsF1("tolerance", ALL_summary_F1)),vjust=0)
lmTolF1plot

bananaTolH3 <- parasiteLoad::bananaPlot(mod = fitTol$H3,
                                        data = ALL_summary_F1,
                                        response = "tolerance",
                                        islog10 = F, group = "Eimeria_species",
                                        cols = c("orange", "darkgreen")) + theme_bw() +
  xlab("Hybrid index") +
  ylab("tolerance index") +
  ggtitle("Tolerance", subtitle = "test for hybrid vigor") +
  theme(axis.title.x = element_text(hjust=1), axis.text = element_text(size=13))

bananaTolH3 

####### Full figure 3 ####### 

Fig3 <- plot_grid(lmResF1plot + theme(legend.position = "none"),
                  bananaResH3 + theme(legend.position = "none"),
                  lmImpF1plot + theme(legend.position = "none"),
                  bananaImpH3 + theme(legend.position = "none"),
                  lmTolF1plot + theme(legend.position = "none"),
                  bananaTolH3 + theme(legend.position = "none"),
                  labels=c("A", "B", "C", "D", "E", "F"), label_size = 15,
                  ncol = 2)

Fig3

pdf(file = "../figures/Fig3.pdf",
    width = 10, height = 12)
Fig3
dev.off()

