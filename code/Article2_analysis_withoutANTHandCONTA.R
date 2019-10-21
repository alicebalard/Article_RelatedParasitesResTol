### Code for data analysis of Article 2 WITHOUT anth or COINF ANIMALs
### August 2019
### Alice Balard

## INFO
# Mouse AA_0088, HI = 0.2
# Mouse AA_0064, HI = 0.08
# Mouse AA_0139, HI = 0.85

#### Load data and functions within ####
source("dataPreparation.R")

######################################
########## Read information ##########
######################################

# full data
rawDF108mice <- DF_all[grep("F0", DF_all$Mouse_genotype),]

# calculate weight retained per day
d <- rawDF108mice[rawDF108mice$dpi == 0., c("weight", "EH_ID")]
names(d)[1] <- "weightdpi0"
rawDF108mice <- merge(d, rawDF108mice)
rawDF108mice$weightRelativeToInfection <- rawDF108mice$weight / rawDF108mice$weightdpi0 * 100

# Rename levels
levels(rawDF108mice$infection_isolate) <- c("Brandenburg139 (E. ferrisi)",
                                            "Brandenburg64 (E. ferrisi)", "Brandenburg88 (E. falciformis)")
# summary data
summaryDF108mice <- ALL_summary_F0
summaryDF108mice$Mouse_genotype <- droplevels(factor(summaryDF108mice$Mouse_genotype))
summaryDF108mice$Mouse_subspecies <- droplevels(factor(summaryDF108mice$Mouse_subspecies))
levels(summaryDF108mice$infection_isolate) <- c("Brandenburg139 (E. ferrisi)",
                                                "Brandenburg64 (E. ferrisi)", "Brandenburg88 (E. falciformis)")
summaryDF108mice$Eimeria_species <- as.factor(summaryDF108mice$Eimeria_species)
summaryDF108mice$Eimeria_species <- relevel(summaryDF108mice$Eimeria_species, "E.ferrisi")

## Set infection group (group 1 had anthelminthics, did not kill worms, stopped after: test effect)
summaryDF108mice$batch <- summaryDF108mice$Exp_ID
levels(summaryDF108mice$batch) <- c("B1", "B2", "B3", "B4", "B5", "B6")

summaryDF108mice$anth<- FALSE
summaryDF108mice$anth[summaryDF108mice$batch %in% "B1"] <- TRUE

# REMOVE animals with anth (N = 108 - 86 = 22)
sumDFnoANTH <- summaryDF108mice[summaryDF108mice$anth == FALSE,]


## Coinfection? (N = 86 -77 = 9)
contaAnimals <- rawDF108mice[rawDF108mice$oocysts.per.tube > 0 & !is.na(rawDF108mice$oocysts.per.tube) &
                rawDF108mice$dpi == 0, "EH_ID"]

sumDFnoANTHnoCONTA <- sumDFnoANTH[!sumDFnoANTH$EH_ID %in% contaAnimals,]

###############
## Resistance ##
################
# Mmd vs Mmm
modResSubsp <- glm.nb(peak.oocysts.per.g.mouse ~ Eimeria_species*Mouse_subspecies, data = sumDFnoANTHnoCONTA)
anova(modResSubsp)
length(na.omit(sumDFnoANTHnoCONTA$peak.oocysts.per.g.mouse))

# plot marginal effects of interaction terms
posx.1 <- c(0.9,1.1, 1.9,2.1)

## FOR PLOT, use Resistance Index 
y = seq(0,1,0.05)
x = -y * 300000 + 300000
plot_model(modResSubsp, type = "int", dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
  scale_color_manual(values = c("blue", "red"),
                     name = "Mouse subspecies",labels = c("Mmd", "Mmm")) +
  ggtitle("Resistance") +
  scale_y_continuous("(predicted) Resistance Index",
                     trans = "reverse",
                     breaks = x,
                     labels = y) +
  xlab("Eimeria species") +
  theme(axis.title.x = element_text(hjust=1), axis.text=element_text(size=13)) +
  geom_text(aes(x=posx.1,y=100000,label=getNs("peak.oocysts.per.g.mouse", sumDFnoANTHnoCONTA, 
                                             "Mouse_subspecies", "Eimeria_species")), vjust=0) 

############
## Impact ##
############
sumDFnoANTHnoCONTA$impact <- sumDFnoANTHnoCONTA$relWL +0.01
length(na.omit(sumDFnoANTHnoCONTA$impact))

modImpSubsp <- survreg(Surv(impact)~Eimeria_species*Mouse_subspecies, data = sumDFnoANTHnoCONTA, dist="weibull")
anova(modImpSubsp) # Eimeria species AND mouse subspecies significant
plot_model(modImpSubsp, type = "int",dot.size = 4, dodge = .5) 
plot_model(modImpSubsp, type = "int",dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
  scale_color_manual(values = c("blue","red"),
                     name = "Mouse subspecies",labels = c("Mmd", "Mmm")) +
  xlab("Eimeria species") +
  ggtitle("Impact on host health") +
  ylab("(predicted) maximum weight loss") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme(axis.title.x = element_text(hjust=1), axis.text=element_text(size=13)) +
  geom_text(aes(x=posx.1,y=0,label=getNs("relWL", sumDFnoANTHnoCONTA,
                                         "Mouse_subspecies", "Eimeria_species")),vjust=0)

###############
## Tolerance ##
###############
sumDFnoANTHnoCONTA$ToleranceIndex <- log10(
  sumDFnoANTHnoCONTA$relWL / sumDFnoANTHnoCONTA$peak.oocysts.per.g.mouse + 1e-8) / (-8)

modTolSubspecies <- lm(ToleranceIndex ~ Eimeria_species*Mouse_subspecies, data = sumDFnoANTHnoCONTA)
anova(modTolSubspecies)
length(na.omit(sumDFnoANTHnoCONTA$ToleranceIndex))


plot_model(modTolSubspecies, type = "int", dot.size = 4, dodge = .5)
plot_model(modTolSubspecies, type = "int", dot.size = 4, dodge = .5) + # mean-value and +/- 1 standard deviation
  scale_color_manual(values = c("blue", "red"),
                     name = "Mouse subspecies",labels = c("Mmd", "Mmm")) +
  xlab("Eimeria species") +
  ylab("(predicted) Tolerance index")+
  ggtitle("Tolerance") +
  theme(axis.title.x = element_text(hjust=1), axis.text = element_text(size=13))+
  geom_text(aes(x=posx.1,y=0.4,label=getNs("ToleranceIndex", sumDFnoANTHnoCONTA,
                                           "Mouse_subspecies", "Eimeria_species")),vjust=0)


# mouse genotype & interactions significant (& EIMERIA)

### Second part: correlation resistance / tolerance
getResistanceIndex <- function(x){
  y = (- x + 300000)/300000
  return(y)}

sumDFnoANTHnoCONTA$ResistanceIndex <- getResistanceIndex(sumDFnoANTHnoCONTA$peak.oocysts.per.g.mouse)

# Calculate mean per group
# Packages we need
library(ggplot2)
library(dplyr)

# Create a group-means data set
gd <- sumDFnoANTHnoCONTA %>% 
  group_by(Mouse_genotype, infection_isolate) %>% 
  summarise(
    ResistanceIndex = mean(ResistanceIndex, na.rm = T),
    ToleranceIndex = mean(ToleranceIndex, na.rm = T)
  )

gd

# define the 8 groups
sumDFnoANTHnoCONTA$group <- paste(sumDFnoANTHnoCONTA$Mouse_genotype, sumDFnoANTHnoCONTA$infection_isolate, sep = "_")

# Plot both data sets
restolplot <- ggplot(sumDFnoANTHnoCONTA, aes(x = ResistanceIndex, y = ToleranceIndex)) +
  geom_smooth(method = "lm", col = "black", alpha = .2, aes(linetype = Eimeria_species)) +
  geom_point(alpha = .4, aes(col = Mouse_genotype, fill = Mouse_genotype, shape = infection_isolate), size = 4) +
  geom_point(data = gd, aes(fill = Mouse_genotype, shape = infection_isolate), size = 10) +
  theme_bw()+
  scale_color_manual(values = c("blue", "cornflowerblue", "red4", "indianred1")) +
  scale_fill_manual(values = c("blue", "cornflowerblue", "red4", "indianred1")) +
  scale_shape_manual(values = c(24,22,21)) +
  ylab(label = "Tolerance index") +
  scale_x_continuous(name = "Resistance index")
restolplot

# https://stats.idre.ucla.edu/r/seminars/interactions-r/
# Test the difference between slopes:
modResTol <- lm(formula = ToleranceIndex ~ ResistanceIndex * Eimeria_species, data = sumDFnoANTHnoCONTA)

summary(modResTol)
#p-value of the t-statistic for the interaction between ResistanceIndex and Eimeria_species: 1.39e-05 ***

# ResistanceIndex:Eimeria_speciesE.falciformis -0.98029    0.27241  -3.599 0.000624 ***
# The interaction ResistanceIndex * Eimeria_species is significant, which suggests
# that the relationship of ToleranceIndex by ResistanceIndex does vary by Eimeria spieces

# Since our goal is to obtain simple slopes of each parasite:
library(lsmeans)
emtrends(modResTol, ~ Eimeria_species, var="ResistanceIndex")
# E.ferrisi                     -0.0806 0.0834 64   -0.247   0.0861
# E.falciformis                 -1.0609 0.2593 64   -1.579  -0.5428
