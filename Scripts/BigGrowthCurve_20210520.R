##############
## GrowthCurve
## Chris Large
##############
## For measuring lag time and maximum growth rate on plate reader data using growthrates
##

## Load in necesary libraries
library(reshape2)
library(tidyr)
library(dplyr)
library(ggplot2)
library(growthrates)
library(ggpubr)

Dir="/Users/livinlrg/Google Drive/Dunham Lab/Postdoc_Brewing/"

GrowthData <- read.csv(
  paste(Dir,
        "GrowthCurve_20210513.csv", sep ="")
  ,sep=",", header=TRUE, check.names = FALSE,stringsAsFactors=FALSE)

GrowthDataMelt <- melt(GrowthData, id.vars = c("Time"))

ParentalStrains <- c("PDB1c1",
                     "PDB1c2",
                     "PDB26c4",
                     "PDB26c7",
                     "PDB26c12")

GrowthDataMelt$ParentalName <- NA

for(i in seq(length(GrowthDataMelt$variable))) {
  for(j in ParentalStrains) {
    if(startsWith(toString(GrowthDataMelt[i,]$variable), j)) {
      GrowthDataMelt[i,]$ParentalName <- j
    }
  }
}

GrowthDataMeltMean <- GrowthDataMelt %>%
  group_by(Time, ParentalName) %>%
  summarize(mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE))

GrowthDataMeltMean$ParentalName <- factor(GrowthDataMeltMean$ParentalName,
                                           levels = c("PDB1c1", "PDB1c2", "PDB26c4", "PDB26c7", "PDB26c12"))

RawGrowthRate <- ggplot(GrowthDataMeltMean) +
  geom_point(aes(x = Time/60, y = mean, color = ParentalName, fill = ParentalName), alpha = 0.5) +
  #geom_pointrange(aes(x = Time/60, y = mean, ymax = mean + sd, ymin = mean - sd, color = ParentalName), 
  #                width = 0.5,
  #                alpha = 0.5,
  #                position = ) +
  geom_line(aes(x = Time/60, y = mean, color = ParentalName)) +
  #geom_smooth(aes(x = Time/60, y = mean, color = ParentalName, fill = ParentalName)) +
  scale_y_continuous(trans='log10', name = expression(paste(Log[10], "(OD600)"))) +
  scale_x_continuous(name = "Time (Hours)") +
  annotate(geom = 'segment',
           y = Inf, yend = Inf,
           x = -Inf, xend = Inf,
           color = "grey80",
           size = 2) +
  annotate(geom = 'segment',
           y = 0, yend = Inf,
           x = Inf, xend = Inf,
           color = "grey80",
           size = 2) +
  labs(color = "Strains", fill = "Strains") +
  theme(line = element_line(color="grey80", size=1),
        axis.line = element_line(color="grey80", size=1),
        axis.line.x = element_line(color="grey80", size=1),
        axis.line.y = element_line(color="grey80", size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(color = "grey80", size = 1),
        axis.ticks.x = element_line(color = "grey80", size = 1),
        text = element_text(color = "black", size=12),
        axis.text = element_text(family = "sans", color = "black", size=8),
#        legend.position="none",
        legend.key = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white"))


ResultsTable <- data.frame()

for(i in levels(factor(GrowthDataMelt$variable))) {
  fit <- fit_easylinear(GrowthDataMelt[GrowthDataMelt$variable == i,]$Time/60,
                        GrowthDataMelt[GrowthDataMelt$variable == i,]$value, h=6)
  print(coef(fit)[3])
  print(plot(fit, log = "y"))

  ResultsTable <- rbind(ResultsTable, data.frame(i, coef(fit)[1],coef(fit)[2],coef(fit)[3], coef(fit)[4], rsquared(fit)))
}

StrainNames <- c("PDB1c1",
                "PDB1c1",
                "PDB1c1",
                "PDB1c2",
                "PDB1c2",
                "PDB1c2",
                "PDB26c4",
                "PDB26c4",
                "PDB26c4",
                "PDB26c7",
                "PDB26c7",
                "PDB26c7",
                "PDB26c12",
                "PDB26c12",
                "PDB26c12")

GrowthDataMelt.Names <- cbind(ResultsTable, StrainNames)

GrowthDataMelt.Names$StrainNames <- factor(GrowthDataMelt.Names$StrainNames,
                                           levels = c("PDB1c1", "PDB1c2", "PDB26c4", "PDB26c7", "PDB26c12"))

GrowthRate <- ggplot(GrowthDataMelt.Names, aes(x = StrainNames,
                                 y = coef.fit..3., color = StrainNames)) + 
  geom_boxplot(width = 0.4, outlier.size = 1) +
  geom_dotplot(data = GrowthDataMelt.Names,
               aes(x = StrainNames,
                   y = coef.fit..3.),
               binaxis = 'y',
               stackratio = 1,
               stackdir='center',
               dotsize = .5,
               colour = "black",
               fill = "grey80",
               alpha = 0.75) +
  scale_y_continuous(name = expression(paste("Growth Rate (OD600/hour)"))) +
  scale_x_discrete(labels = c("PDB1\nClone 1",
                              "PDB1\nClone 2",
                              "PDB26\nClone 4",
                              "PDB26\nClone 7",
                              "PDB26\nClone 12")) +
  annotate(geom = 'segment',
           y = Inf, yend = Inf,
           x = -Inf, xend = Inf,
           color = "grey80",
           size = 2) +
  annotate(geom = 'segment',
           y = -Inf, yend = Inf,
           x = Inf, xend = Inf,
           color = "grey80",
           size = 2) +
  theme(line = element_line(color="grey80", size=1),
        axis.line = element_line(color="grey80", size=1),
        axis.line.x = element_line(color="grey80", size=1),
        axis.line.y = element_line(color="grey80", size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.ticks.y = element_line(color = "grey80", size = 1),
        axis.ticks.x = element_line(color = "grey80", size = 1),
        text = element_text(color = "black", size=12),
        axis.text = element_text(family = "sans", color = "black", size=8),
        legend.position="none",
        panel.background = element_rect(fill = "white"))
  
kruskal.test(coef.fit..3. ~ StrainNames,data = GrowthDataMelt.Names)

pdf(paste('/Users/livinlrg/Google\ Drive/PostdocBrewing_Manuscript/Figures/GrowthRate_20210521',".pdf", sep = ""),height=4, width=8, title="Post-Doc_Brewing")
print(ggdraw() +
        draw_plot(RawGrowthRate, x = 0, y = 0, height = 1, width = 0.60) +
        draw_plot(GrowthRate, x = 0.60, y = 0, height = 1, width = 0.40) +
        draw_plot_label(label = c("A", "B"), size = 12,
                        x = c(0, 0.60), y = c(1, 1)))
dev.off()
      

