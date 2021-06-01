##############
## Flocculation
## Chris Large
##############

## Load in necesary libraries
library(ggplot2)
library(colorspace) # for darken()
library(cowplot)
library(extrafont)
library(dplyr)
library(reshape2)
library(tidyr)
loadfonts()

FloccData <- read.table('/Users/livinlrg/Google\ Drive/Dunham\ Lab/Postdoc_Brewing/Flocculation/Aggregate_ForR.csv',
                        sep=',',
                        head=TRUE)

FloccData.melt <- melt(FloccData, na.rm = TRUE)

SelectedData = c("PDB1c1","PDB1c2","PDB26c1",
                 "PDB26c7","PDB26c12","PDB26c2","PDB26c4","PDB26c6","PDB26c23",
                 "PositiveControl")

FloccData.melt.select <- FloccData.melt[which(FloccData.melt$Strain %in% SelectedData),]

FloccPlot <- ggplot(FloccData.melt.select, aes(x = Strain, y = value)) +
  geom_boxplot(width = 0.4,
               outlier.size = 1,
               outlier.alpha = 0) +
  geom_dotplot(binaxis = 'y',
               stackratio = 1,
               stackdir='center',
               dotsize = .5,
               colour = "black",
               fill = "grey80",
               alpha = 0.75) +
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
  scale_y_continuous(name = "Percent Flocculant",
                     breaks = c(0.0, 20, 40, 60, 80, 100),
                     minor_breaks = seq(0,100,2.5),
                     limits = c(35,100)) +
  scale_x_discrete(limits = SelectedData,
                   labels = c("PDB1\nClone 1",
                              "PDB1\nClone 2",
                              "PDB26\nClone 1",
                              "PDB26\nClone 7",
                              "PDB26\nClone 12",
                              "PDB26\nClone 2",
                              "PDB26\nClone 4",
                              "PDB26\nClone 6",
                              "PDB26\nClone 23",
                              "Positive\nControl")) +
  theme(line = element_line(color="grey80", size=1),
        axis.line = element_line(color="grey80", size=1),
        axis.line.x.bottom = element_line(color="grey80", size=1),
        axis.line.y.left = element_line(color="grey80", size=1),
        axis.line.y.right = element_line(color="grey80", size=1),
        axis.line.x.top = element_line(color="grey80", size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey80", size=.25),
        panel.grid.minor.y = element_line(color="grey80", size=.25),
        axis.ticks.y = element_line(color = "grey80", size = 1),
        axis.ticks.x = element_line(color = "grey80", size = 1),
        text = element_text(family = "sans", color = "black", size=12),
        axis.text = element_text(family = "sans", color = "black", size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        panel.background = element_rect(fill = "white"))

kruskal.test(value ~ Strain,data = FloccData.melt.select)

# Test between mitotic recombination and not
mann (FloccData.melt.select[FloccData.melt.select$Strain == "PDB1c1" |
                                  FloccData.melt.select$Strain == "PDB1c2",]$value,
            FloccData.melt.select[FloccData.melt.select$Strain == "PDB26c2" |
                                  FloccData.melt.select$Strain == "PDB26c4" |
                                  FloccData.melt.select$Strain == "PDB26c6" |
                                  FloccData.melt.select$Strain == "PDB26c23",]$value)

# Test between BAT1 loss and not
wilcox.test(FloccData.melt.select[FloccData.melt.select$Strain == "PDB1c1" |
                                    FloccData.melt.select$Strain == "PDB1c2",]$value,
            FloccData.melt.select[FloccData.melt.select$Strain == "PDB26c7" |
                                    FloccData.melt.select$Strain == "PDB26c12",]$value)


FloccPlot <- ggplot(FloccData.melt.select, aes(x = Strain, y = value)) +
  geom_col(position = position_dodge()) +
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
  scale_y_continuous(name = "Percent Flocculant",
                     breaks = c(0.0, 20, 40, 60, 80, 100),
                     minor_breaks = seq(0,100,2.5),
                     limits = c(0,100)) +
  scale_x_discrete(limits = SelectedData,
                   labels = c("PDB1\nClone 1",
                              "PDB1\nClone 2",
                              "PDB26\nClone 1",
                              "PDB26\nClone 7",
                              "PDB26\nClone 12",
                              "PDB26\nClone 2",
                              "PDB26\nClone 4",
                              "PDB26\nClone 6",
                              "PDB26\nClone 23",
                              "Positive\nControl")) +
  theme(line = element_line(color="grey80", size=1),
        axis.line = element_line(color="grey80", size=1),
        axis.line.x.bottom = element_line(color="grey80", size=1),
        axis.line.y.left = element_line(color="grey80", size=1),
        axis.line.y.right = element_line(color="grey80", size=1),
        axis.line.x.top = element_line(color="grey80", size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey80", size=.25),
        panel.grid.minor.y = element_line(color="grey80", size=.25),
        axis.ticks.y = element_line(color = "grey80", size = 1),
        axis.ticks.x = element_line(color = "grey80", size = 1),
        text = element_text(family = "sans", color = "black", size=12),
        axis.text = element_text(family = "sans", color = "black", size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        panel.background = element_rect(fill = "white"))

FloccData.melt.select.summary <- FloccData.melt.select %>% # the names of the new data frame and the data frame to be summarised
  group_by(Strain) %>%   # the grouping variable
  summarise(mean = mean(value),  # calculates the mean of each group
            sd = sd(value), # calculates the standard deviation of each group
            n = n(),  # calculates the sample size per group
            SE = sd(value)/sqrt(n())) # calculates the standard error of each group

FloccPlot <- ggplot() +
  geom_bar(data = FloccData.melt.select.summary,
           aes(x=Strain, y=mean),
           stat="identity",
           fill = "grey80",
           color = "black",
           width = 0.75) +
  geom_errorbar(data = FloccData.melt.select.summary, aes(x=Strain, 
                                                          ymin = mean - sd,
                                                          ymax = mean + sd),
           stat="identity",
           width = 0.25) +
  geom_dotplot(data = FloccData.melt.select, aes(x = Strain, y = value),
               binaxis = 'y',
               stackratio = 1,
               stackdir='center',
               dotsize = .5,
               colour = "black",
               fill = "grey80",
               alpha = 0.75) +
  annotate(geom = 'segment',
           y = Inf, yend = Inf,
           x = -Inf, xend = Inf,
           color = "black",
           size = 2) +
  annotate(geom = 'segment',
           y = -Inf, yend = Inf,
           x = Inf, xend = Inf,
           color = "black",
           size = 2) +
  scale_y_continuous(name = "Percent Flocculant",
                     breaks = seq(0,100,10),
                     minor_breaks = seq(0,100,2.5),
                     limits = c(0,100),
                     expand = expansion(mult = c(0, .05))) +
  scale_x_discrete(limits = SelectedData,
                   labels = c("PDB1\nClone 1",
                              "PDB1\nClone 2",
                              "PDB26\nClone 1",
                              "PDB26\nClone 7",
                              "PDB26\nClone 12",
                              "PDB26\nClone 2",
                              "PDB26\nClone 4",
                              "PDB26\nClone 6",
                              "PDB26\nClone 23",
                              "Positive\nControl")) +
  theme(line = element_line(color="black", size=1),
        axis.line = element_line(color="black", size=1),
        axis.line.x.bottom = element_line(color="black", size=1),
        axis.line.y.left = element_line(color="black", size=1),
        axis.line.y.right = element_line(color="black", size=1),
        axis.line.x.top = element_line(color="black", size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks.y = element_line(color = "black", size = 1),
        axis.ticks.x = element_line(color = "black", size = 1),
        text = element_text(family = "sans", color = "black", size=12),
        axis.text = element_text(family = "sans", color = "black", size=12),
        axis.title.x = element_blank(),
        legend.position="none",
        panel.background = element_rect(fill = "white"))
