## 
## Initialize
##

library(ggplot2)
library(colorspace) # for darken()
library(cowplot)

ChromMidPointList = c(
  115109,
  636810,
  1201712,
  2125988.5,
  3180392,
  3603909.5,
  4284460,
  5111251.5,
  5612517,
  6205336.5,
  6911620,
  7784116.5,
  8785420.5,
  9639802.5,
  10577614.5,
  11597293)

ChromRunLengthList <- c(
  0,
  230218,
  1043402,
  1360022,
  2891955,
  3468829,
  3738990,
  4829930,
  5392573,
  5832461,
  6578212,
  7245028,
  8323205,
  9247636,
  10031969,
  11123260,
  12071326)

ChromRomanList <- c("I",
                    "II",
                    "III",
                    "IV",
                    "V",
                    "VI",
                    "VII",
                    "VIII",
                    "IX",
                    "X",
                    "XI",
                    "XII", 
                    "XIII",
                    "XIV", 
                    "XV", 
                    "XVI")


##
## Load Data
##

Dir="/Users/livinlrg/Desktop/dunham/Cris_L/PostDocBrewing/WorkDirectory/Clones"

#P1c1 <- read.table(paste(Dir,'/PDB1_c1/PDB1_c1.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)
#P1c2 <- read.table(paste(Dir,'/PDB1_c2/PDB1_c2.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)

ClonesList <- c(1,10,11,12,13,2,20,23,3,4,5,6,7,9)

Chrom=16

for(i in ClonesList) {
  #assign(paste("P48c", i, sep=""),
  #       read.table(paste(Dir,'/PDB48_c',i,'/PDB48_c',i,'.ADALL_SNPS.txt',sep = ""), sep='\t', head=FALSE))
  assign(paste("subP48c", i, sep=""),
         subset(get(paste("P48c",i,sep="")),V1==Chrom))
}

subP1c1 <- subset(P1c1,V1==Chrom)
subP1c2 <- subset(P1c2,V1==Chrom)

PlotFunc48 <- function(CloneValues, Number) {
  print(Number)
  return(ggplot() + 
    geom_point(data=get(paste("subP48c",Number,sep="")),
               aes(x=get(paste("subP48c",Number,sep=""))$V2,
                   y=get(paste("subP48c",Number,sep=""))$V3,
                   color="End", fill="End",
                   alpha=0.5)) +
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
    scale_color_manual(
      name = NULL,
      values = darken(c("Start"='#33b4ff',
                        "End"='#e69f00')), 0.3) +
    scale_fill_manual(
      name = NULL,
      values = c("Start"='#33b4ff',
                 "End"='#e69f00')) +
    scale_y_continuous(name = "Allele Frequency",
                       expand = c(0.025, 0),
                       breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) +
    scale_x_continuous(name = "Genome Position",
                       expand = c(0.0125, 0.0125),
                       breaks = ChromMidPointList[Chrom],
                       labels = c(paste("PDB26 Clone ", Number,sep=""))) +
    theme(line = element_line(color="grey80", size=1),
          axis.line = element_line(color="grey80", size=1),
          axis.line.x = element_line(color="grey80", size=1),
          axis.line.y = element_line(color="grey80", size=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.y = element_line(color = "grey80", size = 1),
          axis.ticks.x = element_blank(),
          text = element_text(color = "black", size=12),
          axis.text.x.bottom = element_text(color = "black", size=12),
          legend.position="none",
          panel.background = element_rect(fill = "white"))
    )
}

for(i in ClonesList) {
  assign(paste("plotP48c",i,sep = ""), PlotFunc48("Hmm", i))
}


plotP1c2 <- ggplot() + 
  geom_point(data=subP1c2,
             aes(x=subP1c2$V2,
                 y=subP1c2$V3,
                 color="End", fill="End",
                 alpha=0.5)) +
  #geom_vline(xintercept = 5411222 + 600000, linetype = "solid", color = "grey80") +
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
  scale_color_manual(
    name = NULL,
    values = darken(c("Start"='#33b4ff',
                      "End"='#e69f00')), 0.3) +
  scale_fill_manual(
    name = NULL,
    values = c("Start"='#33b4ff',
               "End"='#e69f00')) +
  scale_y_continuous(name = "Allele Frequency",
                     expand = c(0.025, 0),
                     breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) +
  scale_x_continuous(name = "Genome Position",
                     expand = c(0.0125, 0.0125),
                     breaks = ChromMidPointList[Chrom],
                     labels = "PDB1 Clone 2") +
  theme(line = element_line(color="grey80", size=1),
        axis.line = element_line(color="grey80", size=1),
        axis.line.x = element_line(color="grey80", size=1),
        axis.line.y = element_line(color="grey80", size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(color = "grey80", size = 1),
        axis.ticks.x = element_blank(),
        text = element_text(color = "black", size=12),
        axis.text.x.bottom = element_text(color = "black", size=12),
        legend.position="none",
        panel.background = element_rect(fill = "white"))

plotP1c1 <- ggplot() + 
  geom_point(data=subP1c1,
             aes(x=subP1c1$V2,
                 y=subP1c1$V3,
                 color="End", fill="End",
                 alpha=0.5)) +
  #geom_vline(xintercept = 5411222 + 600000, linetype = "solid", color = "grey80") +
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
  scale_color_manual(
    name = NULL,
    values = darken(c("Start"='#33b4ff',
                      "End"='#e69f00')), 0.3) +
  scale_fill_manual(
    name = NULL,
    values = c("Start"='#33b4ff',
               "End"='#e69f00')) +
  scale_y_continuous(name = "Allele Frequency",
                     expand = c(0.025, 0),
                     breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) +
  scale_x_continuous(name = "Genome Position",
                     expand = c(0.0125, 0.0125),
                     breaks = ChromMidPointList[Chrom],
                     labels = "PDB1 Clone 1") +
  theme(line = element_line(color="grey80", size=1),
        axis.line = element_line(color="grey80", size=1),
        axis.line.x = element_line(color="grey80", size=1),
        axis.line.y = element_line(color="grey80", size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(color = "grey80", size = 1),
        axis.ticks.x = element_blank(),
        text = element_text(color = "black", size=12),
        axis.text.x.bottom = element_text(color = "black", size=12),
        legend.position="none",
        panel.background = element_rect(fill = "white"))


pdf(paste('/Volumes/GoogleDrive/My\ Drive/PostdocBrewing_Manuscript/Figures/Clones_Chr', ChromRomanList[Chrom], '_AF_Plot_20200620.pdf', sep = ''),height=12, width=12, title="Postdoc_Brewing")
ggdraw() +
  draw_plot(plotP1c1, x = 0,    y = 0.75, height = 0.25, width = 0.25) +
  draw_plot(plotP1c2, x = 0.25, y = 0.75, height = 0.25, width = 0.25) +
  draw_plot(PlotFunc48(plotP48c1, 1), x = 0.5,  y = 0.75, height = 0.25, width = 0.25) +
  draw_plot(PlotFunc48(plotP48c2, 2), x = 0.75, y = 0.75, height = 0.25, width = 0.25) +
  draw_plot(PlotFunc48(plotP48c3, 3), x = 0,  y = 0.5, height = 0.25, width = 0.25) +
  draw_plot(PlotFunc48(plotP48c4, 4), x = 0.25, y = 0.5, height = 0.25, width = 0.25) +
  draw_plot(PlotFunc48(plotP48c5, 5), x = 0.5,    y = 0.5, height = 0.25, width = 0.25) +
  draw_plot(PlotFunc48(plotP48c6, 6), x = 0.75, y = 0.5, height = 0.25, width = 0.25) +
  draw_plot(PlotFunc48(plotP48c7, 7), x = 0,  y = 0.25, height = 0.25, width = 0.25) +
  draw_plot(PlotFunc48(plotP48c9, 9), x = 0.25, y = 0.25, height = 0.25, width = 0.25) +
  draw_plot(PlotFunc48(plotP48c10, 10), x = 0.5, y = 0.25, height = 0.25, width = 0.25) +
  draw_plot(PlotFunc48(plotP48c11, 11), x = 0.75,    y = 0.25, height = 0.25, width = 0.25) +
  draw_plot(PlotFunc48(plotP48c12, 12), x = 0, y = 0, height = 0.25, width = 0.25) +
  draw_plot(PlotFunc48(plotP48c13, 13), x = 0.25,  y = 0, height = 0.25, width = 0.25) +
  draw_plot(PlotFunc48(plotP48c20, 20), x = 0.5,    y = 0, height = 0.25, width = 0.25) +
  draw_plot(PlotFunc48(plotP48c23, 23), x = 0.75, y = 0, height = 0.25, width = 0.25)
dev.off()