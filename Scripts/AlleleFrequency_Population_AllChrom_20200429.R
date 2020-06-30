

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

Dir="/Users/livinlrg/Desktop/dunham/Cris_L/PostDocBrewing/WorkDirectory/"

PDB1 <- read.table(paste(Dir, 'Populations/PDB1/PDB1.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)
PDB9 <- read.table(paste(Dir, 'Populations/PDB9/PDB9.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)
PDB24 <- read.table(paste(Dir, 'Populations/PDB24/PDB24.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)
PDB34 <- read.table(paste(Dir, 'Populations/PDB34/PDB34.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)
PDB48 <- read.table(paste(Dir, 'Populations/PDB48/PDB48.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)

PD68 <- read.table(paste(Dir, 'Populations/PD68/PD68.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)
PD140 <- read.table(paste(Dir, 'Populations/PD140/PD140.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)

E01 <- read.table(paste(Dir, 'Populations/E01_S17/E01_S17.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)
E03 <- read.table(paste(Dir, 'Populations/E03_S18/E03_S18.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)
E07 <- read.table(paste(Dir, 'Populations/E07_S20/E07_S20.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)
E10 <- read.table(paste(Dir, 'Populations/E10_S22/E10_S22.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)
E08 <- read.table(paste(Dir, 'Populations/E08_S22/E08_S22.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)
E09 <- read.table(paste(Dir, 'Populations/E09_S23/E09_S23.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)

WLP001 <- read.table(paste(Dir, 'Populations/WLP001/WLP001_Preiss.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)
Drakes <- read.table(paste(Dir, 'Populations/YMD3980_Pop/YMD3980_Pop.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)
RedCircle <- read.table(paste(Dir, 'Populations/YMD4251/YMD4251.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)

PlotFunc <- function(StartData, EndData, Chrom, NameOfPlot) {
  return(ggplot() + 
    geom_point(data = StartData,
                aes(x = StartData$V2,
                    y = StartData$V3,
                    color = "Start",
                    fill = "Start"),
                alpha = 0.5, size = 1) +
    geom_point(data = EndData,
                aes(x = EndData$V2,
                    y = EndData$V3,
                    color = "End",
                    fill = "End"),
                alpha = 0.5, size = 1) +
    #  geom_jitter(data = CloneData,
    #              aes(x = CloneData$V2,
    #                  y = CloneData$V3,
    #                  color = "Clone",
    #                  fill = "Clone"),
    #              alpha = 0.5) +  
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
        values = darken(c("End"='#e69f00',
                          "Start"='#33b4ff')), 0.3) +
      scale_fill_manual(
        name = NULL,
        values = c("End"='#e69f00',
                   "Start"='#33b4ff')) +
    scale_y_continuous(name = "Allele Frequency",
                       expand = c(0.025, 0),
                       breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) +
    scale_x_continuous(name = paste(NameOfPlot,"\nGenome Position", sep=""),
                       expand = c(0.025, 0.025),
                       limits = c(ChromRunLengthList[as.numeric(Chrom)],
                                  ChromRunLengthList[as.numeric(Chrom) + 1])
#                       ,
#                       breaks = ChromMidPointList[as.numeric(Chrom)],
#                       labels = paste("Chromosome: ",ChromRomanList[as.numeric(Chrom)],sep = "")
    ) +
    theme(line = element_line(color="grey80", size=1),
          axis.line = element_line(color="grey80", size=1),
          axis.line.x = element_line(color="grey80", size=1),
          axis.line.y = element_line(color="grey80", size=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.y = element_line(color = "grey80", size = 1),
          axis.ticks.x = element_blank(),
          axis.title.y.left = element_blank(),
          text = element_text(color = "black", size=12),
          axis.text = element_text(family = "sans", color = "black", size=8),
          legend.position="none",
          axis.text.x.bottom = element_blank(),
          panel.background = element_rect(fill = "white")))
}

for(i in seq(1:16)) {
  pdf(paste('/Volumes/GoogleDrive/My\ Drive/PostdocBrewing_Manuscript/Figures/Population_AF_Chr',i,".pdf", sep = ""),height=4, width=16, title="Post-Doc_Brewing")
  print(ggdraw() +
    draw_plot(PlotFunc(PDB1, PDB48, i, "PDB Rep. 1"), x = 0, y = 0, height = 1, width = 1/7) +
    draw_plot(PlotFunc(PD68, PD140, i, "PDB Rep. 2"), x = 1/7, y = 0, height = 1, width = 1/7) +
    draw_plot(PlotFunc(E01, E03, i, "BreweryX Rep. 1"), x = 2/7, y = 0, height = 1, width = 1/7) +
    draw_plot(PlotFunc(E07, E10, i, "BreweryX Rep. 2"), x = 3/7, y = 0, height = 1, width = 1/7) +
    draw_plot(PlotFunc(E01, E09, i, "BreweryX Rep. 3"), x = 4/7, y = 0, height = 1, width = 1/7) +
    draw_plot(PlotFunc(WLP001, Drakes, i, "Drake's"), x = 5/7, y = 0, height = 1, width = 1/7) +
    draw_plot(PlotFunc(WLP001, RedCircle, i, "Red Circle"), x = 6/7, y = 0, height = 1, width = 1/7))
  dev.off()
}

