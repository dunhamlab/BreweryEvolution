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

SubSetOn = '8'

PDB1 <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Populations/PDB1/PDB1.ADALL_SNPS.txt', sep='\t', head=FALSE)
PDB9 <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Populations/PDB9/PDB9.ADALL_SNPS.txt', sep='\t', head=FALSE)
PDB24 <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Populations/PDB24/PDB24.ADALL_SNPS.txt', sep='\t', head=FALSE)
PDB34 <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Populations/PDB34/PDB34.ADALL_SNPS.txt', sep='\t', head=FALSE)
PDB48 <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Populations/PDB48/PDB48.ADALL_SNPS.txt', sep='\t', head=FALSE)

PDB48_c2 <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Clones/PDB48_c2/PDB48_c2.ADALL_SNPS.txt', sep='\t', head=FALSE)
PDB48_c13 <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Clones/PDB48_c13/PDB48_c13.ADALL_SNPS.txt', sep='\t', head=FALSE)
PDB48_c23 <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Clones/PDB48_c23/PDB48_c23.ADALL_SNPS.txt', sep='\t', head=FALSE)

PD68 <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Populations/PD68/PD68.ADALL_SNPS.txt', sep='\t', head=FALSE)
PD140 <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Populations/PD140/PD140.ADALL_SNPS.txt', sep='\t', head=FALSE)

E01 <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Populations/E01_S17/E01_S17.ADALL_SNPS.txt', sep='\t', head=FALSE)
E03 <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Populations/E03_S18/E03_S18.ADALL_SNPS.txt', sep='\t', head=FALSE)
E07 <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Populations/E07_S20/E07_S20.ADALL_SNPS.txt', sep='\t', head=FALSE)
E10 <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Populations/E10_S22/E10_S22.ADALL_SNPS.txt', sep='\t', head=FALSE)
E08 <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Populations/E08_S22/E08_S22.ADALL_SNPS.txt', sep='\t', head=FALSE)
E09 <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Populations/E09_S23/E09_S23.ADALL_SNPS.txt', sep='\t', head=FALSE)

WLP001 <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Populations/WLP001/WLP001_Preiss.ADALL_SNPS.txt', sep='\t', head=FALSE)
Drakes <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Populations/YMD3980_Pop/YMD3980_Pop.ADALL_SNPS.txt', sep='\t', head=FALSE)
RedCircle <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Populations/YMD4251/YMD4251.ADALL_SNPS.txt', sep='\t', head=FALSE)


PD68_CN <- read.table('/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Populations/PD68/CNV_new_1000bp/PD68_1000bp_norm.All.wig', sep='\t', head=FALSE)

subPDB1 <- subset(PDB1,V1==SubSetOn)
subPDB9 <- subset(PDB9,V1==SubSetOn)
subPDB24 <- subset(PDB24,V1==SubSetOn)
subPDB34 <- subset(PDB34,V1==SubSetOn)
subPDB48 <- subset(PDB48,V1==SubSetOn)

subPDB48_c2 <- subset(PDB48_c2,V1==SubSetOn)
subPDB48_c13 <- subset(PDB48_c13,V1==SubSetOn)
subPDB48_c23 <- subset(PDB48_c23,V1==SubSetOn)

subPD68 <- subset(PD68,V1=='8')
subPD140 <- subset(PD140,V1=='8')

subE01 <- subset(E01,V1=='8')
subE03 <- subset(E03,V1=='8')
subE07 <- subset(E07,V1=='8')
subE09 <- subset(E09,V1=='8')
subE10 <- subset(E10,V1=='8')

subDrakes <- subset(Drakes, V1=='8')
subRedCircle <- subset(RedCircle, V1=='8')
subWLP001 <- subset(WLP001, V1 =='8')

subPD68_CN <- subset(PD68_CN, V1 == '8')

GeneList <- read.table("/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/ChromosomeRegion_AllGenes_ChrVIII.tsv",
                       sep = '\t',
                       stringsAsFactors = FALSE,
                       col.names = c("Chrom",
                                     "Gene",
                                     "Type",
                                     "CommonName",
                                     "Start",
                                     "End",
                                     "Strand"))

Mutations <- read.table("/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Clones/PDB1_c1_by_PDB48_c2_GATK.vcf.table",
                        sep = '\t',
                        header = TRUE)

Mutations_Het <- read.table("/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/Populations/PDB1/PDB1.comb.filt.grep.chrVIII.table",
                        sep = '\t',
                        header = TRUE)

for(i in seq(1, length(GeneList[,1]))) {
  if(GeneList[i,]$CommonName == "") {
    GeneList[i,]$CommonName = GeneList[i,]$Gene
  }
}

GeneList <- GeneList[GeneList$Type == "ORF",]

##
## Allele Time-Series plots
##
AFTimeSeriesPlot <- ggplot() + 
  geom_point(data=subPDB1, aes(x=subPDB1$V2, y=subPDB1$V3, color="A", fill="A", alpha=0.5))+
  geom_point(data=subPDB9, aes(x=subPDB9$V2 + 600000, y=subPDB9$V3, color="B", fill="B", alpha=0.5))+
  geom_point(data=subPDB24, aes(x=subPDB24$V2 + 1200000, y=subPDB24$V3, color="C", fill="C", alpha=0.5))+
  geom_point(data=subPDB34, aes(x=subPDB34$V2 + 1800000, y=subPDB34$V3, color="D", fill="D", alpha=0.5))+
  geom_point(data=subPDB48, aes(x=subPDB48$V2 + 2400000, y=subPDB48$V3, color="E", fill="E", alpha=0.5))+
  #geom_vline(xintercept = 4811222, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 5411222, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 5411222 + 600000, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 5411222 + 1200000, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 5411222 + 1800000, linetype = "solid", color = "grey80") +
  #geom_vline(xintercept = 5392396 + 2400000 + 18826, linetype = "solid", color = "grey80") +
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
  scale_color_manual(values = darken(
    c("A"='#0066FF',
      "B"='#0080BF',
      "C"='#009980',
      "D"='#00B240',
      "E"='#00CC00')), 0.3) +
  scale_fill_manual(values = c("A"='#0066FF', "B"='#0080BF', "C"='#009980', "D"='#00B240', "E"='#00CC00')) +
  scale_y_continuous(name = "Allele Frequency",
                     expand = c(0.025, 0),
                     breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) +
  scale_x_continuous(name = "Genome Position",
                     expand = c(0.0125, 0.0125),
                     breaks = c(5111222,
                                5711222,
                                6311222,
                                6911222,
                                7511222),
                     labels = c("Repitch 0",
                                "Repitch 6",
                                "Repitch 15",
                                "Repitch 19",
                                "Repitch 26")) +
  theme(line = element_line(color="grey80", size=1),
        axis.line = element_line(color="grey80", size=1),
        axis.line.x = element_line(color="grey80", size=1),
        axis.line.y = element_line(color="grey80", size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(color = "grey80", size = 1),
        axis.ticks.x = element_blank(),
        text = element_text(family = "sans", color = "black", size=12),
        axis.text = element_text(family = "sans", color = "black", size=12),
        legend.position="none",
        panel.background = element_rect(fill = "white"))

PointSize = 1

##
## Clone plots
##

ClonePlot <- ggplot() +
  geom_rect(aes(xmin = 4830616 + 465717 + 1200000,
                xmax = 4830616 + 525315 + 1200000,
                ymin = -Inf,
                ymax = Inf), color = "grey", alpha = 0.2) +
  geom_point(data = subPDB1,
              aes(x = subPDB1$V2,
                  y = subPDB1$V3,
                  color = "Start",
                  fill = "Start"),
              alpha = 0.5,
              size = PointSize) +
  geom_point(data = subPDB1,
             aes(x = subPDB1$V2 + 600000,
                 y = subPDB1$V3,
                 color = "Start",
                 fill = "Start"),
             alpha = 0.5,
             size = PointSize) +
  geom_point(data = subPDB1,
             aes(x = subPDB1$V2 + 1200000,
                 y = subPDB1$V3,
                 color = "Start",
                 fill = "Start"),
             alpha = 0.5,
             size = PointSize) +
  geom_point(data = subPDB48_c2,
              aes(x = subPDB48_c2$V2,
                  y = subPDB48_c2$V3,
                  color = "Clone",
                  fill = "Clone"),
              alpha = 0.5,
             size = PointSize) +  
  geom_point(data = subPDB48_c13,
             aes(x = subPDB48_c13$V2 + 600000,
                 y = subPDB48_c13$V3,
                 color = "Clone",
                 fill = "Clone"),
             alpha = 0.5,
             size = PointSize) +  
  geom_point(data = subPDB48_c23,
             aes(x = subPDB48_c23$V2 + 1200000,
                 y = subPDB48_c23$V3,
                 color = "Clone",
                 fill = "Clone"),
             alpha = 0.5,
             size = PointSize) +  
  #geom_vline(xintercept = 4811222, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 5411222, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 5411222 + 600000, linetype = "solid", color = "grey80") +
  #geom_vline(xintercept = 5411222 + 1200000, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 4830616 + 465717 + 1200000, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 5198744 + 600000, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 173562 + 4829930, linetype = "dashed", color = "black") +
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
    values = darken(c("Clone" = '#00CC00',
                      "Start"='#0066FF',
                      "End"='#b3b3b3')), 0.3) +
  scale_fill_manual(
    name = NULL,
    values = c("Clone" = '#00CC00',
               "Start"='#0066FF',
               "End"='#b3b3b3')) +
  scale_y_continuous(name = "Allele Frequency",
                     expand = c(0.025, 0),
                     breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) +
  scale_x_continuous(name = "Genome Position",
                     expand = c(0.0125, 0.0125),
                     breaks = c(ChromMidPointList[as.numeric(SubSetOn)],
                                ChromMidPointList[as.numeric(SubSetOn)] + 600000,
                                ChromMidPointList[as.numeric(SubSetOn)] + 1200000),
                     labels = c("Clone 2",
                     "Clone 13",
                     "Clone 23")) +
  theme(line = element_line(color="grey80", size=1),
        axis.line = element_line(color="grey80", size=1),
        axis.line.x = element_line(color="grey80", size=1),
        axis.line.y = element_line(color="grey80", size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(color = "grey80", size = 1),
        axis.ticks.x = element_blank(),
        text = element_text(family = "sans", color = "black", size=12),
        axis.text = element_text(family = "sans", color = "black", size=12),
        legend.position="none",
        panel.background = element_rect(fill = "white"))

##
## Replicate Plots
##

ReplicationPlots <- ggplot() + 
  geom_point(data = subPD68,
             aes(x = subPD68$V2,
                 y = subPD68$V3,
                 color = "Clone",
                 fill = "Clone"),
             alpha = 0.5,
             size = PointSize) +
  geom_point(data = subE01,
             aes(x = subE01$V2 + 600000,
                 y = subE01$V3,
                 color = "Clone",
                 fill = "Clone"),
             alpha = 0.5,
             size = PointSize) +  
  geom_point(data = subE01,
             aes(x = subE01$V2 + 1200000,
                 y = subE01$V3,
                 color = "Clone",
                 fill = "Clone"),
             alpha = 0.5,
             size = PointSize) + 
  geom_point(data = subE07,
             aes(x = subE07$V2 + 1800000,
                 y = subE07$V3,
                 color = "Clone",
                 fill = "Clone"),
             alpha = 0.5,
             size = PointSize) +  
  geom_point(data = subWLP001,
             aes(x = subWLP001$V2 + 2400000,
                 y = subWLP001$V3,
                 color = "Clone",
                 fill = "Clone"),
             alpha = 0.5,
             size = PointSize) +
  geom_point(data = subWLP001,
             aes(x = subWLP001$V2 + 3000000,
                 y = subWLP001$V3,
                 color = "Clone",
                 fill = "Clone"),
             alpha = 0.5,
             size = PointSize) +
  geom_point(data = subPD140,
             aes(x = subPD140$V2,
                 y = subPD140$V3,
                 color = "Start",
                 fill = "Start"),
             alpha = 0.5,
             size = PointSize) +
  geom_point(data = subE03,
             aes(x = subE03$V2 + 600000,
                 y = subE03$V3,
                 color = "Start",
                 fill = "Start"),
             alpha = 0.5) +
  geom_point(data = subE09,
             aes(x = subE09$V2 + 1200000,
                 y = subE09$V3,
                 color = "Start",
                 fill = "Start"),
             alpha = 0.5,
             size = PointSize) + 
  geom_point(data = subE10,
             aes(x = subE10$V2 + 1800000,
                 y = subE10$V3,
                 color = "Start",
                 fill = "Start"),
             alpha = 0.5,
             size = PointSize) +  
  geom_point(data = subDrakes,
             aes(x = subDrakes$V2 + 2400000,
                 y = subDrakes$V3,
                 color = "Start",
                 fill = "Start"),
             alpha = 0.5,
             size = PointSize) +
  geom_point(data = subRedCircle,
             aes(x = subRedCircle$V2 + 3000000,
                 y = subRedCircle$V3,
                 color = "Start",
                 fill = "Start"),
             alpha = 0.5,
             size = PointSize) +
  #geom_vline(xintercept = 4811222, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 5411222, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 5411222 + 600000, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 5411222 + 1200000, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 5411222 + 1800000, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 5411222 + 2400000, linetype = "solid", color = "grey80") +
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
                     breaks = c(ChromMidPointList[as.numeric(SubSetOn)],
                                ChromMidPointList[as.numeric(SubSetOn)] + 600000,
                                ChromMidPointList[as.numeric(SubSetOn)] + 1200000,
                                ChromMidPointList[as.numeric(SubSetOn)] + 1800000,
                                ChromMidPointList[as.numeric(SubSetOn)] + 2400000,
                                ChromMidPointList[as.numeric(SubSetOn)] + 3000000),
                     labels = c("PDB Rep. 2",
                                "Elysian Rep. 1",
                                "Elysian Rep. 2",
                                "Elysian Rep. 3",
                                "Drake's",
                                "Red Circle")) +
  theme(line = element_line(color="grey80", size=1),
        axis.line = element_line(color="grey80", size=1),
        axis.line.x = element_line(color="grey80", size=1),
        axis.line.y = element_line(color="grey80", size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(color = "grey80", size = 1),
        axis.ticks.x = element_blank(),
        text = element_text(family = "sans", color = "black", size=12),
        axis.text = element_text(family = "sans", color = "black", size=12),
        legend.position="none",
        #        axis.text.x.bottom = element_blank(),
        #        axis.title.x.bottom = element_blank(),
        panel.background = element_rect(fill = "white"))

##
## PD68/140
##

PDPlots <- ggplot() + 
  geom_rect(aes(xmin = 4830616 + 465717,
                xmax = 4830616 + 525315,
                ymin = -Inf,
                ymax = Inf), color = "grey", alpha = 0.2) +
  geom_point(data = subPD68,
             aes(x = subPD68$V2,
                 y = subPD68$V3,
                 color = "Start",
                 fill = "Start"),
             alpha = 0.5) +
  geom_point(data = subPD140,
             aes(x = subPD140$V2,
                 y = subPD140$V3,
                 color = "End",
                 fill = "End"),
             alpha = 0.5) +  
  geom_vline(xintercept = 4830616, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 5392402, linetype = "longdash", color = "grey80") +
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
  scale_x_continuous(name = "Genome Position",
                     expand = c(0, 0),
                     breaks = ChromMidPointList[as.numeric(SubSetOn)],
                     labels = "Postdoc Brewing Co. Replicate 2") +
  theme(line = element_line(color="grey80", size=1),
        axis.line = element_line(color="grey80", size=1),
        axis.line.x = element_line(color="grey80", size=1),
        axis.line.y = element_line(color="grey80", size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(color = "grey80", size = 1),
        axis.ticks.x = element_blank(),
        text = element_text(family = "sans", color = "black", size=12),
        axis.text = element_text(family = "sans", color = "black", size=12),
        legend.position="none",
                axis.text.x.bottom = element_blank(),
                axis.title.x.bottom = element_blank(),
        panel.background = element_rect(fill = "white"))

##
## PD68/140 Copy Number
##

PD_CN_Plots <- ggplot() + 
  geom_rect(aes(xmin = 4830616 + 465717,
                xmax = 4830616 + 525315,
                ymin = -Inf,
                ymax = Inf), color = "grey", alpha = 0.2) +
  geom_point(data = subPD68_CN,
             aes(x = subPD68_CN$V2,
                 y = ((subPD68_CN$V3) / (subPD68_CN$V4)) * 4,
                 color = "Start",
                 fill = "Start"),
             alpha = 0.5) +
  geom_vline(xintercept = 4830616, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 5392402, linetype = "longdash", color = "grey80") +
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
    values = darken(c("End"='#33b4ff',
                      "Start"='#33b4ff')), 0.3) +
  scale_fill_manual(
    name = NULL,
    values = c("End"='#e69f00',
               "Start"='#33b4ff')) +
  scale_y_continuous(name = "Copy Number\n",
                     expand = c(0.025, 0),
                     breaks = c(0.0, 2.0, 4.0, 6.0),
                     limits = c(0.0, 6.5)) +
  scale_x_continuous(name = "Genome Position",
                     expand = c(0, 0),
                     breaks = ChromMidPointList[as.numeric(SubSetOn)],
                     labels = "Postdoc Brewing Co. Replicate 2") +
  theme(line = element_line(color="grey80", size=1),
        axis.line = element_line(color="grey80", size=1),
        axis.line.x = element_line(color="grey80", size=1),
        axis.line.y = element_line(color="grey80", size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(color = "grey80", size = 1),
        axis.ticks.x = element_blank(),
        text = element_text(family = "sans", color = "black", size=12),
        axis.text = element_text(family = "sans", color = "black", size=12),
        legend.position="none",
        #        axis.text.x.bottom = element_blank(),
        #        axis.title.x.bottom = element_blank(),
        panel.background = element_rect(fill = "white"))

##
## Gene Plot
##
GenePlot <- ggplot() + 
  geom_rect(aes(xmin = 465717,
                xmax = 525315,
                ymin = -Inf,
                ymax = Inf), color = "grey", fill = "grey", alpha = 0.2) +
  geom_hline(yintercept = 3.5) +
  geom_hline(yintercept = -3.5) +
  geom_label_repel(
    aes(x = (GeneList$Start + GeneList$End)/2, y = ifelse(GeneList$Strand > 0, GeneList$Strand / .222222222222, GeneList$Strand / .222222222222),
        label = GeneList$CommonName),
    color = "black",
    size = 7/.pt,
    direction = "x",
    nudge_y = ifelse(GeneList$Strand > 0, 4, -4),
    segment.alpha = 0.25,
    min.segment.length = 0.1,
    seed = 7654) +
  geom_rect(aes(xmin = GeneList$Start,
                xmax = GeneList$End,
                ymin = ifelse(GeneList$Strand > 0, 2.5, -2.5),
                ymax = ifelse(GeneList$Strand > 0, 4.5, -4.5)),
            color = "black",
            fill = "white") +
  geom_point(aes(x = (GeneList$Start + GeneList$End)/2,
                 y = ifelse(GeneList$Strand > 0, GeneList$Strand / .222222222222, GeneList$Strand / .222222222222)),
             color = "darkgrey",
             fill = "darkgrey") +
  geom_point(aes(x = Mutations_Het$POS,
                y = 0),
            shape = "*",
            size = 12,
            color = "darkblue",
#            width = 0,
#            height = 1,
            alpha = 0.5) + 
  geom_point(aes(x = Mutations$POS,
                y = 0),
            shape = "*",
            size = 12,
            color = "red") +
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
    values = darken(c("End"='#33b4ff',
                      "Start"='#33b4ff')), 0.3) +
  scale_fill_manual(
    name = NULL,
    values = c("End"='#e69f00',
               "Start"='#33b4ff')) +
  scale_y_continuous(name = "Candidate Region",
                     limits = c(-10, 10),
                     breaks = c(-3.5, 3.5),
                     labels = c(" -"," +")) +
  xlim(465717 - 1000, 525315 + 3500) +
  theme(line = element_line(color="grey80", size=1),
        axis.line = element_line(color="grey80", size=1),
        axis.line.x = element_line(color="grey80", size=1),
        axis.line.y = element_line(color="grey80", size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(color="grey80", size=1),
        axis.ticks.x = element_blank(),
        text = element_text(family = "sans", color = "black", size=12),
        axis.text = element_text(family = "sans", color = "black", size=12),
        legend.position="none",
        axis.text.x.bottom = element_blank(),
        axis.title.x.bottom = element_blank(),
        axis.text.y = element_text(color = "black", size=20, hjust = 0.5, vjust = 0.5),
        panel.background = element_rect(fill = "white"))


## Composition
pdf('//Volumes/GoogleDrive/My\ Drive/PostdocBrewing_Manuscript/Figures/AlleleFrequencyPlot_20200622.pdf',height=12, width=16, title="Post-Doc_Brewing")
ggdraw() +
  draw_plot(AFTimeSeriesPlot, x = 0, y = 0.60, height = 0.40, width = 1.0) +
  draw_plot(ClonePlot, x = 0.625, y = 0.2, height = 0.38, width = 0.375) +
  draw_plot(ReplicationPlots, x = 0, y = 0.2, height = 0.38, width = 0.625) +
  draw_plot(GenePlot, x = 0.0, y = 0.0, height = 0.18, width = 1.0) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                  x = c(0, 0, 0.625, 0.0), y = c(1, 0.6, 0.6, 0.20))
dev.off()

