## 
## Initialize
##

library(ggplot2)
library(colorspace) # for darken()
library(cowplot)
library(extrafont)
loadfonts()

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

SubSetOn = '5'

PDB1 <- read.table(paste(Dir, 'Populations/PDB1/CNV_new_1000bp/PDB1_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)
PDB9 <- read.table(paste(Dir, 'Populations/PDB9/CNV_new_1000bp/PDB9_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)
PDB24 <- read.table(paste(Dir, 'Populations/PDB24/CNV_new_1000bp/PDB24_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)
PDB34 <- read.table(paste(Dir, 'Populations/PDB34/CNV_new_1000bp/PDB34_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)
PDB48 <- read.table(paste(Dir, 'Populations/PDB48/CNV_new_1000bp/PDB48_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)

PD68 <- read.table(paste(Dir, 'Populations/PD68/CNV_new_1000bp/PD68_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)
PD140 <- read.table(paste(Dir, 'Populations/PD140/CNV_new_1000bp/PD140_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)

E01 <- read.table(paste(Dir, 'Populations/E01_S17/CNV_new_1000bp/E01_S17_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)
E03 <- read.table(paste(Dir, 'Populations/E03_S18/CNV_new_1000bp/E03_S18_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)
E07 <- read.table(paste(Dir, 'Populations/E07_S20/CNV_new_1000bp/E07_S20_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)
E10 <- read.table(paste(Dir, 'Populations/E10_S22/CNV_new_1000bp/E10_S22_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)
E08 <- read.table(paste(Dir, 'Populations/E08_S22/CNV_new_1000bp/E08_S22_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)
E09 <- read.table(paste(Dir, 'Populations/E09_S23/CNV_new_1000bp/E09_S23_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)

WLP001 <- read.table('/Users/livinlrg/Desktop/dunham/Cris_L/BeerID/WorkDirectory/Preiss_Dataset/SRR7406282_1000bp_norm.wig', sep='\t', head=FALSE)
Drake <- read.table(paste(Dir, 'Populations/YMD3980_Pop/CNV_new_1000bp/YMD3980_Pop_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)
RedCircle <- read.table(paste(Dir, 'Populations/YMD4251/CNV_new_1000bp/YMD4251_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)

subPDB1 <- subset(PDB1,V1==SubSetOn)
subPDB9 <- subset(PDB9,V1==SubSetOn)
subPDB24 <- subset(PDB24,V1==SubSetOn)
subPDB34 <- subset(PDB34,V1==SubSetOn)
subPDB48 <- subset(PDB48,V1==SubSetOn)

subPD68 <- subset(PD68,V1==SubSetOn)
subPD140 <- subset(PD140,V1==SubSetOn)

subE01 <- subset(E01,V1==SubSetOn)
subE03 <- subset(E03,V1==SubSetOn)
subE07 <- subset(E07,V1==SubSetOn)
subE10 <- subset(E10,V1==SubSetOn)
subE08 <- subset(E08,V1==SubSetOn)
subE09 <- subset(E09,V1==SubSetOn)

subWLP001 <- subset(WLP001,V1==SubSetOn)
subDrake <- subset(Drake,V1==SubSetOn)
subRedCircle <- subset(RedCircle,V1==SubSetOn)

PD68_AF <- read.table(paste(Dir, 'Populations/PD68/PD68.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)
PD140_AF <- read.table(paste(Dir, 'Populations/PD140/PD140.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)
  
PDB1_AF <- read.table(paste(Dir, 'Populations/PDB1/PDB1.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)
PDB48_AF <- read.table(paste(Dir, 'Populations/PDB48/PDB48.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)

PDB48_c1 <- read.table(paste(Dir, 'Clones/PDB48_c1/PDB48_c1.ADALL_SNPS.txt', sep = ""), sep='\t', head=FALSE)

subPD68_AF <- subset(PD68_AF,V1==SubSetOn)
subPD140_AF <- subset(PD140_AF,V1==SubSetOn)
subPDB1_AF <- subset(PDB1_AF,V1==SubSetOn)
subPDB48_AF <- subset(PDB48_AF,V1==SubSetOn)

subPDB48_c1 <- subset(PDB48_c1,V1==SubSetOn)

GeneList <- read.table("/Volumes/dunhamlab2/Cris_L/PostDocBrewing/WorkDirectory/ChromosomeRegion_AllGenes_AllChr.tsv",
                       sep = '\t',
                       stringsAsFactors = FALSE,
                       col.names = c("Chrom",
                                     "Gene",
                                     "Type",
                                     "CommonName",
                                     "Start",
                                     "End",
                                     "Strand"))

Shift = 600000

TimeCoursePlot <- ggplot() +
  geom_point(data=subPDB1, aes(x=V2 + (Shift*0), y=V3, color="A", fill="A", alpha=0.5))+
  geom_point(data=subPDB9, aes(x=V2 + (Shift*1), y=V3, color="B", fill="B", alpha=0.5))+
  geom_point(data=subPDB24, aes(x=V2 + (Shift*2), y=V3, color="C", fill="C", alpha=0.5))+
  geom_point(data=subPDB34, aes(x=V2 + (Shift*3), y=V3, color="D", fill="D", alpha=0.5))+
  geom_point(data=subPDB48, aes(x=V2 + (Shift*4), y=V3, color="E", fill="E", alpha=0.5))+
  geom_vline(xintercept = 2916956 + (Shift*1) - 42000, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 2916956 + (Shift*2) - 42000, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 2916956 + (Shift*3) - 42000, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 2916956 + (Shift*4) - 42000, linetype = "solid", color = "grey80") +
  annotate("point", x = 3175000, y = median(subPDB1$V3), color = "black", size = 3) +
#  annotate("text", x = 3175000  + (Shift * 0), y = 1.6, family = "sans", label = round(median(subPDB1$V3),2), color = "black") +
  annotate("point", x = 3175000 + Shift, y = median(subPDB9$V3), color = "black", size = 3) +
#  annotate("text", x = 3175000  + (Shift * 1), y = 1.6, family = "sans", label = round(median(subPDB9$V3),2), color = "black") +
  annotate("point", x = 3175000 + (Shift * 2), y = median(subPDB24$V3), color = "black", size = 3) +
#  annotate("text", x = 3175000  + (Shift * 2), y = 1.6, family = "sans", label = round(median(subPDB24$V3),2), color = "black") +
  annotate("point", x = 3175000 + (Shift * 3), y = median(subPDB34$V3), color = "black", size = 3) +
#  annotate("text", x = 3175000  + (Shift * 3), y = 1.6, family = "sans", label = round(median(subPDB34$V3),2), color = "black") +
  annotate("point", x = 3175000 + (Shift * 4), y = median(subPDB48$V3), color = "black", size = 3) +
#  annotate("text", x = 3175000  + (Shift * 4), y = 1.6, family = "sans", label = round(median(subPDB48$V3),2), color = "black") +
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
  scale_y_continuous(name = "Copy Number",
                     expand = c(0.025, 0),
                     breaks = c(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0),
                     limits = c(1.5, 5)) +
  scale_x_continuous(name = "Genome Position",
                     expand = c(0.0125, 0.0125),
                     breaks = c(3177956,
                                3177956+Shift,
                                3177956+Shift*2,
                                3177956+Shift*3,
                                3177956+Shift*4),
                     labels = c("Repitch 1\n\n",
                                "Repitch 6\n\n",
                                "Repitch 14\n\n",
                                "Repitch 18\n\n",
                                "Repitch 25\n\n")) +
  scale_color_manual(values = darken(
    c("A"='#0066FF',
      "B"='#0080BF',
      "C"='#009980',
      "D"='#00B240',
      "E"='#00CC00')), 0.3) +
  scale_fill_manual(values = c("A"='#0066FF',
                               "B"='#0080BF',
                               "C"='#009980',
                               "D"='#00B240',
                               "E"='#00CC00')) +
  theme(line = element_line(color="grey80", size=1),
        axis.line = element_line(color="grey80", size=1),
        axis.line.x.bottom = element_line(color="grey80", size=1),
        axis.line.y.left = element_line(color="grey80", size=1),
        axis.line.y.right = element_line(color="grey80", size=1),
        axis.line.x.top = element_line(color="grey80", size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(color = "grey80", size = 1),
        axis.ticks.x = element_blank(),
        text = element_text(family = "sans", color = "black", size=12),
        axis.text = element_text(family = "sans", color = "black", size=12),
        legend.position="none",
        panel.background = element_rect(fill = "white"))



WholeGenomePlot <- ggplot() +
  geom_vline(xintercept = 230218, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 1043402, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 1360022, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 2891955, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 3468829, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 3738990, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 4829930, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 5392573, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 5832461, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 6578212, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 7245028, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 8323205, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 9247636, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 10031969, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 11123260, linetype = "solid", color = "grey80") +
  geom_jitter(data=PD68,
              aes(x=V2,
                  y=V3,
                  shape='1',
                  size='1',                  
                  color = ifelse(V1 %% 2 == 0, "A", "B"),
                  fill = ifelse(V1 %% 2 == 0, "A", "B")),
              alpha = 0.5) +
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
    values=darken(c("A"='#33b4ff',
                    "B"='#0072b2',
                    "C"='#e69f00',
                    "D"='#ffbe33')), 0.3) +
  scale_fill_manual(
    name = NULL,
    values = c("A"='#33b4ff',
               "B"='#0072b2',
               "C"='#e69f00',
               "D"='#ffbe33')) +
  scale_y_continuous(name = "Copy Number",
                     expand = c(0.025, 0),
                     breaks = c(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0),
                     limits = c(0, 7.5)) +
  scale_x_continuous(name = "Genome Position",
                     expand = c(0.0, 0),
                     limits = c(-100000, 12171326),
                     breaks = c(-100000,
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
                                11597293),
                     labels = c("Chr.",
                                "I",
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
                                "XVI")) +
  theme(line = element_line(color="grey80", size=1),
        axis.line = element_line(color="grey80", size=1),
        axis.line.x.bottom = element_line(color="grey80", size=1),
        axis.line.y.left = element_line(color="grey80", size=1),
        axis.line.y.right = element_line(color="grey80", size=1),
        axis.line.x.top = element_line(color="grey80", size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(color = "grey80", size = 1),
        axis.ticks.x = element_blank(),
        text = element_text(family = "sans", color = "black", size=12),
        axis.text = element_text(family = "sans", color = "black", size=12),
        legend.position="none",
        panel.background = element_rect(fill = "white"))

ReplicatePlot <- ggplot() +
  geom_point(data=subPD68, aes(x=V2 + (Shift*0), y=V3, color="A", fill="A", alpha=0.5))+
  geom_point(data=subPD140, aes(x=V2 + (Shift*1), y=V3, color="D", fill="D", alpha=0.5))+
  geom_point(data=subE01, aes(x=V2 + (Shift*2), y=V3, color="A", fill="A", alpha=0.5))+
  geom_point(data=subE03, aes(x=V2 + (Shift*3), y=V3, color="D", fill="D", alpha=0.5))+
  geom_point(data=subWLP001, aes(x=V2 + (Shift*4), y=V3, color="A", fill="A", alpha=0.5))+
  geom_point(data=subDrake, aes(x=V2 + (Shift*5), y=V3, color="D", fill="D", alpha=0.5))+
  geom_point(data=subRedCircle, aes(x=V2 + (Shift*5), y=V3, color="D", fill="D", alpha=0.5))+
#  geom_point(data=subRedCircle, aes(x=V2 + (Shift*6), y=V3, color="D", fill="D", alpha=0.5))+
  scale_color_manual(
    name = NULL,
    values=darken(c("A"='#33b4ff',
                    "D"='#ffbe33')), 0.3) +
  scale_fill_manual(
    name = NULL,
    values = c("A"='#33b4ff',
               "D"='#ffbe33')) +
  geom_vline(xintercept = 2916956 + (Shift*1) - 42000, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 2916956 + (Shift*2) - 42000, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 2916956 + (Shift*3) - 42000, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 2916956 + (Shift*4) - 42000, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 2916956 + (Shift*5) - 42000, linetype = "solid", color = "grey80") +
  annotate("point", x = 3175000, y = median(subPD68$V3), color = "black", size = 3) +
#  annotate("text", x = 3175000, y = 1.6, family = "sans", label = round(median(subPD68$V3),2), color = "black") +
  annotate("point", x = 3175000 + Shift, y = median(subPD140$V3), color = "black", size = 3) +
#  annotate("text", x = 3175000  + (Shift * 1), y = 1.6, family = "sans", label = round(median(subPD140$V3),2), color = "black") +
  annotate("point", x = 3175000 + (Shift * 2), y = median(subE01$V3), color = "black", size = 3) +
#  annotate("text", x = 3175000  + (Shift * 2), y = 1.6, family = "sans", label = round(median(subE01$V3),2), color = "black") +
  annotate("point", x = 3175000 + (Shift * 3), y = median(subE03$V3), color = "black", size = 3) +
#  annotate("text", x = 3175000  + (Shift * 3), y = 1.6, family = "sans", label = round(median(subE03$V3),2), color = "black") +
  annotate("point", x = 3175000 + (Shift * 4), y = median(subWLP001$V3), color = "black", size = 3) +
#  annotate("text", x = 3175000  + (Shift * 4), y = 1.6, family = "sans", label = round(median(subWLP001$V3),2), color = "black") +
  annotate("point", x = 3175000 + (Shift * 5), y = median(subDrake$V3), color = "black", size = 3) +
#  annotate("text", x = 3175000  + (Shift * 5), y = 1.6, family = "sans", label = round(median(subDrake$V3),2), color = "black") +
  annotate("point", x = 3175000 + (Shift * 6), y = median(subRedCircle$V3), color = "black", size = 3) +
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
  scale_y_continuous(name = "Copy Number",
                     expand = c(0.025, 0),
                     breaks = c(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0),
                     limits = c(1.5, 5)) +
  scale_x_continuous(name = "Genome Position",
                     expand = c(0.0125, 0.0125),
                     breaks = c(3177956,
                                3177956+Shift,
                                3177956+Shift*2,
                                3177956+Shift*3,
                                3177956+Shift*4,
                                3177956+Shift*5
#                                ,3177956+Shift*6
                                ),
                     labels = c("Repitch 1\n\n",
                                "Repitch 29\n\n",
                                "Repitch 1\n\n",
                                "Repitch 15\n\n",
                                "WLP001\n\n",
                                "Repitch 24\n\n"
                                ,"Repitch 36\n\n"
                                )) +
  #coord_cartesian(clip = "off") +
  theme(line = element_line(color="grey80", size=1),
        axis.line = element_line(color="grey80", size=1),
        axis.line.x.bottom = element_line(color="grey80", size=1),
        axis.line.y.left = element_line(color="grey80", size=1),
        axis.line.y.right = element_line(color="grey80", size=1),
        axis.line.x.top = element_line(color="grey80", size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(color = "grey80", size = 1),
        axis.ticks.x = element_blank(),
        text = element_text(family = "sans", color = "black", size=12),
        axis.text = element_text(family = "sans", color = "black", size=12),
        legend.position="none",
        panel.background = element_rect(fill = "white"))

PointSize = 2

AFPlot <- ggplot() + 
  geom_point(data = subPDB1_AF,
             aes(x = V2,
                 y = V3,
                 color = "Clone",
                 fill = "Clone"),
             alpha = 0.5,
             size = PointSize) +
  geom_point(data = subPDB48_c1,
             aes(x = V2 + 600000,
                 y = V3,
                 color = "Clone",
                 fill = "Clone"),
             alpha = 0.5,
             size = PointSize) +  
  geom_point(data = subPDB48_AF,
             aes(x = V2,
                 y = V3,
                 color = "Start",
                 fill = "Start"),
             alpha = 0.5,
             size = PointSize) +
  #geom_point(data = subE03,
  #           aes(x = V2 + 600000,
  #               y = V3,
  #               color = "Start",
  #               fill = "Start"),
  #           alpha = 0.5) +
  #geom_vline(xintercept = 4811222, linetype = "longdash", color = "grey80") +
  geom_vline(xintercept = 2916956 + (Shift*1) - 42000, linetype = "solid", color = "grey80") +
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
                                ChromMidPointList[as.numeric(SubSetOn)] + 600000),
                     labels = c("\n\n\n\n\n",
                                "\n\n\n\n\n")) +
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

pdf('//Volumes/GoogleDrive/My\ Drive/PostdocBrewing_Manuscript/Figures/CopyNumberPlot_20200617.pdf',height=12, width=12, title="Post-Doc_Brewing")
ggdraw() +
  draw_plot(WholeGenomePlot, x = 0, y = 0.70, height = 0.30, width = 1.0) +
  draw_plot(TimeCoursePlot, x = 0, y = 0.4, height = 0.28, width = 0.45) +
  draw_plot(ReplicatePlot, x = 0.45, y = 0.4, height = 0.28, width = 0.55) +
  draw_plot(AFPlot, x = 0.0, y = 0, height = 0.36, width = 0.65) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                  x = c(0, 0, 0.45, 0.0), y = c(1, 0.7, 0.7, 0.40))
dev.off()

