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


## Location of data files from plate reader
Dir="/Users/livinlrg/Google Drive/Dunham Lab/Postdoc_Brewing/PlateReader/"

growthData_06_08_YEPD_WORT <- read.csv(
  paste(Dir,
        "2018_05_08_GrowthCurve_YEPD_WORT/",
        "2018_06_08_GrowthCurve_YEPD_WORT.csv", sep ="")
  ,sep=",", header=TRUE, check.names = FALSE,stringsAsFactors=FALSE)

growthData_05_09_EtOH <- read.csv(
  paste(Dir,
        "2018_05_09_PDB_EtOh/",
        "2018_05_09_PDB_EtOh_Formatted.csv", sep ="")
  ,sep=",", header=TRUE, check.names = FALSE,stringsAsFactors=FALSE)

growthData_05_21_Cold <- read.csv(
  paste(Dir,
        "2018_05_21_OutOfCold/",
        "2018_06_21_OutOfCold_48hrRun_All.csv", sep ="")
  ,sep=",", header=TRUE, check.names = FALSE,stringsAsFactors=FALSE)

growthData_06_05_EtOH <- read.csv(
  paste(Dir,
        "2018_06_05_TestForEthanolInWort/",
        "2018_06_05_TestForEthanolInWort.csv", sep ="")
  ,sep=",", header=TRUE, check.names = FALSE,stringsAsFactors=FALSE)

growthData_06_13_EtOH <- read.csv(
  paste(Dir,
        "2018_06_13_TestForEthanolInWort/",
        "2018_06_13_TestForEthanolInWort.csv", sep ="")
  ,sep=",", header=TRUE, check.names = FALSE,stringsAsFactors=FALSE)

growthData_08_01_Cold <- read.csv(
  paste(Dir,
        "2018_08_01_Biological_Out_Of_Cold/",
        "2018_08_01_Biological_Out_Of_Cold.csv", sep ="")
  ,sep=",", header=TRUE, check.names = FALSE,stringsAsFactors=FALSE)

GrowthList <- c('growthData_06_08_YEPD_WORT',
                'growthData_05_09_EtOH',
                'growthData_05_21_Cold',
                'growthData_06_05_EtOH',
                'growthData_06_13_EtOH',
                'growthData_08_01_Cold')

## Function to melt the data to time as a column value with OD600 values
AssignVariables <- function(GrowthData) {
  GrowthData.melt <- melt(GrowthData,id.vars = c("Strain","Replicate","Well","Condition"))
  colnames(GrowthData.melt)[5:6] <- c("Time","OD600")
  
  GrowthData.melt$Time <- as.numeric(as.character(GrowthData.melt$Time))
  GrowthData.melt <- na.omit(GrowthData.melt)
  
  GrowthData.melt <- GrowthData.melt[which(GrowthData.melt$Strain != "Blank"),]
  return(GrowthData.melt)
}

## Function to loop through the samples and replictates to generate linear fits to the exponential phase
GrowthCurve <- function(GrowthData.melt, growthData) {
  ResultsTable <- data.frame()
  for(i in levels(factor(GrowthData.melt$Strain))) {
    for(j in levels(factor(GrowthData.melt[GrowthData.melt$Strain == i & GrowthData.melt$Condition == '0',]$Replicate))) {
      GrowthData.melt.update <- GrowthData.melt[GrowthData.melt$Strain == i & GrowthData.melt$Condition == 0 & GrowthData.melt$Replicate == j,c(1,2,4,5,6)]
      
      #print(summary(GrowthData.melt.update$OD600))
      print(paste(i, " ", j, sep=""))
      
      print(max(GrowthData.melt.update$OD600))
      
      fit <- fit_easylinear(GrowthData.melt.update$Time, GrowthData.melt.update$OD600)
      print(coef(fit)[3])
      print(plot(fit, log = "y"))

      ResultsTable <- rbind(ResultsTable, data.frame(i, j, coef(fit)[1],coef(fit)[2],coef(fit)[3], coef(fit)[4], rsquared(fit)))
     }
  }
  return(ResultsTable)
}

## Cycles through the experiments
for(i in GrowthList) {
  print(i)
  assign(paste(i, ".table", sep=""),GrowthCurve(AssignVariables(get(i)),get(i)))
}

## Gives a new column with the sample name
growthData_05_09_EtOH.table$Sample <- 'EtOH_05_09'
growthData_06_05_EtOH.table$Sample <- 'EtOH_06_05'
growthData_06_13_EtOH.table$Sample <- 'EtOH_06_13'
growthData_06_08_YEPD_WORT.table$Sample <- 'EtOH_06_08'
growthData_05_21_Cold.table$Sample <- 'Cold_05_21'
growthData_08_01_Cold.table$Sample <- 'Cold_08_01'

## Combines the samples that are conducted simiarily on the same samples
Combined.Table <- rbind(
      growthData_06_05_EtOH.table,
      growthData_06_13_EtOH.table,
      growthData_06_08_YEPD_WORT.table[growthData_06_08_YEPD_WORT.table$i == "PDB1_c1" | 
                                         growthData_06_08_YEPD_WORT.table$i == "PDB1_c2" | 
                                         growthData_06_08_YEPD_WORT.table$i == "PDB48_c1" | 
                                         growthData_06_08_YEPD_WORT.table$i == "PDB48_c2" | 
                                         growthData_06_08_YEPD_WORT.table$i == "PDB48_c6",])

pdf('//Volumes/GoogleDrive/My\ Drive/PostdocBrewing_Manuscript/Figures/GrowthRates_20200619.pdf',height=4, width=6, title="Post-Doc_Brewing")
ggplot(Combined.Table, aes(x=i, y=coef.fit..3.) ) +
  geom_boxplot(width = 0.4, outlier.size = 1) +
  geom_dotplot(binaxis = 'y',
               stackratio=1,
               stackdir='center',
               dotsize = 0.25) +
  scale_y_continuous(name = expression(paste("Growth Rate (", min^{-1}, ")")),
                     limits = c(0.0014, 0.00185)) +
  scale_x_discrete(labels = c("PDB1\nClone 1\nAncestor",
                              "PDB1\nClone 2\nAncestor",
                              "PDB26\nClone 1\nCNV Chr. V",
                              "PDB26\nClone 2\nLOH Chr. VIII",
                              "PDB26\nClone 6\nLOH Chr. VIII")) +
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
        panel.background = element_rect(fill = "white"))
dev.off()

## Tests for any differences
kruskal.test(coef.fit..3. ~ i,data = Combined.Table)

## Tests for specific differences
wilcox.test(Combined.Table[Combined.Table$i == "PDB1_c1",]$coef.fit..3.,
       Combined.Table[Combined.Table$i == "PDB48_c2",]$coef.fit..3.)

