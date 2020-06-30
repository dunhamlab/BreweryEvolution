
library(ggplot2)
library(zoo)

PDB1 <- read.table('/Users/livinlrg/Desktop/dunham/Cris_L/PostDocBrewing/WorkDirectory/Populations/PDB1/PDB1.ADALL_SNPS.txt', sep='\t', head=FALSE)
PDB9 <- read.table('/Users/livinlrg/Desktop/dunham/Cris_L/PostDocBrewing/WorkDirectory/Populations/PDB9/PDB9.ADALL_SNPS.txt', sep='\t', head=FALSE)
PDB24 <- read.table('/Users/livinlrg/Desktop/dunham/Cris_L/PostDocBrewing/WorkDirectory/Populations/PDB24/PDB24.ADALL_SNPS.txt', sep='\t', head=FALSE)
PDB34 <- read.table('/Users/livinlrg/Desktop/dunham/Cris_L/PostDocBrewing/WorkDirectory/Populations/PDB34/PDB34.ADALL_SNPS.txt', sep='\t', head=FALSE)
PDB48 <- read.table('/Users/livinlrg/Desktop/dunham/Cris_L/PostDocBrewing/WorkDirectory/Populations/PDB48/PDB48.ADALL_SNPS.txt', sep='\t', head=FALSE)

PD68 <- read.table('/Users/livinlrg/Desktop/dunham/Cris_L/PostDocBrewing/WorkDirectory/Populations/PD68/PD68.ADALL_SNPS.txt', sep='\t', head=FALSE)
PD140 <- read.table('/Users/livinlrg/Desktop/dunham/Cris_L/PostDocBrewing/WorkDirectory/Populations/PD140/PD140.ADALL_SNPS.txt', sep='\t', head=FALSE)

PDB48_c2 <- read.table('/Users/livinlrg/Desktop/dunham/Cris_L/PostDocBrewing/WorkDirectory/Clones/PDB48_c2/PDB48_c2.ADALL_SNPS.txt', sep='\t', head=FALSE)
PDB48_c13 <- read.table('/Users/livinlrg/Desktop/dunham/Cris_L/PostDocBrewing/WorkDirectory/Clones/PDB48_c13/PDB48_c13.ADALL_SNPS.txt', sep='\t', head=FALSE)
PDB48_c20 <- read.table('/Users/livinlrg/Desktop/dunham/Cris_L/PostDocBrewing/WorkDirectory/Clones/PDB48_c20/PDB48_c20.ADALL_SNPS.txt', sep='\t', head=FALSE)
PDB48_c23 <- read.table('/Users/livinlrg/Desktop/dunham/Cris_L/PostDocBrewing/WorkDirectory/Clones/PDB48_c23/PDB48_c23.ADALL_SNPS.txt', sep='\t', head=FALSE)

WLP001 <- read.table('/Users/livinlrg/Desktop/dunham/Cris_L/PostDocBrewing/WorkDirectory/Populations/WLP001/WLP001_Preiss.ADALL_SNPS.txt', sep='\t', head=FALSE)

## Positions that reach 0.5 freqeuncy between the four clones
PDB48_clones <- Reduce(function(x, y)
  merge(x, y, 
        by=c('V1','V2'),
        all=FALSE),
  list(PDB48_c2,
       PDB48_c13,
       PDB48_c20))

subPDB48_clones <- subset(PDB48_clones,V1=='8' &
                            V2 > 5240000 &
                            V2 < 5352500)

subPDB1 <- subset(PDB1,V1=='8' &
                    V2 > 5240000 &
                    V2 < 5352500)

## Plot every clone in region
ggplot() +
  geom_point(aes(x = subPDB48_clones$V2, y = subPDB48_clones$V3.x), color = "red", fill = "red") +
  geom_point(aes(x = subPDB48_clones$V2, y = subPDB48_clones$V3.y), color = "green", fill = "green") +
  geom_point(aes(x = subPDB48_clones$V2, y = subPDB48_clones$V3), color = "blue", fill = "blue") +
  theme_minimal()

## Plot average of these clones
plot(x = subPDB48_clones$V2,
     y = (subPDB48_clones$V3.x + subPDB48_clones$V3.y + subPDB48_clones$V3)/3)

## Create df of the average
subPDB48_clones_Values <- data.frame(subPDB48_clones$V1,
                             subPDB48_clones$V2,
                             (subPDB48_clones$V3.x + subPDB48_clones$V3.y + subPDB48_clones$V3)/3)

## Use where the values are within the 0.5 range as a selection for the best positions to estimate proportion
PDB48_r8 <- PDB48_c23[PDB48_c23$V2 %in% subPDB48_clones[0.54 > subPDB48_clones_Values[,3] &
                                                  subPDB48_clones_Values[,3] > 0.46,]$V2,]

plot(PDB48_r8$V2, PDB48_r8$V3)

## Function 
ValuesOfAF <- function(Clone) {
  return(Clone[(Clone$V2 %in% (subPDB48_clones[0.54 > subPDB48_clones_Values[,3] & subPDB48_clones_Values[,3] > 0.46,]$V2))
               & (Clone$V2 %in% subPDB1[subPDB1$V3 > 0.55,]$V2 
               | Clone$V2 %in% subPDB1[subPDB1$V3 < 0.45,]$V2),])
}

## Assume AF going to 0.5 to determine frequency in population
## 0.5 - 0.75 * 4 = 1
## 0.5 - 0.5 * 4 = 0
ThingsToMeasure <- c("PDB1", "PDB9", "PDB24", "PDB34", "PDB48", "PDB48_c2", "PDB48_c13", "PDB48_c20", "PDB48_c23", "WLP001","PD68", "PD140")

for(i in ThingsToMeasure) {
  plot(ValuesOfAF(get(i))$V2, ValuesOfAF(get(i))$V3, xlab = i)
  print(i)
  print(mean(ValuesOfAF(get(i))$V3))
  print(1 - (( (abs(0.5 - mean(ValuesOfAF(get(i))$V3[ValuesOfAF(get(i))$V3 > 0.50]))+
                abs(0.5 - mean(ValuesOfAF(get(i))$V3[ValuesOfAF(get(i))$V3 < 0.50])))/2) *4))
}

PDB_TimeCourse <- c("PDB1", "PDB9", "PDB24", "PDB34", "PDB48")
PDB_Values <- data.frame(c(0, 18, 45, 57, 78),c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0))
j = 1

for(i in PDB_TimeCourse) {
  PDB_Values[j,2] <- 1 - (( (abs(0.5 - mean(ValuesOfAF(get(i))$V3[ValuesOfAF(get(i))$V3 > 0.50]))+
                  abs(0.5 - mean(ValuesOfAF(get(i))$V3[ValuesOfAF(get(i))$V3 < 0.50])))/2) *4)
  j <- j + 1
}

j = 1
for(i in PDB_Values[,2]) {
  PDB_Values[j,3] <-(log(i/(1-i)))
  j <- j + 1
}

plot(x = PDB_Values[,1], y=PDB_Values[,3], xlab = "Generations", ylab = "Ln(A/a)")
abline(lm(PDB_Values[,3] ~ PDB_Values[,1]), col = "red")
exp(lm(PDB_Values[,3] ~ PDB_Values[,1])$coefficients[2]) - 1

### Anueploidy
Dir="/Users/livinlrg/Desktop/dunham/Cris_L/PostDocBrewing/WorkDirectory/"
PDB1 <- read.table(paste(Dir, 'Populations/PDB1/CNV_new_1000bp/PDB1_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)
PDB9 <- read.table(paste(Dir, 'Populations/PDB9/CNV_new_1000bp/PDB9_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)
PDB24 <- read.table(paste(Dir, 'Populations/PDB24/CNV_new_1000bp/PDB24_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)
PDB34 <- read.table(paste(Dir, 'Populations/PDB34/CNV_new_1000bp/PDB34_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)
PDB48 <- read.table(paste(Dir, 'Populations/PDB48/CNV_new_1000bp/PDB48_1000bp_norm.wig', sep = ""), sep='\t', head=FALSE)

PDB_Values <- data.frame(c(0, 18, 45, 57, 78),c(0, 0, 0, 0, 0), c(0, 0, 0, 0, 0))

PDB_Values[1,2] <- mean(PDB1[PDB1$V1 == 5,]$V3) -3
PDB_Values[2,2] <- mean(PDB9[PDB9$V1 == 5,]$V3) -3
PDB_Values[3,2] <- mean(PDB24[PDB24$V1 == 5,]$V3) -3
PDB_Values[4,2] <- mean(PDB34[PDB34$V1 == 5,]$V3) -3
PDB_Values[5,2] <- mean(PDB48[PDB48$V1 == 5,]$V3) -3

j = 1
for(i in PDB_Values[,2]) {
  PDB_Values[j,3] <-(log(i/(1-i)))
  j <- j + 1
}

plot(x = PDB_Values[,1], y=PDB_Values[,3], xlab = "Generations", ylab = "Ln(A/a)")
abline(lm(PDB_Values[,3] ~ PDB_Values[,1]), col = "red")
exp(lm(PDB_Values[,3] ~ PDB_Values[,1])$coefficients[2]) - 1

