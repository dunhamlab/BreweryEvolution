##############
## FlavorProfillingAnalysis
## Chris Large
##############
## For measuring differences between panelist response to a beer sensory panel
##

## Load in necesary library
library(ggplot2)

## Load in data
flavorData <- read.delim('/Users/livinlrg/Google\ Drive/Dunham\ Lab/Postdoc_Brewing/FlavorProfiling/Beer_ScoreSheet_Reformatted.txt', 
                         sep='\t', head=TRUE, na.strings = "")

CharecteristicLevels <- levels(flavorData$Charecteristic)

flavorDataNAOmit<- flavorData[complete.cases(flavorData$Charecteristic),]

## Test Case
Astringency <- flavorDataNAOmit[flavorDataNAOmit$Charecteristic=="Astringency",]
AstringencyNAOmit<- Astringency[complete.cases(Astringency$Score_1),]

AstringencyNAOmit$Score_1 <- as.numeric(as.character(AstringencyNAOmit$Score_1))

AstringencyKruskal <- kruskal.test(Score_1 ~ Beer, data = AstringencyNAOmit)

TestCharecteristics <- c("Esters", "Phenols", "Alcohol", 
                         "Sweetness", "Phenols", "Harshness", 
                         "Warmth", "Acidity", "Hops",
                         "Bitterness", "Malt", "Clarity",
                         "Head Size", "Head Retention ")

## Loops through all mentioned test characteristics and does a kruskal.test
for(i in TestCharecteristics) {
  Charecteristic <- flavorDataNAOmit[flavorDataNAOmit$Charecteristic==i,]
  CharecteristicNAOmit<- Charecteristic[complete.cases(Charecteristic$Score_1),]
  CharecteristicNAOmit$Score_1 <- as.numeric(as.character(CharecteristicNAOmit$Score_1))
  print(i)
  print(kruskal.test(Score_1 ~ Beer, data = CharecteristicNAOmit))
  print(mean(CharecteristicNAOmit[CharecteristicNAOmit$Beer=="1",]$Score_1))
  print(mean(CharecteristicNAOmit[CharecteristicNAOmit$Beer=="2",]$Score_1))
  print(mean(CharecteristicNAOmit[CharecteristicNAOmit$Beer=="3",]$Score_1))
}

Flavors <- c("Grape", "Fruity","Dried Fruit", "Banana",
             "Citrus", "Stone Fruit", "Berry")

for(i in Flavors) {
  print(i)
  print(sum(flavorDataNAOmit[flavorDataNAOmit$Charecteristic==i,]$Score_1 == "Yes",na.rm=TRUE))
}
