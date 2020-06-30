setwd("/Users/Chris/Google Drive/Dunham Lab/Postdoc_Brewing/Mutations/")

MutationData <- read.csv("SNP_INDEL_ALL.csv",sep=",", header=TRUE, stringsAsFactors = FALSE)

Clones <- levels(as.factor(MutationData$Sample))

OverlapMatrix <- matrix(, nrow = length(Clones), ncol = length(Clones))
iValue = 1
jValue = 1
for(i in Clones) {
  for(j in Clones) {
    print(paste(i,j,iValue,jValue))
    OverlapMatrix[iValue,jValue] = length(c(MutationData$Gene[which(MutationData$Sample==i)],
                                            MutationData$Gene[which(MutationData$Sample==j)])) - length(unique(c(MutationData$Gene[which(MutationData$Sample==i)],
                                                                                                                 MutationData$Gene[which(MutationData$Sample==j)])))
    jValue = jValue + 1
  }
  jValue = 1
  iValue = iValue + 1
}

colnames(OverlapMatrix) = Clones
rownames(OverlapMatrix) = Clones

write.csv(OverlapMatrix,file="OverlapMatrix.csv")
