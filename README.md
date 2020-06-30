# Brewery Evolution

These are the scripts and data that were used to generate the figures presented in:

** Genomic stability and adaptation of beer brewing yeasts during serial repitching in the brewery
* * Christopher R. L. Large, Noah A Hanson, Andreas Tsouris, Omar Abou Saada, Jirasin Koonthongkaew, Yoichi Toyokawa, Tom Schmidlin, Daniela A Moreno-Habel, Hal McConnellogue, Richard Preiss, Hiroshi Takagi, Joseph Schacherer, Maitreya J Dunham

doi: https://doi.org/10.1101/2020.06.26.166157

## Table of Contents
1. Data
- (A) Plate Reader
- (B) Settling Assay
- (C) Sensory Data
2. Scripts
- (A) Alignment and Variant Calling
- (B) Copy Number Analysis
- (C) Allele Frequency


# Scripts
Please note that in many instances, the scripts outlined here are written for my computing cluster and will require some retooling if they are to be adapted for other uses. Please bear with my occasional hard coding of directories. Hopefully, they can serve as inspiration for further studies. 

## (A) Alignment and Variant Calling
###### runSeqAlignVariantCall_20200504.sh
- First, the reads mentioned in the paper were aligned to the SacCer3 genome from SGD. They were then preprocessed with GATK and picard tools
- Second, variant calls for the ancestor were generated using Samtools and Freebayes.
- Third, with the ancestor variant calls in hand, the evolved samples were called with LoFreq in paired mode and the variant calls from Samtools and Freebayes were compared to the ancestor.
- Fourth, the variants were filetered and annotated using: yeast_annotation_chris_edits_20170925.py
- Finally, the variants found in the clones were comapred using MakeOverlapMatrix.R to generate an overlap matrix.

## (B) Copy Number Analysis
######  runCNScripts_Populations_20200310.sh
- First, the alignments from part (A) were measured for their genome wide depth of coverage using GATK/2. Please note that there is a new tool in GATK/4 that serves the same function. Next, 1000-bp sliding windows of coverage were generated with the IGVtools command line implementation. 
Copy_Number_Plot_Text_20200619.R
######  Copy_Number_Population_AllChrom_20200609.R


## (C) Allele Frequency
######  runVariantCall_GATK_Populations_20200310.sh >
######  VaraiantTableParse_GATK_BAF_20190206.jar
######  AlleleFrequency_Population_AllChrom_20200429.R
######  Allele_Frequency_Clones_Plot_2020428.R




