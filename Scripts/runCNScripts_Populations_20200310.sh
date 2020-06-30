#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Cris_L/PostDocBrewing
#$ -e /net/dunham/vol2/Cris_L/PostDocBrewing/Scripts/ErrStd/20200109/
#$ -o /net/dunham/vol2/Cris_L/PostDocBrewing/Scripts/OutStd/20200109/
#$ -N BamToCN
#$ -l mfree=4G

## CNV pipeline for figuring out CN from bam files: uses both wig file from igvtools and mpileup

module load modules modules-init modules-gs
module load python/2.7.3 
module load numpy/1.7.0 #latest, needs at least python 2.7.3
module load java/7u17 java_jdk/7u17
module load GATK/2.6.5
IGVTOOLS=/net/dunham/vol2/Caiti/hybrid_seq/IGVTools/igvtools.jar

#cer
SAMPLE=$1 #sample prefix (ex: Sample-01)
SIZE=$2
WORKDIR=/net/dunham/vol2/Cris_L/PostDocBrewing/WorkDirectory/Populations #CHANGE
BAMDIR=${WORKDIR}/${SAMPLE} #CHANGE
CNDIR=${WORKDIR}/${SAMPLE}/CNV_new_${SIZE}bp #CHANGE
SCRIPTS=/net/dunham/vol2/Caiti/scripts
REF=/net/dunham/vol2/Caiti/reference_seq/sacCer3.fasta
#ANC=$3
CERORF=/net/dunham/vol2/Caiti/reference_seq/cer_homolog_coordinates_filt.txt

cd ${WORKDIR}/${SAMPLE}
mkdir -p CNV_new_${SIZE}bp

## Get depth of coverage info (old version ran this on earlier bam). Ignores mito
java -Xmx2g -jar $GATK_DIR/GenomeAnalysisTK.jar -T DepthOfCoverage \
	-R ${REF} -I ${BAMDIR}/${SAMPLE}_comb_R1R2.RG.MD.realign.sort.bam \
	-o ${CNDIR}/${SAMPLE}_comb_R1R2.dedup.RG.sort.bam_DOC \
	-XL chrM -omitBaseOutput -omitLocusTable -omitIntervals -rf BadCigar

## Make wig file (can also run igvtools directly, but java allows mem management)
## Can change window size
## NOTE: RUN ON ANCESTRAL FIRST
java -Xmx2g -Djava.awt.headless=true -jar $IGVTOOLS count -w ${SIZE} --minMapQuality 30 \
	${BAMDIR}/${SAMPLE}_comb_R1R2.RG.MD.realign.sort.bam \
	${CNDIR}/${SAMPLE}_${SIZE}bp.wig ${REF}

java -Xmx2g -Djava.awt.headless=true -jar $IGVTOOLS count -w ${SIZE} --minMapQuality 0 \
        ${BAMDIR}/${SAMPLE}_comb_R1R2.RG.MD.realign.sort.bam \
        ${CNDIR}/${SAMPLE}_${SIZE}bp.All.wig ${REF}

python /net/dunham/vol2/Cris_L/BeerID/Scripts/wigNormalizedToAverageReadDepth_MapQ_ForPlot.py       ${CNDIR}/${SAMPLE}_comb_R1R2.dedup.RG.sort.bam_DOC.sample_summary         ${CNDIR}/${SAMPLE}_${SIZE}bp.wig    ${CNDIR}/${SAMPLE}_${SIZE}bp.All.wig     4 ${CNDIR}/${SAMPLE}_${SIZE}bp_norm.wig

python /net/dunham/vol2/Cris_L/BeerID/Scripts/wigNormalizedToAverageReadDepth_MapQ_ForPlot.py \
        ${CNDIR}/${SAMPLE}_comb_R1R2.dedup.RG.sort.bam_DOC.sample_summary \
        ${CNDIR}/${SAMPLE}_${SIZE}bp.wig \
	${CNDIR}/${SAMPLE}_${SIZE}bp.All.wig
	4 \
        ${CNDIR}/${SAMPLE}_${SIZE}bp_norm.wig

module load gcc/8.1.0
module load R/latest

Rscript /net/dunham/vol2/Cris_L/PostDocBrewing/Scripts/PlotCopyNumber_OneSample_20200106.R ${WORKDIR} ${SAMPLE}

