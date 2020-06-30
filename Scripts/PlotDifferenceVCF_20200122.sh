#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Cris_L/BeerID/
#$ -o /net/dunham/vol2/Cris_L/BeerID/Scripts/OutStd/20200122/
#$ -e /net/dunham/vol2/Cris_L/BeerID/Scripts/ErrStd/20200122/
#$ -N Filter
#$ -l mfree=4G
#$ -l h_rt=2:0:0:0

module load modules modules-init modules-gs
module load java/8u25
module load bwa/0.7.15
module load picard/2.9.0
module load samtools/1.9
module load bedtools/2.28.0

SAMPLE_MAIN=$1
SAMPLE_FILT=$2

WORK_DIR=/net/dunham/vol2/Cris_L/BeerID/WorkDirectory_Combined # Give the desired location of intermediate output files.
SCRIPT_DIR=/net/dunham/vol2/Cris_L/BeerID/Scripts
REF=/net/dunham/vol2/Cris_L/ReferenceGenome/Reference_2Micron/sacCer3_2micron.fasta

OUT_NAME=CerLager_20191217
GATK_DIR=/net/dunham/vol2/Cris_L/Software/gatk-4.1.3.0
GATK_LIST=${SCRIPT_DIR}/GATK_Style_Interval_50kb.list

## Pull out variants from the large snp file 
${GATK_DIR}/gatk SelectVariants \
        -R ${REF} \
        --variant ${WORK_DIR}/${OUT_NAME}.snp.filt.vcf \
        --remove-unused-alternates true \
        --sample-name ${SAMPLE_MAIN} \
	-O ${WORK_DIR}/SubSample/${SAMPLE_MAIN}.vcf

## Extracts just variable sites
awk '{ if( $5 != ".") { print } }' < ${WORK_DIR}/SubSample/${SAMPLE_MAIN}.vcf > ${WORK_DIR}/SubSample/${SAMPLE_MAIN}.filt.vcf

# Filters by 'ancestor'
bedtools intersect \
	-header \
	-v \
	-a ${WORK_DIR}/SubSample/${SAMPLE_MAIN}.filt.vcf \
	-b ${WORK_DIR}/SubSample/${SAMPLE_FILT}.filt.vcf > ${WORK_DIR}/SubSample/${SAMPLE_MAIN}.filt.intersect.vcf

IGVTOOLS=/net/dunham/vol2/Caiti/hybrid_seq/IGVTools/igvtools.jar
java -Xmx2g -Djava.awt.headless=true -jar $IGVTOOLS index ${WORK_DIR}/SubSample/${SAMPLE_MAIN}.filt.intersect.vcf

# Variants to table
${GATK_DIR}/gatk VariantsToTable \
        -R $REF \
        -V ${WORK_DIR}/SubSample/${SAMPLE_MAIN}.filt.intersect.vcf \
        -F CHROM -F POS -GF AD \
        -O ${WORK_DIR}/SubSample/${SAMPLE_MAIN}.filt.intersect.table

# Move to allele frequency plot
java -Xmx2g -jar /net/dunham/vol2/Cris_L/Scripts/VaraiantTableParse_GATK_AF_10Reads_20190227.jar \
	${WORK_DIR}/SubSample/${SAMPLE_MAIN}.filt.intersect.table \
	${SCRIPT_DIR}/ChromSizes_BeerID.txt \
	${WORK_DIR}/SubSample/${SAMPLE_MAIN}_${SAMPLE_FILT}_

module load gcc/8.1.0
module load R/3.5.1

## Plot allele frequency REF is 1.0 and Alt is 0.0
Rscript ${SCRIPT_DIR}/PlotAlleleFreq_OneSample_SansDarken_20200122.R \
        ${WORK_DIR}/SubSample \
        ${SAMPLE_MAIN}_${SAMPLE_FILT}_${SAMPLE_MAIN} ## Sample name

