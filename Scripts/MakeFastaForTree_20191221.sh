#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Cris_L/BeerID/
#$ -o /net/dunham/vol2/Cris_L/BeerID/Scripts/OutStd/20191203/
#$ -e /net/dunham/vol2/Cris_L/BeerID/Scripts/ErrStd/20191203/
#$ -N FastaForTree
#$ -l mfree=8G
#$ -l h_rt=10:0:0:0

module load modules modules-init modules-gs
module load java/8u25
module load bwa/0.7.15
module load GATK/4.1.1.0
module load picard/2.9.0
module load samtools/1.9
module load htslib/1.9
module load bcftools/1.9
module load tabix/0.2.6

## Pulls relevant info from the sample sheet that has been submitted with runQsub_SampleList.py
SAMPLE=$1
LIBRARY_NAME=$2
COMMON_NAME=$3
PLOIDY=$4
DATASET=$5
WORK_DIR=/net/dunham/vol2/Cris_L/BeerID/WorkDirectory/Fasta/AmericanBeer_NoBE051_20200430
FASTQ_DIR=/net/dunham/vol2/Cris_L/BeerID/Fastqs/${DATASET}_Dataset
## Aligns to the S288C reference genome with the 2-Micron plasmid added
REF=/net/dunham/vol2/Cris_L/ReferenceGenome/Reference_2Micron/sacCer3_2micron.50kb.fasta

VCF=/net/dunham/vol2/Cris_L/BeerID/WorkDirectory_Combined/AmericanBeer_NoBE051_20200430.snp.filt.BE051_filt.vcf

cd ${WORK_DIR}

## Update the refernece with the reference allele at heterozygous sites
bcftools consensus -f ${REF} \
	--haplotype R \
	-s ${COMMON_NAME} \
	-o ${WORK_DIR}/${COMMON_NAME}.ref.fasta \
	${VCF}.gz

## Update the reference with the alternate allele at heterozyous sites
bcftools consensus -f ${REF} \
        --haplotype A \
        -s ${COMMON_NAME} \
        -o ${WORK_DIR}/${COMMON_NAME}.alt.fasta \
        ${VCF}.gz

grep -v '>' ${WORK_DIR}/${COMMON_NAME}.ref.fasta | tr -d '\n' > ${WORK_DIR}/${COMMON_NAME}.ref.noArrow.fasta

grep -v '>' ${WORK_DIR}/${COMMON_NAME}.alt.fasta | tr -d '\n' > ${WORK_DIR}/${COMMON_NAME}.alt.noArrow.fasta

echo \>${COMMON_NAME} | cat - ${WORK_DIR}/${COMMON_NAME}.ref.noArrow.fasta ${WORK_DIR}/${COMMON_NAME}.alt.noArrow.fasta > ${COMMON_NAME}.cat.fasta

rm ${WORK_DIR}/${COMMON_NAME}.ref.fasta
rm ${WORK_DIR}/${COMMON_NAME}.alt.fasta
rm ${WORK_DIR}/${COMMON_NAME}.ref.noArrow.fasta
rm ${WORK_DIR}/${COMMON_NAME}.alt.noArrow.fasta

