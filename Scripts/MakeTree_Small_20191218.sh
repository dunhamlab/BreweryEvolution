#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Cris_L/BeerID/
#$ -o /net/dunham/vol2/Cris_L/BeerID/Scripts/OutStd/20191211/
#$ -e /net/dunham/vol2/Cris_L/BeerID/Scripts/ErrStd/20191211/
#$ -N Tree
#$ -l mfree=4G
#$ -pe serial 20
#$ -l h_rt=10:0:0:0

module load modules modules-init modules-gs
module load java/8u25
module load bwa/0.7.15
module load picard/2.9.0
module load samtools/1.9

module load htslib/1.9
module load bcftools/1.9

WORK_DIR=/net/dunham/vol2/Cris_L/BeerID/WorkDirectory_Combined # Give the desired location of intermediate output files.
SCRIPT_DIR=/net/dunham/vol2/Cris_L/BeerID/Scripts
FASTA_DIR=/net/dunham/vol2/Cris_L/BeerID/WorkDirectory/Fasta/AmericanBeer_NoBE051_20200430
REF=/net/dunham/vol2/Cris_L/ReferenceGenome/Reference_2Micron/sacCer3_2micron.fasta

SNPEFF_DIR=/net/dunham/vol2/Cris_L/Software/snpEff
VCF2PHYLIP=/net/dunham/vol2/Cris_L/Software/vcf2phylip
OUT_NAME=AmericanBeer_NoBE051_20200430
GATK_DIR=/net/dunham/vol2/Cris_L/Software/gatk-4.1.3.0

IQTREE=/net/dunham/vol2/Cris_L/Software/iqtree-2.0-rc1-Linux

cd ${WORK_DIR}

sed -e '$s/$/\n/' -s ${FASTA_DIR}/*.cat.fasta > ${WORK_DIR}/${OUT_NAME}.fasta
sed -e '$s/^$//' ${WORK_DIR}/${OUT_NAME}.fasta > ${WORK_DIR}/${OUT_NAME}.fasta.temp
mv ${WORK_DIR}/${OUT_NAME}.fasta.temp ${WORK_DIR}/${OUT_NAME}.fasta
sed -i 's/*/N/g' ${WORK_DIR}/${OUT_NAME}.fasta

rm ${WORK_DIR}/${OUT_NAME}.fasta.temp

${IQTREE}/bin/iqtree \
	-s ${WORK_DIR}/${OUT_NAME}.fasta \
	-m GTR+F+R4 \
	-B 1000 \
	--seqtype DNA \
	-T 20
