#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Cris_L/BeerID/
#$ -o /net/dunham/vol2/Cris_L/BeerID/Scripts/OutStd/20191217/
#$ -e /net/dunham/vol2/Cris_L/BeerID/Scripts/ErrStd/20191217/
#$ -N Combine
#$ -l mfree=180G
#$ -l h_rt=20:0:0:0

module load modules modules-init modules-gs
module load java/8u25
module load bwa/0.7.15
module load picard/2.9.0
module load samtools/1.9
module load tabix/0.2.6

WORK_DIR=/net/dunham/vol2/Cris_L/BeerID/WorkDirectory_Combined # Give the desired location of intermediate output files.
SCRIPT_DIR=/net/dunham/vol2/Cris_L/BeerID/Scripts
REF=/net/dunham/vol2/Cris_L/ReferenceGenome/Reference_2Micron/sacCer3_2micron.fasta
FILE_LIST=/net/dunham/vol2/Cris_L/BeerID/Scripts/$1
SNPEFF_DIR=/net/dunham/vol2/Cris_L/Software/snpEff
VCF2PHYLIP=/net/dunham/vol2/Cris_L/Software/vcf2phylip
IQTREE=/net/dunham/vol2/Cris_L/Software/iqtree-2.0-rc1-Linux
OUT_NAME=$2
GATK_DIR=/net/dunham/vol2/Cris_L/Software/gatk-4.1.3.0
GATK_LIST=${SCRIPT_DIR}/GATK_Style_Interval.list

${GATK_DIR}/gatk GenomicsDBImport \
	--java-options "-Xmx165g -Xms165g" \
	-R ${REF} \
	--variant ${FILE_LIST} \
	-L ${GATK_LIST} \
	--genomicsdb-workspace-path ${WORK_DIR}/${OUT_NAME} \
	--tmp-dir=${WORK_DIR}/

echo ImportFinished

${GATK_DIR}/gatk GenotypeGVCFs \
	--java-options "-Xmx165g -Xms165g" \
	-R ${REF} \
	--variant gendb://${WORK_DIR}/${OUT_NAME}  \
	-O ${WORK_DIR}/${OUT_NAME}.vcf

${GATK_DIR}/gatk SelectVariants \
	-R ${REF} \
	--variant ${WORK_DIR}/${OUT_NAME}.vcf \
	--select-type-to-include SNP \
	-O ${WORK_DIR}/${OUT_NAME}.snp.vcf

${GATK_DIR}/gatk SelectVariants \
        -R ${REF} \
        --variant ${WORK_DIR}/${OUT_NAME}.vcf \
        --select-type-to-include INDEL \
        -O ${WORK_DIR}/${OUT_NAME}.indel.vcf

${GATK_DIR}/gatk VariantFiltration \
	-V ${WORK_DIR}/${OUT_NAME}.indel.vcf \
	--filter-name "QD_Filter" \
	--filter-expression "QD < 2.0" \
	--filter-name "FS_Filter" \
	--filter-expression "FS > 200.0" \
	--filter-name "SOR_Filter" \
	--filter-expression "SOR > 10.0" \
	--filter-name "ReadPosRankSum_Filter" \
	--filter-expression "ReadPosRankSum < -20.0" \
        --genotype-filter-expression "DP<10" \
        --genotype-filter-name "DP_Filter" \
	--verbosity ERROR \
	-O ${WORK_DIR}/${OUT_NAME}.indel.prefilt.vcf

${GATK_DIR}/gatk VariantFiltration \
        --java-options "-Xmx165g -Xms165g" \
	-R ${REF} \
	-V ${WORK_DIR}/${OUT_NAME}.vcf \
	--filter-name "QD_Filter" \
	--filter-expression "QD < 2.0" \
        --filter-name "FS_Filter" \
        --filter-expression "FS > 60.0" \
        --filter-name "SOR_Filter" \
        --filter-expression "SOR > 3.0" \
        --filter-name "MQ_Filter" \
        --filter-expression "MQ < 40.0" \
        --filter-name "MQRankSum_Filter" \
        --filter-expression "MQRankSum < -12.5" \
        --filter-name "ReadPosRankSum_Filter" \
        --filter-expression "ReadPosRankSum < -8.0" \
	--genotype-filter-expression "DP<10" \
	--genotype-filter-name "DP_Filter" \
	--verbosity ERROR \
	-O ${WORK_DIR}/${OUT_NAME}.snp.prefilt.vcf

${GATK_DIR}/gatk SelectVariants \
        --java-options "-Xmx165g -Xms165g" \
        -R ${REF} \
        --variant ${WORK_DIR}/${OUT_NAME}.snp.prefilt.vcf \
        --exclude-filtered true \
	--set-filtered-gt-to-nocall \
        --select-type-to-include SNP \
	--remove-unused-alternates true \
        -O ${WORK_DIR}/${OUT_NAME}.snp.filt.vcf

${GATK_DIR}/gatk SelectVariants \
        --java-options "-Xmx165g -Xms165g" \
        -R ${REF} \
        --variant ${WORK_DIR}/${OUT_NAME}.indel.prefilt.vcf \
        --exclude-filtered true \
        --set-filtered-gt-to-nocall \
        -O ${WORK_DIR}/${OUT_NAME}.indel.filt.vcf

${GATK_DIR}/gatk MergeVcfs \
	-I ${WORK_DIR}/${OUT_NAME}.snp.filt.vcf \
	-I ${WORK_DIR}/${OUT_NAME}.indel.filt.vcf \
	-O ${WORK_DIR}/${OUT_NAME}.comb.filt.vcf

bgzip -c ${WORK_DIR}/${OUT_NAME}.comb.filt.vcf > ${WORK_DIR}/${OUT_NAME}.comb.filt.vcf.gz
bgzip -c ${WORK_DIR}/${OUT_NAME}.snp.filt.vcf > ${WORK_DIR}/${OUT_NAME}.snp.filt.vcf.gz

tabix -f ${WORK_DIR}/${OUT_NAME}.comb.filt.vcf.gz
tabix -f ${WORK_DIR}/${OUT_NAME}.snp.filt.vcf.gz

rm ${WORK_DIR}/${OUT_NAME}.snp.prefilt.vcf
rm ${WORK_DIR}/${OUT_NAME}.indel.prefilt.vcf
