#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Cris_L/PostDocBrewing/
#$ -o /net/dunham/vol2/Cris_L/PostDocBrewing/Scripts/OutStd/20190502/
#$ -e /net/dunham/vol2/Cris_L/PostDocBrewing/Scripts/ErrStd/20190502/
#$ -N Align_PostDoc
#$ -l mfree=4G
#$ -l h_rt=10:0:0:0

## SNP calling and alignment pipeline for PDB data
## Chris Large and Caiti S. Heil
## Uses the recommended SNP calling pipeline from Samtools
## Then filters based on the ancestral sequence

module load modules modules-init modules-gs
module load java/8u25
module load zlib/latest
module load bwa/latest
module load htslib/1.9
module load samtools/latest
module load picard/2.6.0
module load GATK/3.7
module load trimmomatic/0.32
module load perl/latest
module load VCFtools/latest
module load bcftools/latest
module load bedtools/latest

SAMPLE=$1 #sample prefix (ex: Sample-01)
WORKDIR=/net/dunham/vol2/Cris_L/PostDocBrewing/WorkDirectory/Populations #CHANGE
REF=/net/dunham/vol2/Caiti/reference_seq/sacCer3.fasta
SCRIPTS=/net/dunham/vol2/Cris_L/Aaron_Reanalyze/Scripts

cd ${WORKDIR}
mkdir -p ${WORKDIR}/${SAMPLE}
cd ${WORKDIR}/${SAMPLE}
java -Xmx2g -jar $GATK_DIR/GenomeAnalysisTK.jar \
	-R ${REF} \
	-T HaplotypeCaller \
	--output_mode EMIT_ALL_CONFIDENT_SITES \
	-I ${WORKDIR}/${SAMPLE}/${SAMPLE}_comb_R1R2.RG.MD.realign.sort.bam \
	-o ${WORKDIR}/${SAMPLE}/${SAMPLE}_GATK_HaplotypeCaller_ALL.vcf

# Extract SNPs from call set.
java -Xmx2g -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R ${REF} \
	-V ${WORKDIR}/${SAMPLE}/${SAMPLE}_GATK_HaplotypeCaller_ALL.vcf \
	--selectTypeToExclude INDEL \
        --selectTypeToExclude MNP \
        --selectTypeToExclude SYMBOLIC \
        --selectTypeToExclude MIXED \
	-o ${WORKDIR}/${SAMPLE}/${SAMPLE}_GATK_HaplotypeCaller_ALL_SNPS.vcf

# Apply filters to call sets.
java -Xmx2g -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R ${REF} \
	-V ${WORKDIR}/${SAMPLE}/${SAMPLE}_GATK_HaplotypeCaller_ALL_SNPS.vcf \
	--filterExpression \
	"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
	--filterName "GATK_Best_Practices_default" \
	-o ${WORKDIR}/${SAMPLE}/${SAMPLE}_GATK_HaplotypeCaller_ALL_SNPS_Tagged.vcf

# Select only variants that have passed the filters.
java -Xmx2g -jar ${GATK_DIR}/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R ${REF} \
	-V ${WORKDIR}/${SAMPLE}/${SAMPLE}_GATK_HaplotypeCaller_ALL_SNPS_Tagged.vcf \
	-select 'vc.isNotFiltered()' \
	-o ${WORKDIR}/${SAMPLE}/${SAMPLE}_GATK_HaplotypeCaller_ALL_SNPS_Tagged_Filtered.vcf

# Variants to table
java -Xmx2g -jar $GATK_DIR/GenomeAnalysisTK.jar \
        -R $REF \
        -T VariantsToTable \
        -V ${WORKDIR}/${SAMPLE}/${SAMPLE}_GATK_HaplotypeCaller_ALL_SNPS_Tagged_Filtered.vcf \
        -F CHROM -F POS -GF AD \
        -o ${WORKDIR}/${SAMPLE}/${SAMPLE}_GATK_HaplotypeCaller_ALL_SNPS_Tagged_Filtered.table
## Need java 8 jdk for this shenanigans
module load java_jdk/8u91

## Changes the table format into allele frequencies for plotting
# double AlleleFreq = (Double.parseDouble(APArray[0]) / Total); // REF is 1 Alt is 0
java -Xmx2g -jar ${SCRIPTS}/VaraiantTableParse_GATK_BAF_20190206.jar \
	${WORKDIR}/${SAMPLE}/${SAMPLE}_GATK_HaplotypeCaller_ALL_SNPS_Tagged_Filtered.table \
	${SCRIPTS}/ChromSizes_PDB.txt \
	${WORKDIR}/${SAMPLE}/

module load gcc/8.1.0
module load R/3.5.1

## Plot allele frequency REF is 1.0 and Alt is 0.0
Rscript ${SCRIPTS}/PlotAlleleFreq_OneSample_20190207.R \
	${WORKDIR} \
	${SAMPLE} ## Sample name

