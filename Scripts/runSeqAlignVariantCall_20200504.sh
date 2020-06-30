#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Cris_L/PostDocBrewing/
#$ -o /net/dunham/vol2/Cris_L/PostDocBrewing/Scripts/OutStd/20200309/
#$ -e /net/dunham/vol2/Cris_L/PostDocBrewing/Scripts/ErrStd/20200309/
#$ -N Align_PostDoc
#$ -l mfree=8G
#$ -l h_rt=36:0:0

## SNP calling and alignment pipeline for PDB data
## Chris Large and Caiti S. Heil
## Uses the recommended SNP calling pipeline from Samtools
## Then filters based on the ancestral sequence

module load modules modules-init modules-gs
module load java/8u25
module load zlib/latest
module load bwa/0.7.15
module load htslib/1.9
module load samtools/1.9
module load picard/2.6.0
module load GATK/3.7
module load perl/latest
module load VCFtools/latest
module load bcftools/1.9
module load bedtools/latest
module load freebayes/latest
module load lofreq/2.1.2

SAMPLE=$1 #sample prefix (ex: Sample-01)
WORKDIR=/net/dunham/vol2/Cris_L/PostDocBrewing/WorkDirectory/Populations #CHANGE
SEQDIR=/net/dunham/vol2/Cris_L/PostDocBrewing/Fastq #CHANGE
SEQID=20200309
REF=/net/dunham/vol2/Cris_L/ReferenceGenome/sacCer3_PD68_mtDNA.fasta
ANNOTATE=/net/dunham/vol2/Cris_L/ReferenceGenome/AnnotationReference
SCRIPTS=/net/dunham/vol2/Cris_L/Aaron_Reanalyze/Scripts
ANC=$2
ANCBAM=/net/dunham/vol2/Cris_L/PostDocBrewing/WorkDirectory/Populations/${ANC}/${ANC}_comb_R1R2.RG.MD.realign.sort.bam
VCFDIR=/net/dunham/vol2/Cris_L/PostDocBrewing/WorkDirectory/Populations/${ANC}

cd ${WORKDIR}

mkdir -p ${WORKDIR}/${SAMPLE}_mtDNA
cd ${WORKDIR}/${SAMPLE}_mtDNA

###Align with bwa
#(>&2 echo ***BWA - mem -R***)
#bwa mem -R '@RG\tID:'${SEQID}'\tSM:'${SAMPLE}'\tLB:1' ${REF} ${SEQDIR}/${SAMPLE}_*R1*.fastq.gz ${SEQDIR}/${SAMPLE}_*R2*.fastq.gz > ${SAMPLE}_R1R2.sam
#
#(>&2 echo ***Samtools - View***)
#samtools view -bS ${SAMPLE}_R1R2.sam -o ${SAMPLE}_R1R2.bam
#
#(>&2 echo ***Samtools - Sort***)
#samtools sort ${SAMPLE}_R1R2.bam -o ${SAMPLE}_R1R2_sort.bam
#
#(>&2 echo ***Samtools - Index***)
#samtools index ${SAMPLE}_R1R2_sort.bam
#
#(>&2 echo ***Samtools - Flagstat***)
#samtools flagstat ${WORKDIR}/${SAMPLE}/${SAMPLE}_R1R2_sort.bam
#
##Remove intermediate files
#rm ${SAMPLE}_R1R2.sam
#rm ${SAMPLE}_R1R2.bam
#
#cd ${WORKDIR}/${SAMPLE}
#
#mkdir -p dup_metrics
#
#(>&2 echo ***Picard - MarkDuplicates***)
#java -Xmx2g -jar $PICARD_DIR/picard.jar MarkDuplicates \
#        INPUT=${SAMPLE}_R1R2_sort.bam \
#        OUTPUT=${SAMPLE}_comb_R1R2.MD.bam \
#        METRICS_FILE=dup_metrics/${SAMPLE}_comb_R1R2.sort_dup_metrics \
#        REMOVE_DUPLICATES=true \
#        VALIDATION_STRINGENCY=LENIENT
#
#(>&2 echo ***Picard - AddOrReplaceReadGroups***)
##Add or replace read groups needs to happen before GATK
#java -Xmx2g -jar $PICARD_DIR/picard.jar AddOrReplaceReadGroups \
#        I=${SAMPLE}_comb_R1R2.MD.bam \
#        O=${SAMPLE}_comb_R1R2.RG.MD.bam \
#        RGID=${SEQID} \
#        RGLB=1 \
#        RGPU=1 \
#        RGPL=illumina \
#        RGSM=${SAMPLE} \
#        VALIDATION_STRINGENCY=LENIENT
#
#(>&2 echo ***Samtools - Sort and Index***)
#samtools sort ${SAMPLE}_comb_R1R2.RG.MD.bam \
#        -o ${SAMPLE}_comb_R1R2.RG.MD.sort.bam
#samtools index ${SAMPLE}_comb_R1R2.RG.MD.sort.bam
#
#(>&2 echo ***GATK - RealingerTargetCreator***)
##GATK Realinger
#java -Xmx2g -jar $GATK_DIR/GenomeAnalysisTK.jar \
#        -T RealignerTargetCreator \
#        -R ${REF} \
#        -I ${SAMPLE}_comb_R1R2.RG.MD.sort.bam \
#        -o ${SAMPLE}_comb_R1R2.bam.intervals
#java -Xmx2g -jar $GATK_DIR/GenomeAnalysisTK.jar \
#        -T IndelRealigner \
#        -R ${REF} \
#        -I ${SAMPLE}_comb_R1R2.RG.MD.sort.bam \
#        -targetIntervals ${SAMPLE}_comb_R1R2.bam.intervals \
#        -o ${SAMPLE}_comb_R1R2.RG.MD.realign.bam
#
#samtools sort ${SAMPLE}_comb_R1R2.RG.MD.realign.bam \
#        -o ${SAMPLE}_comb_R1R2.RG.MD.realign.sort.bam
#samtools index ${SAMPLE}_comb_R1R2.RG.MD.realign.sort.bam
#
cd ${WORKDIR}/${SAMPLE}/

###pileup
###note: when running for ancestor or control, use varFilter -d1 
##(>&2 echo ***BCFtools - Pileup***)
##bcftools mpileup --ignore-RG -Ou -ABf ${REF} ${SAMPLE}_comb_R1R2.RG.MD.realign.sort.bam | bcftools call -vmO v -o ${SAMPLE}_samtools_AB.vcf
#
#freebayes -f ${REF} \
#        --pooled-discrete --pooled-continuous --report-genotype-likelihood-max --allele-balance-priors-off --min-alternate-fraction 0.1 \
#        ${SAMPLE}_comb_R1R2.RG.MD.realign.sort.bam > ${SAMPLE}_freebayes_BCBio.vcf
#
## If ancestor, comment out below section
#(>&2 echo ***LoFreq - Somatic***)
lofreq somatic -n ${ANCBAM} -t ${WORKDIR}/${SAMPLE}/${SAMPLE}_comb_R1R2.RG.MD.realign.sort.bam -f ${REF} \
       -o ${SAMPLE}_lofreq_

bgzip -d ${SAMPLE}_lofreq_somatic_final.snvs.vcf.gz
bgzip -d ${SAMPLE}_lofreq_tumor_relaxed.vcf.gz
bgzip -d ${SAMPLE}_lofreq_normal_relaxed.vcf.gz

#(>2 echo ***Bedtools - Intersect***)
#bedtools intersect -v -header \
#        -a ${WORKDIR}/${SAMPLE}/${SAMPLE}_samtools_AB.vcf \
#        -b ${VCFDIR}/${ANC}_samtools_AB.vcf \
#        > ${WORKDIR}/${SAMPLE}/${SAMPLE}_samtools_AB_AncFiltered.vcf

bedtools intersect -v -header \
        -a ${WORKDIR}/${SAMPLE}/${SAMPLE}_freebayes_BCBio.vcf \
        -b ${VCFDIR}/${ANC}_freebayes_BCBio.vcf \
        > ${WORKDIR}/${SAMPLE}/${SAMPLE}_freebayes_BCBio_AncFiltered.vcf

bedtools intersect -v -header \
        -a ${WORKDIR}/${SAMPLE}/${SAMPLE}_lofreq_tumor_relaxed.vcf \
        -b ${WORKDIR}/${SAMPLE}/${SAMPLE}_lofreq_normal_relaxed.vcf \
        > ${WORKDIR}/${SAMPLE}/${SAMPLE}_lofreq_tumor_relaxed_AncFiltered.vcf

##(>2 echo ***BCFtools - Filter***)
#bcftools filter -O v -o ${SAMPLE}_samtools_AB_AncFiltered.filt.vcf \
#        -i 'MQ>30 & QUAL>50 & DP>40 & (DP4[2]+DP4[3])>4 & (DP4[0]+DP4[2])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.01 & (DP4[1]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.01' \
#        ${SAMPLE}_samtools_AB_AncFiltered.vcf

bcftools filter -O v -o ${SAMPLE}_freebayes_BCBio_AncFiltered.filt.vcf \
        -i 'MQM>30 & QUAL>20 & INFO/DP>40 & (SAF+SAR)>4 & (SRF+SAF)/(INFO/DP)>0.01 & (SRR+SAR)/(INFO/DP)>0.01' \
        ${SAMPLE}_freebayes_BCBio_AncFiltered.vcf

bcftools filter -O v -o ${SAMPLE}_lofreq_tumor_relaxed_AncFiltered.filt.vcf \
       -i 'QUAL>20 & DP>20 & (DP4[2]+DP4[3])>4 & (DP4[0]+DP4[2])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.01 & (DP4[1]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.01' \
       ${SAMPLE}_lofreq_tumor_relaxed_AncFiltered.vcf

#bedtools intersect -v -header \
#        -a ${SAMPLE}_samtools_AB_AncFiltered.filt.vcf \
#        -b ${SAMPLE}_freebayes_BCBio_AncFiltered.filt.vcf ${WORKDIR}/${SAMPLE}/${SAMPLE}_lofreq_tumor_relaxed_AncFiltered.filt.vcf \
#        > ${SAMPLE}_samtools_AB_AncFiltered.filt.noOverlap.vcf

bedtools intersect -v -header \
        -a ${SAMPLE}_freebayes_BCBio_AncFiltered.filt.vcf \
        -b ${WORKDIR}/${SAMPLE}/${SAMPLE}_lofreq_tumor_relaxed_AncFiltered.vcf \
        > ${SAMPLE}_freebayes_BCBio_AncFiltered.filt.noOverlap.vcf

(>2 echo ***Annotate***)
python ${SCRIPTS}/yeast_annotation_chris_edits_20170925.py \
        -f ${WORKDIR}/${SAMPLE}/${SAMPLE}_freebayes_BCBio_AncFiltered.filt.noOverlap.vcf \
        -s ${ANNOTATE}/orf_coding_all_R64-1-1_20110203.fasta \
        -n ${ANNOTATE}/saccharomyces_cerevisiae_R64-1-1_20110208.gff.filtered \
        -g ${ANNOTATE}/S288C_reference_sequence_R64-1-1_20110203.fsa

#python ${SCRIPTS}/yeast_annotation_chris_edits_20170925.py \
#        -f ${WORKDIR}/${SAMPLE}/${SAMPLE}_samtools_AB_AncFiltered.filt.noOverlap.vcf \
#        -s ${ANNOTATE}/orf_coding_all_R64-1-1_20110203.fasta \
#        -n ${ANNOTATE}/saccharomyces_cerevisiae_R64-1-1_20110208.gff.filtered \
#        -g ${ANNOTATE}/S288C_reference_sequence_R64-1-1_20110203.fsa

python ${SCRIPTS}/yeast_annotation_chris_edits_20170925.py \
        -f ${WORKDIR}/${SAMPLE}/${SAMPLE}_lofreq_tumor_relaxed_AncFiltered.filt.vcf \
        -s ${ANNOTATE}/orf_coding_all_R64-1-1_20110203.fasta \
        -n ${ANNOTATE}/saccharomyces_cerevisiae_R64-1-1_20110208.gff.filtered \
        -g ${ANNOTATE}/S288C_reference_sequence_R64-1-1_20110203.fsa

awk '$8 = $8 FS "NA NA"' ${SAMPLE}_lofreq_tumor_relaxed_AncFiltered.filt.vcf.annotated > ${SAMPLE}_lofreq_tumor_relaxed_AncFiltered.filt.twoCol.vcf.annotated

# Remove
rm ${SAMPLE}_R1R2.bam.intervals
rm ${SAMPLE}_R1R2.MD.bam
rm ${SAMPLE}_R1R2.RG.MD.bam
rm ${SAMPLE}_R1R2.RG.MD.realign.bam
rm ${SAMPLE}_R1R2.RG.MD.realign.bai
rm ${SAMPLE}_R1R2_sort.bam
rm ${SAMPLE}_R1R2_sort.bam.bai
