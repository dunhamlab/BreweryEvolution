#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Cris_L/BeerID/
#$ -o /net/dunham/vol2/Cris_L/BeerID/Scripts/OutStd/20191022/
#$ -e /net/dunham/vol2/Cris_L/BeerID/Scripts/ErrStd/20191022/
#$ -N AlignBam
#$ -l mfree=16G
#$ -l h_rt=10:0:0:0

module load modules modules-init modules-gs
module load java/8u25
module load bwa/0.7.15
module load GATK/4.1.1.0
module load picard/2.9.0
module load samtools/1.9

## Pulls relevant info from the sample sheet that has been submitted with runQsub_SampleList.py
SAMPLE=$1
LIBRARY_NAME=$2
COMMON_NAME=$3
DATASET=$5
WORK_DIR=/net/dunham/vol2/Cris_L/BeerID/WorkDirectory/${DATASET}_Dataset
FASTQ_DIR=/net/dunham/vol2/Cris_L/BeerID/Fastqs/${DATASET}_Dataset
## Aligns to the S288C reference genome with the 2-Micron plasmid added
REF=/net/dunham/vol2/Cris_L/ReferenceGenome/Reference_2Micron/sacCer3_2micron.fasta

cd ${WORK_DIR}

## Puts all adapter trim info in a joint file in the work directory
mkdir -p AdapterTrimMetrics

## Prints to ErrStd file what the sample name is
(>&2 echo ${SAMPLE})

## Fastq -> Unaligned Bam
java -Xmx8G -jar ${PICARD_DIR}/picard.jar FastqToSam \
	FASTQ=${FASTQ_DIR}/${SAMPLE}_1.fastq.gz \
	FASTQ2=${FASTQ_DIR}/${SAMPLE}_2.fastq.gz \
	OUTPUT=${WORK_DIR}/${SAMPLE}_FastqToSam.bam \
	READ_GROUP_NAME=${SAMPLE} \
	SAMPLE_NAME=${COMMON_NAME} \
	LIBRARY_NAME=${LIBRARY_NAME} \
	PLATFORM=illumina

## Mark adapters in unaligned Bam
java -Xmx8G -jar ${PICARD_DIR}/picard.jar MarkIlluminaAdapters \
	I=${WORK_DIR}/${SAMPLE}_FastqToSam.bam \
	O=${WORK_DIR}/${SAMPLE}_FastqToSam_AdaptersMarked.bam \
	M=${WORK_DIR}/AdapterTrimMetrics/${SAMPLE}_AdapterMetrics.txt \
	ADAPTERS=NEXTERA_V2 \
	TMP_DIR=/net/dunham/vol2/Cris_L/BeerID/WorkDirectory

## SamToFastq to convert unmapped bam into fastq
java -Xmx8G -jar ${PICARD_DIR}/picard.jar SamToFastq \
	I=${WORK_DIR}/${SAMPLE}_FastqToSam_AdaptersMarked.bam \
	FASTQ=${SAMPLE}_AdaptersMarked_Interleaved.fastq \
	CLIPPING_ATTRIBUTE=XT \
	CLIPPING_ACTION=2 \
	INTERLEAVE=true \
	NON_PF=true \
	TMP_DIR=/net/dunham/vol2/Cris_L/BeerID/WorkDirectory

## Map the interleaved reads
bwa mem -M -p ${REF} ${WORK_DIR}/${SAMPLE}_AdaptersMarked_Interleaved.fastq > ${WORK_DIR}/${SAMPLE}_AdaptersMarked_Mapped.bam

## Merge the unmapped bam meta-data of the reads with the aligned bam
java -Xmx8G -jar ${PICARD_DIR}/picard.jar MergeBamAlignment \
	ALIGNED_BAM=${WORK_DIR}/${SAMPLE}_AdaptersMarked_Mapped.bam \
	UNMAPPED_BAM=${WORK_DIR}/${SAMPLE}_FastqToSam.bam \
	OUTPUT=${WORK_DIR}/${SAMPLE}_AdaptersMarked_Mapped_Merged.bam \
	R=${REF} CREATE_INDEX=true ADD_MATE_CIGAR=true \
	CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
	INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
	PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
	TMP_DIR=/net/dunham/vol2/Cris_L/BeerID/WorkDirectory

## Removes intermediate files
rm ${WORK_DIR}/${SAMPLE}_FastqToSam.bam
rm ${WORK_DIR}/${SAMPLE}_FastqToSam_AdaptersMarked.bam
rm ${WORK_DIR}/${SAMPLE}_AdaptersMarked_Interleaved.fastq
rm ${WORK_DIR}/${SAMPLE}_AdaptersMarked_Mapped.bam

