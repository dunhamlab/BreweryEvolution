#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Cris_L/PostDocBrewing/WorkDirectory
#$ -e /net/dunham/vol2/Cris_L/PostDocBrewing/WorkDirectory/BeerStrains_FromOmar/StdErr
#$ -o /net/dunham/vol2/Cris_L/PostDocBrewing/WorkDirectory/BeerStrains_FromOmar/StdOut
#$ -N AlignReads
#$ -l mfree=16G

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
module load minimap2/2.17

SAMPLE=$1 #sample prefix (ex: Sample-01)
DIR=/net/dunham/vol2/Cris_L/PostDocBrewing/WorkDirectory/BeerStrains_FromOmar
WORKDIR=${DIR}/Alignments/$1 #CHANGE
SEQDIR=${DIR}/PolishedAssembly/PDB1_c1/${1}.fasta #CHANGE
SEQID=PostDocBrewing_20200210
REF=/net/dunham/vol2/Cris_L/ReferenceGenome/sacCer3.fasta
REF_ANN=/net/dunham/vol2/Cris_L/ReferenceGenome/S288C_reference_genome_R64-2-1_20150113/orf_coding_all_R64-2-1_20150113.fasta

cd ${WORKDIR}

minimap2 -Pa ${SEQDIR} ${REF_ANN} > ${SAMPLE}.RefAnnotations.sam

(>&2 echo ***Samtools - View***)
samtools view -bS ${SAMPLE}.RefAnnotations.sam -o ${SAMPLE}.RefAnnotations.bam

(>&2 echo ***Samtools - Sort***)
samtools sort ${SAMPLE}.RefAnnotations.bam -o ${SAMPLE}.RefAnnotations.sort.bam

(>&2 echo ***Samtools - Index***)
samtools index ${SAMPLE}.RefAnnotations.sort.bam

rm ${SAMPLE}.RefAnnotations.sam
rm ${SAMPLE}.RefAnnotations.bam



