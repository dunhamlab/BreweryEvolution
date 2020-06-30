#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Cris_L/PostDocBrewing/WorkDirectory
#$ -e /net/dunham/vol2/Cris_L/PostDocBrewing/WorkDirectory/BeerStrains_FromOmar/StdErr
#$ -o /net/dunham/vol2/Cris_L/PostDocBrewing/WorkDirectory/BeerStrains_FromOmar/StdOut
#$ -N PolishAssembly
#$ -pe serial 8
#$ -l mfree=4G
#$ -l h_rt=10:0:0:0

module load modules modules-init modules-gs
module load R/3.6.1
module load minimap2/2.17
module load htslib/1.9
module load samtools/1.9
module load bwa/0.7.17
module load racon/1.4.13
module load pilon/1.23

module load python/3.6.4
module load numpy/1.17.3
module load tensorflow/1.15.2
module load keras/2.3.1
module load pysam/0.15.2
module load cffi/1.13.2
module load h5py/2.9.0
module load intervaltree/3.0.2
module load ont-fast5-api/3.0.2
module load parasail/1.2
module load mappy/2.17
module load medaka/0.11.5

PILON_DIR=/net/gs/vol3/software/modules-sw/pilon/1.23/Linux/RHEL6/x86_64
REF=/net/dunham/vol2/Cris_L/ReferenceGenome/sacCer3.fasta
SAMPLE=$1

DIR=/net/dunham/vol2/Cris_L/PostDocBrewing/WorkDirectory/BeerStrains_FromOmar

mkdir ${DIR}/PolishedAssembly/${SAMPLE}
WORKDIR=${DIR}/PolishedAssembly/${SAMPLE} #CHANGE
SEQDIR=${DIR}/minIONReads #CHANGE

NPROC=$(nproc)
DRAFT=${DIR}/deNovo/${SAMPLE}.dmo.lay.utg.fasta
OUTDIR=${DIR}/PolishedAssembly/${SAMPLE}

#(>&2 echo *** Minimap2 ***)
minimap2 -t ${NPROC} -x map-ont \
	${DRAFT} \
	${SEQDIR}/porechop_${SAMPLE}.fastq.gz > ${WORKDIR}/${SAMPLE}_ont.paf

##(>&2 echo *** Racon ***)
racon \
	-m 8 -x -6 -g -8 -w 500 \
	${SEQDIR}/porechop_${SAMPLE}.fastq.gz \
	${WORKDIR}/${SAMPLE}_ont.paf \
	${DRAFT} > ${OUTDIR}/${SAMPLE}.racon.fasta

bwa index ${OUTDIR}/${SAMPLE}.racon.fasta

##(>&2 echo *** Medaka ***)
medaka_consensus \
	-i ${SEQDIR}/porechop_${SAMPLE}.fastq.gz \
	-d ${OUTDIR}/${SAMPLE}.racon.fasta \
	-o ${OUTDIR} \
	-t ${NPROC} \
	-m r941_min_high_g344

#(>&2 echo *** Bwa - mem ***)
mkdir ${OUTDIR}/pilon

bwa index ${OUTDIR}/consensus.fasta

bwa mem -t ${NPROC} \
	${OUTDIR}/consensus.fasta \
	${DIR}/illuminaReads/${SAMPLE}_R1.fastq.gz \
	${DIR}/illuminaReads/${SAMPLE}_R2.fastq.gz \
	| samtools view - -Sb \
	| samtools sort - -@8 -o ${OUTDIR}/pilon/${SAMPLE}_pilon1.sorted.bam

samtools index ${OUTDIR}/pilon/${SAMPLE}_pilon1.sorted.bam

#(>&2 echo *** Pilon ***)
/bin/env java -Xmx16G -jar ${PILON_DIR}/pilon-1.23.jar \
	--genome ${OUTDIR}/consensus.fasta \
	--fix all \
	--changes \
	--frags ${OUTDIR}/pilon/${SAMPLE}_pilon1.sorted.bam \
	--threads ${NPROC} \
	--output ${OUTDIR}/pilon/${SAMPLE}_pilon1

bwa index ${OUTDIR}/pilon/${SAMPLE}_pilon1.fasta

#(>&2 echo *** Bwa - mem 2 ***)
bwa mem -t ${NPROC} \
        ${OUTDIR}/pilon/${SAMPLE}_pilon1.fasta \
        ${DIR}/illuminaReads/${SAMPLE}_R1.fastq.gz \
        ${DIR}/illuminaReads/${SAMPLE}_R2.fastq.gz \
        | samtools view - -Sb \
        | samtools sort - -@8 -o ${OUTDIR}/pilon/${SAMPLE}_pilon2.sorted.bam

samtools index ${OUTDIR}/pilon/${SAMPLE}_pilon2.sorted.bam

#(>&2 echo *** Pilon 2 ***)
/bin/env java -Xmx16G -jar ${PILON_DIR}/pilon-1.23.jar \
        --genome ${OUTDIR}/pilon/${SAMPLE}_pilon1.fasta \
        --fix all \
        --changes \
        --frags ${OUTDIR}/pilon/${SAMPLE}_pilon2.sorted.bam \
        --threads ${NPROC} \
        --output ${OUTDIR}/pilon/${SAMPLE}_pilon2

bwa index ${OUTDIR}/pilon/${SAMPLE}_pilon2.fasta

#(>&2 echo *** Bwa - mem 3 ***)
bwa mem -t ${NPROC} \
        ${OUTDIR}/pilon/${SAMPLE}_pilon2.fasta \
        ${DIR}/illuminaReads/${SAMPLE}_R1.fastq.gz \
        ${DIR}/illuminaReads/${SAMPLE}_R2.fastq.gz \
        | samtools view - -Sb \
        | samtools sort - -@8 -o ${OUTDIR}/pilon/${SAMPLE}_pilon3.sorted.bam

samtools index ${OUTDIR}/pilon/${SAMPLE}_pilon3.sorted.bam

#(>&2 echo *** Pilon 3 ***)
/bin/env java -Xmx16G -jar ${PILON_DIR}/pilon-1.23.jar \
        --genome ${OUTDIR}/pilon/${SAMPLE}_pilon2.fasta \
        --fix all \
        --changes \
        --frags ${OUTDIR}/pilon/${SAMPLE}_pilon3.sorted.bam \
        --threads ${NPROC} \
        --output ${OUTDIR}/pilon/${SAMPLE}_pilon3


