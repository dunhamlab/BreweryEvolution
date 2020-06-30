#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Cris_L/BeerID/
#$ -o /net/dunham/vol2/Cris_L/BeerID/Scripts/OutStd/20191024/
#$ -e /net/dunham/vol2/Cris_L/BeerID/Scripts/ErrStd/20191024/
#$ -N Variant
#$ -pe serial 4
#$ -l mfree=16G
#$ -l h_rt=15:0:0:0

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
PLOIDY=$4
DATASET=$5
WORK_DIR=/net/dunham/vol2/Cris_L/BeerID/WorkDirectory/${DATASET}_Dataset
FASTQ_DIR=/net/dunham/vol2/Cris_L/BeerID/Fastqs/${DATASET}_Dataset
## Aligns to the S288C reference genome with the 2-Micron plasmid added
REF=/net/dunham/vol2/Cris_L/ReferenceGenome/Reference_2Micron/sacCer3_2micron.fasta
SAMPLE_SHEET=/net/dunham/vol2/Cris_L/BeerID/Scripts/SampleSheet/SampleSheet_4251.txt

cd ${WORK_DIR}

## Test whether there is avaliable ploidy information from the sample sheet.
# If not, then set the ploidy designation to 2
if test ${PLOIDY} = '0'; then 
	PLOIDY=4
else
	echo ${PLOIDY}
fi
##
#### Casts read_group files into array from the sample sheet using the sample name
### In effect, it will dynamically be able to grab all of the read_group files for a given sample name
readarray -t bam_array < <(grep $LIBRARY_NAME $SAMPLE_SHEET | awk -v var="${WORK_DIR}/" '{print "INPUT="var$1"_AdaptersMarked_Mapped_Merged.bam"}')
#
#### Use bam_array to merge and MarkDuplicates for a given sample
java -Xmx8G -jar ${PICARD_DIR}/picard.jar MarkDuplicates \
	$(echo ${bam_array[@]}) \
	OUTPUT=${WORK_DIR}/${LIBRARY_NAME}_AdaptersMarked_MD.bam \
	CREATE_INDEX=true \
	METRICS_FILE=${WORK_DIR}/${LIBRARY_NAME}_MD_Metrics.txt \
	TMP_DIR=${WORK_DIR}

java -Xmx8G -jar ${PICARD_DIR}/picard.jar SortSam \
	INPUT=${WORK_DIR}/${LIBRARY_NAME}_AdaptersMarked_MD.bam \
	OUTPUT=${WORK_DIR}/${LIBRARY_NAME}_AdaptersMarked_MD_sort.bam \
	SORT_ORDER=coordinate \
	CREATE_INDEX=true \
	TMP_DIR=${WORK_DIR}

rm ${WORK_DIR}/${LIBRARY_NAME}_AdaptersMarked_MD.bam

qsub /net/dunham/vol2/Cris_L/BeerID/Scripts/runFlagstat_DepthOfCoverage_20190909.sh $1 $2 $3 $4 $5

${GATK_DIR}/gatk HaplotypeCaller \
    --java-options "-Xmx8g -Xms8g" \
	--native-pair-hmm-threads 4 \
	-R ${REF} \
	-I ${WORK_DIR}/${LIBRARY_NAME}_AdaptersMarked_MD_sort.bam \
	-ERC GVCF \
	-ploidy ${PLOIDY} \
	-O ${WORK_DIR}/${LIBRARY_NAME}.g.vcf

