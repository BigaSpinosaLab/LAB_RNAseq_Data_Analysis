#!/bin/bash

#SBATCH --job-name=PDOFetal_STAR
#SBATCH --partition=long
#SBATCH --cpus-per-task=3 
#SBATCH --mem=12G
#SBATCH --nodes=1  
#SBATCH --output=logs/STAR_alignment.out
#SBATCH --error=logs/STAR_alignment.err
#SBATCH --array=1-21%4

# REMARK!!!! Adapt the number of array tasks to the number of samples i.e. if you have
# 21 samples to trim (42 fastq files) you need to specify 21 as indicated
# NOTE: %4 means that only two tasks will be sent to cluster execution simultaneously

#=========================
# User defined parameters: relevant paths
#=========================

# STAR alignment requires a genome INDEX (previously computed) - see STAR_Build_Index script.
# STAR default parameters are optimized for mammalian genomes, however specific parameters
# can be tuned. Check below section (Tool Parametrization) as an example.

# This script is assuming you are dealing with paired-end data, if this
# is not the case, you have to adapt the trimming tool execution 
# (Go to > Execution section). Additionally, check fastq extensions match with your data

#=========================
# General configuration: paths and singularity images binding
#=========================

# Root directory
ROOTDIR=$(sed -n 6p RNAseq_User_defined_parameters.txt)
# Project working directory
WKD=$ROOTDIR/$(sed -n 12p RNAseq_User_defined_parameters.txt)

START=$(date +%s)
# Enable Singularity image to look into the general path (equivalent to -B)
export SINGULARITY_BIND=$ROOTDIR 
# Path to images folder in cluster
IMAGES_PATH=$ROOTDIR"/images"
# Path to databases folder in cluster
DB_PATH=$ROOTDIR"/db_files"

# Folder where data to be aligned is located (expected 'trimmed_data' or 'raw_data')
DATA=$WKD'/trimmed_data'
# Folder where index is stored
INDEX=$ROOTDIR/$(sed -n 15p User_defined_parameters.txt)
# Folder where STAR alignment results will be stored
OUT=$WKD'/STAR_align/Other_results'
# Folder where BAM files will be finally stored
OUTBAM=$WKD'/STAR_align/BAM'

#=========================
# Singularity image and Tool Parametrization
#=========================

# Link to STAR manual
# https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

# Specify image/s name to be used (tool-related)
STAR='star_2.7.8a.sif'  #This image inludes STAR 2.7.8
SAMTOOLS='samtools_v1.15.sif' # This image includes SAMTOOLS 1.15  

# Specify any particular tool parameters
# Number of threads                         
T='8' 

# Number of allowed mismatches 
MISMATCH='0' # Alignment with 0 mismatch allowed

# Number of multimappers allowed (1==unique mapping reads)
MULTIMAP="1"

# Other options that may be of interest
# Redefined the 'short reads' consideration (by default is 66%)
# Included in the third round of alignment
# SHORT='0.25'

  # If you use it, you have to include next parameters in STAR execution
#--outFilterScoreMinOverLread $SHORT  \
#--outFilterMatchNminOverLread $SHORT  \


#--outReadsUnmapped Fastx \  # Include in STAR execution if you want to store unmapped reads

#=========================
# Execution: STAR alignment
#=========================

# Command file preparation: to execute batch array

for FILENAME in $DATA/*_R1_001.fastq.gz
do
    NAME=${FILENAME%_R1_001.fastq.gz}
    SAMPLE=$(basename $NAME)

    # Forward and Reverse Reads for that sample
	  READ1=$NAME'_R1_001.fastq.gz'
	  READ2=$NAME'_R2_001.fastq.gz'

    # Construct the full execution command for STAR alignment
    echo singularity exec $IMAGES_PATH/$STAR STAR --runThreadN $T \
                                                --genomeDir $INDEX \
                                                --genomeLoad LoadAndKeep \
						                                    --limitBAMsortRAM 12000000000 \
                                                --readFilesIn $READ1 $READ2 \
                                                --readFilesCommand zcat \
                                                --outSAMtype BAM SortedByCoordinate \
                                                --outFilterMismatchNmax $MISMATCH \
                                                --outFilterMultimapNmax $MULTIMAP \
                                                --quantMode GeneCounts  \
                                                --outFileNamePrefix $OUT/$SAMPLE'_'

done > $WKD'/scripts/cmds/STAR_alignment.cmd'

# Execute commant in batch array

echo "-----------------------------------------"
echo "Starting Alignment to reference genome: Loading genome index"
echo "-----------------------------------------"

singularity exec $IMAGES_PATH/$STAR STAR --genomeDir $INDEX

## --genomeLoad LoadAndExit

DATE=$(date +%m-%d-%Y--%T)
echo "  Samples alignment in array mode: $DATE"
echo " "

SEEDFILE=$WKD'/scripts/cmds/STAR_alignment.cmd'
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
eval $SEED

DATE=$(date +%m-%d-%Y--%T)
echo "  All samples aligned: $DATE"


# Clean up: remove index, move BAM files and create respective bam index

echo "  Remove index from memory"
echo " "
singularity exec $IMAGES_PATH/$STAR STAR --genomeDir $INDEX --genomeLoad Remove

echo " Move BAM files to specific folder"
echo " "

for BAM in $OUT/*.bam
do
	mv $BAM $OUTBAM
done

echo " Create BAM index"
for BAM in $OUTBAM/*.bam
do
	singularity exec $IMAGES_PATH/$SAMTOOLS samtools index $BAM
done

#=========================
# End Preprocessing
#=========================

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'STAR alignment completed'
echo "Processing Time: $DIFF seconds"
