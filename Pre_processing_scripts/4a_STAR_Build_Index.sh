#!/bin/bash

#SBATCH --job-name=index_STAR 
#SBATCH --partition=long
#SBATCH --cpus-per-task=1
#SBATCH --mem=42G
#SBATCH --nodes=1  
#SBATCH --output=logs/index_STAR.out
#SBATCH --error=logs/index_STAR.err

#=========================
# User defined parameters: relevant paths
#=========================

# NOTE: You can skip this script if you already have an index built in a previous
# project. You should built it again if you use a different STAR tool version, 
# different read length or different genome assembly

#=========================
# General configuration: paths and singularity images binding
#=========================

# Root directory in the cluster 
ROOTDIR=$(sed -n 6p RNAseq_User_defined_parameters.txt)
# Project working directory. STAR index will be stored there
WKD=$ROOTDIR/$(sed -n 12p RNAseq_User_defined_parameters.txt)
# GTF annotation reference genome
GTF=$ROOTDIR/$(sed -n 19p RNAseq_User_defined_parameters.txt)
# FASTA sequence reference genome
FASTA=$ROOTDIR/$(sed -n 23p RNAseq_User_defined_parameters.txt)
# ReadLength-1 for sjdbOverhand 
LENGTH=$(sed -n 27p RNAseq_User_defined_parameters.txt) 

START=$(date +%s)
# Enable Singularity image to look into the general path (equivalent to -B)
export SINGULARITY_BIND=$ROOTDIR 
# Path to images folder in cluster
IMAGES_PATH=$ROOTDIR"/images"
# Path to databases folder in cluster
DB_PATH=$ROOTDIR"/db_files"

# Folder where to store the STAR index
INDEX=$WKD'/STAR_align/Index'

#=========================
# General configuration: paths and singularity images binding
#=========================

# STAR index SHALL be re-computed in following cases:
#	 .- different reference genome
#	 .- different STAR version
#	 .- your reads (maximum length) considerably differ (i.e. 50bp vs 150bp). In 
#	    case you are interested in splicing junctions, you SHOULD always re-build
# 		the index with the proper length

# Link to STAR manual
# https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

# Specify image/s name to be used (tool-related)
STAR='star_2.7.8a.sif'  #This image inludes START 2.7.8

# Specify any particular tool parameters
# Number of threads                         
T='8' 

#=========================
# Execution: Build Genome Index
#=========================

DATE=$(date +%m-%d-%Y--%T)
echo "Starting building Genome Index: $DATE"
echo ''

# NOTE: Include this parameter if using gff3 instead of gtf annotations
#--sjdbGTFtagExonParentTranscript Parent

singularity exec $IMAGES_PATH/$STAR STAR --runMode genomeGenerate --genomeDir $INDEX --genomeFastaFiles $FASTA --sjdbOverhang $LENGTH --sjdbGTFfile $GTF --runThreadN $T 

#=========================
# End Preprocessing
#=========================

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'Genome Index Built' 
echo "Processing Time: $DIFF seconds"
