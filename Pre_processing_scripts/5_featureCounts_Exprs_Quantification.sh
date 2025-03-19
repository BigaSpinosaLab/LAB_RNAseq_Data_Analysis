#!/bin/bash

#SBATCH --job-name=fCounts_PDOs 
#SBATCH --partition=long
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=2 
#SBATCH --mem=8G
#SBATCH --output=logs/fCounts_quant.out
#SBATCH --error=logs/fCounts_quant.err

#=========================
# User defined parameters: relevant paths
#=========================

# NOTE: Check tool parametrization in case you need to change anything (i.e. library type)

#=========================
# General configuration
#=========================

# Root directory
ROOTDIR=$(sed -n 6p RNAseq_User_defined_parameters.txt)
# Project working directory
WKD=$ROOTDIR/$(sed -n 12p RNAseq_User_defined_parameters.txt)
# Project name
PROJECT=$(sed -n 9p RNAseq_User_defined_parameters.txt)
# GTF annotation reference genome
GTF=$ROOTDIR/$(sed -n 19p RNAseq_User_defined_parameters.txt)

START=$(date +%s)
# Enable Singularity image to look into the general path (equivalent to -B)
export SINGULARITY_BIND=$ROOTDIR 
# Path to images folder in cluster
IMAGES_PATH=$ROOTDIR"/images"
# Path to databases folder in cluster
DB_PATH=$ROOTDIR"/db_files"

# Folder where BAM files are located
DATA=$WKD'/STAR_align/BAM'
# Folder where featureCounts output should be stored
OUT=$WKD'/quantification'

#=========================
# Singularity image and Tool Parametrization
#=========================

# Link to featureCounts tutorial
# http://subread.sourceforge.net/featureCounts.html

# Specify image/s name to be used (tool-related)
FCOUNTS='featureCounts_subRead_v2.0.1.sif'  #This image inludes featureCounts from subRead v2.0.1

# Specify any particular tool parameters
# Number of threads                         
T='6' 

# Stranded library
# 0 is for unstranded, 1 for stranded and 2 for reversely stranded
# Typically = 2, if succesfully assigned reads is extremely low, switch to 1
# If single-end reads => typically unstranded
ST='2'                                          

# Feature type (-t)
ft='exon'   #exon is by default
#ft='gene'   #gene feature to account for intronic reads

# Attribute type to group features based on GTF annotations (gene id by default)
at='gene_id'

# Other interesting parameters:
# -p only applicable for paired-end: fragments to be counted instead of reads.
# -C for NOT counting chimeric reads (although already discarded in STAR, by default)
# -B only count read pairs that have both ends aligned
# -M if multi-mappers have to be also counted

#=========================
# Execution: gene quantification with featureCounts
#=========================

# Listing all BAM files from which you wanna quantify genes
BAM=$(for f in $DATA/*.bam; do printf '%s ' $f; done)

echo ''
echo 'Begin Read Quantification with featureCounts....................... '`date`
echo ''

singularity exec $IMAGES_PATH/$FCOUNTS featureCounts  -T $T -s $ST -t $ft -g $at -a $GTF -pCB -o $OUT/$PROJECT'_exprs_quantification.txt' $BAM

#=========================
# End Preprocessing
#=========================

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Processing Time: $DIFF seconds"