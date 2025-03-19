#!/bin/bash

################################################################################
##       Create Project directory tree for an standard RNA-seq analysis with STAR aligner
################################################################################

# Project directory tree consists of (inside main project folder):
#   
#    .- raw_data: Includes raw data (FASTQ files). Subdirs: FASTQC and MULTIQC
#    .- trimmed_data: Includes trimmed data (from i.e. trimgalore). Subdirs: FASTQC and MULTIQC
#    .- STAR_align: Includes output from STAR alignment. Subdirs:  Index, BAM, Other_results
#    .- quantification: Includes gene quantification from featureCounts
#    .- scripts: Includes all scripts used to analyze the data
#    .- summary: Includes an html report generated with MultiQC with a summary of all
#                previous steps

# IMPORTANT REMARK: It is assumed that an initial raw_data folder has been previously
#         created during the data download process.

# NOTE: This is the main set of subfolders. Adapt this script for additional subfolders

#=========================
# User defined parameters
#=========================

ROOTDIR=$(sed -n 6p RNAseq_User_defined_parameters.txt)
WKD=$ROOTDIR/$(sed -n 12p RNAseq_User_defined_parameters.txt)

mkdir -p $WKD/raw_data/{FASTQC,MULTIQC}
mkdir -p $WKD/{trimmed_data/{FASTQC,MULTIQC},STAR_align/{Index,BAM,Other_results},quantification,scripts/{logs,cmds},summary}
