#!/bin/bash

#SBATCH --job-name=Summary_RNAseq
#SBATCH --partition=long
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1  
#SBATCH --output=logs/Summary_RNAseq.out
#SBATCH --error=logs/Summary_RNAseq.err

#=========================
# User defined parameters: relevant paths
#=========================

# Nothing to define

#=========================
# General configuration
#=========================

START=$(date +%s)
# Enable Singularity image to look into the general path (equivalent to -B)
export SINGULARITY_BIND=$ROOTDIR 
# Path to images folder in cluster
IMAGES_PATH=$ROOTDIR"/images"
# Path to databases folder in cluster
DB_PATH=$ROOTDIR"/db_files"

# Location where to add the summary output results
OUTMULTIQC=$WKD/'summary'

#=========================
# Singularity image and Tool Parametrization
#=========================

# Specify image/s name to be used (tool-related)
MULTIQC='multiqc_v1.12.sif' # This images includes MultiQC 1.12

#=========================
# Execution: Create a summary file gathering all results
#=========================

echo $WKD/'raw_data' > $OUTMULTIQC/'results.paths.txt'
echo $WKD/'STAR_align/Other_results' >> $OUTMULTIQC/'results.paths.txt'
echo $WKD/'quantification' >> $OUTMULTIQC/'results.paths.txt'


 # Build the final summary report

echo ''
echo 'Run MULTIQC.............................. '`date`
echo ''

cd $OUTMULTIQC

singularity exec $IMAGES_PATH/$MULTIQC multiqc --file-list results.paths.txt \
                  --ignore '*_ReadsPerGene.out.tab' --exclude 'snippy' \
                  --title $PROJECT': Summary report transcriptome profiling (RNAseq data)' \
                  --filename $PROJECT'_Summary_Transcriptome_Profiling.html'

#=========================
# End Preprocessing
#=========================

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Processing Time: $DIFF seconds"
