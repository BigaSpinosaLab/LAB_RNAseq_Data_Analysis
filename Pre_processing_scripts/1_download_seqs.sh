#!/bin/bash

#SBATCH --job-name=Download_Files 
#SBATCH --partition=long
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1  
#SBATCH --output=logs/download_files.out
#SBATCH --error=logs/download_files.err

ROOTDIR=$(sed -n 6p RNAseq_User_defined_parameters.txt)
WKD=$ROOTDIR/$(sed -n 12p RNAseq_User_defined_parameters.txt)

cd $WKD'/raw_data'

while read file
do
  #Download sequences from CRG server
  wget --user SPECIFY_USER --password SPECIFY_PWD http://seq.crg.es/download/external/SPECIFY_PATH/$file

done < $WKD'/scripts/CRG_files_example.txt'