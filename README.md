# RNA-seq data analysis
This repository includes scripts for analyzing 

the bulk RNAseq (both datasets included in the paper) and ChIPseq IkBa data analysis included in Alvarez-Villanueva et al. (Under peer-review). All scripts include comments so they are self-explanatory.

The repository is organized in the following subfolders:

Obtained raw counts matrix was imported into R 

## RNAseq data analysis folder

Scripts required to reproduce the complete RNAseq data analysis, specifically:

- Data preprocessing: to obtain a raw expression matrix from FASTQ files. Scripts from 0 to 7.
- Downstream analysis: to conduct differential expression analysis and functional analysis (GSEA and Overrepresentation analysis). Scripts 8 and 12. NOTE: Although these scripts are particularized for RNAseq inducible IkBa dataset, they have also been applied for the RNAseq Knock-in dataset. Experimental design has been adapted for the corresponding comparisons.




# Proposed methods

Quality control was performed on raw data with FASTQC tool (v0.11.9). Raw reads were trimmed to remove adapters presence with Trimgalore (v0.6.6). Default parameters were used except for a minimum quality of 15 (Phred score) and an adapter removal stringency of 3bp overlap. Trimmed reads were aligned to reference genome with STAR aligner tool (v2.7.8). STAR was executed with default parameters except for the number of allowed mismatches which was set to 1. Required genome index was built with corresponding [Specify org assembly] gtf and fasta files retrieved from Ensembl (http://ftp.ensembl.org/pub/release-XXX/ Specify the specific release). Obtained BAM files with uniquely mapped reads were considered for further analysis. Raw gene expression was quantified using featureCounts tool from subRead software (v2.0.1) with exon as feature.

References for required tools
FASTQC
Trimgalore
STAR aligner
featureCounts
