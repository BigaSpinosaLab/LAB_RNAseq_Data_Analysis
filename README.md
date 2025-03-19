# RNA-seq data analysis
This repository includes scripts for pre-processing RNAseq data to obtain a raw expression matrix from a set of FASTQ files using a splice-aware aligner tool. The final raw expression matrix can be imported into i.e. R for downstream statistical analyisis.

Considerations to bear in mind:

- All scripts include comments so they are self-explanatory. Nevertheless, below is a brief explanation of how scripts are structured.
- Scripts are prepared to be executed in a HPC environment with SLURM. Even in a similar context, the SLURM directive should be adapted to the specific HPC used and user requirements.
- All scripts required a singularity image containing a specific tool to be run. They can be obtained by executing a singularity recipe or pulling an existing Docker container.
- There is one specific script per analysis step and all of them should be sequentially run. With SLURM, it is possible to run them sequentially using the parameter --dependency when submitting a job to the HPC.

## General Pipeline and scripts structure 

The scripts folder contains the folder bash scripts to complete this part of the analysis. There is one script per specific step in the pipeline, specifically:

- 0_Create_dir_tree_RNAseq_STAR.sh: used to create the proper directory tree for storing intermediate and final results. Rest of scripts rely on this structure.
- 1_download_seqs.sh: used to download files from a sequencing platform (particularized to CRG seq platfom at PRBB, Barcelona).
- 2_Quality_check_data.sh: used to run FASTQ quality check with basic metrics. This can be used for raw and trimmed data.
- 3_Trimming_trimgalore.sh: used to run reads trimming. This is optional depending on the quality/adapters presence in your reads.
- 4a_STAR_Build_Index.sh: used to build the required STAR index for reads alignment. This is optional since you may already have an index.
- 4b_STAR_alignment.sh: used to align FASTQ files with STAR aligner tool.
- 5_featureCounts_Exprs_Quantification.sh: used to quantify raw counts per gene in your annotation file per sample under test.
- 6_Summary_RNAseq_profiling.sh: used to generate a summary report including the results per each step.

Importantly, there is a txt file 'RNAseq_User_defined_parameters.txt' that MUST be adapted for each RNAseq data analysis project. It basically contains the required data paths.

### Scripts structure

All scripts follow the same structure with following sections:

- User defined parameters. This sections specifies the required parameters to define. Most of the scripts do not require anything.
- General configuration. This defines correct paths and sets singularity binding. User do not have to change anything in this part.
- Singularity image and tool parametrization. Definition of the singularity image to be used and the tools parameters to consider.
- Execution. Execution of the corresponding tool using previous information.

## Proposed methods section

Following paragraph can be used in a methods section to explain RNA-seq data pre-processing analysis. Information can be adapted/customized in any case.

"Quality control was performed on raw data with FASTQC tool (v0.11.9). Raw reads were trimmed to remove adapters presence with Trimgalore (v0.6.6). Default parameters were used except for a minimum quality of 15 (Phred score) and an adapter removal stringency of 3bp overlap. Trimmed reads were aligned to reference genome with STAR aligner tool (v2.7.8). STAR was executed with default parameters except for the number of allowed mismatches which was set to 1. Required genome index was built with corresponding [*Specify org assembly*] gtf and fasta files retrieved from Ensembl (*Specify the specific release XXX* i.e. http://ftp.ensembl.org/pub/release-XXX/). Obtained BAM files with uniquely mapped reads were considered for further analysis. Raw gene expression was quantified using featureCounts tool from subRead software (v2.0.1) with exon as feature."

References for required tools:

- FASTQC: Andrews S.,FASTQC: a quality control tool for high throughput sequence data. https://github.com/s-andrews/FastQC
- Trimgalore: Felix Krueger, Frankie James, Phil Ewels, Ebrahim Afyounian, & Benjamin Schuster-Boeckler. (2021). FelixKrueger/TrimGalore: v0.6.7 - DOI via Zenodo (0.6.7). Zenodo. https://doi.org/10.5281/zenodo.5127899
- STAR aligner: Alexander Dobin, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, Thomas R. Gingeras, STAR: ultrafast universal RNA-seq aligner, Bioinformatics, Volume 29, Issue 1, January 2013, Pages 15–21, https://doi.org/10.1093/bioinformatics/bts635
- featureCounts: Liao, Y., Smyth, G. K. & Shi, W. featureCounts: an efficient general purpose program for assigning sequence reads  to genomic features. Bioinformatics 30, 923–930 (2014).
