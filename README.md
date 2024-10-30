# Computational analyses for: SMG1:SMG8:SMG9-complex integrity maintains robustness of nonsense-mediated mRNA decay 
___
This repository contains the codes, scripts and log files for the high-throughput sequencing analyses of the project: <br/> [__SMG1:SMG8:SMG9-complex integrity maintains robustness of nonsense-mediated mRNA decay__](https://doi.org/10.1101/2024.04.15.589496) <br/> (available as bioRxiv preprint)

Preprint:  
[![DOI](https://img.shields.io/badge/DOI-10.1101%2F2024.04.15.589496%20-red)](https://doi.org/10.1101/2024.04.15.589496) 

## Graphical abstract

<img src="https://github.com/boehmv/2024_SMG189/blob/main/doc/Kueckelmann_2024_Model.png?raw=false" max-height="300">

## Scope
This repository primarily aims to provide transparent insight into the high-throughput analysis steps used in the study of SMG8- and/or SMG9-KO, SMG8-deltaKID and SMG1i treatment in human cells (all high-throughput data obtained from colon cancer HCT116 cell line). :warning: **NOTE:** The complete pipeline is currently not optimized to run on different computing infrastructures in a standardized/portable manner. This means that all required packages have to be installed manually and configured accordingly to reproduce the results.

## Features / Requirements
* Complete analysis of multiple RNA-Seq datasets (provided in FASTQ format; see [here](https://github.com/boehmv/2024_SMG189/blob/main/metadata/SMG189_datasets.csv) for dataset overview and [here](https://github.com/boehmv/2024_SMG189/blob/main/metadata/SMG189_samples.csv) for individual sample identification), mapped to [Gencode v42](https://www.gencodegenes.org/human/release_42.html) / GRCh38.primary_assembly supplemented with SIRVomeERCCome (from Lexogen; [download](https://www.lexogen.com/wp-content/uploads/2018/08/SIRV_Set3_Sequences_170612a-ZIP.zip)) using STAR, followed by transcript quantification using Salmon in mapping-based mode with a decoy-aware transcriptome index and the options --numGibbsSamples 30 --useVBOpt --gcBias --seqBias, finished with analyses of differential gene expression (DGE) via DESeq2 and differential transcript expression (DTE) via Swish (pre-revision) or edgeR (post-revision).
* The main Bash script [CRSA_V009.sh](https://github.com/boehmv/2024_SMG189/blob/main/scripts/CRSA_V009.sh) or [CRSA_V010.sh](https://github.com/boehmv/2024_SMG189/blob/main/scripts/CRSA_V010.sh) runs the complete pipeline or individual modules using the options (see CRSA_V009.sh -h) and requires a design file specifying the following:
  * reference type (gencode.v42.SIRVomeERCCome was used in this study)
  * sequencing design (single- or paired-end reads)
  * study name
  * folder locations (srvdir for raw file locations, mydir for analyses output)
  * location of the experiment file which specifies sample IDs and condition
* Please see the provided [design.txt](https://github.com/boehmv/2024_SMG189/blob/main/data/design.txt) file example for more information concerning this design file. An example for the tab-delimited [experiment.txt](https://github.com/boehmv/2024_SMG189/blob/main/data/experiment.txt) file is provided as well. Please see the comments in CRSA_V009.sh or CRSA_V010.sh for further instructions 
* To run/reproduce the complete analysis script, many modules require specific tools. Please make sure you have the following tools installed and configured if required:
  * [STAR](https://github.com/alexdobin/STAR) - version 2.7.10b was used for the analyses - with genome indices generated using GRCh38.primary.SIRVomeERCCome.fa and gencode.v42.SIRVomeERCCome.annotation.gtf (both reference files can be found [here](https://uni-koeln.sciebo.de/s/RFID1U3YYBZmkkE)). The following code was used for genome index generation: 
  ```
  STAR   --runMode genomeGenerate   --runThreadN 15   --genomeDir /home/volker/reference/gencode.v42.SIRVomeERCCome  --genomeFastaFiles /home/volker/reference/Gencode/GRCh38.primary.SIRVomeERCCome.fa --sjdbGTFfile /home/volker/reference/Gencode/gencode.v42.SIRVomeERCCome.annotation.gtf   --sjdbOverhang 99
  ```
  * [Alfred](https://github.com/tobiasrausch/alfred) - version v0.2.6 was used for the analyses
  * [samtools](http://www.htslib.org/) - version 1.16.1 (using htslib 1.16) was used for the analyses
  * [IGV tools](http://software.broadinstitute.org/software/igv/download) - version 2.14.1 or 2.17.2 was used for the analyses - make sure you have the gencode.v42.SIRVomeERCCome.chrom.sizes file (can be found [here](https://uni-koeln.sciebo.de/s/RFID1U3YYBZmkkE)) located in /PATH/TO/IGV/lib/genomes
  * [Salmon](https://github.com/COMBINE-lab/salmon) - version v1.9.0 was used for the analyses - with an index generated using gentrome.v42.SIRV.ERCC.fa.gz and decoys.txt (can be found [here](https://uni-koeln.sciebo.de/s/RFID1U3YYBZmkkE)). A separate conda environment was created for Salmon. The following code was used for index generation: 
  ```
  salmon index -t /home/volker/reference/Gencode/gentrome.v42.SIRV.ERCC.fa.gz -d /home/volker/reference/Gencode/decoys.txt -p 12 -i /home/volker/reference/Transcriptome/gencode.v42.SIRVomeERCCome --gencode
  ```
  * [DESeq2](https://github.com/mikelove/DESeq2) - version 1.40.1 was used for the analyses. The tx2gene file used for the analyses can be found [here](https://uni-koeln.sciebo.de/s/RFID1U3YYBZmkkE) 
  * [FastQC](https://github.com/s-andrews/FastQC) - version 0.11.9 was used for the analyses
  * [MultiQC](https://github.com/ewels/MultiQC) - version v1.14  was used for the analyses
* Additionally, many analyses were run using a plethora of R packages (including swish, edgeR, ...), please see the session info for the individual R scripts for more information.
* All analyses were performed on a 16-core (2x Intel(R) Xeon(R) CPU E5-2687W v2 @ 3.40GHz) workstation with 128 GB RAM running Ubuntu 22.04.2 LTS
* Please make sure to change installation and file paths in the respective scripts to match your local environment

## Individual scripts
The specialized scripts called by the main CRSA_V009.sh (pre-revision) or CRSA_V010.sh (post-revision) script can be found [here](https://github.com/boehmv/2024_SMG189/tree/main/scripts).
The main R script to produce the RNA-Seq-based figures can be found [here](https://github.com/boehmv/2024_SMG189/tree/main/scripts/SMG89_Revision_Analysis.R). This script uses the "SMG89_datasets.csv" file in the same folder to load DESeq2 or edgeR output data.
The main R script to produce the proteomics-based figures can be found [here](https://github.com/boehmv/2024_SMG189/tree/main/scripts/SMG89_Proteomics_Analysis.R). This script uses the "SMG89_proteomics.xlsx" file in the data folder.

## Note on DTE analyses
We have previously used swish to perform DTE analyses, but switched now to edgeR. That is why e.g. the absolute numbers of differentially expressed transcripts in the bioRxiv version will differ from the one in the final version of the manuscript. The main reason was that we made use of the edgeR function to better control “read-to-transcript ambiguity” and false discovery rate (see: https://doi.org/10.1093/nar/gkad1167).

## Output data 
The required Salmon, DESeq2 and edgeR output data to re-run most of the analyses can be found on [here](https://doi.org/10.5281/zenodo.14003499)
Several helper files, Log files and QC-related data, as well as an "ready-to-load-in-R" 2024-10-28_SMG189_datasources.rds file can be found there as well.

## Feedback / Questions
Feedback is welcome! For any question, please email: boehmv@uni.koeln.de or [create an issue](https://github.com/boehmv/2024_SMG189/issues)

## Citation
### Journal article
TBD

### bioRxiv preprint
Sabrina Kueckelmann, Sophie Theunissen, Jan-Wilm Lackmann, Marek Franitza, Kerstin Becker, Volker Boehm, Niels H. Gehring (2024) __SMG1:SMG8:SMG9-complex integrity maintains robustness of nonsense-mediated mRNA decay__. 
bioRxiv 2024.04.15.589496; doi: [https://doi.org/10.1101/2024.04.15.589496](https://doi.org/10.1101/2024.04.15.589496)
