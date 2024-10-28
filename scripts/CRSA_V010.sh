#!/bin/bash

# Title: Complete RNA-Seq Analysis [CRSA] (Version 0.1.0)
# Objective: Run various tools to map FASTQ format RNA-Seq data (via STAR), count (via salmon) and perform downstream analyses (such as DGE, DTE, DTU, AS, IR), as well as QC (via e.g. FastQC and multiQC)
# Created by: boehmv (Volker Böhm; boehmv@uni-koeln.de)

########################
#			#
# Usage:		#
# CRSA_V010.sh -h	#
#			#
########################

# Version 0.1.0 (started on 10.07.24):
# - Changed Option to allow mapping with Salmon (with 30 Gibbs samples) against de Novo NMDRegHumanTxome transcriptome indices (-N)
# - Can use edgeR for DTE analysis on both GENCODE reference and NMDRegHumanTxome (-E)

# Version 0.0.9 (started on 15.08.23):
# - *** NEW DEFAULT FOR SALMON *** Bootstrap Salmon - uses 30 Gibbs samples -> allows DGE, DTE, DTU with Swish using inferential replicates
# - Swish analysis option available (-S)
# - Combined DESeq2 DGE and DTE into one script (-D) option only
# - Slightly different structure of DESeq2 and Swish outputs (DGE, DTE, DTU subfolders)
# - Option to allow mapping with Salmon (with 30 Gibbs samples) against StringTie-generated *NMD_HCT116_gencode.v42* indices (-N)
# - Option to quantify reads mapping to rRNA (via Bowtie2) (-r)

# Version 0.0.8 (24.10.22):
# - Enable usage of newest GENCODE annotations (v42) 
# - GENCODE files + "old" SIRVomeERCCome used to prepare new annotations, indices, etc.
# - Updated Salmon (from v1.3 -> v.1.9) -> new environment (!)
# - New Salmon transcriptome index (see /home/volker/reference/Gencode/v42/Prepare_Salmon_Index.sh)
# - Updated STAR (from 2.7.3a -> 2.7.10b) via git
# - Updated alfred (from 0.2.1 -> 0.2.6) via conda
# - Updated IGV (from 2.8.0 -> 2.14.1)
# - New STAR index generated
 
# Version 0.0.7.beta (08.02.22):
# - Enable -T option for DESeq2 DTE analysis

# Version 0.0.6 (04.08.20):
# - Use design file (design.txt) for study details, initialized with '-f /path/to/file' option (MANDATORY!). This allows to have one universal CRSA script for all analyses
# - Improved analyses and plotting scripts (DESeq2, ISAR and AS; plus combined analysis)
# - Updated alfred (v0.2.1) and salmon (1.3.0)
# - Newly generated salmon index with current version (1.3.0): /home/volker/reference/Transcriptome/gencode.v33.SIRVomeERCCome
#   - Code: salmon index -t /home/volker/reference/Gencode/gentrome.v33.SIRV.ERCC.fa.gz -d /home/volker/reference/Gencode/decoys.txt -p 12 -i /home/volker/reference/Transcriptome/gencode.v33.SIRVomeERCCome --gencode
# - Minor bugfixes

# Version 0.0.5 (06.05.20):
# - Changes in the leafcutter analysis: run leafcutter R analysis script -> produces volcano plots, calculates 5'ss and 3'ss MaxEnt scores and produces boxplots and volcano plots with these scores
# - IRFinder updated to version 1.2.6, change from BAM to FASTQ mode since this removes having to re-sort the BAM files
# - added option for single or paired end reads - changes in STAR, Salmon and IRFinder scripts

# Version 0.0.4 (26.04.20):
# - Minor changes to DESeq2 and ISAR R scripts, they now by default produce comparison plots (heatmap, barcode plots)
# - Error-resistant execution of g:profiler after DESeq2 analysis with static and interactive output
# - Integrate MultiQC as additional option (-Q), which scans for FastQC, Salmon and STAR files
# - [Optional] generate FastQC in FASTQ symlink step

# Version 0.0.3 (03.03.20):
# - Use Options via getopts (see below)
# - Better arrangement of the output to screen/logfile
# - Got rid of the "myprefix" and range variables, instead use experiment file directly via samples array
# - Updated STAR (2.7.3a) via download and stored @ /opt/STAR
# - Updated Alfred (0.1.19) via conda update
# - Updated IGV (2.8.0) via download and stored @  ~/Tools
# - Updated Salmon (1.1.0) via conda create -n salmon salmon
# - Use gencode.v33 (primary and transcripts) instead of plain ensembl
# - For human RNA-Seq data; Use SIRV and ERCC spike-ins by default (SIRVomeERCCome from SIRV_Set3_Sequences_170612a, modified!)
#   - SIRVomeERCCome annotation (GTF) modified:
#     - Change ERCC transcript_id (e.g. "DQ516748") to the ERCC Name (e.g. "ERCC-00046")
#     - Add "gene_type" identifier "spike_in" and "gene_name" identifier (= gene_id; e.g. "SIRV7")
#   - SIRV_ERCC fasta file modified:
#     - Use gffread to get SIRV individual transcript fasta sequences (needed for Salmon)
#     - Combine individual SIRV with ERCCs, but REMOVE "ERCC-00007", "ERCC-00018", "ERCC-00023" and "ERCC-00128" from "SIRV_ERCC_2.fa", as they are not annotated and will be problematic in the Salmon/ISAR steps
# - New concatenated fasta and gtf files in /home/volker/reference/Gencode/*SIRVomeERCCome*
# - New STAR genome index (/home/volker/reference/gencode.v33.SIRVomeERCCome)
# - New Salmon index (/home/volker/reference/Transcriptome/gentrome.v33.SIRVomeERCCome_index) using Decoys and --gencode
# - New DESeq2 tx2gene file (from gencode.v33.SIRVomeERCCome)
# - Adjusted ISAR conditions allowing spike_ins quantification
# - New Leafcutter exon_IDs based on gencode.v33.SIRVomeERCCome
# - New IRFinder REF based on gencode.v33.SIRVomeERCCome


# Version 0.0.2 (19.10.19):
# - More verbose on each step, giving version, etc.
# - Store BAM files by default on /srv2 (HDD) to save disk space on /home (SSD)
# - Remove bedGraph files to save disk space
# - Massively updated IRFinder script: Remove aligned files to safe disk space, use DESeq2 for IR quantification, get meaningful output files (xlsx, volcano plot, etc.)
# - DESeq2, leafcutter and IRFinder now work with updated "combine".py scripts - giving formatted final output excel files (in respective /Combined folder)
# - Updated leafcutter procedure: modified the reference to include ensembl_id, modified the leafcutter R script to output gene_ids, keep alternative genes, modified the downstream scripts

#################################################
#						#
# CRSA getopts: define and get options		#
#						#
#################################################

if [ $# -eq 0 ]
then
  echo "Missing options!"
  echo "(run $0 -h for help)"
  echo ""
  exit 0
fi

ECHO="false"

while getopts "ahMCNDSEIRLirQf:" OPTION; do
  case $OPTION in

    a)
    STAR_enable="true"
    Salmon_enable="true"
    Salmon_NMD_enable="true"
    DESeq2_enable="true"
    Swish_enable="true"
    edgeR_enable="true"
    ISAR_enable="true"
    ISAR_NMD_enable="true"
    leafcutter_enable="true"
    IRFinder_enable="true"
    rRNA_quant_enable="true"
    QC_enable="true"
    ;;

    M)
    STAR_enable="true"
    ;;

    C)
    Salmon_enable="true"
    ;;
    
    N)
    Salmon_NMD_enable="true"
    ;;

    D)
    DESeq2_enable="true"
    ;;    
    
    S)
    Swish_enable="true"
    ;;    
    
    E)
    edgeR_enable="true"
    ;;  

    I)
    ISAR_enable="true"
    ;;
    
    R)
    ISAR_NMD_enable="true"
    ;;

    L)
    leafcutter_enable="true"
    ;;

    i)
    IRFinder_enable="true"
    ;;
    
    r)
    rRNA_quant_enable="true"
    ;;

    Q)
    QC_enable="true"
    ;;

    f)
    design_file=$OPTARG
    if [ -e "${design_file}" ]; then

	# Read in design file and variables therein
	source $design_file

	# Check if variables are defined
	if [[ -z $reference || -z $seq_design || -z $myname || -z $srvdir || -z $mydir || -z $experiment_file ]]; then
		echo "One or more variable are undefined"
	else
		ECHO="true"
	fi
    else
	echo "Design file ${design_file} does not exist"
    fi    
    ;;

    h)
    echo "Usage:"
    echo "CRSA_V010.sh -f /path/to/design.txt [-OPTION] "
    echo ""
    echo "   -f path/to/design.txt	MANDATORY: give this path to the design file; please see example design file for guidance"
    echo "   -a     to execute the full analysis (all other options enabled)"
    echo "   -h     help (this output)"
    echo "   -M     to execute STAR"
    echo "   -C     to execute Salmon - 30 Gibbs Sample"
    echo "   -N     to execute Salmon - NMDRegHumanTxome assembly"
    echo "   -D     to execute DESeq2"
    echo "   -S     to execute Swish-based DGE, DTE and DTU analysis"
    echo "   -E     to execute edgeR-based DTE analysis"
    echo "   -I     to execute ISAR"
    echo "   -R     to execute ISAR - NMDRegHumanTxome assembly"
    echo "   -L     to execute leafcutter"
    echo "   -i     to execute IRFinder"
    echo "   -r     to execute rRNA quantification with Bowtie2"
    echo "   -Q     to execute QC"
    exit 0
    ;;

  esac
done

########################################
#					#
#     ---- Start of script ----	#
#					#
########################################

# Only proceed if acceptable options were given (NOTE: this is one big 'if' command until the end of the script!)

if [ $ECHO = "true" ]
then

# Get current date
date=$(date "+%Y-%m-%d")

# Define log file and redirect stdout and stderr to this file
if [ ! -d "${mydir}/Logs/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${mydir}/Logs/"
  mkdir ${mydir}/Logs/
fi
log_file="${mydir}/Logs/log_$date"
exec &> >(tee -a "$log_file")

##############################################################################
# Give absolute paths to local folders of STAR, leafcutter and IRFinder 	#
# (change according to installation path)					#
##############################################################################

STAR_dir="/home/volker/reference"
leafcutter_dir="/home/volker/Tools/leafcutter"
IRFinder_dir="/home/volker/Tools/IRFinder"
DESeq2_dir="/home/volker/Tools/DESeq2"
Swish_dir="/home/volker/Tools/Swish"
edgeR_dir="/home/volker/Tools/edgeR"
ISAR_dir="/home/volker/Tools/ISAR"
bowtie_rRNA_dir="/home/volker/reference/Bowtie2/rRNA_all/rRNA_all"

# State Script Version
echo ""
echo "##########################################################"
echo "## Complete RNA-Seq Analysis [CRSA] Version 0.1.0 	##"
echo "##########################################################"
echo ""

# Show design parameters in log
echo "######################"
echo "## General settings ##"
echo "######################"
echo ""
echo "Reference used: $reference"
echo "Sequencing design (paired/single): $seq_design"
echo "Name of the study: $myname"
echo "Location of files on server: $srvdir"
echo "Location of project folder: $mydir"
echo "Experiment file used for this analysis: $experiment_file"
echo ""

# Declare samples array
declare -a samples

# Load experiment file, parse for samples and save them into the array
let p=0
while read -r f1 f2; do
  samples[p]="${f1}"
  ((++p))
done < $experiment_file

# Declare condition arrays
declare -a cond

# Load experiment file, parse for conditions and save unique conditions into array. Note: Condition1 is always "control"
let i=1
cond[0]="control"
while read -r f1 f2; do
  if [[ " ${cond[*]} " == *"$f2"* ]];then
  continue
else
  cond[i]="${f2}"
  ((++i))
fi
done < $experiment_file

# Declare individual condition arrays
arr_length="$((${#cond[@]}-1))"
for i in $( eval echo {0..${arr_length}} );do
  declare -a cond${i}
done

# Load experiment file again, parse for conditions and save filenames into condition-specific arrays.
while read -r f1 f2; do
  for i in $( eval echo {0..${arr_length}} );do
  if [[ "$f2" == "${cond[i]}" ]];then
    eval cond${i}[cond${i}count]="${f1}"
    ((++cond${i}count))
  fi
done
done < $experiment_file

# State the conditions and samples for this analysis
echo "############################"
echo "## Conditions and samples ##"
echo "############################"
echo ""

arr_length="$((${#cond[@]}-1))"
for i in $( eval echo {0..${arr_length}} );do
  echo -e "cond${i} \t ${cond[i]} \t $(eval echo \${cond$i[*]})"
done

#Generate folders if necessary

# Generate experiment folder in srvdir
if [ ! -d "${srvdir}/${myname}/" ]; then
  mkdir ${srvdir}/${myname}/
fi

# Generate srv folder
if [ ! -d "${srvdir}/" ]; then
  mkdir ${srvdir}/
fi

# Generate srv myname folder
if [ ! -d "${srvdir}/${myname}/" ]; then
  mkdir ${srvdir}/${myname}/
fi

# Generate srv BAM folder
if [ ! -d "${srvdir}/${myname}/BAM/" ]; then
  mkdir ${srvdir}/${myname}/BAM/
fi

# Generate home BAM folder
if [ ! -d "${mydir}/BAM/" ]; then
  mkdir ${mydir}/BAM/
fi

# Generate Tracks folder for easy visualization
if [ ! -d "${mydir}/Tracks/" ]; then
  mkdir ${mydir}/Tracks/
fi

# Create complete sample text file
> ${mydir}/Samples.txt
echo -e "sample \t condition" >> ${mydir}/Samples.txt
cat $experiment_file >> ${mydir}/Samples.txt

# Return state of analysis options
echo ""
echo "######################"
echo "## Analysis options ##"
echo "######################"
echo ""
echo "STAR_enable (-M) = ${STAR_enable}"
echo "Salmon_enable (-C) = ${Salmon_enable}"
echo "Salmon_NMD_enable (-N) = ${Salmon_NMD_enable}"
echo "DESeq2_enable (-D) = ${DESeq2_enable}"
echo "Swish_enable (-S) = ${Swish_enable}"
echo "edgeR_enable (-E) = ${edgeR_enable}"
echo "ISAR_enable (-I) = ${ISAR_enable}"
echo "ISAR_NMD_enable (-R) = ${ISAR_NMD_enable}"
echo "leafcutter_enable (-L) = ${leafcutter_enable}"
echo "IRFinder_enable (-i) = ${IRFinder_enable}"
echo "rRNA_quant_enable (-r) = ${rRNA_quant_enable}"
echo "QC_enable (-Q) = ${QC_enable}"
echo ""

######
######################
#
# Perform STAR mapping, indexing, bedGraph (via Alfred) and tracks (via IGV) generation
#
######################
######

if  [ "$STAR_enable" = "true" ]; then

echo ""
echo "#########################################################"
echo "## Perform STAR mapping, indexing and track generation ##"
echo "#########################################################"
echo ""
echo "STAR directory: $STAR_dir"
echo -n "STAR version: "
STAR --version
alfred --version

# State IGV tools version
/home/volker/Tools/IGV_2.14.1/igvtools version


for i in "${samples[@]}"; do
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo ""
  echo "$mytime Start STAR mapping of sample ${i}"
  echo ""

  # Check Sequencing design and designate variable accordingly
  if  [ "$seq_design" = "paired" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime sample ${i} paired reads detected"
    seq_reads="$mydir/FASTQ/${i}_1.fq.gz $mydir/FASTQ/${i}_2.fq.gz"
  elif  [ "$seq_design" = "single" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime sample ${i} single reads detected"
    seq_reads="$mydir/FASTQ/${i}.fq.gz"
  fi

  # Generate individual folder
  if [ ! -d "${srvdir}/${myname}/BAM/${i}/" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime Make directory ${srvdir}/${myname}/BAM/${i}/"
    mkdir ${srvdir}/${myname}/BAM/${i}/
  fi

  # Check if BAM file was already generated
  if [ -e "${srvdir}/${myname}/BAM/${i}/Aligned.sortedByCoord.out.bam" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime ${i} BAM file already generated"
  else
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime ${i} mapping in process"
    # Modified conditions
    STAR \
    --runThreadN 15 \
    --quantMode GeneCounts \
    --genomeDir $STAR_dir/$reference \
    --genomeLoad NoSharedMemory \
    --readFilesIn $seq_reads \
    --readFilesCommand gunzip -c \
    --outReadsUnmapped Fastx \
    --outSJfilterOverhangMin 15 15 15 15 \
    --alignSJoverhangMin 15 \
    --alignSJDBoverhangMin 10 \
    --outFilterMultimapNmax 20 \
    --outFilterScoreMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.05 \
    --outFilterMatchNminOverLread 0.7 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --chimSegmentMin 15 \
    --chimScoreMin 15 \
    --chimScoreSeparation 10 \
    --chimJunctionOverhangMin 15 \
    --twopassMode Basic \
    --chimOutType SeparateSAMold \
    --alignSoftClipAtReferenceEnds No \
    --outSAMattributes NH HI AS nM NM MD jM jI XS \
    --outFileNamePrefix ${srvdir}/${myname}/BAM/${i}/ \
    --outSAMtype BAM SortedByCoordinate
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime ${i} mapping complete"
  fi

  # Index BAM file if necessary
  if [ -e "${srvdir}/${myname}/BAM/${i}/Aligned.sortedByCoord.out.bam.bai" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime ${i} BAM file already indexed"
  else
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime Index ${srvdir}/${myname}/BAM/${i}"

    samtools index \
    ${srvdir}/${myname}/BAM/${i}/Aligned.sortedByCoord.out.bam \
    ${srvdir}/${myname}/BAM/${i}/Aligned.sortedByCoord.out.bam.bai
  fi

  # Generate symbolic bam link for easy access
  if [ -e "$mydir/BAM/${i}.bam" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime ${i} BAM file link already created"
  else
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime ${i} BAM file symlink created"
    ln -s ${srvdir}/${myname}/BAM/${i}/Aligned.sortedByCoord.out.bam $mydir/BAM/${i}.bam
  fi

  # Generate symbolic index link for easy access
  if [ -e "$mydir/BAM/${i}.bam.bai" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime ${i} BAM file index link already created"
  else
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime ${i} BAM file index symlink created"
    ln -s ${srvdir}/${myname}/BAM/${i}/Aligned.sortedByCoord.out.bam.bai $mydir/BAM/${i}.bam.bai
  fi

  # Generate bedGraph using Alfred and convert to tdf using IGVtools
  if [ -e "$mydir/Tracks/${i}.tdf" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime ${i} BAM file already converted to track"
  else
    if [ ! -d "${mydir}/Tracks/BedGraph/" ]; then
      mytime=$(date "+%Y-%m-%d %H:%M:%S")
      echo "$mytime Make directory ${mydir}/Tracks/BedGraph/"
      mkdir ${mydir}/Tracks/BedGraph/
    fi

    # Use Alfred to generate bedGraphs
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime ${i} create BedGraph"
    alfred tracks -r 1 -o $mydir/Tracks/BedGraph/${i}.bedGraph.gz $mydir/BAM/${i}.bam
    cd /home/volker/Tools/IGV_2.14.1/

    # Convert to tdf
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime ${i} convert to TDF"
    /home/volker/Tools/IGV_2.14.1/igvtools totdf $mydir/Tracks/BedGraph/${i}.bedGraph.gz $mydir/Tracks/${i}.tdf $reference
    # Remove bedGraph files to save disk space
    #rm $mydir/Tracks/BedGraph/${i}.bedGraph.gz
  fi
done
fi

######
######################
#
# Perform salmon analysis
#
######################
######

if  [ "$Salmon_enable" = "true" ]; then

source /home/volker/miniconda3/bin/activate salmon_recent

echo ""
echo "#############################"
echo "## Perform salmon bootstrap analysis ##"
echo "#############################"
echo ""
echo -n "Salmon version: "
salmon -v
echo ""


# Iterate over fastq files and quantify with Salmon

for i in "${samples[@]}"; do
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Start Salmon analysis ${i}"

  if [ ! -d "${mydir}/Salmon/" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime Make directory ${mydir}/Salmon/"
    mkdir ${mydir}/Salmon/
  fi

  if [ -d "$mydir/Salmon/${i}" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime Sample ${i} already analyzed"
  else

    # Paired end reads
   if  [ "$seq_design" = "paired" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime Paired-end sample ${i} quant in progress"
    salmon quant -i /home/volker/reference/Transcriptome/$reference -l A \
    -1 $mydir/FASTQ/${i}_1.fq.gz \
    -2 $mydir/FASTQ/${i}_2.fq.gz \
    -p 15 \
    --numGibbsSamples 30 \
    --useVBOpt \
    --gcBias \
    --seqBias \
    --validateMappings \
    -o $mydir/Salmon/${i}

  # Single end reads
  elif  [ "$seq_design" = "single" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime Single-end sample ${i} quant in progress"
    salmon quant -i /home/volker/reference/Transcriptome/$reference -l A \
    -r $mydir/FASTQ/${i}.fq.gz \
    -p 15 \
    --numGibbsSamples 30 \
    --useVBOpt \
    --gcBias \
    --seqBias \
    --validateMappings \
    -o $mydir/Salmon/${i}
  fi
fi
done

conda deactivate

fi


if  [ "$Salmon_NMD_enable" = "true" ]; then

source /home/volker/miniconda3/bin/activate salmon_recent

echo ""
echo "#############################"
echo "## Perform salmon NMD assembly analysis ##"
echo "#############################"
echo ""
echo -n "Salmon version: "
salmon -v
echo ""


# Iterate over fastq files and quantify with Salmon

for i in "${samples[@]}"; do
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Start Salmon analysis ${i}"

  if [ ! -d "${mydir}/Salmon_NMD/" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime Make directory ${mydir}/Salmon_NMD/"
    mkdir ${mydir}/Salmon_NMD/
  fi

  if [ -d "$mydir/Salmon_NMD/${i}" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime Sample ${i} already analyzed"
  else

    # Paired end reads
   if  [ "$seq_design" = "paired" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime Paired-end sample ${i} quant in progress"
    salmon quant -i /home/volker/reference/Transcriptome/NMDRegHumanTxome.v1.1.SIRVomeERCCome -l A \
    -1 $mydir/FASTQ/${i}_1.fq.gz \
    -2 $mydir/FASTQ/${i}_2.fq.gz \
    -p 15 \
    --numGibbsSamples 30 \
    --useVBOpt \
    --gcBias \
    --seqBias \
    --validateMappings \
    -o $mydir/Salmon_NMD/${i}

  # Single end reads
  elif  [ "$seq_design" = "single" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime Single-end sample ${i} quant in progress"
    salmon quant -i /home/volker/reference/Transcriptome/NMDRegHumanTxome.v1.1.SIRVomeERCCome -l A \
    -r $mydir/FASTQ/${i}.fq.gz \
    -p 15 \
    --numGibbsSamples 30 \
    --useVBOpt \
    --gcBias \
    --seqBias \
    --validateMappings \
    -o $mydir/Salmon_NMD/${i}
  fi
fi
done

conda deactivate

fi

######
####################################
#
# Perform DESeq2 analysis
#
####################################
######

if  [ "$DESeq2_enable" = "true" ]; then

echo ""
echo "#############################"
echo "## Perform DESeq2 analysis ##"
echo "#############################"
echo ""
echo "DESeq2 directory: $DESeq2_dir"
echo ""

	# Generate DESeq2 folder if necessary
	if [ ! -d "${mydir}/DESeq2/" ]; then
	  mytime=$(date "+%Y-%m-%d %H:%M:%S")
	  echo "$mytime Make directory ${mydir}/DESeq2/"
	  mkdir ${mydir}/DESeq2/
	fi
	
	# Run DESeq2 script using mydir as arguments
	$DESeq2_dir/DESeq2_${reference}.R ${mydir} 


fi

######
####################################
#
# Perform Swish DGE, DTE and DTU analysis
#
####################################
######

if  [ "$Swish_enable" = "true" ]; then

echo ""
echo "#################################"
echo "## Perform Swish DGE, DTE and DTU analysis ##"
echo "#################################"
echo ""
echo "Swish directory: $Swish_dir"
echo ""

# Generate Swish folder if necessary
if [ ! -d "${mydir}/Swish/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${mydir}/Swish/"
  mkdir ${mydir}/Swish/
fi

##
# Run self-made Swish scripts using mydir as positional arguments
##

# DGE
$Swish_dir/swish_DGE_${reference}.R ${mydir} 

# DTE
$Swish_dir/swish_DTE_${reference}.R ${mydir} 

# DTU
$Swish_dir/swish_DTU_${reference}.R ${mydir} 

fi

######
####################################
#
# Perform edgeR DTE analysis
#
####################################
######

if  [ "$edgeR_enable" = "true" ]; then

echo ""
echo "#################################"
echo "## Perform edgeR DTE analysis ##"
echo "#################################"
echo ""
echo "edgeR directory: $edgeR_dir"
echo ""

# Generate Swish folder if necessary
if [ ! -d "${mydir}/edgeR/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${mydir}/edgeR/"
  mkdir ${mydir}/edgeR/
fi

##
# Run self-made edgeR scripts using mydir as positional arguments
##

# DTE reference
$edgeR_dir/edgeR_DTE_${reference}.R ${mydir} 

	if [ ! -d "${mydir}/Salmon_NMD/" ]; then
	  mytime=$(date "+%Y-%m-%d %H:%M:%S")
	  echo "$mytime Salmon not performed on NMDRegHumanTxome"
	  
	  else
	  # DTE NMDRegHumanTxome
	  $edgeR_dir/edgeR_DTE_NMDRegHumanTxome.R ${mydir} 
	fi
fi

######
####################################
#
# Perform Isoform Switch AnalyseR (ISAR) analysis
#
####################################
######

if  [ "$ISAR_enable" = "true" ]; then

echo ""
echo "###########################"
echo "## Perform ISAR analysis ##"
echo "###########################"
echo ""
echo "IsoformSwitchAnalyzeR directory: $ISAR_dir"

# Generate ISAR folder if necessary
if [ ! -d "${mydir}/ISAR/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${mydir}/ISAR/"
  mkdir ${mydir}/ISAR/
fi

# Organize conditions from array in a single string
cond_string=$(printf "%s," "${cond[@]}" | cut -d "," -f 1-${#cond[@]})

# Run self-made ISAR script using mydir and the cond_string as positional arguments
$ISAR_dir/ISAR_${reference}.R ${mydir} ${cond_string}

# Generate ISAR-specific condition table for running the comparison script
> ${mydir}/ISAR/ConditionTable.txt
echo -e "cond\tpath" >> ${mydir}/ISAR/ConditionTable.txt
echo -e "${cond[@]:1}\t${mydir}/ISAR/SwitchList_filt_Analyzed.csv" >> ${mydir}/ISAR/ConditionTable.txt

# Run ISAR Comparison script
$ISAR_dir/ISAR_Comparison.R ${mydir} ${mydir}/ISAR/ConditionTable.txt
fi

######
####################################
#
# Perform NMDRegHumanTxome - Isoform Switch AnalyseR (ISAR) analysis
#
####################################
######

if  [ "$ISAR_NMD_enable" = "true" ]; then

echo ""
echo "###############################################"
echo "## Perform NMDRegHumanTxome - ISAR analysis ##"
echo "###############################################"
echo ""
echo "IsoformSwitchAnalyzeR directory: $ISAR_dir"

# Generate ISAR folder if necessary
if [ ! -d "${mydir}/ISAR_NMD/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${mydir}/ISAR_NMD/"
  mkdir ${mydir}/ISAR_NMD/
fi

# Organize conditions from array in a single string
cond_string=$(printf "%s," "${cond[@]}" | cut -d "," -f 1-${#cond[@]})

# Run self-made ISAR script using mydir and the cond_string as positional arguments
$ISAR_dir/ISAR_NMDRegHumanTxome.R ${mydir} ${cond_string}

# Generate ISAR-specific condition table for running the comparison script
> ${mydir}/ISAR_NMD/ConditionTable.txt
echo -e "cond\tpath" >> ${mydir}/ISAR_NMD/ConditionTable.txt
echo -e "${cond[@]:1}\t${mydir}/ISAR_NMD/SwitchList_filt_Analyzed.csv" >> ${mydir}/ISAR_NMD/ConditionTable.txt

# Run ISAR Comparison script
# $ISAR_dir/ISAR_Comparison.R ${mydir} ${mydir}/ISAR/ConditionTable.txt
fi

######
####################################
#
# Perform leafcutter analysis
#
####################################
######

if  [ "$leafcutter_enable" = "true" ]; then

echo ""
echo "######################################################"
echo "## Perform leafcutter alternative splicing analysis ##"
echo "######################################################"
echo ""
echo "Leafcutter directory: $leafcutter_dir"

# Step 1: Perform the junction counting step

# Generate necessary folders
if [ ! -d "${mydir}/leafcutter/" ]; then
  mkdir ${mydir}/leafcutter/
fi
if [ ! -d "${mydir}/leafcutter/Juncfile" ]; then
  mkdir ${mydir}/leafcutter/Juncfile
fi

# Iterate over bamfiles
for bamfile in `ls ${mydir}/BAM/*.bam`
do

  # Get simplified filenames
  filename=$(basename -- "$bamfile")
  filename=${filename%%.*}

  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo $mytime Converting $filename to $filename.junc

  # Check if BAM file was already counted
  if [ -e "${mydir}/leafcutter/Juncfile/$filename.junc" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime BAM file already converted to junc file"

    # Count junctions
  else
  	/home/volker/Tools/regtools/build/regtools junctions extract -a 8 -m 50 -M 500000 -s XS $bamfile -o ${mydir}/leafcutter/Juncfile/$filename.junc

  fi
done

# Loop over all conditions vs control and finish step 1-3
for j in $( eval echo {1..${arr_length}} );do

  # Define current analysis (control vs ?)
  myexp="${cond[0]}_vs_$(eval echo \${cond[$j]})"

  # Generate current analysis folder
  if [ ! -d "${mydir}/leafcutter/$myexp/" ]; then
    mkdir ${mydir}/leafcutter/$myexp/
  fi

  # Generate text files for each analysis
  > ${mydir}/leafcutter/Juncfile/${date}_juncfile_$myexp.txt
  > ${mydir}/leafcutter/$myexp/Groups_$myexp.txt

  # Iterate over bamfiles and fill essential information in text files
  for bamfile in `ls ${mydir}/BAM/*.bam`
  do

    # Get simplified filenames
    filename=$(basename -- "$bamfile")
    filename=${filename%%.*}

    # Output juncfile information to necessary textfiles if sample belongs to respective conditions
    if [[ " ${cond0[@]} " == *"$filename"* ]];then
      echo -e $filename'\t'${cond[0]} >> ${mydir}/leafcutter/$myexp/Groups_$myexp.txt
      echo ${mydir}/leafcutter/Juncfile/$filename.junc >> ${mydir}/leafcutter/Juncfile/${date}_juncfile_$myexp.txt
    elif [[ " $(eval echo \${cond$j[*]}) " == *"$filename"* ]];then
      echo -e $filename'\t'$(eval echo \${cond[$j]}) >> ${mydir}/leafcutter/$myexp/Groups_$myexp.txt
      echo ${mydir}/leafcutter/Juncfile/$filename.junc >> ${mydir}/leafcutter/Juncfile/${date}_juncfile_$myexp.txt
    fi
  done
  
  # Step2

  source /home/volker/miniconda3/bin/activate activate Python27

  cd ${mydir}/leafcutter/$myexp/

  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo $mytime Leafcutter analysis Step2 $myexp
  
  # Check if Step2 is already done
  if [ -e "${mydir}/leafcutter/$myexp/${myexp}_perind_numers.counts.gz" ]; then
    echo "Step2 already completed"
  else
    python ${leafcutter_dir}/clustering/leafcutter_cluster_regtools.py -j ${mydir}/leafcutter/Juncfile/${date}_juncfile_$myexp.txt -m 50 -o $myexp -l 500000
  fi  

  conda deactivate

  # Step3
  
  source /home/volker/miniconda3/bin/activate activate leafcutter

  cd ${mydir}/leafcutter/$myexp/

  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo $mytime Leafcutter analysis Step3 $myexp

  # Check if reference exon file exists
  if [ -e "${leafcutter_dir}/reference/${reference}_exons_IDs.txt.gz" ]; then
    echo "Reference ${reference}_exons_IDs.txt.gz used for the leafcutter analysis"
  else
    ${leafcutter_dir}/scripts/gtf_to_exons.R /home/volker/reference/${reference}.annotation.gtf.gz ${leafcutter_dir}/reference/${reference}_exons_IDs.txt.gz
  fi
  
  # Check if Step 3 output files exist
  if [ -e "${mydir}/leafcutter/$myexp/${myexp}_cluster_significance.txt" ]; then
    echo "Step 3 leafcutter_ds Analysis already performed"
  else
    # Run the leafcutter R script
  ${leafcutter_dir}/scripts/leafcutter_ds.R --num_threads 10 --min_samples_per_intron=2 --min_samples_per_group=2 --output_prefix=${mydir}/leafcutter/$myexp/$myexp --exon_file=${leafcutter_dir}/reference/${reference}_exons_IDs.txt.gz ${mydir}/leafcutter/$myexp/${myexp}_perind_numers.counts.gz ${mydir}/leafcutter/$myexp/Groups_$myexp.txt
  fi  
  
  # Run leafviz
  
  ${leafcutter_dir}/leafviz/prepare_results.R \
  --meta_data_file ${mydir}/leafcutter/$myexp/Groups_$myexp.txt \
  -f 0.0001 \
  --code $myexp \
  -o ${mydir}/leafcutter/$myexp/${myexp}.RData \
  ${mydir}/leafcutter/$myexp/${myexp}_perind_numers.counts.gz \
  ${mydir}/leafcutter/$myexp/${myexp}_cluster_significance.txt \
  ${mydir}/leafcutter/$myexp/${myexp}_effect_sizes.txt \
  ${leafcutter_dir}/reference/${reference} 
  
  
  conda deactivate
  
  # Run the leafcutter combination R script to combine effect size and cluster significance and run analyses
  ${leafcutter_dir}/leafcutter_analysis_${reference}.R ${mydir} ${myexp}
  
  
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime All leafcutter steps for ${myexp} finished successfully"

done

# Organize conditions from array in a single string
cond_string=$(printf "%s," "${cond[@]}" | cut -d "," -f 1-${#cond[@]})

# Generate leafcutter-specific condition table for running the comparison script
> ${mydir}/leafcutter/ConditionTable.txt
echo -e "cond\tpath" >> ${mydir}/leafcutter/ConditionTable.txt
for i in "${cond[@]}"; do
	if [ ${i} != "control" ];then
		echo -e "${i}\t${mydir}/leafcutter/control_vs_${i}/control_vs_${i}_final.csv" >> ${mydir}/leafcutter/ConditionTable.txt
	fi
done	

# Run leafcutter Comparison script
${leafcutter_dir}/leafcutter_combined.R ${mydir} 

mytime=$(date "+%Y-%m-%d %H:%M:%S")
echo "$mytime All leafcutter steps finished successfully"
fi

######
####################################
#
# Perform combined analysis of DEG, DTU and AS
#
####################################
######

if  [ "$DESeq2_enable" = "true" ] && [ "$ISAR_enable" = "true" ] && [ "$leafcutter_enable" = "true" ]; then
	/home/volker/Tools/Combined_analysis_DEG_DTU_AS.R ${mydir}
fi

######
####################################
#
# Perform IRFinder analysis
#
####################################
######

if  [ "$IRFinder_enable" = "true" ]; then

echo ""
echo "###############################"
echo "## Perform IRFinder analysis ##"
echo "###############################"
echo ""
echo "IRFinder directory: $IRFinder_dir"

if [ ! -d "${mydir}/IRFinder/" ]; then
  mkdir ${mydir}/IRFinder/
fi

# Increase limit of simultaneously opened files
ulimit -n 20000
echo "ulimit is:"
ulimit -n


# Step1: Run IRFinder in FASTQ mode

mytime=$(date "+%Y-%m-%d %H:%M:%S")
echo $mytime IRFinder Step1 $myname

# Iterate over individual files
for i in "${samples[@]}"; do
mytime=$(date "+%Y-%m-%d %H:%M:%S")
echo $mytime Running IRFinder on ${i}

# Generate individual folder if neceassry
if [ ! -d "${mydir}/IRFinder/${i}/" ]; then
  mkdir ${mydir}/IRFinder/${i}/
fi
if [ -e "${mydir}/IRFinder/${i}/IRFinder-IR-dir.txt" ]; then
  echo "IRFinder analysis already performed for ${i}"
elif [ -e "${mydir}/IRFinder/${i}/IRFinder-IR-nondir.txt" ]; then
  echo "IRFinder analysis already performed for ${i}"
else

  # Check Sequencing design and designate variable accordingly
  if  [ "$seq_design" = "paired" ]; then
    seq_reads="$mydir/FASTQ/${i}_1.fq.gz $mydir/FASTQ/${i}_2.fq.gz"
  elif  [ "$seq_design" = "single" ]; then
    seq_reads="$mydir/FASTQ/${i}.fq.gz"
  fi

  # Run IRFinder using 20 threads
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Run IRFinder on ${i}"
  IRFinder -t 20 -M 60000 -r ${IRFinder_dir}/REF/$reference -d ${mydir}/IRFinder/${i} $seq_reads
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime IRFinder run on ${i} completed"

  # Remove output alignment files to safe disk space and as they are useless for the DESeq2 analysis
  rm ${mydir}/IRFinder/${i}/Unsorted.bam
fi
done

# Step2: Run differential intron retention analysis using DESeq2 -

mytime=$(date "+%Y-%m-%d %H:%M:%S")
echo $mytime IRFinder Step2 $myname

# Generate text files for each analysis
> ${mydir}/IRFinder/filePaths.txt
> ${mydir}/IRFinder/experiment.txt
echo -e "SampleNames"'\t'"Condition" >> ${mydir}/IRFinder/experiment.txt

# Iterate over bamfiles and fill essential information in text files
for j in $( eval echo {0..${arr_length}} );do
  for bamfile in `ls ${mydir}/BAM/*.bam`
  do

    # Get simplified filenames
    filename=$(basename -- "$bamfile")
    filename=${filename%%.*}

    # Output filePaths and experiment information to necessary textfiles if sample belongs to respective conditions
    if [[ " $(eval echo \${cond$j[*]}) " == *"$filename"* ]];then
      if [ -e "${mydir}/IRFinder/${i}/IRFinder-IR-dir.txt" ]; then
        echo -e $filename'\t'$(eval echo \${cond[$j]}) >> ${mydir}/IRFinder/experiment.txt
        echo ${mydir}/IRFinder/$filename/IRFinder-IR-dir.txt >> ${mydir}/IRFinder/filePaths.txt
      elif [ -e "${mydir}/IRFinder/${i}/IRFinder-IR-nondir.txt" ]; then
        echo -e $filename'\t'$(eval echo \${cond[$j]}) >> ${mydir}/IRFinder/experiment.txt
        echo ${mydir}/IRFinder/$filename/IRFinder-IR-nondir.txt >> ${mydir}/IRFinder/filePaths.txt
      fi
    fi
  done
done

# Organize conditions from array in a single string
cond_string=$(printf "%s," "${cond[@]}" | cut -d "," -f 1-${#cond[@]})

# Run R script for IRFinder analysis
mytime=$(date "+%Y-%m-%d %H:%M:%S")
echo "$mytime Run IRFinder main R script"

${IRFinder_dir}/IRFinder.R ${mydir} ${cond_string}

mytime=$(date "+%Y-%m-%d %H:%M:%S")
echo "$mytime IRFinder main R script completed"

# Run python script to get combined xlsx output
mytime=$(date "+%Y-%m-%d %H:%M:%S")
echo "$mytime Run IRFinder combination script"

python3 ${IRFinder_dir}/xlsx_combine_IRFinder.py ${mydir}/IRFinder/Combined ${mydir}/IRFinder/Combined/${myname}.IRFinder.out.final.xlsx ${cond_string}

mytime=$(date "+%Y-%m-%d %H:%M:%S")
echo $mytime All IRFinder steps finished successfully

####################

mytime=$(date "+%Y-%m-%d %H:%M:%S")
echo $mytime All analysis steps completed

fi

######
####################################
#
# Perform rRNA quantification analysis
#
####################################
######

if  [ "$rRNA_quant_enable" = "true" ]; then

echo ""
echo "##############################"
echo "## Perform rRNA quantification analysis ##"
echo "##############################"
echo ""

# Generate QC folder if necessary
if [ ! -d "${mydir}/QC/rRNA" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${mydir}/QC/rRNA"
  mkdir -p ${mydir}/QC/rRNA
fi


# Execute bowtie2 mapping to rRNA

  source /home/volker/miniconda3/bin/activate activate bowtie2
  

for i in "${samples[@]}"; do

	if [ -e "$mydir/QC/rRNA/${i}_bowtie2_rRNA.log" ]; then
	   	 	mytime=$(date "+%Y-%m-%d %H:%M:%S")
	   		echo "$mytime Sample ${i} already analyzed with bowtie2"
	   	else
	   		mytime=$(date "+%Y-%m-%d %H:%M:%S")
	   		echo "$mytime start rRNA quantification for ${i}"
	   		(bowtie2 -p 15 --no-unal --local -x $bowtie_rRNA_dir -1 ${mydir}/FASTQ/${i}_1.fq.gz -2 ${mydir}/FASTQ/${i}_2.fq.gz -S ${mydir}/FASTQ/${i}_forREMOVAL.bam)  2> ${mydir}/QC/rRNA/${i}_bowtie2_rRNA.log

		rm ${mydir}/FASTQ/${i}_forREMOVAL.bam
	fi

done

  conda deactivate
  
fi

######
####################################
#
# Perform multiQC analysis
#
####################################
######

if  [ "$QC_enable" = "true" ]; then

echo ""
echo "##############################"
echo "## Perform MultiQC analysis ##"
echo "##############################"
echo ""

# Generate QC folder if necessary
if [ ! -d "${mydir}/QC/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${mydir}/QC/"
  mkdir ${mydir}/QC/
fi

# Generate QC qualimap folder if necessary
if [ ! -d "${mydir}/QC/QC_qualimap/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${mydir}/QC/QC_qualimap/"
  mkdir ${mydir}/QC/QC_qualimap/
fi

###
# Run Qualimap 
###

# Unset DISPLAY to prevent X11 error
unset DISPLAY

source /home/volker/miniconda3/bin/activate activate qualimap

for i in "${samples[@]}"; do

  if [ -e "$mydir/BAM/${i}.bam" ]; then
	  
	   if [ ! -d "${mydir}/QC/QC_qualimap" ]; then
		  mkdir ${mydir}/QC/QC_qualimap/
	   fi

	    # Paired end reads
	   if  [ "$seq_design" = "paired" ]; then
	   	if [ -e "$mydir/QC/QC_qualimap/${i}/qualimapReport.html" ]; then
   	 	mytime=$(date "+%Y-%m-%d %H:%M:%S")
   		echo "$mytime Sample ${i} already analyzed with Qualimap"
 			else
		    	mytime=$(date "+%Y-%m-%d %H:%M:%S")
		    	echo "$mytime Paired-end sample ${i} Qualimap analysis"
			    /home/volker/Tools/qualimap_v2.3/qualimap \
			    rnaseq \
			    --java-mem-size=20G \
			    -bam $mydir/BAM/${i}.bam \
			    -gtf /home/volker/reference/Gencode/${reference}.annotation.gtf \
			    -p 'strand-specific-reverse' \
			    -pe \
			    -s \
			    -outdir $mydir/QC/QC_qualimap/${i}
	        fi

	  # Single end reads
	  elif  [ "$seq_design" = "single" ]; then
		if [ -e "$mydir/QC/QC_qualimap/${i}/qualimapReport.html" ]; then
	   	mytime=$(date "+%Y-%m-%d %H:%M:%S")
	   	echo "$mytime Sample ${i} already analyzed with Qualimap"
	 		else
			mytime=$(date "+%Y-%m-%d %H:%M:%S")
			echo "$mytime Single-end sample ${i} Qualimap analysis"
			    /home/volker/Tools/qualimap_v2.3/qualimap \
			    rnaseq \
			    --java-mem-size=20G \
			    -bam $mydir/BAM/${i}.bam \
			    -gtf /home/volker/reference/Gencode/${reference}.annotation.gtf \
			    -p 'strand-specific-reverse' \
			    -outdir $mydir/QC/QC_qualimap/${i}
		fi
	  fi
	  
  fi        
done

conda deactivate

# Execute multiQC

  source /home/volker/miniconda3/bin/activate activate multiqc
  
  multiqc --version

  multiqc ${mydir}/FASTQ ${mydir}/Salmon ${srvdir}/${myname} ${mydir}/QC/QC_qualimap -v -f -o ${mydir}/QC/ --sample-names ${mydir}/Samples.txt

# Check for Salmon -> perform Salmon specific multiQC for parsing in R
	if [ -d "$mydir/Salmon" ]; then

			multiqc ${mydir}/Salmon -f -o ${mydir}/QC/QC_salmon --sample-names ${mydir}/Samples.txt

	fi
	
# Check for Salmon_NMDRegHumanTxome -> perform Salmon specific multiQC for parsing in R
	if [ -d "$mydir/Salmon" ]; then

			multiqc ${mydir}/Salmon_NMD -f -o ${mydir}/QC/QC_salmon_NMD --sample-names ${mydir}/Samples.txt

	fi	
	
# Check for rRNA quantification -> perform rRNA specific multiQC for parsing in R
	if [ -d "$mydir/QC/rRNA" ]; then

			multiqc ${mydir}/QC/rRNA -f -o ${mydir}/QC/QC_rRNA --sample-names ${mydir}/Samples.txt

	fi

  conda deactivate
  
fi

fi
