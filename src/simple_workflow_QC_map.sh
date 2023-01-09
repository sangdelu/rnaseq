#!/usr/bin/env bash
# RNA-Seq analysis workflow
#
###################################
# Important metadata about raw sequencing data, that you should 
# ask data-providers:
# 1. Sequencing platform/technique.
# 2. Sequencing direction/strand.
###################################

# Necessary tools/packages:
# FastQC, MultiQC
# STAR
# R-DESeq2
#
# R-ggplot2 or other plot packages
# Trimmomatic or other trimmers
###################################

#### Data security ####
RawDataDir="/.../"
# Before start, make sure the raw data is UNWRITABLE #
# DO NEVER overwrite the raw data #

# Other file paths
LinkDataDir="/.../"
TempDir="/.../"
TrimmDir=$LinkDataDir"trimmed/"
OutputDir="/.../"
STARIndex="/ref/.../Ensembl/STAR_idx"
DeSeq2Code="/.../Deseq2.R"
STARTransCode="/.../STAR_transpose.py"
SpliceJunctionGTF="/ref/.../*.gtf"
CountMatrix=$OutputDir"countmatrix.csv"

# If you have these directories already created, unable these three lines.
mkdir $LinkDataDir
mkdir $TempDir
mkdir $OutputDir
##########################################

# 0. Create file links and change their names:
###### You need to name your links as following: ####
# A file name should be:
# Single strand: SampleName_SampleNumber_L*.fastq.gz
# Paired-end: SampleName_SampleNumber_<1 or 2>.fastq.gz <- Illumina standard
# Format file name while linking. If one sample was split, merge them into one.
# Think twice about raw data path before start. Modify this loop according to how your data was grouped.
for x in $(ls $RawDataDir)
do
    ln -s $RawDataDir$x $LinkDataDir$x
done

# 1. Quatlity control
# Depending on sequencing technique, you might need to modify the raw data 
# (e.g. delete few bases in tail/head of each sequence).

# FastQC
mkdir $LinkDataDir"FastQC"

for x in $(ls $LinkDataDir*fq.gz)
do
    fastqc -o $LinkDataDir"FastQC" --dir $TempDir --thread 16 --noextract $x
done

# MultiQC
multiqc --dirs $LinkDataDir"FastQC" --filename "RawSeq_MultiQC_res" --outdir $LinkDataDir
# Check point 1 - Read mQC report directly.
# Trimming if necessary.
 mkdir trimmed
for x in *_1.fq.gz
do
    bn=$(basename -s _1.fq.gz $x)
    java -jar /opt/bioinf/Trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 $x $bn"_2.fq.gz" "trimmed/"$bn"_trim_1.fq.gz" /dev/null "trimmed/"$bn"_trim_2.fq.gz" /dev/null TRAILING:20 SLIDINGWINDOW:4:20 ILLUMINACLIP:/opt/bioinf/Trimmomatic/Trimmomatic-0.38/adapters/NexteraPE-PE.fa:2:30:10:1:true MINLEN:40
done

# Re-do fastQC after trimming, and check

for x in $(ls $TrimmDir*fq.gz)
do
    fastqc -o $TrimmDir"FastQC" --dir $TempDir --thread 16 --noextract $x
done

# MultiQC
multiqc --dirs $TrimmDir"FastQC" --filename "Trimmed_MultiQC_res" --outdir $TrimmDir


# 2. Alignment & Count
# STAR
# Don't forget to check STAR's version.
STAR --version
mkdir $OutputDir"STAR"
# Make sure the version of .GTF file is the newest.
# If you have a STAR index, which was built from GTF file,
# you could load it first with '-genomeLoad' and then start loop. It saves time.
# Input Single-end:readA,readB,readC
# Input paired-end: readA1,readB1,readC1 readA2,readB2,readC2
# A file name should be: SampleName_SampleNumber_L*_<1 or 2>.fastq.gz <- Illumina standard
for fq in $(ls $TrimmDir*_1.fq.gz)
do
    bn=$(basename -s _1.fq.gz $fq)
    mkdir $OutputDir"STAR/"$bn"/"
    # Modify this line according to your data type.
    STAR --genomeDir $STARIndex --readFilesIn $TrimmDir$bn"_1.fq.gz" $TrimmDir$bn"_2.fq.gz" --sjdbGTFfile $SpliceJunctionGTF --runThreadN 16 --readFilesCommand zcat --outFileNamePrefix $OutputDir"STAR/"$bn"/" --outSAMmode None --quantMode GeneCounts
    # Set '-outSAMmode None' when you do not need .sam file.
done

# Check point 2:
# Read the report from multiqc
multiqc --module star --dirs $OutputDir"STAR/" -dd 2 --filename "STAR_MultiQC_res.log" --outdir $OutputDir"STAR/"

# 3. Read-counts data transpose, Matrix generation
# Automatic generation with Python 3
# Hint: choose unstrand count automatically. If you need other strands, please modify the code.
python3 STAR_transpose.py $OutputDir"STAR" $CountMatrix
# The matrix should have sample names as column names, gene names as row names.
