#!/usr/bin/env bash
# RNA-Seq analysis workflow

source ./rna_seq.config 


###########################################################
# Not really necessary to change following scripts.       #
###########################################################
mkdir $LinkDataDir
mkdir $TempDir
mkdir $OutputDir
mkdir $OutputDir"qc"

# 0. Create file links and change their names:
for x in $(find $RawDataDir -iname *$FileSuffix)
do
    bn=$(basename $x)
    ln -s $x $LinkDataDir$bn
done

# 1. Quatlity control of raw data
# FastQC
mkdir $LinkDataDir"FastQC"

for x in $(ls $LinkDataDir*$FileSuffix)
do
    fastqc -o $LinkDataDir"FastQC" --dir $TempDir --thread 16 --noextract $x
done

mkdir $OutputDir"qc/pre"
# MultiQC
multiqc --dirs $LinkDataDir"FastQC" --filename "RawSeq_MultiQC_res" --outdir $OutputDir"qc/pre/"
echo -n "Please find the multiQC report of raw data in "$OutputDir"qc/pre ."
#echo -n "Is it necessary to trimm raw data?[yes/no]"
#read totrimm
#while [ $totrimm != "yes" -o $totrimm != "no" ];
#do
#    echo -n "Please type in yes or no."
#    read totrimm
#done

if  [ $totrimm == "yes" ]; then
    echo -n "Data will be trimmed."
    mkdir trimmed
    for x in *$FileSuffix_1
    do
        bn=$(basename -s $FileSuffix_1 $x)
        java -jar /opt/bioinf/Trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 $x $bnFileSuffix_2 "trimmed/"$bn"_trim"$FileSuffix_1 /dev/null "trimmed/"$bn"_trim"$FileSuffix_2 /dev/null TRAILING:20 SLIDINGWINDOW:4:20 ILLUMINACLIP:$adapter_file:2:30:10:1:true MINLEN:40
    done

# Re-do fastQC after trimming, and check
    for x in $(ls $TrimmDir*$FileSuffix)
    do
        fastqc -o $TrimmDir"FastQC" --dir $TempDir --thread 16 --noextract $x
    done

# MultiQC
    mkdir $OutputDir"qc/post"
    multiqc --dirs $TrimmDir"FastQC" --filename "Trimmed_MultiQC_res" --outdir $OutputDir"qc/post/"
    echo -n "Please find the multiQC report of trimmed data in "$OutputDir"qc/post ."
elif [ $totrimm == "no" ]; then
    echo -n "Data won't be trimmed."
    TrimmDir=$LinkDataDir
fi


# 2. Alignment & Count
# STAR
echo "The version of STAR is: "
$STARPath --version
mkdir $OutputDir"align"
# Input Single-end:readA,readB,readC
# Input paired-end: readA1,readB1,readC1 readA2,readB2,readC2
echo -n "Start mapping with STAR..."

if  [ $SequencingMode == "paired" ]; then
    for fq in $(ls $TrimmDir*$FileSuffix_1)
    do
        bn=$(basename -s $FileSuffix_1 $fq)
        mkdir $OutputDir"align/"$bn"/"
        # Modify this line according to your data type.
        $STARPath --genomeDir $STARIndex --readFilesIn $TrimmDir$bn$FileSuffix_1 $TrimmDir$bn$FileSuffix_2 --sjdbGTFfile $SpliceJunctionGTF --runThreadN 16 --readFilesCommand zcat --outFileNamePrefix $OutputDir"align/"$bn"/" --outSAMmode None --quantMode GeneCounts
    done
elif [ $SequencingMode == "single" ]; then
    for fq in $(ls $TrimmDir*$FileSuffix)
    do
        bn=$(basename -s $FileSuffix $fq)
        mkdir $OutputDir"align/"$bn"/"
        # Modify this line according to your data type.
        $STARPath --genomeDir $STARIndex --readFilesIn $TrimmDir$bn$FileSuffix --sjdbGTFfile $SpliceJunctionGTF --runThreadN 16 --readFilesCommand zcat --outFileNamePrefix $OutputDir"align/"$bn"/" --outSAMmode None --quantMode GeneCounts
    done
fi

multiqc --module star --dirs $OutputDir"align/" -dd 2 --filename "STAR_MultiQC_res.log" --outdir $OutputDir"align/"

# Check point 2:
grep "Uniquely mapped reads %" $OutputDir"align/"*"/Log.final.out" > $OutputDir"align/uniquely.report"
echo "Please check the number of uniquely mapped reads in "$OutputDir"align/uniquely.report"
echo ""

# 3. Read-counts data transpose, Matrix generation
# Automatic generation with Python 3
echo -n "Generating count table..."
python3 $STARTransCode $OutputDir"align" $CountMatrix $Strandness

