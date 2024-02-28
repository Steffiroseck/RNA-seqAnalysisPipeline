#!/bin/bash

SECONDS=0

# Software Requirements
# FastQC v0.12.1 
# cutadapt version 1.15
# TrimmomaticPE
# HISAT2 version 2.1.0
# STAR
# featureCounts v2.0.1
# StringTie 2.2.1

# Path to the working directory where fastq files are present. In this example, the samples are present in a sub-directory of sheep named 24_samples. Hence, the current working directory is set.
cwd = "/home/afbi-roses@nigov.net/Sheep/"
cd $cwd

# Create directories to store outputs and give full permission to those folders for the output files to be saved.
echo "Creating directories to store outputs!"
mkdir 1.fatsqc
mkdir 2.trimmomatic
mkdir 3.fastqc.after.trimmomatic
mkdir 4.hisat2
mkdir 4.star
mkdir 5.featurecounts
mkdir logs

# Step 1: Quality Control
# Run FASTQC (2 minutes per sample)

echo "Starting FastQC"
for file in 24_Samples/*R1_001.fastq.gz; do
    prefix="${file%R1_001.fastq,gz}"
    reverse="${file%R1_001.fastq.gz}R2_001.fastq.gz"
    fastqc -1 $file -2 $reverse -t 15 -o 1.fastqc
done
echo "FastQC finished succesfully!"

# Step 2: Adapter trimming
# Run Trimmomatic

echo "Starting Trimmomatic"
threads="8" #higher the number faster the speed; but make sure you have enough CPU support
for R1 in 24_Samples/*_R1_001.fastq.gz ; do
    R2="${R1%_R1_001.fastq.gz}_R2_001.fastq.gz"
    sample=$(echo $R1|sed 's/_R1_001.fastq.gz//'|sed 's/24_Samples\///'); 
    TrimmomaticPE -threads $threads "$R1" "$R2" 2.trimmomatic/"${sample}_R1_paired.fastq.gz" 2.trimmomatic/"${sample}_R1_unpaired.fastq.gz"  \
    2.trimmomatic/"${sample}_R2_paired.fastq.gz" 2.trimmomatic/"${sample}_R2_unpaired.fastq.gz" ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 
done
echo "Trimmomatic run finished succesfully!"

# Step3: Fastqc on adapter trimmed reads

echo "Starting FastQC on trimmed reads"
for file in 2.trimmomatic/*_R1_paired.fastq.gz; do
    prefix="${file%_R1_paired.fastq.gz}"
    reverse="${file%_R1_paired.fastq.gz}_R2_paired.fastq.gz"
    fastqc -1 $file -2 $reverse -t 15 -o 3.fastqc.after.trimmomatic
done
echo "FastQC on trimmed reads finished succesfully!"

# Step 4: Alignment to the reference genome using HISAT2
# 4.1: Build reference genome index. The genome should be in unzipped format or else it wont work

echo "Starting to build the genome index files"
refgenome="Reference_geneome/ARS-UI_Ramb_v3.0/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna"
hisat2-build  $refgenome $refgenome

# 4.2 Now Align to the reference genome

echo "Index files generated. Now performing the alignment. This might take a while.."
for R1 in 2.trimmomatic/*_R1_paired.fastq.gz; do
    R2=$(echo $R1| sed 's/_R1_/_R2_/'); 
    sample=$(echo $R1|sed 's/_R1_paired.fastq.gz//'|sed 's/2.trimmomatic\///'); 
    hisat2 -q --time --novel-splicesite-outfile 4.hisat2/$sample.tsv --summary-file 4.hisat2/$sample.summary.txt \
    --met-file 4.hisat2/$sample.met.txt --threads $threads -x Reference_geneome/ARS-UI_Ramb_v3.0/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna \
    -1 $R1 -2 $R2 | tee >(samtools flagstat - > 4.hisat2/$sample.flagstat) \
    | samtools sort -O BAM | tee 4.hisat2/$sample.bam \
    | samtools index - 4.hisat2/$sample.bam.bai &> 4.hisat2/$sample.hisat2Log.txt;
done
echo "Alignment finished succesfully!"

# Step 5: Alignment using STAR
# 5.1 Build genome index files. You need to create a directory to store the index files or if you wish you can write it in the cwd

fastafile="/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0/GCF_016772045.2_ARS-UI_Ramb_v3.0_genomic.fna"
gtffile="/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0/genomic.gtf"
threads=8

STAR --runThreadN $threads --runMode genomeGenerate --genomeDir StarIndex/ --genomeFastaFiles $fastafile --sjdbGTFfile $gtffile --sjdbOverhang 100

# 5.2 After indexing, we go for STAR alignment.The input files of STAR can be single-end or pair-end fastq files. 
# If the annotation file (.gtf file) is provided, the accuracy of alignment can be increased.

indexDir='/mnt/sda1/RNA/40-815970407/Sheep/StarIndex/'
threads=8

for R1 in 2.trimmomatic/*_R1_paired.fastq.gz; do
    R2=$(echo $R1| sed 's/_R1_/_R2_/'); 
    sample=$(echo $R1|sed 's/_R1_paired.fastq.gz//'|sed 's/2.trimmomatic\///'); 
    STAR --runThreadN $threads --genomeDir $indexDir  --readFilesIn $R1 $R2 \
    --outFileNamePrefix 4.star/$sample  --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c --limitBAMsortRAM 1756793772
done

# Step 6 Generate read counts matrix using featureCounts
# When you want to analyze the data for differential gene expression analysis, it would be convenient to have counts for all samples in a single file (gene count matrix).

echo "Generating the read counts matrix.."
gtffile="/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0/genomic.gtf"

featureCounts -T 8 -t 'gene' -g 'gene_id' -f -a $gtffile -o 5.featurecounts/Lambs.featurecounts.hisat2 4.hisat2/*.bam
featureCounts -T 8 -t 'gene' -g 'gene_id' -f -a $gtffile -o 5.featurecounts/Lambs.featurecounts.star 4.star/*.bam

## Since the featureCounts output has additional columns with information about genomic coordinates, gene length etc., 
## we can use the cut command to select only those columns that you are interested in. Columns 1 and sample wise counts columns

cut -f1,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32 5.featurecounts/Lambs.featurecounts.hisat2 > 5.featurecounts/Lambs.featurecounts.hisat2.Rmatrix
cut -f1,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32 5.featurecounts/Lambs.featurecounts.star > 5.featurecounts/Lambs.featurecounts.star.Rmatrix
echo "Read count matrix generated!"
echo "Have a look at the counts matrix. You may change the column names if needed.."

# After these step, try to remove extra texts in the column names (eg..bam in the sample names)

# Step 6 Multiqc Report generation
echo "Generating Multiqc report.."
multiqc -o 9.multiqc .
echo "Report generated succesfully!"

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
