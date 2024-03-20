# RNASequencingAnalysisPipeline
RNA-Seq analysis pipeline for sheep (Ovis aries) samples. The RNA-Sequencing on 24 samples were performed on Illumina Platform and reads were generated in fastq format. Once the fastq data is downloaded, then the below pipeline have been used for generating the biological information.

# Pre-requisites
Make sure you have installed the necessary softwares required to run the pipeline. These are mentioned in detail below.

**1. FastQC:** To generate quality reports for the raw fastq files coming from the sequencing platforms.
  ```
  conda install bioconda::fastqc
  ```
   or
   ```
   sudo apt-install fastqc
  ```
   or install from source
   ```
  wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
  unzip fastqc_v0.12.1.zip
  cd FastQC/
  chmod 755 fastqc
  ./fastqc --help

  ```
   
**2. Trimmomatic:** 
  ```
  conda install bioconda::trimmomatic
  ```
  or install from source. Refer to the link below.
  ```
  http://www.usadellab.org/cms/?page=trimmomatic
  ```

**3. Hisat2:** 
  ```
  sudo apt-get update
  sudo apt-get -y install hisat2
  ```
  or you can simply download from the source
  ```
  curl -s https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download
  unzip hisat2-2.2.1-Linux_x86_64.zip
  cd hisat2-2.2.1
  ./hisat2 -h
  ```
  **5. featureCounts:**
  The installation instructions can be found at : https://subread.sourceforge.net/featureCounts.html

After you have installed the necessary software, please proceed to running the actual pipeline. We would run them in the following order:
1. bash V1.0.RNASeqPipeline.sh or ./V1.0.RNASeqPipeline.sh
2. Rscript V1.0_Deseq2.R
3. Rscript V1.0_GO_KEGG_enrichment_analysis.R
4. Rscript V1.0_wgcna.R
5. Rscript V1.0_Stringdb.R
