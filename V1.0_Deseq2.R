# This analysis takes the counts data generated from featureCounts and the metadata with the phenotype information as the input files.
# So make sure you have both files in the correct directories. In this scenario, I have copied the metadata to the main working directory (/mnt/sda1/RNA/40-815970407/Sheep)

# Load all the necessary libraries. If not installed, install using either:
install.packages("package") or
# to install biconductor packages, 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install(c("clusterProfiler", "AnnotationHub"))

library(DESeq2)
library(dplyr)
library(tidyverse)
library(EnhancedVolcano)

# set the working directory
setwd("/mnt/sda1/RNA/40-815970407/Sheep")

# create a directory to store DESeq2 outputs. Since I am pre-filtering reads in coming steps, 
# I am creating directories to keep stringent results (where pre-filtering is done), and relaxed ( where no pre-filtering applied).
system("mkdir 6.deseq2")
#system("mkdir 6.deseq2/stringent")
#system("mkdir 6.deseq2/relaxed")

# Read the raw count data generated from featureCounts
countData<-read.csv("5.featurecounts/Lambs.featurecounts.hisat2.Rmatrix",sep="\t", header=T, check.names=F)
# run the below step if you want to remove any of the samples with poor mapping rates, as including this might induce noise in the deseq2 results
countData<-countData[ , !names(countData) %in% c("7085","7073")]# Remove 7085 and 7073 as they had poor mapping rates

# Remove the .bam, control, low, medium and high from the column names
colnames(countData)<-gsub("Control","",colnames(countData))
colnames(countData)<-gsub("Low","",colnames(countData))
colnames(countData)<-gsub("Medium","",colnames(countData))
colnames(countData)<-gsub("High","",colnames(countData))

orig_names <- names(countData) # keep a back-up copy of the original names
geneID <- countData$Geneid# Convert count data to a matrix of appropriate form that DEseq2 can read
countData <- as.matrix(countData[ , -1]) #removing first column geneID from the table
# make sure the rownames are gene ids and first column onwards should be samples. any other columns should be removed.otherwise deseq2 will throw error.
sampleIndex <- colnames(countData)
countData <- as.matrix(countData[,sampleIndex])
rownames(countData) <- geneID
head(countData)

# check library size (total mapped reads)
colSums(countData[,2:ncol(countData)])

# Read the metadata file
metaData <-read.csv("metadata_with_methaneinfoadded_metadata.csv",sep=",",header=T)
#metaData1 = metaData[!grepl('Low|Medium', metaData$ID),]# To remove any rows with Low or medium in the ID column. (Only for test case)
#metaData = metaData1
dim(metaData)
head(metaData)
rownames(metaData) <- metaData$ID
metaData$ID <- factor(metaData$ID)
rownames(metaData)<-gsub("[a-zA-Z ]", "", rownames(metaData))
head(metaData)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
countData <- countData[,unique(rownames(metaData))]
all(colnames(countData) == rownames(metaData))

# Users should center and scale numeric variables in the design to improve GLM convergence if using numeric variables (recommended by deseq2 vignette)
# metaData$CH4production=scale(metaData$CH4production, center=TRUE) 
deseq2Data <- DESeqDataSetFromMatrix(countData=countData, colData=metaData, design= ~CH4production)

# Stringent approach where we keep only rows that have at least 10 reads total
keep <- rowSums(counts(deseq2Data)) >= 10
deseq2Data <- deseq2Data[keep,] #22,750 remaining

deseq2Data <- DESeq(deseq2Data)


# loop through results and extract significant DEGs for each model term
# specify the cut-offs for pval and lfc in the below variables.
# make sure to change the filenames with the cutoff values before saving the deg file (Line 91)

pval = 0.05
lfc = 0.584
results = resultsNames(deseq2Data)
upresultstable = matrix(nrow = length(results), ncol = 1, dimnames = list(results,"upDEGs"))
downresultstable = matrix(nrow = length(results), ncol = 1, dimnames = list(results,"downDEGs"))

for(i in 1:length(results)){

  res = results(deseq2Data, 
                name = results[i])
  resorder <- res[order(res$padj),]
  upDEGs = (length(na.omit(which(res$padj<pval & res$log2FoldChange > lfc))))
  downDEGs = (length(na.omit(which(res$padj<pval & res$log2FoldChange < -lfc))))
  resSig = subset(resorder, padj < pval & log2FoldChange > lfc | padj < pval & log2FoldChange < -lfc)
  write.csv(resSig , file=paste0("6.deseq2/",results[i],".0.05P.0LFC.updownDEGs.csv"), row.names = T)
  upresultstable[results[i],"upDEGs"] = upDEGs
  downresultstable[results[i],"downDEGs"] = downDEGs 
}

# fold change = 0.5, Log2FC= -1 (2 fold decrease)
# fc = 1.4, log2fc = 0.5 (1.5 times higher expression)
# log2fc =0 means no change

# these next R scripts are for a variety of visualization, QC and other plots to
# get a sense of what the RNAseq data looks like based on DESEq2 analysis
#
# 1) MA plot
# 2) PCA plot
# 3) Volcano plot
#
 
# 1. MA plot
pdf("6.deseq2/MAplot.pdf")
plotMA(res,main = "MA plot", alpha=0.1,colNonSig = "black", colSig = "red",colLine = "lightblue")
dev.off()

# 2. Volcano plot
# change the pvalue and fccutoffs as you wish (lines 124 and 125)
pdf("6.deseq2/EnhancedVolcano.pdf", width=14, height=9)
EnhancedVolcano(res,
	lab=rownames(res),
	x='log2FoldChange',
	y='padj',
	title='Volcano plot',
	caption = 'FC cutoff, 0.584; p-value cutoff, 0.05',
	pCutoff=0.05,
	FCcutoff=0.584,
	pointSize = 3.0)
dev.off()

##labSize=2,legendPosition = 'top',legendLabSize = 12 for volcanoplot if need to add gene names

# Extract the rawcounts for these significant genes (for WGCNA)
#required_df <- countData[rownames(countData) %in% rownames(resSig1),]
#write.table(as.data.frame(required_df), '6.deseq2/Control.vs.High.sig.genes.raw.counts.csv',quote=F, row.names=TRUE)
