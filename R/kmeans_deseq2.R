# Load libraries
library(DESeq2)
library(dplyr)
library(pheatmap)
library(tidyverse)
library(ggbeeswarm)
library(genefilter)

#set the working directory
setwd("/mnt/sda1/RNA/40-815970407/Sheep")

# Load the metadata file
metaData <-read.csv("metadata_with_methaneinfoadded_metadata.csv",sep=",",header=T)
kmeans_group <- kmeans(metaData$CH4production, centers = 2)
metaData$kmeansgroup <- kmeans_group$cluster
metaData<-metaData %>% mutate(kmeansgroup = dplyr::case_when(kmeansgroup== 1 ~ "High", kmeansgroup== 2 ~ "Low"))

# removing unwanted info
metaData$ID<-gsub("Control","",metaData$ID)
metaData$ID<-gsub("Low","",metaData$ID)
metaData$ID<-gsub("Medium","",metaData$ID)
metaData$ID<-gsub("High","",metaData$ID)
metaData$kmeansconcat<-paste(metaData$kmeansgroup,metaData$ID,sep="")
rownames(metaData)<-metaData$kmeansconcat
metaData$kmeansconcat<-NULL
metaData %>% 
  write_csv("6.deseq2/kmeans_metadata_with_methaneinfo_added.csv")

# Based on this kmeans grouping, the counts data need to be transformed in such a way that sample names should include information of high or low methane production info.
# example of counts data should look like:
# Geneid Low7050 Low7066 Low6914 Low6976 Low6968 Low6768 Low7220 Low7070
# LOC114110831       2       0       8       0       0       0       2       2
# LOC114112203       0       0       6       3       2       0       3       0
# LOC114110836       0       0       0       0       0       0       0       0
# LOC101102048       0       0       3       0       5       4       0       1
# LOC114113923      33      36      20      46      43      24      22      16
# LOC114112459      31      15      49      35      30      32      48      16

# Read the counts data
countData<-read.csv("6.deseq2/kmeans_LambAllSamples.featureCounts.RmatrixNew.csv")
colnames(countData)<-gsub(".bam","",colnames(countData))
rownames(countData) <- countData$Geneid
countData$Geneid<-NULL

# check library size (total mapped reads)
colSums(countData[,2:ncol(countData)])

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
countData <- countData[,unique(rownames(metaData))]
all(colnames(countData) == rownames(metaData))

deseq2Data <- DESeqDataSetFromMatrix(countData=countData, colData=metaData, design= ~kmeansgroup)

deseq2Data <- deseq2Data[rowSums(assay(deseq2Data)) > 0, ]
keep <- rowSums(counts(deseq2Data)) >= 4
deseq2Data <- deseq2Data[keep,]
dim(deseq2Data)#23871 genes

deseq2Data <- DESeq(deseq2Data)

results = resultsNames(deseq2Data)
upresultstable = matrix(nrow = length(results), ncol = 1, dimnames = list(results,"upDEGs"))
downresultstable = matrix(nrow = length(results), ncol = 1, dimnames = list(results,"downDEGs"))

for(i in 1:length(results)){

  res = results(deseq2Data, 
                name = results[i])
  resorder <- res[order(res$padj),]
  upDEGs = (length(na.omit(which(res$padj<0.1 & res$log2FoldChange > 0.584))))
  downDEGs = (length(na.omit(which(res$padj<0.1 & res$log2FoldChange < -0.584))))
  resSig = subset(resorder, padj < 0.1 & log2FoldChange > 0.584 | padj < 0.1 & log2FoldChange < -0.584)
  write.csv(resSig , file=paste0("6.deseq2/",results[i],".0.1p.lfc1.updownDEGs.csv"), row.names = T)
  upresultstable[results[i],"upDEGs"] = upDEGs
  downresultstable[results[i],"downDEGs"] = downDEGs 
}

# Plot gene counts for the significant genes

gene_set <- c("LOC114114465", "FAM3B", "NME4", "ATP6V0A4")
names(gene_set) <- gene_set
df <- lapply(gene_set, \(x) {
    y <- plotCounts(deseq2Data, x, c("kmeansgroup"), normalized=TRUE,returnData=TRUE)
    y$feature <- x
    return(y)
})

df <- do.call(rbind, df)
pdf("6.deseq2/kmeansgroup-Genecounts_of_sig_genes.pdf")
ggplot(df, aes(x=kmeansgroup, y=count, color=kmeansgroup)) +
  geom_jitter(
    position=position_jitterdodge(dodge.width=0.75),
    size=0.75) +
  facet_grid(feature~.) +
  scale_y_continuous(limits = c(0, 600))
  dev.off()

# Plot gene counts individually
pdf("6.deseq2/FAM3Bgene.pdf")
geneCounts <- plotCounts(deseq2Data, gene = "FAM3B", intgroup = c("kmeansgroup"),
                         returnData = TRUE, normalized=TRUE)
ggplot(geneCounts, aes(x = kmeansgroup, y = count, color = kmeansgroup)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()

pdf("6.deseq2/LOC114114465gene.pdf")
geneCounts <- plotCounts(deseq2Data, gene = "LOC114114465", intgroup = c("kmeansgroup"),
                         returnData = TRUE, normalized=TRUE)
ggplot(geneCounts, aes(x = kmeansgroup, y = count, color = kmeansgroup)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()

pdf("6.deseq2/NME4gene.pdf")
geneCounts <- plotCounts(deseq2Data, gene = "NME4", intgroup = c("kmeansgroup"),
                         returnData = TRUE, normalized=TRUE)
ggplot(geneCounts, aes(x = kmeansgroup, y = count, color = kmeansgroup)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()

pdf("6.deseq2/ATP6V0A4gene.pdf")
geneCounts <- plotCounts(deseq2Data, gene = "ATP6V0A4", intgroup = c("kmeansgroup"),
                         returnData = TRUE, normalized=TRUE)
ggplot(geneCounts, aes(x = kmeansgroup, y = count, color = kmeansgroup)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
dev.off()

# For comparisons of gene counts between raw and normalized dataset

#Heatmaps are a great way to look at gene counts. To do that, we can use a function in the pheatmap package.
#Next, we can select a subset of genes to plot. Although we could plot all  genes, let’s choose the 20 genes with the largest positive log2fold change.

vsd <- varianceStabilizingTransformation(deseq2Data)

pdf("6.deseq2/kmeans_pheatmap_top20genes_with_highest_variance.pdf")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("kmeansgroup")])
rownames(anno)<-colnames(mat)
pheatmap(mat, annotation_col = anno)
dev.off()
