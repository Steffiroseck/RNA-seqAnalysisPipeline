# This R script does a differential expression analysis on the two important modules with highest correlations to methane production.
# yellow green and palevioletred3 modules

library(DESeq2)
library(stringr)
library(clusterProfiler)
library(pathview)
library(stringr)
library(AnnotationHub)
library(ggridges)
library(enrichplot)

# set the working directory
setwd("/mnt/sda1/00_fastq/Sheep/")
# create a directory to store the results
system("mkdir 8.wgcna.enrichments")

# Load the count data and metadata
countData<-read.csv("5.featurecounts/Lambs.featurecounts.hisat2.Rmatrix",sep="\t", header=T, check.names=F)
# run the below step if you want to remove any of the samples with poor mapping rates, as including this might induce noise in the deseq2 results
countData<-countData[ , !names(countData) %in% c("7085","7073")]# Remove 7085 and 7073 as they had poor mapping rates

# Remove the .bam, control, low, medium and high from the column names
colnames(countData)<-gsub(".bam","",colnames(countData))
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

# extract genes from the inetersting modules
yellowgreen = geneInfo %>% filter(moduleColor %in% "yellowgreen")

# extract the raw counts for these genes in the modules
ygLfcRaw <- countData[rownames(countData) %in% yellowgreen$Genes,]

# removing "LOC" info
rownames(ygLfcRaw) = stringr::str_remove(rownames(ygLfcRaw), "LOC")

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
ygLfcRaw <- ygLfcRaw[,unique(rownames(metaData))]
all(colnames(ygLfcRaw) == rownames(metaData))

# Scale the variables if its needed
metaData$CH4production = scale(metaData$CH4production, center=TRUE)
deseq2Data <- DESeqDataSetFromMatrix(countData=ygLfcRaw, colData=metaData, design= ~CH4production)

deseq2Data <- DESeq(deseq2Data)

#loop through results and extract significant DEGs for each model term
# speify the cut-offs for pval and lfc in the below variables.
# make sure to change the filenames with the cutoff values before saving the deg file (Line 60)

pval = 0.1
lfc = 0
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
  #write.csv(resSig , file=paste0("7.wgcna/deseq2.YG.PVR/",results[i],".0.05P.0LFC.updownDEGs.csv"), row.names = T)
  upresultstable[results[i],"upDEGs"] = upDEGs
  downresultstable[results[i],"downDEGs"] = downDEGs 
}

########################################
# Enrichment Analysis of the genes
########################################

# extract annotation DB of sheep

ah <- AnnotationHub()
AnnotationHub::query(ah, c("Ovis", "aries"))
Oaries <- ah[["AH111978"]]
columns(Oaries)

# Taking the res variable from the above deseq2 as it has all the genes with corresponding pvalue and log2fc values
res$GeneID = rownames(res)
YgEntrez <- AnnotationDbi::select(Oaries, keys =  rownames(res),
  columns = c('ENTREZID','GENENAME'), keytype = 'SYMBOL') # Oaries is the annotation db created from the geneset_enrichment_analysis.R code

# Now replace the NA values in entrezid column with the values from first column. This is done because many values were not converted to ENTREZ IDs but already has it from GTF file. We will retain those.
YgEntrez$ENTREZID <- ifelse(is.na(YgEntrez$ENTREZID), YgEntrez$SYMBOL, YgEntrez$ENTREZID)
colnames(YgEntrez) = c("GeneID", "ENTREZID")
Yg_with_entrez = merge(as.data.frame(res), YgEntrez, by = "GeneID") #All entrezIDs have been retrieved for the 36 genes with p<0.1 and lfc=0 threshold

# Extract entrezIDs for the res
res_allgenes<-res
res_allgenes$GeneID<-rownames(res)

# removing "LOC" info
res_allgenes$GeneID <- stringr::str_remove(res_allgenes$GeneID, "LOC")
rownames(res_allgenes) = stringr::str_remove(rownames(res_allgenes), "LOC")

res_allgenesEntrez <- AnnotationDbi::select(Oaries, keys =  res_allgenes$GeneID,
  columns = c('ENTREZID','GENENAME'), keytype = 'SYMBOL')

# Now replace the NA values in entrezid column with the values from first column. This is done because many values were not converted to ENTREZ IDs but already has it from GTF file. We will retain those.
res_allgenesEntrez$ENTREZID <- ifelse(is.na(res_allgenesEntrez$ENTREZID), res_allgenesEntrez$SYMBOL, res_allgenesEntrez$ENTREZID)
colnames(res_allgenesEntrez) = c("GeneID", "ENTREZID")
res_allgenes_with_entrez = merge(as.data.frame(res_allgenes), res_allgenesEntrez, by = "GeneID") #All entrezIDs have been retrieved for the DEGs

#####################################################
# GO AND KEGG ENRICHMENTS
#####################################################

# In order to asses functional enrichment, DE gene list must be annotated in Entrez IDs:
# sort the list in decreasing order (required for clusterProfiler)
res = sort(res, decreasing = TRUE)
res_Genes<-Yg_with_entrez$ENTREZID
res_Genes = sort(res_Genes, decreasing = TRUE)

ans.go <- enrichGO(gene = res_Genes, ont = "ALL",
                   OrgDb = Oaries,
                   readable=TRUE,
                   pvalueCutoff = 0.05)
tab.go <- as.data.frame(ans.go)
write.csv(tab.go,"8.geneset.enrichments/WGCNA_Yg_Pvr_GO_enrichments.csv")

ans.kegg <- enrichKEGG(gene = res_Genes,
                       organism = 'oas',
                       pvalueCutoff = 0.05)
tab.kegg <- as.data.frame(ans.kegg)
write.csv(tab.kegg,"8.geneset.enrichments/WGCNA_Yg_Pvr_KEGG_enrichments.csv")

#####################################################
# VISUALIZATIONS OF GO AND KEGG ENRICHMENTS
#####################################################

pdf("8.geneset.enrichments/WGCNA_Yg_Pvr_GO-Barplot.pdf", height=14)
barplot(ans.go, showCategory=20)
dev.off()

pdf("8.geneset.enrichments/WGCNA_Yg_Pvr_upsetplot_GO.pdf")
upsetplot(ans.go, showCategory=100)
dev.off()

pdf("8.geneset.enrichments/WGCNA_Yg_Pvr_upsetplot_kegg.pdf")
upsetplot(ans.kegg, showCategory=30)
dev.off()

pdf("8.geneset.enrichments/WGCNA_Yg_Pvr_emapplot_kegg.pdf", width=12)
emapplot(pairwise_termsim(ans.kegg))
dev.off()

# In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories,  cnetplot function will extract the complex association between genes and pathways.
pdf("8.geneset.enrichments/WGCNA_Yg_Pvr_cnetplot_kegg.pdf", width=12)
cnetplot(ans.kegg, categorySize="pvalue", foldChange=res_Genes)
dev.off()

#####################################################
# PATHWAY ANALYSIS - RETRIEVING PATHWAY IMAGES
#####################################################

#preparing tables for pathway analysis
de=data.frame(Yg_with_entrez$ENTREZID,Yg_with_entrez$log2FoldChange)
dl=data.frame(Yg_with_entrez$GeneID,Yg_with_entrez$log2FoldChange)
head(de)
head(dl)

geneList = de[,2]
names(geneList) = as.character(de[,1])
geneList = sort(geneList, decreasing = TRUE)
head(geneList)#geneList has entrezids and FC values
gene <- names(geneList)

# Retrieveing the pathway images
pathwayids=ans.kegg$ID   #make a vector of Pathway ids
keggspecies="oas"

x <- pathview(gene.data  = geneList,
              pathway.id = pathwayids,
              species    = keggspecies,
              gene.idtype = "KEGG",
              limit      = list(gene=max(abs(geneList)), cpd=1),
             #kegg.dir="8.geneset.enrichments/wgcna_YG_PVR_pathview")

# The .pathview images will be generated in the current directory, whereas .xml and original kegg images will be in the 8.geneset.enrichments/wgcna_YG_PVR_pathview folder.
