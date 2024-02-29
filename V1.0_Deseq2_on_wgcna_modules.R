# extract genes from the inetersting modules
yg_pvr = geneInfo %>% filter(moduleColor %in% "yellowgreen" | moduleColor %in% "palevioletred3")
yellowgreen = geneInfo %>% filter(moduleColor %in% "yellowgreen")
palevioletred3 = geneInfo %>% filter(moduleColor %in% "palevioletred3")

# extract the raw counts for these genes in the modules
ygLfcRaw <- countData[rownames(countData) %in% yellowgreen$Genes,]
pvrLfcRaw <- countData[rownames(countData) %in% palevioletred3$Genes,]
yg_pvr_Raw = rbind(ygLfcRaw,pvrLfcRaw)

# removing "LOC" info
rownames(yg_pvr_Raw) = stringr::str_remove(rownames(yg_pvr_Raw), "LOC")

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
yg_pvr_Raw <- yg_pvr_Raw[,unique(rownames(metaData))]
all(colnames(yg_pvr_Raw) == rownames(metaData))

deseq2Data <- DESeqDataSetFromMatrix(countData=yg_pvr_Raw, colData=metaData, design= ~CH4production)

# Stringent approach where we keep only rows that have at least 10 reads total
keep <- rowSums(counts(deseq2Data)) >= 10
deseq2Data <- deseq2Data[keep,] #240 remaining

deseq2Data <- DESeq(deseq2Data)


#loop through results and extract significant DEGs for each model term
# speify the cut-offs for pval and lfc in the below variables.
# make sure to change the filenames with the cutoff values before saving the deg file (Line 60)

# STRINGENT - with reads removed in steps 32 to 34
# RELAXED - with reads not removed

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
  write.csv(resSig , file=paste0("7.wgcna/",results[i],".0.1P.0LFC.updownDEGs_STRINGENT.csv"), row.names = T)
  upresultstable[results[i],"upDEGs"] = upDEGs
  downresultstable[results[i],"downDEGs"] = downDEGs 
}

pval = 0.1
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
  write.csv(resSig , file=paste0("7.wgcna/",results[i],".0.1P.0.584LFC.updownDEGs_RELAXED.csv"), row.names = T)
  upresultstable[results[i],"upDEGs"] = upDEGs
  downresultstable[results[i],"downDEGs"] = downDEGs 
}

# Enrichment Analysis of the genes
# Taking the resSig variable from the above deseq2 (one with p<0.1 and lfc=0 stringent approach)
resSig$GeneID = rownames(resSig)
YgPvrEntrez <- AnnotationDbi::select(Oaries, keys =  rownames(resSig),
  columns = c('ENTREZID','GENENAME'), keytype = 'SYMBOL') # Oaries is the annotation db created from the geneset_enrichment_analysis.R code

# Now replace the NA values in entrezid column with the values from first column. This is done because many values were not converted to ENTREZ IDs but already has it from GTF file. We will retain those.
YgPvrEntrez$ENTREZID <- ifelse(is.na(YgPvrEntrez$ENTREZID), YgPvrEntrez$SYMBOL, YgPvrEntrez$ENTREZID)
colnames(YgPvrEntrez) = c("GeneID", "ENTREZID")
YgPvr_with_entrez = merge(as.data.frame(resSig), YgPvrEntrez, by = "GeneID") #All entrezIDs have been retrieved for the 36 genes with p<0.1 and lfc=0 threshold

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

# In order to asses functional enrichment, both DE gene list and gene universe must be annotated in Entrez IDs:

res_universe<-res_allgenes_with_entrez$ENTREZID
# omit any NA values 
res_universe<-na.omit(res_universe)
# sort the list in decreasing order (required for clusterProfiler)
res_universe = sort(res_universe, decreasing = TRUE)
res_sigGenes<-YgPvr_with_entrez$ENTREZID
res_sigGenes = sort(res_sigGenes, decreasing = TRUE)

ans.go <- enrichGO(gene = res_sigGenes, ont = "ALL",
                   OrgDb = Oaries,
                   universe = res_universe,
                   readable=TRUE,
                   pvalueCutoff = 0.05)
tab.go <- as.data.frame(ans.go)
write.csv(tab.go,"8.geneset.enrichments/WGCNA_Yg_Pvr_GO_enrichments.csv")

ans.kegg <- enrichKEGG(gene = res_sigGenes,
                       organism = 'oas',
                       universe = res_universe,
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
cnetplot(ans.kegg, categorySize="pvalue", foldChange=res_sigGenes)
dev.off()

#####################################################
# PATHWAY ANALYSIS - RETRIEVING PATHWAY IMAGES
#####################################################

#preparing tables for pathway analysis
de=data.frame(YgPvr_with_entrez$ENTREZID,YgPvr_with_entrez$log2FoldChange)
dl=data.frame(YgPvr_with_entrez$GeneID,YgPvr_with_entrez$log2FoldChange)
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
              limit      = list(gene=max(abs(geneList)), cpd=1))
