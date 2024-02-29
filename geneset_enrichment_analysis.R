# https://rpubs.com/jrgonzalezISGlobal/enrichment

# Load the necessary libraries

library(stringr)
library(clusterProfiler)
library(pathview)
library(stringr)
library(AnnotationHub)
library(ggridges)
library(enrichplot)

# extract annotation DB of sheep

ah <- AnnotationHub()
AnnotationHub::query(ah, c("Ovis", "aries"))
Oaries <- ah[["AH111978"]]
columns(Oaries) # see available identifiers that can be used in this package

# create geneList, which is a list of all genes from deseq2 with p<0.1 and lfc =0
# The variable res is the result from deseq2 and resSig contains all the DEGs

# The deseq2 results contains a mix of gene symbols and gene IDs (entrez IDs). For the enrichment analysis, it is advisable to use ensembl IDs or ENTREZ Ids. Here, we will be using Entrez IDs.
# So, the genes with gene symbols need to be converted to their entrezids. The entrez ids were retrieved from the gtf file which can be used to get entrez ids for our significant gene list. This is one way
# There is another way to retrieve the entrez ids using a function in the clusterProfiler package.

# 1. method 1 to retrieve entrez ids
entrez<-read.csv("/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0/EntrezIDsfromGtf",sep="\t", header=TRUE,row.names=NULL)#35105
# modify column names
colnames(entrez) <- colnames(entrez)[2:ncol(entrez)]
#drop last column
entrez <- entrez[1:(ncol(entrez)-1)]
# # Some of the GeneIDs and symbols contain "LOC" info, we need to remove it to merge it with the res and resSig files from deseq2.
entrez$GeneID <- stringr::str_remove(entrez$GeneID, "LOC")
entrez$GeneSymbol <- stringr::str_remove(entrez$GeneSymbol, "LOC")

res_allgenes<-res
res_allgenes$GeneID<-rownames(res)
res_degs<-resSig
res_degs$GeneID<-rownames(resSig)
res_degs$GeneSymbol<-rownames(resSig)

# removing "LOC" info
res_allgenes$GeneID <- stringr::str_remove(res_allgenes$GeneID, "LOC")
rownames(res_allgenes) = stringr::str_remove(rownames(res_allgenes), "LOC")
res_degs$GeneID <- stringr::str_remove(res_degs$GeneID, "LOC")
rownames(res_degs) = stringr::str_remove(rownames(res_degs), "LOC")

###### res_allgenes_entrez<-merge(as.data.frame(res_allgenes), entrez, by="GeneID") # This significantly reduced the number of mapping gene IDS from gtf file, 35,117 genes were reduced to 16,120
###### res_degs_entrez<-merge(as.data.frame(res_degs), entrez, by="GeneSymbol") # same here
# Therefore trying another approach to retrieve the ENTREZ IDs

# method 2 to retireve entrez ids
# get entrezIds for the genes with symbols
res_allgenesEntrez <- AnnotationDbi::select(Oaries, keys =  res_allgenes$GeneID,
  columns = c('ENTREZID'), keytype = 'SYMBOL')

res_sigGenesEntrez <- AnnotationDbi::select(Oaries, keys = res_degs$GeneID,
  columns = c('ENTREZID'), keytype = 'SYMBOL')

# Now replace the NA values in entrezid column with the values from first column. This is done because many values were not converted to ENTREZ IDs but already has it from GTF file. We will retain those.
res_allgenesEntrez$ENTREZID <- ifelse(is.na(res_allgenesEntrez$ENTREZID), res_allgenesEntrez$SYMBOL, res_allgenesEntrez$ENTREZID)
colnames(res_allgenesEntrez) = c("GeneID", "ENTREZID")
res_allgenes_with_entrez = merge(as.data.frame(res_allgenes), res_allgenesEntrez, by = "GeneID") #All entrezIDs have been retrieved for the DEGs (except for 2 genes)

res_sigGenesEntrez$ENTREZID <- ifelse(is.na(res_sigGenesEntrez$ENTREZID), res_sigGenesEntrez$SYMBOL, res_sigGenesEntrez$ENTREZID)
colnames(res_sigGenesEntrez) = c("GeneID", "ENTREZID")
res_degs_with_entrez = merge(as.data.frame(res_degs), res_sigGenesEntrez, by = "GeneID") #All entrezIDs have been retrieved for the DEGs (except for 2 genes)
# write.csv(res_degs_with_entrez,"6.deseq2/CH4production_DEGs_mapped_to_Entrezids.csv")

#####################################################
# GO AND KEGG ENRICHMENTS
#####################################################

# In order to asses functional enrichment, both DE gene list and gene universe must be annotated in Entrez IDs:

res_universe<-res_allgenes_with_entrez$ENTREZID
# omit any NA values 
res_universe<-na.omit(res_universe)
# sort the list in decreasing order (required for clusterProfiler)
res_universe = sort(res_universe, decreasing = TRUE)
res_sigGenes<-res_degs_with_entrez$ENTREZID
res_sigGenes = sort(res_sigGenes, decreasing = TRUE)

ans.go <- enrichGO(gene = res_sigGenes, ont = "ALL",
                   OrgDb =Oaries,
                   universe = res_universe,
                   readable=TRUE,
                   pvalueCutoff = 0.05)
tab.go <- as.data.frame(ans.go)
write.csv(tab.go,"8.geneset.enrichments/GO_enrichments.csv")

ans.kegg <- enrichKEGG(gene = res_sigGenes,
                       organism = 'oas',
                       universe = res_universe,
                       pvalueCutoff = 0.05)
tab.kegg <- as.data.frame(ans.kegg)
write.csv(tab.kegg,"8.geneset.enrichments/KEGG_enrichments.csv")

#####################################################
# VISUALIZATIONS OF GO AND KEGG ENRICHMENTS
#####################################################

pdf("8.geneset.enrichments/relaxed.nonscaled/GO-Barplot.pdf", height=14)
barplot(ans.go, showCategory=20)
dev.off()

pdf("8.geneset.enrichments/relaxed.nonscaled/upsetplot_kegg.pdf")
upsetplot(ans.kegg, showCategory=18)
dev.off()

pdf("8.geneset.enrichments/relaxed.nonscaled/emapplot_kegg.pdf", width=12)
emapplot(pairwise_termsim(ans.kegg))
dev.off()

# In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories,  cnetplot function will extract the complex association between genes and pathways.
pdf("8.geneset.enrichments/relaxed.nonscaled/cnetplot_kegg.pdf", width=12)
cnetplot(ans.kegg, categorySize="pvalue", foldChange=res_sigGenes)
dev.off()

#####################################################
# PATHWAY ANALYSIS - RETRIEVING PATHWAY IMAGES
#####################################################

#preparing tables for pathway analysis
de=data.frame(res_degs_with_entrez$ENTREZID,res_degs_with_entrez$log2FoldChange)
dl=data.frame(res_degs_with_entrez$GeneID,res_degs_with_entrez$log2FoldChange)
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

# Another method if kegg org annotations not working or kegg species invalid error occurs
# preparing korg for oas as the above code does not give full info of pathways (https://support.bioconductor.org/p/9146074/)
data(korg, package="pathview")
head(korg)
korg[korg[,3]=="oas",]
# Next create your own, single line, korg object.
korg <- cbind("ktax.id" = "T03117", "tax.id" = "9940", "kegg.code" = "oas",
               "scientific.name" = "Ovis aries", "common.name" = "sheep",
               "entrez.gnodes" = "1", "kegg.geneid" = "101112667", "ncbi.geneid" = "101112667",
               "ncbi.proteinid" = "XP_027833621", "uniprot" = "A0A6P7ER41")

pv.out <- pathview(gene.data =geneList, pathway.id = pathwayids,
                    species = "oas", out.suffix = "oas",
                    kegg.native = TRUE)
