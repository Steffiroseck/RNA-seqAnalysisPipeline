library(AnnotationHub)
library(stringr)
library(dplyr)

hub <- AnnotationHub()
query(hub, c("ovis","orgdb"))
oaries <- hub[["AH111977"]]
length(keys(oaries))#34775
sum(is.na(mapIds(oaries, keys(z), "ENTREZID", "GID"))) #how many have no symbol? 268
select(oaries, head(keys(z)), "SYMBOL", "GID") # what do those symbols look like

oaries_info<- select(oaries, keys(oaries), c("SYMBOL","ALIAS","GENENAME","ENSEMBL"))

LOCID <- oaries_info[which(
  str_starts(oaries_info$SYMBOL, "LOC") == T),]

# we will use the hub genes identified from wgcna.R script
intra_modular_analysis.hubgene$SYMBOL<- rownames(intra_modular_analysis.hubgene)
# remove LOC in GeneID
#intra_modular_analysis.hubgene$GeneID<-gsub("LOC","",as.character(intra_modular_analysis.hubgene$GeneID))

merge <- merge(oaries_info,intra_modular_analysis.hubgene, by ="SYMBOL")
write.csv(merge,"7.wgcna/hubgenes.with.genenamesinfo.csv",row.names=F)

