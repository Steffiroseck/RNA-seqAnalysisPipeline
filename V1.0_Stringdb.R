# Instal library stringdb
#install.packages("https://www.bioconductor.org/packages/release/bioc/src/contrib/STRINGdb_2.14.0.tar.gz", repos = NULL, type="source")
# Load the libraries
 library(STRINGdb)


system("mkdir 9.NetworkAnalysis")

# load string database
#9940 for sheep
string_db <- STRINGdb$new(version = "11.5", species = 9940,score_threshold=400, protocol="https")
class(string_db)

# map the gene names to the STRING database identifiers using the map method. The map function adds an additional column with STRING identifiers to the dataframe that is passed as first parameter.
data<-read.csv("6.deseq2/CH4production_DEGs_mapped_to_Entrezids.csv",sep=",",row.names=NULL)
rownames(data)<-data$ENTREZID
data$gene<-rownames(data)
data_mapped <- string_db$map(data, "gene", removeUnmappedRows = TRUE )
dim(data_mapped)
hits <- data_mapped$STRING_id
string_db$plot_network(hits)

# PAYLOAD MECHANISM
# filter by p-value and add a color column
# (i.e. green down-regulated gened and red for up-regulated genes)
data_mapped_sig <- string_db$add_diff_exp_color(subset(data_mapped, log10(padj) >= -log10(0.01) | abs(log2FoldChange) >= 0),
                                                        logFcColStr="log2FoldChange" )
head(data_mapped_sig)
table(data_mapped_sig$color)

# post payload information to the STRING server
payload_id <- string_db$post_payload( data_mapped_sig$STRING_id,
                                      colors=data_mapped_sig$color )

# display a STRING network png with the "halo"
pdf("9.NetworkAnalysis/String_network.pdf")
string_db$plot_network(hits, payload_id=payload_id)
dev.off()

#clustering
clustersList <- string_db$get_clusters(data_mapped$STRING_id)

# plot first 4 clusters
par(mfrow=c(2,2))
for(i in seq(1:4)){
   string_db$plot_network(clustersList[[i]])
   }

# perform GO and KEGG pathway enrichment analysis
#GO enrichment analysis
enrichmentGO <- string_db$get_enrichment(hits, category = "Process")
head(enrichmentGO)

#KEGG enrichment analysis
enrichmentKEGG <- string_db$get_enrichment( hits, category = "KEGG" )
head(enrichmentKEGG)
