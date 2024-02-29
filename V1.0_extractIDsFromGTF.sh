# extract the entrez ids from the gtf and gff files (downloaded from NCBI)
gffpath="/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0"
gfffile="/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0/genomic.gff"
gtffile="/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0/genomic.gtf"

# From the GTF file
cat $gtffile | awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"}' \
| sed 's/gene_id "//' | sed 's/gene_name "//' | sed 's/"//g' | sed 's/ //g' > $gffpath/GeneIDsfromGTF

sed 's/[db_xrefGeneID:,]//g' $gffpath/GeneIDsfromGTF > $gffpath/EntrezIDsfromGtf
sed -i '1i GeneSymbol\tGeneID' $gffpath/EntrezIDsfromGtf

# From the GFF file
zgrep 'GeneID' $gfffile \
  | cut -f9 | perl -pe 's/ID.*(GeneID:\d+).*gene=([^;]*).*/\1\t\2/g' \
  | sort -u > $gffpath/extracted_entrezids

# remove GeneID from all the lines and add column names
sed 's/GeneID://g' $gffpath/extracted_entrezids > $gffpath/entrezids 
sed  -i '1i GeneID\tGeneName' $gffpath/entrezids # 35106 total

# getting description
grep 'GeneID' $gfffile \
  | cut -f9 | perl -pe 's/ID.*(GeneID:\d+).*description=([^;]*).*/\1\t\2/g' \
  | sort -u > $gffpath/extracted_entrezids_with_description
  
# remove GeneID from all the lines and add column names
sed 's/GeneID://g' $gffpath/extracted_entrezids_with_description > $gffpath/entrezids_with_description
sed  -i '1i GeneID\tdescription' $gffpath/entrezids_with_description

rm extracted_entrezids_with_description
rm extracted_entrezids

# extract the gene ontology ids and gene ids from the .gaf file (downloaded from NCBI)
awk -F '\t| |"' -v OFS="\t" '{print $2,$5}' $gffpath/GCF_016772045.2-RS_2023_10_gene_ontology.gaf > $gffpath/entrezids_GOids #93113 IDs only have GO Ids

# combine the GOids and entrezIds file, because the GO ID file does not have all the entrezIDs of rambouillet genome
awk 'NR==FNR {h[$1] = $2; next} {print $1"\t"$2"\t"h[$1]}' $gffpath/entrezids_GOids $gffpath/entrezids > $gffpath/Full_entrezids_with_GOids

gtf2bed < genomic.gtf | cut -f4 > EntrzGeneIDs
