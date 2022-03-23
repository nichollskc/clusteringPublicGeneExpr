library(biomaRt)

mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
# Around 100 refseqs have multiple genes...
biomart_df = getBM(attributes=c( "hgnc_symbol", 'refseq_mrna'), mart=mart)
biomart_df$Gene.Symbol = biomart_df$hgnc_symbol
save(biomart_df, "data/datasets/biomart.Rdata")
