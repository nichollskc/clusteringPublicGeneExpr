logfc_summaries = lapply(snakemake@input, function(x) read.csv(x, sep='\t'))
print(sapply(logfc_summaries, colnames))

shared_genes = Reduce(intersect, sapply(logfc_summaries, colnames))
print("Genes in each dataset:")
print(sapply(logfc_summaries, ncol))
print(paste0("Shared genes: ", length(shared_genes)))

write.csv(data.frame(SYMBOL=shared_genes), snakemake@output[[1]], row.names=FALSE)
