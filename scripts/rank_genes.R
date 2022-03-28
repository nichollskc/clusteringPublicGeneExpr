library(dplyr)
library(tibble)
library(tidyr)

full_dataset = paste0(snakemake@wildcards[["dataset"]], snakemake@wildcards[["normalisation"]])
expression = read.csv(snakemake@input[["expression"]], check.names=FALSE, sep='\t', row.names=NULL)

probe_to_gene = read.csv(snakemake@input[["probe_to_gene"]], check.names=FALSE, sep='\t', row.names=NULL)
probe_values = inner_join(expression, probe_to_gene, by="probe_id")

gene_means = probe_values %>%
    group_by(gene) %>%
    summarise(across(where(is.numeric), mean)) %>%
    column_to_rownames("gene")

write.table(gene_means, snakemake@output[["gene_means"]], sep='\t', col.names=TRUE, row.names=TRUE)

gene_ranks = gene_means %>%
    mutate(across(everything(), function(x) rank(x, ties.method="average")))

write.table(gene_ranks, snakemake@output[["ranks"]], sep='\t', col.names=TRUE, row.names=TRUE)
