library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)

read_gene_list <- function(gene_list_file) {
    genes_table = read.table(gene_list_file, header=1, sep=',')
    return (genes_table$SYMBOL)
}

std <- function(x) sd(x)/sqrt(length(x))

calculate_average_ranks <- function(genelist_name, ranks, full_dataset, sample_info) {
    dir_genelists = "~/rds/rds-cew54-basis/People/KATH/publicGeneExpr/"
    genelist_file = paste0(dir_genelists, "data/pathways/processed/", genelist_name, ".csv")
    genelist <- read_gene_list(genelist_file)

    print(paste0("Total gene list length: ", length(genelist)))
    ranks_in_signature = ranks[rownames(ranks) %in% genelist, ]
    print(paste0("Genes also found in dataset: ", nrow(ranks_in_signature)))

    ranks_pivoted = pivot_longer(ranks_in_signature %>% rownames_to_column("gene"), !gene, names_to="neat_sample_id")
    print(head(ranks_pivoted))
    print(head(sample_info))
    pivoted_ranks_with_info = inner_join(ranks_pivoted, sample_info, by='neat_sample_id')
    print(head(pivoted_ranks_with_info))
    pivoted_ranks_with_info$value = pivoted_ranks_with_info$value / nrow(ranks)

    print(head(pivoted_ranks_with_info))

    g = ggplot(pivoted_ranks_with_info, aes(x=disease, y=value, colour=disease)) + geom_boxplot() + facet_wrap("gene")
    ggsave(g, filename=paste0("plots/ranks_", full_dataset, "_", genelist_name, ".png"), width=20, height=20)

    write.table(pivoted_ranks_with_info, snakemake@output[["ranks"]], sep='\t', col.names=TRUE, row.names=FALSE)

    average_ranks_individuals = pivoted_ranks_with_info %>%
        select(neat_sample_id, gene, value) %>%
        group_by(neat_sample_id) %>%
        summarise(mean_ranks = mean(value),
                  se_ranks = std(value),
                  .groups="keep") %>%
        arrange(neat_sample_id)

    write.table(average_ranks_individuals, snakemake@output[["mean_ranks"]], sep='\t', col.names=TRUE, row.names=FALSE)
}

full_dataset = paste0(snakemake@wildcards[["dataset"]], snakemake@wildcards[["normalisation"]])

sample_info = read.csv(snakemake@input[["sample_info"]], check.names=FALSE, sep='\t', row.names=NULL)
ranks = read.csv(snakemake@input[["ranks"]], check.names=FALSE, sep='\t', row.names=1)
calculate_average_ranks(snakemake@wildcards[["genelist_name"]], ranks, full_dataset, sample_info)
