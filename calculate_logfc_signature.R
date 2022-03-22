library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)

read_gene_list <- function(gene_list_file) {
    genes_table = read.table(gene_list_file, header=1, sep=',')
    return (genes_table$SYMBOL)
}

std <- function(x) sd(x)/sqrt(length(x))

calculate_average_logfc <- function(gene_means, full_dataset, sample_info, genelist_name) {
    hc_medians = gene_means %>%
        rownames_to_column("neat_sample_id") %>%
        inner_join(sample_info, by="neat_sample_id") %>%
        filter(disease == "HC") %>%
        select(where(is.numeric)) %>%
        summarise(across(everything(), median))

    logfc = data.frame(mapply('-', gene_means, hc_medians))
    colnames(logfc) = colnames(gene_means)
    rownames(logfc) = rownames(gene_means)
    write.table(logfc, snakemake@output[["logfc"]], sep='\t', col.names=TRUE, row.names=FALSE)

    logfc_pivoted = pivot_longer(logfc %>% rownames_to_column("neat_sample_id"), !neat_sample_id, names_to="gene")
    pivoted_logfc_with_info = inner_join(logfc_pivoted, sample_info, by='neat_sample_id')

    g = ggplot(pivoted_logfc_with_info, aes(x=disease, y=value, colour=disease)) + geom_boxplot() + facet_wrap("gene")
    ggsave(g, filename=paste0("plots/logfc_", full_dataset, "_", genelist_name, ".png"), width=20, height=20)

    average_logfc_individuals = pivoted_logfc_with_info %>%
        select(neat_sample_id, gene, value) %>%
        group_by(neat_sample_id) %>%
        summarise(mean_logfc = mean(value),
                  se_logfc = std(value),
                  .groups="keep") %>%
        arrange(neat_sample_id)

    return(average_logfc_individuals)
}

calculate_logfc_signature <- function(genelist_name, expression, full_dataset, sample_info) {
    dir_genelists = "~/rds/rds-cew54-basis/People/KATH/publicGeneExpr/"
    genelist_file = paste0(dir_genelists, "data/pathways/processed/", genelist_name, ".csv")
    genelist <- read_gene_list(genelist_file)

    probe_to_gene = read.csv(snakemake@input[["probe_to_gene"]], check.names=FALSE, sep='\t', row.names=NULL)
    print(paste0("Total gene list length: ", length(genelist)))
    probes_in_list <- probe_to_gene[probe_to_gene$gene %in% genelist, ]
    probe_values_in_list = inner_join(expression, probes_in_list, by="probe_id")
    print(paste0("Genes also found in dataset: ", length(unique(probe_values_in_list$gene))))

    gene_means = probe_values_in_list %>%
        group_by(gene) %>%
        summarise(across(where(is.numeric), mean)) %>%
        column_to_rownames("gene") %>%
        t() %>%
        as.data.frame()

    gene_means_pivoted = pivot_longer(gene_means %>% rownames_to_column("neat_sample_id"), !neat_sample_id, names_to="gene")
    pivoted_with_info = inner_join(gene_means_pivoted, sample_info, by='neat_sample_id')

    g = ggplot(pivoted_with_info, aes(x=disease, y=value, colour=disease)) + geom_boxplot() + facet_wrap("gene")
    ggsave(g, filename=paste0("plots/gene_means_", full_dataset, "_", genelist_name, ".png"), width=20, height=20)

    write.table(gene_means, snakemake@output[["gene_means"]], sep='\t', col.names=TRUE, row.names=FALSE)

    average_logfc_individuals = calculate_average_logfc(gene_means, full_dataset, sample_info, genelist_name)
    write.table(average_logfc_individuals, snakemake@output[["average_logfc"]], sep='\t', col.names=TRUE, row.names=FALSE)
}

full_dataset = paste0(snakemake@wildcards[["dataset"]], snakemake@wildcards[["extra"]])

expression = read.csv(snakemake@input[["expression"]], check.names=FALSE, sep='\t', row.names=NULL)
sample_info = read.csv(snakemake@input[["sample_info"]], check.names=FALSE, sep='\t', row.names=NULL)
calculate_logfc_signature(snakemake@wildcards[["genelist_name"]], expression, full_dataset, sample_info)
