library(biomaRt)
library(dplyr)
library(ggplot2)
library(impute)
library(stringr)
library(tibble)
library(tidyr)
library(vsn)

build_sample_info <- function() {
    sample_info = read.csv("data/datasets/smith/sample_info.txt", sep='\t')
    sample_info$disease = sapply(sample_info$Factor.Value.DISEASESTATE., function(x) {
        result = ""
        if (x == "normal") result = "HC"
        if (x == "Microscopic Polyangiitis") result = "VMP"
        if (x == "Systemic lupus erythematosus") result = "SLE"
        if (x == "Wegener's Granulomatosis") result = "VWG"
        return (result)
    })
    sample_info$patient_num = sapply(sample_info$Source.Name, function(x) substring(x, 2))
    sample_info$extra = sapply(sample_info$Normalization.Name, function(x) str_split(x, ":")[[1]][2])
    sample_info$dataset = "SMITH"
    sample_info$patient_id = paste0(sample_info$disease, "_", sample_info$dataset, "_", sample_info$patient_num)
    sample_info$neat_sample_id = paste0(sample_info$patient_id, "_", sample_info$extra)
    sample_info$column_name = gsub(":", ".", sample_info$Normalization.Name)

    empty_pbmcs = read.table("data/datasets/smith/PBMC_norm.txt", sep='\t', row.names=1, header=TRUE, nrow=1)
    sample_info = sample_info[sample_info$column_name %in% colnames(empty_pbmcs), ]
    write.csv(sample_info, "data/datasets/smith/neat_sample_info.csv")
}

calculate_average_logfc <- function(gene_means, name) {
    hc_medians = filter(gene_means, disease == "HC") %>% summarise(across(starts_with("GENE_"), median))

    logfc = data.frame(mapply('-', gene_means %>% select(starts_with("GENE_")), hc_medians))
    logfc = cbind(logfc, select(gene_means, !starts_with("GENE_")))
    colnames(logfc) = colnames(gene_means)
    rownames(logfc) = rownames(gene_means)

    to_plot = pivot_longer(logfc, cols=starts_with("GENE_"), names_to="GENE", names_prefix="GENE_")
    g = ggplot(to_plot, aes(x=GENE, y=value, color=disease))  + geom_point() + facet_wrap("~disease")
    ggsave(paste0("plots/smith", name, ".png"), g)

    average_logfc_individuals = to_plot %>%
        select(neat_sample_id, GENE, value) %>%
        group_by(neat_sample_id) %>%
        summarise(mean_logfc = mean(value),
                  se_logfc = std(value),
                  .groups="keep") %>%
        arrange(neat_sample_id)

    return(average_logfc_individuals)
}

read_gene_list <- function(gene_list_file) {
    genes_table = read.table(gene_list_file, header=1, sep=',')
    return (genes_table$SYMBOL)
}

std <- function(x) sd(x)/sqrt(length(x))

calculate_ranked_data_smith <- function(genelist, processed_data, biomart_df, col_to_reporter, name) {
    dir_genelists = "~/rds/rds-cew54-basis/People/KATH/publicGeneExpr/"
    genelist_file = paste0(dir_genelists, "data/pathways/processed/", genelist, ".csv")
    gene_list <- read_gene_list(genelist_file)

    print(paste0("Total gene list length: ", length(gene_list)))
    restricted_biomart = biomart_df %>% filter(hgnc_symbol %in% gene_list, refseq_mrna != "")
    combined = merge(restricted_biomart, col_to_reporter, by.y="Reporter.Database.Entry.refseq.", by.x="refseq_mrna")

    restricted_processed_data = processed_data[rownames(processed_data) %in% combined$Reporter.Name, ]
    restricted_processed_data$Reporter.Name = rownames(restricted_processed_data)

    probe_values = merge(restricted_processed_data,
                         combined[, c("Gene.Symbol", "Reporter.Name")],
                         by="Reporter.Name")
    print(paste0("Genes also found in dataset: ", length(unique(probe_values$Gene.Symbol))))

    # Impute missing values
    values_only = probe_values %>% select(contains("PBMC")) %>% apply(FUN=as.numeric, MARGIN=1:2)
    values_only = data.frame(impute.knn(values_only)$data)
    values_only$Gene.Symbol = probe_values$Gene.Symbol

    gene_means = values_only %>%
        mutate(across(!starts_with("Gene.Symbol"), as.numeric)) %>%
        group_by(Gene.Symbol) %>%
        summarise(across(!starts_with("Gene.Symbol"), mean)) %>%
        column_to_rownames(var = "Gene.Symbol") %>%
        t() %>%
        data.frame()
    gene_names = paste0("GENE_", colnames(gene_means))
    colnames(gene_means) = gene_names
    gene_means$column_name = rownames(gene_means)
    small_sample_info = sample_info[!duplicated(sample_info$neat_sample_id),
                                    c("patient_id", "disease", "dataset", "column_name", "neat_sample_id")]
    gene_means = merge(gene_means, small_sample_info, by="column_name") %>%
        select(starts_with("GENE_"), everything())

    to_plot = pivot_longer(gene_means, cols=starts_with("GENE_"), names_to="GENE", names_prefix="GENE_")
    g = ggplot(to_plot, aes(x=GENE, y=value, color=disease))  + geom_point() + facet_wrap("~disease")
    ggsave(paste0("plots/smith_logtransformed_", name, ".png"), g)

    average_logfc_individuals = calculate_average_logfc(gene_means, name)
    write.csv(average_logfc_individuals, paste0("data/datasets/smith/", name, ".csv"))

    return(average_logfc_individuals)
}

apply_vsn_transformation <- function(processed_data, name) {
    g = meanSdPlot(as.matrix(processed_data))$gg
    ggsave(paste0("plots/smith_mean_variance_ranks_log_", name, ".png"), g)
    g = meanSdPlot(as.matrix(processed_data), ranks=FALSE)$gg
    ggsave(paste0("plots/smith_mean_variance_log_", name, ".png"), g)

    processed_data = processed_data %>% mutate(across(everything(), function(x) 2^x)) %>% as.matrix()

    vsn_smith = vsn2(as.matrix(processed_data))

    g = meanSdPlot(vsn_smith)$gg
    ggsave(paste0("plots/smith_mean_variance_vsn_ranks_", name, ".png"), g)
    g = meanSdPlot(vsn_smith, ranks=FALSE)$gg
    ggsave(paste0("plots/smith_mean_variance_vsn_", name, ".png"), g)

    df = data.frame(as.matrix(vsn_smith))
    colnames(df) = colnames(processed_data)
    rownames(df) = rownames(processed_data)
    return(df)
}

pbmcs = read.table("data/datasets/smith/PBMC_norm.txt", sep='\t', row.names=1, header=TRUE)
pbmcs = pbmcs[-c(1), ] %>% mutate(across(everything(), as.numeric))

#vsn_smith = apply_vsn_transformation(pbmcs, "")
#write.csv(vsn_smith, "data/datasets/smith/after_vsn.csv")
vsn_smith = read.csv("data/datasets/smith/after_vsn.csv", row.names=1)

#mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
# Around 100 refseqs have multiple genes...
#biomart_df = getBM(attributes=c( "hgnc_symbol", 'refseq_mrna'), mart=mart)
#biomart_df$Gene.Symbol = biomart_df$hgnc_symbol
#save(biomart_df, "data/datasets/smith/biomart.Rdata")
load("data/datasets/smith/biomart.Rdata")

col_to_reporter = read.csv("data/datasets/smith/A-MEXP-1152.adf.txt", skip=18, sep='\t')

sample_info = read.csv("data/datasets/smith/neat_sample_info.csv")

calculate_ranked_data_smith("ifn", pbmcs, biomart_df, col_to_reporter, name="ifn")
calculate_ranked_data_smith("ifn", vsn_smith, biomart_df, col_to_reporter, name="ifn_vsn")
