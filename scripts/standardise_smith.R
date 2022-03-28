library(dplyr)
library(impute)
library(stringr)
library(tibble)
library(tidyr)

source("scripts/utils.R")

build_sample_info <- function() {
    sample_info = read.csv("data/datasets/smith/sample_info.txt", sep='\t', check.names=FALSE)
    sample_info$disease = sapply(sample_info[["Factor Value[DISEASESTATE]"]], function(x) {
        result = ""
        if (x == "normal") result = "HC"
        if (x == "Microscopic Polyangiitis") result = "VMP"
        if (x == "Systemic lupus erythematosus") result = "SLE"
        if (x == "Wegener's Granulomatosis") result = "VWG"
        return (result)
    })
    sample_info$patient_num = sapply(sample_info[["Normalization Name"]], function(x) str_sub(str_split(x, ":")[[1]][1], 2, -1))
    sample_info$extra = sapply(sample_info[["Normalization Name"]], function(x) str_split(x, ":")[[1]][2])
    sample_info$dataset = "SMITH"
    sample_info$patient_id = paste0(sample_info$disease, "_", sample_info$dataset, "_", sample_info$patient_num)
    sample_info$neat_sample_id = paste0(sample_info$patient_id, "_", sample_info$extra)
    sample_info$column_name = sample_info[["Normalization Name"]]
    sample_info = sample_info[, c("neat_sample_id", "patient_id", "dataset", "extra", "disease", "column_name")]
    sample_info = distinct(sample_info, column_name, .keep_all=TRUE)

    return (sample_info)
}

build_expression <- function() {
    expression = read.csv("data/datasets/smith/PBMC_norm.txt", sep='\t', check.names=FALSE)
    # Trim first row
    expression = expression[-c(1), ]
    colnames(expression)[colnames(expression) == 'Normalization REF'] <- 'probe_id'
    expression = expression %>% mutate(across(-c(probe_id), as.numeric))

    # Remove rows that are entirely empty
    expression = expression[rowSums(is.na(expression)) != ncol(expression) - 1, ]
    return (expression)
}

impute_missing <- function(expression) {
    imputed = data.frame(impute.knn(as.matrix(expression[, -c(1)]))$data)
    imputed = cbind(probe_id = expression$probe_id,
                    imputed)
    rownames(imputed) = NULL

    return (imputed)
}

build_probe_to_gene <- function() {
    load("data/datasets/biomart.Rdata")
    biomart_df = biomart_df[, c("refseq_mrna", "Gene.Symbol")]
    # Remove rows containing blanks
    biomart_df = biomart_df[!apply(biomart_df == "", MARGIN=1, any), ]

    col_to_reporter = read.csv("data/datasets/smith/A-MEXP-1152.adf.txt", skip=18, sep='\t')
    col_to_reporter = col_to_reporter[, c("Reporter.Name", "Reporter.Database.Entry.refseq.")]
    colnames(col_to_reporter) = c("probe_id", "refseq_mrna")
    col_to_reporter = col_to_reporter[col_to_reporter$refseq_mrna != "", ]

    probe_to_gene = inner_join(col_to_reporter, biomart_df)
    probe_to_gene = probe_to_gene[, c("probe_id", "Gene.Symbol")]
    colnames(probe_to_gene) = c("probe_id", "gene")

    return (probe_to_gene)
}

expression = build_expression()
sample_info = build_sample_info()
# Use neat_sample_id as column names
colnames(expression)[-c(1)] <- sapply(colnames(expression)[-c(1)],
                                      function(x) unique(sample_info[sample_info$column_name == x, ]$neat_sample_id))

print(head(expression))
print(head(sample_info))
# Restrict samples to those in expression df
sample_info = restrict_to_present_samples(sample_info, expression)

# No longer need this column - was only used to match up expression and sample_info
sample_info$column_name = NULL

# Smith dataset has missing values
expression = impute_missing(expression)

probe_to_gene = build_probe_to_gene()

print(head(expression))
print(head(sample_info))
print(head(probe_to_gene))

print("Labels should match up:")
print(all(sample_info$neat_sample_id ==
          colnames(expression)[-c(1)]))
print(all(length(sample_info$neat_sample_id) ==
          length(colnames(expression)[-c(1)])))

print("neat_sample_id should have unique values")
print(length(unique(sample_info$neat_sample_id)) == length(sample_info$neat_sample_id))

write.table(sample_info, "data/datasets/smith/sample_info.tsv", sep='\t', row.names=FALSE, col.names=TRUE)
write.table(expression, "data/datasets/smith/expression.tsv", sep='\t', row.names=FALSE, col.names=TRUE)
write.table(probe_to_gene, "data/datasets/smith/probe_to_gene.tsv", sep='\t', row.names=FALSE, col.names=TRUE)
