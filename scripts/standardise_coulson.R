library(stringr)
library(dplyr)
library(tidyr)
library(tibble)

extract_gene_name <- function(x) {
    strsplit(x, " // ")[[1]][2]
}

create_gene_rows <- function(row) {
    unique_gene_names = unique(sapply(strsplit(row['gene_assignment'], " /// ")[[1]],
                                      extract_gene_name))
    if (length(unique_gene_names) > 0) {
        df = data.frame(probe_id=row['probeset_id'],
                        gene=unique_gene_names,
                        row.names=NULL) 
    } else {
        df = data.frame()
    }
    return (df)
}

build_sample_info <- function() {
    sample_info = read.csv("data/datasets/coulson/sample_info.txt", sep='\t', check.names=FALSE)

    parts = sapply(sample_info[["Characteristics[individual]"]], function(x) {
           parts = strsplit(x, "_")[[1]]
           disease = parts[1]
           if (disease == "control") {
               disease = "HC"
           }
           subtype = disease
           if (parts[1] == "SLE" && parts[2] == "12") {
               patient_num = "122"
           } else if (parts[1] == "T1D" && length(parts) == 3) {
               patient_num = parts[3]
               subtype = paste0("T1D", parts[2])
           } else {
               patient_num = parts[2]
           }

           return (c(disease, patient_num, subtype))
    })

    parts_df = data.frame(t(parts))
    colnames(parts_df) = c("disease", "patient_num", "subtype")

    sample_info$disease = parts_df$disease
    sample_info$patient_num = parts_df$patient_num
    sample_info$subtype = parts_df$subtype
    sample_info$dataset = "COULSON"

    sample_info$neat_sample_id = paste0(sample_info$disease, "_", sample_info$dataset, "_", sample_info$patient_num)
    sample_info$patient_id = sample_info$neat_sample_id
    sample_info$column_name = sample_info[["Assay Name"]]
    sample_info = sample_info[, c("neat_sample_id", "patient_id", "subtype", "dataset", "disease", "column_name")]
    sample_info = distinct(sample_info, column_name, .keep_all=TRUE)

    sample_info = sample_info[sample_info$disease %in% c("SLE", "HC", "T1D"), ]

    return (sample_info)
}

build_expression <- function() {
    expression = read.csv("data/datasets/coulson/ifn_signature.data_matrix.tab", sep='\t', check.names=FALSE)

    # Trim first row
    expression = expression[-c(1), ]
    colnames(expression)[colnames(expression) == 'Hybridization REF'] <- 'probe_id'
    expression = expression %>% mutate(across(-c(probe_id), as.numeric))

    return (expression)
}

build_probe_to_gene <- function() {
    affy = read.csv("data/datasets/coulson/HuGene-1_1-st-v1.na36.hg19.probeset.csv", sep=',', skip=22)
    all_gene_names = do.call(rbind, apply(affy , MARGIN=1, FUN=create_gene_rows)) 
    probe_to_gene = na.omit(all_gene_names) 

    return (probe_to_gene)
}

sample_info = build_sample_info()
expression = build_expression()
probe_to_gene = build_probe_to_gene()

restricted = expression[, sample_info$column_name]
colnames(restricted) = sample_info$neat_sample_id
expression = cbind(probe_id = expression$probe_id,
                   restricted)

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

write.table(sample_info, "data/datasets/coulson/sample_info.tsv", sep='\t', row.names=FALSE, col.names=TRUE)
write.table(expression, "data/datasets/coulson/expression.tsv", sep='\t', row.names=FALSE, col.names=TRUE)
write.table(probe_to_gene, "data/datasets/coulson/probe_to_gene.tsv", sep='\t', row.names=FALSE, col.names=TRUE)
