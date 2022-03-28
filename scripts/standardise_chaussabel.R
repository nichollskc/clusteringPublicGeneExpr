library(dplyr)
library(stringr)
library(tibble)
library(tidyr)

create_gene_rows <- function(row) {
    unique_gene_names = unique(strsplit(row['Gene.Symbol'], " /// ")[[1]])
    if (length(unique_gene_names) > 0) {
        df = data.frame(probe_id=row['Probe.Set.ID'],
                        gene=unique_gene_names,
                        row.names=NULL) 
    } else {
        df = data.frame()
    }
    return (df)
}

# Function to allow named capture groups
re.capture = function(pattern, string, ...) {
    rex = list(src=string,
               result=regexpr(pattern, string, perl=TRUE, ...),
               names=list())

    for (.name in attr(rex$result, 'capture.name')) {
        rex$names[[.name]] = substr(rex$src,
                                    attr(rex$result, 'capture.start')[,.name],
                                    attr(rex$result, 'capture.start')[,.name]
                                    + attr(rex$result, 'capture.length')[,.name]
                                    - 1)
    }

    return(rex)
}

try_extract_sample_info <- function(sample_id, pattern, disease) {
    match_result = re.capture(pattern, sample_id)
    if (match_result$result[1] != -1) {
        patient_num = match_result$names$pat_num
        extra_info = match_result$names$extra
        AorB = match_result$names$AB
        result = c(sample_id,
                   disease,
                   patient_num,
                   extra_info,
                   AorB)
    } else {
        result = NULL
    }
    return (result)
}

extract_sample_id <- function(sample_id) {
    this_sample_info = NULL

    if (is.null(this_sample_info)) {
        this_sample_info = try_extract_sample_info(sample_id, "^SLE(?<pat_num>[0-9]+)_?(?<extra>.*?)_?(?<AB>A|B|)$", "SLE")
        if (! is.null(this_sample_info)) {
            this_sample_info[4] = gsub("/", ".", this_sample_info[4])
        }
    }
    if (is.null(this_sample_info)) {
        this_sample_info = try_extract_sample_info(sample_id, "^Lupus_(?<pat_num>[0-9]+|pre)(?<extra>)(?<AB>A|B)", "SLE")
    }
    if (is.null(this_sample_info)) {
        this_sample_info = try_extract_sample_info(sample_id, "^sys(?<pat_num>[0-9]+)(?<AB>[AB])[ \\.]?(?<extra>[0-9]{6}|)", "JIA")
    }
    if (is.null(this_sample_info)) {
        this_sample_info = try_extract_sample_info(sample_id, "^RET[ \\.]?(?<pat_num>[0-9]+)(?<extra>Froz|)(?<AB>A|B|)", "RET")
    }
    if (is.null(this_sample_info)) {
        this_sample_info = try_extract_sample_info(sample_id, "^(?i)Mel[ \\.](?<pat_num>[0-9]+)[ \\.](?<extra>pre).*U133[ \\.]?(?<AB>A|B)", "MEL")
    }
    if (is.null(this_sample_info)) {
        this_sample_info = try_extract_sample_info(sample_id, "^(?i)INF(?<pat_num>[0-9]+)(?<AB>A|B|).*Ecoli(?<extra>)", "InfECOLI")
    }
    if (is.null(this_sample_info)) {
        this_sample_info = try_extract_sample_info(sample_id, "^(?i)INF(?<pat_num>[0-9]+)(?<AB>A|B|).*STAPH(?<extra>)", "InfSTAPH")
    }
    if (is.null(this_sample_info)) {
        this_sample_info = try_extract_sample_info(sample_id, "^(?i)IDDM[\\. ]?(?<pat_num>[0-9]*)(?<extra>)[ \\.]?(?<AB>A|B)", "T1D")
    }
    if (is.null(this_sample_info)) {
        this_sample_info = try_extract_sample_info(sample_id, "^X?003.*ACC[\\.-](?<pat_num>[0-9]+)(?<extra>F|)[\\._](?:U133)?(?<AB>B|)", "MEL")
    }
    if (is.null(this_sample_info)) {
        this_sample_info = try_extract_sample_info(sample_id, "^(?i)(?:PBMCS__|)H[ \\.-]?(?<pat_num>[0-9]+|[A-Z]{2})_?(?<extra>[0-9]{6}|)(?:Healthy_U133)?(?<AB>B|)", "HC")
    }

    if (is.null(this_sample_info)) {
        this_sample_info = c(sample_id, "", "", "", "")
    }
    return(this_sample_info)
}

add_detailed_sample_info <- function(sample_df, full_ids) {
    res = lapply(full_ids, extract_sample_id)
    df = data.frame(do.call(rbind, res))
    colnames(df) = c("old_sample_id", "disease", "patient_num", "extra", "AorB")
    df$dataset = "CHA"
    df$patient_id = paste0(df$disease, "_", df$dataset, "_", df$patient_num)
    df$neat_sample_id = paste0(df$patient_id, "_", df$extra, "_", df$AorB)

    return(cbind(sample_df, df))
}

build_sample_info <- function() {
    sample_info = read.csv("data/datasets/chaussabel/sample_info.txt", sep='\t')
    sample_info = add_detailed_sample_info(sample_info, sample_info[, "Scan.Name"])
    sample_info$dataset = sapply(sample_info$Derived.Array.Data.Matrix.File, function(x) {
                                 if (x == "E-GEOD-11907-processed-data-1673830054.txt") {
                                     dataset = "CHA"
                                 } else if (x == "E-GEOD-11907-processed-data-1673830055.txt") {
                                     dataset = "CHB"
                                 } else {
                                     dataset = "CH"
                                 }
                                 return (dataset)
                   })
    sample_info$neat_sample_id = paste0(sample_info$disease, "_", sample_info$dataset, "_", sample_info$patient_num, "_", sample_info$extra)
    sample_info$column_name = sample_info$old_sample_id
    sample_info = sample_info[, c("neat_sample_id", "patient_id", "dataset", "extra", "disease", "column_name")]

    return (sample_info)
}

build_expression <- function(filename) {
    expression = read.csv(filename, sep='\t', row.names=1, check.names=FALSE)
    expression = expression[, c(FALSE, TRUE)]
    expression = expression[-c(1), ]
    expression = expression %>% mutate(across(everything(), as.numeric))
    expression = expression %>% mutate(across(everything(), log2))

    expression = cbind(probe_id = rownames(expression),
                       expression)
    rownames(expression) = NULL

    return (expression)
}

build_sample_expression_subset <- function(sample_info, filename) {
    rownames(sample_info) = sample_info$column_name

    expression_subset = build_expression(filename)
    sample_info_subset = sample_info[colnames(expression_subset)[-c(1)], ]
    sample_info_subset[["column_name"]] = NULL
    rownames(sample_info_subset) = NULL
    colnames(expression_subset)[-c(1)] = sample_info_subset$neat_sample_id

    print(head(expression_subset))
    print(head(sample_info_subset))

    return(list(sample_info = sample_info_subset,
                expression = expression_subset))
}

build_probe_to_gene <- function() {
    affy = read.csv("data/datasets/HG-U133_Plus_2.na36.annot.csv", skip=25)
    probe_to_gene = do.call(rbind, apply(affy , MARGIN=1, FUN=create_gene_rows))
    # Replace --- with NA
    probe_to_gene = mutate(probe_to_gene, gene = na_if(gene, "---"))
    # Remove any rows with NA
    probe_to_gene = probe_to_gene[!(rowSums(is.na(probe_to_gene))), ]

    return (probe_to_gene)
}

sample_info = build_sample_info()

chaussabelA = build_sample_expression_subset(sample_info, "data/datasets/chaussabel/E-GEOD-11907-processed-data-1673830054.txt")
chaussabelB = build_sample_expression_subset(sample_info, "data/datasets/chaussabel/E-GEOD-11907-processed-data-1673830055.txt")

print("Labels should match up:")
print(all(chaussabelA$sample_info$neat_sample_id ==
          colnames(chaussabelA$expression)[-c(1)]))
print(all(length(chaussabelA$sample_info$neat_sample_id) ==
          length(colnames(chaussabelA$expression)[-c(1)])))
print(all(chaussabelB$sample_info$neat_sample_id ==
          colnames(chaussabelB$expression)[-c(1)]))
print(all(length(chaussabelB$sample_info$neat_sample_id) ==
          length(colnames(chaussabelB$expression)[-c(1)])))

print("neat_sample_id should have unique values")
print(length(unique(chaussabelA$sample_info$neat_sample_id)) == length(chaussabelA$sample_info$neat_sample_id))
print(length(unique(chaussabelB$sample_info$neat_sample_id)) == length(chaussabelB$sample_info$neat_sample_id))

write.table(chaussabelA$sample_info, "data/datasets/chaussabelA/sample_info.tsv", sep='\t', row.names=FALSE, col.names=TRUE)
write.table(chaussabelA$expression, "data/datasets/chaussabelA/expression.tsv", sep='\t', row.names=FALSE, col.names=TRUE)
write.table(chaussabelB$sample_info, "data/datasets/chaussabelB/sample_info.tsv", sep='\t', row.names=FALSE, col.names=TRUE)
write.table(chaussabelB$expression, "data/datasets/chaussabelB/expression.tsv", sep='\t', row.names=FALSE, col.names=TRUE)

probe_to_gene = build_probe_to_gene()
print(head(probe_to_gene))
write.table(probe_to_gene, "data/datasets/chaussabelA/probe_to_gene.tsv", sep='\t', row.names=FALSE, col.names=TRUE)
write.table(probe_to_gene, "data/datasets/chaussabelB/probe_to_gene.tsv", sep='\t', row.names=FALSE, col.names=TRUE)
