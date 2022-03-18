library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)
library(tidyr)

library(vsn)

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
    sample_info = NULL

    if (is.null(sample_info)) {
        sample_info = try_extract_sample_info(sample_id, "^SLE(?<pat_num>[0-9]+)_?(?<extra>.*?)_?(?<AB>A|B|)$", "SLE")
        if (! is.null(sample_info)) {
            sample_info[4] = gsub("/", ".", sample_info[4])
        }
    }
    if (is.null(sample_info)) {
        sample_info = try_extract_sample_info(sample_id, "^Lupus_(?<pat_num>[0-9]+|pre)(?<extra>)(?<AB>A|B)", "SLE")
    }
    if (is.null(sample_info)) {
        sample_info = try_extract_sample_info(sample_id, "^sys(?<pat_num>[0-9]+)(?<AB>[AB])[ \\.]?(?<extra>[0-9]{6}|)", "JIA")
    }
    if (is.null(sample_info)) {
        sample_info = try_extract_sample_info(sample_id, "^RET[ \\.]?(?<pat_num>[0-9]+)(?<extra>Froz|)(?<AB>A|B|)", "RET")
    }
    if (is.null(sample_info)) {
        sample_info = try_extract_sample_info(sample_id, "^(?i)Mel[ \\.](?<pat_num>[0-9]+)[ \\.](?<extra>pre).*U133[ \\.]?(?<AB>A|B)", "MEL")
    }
    if (is.null(sample_info)) {
        sample_info = try_extract_sample_info(sample_id, "^(?i)INF(?<pat_num>[0-9]+)(?<AB>A|B|).*Ecoli(?<extra>)", "InfECOLI")
    }
    if (is.null(sample_info)) {
        sample_info = try_extract_sample_info(sample_id, "^(?i)INF(?<pat_num>[0-9]+)(?<AB>A|B|).*STAPH(?<extra>)", "InfSTAPH")
    }
    if (is.null(sample_info)) {
        sample_info = try_extract_sample_info(sample_id, "^(?i)IDDM[\\.]?(?<pat_num>[0-9]*)(?<extra>)[ \\.]?(?<AB>A|B)", "T1D")
    }
    if (is.null(sample_info)) {
        sample_info = try_extract_sample_info(sample_id, "^X?003.*ACC[\\.-](?<pat_num>[0-9]+)(?<extra>F|)[\\._](?:U133)?(?<AB>B|)", "MEL")
    }
    if (is.null(sample_info)) {
        sample_info = try_extract_sample_info(sample_id, "^(?i)(?:PBMCS__|)H[ \\.-]?(?<pat_num>[0-9]+|[A-Z]{2})_?(?<extra>[0-9]{6}|)(?:Healthy_U133)?(?<AB>B|)", "HC")
    }

    if (is.null(sample_info)) {
        sample_info = c(sample_id, "", "", "", "")
    }
    return(sample_info)
}

add_detailed_sample_info <- function(sample_df, full_ids) {
    res = lapply(full_ids, extract_sample_id)
    df = data.frame(do.call(rbind, res))
    colnames(df) = c("old_sample_id", "disease", "patient_num", "extra", "AorB")
    df$dataset = "CHA"
    df$patient_id = paste0(df$disease, "_", df$dataset, "_", df$patient_num)
    df$neat_sample_id = paste0(df$patient_id, "_", df$extra, "_", df$AorB)

    #n_occur <- data.frame(table(df$neat_sample_id))
    #df[df$neat_sample_id %in% n_occur[n_occur$Freq > 1, "Var1"], ]

    return(cbind(sample_df, df))
}

fix_chaussabel_columns <- function(processed_data) {
    processed_data = processed_data[, c(FALSE, TRUE)]
    colnames(processed_data) = gsub("\\.1", "", colnames(processed_data))
    processed_data = processed_data[-c(1), ]
    processed_data = processed_data %>% mutate(across(everything(), as.numeric))
    return (processed_data)
}

read_gene_list <- function(gene_list_file) {
    genes_table = read.table(gene_list_file, header=1, sep=',')
    return (genes_table$SYMBOL)
}

std <- function(x) sd(x)/sqrt(length(x))

calculate_average_logfc <- function(gene_means, name) {
    hc_medians = filter(gene_means, disease == "HC") %>% summarise(across(starts_with("GENE_"), median))

    logfc = data.frame(mapply('-', gene_means %>% select(starts_with("GENE_")), hc_medians))
    logfc = cbind(logfc, select(gene_means, !starts_with("GENE_")))
    colnames(logfc) = colnames(gene_means)
    rownames(logfc) = rownames(gene_means)

    to_plot = pivot_longer(logfc, cols=starts_with("GENE_"), names_to="GENE", names_prefix="GENE_")
    g = ggplot(to_plot, aes(x=GENE, y=value, color=disease))  + geom_point() + facet_wrap("~disease")
    ggsave(paste0("plots/chaussabel_", name, ".png"), g)

    average_logfc_individuals = to_plot %>%
        select(neat_sample_id, GENE, value) %>%
        group_by(neat_sample_id) %>%
        summarise(mean_logfc = mean(value),
                  se_logfc = std(value),
                  .groups="keep") %>%
        arrange(neat_sample_id)

    return(average_logfc_individuals)
}

apply_vsn_transformation <- function(processed_data, name) {
    g = meanSdPlot(as.matrix(processed_data %>% mutate(across(everything(), function(x) log2(x))) ))$gg
    ggsave(paste0("plots/chaussabel_mean_variance_ranks_log_", name, ".png"), g)
    g = meanSdPlot(as.matrix(processed_data %>% mutate(across(everything(), function(x) log2(x))) ), ranks=FALSE)$gg
    ggsave(paste0("plots/chaussabel_mean_variance_log_", name, ".png"), g)

    after_vsn = vsn2(as.matrix(processed_data))

    g = meanSdPlot(after_vsn)$gg
    ggsave(paste0("plots/chaussabel_mean_variance_ranks_vsn_", name, ".png"), g)
    g = meanSdPlot(after_vsn, ranks=FALSE)$gg
    ggsave(paste0("plots/chaussabel_mean_variance_vsn_", name, ".png"), g)

    df = data.frame(as.matrix(after_vsn))
    colnames(df) = colnames(processed_data)
    rownames(df) = rownames(processed_data)
    return (df)
}

calculate_ranked_data_chaussabel <- function(genelist, processed_data, name) {
    dir_genelists = "~/rds/rds-cew54-basis/People/KATH/publicGeneExpr/"
    genelist_file = paste0(dir_genelists, "data/pathways/processed/", genelist, ".csv")
    gene_list <- read_gene_list(genelist_file)

    processed_data$Probe.Set.ID = rownames(processed_data)

    affy = read.csv("data/datasets/HG-U133_Plus_2.na36.annot.csv", skip=25)
    print(paste0("Total gene list length: ", length(gene_list)))
    probes_in_list <- affy[affy$Gene.Symbol %in% gene_list, c("Gene.Symbol", "Probe.Set.ID")]
    probe_values_in_list = merge(processed_data, probes_in_list, by="Probe.Set.ID")
    print(paste0("Genes also found in dataset: ", length(unique(probe_values_in_list$Gene.Symbol))))

    gene_means = probe_values_in_list %>%
        select(!starts_with("Probe.Set.ID")) %>%
        mutate(across(!starts_with("Gene.Symbol"), as.numeric)) %>%
        mutate(across(!starts_with("Gene.Symbol"), log2)) %>%
        group_by(Gene.Symbol) %>%
        summarise(across(!starts_with("Gene.Symbol"), mean)) %>%
        column_to_rownames(var = "Gene.Symbol") %>%
        t()
    gene_names = paste0("GENE_", colnames(gene_means))
    colnames(gene_means) = gene_names
    gene_means = add_detailed_sample_info(gene_means, rownames(gene_means))
    to_plot = pivot_longer(gene_means, cols=starts_with("GENE_"), names_to="GENE", names_prefix="GENE_")
    g = ggplot(to_plot, aes(x=GENE, y=value, color=disease))  + geom_point() + facet_wrap("~disease")
    ggsave(paste0("plots/chaussabel_logtransformed_", name, ".png"), g)

    average_logfc_individuals = calculate_average_logfc(gene_means, name)
    write.csv(average_logfc_individuals, paste0("data/datasets/chaussabel/", name, ".csv"))

    return(average_logfc_individuals)
}

compare_datasets <- function(suffix="") {
    ifnA = read.csv(paste0("data/datasets/chaussabel/ifn_A", suffix, ".csv"))
    ifnB = read.csv(paste0("data/datasets/chaussabel/ifn_B", suffix, ".csv"))
    ifnSMITH = read.csv(paste0("data/datasets/smith/ifn", suffix, ".csv"))
    ifnA$dataset = "CHA_A"
    ifnB$dataset = "CHA_B"
    ifnSMITH$dataset = "SMITH"

    both = rbind(ifnA, ifnB, ifnSMITH)
    both$disease = sapply(both$neat_sample_id, function(x) str_split(string=x, pattern="_")[[1]][1])
    g = ggplot(both, aes(x=dataset, y=mean_logfc, colour=disease)) + geom_point() + facet_wrap("~ disease") + geom_boxplot()
    ggsave(paste0("plots/all_compared", suffix, ".png"), g, width = 7, height = 8, dpi = 300, units = "in")
    g = ggplot(both, aes(x=dataset, y=se_logfc, colour=disease)) + geom_point() + facet_wrap("~ disease") + geom_boxplot()
    ggsave(paste0("plots/all_compared_se", suffix, ".png"), g, width = 7, height = 8, dpi = 300, units = "in")
}

full_sample_info = read.csv("data/datasets/chaussabel/sample_info.txt", sep='\t')
full_sample_info = add_detailed_sample_info(full_sample_info, full_sample_info[, "Scan.Name"])
full_sample_info$dataset = ifelse(full_sample_info$Derived.Array.Data.Matrix.File == "E-GEOD-11907-processed-data-1673830054.txt", "CHA", "CHB")
full_sample_info$neat_sample_id = paste0(full_sample_info$disease, "_", full_sample_info$dataset, "_", full_sample_info$patient_num, "_", full_sample_info$extra)

dfA = read.csv("data/datasets/chaussabel/E-GEOD-11907-processed-data-1673830054.txt",
               sep='\t', row.names=1)
dfA = fix_chaussabel_columns(dfA)
dfB = read.csv("data/datasets/chaussabel/E-GEOD-11907-processed-data-1673830055.txt",
               sep='\t', row.names=1)
dfB = fix_chaussabel_columns(dfB)

processed_samples = rbind(data.frame(column_name = colnames(dfA),
                                     file = "E-GEOD-11907-processed-data-1673830054.txt"),
                          data.frame(column_name = colnames(dfB),
                                     file = "E-GEOD-11907-processed-data-1673830055.txt")) %>%
    filter(!grepl("[.]1$", column_name))
processed_samples = add_detailed_sample_info(processed_samples, processed_samples[, "column_name"])
processed_samples$dataset = ifelse(processed_samples$file == "E-GEOD-11907-processed-data-1673830054.txt", "CHA", "CHB")
processed_samples$neat_sample_id = paste0(processed_samples$disease, "_", processed_samples$dataset, "_", processed_samples$patient_num, "_", processed_samples$extra)

# There are just 6 samples not found in one of the processed files
combined = merge(processed_samples, full_sample_info, by="neat_sample_id")
write.table(combined, "data/datasets/chaussabel/neat_sample_info.txt")

dfA_vsn = apply_vsn_transformation(dfA, "A")
dfB_vsn = apply_vsn_transformation(dfB, "B")
write.csv(dfA_vsn, "data/datasets/chaussabel/after_vsnA.csv")
write.csv(dfB_vsn, "data/datasets/chaussabel/after_vsnB.csv")
#dfA_vsn = read.csv("data/datasets/chaussabel/after_vsnA.csv", row.names=1)
#dfB_vsn = read.csv("data/datasets/chaussabel/after_vsnB.csv", row.names=1)
calculate_ranked_data_chaussabel("ifn", dfA_vsn, "ifn_A_vsn")
calculate_ranked_data_chaussabel("ifn", dfB_vsn, "ifn_B_vsn")

calculate_ranked_data_chaussabel("ifn", dfA, "ifn_A")
calculate_ranked_data_chaussabel("ifn", dfB, "ifn_B")

#compare_datasets("")
#compare_datasets("_vsn")
