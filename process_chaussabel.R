library(stringr)

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
        sample_info = try_extract_sample_info(sample_id, "^(?i)IDDM.(?<pat_num>[0-9]+)(?<extra>)[ \\.]?(?<AB>A|B)", "T1D")
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

add_detailed_sample_info <- function(sample_df, col_name) {
    res = lapply(sample_df[, col_name], extract_sample_id)
    df = data.frame(do.call(rbind, res))
    colnames(df) = c("old_sample_id", "disease", "patient_num", "extra", "AorB")
    df$neat_sample_id = paste0(df$disease, "_", df$patient_num, "_", df$extra, "_", df$AorB)

    #n_occur <- data.frame(table(df$neat_sample_id))
    #df[df$neat_sample_id %in% n_occur[n_occur$Freq > 1, "Var1"], ]

    return(cbind(sample_df, df))
}

full_sample_info = read.csv("data/datasets/chaussabel/sample_info.txt", sep='\t')
full_sample_info = add_detailed_sample_info(full_sample_info, "Scan.Name")

dfA = read.csv("data/datasets/chaussabel/E-GEOD-11907-processed-data-1673830054.txt",
               sep='\t', row.names=1)
dfA = apply(dfA, MARGIN=1:2, FUN=as.numeric)
dfB = read.csv("data/datasets/chaussabel/E-GEOD-11907-processed-data-1673830055.txt",
               sep='\t', row.names=1)
dfB = apply(dfB, MARGIN=1:2, FUN=as.numeric)

processed_samples = rbind(data.frame(column_name = colnames(dfA),
                                     file = "E-GEOD-11907-processed-data-1673830054.txt"),
                          data.frame(column_name = colnames(dfB),
                                     file = "E-GEOD-11907-processed-data-1673830055.txt")) %>%
    filter(!grepl("[.]1$", column_name))
processed_samples = add_detailed_sample_info(processed_samples, "column_name")

# There are just 6 samples not found in one of the processed files
combined = merge(processed_samples, full_sample_info, by="neat_sample_id")

values = df[rownames(df) %in% yme, sample_info$df_colname]
values = apply(values, MARGIN=1:2, FUN=as.numeric)

hc_values = df[rownames(df) %in% found$Probe.Set.ID, sample_info[sample_info$disease == "HC", ]$df_colname]
hc_values = apply(hc_values, MARGIN=1:2, FUN=as.numeric)
hc_median = apply(hc_values, MARGIN=1, FUN=median)
