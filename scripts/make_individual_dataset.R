library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)
library(tidyr)

filenames = snakemake@input
variables = str_match(filenames, "-(\\w*)\\.tsv")[, 2]

fetch_average_logfc <- function(filename) {
    df = read.csv(filename, sep='\t', row.names=1)
    label = str_match(filename, "-(\\w*)\\.tsv")[, 2]

    colnames(df) = paste0(colnames(df), "_", label)
    return(df)
}

logfc_data = lapply(filenames, fetch_average_logfc)
combined = do.call(cbind, logfc_data)

file_prefix = str_match(snakemake@output[[1]], "[\\w\\/]*_logfc")[, 1]
construct_filename <- function(signature_name, novar=FALSE, filetype="Data") {
    return(paste0(file_prefix,
                  ifelse(novar, "_novar", ""),
                  "/",
                  signature_name,
                  "/obs",
                  filetype,
                  ".tsv"))
}

is_avg_logfc_col = grepl("mean_logfc", colnames(combined))
obsData = combined[, is_avg_logfc_col]
write.table(obsData, construct_filename("signatures_v2"), sep='\t', row.names=TRUE, col.names=TRUE)
write.table(obsData, construct_filename("signatures_v2", novar=TRUE), sep='\t', row.names=TRUE, col.names=TRUE)

is_se_logfc_col = grepl("se_logfc", colnames(combined))
obsVars = combined[, is_se_logfc_col] ** 2
write.table(obsVars, construct_filename("signatures_v2", filetype="Vars"), sep='\t', row.names=TRUE, col.names=TRUE)
obsVars_novar = apply(obsVars, function(x) 1e-8, MARGIN=1:2)
write.table(obsVars_novar, construct_filename("signatures_v2", filetype="Vars", novar=TRUE), sep='\t', row.names=TRUE, col.names=TRUE)

write.table(obsData[, 2:sum(is_avg_logfc_col)], construct_filename("signatures_v3"), sep='\t', row.names=TRUE, col.names=TRUE)
write.table(obsVars[, 2:sum(is_avg_logfc_col)], construct_filename("signatures_v3", filetype="Vars"), sep='\t', row.names=TRUE, col.names=TRUE)

write.table(obsData[, 2:sum(is_avg_logfc_col)], construct_filename("signatures_v3", novar=TRUE), sep='\t', row.names=TRUE, col.names=TRUE)
write.table(obsVars_novar[, 2:sum(is_avg_logfc_col)], construct_filename("signatures_v3", filetype="Vars", novar=TRUE), sep='\t', row.names=TRUE, col.names=TRUE)

for (i in 1:length(variables)) {
    write.table(obsData[ , i, drop=FALSE], construct_filename(variables[[i]]), sep='\t', row.names=TRUE, col.names=TRUE)
    write.table(obsVars[ , i, drop=FALSE], construct_filename(variables[[i]], filetype="Vars"), sep='\t', row.names=TRUE, col.names=TRUE)

    write.table(obsData[ , i, drop=FALSE], construct_filename(variables[[i]], novar=TRUE), sep='\t', row.names=TRUE, col.names=TRUE)
    write.table(obsVars_novar[ , i, drop=FALSE], construct_filename(variables[[i]], filetype="Vars", novar=TRUE), sep='\t', row.names=TRUE, col.names=TRUE)
}
