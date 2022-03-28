library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(vsn)

plot_meansd <- function(data, name, extra_name) {
    g = meanSdPlot(as.matrix(data))$gg
    ggsave(paste0("plots/", name, "_mean_variance_ranks", extra_name, ".png"), g)
    g = meanSdPlot(as.matrix(data), ranks=FALSE)$gg
    ggsave(paste0("plots/", name, "_mean_variance", extra_name, ".png"), g)
}

apply_vsn_transformation <- function(expression, name) {
    plot_meansd(expression, name, "_log")

    unlogged_expression = expression %>% mutate(across(everything(), function(x) 2^x)) %>% as.matrix()
    plot_meansd(unlogged_expression, name, "_raw")

    vsn_expression = vsn2(as.matrix(unlogged_expression), subsample=min(nrow(unlogged_expression), 30000))
    plot_meansd(vsn_expression, name, "_vsn")

    df = data.frame(as.matrix(vsn_expression))
    colnames(df) = colnames(expression)
    rownames(df) = rownames(expression)
    return(df)
}

name = snakemake@wildcards[["dataset"]]
filename = snakemake@input[[1]]
output_file = snakemake@output[[1]]
expression = read.csv(filename, row.names=1, check.names=FALSE, sep='\t')

vsn_expression = apply_vsn_transformation(expression, name)
vsn_expression = cbind(probe_id=rownames(vsn_expression),
                       vsn_expression)
write.table(vsn_expression, output_file, row.names=FALSE, col.names=TRUE, sep='\t')
