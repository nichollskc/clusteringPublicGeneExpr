library(ggplot2)
library(stringr)

logfc_summaries = lapply(snakemake@input, function(x) read.csv(x, sep='\t'))
combined = do.call(rbind, logfc_summaries)
combined$disease = sapply(combined$neat_sample_id, function(x) str_split(string=x, pattern="_")[[1]][1])
combined$dataset = sapply(combined$neat_sample_id, function(x) str_split(string=x, pattern="_")[[1]][2])

g = ggplot(combined, aes(x=dataset, y=mean_logfc, colour=disease)) + geom_point() + facet_wrap("~ disease") + geom_boxplot()
ggsave(snakemake@output[[1]], g, width = 7, height = 8, dpi = 300, units = "in")
g = ggplot(combined, aes(x=dataset, y=se_logfc, colour=disease)) + geom_point() + facet_wrap("~ disease") + geom_boxplot()
ggsave(snakemake@output[[2]], g, width = 7, height = 8, dpi = 300, units = "in")
