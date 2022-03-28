library(ggplot2)
library(stringr)

summary_name = snakemake@wildcards[["summary"]]

summaries = lapply(snakemake@input, function(x) read.csv(x, sep='\t', col.names=c("neat_sample_id", "mean_summary", "se_summary")))
combined = do.call(rbind, summaries)
combined$disease = sapply(combined$neat_sample_id, function(x) str_split(string=x, pattern="_")[[1]][1])
combined$dataset = sapply(combined$neat_sample_id, function(x) str_split(string=x, pattern="_")[[1]][2])

g = ggplot(combined, aes(x=dataset, y=mean_summary, colour=disease)) + geom_point() + facet_wrap("~ disease") + geom_boxplot() + labs(y=summary_name)
ggsave(snakemake@output[[1]], g, width = 7, height = 8, dpi = 300, units = "in")
g = ggplot(combined, aes(x=dataset, y=se_summary, colour=disease)) + geom_point() + facet_wrap("~ disease") + geom_boxplot() + labs(y=paste0("SE ", summary_name))
ggsave(snakemake@output[[2]], g, width = 7, height = 8, dpi = 300, units = "in")
