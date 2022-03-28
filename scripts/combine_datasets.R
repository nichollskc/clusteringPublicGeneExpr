library(ggplot2)
library(stringr)

print(snakemake@input)
logfc_summaries = lapply(snakemake@input, function(x) read.csv(x, sep='\t', row.names=1))

csv_info = str_match(snakemake@input, 'data/datasets/(\\w+)/.*-(\\w+)\\.')
datasets = unique(csv_info[, 2])
signatures = unique(csv_info[, 3])
print(datasets)
print(signatures)

logfc_summaries_renamed = lapply(signatures, function(sig_name) {
    lapply(logfc_summaries[csv_info[, 3] == sig_name],
           function(df) {
                colnames(df) = paste0(colnames(df), "_", sig_name)
                df
           })
})
logfc_summaries_all_datasets = lapply(logfc_summaries_renamed,
                                      function(x) do.call(rbind, x))
combined = do.call(cbind, logfc_summaries_all_datasets)
print(head(combined))

write.csv(combined, "combined.csv")

library(ggplot2)
library(stringr)
combined$disease = sapply(rownames(combined), function(x) str_split(string=x, pattern="_")[[1]][1])
combined$dataset = sapply(rownames(combined), function(x) str_split(string=x, pattern="_")[[1]][2])
pca = prcomp(combined[combined$dataset != "SMITH", seq(2, 12, 2)])

with_pca_no_smith = cbind(combined[combined$dataset != "SMITH", ], pca$x)
ggplot(with_pca, aes(x=PC1, y=PC2, colour=dataset)) + geom_point(alpha=0.5)
ggplot(with_pca_no_smith, aes(x=PC1, y=PC2, colour=dataset)) + geom_point(alpha=0.5)
ggplot(with_pca_no_smith, aes(x=PC1, y=PC2, colour=disease)) + geom_point(alpha=0.5)

pca = prcomp(combined[combined$dataset %in% c("CHA", "COULSON"), seq(1, 12, 2)])
with_pca_2 = cbind(combined[combined$dataset %in% c("CHA", "COULSON"), ], pca$x)
g = ggplot(with_pca_2, aes(x=PC1, y=PC2, colour=disease)) + geom_point(alpha=0.5, aes(shape=dataset, size=total_se)) + theme(legend.position="bottom")
ggsave("plots/pca_2.png", g, width=10, height=10)

pca = prcomp(combined[combined$dataset %in% c("CHA", "CHB", "COULSON"), seq(1, 12, 2)])
with_pca_3 = cbind(combined[combined$dataset %in% c("CHA", "CHB", "COULSON"), ], pca$x)
g = ggplot(with_pca_3, aes(x=PC1, y=PC2, colour=disease)) + geom_point(alpha=0.5, aes(shape=dataset, size=total_se)) + theme(legend.position="bottom")
ggsave("plots/pca_3.png", g, width=10, height=10)

pca = prcomp(combined[, seq(1, 12, 2)])
with_pca = cbind(combined, pca$x)
g = ggplot(with_pca, aes(x=PC1, y=PC2, colour=disease)) + geom_point(alpha=0.5, aes(shape=dataset, size=total_se)) + theme(legend.position="bottom")
ggsave("plots/pca_4.png", g, width=10, height=10)
