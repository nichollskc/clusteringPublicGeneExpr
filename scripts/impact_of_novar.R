library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)

palette <- c("#77AADD", "#000000", "#9E0142", "#D53E4F", "#F46D43", "#FEE08B", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2", "#FDAE61", "#E6F598", "#771155", "#AA4488", "#CC99BB", "#114477", "#774411", "#EEEEEE", "#117777", "#117744", "#44AA77", "#88CCAA", "#777711", "#44AAAA", "#AAAA44", "#77CCCC", "#DDDD77", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

get_var_info <- function(dataset) {
  obsData = read.csv(paste0("data/processed/", dataset, "_vsn_mad_logfc/signatures_v2/obsData.tsv"),
                     sep='\t', row.names=1)
  obsVars = read.csv(paste0("data/processed/", dataset, "_vsn_mad_logfc/signatures_v2/obsVars.tsv"),
                     sep='\t', row.names=1)
  scaled = DPMUnc::scale_data(obsData, obsVars)

  all_calls = readRDS(paste0("plots/", dataset, "/all_calls.rds"))
  df = data.frame(name=colnames(all_calls),
                  signature = sapply(colnames(all_calls), function(x) strsplit(x, '\\.')[[1]][2]),
                  novar=ifelse(grepl("novar", colnames(all_calls)), "novar", "withvar"),
                  method=sapply(colnames(all_calls), function(x) strsplit(x, '__')[[1]][1]),
                  K=apply(all_calls, 2, function(x) length(unique(x))))
  reduced_df = df %>%
    filter(method == "DPMUnc" | novar == "") %>%
    select(!name)

  all_ari = read.csv(paste0("plots/", dataset, "/ari_all_calls.csv"), row.names=1)
  reduced_ari = all_ari[rownames(reduced_df), rownames(reduced_df)]

  return(list(dataset=dataset,
              obsData=obsData, obsVars=obsVars, scaled=scaled,
              runs_df=reduced_df, ari=reduced_ari))
}

datasets = lapply(c("smith", "chaussabelA", "chaussabelB", "coulson"),
                  get_var_info)

ari_vs_vars = lapply(datasets, function(x) {
  result = x$runs_df %>%
    filter(!grepl("signature", signature), method=="DPMUnc") %>%
    select(novar, signature, K) %>%
    pivot_wider(names_from=c(novar), values_from=c(K), names_prefix="K_")
  result$ari = sapply(result$signature, function(y) {
    runs = rownames(x$runs_df %>% filter(signature == y))
    return(x$ari[runs[[1]], runs[[2]]])
  })

  result$obsvars_median = apply(x$obsVars, 2, median)
  result$scaledvars_median = apply(x$scaled$vars, 2, median)
  result$dataset = x$dataset
  return(result)
})

ari_vs_var_combined = do.call(rbind, ari_vs_vars)

g = ggplot(ari_vs_var_combined, aes(x=scaledvars_median, y=ari)) +
  geom_point(aes(colour=signature, shape=dataset), size=4, alpha=0.7) +
  geom_smooth(method="lm") +
  scale_color_manual(values=palette[seq(1,length.out=6, by=3)]) +
  scale_shape_manual(values=15:18) +
  stat_cor(label.x=0.5,
           aes(label = paste(..r.label.., ..p.label.., sep=', ')),
           method = "pearson",
           output.type = "text") +
  labs(title = "Uncertainty vs impact of ignoring uncertainty",
       x = "Median uncertainty (after scaling)",
       y = "Similarity (ARI) between DPMUnc run with and without uncertainty")
ggsave("plots/uncertainty_impact.png", g)
