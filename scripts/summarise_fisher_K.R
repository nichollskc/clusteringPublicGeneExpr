library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(tibble)

library(ComplexHeatmap)

construct_df_dataset <- function(dataset) {
  calls = read.table(paste0("plots/", dataset, "/calls_and_sample_info.csv"),
                     sep=',', header=TRUE)
  K = calls %>% select(contains("__")) %>%
    apply(2, function(x) length(unique(x))) %>%
    data.frame() %>%
    setNames("K") %>%
    rownames_to_column(var="method_sig")
  fisher = read.table(paste0("plots/", dataset, "/fisher_pvals.csv"),
                      sep=',', header=TRUE)

  df = merge(K, fisher, by.x = "method_sig", by.y="X")
  str_to_delete = paste0("__", dataset, "_vsn_mad_logfc")
  df$method = str_replace(str_replace(df$method_sig, str_to_delete, ""), "\\..*", "")
  df$signature = str_replace(df$method_sig, ".*\\.", "")

  trait_cols = intersect(c("HC", "disease", "subtype"), colnames(fisher))

  df = df %>%
    pivot_longer(cols = trait_cols, names_to="trait", values_to="fisher") %>%
    mutate(fisher_log = log10(fisher),
           dataset = dataset)
  return (df)
}

datasets = c("smith", "coulson1", "chaussabelA", "chaussabelB")
combined = do.call(rbind, lapply(datasets,
                                 construct_df_dataset)) %>%
  filter(trait == "disease") %>%
  mutate(fisher_q = p.adjust(fisher, method='BY'),
         fisher_log_q = log10(fisher_q),
         q_is_sig = fisher_q < 0.01)

disease_log_fishers = lapply(datasets,
       function(x) {
          combined %>%
            select(trait, method, signature, fisher_log_q, dataset) %>%
            pivot_wider(names_from = method, values_from = fisher_log_q) %>%
            filter(dataset == x) %>%
            select(unique(combined$method), signature) %>%
            column_to_rownames("signature")
       })

ht_list = NULL  ## Heatmap(...) + NULL gives you a HeatmapList object
for(i in 1:4) {
  ht_list = ht_list %v% Heatmap(disease_log_fishers[[i]],
                                col=c("blue", "white"),
                                cluster_rows = FALSE,
                                row_title = datasets[[i]],
                                cluster_columns = FALSE,
                                name = "log10 qval")
}

comparison = combined %>%
  select(trait, method, signature, fisher_log, dataset) %>%
  pivot_wider(names_from = method, values_from = fisher_log) %>%
  filter(trait == "disease") %>%
  mutate(DP_less_DPZ = DPMUnc < DPMUnc_novar,
         DP_minus_DPZ = DPMUnc - DPMUnc_novar)
table(comparison$DP_less_DPZ)
table(comparison$DP_less_DPZ, comparison$dataset)

K_counts = lapply(datasets,
                             function(x) {
                               combined %>%
                                 select(trait, method, signature, K, dataset) %>%
                                 pivot_wider(names_from = method, values_from = K) %>%
                                 filter(dataset == x) %>%
                                 select(unique(combined$method), signature) %>%
                                 column_to_rownames("signature")
                             })

ht_list2 = NULL  ## Heatmap(...) + NULL gives you a HeatmapList object
for(d in 1:4) {
  ht_list2 = ht_list2 %v% Heatmap(K_counts[[d]],
                                col=circlize::colorRamp2(c(1, 2, 3, 10, 40),
                                                         c("black", "grey", "white", "greenyellow", "green4")),
                                cell_fun = function(j, i, x, y, width, height, fill) {
                                    if (fill == "#000000FF") {
                                      textcol = "#FFFFFFFF"
                                    } else {
                                      textcol = "#000000FF"
                                    }
                                    grid.text(K_counts[[d]][i, j], x, y, gp = gpar(fontsize = 10,
                                                                                   col = textcol))
                                },
                                cluster_rows = FALSE,
                                row_title = datasets[[d]],
                                cluster_columns = FALSE,
                                name = "K")
}

ggplot(combined, aes(x=K)) + geom_histogram() + facet_grid( method ~ .)
ggsave("plots/K_histogram.png", width=4, height=6, units="in")

plots_dir = "plots/"

png(paste0(plots_dir, "/fisher_q_values.png"), width=5, height=8, units="in", res=1200)
draw(ht_list)
dev.off()

png(paste0(plots_dir, "/K_heatmap.png"), width=5, height=8, units="in", res=1200)
draw(ht_list2)
dev.off()
