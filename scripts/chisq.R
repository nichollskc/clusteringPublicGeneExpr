library(pheatmap)
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(tibble)
library(R.cache)

generate_balanced_colours <- function(obj) {
  paletteLength = 100
  customColours <- colorRampPalette(c("firebrick3", "white", "navy"))(paletteLength)
  customBreaks <- c(seq(min(obj), 0, length.out=ceiling(paletteLength/2) + 1),
                    seq(max(obj)/paletteLength, max(obj), length.out=floor(paletteLength/2)))
  return(list("colours"=customColours, "breaks"=customBreaks))
}

raw_fisher_pval <- function(tab, B) {
    fisher = fisher.test(tab, simulate.p.value=TRUE, B=B)
    return (fisher)
}

calc_fisher_pval <- addMemoization(raw_fisher_pval)

get_chisq_pvalue <- function(df, var1, var2) {
    tab = table(df[, c(var1, var2)])
    props = tab / rowSums(tab)

    print(tab)
    if (length(unique(df[[var1]])) < 2) {
        pval = 1
        print("Can't apply chisq.test since only one cluster")
    } else {
        chi_res = chisq.test(tab, simulate.p.value=TRUE, B=1000000)
        print(chi_res)
        pval = chi_res$p.value
    }

    return (list(signature=var1,
                 pval=pval,
                 freq_counts=tab,
                 proportions=props))
}

get_fisher_pvalue <- function(df, var1, var2) {
    tab = table(df[, c(var1, var2)])
    props = prop.table(tab, 2)
    print("Proportions")
    print(props)

    print(tab)
    if (length(unique(df[[var1]])) < 2) {
        pval = 1
        print("Can't apply fisher test since only one cluster")
    } else {
        fisher = calc_fisher_pval(tab, 1e8)
        print(fisher)
        pval = fisher$p.value
    }

    return (list(signature=var1,
                 pval=pval,
                 freq_counts=tab,
                 proportions=props))
}

heatmap_proportions <- function(enrichment_res) {
    pheatmap(enrichment_res$proportions,
             display_numbers=enrichment_res$freq_counts,
             color=colorRampPalette((RColorBrewer::brewer.pal(n = 7,
                                                              name = "Blues")))(100),
             breaks=seq(0, 1, length.out=101),
             number_color = "#888888",
             cluster_rows=FALSE,
             cluster_cols=FALSE,
             legend=FALSE,
             main=paste0(enrichment_res$signature, " (p:", signif(enrichment_res$pval, digits=2), ")"))[[4]]
}

heatmap_mean_by_cluster <- function(variable_name, combined_df, col_start="mean_logfc_") {
  clusterMeans = combined_df %>%
    group_by(get(variable_name)) %>%
    select(starts_with(col_start)) %>%
    summarise(across(everything(), mean)) %>%
    select(-c(1)) %>%
    rename_with(~ gsub(col_start, "", .x))
  print(col_start)
  print(clusterMeans)
  customColours = generate_balanced_colours(clusterMeans)
  pheatmap(clusterMeans,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = customColours$colours,
           breaks = customColours$breaks)[[4]]
}

plot_props_obs <- function(signature, combined_df, disease_var="disease") {
#  chisq_res = get_chisq_pvalue(combined_df, signature, disease_var)
  fisher_res = get_fisher_pvalue(combined_df, signature, disease_var)
  props_heatmap = heatmap_proportions(fisher_res)
  clustermeans_heatmap = heatmap_mean_by_cluster(signature, combined_df)

  all_grobs = c(props_heatmap$grobs, list(textGrob("Counts/proportions by cluster")),
                clustermeans_heatmap$grobs, list(textGrob("mean logfc by cluster")))

  g = grid.arrange(grobs=all_grobs, widths=c(10,1,10,1,1.5), heights=c(1,1,10,6),
               layout_matrix=rbind(c(1,1,1,1,1),
                                   c(5,NA,11,NA,NA),
                                   c(2,4,6,8,9),
                                   c(3,NA,7,NA,NA)))
  ggsave(paste0(plots_dir, "/fisher_", disease_var, "_", signature, ".png"), g)

  return(fisher_res$pval)
}

plot_props_obs_latents <- function(signature, combined_df, disease_var="disease") {
  chisq_res = get_chisq_pvalue(combined_df, signature, disease_var)
  props_heatmap = heatmap_proportions(chisq_res)
  clustermeans_heatmap = heatmap_mean_by_cluster(signature, combined_df)

  latents_file = paste0("results/",
                        gsub("\\.", "/", str_split(signatures[[1]], "__")[[1]][2]),
                        "/meanLatents.csv")
  latentMeans = read.csv(latents_file, row.names=1)
  latentMeans = latentMeans %>%
    rownames_to_column("neat_sample_id") %>%
    rename_with(~ gsub("mean", "latent_mean", .x))

  df_with_latents = combined_df %>%
    rownames_to_column("neat_sample_id") %>%
    merge(latentMeans, by="neat_sample_id") %>%
    column_to_rownames("neat_sample_id")
  print(head(df_with_latents))

  latentclustermeans_heatmap = heatmap_mean_by_cluster(signature, df_with_latents, col_start="latent_mean_logfc_")

  with_latents = c(props_heatmap$grobs, list(textGrob("Counts/proportions by cluster")),
                   clustermeans_heatmap$grobs, list(textGrob("mean logfc by cluster")),
                   latentclustermeans_heatmap$grobs, list(textGrob("mean latent logfc by cluster")))

  g = grid.arrange(grobs=with_latents, widths=c(10,1,10,1,1.5,(5/3) * (ncol(latentMeans) - 1),1,1.5),
                   heights=c(1,1,10,6),
               layout_matrix=rbind(c(1,1,1,1,1,1,1,1),
                                   c(5,NA,10,NA,NA,15,NA,NA),
                                   c(2,4,6,8,9,11,13,14),
                                   c(3,NA,7,NA,NA,12,NA,NA)))
  ggsave(paste0(plots_dir, "/chisq_", disease_var, "_", signature, ".png"), g)

  return(chisq_res$pval)
}

# Create HC variable
sample_info = read.csv(snakemake@input[["sample"]], sep='\t') %>%
  mutate(HC = case_when(disease == "HC" ~ "HC",
                        TRUE ~ "nonHC"))
if (! "subtype" %in% colnames(sample_info)) {
    sample_info["subtype"] = sample_info["disease"]
}
all_calls = readRDS(snakemake@input[["calls"]])

obsData = read.csv(snakemake@input[["obs"]], sep='\t', row.names=1) %>%
  rownames_to_column("neat_sample_id")
obsVars = read.csv(snakemake@input[["var"]], sep='\t', row.names=1)

df = data.frame(name=colnames(all_calls),
                signature = sapply(colnames(all_calls), function(x) strsplit(x, '\\.')[[1]][2]),
                novar=ifelse(grepl("novar", colnames(all_calls)), "novar", ""),
                method=sapply(colnames(all_calls), function(x) strsplit(x, '__')[[1]][1]),
                K=apply(all_calls, 2, function(x) length(unique(x))))

combined = all_calls %>%
  select(rownames(df)) %>%
  rownames_to_column("neat_sample_id") %>%
  merge(sample_info[, c("neat_sample_id", "disease", "HC", "subtype")], by="neat_sample_id") %>%
  merge(obsData, by="neat_sample_id") %>%
  column_to_rownames("neat_sample_id")

write.csv(combined, snakemake@output[["combined"]])

signatures = rownames(df)
plots_dir = dirname(snakemake@output[["plot"]])

pvals = data.frame("HC" = rep(1, length(signatures)),
                   "disease" = rep(1, length(signatures)),
                   row.names=signatures)
if (! all(combined["subtype"] == combined["disease"])) {
    pvals["subtype"] = 1
}

print(pvals)

#print("Sig and trait")
#print(snakemake@wildcards[["signature"]])
#print(snakemake@wildcards[["trait"]])
#plot_props_obs(snakemake@wildcards[["signature"]],
#               combined,
#               snakemake@wildcards[["trait"]])
#
print(pvals)
for (trait in colnames(pvals)) {
  for (signature in rownames(pvals)) {
    print(paste(trait, signature))
    pvals[signature, trait] = plot_props_obs(signature, combined, trait)
  }
  print(pvals)
}

write.csv(pvals, snakemake@output[["pvals"]])
