source("scripts/utils.R")
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(mcclust)
library(mclust)
library(pheatmap)
library(stringr)

compare_versions <- function(psm_results, base_run, plots_dir) {
    base_calls = psm_results[[base_run]]$calls$cl
    hclust.comp = psm_results[[base_run]]$hclust.comp

    print(lapply(psm_results, function(x) sort(unique(x$calls$cl))))
    all_calls = do.call(cbind, lapply(psm_results, function(x) adjust_labels_B_to_match_A(base_calls, x$calls$cl)))
    colnames(all_calls) = lapply(psm_results, function(x) str_replace(x$name, "/", "."))
    print("RDS files")
    print(colnames(all_calls))

    all_calls_mclust = do.call(cbind, lapply(psm_results, function(x) adjust_labels_B_to_match_A(base_calls, x$mclust_calls)))
    colnames(all_calls_mclust) = paste0("mclust__", colnames(all_calls))
    all_calls_mclust = all_calls_mclust %>% data.frame() %>% select(-contains("novar"))

    all_calls_kmeans = do.call(cbind, lapply(psm_results, function(x) adjust_labels_B_to_match_A(base_calls, x$kmeans_calls$cluster)))
    colnames(all_calls_kmeans) = paste0("kmeans__", colnames(all_calls))
    all_calls_kmeans = all_calls_kmeans %>% data.frame() %>% select(-contains("novar"))

    colnames(all_calls) = paste0("DPMUnc__", colnames(all_calls))
    all_calls = cbind(all_calls, all_calls_mclust, all_calls_kmeans)
    print(apply(all_calls, MARGIN=1, unique))

    print(all_calls)

    map_to_call_counts <- function(calls) {
        call_labels = paste0(LETTERS[1:max(calls)], " (", table(calls), ")")
        return(plyr::mapvalues(calls, from=1:max(calls), to=call_labels))
    }
    cell_labels = apply(all_calls, MARGIN=2, map_to_call_counts)

    obsData = read.table(snakemake@input[["obs"]],
                         header=1, row.names=1)
    obsVars = read.table(snakemake@input[["var"]],
                         header=1, row.names=1)

    calls_heatmap = pheatmap(all_calls,
                             display_numbers = cell_labels,
                             fontsize_number = 4,
                             number_color = "#CCCCCC",
                             show_rownames = TRUE,
                             show_colnames = TRUE,
                             color = palette[1:max(all_calls)],
                             breaks = seq(0.5, max(all_calls) + 0.5, by=1),
                             cluster_col = FALSE,
                             cluster_row = hclust.comp,
                             fontsize_row = 6,
                             width=8,
                             height=14,
                             filename=paste0(plots_dir, "/all_calls_heatmap.png"))
    print("Done calls heatmap")

    generate_balanced_colours <- function(obj) {
        paletteLength = 100
        customColours <- colorRampPalette(c("firebrick3", "white", "navy"))(paletteLength)
        customBreaks <- c(seq(min(obj), 0, length.out=ceiling(paletteLength/2) + 1),
                          seq(max(obj)/paletteLength, max(obj), length.out=floor(paletteLength/2)))
        return(list("colours"=customColours, "breaks"=customBreaks))
    }

    customColours = generate_balanced_colours(obsData)
    annotations = list(colors=list())
    annotations$ann = data.frame(all_calls)
    for (c in colnames(all_calls)) {
        annotations$colors[[c]] = palette[1:max(all_calls)]
    }
    print(annotations)
    print("About to do obs heatmap")
    obs_heatmap = pheatmap(obsData,
                           cluster_rows = hclust.comp,
                           cluster_col = FALSE,
                           color = customColours$colours,
                           annotation_colors = annotations$colors,
                           annotation_row = annotations$ann,
                           fontsize_col = 8,
                           fontsize_row = 6,
                           width=9,
                           height = 14,
                           breaks = customColours$breaks,
                           filename=paste0(plots_dir, "/all_calls_obs_heatmap.png"))

    customColours = generate_balanced_colours(obsVars)
    obs_heatmap = pheatmap(obsVars,
                           cluster_rows = hclust.comp,
                           cluster_col = FALSE,
                           color = customColours$colours,
                           annotation_colors = annotations$colors,
                           annotation_row = annotations$ann,
                           fontsize_col = 8,
                           fontsize_row = 6,
                           width=9,
                           height = 14,
                           breaks = customColours$breaks,
                           filename=paste0(plots_dir, "/all_calls_obs_vars_heatmap.png"))

    print("About to save rds")
    saveRDS(annotations$ann, file=paste0(plots_dir, "/all_calls.rds"))
    print("Saved RDS")

    plot_ari_for_datasets(all_calls, plots_dir, width=4000, height=4000)
}

rds_files = snakemake@input[grepl(".rds", snakemake@input)]
psm_results = list()
for (rds_file in rds_files) {
    print(rds_file)
    dataset = str_match(rds_file, "\\w*/([/\\w]*)/\\w*")[, 2]
    print(dataset)
    psm_results[[str_replace(dataset, "/", ".")]] = readRDS(rds_file)
}

base_run = 1

compare_versions(psm_results, base_run, dirname(snakemake@output[[1]]))
