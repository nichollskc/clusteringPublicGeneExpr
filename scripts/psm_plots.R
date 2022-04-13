source("scripts/utils.R")
library(cluster)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(mcclust)
library(mclust)
library(pheatmap)

psm_plots <- function(dataset, datasets, name, focus_dataset=NULL) {
    obsData = read.table(paste0(dataset, "/obsData.tsv"),
                         header=1, row.names=1)
    obsVars = read.table(paste0(dataset, "/obsVars.tsv"),
                         header=1, row.names=1)

    result = calc_psms(datasets)
    bigpsm = result$bigpsm
    psms = result$psms

    print(isSymmetric(bigpsm))
    print(max(bigpsm))
    print(min(bigpsm))
    print(diag(bigpsm))
    rownames(bigpsm) = rownames(obsData)
    colnames(bigpsm) = rownames(obsData)

    calls=maxpear(bigpsm, method="comp") ## calls

    # MClust solution
    print("Calculating mclust solution")
    mclust_solution <- calc_mclust(obsData)
    print(mclust_solution)

    scaled_obsData = DPMUnc::scale_data(obsData, obsVars)$data
    print("Calculating mclust solution")
    mclust_solution_scaled <- calc_mclust(scaled_obsData)
    print(mclust_solution_scaled)

    gap_stat <- clusGap(obsData,
                        FUN = kmeans,
                        nstart = 25,
                        K.max = 20,
                        B = 50)
    print(gap_stat)
    best_K = which.max(gap_stat$Tab[, 3])
    kmeans_solution = kmeans(obsData, centers=best_K, nstart=25)

    gap_stat <- clusGap(scaled_obsData,
                        FUN = kmeans,
                        nstart = 25,
                        K.max = 20,
                        B = 50)
    print(gap_stat)
    best_K = which.max(gap_stat$Tab[, 3])
    kmeans_solution_scaled = kmeans(scaled_obsData, centers=best_K, nstart=25)

    annotations = get_ann_colors(calls$cl, mclust_solution$classification, obsData)

    # Same hclust calculation that maxpear does prior to choosing optimal cut point
    hclust.comp <- hclust(as.dist(1 - bigpsm), method = "complete")
    psm_heatmap = pheatmap(bigpsm,
                           show_rownames = TRUE,
                           show_colnames = FALSE,
                           cluster_rows = hclust.comp,
                           cluster_cols = hclust.comp,
                           annotation_names_row = FALSE,
                           treeheight_col=0,
                           fontsize_row=6,
                           annotation_row = annotations$ann,
                           annotation_col = annotations$ann,
                           color=colorRampPalette((RColorBrewer::brewer.pal(n = 7,
                                                                            name = "Blues")))(100),
                           annotation_colors = annotations$colors,
                           width=20,
                           height=14,
                           filename=paste0("plots/", name, "/psm_heatmap.png"))

    generate_balanced_colours <- function(obj) {
        paletteLength = 100
        customColours <- colorRampPalette(c("firebrick3", "white", "navy"))(paletteLength)
        customBreaks <- c(seq(min(obj), 0, length.out=ceiling(paletteLength/2) + 1),
                          seq(max(obj)/paletteLength, max(obj), length.out=floor(paletteLength/2)))
        return(list("colours"=customColours, "breaks"=customBreaks))
    }

    customColours = generate_balanced_colours(obsData)
    obs_heatmap = pheatmap(obsData,
                           clustering_method="complete",
                           cluster_rows = hclust.comp,
                           annotation_colors = annotations$colors,
                           annotation_row = annotations$ann,
                           cluster_col = FALSE,
                           color = customColours$colours,
                           fontsize_col = 8,
                           fontsize_row = 6,
                           width=9,
                           height = 14,
                           breaks = customColours$breaks,
                           filename=paste0("plots/", name, "/obs_heatmap.png"))

    if(grepl("novar", name)) {
        min_val = min(obsVars)
        max_val = max(obsVars)
        customColours = generate_balanced_colours(c(min_val / 10, max_val * 10))
    } else {
        customColours = generate_balanced_colours(obsVars)
    }
    obs_vars_heatmap = pheatmap(obsVars,
                                clustering_method="complete",
                                cluster_rows = hclust.comp,
                                cluster_cols = FALSE,
                                annotation_colors = annotations$colors,
                                annotation_row = annotations$ann,
                                color = customColours$colours,
                                fontsize_col = 8,
                                fontsize_row = 6,
                                width=9,
                                height = 14,
                                breaks = customColours$breaks,
                                filename=paste0("plots/", name, "/obs_vars_heatmap.png"))

    print("Entries from different PSMs")
    print(sapply(psms, function(x) x[4, 5]))
    heatmaps = lapply(psms, function(x) pheatmap(x,
                                                 legend=FALSE,
                                                 color=colorRampPalette((RColorBrewer::brewer.pal(n = 7,
                                                                                                  name = "Blues")))(100),
                                                 cluster_rows = hclust.comp,
                                                 cluster_cols = hclust.comp,
                                                 treeheight_col=0,
                                                 border_color = NA,
                                                 treeheight_row=0))
    for (i in 1:10) {
        ggsave(paste0("plots/", name, "/psm_heatmap_seed_", i, ".png"), heatmaps[[i]])
    }
    heatmap_grid = grid.arrange(grobs=lapply(heatmaps, function(x) x[[4]]))
    ggsave(paste0("plots/", name, "/psm_heatmap_grid.png"), heatmap_grid)

    customColours = generate_balanced_colours(bigpsm - psms[[1]])
    diff_heatmaps = lapply(psms, function(x) pheatmap(bigpsm - x,
                                                 legend=FALSE,
                                                 show_colnames = FALSE,
                                                 show_rownames = FALSE,
                                                 color=customColours$colours,
                                                 breaks=customColours$breaks,
                                                 cluster_rows = hclust.comp,
                                                 cluster_cols = hclust.comp,
                                                 treeheight_col=0,
                                                 border_color = NA,
                                                 treeheight_row=0)[[4]])
    colours = customColours$colours[c(100, 76, 51, 26, 1)]
    breaks = sapply(customColours$breaks[c(101, 76, 51, 26, 1)], function(x) signif(x, 3))
    dummy_plot = ggplot(data.frame(x=as.factor(breaks))) +
        scale_color_manual(values=colours, breaks=breaks) +
        labs(color="Overall PSM - Seed PSM") +
        geom_point(aes(x, x, color=x))
    dummy_grobs <- ggplot_gtable(ggplot_build(dummy_plot))
    leg_index <- which(sapply(dummy_grobs$grobs, function(x) x$name) == "guide-box")
    legend = dummy_grobs$grobs[[leg_index]]

    print(class(legend))
    diff_heatmaps[[length(diff_heatmaps) + 1]] = legend
    heatmap_grid = grid.arrange(grobs=diff_heatmaps)
    ggsave(paste0("plots/", name, "/psm_heatmap_diff_grid.png"), heatmap_grid)

    individual_calls = do.call(cbind, lapply(psms, function(x) adjust_labels_B_to_match_A(calls$cl, maxpear(x, method="comp")$cl)))
    all_calls = data.frame(cbind(calls$cl, individual_calls))
    colnames(all_calls) = c("Overall", paste("Seed", 1:length(datasets)))

    print(all_calls)

    map_to_call_counts <- function(calls) {
        call_labels = paste0(LETTERS[1:max(calls)], " (", table(calls), ")")
        return(plyr::mapvalues(calls, from=1:max(calls), to=call_labels))
    }
    cell_labels = apply(all_calls, MARGIN=2, map_to_call_counts)

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
                             filename=paste0("plots/", name, "/calls_heatmap.png"))
    generate_balanced_colours <- function(obj) {
        paletteLength = 100
        customColours <- colorRampPalette(c("firebrick3", "white", "navy"))(paletteLength)
        customBreaks <- c(seq(min(obj), 0, length.out=ceiling(paletteLength/2) + 1),
                          seq(max(obj)/paletteLength, max(obj), length.out=floor(paletteLength/2)))
        return(list("colours"=customColours, "breaks"=customBreaks))
    }

    return(list(name=name, bigpsm=bigpsm, calls=calls, hclust.comp=hclust.comp, mclust_calls=mclust_solution$classification,
                mclust_scaled_calls=mclust_solution_scaled$classification, kmeans_calls=kmeans_solution, kmeans_scaled_calls=kmeans_solution_scaled))
}

dataset_name = snakemake@wildcards[["dataset"]]
dataset_folders = snakemake@params[["dataset_folders"]]
obsfile = snakemake@input[["obs"]]
resfile = snakemake@output[["result"]]

print(dataset_folders)
print(class(dataset_folders))
dataset_folders = as.character(dataset_folders)
print(class(dataset_folders))

result = psm_plots(dirname(obsfile), dataset_folders, name=dataset_name)
saveRDS(result, paste0("plots/", dataset_name, "/psm_data.rds"))
