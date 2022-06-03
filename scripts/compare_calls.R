source("scripts/utils.R")
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(mcclust)
library(mclust)
library(pheatmap)
library(stringr)

compare_versions <- function(psm_results, base_run, plots_dir, sample_info) {
    base_calls = psm_results[[base_run]]$calls$cl
    hclust.comp = psm_results[[base_run]]$hclust.comp

    #print(lapply(psm_results, function(x) sort(unique(x$calls$cl))))
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
    #print(all_calls_kmeans)

    colnames(all_calls) = paste0("DPMUnc__", colnames(all_calls))
    all_calls = cbind(all_calls, all_calls_mclust, all_calls_kmeans)
    #print(apply(all_calls, MARGIN=1, unique))

    colnames(all_calls) = str_replace(colnames(all_calls), "__.*_vsn_mad_logfc", "")
    colnames(all_calls) = str_replace(colnames(all_calls), "DPMUnc_novar", "DPMZeroUnc")

    #print(all_calls)

    map_to_call_counts <- function(calls) {
        print(table(factor(calls, levels=1:max(all_calls))))
        call_labels = paste0(LETTERS[1:max(all_calls)],
                             " (",
                             table(factor(calls, levels=1:max(all_calls))),
                             ")")
        return(plyr::mapvalues(calls, from=1:max(all_calls), to=call_labels))
    }
    cell_labels = apply(all_calls, MARGIN=2, map_to_call_counts)

    obsData = read.table(snakemake@input[["obs"]],
                         header=1, row.names=1)
    obsVars = read.table(snakemake@input[["var"]],
                         header=1, row.names=1)

    generate_balanced_colours <- function(obj) {
        paletteLength = 100
        customColours <- colorRampPalette(c("firebrick3", "white", "navy"))(paletteLength)
        customBreaks <- c(seq(min(obj), 0, length.out=ceiling(paletteLength/2) + 1),
                          seq(max(obj)/paletteLength, max(obj), length.out=floor(paletteLength/2)))
        return(list("colours"=customColours, "breaks"=customBreaks))
    }

    annotations = list(colors=list())
    annotations$ann = data.frame(all_calls)
    colnames(annotations$ann) = colnames(all_calls) %>% str_replace("__with_subtypes_noblood_noGA_001", "")
    for (c in colnames(annotations$ann)) {
        annotations$colors[[c]] = palette[1:max(all_calls)]
    }
    print(annotations)

    palette <- c("#77AADD", "#CCCCCC", "#D53E4F", "#FEE08B", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2", "#FDAE61", "#E6F598", "#771155", "#AA4488", "#CC99BB", "#114477", "#774411", "#EEEEEE", "#117777", "#117744", "#44AA77", "#88CCAA", "#777711", "#44AAAA", "#AAAA44", "#77CCCC", "#DDDD77", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    print("Reduced calls")
    print(all_calls %>% select(contains("signatures_v2")))
    ch_calls_df = all_calls %>% select(contains("signatures_v2"))
    colnames(ch_calls_df) = str_replace(colnames(ch_calls_df), ".signatures_v2", "")
    ch_calls = Heatmap(ch_calls_df,
                       show_heatmap_legend = FALSE,
                       name = "calls",
                       cluster_rows = hclust.comp,
                       show_row_names = FALSE,
                       column_names_gp = gpar(fontsize = 10),
                       cluster_columns = FALSE,
                       col = palette[1:max(all_calls)],
#                      cell_fun = function(j, i, x, y, width, height, fill) {
#                          grid.text(substr(cell_labels[i, j], 1, 2), x, y, gp = gpar(fontsize = 10))
#                      },
                       column_title = "Clusters",
                       width = unit(1.2, "in"),
                       height = unit(8, "in"))

    if ("subtype" %in% colnames(sample_info)) {
        sample_info["disease"] = sample_info["subtype"]
    }
    ch_info = Heatmap(sample_info[, c("disease"), drop=FALSE],
                      name = "disease",
                      cluster_rows = hclust.comp,
                      col = cbbPalette[1:length(unique(sample_info$disease))],
                      show_row_names = FALSE,
                      column_names_gp = gpar(fontsize = 8),
                      width = unit(0.2, "in"),
                      height = unit(8, "in"))

    ch_obs = Heatmap(obsData,
                     name = "components",
                     cluster_rows = hclust.comp,
                     show_row_names = FALSE,
                     column_names_gp = gpar(fontsize = 8),
                     cluster_columns = FALSE,
                     col = circlize::colorRamp2(c(-max(abs(obsData)), 0, max(abs(obsData))),
                                                c("firebrick3", "white", "navy")),
                     heatmap_legend_param = list(direction = "horizontal"),
                     column_title = "Basis components",
                     width = unit(1.5, "in"),
                     height = unit(8, "in"))

    ch_vars = Heatmap(obsVars,
                     name = "uncertainty",
                     cluster_rows = hclust.comp,
                     cluster_columns = FALSE,
                       show_row_names = FALSE,
                     column_names_gp = gpar(fontsize = 8),
                     row_names_gp = gpar(fontsize = 10),
                     col = circlize::colorRamp2(c(0, max(abs(obsVars))),
                                                c("white", "green4")),
                     heatmap_legend_param = list(direction = "horizontal"),
                     column_title = "Uncertainty",
                     width = unit(1.5, "in"),
                     height = unit(8, "in"))

    png(paste0(plots_dir, "/combined_heatmap.png"), width=7.5, height=12, units="in", res=1200)
    draw(ch_calls + ch_info + ch_obs + ch_vars, merge_legend=TRUE, heatmap_legend_side = "bottom")
    dev.off()

    ch_psm = Heatmap(psm_results[[base_run]]$bigpsm,
                     name = "PSM",
                     cluster_columns = hclust.comp,
                     cluster_rows = hclust.comp,
                     show_column_names = FALSE,
                     row_names_gp = gpar(fontsize = 10),
                     col = circlize::colorRamp2(c(0, 1),
                                                c("white", "blue4")),
                     width = unit(8, "in"),
                     height = unit(8, "in"))

    png(paste0(plots_dir, "/psm_with_calls.png"), width=14, height=10, units="in", res=1200)
    draw(ch_calls + ch_psm)
    dev.off()

    calls_heatmap = pheatmap(all_calls,
                             display_numbers = cell_labels,
                             fontsize_number = 4,
                             number_color = "#333333",
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

    customColours = generate_balanced_colours(obsData)
    print("About to do obs heatmap")
    obs_heatmap = pheatmap(obsData,
                           cluster_rows = hclust.comp,
                           cluster_col = FALSE,
                           color = customColours$colours,
                           breaks = customColours$breaks,
                           annotation_colors = annotations$colors,
                           annotation_row = annotations$ann,
                           fontsize_col = 8,
                           fontsize_row = 6,
                           width=9,
                           height = 14,
                           filename=paste0(plots_dir, "/all_calls_obs_heatmap.png"))

    customColours = generate_balanced_colours(obsVars)
    obs_heatmap = pheatmap(obsVars,
                           cluster_rows = hclust.comp,
                           cluster_col = FALSE,
                           color = customColours$colours,
                           breaks = customColours$breaks,
                           annotation_colors = annotations$colors,
                           annotation_row = annotations$ann,
                           fontsize_col = 8,
                           fontsize_row = 6,
                           width=9,
                           height = 14,
                           filename=paste0(plots_dir, "/all_calls_obs_vars_heatmap.png"))

    obs_heatmap = pheatmap(log10(obsVars),
                           cluster_rows = hclust.comp,
                           cluster_col = FALSE,
                           color = colorRampPalette(c("white", "navy"))(50),
                           annotation_colors = annotations$colors,
                           annotation_row = annotations$ann,
                           fontsize_col = 8,
                           fontsize_row = 6,
                           width=9,
                           height = 14,
                           filename=paste0(plots_dir, "/all_calls_obs_vars_log_heatmap.png"))

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

sample_info = read.table(snakemake@input[["sample"]], header=TRUE)

compare_versions(psm_results, base_run, dirname(snakemake@output[[1]]), sample_info)
