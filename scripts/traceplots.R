library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

read_numeric=function(dirs, filename) {
  files=file.path(dirs,filename)
  use=file.exists(files)
  values=lapply(which(use), function(i) {
    message(dirs[i])
    scan(files[i],what="", sep=",") %>% as.numeric()
  })
  nmin=sapply(values,length) %>% min()
  values %<>% lapply(., "[", 1:nmin) %>% do.call("cbind",.)
  colnames(values)=dirs[use]
  values
}

get_only_legend <- function(plot) {
  plot_table <- ggplot_gtable(ggplot_build(plot))
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")
  legend <- plot_table$grobs[[legend_plot]]
  return(legend)
}

calculate_quantiles <- function(datasets, filename, block_size=100) {
  data_mat = read_numeric(datasets, filename)
  df = data.frame(data_mat)
  colnames(df) = paste0("seed", 1:ncol(df))
  df["block"] = rep(1:(nrow(df) / block_size), each=block_size)
  df = df%>% tidyr::pivot_longer(cols=1:(ncol(df) - 1))

  quantiled = df %>%
    group_by(block, name) %>%
    summarise(value = quantile(value, c(0.1, 0.25, 0.5, 0.75, 0.9)),
              quantile = c(0.1, 0.25, 0.5, 0.75, 0.9),
              count = n())

  write.table(quantiled, "quantiled.csv")
  return(quantiled)
}

median_traceplot <- function(label, quantiled) {

  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  if (length(unique(quantiled$name)) > 8) {
      # Using RColorBrewer::brewer.pal(10, "Spectral")
      cbbPalette <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")
  }

  g = ggplot(quantiled, aes(x=block, y=value, colour=name)) +
    geom_vline(xintercept=max(quantiled["block"])/ 2, colour="red") +
    scale_color_manual(values=cbbPalette) +
    labs(y=label) +
    geom_line(data=subset(quantiled,quantile==0.5), alpha=0.5)

  return(g)
}

quantile_traceplot <- function(label, quantiled) {

  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  if (length(unique(quantiled$name)) > 8) {
      # Using RColorBrewer::brewer.pal(10, "Spectral")
      cbbPalette <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")
  }

  g = ggplot(quantiled, aes(x=block, y=value, colour=name)) +
    geom_vline(xintercept=max(quantiled["block"])/ 2, colour="red") +
    scale_color_manual(values=cbbPalette) +
    labs(y=label) +
    geom_line(data=subset(quantiled,quantile==0.5)) +
    geom_line(data=subset(quantiled,quantile==0.1), alpha = 0.2, linetype='dashed') +
    geom_line(data=subset(quantiled,quantile==0.9), alpha = 0.2, linetype='dashed') +
    geom_line(data=subset(quantiled,quantile==0.25), alpha = 0.2) +
    geom_line(data=subset(quantiled,quantile==0.75), alpha = 0.2)

  return(g)
}

quantile_traceplots_dataset <- function(dataset, datasets, block_size=100) {
  quant_alpha = calculate_quantiles(datasets, "alpha.csv")
  quant_K = calculate_quantiles(datasets, "K.csv")
  quant_latent = calculate_quantiles(datasets, "pLatentsGivenClusters.csv")

  g_alpha = median_traceplot("alpha", quant_alpha)
  g_K = median_traceplot("K", quant_K)
  g_latent = median_traceplot("pLatents", quant_latent)

  g = grid.arrange(g_alpha + theme(legend.position = "none"),
                   g_K + theme(legend.position = "none"),
                   g_latent + theme(legend.position = "none"),
                   get_only_legend(g_latent),
                   top = textGrob(paste0("Median traceplot - ", dataset),gp=gpar(fontsize=20,font=3)),
                   ncol=4,
                   widths=c(2,2,2,1))
  ggsave(paste0("plots/", dataset, "/trace_medians.png"), g, width=12, height=7, units="in")

  g_alpha = quantile_traceplot("alpha", quant_alpha)
  g_K = quantile_traceplot("K", quant_K)
  g_latent = quantile_traceplot("pLatents", quant_latent)

  g = grid.arrange(g_alpha + theme(legend.position = "none"),
                   g_K + theme(legend.position = "none"),
                   g_latent + theme(legend.position = "none"),
                   get_only_legend(g_latent),
                   top = textGrob(paste0("Quantile traceplot - ", dataset),gp=gpar(fontsize=20,font=3)),
                   ncol=4,
                   widths=c(2,2,2,1))
  ggsave(paste0("plots/", dataset, "/trace_quantiles.png"), g, width=12, height=7, units="in")

  return(g)
}

dataset_name = snakemake@wildcards[["dataset"]]
dataset_folders = snakemake@params[["dataset_folders"]]
print(dataset_folders)
print(class(dataset_folders))
dataset_folders = as.character(dataset_folders)
print(class(dataset_folders))
quantile_traceplots_dataset(dataset_name, dataset_folders)
