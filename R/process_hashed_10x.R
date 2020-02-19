library(assertthat)
library(Seurat)
library(dplyr)  # inline modification of matrices
library(cowplot)  # Pretty plots
library(ggplot2)  # Pretty plots
library(grid)  #  for plotting multiple plots in one frame
library(gridExtra)  #  for plotting multiple plots in one frame
library(scales)  # to access break formatting functions
 
# Set your working directory
setwd('/data')
source('/scripts/identify_hto_clusters.R')
source('/scripts/demux_HTOs_by_inflexions.R')
source('/scripts/generate_profile_plot.R')
source('/scripts/volcanize_seurat_markers.R')

pool_map <- read.table('pool_map.tsv', header=T, sep='\t', stringsAsFactors = F)

robjs <- list()
for (pool in unique(pool_map$pool)){
  robj1 <- paste0(pool, '_scTransformed_processed.Robj')
  robj2 <- paste0(pool, '_scTransformed.Robj')
  if (file.exists(robj1)) {
    robjs[[pool]] <- robj1
  } else if(file.exists(robj2)) {
    robjs[[pool]] <- robj2
  } else{
    print(paste0('Could not find ', robj1, ' or ', robj2, '.'))
    quit(save = 'no', status = 1)
  }
}

sobjs <- list()
for (pool in names(robjs)){
  load(robjs[[pool]])
  sobjs[[pool]] <- get(pool)
  rm(list=pool)
}

if (file.exists('inflexions.tsv')){
  inflexions_table <- read.table('inflexions.tsv', sep='\t', header=T, stringsAsFactors=F)
} else {
  samples_to_process <- names(sobjs)
  inflexions <- list()
}

for (pool in samples_to_process){
  assert_that(sum(!as.vector(pool_map[pool_map$pool == pool, 'tag']) %in% rownames(sobjs[[pool]]@assays$Hashes))==0, msg=paste0('Tags reported for ', pool, ' were not found in the sobj.'))
}

results <- list()
for (pool in samples_to_process){
  sample_names <- as.vector(pool_map[pool_map$pool == pool, 'sample_name'])
  names(sample_names) <- as.vector(pool_map[pool_map$pool == pool, 'tag'])
  
  if (!all(sort(as.vector(pool_map[pool_map$pool == pool, 'tag']))==sort(rownames(sobjs[[pool]]@assays$Hashes)))) {
    keep_features = c(rownames(sobjs[[pool]]@assays$RNA), as.vector(pool_map[pool_map$pool == pool, 'tag']))
    sobjs[[pool]] <- subset(sobjs[[pool]], features = keep_features)
  }
  
  if (is.null(inflexions[[pool]])){
    results[[pool]] <- demux_by_inflexions(
      sobj=sobjs[[pool]],
      sample_names=sample_names
    )
  } else {
    results[[pool]] <- demux_by_inflexions(
      sobj=sobjs[[pool]],
      sample_names=sample_names, 
      inflexions = inflexions[[pool]]
    )
  }
  
  nplots <- length(results[[pool]]$background_plots)
  
  png(paste0(pool, '_raw_hto_pairplots.png'), width = nplots*750, height = nplots*750, units = 'px')
  print(plot_grid(plotlist=unlist(results[[pool]][['raw_plots']], recursive = FALSE), ncol=nplots))
  dev.off()
  png(paste0(pool, '_scaled_hto_pairplots.png'), width = nplots*750, height = nplots*750, units = 'px')
  print(plot_grid(plotlist=unlist(results[[pool]][['scaled_plots']], recursive = FALSE), ncol=nplots))
  dev.off()
  png(paste0(pool, '_hto_background_profiles.png'), width = 1500, height = 750*max(1, floor(nplots/2)), units = 'px')
  print(plot_grid(plotlist=results[[`pool`]][['background_plots']], ncol=2))
  dev.off()
  inflexions[[pool]] <- results[[pool]]$inflexions
}

for (pool in names(results)){
  sobjs[[pool]]@meta.data$sample_name <- results[[pool]]$cell_ids

  sample_names <- as.vector(pool_map[pool_map$pool == pool, 'sample_name'])
  n = length(sample_names)
  
  Idents(sobjs[[pool]]) <- sobjs[[pool]]@meta.data$sample_name
  png(paste0(pool, '_ridgeplot_inhouse.png'), width = 1500, height = ceiling(n/2)*750, units = 'px')
  print(RidgePlot(sobjs[[pool]], assay = "Hashes", features = rownames(sobjs[[pool]][["Hashes"]]), ncol = 2))
  dev.off()
  
  sobjs[[pool]] <- HTODemux(sobjs[[pool]], assay = "Hashes", positive.quantile = 0.99)
  png(paste0(pool, '_ridgeplot_HTODemux.png'), width = 1500, height = ceiling(n/2)*750, units = 'px')
  print(RidgePlot(sobjs[[pool]], assay = "Hashes", features = rownames(sobjs[[pool]][["Hashes"]]), ncol = 2))
  dev.off()
  if('seurat_clusters' %in% colnames(sobjs[[pool]]@meta.data)) {
    Idents(sobjs[[pool]]) <- sobjs[[pool]]@meta.data$seurat_clusters  
  } else {
    Idents(sobjs[[pool]]) <- sobjs[[pool]]@meta.data$sample_name
  }
  assign(pool, sobjs[[pool]])
  save(list=pool, file=gsub('.Robj', '_HTODemuxed.Robj', robjs[[pool]]))
  rm(list=pool)
}


