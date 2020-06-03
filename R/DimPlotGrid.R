suppressPackageStartupMessages({
  if (!'package:cowplot' %in% search()) library(cowplot)
  if (!'package:ggplot2' %in% search()) library(ggplot2)
  if (!'package:Seurat' %in% search()) library(Seurat)
})

DimPlotGrid <- function(sobj, split.by, prefix=NULL, ncol=5, cluster_col='red', other_col='lightgrey') {
    ## A function to plot a grid of dimplots where each member in the grid is one cluster colored by 
    ## `cluster_col` and the rest are colored `other_col`. The top left of the grid is the labeled DimPlot
    ## 
    ##
    ## Parameters:
    ##     sobj : The Seurat object (REQUIRED)
    ##     split.by : The metadata column to split the grid on (REQUIRED)
    ##     prefix : A character prefix for the clusters (DEFAULT=NULL)
    ##     ncol : The number of columns in the final grid (DEFAULT=5)
    ##     cluster_col : The color for the highlighted cluster per grid point  (DEFAULT=red)
    ##     other_col : The color for the background cells per grid point (DEFAULT=lightgrey)
    ##
    ## Return value
    ##     The result of running `plot_grid` on all plots with ncol=`ncol`
    ##    
    clusters <- as.vector(unique(sobj@meta.data[[split.by]]))
    clusters <- clusters[order(as.numeric(clusters))]
    cluster_cells <- lapply(clusters, 
                            function(x){rownames(sobj@meta.data)[sobj@meta.data[[split.by]]==x]})
    if (!is.null(prefix)){
        names(cluster_cells) <- paste(prefix, clusters, sep='')
    }
    
    plots <- list()
    plots[['Dimplot']] <- DimPlot(sobj, group.by=split.by, label=T) + NoLegend()
    for (cluster_name in names(cluster_cells)) {
        plots[[cluster_name]] <- DimPlot(sobj, cells.highlight = cluster_cells[[cluster_name]]) + 
                                        scale_color_manual(labels=c('other', cluster_name), 
                                                           values=c(other_col, cluster_col)) +
                                        NoLegend() +
                                        theme(axis.title.x=element_blank(),
                                              axis.text.x=element_blank(),
                                              axis.ticks.x=element_blank(),
                                              axis.title.y=element_blank(),
                                              axis.text.y=element_blank(),
                                              axis.ticks.y=element_blank()) +
                                        ggtitle(cluster_name)
    }
    plot_grid(plotlist=plots, ncol=ncol)
}