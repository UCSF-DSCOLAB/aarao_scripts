library(Seurat)
library(ggplot2)
library(gridExtra)

TriPlot <- function(sobj, features, reduction.use="umap", group.by="seurat_clusters") {
  plots <- list()
  layout = rbind(c(1,1,2,2),
                 c(1,1,2,2),
                 c(3,3,3,3),
                 c(3,3,3,3))
  for (f in features){
    tmp_plots <- list(p1=DimPlot(sobj,
                                 group.by=group.by, 
                                 reduction=reduction.use, 
                                 label=TRUE) + NoLegend(),
                      p2=FeaturePlot(sobj,
                                     features=f,
                                     reduction=reduction.use),
                      p3=VlnPlot(sobj, 
                                 features=f,
                                 group.by=group.by) + NoLegend()
    )
    plots[[f]] <- grid.arrange(grobs=tmp_plots, layout_matrix=layout)
  }
  plots
}
