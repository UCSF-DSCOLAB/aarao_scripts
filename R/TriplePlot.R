suppressPackageStartupMessages({
  if (!'package:Seurat' %in% search()) library(Seurat)
  if (!'package:ggplot2' %in% search()) library(ggplot2)
  if (!'package:gridExtra' %in% search()) library(gridExtra)
})

TriPlot <- function(sobj, features, reduction.use="umap", group.by="seurat_clusters", jitter=TRUE) {
  if (jitter){
    pt.size = 1
  } else {
    pt.size = 0
  }
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
                                 group.by=group.by,
                                 pt.size = pt.size) + NoLegend()
    )
    plots[[f]] <- grid.arrange(grobs=tmp_plots, layout_matrix=layout)
  }
  plots
}
