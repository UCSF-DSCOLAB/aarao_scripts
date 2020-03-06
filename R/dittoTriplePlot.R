library(dittoSeq)
library(gridExtra)

dittoTriPlot <- function(sobj, features, reduction.use="umap", group.by="seurat_clusters") {
  plots <- list()
  layout = rbind(c(1,1,2,2),
                 c(1,1,2,2),
                 c(3,3,3,3),
                 c(3,3,3,3))
  for (f in features){
    tmp_plots <- list(p1=dittoDimPlot(group.by, 
                                      sobj, 
                                      reduction.use=reduction.use, 
                                      labels.repel=T, 
                                      do.label=F, 
                                      legend.show=F),
                      p2=dittoDimPlot(f, 
                                      sobj, 
                                      reduction.use=reduction.use),
                      p3=dittoBoxPlot(f, 
                                      sobj, 
                                      group.by=group.by, 
                                      legend.show = F)
    )
    plots[[f]] <- grid.arrange(grobs=tmp_plots, layout_matrix=layout)
  }
  plots
}