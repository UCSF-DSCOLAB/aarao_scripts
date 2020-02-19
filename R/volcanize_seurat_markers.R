library(ggrepel)

volcanize_from_FindAllMarkers <- function(x, pos_group, label_n=10, logFC_threshold=0.4) {
  x <- transform(x, 
                 plmin_log_FC=ifelse(cluster==pos_group, avg_logFC, -avg_logFC)
  )
  top_genes <- (x %>% group_by(cluster) %>% filter(avg_logFC > logFC_threshold) %>% top_n(label_n, wt=avg_logFC))$gene
  x <- transform(x, col=factor(ifelse(gene %in% top_genes, 2, ifelse(avg_logFC>logFC_threshold, 1, 0)), levels=c(0,1,2)))  
  
  plot <- ggplot(x, aes(plmin_log_FC, -log10(p_val))) + 
    geom_point(aes(col=col)) + 
    scale_color_manual(values=c("black", "red", "red"), drop=F) + 
    geom_text_repel(data=filter(x, col=='2'), aes(label=gene), size=5, color='red') +
    theme_minimal() +
    theme(legend.position = "none") +
    xlab('logFC') +
    xlim(-max(x$avg_logFC),max(x$avg_logFC)) +
    ggtitle(label = paste0('Upregulated in ', pos_group)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  if (any(x$p_val == 0)) {
    pbuild <- ggplot_build(plot = plot)
    y.min <- pbuild$layout$panel_params[[1]]$y.range[1]
    y.max <- pbuild$layout$panel_params[[1]]$y.range[2]
    
    plot <- plot + scale_y_continuous(breaks=c(pbuild$layout$panel_params[[1]]$y.major_source, y.max),
                                      labels=c(pbuild$layout$panel_params[[1]]$y.labels, expression(infinity))) +
      theme(axis.text.y = element_text(size = c(rep(15, length(pbuild$layout$panel_params[[1]]$y.major_source)), 25)),
            axis.text.x = element_text(size = rep(15, length(pbuild$layout$panel_params[[1]]$y.major_source))),
            plot.margin=unit(c(30, 30, 5.5, 5.5), "points"))
  }
  return(plot)
}

volcanize_from_FindMarkers <- function(x, pos_group='Positive_group', label_n=10, logFC_threshold=0.4) {
  x <- transform(x, 
                 regulation=ifelse(avg_logFC<0, 'DOWN', 'UP'),
                 abs_log_FC = abs(avg_logFC)
  )
  x$gene <- rownames(x)
  
  top_genes <- (x %>% group_by(regulation) %>% filter(abs_log_FC > logFC_threshold) %>% top_n(label_n, wt=abs_log_FC))$gene
  
  x <- transform(x, col=factor(ifelse(gene %in% top_genes, 2, ifelse(abs_log_FC>logFC_threshold, 1, 0)), levels=c(0,1,2)))  
  
  plot <- ggplot(x, aes(avg_logFC, -log10(p_val))) + 
    geom_point(aes(col=col)) + 
    scale_color_manual(values=c("black", "red", "red"), drop=F) + 
    geom_text_repel(data=filter(x, col=='2'), aes(label=gene), size=5, color='red') +
    theme_minimal() +
    theme(legend.position = "none") +
    xlab('logFC') +
    xlim(-max(x$avg_logFC),max(x$avg_logFC)) +
    ggtitle(label = paste0('Upregulated in ', pos_group)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  if (any(x$p_val == 0)) {
    pbuild <- ggplot_build(plot = plot)
    y.min <- pbuild$layout$panel_params[[1]]$y.range[1]
    y.max <- pbuild$layout$panel_params[[1]]$y.range[2]
    
    plot <- plot + scale_y_continuous(breaks=c(pbuild$layout$panel_params[[1]]$y.major_source, y.max),
                                      labels=c(pbuild$layout$panel_params[[1]]$y.labels, expression(infinity))) +
      theme(axis.text.y = element_text(size = c(rep(15, length(pbuild$layout$panel_params[[1]]$y.major_source)), 25)),
            axis.text.x = element_text(size = rep(15, length(pbuild$layout$panel_params[[1]]$y.major_source))),
            plot.margin=unit(c(30, 30, 5.5, 5.5), "points"))
  }
  return(plot)
}
