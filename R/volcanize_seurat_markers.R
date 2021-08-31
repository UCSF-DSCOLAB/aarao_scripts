suppressPackageStartupMessages({
  if (!'package:assertthat' %in% search()) library(assertthat)
  if (!'package:dplyr' %in% search()) library(dplyr)
  if (!'package:ggplot2' %in% search()) library(ggplot2)
  if (!'package:ggrepel' %in% search()) library(ggrepel)
})

# A function to Plot a volcano of a "test" cluster vs an "other" one, given the results 
# of Seurat FindAllmarkers. "test" vs "other" means that genes with positive Fold Change
# in the plot are UPREGULATED in "test".
#
# For best results and an actual volcano, run with logfc.threshold=0.001
# Parameters are
# df                 : The input data frame created using FindAllMarkers
# pos_cluster        : The name within the `cluster` column to use as the "test"
# neg_cluster        : The name within the `cluster` column to use as the "other"
# label_n            : The number of significant hits to label on both sides (default=10)
# logFC_threshold    : The logFC threshold above which we will color the point red (default=0.4)
# logFC_colname      : The column name in the df for logFC (default='avg_log2FC')
# p_threshold        : The threshold for significance (default=0.05)
# p_colname          : The column name in the df for p (or p adjusted) value (default='p_val_adj')
# logp_visual_cutoff : Values "higher" than this value will be set to infinity for visual clarity (default=NULL)
# sig_color          : Color to use for significant hits i.e. lFC > logFC_threshold and p < p_threshold (default=red')
# background_color   : Color to use for the background points (default='lightgrey')
# text_color         : Color to use for labeled genes (default='black')
# point_size         : Point size for ggplot2 (default=0.5)
volcanize_from_FindAllMarkers <- function(df, pos_cluster, neg_cluster, label_n=10, 
                                          logFC_threshold=0.4, logFC_colname='avg_log2FC', 
                                          p_threshold=0.05, p_colname='p_val_adj', logp_visual_cutoff=NULL,
                                          sig_color='red',background_color='lightgrey',  
                                          text_color='black', point_size=0.5) {
  
  p_visual_cutoff = ifelse(is.null(logp_visual_cutoff), 0, 10**(-logp_visual_cutoff))

  temp <- assert_that(all(sign(df[[logFC_colname]])==1), 
                   msg = paste0("Not all LogFC values in ", logFC_colname, "are positive. Rerun `FindAllMarkers` with `only.pos=TRUE`"))
  x <- df %>% 
          select(!!sym(logFC_colname), !!sym(p_colname), gene, cluster) %>%
          rename(logFC = !!sym(logFC_colname),
                 p_val = !!sym(p_colname)) %>% 
          filter(cluster %in% c(pos_cluster, neg_cluster)) %>%
          mutate(plmin_log_FC=ifelse(cluster==pos_cluster, 
                                     logFC, 
                                     -logFC))
  top_genes <- x %>% 
                  group_by(cluster) %>% 
                  filter(logFC >= logFC_threshold & p_val <= p_threshold) %>% 
                  top_n(label_n, wt=logFC) %>%
                  pull(gene)
  
  x <- x %>% 
        mutate(col=factor(ifelse(gene %in% top_genes, 
                                    2, 
                                    ifelse(logFC>logFC_threshold & p_val <= p_threshold,
                                           1, 
                                           0)), 
                             levels=c(0,1,2)))
  
  # Now conduct the visual p pvalue adjustment
  x <- x %>% 
    mutate(p_val = ifelse(p_val <= p_visual_cutoff, 0, p_val))
  
  plot <- ggplot(x, aes(x=plmin_log_FC, y=-log10(p_val))) + 
              geom_point(aes(col=col), size=point_size) + 
              scale_color_manual(values=c(background_color, 
                                          sig_color, 
                                          sig_color), drop=F) + 
              geom_text_repel(data=filter(x, col=='2'), 
                              aes(label=gene), 
                              size=5, 
                              color=text_color) +
              theme_minimal() +
              theme(legend.position = "none") +
              xlab(logFC_colname) +
              ylab(paste0("-log10(", p_colname,")")) +
              xlim(-max(x$logFC),
                   max(x$logFC)) +
              ggtitle(label = paste0('Volcano of ', pos_cluster, ' vs. ', neg_cluster)) +
              theme(plot.title = element_text(hjust = 0.5))
      
  if (any(x$p_val == 0)) {
    pbuild <- ggplot_build(plot = plot)
    y.min <- pbuild$layout$panel_params[[1]]$y.range[1]
    y.max <- pbuild$layout$panel_params[[1]]$y.range[2]
    
    plot <- plot + scale_y_continuous(breaks=c(pbuild$layout$panel_params[[1]]$y$minor_breaks, y.max),
                                      labels=c(pbuild$layout$panel_params[[1]]$y$minor_breaks, expression(infinity))) +
      theme(axis.text.y = element_text(size = 15),
            axis.text.x = element_text(size = 15),
            plot.margin=unit(c(30, 30, 5.5, 5.5), "points"))
  }
  return(plot)
}


# This literally just modifies the input dataframe to look like the results of FindAllMarkers 
# and calls the above function
volcanize_from_FindMarkers <- function(x, pos_group_name='Positive_group', neg_group_name='Negative_group', ...) {
    x <- df %>%
            mutate(cluster=ifelse(!!sym(logFC_colname) < 0, neg_group_name, pos_group_name),
                   !!sym(logFC_colname) := abs(!!sym(logFC_colname)))
    return(volcanize_from_FindAllMarkers(x, pos_group=pos_group_name, neg_group=neg_group_name, ...))
}  
