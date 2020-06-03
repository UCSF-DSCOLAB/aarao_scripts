suppressPackageStartupMessages({
  if (!'package:cowplot' %in% search()) library(cowplot)
  if (!'package:ggplot2' %in% search()) library(ggplot2)
  if (!'package:RColorBrewer' %in% search()) library(RColorBrewer)
  if (!'package:reshape2' %in% search()) library(reshape2)
  if (!'package:Seurat' %in% search()) library(Seurat)
})

plot_cluster_diffs_by_sample <- function(sobj, ann.col='SAMPLE.by.SNPs', 
                                         ident.col='seurat_clusters', palette='Dark2', 
                                         ann.col.prefix='S', ident.col.prefix='C') {
  ## A function to plot and tabulate  the fractions of cells belonging to `ann.col` in 
  ## `ident.col` and vice versa. 
  ## 
  ##
  ## Parameters:
  ##     sobj : The Seurat object (REQUIRED)
  ##     ann.col : The sample annotations of interest. Must exist in the metadata of the object
  ##                  (DEFAULT=SAMPLE.by.SNPs), 
  ##     ident.col : The identities for the object. Must exist in the metadata of the object
  ##                   (DEFAULT=seurat_clusters)
  ##     palette : Custom color paletted for output plots. CURRENTLY NOT IMPLEMENTED. (DEFAULT=Dark2)
  ##     ann.col.prefix : character prefix for ann.col if any sample name is strictly   
  ##                         non-numeric (DEFAULT=S)
  ##     ident.col.prefix : character prefix for ident.col if any cluster name is strictly   
  ##                         non-numeric (DEFAULT=C)
  ##
  ## Return value
  ##     A list with the following elements
  ##         plots : A list of the following plots
  ##             dimplot:
  ##             ann_fractions_per_ident_hist : A histogram of
  ##             ident_fractions_per_ann_hist : A histogram of 
  ##             ident_fractions_per_ann_line : A line plot of         
  ##         fraction_table: A tibble with counts and percentages of idents per ann and vice 
  ##                         versa
  ##

  annotations <- levels(as.factor(sobj@meta.data[[ann.col]]))
  suppressWarnings({
    if (any(!is.na(as.numeric(annotations)))){
      cat(paste0('Not all annotations are strictly non-numeric. Prefixing with `', ann.col.prefix, 
                 '`.\n'))
      sobj@meta.data[[ann.col]] <- gsub('^', ann.col.prefix, 
                                        as.vector(sobj@meta.data[[ann.col]]))
      annotations <- gsub('^', ann.col.prefix, as.vector(annotations))
    }
  })

  identities <- levels(as.factor(sobj@meta.data[[ident.col]]))
  suppressWarnings({
    if (any(!is.na(as.numeric(identities)))){
      cat(paste0('Not all identities are strictly non-numeric. Prefixing with `', ident.col.prefix, 
                 '`.\n'))
      sobj@meta.data[[ident.col]] <- gsub('^', ident.col.prefix, 
                                          as.vector(sobj@meta.data[[ident.col]]))
      identities <- gsub('^', ident.col.prefix, as.vector(identities))
    }
  })

  # TODO Add custom palettes
  colors <- palette
  #try(colors <- brewer.pal(name = palette, n=length(annotations)), silent = TRUE)
  
  #if (is.null(colors)){
  #  stop("Palette `", palette, "` supports only (", brewer.pal.info[palette, 'maxcolors'], ") colors but you have (", length(annotations), ") annotations. please choose a different palette.")
  #}
  
  #names(colors) <- annotations
  outs <- list()
  outs[['plots']] <- list()
  outs[['plots']][['dimplot']] <- DimPlot(sobj, split.by = ann.col) + 
                                      guides(colour = guide_legend(override.aes = list(colour='white'))) +
                                      theme(legend.title = element_blank(), 
                                            legend.key=element_blank(), 
                                            legend.text = element_blank())
  
  outs[['fraction_table']] <- sobj@meta.data %>% 
      count((!!sym(ident.col)), (!!sym(ann.col))) %>%
      group_by((!!sym(ident.col))) %>%
      mutate(ident_size=sum(n)) %>% 
      mutate(fraction_of_ident = n/ident_size) %>%
      group_by((!!sym(ann.col))) %>%
      mutate(ann_size=sum(n)) %>% 
      mutate(fraction_of_ann = n/ann_size) %>%
      arrange(-ident_size, .by_group = TRUE) %>%
      mutate(idx=row_number())
  outs[['plots']][['ann_fractions_per_ident_hist']] <- ggplot(outs[['fraction_table']], 
                                                                   aes_string(y="fraction_of_ident", 
                                                                              x=ident.col,
                                                                              fill=ann.col)) +
                                                              geom_col(position = "dodge") +
                                                              theme(legend.key.size = unit(0.3, "cm"))
  
  # convenience
  temp <- outs[['fraction_table']]
  outs[['plots']][['ident_fractions_per_ann_hist']] <- ggplot(temp, aes_string(y="fraction_of_ann", x=ann.col, fill=ident.col)) +
                                                              geom_col(position = "dodge") +
                                                              theme(legend.key.size = unit(0.3, "cm"), 
                                                                    axis.text.x=element_text(angle=90, hjust = 1))

  outs[['plots']][['ident_fractions_per_ann_line']] <- ggplot(temp, aes_string(x="idx", y="fraction_of_ann", col=ann.col)) +
                                                              geom_line() + 
                                                              geom_point(shape=1, size=3) +
                                                              #scale_color_brewer(palette=palette)
                                                              xlab(ident.col) +
                                                              ylab('percent') +
                                                              scale_x_continuous(breaks=temp[temp[[ann.col]]==temp[[ann.col]][1], 'idx', drop=T], 
                                                                                 labels=as.vector(temp[temp[[ann.col]]==temp[[ann.col]][1], ident.col, drop=T]))
  outs
}