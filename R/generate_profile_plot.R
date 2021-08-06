suppressPackageStartupMessages({
  if (!'package:assertthat' %in% search()) library(assertthat)
  if (!'package:cowplot' %in% search()) library(cowplot)
  if (!'package:grid' %in% search()) library(grid)         #  for plotting multiple plots in one frame
  if (!'package:gridExtra' %in% search()) library(gridExtra)    #  for plotting multiple plots in one frame
  if (!'package:RColorBrewer' %in% search()) library(RColorBrewer)
  if (!'package:Seurat' %in% search()) library(Seurat)
})

plot_all_profiles <- function(sobj, out_prefix, plot_extras=FALSE, scatter_color='red', other_assay_names=c('IDX', 'ADT')) {
    plots <- list()
    if (scatter_color %in% colnames(sobj@meta.data)){
      plots[[1]] <- generate_profile_plot(sobj,
                                          feature1 = "nCount_RNA",
                                          feature2 = "percent.mt",
                                          feature1_binwidth=100,
                                          feature2_binwidth=0.1,
                                          scatter_color=scatter_color,
                                          plot_legend=TRUE)
    } else {
      plots[[1]] <- generate_profile_plot(sobj,
                                          feature1 = "nCount_RNA",
                                          feature2 = "percent.mt",
                                          feature1_binwidth=100,
                                          feature2_binwidth=0.1,
                                          scatter_color=scatter_color)      
    }
    plots[[2]] <- generate_profile_plot(sobj,
                                        feature1 = "nCount_RNA",
                                        feature2 = "percent.ribo",
                                        feature1_binwidth=100,
                                        feature2_binwidth=0.1,
                                        scatter_color=scatter_color)
    plots[[3]] <- generate_profile_plot(sobj,
                                        feature1 = "nCount_RNA",
                                        feature2 = "nFeature_RNA",
                                        feature1_binwidth=100,
                                        feature2_binwidth=100,
                                        scatter_color=scatter_color)

    plots[[4]] <- generate_profile_plot(sobj,
                                        feature1 = "percent.ribo",
                                        feature2 = "percent.mt",
                                        feature1_binwidth=0.1,
                                        feature2_binwidth=0.1,
                                        scatter_color=scatter_color)
    plots[[5]] <- generate_profile_plot(sobj,
                                        feature1 = "nFeature_RNA",
                                        feature2 = "percent.mt",
                                        feature1_binwidth=100,
                                        feature2_binwidth=0.1,
                                        scatter_color=scatter_color)
    plots[[6]] <- generate_profile_plot(sobj,
                                        feature1 = "percent.ribo",
                                        feature2 = "nFeature_RNA",
                                        feature1_binwidth=0.1,
                                        feature2_binwidth=100,
                                        scatter_color=scatter_color)
    df <- data.frame(
        cell_counts=seq(0, 1.01, 0.1)*dim(sobj@meta.data)[1],
        percent.mt=quantile(sobj@meta.data[["percent.mt"]], seq(0, 1.01, 0.1)),
        percent.ribo=quantile(sobj@meta.data[["percent.ribo"]], seq(0, 1.01, 0.1)),
        nFeature_RNA=quantile(sobj@meta.data[["nFeature_RNA"]], seq(0, 1.01, 0.1)),
        nCount_RNA=quantile(sobj@meta.data[["nCount_RNA"]], seq(0, 1.01, 0.1)),
        row.names=seq(0, 1.01, 0.1)
        )
    for (assay_name in other_assay_names){
      if (assay_name %in% names(sobj@assays)){
        plots[[length(plots)+1]] <- generate_profile_plot(sobj,
                                                          feature1 = paste0("nCount_", assay_name),
                                                          feature2 = "nCount_RNA",
                                                          feature1_binwidth=100,
                                                          feature2_binwidth=100,
                                                          scatter_color=scatter_color)
        plots[[length(plots)+1]] <- generate_profile_plot(sobj,
                                                          feature1 = paste0("nCount_", assay_name),
                                                          feature2 = "percent.mt",
                                                          feature1_binwidth=100,
                                                          feature2_binwidth=0.1,
                                                          scatter_color=scatter_color)
        df[paste0("nCount_", assay_name)] <- quantile(sobj@meta.data[[paste0("nCount_", assay_name)]], seq(0, 1.01, 0.1))
      }
    }

    write.table(format(df, digits=2), 
                file=paste0(out_prefix, ".tsv"), 
                row.names=T, 
                col.names=T, 
                quote=F, 
                sep="\t")

    pdf(paste0(out_prefix, ".pdf"), width = 21, height = 7 * length(plots)%/%3, onefile=TRUE)
    print(CombinePlots(plots = plots, ncol=3))
    if (plot_extras) {
      print(generate_inset_histogram(sobj=sobj, feature = "nCount_RNA",
                                     binwidth=100))
      print(generate_inset_histogram(sobj=sobj, feature = "nFeature_RNA",
                                     binwidth=100))
      for (assay_name in other_assay_names){
        if (assay_name %in% names(sobj@assays)){
          print(generate_inset_histogram(sobj=sobj, feature = paste0("nCount_", assay_name),
                                         binwidth=100))
        }
      }
      print(generate_inset_histogram(sobj=sobj, feature = "percent.mt",
                                     binwidth=0.1,
                                     inset_cutoff_percentile = 0.9))
      print(generate_inset_histogram(sobj=sobj, feature = "percent.ribo",
                                     binwidth=0.1,
                                     inset_cutoff_percentile = 0.9))
    }
    dev.off()
}


generate_profile_plot <- function(sobj, feature1, feature2, feature1_binwidth=100, feature2_binwidth=100, 
                                  visual_outlier_cutoff1=0.999, visual_outlier_cutoff2=0.999, scatter_color='red',
                                  plot_legend=FALSE) {
  suppmsg <- assert_that(feature1 %in% colnames(sobj@meta.data), msg=paste0(feature1, ' was not present in the metadata of sobj'))
  suppmsg <- assert_that(feature2 %in% colnames(sobj@meta.data), msg=paste0(feature2, ' was not present in the metadata of sobj'))
  suppmsg <- assert_that(0 < visual_outlier_cutoff1 && visual_outlier_cutoff1 <=1.0, msg='visual_outlier_cutoff1 must be in the range (0,1]')
  suppmsg <- assert_that(0 < visual_outlier_cutoff2 && visual_outlier_cutoff2 <=1.0, msg='visual_outlier_cutoff2 must be in the range (0,1]')
  
  lay <- rbind(c(1,  1,2,2,2,2),
               c(1,  1,2,2,2,2),
               c(1,  1,2,2,2,2),
               c(NA,NA,3,3,3,3),
               c(NA,NA,3,3,3,3))
  
  lims = as.vector(
    c(quantile(sobj@meta.data[[feature1]], visual_outlier_cutoff1), 
      quantile(sobj@meta.data[[feature2]], visual_outlier_cutoff2)))
  xticks <- as.vector(quantile(sobj@meta.data[[feature1]], seq(0, max(0.9, visual_outlier_cutoff1), 0.1)))
  if (xticks[length(xticks)] != lims[1]) {
    xticks <- c(xticks, lims[1])
  }
  yticks <- as.vector(quantile(sobj@meta.data[[feature2]], seq(0, max(0.9, visual_outlier_cutoff2), 0.1)))
  if (yticks[length(yticks)] != lims[2]) {
    yticks <- c(yticks, lims[2])
  }

  main <- ggplot(sobj@meta.data, aes_string(x=feature1, y=feature2)) 
   
  if (scatter_color %in% colnames(sobj@meta.data)) {
    custom_colors = brewer.pal(9, 'Set1')[c(1:2, 9, 3:8)]
    main <- main +
      geom_point(aes_string(col=scatter_color, shape=scatter_color), size=0.5) + 
      scale_color_manual(values = custom_colors) +
      guides(colour = guide_legend(override.aes = list(size=3)))
  } else {
    main <- main +
      geom_point(aes(col=scatter_color), size=0.5)
  }

  main <- main +
    xlim(NA, lims[1]) +
    ylim(NA, lims[2]) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()
    )

  if (plot_legend){
    main <- main +
      theme(legend.position="top")
  } else {
    main <- main +
      NoLegend()
  }
  
  y_hist <- ggplot(sobj@meta.data, aes_string(x=feature2)) + 
    geom_histogram(aes(col='red', fill='red'), binwidth=feature2_binwidth) +
    theme_bw() +
    theme(axis.title.x=element_blank(), 
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank()) +
    scale_x_continuous(limits=c(NA, lims[2]),
                       sec.axis = sec_axis(trans = ~., 
                                           breaks=yticks, labels=NULL)) +
    coord_flip() + 
    scale_y_reverse()
  
  if (plot_legend){
    y_hist <- y_hist +
      guides(colour = guide_legend(override.aes = list(colour='white', fill='white')),
             fill = guide_legend(override.aes = list(fill='white'))) +
      theme(legend.position="top", 
            legend.title = element_blank(), 
            legend.key=element_blank(), 
            legend.text = element_blank()
      )
  } else {
    y_hist <- y_hist +
      NoLegend()
  }
  
  
  x_hist <- ggplot(sobj@meta.data, aes_string(x=feature1)) + 
    geom_histogram(aes(col='red', fill='red'), binwidth=feature1_binwidth) +
    theme_bw() +
    theme(axis.title.y=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) +
    scale_x_continuous(limits=c(NA, lims[1]),
                       sec.axis = sec_axis(trans = ~., 
                                           breaks=xticks, labels=NULL)) +
    NoLegend() +
    scale_y_reverse()
  empty_plot <- ggplot() + theme_void()
  grid.arrange(grobs=list(y_hist, main, x_hist), layout_matrix = lay)
}


generate_inset_histogram <- function(sobj, feature, binwidth=100, inset_cutoff_percentile=0.5, visual_outlier_cutoff=0.999) {
  suppmsg <- assert_that(feature %in% colnames(sobj@meta.data), msg=paste0(feature, ' was not present in the metadata of sobj'))
  suppmsg <- assert_that(0 < inset_cutoff_percentile && inset_cutoff_percentile <=1.0, msg='inset_cutoff_percentile must be in the range (0,1]')
  suppmsg <- assert_that(0 < visual_outlier_cutoff && visual_outlier_cutoff <=1.0, msg='visual_outlier_cutoff must be in the range (0,1]')
  
  inset_lim <- quantile(sobj@meta.data[[feature]], inset_cutoff_percentile)[[1]]
  visual_outlier_lim <- quantile(sobj@meta.data[[feature]], visual_outlier_cutoff)[[1]]
  
  sec_xticks <- as.vector(quantile(sobj@meta.data[[feature]], seq(0, max(0.9, visual_outlier_cutoff), 0.1)))
  if (sec_xticks[length(sec_xticks)] != visual_outlier_lim){
    sec_xticks <- c(sec_xticks, visual_outlier_lim)
  }
  
  main_plot <- ggplot(sobj@meta.data, aes_string(x=feature)) + 
               geom_histogram(aes(col='red', fill='red'), binwidth=binwidth) +
               theme_bw()
  inset_hist <- main_plot +
                geom_vline(xintercept=inset_lim, lty=2) +
                theme(axis.title.x=element_blank(),
                      axis.text.x = element_text(angle = 45, hjust = 1),
                      axis.title.y=element_blank()) +
                scale_x_continuous(limits=c(NA, visual_outlier_lim),
                                   sec.axis=sec_axis(trans=~., breaks=sec_xticks, labels=NULL)) +
                NoLegend()
  outer_hist <- main_plot + 
                scale_x_continuous(limits=c(NA, inset_lim)) +
                NoLegend()

  p <- ggdraw() +
       draw_plot(outer_hist) +
       draw_plot(inset_hist, x = 0.55, y = .55, width = .4, height = .4)
  return(p)
}
