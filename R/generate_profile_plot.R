library(assertthat)
generate_profile_plot <- function(sobj, feature1, feature2, feature1_binwidth=100, feature2_binwidth=100, visual_outlier_cutoff1=0.999, visual_outlier_cutoff2=0.999) {
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

  main <- ggplot(sobj@meta.data, aes_string(x=feature1, y=feature2)) + 
    geom_point(aes(col='red'), size=0.5) + 
    xlim(NA, lims[1]) +
    ylim(NA, lims[2]) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()
    ) +
    NoLegend()
  
  y_hist <- ggplot(sobj@meta.data, aes_string(x=feature2)) + 
    geom_histogram(aes(col='red', fill='red'), binwidth=feature2_binwidth) +
    theme_bw() +
    theme(axis.title.x=element_blank(), 
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank()) +
    scale_x_continuous(limits=c(NA, lims[2]),
                       sec.axis = sec_axis(trans = ~., 
                                           breaks=yticks, labels=NULL)) +
    NoLegend() + 
    coord_flip() + 
    scale_y_reverse()
  
  
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
