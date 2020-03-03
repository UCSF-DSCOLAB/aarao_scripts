library(assertthat) # Assertions
library(Seurat)  # Seurat
library(ggplot2)  # Pretty plots
library(cowplot)  # Grid plot
library(grid)  #  for plotting multiple plots in one frame
library(gridExtra)  #  for plotting multiple plots in one frame
library(scales)  # to access break formatting functions
library(RColorBrewer)  # For colors

demux_by_inflexions <- function(sobj, sample_names, inflexions=NULL, visual_outlier_frac=0.001, 
                                assay_name="Hashes", slot='scale.data', dens_windows=1000, 
                                num_inflexions=5, sample_colors=NULL){
  ## A function to identify HTO clusters based the inflexion point in the density plot of the
  ## HTO count distribution
  ##
  ## Parameters:
  ##     sobj : The Seurat object (REQUIRED)
  ##     sample_names : A named chr or list with names as the HTO and values 
  ##         as the sample names (REQUIRED)
  ##     inflexions : A named chr or list with names as the HTO and values 
  ##         are the inflexion points to use. (DEFAULT=NULL)
  ##     visual_outlier_frac : For plotting raw and scaled plots, what is the
  ##         fraction of highed HTO counts that must be dropped for visual
  ##         clarity (DEFAULT=0.001)
  ##     assay_name : The assay within `sobj` that contains the HTO sounts
  ##         (DEFAULT="Hashes")
  ##     slot : The slow used to find inflexions. One of counts, data, or scale.data
  ##         (DEFAULT=scale.data)
  ##     dens_windows : The number of windows used while calculating density. Higher values 
  ##         give finer detail. (DEFAULT=1000)
  ##     num_inflexions : The number of inflexions we will report at max (DEFAULT=5)
  ##     sample_colors : a named chr or list of colors corresponding to each
  ##         sample name (DEFAULT=NULL)
  ## 
  ## Return value
  ##     A list with the following elements
  ##         background_plots : A list of background plots showing the 
  ##             distribution of each HTO counts in the background set and 
  ##             the fitted Negative Binomial distribution.
  ##         inflexions : A list with the inflexion used for each HTO
  ##         secondary_inflexions : A list of addditional inflexions for each HTO
  ##         raw_plots : A list of plots or raw counts in normal 
  ##             (below diagonal) and log10 (above diagonal) scale for each 
  ##             HTO pair.
  ##         scaled_plots : A list of plots of scaled HTO counts in a similar
  ##             form as `raw_plots` but only below the diagonal.
  ##         cell_ids : a vector of identities, one for each cell in the input assay,
  ##             in the same cell order as the assay provided.
  ##         summary : A summary of the number and percent of cells identified as either an
  ##             expected sample, or as a MULTIPLET, or as NEGATIVE for all HTOs.
  ##         colors : named chr of colors used in the plots. Can be used to
  ##             tweak colors according to your needs.
  
  assert_that(num_inflexions > 0, msg="num_inflexions must be > 0")
  results <- list(
    background_plots = NULL,
    inflexions = NULL,
    secondary_inflexions = NULL,
    raw_plots = NULL,
    scaled_plots = NULL,
    cell_ids = NULL,
    summary = NULL,
    colors = NULL
  )
  assay_obj <- GetAssay(sobj, assay=assay_name)
  suppmsg <- assert_that(!is.null(assay_obj), msg=paste0('The Seurat object must have an assay named `', assay_name, '`'))
  
  scaled_hashes <- as.data.frame(t(as.matrix(GetAssayData(assay_obj, slot = 'scale.data'))))
  raw_hashes <- as.data.frame(t(as.matrix(GetAssayData(assay_obj, slot = 'counts'))))

  hto_names <- names(sample_names)
  l <- length(hto_names)
  suppmsg <- assert_that(sum(sort(hto_names)!=sort(colnames(raw_hashes)))==0, msg='Not all hashes in sobj are accounted for in sample_names')
  
  scaled_hashes <- scaled_hashes[, hto_names, drop=F]
  raw_hashes <- raw_hashes[, hto_names, drop=F]

  raw_limits <- lapply(raw_hashes, function(x){
    c(NA, quantile(x, 1-visual_outlier_frac)[[1]])
  })
  scaled_limits <- lapply(scaled_hashes, function(x){
    c(NA, quantile(x, 1-visual_outlier_frac)[[1]])
  })
  
  
  if (is.null(sample_colors)){
    sample_colors <- c(brewer.pal(9, 'Set1'), brewer.pal(12, 'Set3'))[c(1:(3+length(sample_names)))]
    names(sample_colors) <- c('NEGATIVE', 'MULTIPLET', 'OTHER', as.vector(sample_names))
  } else{
    suppmsg <- assert_that(all(c('NEGATIVE', 'MULTIPLET', 'OTHER', as.vector(sample_names)) %in% names(sample_colors)), 
                           msg=paste(c('`sample_colors` must have a colors for each of', 'NEGATIVE', 'MULTIPLET', 
                                    'OTHER', sample_names), collapse=" "))
  }
  
  if (is.null(inflexions)){
    results[['inflexions']] <- list()
  } else {
    # This is a valid check since we've previously ensured sample_names is valid.
    suppmsg <- assert_that(sum(sort(names(inflexions)) != sort(names(sample_names)))==0, msg='Inflexion names does not match the HTO names in sobj')
    results[['inflexions']] <- inflexions
  }
  
  hashes <- as.data.frame(t(as.matrix(GetAssayData(assay_obj, slot = slot))))[, hto_names, drop=F]
  background_limits <- lapply(hashes, function(x){
    c(NA, quantile(x, 1-visual_outlier_frac)[[1]])
  })
  
  all_inflexions <- list()
  results[['background_plots']] <- list()
  for (sample in names(sample_names)) {
    if (is.null(results[['inflexions']][[`sample`]])) {
      dens <- density(hashes[[`sample`]], n = dens_windows)
      maxima <- which.max(dens$y)
      #all_inflexions[[sample]] = data.frame(row.names = c('x', 'y'))
      all_inflexions[[sample]] = c()
      prev = NULL
      #pi_id=1
      for (i in c(maxima:length(dens$x))){
        if (is.null(prev)){
          prev = sign((dens$y[i]-dens$y[i-1])/(dens$x[i]-dens$x[i-1]))
          next
        }
        if (length(all_inflexions[[sample]]) == num_inflexions || i==length(dens$x)) {
          suppmsg <- assert_that(!is.null(all_inflexions[[sample]][1]), 
                                 msg=paste0('Could not identify an inflexion for ', sample, 
                                            '. Consider increasing inflexion_searchspace.'))
          results[['inflexions']][[sample]] <- all_inflexions[[sample]][1]
          results[['secondary_inflexions']][[sample]] <- all_inflexions[[sample]][c(2:num_inflexions)][!is.na(all_inflexions[[sample]][c(2:num_inflexions)])]
          break
        }
        s = sign((dens$y[i]-dens$y[i-1])/(dens$x[i]-dens$x[i-1]))
        if (s != prev){
          # This is an inflexion. i-1 was the point of inflexion
          if (s == 1){
            # This is a positive inflexion and that's what we need
            all_inflexions[[sample]] = c(all_inflexions[[sample]], dens$x[i-1])
            #all_inflexions[[sample]][[pi_id]] = c(dens$x[i-1], dens$y[i-1])
            #pi_id <- pi_id + 1
          }
          prev = s
        }
      }
    }
    results[['background_plots']][[sample]] <- ggplot(hashes, aes_string(x=paste0("`", sample, "`"))) + 
                                               geom_density() +
                                               geom_vline(xintercept=results[['inflexions']][[sample]], lty=2) +
                                               xlim(background_limits[[`sample`]])
  }
  results[['inflexions']] <- unlist(results[['inflexions']])
  
  clusters <- data.frame(row.names=rownames(hashes))
  for (sample in names(sample_names)) {
    clusters[sample] <- hashes[`sample`] >= results[['inflexions']][[`sample`]]
    clusters[clusters[[`sample`]], `sample`] <- sample_names[`sample`]
  }
  
  clusters$sample_name <- apply(clusters, 1, function(x) {
    ifelse(sum(x!=FALSE)==0, 'NEGATIVE', ifelse(sum(x!=FALSE)==1, x[x!=FALSE], 'MULTIPLET'))
    })
  
  raw_hashes$sample_name <- factor(clusters$sample_name, levels=c(sample_names, 'NEGATIVE', 'MULTIPLET'))
  scaled_hashes$sample_name <- factor(clusters$sample_name, levels=c(sample_names, 'NEGATIVE', 'MULTIPLET'))
  
  raw_plots <- rep(list(list()), l)
  scaled_plots <- rep(list(list()), l)
  for (i in c(1:l)){
    for (j in c(i:l)){
      if (i == j){
        raw_plots[[i]][[i]] = ggplot(data.frame()) + geom_blank()
        scaled_plots[[i]][[i]] = ggplot(data.frame()) + geom_blank()
        next
      }
      
      temp <- as.data.frame(cbind(raw_hashes[, c(hto_names[i], hto_names[j]), drop=F], 
                                  scaled_hashes[, c(hto_names[i], hto_names[j], 'sample_name'), drop=F]))
      colnames(temp) <- c(paste0('RAW_', c(hto_names[i], hto_names[j])), 
                          paste0('SCALED_', c(hto_names[i], hto_names[j])), 'sample_name')
      temp$sample_name <- as.vector(temp$sample_name)
      temp[!temp$sample_name%in%c('NEGATIVE', 'MULTIPLET', sample_names[i], sample_names[j]), 'sample_name'] <- 'OTHER'
      temp$sample_name <- factor(temp$sample_name, 
                             levels=c('NEGATIVE', 'MULTIPLET', sample_names[i], sample_names[j], 'OTHER'))
      
      temp_table <- table(temp$sample_name)
      temp_table <- round(temp_table*100/sum(temp_table), 2)
      temp_table <- sapply(names(temp_table), function(x){paste0(x, ": ", temp_table[x], "%")})
      
      s1 <- hto_names[[i]]
      s2 <- hto_names[[j]]

      raw_plots[[j]][[i]] <- ggplot(temp, aes_string(x=paste0("`RAW_", s1, "`"), 
                                                     y=paste0("`RAW_", s2, "`"), 
                                                     col="sample_name", fill="sample_name")) + 
                             geom_point(size=2, shape=23) +
                             scale_x_continuous(limits=c(raw_limits[[s1]][1], raw_limits[[s1]][2])) +
                             scale_y_continuous(limits=c(raw_limits[[s2]][1], raw_limits[[s2]][2])) +
                             scale_color_manual(values=sample_colors, labels=temp_table, drop=F) +
                             scale_fill_manual(values=sample_colors, labels=temp_table, drop=F) +
                             theme_minimal()
      raw_plots[[i]][[j]] <- ggplot(temp, aes_string(x=paste0("`RAW_", s1, "`"), 
                                                     y=paste0("`RAW_", s2, "`"), 
                                                     col="sample_name", fill="sample_name")) + 
                             geom_point(size=2, shape=23) +
                             scale_x_log10(limits=c(raw_limits[[s1]][1], raw_limits[[s1]][2])) +                     
                             scale_y_log10(limits=c(raw_limits[[s2]][1], raw_limits[[s2]][2])) +
                             scale_color_manual(values=sample_colors, labels=temp_table, drop=F) +
                             scale_fill_manual(values=sample_colors, labels=temp_table, drop=F) +
                             theme_minimal()
      
      scaled_plots[[j]][[i]] <- ggplot(temp, aes_string(x=paste0("`SCALED_", s1, "`"), 
                                                        y=paste0("`SCALED_", s2, "`"), 
                                                        col="sample_name", fill="sample_name")) + 
                             geom_point(size=2, shape=23) +
                             scale_x_continuous(limits=c(scaled_limits[[s1]][1], scaled_limits[[s1]][2])) +
                             scale_y_continuous(limits=c(scaled_limits[[s2]][1], scaled_limits[[s2]][2])) +
                             scale_color_manual(values=sample_colors, labels=temp_table, drop=F) +
                             scale_fill_manual(values=sample_colors, labels=temp_table, drop=F) +
                             theme_minimal()
      scaled_plots[[i]][[j]] = ggplot(data.frame()) + geom_blank()
    }
  }
  print("Demuxed the following identities")
  results[['summary']] <- data.frame(rbind(table(raw_hashes$sample_name), round(100*table(raw_hashes$sample_name)/length(raw_hashes$sample_name), 2)), row.names=c('counts', 'pct'))
  print(results[['summary']])
  print('You can plot a grid with: ')
  print(paste0("plot_grid(plotlist=unlist(results[['raw_plots']], recursive = FALSE), ncol=", l, ")"))
  print(paste0("plot_grid(plotlist=unlist(results[['scaled_plots']], recursive = FALSE), ncol=", l, ")"))
  print(paste0("plot_grid(plotlist=results[['background_plots']], ncol=2)"))
  
  results[['raw_plots']] <- raw_plots
  results[['scaled_plots']] <- scaled_plots
  
  results[['cell_ids']] <- clusters$sample_name
  results[['colors']] <- sample_colors
  
  return(results)
}
