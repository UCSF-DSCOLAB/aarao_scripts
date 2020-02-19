library(assertthat) # Assertions
library(Seurat)  # Seurat
library(ggplot2)  # Pretty plots
library(cowplot)  # Grid plot
library(grid)  #  for plotting multiple plots in one frame
library(gridExtra)  #  for plotting multiple plots in one frame
library(scales)  # to access break formatting functions
library(fitdistrplus)  # To fit the negbinom distro
library(RColorBrewer)  # For colors

identify_hto_clusters <- function(sobj, sample_names, quantiles=NULL, outlier_quantile=0.995, 
                                  background_quantile=0.95, visual_outlier_frac=0.001, assay_name="Hashes",
                                  sample_colors=NULL){
  ## A function to identify HTO clusters based on a negative binomial 
  ## background distribution of HTO counts
  ##
  ## Parameters:
  ##     sobj : The Seurat object (REQUIRED)
  ##     sample_names : A named chr or list with names as the HTO and values 
  ##         as the sample names (REQUIRED)
  ##     quantiles : A list of quantiles to use for identifying the 
  ##         outliers and background. A copy of the automatically generated
  ##         one is returned when run for the first time. This can be tweaked
  ##         and used to rerun with tweaked parameters. (DEFAULT=NULL)
  ##     outlier_quantile : The HTO count quantile ABOVE which we will 
  ##         identify a cell as being an "outlier" for that HTO 
  ##         (DEFAULT=0.995)
  ##     background_quantile : The non-outlier HTO count quantile BELOW which
  ##         we will define a cell as having "background" expression of the 
  ##         HTO (DEFAULT=0.95)
  ##     visual_outlier_frac : For plotting raw and scaled plots, what is the
  ##         fraction of highed HTO counts that must be dropped for visual
  ##         clarity (DEFAULT=0.001)
  ##     assay_name : The assay within `sobj` that contains the HTO sounts
  ##         (DEFAULT="Hashes")
  ##     sample_colors : a named chr or list of colors corresponding to each
  ##         sample name (DEFAULT=NULL)
  ## 
  ## Return value
  ##     A list with the following elements
  ##         background_plots : A list of background plots showing the 
  ##             distribution of each HTO counts in the background set and 
  ##             the fitted Negative Binomial distribution.
  ##         quantiles : A list of quantile values used to identify the 
  ##             outlier and background values of reach HTO
  ##         raw_plots : A list of plots or raw counts in normal 
  ##             (below diagonal) and log10 (above diagonal) scale for each 
  ##             HTO pair.
  ##         scaled_plots : A list of plots of scaled HTO counts in a similar
  ##             form as `raw_plots` but only below the diagonal.
  ##         clusters : a vector of cluster identities per cell in the same
  ##             order as the hashes in the input assay.
  ##         colors : named chr of colors used in the plots. Can be used to
  ##             tweak colors according to your needs.
  
  results <- list(
    background_plots = NULL,
    quantiles = NULL,
    raw_plots = NULL,
    scaled_plots = NULL,
    clusters = NULL,
    colors = NULL
  )
  assay_obj <- GetAssay(sobj, assay=assay_name)
  suppmsg <- assert_that(!is.null(assay_obj), msg=paste0('The Seurat object must have an assay named `', assay_name, '`'))
  
  scaled_hashes <- as.data.frame(t(as.matrix(assay_obj@scale.data)))
  raw_hashes <- as.data.frame(t(as.matrix(assay_obj@counts)))

  hto_names <- names(sample_names)
  l <- length(hto_names)
  if (sum(hto_names %in% colnames(scaled_hashes)) != length(hto_names)){
    warning("Not all HTO names provided in `sample_names` are present in the assay object. Using the intersection.")
    hto_names <- hto_names[hto_names %in% colnames(scaled_hashes)]
    assert_that(length(hto_names) >= 2, msg='Require at least 2 HTOs to continue')
    l <- length(hto_names)
    sample_names <- sample_names[hto_names]
  }
  
  scaled_hashes <- scaled_hashes[, hto_names, drop=F]
  colnames(scaled_hashes) <- sample_names
  raw_hashes <- raw_hashes[, hto_names, drop=F]
  colnames(raw_hashes) <- sample_names
  
  if (is.null(sample_colors)){
    sample_colors <- c(brewer.pal(9, 'Set1'), brewer.pal(12, 'Set3'))[c(1:(3+length(sample_names)))]
    names(sample_colors) <- c('NEGATIVE', 'MULTIPLET', 'OTHER', as.vector(sample_names))
  } else{
    assert_that(all(c('NEGATIVE', 'MULTIPLET', 'OTHER', as.vector(sample_names)) %in% names(sample_colors)), 
                msg=paste(c('`sample_colors` must have a colors for each of', 'NEGATIVE', 'MULTIPLET', 
                            'OTHER', sample_names), collapse=" "))
  }
  
  
  if (!(is.null(outlier_quantile)&is.null(quantiles))){
    background_plots <- list()
    if (!is.null(quantiles)){
      expected_names <- c(paste(sample_names, 'outlier', sep='.')) 
      if (!all(expected_names %in% names(unlist(quantiles)))){
        assert_that(
          !is.null(outlier_quantile), 
          msg = paste0('If a Fixed value for all quantile outliers is not provided, `outlier_quantile` must be non-NULL'))
        for (i in sample_names){
          quantiles[[i]][['outlier']] <- quantile(raw_hashes[[i]], outlier_quantile)[[1]]
        }
      }
      expected_names <- c(paste(sample_names, 'background', sep='.'))
      if (!all(expected_names %in% names(unlist(quantiles)))){
        assert_that(
          !is.null(background_quantile), 
          msg = paste0('If a Fixed value for all quantile backgrounds is not provided, `background_quantile` must be non-NULL'))
      }
    } else { # outlier_quantile is not NULL
      assert_that(
        !is.null(background_quantile), 
        msg = paste0('`background_quantile` cannot be NULL if quantiles is not provided but ',
                     '`outlier_quantile_threshold` is'))
      quantiles <- lapply(raw_hashes, function(x) {
        list(outlier=quantile(x, outlier_quantile)[[1]])
        })
    }
    
    for (i in sample_names){
      # Drop outliers with super high values
      hashes <- raw_hashes[raw_hashes[i] < quantiles[[i]][['outlier']],,drop=F]
      
      
      # Get the cells that are NOT in the top quantile (i.e. get the background)
      #quantiles[[i]][['confidence']] <- quantile(hashes[[i]], confidence_quantile)
      #hashes <- raw_hashes[raw_hashes[i] < quantiles[[i]][['confidence']],]
      
      # Discard top 0.5% as outliers
      #hashes <- hashes[hashes[i] < quantile(hashes[[i]], 0.995),]

      # Get the fit
      fit <- fitdist(hashes[[i]], "nbinom")
      # Get a line corresponding to the fit.
      fitD <- dnbinom(0:max(hashes[[i]]), size=fit$estimate['size'], mu=fit$estimate['mu'])
      # Get the background threshold quantile value
      if (is.null(quantiles[[i]][['background']])){
        quantiles[[i]][['background']] <- qnbinom(background_quantile, 
                                                  size=fit$estimate['size'], 
                                                  mu=fit$estimate['mu'])
      }
      
      # Create a data frame to make ggplot happy
      fitD <- data.frame(fit=fitD, idx=c(1:length(fitD)))
      fitD[['col']] = fitD$idx >= quantiles[[i]][['background']]
    
      background_plots[[i]] <- ggplot(hashes) + 
                    geom_histogram(aes_string(x=paste0("`", i, "`"), y = "..density.."), 
                                   binwidth=density(hashes[[i]])$bw) + 
                    geom_area(data=fitD, aes(x=idx, y=fit, col=col, fill=col), size=1.5) +
                    scale_fill_manual("Prediction", values=alpha(c('red', 'blue'), 0.25)) +
                    scale_color_manual("Prediction", values=alpha(c('red', 'blue'), 0.25)) +
                    geom_vline(xintercept=quantiles[[i]][['background']], lty=2) + 
                    annotate("text", x=quantiles[[i]][['background']], y=max(fitD$fit)/2, 
                             label = paste0('Threshold=', quantiles[[i]][['background']]), 
                             angle = 90, vjust=-1.5)
    }
  }
  
  clusters <- data.frame(row.names=rownames(raw_hashes))
  for (i in sample_names){
    clusters[i] <- raw_hashes[i] >= quantiles[[i]][['background']]
    clusters[clusters[[i]], i] <- i
  }
  
  clusters$cluster <- apply(clusters, 1, function(x) {
    ifelse(sum(x!=FALSE)==0, 'NEGATIVE', ifelse(sum(x!=FALSE)==1, x[x!=FALSE], 'MULTIPLET'))
    })
  
  raw_limits <- lapply(raw_hashes, function(x){
    c(NA, quantile(x, 1-visual_outlier_frac)[[1]])
    })
  scaled_limits <- lapply(scaled_hashes, function(x){
    c(NA, quantile(x, 1-visual_outlier_frac)[[1]])
    })
  
  raw_hashes$cluster <- factor(clusters$cluster, levels=c(sample_names, 'NEGATIVE', 'MULTIPLET'))
  scaled_hashes$cluster <- factor(clusters$cluster, levels=c(sample_names, 'NEGATIVE', 'MULTIPLET'))
  
  
  raw_plots <- rep(list(list()), l)
  scaled_plots <- rep(list(list()), l)
  for (i in c(1:l)){
    for (j in c(i:l)){
      if (i == j){
        raw_plots[[i]][[i]] = ggplot(data.frame()) + geom_blank()
        scaled_plots[[i]][[i]] = ggplot(data.frame()) + geom_blank()
        next
      }
      
      temp <- as.data.frame(cbind(raw_hashes[, c(sample_names[i], sample_names[j]), drop=F], 
                                  scaled_hashes[, c(sample_names[i], sample_names[j], 'cluster'), drop=F]))
      colnames(temp) <- c(paste0('RAW_', c(sample_names[i], sample_names[j])), paste0('SCALED_', c(sample_names[i], sample_names[j])), 'cluster')
      temp$cluster <- as.vector(temp$cluster)
      temp[!temp$cluster%in%c('NEGATIVE', 'MULTIPLET', sample_names[i], sample_names[j]), 'cluster'] <- 'OTHER'
      temp$cluster <- factor(temp$cluster, 
                             levels=c('NEGATIVE', 'MULTIPLET', sample_names[i], sample_names[j], 'OTHER'))
      
      temp_table <- table(temp$cluster)
      temp_table <- round(temp_table*100/sum(temp_table), 2)
      temp_table <- sapply(names(temp_table), function(x){paste0(x, ": ", temp_table[x], "%")})
      
      s1 <- sample_names[[i]]
      s2 <- sample_names[[j]]

      raw_plots[[j]][[i]] <- ggplot(temp, aes_string(x=paste0("`RAW_", sample_names[i], "`"), 
                                                     y=paste0("`RAW_", sample_names[j], "`"), 
                                                     col="cluster", fill="cluster")) + 
                             geom_point(size=2, shape=23) +
                             scale_x_continuous(limits=c(raw_limits[[s1]][1], raw_limits[[s1]][2])) +
                             scale_y_continuous(limits=c(raw_limits[[s2]][1], raw_limits[[s2]][2])) +
                             scale_color_manual(values=sample_colors, labels=temp_table, drop=F) +
                             scale_fill_manual(values=sample_colors, labels=temp_table, drop=F) +
                             theme_minimal()
      raw_plots[[i]][[j]] <- ggplot(temp, aes_string(x=paste0("`RAW_", sample_names[i], "`"), 
                                                     y=paste0("`RAW_", sample_names[j], "`"), 
                                                     col="cluster", fill="cluster")) + 
                             geom_point(size=2, shape=23) +
                             scale_x_log10(limits=c(raw_limits[[s1]][1], raw_limits[[s1]][2])) +                     
                             scale_y_log10(limits=c(raw_limits[[s2]][1], raw_limits[[s2]][2])) +
                             scale_color_manual(values=sample_colors, labels=temp_table, drop=F) +
                             scale_fill_manual(values=sample_colors, labels=temp_table, drop=F) +
                             theme_minimal()
      
      scaled_plots[[j]][[i]] <- ggplot(temp, aes_string(x=paste0("`SCALED_", sample_names[i], "`"), 
                                                        y=paste0("`SCALED_", sample_names[j], "`"), 
                                                        col="cluster", fill="cluster")) + 
                             geom_point(size=2, shape=23) +
                             scale_x_continuous(limits=c(scaled_limits[[s1]][1], scaled_limits[[s1]][2])) +
                             scale_y_continuous(limits=c(scaled_limits[[s2]][1], scaled_limits[[s2]][2])) +
                             scale_color_manual(values=sample_colors, labels=temp_table, drop=F) +
                             scale_fill_manual(values=sample_colors, labels=temp_table, drop=F) +
                             theme_minimal()
      scaled_plots[[i]][[j]] = ggplot(data.frame()) + geom_blank()
    }
  }
  print('You can plot a grid with ')
  print(paste0("plot_grid(plotlist=unlist(results[['raw_plots']], recursive = FALSE), ncol=", l, ")"))
  print(paste0("plot_grid(plotlist=unlist(results[['scaled_plots']], recursive = FALSE), ncol=", l, ")"))
  print(paste0("plot_grid(plotlist=results[['background_plots']], ncol=2)"))
  results[['background_plots']] <- background_plots
  results[['raw_plots']] <- raw_plots
  results[['scaled_plots']] <- scaled_plots
  results[['quantiles']] <- quantiles
  results[['clusters']] <- clusters$cluster
  results[['colors']] <- sample_colors
  return(results)
}