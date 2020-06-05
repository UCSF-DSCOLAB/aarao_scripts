#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(assertthat)
  library(dplyr)
  library(edgeR)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(ggfortify)
  library(gridExtra)
  library(limma)
  library(methods)
  library(reshape2)
  library(RColorBrewer)
})

# poor vs rich
# logFc < 0 -> poor/rich < 0 -> upregulated in poor
# logFc > 0 -> poor/rich > 0 -> downregulated in poor


# Function definitions

args = list(
  COUNTS_TSV=NA,  # Directory with folders containing fastqs
  CLASSES_TSV=NA,  # Samples to pull from the directory
  OUT_FOLDER=".",  # Where to write results
  MODEL_COVARIATE=c(),
  ANNOT_COLS_TSV=NULL,
  STRING_FOLDER=paste(Sys.getenv("KLAB"), "/ipi/data/databases/string", sep=""),
  RSCRIPTS_DIR=paste(Sys.getenv("SCRIPTS"), "R", sep="/"),  # Where to find aux scripts
  RUN_STRING_DE=TRUE,
  RUN_STRING_PCA=TRUE,
  RUN_DAVID_DE=TRUE,
  RUN_DAVID_PCA=TRUE,
  MOST_VARIABLE_CUTOFF=c(200, 500, 1000)
  )

argsClasses = list(
  COUNTS_TSV=as.character,
  CLASSES_TSV=as.character,
  OUT_FOLDER=as.character,
  MODEL_COVARIATE=as.character,
  ANNOT_COLS_TSV=as.character,
  STRING_FOLDER=as.character,
  RSCRIPTS_DIR=as.character,
  RUN_STRING_DE=as.logical,
  RUN_STRING_PCA=as.logical,
  RUN_DAVID_DE=as.logical,
  RUN_DAVID_PCA=as.logical,
  MOST_VARIABLE_CUTOFF=as.numeric
)

user_args <- commandArgs(trailingOnly=TRUE)
parsed_user_args <- c()
for (i in user_args){
  key_val = strsplit(i, split="=", fixed=TRUE)[[1]]
  assert_that(length(key_val) == 2, msg=paste0("Invalid option `", i, "`. ",
                                               "Args must be of the form ",
                                               "KEY=VAL ."))
  if (!key_val[1] %in% names(args)){
    print(paste0("WARNING: Invalid key provided : ", i, ""))
    next
  } else {
    if (key_val[1] %in% parsed_user_args){
      args[[key_val[1]]] = c(args[[key_val[1]]], argsClasses[[key_val[1]]](key_val[2]))
    } else {
      args[[key_val[1]]] = argsClasses[[key_val[1]]](key_val[2])
      parsed_user_args <- c(parsed_user_args, key_val[1])
    }
  }
}

cat("Running with arguments:\n")
cat(paste0("COUNTS_TSV : ", args$COUNTS_TSV, " (", mode(args$COUNTS_TSV), ")\n"))
cat(paste0("CLASSES_TSV : ", args$CLASSES_TSV, " (", mode(args$CLASSES_TSV), ")\n"))
cat(paste0("OUT_FOLDER  : ", args$OUT_FOLDER, " (", mode(args$OUT_FOLDER), ")\n"))
cat(paste("MODEL_COVARIATE : ", args$MODEL_COVARIATE, " (", mode(args$MODEL_COVARIATE), ")\n", sep='', collapse=''))
cat(paste0("ANNOT_COLS_TSV : ", args$ANNOT_COLS_TSV, " (", mode(args$ANNOT_COLS_TSV), ")\n"))
cat(paste0("STRING_FOLDER : ", args$STRING_FOLDER, " (", mode(args$STRING_FOLDER), ")\n"))
cat(paste0("RSCRIPTS_DIR : ", args$RSCRIPTS_DIR, " (", mode(args$RSCRIPTS_DIR), ")\n"))
cat(paste0("RUN_STRING_DE : ", args$RUN_STRING_DE, " (", mode(args$RUN_STRING_DE), ")\n"))
cat(paste0("RUN_STRING_PCA : ", args$RUN_STRING_PCA, " (", mode(args$RUN_STRING_PCA), ")\n"))
cat(paste0("RUN_DAVID_DE : ", args$RUN_DAVID_DE, " (", mode(args$RUN_DAVID_DE), ")\n"))
cat(paste0("RUN_DAVID_PCA : ", args$RUN_DAVID_PCA, " (", mode(args$RUN_DAVID_PCA), ")\n"))
cat(paste("MOST_VARIABLE_CUTOFF : ", args$MOST_VARIABLE_CUTOFF, " (", mode(args$MOST_VARIABLE_CUTOFF), ")\n", sep='', collapse=''))

if (is.na(args$COUNTS_TSV) | is.na(args$CLASSES_TSV)) {
  stop(paste0("Need COUNTS_TSV=<COUNTS.tsv> and CLASSES_TSV=<CLASSES.tsv> to continue. ",
              "Other options include OUT_FOLDER, BATCH_CORRECT, ANNOT_COLS_TSV, ",
              "RSCRIPTS_DIR, and STRING_FOLDER ."),
       call.=FALSE)
}

source(paste0(args$RSCRIPTS_DIR, '/coolmap2_heatmap3.R'))
source(paste0(args$RSCRIPTS_DIR, '/biomart_operations.R'))
source(paste0(args$RSCRIPTS_DIR, '/string_operations.R'))

cat ("Reading counts.....")
counts <- read.table(args$COUNTS_TSV, sep='\t', header=1, row.names=1)
cat ("Done\n")

cat ("Reading classes.....")
groups <- read.table(args$CLASSES_TSV, sep='\t', header=1, row.names=1, stringsAsFactors = T, na.strings = 'unknown')
cat ("Done\n")

if (!is.null(args$ANNOT_COLS_TSV)){
  cat ("Reading Additional annotation columns.....")
  annot_cols <- read.table(args$ANNOT_COLS_TSV, sep='\t', 
                           header=1, row.names=1, 
                           stringsAsFactors = FALSE, 
                           quote="")
  cat ("Done\n")
} else {
  annot_cols <- data.frame()
}


if (is.numeric(groups$plate)){
  cat(paste0('`plate` in <GROUPS.TSV> cannot be numerical. Renaming to ',
             'prefix with `plate` (i.e. `1` becomes `plate1`). This ',
             'is the expected format.\n'))
  groups$plate <- paste0('plate', groups$plate)
}

if (!'group' %in% colnames(groups)) {
    stop("<GROUPS.TSV> must have a column named group", call.=FALSE)
} else if (is.numeric(groups$group)){
  cat(paste0('`group` in <GROUPS.TSV> cannot be numerical. Renaming to ',
             'prefix with `g` (i.e. `1` becomes `g1`).\n'))
  groups$group <- as.factor(paste0('g', groups$group))
}

if (args$RUN_STRING_DE|args$RUN_STRING_PCA){
  cat ("Populating String DB.....")
  string_db <- setup_string(args$STRING_FOLDER)
  cat ("Done\n")
}

if (length(rownames(groups)) != length(colnames(counts)) ||
        !all(rownames(groups) %in% colnames(counts)) ||
        !all(colnames(counts) %in% rownames(groups))) {
    warning("<COUNTS.TSV> and <GROUPS.TSV> do not have the same columns. Continuing with the intersection", call.=FALSE)
    common_cols <- intersect(colnames(counts), rownames(groups))
} else {
    common_cols <- rownames(groups)
}

groups <- groups[common_cols,, drop=F]
groups <- groups[order(groups$group, rownames(groups)),, drop=F]
counts <- counts[, rownames(groups)]


if (!all(rownames(annot_cols) %in% colnames(groups))) {
  nc <- rownames(annot_cols)[!rownames(annot_cols) %in% colnames(groups)]
  warning(paste0("The following fields in <ANNOT_COLS.TSV> are not in ",
                 "<GROUPS.TSV>: ", 
                 paste(nc, collapse=", "),
                 " ....Discarding missing columns"),
                 call.=FALSE)
  cc <- rownames(annot_cols)[rownames(annot_cols) %in% colnames(groups)]
  annot_cols <- annot_cols[cc, , drop=F]
}

# Add IPI batch annotations into the groups dataframe
fixed_levels <- list()
fixed_levels[['sequencer']] <- c('hiseq', 'novaseq')
groups$sequencer <- factor(x=sapply(groups$plate, function(x){if (as.numeric(gsub('plate', '', x)) > 25) 'novaseq' else 'hiseq'}), levels=fixed_levels[['sequencer']])
fixed_levels[['washes']] <- c('extraWash', 'noExtraWash')
groups$washes <- factor(x=sapply(groups$plate, function(x){if (as.numeric(gsub('plate', '', x)) > 8) 'extraWash' else 'noExtraWash'}), levels=fixed_levels[['washes']])

# Add proper levels to EHK and plates
groups$EHK <- factor(as.vector(groups$EHK),
                     levels = as.vector(c(1:10)))

fixed_levels[['plate']] <- as.vector(levels(groups$plate))
fixed_levels[['plate']] <- fixed_levels[['plate']][order(as.numeric(gsub("plate", "", fixed_levels[['plate']])))]
groups$plate <- factor(as.vector(groups$plate), levels = fixed_levels[['plate']])

cat('\nThe plate distribution across plates is seen to be:\n')
print(table(groups[,c('sequencer', 'plate')]))
cat('\nThe wash distribution across plates is seen to be:\n')
print(table(groups[,c('washes', 'plate')]))

# See https://support.bioconductor.org/p/36029/#100297
# See https://support.bioconductor.org/p/66251/#66252
# This is how Gordon Smyth suggests you handle batch effect in the DE
design_formula = '~0+group'
factor_covariates = c('group')

for (m in  args$MODEL_COVARIATE) {
  if (!m%in%colnames(groups)){
    cat(paste0('Cannot add ', m, ' to the model since it is not in the groups dataframe.\n'))
    next
  }
  
  if (length(unique(groups[[m]])) > 1) {
    if (is.factor(groups[[m]])){
      factor_covariates = c(factor_covariates, m)
    }
    design_formula <- paste0(design_formula, '+', m)
    useSequencerToBC <- TRUE
    cat(paste0('Added ', m, ' to the model.\n'))
  } else {
    cat(paste0('Covariate ', m, ' is constant across all samples. Skipping.\n'))
  }
}

design <- model.matrix(as.formula(design_formula), data=groups) #limma

for (l in levels(groups$group)) {
  colnames(design) <- gsub(paste0("^", "group", l, "$"), l, colnames(design))
}
cat("\nUsing design Model: \n")
print(design)

logcpm_matrices <- list()
dge_list <- DGEList(counts=counts, group=groups$group)
logcpm_matrices[['unnormalized']] <- cpm(dge_list, log=TRUE)
cat ("Made DGElist.\n")
dge_list <- dge_list[ rowSums(dge_list$counts) > 60, keep.lib.size=FALSE ]
dge_list <- calcNormFactors(dge_list, method="TMM")
logcpm_matrices[['TMMnormalized']] <- cpm(dge_list, log=TRUE)
cat ("Normalization done.\n")

v <- voom(dge_list, design)
fit <- lmFit(v, design)
cat ("Limma fit done.\n")

i=1
conts <- c()
cnames <- c()
combs <- combn(levels(groups$group), 2, simplify=F)
for (comb in combs){
    conts[i] = paste(comb, collapse=" - ")
    cnames[i] = paste(comb, collapse="_vs_")
    i <- i+1
}

cont_matrix <- makeContrasts(contrasts=conts, levels = design )#limma
colnames(cont_matrix) <- cnames
cat ("Using contrast matrix:\n")  #limma
print(cont_matrix)

fit2 <- contrasts.fit(fit, cont_matrix)  #limma
fit2 <- eBayes(fit2, robust=TRUE)  #limma
cat ("Fit2 done.\n")


# Setup the colors for the plots
cols <- list(
  group = brewer.pal(max(length(levels(groups$group)), 3), "Accent"),
  washes = c('black', 'lightgrey'),
  sequencer = c('black', 'lightgrey')
)

numeric_cols <- list()

for (annot in rownames(annot_cols)) {
  if (annot %in% names(cols)){
    next
  }
  
  if (annot_cols[annot, 'type'] == "factor"){
    num_vals <- length(levels(groups[[annot]]))
  } else if (annot_cols[annot, 'type'] == "numeric") {
    if (!annot_cols[annot, 'cmap'] %in% rownames(brewer.pal.info)) {
      stop(paste0('Currently only color palettes from brewer are allowed for numerical annotations',
                  'Choose from one of ', paste(rownames(brewer.pal.info), collapse=", ")))
    }
    numeric_cols[[annot]] <- list(palette=annot_cols[annot, 'cmap'])
    if (endsWith(annot, "_frac")){
      numeric_cols[[annot]][["low"]] <- 0
      numeric_cols[[annot]][["high"]] <- 1
    } else if (endsWith(annot, "_pct")){
      numeric_cols[[annot]][["low"]] <- 0
      numeric_cols[[annot]][["high"]] <- 100
    } else {
      numeric_cols[[annot]][["low"]] <- min(groups[[annot]], na.rm = T)
      numeric_cols[[annot]][["high"]] <- max(groups[[annot]], na.rm = T)
    }
    num_vals <- 10
  } else{
    assert_that(FALSE, msg="Cannot handle non-factor/numeric types for annot right now.")
  }
  
  if(annot_cols[annot, 'cmap'] %in% rownames(brewer.pal.info)) {
    cols[[annot]] = brewer.pal(max(num_vals, 3), annot_cols[annot, 'cmap'])  
  } else {
    # This assumes the user has passed in a function that can be evaluated to
    # generate a palette
    pal = eval(parse(text=annot_cols[annot, 'cmap']))
    cols[[annot]] = pal(num_vals)
  }
} # for (annot in rownames(annot_cols))

labs <- c()
labs.colnames <- c()

legend.text <- c()
legend.fill <- c()

for (annot_group in names(cols)) {
  if (annot_group %in% names(numeric_cols)) {
    pct <- (numeric_cols[[annot_group]][['high']]-numeric_cols[[annot_group]][['low']])/10
    vals <- sapply(groups[[annot_group]], function(x){min((x-numeric_cols[[annot_group]][['low']])%/%pct+1, 10)})
    vals <- cols[[annot_group]][vals]
    vals[is.na(vals)] <- 'black'
    labs <- cbind(labs, vals)
    labs.colnames <- c(labs.colnames, annot_group)
    
    legend.text <- c(legend.text, c("NA", "low", "medium", "high"), "")
    legend.fill <- c(legend.fill, 'grey', cols[[annot_group]][c(1,5, 10)], 'white')
    
  } else {
    labs <- cbind(labs, cols[[annot_group]][groups[[annot_group]]])
    labs.colnames <- c(labs.colnames, annot_group)
  
    legend.text <- c(legend.text, levels(groups[[annot_group]]), "")
    
    if (annot_group %in% c('plate', 'washes', 'sequencer')) {
      # We specifically use the levels since we've manually set them in the order we care
      legend.fill <- c(legend.fill, cols[[annot_group]][factor(levels(groups[[annot_group]]), levels=fixed_levels[[annot_group]])], 'white')
    } else if (annot_group == 'EHK') {
      # Has to be numeric for this to be correct)
      legend.fill <- c(legend.fill, cols[[annot_group]][as.numeric(levels(groups[[annot_group]]))], 'white')
    } else {
      legend.fill <- c(legend.fill, cols[[annot_group]][as.factor(levels(groups[[annot_group]]))], 'white')  
    }
  } # if (annot_group %in% names(numeric_cols)) / else
} # for (annot_group in names(cols)) 

colnames(labs) <- labs.colnames
rownames(labs) <- rownames(groups)

hclust.ward = function(d) hclust(d,method="ward.D2")

cat ("Setup all color palettes.\n")

for (nm in names(logcpm_matrices)) {
  nmlcpm <- logcpm_matrices[[nm]]
  
  cat(paste0('Plotting raw pca with ', nm, ' cpms\n'))
  pcs <- prcomp(t(nmlcpm))
  ap <- list()
  for (annot_name in names(cols)) {
    if (annot_name == 'group'){
      ap[[annot_name]] <- autoplot(pcs, data=groups, shape=annot_name, col=annot_name, size=5, alpha=0.75) + 
        geom_text_repel(aes(label=rownames(groups))) + 
        scale_color_manual(values=cols[[annot_name]]) +
        theme_bw()
    } else if (annot_name %in% names(numeric_cols)){
      ap[[annot_name]] <- autoplot(pcs, data=groups, shape='group', col=annot_name, size=5, alpha=0.75) + 
        scale_color_distiller(palette=numeric_cols[[annot_name]][["palette"]], limits = c(numeric_cols[[annot_name]][["low"]], numeric_cols[[annot_name]][["high"]])) +
        theme_bw()
    } else {
      ap[[annot_name]] <- autoplot(pcs, data=groups, shape='group', col=annot_name, size=5, alpha=0.75) + 
        scale_color_manual(values=cols[[annot_name]]) +
        theme_bw()
    }
  }
  
  png(filename=paste(args$OUT_FOLDER, paste('voom_', gsub('.{6}$', '', nm), '_pca.png',sep=''), sep='/'), width = 30, height =  ceiling(length(ap)/3) *10 , units = "in", res = 150)
  grid.arrange(grobs=ap, ncol = 3)
  dev.off()
  
  x <- as.data.frame(pcs$x)
  x$group = groups$group
  y <- melt(x, id.vars = 'group', variable.name = 'PC')
  
  png(filename=paste(args$OUT_FOLDER, paste('voom_', gsub('.{6}$', '', nm), '_pca_boxplots.png',sep=''), sep='/'), width = 30, height =  ceiling(length(ap)/3) *10 , units = "in", res = 150)
  print(ggplot(y, aes(x=group, y=value)) + geom_boxplot() + stat_compare_means(method = "t.test") + facet_wrap(~PC, ncol=5))
  dev.off()
  
  write_text <- c()
  for (pc in c(1, 2)) {
    best_correlation <- sort(apply(nmlcpm, 1, function(x) {cor(x, pcs$x[,pc], method="spearman")}))
    write_text <- c(write_text, paste0("Top 20 negatively correlated genes with PC", pc))
    write_text <- c(write_text, "Gene\tr\tr^2")
    for (i in names(best_correlation[c(1:20)])){
      write_text <- c(write_text, paste(i, round(best_correlation[i], 2), round(best_correlation[i]**2, 2), sep='\t'))
    }
    write_text <- c(write_text, '')
    write_text <- c(write_text, paste0("Top 20 positively correlated genes with PC", pc))
    write_text <- c(write_text, "Gene\tr\tr^2")
    for (i in rev(names(best_correlation[c((length(best_correlation)-(20-1)):length(best_correlation))]))){
      write_text <- c(write_text, paste(i, round(best_correlation[i], 2), round(best_correlation[i]**2, 2), sep='\t'))
    }
    write_text <- c(write_text, '')
  } # for (pc in c(1, 2))
  fileConn<-file(paste(args$OUT_FOLDER, paste('voom_', gsub('.{6}$', '', nm), '_pca_correlating_genes.txt',sep=''), sep='/'))
  writeLines(write_text, fileConn)
  close(fileConn)
  
  for (pc in c(1, 2)) {
    x <- data.frame(weight=pcs$rotation[, pc][order(-abs(pcs$rotation[, pc]))][c(1:200)])
    x$ensembl <- rownames(x)
    x <- ensg_to_hugo(x, 
                      ensembl_col = 'ensembl',
                      hugo_col = 'hugo',
                      fill = FALSE,
                      handle_duplicates = 'first')
    write.table(x, 
                file=paste(args$OUT_FOLDER, paste0('voom_', gsub('.{6}$', '', nm), '_pca_top200_weights.tsv'), sep='/'),
                row.names = T, 
                col.names=T, 
                sep="\t", 
                quote=F)
    if (args$RUN_STRING_PCA){
      cat(paste0('Running STRING on pc ', pc, '\n'))
      run_string(paste(args$OUT_FOLDER, paste0('voom_', gsub('.{6}$', '', nm), '_PC', pc), sep="/"),
                        top_table=x, 
                        string_db=string_db, 
                        usesig=FALSE,
                        logFcColStr="weight",
                        trend_col="weightage",
                        trend_classes=c('POS', 'NEG'))
    }
  } 
  
  for (mvc in args$MOST_VARIABLE_CUTOFF) {
    cat(paste0('Plotting unbiased heatmap with ', nm, ' cpms across ', mvc, ' variable genes.\n'))
    top_most_variable_genes = rownames(fit2$coefficients)[order(-fit2$sigma)][c(0:mvc)]
    png(filename=paste(args$OUT_FOLDER, paste('voom_', 
                                              gsub('.{6}$', '', nm), 
                                              '_unbiased_top', 
                                              mvc, 
                                              '_heatmap.png',
                                              sep=''), 
                       sep='/'), 
        width = 20, 
        height = 15, units = "in", res = 150)
    coolmap.2(nmlcpm[top_most_variable_genes,], 
              ColSideColors=labs, 
              margins=c(12, 25),
              cbar_extreme=2,
              ColSideColorSeparator=TRUE,
              cexCol=1)
    legend("topright",
           legend=legend.text,
           fill=legend.fill,
           border=FALSE,
           bty="n",
           y.intersp = 0.7,
           cex=0.7,
           ncol=1)
    dev.off()
  } # for (mvc in args$MOST_VARIABLE_CUTOFF)
} # for (nm in names(logcpm_matrices))

for (cname in cnames) {
  top_table = topTable(fit2, coef=cname, sort="p", number=Inf)
  top_table$col <- 0
  top_table$gene <- rownames(top_table)
  top_table <- top_table[order(top_table$P.Value),]
  top_100 <- c(rownames(head(top_table[top_table$logFC<0 & top_table$P.Value <= 0.005,], 100)), rownames(head(top_table[top_table$logFC>0 & top_table$P.Value <= 0.005,], 100)))
  #top_n <- c(rownames(head(top_table[top_table$logFC<0 & top_table$P.Value <= 0.005,], 25)), rownames(head(top_table[top_table$logFC>0 & top_table$P.Value <= 0.005,], 25)))
  #top_n <- c(rownames(head(top_table[top_table$logFC<0 & top_table$P.Value <= 0.005,], 50)), rownames(head(top_table[top_table$logFC>0 & top_table$P.Value <= 0.005,], 50)))
  top_n <-rownames(top_table[top_table$P.Value <= 0.005,])
  top_table[rownames(top_table[top_table$P.Value<0.05,]), "col"] = 1
  top_table[top_100, "col"] = 2
  top_table$col <- as.factor(top_table$col)
  
  write.table(top_table, file=paste(args$OUT_FOLDER, paste('voom_', cname, '.tsv',sep=''), sep='/'), sep='\t', row.names=T, col.names=T, quote=F)
  p <- ggplot(top_table, aes(logFC, -log10(P.Value))) + 
          geom_point(aes(col=col)) + 
          scale_color_manual(values=c("black", "red", "red")) + 
          theme(legend.position = "none") + 
          geom_text_repel(data=filter(top_table, col=='2'), aes(label=gene), size=2, color='black')
  png(filename=paste(args$OUT_FOLDER, paste('voom_', cname, '.png',sep=''), sep='/'), width = 10, height = 10, units = "in", res = 150)
  print(p)
  dev.off()
  
  if (args$RUN_STRING_DE){
    run_string(paste(args$OUT_FOLDER, 
                     paste0(cname, '_DE'), 
                     sep="/"),
               top_table, 
               string_db, 
               usesig=FALSE,
               logFcColStr="logFC",
               trend_col="regulation",
               trend_classes=c('UP', 'DOWN'))
  }

  for (nm in names(logcpm_matrices)) {
    nmlcpm <- logcpm_matrices[[nm]]
    cat(paste0('Plotting raw heatmap with ', nm, ' cpms\n'))
    
    png(filename=paste(args$OUT_FOLDER, paste('voom_', cname, '_', gsub('.{6}$', '', nm), '_heatmap.png',sep=''), sep='/'), width = 20, height = 15, units = "in", res = 150)
    coolmap.2(nmlcpm[top_n,], 
              ColSideColors=labs, 
              margins=c(12, 25),
              cbar_extreme=2)  
    legend("topright",
           legend=legend.text,
           fill=legend.fill,
           border=FALSE,
           bty="n",
           y.intersp = 0.7,
           cex=0.7)
    dev.off()
  } # for (nm in names(logcpm_matrices))
} # for (cname in cnames)

