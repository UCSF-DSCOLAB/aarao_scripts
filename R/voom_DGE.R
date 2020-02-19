#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(ggfortify)
  library(gridExtra)
  library(limma)
  library(methods)
  library(RColorBrewer)
  library(STRINGdb)
})

# poor vs rich
# logFc < 0 -> poor/rich < 0 -> upregulated in poor
# logFc > 0 -> poor/rich > 0 -> downregulated in poor


string_folder <- paste0(Sys.getenv('KLAB'), '/ipi/data/databases/string')
source(paste0(Sys.getenv('SCRIPTS'), '/R/coolmap2_heatmap3.R'))
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2 || length(args) > 5) {
  stop("Need 2 or more arguments <COUNTS.TSV>, <CLASSES.TSV>, [<OUT_FOLDER>, <BATCH_CORRECT>, <ANNOT_COLS.TSV>] in that order", call.=FALSE)
} else if (length(args) == 2) {
  args[3] <- "."
  batchCorrect <- FALSE
  annot_cols <- NULL
} else if (length(args) == 3) {
    args[3] <- gsub('/$', '', args[3])
    batchCorrect <- FALSE
    annot_cols <- NULL
} else if (length(args) == 4) {
    args[3] <- gsub('/$', '', args[3])
    batchCorrect <- as.logical(args[4])
    annot_cols <- NULL
} else {
  args[3] <- gsub('/$', '', args[3])
  batchCorrect <- as.logical(args[4])
  annot_cols <- args[5]
}


print ("Reading counts.....")
counts <- read.table(args[1], sep='\t', header=1, row.names=1)
print ("Done")

print ("Reading classes.....")
groups <- read.table(args[2], sep='\t', header=1, row.names=1, colClasses=c('factor'))
print ("Done")

if (!is.null(annot_cols)){
  print ("Reading Additional annotation columns.....")
  annot_cols <- read.table(annot_cols, sep='\t', header=1, row.names=1, stringsAsFactors = FALSE)
  print ("Additional annotation columns read")
} else {
  annot_cols <- data.frame()
}

#if (length(groups) != 1 || colnames(groups) != 'group') {
#    stop("<GROUPS.TSV> must have only one column named group", call.=FALSE)
#}


suppressWarnings({
    if (any(sapply(levels(groups$plate), function(x){!is.na(as.numeric(x))}))){
        stop("Plates in <GROUPS.TSV> cannot be numbers. Consider renaming `1` (say) to `group_1`", call.=FALSE)
    }
})

suppressWarnings({
    if (any(sapply(levels(groups$group), function(x){!is.na(as.numeric(x))}))){
        stop("Groups in <GROUPS.TSV> cannot be numbers. Consider renaming `1` (say) to `plate1`", call.=FALSE)
    }
})


print ("Populating String DB.....")
string_db <- STRINGdb$new(version="10", species=9606, score_threshold=0, input_directory=string_folder)
temp <- string_db$get_annotations()
print ("Done")


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
  warning("The following fields in <ANNOT_COLS.TSV> are not in <GROUPS.TSV>. Discarding missing columns", call.=FALSE)
  cc <- rownames(annot_cols)[rownames(annot_cols) %in% colnames(groups)]
  annot_cols <- annot_cols[cc, , drop=F]
}


y = DGEList(counts=counts, group=groups$group)
cpm <- cpm(y)
lcpm <- cpm(y, log=TRUE)
print ("made DGElist...")
y = y[ rowSums(y$counts) > 60, keep.lib.size=FALSE ]
y <- calcNormFactors(y, method="TMM")
cpm2 <- cpm(y)
lcpm2 <- cpm(y, log=TRUE)

print ("Normalization done")

# Add batch annotations into the groups dataframe
fixed_levels <- list()
fixed_levels[['sequencer']] <- c('hiseq', 'novaseq')
groups$sequencer <- factor(x=sapply(groups$plate, function(x){if (as.numeric(gsub('plate', '', x)) > 25) 'novaseq' else 'hiseq'}), levels=fixed_levels[['sequencer']])
fixed_levels[['washes']] <- c('extraWash', 'noExtraWash')
groups$washes <- factor(x=sapply(groups$plate, function(x){if (as.numeric(gsub('plate', '', x)) > 8) 'extraWash' else 'noExtraWash'}), levels=fixed_levels[['washes']])

# histology may or may not be known
if (!'histology' %in% colnames(groups)){
  groups$histology <- 'unknown'
}

# Add proper levels to EHK and plates
groups$EHK <- factor(as.vector(groups$EHK),
                     levels = as.vector(c(1:10)))


fixed_levels[['plate']] <- as.vector(levels(groups$plate))
fixed_levels[['plate']] <- fixed_levels[['plate']][order(as.numeric(gsub("plate", "", fixed_levels[['plate']])))]
groups$plate <- factor(as.vector(groups$plate), levels = fixed_levels[['plate']])

print('The plate distribution across plates is seen to be:')
print(table(groups[,c('sequencer', 'plate')]))
print('The wash distribution across plates is seen to be:')
print(table(groups[,c('washes', 'plate')]))

useSequencerToBC <- FALSE

# See https://support.bioconductor.org/p/36029/#100297
# This is how Gordon Smyth suggests you handle batch effect in the DE
design_formula = '~0+group'
covariates = c()

if (batchCorrect) {
  if (!'plate' %in% colnames(groups)) {
    stop('Cannot run a DGE without plates to batch correct on.')
  }
  design_formula <- paste0(design_formula, '+plate')
  covariates = c(covariates, 'plates')
  useSequencerToBC <- FALSE
  temp <- colSums(apply(table(groups[,c('sequencer', 'plate')]), 1, as.logical))
  if (any(temp < 2)){
    print(paste0('Even though multiple sequencers were detected, they won\'t ',
                 'be used for batch correction since it only corresponds to ',
                 'one plate or fewer.'))
  } else {
    covariates = c(covariates, 'sequencer')
    design_formula <- paste0(design_formula, '+sequencer')
    useSequencerToBC <- TRUE
  }
  #temp <- colSums(apply(table(groups[,c('washes', 'plate')]), 1, as.logical))
  #if (any(temp < 2)){
  #  print(paste0('Even though multiple sequencers were detected, they won\'t ',
  #               'be used for batch correction since it only corresponds to ',
  #               'one plate or fewer.'))
  #} else {
  #  covariates = c(covariates, 'sequencer')
  #  design_formula <- paste0(design_formula, '+sequencer')
  #  useSequencerToBC <- TRUE
  #}
  print(paste0('Batch correcting with [', paste(covariates, collapse=', ') ,'] as covariates'))
} else {
  print('Foregoing Batch correction')
}

design <- model.matrix(as.formula(design_formula), data=groups) #limma

colnames(design) <- gsub("washes", "", gsub("sequencer", "",  gsub("plateplate","plate", gsub("group", "", colnames(design)))))
print ("Made Design Model")
v <- voom(y, design)
fit <- lmFit(v, design)
print ("Limma fit done")

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
print ("Made contrast matrix")  #limma

fit2 <- contrasts.fit(fit, cont_matrix)  #limma
fit2 <- eBayes(fit2, robust=TRUE)  #limma
print ("Fit2 done")

pal <- colorRampPalette(c("light grey", "Dark Blue"))

cols <- list(
  group = brewer.pal(max(length(levels(groups$group)), 3), "Accent"),
  plate = rainbow(length(levels(groups$plate))),
  EHK = pal(10),
  washes = c('black', 'lightgrey'),
  sequencer = c('black', 'lightgrey'),
  histology = brewer.pal(max(length(levels(groups$histology)), 3), "Spectral")
)

for (annot in rownames(annot_cols)){
  if (annot %in% names(cols)){
    next
  }
  if (annot_cols[annot, 'cmap'] == 'rainbow'){
    cols[[annot]] = rainbow(length(levels(groups[[annot]])))  
  } else {
    cols[[annot]] = brewer.pal(max(length(levels(groups[[annot]])), 3), annot_cols[annot, 'cmap'])  
  }
}


labs <- c()
labs.colnames <- c()

legend.text <- c()
legend.fill <- c()

for (annot_group in names(cols)){
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
}

colnames(labs) <- labs.colnames

design <- model.matrix(~0+group, data=groups) #limma
colnames(design) <- gsub("group", "", colnames(design))
hclust.ward = function(d) hclust(d,method="ward.D2")

for (cname in cnames){
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
  write.table(top_table, file=paste(args[3], paste('voom_', cname, '.tsv',sep=''), sep='/'), sep='\t', row.names=T, col.names=T, quote=F)
  p <- ggplot(top_table, aes(logFC, -log10(P.Value))) + geom_point(aes(col=col)) + scale_color_manual(values=c("black", "red", "red")) + theme(legend.position = "none") + geom_text_repel(data=filter(top_table, col=='2'), aes(label=gene), size=2, color='black')
  png(filename=paste(args[3], paste('voom_', cname, '.png',sep=''), sep='/'), width = 10, height = 10, units = "in", res = 150)
  print(p)
  dev.off()
  
  # Do the string stuff
  top_table_mapped <- string_db$map(top_table, "gene", removeUnmappedRows = TRUE )
  hits <- top_table_mapped$STRING_id[1:200]
  top_table_mapped_pval05 <- string_db$add_diff_exp_color( subset(top_table_mapped, P.Value<0.05), logFcColStr="logFC" )
  # post payload information to the STRING server
  payload_id <- string_db$post_payload( top_table_mapped_pval05$STRING_id, colors=top_table_mapped_pval05$color )
  # display a STRING network png with the "halo"
  png(filename=paste(args[3], paste('voom_', cname, '_top200STRINGdb.png',sep=''), sep='/'), width = 10, height = 10, units = "in", res = 150)
  string_db$plot_network( hits, payload_id=payload_id )
  dev.off()
  
  hits <- top_table_mapped_pval05$STRING_id
  enrichmentGO <- string_db$get_enrichment( hits, category = "Process", methodMT = "fdr", iea = TRUE )
  enrichmentKEGG <- string_db$get_enrichment( hits, category = "KEGG", methodMT = "fdr", iea = TRUE )
  write.table(enrichmentGO, file=paste(args[3], paste('voom_', cname, '_GOenrichment.tsv',sep=''), sep='/'), sep='\t', row.names=F, col.names=T, quote=F)
  write.table(enrichmentKEGG, file=paste(args[3], paste('voom_', cname, '_KEGGenrichment.tsv',sep=''), sep='/'), sep='\t', row.names=F, col.names=T, quote=F)
  
  enrichmentGO20 <- head(enrichmentGO, n=20) 
  top_annotations <- string_db$annotations[string_db$annotations$term_id %in% enrichmentGO20$term_id,]
  top_annotations <- top_annotations[top_annotations$STRING_id %in% top_table_mapped_pval05$STRING_id,]
  top_annotations$gene <- sapply(top_annotations$STRING_id, function(x) {paste(top_table_mapped_pval05[top_table_mapped_pval05$STRING_id==x, 'gene'], collapse=',')})
  top_annotations$description <- sapply(top_annotations$term_id, function(x) {enrichmentGO20[enrichmentGO20$term_id==x, 'term_description']})
  top_annotations$regulation <- sapply(top_annotations$gene, function(x) {ifelse(!(x %in% top_table_mapped_pval05$gene), 'UNKNOWN', ifelse(top_table_mapped_pval05[top_table_mapped_pval05$gene==x, 'logFC'] > 0, 'UP', 'DOWN'))})
  top_annotations <- top_annotations[, c('gene', 'regulation', 'term_id', 'description')]
  top_annotations <- top_annotations[order(top_annotations$term_id, top_annotations$regulation),]
  write.table(top_annotations, file=paste(args[3], paste('voom_', cname, '_GOenrichment_genes.tsv',sep=''), sep='/'), sep='\t', row.names=F, col.names=T, quote=F)

  enrichmentKEGG20 <- head(enrichmentKEGG, n=20) 
  top_annotations <- string_db$annotations[string_db$annotations$term_id %in% enrichmentKEGG20$term_id,]
  top_annotations <- top_annotations[top_annotations$STRING_id %in% top_table_mapped_pval05$STRING_id,]
  top_annotations$gene <- sapply(top_annotations$STRING_id, function(x) {paste(top_table_mapped_pval05[top_table_mapped_pval05$STRING_id==x, 'gene'], collapse=',')})
  top_annotations$description <- sapply(top_annotations$term_id, function(x) {enrichmentKEGG20[enrichmentKEGG20$term_id==x, 'term_description']})
  top_annotations$regulation <- sapply(top_annotations$gene, function(x) {ifelse(!(x %in% top_table_mapped_pval05$gene), 'UNKNOWN', ifelse(top_table_mapped_pval05[top_table_mapped_pval05$gene==x, 'logFC'] > 0, 'UP', 'DOWN'))})
  top_annotations <- top_annotations[, c('gene', 'regulation', 'term_id', 'description')]
  top_annotations <- top_annotations[order(top_annotations$term_id, top_annotations$regulation),]
  write.table(top_annotations, file=paste(args[3], paste('voom_', cname, '_KEGGenrichment_genes.tsv',sep=''), sep='/'), sep='\t', row.names=F, col.names=T, quote=F)
  
  if (batchCorrect){
    print('Creating a batch corrected matrix for PCA and heatmaps')
    if (useSequencerToBC){
      x <- removeBatchEffect(lcpm, batch=groups$plate, batch2=groups$sequencer, design=design)
    } else {
      x <- removeBatchEffect(lcpm, batch=groups$plate, design=design)
    }
    print('Plotting Batch corrected heatmap')
    png(filename=paste(args[3], paste('voom_', cname, '_bc_heatmap.png',sep=''), sep='/'), width = 10, height = 10, units = "in", res = 150)
    #coolmap.2(x[top_n,], linkage.row='ward', linkage.col='ward', ColSideColors=labs, margins=c(12, 12))
    #coolmap.2(x[top_n,], hclust=hclust.ward, ColSideColors=labs, margins=c(12, 12))
    coolmap.2(x[top_n,], ColSideColors=labs, margins=c(12, 12))
    legend("topright",
           legend=c(as.character(levels(groups$group)), "", as.character(levels(groups$plate)), "", c('hiseq', 'novaseq')),
           fill=c(group_cols[as.factor(levels(groups$group))], 'white', plate_cols[as.factor(levels(groups$plate))], 'white', 'black', 'lightgrey'),
           border=FALSE,
           bty="n",
           y.intersp = 0.7,
           cex=0.7)
    dev.off()
    print('Plotting Batch corrected pca')
    png(filename=paste(args[3], paste('voom_', cname, '_bc_pca.png',sep=''), sep='/'), width = 10, height = 10, units = "in", res = 150)
    #ap <- autoplot(prcomp(t(x), center=TRUE, scale.=TRUE), data=groups, shape='group', col='plate', size=5, alpha=0.5)
    ap <- autoplot(prcomp(t(x)), data=groups, shape='group', col='plate', size=5, alpha=0.5)
    print(ap)
    dev.off()
  }
  
  for (nm in c('normalized', 'unnormalized')){
    print(paste0('Plotting raw heatmap with ', nm, ' cpms'))
    
    if (nm == 'normalized'){
      nmlcpm <- lcpm2
    } else {
      nmlcpm <- lcpm
    }
    
    png(filename=paste(args[3], paste('voom_', cname, '_nobc_', gsub('.{6}$', '', nm), '_heatmap.png',sep=''), sep='/'), width = 20, height = 15, units = "in", res = 150)
    coolmap.2(nmlcpm[top_n,], ColSideColors=labs, margins=c(12, 25))  
    legend("topright",
           legend=legend.text,
           fill=legend.fill,
           border=FALSE,
           bty="n",
           y.intersp = 0.7,
           cex=0.7)
    dev.off()
    
    print(paste0('Plotting raw pca with ', nm, ' cpms'))
    pcs <- prcomp(t(nmlcpm))
    ap <- list()
    for (annot_name in names(cols)){
      if (annot_name == 'group'){
        ap[[annot_name]] <- autoplot(pcs, data=groups, shape=annot_name, col=annot_name, size=5, alpha=0.75) + geom_text_repel(aes(label=rownames(groups))) + scale_color_manual(values=cols[[annot_name]])
      } else{
        ap[[annot_name]] <- autoplot(pcs, data=groups, shape='group', col=annot_name, size=5, alpha=0.75) + scale_color_manual(values=cols[[annot_name]])
      }
    }
    
    png(filename=paste(args[3], paste('voom_', cname, '_nobc_', gsub('.{6}$', '', nm), '_pca.png',sep=''), sep='/'), width = 30, height =  ceiling(length(ap)/3) *10 , units = "in", res = 150)
    grid.arrange(grobs=ap, ncol = 3)
    dev.off()
    
    write_text <- c()
    for (pc in c(1, 2)){
      top_most_variable <- names(sort(apply(nmlcpm, 1, var), decreasing = TRUE)[c(1:5000)])
      top_most_variable <- t(lcpm2[top_most_variable, ])
      best_correlation <- sort(apply(top_most_variable, 2, function(x) {cor(x, pcs$x[,pc])}))
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
    }
    fileConn<-file(paste(args[3], paste('voom_', cname, '_nobc_', gsub('.{6}$', '', nm), '_pca_correlating_genes.txt',sep=''), sep='/'))
    writeLines(write_text, fileConn)
    close(fileConn)
  }
}

