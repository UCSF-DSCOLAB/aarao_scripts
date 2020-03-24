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
  library(STRINGdb)
})

# poor vs rich
# logFc < 0 -> poor/rich < 0 -> upregulated in poor
# logFc > 0 -> poor/rich > 0 -> downregulated in poor


args = list(
  COUNTS_TSV=NA,  # Directory with folders containing fastqs
  CLASSES_TSV=NA,  # Samples to pull from the directory
  OUT_FOLDER=".",  # Where to write results
  BATCH_CORRECT=FALSE,
  ANNOT_COLS_TSV=NULL,
  STRING_FOLDER=paste(Sys.getenv("KLAB"), "/ipi/data/databases/string", sep=""),
  RSCRIPTS_DIR=paste(Sys.getenv("SCRIPTS"), "R", sep="/"),  # Where to find aux scripts
  RUN_STRING_DE=TRUE,
  RUN_STRING_PCA=TRUE,
  MOST_VARIABLE_CUTOFF=1000
  )

argsClasses = list(
  COUNTS_TSV=as.character,
  CLASSES_TSV=as.character,
  OUT_FOLDER=as.character,
  BATCH_CORRECT=as.logical,
  ANNOT_COLS_TSV=as.character,
  STRING_FOLDER=as.character,
  RSCRIPTS_DIR=as.character,
  RUN_STRING_DE=as.logical,
  RUN_STRING_PCA=as.logical,
  MOST_VARIABLE_CUTOFF=as.numeric
)

UserArgs = commandArgs(trailingOnly=TRUE)

for (i in UserArgs){
  key_val = strsplit(i, split="=", fixed=TRUE)[[1]]
  assert_that(length(key_val) == 2, msg=paste0("Invalid option `", i, "`. ",
                                               "Args must be of the form ",
                                               "KEY=VAL ."))
  if (!key_val[1] %in% names(args)){
    print(paste0("WARNING: Invalid key provided : ", i, ""))
    next
  } else {
    args[[key_val[1]]] = argsClasses[[key_val[1]]](key_val[2])
  }
}

print("Running with arguments:")
print(paste0("COUNTS_TSV : ", args$COUNTS_TSV, " (", mode(args$COUNTS_TSV), ")"))
print(paste0("CLASSES_TSV : ", args$CLASSES_TSV, " (", mode(args$CLASSES_TSV), ")"))
print(paste0("OUT_FOLDER  : ", args$OUT_FOLDER, " (", mode(args$OUT_FOLDER), ")"))
print(paste0("BATCH_CORRECT : ", args$BATCH_CORRECT, " (", mode(args$BATCH_CORRECT), ")"))
print(paste0("ANNOT_COLS_TSV : ", args$ANNOT_COLS_TSV, " (", mode(args$ANNOT_COLS_TSV), ")"))
print(paste0("STRING_FOLDER : ", args$STRING_FOLDER, " (", mode(args$STRING_FOLDER), ")"))
print(paste0("RSCRIPTS_DIR : ", args$RSCRIPTS_DIR, " (", mode(args$RSCRIPTS_DIR), ")"))
print(paste0("RUN_STRING_DE : ", args$RUN_STRING_DE, " (", mode(args$RUN_STRING_DE), ")"))
print(paste0("RUN_STRING_PCA : ", args$RUN_STRING_PCA, " (", mode(args$RUN_STRING_PCA), ")"))
print(paste0("MOST_VARIABLE_CUTOFF : ", args$MOST_VARIABLE_CUTOFF, " (", mode(args$MOST_VARIABLE_CUTOFF), ")"))


if (is.na(args$COUNTS_TSV) | is.na(args$CLASSES_TSV)) {
  stop(paste0("Need COUNTS_TSV=<COUNTS.tsv> and CLASSES_TSV=<CLASSES.tsv> to continue. ",
              "Other options include OUT_FOLDER, BATCH_CORRECT, ANNOT_COLS_TSV, ",
              "RSCRIPTS_DIR, and STRING_FOLDER ."),
       call.=FALSE)
}

source(paste0(args$RSCRIPTS_DIR, '/coolmap2_heatmap3.R'))

suppmsg <- assert_that(is.logical(args$BATCH_CORRECT), msg="BATCH_CORRECT must be one of TRUE or FALSE")
suppmsg <- assert_that(is.logical(args$RUN_STRING_DE), msg="RUN_STRING_DE must be one of TRUE or FALSE")
suppmsg <- assert_that(is.logical(args$RUN_STRING_PCA), msg="RUN_STRING_PCA must be one of TRUE or FALSE")

print ("Reading counts.....")
counts <- read.table(args$COUNTS_TSV, sep='\t', header=1, row.names=1)
print ("Done")

print ("Reading classes.....")
groups <- read.table(args$CLASSES_TSV, sep='\t', header=1, row.names=1, stringsAsFactors = T)
print ("Done")

if (!is.null(args$ANNOT_COLS_TSV)){
  print ("Reading Additional annotation columns.....")
  annot_cols <- read.table(args$ANNOT_COLS_TSV, sep='\t', header=1, row.names=1, stringsAsFactors = FALSE)
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


if (args$RUN_STRING_DE|args$RUN_STRING_PCA){
  print ("Populating String DB.....")
  string_db <- STRINGdb$new(version="10", species=9606, score_threshold=0, input_directory=args$STRING_FOLDER)
  temp <- string_db$get_annotations()
  print ("Done")
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
  groups$histology <- as.factor('unknown')
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

if (args$BATCH_CORRECT) {
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
} # if (args$BATCH_CORRECT)/else

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

temp <- sort(fit2$sigma, decreasing = T)[args$MOST_VARIABLE_CUTOFF]
top_most_variable_genes = rownames(fit2$coefficients)[fit2$sigma >= temp]

pal <- colorRampPalette(c("light grey", "Dark Blue"))

cols <- list(
  group = brewer.pal(max(length(levels(groups$group)), 3), "Accent"),
  plate = rainbow(length(levels(groups$plate))),
  EHK = pal(10),
  washes = c('black', 'lightgrey'),
  sequencer = c('black', 'lightgrey'),
  histology = brewer.pal(max(length(levels(groups$histology)), 3), "Spectral")
)

num_cols <- list()

for (annot in rownames(annot_cols)) {
  if (annot %in% names(cols)){
    next
  }
  suppmsg <- assert_that(annot_cols[annot, 'cmap'] %in% c('rainbow', rownames(brewer.pal.info)))
  
  if (annot_cols[annot, 'type'] == "factor"){
    num_vals <- length(levels(groups[[annot]]))
  } else if (annot_cols[annot, 'type'] == "numeric") {
    
    num_cols[[annot]] <- list(palette=annot_cols[annot, 'cmap'])
    if (endsWith(annot, "_frac")){
      num_cols[[annot]][["low"]] <- 0
      num_cols[[annot]][["high"]] <- 1
    } else if (endsWith(annot, "_pct")){
      num_cols[[annot]][["low"]] <- 0
      num_cols[[annot]][["high"]] <- 100
    } else {
      num_cols[[annot]][["low"]] <- min(groups[[annot]])
      num_cols[[annot]][["high"]] <- max(groups[[annot]])
    }
    num_vals <- 10
  } else{
    assert_that(FALSE, msg="Cannot handle non-factor/numeric types for annot right now.")
  }
  
  if (annot_cols[annot, 'cmap'] == 'rainbow'){
    cols[[annot]] = rainbow(num_vals)  
  } else {
    cols[[annot]] = brewer.pal(max(num_vals, 3), annot_cols[annot, 'cmap'])  
  }
} # for (annot in rownames(annot_cols))

labs <- c()
labs.colnames <- c()

legend.text <- c()
legend.fill <- c()

for (annot_group in names(cols)) {
  if (annot_group %in% names(num_cols)) {
    pct <- (num_cols[[annot_group]][['high']]-num_cols[[annot_group]][['low']])/10
    vals <- sapply(groups[[annot_group]], function(x){min((x-num_cols[[annot_group]][['low']])%/%pct+1, 10)})
    labs <- cbind(labs, cols[[annot_group]][vals])
    labs.colnames <- c(labs.colnames, annot_group)
    
    legend.text <- c(legend.text, c("low", "medium", "high"), "")
    legend.fill <- c(legend.fill, cols[[annot_group]][c(1,5, 10)], 'white')
    
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
  } # if (annot_group %in% names(num_cols)) / else
} # for (annot_group in names(cols)) 

colnames(labs) <- labs.colnames
rownames(labs) <- rownames(groups)

design <- model.matrix(~0+group, data=groups) #limma
colnames(design) <- gsub("group", "", colnames(design))
hclust.ward = function(d) hclust(d,method="ward.D2")

string_processing_DE <- function(cname, top_table, string_db) {
  # Do the string stuff
  top_table_mapped <- string_db$map(top_table, "gene", removeUnmappedRows = TRUE )
  dropped_genes <- top_table[!top_table$gene%in%top_table_mapped$gene, 'gene']
  print("The following genes were dropped from the STRING analysis:")
  cat(dropped_genes, sep="\n")
  hits <- top_table_mapped$STRING_id[1:200]
  top_table_mapped_pval05 <- string_db$add_diff_exp_color( subset(top_table_mapped, P.Value<0.05), logFcColStr="logFC" )
  # post payload information to the STRING server
  payload_id <- string_db$post_payload( top_table_mapped_pval05$STRING_id, colors=top_table_mapped_pval05$color )
  # display a STRING network png with the "halo"
  png(filename=paste(args$OUT_FOLDER, paste('voom_', cname, '_DE_top200STRINGdb.png',sep=''), sep='/'), width = 10, height = 10, units = "in", res = 150)
  string_db$plot_network( hits, payload_id=payload_id )
  dev.off()
  
  hits <- top_table_mapped_pval05$STRING_id
  enrichmentGO <- string_db$get_enrichment( hits, category = "Process", methodMT = "fdr", iea = TRUE )
  enrichmentKEGG <- string_db$get_enrichment( hits, category = "KEGG", methodMT = "fdr", iea = TRUE )
  write.table(enrichmentGO, file=paste(args$OUT_FOLDER, paste('voom_', cname, '_DE_GOenrichment.tsv',sep=''), sep='/'), sep='\t', row.names=F, col.names=T, quote=F)
  write.table(enrichmentKEGG, file=paste(args$OUT_FOLDER, paste('voom_', cname, '_DE_KEGGenrichment.tsv',sep=''), sep='/'), sep='\t', row.names=F, col.names=T, quote=F)
  
  enrichmentGO20 <- head(enrichmentGO, n=20) 
  top_annotations <- string_db$annotations[string_db$annotations$term_id %in% enrichmentGO20$term_id,]
  top_annotations <- top_annotations[top_annotations$STRING_id %in% top_table_mapped_pval05$STRING_id,]
  top_annotations$gene <- sapply(top_annotations$STRING_id, function(x) {paste(top_table_mapped_pval05[top_table_mapped_pval05$STRING_id==x, 'gene'], collapse=',')})
  top_annotations$description <- sapply(top_annotations$term_id, function(x) {enrichmentGO20[enrichmentGO20$term_id==x, 'term_description']})
  top_annotations$regulation <- sapply(top_annotations$gene, function(x) {ifelse(!(x %in% top_table_mapped_pval05$gene), 'UNKNOWN', ifelse(top_table_mapped_pval05[top_table_mapped_pval05$gene==x, 'logFC'] > 0, 'UP', 'DOWN'))})
  top_annotations <- top_annotations[, c('gene', 'regulation', 'term_id', 'description')]
  top_annotations <- top_annotations[order(top_annotations$term_id, top_annotations$regulation),]
  write.table(top_annotations, file=paste(args$OUT_FOLDER, paste('voom_', cname, '_DE_GOenrichment_genes.tsv',sep=''), sep='/'), sep='\t', row.names=F, col.names=T, quote=F)

  enrichmentKEGG20 <- head(enrichmentKEGG, n=20) 
  top_annotations <- string_db$annotations[string_db$annotations$term_id %in% enrichmentKEGG20$term_id,]
  top_annotations <- top_annotations[top_annotations$STRING_id %in% top_table_mapped_pval05$STRING_id,]
  top_annotations$gene <- sapply(top_annotations$STRING_id, function(x) {paste(top_table_mapped_pval05[top_table_mapped_pval05$STRING_id==x, 'gene'], collapse=',')})
  top_annotations$description <- sapply(top_annotations$term_id, function(x) {enrichmentKEGG20[enrichmentKEGG20$term_id==x, 'term_description']})
  top_annotations$regulation <- sapply(top_annotations$gene, function(x) {ifelse(!(x %in% top_table_mapped_pval05$gene), 'UNKNOWN', ifelse(top_table_mapped_pval05[top_table_mapped_pval05$gene==x, 'logFC'] > 0, 'UP', 'DOWN'))})
  top_annotations <- top_annotations[, c('gene', 'regulation', 'term_id', 'description')]
  top_annotations <- top_annotations[order(top_annotations$term_id, top_annotations$regulation),]
  write.table(top_annotations, file=paste(args$OUT_FOLDER, paste('voom_', cname, '_DE_KEGGenrichment_genes.tsv',sep=''), sep='/'), sep='\t', row.names=F, col.names=T, quote=F)
} # string_processing_DE


string_processing_PCA <- function(cname, pca_weights_table, pc, string_db) {
  # Do the string stuff
  pca_weights_table_mapped <- string_db$map(pca_weights_table, "gene", removeUnmappedRows = TRUE )
  dropped_genes <- pca_weights_table[!pca_weights_table$gene%in%pca_weights_table_mapped$gene, 'gene']
  print("The following genes were dropped from the STRING analysis:")
  cat(dropped_genes, sep="\n")
  hits <- pca_weights_table_mapped$STRING_id
  pca_weights_table_mapped <- string_db$
    add_diff_exp_color(pca_weights_table_mapped, logFcColStr="weight" )
  # post payload information to the STRING server
  payload_id <- string_db$post_payload( pca_weights_table_mapped$STRING_id, colors=pca_weights_table_mapped$color )
  # display a STRING network png with the "halo"
  png(filename=paste(args$OUT_FOLDER, paste('voom_', cname, '_PC', pc, '_top200STRINGdb.png',sep=''), sep='/'), width = 10, height = 10, units = "in", res = 150)
  string_db$plot_network( hits, payload_id=payload_id )
  dev.off()
  
  enrichmentGO <- string_db$get_enrichment( hits, category = "Process", methodMT = "fdr", iea = TRUE )
  enrichmentKEGG <- string_db$get_enrichment( hits, category = "KEGG", methodMT = "fdr", iea = TRUE )
  write.table(enrichmentGO, file=paste(args$OUT_FOLDER, paste('voom_', cname, '_PC', pc, '_GOenrichment.tsv',sep=''), sep='/'), sep='\t', row.names=F, col.names=T, quote=F)
  write.table(enrichmentKEGG, file=paste(args$OUT_FOLDER, paste('voom_', cname, '_PC', pc, '_KEGGenrichment.tsv',sep=''), sep='/'), sep='\t', row.names=F, col.names=T, quote=F)
  
  enrichmentGO20 <- head(enrichmentGO, n=20) 
  top_annotations <- string_db$annotations[string_db$annotations$term_id %in% enrichmentGO20$term_id,]
  top_annotations <- top_annotations[top_annotations$STRING_id %in% pca_weights_table_mapped$STRING_id,]
  top_annotations$gene <- sapply(top_annotations$STRING_id, function(x) {paste(pca_weights_table_mapped[pca_weights_table_mapped$STRING_id==x, 'gene'], collapse=',')})
  top_annotations$description <- sapply(top_annotations$term_id, function(x) {enrichmentGO20[enrichmentGO20$term_id==x, 'term_description']})
  top_annotations$weightage <- sapply(top_annotations$gene, function(x) {ifelse(!(x %in% pca_weights_table_mapped$gene), 'UNKNOWN', ifelse(pca_weights_table_mapped[pca_weights_table_mapped$gene==x, 'weight'] > 0, 'POS', 'NEG'))})
  top_annotations <- top_annotations[, c('gene', 'weightage', 'term_id', 'description')]
  top_annotations <- top_annotations[order(top_annotations$term_id, top_annotations$weightage),]
  write.table(top_annotations, file=paste(args$OUT_FOLDER, paste('voom_', cname, '_PC', pc, '_GOenrichment_genes.tsv',sep=''), sep='/'), sep='\t', row.names=F, col.names=T, quote=F)
  
  enrichmentKEGG20 <- head(enrichmentKEGG, n=20) 
  top_annotations <- string_db$annotations[string_db$annotations$term_id %in% enrichmentKEGG20$term_id,]
  top_annotations <- top_annotations[top_annotations$STRING_id %in% pca_weights_table_mapped$STRING_id,]
  top_annotations$gene <- sapply(top_annotations$STRING_id, function(x) {paste(pca_weights_table_mapped[pca_weights_table_mapped$STRING_id==x, 'gene'], collapse=',')})
  top_annotations$description <- sapply(top_annotations$term_id, function(x) {enrichmentKEGG20[enrichmentKEGG20$term_id==x, 'term_description']})
  top_annotations$weightage <- sapply(top_annotations$gene, function(x) {ifelse(!(x %in% pca_weights_table_mapped$gene), 'UNKNOWN', ifelse(pca_weights_table_mapped[pca_weights_table_mapped$gene==x, 'weight'] > 0, 'POS', 'NEG'))})
  top_annotations <- top_annotations[, c('gene', 'weightage', 'term_id', 'description')]
  top_annotations <- top_annotations[order(top_annotations$term_id, top_annotations$weightage),]
  write.table(top_annotations, file=paste(args$OUT_FOLDER, paste('voom_', cname, '_PC', pc, '_KEGGenrichment_genes.tsv',sep=''), sep='/'), sep='\t', row.names=F, col.names=T, quote=F)
} # string_processing_PCA


for (nm in c('normalized', 'unnormalized')) {
  print(paste0('Plotting unbiased heatmap with ', nm, ' cpms'))
  
  if (nm == 'normalized'){
    nmlcpm <- lcpm2
  } else {
    nmlcpm <- lcpm
  }
  
  png(filename=paste(args$OUT_FOLDER, paste('voom_', gsub('.{6}$', '', nm), '_unbiased_top', args$MOST_VARIABLE_CUTOFF, '_heatmap.png',sep=''), sep='/'), width = 20, height = 15, units = "in", res = 150)
  coolmap.2(nmlcpm[top_most_variable_genes,], 
            ColSideColors=labs, 
            margins=c(12, 25))  
  legend("topright",
         legend=legend.text,
         fill=legend.fill,
         border=FALSE,
         bty="n",
         y.intersp = 0.7,
         cex=0.7)
  dev.off()
} # for (nm in c('normalized', 'unnormalized')) 

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
    string_processing_DE(cname, top_table, string_db)
  }

  if (args$BATCH_CORRECT) {
    print('Creating a batch corrected matrix for PCA and heatmaps')
    if (useSequencerToBC){
      x <- removeBatchEffect(lcpm, batch=groups$plate, batch2=groups$sequencer, design=design)
    } else {
      x <- removeBatchEffect(lcpm, batch=groups$plate, design=design)
    }
    print('Plotting Batch corrected heatmap')
    png(filename=paste(args$OUT_FOLDER, paste('voom_', cname, '_bc_heatmap.png',sep=''), sep='/'), width = 10, height = 10, units = "in", res = 150)
    #coolmap.2(x[top_n,], linkage.row='ward', linkage.col='ward', ColSideColors=labs, margins=c(12, 12))
    #coolmap.2(x[top_n,], hclust=hclust.ward, ColSideColors=labs, margins=c(12, 12))
    coolmap.2(x[top_n,], ColSideColors=labs, margins=c(12, 12))
    legend("topright",
           legend=c(as.character(levels(groups$group)), 
                    "", 
                    as.character(levels(groups$plate)), 
                    "", 
                    c('hiseq', 'novaseq')),
           fill=c(group_cols[as.factor(levels(groups$group))], 
                  'white', 
                  plate_cols[as.factor(levels(groups$plate))], 
                  'white', 
                  'black', 'lightgrey'),
           border=FALSE,
           bty="n",
           y.intersp = 0.7,
           cex=0.7)
    dev.off()
    print('Plotting Batch corrected pca')
    png(filename=paste(args$OUT_FOLDER, paste('voom_', cname, '_bc_pca.png',sep=''), sep='/'), width = 10, height = 10, units = "in", res = 150)
    #ap <- autoplot(prcomp(t(x), center=TRUE, scale.=TRUE), data=groups, shape='group', col='plate', size=5, alpha=0.5)
    ap <- autoplot(prcomp(t(x)), data=groups, shape='group', col='plate', size=5, alpha=0.5)
    print(ap)
    dev.off()
  } # if (args$BATCH_CORRECT)
  
  for (nm in c('unnormalized', 'normalized')) {
    print(paste0('Plotting raw heatmap with ', nm, ' cpms'))
    
    if (nm == 'normalized'){
      nmlcpm <- lcpm2
    } else {
      nmlcpm <- lcpm
    }
    
    png(filename=paste(args$OUT_FOLDER, paste('voom_', cname, '_', gsub('.{6}$', '', nm), '_heatmap.png',sep=''), sep='/'), width = 20, height = 15, units = "in", res = 150)
    coolmap.2(nmlcpm[top_n,], 
              ColSideColors=labs, 
              margins=c(12, 25))  
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
        ap[[annot_name]] <- autoplot(pcs, data=groups, shape=annot_name, col=annot_name, size=5, alpha=0.75) + 
                                geom_text_repel(aes(label=rownames(groups))) + 
                                scale_color_manual(values=cols[[annot_name]])
      } else if (annot_name %in% names(num_cols)){
        ap[[annot_name]] <- autoplot(pcs, data=groups, shape='group', col=annot_name, size=5, alpha=0.75) + 
                                scale_color_distiller(palette=num_cols[[annot_name]][["palette"]], limits = c(num_cols[[annot_name]][["low"]], num_cols[[annot_name]][["high"]]))
      } else {
        ap[[annot_name]] <- autoplot(pcs, data=groups, shape='group', col=annot_name, size=5, alpha=0.75) + 
                                scale_color_manual(values=cols[[annot_name]])
      }
    }
    
    png(filename=paste(args$OUT_FOLDER, paste('voom_', cname, '_', gsub('.{6}$', '', nm), '_pca.png',sep=''), sep='/'), width = 30, height =  ceiling(length(ap)/3) *10 , units = "in", res = 150)
    grid.arrange(grobs=ap, ncol = 3)
    dev.off()

    x <- as.data.frame(pcs$x)
    x$group = groups$group
    y <- melt(x, id.vars = 'group', variable.name = 'PC')

    png(filename=paste(args$OUT_FOLDER, paste('voom_', cname, '_', gsub('.{6}$', '', nm), '_pca_boxplots.png',sep=''), sep='/'), width = 30, height =  ceiling(length(ap)/3) *10 , units = "in", res = 150)
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
    fileConn<-file(paste(args$OUT_FOLDER, paste('voom_', cname, '_', gsub('.{6}$', '', nm), '_pca_correlating_genes.txt',sep=''), sep='/'))
    writeLines(write_text, fileConn)
    close(fileConn)
    
    for (pc in c(1, 2)) {
      x <- data.frame(weight=pcs$rotation[, pc][order(-abs(pcs$rotation[, pc]))][c(1:200)])
      x$gene <- rownames(x)
      write.table(x, 
                  file=paste(args$OUT_FOLDER, paste('voom_', cname, '_', gsub('.{6}$', '', nm), '_pca_top200_weights.tsv',sep=''), sep='/'),
                  row.names = T, 
                  col.names=T, 
                  sep="\t", 
                  quote=F)
      if (args$RUN_STRING_PCA){
        string_processing_PCA(paste(cname, gsub('.{6}$', '', nm), sep="_"), x, pc, string_db)
      }
    } 
  } # for (nm in c('unnormalized', 'normalized')) 
} # for (cname in cnames)

