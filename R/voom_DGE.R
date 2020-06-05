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

# Remember, this is poor vs rich: So
# logFc < 0 -> poor/rich < 0 -> upregulated in poor
# logFc > 0 -> poor/rich > 0 -> downregulated in poor

# TODO replace hugo and ensembl with args$SYMBOL_COL and args$GENE_ID_COL

# Colors used by hex
# #2B2D2F = Dark Grey
# #A9A9A9 = A lighter dark grey
# #FF0000 = Red
# #FFFFFF = White
# #000000 = Black

# Function definitions

args = list(
  COUNTS_TSV=NA,  # Directory with folders containing fastqs
  CLASSES_TSV=NA,  # Metadata for the samples
  COUNTS_TYPE='counts', # counts or tpm
  GROUP_COL='group', # Column to group all data by
  GENE_TO_SYMBOL_TSV=NA,
  OUT_FOLDER=".",  # Where to write results
  MODEL_COVARIATE=c(),
  ANNOT_COLS_TSV=NULL,
  STRING_FOLDER=paste(Sys.getenv("KLAB"), "/ipi/data/databases/string", sep=""),
  RSCRIPTS_DIR=paste(Sys.getenv("SCRIPTS"), "R", sep="/"),  # Where to find aux scripts
  RUN_STRING_DE=TRUE,
  RUN_STRING_PCA=TRUE,
  RUN_DAVID_DE=TRUE,
  RUN_DAVID_PCA=TRUE,
  DAVID_AUTH_EMAIL=NA,
  MOST_VARIABLE_CUTOFF=c(200, 500, 1000),
  MOST_VARIABLE_METHOD="var",
  FDR_THRESHOLD=0.05,
  GENESETS=c()
)

argsClasses = list(
  COUNTS_TSV=as.character,
  CLASSES_TSV=as.character,
  COUNTS_TYPE=as.character,
  GROUP_COL=as.character,
  GENE_TO_SYMBOL_TSV=as.character,
  OUT_FOLDER=as.character,
  MODEL_COVARIATE=as.character,
  ANNOT_COLS_TSV=as.character,
  STRING_FOLDER=as.character,
  RSCRIPTS_DIR=as.character,
  RUN_STRING_DE=as.logical,
  RUN_STRING_PCA=as.logical,
  RUN_DAVID_DE=as.logical,
  RUN_DAVID_PCA=as.logical,
  DAVID_AUTH_EMAIL=as.character,
  MOST_VARIABLE_CUTOFF=as.numeric,
  MOST_VARIABLE_METHOD=as.character,
  FDR_THRESHOLD=as.numeric,
  GENESETS=as.character
)

user_args <- commandArgs(trailingOnly=TRUE)
# user_args = strsplit("COUNTS_TSV=/Users/arjunarkalrao/projects/gyn_indication/poor_rich_dge/05_04_2020/live/live_valid_samples_counts.tsv CLASSES_TSV=/Users/arjunarkalrao/projects/gyn_indication/poor_rich_dge/05_04_2020/live/live_valid_samples.tsv OUT_FOLDER=/Users/arjunarkalrao/projects/gyn_indication/poor_rich_dge/05_04_2020/live/ ANNOT_COLS_TSV=/Users/arjunarkalrao/projects/gyn_indication/poor_rich_dge/05_04_2020/important_fields.tsv STRING_FOLDER=/Users/arjunarkalrao/important_data/databases/string RSCRIPTS_DIR=/Users/arjunarkalrao/scripts_bk/R RUN_STRING_DE=FALSE RUN_STRING_PCA=FALSE MOST_VARIABLE_CUTOFF=200 MOST_VARIABLE_CUTOFF=500 MOST_VARIABLE_CUTOFF=1000 DAVID_AUTH_EMAIL=arjunarkal.rao@ucsf.edu", " ")[[1]]
# user_args = strsplit("COUNTS_TSV=/Users/arjunarkalrao/projects/gyn_indication/poor_rich_dge/05_04_2020/tcell/tcell_valid_samples_counts.tsv CLASSES_TSV=/Users/arjunarkalrao/projects/gyn_indication/poor_rich_dge/05_04_2020/tcell/tcell_valid_samples.tsv OUT_FOLDER=/Users/arjunarkalrao/projects/gyn_indication/poor_rich_dge/05_04_2020/tcell/ ANNOT_COLS_TSV=/Users/arjunarkalrao/projects/gyn_indication/poor_rich_dge/05_04_2020/important_fields.tsv STRING_FOLDER=/Users/arjunarkalrao/important_data/databases/string RSCRIPTS_DIR=/Users/arjunarkalrao/scripts_bk/R RUN_STRING_DE=FALSE RUN_STRING_PCA=FALSE MOST_VARIABLE_CUTOFF=200 MOST_VARIABLE_CUTOFF=500 MOST_VARIABLE_CUTOFF=1000 DAVID_AUTH_EMAIL=arjunarkal.rao@ucsf.edu", " ")[[1]]
# user_args <- strsplit("COUNTS_TSV=/Users/arjunarkalrao/projects/lung_indication/risk_dge/05_27_2019/live/live_valid_samples_counts.tsv CLASSES_TSV=/Users/arjunarkalrao/projects/lung_indication/risk_dge/05_27_2019/live/live_valid_samples.tsv GENE_TO_SYMBOL_TSV=/Users/arjunarkalrao/projects/lung_indication/risk_dge/05_27_2019/live/ensembl_to_symbol.tsv GROUP_COL=risk_group OUT_FOLDER=/Users/arjunarkalrao/projects/lung_indication/risk_dge/05_27_2019/live/ STRING_FOLDER=/Users/arjunarkalrao/important_data/databases/string RSCRIPTS_DIR=/Users/arjunarkalrao/scripts_bk/R MOST_VARIABLE_CUTOFF=200 MOST_VARIABLE_CUTOFF=500 MOST_VARIABLE_CUTOFF=1000 RUN_STRING_DE=FALSE RUN_STRING_PCA=FALSE RUN_DAVID_PCA=FALSE RUN_DAVID_DE=FALSE DAVID_AUTH_EMAIL=arjunarkal.rao@ucsf.edu", " ")[[1]]

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
cat(paste0("COUNTS_TSV          : ", args$COUNTS_TSV, " (", mode(args$COUNTS_TSV), ")\n"))
cat(paste0("CLASSES_TSV         : ", args$CLASSES_TSV, " (", mode(args$CLASSES_TSV), ")\n"))
cat(paste0("COUNTS_TYPE         : ", args$COUNTS_TYPE, " (", mode(args$COUNTS_TYPE), ")\n"))
cat(paste0("GROUP_COL           : ", args$GROUP_COL, " (", mode(args$GROUP_COL), ")\n"))
cat(paste0("GENE_TO_SYMBOL_TSV  : ", args$GENE_TO_SYMBOL_TSV, " (", mode(args$GENE_TO_SYMBOL_TSV), ")\n"))
cat(paste0("OUT_FOLDER          : ", args$OUT_FOLDER, " (", mode(args$OUT_FOLDER), ")\n"))
cat(paste("MODEL_COVARIATE      : ", args$MODEL_COVARIATE, " (", mode(args$MODEL_COVARIATE), ")\n", sep='', collapse=''))
cat(paste0("ANNOT_COLS_TSV      : ", args$ANNOT_COLS_TSV, " (", mode(args$ANNOT_COLS_TSV), ")\n"))
cat(paste0("STRING_FOLDER       : ", args$STRING_FOLDER, " (", mode(args$STRING_FOLDER), ")\n"))
cat(paste0("RSCRIPTS_DIR        : ", args$RSCRIPTS_DIR, " (", mode(args$RSCRIPTS_DIR), ")\n"))
cat(paste0("RUN_STRING_DE       : ", args$RUN_STRING_DE, " (", mode(args$RUN_STRING_DE), ")\n"))
cat(paste0("RUN_STRING_PCA      : ", args$RUN_STRING_PCA, " (", mode(args$RUN_STRING_PCA), ")\n"))
cat(paste0("RUN_DAVID_DE        : ", args$RUN_DAVID_DE, " (", mode(args$RUN_DAVID_DE), ")\n"))
cat(paste0("RUN_DAVID_PCA       : ", args$RUN_DAVID_PCA, " (", mode(args$RUN_DAVID_PCA), ")\n"))
cat(paste0("DAVID_AUTH_EMAIL    : ", args$DAVID_AUTH_EMAIL, " (", mode(args$DAVID_AUTH_EMAIL), ")\n"))
cat(paste("MOST_VARIABLE_CUTOFF : ", args$MOST_VARIABLE_CUTOFF, " (", mode(args$MOST_VARIABLE_CUTOFF), ")\n", sep='', collapse=''))
cat(paste("MOST_VARIABLE_METHOD : ", args$MOST_VARIABLE_METHOD, " (", mode(args$MOST_VARIABLE_METHOD), ")\n", sep='', collapse=''))
cat(paste0("FDR_THRESHOLD       : ", args$FDR_THRESHOLD, " (", mode(args$FDR_THRESHOLD), ")\n"))
cat(paste("GENESETS             : ", args$GENESETS, " (", mode(args$GENESETS), ")\n", sep='', collapse=''))

if (endsWith(args$OUT_FOLDER, '/')) {
  args$OUT_FOLDER = gsub('/$', '', args$OUT_FOLDER)
}

if (is.na(args$COUNTS_TSV) | is.na(args$CLASSES_TSV)) {
  stop(paste0("Need COUNTS_TSV=<COUNTS.tsv> and CLASSES_TSV=<CLASSES.tsv> to continue. ",
              "Other options include OUT_FOLDER, BATCH_CORRECT, ANNOT_COLS_TSV, ",
              "RSCRIPTS_DIR, and STRING_FOLDER ."),
       call.=FALSE)
}

if (! all(args$MOST_VARIABLE_METHOD %in% c('var', 'IQR', 'IQV', 'CQV'))){
  stop("MOST_VARIABLE_METHOD must be one (or more) of `var` (standard variance), `IQR` (Inter Quartile Vange), `IQV` (Inter Quartile variance), or `CQV` (Coefficient of Quartile Variance).",
       call.=FALSE)
}


source(paste0(args$RSCRIPTS_DIR, '/coolmap2_heatmap3.R'))
source(paste0(args$RSCRIPTS_DIR, '/biomart_operations.R'))
source(paste0(args$RSCRIPTS_DIR, '/string_operations.R'))
source(paste0(args$RSCRIPTS_DIR, '/david_operations.R'))

cat ("Reading counts.....")
counts <- read.table(args$COUNTS_TSV, 
                     sep='\t', 
                     header=1, 
                     row.names=1)
cat ("Done\n")

cat ("Reading classes.....")
classes <- read.table(args$CLASSES_TSV, 
                      sep='\t', 
                      header=1, 
                      row.names=1, 
                      stringsAsFactors = T,
                      na.strings = 'unknown')
cat ("Done\n")

if (is.na(args$GENE_TO_SYMBOL_TSV)){
  cat ('Getting Gene Symbols for all genes')
  mart <- setup_biomart(biomart="ensembl", 
                        dataset = "hsapiens_gene_ensembl")
  
  ensembl_to_symbol <- data.frame(ensembl=rownames(counts), stringsAsFactors = F)
  ensembl_to_symbol <- ensg_to_hugo(mart, ensembl_to_symbol, fill = T)
} else {
  cat (paste0('Reading Gene Symbols  for all genes from ', args$GENE_TO_SYMBOL_TSV, '\n'))
  ensembl_to_symbol <- read.table(args$GENE_TO_SYMBOL_TSV,  
                                  sep='\t',
                                  header=1,
                                  stringsAsFactors = FALSE)
  if (!all(c('hugo', 'ensembl') %in% colnames(ensembl_to_symbol))){
    stop("Can't use provided GENE_TO_SYMBOL_TSV since it does not contain the columns `hugo` and `ensembl`")
  }
}
assign('ensembl_to_symbol', ensembl_to_symbol)

add_gene_symbol <- function(df, df_ensembl_col='ensembl') {
  if (is.null(ensembl_to_symbol)){
    stop("Can't run this without setting ensembl_to_symbol")
  }
  
  temp_rownames <- rownames(df)
  df$TEMP_ROW_NAMES <- rownames(df)
  df <- merge(df, 
              ensembl_to_symbol, 
              by.x=df_ensembl_col, 
              by.y='ensembl',
              all.x=T,
              all.y=FALSE)
  rownames(df) <- df$TEMP_ROW_NAMES
  df$TEMP_ROW_NAMES <- NULL
  df <- df[temp_rownames, ]
  df$hugo[is.na(df$hugo)] <- df$ensembl[is.na(df$hugo)]
  df$hugo[df$hugo==''] <- df$ensembl[df$hugo=='']
  df
}

if (!is.null(args$ANNOT_COLS_TSV)){
  cat ("Reading Annotation columns.....")
  annot_cols <- read.table(args$ANNOT_COLS_TSV, sep='\t', 
                           header=1, row.names=1, 
                           stringsAsFactors = FALSE, 
                           quote="")
  cat ("Done\n")
  if (args$GROUP_COL %in% rownames(annot_cols)){
    cat(paste0('Dropping ', args$GROUP_COL, ' from ANNOT_COLS.tsv as it is the main grouping var.\n'))
    annot_cols <- annot_cols[rownames(annot_cols)!=args$GROUP_COL, ]
  }
} else {
  annot_cols <- data.frame()
}


if (is.numeric(classes$plate)){
  cat(paste0('`plate` in <CLASSES.TSV> cannot be numerical. Renaming to ',
             'prefix with `plate` (i.e. `1` becomes `plate1`). This ',
             'is the expected format.\n'))
  classes$plate <- paste0('plate', classes$plate)
}

if (!args$GROUP_COL %in% colnames(classes)) {
  stop(paste0("<CLASSES.TSV> must have a column named ", args$GROUP_COL), call.=FALSE)
} else if (is.numeric(classes[[args$GROUP_COL]])){
  cat(paste0('`', args$GROUP_COL, '` in <CLASSES.TSV> cannot be numerical. Renaming to ',
             'prefix with `g` (i.e. `1` becomes `g1`).\n'))
  classes[[args$GROUP_COL]] <- as.factor(paste0('g', classes[[args$GROUP_COL]]))
}

if (args$RUN_STRING_DE|args$RUN_STRING_PCA){
  cat ("Populating String DB.....")
  string_db <- setup_string(args$STRING_FOLDER)
  cat ("Done\n")
}

if (args$RUN_DAVID_DE|args$RUN_DAVID_PCA){
  suppmsg <- assert_that(!is.na(args$DAVID_AUTH_EMAIL), msg = 'Need DAVID email auth to run DAVID')
  cat ("Authenticating with DAVID.....")
  david_db <- setup_david(email=args$DAVID_AUTH_EMAIL)
  cat ("Done\n")
}

cat('The distribution of the grouping variable is\n')
print(table(classes[[args$GROUP_COL]]))

# Drop any items in classes if the grouping var is NA
if (any(is.na(classes[[args$GROUP_COL]]))) {
  drop_rownames = rownames(classes)[is.na(classes[[args$GROUP_COL]])]
  warning(paste0("The following samples in <CLASSES.TSV> are NA in the ",
                 "grouping variable `", args$GROUP_COL, "`: ", 
                 paste(drop_rownames, collapse=", "),
                 " ....Discarding NA group samples"),
          call.=FALSE)
  classes <- classes[!rownames(classes) %in% drop_rownames, ]
}

# Drop any singleton samples in classes
if (any(table(classes[[args$GROUP_COL]])<2)) {
  drop_cats <- names(table(classes[[args$GROUP_COL]])[table(classes[[args$GROUP_COL]])<2])
  drop_rownames = rownames(classes)[classes[[args$GROUP_COL]]%in%drop_cats]
  warning(paste0("The following samples in <CLASSES.TSV> are singletons in the ",
                 "grouping variable `", args$GROUP_COL, "`: ", 
                 paste(drop_rownames, collapse=", "),
                 " ....Discarding NA group samples"),
          call.=FALSE)
  classes <- classes[!rownames(classes) %in% drop_rownames, ]
}
classes[[args$GROUP_COL]] <- as.factor(as.vector(classes[[args$GROUP_COL]]))

cat('The final distribution of the grouping variable is\n')
print(table(classes[[args$GROUP_COL]]))

if (length(rownames(classes)) != length(colnames(counts)) ||
    !all(rownames(classes) %in% colnames(counts)) ||
    !all(colnames(counts) %in% rownames(classes))) {
  warning("<COUNTS.TSV> and <CLASSES.TSV> do not have the same columns. Continuing with the intersection", call.=FALSE)
  common_cols <- intersect(colnames(counts), rownames(classes))
} else {
  common_cols <- rownames(classes)
}

classes <- classes[common_cols,, drop=F]
classes <- classes[order(classes[[args$GROUP_COL]], rownames(classes)),, drop=F]
counts <- counts[, rownames(classes)]


if (!all(rownames(annot_cols) %in% colnames(classes))) {
  nc <- rownames(annot_cols)[!rownames(annot_cols) %in% colnames(classes)]
  warning(paste0("The following fields in <ANNOT_COLS.TSV> are not in ",
                 "<CLASSES.TSV>: ", 
                 paste(nc, collapse=", "),
                 " ....Discarding missing columns"),
          call.=FALSE)
  cc <- rownames(annot_cols)[rownames(annot_cols) %in% colnames(classes)]
  annot_cols <- annot_cols[cc, , drop=F]
}

# Add IPI batch annotations into the classes dataframe
annotations_added <- c()
if ('plate' %in% colnames(classes)){
  fixed_levels <- list()
  # Sequencer
  fixed_levels[['sequencer']] <- c('hiseq', 'novaseq')
  classes$sequencer <- factor(x=sapply(classes$plate, function(x){if (as.numeric(gsub('plate', '', x)) > 25) 'novaseq' else 'hiseq'}), levels=fixed_levels[['sequencer']])
  # Washes
  fixed_levels[['washes']] <- c('extraWash', 'noExtraWash')
  classes$washes <- factor(x=sapply(classes$plate, function(x){if (as.numeric(gsub('plate', '', x)) > 8) 'extraWash' else 'noExtraWash'}), levels=fixed_levels[['washes']])
  # Fix Plate to be ordered correctly as a factor
  fixed_levels[['plate']] <- as.vector(levels(classes$plate))
  fixed_levels[['plate']] <- fixed_levels[['plate']][order(as.numeric(gsub("plate", "", fixed_levels[['plate']])))]
  classes$plate <- factor(as.vector(classes$plate), levels = fixed_levels[['plate']])
  
  annotations_added <- c('sequencer', 'washes')
  
  cat('\nThe plate distribution across plates is seen to be:\n')
  print(table(classes[,c('sequencer', 'plate')]))
  cat('\nThe wash distribution across plates is seen to be:\n')
  print(table(classes[,c('washes', 'plate')]))
  
}

if ('EHK' %in% colnames(classes)){
  # Add proper levels to EHK and plates
  classes$EHK <- factor(as.vector(classes$EHK),
                        levels = as.vector(c(1:10)))
}

# See https://support.bioconductor.org/p/36029/#100297
# See https://support.bioconductor.org/p/66251/#66252
# This is how Gordon Smyth suggests you handle batch effect in the DE
design_formula = paste0('~0+', args$GROUP_COL)
factor_covariates = c(args$GROUP_COL)

for (m in  args$MODEL_COVARIATE) {
  if (!m%in%colnames(classes)){
    cat(paste0('Cannot add ', m, ' to the model since it is not in <CLASSES.TSV>.\n'))
    next
  }
  
  if (length(unique(classes[[m]])) > 1) {
    if (is.factor(classes[[m]])){
      factor_covariates = c(factor_covariates, m)
    }
    design_formula <- paste0(design_formula, '+', m)
    cat(paste0('Added ', m, ' to the model.\n'))
  } else {
    cat(paste0('Covariate ', m, ' is constant across all samples. Skipping.\n'))
  }
}

design <- model.matrix(as.formula(design_formula), data=classes) #limma

for (l in levels(classes[[args$GROUP_COL]])) {
  colnames(design) <- gsub(paste0("^", args$GROUP_COL, l, "$"), l, colnames(design))
}
cat("\nUsing design Model: \n")
print(as.data.frame(design))

if (args$COUNTS_TYPE == 'tpm') {
  unnorm_string = "unnormalized"
  prior_count = 1e-6
  norm_technique = 'upperquartile'
  min_count = 2
  min_total_count = 5
  trend = TRUE
} else if (args$COUNTS_TYPE == 'counts') {
  unnorm_string = "unnormalized"
  prior_count = 2
  norm_technique = 'TMM'
  min_count = 10
  min_total_count = 15
  trend = FALSE
} else {
  cat(paste0('WARNING: Received count type: ', args$COUNTS_TYPE, '.\n',
             'Assuming data was prenormalized. Setting trend=FALSE and skipping normalization.\n'))
  norm_technique = 'skip'
  unnorm_string = args$COUNTS_TYPE
  trend = FALSE
}

logcpm_matrices <- list()
dge_list <- DGEList(counts=counts, group=classes[[args$GROUP_COL]])
cat ("Made DGElist.\n")
if (norm_technique != 'skip'){
  keep_genes <- filterByExpr(dge_list, 
                             group = dge_list$samples$group,
                             min.count=min_count,
                             min.total.count=min_total_count)
  dge_list <- dge_list[keep_genes, , keep.lib.size=FALSE ]
  logcpm_matrices[[unnorm_string]] <- cpm(dge_list, log=TRUE, prior.count=prior_count)
  dge_list <- calcNormFactors(dge_list, method=norm_technique)
  logcpm_matrices[[paste0(norm_technique, 'normalized')]] <- cpm(dge_list, log=TRUE, prior.count=prior_count)
  cat ("Normalization done.\n")
} else {
  logcpm_matrices[[unnorm_string]] <- cpm(dge_list, log=TRUE, prior.count=prior_count)  
}


v <- voom(dge_list, design)
fit <- lmFit(v, design)
cat ("Limma fit done.\n")

i=1
conts <- c()
cnames <- c()
combs <- combn(levels(classes[[args$GROUP_COL]]), 2, simplify=F)
for (comb in combs){
  conts[i] = paste(comb, collapse=" - ")
  cnames[i] = paste(comb, collapse="_vs_")
  i <- i+1
}

if (length(levels(classes[[args$GROUP_COL]])) > 2){
  lvls = levels(classes[[args$GROUP_COL]])
  for (l in lvls){
    conts[i] = paste(c(l, lvls[lvls!=l]), collapse=' - ')
    cnames[i] = paste0(l, '_vs_rest')
    i <- i+1
  }
}

cont_matrix <- makeContrasts(contrasts=conts, levels = design ) #limma
colnames(cont_matrix) <- cnames
cat ("Using contrast matrix:\n")  #limma
print(cont_matrix)

fit2 <- contrasts.fit(fit, cont_matrix)  #limma
fit2 <- eBayes(fit2, robust=TRUE, trend=trend)  #limma
cat ("Fit2 done.\n")

# Setup the colors for the plots
cols <- list()
cols[[args$GROUP_COL]] <- brewer.pal(max(length(levels(classes[[args$GROUP_COL]])), 3), "Accent")
for (l in annotations_added) {
  # Currently all our added annotations are boolean so grey/red works
  cols[[l]] <- c('#A9A9A9',     # darkgrey 
                 '#FF0000'       # red
  )
}

numeric_cols <- list()

for (annot in rownames(annot_cols)) {
  if (annot %in% names(cols)){
    next
  }
  add_UK = FALSE
  if (annot_cols[annot, 'type'] == "factor"){
    if (any(is.na(classes[[annot]]))){
      classes[[annot]] <- factor(classes[[annot]], 
                                 levels=c(levels(as.factor(classes[[annot]])), 'UNKNOWN'))
      classes[[annot]][is.na(classes[[annot]])] <- 'UNKNOWN'
      num_vals <- length(levels(classes[[annot]])) - 1
      add_UK = TRUE
    } else {
      classes[[annot]] <- as.factor(classes[[annot]])  
      num_vals <- length(levels(classes[[annot]]))
    }
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
      numeric_cols[[annot]][["low"]] <- min(classes[[annot]], na.rm = T)
      numeric_cols[[annot]][["high"]] <- max(classes[[annot]], na.rm = T)
    }
    num_vals <- 10
  } else {
    assert_that(FALSE, msg="Cannot handle non-factor/numeric types for annot right now.")
  }
  
  if(annot_cols[annot, 'cmap'] %in% rownames(brewer.pal.info)) {
    cols[[annot]] = brewer.pal(max(num_vals, 3), annot_cols[annot, 'cmap'])
    if (num_vals < length(cols[[annot]])){
      cols[[annot]] <- cols[[annot]][1:num_vals]
    }
  } else {
    # This assumes the user has passed in a function that can be evaluated to
    # generate a palette
    pal = eval(parse(text=annot_cols[annot, 'cmap']))
    cols[[annot]] = pal(num_vals)
  }
  if (add_UK){
    cols[[annot]] <- c(cols[[annot]], '#000000')
  }
} # for (annot in rownames(annot_cols))

labs <- c()
labs.colnames <- c()

legend.text <- c()
legend.text.col <- c()
legend.text.font <- c()
legend.fill <- c()


for (annot in names(cols)) {
  if (annot %in% names(numeric_cols)) {
    pct <- (numeric_cols[[annot]][['high']]-numeric_cols[[annot]][['low']])/10
    vals <- sapply(classes[[annot]], function(x){min((x-numeric_cols[[annot]][['low']])%/%pct+1, 10)})
    vals <- cols[[annot]][vals]
    vals[is.na(vals)] <- '#000000'
    labs <- cbind(labs, vals)
    labs.colnames <- c(labs.colnames, annot)
    
    legend.text <- c(legend.text, 
                     annot,
                     "NA", 
                     c("low", "medium", "high"))
    legend.fill <- c(legend.fill, 
                     '#FFFFFF',
                     '#000000', 
                     cols[[annot]][c(1,5, 10)])
    legend.text.col <- c(legend.text.col,
                         '#000000',
                         '#2B2D2F',
                         rep('#2B2D2F', 3)
    )
    legend.text.font <- c(legend.text.font,
                          4,
                          1,
                          rep(1, 3)
    )
  } else {
    labs <- cbind(labs, cols[[annot]][classes[[annot]]])
    labs.colnames <- c(labs.colnames, 
                       annot)
    
    legend.text <- c(legend.text, 
                     annot,
                     levels(classes[[annot]]))
    legend.text.col <- c(legend.text.col, 
                         "#000000",
                         rep('#2B2D2F', length(levels(classes[[annot]]))))
    legend.text.font <- c(legend.text.font, 
                          4,
                          rep(1, length(levels(classes[[annot]]))))
    
    if (annot %in% c('plate', 'washes', 'sequencer')) {
      # We specifically use the levels since we've manually set them in the order we care
      legend.fill <- c(legend.fill, 
                       '#FFFFFF',
                       cols[[annot]][factor(levels(classes[[annot]]), 
                                            levels=fixed_levels[[annot]])])
      
    } else if (annot == 'EHK') {
      # Has to be numeric for this to be correct)
      legend.fill <- c(legend.fill, 
                       '#FFFFFF',
                       cols[[annot]][as.numeric(levels(classes[[annot]]))])
    } else {
      legend.fill <- c(legend.fill, 
                       '#FFFFFF',
                       # Need to do the levels hack so the levels don't get resorted
                       cols[[annot]][factor(levels(classes[[annot]]), levels=as.vector(levels(classes[[annot]])))])
    }
  } # if (annot %in% names(numeric_cols)) / else
} # for (annot in names(cols)) 

colnames(labs) <- labs.colnames
rownames(labs) <- rownames(classes)
cat ("Setup all color palettes.\n")

hclust.ward = function(d) hclust(d,method="ward.D2")

for (nm in names(logcpm_matrices)) {
  prefix_nm=file.path(args$OUT_FOLDER, args$GROUP_COL, nm)
  prefix_pca=file.path(prefix_nm, 'pca')
  dir.create(prefix_pca, showWarnings = FALSE, recursive = TRUE)
  
  nmlcpm <- logcpm_matrices[[nm]]
  
  cat(paste0('Plotting raw pca with ', nm, ' cpms\n'))
  pcs <- prcomp(t(nmlcpm))
  ap <- list()
  for (annot_name in names(cols)) {
    if (annot_name == args$GROUP_COL){
      ap[[annot_name]] <- autoplot(pcs, data=classes, shape=annot_name, col=annot_name, size=5, alpha=0.75) + 
        geom_text_repel(aes(label=rownames(classes))) + 
        scale_color_manual(values=cols[[annot_name]]) +
        theme_bw()
    } else if (annot_name %in% names(numeric_cols)){
      ap[[annot_name]] <- autoplot(pcs, data=classes, shape=args$GROUP_COL, col=annot_name, size=5, alpha=0.75) + 
        scale_color_distiller(palette=numeric_cols[[annot_name]][["palette"]], 
                              limits = c(numeric_cols[[annot_name]][["low"]], 
                                         numeric_cols[[annot_name]][["high"]]),
                              na.value='#000000') +
        theme_bw()
    } else {
      ap[[annot_name]] <- autoplot(pcs, data=classes, shape=args$GROUP_COL, col=annot_name, size=5, alpha=0.75) + 
        scale_color_manual(values=cols[[annot_name]]) +
        theme_bw()
    }
  }
  
  png(filename=file.path(prefix_pca, paste('voom_', gsub('.{6}$', '', nm), '_pca.png',sep='')), width = 30, height =  ceiling(length(ap)/3) *10 , units = "in", res = 150)
  grid.arrange(grobs=ap, ncol = 3)
  dev.off()
  
  x <- as.data.frame(pcs$x)
  x[[args$GROUP_COL]] = classes[[args$GROUP_COL]]
  y <- melt(x, id.vars = args$GROUP_COL, variable.name = 'PC')
  
  png(filename=file.path(prefix_pca, paste('voom_', gsub('.{6}$', '', nm), '_pca_boxplots.png',sep='')), width = 30, height =  ceiling(length(ap)/3) *10 , units = "in", res = 150)
  print(ggplot(y, aes_string(x=args$GROUP_COL, y="value", fill=args$GROUP_COL)) + 
          geom_boxplot() + 
          stat_compare_means(method = "t.test") + 
          facet_wrap(~PC, ncol=5))
  dev.off()
  
  for (pc in c(1, 2)) {
    prefix_pc_x=file.path(prefix_pca, paste0('pc_', pc))
    dir.create(prefix_pc_x, showWarnings = FALSE)
    
    best_correlation <- data.frame(Sr=sort(apply(nmlcpm, 1, function(x) {cor(x, pcs$x[,pc], method="spearman")})))
    best_correlation$ensembl <- rownames(best_correlation)
    best_correlation$trend <- ifelse(best_correlation$Sr>0, 'POS', 'NEG')
    best_correlation <- best_correlation %>% group_by(trend) %>% top_n(20, wt=abs(Sr))
    best_correlation[['Sr_sqr']]  <-  round(best_correlation$Sr ** 2, 2)
    best_correlation[['Sr']]  <-  round(best_correlation$Sr, 2)
    best_correlation <- add_gene_symbol(best_correlation)
    best_correlation <- best_correlation[, c('hugo', 'ensembl', 'Sr', 'Sr_sqr', 'trend')]
    write.table(best_correlation, 
                file=file.path(prefix_pc_x, 
                               paste0('voom_', gsub('.{6}$', '', nm), '_PC', pc, '_correlating_genes.tsv')),
                sep='\t', quote=F, row.names = F, col.names = T)
    
    top_rotations <- data.frame(weight=round(pcs$rotation[, pc][order(-abs(pcs$rotation[, pc]))][c(1:200)], 4))
    top_rotations$ensembl <- rownames(top_rotations)
    top_rotations <- add_gene_symbol(top_rotations)
    top_rotations$trend <- ifelse(top_rotations$weight>0, 'POS', 'NEG')
    top_rotations <- top_rotations[, c("ensembl", "hugo", "weight", "trend")]
    write.table(top_rotations, 
                file=file.path(prefix_pc_x, 
                               paste0('voom_', gsub('.{6}$', '', nm), '_PC', pc, '_top200_weights.tsv')),
                row.names = F, 
                col.names=T, 
                sep="\t", 
                quote=F)
    
    if (args$RUN_STRING_PCA){
      cat(paste0('Running STRING on pc ', pc, '\n'))
      try(run_string(file.path(prefix_pc_x, paste0('voom_', gsub('.{6}$', '', nm), '_PC', pc)),
                     top_table=top_rotations, 
                     string_db=string_db, 
                     usesig=FALSE,
                     logFcColStr="weight",
                     trend_col="weightage",
                     trend_classes=c('POS', 'NEG')))
    }
    if (args$RUN_DAVID_PCA){
      cat(paste0('Running DAVID on pc ', pc, '\n'))
      try(run_david(file.path(prefix_pc_x, paste0('voom_', gsub('.{6}$', '', nm), '_PC', pc)),
                    top_table=top_rotations, 
                    david_db=david_db, 
                    genelist_prefix = paste0(nm, '_', pc),
                    logFcColStr="weight",
                    trend_col="weightage",
                    trend_classes=c('POS', 'NEG')))
    }
  } 
  
  for (mvm in args$MOST_VARIABLE_METHOD) {
    if (mvm == 'var') {
      ranked_variance = rownames(fit2$coefficients)[order(-fit2$sigma)][1:max(args$MOST_VARIABLE_CUTOFF)]  
    } else if (mvm == 'IQR') {
      cat('Calculating IQRs\n')
      iqr = apply(nmlcpm, 1, IQR)
      ranked_variance = names(sort(iqr, decreasing = T)[1:max(args$MOST_VARIABLE_CUTOFF)])
    } else if (mvm == 'CQV') {
      cat('Calculating CQVs\n')
      cqv = apply(nmlcpm, 1, function(x){
        q1 = quantile(x, 0.25)
        q3 = quantile(x, 0.75)
        100 * (q3-q1) / (q3+q1)
        })
      ranked_variance = names(sort(cqv, decreasing = T)[1:max(args$MOST_VARIABLE_CUTOFF)])
    } else if (mvm == 'IQV') {
      iqv = apply(nmlcpm, 1, function(x){
        q1 = quantile(x, 0.25)
        q3 = quantile(x, 0.75)
        var(x[x >= q1 & x <= q3])
      })
      ranked_variance = names(sort(iqv, decreasing = T)[1:max(args$MOST_VARIABLE_CUTOFF)])
    } else {
      cat(paste0('ERROR: Unimplemented MOST_VARIABLE_METHOD encountered: ', mvm, '\n'))
      next
    }
    
    prefix_ubh=file.path(prefix_nm, 'unbiased_heatmaps', mvm)
    dir.create(prefix_ubh, showWarnings = FALSE, recursive = TRUE)
    
    for (mvc in args$MOST_VARIABLE_CUTOFF) {
      cat(paste0('Plotting unbiased heatmap with ', nm, ' cpms across ', mvc, ' variable genes using the ', mvm, ' method.\n'))
      top_most_variable_genes = ranked_variance[1:mvc]
      png(filename=file.path(prefix_ubh, paste('voom_', 
                                               gsub('.{6}$', '', nm), 
                                               '_unbiased_top', 
                                               mvc,
                                               '_',
                                               mvm,
                                               '_heatmap.png',
                                               sep='')), 
          width = 20, height = 15, units = "in", res = 150)
      temp_rownames <- ensembl_to_symbol[ensembl_to_symbol$ensembl %in% top_most_variable_genes, ]
      rownames(temp_rownames) <- temp_rownames$ensembl
      temp_rownames <- temp_rownames[top_most_variable_genes, ]
      temp <- nmlcpm[top_most_variable_genes,]
      rownames(temp) <- temp_rownames$hugo
      x <- coolmap.2(temp, 
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
             ncol=1, 
             text.col = legend.text.col,
             text.font = legend.text.font)
      dev.off()
      png(filename=file.path(prefix_ubh, paste('voom_', 
                                               gsub('.{6}$', '', nm), 
                                               '_unbiased_top', 
                                               mvc,
                                               '_',
                                               mvm,
                                               '_rowdendro.png',
                                               sep='')),
          width = 25, height = 5, units = "in", res = 150)
      plot(x$rowDendrogram)
      dev.off()
      writeLines(rownames(temp)[order.dendrogram(x$rowDendrogram)], 
                 con=file.path(prefix_ubh, paste('voom_', 
                                                 gsub('.{6}$', '', nm), 
                                                 '_unbiased_top', 
                                                 mvc,
                                                 '_',
                                                 mvm,
                                                 '_rowdendro.list',
                                                 sep='')))
    } # for (mvc in args$MOST_VARIABLE_CUTOFF)
  } # for (mvm in args$MOST_VARIABLE_METHOD)
  
  for (gs in args$GENESETS) {
    gs_name = strsplit(basename(gs), '[.]')[[1]][1]
    
    gs_genes = readLines(gs)
    gs_genes = gs_genes[gs_genes %in% rownames(nmlcpm)]
    
    prefix_gs=file.path(prefix_nm, 'unbiased_heatmaps', gs_name)
    dir.create(prefix_gs, showWarnings = FALSE, recursive = TRUE)
    
    cat(paste0('Plotting unbiased heatmap using geneset ', gs_name, '.\n'))
    png(filename=file.path(prefix_gs, paste('voom_', 
                                            gsub('.{6}$', '', nm), 
                                            '_unbiased_',
                                            gs_name,
                                            '_heatmap.png',
                                            sep='')), 
        width = 20, height = 15, units = "in", res = 150)
    temp_rownames <- ensembl_to_symbol[ensembl_to_symbol$ensembl %in% gs_genes, ]
    rownames(temp_rownames) <- temp_rownames$ensembl
    temp_rownames <- temp_rownames[gs_genes, ]
    temp <- nmlcpm[gs_genes,]
    rownames(temp) <- temp_rownames$hugo
    x <- coolmap.2(temp, 
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
           ncol=1, 
           text.col = legend.text.col,
           text.font = legend.text.font)
    dev.off()
    png(filename=file.path(prefix_gs, paste('voom_', 
                                            gsub('.{6}$', '', nm), 
                                            '_unbiased_',
                                            gs_name,
                                            '_rowdendro.png',
                                            sep='')),
        width = 25, height = 5, units = "in", res = 150)
    plot(x$rowDendrogram)
    dev.off()
    writeLines(rownames(temp)[order.dendrogram(x$rowDendrogram)], 
               con=file.path(prefix_gs, paste('voom_', 
                                              gsub('.{6}$', '', nm), 
                                              '_unbiased_',
                                              gs_name,
                                              '_rowdendro.list',
                                              sep='')))
  } # for (gs in args$GENESETS)
} # for (nm in names(logcpm_matrices))

top_tables = list()

n=50
for (cname in cnames) {
  prefix_dge=file.path(args$OUT_FOLDER, args$GROUP_COL, 'DGE', cname)
  dir.create(prefix_dge, showWarnings = FALSE, recursive = TRUE)
  
  cat(paste0('Processing contrast ', cname, '\n'))
  top_tables[[cname]] = topTable(fit2, coef=cname, sort="p", number=Inf)
  top_tables[[cname]]$col <- 0
  top_tables[[cname]]$ensembl <- rownames(top_tables[[cname]])
  top_tables[[cname]] <- add_gene_symbol(top_tables[[cname]])
  top_tables[[cname]] <- top_tables[[cname]][order(top_tables[[cname]]$logFC),]
  top_n <- c(rownames(head(top_tables[[cname]][top_tables[[cname]]$logFC<0 & top_tables[[cname]]$adj.P.Val <= args$FDR_THRESHOLD,], n)), 
             rownames(tail(top_tables[[cname]][top_tables[[cname]]$logFC>0 & top_tables[[cname]]$adj.P.Val <= args$FDR_THRESHOLD,], n)))
  top_tables[[cname]][rownames(top_tables[[cname]][top_tables[[cname]]$P.Value <= 0.01,]), "col"] = 1
  top_tables[[cname]][rownames(top_tables[[cname]][top_tables[[cname]]$adj.P.Val<=args$FDR_THRESHOLD,]), "col"] = 2
  top_tables[[cname]]$col <- factor(top_tables[[cname]]$col, levels=c(0,1,2))
  
  write.table(top_tables[[cname]], 
              file=file.path(prefix_dge, paste('voom_', cname, '.tsv',sep='')), 
              sep='\t', 
              row.names=F, 
              col.names=T, 
              quote=F)
  p <- ggplot(top_tables[[cname]], aes(logFC, -log10(P.Value))) + 
    geom_point(aes(col=col)) + 
    scale_color_manual(name='Significance',
                       values=c("grey", "black", "red"),
                       labels=c('p > 0.01', 'p <= 0.01', paste0('fdr <= ', args$FDR_THRESHOLD)),
                       drop=F) + 
    theme(legend.position = "none") + 
    theme_bw() +
    geom_text_repel(data=top_tables[[cname]][top_n,], 
                    aes(label=hugo), 
                    size=2, 
                    color='black')
  png(filename=file.path(prefix_dge, paste('voom_', cname, '.png',sep='')), width = 10, height = 10, units = "in", res = 150)
  print(p)
  dev.off()
  
  top_tables[[cname]] = top_tables[[cname]][top_tables[[cname]]$adj.P.Val<args$FDR_THRESHOLD, ]
  if(dim(top_tables[[cname]])[1] == 0) {
    cat(paste0('WARNING: No genes passed the fdr threshold of ', args$FDR_THRESHOLD, '.\n'))
    next
  }
  
  if (args$RUN_STRING_DE){
    try(run_string(file.path(prefix_dge, 
                             paste0(cname, '_DE')),
                   top_table = top_tables[[cname]], 
                   string_db = string_db, 
                   usesig=FALSE,
                   logFcColStr="logFC",
                   trend_col="regulation",
                   trend_classes=c('UP', 'DOWN')))
  }
  if (args$RUN_DAVID_DE){
    try(run_david(file.path(prefix_dge, 
                            paste0(cname, '_DE')),
                  top_table = top_tables[[cname]], 
                  david_db = david_db, 
                  genelist_prefix = cname,
                  logFcColStr="logFC",
                  trend_col="regulation",
                  trend_classes=c('UP', 'DOWN')))
  }
  
  for (nm in names(logcpm_matrices)) {
    nmlcpm <- logcpm_matrices[[nm]]
    cat(paste0('Plotting raw heatmap with ', nm, ' cpms\n'))
    
    png(filename=file.path(prefix_dge, paste('voom_', cname, '_', gsub('.{6}$', '', nm), '_heatmap.png',sep='')), width = 20, height = 15, units = "in", res = 150)
    
    temp <- nmlcpm[rownames(top_tables[[cname]]),]
    rownames(temp) <- top_tables[[cname]]$hugo
    x <- coolmap.2(temp, 
                   ColSideColors=labs, 
                   margins=c(12, 25),
                   cbar_extreme=2)  
    legend("topright",
           legend=legend.text,
           fill=legend.fill,
           border=FALSE,
           bty="n",
           y.intersp = 0.7,
           cex=0.7,
           ncol=1, 
           text.col = legend.text.col,
           text.font = legend.text.font)
    dev.off()
    png(filename=file.path(prefix_dge, paste('voom_', cname, '_', gsub('.{6}$', '', nm), '_rowdendro.png',sep='')), width = 25, height = 5, units = "in", res = 150)
    plot(x$rowDendrogram)
    dev.off()
    writeLines(rownames(temp)[order.dendrogram(x$rowDendrogram)], 
               con=file.path(prefix_dge, paste('voom_', cname, '_', gsub('.{6}$', '', nm), '_rowdendro.list',sep='')))
  } # for (nm in names(logcpm_matrices))
} # for (cname in cnames)

n=100
if (length(top_tables) > 1) {
  prefix_ava=file.path(args$OUT_FOLDER, args$GROUP_COL, 'DGE', 'all_vs_all')
  dir.create(prefix_ava, showWarnings = FALSE, recursive = TRUE)
  
  top_tables_rest <- top_tables[grepl('_vs_rest', names(top_tables))]
  top_n_per <- lapply(
    top_tables_rest,
    function(x){
      (x %>% 
         filter(logFC>0) %>%
         top_n(n, wt=abs(logFC)))$ensembl
    }
  )
  top_n_per <- add_gene_symbol(data.frame(ensembl=unname(unlist(top_n_per)), stringsAsFactors = F))
  
  for (nm in names(logcpm_matrices)) {
    nmlcpm <- logcpm_matrices[[nm]]
    cat(paste0('Plotting raw heatmap with ', nm, ' cpms\n'))
    
    png(filename=file.path(prefix_ava, paste('voom_all_vs_all_', gsub('.{6}$', '', nm), '_heatmap.png',sep='')), width = 20, height = 15, units = "in", res = 150)
    
    temp <- nmlcpm[top_n_per$ensembl,]
    rownames(temp) <- top_n_per$hugo
    x <- coolmap.2(temp, 
                   ColSideColors=labs, 
                   margins=c(12, 25),
                   cbar_extreme=2)  
    legend("topright",
           legend=legend.text,
           fill=legend.fill,
           border=FALSE,
           bty="n",
           y.intersp = 0.7,
           cex=0.7,
           ncol=1, 
           text.col = legend.text.col,
           text.font = legend.text.font)
    dev.off()
    png(filename=file.path(prefix_ava, paste('voom_all_vs_all_', gsub('.{6}$', '', nm), '_rowdendro.png',sep='')), width = 25, height = 5, units = "in", res = 150)
    plot(x$rowDendrogram)
    dev.off()
    writeLines(rownames(temp)[order.dendrogram(x$rowDendrogram)], 
               con=file.path(prefix_ava, paste('voom_all_vs_all_', gsub('.{6}$', '', nm), '_rowdendro.list',sep='')))
    
  } # for (nm in names(logcpm_matrices))
}
