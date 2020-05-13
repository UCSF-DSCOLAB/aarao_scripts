library(STRINGdb)

setup_string <- function(string_folder, version="10", species=9606, score_threshold=0){
  string_db <- STRINGdb$new(version=version, 
                            species=species, 
                            score_threshold=score_threshold, 
                            input_directory=string_folder)
  temp <- string_db$get_annotations()
  return(string_db)
}

run_string <- function(out_prefix, top_table, string_db, gene_col='ensembl',
                       symbol_col='hugo', usesig=TRUE, sigval=0.05, 
                       logFcColStr="logFC", trend_col='regulation',
                       trend_classes=c('UPREG', 'DOWNREG'),
                       top_n=200, enrichment_n=20, num_attempts=10) {
  print("Beginning a STRING analysis:")
  cat("Attempting to connect to STRING....")
  attempts = 0
  while (attempts < num_attempts) {
    top_table_mapped <- NULL
    try(top_table_mapped <- string_db$map(top_table, 
                                          gene_col, 
                                          removeUnmappedRows=TRUE),
        silent = TRUE
    )
    if (!is.null(top_table_mapped)){
      break
    } else {
      attempts <- attempts + 1
    }
  }
 
  if (is.null(top_table_mapped)){
    stop(paste0('ERROR: Could not reach STRING after ', num_attempts, ' attempts'))
  } else {
    cat("Done.\n")
  }
  
  dropped_genes <- top_table[!top_table[, gene_col, drop=T]%in%top_table_mapped[, gene_col, drop=T], gene_col]
  print("The following genes were dropped from the STRING analysis:")
  cat(dropped_genes, sep="\n")
  if (usesig){
    sigstring <- 'significant'
    top_table_mapped <- subset(top_table_mapped, P.Value<sigval)
  } else {
    sigstring <- 'all'
  }
  hits <- top_table_mapped$STRING_id[1:top_n]
  hits <- hits[!is.na(hits)]

  top_table_mapped <- string_db$add_diff_exp_color(top_table_mapped, 
                                                   logFcColStr=logFcColStr)
  # post payload information to the STRING server
  cat("Attempting to post a payload to STRING....")
  attempts = 0
  while (attempts < num_attempts) {
    payload_id <- NULL
    try(payload_id <- string_db$post_payload(top_table_mapped$STRING_id, 
                                             colors=top_table_mapped$color),
        silent = TRUE
    )
    if (!is.null(payload_id)){
      break
    } else {
      attempts <- attempts + 1
    }
  }
 
  if (is.null(payload_id)){
    stop(paste0('ERROR: Could not get a payload ID after ', num_attempts, ' attempts'))
  } else {
    cat("Done.\n")
  }
  
  # display a STRING network png with the "halo"
  png(filename=paste0(out_prefix, '_top', top_n, sigstring, '_STRINGdb.png',sep=''), 
                      width = 10, height = 10, units = "in", res = 150)
  string_db$plot_network(hits, payload_id=payload_id)
  dev.off()
  
  cat("Attempting to get GO and KEGG enrichments from STRING....")
  attempts = 0
  enrichments <- list()
  while (attempts < num_attempts) {
    payload_id <- NULL
    if (is.null(enrichments[['GO']])){
      try(enrichments[['GO']] <- string_db$get_enrichment(hits, category="Process", 
                                                          methodMT="fdr", iea=TRUE),
          silent = TRUE
      )
    } else if (is.null(enrichments[['KEGG']])){
      try(enrichments[['KEGG']] <- string_db$get_enrichment(hits, category="KEGG", 
                                                            methodMT="fdr", iea=TRUE),
          silent = TRUE
      )
    } else {
      break
    }
    attempts <- attempts + 0.5
  }
 
  if (is.null(enrichments[['GO']]) || is.null(enrichments[['KEGG']])){
    stop(paste0('ERROR: Could not get GO and KEGG enrichments after ', num_attempts, ' attempts'))
  } else {
    cat("Done.\n")
  }

  for (enrichment_name in names(enrichments)){
    enrichment <- enrichments[[enrichment_name]]
    write.table(enrichment, 
                file=paste0(out_prefix, '_', enrichment_name, 'enrichment', '.tsv'), 
                sep='\t', row.names=F, col.names=T, quote=F)

    enrichment_red <- head(enrichment, n=enrichment_n) 
    top_annotations <- string_db$annotations[string_db$annotations$term_id %in% enrichment_red$term_id,]
    top_annotations <- top_annotations[top_annotations$STRING_id %in% top_table_mapped$STRING_id,]
    top_annotations$gene <- sapply(top_annotations$STRING_id, function(x) {paste(top_table_mapped[top_table_mapped$STRING_id==x, gene_col], collapse=',')})
    if (!is.null(symbol_col) && symbol_col %in% colnames(top_table_mapped)){
      top_annotations$symbol <- sapply(top_annotations$STRING_id, function(x) {paste(top_table_mapped[top_table_mapped$STRING_id==x, symbol_col], collapse=',')})
    } else {
      top_annotations$symbol <- top_annotations$gene
    }
    top_annotations$description <- sapply(top_annotations$term_id, function(x) {enrichment_red[enrichment_red$term_id==x, 'term_description']})
  
    top_annotations[trend_col] <- sapply(top_annotations$gene, 
                                         function(x) {
                                              ifelse(!(x %in% top_table_mapped[[gene_col]]), 
                                                     'UNKNOWN', 
                                                     ifelse(top_table_mapped[top_table_mapped[[gene_col]]==x, logFcColStr] > 0, 
                                                            trend_classes[1], 
                                                            trend_classes[2]))
                                         })

    top_annotations <- top_annotations[, c('gene', 'symbol', trend_col, 'term_id', 'description')]
    top_annotations$term_id <- factor(top_annotations$term_id, levels=enrichment_red$term_id)
    top_annotations <- top_annotations[order(top_annotations$term_id, top_annotations$weightage),]
    write.table(top_annotations, 
                file=paste0(out_prefix,'_', enrichment_name, 'enrichment', '_top', enrichment_n, 'term_genes.tsv'), 
                sep='\t', row.names=F, col.names=T, quote=F)
  }
}