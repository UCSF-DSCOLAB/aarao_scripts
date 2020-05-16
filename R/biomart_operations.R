library(biomaRt)

setup_biomart  <- function(biomart, dataset, num_attempts=10){
  attempts = 0
  while (attempts < num_attempts) {
    human <- NULL
    try(mart <- useMart(biomart = biomart, 
                        dataset = dataset),
        silent = TRUE
    )
    if (!is.null(mart)){
      break
    } else {
      attempts <- attempts + 1
    }
  }
  mart
}


getBM_or_bust <- function(mart, attributes, filters, values, uniqueRows = T, curr_iter=1, max_attempts=10) {
  if (curr_iter > max_attempts){
    stop(paste0('Could not reach biomart after ', max_attempts, ' attempts'))
  }
  temp <- NULL
  try(temp <- getBM(mart = mart,
                    attributes = attributes, 
                    filters = filters,
                    values = values,
                    uniqueRows = uniqueRows))
  if (is.null(temp)) {
    getBM_or_bust(mart = mart,
                  attributes = attributes, 
                  filters = filters,
                  values = values,
                  uniqueRows = uniqueRows,
                  curr_iter=curr_iter+1,
                  max_attempts=max_attempts)
  } else {
    temp
  }
}

iterative_getBM <- function(mart, attributes, filters, values, uniqueRows = T, max_records_per_iter=500, max_attempts_per_iter=10) {
  # Biomart gets annoyed if you submit too many genes so chunking it is better
  out <- NULL
  values = split(values, ceiling(seq_along(values)/max_records_per_iter))
  cat(paste0('Created ', length(values), ' chunks.\n'))
  for (i in 1:length(values)){
    cat(paste0('Processing chunk ', i, '\n'))
    temp_values <- values[[i]]
    temp <- getBM_or_bust(mart = human,
                          attributes = attributes, 
                          filters = filters,
                          values = temp_values,
                          uniqueRows = uniqueRows,
                          max_attempts = max_attempts_per_iter)
    if (is.null(out)){
      out = temp
    } else {
      out <- rbind(out, temp)
    }
  }
  out
}

ensg_to_hugo <- function(mart, 
                         input_df, 
                         ensembl_col='ensembl', 
                         hugo_col='hugo', 
                         fill=TRUE,
                         handle_duplicates='first',
                         max_records_per_iter=500,
                         max_attempts_per_iter=10){
  
  if (!handle_duplicates %in% c('first', 'merge')){
    stop("handle_duplicates must be one of 'first' or 'merge'")
  }
  if (!fill %in% c(TRUE, FALSE)){
    stop("fill must be one of TRUE or FALSE")
  }

  ensembl <- iterative_getBM(mart = human,
                             attributes = c("hgnc_symbol","ensembl_gene_id"), 
                             filters = "ensembl_gene_id", 
                             values = input_df[, ensembl_col, drop=T], 
                             uniqueRows = T,
                             max_records_per_iter=max_records_per_iter,
                             max_attempts_per_iter =  max_attempts_per_iter)

  if (handle_duplicates == 'merge'){
    ensembl$hgnc_symbol <- sapply(ensembl$ensembl_gene_id, function(x){
      temp = unique(ensembl[ensembl$ensembl_gene_id==x, 'hgnc_symbol', drop=T])
      temp = temp[temp!=""]
      paste(temp, collapse='__')
    })
  }
  # else 'first'
  ensembl <- ensembl[!duplicated(ensembl$ensembl_gene_id),]
  
  original_order <- rownames(input_df)
  input_df$Row.names <- rownames(input_df)
  input_df <- transform(
    merge(input_df,
          ensembl,
          by.x=ensembl_col,
          by.y='ensembl_gene_id',
          all.x=TRUE),
    row.names=Row.names,
    Row.names=NULL)
  input_df[hugo_col] = input_df$hgnc_symbol
  input_df$hgnc_symbol = NULL
    
  input_df <- input_df[original_order,, drop=F]

  if (fill){
    input_df[is.na(input_df[[hugo_col]]), hugo_col] = input_df[is.na(input_df[[hugo_col]]), ensembl_col]
    input_df[input_df[[hugo_col]]=="", hugo_col] = input_df[input_df[[hugo_col]]=="", ensembl_col]
  }
  input_df
}
