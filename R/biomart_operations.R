library(biomaRt)

ensg_to_hugo <- function(input_df, 
                         ensembl_col='ensembl', 
                         hugo_col='hugo', 
                         fill=TRUE,
                         handle_duplicates='first',
                         num_attempts=10){
  
  if (!handle_duplicates %in% c('first', 'merge')){
    stop("handle_duplicates must be one of 'first' or 'merge'")
  }
  if (!fill %in% c(TRUE, FALSE)){
    stop("fill must be one of TRUE or FALSE")
  }
  # Load human ensembl attributes
  attempts = 0
  while (attempts < num_attempts) {
    human <- NULL
    try(human <- useEnsembl(biomart="ensembl", 
                            dataset = "hsapiens_gene_ensembl",
                            mirror = "useast"),
        silent = TRUE
    )
    if (!is.null(human)){
      break
    } else {
      attempts <- attempts + 1
    }
  }
 
   if (is.null(human)){
    stop(paste0('ERROR: Could not load ensembl biomart after ', num_attempts, ' attempts'))
  }
  
  attempts = 0
  while (attempts < num_attempts) {
    ensembl <- NULL
    try(ensembl <- getBM(mart = human,
                         attributes = c("hgnc_symbol","ensembl_gene_id"), 
                         filters = "ensembl_gene_id", 
                         values = input_df[, ensembl_col, drop=T], 
                         uniqueRows = T),
        silent = TRUE
    )
    if (!is.null(ensembl)){
      break
    } else {
      attempts <- attempts + 1
    }
  }
  if (is.null(ensembl)){
    stop(paste0('ERROR: Could not convert symbols after ', num_attempts, ' attempts'))
  }
  
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
    input_df[hugo_col] <- sapply(rownames(input_df), function(x){
      if (is.na(input_df[x, hugo_col]) || input_df[x, hugo_col] == ""){
        input_df[x, ensembl_col]
      } else {
        input_df[x, hugo_col]
      }
    })
  }
  input_df
}
