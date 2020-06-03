suppressPackageStartupMessages({
  library(rlang)
  if (!'package:Seurat' %in% search()) library(Seurat)
})


seurat_to_loupe <- function(sobj, 
                            file_prefix, 
                            reductions=c('umap'),
                            metadata_cols=NULL,
                            barcode_suffix=NULL){
  ## A function to extract data from a seurat object and bump it into files that can be loaded into
  ## the 10x loupe cell browser. The files that will be written include
  ##
  ##     <FILE_PREFIX>_metadata.csv 
  ##         A comma separated table containing metadata from the seurat object
  ## 
  ##     <FILE_PREFIX>_<REDUCTION>_coords.csv 
  ##         A comma separated table with the first and second components of a dimensional 
  ##         reduction. One such file will be generated for each specified reduction.
  ##
  ## Parameters:
  ##     sobj : The Seurat object (REQUIRED)
  ##     file_prefix: The prefix to all output files. (REQUIRED)
  ##     reductions: The reductions to write to file. The reduction must be in the reductions slot 
  ##                 of sobj. (DEFAULT='umap')
  ##     metadata_cols: Specific columns in the metadata to write. (DEFAULT=all columns)
  ##
  ## Return value
  ##     None
  ##
  if (is.null(sobj)){
    stop('sobj cannot be null')
  }
  if (is.null(file_prefix)){
    stop('file_prefix cannot be null')
  }
  
  metadata_cols = metadata_cols %||% colnames(sobj@meta.data)
  
  drop_metadata_cols = metadata_cols[!metadata_cols%in%colnames(sobj@meta.data)]
  for (dmc in drop_metadata_cols){
    cat(paste0('Dropping missing metadata column : ', dmc, '\n'))
  }
  metadata_cols = metadata_cols[metadata_cols%in%colnames(sobj@meta.data)]
  
  if (length(metadata_cols) == 0){
    cat('WARNING: None of the provided metadata columns are in the object. Not generating a metadata csv')
  } else {
    metadata = sobj@meta.data[, metadata_cols]
    if ('barcode' %in% metadata_cols) {
      cat('Using barcode column in the metadata table. \n')
    } else {
      cat('Using cell names as the Barcodes. \n')
        metadata$barcode <- rownames(metadata)
    }
    if (!is.null(barcode_suffix)) {
      cat('Appending `', barcode_suffix, '` to all barcodes. \n')
      metadata$barcode <- paste0(metadata$barcode, barcode_suffix)
    }
    metadata = metadata[, c('barcode', metadata_cols[metadata_cols!='barcode'])]
    write.table(metadata,
                file=paste0(file_prefix, '_metadata.csv'),
                sep=',',
                row.names=F,
                col.names=T,
                quote=F)
  }

  for (redn in reductions){
    if (!redn %in% names(sobj@reductions)){
      cat(paste0('WARNING: Could not find reduction `', redn, '` in sobj. Skipping.\n'))
    } else {
      cell_embeddings <- as.data.frame(sobj@reductions[[redn]]@cell.embeddings[, c(1, 2)])
      cc <- colnames(cell_embeddings)
      cell_embeddings$barcode <- metadata$barcode
      cell_embeddings <- cell_embeddings[, c('barcode', cc)]
      write.table(cell_embeddings,
                  file=paste0(file_prefix, '_', redn,'_coords.csv'),
                  sep=',',
                  row.names=F,
                  col.names=T,
                  quote=F)
    }
  }
}
