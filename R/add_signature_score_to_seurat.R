get_genes <- function(genes, sobj) {
  ###
  # DESCRIPTION
  # Get a dataframe of gene values per cell from the input seurat object
  #
  # INPUTS
  # genes: A character vector containing gene names
  # sobj: The Seurat object
  # OUTPUT
  # A dataframe of input genes as rows and all cells as columns
  ###
  gene_idx <- sapply(genes, function(x) { match(x, rownames(sobj@data)) })
  if (sum(is.na(gene_idx)) > 0) {
    print("The following genes are not in the gene list.")
    print(names(gene_idx[is.na(gene_idx)]))
    return(1)
  }
  genes <- as.data.frame(as.matrix(sobj@data[genes,]))
}


get_genes3 <- function(genes, sobj) {
  ###
  # DESCRIPTION
  # Get a dataframe of gene values per cell from the input seurat3 object
  #
  # INPUTS
  # genes: A character vector containing gene names
  # sobj: The Seurat object
  # OUTPUT
  # A dataframe of input genes as rows and all cells as columns
  ###
  gene_idx <- sapply(genes, function(x) { match(x, rownames(sobj)) })
  if (sum(is.na(gene_idx)) > 0) {
    print("The following genes are not in the gene list.")
    print(names(gene_idx[is.na(gene_idx)]))
    return(1)
  }
  genes <- as.data.frame(as.matrix(sobj@assays[[sobj@active.assay]]@data[genes,]))
}


saturate <- function(vec, sat=0, binary=FALSE){
  ###
  # DESCRIPTION
  # A Function to convert a vector of scores into a saturated vectore of scores. A saturated vector is one where all values below the
  # provided "saturation" (percentile of data) are set to 0. If the binary flag is specified, all values greater than or equal to the
  # saturation will be set to 1.
  #
  # INPUTS
  # vec: A numeric vector of scores
  # sat: A value to saturate the vector at (float (0.0-1.0) or percent (1.0-100.0))
  # binary: A flag to indicate if we should make the output vector a binary one.
  #
  # OUTPUT
  # A vector of saturated scores
  ###
  sat = if (sat > 1.0) sat/100 else sat
  z <- quantile(vec, sat)
  for (i in 1:length(vec)){
    if (vec[i] < z) {
      vec[i] = 0
    } else if(binary) {
      vec[i] = 1
    }
  }
  vec
}




# Example usage
if (FALSE){
  # A list of genes in the signature
  genes = c("SPP1", "APOE", "SEPP1", "RNASE1", "APOC1", "C1QB", "C1QA", "C1QC", "CD14", "IFI27",
            "FOLR2", "HAMP", "GPNMB", "NUPR1", "FCGR3A", "PLTP", "CCL3", "A2M", "CTSB", "CCL2",
            "TREM2", "CTSL", "VSIG4", "MARCO", "CD63", "TMEM176B", "MARCKS", "CTSD", "LGALS1",
            "PSAP", "MS4A7")

  # Get the cell values from the Seurat object
  gene_df <- get_genes(genes, sobj)

  # Define the score vector. In this case, we call the score for a cell the mean value of all genes in the vector
  score <- colMeans(gene_df)

  # Code for if you want to define a custom score for the data
  # for example, median instead of mean
  # score <- sapply(gene_df, function(x) {median(x)})

  # Add score columns to the metadata column of the Seurat object
  sobj@meta.data$score50 <- saturate(vec=score, sat=0.5, binary=TRUE)    # Saturate at 50%
  sobj@meta.data$score75 <- saturate(vec=score, sat=0.75, binary=TRUE)   # Saturate at 75%
  sobj@meta.data$score90 <- saturate(vec=score, sat=0.90, binary=TRUE)   # Saturate at 90%

  # Make a featureplot
  FeaturePlot(object = sobj,
              features.plot = c("score50", "score75", "score90"),
              cols.use=c("light grey", "red"))
}