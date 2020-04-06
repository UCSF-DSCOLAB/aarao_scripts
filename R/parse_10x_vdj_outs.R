parse_tcr_clonotype <- function(tcr_outs_folder){
    tcr <- read.csv(paste(tcr_outs_folder,"filtered_contig_annotations.csv", sep=""))

    # Remove the -1 at the end of each barcode.
    # Subsets so only the first line of each barcode is kept,
    # as each entry for given barcode will have same clonotype.
    tcr$barcode <- gsub("-1", "", tcr$barcode)
    tcr <- tcr[!duplicated(tcr$barcode), ]

    # Only keep the barcode and clonotype columns.
    # We'll get additional clonotype info from the clonotype table.
    tcr <- tcr[,c("barcode", "raw_clonotype_id")]
    names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"

    # Clonotype-centric info.
    clono <- read.csv(paste(tcr_outs_folder,"clonotypes.csv", sep=""))

    # Slap the AA sequences onto our original table by clonotype_id.
    tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])

    # Reorder so barcodes are first column and set them as rownames.
    tcr <- tcr[, c(2,1,3)]
    rownames(tcr) <- tcr[,1]
    tcr[,1] <- NULL
    colnames(tcr) <- c("TCR.clonotype_id", "TCR.cdr3s_aa")
    tcr
}


# Parse the results of a 10x VDJ IG sequencing
parse_bcr_clonotype <- function(bcr_outs_folder){
    bcr <- read.csv(paste(bcr_outs_folder,"filtered_contig_annotations.csv", sep=""))

    # Remove the -1 at the end of each barcode.
    # Subsets so only the first line of each barcode is kept,
    # as each entry for given barcode will have same clonotype.
    bcr$barcode <- gsub("-1", "", bcr$barcode)
    bcr <- bcr[!duplicated(bcr$barcode), ]

    # Only keep the barcode and clonotype columns.
    # We'll get additional clonotype info from the clonotype table.
    bcr <- bcr[,c("barcode", "raw_clonotype_id")]
    names(bcr)[names(bcr) == "raw_clonotype_id"] <- "clonotype_id"

    # Clonotype-centric info.
    clono <- read.csv(paste(bcr_outs_folder,"clonotypes.csv", sep=""))

    # Slap the AA sequences onto our original table by clonotype_id.
    bcr <- merge(bcr, clono[, c("clonotype_id", "cdr3s_aa")])

    # Reorder so barcodes are first column and set them as rownames.
    bcr <- bcr[, c(2,1,3)]
    rownames(bcr) <- bcr[,1]
    bcr[,1] <- NULL
    colnames(bcr) <- c("BCR.clonotype_id", "BCR.cdr3s_aa")
    bcr
}