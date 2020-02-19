#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(beeswarm)
  library(dplyr)
  library(edgeR)
  library(gplots)
  library(ggplot2)
  library(ggrepel)
  library(limma)
  library(methods)
  library(RColorBrewer)
})

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 3 || length(args) > 4) {
  stop("Need 3 or 4 arguments <COUNTS.TSV>, <CLASSES.TSV>, <GENES_OF_INTEREST.LIST> [<OUT_FOLDER>] in that order", call.=FALSE)
} else if (length(args) == 3) {
  args[4] <- "."
} else {
    args[4] <- gsub('/$', '', args[3])
}

print ("Reading counts.....")
counts <- read.table(args[1], sep='\t', header=1, row.names=1)
print ("counts read")

print ("Reading classes.....")
groups <- read.table(args[2], sep='\t', header=1, row.names=1, colClasses=c('factor'))
print ("classes read")

print ("Reading genes.....")
genes_of_interest <- read.table(args[3], sep='\t', header=FALSE, row.names=NULL, colClasses=c('character'))
colnames(genes_of_interest) <- "name"
print ("classes read")

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

# See https://support.bioconductor.org/p/76837/
# Suggests this method to remove batch effects before heatmap generation is sound. 

y = DGEList(counts=counts, group=groups$group)
print ("made DGElist...")
y = y[ rowSums(y$counts) > 60, keep.lib.size=FALSE ]
y <- calcNormFactors(y)
print ("Normalization done")
design <- model.matrix(~0+group, data=groups) #limma
colnames(design) <- gsub("group", "", colnames(design))
print ("Made Design Model")

logCPM <- cpm(y, log=TRUE, prior.count=0.1)
logCPM_no_batch <- removeBatchEffect(logCPM, batch=groups$plate, design=design)

gene_mats = list(
    uncorrected =  logCPM[genes_of_interest$name,],
    corrected =  logCPM_no_batch[genes_of_interest$name,]
  )


if (length(genes_of_interest$name) < 4) {  
  for (cc in names(gene_mats)) {
    gene_mats[[cc]] =  cbind(t(gene_mats[[cc]]), groups[, 'group', drop=FALSE])
    png(filename=paste(args[4], paste(cc, '_expressions.png',sep=''), sep='/'), width = 10, height = 10, units = "in", res = 300)
    par(mfcol=c(length(genes_of_interest$name), 1))
    for (gene in genes_of_interest$name){
      beeswarm(as.formula(paste0(gene, '~group')),data=gene_mats[[cc]], pch = 16, col = rainbow(8),
               main=paste0("Swarm for ", gene, "(", cc, ")"),  xlab="Group", 
               ylab="Mean Centered TPM")
    }
    dev.off()
  }
} else {
  colors = sapply(groups$group, function(x) {rainbow(8)[x]})
  for (cc in names(gene_mats)) {
    png(filename=paste(args[4], paste(cc, '_expressions.png',sep=''), sep='/'), width = 10, height = 10, units = "in", res = 300)
    par(mar=c(7,4,4,5)+0.1)
    heatmap.2(gene_mats[[cc]], ColSideColors=colors, scale="row", col=brewer.pal(11,"RdBu"), Rowv=NA, trace="none", margins=c(12,5), cexRow=1)
    dev.off()
  }
}