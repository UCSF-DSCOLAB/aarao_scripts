suppressPackageStartupMessages({
  if (!'package:assertthat' %in% search()) library(assertthat)
  if (!'package:cowplot' %in% search()) library(cowplot)
  if (!'package:DoubletFinder' %in% search()) library(DoubletFinder)
  if (!'package:Seurat' %in% search()) library(Seurat)
})

runDoubletFinder <- function(sobj, doublet_rate, out_prefix, num_cores=1, annot="SCT_snn_res.0.8", PCs=1:30, sct=TRUE, pN=0.25) {
    sweep.list <- paramSweep_v3(sobj, PCs=PCs, sct=sct, num.cores=num_cores)
    sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
    png(paste0(out_prefix, '_doubletfinder_bcmvn.png'), height=500, width=500, units='px')
    bcmvn <- find.pK(sweep.stats)
    dev.off()
    write.table(bcmvn,
                file=paste0(out_prefix, '_doubletfinder_bcmvn.tsv'),
                sep='\t',
                row.names=T,
                col.names=T,
                quote=F)
    homotypic.prop <- modelHomotypic(as.vector(sobj@meta.data[[annot]]))
    nExp_poi <- round(doublet_rate*(dim(sobj)[1]))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    # Kinda invariant.... use default
    pN <- pN
    # Maxima of the BCmvn distribution
    pK <- as.numeric(as.vector(bcmvn[which.max(bcmvn$BCmetric), 'pK']))
    sobj <- doubletFinder_v3(sobj,
                             PCs = 1:30,
                             pN = pN,
                             pK = pK,
                             nExp = nExp_poi.adj,
                             reuse.pANN = FALSE,
                             sct = TRUE)
    sobj$DF.DROPLET.TYPE <- sobj@meta.data[[paste0('DF.classifications_', pN, '_', pK, '_', nExp_poi.adj)]]
    sobj
}