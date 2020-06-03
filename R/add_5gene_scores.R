library(Seurat)
library(cowplot)
require(reshape2)

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

gene_sigs <- list(
    five_gene_sigs =  list(
        Mono_MACs_DCs=c('CD14', 'VCAN', 'LYZ', 'MS4A7', 'HLA-DRB1'),
        Bcell_plasmacell=c('MS4A1', 'CD79A', 'IGHA1', 'JCHAIN', 'HLA-DQA1'),
        Neuts=c('CSF3R', 'FCGR3B', 'CXCL8', 'G0S2', 'S100A8'),
        Platlet=c('PPBP', 'PF4', 'RGS10', 'CLEC1B', 'RGS18'),
        RBC=c('HBB', 'HBA2', 'HBA1', 'HBD', 'HBM'),
        Tcell_NKcell=c('CD3E', 'CD3D', 'IL7R', 'PRF1', 'GNLY')
    ),
    gtfive_gene_sigs = list(
      Neut_gt5=c("S100A8", "IFITM2", "S100A9", "S100A12", "S100A6", "IFITM3", "S100A11", "NEAT1", "TYROBP", "H3F3A", "FTH1", "ITM2B", "H3F3B", "S100A4", "IFITM1", "MNDA", "FTL", "GABARAP", "SRGN", "FCER1G", "MYL6", "ATP5E"),
      Platelet_gt5=c("PPBP", "GPX1", "PF4", "TUBB1", "TAGLN2", "RGS18", "SDPR", "SPARC", "HIST1H2AC", "OAZ1", "OST4", "TMSB4X", "CLU", "GNG11", "SH3BGRL3", "C6orf25", "CCL5", "HLA-E", "NRGN", "RGS10", "H3F3A", "ITM2B", "SERF2", "MYL6", "ACTB", "MYL12A"),
      Mono_Mac_gt5=c("S100A8", "S100A9", "LYZ", "S100A6", "CD74", "FCN1", "CST3", "FTL", "S100A4", "TMSB10", "TYROBP", "S100A10", "HLA-DRA", "VCAN", "IFITM3", "PSAP", "LGALS1", "VIM", "CTSS", "HLA-DRB1", "S100A11", "GRN", "TYMP", "AIF1", "IFI30", "FCER1G", "GAPDH", "GPX1", "TSPO", "SH3BGRL3", "ATP5G2", "ITGB2", "GABARAP", "SERF2"),
      T_NK_gt5=c("EEF1A1", "IL32", "CD3E", "TMSB4X", "HLA-B", "CD52", "PTPRCAP", "HLA-A", "TPT1", "TXNIP", "B2M", "EEF1B2", "HNRNPA1", "HCST", "HSPA8", "AC090498.1", "ZFP36L2", "PPIA", "HLA-C", "EEF2", "NPM1", "PFN1", "IFITM1", "DDX5", "IL2RG", "CALM1", "ACTG1"),
      RBC_gt5=c("HBB", "HBA2", "HBA1", "SLC25A37", "SLC25A39"),
      B_Plasma_gt5=c("EEF1A1", "FAU", "EEF1B2")
    )
)


for (i in Sys.glob("*/*_processed.RData")){
    j <- strsplit(i, '/')[[1]][1]
    load(i)
    sobj <- get(j)
    rm(list=j)
    

    for (gs in names(gene_sigs)){
        plots1 <- list()
        plots2 <- list()

        plots1[['clusters']] <- print(DimPlot(object=sobj,
                                             group.by='seurat_clusters') + NoLegend())
        plots2[['clusters']] <- print(DimPlot(object=sobj,
                                             group.by='seurat_clusters') + NoLegend())
        for (sig in names(gene_sigs[[gs]])){
            for (nb in c(50, 25, 12, 6)) {
                try(sobj <- AddModuleScore(sobj, 
                                           features=list(gene_sigs[[gs]][[sig]]), 
                                           name=sig, 
                                           nbin = nb))
                if (paste0(sig, '1') %in% colnames(sobj@meta.data)) {
                    break
                }
            }
            if (!paste0(sig, '1') %in% colnames(sobj@meta.data)) {
                stop()
            }
            
            sobj@meta.data[sig] <- sobj@meta.data[paste0(sig, '1')]
            sobj@meta.data[paste0(sig, '1')] <- NULL
            sobj@meta.data[[paste0(sig, '_75sat')]] <- saturate(sobj@meta.data[[sig]], sat=0.75, binary=T)
            plots1[[sig]] <- print(FeaturePlot(sobj,
                                              features=sig))
            plots2[[sig]] <- print(FeaturePlot(sobj,
                                              features=paste0(sig, '_75sat')))

        }

        png(filename=file.path(j, paste(j, gs, 'umap.png', sep='_')), width = 15, height = 15, units = "in", res = 300)
        print(plot_grid(plotlist = plots1, ncol=4))
        dev.off()

        png(filename=file.path(j, paste(j, gs, '_sat_umap.png', sep='_')), width = 15, height = 15, units = "in", res = 300)
        print(plot_grid(plotlist = plots2, ncol=4))
        dev.off()

        y <- reshape2::melt(sobj@meta.data[, names(gene_sigs[[gs]])])

        png(file.path(j, paste(j, gs, 'profiles.png', sep='_')), width = 15, height = 15, units = "in", res = 300)
        print(ggplot(y) + 
                  geom_density(aes(x=value, col=variable)) + 
                  facet_wrap(~variable))
        dev.off()
    }
}