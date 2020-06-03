suppressPackageStartupMessages({
  library(assertthat)
  library(Seurat)
  library(dplyr)    # inline modification of matrices
  library(cowplot)  # Pretty plots
  library(ggplot2)  # Pretty plots
  library(scales)   # to access break formatting functions
  library(yaml)     # For reading cutoffs
})

# Set a seed so results are reproducible
seed = 21212

defaultArgs = list(
  SAMPLE_YML=list(
    value=NA,
    class=as.character,
    description="Samples to pull from the directory."
  ),
  OUT_FOLDER=list(
    value=getwd(),
    class=as.character,
    description="Where to write results."
  ),
  RSCRIPTS_DIR=list(
    value=paste(Sys.getenv("SCRIPTS"), "R", sep="/"),
    class=as.character,
    description="Where to find aux scripts"
  ),
  GENESET_DIR=list(
    value=paste(Sys.getenv("KLAB"), "ipi/data/refs/10x/genesets", sep="/"),
    class=as.character,
    description="Where to find genesets for ribo and cell cycling."
  ),
  SPECIES=list(
    value="human",
    class=as.character,
    description="Species to use for cell cycling and ribo. Must be `human` or `mouse` (case-sensitive)."
  ),
  RERUN_STAGE=list(
    value=FALSE,
    class=as.logical,
    description="Should we rerun the last run stage?"
  ),
  AB_ASSAY_NAME=list(
    value="IDX",
    class=as.character,
    description="What should the antibody capture assay (if any) be called?"
  )
)

print_help <- function(defaultArgs) {
  helptext = "
This is an R Script to process a 10X samples from raw feature/barcode matrices.

USAGE:
  
  Rscript process_10x_with_Seurat.R SAMPLE_YML=/path/to/samples.yml \
                                    OPTIONAL_ARG1=OPTIONAL_VAL1 \
                                    OPTIONAL_ARG2=OPTIONAL_VAL2 \
                                              ....
                                    OPTIONAL_ARGn=OPTIONAL_VALn \

Required Arguments:
SAMPLE_YML : "
  cat(helptext)
  cat(paste0(defaultArgs[['SAMPLE_YML']][['description']], 
      "\n\n",
      "Optional Arguments:\n"))
  for (arg in names(defaultArgs)){
    if (arg == 'SAMPLE_YML'){
      next
    }
    cat(paste0(arg, 
               " : ",
               defaultArgs[[arg]][['description']],
               " [default=",
               defaultArgs[[arg]][['value']],
               "]\n"))
  }
  quit(save="no")
}

user_args = commandArgs(trailingOnly=TRUE)

args <- list()
for (i in user_args){
  key_val = strsplit(i, split="=", fixed=TRUE)[[1]]
  if (length(key_val) != 2) {
    if (length(key_val) == 1 & key_val %in% c('-h', '--help')){
        print_help(defaultArgs)
      } else {
        stop(paste0("Invalid option `", i, "`. Args must be of the form KEY=VAL."))
      }
  }
  if (!key_val[1] %in% names(defaultArgs)){
    print(paste0("WARNING: Invalid key provided : ", i, ""))
    next
  } else {
    if (key_val[1] %in% names(args)){
      args[[key_val[1]]] = c(args[[key_val[1]]], defaultArgs[[key_val[1]]][['class']](key_val[2]))
    } else {
      args[[key_val[1]]] = defaultArgs[[key_val[1]]][['class']](key_val[2])
    }
  }
}

for (i in names(defaultArgs)) {
  if (i == 'SAMPLE_YML') {
    next
  }
  if (!i %in% names(args)){
    args[[i]] <- defaultArgs[[i]][['value']]
  }
}

print("Running with arguments:")
print(paste0("SAMPLE_YML     : ", args$SAMPLE_YML))
print(paste0("OUT_FOLDER     : ", args$OUT_FOLDER))
print(paste0("RSCRIPTS_DIR   : ", args$RSCRIPTS_DIR))
print(paste0("GENESET_DIR    : ", args$GENESET_DIR))
print(paste0("SPECIES        : ", args$SPECIES))
print(paste0("RERUN_STAGE    : ", args$RERUN_STAGE))
print(paste0("AB_ASSAY_NAME  : ", args$AB_ASSAY_NAME))

if (is.null(args$SAMPLE_YML)) {
  print_help(defaultArgs)
} else if (!file.exists(args$SAMPLE_YML)){
  stop("SAMPLE_YML did not exist.")
} else {
  samples <- read_yaml(args$SAMPLE_YML)
}

species_args <- list()
if (args$SPECIES == "human"){
  species_args$vars_to_regress <- c("percent.mt", "percent.ribo", "S.Score", "G2M.Score")
  species_args$mito_regex <- "^MT-"
  species_args$reference_dir <- paste(args$GENESET_DIR,"GRCh38", sep="/")
  species_args$CD45 <- "PTPRC"
} else if (args$SPECIES == "mouse"){
  species_args$vars_to_regress <- c("percent.mt", "percent.ribo", "percent.cc")
  species_args$mito_regex <- "^mt-"
  species_args$reference_dir <- paste(args$GENESET_DIR,"GRCm38", sep="/")
  species_args$CD45 <- "Ptprc"
} else {
  stop("SPECIES must be one of `human` or `mouse`", call.=FALSE)
}

if (!args$RERUN_STAGE %in% c(TRUE, FALSE)){
  stop("RERUN_STAGE must be one of `TRUE` or `FALSE`", call.=FALSE)
}

setwd(args$OUT_FOLDER)

suppressPackageStartupMessages({
  source(paste(args$RSCRIPTS_DIR, "identify_hto_clusters.R", sep="/"))
  source(paste(args$RSCRIPTS_DIR, "generate_profile_plot.R", sep="/"))
  source(paste(args$RSCRIPTS_DIR, "demux_HTOs_by_inflexions.R", sep="/"))
  source(paste(args$RSCRIPTS_DIR, "parse_10x_vdj_outs.R", sep="/"))
  source(paste(args$RSCRIPTS_DIR, "TriplePlot.R", sep="/"))
})

remove_rplots <- TRUE
if (file.exists("Rplots.pdf")){
  # Don't delete an existing file. Just delete if we create it.
  remove_rplots <- FALSE
}

sobjs <- list()
for (i in names(samples)){
  # Reset the value with each sample
  rerun_stage <- args$RERUN_STAGE
  if (!file.exists(paste0(i, "/", i, "_scTransformed_processed.RData")) | rerun_stage) {
    if (file.exists(paste0(i, "/", i, "_scTransformed_processed.RData"))) {
      print(paste0("Found `", i, "_scTransformed_processed.RData` but regenerating it as per request."))
      rerun_stage <- FALSE
    }
    if (!file.exists(paste0(i, "/", i, "_scTransformed.RData")) | rerun_stage) {
      if (file.exists(paste0(i, "/", i, "_scTransformed.RData"))) {
        print(paste0("Found `", i, "_scTransformed.RData` but regenerating it as per request."))
        rerun_stage <- FALSE
      }
      if (!file.exists(paste0(i, "/", i, "_filtered_singlets.RData")) | rerun_stage) {
        if (file.exists(paste0(i, "/", i, "_filtered_singlets.RData"))) {
          print(paste0("Found `", i, "_filtered_singlets.RData` but regenerating it as per request."))
          rerun_stage <- FALSE
        }
        if (!file.exists(paste0(i, "/", i, "_filtered.RData")) | rerun_stage) {
          if (file.exists(paste0(i, "/", i, "_filtered.RData"))) {
            print(paste0("Found `", i, "_filtered.RData` but regenerating it as per request."))
            rerun_stage <- FALSE
          }
          if (!file.exists(paste0(i, "/", i, "_raw.RData")) | rerun_stage) {
            if (file.exists(paste0(i, "/", i, "_raw.RData"))) {
              print(paste0("Found `", i, "_raw.RData` but regenerating it as per request."))
              rerun_stage <- FALSE
            }
            print(paste0("Generating `", i, "_raw.RData` ..."))
            datadir = file.path(samples[[i]][['datadir']], "raw_feature_bc_matrix")
            if (!dir.exists(datadir)) {
              cat(paste0("Could not find a directory within `datadir` for (", i, "). ",
                         "Tried to access (", datadir, ").\n"))
              next
            }

            metadata <- list()
            metadata[['tcr']] = data.frame()
            if (!is.null(samples[[i]][['tcrdir']])){
              tcrdir <- file.path(samples[[i]][['tcrdir']])
              if(dir.exists(tcrdir)){
                  cat(paste0("Reading TCR information from ", tcrdir, '\n'))
                  metadata[['tcr']] <- parse_tcr_clonotype(tcrdir)
                } else {
                  cat(paste0("`tcrdir` specified for ", i, " was not found at the ",
                             "specified location ", tcrdir, '\n'))
                }
            }

            metadata[['bcr']] = data.frame()
            if (!is.null(samples[[i]][['bcrdir']])){
              bcrdir <- file.path(samples[[i]][['bcrdir']])
              if(dir.exists(bcrdir)){
                  cat(paste0("Reading BCR information from ", bcrdir, '\n'))
                  metadata[['bcr']] <- parse_bcr_clonotype(bcrdir)
                } else {
                  cat(paste0("`bcrdir` specified for ", i, " was not found at the ",
                             "specified location ", bcrdir, '\n'))
                }
            }

            for (md in names(samples[[i]][['metadata']])){
              fn = file.path(samples[[i]][['metadata']][[md]])
              if (file.exists(fn)){
                cat(paste0('Reading metadata file: ' ,md, '\n'))
                metadata[[md]] <- read.table(fn,
                                             sep = "\t",
                                             header = TRUE,
                                             row.names=1,
                                             stringsAsFactors = TRUE)
                rownames(metadata[[md]]) <- gsub("-1$", '', rownames(metadata[[md]]))
                if (dim(metadata[[md]])[2] < 1){
                  stop(paste0('Something wrong with metadata file ', md, '. Found 0 columns... Is the file a tsv?'))
                }
              } else {
                cat(paste0("`metadata` file ", md, " specified for ", i, 
                           " was not found at the specified location ", 
                           fn, '\n'))
              }
            }

            merged_metadata <- data.frame()
            for (md in names(metadata)){
              merged_metadata <- merge(merged_metadata, metadata[[md]], by=0, all.x=TRUE, all.y=TRUE)
              rownames(merged_metadata) <- merged_metadata$Row.names
              merged_metadata$Row.names <- NULL
            }

            if (dim(merged_metadata)[1] == 0) {
              # Empty data frame. NULL is the default for metadata in a seurat object
              merged_metadata = NULL
            }

            # Create the fodler if it doesn't exist
            dir.create(i, showWarnings=FALSE)

            data <- Read10X(data.dir=datadir)
            # newer versions of Seurat
            if (endsWith(x=colnames(data)[1], "-1")){
              rownames(merged_metadata) <- paste0(rownames(merged_metadata), '-1')
            }
            if (typeof(data) == "list"){
              sobjs[[i]] <- CreateSeuratObject(counts = data$`Gene Expression`,
                                               project = i,
                                               min.cells = 3,
                                               min.features = 100,
                                               meta.data=merged_metadata)
              ab_capture <- CreateAssayObject(counts = data$`Antibody Capture`)
              sobjs[[i]][[args$AB_ASSAY_NAME]] <- subset(ab_capture, cells=colnames(sobjs[[i]]))
              rm(list=c("ab_capture", "data"))
            } else {
              sobjs[[i]] <- CreateSeuratObject(counts = data,
                                               project = i,
                                               min.cells = 3,
                                               min.features = 100,
                                               meta.data=merged_metadata)
            }
            # store mitochondrial percentage in object metadata
            sobjs[[i]] <- PercentageFeatureSet(sobjs[[i]],
                                               pattern = species_args$mito_regex,
                                               col.name = "percent.mt")
            #Store ribosomal percentage in the object metadata
            ribo_genes <- read.table(paste(species_args$reference_dir,
                                         "ribo_genes.tsv", sep="/"),
                                     sep = "\t", header=TRUE,
                                     stringsAsFactors = FALSE)
            ribo_genes <- ribo_genes[ribo_genes[["HUGO"]] %in% rownames(sobjs[[i]]), ]
            sobjs[[i]] <- PercentageFeatureSet(sobjs[[i]],
                                               features = ribo_genes[["HUGO"]],
                                               col.name = "percent.ribo")

            #Store cell cycle state in the object metadata
            cc_genes <- read.table(paste(species_args$reference_dir,
                                         "cell_cycle_genes.tsv", sep="/"),
                                   sep = "\t", header=TRUE,
                                   stringsAsFactors = FALSE)

            if (args$SPECIES == "human"){
              sobjs[[i]] <- CellCycleScoring(sobjs[[i]],
                                             s.features = cc_genes[cc_genes$stage=="G1-S", "HUGO"],
                                             g2m.features = cc_genes[cc_genes$stage=="G2-M", "HUGO"],
                                             nbin = 12)
            } else if (args$SPECIES == "mouse") {

              cc_genes <- cc_genes[cc_genes[["HUGO"]] %in% rownames(sobjs[[i]]), ]
              sobjs[[i]] <- PercentageFeatureSet(sobjs[[i]],
                                                 features = cc_genes[["HUGO"]],
                                                 col.name = "percent.cc")
            }

            plot_all_profiles(sobjs[[i]],
                              out_prefix=paste0(i, "/", i, "_dualscatter_pre"),
                              plot_extra=FALSE,
                              ab_assay_name=args$AB_ASSAY_NAME)

            cutoffs = list()
            for (hl in c('high', 'low')){
              for (feature in c('percent.mt', 'percent.ribo')){
                cutoffs[[paste0(feature, '.', hl)]] = NA
              }
              for (assay in Assays(sobjs[[i]])){
                cutoffs[[paste0('nCount_', assay, '.', hl)]] = NA
                cutoffs[[paste0('nFeature_', assay, '.', hl)]] = NA
              }
            }

            if ("DROPLET.TYPE" %in% colnames(sobjs[[i]]@meta.data)) {
              # We have passed in Demuxlet/Freemuxlet info
              plot_all_profiles(sobjs[[i]],
                                out_prefix=paste0(i, "/", i, "_droplet_dualscatter_pre"),
                                plot_extra=FALSE,
                                scatter_color="DROPLET.TYPE",
                                ab_assay_name=args$AB_ASSAY_NAME)
              cutoffs[["DROPLET.TYPE.keep"]] = NA
              cutoffs[["BEST.GUESS.keep"]] = NA
            }

            write_yaml(cutoffs, file=paste0(i, "/cutoffs.yml"))

            assign(i, sobjs[[i]])
            save(list=i, file=paste0(i, "/", i, "_raw.RData"))
            next
          } else {
            print(paste0("Using existing `", i, "_raw.RData`."))
            load(paste0(i, "/", i, "_raw.RData"))
            sobjs[[i]] <- get(i)
          }  # END STAGE 1
          rm(list=i)
          print(paste0("Generating `", i, "_filtered.RData` ..."))
          if (!file.exists(paste0(i, "/cutoffs.yml"))) {
            print("Cannot continue without a cutoffs.yml file in the output directory")
            next
          }
          cutoffs <- read_yaml(paste0(i, "/cutoffs.yml"))
          # Filter cells with high mito content (dying/dead cells)
          keep_rownames = rep(TRUE, dim(sobjs[[i]]@meta.data)[1])
          total_cells = c(dim(sobjs[[i]]@meta.data)[1], 100)
          additional_text = ""

          for (filter_key in names(cutoffs)){
            filter_feature = strsplit(filter_key, split = ".",
                                      fixed = TRUE)[[1]]
            suppmsg <- assert_that(length(filter_feature) >= 2,
                                   msg=paste0("Names in cutoffs.yml must ",
                                              "be of the form `feature.high`, ",
                                              "`feature.low`, or `feature.keep`. ",
                                              "Got ", filter_key))
            filter_cat <- filter_feature[length(filter_feature)]
            filter_feature <- paste(filter_feature[1:(length(filter_feature)-1)],
                                    collapse=".")
            suppmsg <- assert_that(filter_cat %in% c("high", "low", "keep"),
                                   msg=paste0("Names in cutoffs.yml must ",
                                              "be of the form `feature.high`, ",
                                              "`feature.low`, or `feature.keep`. ",
                                              "Got ", filter_key))
            if (!filter_feature %in% colnames(sobjs[[i]]@meta.data)){
              print(paste0("Could not find ", filter_feature, " in the metadata ",
                           "for ", i))
              next
            }
            if(!is.na(cutoffs[[filter_key]])){
              if (filter_cat == "high") {
                keep_rownames <- keep_rownames & sobjs[[i]]@meta.data[[filter_feature]] <= cutoffs[[filter_key]]
                filter_cat_text <- " > "
              } else if (filter_cat == "low") {
                keep_rownames <- keep_rownames & sobjs[[i]]@meta.data[[filter_feature]] >= cutoffs[[filter_key]]
                filter_cat_text <- " < "
              } else {
                keep_rownames <- keep_rownames &
                  sobjs[[i]]@meta.data[[filter_feature]] %in% strsplit(cutoffs[[filter_key]], ',')[[1]]
                filter_cat_text <- " in "
              }

              print(paste0("Dropping ", additional_text,
                           total_cells[1] - sum(keep_rownames),
                           " cells (",
                           round((total_cells[1] - sum(keep_rownames))/length(keep_rownames) * 100, 2),
                           "%) for having `", filter_feature,
                           filter_cat_text,
                           cutoffs[[filter_key]], "`."))
              additional_text = "an additional "
              total_cells = c(sum(keep_rownames),
                              round(sum(keep_rownames)/length(keep_rownames) * 100, 2))
            }
          }

          print(paste0("Dropping a total of ", sum(keep_rownames==FALSE),
                       " cells (",
                       round(sum(keep_rownames==FALSE)/length(keep_rownames) *100, 2),
                       "%)."))
          print(paste0("Retaining a total of ", sum(keep_rownames),
                       " cells (",
                       round(sum(keep_rownames)/length(keep_rownames) *100, 2),
                       "%)."))

          sobjs[[i]] <- subset(sobjs[[i]],
                               cells = rownames(sobjs[[i]]@meta.data[keep_rownames, ]))
          keep_features <- c()
          for (assay in names(sobjs[[i]]@assays)){
            if (assay == "RNA"){
              next
            }
            keep_features <- c(keep_features, rownames(sobjs[[i]]@assays[[assay]]@counts))
          }
          keep_genes <- rownames(sobjs[[i]]@assays$RNA@counts)[Matrix::rowSums(sobjs[[i]]@assays$RNA@counts>0)>3]
          print(paste0("Dropping ",
                       (length(rownames(sobjs[[i]]@assays$RNA@counts)) -
                        length(keep_genes)),
                       " genes for being in fewer than 3 cells"))
          sobjs[[i]] <- subset(sobjs[[i]], features = c(keep_genes, keep_features))

          plot_all_profiles(sobjs[[i]],
                            out_prefix=paste0(i, "/", i, "_dualscatter_post"),
                            plot_extra=TRUE,
                            ab_assay_name=args$AB_ASSAY_NAME)

          if ("DROPLET.TYPE" %in% colnames(sobjs[[i]]@meta.data)) {
            # We have passed in Demuxlet/Freemuxlet info
            plot_all_profiles(sobjs[[i]],
                              out_prefix=paste0(i, "/", i, "_droplet_dualscatter_post"),
                              plot_extra=TRUE,
                              scatter_color="DROPLET.TYPE",
                              ab_assay_name=args$AB_ASSAY_NAME)
            sobjs[[i]]@meta.data$SAMPLE.by.SNPs = as.factor(sapply(as.vector(sobjs[[i]]@meta.data$BEST.GUESS), function(x) {
              temp <- paste(sort(unique(strsplit(x, ',')[[1]])), collapse='_')
              }))
          }

          # Now that we"ve removed bogus cells, let"s normalize the counts (if present)
          if (args$AB_ASSAY_NAME %in% names(sobjs[[i]]@assays)){
            sobjs[[i]] <- NormalizeData(sobjs[[i]], assay = args$AB_ASSAY_NAME,
                                        normalization.method = "CLR")
            sobjs[[i]] <- ScaleData(sobjs[[i]], assay = args$AB_ASSAY_NAME)
          }
          assign(i, sobjs[[i]])
          save(list=i, file=paste0(i, "/", i, "_filtered.RData"))
          next
        } else {
          print(paste0("Using existing `", i, "_filtered.RData`."))
          load(paste0(i, "/", i, "_filtered.RData"))
          sobjs[[i]] <- get(i)
        }  # END STAGE 2
        rm(list=i)
        if (args$AB_ASSAY_NAME %in% names(sobjs[[i]]@assays)){
          print(paste0("Generating `", i, "_filtered_singlets.RData` ..."))
          if (!file.exists(paste0(i, "/IDX_map.tsv"))){
            print("Cannot continue without an IDX_map.tsv file")
            next
          }
          IDX_map <- read.table(paste0(i, "/IDX_map.tsv"), sep="\t",
                                header=TRUE, row.names=1, stringsAsFactors=FALSE)
          # Hacky but it works ¯\_(ツ)_/¯
          IDX_map <- sapply(colnames(IDX_map), function(x) {
                              y = IDX_map[,x, drop=T]
                              names(y) = rownames(IDX_map)
                              y
                              }, simplify = FALSE, USE.NAMES = TRUE)
          sample_names <- names(IDX_map$sample_name)
          sample_names <- sample_names[!sample_names %in% c("NEGATIVE", "MULTIPLET", "OTHER")]
          suppmgg <- assert_that(all(sample_names %in% rownames(sobjs[[i]]@assays[[args$AB_ASSAY_NAME]])),
                                 msg="Not all hashtags in IDX_map.tsv are in the seurat object")
          keep_features <- c(rownames(sobjs[[i]]@assays$RNA),
                             sample_names)
          sobjs[[i]] <- subset(sobjs[[i]], features = keep_features)

          if ("threshold" %in% names(IDX_map)){
            inflexions <- IDX_map$threshold[sample_names]
            if (any(is.na(inflexions))){
              inflexions <- NULL
            } else {
              print("Using user-provided thresholds for demultiplexing.")
            }
          } else {
            inflexions <- NULL
          }

          if ("color" %in% names(IDX_map)){
            sample_colors <- IDX_map$color
            suppmsg <- assert_that(all(c("NEGATIVE", "MULTIPLET", "OTHER") %in% names(sample_colors)),
                                   msg=paste0("Cannot use user-provided colors if values ",
                                              "for NEGATIVE,  MULTIPLET, and OTHER are not ",
                                              "provided"))
            if (any(is.na(sample_colors))){
              sample_colors <- NULL
            } else {
              #TODO
              sample_colors <- NULL
              print("Can't use user-provided colors for plotting yet.")
              #print("Using user-provided colors for plotting.")
            }
          } else {
            sample_colors <- NULL
          }

          demux_results <- demux_by_inflexions(sobjs[[i]],
                                               sample_names=IDX_map$sample_name[sample_names],
                                               inflexions=inflexions,
                                               sample_colors=sample_colors,
                                               assay_name=args$AB_ASSAY_NAME)

          nplots <- length(demux_results$background_plots)

          png(paste0(i, "/", i, "_raw_hto_pairplots.png"), width=nplots*750,
              height=nplots*750, units="px")
          print(plot_grid(plotlist=unlist(demux_results[['raw_plots']], recursive = FALSE),
                          ncol=nplots))
          dev.off()
          png(paste0(i, "/", i, "_scaled_hto_pairplots.png"), width=nplots*750,
              height=nplots*750, units="px")
          print(plot_grid(plotlist=unlist(demux_results[['scaled_plots']], recursive=FALSE),
                          ncol=nplots))
          dev.off()
          png(paste0(i, "/", i, "_hto_background_profiles.png"), width=1500,
              height=750*max(1, floor(nplots/2)), units="px")
          print(plot_grid(plotlist=demux_results[['background_plots']], ncol=2))
          dev.off()

          sobjs[[i]]@meta.data$SAMPLE.by.ABs <- factor(demux_results$cell_ids,
                                                     levels=c(unname(sort(IDX_map$sample_name)), 'MULTIPLET', 'NEGATIVE'))

          drf <- data.frame(sample_name=unlist(IDX_map$sample_name[sample_names]),
                  threshold=demux_results[['inflexions']])
          write.table(drf, file=paste0(i, "/IDX_map_generated.tsv"),
                      quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

          capture.output(demux_results[['secondary_inflexions']],
                         file=paste0(i, "/IDX_additional_inflexions.txt"))

          write.table(demux_results[["summary"]], file=paste0(i, "/IDX_summary.tsv"),
                      quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")

          Idents(sobjs[[i]]) <- sobjs[[i]]@meta.data$SAMPLE.by.ABs
          png(paste0(i, "/", i, '_ridgeplot.png'), width=1500,
              height=ceiling(length(sample_names)/2)*750, units = 'px')
          print(RidgePlot(sobjs[[i]],
                          assay=args$AB_ASSAY_NAME,
                          features=sample_names,
                          ncol = 2))
          dev.off()
          sobjs[[i]] <- subset(sobjs[[i]],
                               cells=colnames(sobjs[[i]])[!sobjs[[i]]$SAMPLE.by.ABs %in% c('NEGATIVE', 'MULTIPLET')])
          sobjs[[i]]$SAMPLE.by.ABs <- factor(sobjs[[i]]$SAMPLE.by.ABs,
                                             levels=unname(sort(IDX_map$sample_name)))
          assign(i, sobjs[[i]])
          save(list=i, file=paste0(i, "/", i, "_filtered_singlets.RData"))
          next
        } else {
          print(paste0("Could not find an assay named ", args$AB_ASSAY_NAME,
                       " in the object. Continuing on to run SCTransform. ",
                       "If you have multiplexed this data and have stored the ",
                       "counts in another assay name, delete the SCTransformed ",
                       "RData rerun this stage with AB_ASSAY_NAME=XXX."))
        }
      } else {
        print(paste0("Using existing `", i, "_filtered_singlets.RData`."))
        load(paste0(i, "/", i, "_filtered_singlets.RData"))
        sobjs[[i]] <- get(i)
      }  # END STAGE 3
      rm(list=i)
      print(paste0("Generating `", i, "_scTransformed.RData` ..."))
      sobjs[[i]] <- SCTransform(sobjs[[i]],
                                vars.to.regress = species_args$vars_to_regress,
                                return.only.var.genes = FALSE,
                                verbose = FALSE)

      # Get the Principal components for the object
      sobjs[[i]] <- RunPCA(sobjs[[i]], verbose = FALSE)
      # Verify the confounding vars have been regressed out
      pdf(paste0(i, "/", i, "_confound_pca.pdf"))
      print(FeaturePlot(sobjs[[i]], reduction = "pca",
                        features=species_args$vars_to_regress))
      dev.off()
      assign(i, sobjs[[i]])
      save(list=i, file=paste0(i, "/", i, "_scTransformed.RData"))
      next
    } else {
      print(paste0("Using existing `", i, "_scTransformed.RData`."))
      load(paste0(i, "/", i, "_scTransformed.RData"))
      sobjs[[i]] <- get(i)
    }  # END STAGE 4
    rm(list=i)
    print(paste0("Generating `", i, "_scTransformed_processed.RData`..."))
    sobjs[[i]] <- RunUMAP(sobjs[[i]],
                     dims = 1:30,  # Num PCs to use
                     n.neighbors = 30,  # Default. Controls how UMAP balances local (low) versus global (large) structure in the data
                     min.dist = 0.3,   # Default. Controls the size of the clusters. Should be smaller than spread
                     spread = 1,  # Default. Controls the inter-cluster distances to some extent. Should be larger than min_dist
                     a = NULL,  # Default. Can be used with b instead of using min.dist/spread
                     b = NULL,  # Default. Can be used with a instead of using min.dist/spread
                     verbose = FALSE)
    pdf(paste0(i, "/", i, "_umap.pdf"))
    print(DimPlot(sobjs[[i]]))
    dev.off()

    # Verify the confounding vars have been regressed out
    pdf(paste0(i, "/", i, "_confound_umap.pdf"))
    print(FeaturePlot(sobjs[[i]], features=species_args$vars_to_regress))
    dev.off()
    # Calculate the neighborhood graph
    sobjs[[i]] <- FindNeighbors(sobjs[[i]],
                          dims = 1:30,  # Num PCs to use
                          k.param = 20,  # k for the knn algorithm
                          verbose = FALSE)
    # Use the neighborhood graph to cluster the data
    sobjs[[i]] <- FindClusters(sobjs[[i]], verbose = FALSE,
                               algorithm = 4)  # Use Leiden

    # View the cluters UMAP
    pdf(paste0(i, "/", i, "_umap.pdf"))
    print(DimPlot(sobjs[[i]], label=TRUE))
    dev.off()

    if ("SAMPLE.by.SNPs" %in% colnames(sobjs[[i]]@meta.data)) {
      # We have passed in Demuxlet/Freemuxlet info
      pdf(paste0(i, "/", i, "_samples_by_SNPs_umap.pdf"))
      print(DimPlot(sobjs[[i]], group.by="SAMPLE.by.SNPs"))
      dev.off()
    }

    if ("SAMPLE.by.ABs" %in% colnames(sobjs[[i]]@meta.data)){
      # View the clusters UMAP
      pdf(paste0(i, "/", i, "_samples_by_ABs_umap.pdf"))
      print(DimPlot(sobjs[[i]], group.by="SAMPLE.by.ABs"))
      dev.off()
    }

    pdf(paste0(i, "/", i, "_immune_triplot_umap.pdf"))
    tryCatch(expr=print(TriPlot(sobjs[[i]], features=species_args$CD45, reduction.use="umap", group.by="seurat_clusters")),
             error=function(e) {print(paste0("WARNING: Could not print CD45(",
                                species_args$CD45,
                                ") expression since it was not found in the object"))
                                ggplot() + theme_void()},
             finally=function(e){dev.off()})
    assign(i, sobjs[[i]])
    save(list=i, file=paste0(i, "/", i, "_scTransformed_processed.RData"))
    next
  } else {
    print(paste0("Nothing more to do for ", i, " since `", i,
                 "_scTransformed_processed.RData` exists."))
  }
}  # END STAGE 5

if (remove_rplots & file.exists("Rplots.pdf")){
  file.remove("Rplots.pdf")
}

uuid <- paste(sample(c(LETTERS,0:9), 10, replace = T), collapse="")
write(paste(names(warnings()), collapse="\n"), file=paste0("run_warnings_", uuid, ".txt"))
print(paste0("Wrote run warnings to run_warnings_", uuid, ".txt"))
