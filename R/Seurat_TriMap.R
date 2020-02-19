RunTRIMAP.Seurat <- function(
  object,
  dims = NULL,
  reduction = 'pca',
  features = NULL,
  assay = 'RNA',
  n.inliers = 10,
  n.outliers = 5,
  n.random = 5,
  distance = 'euclidean',
  weight.adj = 500.0,
  learning.rate = 1000,
  n.iters = 400,
  apply.pca = TRUE,
  opt.method = 'dbd',
  seed.use = 42,
  verbose = TRUE,
  reduction.name = 'trimap',
  reduction.key = 'TriMap_',
  ...
) {
  Seurat:::CheckDots(...)
  if (sum(c(is.null(x = dims), is.null(x = features))) != 1) {
    stop("Please specify only one of the following arguments: dims, or features")
  }
  if (!is.null(x = features)) {
    data.use <- t(x = GetAssayData(object = object, slot = 'data', assay = assay)[features, ])
  } else if (!is.null(x = dims)) {
    data.use <- Embeddings(object[[reduction]])[, dims]
    assay <- DefaultAssay(object = object[[reduction]])
    apply.pca = FALSE
  } else {
    stop("Please specify either dims or features")
  }
  object[[reduction.name]] <- RunTRIMAP.default(
    object = data.use,
    assay = assay,
    n.inliers = n.inliers,
    n.outliers = n.outliers,
    n.random = n.random,
    distance = distance,
    weight.adj = weight.adj,
    learning.rate = learning.rate,
    n.iters = n.iters,
    apply.pca = apply.pca,
    opt.method = opt.method,
    seed.use = seed.use,
    reduction.key = reduction.key,
    verbose = verbose
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}


RunTRIMAP.default <- function(
  object,
  assay = NULL,
  n.inliers = 10,
  n.outliers = 5,
  n.random = 5,
  distance = 'euclidean',
  weight.adj = 500.0,
  learning.rate = 1000,
  n.iters = 400, 
  apply.pca = TRUE,
  opt.method = 'dbd',
  seed.use = 42,
  reduction.key = 'TriMap_',
  verbose = TRUE,
  ...
) {
  Seurat:::CheckDots(...)
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (!py_module_available(module = 'trimap')) {
    stop("Cannot find trimap, please install through pip (e.g. pip install trimap).")
  }
  if (!is.null(x = seed.use)) {
    py_set_seed(seed = seed.use)
  }
  if (typeof(x = n.epochs) == "double") {
    n.epochs <- as.integer(x = n.iters)
  }
  trimap_import <- import(module = "trimap", delay_load = TRUE)
  trimap <- trimap_import$TRIMAP(
    n_inliers = as.integer(x = n.inliers),
    n_outliers = as.integer(x = n.outliers),
    n_random = as.integer(x = n.random),
    distance = distance,
    weight_adj = as.double(x = weight.adj),
    lr = as.double(x = learning.rate),
    n_iters = n.epochs,
    apply_pca = apply.pca,
    opt_method = opt.method,
    return_seq = FALSE,
    verbose = verbose
  )
  trimap.output <- trimap$fit_transform(as.matrix(x = object))
  
  colnames(x = trimap.output) <- paste0(reduction.key, 1:ncol(x = trimap.output))
  if (inherits(x = object, what = 'dist')) {
    rownames(x = trimap.output) <- attr(x = object, "Labels")
  } else {
    rownames(x = trimap.output) <- rownames(x = object)
  }
  trimap.reduction <- CreateDimReducObject(
    embeddings = trimap.output,
    key = reduction.key,
    assay = assay
  )
  return(trimap.reduction)
}
