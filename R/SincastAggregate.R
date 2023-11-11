generate.one.pseudo.bulk <- function(data, n.pool) {
  idx <- sample(
    x = 1:ncol(data),
    size = n.pool,
    replace = T
  )
  rowMeans(data[, idx, drop = F])
}

SincastAggregate <- function(object,
                             assay = NULL,
                             layer = "counts",
                             features = NULL,
                             group.by = "ident",
                             sample.method = c("equal.size", "equal.prop"),
                             n.pool = 15,
                             size.factor = 1,
                             pool.factor = NULL,
                             sep = ".", ...) {
  # Get defaults assay.
  if (is.null(assay)) {
    assay <- DefaultAssay(object)
  }

  # Get default features to aggregate.
  if (is.null(features)) {
    features <- rownames(object[[assay]])
  }

  # Get default group.
  if (group.by == "ident") {
    group <- object@active.ident
  } else {
    group <- object@meta.data[, group.by]
  }
  group.name <- unique(group) %>% sort()

  # Check sample.method format.
  sample.method <- match.arg(sample.method)

  # Store pseudo-bulk.
  out <- list()
  agg.label <- c()

  # Aggregation.
  for (g in group.name) {
    message(paste("\r Now aggregate", g), appendLF = F)

    # Calculate how many cells are in cluster g.
    n.c <- sum(group == g)

    # Calculate how many to generate for cluster g.
    if (sample.method == "equal.size") {
      n.gen <- ceiling(size.factor)
    } else {
      n.gen <- ceiling(n.c * size.factor)
    }

    # Calculate how many to pool for cluster g.
    if (!is.null(pool.factor)) {
      n.pool <- ceiling(n.c * pool.factor)
    }

    # Subset data.
    tmp.data <- LayerData(
      object = object,
      layer = layer,
      assay = assay
    )[, group == g]

    # Generate pseudo-bulks.
    out[[g]] <- replicate(n = n.gen, expr = generate.one.pseudo.bulk(tmp.data, n.pool))
    colnames(out[[g]]) <- paste(g, sep, 1:n.gen, sep = "")
    agg.label <- c(agg.label, rep(g, n.gen))
  }

  # Collapse the list
  out <- unname(out)
  out <- do.call(cbind, out)
  meta.data <- data.frame(agg.label)
  rownames(meta.data) <- colnames(out)

  # Create Seurat object
  SincastPseudoBulk <- Seurat::CreateSeuratObject(
    counts = Matrix::Matrix(out),
    assay = assay,
    meta.data = meta.data,
    project = "Sincast", ...
  )

  Idents(SincastPseudoBulk) <- "agg.label"
  SincastPseudoBulk
}
