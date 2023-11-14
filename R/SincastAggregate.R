# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Not exported satellite functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
generate.one.pseudo.bulk <- function(data, n.pool, aggregate.method) {
  idx <- sample(
    x = 1:ncol(data),
    size = n.pool,
    replace = T
  )
  if (aggregate.method == "addition") {
    rowSums(data[, idx, drop = F])
  } else {
    rowMeans(data[, idx, drop = F])
  }
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Aggregate single cells to generate pseudobulk samples.
#'
#' To be added.
#'
#' @param object A \code{Seurat} object
#' @param assay To be added.
#' @param layer To be added.
#' @param features To be added.
#' @param group.by To be added.
#' @param sample.method To be added.
#' @param n.pool To be added.
#' @param size.factor To be added.
#' @param pool.factor To be added.
#' @param sep To be added.
#' @param replace Whether to replace the existing 'pseudobulk' assay.
#'
#' @return An updated \code{Seurat} object with an additional \code{SincastAssays}
#' object in the \code{misc} slot. Aggregated pseudobulk samples are stored in
#' the Sincast assay \code{pseudobulk} embedded in \code{Seurat::MISC(object, slot = 'SincastAssays')}.
#'
#' @family SincastAssays related methods
#'
#' @export
#' @rdname SincastAggregate
#' @aliases Sincast, SincastAssays, Seurat
setGeneric("SincastAggregate", function(object,
                                        assay = NULL,
                                        layer = "counts",
                                        features = NULL,
                                        group.by = "ident",
                                        sample.method = c("equal.size", "equal.prop"),
                                        aggregate.method = c("addition", "average"),
                                        n.pool = 15,
                                        size.factor = 1,
                                        pool.factor = NULL,
                                        sep = ".", replace = c("None", "All", "invalid"), ...) {
  standardGeneric("SincastAggregate")
})

#' @rdname SincastAggregate
setMethod("SincastAggregate", "Seurat", function(object,
                                                 assay = NULL,
                                                 layer = "counts",
                                                 features = NULL,
                                                 group.by = "ident",
                                                 sample.method = c("equal.size", "equal.prop"),
                                                 aggregate.method = c("addition", "average"),
                                                 n.pool = 15,
                                                 size.factor = 1,
                                                 pool.factor = NULL,
                                                 sep = ".", replace = FALSE, ...) {

  object <- Sincast::CreateSincastAssays(object)

  # Check the existence of the 'pseudobulk' assay.
  SincastAssays <- Seurat::Misc(object, slot = "SincastAssays")

  if(is(SincastAssays, 'SincastAssays')){
    if(!replace & !is.null(SincastAssays@pseudobulk)){
      stop("A 'pseudobulk‘ assay already exists, set replace = T to enforce a replacement.")
    }

  }else{
    if(!replace){
      stop("Invalid 'SincastAssays' object, set replace = T to enforce a replacement.")
    }else{
      SincastAssays <- new('SincastAssays')
    }

  }



  # Get defaults assay.
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
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

  # Check aggregate.method format.
  aggregate.method <- match.arg(aggregate.method)

  # Store pseudo-bulk.
  out <- list()
  agg.label <- c()

  # Aggregation.
  for (g in group.name) {
    message(paste("\rSincastAggregate: Now aggregate ", g), appendLF = F)

    # Calculate how many cells are in cluster g.
    n.c <- sum(group == g)

    # Calculate how many to generate for cluster g.
    if (sample.method == "equal.size") {
      n.gen <- round(size.factor) + 1
    } else {
      n.gen <- round(n.c * size.factor) + 1
    }

    # Calculate how many to pool for cluster g.
    if (!is.null(pool.factor)) {
      n.pool <- round(n.c * pool.factor) + 1
    }

    # Subset data.
    tmp.data <- SeuratObject::LayerData(
      object = object,
      layer = layer,
      assay = assay
    )[, group == g]

    # Generate pseudo-bulks.
    out[[g]] <- replicate(n = n.gen, expr = generate.one.pseudo.bulk(tmp.data, n.pool, aggregate.method))
    colnames(out[[g]]) <- paste(g, sep, 1:n.gen, sep = "")
    agg.label <- c(agg.label, rep(g, n.gen))
  }
  message("\n", rep("-", 50))

  # Collapse the list.
  out <- unname(out)
  out <- do.call(cbind, out)
  meta.data <- data.frame(agg.label)
  rownames(meta.data) <- colnames(out)


  # Generate a token for this aggregation run.
  SincastToken <- GenerateSincastToken()
  # Create Seurat object: SincastPseudoBulk to store the result.
  suppressWarnings(
    SincastPseudoBulk <- Seurat::CreateSeuratObject(
      counts = Matrix::Matrix(out),
      assay = assay,
      meta.data = meta.data,
      project = SincastToken@timestamp, ...
    )
  )
  Seurat::Idents(SincastPseudoBulk) <- "agg.label" # Add aggregated annotation.
  Seurat::Misc(SincastPseudoBulk, slot = "SincastToken") <- SincastToken # Add token.

  # Generate command log.
  SincastPseudoBulk@commands$CreateSincastAssays <- deparse(match.call())

  # Append the "pseudobulk" assay by "SincastPseudoBulk".
  message("SincastAggregate: Append the 'pseudobulk' assay in the 'SincastAssays' object by a new aggregation result.")
  suppressWarnings(Seurat::Misc(object, slot = "SincastAssays")@pseudobulk <- SincastPseudoBulk)

  # Summary
  message(
    "SincastAggregate: 'SincastAssays' object contains ",
    length(pseudobulk)," valid pseudobulk objects and ",
    length(Seurat::Misc(object, slot = "SincastAssays")@imputation), " valid imputation objects."
  )

  object
})
