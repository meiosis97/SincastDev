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
#' @param cells To be added.
#' @param features To be added.
#' @param group.by To be added.
#' @param sample.method To be added.
#' @param n.pool To be added.
#' @param size.factor To be added.
#' @param pool.factor To be added.
#' @param sep To be added.
#' @param replace Whether to replace the existing 'pseudobulk' assay.
#'
#' @return A \code{Seurat} object with updated \code{Sincast pseudobulk} assay.
#'
#' @export
#' @rdname SincastAggregate
#' @aliases Sincast, SincastAssays, Seurat
setGeneric("SincastAggregate", function(object,
                                        assay = NULL,
                                        layer = "counts",
                                        cells = NULL,
                                        features = NULL,
                                        group.by = "ident",
                                        sample.method = c("equal.size", "equal.prop"),
                                        aggregate.method = c("addition", "average"),
                                        n.pool = 15,
                                        size.factor = 1,
                                        pool.factor = NULL,
                                        sep = ".", replace = FALSE, ...) {
  standardGeneric("SincastAggregate")
})

#' @rdname SincastAggregate
setMethod("SincastAggregate", "Seurat", function(object,
                                                 assay = NULL,
                                                 layer = "counts",
                                                 cells = NULL,
                                                 features = NULL,
                                                 group.by = "ident",
                                                 sample.method = c("equal.size", "equal.prop"),
                                                 aggregate.method = c("addition", "average"),
                                                 n.pool = 15,
                                                 size.factor = 1,
                                                 pool.factor = NULL,
                                                 sep = ".", replace = FALSE, ...) {

  # Check the validity of the "Sincast" object in "Seurat"'s "misc" slot.
  test.SincastObject <- Sincast::CheckSincastObject(object, complete = FALSE)

  # If the 'Sincast' object is missing or unrecognized, create an empty 'Sincast' object.
  if (test.SincastObject != "Valid") {
    message("SincastAggregate: Create an empty 'Sincsat' object.")
    SincastObject <- Sincast::CreateSincastObject()
  } else {
    SincastObject <- Sincast::GetSincastObject(object)
  }

  # If the 'Sincast' object is not of the class 'Sincast', either replace it or return an error.
  if (!replace & test.SincastObject == "Wrong class") {
    stop(
      "SincastAggregate: Set 'replace = T' to enforce a replacement for the unrecorgnized 'Sincast' object."
    )
  }

  SincastAssays <- SincastObject@SincastAssays

  # If a 'pseudobulk' assay already exists, either replace it or return an error.
  if (!is.null(SincastAssays@pseudobulk)) {
    if(!replace){
      stop(
        "SincastAggregate: A 'pseudobulk‘ assay already exists, set replace = T to enforce a replacement."
      )
    }else{
      message("SincastAggregate: A 'pseudobulk‘ assay already exists. Will be replaced as 'replace = T'.")
    }
  }

  # Get defaults assay.
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }

  # Get features to aggregate.
  if(!is.null(features)) object <- object[features, ]

  # Get cells to aggregate.
  if(!is.null(cells)) object <- object[,cells]

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

    # Generate pseudobulk.
    out[[g]] <- replicate(
      n = n.gen,
      expr = generate.one.pseudo.bulk(tmp.data, n.pool, aggregate.method)
    )
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
  SincastToken <- GenerateSincastToken(by = "SincastAggregate")
  ids <- c(names(SincastObject@SincastToken@command), SincastToken@id)
  SincastToken@command <- c(SincastObject@SincastToken@command, deparse(match.call()))
  names(SincastToken@command) <- ids

  # Create Seurat object: SincastPseudoBulk to store the result.
  suppressWarnings(
    SincastPseudobulk <- Seurat::CreateSeuratObject(
      counts = Matrix::Matrix(out),
      assay = assay,
      meta.data = meta.data,
      project = SincastToken@id, ...
    )
  )
  Seurat::Idents(SincastPseudobulk) <- "agg.label" # Add aggregated annotation.
  Seurat::Misc(SincastPseudobulk, slot = "SincastToken") <- SincastToken # Add token.

  # Add or replace the pseudobulk assay by SincastPseudobulk.
  message(
    "SincastAggregate: Add or replace the 'pseudobulk' assay in the 'SincastAssays'",
    "object by a new aggregation result."
  )

  SincastAssays@pseudobulk <- SincastPseudobulk
  SincastObject@SincastAssays <- SincastAssays
  Sincast::GetSincastObject(object) <- SincastObject

  # Calculate sparsity
  sparsity.before <-  1-mean(object$nFeature_RNA/nrow(object))
  sparsity.after <- mean(out == 0)

  message(
    "SincastAggregate: Sparsity before aggregation: ", round(sparsity.before,3),
    "; After aggregation: ", round(sparsity.after,3)
  )

  object
})
