# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Not exported satellite functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GeneratePseudobulk <- function(data, n.pool, aggregate.method) {
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
#' Perform \code{Sincast} aggregation.
#'
#' @param object A \code{Sincast} object.
#' @param assay Which \code{Seurat} assay to use. Default is the default \code{Seurat} assay.
#' @param layer Which \code{Seurat} layer to use. Default is the \code{counts} layer.
#' @param features Features to analyze. Default is all features in the assay.
#' @param group.by To be added.
#' @param sample.method To be added.
#' @param n.pool To be added.
#' @param size.factor To be added.
#' @param pool.factor To be added.
#' @param sep To be added.
#' @param replace Logical; if TRUE, replace the existing \code{pseudobulk} assay.
#'
#' @return A \code{Sincast} object with updated \code{pseudobulk} assay.
#'
#' @seealso [SincastImpute()]
#'
#' @export
#' @name SincastAggregate
#' @rdname SincastAggregate
#' @aliases Sincast, SincastAssays
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
                                        sep = "_", replace = FALSE, ...) {
  standardGeneric("SincastAggregate")
})

#' @name SincastAggregate
#' @rdname SincastAggregate
setMethod("SincastAggregate", "Sincast", function(object,
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
                                                  sep = "_", replace = FALSE, ...) {
  # Check the validity of the Sincast object.
  Sincast::CheckSincastObject(object, complete = FALSE, test = FALSE)

  # If a pseudobulk assay already exists, either replace it or return an error.
  if (!is.null(
    Sincast::GetSincastAssays(object, assay = "pseudobulk")
  )) {
    if (!replace) {
      stop(
        "SincastAggregate: A pseudobulk assay already exists, set replace = T to enforce a replacement."
      )
    } else {
      message("SincastAggregate: A pseudobulk assay already exists. Will be replaced as 'replace = T'.")
    }
  }

  # Get the original Seurat object.
  original <- Sincast::GetSincastAssays(object, assay = "original")

  # Get the default assay.
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(original)
  }else{
    Seurat::DefaultAssay(original) <- assay
  }

  # Get default features to aggregate.
  if (!is.null(features)) original <- original[features, ]

  # Get the default group.
  if (group.by == "ident") {
    group <- slot(original, name = "active.ident")
  } else {
    group <- slot(original, name = "meta.data")[, group.by]
  }
  group.name <- unique(group) %>% sort()

  # Check sample.method format.
  sample.method <- match.arg(sample.method)

  # Check aggregate.method format.
  aggregate.method <- match.arg(aggregate.method)

  # Store pseudobulk.
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
    tmp.data <- SeuratObject::GetAssayData(
      object = original,
      layer = layer
    )[, group == g, drop = F]

    # Generate pseudobulk.
    out[[g]] <- replicate(
      n = n.gen,
      expr = GeneratePseudobulk(tmp.data, n.pool, aggregate.method)
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

  # Calculate sparsity
  sparsity.before <- mean(
    SeuratObject::GetAssayData(
      object = original,
      layer = layer
    ) == 0
  )
  sparsity.after <- mean(out == 0)

  # Generate a token for the original data.
  if (is.null(
    Seurat::Misc(original, slot = "SincastToken")
  )) {
    SincastToken <- GenerateSincastToken()
    Seurat::Misc(
      Sincast::GetSincastAssays(object, assay = "original"),
      slot = "SincastToken"
    ) <- SincastToken
  }else{
    SincastToken <- Seurat::Misc(original, slot = "SincastToken")
  }

  # Generate a token for this aggregation run.

  SincastToken <- GenerateSincastToken(by = "SincastAggregate",
                                       command = deparse(match.call()),
                                       extend = SincastToken@command)

  # Create a new Seurat object to store the result.
  suppressWarnings(
    out <- Seurat::CreateSeuratObject(
      counts = Matrix::Matrix(out),
      assay = assay,
      meta.data = meta.data,
      project = SincastToken@id, ...
    )
  )
  Seurat::Idents(out) <- "agg.label" # Add aggregated annotation.

  # Write the summary information
  object@summary@summary["pseudobulk",] <-  c(assay, layer, nrow(out), ncol(out),
                                              NA, round(sparsity.before,3), round(sparsity.after,3))
  SincastToken@summary <- object@summary@summary["pseudobulk",]

  # Add token to the new Seurat object,
  Seurat::Misc(out, slot = "SincastToken") <- SincastToken

  # Add or replace the pseudobulk assay by the new Seurat object.
  Sincast::GetSincastAssays(object, assay = "pseudobulk") <- out

  message(
    "SincastAggregate: Sparsity before aggregation: ", round(sparsity.before, 3),
    "; After aggregation: ", round(sparsity.after, 3)
  )

  object
})
