# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SincastToken
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' An S4 class of \code{SincastToken} object.
#'
#' To be added.
#'
#' @slot id A 16 character unique identifier of a particular \code{Sincast} result.
#' @slot timestamp A string recording the time at which the result was generated.
#' @slot by A string recording the function by which this \code{Sincast} result was generated.
#' @slot command A list recording \code{Sincast} command history.
#' @slot summary A one row \code{data.frame} storing summary information.
#'
#' @family Sincast classes
#'
#' @name SincastToken-class
#' @rdname SincastToken-class
SincastToken <- setClass(
  Class = "SincastToken",
  slots = list(
    id = "character",
    timestamp = "POSIXct",
    by = "character",
    command = "list",
    summary = "data.frame"
  )
)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SincastSummary
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' An S4 class of \code{SincastSummary} object.
#'
#' To be added.
#'
#' @slot summary A \code{data.frame} summarize \code{Sincast} results.
#'
#' @family Sincast classes
#'
#' @name SincastSummary-class
#' @rdname SincastSummary-class
SincastSummary <- setClass(
  Class = "SincastSummary",
  slots = list(
    summary = "data.frame",
    active.assay = "character",
    active.atlas = "character"
  )
)

#' SincastSummary Object Validity
#'
#' Validation of \code{SincastSummary} objects is handled by \code{\link[methods]{validObject}}.
#'
#' @name SincastSummary-validity
#' @rdname SincastSummary-class
setValidity(
  Class = "SincastSummary",
  method = function(object) {
    out <- NULL

    if (any(rownames(object@summary) != c(
      "pseudobulk", "imputation",
      "original.atlas", "pseudobulk.atlas",
      "imputation.atlas"
    ))) {
      out <- c(out, "SincastSummary table has unknow row names.")
    }

    if (any(colnames(object@summary) != c(
      "assay", "layer", "nfeatures",
      "nsamples", "n.components", "var.explained",
      "sparsity.before", "sparsity.after"
    ))) {
      out <- c(out, "SincastSummary table has unknow column names.")
    }


    if (
      !IsCharacterVector(
        object@active.assay,
        n = 1, in.which = c("original", "imputation", "pseudobulk")
      )) {
      out <- c(out, "active.assay should be either one of 'original', 'imputation' or 'pseudobulk'.")
    }

    if (
      !IsCharacterVector(
        object@active.atlas,
        n = c(0,1), in.which = c("original", "imputation", "pseudobulk")
      )) {
      out <- c(out, "active.atlas should be either one of 'original', 'imputation' or 'pseudobulk'.")
    }

    if (is.null(out)) out <- TRUE

    out
  }
)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SincastAssays
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' A S4 class to store \code{Sincast} aggregation or imputation results.
#'
#' To be added.
#'
#' @slot original A \code{Seurat} object storing the original single cell data.
#' @slot pseudobulk A \code{Seurat} object storing aggregated pseudobulk.
#' @slot imputation A \code{Seurat} object storing imputed single cells.
#'
#' @family Sincast classes
#'
#' @name SincastAssays-class
#' @rdname SincastAssays-class
SincastAssays <- setClass(
  Class = "SincastAssays",
  slots = list(
    original = "Seurat",
    pseudobulk = "ANY",
    imputation = "ANY"
  )
)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SincastAtlas
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' A S4 class to store \code{Sincast} pseudobulk or imputation atlas
#'
#' To be added.
#'
#' @slot original A \code{Seurat} object storing the atlas constructed on the original data.
#' @slot pseudobulk A \code{Seurat} object storing the pseudobulk atlas.
#' @slot imputation A \code{Seurat} object storing the imputation atlas.
#'
#' @family Sincast classes
#'
#' @name SincastAtlas-class
#' @rdname SincastAtlas-class
SincastAtlas <- setClass(
  Class = "SincastAtlas",
  slots = list(
    original = "ANY",
    pseudobulk = "ANY",
    imputation = "ANY"
  )
)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Sincast
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' An S4 class of \code{Sincast} object.
#'
#' To be added.
#'
#' @slot SincastAssays A S4 class to store \code{Sincast} aggregation or imputation results.
#' @slot SincastAtlas A S4 class to store \code{Sincast} pseudobulk or imputation atlas.
#' @slot Summary A \code{SincastSummary} object.
#'
#' @family Sincast classes
#'
#' @name Sincast-class
#' @rdname Sincast-class
Sincast <- setClass(
  Class = "Sincast",
  slots = list(
    SincastAssays = "SincastAssays",
    SincastAtlas = "SincastAtlas",
    summary = "SincastSummary"
  )
)
