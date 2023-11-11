# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Sincast object, could be used
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Sincast <- setClass(
#   Class = "Sincast",
#   contains = "Seurat"
# )
#
# setAs("Seurat", "Sincast", function(from, to) {
#   arguments <- paste(slotNames(from), "=from@", slotNames(from), sep = "", collapse = ",")
#   text2expr <- paste("Sincast(", arguments, ")", sep = "")
#   eval(parse(text = text2expr))
# })

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 class definition
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' An S4 class to store \code{Sincast} aggregation or imputation results.
#'
#' To be added.
#'
#' @slot pseudobulk A list of \code{Seurat} object storing aggregated pseudobulk.
#' @slot imputation A list of \code{Seurat} object storing imputed single cells.
#'
#' @family SincastAssays related methods
#'
#' @name SincastAssays-class
#' @rdname SincastAssays-class
#' @aliases Sincast, SincastAssays
SincastAssays <- setClass(
  Class = "SincastAssays",
  slots = list(
    pseudobulk = "list",
    imputation = "list"
  )
)

#' SincastAssays Object Validity
#'
#' Validation of \code{SincastAssays} objects is handled by \code{\link[methods]{validObject}}.
#'
#' @name SincastAssays-validity
#' @rdname SincastAssays-class
#' @aliases Sincast, SincastAssays
setValidity(
  Class = "SincastAssays",
  method = function(object) {
    if (length(slot(object, "pseudobulk"))) {
      test <- sapply(slot(object, "pseudobulk"), function(x) is(x, "Seurat"))
      if (!all(test)) {
        return("'pseudobulk' must be a list of Seurat objects.")
      }
    }

    if (length(slot(object, "imputation"))) {
      test <- sapply(slot(object, "imputation"), function(x) is(x, "Seurat"))
      if (!all(test)) {
        return("'imputation' must be a list of Seurat objects.")
      }
    }
  }
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Add a \code{SincastAssay} object to Seurat's \code{misc} slots
#'
#' To be added.
#'
#' @param object A \code{Seurat} object
#'
#' @return An updated \code{Seurat} object with an additional \code{SincastAssays} object in the \code{misc} slot.
#'
#' @family SincastAssays related methods
#'
#' @export
#' @rdname AddSincast
setGeneric("AddSincast", function(object, ...) standardGeneric("AddSincast"))

#' @rdname AddSincast
#' @aliases Sincast, SincastAssays, Seurat
setMethod("AddSincast", "Seurat", function(object) {
  if (is.null(Seurat::Misc(object, slot = "SincastAssays"))) {
    Seurat::Misc(object, slot = "SincastAssays") <- new("SincastAssays")
  } else if (is(Seurat::Misc(object, slot = "SincastAssays"), "SincastAssays")) {
    message("SincastAssays object already exists, check if it's valid.")

    if (
      validObject(
        Seurat::Misc(object, slot = "SincastAssays"),
        test = T
      ) != TRUE
    ) {
      message("Replace the bad Assay.")
      Seurat::Misc(object, slot = "SincastAssays") <- new("SincastAssays")
    }
  } else {
    warnings("Replacing an unknown SincastAssays slot.")
    Seurat::Misc(object, slot = "SincastAssays") <- new("SincastAssays")
  }

  object
})
