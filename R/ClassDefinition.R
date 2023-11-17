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
#'
#' @family Sincast classes
#'
#' @name SincastToken-class
#' @rdname SincastToken-class
#' @aliases Sincast
SincastToken <- setClass(
  Class = "SincastToken",
  slots = list(
    id = "character",
    timestamp = "character",
    by = "character",
    command = 'list'
  )
)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SincastAssays
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NullSeurat <- setClassUnion(
  name = "NullSeurat",
  members = c("NULL", "Seurat")
)


#' A S4 class to store \code{Sincast} aggregation or imputation results.
#'
#' To be added.
#'
#' @slot pseudobulk A list of \code{Seurat} object storing aggregated pseudobulk.
#' @slot imputation A list of \code{Seurat} object storing imputed single cells.
#'
#' @family Sincast classes
#'
#' @name SincastAssays-class
#' @rdname SincastAssays-class
#' @aliases Sincast, SincastAssays
SincastAssays <- setClass(
  Class = "SincastAssays",
  slots = list(
    pseudobulk = "NullSeurat",
    imputation = "NullSeurat"
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
    test <- function(assay) {
      out <- "Valid"
      if (is.null(assay)) {
        out <- "Empty"
      } else {
        SincastToken <- Seurat::Misc(assay, slot = "SincastToken")
        if (!is(SincastToken, "SincastToken")) out <- "Sincast token is either missing or invalid."
      }
      out
    }

    test.result <- c(pseudobulk = "Valid", imputation = "Valid")
    test.result["pseudobulk"] <- test(object@pseudobulk)
    test.result["imputation"] <- test(object@imputation)


    test.result
  }
)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Sincast
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' An S4 class of \code{Sincast} object.
#'
#' To be added.
#'
#' @slot SincastAssays A S4 class to store \code{Sincast} aggregation or imputation results.
#' @slot SincastToken A S4 class to mark and time stamp the \code{Sincast} object.
#'
#' @family Sincast classes
#'
#' @name Sincast-class
#' @rdname Sincast-class
#' @aliases Sincast, SincastAssays
Sincast <- setClass(
  Class = "Sincast",
  slots = list(
    SincastAssays = "SincastAssays",
    SincastToken = "SincastToken"
  )
)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SincastSeurat
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' An S4 class of \code{SincastSeurat} object, extended based on \code{Seurat}.
#'
#' Inherit all slots from a \code{Seurat} object (See \code{\link[Seurat]{Seurat-class}}),
#' plus an additional slot \code{Sincast} storing a \code{Sincast}
#' object (See \code{\link[Sincast]{Sincast-class}}).
#'
#' @slot Sincast An S4 class of \code{Sincast} object.
#'
#' @family Sincast classes
#'
#' @seealso [as.SincastSeurat()] for converting a \code{Seurat} object to \code{Sincast}.
#'
#' @name SincastSeurat-class
#' @rdname SincastSeurat-class
#' @aliases Sincast, Seurat
SincastSeurat <- setClass(
  Class = "SincastSeurat",
  contains = "Seurat",
  slot = list(Sincast = 'Sincast')
)
