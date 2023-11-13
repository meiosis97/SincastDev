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
  method = function(SincastAssays) {
    test <- function(i, assay) {
      # Check whether a element is a "Seurat" object.
      out <- "valid"
      if (!is(assay[[i]], "Seurat")) out <- "Not a seurat object"

      # Check whether the Seurat object contains a 'SincastToken' object.
      if (out == "valid") {
        SincastToken <- Seurat::Misc(assay[[i]], slot = "SincastToken")
        if (!is(SincastToken, "SincastToken")) out <- "Sincast token is either missing or invalid"
      }
      # Check whether the 'SincastToken' object matches with the name of the element.
      if (out == "valid") {
        if (is.null(names(assay[i]))) {
          out <- "Object is not named"
        } else if (names(assay[i]) != SincastToken@id) {
          out <- "Sincast token and the name of the object do not match"
        }
      }
      out
    }

    test.result <- list(pseudobulk = NULL, imputation = NULL)
    n.pseudobulk <- length(SincastAssays@pseudobulk)
    if (n.pseudobulk) {
      test.result$pseudobulk <- sapply(1:n.pseudobulk,
        test,
        assay = SincastAssays@pseudobulk
      )
    }

    n.imputation <- length(SincastAssays@imputation)
    if (n.imputation) {
      test.result$imputation <- sapply(1:n.imputation,
        test,
        assay = SincastAssays@imputation
      )
    }

    test.result
  }
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CreateSincastAssays
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Add a \code{SincastAssays} object to \code{Seurat}'s \code{misc} slot.
#'
#' To be added.
#'
#' @param object A \code{Seurat} object.
#' @param remove.invalid Logical. Whether to check the validity of the existing \code{SincastAssays} object, and
#' remove the invalid elements (e.g, non-\code{Seurat} object) identified in \code{Sincast} assays.
#'
#' @return An updated \code{Seurat} object with an additional \code{SincastAssays} object in the \code{misc} slot.
#'
#' @family SincastAssays related methods
#'
#' @export
#' @rdname CreateSincastAssays
#' @aliases Sincast, SincastAssays, Seurat
setGeneric("CreateSincastAssays", function(object, remove.invalid = TRUE, ...) standardGeneric("CreateSincastAssays"))

#' @rdname CreateSincastAssays
setMethod("CreateSincastAssays", "Seurat", function(object, remove.invalid = TRUE) {
  if (is.null(Seurat::Misc(object, slot = "SincastAssays"))) {
    Seurat::Misc(object, slot = "SincastAssays") <- new("SincastAssays")
    Seurat::Misc(object, slot = "SincastToken") <- GenerateSincastToken()
  } else {
    message("CreateSincastAssays: 'SincastAssays' object already exists. Check if it's valid.")
    object <- CleanSincastAssays(object, remove.invalid = remove.invalid)

  }

  object
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CleanSincastAssays
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Remove unwanted or invalid elements (i.e, non-\code{Seurat} object) in \code{Sincast} assays
#' \code{pseudobulk} and \code{imputation} embedded in \code{Seurat::MISC(object, slot = "SincastAssays")}.
#'
#' To be added.
#'
#' @param object A \code{Seurat} object with an additional \code{SincastAssays} object in the \code{misc} slot.
#' @param remove.id.pseudobulk A vector of ids of the elements to be removed from the \code{pseduobulk} assay.
#' @param remove.id.imputation A vector of ids of the elements to be removed from the \code{Imputation} assay.
#' @param clean.up Either "All", "pseudobulk" or "imputation" indicating whether to clean up
#' the \code{SincastAssays} object or clean up a specific \code{Sincast} assay as proposed.
#' @param remove.invalid Logical. Whether to check the validity of the existing \code{SincastAssays} object, and
#' remove the invalid elements (e.g, non-\code{Seurat} object) identified in \code{Sincast} assays.
#'
#' @return An updated \code{Seurat} object with an additional \code{SincastAssays} object in the \code{misc} slot.
#'
#' @family SincastAssays related methods
#'
#' @export
#' @rdname CleanSincastAssays
#' @aliases Sincast, SincastAssays, Seurat
setGeneric("CleanSincastAssays", function(object,
                                          remove.id.pseudobulk = NULL,
                                          remove.id.imputation = NULL,
                                          clean.up = c("None", "All", "pseudobulk", "imputation"),
                                          remove.invalid = TRUE,
                                          ...) {
  standardGeneric("CleanSincastAssays")
})

#' @rdname CleanSincastAssays
setMethod("CleanSincastAssays", "Seurat", function(object,
                                                   remove.id.pseudobulk = NULL,
                                                   remove.id.imputation = NULL,
                                                   clean.up = c("None", "All", "pseudobulk", "imputation"),
                                                   remove.invalid = TRUE) {
  SincastAssays <- Seurat::Misc(object, slot = "SincastAssays")
  ret.summary <- FALSE

  if (is.null(SincastAssays)) {
    if (remove.invalid) {
      message("CleanSincastAssays: 'SincastAssays' object does not exists. Add a 'SincastAssays' object.")
      SincastAssays <- new("SincastAssays")
    } else {
      message("CleanSincastAssays: 'SincastAssays' object does not exists. Add it by setting 'remove.invalid = TRUE'.")
    }
  } else if (is(SincastAssays, "SincastAssays")) {
    # Remove elements in the pseudobulk assay.
    if (!is.null(remove.id.pseudobulk)) {
      remove.id.pseudobulk <- intersect(
        remove.id.pseudobulk,
        1:length(SincastAssays@pseudobulk)
      )
      if (length(remove.id.pseudobulk)) {
        message(
          "CleanSincastAssays: Remove elements",
          remove.id.pseudobulk %>% paste(collapse = ","), " in the 'pseudobulk' assay."
        )
        SincastAssays@pseudobulk <-
          SincastAssays@pseudobulk[-remove.id.pseudobulk]
      }
    }

    # Remove elements in the imputation assay.
    if (!is.null(remove.id.imputation)) {
      remove.id.imputation <- intersect(
        remove.id.imputation,
        1:length(SincastAssays@imputation)
      )
      if (length(remove.id.imputation)) {
        message(
          "CleanSincastAssays: Remove elements",
          remove.id.imputation %>% paste(collapse = ","), " in the 'imputation' assay."
        )
        SincastAssays@imputation <-
          SincastAssays@imputation[-remove.id.imputation]
      }
    }

    # Clean up a assay (assays)
    if (!is.null(clean.up)) {
      # Check the validity of clean.up
      clean.up <- match.arg(clean.up)
      if (clean.up == "All") {
        SincastAssays <- new("SincastAssays")
      } else if (clean.up == "pseudobulk") {
        SincastAssays@pseudobulk <- list()
      } else if (clean.up == "imputation") {
        SincastAssays@imputation <- list()
      }
    }

    # Check and remove invalid elements in Sincast assays.
    test.result <- validObject(
      SincastAssays,
      test = T
    )
    all.valid <- TRUE
    invalid.pseudobulk <- which(test.result$pseudobulk != "valid")
    invalid.imputation <- which(test.result$imputation != "valid")
    n.valid.pseudobulk <- sum(test.result$pseudobulk == "valid")
    n.valid.imputation <- sum(test.result$imputation == "valid")

    if (length(invalid.pseudobulk)) {
      if(remove.invalid){
        message(
          "CleanSincastAssays: Remove invalid elements ",
          invalid.pseudobulk %>% paste(collapse = ","), " in the 'pseudobulk' assay."
        )
        SincastAssays@pseudobulk <-
          SincastAssays@pseudobulk[-invalid.pseudobulk]
      }else{
        message(
          "CleanSincastAssays: ", invalid.pseudobulk %>% paste(collapse = ","),
          " in the 'pseudobulk' assay are invalid. Set 'remove.invalid = TRUE' to remove."
        )
      }
      all.valid <- FALSE
    }

    if (length(invalid.imputation)) {
      if(remove.invalid){
        message(
          "CleanSincastAssays: Remove invalid elements ",
          invalid.imputation %>% paste(collapse = ","), " in the 'imputation' assay."
        )
        SincastAssays@imputation <-
          SincastAssays@imputation[-invalid.imputation]
      }else{
        message(
          "CleanSincastAssays: ", invalid.imputation %>% paste(collapse = ","),
          " in the 'imputation' assay are invalid. Set 'remove.invalid = TRUE' to remove."
        )
      }
      all.valid <- FALSE
    }
    if (all.valid) {
      message("CleanSincastAssays: Valid.")
    }
    ret.summary <- TRUE

  } else {
    if (remove.invalid) {
      message("CleanSincastAssays: Replacing an unknown 'SincastAssays' object.")
      SincastAssays <- new("SincastAssays")
    } else {
      message("CleanSincastAssays: There is an unknown ’SincastAssays‘ object. Remove it by setting 'remove.invalid = TRUE'.")
    }
  }

  # Update
  suppressWarnings(Seurat::Misc(object, slot = "SincastAssays") <- SincastAssays)
  SincastToken <- Seurat::Misc(object, slot = "SincastToken")
  if (is.null(SincastToken)) {
    SincastToken <- GenerateSincastToken()
  } else if (!is(SincastToken, "SincastToken")) {
    SincastToken <- GenerateSincastToken()
  }
  suppressWarnings(Seurat::Misc(object, slot = "SincastToken") <- SincastToken)

  if (ret.summary) {
    message(
      "CleanSincastAssays: 'SincastAssays' object contains ", n.valid.pseudobulk,
      " valid pseudobulk objects and ", n.valid.imputation, " valid imputation objects."
    )
  }

  object
})
