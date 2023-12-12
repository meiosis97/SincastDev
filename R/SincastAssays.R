#' Check the validity of the \code{SincastAssays} object of \code{Sincast}.
#'
#' To be added.
#'
#' @param object A \code{SincastAssays} object.
#' @param test Logical; if TRUE (the default) and validity fails, the function returns a
#'  vector of strings describing the problems. If test is FALSE validity failure generates an error.
#' @param silent Logical; if TRUE, suppress all messages.
#'
#' @return A vector of strings describing the problems.
#'
#' @family SincastAssays related methods
#'
#' @export
#' @name CheckSincastAssays
#' @rdname CheckSincastAssays
setGeneric("CheckSincastAssays", function(object, test = TRUE, slient = FALSE, ...) {
  standardGeneric("CheckSincastAssays")
})

#' @rdname CheckSincastAssays
setMethod("CheckSincastAssays", "SincastAssays", function(object, test = TRUE, silent = FALSE, ...) {
  # A single test.
  check <- function(x) {
    out <- "Valid"
    if (is.null(x)) {
      out <- "Empty"
    } else if (class(x) != "Seurat") {
      out <- "Not a 'Seurat' object."
    } else {
      SincastToken <- Seurat::Misc(x, slot = "SincastToken")
      if (!is(SincastToken, "SincastToken")) out <- "Sincast token is either missing or invalid"
    }
    out
  }

  test.result <- c(pseudobulk = "Valid", imputation = "Valid")
  test.result["original"] <- check(object@original)
  test.result["pseudobulk"] <- check(object@pseudobulk)
  test.result["imputation"] <- check(object@imputation)

  # Create an error message.
  if (!all(test.result %in% c("Valid", "Empty"))) {
    problem <- paste(
      "CheckSincastAssays: ",
      "\n original Seurat: ", test.result["original"],
      "\n pseudobulk assay: ", test.result["pseudobulk"],
      "\n imputation assay: ", test.result["imputation"],
      collapse = ""
    )
  } else {
    problem <- NULL
  }

  # Print the error messages
  if (!is.null(problem)) {
    if (!test) stop(problem) else if (!silent) message(problem)
  }

  test.result
})

#' @rdname CheckSincastAssays
setMethod("CheckSincastAssays", "Sincast", function(object, test = TRUE, silent = FALSE, ...) {
  # Check the validity of the Sincast object.
  Sincast::CheckSincastObject(object, complete = FALSE, test = FALSE)

  object <- Sincast::GetSincastAssays(object)
  Sincast::CheckSincastAssays(
    object = object, test = test, silent = silent, ...  )
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GetSincastAssays
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Extract the \code{SincastAssays} object from a \code{Sincast} object.
#'
#' To be added.
#'
#' @param object A \code{Sincast} object.
#' @param assay Either \code{"all"}, \code{"original"}, \code{"pseudobulk"} or
#' \code{"imputation"} indicating which specific \code{Sincast} assay to extract.
#'
#' @return Depending on the \code{assay} argument, can be either type of a \code{Sicast} assay,
#' or a \code{SincastAssays} object containing all the types.
#'
#' @family SincastAssays related methods
#'
#' @export
#' @name GetSincastAssays
#' @rdname GetSincastAssays
#' @aliases Sincast, SincastAssays
setGeneric("GetSincastAssays", function(object, assay = c(
                                          "all", "original",
                                          "pseudobulk", "imputation"
                                        ), ...) {
  standardGeneric("GetSincastAssays")
})

#' @rdname GetSincastAssays
setMethod("GetSincastAssays", "Sincast", function(object,
                                                  assay = c(
                                                    "all", "original",
                                                    "pseudobulk", "imputation"
                                                  ), ...) {
  # Check the validity of the Sincast object.
  test.SincastObject <- Sincast::CheckSincastObject(object, complete = FALSE)

  # If the Sincast object is missing, or invalid, return a NULL
  if (test.SincastObject != "Valid") {
    out <- NULL
  } else {
    assay <- match.arg(assay)

    # Get all assays.
    out <- object@SincastAssays

    # Get the original Seurat.
    if (assay == "original") out <- out@original

    # Get the pseudobulk assay.
    if (assay == "pseudobulk") out <- out@pseudobulk

    # Get the imputation assay.
    if (assay == "imputation") out <- out@imputation
  }

  out
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GetSincastAssays<-
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Setter function for \code{GetSincastAssays}.
#'
#' To be added.
#'
#' @param object A \code{Sincast} object.
#' @param assay Either \code{"all"}, \code{"original"}, \code{"pseudobulk"} or
#' \code{"imputation"} indicating which specific \code{Sincast} assay to modify.
#'
#' @return A \code{Sincast} object with updated \code{Sincast} assays.
#'
#' @family SincastAssays related methods
#'
#' @export
#' @name GetSincastAssays
#' @rdname GetSincastAssays
setGeneric("GetSincastAssays<-", function(object, assay = c(
                                            "all", "original",
                                            "pseudobulk", "imputation"
                                          ), value, ...) {
  standardGeneric("GetSincastAssays<-")
})

#' @rdname GetSincastAssays
setMethod("GetSincastAssays<-", "Sincast", function(object, assay = c(
                                                      "all", "original",
                                                      "pseudobulk", "imputation"
                                                    ), value, ...) {
  # Check the validity of the Sincast object.
  Sincast::CheckSincastObject(object, complete = FALSE, test = FALSE)

  assay <- match.arg(assay)

  # Assign all assays.
  if (assay == "both") object@SincastAssays <- value

  # Assign the original Seurat.
  if (assay == "original") object@SincastAssays@original <- value

  # Assign the pseudobulk assay.
  if (assay == "pseudobulk") object@SincastAssays@pseudobulk <- value

  # Assign the imputation assay.
  if (assay == "imputation") object@SincastAssays@imputation <- value


  # Check the validity of the Sincast assays.
  test.SincastAssays <- Sincast::CheckSincastAssays(object@SincastAssays, silent = TRUE)

  if (assay == "all" & !all(test.SincastAssays %in% c("Valid", "Empty"))) {
    warning("At least one Sincast Assay is invalid.")
  }

  if(assay != "all") {
    if(test.SincastAssays[assay] != "Valid"){
      warning("Invalid ", assay, " assay provided.")
    }
  }

  object
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CleanSincastAssays
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Remove unwanted or invalid (i.e, non-\code{Seurat} object) in \code{Sincast} assays
#' (\code{pseudobulk} and \code{imputation}) embedded in \code{Sincast}.
#'
#' To be added.
#'
#' @param object A \code{Sincast} object.
#' @param clean.up Either \code{"None"}, \code{"All"}, \code{"pseudobulk"} or \code{"imputation"} indicating whether to clean up
#' the \code{SincastAssays} object or clean up a specific \code{Sincast} assay as proposed.
#' @param remove.invalid Logical; if TRUE, check the validity of the existing \code{Sincast} assays and
#' remove invalid ones (e.g, \code{Seurat} objects not generated by \code{Sincast} functions).
#'
#' @return A \code{Sincast} object with updated \code{Sincast} assays.
#'
#' @family SincastAssays related methods
#'
#' @export
#' @name CleanSincastAssays
#' @rdname CleanSincastAssays
setGeneric("CleanSincastAssays", function(object,
                                          clean.up = c("none", "all", "pseudobulk", "imputation"),
                                          remove.invalid = TRUE,
                                          ...) {
  standardGeneric("CleanSincastAssays")
})

#' @rdname CleanSincastAssays
setMethod("CleanSincastAssays", "SincastAssays", function(object,
                                                    clean.up = c("none", "all", "pseudobulk", "imputation"),
                                                    remove.invalid = TRUE, ...) {
  test.SincastAssays <- Sincast::CheckSincastAssays(object, silent = TRUE)
  message(
    "CleanSincastAssays: Before clean up: ",
    "\n pseudobulk assay: ", test.SincastAssays["pseudobulk"],
    "\n imputation assay: ", test.SincastAssays["imputation"]
  )

  # Clean up a assay (assays)
  if (!is.null(clean.up)) {
    # Check the validity of clean.up
    clean.up <- match.arg(clean.up)
    if (clean.up == "all") {
      object <- new("SincastAssays")
    } else if (clean.up == "pseudobulk") {
      object@pseudobulk <- NULL
    } else if (clean.up == "imputation") {
      object@imputation <- NULL
    }
  }

  # Check and remove invalid elements in Sincast assays.
  test.SincastAssays <- Sincast::CheckSincastAssays(object, silent = TRUE)
  is.pseudobulk.valid <- test.SincastAssays["pseudobulk"] %in% c("Valid", "Empty")
  is.imputation.valid <- test.SincastAssays["imputation"] %in% c("Valid", "Empty")

  # Clean up the pseudobulk assay
  if (!is.pseudobulk.valid) {
    if (remove.invalid) {
      message(
        "CleanSincastAssays: Remove invalid 'pseudobulk' assay  as 'remove.invalid = TRUE'."
      )
      object@pseudobulk <- NULL
      test.SincastAssays["pseudobulk"] <- "Empty"
    } else {
      message(
        "CleanSincastAssays: ", "the 'pseudobulk' assay is invalid. Consider replace it by setting 'remove.invalid = TRUE'."
      )
    }
  }

  # Clean up the imputation assay
  if (!is.imputation.valid) {
    if (remove.invalid) {
      message(
        "CleanSincastAssays: Remove invalid 'imputation' assay as 'remove.invalid = TRUE'."
      )
      object@imputation <- NULL
      test.SincastAssays["imputation"] <- "Empty"
    } else {
      message(
        "CleanSincastAssays: ", "the 'imputation' assay is invalid. Consider replace it by setting 'remove.invalid = TRUE'."
      )
    }
  }

  message(
    "CleanSincastAssays: After clean up: ",
    "\n pseudobulk assay: ", test.SincastAssays["pseudobulk"],
    "\n imputation assay: ", test.SincastAssays["imputation"]
  )

  object
})

#' @rdname CleanSincastAssays
setMethod("CleanSincastAssays", "Sincast", function(object,
                                                          clean.up = c("none", "all", "pseudobulk", "imputation"),
                                                          remove.invalid = TRUE, ...) {
  # Check the validity of the Sincast object.
  Sincast::CheckSincastObject(object, complete = FALSE, test = FALSE)

  SincastAssays <- Sincast::GetSincastAssays(object)
  SincastAssays <- Sincast::CleanSincastAssays(object = SincastAssays,
                              clean.up = clean.up,
                              remove.invalid = remove.invalid, ...)
  Sincast::GetSincastAssays(object) <- SincastAssays
  object
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# show.SincastAssays
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setMethod("show", "SincastAssays", function(object) {
  # Convert Seurat show
  out <- capture.output(object@original)
  out[1] <- "An object of class SincastAssays"
  out <- paste(out, collapse = "\n")

  test.SincastAssays <- Sincast::CheckSincastAssays(object, silent = TRUE)

  out <- paste(out, "\nSincast assay:")
  if(test.SincastAssays["pseudobulk"] == "Empty" &
     test.SincastAssays["imputation"] == "Empty"){
    out <- paste(out, "\n", " No other Sincast assay present.", sep = "")
  }

  # Print information for Sincast's pseudobulk assay.
  if (test.SincastAssays["pseudobulk"] == "Valid") {
    summary <- Seurat::Misc(object@pseudobulk, "SincastToken")@summary
    out <- paste(out, "\n", " pseudobulk assay: (",
      "Computed using ", summary[,"assay"], "-", summary[,"layer"], ", ",
      summary[,"nfeatures"], " features, ",
      summary[,"nsamples"], " Samples, ",
      summary[,"sparsity.before"], " sparsity before, ",
      summary[,"sparsity.after"], " sparsity after)",
      sep = ""
    )
  } else if(test.SincastAssays["pseudobulk"] != "Empty"){
    out <- paste(out, "\n", " pseudobulk assay: (", test.SincastAssays["pseudobulk"], ")", sep = "")
  }

  # Print information for Sincast's imputation assay.
  if (test.SincastAssays["imputation"] == "Valid") {
    summary <- Seurat::Misc(object@imputation, "SincastToken")@summary
    out <- paste(out, "\n", " imputation assay: (",
      "Computed using ", summary[,"assay"], "-", summary[,"layer"], ", ",
      summary[,"nfeatures"], " features, ",
      summary[,"nsamples"], " Samples, ",
      summary[,"sparsity.before"], " sparsity before, ",
      summary[,"sparsity.after"], " sparsity after)",
      sep = ""
    )
  } else if(test.SincastAssays["imputation"] != "Empty"){
    out <- paste(out, "\n", " imputation assay: (", test.SincastAssays["imputation"], ")", sep = "")
  }

  cat(out)
})
