#' Check the validity of the \code{SincastAtlas} object of \code{Sincast}.
#'
#' To be added.
#'
#' @param object A \code{SincastAtlas} object.
#' @param test Logical; if TRUE (the default) and validity fails, the function returns a
#'  vector of strings describing the problems. If test is FALSE validity failure generates an error.
#' @param silent Logical; if TRUE, suppress all messages.
#'
#' @return A vector of strings describing the problems.
#'
#' @family SincastAtlas related methods
#'
#' @export
#' @name CheckSincastAtlas
#' @rdname CheckSincastAtlas
#' @aliases Sincast, SincastAtlas
setGeneric("CheckSincastAtlas", function(object, test = TRUE, slient = FALSE, ...) {
  standardGeneric("CheckSincastAtlas")
})

#' @rdname CheckSincastAtlas
setMethod("CheckSincastAtlas", "SincastAtlas", function(object, test = TRUE, silent = FALSE, ...) {
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
      "CheckSincastAtlas: ",
      "\n original atlas: ", test.result["original"],
      "\n pseudobulk atlas: ", test.result["pseudobulk"],
      "\n imputation atlas: ", test.result["imputation"],
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

#' @rdname CheckSincastAtlas
setMethod("CheckSincastAtlas", "Sincast", function(object, test = TRUE, silent = FALSE, ...) {
  # Check the validity of the Sincast object.
  Sincast::CheckSincastObject(object, complete = FALSE, test = FALSE)

  object <- Sincast::GetSincastAtlas(object)
  Sincast::CheckSincastAtlas(
    object = object, test = test, silent = silent, ...  )
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GetSincastAtlas
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Extract the \code{SincastAtlas} object from a \code{Sincast} object.
#'
#' To be added.
#'
#' @param object A \code{Sincast} object.
#' @param atlas Either \code{"all"}, \code{"original"}, \code{"pseudobulk"} or
#' \code{"imputation"} indicating which specific \code{Sincast} atlas to extract.
#'
#' @return Depending on the \code{atlas} argument, can be either type of a \code{Sicast} atlas,
#' or a \code{SincastAtlas} object containing all the types.
#'
#' @family SincastAtlas related methods
#'
#' @export
#' @name GetSincastAtlas
#' @rdname GetSincastAtlas
#' @aliases Sincast, SincastAtlas
setGeneric("GetSincastAtlas", function(object, atlas = c(
  "all", "original",
  "pseudobulk", "imputation"
), ...) {
  standardGeneric("GetSincastAtlas")
})

#' @rdname GetSincastAtlas
setMethod("GetSincastAtlas", "Sincast", function(object,
                                                  atlas = c(
                                                    "all", "original",
                                                    "pseudobulk", "imputation"
                                                  ), ...) {
  # Check the validity of the Sincast object.
  test.SincastObject <- Sincast::CheckSincastObject(object, complete = FALSE)

  # If the Sincast object is missing, or invalid, return a NULL
  if (test.SincastObject != "Valid") {
    out <- NULL
  } else {
    atlas <- match.arg(atlas)

    # Get all atlas.
    out <- object@SincastAtlas

    # Get the original Seurat.
    if (atlas == "original") out <- out@original

    # Get the pseudobulk atlas.
    if (atlas == "pseudobulk") out <- out@pseudobulk

    # Get the imputation atlas.
    if (atlas == "imputation") out <- out@imputation
  }

  out
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GetSincastAtlas<-
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Setter function for \code{GetSincastAtlas}.
#'
#' To be added.
#'
#' @param object A \code{Sincast} object.
#' @param atlas Either \code{"all"}, \code{"original"}, \code{"pseudobulk"} or
#' \code{"imputation"} indicating which specific \code{Sincast} atlas to modify.
#'
#' @return A \code{Sincast} object with updated \code{Sincast} atlas.
#'
#' @family SincastAtlas related methods
#'
#' @export
#' @name GetSincastAtlas
#' @rdname GetSincastAtlas
#' @aliases Sincast, SincastAtlas
setGeneric("GetSincastAtlas<-", function(object, atlas = c(
  "all", "original",
  "pseudobulk", "imputation"
), value, ...) {
  standardGeneric("GetSincastAtlas<-")
})

#' @rdname GetSincastAtlas
setMethod("GetSincastAtlas<-", "Sincast", function(object, atlas = c(
  "all", "original",
  "pseudobulk", "imputation"
), value, ...) {
  # Check the validity of the Sincast object.
  Sincast::CheckSincastObject(object, complete = FALSE, test = FALSE)

  atlas <- match.arg(atlas)

  # Assign all atlas.
  if (atlas == "all") object@SincastAtlas <- value

  # Assign the original Seurat.
  if (atlas == "original") object@SincastAtlas@original <- value

  # Assign the pseudobulk atlas.
  if (atlas == "pseudobulk") object@SincastAtlas@pseudobulk <- value

  # Assign the imputation atlas.
  if (atlas == "imputation") object@SincastAtlas@imputation <- value


  # Check the validity of the Sincast atlas.
  test.SincastAtlas <- Sincast::CheckSincastAtlas(object@SincastAtlas, silent = TRUE)

  if (atlas == "all" & !all(test.SincastAtlas %in% c("Valid", "Empty"))) {
    warning("At least one Sincast atlas is invalid.")
  }

  if (atlas != "all") {
    if (test.SincastAtlas[atlas] != "Valid"){
      warning("Invalid ", atlas, " atlas provided.")
    }
  }

  object
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CleanSincastAtlas
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Remove unwanted or invalid (i.e, non-\code{Seurat} object) \code{Sincast} atlas
#' (\code{pseudobulk} and \code{imputation}) embedded in \code{Sincast}.
#'
#' To be added.
#'
#' @param object A \code{Sincast} object.
#' @param clean.up Either \code{"none"}, \code{"all"}, \code{original},
#' \code{"pseudobulk"} or \code{"imputation"} indicating whether to clean up
#' the \code{SincastAtlas} object or clean up a specific \code{Sincast} atlas as proposed.
#' @param remove.invalid Logical; if TRUE, check the validity of the existing \code{Sincast} atlas and
#' remove invalid ones (e.g, \code{Seurat} objects not generated by \code{Sincast} functions).
#'
#' @return A \code{Sincast} object with updated \code{Sincast} atlas
#'
#' @family SincastAtlas related methods
#'
#' @export
#' @name CleanSincastAtlas
#' @rdname CleanSincastAtlas
#' @aliases Sincast, SincastAtlas
setGeneric("CleanSincastAtlas", function(object,
                                          clean.up = c("none", "all", "original", "pseudobulk", "imputation"),
                                          remove.invalid = TRUE,
                                          ...) {
  standardGeneric("CleanSincastAtlas")
})

#' @rdname CleanSincastAtlas
setMethod("CleanSincastAtlas", "SincastAtlas", function(object,
                                                          clean.up = c("none", "all", "original", "pseudobulk", "imputation"),
                                                          remove.invalid = TRUE, ...) {
  test.SincastAtlas <- Sincast::CheckSincastAtlas(object, silent = TRUE)
  message(
    "CleanSincastAtlas: Before clean up: ",
    "\n original atlas: ", test.SincastAtlas["original"],
    "\n pseudobulk atlas: ", test.SincastAtlas["pseudobulk"],
    "\n imputation atlas: ", test.SincastAtlas["imputation"]
  )

  # Clean up a atlas
  if (!is.null(clean.up)) {
    # Check the validity of clean.up
    clean.up <- match.arg(clean.up)
    if (clean.up == "all") {
      object <- new("SincastAtlas")
    } else if (clean.up == "original") {
      object@original <- NULL
    } else if (clean.up == "pseudobulk") {
      object@pseudobulk <- NULL
    } else if (clean.up == "imputation") {
      object@imputation <- NULL
    }
  }

  # Check and remove invalid elements in Sincast atlas.
  test.SincastAtlas <- Sincast::CheckSincastAtlas(object, silent = TRUE)
  is.original.valid <- test.SincastAtlas["original"] %in% c("Valid", "Empty")
  is.pseudobulk.valid <- test.SincastAtlas["pseudobulk"] %in% c("Valid", "Empty")
  is.imputation.valid <- test.SincastAtlas["imputation"] %in% c("Valid", "Empty")

  # Clean up the original atlas
  if (!is.original.valid) {
    if (remove.invalid) {
      message(
        "CleanSincastAtlas: Remove invalid original atlas  as 'remove.invalid = TRUE'."
      )
      object@original <- NULL
      test.SincastAtlas["original"] <- "Empty"
    } else {
      message(
        "CleanSincastAtlas: ", "the original atlas is invalid. Consider replace it by setting 'remove.invalid = TRUE'."
      )
    }
  }

  # Clean up the pseudobulk atlas
  if (!is.pseudobulk.valid) {
    if (remove.invalid) {
      message(
        "CleanSincastAtlas: Remove invalid pseudobulk atlas  as 'remove.invalid = TRUE'."
      )
      object@pseudobulk <- NULL
      test.SincastAtlas["pseudobulk"] <- "Empty"
    } else {
      message(
        "CleanSincastAtlas: ", "the pseudobulk atlas is invalid. Consider replace it by setting 'remove.invalid = TRUE'."
      )
    }
  }

  # Clean up the imputation atlas
  if (!is.imputation.valid) {
    if (remove.invalid) {
      message(
        "CleanSincastAtlas: Remove invalid imputation atlas as 'remove.invalid = TRUE'."
      )
      object@imputation <- NULL
      test.SincastAtlas["imputation"] <- "Empty"
    } else {
      message(
        "CleanSincastAtlas: ", "the imputation atlas is invalid. Consider replace it by setting 'remove.invalid = TRUE'."
      )
    }
  }

  message(
    "CleanSincastAtlas: After clean up: ",
    "\n original atlas: ", test.SincastAtlas["original"],
    "\n pseudobulk atlas: ", test.SincastAtlas["pseudobulk"],
    "\n imputation atlas: ", test.SincastAtlas["imputation"]
  )

  object
})

#' @rdname CleanSincastAtlas
setMethod("CleanSincastAtlas", "Sincast", function(object,
                                                    clean.up = c("none", "all", "original", "pseudobulk", "imputation"),
                                                    remove.invalid = TRUE, ...) {
  # Check the validity of the Sincast object.
  Sincast::CheckSincastObject(object, complete = FALSE, test = FALSE)

  SincastAtlas <- Sincast::GetSincastAtlas(object)
  SincastAtlas <- Sincast::CleanSincastAtlas(object = object,
                              clean.up = clean.up,
                              remove.invalid = remove.invalid, ...)
  Sincast::GetSincastAtlas(object) <- SincastAtlas
  object
})
