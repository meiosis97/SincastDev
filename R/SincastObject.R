# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CreateSincastObject
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Create an empty \code{Sincast} object.
#'
#' To be added.
#'
#' @slot by A string recording the function by which \code{CreateSincastObject} was called.
#' @slot command A string recording \code{Sincast} command history.
#'
#' @return An empty \code{Sincast} object.
#'
#' @family SincastObject related methods
#'
#' @export
#' @rdname CreateSincastObject
#' @aliases Sincast, SincastObject
CreateSincastObject <- function(by = "CreateSincastObject", command = deparse(match.call())) {
  # Generate Sincast assays.
  SincastAssays <- new("SincastAssays")
  # Generate a Sincast token.
  SincastToken <- GenerateSincastToken(by = by)
  ids <- c(names(SincastToken@command), SincastToken@id)
  SincastToken@command <- c(SincastToken@command, command)
  names(SincastToken@command) <- ids

  # Generate a Sincast object.
  SincastObject <- new("Sincast",
    SincastAssays = SincastAssays,
    SincastToken = SincastToken
  )

  SincastObject
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CheckSincastObject
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Check the validity of the \code{Sincast} object in \code{Seurat}'s \code{misc} slot.
#'
#' To be added.
#'
#' @param object A \code{Seurat} object.
#' @param test Logical; if TRUE (the default) and validity fails, the function returns a
#'  vector of strings describing the problems. If test is FALSE validity failure generates an error.
#' @param complete Logical; if TRUE, call validity check for each slot of the
#' existing \code{Sincast} object.
#' @param silent Logical; if TRUE, suppress all messages.
#'
#' @return Logical indicates whether \code{Seurat}'s \code{misc} slot
#' contains a valid \code{Sincast} object.
#'
#' @family SincastObject related methods
#'
#' @export
#' @rdname CheckSincastObject
#' @aliases Sincast, SincastObject, Seurat
setGeneric("CheckSincastObject", function(object, test = TRUE,
                                          complete = TRUE, slient = FALSE, ...) {
  standardGeneric("CheckSincastObject")
})

#' @rdname CheckSincastObject
setMethod("CheckSincastObject", "Seurat", function(object, test = TRUE,
                                                   complete = TRUE, silent = FALSE, ...) {
  if ("Sincast" %in% slotNames(object)) {
    SincastObject <- object@Sincast
  } else {
    SincastObject <- NULL
  }

  problem <- NULL
  test.SincastObject <- "Valid"
  test.SincastAssays <- NULL

  if (is.null(SincastObject)) {
    problem <- "CheckSincastObject: 'Sincast' object doesn't exist."
    test.SincastObject <- "Missing object"
  } else if (is(SincastObject, "SincastObject")) {
    problem <- "CheckSincastObject: Unrecroglized 'Sincast' object."
    test.SincastObject <- "Wrong class"
  } else if (complete) {
    test.SincastAssays <- validObject(SincastObject@SincastAssays, test = TRUE)

    if (!silent) {
      message(
        "CheckSincastObject: A 'Sincast' object exists. Check the validity of each slot."
      )
    }

    if (!all(test.SincastAssays %in% c("Valid", "Empty"))) {
      problem <- paste(
        "CheckSincastObject: Check 'SincastAssays': ",
        "\n \t 'pseudobulk' assay: ", test.SincastAssays["pseudobulk"],
        "; 'imputation' assay: ", test.SincastAssays["imputation"],
        collapse = ""
      )
    }
  }

  # Print error messages.
  if (!is.null(problem)) {
    if (!test) strop(problem) else if (!silent) message(problem)
  }

  attr(test.SincastObject, "Sincastassays") <- test.SincastAssays
  test.SincastObject
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GetSincastObject
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Extract the \code{Sincast} object from \code{Seurat}'s \code{misc} slot.
#'
#' To be added.
#'
#' @param object A \code{Seurat} object.
#'
#' @return A \code{Sincast} object.
#'
#' @family SincastObject related methods
#'
#' @export
#' @rdname GetSincastObject
#' @aliases Sincast, SincastObject, Seurat
setGeneric("GetSincastObject", function(object, ...) {
  standardGeneric("GetSincastObject")
})

#' @rdname GetSincastObject
setMethod("GetSincastObject", "Seurat", function(object, ...) {
  # Check the validity of the "Sincast" object.
  test.SincastObject <- Sincast::CheckSincastObject(object, complete = FALSE)

  # If the "Sincast" object is missing, or invalid, return a NULL
  if (test.SincastObject != "Valid") {
    out <- NULL
  } else {
    out <- object@Sincast
  }

  out
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GetSincastObject<-
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Setter function for \code{GetSincastObject}.
#'
#' To be added.
#'
#' @param object A \code{Seurat} object.
#'
#' @return A \code{Seurat} object with updated \code{Sincast} object
#'
#' @family SincastObject related methods
#'
#' @export
#' @rdname GetSincastObject
#' @aliases Sincast, SincastObject, Seurat
setGeneric("GetSincastObject<-", function(object, value, ...) {
  standardGeneric("GetSincastObject<-")
})

#' @rdname GetSincastObject
setMethod("GetSincastObject<-", "Seurat", function(object, value, ...) {
  if (!"Sincast" %in% slotNames(object)) {
    message("GetSincastObject: Convert to a SincastSeurat object.")
    object <- as.SincastSeurat(object)
  }
  object@Sincast <- value

  object
})
