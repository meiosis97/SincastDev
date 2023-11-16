# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CreateSincastObject
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Create an empty \code{Sincast} object.
#'
#' To be added.
#'
#' @return An empty \code{Sincast} object.
#'
#' @family Sincast related methods
#'
#' @export
#' @rdname CreateSincastObject
#' @aliases Sincast, SincastAssays
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
# AddSincastObject
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Add a \code{Sincast} object to \code{Seurat}'s \code{misc} slot.
#'
#' To be added.
#'
#' @param object A \code{Seurat} object.
#' @param replace Logical; if \code{TRUE}, replace the existing \code{Sincast} object in the \code{misc} slot.
#'
#' @return An updated \code{Seurat} object with an additional \code{Sincast} object in the \code{misc} slot.
#'
#' @family Sincast related methods
#'
#' @export
#' @rdname AddSincastObject
#' @aliases Sincast, SincastAssays, Seurat
setGeneric("AddSincastObject", function(object,
                                        replace = FALSE, ...) {
  standardGeneric("AddSincastObject")
})

#' @rdname AddSincastObject
setMethod("AddSincastObject", "Seurat", function(object,
                                                 replace = FALSE, ...) {
  SincastObject <- Seurat::Misc(object, slot = "Sincast")

  if (is.null(SincastObject)) {
    # Generate a Sincast object,
    SincastObject <- Sincast::CreateSincastObject(
      by = "AddSincastObject",
      command = deparse(match.call())
    )
  } else {
    message(
      "AddSincastObject: 'Sincast' object already exists. Check the validity of the existing object."
    )
    CheckSincastObject(object)

    if (replace) {
      message(
        "AddSincastObject: Replacing the existing 'Sincast' object as replace = TRUE."
      )
      # Generate a Sincast object,
      SincastObject <- Sincast::CreateSincastObject(
        by = "AddSincastObject",
        command = deparse(match.call())
      )
    } else {
      message(
        "AddSincastObject: A 'Sincast‘ object already exists, set 'replace = T' to enforce a replacement,",
        "or use 'CleanSincastAssays' function to modify 'SincastAssays' specifically."
      )
    }
  }

  suppressWarnings(
    Seurat::Misc(object, slot = "Sincast") <- SincastObject
  )

  object
})


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
#' @aliases Sincast, SincastAssays, Seurat
setGeneric("CheckSincastObject", function(object, test = TRUE,
                                          complete = TRUE, slient = FALSE, ...) {
  standardGeneric("CheckSincastObject")
})

#' @rdname CheckSincastObject
setMethod("CheckSincastObject", "Seurat", function(object, test = TRUE,
                                                   complete = TRUE, silent = FALSE, ...) {

  SincastObject <- Seurat::Misc(object, slot = "Sincast")
  missing.object <- FALSE
  wrong.class <- FALSE
  problem <- NULL
  test.SincastAssays <- NULL

  if (is.null(SincastObject)) {
    problem <- "'Sincast object' doesn't exist."
    missing.object <- TRUE
  } else if (is(SincastObject, "SincastObject")) {
    problem <- "Unrecroglized 'Sincast' object."
    wrong.class <- TRUE
  } else if (complete) {
    test.SincastAssays <- validObject(SincastObject@SincastAssays, test = TRUE)

    if (!silent) {
      message(
        "CheckSincastObject: A 'Sincast' object exists. Check the validity of each slot."
      )
    }

    if(!all(test.SincastAssays %in% c("Valid", "Empty"))){
      problem <- paste(
        "CheckSincastObject: Check 'SincastAssays': ",
        "\n \t 'pseudobulk' assay: ", test.SincastAssays["pseudobulk"],
        "; 'imputation' assay: ", test.SincastAssays["imputation"], collapse = ''
      )
    }

  }

  # Print error messages.
  if (!is.null(problem)) {
    if (!test) strop(problem) else if (!silent) message(problem)
  }

  out <- c(missing.object = missing.object, wrong.class = wrong.class)
  attr(out, 'Sincastassays') <- test.SincastAssays
  out
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
#' @aliases Sincast, SincastAssays, Seurat
setGeneric("GetSincastObject", function(object, ...) {
  standardGeneric("GetSincastObject")
})

#' @rdname GetSincastObject
setMethod("GetSincastObject", "Seurat", function(object, ...) {
  # Check the validity of the "Sincast" object in "Seurat"'s "misc" slot.
  test.SincastObject <- Sincast::CheckSincastObject(object, complete = FALSE)

  # If the "Sincast" object is missing, or invalid, return a NULL
  if (any(test.SincastObject)) {
    out <- NULL
  } else {
    out <- Seurat::Misc(object, slot = "Sincast")
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
#' @aliases Sincast, SincastAssays, Seurat
setGeneric("GetSincastObject<-", function(object, value, ...) {
  standardGeneric("GetSincastObject<-")
})

#' @rdname GetSincastObject
setMethod("GetSincastObject<-", "Seurat", function(object, value, ...) {
  suppressWarnings(Seurat::Misc(object, slot = "Sincast") <- value)

  object
})
