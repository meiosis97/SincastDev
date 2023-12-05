# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# as.Sincast
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Convert a \code{Seurat} object to a \code{Sincast} object.
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
#' @name as.Sincast
#' @rdname as.Sincast
#' @aliases Sincast, SincastObject, Seurat
setGeneric("as.Sincast", function(object, ...) {
  standardGeneric("as.Sincast")
})

#' @rdname as.Sincast
setMethod("as.Sincast", "Seurat", function(object, ...) {
  # Generate a Sincast token.
  Seurat::Misc(object, slot = "SincastToken") <-
    GenerateSincastToken(by = "as.Sincast", command = deparse(match.call()))

  # Generate a SincastAssays object.
  object <- new("SincastAssays", original = object)

  # Generate a Sincast object.
  object <- new("Sincast", SincastAssays = object,
                summary = new("SincastSummary"))

  object
})

#' @rdname as.Sincast
setMethod("as.Sincast", "Sincast", function(object, ...) {
  Sincast::CheckSincastObject(object, complete = FALSE)
  object
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CheckSincastObject
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Check the validity of the \code{Sincast} object.
#'
#' To be added.
#'
#' @param object A \code{Sincast} object.
#' @param test Logical; if TRUE (the default) and validity fails, the function returns a
#'  vector of strings describing the problems. If test is FALSE validity failure generates an error.
#' @param complete Logical; if TRUE, call validity check for each slot of the
#' existing \code{Sincast} object.
#' @param silent Logical; if TRUE, suppress all messages.
#'
#' @return A vector of strings describing the problems.
#'
#' @family SincastObject related methods
#'
#' @export
#' @name CheckSincastObject
#' @rdname CheckSincastObject
#' @aliases Sincast, SincastObject
setGeneric("CheckSincastObject", function(object, test = TRUE,
                                          complete = TRUE, slient = FALSE, ...) {
  standardGeneric("CheckSincastObject")
})


#' @rdname CheckSincastObject
setMethod("CheckSincastObject", "Sincast", function(object, test = TRUE,
                                                   complete = TRUE, silent = FALSE, ...) {
  problem <- NULL
  test.SincastObject <- "Valid"
  test.SincastAssays <- NULL
  test.SincastAtlas <- NULL

  if (FALSE) {
    # %%%%%%%%%% TO BE ADDED: TESTS FOR THE MAIN OBJECT %%%%%%%%% #

  } else if (complete) {

    if (!silent) {
      message(
        "CheckSincastObject: Check the validity of each slot of the Sincast object, as 'complete = TRUE'."
      )
    }

    #  Check the validity of SincastAssays.
    test.SincastAssays <- Sincast::CheckSincastAssays(object@SincastAssays, silent = TRUE)

    if (!all(test.SincastAssays %in% c("Valid", "Empty"))) {
      problem <- paste(
        "CheckSincastObject: Check SincastAssays: ",
        "\n original Seurat: ", test.SincastAssays["original"],
        "\n pseudobulk assay: ", test.SincastAssays["pseudobulk"],
        "\n imputation assay: ", test.SincastAssays["imputation"],
        collapse = ""
      )
    }

    #  Check the validity of SincastAtlas
    test.SincastAtlas <- Sincast::CheckSincastAtlas(object@SincastAtlas, silent = TRUE)

    if (!all(test.SincastAtlas %in% c("Valid", "Empty"))) {
      problem <- paste(
        "CheckSincastObject: Check SincastAtlas: ",
        "\n original atlas: ", test.SincastAtlas["original"],
        "\n pseudobulk atlas: ", test.SincastAtlas["pseudobulk"],
        "\n imputation atlas: ", test.SincastAtlas["imputation"],
        collapse = ""
      )
    }
  }

  # Print error messages.
  if (!is.null(problem)) {
    if (!test) stop(problem) else if (!silent) message(problem)
  }

  attr(test.SincastObject, "Sincastassays") <- test.SincastAssays
  attr(test.SincastObject, "SincastAtlas") <- test.SincastAtlas
  test.SincastObject

})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# show.Sincast
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setMethod("show", "Sincast", function(object) {
  Sincast::CheckSincastObject(object, complete = FALSE)
  out <- capture.output(object@SincastAssays)
  out[1] <- "An object of class Sincast"
  out <- paste(out, collapse = "\n")
  cat(out)
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# generics
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setMethod(
  f = '[[',
  signature = c(
    x = 'Sincast',
    i = 'character',
    j = 'missing'
  ),
  definition = function(x, i, ...) {
    Sincast::GetSincastAssays(x, assay = i, ...)
  }
)

setMethod(
  f = '[[',
  signature = c(
    x = 'Sincast',
    i = 'character',
    j = 'character'
  ),
  definition = function(x, i, j = c("assay", "atlas"),...) {
    j <- match.arg(j)
    if(j == "assay") Sincast::GetSincastAssays(x, assay = i, ...)
    if(j == "atlas") Sincast::GetSincastAtlas(x, atlas = i, ...)
  }
)

setMethod(
  f = '[[<-',
  signature = c(
    x = 'Sincast',
    i = 'character',
    j = 'missing',
    value = "ANY"
  ),
  definition = function(x, i, j, ..., value) {
    Sincast::GetSincastAssays(x, assay = i, ...) <- value
    x
  }
)


setMethod(
  f = '[[<-',
  signature = c(
    x = 'Sincast',
    i = 'character',
    j = 'character',
    value = "ANY"
  ),
  definition = function(x, i, j = c("assay", "atlas"),..., value) {
    j <- match.arg(j)
    if(j == "assay") Sincast::GetSincastAssays(x, assay = i, ...) <- value
    if(j == "atlas") Sincast::GetSincastAtlas(x, atlas = i, ...) <- value
    x
  }
)
