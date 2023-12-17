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

  object@summary@active.assay <- "original"

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


  # %%%%%%%%%% TO BE ADDED: TESTS FOR THE MAIN OBJECT %%%%%%%%% #
  test.SincastSummary <- methods::validObject(object@summary, test = TRUE)
  if(!is.logical(test.SincastSummary)){
    # Update error message
    tmp <- paste(test.SincastSummary, collapse = "\n  ")
    problem <- paste("CheckSincastObject (SincastSummary): ",
                     tmp, sep = "")

    # Update return result.
    tmp <- paste(test.SincastSummary, collapse = "; ")
    test.SincastObject <- paste("SincastSummary:", tmp)

  }

  if (complete) {
    tmp <- "CheckSincastObject: Check the validity of each slot of the Sincast object, as 'complete = TRUE'."
    problem <- if(is.null(problem)) tmp else paste(problem, tmp, sep = "\n")

    #  Check the validity of SincastAssays.
    test.SincastAssays <- Sincast::CheckSincastAssays(object@SincastAssays, silent = TRUE)

    if (!all(test.SincastAssays %in% c("Valid", "Empty"))) {
      tmp <- paste(
        "CheckSincastObject (SincastAssays): original Seurat: ", test.SincastAssays["original"],
        "\n  pseudobulk assay: ", test.SincastAssays["pseudobulk"],
        "\n  imputation assay: ", test.SincastAssays["imputation"],
        sep = ""
      )
      problem <- paste(problem, tmp, sep = "\n")
    }

    #  Check the validity of SincastAtlas
    test.SincastAtlas <- Sincast::CheckSincastAtlas(object@SincastAtlas, silent = TRUE)

    if (!all(test.SincastAtlas %in% c("Valid", "Empty"))) {
      tmp <- paste(
        "CheckSincastObject (SincastAtlas): original atlas: ", test.SincastAtlas["original"],
        "\n  pseudobulk atlas: ", test.SincastAtlas["pseudobulk"],
        "\n  imputation atlas: ", test.SincastAtlas["imputation"],
        sep = ""
      )
      problem <- paste(problem, tmp, sep = "\n")
    }
  }

  # Print error messages.
  if (!is.null(problem)) {
    if (!test){
      problem <- gsub("\n", "\n  ", problem)
      stop(problem)
    }else if (!silent){
      message(problem)
    }
  }

  attr(test.SincastObject, "Sincastassays") <- test.SincastAssays
  attr(test.SincastObject, "SincastAtlas") <- test.SincastAtlas
  test.SincastObject

})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# show.Sincast
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setMethod("show", "Sincast", function(object) {
  # Shoe the Sincast assay.
  Sincast::CheckSincastObject(object, complete = FALSE)
  out <- capture.output(object@SincastAssays)
  out[1] <- "An object of class Sincast"

  # Show the active assay.
  test.SincastAssays <- Sincast::CheckSincastAssays(object, silent = TRUE)
  active.assay <- summary(object)@active.assay
  idx <- which(out == "Sincast assay:")
  if(length(active.assay)){
    if(test.SincastAssays[active.assay] == "Valid"){
      out[idx] <- paste(out[idx], "\n", " Active assay: ", active.assay, sep = "")
    }else{
      out[idx] <- paste(out[idx], "\n", " Active assay: none", sep = "")
    }
  }

  out <- paste(out, collapse = "\n")

  # Shoe the Sincast atlas.
  out2 <- capture.output(object@SincastAtlas)
  # Show the active atlas
  test.SincastAtlas <- Sincast::CheckSincastAtlas(object, silent = TRUE)
  active.atlas <- summary(object)@active.atlas
  if(length(active.atlas)){
    if(test.SincastAtlas[active.atlas] == "Valid"){
      out2[1] <- paste("Sincast atlas:", "\n", " Active atlas: ", active.atlas, sep = "")
    }else{
      out2[1] <- paste("Sincast atlas:", "\n", " Active atlas: none", sep = "")
    }
  }
  out2 <- paste(out2, collapse = "\n")

  out <- paste(out, out2, sep = "\n")

  cat(out)
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# generics
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setMethod(
  f = '[',
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
  f = '[',
  signature = c(
    x = 'Sincast',
    i = 'missing',
    j = 'missing'
  ),
  definition = function(x, i, ...) {
    Sincast::GetSincastAssays(x, ...)
  }
)


setMethod(
  f = '[<-',
  signature = c(
    x = 'Sincast',
    i = 'missing',
    j = 'missing',
    value = "ANY"
  ),
  definition = function(x, i, ..., value) {
    Sincast::GetSincastAssays(x, assay = i, ...) <- value
    x
  }
)

setMethod(
  f = '[<-',
  signature = c(
    x = 'Sincast',
    i = 'missing',
    j = 'missing',
    value = "SincastAssays"
  ),
  definition = function(x, i, ..., value) {
    Sincast::GetSincastAssays(x, ...) <- value
    x
  }
)

setMethod(
  f = '[[',
  signature = c(
    x = 'Sincast',
    i = 'character',
    j = 'missing'
  ),
  definition = function(x, i, ...) {
    Sincast::GetSincastAtlas(x, atlas = i, ...)
  }
)

setMethod(
  f = '[[',
  signature = c(
    x = 'Sincast',
    i = 'missing',
    j = 'missing'
  ),
  definition = function(x, i, ...) {
    Sincast::GetSincastAtlas(x, ...)
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
  definition = function(x, i, ..., value) {
    Sincast::GetSincastAtlas(x, atlas = i, ...) <- value
    x
  }
)

setMethod(
  f = '[<-',
  signature = c(
    x = 'Sincast',
    i = 'missing',
    j = 'missing',
    value = "SincastAtlas"
  ),
  definition = function(x, i, ..., value) {
    Sincast::GetSincastAtlas(x, ...) <- value
    x
  }
)



# setMethod(
#   f = '[[<-',
#   signature = c(
#     x = 'Sincast',
#     i = 'character',
#     j = 'character',
#     value = "ANY"
#   ),
#   definition = function(x, i, j = c("assay", "atlas"),..., value) {
#     j <- match.arg(j)
#     if(j == "assay") Sincast::GetSincastAssays(x, assay = i, ...) <- value
#     if(j == "atlas") Sincast::GetSincastAtlas(x, atlas = i, ...) <- value
#     x
#   }
# )
#
#
#
# setMethod(
#   f = '[[',
#   signature = c(
#     x = 'Sincast',
#     i = 'character',
#     j = 'character'
#   ),
#   definition = function(x, i, j = c("assay", "atlas"),...) {
#     j <- match.arg(j)
#     if(j == "assay") Sincast::GetSincastAssays(x, assay = i, ...)
#     if(j == "atlas") Sincast::GetSincastAtlas(x, atlas = i, ...)
#   }
# )
