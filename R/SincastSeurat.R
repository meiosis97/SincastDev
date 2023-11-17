# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# "as" method for S4 class "Seurat" and "SincastSeurat"
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# As("Seurat", "SincastSeurat")
setAs("Seurat", "SincastSeurat", function(from, to) {
  SincastObject <- CreateSincastObject(by = "as", command = deparse(match.call()))
  arguments <- paste(slotNames(from), "=from@", slotNames(from), sep = "", collapse = ",")
  arguments <- paste(arguments, "Sincast=SincastObject", sep = ",")
  text2expr <- paste("SincastSeurat(", arguments, ")", sep = "")
  eval(parse(text = text2expr))
})


# As("SincastSeurat", "Seurat")
setAs("SincastSeurat", "Seurat", function(from, to) {
  arguments <- paste(slotNames(to), "=from@", slotNames(to), sep = "", collapse = ",")
  text2expr <- paste("new('Seurat',", arguments, ")", sep = "")
  eval(parse(text = text2expr))
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# as.SincastSeurat
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Convert a \code{Seurat} object to a \code{SincastSeurat} object.
#'
#' To be added.
#'
#' @param object A \code{Seurat} or a \code{SincastSeurat} object.
#'
#' @return A \code{SincastSeurat} object.
#'
#' @family Sincast related methods
#'
#' @export
#' @name as.SincastSeurat
#' @rdname as.SincastSeurat
#' @aliases Sincast, Seurat
setGeneric("as.SincastSeurat", function(object, ...) {
  standardGeneric("as.SincastSeurat")
})

setMethod("as.SincastSeurat", "Seurat", function(object, ...) {
  if(class(object) == 'Seurat'){
    object <- as(object, "SincastSeurat")
    object@Sincast <- CreateSincastObject(by = "as.SincastSeurat",
                                          command = deparse(match.call()))

  }else if(class(object) != "SincastSeurat"){
    stop("Method not defined for class ", class(object))
  }

  object

})
