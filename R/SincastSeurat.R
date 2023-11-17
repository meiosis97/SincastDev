setAs("Seurat", "SincastSeurat", function(from, to) {
  arguments <- paste(slotNames(from), "=from@", slotNames(from), sep = "", collapse = ",")
  text2expr <- paste("SincastSeurat(", arguments, ")", sep = "")
  eval(parse(text = text2expr))
})

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
#' @rdname as.SincastSeurat
#' @aliases Sincast, SincastAssays, Seurat
setGeneric("as.SincastSeurat", function(object, ...) {
  standardGeneric("as.SincastSeurat")
})

setMethod("as.SincastSeurat", "Seurat", function(object, ...) {
  if(class(object) == 'Seurat'){
    object <- as(object, "SincastSeurat")
    object@Sincast <- CreateSincastObject()

  }else if(class(object) == "SincastSeurat"){
    if(is.null(object@Sincast)) object@Sincast <- CreateSincastObject()

  }else{
    stop("Method not defined for class ", class(object))
  }

  object

})
