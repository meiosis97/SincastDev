# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# "as" method for S4 class "Seurat" and "SincastSeurat"
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# As("Seurat", "SincastSeurat")
setAs("Seurat", "SincastSeurat", function(from, to) {
  arguments <- paste(slotNames(from), "=from@", slotNames(from), sep = "", collapse = ",")
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
#' @family SincastSeurat related methods
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


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# show.SincastSeurat
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setMethod("show", "SincastSeurat", function(object) {
  # Convert Seurat show
  out <- capture.output(as(object, "Seurat"))
  out[1] <- "An object of class SincastSeurat"
  out <- paste(out, collapse = "\n")


  SincastAssays <- object@Sincast@SincastAssays
  test.SincastAssays <- validObject(SincastAssays, test = "TRUE")

  sparsity <-  1-mean(object$nFeature_RNA/nrow(object))
  out <- paste(out,"\nSincast assay: Original data sparisty:", round(sparsity,3))

  # Print information for "Sincast"'s "pseudobulk" assay.
  if(test.SincastAssays["pseudobulk"] == "Valid"){
    n.features <- nrow(SincastAssays@pseudobulk)
    n.samples <- ncol(SincastAssays@pseudobulk)
    sparsity <-  1-mean(SincastAssays@pseudobulk$nFeature_RNA/n.features)
    out <- paste(out, "\n", " pseudobulk assay: (", n.features, " features, ",
                 n.samples, " Samples, ", round(sparsity,3) , " sparsity)", sep = "")
  }else{
    out <- paste(out, "\n", " pseudobulk assay: (", test.SincastAssays["pseudobulk"], ")", sep = "")
  }

  # Print information for "Sincast"'s "imputation" assay.
  if(test.SincastAssays["imputation"] == "Valid"){
    n.features <- nrow(SincastAssays@imputation)
    n.samples <- ncol(SincastAssays@imputation)
    sparsity <-  1-mean(SincastAssays@imputation$nFeature_RNA/n.features)
    out <- paste(out, "\n", " imputation assay: (", n.features, " features, ",
                 n.samples, " Samples, ", round(sparsity,3) , " sparsity)", sep = "")
  }else{
    out <- paste(out, "\n", " imputation assay: (", test.SincastAssays["imputation"], ")", sep = "")
  }

  cat(out)
})
