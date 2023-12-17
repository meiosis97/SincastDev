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
      "CheckSincastAtlas: original atlas: ", test.result["original"],
      "\n  pseudobulk atlas: ", test.result["pseudobulk"],
      "\n  imputation atlas: ", test.result["imputation"],
      collapse = ""
    )
  } else {
    problem <- NULL
  }

  # Print the error messages
  if (!is.null(problem)) {
    if (!test){
      problem <- gsub("\n", "\n  ", problem)
      stop(problem)
    }else if (!silent){
      message(problem)
    }
  }

  test.result
})

#' @rdname CheckSincastAtlas
setMethod("CheckSincastAtlas", "Sincast", function(object, test = TRUE, silent = FALSE, ...) {
  # Check the validity of the Sincast object.
  Sincast::CheckSincastObject(object, complete = FALSE, test = FALSE)

  object <- Sincast::GetSincastAtlas(object)
  Sincast::CheckSincastAtlas(
    object = object, test = test, silent = silent, ...
  )
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
    warning("GetSincastAtlas: At least one Sincast atlas is invalid.")
  }

  if (atlas != "all") {
    if (test.SincastAtlas[atlas] != "Valid") {
      warning("GetSincastAtlas: Invalid ", atlas, " atlas provided.")
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
#' @return A \code{Sincast} object with updated \code{Sincast} atlas.
#'
#' @family SincastAtlas related methods
#'
#' @export
#' @name CleanSincastAtlas
#' @rdname CleanSincastAtlas
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
    "\n  original atlas: ", test.SincastAtlas["original"],
    "\n  pseudobulk atlas: ", test.SincastAtlas["pseudobulk"],
    "\n  imputation atlas: ", test.SincastAtlas["imputation"]
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
    "\n  original atlas: ", test.SincastAtlas["original"],
    "\n  pseudobulk atlas: ", test.SincastAtlas["pseudobulk"],
    "\n  imputation atlas: ", test.SincastAtlas["imputation"]
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
  SincastAtlas <- Sincast::CleanSincastAtlas(
    object = object,
    clean.up = clean.up,
    remove.invalid = remove.invalid, ...
  )
  Sincast::GetSincastAtlas(object) <- SincastAtlas
  object
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# BuildSincastAtlas
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Build a \code{Sincast} atlas, which is a \code{Seurat} object with a rank-transformed reference
#' data, as well a pca dimensional reduction performed on the data.
#'
#' To be added.
#'
#' @param object A \code{Sincast} object.
#' @param cells Cells to be included in the atlas. Can be indexes, cell names or logical.
#' @param sincast.assay Either one of the "pseudobulk", "imputation" or "original" indicating the
#' type of atlas to build. Default is the active \code{Sincast} assay (Can be check by calling the \code{Sincast} object).
#' @param seurat.assay  Name of the \code{Seurat} assay rank transformation, data scaling and PCA are being run on.
#' @param layer Name of the \code{Seurat} layer rank transformation, data scaling and PCA are being run on.
#' @param reference.features Features on which the atlas is being built. Can be indexes, feature names, logical
#' or a single character corresponding to a logical attribute in \code{Seurat} assay's feature metadata. If the data
#' does not contain any of the given reference features, the missing features will be set to zero in the atlas construction.
#' @param query.features Must be a character vector of feature names. If provided, the atlas will be build on the
#' intersection between reference.features and query.features.
#' @param ScaleData.control A list of arguments to the \code{Seurat} \code{\link[Seurat]{ScaleData}} function call.
#' @param RunPCA.control A list of arguments to the \code{Seurat} \code{\link[Seurat]{RunPCA}} function call.
#' @param replace Logical; if TRUE, replace the existing atlas corresponding to the \code{sincast.assay}.
#'
#' @return A \code{Sincast} object with updated \code{Sincast} atlas.
#'
#' @family SincastAtlas related methods
#'
#' @export
#' @name BuildSincastAtlas
#' @rdname BuildSincastAtlas
setGeneric("BuildSincastAtlas", function(object,
                                         cells = NULL,
                                         sincast.assay = NULL,
                                         seurat.assay = NULL,
                                         layer = "counts",
                                         reference.features = NULL,
                                         query.features = NULL,
                                         ScaleData.control = list(),
                                         RunPCA.control = list(),
                                         replace = FALSE,
                                         ...) {
  standardGeneric("BuildSincastAtlas")
})

#' @rdname BuildSincastAtlas
setMethod("BuildSincastAtlas", "Sincast", function(object,
                                                   cells = NULL,
                                                   sincast.assay = NULL,
                                                   seurat.assay = NULL,
                                                   layer = "counts",
                                                   reference.features = NULL,
                                                   query.features = NULL,
                                                   ScaleData.control = list(),
                                                   RunPCA.control = list(),
                                                   replace = FALSE,
                                                   ...) {
  # Check the validity of the Sincast object.
  Sincast::CheckSincastObject(object, complete = FALSE, test = FALSE)

  # Get the Sincast assay.
  if (is.null(sincast.assay)) {
    active.assay <- summary(object)@active.assay
    test.SincastAssays <- Sincast::CheckSincastAssays(object, silent = TRUE)

    if (length(active.assay) & test.SincastAssays[active.assay] == "Valid") {
      sincast.assay <- active.assay
    } else {
      stop("BuildSincastAtlas: No valid active Sincast assay found.")
    }
  } else if (!IsCharacterVector(sincast.assay, 1, c("original", "imputation", "pseudobulk"))) {
    stop("BuildSincastAtlas: sincast.assay should be either one of 'original', 'imputation' or 'pseudobulk'.")
  }

  # If a desired atlas already exists, either replace it or return an error.
  if (!is.null(
    Sincast::GetSincastAtlas(object, atlas = sincast.assay)
  )) {
    if (!replace) {
      stop(
        "BuildSincastAtlas: An ", sincast.assay, " atlas already exists, set 'replace = T' to enforce a replacement."
      )
    } else {
      message("BuildSincastAtlas: An ", sincast.assay, " atlas already exists. Will be replaced as 'replace = T'.")
    }
  }

  seurat.object <- Sincast::GetSincastAssays(object = object, assay = sincast.assay)

  # Get the Seurat assay.
  if (is.null(seurat.assay)) {
    seurat.assay <- Seurat::DefaultAssay(seurat.object)
  } else {
    Seurat::DefaultAssay(seurat.object) <- seurat.assay
  }

  # Get cells.
  if (is.null(cells)){
    cells <- colnames(SeuratObject::LayerData(seurat.object,
                                              layer = layer))
  }

  # Get features.
  data.features <- rownames(SeuratObject::LayerData(seurat.object,
                                                   layer = layer))

  # Get the reference features
  if (is.null(reference.features)) {
    if (!is.null(Seurat::VariableFeatures(seurat.object))) {
      message("BuildSincastAtlas: reference.features not provided, set to the most variable features.")
      reference.features <- Seurat::VariableFeatures(seurat.object)
    } else {
      message("BuildSincastAtlas: reference.features not provided, set to all the features.")
      reference.features <- data.features
    }
  } else if (IsCharacterVector(reference.features, 1)) {
    feature.metadata <- Seurat::GetAssay(seurat.object)@metadata
    if (!reference.features %in% colnames(feature.metadata)) {
      stop("BuildSincastAtlas: reference.features is set to a single string that cannot be found in the feature metadata.")
    } else if (!is.logical(feature.metadata[, reference.features])) {
      stop("BuildSincastAtlas: reference.features is not a logical attribute in the feature metadata.")
    }
    reference.features <- rownames(feature.metadata)[feature.metadata[, reference.features]]
  } else if (IsCharacterVector(reference.features)) {
    missing.features <- reference.features[!reference.features %in% data.features]
    if (length(missing.features)) {
      warning("BuildSincastAtlas: ", paste(missing.features, collapse = ","), " are not found in the data, ",
      "they will be assigned with values of zero in the atlas.")
    }
  } else if (IsNumericVector(reference.features) | IsLogicalVector(reference.features, n = nrow(out))) {
    reference.features <- data.features[reference.features]
  } else {
    stop("BuildSincastAtlas: Invalid reference features.")
  }

  # Get the query features.
  if (IsCharacterVector(query.features)) {
    message(
      "BuildSincastAtlas: query.features provided, build the atlas based on the intersection between ",
      "the query.features and reference.features."
    )
    reference.features <- intersect(reference.features, query.features)
    if (length(reference.features)) {
      message("BuildSincastAtlas: ", length(reference.features), " common features found.")
    } else {
      stop("BuildSincastAtlas: No common features found.")
    }
  } else if(!is.null(query.features)){
    stop("BuildSincastAtlas: Invalid query features.")
  }

  # Subset data by cells.
  seurat.object <- seurat.object[,cells]

  # Subset data by the reference features.
  common.features <- intersect(data.features, reference.features)
  missing.features <- reference.features[!reference.features %in% data.features]
  seurat.object <- seurat.object[common.features,]
  out <- SeuratObject::LayerData(seurat.object, layer = layer)

  if (length(missing.features)) {
    tmp <- matrix(0, nrow = length(missing.features), ncol = ncol(out))
    rownames(tmp) <- missing.features
    out <- rbind(out[common.features, ], tmp)
  } else {
    out <- out[common.features, ]
  }

  message("BuildSincastAtlas: Ready to build an ", sincast.assay, "-", seurat.assay, "-",
          layer, " atlas containing ", nrow(out), " features, ", ncol(out), " cells.")
  message(rep("-", 100))

  # Rank transformation
  message("  Now perform rank transformation.")
  out <- Sincast::RankTrans(out)

  # Calculate sparisity
  sparsity <- round(mean(out == 0), 3)

  # Get metadata.
  meta.data <- seurat.object@meta.data

  # Generate a token for the this atlas construction run.
  SincastToken <- Seurat::Misc(seurat.object, slot = "SincastToken")

  if(is.null(SincastToken) | !is(SincastToken, "SincastToken")){
    # Generate a token for this imputation run.
    SincastToken <- GenerateSincastToken(
      by = "BuildSincastAtlas",
      command = deparse(match.call())
    )
  }else{
    SincastToken <- GenerateSincastToken(
      by = "SincastAggregate",
      command = deparse(match.call()),
      extend = SincastToken@command
    )
  }

  # Create a new Seurat object to store the result.
  suppressWarnings(
    out <- Seurat::CreateSeuratObject(
      counts = Matrix::Matrix(out),
      assay = seurat.assay,
      meta.data = meta.data,
      project = SincastToken@id, ...
    )
  )
  Seurat::Idents(out) <- Seurat::Idents(seurat.object)

  # Set layer data
  SeuratObject::LayerData(out, layer = "data") <-
    SeuratObject::LayerData(out, layer = "counts")

  # Scale data
  message("  Now scale the rank transformed data.")
  if(!is.list(ScaleData.control)){
    warning("BuildSincastAtlas: ScaleData.control is not a list. Using the default control.")
    ScaleData.control <- list()
  }
  if("object" %in% names(ScaleData.control)){
    ScaleData.control <- ScaleData.control[names(ScaleData.control)!="object"]
  }
  if("features" %in% names(RunPCA.control)){
    RunPCA.control <- RunPCA.control[names(RunPCA.control)!="features"]
  }
  if(!"do.scale" %in% names(ScaleData.control)){
    ScaleData.control[["do.scale"]] <- FALSE
  }
  if(!"verbose" %in% names(ScaleData.control)){
    ScaleData.control[["verbose"]] <- FALSE
  }
  if(length(ScaleData.control)){
    args <- paste(names(ScaleData.control), "=ScaleData.control$",
                  names(ScaleData.control), sep = "", collapse = ",")
    args <- paste("Seurat::ScaleData(object = out,", args, ")", sep = "")
    out <- eval(parse(text = args))
  }else{
    out <- Seurat::ScaleData(object = out)
  }

  # PCA
  message("  Now perform principal component analysis.")
  if(!is.list(RunPCA.control)){
    warning("BuildSincastAtlas: RunPCA.control is not a list. Using the default control.")
    RunPCA.control <- list()
  }
  if("object" %in% names(RunPCA.control)){
    RunPCA.control <- RunPCA.control[names(RunPCA.control)!="object"]
  }
  if("features" %in% names(RunPCA.control)){
    RunPCA.control <- RunPCA.control[names(RunPCA.control)!="features"]
  }
  if(!"verbose" %in% names(RunPCA.control)){
    RunPCA.control[["verbose"]] <- FALSE
  }
  if(length(RunPCA.control)){
    args <- paste(names(RunPCA.control), "=RunPCA.control$",
                  names(RunPCA.control), sep = "", collapse = ",")
    args <- paste("Seurat::RunPCA(object = out, features = reference.features, ",
                  args, ")", sep = "")
    out <- eval(parse(text = args))
  }else{
    out <- Seurat::RunPCA(object = out, features = reference.features)
  }
  npcs <- ncol(Seurat::Embeddings(out, reduction = "pca"))
  tot.var <- sum(SeuratObject::LayerData(out,layer = "scale.data")^2) / max(1, ncol(out) - 1)
  var.explained <- round(sum(Seurat::Stdev(out)^2) / tot.var,3)

  # Remove duplicated data
  SeuratObject::LayerData(out, layer = "data") <- NULL

  # Write the summary information
  idx <- paste(sincast.assay, ".atlas", sep = "")
  object@summary@summary[idx, ] <- c(
    seurat.assay, layer, nrow(out), ncol(out),
    npcs, var.explained, NA, sparsity
  )
  object@summary@active.atlas <- sincast.assay
  SincastToken@summary <- object@summary@summary[idx, ]

  # Add token to the new Seurat object,
  Seurat::Misc(out, slot = "SincastToken") <- SincastToken

  # Add or replace the desired atlas by the new atlas.
  Sincast::GetSincastAtlas(object, atlas = sincast.assay) <- out

  message(rep("-", 100))
  message(
    "BuildSincastAtlas: Built an atlas of ", npcs, " PCs, explained ",
    var.explained * 100, " percent of variance."
  )

  object

})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# show.SincastAtlas
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setMethod("show", "SincastAtlas", function(object) {
  # Convert Seurat show
  out <- "An object of class SincastAtlas"
  out <- paste(out, collapse = "\n")

  test.SincastAtlas <- Sincast::CheckSincastAtlas(object, silent = TRUE)

  if(test.SincastAtlas["pseudobulk"] == "Empty" &
     test.SincastAtlas["imputation"] == "Empty"){
    out <- paste(out, "\n", " No Sincast atlas present", sep = "")
  }

  # Print information for Sincast's pseudobulk atlas
  if (test.SincastAtlas["pseudobulk"] == "Valid") {
    summary <- Seurat::Misc(object@pseudobulk, "SincastToken")@summary
    out <- paste(out, "\n", " pseudobulk atlas: (",
                 "Computed using ", summary[,"assay"], "-", summary[,"layer"], ", ",
                 summary[,"nfeatures"], " features, ",
                 summary[,"nsamples"], " Samples, ",
                 summary[,"n.components"], " components, ",
                 summary[,"var.explained"], " percent variance explained, ",
                 summary[,"sparsity.after"], " sparsity)",
                 sep = ""
    )
  } else if(test.SincastAtlas["pseudobulk"] != "Empty"){
    out <- paste(out, "\n", " pseudobulk atlas: (", test.SincastAtlas["pseudobulk"], ")", sep = "")
  }

  # Print information for Sincast's imputation atlas
  if (test.SincastAtlas["imputation"] == "Valid") {
    summary <- Seurat::Misc(object@imputation, "SincastToken")@summary
    out <- paste(out, "\n", " imputation atlas: (",
                 "Computed using ", summary[,"assay"], "-", summary[,"layer"], ", ",
                 summary[,"nfeatures"], " features, ",
                 summary[,"nsamples"], " Samples, ",
                 summary[,"n.components"], " components, ",
                 summary[,"var.explained"], " percent variance explained, ",
                 summary[,"sparsity.after"], " sparsity)",
                 sep = ""
    )
  } else if(test.SincastAtlas["imputation"] != "Empty"){
    out <- paste(out, "\n", " imputation atlas: (", test.SincastAtlas["imputation"], ")", sep = "")
  }

  cat(out)
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# AtlasPlot
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Plot \code{Sincast} atlas.
#'
#' To be added.
#'
#' @param object A \code{Sincast} atlas, which suppose to be a \code{Seurat} object
#' generated by \code{\link[Sincast]{BuildSincastAtlas}}.
#' @param atlas Either one of the "pseudobulk", "imputation" or "original" indicating the
#' type of atlas to be plotted Default is the active \code{Sincast} atlas (Can be check by calling the \code{Sincast} object).
#' @param dims Dimensions to plot, must be a three-length numeric vector specifying x-, y- and z-dimensions.
#' @param cells Cells to plot.
#' @param color.by Color cells by which feature or meta.data attribute.
#' @param colors Color Scheme of col.by. Should be a named vector using color
#'  codes as values and labels of \code{color.by} as names.
#' @param anno.by Additional annotation of cells.
#'
#' @return A \code{plotly} object illustrating 3d principal components.
#'
#' @family Sincast plot methods
#'
#' @export
#' @name AtlasPlot
#' @rdname AtlasPlot
setGeneric("AtlasPlot", function(object,
                                      atlas = NULL,
                                      dims = 1:3,
                                      cells = NULL,
                                      color.by = "ident",
                                      colors = NULL,
                                      anno.by = NULL, ...) {
  standardGeneric("AtlasPlot")
})


#' @name AtlasPlot
#' @rdname AtlasPlot
setMethod("AtlasPlot", "Seurat", function(object,
                                               dims = 1:3,
                                               cells = NULL,
                                               color.by = "ident",
                                               colors = NULL,
                                               anno.by = NULL, ...) {
  if (!all(is.numeric(dims)) | length(dims) != 3) {
    warning("dims must be a three-length numeric vector.")
  }

  if (!is.null(cells)) object <- object[cells, ]

  dcs <- Seurat::Embeddings(object, reduction = "pca")
  dcs <- dcs[, dims]
  axis.labels <- colnames(dcs)
  colnames(dcs) <- c("x", "y", "z")


  # Generate color
  color <- NULL
  if(! IsCharacterVector(color.by, n= 1) ){
    warning("AtlasPlot: ", "color.by should be a single character.")
    color.by <- "ident"
  }

  if (color.by == "ident") {
    color <- Seurat::Idents(object)
  } else if (color.by %in% rownames(object)) {
    color <- suppressWarnings( SeuratObject::LayerData(object)[color.by, ] )
  } else if (color.by %in% colnames(object@meta.data)) {
    color <- object@meta.data[, color.by]
  } else {
    warning("AtlasPlot: ", color.by, " was not found in neither the data nor the metadata.")
  }

  # Generate annotation
  anno <- NULL

  if (!is.null(anno.by) & !IsCharacterVector(anno.by) ) {
    warning("AtlasPlot: ", "anno.by should be a character vector.")
    anno.by <- NULL
  }

  for (i in anno.by) {
    if (i == "ident") {
      tmp <- paste(i, Seurat::Idents(object), sep = ":")
      anno <- paste(anno, tmp, sep = "\n")
    } else if (i %in% rownames(object)) {
      tmp <- suppressWarnings( paste(i, SeuratObject::LayerData(object)[i, ], sep = ":") )
      anno <- paste(anno, tmp, sep = "\n")
    } else if (i %in% colnames(object@meta.data)) {
      tmp <- paste(i, object@meta.data[, i], sep = ":")
      anno <- paste(anno, tmp, sep = "\n")
    } else {
      warning("AtlasPlot: ", i, " was not found in either the data nor the metadata.")
    }
  }

  # Generate plot
  fig <- plot_ly() %>%
    add_trace(
      data = data.frame(dcs), type = "scatter3d",
      mode = "markers",
      x = ~x, y = ~y, z = ~z,
      color = color,
      colors = colors,
      text = anno, ...
    ) %>%
    layout(
      legend = list(orientation = "h"),
      scene = list(
        xaxis = list(title = axis.labels[1]),
        yaxis = list(title = axis.labels[2]),
        zaxis = list(title = axis.labels[3])
      )
    )

  suppressMessages(suppressWarnings(fig))
})


#' @name AtlasPlot
#' @rdname AtlasPlot
setMethod("AtlasPlot", "Sincast", function(object,
                                                atlas = NULL,
                                                dims = 1:3,
                                                cells = NULL,
                                                color.by = "ident",
                                                colors = NULL,
                                                anno.by = NULL, ...) {

  # Get the Sincast assay.
  if (is.null(atlas)) {
    active.atlas <- summary(object)@active.atlas
    test.SincastAtlas<- Sincast::CheckSincastAtlas(object, silent = TRUE)

    if (length(active.atlas) & test.SincastAtlas[active.atlas] == "Valid") {
      atlas <- active.atlas
    } else {
      stop("AtlasPlot: No valid active Sincast atlas found.")
    }
  } else if (!IsCharacterVector(atlas, 1, c("original", "imputation", "pseudobulk"))) {
    stop("AtlasPlot: atlas should be either one of 'original', 'imputation' or 'pseudobulk'.")
  }


  # Check the validity of the Sincast object.
  Sincast::CheckSincastObject(object, complete = FALSE, test = FALSE)

  object <- Sincast::GetSincastAtlas(object, atlas = atlas)
  Sincast::AtlasPlot(
    object = object, dims = dims, cells = cells,
    color.by = color.by, colors = colors,
    anno.by = anno.by, ...
  )
})

