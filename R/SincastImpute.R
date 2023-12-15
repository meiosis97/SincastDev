# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Not exported satellite functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pdist <- function(tmat) {
  # @param tmat A non-negative matrix with samples by features
  # @reference http://r.789695.n4.nabble.com/dist-function-in-R-is-very-slow-td4738317.html
  mtm <- Matrix::tcrossprod(tmat)
  sq <- rowSums(tmat^2)
  out0 <- outer(sq, sq, "+") - 2 * mtm
  out0[out0 < 0] <- 0

  sqrt(out0)
}


ScaleDistance <- function(dist) {
  for (i in 1:nrow(dist)) {
    r <- sort(dist[i, ], partial = 2)[2]
    dist[i, ] <- dist[i, ] - r
    dist[i, ][dist[i, ] <= 0] <- 0
  }
  dist
}


Symmetrization <- function(aff, norm) {
  if (norm == "probabilistic") {
    aff <- aff + t(aff) - aff * t(aff)
  }

  if (norm == "average") {
    aff <- (aff + t(aff)) / 2
  }
  aff
}


Laplacian <- function(aff) {
  Z <- rowSums(aff)
  Z.mat <- tcrossprod(Z)
  aff <- aff / Z.mat
}


FindSigma <- function(dk, a, k) {
  lower <- 0
  upper <- Inf
  cur <- dk[k]
  while (T) {
    psum <- sum(exp(-0.5 * (dk / cur)^2))
    if (psum > a) {
      upper <- cur
      cur <- (lower + cur) / 2
    } else if (psum < a) {
      if (is.infinite(upper)) {
        lower <- cur
        cur <- 2 * cur
      } else {
        lower <- cur
        cur <- (upper + cur) / 2
      }
    }
    if (abs(psum - a) < 1e-5) break
  }
  cur
}


QQSCale <- function(Y, X){
  n <- ncol(Y)
  if(n > 1000){
    q <- seq(0,1,length.out = n/10)
  }

  for(i in 1:nrow(Y)){
    if(n > 1000){
      y <- quantile(Y[i,], q)
      x <- quantile(X[i,], q)
    }else{
      y <- sort(Y[i,])
      x <- sort(X[i,])
    }
    mod <-  lm(y~0+x)
    X[i,] <- as.numeric( X[i,] * coef(mod)[1] )
  }
  X
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Low rank approximation.
#'
#' Perform the scRNA-seq low rank approximation method proposed by
#' \href{https://www.nature.com/articles/s41467-021-27729-z}{George C. Linderman (2022, nat communication)}.
#' The major contribution of the method is that it can perform imputation while preserve biological zeros.
#'
#' @param x a non-negative feature-by-cell matrix.
#' @param k Number of singular values requested for single value decomposition. If \code{estimate.rank} is set to TRUE, default is 100, otherwise 50.
#' @param estimate.rank Logical; If TRUE, the rank of \code{x} will be estimated by the approach proposed by Linderman et al (2022),
#' which will be less or equal to \code{k}. Otherwise the rank will be set to \code{k}. Also note that if TRUE, \code{k} should be a numeric larger than 100.
#' @param do.rsvd Logical; If TRUE, perform random single value decomposition.
#' @param seed The random seed for random SVD if do.rsvd = TRUE.
#' @param q If a svd-approximation of a feature contains negative values, calculate the qth smallest negative value. Any imputed expression that are lower
#' than the absolute of this value is set to zero.
#'
#' @return A low rank approximation of \code{x}
#'
#' @export
#'
#' @name LowRankApprox
#' @rdname LowRankApprox
LowRankApprox <-  function(x,
                           k = NULL,
                           estimate.rank = TRUE,
                           do.rsvd = TRUE,
                           seed = 521626,
                           q = 0.9, ...){
  x <- t(x)

  if(estimate.rank){
    # Number of components to calculate.
    if(is.null(k)){
      k <- 100
    }else if(!is.numeric(k)){
      warning("k is not a numeric larger than 100, reset k to 100.")
      k <- 100
    }else if(k < 100){
      warning("ran' is not a numeric larger than 100, reset k to 100.")
      k <- 100
    }

  }else{
    # Number of components to calculate.
    if(is.null(k)){
      k <- 50
    }else if(!is.numeric(k)){
      warning("k is not a numeric, reset k to 50.")
      k <- 50
    }

  }

  # Single Value decomposition.
  if(do.rsvd){
    set.seed(seed)
    svd.x <- rsvd::rsvd(x, k, ...)
  }else{
    svd.x <- RSpectra::svds(x, k, ...)
  }

  # Determine the number of rank to use if rank is not given.
  if(estimate.rank){
    s <- abs(diff(svd.x$d))
    mu <- mean(s[(0.8*k):k-1])
    sigma <- sd(s[(0.8*k):k-1])
    sk <- mu + 6*sigma
    rank <- which(s < sk)[1]
    if(length(rank)==0){
      warning(rank, " single values could be insufficient to approximate the data.")
      rank <- k
    }
  }else{
    rank <- k
  }

  # Low rank approximation
  x.lra <- tcrossprod(x %*% svd.x$v[,1:rank], svd.x$v[,1:rank])

  # Restore biological zeros.
  x.lra <- apply(x.lra, 2, function(x){
    x[x<abs(quantile(x[x<0], 1-q))] <- 0
    x
  } )

  dimnames(x.lra) <- dimnames(x)
  t(x.lra)

}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SincastImpute
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Sincast imputation
#'
#' Perform Sincast imputation.
#'
#' @param object A \code{Sincast} object.
#' @param assay Which \code{Seurat} assay to use. Default is the default \code{Seurat} assay.
#' @param features Features to impute. Default is all features in the assay.
#' @param npcs How many principal components to compute on the query data. The PCs will be used to construct the knn graph.
#' @param t Diffusion time, or the power of Markov transition matrix.
#' @param k The number of neighbors used to infer adaptive Gaussian kernel for each cell.
#' @param knn A cell can only be connected to knn neighbors when t is set 1.
#' @param do.umap.dist Logical; if TRUE, scale Euclidean distances to distances beyond nearest neighbors as in the UMAP algorithm.
#' @param a Default: log2(k) and log(k/(k-1)) when setting umap.dist to TRUE and FALSE respectively.
#'        a < 1 represents the probability of a cell communicating with its kth nearest neighbor.
#'        a > 1 represents the sum of probabilities of a cell communicating with its k nearest neighbors.
#' @param do.laplacian Logical; if TRUE, perform Laplacian normalization on the affinity matrix.
#' @param norm How to symmetrize the affinity matrix. Default is Probabilistic t-norm. The other option is 'average'.
#' @param do.post.scale Logical; if TRUE, perform post-imputation scaling, which match the scale of \code{Sincast} imputed data
#' with the data imputed by low rank approximation (See \code{\link[Sincast]{LowRankApprox}}). After scale matching,
#' the average of the low rank approximation and the scaled \code{Sincast} imputation will be returned as the final imputed
#' data.
#' @param preserve.zero Logical; if TRUE, perform zero-preserving low rank approximation on the \code{Sincast} imputed data.
#' (See \code{\link[Sincast]{LowRankApprox}}).
#' @param lra.control A list of controls for \code{\link[Sincast]{LowRankApprox}}.
#' @param ret.graph Logical; if TRUE, return the diffusion operator and store it in the \code{graph} slot of the imputation assay.
#' @param ndcs Calculate \code{ndcs} number of diffusion components by eigen decomposing the diffusion operator.
#'        Resulting cell embedding is stored in the \code{reduction} slot of the imputation assay.
#' @param replace Logical; if TRUE, replace the existing \code{pseudobulk} assay.
#'
#' @return A \code{Sincast} object with updated \code{imputation} assay.
#'
#' @seealso \code{\link[Sincast]{SincastAggregate}()}, \code{\link[Sincast]{LowRankApprox}()}
#'
#' @export
#' @name SincastImpute
#' @rdname SincastImpute
setGeneric("SincastImpute", function(object,
                                     assay = NULL,
                                     features = NULL,
                                     npcs = NULL,
                                     t = 3,
                                     k = 30,
                                     knn = NULL,
                                     do.umap.dist = TRUE,
                                     a = NULL,
                                     do.laplacian = TRUE,
                                     norm = c("probabilistic", "average"),
                                     do.post.scale = TRUE,
                                     preserve.zero = TRUE,
                                     lra.control = list(),
                                     ret.graph = TRUE,
                                     ndcs = 10,
                                     replace = FALSE, ...) {
  standardGeneric("SincastImpute")
})



#' @name SincastImpute
#' @rdname SincastImpute
setMethod("SincastImpute", "Sincast", function(object,
                                               assay = NULL,
                                               features = NULL,
                                               npcs = NULL,
                                               t = 3,
                                               k = 30,
                                               knn = NULL,
                                               do.umap.dist = TRUE,
                                               a = NULL,
                                               do.laplacian = TRUE,
                                               norm = c("probabilistic", "average"),
                                               do.post.scale = TRUE,
                                               preserve.zero = TRUE,
                                               lra.control = list(),
                                               ret.graph = TRUE,
                                               ndcs = 10,
                                               replace = FALSE, ...) {
  # Check the validity of the Sincast object.
  Sincast::CheckSincastObject(object, complete = FALSE, test = FALSE)

  # If a imputation assay already exists, either replace it or return an error.
  if (!is.null(
    Sincast::GetSincastAssays(object, assay = "imputation")
  )) {
    if (!replace) {
      stop(
        "SincastImpute: An imputation assay already exists, set 'replace = T' to enforce a replacement."
      )
    } else {
      message("SincastImpute: An imputation assay already exists. Will be replaced as 'replace = T'.")
    }
  }

  # Get the original Seurat object.
  original <- Sincast::GetSincastAssays(object, assay = "original")

  # Get the default assay.
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(original)
  } else {
    Seurat::DefaultAssay(original) <- assay
  }

  # Check norm format.
  norm <- match.arg(norm)

  # Get pcs.
  pcs <- Seurat::Embeddings(original, reduction = "pca")
  cells <- rownames(pcs)
  if (!is.null(npcs)) {
    pcs <- pcs[, 1:npcs]
  }
  if( Seurat::Reductions(original, "pca")@assay.used != assay ){
    stop("SincastImpute: PCA was calculated on another assay other than ", assay)
  }
  message("SincastImpute: ready to impute ", length(cells), " cells that are stored in ", assay, " assay, data layer.")

  # Subset the original Seurat object.
  original <- original[, cells]
  if (!is.null(features)) original <- original[features, ]

  # Construct affinity matrix.
  message("Now construct affinity matrix.")
  message("\t Calcualting distance.")
  dist <- pdist(pcs)

  # Adjust bandwidth.
  if (is.null(a)) {
    if (do.umap.dist) a <- log2(k) else a <- log(k / (k - 1))
  }
  if (is.null(knn)) knn <- k

  message("\t Scaling distance.")
  if (do.umap.dist) dist <- ScaleDistance(dist)

  message("\t Calculating band width.")
  sigma <- c()
  for (i in 1:length(cells)) {
    dk <- sort(dist[i, ])[2:(k + 1)]
    if (a < 1) {
      sigma[i] <- dk[k] / sqrt(-2 * log(a))
    } else if (a > 1) {
      sigma[i] <- FindSigma(dk, a, k)
    } else {
      sigma[i] <- Inf
    }
  }
  names(sigma) <- cells

  message("\t Calculating affinities.")
  aff <- exp(-0.5 * (dist / sigma)^2)
  for (i in 1:length(cells)) aff[i, ][rank(-aff[i, ]) > knn] <- 0
  aff <- Matrix::Matrix(aff)

  message("\t Symmetrization.")
  aff <- Symmetrization(aff, norm = norm)

  # Laplacian normalization
  if (do.laplacian) {
    message("\t Laplacian normalization.")
    aff <- Laplacian(aff)
  }
  message("Construct affinity matrix: Done.")

  message("Now impute.")
  message("\t Constructing diffusion operator.")
  p <- aff / rowSums(aff)

  message("\t Diffusing.")
  out <- SeuratObject::GetAssayData(
    object = original,
    layer = "data"
  )

  for (i in 1:t) {
    out <- tcrossprod(out, t(p))

  }
  out <- as.matrix(out)

  if(preserve.zero){
    message("\t Zero-preserving low rank approximation.")
    if(!is.list(lra.control)){
      warning("lra.control is not a list. Using the default control.")
      lra.control <- list()
    }
    lra.control[["x"]] <- out
    out <- do.call(Sincast::LowRankApprox, lra.control)
  }

  if(do.post.scale){
    message("\t Post-scaling.")
    if(!is.list(lra.control)){
      warning("lra.control is not a list. Using the default control.")
      lra.control <- list()
    }
    lra.control[["x"]] <- SeuratObject::GetAssayData(
      object = original,
      layer = "data"
    )
    lra.original <- do.call(Sincast::LowRankApprox, lra.control)
    out <- QQSCale(lra.original, out)
    out <- (lra.original + out)/2
  }
  message("Finish impute")

  # Get metadata.
  meta.data <- original@meta.data

  # Calculate sparsity
  sparsity.before <- mean(
    SeuratObject::GetAssayData(
      object = original,
      layer = "data"
    ) == 0
  )
  sparsity.after <- mean(out == 0)

  # Generate a token for the original data.
  if (is.null(
    Seurat::Misc(original, slot = "SincastToken")
  )) {
    SincastToken <- GenerateSincastToken()
    Seurat::Misc(
      Sincast::GetSincastAssays(object, assay = "original"),
      slot = "SincastToken"
    ) <- SincastToken
  } else {
    SincastToken <- Seurat::Misc(original, slot = "SincastToken")
  }

  # Generate a token for this imputation run.
  SincastToken <- GenerateSincastToken(
    by = "SincastImpute",
    command = deparse(match.call()),
    extend = SincastToken@command
  )

  # Create a new Seurat object to store the result.
  suppressWarnings(
    out <- Seurat::CreateSeuratObject(
      counts = Matrix::Matrix(out),
      assay = assay,
      meta.data = meta.data,
      project = SincastToken@id, ...
    )
  )
  Seurat::Idents(out) <- Seurat::Idents(original)

  # Write the summary information
  object@summary@summary["imputation", ] <- c(
    assay, "data", nrow(out), ncol(out),
    NA, round(sparsity.before, 3), round(sparsity.after, 3)
  )
  SincastToken@summary <- object@summary@summary["imputation", ]

  # Add token to the new Seurat object,
  Seurat::Misc(out, slot = "SincastToken") <- SincastToken

  # Return the diffusion operator and diffusion map.
  if (ret.graph) {
    graph.names <- paste("Sincast", assay, "DO", sep = "_")
    out[[graph.names]] <- as.Graph(p)
  }

  message("Constructing diffusion map.")
  if (!is.null(ndcs)) {
    s <- RSpectra::eigs(p, k = ndcs + 1, method = "LR")
    s$values <- Re(s$values)
    s$vectors <- Re(s$vectors)
    dc <- sweep(s$vectors, 2, s$values^t, "*")
    dc <- dc[, -1]
    rownames(dc) <- cells
    colnames(dc) <- paste("DC", 1:ndcs, sep = "_")
    out[["Sincast_dm"]] <- Seurat::CreateDimReducObject(
      embeddings = dc,
      assay = assay
    )
  }

  # Add or replace the imputation assay by the new Seurat object.
  Sincast::GetSincastAssays(object, assay = "imputation") <- out

  message(
    "SincastImpute: Sparsity before imputation: ", round(sparsity.before, 3),
    "; After imputation: ", round(sparsity.after, 3)
  )

  object
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ImputationPlot
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Plot the diffusion components computed on the diffusion operator used for \code{Sincast}
#' imputation
#'
#' To be added.
#'
#' @param object \code{Sincast}'s imputation assay, which suppose to be a \code{Seurat} object
#'  generated by \code{\link[Sincast]{SincastImpute}}.
#' @param dims Dimensions to plot, must be a three-length numeric vector specifying x-, y- and z-dimensions.
#' @param cells Cells to plot.
#' @param color.by Color query cells by which feature or meta.data attribute.
#' @param colors Color Scheme of col.by. Should be a named vector using color
#'  codes as values and labels of \code{color.by} as names.
#' @param anno.by Additional annotation of query cells.
#'
#' @return A \code{plotly} object.
#'
#' @seealso \code{\link[Sincast]{SincastAggregate}()}
#'
#' @family Sincast plot methods
#'
#' @export
#' @name ImputationPlot
#' @rdname ImputationPlot
setGeneric("ImputationPlot", function(object,
                                      dims = 1:3,
                                      cells = NULL,
                                      color.by = "ident",
                                      colors = NULL,
                                      anno.by = NULL, ...) {
  standardGeneric("ImputationPlot")
})


#' @name ImputationPlot
#' @rdname ImputationPlot
setMethod("ImputationPlot", "Seurat", function(object,
                                               dims = 1:3,
                                               cells = NULL,
                                               color.by = "ident",
                                               colors = NULL,
                                               anno.by = NULL, ...) {
  if (!all(is.numeric(dims)) | length(dims) != 3) {
    warning("dims must be a three-length numeric vector.")
  }

  if (!is.null(cells)) object <- object[cells, ]

  dcs <- Seurat::Embeddings(object, reduction = "Sincast_dm")
  dcs <- dcs[, dims]
  axis.labels <- colnames(dcs)
  colnames(dcs) <- c("x", "y", "z")


  # Generate color
  color <- NULL
  if (!is.character(color.by) | length(color.by) != 1) {
    warning("ImputationPlot: ", "color.by should be a single character.")
    color.by <- "ident"
  }

  if (color.by == "ident") {
    color <- Seurat::Idents(object)[]
  } else if (color.by %in% rownames(object)) {
    color <- suppressWarnings( Seurat::GetAssayData(object)[color.by, ] )
  } else if (color.by %in% colnames(object@meta.data)) {
    color <- object@meta.data[, color.by]
  } else {
    warning("ImputationPlot: ", color.by, " was not found in neither the data nor the metadata.")
  }

  # Generate annotation
  anno <- NULL

  if (!is.null(anno.by) & !all(is.character(anno.by))) {
    warning("ImputationPlot: ", "anno.by should be a character vector.")
    anno.by <- NULL
  }

  for (i in anno.by) {
    if (i == "ident") {
      tmp <- paste(i, Seurat::Idents(object)[], sep = ":")
      anno <- paste(anno, tmp, sep = "\n")
    } else if (i %in% rownames(object)) {
      tmp <- suppressWarnings( paste(i, Seurat::GetAssayData(object)[i, ], sep = ":") )
      anno <- paste(anno, tmp, sep = "\n")
    } else if (i %in% colnames(object@meta.data)) {
      tmp <- paste(i, object@meta.data[, i], sep = ":")
      anno <- paste(anno, tmp, sep = "\n")
    } else {
      warning("ImputationPlot: ", i, " was not found in either the data nor the metadata.")
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


#' @name ImputationPlot
#' @rdname ImputationPlot
setMethod("ImputationPlot", "Sincast", function(object,
                                                dims = 1:3,
                                                cells = NULL,
                                                color.by = "ident",
                                                colors = NULL,
                                                anno.by = NULL, ...) {
  # Check the validity of the Sincast object.
  Sincast::CheckSincastObject(object, complete = FALSE, test = FALSE)

  object <- Sincast::GetSincastAssays(object, "imputation")
  Sincast::ImputationPlot(
    object = object, dims = dims, cells = cells,
    color.by = color.by, colors = colors,
    anno.by = anno.by, ...
  )
})
