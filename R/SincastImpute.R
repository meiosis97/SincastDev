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


MedianScale <- function(X, Y) {
  scale.factor <- apply(X, 1, function(x) median(x[x != 0])) /
      apply(replace(Y, X == 0, NA), 1, function(y) median(y, na.rm = T))

  Y <- Y * scale.factor
  Y[is.na(Y)] <- 0
  Y
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

PostScale <- function(before, after){
  after <- as.matrix(after)
  N <- ncol(before)
  G <- nrow(before)
  q <- ppoints(N) %>% qnorm()
  w <- rowMeans(before!=0)

  # Gene-wise mean and variance estimate on expressed genes.
  message("\t Genewise Mean and Variance estimation on imputed data.")
  mu <-c()
  v <- c()
  for(i in 1:G){
    qx <- sort(after[i,])
    if(sum(qx > 0)){
      lmod <- lm(qx~q, subset = qx > 0)
      mu[i] <- lmod$coefficients[1]
      v[i] <- lmod$coefficients[2]^2
    }else{
      mu[i] <- 0
      v[i] <- 0
    }
  }
  param <- data.frame(mu = mu, v = v, w = w)
  message("\t Done")

  # Estimate mean and variance trend.
  param$log.mu <- log(param$mu + min(param$mu[param$mu>0]))
  param$log.v <- log(param$mu + min(param$v[param$v>0]))

  message("\t Now perform GAM fit.")
  k <- 3
  gam.mod <- mgcv::gam(log.v~s(log.mu, k = k, bs = 'cr'), weights = w, data = param)
  # while(mgcv::k.check(gam.mod)[,4] < 0.05){
  #   k <- k + 2
  #   gam.mod <- mgcv::gam(log.v~s(log.mu, k = k, bs = 'cr'), weights = w, data = param)
  # }
  message("\t Finish regress. The basis dimension is ", k)

  p <- ggplot(data = param) + geom_point(aes(log.mu, log.v, col = w))+
    scale_color_continuous('Zero proportion') +
    geom_path(aes(log.mu, predict(gam.mod, data.frame(log.mu))),size = 1.5, linetype = 'dashed') +
    theme_bw() + xlab('Log-Mean') + ylab('Log-Variance') + theme(text = element_text(size=15))
  print(p)

  message("\t Perform observation-wise variance estimation.")
  sigma <- log(after +  min(param$mu[param$mu>0]))
  for(i in 1:N){
    sigma[,i] <- predict(gam.mod, data.frame(log.mu = sigma[,i])) %>% exp() -  min(param$v[param$v>0])
  }
  message("\t Done")

  # Shrink back
  e <- (before-after)^2
  lambda<- colMeans(e/(e+sigma),na.rm = T)

  sweep(before,2,lambda,"*") +  sweep(after,2,1-lambda,"*")

}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported
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
#' @param ret.graph Logical; if TRUE, return the diffusion operator and store it in the \code{graph} slot of the imputation assay.
#' @param ndcs Calculate \code{ndcs} number of diffusion components by eigen decomposing the diffusion operator.
#'        Resulting cell embedding is stored in the \code{reduction} slot of the imputation assay.
#' @param replace Logical; if TRUE, replace the existing \code{pseudobulk} assay.
#'
#' @return A \code{Sincast} object with updated \code{imputation} assay.
#'
#' @seealso [SincastAggregate()]
#'
#' @export
#' @name SincastImpute
#' @rdname SincastImpute
#' @aliases Sincast, SincastAssays
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
  message("Ready to impute ", length(cells), " cells that are stored in ", assay, " pca.")

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

  message("\t Calculating band width")
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

  message("\t Scaling.")
  out <- MedianScale(
    SeuratObject::GetAssayData(
      object = original,
      layer = "data"
    )  %>% as.matrix(), out
  )
  message("Finish impute")

  message("Post-scaling.")
  # out <- PostScale(
  #   SeuratObject::GetAssayData(
  #     object = original,
  #     layer = "data"
  #   ), out
  # )
  message("Finish post-scaling.")

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
    graph.names <- paste(assay, "Sincast", sep = "_")
    out[[graph.names]] <- as.Graph(p)
  }

  message("Constructing diffusion map.")
  if (!is.null(ndcs)) {
    s <- eigs(p, k = ndcs + 1, method = "LR")
    s$values <- Re(s$values)
    s$vectors <- Re(s$vectors)
    dc <- sweep(s$vectors, 2, s$values^t, "*")
    dc <- dc[, -1]
    rownames(dc) <- cells
    colnames(dc) <- paste("DC", 1:ndcs, sep = "_")
    out[["dm_Sincast"]] <- Seurat::CreateDimReducObject(
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
#' @seealso [SincastImpute()]
#'
#' @family Sincast plot methods
#'
#' @export
#' @name ImputationPlot
#' @rdname ImputationPlot
#' @aliases Sincast, SincastAssays, Seurat
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

  dcs <- Seurat::Embeddings(object, reduction = "dm_Sincast")
  dcs <- dcs[, dims]
  axis.labels <- colnames(dcs)
  colnames(dcs) <- c("x", "y", "z")


  # Generate color
  color <- NULL
  if (!is.character(color.by) | length(color.by) != 1) {
    warning("color.by should be a single character.")
  }

  if (color.by == "ident") {
    color <- Seurat::Idents(object)[]
  } else if (color.by %in% rownames(object)) {
    color <- suppressWarnings( Seurat::GetAssayData(object)[color.by, ] )
  } else if (color.by %in% colnames(object@meta.data)) {
    color <- object@meta.data[, color.by]
  } else {
    warning(color.by, " was not found in neither the data nor the metadata.")
  }

  # Generate annotation
  anno <- NULL

  if (!is.null(anno.by) & !all(is.character(anno.by))) {
    stop("anno.by should be a character vector.")
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
      warning(i, " was not found in either the data nor the metadata.")
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
