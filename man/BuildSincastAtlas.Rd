% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SincastAtlas.R
\name{BuildSincastAtlas}
\alias{BuildSincastAtlas}
\alias{BuildSincastAtlas,Sincast-method}
\title{Build a \code{Sincast} atlas, which is a \code{Seurat} object with a rank-transformed reference
data, as well a pca dimensional reduction performed on the data.}
\usage{
BuildSincastAtlas(
  object,
  cells = NULL,
  sincast.assay = NULL,
  seurat.assay = NULL,
  layer = "counts",
  reference.features = NULL,
  query.features = NULL,
  ScaleData.control = list(),
  RunPCA.control = list(),
  replace = FALSE,
  ...
)

\S4method{BuildSincastAtlas}{Sincast}(
  object,
  cells = NULL,
  sincast.assay = NULL,
  seurat.assay = NULL,
  layer = "counts",
  reference.features = NULL,
  query.features = NULL,
  ScaleData.control = list(),
  RunPCA.control = list(),
  replace = FALSE,
  ...
)
}
\arguments{
\item{object}{A \code{Sincast} object.}

\item{cells}{Cells to be included in the atlas. Can be indexes, cell names or logical.}

\item{sincast.assay}{Either one of the "pseudobulk", "imputation" or "original" indicating the
type of atlas to build. Default is the active \code{Sincast} assay (Can be check by calling the \code{Sincast} object).}

\item{seurat.assay}{Name of the \code{Seurat} assay rank transformation, data scaling and PCA are being run on.}

\item{layer}{Name of the \code{Seurat} layer rank transformation, data scaling and PCA are being run on.}

\item{reference.features}{Features on which the atlas is being built. Can be indexes, feature names, logical
or a single character corresponding to a logical attribute in \code{Seurat} assay's feature metadata. If the data
does not contain any of the given reference features, the missing features will be set to zero in the atlas construction.}

\item{query.features}{Must be a character vector of feature names. If provided, the atlas will be build on the
intersection between reference.features and query.features.}

\item{ScaleData.control}{A list of arguments to the \code{Seurat} \code{\link[Seurat]{ScaleData}} function call.}

\item{RunPCA.control}{A list of arguments to the \code{Seurat} \code{\link[Seurat]{RunPCA}} function call.}

\item{replace}{Logical; if TRUE, replace the existing atlas corresponding to the \code{sincast.assay}.}
}
\value{
A \code{Sincast} object with updated \code{Sincast} atlas.
}
\description{
To be added.
}
\seealso{
Other SincastAtlas related methods: 
\code{\link{CheckSincastAtlas}()},
\code{\link{CleanSincastAtlas}()},
\code{\link{GetSincastAtlas}()}
}
\concept{SincastAtlas related methods}
