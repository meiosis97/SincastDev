% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SincastAtlas.R
\name{CleanSincastAtlas}
\alias{CleanSincastAtlas}
\alias{CleanSincastAtlas,SincastAtlas-method}
\alias{CleanSincastAtlas,Sincast-method}
\title{Remove unwanted or invalid (i.e, non-\code{Seurat} object) \code{Sincast} atlas
(\code{pseudobulk} and \code{imputation}) embedded in \code{Sincast}.}
\usage{
CleanSincastAtlas(
  object,
  clean.up = c("none", "all", "original", "pseudobulk", "imputation"),
  remove.invalid = TRUE,
  ...
)

\S4method{CleanSincastAtlas}{SincastAtlas}(
  object,
  clean.up = c("none", "all", "original", "pseudobulk", "imputation"),
  remove.invalid = TRUE,
  ...
)

\S4method{CleanSincastAtlas}{Sincast}(
  object,
  clean.up = c("none", "all", "original", "pseudobulk", "imputation"),
  remove.invalid = TRUE,
  ...
)
}
\arguments{
\item{object}{A \code{Sincast} object.}

\item{clean.up}{Either \code{"none"}, \code{"all"}, \code{original},
\code{"pseudobulk"} or \code{"imputation"} indicating whether to clean up
the \code{SincastAtlas} object or clean up a specific \code{Sincast} atlas as proposed.}

\item{remove.invalid}{Logical; if TRUE, check the validity of the existing \code{Sincast} atlas and
remove invalid ones (e.g, \code{Seurat} objects not generated by \code{Sincast} functions).}
}
\value{
A \code{Sincast} object with updated \code{Sincast} atlas.
}
\description{
To be added.
}
\seealso{
Other SincastAtlas related methods: 
\code{\link{BuildSincastAtlas}()},
\code{\link{CheckSincastAtlas}()},
\code{\link{GetSincastAtlas}()}
}
\concept{SincastAtlas related methods}
