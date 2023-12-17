#' Check if a \code{R} variable is a numerical vector.
#'
#' To be added.
#'
#' @param x a \code{R} code
#' @param n Check if \code{x} is a numerical vector of length \code{n}.
#' @param in.which A numerical vector. Check if \code{x} is a subset of \code{in.which}.
#'
#' @return Logical indicates the test result.
#'
#' @export
#' @name IsNumericVector
#' @rdname IsNumericVector
IsNumericVector <- function(x, n = NULL, in.which = NULL){
  out <- TRUE

  if(!is.vector(x)|!is.numeirc(x)){
    out <- FALSE
  } else {
    if(!is.null(n)){
      if(!length(x) %in% n){
        out <- FALSE
      }
    }

    if(!is.null(in.which)){
      if(!all(x %in% in.which)){
        out <- FALSE
      }
    }
  }
  out
}

#' Check if a \code{R} variable is a character vector.
#'
#' To be added.
#'
#' @param x a \code{R} code
#' @param n Check if \code{x} is a character vector of length \code{n}.
#' @param in.which A character vector. Check if \code{x} is a subset of \code{in.which}.
#'
#' @return Logical indicates the test result.
#'
#' @export
#' @name IsCharacterVector
#' @rdname IsCharacterVector
IsCharacterVector <- function(x, n = NULL, in.which = NULL){
  out <- TRUE

  if(!is.vector(x)|!is.character(x)){
    out <- FALSE
  } else {
    if(!is.null(n)){
      if(!length(x) %in% n){
        out <- FALSE
      }
    }

    if(!is.null(in.which)){
      if(!all(x %in% in.which)){
        out <- FALSE
      }
    }
  }
  out
}

#' Check if a \code{R} variable is a logical vector.
#'
#' To be added.
#'
#' @param x a \code{R} code
#' @param n Check if \code{x} is a logical vector of length \code{n}.
#' @param in.which A logical vector. Check if \code{x} is a subset of \code{in.which}.
#'
#' @return Logical indicates the test result.
#'
#' @export
#' @name IsLogicalVector
#' @rdname IsLogicalVector
IsLogicalVector <- function(x, n = NULL, in.which = NULL){
  out <- TRUE

  if(!is.vector(x)|!logical(x)){
    out <- FALSE
  } else {
    if(!is.null(n)){
      if(!length(x) %in% n){
        out <- FALSE
      }
    }

    if(!is.null(in.which)){
      if(!all(x %in% in.which)){
        out <- FALSE
      }
    }
  }
  out
}

#' Perform rank transformation
#'
#' To be added.
#'
#' @param x a sample by feature matrix.
#'
#' @return A rank transformed data matrix.
#'
#' @export
#' @name RankTrans
#' @rdname RankTrans
RankTrans <- function(x){
  (apply(x,2,rank, ties.method = 'min')-1)/(nrow(x) - 1)
}
