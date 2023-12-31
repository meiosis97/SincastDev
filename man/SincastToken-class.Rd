% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClassDefinition.R
\docType{class}
\name{SincastToken-class}
\alias{SincastToken-class}
\alias{SincastToken}
\title{An S4 class of \code{SincastToken} object.}
\description{
To be added.
}
\section{Slots}{

\describe{
\item{\code{id}}{A 16 character unique identifier of a particular \code{Sincast} result.}

\item{\code{timestamp}}{A string recording the time at which the result was generated.}

\item{\code{by}}{A string recording the function by which this \code{Sincast} result was generated.}

\item{\code{command}}{A list recording \code{Sincast} command history.}

\item{\code{summary}}{A one row \code{data.frame} storing summary information.}
}}

\seealso{
Other Sincast classes: 
\code{\link{Sincast-class}},
\code{\link{SincastAssays-class}},
\code{\link{SincastAtlas-class}},
\code{\link{SincastSummary-class}}
}
\concept{Sincast classes}
