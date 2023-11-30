# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# initialize.SincastSummary
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setMethod("initialize", "SincastSummary", function(.Object, ..., summary) {
  rname <- c(
    "pseudobulk", "imputation", "original.atlas",
    "pseudobulk.atlas", "imputation.atlas"
  )
  cname <- c(
    "assay", "layer", "nfeatures", "nsamples",
    "ncomponents", "sparsity before", "sparsity after"
  )

  .Object@summary <- data.frame(matrix(
    nrow = length(rname),
    ncol = length(cname), dimnames = list(rname, cname)
  ))

  .Object
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# show.SincastSummary
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setMethod("show", "SincastSummary", function(object) {
  cat("Summary of Sincast results\n")
  print(object@summary)
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# summary.Sincast
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setMethod("summary", "Sincast", function(object) {
  show(object@summary)
})

