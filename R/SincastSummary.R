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
    "n.components", "var.explained", "sparsity.before", "sparsity.after"
  )

  .Object@summary <- data.frame(matrix(
    nrow = length(rname),
    ncol = length(cname), dimnames = list(rname, cname)
  ))

  .Object@active.assay <- new("character")
  .Object@active.atlas <- new("character")

  .Object
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# show.SincastSummary
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setMethod("show", "SincastSummary", function(object) {
  cat("Summary of Sincast results\n")
  summary <- object@summary
  summary[is.na(summary)] <- "-"
  print(summary)
})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# summary.Sincast
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setMethod("summary", "Sincast", function(object) {
  object@summary
})

