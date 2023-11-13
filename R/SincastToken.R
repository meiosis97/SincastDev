# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 class definition
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SincastToken <- setClass(
  Class = "SincastToken",
  slots = list(
    id = "character",
    timestamp = "character"
  )
)

GenerateSincastToken <- function(){
  id <- paste(sample(c(rep(0:9, 3), letters, LETTERS), 16, TRUE), collapse = "")
  timestamp <- timestamp(prefix = NULL, suffix = NULL, quiet = T)
  new('SincastToken', id = id, timestamp = timestamp)
}
