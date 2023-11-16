GenerateSincastToken <- function(by = "", command = list()){
  id <- paste(sample(c(rep(0:9, 3), letters, LETTERS), 16, TRUE), collapse = "")
  timestamp <- timestamp(prefix = NULL, suffix = NULL, quiet = T)
  new('SincastToken', id = id, timestamp = timestamp, by = by, command = command)
}

