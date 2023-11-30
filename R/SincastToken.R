GenerateSincastToken <- function(by = "GenerateSincastToken", command = deparse(match.call()), extend = NULL){
    id <- paste(sample(c(rep(0:9, 3), letters, LETTERS), 16, TRUE), collapse = "")
    timestamp <- timestamp(prefix = NULL, suffix = NULL, quiet = T)
    if(is.null(extend)){
      command <- list(command)
      names(command) <- id

    }else if(SeuratObject::IsNamedList(extend)){
      ids <- c(names(extend),id)
      command <- c(extend, command)
      names(command) <- ids

    }else{
      stop("Invalid Extension.")
    }

    new("SincastToken", id = id, timestamp = timestamp, by = by, command = command)
}

