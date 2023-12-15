GenerateSincastToken <- function(by = "GenerateSincastToken", command = deparse(match.call()), extend = NULL){
    id <- paste(sample(c(rep(0:9, 3), letters, LETTERS), 16, TRUE), collapse = "")
    timestamp <- Sys.time()
    if(is.null(extend)){
      command <- list(command)
      names(command) <- id

    }else if(SeuratObject::IsNamedList(extend)){
      extend[[id]] <- command
      command <- extend

    }else{
      stop("Invalid Extension.")
    }

    new("SincastToken", id = id, timestamp = timestamp, by = by, command = command)
}

