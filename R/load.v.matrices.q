load.v.matrices <- function(family.names) {
  v.matrix.files <- multic.system("ls v.matrix.*.log")

  if(length(v.matrix.files) > 0) {
    null.index <- match("v.matrix.null.log", v.matrix.files)
    v.matrix.files <- c(v.matrix.files[null.index],
                        v.matrix.files[-null.index])
    
    family.count <- length(family.names)
    v.matrices <- sapply(v.matrix.files, load.v.matrix.file, family.count)
    
    ## Generate the loci name labels
    if(length(v.matrix.files) > 1) {
      loci.names <- matrix(multic.strsplit(v.matrix.files[-1], sep = "."),
                           byrow = T, nrow = length(v.matrix.files) - 1)
      loci.names <- loci.names[, -c(1, 2, ncol(loci.names))]
      loci.names <- apply(loci.names, 1, paste, collapse = '.')
    }else {
      loci.names <- NULL
    }
    
    dimnames(v.matrices) <- list(family.names, c("null", loci.names))
    
    return (v.matrices)
  }

  return ("No V matrix data was calculated for this multic object.")
}

load.v.matrix.file <- function(file.name, family.count) {
  v.matrix.file <- .Call("loadVMatrixFile",
                         as.character(file.name),
                         as.integer(family.count),
                         PACKAGE = "multic"
                         )
  
  return (v.matrix.file)
}


