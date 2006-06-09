load.y.beta.diffs <- function(family.names) {
  y.beta.diff.files <- multic.system("ls y.beta.diff.*.log")

  if(length(y.beta.diff.files) > 0) {
    null.index <- match("y.beta.diff.null.log", y.beta.diff.files)
    y.beta.diff.files <- c(y.beta.diff.files[null.index],
                           y.beta.diff.files[-null.index])
    
    family.count <- length(family.names)
    y.beta.diffs <- sapply(y.beta.diff.files, load.y.beta.file, family.count)

    ## Generate the loci name labels
    if(length(y.beta.diff.files) > 1) {
      loci.names <- matrix(multic.strsplit(y.beta.diff.files[-1], sep = "."),
                           byrow = T, nrow = length(y.beta.diff.files) - 1)
      loci.names <- loci.names[, -c(1, 2, 3, ncol(loci.names))]
      loci.names <- apply(loci.names, 1, paste, collapse = '.')
    }else {
      loci.names <- NULL
    }
      
    dimnames(y.beta.diffs) <- list(family.names, c("null", loci.names))
    
    return (y.beta.diffs)
  }

  return ("No Y-beta data was calculated for this multic object.")
}

load.y.beta.file <- function(file.name, family.count) {
  y.beta.diff.file <- .Call("loadYBetaDiffFile",
                            as.character(file.name),
                            as.integer(family.count),
                            PACKAGE = "multic"
                            )
  
  return (y.beta.diff.file)
}
