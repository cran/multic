print.multic <- function(x, ...) {
  if(!is.multic(x)) {
    stop("\n", substitute(x), "is not a multic object\n")
  }  
  multic.obj <- x

  ## Define the significant digits amount
  significant.digits <- 5
  
  ## Print the way the multic function was called that created multic.obj
  cat("Call:\n")
  dput(multic.obj$call)
  cat("\n")
  
  ## multi/long
  if(multic.obj$metadata$repeat.count == 1) {
    cat("Multivariate analysis counts:\n")
  }else {
    cat("Longitudinal analysis counts (", multic.obj$metadata$repeat.count,
        "time points ):\n")
  }

  ## Print the pedigree information  
  cat("\n")
  cat("  Pedigree Information\n")
  cat("  --------------------\n")
  print(multic.obj$counts[1:5])

  ## Print the complete count information
  cat("\n")
  cat("  Complete Phenotype Count Information\n")
  cat("  ------------------------------------\n")
  print(multic.obj$counts[6:8])

  ## Print other count information  
  cat("\n")
  cat("  Other Information\n")
  cat("  -----------------\n")
  print(multic.obj$counts[9:11])
  cat("\n")

  if(is.null(multic.obj$metadata$mloci.out)) {
    cat("No marker data was used in this analysis to generate lod scores.\n")
  }else {
    ## Get the max log odd score and the locus at which it occured
    max.lod.and.locus <- get.max.lod.and.locus(multic.obj)
    max.lod.score <- max.lod.and.locus$max.lod.score
    max.lod.locus <- max.lod.and.locus$max.lod.locus
    max.lod.centimorgan <- max.lod.and.locus$max.lod.centimorgan
    
    multiple.max <- FALSE
    if(length(max.lod.locus) > 1) {
      max.lod.locus <- max.lod.locus[1]
      multiple.max <- TRUE
    }
    if(length(max.lod.centimorgan) > 1) {
      max.lod.centimorgan <- max.lod.centimorgan[1]
      multiple.max <- TRUE
    }
    
    cat("Maximum lod score:", signif(max.lod.score, significant.digits),
        "\n         at locus:", max.lod.locus,
        "\nat positions (cM):", max.lod.centimorgan, "\n")
    if(multiple.max) {
      cat('Note: Multiple loci generated the same maximum lod score.\n')
    }
  }
  
  invisible()
}
