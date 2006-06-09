print.summary.multic <- function(x, ...) {
  if(oldClass(x) != "summary.multic") {
    stop("\n", substitute(x), "is not a summary.multic object\n")
  }
  summary.multic.object <- x
  
  significant.digits <- 5

  cat("Call:\n")
  print(summary.multic.object$call)
  cat("\n")

  if(is.null(summary.multic.object$call$mloci.out)) {
    cat("No marker data was used in this analysis to generate lod scores.\n")
  }else {
    max.lod.score <- summary.multic.object$max.lod.score
    max.lod.locus <- summary.multic.object$max.lod.locus
    max.lod.centimorgan <- summary.multic.object$max.lod.centimorgan
    
    ## Print the max lod score, its locus, and its centimorgan distance
    cat("Maximum lod score:", signif(max.lod.score, significant.digits),
        "\n         at locus:", max.lod.locus,
        "\nat positions (cM):", max.lod.centimorgan, "\n\n")
    
    ## Print the top n families at the max lod score
    if(!is.null(summary.multic.object$top.n.families)) {
      cat("The top", summary.multic.object$n, "families and their",
          "lod score contributions","\nto the maximum lod score are:\n")
      top.n.families <- summary.multic.object$top.n.families
      top.n.families[, 2] <- signif(top.n.families[, 2], significant.digits)
      print(top.n.families[, 2])
      cat("\n")
    } else {
      cat("Since multic was run with calc.fam.log.liks = F (default),",
          "\nthe top families and their lod scores have not been calculated.",
          "\n\n")
    }
    
    ## Print the centimorgan that were contiguous to the peak and produced lod
    ## scores at least the maximum lod score - 1
    cat("The minimum and maximum positions (cM) that produced a lod score",
        "\ngreater than the maximum - 1 (",
        signif(max.lod.score, significant.digits),
        "- 1 ) \nand are contiguous to", max.lod.locus, "are:\n")
    if(is.null(summary.multic.object$centimorgan.close.to.peak)) {
      print(NULL)
    }else {
      range.count <- length(summary.multic.object$centimorgan.close.to.peak)
      cat(c(summary.multic.object$centimorgan.close.to.peak[1],
            summary.multic.object$centimorgan.close.to.peak[range.count]))
    }
    cat("\n")
  }
}
