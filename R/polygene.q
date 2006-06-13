polygene <- function(multic.obj) {
  significant.digits <- 4

  ## Print the pedigree information  
  cat("\n")
  cat("  Pedigree Information\n")
  cat("  --------------------\n")
  print(multic.obj$counts[1:5])

  ## Print the complete count information
  cat("\n")
  cat("  Complete-Count Information\n")
  cat("  --------------------------\n")
  print(multic.obj$counts[6:8])

  ## Print other count information  
  cat("\n")
  cat("  Other Information\n")
  cat("  -----------------\n")
  print(multic.obj$counts[9:11])

  ## Print the mean, std err, min, and max for the traits and covariates
  cat("\n")
  cat("  Descriptive Statistics for the Variables\n")
  cat("  ----------------------------------------\n")
  print(signif(multic.obj$descriptives, significant.digits))

  ## Create local variables holding the number and names of the traits and
  ## covariates for use in printing the fixed and random effects
  trait.count <- multic.obj$metadata$trait.count
  covariate.count <- multic.obj$metadata$covariate.count
  trait.names <- multic.obj$metadata$trait.names
  covariate.names <- multic.obj$metadata$covariate.names

  ## Print the fixed effects under the null hypothesis
  cat("\n")
  cat("  Covariate coefficients\n")
  cat("  ----------------------\n")
  trait.mean.index <- 1

  polygene.fixed.effects <-
    matrix(multic.obj$fixed.effects[, , 1],
           nrow = trait.count * (1 + covariate.count),
           dimnames = dimnames(multic.obj$fixed.effects)[-3])
  for(i in 1:trait.count) {
    cat("Trait", i, "(", trait.names[i], "):\n")
    ## It doesn't make sense to list the trait name as a covariate
    ## coefficient, so change the name here to "intercept"
    dimnames(polygene.fixed.effects)[[1]][trait.mean.index] <- "(Intercept)"
    print(signif(polygene.fixed.effects[seq(from=trait.mean.index,
                                            length=covariate.count + 1), ,
                                        drop = FALSE],
                 significant.digits)
          )
    trait.mean.index <- trait.mean.index + covariate.count + 1
    cat("\n")
  }

  ## Print the random effects under the null hypothesis
  cat("  Variance Components\n")
  cat("  -------------------\n")
  cat("Polygenic:\n")
  print(signif(multic.obj$polygenic[, , 1]), significant.digits)
  cat("\n")
  
  cat("Environmental (non-genetic component):\n")
  print(signif(multic.obj$environmental[, , 1], significant.digits))
  cat("\n")

  cat("  Proportion of Variance due to the Covariates\n")
  cat("  --------------------------------------------\n")
  cat("R.sq: ")
  R.sq <- multic.obj$R.sq
  if(length(R.sq) == 1 && is.numeric(R.sq)) {
    cat(signif(R.sq, significant.digits), "\n", sep = "")
  } else {
    R.sq[2] <- paste("      ", R.sq[2], sep = "")
    cat(R.sq, sep = "\n")
  }

  invisible()  
}
