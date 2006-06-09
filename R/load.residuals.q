#############################################################################
## Title: load.residuals
## Description: function to read the v.matrix.*.log and y.beta.diff.*.log
##              files and generate the residuals
## Author: Eric Lunde, 9-21-2004
## Updates:
## 9-21-2004, Eric Lunde, When running only the polygene model, there are no
##            v.matrix.*.log nor y.beta.diff.*.log files (except for the null
##            files).  I've added logic to use the alternative file names
##            only when we are also calculating the sporadic models as well.
#############################################################################
load.residuals <- function(famid, id, trait.count, residual.placement) {
  ## browser()
  family.sizes <- table(famid)
  unique.famid <- unique(famid)
  family.sizes <- family.sizes[match(unique.famid, names(family.sizes))]
  
  families.count <- length(family.sizes)
  people.count <- length(id)
  
  ## Get the names of the y.beta.diff.*.log files and put the
  ## y.beta.diff.null.log first
  y.betas <- multic.system('ls y.beta.diff.*.log')
  if(length(y.betas) == 0) {
    cat(paste("- No y.beta.diff.*.log files were found\n",
              "- No residual information will be loaded.\n"))
    return (NULL)
  } else {
    y.betas <- c(y.betas[length(y.betas)],
                 y.betas[-length(y.betas)])
  }
  
  ## Get the names of the v.matrix.*.log files and put the
  ## v.matrix.null.log first
  v.mats <- multic.system('ls v.matrix.*.log')
  if(length(v.mats) == 0) {
    cat(paste("- No v.matrix.*.log files were found\n",
              "- No residual information will be loaded.\n", sep = ""))
    return (NULL)
  } else {
    v.mats <- c(v.mats[length(v.mats)],
                v.mats[-length(v.mats)])
  }

  if(length(v.mats) != length(y.betas)) {
    stop("length(v.mats) (", length(v.mats),
         ") != length(y.betas) (", length(y.betas), ")\n",
         "This is an error.  load.residuals.q key 13")
  }

  number.of.file.pairs <- length(v.mats)

  ## missing.value is used to determine which residual locations were not
  ## written to.
  missing.value <- -9.1
  marginal.resids.input <- rep(-9.1, people.count * number.of.file.pairs
                               * trait.count)
  overall.resids.input <- rep(-9.1, people.count * number.of.file.pairs
                              * trait.count)
  
  load.residuals.result <-
    .C("loadResiduals",
       as.character(y.betas),
       as.character(v.mats),
       as.integer(number.of.file.pairs),
       as.integer(family.sizes),
       as.integer(families.count),
       marginal.resids = marginal.resids.input,
       ##rep(-9.1, people.count * number.of.file.pairs * trait.count),
       ##       numeric(people.count * number.of.file.pairs * trait.count),
       ##       marginal.resids = input.marg,#numeric(sum(residual.placement))
       overall.resids = overall.resids.input,
       ##rep(-9.1, people.count * number.of.file.pairs * trait.count),
       ##       numeric(people.count * number.of.file.pairs * trait.count),
       ##       overall.resids = numeric(sum(residual.placement)),
       as.integer(people.count * number.of.file.pairs),
       PACKAGE = "multic")
 
  marginal.resid <- array(load.residuals.result$marginal.resids,
                          dim = c(people.count * number.of.file.pairs,
                            trait.count))
  overall.resid <- array(load.residuals.result$overall.resids,
                         dim = c(people.count * number.of.file.pairs,
                           trait.count))
  
  ## Combine the marginal and overall into one array
  residuals <- array(0, dim = c(people.count * number.of.file.pairs,
                          trait.count * 2))
  residuals[, 2 * (1:trait.count - 1) + 1] <- marginal.resid
  residuals[, 2 * 1:trait.count] <- overall.resid
  residuals.dimnames <- character(trait.count * 2)
  residuals.dimnames[2 * (1:trait.count - 1) + 1] <-
    paste("marg", 1:trait.count, sep = "")
  residuals.dimnames[2 * 1:trait.count] <-
    paste("over", 1:trait.count, sep = "")
  dimnames(residuals) <- list(NULL, residuals.dimnames)

  residuals.assigned <- apply(residuals != missing.value, 2, sum)
  equal.to.first.value <- residuals.assigned == residuals.assigned[1]
  if( !all(equal.to.first.value) ) {
    warning("The number of marginal and overall residuals calculated is not",
            " consistant.\nThe results may be inaccurate.\n",
            "load.residuals.q key 103\n")
  }

  ## Remove rows based on whether the first column has a missing value or not
  valued.residuals <- residuals[residuals[,1] != missing.value, ]

  ## If there were errors in the residual calculation and some missing values
  ## (-9.1) were not removed, make them NA's now
  valued.residuals[valued.residuals == missing.value] <- NA

  ## Replace the residuals in the correct order
  residuals[residual.placement, ] <- valued.residuals
  residuals[!residual.placement, ] <- NA
  
  ## Get the marker names from the y.beta.diff.*.log files
  null.y.beta <- y.betas[1]
  null.y.beta.split <- multic.strsplit(null.y.beta, sep = '.')
  null.y.beta.name <-
    null.y.beta.split[-match(c('y', 'beta', 'diff', 'log'),
                             null.y.beta.split)]

  ## Commented out on 9-21-2004 by Eric Lunde, the functionality was replaced
  ## by get.cM.from.v.mat
  ##  alternative.y.betas <- y.betas[-1]
  ##  if(length(alternative.y.betas) != 0) {
  ##    alternative.y.betas.split <-
  ##      apply(array(alternative.y.betas), 1, multic.strsplit, sep = '.')
  ##    alternative.y.betas.names <-
  ##      apply(alternative.y.betas.split[-match(c('y', 'beta', 'diff', 'log'),
  ##                                             alternative.y.betas.split), ],
  ##            2, paste, collapse = '.')
  ##  }
  
  ## This section of code creates the famid column of the residual
  ## data.frame.  We need to repeat each individual's family id 
  full.famid <- array(famid, dim = c(length(famid), number.of.file.pairs))
  full.famid <- as.vector(t(full.famid))
  
  cM <- c(null.y.beta.name, get.cM.from.v.mat(v.mats[-1]))
  repeat.values <- as.vector(apply(array(family.sizes), 1, rep, length(cM)))
  full.cM <- rep(rep(cM, length(unique.famid)), repeat.values)

  id.max.length <- 32
  full.id <- rep(paste(rep(" ", id.max.length), collapse = ""),
                 length(id) * number.of.file.pairs)

  full.id.results <- .C("fullId",
                        as.character(id),
                        as.integer(length(id)),
                        as.integer(family.sizes),
                        as.integer(length(family.sizes)),
                        as.integer(number.of.file.pairs),
                        full.id = as.character(full.id),
                        ## full.id = integer(length(id) * number.of.file.pairs),
                        as.integer(length(id) * number.of.file.pairs),
                        as.integer(id.max.length),
                        PACKAGE = "multic"
                        )
  
  ##  full.id <- full.id.results$full.id

  residuals <- data.frame(famid = full.famid, cM = full.cM,
                          id = full.id.results$full.id, round(residuals, 4))

  return (residuals)
}

get.cM.from.v.mat <- function(v.mats) {
  if(length(v.mats) == 0) {
    return (NULL)
  }else if(length(v.mats) > 1) {
    v.mat <- v.mats[1]
  }else {
    v.mat <- v.mats
  }
  
  ##  browser()

  tokens <- multic.strsplit(v.mat, sep = '.')
  if(substring(v.mat, 9, 13) == '.ibd.') {
    if(length(tokens) == 5) {
      desired.tokens <- 4
    }else {
      stop("\nThe v.matrix file '", v.mat, "', does not have 5 tokens.\n",
           "length(tokens) = ", length(tokens), ".\n load.residuals.q key 74")
    }
  }else {  
    if(length(tokens) == 6) {
      desired.tokens <- 5
    }else if(length(tokens) == 7) {
      desired.tokens <- c(5, 6)
    }else {
      stop("\nThe v.matrix file '", v.mat, "', does not have 6 nor 7 tokens.\n",
           "length(tokens) = ", length(tokens), ".\n load.residuals.q key 75")
    }
  }

  all.tokens <- apply(array(v.mats), 1, multic.strsplit, sep = '.')
  tokens.to.be.pasted <- matrix(all.tokens[desired.tokens, ],
                                nrow = length(desired.tokens))
  pasted.tokens <- apply(tokens.to.be.pasted, 2, paste,
                         collapse = '.')
  
  return (pasted.tokens)
}
