##############################################################################
## Title: load.effects
## Description: load.effects is responsible for allocating the memory needed
##              by a compiled C++ function, calling that function to fill that
##              memory with the values from the multic output file
##              'summary.log', modifying the data structure's dimensions and
##              names, and then finally returning the data in a list object.
## Input: fixed.effects.count - fixed.effects.count is a scalar value
##                              representing the number of fixed effects that
##                              the last run of multic output.
##        random.effects.count - random.effects.count is a scalar value
##                               representing the number of random effects per
##                               effect (polygene, major gene) that the last
##                               run of multic output.
##        locus.count - locus.count is the number of loci that the last run of
##                      mutlic calculated with (the number of alternative
##                      hypotheses plus the null hypothesis)
##        individual.count - the number of total individuals in this analysis.
##                           This is used for the calculating the degrees of
##                           freedom for the p-value for the fixed effects
##        covariate.count - the number of covariates.  This is used for the
##                          calculating the degrees of freedom for the p-value
##                          for the fixed effects
## Ouput: A list object containing the values from summary.log concerning the
##        fixed effects, polygenic, major.gene1, environmental effects, and
##        log likelihood data.
## Side Effects: Since the C++ function is hardcoded to read 'summary.log' the
##               program will crash if that file is not in the same directory
##               as the Splus session.  If the file fails to open, an
##               appropriate (and hopefully, helpful) error message will be
##               printed to standard error.
## Author: Eric Lunde 10-26-03
## Updates: (Date, Modified By, Modification Description)
## 11-28-03, Eric Lunde, 11-28-03, I removed the parameters trait.count and
##           covariate.count and replaced them with fixed.effects.count and
##           random.effects.count.  I made a decision that these effect count
##           values were the important (useful) values.  So they are
##           calculated early in the program and the values are passed through
##           rather than each function calculating the values again and again.
## 02-20-2004, EricLunde, I moved much of the additional calculations (like
##             t.value, p.value, Wald score, etc) from multic into this
##             function.  This helps to make the multic function more
##             logically designed.
## 09-21-2004, Eric Lunde, I added the parameter 'run.alternative.hyps' to
##             determine if 'load.effects' was being called from polygene or
##             multic.  This way, we only load "distances" during a multic
##             call and not a polygene call.
## 03-31-2005, Eric Lunde, The number of environmental variance (covariance)
##             values no longer is dependent upon whether multic is run in
##             multivariate or longitudinal mode (it now always operates as if
##             in multivariate mode).  Therefore I've removed the
##             environment.count variable and replaced it with the
##             random.effects.count which now always is the correct value.
##############################################################################
load.effects <- function(fixed.effects.count, random.effects.count,
                         locus.count, individual.count,
                         trait.count, covariate.count, repeat.count,
                         run.alternative.hyps) {
  ## If we do not have access to summary.log (we cannot open it), set effects
  ## to null. This is a very important step.  If we cannot access summary.log,
  ## we do not acquire the names of the effects and the heart of the multic
  ## output.
  if( !file.exists("summary.log")) {
    cat(paste("\nThe multic output file 'summary.log' is not present.\n",
              "The multic object will not contain data concerning the fixed\n",
              "effects, polygenic, major gene, environmental, or log ",
              "likelihood effects.\n",
              "Also, the names of the effects will not be present.\n",
              sep=""))
    return(NULL)
  }
  
  ## Verify that the arguments are of the correct type and range
  fixed.effects.count <- as.integer(fixed.effects.count)
  if(fixed.effects.count < 0) {
    stop(paste("fixed.effects.count has the value '", fixed.effects.count,
               "', but must be positive.  load.effects.s key 39"))
  }
  
  random.effects.count <- as.integer(random.effects.count)
  if(random.effects.count < 0) {
    stop(paste("random.effects.count has the value '", random.effects.count,
               "', but must be positive.  load.effects.s key 45"))
  }
  
  locus.count <- as.integer(locus.count)  
  if(locus.count < 0) {
    stop(paste("locus.count as the value '", locus.count, "', but must be",
               "positive.  load.effects.s key 57"))
  }

  ## Create the space for the data from the file 'summary.log' to be placed
  loci.names <- rep("                        ", locus.count)
  variable.names <- rep("                        ",
                        fixed.effects.count + (6 * random.effects.count))
  
  fixed.effects <- numeric(fixed.effects.count * 2 * locus.count)
  
  polygenic <- numeric(random.effects.count * 2 * locus.count)  
  major.gene1 <- numeric(random.effects.count * 2 * locus.count)
  environmental <- numeric(random.effects.count * 2 * locus.count)
  sibling.sibling <- numeric(random.effects.count * 2 * locus.count)
  parent.parent <- numeric(random.effects.count * 2 * locus.count)
  parent.offspring <- numeric(random.effects.count * 2 * locus.count)
  
  log.likelihoods <- numeric(locus.count)
  log.lik.status <- rep("                        ", locus.count)

  wald.data.names <- c("Wald", "W.p.value")
  herit.data.names <- c("h^2", "se.h^2", "h.p.value")

  ## Call the C++ routine to load the effects from 'summary.log'
  effects <- .C("loadEffects",
                loci.names=as.character(loci.names),
                variable.names=as.character(variable.names),
                fixed.effects=as.numeric(fixed.effects),
                polygenic=as.numeric(polygenic),
                major.gene1=as.numeric(major.gene1),
                environmental=as.numeric(environmental),
                sibling.sibling=as.numeric(sibling.sibling),
                parent.parent=as.numeric(parent.parent),
                parent.offspring=as.numeric(parent.offspring),
                log.likelihoods=as.numeric(log.likelihoods),
                log.lik.status=as.character(log.lik.status),
                fixed.effects.count=as.integer(fixed.effects.count),
                random.effects.count=as.integer(random.effects.count),
                locus.count=as.integer(locus.count),
                PACKAGE = "multic")

  loci.names <- effects$loci.names
  variable.names <- effects$variable.names
  fixed.effects.names <- variable.names[1:fixed.effects.count]

  fixed.effects <- effects$fixed.effects
  fixed.effects <- array(fixed.effects, c(fixed.effects.count, 2, locus.count))
  dimnames(fixed.effects) <- list(fixed.effects.names,
                                  c("Estimate", "Std.err"), loci.names)

  ## Add the t and p-values to the fixed.effects
  temp.fixed.effects <- array(NA, dim=c(dim(fixed.effects)[1],
                                   dim(fixed.effects)[2] + 2,
                                   dim(fixed.effects)[3]),
                              dimnames=list(dimnames(fixed.effects)[[1]],
                                c(dimnames(fixed.effects)[[2]],
                                  "t.value", "p.value"),
                                dimnames(fixed.effects)[[3]]))
  temp.fixed.effects[, 1:2, ] <- fixed.effects
  fixed.effects <- temp.fixed.effects

  ## Don't round p and t values per Mariza request August 2015
  fixed.effects[, "t.value", ] <- fixed.effects[, 1, ] / fixed.effects[, 2, ]
#    round(fixed.effects[, 1, ] / fixed.effects[, 2, ], 4)
  fixed.effects[, "p.value", ] <- 2 * (1 - pt(abs(fixed.effects[, 3, ]),
                      individual.count - (1 + covariate.count)))
#    round(2 * (1 - pt(abs(fixed.effects[, 3, ]),
#                      individual.count - (1 + covariate.count))), 4)
  
  ## Prepare the polygenic data
  polygenic <- effects$polygenic
  start.index <- fixed.effects.count + 1
  end.index <- start.index + random.effects.count - 1
  polygenic.names <- variable.names[start.index:end.index]
  polygenic <- array(polygenic, c(random.effects.count, 2, locus.count))
  dimnames(polygenic) <- list(polygenic.names,
                              c("Estimate", "Std.err"), loci.names)
  ## 69 and 73 are the integers chosen to represent a fixed externally and a
  ## fixed internally (respectively) standard error numerically
  polygenic[, 2, ][polygenic[, 2, ] == 69] <- NA
  polygenic[, 2, ][polygenic[, 2, ] == 73] <- NA
  
  ## Add Wald test, p-value, h2 and se.h2 to the polygenic compenents
  temp.polygenic <- array(NA, dim=c(dim(polygenic)[1],
                                dim(polygenic)[2] + 5,
                                dim(polygenic)[3]),
                          dimnames=list(dimnames(polygenic)[[1]],
                            c(dimnames(polygenic)[[2]],
                              wald.data.names, herit.data.names),
                            dimnames(polygenic)[[3]]))
  temp.polygenic[, 1:2, ] <- polygenic
  polygenic <- temp.polygenic
  
  ## Calculate the Wald score
  polygenic[, 3, ] <- (polygenic[, 1, ] / polygenic[, 2, ]) ^ 2
  
  ## Calculate the Wald score's p-value
  ## Don't round p-values per Mariza's request, August 2015
  polygenic[, 4, ] <- pchisq.mix(polygenic[, 3, ],
                                 trait.count = trait.count,
                                 lod = F)

#  polygenic[, 4, ] <- round(pchisq.mix(polygenic[, 3, ],
#                                       trait.count = trait.count,
#                                       lod = F),
#                            4)
  
  ## Prepare the major.gene1 data
  major.gene1 <- effects$major.gene1
  start.index <- end.index + 1
  end.index <- start.index + random.effects.count - 1
  major.gene1.names <- variable.names[start.index:end.index]
  major.gene1 <- array(major.gene1, c(random.effects.count, 2, locus.count))
  dimnames(major.gene1) <- list(major.gene1.names,
                                c("Estimate", "Std.err"), loci.names)
  ## 69 and 73 are the integers chosen to represent a fixed externally and a
  ## fixed internally (respectively) standard error numerically
  major.gene1[, 2, ][major.gene1[, 2, ] == 69] <- NA
  major.gene1[, 2, ][major.gene1[, 2, ] == 73] <- NA

  ## Add Wald test, p-value, h2 and se.h2 to the major.gene1 compenents
  temp.major.gene1 <- array(NA, dim=c(dim(major.gene1)[1],
                                  dim(major.gene1)[2] + 5,
                                  dim(major.gene1)[3]),
                            dimnames=list(dimnames(major.gene1)[[1]],
                              c(dimnames(major.gene1)[[2]],
                                wald.data.names, herit.data.names),
                              dimnames(major.gene1)[[3]]))
  temp.major.gene1[, 1:2, ] <- major.gene1
  major.gene1 <- temp.major.gene1
  
  ## Calculate the Wald score
  major.gene1[, 3, ] <- (major.gene1[, 1, ] / major.gene1[, 2, ]) ^ 2
  
  ## Calculate the Wald score's p-value
  major.gene1[, 4, ] <- pchisq.mix(major.gene1[, 3, ],
                                   trait.count = trait.count,
                                   lod = F)
 
#  major.gene1[, 4, ] <- round(pchisq.mix(major.gene1[, 3, ],
#                                         trait.count = trait.count,
#                                         lod = F),
#                              4)
  
  ## Prepare the environmental data
  environmental <- effects$environmental
  start.index <- end.index + 1
  end.index <- start.index + random.effects.count - 1
  environmental.names <- variable.names[start.index:end.index]
  environmental <- array(environmental, c(random.effects.count, 2,
                                          locus.count))
  dimnames(environmental) <- list(environmental.names,
                                  c("Estimate", "Std.err"), loci.names)
  ## 69 and 73 are the integers chosen to represent a fixed externally and a
  ## fixed internally (respectively) standard error numerically
  environmental[, 2, ][environmental[, 2, ] == 69] <- NA
  environmental[, 2, ][environmental[, 2, ] == 73] <- NA
  
  ## Add Wald test, p-value, h2 and se.h2 to the environmental compenents
  temp.environmental <- array(NA, dim=c(dim(environmental)[1],
                                    dim(environmental)[2] + 2,
                                    dim(environmental)[3]),
                              dimnames=list(dimnames(environmental)[[1]],
                                c(dimnames(environmental)[[2]],
                                  wald.data.names),
                                  dimnames(environmental)[[3]]))
  temp.environmental[, 1:2, ] <- environmental
  environmental <- temp.environmental
  
  ## Calculate the Wald score
  environmental[, 3, ] <- (environmental[, 1, ] / environmental[, 2, ]) ^ 2
  
  ## Calculate the Wald score's p-value
  environmental[, 4, ] <- pchisq.mix(environmental[, 3, ],
                                     trait.count = trait.count,
                                     lod = F)
 
#  environmental[, 4, ] <- round(pchisq.mix(environmental[, 3, ],
#                                           trait.count = trait.count,
#                                           lod = F),
#                                4)
  
  ## Prepare the sibling.sibling data
  sibling.sibling <- effects$sibling.sibling
  start.index <- end.index + 1
  end.index <- start.index + random.effects.count - 1
  sibling.sibling.names <- variable.names[start.index:end.index]
  sibling.sibling <- array(sibling.sibling, c(random.effects.count, 2,
                                              locus.count))
  dimnames(sibling.sibling) <- list(sibling.sibling.names,
                                    c("Estimate", "Std.err"), loci.names)
  ## 69 and 73 are the integers chosen to represent a fixed externally and a
  ## fixed internally (respectively) standard error numerically
  sibling.sibling[, 2, ][sibling.sibling[, 2, ] == 69] <- NA
  sibling.sibling[, 2, ][sibling.sibling[, 2, ] == 73] <- NA
  
  ## Add Wald test, p-value, h2 and se.h2 to the sibling.sibling compenents
  temp.sibling.sibling <- array(NA, dim=c(dim(sibling.sibling)[1],
                                     dim(sibling.sibling)[2] + 2,
                                     dim(sibling.sibling)[3]),
                                dimnames=list(dimnames(sibling.sibling)[[1]],
                                  c(dimnames(sibling.sibling)[[2]],
                                    wald.data.names),
                                  dimnames(sibling.sibling)[[3]]))
  temp.sibling.sibling[, 1:2, ] <- sibling.sibling
  sibling.sibling <- temp.sibling.sibling
  
  ## Calculate the Wald score
  sibling.sibling[, 3, ] <- (sibling.sibling[, 1, ]
                             / sibling.sibling[, 2, ]) ^ 2
  
  ## Calculate the Wald score's p-valu
  sibling.sibling[, 4, ] <- pchisq.mix(sibling.sibling[, 3, ],
                                       trait.count = trait.count,
                                       lod = F)

#  sibling.sibling[, 4, ] <- round(pchisq.mix(sibling.sibling[, 3, ],
#                                             trait.count = trait.count,
#                                             lod = F),
#                                  4)  

  ## Prepare the parent.parent data
  parent.parent <- effects$parent.parent
  start.index <- end.index + 1
  end.index <- start.index + random.effects.count - 1
  parent.parent.names <- variable.names[start.index:end.index]
  parent.parent <- array(parent.parent, c(random.effects.count, 2,
                                          locus.count))
  dimnames(parent.parent) <- list(parent.parent.names,
                                  c("Estimate", "Std.err"), loci.names)
  ## 69 and 73 are the integers chosen to represent a fixed externally and a
  ## fixed internally (respectively) standard error numerically
  parent.parent[, 2, ][parent.parent[, 2, ] == 69] <- NA
  parent.parent[, 2, ][parent.parent[, 2, ] == 73] <- NA
  
  ## Add Wald test, p-value, h2 and se.h2 to the parent.parent compenents
  temp.parent.parent <- array(NA, dim=c(dim(parent.parent)[1],
                                   dim(parent.parent)[2] + 2,
                                   dim(parent.parent)[3]),
                              dimnames=list(dimnames(parent.parent)[[1]],
                                c(dimnames(parent.parent)[[2]],
                                  wald.data.names),
                                dimnames(parent.parent)[[3]]))
  temp.parent.parent[, 1:2, ] <- parent.parent
  parent.parent <- temp.parent.parent
  
  ## Calculate the Wald score
  parent.parent[, 3, ] <- (parent.parent[, 1, ] / parent.parent[, 2, ]) ^ 2
  
  ## Calculate the Wald score's p-value
  parent.parent[, 4, ] <- pchisq.mix(parent.parent[, 3, ],
                                     trait.count = trait.count,
                                     lod = F)

#  parent.parent[, 4, ] <- round(pchisq.mix(parent.parent[, 3, ],
#                                           trait.count = trait.count,
#                                           lod = F),
#                                4)  
  
  ## Prepare the parent.offspring data
  parent.offspring <- effects$parent.offspring
  start.index <- end.index + 1
  end.index <- start.index + random.effects.count - 1
  parent.offspring.names <- variable.names[start.index:end.index]
  parent.offspring <- array(parent.offspring, c(random.effects.count, 2,
                                                locus.count))
  dimnames(parent.offspring) <- list(parent.offspring.names,
                                     c("Estimate", "Std.err"), loci.names)
  ## 69 and 73 are the integers chosen to represent a fixed externally and a
  ## fixed internally (respectively) standard error numerically
  parent.offspring[, 2, ][parent.offspring[, 2, ] == 69] <- NA
  parent.offspring[, 2, ][parent.offspring[, 2, ] == 73] <- NA
  
  ## Add Wald test, p-value, h2 and se.h2 to the parent.offspring compenents
  temp.parent.offspring <- array(NA, dim=c(dim(parent.offspring)[1],
                                      dim(parent.offspring)[2] + 2,
                                      dim(parent.offspring)[3]),
                                 dimnames=list(dimnames(parent.offspring)[[1]],
                                   c(dimnames(parent.offspring)[[2]],
                                     wald.data.names),
                                   dimnames(parent.offspring)[[3]]))
  temp.parent.offspring[, 1:2, ] <- parent.offspring
  parent.offspring <- temp.parent.offspring
  
  ## Calculate the Wald score
  parent.offspring[, 3, ] <- (parent.offspring[, 1, ]
                              / parent.offspring[, 2, ]) ^ 2
  
  ## Calculate the Wald score's p-value
  parent.offspring[, 4, ] <- pchisq.mix(parent.offspring[, 3, ],
                                              trait.count = trait.count,
                                              lod = F)
                             
#  parent.offspring[, 4, ] <- round(pchisq.mix(parent.offspring[, 3, ],
#                                              trait.count = trait.count,
#                                              lod = F),
#                                   4)  

  if(run.alternative.hyps) {
    distance <- unlist(sapply(loci.names, get.dist))
  }else {
    distance <- NA
  }
  log.likelihood <- effects$log.likelihoods
  log.lik.status <- effects$log.lik.status
  ## I received different values for lod.score and summing the values
  ## directly from the fam.log.liks.  The is due to the loss of precision
  ## in the log.liklihood values calculated and printed to summary.log.  The
  ## lod scores in multic are calculated from these summary.log values.  With
  ## all that said, it is only a loss of precision at about the
  ## one-tenthousandth place.  Eric Lunde 2005-06-06
  lod.score <- -2 * (log.likelihood[1] - log.likelihood) / 4.6
  lod.score[1] <- NA

  ## Mariza no longer wants the p-values rounded. Pat Votruba 2015-08-05
  #p.value <- round(pchisq.mix(lod.score, trait.count = trait.count), 4)
  p.value <- pchisq.mix(lod.score, trait.count = trait.count)
  
  log.liks <- data.frame(log.likelihood, distance, log.lik.status,
                         lod.score = round(lod.score, digits = 4), p.value,
                         row.names=loci.names)

  ## Order the log likelihoods based on centimorgan value
  log.liks <- rbind(log.liks[1, ],
                    log.liks[-1, ][order(as.numeric(log.liks[-1, 2])), ])

  effects.list <- list(fixed.effects = fixed.effects,
                       polygenic = polygenic,
                       major.gene1 = major.gene1,
                       environmental = environmental,
                       sibling.sibling = sibling.sibling,
                       parent.parent = parent.parent,
                       parent.offspring = parent.offspring,
                       log.liks = log.liks,
                       polygenic.names = polygenic.names,
                       major.gene1.names = major.gene1.names,
                       environmental.names = environmental.names,
                       sibling.sibling.names = sibling.sibling.names,
                       parent.parent.names = parent.parent.names,
                       parent.offspring.names = parent.offspring.names)
  
  return (effects.list)
}

get.dist <- function(ibd.name) {
  tokens <- multic.strsplit(ibd.name, sep = ".")
  tokens.length <- length(tokens)
  if(tokens[1] == "ibd" && tokens.length > 2) {
    dist <- as.numeric(paste(tokens[4:tokens.length], collapse = "."))   
  }else if(tokens[1] == "mibd") {
    dist <- as.numeric(paste(tokens[3:tokens.length], collapse = "."))   
  }else {
    dist <- NA
  }
  
  ##  if(tokens.length > 2) {
  ##    dist <- as.numeric(paste(tokens[-(1:(length(tokens)-2))], collapse = "."))
  ##  } else {
  ##    dist <- NA
  ##  }
  
  return (dist)
}
