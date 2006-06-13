# function: gene.eff.s
#
# The goal of this function is to estimate the overall contribution of two or
# more genes to a specific trait or the contribution of a specific gene/marker
# to two or more traits.
#
#  Required:
#  MG = vector of variance estimates.
#  SE = vector of the standard errors of the elements of MG
#
#  Optional:
#  combine = number of traits to combine (look at all possible)
#  normal = F   outputs only the chisquare results
#         = T   outputs normal and chisquare results
#
#  Example:
#
#  mg <- c(10.89, 1.82, 13.22, 9.93)
#  SE <- c(8.94, 11.86, 7.05, 7.67)
#  n.traits = 4 (number of traits) 
#  combine = 2
#
#  gene.eff(mg, SE, combine=2)

gene.eff <- function(mg, se, combine=2, normal=F){

  if (length(mg) != length(se))
    stop("Check the length of mg and se must be the same")

  if (any(mg < 0)) stop("Cannot have negative variances")
  if (any(se < 0)) stop("Cannot have negative standard errors")
  
  #  The number of markers/traits used
  n.traits <- length(mg)

  # Possible number of combinations
  n.comb <- choose(n.traits, combine)
  all.comb <- matrix(subsets(n.traits, combine), ncol=combine)
 
  # Calculate the statistic T
  #  statNormal  normal distribution
  #  statChi chisquare distribution
  #  pNormal  p-value, normal distribution (half-normal p-value = 2*p.normal)
  #  pChi  p-value, chisquare distribution
  
  tmp <- matrix(mg[all.comb],ncol=combine)/matrix(se[all.comb],ncol=combine)

  statNormal <- rowSums(tmp)/sqrt(combine)
  pNormal  <- 1 - pnorm(statNormal)
  statChi <- (rowSums(tmp)/sqrt(combine))^2  
  pChi  <- pchisq.mix(statChi, 1, lod=F)

  unique.names <- apply(all.comb, 1, paste, collapse=', ')

  if (normal == T){
    results <- data.frame(comb=unique.names, statNormal=statNormal, pNormal=pNormal,
                statChi=statChi, pChi=pChi)
  }
  else { 
    results <- data.frame(comb=unique.names, statChi=statChi, pChi=pChi)
  }
  
  return(results)

}


