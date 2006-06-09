#############################################################################
## get.top.n.families
## Report the top n family identifiers at the max lod score.  If there is
## more than one max lod score, the first of the cM is chosen.  There are no
## plans to be able to choose which max lod score cM to choose.
## Eric Lunde
#############################################################################
get.top.n.families <- function(multic.obj, n) {
  if(n < 1) {
    stop(paste("n must be positive, n = ", n, "\n\n", sep=""))
  }

  max.lod.locus <- get.max.lod.and.locus(multic.obj)$max.lod.locus
  ## Ideally, max.lod.locus is a single value.  When it is not, we need to
  ## pare it down to only the first max value.  Choosing the first max value
  ## is completely arbitrary.  We will still report all the cM that
  ## have that LOD, but the == below needs only one file.
  max.lod.locus <- max.lod.locus[1]
  
  ## Get all of the family data from the specified locus
  all.fam.at.max.locus <-
    multic.obj$fam.log.liks[, , dimnames(multic.obj$fam.log.liks)[[3]]
                            == max.lod.locus]

  ## Get the lod odd scores. Since the order function sorts in increasing
  ## order, I negate all values because I want the n largest scores.  These
  ## values will then be the at the beginning of the vector after sorting
  lods <- -all.fam.at.max.locus[, 2]
  
  ## Get the n best family log odd scores
  n.best.families <- all.fam.at.max.locus[order(lods), ][1:n, ]

  return(n.best.families)
}
