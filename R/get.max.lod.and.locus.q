get.max.lod.and.locus <- function(multic.obj) {  
  ## Gain direct acces to the log likelihoods, log odd scores, and locus names
  log.liks <- multic.obj$log.liks
  lods <- log.liks[, 4]
  locus.names <- dimnames(log.liks)[[1]]

  if(locus.names[1] == "null") {
    lods[1] <- -Inf
  }

  ## Find the maximum lod score
  max.lod.score <- max(lods, na.rm=TRUE)

  locus.distances <- log.liks[, 2]

  ## Finally, get the name of the locus with the maximum lod score and report
  ## it to the user.
  max.lod.locus <- locus.names[lods == max.lod.score]

  ## Get the centimorgan distance (if availiable)
  max.lod.centimorgan <- locus.distances[lods == max.lod.score]

  return (list(max.lod.score=max.lod.score,
               max.lod.locus=max.lod.locus,
               max.lod.centimorgan=max.lod.centimorgan))
}
