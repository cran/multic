## Calculate the centiMorgan around the location that produced the max lod
## score that also have lod scores greater than or equal to the max lod
## score - 1
get.cm.close.to.peak <- function(lod.scores, centimorgan, epsilon = 1) {
  ## If there is no distance (centimorgan) data, this calculation is worthless
  ## This is a test to see if there is any distance data

  ## Force that centimorgan and lod.scores are numeric.  Previous multic
  ## object versions stored these as character, thus throwing off the
  ## sorting.
  centimorgan <- as.numeric(centimorgan)
  lod.scores <- as.numeric(lod.scores)

  all.na.centimorgan <- all(apply(array(centimorgan), 1, is.na))
  if(all.na.centimorgan) {
    warning(paste("There is no centimorgan data to determine a range around",
                  "the max lod score."))
    return (NULL)
  }

  ## Get the ascending order of centimorgan values
  ascending.order <- order(centimorgan)

  ## Sort the centimorgan (and the lod scores they produce) in increasing order
  centimorgan <- centimorgan[ascending.order]
  lod.scores <- lod.scores[ascending.order]

  ## Get the maximum lod odds score and the centimorgan the index of the
  ## centimorgan that generates it.
  max.lod.score <- max(lod.scores)
  index.of.max.lod.score <- match(max.lod.score, lod.scores)

  ## centimorgan.of.interest is all the centimorgan that produced a lod score
  ## greater than or equal to the max lod score - epsilon
  centimorgan.of.interest <- centimorgan[lod.scores >=
                                         (max.lod.score - epsilon)]

  ## clusters.above.threshold is the groupings of indices of centimorgan that
  ## produced lod scores above the max lod score - 1
  ## The indices of clusters.above.threshold that hold the value NA correspond
  ## to the indices of centimorgan values that did not produce at least a lod
  ## score of the max lod score - epsilon 
  clusters.above.threshold <- match(centimorgan, centimorgan.of.interest)

  ## first.NA.greater is the first index of an NA (centimorgan that is larger
  ## than the location that produced the max that did not produce a lod score
  ## bigger than the max lod score - epsilon
  first.NA.greater <-
    match(NA, clusters.above.threshold[index.of.max.lod.score +
                                       (1 : (length(centimorgan) -
                                             index.of.max.lod.score)) ])

  biggest.index.of.interest <- ifelse(is.na(first.NA.greater),
                                      length(centimorgan),
                                      index.of.max.lod.score
                                      + first.NA.greater - 1)

  first.NA.less <-
    match(NA, clusters.above.threshold[(index.of.max.lod.score - 1) : 1])

  lowest.index.of.interest <- ifelse(is.na(first.NA.less),
                                     1,
                                     index.of.max.lod.score
                                     - first.NA.less + 1)

  ## return the centimorgan that were adjacent to the peak and also produced
  ## a lod score of above (or equal to) the max lod score - epsilon
  return (centimorgan[lowest.index.of.interest : biggest.index.of.interest])
}

summary.multic <- function(object, ...) {
  if(!is.multic(object)) {
    stop("\n", substitute(object), "is not a multic object\n")
  }
  multic.object <- object
  
  ## n is a variable used to determine the top n families.  I put this
  ## assignment up here so we can change it easier later.
  n <- 5

  call <- multic.object$call

  if(is.null(multic.object$metadata$mloci.out)) {
    summary <- list(call = call)
  }else {
    max.lod.and.locus <- get.max.lod.and.locus(multic.object)
    
    max.lod.score <- max.lod.and.locus$max.lod.score
    max.lod.locus <- max.lod.and.locus$max.lod.locus
    max.lod.centimorgan <- max.lod.and.locus$max.lod.centimorgan
    
    top.n.families <- NULL
    if(multic.object$metadata$calc.fam.log.liks) {
      top.n.families <- get.top.n.families(multic.object, n)
    }
      
    ## Get the centimorgan that produced values above the maximum lod score - 1
    ## that are also contiguous to the centimorgan that produced the maximum
    lod.scores <- multic.object$log.liks$lod.score[-1]
    centimorgan <- multic.object$log.liks$distance[-1]
    
    centimorgan.close.to.peak <- get.cm.close.to.peak(lod.scores, centimorgan)
    
    ## Build summary.multic.object
    summary <- list(call = call,
                    max.lod.score = max.lod.score,
                    max.lod.locus = max.lod.locus,
                    max.lod.centimorgan = max.lod.centimorgan,
                    n = n,
                    top.n.families = top.n.families,
                    centimorgan.close.to.peak = centimorgan.close.to.peak
                    )
  }
  
  oldClass(summary) <- 'summary.multic'
  return (summary)
}
