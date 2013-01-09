plot.multic <- function(x, ylim=NULL, xlim=NULL, lty=1, col=1, lwd=1, ...) {
  if(!is.multic(x)) {
    stop("\n", substitute(x), "is not a multic object\n")
  }
  multic.object <- x
  
  cM <- as.numeric(multic.object$log.liks$distance[-1])
  lod.scores <- as.numeric(multic.object$log.liks$lod.score[-1])

  if(length(cM) == 0) {
    stop(paste("\nThere are no lod scores to plot because",
               "there were no alternative\nhypotheses calculated. ",
               "You must specify an mloci.out file to the multic",
               "\ncall to calculate alternative hypotheses.\n"))
  } else if(length(lod.scores) == 0) {
    stop(paste("\nThere are no lod scores to plot because",
               "there were no alternative\nhypotheses calculated. ",
               "You must specify an mloci.out file to the multic",
               "\ncall to calculate alternative hypotheses.\n"))
  } else if(sum(is.na(cM)) > 0) {
    stop(paste("\nThere are some lod scores that have no distance ",
               "measurement ",
               "(possibly due\nto ibd names without a distance component).  ",
               substitute(x), " cannot be\nplotted with missing distance ",
               "data.  Look at ", substitute(x), "$log.liks to\ndetermine ",
               "which ibds are causing this error.\nplot.multic.q key 25",
               sep = ""))
  }
  
  max.lod.score <- max(lod.scores)
  min.lod.score <- min(lod.scores)

  if(is.null(xlim))  xlim <- c(cM[1], cM[length(cM)])

  y.min <- 0
  y.min <- ifelse(min.lod.score >= 0, 0, min.lod.score - 1)

  if(is.null(ylim))  ylim <- c(y.min, max.lod.score + 1)

  ## Draw the actual curve of lod scores
  plot(x = cM,
       y = lod.scores,
       ...,
       type = 'l',
       lty = lty,
       lwd = lwd,
       col=col,
       xlim = xlim,
       ylim = ylim,
       xlab = "cM",
       ylab = "LOD")
  
  if(y.min < 0) {
    abline(h = 0)
  }  
  
  invisible()
}
