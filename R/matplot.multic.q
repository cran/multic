
matplot.multic <- function(multic.objects, xlab='cM', ylab='LOD', xlim=NULL, ylim=NULL, lty=1, col=1, main='', ...) {

  n.obj <- length(multic.objects)
  results <- NULL

  for(mj in 1:n.obj) {

    multic.object <- multic.objects[[mj]]
    ### code from plot.multic ###
    obj.lod  <-  cbind(as.numeric(multic.object$log.liks$distance),
                       as.numeric(multic.object$log.liks$lod.score),
                       rep(mj,length(multic.object$log.liks$distance)))

    results <- rbind(obj.lod, results)
  }

  results <- results[!is.na(results[,1]),]
  results <- results[order(results[,3],results[,1]),]
  
  cM  <-  results[,1]
  lod.scores <- results[,2]
  group <- results[,3]

  ## change all negative lod values to 0
  lod.scores <- ifelse(lod.scores<0,0,lod.scores)
  
  if(is.null(xlim)) {
    xlim  <-  range(cM, na.rm=T)
  }
  
  if(is.null(ylim)) {
    max.lod.score  <-  max(lod.scores, na.rm=T)
    ylim  <-  c(0, max.lod.score + 1)
  }
  if(length(col)==1) {
    col  <-  1:n.obj
  }
  if(length(lty)==1) {
    lty  <-  1:n.obj
  }
  
  ## Draw the actual curve of lod scores  
  plot(x= cM, y= lod.scores, type='n', xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)
  for(i in 1:n.obj) {
    ok <- group==i
    lines(cM[ok], lod.scores[ok], col=col[i], lty=lty[i], ...)
  }
  

  invisible(results)
} 

