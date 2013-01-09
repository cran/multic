
########################################################################
# Function: addGE
# Purpose:  Assess whether running "multivariate" multic is worthwhile
#           by looking at various combinations of the univariate results.
#
# Returns:  p-value, -log10pvalue, position, and chromosome
#           Only done where all variances/standard Errors are NOT NA
#
# Author: EJA (using functions created by Mariza de Andrade, Curt Olswold)
# Created: 7/10/2005
#
# Parameters:
#   multic.objs: List of multic objects for combining gene effects
#   combine:     Indicate if how many of the multic.objs should be combined
#   plotit:      Plot position vs LOD values
#   ibd.dist:    Distance information
#                Default is to use location from multic object
#   statistic:   Use LRT or Wald test statistic (default=LRT)
#   legend:      Automatically include legend (default=T)
#
# test/example:
# x.2 <- addGE(list(fit.x1, fit.x2, fit.x3), combine=2, plotit=T)
#
########################################################################

addGE <- function(multic.objs, combine=2, plotit=FALSE, ibd.dist, statistic=c("lrt","wald"),
                  legend=TRUE, ylim=NULL, ...) {


  statistic <- match.arg(statistic)
  
  if(missing(ibd.dist)) cM <- as.numeric(multic.objs[[1]]$log.liks$distance[-1]) else
    cM <- ibd.dist
    
  ## check to see that all objects have the same length
  ## if using mibd files
  
  n.traits <- length(multic.objs)
  n.markers <- length(cM)
  n.comb <- choose(n.traits, combine)
  check.lengths <- any(unlist(lapply(multic.objs,
                       function(x) length(as.numeric(x$log.liks$distance[-1])))) != n.markers)

  if(check.lengths) stop('You have multic objects using different numbers of markers\n or your ibd distance file is the wrong length')



  ###################################################################################
  ###################################################################################

  if(statistic=="lrt") {

     constant <- 2/log10(exp(1))

     ## extract all the LOD estimates (except the first one which is NA)
     LOD <- matrix(unlist(lapply(multic.objs,
                                 function(x) x$log.liks$lod.score[-1])), ncol=n.traits) 
     LRT <- LOD*constant
  
     ## remove rows with NAs
     rowNA <- as.logical(rowSums(is.na(LRT)))
     LRT <- LRT[!rowNA,]
     cM <- cM[!rowNA]
     n.markers <- length(cM)
  
     if(nrow(LRT) ==1) stop('You only have 1 marker with non-missing results')

     ## COMBINE BY 2, BY 3, ETC.
     pullfunc.lrt <- function(x, n.traits, combine) {
       gene.eff.lrt <- function(lrt, combine = 2) {
         ##  The number of markers/traits used
         n.traits <- length(lrt)
         ## Possible number of combinations
         n.comb <- choose(n.traits, combine)
         all.comb <- matrix(subsets(n.traits, combine), ncol = combine)

         ## Calculate the LRT statistic 
                                        #  pChi  p-value, chisquare distribution
         tmp <- matrix(lrt[all.comb], nrow=n.comb)
         statChi <- rowSums(tmp)
         pChi <- pchisq.mix(statChi, 1, lod = F)
         unique.names <- apply(all.comb, 1, paste, collapse = ", ")
         results <- data.frame(comb = unique.names, statChi = statChi,
                               pChi = pChi)
         return(results)
       }

       tmp <- gene.eff.lrt(x[1:n.traits], combine=combine)
       return(tmp)
     }

     ## a list of matrices
     GEresults <- apply(LRT, 1, pullfunc.lrt, n.traits=n.traits, combine=combine) 

     ## create overall results each as vectors
     comb <- unlist(lapply(GEresults, function(x) x[,1]))
     statChi <- unlist(lapply(GEresults, function(x) x[,2]))
     pChi <- unlist(lapply(GEresults, function(x) x[,3]))
     cMrep <- rep(cM, rep(n.comb, n.markers))
     lod <- -log10(pChi)

     ge.all <- data.frame(comb, statChi, pChi, cM=cMrep, lod, method=statistic)
   }

  ###################################################################################
  ###################################################################################

  if(statistic=="wald"){
    ## extract all the Major Gene variance estimates
    MG <- matrix(unlist(lapply(multic.objs,
                               function(x) x$major.gene1[1,1,] [-1])), ncol=n.traits)
  
    ## extract all the SE estimates for the MG
    SE <- matrix(unlist(lapply(multic.objs,
                               function(x) x$major.gene1[1,2,] [-1])), ncol=n.traits)
  
    ## remove rows with NAs
    rowNA <- as.logical(rowSums(is.na(MG)) + rowSums(is.na(SE)))
    MG <- MG[!rowNA,]
    SE <- SE[!rowNA,]
    cM <- cM[!rowNA]
    n.markers <- length(cM)
  
    if(nrow(MG) ==1) stop('You only have 1 marker with non-missing results')

    ## COMBINE BY 2, BY 3, ETC.
    both <- cbind(MG,SE)

    pullfunc <- function(x, n.traits, combine) {
      tmp <- gene.eff(x[1:n.traits], x[(n.traits+1):(2*n.traits)], combine=combine)
      return(tmp)
    }

    ## a list (length n.markers) of matrices (columns=comb, statChi, pChi) 
    GEresults <- apply(both, 1, pullfunc, n.traits=n.traits, combine=combine) 

    ## create overall results each as vectors
    comb <- unlist(lapply(GEresults, function(x) x[,1]))
    statChi <- unlist(lapply(GEresults, function(x) x[,2]))
    pChi <- unlist(lapply(GEresults, function(x) x[,3]))
    cMrep <- rep(cM, rep(n.comb, n.markers))
    lod <- -log10(pChi)

    ge.all <- data.frame(comb, statChi, pChi, cM=cMrep, lod, method=statistic)
  }

  ###################################################################################
  ###################################################################################

  if(plotit==TRUE) {
    if(!is.null(ylim)) my.ylim <- ylim
    if(is.null(ylim)) my.ylim <- c(0, max(c(3,lod),na.rm=T) + .5)
    
    ## Plot the -log10 p-values
    plot(x=cMrep, y=lod, type='n', ylab=paste('-log10(',casefold(statistic,upper=T),'P-Value)'),
         xlab=' Position', ylim=my.ylim, ...)
    for(i in unique(comb)) {
      ok <- comb == i
      code <- match(i,unique(comb))
      lines(x=cMrep[ok], y=lod[ok], lty=code, lwd=2, col=code)
    }

    if(legend==TRUE){

      legend(x=min(cMrep), y=par()$usr[4],
             lty=1:n.comb, col=1:n.comb, legend=unique(comb),
             ncol=n.comb,lwd=2)
    }

  }
  ###################################################################################
  ###################################################################################
  
  invisible(ge.all)

} ### end multicGeneEffects


