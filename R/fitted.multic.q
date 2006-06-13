fitted.multic <- function(object, collapse.family = TRUE, position = 0, ...)
{
  ## STILL NEED TO ADD IN POSITION
  if(is.character(object$residuals)) {
    stop("\n", substitute(object), " does not contain residual data.\n",
         "Please re-run multic with \"calc.residuals = TRUE\" to calculate ",
         "them.\n",
         "fitted.multic.q key 7\n")
  }

  traits <- object$metadata$traits
  resid.raw <- object$residuals
  trait.names <- object$metadata$trait.names
  n.traits <- length(trait.names)
  n.subj <- nrow(traits)
  n.fam <- nrow(resid.raw)
  
  ## manipulate residuals so that they are in the right format
  complete.resid <- matrix(NA, ncol=n.traits, nrow=n.subj,
                           dimnames=dimnames(traits))
  for(i in 1:n.traits) {
    complete.resid[!is.na(traits[, 1]),i] <-
      # This [, 1] gets to the null hyp, it will have to be more generic.
      unlist(lapply(resid.raw[, position + 1], function(x,i,n.traits) 
                    matrix(x, ncol=n.traits)[,i],i,n.traits))
  }
  
  fitted <- traits - complete.resid  ## y-yhat=resid so y-resid=yhat

  ## now collapse by family
  if(collapse.family) {
    fam.id <- multic.strings.split(dimnames(complete.resid)[[1]],'-')[, 1]
    fitted2 <- matrix(NA, ncol=n.traits, nrow=n.fam, dimnames=
                      list(dimnames(resid.raw)[[1]], trait.names))

    for(i in 1:n.traits){
      fitted2[,i] <-
        sapply(split(fitted[,i], fam.id), mean, na.rm=T) [unique(fam.id)]
    }
    fitted <- fitted2
  }
  
  return(fitted)
  
}
