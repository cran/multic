residuals.multic <- function(object,
                             type = c("working", "standardized", "Q1", "Q2"),
                             collapse.family=TRUE, position=0, ...)
{
  ## STILL NEED TO ADD IN POSITION
  ## could also add an option, collapse residuals using 'mean' or other
  ## function

  if( length(find("gchol")) == 0 ) {
    stop("\nThe kinship2 function \"gchol\" is not available.\n",
         "You must have the kinship2 package installed.",
         "\n(Splus 7 includes it.\n")
  }
  
  if(is.character(object$residuals)) {
    stop("\n", substitute(object), " does not contain residual data.\n",
         "Please re-run multic with \"calc.residuals = TRUE\" to calculate ",
         "them.\n",
         "residuals.multic.q key 13\n")
  }

  type <- match.arg(type)
  
  ## extract basic information.  Need residuals only for one location.
  resid.raw <- object$residuals
  v.raw <- object$v.matrices
  traits <- object$metadata$traits
  id.names <- dimnames(traits)[[1]]
  trait.names <- object$metadata$trait.names
  loci.names <- row.names(object$metadata$iterations)
  n.traits <- length(trait.names)
  n.subj <- nrow(traits)
  n.fam <- nrow(resid.raw)
  
  ## first do the calculations by individual.  
  ## Note that it should be the same length as the original data.
  if( !collapse.family ) {
    r <- switch(type, working=
                {
                  complete.resid <- matrix(NA, ncol=n.traits, nrow=n.subj,
                                           dimnames=dimnames(traits))
                  for(i in 1:n.traits) {
                    complete.resid[!is.na(traits[, 1]),i] <- unlist(lapply(resid.raw,
                                                                           function(x,i,n.traits) matrix(x, ncol=n.traits)[,i],i,n.traits))}
                  complete.resid
                },
                standardized=
                {
                  ## Standardized Residuals
                  stdfn <- function(j,x,y,n.traits) {
                    n.ped <- nrow(x[[j]])/n.traits
                    if(n.ped==0) return(NULL)
                    else {
                      results <- matrix(NA, ncol=n.traits, nrow=nrow(x[[j]])/n.traits) 
                      for(i in 1:n.traits) {
                        beg <- (i-1)*n.ped + 1
                        end <- (i-1)*n.ped + n.ped
                        res <- y[[j]][beg:end]
                        vmat <- x[[j]][beg:end,beg:end, drop = FALSE]
                        gcholmat <- kinship2::gchol(vmat)
                        stdres <- solve(gcholmat,res, full=FALSE)
                        results[,i] <- stdres 
                      }
                      return(results)
                    }
                  }
                  
                  std.resid.raw <- lapply(1:nrow(resid.raw), stdfn,
                                          x=v.raw, y=resid.raw, n.traits=n.traits)

                  std.resid.raw <- std.resid.raw[sapply(std.resid.raw, function(x)!is.null(x))]

                  ## return a matrix that has the same length as the original data
                  complete.std.resid <- matrix(NA, ncol=n.traits, nrow=n.subj,
                                               dimnames=dimnames(traits))
                  for(i in 1:n.traits) {
                    ## STILL NEED TO DEAL with families where everyone is missing!!
                    complete.std.resid[!is.na(traits[, 1]),i] <- unlist(lapply(std.resid.raw,
                                                                               function(x,i,n.traits) matrix(x, ncol=n.traits)[,i],i,n.traits))}
                  complete.std.resid
                },
                Q1=stop('Q1 is only defined by family'),
                Q2=stop('Q2 is only defined by family')
                )
    return(r)
  }
  
  ## Now do the calculations by family (averaging within a family)
  if(collapse.family) {
    r <- switch(type, working=
                matrix(unlist(lapply(resid.raw,function(x,n.traits)
                                     colMeans(matrix(x, ncol=n.traits),na.rm=TRUE),n.traits)),
                       ncol=n.traits,byrow=TRUE),
                standardized=
                {
                  ## Standardized Residuals
                  stdfn <- function(j,x,y,n.traits) {
                    n.ped <- nrow(x[[j]])/n.traits
                    if(n.ped==0) return(matrix(NA,ncol=n.traits,nrow=1))
                    else {
                      results <- matrix(NA, ncol=n.traits, nrow=nrow(x[[j]])/n.traits) 
                      for(i in 1:n.traits) {
                        beg <- (i-1)*n.ped + 1
                        end <- (i-1)*n.ped + n.ped
                        res <- y[[j]][beg:end]
                        vmat <- x[[j]][beg:end,beg:end, drop = FALSE]
                        gcholmat <- kinship2::gchol(vmat)
                        stdres <- solve(gcholmat,res, full = FALSE)
                        results[,i] <- stdres 
                      }
                      return(results)
                    }
                  }
                  
                  std.resid.raw <- lapply(1:nrow(resid.raw), stdfn,
                                          x=v.raw, y=resid.raw, n.traits=n.traits)

                  ## create a mean value for each trait and each family
                  matrix(unlist(lapply(std.resid.raw,colMeans,na.rm=TRUE)),
                         ncol=n.traits,byrow=TRUE)
                },
                Q1={
                  ## Q1 function
                  q1fn <- function(j,R,V,n.traits) {
                    n.ped <- nrow(V[[j]])/n.traits

                    if(n.ped==0) return(matrix(NA,ncol=n.traits,nrow=1))
                    else {
                      results <- matrix(NA, ncol=n.traits, nrow=1) 
                      for(i in 1:n.traits) {
                        beg <- (i-1)*n.ped + 1
                        end <- (i-1)*n.ped + n.ped
                        res <- R[[j]][beg:end]
                        vmat <- V[[j]][beg:end,beg:end, drop = FALSE]
                        gcholmat <- kinship2::gchol(vmat)
                        Q <- t(res)%*%solve(gcholmat)%*%res
                        ans <- sqrt(2*Q) - sqrt(2*n.ped - 1)
                        results[,i] <- ans
                      }
                      return(results)
                    }
                  }
                  
                  Q1 <- lapply(1:nrow(resid.raw), q1fn,
                               V=v.raw, R=resid.raw, n.traits=n.traits)

                  ## create a mean value for each trait and each family
                  matrix(unlist(Q1), ncol=n.traits,byrow=TRUE)

                },
                Q2={
                  ## Q2 function
                  q2fn <- function(j,R,V,n.traits) {
                    n.ped <- nrow(V[[j]])/n.traits

                    if(n.ped==0) return(matrix(NA,ncol=n.traits,nrow=1))
                    else {
                      results <- matrix(NA, ncol=n.traits, nrow=1) 
                      for(i in 1:n.traits) {
                        beg <- (i-1)*n.ped + 1
                        end <- (i-1)*n.ped + n.ped
                        res <- R[[j]][beg:end]
                        vmat <- V[[j]][beg:end,beg:end, drop = FALSE]
                        gcholmat <- kinship2::gchol(vmat)
                        Q <- t(res)%*%solve(gcholmat)%*%res
                        ans <- sqrt(2*Q) - sqrt(2*n.ped - 1)
                        ans <- ((Q/n.ped)^(1/3) - 1 + 2/(9*n.ped)) * sqrt(9*(n.ped/2))
                        results[,i] <- ans
                      }
                      return(results)
                    }
                  }
                  
                  Q2 <- lapply(1:nrow(resid.raw), q2fn,
                               V=v.raw, R=resid.raw, n.traits=n.traits)

                  ## create a mean value for each trait and each family
                  matrix(unlist(Q2), ncol=n.traits,byrow=TRUE)

                }
                )
    return(r) 
  }
}
