###############################################################################
## Title: load.var.sandwich
## Description: load.var.sandwich is responsible for allocating memory, filling
##              that memory with values read from the multic output file
##              'varCovar.log', and adjusting the dimensions to accurately
##              model the data.  varCovar.log holds the values of the variance
##              covariance matrix used in multic.
## Input: random.effects.count - random.effects.count is a scalar value
##                               representing the number of random effects per
##                               effect (polygene and major gene) that the
##                               last run of multic output.
##        alt.hyp.count - alt.hyp.count is the number of alternative hypotheses
##                        that the last run of mutlic calculated.  The null
##                        hypothesis is not included in this value.
##        random.effects.names - a character vector containing the names
##                               (dimnames or labels) for the first two
##                               dimensions of var.sandwich.
##        ibd.names - a character vector containing the names (dimnames or
##                    labels) of the names of the ibd files that will be the
##                    names of the third dimension of var.sandwich.
## Ouput: var.sandwich - var.sandwich is a three dimensional array where each index
##                    of the third dimension contains the two dimensional
##                    variance covariance information from the null
##                    hypothesis each of the alternative hypotheses.
## Side Effects: Since the file name 'varCovar.log' is hardcoded into the C++
##               routine that extracts the variance covariance matrix data, if
##               that file cannot be opened, the C++ routine will cause the
##               program (and the Splus session) to terminate.
##
##               As of 02-07-2004, if the file 'invExpSecDerRandom.log' doesn't
##               exist, a NULL object will be returned, the program will not
##               terminate.
## Author: Eric Lunde 10-27-03
## 11-28-03, Eric Lunde, 11-28-03, I removed the parameters trait.count and
##           replaced it with random.effects.count.  I made a decision that
##           this effect count value was the important (useful) value.  So it is
##           calculated early in the program and the value is passed through,
##           rather than each function calculating the value again and again.
## 02-07-2004, Eric Lunde, I added two more input parameters:
##             fixed.effects.names and ibd.names.  They are the dimnames for
##             the return value var.sandwich.  Also I added a file
##             existing test at the beginning to see if the file
##             'varCovar.log' exists or not.
## 03-31-2005, Eric Lunde, The number of environmental variance (covariance)
##             values no longer is dependent upon whether multic is run in
##             multivariate or longitudinal mode (it now always operates as if
##             in multivariate mode).  Therefore I've removed the
##             environment.count variable and replaced it with the
##             random.effects.count which now always is the correct value.
###############################################################################
load.var.sandwich <- function(random.effects.count, alt.hyp.count,
                              random.effects.names, ibd.names,
                              constraints) {
  if( !file.exists("varCovar.log") ) {
    cat(paste("- 'varCovar.log' was not found\n",
              "- No robust variance covariance information will be loaded.\n",
              sep=""))
    return ("No robust variance covariance data was calculated for this multic object.")
  }
  
  
  ## Make sure the parameters are in their correct ranges
  random.effects.count <- as.integer(random.effects.count)
  if(random.effects.count < 0) {
    stop(paste("random.effects.count must be non-negative.  Current value is",
               random.effects.count, "\nload.var.sandwich.s key 35"))
  } 

  alt.hyp.count <- as.integer(alt.hyp.count)
  if(alt.hyp.count < 0) {
    stop(paste("alt.hyp.count must be non-negative.  Current value is",
               alt.hyp.count, "\nload.var.sandwich.s key 41"))
  }

  ## n.vc is the number of variance components (including sib, par, and
  ## par-off)
  n.vc <- sum(constraints[-1] != "F")
  is.poly.fixed <- as.integer(constraints[2] == "F")
  is.mg1.fixed <- as.integer(constraints[3] == "F")

  var.sandwich <- numeric( (n.vc * random.effects.count)
                          * (n.vc * random.effects.count)
                          * (alt.hyp.count + 1) )
  
  c.result <- .C("loadVarSandwich",
                 var.sandwich=as.numeric(var.sandwich),
                 random.effects.count=as.integer(random.effects.count),
                 alt.hyp.count=as.integer(alt.hyp.count),
                 as.integer(n.vc),
                 as.integer(is.poly.fixed),
                 as.integer(is.mg1.fixed),
                 PACKAGE = "multic")
  
  ## Adjust the dimensions of var.sandwich to accurately represent the data
  var.sandwich <-
    array(c.result$var.sandwich,
          dim=c( (n.vc * random.effects.count),
            (n.vc * random.effects.count),
            (alt.hyp.count + 1) ))
  
  ## Remove all major gene covariances in the null hypothesis
  var.sandwich[,,1][var.sandwich[,,1] == -9] <- NA
  
  dimnames(var.sandwich) <- list(random.effects.names, random.effects.names,
                              ibd.names)

  ## The original multic placed the major gene data before the polygene data.
  ## To uniform that data, we swap the first columns and rows, if needed.
  if(!is.poly.fixed && !is.mg1.fixed) {
    tmp.data <- var.sandwich[1, , -1]
    var.sandwich[1, , -1] <- var.sandwich[2, , -1]
    var.sandwich[2, , -1] <- tmp.data
    
    tmp.data <- var.sandwich[, 1, -1]
    var.sandwich[, 1, -1] <- var.sandwich[, 2, -1]
    var.sandwich[, 2, -1] <- tmp.data
  }
  
  return (var.sandwich)
}
