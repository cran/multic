###############################################################################
## Title: load.inv.exp.sec.der.random
## Description: load.inv.exp.sec.der.random allocates memory to hold the values
##              of the inverse of the expected second derivative matrix created
##              in multic.  It receives the information via the multic output
##              file 'invExpSecDerRandom.log'.
## Input: random.effects.count - random.effects.count is a scalar value
##                               representing the number of random effects per
##                               effect (polygene, major gene) that the last
##                               run of multic output.
##        alt.hyp.count - alt.hyp.count is the number of alternative hypotheses
##                        that the last run of mutlic calculated.
##        random.effects.names - a character vector containing the names
##                               (dimnames or labels) for the first two
##                               dimensions of inv.exp.sec.der.random.
##        ibd.names - a character vector containing the names (dimnames or
##                    labels) of the names of the ibd files that will be the
##                    names of the third dimension of inv.exp.sec.der.random.
## Ouput: inv.exp.sec.der.random - inv.exp.sec.der.random is a three
##                                 dimensional array that contains the filled
##                                 values of the inverse of the second
##                                 derivative under the null hypothesis and
##                                 each of the alternative  hypotheses.
## Side Effects: load.inv.exp.sec.der.random calls a compiled C/C++ routine
##               named 'loadInvExpSecDerRandom'.  If the file to be read is not
##               readable or not present, Splus will terminate.
##
##               As of 02-07-2004, if the file 'invExpSecDerRandom.log' doesn't
##               exist, a NULL object will be returned, the program will not
##               terminate.
## Author: Eric Lunde 10-24-03
## Updates: (Date, Modified By, Modification Description)
## 11-28-03, Eric Lunde, 11-28-03, I removed the parameters trait.count and
##           replaced it with random.effects.count.  I made a decision that
##           this effect count value was the important (useful) value.  So it
##           is calculated early in the program and the value is passed
##           through, rather than each function calculating the value again
##           and again.
## 02-07-2004, Eric Lunde, I added two more input parameters:
##             fixed.effects.names and ibd.names.  They are the dimnames for
##             the return value inv.exp.sec.der.random.  Also I added a file
##             existing test at the beginning to see if the file
##             'invExpSecDerRandom.log' exists or not.
## 03-31-2005, Eric Lunde, The number of environmental variance (covariance)
##             values no longer is dependent upon whether multic is run in
##             multivariate or longitudinal mode (it now always operates as if
##             in multivariate mode).  Therefore I've removed the
##             environment.count variable and replaced it with the
##             random.effects.count which now always is the correct value.
###############################################################################
load.inv.exp.sec.der.random <- function(random.effects.count, alt.hyp.count,
                                        random.effects.names, ibd.names) {
  if( !file.exists("invExpSecDerRandom.log") ) {
    cat(paste("\nThe multic output file 'invExpSecDerRandom.log' is not ",
              "present.\n",
              "The multic object will not contain data concerning the inverse",
              " of the expected second derivative for the random effects.\n",
              sep=""))
    return (NULL)
  }
  
  ## Verify that the arguments are of the correct type and range
  random.effects.count <- as.integer(random.effects.count)
  if(random.effects.count < 0) {
    stop(paste("random.effects.count has the value '", random.effects.count,
               "', but must be positive.  load.effects.s key 45"))
  }
  
  alt.hyp.count <- as.integer(alt.hyp.count)
  if(alt.hyp.count < 0) {
    stop(paste("alt.hyp.count must be positive.  Current value is  -",
               alt.hyp.count, "\nload.inv.exp.sec.der.random.s key 40"))
  }

  ## inv.exp.sec.der.random is the variable that will hold the inverse of the
  ## expected second derivative for the null hypothesis and all of the
  ## alternative hypotheses
  inv.exp.sec.der.random <-
    numeric( (3 * random.effects.count)
            * (3 * random.effects.count)
            * (alt.hyp.count + 1 ) )

  c.result <- .C("loadInvExpSecDerRandom",
                 inv.exp.sec.der.random=as.numeric(inv.exp.sec.der.random),
                 random.effects.count=as.integer(random.effects.count),
                 alt.hyp.count=as.integer(alt.hyp.count),
                 PACKAGE ="multic" )

  ## Adjust the dimensions of inv.exp.sec.der.random
  inv.exp.sec.der.random <-
    array(c.result$inv.exp.sec.der.random,
          dim=c( (3 * random.effects.count),
            (3 * random.effects.count),
            (alt.hyp.count + 1 ) ))
  
  ## Eliminate the major gene variances in the null hypothesis.  The '-9' is
  ## the value that C++ used to indicate a missing value.  Also eliminate the
  ## variances from effects that have been fixed internally or constrained
  ## covariances.
  inv.exp.sec.der.random[inv.exp.sec.der.random == -9] <- NA
  
  dimnames(inv.exp.sec.der.random) <- list(random.effects.names,
                                           random.effects.names,
                                           ibd.names)

  return (inv.exp.sec.der.random)
}
