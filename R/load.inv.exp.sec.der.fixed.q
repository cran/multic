###############################################################################
## Title: load.inv.exp.sec.der.fixed
## Description: load.inv.exp.sec.der.fixed is responsible for collecting the
##              informtion contained in the inverse of the expected second
##              derivative matrix for the fixed effects.  It receives the
##              information via the multic output file 'invExpSecDerFixed.log'.
## Input: fixed.effects.count - fixed.effects.count is a scalar value
##                              representing the number of fixed effects that
##                              the last run of multic output.
##        locus.count - locus.count is the number of loci that the last run of
##                      mutlic calculated with (the number of alternative
##                      hypotheses plus the null hypothesis)
##        fixed.effects.names - a character vector containing the names
##                              (dimnames or labels) for the first two
##                              dimensions of inv.exp.sec.der.fixed
##        ibd.names - a character vector containing the names (dimnames or
##                    labels) of the names of the ibd files that will be the
##                    names of the third dimension of inv.exp.sec.der.fixed.
## Ouput: A three dimensional array of numeric values representing the inverse
##        of the expected second derivative for the fixed effects for each of
##        the loci presented to multic.
## Side Effects: load.inv.exp.sec.der.fixed calls a compiled C/C++ routine
##               named 'loadInvExpSecDerRandom'.  If the file to be read is not
##               readable or not present, Splus will terminate.
##
##               As of 02-07-2004, if the file 'invExpSecDerFixed.log' doesn't
##               exist, a NULL object will be returned, the program will not
##               terminate.
## Author: Eric Lunde 10-24-03
## Updates: (Date, Modified By, Modification Description)
## 11-28-03, Eric Lunde, 11-28-03, I removed the parameters trait.count and
##           covariate.count and replaced them with fixed.effects.count.  I
##           made a decision that this fixed effect count value was the more
##           important (useful) value.  So it was calculated early in the
##           program and passed through, rather than each function calculating
##           the value again and again.
## 02-07-2004, Eric Lunde, I added two more input parameters:
##             fixed.effects.names and ibd.names.  They are the dimnames for
##             the return value inv.exp.sec.der.fixed.  Also I added a file
##             existing test at the beginning to see if the file
##             'invExpSecDerFixed.log' exists or not.
###############################################################################
load.inv.exp.sec.der.fixed <- function(fixed.effects.count, loci.count,
                                       fixed.effects.names, ibd.names) {
  if( !file.exists('invExpSecDerFixed.log') ) {
    cat(paste("\nThe multic output file 'invExpSecDerFixed.log' is not ",
              "present.\n",
              "The multic object will not contain data concerning the inverse",
              " of the expected second derivative for the fixed effects.\n",
              sep=""))
    return (NULL)
  }
  
  ## Verify that the arguments are of the correct type and range
  fixed.effects.count <- as.integer(fixed.effects.count)
  if(fixed.effects.count < 0) {
    stop(paste("fixed.effects.count has the value '", fixed.effects.count,
               "', but must be positive. load.inv.exp.sec.der.fixed.s key 35",
               sep=""))
  }
  
  loci.count <- as.integer(loci.count)
  if(loci.count < 0) {
    stop(paste("loci.count has the value '", loci.count, "', but must be ",
               "positive.  load.inv.exp.sec.der.fixed.s key 42", sep=""))
  }
  
  inv.exp.sec.der.fixed <- read.table('invExpSecDerFixed.log')

  # Adjust the dimensions of inv.exp.sec.dir.fixed
  inv.exp.sec.der.fixed <- array(inv.exp.sec.der.fixed[[1]],
                                 c(fixed.effects.count, fixed.effects.count,
                                   loci.count))
  
  dimnames(inv.exp.sec.der.fixed) <- list(fixed.effects.names,
                                          fixed.effects.names,
                                          ibd.names)
  
  return (inv.exp.sec.der.fixed)
}
