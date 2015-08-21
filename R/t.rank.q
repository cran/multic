## t.rank is a transformation function.  It is intended to be wrapped around
## the formula components in a call to multic.
##
## Example:
## multic(t.rank(trig) ~ age + bmi, ...)
##
## t.rank was originally written by John Paul Bida
## $Id: t.rank.q,v 1.2 2014/03/14 17:59:44 m015733 Exp $
tRank <- function(x) {
  t1 <- x
  t2 <- x[!is.na(x)]
  r <- qnorm(rank(t2)/(length(rank(t2))+1))
  t1[!is.na(t1)] <- r
  t1
}

## Method renamed to tRank because it's not S3.
## Keep alias for backward compatibility.
t.rank <- tRank
