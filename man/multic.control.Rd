\name{multic.control}
\alias{multic.control}
\title{
  Set control parameters for multic 
}
\description{
Allows users to alter the default behavior of multic
}
\usage{
multic.control(epsilon = 1e-5,
               max.iterations = 50,
               boundary.fix = TRUE,
               constraints = c("E", "E", "E", "E", "F", "F", "F"),
               initial.values = NULL,
               save.output.files = FALSE,
               method = c("multic", "leastsq", "maxfun", "emvc"),
               calc.fam.log.liks = FALSE,
               calc.residuals = FALSE,
               keep.input = calc.residuals)
}
\arguments{
\item{epsilon}{
a numeric value specifying the convergence threshold.  When the
difference of an iteration's loglikelihood and the previous
iteration's loglikelihood are less than
\code{epsilon}, the value has 
"converged".
}
\item{max.iterations}{
an integer value specifying the maximum number of iterations
\code{multic} 
will take to converge during the polygenic and sporadic model
calculations.
}
\item{boundary.fix }{
logical flag: if \code{TRUE}, then the variances
generated will be fixed to \code{0}
and no longer estimated when they become less than
\code{0.00001} (\code{1e-5}).
}
\item{constraints}{
a character vector of length seven (7) specifying the constraints on
the random effects variance components.  Each
value of the vector needs to be either
\code{"E"} - `E'stimate the variance and 
covariance, \code{"C"} - estimate the variance
and `C'onstrain the covariance, 
or \code{"F"} - `F'ix the variance and
covariance to \code{0}.  Each index of 
constraints corresponds to (in this exact order) mu, polygene, major
gene, environment, sibling-sibling, parent-parent, and
parent-offspring.
}
\item{initial.values }{
numeric vector: use the specified initial values instead of
calculating them automatically.  This vector has a very specific
length and order.  If n is the number of traits and m is ( n + (n-1) +
(n-2) + ... + 1 ), then the length must be n + 6 * m.  So for one
trait (univariate), the length must be 7, for two traits (bivariate),
20, and so on. The position of the values in the vector is important
as well.  The first n terms are the mu starting values.  The next
starting values come in chunks of m.  The next m values are the
polygenic starting values, followed by major gene, environmental,
sibling-sibling, parent-parent, and parent-offspring starting values.
The metadata\$null.initial.values contains the placement of the
starting values.  You can use this to verify your order is correct.
}
\item{save.output.files }{
logical flag: if \code{TRUE}, then the multiple
temporary output files \code{multic}
generates are not removed.  This is mostly for debugging purposes and
is very likely to be not useful to the user community.
}
\item{method}{
a character value specifying the method to use in fitting the model.
Possible values include \code{"multic"}
(default), \code{"leastsq"},
\code{"maxfun"}, and  
\code{"emvc"} (all case insensitive).
}
\item{calc.fam.log.liks }{
logical flag: if \code{TRUE}, then the family
log likelihoods will be 
returned in the \code{multic} object.  WARNING:
This significantly increases 
the size of the returned \code{multic} object.
}
\item{calc.residuals }{
logical flag: if \code{TRUE}, then the
residuals will be calculated and Y 
beta differences and V matrix data will be returned in the
\code{multic} 
object.  WARNING: This dramatically increases the size of the returned
\code{multic} object.
}
\item{keep.input }{
logical flag: if \code{TRUE}, then the traits
and covariates will be 
saved in the \code{metdata} list of the
\code{multic} 
object.  Since the input is 
needed during special residual calculations, its default value is that
of \code{calc.residuals}.
}
}
\value{
a list that is designed to be supplied as a control argument to
\code{multic}. The values for
\code{multic.control} can 
be supplied directly in a 
call to \code{multic} (via the
\code{...} parameter).  These values are then 
filtered through \code{multic.control} inside \code{multic}.
}
\seealso{
\code{\link{multic}},
\code{\link{multic.object}}
}
\examples{
\dontrun{
## The following calls to multic are equivalent 
multic(formula, data, control = multic.control(calc.fam.log.liks = TRUE,
                                               calc.residuals = TRUE))  
multic(formula, data, calc.fam.log.liks = TRUE, calc.residuals = TRUE)
}
}
% docclass is function
% Converted by Sd2Rd version 52193.
