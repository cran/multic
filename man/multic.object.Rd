\name{multic.object}
\alias{multic.object}
\title{
  a multic object 
}
\description{
Object of class \code{"multic"} returned from
the function \code{multic}. 
}
\section{Generation}{
This class of objects is returned by the
\code{multic} function to represent a fitted
variance components model.
}
\section{Methods}{
Objects of this class have methods for the functions
\code{polygene},
\code{print},
\code{plot},
\code{fitted},
\code{residuals},
and \code{summary}
}
\arguments{
\item{fam.log.liks }{
the log likelihoods and lod scores for each family at each marker
(including the null hypothesis).
\code{fam.log.liks} is a 3-dimensional 
matrix.  The first dimension is indexed by the family identifiers.
The second dimension is indexed by the words
\code{"log.lik"} and
\code{"lod.score"}.  The third dimension is
indexed by the word \code{"null"} and the 
names of the marker file names.  To calculate the family log
likelihoods, \code{calc.fam.log.liks = TRUE}
must be passed 
to \code{multic} via the
\code{...} parameter or a 
\code{multic.control} object.  If 
\code{fam.log.liks} are not calculated, then
\code{fam.log.liks} is a character 
vector providing instructions how to calculate the values.
}
\item{fixed.effects }{
the estimate, standard error, t value, and p value of the fixed
effects for the traits and covariates for the null hypothesis and each
marker.  \code{fixed.effects} is a 3-dimensional
matrix.  The first 
dimension is indexed by the trait and covariate names.  The second
dimension is indexed by the words
\code{"Estmate"},
\code{"Std.err"},
\code{"t.value"}, and 
\code{"p.value"}.  The third dimension is
indexed by the word \code{"null"} and the 
marker file names.
}
\item{polygenic}{
the estimate, standard error, Wald score, Wald score P-value,
heritabilty estimate, standard error of the heritabilty
estimate, and heritably estimate P-value for the variance and
covariance of the polygenic effect of the formula
for the null hypothesis and each marker.
\code{polygenic} is 
a 3-dimensional matrix.  The first dimension is indexed by the letter
\code{"s"} followed by a
\code{1}, \code{2},
etc. for the variance of the first trait, 
second trait, and so on or \code{12},
\code{13}, \code{23},
etc. for the covariance between the 
first and second traits, first and third traits, second and third
traits, and so on.  The second dimension is indexed by the words
\code{"Estimate"},
\code{"Std.err"},
\code{"Wald"},
\code{"W.p.value"},
\code{"h^2"},
\code{"se.h^2"}, and   
\code{"h.p.value"}.  The third dimension is
indexed by the word \code{"null"} and 
the marker file names.
}
\item{major.gene1 }{
the estimate, standard error, Wald score, Wald score P-value,
heritabilty estimate, standard error of the heritabilty
estimate, and heritably estimate P-value for the variance and
covariance of the major gene effect of formula
for the null hypothesis and each marker.
\code{major.gene1} is 
a 3-dimensional matrix.  The first dimension is indexed by the letters
\code{"mg"} followed by a
\code{1}, \code{2},
etc. for the variance of the first trait, 
second trait, and so on or \code{12},
\code{13}, \code{23},
etc. for the covariance between the 
first and second traits, first and third traits, second and third
traits, and so on.  The second dimension is indexed by the words
\code{"Estimate"},
\code{"Std.err"},
\code{"Wald"},
\code{"W.p.value"},
\code{"h^2"},
\code{"se.h^2"}, and   
\code{"h.p.value"}.  The third dimension is
indexed by the word \code{"null"} and 
the marker file names.
}
\item{environmental}{
the estimate, standard error, Wald score, and Wald score P-value
for the variance and covariance of the environmental effect of formula
for the null hypothesis and each marker.  environmental is
a 3-dimensional matrix.  The first dimension is indexed by the letter
\code{"e"} followed by a
\code{1}, \code{2},
etc. for the variance of the first trait, 
second trait, and so on or \code{12},
\code{13}, \code{23},
etc. for the covariance between the 
first and second traits, first and third traits, second and third
traits, and so on.  The second dimension is indexed by the words
\code{"Estimate"},
\code{"Std.err"},
\code{"Wald"}, and
\code{"W.p.value"}.  The third dimension is 
indexed by the word \code{"null"} and the
marker file names. 
}
\item{sibling.sibling }{
the estimate, standard error, Wald score, and Wald score P-value
for the variance and covariance of the sibling to sibling effect of formula
for the null hypothesis and each marker.
\code{sibling.sibling} is 
a 3-dimensional matrix.  The first dimension is indexed by the letters
\code{"sib"} followed by a
\code{1}, \code{2},
etc. for the variance of the first trait, 
second trait, and so on or \code{12},
\code{13}, \code{23},
etc. for the covariance between the 
first and second traits, first and third traits, second and third
traits, and so on.  The second dimension is indexed by the words
\code{"Estimate"},
\code{"Std.err"},
\code{"Wald"}, and
\code{"W.p.value"}.  The third dimension is 
indexed by the word \code{"null"} and the
marker file names.  To receive 
valuable data, the 5th member of
\code{constraints} in the
\code{multic.control} 
object must be set to not \code{"F"} (fixed).
}
\item{parent.parent }{
the estimate, standard error, Wald score, and Wald score P-value
for the variance and covariance of the parent to parent effect of formula
for the null hypothesis and each marker.
\code{parent.parent} is 
a 3-dimensional matrix.  The first dimension is indexed by the letter
\code{"p"} followed by a
\code{1}, \code{2},
etc. for the variance of the first trait, 
second trait, and so on or \code{12},
\code{13}, \code{23},
etc. for the covariance between the 
first and second traits, first and third traits, second and third
traits, and so on.  The second dimension is indexed by the words
\code{"Estimate"},
\code{"Std.err"},
\code{"Wald"}, and
\code{"W.p.value"}.  The third dimension is 
indexed by the word \code{"null"} and the
marker file names.  To receive 
valuable data, the 6th member of
\code{constraints} in the
\code{multic.control} 
object must be set to not \code{"F"} (fixed).
}
\item{parent.offspring }{
the estimate, standard error, Wald score, and Wald score P-value
for the variance and covariance of the parent to offspring effect of formula
for the null hypothesis and each marker.
\code{parent.offspring} is 
a 3-dimensional matrix.  The first dimension is indexed by the letter
\code{"q"} followed by a
\code{1}, \code{2},
etc. for the variance of the first trait, 
second trait, and so on or \code{12},
\code{13}, \code{23},
etc. for the covariance between the 
first and second traits, first and third traits, second and third
traits, and so on.  The second dimension is indexed by the words
\code{"Estimate"},
\code{"Std.err"},
\code{"Wald"}, and
\code{"W.p.value"}.  The third dimension is 
indexed by the word \code{"null"} and the
marker file names.  To receive 
valuable data, the 7th member of
\code{constraints} in the
\code{multic.control}
object must be set to not \code{"F"} (fixed).
}
\item{log.liks }{
the log likelihood, centimorgan distance, log likelihood status, and
lod score and P-value for the null hypothesis and each marker.
\code{log.liks} is a
\code{data.frame}.  The row names are
\code{"null"} and the markder 
file names.  The column names are
\code{"log.likelihood"},
\code{"distance"}, 
\code{"log.lik.status"},
\code{"lod.score"}, and
\code{"p.value"}.  The log likelihood 
status represents whether the log likelihood converged before the
maximum interations allowed or not and have the values of either
\code{"converg"} or
\code{"non-converg"}. 
}
\item{var.fixed }{
the variance of the fixed effects of the traits and covariates for the
null hypothesis and each marker.
\code{var.fixed} is a 3-dimensional 
matrix.  The first and second dimensions are indexed by the trait and
covariate names.  The third dimension is indexed by the word
\code{"null"} and the 
marker file names.
}
\item{var.random }{
the variance of the polygenic, major gene, and
environmental effects for the null hypothesis and each marker.
\code{var.random} is a 3-dimensional matrix.
The first and second dimensions 
are indexed as described by the polygenic, major.gene1, and
environmental components above.  The third dimension is indexed by the
word \code{"null"} and the marker file names.
}
\item{var.sandwich }{
a more precise variance estimator after using a sandwich estimator
approach.  This is only calculated if the multic object represents a
univariate model.  \code{var.sandwich} is a
3-dimensional matrix.  The first 
and second dimensions are indexed by
\code{"s1"},
\code{"mg1"}, and
\code{"e1"}.  The third 
dimension is indexed by the word \code{"null"}
and the marker file names. 
}
\item{cors}{
the Pearson, Spearman, genetic, environmental, and phenotypic
correlations.  \code{cors} is a list made up of
the components \code{"pearson"}, 
\code{"spearman"},
\code{"genetic"},
\code{"environment"}, and
\code{"phenotype"}.  Both
\code{"pearson"} 
and \code{"spearman"} are their respective
correlations between the traits and 
covariates.  They are 2-dimensional matrices indexed by the trait and
covariate names.  \code{"genetic"},
\code{"environment"}, and
\code{"phenotype"} are the 
respective correlations between the
\code{polygenic} and
\code{environmenal} 
estimates.  They are 2 dimensional matrices.  The first dimension is
indexed by the word \code{"null"} and the
marker file names.  The second 
dimension is indexed as described by the covariance portions of the
\code{polygenic} and
\code{environmenal} components above. 
}
\item{v.matrices }{
the variance-covariance matrix of the trait (y) that incorporates the
polygenic, major gene, shared common environment, and error matrices.
\code{v.matrices} is a 2-dimensional matrix.
The first dimension is indexed 
by the family identifier (\code{famid}) values.
The second dimension is 
indexed by the word \code{"null"} and the
marker file names.  Currently, 
there are no individual identifiers on each of the V matrices.  If the
V matrices are not calculated, then
\code{v.matrices} is a character vector 
providing instructions how to calculate the values. 
}
\item{residuals}{
the observed values minus the fitted values of the trait (y) divided by
the square root of the V matrix for each family.  If the residuals are
not calculated, then \code{residuals} is a
character vector providing 
instructions how to calculate the values.
}
\item{descriptives}{
the total individuals used, mean, standard deviation, minimum,
maximum, kurtosis, and skewness for each trait and covariate.
}
\item{counts}{
various counts of the total number of pedigrees, people, females,
males, and so on.  This is mostly for passing data for
\code{print} and 
\code{summary} to display and is very likely to
be not useful to the user community.
}
\item{call}{
how \code{multic} was called.  call is a call
object. 
}
\item{R.sq }{
the proportion of variance due to the covariates.
}
\item{metadata}{
a list of useful data like \code{start.time},
\code{finish.time},
\code{call},
\code{epsilon}, 
\code{trait.count},
\code{iterations},
\code{null.initial.values},
\code{method}, etc. 
}
}
\seealso{
\code{\link{multic}}
}
\keyword{classes}
\keyword{methods}
% docclass is function
% Converted by Sd2Rd version 1.21.
