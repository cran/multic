\name{multic}
\alias{multic}
\title{
  Create a multic object
}
\description{
Calculate the polygenic and major gene models for quantitative trait
linkage analysis using variance components approach.
}
\usage{
multic(formula, data = sys.parent(), famid, id, dadid, momid, sex,
       mloci.out = NULL, share.out = "kinship", longitudinal = FALSE,
       subset = NULL, ascertainment = NULL,
       control = multic.control(...), ...)
}
\arguments{
\item{formula}{
a formula object, with the traits on the
left of a \code{~} (tilde) operator 
and the covariates, separated by \code{+}
operators, on the right.  The 
traits may be a single numeric vector or a matrix.  Commonly, traits
are aggragated together using the \code{cbind}
command.  See the Examples 
section for examples. 
}
\item{famid}{
integer, numeric, or character vector specifying each individual's
family identifier.  Members of the same family must have the same
\code{famid} and each family must have a unique
\code{famid}.  Any missing data 
will result in an error message and the termination of
\code{multic}. 
}
\item{id}{
integer, numeric, or character vector specifying each individual's
identifier.  Members of the same family must have a unique
\code{id} within 
the family.  \code{id} does not have to be
universally unique among all 
individuals.  Any missing data will result in an error
message and the termination of \code{multic}. 
}
\item{dadid}{
integer, numeric, or character vector specifying each individual's
father identifier.  This father identifier must have the same
\code{famid} as 
the individual.  Any missing data will result in an
error message and the termination of \code{multic}. 
}
\item{momid}{
integer, numeric, or character vector specifying each individual's
mother identifier.  This mother identifier must have the same
\code{famid} as 
the individual.  Any missing data will result in an
error message and the termination of
\code{multic}. 
}
\item{sex}{
integer, numeric, or character vector specifying each individual's
sex.  Acceptable forms of sex-coding are
\code{"M"}, \code{"m"}, or
\code{1} for male and
\code{"F"}, \code{"f"}, or
\code{2} for female. Any missing data will 
result in an error message and the termination of
\code{multic}.
}
\item{data}{
\code{data.frame} in which to interpret the
variables named in \code{formula}, 
\code{famid}, \code{id},
\code{dadid},
\code{momid},
\code{sex},
\code{subset}, and
\code{ascertainment} arguments.  If 
\code{data} is missing, the 
variables in \code{formula} should be in the
search path. 
}
\item{mloci.out }{
a character value specifying a path to an mloci.out (or
similarly formatted) file.  Specifying a non-empty mloci.out file will
allow multic to calculate sporadic models using the ibd (identity by
decent) information in the mloci.out file.  Due to the general size of
mloci.out, it is often stored in .gz format.
\code{multic} will manage this
for the user.  Whether the user specifies an mloci.out file with a .gz
suffix or not will not effect how \code{multic}
operates on the file.  See
\code{solar2mloci} for more details.
}
\item{share.out }{
a character value specifying a path to a share.out (or
similarly formatted) file.  This file contains the amount of
genetic material shared between each family member pair based on
family structure only.  It also contains boolean values to indicate
whether two familiy members have a sibling-sibling, parent-parent, or
parent-offspring relationships.  Due to the general size of
share.out, it is often stored in .gz format.
\code{multic} will manage this
for the user.  Whether the user specifies the file wiht a .gz suffix
format or not will not effect how \code{multic}
operates on the file.  See 
\code{phi2share} for more details.
}
\item{longitudinal}{
logical flag: if \code{TRUE}, then fomula will
be interpreted as a 
longitudinal model.  In this case, the fomula argument requires
special formatting as described.  The number of traits on the
left side of the \code{~} (tilde) is the number
of time-points for \code{multic} to 
analyze.  The number of covariates on the right side of the
\code{~} (tilde) 
must be a multiple of the number of traits on the left side.  That
multiple is the amount of covariates to analyze at each time-point.
All covariates for the first time point must be specified before any
of the second, all second before any third, etc.  See the Examples
section for examples.
}
\item{subset}{
a logical vector specifying which subset of the rows in
\code{data} to 
use in the fit.
}
\item{ascertainment}{
vector specifying each individual's ascertainment (effected) status. 
Acceptable forms of ascertainment are \code{T},
\code{TRUE}, or \code{1}
for a proband (effected) and \code{F},
\code{FALSE}, or
\code{0} for a non-proband (non-effected 
person).
}
\item{control}{
\code{list} of iteration and algorithmic
constants. See 
\code{multic.control} for their names and
default values. These can also be given directly as arguments to
\code{multic} itself, instead of through a
\code{multic.control} object.  If control 
is specified, the arguments specified in
\code{...} will not be used. 
}
\item{\dots }{
further arguments passed to
\code{multic.control} to alter
\code{multic}'s default behavior.
}
}
\value{
an object of class \code{"multic"}.  See
\code{multic.object} for more details.
}
\section{Side Effects}{
Many temporary files are created during
\code{multic}'s execution.  These 
files are deleted afterwards (by default).  If they are not deleted
(due to a crash or some other unexpected action), use the included
function clean() to delete them.  Also,
\code{multic} copies, gunzip's, and 
removes the copies of share.out and mloci.out (if specified).
}
\details{
See the technical report.
}
\section{References}{
Amos, C. I. (1994). "Robust variance-components approach for assessing
genetic linkage in pedigrees." American Journal of Human Genetics
54(3): 535-543. 

Almasy, L. and J. Blangero (1998). "Multipoint quantitative-trait
linkage analysis in general pedigrees." American Journal of Human
Genetics 62(5): 1198-1211. 
}
\seealso{
\code{\link{multic.object}},
\code{\link{multic.control}},
\code{\link{phi2share}},
\code{\link{solar2mloci}},
\code{\link{solar2multic}},
\code{\link{sw2mloci}}
}
\examples{
\dontrun{
# Call multic with a univariate formula with two covariates and no
# markers (no mloci.out argument).
fit.ibd.uni <- multic(k.trig ~ sex.x + agexam,
                      data = ped.phen.data,
                      famid, id, fa, mo, sex.x,
                      share.out = 'multicInput/share.out')

# Call multic with a bivariate formula with three covariates, no
# markers (no mloci.out argument), and calculate the family log
# likelihoods.
fit.ibd.bi <- multic(cbind(k.trig, k.chol) ~ sex.x + agexam + agexam2,
                     data = ped.phen.data,
                     famid, id, fa, mo, sex.x,
                     share.out = 'multicInput/share.out',
                     calc.fam.log.liks = TRUE)

# Call multic with a longitudinal formula with six covariates letting
# the kinship library calculate the share.out argument.
long.fit <- multic(cbind(sbpA, sbpB, sbpC) ~
                     sexA + ageA + bmiA + generA + ageAg + smkA +
                     sexB + ageB + bmiB + generB + ageBg + smkB +
                     sexC + ageC + bmiC + generC + ageCg + smkC,
                   data = long.data,
                   famid, id, dadid, momid, sex,
                   longitudinal = TRUE)
}
}
% docclass is function
% Converted by Sd2Rd version 52193.
