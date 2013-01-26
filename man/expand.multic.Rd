\name{expand.multic}
\alias{expand.multic}
\title{
  Create "bootstrapped" versions of mloci.out and share.out for multic
}
\description{
expand.multic is a utility function to create "bootstrap"ed versions
of mloci.out and share.out
}
\usage{
expand.multic(famids, mloci.out=NULL, share.out=NULL)
}
\arguments{
\item{famids}{
famids is a character or integer vector that specifies the family
order in a "bootstrapped" fashon.  Each index of famids is the famid
(family identifier) from the original dataset not the index of the
family.  An example famids argument would be
\code{famids <- sample(famid, length(unique(famid)), replace = TRUE)}.
IMPORTANT NOTE: This sequence of famids must be the same as that passed to
expand.data.  If they are not, the dataset and the external data
will not match.
}
\item{mloci.out }{
a character value specifying the name of an mloci.out file.  This file
needs to have the famid portion (i.e., the characters before the
hyphen [-]) of the unique id for each entry.
}
\item{share.out }{
a character value specifying the name of an share.out file.  This file
needs to have the famid portion (i.e., the characters before the
hyphen [-]) of the unique id for each entry.
}
}
\value{
a list of two elements.  The first is the name of the new mloci.out
file.  The second element is the name of the new share.out.  Either
element may be \code{NULL} if the respective
input was \code{NULL}.
}
\section{Side Effects}{
the output files are created in the current directory.  If either of
the input files (mloci.out or share.out) were gzip'ed, expand.multic
will gunzip them.  Currently, this is done in their own directory.
However, in the future, this can be done in a temporary.  Also, a
directory named "loci" is temporarily craeted to hold split mloci.out
file.
}
\seealso{
\code{\link{expand.data}}
}
\examples{
\dontrun{
famids <- sample(famid, length(unique(famid)), replace = TRUE)
new.files <- expand.multic(famids, "input/mloci.out", "input/share.out")
mult.obj <- multic( -- your formula, data, famid, etc. here --
                   mloci.out = new.files$new.mloci.out,
                   share.out = new.files$new.share.out)
}
}
% docclass is function
% Converted by Sd2Rd version 52193.
