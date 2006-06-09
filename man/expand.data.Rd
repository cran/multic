\name{expand.data}
\alias{expand.data}
\title{
  Create a "bootstrapped" version of a dataset to be used in multic.
}
\description{
When using multic to bootstrap over families, an appropriate data set
is needed.  By providing a random set of famids, expand.data creates
such a dataset.
}
\usage{
expand.data(famids, d.frame)
}
\arguments{
\item{famids}{
famids is a character or integer vector that specifies the family
order in a "bootstrapped" fashon.  Each index of famids is the famid
(family identifier) from the original dataset not the index of the
family.  An example famids argument would be
\code{famids <- sample(famid, length(unique(famid)), replace = TRUE)}.
IMPORTANT NOTE: This sequence of famids must be the same as that passed to
expand.multic.  If they are not, the dataset and the external data
will not match.
}
\item{d.frame }{
the data.frame that holds the family structure and phenotype data.
This should be the dataset that was used to sample famid.
}
}
\value{
a data.frame that contains the bootstrapped version of the input dataset
}
\seealso{
\code{\link{expand.multic}}
}
\examples{
\dontrun{
famids <- sample(famid, length(unique(famid)), replace = TRUE)
expaned.ped.phen <- expand.data(famids, ped.phen)
}
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
