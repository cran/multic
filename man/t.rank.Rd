\name{t.rank}
\alias{t.rank}
\title{
  Empirical Normal Quantile transformation
}
\description{
Transform a vector using the empirical normal quantile distribution (also
called the van der Waerden normal scores).
}
\usage{
t.rank(x)
}
\arguments{
\item{x}{
 a numerical vector 
}
}
\value{
A numerical vector that are approximately normally distributed.
}
\section{References}{
Lehman E.L. 1975.  Nonparametrics: Statistical Methods Based On
Ranks. Holden-Day, Inc., San Francisco, CA (page 97)
}
\seealso{

}
\examples{
\dontrun{
stem(t.rank(1:100))
}
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
