\name{tRank}
\alias{tRank}
\title{
  Empirical Normal Quantile transformation
}
\description{
Transform a vector using the empirical normal quantile distribution (also
called the van der Waerden normal scores).
}
\usage{
tRank(x)
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
\examples{
\dontrun{
stem(tRank(1:100))
}
}
% docclass is function
% Converted by Sd2Rd version 52193.
