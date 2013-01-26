\name{summary.multic}
\alias{summary.multic}
\title{
  Use print() on a multic object
}
\description{
Produces a summary of a fitted multic object.
}
\usage{
\method{summary}{multic}(object, ...)
}
\arguments{
\item{object}{
an object of class "multic", usually a result of a call to "multic".
}
\item{\dots }{
further arguments passed to or from other methods.
Currently \code{...} only exists to pass 'R CMD
check' tests.
}
}
\value{
a list is returned with the following components.
<s-arg name="call">
as contained on object
</s-arg>
<s-arg name="max.lod.score">
maximum lod score from the multic object
 </s-arg>
<s-arg name="max.lod.locus">
IBD file that corresponds to the maximum lod score
 </s-arg>
<s-arg name="max.lod.centimorgan">
position where the maximum lod score was found
 </s-arg>
<s-arg name="n">
default=5. Number of families to examine for the top families
 </s-arg>
<s-arg name="top.n.families">
ids of the top n families as they contribute to the maximum lod score
 </s-arg>
<s-arg name="centimorgan.close.to.peak">
The minimum and maximum positions (cM) that produced a lod score 
greater than the "maximum - 1" and are contiguous to max.lod.locus.
 </s-arg>
}
\details{
This is a method for the function summary() for objects inheriting from
class multic.  See "summary" for the general behavior of this function.
 
}
\seealso{
\code{\link{multic}}, \code{\link{multic.object}}
}
\examples{
\dontrun{summary(fit)}
}
\keyword{multic, summary}
% docclass is function
% Converted by Sd2Rd version 52193.
