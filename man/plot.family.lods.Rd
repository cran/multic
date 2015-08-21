\name{plotFamilyLods}
\alias{plotFamilyLods}
\title{
  Plot family contribution to peak LOD
}
\description{
Plot the top N families that contribute to the peak LOD score as obtained
from a multic object
}
\usage{
plotFamilyLods(x, npeakfams=5, title=NULL, title.cex=0.75, legend=T,
               legend.loc=c("left", "middle", "right", "extra"),
               type=c("top", "total", "proportion", "all"), ...)
}
\arguments{
\item{x}{
a \code{multic} object
}
\item{npeakfams}{
number of top families to plot
}
\item{title}{
title for the plot
}
\item{title.cex }{
character size for the title
}
\item{legend}{
logical (default = TRUE) indicating whether a legend is displayed or not
}
\item{legend.loc }{
character string indicating the general legend location, if legend =T.
"left" indicates the top left, "right" indicates the top right, "middle"
indicates the upper middle, and "extra" indicates the legend should go on a
separate plot.
}
\item{type}{
character string indicating which plot should be shown.  "top" shows the
top npeakfams families, "total" shows the top families plus the overall
total, "proportion" shows the percentage for the top families of the total
LOD score, and "all" shows all 3 plots.
}
\item{\dots }{
additional parameters to alter the default behavior of the plot.
}
}
\value{
the lod scores for the top families
}
\details{
Based on the peak total lod score for the multic object, determine which
families contribute most to that peak.  Then display their contribution
across the entire area that was used to fit the multic object.
}
\note{
The multic object must have been fit with the option
"calc.fam.log.liks=TRUE".  Additionally, the function will not work on
polygenic multic objects.
}
\seealso{
\code{\link{multic.object}}, 
\code{\link{multic}}
}
\examples{
\dontrun{
mult10 <- multic(sys.avg ~ sex + agexam + agexam^2, data=d10,
                 famid=famid, id=id, dadid=fa, momid=mo, sex=sex,
                 mloci='multicInput/mloci.out.gz', share='multicInput/share.out.gz',
                 calc.fam.log.liks=T)

plotFamilyLods(mult10, npeakfams=3, plot="total")
}
}
% docclass is function
% Converted by Sd2Rd version 52193.
