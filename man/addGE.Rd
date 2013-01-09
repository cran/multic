\name{addGE}
\alias{addGE}
\title{
  Assess combinations of univariate multic objects
}
\description{
Determine whether there is any evidence that running a multivariate multic
model will significantly increase the evidence of a genetic effect.
}
\usage{
addGE(multic.objs, combine=2, plotit=FALSE, ibd.dist, statistic=c("lrt",
"wald"), legend=TRUE, ylim=NULL, ...)
}
\arguments{
\item{multic.objs }{
A list of 1-trait multic objects. 
}
\item{combine}{
Indicate how many traits should be examined together.  The program will
then look at all N traits choose 'combine'.
}
\item{plotit}{
Logical, default=FALSE.  If TRUE, a LOD plot is generated with a separate
line for each combination of traits.
}
\item{ibd.dist }{
The default is to use the distances from the first multic object.  This
options allows the user to provide a different set of distances.
}
\item{statistic}{
Character, default="lrt".  This determines whether the Wald statistic
(MG/SE) or the LRT is used when combining the traits. 
}
\item{legend}{
Logical, default=TRUE.  If TRUE and if plotit=TRUE then a legend is
automatically provided.  
}
\item{ylim}{
The extremes of the range of the y-axis to pass to the plot function.
}
\item{\dots }{
Allows for graphical parameters to be passed to the plot function (only
applicable when plotit=TRUE).
}
}
\value{
A data frame is returned if the function is assigned to new object.
Included are the various combinations (listed in order 1-N), the Chi-square
statistic, the p-value, the distance, and the LOD score.
}
\section{Side Effects}{
If plotit=T, a plot is generated on the current graphical device.
}
\section{References}{
M. de Andrade, C. Olswold, J.P. Slusser, L.A.Tordsen, E.J. Atkinson,
K.G. Rabe, and S.L.Slager.
Identification of genes involved in alcohol consumption and cigarette
smoking. BMC Genetics, 6:S112, 2005.
}
\seealso{
multic, gene.eff
}
\examples{
\dontrun{
add2 <- addGE(list(bmi10, dia10, sys10), combine = 2, plotit = T, ylim=c(0,8), legend=F)
add3 <- addGE(list(bmi10, dia10, sys10), combine = 3, plotit = F)
lines(add3$cM, add3$lod, col=4, lwd=2, lty=4)
key(corner=c(0,1), lines=list(lwd=2, col=1:4, lty=1:4),
    text=list(c('BMI-Dia','BMI-Sys','Dia-Sys','BMI-Dia-Sys'), col=1:4))
}
}
% docclass is function
% Converted by Sd2Rd version 52193.
