\name{sw2mloci}
\alias{sw2mloci}
\title{
  Convert SimWalk IBD files into a single mloci.out.
}
\description{
sw2mloci converts all IBD files in a given directory into a sinlge
mloci.out needed my multic, altering the centimorgan values if a map
argument is provided.
}
\usage{
sw2mloci(directory, map="", output.directory=".")
}
\arguments{
\item{directory}{
a character object specifying the name of the directory that contains
the SimWalk IBD files.  This can be an absolute or relative path.
}
\item{map}{
a character object specifying the name of a .map file to be used to
modify the centimorgan values in the mloci.out file.  This can be an
absolute or relative path and does not have to be in the same
directory as the parameter \code{directory}.
}
\item{output.directory }{
a character object specifying the name of the directory to put the
finished mloci.out.gz.  This can be an absolute or relative path.  If
the directory does not exist, it will be created.
}
}
\value{
a character object specifying the name of the file created.  In
general, this will be "mloci.out.gz".
}
\section{Side Effects}{
If a file named "mloci.out" or "mloci.out.gz" already exist in the
current directory, \code{sw2mloci} will move
them to "mloci.out.before" or "mloci.out.before.gz" respectively
before doing any calculations.
\code{sw2mloci} also copies the IBD files and
map file (if it is specified) to a temp space.  This is done to bypass
any write permission issues.  This temp space is deleted when the
function is finished.  It also creates a temp space to hold the
intermediate mibd files.  These also will be deleted at the end of the
function.
}
\seealso{
There are similar functions to deal with SOLAR mibds, see 
\code{\link{phi2share}},
\code{\link{solar2mloci}}, and
\code{\link{solar2multic}}.
}
\examples{
\dontrun{
sw2mloci("../otherInput/sw18", "../otherInput/sw18/c18.map")

sw2mloci("sw18")

sw2mloci(".", "sw18/c18.map", output.directory = "multicInput")
}
}
\keyword{}
% docclass is function
% Converted by Sd2Rd version 1.21.
