\name{phi2share}
\alias{phi2share}
\title{
  Convert a SOLAR-formatted phi2 file into a multic-formatted share.out file.
}
\description{
\code{phi2share} reads in the specified files
and generates a multic-formatted share.out file.  share.out contains
unique identifiers, expected shared genetic material, and sibling,
spousal, and parent-offspring true/false values.
}
\usage{
phi2share(phi2, pedigree.file, pedindex.out, pedindex.cde,
          output.directory=".")
}
\arguments{
\item{phi2}{
a character value to specify the location of a SOLAR-formatted phi2
(or similarly formatted) file.  Due to the general size of a typical
phi2 file, it is
often stored in .gz format.  \code{phi2share}
will manage this for the user. 
Whether the user specifies the file with a .gz suffix or not will not
effect how \code{phi2share} operates on the
file.  
}
\item{pedigree.file }{
a character value to specify the location of a .ped (or
similarly formatted) file.  This file must have a header of
\code{famid}, \code{id},
\code{fa}, \code{mo},
and \code{sex} (case insensitive).  The file
must also be comma separated. 
}
\item{pedindex.out }{
a character value to specify the location of a pedindex.out file.
This file must be the same that was output from SOLAR.  It provides a
mapping between the sequential number system assigned by SOLAR and the
original family and individual identifiers.
}
\item{pedindex.cde }{
a character value to specify the location of a pedindex.cde file.
This file must be the same that was output from SOLAR.  This file
describes how pedindex.out is organized.  This is necessary to read
pedindex.out correctly.
}
\item{output.directory }{
a character value specifying which directory to place the output
share.out.  If \code{output.directory}
(including supporting path) does not exist yet, 
it will be created.  The default directory is the current directory.
}
}
\section{Side Effects}{
\code{phi2share} creates a local copy of,
gunzip's, and removes the copy of \code{phi2}.
It also will overwrite share.out and share.out.gz
if they exist in \code{output.directory}. 
}
\seealso{
\code{\link{solar2mloci}},
\code{\link{solar2multic}}
}
\examples{
\dontrun{
phi2share(phi2 = "phi2",
          pedigree.file = "chrom18.ped",
          pedindex.out = "pedindex.out",
          pedindex.cde = "pedindex.cde",
          output.directory = "multicInput")
phi2share("solarOutput/phi2.gz",
          "solarOutput/chrom18.ped",
          "solarOutput/pedindex.out",
          "solarOutput/pedindex.cde")
}
}
% docclass is function
% Converted by Sd2Rd version 52193.
