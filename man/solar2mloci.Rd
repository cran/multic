\name{solar2mloci}
\alias{solar2mloci}
\title{
  Convert a directory of SOLAR-formatted ibd and/or mibd files into a
multic-formatted mloci.out
}
\description{
\code{solar2mloci} reads all of the ibd and mibd files in the given directory,
and creates mloci.out.gz in the specified output directory.
}
\usage{
solar2mloci(directory, phi2, pedindex.out, pedindex.cde,
          ibd.dist = NULL, output.directory = ".",
          delete.fixed.dir = TRUE)
}
\arguments{
\item{directory}{
charater value specifying a path to a directory of SOLAR-formatted
ibd and/or mibd files.  These files are often kept in .gz format.
\code{solar2mloci} will manage this for the user.
}
\item{phi2}{
charater value specifying a path to a SOLAR-formatted phi2 file.
Due to the general size of a typical phi2 file, it is often stored in
.gz format.  \code{solar2mloci} will manage this
for the user.  Whether the user specifies the file with a .gz suffix
or not will not effect how \code{solar2mloci}
operates on the file.
}
\item{pedindex.out }{
character value specifying a path to a SOLAR-formatted pedindex.out
file.  This must be the same file that was output from SOLAR.  It
provides a mapping between the sequential number system assigned by
SOLAR and the original family and individual identifiers.
}
\item{pedindex.cde }{
character value specifying a path to a SOLAR-formatted pedindex.cde
file.  This must be the same file that was output from SOLAR.  This
file describes the format of pedindex.out.  This is necessary to
read pedindex.out correctly.
}
\item{ibd.dist }{
charater value specifying a path to a SOLAR-formatted .dist file that
maps the character marker names of ibd files to their corresponding
numeric centimorgan values.
}
\item{output.directory }{
charater value specifying a path to a directory where the output file
mloci.out.gz will be created. If the directory (including supporting
path) does not exist yet, it will be created.  Also,
\code{solar2mloci} will overwrite mloci.out and
mloci.out.gz if they exist in
\code{output.directory}. 
}
\item{delete.fixed.dir }{
logical flag: if \code{TRUE}, delete the
temporary directory used 
to "fix" the SOLAR-formatted ibds and mibd files.  This is mostly for
debugging 
purposes and is very likely to be not useful to the user community.
}
}
\section{Side Effects}{
Due to write permissions possibly not allowing the user to gunzip
and create files in the specified
\code{directory}, 
\code{solar2mloci} first copies 
\code{directory} and
\code{phi2} to the current directory.
\code{solar2mloci} then 
creates a temporary directory to hold the "fixed," intermediate files
that will be deleted (by default).  Also,
\code{solar2mloci} will overwrite 
mloci.out and mloci.out.gz if they exist in
\code{output.directory}.
}
\seealso{
\code{\link{phi2share}},
\code{\link{solar2multic}}
}
\examples{
\dontrun{
solar2mloci(directory = "ibddir", phi2 = "phi2",
          pedindex.out = "pedindex.out", pedindex.cde = "pedindex.cde",
          ibd.dist = "solar.dist", output.directory = "multicInput",
          delete.fixed.dir = FALSE)

solar2mloci("mibds/chrm10", "solarOutput/phi2.gz",
          "solarOutput/pedindex.out", "solarOutput/pedindex.cde") 
}
}
% docclass is function
% Converted by Sd2Rd version 52193.
