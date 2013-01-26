\name{solar2multic}
\alias{solar2multic}
\title{
  Convert SOLAR-formatted output into multic-formatted mloci.out and
share.out
}
\description{
\code{solar2multic} is a utility function to
convert the ibd and mibd files  (identity by 
descent) files created by SOLAR into the multic input file mloci.out
and convert the phi2 created by SOLAR into the multic input file
share.out.
}
\usage{
solar2multic(phi2, pedigree.file, pedindex.out, pedindex.cde,
             ibd.directory, ibd.dist = NULL, output.directory = ".", 
             delete.fixed.dir = TRUE)
}
\arguments{
\item{phi2}{
character value specifying a path to a SOLAR-formatted phi2 file.
Due to the general size of a typical phi2 file, it is often stored in
.gz format. 
\code{solar2multic} will manage this for the
user. Whether the user specifies 
the file with a .gz suffix or not will not effect how
\code{solar2multic} 
operates on the file.  
}
\item{pedigree.file }{
character value specifying a path to a SOLAR-formatted pedigree
structure file (.ped).  This file must have a header of
\code{famid}, \code{id},
\code{fa}, \code{mo},
and \code{sex} (case insensitive).  The file
must also be comma separated. 
}
\item{pedindex.out }{
character value specifying a path to a SOLAR-formatted pedindex.out
file.  This file must be the same that was output from SOLAR.  It
provides a mapping between the sequential number system assigned by
SOLAR and the original family and individual identifiers.
}
\item{pedindex.cde }{
character value specifying a path to a SOLAR-formatted pedindex.cde
file.  This file must be the same that was output from SOLAR.  This
file describes how pedindex.out is organized.  This is necessary to
read pedindex.out correctly.
}
\item{ibd.directory }{
character value specifying a path to a directory containing
SOLAR-formatted ibd and/or mibd files.
}
\item{ibd.dist }{
charater value specifying a path to a SOLAR-formatted .dist file that
maps the character marker names to numeric centimorgan values.
}
\item{output.directory }{
character value specifying a path to a directory that the output files
(mloci.out and share.out) will be placed.  If any of the specified
directory path does not exist,
\code{solar2multic} will create the necessary 
directories.
}
\item{delete.fixed.dir }{
logical flag: if \code{TRUE} (default), then the
temporary directory that is created to hold intermediate files is
deleted. 
}
}
\section{Side Effects}{
Due to write permissions possibly not allowing the user to gunzip
and create files in the specified
\code{directory}, 
\code{solar2multic} first copies 
\code{directory} and
\code{phi2} to the current directory.
\code{solar2multic} then 
creates a temporary directory to hold the "fixed," intermediate files
that will be deleted (by default).  Also,
\code{solar2multic} will overwrite 
mloci.out, mloci.out.gz, share.out and share.out.gz if they exist in
\code{output.directory}.
}
\seealso{
\code{\link{solar2mloci}},
\code{\link{phi2share}}
}
\examples{
\dontrun{
solar2multic(phi2 = "phi2.gz",
             pedigree.file = "chrom18.ped",
             pedindex.out = "pedindex.out",
             pedindex.cde = "pedindex.cde",
             ibd.directory = "mibddir",
             output.directory = "multicInput",
             delete.fixed.dir = FALSE)

solar2multic("solarOutput/phi2",
             "solarOutput/chrom10.ped",
             "solarOutput/pedindex.out",
             "solarOutput/pedindex.cde",
             "mibds/chrm10"))
}
}
% docclass is function
% Converted by Sd2Rd version 52193.
