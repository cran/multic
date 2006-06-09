#############################################################################
## mlociCut is a function to make a smaller version of an mloci.out file by
## keeping only centimorgan values within a specified "region" of the
## original mloci.out.

## Eric Lunde, 2006-04-03
## $Id: mlociCut.q,v 1.2 2006/04/03 16:05:06 lunde Exp $
#############################################################################
mlociCut <- function(mloci.out = 'mloci.out', region, mloci.cut='mloci.cut',
                     output = TRUE)
  ## Possible signature
  ## mlociCut <- function(mloci.out = 'mloci.out', from, to,
  ##                      mloci.cut='mloci.cut', output = T,
  ##                      by.index = FALSE)
  ## from, to, by, length, cM=T?F
{  
  ## this fixes it for now, but is over a region
  from<-min(region)
  to<-max(region)
  
  ## gunzip and rm mloci and outf respectively, if necessary
  if(output) {
    cat(paste("gunzip'ing", mloci.out, "\n"))
  }
  gunzip(mloci.out)
  multic.system(paste("rm -f", mloci.cut))

  ## read mloci and determine line number of "#"'s and the loci length
  if(output) {
    cat(paste("Cutting ", mloci.out, "...\n", sep = ""))
  }
  lines <- multic.system(paste("cat", mloci.out))
  pound.lines <- grep("^#", lines)
  if(length(pound.lines) > 1) {
    loci.length <- pound.lines[2] - 1
  }else {
    loci.length <- length(pound.lines)
  }

  ## for each locus, determine if it is in the range desired, if so, print
  ## the locus
  for(pound.line in pound.lines) {
    line <- lines[pound.line]
    distance = as.numeric(paste(multic.strsplit(line, sep = ".")[-(1:2)],
                                collapse = '.'))
    if(from <= distance & distance <= to) {
      cat(lines[seq(pound.line, length = loci.length)], file = mloci.cut,
          append = TRUE, sep = "\n")
    }
  }

  if(output) {
    cat(paste("gzip'ing", mloci.out, "and", mloci.cut, "\n"))
  }
  gzip(mloci.out)
  gzip(mloci.cut)
  invisible()
}
