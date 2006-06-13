#############################################################################
## Title: expand.multic
## Description: expand.multic is a utility function to create "bootstrap"ed
##              versions of mloci.out and share.out.
## Input: mloci.out - character variable specifying the path to an mloci.out
##                    (or similarly formatted) file.
##        share.out - character variable specifying the path to an share.out
##                    (or similarly formatted) file.
## Output: files - list of two character values specifying the names of the
##                 "expand"ed mloci.out and share.out
## Side Effects: "expand"ed versions of the input mloci.out and share.out
##               files are created in the current directory.
##
##               A directory named "loci" is used as temporary storage
##               for the split mloci.out file.
## Author: Eric Lunde 2005-08-16
#############################################################################
expand.multic <- function(famids, mloci.out = NULL, share.out = NULL) {
 
  new.mloci.out <- expand.mloci(famids, mloci.out)
  new.share.out <- expand.share(famids, share.out)
  
  files <- list(new.mloci.out = new.mloci.out, new.share.out = new.share.out)
  
  return (files)
}

expand.mloci <- function(famids, file) {
  if(is.null(file)) {
    return (NULL)
  }

  ## Extract only the base name of the file (remove any directory info) and
  ## remove the .gz if it exists
  base.file <-
    get.basename(file, keep.extension = "gz" != get.extension(file))

  new.file <- paste(base.file, "expanded", sep = ".")
  multic.system(paste('rm -f', new.file))

  ## Copy and gunzip (if necessary)
  file <- gunzip(file, copy = TRUE)
  
  loci.names <- mloci.split(file)
  
  sapply(loci.names, expand.locus, famids, new.file)
    
  ## delete loci dir and copied file
  multic.system('rm -rf loci')
  multic.system(paste('rm -f', file))

  return (new.file)
}

expand.share <- function(famids, file) {
  if(is.null(file)) {
    return (NULL)
  }

  ## Extract only the base name of the file (remove any directory info) and
  ## remove the .gz if it exists
  base.file <-
    get.basename(file, keep.extension = "gz" != get.extension(file))
  
  new.file <- paste(base.file, "expanded", sep = ".")
  multic.system(paste('rm -f', new.file))

  ## Copy and gunzip (if necessary)
  file <- gunzip(file, copy = TRUE)
  
  expand.locus(file, famids, new.file)
  
  ## delete the copied file
  multic.system(paste('rm -f', file))

  return (new.file)
}

expand.locus <- function(locus.name, famids, new.file) {
  lines <- readLines(locus.name, n = -1)
  ## if first line has # repeat it.
  ## This was added because the majority of the functionality for "expand"ing
  ## a locus is the same for a share.out, except for the "#" header.
  if(substring(lines[1], 1, 1) == "#") {
    cat(lines[1], "\n", file = new.file, append = TRUE)
  }  
  
  expand.famid <- function(famid, lines, new.file) {
    ## grep for famids
    ## This could be sped-up a little by using substring to get the
    ## beginning of each string.

    ## Should I do some error checking to see if the famid is even a valid
    ## one (one that exists in the file)?
    famid.lines <- lines[grep(paste("^", famid, "-", sep = ""), lines)]
    if(length(famid.lines) > 0) {
      cat(famid.lines, sep = "\n", file = new.file, append = TRUE)
    }
    invisible()
  }

  ## Make sure all famids supplied are valid.  If not, produce error message
  regexpr.match <- regexpr("^[0-9a-zA-Z]+", lines)
  unique.famids <- 
    unique(substring(lines, regexpr.match,
                     attr(regexpr.match, "match.length")))
  matches <- match(famids, unique.famids)
  if(any(is.na(matches))) {
    ## Ask for Beth's help about this topic.  Can I have this printed to the
    ## warnings only once. or without "FUN(...X.sub.i...., ..."?
    warning("The famids ", famids[is.na(matches)], " were not found")
  }
  
  sapply(famids, expand.famid, lines, new.file)
  invisible()
}

expand.data <- function(famids, d.frame) {
  ## Make sure all famids supplied are valid.  If not, produce error message
  matches <- match(famids, unique(d.frame$famid))
  if(any(is.na(matches))) {
    warning("The famids ", paste(famids[is.na(matches)], collapse = " "),
            " were not found")
  }

  famid.data <- function(famid, d.frame) {
    return (d.frame[d.frame$famid == famid, ])
  }

  ## I'd like to find a faster solution to this problem, but this
  ## sledgehammer gets the job done.  About 45 seconds for 100 families.
  result <- NULL
  for(famid in famids) {
    result <- rbind(result, famid.data(famid, d.frame))
  }
  return (result)
}

get.basename <- function(file, keep.extension = TRUE) {  
  ## Extract only the base name of the file (remove any directory info)
  split.file <- multic.strsplit(file, sep = "/")
  basename <- split.file[length(split.file)]
  if( !keep.extension ) {
    basename.split <- multic.strsplit(basename, sep = ".")
    basename <- paste(basename.split[-length(basename.split)],
                      collapse = ".")
  }
  return (basename)
}
