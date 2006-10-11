#############################################################################
## sw2mloci is a function to create a multic input file mloci.out from 
## SimWalk IBD files.  It takes as arguments the name (path) of the directory
## that holds the IBD files and an optional .amp file argument in map.
## It creates a file named mloci.out.gz in the current directory.  If
## mloci.out.gz or mloci.out are already present, sw2mloci will move them to
## mloci.out.before.gz or mloci.out.before respectively.
#############################################################################
sw2mloci <- function(directory, map = "", output.directory = ".") {

  stop("The sw2mloci function does not currently work properly.  ",
       "Look for bug fixes in a future multic release.\n")
  
  directory <- as.character(directory)
  output.directory <- as.character(output.directory)

  ## Check to see if all the input file names represent existing files.
  if(!multic.is.dir(directory)) {
    stop("\nDirectory ", directory, " does not exist.\n",
         "sw2mloci.q key 6")
  }

  if( map != "" ) {
    map <- as.character(map)
    if( !file.exists(map) ) {
      stop(paste("\nThe map file ", map, " does not exist.\n",
                 "sw2mloci.q key 13", sep = ""))
    }
  }

  ## Get the names of all the files in the specified directory.
  all.files <- if(using.R()) {
    list.files(directory)
  } else {
    files.in.dir(directory)
  }

  ## Narrow the files down to only the IBD files.
  ibd.regexpr <- "^IBD-.*"
  ibd.files <- all.files[regexpr(ibd.regexpr, all.files) > 0]

  ## If there are no IBD files, exit.
  if(length(ibd.files) == 0) {
    stop("\nThere are no IBD-* files in ", directory, "\nsw2mloci.q key 25")
  }

  ## Move sw2mloci's output files, if they exist, so it doesn't overwrite
  ## them.
  if(file.exists("mloci.out")) {
    cat("mloci.out exists, moving to mloci.before ...\n")
    multic.system("mv -f mloci.out mloci.out.before")
  }
  if(file.exists("mloci.out.gz")) {
    cat("mloci.out.gz exists, moving to mloci.out.before.gz ...\n")
    multic.system("mv -f mloci.out.gz mloci.out.before.gz")
  }

  ## We only need to move the target mloci in the output dir.
  target.mloci <- paste(output.directory, "mloci.out", sep = "/")
  if(file.exists(paste(target.mloci, "gz", sep = "."))) {
    cat(target.mloci, ".gz exists, moving to ", target.mloci,
        ".before.gz ...\n", sep = "")
    multic.system(paste("mv -f ", target.mloci, ".gz ", target.mloci,
                        ".before.gz", sep = ""))
  }

  ## Keep the original directory name, but begin using a tempoaray directory.
  ## Copy the IBD and .map files to the temp space.
  directory.orig <- directory
  directory <- tempfile("sw2mloci.")
  copy.files.to.tmp.space(directory,
                          paste(directory.orig, ibd.files, sep = "/",
                                collapse = " "), 
                          map)

  ## gunzip the files as necessary (keeping the names of the new files).
  gziped.files <- ibd.files[regexpr("\\.gz$", ibd.files) > 0]
  gunziped.files <- NULL
  if(length(gziped.files) > 0) {
    cat("gunzip'ing IBD files ...\n")
    gunziped.files <- sapply(paste(directory, gziped.files, sep = "/"),
                             gunzip)
  }

  ## Refresh the list after files have been gunzip'ed.
  all.files <- if(using.R()) {
    list.files(directory)
  } else {
    files.in.dir(directory)
  }
  ibd.files <- paste(directory,
                     all.files[regexpr(ibd.regexpr, all.files) > 0],
                     sep = "/")

  ## Call the C routines to create the temporary mibd.* files and create
  ## mloci.
  .C("sw2lociFiles", as.character(ibd.files),
     as.integer(length(ibd.files)),
     as.character(map),
     PACKAGE = "multic")

  ## Compress (gzip) mloci.out.
  cat("gzip'ing mloci.out\n")
  output.file <- gzip("mloci.out")

  ## Check to see if all the output.directory exists.
  if(!multic.is.dir(output.directory)) {
    result <- multic.mkdir(output.directory)
    if( !result ) {
      warning(output.directory, " could ne be created or writen to.\n",
              output.file, " will be in the current directory.\n",
              "swm2loci.q key 97")
      output.directory = "."
    }
  }
  if(output.directory != ".") {
    multic.system(paste("mv", output.file, output.directory))
    output.file <- paste(output.directory, output.file, sep = "/")
  }

  ## Remove the temporary directory.
  multic.system(paste("rm -r", directory))
  
  return(output.file)
}

copy.files.to.tmp.space <- function(tmp.dir, ibd.files, map.file) {
  copy.command <- paste("mkdir", tmp.dir, ";",
                        "cp -r", ibd.files, map.file, tmp.dir)
  multic.system(copy.command)
}
