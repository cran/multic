#############################################################################
## sw2mloci is a function to create a multic input file mloci.out from 
## SimWalk IBD files.  It takes as arguments the name (path) of the directory
## that holds the IBD files and an optional .map file argument in map.
## It creates a file named mloci.out.gz in the current directory.  If
## mloci.out.gz or mloci.out are already present, sw2mloci will move them to
## mloci.out.before.gz or mloci.out.before respectively.
#############################################################################
sw2mloci <- function(ibd.directory, map="", output.directory=".",
                     famid=NULL, id=NULL, dadid=NULL, momid=NULL,
                     directory=NULL) {

  ## This function requires R or Splus 8
  if (!is.R() && version$major < 8) {
    stop("sw2mloci function requires R or Splus 8 or greater.\n")
  }

  ## Check for valid famid, id, dadid, and momid (needed for kinship file.)
  if (is.null(famid) || is.null(id) || is.null(dadid) || is.null(momid)) {
    stop("\nInsufficient pedigree information passed.\n",
         "Make sure that famid, id, dadid, and momid arguments",
         "are specified.\nsw2mloci.q key ??")
  }

  ## "directory" was renamed to "ibd.directory to be consistent with
  ## solar2multic.  But, we want to stay backward compatible.
  if (is.null(ibd.directory) && !is.null(directory)) {
    ibd.directory <- directory
  }
  
  ibd.directory <- as.character(ibd.directory)
  output.directory <- as.character(output.directory)

  ## Check to see if all the input file names represent existing files.
  if(!multic.is.dir(ibd.directory)) {
    stop("\nDirectory ", ibd.directory, " does not exist.\n",
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
  all.files <- list.files(ibd.directory)

  ## Narrow the files down to only the IBD files.
  ibd.regexpr <- "^IBD-.*"
  ibd.files <- all.files[regexpr(ibd.regexpr, all.files) > 0]

  ## If there are no IBD files, exit.
  if(length(ibd.files) == 0) {
    stop("\nThere are no IBD-* files in ", ibd.directory, "\nsw2mloci.q key 25")
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
  directory.orig <- ibd.directory
  ibd.directory <- tempfile("sw2mloci.")
  copy.files.to.tmp.space(ibd.directory,
                          paste(directory.orig, ibd.files, sep = "/",
                                collapse = " "), 
                          map)

  ## gunzip the files as necessary (keeping the names of the new files).
  gziped.files <- ibd.files[regexpr("\\.gz$", ibd.files) > 0]
  gunziped.files <- NULL
  if(length(gziped.files) > 0) {
    cat("gunzip'ing IBD files ...\n")
    gunziped.files <- sapply(paste(ibd.directory, gziped.files, sep = "/"),
                             gunzip)
  }

  ## Refresh the list after files have been gunzip'ed.
  all.files <- list.files(ibd.directory)
 
  ibd.files <- paste(ibd.directory,
                     all.files[regexpr(ibd.regexpr, all.files) > 0],
                     sep = "/")

  ## Call the C routines to create the temporary mibd.* files and create
  ## mloci.
  .C("sw2lociFiles", as.character(ibd.files),
     as.integer(length(ibd.files)),
     as.character(map),
     PACKAGE = "multic")

  ## The mloci.out created by the C++ function is incorrect.
  ## We need to fix it by getting the id's in the proper columns
  ## and merging with a kinship file to make sure that all combinations
  ## are included.
  mergeWithKinship(famid=famid, id=id, dadid=dadid, momid=momid)
  
  ## Compress (gzip) mloci.out.
  cat("gzip'ing mloci.out\n")
  output.file <- gzip("mloci.out")

  ## Check to see if all the output.directory exists.
  if(!multic.is.dir(output.directory)) {
    result <- dir.create(output.directory)
    if( !result ) {
      warning(output.directory, " could ne be created or written to.\n",
              output.file, " will be in the current directory.\n",
              "swm2mloci.q key 97")
      output.directory = "."
    }
  }
  if(output.directory != ".") {
    multic.system(paste("mv", output.file, output.directory))
    output.file <- paste(output.directory, output.file, sep = "/")
  }

  ## Remove the temporary directory.
  multic.system(paste("rm -r", ibd.directory))
  
  return(output.file)
}

###
# The IBD files likely don't contain all pedigree combinations.
# Create a kinship file to generate all id pairs and use the phi2
# value there if it's not in the IBD files.
###
mergeWithKinship <- function (famid, id, dadid, momid) {
  ## args were checked for NULL in parent function
  
  ## Make a kinship file and read it back in.
  ## Note:  "share.out" is named kinship if we don't have mloci
  kinTry <- try(make.kinship.file(famid, id, dadid, momid, "kinship", FALSE, NULL))
  if (inherits(kinTry, "try-error")) {
    traceback()
    stop("\nmakeBlockIdPair failed because insufficient memory could be allocated.\n",
         "Try running sw2mloci on a computer with more memory.  Or if you are\n",
         "using R and have access to Splus, try using the Splus version of multic.\n")
  }
  kinship.matrix <- read.table("kinship", as.is=TRUE)[, 1:3]
  kinship.matrix <- expandMibd(kinship.matrix)
  dimnames(kinship.matrix)[[2]] <- c("famid", "id1", "id2", "phi2.kin")

  ######
  # mloci.out has all the mibd file information concatenated into one file.
  # It requires some processing to separate them into one table for each mibd.
  # Each new mibd starts with a comment line (e.g. "# mibd.22.21.406")
  ######

  ## Read all the lines of mloci.out out as vector of strings
  ## and delete existing file (will write out corrected version).
  mloci.out.file <- readLines("mloci.out", n = -1)
  multic.system("rm -f mloci.out") # replace system call when Splus 8 is released?

  ## Get rid of the blank lines
  mloci.out.file <- mloci.out.file[mloci.out.file != ""]

  # Get the indeces of the mibd names and add marker for the end of file
  mibd.index <- grep("#", mloci.out.file)
  mibd.index <- c(mibd.index, (length(mloci.out.file) + 1))

  cat("Merging the mibd's with kinship file.\n",
      "Please be patient.  This could take a while.")
  
  ## Example:  First mibd records are from 2 through 86
  # mibd.index
  # [1]   1  87 173 259 345 431 517 603 689 775 861

  for (ii in 1:(length(mibd.index) - 1)) {
    
    start.index <- mibd.index[ii] + 1
    stop.index <- mibd.index[ii + 1] - 1

    #cat("\nProcessing mibd", ii, "start index =", start.index, "stop index =", stop.index, "\n")
    #cat(mloci.out.file[mibd.index[ii]])
    
    this.mibd <- mloci.out.file[start.index:stop.index]

    #browser()
    
    ## mibd is all in text, so we have to convert it to a matrix
    this.mibd <- strsplit(this.mibd, "\t")
    this.mibd <- do.call("rbind", this.mibd)
    this.mibd <- data.frame(I(this.mibd[, 1]), I(this.mibd[, 2]),
                            as.numeric(this.mibd[, 3]))

    ## Assign the dimnames (and rewrite C++ function so that last col is not created)
    dimnames(this.mibd)[[2]] <- c("id1", "id2", "phi2")

    ## Multiple phi2 by 2 so that it really is "2 phi"
    this.mibd$phi2 <- 2 * this.mibd$phi2

    ## Expand famid to its own column and make sure id's are in correct column
    this.mibd <- expandMibd(this.mibd)
    this.mibd <- sortIbdCols(this.mibd) # puts ids in correct column

    ## Now, we want to add the phi2 values to the kinship.file data.frame.
    this.mibd <- merge(kinship.matrix, this.mibd, all.x=TRUE)

    ## New entries in phi2 column are NA.  Need to pull from kinship.
    na.idx <- is.na(this.mibd[, "phi2"])
    this.mibd[na.idx, "phi2"] <- this.mibd[na.idx, "phi2.kin"]

    ## Drop the phi2.kin column
    col.idx <- dimnames(this.mibd)[[2]] != "phi2.kin"
    this.mibd <- this.mibd[, col.idx]

    # Doesn't come out sorted in Splus
    idx <- order(this.mibd[, 1], this.mibd[, 2], this.mibd[, 3])
    this.mibd <- this.mibd[idx, ]
    
    ## Collapse famid back into the id columns
    this.mibd <- collapseMibd(this.mibd)

    ## Write out mibd header info and mibd table
    write(mloci.out.file[mibd.index[ii]], file="mloci.out", append=TRUE)
    write.table(this.mibd, file="mloci.out", quote=FALSE,
                row.names=FALSE, col.names=FALSE, sep="\t", append=TRUE)

    cat(".")
    if (ii %% 50 ==0) {
      cat("\n")
    }
  } ## end for all mibd's in mloci.out file
  cat("\n")
}

#####
# Function to unpaste famid from person id info and put
# it into a new column and returns matrix.
#####
expandMibd <- function (mibd) {
  famid.and.id1.col <- do.call("rbind", strsplit(mibd[, 1], "-"))
  id2.col <- do.call("rbind", strsplit(mibd[, 2], "-"))[, 2]

  phi2 <- NULL
  if (ncol(mibd) > 2) {
    phi2 <- mibd[, 3]
  }
  
  new.dimnames <- c("famid", dimnames(mibd)[[2]])

  # Make sure it handles mibd with only 2 cols
  numCol <- ncol(mibd)
  mibd <- cbind(as.numeric(famid.and.id1.col[, 1]),
                as.numeric(famid.and.id1.col[, 2]),
                as.numeric(id2.col)
                )
  # We may or may not have a phi2 column
  if (!is.null(phi2)) {
    mibd <- cbind(mibd, as.numeric(phi2))
  }
  dimnames(mibd)[[2]] <- new.dimnames
  return(mibd)
}

collapseMibd <- function (mibd) {
  if (dimnames(mibd)[[2]][1] != "famid" && dimnames(mibd)[[2]][1] != "FAMILY") {
    stop("ERROR: collapseMibd can't find famid column to collapse.\n")
  }
  # remove "famid" from dimnames
  new.dimnames <- dimnames(mibd)[[2]][-1]

  # Now, we need to paste the famid to the next to columns separated by "-"
  id1 <- paste(mibd[, 1], mibd[, 2], sep="-")
  id2 <- paste(mibd[, 1], mibd[, 3], sep="-")

  # Do we want a matrix or data frame?  matrix seems to work fine
  mibd <- data.frame(I(id1), I(id2), as.numeric(mibd[, 4]))
  #mibd <- cbind(id1, id2, as.numeric(mibd[, 4:ncol(mibd)]))

  dimnames(mibd)[[2]] <- new.dimnames
  
  return(mibd)
}

# Input matrix must be numeric
sortIbdCols <- function (mibd) {
  need.to.swap <- mibd[, 2] > mibd[, 3]
  for (i in 1:nrow(mibd)) {
    if (need.to.swap[i]) {
      # swap columns 2 and 3
      tmp <- mibd[i, 2]
      mibd[i, 2] <- mibd[i, 3]
      mibd[i, 3] <- tmp
    }
  }
  return(mibd)
}

copy.files.to.tmp.space <- function(tmp.dir, ibd.files, map.file) {
  copy.command <- paste("mkdir", tmp.dir, ";",
                        "cp -r", ibd.files, map.file, tmp.dir)
  multic.system(copy.command)
}
