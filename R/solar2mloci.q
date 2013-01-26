solar2mloci <- function(directory, phi2, pedindex.out, pedindex.cde,
                      ibd.dist = NULL, output.directory = ".",
                      delete.fixed.dir = TRUE) {
  ## Test directory for a supplied argument, existence, and (m)ibd contents
  if(missing(directory)) {
    stop("\nargument \"directory\" is missing with no default\n",
         "solar2mloci.q key 17")
  }
  if(!file.exists(directory)) {
    stop("\ndirectory = \"", directory, "\" does not exist.\n",
         "solar2mloci.q key 12")
  }
  file.list <- multic.system(paste("ls -1 ", directory, "/*", sep = ""))
  all.ibd.files <- file.list[grep(paste(directory, "m?ibd\\..*", sep = "/"),
                                  file.list)]
  if(length(all.ibd.files) == 0) {
    stop("\ndirectory = \"", directory,
         "\" contains no ibd or mibd files.\n",
         "solar2mloci.q key 15")
  }


  ## Test phi2 for a supplied argument
  if(missing(phi2)) {
    stop("\nargument \"phi2\" is missing with no default\n",
         "solar2mloci.q key 17")
  }

  ## Test pedindex.out for a supplied argument and existence
  if(missing(pedindex.out)) {
    stop("\nargument \"pedindex.out\" is missing with no default\n",
         "solar2mloci.q key 17")
  }
  if(!file.exists(pedindex.out)) {
    stop("\npedindex.out = \"", pedindex.out, "\" does not exist.\n",
         "solar2mloci.q key 21")
  }

  ## Test pedindex.cde for a supplied argument and existence
  if(missing(pedindex.cde)) {
    stop("\nargument \"pedindex.cde\" is missing with no default\n",
         "solar2mloci.q key 26")
  }
  if(!file.exists(pedindex.cde)) {
    stop("\npedindex.cde = \"", pedindex.cde, "\" does not exist.\n",
         "solar2mloci.q key 30")
  }

  ## if provided, test ibd.dist for existence
  if(!is.null(ibd.dist)) {
    if(!file.exists(ibd.dist)) {
      stop("\nibd.dist = \"", ibd.dist, "\" does not exist.\n",
           "solar2mloci.q key 38")
    }
  }

  ## Test phi2 for existence (and copy and gunzip).  This is performed after
  ## all the other tests because if this works and some other file is missing
  ## the "local" copy of phi2 still remains.  I don't want to create it until
  ## I'm sure it is needed. - Eric Lunde 2005-10-13
  phi2 <- gunzip(phi2, copy = T)

  ## Copy directory if it is not the current directory.
  ## Think about if the dir is a parent dir

  ## Copy directory
  if("." != directory) {
    local.dir.name <- paste("local", get.basename(directory), sep = "")

    ## If local* already exists, stop and tell the user what is going on.
    if(multic.is.dir(local.dir.name)) {
      remove.file(phi2)
      stop("\nThe directory \"", local.dir.name, "\" exists.  solar2mloci ",
           "needs this\ndirectory for temporary storage. Please rename it ",
           "if it contains data you\nwish to keep, or delete it.\n",
           "solar2mloci.q key 73")
    }
    multic.system(paste("mkdir", local.dir.name))
    multic.system(paste("cp", paste(all.ibd.files, collapse = " "),
                        local.dir.name))

    ## Use on.exit-type function to delete dir.
    on.exit(multic.system(paste("rm -r", local.dir.name)))
  } else {
    local.dir.name <- "."
  }

  ## Save original directory
  directory.orig <- directory

  ## The rest of the program relies on directory, so save local.dir.name
  ## on to it.
  directory <- local.dir.name
  
  ## Clear the space for new mloci.out or mloci.out.gz
  remove.file(paste(output.directory, "/mloci.out", sep = ""))
  remove.file(paste(output.directory, "/mloci.out.gz", sep = ""))
  
  fixed.ibd.directory <- paste(directory, "_fixed", sep = "")
  fixibd(directory, fixed.ibd.directory, phi2, directory.orig)
  
  ## Get the names of all the (m)ibd files.  
  all.file.names <- multic.system(paste('ls', fixed.ibd.directory))
  
  ## for file names that begin with "mibd."
  ## make a data frame of "mibd", "15", "10" and convert the numeric ones to
  ## numeric and reorder, send that reordered list as ibd.names
  mibd.names <- all.file.names[substring(all.file.names, 1, 5) == "mibd."]
  ibd.names <- all.file.names[substring(all.file.names, 1, 4) == "ibd."]
  
  if(length(mibd.names) != 0) {
    ## In general, mibds will be of the format "mibd.*.*".  Sometimes, they
    ## will be in the format "mibd.*.*.*".  And it even possible that both
    ## formats will be used.  In this last case, we need to determine which
    ## files are in the first format and convert them to the second format
    ## for consistancy.  needs.dot.zero determines if the conversions are
    ## necessary and normalize.mibds does the actual conversion.
    needs.dot.zero <-
      unlist(lapply(sapply(mibd.names,
                           multic.strsplit, sep = '.'), length)) == 3
    
    if(any(needs.dot.zero) && !all(needs.dot.zero)) {
      mibd.names <- sapply(mibd.names, normalize.mibd, fixed.ibd.directory)
    }
    
    mibd.mat <- matrix(multic.strsplit(mibd.names, sep = "."),
                       byrow = T, nrow = length(mibd.names))
    if(ncol(mibd.mat) == 4) {
      mibd.mat[,3] <- paste(mibd.mat[,3], mibd.mat[,4], sep = '.')
      mibd.mat <- mibd.mat[, -4, drop = FALSE]
    }
    
    ## Order the centimorgan measure
    numeric.values <- as.numeric(mibd.mat[,3])
    numeric.values <- numeric.values[order(numeric.values)]
    temp.order <- match(numeric.values, as.numeric(mibd.mat[,3]))
    mibd.mat <- mibd.mat[temp.order, , drop = FALSE]
    
    ## Order the chromosome numbers
    ##mibd.mat[,2] <- as.numeric(mibd.mat[,2])
    ##mibd.mat <- mibd.mat[order(mibd.mat[,2]), ]
    mibd.mat <- mibd.mat[order(as.numeric(mibd.mat[,2])), , drop = FALSE]
    
    mibd.mat[, 2] <- trim(mibd.mat[, 2])
    mibd.mat[, 3] <- trim(mibd.mat[, 3])
    
    mibd.names <- apply(mibd.mat, 1, paste, collapse = '.')
  }
  
  if( !is.null(ibd.dist) && ibd.dist != "" && length(ibd.names) != 0) {
    ibd.names.orig <- ibd.names
    ibd.names <-
      data.frame(matrix(multic.strsplit(ibd.names, sep = "."),
                        ncol = 2, byrow = TRUE))
    names(ibd.names) <- c("ibd", "marker")

    chrm.number <- scan(ibd.dist, n = 1, quiet = TRUE)

    marker.map <- read.table(ibd.dist, skip = 1, row.names = NULL)
    names(marker.map) <- c("marker", "cM")

    match.order <- match(ibd.names$marker, marker.map$marker)

    non.matched.ibds <- 
      if(nrow(ibd.names[is.na(match.order), ]) > 0) {
        apply(ibd.names[is.na(match.order), ], 1, paste,
              collapse = ".")
      } else {
        NULL
      }
    matched.ibds <- 
      if(nrow(ibd.names[na.omit(match.order), ]) > 0) {
        apply(ibd.names[na.omit(match.order), ], 1, paste,
              collapse = ".")
      } else {
        NULL
      }

    if(any(is.na(match.order))) {
      warning("The ibd files:\n", paste(non.matched.ibds, collapse = ", "),
              "\ndid not find a match in ", ibd.dist, "\n",
              "solar2mloci.q key 168")
    }
    marker.map <- marker.map[match.order, ]

    ## Replace any matched ibd marker with its cM value
    new.ibd.names <- paste(ibd.names$ibd, ibd.names$marker, chrm.number,
                           marker.map$cM, sep = ".")

    ibd.names <- new.ibd.names
    if(any(is.na(match.order))) {
      ibd.names[is.na(match.order)] <- non.matched.ibds
    }

    ## Rename the appropriate files
    if(any(!is.na(match.order))) {
      multic.system(paste("mv ", fixed.ibd.directory, "/",
                          ibd.names.orig[!is.na(match.order)],
                          " ", fixed.ibd.directory, "/",
                          ibd.names[!is.na(match.order)],
                          sep = "", collapse = ";"))
    }
  }
  
  ibd.names <- c(ibd.names, mibd.names)
  
  ## Get the number of (m)ibd files.
  ibd.count <- length(ibd.names)

  ## If the output directory is not the current directory and the directory
  ## does not exist, create it.
  if( output.directory != '.' ) {
    if( !file.exists(output.directory) ) {
      multic.system(paste('mkdir -p', output.directory))
    }
  }

  ids <- read.pedindex.out(pedindex.out, pedindex.cde)

  .C("solar2mloci",
     directory = as.character(fixed.ibd.directory),
     ibd.names = as.character(ibd.names),
     ibd.count = as.integer(ibd.count),
     uniqueIds = as.character(ids[,4]),
     gzipWhenComplete = as.integer(!delete.fixed.dir),
     PACKAGE = "multic"
     )

  gzip("mloci.out")

  ## If the output directory is not the current directory, move mloci.out
  ## to that output directory.
  if( output.directory != '.' ) {
    multic.system(paste('mv -f mloci.out.gz ', output.directory, sep = ""))
  }

  ## Delete the 'fixed' directory
  if(delete.fixed.dir) {
    multic.system(paste('rm -rf', fixed.ibd.directory))
  }

  remove.file(phi2)

  invisible()
}

normalize.mibd <- function(mibd.name, directory) {
  split.name <- multic.strsplit(mibd.name, sep='.')
  if(length(split.name) == 3) {
    new.name <- paste(mibd.name, ".0", sep = "")
    multic.system(paste("mv", paste(directory, mibd.name, sep = "/"),
                        paste(directory, new.name, sep = "/")))
    return (new.name)
  }
  return(mibd.name)
}
