## This file brings together the functions that were different enough between
## R and Splus that they needed to be encapsulated in a universal function.

## This file also brings together the general purpose and Unix-like commands

remove.file <- function(file.name) {
  if(file.exists(file.name)) {
    multic.system(paste('rm', file.name))
  }

  invisible()
}

clean <- function() {
  remove.temp.files()
  multic.system("rm -f local*")
  invisible()
}

remove.temp.files <- function() {
  ## If these files exist, delete them
  singular.files <- c("fam.lik",
                      "fam.lik0",
                      "fam.likA",
                      "fort.12",
                      "invExpSecDerFixed.log",
                      "invExpSecDerRandom.log",
                      "initial.values",
                      "null.initial.values",
                      "alt.initial.values",
                      "least.log",
                      "multic.env",
                      "multic.lik",
                      "multic.log",
                      "multic.mg1",
                      "multic.mg2",
                      "multic.mu",
                      "multic.out",
                      "multic.po",
                      "multic.poly",
                      "multic.pp",
                      "multic.sib",
                      "multiclookup.log",
                      "summary.log",
                      "varCovar.log",
                      "tempmloci.out",
                      "loci.out",
                      "iterations.log")
  
  multic.system(paste(c('rm -f', singular.files, 'v.matrix.*.log',
               'y.beta.diff.*.log'), collapse = ' '))  
  multic.system('rm -rf loci')
  
  invisible()
}

gunzip <- function(file, copy = FALSE) {
  if(file.exists(file)) {
    if(copy) {
      basename <- get.basename(file)
      
      local.file <- paste("local", basename, sep = "")
      multic.system(paste("cp -f", file, local.file))
      file <- local.file
    }
  }else {
    first.try <- file
    if("gz" == get.extension(file)) {
      file <- substring(file, 1, nchar(file) - 3)
    } else {
      file <- paste(file, ".gz", sep = "")
    }
    
    if(file.exists(file)) {
      if(copy) {
        basename <- get.basename(file)
        
        local.file <- paste("local", basename, sep = "")
        multic.system(paste("cp -f", file, local.file))
        file <- local.file
      }
    } else {
      stop(paste("\nNeither '", first.try, "' nor '", file, "' exist.\n",
                 "multic.q key 263", sep = ""))
    }
  }
  
  if("gz" == get.extension(file)) {
    file <- substring(file, 1, nchar(file) - 3)
  }
  multic.system(paste("gunzip -qf", file))
  
  return (file)
}

gzip <- function(file) {
  if(!file.exists(file)) {
    warning(file, ": No such file or directory.\n", "multic.q key 1617")
    return(invisible(file))
  }

  if("gz" == get.extension(file)) {
    warning(file, " already has a .gz suffix -- unchanged.\n",
            "multic.q key 1623")
    return(invisible(file))
  }

  multic.system(paste("gzip -f", file))
  gzip.file <- paste(file, "gz", sep = ".")

  invisible(gzip.file)
}

## Given a file name, return the characters following the last . (dot)
get.extension <- function(file) {
  ## This was first implemented without the get.basename step, but that
  ## caused issues with a name like `projects/multic.test/file'
  basename <- get.basename(file)
  split.file <- multic.strsplit(basename, sep = ".")
  extension <- split.file[length(split.file)]
  return (extension)
}

trim <- function(str) {
  result <- sub("[ \t\n]*$", "", str)
  result <- sub("^[ \t\n]*", "", result)

  return (result)
}

recover.multic <- function() {
  print("please implement recover.multic")
}

## Shouldn't need as of Splus 8, but we'll keep it around for now.
multic.strsplit <- function(str, sep = " ") {
  result <- NULL
  if(!missing(sep)){
    result <- unlist(strsplit(str, sep, fixed = TRUE))
  }

  result <- result[result != ""]
  
  return (result)
}

multic.strings.split <- function(strings, sep = " ")  {
  split.strings <- t(sapply(strings, multic.strsplit, sep))
  return (split.strings)
}

multic.system <- function(command) {
  result <- system(command, intern = TRUE)
  return (result)
}

multic.rownames <- function(x, do.NULL = TRUE, prefix = "row") {
  result <- NULL
  result <- rownames(x, do.NULL, prefix)

  return (result)
}

multic.colnames <- function(x, do.NULL = TRUE, prefix = "col") {
  result <- colnames(x, do.NULL, prefix)
  return (result)
}

multic.is.dir <- function(dir) {
  result <- file.info(dir)[,2]
  if(is.na(result)) {
    result <- FALSE
  }
  return (result)
}

multic.kurtosis <- function (x, na.rm = FALSE) {
  if (na.rm) 
    x <- x[!is.na(x)]
  result <- sum((x - mean(x))^4)/(length(x) * var(x)^2) - 3
  return (result)
}

multic.skewness <- function (x, na.rm = FALSE) {  
  std.dev <- sd

  if (na.rm) 
    x <- x[!is.na(x)]
  result <- sum((x - mean(x))^3)/(length(x) * std.dev(x)^3)
  return (result)
}
