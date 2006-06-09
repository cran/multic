phi2share <- function(phi2, pedigree.file, pedindex.out, pedindex.cde,
                      output.directory = '.') {
  ## Test phi2 for a supplied argument
  if(missing(phi2)) {
    stop("\nargument \"phi2\" is missing with no default\n",
         "phi2share.q key 6")
  }

  ## Test pedigree.file for a supplied argument and existence
  if(missing(pedigree.file)) {
    stop("\nargument \"pedigree.file\" is missing with no default\n",
         "phi2share.q key 14")
  }
  if(!file.exists(pedigree.file)) {
    stop("\npedigree.file = \"", pedigree.file, "\" does not exist.\n",
         "phi2share.q key 18")
  }

  ## Test pedindex.out for a supplied argument and existence
  if(missing(pedindex.out)) {
    stop("\nargument \"pedindex.out\" is missing with no default\n",
         "phi2share.q key 24")
  }
  if(!file.exists(pedindex.out)) {
    stop("\npedindex.out = \"", pedindex.out, "\" does not exist.\n",
         "phi2share.q key 28")
  }

  ## Test pedindex.cde for a supplied argument and existence
  if(missing(pedindex.cde)) {
    stop("\nargument \"pedindex.cde\" is missing with no default\n",
         "phi2share.q key 34")
  }
  if(!file.exists(pedindex.cde)) {
    stop("\npedindex.cde = \"", pedindex.cde, "\" does not exist.\n",
         "phi2share.q key 38")
  }

  ## 1) unpack phi2 (if necessary)
  phi2 <- gunzip(phi2, copy = TRUE)

  ## Read the pedigree file
  cat("Loading ", pedigree.file, "...\n", sep = "")
  ped <- read.table(pedigree.file, header = TRUE, sep = ",")
  
  dimnames(ped)[[2]] <- casefold(dimnames(ped)[[2]])

  if(is.null(ped$famid)) {
    remove.file(phi2);
    stop(paste("\nThe pedigree file '", pedigree.file, "' has no 'famid' ",
               "column.", "\nphi2share.q key 62", sep = ""))
  }

  ## Clear space for share.out or share.out.gz
  remove.file(paste(output.directory, "/share.out", sep = ""))
  remove.file(paste(output.directory, "/share.out.gz", sep = ""))

  ## Read phi2  
  data <- read.table(phi2, header = FALSE, sep = "")
  
  ## 3) swap first two columns of phi2
  data <- cbind( data[, 2], data[, 1], data[, 3:ncol(data)] )
  
  ## 4) order the table based on the sorted first two columns
  ## This ordering by the first column not only sorts the data based on the
  ## first column, but the data is also sorted on the second column as well
  ## because of how the original phi2 file was ordered.
  data <- data[ order(data[, 1]), ]
  
  ## 5) remove self reference rows
  data <- data[ data[,1] != data[,2], ]
  
  ## 6) load pedindex.out (pedindex.cde is needed to do this)
  cat("Loading ", pedindex.out, "...\n", sep = "")
  ids <- read.pedindex.out(pedindex.out, pedindex.cde)
  
  ## 7) convert all ids in first two columns of table to original ids using
  ##    pedindex.out
  unique.ids <- array(character(length(data[, 1]) * 2),
                      dim = c(length(data[, 1]), 2))
  unique.ids[, 1] <- ids[match(data[, 1], ids[, 1]), 4]
  unique.ids[, 2] <- ids[match(data[, 2], ids[, 1]), 4]

  ## 8) generate parent-parent, sib-sib and parent-offspring info
  cat("Calculating sibling-sibling...\n")
  siblings <- is.sibling(unique.ids[, 1], unique.ids[, 2], ped)
  cat("Calculating spouse-spouse...\n")
  spouses <- is.spouse(unique.ids[, 1], unique.ids[, 2], ped)
  cat("Calculating parent-offspring...\n")
  parent.offsprings <-
    is.parent.offspring(unique.ids[, 1], unique.ids[, 2], ped)
  
  ## 9) output ids, common genetics, and pp, ss, po info to SHARE.OUT
  share.out <- cbind(unique.ids[, 1], unique.ids[, 2], data[, 3],
                     siblings, spouses, parent.offsprings)

  if(using.R()) {
    write.table(share.out, 'share.out', row.names = FALSE, col.names = FALSE,
                sep = "\t", quote = FALSE)
  } else {
    write.table(share.out, 'share.out', dimnames.write = FALSE, sep = "\t")
  }

  ## If the output directory is not the current directory and the directory
  ## does not exist, create it.
  if( output.directory != '.' ) {
    if( !file.exists(output.directory) ) {
      multic.system(paste('mkdir -p ', output.directory, sep = ""))
    }
    ## If the output directory is not the current directory, move share.out
    ## to that output directory.
    multic.system(paste('mv share.out ', output.directory, sep = ""))
  }

   ##  unix(paste("gzip", phi2))
  ## The abov line was replaced after "copy = T" was added to the gunzip
  ## command at the beginning of the function.
  remove.file(phi2);
  gzip(paste(output.directory, "share.out", sep = "/"))
  
  invisible()
}

read.pedindex.out <- function(pedindex.out, pedindex.cde) {
  ## Read pedindex.cde for the locations of sequential and original
  ## identifiers
  if( !file.exists(pedindex.cde) ) {
    stop(paste("\nThe phi2 file '", pedindex.cde, "' does not exist.\n",
               "phi2share.q key 57", sep = ""))
  }

  if(using.R()) {
    ## Add code to dynamically calcuate the 51 below (the lengths of the
    ## lines).
    pedindex.cde.data <-
      read.fwf(pedindex.cde, widths = c(3, 51), skip = 1)
    pedindex.cde.data[, 2] <- as.character(pedindex.cde.data[, 2])
  } else {
    pedindex.cde.data <- importData(pedindex.cde, type='ASCII')
  }

  ## famid.index is the row of pedindex.cde that holds the id "FAMILY" in
  ## the second column
  famid.index <- pmatch('FAMILY', pedindex.cde.data[, 2])
  if( is.na(famid.index) ) {
    famid.index <- pmatch('PEDIGREE', pedindex.cde.data[, 2])
    if( is.na(famid.index) ) {
      stop(paste("\nThe file '", pedindex.cde, "' does not contain a column",
                 "name of 'FAMILY' or 'PEDIGREE'.\nphi2share.q key 71",
                 sep = ""))
    }
  }
  
  ## id.index is the row of pedindex.cde that holds the id "ID" in the
  ## second column
  id.index <- pmatch('ID', pedindex.cde.data[, 2])
  if( is.na(id.index) ) {
    stop(paste("\nThe file '", pedindex.cde, "' does not contain a column",
               "name of 'ID'.\nphi2share.q key 66", sep = ""))
  }

  ## Calculate the sum of adding all the rows of column 1 up to (but not
  ## including) the row that holds the id "PEDIGREE"
  famid.start <- sum(pedindex.cde.data[1:(famid.index-1), 1])
  
  ## Find how many characters are in the field labeled by "PEDIGREE"
  famid.len <- pedindex.cde.data[famid.index, 1]
    
  ## Calculate the sum of adding all the rows of column 1 up to (but not
  ## including) the row that holds the id "ID"
  origIdStart <- sum(pedindex.cde.data[seq(1, id.index-1), 1])
  
  ## Find how many characters are in the field labeled by "ID"
  origIdLen <- pedindex.cde.data[id.index, 1]
  
  ## ibdid.index is the row of pedindex.cde that holds the id "IBDID" in the
  ## second column
  ibdid.index <- pmatch('IBDID', pedindex.cde.data[, 2])
  if( is.na(id.index) ) {
    stop(paste("\nThe file '", pedindex.cde, "' does not contain\na column ",
               "name of 'IBDID'.\nphi2share.q key 81", sep = ""))
  }
  
  ## Find how many characters are in the field labeled by "IBDID"
  seqIdLen <- pedindex.cde.data[1, ibdid.index]
  
  ## Read pedindex.out for the sequential and original identifiers
  form <- paste("%", seqIdLen, "f,%", famid.start-seqIdLen, "*,%",
                famid.len, "s,%",
                origIdStart-(famid.start + famid.len), "*,%",
                origIdLen, "s", sep='')
  if(using.R()) {
    ids <- read.fwf(pedindex.out,
                    c(seqIdLen, famid.start-seqIdLen, famid.len,
                      origIdStart - (famid.start + famid.len), origIdLen ))
    ids <- ids[, c(1, 3, 5)]
    ids[, 2] <- as.character(ids[, 2])
    ids[, 3] <- as.character(ids[, 3])
  } else {
    ids <- importData(pedindex.out, type='FASCII', format=form)
  }

  ## Remove leading and trailing spaces
  ids[, 2] <- trim(ids[, 2])
  ids[, 3] <- trim(ids[, 3])

  ## Combine the famid and original id to create a universally unique id
  ids[, 4] <- paste(ids[, 2], ids[, 3], sep = "-")

  return (ids)
}

is.parent.offspring <- function(x, y, ped) {
  unique.ids <- paste(ped$famid, ped$id, sep = '-')
  
  x.index <- match(x, unique.ids)
  y.index <- match(y, unique.ids)
  
  if(any(ped$famid[x.index] != ped$famid[y.index])) {
    error.index <- match(TRUE, ped$famid[x.index] != ped$famid[y.index])
    stop(paste("\nThe person combination, ", ped$id[x.index][error.index],
               " (from family ", ped$famid[x.index][error.index], "), and ",
               ped$id[y.index][error.index], " (from family ",
               ped$famid[y.index][error.index],
               ")\nare not in the same family.\n",
               "They should either be in the same family, ",
               "or not a combination in the\n",
               "share.out file.\nphi2share.q key 111", sep = ""))
  }

  ## Is x the father of y
  x.is.father.of.y <- ped$id[x.index] == ped$fa[y.index]
  x.is.mother.of.y <- ped$id[x.index] == ped$mo[y.index]
  
  y.is.father.of.x <- ped$id[y.index] == ped$fa[x.index]
  y.is.mother.of.x <- ped$id[y.index] == ped$mo[x.index]

  return (as.integer(x.is.father.of.y
                     | x.is.mother.of.y
                     | y.is.father.of.x
                     | y.is.mother.of.x))
}

is.sibling <- function(x, y, ped) {
  unique.ids <- paste(ped$famid, ped$id, sep = '-')
  
  x.index <- match(x, unique.ids)
  y.index <- match(y, unique.ids)
  if(any(ped$famid[x.index] != ped$famid[y.index])) {
    error.index <- match(TRUE, ped$famid[x.index] != ped$famid[y.index])
    stop(paste("\nThe person combination, ", ped$id[x.index][error.index],
               " (from family ", ped$famid[x.index][error.index], "), and ",
               ped$id[y.index][error.index], " (from family ",
               ped$famid[y.index][error.index],
               ") are not in the same family.\n",
               "They should either be in the same family, ",
               "or not a combination in the\n",
               "share.out file.\nphi2share.q key 131", sep = ""))
  }

  father.of.x <- ped$fa[x.index]
  ## change the father id of x to -1 if the father id is 0
  father.of.x <- ifelse(father.of.x == 0, -1, father.of.x)
  father.of.y <- ped$fa[y.index]
  
  mother.of.x <- ped$mo[x.index]
  ## change the mother id of x to -1 if the father id is 0
  mother.of.x <- ifelse(mother.of.x == 0, -1, mother.of.x)
  mother.of.y <- ped$mo[y.index]
  
  return ( as.integer((father.of.x == father.of.y)
                      & (mother.of.x == mother.of.y)) )
}

is.spouse <- function(x, y, ped) {
  ## This family.sizes setup and stuff does not have to be updated to using
  ## get.family.sizes, because it will never have to worry about
  ## bootstrap'ed data.  Eric Lunde, 2005-09-07
  family.sizes <- table(ped$famid)
  family.sizes <- family.sizes[match(unique(ped$famid),
                                     as.character(names(family.sizes)))]
  unique.famids.count <- length(family.sizes)
  
  unique.ids <- paste(ped$famid, ped$id, sep = '-')
  unique.fa <- paste(ped$famid, ped$fa, sep = '-')
  unique.mo <- paste(ped$famid, ped$mo, sep = '-')
  
  x.index <- match(x, unique.ids)
  y.index <- match(y, unique.ids)
  
  if(any(ped$famid[x.index] != ped$famid[y.index])) {
    error.index <- match(TRUE, ped$famid[x.index] != ped$famid[y.index])
    stop(paste("\nThe person combination, ", ped$id[x.index][error.index],
               " (from family ", ped$famid[x.index][error.index], "), and ",
               ped$id[y.index][error.index], " (from family ",
               ped$famid[y.index][error.index],
               ") are not in the same family.\n",
               "They should either be in the same family, ",
               "or not a combination in the\n",
               "share.out file.\nphi2share.q key 176", sep = ""))
  }

  result <- .C('isSpouse',
               as.character(unique.ids[x.index]),
               as.character(unique.ids[y.index]),
               as.integer(length(x)),
               as.character(unique.fa),
               as.character(unique.mo),
               as.integer(length(ped$fa)),
               as.integer(family.sizes),
               as.integer(length(family.sizes)),
               result = integer(length(x)),
               PACKAGE = "multic"
               )$result
  
  return (result)
}
