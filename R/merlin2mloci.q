##
# Converts a merlin merlin.ibd (or merlin.kin) file to a multic mloci.out file
#
merlin2mloci <- function (merlin.ibd.file=NULL, merlin.kin.file=NULL, output.directory=".") {
  
  if (!file.exists(output.directory)) {
    stop ("ERROR in merlin2mloci:  output.directory ", output.directory, " not found.\n")
  }

  mloci.out <- paste(output.directory, "mloci.out", sep="/")

  if (file.exists(mloci.out)) {
    paste(mloci.out, " already exists.  Renaming to ", mloci.out, ".old\n", sep="")
    file.copy(mloci.out, paste(mloci.out, ".old", sep=""), overwrite=TRUE)
    file.remove(mloci.out)
  }

  ## TODO:  Sort merlin input file with Unix system call?
  
  if (is.null(merlin.ibd.file) && is.null(merlin.kin.file)) {
    stop ("ERROR in merlin2mloci:  Must specify either a merlin.ibd or merlin.kin file.\n");
    traceback()
  }
  if (!is.null(merlin.kin.file)) {
    merlin.kin <- read.table(merlin.kin.file, header=TRUE)

    ## TODO:  Check dimnames[[2]]
  }
  
  ## We need to create the merlin.kin table from the merlin.ibd file
  else {
    merlin.ibd <- read.table("merlin.ibd", header=TRUE)

    ## TODO:  Check dimnames[[2]]

    ## The kinship coefficient (phi) is 0.25*P1 + 0.5*P2 in merlin.ibd file
    ## This variable is named in all caps so that the dimname mathese the merlin.kin file format 
    KINSHIP <- 0.25 * merlin.ibd[,"P1"] + 0.5 * merlin.ibd[,"P2"]
    merlin.kin <- cbind(merlin.ibd[, 1:4], KINSHIP)

    # These objects are large and we don't need them anymore
    rm (merlin.ibd, KINSHIP)
  }

                                        # Make the kinship matrix before this?
                                        #id1 <- paste(famid, merlin.ibd[,"ID1"], sep="-")

  ## Weed out the duplicates
  merlin.kin <- merlin.kin[merlin.kin[, "ID1"] != merlin.kin[, "ID2"], ]

#  if (fix.cols) {
    ## Looping in R is way too slow.  Wrote a C++ function instead
#    for (ii in 1:nrow(merlin.kin)) {
     
#    if (merlin.kin[ii, "ID1"] > merlin.kin[ii, "ID2"]) {
#        tmp <- merlin.kin[ii, "ID1"]
#        merlin.kin[ii, "ID1"] <- merlin.kin[ii, "ID2"]
#        merlin.kin[ii, "ID2"] <- tmp
#      }
#    }
##    dyn.load("sortIdPairs.so")
#    .Call("sortIdPairs", merlin.kin[, "ID1"], merlin.kin[, "ID2"])
    
#  }

  ## sortIdPairs will switch the id pair order if ID1 > ID2
  #dyn.load("sortIdPairs.so")
  .Call("sortIdPairs", merlin.kin[, "ID1"], merlin.kin[, "ID2"], PACKAGE = "multic")
  
  ## Grab the marker names
  mibd.nums <- unique(merlin.kin[, "MARKER"])

  ## mloci.out needs 2 * phi
  merlin.kin[, "KINSHIP"] <- 2 * merlin.kin[, "KINSHIP"]
    
  for (ii in 1:length(mibd.nums)) {
    this.mibd <- merlin.kin[merlin.kin[, "MARKER"] == mibd.nums[ii], ]

    ## Drop MARKER column (4th col)
    this.mibd <- this.mibd[, -4]

    ## Sort this.mibd by family, id1, and id2
    idx <- order(this.mibd[, "FAMILY"], this.mibd[, "ID1"], this.mibd[, "ID2"])
    this.mibd <- this.mibd[idx, ]

    ## Paste FAMILY and ID's together delimited with '-'
    this.mibd <- collapseMibd(this.mibd)

    ## Write (append) it to the mloci.out file
    mibd.comment <- paste("# mibd.", mibd.nums[ii], sep="")

    cat("Writing ", mibd.comment, " to ", mloci.out, "\n", sep="")
    
    ## TODO:  Write to tempfiles and then use Unix cat
    write (mibd.comment, file=mloci.out, append=TRUE)
    write.table(this.mibd, file=mloci.out, append=TRUE, quote=FALSE,
                sep="\t", row.names=FALSE, col.names=FALSE)
  } # end foreach mibd
  
} # end function merlin2mloci
