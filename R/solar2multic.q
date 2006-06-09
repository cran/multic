##############################################################################
## Title: solar2multic
## Description: solar2multic is a utility function to convert the (m)ibd
##              (identity by descent) files created by SOLAR into the multic
##              input file 'mloci.out' and convert the 'phi2.gz' created by
##              SOLAR into the multic input file 'share.out'.
## Input: phi2 - character variable specifying the name of the phi2 file
##               (the phi2 file is preferred to be in 'gzip' format)
##        pedigree.file - character variable specifying the name of the file
##                        used to describe the pedigree structure.  The file
##                        needs to be a comma value seperated file with a
##                        header of 'famid', 'id', 'fa', 'mo', and 'sex' in no
##                        particular order.
##        pedindex.out - character variable specifying the name of the
##                       pedindex.out file that was output from SOLAR.
##        pedindex.cde - character variable specifying the name of the
##                       pedindex.cde file that was output from SOLAR.
##        ibd.directory - character variable specifying the name of the
##                        directory that contains the (m)ibd files (the (m)ibd
##                        files are preferred to be in 'gzip' format)
##        ibd.dist - character variable specifying a path to a .dist file
##                   created by SOLAR.  This file links the names of markers
##                   from ibd files to their cM value.
##        output.directory - character variable specifying the name of the
##                           directory that the output files ('mloci.out' and
##                           'share.out') will be moved to after they've been
##                           created
## Output: NONE.
## Side Effects: The (m)ibd files need to be 'fixed' before we can compile
##               them into an 'mloci.out' file.  So we place the 'fix'ed
##               (m)ibd file in a directory with the same name as the one
##               specified in ibd.directory but with an added "_fixed" to the
##               directory name containing the 'fix'ed (m)ibds.
##
##               During the 'fix'ing process, a file named 'temp_ibd' is
##               created and is deleted at the end.
##
##               In the creation of 'share.out' many temporary files (e.g.
##               'tempshare.out, one.phi2, phi2.sorted') are generated and
##               deleted at the end.
##
##               The final files 'mloci.out' and 'share.out' are created in
##               the same directory that the Splus session is in.  At the end
##               of solar2multic, those files are moved using the UNIX command
##               'mv' to the specified output directory.  This will cause an
##               issue of files overwriting each other if we ever want to run
##               more than one solar2multic conversion simultaneously.
##
##               If any of the input files are not present, solar2multic will
##               print an error message and quit execution.
## Author: Eric Lunde 04/05/2004
## Updates (Date, Modified By, Changes made) :
## 04/20/2004, Eric Lunde, I incorporated the new phi2share Splus function
##                         into solar2multic.  I also added pedindex.out and
##                         pedindex.cde to the parameter list of solar2multic.
##############################################################################
solar2multic <- function(phi2, pedigree.file, pedindex.out,
                         pedindex.cde, ibd.directory, ibd.dist = NULL,
                         output.directory = ".",
                         delete.fixed.dir = TRUE)
{
  ## Check the inputs to see if there are in the correct (character) mode
  if( !is.character(pedigree.file) ) {
    if( is.null(pedigree.file) ) {
      stop(paste("\nThe pedigree file name '", pedigree.file,
                 "' was given a valueless (NULL) input.\n",
                 "solar2multic.q key 11", sep = ""))
    }
    warning(paste("The pedigree file name '", pedigree.file,
                  "' was not in character format.\n",
                  "solar2multic.q key 4", sep = ""))
    pedigree.file <- as.character(pedigree.file)
  }  
  if( !is.character(pedindex.out) ) {
    if( is.null(pedindex.out) ) {
      stop(paste("\nThe pedindex.out file name '", pedindex.out,
                 "' was given a valueless (NULL) input.\n",
                 "solar2multic.q key 66", sep = ""))
    }
    warning(paste("The pedindex.out file name '", pedindex.out,
                  "' was not in character format.\n",
                  "solar2multic.q key 70", sep = ""))
    pedindex.out <- as.character(pedindex.out)
  }  
  if( !is.character(pedindex.cde) ) {
    if( is.null(pedindex.cde) ) {
      stop(paste("\nThe pedindex.cde file name '", pedindex.cde,
                 "' was given a valueless (NULL) input.\n",
                 "solar2multic.q key 77", sep = ""))
    }
    warning(paste("The pedindex.cde file name '", pedindex.cde,
                  "' was not in character format.\n",
                  "solar2multic.q key 81", sep = ""))
    pedindex.cde <- as.character(pedindex.cde)
  }  
  if( !is.character(pedindex.cde) ) {
    if( is.null(pedindex.cde) ) {
      stop(paste("\nThe pedindex.cde file name '", pedindex.cde,
                 "' was given a valueless (NULL) input.\n",
                 "solar2multic.q key 77", sep = ""))
    }
    warning(paste("The pedindex.cde file name '", pedindex.cde,
                  "' was not in character format.\n",
                  "solar2multic.q key 81", sep = ""))
    pedindex.cde <- as.character(pedindex.cde)
  }  
  if( !missing(ibd.dist) && !is.character(ibd.dist) ) {
    if( is.null(ibd.dist) ) {
      stop(paste("\nThe ibd distance file name '", ibd.dist,
                 "' was given a valueless (NULL) input.\n",
                 "solar2multic.q key 111", sep = ""))
    }
    warning(paste("The ibd distance file name '", ibd.dist,
                  "' was not in character format.\n",
                  "solar2multic.q key 115", sep = ""))
    ibd.dist <- as.character(ibd.dist)
  }  
  if( !is.character(phi2) ) {
    if( is.null(phi2) ) {
      stop(paste("\nThe phi2 file name '", phi2,
                 "' was given a valueless (NULL) input.\n",
                 "solar2multic.q key 33", sep = ""))
    }
    warning(paste("The phi2 file name '", phi2,
                  "' was not in character format.\n",
                  "solar2multic.q key 16", sep = ""))
    phi2 <- as.character(phi2)
  }  
  if( !is.character(output.directory) ) {
    if( is.null(output.directory) ) {
      stop(paste("\nThe output directory name '", output.directory,
                 "' was given a valueless (NULL) input.\n",
                 "solar2multic.q key 42", sep = ""))
    }
    warning(paste("The output directory name '", output.directory,
                  "' was not in character format.\n",
                  "solar2multic.q key 46", sep = ""))
    output.directory <- as.character(output.directory)
  }
  
  ## Check to see if all the input file names represent existing files.
  if( !file.exists(pedigree.file) ) {
    stop(paste("\nThe pedigree file '", pedigree.file, "' does not exist.\n",
               "solar2multic.q key 24", sep = ""))
  }
  if( !file.exists(pedindex.out) ) {
    stop(paste("\nThe pedindex.out file '", pedindex.out,
               "' does not exist.\nsolar2multic.q key 125", sep = ""))
  }
  if( !file.exists(pedindex.cde) ) {
    stop(paste("\nThe pedindex.cde file '", pedindex.cde,
               "' does not exist.\nsolar2multic.q key 129", sep = ""))
  }
  if( !file.exists(ibd.directory) ) {
    stop(paste("\nThe ibd directory '", ibd.directory, "' does not exist.\n",
               "solar2multic.q key 29", sep = ""))
  }

  ## Both solar2mloci and phi2share worry about this.  I don't think
  ## solar2multic has to.
  if(F) {
    ## Unpack phi2 (if necessary)
    phi2 <- gunzip(phi2)
  }
  
  ## Then run the solar2mloci function on the fixed directory
  solar2mloci(ibd.directory, phi2, pedindex.out, pedindex.cde,
            ibd.dist = ibd.dist,
            delete.fixed.dir = delete.fixed.dir)
  
  ## Create share.out from phi2 and the pedigree file 
  ##  create.share.out(phi2, pedigree.file)
  phi2share(phi2, pedigree.file, pedindex.out, pedindex.cde)

  ## If the output directory is not the current directory and the directory
  ## does not exist, create it.
  if( output.directory != '.' ) {
    if( !file.exists(output.directory) ) {
      multic.system(paste('mkdir -p', output.directory))
    }
    ## If the output directory is not the current directory, move share.out
    ## and mloci.out to that output directory.
    multic.system(paste('mv share.out.gz mloci.out.gz', output.directory))
  }

  cat(paste(output.directory, "/mloci.out.gz and ", output.directory,
            "/share.out.gz\nwere created successfully.\n", sep = ""))

  invisible()
}
