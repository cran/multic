fixibd <- function(directory.name, fixed.directory.name, phi2.file,
                   directory.name.to.print) {
  ## Check to make sure the file names are in character mode.
  if( !is.character(directory.name) ) {
    warning("The directory name was not in character format.\n",
            "create.mloci.q key 4")
    directory.name <- as.chacter(directory.name)
  }  
  if( !is.character(phi2.file) ) {
    warning("The phi2 file was not in character format.\n",
            "create.mloci.q key 10")
    phi2.file <- as.chacter(phi2.file)
  }  
  if( !is.character(fixed.directory.name) ) {
    warning("The fixed directory name was not in character format.\n",
            "create.mloci.q key 16")
    fixed.directory.name <- as.chacter(fixed.directory.name)
  }

  ## If the output directory already exists, delete it and its contents
  if( file.exists(fixed.directory.name) ) {
    multic.system(paste('rm -rf', fixed.directory.name))
  }

  .C("fixibd",
     directory.name = as.character(directory.name),
     fixed.directory.name = as.character(fixed.directory.name),
     phi2.file = as.character(phi2.file),
     directory.name.to.print = as.character(directory.name.to.print),
     PACKAGE = "multic"
     )

  invisible()
}
