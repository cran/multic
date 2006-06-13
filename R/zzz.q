.First.lib <- function(lib, pkg) {
  ## use cat() to include messages for the user 
  ## when loading the library
  no.kinship.lib.string <-
    paste("The kinship library could not be loaded.\n",
          "Not all of multic's functionality will be present.\n", sep = "")

  if(using.R()) {
    library.dynam("multic", pkg, lib)
    if( !require("kinship") ) {
      cat(no.kinship.lib.string)
    }
  } else {
    if(version$major < 7) {
      cat("loading kinship library\n")
      library("kinship")
      ## I would like to print a message to the user showing the
      ## no.kinship.lib.string, but if kinship cannot be loaded, .First.lib
      ## is exited immediately and I cannot display the message like I can
      ## in R.
      cat("\nIf the kinship library could not be loaded\n",
          "(you would have seen an error message),\n",
          "Not all of multic's functionality will be present.\n",sep = "")
    }
  }
}

## The S-PLUS version of .Last.lib requres section, .data, and where, the R
## version does not.  They are added with defaults to work in both
## environments.
.Last.lib <- function(libpath, section = NULL, .data = NULL, where = NULL) {
  if(using.R()) {
    library.dynam.unload("multic", libpath)
  }
}
