.onLoad <- function(lib, pkg) {
    library.dynam("multic", pkg, lib)
}

.onUnload <- function(libpath, section = NULL, .data = NULL, where = NULL) {
    library.dynam.unload("multic", libpath)
}
