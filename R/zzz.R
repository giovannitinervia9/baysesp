.onLoad <- function(libname, pkgname) {
  # register S7 methods, including those for generics of other packages
  S7::methods_register()
}
