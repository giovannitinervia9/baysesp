.onLoad <- function(libname, pkgname) {
  # registra i metodi S7, compresi quelli per generici di altri pacchetti
  S7::methods_register()
}
