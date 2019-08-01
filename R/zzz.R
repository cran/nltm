.onUnload <- function (libpath) {
  library.dynam.unload("nltm", libpath)
}
