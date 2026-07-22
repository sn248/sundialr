#.onAttach <- function(...){
#  packageStartupMessage(
#    "SUNDIALS - Copyright (c) 2002-2020, Lawrence Livermore National Security and Southern Methodist University.\nAll rights reserved.")
#  }

# Register the sundialr C API (see src/sundialr_capi.cpp, src/capi_register.cpp)
# so other packages can reach CVODE via R_GetCCallable("sundialr", ...). Rcpp
# owns R_init_sundialr, so registration is done here instead. `.register_capi`
# is the dot-prefixed, unexported wrapper generated for register_capi().
.onLoad <- function(libname, pkgname) {
  .Call("_sundialr_register_capi", PACKAGE = "sundialr")
}
