library(Rcpp)
library(RcppArmadillo)
# need to source mep.cpp in the same dir with sourceCpp("mep.cpp")


run.mep <- function(model, params) {

  mep.res <- mep(model, params$beta, params$damp, params$max.iter, params$dlb, params$dub, params$epsil, params$fix.flux, params$fflux.id, params$fflux.mean, params$fflux.var)
  
  model$mep.res <- mep.res

  model
}
