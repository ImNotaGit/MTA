library(Matrix)
library(Rcplex)
library(parallel)

sample.model <- function(model, params) {
  warmup.pnts <- sample.warmup.pnts(model, params$n.warmup, 1)
  state <- model_sampling_ACHR(model, warmup.pnts, params$n.sampl.state, params$steps.per.pnt)$state
  model_sampling_ACHR(model, warmup.pnts, params$n.sampl, params$steps.per.pnt, state)$sample_points
}

sample.warmup.pnts <- function(model, n) {
  orth.pnts <- get.orth.pnts(model, n)
  rand.pnts <- get.rand.pnts(model, n)
  n.rxns <- length(model$rxns)
  r <- rep(runif(n), each=n.rxns)
  dim(r) <- c(n.rxns, n)
  orth.pnts*r + rand.pnts*(1-r)
}

get.orth.pnts <- function(model, n) {
  n.rxns <- length(model$rxns)
  mat <- cbind(Diagonal(n.rxns), Diagonal(n.rxns, x=-1))
  if (n<=2*n.rxns) {
    mat <- mat[, sample(2*n.rxns, n)]
  } else {
    mat <- cbind(mat[, sample(2*n.rxns)], mat[, sample(2*n.rxns, n-2*n.rxns, replace=TRUE)])
  }
  apply(mat, 2, get.opt.pnt, model=model)
}

get.rand.pnts <- function(model, n) {
  n.rxns <- length(model$rxns)
  cs <- runif(n.rxns*n) - 0.5
  dim(cs) <- c(n.rxns, n)
  apply(cs, 2, get.opt.pnt, model=model)
}

get.opt.pnt <- function(model, c) {
  cvec <- c / sqrt(sum(c^2))
  objsense <- "max"
  Amat <- rbind(model$S, model$S)
  bvec <- c(model$rowlb, model$rowub)
  sense <- rep(c("G","L"), c(length(model$rowlb), length(model$rowub)))
  lb <- model$lb
  ub <- model$ub
  
  res <- Rcplex(cvec=cvec, objsense=objsense, Amat=Amat, bvec=bvec, sense=sense, lb=lb, ub=ub, control=list(maxcalls=Inf))
  if (res$status!=1) stop("Failed running LP. Solver status: ", res$status, ".\n")
  res$xopt
}

