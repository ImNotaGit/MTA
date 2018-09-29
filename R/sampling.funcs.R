library(Matrix)
library(MASS)
library(Rcplex)
library(parallel)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("achr.cpp")

sample.model <- function(model, params) {
  warmup.pnts <- sample.warmup.pnts(model, params$n.warmup)
  centr.pnt <- rowMeans(warmup.pnts)
  init.stat <- list(centr.pnt=centr.pnt, prev.pnt=centr.pnt, n.tot.steps=0)
  stat <- achr(model, init.stat, warmup.pnts, params$n.burnin, params$steps.per.pnt)$stat
  res <- achr(model, stat, warmup.pnts, params$n.sampl, params$steps.per.pnt)
  # model is an environment, modified in place outsite this function
  model$sample <- res
  model$v.ref <- rowMeans(res$sampl.pnts)
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
  cl <- makeCluster(detectCores(), type="FORK")
  parApply(cl, mat, 2, get.opt.pnt, model=model)
  stopCluster(cl)
}

get.rand.pnts <- function(model, n) {
  n.rxns <- length(model$rxns)
  cs <- runif(n.rxns*n) - 0.5
  dim(cs) <- c(n.rxns, n)
  cl <- makeCluster(detectCores(), type="FORK")
  parApply(cl, cs, 2, get.opt.pnt, model=model)
  stopCluster(cl)
}

get.opt.pnt <- function(model, c) {
  cvec <- c / norm(c,"2")
  objsense <- "max"
  Amat <- rbind(model$S, model$S)
  bvec <- c(model$rowlb, model$rowub)
  sense <- rep(c("G","L"), c(length(model$rowlb), length(model$rowub)))
  lb <- model$lb
  ub <- model$ub
  
  res <- Rcplex(cvec=cvec, objsense=objsense, Amat=Amat, bvec=bvec, sense=sense, lb=lb, ub=ub)
  if (res$status!=1) warning("Potential problems running LP. Solver status: ", res$status, ".\n")
  res$xopt
}

