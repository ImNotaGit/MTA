library(Matrix)
library(MASS)
library(Rcplex.my)
library(parallel)
library(Rcpp)
library(RcppArmadillo)

# need to source achr.cpp in the same dir with sourceCpp("achr.cpp")

sample.model <- function(model, params) {
  if ("sampl" %in% ls(model)) { # ls(model) works for model as either an environment or a list; but here we intend it to be an environment
    cat("Will use the warmup points and status stored in the model.\n")
    cat("Will sample", params$n.sampl, "points.\n")
    res <- achr(model, model$sampl$stat, model$sampl$warmup.pnts, params$n.sampl, params$steps.per.pnt)
    model$sampl$stat <- res$stat
    model$sampl$pnts <- cbind(model$sampl$pnts, res$sampl.pnts)
    model$sampl$mean.rng <- c(params$n.burnin+1, ncol(model$sampl$pnts))
    model$sampl$mean <- rowMeans(model$sampl$pnts[, model$sampl$mean.rng[1]:model$sampl$mean.rng[2]])
  } else {
    warmup.pnts <- sample.warmup.pnts(model, params$n.warmup, params$ncores)
    centr.pnt <- rowMeans(warmup.pnts)
    init.stat <- list(centr.pnt=centr.pnt, prev.pnt=centr.pnt, n.tot.steps=0)
    cat("Will sample", params$n.sampl, "points after", params$n.burnin, "burn-in points.\n")
    res <- achr(model, init.stat, warmup.pnts, params$n.burnin+params$n.sampl, params$steps.per.pnt)
    model$sampl <- list()
    model$sampl$warmup.pnts <- warmup.pnts
    model$sampl$stat <- res$stat
    model$sampl$pnts <- res$sampl.pnts
    model$sampl$mean.rng <- c(params$n.burnin+1, ncol(model$sampl$pnts))
    model$sampl$mean <- rowMeans(model$sampl$pnts[, model$sampl$mean.rng[1]:model$sampl$mean.rng[2]])
  }

  model
}

sample.warmup.pnts <- function(model, n, ncores) {
  n.rxns <- ncol(model$S)
  if (n<2*n.rxns) {
    n <- 2*n.rxns
    warning(sprintf("#{warmup points} should be at least 2*#{reactions}=%d.\n", 2*n.rxns))
  }
  cat("Will generate", n, "warmup points.\n")
  cat("Begin generating warmup points...\n")
  orth.pnts <- get.orth.pnts(model, n, ncores)
  rand.pnts <- get.rand.pnts(model, n, ncores)
  r <- rep(runif(n), each=n.rxns)
  dim(r) <- c(n.rxns, n)
  res <- orth.pnts*r + rand.pnts*(1-r)
  cat("Finished generating warmup points.\n")
  res
}

get.orth.pnts <- function(model, n, ncores) {
  n.rxns <- ncol(model$S)
  mat <- cbind(Diagonal(n.rxns), Diagonal(n.rxns, x=-1))
  if (n<=2*n.rxns) {
    mat <- mat[, sample(2*n.rxns, n)]
  } else {
    mat <- cbind(mat[, sample(2*n.rxns)], mat[, sample(2*n.rxns, n-2*n.rxns, replace=TRUE)])
  }
  cl <- makeCluster(ncores, type="FORK")
  res <- parApply(cl, mat, 2, get.opt.pnt, model=model)
  stopCluster(cl)
  res
}

get.rand.pnts <- function(model, n, ncores) {
  n.rxns <- ncol(model$S)
  cs <- runif(n.rxns*n) - 0.5
  dim(cs) <- c(n.rxns, n)
  cl <- makeCluster(ncores, type="FORK")
  res <- parApply(cl, cs, 2, get.opt.pnt, model=model)
  stopCluster(cl)
  res
}

get.opt.pnt <- function(model, c) {
  cvec <- c / norm(c,"2")
  objsense <- "max"
  Amat <- rbind(model$S, model$S)
  bvec <- c(model$rowlb, model$rowub)
  sense <- rep(c("G","L"), c(length(model$rowlb), length(model$rowub)))
  lb <- model$lb
  ub <- model$ub
  
  res <- rcplex(cvec=cvec, objsense=objsense, Amat=Amat, bvec=bvec, sense=sense, lb=lb, ub=ub, control=list(trace=0, maxcalls=5000, tilim=120, threads=1))
  res[[1]]$xopt
}

