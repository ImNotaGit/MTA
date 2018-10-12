library(Matrix)
library(MASS)
library(Rcplex.my)
library(parallel)
library(Rcpp)
library(RcppArmadillo)

# need to source achr.cpp in the same dir with sourceCpp("achr.cpp")

sample.model <- function(model, params) {
  if ("sample" %in% ls(model)) {
    cat("Will use the warmup points and status stored in the model.\n")
    cat("Will sample", params$n.sampl, "points after", params$n.burnin, "burn-in points.\n")
    res <- achr(model, model$sample$stat, model$sample$warmup.pnts, params$n.burnin+params$n.sampl, params$steps.per.pnt)
    # model is an environment, modified in place outside this function
    model$sample$stat <- res$stat
    model$sample$sampl.pnts <- cbind(model$sample$sampl.pnts, res$sampl.pnts)
    model$sample$using.last.n.pnts <- params$n.sampl
    model$sample$v.ref <- rowMeans(res$sampl.pnts[, (params$n.burnin+1):params$n.sampl])
  } else {
    warmup.pnts <- sample.warmup.pnts(model, params$n.warmup)
    centr.pnt <- rowMeans(warmup.pnts)
    init.stat <- list(centr.pnt=centr.pnt, prev.pnt=centr.pnt, n.tot.steps=0)
    cat("Will sample", params$n.sampl, "points after", params$n.burnin, "burn-in points.\n")
    res <- achr(model, init.stat, warmup.pnts, params$n.burnin+params$n.sampl, params$steps.per.pnt)
    # model is an environment, modified in place outside this function
    model$sample <- res
    model$sample$warmup.pnts <- warmup.pnts
    model$sample$using.last.n.pnts <- params$n.sampl
    model$sample$v.ref <- rowMeans(res$sampl.pnts[, (params$n.burnin+1):params$n.sampl])
  }
}

sample.warmup.pnts <- function(model, n) {
  n.rxns <- length(model$rxns)
  if (n<2*n.rxns) {
    n <- 2*n.rxns
    warning(sprintf("#{warmup points} should be at least 2*#{reactions}=%d.\n", 2*n.rxns))
  }
  cat("Will generate", n, "warmup points.\n")
  cat("Begin generating warmup points...\n")
  orth.pnts <- get.orth.pnts(model, n)
  rand.pnts <- get.rand.pnts(model, n)
  r <- rep(runif(n), each=n.rxns)
  dim(r) <- c(n.rxns, n)
  res <- orth.pnts*r + rand.pnts*(1-r)
  cat("Finished generating warmup points.\n")
  res
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
  res <- parApply(cl, mat, 2, get.opt.pnt, model=model)
  stopCluster(cl)
  res
}

get.rand.pnts <- function(model, n) {
  n.rxns <- length(model$rxns)
  cs <- runif(n.rxns*n) - 0.5
  dim(cs) <- c(n.rxns, n)
  cl <- makeCluster(detectCores(), type="FORK")
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
  
  res <- Rcplex(cvec=cvec, objsense=objsense, Amat=Amat, bvec=bvec, sense=sense, lb=lb, ub=ub, control=list(trace=0, maxcalls=5000, tilim=120, threads=1))
  if (res$status!=1) warning("Potential problems running LP. Solver status: ", res$status, ".\n")
  res$xopt
}

