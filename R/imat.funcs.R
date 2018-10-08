library(Matrix)
library(Rcplex.my)
source("utils.R")
source("sampling.funcs.R")

imat.pars <- list(flux.act=1, flux.inact=0.1, flux.bound=1000)
milp.pars <- list(trace=1, maxcalls=5000, tilim=120, nodesel=0)
sampl.pars <- list(n.warmup=5000, n.burnin=1000, n.sampl=2000, steps.per.pnt=400)

imat <- function(model, expr, imat.params=imat.pars, milp.params=milp.pars, sampl.params=sampl.pars) {
  
  # model as environment
  model <- as.environment(model)

  # formulate iMAT model
  imat.model <- form.imat(model, expr, imat.params)

  # run the iMAT MILP
  run.imat(imat.model, milp.params) # modify imat.model in place
  
  # update the original metabolic model based on iMAT result
  update.model(model, imat.model, imat.params) # update model in place

  # sample the metabolic model to get the fluxes of the reference state
  sample.model(model, sampl.params) # update model in place
  
  # force close CPLEX
  Rcplex.close()

  # model back to list
  as.list(model)
}

form.imat <- function(model, expr, params) {

  # reaction data
  rxns.int.raw <- discrt.genes2rxns(expr, "stat", model)
  rxns.int <- rxns.int.raw
  # remove integers for dead end rxns
  rxns.int[model$lb==0 & model$ub==0] <- 0

  n.mets <- nrow(model$S)
  n.rxns <- ncol(model$S)

  # 1. Active reactions: specify the y+ indicator variables, representing activation in the forward direction (i.e. v>flux.act)
  rxns.act <- which(rxns.int==1)
  n.act <- length(rxns.act)
  m1 <- sparseMatrix(1:n.act, rxns.act, dims=c(n.act, n.rxns))
  m2 <- Diagonal(n.act, x=(-params$flux.act-params$flux.bound))
  S <- rbind(cbind(model$S, sparseMatrix(NULL, NULL, dims=c(n.mets, n.act))), cbind(m1, m2))

  # 2. Reversible active reactions: for those reversible ones among the active reactions, specify the extra y- indicator variables, representing activation in the backward direction (i.e. v<-flux.act)
  # thus, an reversible active reaction has both the y+ and y- indicator variables, because it can be active in either direction (but never both, i.e. 1 XOR 2)
  rxns.act.rev <- which(rxns.int==1 & model$lb<0)
  n.act.rev <- length(rxns.act.rev)
  m1 <- sparseMatrix(1:n.act.rev, rxns.act.rev, dims=c(n.act.rev, ncol(S)))
  m2 <- Diagonal(n.act.rev, x=params$flux.act+params$flux.bound)
  S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(nrow(S), n.act.rev))), cbind(m1, m2))

  # 3. Inactive reactions: specify the y0 indicator variables
  # 3a. specify inactivation in the forward direction (i.e. v<flux.inact)
  rxns.inact <- which(rxns.int==-1)
  n.inact <- length(rxns.inact)
  m1 <- sparseMatrix(1:n.inact, rxns.inact, dims=c(n.inact, ncol(S)))
  m2 <- Diagonal(n.inact, x=params$flux.bound-params$flux.inact)
  # 3b. for those reversible inactive reactions, need to further specify inactivation in the backward direction (i.e. v>-flux.inact)
  # note that a reversible inactive reaction has only one y0 indicator variable, because for these reactions we want -flux.inact<v<flux.inact (3a AND 3b) 
  rxns.inact.rev <- which(rxns.int==-1 & model$lb<0)
  n.inact.rev <- length(rxns.inact.rev)
  m3 <- sparseMatrix(1:n.inact.rev, rxns.inact.rev, dims=c(n.inact.rev, ncol(S)))
  m4 <- sparseMatrix(1:n.inact.rev, match(rxns.inact.rev, rxns.inact), x=params$flux.inact-params$flux.bound, dims=c(n.inact.rev, n.inact))
  S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(nrow(S), n.inact))), cbind(m1, m2), cbind(m3, m4))

  # other parameters
  n <- n.act + n.act.rev + n.inact + n.inact.rev
  rowlb <- c(model$b, rep(-params$flux.bound, n))
  rowub <- c(model$b, rep(params$flux.bound, n))
  n <- ncol(S) - n.rxns
  c <- rep(c(0,1), c(n.rxns, n))
  vtype <- ifelse(c==1, "I", "C")
  lb <- c(model$lb, rep(0, n))
  ub <- c(model$ub, rep(1, n))
  var.ind <- rep(c("v","y+","y-","y0"), c(n.rxns, n.act, n.act.rev, n.inact)) # iMAT variable type indicators (v: fluxex; y+/-/0: indicator variables)

  # return iMAT model as an environment
  as.environment(list(genes.int=expr, rxns.int.raw=rxns.int.raw, rxns.int=rxns.int,
                      rxns.act=rxns.act, rxns.act.rev=rxns.act.rev, rxns.inact=rxns.inact, rxns.inact.rev=rxns.inact.rev, var.ind=var.ind,
                      c=c, S=S, rowlb=rowlb, rowub=rowub, lb=lb, ub=ub, vtype=vtype))
}

run.imat <- function(model, params) {
  cvec <- model$c
  objsense <- "max"
  Amat <- rbind(model$S, model$S)
  bvec <- c(model$rowlb, model$rowub)
  sense <- rep(c("G","L"), c(length(model$rowlb), length(model$rowub)))
  lb <- model$lb
  ub <- model$ub
  vtype <- model$vtype
  
  res <- Rcplex(cvec=cvec, objsense=objsense, Amat=Amat, bvec=bvec, sense=sense, lb=lb, ub=ub, vtype=vtype, control=params)
  if (!res$status %in% c(101,102)) stop("iMAT: Potential problems running MILP. Solver status: ", res$status, ".\n")

  model$milp.out <- res # model is an environment, so will be modified outside this function
}

update.model <- function(model, imat.res, params) {

  yp <- imat.res$milp.out$xopt[imat.res$var.ind=="y+"]
  ym <- imat.res$milp.out$xopt[imat.res$var.ind=="y-"]
  y0 <- imat.res$milp.out$xopt[imat.res$var.ind=="y0"]

  fw <- imat.res$rxns.act[yp==1]
  bk <- imat.res$rxns.act.rev[ym==1]
  inact <- imat.res$rxns.inact[y0==1]
  inact.rev <- intersect(inact, imat.res$rxns.inact.rev)
  
  # update model: model is an environment, so will be modified outside this function
  model$lb[fw] <- params$flux.act
  model$ub[bk] <- -params$flux.act
  model$ub[inact] <- params$flux.inact
  model$lb[inact.rev] <- -params$flux.inact
}

