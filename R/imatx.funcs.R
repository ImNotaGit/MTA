library(Matrix)
library(data.table)
library(Rcplex.my)
# need to source utils.R
# need to source sampling.funcs.R

imat.pars <- list(flux.act=1, flux.inact=0.1, flux.delta.rel=0.1, flux.delta=0.1, flux.bound=1000)
milp.pars <- list(trace=1, nodesel=0, solnpoolagap=0, solnpoolgap=0, solnpoolintensity=2, n=1)
sampl.pars <- list(n.warmup=5000, n.burnin=1000, n.sampl=2000, steps.per.pnt=400, ncores=1L)

imatx <- function(model, expr, dflux, imat.params=imat.pars, milp.params=milp.pars, sampl.params=sampl.pars) {
  
  # formulate iMAT model
  imat.model <- form.imatx(model, expr, dflux, imat.params)

  # run the iMAT MILP
  run.imat(imat.model, milp.params) # modify imat.model in place
  
  # update the original metabolic model based on iMAT result
  res.model <- update.model.imatx(imat.model, sol=1, imat.params)

  if (!is.null(sampl.params)) {
    # sample the metabolic model to get the fluxes of the reference state
    sample.model(res.model, sampl.params) # update res.model in place
  }
  
  # close CPLEX
  Rcplex.close()

  # return
  res.model
}

imatx2steps <- function(model, expr, dflux, imat.params=imat.pars, milp.params=milp.pars, sampl.params=sampl.pars) {
  
  # formulate iMAT model step 1 (i.e. the original imat)
  imat.model1 <- form.imat(model, expr, imat.params)
  # run the iMAT MILP
  run.imat(imat.model1, milp.params) # modify imat.model in place
  
  # update the original metabolic model based on iMAT result
  imat.model2 <- update.model.imat(model, imat.model1, sol=1, imat.params)
  # formulate iMAT model step 2 (xde) upon the updated model from the above step
  imat.model2 <- form.imat.xde(imat.model2, dflux, imat.params)
  
  # run the iMAT MILP
  run.imat(imat.model2, milp.params) # modify imat.model2 in place
  
  # update the original metabolic model based on iMAT result
  res.model <- update.model.imat.xde(imat.model2, sol=1, imat.params)
  
  if (!is.null(sampl.params)) {
    # sample the metabolic model to get the fluxes of the reference state
    sample.model(res.model, sampl.params) # update model in place
  }
  
  # close CPLEX
  Rcplex.close()
  
  # return
  res.model
}

form.imat <- function(model, expr, params) {

  # reaction data
  rxns.int.raw <- exprs2rxns(expr, 0, model)
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
  rowlb <- c(model$rowlb, rep(-params$flux.bound, n)) # originally, it was c(model$b, ...)
  rowub <- c(model$rowub, rep(params$flux.bound, n)) # originally, it was c(model$b, ...)
  n <- ncol(S) - n.rxns
  c <- rep(c(0,1), c(n.rxns, n))
  vtype <- ifelse(c==1, "I", "C")
  lb <- c(model$lb, rep(0, n))
  ub <- c(model$ub, rep(1, n))
  var.ind <- rep(c("v","y+","y-","y0"), c(n.rxns, n.act, n.act.rev, n.inact)) # iMAT variable type indicators (v: fluxex; y+/-/0: indicator variables)

  # return iMAT model as an environment
  as.environment(list(irxn.ids=model$irxn.ids, # if not exist, will be NULL
                      genes.int=expr, rxns.int.raw=rxns.int.raw, rxns.int=rxns.int,
                      rxns.act=rxns.act, rxns.act.rev=rxns.act.rev, rxns.inact=rxns.inact, rxns.inact.rev=rxns.inact.rev, var.ind=var.ind,
                      c=c, S=S, rowlb=rowlb, rowub=rowub, lb=lb, ub=ub, vtype=vtype))
}

form.imat.xde0 <- function(model, i1, i2, df, rr, params) {
  S <- model$S
  # z+
  if (df>0) {
    S <- rbind(cbind(S, sparseMatrix(NULL,NULL,dims=c(nrow(S),1))),
               sparseMatrix(rep(1,3), c(i1, i2, ncol(S)+1), x=c(1, params$flux.delta.rel-1, (2-params$flux.delta.rel)*params$flux.bound), dims=c(1,ncol(S)+1)))
  } else if (df<0) {
    S <- rbind(cbind(S, sparseMatrix(NULL,NULL,dims=c(nrow(S),1))),
               sparseMatrix(rep(1,3), c(i1, i2, ncol(S)+1), x=c(params$flux.delta.rel-1, 1, (2-params$flux.delta.rel)*params$flux.bound), dims=c(1,ncol(S)+1)))
  }
  if (df!=0) { # for now I haven't implemented df==0, so when df==0 should do nothing
    model$rowlb <- c(model$rowlb, (params$flux.delta.rel-2)*params$flux.bound)
    model$rowub <- c(model$rowub, (2-params$flux.delta.rel)*params$flux.bound-params$flux.delta)
    model$lb <- c(model$lb, 0)
    model$ub <- c(model$ub, 1)
    model$c <- c(model$c, 1)
    model$vtype <- c(model$vtype, "I")
    model$var.ind <- c(model$var.ind, "z+")
  }
  # reversible reactions
  if (rr) {
    if (df>0) {
      # additional z+
      S <- rbind(S, sparseMatrix(rep(1,3), c(i1, i2, ncol(S)), x=c(1, 1-params$flux.delta.rel, (params$flux.delta.rel-2)*params$flux.bound), dims=c(1,ncol(S))))
      # z-
      S <- rbind(cbind(S, sparseMatrix(NULL,NULL,dims=c(nrow(S),1))),
                 sparseMatrix(rep(1,3), c(i1, i2, ncol(S)+1), x=c(1, 1-params$flux.delta.rel, (2-params$flux.delta.rel)*params$flux.bound), dims=c(1,ncol(S)+1)))
      S <- rbind(S, sparseMatrix(rep(1,3), c(i1, i2, ncol(S)), x=c(1, params$flux.delta.rel-1, (params$flux.delta.rel-2)*params$flux.bound), dims=c(1,ncol(S))))
    } else if (df<0) {
      # additional z+
      S <- rbind(S, sparseMatrix(rep(1,3), c(i1, i2, ncol(S)), x=c(1-params$flux.delta.rel, 1, (params$flux.delta.rel-2)*params$flux.bound), dims=c(1,ncol(S))))
      # z-
      S <- rbind(cbind(S, sparseMatrix(NULL,NULL,dims=c(nrow(S),1))),
                 sparseMatrix(rep(1,3), c(i1, i2, ncol(S)+1), x=c(1-params$flux.delta.rel, 1, (2-params$flux.delta.rel)*params$flux.bound), dims=c(1,ncol(S)+1)))
      S <- rbind(S, sparseMatrix(rep(1,3), c(i1, i2, ncol(S)), x=c(params$flux.delta.rel-1, 1, (params$flux.delta.rel-2)*params$flux.bound), dims=c(1,ncol(S))))
    }
    if (df!=0) { # for now I haven't implemented df==0, so when df==0 should do nothing
      # (z+) + (z-) = 1
      S <- rbind(S, sparseMatrix(rep(1,2), c(ncol(S)-1, ncol(S)), dims=c(1,ncol(S))))
      model$rowlb <- c(model$rowlb, c((params$flux.delta.rel-2)*params$flux.bound+params$flux.delta, (params$flux.delta.rel-2)*params$flux.bound, (params$flux.delta.rel-2)*params$flux.bound+params$flux.delta), 0)
      model$rowub <- c(model$rowub, c((2-params$flux.delta.rel)*params$flux.bound, (2-params$flux.delta.rel)*params$flux.bound-params$flux.delta, (2-params$flux.delta.rel)*params$flux.bound), 1)
      model$lb <- c(model$lb, 0)
      model$ub <- c(model$ub, 1)
      model$c <- c(model$c, 1)
      model$vtype <- c(model$vtype, "I")
      model$var.ind <- c(model$var.ind, "z-")
    }
  }
  model$S <- S
  model
}

form.imat.xde <- function(model, dflux, params) {
  # model can be a "raw" model, or can be output from form.imat()
  # model will be converted to an environment if it is not, and will be modified in place by form.imat.xde0()
  if (!is.environment(model)) model <- as.environment(model)
  irxns <- model$irxn.ids
  if (is.list(dflux)) dflux <- lapply(dflux, function(x) x[irxns]) else dflux <- dflux[irxns]
  n <- length(irxns)
  nc <- sum(model$vtype=="C") - n
  
  if (!is.list(dflux)) {
    for (i in 1:n) {
      rr <- model$lb[irxns[i]]<0
      i1 <- irxns[i]
      i2 <- nc + i
      df <- dflux[i]
      form.imat.xde0(model, i1, i2, df, rr, params)
    }
  }
  if (is.list(dflux) && length(dflux)==3) {
    for (i in 1:n) {
      rr <- model$lb[irxns[i]]<0
      # 1--2
      i1 <- irxns[i]
      i2 <- nc - n + i
      df <- dflux$de12[i]
      form.imat.xde0(model, i1, i2, df, rr, params)
      # 2--3
      i1 <- nc - n + i
      i2 <- nc + i
      df <- dflux$de23[i]
      form.imat.xde0(model, i1, i2, df, rr, params)
      # 1--3
      i1 <- irxns[i]
      i2 <- nc + i
      df <- dflux$de13[i]
      form.imat.xde0(model, i1, i2, df, rr, params)
    }
  }
  model
}

form.imatx <- function(model, expr, dflux, params) {
  
  # original imat model
  imat.model <- form.imat(model, expr, params)
  
  # DE model
  form.imat.xde(imat.model, dflux, params) # modify imat.model in place
  
  imat.model
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
  if ("n" %in% names(params)) {
    n <- params$n
    params$n <- NULL
  } else n <- 1
  
  res <- Rcplex(cvec=cvec, objsense=objsense, Amat=Amat, bvec=bvec, sense=sense, lb=lb, ub=ub, vtype=vtype, control=params, n=n)

  if (!is.na(n) && n==1) res <- list(res)
  res <- lapply(res, function(x) x[c("xopt","obj","status")])
  stats <- sapply(res, function(x) x$status)
  stats <- unique(stats[!stats %in% c(101,102,128,129,130)])

  if (length(stats)>0) {
    model$milp.out <- res
    stop("iMAT: Potential problems running MILP. Please check before going on. Solver status: ", paste(stats, collapse=", "), ".\n")
  }

  model$milp.out <- res

  model
}

update.model.imat <- function(model, imat.res, sol=1, params) {
  
  yp <- imat.res$milp.out[[sol]]$xopt[imat.res$var.ind=="y+"]
  ym <- imat.res$milp.out[[sol]]$xopt[imat.res$var.ind=="y-"]
  y0 <- imat.res$milp.out[[sol]]$xopt[imat.res$var.ind=="y0"]
  
  fw <- imat.res$rxns.act[yp==1]
  bk <- imat.res$rxns.act.rev[ym==1]
  inact <- imat.res$rxns.inact[y0==1]
  inact.rev <- intersect(inact, imat.res$rxns.inact.rev)
  
  # update model
  model$lb[fw] <- params$flux.act
  model$ub[bk] <- -params$flux.act
  model$ub[inact] <- params$flux.inact
  model$lb[inact.rev] <- -params$flux.inact
  model$c <- rep(0, length(model$c))
  model$vtype <- rep("C", length(model$c))
  model$var.ind <- rep("v", length(model$c))
  model
}

update.model.imat.xde <- function(imat.res, sol=1, params) {
  
  # extract the model from imat.res
  # rows and cols to keep from imat.res$S: the "v" part and the z==1 part.
  tmp <- imat.res$S[, (imat.res$var.ind %in% c("z+","z-") & imat.res$milp.out[[sol]]$xopt==0) | imat.res$var.ind %in% c("y+","y-","y0")]
  rind <- rowSums(tmp)==0
  cind <- imat.res$var.ind=="v"
  x <- rowSums(imat.res$S[rind, imat.res$var.ind %in% c("z+","z-") & imat.res$milp.out[[sol]]$xopt==1])
  res <- subset.model(imat.res, rind, cind)
  res$rowlb <- res$rowlb - x
  res$rowub <- res$rowub - x

  res
}

update.model.imatx <- function(imat.res, sol=1, params) {
  
  ## the part for imat.de: extract the model from imat.res
  res <- update.model.imat.xde(imat.res, sol, params)
  
  ## the part for the original iMAT
  res <- update.model.imat(res, imat.res, sol, params)
}

