library(Matrix)
library(data.table)
library(Rcplex.my)
# need to source utils.R
# need to source sampling.funcs.R


##### ----- parameters ----- #####

imat.pars <- list(flux.act=1, flux.inact=0.1, flux.delta.rel=0, flux.delta=0.1, flux.bound=1000)
milp.pars <- list(trace=1, nodesel=0, solnpoolagap=0, solnpoolgap=0, solnpoolintensity=2, n=1)
sampl.pars <- list(n.warmup=5000, n.burnin=1000, n.sampl=2000, steps.per.pnt=400, ncores=1L)
mep.pars <- list(beta=1e9, damp=0.9, max.iter=2000, dlb=1e-50, dub=1e50, epsil=1e-6, fix.flux=FALSE, fflux.id=0, fflux.mean=0, fflux.var=0)


##### ----- different versions of iMAT ----- #####

imat <- function(model, expr, imat.params=imat.pars, milp.params=milp.pars, sampl.params=sampl.pars) {
  
  # formulate iMAT model
  imat.model <- form.imat(model, expr, imat.params)

  # run the iMAT MILP
  run.imat(imat.model, milp.params) # modify imat.model in place
  
  # update the original metabolic model based on iMAT result
  res.model <- update.model.imat(model, imat.model, sol=0, imat.params)

  if (!is.null(sampl.params)) {
    # sample the metabolic model to get the fluxes of the reference state
    sample.model(res.model, sampl.params) # update model in place
  }
  
  # close CPLEX
  Rcplex.close()

  # return
  res.model
}

imat.mep <- function(model, expr, imat.params=imat.pars, milp.params=milp.pars, mep.params=mep.pars, nc=1L) {
  
  # formulate iMAT model, return an environment
  imat.model <- form.imat(model, expr, imat.params)

  # run the iMAT MILP
  run.imat(imat.model, milp.params) # imat.model modified in place since it's an environment

  # update the original metabolic model based on iMAT result
  res.model <- update.model.imat(model, imat.model, sol=0, imat.params)

  # process model for MEP
  cat("Preparing for MEP...\n")
  preprocess.model(res.model, nc=nc) # res.model modified in place since it's an environment

  # close CPLEX
  Rcplex.close()

  # run MEP to get the fluxes of the reference state
  run.mep(res.model, mep.params) # results added to res.model in place since it's an environment
  
  # return
  res.model
}

imatx <- function(model, expr, dflux, imat.params=imat.pars, milp.params=milp.pars, sampl.params=sampl.pars) {

  # formulate iMAT model
  imat.model <- form.imat(model, expr, imat.params)
  imat.model <- form.imat.xde(imat.model, dflux, imat.params)
  
  # run the iMAT MILP
  run.imat(imat.model, milp.params) # modify imat.model in place
  
  # update the original metabolic model based on iMAT result
  res.model <- update.model.imat(model, imat.model, sol=0, imat.params)
  
  if (!is.null(sampl.params)) {
    # sample the metabolic model to get the fluxes of the reference state
    sample.model(res.model, sampl.params) # update res.model in place
  }
  
  # close CPLEX
  Rcplex.close()
  
  # return
  res.model
}

imatx2steps <- function(model, expr, dflux, imat.params=imat.pars, milp.params1=milp.pars, milp.params2=milp.pars, sampl.params=sampl.pars) {
  
  # formulate iMAT model step 1 (i.e. the original imat)
  imat.model1 <- form.imat(model, expr, imat.params)
  # run the iMAT MILP
  run.imat(imat.model1, milp.params1) # modify imat.model in place
  
  # update the original metabolic model based on iMAT result
  imat.model2 <- update.model.imat(model, imat.model1, sol=0, imat.params)
  # formulate iMAT model step 2 (xde) upon the updated model from the above step
  imat.model2 <- form.imat.xde(imat.model2, dflux, imat.params)
  
  # run the iMAT MILP
  run.imat(imat.model2, milp.params2) # modify imat.model in place
  
  # update the original metabolic model based on iMAT result
  res.model <- update.model.imat(model, imat.model2, sol=0, imat.params)
  
  if (!is.null(sampl.params)) {
    # sample the metabolic model to get the fluxes of the reference state
    sample.model(res.model, sampl.params) # update model in place
  }
  
  # close CPLEX
  Rcplex.close()
  
  # return
  res.model
}

imat2 <- function(model, expr1, expr2, dflux, imat.params=imat.pars, milp.params=milp.pars, sampl.params=sampl.pars) {
  
  # formulate iMAT model
  # original imat model for expr1
  model1 <- form.imat(model, expr1, imat.params)
  # original imat model for expr2
  model2 <- form.imat(model, expr2, imat.params)
  # DE model
  imat.model <- form.imat.de(model1, model2, dflux, imat.params)

  # run the iMAT MILP
  run.imat(imat.model, milp.params) # modify imat.model in place
  
  # update the original metabolic model based on iMAT result
  res.model <- update.model.imat(model, imat.model, sol=0, imat.params)
  
  if (!is.null(sampl.params)) {
    # sample the metabolic model to get the fluxes of the reference state
    sample.model(res.model, sampl.params) # update res.model in place
  }
  
  # close CPLEX
  Rcplex.close()
  
  # return
  res.model
}

imat2steps <- function(model, expr1, expr2, dflux, imat.params=imat.pars, milp.params1=milp.pars, milp.params2=milp.pars, sampl.params=sampl.pars) {
  
  # formulate iMAT model step 1 (i.e. the original imat)
  imat.model1 <- form.imat(model, expr1, imat.params)
  imat.model2 <- form.imat(model, expr2, imat.params)
  # run the iMAT MILP
  run.imat(imat.model1, milp.params1) # modify imat.model in place
  run.imat(imat.model2, milp.params1) # modify imat.model in place
  # update the original metabolic model based on iMAT result
  model1 <- update.model.imat(model, imat.model1, sol=0, imat.params)
  model2 <- update.model.imat(model, imat.model2, sol=0, imat.params)

  # formulate iMAT model step 2 (de) upon the updated model from the above step
  imat.model.de <- form.imat.de(model1, model2, dflux, imat.params)
  # run the iMAT MILP
  run.imat(imat.model.de, milp.params2) # modify imat.model in place
  
  # update model (model1+model2) based on iMAT result
  res.model <- update.model.imat(model, imat.model.de, sol=0, imat.params) # here model is actually not used
  
  if (!is.null(sampl.params)) {
    # sample the metabolic model to get the fluxes of the reference state
    sample.model(res.model, sampl.params) # update model in place
  }
  
  # close CPLEX
  Rcplex.close()
  
  # return
  res.model
}


##### ----- helper functions for the individual internal steps of iMAT ----- #####

form.imat <- function(model, expr, params) {
  # formulate the original iMAT model
  
  # reaction data
  rxns.int.raw <- exprs2rxns(expr, 0, model)
  rxns.int <- rxns.int.raw
  # remove integers for dead end rxns
  rxns.int[model$lb==0 & model$ub==0] <- 0
  
  n.mets <- nrow(model$S)
  n.rxns <- ncol(model$S)
  S <- model$S
  
  # 1. Active reactions: specify the y+ indicator variables, representing activation in the forward direction (i.e. v>flux.act)
  rxns.act <- which(rxns.int==1)
  n.act <- length(rxns.act)
  if (n.act!=0) {
    m1 <- sparseMatrix(1:n.act, rxns.act, dims=c(n.act, n.rxns))
    m2 <- Diagonal(n.act, x=(-params$flux.act-params$flux.bound))
    S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(n.mets, n.act))), cbind(m1, m2))
  }
  
  # 2. Reversible active reactions: for those reversible ones among the active reactions, specify the extra y- indicator variables, representing activation in the backward direction (i.e. v<-flux.act)
  # thus, an reversible active reaction has both the y+ and y- indicator variables, because it can be active in either direction (but never both, i.e. 1 XOR 2)
  rxns.act.rev <- which(rxns.int==1 & model$lb<0)
  n.act.rev <- length(rxns.act.rev)
  if (n.act.rev!=0) {
    m1 <- sparseMatrix(1:n.act.rev, rxns.act.rev, dims=c(n.act.rev, ncol(S)))
    m2 <- Diagonal(n.act.rev, x=params$flux.act+params$flux.bound)
    S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(nrow(S), n.act.rev))), cbind(m1, m2))
  }
  
  # 3. Inactive reactions: specify the y0 indicator variables
  # 3a. specify inactivation in the forward direction (i.e. v<flux.inact)
  rxns.inact <- which(rxns.int==-1)
  n.inact <- length(rxns.inact)
  if (n.inact!=0) {
    m1 <- sparseMatrix(1:n.inact, rxns.inact, dims=c(n.inact, ncol(S)))
    m2 <- Diagonal(n.inact, x=params$flux.bound-params$flux.inact)
    S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(nrow(S), n.inact))), cbind(m1, m2))
  }
  # 3b. for those reversible inactive reactions, need to further specify inactivation in the backward direction (i.e. v>-flux.inact)
  # note that a reversible inactive reaction has only one y0 indicator variable, because for these reactions we want -flux.inact<v<flux.inact (3a AND 3b) 
  rxns.inact.rev <- which(rxns.int==-1 & model$lb<0)
  n.inact.rev <- length(rxns.inact.rev)
  if (n.inact.rev!=0) {
    m3 <- sparseMatrix(1:n.inact.rev, rxns.inact.rev, dims=c(n.inact.rev, ncol(S)-n.inact))
    m4 <- sparseMatrix(1:n.inact.rev, match(rxns.inact.rev, rxns.inact), x=params$flux.inact-params$flux.bound, dims=c(n.inact.rev, n.inact))
    S <- rbind(S, cbind(m3, m4))
  }
  
  # other parameters
  n <- n.act + n.act.rev + n.inact + n.inact.rev
  rowlb <- c(model$rowlb, rep(-params$flux.bound, n))
  rowub <- c(model$rowub, rep(params$flux.bound, n))
  n <- ncol(S) - n.rxns
  c <- rep(c(0,1), c(n.rxns, n))
  vtype <- ifelse(c==1, "I", "C")
  lb <- c(model$lb, rep(0, n))
  ub <- c(model$ub, rep(1, n))
  var.ind <- rep(c("v","y+","y-","y0"), c(n.rxns, n.act, n.act.rev, n.inact)) # iMAT variable type indicators (v: fluxex; y+/-/0: indicator variables)
  
  # return iMAT model as an environment; set class as "imat.model" so later functions can identify imat models based on this (but currently not in the canonical S3 OOP way)
  res.model <- as.environment(list(irxn.ids=model$irxn.ids, # if not exist, will be NULL
                      genes.int=expr, rxns.int=rxns.int,
                      rxns.act=rxns.act, rxns.act.rev=rxns.act.rev, rxns.inact=rxns.inact, rxns.inact.rev=rxns.inact.rev, var.ind=var.ind,
                      c=c, S=S, rowlb=rowlb, rowub=rowub, lb=lb, ub=ub, vtype=vtype))
  class(res.model) <- "imat.model"
  res.model
}

form.imat.de0 <- function(model, i1, i2, df, rr, params) {
  # a helper function to formulate the iMAT DE model

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
      model$rowlb <- c(model$rowlb,
                       c((params$flux.delta.rel-2)*params$flux.bound+params$flux.delta,
                         (params$flux.delta.rel-2)*params$flux.bound,
                         (params$flux.delta.rel-2)*params$flux.bound+params$flux.delta),
                       0)
      model$rowub <- c(model$rowub,
                       c((2-params$flux.delta.rel)*params$flux.bound,
                         (2-params$flux.delta.rel)*params$flux.bound-params$flux.delta,
                         (2-params$flux.delta.rel)*params$flux.bound),
                       1)
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
  # formulate the iMAT DE model for the DE among cells for a multi-cellular model
  
  # model can be a "raw" or "updated" model (i.e. from update.model.imat), or can be an imat model (from form.imat)
  # if model is not an output from form.imat, then need to initiate a few things
  res.model <- as.environment(model)
  if (!"imat.model" %in% class(model)) {
    res.model$c <- rep(0, ncol(model$S))
    res.model$vtype <- rep("C", ncol(model$S))
    res.model$var.ind <- rep("v", ncol(model$S))
  }
  irxns <- res.model$irxn.ids
  if (is.list(dflux)) dflux <- lapply(dflux, function(x) x[irxns]) else dflux <- dflux[irxns]
  n <- length(irxns)
  nc <- sum(res.model$vtype=="C") - n
  
  if (!is.list(dflux)) {
    for (i in 1:n) {
      rr <- model$lb[irxns[i]]<0
      i1 <- irxns[i]
      i2 <- nc + i
      df <- dflux[i]
      form.imat.de0(res.model, i1, i2, df, rr, params) # modify res.model in place
    }
  }
  if (is.list(dflux) && length(dflux)==3) {
    for (i in 1:n) {
      rr <- model$lb[irxns[i]]<0
      # 1--2
      i1 <- irxns[i]
      i2 <- nc - n + i
      df <- dflux$de12[i]
      form.imat.de0(res.model, i1, i2, df, rr, params) # modify res.model in place
      # 2--3
      i1 <- nc - n + i
      i2 <- nc + i
      df <- dflux$de23[i]
      form.imat.de0(res.model, i1, i2, df, rr, params) # modify res.model in place
      # 1--3
      i1 <- irxns[i]
      i2 <- nc + i
      df <- dflux$de13[i]
      form.imat.de0(res.model, i1, i2, df, rr, params) # modify res.model in place
    }
  }
  res.model
}

form.imat.de <- function(model1, model2=model1, dflux, params) {
  # formulate the iMAT DE model for the DE between two models (i.e. two samples)
  
  # model1 and model2 can be "raw" or "updated" models (i.e. from update.model.imat), or can be an imat model (from form.imat)
  res.model <- as.environment(c.model(model1, model2))
  # if model1 or model2 is not an output from form.imat, then need to initiate a few things; otherwise, need to fix a few things not handled properly by c.model()
  if (!"imat.model" %in% class(model1)) {
    res.model$c <- rep(0, ncol(res.model$S))
    res.model$vtype <- rep("C", ncol(res.model$S))
    res.model$var.ind <- c(rep("v_1", ncol(model1$S)), rep("v_2", ncol(model2$S)))
  } else {
    nv1 <- sum(model1$var.ind=="v")
    res.model$rxns.act <- c(model1$rxns.act, model2$rxns.act+nv1)
    res.model$rxns.act.rev <- c(model1$rxns.act.rev, model2$rxns.act.rev+nv1)
    res.model$rxns.inact <- c(model1$rxns.inact, model2$rxns.inact+nv1)
    res.model$rxns.inact.rev <- c(model1$rxns.inact.rev, model2$rxns.inact.rev+nv1)
  }
  nc1 <- ncol(model1$S)

  for (i in 1:length(dflux)) {
    rr <- model1$lb[i]<0 # use model1 is OK since even if it's updated by iMAT (first step) the reversibilities of reactions don't change
    form.imat.de0(res.model, i, nc1+i, dflux[i], rr, params)
  }
  
  res.model
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
  
  if (!is.null(names(res))) res <- list(res) # if only 1 solution, wrap it in a list
  res <- lapply(res, function(x) x[c("xopt","obj","status")])
  model$milp.out <- res
  model
}

get.imat.xopt <- function(imat.res, sol=0) {
  # if sol==0, pool all solutions and take the "consensus", otherwise use the one solution specified
  if (sol==0) {
    xopt <- do.call(cbind, lapply(imat.res$milp.out, function(x) {
      if (x$status %in% c(101,102,128,129,130)) x$xopt else NULL
    }))
    if (is.null(xopt)) stop("run.imat: all MILP solutions may contain issues.\n")
    xopt <- rowMeans(xopt)
    xopt[imat.res$vtype=="I"] <- as.numeric(xopt[imat.res$vtype=="I"]>=0.667)
  } else {
    if (!imat.res$milp.out[[sol]]$status %in% c(101,102,128,129,130)) stop("run.imat: the selected MILP solution may contain issues.\n")
    xopt <- imat.res$milp.out[[sol]]$xopt
  }
  xopt
}

update.model.imat <- function(model, imat.res, sol=0, params) {
  
  xopt <- get.imat.xopt(imat.res, sol)
  
  # the de part of the model
  if ("z+" %in% imat.res$var.ind) {
    # if imat.res is from a de model, then obtain the resulting model from imat.res
    # rows and cols to keep from imat.res$S: the "v" part and the z==1 part.
    tmp <- imat.res$S[, (imat.res$var.ind %in% c("z+","z-") & xopt==0) | imat.res$var.ind %in% c("y+","y-","y0")]
    rind <- apply(tmp, 1, function(x) all(x==0))
    cind <- imat.res$var.ind=="v"
    x <- rowSums(imat.res$S[rind, imat.res$var.ind %in% c("z+","z-") & xopt==1])
    res.model <- as.environment(subset.model(imat.res, rind, cind))
    res.model$rowlb <- res.model$rowlb - x
    res.model$rowub <- res.model$rowub - x
  } else {
    # if imat.res if not from a de model, then the resulting model will be based on model
    res.model <- as.environment(model)
    res.model$milp.out <- imat.res$milp.out # add the milp output
  }
  
  # the original imat part of the model
  if ("y+" %in% imat.res$var.ind || "y+_1" %in% imat.res$var.ind) {
    yp <- xopt[imat.res$var.ind %in% c("y+","y+_1","y+_2")]
    ym <- xopt[imat.res$var.ind %in% c("y-","y-_1","y-_2")]
    y0 <- xopt[imat.res$var.ind %in% c("y0","y0_1","y0_2")]
    
    fw <- imat.res$rxns.act[yp==1]
    bk <- imat.res$rxns.act.rev[ym==1]
    inact <- imat.res$rxns.inact[y0==1]
    inact.rev <- intersect(inact, imat.res$rxns.inact.rev)
    
    # update model
    res.model$lb[fw] <- params$flux.act
    res.model$ub[bk] <- -params$flux.act
    res.model$ub[inact] <- params$flux.inact
    res.model$lb[inact.rev] <- -params$flux.inact
  }
  res.model
}


