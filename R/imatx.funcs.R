library(Matrix)
library(data.table)
library(Rcplex.my)
# need to source utils.R
# need to source sampling.funcs.R

imat.pars <- list(flux.act=1, flux.inact=0.1, flux.delta.rel=0, flux.delta=0.1, flux.bound=1000)
milp.pars <- list(trace=1, nodesel=0, solnpoolagap=0, solnpoolgap=0, solnpoolintensity=2, n=1)
sampl.pars <- list(n.warmup=5000, n.burnin=1000, n.sampl=2000, steps.per.pnt=400, ncores=1L)

imat <- function(model, expr, imat.params=imat.pars, milp.params=milp.pars, sampl.params=sampl.pars) {
  
  # model as environment
  model <- as.environment(model)

  # formulate iMAT model
  imat.model <- form.imat(model, expr, imat.params)

  # run the iMAT MILP
  run.imat(imat.model, milp.params) # modify imat.model in place
  
  # update the original metabolic model based on iMAT result
  update.model(model, imat.model, sol=1, imat.params) # update model in place

  if (!is.null(sampl.params)) {
    # sample the metabolic model to get the fluxes of the reference state
    sample.model(model, sampl.params) # update model in place
  }
  
  # close CPLEX
  Rcplex.close()

  # return
  model
}

form.imatx <- function(model, expr, params) {

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
  as.environment(list(irxn.ids=model$irxn.ids, genes.int=expr, rxns.int.raw=rxns.int.raw, rxns.int=rxns.int,
                      rxns.act=rxns.act, rxns.act.rev=rxns.act.rev, rxns.inact=rxns.inact, rxns.inact.rev=rxns.inact.rev, var.ind=var.ind,
                      c=c, S=S, rowlb=rowlb, rowub=rowub, lb=lb, ub=ub, vtype=vtype))
}

form.imatx.de2 <- function(model, dflux, params) {
  
  irxns <- model$irxn.ids
  dflux <- dflux[irxns]
  nc2 <- length(irxns)
  nc1 <- ncol(model$S) - nc2
  
  # 1. Changed reactions in cell2 vs cell1: specify the z+ indicator variables representing v1 <= (1-delta)*v2 when rxns.int>0 (i.e. upregulated in 2 vs 1) and v2 <= (1-delta)*v1 when rxns.int<0 (i.e. downregulated in 2 vs 1)
  rxns.d <- which(dflux!=0)
  n.d <- length(rxns.d)
  m1 <- sparseMatrix(1:n.d, irxns[rxns.d], x=ifelse(dflux[rxns.d]>0, 1, params$flux.delta.rel-1), dims=c(n.d, nc1))
  m2 <- sparseMatrix(1:n.d, rxns.d, x=ifelse(dflux[rxns.d]<0, 1, params$flux.delta.rel-1), dims=c(n.d, nc2))
  m3 <- Diagonal(n.d, x=(2-params$flux.delta.rel)*params$flux.bound)
  S <- cbind(m1, m2, m3)
  
  # 1.1 For those reversible ones among #1:
  rxns.drev <- which(dflux!=0 & model1$lb[irxns]<0)
  n.drev <- length(rxns.drev)
  # 1.1a. need to add additional constraint involving z+, representing v1 >= -(1-delta)*v2 (when 2>1) or v2 >= -(1-delta)*v1 (when 2<1)
  m1 <- sparseMatrix(1:n.drev, irxns[rxns.drev], x=ifelse(dflux[rxns.drev]>0, 1, 1-params$flux.delta.rel), dims=c(n.drev, nc1))
  m2 <- sparseMatrix(1:n.drev, rxns.drev, x=ifelse(dflux[rxns.drev]<0, 1, 1-params$flux.delta.rel), dims=c(n.drev, nc2))
  m3 <- sparseMatrix(1:n.drev, match(rxns.drev, rxns.d), x=(params$flux.delta.rel-2)*params$flux.bound, dims=c(n.drev, n.d))
  S <- rbind(S, cbind(m1, m2, m3))
  # 1.1b. need to specify the extra z- indicator variables, representing the alternative case of (1-delta)*v2 <= v1 <= -(1-delta)*v2 (when 2>1) or (1-delta)*v1 <= v2 <= -(1-delta)*v1 (when 2<1)
  # thus, an reversible active reaction has both the z+ and z- indicator variables, at most one of z+ and z- can be 1
  m1 <- sparseMatrix(1:n.drev, irxns[rxns.drev], x=ifelse(dflux[rxns.drev]>0, 1, 1-params$flux.delta.rel), dims=c(n.drev, nc1))
  m2 <- sparseMatrix(1:n.drev, rxns.drev, x=ifelse(dflux[rxns.drev]<0, 1, 1-params$flux.delta.rel), dims=c(n.drev, nc2))
  m3 <- sparseMatrix(NULL, NULL, dims=c(n.drev, n.d))
  m4 <- Diagonal(n.drev, x=(2-params$flux.delta.rel)*params$flux.bound)
  m1a <- sparseMatrix(1:n.drev, irxns[rxns.drev], x=ifelse(dflux[rxns.drev]>0, 1, params$flux.delta.rel-1), dims=c(n.drev, nc1))
  m2a <- sparseMatrix(1:n.drev, rxns.drev, x=ifelse(dflux[rxns.drev]<0, 1, params$flux.delta.rel-1), dims=c(n.drev, nc2))
  S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(nrow(S), n.drev))), cbind(m1,m2,m3,m4), cbind(m1a,m2a,m3,-m4))
  # 1.1c. forcing (z+) + (z-) <=1 for the reversible reactions, otherwise there can be cases where both z+ and z- = 1 due to numerical errors (?)
  m1 <- sparseMatrix(NULL, NULL, dims=c(n.drev, nc1+nc2))
  m2 <- sparseMatrix(1:n.drev, match(rxns.drev, rxns.d), dims=c(n.drev, n.d))
  m3 <- Diagonal(n.drev)
  S <- rbind(S, cbind(m1, m2, m3))
  
  # the final S matrix
  S <- rbind(cbind(model$S, sparseMatrix(NULL, NULL, dims=c(nrow(model$S), ncol(S)-ncol(model$S)))),
             S)
  
  # other parameters
  rowlb <- c(model$rowlb,
             rep((params$flux.delta.rel-2)*params$flux.bound, n.d),
             rep((params$flux.delta.rel-2)*params$flux.bound+params$flux.delta, n.drev),
             rep((params$flux.delta.rel-2)*params$flux.bound, n.drev),
             rep((params$flux.delta.rel-2)*params$flux.bound+params$flux.delta, n.drev),
             rep(0,n.drev))
  rowub <- c(model$rowub,
             rep((2-params$flux.delta.rel)*params$flux.bound-params$flux.delta, n.d),
             rep((2-params$flux.delta.rel)*params$flux.bound, n.drev),
             rep((2-params$flux.delta.rel)*params$flux.bound-params$flux.delta, n.drev),
             rep((2-params$flux.delta.rel)*params$flux.bound, n.drev),
             rep(1,n.drev))
  n <- ncol(S) - ncol(model$S)
  lb <- c(model$lb, rep(0, n))
  ub <- c(model$ub, rep(1, n))
  c <- c(model$c, rep(1, n))
  vtype <- c(model$vtype, rep("I", n))
  var.ind <- c(model$var.ind, rep(c("z+","z-"), c(n.d, n.drev)))
  res <- list(dflux=dflux, rxns.diff=rxns.d, rxns.diff.rev=rxns.drev,
              rxns.act=model$rxns.act, rxns.act.rev=model$rxns.act.rev,
              rxns.inact=model$rxns.inact, rxns.inact.rev=model$rxns.inact.rev, var.ind=var.ind,
              c=c, S=S, rowlb=rowlb, rowub=rowub, lb=lb, ub=ub, vtype=vtype)
  
  # return iMAT model as an environment
  as.environment(res)
}

form.imatx2 <- function(model, expr, dflux, params) {
  
  # original imat model for expr
  model <- form.imatx(model, expr, params)

  # DE model
  form.imatx.de2(model, dflux, params)
}

form.imatx3 <- function(model, expr, df12, df23, df13, params) {
  
  # original imat model for expr
  model <- form.imatx(model, expr, params)
  
  # DE model
  form.imatx.de3(model, df12, df23, df13, params)
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

  if (n==1) res <- list(res)
  res <- lapply(res, function(x) x[c("xopt","obj","status")])
  stats <- sapply(res, function(x) x$status)
  stats <- unique(stats[!stats %in% c(101,102,129,130)])

  if (length(stats)>0) {
    model$milp.out <- res
    stop("iMAT: Potential problems running MILP. Please check before going on. Solver status: ", paste(stats, collapse=", "), ".\n")
  }

  model$milp.out <- res

  model
}

update.model <- function(model, imat.res, sol=1, params) {
  
  if ("dflux" %in% ls(imat.res)) {
    # join two `model` to create the bi-model here
  }
  ## the part for the original iMAT
  yp <- imat.res$milp.out[[sol]]$xopt[imat.res$var.ind %in% c("y+","y+_1","y+_2")]
  ym <- imat.res$milp.out[[sol]]$xopt[imat.res$var.ind %in% c("y-","y-_1","y-_2")]
  y0 <- imat.res$milp.out[[sol]]$xopt[imat.res$var.ind %in% c("y0","y0_1","y0_2")]
  
  fw <- imat.res$rxns.act[yp==1]
  bk <- imat.res$rxns.act.rev[ym==1]
  inact <- imat.res$rxns.inact[y0==1]
  inact.rev <- intersect(inact, imat.res$rxns.inact.rev)
  
  # update model
  model$lb[fw] <- params$flux.act
  model$ub[bk] <- -params$flux.act
  model$ub[inact] <- params$flux.inact
  model$lb[inact.rev] <- -params$flux.inact
  
  ## the part for imat.de, if applicable
  if ("dflux" %in% ls(imat.res)) {
    dflx <- imat.res$dflux[imat.res$rxns.diff]
    zp <- imat.res$milp.out[[sol]]$xopt[imat.res$var.ind=="z+"]
    up.zp <- imat.res$rxns.diff[zp==1 & dflx>0]
    n.up.zp <- length(up.zp)
    up.zp.rev <- intersect(up.zp, imat.res$rxns.diff.rev)
    n.up.zp.rev <- length(up.zp.rev)
    dn.zp <- imat.res$rxns.diff[zp==1 & dflx<0]
    n.dn.zp <- length(dn.zp)
    dn.zp.rev <- intersect(dn.zp, imat.res$rxns.diff.rev)
    n.dn.zp.rev <- length(dn.zp.rev)
    
    dflx <- imat.res$dflux[imat.res$rxns.diff.rev]
    zm <- imat.res$milp.out[[sol]]$xopt[imat.res$var.ind=="z-"]
    up.zm <- imat.res$rxns.diff.rev[zm==1 & dflx>0]
    n.up.zm <- length(up.zm)
    dn.zm <- imat.res$rxns.diff.rev[zm==1 & dflx<0]
    n.dn.zm <- length(dn.zm)
    
    # update model
    nc1
    nc2
    S <- rind(
      cbind(sparseMatrix(1:n.up.zp, up.zp, dims=c(n.up.zp, nc1)), sparseMatrix(1:n.up.zp, up.zp, x=params$flux.delta.rel-1, dims=c(n.up.zp, nc2))),
      cbind(sparseMatrix(1:n.dn.zp, dn.zp, x=params$flux.delta.rel-1, dims=c(n.dn.zp, nc1)), sparseMatrix(1:n.dn.zp, dn.zp, dims=c(n.up.zp, nc2))),
      cbind(sparseMatrix(1:n.up.zp.rev, up.zp.rev, dims=c(n.up.zp.rev, nc1)), sparseMatrix(1:n.up.zp.rev, up.zp.rev, x=1-params$flux.delta.rel, dims=c(n.up.zp.rev, nc2))),
      cbind(sparseMatrix(1:n.dn.zp.rev, dn.zp.rev, x=1-params$flux.delta.rel, dims=c(n.dn.zp.rev, nc1)), sparseMatrix(1:n.dn.zp.rev, dn.zp.rev, dims=c(n.dn.zp.rev, nc2))),
      cbind(sparseMatrix(1:n.up.zm, up.zm, dims=c(n.up.zm, nc1)), sparseMatrix(1:n.up.zm, up.zm, x=1-params$flux.delta.rel, dims=c(n.up.zm, nc2))),
      cbind(sparseMatrix(1:n.dn.zm, dn.zm, x=1-params$flux.delta.rel, dims=c(n.dn.zm, nc1)), sparseMatrix(1:n.dn.zm, dn.zm, dims=c(n.dn.zm, nc2))),
      cbind(sparseMatrix(1:n.up.zm, up.zm, dims=c(n.up.zm, nc1)), sparseMatrix(1:n.up.zm, up.zm, x=params$flux.delta.rel-1, dims=c(n.up.zm, nc2))),
      cbind(sparseMatrix(1:n.dn.zm, dn.zm, x=params$flux.delta.rel-1, dims=c(n.dn.zm, nc1)), sparseMatrix(1:n.dn.zm, dn.zm, dims=c(n.dn.zm, nc2)))
    )
    model$S <- rbind(model$S, S)
    model$rowlb <- c(model$rowlb,
                     rep((params$flux.delta.rel-2)*params$flux.bound, n.up.zp+n.dn.zp),
                     rep(params$flux.delta, n.up.zp.rev+n.dn.zp.rev),
                     rep((params$flux.delta.rel-2)*params$flux.bound, n.up.zm+n.dn.zm),
                     rep(params$flux.delta, n.up.zm+n.dn.zm))
    model$rowub <- c(model$rowub,
                     rep(-params$flux.delta, n.up.zp+n.dn.zp),
                     rep((2-params$flux.delta.rel)*params$flux.bound, n.up.zp.rev+n.dn.zp.rev),
                     rep(-params$flux.delta, n.up.zm+n.dn.zm),
                     rep((2-params$flux.delta.rel)*params$flux.bound, n.up.zm+n.dn.zm))
  }
  
  model
}

