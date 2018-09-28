library(Matrix)
library(Rcplex)
library(data.table)
library(parallel)

mta.params <- list(v.min=-50, v.max=50, v.min.c=-1000, v.max.c=1000, alpha=0.9, epsil=0.01)
miqp.params <- list(trace=1, tilim=120)
 
mta <- function(model, v.ref, dflux, del, mta.params=mta.params, miqp.params=miqp.params) {
  
  # formulate MTA model
  mta.model <- form.mta(model, v.ref, dflux, mta.params)

  # run the MTA MIQP
  run.mta(mta.model, del, miqp.params)
}


form.mta <- function(model, v.ref, dflux, params) {
  
  n.mets <- nrow(model$S)
  n.rxns <- ncol(model$S)

  # constraint matrix
  ## for reactions to be changed in the forward direction
  rxns.fw <- which(v.ref>=0 & dflux==1 | v.ref<0 & dflux==-1)
  n.fw <- length(rxns.fw)
  m1 <- sparseMatrix(1:n.fw, rxns.fw, dims=c(2*n.fw, n.rxns))
  m2 <- rbind(cbind(Diagonal(n.fw, x=(-v.ref[rxns.fw]-params$epsil)), Diagonal(n.fw, x=(-params$v.min))),
              cbind(Diagonal(n.fw), Diagonal(n.fw)))
  S <- rbind(cbind(model$S, sparseMatrix(NULL, NULL, dims=c(n.mets, 2*n.fw))), cbind(m1, m2))
  ## for reactions to be changed in the backward direction
  rxns.bk <- which((v.ref<=0 & dflux==1 | v.ref>0 & dflux==-1) & !(v.ref==0 & model$lb==0))
  n.bk <- length(rxns.bk)
  m1 <- sparseMatrix(1:n.bk, rxns.bk, dims=c(2*n.bk, n.rxns+2*n.fw))
  m2 <- rbind(cbind(Diagonal(n.bk, x=(-v.ref[rxns.bk]+params$epsil)), Diagonal(n.bk, x=(-params$v.max))),
              cbind(Diagonal(n.bk), Diagonal(n.bk)))
  S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(n.mets+2*n.fw, 2*n.bk))), cbind(m1, m2))
  
  # constraints
  rowlb <- c(model$rowlb, rep(0:1, each=n.fw), rep(c(params$v.min.c,1), each=n.bk))
  rowub <- c(model$rowub, rep(c(params$v.max.c,1), each=n.fw), rep(0:1, each=n.bk))
  lb <- c(model$lb, rep(0, 2*n.fw+2*n.bk))
  ub <- c(model$ub, rep(1, 2*n.fw+2*n.bk))

  # objective function and others
  n <- ncol(S)
  rxns.st <- which(dflux==0 & model$c!=1)
  tmp <- rep(0, n)
  tmp[rxns.st] <- 2*(1-params$alpha)
  F <- Diagonal(x=vec)
  c <- rep(c(0, params$alpha/2, 0, params$alpha/2), c(n.rxns+n.fw, n.fw, n.bk, n.bk))
  c[rxns.st] <- -2*(1-params$alpha)*v.ref[rxns.st]
  vtype <- rep(c("C","I"), c(n.rxns, n-n.rxns))

  # return MTA model
  list(v.ref=v.ref, dflux=dflux,
       rxns.fw=rxns.fw, rxns.bk=rxns.bk, rxns.st=rxns.st,
       c=c, F=F, S=S, rowlb=rowlb, rowub=rowub, lb=lb, ub=ub, vtype=vtype))
}

run.mta <- function(model, del, params) {
  
  names(del) <- del
  res <- lapply(del, function(i) {
    miqp.res <- run.miqp(model, i, params)
    analyz.mta.res(model, miqp.res)
  })
  rbindlist(res, idcol="del.rxn")
}

run.miqp <- function(model, del, params) {

  cvec <- model$c
  Qmat <- model$F
  objsense <- "min"
  Amat <- rbind(model$S, model$S)
  bvec <- c(model$rowlb, model$rowub)
  sense <- rep(c("G","L"), c(length(model$rowlb), length(model$rowub)))
  lb <- model$lb
  lb[del] <- 0
  ub <- model$ub
  ub[del] <- 0
  vtype <- model$vtype
  
  res <- Rcplex(cvec=cvec, Qmat=Qmat, objsense=objsense, Amat=Amat, bvec=bvec, sense=sense, lb=lb, ub=ub, vtype=vtype, control=params, n=1)
  if (res$status!=101) warning("MTA: Potential problem running MIQP. Solver status: ", res$status, ".\n")

  res
}

analyz.mta.res <- function(model, miqp.res) {
    
  v <- miqp.res$xopt
  v0 <- model$v.ref
  fw <- (1:length(v)) %in% setdiff(model$rxns.fw, model$rxns.bk)
  bk <- (1:length(v)) %in% setdiff(model$rxns.bk, model$rxns.fw)
  fw.or.bk <- (1:length(v)) %in% intersect(model$rxns.bk, model$rxns.fw) # these are those reactions intended to change in either direction with v.ref=0, thus these are always successful 
  
  # reactions intended to change: overdone (thus regared as failed)
  fw.overdo <- fw & v0<0 & v>(-v0)
  bk.overdo <- bk & v0>0 & v<(-v0)
  # reactions intended to change: successful
  fw.yes <- fw & v>v0 & !fw.overdo
  bk.yes <- bk & v<v0 & !bk.overdo
  # reactions intended to change: failed
  fw.no <- fw & v<v0
  bk.no <- bk & v>v0

  # calculate mta score
  # reactions intended to change
  s.ch <- sum(abs(v[fw.yes|bk.yes|fw.or.bk]-v0[fw.yes|bk.yes|fw.or.bk])) -
          sum(abs(v[fw.no|bk.no]-v0[fw.no|bk.no])) - 
          sum(abs(abs(v[fw.overdo|bk.overdo])-abs(v0[fw.overdo|bk.overdo])))
  # reactions intended to stay steady
  dv.st <- v[model$rxns.st]-v0[model$rxns.st]
  s.st <- sum(abs(dv.st))
  s <- d.ch/d.st

  # return
  data.table(solv.stat=miqp.res$status, obj.opt=miqp.res$obj, v.opt=v, rxns.change.yes=which(fw.yes|bk.yes|fw.or.bk), rxns.change.no=which(fw.no|bk.no), rxns.change.overdo=which(fw.overdo|bk.overdo), dv.rxns.steady=list(v[model$rxns.st]-v0[model$rxns.st]), score.change=s.ch, score.steady=s.st, mta.score=s)
}