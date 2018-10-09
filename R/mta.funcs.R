library(Matrix)
library(data.table)
library(Rcplex.my)
library(parallel)
# need to source utils.R

mta.pars <- list(v.min=-50, v.max=50, v.min.c=-1000, v.max.c=1000, alpha=0.9, epsil=0.01)
miqp.pars <- list(trace=0, tilim=120, threads=1)
 
mta <- function(model, v.ref, dflux, del="default", mta.params=mta.pars, miqp.params=miqp.pars) {
  
  # formulate MTA model
  mta.model <- form.mta(model, v.ref, dflux, mta.params)

  # run the MTA MIQP
  if (del=="default") del <- 0:ncol(model$S)
  res <- run.mta(mta.model, del, miqp.params)
  res[, del.rxn:=as.integer(del.rxn)]
  res[, genes:=list(rxns2genes(del.rxn, model))]
  rbind(res[del.rxn==0], res[del.rxn!=0][order(-mta.score)])
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
  F <- .sparseDiagonal(x=tmp)
  c <- rep(c(0, params$alpha/2, 0, params$alpha/2), c(n.rxns+n.fw, n.fw, n.bk, n.bk))
  c[rxns.st] <- -2*(1-params$alpha)*v.ref[rxns.st]
  vtype <- rep(c("C","I"), c(n.rxns, n-n.rxns))

  # return MTA model
  list(v.ref=v.ref, dflux=dflux,
       rxns.fw=rxns.fw, rxns.bk=rxns.bk, rxns.st=rxns.st,
       c=c, F=F, S=S, rowlb=rowlb, rowub=rowub, lb=lb, ub=ub, vtype=vtype)
}

run.mta <- function(model, del, params) {
  
  names(del) <- del
  res <- mclapply(del, function(i) {
    miqp.res <- run.miqp(model, i, params)
    analyz.mta.res(model, miqp.res)
  }, mc.cores=detectCores())

  # close CPLEX
  Rcplex.close()

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
  lb[del] <- 0 # if del==0, nothing will be changed to lb, meaning do not delete any reaction (the control)
  ub <- model$ub
  ub[del] <- 0 # if del==0, nothing will be changed to ub, meaning do not delete any reaction (the control)
  vtype <- model$vtype
  
  res <- Rcplex(cvec=cvec, Qmat=Qmat, objsense=objsense, Amat=Amat, bvec=bvec, sense=sense, lb=lb, ub=ub, vtype=vtype, control=params)
  if (!res$status %in% c(101,102)) warning("MTA: Potential problem running MIQP. Solver status: ", res$status, ".\n")

  res
}

analyz.mta.res <- function(model, miqp.res) {
  
  v0 <- model$v.ref
  n <- length(v0)
  v <- miqp.res$xopt[1:n]
  fw <- (1:n) %in% setdiff(model$rxns.fw, model$rxns.bk)
  bk <- (1:n) %in% setdiff(model$rxns.bk, model$rxns.fw)
  fw.or.bk <- (1:n) %in% intersect(model$rxns.bk, model$rxns.fw) # these are those reactions intended to change in either direction with v.ref=0
  
  # reactions intended to change: overdone (thus regared as failed)
  fw.overdo <- fw & v0<0 & v>(-v0)
  bk.overdo <- bk & v0>0 & v<(-v0)
  # reactions intended to change: successful
  fw.yes <- fw & v>v0 & !fw.overdo
  bk.yes <- bk & v<v0 & !bk.overdo
  fw.or.bk.yes <- fw.or.bk & v!=0
  # reactions intended to change: failed
  fw.no <- fw & v<v0
  bk.no <- bk & v>v0
  fw.or.bk.no <- fw.or.bk & v==0

  # calculate mta score
  # reactions intended to change
  yes <- which(fw.yes|bk.yes|fw.or.bk.yes)
  no <- which(fw.no|bk.no|fw.or.bk.no)
  overdo <- which(fw.overdo|bk.overdo)
  adv.yes <- abs(v[yes]-v0[yes])
  adv.no <- abs(v[no]-v0[no])
  adv.overdo <- abs(abs(v[overdo])-abs(v0[overdo]))
  s.ch <- sum(adv.yes)-sum(adv.no)-sum(adv.overdo)
  # reactions intended to stay steady
  adv.st <- abs(v[model$rxns.st]-v0[model$rxns.st])
  s.st <- sum(adv.st)
  s <- s.ch/s.st

  # return
  data.table(solv.stat=miqp.res$status, obj.opt=miqp.res$obj, v.opt=list(v), rxns.change.yes=list(yes), rxns.change.no=list(no), rxns.change.overdo=list(overdo), advs.change.yes=list(adv.yes), advs.change.no=list(adv.no), advs.change.overdo=list(adv.overdo), advs.steady=list(adv.st), score.change=s.ch, score.steady=s.st, mta.score=s)
}