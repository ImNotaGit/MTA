library(Matrix)
library(data.table)
library(Rcplex.my)
library(parallel)
# need to source utils.R

mta.pars <- list(v.min=-50, v.max=50, v.min.c=-1000, v.max.c=1000, alpha=0.9, epsil=0.01)
lp.pars <- list(trace=0, maxcalls=5000, tilim=120, threads=1)
 
mta <- function(model, v.ref, dflux, del="default", mta.params=mta.pars, lp.params=lp.pars) {
  
  # formulate MTA model
  mta.model <- form.mta(model, v.ref, dflux, mta.params)

  # run the MTA MIQP
  if (length(del)==1 && del=="default") del <- 0:ncol(model$S)
  res <- run.mta(mta.model, del, lp.params, ncores=detectCores())
  res[, del.rxn:=as.integer(del.rxn)]
  res[, genes:=list(rxns2genes(del.rxn, model))]
  rbind(res[del.rxn==0], res[del.rxn!=0][order(-score.adj)])
}

form.mta <- function(model, v.ref, dflux, params) {
  
  n.mets <- nrow(model$S)
  n.rxns <- ncol(model$S)

  # constraint matrix
  ## for reactions meant to remain steady
  rxns.st <- which(dflux==0 & model$c!=1)
  n.st <- length(rxns.st)
  m1 <- rbind(sparseMatrix(1:n.st, rxns.st, x=1, dims=c(n.st, n.rxns)), sparseMatrix(1:n.st, rxns.st, x=-1, dims=c(n.st, n.rxns)))
  m2 <- rbind(Diagonal(n.st), Diagonal(n.st))
  S <- rbind(cbind(model$S, sparseMatrix(NULL, NULL, dims=c(n.mets, n.st))),
             cbind(m1, m2))
  ## for **reversible** reactions that is meant to have reduced fluxes (i.e. these have the potential to "overshoot", thus need to be treated specially)
  rxns.x <- which(model$lb<0 & dflux==-1)
  n.x <- length(rxns.x)
  m1 <- rbind(sparseMatrix(1:n.x, rxns.x, x=1, dims=c(n.x, n.rxns)), sparseMatrix(1:n.x, rxns.x, x=-1, dims=c(n.x, n.rxns)))
  m2 <- rbind(Diagonal(n.x), Diagonal(n.x))
  S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(n.mets+2*n.st, n.x))),
             cbind(m1, sparseMatrix(NULL, NULL, dims=c(2*n.x, n.st)), m2))

  # constraints
  rowlb <- c(model$rowlb, c(v.ref[rxns.st], -v.ref[rxns.st]), rep(0, 2*n.x))
  rowub <- c(model$rowub, rep(Inf, 2*(n.st+n.x)))
  lb <- c(model$lb, rep(0, n.st+n.x))
  ub <- c(model$ub, rep(Inf, n.st+n.x))

  # objective function and others
  ## reactions meant to change in the forward direction, excluding those in rxns.x (i.e. w/o the potential to "overshoot")
  rxns.fw0 <- v.ref>0 & dflux==1 | v.ref==0 & dflux==1 & model$lb>=0
  ## reactions meant to change in the backward direction, excluding those in rxns.x (i.e. w/o the potential to "overshoot")
  rxns.bk0 <- v.ref<0 & dflux==1 | v.ref>0 & dflux==-1 & model$lb>=0
  tmp <- rep(0, n.rxns)
  n.ch <- sum(dflux!=0)
  tmp[rxns.fw0] <- -1/n.ch
  tmp[rxns.bk0] <- 1/n.ch
  c <- c(tmp, rep(1/n.st, n.st), rep(1/n.ch, n.x))

  # for record
  ## reactions meant to change in the forward direction
  fw <- v.ref>0 & dflux==1 | v.ref<0 & dflux==-1
  ## reactions meant to change in the backward direction
  bk <- v.ref<0 & dflux==1 | v.ref>0 & dflux==-1
  ## **reversible** reactions with v.ref==0 and dflux==1: it's fine for them to change either forward or backward; Note that I did not include these in the objective function: they would require maximize |v|, which is slightly tricker; I simply neglect such cases and adjust for them later in the scores
  fw.or.bk <- v.ref==0 & dflux==1 & model$lb<0

  # return MTA model
  list(v.ref=v.ref, dflux=dflux,
       fw=fw, bk=bk, fw.or.bk=fw.or.bk, rxns.st=rxns.st, n.ch=n.ch,
       c=c, S=S, rowlb=rowlb, rowub=rowub, lb=lb, ub=ub)
}

run.mta <- function(model, del, params, ncores=detectCores()) {
  
  names(del) <- del
  res <- mclapply(del, function(i) {
    lp.res <- run.lp(model, i, params)
    analyz.mta.res(model, lp.res)
  }, mc.cores=ncores)

  # close CPLEX
  Rcplex.close()

  rbindlist(res, idcol="del.rxn")
}

run.lp <- function(model, del, params) {

  cvec <- model$c
  objsense <- "min"
  Amat <- rbind(model$S, model$S)
  bvec <- c(model$rowlb, model$rowub)
  sense <- rep(c("G","L"), c(length(model$rowlb), length(model$rowub)))
  lb <- model$lb
  lb[del] <- 0 # if del==0, nothing will be changed to lb, meaning do not delete any reaction (the control)
  ub <- model$ub
  ub[del] <- 0 # if del==0, nothing will be changed to ub, meaning do not delete any reaction (the control)
  if ("n" %in% names(params)) {
    n <- params$n
    params$n <- NULL
  } else n <- 1
  
  tryCatch(
    {
      res <- Rcplex(cvec=cvec, objsense=objsense, Amat=Amat, bvec=bvec, sense=sense, lb=lb, ub=ub, control=params, n=n)
      if (res$status!=1) warning("MTA: Potential problem running LP for del=", del, ". Solver status: ", res$status, ".\n")
      # Not sure why the returned res$obj are NaNs... just manually recover obj
      if (is.na(res$obj)) res$obj <- sum(model$c * res$xopt)
      res[c("xopt","obj","status")]
    },
    error=function(e) {
      warning("MTA: Failed running LP for del=", del, ". Message: ", e, "NA returned.\n")
      list(xopt=NA, obj=NA, status=NA)
    }
  )
}

analyz.mta.res <- function(model, lp.res) {
 
  v0 <- model$v.ref
  v <- lp.res$xopt[1:length(v0)]
  fw <- model$fw
  bk <- model$bk
  fw.or.bk <- model$fw.or.bk
  
  # reactions intended to change: overshoot (thus regarded as failed)
  fw.overdo <- fw & v0<0 & v>(-v0)
  bk.overdo <- bk & v0>0 & v<(-v0)
  overdo <- which(fw.overdo|bk.overdo)
  # reactions intended to change: successful
  fw.yes <- fw & v>v0 & !fw.overdo
  bk.yes <- bk & v<v0 & !bk.overdo
  fw.or.bk.yes <- fw.or.bk & v!=0
  yes <- which(fw.yes|bk.yes|fw.or.bk.yes)
  # reactions intended to change: failed
  fw.no <- fw & v<v0
  bk.no <- bk & v>v0
  fw.or.bk.no <- fw.or.bk & v==0
  no <- which(fw.no|bk.no|fw.or.bk.no)

  # absolute differences between v and v.ref for the different sets of reactions
  adv.yes <- abs(v[yes]-v0[yes])
  adv.no <- abs(v[no]-v0[no])
  adv.overdo <- abs(abs(v[overdo])-abs(v0[overdo]))
  adv.st <- abs(v[model$rxns.st]-v0[model$rxns.st]) # reactions meant to remain steady

  # raw score
  s.raw <- -lp.res$obj
  # adjusted score
  s.adj <- s.raw + (sum(v0[bk])-sum(v0[fw]))/model$n.ch + sum(v[fw.or.bk])/model$n.ch
  # raw score recalculated
  #xxxx
  # original mta score
  s.ch <- sum(adv.yes)-sum(adv.no)-sum(adv.overdo)
  s.st <- sum(adv.st)
  s.mta <- s.ch/s.st
  # alternative score (approximately equivalent to adjusted score??)
  s.alt <- s.ch/(length(yes)+length(no)+length(overdo)) - s.st/length(model$rxns.st)

  # return
  data.table(solv.stat=lp.res$status, v.opt=list(v), v.opt.full=list(lp.res$xopt), rxns.change.yes=list(yes), rxns.change.no=list(no), rxns.change.overdo=list(overdo), advs.change.yes=list(adv.yes), advs.change.no=list(adv.no), advs.change.overdo=list(adv.overdo), advs.steady=list(adv.st), score.raw=s.raw, score.adj=s.adj, score.mta=s.mta)
}