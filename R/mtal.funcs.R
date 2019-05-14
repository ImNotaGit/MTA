library(data.table)
library(Matrix)
library(parallel)
library(Rcplex.my)

lp.pars <- list(trace=0, maxcalls=5000, tilim=120, threads=1)

 
mtal <- function(model, v.ref, dflux, del="default", lp.params=lp.pars, nc=1L) {
  
  # formulate MTA model
  mtal.model <- form.mtal(model=model, v.ref=v.ref, dflux=dflux)

  # run MTA
  if (length(del)==1 && del=="default") del <- 0:ncol(model$S)
  run.mtal(model=mtal.model, del=del, params=lp.params, nc=nc)
}

form.mtal <- function(model, v.ref, dflux) {
  
  n.mets <- nrow(model$S)
  n.rxns <- ncol(model$S)

  # constraint matrix
  ## for reactions meant to remain steady
  st <- which(dflux==0 & model$c!=1)
  n.st <- length(st)
  S <- rbind(
    cbind( model$S,                                              sparseMatrix(NULL, NULL, dims=c(n.mets, 2*n.st)) ),
    cbind( sparseMatrix(1:n.st, st, x=1, dims=c(n.st, n.rxns)),  Diagonal(n.st),  Diagonal(n.st, -1) )
  )
  ## for **reversible** reactions that is meant to have reduced fluxes (i.e. these have the potential to "overshoot", thus need to be treated specially)
  rvdn.b <- model$lb<0 & dflux<0
  rvdn <- which(rvdn.b)
  n.rvdn <- length(rvdn)
  if (n.rvdn>0) {
    S <- rbind(
      cbind( S,                                                                                                             sparseMatrix(NULL, NULL, dims=c(n.mets+n.st, 2*n.rvdn)) ),
      cbind( sparseMatrix(1:n.rvdn, rvdn, x=1, dims=c(n.rvdn, n.rxns)),  sparseMatrix(NULL, NULL, dims=c(n.rvdn, 2*n.st)),  Diagonal(n.rvdn),  Diagonal(n.rvdn, -1) )
    )
  }
  
  # constraints
  rowlb <- c(model$rowlb, v.ref[st], rep(0, n.rvdn))
  rowub <- c(model$rowub, v.ref[st], rep(0, n.rvdn))
  lb <- c(model$lb, rep(0, 2*n.st+2*n.rvdn))
  ub <- c(model$ub, rep(Inf, 2*n.st+2*n.rvdn))

  # objective function and others
  ## reactions meant to change in the forward direction, excluding those in rvdn (i.e. w/o the potential to "overshoot")
  fw0.b <- v.ref>0 & dflux>0 | v.ref==0 & dflux>0 & model$lb>=0
  ## reactions meant to change in the backward direction, excluding those in rvdn (i.e. w/o the potential to "overshoot")
  bk0.b <- v.ref<0 & dflux>0 | v.ref>0 & dflux<0 & model$lb>=0
  w <- abs(dflux) / sum(abs(dflux), na.rm=TRUE) # weight
  c <- c(ifelse(fw0.b, -w, ifelse(bk0.b, w, 0)), rep(1/n.st, 2*n.st), w[rvdn], w[rvdn])
  c[is.na(c)] <- 0

  # things to keep for downstream analysis of MTA score
  ## reactions meant to change in the forward direction
  fw <- which(v.ref>0 & dflux>0 | v.ref<0 & dflux<0)
  ## reactions meant to change in the backward direction
  bk <- which(v.ref<0 & dflux>0 | v.ref>0 & dflux<0)
  ## **reversible** reactions with v.ref==0 and dflux==1: it's fine for them to change either forward or backward; Note that I did not include these in the objective function: they would require maximize |v|, which is slightly tricker; I simply neglect such cases and adjust for them later in the scores
  fw.or.bk <- which(v.ref==0 & dflux>0 & model$lb<0)

  # return MTA model
  list(v.ref=v.ref, dflux=dflux, w=w,
       fw0.b=fw0.b, bk0.b=bk0.b, rvdn.b=rvdn.b, fw=fw, bk=bk, fw.or.bk=fw.or.bk, st=st, n.st=n.st,
       c=c, S=S, rowlb=rowlb, rowub=rowub, lb=lb, ub=ub)
}

run.mtal <- function(model, del, params, nc) {
  
  #x0 <- run.lp(model, 0, NULL, params)$xopt # warm start solution
  #if (length(x0)==1 && is.na(x0)) x0 <- NULL
  names(del) <- del
  res <- mclapply(del, function(i) {
    lp.res <- run.lp(model=model, del=i, x0=NULL, params=params)
    analyz.mtal.res(model=model, lp.res=lp.res)
  }, mc.cores=nc)

  # close CPLEX
  Rcplex.close()

  res <- rbindlist(res, idcol="del.rxn")
  res[, del.rxn:=as.integer(del.rxn)]

  # give a summary here of (potential) optimization issues if any
  e <- res[solv.stat!="1", .(del.rxn, solv.stat)]
  if (nrow(e)>0) {
    cat("MTA: Potential problem or failed running LP for the case(s) below:\n")
    print(e)
  }
  
  #res <- rbind(res[del.rxn==0], res[del.rxn!=0][order(-score.mta)])
  res
}

run.lp <- function(model, del, x0, params) {

  cvec <- model$c
  Amat <- rbind(model$S, model$S)
  bvec <- c(model$rowlb, model$rowub)
  sense <- rep(c("G","L"), c(length(model$rowlb), length(model$rowub)))
  lb <- model$lb
  lb[del] <- 0 # if del==0, nothing will be changed to lb, meaning do not delete any reaction (the control)
  ub <- model$ub
  ub[del] <- 0 # if del==0, nothing will be changed to ub, meaning do not delete any reaction (the control)
  
  tryCatch(
    {
      res <- Rcplex(cvec=cvec, objsense="min", Amat=Amat, bvec=bvec, sense=sense, lb=lb, ub=ub, x0=x0, control=params)
      list(xopt=res$xopt, obj=ifelse(is.na(res$obj), sum(model$c*res$xopt, na.rm=TRUE), res$obj), status=as.character(res$status))
    },
    error=function(e) list(xopt=NA, obj=NA, status=as.character(e))
  )
}

analyz.mtal.res <- function(model, lp.res) {

  if (is.na(lp.res$xopt) && is.na(lp.res$obj)) {
    return(data.table(solv.stat=lp.res$status, v.opt=NA, v.opt.full=NA, rxns.change.yes=NA, rxns.change.no=NA, advs.change.yes=NA, advs.change.no=NA, advs.steady=NA, score.raw=NA, score.adj=NA, score.mta=NA, score.mta.uw=NA))
  }
 
  v0 <- model$v.ref
  v <- lp.res$xopt[1:length(v0)]
  fw0 <- model$fw0.b
  bk0 <- model$bk0.b
  rvdn <- model$rvdn.b
  
  # reactions intended to change: successful
  fw0.yes <- fw0 & v>v0
  bk0.yes <- bk0 & v<v0
  rvdn.yes <- rvdn & abs(v)<abs(v0)
  yes <- which(fw0.yes|bk0.yes|rvdn.yes)

  # reactions intended to change: failed
  fw0.no <- fw0 & v<v0
  bk0.no <- bk0 & v>v0
  rvdn.no <- rvdn & abs(v)>abs(v0)
  no <- which(fw0.no|bk0.no|rvdn.no)

  # absolute differences between v and v.ref for the different sets of reactions
  ## reactions intended to change
  adv.yes <- ifelse(yes %in% which(rvdn.yes), abs(v0[yes])-abs(v[yes]), abs(v[yes]-v0[yes]))
  adv.no <- ifelse(no %in% which(rvdn.no), abs(v[no])-abs(v0[no]), abs(v[no]-v0[no]))
  ## reactions intended to remain steady
  adv.st <- abs(v[model$st]-v0[model$st])

  # score for reactions intended to change
  s.ch <- sum(adv.yes * model$w[yes]) + sum(v[model$fw.or.bk] * model$w[model$fw.or.bk]) - sum(adv.no * model$w[no]) # score.yes - score.no
  # score for reactions intended to remain steady
  s.st.uw <- sum(adv.st) # un-weighted
  s.st <- s.st.uw / model$n.st

  # raw score: just the (negated) optimal objective value
  s.raw <- -lp.res$obj
  # adjusted score
  s.adj <- s.raw + sum(v0[model$bk] * model$w[model$bk]) - sum(v0[model$fw] * model$w[model$fw]) + sum(v[model$fw.or.bk] * model$w[model$fw.or.bk])
  # recalculated score (should be the same as the adjusted score) ## yes, same
  #s.rec <- s.ch - s.st
  # ratio score like the original mta
  s.mta <- s.ch / s.st
  s.mta.uw <- (sum(adv.yes) + sum(v[model$fw.or.bk]) - sum(adv.no)) / s.st.uw # un-weighted

  # return
  data.table(solv.stat=lp.res$status, v.opt=list(v), v.opt.full=list(lp.res$xopt), rxns.change.yes=list(yes), rxns.change.no=list(no), advs.change.yes=list(adv.yes), advs.change.no=list(adv.no), advs.steady=list(adv.st), score.raw=s.raw, score.adj=s.adj, score.mta=s.mta, score.mta.uw=s.mta.uw)
}