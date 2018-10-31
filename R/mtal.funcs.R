library(Matrix)
library(data.table)
library(Rcplex.my)
library(parallel)

lp.pars <- list(trace=0, maxcalls=5000, tilim=120, threads=1)
#nc <- detectCores() # can be unreliable??
nc <- 4L
 
mtal <- function(model, v.ref, dflux, del="default", lp.params=lp.pars, ncores=nc) {
  
  # formulate MTA model
  mtal.model <- form.mtal(model, v.ref, dflux)

  # run the MTA MIQP
  if (length(del)==1 && del=="default") del <- 0:ncol(model$S)
  run.mtal(mtal.model, del, lp.params, ncores)
}

form.mtal <- function(model, v.ref, dflux) {
  
  n.mets <- nrow(model$S)
  n.rxns <- ncol(model$S)

  # constraint matrix
  ## for reactions meant to remain steady
  st <- which(dflux==0 & model$c!=1)
  n.st <- length(st)
  m1 <- rbind(sparseMatrix(1:n.st, st, x=1, dims=c(n.st, n.rxns)), sparseMatrix(1:n.st, st, x=-1, dims=c(n.st, n.rxns)))
  m2 <- rbind(Diagonal(n.st), Diagonal(n.st))
  S <- rbind(cbind(model$S, sparseMatrix(NULL, NULL, dims=c(n.mets, n.st))),
             cbind(m1, m2))
  ## for **reversible** reactions that is meant to have reduced fluxes (i.e. these have the potential to "overshoot", thus need to be treated specially)
  rvdn.b <- model$lb<0 & dflux<0
  rvdn <- which(rvdn.b)
  n.rvdn <- length(rvdn)
  m1 <- rbind(sparseMatrix(1:n.rvdn, rvdn, x=1, dims=c(n.rvdn, n.rxns)), sparseMatrix(1:n.rvdn, rvdn, x=-1, dims=c(n.rvdn, n.rxns)))
  m2 <- rbind(Diagonal(n.rvdn), Diagonal(n.rvdn))
  S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(n.mets+2*n.st, n.rvdn))),
             cbind(m1, sparseMatrix(NULL, NULL, dims=c(2*n.rvdn, n.st)), m2))

  # constraints
  rowlb <- c(model$rowlb, c(v.ref[st], -v.ref[st]), rep(0, 2*n.rvdn))
  rowub <- c(model$rowub, rep(Inf, 2*(n.st+n.rvdn)))
  lb <- c(model$lb, rep(0, n.st+n.rvdn))
  ub <- c(model$ub, rep(Inf, n.st+n.rvdn))

  # objective function and others
  ## reactions meant to change in the forward direction, excluding those in rvdn (i.e. w/o the potential to "overshoot")
  fw0.b <- v.ref>0 & dflux>0 | v.ref==0 & dflux>0 & model$lb>=0
  ## reactions meant to change in the backward direction, excluding those in rvdn (i.e. w/o the potential to "overshoot")
  bk0.b <- v.ref<0 & dflux>0 | v.ref>0 & dflux<0 & model$lb>=0
  w <- abs(dflux) / sum(abs(dflux), na.rm=TRUE) # weight
  w[is.na(w)] <- 0
  c <- c(ifelse(fw0.b, -w, ifelse(bk0.b, w, 0)), rep(1/n.st, n.st), w[rvdn.b])

  # things to keep for downstream analysis of MTA score
  ## reactions meant to change in the forward direction
  fw <- which(v.ref>0 & dflux>0 | v.ref<0 & dflux<0)
  ## reactions meant to change in the backward direction
  bk <- which(v.ref<0 & dflux>0 | v.ref>0 & dflux<0)
  ## **reversible** reactions with v.ref==0 and dflux==1: it's fine for them to change either forward or backward; Note that I did not include these in the objective function: they would require maximize |v|, which is slightly tricker; I simply neglect such cases and adjust for them later in the scores
  fw.or.bk <- which(v.ref==0 & dflux>0 & model$lb<0)
  ## number reactions meant to change
  n.ch <- sum(dflux!=0, na.rm=TRUE)

  # return MTA model
  list(v.ref=v.ref, dflux=dflux, w=w,
       fw0.b=fw0.b, bk0.b=bk0.b, rvdn.b=rvdn.b, fw=fw, bk=bk, fw.or.bk=fw.or.bk, st=st, n.ch=n.ch, n.st=n.st,
       c=c, S=S, rowlb=rowlb, rowub=rowub, lb=lb, ub=ub)
}

run.mtal <- function(model, del, params, ncores) {
  
  names(del) <- del
  res <- mclapply(del, function(i) {
    lp.res <- run.lp(model, i, params)
    analyz.mtal.res(model, lp.res)
  }, mc.cores=ncores)

  # close CPLEX
  Rcplex.close()

  res <- rbindlist(res, idcol="del.rxn")

  # if parallelling, the warning messages from each core won't show up, so give a summary here if any
  e <- res[is.na(solv.stat) | solv.stat!="1"]
  if (nrow(e)>0) {
    for (i in 1:nrow(e)) {
      warning("MTA: Potential problem or failed running LP for del=", e[i, del.rxn], ". Solver status: ", e[i, solv.stat], ".\n")
    }
  }

  res[, del.rxn:=as.integer(del.rxn)] # just in case later we need to use rxns2genes, which requires numeric rxn indeces
  rbind(res[del.rxn==0], res[del.rxn!=0][order(-score.adj)])
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
      res$status <- as.character(res$status)
      if (res$status!="1") warning("MTA: Potential problem running LP for del=", del, ". Solver status: ", res$status, ".\n")
      # Not sure why the returned res$obj are NaNs... just manually recover obj
      if (is.na(res$obj)) res$obj <- sum(model$c * res$xopt)
      res[c("xopt","obj","status")]
    },
    error=function(e) {
      warning("MTA: Failed running LP for del=", del, ". Message: ", e, "NA returned.\n")
      list(xopt=NA, obj=NA, status=as.character(e))
    }
  )
}

analyz.mtal.res <- function(model, lp.res) {

  if (is.na(lp.res$xopt) && is.na(lp.res$obj)) {
    return(data.table(solv.stat=lp.res$status, v.opt=NA, v.opt.full=NA, rxns.change.yes=NA, rxns.change.no=NA, advs.change.yes=NA, advs.change.no=NA, advs.steady=NA, score.raw=NA, score.adj=NA, score.rec=s.rec, score.mta=NA))
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

  s.yes <- sum(adv.yes)
  s.no <- sum(adv.no)
  s.st <- sum(adv.st)

  # raw score: just the (negated) optimal objective value
  s.raw <- -lp.res$obj
  # adjusted score
  s.adj <- s.raw + (sum(v0[model$bk]) - sum(v0[model$fw]))/model$n.ch + sum(v[model$fw.or.bk])/model$n.ch
  # recalculated score (should be the same as the adjusted score)
  s.rec <- (s.yes - s.no + sum(v[model$fw.or.bk]))/model$n.ch - s.st/model$n.st
  # ratio score like the original mta
  s.mta <- (s.yes-s.no)/s.st

  # return
  data.table(solv.stat=lp.res$status, v.opt=list(v), v.opt.full=list(lp.res$xopt), rxns.change.yes=list(yes), rxns.change.no=list(no), advs.change.yes=list(adv.yes), advs.change.no=list(adv.no), advs.steady=list(adv.st), score.raw=s.raw, score.adj=s.adj, score.rec=s.rec, score.mta=s.mta)
}