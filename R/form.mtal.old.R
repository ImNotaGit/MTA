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
  if (n.rvdn>0) {
    m1 <- rbind(sparseMatrix(1:n.rvdn, rvdn, x=1, dims=c(n.rvdn, n.rxns)), sparseMatrix(1:n.rvdn, rvdn, x=-1, dims=c(n.rvdn, n.rxns)))
    m2 <- rbind(Diagonal(n.rvdn), Diagonal(n.rvdn))
    S <- rbind(cbind(S, sparseMatrix(NULL, NULL, dims=c(n.mets+2*n.st, n.rvdn))),
               cbind(m1, sparseMatrix(NULL, NULL, dims=c(2*n.rvdn, n.st)), m2))
  }
  
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
  c <- c(ifelse(fw0.b, -w, ifelse(bk0.b, w, 0)), rep(1/n.st, n.st), w[rvdn])
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