library(gurobi)

lp.pars2 <- list(OutputFlag=0, TimeLimit=120, Threads=1)

mtal2 <- function(model, v.ref, dflux, del="default", lp.params=lp.pars2, ncores=nc) {
  
  # formulate MTA model
  mtal.model <- form.mtal(model, v.ref, dflux)

  # run the MTA MIQP
  if (length(del)==1 && del=="default") del <- 0:ncol(model$S)
  run.mtal2(mtal.model, del, lp.params, ncores)
}

run.mtal2 <- function(model, del, params, ncores) {
  
  #x0 <- run.lp(model, 0, NULL, params)$xopt # warm start solution
  #if (length(x0)==1 && is.na(x0)) x0 <- NULL
  names(del) <- del
  res <- mclapply(del, function(i) {
    lp.res <- run.lp2(model, i, NULL, params)
    analyz.mtal.res(model, lp.res)
  }, mc.cores=ncores)

  # close CPLEX
  Rcplex.close()

  res <- rbindlist(res, idcol="del.rxn")

  # if parallelling, the warning messages from each core won't show up, so give a summary here if any
  e <- res[is.na(solv.stat) | solv.stat!="OPTIMAL"]
  if (nrow(e)>0) {
    for (i in 1:nrow(e)) {
      warning("MTA: Potential problem or failed running LP for del=", e[i, del.rxn], ". Solver status: ", e[i, solv.stat], ".\n")
    }
  }

  res[, del.rxn:=as.integer(del.rxn)] # just in case later we need to use rxns2genes, which requires numeric rxn indeces
  rbind(res[del.rxn==0], res[del.rxn!=0][order(-score.adj)])
}

run.lp2 <- function(model, del, x0, params) {

  gm <- list()
  gm$obj <- model$c
  gm$sense <- "min"
  gm$A <- rbind(model$S, model$S)
  gm$rhs <- c(model$rowlb, model$rowub)
  gm$sense <- rep(c(">","<"), c(length(model$rowlb), length(model$rowub)))
  gm$lb <- model$lb
  gm$lb[del] <- 0 # if del==0, nothing will be changed to lb, meaning do not delete any reaction (the control)
  gm$ub <- model$ub
  gm$ub[del] <- 0 # if del==0, nothing will be changed to ub, meaning do not delete any reaction (the control)
  
  tryCatch(
    {
      res <- gurobi(gm, params)
      if (res$status!="OPTIMAL") warning("MTA: Potential problem running LP for del=", del, ". Solver status: ", res$status, ".\n")
      # Not sure why the returned res$obj are NaNs... just manually recover obj
      if (is.na(res$objval)) res$objval <- sum(model$c * res$x)
      res <- res[c("x","objval","status")]
      names(res) <- c("xopt", "obj", "status")
      res
    },
    error=function(e) {
      warning("MTA: Failed running LP for del=", del, ". Message: ", e, "NA returned.\n")
      list(xopt=NA, obj=NA, status=as.character(e))
    }
  )
}
