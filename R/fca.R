prepare.model.for.fca <- function(model) {
  # prepare metabolic model for flux coupling analysis
  # first convert each reversible rxn in model into two irreversible rxns representing the positive and negative halves, both with lb==0
  # then remove blocked reactions
  
  # convert reversible reactions
  model <- as.list(model)
  rev.idx <- which(model$lb<0)
  model$S <- cbind(model$S, -model$S[, rev.idx])
  rev.lb <- model$lb[rev.idx]
  model$lb[rev.idx] <- 0
  model$lb <- c(model$lb, rep(0, length(rev.idx)))
  model$ub <- c(model$ub, -rev.lb)
  model$rxns <- c(model$rxns, paste0(model$rxns[rev.idx],"_REV"))
  model$rxnNames <- c(model$rxnNames, paste0(model$rxnNames[rev.idx],"_REV"))
  model$rules <- c(model$rules, model$rules[rev.idx])
  # for now, not updating model$rxnGeneMat
  
  # remove blocked reactions
  ubs <- sapply(1:ncol(model$S), get.opt.flux, model=model)
  idx <- which(ubs>0)
  model <- subset.model(model, 1:nrow(model$S), idx)
}

fca <- function(model, rxn0, rxn1) {
  # flux coupling analysis, just for two reactions; the formalization is rxn1/rxn0 (i.e. normalized by rxn0)

  # formalize fca model
  model$lb <- rep(0, length(model$lb))
  ub <- model$ub
  model$ub <- rep(1e100, length(model$ub))
  model$lb[rxn0] <- 1
  model$ub[rxn0] <- 1
  ## transport reactions
  trxns <- which(apply(model$S, 2, function(x) all(x>=0) || all(x<=0)))
  n.trxns <- length(trxns)
  nr <- nrow(model$S)
  nc <- ncol(model$S)
  model$S <- rbind(cbind(model$S,                                             rep(0, nr)),
                   cbind(sparseMatrix(1:n.trxns, trxns, dims=c(n.trxns, nc)), -ub[trxns]))
  
  model$lb <- c(model$lb, 0)
  model$ub <- c(model$ub, 1e100)
  model$rowlb <- c(model$rowlb, rep(-1e100, n.trxns))
  model$rowub <- c(model$rowub, rep(0, n.trxns))
  
  # run
  rmax <- get.opt.flux(model, rxn1, dir="max")
  rmin <- get.opt.flux(model, rxn1, dir="min")
  c(rmin, rmax)
}
