library(data.table)
library(Matrix)
library(stringr)
library(parallel)
library(Rcplex.my)


#### ---- metabolic model utils ----

exprs2rxns <- function(vec, type=0, model, discrt=TRUE, na.replace=TRUE) {
  # map a vector of expression values either meant to be levels or direction of changes of genes to those of the reactions in the metabolic model
  # vec is a named vector of, names being gene symbol; or if it's unnamed and length being model$genes, assume it's already in the same order as model$genes
  # type==0 for vec as gene levels, type==1 for vec as directions of changes
  # discrt==TRUE for returning a flux vector or -1/0/1, otherwise the result will be continuous values
  # NA's will be replaced by zeros if na.replace==TRUE

  if (is.null(names(vec)) && length(vec)==length(model$genes)) x <- vec else x <- vec[model$genes]
  if (discrt) x <- sign(x)
  if (type==0) {
    `&` <- function(a,b) {
      #if (isTRUE(is.na(a) && b<0)) return(b)
      #if (isTRUE(is.na(b) && a<0)) return(a) # if one is NA and the other <0, for sure the result is the <0 value; all other NA cases are undetermined and NA will be returned
      return(min(a,b))
    }
    `|` <- function(a,b) {
      #if (isTRUE(is.na(a) && b>0)) return(b)
      #if (isTRUE(is.na(b) && a>0)) return(a) # if one is NA and the other >0, for sure the result is the >0 value; all other NA cases are undetermined and NA will be returned
      return(max(a,b))
    }
  } else if (type==1) {
    `&` <- function(a,b) { # if one is NA and the other is 0, for sure the result is 0; all other NA cases are undetermined and NA will be returned
      if (isTRUE(sign(a)!=sign(b) || a==0 || b==0)) return(0)
      return(min(a,b))
    }
    `|` <- function(a,b) { # all NA cases are undetermined and NA will be returned
      if (isTRUE(sign(a)==sign(b))) return(max(a,b))
      return(abs(sign(a)+sign(b))*(a+b))
    }
  }
  res <- sapply(model$rules, function(i) eval(parse(text=i)))
  if (na.replace) res[is.na(res)] <- 0
  if (discrt) res <- as.integer(res)
  unname(res)
}

rxns2genes <- function(vec, model) {
  # given a numerical vector corresponding to the reaction indeces in the model, map each of them to gene names, return as a list; NA will be returned for reaction indeces outside the proper range (including 0)
  vec[vec==0] <- NA
  res <- lapply(str_extract_all(model$rules[vec], "[1-9][0-9]*"), function(x) unique(model$genes[as.numeric(x)]))
  #if (length(res)==1) res <- res[[1]]
  
  res
}

genes2rxns <- function(genes, type=0, model) {
  # given a vector of gene symbols, for each of them map to reactions (rxn indeces in the model), return as a list
  # type==0 for any reactions involving the gene; type==1 for reactions where the gene is essential (corresponding to that if the gene level is low, then the reaction flux should be low)
  if (type==0) {
    r2g <- rxns2genes(1:length(model$rules), model)
    res <- lapply(genes, function(gi) which(sapply(r2g, function(gns) gi %in% gns)))
  } else if (type==1) {
    gind <- match(genes, model$genes)
    res <- lapply(gind, function(gi) {
      tmp <- rep(0, length(model$genes))
      tmp[gi] <- -1
      which(exprs2rxns(tmp, 0, model)==-1)
    })
  }

  res
}

rxns2mets <- function(vec, model) {
  # given a vector of reaction indices, map each of them to metabolite indices
  res <- lapply(vec, function(i) which(model$S[,i]!=0))
  names(res) <- vec
  res
}

mets2rxns <- function(vec, model) {
  # given a vector of metabolite indices, map each of them to reaction indices
  res <- lapply(vec, function(i) which(model$S[i,]!=0))
  names(res) <- vec
  res
}

get.rxn.equation <- function(vec, model) {
  # given a numerical vector corresponding to the reaction indeces in the model, return string vector containing the equations of the corresponding reactions
  sapply(vec, function(i) {
    x <- model$S[, i]
    rs <- paste(trimws(paste(ifelse(x[x<0]==-1,"",-x[x<0]), model$mets[x<0])), collapse=" + ")
    ps <- paste(trimws(paste(ifelse(x[x>0]==1,"",x[x>0]), model$mets[x>0])), collapse=" + ")
    if (model$lb[i]>=0) arrow <- "-->" else arrow <- "<==>"
    paste(rs, arrow, ps)
  })
}

get.neighborhood <- function(model, ids, order=1, type="rxn", exclude.mets.default="^h[[_].\\]?|^oh1[[_].\\]?|^h2o[[_].\\]?|^atp[[_].\\]?|^adp[[_].\\]?|^pi[[_].\\]?|^ppi[[_].\\]?|^coa[[_].\\]?|^o2[[_].\\]?|^co2[[_].\\]?|^nadp[[_].\\]?|^nadph[[_].\\]?|^nad[[_].\\]?|^nadh[[_].\\]?|^fad[[_].\\]?|^fadh2[[_].\\]?|^na1[[_].\\]?|^so4[[_].\\]?|^nh4[[_].\\]?|^cl[[_].\\]?", exclude.mets.degree=nrow(model$S), exclude.mets=NULL) {
  # given the IDs of a set of rxns or mets (specified by type), return the IDs of the rxns or mets with distance<=order from each of the given ones (as a list). by default order=1 means the rxns sharing a met or the mets within the same rxn.
  # exclude.mets.default: regex of some high degree metabolites to be excluded by default when determining the neighborhood; the regex works for mets formats like "h[c]" and "h_c"
  # exclude.mets.degree: exclude metabolites with degree greater than this
  # exclude.mets are IDs (indices) of other mets to be excluded
  
  library(igraph)
  s <- model$S
  # exclude metabolites
  tmp <- graph.incidence(s)
  tmp <- bipartite.projection(tmp, which="false") # projected network of metabolites
  deg <- igraph::degree(tmp) # make sure use degree from igraph
  exclude.mets <- c(exclude.mets, which(deg>exclude.mets.degree))
  if (!(is.null(exclude.mets.default) || exclude.mets.default=="")) {
    emd <- grep(exclude.mets.default, model$mets)
    exclude.mets <- unique(c(exclude.mets, emd))
  }
  cat("The following metabolites are excluded:\n")
  cat(paste(model$mets[exclude.mets], collapse=", "), "\n")
  s[exclude.mets, ] <- 0
  # create bipartite graph between mets and rxns
  gb <- graph.incidence(s)
  if (type=="rxn") {
    # project bipartite graph into the graph of rxns
    gp <- bipartite.projection(gb, which="true")
  } else if (type=="met") {
    # project bipartite graph into the graph of mets
    gp <- bipartite.projection(gb, which="false")
  }
  # get the neighborhood
  res <- lapply(ego(gp, order=order, nodes=ids), as.vector)
  names(res) <- ids
  res
}

subset.model <- function(model, i, j) {
  # subset model like a matrix, i for metabolites, j for reactions
  # the approach below will not work properly if e.g. the number of reactions happens to be equal to the number of metabolites, of any of them happens to be equal to the number of genes.
  m <- nrow(model$S)
  n <- ncol(model$S)
  p <- length(model$genes)
  if (m==n || n==p || m==p) stop("subset.model: sorry this case won't work.\n")
  res <- lapply(model, function(x) {
    if (length(dim(x))==2) {
      if (nrow(x)==m) xx <- x[i, ] else xx <- x
      if (ncol(x)==n) xx <- xx[, j]
    } else if (length(x)==m) {
      xx <- x[i]
    } else if (length(x)==n) {
      xx <- x[j]
    } else {
      xx <- x
    }
    xx
  })
}

c.model <- function(model1, model2) {
  # simply concatenate two models, for now, ignore rules, grRules, rxnGeneMat
  xs <- setdiff(intersect(names(model1), names(model2)), c("rules","grRules","rxnGeneMat"))
  names(xs) <- xs
  res <- lapply(xs, function(i) {
    x1 <- model1[[i]]
    x2 <- model2[[i]]
    if (length(dim(x1))==2) {
      xx <- rbind(cbind(x1, sparseMatrix(NULL, NULL, dims=c(nrow(x1), ncol(x2)))),
                  cbind(sparseMatrix(NULL, NULL, dims=c(nrow(x2), ncol(x1))), x2))
    } else if (is.vector(x1)) {
      if (is.numeric(x1) || i=="vtype" || i=="metFormulas" || i=="csense") {
        xx <- c(x1, x2)
      } else if (i=="description") {
        xx <- paste(x1, x2)
      } else xx <- c(paste0(x1,"_1"), paste0(x2,"_2"))
    }
    xx
  })
}

get.opt.flux <- function(model, i, coef=1, dir="max", ko=NULL, keep.xopt=FALSE, nc=1L, na=FALSE) {
  # get the max or min flux of the i'th reaction in the model
  # `i` can also be a vector of multiple reaction indices, with coef being their coefficients, then the corresponding linear objective function will be optimized
  # to knockout reaction(s), pass KO as their indices
  # if keep.xopt, the optimal xopt vector will also be returned
  # if na, then will return NA if solver status is not 1 (i.e. optimal), otherwise will keep the original value
  cvec <- rep(0, ncol(model$S))
  cvec[i] <- coef
  objsense <- dir
  Amat <- rbind(model$S, model$S)
  bvec <- c(model$rowlb, model$rowub)
  sense <- rep(c("G","L"), c(length(model$rowlb), length(model$rowub)))
  lb <- model$lb
  ub <- model$ub
  lb[ko] <- 0
  ub[ko] <- 0
  
  res <- Rcplex(cvec=cvec, objsense=objsense, Amat=Amat, bvec=bvec, sense=sense, lb=lb, ub=ub, control=list(trace=0, maxcalls=10000, tilim=120, threads=nc))
  if (res$status!=1) {
    warning("Potential problems running LP. Solver status: ", res$status, ".\n")
    if (na) {
      res$obj <- NA
      res$xopt <- NA
    }
  }
  if (keep.xopt) {
    return(list(obj=res$obj, xopt=res$xopt))
  } else return(res$obj)
}

set.lb.biomass <- function(model, value=0.1, of.max=TRUE, rgx="biomass") {
  i <- grep(rgx, model$rxns, ignore.case=TRUE)
  if (length(i)==0) stop("No biomass reaction mapped by regex.\n")
  if (length(i)>1) stop("Multiple reactions mapped by regex:\n", paste(model$rxns[i], collapse="\n"))
  if (of.max) {
    max.val <- get.opt.flux(model, i)
    model$lb[i] <- value * max.val
  } else {
    model$lb[i] <- value
  }
  model
} 

rm.rxns <- function(model, vec) {
  # remove reactions from model
  # vec is a vector containing the indices of reactions to remove, or a logical vector with length equal to the number of reactions where the reactions to be removed is represented by TRUE
  if (is.numeric(vec)) {
    keeps <- setdiff(1:ncol(model$S), vec)
  } else if (is.logical(vec) && length(vec)==ncol(model$S)) {
    keeps <- !vec
  }

  model$S <- model$S[, keeps]
  model$lb <- model$lb[keeps]
  model$ub <- model$ub[keeps]

  if ("rxns" %in% names(model)) model$rxns <- model$rxns[keeps]
  if ("rxnNames" %in% names(model)) model$rxnNames <- model$rxnNames[keeps]
  if ("rxnConfidenceScores" %in% names(model)) model$rxnConfidenceScores <- model$rxnConfidenceScores[keeps]
  if ("rxnECNumbers" %in% names(model)) model$rxnECNumbers <- model$rxnECNumbers[keeps]
  if ("rxnReferences" %in% names(model)) model$rxnReferences <- model$rxnReferences[keeps]
  if ("subSystems" %in% names(model)) model$subSystems <- model$subSystems[keeps]
  if ("rev" %in% names(model)) model$rev <- model$rev[keeps]
  if ("c" %in% names(model)) model$c <- model$c[keeps]
  if ("int.vars" %in% names(model)) model$int.vars <- model$int.vars[keeps]
  if ("rules" %in% names(model)) {
    model$rules <- model$rules[keeps]
    # remove the genes that are no longer present in the model
    gns.keep <- as.numeric(unique(unlist(str_extract_all(model$rules, "[0-9]+"))))
    gns.rm <- setdiff(1:length(model$genes), gns.keep)
    model$genes[gns.rm] <- NA # set to NA to avoid disrupting the indices of rest of the genes
    if ("gene.ids" %in% names(model)) model$gene.ids[gns.rm] <- NA
  }
  if ("grRules" %in% names(model)) model$grRules <- model$grRules[keeps]
  if ("rxnGeneMat" %in% names(model)) model$rxnGeneMat <- model$rxnGeneMat[keeps, ]

  # remove the metabolites that are no longer present in the model
  mkeeps <- apply(model$S, 1, function(x) any(x!=0))
  model$S <- model$S[mkeeps, ]
  if ("mets" %in% names(model)) model$mets <- model$mets[mkeeps]
  if ("metNames" %in% names(model)) model$metNames <- model$metNames[mkeeps]
  if ("metFormulas" %in% names(model)) model$metFormulas <- model$metFormulas[mkeeps]
  if ("metCharges" %in% names(model)) model$metCharges <- model$metCharges[mkeeps]
  if ("met.production" %in% names(model)) model$met.production <- model$met.production[mkeeps]
  if ("rowlb" %in% names(model)) model$rowlb <- model$rowlb[mkeeps]
  if ("rowub" %in% names(model)) model$rowub <- model$rowub[mkeeps]
  if ("b" %in% names(model)) model$b <- model$b[mkeeps]
  if ("csense" %in% names(model)) model$csense <- model$csense[mkeeps]

  # if resulting in metabolites that are only produced or only consumed, give a warning
  pr <- apply(model$S, 1, function(x) all(x>=0 & model$lb>=0 | x<=0 & model$ub<=0))
  co <- apply(model$S, 1, function(x) all(x>=0 & model$ub<=0 | x<=0 & model$lb>=0))
  if (any(pr|co)) warning("The returned model contains metabolites that are produced or consumed only. Please manually check and fix this.\n")

  model
}

preprocess.model <- function(model, nc=1L) {
  # preprocess metabolic model by removing the fixed reactions (i.e. those in practise must carry fixed fluxes, including zero, based on the model constraints)
  ub <- unlist(mclapply(1:ncol(model$S), get.opt.flux, model=model, dir="max", mc.cores=nc))
  ub <- pmin(model$ub, ub)
  lb <- unlist(mclapply(1:ncol(model$S), get.opt.flux, model=model, dir="min", mc.cores=nc))
  lb <- pmax(model$lb, lb)
  Rcplex.close()

  `%gt%` <- function(x, y) x-y > sqrt(.Machine$double.eps) # "substantially" greater than, in the way of all.equal
  nfx <- ub %gt% lb # non-fixed cases
  fx <- !nfx # all other cases: we expect that in all these cases, ub equals or nearly equals lb
  if ("b" %in% names(model)) model$b <- as.vector(model$b - model$S[, fx] %*% lb[fx]) # for the "fixed" cases, we expect ub equals or nearly equals lb, so just use lb here to correct for b
  if ("rowlb" %in% names(model)) model$rowlb <- as.vector(model$rowlb - model$S[, fx] %*% lb[fx])
  if ("rowub" %in% names(model)) model$rowub <- as.vector(model$rowub - model$S[, fx] %*% lb[fx])
  model <- rm.rxns(model, fx)
  cat(sum(fx), "fixed reactions removed. Now the model contains", ncol(model$S), "reactions and", nrow(model$S), "metabolites.\n")
  
  model
}

get.diff.flux <- function(imat.model0, imat.model1, use.sample=TRUE, sample.range=NULL, rxns="all", nc=1L, padj.cutoff=0.01, r.cutoff=0.1, diff.med.cutoff=0.1) {
  # do differential flux analysis: imat.model1 compared to imat.model0
  # if use.sample, use the sampled flux distributions saved in the imat.models, by default, use 1001:end sample points; or specify range of sample points for model0 and model1 respectively in a list
  # is use.sample==FALSE, then deterimine diff flux by obtaining the min and max fluxes of each reaction for the two models, and regard reactions with either increased lb or increased ub as upregulated and vice versa; ambiguous cases (e.g. lowered lb and increased ub) are regarded as non-differential
  # by default, rxns="all" means performing the analysis across all reactions; or specify reaction indices; nc is the number of cores for paralleling across the reactions
  # padj.cutoff and r.cutoff and log.fc.cutoff are used to determine the significantly changed reactions when use.sample==TRUE

  if (length(rxns)==1 && rxns=="all") rxns <- 1:length(imat.model0$rxns)
  lb <- pmin(imat.model0$lb[rxns], imat.model1$lb[rxns])
  ub <- pmax(imat.model0$ub[rxns], imat.model1$ub[rxns])
  
  if (use.sample) {
    if (is.null(sample.range)) {
      sr0 <- 1001:ncol(imat.model0$sampl$pnts)
      sr1 <- 1001:ncol(imat.model1$sampl$pnts)
    } else {
      sr0 <- sample.range[[1]]
      sr1 <- sample.range[[2]]
    }
    dflux.test <- function(s0, s1, lb, ub) {
      # run wilcox test
      tryCatch({
        wilcox.res <- wilcox.test(s0, s1)
        # p value
        wilcox.p <- wilcox.res$p.value
        # effect size for unpaired test: rank biserial correlation
        wilcox.r <- unname(1 - 2 * wilcox.res$statistic / (length(s0)*length(s1)))
        # another effect size measure: "normalized" difference of median fluxes
        m0 <- median(s0)
        m1 <- median(s1)
        dm <- m1-m0
        data.table(lb0=min(s0), ub0=max(s0), med0=m0,
                   lb1=min(s1), ub1=max(s1), med1=m1,
                   diff.med=dm, r=wilcox.r, pval=wilcox.p)
      }, error=function(e) {
        data.table(lb0=NA, ub0=NA, med0=NA, lb1=NA, ub1=NA, med1=NA, diff.med=NA, r=NA, pval=NA)
      })
    }
    res <- rbindlist(mclapply(rxns, function(i) dflux.test(imat.model0$sampl$pnts[i, sr0], imat.model1$sampl$pnts[i, sr1], lb, ub), mc.cores=nc))
    res[, padj:=p.adjust(pval, method="BH")]
    res <- cbind(data.table(id=rxns, rxn=imat.model0$rxns[rxns]), res)
    res <- res[order(-abs(diff.med), -abs(r), padj, pval)]
    # add summary of flux differences: positive value means flux value changes towards the positive side, vice versa; 0 means unchanged
    res[, dir:=ifelse(!(padj<padj.cutoff & abs(r)>quantile(abs(r), r.cutoff, na.rm=TRUE) & abs(diff.med)>quantile(abs(diff.med), diff.med.cutoff, na.rm=TRUE)), 0, ifelse(r>0, 1, -1))]
  } else {
    ub0 <- unlist(mclapply(rxns, get.opt.flux, model=imat.model0, dir="max", mc.cores=nc))
    lb0 <- unlist(mclapply(rxns, get.opt.flux, model=imat.model0, dir="min", mc.cores=nc))
    ub1 <- unlist(mclapply(rxns, get.opt.flux, model=imat.model1, dir="max", mc.cores=nc))
    lb1 <- unlist(mclapply(rxns, get.opt.flux, model=imat.model1, dir="min", mc.cores=nc))
    Rcplex.close()
    m0 <- (ub0+lb0)/2
    m1 <- (ub1+lb1)/2
    dm <- m1-m0
    res <- data.table(id=rxns, rxn=imat.model0$rxns[rxns],
                      lb0=lb0, ub0=ub0, med0=m0,
                      lb1=lb1, ub1=ub1, med1=m1,
                      diff.med=dm)
    res <- res[order(-abs(diff.med))]
    # add summary of flux differences: positive value means flux value changes towards the positive side, vice versa; 0 means unchanged
    res[, dir:=ifelse(ub1>ub0 & lb1>lb0, 3,
               ifelse(ub1<ub0 & lb1<lb0, -3,
               ifelse(ub1>ub0 & lb1==lb0 | ub1==ub0 & lb1>lb0, 2,
               ifelse(ub1<ub0 & lb1==lb0 | ub1==ub0 & lb1<lb0, -2,     
               ifelse(med1>med0, 1,
               ifelse(med1<med0, -1, 0))))))]
  }
  
  res
}

get.diff.flux.x2 <- function(imat.model, c0="cell2", c1="cell1", use.sample=TRUE, sample.range=NULL, nc=1L, padj.cutoff=0.01, r.cutoff=0.1, diff.med.cutoff=0.1) {
  # diff.flux for bicellular model, between the two cells, i.e. c1 (_cell1) vs c0 (_cell2)
  # will just re-use get.diff.flux; but perform analysis always for all reactions that are not exlusively in the extracellular space
  # in the result, reaction id correspond to those in the single cellular model
  r0 <- grep(c0, imat.model$rxns)
  r1 <- grep(c1, imat.model$rxns)
  if (c1=="cell1") rxn.ids <- r1 else rxn.ids <- r0
  rxns <- stringr::str_sub(imat.model$rxns[rxn.ids],1,-7)
  im0 <- list(rxns=rxns, lb=imat.model$lb[r0], ub=imat.model$ub[r0], sampl=list(pnts=imat.model$sampl$pnts[r0,]))
  im1 <- list(rxns=rxns, lb=imat.model$lb[r1], ub=imat.model$ub[r1], sampl=list(pnts=imat.model$sampl$pnts[r1,]))
  res <- get.diff.flux(im0, im1, use.sample, sample.range, rxns="all", nc, padj.cutoff, r.cutoff, diff.med.cutoff)
  res[, id:=rxn.ids[id]]
}

get.diff.comb.flux <- function(imat.model0, imat.model1, use.sample=TRUE, sample.range=NULL, rxns, nc=1L, padj.cutoff=0.01, r.cutoff=0.1, diff.med.cutoff=0.1) {
  # do differential flux analysis for the linear combination of fluxes of rxns; rxns is a list, each element is a named vector, the name being the rxn indices to combine, and the values being the coefficient of the linear combination; this function will do diff flux for each case corresponding to the elements of rxns.
  # if use.sample, use the sampled flux distributions saved in the imat.models, by default, use 1001:end sample points; or specify range of sample points for model0 and model1 respectively in a list
  # is use.sample==FALSE, then deterimine diff flux by obtaining the min and max fluxes through each metabolite for the two models, and regard metabolites with either increased lb or increased ub as upregulated and vice versa; ambiguous cases (e.g. lowered lb and increased ub) are regarded as non-differential
  # by default, mets="all" means performing the analysis across all metabolites; or specify metabolite indices; nc is the number of cores for paralleling across the metabolites
  # padj.cutoff and r.cutoff and log.fc.cutoff are used to determine the significantly changed metabolites when use.sample==TRUE
  if (use.sample) {
    if (is.null(sample.range)) {
      sr0 <- 1001:ncol(imat.model0$sampl$pnts)
      sr1 <- 1001:ncol(imat.model1$sampl$pnts)
    } else {
      sr0 <- sample.range[[1]]
      sr1 <- sample.range[[2]]
    }
    samp0 <- sapply(rxns, function(x) {
      i <- as.integer(names(x))
      colSums(imat.model0$sampl$pnts[i,sr0]*x)
    })
    samp1 <- sapply(rxns, function(x) {
      i <- as.integer(names(x))
      colSums(imat.model1$sampl$pnts[i,sr1]*x)
    })
    dflux.test <- function(s0, s1) {
      # run wilcox test
      tryCatch({
        wilcox.res <- wilcox.test(s0, s1)
        # p value
        wilcox.p <- wilcox.res$p.value
        # effect size for unpaired test: rank biserial correlation
        wilcox.r <- unname(1 - 2 * wilcox.res$statistic / (length(s0)*length(s1)))
        # another effect size measure: difference of median fluxes
        m0 <- median(s0)
        m1 <- median(s1)
        dm <- m1-m0
        data.table(lb0=min(s0), ub0=max(s0), med0=m0,
                   lb1=min(s1), ub1=max(s1), med1=m1,
                   diff.med=dm, r=wilcox.r, pval=wilcox.p)
      }, error=function(e) {
        data.table(lb0=NA, ub0=NA, med0=NA, lb1=NA, ub1=NA, med1=NA, diff.med=NA, r=NA, pval=NA)
      })
    }
    res <- rbindlist(mclapply(1:length(rxns), function(i) dflux.test(samp0[,i], samp1[,i]), mc.cores=nc))
    res[, padj:=p.adjust(pval, method="BH")]
    res <- cbind(data.table(id=ifelse(is.null(names(rxns)), 1:length(rxns), names(rxns))), res)
    res <- res[order(-abs(diff.med), -abs(r), padj, pval)]
    # add summary of flux differences: 1 means increased flux through the metabolite, vice versa; 0 means unchanged
    res[, dir:=ifelse(!(padj<padj.cutoff & abs(r)>quantile(abs(r), r.cutoff, na.rm=TRUE) & abs(diff.med)>quantile(abs(diff.med), diff.med.cutoff, na.rm=TRUE)), 0, ifelse(r>0, 1, -1))]
  } else {
    ids <- lapply(rxns, function(x) as.integer(names(x)))
    ub0 <- abs(unlist(mcmapply(get.opt.flux, ids, rxns, MoreArgs=list(model=imat.model0, dir="max"), mc.cores=nc)))
    lb0 <- abs(unlist(mcmapply(get.opt.flux, ids, rxns, MoreArgs=list(model=imat.model0, dir="min"), mc.cores=nc)))
    ub1 <- abs(unlist(mcmapply(get.opt.flux, ids, rxns, MoreArgs=list(model=imat.model1, dir="max"), mc.cores=nc)))
    lb1 <- abs(unlist(mcmapply(get.opt.flux, ids, rxns, MoreArgs=list(model=imat.model1, dir="min"), mc.cores=nc)))
    Rcplex.close()
    m0 <- (ub0+lb0)/2
    m1 <- (ub1+lb1)/2
    dm <- m1-m0
    res <- data.table(id=1:length(rxns), rxn=names(rxns),
                      lb0=lb0, ub0=ub0, med0=m0,
                      lb1=lb1, ub1=ub1, med1=m1,
                      diff.med=dm)
    res <- res[order(-abs(diff.med))]
    # add summary of flux differences: positive value means increased flux through the metabolite, vice versa; 0 means unchanged
    res[, dir:=ifelse(ub1>ub0 & lb1>lb0, 3,
               ifelse(ub1<ub0 & lb1<lb0, -3,
               ifelse(ub1>ub0 & lb1==lb0 | ub1==ub0 & lb1>lb0, 2,
               ifelse(ub1<ub0 & lb1==lb0 | ub1==ub0 & lb1<lb0, -2,     
               ifelse(med1>med0, 1,
               ifelse(med1<med0, -1, 0))))))]
  }
}

get.diff.flux.by.met <- function(imat.model0, imat.model1, use.sample=TRUE, sample.range=NULL, mets="all", nc=1L, padj.cutoff=0.01, r.cutoff=0.1, diff.med.cutoff=0.1) {
    # do differential flux analysis for the flux through each metabolite (i.e. either the production or consumption, they should be the same in magnitude if S*v=0; below I am using the production flux): imat.model1 compared to imat.model0
    # if use.sample, use the sampled flux distributions saved in the imat.models, by default, use 1001:end sample points; or specify range of sample points for model0 and model1 respectively in a list
    # is use.sample==FALSE, then deterimine diff flux by obtaining the min and max fluxes through each metabolite for the two models, and regard metabolites with either increased lb or increased ub as upregulated and vice versa; ambiguous cases (e.g. lowered lb and increased ub) are regarded as non-differential
    # by default, mets="all" means performing the analysis across all metabolites; or specify metabolite indices; nc is the number of cores for paralleling across the metabolites
    # padj.cutoff and r.cutoff and log.fc.cutoff are used to determine the significantly changed metabolites when use.sample==TRUE
  if (length(mets)==1 && mets=="all") mets <- 1:length(imat.model0$mets)
  if (use.sample) {
    if (is.null(sample.range)) {
      sr0 <- 1001:ncol(imat.model0$sampl$pnts)
      sr1 <- 1001:ncol(imat.model1$sampl$pnts)
    } else {
      sr0 <- sample.range[[1]]
      sr1 <- sample.range[[2]]
    }
    samp0 <- apply(imat.model0$S[mets,,drop=FALSE], 1, function(x) {
      tmp <- imat.model0$sampl$pnts[,sr0]*x
      apply(tmp, 2, function(y) sum(y[y>0]))
    })
    samp1 <- apply(imat.model1$S[mets,,drop=FALSE], 1, function(x) {
      tmp <- imat.model1$sampl$pnts[,sr1]*x
      apply(tmp, 2, function(y) sum(y[y>0]))
    })
    dflux.test <- function(s0, s1) {
      # run wilcox test
      tryCatch({
        wilcox.res <- wilcox.test(s0, s1)
        # p value
        wilcox.p <- wilcox.res$p.value
        # effect size for unpaired test: rank biserial correlation
        wilcox.r <- unname(1 - 2 * wilcox.res$statistic / (length(s0)*length(s1)))
        # another effect size measure: difference of median fluxes
        m0 <- median(s0)
        m1 <- median(s1)
        dm <- m1-m0
        data.table(lb0=min(s0), ub0=max(s0), med0=m0,
                   lb1=min(s1), ub1=max(s1), med1=m1,
                   diff.med=dm, r=wilcox.r, pval=wilcox.p)
      }, error=function(e) {
        data.table(lb0=NA, ub0=NA, med0=NA, lb1=NA, ub1=NA, med1=NA, diff.med=NA, r=NA, pval=NA)
      })
    }
    res <- rbindlist(mclapply(1:length(mets), function(i) dflux.test(samp0[,i], samp1[,i]), mc.cores=nc))
    res[, padj:=p.adjust(pval, method="BH")]
    res <- cbind(data.table(id=mets, met=imat.model0$mets[mets]), res)
    res <- res[order(-abs(diff.med), -abs(r), padj, pval)]
    # add summary of flux differences: 1 means increased flux through the metabolite, vice versa; 0 means unchanged
    res[, dir:=ifelse(!(padj<padj.cutoff & abs(r)>quantile(abs(r), r.cutoff, na.rm=TRUE) & abs(diff.med)>quantile(abs(diff.med), diff.med.cutoff, na.rm=TRUE)), 0, ifelse(r>0, 1, -1))]
  } else {
    # here there is a problem... I won't be able to know which reactions to use to represent the flux through a metabolite in the cases involving multiple reversible reactions, e.g. (1) x->y, (2) z<=>x, (3) x<=>, here for x the fluxes will sum to 0 if S*v=0, but to represent the flux through x, in general I don't know whether I should use (1)-(2) or (1)+(3) or just (1), so I omit these cases
    tmp <- apply(imat.model0$S[mets,], 1, function(x) sum(model$lb[x!=0]<0)>1)
    if (any(tmp)) warning("Cannot decide the flux through some of the metabolites, these metabolites are omitted.\n")
    mets <- mets[!tmp]
    ids <- apply(imat.model0$S[mets,], 1, function(x) which(x>0))
    coefs <- apply(imat.model0$S[mets,], 1, function(x) x[x>0])
    ub0 <- abs(unlist(mcmapply(get.opt.flux, ids, coefs, MoreArgs=list(model=imat.model0, dir="max"), mc.cores=nc)))
    lb0 <- abs(unlist(mcmapply(get.opt.flux, ids, coefs, MoreArgs=list(model=imat.model0, dir="min"), mc.cores=nc)))
    ub1 <- abs(unlist(mcmapply(get.opt.flux, ids, coefs, MoreArgs=list(model=imat.model1, dir="max"), mc.cores=nc)))
    lb1 <- abs(unlist(mcmapply(get.opt.flux, ids, coefs, MoreArgs=list(model=imat.model1, dir="min"), mc.cores=nc)))
    Rcplex.close()
    m0 <- (ub0+lb0)/2
    m1 <- (ub1+lb1)/2
    dm <- m1-m0
    res <- data.table(id=mets, met=imat.model0$mets[mets],
                      lb0=lb0, ub0=ub0, med0=m0,
                      lb1=lb1, ub1=ub1, med1=m1,
                      diff.med=dm)
    res <- res[order(-abs(diff.med))]
    # add summary of flux differences: positive value means increased flux through the metabolite, vice versa; 0 means unchanged
    res[, dir:=ifelse(ub1>ub0 & lb1>lb0, 3,
               ifelse(ub1<ub0 & lb1<lb0, -3,
               ifelse(ub1>ub0 & lb1==lb0 | ub1==ub0 & lb1>lb0, 2,
               ifelse(ub1<ub0 & lb1==lb0 | ub1==ub0 & lb1<lb0, -2,     
               ifelse(med1>med0, 1,
               ifelse(med1<med0, -1, 0))))))]
  }
}

get.diff.flux.by.met.x2 <- function(imat.model, c0="cell2", c1="cell1", use.sample=TRUE, sample.range=NULL, nc=1L, padj.cutoff=0.01, r.cutoff=0.1, diff.med.cutoff=0.1) {
  # diff.flux.by.met for bicellular model, between the two cells, i.e. c1 (_cell1) vs c0 (_cell2)
  # will just re-use get.diff.flux.by.met; but perform analysis always for all intracellular metabolites
  # in the result, metabolite id correspond to those in the single cellular model
  m0 <- grep(c0, imat.model$mets)
  m1 <- grep(c1, imat.model$mets)
  if (c1=="cell1") met.ids <- m1 else met.ids <- m0
  mets <- stringr::str_sub(imat.model$mets[met.ids],1,-7)
  im0 <- list(mets=mets, S=imat.model$S[m0,], sampl=list(pnts=imat.model$sampl$pnts))
  im1 <- list(mets=mets, S=imat.model$S[m1,], sampl=list(pnts=imat.model$sampl$pnts))
  res <- get.diff.flux.by.met(im0, im1, use.sample, sample.range, mets="all", nc, padj.cutoff, r.cutoff, diff.med.cutoff)
  res[, id:=met.ids[id]]
}

check.diff.flux.of.met <- function(dflux.res, met.ids, model) {
  # given dflux.res and a set of metabolite indices in met.ids, for each metabolite get the diff.flux result of the reactions associated with it
  names(met.ids) <- model$mets[met.ids]
  lapply(met.ids, function(i) {
    tmp <- sign(model$S[i,])
    rxn.ids <- which(tmp!=0)
    dirs <- tmp[rxn.ids]
    tmp <- dflux.res[match(rxn.ids, id), .(diff.med=med1-med0, dir)]
    res <- cbind(data.table(met=model$mets[i], rxn.id=rxn.ids, met.dir=dirs), tmp, data.table(subsystem=model$subSystems[rxn.ids], equation=get.rxn.equation(rxn.ids,model)))
    res[order(-met.dir, -diff.med)]
  })
}

get.dflux.subnetwork <- function(dflux.res, model, dflux.cutoff=1, exclude.mets="^h[[_].\\]?|^oh1[[_].\\]?|^h2o[[_].\\]?|^atp[[_].\\]?|^adp[[_].\\]?|^pi[[_].\\]?|^ppi[[_].\\]?|^coa[[_].\\]?|^o2[[_].\\]?|^co2[[_].\\]?|^nadp[[_].\\]?|^nadph[[_].\\]?|^nad[[_].\\]?|^nadh[[_].\\]?|^fad[[_].\\]?|^fadh2[[_].\\]?|^na1[[_].\\]?|^so4[[_].\\]?|^nh4[[_].\\]?|^cl[[_].\\]?") {
  # from the result of differential flux analysis with get.diff.flux, identify all the subnetworks (of >2 reactions) with a consistent direction of flux difference (i.e. all increase or all decrease).
  # dflux.cutoff: used to determine the reactions with differential fluxes
  # exclude.mets: regex of some high degree metabolites to be excluded; the regex works for mets formats like "h[c]" and "h_c"
  
  library(igraph)
  s <- model$S
  # exclude metabolites
  if (!(is.null(exclude.mets) || exclude.mets=="")) {
    exclude.mets <- grep(exclude.mets, model$mets)
    cat("The following metabolites are excluded:\n")
    cat(paste(model$mets[exclude.mets], collapse=", "), "\n")
    s[exclude.mets, ] <- 0
  }
  # create bipartite graph between mets and rxns
  gb <- graph.incidence(s)
  # project bipartite graph into the graph of rxns
  gp <- bipartite.projection(gb, which="true")
  V(gp)$name <- 1:length(V(gp)) # set vertex names to the rxn ids, since the induced_subgraph below will re-number the vertices and then the rxn ids will be lost
  # for reactions with increased fluxes
  df.rxns <- dflux.res[dir>=dflux.cutoff, id]
  gp1 <- induced_subgraph(gp, df.rxns)
  subn.pos <- groups(components(gp1))
  subn.pos <- cbind(data.table(df.dir=1), rbindlist(lapply(subn.pos, function(x) data.table(size=length(x), rxn.ids=list(as.numeric(x))))))[size>1][order(-size)]
  # for reactions with decreased fluxes
  df.rxns <- dflux.res[-dir>=dflux.cutoff, id]
  gp1 <- induced_subgraph(gp, df.rxns)
  subn.neg <- groups(components(gp1))
  subn.neg <- cbind(data.table(df.dir=-1), rbindlist(lapply(subn.neg, function(x) data.table(size=length(x), rxn.ids=list(as.numeric(x))))))[size>1][order(-size)]
  # combine
  subn <- rbind(subn.pos, subn.neg)
  # add information on subsystems and metabolites
  subn[, rxns:=lapply(rxn.ids, function(x) model$rxns[x])]
  subn[, subsystems:=lapply(rxn.ids, function(x) unique(model$subSystems[x]))]
  subn[, met.ids:=lapply(rxn.ids, function(x) which(Matrix::rowSums(abs(s[,x,drop=FALSE]))!=0))]
  subn[, mets:=lapply(met.ids, function(x) model$mets[x])]
  subn
}

get.flux.diversion <- function(dflux.res, model, dflux.cutoff=1, exclude.mets="^h[[_].\\]?|^oh1[[_].\\]?|^h2o[[_].\\]?|^co2[[_].\\]?|^na1[[_].\\]?|^so4[[_].\\]?|^nh4[[_].\\]?|^cl[[_].\\]?") {
  # from the result of differential flux analysis with get.diff.flux, identify the flux differences of different directions (i.e. include both increase and decrease) associated with "branching point" metabolites (i.e. metabolites involved in >= 3 reactions): these cases reflect the metabolic flux diversion between two conditions.
  # dflux.cutoff: used to determine the reactions with differential fluxes
  # exclude.mets: regex of some high degree metabolites to be excluded; the regex works for mets formats like "h[c]" and "h_c"
  
  df.rxns <- dflux.res[abs(dir)>=dflux.cutoff, id]
  s <- model$S
  s[s!=0] <- 1
  # "branching point" metabolites associated with diff flux
  div.mets <- which(Matrix::rowSums(s)>=3 & Matrix::rowSums(s[,df.rxns,drop=FALSE])!=0)
  # exclude metabolites
  if (!(is.null(exclude.mets) || exclude.mets=="")) {
    exclude.mets <- grep(exclude.mets, model$mets)
    cat("The following metabolites are excluded:\n")
    cat(paste(model$mets[exclude.mets], collapse=", "), "\n")
    div.mets <- setdiff(div.mets, exclude.mets)
  }
  names(div.mets) <- div.mets
  
  res <- rbindlist(lapply(div.mets, function(x) {
    x.df.rxns <- df.rxns[s[x, df.rxns]!=0]
    dflux.res[id %in% x.df.rxns, .(rxn.id=id, rxn, met.dir=sign(model$S[x,id]), df.dir=dir, subsystem=model$subSystems[id], equation=get.rxn.equation(id,model))][uniqueN(sign(df.dir))>1][order(-met.dir, -df.dir)]
  }), idcol="met.id")
  res[, met.id:=as.numeric(met.id)]
  res[, met:=model$mets[met.id]]
  setcolorder(res, c("met.id","met","rxn.id","rxn","met.dir","df.dir","subsystem","equation"))
  res
}

subsystems2gsets <- function(model, by=c("rxn","met"), exclude.mets="^h[[_].\\]?|^oh1[[_].\\]?|^h2o[[_].\\]?|^atp[[_].\\]?|^adp[[_].\\]?|^pi[[_].\\]?|^ppi[[_].\\]?|^coa[[_].\\]?|^o2[[_].\\]?|^co2[[_].\\]?|^nadp[[_].\\]?|^nadph[[_].\\]?|^nad[[_].\\]?|^nadh[[_].\\]?|^fad[[_].\\]?|^fadh2[[_].\\]?|^na1[[_].\\]?|^so4[[_].\\]?|^nh4[[_].\\]?|^cl[[_].\\]?", exclude.mets.degree=nrow(model$S), name="subSystems") {
  # create a list of reaction or metabolite ("by") sets from the "subSystems" field of a metabolic model
  # if by metabolite, exclude.mets are the regex of high-degree mets to exclude, also can set degree cutoff with exclude.mets.degree (metabolites with degree higher than this will be excluded); the union of these two will be excluded
  by <- match.arg(by)
  if (is.null(model[[name]])) stop("subSystems not in model.\n")
  tmp <- data.table(path=model$subSystems, rxn.id=as.character(1:length(model$subSystems)))
  tmp <- tmp[!is.na(path), .(rxn.id=list(rxn.id)), by=path]
  gsets <- tmp$rxn.id
  names(gsets) <- tmp$path
  if (by=="rxn") return(gsets)
  # if by reaction, exclude metabolites
  library(igraph)
  tmp <- graph.incidence(model$S)
  tmp <- bipartite.projection(tmp, which="false") # projected network of metabolites
  deg <- igraph::degree(tmp) # make sure use degree from igraph
  mets.rm <- unique(c(grep(exclude.mets, model$mets), which(deg>exclude.mets.degree)))
  cat("The following metabolites are excluded:\n")
  cat(paste(model$mets[mets.rm], collapse=", "), "\n")
  
  gsets <- lapply(gsets, function(x) {
    mets <- rxns2mets(as.integer(x), model)
    mets <- unique(unlist(mets))
    as.character(setdiff(mets, mets.rm))
  })
}

pathway.gsea <- function(dflux.res, pathways, value.name="r", id.name="id") {
  # metabolic pathway enrichment with gsea, from the result of differential flux analysis with get.diff.flux or get.diff.flux.by.met
  # value.name: the variable name in dflux.res for the measure of flux difference
  # id.name: the variable name in dflux.res for the reaction/metabolite id
  library(fgsea)
  vec <- dflux.res[[value.name]]
  names(vec) <- dflux.res[[id.name]]
  res <- fgsea(pathways, vec, nperm=1e4)
  res <- res[order(padj, pval)]
}


#### ---- visualization ----

plot.model <- function(model, rxn.ids, fluxes=rep(1, length(rxn.ids)), dfluxes=rep(0, length(rxn.ids)), met.ids=1:length(model$mets), exclude.mets="^h[[_].\\]?|^oh1[[_].\\]?|^h2o[[_].\\]?|^atp[[_].\\]?|^adp[[_].\\]?|^pi[[_].\\]?|^ppi[[_].\\]?|^coa[[_].\\]?|^o2[[_].\\]?|^co2[[_].\\]?|^nadp[[_].\\]?|^nadph[[_].\\]?|^nad[[_].\\]?|^nadh[[_].\\]?|^fad[[_].\\]?|^fadh2[[_].\\]?|^na1[[_].\\]?|^so4[[_].\\]?|^nh4[[_].\\]?|^cl[[_].\\]?", dup.mets=exclude.mets, use.flux=c("dflux","flux"), use=c("both","color","width"), cols=c("green4","grey","red3"), sizes=c(0.5,5), layout=c("neato","fdp","dot","circo","twopi"), margins=c(150,150,150,150)) {
  # model: the base metabolic model
  # rxn.ids: IDs of the reactions to plot
  # fluxes: the flux values corresponding to the reactions in rxn.ids
  # dfluxes: the values of flux changes corresponding to the reactions in rxn.ids
  # met.ids: IDs of the metabolites to include in the plot
  # exclude.mets: regex for names of metabolites (as in model$mets) to be excluded from the plots; the default regex works for some of the high-degree mets in recon1 and iMM1415
  # dup.mets: after keeping the mets in mets.ids and excluding those in exclude.mets, for the remaining mets, use dup.mets to specify regex of mets to be plot as separate nodes for each reaction, when they are recurrent in multiple reactions; the default is the same as exclude.mets, so these will be excluded; to duplicate these instead of removing them, set exclude.mets to NULL
  # use.flux: choose to plot fluxes or dfluxes
  # use: use line color, or line width, or both to represent the flux or dflux values
  # cols: a vector of length 3, colors for decreased, unchanged, and increased reactions respectively if plotting dfluxes; if plotting fluxes, the 2nd and 3rd colors will be used for low and high fluxes respectively; if not using colors, the 2nd color will be used for all reactions
  # sizes: a vector of length 2, range of line widths
  # layout: graph layout for hyperdraw::graphLayout
  # margins: plot margins on the up, bottom, left, and right
  library(hypergraph)
  library(hyperdraw)
  library(RColorBrewer)
  
  use.flux <- match.arg(use.flux)
  use <- match.arg(use)
  layout <- match.arg(layout)

  # build hyperedges from rxns
  rxns <- model$rxns[rxn.ids]
  met.ids <- setdiff(met.ids, grep(exclude.mets, model$mets))
  mets <- model$mets[met.ids]
  md.ids <- intersect(met.ids, grep(dup.mets, model$mets))
  hypeds <- lapply(1:length(rxn.ids), function(i) {
    x <- model$S[met.ids, rxn.ids[i]]
    mi <- x!=0
    mets.i <- ifelse(met.ids[mi] %in% md.ids, paste0(mets[mi],i), mets[mi])
    rs <- mets.i[x[mi]<0] # reactants
    ps <- mets.i[x[mi]>0] # products
    if (length(rs)==0) rs <- paste0("EX_",ps)
    if (length(ps)==0) ps <- paste0("EX_",rs)
    # by default, create a hyperedge for a reaction, with arrows pointing from reactants to products:
    res <- DirectedHyperedge(rs, ps, label=rxns[i])
    # when plotting flux, if a reversible reaction is going backwards (i.e. flux<0), draw the arrow in the corresponding direction:
    if (use.flux=="flux" && fluxes[i]<0) res <- DirectedHyperedge(ps, rs, label=rxns[i])
    # when plotting dflux, (only) for the reversible reactions, draw the direction of the arrow according to the direction of dflux (also note: when using color, these will always be plotted in red or the color representing increase):
    if (use.flux=="dflux" && model$lb[rxn.ids][i]<0 && dfluxes[i]<0) res <- DirectedHyperedge(ps, rs, label=rxns[i])
    res
  })
  hypeds <- hypeds[!sapply(hypeds, is.null)]
  # now arrows for reversible reactions have been sorted out, we update flux and dflux values for plotting
  fluxes <- abs(fluxes)
  dfluxes[model$lb[rxn.ids]<0] <- abs(dfluxes[model$lb[rxn.ids]<0])

  # if plotting flux:
  if (use.flux=="flux") {
    v <- fluxes
    # trimming large flux values to facilitate plotting
    v0 <- median(v) + 2*mad(v)
    v[v>v0] <- v0
    # line colors:
    if (use %in% c("color","both")) {
      nc <- uniqueN(c(0,v))
      if (nc==1) {
        cols <- rep(cols[2],length(v))
      } else {
        cols <- colorRampPalette(cols[c(2,3)])(nc)
        bid <- as.numeric(cut(c(0,v), nc))
        cols <- cols[bid[-1]]
      }
    } else cols <- rep(cols[2],length(v))
    # line widths:
    if (use %in% c("width","both")) lwds <- v / (max(v)-min(v)) * diff(sizes) - min(v) + sizes[1] else lwds <- rep(2,length(v))
  }
  # if plotting dflux:
  if (use.flux=="dflux") {
    v <- dfluxes
    # if dfluxes values do not range from -1 and 1, trimm large positive and negative dflux values to facilitate plotting
    if (any(v>1) || any(v< -1)) {
      v0 <- median(v[v>0]) + 2*mad(v[v>0])
      v[v>v0] <- v0
      v0 <- median(v[v<0]) - 2*mad(v[v<0])
      v[v<v0] <- v0
    }
    # line colors:
    if (use %in% c("color","both")) {
      unqv <- unique(c(0,v))
      nc1 <- sum(unqv>=0)
      nc2 <- sum(unqv<=0)
      idx1 <- v>=0
      idx2 <- v<0
      if (nc1==1) {
        c1 <- rep(cols[2], sum(idx1))
      } else {
        c1 <- colorRampPalette(cols[2:3])(nc1)
        bid <- as.numeric(cut(c(0,v[idx1]), nc1))
        c1 <- c1[bid[-1]]
      }
      if (nc2==1) {
        c2 <- rep(cols[2], sum(idx2))
      } else {
        c2 <- colorRampPalette(cols[1:2])(nc2)
        bid <- as.numeric(cut(c(0,v[idx2]), nc2))
        c2 <- c2[bid[-1]]
      }
      cols <- c()
      cols[idx1] <- c1
      cols[idx2] <- c2
    } else cols <- rep(cols[2],length(v))
    # line widths:
    if (use %in% c("width","both")) {
      v <- abs(v)
      lwds <- v / (max(v)-min(v)) * diff(sizes) - min(v) + sizes[1]
    } else lwds <- rep(2,length(v))
  }

  # build graph object
  node.names <- unique(unlist(lapply(hypeds, function(x) c(x@head, x@tail))))
  hg <- Hypergraph(node.names, hypeds)
  testbph <- graphBPH(hg)
  my.graph <- graphLayout(testbph, layoutType=layout)
  # various plot parameters
  nodeDataDefaults(my.graph, "shape") <- "box"
  nodeDataDefaults(my.graph, "margin") <- 'unit(3, "mm")'  
  edgeDataDefaults(my.graph, "lwd") <- 2
  graphDataDefaults(my.graph, "arrowLoc") <- "end"
  # set line widths and colors
  for (i in 1:length(rxn.ids)) {
    rxn <- rxns[i]
    lwd <- as.character(lwds[i])
    col <- cols[i]
    lapply(my.graph@edgeNodeIO$outgoing[[rxn]], function(x) edgeData(my.graph, rxn, x, "lwd") <- lwd)
    lapply(my.graph@edgeNodeIO$incoming[[rxn]], function(x) edgeData(my.graph, x, rxn, "lwd") <- lwd)
    lapply(my.graph@edgeNodeIO$outgoing[[rxn]], function(x) edgeData(my.graph, rxn, x, "color") <- col)
    lapply(my.graph@edgeNodeIO$incoming[[rxn]], function(x) edgeData(my.graph, x, rxn, "color") <- col)
  }
  # plot margins
  my.graph@graph@boundBox@upRight@y <- my.graph@graph@boundBox@upRight@y + margins[1] # top
  my.graph@graph@boundBox@botLeft@y <- my.graph@graph@boundBox@botLeft@y - margins[2] # bottom
  my.graph@graph@boundBox@botLeft@x <- my.graph@graph@boundBox@botLeft@x - margins[3] # left 
  my.graph@graph@boundBox@upRight@x <- my.graph@graph@boundBox@upRight@x + margins[4] # right  
  
  # plot
  plot(my.graph)
  #return(my.graph)
}

check.mets.plot <- function(model, fluxes=rep(1, length(rxn.ids)), dfluxes=rep(0, length(rxn.ids)), met.ids=1:length(model$mets), exclude.mets="^h[[_].\\]?|^oh1[[_].\\]?|^h2o[[_].\\]?|^atp[[_].\\]?|^adp[[_].\\]?|^pi[[_].\\]?|^ppi[[_].\\]?|^coa[[_].\\]?|^o2[[_].\\]?|^co2[[_].\\]?|^nadp[[_].\\]?|^nadph[[_].\\]?|^nad[[_].\\]?|^nadh[[_].\\]?|^fad[[_].\\]?|^fadh2[[_].\\]?|^na1[[_].\\]?|^so4[[_].\\]?|^nh4[[_].\\]?|^cl[[_].\\]?", dup.mets=exclude.mets, use.flux=c("dflux","flux"), use=c("both","color","width"), cols=c("green4","grey","red3"), sizes=c(0.5,5), layout=c("neato","fdp","dot","circo","twopi")) {
  
  library(visNetwork)
  
  s <- model$S
  s[grep(exclude.mets, model$mets), ] <- 0
  rxn.ids <- unique(unlist(mets2rxns(met.ids, model)))
  # network data
  tmp <- lapply(rxn.ids, function(i) {
    from <- which(s[,i]<0)
    to <- which(s[,i]>0)
    # edge info
    ed <- as.data.table(expand.grid(from=from, to=to))[, c("title","arrows","smooth"):=list(model$rxns[i], "to", TRUE)]
    # if only 2 or 3 reactants/products, also pull them together by edges
    if (length(from)==2) ed <- rbind(ed, data.table(from=from[1], to=from[2], title=model$rxns[i], arrows=NA, smooth=TRUE))
    if (length(from)==3) ed <- rbind(ed, as.data.table(expand.grid(from=from, to=from))[c(2,3,6)][, c("title","arrows","smooth"):=list(model$rxns[i], NA, TRUE)])
    if (length(to)==2) ed <- rbind(ed, data.table(from=to[1], to=to[2], title=model$rxns[i], arrows=NA, smooth=TRUE))
    if (length(to)==3) ed <- rbind(ed, as.data.table(expand.grid(from=to, to=to))[c(2,3,6)][, c("title","arrows","smooth"):=list(model$rxns[i], NA, TRUE)])
    # node info
    nd <- data.table(id=c(from,to), label=model$mets[c(from,to)])
    list(ed=ed, nd=nd)
  })
  # collect all edges
  eds <- rbindlist(lapply(tmp, function(x) x$ed))
  # collapse
  # collect all nodes
  nds <- rbindlist(lapply(tmp, function(x) x$nd))
  nds <- unique(nds)
  nds[, group:=ifelse(id %in% met.ids, "main", "adjacent")]
  visNetwork(nds, eds) %>% visIgraphLayout(layout="layout_with_fr")
}

#### ---- data preprocessing for iMAT and MTA ----

get.dflux.for.mta <- function(de.res, topn=Inf, padj.cutoff=1.1, model, discrt=TRUE, na.replace=TRUE, reverse.de=TRUE) {
  # get the directions of reaction changes as input for MTA
  # de.res: DE result as output from de(), at most topn (can be Inf) genes in either direction with fdr<0.1 are selected then mapped to reaction changes, and the resulting vector is returned; the number of changed reactions will be printed
  # NOTE that by default (reverse.de==TRUE), dflux seeks to REVERSE the DE changes!
  # if na.replace, NA's will be replaced by zeros
  # as a rule of thumb, try different topn values so that about 100 changed reactions in either direction are obtained (for MIQP only)
  
  # get gene DE vector
  #ex.genes <- unique(model.data$genes[table(unlist(str_extract_all(model.data$rules, "[1-9][0-9]*")))>10]) # genes mapped to > 10 reactions are excluded, this was in the original MATLAB code but not actually used.
  #de.res <- de.res[id %in% model.data$genes & !id %in% ex.genes]
  de.res <- de.res[id %in% model$genes]
  de.res[, padj:=p.adjust(pval, method="BH")]
  gns.ch <- c(de.res[padj<padj.cutoff & log.fc<0][order(log.fc)][1:min(.N,topn), id], de.res[padj<padj.cutoff & log.fc>0][order(-log.fc)][1:min(.N,topn), id])
  de.res[, df:=ifelse(id %in% gns.ch, log.fc, 0)]
  df <- de.res$df
  names(df) <- de.res$id
  
  # map to reactions
  vec <- exprs2rxns(df, type=1, model=model, discrt=discrt, na.replace=na.replace)
  npos <- sum(vec>0, na.rm=TRUE)
  nneg <- sum(vec<0, na.rm=TRUE)
  if (is.infinite(topn)) {
    cat(sprintf("All the DE genes with p.adj<%g in the model selected, mapped to %d up-regulated reactions and %d down-regulated reactions.\n", padj.cutoff, npos, nneg))
  } else cat(sprintf("At most top %g DE genes with p.adj<%g in the model selected, mapped to %d up-regulated reactions and %d down-regulated reactions.\n", topn, padj.cutoff, npos, nneg))

  if (reverse.de) {
    cat("reverse.de==TRUE, will return dflux that is to reverse the DE changes.\n")
    vec <- -vec
  }
  vec
}

discrt.exprs.for.imat <- function(dat, q.lo=0.25, q.hi=0.75, na.replace=TRUE, f=mean, model) {
  # produce input vector for iMAT:
  # average across all samples in the data into a single vector of expression values of all genes, then select only those genes in the model, then discretize the expression values into low (-1L), medium (0L), and high (1L), missing genes (i.e. model genes that are not in the expression data) will be NA's.

  if (class(dat)=="ExpressionSet") {
    mat <- exprs(dat)
    rownames(mat) <- fData(dat)$Gene.symbol
  } else if (is.matrix(dat)) {
    mat <- dat
  } else if (is.vector(dat)) vec <- dat
  
  if (!exists("vec")) vec <- apply(mat, 1, f, na.rm=TRUE)
  vec <- vec[model$genes]
  na.idx <- is.na(vec)
  cat(sprintf("%d model genes not in the expression data, ", sum(na.idx)))
  qlo <- quantile(vec, q.lo, na.rm=TRUE)
  qhi <- quantile(vec, q.hi, na.rm=TRUE)
  vec <- ifelse(vec<qlo, -1L, ifelse(vec>qhi, 1L, 0L))
  if (na.replace) {
    vec[na.idx] <- 0L
    cat("these gene values are set to 0.\n")
  } else cat("these genes have NA values.\n")
  unname(vec)
}

make.ortho.dflux <- function(seed=0, x) {
  # generate random orthogonal dflux vectors to a given dflux vector, x
  
  # shuffle x
  set.seed(seed); y <- sample(x)
  # current inner product
  k <- sum(x*y)
  # if k happens to be 0, simply return y
  if (k==0) return(y)
  # if k is odd, we randomly swap a pair of y's where y0x!0 and y!0x0, so that we reduce one product-is-zero case; if it happens that no such case exists, we randomly swap a pair of y's where y0x0 and y!0x!0, so that we increase on product-is-zero case
  if (k%%2==1) {
    if (any(y==0 & x!=0)) {
      set.seed(seed); p1 <- sample(which(y==0 & x!=0), 1)
      set.seed(seed); p2 <- sample(which(y!=0 & x==0), 1)
    } else {
      set.seed(seed); p1 <- sample(which(y==0 & x==0), 1)
      set.seed(seed); p2 <- sample(which(y!=0 & x!=0), 1)
    }
    tmp <- y[p1]
    y[p1] <- y[p2]
    y[p2] <- tmp
  }
  k <- sum(x*y)
  # now the k should become even; randomly "flip" k/2 product-non-zero positions in the proper direction
  set.seed(seed); fp <- sample(which(x*y==sign(k)), abs(k)/2)
  y[fp] <- -y[fp]
  # then randomly "flip" equal number of positions of product +1 and -1
  set.seed(seed); nf <- sample(0:(sum(x*y!=0)/2), 1)
  set.seed(seed); fpp <- sample(which(x*y==1), nf)
  set.seed(seed); fpn <- sample(which(x*y==-1), nf)
  fp <- c(fpp, fpn)
  y[fp] <- -y[fp]
  y
}



