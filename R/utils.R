library(data.table)
library(stringr)
library(parallel)


#### ---- metabolic model utils ----

library(Rcplex.my)

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
      if (isTRUE(is.na(a) && b<0)) return(b)
      if (isTRUE(is.na(b) && a<0)) return(a) # if one is NA and the other <0, for sure the result is the <0 value; all other NA cases are undetermined and NA will be returned
      return(min(a,b))
    }
    `|` <- function(a,b) {
      if (isTRUE(is.na(a) && b>0)) return(b)
      if (isTRUE(is.na(b) && a>0)) return(a) # if one is NA and the other >0, for sure the result is the >0 value; all other NA cases are undetermined and NA will be returned
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

get.rxn.equation <- function(vec, model) {
  # given a numerical vector corresponding to the reaction indeces in the model, return string vector containing the equations of the corresponding reactions
  sapply(vec, function(i) {
    x <- model$S[, i]
    rs <- paste(model$mets[x<0], collapse=" + ")
    ps <- paste(model$mets[x>0], collapse=" + ")
    if (model$lb[i]>=0) arrow <- "-->" else arrow <- "<==>"
    paste(rs, arrow, ps)
  })
}

get.opt.flux <- function(model, i, coef=1, dir="max", ko=NULL, keep.xopt=FALSE, nc=1L) {
  # get the max or min flux of the i'th reaction in the model
  # `i` can also be a vector of multiple reaction indices, with coef being their coefficients, then the corresponding linear objective function will be optimized
  # to knockout reaction(s), pass KO as their indices
  # if keep.xopt, the optimal xopt vector will also be returned
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
  if (res$status!=1) warning("Potential problems running LP. Solver status: ", res$status, ".\n")
  if (keep.xopt) {
    return(list(obj=res$obj, xopt=res$xopt))
  } else return(res$obj)
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


#### ---- visualization ----

library(hypergraph)
library(hyperdraw)
library(RColorBrewer)

plot.model <- function(model, rxn.ids=1:length(model$rxns), fluxes=rep(1, length(rxns)), dfluxes=rep(0, length(rxns)), met.ids=1:length(model$mets), mets.dup=c(), layout="dot", mode=0, lwd.rng=c(0.5,10), cols=c("green4","grey","red3"), plt.margins=c(150,150,150,150)) {
  # model: the base metabolic model
  # rxn.ids: IDs of the reactions to plot
  # fluxes and dfluxes: the reference fluxes and flux changes of the reactions in rxn.ids
  # met.ids: IDs of the metabolites to plot
  # mets.dup: names of metabolites (as in model$mets) to be plot as separate nodes for each reaction, when they are recurrent in multiple reactions
  # layout
  # mode 0: line width represents dfluxes, line color represents the direction of dfluxes; mode 1: line width represents fluxes, and line color represents dfluxes
  # lwd.rng: range of line widths to use
  # cols: color values for decreased, unchanged, and increased fluxes (decrease/increase in terms of absolute value)
  # plt.margins: plot margins on the up, bottome, left, and right
  
  # build hyperedges from rxns
  rxns <- model$rxns[rxn.ids]
  mets <- model$mets[met.ids]
  hyperedges <- lapply(1:length(rxn.ids), function(i) {
    rxn.id <- rxn.ids[i]
    rxn <- rxns[i]
    flux <- fluxes[i]
    x <- model$S[met.ids, rxn.id]
    reactants <- mets[x<0]
    mdup <- reactants %in% mets.dup
    reactants[mdup] <- paste0(reactants[mdup], i)
    products <- mets[x>0]
    mdup <- products %in% mets.dup
    products[mdup] <- paste0(products[mdup], i)
    if (flux>=0) {
      res <- DirectedHyperedge(reactants, products, label=rxn)
    } else {
      res <- DirectedHyperedge(products, reactants, label=rxn)
    }
    res
  })
  
  # build graph object
  node.names <- unique(unlist(lapply(hyperedges, function(x) c(x@head, x@tail))))
  hg <- Hypergraph(node.names, hyperedges)
  testbph <- graphBPH(hg)
  my.graph <- graphLayout(testbph, layoutType=layout)
  
  nodeDataDefaults(my.graph, "shape") <- "box"
  nodeDataDefaults(my.graph, "margin") <- 'unit(3, "mm")'  
  edgeDataDefaults(my.graph, "lwd") <- 1
  graphDataDefaults(my.graph, "arrowLoc") <- "end"
  
  df <- ifelse(sign(fluxes)*sign(dfluxes)>=0, abs(dfluxes), -abs(dfluxes)) / abs(fluxes)
  df[is.nan(df)] <- 0
  df.med <- median(df)
  df.3mad <- 3*median(abs(df-df.med))
  df.ub <- df.med+df.3mad
  df.lb <- df.med-df.3mad
  df[df>df.ub] <- df.ub
  df[df<df.lb] <- df.lb
  if (mode==0) {
    # line widths represents dfluxes, line color means the direction of dfluxes
    ## line widths
    lwds <- abs(df) / median(abs(df)) * (lwd.rng[1]+lwd.rng[2])/2
    lwds[is.nan(lwds)] <- (lwd.rng[1]+lwd.rng[2])/2
    lwds[lwds<lwd.rng[1]] <- lwd.rng[1]
    lwds[lwds>lwd.rng[2]] <- lwd.rng[2]
    ## line colors
    cols <- ifelse(df>0, cols[3], ifelse(df<0, cols[1], cols[2]))
  } else if (mode==1) {
    # line widths means fluxes, and line colors means dfluxes
    ## line widths
    lwds <- abs(fluxes) / max(abs(fluxes)) * lwd.rng[2]
    lwds[lwds<lwd.rng[1]] <- lwd.rng[1]
    ## line colors
    if (all(dfluxes==0)) {
      cols <- rep(cols[2], length(dfluxes))
    } else {
      ncols <- max(sum(df>0), sum(df<0))
      cols <- colorRampPalette(cols)(2*ncols+3)
      tmp <- as.numeric(cut(c(0, df), ncols+1))
      delt <- ncols+2 - tmp[1]
      cols <- cols[tmp[-1]+delt]
    }
  }
  # set line widths and colors
  for (i in 1:length(rxns)) {
    rxn <- rxns[i]
    lwd <- lwds[i]
    col <- cols[i]
    lapply(my.graph@edgeNodeIO$outgoing[[rxn]], function(x) edgeData(my.graph, rxn, x, "lwd") <- as.character(lwd))
    lapply(my.graph@edgeNodeIO$incoming[[rxn]], function(x) edgeData(my.graph, x, rxn, "lwd") <- as.character(lwd))
    lapply(my.graph@edgeNodeIO$outgoing[[rxn]], function(x) edgeData(my.graph, rxn, x, "color") <- col)
    lapply(my.graph@edgeNodeIO$incoming[[rxn]], function(x) edgeData(my.graph, x, rxn, "color") <- col)
  }
  
  # set plot margins
  my.graph@graph@boundBox@botLeft@y <- my.graph@graph@boundBox@botLeft@y-plt.margins[2] #bottom
  my.graph@graph@boundBox@botLeft@x <- my.graph@graph@boundBox@botLeft@x-plt.margins[3] #left 
  my.graph@graph@boundBox@upRight@y <- my.graph@graph@boundBox@upRight@y+plt.margins[1] #top
  my.graph@graph@boundBox@upRight@x <- my.graph@graph@boundBox@upRight@x+plt.margins[4] #right  
  
  plot(my.graph)
  
  return(my.graph)
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

discrt.exprs.for.imat <- function(dat, q.lo=0.25, q.hi=0.75, na.replace=TRUE, model) {
  # produce input vector for iMAT:
  # average across all samples in the data into a single vector of expression values of all genes, then select only those genes in the model, then discretize the expression values into low (-1L), medium (0L), and high (1L), missing genes (i.e. model genes that are not in the expression data) will be NA's.

  if (class(dat)=="ExpressionSet") {
    mat <- exprs(dat)
    rownames(mat) <- fData(dat)$Gene.symbol
  } else if (is.matrix(dat)) mat <- dat
  
  vec <- rowMeans(mat)
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



