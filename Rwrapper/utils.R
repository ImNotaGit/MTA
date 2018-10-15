library(data.table)
library(stringr)

user <- system("whoami", intern=TRUE)
if (user=="fountain") {
  mta.path <- "/home/fountain/Documents/Projects/ourtools/MTA"
  cplex.path <- "/home/fountain/ibm/ILOG/CPLEX_Studio128/cplex/matlab/x86-64_linux"
  tomlab.path <- "/home/fountain/tomlab"
  gurobi.path <- "/home/fountain/gurobi800/linux64/matlab"
  cobra.path <- "/home/fountain/biocomp_tools/cobratoolbox"
} else if (user=="kycheng") {
  mta.path <- "/fs/cbcb-scratch/kycheng/MTA"
  cplex.path <- "/cbcb/project2-scratch/kycheng/tools/ibm/ILOG/CPLEX_Studio128/cplex/matlab/x86-64_linux"
  tomlab.path <- "/cbcb/project2-scratch/kycheng/tools/tomlab"
  gurobi.path <- "/cbcb/project2-scratch/kycheng/tools/gurobi800/linux64/matlab"
  cobra.path <- "/cbcb/project2-scratch/kycheng/tools/cobratoolbox"
}
rm(user)


# Recon1 model data
#library(R.matlab)
#recon <- readMat("/home/fountain/Documents/Projects/ourtools/MTA/models/recon.mat")
#genes <- unname(unlist(recon$model[[24]])) # unique gene symbols
#genes.map <- as.vector(recon$model[[22]]) # mapping between all model genes and unique genes
#model.gene.symbols <- genes[genes.map]
#model.rxn.rules <- sapply(recon$model[[8]], function(x) { # rules mapping genes to reactions
#  x <- unlist(x)
#  if (length(x)==0) x <- "0"
#  str_replace_all(x, "x\\(([0-9]+)\\)", "x\\[\\1\\]")
#})
#model <- list(genes=model.gene.symbols, rules=model.rxn.rules)
#save(model, file="recon.data.RData")

# ElegCyc model data
#load("/home/fountain/Documents/Projects/ourtools/MTA/R/ElegCyc.RData")
#model <- elegcyc[c("genes","rules")]
#save(model, file="elegcyc.data.RData")


use.model <- function(model="recon") {
  if (model=="recon") {
    load(file.path(mta.path, "Rwrapper", "recon.data.RData"), envir=.GlobalEnv)
    model.mat <<- file.path(mta.path, "models", "recon.mat")
  } else if (model=="elegcyc") {
    load(file.path(mta.path, "Rwrapper", "elegcyc.data.RData"), envir=.GlobalEnv)
    model.mat <<- file.path(mta.path, "models", "elegcyc.mat")
  }
}


#### ---- matlab utils ----

write.mat.vec <- function(vec, mat.file, mat.var, by.row=FALSE) {
  # write a vector to .mat file, with the matlab variable (object) name being mat.var
  # by.row: whether to write as a row vector (TRUE) or a column vector (FALSE)
  library(R.matlab)

  if (by.row) {
    dim(vec) <- c(1, length(vec))
  } else vec <- as.vector(vec)
  vec <- unname(vec)
  
  if (!endsWith(mat.file, ".mat")) mat.file <- paste0(mat.file, ".mat")
  tmp <- list(mat.file, vec)
  names(tmp) <- c("con", mat.var)
  do.call(writeMat, tmp)
  cat(sprintf("Written to %s.\n", mat.file))
}

exists.matlab.var <- function(server=matlab, var, exit=FALSE) {
  # check whether a variable exists in the MATLAB workspace; when exit=TRUE and the variable does not exist, stop.
  
  if (!isOpen(server)) stop("The MATLAB server is not open.\n")
  evaluate(server, paste0("tmp = ismember('", var, "', who);"))
  tmp <- getVariable(server, "tmp")
  evaluate(server, "clear tmp")
  res <- as.logical(tmp[[1]])
  if (!res && exit) stop(sprintf("Variable by the name %s not found.\n", var))
  return(res)
}

#### ---- metabolic model utils ----

discrt.genes2rxns <- function(vec, type="stat", model.data=model) {
  # map a vector of values -1/0/1 either meant to be levels or direction of changes of genes to those of the reactions in the metabolic model
  # vec is a named vector of -1/0/1, names being gene symbol
  # type=="stat" for vec as gene levels, type=="diff" for vec as directions of changes

  x <- vec[model.data$genes]
  if (type=="stat") {
    `&` <- min
    `|` <- max
  } else if (type=="diff") {
    `&` <- function(a,b) ifelse(a==b, a, 0)
    `|` <- function(a,b) ifelse(a==b, a, a+b)
  }
  res <- sapply(model.data$rules, function(i) eval(parse(text=i)))
  res[is.na(res)] <- 0 # NA's will be replaced by zeros!
  as.integer(unname(res))
}

rxns2genes <- function(vec, model.data=model) {
  # given a numerical vector corresponding to the reaction indeces in the model, map each of them to gene names, return as a list (or a simple vector for a single reaction); NA will be returned for reaction indeces outside the proper range (including 0)
  vec[vec==0] <- NA
  res <- lapply(str_extract_all(model.data$rules[vec], "[1-9][0-9]*"), function(x) unique(model.data$genes[as.numeric(x)]))
  if (length(res)==1) res <- res[[1]]
  return(res)
}


#### ---- data preprocessing and matlab input file generating ----

library(GEOquery)

prep.data <- function(dat, log="default", norm.method="loess") {
  # prepare expression data: transformation, normalization, etc.
  library(affy)
  
  if (class(dat)=="ExpressionSet") {
    # for featureData, fix potential improper column names so that later limma::topTable can use them
    fvarLabels(dat) <- make.names(fvarLabels(dat))
    mat <- exprs(dat)
  } else if (is.matrix(dat)) mat <- dat
  
  # log2 transform
  if (log=="default") {
    # following the method as in GEO2R, for microarray data
    qx <- as.numeric(quantile(mat, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm=TRUE))
    log <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  }
  if (log) {
    nlt0 <- sum(mat<0)
    if (nlt0>0) {
      warning(sprintf("Log-transformation error: there are %d negative values in the data, the data may be already on log-scale.\nHere is a summary of the data:\n", nlt0))
      summary(as.vector(mat))
      stop()
    }
    mat <- log2(mat+1)
    cat("log2-transformation performed.\n")
  } else cat("log2-transformation NOT performed.\n")

  # normalization
  if (norm.method=="loess") {
    nna <- sum(is.na(mat))
    if (nna>0) {
      stop(sprintf("Loess normalization error: there are %d NA/NaN's in the data.\n", nna))
    } else mat <- normalize.loess(mat, log.it=FALSE)
  } else if (norm.method=="quantile") {
    mat <- limma::normalizeQuantiles(mat)
  } else cat("Normalization NOT performed.\n")
  
  # return
  if (class(dat)=="ExpressionSet") {
    exprs(dat) <- mat
    return(dat)
  } else if (is.matrix(dat)) return(mat)
}

de <- function(dat, pheno, model="~.", coef, robust=FALSE, trend=FALSE) {
  # differential expression analysis
  library(limma)
  
  if (class(dat)=="ExpressionSet") {
    mat <- exprs(dat)
    rownames(mat) <- fData(dat)$Gene.symbol
  } else if (is.matrix(dat)) mat <- dat
  
  design <- model.matrix(as.formula(model), pheno)
  fit <- lmFit(mat, design)
  fit <- eBayes(fit, robust=robust, trend=trend)
  coef <- as.data.table(topTable(fit, coef=coef, number=Inf, sort.by="none"))
  if (length(coef)==6) coef <- cbind(data.table(id=rownames(mat)), coef)
  setnames(coef, c("id","log.fc","ave.expr","t","pval","padj","B"))
  coef[order(-B)]
}

discrt.dflux.for.mta <- function(de.res, topn=100, padj.cutoff=0.1, model.data=model) {
  # get the directions of reaction changes as input for MTA
  # de.res: DE result as output from de(), at most topn genes in either direction with fdr<0.1 are selected then mapped to reaction changes, and the resulting vector is returned; the number of changed reactions will be printed
  # NOTE that by default, dflux seeks to REVERSE the DE changes!
  # as a rule of thumb, try different topn values so that about 100 changed reactions in either direction are obtained
  
  # get gene DE vector
  #ex.genes <- unique(model.data$genes[table(unlist(str_extract_all(model.data$rules, "[1-9][0-9]*")))>10]) # genes mapped to > 10 reactions are excluded, this was in the original MATLAB code but not actually used.
  #de.res <- de.res[id %in% model.data$genes & !id %in% ex.genes]
  de.res <- de.res[id %in% model.data$genes]
  de.res[, padj:=p.adjust(pval, method="BH")]
  de.res[, df:=as.integer(sign(log.fc))]
  de.res <- rbind(de.res[padj<padj.cutoff & log.fc<0][order(log.fc)][1:min(.N,topn)], de.res[padj<padj.cutoff & log.fc>0][order(-log.fc)][1:min(.N,topn)])
  df <- de.res$df
  names(df) <- de.res$id
  
  # map to reactions
  vec <- discrt.genes2rxns(df, type="diff", model.data=model.data)
  npos <- sum(vec==1L)
  nneg <- sum(vec==-1L)
  cat(sprintf("Top %d DE genes in the model selected, mapped to %d up-regulated reactions and %d down-regulated reactions.\n", topn, npos, nneg))
  
  return(-vec)
}

discrt.exprs.for.imat <- function(dat, q.lo=0.25, q.hi=0.75, model.data=model) {
  # produce input vector for iMAT:
  # average across all samples in the data into a single vector of expression values of all genes, then select only those genes in the model, then discretize the expression values into low (-1L), medium (0L), and high (1L), with missing genes (i.e. model genes that are not in the expression data) set to 0L.
  
  if (class(dat)=="ExpressionSet") {
    mat <- exprs(dat)
    rownames(mat) <- fData(dat)$Gene.symbol
  } else if (is.matrix(dat)) mat <- dat
  
  vec <- rowMeans(mat)
  vec <- vec[model.data$genes]
  na.idx <- is.na(vec)
  cat(sprintf("%d model genes not in the expression data.\n", sum(na.idx)))
  qlo <- quantile(vec, q.lo, na.rm=TRUE)
  qhi <- quantile(vec, q.hi, na.rm=TRUE)
  vec <- ifelse(vec<qlo, -1L, ifelse(vec>qhi, 1L, 0L))
  vec[na.idx] <- 0L # NA's will be replaced by zeros!
  unname(vec)
}


#### ---- R wrapper of matlab MTA ----

library(R.matlab)

start.matlab <- function(name="matlab", verbose=-2, ...) {
  # start and connect to matlab server
  
  Matlab$startServer(...)
  assign(name, Matlab(), envir=.GlobalEnv)
  setVerbose(get(name), verbose)
  is.open <- open(get(name))
  if (!is.open) stop("Failed to connect to the MATLAB server.\n")
  #print(get(name))
}

stop.matlab <- function(server=matlab, save.log=NULL) {
  # close matlab server; save the workspace if save.log is a character string as the file name
  
  if (!is.null(save.log)) {
    if (!endsWith(save.log, ".mat")) save.log <- paste0(save.log, ".mat")
    evaluate(server, paste("save", save.log))
    cat(sprintf("MATLAB workspace saved to %s.\n", save.log))
  }
  close(server)
  rm(list=as.character(substitute(server)), envir=.GlobalEnv)
}

init.matlab.for.mta <- function(server=matlab, cplex=cplex.path, tomlab=tomlab.path, gurobi=gurobi.path, cobra=NA, mta=mta.path) {
  # initialize the paths in MATLAB for running MTA; to disable one solver, set to NA; cobratoolbox by default disabled, set it as cobra=cobra.path to enable it.
  
  if (!isOpen(server)) stop("The MATLAB server is not open.\n")
  # add cplex path
  if (!is.na(cplex)) evaluate(server, paste("addpath", cplex))
  # initiate tomlab
  if (!is.na(tomlab)) {
    evaluate(server, paste("cd", tomlab))
    evaluate(server, "startup")
  }
  # initiate gurobi
  if (!is.na(gurobi)) {
    evaluate(server, paste("cd", gurobi))
    evaluate(server, "gurobi_setup")
  }
  # initiate cobratoolbox
  if (!is.na(cobra)) {
    evaluate(server, paste("cd", cobra))
    evaluate(server, "initCobraToolbox")
  }
  # add MTA paths
  evaluate(server, paste("addpath", file.path(mta, "matlab", "iMAT")))
  evaluate(server, paste("addpath", file.path(mta, "matlab", "MTA")))
  evaluate(server, paste("addpath", file.path(mta, "matlab", "models")))
  # return to current dir
  evaluate(server, paste("cd", getwd()))
}

load.model <- function(model=model.mat, var="model", server=matlab) {
  # load metabolic model into MATLAB
  
  if (!isOpen(server)) stop("The MATLAB server is not open.\n")
  evaluate(server, paste0(var, "=load('", model, "');"))
  evaluate(server, sprintf("tmp=fieldnames(%s); %s=%s.(tmp{1}); clear tmp", var,var,var))
  cat(sprintf("Loaded %s into MATLAB as variable: %s\n", model, var))
}

imat <- function(expr, vout1="imat_v", vout2="imat_sampl_pnts", model="model", milp.solver="cplex", lp.solver="cplex", server=matlab) {
  # run iMAT in MATLAB
  # expr is the discretized gene expression vector as iMAT input; if it's given as character, then the variable with the same name in MATLAB is used; otherwize it should be an R object generated by the function discrt.exprs.for.imat, and is passed into MATLAB
  # by default, the MATLAB variable name for the metabolic model is called "model"
  # specify the MATLAB variable names for the two outputs in vout1 and vout2
  
  if (!isOpen(server)) stop("The MATLAB server is not open.\n")
  exists.matlab.var(server=server, model, exit=TRUE)
  
  if (is.character(expr)) {
    exists.matlab.var(server=server, expr, exit=TRUE)
    expr.n <- expr
  } else {
    expr.n <- str_replace_all(as.character(substitute(expr)), "\\.", "_") # the variable name to be used in MATLAB
    tmp <- list(server, expr)
    names(tmp) <- c("this", expr.n)
    do.call(setVariable, tmp)
  }
  
  # run iMAT
  evaluate(server, sprintf("[%s, %s] = sampleiMAT(%s, %s, '%s', '%s');", vout2, vout1, model, expr.n, milp.solver, lp.solver))
  
  # retrieve results from MATLAB, need to do it separately otherwise the order will not keep as is.
  res1 <- getVariable(server, vout1)
  res2 <- getVariable(server, vout2)
  list(v=as.vector(res1[[1]]), sampl.pnts=res2[[1]])
}

mta <- function(v.ref, dflux, del=NULL, vout="mta_res", model="model", solver="cplex", server=matlab) {
  # run MTA in MATLAB
  # usage on the arguments is similar to the function imat
  # v.ref: source state fluxes; dflux: intented flux changes; del: the reaction-deletion(s) to screen
  # by default (when del is NULL), deletions of all reactions in the model will be screened, in addition to the control (i.e. no reaction is deleted, represented by del=0)
  
  if (!isOpen(server)) stop("The MATLAB server is not open.\n")
  exists.matlab.var(server=server, model, exit=TRUE)
  
  if (is.character(v.ref)) {
    exists.matlab.var(server=server, v.ref, exit=TRUE)
    v.ref.n <- v.ref
  } else {
    v.ref.n <- str_replace_all(as.character(substitute(v.ref)), "\\.", "_")
    tmp <- list(server, v.ref)
    names(tmp) <- c("this", v.ref.n)
    do.call(setVariable, tmp)
  }
  
  if (is.character(dflux)) {
    exists.matlab.var(server=server, dflux, exit=TRUE)
    dflux.n <- dflux
  } else {
    dflux.n <- str_replace_all(as.character(substitute(dflux)), "\\.", "_")
    tmp <- list(server, dflux)
    names(tmp) <- c("this", dflux.n)
    do.call(setVariable, tmp)
  }
  
  if (is.null(del)) {
    evaluate(server, sprintf("del = 0:length(%s.rxns);", model)) # 0 and all reactions
    del.n <- "del" # use the variable name "del" in MATLAB
  } else if (is.character(del)) {
    exists.matlab.var(server=server, del, exit=TRUE)
    del.n <- del
  } else {
    tmp <- substitute(del)
    del.n <- ifelse(is.numeric(tmp) || length(as.character(tmp))>1, "del", str_replace_all(as.character(tmp), "\\.", "_")) # if del not assigned from another variable but directly as expression, use the variable name "del" in MATLAB
    tmp <- list(server, del)
    names(tmp) <- c("this", del.n)
    do.call(setVariable, tmp)
  }
  
  # run MTA
  evaluate(server, sprintf("%s = MTA(%s, %s, %s, %s, '%s');", vout, model, v.ref.n, dflux.n, del.n, solver))

  # retrieve results from MATLAB
  tmp <- getVariable(server, vout)[[1]]
  res <- rbindlist(apply(tmp, 3, function(x) {
    x <- as.data.table(lapply(x, function(i) {
      i <- as.vector(i)
      if (length(i)>1) i <- list(i)
      i
    }))
    setnames(x, dimnames(tmp)[[1]])
    x
  }))
  return(res)
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

