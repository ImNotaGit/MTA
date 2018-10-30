library(data.table)
library(stringr)


#### ---- metabolic model utils ----

discrt.genes2rxns <- function(vec, type=0, model, na.replace=TRUE) {
  # map a vector of values -1/0/1 either meant to be levels or direction of changes of genes to those of the reactions in the metabolic model
  # vec is a named vector of -1/0/1, names being gene symbol; or if it's unnamed and length being model$genes, assume it's already in the same order as model$genes
  # type==0 for vec as gene levels, type==1 for vec as directions of changes
  # NA's will be replaced by zeros if na.replace==TRUE

  if (is.null(names(vec)) && length(vec)==length(model$genes)) x <- vec else x <- vec[model$genes]
  if (type==0) {
    `&` <- function(a,b) ifelse(is.na(a) & b==-1 | is.na(b) & a==-1, -1, pmin(a,b)) # if one is NA and the other is -1, for sure the result is -1; all other NA cases are undetermined and NA will be returned
    `|` <- function(a,b) ifelse(is.na(a) & b==1 | is.na(b) & a==1, 1, pmax(a,b)) # if one is NA and the other is 1, for sure the result is 1; all other NA cases are undetermined and NA will be returned
  } else if (type==1) {
    `&` <- function(a,b) ifelse(a!=b | a==0 | b==0, 0, a) # if one is NA and the other is 0, for sure the result is 0; all other NA cases are undetermined and NA will be returned
    `|` <- function(a,b) ifelse(a==b, a, a+b) # all NA cases are undetermined and NA will be returned
  }
  res <- sapply(model$rules, function(i) eval(parse(text=i)))
  if (na.replace) res[is.na(res)] <- 0
  as.integer(unname(res))
}

rxns2genes <- function(vec, model) {
  # given a numerical vector corresponding to the reaction indeces in the model, map each of them to gene names, return as a list; NA will be returned for reaction indeces outside the proper range (including 0)
  vec[vec==0] <- NA
  res <- lapply(str_extract_all(model$rules[vec], "[1-9][0-9]*"), function(x) unique(model$genes[as.numeric(x)]))
  #if (length(res)==1) res <- res[[1]]
  return(res)
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
      which(discrt.genes2rxns(tmp, 0, model)==-1)
    })
  }
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
  head(topTable(fit, coef=coef, number=Inf))
  res <- as.data.table(topTable(fit, coef=coef, number=Inf))
  setnames(res, c("id","log.fc","ave.expr","t","pval","padj","B"))
  res
}

discrt.dflux.for.mta <- function(de.res, topn=Inf, padj.cutoff=1.1, model, na.replace=TRUE, reverse.de=TRUE) {
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
  de.res[, df:=as.integer(sign(log.fc))]
  de.res <- rbind(de.res[padj<padj.cutoff & log.fc<0][order(log.fc)][1:min(.N,topn)], de.res[padj<padj.cutoff & log.fc>0][order(-log.fc)][1:min(.N,topn)])
  df <- de.res$df
  names(df) <- de.res$id
  
  # map to reactions
  vec <- discrt.genes2rxns(df, type=1, model=model, na.replace=na.replace)
  npos <- sum(vec==1L, na.rm=TRUE)
  nneg <- sum(vec==-1L, na.rm=TRUE)
  cat(sprintf("Top %g DE genes in the model selected, mapped to %d up-regulated reactions and %d down-regulated reactions.\n", topn, npos, nneg))

  if (reverse.de) {
    cat("reverse.de==TRUE, will return dflux that is to reverse the DE changes.\n")
    vec <- -vec
  }
  vec
}

discrt.exprs.for.imat <- function(dat, q.lo=0.25, q.hi=0.75, model) {
  # produce input vector for iMAT:
  # average across all samples in the data into a single vector of expression values of all genes, then select only those genes in the model, then discretize the expression values into low (-1L), medium (0L), and high (1L), missing genes (i.e. model genes that are not in the expression data) will be NA's.

  if (class(dat)=="ExpressionSet") {
    mat <- exprs(dat)
    rownames(mat) <- fData(dat)$Gene.symbol
  } else if (is.matrix(dat)) mat <- dat
  
  vec <- rowMeans(mat)
  vec <- vec[model$genes]
  na.idx <- is.na(vec)
  cat(sprintf("%d model genes not in the expression data.\n", sum(na.idx)))
  qlo <- quantile(vec, q.lo, na.rm=TRUE)
  qhi <- quantile(vec, q.hi, na.rm=TRUE)
  vec <- ifelse(vec<qlo, -1L, ifelse(vec>qhi, 1L, 0L))
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

