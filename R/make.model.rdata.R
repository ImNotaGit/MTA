library(R.matlab)
library(stringr)
library(rlist)
library(sybilSBML)


### ---- Recon 1 ----

# this file is from Noam, probably contains an old version of the model
tmp <- readMat("../models/recon.mat")

# tmp[[1]] contains the model but as an array, each element of the array being a field of the MATLAB struct... so need to use apply
# for some weird reason each element of the array is a deeply nested list, so use rlist::list.flatten to make each of them a simple list of only one level
# cannot use unlist here because it can flattern things directly to vectors and discarding NULLs/empty values.
recon1 <- apply(tmp[[1]], 1, function(x) list.flatten(x))
# then manually unlist each of the simple lists, so that NULLs/empty values are not lost
recon1 <- lapply(recon1, function(x) {
  # if of length 1, just extract x[[1]], simplify to vector when appropriate
  if (length(x)==1) {
    res <- x[[1]]
    if (is.matrix(res) && (nrow(res)==1 || ncol(res)==1)) res <- as.vector(res)
  } else { # if length >1, carefully unlist into a vector so that NULLs/empty values are not lost
    res <- unname(sapply(x, function(y) {
      y <- as.vector(y) # y can be a (single-element or empty) matrix/array, so as.vector
      if (length(y)==0 || is.null(y) || y=="") y <- NA
      y
    }))
  }
  return(res)
})

# change the "genes" field
recon1$gene.ids <- recon1$genes
recon1$genes <- recon1$genes.unique.names[recon1$genes.unique.map]
recon1$genes.unique <- NULL
recon1$genes.unique.map <- NULL
recon1$genes.unique.names <- NULL
# 2019.11.12: I find that there are some non-standard gene symbols (maybe outdated symbols?), e.g. ASS1 gene is written as "ASS", so I recreated gene symbols from gene ids
gid <- str_split(recon1$gene.ids, "\\.", simplify=TRUE)[,1]
library("org.Hs.eg.db")
mapp <- select(org.Hs.eg.db, keys=gid, columns=c("ENTREZID","SYMBOL"), keytype="ENTREZID") # many:1 mapping (with NA's)
all(gid==mapp$ENTREZID) # TRUE
recon1$genes <- ifelse(is.na(mapp$SYMBOL), recon1$gene.ids, mapp$SYMBOL)

# format the "rules" field
recon1$rules <- unname(sapply(recon1$rules, function(x) {
  if (is.na(x)) x <- "0"
  str_replace_all(x, "x\\(([0-9]+)\\)", "x\\[\\1\\]")
}))

# change model description
recon1$description <- "Recon 1"

save(recon1, file="Recon1.RData")


# this is another version of recon 1, downloaded from BiGG (2019.11.12)
tmp <- readMat("../models/BiGG/RECON1.mat")

# tmp[[1]] contains the model but as an array, each element of the array being a field of the MATLAB struct... so need to use apply
# for some weird reason each element of the array is a deeply nested list, so use rlist::list.flatten to make each of them a simple list of only one level
# cannot use unlist here because it can flattern things directly to vectors and discarding NULLs/empty values.
recon1 <- apply(tmp[[1]], 1, function(x) list.flatten(x))
# then manually unlist each of the simple lists, so that NULLs/empty values are not lost
recon1 <- lapply(recon1, function(x) {
  # if of length 1, just extract x[[1]], simplify to vector when appropriate
  if (length(x)==1) {
    res <- x[[1]]
    if (is.matrix(res) && (nrow(res)==1 || ncol(res)==1)) res <- as.vector(res)
  } else { # if length >1, carefully unlist into a vector so that NULLs/empty values are not lost
    res <- unname(sapply(x, function(y) {
      y <- as.vector(y) # y can be a (single-element or empty) matrix/array, so as.vector
      if (length(y)==0 || is.null(y) || y=="") y <- NA
      y
    }))
  }
  return(res)
})

# gene symbols
recon1$gene.ids <- recon1$genes
gid <- str_split(recon1$gene.ids, "_", simplify=TRUE)[,1]
library("org.Hs.eg.db")
mapp <- select(org.Hs.eg.db, keys=gid, columns=c("ENTREZID","SYMBOL"), keytype="ENTREZID") # many:1 mapping (with NA's)
all(gid==mapp$ENTREZID) # TRUE
recon1$genes <- ifelse(is.na(mapp$SYMBOL), recon1$gene.ids, mapp$SYMBOL)

# rules mapping genes to reactions
rules <- sapply(recon1$grRules, function(x) {
  if (is.na(x)) {
    x <- "0"
  } else {
    x <- str_replace_all(x, "[0-9]+_AT[0-9]+", function(x) paste0("x[", match(x, recon1$gene.ids), "]"))
    x <- str_replace_all(x, "and", "&")
    x <- str_replace_all(x, "or", "\\|")
  }
  x
})
recon1$rules <- unname(rules)
recon1$rowlb <- recon1$b # all 0
recon1$rowub <- recon1$b # all 0
recon1$S <- Matrix(recon1$S, sparse=TRUE)
recon1$lb[recon1$lb==-999999] <- -1000
recon1$ub[recon1$ub==999999] <- 1000

save(recon1, file="Recon1.BiGG.20191112.RData")


### ---- Recon 2.2 ----

validateSBMLdocument("../models/Recon2.2.xml")
tmp <- readSBMLmod("../models/Recon2.2.xml")

recon2.2 <- list()
recon2.2$description <- tmp@mod_desc
recon2.2$modelName <- tmp@mod_name
recon2.2$modelVersion <- tmp@version
recon2.2$comps <- tmp@mod_compart
recon2.2$mets <- tmp@met_id
recon2.2$metNames <- tmp@met_name
recon2.2$metCharges <- tmp@met_attr$charge
recon2.2$metFormulas <- tmp@met_attr$chemicalFormula
recon2.2$rxns <- tmp@react_id
recon2.2$rxnNames <- tmp@react_name
recon2.2$S <- tmp@S
recon2.2$lb <- tmp@lowbnd
recon2.2$ub <- tmp@uppbnd
recon2.2$c <- tmp@obj_coef
recon2.2$rowlb <- rep(0, length(tmp@met_id))
recon2.2$rowub <- rep(0, length(tmp@met_id))
recon2.2$b <- rep(0, length(tmp@met_id))
recon2.2$rxnGeneMat <- tmp@rxnGeneMat
recon2.2$subSystems <- tmp@subSys
recon2.2$rules <- tmp@gprRules
# gene symbols
library(biomaRt)
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
tmp1 <- getBM(attributes=c("hgnc_id", "hgnc_symbol"), filters="hgnc_id", values=tmp@allGenes, mart=mart)
anyDuplicated(tmp1$hgnc_id) # 0
setdiff(tmp@allGenes, tmp1$hgnc_id) # some bad IDs, and others somehow not matched by biomart, manually map them to gene symbols based on the HGNC database
tmp1 <- rbind(tmp1, data.frame(hgnc_id=c("HGNC:HGNC:987","HGNC:HGNC:2898","HGNC:8780","HGNC:8789","HGNC:8790","HGNC:8903","HGNC:14423","HGNC:11240"), hgnc_symbol=c("BCKDHB","DLD","PDE4A","PDE6G","PDE6H","PGLS","RDH8","SPHK1")))
all(tmp@allGenes %in% tmp1$hgnc_id) # now TRUE
tmp1 <- as.data.table(tmp1)
recon2.2$genes <- tmp1[match(tmp@allGenes, hgnc_id), hgnc_symbol]

save(recon2.2, file="recon2.2.RData")


### ---- Recon 3D ----

tmp <- readMat("../models/Recon3D_301/Recon3DModel_301.mat")

# tmp[[1]] contains the model but as an array, each element of the array being a field of the MATLAB struct... so need to use apply
# for some weird reason each element of the array is a deeply nested list, so use rlist::list.flatten to make each of them a simple list of only one level
# cannot use unlist here because it can flattern things directly to vectors and discarding NULLs/empty values.
recon3d <- apply(tmp[[1]], 1, function(x) list.flatten(x))
# then manually unlist each of the simple lists, so that NULLs/empty values are not lost
recon3d <- lapply(recon3d, function(x) {
  # if of length 1, just extract x[[1]], simplify to vector when appropriate
  if (length(x)==1) {
    res <- x[[1]]
    if (is.matrix(res) && (nrow(res)==1 || ncol(res)==1)) res <- as.vector(res)
  } else { # if length >1, carefully unlist into a vector so that NULLs/empty values are not lost
    res <- unname(sapply(x, function(y) {
      y <- as.vector(y) # y can be a (single-element or empty) matrix/array, so as.vector
      if (length(y)==0 || is.null(y) || y=="") y <- NA
      y
    }))
  }
  return(res)
})

# gene symbols
recon3d$gene.ids <- recon3d$genes
gid <- str_split(recon3d$gene.ids, "\\.", simplify=TRUE)[,1]
library("org.Hs.eg.db")
mapp <- select(org.Hs.eg.db, keys=gid, columns=c("ENTREZID","SYMBOL"), keytype="ENTREZID") # many:1 mapping (with NA's)
all(gid==mapp$ENTREZID) # TRUE
recon3d$genes <- ifelse(is.na(mapp$SYMBOL), recon3d$gene.ids, mapp$SYMBOL)

# format the "rules" field
recon3d$rules <- unname(sapply(recon3d$rules, function(x) {
  if (is.na(x)) x <- "0"
  str_replace_all(x, "x\\(([0-9]+)\\)", "x\\[\\1\\]")
}))

save(recon3d, file="Recon3D.RData")


### ---- ElegCyc ----

validateSBMLdocument("../models/Celegans_2016_2.xml")
tmp <- readSBMLmod("../models/Celegans_2016_2.xml")

elegcyc <- list()
elegcyc$description <- tmp@mod_desc
elegcyc$modelName <- tmp@mod_name
elegcyc$modelVersion <- tmp@version
elegcyc$comps <- tmp@mod_compart
elegcyc$mets <- tmp@met_id
elegcyc$metNames <- tmp@met_name
elegcyc$metCharges <- tmp@met_attr$charge
elegcyc$metFormulas <- tmp@met_attr$chemicalFormula
elegcyc$rxns <- tmp@react_id
elegcyc$rxnNames <- tmp@react_name
elegcyc$S <- tmp@S
elegcyc$lb <- tmp@lowbnd
elegcyc$ub <- tmp@uppbnd
elegcyc$c <- tmp@obj_coef
elegcyc$rowlb <- rep(0, length(tmp@met_id))
elegcyc$rowub <- rep(0, length(tmp@met_id))
elegcyc$b <- rep(0, length(tmp@met_id))
elegcyc$rxnGeneMat <- tmp@rxnGeneMat
# readSBMLmod() does not construct tmp@gprRules properly!!!
rls <- tmp@gpr
rls <- str_replace_all(rls, " and ", " & ")
rls <- str_replace_all(rls, " or ", " | ")
rls[rls==""] <- "0"
gns <- tmp@allGenes
for (i in 1:length(gns)) {
  rls <- str_replace_all(rls, fixed(gns[i]), paste0("x[",i,"]"))
}
elegcyc$rules <- rls
library(biomaRt)
mart <- useMart(biomart="ensembl", dataset="celegans_gene_ensembl")
tmp1 <- getBM(attributes=c("external_transcript_name", "external_gene_name"), filters="external_transcript_name", values=gns, mart=mart)
gns.map <- tmp1$external_gene_name
names(gns.map) <- tmp1$external_transcript_name
gns <- ifelse(gns %in% names(gns.map), gns.map[gns], gns)
elegcyc$genes <- gns
elegcyc$subSystems <- tmp@subSys

# note: the metabolites that are not involved in any reactions was automatically removed by readSBMLmod(). folate polyglutamate (M_Folatepolyglutamate_n_c or .._m) was also automatically removed as it is a polymer, and it is only involved in reactions where the polymer is elongated or shortened by monomer unit; in such reactions polymers of different lengths (i.e. the reactant and the product) are represented by this same identity, thus amounting to no net change in this identity whatsoever.

save(elegcyc, file="ElegCyc.cleaned.RData")

# 2. from .mat

tmp <- readMat("../models/elegcyc.mat")

# tmp[[1]] contains the model but as an array, each element of the array being a field of the MATLAB struct... so need to use apply
# for some weird reason each element of the array is a deeply nested list, so use rlist::list.flatten to make each of them a simple list of only one level
# cannot use unlist here because it can flattern things directly to vectors and discarding NULLs/empty values.
elegcyc <- apply(tmp[[1]], 1, function(x) list.flatten(x))
# then manually unlist each of the simple lists, so that NULLs/empty values are not lost
elegcyc <- lapply(elegcyc, function(x) {
  # if of length 1, just extract x[[1]], simplify to vector when appropriate
  if (length(x)==1) {
    res <- x[[1]]
    if (is.matrix(res) && (nrow(res)==1 || ncol(res)==1)) res <- as.vector(res)
  } else { # if length >1, carefully unlist into a vector so that NULLs/empty values are not lost
    res <- unname(sapply(x, function(y) {
      y <- as.vector(y) # y can be a (single-element or empty) matrix/array, so as.vector
      if (length(y)==0 || is.null(y) || y=="") y <- NA
      y
    }))
  }
  return(res)
})
elegcyc$modelID <- NULL
elegcyc$metNotes <- NULL
elegcyc$rxnNotes <- NULL
elegcyc$proteins <- NULL
all(elegcyc$genes==elegcyc$geneNames) # TRUE
elegcyc$geneNames <- NULL
elegcyc$genes <- ifelse(elegcyc$genes %in% names(gns.map), gns.map[elegcyc$genes], elegcyc$genes)
# format the "rules" field
elegcyc$rules <- unname(sapply(elegcyc$rules, function(x) {
  if (is.na(x)) x <- "0"
  str_replace_all(x, "x\\(([0-9]+)\\)", "x\\[\\1\\]")
}))
elegcyc$rowlb <- rep(0, length(elegcyc$mets))
elegcyc$rowub <- rep(0, length(elegcyc$mets))
elegcyc$S <- Matrix(elegcyc$S, sparse=TRUE)


save(elegcyc, file="ElegCyc.RData")


### ---- iMM1415 (mice) ----

# from .mat

tmp <- readMat("../models/BiGG/iMM1415.mat")

# tmp[[1]] contains the model but as an array, each element of the array being a field of the MATLAB struct... so need to use apply
# for some weird reason each element of the array is a deeply nested list, so use rlist::list.flatten to make each of them a simple list of only one level
# cannot use unlist here because it can flattern things directly to vectors and discarding NULLs/empty values.
iMM1415 <- apply(tmp[[1]], 1, function(x) list.flatten(x))
# then manually unlist each of the simple lists, so that NULLs/empty values are not lost
iMM1415 <- lapply(iMM1415, function(x) {
  # if of length 1, just extract x[[1]], simplify to vector when appropriate
  if (length(x)==1) {
    res <- x[[1]]
    if (is.matrix(res) && (nrow(res)==1 || ncol(res)==1)) res <- as.vector(res)
  } else { # if length >1, carefully unlist into a vector so that NULLs/empty values are not lost
    res <- unname(sapply(x, function(y) {
      y <- as.vector(y) # y can be a (single-element or empty) matrix/array, so as.vector
      if (length(y)==0 || is.null(y) || y=="") y <- NA
      y
    }))
  }
  return(res)
})

# modify the genes field
gid <- iMM1415$genes # gene entrez ids
# gene ids to gene symbols
library("org.Mm.eg.db")
mapp <- select(org.Mm.eg.db, keys=gid, columns=c("ENTREZID","SYMBOL"), keytype="ENTREZID") # 1:1 mapping (with NA's)
all(gid==mapp$ENTREZID) # TRUE
iMM1415$gene.ids <- gid
iMM1415$genes <- ifelse(is.na(mapp$SYMBOL), gid, mapp$SYMBOL)
# rules mapping genes to reactions
rules <- sapply(iMM1415$grRules, function(x) {
  if (is.na(x)) {
    x <- "0"
  } else {
    x <- str_replace_all(x, "[0-9]+", function(x) paste0("x[", match(x, gid), "]"))
    x <- str_replace_all(x, "and", "&")
    x <- str_replace_all(x, "or", "\\|")
  }
})
iMM1415$rules <- unname(rules)
iMM1415$rowlb <- iMM1415$b # all 0
iMM1415$rowub <- iMM1415$b # all 0
iMM1415$S <- Matrix(iMM1415$S, sparse=TRUE)
iMM1415$csense <- str_split(iMM1415$csense, "")[[1]]

save(iMM1415, file="iMM1415.RData")


