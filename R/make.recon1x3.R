# create a tri-cellular model by connecting the extracellular compartment of three Recon 1 models.
# the cell-specific model$genes/mets/rxns will be identified by adding suffixes "_cell1", "_cell2" or "_cell3"
library(stringr)
source("utils.R")
load("Recon1.RData")

# mets in the extracellular space (emet.ids)
tmp <- grepl("\\[e\\]$", recon1$mets)
emet.ids <- which(tmp)
# and the rest (i.e. intracellular)
imet.ids <- which(!tmp)

# exclusively extracellular reactions, these include the extracellular boundary reactions (those with names EX_.*, and equations like "xxx[e]-->") and those involving only extracellular metabolites
tmp <- apply(recon1$S, 2, function(x) all(which(x!=0) %in% emet.ids))
erxn.ids <- which(tmp)
# and the rest
not.erxn.ids <- which(!tmp)

# ---- create combined model ----

# 1. reactions
# concatenate: (a) the current reactions with (b) an extra set of reactions for those in the second and third cell; the not.erxn.ids part of (a) is used to represent those in the first cell, and this part is simply duplicated to obtain (b); the erxn.ids part of (a) is common to all cells
tmpf <- function(a, b=not.erxn.ids) {
  res <- c(recon1[[a]], paste0(recon1[[a]][b], "_cell2"), paste0(recon1[[a]][b], "_cell3"))
  res[b] <- paste0(res[b], "_cell1")
  res
}
rxns <- tmpf("rxns")
rxnNames <- tmpf("rxnNames")
subsys <- tmpf("subSystems")
subsys[subsys %in% c("NA_cell1","NA_cell2","NA_cell3")] <- NA
lb <- c(recon1$lb, rep(recon1$lb[not.erxn.ids], 2))
ub <- c(recon1$ub, rep(recon1$ub[not.erxn.ids], 2))
cvec <- c(recon1$c, rep(recon1$c[not.erxn.ids], 2)) # the default cvec (objective function coefficient)

# 2. metabolites
# concatenate: (a) the current metabolites with (b) an extra set of metabolites for those in the second cell; the imet.ids part of (a) is used to represent those in the first cell, and this part is simply duplicated to obtain (b); the emet.ids part of (a) is common to both cells
mets <- tmpf("mets", imet.ids)
metNames <- tmpf("metNames", imet.ids)
metFor <- c(recon1$metFormulas, rep(recon1$metFormulas[imet.ids], 2))
rowlb <- c(recon1$rowlb, rep(recon1$rowlb[imet.ids], 2))
rowlb[emet.ids] <- 3*rowlb[emet.ids] # x3 since these are produced/consumed by all three cells
rowub <- c(recon1$rowub, rep(recon1$rowub[imet.ids], 2))
rowub[emet.ids] <- 3*rowub[emet.ids] # x3 since these are produced/consumed by all three cells
b <- c(recon1$b, rep(recon1$b[imet.ids], 2))
b[emet.ids] <- 3*b[emet.ids] # x3 since these are produced/consumed by all three cells

# 3. S matrix
# the order of reactions is as above; for the (a) part, the matrix is just the original S unchanged but added more rows of zeros for the intracellular metabolites of the second and third cell; for the (b) part, the matrix is based on the not.erxn.ids columns of the original S matrix just with the rows corresponding to imet.ids "cut and pasted" to bottom
# this is the (b) part:
tmp <- recon1$S[, not.erxn.ids]
tmp[imet.ids, ] <- 0
tmp <- cbind(rbind(tmp, recon1$S[imet.ids, not.erxn.ids], matrix(0,nrow=length(imet.ids),ncol=length(not.erxn.ids))), rbind(tmp, matrix(0,nrow=length(imet.ids),ncol=length(not.erxn.ids)), recon1$S[imet.ids, not.erxn.ids]))
# cbind (a) and (b)
S <- cbind(rbind(recon1$S, matrix(0,nrow=length(imet.ids)*2,ncol=ncol(recon1$S))), tmp)

# 4. genes
# just duplicated everything since each cell has its own genes
genes <- c(paste0(recon1$genes, "_cell1"), paste0(recon1$genes, "_cell2"), paste0(recon1$genes, "_cell3"))
gene.ids <- c(paste0(recon1$gene.ids, "_cell1"), paste0(recon1$gene.ids, "_cell2"), paste0(recon1$gene.ids, "_cell3"))

# 5. rules (not using "grRules" or "rxnGeneMat", and these are also not included in the resulting model)
# just need to correct the gene indices for the second and third cell
tmp <- recon1$rules[not.erxn.ids]
tmp <- str_replace_all(tmp, "[0-9]+", function(x) {
  x <- as.integer(x)
  if (x==0) 0 else x+length(recon1$genes)
})
rules <- c(recon1$rules, tmp)
tmp <- str_replace_all(tmp, "[0-9]+", function(x) {
  x <- as.integer(x)
  if (x==0) 0 else x+length(recon1$genes)
})
rules <- c(rules, tmp)
# if any genes responsible for erxn.ids, then need to recreate rules for these reactions as "(rules.cell1) | (rules.cell2) | (rules.cell3)"
tmp <- rxns2genes(erxn.ids, recon1)
dupg.rxn.ids <- erxn.ids[sapply(tmp, length)!=0]
tmp <- str_replace_all(rules[dupg.rxn.ids], "[0-9]+", function(x) {
  x <- as.integer(x)
  if (x==0) 0 else x+length(recon1$genes)
})
rules[dupg.rxn.ids] <- paste0("(",rules[dupg.rxn.ids],") | (",tmp,")")
tmp <- str_replace_all(tmp, "[0-9]+", function(x) {
  x <- as.integer(x)
  if (x==0) 0 else x+length(recon1$genes)
})
rules[dupg.rxn.ids] <- paste0(rules[dupg.rxn.ids], " | (",tmp,")")


# resulting model
recon1x3 <- list(description="Recon 1 Tricellular", rxns=rxns, rxnNames=rxnNames, subSystems=subsys, mets=mets, metNames=metNames, metFormulas=metFor, lb=lb, ub=ub, c=cvec, rowlb=rowlb, rowub=rowub, b=b, S=S, genes=genes, gene.ids=gene.ids, rules=rules, imet.ids=imet.ids, emet.ids=emet.ids, irxn.ids=not.erxn.ids, erxn.ids=erxn.ids)

save(recon1x3, file="Recon1x3.RData")



