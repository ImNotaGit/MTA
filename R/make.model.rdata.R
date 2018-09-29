library(R.matlab)
library(stringr)

tmp <- readMat("../matlab/models/recon.mat")

recon1 <- apply(tmp[[1]], 1, function(x) unname(unlist(x)))
recon1$S <- recon1$S[[1]]
recon1$rxnGeneMat <- recon1$rxnGeneMat[[1]]
recon1$gene.ids <- recon1$genes
recon1$genes <- recon1$genes.unique.names[recon1$genes.unique.map]
recon1$genes.unique <- NULL
recon1$genes.unique.map <- NULL
recon1$genes.unique.names <- NULL
recon1$rules <- unname(sapply(recon1$rules, function(x) {
  if (length(x)==0) x <- "0"
  str_replace_all(x, "x\\(([0-9]+)\\)", "x\\[\\1\\]")
}))
recon1$description <- "Recon 1"

save(recon1, file="Recon1.RData")
