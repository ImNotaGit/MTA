library(R.matlab)
library(stringr)

tmp <- readMat("../matlab/models/recon.mat")

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

# format the "rules" field
recon1$rules <- unname(sapply(recon1$rules, function(x) {
  if (is.na(x)) x <- "0"
  str_replace_all(x, "x\\(([0-9]+)\\)", "x\\[\\1\\]")
}))

# change model description
recon1$description <- "Recon 1"

save(recon1, file="Recon1.RData")
