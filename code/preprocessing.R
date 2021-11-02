# Functions for pre-processing of the ST data


# Remove samples with mostly zero genes
# Filter out samples with more than zero.cutoff*100 zero expressed genes
selectSamples <- function(data.tbl, zero.cutoff=0.99) {
  print("Filter samples")
  print(sprintf("Number of samples %d", dim(data.tbl)[2]))
  keep.idx = colSums(data.tbl==0)<=(zero.cutoff*dim(data.tbl)[1])
  data.tbl = data.tbl[,keep.idx]
  print(sprintf("Remaining number of samples %d", dim(data.tbl)[2]))
  return(list("data"=data.tbl, "idx"=keep.idx))
}


# Remove genes with mostly zeros and low variances
# Filter out genes with more than zero.cutoff*100% zero expressed samples
# Filter our genes with lowest var.cutoff*100% variances
selectGenes <- function(data.tbl, zero.cutoff=0.5, var.cutoff=0.1) {
  print("Filter genes")
  print(sprintf("Number of genes %d", dim(data.tbl)[1]))
  keep.idx1 = rowSums(data.tbl==0)<=(zero.cutoff*dim(data.tbl)[2])
  data.tbl = data.tbl[keep.idx1,]
  keep.idx = keep.idx1
  print(sprintf("Remaining number of genes %d", dim(data.tbl)[1]))
  vars = apply(data.tbl, 1, var)
  if (var.cutoff > 0) {
    keep.idx2 = vars>=quantile(vars, var.cutoff)
    data.tbl = data.tbl[keep.idx2,]
    keep.idx = keep.idx[keep.idx2]
  }
  print(sprintf("Remaining number of genes %d", dim(data.tbl)[1]))
  return(list("data"=data.tbl, "idx"=keep.idx))
}


# Normalize the library size of each sample to sum up to scale
NormalizeLibrarySize <- function(data.tbl, scale=10000) {
  return(sweep(data.tbl, 2, colSums(data.tbl), FUN="/") * scale)
}









