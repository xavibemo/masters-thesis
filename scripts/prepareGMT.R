## LIBRARIES
library(fgsea)


#' @title Matrix to list transformation
#' Adapted from https://biostatsquid.com/fgsea-tutorial-gsea/
#' 
#' @param pws matrix like object
#' 
#' @returns transformed list

matrix_to_list <- function(pws)
{
  pws.l <- list()
  
  for (pw in colnames(pws)) 
  {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}


#' @title Background preparation for GSEA analysis
#' Adapted from https://biostatsquid.com/fgsea-tutorial-gsea/
#' 
#' @param gmt_file path to .gmt file
#' @param genes_in_data array-like object of all genes in the study
#' 
#' @returns a final list of all gene sets and the correct gene backgrounds

prepareGMT <- function(gmt_file, genes_in_data)
{
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  
  for (i in 1:dim(mat)[2])
  {
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  # Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  # filter for gene sets with more than 5 genes annotated
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] 
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  return(final_list)
}