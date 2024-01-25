#' @title TPM Transformation
#' @author Xavier Benedicto Molina
#' Adapted from https://support.bioconductor.org/p/91218/.
#' 
#' @param counts matrix like object containing raw counts
#' @param length array like object containing transcript lengthts
#' 
#' @returns matrix like object with TPMs

TPMmatrix <- function(counts, length) 
{
  x <- counts / length
  return( t( t(x)*1e6 / colSums(x) ) )
}
