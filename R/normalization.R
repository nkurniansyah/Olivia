#'Median Normalization
#'
#' Median normalization is using median value for each library so they have equal size
#'
#' @param count_matrix  A matrix of gene counts (possibly transformed). rows are genes, columns are individuals
#' @return median_norm  A matrix of gene counts  after normalization
#' @examples
#' data(rnaseq_count_matrix)
#' rnaseq_count_matrix<- rnaseq_count_matrix[rowSums(rnaseq_count_matrix)>0,]
#' median_normalization(rnaseq_count_matrix)
#' @export
median_normalization <- function(count_matrix){
  median_norm <- t(t(count_matrix)/(colSums(count_matrix))*median(colSums(count_matrix)))
  return(median_norm)
}