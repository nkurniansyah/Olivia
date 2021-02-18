
#' Applying log transform on a matrix of gene counts
#'
#' @param count_matrix  A matrix of gene counts (possibly transformed). rows are genes, columns are individuals
#' @param transform One of the transformations log_replace_half_min, log_add_min, log_add_0.5
#' @return  A matrix of gene counts after transformation
#' @examples
#' data(rnaseq_count_matrix)
#' rnaseq_count_matrix<- rnaseq_count_matrix[rowSums(rnaseq_count_matrix)>0,]
#' log_transform_count(count_matrix=rnaseq_count_matrix, transform = "log_replace_half_min")
#' @export
#'
#'

log_transform_count <- function(count_matrix, transform = "log_replace_half_min"){
  if (is.null(transform)){
    message("No transformation of gene counts requested")
    return(count_matrix)
  }
  if (!is.element(transform, c("log_replace_half_min", "log_add_min", "log_add_0.5"))){
    stop("Requested transformation not allowed.
         Allowed transformation names are log_replace_half_min, log_add_min, log_add_0.5")
  }
  if (transform == "log_replace_half_min") {
    return(log_replace_half_min(count_matrix))
  }

  if (transform == "log_add_min"){
    return(log_add_min(count_matrix))
  }

  if (transform == "log_add_0.5"){
    return(log_add_0.5(count_matrix))
  }


  }




#' Log replace half min
#'
#' Type of log transform, replacing zero with minimum value for a gene divided by 2
#'
#' @param count_matrix  A matrix of gene counts (possibly transformed). rows are genes, columns are individuals
#'
#' @return A matrix of gene counts after transformation
#' @examples
#' data(rnaseq_count_matrix)
#' rnaseq_count_matrix<- rnaseq_count_matrix[rowSums(rnaseq_count_matrix)>0,]
#' log_replace_half_min(count_matrix=rnaseq_count_matrix)
#' @export

#

log_replace_half_min<- function(count_matrix){
  imputed_mat <- t(apply(count_matrix,1,function(x){x[x==0] <- min(x[x>0]/2);x}))
  imputed_mat <- log2(imputed_mat)
  return(imputed_mat)
}


#' Log add half
#'
#' Type of log transform, adding 0.5 into entire gene counts matrix
#'
#' @param count_matrix  A matrix of gene counts (possibly transformed). rows are genes, columns are individuals
#'
#' @return A matrix of gene expression counts after transformation
#' @examples
#' data(rnaseq_count_matrix)
#' rnaseq_count_matrix<- rnaseq_count_matrix[rowSums(rnaseq_count_matrix)>0,]
#' log_add_0.5(count_matrix=rnaseq_count_matrix)
#' @export

log_add_0.5<- function(count_matrix){
  imputed_mat <- log2(count_matrix + 0.5)
  return(imputed_mat)
}


#' Log add minimum value
#'
#' Type of log transform, replacing zero with minimum values for each genes
#'
#' @param count_matrix  A matrix of gene counts (possibly transformed). rows are genes, columns are individuals
#'
#' @return A matrix of gene expression counts after transformation
#' @examples
#' data(rnaseq_count_matrix)
#' rnaseq_count_matrix<- rnaseq_count_matrix[rowSums(rnaseq_count_matrix)>0,]
#' log_add_min(count_matrix=rnaseq_count_matrix)
#' @export

log_add_min <-  function(count_matrix){
  imputed_mat = t(apply(count_matrix,1,function(x){x= x + min(x[x>0]/2);x}))
  imputed_mat <- log2(imputed_mat)
  return(imputed_mat)
}
