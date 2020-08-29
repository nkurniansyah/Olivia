

#' Title apply log transform on a matrix of transcript counts
#'
#' @param count_matrix : A p x n matrix of gene expression counts (possibly transformed). p are genes, n are individuals. Rownames are gene names.
#' @param log_transform_type  : One of the transformations log_replace_half_min, log_add_min, log_add_0.5, or NULL 
#' @return  A p x n matrix of gene expression counts after transformation
#' @export
#'
#' @examples
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