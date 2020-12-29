

#' Title Normalize a matrix of transcript counts
#'
#' @param count_matrix : A p x n matrix of gene expression counts (possibly transformed). p are genes, n are individuals. Rownames are gene names.
#' @param normalized_type  : One of the normalization method median_normalization, SizeFactor, TMM, or NULL 
#' @return  A p x n matrix of gene expression counts after normalization
#' @export
#' 

normalize_trancript_count <- function(count_matrix,covariates_string=NULL, outcome=NULL, phenotype=NULL, normalized_type = "median_normalization"){
  if (is.null(median_normalization)){
    message("No normalization method of gene counts requested")
    return(count_matrix)
  }
  if (!is.element(normalized_type, c("median_normalization", "SizeFactor", "TMM"))){
    stop("Requested normalization not allowed. 
         Allowed normalization names are median_normalization, SizeFactor, TMM")
  }
  
  if (normalized_type == "median_normalization") {
    message("Applying median normalization  ...")
    
    return(median_normalization(count_matrix))
  }
  
  if (normalized_type == "SizeFactor"){
    message("Applying size factor ...")
    return(SizeFactor(count_matrix,covariates_string,outcome,phenotype))
  }
  
  if (normalized_type == "TMM"){
    message("Applying TMM ...")
    
    return(TMM(count_matrix))
  }
  
}